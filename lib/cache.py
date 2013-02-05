import sys
import commands
import re
import os
import shutil
import hashlib
import string
import tempfile
import glob

import Queue
import threading

from cogent.db.ensembl import Species, Genome, Compara

from lib.progress import Progress
from lib.tools import Prank


class EnsemblInfo(object) :
    def __init__(self) :
        pass

    def _convert_to_range(self, releases) :
        releases.sort()
        return "%d-%d" % (releases[0], releases[-1])

    def _get_latest_release_versions(self) :
        showdb = "mysql -h ensembldb.ensembl.org -u anonymous -P5306 -B -e 'SHOW DATABASES;'"
        stat,output = commands.getstatusoutput(showdb)

        if stat != 0 :
            print >> sys.stderr, "Error: could not run \"%s\"" % showdb
            sys.exit(-1)

        dbpat = re.compile("^(.*)_core_(\d+)_.*")
        db2rel = {}

        for dbdesc in output.split('\n') :
            if "core" in dbdesc :
                m = dbpat.match(dbdesc)
                if (m != None) and (len(m.groups()) == 2) :
                    dbname,dbrel = m.groups()
                    if dbname not in db2rel :
                        db2rel[dbname] = []
                    db2rel[dbname].append(int(dbrel))

        for dbname in db2rel :
            db2rel[dbname] = self._convert_to_range(db2rel[dbname])

        return db2rel

    def print_species_table(self) :
        db2rel = self._get_latest_release_versions()

        print "Name".rjust(30), "Common Name".rjust(20), "Releases".rjust(12)
        print "-" * 64

        for name in Species.getSpeciesNames() :
            dbprefix = Species.getEnsemblDbPrefix(name)
            if db2rel.has_key(dbprefix) :
                print name.rjust(30), 
                print Species.getCommonName(name).rjust(20), 
                print db2rel.get(dbprefix, "???").rjust(12)


class TranscriptCache(object) :
    file_manifest = 'manifest'
    file_prefix = 'genefamily_'

    def __init__(self, workingdir, tmpdir, species, release) :
        self.stop = False
        self.workingdir = workingdir
        self.tmpdir = tmpdir
        self.species = species
        self.release = release
        self.account = None
        self.basedir = os.path.join(self.workingdir, str(release), species)

        self._check_directory(self.tmpdir, create=True)
        self._check_directory(self.workingdir, create=True)
        self._check_directory(os.path.join(self.workingdir, str(release)), create=True)
        self._check_directory(os.path.join(self.workingdir, str(release), species), create=True)

        self.genes = set()
        self._verify_manifest()

        # alignments are handled in one thread
        self.alignment_queue = Queue.Queue()
        self.alignment_thread = threading.Thread(target=self._consume_alignment_queue)
        self.alignment_thread.daemon = True
        self.alignment_thread.start()

        # there is a separate thread to write the manifest and files
        self.manifest_queue = Queue.Queue()
        self.manifest_thread = threading.Thread(target=self._consume_manifest_queue)
        self.manifest_thread.daemon = True
        self.manifest_thread.start()

    def shutdown(self) :
        self.stop = True

    def _consume_alignment_queue(self) :
        while not self.stop :
            fname = self.alignment_queue.get()
            self._align(fname)
            self.alignment_queue.task_done()

    def _consume_manifest_queue(self) :
        while not self.stop :
            data = self.manifest_queue.get()
           
            if len(data) != 3 :
                print >> sys.stderr, "Error: manifest queue contained %s" % str(data)
                sys.exit(-1)

            fname = self._write_file(data[0], data[1])

            if data[2] :
                self._add_to_alignment_queue(fname)

            self.manifest_queue.task_done()

    def _md5(self, s) :
        return hashlib.md5(s).hexdigest()
        
    def _contents(self, fname) :
        s = ""
        f = open(fname)
        
        for line in f :
            s += line

        f.close()
        return s

    def _write_file(self, filename, filecontents) :
        manifestdata = self._md5(filecontents)
        
        # write md5 hash to the end of the manifest file
        f = open(self.manifest_name, 'a')
        print >> f, "%s %s" % (filename, manifestdata)
        f.close()

        # write the actual file, we don't care if this gets interrupted
        # because then it will be discovered by recalculating the hash
        f = open(os.path.join(self.basedir, filename), 'w')
        f.write(filecontents)
        f.close()

        print "Info: written %s" % filename
        return os.path.join(self.basedir, filename)

    def _add_to_manifest_queue(self, filename, filecontents, align) :
        self.manifest_queue.put((filename, filecontents, align))

    def _add_to_alignment_queue(self, filename) :
        self.alignment_queue.put(filename)

    def _check_directory(self, dirname, create=False) :
        if not os.path.exists(dirname) :
            if create :
                try :
                    os.makedirs(dirname)
                    print "Info: created %s" % dirname
                except OSError, ose :
                    print >> sys.stderr, "Error: %s" % str(ose)
                    sys.exit(-1)
            else :
                print >> sys.stderr, "Error: '%s' does not exist." % dirname
                sys.exit(-1)

        elif not os.path.isdir(dirname) :
            print >> sys.stderr, "Error: '%s' exists, but is not a directory!" % dirname
            sys.exit(-1)

        else :
            pass
        
    def _reset_cache(self) :
        #shutil.rmtree(self.basedir, ignore_errors=True)

        #try:
        #    os.rmdir(os.path.dirname(self.basedir))
        #except OSError, ose :
        #    pass

        for fname in glob.glob(os.path.join(self.basedir, '*')) :
            os.remove(fname)

        open(self.manifest_name, 'w').close()

    def _swap_dirname(self, fname, directory=None) :
        return os.path.join(self.basedir if not directory else directory, os.path.basename(fname))

    def _file_md5(self, fname) :
        return self._md5(self._contents(fname))

    def _file_md5_and_genes(self, fname) :
        md = md5.new()
        tmp = []
        f = open(fname)

        for line in f :
            md.update(line)
            
            if line.startswith('>') :
                tmp.append(line.strip()[1:])

        f.close()
        return md.hexdigest(), tmp

    @property
    def manifest_name(self) :
        return os.path.join(self.basedir, type(self).file_manifest)

    # TODO
    # verify that files match md5s in manifest
    # figure out which genes cds has been downloaded and stored
    # figure out which alignments have not been done and add to the prank queue
    def _verify_manifest(self) :
        pat = re.compile("^(.+) ([" + string.ascii_lowercase + string.digits + "]{32})$")
        f2md5 = {}

        # read manifest into memory
        linenum = 0

        try :
            f = open(self.manifest_name)
        except IOError, ioe :
            print "Info: creating new manifest file..."
            open(self.manifest_name, 'w').close()
            return

        for line in f :
            linenum += 1
            line = line.strip()

            if line == '' :
                continue

            result = pat.match(line)

            if not result :
                print >> sys.stderr, "Warning: line %d in %s appears to be corrupt..." % (linenum, f.name)
                continue

            f2md5[result.group(1)] = result.group(2)

        f.close()

        # check md5 sums
        correct_md5 = []
        bad_md5 = []
        for fname in f2md5 :
            try :
                if self._file_md5(os.path.join(self.basedir, fname)) == f2md5[fname] :
                    correct_md5.append(f2md5[fname])
                else :
                    bad_md5.append(f2md5[fname])

            except IOError, ioe :
                bad_md5.append(f2md5[fname])

        # TODO
        # based on what the file is, do file specific checks
        # ie: if the data is present, but the alignment is not, then queue the alignment
        print "Info: TODO file type specific checks..."

    def _random_filename(self) :
        return os.path.basename(tempfile.mktemp(prefix=type(self).file_prefix, dir=self.basedir))

    def _align(self, infile) :
        outfile = self._swap_dirname(infile, directory=self.tmpdir)
        
        outfiles = Prank(self.tmpdir).align(infile, outfile)
        
        for f in outfiles :
            self._add_to_manifest_queue(os.path.basename(f), self._contents(f), False)

        print "Info: aligned sequences in %s" % infile

    def build(self, resume=True) :
        if not resume :
            self._reset_cache()
        else :
            print "Info: resuming..."

        print "Info: enumerating gene families in %s release %d" % (self.species, self.release)
        genome = Genome(self.species, Release=self.release, account=self.account)
        compara = Compara([self.species], Release=self.release, account=self.account)

        for gene in genome.getGenesMatching(BioType='protein_coding') :
            stableid = gene.StableId.lower()

            # ignore genes that have already been seen as members of
            # gene families
            if stableid in self.genes :
                continue

            self.genes.add(stableid)

            # get cds sequences of any paralogs
            paralogs = compara.getRelatedGenes(StableId=stableid, Relationship='within_species_paralog')

            paralog_seqs = []
            if paralogs is None :
                paralog_seqs.append((stableid, gene.getLongestCdsTranscript().Cds))
            else :
                for paralog in paralogs.Members :
                    paralog_id = paralog.StableId.lower()
                    self.genes.add(paralog_id)
                    paralog_seqs.append((paralog_id, paralog.getLongestCdsTranscript().Cds))

            # write md5 to the manifest and save to disk
            fname = self._random_filename()
            fcontents = ""
            for geneid,cds in paralog_seqs :
                fcontents += (">%s\n%s\n" % (geneid, cds))

            self._add_to_manifest_queue(fname, fcontents, len(paralog_seqs) >= 2)

            print "%s - %d genes in family (%s)" % (stableid, len(paralog_seqs), fname)

            # check to see if we have been told to stop
            if self.stop :
                break

        print "Killed by user" if self.stop else "Database downloaded..."

