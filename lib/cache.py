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

from cogent.parse.fasta import MinimalFastaParser
from cogent.db.ensembl import Species, Genome, Compara, HostAccount
from cogent.db.ensembl.database import Database

from lib.progress import Progress
from lib.tools import Prank, PrankError


class EnsemblInfo(object) :
    def __init__(self, options) :
        self.db_host = options['db-host']
        self.db_port = options['db-port']
        self.db_user = options['db-user']
        self.db_pass = options['db-pass']

    def _convert_to_range(self, releases) :
        releases.sort()
        return "%d-%d" % (releases[0], releases[-1])

    def _get_latest_release_versions(self, only_latest=False) :
        #showdb = "mysql -h ensembldb.ensembl.org -u anonymous -P5306 -B -e 'SHOW DATABASES;'"
        #showdb = "mysql -h mysql.ebi.ac.uk -u anonymous -P4157 -B -e 'SHOW DATABASES;'"
        passwd = "" if self.db_pass == "" else "-p %s" % self.db_pass
        showdb = "mysql -h %s -u %s -P %d %s -B -e 'SHOW DATABASES;'" % (self.db_host, self.db_user, self.db_port, passwd)
        
        stat,output = commands.getstatusoutput(showdb)

        if stat != 0 :
            print >> sys.stderr, "Error: could not run \"%s\"" % showdb
            sys.exit(-1)

        dbpat = re.compile("^(.*)_core_(\d+_)?(\d+)_.+")
        db2rel = {}

        for dbdesc in output.split('\n') :
            if "core" in dbdesc :
                m = dbpat.match(dbdesc)
                if (m != None) and (len(m.groups()) == 3) :
                    dbname,chaff,dbrel = m.groups()
                    if dbname not in db2rel :
                        db2rel[dbname] = []
                    if chaff is None :
                        db2rel[dbname].append(int(dbrel))
                    else :
                        # in the case of the ensembl-metazoa species
                        db2rel[dbname].append(int(chaff[:-1]))
                #else :
                #    if "core" in dbdesc :
                #        print "rejected", dbdesc

        for dbname in db2rel :
            if only_latest :
                db2rel[dbname] = max(db2rel[dbname])
            else :
                db2rel[dbname] = self._convert_to_range(db2rel[dbname])

        return db2rel

    def get_latest_release(self, species) :
        db2rel = self._get_latest_release_versions(only_latest=True)
        try :
            return db2rel.get(Species.getEnsemblDbPrefix(species))
        except :
            pass
            
        tokens = map(lambda x : x.lower(), species.split())
        dbname = '_'.join(tokens)
 
        return db2rel.get(dbname, -1)

    def print_species_table(self) :
        db2rel = self._get_latest_release_versions()
        printed = {}

        print "Name".rjust(32), "Common Name".rjust(18), "Releases".rjust(12)
        print "-" * 64

        for name in Species.getSpeciesNames() :
            dbprefix = Species.getEnsemblDbPrefix(name)
            if db2rel.has_key(dbprefix) :
                print name.rjust(32), 
                print Species.getCommonName(name).rjust(18), 
                print db2rel.get(dbprefix).rjust(12)
                printed[dbprefix] = True

        print "-" * 64

        for dbprefix in sorted(db2rel.keys()) :
            if dbprefix not in printed :
                print dbprefix.capitalize().replace('_', ' ').rjust(32),
                print "undefined".rjust(18),
                print db2rel.get(dbprefix).rjust(12),
                print dbprefix.rjust(30)


class TranscriptCache(object) :
    file_manifest = 'manifest'
    file_prefix = 'genefamily_'
    queue_timeout = 1

    def __init__(self, options) :
        self.stop = False
        self.download_complete = False
        self.alignments_complete = False
        self.workingdir = options['workingdir']
        self.tmpdir = options['tmpdir']
        self.species = options['species']
        self.release = options['release']
        self.account = HostAccount(options['db-host'], options['db-user'], options['db-pass'], port=options['db-port'])
        self.basedir = os.path.join(self.workingdir, str(self.release), self.species)
        self.resume = options['resume']
        self.database = options['database'] # this is necessary for a hack used later

        self._is_valid_species()
        #Species.amendSpecies("Tribolium castaneum", "T.castaneum")

        self._check_directory(self.tmpdir, create=True)
        self._check_directory(self.workingdir, create=True)
        self._check_directory(os.path.join(self.workingdir, str(self.release)), create=True)
        self._check_directory(os.path.join(self.workingdir, str(self.release), self.species), create=True)

        self.genes = set()

        # alignments are handled in one thread
        self.alignment_queue = Queue.Queue()
        self.alignment_thread = threading.Thread(target=self._consume_alignment_queue)
        self.alignment_thread.daemon = True

        # there is a separate thread to write the manifest and files
        self.manifest_queue = Queue.Queue()
        self.manifest_thread = threading.Thread(target=self._consume_manifest_queue)
        self.manifest_thread.daemon = True

        if not self.resume :
            self._reset_cache()
        else :
            self._verify_manifest()

        self.alignment_thread.start()
        self.manifest_thread.start()

    def shutdown(self) :
        self.stop = True

    # this is a hack, but will often work
    def latin2common(self, latin_name) :
        tokens = latin_name.lower().split()
        return tokens[0][0].capitalize() + "." + tokens[1]

    def _is_valid_species(self) :
        try :
            Species.getCommonName(self.species)
            return

        except :
            if 'collection' in self.species.lower() :
                print >> sys.stderr, "Error: species containing the word 'collection' do not work..."
                sys.exit(-1)

        common = self.latin2common(self.species)
        print >> sys.stderr, "Info: Requested species '%s' not found in pycogent, adding as '%s' (this might not work...)" % (self.species, common)
        
        Species.amendSpecies(self.species, common)

    def _consume_alignment_queue(self) :
        while not self.stop :
            try :
                fname = self.alignment_queue.get(timeout=type(self).queue_timeout)
            
            except Queue.Empty :
                if self.download_complete :
                    break
                continue

            self._align(fname)
            self.alignment_queue.task_done()

        self.alignments_complete = True
        print "Info: alignment thread finished"

    def _consume_manifest_queue(self) :
        while not self.stop :
            try :
                data = self.manifest_queue.get(timeout=type(self).queue_timeout)
           
            except Queue.Empty :
                if self.alignments_complete :
                    break
                continue

            if len(data) != 3 :
                print >> sys.stderr, "Error: manifest queue contained %s" % str(data)
                sys.exit(-1)

            fname = self._write_file(data[0], data[1])

            if data[2] :
                self._add_to_alignment_queue(fname)

            self.manifest_queue.task_done()

        print "Info: manifest thread finished"

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
        os.fsync(f)
        f.close()

        # write the actual file, we don't care if this gets interrupted
        # because then it will be discovered by recalculating the hash
        f = open(os.path.join(self.basedir, filename), 'w')
        f.write(filecontents)
        os.fsync(f)
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

    def _verify_manifest(self) :
        manifest_pat = re.compile("^(.+) ([" + string.ascii_lowercase + string.digits + "]{32})$") # filename + md5
        family_pat = re.compile("^" + type(self).file_prefix + "[" + string.ascii_letters + string.digits + "]{6}$") # filename for CDS data
        f2md5 = {}
        linenum = 0

        # 1. parse manifest file 
        #    store copy in memory in f2md5
        try :
            f = open(self.manifest_name)
        except IOError, ioe :
            return

        print "Info: verifying manifest file..."

        for line in f :
            linenum += 1
            line = line.strip()

            if line == '' :
                continue

            result = manifest_pat.match(line)

            if not result :
                print >> sys.stderr, "Warning: line %d in %s appears to be corrupt..." % (linenum, f.name)
                continue

            f2md5[result.group(1)] = result.group(2)

        f.close()

        # 2. check md5 sums of gene families containing
        #    CDS data
        good_gene_families = []
        
        for fname in f2md5 :
            if not family_pat.match(fname) :
                continue

            try :
                if self._file_md5(os.path.join(self.basedir, fname)) == f2md5[fname] :
                    good_gene_families.append(fname)
                else :
                    print >> sys.stderr, "Info: md5 for %s was bad..." % fname

            except IOError, ioe :
                print >> sys.stderr, "Info: %s not found!" % fname
                continue

        # 3. for all the gene family files with good md5s
        #    if there is more than one sequence, check that 
        #    alignments are present and md5s match
        good_alignments = []

        for fname in good_gene_families :
            if self._count_sequences(os.path.join(self.basedir, fname)) < 2 :
                continue

            alignment_files = map(lambda x: fname + x, [".1.dnd", ".2.dnd", ".nuc.1.fas", ".nuc.2.fas", ".pep.1.fas", ".pep.2.fas"])
            realign = False

            for align_fname in alignment_files :
                try :
                    if self._file_md5(os.path.join(self.basedir, align_fname)) != f2md5[align_fname] : # bad md5
                        print >> sys.stderr, "Info: md5 for %s was bad, queuing for alignment..." % align_fname
                        realign = True
                        break
                
                except IOError, ioe : # no file to open
                    print >> sys.stderr,  "Info: %s alignment files missing, queuing for alignment..." % fname
                    realign = True
                    break

                except KeyError, ke : # no mention in the manifest
                    print >> sys.stderr, "Info: md5 for %s not present in manifest, queuing for alignment..." % align_fname
                    realign = True
                    break

            if realign :
                self._add_to_alignment_queue(os.path.join(self.basedir, fname))
            else :
                good_alignments += alignment_files


        good_files = good_gene_families + good_alignments
        # 4. rewrite manifest
        #    XXX what if the program is killed during?
        f = open(self.manifest_name, 'w')

        for fname in good_files :
            if family_pat.match(fname) :
                self._add_genes(os.path.join(self.basedir, fname))
            print >> f, "%s %s" % (fname, f2md5[fname])

        f.close()

        # 5. unlink chaff
        for fname in glob.glob(os.path.join(self.basedir, '*')) :
            basename = os.path.basename(fname)

            if basename == type(self).file_manifest :
                continue

            if basename not in good_files :
                print >> sys.stderr, "Info: removing %s ..." % basename
                os.remove(fname)

        print "Info: verification complete (%d genes in %d families)" % (len(self.genes), len(good_gene_families))

    def _count_sequences(self, fname) :
        count = 0
        for label,seq in MinimalFastaParser(open(fname)) :
            count += 1
        return count

    def _add_genes(self, fname) :
        f = open(fname)

        for line in f :
            if line.startswith('>') :
                self.genes.add(line.strip()[1:])

        f.close()

    def _random_filename(self) :
        return os.path.basename(tempfile.mktemp(prefix=type(self).file_prefix, dir=self.basedir))

    def _align(self, infile) :
        outfile = self._swap_dirname(infile, directory=self.tmpdir)
        
        print "Prank: aligning %s ..." % (os.path.basename(infile))

        try :
            outfiles = Prank(self.tmpdir).align(infile, outfile)
        
        except PrankError, pe :
            print >> sys.stderr, "Error: Prank died on %s ..." % (os.path.basename(infile))
            return

        for f in outfiles :
            self._add_to_manifest_queue(os.path.basename(f), self._contents(f), False)

        print "Info: aligned sequences in %s" % os.path.basename(infile)

    def build(self) :
        if self.resume :
            print "Info: resuming..."

        print "Info: enumerating gene families in %s release %d" % (self.species, self.release)
        genome = Genome(self.species, Release=self.release, account=self.account)
        compara = Compara([self.species], Release=self.release, account=self.account)



        # DON'T DO THIS AT HOME!
        #
        # what happens is it searches for compara databases, but unfortunately finds more than one
        # in this situation pycogent just connects to the first one, which is always compara_bacteria
        # so one solution is to dig through all the compara objects internals to provide a connection
        # to the correct database ... obviously not the best solution, but at 6 lines of code definitely 
        # the shortest ;-P
        #
        if self.database == 'ensembl-metazoa' :
            from cogent.db.ensembl.host import DbConnection
            from cogent.db.ensembl.name import EnsemblDbName
            import sqlalchemy

            new_db_name = EnsemblDbName(compara.ComparaDb.db_name.Name.replace('bacteria', 'metazoa'))
            compara.ComparaDb._db = DbConnection(account=self.account, db_name=new_db_name)
            compara.ComparaDb._meta = sqlalchemy.MetaData(compara.ComparaDb._db)
        # end of DON'T DO THIS AT HOME!



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

            print "Info: %s - %d genes in family (%s)" % (stableid, len(paralog_seqs), fname)

            # check to see if we have been told to stop
            if self.stop :
                break

            # XXX DEBUG
            if len(self.genes) > 10 :
                break


        if self.stop :
            print "Info: killed by user..."
        else :
            print "Info: download complete..."
            self.download_complete = True
            self.alignment_thread.join()
            self.manifest_thread.join()
            print "Info: done!"

