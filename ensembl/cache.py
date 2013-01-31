import sys
import commands
import re
import os
import shutil
import md5

from cogent.db.ensembl import Species, Genome, Compara

from ensembl.progress import Progress


class EnsemblCache(object) :
    def __init__(self, workingdir, tmpdir) :
        self.stop = False
        self.workingdir = workingdir
        self.tmpdir = tmpdir

        self._check_directory(self.workingdir, create=True)
        self._check_directory(self.tmpdir, create=True)

        self.species = None
        self.release = None
        self.account = None
        self.basedir = None
        self.genes = set()

    def shutdown(self) :
        self.stop = True

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

    def _check_directory(self, dirname, create=False, silent=False) :
        if not os.path.exists(dirname) :
            if create :
                #print >> sys.stderr, "Info: '%s' does not exist, creating..." % dirname
                try :
                    os.makedirs(dirname)
                except OSError, ose :
                    print >> sys.stderr, str(ose)
                    return False
            else :
                if not silent :
                    print >> sys.stderr, "Error: '%s' does not exist." % dirname
                return False

        elif not os.path.isdir(dirname) :
            if not silent :
                print >> sys.stderr, "Error: '%s' exists, but is not a directory!" % tmp
            return False

        else :
            return True
        
    def _check_cache_directory(self, species, release) :
        tmp = self.workingdir + os.sep + str(release)
        self._check_directory(tmp, create=True)

        tmp += (os.sep + Species.getCommonName(species))
        self._check_directory(tmp, create=True)

        return tmp

    def _kill_cache(self, species, release) :
        tmp = self.workingdir + os.sep + str(release)
        if self._check_directory(tmp, silent=True) :
            tmp += (os.sep + Species.getCommonName(species))
            if self._check_directory(tmp, silent=True) :
                shutil.rmtree(tmp, ignore_errors=True)

    def _file_md5(self, fname) :
        md = md5.new()
        f = open(fname)

        for line in f :
            md.update(line)

        f.close()
        return md.hexdigest()

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

    def _md5(self, s) :
        return md5.new(s).hexdigest()
        
    def _write_file_and_manifest(self, basedir, filename, filecontents) :
        manifestdata = self._md5(filecontents)
        
        # write md5 hash to the end of the manifest file
        f = open(basedir + os.sep + "manifest", 'a')
        print >> f, "%s %s" % (filename, manifestdata)
        f.close()        

        # write the actual file, we don't care if this gets interrupted
        # because then it will be discovered by recalculating the hash
        f = open(basedir + os.sep + filename, 'w')
        f.write(filecontents)
        f.close()

    def _verify_manifest(self) :
        f2md5 = {}
        
        # read manifest file
        linenum = 0
        f = open(self.basedir + os.sep + "manifest")

        for line in f :
            linenum += 1
            line = line.strip()
            
            if line == '' :
                continue
            
            data = line.split()

            if len(data) != 2 :
                print >> sys.stderr, "Warning: line %d in %s appears to be corrupt..." % (linenum, f.name)
                continue

            f2md5[data[0]] = data[1]

        f.close()

        # check md5 sums
        p = Progress("Verifying", len(f2md5))
        p.start()

        for fname in f2md5 :
            md5hash,genelist = self._file_md5_and_genes(self.basedir + os.sep + fname)
            if md5hash == f2md5[fname] :
                self.genes.update(genelist)
            p.increment()

        p.end()

    def build_cache(self, species, release, resume=True) :

        self.species = species
        self.release = release
        self.basedir = self.workingdir + os.sep + str(self.release) + os.sep + self.species

        if not resume :
            self._kill_cache(species, release)

        self._check_cache_directory(species, release)
        
        if resume :
            self._verify_manifest()

        genome = Genome(self.species, Release=self.release, account=self.account)
        compara = Compara([self.species], Release=self.release, account=self.account)

        gene_families = []

        for gene in genome.getGenesMatching(BioType='protein_coding') :
            gene_id = gene.StableId.lower()

            # ignore genes that have already been seen as members of
            # gene families
            if gene_id in self.genes :
                continue

            print gene_id
            self.genes.add(gene_id)

            #paralogs = compara.getRelatedGenes(StableId=gene_id, Relationship='within_species_paralog')

            #if paralogs is None :
            #    gene_families.append([gene_id])
            #else :
            #    gene_family = []

            #    for paralog in paralogs.Members :
            #        paralog_id = paralog.StableId.lower()

            #        print "\t" + paralog_id
            #        self.genes.add(paralog_id)

            #        gene_family.append(paralog_id)

            #    gene_families.append(gene_family)

            filecontents = ">%s\n%s\n" % (gene_id, gene.getLongestCdsTranscript().ProteinSeq)
            self._write_file_and_manifest(self.basedir, gene_id + ".fa", filecontents)

            if self.stop :
                break

        # write out for testing
        print "\n\nwriting %d genes in %d gene families..." % (len(self.genes), len(gene_families))

        f = open(self.basedir + os.sep + 'genes.txt', 'w')
        g = open(self.basedir + os.sep + 'genefamilies.txt', 'w')

        for fam in gene_families :
            for gene in fam :
                print >> f, gene
                print >> g, gene,
            print >> g, ""

        g.close()
        f.close()

