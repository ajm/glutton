import sys
import commands
import re
import os

from cogent.db.ensembl import Species

class EnsemblCache(object) :
    def __init__(self, workingdir, tmpdir) :
        self.stop = False
        self.workingdir = workingdir
        self.tmpdir = tmpdir

        self._check_directory(self.workingdir, create=True)
        self._check_directory(self.tmpdir, create=True)

    def stop(self) :
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

    def _check_directory(self, dirname, create=False) :
        if not os.path.exists(dirname) :
            if create :
                print >> sys.stderr, "Info: '%s' does not exist, creating..." % dirname
                try :
                    os.makedirs(dirname)
                except OSError, ose :
                    print >> sys.stderr, str(ose)
                    return False
            else :
                print >> sys.stderr, "Error: '%s' does not exist." % dirname
                return False

        elif not os.path.isdir(dirname) :
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

    def build_cache(self, species, release, resume=True) :
        basedir = self._check_cache_directory(species, release)

        print basedir
