import commands
import re

from cogent.db.ensembl import Species, Genome, Compara, HostAccount
from cogent.db.ensembl.database import Database

from lib.base import Base
from lib.genefamily import GeneFamily

class EnsemblDbInfo(object) :
    def __init__(self, db_name, low_release, high_release) :
        self.db_name = db_name
        self.latin_name = self._db2latin(self.db_name)
        self.common_name = self._latin2common(self.latin_name)
        self.low_release = low_release
        self.high_release = high_release
        self.release_str = "%d-%d" % (self.low_release, self.high_release)

        #print "%s %s %s" % (self.db_name, self.latin_name, self.common_name)

    def _db2latin(self, db_name) :
        tmp = Species.getSpeciesName(db_name)
        
        if tmp != 'None' :
            return tmp

        return db_name.capitalize().replace('_', ' ')

    def _latin2common(self, latin_name) :
        try :
            return Species.getCommonName(latin_name)
        except :
            pass

        tokens = latin_name.lower().split()

        if len(tokens) == 2 :
            return tokens[0][0].capitalize() + "." + tokens[1]

        raise Exception("Bad latin name: %s" % latin_name)

    def table_str(self, latin_width, common_width, release_width) :
        return self.latin_name.rjust(latin_width) + \
               self.common_name.rjust(common_width) + \
               self.release_str.rjust(release_width)
    
    def __str__(self) :
        return "%s %s %s %s" % (self.latin_name, self.common_name, self.release_str, self.db_name)

class EnsemblInfo(object) :
    def __init__(self, options) :
        self.db_host = options['db-host']
        self.db_port = options['db-port']
        self.db_user = options['db-user']
        self.db_pass = options['db-pass']

        self.verbose = options['verbose']

        self.databases = self._get_databases()

    def _convert_to_range(self, releases) :
        releases.sort()
        return "%d-%d" % (releases[0], releases[-1])

    def _get_databases(self) :
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

        databases = {}

        for dbname in db2rel :
            try :
                databases[dbname] = EnsemblDbInfo(dbname, min(db2rel[dbname]), max(db2rel[dbname]))
                
                # add to pycogent as well
                if Species.getSpeciesName(databases[dbname].latin_name) == 'None' :
                    if self.verbose :
                        print >> sys.stderr, "Info: adding '%s' to pycogent" % databases[dbname].latin_name
                    Species.amendSpecies(databases[dbname].latin_name, databases[dbname].common_name)
                    
                    #print >> sys.stderr, "\t" + Species.getCommonName(databases[dbname].latin_name)
                    #print >> sys.stderr, "\t" + Species.getEnsemblDbPrefix(databases[dbname].latin_name)
            except :
                if self.verbose :
                    print >> sys.stderr, "Info: rejected '%s'" % dbname

        return databases

    def get_latest_release(self, species) :
        try :
            return self.databases.get(Species.getEnsemblDbPrefix(species)).high_release
        except :
            return -1

    def _calc_rjust(self, title, variable) :
        return len(sorted([title] + map(lambda x: getattr(x, variable), self.databases.values()), key=len, reverse=True)[0]) + 2

    def print_species_table(self) :
        l_len = self._calc_rjust("Name", "latin_name")
        c_len = self._calc_rjust("Common name", "common_name")
        r_len = self._calc_rjust("Releases", "release_str")

        print "Name".rjust(l_len) + "Common Name".rjust(c_len) + "Releases".rjust(r_len)
        print "-" * (l_len + c_len + r_len)

        for name in Species.getSpeciesNames() :
            try :
                print self.databases[Species.getEnsemblDbPrefix(name)].table_str(l_len, c_len, r_len)
            except KeyError, ke :
                pass

    def is_valid_species(self, species) :
        try :
            return self.databases.has_key(Species.getEnsemblDbPrefix(species))
        except :
            return False

    def is_valid_release(self, species, release) :
        if not self.is_valid_species(species) :
            return False

        tmp = self.databases[Species.getEnsemblDbPrefix(species)]

        return (release >= tmp.low_release) and (release <= tmp.high_release)

class EnsemblDownloader(Base) :
    def __init__(self, opt) :
        super(EnsemblDownloader, self).__init__(opt)
        
        self.species = opt['species']
        self.release = opt['release']
        self.account = HostAccount(
                            opt['db-host'], 
                            opt['db-user'], 
                            opt['db-pass'], 
                            port=opt['db-port']
                         )
        self.database = opt['database']
        self.genes = set()

    def set_already_downloaded(self, genes) :
        self.genes.update(genes)

    def __iter__(self) :
        return self.genefamilies()

    def genefamilies(self) :
        genome = Genome(self.species, Release=self.release, account=self.account)
        compara = Compara([self.species], Release=self.release, account=self.account)

        self.warn("current version only works with species in ensembl and ensembl-metazoa")

        # DON'T TRY THIS AT HOME!
        #
        # what happens is it searches for compara databases, but unfortunately finds more than one
        # in this situation pycogent just connects to the first one, which is always compara_bacteria
        # so one solution is to dig through all the compara objects internals to provide a connection
        # to the correct database ... obviously not the best solution, but at 6 lines of code definitely 
        # the shortest ;-P
        #
        if self.database == 'ensembl-genomes' :
            from cogent.db.ensembl.host import DbConnection
            from cogent.db.ensembl.name import EnsemblDbName
            import sqlalchemy

            new_db_name = EnsemblDbName(compara.ComparaDb.db_name.Name.replace('bacteria', 'metazoa'))
            compara.ComparaDb._db = DbConnection(account=self.account, db_name=new_db_name)
            compara.ComparaDb._meta = sqlalchemy.MetaData(compara.ComparaDb._db)
        # end of DON'T TRY THIS AT HOME!

        for gene in genome.getGenesMatching(BioType='protein_coding') :
            stableid = gene.StableId.lower()

            # ignore genes that have already been seen as members of
            # gene families
            if stableid in self.genes :
                continue

            self.genes.add(stableid)

            # get cds sequences of any paralogs
            paralogs = compara.getRelatedGenes(StableId=stableid, Relationship='within_species_paralog')


            gf = GeneFamily()

            if paralogs is None :
                gf[stableid] = gene.getLongestCdsTranscript().Cds
            else :
                for paralog in paralogs.Members :
                    paralog_id = paralog.StableId.lower()
                    self.genes.add(paralog_id)
                    gf[paralog_id] = paralog.getLongestCdsTranscript().Cds

            yield gf

