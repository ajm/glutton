from sys import exit

from glutton.utils import get_log
from glutton.ensembl_sql import ensembl_sql_hosts, get_latest_release_sql, find_database_for_species


def get_missing_info(species, release) :
    global ensembl_sql_hosts

    if not release :
        release = get_latest_release_sql(species)

    genome_db_id, db_host, db_name = find_database_for_species(species, release)

    project_name = 'ensembl'
    for i in ('metazoa', 'bacteria', 'fungi', 'protists', 'plants') :
        if i in db_name :   
            project_name = i
            break

    get_log().info("%s is found in %s" % (species, project_name))

    return release, project_name, ensembl_sql_hosts[db_host]

def get_all_species_pycogent(db, suppress) :
    log = get_log()

    try :
        from cogent.db.ensembl import Species

    except ImportError :
        log.fatal("importing pycogent failed, exiting...")
        exit(1)

    log.warning("pycogent cannot differentiate between ensembl projects, listing everything...")

    # annoying that i cannot see another way to programmatically get
    # a list of supported species
    return [ (i.split()[-1], 'unknown') for i in str(Species).split('\n')[3:-1] ]

# db_name has to be one of (metazoa, protists, plants, bacteria) or empty string for ensembl main
def download_database_pycogent(species, release, nucleotide=False) :
    log = get_log()

    try :
        import cogent
        from cogent.db.ensembl import Species, Genome, Compara, HostAccount
        from cogent.db.ensembl.database import Database

    except ImportError :
        log.fatal("importing pycogent failed, exiting...")
        exit(1)

    if cogent.version_info != (1,5,3) :
        log.warning("only tested with pycogent version 1.5.3 (you are running %s)" % cogent.version)


    release, db_name, db_details = get_missing_info(species, release)

    account = HostAccount(
                db_details['hostname'],
                db_details['username'],
                db_details['password'],
                port=db_details['port'])

    if Species.getSpeciesName(species) == 'None' : # this is not an error, it returns the string "None"
        log.warning("%s not found in pycogent, attempting to add it manually" % species)
        Species.amendSpecies(species.capitalize().replace('_', ' '), species)

    genome = Genome(species, Release=release, account=account)
    compara = Compara([species], Release=release, account=account)



    # DON'T TRY THIS AT HOME!
    #
    # what happens is it searches for compara databases, but unfortunately finds more than one
    # in this situation pycogent just connects to the first one, which is always compara_bacteria
    # so one solution is to dig through all the compara objects internals to provide a connection
    # to the correct database ... obviously not the best solution, but at 6 lines of code definitely
    # the shortest ;-P
    #
    if db_name not in ('main', 'bacteria') :
        log.warning("accessing compara from pycogent with species outside of ensembl-main and ensembl-bacteria is problematic, attempting to patch...")

        from cogent.db.ensembl.host import DbConnection
        from cogent.db.ensembl.name import EnsemblDbName
        import sqlalchemy

        new_db_name = EnsemblDbName(compara.ComparaDb.db_name.Name.replace('bacteria', db_name))
        compara.ComparaDb._db = DbConnection(account=account, db_name=new_db_name)
        compara.ComparaDb._meta = sqlalchemy.MetaData(compara.ComparaDb._db)
    # end of DON'T TRY THIS AT HOME!



    genes = set()
    families = []

    for gene in genome.getGenesMatching(BioType='protein_coding') :
        stableid = gene.StableId

        # ignore genes that have already been seen as members of other gene families
        if stableid in genes :
            continue

        genes.add(stableid)

        paralogs = compara.getRelatedGenes(StableId=stableid, Relationship='within_species_paralog')
        
        current = []
        
        if paralogs is None :
            log.info("downloading %s" % stableid)

            current.append((stableid, str(gene.CanonicalTranscript.Cds) if nucleotide else str(gene.CanonicalTranscript.ProteinSeq)))

        else :
            for paralog in paralogs.Members :
                paralogid = paralog.StableId

                log.info("downloading %s" % paralogid)

                genes.add(paralogid)

                try :
                    current.append((paralogid, str(paralog.CanonicalTranscript.Cds) if nucleotide else str(paralog.CanonicalTranscript.ProteinSeq)))
                
                except AttributeError :
                    log.fatal("pycogent did not find a canonical transcript for %s" % paralogid)
                    exit(1)

        #print ','.join([ i for i,j in current ])
        families.append(current)

        log.info("")        

    return families

