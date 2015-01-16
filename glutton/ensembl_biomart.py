import urllib
import urllib2

from xml.dom.minidom import parse
from sys import exit, stderr
from collections import defaultdict

from Bio import SeqIO

from glutton.utils import get_log


# 
# biomart requires slightly different values for virtualSchemaName
# depending on the installation
#
# biomart.org works with it set to "default", but ensembl.org or
# metazoa.ensembl.org requires something like "metazoa_mart_24"
# or else it does not return anything...
#



TIMEOUT = 30
#URL = 'http://www.biomart.org/biomart/martservice'
#URL = 'http://metazoa.ensembl.org/biomart/martservice'

def get_URL(database_name) :
    return "http://%s.ensembl.org/biomart/martservice" % (database_name if database_name != 'ensembl' else 'www')

def get_marts(database_name) :
    global TIMEOUT

    if database_name :
        marts = { database_name : None }
    else :
        marts = {   
                'ensembl'   : None,
                'metazoa'   : None,
                'plants'    : None,
                'protists'  : None,
                'fungi'     : None,
                'bacteria'  : None  }

    log = get_log()
    query = get_URL(database_name) + '?type=registry'
    
    try :
        dom = parse(urllib2.urlopen(query, timeout=TIMEOUT))
    
    except urllib2.URLError, urle :
        log.fatal("biomart registry query timed out")
        exit(1)

    for ele in dom.getElementsByTagName("MartURLLocation") :
        schema_name =   ele.attributes["database"].value
        db_display =    ele.attributes["displayName"].value
        mart_name =     ele.attributes["name"].value

        log.debug("%s %s %s" % (schema_name, db_display, mart_name))

        # this messes things up
        if mart_name == 'ENSEMBL_MART_PLANT' :
            continue

        for k in marts :
            if schema_name.startswith(k + '_mart_') :
                marts[k] = (mart_name, schema_name, int(schema_name.split('_')[-1]))
                break


    return marts

def get_latest_release_biomart(species, database_name) :
    mart_name, schema_name, release, db_name, species_desc = get_all_species(get_marts(database_name), database_name)[species]
    
    return release

def get_all_species(marts, database_name) :
    global TIMEOUT

    query = get_URL(database_name) + '?type=datasets&mart=%s'
    species = {}

    log = get_log()

    for k in marts :
        if not marts[k] :
            continue

        mart_name,schema_name,release = marts[k]

        try :
            f = urllib2.urlopen(query % mart_name, timeout=TIMEOUT)

        except urllib2.URLError, ue :
            log.fatal("biomart dataset listing timed out")
            exit(1)

        for line in f : 
            line = line.strip()
            if not line :
                continue

            data = line.split('\t')

            db_name = data[1]
            description = data[2]

            if ' genes ' in description :
                species_name = description.split(' genes ')[0].lower().replace(' ', '_')
            else :
                species_name = db_name.split('_')[0]

            species[species_name] = (mart_name, schema_name, release, db_name, description)

        f.close()

    return species

def get_all_species_biomart(database_name, suppress) :
    return [ (i, str(j[2])) for i,j in get_all_species(get_marts(database_name), database_name).items() ]

# check release is the most recent one and raise 
# a ReleaseNotFoundException if it is not the most recent one
def download_database_biomart(species, release, database_name='ensembl', nucleotide=True) :
    
    mart_name, schema_name, mart_release, table_name, species_desc = get_all_species(get_marts(database_name), database_name)[species]

    if mart_release != release :
        get_log().fatal("requested release (%d) does not match current biomart release (%d)" % (release, mart_release))
        exit(1)

    seq = get_sequences(species, database_name, schema_name, table_name, nucleotide)
    homo = get_homology_info(species, database_name, schema_name, table_name, nucleotide)

    return group_into_families(seq, homo)

def group_into_families(peptides, homologies) :
    seen = set()
    all_families = []

    for pep in peptides :
        if pep in seen :
            continue

        fam = []

        if homologies[pep] :
            for pepid in homologies[pep] :
                fam.append((pepid, peptides[pepid]))
                seen.add(pepid)
        else :
            fam.append((pep, peptides[pep]))

        all_families.append(fam)

    return all_families

def get_sequences(species, database_name, schema_name, table_name, nucleotide) :
    global TIMEOUT
    
    log = get_log()
    log.debug("getting sequences...")

    payload_sequences = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="%s" formatter="FASTA" header="0" uniqueRows="0" count="" datasetConfigVersion="0.7">
    <Dataset name="%s" interface="default">
        <Attribute name="ensembl_gene_id" />
        <Attribute name="%s" />
    </Dataset>
</Query>"""

    # e.g. 'dipodomys_ordii': ('dordii_gene_ensembl', 'Dipodomys ordii genes (dipOrd1)')
    #
    seq_type = 'cdna' if nucleotide else 'peptide'
    query = payload_sequences % (schema_name if database_name != 'ensembl' else 'default', table_name, seq_type)
    params = urllib.urlencode({ 'query' : query })

    sequences = {}

    try :
        f = urllib2.urlopen(get_URL(database_name), params)

    except urllib2.URLError, ue :
        log.fatal("biomart sequence query error (%s)" % str(ue))
        exit(1)


    stderr.write("\r[downloading %s] got %d sequences " % ("cDNA" if nucleotide else "protein", len(sequences)))

    for s in SeqIO.parse(f, 'fasta') :
        sequences[s.id] = str(s.seq)
        stderr.write("\r[downloading %s] got %d sequences " % ("cDNA" if nucleotide else "protein", len(sequences)))

    f.close()

    stderr.write("\r[downloading %s] got %d sequences\n" % ("cDNA" if nucleotide else "protein", len(sequences)))

    return sequences

def get_homology_info(species, database_name, schema_name, table_name, nucleotide) :
    global TIMEOUT

    log = get_log()
    log.debug("getting homology...")

    payload_homology = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="%s" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.7">
    <Dataset name="%s" interface="default">
        <Attribute name="ensembl_gene_id" />
        <Attribute name="%s" />
    </Dataset>
</Query>"""

    # is this always correct?
    short_species_name = table_name.split('_')[0]

    # this does not seem like a great idea, but until it breaks
    # it will do...
    if '_eg_' in table_name :
        paralog_name = short_species_name + '_eg_paralog_gene'
    else :
        paralog_name = short_species_name + '_paralog_ensembl_gene'

    query2 = payload_homology % (schema_name, table_name, paralog_name)
    params2 = urllib.urlencode({ 'query' : query2 })

    # make api call
    try :
        f = urllib2.urlopen(get_URL(database_name), params2)

    except urllib2.URLError, ue :
        log.fatal("biomart homology query error (%s)" % str(ue))
        exit(1)

    # parse
    homologies = defaultdict(set)

    stderr.write("\r[downloading homology] got 0 records ")
    count = 0

    for line in f :
        try :
            x,y = line.strip().split()

        except ValueError :
            continue

        count += 1
        stderr.write("\r[downloading homology] got %d records " % count)

        homologies[x].add(x)
        homologies[y].add(y)
        homologies[x].add(y)
        homologies[y].add(x)

    f.close()

    stderr.write("\r[downloading homology] got %d records\n" % count)
    
    return homologies


if __name__ == '__main__' :
    from glutton.utils import setup_logging

    setup_logging()

    #print get_all_species(get_marts())

    download_database_biomart('homo_sapiens')
    #download_database_biomart('tribolium_castaneum')

