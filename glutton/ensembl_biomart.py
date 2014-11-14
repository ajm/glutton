import urllib
import urllib2

from xml.dom.minidom import parse
from sys import exit, stderr

from Bio import SeqIO

from glutton.utils import get_log


TIMEOUT = 30
URL = 'http://www.biomart.org/biomart/martservice'

def get_marts() :
    global TIMEOUT, URL

    marts = {   'ensembl'   : None,
                'metazoa'   : None,
                'plants'    : None,
                'protists'  : None,
                'fungi'     : None,
                'bacteria'  : None  }

    log = get_log()
    query = URL + '?type=registry'
    
    try :
        dom = parse(urllib2.urlopen(query, timeout=TIMEOUT))
    
    except urllib2.URLError, urle :
        log.fatal("biomart registry query timed out")
        exit(1)

    for ele in dom.getElementsByTagName("MartURLLocation") :
        db_name =       ele.attributes["database"].value
        db_display =    ele.attributes["displayName"].value
        mart_name =     ele.attributes["name"].value

        for k in marts :
            if db_name.startswith(k + '_mart_') :
                marts[k] = mart_name
                break


    return marts

def guess_latest_release(mart_name) :
    return int(mart_name.split('_')[-1])

def get_latest_release_biomart(species) :
    mart_name, db_name, species_desc = get_all_species(get_marts())[species]
    
    return guess_latest_release(mart_name)

def get_all_species(marts) :
    global TIMEOUT, URL

    query = URL + '?type=datasets&mart=%s'
    species = {}

    log = get_log()

    for k in marts :
        #print k

        if not marts[k] :
            continue

        try :
            f = urllib2.urlopen(query % marts[k], timeout=TIMEOUT)

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

            species[species_name] = (marts[k], db_name, description)

            #print "\t", species_name

        f.close()

    return species

def get_all_species_biomart(db, suppress) :
    marts = dict([ (i,j) for i,j in get_marts().items() if i.startswith(db) ])

    if len(marts) != 1 :
        get_log().fatal('internal error, multiple marts found for database (%s)' % db)
        exit(1)

    return [ (i, 'latest only') for i in get_all_species(marts).keys() ]

# check release is the most recent one and raise 
# a ReleaseNotFoundException if it is not the most recent one
def download_database_biomart(species, release, nucleotide=True) :
    # e.g. 'dipodomys_ordii': ('dordii_gene_ensembl', 'Dipodomys ordii genes (dipOrd1)')
    #
    mart_name, db_name, species_desc = get_all_species(get_marts())[species]

    biomart_release = guess_latest_release(mart_name)
    if biomart_release != release :
        get_log().fatal("requested release (%d) does not match current biomart release (%d)" % (release, biomart_release))
        exit(1)

    seq = get_sequences(species, db_name, nucleotide)
    homo = get_homology_info(species, db_name, nucleotide)

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

def get_sequences(species, db_name, nucleotide) :
    global TIMEOUT, URL    
    
    get_log().debug("getting sequences...")

    payload_sequences = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="FASTA" header="0" uniqueRows="0" count="" datasetConfigVersion="0.7">
    <Dataset name="%s" interface="default">
        <Attribute name="ensembl_gene_id" />
        <Attribute name="%s" />
    </Dataset>
</Query>"""

    # e.g. 'dipodomys_ordii': ('dordii_gene_ensembl', 'Dipodomys ordii genes (dipOrd1)')
    #
    seq_type = 'cdna' if nucleotide else 'peptide'
    query = payload_sequences % (db_name, seq_type)
    params = urllib.urlencode({ 'query' : query })

    sequences = {}

    try :
        f = urllib2.urlopen(URL, params)

    except urllib2.URLError, ue :
        log.fatal("biomart sequence query error (%s)" % str(ue))
        exit(1)


    count = 0

    for s in SeqIO.parse(f, 'fasta') :
        sequences[s.id] = s.seq
        #print s.id, s.seq

        count += 1

        if (count % 100) == 0 :
            print "downloaded %d genes..." % count

    f.close()

    return sequences

def get_homology_info(species, db_name, nucleotide) :
    global TIMEOUT, URL

    get_log().debug("getting homology...")

    payload_homology = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.7">
    <Dataset name="%s" interface="default">
        <Attribute name="ensembl_gene_id" />
        <Attribute name="%s" />
    </Dataset>
</Query>"""

    # is this always correct?
    short_species_name = db_name.split('_')[0]

    # this does not seem like a great idea, but until it breaks
    # it will do...
    if '_eg_' in db_name :
        paralog_name = short_species_name + '_eg_paralog_gene'
    else :
        paralog_name = short_species_name + '_paralog_ensembl_gene'

    query2 = payload_homology % (db_name, paralog_name)
    params2 = urllib.urlencode({ 'query' : query2 })

    # make api call
    try :
        f = urllib2.urlopen(URL, params2)

    except urllib2.URLError, ue :
        log.fatal("biomart homology query error (%s)" % str(ue))
        exit(1)

    # parse
    homologies = defaultdict(set)

    for line in f :
        try :
            x,y = line.strip().split()

        except ValueError :
            continue

        homologies[x].add(x)
        homologies[y].add(y)
        homologies[x].add(y)
        homologies[y].add(x)

    f.close()
    
    return homologies


if __name__ == '__main__' :
    from glutton.utils import setup_logging

    setup_logging()

    #print get_all_species(get_marts())

    #download_database_biomart('homo_sapiens')
    download_database_biomart('tribolium_castaneum')

