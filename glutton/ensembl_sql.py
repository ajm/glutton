from sqlalchemy import create_engine
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.sql import text

from time import time
from sys import stderr, exit
from collections import defaultdict

import itertools 
import re

from glutton.utils import get_log


DEBUG = False
#DEBUG = True

ensembl_sql_hosts = {
            'ensembl' : {
                    'username' : 'anonymous',
                    'password' : '',
                    'hostname' : 'ensembldb.ensembl.org',
                    'port'     : 5306
                },
            'ensembl-genomes' : {
                    'username' : 'anonymous',
                    'password' : '',
                    'hostname' : 'mysql-eg-publicsql.ebi.ac.uk',
                    'port'     : 4157
                }
        }

def custom_database(hostname, port, username, password) :
    global ensembl_sql_hosts

    ensembl_sql_hosts.clear()

    ensembl_sql_hosts['user'] = {}

    ensembl_sql_hosts['user']['hostname'] = hostname
    ensembl_sql_hosts['user']['port']     = port
    ensembl_sql_hosts['user']['username'] = username
    ensembl_sql_hosts['user']['password'] = password

def get_all_sequences_SQL(genome_db_id, release, nucleotide=False) :
    source_name = "ENSEMBLTRANS" if nucleotide else "ENSEMBLPEP"
 
    if release < 76 :
        return """
    SELECT m.member_id, m2.stable_id, s.sequence
    FROM sequence s 
    JOIN member m 
        USING (sequence_id)
    JOIN member m2 
        ON m2.canonical_member_id=m.member_id 
    WHERE m2.source_name="ENSEMBLGENE"
        AND m.source_name="%s"
        AND m2.genome_db_id=%d""" % (source_name, genome_db_id)
    else :
        return """
    SELECT gm.gene_member_id, gm.stable_id, s.sequence
    FROM sequence s 
    JOIN seq_member sm 
        USING (sequence_id)
    JOIN gene_member gm
        ON sm.seq_member_id=gm.canonical_member_id
    WHERE sm.source_name="%s"
        AND gm.genome_db_id=%d""" % (source_name, genome_db_id)
 
def get_all_homology_SQL(genome_db_id, release) :
    if release < 76 :
        return """
    SELECT h.homology_id, h.description, GROUP_CONCAT(hm.peptide_member_id) as member_list 
        FROM homology h 
    JOIN homology_member hm 
        USING (homology_id) 
    JOIN method_link_species_set m 
        ON h.method_link_species_set_id=m.method_link_species_set_id AND h.description="within_species_paralog" 
    JOIN species_set s 
        USING (species_set_id) 
    WHERE s.genome_db_id=%d 
    GROUP BY h.homology_id""" % genome_db_id
    else :
        return """
    SELECT h.homology_id, h.description, GROUP_CONCAT(hm.gene_member_id) as member_list 
        FROM homology h 
    JOIN homology_member hm 
        USING (homology_id) 
    JOIN method_link_species_set m 
        ON h.method_link_species_set_id=m.method_link_species_set_id AND h.description="within_species_paralog" 
    JOIN species_set s 
        USING (species_set_id) 
    WHERE s.genome_db_id=%d 
    GROUP BY h.homology_id""" % genome_db_id


# need to include '' as a wildcard
def invalid_ensembl_db(name) :
    return name not in ['', 'ensembl', 'metazoa', 'bacteria', 'fungi', 'protists', 'plants']


class SpeciesNotFoundError(Exception) :
    pass

class ReleaseNotFoundError(Exception) :
    pass

class SQLQueryError(Exception) :
    pass

class NoResultsError(Exception) :
    pass


def make_connection(user, password, host, port, db="", echo=False) :
    e = create_engine('mysql://%s:%s@%s:%d/%s' % (user, password, host, port, db), echo=echo)
    try :
        return e.connect()
    
    except SQLAlchemyError, sql :
        raise SQLQueryError(sql.message)

def make_connection_dict(d, db="", echo=False) :
    return make_connection(d['username'], d['password'], d['hostname'], d['port'], db, echo)

# transform a list of integers into a 'range string'
# e.g. [1,2,3,5,6,8] -> "1-3,5-6,8"
def list2rangestr(l) :
    numbers = sorted(l)
    ranges = []

    begin = last = numbers.pop(0)
    
    while numbers :
        current = numbers.pop(0)

        if (current - 1) != last :
            ranges.append((begin, last))
            begin = current 

        last = current
    
    ranges.append((begin, last))

    return ','.join(["%d" % r1 if r1 == r2 else "%d-%d" % (r1, r2) for r1,r2 in ranges])

# transform a 'range string', i.e. a string of the form "X-Y,Z,A-B,"
# into a list of integers in ascending order
# intervals X-Y are assumed to be inclusive
def rangestr2list(rangestr) :
    return sorted(itertools.chain.from_iterable([ 
                    range(*[int(j)+ind for ind,j in enumerate(i.split('-'))]) 
                        if '-' in i else [int(i)] 
                            for i in rangestr.split(',') ]
                ))

# get all versions of the compara database at the db hosts
# listed in ensembl_sql_hosts
def get_compara_versions() :
    global ensembl_sql_hosts
    
    version_table = defaultdict(list)
    db_table = {}

    for hostkey in ensembl_sql_hosts :
        c = make_connection_dict(ensembl_sql_hosts[hostkey])

        result = c.execute(text('show databases'))

        for r in result :
            db_name = r[0]

            is_compara = re.match("ensembl_compara_(\d+)|ensembl_compara_(\w+)_(\d+)_(\d+)", db_name)

            if is_compara :
                compara = is_compara.groups()
                if compara[0] :
                    name = 'ensembl'
                    version = int(compara[0])
                else :
                    name = compara[1]
                    version = int(compara[2])

                # i have no idea what pan_homology is...
                if name != 'pan_homology' :
                    version_table[name].append(version)
                    db_table[(name, version)] = (hostkey, db_name)

    for i in version_table :
        version_table[i] = list2rangestr(version_table[i])

    return version_table, db_table

# get all species listed in a specific database, db_name
# e.g. ('ensembl-genomes', 'ensembl_compara_metazoa_19_75')
def get_compara_species(db_host, db_name) :
    global ensembl_sql_hosts

    species_table = {}
    c = make_connection_dict(ensembl_sql_hosts[db_host], db_name)

    result = c.execute(text("select genome_db_id, name, assembly, genebuild from genome_db"))

    for r in result :
        genome_db_id,name,assembly,genebuild = r
        species_table[name] = (genome_db_id, assembly, genebuild)

    result.close()

    return species_table

def get_all_species_sql(db, suppress) :
    return get_species_versions(db_name=db, suppress=suppress).items()

# return dictionary of species name -> range string
# parameter db_name can be used to limit the enumeration to
# e.g. 'metazoa'
def get_species_versions(db_name="", species=None, human_readable=True, suppress=None) :
    
    if invalid_ensembl_db(db_name) :
        raise BadEnsemblDBNameError("%s not a valid name" % db_name)
    
    version_table, db_table = get_compara_versions()
    species_version_table = defaultdict(list)

#    print db_table
#    print "\n\n"
#    print sorted(db_table, key=lambda x : x[1], reverse=True)
#    exit(0)

    for db in sorted(db_table, key=lambda x : x[1], reverse=True) :
#        print db
        
        if not db[0].startswith(db_name) :
            continue

        species_table = get_compara_species(*db_table[db])

        if species :
            if species in species_table :
                db_name = db[0] # this will speed up searches for specific species
                species_version_table[species].append(db[1])
        else :
            for s in species_table :
                species_version_table[s].append(db[1])

    # suppress releases older than a certain version
    # but don't do this by default
    if suppress :
        keys_to_delete = []
        for i in species_version_table :
            species_version_table[i] = [ x for x in species_version_table[i] if not x < suppress ]

            if len(species_version_table[i]) == 0 :
                keys_to_delete.append(i)

        for i in keys_to_delete :
            del species_version_table[i]

    if human_readable :
        for i in species_version_table :
            species_version_table[i] = list2rangestr(species_version_table[i])

    return dict(species_version_table)

def test_species_listing() :
    version_table, db_table = get_compara_versions()

    for db in db_table :
        try :
            print "trying %s %d -" % db,
            species_table = get_compara_species(*db_table[db])
            print "got %d species" % len(species_table)
        except :
            print "error"

# find out what the version of the compara database is called
# for this species and release number
#
# returns a 3-tuple
#   genome_db_id - internal species identity in this version of ensembl
#   db_host - key into ensembl_sql_hosts specifying connection information
#   db_name - name of database table to read in compara for this species + release
#
# XXX
# 
# some species are in multiple databases, e.g. drosophila melanogaster is in both
# ensembl-main and ensembl-metazoa, we need to make sure that the non-ensembl-main
# one is always found first
# 
def find_database_for_species(species, release, database_name) :
    version_table, db_table = get_compara_versions()

    for db in [ (name,version) for name,version in db_table if version == release ] :
        if not db[0].startswith(database_name) :
            continue
        
        species_table = get_compara_species(*db_table[db])

        if species in species_table :
            db_host, db_name = db_table[db]
            return species_table[species][0], db_host, db_name

    raise ReleaseNotFoundError("could not find release %d for %s" % (release, species))

def perform_query(connection, query) :
    try :
        return connection.execute(text(query)).fetchall()

    except SQLAlchemyError, sae :
        raise SQLQueryError(str(sae))

def group_into_families(peptides, homologies) :
    seen = set()
    all_families = []

    for pep in peptides :
        if pep in seen :
            continue

        fam = []

        if homologies[pep] :
            for pepid in homologies[pep] :
                fam.append(peptides[pepid])
                seen.add(pepid)
        else :
            fam.append(peptides[pep])

        all_families.append(fam)

    return all_families

def get_canonical_sequences(connection, species, release, nucleotide, genome_db_id) :
    # get a complete listing of canonical peptides for species
    raw_results = perform_query(connection, get_all_sequences_SQL(genome_db_id, release, nucleotide))

    id2peptide = dict([ (r[0], (r[1], r[2])) for r in raw_results ])
    
    if DEBUG :
        with open("%s_%d_%s_raw.txt" % (species, release, "nuc" if nucleotide else "pep"), 'w') as f :
            for r in raw_results :
                print >> f, r

    return id2peptide

def get_homology_information(connection, species, release, genome_db_id) :
    # get homology information about peptide sequences 
    raw_results = perform_query(connection, get_all_homology_SQL(genome_db_id, release))

    homologies = defaultdict(set)

    for r in raw_results :
        x,y = [ int(i) for i in r[2].split(',') ]
        homologies[x].add(x)
        homologies[y].add(y)
        homologies[x].add(y)
        homologies[y].add(x)

    if DEBUG :
        with open("%s_%d_homology_raw.txt" % (species, release), 'w') as f :
            for r in raw_results :
                print >> f, r

    return homologies

def get_latest_release_sql(species, database_name) :
    try :
        return max(get_species_versions(species=species, db_name=database_name, human_readable=False)[species])

    except KeyError :
        raise SpeciesNotFoundError("%s not found in any ensembl database" % species)

# download a listing of the canonical peptides and homology relations
# for within_species_paralogs
# e.g. tribolium_castineum, 75
def download_database_sql(species, release=None, database_name='ensembl', nucleotide=False) :
    global ensembl_sql_hosts
    global DEBUG

    # if release is not specified we need to find the latest release
    if not release :
        if DEBUG :
            print >> stderr, "release not specified, finding latest version..."

        release = get_latest_release(species)

    if DEBUG :
        print >> stderr, "downloading %s/%d -" % (species, release),

    genome_db_id, db_host, db_name = find_database_for_species(species, release, database_name)

    #print db_host, db_name

    if DEBUG :
        print >> stderr, "genome_db_id =", genome_db_id 

    connection = make_connection_dict(ensembl_sql_hosts[db_host], db_name)

    # a user would want to use the ensembl-genomes release number which differs
    # by 53 compared to the ensembl main version number, so increase it here to 
    # allow the correct query for the appropriate db schema 
    if db_host == 'ensembl-genomes' :
        release += 53

    id2peptide = get_canonical_sequences(connection, species, release, nucleotide, genome_db_id)
    homologies = get_homology_information(connection, species, release, genome_db_id)

    if not id2peptide or not homologies :
        raise NoResultsError("no results returned - %s" % "nucleotide sequences are unavailable in ensembl-compara outside of most of the species in ensembl-main, use 'biomart' instead of 'sql'" if nucleotide else "maybe try again later?")

    return group_into_families(id2peptide, homologies)

def test_download_database(species) :
    for release in rangestr2list('15-23') :
        download_database(species, release)


#if __name__ == '__main__' :
#    test_species_listing()
#    exit(0)

#    e = EnsemblDownloader()

    #for i in e.get_all_species(db='protists') :
    #    print i

#    families = e.download('homo_sapiens', 77, nucleotide=False)


#version_table, db_table = get_compara_versions()
#print version_table

#test_download_database('tribolium_castaneum')
#exit(0)


# get all versions of ensembl + release numbers
#
# version_table is dict for printing (human readable) release numbers
#   key can be 'ensembl', 'metazoa', 'fungi', 'plants', 'bacteria'
#   value is range string in the form 'X-Y,A,B-C'
#
# db_table maps a 2-tuple of (key, release) to the host + database name
#   key can be 'ensembl', 'metazoa', 'fungi', 'plants', 'bacteria'
#version_table, db_table = get_compara_versions()
#print version_table

# get all species in version 76 of ensembl metazoa
# 
# species_table is a dict 
#   key = species name
#   value = (genome_db_id, assembly, genebuild)
#species_table = get_compara_species(*db_table[('metazoa', 76)])
#print species_table

# get all species + release range strings from specified database - all releases
# returns dict 
#   key = species name
#   value = range string of releases
#
# XXX may not that useful because there is so much crap in the early releases...
#       but can do get_species_versions('metazoa')['tribolium_castaneum'] to get all releases
#species_table = get_species_versions('metazoa')
#print species_table

# download species database for a given release
# XXX still working on...
#download_database('tribolium_castaneum', 76)

