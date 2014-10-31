import sys
import json
import threading
import collections

from os.path import isfile, join

from glutton.utils import get_log, md5, check_dir


PARAM_FILE  = 'parameters.json'
CONTIG_FILE = 'contigs.json'
BLAST_FILE  = 'blastx.json'
PAGAN_FILE  = 'pagan.json'

QUERY_ID = 'query'


def do_locking(fn) :
    def thread_safe(*args) :
        args[0].lock.acquire()
        ret = fn(*args)
        args[0].lock.release()

        return ret

    return thread_safe

class GluttonException(Exception) :
    pass

class GluttonInformation(object) :
    def __init__(self, db, contigs_fname, alignments_dir) :
        self.db = db
        self.contigs_fname = contigs_fname
        self.directory = alignments_dir

        check_dir(self.directory)

        self.log = get_log()
        self.lock = threading.Lock()

        # the alignment procedure can take a long time, so everything needs to be 
        # restartable, in addition - if we restart it then we need to be sure that 
        # the parameters used are the same, i.e.: same reference database and maybe
        # more
        #
        # blastx results are tricky because another object is doing that, so we need
        # to query that object to get the information
        #
        # pagan results are easier because that is what this object is doing...
        #
        self.params = {}
        self.contig_query_map = {}          # contig id -> query id
        self.query_gene_map = {}            # query id -> gene id
        self.genefamily_filename_map = {}   # gene family id -> filename

        self.read_progress_files()

        self.check_params()

    # files used for recording the progress are just dirty globals
    # at the moment these properties can be in place of a decent solution
    # for now...
    @property
    def parameter_filename(self) :
        global PARAM_FILE
        return join(self.directory, PARAM_FILE)
   
    @property
    def contig_filename(self) :
        global CONTIG_FILE
        return join(self.directory, CONTIG_FILE)

    @property
    def blastx_filename(self) :
        global BLAST_FILE
        return join(self.directory, BLAST_FILE)

    @property
    def pagan_filename(self) :
        global PAGAN_FILE
        return join(self.directory, PAGAN_FILE)

    def flush(self) :
        self.log.info("flushing data to disk...")
        self.write_progress_files()
        self.log.info("done")

    # functions and utilities to read and write the different progress files
    #
    def _load(self, fname) :
        if isfile(fname) :
            self.log.info("found progress file %s ..." % fname)
            return json.loads(open(fname).read())

        return {}

    def _dump(self, fname, data) :
        if data :
            open(fname, 'w').write(json.dumps(data))

    def read_progress_files(self) :
        self.params                     = self._load(self.parameter_filename)
        self.contig_query_map           = self._load(self.contig_filename)
        self.query_gene_map             = self._load(self.blastx_filename)
        self.genefamily_filename_map    = self._load(self.pagan_filename)

        if self.query_gene_map :
            self.log.info("read %d blast results" % len(self.query_gene_map))

        if self.genefamily_filename_map :
            self.log.info("read %d pagan results" % len(self.genefamily_filename_map))

    @do_locking
    def write_progress_files(self) :
        self._dump(self.parameter_filename, self.params)
        self._dump(self.contig_filename,    self.contig_query_map)
        self._dump(self.blastx_filename,    self.query_gene_map)
        self._dump(self.pagan_filename,     self.genefamily_filename_map)

    # related to database parameters
    #
    def get_params(self) :
        p = {}

        p['db_species']  = self.db.species
        p['db_release']  = self.db.release
        p['db_filename'] = self.db.filename
        p['db_checksum'] = self.db.checksum

        p['contig_filename'] = self.contigs_fname
        p['contig_checksum'] = md5(self.contigs_fname)

        return p

    @do_locking
    def check_params(self) :
        db_params = self.get_params()

        if not self.params :
            self.params = db_params
            return
        
        if self._not_same_db(db_params) :
            self.log.fatal("alignment has been resumed with a different reference!")
            self.log.fatal("\toriginal = %s" % self._param_str(self.params))
            self.log.fatal("\tcurrent  = %s" % self._param_str(db_params))
            sys.exit(1)

    def _param_str(self, p) :
        return "%s/%d %s md5=%s" % (p['db_species'], p['db_release'], p['contig_filename'], p['contig_checksum'])

    def _not_same_db(self, par) :
        for key in self.params :
            if self.params[key] != par[key] :
                return True
        return False

    # contig to query ids are only put/get
    #   initialised?/put/get
    def initialised_contig2query(self) :
        return self.contig2query == {}

    @do_locking
    def put_contig2query(self, contig_id, query_id) :
        self.contig_query_map[contig_id] = query_id

    @do_locking
    def get_contig2query(self, contig_id) :
        global QUERY_ID

        try :
            return self.contig_query_map[contig_id]

        except KeyError :
            pass
        
        # well... this makes we queasy...
        #   if there is no attribute in this class called something, then
        #   just create it and initialise it to a sensible value
        if not hasattr(self, 'query_id_counter') :
            self.query_id_counter = 1 + max([0] + [ int(i[len(QUERY_ID):]) for i in self.contig_query_map.values() ])

        new_query_id = "%s%d" % (QUERY_ID, self.query_id_counter)
        self.query_id_counter += 1

        self.contig_query_map[contig_id] = new_query_id

        return new_query_id

    # query id to gene id
    #   get/update/in
    def get_query2gene(self, query_id) :
        return self.query_gene_map[query_id]

    def iter_query2gene(self) :
        for i in self.query_gene_map :
            if self.query_gene_map[i] != 'FAIL' :
                yield i

    @do_locking
    def update_query2gene(self, new_dict) :
        self.query_gene_map.update(new_dict)

    def in_query2gene(self, query_id) :
        return query_id in self.query_gene_map

    # genefamily id to filename or FAIL
    #   put/get/fail/in
    @do_locking
    def put_genefamily2filename(self, genefamily_id, filename='FAIL') :
        self.genefamily_filename_map[genefamily_id] = filename

    def get_genefamily2filename(self, genefamily_id) :
        return self.genefamily_filename_map[genefamily_id]

    def in_genefamily2filename(self, genefamily_id) :
        return genefamily_id in self.genefamily_filename_map

    def len_genefamily2filename(self) :
        return len(self.genefamily_filename_map)

    def build_genefamily2contigs(self) :
        genefamily_contig_map = collections.defaultdict(list)

        for i in self.query_gene_map :
            if self.query_gene_map[i] == 'FAIL' :
                continue

            genefamily_contig_map[self.db.get_familyid_from_geneid(self.query_gene_map[i])].append(i)
        
        return genefamily_contig_map

    # aggregate actions
    # 
    def pending_queries(self) :
        return [ i for i in self.contig_query_map.values() if i not in self.query_gene_map ]

    def alignments_complete(self) :
        genefamily_contig_map = self.build_genefamily2contigs()

        for i in genefamily_contig_map :
            if i not in self.genefamily_filename_map :
                return False

        return True

if __name__ == '__main__' :
    
    from glutton.db import GluttonDB

    db = GluttonDB('tc23.glt')
    gi = GluttonInformation(db, 'queries_test.fasta', './alignment_test')

