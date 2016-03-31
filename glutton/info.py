import json
import threading
import collections

from sys import stderr, exit
from os.path import isfile, join, abspath

from glutton.db import GluttonDB
from glutton.utils import get_log, md5, check_dir
from glutton.table import pretty_print_table


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

GluttonSample = collections.namedtuple('GluttonSample', ['id','species','contigs','checksum','bam'])

class GluttonJSON(object) :
    def __init__(self) :
        self.log = get_log()

    def load(self, fname) :
        if isfile(fname) :
            self.log.info("found progress file %s ..." % fname)
            return json.loads(open(fname).read())

        return {}

    def dump(self, fname, data) :
        if data :
            open(fname, 'w').write(json.dumps(data, sort_keys=True, indent=4, separators=(',', ': ')))

class GluttonParameters(GluttonJSON) :
    def __init__(self, project_dir, create=False) :
        self.directory = project_dir
        self.create = create
        check_dir(self.directory, create=self.create)

        self.log = get_log()

        self.params = self.load(self.parameter_filename)

        if not self.params :
            self.params = { 'db_species'    : None,
                            'db_release'    : None,
                            'db_filename'   : None,
                            'db_checksum'   : None,
                            'samples'       : {} }

    @property
    def parameter_filename(self) :
        global PARAM_FILE
        return join(self.directory, PARAM_FILE)

    def get_sample_ids(self) :
        return self.params['samples'].keys()
    
    def add(self, contigfile, sampleid, species, bamfile=None, assembler=None, copy=False) :
        self.log.info("adding %s (%s, %s, %s, %s)" % (sampleid, species, contigfile, bamfile, assembler))
        contigfile = self.copy(contigfile) if copy else abspath(contigfile)
        bamfile = self.copy(bamfile) if copy else abspath(bamfile) if bamfile else bamfile 

        self.params['samples'][sampleid] = { 'contigs'   : contigfile,
                                             'checksum'  : md5(contigfile),
                                             'species'   : species,
                                             'bam'       : bamfile,
                                             'assembler' : assembler }

    def copy(self, src) :
        data_dir = join(self.directory, 'data')
        check_dir(data_dir, create=True)
        dst = join(data_dir, basename(src))
        self.log("cp %s %s" % (src, dst))
        shutil.copy(src, dst)
        return dst

    def remove(self, sampleid) :
        self.log.info("removing %s" % sampleid)
        del self.params['samples'][sampleid]

    def list(self) :
        if not self.params['db_filename'] :
            print "No reference used..."
        else :
            pretty_print_table(('SPECIES','RELEASE','FILENAME','CHECKSUM'), \
                [(self.params['db_species'], self.params['db_release'], self.params['db_filename'], self.params['db_checksum'])])
        
        if self.count() == 0 :
            print "No samples found..."
        else :
            pretty_print_table(('SAMPLE','SPECIES','CONTIGS','MD5','BAM'), \
                [ (k,v['species'], v['contigs'], v['checksum'], v['bam']) for k,v in self.params['samples'].iteritems() ])

    def count(self) :
        return len(self.params['samples'])

    def flush(self) :
        self.dump(self.parameter_filename, self.params)

    def get_species(self, id) :
        return self.params['samples'][id]['species']

    def get_checksum(self, id) :
        return self.params['samples'][id]['checksum']

    def get_contigs(self, id) :
        return self.params['samples'][id]['contigs']

    def get_bam(self, id) :
        return self.params['samples'][id]['bam']

    def all(self) :
        for k,v in self.params['samples'].iteritems() :
            yield GluttonSample(id=k, contigs=v['contigs'], species=v['species'], bam=v['bam'], checksum=v['checksum'])

class GluttonInformation(GluttonJSON) :
    def __init__(self, alignments_dir, parameters, db) :
        self.directory = alignments_dir
        self.params = parameters
        self.db = db

        check_dir(self.directory)

        self.log = get_log()
        self.lock = threading.RLock() # a single function requires this be an RLock over a Lock

        # the alignment procedure can take a long time, so everything needs to be 
        # restartable, in addition - if we restart it then we need to be sure that 
        # the parameters used are the same, i.e.: same reference database etc etc

        self.contig_query_map = {}          # file id -> contig id -> query id (file id is provided by the user, called a 'label')
        self.query_gene_map = {}            # query id -> (gene id, +/-) or None
        self.genefamily_filename_map = {}   # gene family id -> filename

        self.read_progress_files()

    @property
    def contig_filename(self) :
        global CONTIG_FILE
        return join(self.directory, CONTIG_FILE)

    @property
    def blast_filename(self) :
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

    def read_progress_files(self) :
        self.contig_query_map           = self.load(self.contig_filename)
        self.query_gene_map             = self.load(self.blast_filename)
        self.genefamily_filename_map    = self.load(self.pagan_filename)

        if self.contig_query_map :
            self.log.info("read %d contig to query id mappings" % sum([ len(self.contig_query_map[label]) for label in self.contig_query_map ]))

        if self.query_gene_map :
            self.log.info("read %d blast results" % len(self.query_gene_map))

        if self.genefamily_filename_map :
            self.log.info("read %d pagan results" % len(self.genefamily_filename_map))

    @do_locking
    def write_progress_files(self) :
        self.dump(self.contig_filename,    self.contig_query_map)
        self.dump(self.blast_filename,     self.query_gene_map)
        self.dump(self.pagan_filename,     self.genefamily_filename_map)

    def _set_id_counter(self) :
        tmp = [0]
        for label in self.contig_query_map :
            for i in self.contig_query_map[label].values() :
                tmp.append(int(i[len(QUERY_ID):]))

        self.query_id_counter = 1 + max(tmp)

    # contig to query ids are only get
    @do_locking
    def get_query_from_contig(self, label, contig_id) :
        global QUERY_ID

        try :
            return self.contig_query_map[label][contig_id]

        except KeyError :
            pass
        
        # well... this makes we queasy...
        #   if there is no attribute in this class called something, then
        #   just create it and initialise it to a sensible value
        if not hasattr(self, 'query_id_counter') :
            self._set_id_counter()

        new_query_id = "%s%d" % (QUERY_ID, self.query_id_counter)
        self.query_id_counter += 1

        if label not in self.contig_query_map :
            self.contig_query_map[label] = {}

        self.contig_query_map[label][contig_id] = new_query_id

        return new_query_id

    # query id to gene id
    #   update
    @do_locking
    def update_query_gene_mapping(self, new_dict) :
        self.query_gene_map.update(new_dict)

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

    # aggregate actions
    #
    @do_locking
    def build_genefamily2contigs(self) :
        genefamily_contig_map = collections.defaultdict(list)

        for i in self.query_gene_map :
            if self.query_gene_map[i] == None :
                continue

            geneid,strand = self.query_gene_map[i]
            genefamily_contig_map[self.db.get_familyid_from_geneid(geneid)].append((i,strand))
        
        return genefamily_contig_map

    @do_locking
    def pending_queries(self) :
        tmp = []

        for label in self.contig_query_map :
            for i in self.contig_query_map[label].values() :
                if i not in self.query_gene_map :
                    tmp.append(i)

        return tmp

    @do_locking
    def num_alignments_not_done(self) :
        genefamily_contig_map = self.build_genefamily2contigs()
        not_done = 0
        failures = 0

        for i in genefamily_contig_map :
            if i not in self.genefamily_filename_map :
                not_done += 1
                continue

            if self.genefamily_filename_map[i] == 'FAIL' :
                failures += 1

        return not_done, failures

    @do_locking
    def alignments_complete(self) :
        genefamily_contig_map = self.build_genefamily2contigs()

        for i in genefamily_contig_map :
            if i not in self.genefamily_filename_map :
                return False

        return True

    # functions used by scaffolder
    #
    @do_locking
    def contig_used(self, contig_id, label) :
        return contig_id in self.contig_query_map[label]

    @do_locking
    def contig_assigned(self, contig_id, label) :
        qid = self.contig_query_map[label][contig_id]

        return self.query_gene_map[qid] != None

    @do_locking
    def query_to_gene(self, query_id) :
        geneid,strand = self.query_gene_map[query_id]
        return geneid

#    @do_locking
#    def contig_aligned(self, contig_id) :
#        qid = self.contig_query_map[contig_id]
#        gid = self.query_gene_map[qid]
#        gfid = self.db.get_genefamily_from_gene(gid)
#        
#        return self.genefamily_filename_map[gfid] != 'FAIL'

    @do_locking
    def get_contig_from_query(self, query_id) :
        # lazy reverse lookup
        if not hasattr(self, 'query_contig_map') :
            self.query_contig_map = {}

            for label in self.contig_query_map :
                cqm = self.contig_query_map[label]
                for contig_id in cqm :
                    self.query_contig_map[cqm[contig_id]] = (contig_id, label)

        if isinstance(query_id, list) :
            return [ self.query_contig_map[i] for i in query_id ]

        return self.query_contig_map[query_id]

