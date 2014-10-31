import os
import collections
import shutil
import json
import sys
import threading

from glutton.db import GluttonDB, GluttonDBFileError
from glutton.localsearch import All_vs_all_search
from glutton.utils import tmpfile, num_threads, get_log, rm_f, check_dir, md5
from glutton.queue import WorkQueue, order_jobs
from glutton.job import PaganJob
from glutton.genefamily import Gene, biopy_to_gene, seqlen

from os.path import isfile

from Bio import SeqIO


PARAM_FILE  = 'parameters.json'
CONTIG_FILE = 'contigs.json'
BLAST_FILE  = 'blastx.json'
PAGAN_FILE  = 'pagan.json'


class Aligner(object) :
    def __init__(self, db, contigs_fname, directory, min_identity=0.9, min_length=1, batch_size=100) :
        self.db = db
        self.contigs_fname = contigs_fname
        self.directory = directory
        self.min_identity = min_identity
        self.min_length = min_length

        check_dir(self.directory, create=True)

        self.search = All_vs_all_search(batch_size)
        self.cleanup_files = []
        self.q = None
        self.lock = threading.Lock()

        self.log = get_log()

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
        self.genefamily_filename_map  = {}  # gene family id -> filename

        self.read_progress_files()

        db_params = self.get_params()

        # if we are starting fresh
        if not self.params :
            self.params = db_params
        # else we are restarting
        else :
            if self._not_same_db(db_params) :
                self.log.fatal("alignment has been resumed with a different reference!")
                self.log.fatal("\toriginal = %s" % self._param_str(self.params))
                self.log.fatal("\tcurrent  = %s" % self._param_str(db_params))
                sys.exit(1)

    # files used for recording the progress are just dirty globals
    # at the moment these properties can be in place of a decent solution
    # for now...
    @property
    def parameter_filename(self) :
        global PARAM_FILE
        return os.path.join(self.directory, PARAM_FILE)
    
    @property
    def contig_filename(self) :
        global CONTIG_FILE
        return os.path.join(self.directory, CONTIG_FILE)

    @property
    def blastx_filename(self) :
        global BLAST_FILE
        return os.path.join(self.directory, BLAST_FILE)

    @property
    def pagan_filename(self) :
        global PAGAN_FILE
        return os.path.join(self.directory, PAGAN_FILE)

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

    def write_progress_files(self) :
        self.lock.acquire()

        self._dump(self.parameter_filename, self.params)
        self._dump(self.contig_filename,    self.contig_query_map)
        self._dump(self.blastx_filename,    self.query_gene_map)
        self._dump(self.pagan_filename,     self.genefamily_filename_map)

        self.lock.release()

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

    def _param_str(self, p) :
        return "%s/%d %s md5=%s" % (p['species'], p['release'], p['filename'], p['checksum'])

    def _not_same_db(self, par) :
        for key in self.params :
            if self.params[key] != par[key] :
                return True
        return False

    def _read_contigs(self, fname) :
        queryid = 0
        contigs = {}
        accepted = 0
        rejected = 0

        for r in SeqIO.parse(fname, 'fasta') :
            if len(r) < self.min_length :
                rejected += 1
                continue

            # assume that if there is anything in contig_query_map then
            # everything is, else exit
            if self.contig_query_map :
                try :
                    qid = self.contig_query_map[r.id]

                except KeyError, ke :
                    self.log.error("It appears that all contigs have query ids, but %s does not..." % r.id)
                    self.log.error("This situation is currently unrecoverable!")
                    sys.exit(1)
            else :
                qid = "query%d" % queryid
                queryid += 1

            contigs[qid] = biopy_to_gene(r, qid)
            accepted += 1

        # if we have not stored the query ids before, then we need to 
        # do so with a lock (writing these dicts to disk uses the lock as well)
        if not self.contig_query_map :
            self.lock.acquire()

            for k in contigs :
                self.contig_query_map[contigs[k].name] = k

            self.lock.release()

        self.log.info("read %d contigs (rejected %d due to length < %d)" % (accepted, rejected, self.min_length))

        return contigs

    def stop(self) :
        self.search.stop()
        self.query_gene_map.update(self.search.get_intermediate_results())
        
        rm_f(self.cleanup_files)

        if self.q :
            self.q.stop()

        self.write_progress_files()

    def align(self) :
        # convert the names of the contigs to something no program can complain about
        # + filter out the ones that could never have a long enough alignment
        contigs = self._read_contigs(self.contigs_fname)

        pending_contigs = [ contigs[i] for i in contigs if i not in self.query_gene_map ]

        self.log.info("%d contigs have not been assigned to genes..." % len(pending_contigs))

        # depending on when the program was terminated this step may be complete or partially
        # complete 
        if pending_contigs :
            db_fname = self.db.extract_all()
            self.cleanup_files.append(db_fname)

            # do an all vs all search of contigs vs database of transcripts
            # return a dict of tmp ids with gene ids
            self.query_gene_map.update(
                self.search.process(
                    db_fname, 
                    #contigs.values(), 
                    pending_contigs,
                    self.min_identity, 
                    self.min_length / 3)
                )

            rm_f(db_fname)


        # use the database to convert the mapping from tmp id -> gene
        # to gene family -> list of tmp ids
        genefamily_contig_map = collections.defaultdict(list)
        assigned_contigs = 0
        for i in self.query_gene_map :
            if self.query_gene_map[i] == 'FAIL' :
                continue

            assigned_contigs += 1

            genefamily_contig_map[self.db.get_familyid_from_geneid(self.query_gene_map[i])].append(contigs[i])


        self.log.info("%d contigs assigned to %d gene families" % (assigned_contigs, len(genefamily_contig_map)))
        self.log.info("(%d have already been run)" % len(self.genefamily_filename_map))
        self.log.info("starting alignments...")

        # queue all the alignments up using a work queue and pagan
        self.q = WorkQueue()

        for famid in self.sort_keys_by_complexity(genefamily_contig_map) :
            # ignore the jobs that have already been run
            if famid in self.genefamily_filename_map :
                continue
            
            # get the alignment and tree from the database
            try :
                alignment = self.db.get_alignment(famid)
                tree = alignment.get_tree()

            except GluttonDBFileError, gde :
                self.log.warn(str(gde))
                continue

            # queue the job
            self.q.enqueue(
                    PaganJob(
                        self.job_callback, 
                        genefamily_contig_map[famid], 
                        famid,
                        alignment, 
                        tree)
                    )
            
        self.log.info("waiting for job queue to drain...")
        self.q.join()

        self.write_progress_files()

    def sort_keys_by_complexity(self, d) :
        return [ k for k,v in sorted(d.items(), reverse=True, key=lambda x : seqlen(x[1])) ]

    def job_callback(self, job) :
        self.log.debug("callback from %s: %s + %s" % (str(job), job.genefamily, ','.join([ i.id for i in job.input ])))
        self.log.debug("alignment file = %s" % job.alignment)

        self.lock.acquire()

        if job.success() :
            dst = tmpfile(directory=self.directory)
            shutil.copyfile(job.alignment, dst)
        
            self.genefamily_filename_map[job.genefamily] = os.path.basename(dst)
        else :
            self.genefamily_filename_map[job.genefamily] = 'FAIL'


        self.lock.release()


if __name__ == '__main__' :
    
    import sys
    import signal
    import os

    from glutton.utils import get_log, glutton_log_defaults, set_threads


    glutton_log_defaults(get_log())

    #set_threads(1)

    db = GluttonDB('tc23.glt')
    contigs = 'queries_test.fasta'
    min_id = 0.7
    min_len = 100
    #batch_size = 10 # 8 threads, 5m 11s
    batch_size = 100 # 8 threads, 4m 49s
    #batch_size = 500 # 8 threads, 4m 45s
    directory = './alignment_test'


    align = Aligner(db, contigs, directory, min_id, min_len, batch_size)
    
    def _cleanup(signal, frame) :
        print >> sys.stderr, "Killed by user, cleaning up..."
        align.stop()
        print >> sys.stderr, "clean up done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    align.align()

    get_log().info("DONE!")

