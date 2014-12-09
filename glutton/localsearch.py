import subprocess
import os
import threading
import sys

from glutton.utils import get_log, tmpfile, openmp_num_threads, rm_f
from glutton.blast import Blast
from glutton.job import BlastJob
from glutton.queue import WorkQueue


class All_vs_all_search(object) :
    def __init__(self, batch_size=100) :
        self.batch_size = batch_size
        self.log = get_log()
        self.cleanup_files = []
        self.gene_assignments = {}
        self.lock = threading.Lock()
        self.q = None

        self.total_jobs = 0
        self.complete_jobs = 0

    def _batch(self, x) :
        tmp = []

        for i in x :
            tmp.append(i)
            
            if len(tmp) == self.batch_size :
                yield tmp
                tmp = []
        
        if not tmp :
            raise StopIteration

        yield tmp

    def process(self, db, queries, nucleotide, min_identity, min_alignment_length) :
        self.min_identity = 100 * min_identity
        self.min_length = min_alignment_length

        # we need to deal with the index files here because 
        # all of the blastx jobs need them
        self.cleanup_files += [db + i for i in [".phr",".pin",".psq"]]

        # creates db + {phr,pin,psq} in same dir as db
        self.log.info("creating blast db...")
        Blast.makedb(db)

        # queue up the jobs
        self.log.info("starting local alignments...")
        self.q = WorkQueue()

        self.total_jobs = len(queries)
        self.complete_jobs = -self.batch_size
        self._progress()

        for query in self._batch(queries) :
            self.q.enqueue(BlastJob(self.job_callback, db, query, 'blastx'))

        self.log.info("waiting for job queue to drain...")
        self.q.join()

        rm_f(self.cleanup_files)

        return self.gene_assignments

    def stop(self) :
        if self.q :
            self.q.stop()

        rm_f(self.cleanup_files)

    def get_intermediate_results(self) :
        return self.gene_assignments

    def _progress(self) :
        self.complete_jobs += self.batch_size
        
        if self.complete_jobs > self.total_jobs :
            self.complete_jobs = self.total_jobs

        sys.stderr.write("\rProgress: %d / %d blastx alignments " % (self.complete_jobs, self.total_jobs))

        if self.complete_jobs == self.total_jobs :
            print >> sys.stderr, "\ndone!"

    def job_callback(self, job) :
        self.log.debug("%d blast results returned" % len(job.results))

        self.lock.acquire()

        self._progress()

        if job.success() :
            for contig,geneid,identity,length in job.results :
                if length < self.min_length or identity < self.min_identity :
                    continue

                self.gene_assignments[contig] = geneid

        for q in job.input :
            if q.id not in self.gene_assignments :
                self.gene_assignments[q.id] = 'FAIL'

        self.lock.release()


if __name__ == '__main__' :

    from glutton.utils import get_log, glutton_log_defaults, set_threads
    from glutton.genefamily import read_alignment_as_genefamily

    glutton_log_defaults(get_log())

    ava = All_vs_all_search()
    tmp = ava.process('tc_test.fasta', read_alignment_as_genefamily('queries_test.fasta', 'test'), 0.8, 100)
    
    for m in tmp :
        print m, tmp[m]

