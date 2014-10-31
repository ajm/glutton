import collections
import shutil
import sys

from glutton.db import GluttonDB, GluttonDBFileError
from glutton.localsearch import All_vs_all_search
from glutton.utils import tmpfile, num_threads, get_log, rm_f, check_dir, md5
from glutton.queue import WorkQueue
from glutton.job import PaganJob
from glutton.genefamily import Gene, biopy_to_gene, seqlen
from glutton.info import GluttonInformation

from os.path import isfile, basename

from Bio import SeqIO


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

        self.log = get_log()

        self.info = GluttonInformation(db, contigs_fname, directory)

    def _read_contigs(self, fname) :
        contigs = {}
        accepted = 0
        rejected = 0

        for r in SeqIO.parse(fname, 'fasta') :
            if len(r) < self.min_length :
                rejected += 1
                continue

            qid = self.info.get_contig2query(r.id)
            
            contigs[qid] = biopy_to_gene(r, qid)
            accepted += 1

        self.log.info("read %d contigs (rejected %d due to length < %d)" % 
                (accepted, rejected, self.min_length))

        return contigs

    def stop(self) :
        self.search.stop()
        self.info.update_query2gene(self.search.get_intermediate_results())
        
        rm_f(self.cleanup_files)

        if self.q :
            self.q.stop()

        self.info.flush()

    def align(self) :
        self.log.info("starting alignment procedure")

        # convert the names of the contigs to something no program can complain about
        # + filter out the ones that could never have a long enough alignment
        contigs = self._read_contigs(self.contigs_fname)

        pending_contigs = [ contigs[i] for i in self.info.pending_queries() ]

        self.log.info("%d contigs have not been assigned to genes..." % len(pending_contigs))

        # depending on when the program was terminated this step may be complete or partially
        # complete 
        if pending_contigs :
            db_fname = self.db.extract_all()
            self.cleanup_files.append(db_fname)

            # do an all vs all search of contigs vs database of transcripts
            # return a dict of tmp ids with gene ids
            self.info.update_query2gene(
                self.search.process(
                    db_fname, 
                    pending_contigs,
                    self.min_identity, 
                    self.min_length / 3)
                )

            rm_f(db_fname)


        # use the database to convert the mapping from tmp id -> gene
        # to gene family -> list of tmp ids
        genefamily_contig_map = self.info.build_genefamily2contigs()
        
        self.log.info("%d contigs assigned to %d gene families" % 
                (sum([ len(i) for i in genefamily_contig_map.values() ]), len(genefamily_contig_map)))
        self.log.info("(%d have already been run)" % self.info.len_genefamily2filename())

        if self.info.len_genefamily2filename() == len(genefamily_contig_map) :
            self.log.info("alignment already done, exiting early...")
            return
        else :
            self.log.info("starting alignments...")


        # queue all the alignments up using a work queue and pagan
        self.q = WorkQueue()

        for famid in self.sort_keys_by_complexity(genefamily_contig_map) :
            # ignore the jobs that have already been run
            if self.info.in_genefamily2filename(famid) :
                continue
            
            # get the alignment and tree from the database
            try :
                alignment = self.db.get_alignment(famid)
                tree = alignment.get_tree()

            except GluttonDBFileError, gde :
                self.log.warn(str(gde))
                continue

            job_contigs = [ contigs[i] for i in genefamily_contig_map[famid] ]

            # queue the job
            self.q.enqueue(
                    PaganJob(
                        self.job_callback, 
                        job_contigs, 
                        famid,
                        alignment, 
                        tree)
                    )
            
        self.log.info("waiting for job queue to drain...")
        self.q.join()

        self.info.flush()

    def sort_keys_by_complexity(self, d) :
        return [ k for k,v in sorted(d.items(), reverse=True, key=lambda x : seqlen(x[1])) ]

    def job_callback(self, job) :
        self.log.debug("callback from %s: %s + %s" % (str(job), job.genefamily, ','.join([ i.id for i in job.input ])))
        self.log.debug("alignment file = %s" % job.alignment)

        if job.success() :
            dst = tmpfile(directory=self.directory)
            shutil.copyfile(job.alignment, dst)
        
            self.info.put_genefamily2filename(job.genefamily, basename(dst))
        else :
            self.info.put_genefamily2filename(job.genefamily)


if __name__ == '__main__' :
    
    import sys
    import signal
    import os

    from glutton.utils import get_log, setup_logging, set_verbosity, set_threads


    set_verbosity(3)
    setup_logging()

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

