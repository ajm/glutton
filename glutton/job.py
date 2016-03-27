import os

from glutton.utils import get_log, tmpfasta, tmpfasta_orfs, tmpfile, rm_f, threadsafe_io, fasta_stats
from glutton.prank import Prank
from glutton.pagan import Pagan
from glutton.blast import Blastx, Tblastx

from abc import abstractmethod
from os.path import basename, isfile, join
from sys import exit

import time


class JobError(Exception) :
    pass

class Job(object) :
    QUEUED,RUNNING,SUCCESS,FAIL,TERMINATED,INTERNAL_ERROR,NOTHING_TO_DO = range(7)

    states = {
        QUEUED          : 'QUEUED',
        RUNNING         : 'RUNNING',
        SUCCESS         : 'SUCCESS',
        FAIL            : 'FAIL',
        TERMINATED      : 'TERMINATED',
        INTERNAL_ERROR  : 'INTERNAL_ERROR',
        NOTHING_TO_DO   : 'NOTHING_TO_DO'
    }

    def __init__(self, callback) :
        self.state = Job.QUEUED
        self.log = get_log()
        self.callback = callback

    def success(self) :
        if self.state not in (Job.SUCCESS, Job.FAIL, Job.TERMINATED, Job.INTERNAL_ERROR) :
            raise JobError('job has not been run')

        return self.state == Job.SUCCESS

    def fail(self) :
        return not self.success()

    def terminated(self) :
        if self.state not in (Job.SUCCESS, Job.FAIL, Job.TERMINATED, Job.INTERNAL_ERROR) :
            raise JobError('job has not been run')
        
        return self.state == Job.TERMINATED

    def start(self) :
        self.state = Job.RUNNING

    def end(self, s) :
        assert s in (Job.SUCCESS, Job.FAIL, Job.TERMINATED, Job.INTERNAL_ERROR), "status should be success, fail or terminated"
        self.state = s

    def run(self) :
        self.start()
        
        ret = self._run()
        #try :
        #    ret = self._run()

        #except Exception, e :
        #    self.log.error(str(e))
        #    self.end(Job.INTERNAL_ERROR)
        #    self.cleanup()
        #    return

        if ret == 0 :
            self.end(Job.SUCCESS)
        elif ret == -2 : # SIGINT = 130
            self.end(Job.TERMINATED)
        else :
            self.end(Job.FAIL)

        if not self.terminated() :
            self.callback(self)

        self.cleanup()

    def state_str(self) :
        return Job.states[self.state]

    # delete only the files the program created
    # the responsibility to delete the input files is for the caller
    # - this is no longer true because the callers are passing data not filenames
    def cleanup(self) :
        for f in self._get_filenames() :
            if f and isfile(f) :
                self.log.debug("deleting %s" % f)
                rm_f(f)

    @abstractmethod
    def _run(self) :
        pass

    @abstractmethod
    def _get_filenames(self) :
        pass

    def __str__(self) :
        return type(self).__name__ #"%s, %s" % (type(self).__name__, self.fname)

class PrankJob(Job) :
    def __init__(self, callback, sequences) :
        super(PrankJob, self).__init__(callback)

        self.sequences = sequences
        self.prank = Prank()

    @property
    def input(self) :
        return self.sequences

    @property
    def tree(self) :
        return self.prank.tree

    @property
    def alignment(self) :
        return self.prank.alignment

    def _get_filenames(self) :
        return [self.infile] + self.prank.output_filenames(self.infile)

    def _run(self) :
        self.infile = tmpfasta(self.sequences)


        start_time = time.time()

        result = self.prank.run(self.infile, self.infile)

        elapsed_time = time.time() - start_time
        q_count, q_sum, q_min, q_max, q_mean, q_sd = fasta_stats(self.infile)

        threadsafe_io('prank_stats.txt', "%d %d %d %d %d %.3f %.3f %d" % \
                                            (result, \
                                             q_count, q_sum, q_min, q_max, q_mean, q_sd, \
                                             elapsed_time))

        return result

class BlastJob(Job) :
    def __init__(self, callback, database, queries, blast_version='blastx') :
        super(BlastJob, self).__init__(callback)

        self.database = database
        self.queries = queries

        assert blast_version in ('blastx', 'tblastx')

        self.blastx = Tblastx() if blast_version == 'tblastx' else Blastx()

    @property
    def input(self) :
        return self.queries

    @property
    def results(self) :
        return self.blastx.results

    def _get_filenames(self) :
        return [self.query_fname, self.out_fname]

    def _run(self) :
        self.query_fname = tmpfasta(self.queries)
        self.out_fname = tmpfile()

        result = self.blastx.run(self.query_fname, self.database, self.out_fname)

        q = dict([ (q.id, len(q)) for q in self.input ])
        for br in self.results :
            threadsafe_io('blastx_stats.txt', "%s %s %.3f %d %d %d %d %d %.3e %d %.3f" % \
                                                                          (br.qseqid, 
                                                                           br.sseqid, 
                                                                           br.pident, 
                                                                           br.length, 
                                                                           br.qstart,
                                                                           br.qend,
                                                                           br.sstart,
                                                                           br.send,
                                                                           br.evalue,
                                                                           q[br.qseqid],
                                                                           ((br.pident / 100.0) * (max(br.qstart, br.qend) - min(br.qstart, br.qend))) / float(q[br.qseqid])))

        return result

class PaganJob(Job) :
    def __init__(self, callback, queries, genefamily_id, alignment, tree, identity, overlap) :
        super(PaganJob, self).__init__(callback)

        self._queries = queries
        self._genefamily = genefamily_id
        self._alignment = alignment
        self._tree = tree
        self.identity = identity
        self.overlap = overlap

        self.pagan = Pagan()

        self.query_fname = None
        self.out_fname = None
        self.alignment_fname = None
        self.tree_fname = None

    @property
    def input(self) :
        return self._queries

    @property
    def genefamily(self) :
        return self._genefamily

    @property
    def nucleotide_alignment(self) :
        return self.pagan.nucleotide_alignment

    @property
    def protein_alignment(self) :
        return self.pagan.protein_alignment

    def _get_filenames(self) :
        #return self.pagan.output_filenames(self.out_fname)
        return [self.query_fname, self.alignment_fname, self.tree_fname] + self.pagan.output_filenames(self.out_fname)

    def _run(self) :

        self.query_fname     = tmpfasta_orfs(self._queries, strand=True)
        #self.query_fname     = tmpfasta(self._queries)
        self.out_fname       = tmpfile()
        self.alignment_fname = tmpfasta(self._alignment) # tmpfasta_kill_n(self._alignment)
        
        self.tree_fname = tmpfile(self._tree) if self._tree else None
        
        
        start_time = time.time()
        
        result = self.pagan.run(self.query_fname, 
                                self.out_fname, 
                                self.alignment_fname, 
                                self.tree_fname,
                                self.identity,
                                self.overlap)
        
        elapsed_time = time.time() - start_time
        q_count, q_sum, q_min, q_max, q_mean, q_sd = fasta_stats(self.query_fname)
        a_count, a_sum, a_min, a_max, a_mean, a_sd = fasta_stats(self.alignment_fname)

        threadsafe_io('pagan_stats.txt', "%s %d %d %d %d %d %.3f %.3f %d %d %d %d %.3f %.3f %d" % \
                                            (self._genefamily, result, \
                                             q_count, q_sum, q_min, q_max, q_mean, q_sd, \
                                             a_count, a_sum, a_min, a_max, a_mean, a_sd, \
                                             elapsed_time))
        
        return result

