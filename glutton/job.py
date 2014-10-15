import os

from glutton.utils import get_log, tmpfasta, tmpfile, rm_f
from glutton.prank import Prank
from glutton.pagan import Pagan
from glutton.blastx import Blastx

from abc import abstractmethod
from os.path import basename, isfile, join

class JobError(Exception) :
    pass

class Job(object) :
    QUEUED,RUNNING,SUCCESS,FAIL,TERMINATED,INTERNAL_ERROR = range(6)

    states = {
        QUEUED          : 'QUEUED',
        RUNNING         : 'RUNNING',
        SUCCESS         : 'SUCCESS',
        FAIL            : 'FAIL',
        TERMINATED      : 'TERMINATED',
        INTERNAL_ERROR  : 'INTERNAL_ERROR'
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
        return self.state == Job.TERMINATED

    def start(self) :
        self.state = Job.RUNNING

    def end(self, s) :
        assert s in (Job.SUCCESS, Job.FAIL, Job.TERMINATED, Job.INTERNAL_ERROR), "status should be success, fail or terminated"
        self.state = s

    def run(self) :
        self.start()
        
        ret = self._run()

#       XXX
#        try :
#            ret = self._run()
#
#        except Exception, e :
#            self.log.error(str(e))
#            self.end(Job.INTERNAL_ERROR)
#            self.cleanup()
#            return

        if ret == 0 :
            self.end(Job.SUCCESS)
        elif ret > 0 : 
            self.end(Job.FAIL)
        else :
            self.end(Job.TERMINATED)

        self.callback(self)
        self.cleanup()

    def state_str(self) :
        return Job.states[self.state]

    # delete only the files the program created
    # the responsibility to delete the input files is for the caller
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

        return self.prank.run(self.infile, self.infile)

class BlastxJob(Job) :
    def __init__(self, callback, database, queries) :
        super(BlastxJob, self).__init__(callback)

        self.database = database
        self.queries = queries

        self.blastx = Blastx()

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

        return self.blastx.run(self.query_fname, self.database, self.out_fname)

class PaganJob(Job) :
    def __init__(self, callback, queries, genefamily_id, alignment, tree) :
        super(PaganJob, self).__init__(callback)

        self._queries = queries
        self._genefamily = genefamily_id
        self._alignment = alignment
        self._tree = tree

        self.pagan = Pagan()

    @property
    def input(self) :
        return self._queries

    @property
    def genefamily(self) :
        return self._genefamily

    @property
    def alignment(self) :
        return self.pagan.alignment

    def _get_filenames(self) :
        return [self.query_fname, self.alignment_fname, self.tree_fname] + self.pagan.output_filenames(self.out_fname)

    def _run(self) :
        self.query_fname        = tmpfasta(self._queries)
        self.out_fname          = tmpfile()
        self.alignment_fname    = tmpfasta(self._alignment)
        
        self.tree_fname = tmpfile(self._tree) if self._tree else None
        

        return self.pagan.run(self.query_fname, 
                              self.out_fname, 
                              self.alignment_fname, 
                              self.tree_fname)

