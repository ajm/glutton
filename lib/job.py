import os

from lib.base import Base
from lib.manifest import Manifest, ManifestError
from lib.prank import Prank

from abc import abstractmethod


class JobError(Exception) :
    pass

class Job(Base) :
    QUEUED,RUNNING,SUCCESS,FAIL = range(4)

    def __init__(self, opt) :
        super(Job, self).__init__(opt)

        self.status = Job.QUEUED

    def success(self) :
        if self.status not in (Job.SUCCESS, Job.FAIL) :
            raise JobError('job has not been run')

        return self.status == Job.SUCCESS

    def fail(self) :
        return not self.success()

    def start(self) :
        self.status = Job.RUNNING

    def end(self, s) :
        assert s in (Job.SUCCESS, Job.FAIL), "status should be success or fail"
        self.status = s

    def run(self) :
        self.start()
        self.end(Job.SUCCESS if self._run() else Job.FAIL)

    @abstractmethod
    def _run(self) :
        pass

    def followup_jobs(self) :
        if self.fail() :
            raise JobError("failed jobs have no followup jobs")

        return self._followup_jobs()

    @abstractmethod
    def _followup_jobs(self) :
        pass

    def __str__(self) :
        return "%s, %s" % (type(self).__name__, self.fname)


class SaveJob(Job) :
    def __init__(self, opt, manifest, fname, fcontents, align=False) :
        super(SaveJob, self).__init__(opt)

        self.manifest = manifest
        self.fname = fname
        self.fcontents = fcontents
        self.align = align

        assert os.path.abspath(self.fname)

    def _run(self) :
        self.manifest.append_to_manifest(self.fname, self.fcontents, create=True)

    def _followup_jobs(self) :
        return [ PrankJob(self.opt, self.manifest, self.fname) ] if self.align else []

class PrankJob(Job) :
    def __init__(self, opt, manifest, fname) :
        super(PrankJob, self).__init__(opt)

        self.manifest = manifest
        self.fname = fname
        self.followup = []

        assert os.path.abspath(self.fname)

    def _run(self) :
        p = Prank(self.fname)

        if p.align() :
            self.followup = p.outfile_filenames()

    def _followup_jobs(self) :
        return [ SaveJob(self.opt, self.manifest, os.path.basename(f), self._contents(f)) for f in self.followup ]

class PaganJob(Job) :
    def __init__(self, opt) :
        super(PaganJob, self).__init__(opt)

    def _run(self) :
        raise NotImplementedError

    def _followup_jobs(self) :
        raise NotImplementedError

class ExonerateJob(Job) :
    def __init__(self, opt) :
        super(ExonerateJob, self).__init__(opt)

    def _run(self) :
        raise NotImplementedError

    def _followup_jobs(self) :
        raise NotImplementedError

