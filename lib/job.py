import os

from lib.base import Base
from lib.manifest import Manifest, ManifestError
from lib.prank import Prank
from lib.pagan import Pagan

from abc import abstractmethod
from os.path import basename, isfile, join

class JobError(Exception) :
    pass

class Job(Base) :
    QUEUED,RUNNING,SUCCESS,FAIL,TERMINATED,INTERNAL_ERROR = range(6)

    states = {
        QUEUED     : 'QUEUED',
        RUNNING    : 'RUNNING',
        SUCCESS    : 'SUCCESS',
        FAIL       : 'FAIL',
        TERMINATED : 'TERMINATED',
        INTERNAL_ERROR : 'INTERNAL_ERROR'
    }

    def __init__(self, opt) :
        super(Job, self).__init__(opt)

        self.state = Job.QUEUED

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
        
        try :
            ret = self._run()
    
        except Exception, e :
            self.error(str(e))
            self.end(Job.INTERNAL_ERROR)
            return

        if ret == 0 :
            self.end(Job.SUCCESS)
        elif ret > 0 : 
            self.end(Job.FAIL)
        else :
            self.end(Job.TERMINATED)

    def state_str(self) :
        return Job.states[self.state]

    @abstractmethod
    def _run(self) :
        pass

    def __str__(self) :
        return "%s, %s" % (type(self).__name__, self.fname)

class PrankJob(Job) :
    def __init__(self, opt, manifest, fname) :
        super(PrankJob, self).__init__(opt)

        self.manifest = manifest
        self.fname = fname

        assert os.path.isabs(self.fname) and os.path.isfile(self.fname)

    def _run(self) :
        p = Prank(self.opt, self.fname)

        ret = p.run()

        if ret == 0 :
            for f in p.output_filenames() :
                self.manifest.append_to_manifest(f, self._contents(f), create=True)

        return ret

class PaganJob(Job) :
    def __init__(self, opt, transcriptome, fname, outdir) :
        super(PaganJob, self).__init__(opt)

        self.fname = fname
        self.transcriptome = transcriptome
        self.outdir = outdir

        assert os.path.isabs(self.fname)  and os.path.isfile(self.fname)
        assert os.path.isabs(self.outdir) and os.path.isdir(self.outdir)

    def _run(self) :
        ret = 1

        try :
            job_list = self.transcriptome.query(self.fname)

        except ExonerateError, ee :
            self.info("%s : %s" % (basename(self.fname), str(ee)))
            return 1


        for index,job in enumerate(job_list) :
            root_fname, num_genes = job
            a_fname = root_fname
            t_fname = None

            self.info("aligning contig in '%s' against gene family in '%s'..." % (basename(self.fname), basename(root_fname)))

            if num_genes > 1 :
                a_fname = root_fname + '.pep.2.fas'
                t_fname = root_fname + '.2.dnd'

            p = Pagan(self.opt,
                  alignment_file=a_fname, 
                  tree_file=t_fname, 
                  query_file=self.fname,
                  out_dir=self.outdir)

            p.run()

            tmp = join(self.outdir, basename(self.fname))
            fas = tmp + '.fas'
            dna = tmp + '.dna.fas'

            if isfile(fas) and isfile(dna) :
                ret = 0
                os.rename(fas, "%s.%d.fas" % (tmp, index))
                os.rename(dna, "%s.%d.dna.fas" % (tmp, index))

                self.info("align '%s' --> '%s' success" % (basename(self.fname), basename(root_fname)))
            else :
                self.info("align '%s' --> '%s' fail" % (basename(self.fname), basename(root_fname)))

        if ret != 0 :
            self._safe_write(join(self.outdir, 'failures.txt'), basename(self.fname))

        return ret

