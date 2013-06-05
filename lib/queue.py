import sys
import threading
import Queue

from multiprocessing import cpu_count

from lib.job import Job, JobError


class WorkQueueError(Exception) :
    pass

class WorkQueue(object):
    def __init__(self, opt, numworkers, qtimeout=1, maxsize=0):
        super(WorkQueue, self).__init__(opt)

        self.q = Queue.Queue(maxsize) # I think I can deadlock this if maxsize != 0 
                                      # (because when jobs end, there may be a number of followup jobs)
        self.workers = self._init_workers(numworkers)
        self.q_timeout = qtimeout
        self.running = False
        self.no_more_jobs = False

    def _introspect_cores(self) :
        return cpu_count()

    def _init_workers(self, numworkers):
        if numworkers == 0 :
            numworkers = self._introspect_cores()

        self.info("WorkQueue: using %d threads" % numworkers)

        tmp = []

        for _ in range(numworkers):
            t = threading.Thread(target=self._consume_queue)
            t.daemon = True
            tmp.append(t)

        return tmp

    def start(self):
        self.running = True
        
        for t in self.workers:
            t.start()
    
        self.info("WorkQueue: started")

    def stop(self):
        self.running = False
        self.q.join()
        
        for t in self.workers:
            t.join()
    
        self.info("WorkQueue: stopped")

    # this tells the queue that there are no more jobs to be added
    # (at least externally, it might be the case that they are added
    #  as follow up jobs, but this will only make the last few jobs
    #  slightly suboptimal)
    def done(self) :
        self.no_more_jobs = True
        self.q.join()

        for t in self.workers :
            t.join()

        self.info("WorkQueue: done")

    def enqueue(self, j):
        self.info("WorkQueue: enqueuing %s" % key)
        
        assert isinstance(j, Job)

        self.q.put(j)

    def _consume_queue(self):
        while self.running :
            try :
                work = self.q.get(timeout=self.q_timeout)
            
            except Queue.Empty, qe:
                if self.no_more_jobs :
                    break

                continue

            self.info("WorkQueue: starting %s" % str(work))

            try :
                work.run()
            
                for job in work.followup_jobs() :
                    self.enqueue(job)

            except JobError, je :
                pass

            self.q.task_done()

            self.info("WorkQueue: completed %s (status: %s)" % (str(work), work.status()))

