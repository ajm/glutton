import sys
import threading
import Queue
import time
import itertools

from glutton.job import Job, JobError
from glutton.utils import get_log, num_threads


class WorkQueueError(Exception) :
    pass

class WorkQueue(object):
    def __init__(self, qtimeout=1, maxsize=0):

        if maxsize == 0 :
            maxsize = num_threads() * 2

        self.log = get_log()

        self.q = Queue.Queue(maxsize)
        self.workers = self._init_workers(num_threads())
        self.q_timeout = qtimeout
        self.running = False
        self.no_more_jobs = False
        
        self.jobs_completed = 0
        self.jobs_counter = itertools.count(start=1)

        self.start()

    def __del__(self) :
        self.stop()

#    def _introspect_cores(self) :
#        return cpu_count()

    def _init_workers(self, numworkers):
#        if numworkers == 0 :
#            numworkers = self._introspect_cores()

        self.log.info("queue using %d thread%s" % (numworkers, "" if numworkers == 1 else "s"))

        tmp = []

        for _ in range(numworkers):
            t = threading.Thread(target=self._consume_queue)
            t.setDaemon(True)
            tmp.append(t)

        return tmp

    def start(self):
        self.running = True
        
        for t in self.workers:
            t.start()
    
        self.log.debug("queue started")

    # block until the currently running jobs complete
    def stop(self):
        self.log.debug("queue stopping...")
        self.running = False
        
        for t in self.workers:
            t.join()
    
        self.log.debug("queue stopped")

    def join(self) :
        # calling join causes the calling thread to try and 
        # acquire the lock in self.q, so poll instead so the
        # user can still do a ctrl-C

        while True :
            if self.q.empty() :
                break

            time.sleep(5)

        self.q.join()

    def size(self) :
        return self.q.qsize()

    # block until the queue is drained
    def done(self) :
        self.no_more_jobs = True

        for t in self.workers :
            t.join()

        self.log.debug("queue drained")

    def enqueue(self, j):
        self.log.debug("enqueuing %s" % str(j))
        
        assert isinstance(j, Job)

        while True :
            try :
                self.q.put(j, timeout=3600)
                break

            except Queue.Full :
                pass

    def _consume_queue(self):
        while self.running :
            try :
                work = self.q.get(timeout=self.q_timeout)
            
            except Queue.Empty, qe:
                if self.no_more_jobs :
                    self.log.debug("no more jobs, exiting...")
                    break

                continue

            self.log.debug("starting %s" % str(work))
            work.run()

            self.log.debug("completed %s %s" % (str(work), work.state_str()))
            self.q.task_done()
            self.jobs_completed = self.jobs_counter.next()

            if work.terminated() :
                self.log.warn("job terminated, thread exiting...")
                break

