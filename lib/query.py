import os
from glob import glob

from cogent.parse.fasta import MinimalFastaParser

from lib.base import Base

class QueryManager(Base) :
    query_prefix = 'query_'

    def __init__(self, opt, fname, outdir, min_length) :
        super(QueryManager, self).__init__(opt)

        assert os.path.isabs(fname)

        self.fname = fname
        self.min_length = min_length
        self.dir = self.tmpdir
        self.outdir = outdir

        self.queries = self._complete_jobs()

    def __iter__(self) :
        return self._queries()

    def _write_query(self, fname, label, seq) :
        self._create_file(fname, ">%s\n%s" % (label, seq), self.dir)

    def _complete_jobs(self) :
        tmp = set()

        # failed jobs
        fail_count = 0
        try :
            f = open(os.path.join(self.outdir, 'failures.txt'))
        
            for line in f :
                line = line.strip()
                if line :
                    tmp.add(line)
                    fail_count += 1

            f.close()
        
        except IOError, ioe :
            pass


        # successful jobs
        success_count = 0
        for f in glob(os.path.join(self.outdir, self.query_prefix + '*.fas')) :

            if f.endswith('.dna.fas') :
                continue

            f = os.path.basename(f)
            f,ext = os.path.splitext(f)
            tmp.add(f)
            success_count += 1

        if (fail_count + success_count) != 0 :
            self.info("resuming - added %d failed job%s" % (fail_count, "" if fail_count == 1 else "s"))
            self.info("resuming - added %d successful job%s" % (success_count, "" if success_count == 1 else "s"))

        return tmp

    def _queries(self) :
        for label,seq in MinimalFastaParser(open(self.fname)) :
            if len(seq) < self.min_length :
                continue

            # exonerate does not like filenames with ':' (it assumes it is a URL)
            # the sequence id can be quite long, so only use the first token and
            # search and replace ':' with '-' just in case
            label = label.split()[0].replace(':', '-')

            query_fname = "%s%s" % (self.query_prefix, label)

            if query_fname in self.queries :
                self.info("skipped %s (already done)" % query_fname)
                continue

            self._write_query(query_fname, label, seq)
            self.queries.add(query_fname)

            yield os.path.join(self.dir, query_fname)

