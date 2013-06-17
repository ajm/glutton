import os
from glob import glob

from cogent.parse.fasta import MinimalFastaParser

from lib.base import Base

class QueryManager(Base) :
    query_prefix = 'query_'

    def __init__(self, opt, fname, min_length) :
        super(QueryManager, self).__init__(opt)

        assert os.path.isabs(fname)

        self.fname = fname
        self.min_length = min_length
        self.dir = self.tmpdir

        self.queries = set()

    def __iter__(self) :
        return self._queries()

    def _write_query(self, fname, label, seq) :
        self._create_file(fname, ">%s\n%s" % (label, seq), self.dir)

    # TODO make restartable
    def _queries(self) :
        for f in glob(os.path.join(self.dir, self.query_prefix + '*')) :
            self.queries.add(f)

        for label,seq in MinimalFastaParser(open(self.fname)) :
            if len(seq) < self.min_length :
                continue

            query_fname = "%s%s" % (self.query_prefix, label.replace(' ', '-'))

            if query_fname in self.queries :
                self.info("skipped %s" % query_fname)
                continue

            self._write_query(query_fname, label, seq)
            self.queries.add(query_fname)

            yield os.path.join(self.dir, query_fname)

