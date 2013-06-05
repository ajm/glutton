import sys
import tempfile

from os.path import join
from glob import glob

from lib.base import Base
from lib.manifest import Manifest, ManifestError
from lib.ensembl import EnsemblDownloader
from lib.job import SaveJob

class Transcriptome(Base) :
    file_prefix = 'paralog_'

    def __init__(self, opt, queue) :
        super(Transcriptome, self).__init__(opt)

        self.options = opt
        self.q = queue

        self.species = opt['species']
        self.release = opt['release']

        self.dir = self._get_dir()

        try :
            self.manifest = Manifest(opt, self.dir, self.file_prefix)
    
        except ManifestError, me :
            self.error(str(me))
            sys.exit(-1)

    def build(self) :
        ensembl = EnsemblDownloader(self.options)

        for gene_family in ensembl :
            self.q.enqueue(
                        SaveJob(
                            self.options, 
                            self.manifest,           # manifest obj
                            self._random_filename(), # save to this file
                            str(gene_family),        # save this contents
                            align=True               # send to prank for alignment
                            )
                        )

        self.q.done()

    def fix(self) :
        pass

    def align(self) :
        pass

    def _get_dir(self) :
        dirname = join(self.workingdir, str(self.release), self.species)
        self._check_dir(self.tmpdir, create=True)
        self._check_dir(dirname, create=True)
        return dirname

    def _random_filename(self) :
        return tempfile.mkstemp(prefix=type(self).file_prefix, dir=self.dir)

    def _genefamilies(self) :
        for i in glob(type(self).file_prefix + '*') :
            if len(basename(i)) == (len(type(self).file_prefix) + 6) : 
                yield i

    def __iter__(self) :
        return self._genefamilies()

    def __str__(self) :
        return "%s: %s-%d" % (type(self).__name__, self.species, self.release)

