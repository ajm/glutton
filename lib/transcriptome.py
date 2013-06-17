import sys
import tempfile
import os

from os.path import join
from glob import glob

from lib.base import Base
from lib.manifest import Manifest, ManifestError
from lib.ensembl import EnsemblDownloader
from lib.job import PrankJob, PaganJob
from lib.exonerate import ExonerateDBBuilder, ExonerateServer
from lib.genefamily import GeneFamily
from lib.query import QueryManager

class Transcriptome(Base) :
    file_prefix = 'paralog_'
    done_filename = 'done'
    db_name = 'db'

    def __init__(self, opt, queue) :
        super(Transcriptome, self).__init__(opt)

        self.q = queue

        self.species = opt['species']
        self.release = opt['release']

        self.dir = self._get_dir()

        self.cancel = False

        self.exonerate = None
        self.gene2file = None

        try :
            self.manifest = Manifest(opt, self.dir, self.file_prefix, self.db_name, skip_checks=True)
            
        except ManifestError, me :
            self.error(str(me))
            sys.exit(1)

        for fname in self.manifest.get_realignments() :
            self.q.enqueue(
                        PrankJob(
                            self.opt,
                            self.manifest,
                            fname
                            )
                        )
        
        # the build may seem complete, but there are alignments missing
        if self._build_complete() :
            if self.manifest.realignments_required() :
                self.warn("download appears to be complete, but realignments were required, " + \
                          "waiting for these to finish before proceeding...")
                self.q.join()

    def __del__(self) :
        if self.exonerate :
            self.exonerate.stop()

    def is_cancelled(self) :
        return self.cancel

    def stop(self) :
        self.info("stopping...")
        self.cancel = True

    def _build_complete(self) :
        return os.path.exists(os.path.join(self.dir, self.done_filename))

    def _build_exonerate_db(self) :
        exonerate = ExonerateDBBuilder(self.opt, self.dir, self.db_name)
        exonerate.build(self)

        for filename in [self.db_name + '.' + ext for ext in ['fa', 'esd', 'esi']] :
            self.manifest.append_to_manifest(filename, None, create=False)

    def build(self) :
        if self._build_complete() and not self.force :
            return

        ensembl = EnsemblDownloader(self.opt)
        ensembl.set_already_downloaded(self.manifest.get_genes())

        for gene_family in ensembl :
            fname = self._random_filename(self.file_prefix, self.dir)
            self.manifest.append_to_manifest(fname, str(gene_family), create=True)

            if self.cancel :
                break

            if len(gene_family) > 1 :
                self.q.enqueue(
                        PrankJob(
                            self.opt, 
                            self.manifest, 
                            fname
                            )
                        )

        self.q.join()
        if self.cancel :
            return

        self._build_exonerate_db()

        # download + alignments are complete, so write 'done' file
        self.manifest.append_to_manifest(self.done_filename, '', create=True)
        self.info("download complete!")

    # TODO check if db was build on present OS
    #      if not rebuild 
    def fix(self) :
        raise NotImplementedError

    def align(self, contig_fname, contig_outdir, min_length=200) :
        self._check_dir(contig_outdir, create=True)
        
        if not self.exonerate :
            self.exonerate = ExonerateServer(self.opt, self.dir, self.db_name)

        if not self.exonerate.started() :
            self.exonerate.start()

        qm = QueryManager(self.opt, contig_fname, min_length)

        count = 0
        for fname in qm :
            self.q.enqueue(
                    PaganJob(
                        self.opt,
                        self,
                        fname,
                        contig_outdir
                        )
                    )

            count += 1
            if count == 10 :
                break

            if self.cancel : # this does nothing because all this part does is load the q
                break

        self.q.join()
        self.exonerate.stop()
        self.info("alignments complete")

    def query(self, query_fname) :
        return self._gene_to_file(self.exonerate.query(query_fname))

    def _gene_to_file(self, gene_name) :
        if not self.gene2file :
            self.gene2file, self.file2gene = self._build_file_mapping()

        fname = self.gene2file[gene_name]
        return fname, len(self.file2gene[fname])

    def _build_file_mapping(self) :
        g2f = {}
        f2g = {}

        for fname in self._genefamilies_files() :
            gf = GeneFamily(fname)
            f2g[fname] = gf
            for g in gf :
                g2f[g] = fname # XXX same gene in multiple files? XXX
        
        return g2f,f2g

    def _get_dir(self) :
        dirname = join(self.workingdir, str(self.release), self.species)
        self._check_dir(self.tmpdir, create=True)
        self._check_dir(dirname, create=True)
        return dirname

    def _genefamilies(self) :
        for i in self._genefamilies_files() :
            yield GeneFamily(i)

    def _genefamilies_files(self) :
        for i in glob(os.path.join(self.dir, self.file_prefix + '*')) :
            if len(os.path.basename(i)) == (len(self.file_prefix) + 6) : 
                yield i
    
    def __iter__(self) :
        return self._genefamilies()

    def __str__(self) :
        return "%s: %s-%d" % (type(self).__name__, self.species, self.release)

