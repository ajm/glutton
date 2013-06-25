import sys
import tempfile
import os
import subprocess

from os.path import join, basename
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
    db_name = 'db'

    def __init__(self, opt, queue, skip_checks=False) :
        super(Transcriptome, self).__init__(opt)

        self.q = queue

        self.species = opt['species']
        self.release = opt['release']

        self.dir = self._get_dir()

        self.cancel = False

        self.exonerate = None
        self.gene2file = None

        self._init_manifest(skip_checks)

    def _init_manifest(self, skip_validation) :
        try :
            self.manifest = Manifest(self.opt, self.dir, self.file_prefix, self.db_name, skip_checks=skip_validation)

        except ManifestError, me :
            self.error(str(me))
            sys.exit(1)

        if skip_validation and not self.manifest.is_complete() :
            if self.force :
                return

            self.error("halting operation due to incomplete transcriptome database" + 
                       "\n\t- use --force to continue anyway or ..." + 
                       "\n\t- finish building the transcriptome with \"%s build -s '%s' -r %d -d '%s'\"" % \
                               (sys.argv[0], self.species, self.release, self.opt['database']))
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
        if self.manifest.is_complete() :
            if self.manifest.realignments_required() :
                self.warn("download appears to be complete, but realignments were required, " + \
                          "waiting for these to finish before proceeding...")
                self.q.join()

    def is_cancelled(self) :
        return self.cancel

    def stop(self) :
        self.info("stopping...")
        self.cancel = True

    def _build_exonerate_db(self) :
        exonerate = ExonerateDBBuilder(self.opt, self.dir, self.db_name)
        exonerate.build(self)

        for filename in exonerate.filenames() : 
            self.manifest.append_to_manifest(filename, None, create=False)

    def build(self) :
        if self.manifest.is_complete() and not self.force :
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
        self.manifest.complete()
        
        self.info("complete!")

    def fix(self) :
        if not self.manifest.is_complete() and not self.force :
            self.warn("transcriptome for %s:%d is not complete, nothing to fix!" % (self.species, self.release))
            return 1

        self.info('rebuilding exonerate database...')
        self._build_exonerate_db()

        self.info("complete!")
        return 0

    def align(self, contig_fname, contig_outdir, min_length=200) :
        self._check_dir(contig_outdir, create=True)
        
        self.exonerate = ExonerateServer(self.opt, self.dir, self.db_name)
        self.exonerate.start()

        qm = QueryManager(self.opt, contig_fname, contig_outdir, min_length)

        for fname in qm :
            self.q.enqueue(
                    PaganJob(
                        self.opt,
                        self,
                        fname,
                        contig_outdir
                        )
                    )

        self.q.join()
        self.exonerate.stop()

        self.info("complete!")

    def query(self, query_fname) :
        return list(set([ self._gene_to_file(i) for i in self.exonerate.query(query_fname) ]))

        # XXX might be more than one
        #return self._gene_to_file(self.exonerate.query(query_fname))

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
        for i in glob(join(self.dir, self.file_prefix + '*')) :
            if len(os.path.basename(i)) == (len(self.file_prefix) + 6) : 
                yield i
    
    def __iter__(self) :
        return self._genefamilies()

    def __str__(self) :
        return "%s: %s-%d" % (type(self).__name__, self.species, self.release)

    def package(self) :
        outfile_name = join(os.getcwd(), "%s_%d.tgz" % (self.species, self.release))

        if os.path.isfile(outfile_name) :
            self.error("file %s already exists..." % outfile_name)
            return 1

        args = ['tar', 'zcf',
                outfile_name,
                join(str(self.release), self.species)]

        self.info("packing %s ..." % outfile_name)

        try :
            subprocess.check_output(args, 
                                    cwd=self.workingdir, 
                                    stderr=open('/dev/null', 'w'))

        except subprocess.CalledProcessError, cpe :
            self.error(str(cpe))
            return 1

        self.info("complete!")

        return 0
        
    def unpackage(self, package_name) :
        args = ['tar', 'xf',
                package_name]

        self.info("unpacking %s ..." % basename(package_name))

        try :
            subprocess.check_output(args, 
                                    cwd=self.workingdir, 
                                    stderr=open('/dev/null', 'w'))

        except subprocess.CalledProcessError, cpe :
            self.error(str(cpe))
            return 1

        self.info("complete!")

        return 0

