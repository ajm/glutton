import collections
import shutil
import sys
import threading

from glutton.db import GluttonDB, GluttonDBError, GluttonDBFileError
from glutton.localsearch import All_vs_all_search
from glutton.utils import tmpfile, num_threads, get_log, rm_f, check_dir, md5
from glutton.queue import WorkQueue
from glutton.job import PaganJob
from glutton.genefamily import Gene, biopy_to_gene, seqlen
from glutton.info import GluttonInformation

from os.path import isfile, basename

from Bio import SeqIO


class Aligner(object) :
    def __init__(self, directory, reference_fname, contig_files, min_length, min_hitidentity, min_hitlength, max_evalue, batch_size, min_alignidentity, min_alignoverlap) :
        self.directory = directory
        self.min_length = min_length # glutton
        self.min_hitidentity = min_hitidentity # blast 
        self.min_hitlength = min_hitlength # blast
        self.max_evalue = max_evalue # blast
        self.min_alignidentity = min_alignidentity # pagan
        self.min_alignoverlap = min_alignoverlap # pagan

        check_dir(self.directory, create=True)

        self.search = All_vs_all_search(batch_size)
        self.cleanup_files = []
        self.q = None

        self.lock = threading.Lock()
        self.complete_jobs = 0
        self.total_jobs = 0

        self.log = get_log()

        self.info = GluttonInformation(directory, reference_fname, contig_files)
        
        self.db = self.info.get_db()    

    def _read_contigs(self) :
        contigs = {}

        for label in self.info.get_sample_ids() :
            accepted = 0
            rejected = 0

            fname = self.info.label_to_contigs(label)

            for r in SeqIO.parse(fname, 'fasta') :
                if len(r) < self.min_length :
                    rejected += 1
                    continue

                qid = self.info.get_query_from_contig(label, r.description)
            
                contigs[qid] = biopy_to_gene(r, qid)
                accepted += 1

            self.log.info("%s: read %d contigs (rejected %d due to length < %d)" % 
                (fname, accepted, rejected, self.min_length))

        return contigs

    def stop(self) :
        self.search.stop()
        self.info.update_query_gene_mapping(self.search.get_intermediate_results())
        
        if self.q :
            self.q.stop()

        rm_f(self.cleanup_files)

        self.info.flush()

    def align(self) :
        self.log.info("starting alignment procedure")

        # convert the names of the contigs to something no program can complain about
        # + filter out the ones that could never have a long enough alignment
        contigs = self._read_contigs()

        pending_contigs = [ contigs[i] for i in self.info.pending_queries() ]

        self.log.info("%d contigs have not been assigned to genes..." % len(pending_contigs))

        # depending on when the program was terminated this step may be complete or partially
        # complete 
        if pending_contigs :
            db_fname = self.db.extract_all()
            self.cleanup_files.append(db_fname)

            # do an all vs all search of contigs vs database of transcripts
            # return a dict of tmp ids with gene ids
            self.info.update_query_gene_mapping(
                self.search.process(
                    db_fname, 
                    pending_contigs,
                    self.db.nucleotide,
                    self.min_hitidentity,
                    self.min_hitlength,
                    self.max_evalue)
                )

            rm_f(db_fname)

        # save intermediate results
        self.info.flush()

        # use the database to convert the mapping from tmp id -> gene
        # to gene family -> list of tmp ids
        genefamily_contig_map = self.info.build_genefamily2contigs()
        
        self.log.info("%d contigs assigned to %d gene families" % 
                (sum([ len(i) for i in genefamily_contig_map.values() ]), len(genefamily_contig_map)))
        self.log.info("(%d have already been run)" % self.info.len_genefamily2filename())

        if self.info.len_genefamily2filename() == len(genefamily_contig_map) :
            self.log.info("alignment already done, exiting early...")
            return
        else :
            self.log.info("starting alignments...")


        # queue all the alignments up using a work queue and pagan
        self.q = WorkQueue()

        self.total_jobs = len(genefamily_contig_map) - self.info.len_genefamily2filename()
        self.complete_jobs = -1
        self._progress()

        for famid in self.sort_keys_by_complexity(genefamily_contig_map) :
            # ignore the jobs that have already been run
            if self.info.in_genefamily2filename(famid) :
                continue

            try :
                # get the alignment and tree from the database
                alignment = self.db.get_alignment(famid)
                tree = alignment.get_tree()

                # get contigs
                job_contigs = [ contigs[i] for i in genefamily_contig_map[famid] ]

                # queue the job
                self.q.enqueue(
                    PaganJob(
                        self.job_callback,
                        job_contigs,
                        famid,
                        alignment,
                        tree,
                        self.min_alignidentity,
                        self.min_alignoverlap)
                    )

                # avoid the split code later in the loop...
                continue

            except GluttonDBError, gde :
                # this means we have never heard of this gene family
                self.log.warn(str(gde))
                continue

            except GluttonDBFileError, gdfe :
                # this means we have heard of the gene family, but the
                # alignment files were missing...
                self.log.warn(str(gdfe))

            # okay, the gene family was not aligned for some reason
            # instead we will split the gene family into constituent genes
            # and handle each one separately...

            self.log.warn("gene family was not aligned, breaking down into separate genes...")
            self.total_jobs += (len(genefamily_contig_map[famid]) - 1)

            for i in genefamily_contig_map[famid] :
                try :
                    geneid = self.info.query_to_gene(i)

                except KeyError : # this should be impossible
                    self.log.warn("no gene assignment for %s" % i)
                    continue

                try :
                    alignment = [ self.db.get_gene(geneid) ]

                except GluttonDBError, gde :
                    self.log.warn(str(gde))
                    continue

                # queue the job
                self.q.enqueue(
                    PaganJob(
                        self.job_callback,
                        [ contigs[i] ],
                        geneid,
                        alignment,
                        None)
                    )


        self.log.debug("waiting for job queue to drain...")
        self.q.join()

        # save all the results again
        self.info.flush()

    def sort_keys_by_complexity(self, d) :
        #return [ k for k,v in sorted(d.items(), reverse=True, key=lambda x : seqlen(x[1])) ]
        return [ k for l,k,v in sorted([ (seqlen(v), k, v) for k,v in d.items() ], reverse=True) ]

    def _progress(self) :
        self.lock.acquire()

        self.complete_jobs += 1

        sys.stderr.write("\rProgress: %d / %d pagan alignments " % (self.complete_jobs, self.total_jobs))
        sys.stderr.flush()

        if self.complete_jobs == self.total_jobs :
            sys.stderr.write("\n")
            sys.stderr.flush()

        self.lock.release()

    def job_callback(self, job) :
        self.log.debug("callback from %s: %s + %s" % (str(job), job.genefamily, ','.join([ i.id for i in job.input ])))
        self.log.debug("protein alignment file = %s" % (job.protein_alignment))

        if self.db.nucleotide :
            self.log.debug("nucleotide alignment file = %s" % (job.nucleotide_alignment))

        self._progress()

        if job.success() :
            dst = tmpfile(directory=self.directory, suffix='.protein')
            dst_base = basename(dst)[:-8]


            self.log.debug("cp %s %s" % (job.protein_alignment, dst))
            shutil.copyfile(job.protein_alignment, dst)
            
            shutil.copyfile(job.query_fname, dst[:-8] + '.queries')
            shutil.copyfile(job.alignment_fname, dst[:-8] + '.reference')
            
            if job.tree_fname :
                shutil.copyfile(job.tree_fname, dst[:-8] + '.tree')

            if self.db.nucleotide :
                self.log.debug("cp %s %s" % (job.nucleotide_alignment, dst[:-8] + '.nucleotide'))
                shutil.copyfile(job.nucleotide_alignment, dst[:-8] + '.nucleotide')
        
            self.info.put_genefamily2filename(job.genefamily, dst_base)
        else :
            self.info.put_genefamily2filename(job.genefamily)

