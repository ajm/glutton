from zipfile import ZipFile, ZIP_DEFLATED
from sys import exit, stderr
from os.path import exists
import tempfile
import os
import logging
import threading
import time
import json

import glutton
from glutton.ensembl_downloader import EnsemblDownloader, EnsemblDownloadError, get_ensembl_download_method
from glutton.genefamily import ensembl_to_glutton, glutton_to_json, json_to_glutton, Gene, GeneFamily, read_alignment_as_genefamily
from glutton.utils import tmpfile, get_log, md5
from glutton.queue import WorkQueue
from glutton.job import PrankJob
from glutton.prank import Prank


class GluttonImportantFileNotFoundError(Exception) :
    pass    

class GluttonDBFileError(Exception) :
    pass

class GluttonDBError(Exception) :
    pass

class GluttonDBBuildError(Exception) :
    pass

metadata_keys = [
    'glutton-version',
    'program-name',
    'program-version',
    'species-name',
    'species-release',
    'download-time',
    'data-file',
    'nucleotide',
    'database-name'
  ]

MANIFEST_FNAME = 'manifest.json'

class GluttonDB(object) :
    def __init__(self, fname=None) :
        self.fname       = fname
        self.compression = ZIP_DEFLATED
        self.metadata    = None
        self.data        = None     # dict of famid -> GeneFamily obj (list of Genes)
        self.seq2famid   = None     # dict of geneid -> famid
        self.dirty       = False
        self.lock        = threading.Lock()
        self.complete_jobs = 0
        self.total_jobs = 0

        self.log = get_log()

        if self.fname :
            self._read()

            if not self.is_complete() :
                self.log.warn("%s is not complete!" % self.fname)

    @property
    def species(self) :
        return self.metadata['species-name']

    @property
    def release(self) :
        return self.metadata['species-release']

    @property
    def nucleotide(self) :
        return self.metadata['nucleotide']

    @property
    def download_time(self) :
        return self.metadata['download-time']

    @property
    def version(self) :
        return self.metadata['glutton-version']

    @property
    def database(self) :
        return self.metadata['database-name']

    @property
    def filename(self) :
        return self.fname

    @property
    def checksum(self) :
        return md5(self.fname)

    def stop(self) :
        if hasattr(self, "q") :
            self.q.stop()

    def flush(self) :
        if self.dirty :
            self._write()

    # read manifest to get data and mapping files
    # read data file to get seqID to gluttonID mapping
    # read mapping file to get gluttonID to file name mapping 
    def _read(self) :
        global MANIFEST_FNAME

        z = ZipFile(self.fname, 'r', compression=self.compression)
    
        def _err(msg) :
            z.close()
            raise GluttonImportantFileNotFoundError(msg)
    
        # without the manifest all is lost
        # we need this to get the names of the other
        # XML files
        if MANIFEST_FNAME not in z.namelist() :
            _err('manifest not found in %s' % self.fname)

        self.metadata = json.load(z.open(MANIFEST_FNAME))
        
        self.log.info("read manifest - created on %s using glutton version %.1f" % \
            (time.strftime('%d/%m/%y at %H:%M:%S', time.localtime(self.download_time)), \
             self.version))

        # the data file is the raw data grouped into gene families
        # when we do a local alignment we need to get the gene id
        # of the best hit and find out which gene family it belongs to 
        if self.metadata['data-file'] not in z.namelist() :
            _err('data file (%s) not found in %s' % (self.metadata['data-file'], self.fname))

        self.data = json_to_glutton(json.load(z.open(self.metadata['data-file'])))
        self.seq2famid = self._create_lookup_table(self.data)

        self.log.info("read %d gene families (%d genes)" % (len(self.data), len(self.seq2famid)))

        z.close()

    def _create_lookup_table(self, families) :
        tmp = {}

        for fid in families :
            for gene in families[fid] :
                tmp[gene.id] = fid
        
        return tmp

    def _valid_manifest(self, m) :
        global metadata_keys

        for k in metadata_keys :
            if k not in m :
                return False
        
        return True

    def _write_to_archive(self, data, zfile, zname) :
        fname = tmpfile()

        f = open(fname, 'w')
        f.write(json.dumps(data))
        f.close()

        zfile.write(fname, arcname=zname)
        os.remove(fname)

    def _write(self) :
        global MANIFEST_FNAME

        assert self._valid_manifest(self.metadata)

        z = ZipFile(self.fname, 'a', compression=self.compression)

        self._write_to_archive(self.metadata,               z, MANIFEST_FNAME)
        self._write_to_archive(glutton_to_json(self.data),  z, self.metadata['data-file'])

        z.close()

        self.dirty = False

    def _default_datafile(self, species, release) :
        return "%s_%d_data.json" % (species, release)

    def build(self, fname, species, release=None, database_name='ensembl', nucleotide=False, download_only=False) :
        self.fname = fname

        # if the name is specified and the file exists, then that means the 
        # download already took place and we should get:
        #   - database_name
        #   - release
        # from the metadata
        if self.fname and exists(self.fname) :
            self.log.info("%s exists, resuming..." % self.fname)
       
        else :
            # release not specified
            if not release :
                self.log.info("release not provided, getting latest release...")
                release = EnsemblDownloader().get_latest_release(species, database_name)
                self.log.info("latest release is %d" % release) 

            # default name if it was not defined
            if not self.fname :
                self.fname = "%s_%d_%s_%s.glt" % (species, release, "nuc" if nucleotide else "pep", get_ensembl_download_method())
                self.log.info("database filename not specified, using '%s'" % self.fname)

            # are we resuming or starting fresh?
            if not exists(self.fname) :
                self.log.info("%s does not exist, starting from scratch..." % self.fname)
                self._initialise_db(species, release, database_name, nucleotide)

        # either way, read contents into memory
        self._read()


        # not really necessary, but check that the species from cli and in the file
        # are the same + nucleotide
        if self.species != species :
            self.log.warn("species from CLI (%s) and glutton file (%s) do not match!" % (species, self.species))

        if release and (self.release != release) :
            self.log.warn("release from CLI (%d) and glutton file (%d) do not match!" % (release, self.release))

        if self.nucleotide != nucleotide :
            self.log.warn("nucleotide/protein from CLI (%s) and glutton file (%s) do not match!" % \
                ("nucleotide" if nucleotide else "protein", "nucleotide" if self.nucleotide else "protein"))


        # no work to do
        if self.is_complete() : 
            self.log.info("%s is already complete!" % self.fname)
            return

        # don't do the analysis, just exit
        if download_only :
            self.log.info("download complete")
            return

        # build db
        self._perform_alignments()
        
        # write to disk
        self._write()

        self.log.info("finished building %s/%s" % (self.species, self.release))

    def _get_unaligned_families(self) :
        unaligned = []

        z = ZipFile(self.fname, 'r', compression=self.compression)
        aligned = set([ i.split('.')[0] for i in z.namelist() if i.endswith('.tree') ])
        z.close()

        for i in self.data :
            if (i not in aligned) and (len(self.data[i]) > 1) :
                unaligned.append(i)

        self.log.info("found %d unaligned gene families" % len(unaligned))

        return unaligned

    def _perform_alignments(self) :
        unaligned = self._get_unaligned_families()

        if not hasattr(self, "q") :
            self.q = WorkQueue()

        self.total_jobs = len(unaligned)
        self.complete_jobs = -1
        self._progress()

        for i in unaligned :
            self.q.enqueue(PrankJob(self.job_callback, self.data[i]))

        self.log.debug("waiting for job queue to drain...")

        self.q.join()

    def _initialise_db(self, species, release, database_name, nucleotide) :
        e = EnsemblDownloader()
        self.log.info("downloading %s/%d" % (species, release))
        
        try :
            self.data = ensembl_to_glutton(e.download(species, release, database_name, nucleotide))

        except EnsemblDownloadError, ede :
            self.log.fatal(ede.message)
            exit(1)



        # drosophila melanogastor - nucleotide - ensembl-main
        # contains transcripts, but not gene families
        count = 0
        for famid in self.data :
            if len(self.data[famid]) == 1 :
                count += 1
        
        if count == len(self.data) :
            raise GluttonDBBuildError("downloaded %d gene families composed of a single gene each ('sql' method will do this on some species that do not contain all transcripts (e.g. drosophila_melanogaster in ensembl-main))" % count)



        self.metadata = {}
        self.metadata['download-time']      = time.time()

        # glutton metadata
        self.metadata['glutton-version']    = glutton.__version__
        self.metadata['program-name']       = Prank().name
        self.metadata['program-version']    = Prank().version
        self.metadata['species-name']       = species
        self.metadata['species-release']    = release
        self.metadata['nucleotide']         = nucleotide
        self.metadata['database-name']      = database_name
        self.metadata['download-method']    = get_ensembl_download_method()

        # other xml files
        self.metadata['data-file']          = self._default_datafile(species, release)
        
        self.dirty = True
        self._write()

    def _progress(self) :
        self.complete_jobs += 1

        stderr.write("\rProgress: %d / %d prank alignments " % (self.complete_jobs, self.total_jobs))
        stderr.flush()

        if self.complete_jobs == self.total_jobs :
            print >> stderr, "\ndone!"

    def job_callback(self, job) :
        self.log.debug("callback from %s - %s" % (str(job), job.input.id))
        self.log.debug("alignment file = %s" % job.alignment)
        self.log.debug("tree file = %s" % job.tree)

        # get a lock on the db
        self.lock.acquire()
        self.dirty = True

        self._progress()

        if job.success() :
            # write alignment file
            # write tree file
            z = ZipFile(self.fname, 'a', compression=self.compression)
            z.write(job.alignment, arcname = self._famid_to_alignment(job.input.id))
            z.write(job.tree,      arcname = self._famid_to_tree(job.input.id))
            z.close()
        else :
            if self.q.running :
                self.log.error("Could not align gene family (%s)" % job.input.id)

        # release lock
        self.lock.release()

    # this is only used by the aligner to give localsearch a file containing 
    # protein sequences
    def extract_all(self) :
        fname = tmpfile()

        with open(fname, 'w') as f :
            for gf in self.data :
                for g in self.data[gf] :
                    print >> f, g.format('protein' if self.nucleotide else 'fasta').rstrip()

        return fname

    def get_familyid_from_geneid(self, geneid) :
        return self.seq2famid[geneid]

    def get_genefamily(self, famid) :
        return self.data[famid]

    # this is quite awkward...
    def get_genename_from_geneid(self, geneid) :
        return dict([ (i.id, i.name) for i in self.data[self.seq2famid[geneid]] ])[geneid]

    def _famid_to_alignment(self, famid) :
        return famid + '.align'

    def _famid_to_tree(self, famid) :
        return famid + '.tree'

    def get_alignment(self, famid) :
        # there are two situations
        # 1. there is only a single gene and therefore no alignment file or tree file
        # 2. there is more than one gene and there should be famid.align and famid.tree in the archive
        # (0. the famid is invalid)

        if not self.data.has_key(famid) :
            raise GluttonDBError("genefamily with id '%s' not found" % famid)

        # situation 1
        if len(self.data[famid]) == 1 :
            return self.data[famid]

        # situation 2
        alignment_fname = self._famid_to_alignment(famid)
        tree_fname = self._famid_to_tree(famid)
        
        self.lock.acquire()

        z = ZipFile(self.fname, 'r')
        listing = z.namelist()

        # check that both files exist
        if alignment_fname not in listing :
            raise GluttonDBFileError("'%s' not found" % alignment_fname)

        if tree_fname not in listing :
            raise GluttonDBFileError("'%s' not found" % tree_fname)

        # read contents
        alignment = read_alignment_as_genefamily(z.open(alignment_fname), famid)
        alignment.set_tree(z.read(tree_fname))

        z.close()
        self.lock.release()


        return alignment

    def is_complete(self) :
        return self.sanity_check(suppress_errmsg=True)

    def sanity_check(self, suppress_errmsg=False, human_readable_summary=False, show_all=False) :
        z = ZipFile(self.fname, 'r')
        listing = z.namelist()
        z.close()

        insane = False

        summary = {
            'num_single_genes'  : 0,
            'num_genes' : 0,
            'error_empty'       : 0,
            'error_noalignment' : 0,
            'error_notree'      : 0,
            'error_nofiles'     : 0
        }

        bad_genefamilies = []


        for famid in self.data :
            gf_size = len(self.data[famid])

            summary['num_genes'] += gf_size

            # err
            if gf_size == 0 :
                summary['error_empty'] += 1

                if not suppress_errmsg :
                    self.log.error("%s contains no genes!" % famid)
                
                insane = True
                continue

            # no extra files
            elif gf_size == 1 :
                summary['num_single_genes'] += 1
                continue

            # .align + .tree expected
            else :
                alignment_fname = self._famid_to_alignment(famid)
                tree_fname = self._famid_to_tree(famid)

                if alignment_fname not in listing :
                    insane = True
                    summary['error_noalignment'] += 1
                    if not suppress_errmsg :
                        self.log.error("%s not found" % alignment_fname)

                if tree_fname not in listing :
                    insane = True
                    summary['error_notree'] += 1
                    if not suppress_errmsg :
                        self.log.error("%s not found" % tree_fname)
                
                if (alignment_fname not in listing) and (tree_fname not in listing) :
                    summary['error_nofiles'] += 1

                if (alignment_fname not in listing) or (tree_fname not in listing) :
                    bad_genefamilies.append(famid)

        if insane :
            if not suppress_errmsg :
                self.log.error("%s failed check" % self.fname)
                exit(1)

        if human_readable_summary :
            print ""
            print "Filename:", self.fname
            print ""
            print "Species:", self.species
            print "Release:", self.release
            print "Database:", self.database
            print "Type:", "nucleotide" if self.nucleotide else "protein"
            print "Downloaded:", time.strftime('%d/%m/%y at %H:%M:%S', time.localtime(self.download_time))
            print "Checksum:", self.checksum
            print ""
            print "Number of genes:", summary['num_genes']
            print "Number of gene families:", len(self.data)
            print "Families containing multiple genes:", len(self.data) - summary['num_single_genes']
            print "Missing alignments:", summary['error_noalignment']
            print "Missing phylogenetic trees:", summary['error_notree']
            print ""
            print "FAIL! (to see missing families use --show option)" if insane else "OK!"
            print ""

            if show_all and bad_genefamilies :
                print "Missing gene families:"
                print ""
                for famid in sorted(bad_genefamilies) :
                    print "\t%s\t(%s)" % (famid, ', '.join([ i.name for i in self.data[famid] ]))
                print ""

        return not insane

    def ls(self) :
        z = ZipFile(self.fname, 'r')
        listing = z.namelist()
        z.close()

        print "\n".join(listing)

    def debug(self) :
        x = self.data['genefamily3859']

        print len(x)
        
        for y in x :
            print y.format('fasta').rstrip()

if __name__ == '__main__' :

    import signal
    import sys
    import os

    from glutton.utils import glutton_log_defaults


    glutton_log_defaults(get_log())

    #gdb = GluttonDB()
    gdb = GluttonDB('tc23.glt')

    def _cleanup(signal, frame) :
        print >> sys.stderr, "Killed by user, cleaning up...",
        gdb.stop()
        print >> sys.stderr, "done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    #gdb.build('tc23.glt', 'tribolium_castaneum')
    
    
    #print gdb.extract_all()
    gdb.sanity_check()
    #gdb.debug()

