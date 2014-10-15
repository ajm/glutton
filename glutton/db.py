from zipfile import ZipFile, ZIP_DEFLATED
from sys import exit
import datetime
import tempfile
import os
import logging
import threading
import hashlib

import glutton
from glutton.xmlparser import GluttonXMLManifest, GluttonXMLData, GluttonXMLMapping
from glutton.ensembl import EnsemblDownloader
from glutton.genefamily import ensembl_to_glutton_internal, Gene, GeneFamily, read_alignment_as_genefamily
from glutton.utils import tmpfile, get_log
from glutton.queue import WorkQueue
from glutton.job import PrankJob
from glutton.prank import Prank


class GluttonImportantFileNotFoundError(Exception) :
    pass    

class GluttonDBFileError(Exception) :
    pass

class GluttonDBError(Exception) :
    pass

metadata_keys = [
    'glutton-version',
    'program-name',
    'program-version',
    'species-name',
    'species-release',
    'download-time',
    'data-file',
#    'mapping-file',
#    'is-complete'
  ]

class GluttonDB(object) :
    def __init__(self, fname=None) :
        self.fname = fname
        self.compression = ZIP_DEFLATED
        self.metadata = None
        self.data = None        # dict of famid -> GeneFamily obj (list of Genes)
#        self.mapping = None     # dict of famid -> alignment file name
        self.seq2famid = None   # dict of geneid -> famid
        self.q = None           # work queue if needed
        self.dirty = False
        self.lock = threading.Lock()

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
    def filename(self) :
        return self.fname

    @property
    def checksum(self) :
        m = hashlib.md5()
        f = open(self.fname)

        for block in f.read(1024) :
            m.update(block)

        f.close()

        return m.hexdigest()

    def stop(self) :
        if self.q :
            self.q.stop()

    def flush(self) :
        if self.dirty :
            self._write()

    # read manifest to get data and mapping files
    # read data file to get seqID to gluttonID mapping
    # read mapping file to get gluttonID to file name mapping 
    def _read(self) :
        z = ZipFile(self.fname, 'r', compression=self.compression)
    
        def _err(msg) :
            z.close()
            raise GluttonImportantFileNotFoundError(msg)
    
        # without the manifest all is lost
        # we need this to get the names of the other
        # XML files
        if 'manifest.xml' not in z.namelist() :
            _err('manifest not found in %s' % self.fname)

        gm = GluttonXMLManifest()
        self.metadata = gm.read(z.open('manifest.xml'))
        
        self.log.info("read manifest - created on %s using glutton version %.1f" % \
            (self.metadata['download-time'].strftime('%d/%m/%y at %H:%M:%S'), \
             self.metadata['glutton-version']))

        # the data file is the raw data grouped into gene families
        # when we do a local alignment we need to get the gene id
        # of the best hit and find out which gene family it belongs to 
        if self.metadata['data-file'] not in z.namelist() :
            _err('data file (%s) not found in %s' % (self.metadata['data-file'], self.fname))

        gd = GluttonXMLData()
        self.data = gd.read(z.open(self.metadata['data-file']))
        self.seq2famid = self._create_lookup_table(self.data)

        self.log.info("read %d gene families (%d genes)" % (len(self.data), len(self.seq2famid)))

        # we need the file mapping between gene family ids
        # and the files output by prank
#        if self.metadata['mapping-file'] not in z.namelist() :
#            _err('map file (%s) not found in %s' % (self.metadata['mapping-file'], self.fname))
#
#        gx = GluttonXMLMapping()
#        self.mapping = gx.read(z.open(self.metadata['mapping-file']))
#
#        self.log.info("read %d gene family to file mappings" % len(self.mapping))

        # XXX
#        print "\n\n"
#        for i in z.namelist() :
#            print i
#        print "\n\n"

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

    def _write_to_archive(self, writer, data, zfile, zname) :
        fname = tmpfile()

        f = open(fname, 'w')
        writer.write(f, data)
        f.close()

        zfile.write(fname, arcname=zname)
        os.remove(fname)

    def _write(self) :
        assert self._valid_manifest(self.metadata)

        z = ZipFile(self.fname, 'a', compression=self.compression)

        self._write_to_archive(GluttonXMLManifest(), self.metadata, z, 'manifest.xml')
        self._write_to_archive(GluttonXMLData(),     self.data,     z, self.metadata['data-file'])
#        self._write_to_archive(GluttonXMLMapping(),  self.mapping,  z, self.metadata['mapping-file'])

        z.close()

        self.dirty = False

    def _default_datafile(self, species, release) :
        return "%s_%d_data.xml" % (species, release)

    def _default_mappingfile(self, species, release) :
        return "%s_%d_mapping.xml" % (species, release)

    def build(self, fname, species, release=None) :

        self.fname = fname

        if not os.path.exists(self.fname) :
            self.log.info("'%s' does not exist, creating..." % self.fname)
            self._initialise_db(species, release)
        else :
            self.log.info("'%s' exists, resuming..." % self.fname)

        self._read()

        if self.is_complete() : 
            self.log.info("'%s' is already complete!" % self.fname)
            return

        self._perform_alignments()
        self._write()

        self.log.info("finished building %s/%s" % \
            (self.metadata['species-name'], self.metadata['species-release']))

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

        if not self.q :
            self.q = WorkQueue()

        for i in unaligned :
            self.q.enqueue(PrankJob(self.job_callback, self.data[i]))

        self.log.info("waiting for job queue to drain...")

        self.q.join()

    def _initialise_db(self, species, release=None) :
        e = EnsemblDownloader()

        if not release :
            self.log.info("release not provided, finding out latest release...")
            release = e.get_latest_release(species)
        
        # default name if it was not defined
        if not self.fname :
            self.fname = "%s_%d.glt" % (species, release)

        self.log.info("downloading %s/%d" % (species, release))
        
        self.data = ensembl_to_glutton_internal(e.download(species, release))

        self.log.info("writing sequences and manifest to %s ..." % (self.fname))

        self.metadata = {}
        self.metadata['download-time'] = datetime.datetime.today()

        # glutton metadata
        self.metadata['glutton-version']    = glutton.__version__
        self.metadata['program-name']       = Prank().name
        self.metadata['program-version']    = Prank().version
        self.metadata['species-name']       = species
        self.metadata['species-release']    = release

        # other xml files
        self.metadata['data-file']          = self._default_datafile(species, release)
#        self.metadata['mapping-file']       = self._default_mappingfile(species, release)
#
#        # there is no work to do for families with only a single
#        # gene, so make a note in the file mapping index        
#        self.mapping = {}
#
#        for i in self.data :
#            if len(self.data[i]) < 2 :
#                self.mapping[i] = 'none'

        self.dirty = True
        self._write()

    def job_callback(self, job) :
        self.log.debug("callback from %s - %s" % (str(job),job.input.id))
        self.log.debug("alignment file = %s" % job.alignment)
        self.log.debug("tree file = %s" % job.tree)

        if job.fail() :
            return

        # get a lock on the db
        self.lock.acquire()
        self.dirty = True

        # write alignment file
        # write tree file
        z = ZipFile(self.fname, 'a', compression=self.compression)
        z.write(job.alignment, arcname = self._famid_to_alignment(job.input.id))
        z.write(job.tree,      arcname = self._famid_to_tree(job.input.id))
        z.close()

        # release lock
        self.lock.release()

    def extract_all(self) :
        fname = tmpfile()

        with open(fname, 'w') as f :
            for gf in self.data :
                for g in self.data[gf] :
                    print >> f, g.format('fasta').rstrip()

        return fname

    def get_familyid_from_geneid(self, geneid) :
        return self.seq2famid[geneid]

    def get_genefamily(self, famid) :
        return self.data[famid]

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

    def sanity_check(self, suppress_errmsg=False) :
        z = ZipFile(self.fname, 'r')
        listing = z.namelist()
        z.close()

        insane = False

        for famid in self.data :
            gf_size = len(self.data[famid])

            # err
            if gf_size == 0 :
                if not suppress_errmsg :
                    self.log.error("%s contains no genes!" % famid)
                
                insane = True
                continue

            # no extra files
            elif gf_size == 1 :
                continue

            # .align + .tree expected
            else :
                alignment_fname = self._famid_to_alignment(famid)
                tree_fname = self._famid_to_tree(famid)
                
                for fname in [ self._famid_to_alignment(famid), self._famid_to_tree(famid) ] :
                    if fname not in listing :
                        if not suppress_errmsg :
                            self.log.error("%s not found" % fname)
                        
                        insane = True

        if insane :
            if not suppress_errmsg :
                self.log.error("%s failed check" % self.fname)
                exit(1)

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

