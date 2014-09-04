from zipfile import ZipFile
import datetime
import tempfile
import os
import logging

import glutton
from glutton.xmlparser import GluttonXMLManifest, GluttonXMLData
from glutton.ensembl import EnsemblDownloader
from glutton.genefamily import ensembl_to_glutton_internal, Gene, GeneFamily
from glutton.utils import tmpfile

class GluttonImportantFileNotFoundError(Exception) :
    def __init__(self, msg) :
        super(Exception, self).__init__(msg)

metadata_keys = [
    'glutton-version',
    'program-name',
    'program-version',
    'species-name',
    'species-release',
    'download-time',
    'data-file',
    'mapping-file',
    'is-complete'
  ]

class GluttonDB(object) :
    def __init__(self, fname=None) :
        self.metadata = {}
        self.gene_families = None
        self.g_gf_mapping = None
        self.gf_file_mapping = None

        self.log = logging.getLogger('glutton')
        
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.log.addHandler(ch)

        if fname :
            self._read(fname)

    # read manifest to get data and mapping files
    # read data file to get seqID to gluttonID mapping
    # read mapping file to get gluttonID to file name mapping 
    def _read(self, fname) :
        z = ZipFile(fname, 'r')
    
        def _err(msg) :
            z.close()
            raise GluttonImportantFileNotFoundError(msg)
    
        # without the manifest all is lost
        # we need this to get the names of the other
        # XML files
        if 'manifest.xml' not in z.namelist() :
            _err('manifest not found in %s' % fname)

        gm = GluttonXMLManifest()
        self.metadata = gm.read(z.open('manifest.xml'))
        
        # the data file is the raw data grouped into gene families
        # when we do a local alignment we need to get the gene id
        # of the best hit and find out which gene family it belongs to 
        if self.metadata['data-file'] not in z.namelist() :
            _err('data file (%s) not found in %s' % (self.metadata['data-file'], fname))

        gd = GluttonXMLData()
        self.gene_families = gd.read(z.open(self.metadata['data-file']))
        self.g_gf_mapping = self._create_gene_genefamily_mapping(self.gene_families)


        # XXX there might be a file mapping between gene family ids
        # and the files output by prank, but these may not exist yet

        z.close()

    def _create_gene_genefamily_mapping(self, families) :
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

    def _write(self, fname) :
        assert self._valid_manifest(self.metadata)

        z = ZipFile(fname, 'w')

        self._write_to_archive(GluttonXMLManifest(), self.metadata, z, 'manifest.xml')
        self._write_to_archive(GluttonXMLData(),     self.data,     z, self.metadata['data-file'])

        z.close()

    def _default_datafile(self, species, release) :
        return "%s_%d_data.xml" % (species, release)

    def _default_mappingfile(self, species, release) :
        return "%s_%d_mapping.xml" % (species, release)

    def build(self, fname, program, species, release=None) :

        if not os.path.exists(fname) :
            self.log.info("'%s' does not exist, creating..." % fname)
            self._initialise_db(fname, program, species, release)

        self._read(fname)

        # continue aligning things

    def _initialise_db(self, fname, program, species, release=None) :
        e = EnsemblDownloader()
        
        release = release if release else e.get_latest_release(species)

        self.data = ensembl_to_glutton_internal(e.download(species, release))

        self.metadata['download-time'] = datetime.datetime.today()

        self.metadata['glutton-version']    = glutton.__version__
        self.metadata['program-name']       = program.name
        self.metadata['program-version']    = program.version
        self.metadata['species-name']       = species
        self.metadata['species-release']    = release
        self.metadata['is-complete']        = False

        self.metadata['data-file']          = self._default_datafile(species, release)
        self.metadata['mapping-file']       = self._default_mappingfile(species, release)
        
        self._write(fname)

if __name__ == '__main__' :

    from collections import namedtuple

    Program = namedtuple('Program', 'name version')

    gdb = GluttonDB()
    gdb.build('tc23.glt', Program('prank', 1.1), 'tribolium_castaneum')
    
    #gdb = GluttonDB('tc23.glt')

