import sys
import os
import re
import string
import hashlib
import threading
import time

from os.path import join, exists
from glob import glob

from cogent.parse.fasta import MinimalFastaParser

from lib.base import Base


class ManifestError(Exception) :
    pass

class Manifest(Base) :
    manifest_filename = 'manifest'
    done_filename = 'done'

    def __init__(self, opt, dir, prefix, db_name, skip_checks=False) :
        super(Manifest, self).__init__(opt)

        self.dir = dir
        self.name = join(self.dir, self.manifest_filename)
        self.db_name = db_name

        self.manifest_pat = re.compile("^(.+) ([%s]{32})$" % (string.ascii_lowercase + string.digits)) # filename + md5
        self.family_pat = re.compile("^%s[%s]{6}$" % (prefix, string.ascii_letters + string.digits + "_"))

        self.lock = threading.Lock()

        self.genes = set()
        self.realign = []

        try :
            if not skip_checks :
                self.genes, self.realign = self._validate()
            else :
                f2md5 = self._load()
                self.genes = set(f2md5.keys())

        except ManifestError, me :
            self._create_manifest()

    def _create_manifest(self) :
        try :
            open(self.name, 'w').close()

        except IOError, ioe :
            raise ManifestError("could not create manifest file in %s: %s" % (self.dir, ioe.strerror))

    def is_empty(self) :
        return len(self.genes) == 0

    # would be nice to set a flag, but then the checks can be skipped
    # and then something else would have to be done
    def is_complete(self) :
        return exists(join(self.dir, self.done_filename))

    def complete(self) :
        self.append_to_manifest(self.done_filename, time.ctime(), create=True)

    def get_genes(self) :
        return self.genes

    def realignments_required(self) :
        return len(self.realign) != 0

    def get_realignments(self) :
        return [os.path.join(self.dir, f) for f in self.realign]

    def append_to_manifest(self, fname, fcontents, create=True) :
        try :
            self._append_to_manifest(fname, fcontents, create)
            return True

        except IOError, ioe :
            self.error(str(ioe))
            sys.exit(1)

    def _append_to_manifest(self, fname, fcontents, create) : 
        if fcontents :
            manifestdata = self._md5(fcontents)
        else :
            manifestdata = self._file_md5(join(self.dir, fname))

        self.lock.acquire()

        f = open(self.name, 'a')
        print >> f, "%s %s" % (os.path.basename(fname), manifestdata)
        os.fsync(f)
        f.close()

        self.lock.release()

        # write the actual file, we don't care if this gets interrupted
        # because then it will be discovered by recalculating the hash
        if create :
            self._create_file(os.path.basename(fname), fcontents, self.dir)

    def _load(self) :
        tmp = {}

        linenum = 0

        try :
            f = open(self.name)

        except IOError, ioe :
            raise ManifestError("could not load %s: %s" % (self.name, ioe.strerror))

        self.info("verifying manifest file...")

        for line in f : 
            linenum += 1  
            line = line.strip()

            if line == '' :
                continue  

            result = self.manifest_pat.match(line)

            if not result :
                self.warn("line %d in %s appears to be corrupt..." % (linenum, f.name))
                continue

            fname = result.group(1)
            md5val = result.group(2)

            tmp[fname] = md5val

        f.close()

        return tmp

    def _check_hashes(self, f2md5) :
        tmp = []
        
        count = 0
        total = float(len(f2md5)) / 100

        for fname in f2md5 :
            if not self.family_pat.match(fname) :
                count += 1
                self.progress('validating gene families', count / total)
                continue  

            try :
                if self._file_md5(join(self.dir, fname)) == f2md5[fname] :
                    if fname.startswith("paralog_") :
                        tmp.append(fname)
                else :
                    self.warn("md5 for %s is bad" % fname)

            except IOError, ioe :
                self.warn("%s not found!" % fname)

            count += 1
            self.progress('validating gene families', count / total)
        
        self.info('validating gene families done.')
        return tmp

    def _check_alignments(self, f2md5, gene_families) :
        aligned = []
        realign = []

        alignment_suffixes = [".1.dnd", ".2.dnd", ".nuc.1.fas", ".nuc.2.fas", ".pep.1.fas", ".pep.2.fas"]

        count = 0
        total = float(len(gene_families)) / 100

        for fname in gene_families :
            count += 1
            self.progress('validating alignments', count / total)

            if self._count_seqs(fname) < 2 :
                continue

            alignment_files = [fname+suffix for suffix in alignment_suffixes]

            for align_fname in alignment_files :
                try :
                    if self._file_md5(join(self.dir, align_fname)) != f2md5[align_fname] : # bad md5
                        self.warn("md5 for %s is bad" % align_fname)
                        realign.append(fname)
                        break

                except IOError, ioe : # no file to open
                    self.warn("alignment files missing for %s" % fname)
                    realign.append(fname)
                    break

                except KeyError, ke : # no mention in the manifest
                    self.warn("md5 for %s not present in manifest" % align_fname)
                    realign.append(fname)
                    break
            else :
                aligned += alignment_files

        self.info('validating alignments done.')

        return aligned, realign

    def _rewrite(self, f2md5, filenames) :
        f = open(self.name + '.tmp', 'w')

        for fname in filenames :
            print >> f, "%s %s" % (fname, f2md5[fname])

        f.close()

        os.rename(f.name, self.name)

    def _cleanup(self, good_files) :
        for fname in glob(join(self.dir, '*')) :
            basename = os.path.basename(fname)

            if basename == self.manifest_filename :
                continue

            if basename == 'done' :
                continue

            if basename.startswith(self.db_name) :
                continue

            if basename not in good_files :
                self.warn("removing %s" % basename)
                os.remove(fname)

    def _build_genelist(self, good_files) :
        tmp = set()
        
        for fname in good_files :
            if self.family_pat.match(fname) :
                tmp.update(self._read_genenames(join(self.dir, fname)))
        
        return tmp

    def _read_genenames(self, fname) :
        f = open(fname)
        tmp = set()

        for line in f :
            if line.startswith('>') :
                tmp.add(line.strip()[1:])

        f.close()

        return tmp

    def _file_md5(self, fname) :
        assert os.path.isabs(fname)
        return self._md5(self._contents(fname))

    def _md5(self, s) :
        return hashlib.md5(s).hexdigest()

    def _count_seqs(self, fname) :
        count = 0

        for label,seq in MinimalFastaParser(open(join(self.dir, fname))) :
            count += 1

        return count

    def _validate(self) :
        # load manifest into dict
        f2md5 = self._load()

        # check that md5 hashes match gene families and then
        # alignments of those families that passed
        gene_families = self._check_hashes(f2md5)
        aligned,realign = self._check_alignments(f2md5, gene_families)

        # rewrite manifest using only files that passed md5
        good_files = gene_families + aligned
        self._rewrite(f2md5, good_files)

        # clean up other files
        self._cleanup(good_files)

        genes = self._build_genelist(good_files)

        self.info("manifest validation complete (%d genes in %d families)" % (len(genes), len(gene_families)))

        return genes, realign

