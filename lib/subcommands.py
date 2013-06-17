import sys
import signal
import os

from lib.ensembl import EnsemblInfo
from lib.local import LocalInfo
from lib.transcriptome import Transcriptome
from lib.queue import WorkQueue


def check_local(opt, fill_in_release=False) :
    return _check_generic(opt, LocalInfo(opt), fill_in_release)

def check_remote(opt, fill_in_release=False) :
    return _check_generic(opt, EnsemblInfo(opt), fill_in_release)

def _check_generic(opt, info, fill_in_release) :
    # check species is species and valid
    if opt['species'] is None :
        print >> sys.stderr, "Error: you must specify a species!"
        return False

    if (opt['release'] is None) and (not fill_in_release) :
        print >> sys.stderr, "Error: you must specify a release!"
        return False

    # check species is valid
    if not info.is_valid_species(opt['species']) :
        print >> sys.stderr, "Error: '%s' is an invalid species" % opt['species']
        return False

    # get latest version
    if opt['release'] is None and fill_in_release :
        opt['release'] = info.get_latest_release(opt['species'])

        if opt['release'] == -1 :
            print >> sys.stderr, "Error: could not find database for \'%s\'" % opt['species']
            return False

        print >> sys.stderr, "Info: release not specified, using release %d" % opt['release']

    else :
        # check release is valid
        if not info.is_valid_release(opt['species'], opt['release']) :
            print >> sys.stderr, "Error: %d is an invalid release for '%s'" % (opt['release'], opt['species'])
            return False

    return True

def build(opt) :
    t = None

    def build_cleanup(signal, frame) :
        if t :
            if not t.is_cancelled() :
                t.stop()
                print >> sys.stderr, "\n[!] Application will shutdown after the next gene family is downloaded or user types '^C' again...\n"
                return

        print >> sys.stderr, "\n[!] Forced exit from user..."
        os._exit(0)

    # list remote targets
    if opt['list'] :
        return list_ensembl(opt, local=False)

    if not check_remote(opt, fill_in_release=True) :
        return 1

    # install clean up
    signal.signal(signal.SIGINT, build_cleanup)

    q = WorkQueue(opt, opt['threads'])
    q.start()

    t = Transcriptome(opt, q)

    if not opt['alignment-only'] :
        t.build()

    q.stop()

    return 0

def fix(opt) :
    if not check_local(opt) :
        return 1

    t = Transcriptome(opt)
    return t.fix()
    
def assemble(opt) :
    raise NotImplementedError

def debug(opt) :
    ei = EnsemblInfo(opt)
    #ei.print_species_table()

    def test(name, rel) :
        print "Is species '%s' valid?" % name,
        print "Yep" if ei.is_valid_species(name) else "Nope" 
        print "Is release %d valid?" % rel,
        print "Yep" if ei.is_valid_release(name, rel) else "Nope"    

    test('pig', 70)
    test('human', 100)
    test('Loxodonta africana', 39)
    test('xxx', 10)   

    return 0

def align(opt) :
    if opt['list'] :
        list_ensembl(opt, remote=False)
        return 0

    if not check_local(opt, fill_in_release=True) :
        return 1

    def align_cleanup(signal, frame) :
        print >> sys.stderr, "\n[!] Forced exit from user..."
        os._exit(0)

    signal.signal(signal.SIGINT, align_cleanup)

    q = WorkQueue(opt, opt['threads'])
    q.start()

    t = Transcriptome(opt, q)
    t.align(opt['contig-file'], opt['contig-outdir'], min_length=opt['contig-minlen'])

    q.stop()

    return 0

def list_ensembl(opt, remote=True, local=True) :
    if remote :
        EnsemblInfo(opt).print_species_table()

    if local :
        LocalInfo(opt).print_species_table()

    return 0

