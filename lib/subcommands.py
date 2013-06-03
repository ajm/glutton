import sys
import signal
import os

from lib.cache import EnsemblInfo, TranscriptCache
from lib.local import LocalInfo

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

    # check species2 if non-None
#    if (opt['species2'] is not None) and (not info.is_valid_species(opt['species2'])) :
#        print >> sys.stderr, "Error: '%s' is an invalid species" % opt['species2']
#        return False

    if opt['release'] is None and fill_in_release :
        # get latest version
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
    cache = None
    
    def build_cleanup(signal, frame) :
        if cache :
            if not cache.stop :
                cache.stop = True
                print >> sys.stderr, "\n[!] Application will shutdown after the next gene family is downloaded or user types '^C' again...\n"
                return

        print >> sys.stderr, "\n[!] Forced exit from user..."
        sys.exit(0)

    # list remote targets
    if opt['list'] :
        return list_ensembl(opt, local=False)

    if not check_remote(opt, fill_in_release=True) :
        return -1

    # install clean up
    signal.signal(signal.SIGINT, build_cleanup)

    cache = TranscriptCache(opt)

    # XXX this is a hack
    # TODO clean up with a proper api to select and deselect different actions
    if opt['alignment-only'] :
        cache.download_complete = True
        cache.join()
        cache.build_exonerate_index()
        print "Info: done!"
    else :
        cache.build()

    return 0

def fix(opt) :
    if not check_local(opt) :
        return -1

    cache = TranscriptCache(opt)
    cache.fix()
    
    print "Info: species: %s, release: %d, OK!" % (opt['species'], opt['release']) 
    return 0

def assemble(opt) :
    print 'assemble'

def debug(opt) :
    ei = EnsemblInfo(opt)
    #ei._get_latest_release_versions2()
    #print ei.get_latest_release2(opt['species'])
    ei.print_species_table()

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
        return -1

    return 0

# TODO list local repositories + status (incomplete, complete)
def list_ensembl(opt, remote=True, local=True) :
    if remote :
        EnsemblInfo(opt).print_species_table()

    if local :
        LocalInfo(opt).print_species_table()

    return 0

