from sys import stderr, exit
import os
import signal

from glutton.utils import get_log
from glutton.ensembl_downloader import EnsemblDownloader
from glutton.table import pretty_print_table
from glutton.db import GluttonDB
from glutton.aligner import Aligner
from glutton.scaffolder import Scaffolder


def list_command(args) :
    log = get_log()
    e = EnsemblDownloader()

    suppression_defaults = {
            'ensembl'   : 70, # errors with 69 and older (missing columns)
            'metazoa'   : 17,
            'protists'  : 17,
            'plants'    : 17,
            'bacteria'  : 17,
        }

    if not args.suppress :
        args.suppress = suppression_defaults[args.database_name]

    log.info("listing species in %s database" % args.database_name)
    log.info("suppressing releases prior to %d" % args.suppress)

    try :
        pretty_print_table(
            ('Species name', 'Releases available'), 
            e.get_all_species(db=args.database_name, 
                              suppress=args.suppress))

    except EnsemblDownloadError, ede :
        log.fatal(ede.message)
        exit(1)

    return 0

def build_command(args) :
    log = get_log()
    gdb = GluttonDB()

    def _cleanup(signal, frame) :
        print >> stderr, "Killed by user, cleaning up..."
        gdb.stop()
        print >> stderr, "done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    gdb.build(args.output, 
              args.species, 
              args.release, 
              not args.protein, 
              args.download_only)

    log.info("built database %s" % gdb.filename)
    
    return 0

def check_command(args) :
    return 0 if GluttonDB(args.gltfile).sanity_check(human_readable_summary=True, suppress_errmsg=True) else 1        

def align_command(args) :
    gdb = GluttonDB(args.reference)
    align = Aligner(gdb, 
                    args.contigs, 
                    args.alignments, 
                    args.identity, 
                    args.length, 
                    args.batch_size)

    def _cleanup(signal, frame) :
        print >> stderr, "Killed by user, cleaning up..."
        align.stop()
        print >> stderr, "clean up done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    align.align()

    return 0

def scaffold_command(args) :
    gdb = GluttonDB(args.reference)
    scaf = Scaffolder(gdb,
                      args.contigs, 
                      args.alignments, 
                      args.scaffolds,
                      args.identity)

    def _cleanup(signal, frame) :
        print >> stderr, "Killed by user, cleaning up..."
        scaf.stop()
        print >> stderr, "clean up done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    scaf.scaffold()

    return 0

