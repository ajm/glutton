from sys import stderr, exit
import os
import signal

from glutton.utils import get_log
from glutton.ensembl_downloader import EnsemblDownloader, EnsemblDownloadError
from glutton.table import pretty_print_table
from glutton.db import GluttonDB, GluttonDBBuildError
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
            'fungi'     : 17,
        }

    if not args.suppress :
        args.suppress = suppression_defaults[args.database_name]

    log.info("listing species in %s database" % args.database_name)
    log.info("suppressing releases prior to %d" % args.suppress)

    try :
        pretty_print_table(
            ('Species name', 'Releases'), 
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
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    try :
        gdb.build(args.output, 
              args.species, 
              args.release, 
              args.database_name,
              True, #not args.protein, 
              args.download_only)

    except GluttonDBBuildError, nmgfe :
        log.fatal(nmgfe.message)
        exit(1)

    log.info("built database %s" % gdb.filename)
    
    return 0

def check_command(args) :
    return 0 if GluttonDB(args.gltfile).sanity_check(human_readable_summary=True, suppress_errmsg=True, show_all=args.show) else 1        

def align_command(args) :
    contigs = zip(args.contigs, args.label, args.species, args.bam) if args.contigs else []

    align = Aligner(args.alignments,
                    args.reference, 
                    contigs, 
                    args.identity, 
                    args.length,
                    args.hitidentity,
                    args.hitlength,
                    args.evalue,
                    args.batchsize)

    def _cleanup(signal, frame) :
        print >> stderr, "Killed by user, cleaning up..."
        align.stop()
        print >> stderr, "clean up done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    align.align()

    return 0

def scaffold_command(args) :
    contigs = zip(args.contigs, args.label, args.species, args.bam) if args.contigs else []

    scaf = Scaffolder(args.alignments, 
                      args.reference,
                      contigs,
                      args.output)

    def _cleanup(signal, frame) :
        print >> stderr, "Killed by user, cleaning up..."
        scaf.stop()
        print >> stderr, "clean up done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    scaf.scaffold()

    return 0

