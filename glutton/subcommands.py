from sys import stderr, exit
import os
import signal

from glutton.utils import get_log
from glutton.ensembl import EnsemblDownloader, SQLQueryError
from glutton.table import pretty_print_table
from glutton.db import GluttonDB
from glutton.aligner import Aligner
from glutton.scaffolder import Scaffolder


def list_command(args) :
    log = get_log()
    e = EnsemblDownloader()

    suppression_defaults = {
            'ensembl'   : 60,
            'metazoa'   : 18,
            'protists'  : 18,
            'plants'    : 18,
            'bacteria'  : 18,
        }

    if not args.suppress :
        args.suppress = suppression_defaults[args.database_name]

    log.info("listing species in %s database" % args.database_name)
    log.info("suppressing releases prior to %d" % args.suppress)

    try :
        pretty_print_table(
            ('Species name', 'Releases available'), 
            e.get_all_species(db=args.database_name, suppress=args.suppress))

    except SQLQueryError, sql :
        log.fatal(str(sql))
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

    gdb.build(args.output, args.species, args.release)

    log.info("built database %s" % gdb.filename)
    
    return 0

def check_command(args) :
    if GluttonDB(args.gltfile).sanity_check() :
        print "OK!"
        return 0

    print "FAIL!"
    return 1

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

