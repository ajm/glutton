from sys import stderr, exit
import os
import signal

from glutton.utils import get_log
from glutton.ensembl import EnsemblDownloader
from glutton.table import pretty_print_table
from glutton.db import GluttonDB


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

    pretty_print_table(
        ('Species name', 'Releases available'), 
        e.get_all_species(db=args.database_name, suppress=args.suppress))

    return 0

def build_command(args) :
    log = get_log()
    gdb = GluttonDB()

    def _cleanup(signal, frame) :
        print >> stderr, "Killed by user, cleaning up...",
        gdb.stop()
        print >> stderr, "done"
        os._exit(0)

    signal.signal(signal.SIGINT, _cleanup)

    gdb.build(args.output, args.species, args.release)

    log.info("built database %s" % gdb.filename)
    
    return 0

def check_command(args) :
    return 0 if GluttonDB(args.gltfile).sanity_check() else 1

def align_command(args) :
    pass

def scaffold_command(args) :
    pass

