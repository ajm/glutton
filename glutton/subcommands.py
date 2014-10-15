from glutton.utils import get_log
from glutton.ensembl import EnsemblDownloader
from glutton.table import pretty_print_table


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

def build_command() :
    pass

def check_command() :
    pass

def align_command() :
    pass

def scaffold_command() :
    pass

