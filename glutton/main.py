from sys import stderr, argv, exit, version_info

if version_info < (2,7) :
    print >> stderr, "Sorry: requires python version 2.7 or greater (but not python 3.x)\n"
    exit(1)

import argparse

import glutton
import glutton.subcommands

from glutton.utils import tmpdir, set_threads, num_threads, set_tmpdir, set_verbosity, glutton_log_defaults, get_log
from glutton.ensembl import custom_database


commands = {
    'list'      : glutton.subcommands.list_command,
    'build'     : glutton.subcommands.build_command,
    'check'     : glutton.subcommands.check_command,
    'align'     : glutton.subcommands.align_command,
    'scaffold'  : glutton.subcommands.scaffold_command
}

def handle_args(args) :
    fmt = lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog,
                                                              max_help_position=50,
                                                              width=110)
    
    parser = argparse.ArgumentParser(
                prog=glutton.__name__, 
                description='Evolutionary transcriptome scaffolding',
                formatter_class=fmt)
    
    def add_generic_options(par) :
        par.add_argument('-t', '--tmpdir', type=str, default=tmpdir(),
                         help='temporary directory')
        par.add_argument('-c', '--threads', type=int, default=num_threads(),
                         help='number of threads')
        par.add_argument('-v', '--verbose',  action='count', default=0, 
                         help='set verbosity, can be set multiple times e.g.: -vvv')
        par.add_argument('-V', '--version',  action='version', version="%(prog)s " + str(glutton.__version__), 
                         help='show version')

    def add_database_options(par) :
        par.add_argument('--database-host', type=str, 
                         help='specify database hostname')
        par.add_argument('--database-port', type=int, default=3306, 
                         help='specify database port')
        par.add_argument('--database-user', type=str, default='anonymous', 
                         help='specify database username')
        par.add_argument('--database-password', type=str, default="", 
                         help='specify database password')

    
    subparsers = parser.add_subparsers(help='subcommands')

    # list options
    parser_list = subparsers.add_parser('list', formatter_class=fmt,
                             help='list ensembl species command')
    parser_list.add_argument('-v', '--verbose',  action='count', default=0,   
                             help='set verbosity, can be set multiple times e.g.: -vvv')
    parser_list.add_argument('-d', '--database-name', default='ensembl',
                             choices=('ensembl', 'metazoa', 'fungi', 'protists', 'plants', 'bacteria'),
                             help='specify main ensembl database or one of the ensembl-genomes databases')
    parser_list.add_argument('--suppress', type=int,
                            help='suppress older releases (default depends on database)')
    add_database_options(parser_list)

    # build options
    parser_build = subparsers.add_parser('build', formatter_class=fmt,
                              help='build database command')
    parser_build.add_argument('-s', '--species', type=str, required=True,
                              help='ensembl species name (see list command output)')
    parser_build.add_argument('-r', '--release', type=int, 
                              help='ensembl release number (default: the latest release)')
    parser_build.add_argument('-o', '--output',  type=str, 
                              help='output filename (default: SPECIES_RELEASE.glt)')

    add_database_options(parser_build)
    add_generic_options(parser_build)

    # check options
    parser_check = subparsers.add_parser('check', 
                              help='check %s database' % glutton.__name__)
    parser_check.add_argument('gltfile')

    # align options
    parser_align = subparsers.add_parser('align', formatter_class=fmt,
                              help='align contigs against reference transcript database')
    parser_align.add_argument('-g', '--reference', type=str, required=True,
                              help='reference database, (normally a .glt file)')
    parser_align.add_argument('-c', '--contigs', type=str, required=True,
                              help='fasta files containing contigs')
    parser_align.add_argument('-o', '--output', type=str, default='./alignment_results',
                              help='output directory to store alignment files')
    parser_align.add_argument('-i', '--identity', type=float, default=0.7,
                              help='minimum protein identity for local alignment step')
    parser_align.add_argument('-l', '--length', type=int, default=200,
                              help='minimum contig length')
    parser_align.add_argument('-b', '--batch-size', type=int, default=100,
                              help='batch size for number of queries for local alignment')

    add_generic_options(parser_align)

    # scaffold options
    # TODO

    return parser.parse_args(args)

def generic_options(args) :
    # database
    if hasattr(args, 'database_host') and args.database_host :
        custom_database(args.database_host, args.database_port, args.database_user, args.database_password)

    # threads
    if hasattr(args, 'threads') :
        set_threads(args.threads)

    # tmpdir
    if hasattr(args, 'tmpdir') :
        try :
            set_tmpdir(args.tmpdir)
        
        except OSError :
            print >> stderr, "ERROR: %s does not exist..." % args.tmpdir 

    # verbosity
    if hasattr(args, 'verbose') :
        set_verbosity(args.verbose)

    glutton_log_defaults(get_log())

def main() :
    args = handle_args(argv[1:])

    generic_options(args)

    return commands[argv[1]](args)


if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> stderr, "Killed by user"
        exit(1)

