from sys import stderr, argv, exit, version_info

if version_info < (2,7) :
    print >> stderr, "Sorry: requires python version 2.7 or greater (but not python 3.x)\n"
    exit(1)

import time
import argparse
import os

import glutton
import glutton.subcommands

from glutton.utils import tmpdir, set_threads, num_threads, set_tmpdir, set_verbosity, setup_logging, get_log, duration_str, check_dir
from glutton.ensembl_sql import custom_database
from glutton.ensembl_downloader import set_ensembl_download_method
from glutton.assembler_output import supported_assemblers


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
    
    parser.add_argument('-V', '--version',  action='version', version="%(prog)s " + str(glutton.__version__),
                     help='show version')

    def add_generic_options(par) :
        par.add_argument('--tmpdir', type=str, default=tmpdir(),
                         help='temporary directory')
        par.add_argument('--threads', type=int, default=num_threads(),
                         help='number of threads')
        par.add_argument('-v', '--verbose',  action='count', default=2, 
                         help='set verbosity, can be set multiple times e.g.: -vvv')

    def add_database_options(par) :
        ensembl_db = ('ensembl', 'metazoa', 'fungi', 'protists', 'plants')

        par.add_argument('-d', '--database-name', default='ensembl', metavar='DB', choices=ensembl_db,
                         help='specify main ensembl database or one of the ensembl-genomes databases, options are %s' % ', '.join(ensembl_db))
        par.add_argument('--database-host', type=str, 
                         help='specify database hostname')
        par.add_argument('--database-port', type=int, default=3306, 
                         help='specify database port')
        par.add_argument('--database-user', type=str, default='anonymous', 
                         help='specify database username')
        par.add_argument('--database-password', type=str, default="", 
                         help='specify database password')

    def add_input_files_options(par) :
        par.add_argument('-c', '--contigs', type=str, action='append',
                         help='fasta file containing contigs')
        par.add_argument('-l', '--label',   type=str, action='append',
                         help='freeform label for file containing contigs (specified with --contigs)')
        par.add_argument('-s', '--species', type=str, action='append',
                         help='species designation for file containing contigs (specified with --contigs)')
        par.add_argument('-b', '--bam', type=str, action='append',
                         help='BAM file of reads mapped to contigs (specified with --contigs)')

    def check_zero_one(v):
        x = float(v)
        if x < 0.0 or x > 1.0 :
            raise argparse.ArgumentTypeError("%s is not in the range [0.0,1.0]" % v)
        return x

    def check_non_negative(v) :
        x = int(v)
        if x < 0 :
            raise argparse.ArgumentTypeError("%s is negative" % v)
        return x
    
    def check_greater_than_zero(v) :
        x = int(v)
        if x < 1 :
            raise argparse.ArgumentTypeError("%s is zero or less" % v)
        return x


    subparsers = parser.add_subparsers(help='subcommands')

    ensembl_methods = ('sql', 'biomart', 'pycogent')

    # list options
    parser_list = subparsers.add_parser('list', formatter_class=fmt,
                             help='list ensembl species')
    parser_list.add_argument('-v', '--verbose',  action='count', default=0,   
                             help='set verbosity, can be set multiple times e.g.: -vvv')
    parser_list.add_argument('-m', '--download-method', default='biomart', metavar='METHOD', choices=ensembl_methods,
                             help='specific download method, options are %s' % ', '.join(ensembl_methods))
    parser_list.add_argument('--suppress', type=int,
                             help='suppress older releases (default depends on database)')
    add_database_options(parser_list)


    # build options
    parser_build = subparsers.add_parser('build', formatter_class=fmt,
                              help='build reference transcript database from ensembl')
    parser_build.add_argument('gltfile', nargs='?')
    parser_build.add_argument('-s', '--species', type=str,
                              help='ensembl species name (see list command output)')
    parser_build.add_argument('-r', '--release', type=int, 
                              help='ensembl release number (default: the latest release)')
    parser_build.add_argument('-o', '--output',  type=str, 
                              help='output filename (default: SPECIES_RELEASE.glt)')
#    parser_build.add_argument('-p', '--protein', action='store_true',
#                              help='download protein sequences instead of CDS sequences')
    parser_build.add_argument('--download-only', action='store_true',
                              help='download sequences and homology information, then exit')
    parser_build.add_argument('-m', '--download-method', default='biomart', metavar='METHOD', choices=ensembl_methods,
                             help='specific download method, options are %s' % ', '.join(ensembl_methods))

    add_database_options(parser_build)
    add_generic_options(parser_build)

    # check options
    parser_check = subparsers.add_parser('check', 
                              help='check %s database for errors/completeness' % glutton.__name__)
    parser_check.add_argument('gltfile')
    parser_check.add_argument('-s', '--show', action='store_true',
                              help='show which gene families have not been aligned')

    
    default_alignment_dir = './alignment_results'
    default_scaffold_dir = './scaffold_results'

    # align options
    parser_align = subparsers.add_parser('align', formatter_class=fmt,
                              help='align contigs against reference transcript database')
    parser_align.add_argument('-g', '--reference', type=str,
                              help='reference database, (normally a .glt file)')
    parser_align.add_argument('-a', '--alignments', type=str, default=default_alignment_dir,
                              help='output directory to store alignment files')
#    parser_align.add_argument(      '--assembler', type=str, default='trinity', metavar='ASSEMBLER', choices=supported_assemblers,
#                              help='assembler used to assemble contigs, options are %s' % ', '.join(supported_assemblers))
    
    parser_align.add_argument('-I', '--hitidentity', type=check_zero_one, default=0.7,
                              help='minimum blast hit identity')
    parser_align.add_argument('-L', '--hitlength', type=check_non_negative, default=100,
                              help='minimum blast hit length')
    parser_align.add_argument('-E', '--evalue', type=float, default=1e-3,
                              help='maximum blast hit E value')

    parser_align.add_argument('-x', '--length', type=check_non_negative, default=200,
                              help='minimum contig length for gene assignment step')
    parser_align.add_argument('-B', '--batchsize', type=check_greater_than_zero, default=100,
                              help='number of queries per batch for local alignment')

    add_input_files_options(parser_align)
    add_generic_options(parser_align)


    # scaffold options
    parser_scaf = subparsers.add_parser('scaffold', formatter_class=fmt,
                             help='scaffold contigs together based on evolutionary alignment results')
    parser_scaf.add_argument('-g', '--reference', type=str,
                             help='reference database, (normally a .glt file)')
    parser_scaf.add_argument('-a', '--alignments', type=str, default=default_alignment_dir,
                             help='directory containing evolutionary alignments')
    parser_scaf.add_argument('-o', '--output', type=str, default=default_scaffold_dir,
                             help='directory to output scaffolded contigs and MSAs')
    parser_scaf.add_argument(      '--assembler', type=str, default='trinity', metavar='ASSEMBLER', choices=supported_assemblers,
                             help='assembler used to assemble contigs, options are %s' % ', '.join(supported_assemblers))
    parser_scaf.add_argument(      '--identity', type=check_zero_one, default=0.5,
                             help='minimum identity in protein space for contig to be included in alignment')
    parser_scaf.add_argument(      '--length', type=check_non_negative, default=100,
                             help='minimum length in nucleotide space for contig to be included in alignment')

    add_input_files_options(parser_scaf)
    add_generic_options(parser_scaf)


    return parser.parse_args(args)

# the reason i have used hasattr is becuase not all the parsers for
# different subcommands have the same options, so they must be tested
# for in order for this to be generic
def generic_options(args) :
    # database
    if hasattr(args, 'database_host') and args.database_host :
        custom_database(args.database_host, 
                        args.database_port, 
                        args.database_user, 
                        args.database_password)

    # threads
    if hasattr(args, 'threads') :
        set_threads(args.threads)

    # tmpdir
    if hasattr(args, 'tmpdir') :
        try :
            check_dir(args.tmpdir, create=True)
            set_tmpdir(args.tmpdir)
        
        except OSError :
            print >> stderr, "ERROR: %s does not exist..." % args.tmpdir 
            exit(1)

    # gltfile
    if hasattr(args, 'gltfile') and args.gltfile :
        if not os.path.isfile(args.gltfile) :
            print >> stderr, "ERROR: %s does not exist..." % args.gltfile
            exit(1)

        args.output = args.gltfile

    if hasattr(args, 'gltfile') and hasattr(args, 'species') :
        if not args.gltfile and not args.species :
            print >> stderr, "ERROR: you must specify either the species or an existing GLT file..."
            exit(1)

    if hasattr(args, 'reference') and args.reference :
        if not os.path.isfile(args.reference) :
            print >> stderr, "ERROR: %s does not exist..." % args.reference
            exit(1)

    # verbosity
    if hasattr(args, 'verbose') :
        set_verbosity(args.verbose)

    # ensembl download method
    if hasattr(args, 'download_method') :
        set_ensembl_download_method(args.download_method)

    # contigs, labels, species, bam files
    # contigs, labels amd bam files must be unique
    #   apart from 'FAKE' bam file
    # length of contigs, label, species must be the same
    # length of bam files is either zero or same as contigs
    if hasattr(args, 'contigs') :
        if args.contigs :
            if not args.label or not args.species :
                print >> stderr, "ERROR: you must specify one --label and one --species argument per --contigs file!"
                exit(1)

            l = len(args.contigs)
            if (l != len(args.label)) or (l != len(args.species)) :
                print >> stderr, "ERROR: you must specify one --label and --species argument per --contigs file! (%d contig files, %d labels, %d species)" % \
                        (l, len(args.label), len(args.species))
                exit(1)

            if (l != len(set(args.label))) :
                print >> stderr, "ERROR: file labels must be unique!"
                exit(1)

            if (args.bam != None) and (len(args.bam) != l) :
                print >> stderr, "ERROR: if you must either specify bam files for all samples or no samples! (%d samples, %d bam files)" % \
                        (l, len(args.bam))
                exit(1)

            if args.bam == None :
                args.bam = [ None ] * l

            bad = False
            for fname in [ f for f in args.contigs + args.bam if f ] :
                if fname != 'FAKE' and not os.path.isfile(fname) :
                    print >> stderr, "ERROR: %s not found" % fname
                    bad = True

            if bad :
                exit(1)

    setup_logging()

def main() :
    args = handle_args(argv[1:])

    generic_options(args)

    start_time = time.time()

    get_log().info("this is glutton version %s" % str(glutton.__version__))

    ret = commands[argv[1]](args)

    get_log().info("%s took %s" % (argv[1], duration_str(time.time() - start_time)))

    return ret

if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> stderr, "Killed by user"
        exit(1)

