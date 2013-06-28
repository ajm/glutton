import sys

if sys.version_info < (2,7) :
    print >> sys.stderr, "Sorry: requires python version 2.7 or greater (but not python 3.x)\n"
    sys.exit(1)

try :
    import cogent
    del cogent

except ImportError :
    print >> sys.stderr, "Error: pycogent is not installed! (tested with version 1.5.3)\n"
    sys.exit(1)

import os
import getopt

import lib.subcommands


databases = {
             'ensembl' : { 
                'db-host' : 'ensembldb.ensembl.org', 
                'db-port' : 5306,
                'db-user' : 'anonymous',
                'db-pass' : '' },
             'ensembl-genomes' : {
                'db-host' : 'mysql.ebi.ac.uk',
                'db-port' : 4157,
                'db-user' : 'anonymous',
                'db-pass' : '' }
             }

commands = {
            'build'     : lib.subcommands.BuildCommand,
            'fix'       : lib.subcommands.FixCommand,
            'assemble'  : lib.subcommands.AssembleCommand,
            'align'     : lib.subcommands.AlignCommand,
            'list'      : lib.subcommands.ListCommand,
            'pack'      : lib.subcommands.PackCommand,
            'unpack'    : lib.subcommands.UnpackCommand,
            'rm'        : lib.subcommands.RmCommand   
        }

def get_default_options() :
    return {
            'database'      : 'ensembl',
            'species'       : None,
            'release'       : None,
            'list'          : None,
            'dbdir'    : os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cache'),
            'tmpdir'        : os.environ.get('TMPDIR', '/tmp'),
            'verbose'       : False,
            'db-host'       : None,
            'db-port'       : None,
            'db-user'       : None,
            'db-pass'       : None,
            'prank'         : 'prank',
            'pagan'         : 'pagan',
            'sga'           : 'sga',
            'exonerate'     : 'exonerate',
            'fastareformat' : 'fastareformat',
            'fasta2esd'     : 'fasta2esd',
            'esd2esi'       : 'esd2esi',
            'threads'       : 1,
            'force'         : False,
            'min-length'    : 200,
            'input-file'    : None,
            'output-dir'    : None,
            'use-exonerate-server' : False
           }

def get_commands() :
    global commands
    tmp = commands.keys() 

    return tmp

def fill_in_database_info(options) :
    global databases

    try :
        options.update(databases[options['database']])

    except KeyError, ke :
        print >> sys.stderr, "Error: '%s' is not a valid database string, valid strings are : %s" % \
                (options['database'], ' '.join(databases.keys()))
        sys.exit(1)

def parse_database_string(db_string, options) :
    db_fields = db_string.split(',')

    if len(db_fields) != 4 :
        print >> sys.stderr, "Error: '%s' is an invalid database string, valid strings are in the form hostname,port,username,password" % db_string
        sys.exit(1)

    # host, port, username, password
    try :
        db_fields[1] = int(db_fields[1])

    except ValueError, ve :
        print >> sys.stderr, "Error: %s is not a valid port number" % db_fields[1]
        sys.exit(1)

    options['db-host'] = db_fields[0]
    options['db-port'] = db_fields[1]
    options['db-user'] = db_fields[2]
    options['db-pass'] = db_fields[3]
    options['database'] = 'user-defined'

def bold(s) :
    return "\033[1m%s\033[0m" % s

def bold_all(l) :
    return map(bold, l)

def quote(s) :
    return "'%s'" % s

def quote_all(l) :
    return map(quote, l)

def list_sentence(l) :
    if len(l) < 2 :
        return "".join(l)
    return "%s and %s" % (', '.join(l[:-1]), l[-1])

def pretty(l) :
    return list_sentence(bold_all(quote_all(map(str, sorted(l)))))

def get_required_programs() :
    return ['esd2esi','exonerate','exonerate-server','fasta2esd','fastareformat','pagan','prank','sga']

def usage() :
    global databases

    options = get_default_options()

    print >> sys.stderr, """Usage: %s command [OPTIONS]

Legal commands are %s (see below for options).
%s assumes that the following programs are installed: %s.

%s options:
    -l      --list                  (list species and release versions from current database)
    -s      --species='species'     (no default, use --list for options with current database, MANDATORY)
    -r      --release='release'     (default = latest)
    -d      --database='db string'  (default = %s, valid arguments : %s)
            --specify-db='db info'  ('db info' is comma separated host, port, username, password)

%s options:
    -l      --list
    -s      --species='species'
    -r      --release='release'

%s options:
    -i      --input-file='file'

%s options:
    -l      --list                  (list downloaded species and release versions)
    -s      --species='species'     (no default, use --list for options locally available, MANDATORY)
    -r      --release='release'     (default = latest available locally)
    -m      --min-length=NUM        (minimum length of contig to align, default = %d)
    -i      --input-file='file'     (input file containing contigs, MANDATORY)
    -o      --output-dir='dir'      (output directory, default = location of contig file)
            --use-exonerate-server

%s options:
            --prank='location'      (default = None, use system-wide version)
            --exonerate='location'  (default = None, use system-wide version)
            --sga='location'        (default = None, use system-wide version)
            --pagan='location'      (default = None, use system-wide version)

    -b      --dbdir='dir'           (default = %s)
    -t      --tmpdir='dir'          (default = %s)
    -c      --threads=NUM           (default = %s)

    -f      --force    
    -v      --verbose
    -h      --help
""" % (
        sys.argv[0], 
        pretty(get_commands()),
        sys.argv[0],
        pretty(get_required_programs()),
        bold('build'),
        pretty([options['database']]),
        pretty(databases.keys()),
        bold('pack'),
        bold('unpack'),
        bold('align'),
        options['min-length'],
        bold('misc'),
        options['dbdir'], 
        options['tmpdir'],
        options['threads'])

def expect_int(parameter, argument) :
    try :
        return int(argument)

    except ValueError, ve :
        print >> sys.stderr, "Error: parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(1)

def expect_file(parameter, argument) :
    if not os.path.isfile(argument) :
        print >> sys.stderr, "Error: file '%s' does not exist" % argument
        sys.exit(1)

    return os.path.abspath(argument)

def expect_dir(parameter, argument) :
    if not os.path.isdir(argument) :
        print >> sys.stderr, "Error: directory '%s' does not exist" % argument
        sys.exit(1)

    return os.path.abspath(argument)

def parse_args(argv) :
    global databases

    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        argv,
                        "s:ld:r:b:t:c:vhfm:i:o:",
                        [   
                            "species=", 
                            "species2=",
                            "release=",
                            "dbdir=",
                            "tmpdir=",
                            "list",
                            "verbose", 
                            "help",
                            "database=",
                            "specify-db=",
                            "prank=",
                            "pagan=",
                            "sga=",
                            "exonerate=",
                            "threads=",
                            "force",
                            "min-length=",
                            "input-file=",
                            "output-dir=",
                            "use-exonerate-server"
                        ]
                    )

    except getopt.GetoptError, err :
        print >> sys.stderr, str(err) + " (see --help for more details)"
        sys.exit(1)

    for o,a in opts :
        if o in ('-h', '--help') :
            usage()
            sys.exit(0)

        elif o in ('-v', '--verbose') :
            options['verbose'] = True

        elif o in ('-l', '--list') :
            options['list'] = True

        elif o in ('-s', '--species') :
            options['species'] = a

        elif o in ('-b', '--dbdir') :
            options['dbdir'] = os.path.abspath(a)

        elif o in ('-t', '--tmpdir') :
            options['tmpdir'] = os.path.abspath(a)

        elif o in ('-r', '--release') :
            options['release'] = expect_int('release', a)

        elif o in ('-d', '--database') :
            options['database'] = a

        elif o in ('--specify-db') :
            parse_database_string(a, options)

        elif o in ('-m', '--min-length') :
            options['min-length'] = expect_int('min-length', a)

        elif o in ('-i', '--input-file') :
            options['input-file'] = expect_file('input-file', a)

        elif o in ('-o', '--output-dir') :
            options['output-dir'] = expect_dir('output-dir', a)

        elif o in ('--prank', '--pagan', '--exonerate', '--sga') :
            program = o[2:]
            options[program] = os.path.join(a, program)

            if program == 'exonerate' :
                for p in ['fastareformat', 'fasta2esd', 'esd2esi'] :
                    options[p] = os.path.join(a, p)

        elif o in ('-c', '--threads') :
            options['threads'] = expect_int("threads", a)

        elif o in ('-f', '--force') :
            options['force'] = True

        elif o in ('--use-exonerate-server') :
            options['use-exonerate-server'] = True

        else :
            assert False, "unhandled option %s" % o

    
    if options['database'] not in databases :
        print >> sys.stderr, "Error: %s is not a predefined database, valid options are: %s" % \
                (quote(options['database']), list_sentence(quote_all(databases.keys())))
        sys.exit(1)

    if options['database'] != 'user-defined' :
        fill_in_database_info(options)

    if options['input-file'] and not options['output-dir']:
        options['output-dir'] = os.path.dirname(options['input-file'])

    return options

def main() :
    global commands

    # a plaster over the command interface
    if (len(sys.argv) < 2) or (sys.argv[1] in ('-h', '--help', 'help')) :
        usage()
        return 0

    options = parse_args(sys.argv[2:])

    if sys.argv[1] in commands :
        return commands[sys.argv[1]](options).run()

    print >> sys.stderr, "Error: command %s not recognised, valid command are %s" % \
            (pretty([sys.argv[1]]), pretty(get_commands()))
    
    return 1

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        pass

