from sys import stderr, argv, exit, version_info

if version_info < (2,7) :
    print >> stderr, "Sorry: requires python version 2.7 or greater (but not python 3.x)\n"
    exit(1)

import os
import getopt
import tempfile

import glutton.subcommands


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
            'build'     : glutton.subcommands.BuildCommand,
            'fix'       : glutton.subcommands.FixCommand,
            'align'     : glutton.subcommands.AlignCommand,
            'list'      : glutton.subcommands.ListCommand,
            'pack'      : glutton.subcommands.PackCommand,
            'unpack'    : glutton.subcommands.UnpackCommand,
            'rm'        : glutton.subcommands.RmCommand,
            'scaffold'  : glutton.subcommands.ScaffoldCommand,
            'check'     : glutton.subcommands.CheckCommand
        }

def get_default_options() :
    return {
            'database'      : 'ensembl',
            'species'       : None,
            'release'       : None,
            'list'          : None,
            'dbdir'         : os.path.join(os.environ.get('HOME'), '.local/glutton'),
            'tmpdir'        : os.environ.get('TMPDIR', '/tmp'),
            'verbose'       : False,
            'db-host'       : None,
            'db-port'       : None,
            'db-user'       : None,
            'db-pass'       : None,
            'prank'         : 'prank',
            'pagan'         : 'pagan',
            'exonerate'     : 'exonerate',
            'fastareformat' : 'fastareformat',
            'fasta2esd'     : 'fasta2esd',
            'esd2esi'       : 'esd2esi',
            'threads'       : 1,
            'force'         : False,
            'min-length'    : 200,
            'min-identity'  : 0.95,
            'scaffold-file' : 'scaffolds.fasta',
            'contig-file'   : None,
            'alignment-dir' : os.path.join(os.getcwd(), 'alignments'),
            'package-file'  : None,
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
        print >> stderr, "Error: '%s' is not a valid database string, valid strings are : %s" % \
                (options['database'], ' '.join(databases.keys()))
        exit(1)

def parse_database_string(db_string, options) :
    db_fields = db_string.split(',')

    if len(db_fields) != 4 :
        print >> stderr, "Error: '%s' is an invalid database string, valid strings are in the form hostname,port,username,password" % db_string
        exit(1)

    # host, port, username, password
    try :
        db_fields[1] = int(db_fields[1])

    except ValueError, ve :
        print >> stderr, "Error: %s is not a valid port number" % db_fields[1]
        exit(1)

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
    return ['esd2esi','exonerate','exonerate-server','fasta2esd','fastareformat','pagan','prank']

def usage() :
    global databases

    options = get_default_options()

    print >> stderr, """Usage: %s command [OPTIONS]

Legal commands are %s (see below for options).
%s assumes that the following programs are installed: %s.

%s options:
    -l      --list                  (list species and release versions from current database)
    -s      --species='species'     (no default, use --list for options with current database, MANDATORY)
    -r      --release='release'     (default = latest)
    -d      --database='db string'  (default = %s, valid arguments : %s)
            --specifydb='db info'   ('db info' is comma separated host, port, username, password)

%s options:
    -l      --list                  (list downloaded species and release versions)
    -s      --species='species'     (no default, use --list for options locally available, MANDATORY)
    -r      --release='release'     (default = latest available locally)
    -m      --minlength=NUM         (minimum length of contig to align, default = %d)
    -i      --contigs='file'        (input file containing contigs, MANDATORY)
    -a      --alignments='dir'      (output directory, default = %s)

%s options:
    -i      --contigs='file'        (input file containing contigs, MANDATORY)
    -a      --alignments='dir'      (output directory of 'align' command)
            --minidentity=FLOAT     (threshold for using alignments, default = %.2f)
            --scaffolds='file'      (output of scaffolder, default = %s)

%s options:
    -l      --list                  (list downloaded species and release versions)
    -s      --species='species'     (no default, use --list for options locally available, MANDATORY)
    -r      --release='release'     (no default, use --list for options locally available, MANDATORY)

%s options:
    -p      --package='file'        (no default, 'file' must come from the pack command)

%s options:
    -l      --list                  (list downloaded species and release versions)
    -s      --species='species'     (no default, use --list for options locally available, MANDATORY)
    -r      --release='release'     (no default, use --list for options locally available, MANDATORY)

%s options:
            --prank='location'      (default = None, use system-wide version)
            --exonerate='location'  (default = None, use system-wide version)
            --pagan='location'      (default = None, use system-wide version)

    -b      --dbdir='dir'           (default = %s)
    -t      --tmpdir='dir'          (default = %s)
    -c      --threads=NUM           (default = %s)

    -f      --force    
    -v      --verbose
    -h      --help
""" % (
        argv[0], 
        pretty(get_commands()),
        argv[0],
        pretty(get_required_programs()),
        bold('build'),
        pretty([options['database']]),
        pretty(databases.keys()),
        bold('align'),
        options['min-length'],
        options['alignment-dir'],
        bold('scaffold'),
        options['min-identity'],
        options['scaffold-file'],
        bold('pack'),
        bold('unpack'),
        bold('rm'),
        bold('misc'),
        options['dbdir'], 
        options['tmpdir'],
        options['threads'])

def expect_type(parameter, argument, this_type) :
    try :
        return this_type(argument)

    except ValueError, ve :
        print >> stderr, "Error: parsing argument for %s: %s\n" % (parameter, str(ve))
        exit(1)

def expect_int(parameter, argument) :
    return expect_type(parameter, argument, int)

def expect_float(parameter, argument) :
    return expect_type(parameter, argument, float)

def expect_file(parameter, argument) :
    if not os.path.isfile(argument) :
        print >> stderr, "Error: file '%s' does not exist" % argument
        exit(1)

    return os.path.abspath(argument)

def expect_dir(parameter, argument, should_exist=True) :
    if should_exist and not os.path.isdir(argument) :
        print >> stderr, "Error: directory '%s' does not exist" % argument
        exit(1)

    return os.path.abspath(argument)

def parse_args(argv) :
    global databases

    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        argv,
                        "ls:r:d:m:i:a:p:b:t:c:fvh",
                        [   
                            "list",
                            "species=", 
                            "release=",
                            "database=",
                            "specifydb=",
                            "minlength=",
                            "contigs=",
                            "alignments=",
                            "minidentity=",
                            "scaffolds=",
                            "package=",
                            "prank=",
                            "exonerate=",
                            "pagan=",
                            "dbdir=",
                            "tmpdir=",
                            "threads=",
                            "force",
                            "verbose", 
                            "help"
                        ]
                    )

    except getopt.GetoptError, err :
        print >> stderr, str(err) + " (see --help for more details)"
        exit(1)

    for o,a in opts :
        if o in ('-h', '--help') :
            usage()
            exit(0)

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

        elif o in ('--specifydb') :
            parse_database_string(a, options)

        elif o in ('-i', '--contigs') :
            options['contig-file'] = expect_file('contigs', a)

        elif o in ('-a', '--alignments') :
            options['alignment-dir'] = expect_dir('alignments', a, should_exist=False)

        elif o in ('--scaffolds') :
            options['scaffold-file'] = a

        elif o in ('-p', '--package') :
            options['package-file'] = expect_file('package', a)

        elif o in ('-m', '--minlength') :
            options['min-length'] = expect_int('min-length', a)

        elif o in ('--prank', '--pagan', '--exonerate') :
            program = o[2:]
            options[program] = os.path.join(a, program)

            if program == 'exonerate' :
                for p in ['fastareformat', 'fasta2esd', 'esd2esi'] :
                    options[p] = os.path.join(a, p)

        elif o in ('-c', '--threads') :
            options['threads'] = expect_int("threads", a)

        elif o in ('-f', '--force') :
            options['force'] = True

        elif o in ('--minidentity') :
            tmp = expect_float('minidentity', a)
            if tmp < 0.0 or tmp > 1.0 :
                print >> stderr, "Error: minidentity must be in range 0.0 - 1.0\n"
                exit(1)
            
            options['min-identity'] = tmp

        else :
            assert False, "unhandled option %s" % o

    
    if options['database'] not in databases :
        print >> stderr, "Error: %s is not a predefined database, valid options are: %s" % \
                (quote(options['database']), list_sentence(quote_all(databases.keys())))
        exit(1)

    if options['database'] != 'user-defined' :
        fill_in_database_info(options)

    # change tmpdir just in case other commands are being run concurrently
    if options['tmpdir'] == os.environ.get('TMPDIR', '/tmp') :
        options['tmpdir'] = tempfile.mkdtemp(prefix='glutton_', dir=options['tmpdir'])

    return options

def main() :
    global commands

    # a plaster over the command interface
    if (len(argv) < 2) or (argv[1] in ('-h', '--help', 'help')) :
        usage()
        return 0

    options = parse_args(argv[2:])

    if argv[1] in commands :
        return commands[argv[1]](options).run()

    print >> stderr, "Error: command %s not recognised, valid command are %s" % \
            (pretty([argv[1]]), pretty(get_commands()))
    
    return 1

if __name__ == '__main__' :
    try :
        exit(main())

    except KeyboardInterrupt :
        print >> stderr, "Killed by user"
        exit(1)
