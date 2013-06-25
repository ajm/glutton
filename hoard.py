import sys

if (sys.version_info.major != 2) and (sys.version_info.minor < 7) :
    print >> sys.stderr, "Sorry: requires python version 2.7 or greater (but not python 3.x)\n"
    sys.exit(1)

try :
    import cogent

except ImportError :
    print >> sys.stderr, "Error: pycogent is not installed! (tested with version 1.5.3)"
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
            'build'     : lib.subcommands.build,
            'fix'       : lib.subcommands.fix,
            'assemble'  : lib.subcommands.assemble,
            'align'     : lib.subcommands.align,
            'debug'     : lib.subcommands.debug,
            'list'      : lib.subcommands.list_ensembl
           }

required_programs = {
            'build'     : ['prank','fastareformat','fasta2esd','esd2esi'],
            'fix'       : ['fastareformat','fasta2esd','esd2esi'],
            'assemble'  : ['sga'],
            'align'     : ['pagan','exonerate-server','exonerate'],
            'debug'     : [],
            'list'      : []
           }

def get_default_options() :
    return {
            'database'      : 'ensembl',
            'species'       : None,
            'release'       : None,
            'list'          : None,
            'workingdir'    : os.path.join(os.getcwd(), 'cache'),
            'tmpdir'        : os.environ.get('TMPDIR', '/tmp'),
            'restart'       : False,
            'verbose'       : False,
            'species2'      : None,
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
            'alignment-only': False,
            'force'         : False,
            'contig-minlen' : 200,
            'contig-file'   : None,
            'contig-outdir' : None
           }

def get_commands() :
    global commands
    tmp = commands.keys() 

    tmp.remove('debug')

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
    tmp = set()

    for i in required_programs :
        tmp.update(required_programs[i])

    return sorted(list(tmp))

def program_exists(programname) :
    for p in os.environ['PATH'].split(os.pathsep) :
        progpath = os.path.join(p, programname)
        if os.path.isfile(progpath) :
            # there may be another executable with the correct
            # permissions lower down in the path, but the shell
            # would not find it, so just return here...
            return os.access(progpath, os.X_OK)

    return False

def check_local_installation(required_programs) :
    bad = False

    # ensure required programs are installed locally
    for prog in required_programs :
        if not program_exists(prog) :
            print >> sys.stderr, "Error: '%s' is not installed." % prog
            bad = True

    if bad :
        sys.exit(1)

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
    -a      --alignment-only        (default = False, assumes --restart is not used)
            --restart               (default = False, start downloading from scratch)
            --specify-db='db info'  ('db info' is comma separated host, port, username, password)

%s options:
    -l      --list                  (list downloaded species and release versions)
    -s      --species='species'     (no default, use --list for options locally available, MANDATORY)
    -r      --release='release'     (default = latest available locally)
    -m      --min-length=NUM        (minimum length of contig to align, default = %d)
    -i      --input='file'          (input file containing contigs, MANDATORY)
    -o      --output='dir'          (output directory, default = location of contig file)

%s options:
            --prank='location'      (default = None, use system-wide version)
            --exonerate='location'  (default = None, use system-wide version)
            --sga='location'        (default = None, use system-wide version)
            --pagan='location'      (default = None, use system-wide version)

    -w      --workingdir='dir'      (default = %s)
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
        bold('align'),
        options['contig-minlen'],
        bold('misc'),
        options['workingdir'], 
        options['tmpdir'],
        options['threads'])

def expect_int(parameter, argument) :
    try :
        return int(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(1)

def parse_args(argv) :
    global databases

    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        argv,
                        "s:ld:r:a:w:t:c:vhfm:i:o:",
                        [   
                            "species=", 
                            "species2=",
                            "release=",
                            "workingdir=",
                            "tmpdir=",
                            "list",
                            "verbose", 
                            "help",
                            "restart",
                            "database=",
                            "specify-db=",
                            "prank=",
                            "pagan=",
                            "sga=",
                            "exonerate=",
                            "threads=",
                            "alignment-only",
                            "force",
                            "min-length=",
                            "input=",
                            "output="
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

        elif o in ('-w', '--workingdir') :
            options['workingdir'] = os.path.abspath(a)

        elif o in ('-t', '--tmpdir') :
            options['tmpdir'] = os.path.abspath(a)

        elif o in ('-r', '--release') :
            options['release'] = expect_int("release", a)

        elif o in ('--restart') :
            options['restart'] = True

        elif o in ('-a', '--alignment-only') :
            options['alignment-only'] = True

        elif o in ('-d', '--database') :
            options['database'] = a

        elif o in ('--specify-db') :
            parse_database_string(a, options)

        elif o in ('-m', '--min-length') :
            options['contig-minlen'] = expect_int("min-length", a)

        elif o in ('-i', '--input') :
            options['contig-file'] = os.path.abspath(a)

        elif o in ('-o', '--output-dir') :
            options['contig-outdir'] = os.path.abspath(a)

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

        else :
            assert False, "unhandled option %s" % o

    
    # check whether option combinations make sense
    if options['alignment-only'] and options['restart'] :
        print >> sys.stderr, "Error: 'alignment-only' and 'restart' cannot be used together, exiting..."
        sys.exit(1)

    if options['database'] not in databases.keys() :
        print >> sys.stderr, "Error: %s is not a predefined database, valid options are: %s" % \
                (pretty([options['database']]), pretty(databases.keys()))
        sys.exit(1)

    if options['database'] != 'user-defined' :
        fill_in_database_info(options)

    if options['contig-file'] and not options['contig-outdir']:
        options['contig-outdir'] = os.path.dirname(options['contig-file'])

    return options

def main() :
    global commands

    # a plaster over the command interface
    if (len(sys.argv) < 2) or (sys.argv[1] in ('-h', '--help', 'help')) :
        usage()
        return 0

    options = parse_args(sys.argv[2:])

    if sys.argv[1] in commands :
        check_local_installation(required_programs[sys.argv[1]])
        return commands[sys.argv[1]](options)

    print >> sys.stderr, "Error: command %s not recognised, valid command are %s" % \
            (pretty([sys.argv[1]]), pretty(get_commands()))
    
    return 1

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        pass

