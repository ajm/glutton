import sys
import os
import getopt

from lib.subcommands import build, fix, assemble, align, debug, list_ensembl
from lib.system import System
from lib.cache import EnsemblInfo

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
            'build'     : build,
            'fix'       : fix,
            'assemble'  : assemble,
            'align'     : align,
            'debug'     : debug,
            'list'      : list_ensembl
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
            'database'   : 'ensembl',
            'species'    : None,
            'release'    : None,
            'list'       : None,
            'workingdir' : os.path.join(os.getcwd(), 'cache'),
            'tmpdir'     : os.environ.get('TMPDIR', '/tmp'),
            'restart'    : False,
            'verbose'    : False,
            'species2'   : None,
            'db-host'    : None,
            'db-port'    : None,
            'db-user'    : None,
            'db-pass'    : None,
            'prank'      : 'prank',
            'pagan'      : 'pagan',
            'sga'        : 'sga',
            'exonerate'  : 'exonerate',
            'threads'    : 1,
            'alignment-only' : False,
            'force'      : False
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
        sys.exit(-1)

def parse_database_string(db_string, options) :
    db_fields = db_string.split(',')

    if len(db_fields) != 4 :
        print >> sys.stderr, "Error: '%s' is an invalid database string, valid strings are in the form hostname,port,username,password" % db_string
        sys.exit(-1)

    # host, port, username, password
    try :
        db_fields[1] = int(db_fields[1])

    except ValueError, ve :
        print >> sys.stderr, "Error: %s is not a valid port number" % db_fields[1]
        sys.exit(-1)

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
    return list_sentence(bold_all(quote_all(sorted(l))))

def get_required_programs() :
    tmp = set()

    for i in required_programs :
        tmp.update(required_programs[i])

    return sorted(list(tmp))

def usage() :
    global databases

    options = get_default_options()

    print >> sys.stderr, """Usage: %s command [OPTIONS]

Legal commands are %s (see below for options).
%s assumes that the following programs are installed: %s.

Mandatory:
    -s      --species='species'     (no default, use --list for options)

%s options:
    -l      --list                  (list species and release versions from current database)
    -d      --database='db string'  (default = '%s', valid arguments : %s)
    -r      --release='release'     (default = latest)
    -a      --alignment-only        (default = False, assumes --restart is not used)
            --restart               (default = False, start downloading from scratch)
            --specify-db='db info'  ('db info' is comma separated host, port, username, password)

%s options:
    -l      --list                  (list downloaded species and release versions)

%s options:
            --prank='location'      (default = None, use system-wide version)
            --exonerate='location'  (default = None, use system-wide version)
            --sga='location'        (default = None, use system-wide version)
            --pagan='location'      (default = None, use system-wide version)

    -w      --workingdir='dir'      (default = %s)
    -t      --tmpdir='dir'          (default = %s)
    -c      --threads=NUM           (default = %s)
    
    -v      --verbose
    -h      --help
""" % (
        sys.argv[0], 
        pretty(get_commands()),
        sys.argv[0],
        pretty(get_required_programs()),
        bold('build'),
        options['database'], 
        pretty(databases.keys()),
        bold('align'),
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
        sys.exit(-1)

def parse_args(argv) :
    global databases

    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        argv,
                        "s:ld:r:a:w:t:c:vhf",
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
                            "force"
                        ]
                    )

    except getopt.GetoptError, err :
        print >> sys.stderr, str(err) + " (see --help for more details)"
        sys.exit(-1)

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

        elif o in ('-o', '--species2') :
            options['species2'] = a

        elif o in ('-w', '--workingdir') :
            options['workingdir'] = a

        elif o in ('-t', '--tmpdir') :
            options['tmpdir'] = a

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

        elif o in ('--prank', '--pagan', '--exonerate', '--sga') :
            options[o[2:]] = a

        elif o in ('-c', '--threads') :
            options['threads'] = expect_int("threads", a)

        elif o in ('-f', '--force') :
            options['force'] = True

        else :
            assert False, "unhandled option %s" % o

    
    # check whether option combinations make sense
    if options['alignment-only'] and options['restart'] :
        print >> sys.stderr, "Error: 'alignment-only' and 'restart' cannot be used together, exiting..."
        sys.exit(-1)

    if options['database'] not in databases.keys() :
        print >> sys.stderr, "Error: %s is not a predefined database, valid options are: %s" % (pretty([options['database']]), pretty(databases.keys()))
        sys.exit(-1)

    if options['database'] != 'user-defined' :
        fill_in_database_info(options)

    return options

def main() :
    global commands

    # a plaster over the command interface
    if sys.argv[1] in ('-h', '--help') :
        usage()
        sys.exit(0)

    options = parse_args(sys.argv[2:])

    if sys.argv[1] in commands :
        System().check_local_installation(required_programs[sys.argv[1]])
        return commands[sys.argv[1]](options)

    print >> sys.stderr, "Error: command %s not recognised, valid command are %s" % \
            (pretty([sys.argv[1]]), pretty(get_commands()))
    return -1

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        pass

