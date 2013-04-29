import sys
import os
import getopt
import signal

from lib.cache import TranscriptCache, EnsemblInfo


cache = None

databases = {'ensembl' : { 
                'db-host' : 'ensembldb.ensembl.org', 
                'db-port' : 5306,
                'db-user' : 'anonymous',
                'db-pass' : '' },
             'ensembl-metazoa' : {
                'db-host' : 'mysql.ebi.ac.uk',
                'db-port' : 4157,
                'db-user' : 'anonymous',
                'db-pass' : '' }
             }

def get_default_options() :
    return {
            'database'   : 'ensembl',
            'species'    : None,
            'release'    : None,
            'list'       : None,
            'workingdir' : os.path.join(os.getcwd(), 'cache'),
            'tmpdir'     : os.environ.get('TMPDIR', '/tmp'),
            'resume'     : False,
            'verbose'    : False,
            'species2'   : None,
            'db-host'    : None,
            'db-port'    : None,
            'db-user'    : None,
            'db-pass'    : None,
            'prank'      : 'prank',
            'prank-threads' : 1,
            'alignment-only' : False
           }

def predefined_databases() :
    global databases
    return ', '.join(map(lambda x : "'%s'" % x, databases.keys()))

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

def clean_up() :
    global cache

    if cache :
        if not cache.stop :
            cache.stop = True
            print >> sys.stderr, "\n[!] Application will shutdown after the next gene family is downloaded or user types '^C' again...\n"
            return

    print >> sys.stderr, "\n[!] Forced exit from user..."
    sys.exit(0)

def handler_sigint(signal, frame) :
    clean_up()

def usage() :
    options = get_default_options()

    print >> sys.stderr, """Usage: %s [OPTIONS]
    -d      --database='db string'  (default = '%s', valid arguments : %s)
    -s      --species='species'     (no default, see output of --list for options)
    -o      --species2='species'    (default = None, only specify if you want orthologs from species2)
    -r      --release='release'     (default = latest)
    -w      --workingdir='dir'      (default = %s)
    -t      --tmpdir='dir'          (default = %s)
    -p      --prank='location'      (default = None, use system-wide version)
    -q      --prank-threads=NUM     (default = %d)
    -a      --alignment-only        (default = %s, assumes --continue)
    -c      --continue              (default = %s, don't start from scratch, but resume download)
    -l      --list
    -x      --specify-db='db info'  ('db info' is comma separated host, port, username, password)
    -v      --verbose
    -h      --help
""" % (sys.argv[0], options['database'], predefined_databases(),
        options['workingdir'], options['tmpdir'], options['prank-threads'], 
        options['alignment-only'], options['resume'])

def expect_int(parameter, argument) :
    try :
        return int(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(-1)

def parse_args() :
    global databases

    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        sys.argv[1:],
                        "s:o:r:hvlw:t:ad:x:p:q:c",
                        [   
                            "species=", 
                            "species2=",
                            "relsease=",
                            "workingdir=",
                            "tmpdir=",
                            "list",
                            "verbose", 
                            "help",
                            "continue",
                            "database=",
                            "specify-db=",
                            "prank=",
                            "prank-threads=",
                            "alignment-only"
                        ]
                    )

    except getopt.GetoptError, err :
        print >> sys.stderr, str(err) + "\n"
        usage()
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

        elif o in ('-c', '--continue') :
            options['resume'] = True

        elif o in ('-a', '--alignment-only') :
            options['alignment-only'] = True
            options['resume'] = True

        elif o in ('-d', '--database') :
            if a not in databases.keys() :
                print >> sys.stderr, "Error: '%s' is not a predefined database, valid strings are : %s" % (a, predefined_databases())
                sys.exit(-1)

            options['database'] = a

        elif o in ('-x', '--specify-db') :
            parse_database_string(a, options)

        elif o in ('-p', '--prank') :
            options['prank'] = a

        elif o in ('-q', '--prank-threads') :
            options['prank-threads'] = expect_int("prank-threads", a)

        else :
            assert False, "unhandled option %s" % o

    return options

def main() :
    global cache

    signal.signal(signal.SIGINT, handler_sigint)
    options = parse_args()

    if options['database'] != 'user-defined' :
        fill_in_database_info(options)

    if options['list'] :
        ei = EnsemblInfo(options)
        ei.print_species_table()
        return 0

    if options['species'] is None :
        print >> sys.stderr, "Error: you must specify a species!"
        return -1

    if options['release'] is None :
        ei = EnsemblInfo(options)
        options['release'] = ei.get_latest_release(options['species'])
        
        if options['release'] == -1 :
            print >> sys.stderr, "Could not find database for \'%s\', run --list command for species names" % options['species']
            return -1

        print >> sys.stderr, "Release not specified, using release %d" % options['release']
        

    cache = TranscriptCache(options)

    if options['resume'] and options['alignment-only'] :
        cache.download_complete = True
        cache.join()
        cache.build_exonerate_index()
        print "Info: done!"
    else :
        cache.build()

    return 0

if __name__ == '__main__' :
    sys.exit(main())

