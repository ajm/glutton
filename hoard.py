import sys
import os
import getopt
import signal

from lib.cache import TranscriptCache, EnsemblInfo

cache = None

def get_default_options() :
    return {
            'verbose'    : False,
            'species'    : None,
            'release'    : None,
            'list'       : None,
            'workingdir' : "cache",
            'tmpdir'     : os.environ.get('TMPDIR', "/tmp"),
            'resume'     : False
           }

def clean_up() :
    global cache

    if cache :
        cache.shutdown()
    else :
        sys.exit(-1)

def handler_sigint(signal, frame) :
    clean_up()

def usage() :
    options = get_default_options()

    print >> sys.stderr, """Usage: %s [OPTIONS]
    -s      --species='species'     (no default)
    -r      --release='release'     (no default)
    -w      --workingdir='dir'      (default = %s)
    -t      --tmpdir='dir'          (default = %s)
    -a      --resume                (default = %s)
    -l      --list
    -v      --verbose
    -h      --help
""" % (sys.argv[0], options['workingdir'], options['tmpdir'], options['resume'])

def expect_int(parameter, argument) :
    try :
        return int(argument)

    except ValueError, ve :
        print >> sys.stderr, "Problem parsing argument for %s: %s\n" % (parameter, str(ve))
        usage()
        sys.exit(-1)

def parse_args() :
    options = get_default_options()

    try :
        opts,args = getopt.getopt(
                        sys.argv[1:],
                        "s:r:hvlw:t:a",
                        [   
                            "species=", 
                            "relsease=",
                            "workingdir=",
                            "tmpdir=",
                            "list",
                            "verbose", 
                            "help",
                            "resume"
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

        elif o in ('-w', '--workingdir') :
            options['workingdir'] = a

        elif o in ('-t', '--tmpdir') :
            options['tmpdir'] = a

        elif o in ('-r', '--release') :
            options['release'] = expect_int("release", a)

        elif o in ('-a', '--resume') :
            options['resume'] = True

        else :
            assert False, "unhandled option %s" % o

    return options

def main() :
    global cache

    signal.signal(signal.SIGINT, handler_sigint)
    options = parse_args()

    if options['list'] :
        ei = EnsemblInfo()
        ei.print_species_table()
        return 0

    cache = TranscriptCache(options['workingdir'], options['tmpdir'], options['species'], options['release'])
    cache.build(resume=options['resume'])

    return 0

if __name__ == '__main__' :
    sys.exit(main())

