import tempfile
import os
import logging
import platform

from sys import exit
from multiprocessing import cpu_count


def tmpfile(contents=None, directory=None) :
    fd,name = tempfile.mkstemp(prefix='glutton', dir=directory)
    os.close(fd)

    if contents :
        with open(name, 'w') as f :
            f.write(contents)

    return name

def set_tmpdir(directory) :
    if not os.path.isdir(directory) :
        raise OSError('directory does not exist')

    tempfile.tempdir = directory

def tmpdir() :
    return tempfile.gettempdir()

def rm(f) :
    try :
        os.remove(f)
        return True

    except OSError :
        return False

def rm_f(f) :
    if isinstance(f, list) :
        for i in f :
            rm(i)
    else :
        rm(f)

def check_dir(d, create=False) :
    if os.path.isdir(d) :
        return

    if os.path.exists(d) :
        get_log().error("'%s' exists, but it is not a directory..." % d)
    else :
        if create :
            try :
                os.mkdir(d)
                get_log().info("created %s ..." % d)
                return

            except OSError, ose :
                get_log().error(str(ose))
        else :
            get_log().error("%s does not exist" % d)
    
    exit(1)

_glutton_threads = 0

def num_threads() :
    global _glutton_threads

    if _glutton_threads == 0 :
        return cpu_count()

    else :
        return _glutton_threads

def is_bad_threading_env() :
    macv = platform.mac_ver()

    # openmp based multithreading is broken in MacOS 10.7 
    if macv[0].startswith('10.7') :
        return True

    return False

def openmp_num_threads() :
    return 1 if is_bad_threading_env() else num_threads()

def set_threads(t) :
    global _glutton_threads
    _glutton_threads = t

def tmpfasta(seq) :
    fname = tmpfile()
    
    with open(fname, 'w') as f :
        if isinstance(seq, list) :
            for s in seq :
                print >> f, s.format('fasta').rstrip()
        else :
            print >> f, seq.format('fasta').rstrip()

    return fname

def get_binary_path(programname) :
    for p in os.environ['PATH'].split(os.pathsep) :
        progpath = os.path.join(p, programname)
        if os.path.isfile(progpath) :
            if os.access(progpath, os.X_OK) :
                return progpath
            break
    
    return None

_glutton_stream_loglevel = logging.ERROR

def set_verbosity(verbosity_level) :
    global _glutton_stream_loglevel

    if verbosity_level == 0 :
        _glutton_stream_loglevel = logging.ERROR
    elif verbosity_level == 1 :
        _glutton_stream_loglevel = logging.WARN
    elif verbosity_level == 2 :
        _glutton_stream_loglevel = logging.INFO
    else :
        _glutton_stream_loglevel = logging.DEBUG

def glutton_log_defaults(log) :
    global _glutton_stream_loglevel

    log.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(_glutton_stream_loglevel)
    ch.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    
    log.addHandler(ch) 


    fh = logging.FileHandler('glutton.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

    log.addHandler(fh)

    return log

def get_log() :
    return logging.getLogger('glutton')


if __name__ == '__main__' :
    print get_binary_path('ls')

    l = get_log()
    l = glutton_log_defaults(l)
    l.error("error")
    l.warn("warn")
    l.info("info")
    l.debug("debug")

