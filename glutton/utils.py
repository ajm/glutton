import tempfile
import os
import logging
import platform
import hashlib
import threading

from sys import exit
from multiprocessing import cpu_count


def tmpfile(contents=None, directory=None, suffix='') :
    fd,name = tempfile.mkstemp(prefix='glutton', dir=directory, suffix=suffix)
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

def _write_six_orfs(fout, seq) :
    for i in range(3) :
        print >> fout, ">%s\n%s" % (seq.id + '_orf' + str(i+1), seq[i:])

    seq.reverse_complement()

    for i in range(3) :
        print >> fout, ">%s\n%s" % (seq.id + '_orf' + str(i+4), seq[i:])

    seq.reverse_complement()

def tmpfasta_orfs(seq) :
    fname = tmpfile()

    with open(fname, 'w') as f :
        if isinstance(seq, list) :
            for s in seq :
                _write_six_orfs(f, s)
        else :
            _write_six_orfs(f, seq)

    return fname

def get_binary_path(programname) :
    for p in os.environ['PATH'].split(os.pathsep) :
        progpath = os.path.join(p, programname)
        if os.path.isfile(progpath) :
            if os.access(progpath, os.X_OK) :
                return progpath
            break
    
    return None

_lock = threading.Lock()

def threadsafe_io(fname, s) :
    global _lock

    _lock.acquire()

    with open(fname, 'a') as f :
        print >> f, s

    _lock.release()

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

def duration_str(seconds) :
    seconds = int(seconds)

    hours = seconds / 3600
    seconds %= 3600
    minutes = seconds / 60
    seconds %= 60

    s = ""
    if hours :
        s += ("%d hour%s " % (hours, "" if hours == 1 else "s"))

    if minutes :
        s += ("%d minute%s " % (minutes, "" if minutes == 1 else "s"))

    s += ("%d second%s" % (seconds, "" if seconds == 1 else "s"))

    return s

def md5(fname) :
    m = hashlib.md5()
    f = open(fname)

    for block in f.read(1024) :
        m.update(block)

    f.close()
    return m.hexdigest()

def log_defaults(log) :
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

def setup_logging() :
    log_defaults(get_log())

def get_log() :
    return logging.getLogger('glutton')


if __name__ == '__main__' :
    print get_binary_path('ls')

    set_verbosity(3)
    setup_logging()

    l = get_log()
    l.error("error")
    l.warn("warn")
    l.info("info")
    l.debug("debug")

