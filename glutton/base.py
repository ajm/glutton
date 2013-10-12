import sys
import os
import subprocess
import threading
import tempfile
import logging


class Base(object) :
    print_lock = threading.Lock()
    safe_write_lock = threading.Lock()
    log = None

    def __init__(self, opt) :
        self.opt = opt
        self.dbdir = opt['dbdir']
        self.tmpdir = opt['tmpdir']
        self.verbose = opt.get('verbose', False)
        self.force = opt.get('force', False)

        if not Base.log :
            Base.log = self._setup_logging()

        self.log = Base.log

    def _setup_logging(self) :
        log = logging.getLogger('glutton')
        log.setLevel(logging.DEBUG)

        fh = logging.FileHandler('glutton.log')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO if self.verbose else logging.WARNING)
        ch.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

        log.addHandler(fh)
        log.addHandler(ch)

        return log

    def _check_dir(self, dirname, create=False) :
        if not os.path.exists(dirname) :
            if create :
                try :
                    os.makedirs(dirname)
                    self.info("created %s" % dirname)
                except OSError, ose :
                    self.error(str(ose))
                    sys.exit(1)
            else :
                self.error("'%s' does not exist." % dirname)
                sys.exit(1)

        elif not os.path.isdir(dirname) :
            self.error("'%s' exists, but is not a directory!" % dirname)
            sys.exit(1)

        else :
            pass
        
        return dirname

    def _check_file(self, fname) :
        return os.path.exists(fname) and os.path.isfile(fname)

    def _swap_dirname(self, fname, directory) :
        return os.path.join(directory, os.path.basename(fname))

    def _change_ext(self, fname, new_ext) :
        return os.path.splitext(fname)[0] + '.' + new_ext

    def _create_file(self, fname, fcontents, dirname) :
        assert not os.path.isabs(fname)

        f = open(os.path.join(dirname, fname), 'w')
        f.write(fcontents)
        os.fsync(f)
        f.close()
        self.info("created %s" % fname)

        return f.name

    def _create_file2(self, fname, fcontents) :
        assert os.path.isabs(fname)

        f = open(fname, 'w')
        f.write(fcontents)
        f.close()
        self.info("created %s" % os.path.basename(fname))

    # XXX this is to be used exceptionally sparingly
    def _safe_write(self, fname, line) :
        assert os.path.isabs(fname)

        self.info("writing '%s' to %s ..." % (line, os.path.basename(fname)))

        self.safe_write_lock.acquire()

        try :
            f = open(fname, 'a')
            print >> f, line
            os.fsync(f)
            f.close()

        except IOError, ioe :
            self.warn(str(ioe))
        finally :
            self.safe_write_lock.release()

    def _contents(self, fname) :
        return open(fname).read()

    def _random_filename(self, prefix, dir) :
        fd,name = tempfile.mkstemp(prefix=prefix, dir=dir)
        os.close(fd)
        return name

    def _p(self, label, s, print_anyway=False, newline=True) :
        if self.verbose or print_anyway :
            Base.print_lock.acquire()
            #print >> sys.stderr, "%s: %s, %s" % (label, type(self).__name__, s)
            sys.stderr.write("%s: %s, %s%s" % (label, type(self).__name__, s, '\n' if newline else ''))
            #sys.stderr.flush()
            Base.print_lock.release()

    def _wrap(self, s) :
        return "%s %s" % (type(self).__name__, s)

    def info(self, s) :
        self.log.info(self._wrap(s))

    def warn(self, s) :
        self.log.warn(self._wrap(s))

    def error(self, s) :
        self.log.error(self._wrap(s))

    def debug(self, s) :
        self.log.debug(self._wrap(s))

    def progress(self, msg, percent) :
        self._p('Progress', "%s: %d%%\r" % (msg, int(percent)), newline=False, print_anyway=True)

    def overwrite(self, tag, msg, nl=False) :
        self._p('\r' + tag, msg, newline=nl, print_anyway=True)

    def rm(self, files) :
        if not isinstance(files, list) :
            files = [files]

        for f in files :
            try :
                os.remove(f)
                self.info("deleted %s" % f)
            except OSError, ose :
                pass

class ToolBase(Base) :
    def __init__(self, opt) :
        super(ToolBase, self).__init__(opt)

        self.binary_location = opt[self.name]

    @property
    def name(self) :
        return type(self).__name__.lower()

    def _execute(self, parameters, expected_outfiles) :
        returncode = 0
        output = ""

        self.info(' '.join(parameters))

        try :
            output = subprocess.check_output(
                                [self.binary_location] + parameters, 
                                stderr=subprocess.STDOUT,
                                close_fds=True
                                )

        except subprocess.CalledProcessError, cpe :
            returncode = cpe.returncode
            output = cpe.output

#        returncode,output = commands.getstatusoutput(' '.join([self.binary_location] + parameters))
#
#        if returncode != 0 :
#            if not self.opt['force'] :
#                self.error("return code = %d\n\n%s\n" % (returncode, output))
#                sys.exit(1)
#        else :
        if returncode == 0 :
            missing = [o for o in expected_outfiles if not self._check_file(o)]
            if len(missing) != 0 :
                returncode = -1
                output = "the following files were missing : %s" % ' '.join(missing)

        if returncode > 0 and not self.force :
            self.error(' '.join([self.binary_location] + parameters))
            self.error("return code = %d\n\n%s\n" % (returncode, output))
            #sys.exit(1)

        if returncode != 0 and self.force :
            self.warn("return code = %d" % (returncode))

        return returncode, output

