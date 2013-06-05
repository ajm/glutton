import sys
import os
import subprocess

class Base(object) :
    def __init__(self, opt) :
        self.workingdir = opt['workingdir']
        self.tmpdir = opt['tmpdir']
        self.verbose = opt.get('verbose', False)

    def _check_dir(self, dirname, create=False) :
        if not os.path.exists(dirname) :
            if create :
                try :
                    os.makedirs(dirname)
                    self.info("created %s" % dirname)
                except OSError, ose :
                    self.error(str(ose))
                    sys.exit(-1)
            else :
                self.error("'%s' does not exist." % dirname)
                sys.exit(-1)

        elif not os.path.isdir(dirname) :
            self.error("'%s' exists, but is not a directory!" % dirname)
            sys.exit(-1)

        else :
            pass
        
        return dirname

    def _check_file(self, fname) :
        return os.path.exists(fname) and os.path.isfile(fname)

    def _swap_dirname(self, fname, directory) :
        return os.path.join(directory, os.path.basename(fname))

    def _create_file(self, fname, fcontents, dirname) :
        f = open(os.path.join(dirname, fname), 'w')
        f.write(fcontents)
        os.fsync(f)
        f.close()

    def _contents(self, fname) :
        return open(fname).read()

    def _p(self, s, print_anyway=False) :
        if self.verbose or print_anyway :
            print >> sys.stderr, s

    def info(self, s) :
        self._p("Info: %s" % s)

    def warn(self, s) :
        self._p("Warn: %s" % s, print_anyway=True)

    def error(self, s) :
        self._p("Error: %s" % s, print_anyway=True)

class ToolBase(Base) :
    def __init__(self, opt) :
        super(ToolBase, self).__init__(opt)

        self.binary_location = self.opt[self.name]

    @property
    def name(self) :
        return type(self).__name__.lower()

    def run(self, parameters, expected_outfiles) :
        DEVNULL = open('/dev/null', 'w')
        returncode = 0
        output = ""

        try :
            subprocess.check_output([self.binary_location] + parameters, stderr=DEVNULL, stdout=DEVNULL)

        except subprocess.CalledProcessError, cpe :
            returncode = cpe.returncode
            output = cpe.output

        finally :
            DEVNULL.close()

        if returncode == 0 :
            missing = [o for o in expected_outfiles if not self.check_file(o)]
            if len(missing) != 0 :
                returncode = -1
                output = "the following files were missing : %s" % ' '.join(missing)

        # if -f is not used then exit, otherwise just continue
        if not self.opt['force'] :
            self.error("%s returncode = %d\n\n%s\n" % (self.name, returncode, output))
            sys.exit(-1)

        return returncode, output

