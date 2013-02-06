import sys
import subprocess
import tempfile

class Prank(object) :
    def __init__(self, tmpdir) :
        self.tmpdir = tmpdir

    def _tmpfilename(self) :
        return tempfile.mktemp(dir=self.tmpdir)

    def align(self, infile, outfile) :
        #command = "prank -d=%s -o=%s -translate -showtree" % (infile, outfile)
        command = ["prank", "-d=%s" % infile, "-o=%s" % outfile, "-translate", "-showtree"]

        #print "\n" + ' '.join(command) + "\n"

        try :
            subprocess.check_call(command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

        except subprocess.CalledProcessError, cpe :
            print >> sys.stderr, "Error: '%s' failed (%s)" % (command, str(cpe))
            sys.exit(-1)


        return map(lambda x: outfile + x, 
                [".1.dnd", ".2.dnd", ".nuc.1.fas", ".nuc.2.fas", ".pep.1.fas", ".pep.2.fas"])

