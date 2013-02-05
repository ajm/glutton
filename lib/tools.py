import sys
import os
import tempfile

class Prank(object) :
    def __init__(self, tmpdir) :
        self.tmpdir = tmpdir

    def _tmpfilename(self) :
        return tempfile.mktemp(dir=self.tmpdir)

    def align(self, infile, outfile) :
        command = "prank -d=%s -o=%s -translate -showtree &> /dev/null" % (infile, outfile)

        if os.system(command) != 0 :
            print >> sys.stderr, "Error: '%s' failed" % command
            sys.exit(0)

        return map(lambda x: outfile + x, 
                [".1.dnd", ".2.dnd", ".nuc.1.fas", ".nuc.2.fas", ".pep.1.fas", ".pep.2.fas"])

