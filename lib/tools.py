import sys
import subprocess
import tempfile

class PrankError(Exception) :
    pass

class Prank(object) :
    def __init__(self, tmpdir, prank_binary) :
        self.tmpdir = tmpdir
        self.binary = prank_binary

    def _tmpfilename(self) :
        return tempfile.mktemp(dir=self.tmpdir)
    
    def align(self, infile, outfile) :
        #command = "prank -d=%s -o=%s -translate -showtree" % (infile, outfile)
        command = [self.binary, "-d=%s" % infile, "-o=%s" % outfile, "-translate", "-showtree"]

        DEVNULL = open('/dev/null', 'w')

        try :
            subprocess.check_call(command, stderr=DEVNULL, stdout=DEVNULL)

        except subprocess.CalledProcessError, cpe :
            #print >> sys.stderr, "Error: '%s' failed (%s)" % (command, str(cpe))
            DEVNULL.close()
            raise PrankError()
        
        DEVNULL.close()

        return map(lambda x: outfile + x, 
                [".1.dnd", ".2.dnd", ".nuc.1.fas", ".nuc.2.fas", ".pep.1.fas", ".pep.2.fas"])

