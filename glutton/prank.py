from glutton.base import ToolBase

class Prank(ToolBase) :
    def __init__(self, opt, infile) :
        super(Prank, self).__init__(opt)
        
        self.infile = infile
        self.outfile = self._swap_dirname(infile, self.tmpdir)

    def output_filenames(self) :
        return [self.outfile + x for x in [".1.dnd", ".2.dnd", ".nuc.1.fas", ".nuc.2.fas", ".pep.1.fas", ".pep.2.fas"]]

    def run(self) :
        parameters = [
                      "-d=%s" % self.infile, 
                      "-o=%s" % self.outfile, 
                      "-translate", 
                      "-showtree"
                     ]

        returncode, output = self._execute(parameters, self.output_filenames())

        return returncode

