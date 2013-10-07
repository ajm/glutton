import os

from lib.base import ToolBase

class Pagan(ToolBase) :
    def __init__(self, opt, alignment_file, tree_file, query_file, out_dir) :
        super(Pagan, self).__init__(opt)

        self.alignment_file = alignment_file
        self.tree_file = tree_file
        self.query_file = query_file
        self.out_file = os.path.join(out_dir, os.path.basename(query_file))

    def output_filenames(self) :
        return []

#        ext = ['.fas','.dna.fas']
#        
#        if self.tree_file :
#            ext.append('.tre')
#
#        return [self.out_file + e for e in ext]

    def run(self) :
        parameters = [
                      "--ref-seqfile",  self.alignment_file,
                      "--queryfile",    self.query_file,
                      "--outfile",      self.out_file,
                      "--fast-placement",
                      "--test-every-terminal-node",
                      "--translate", 
                      "--find-best-orf",
                      "--threads", "1"
                     ]

        if self.tree_file :
            parameters.append("--ref-treefile")
            parameters.append(self.tree_file)
        else :
            parameters.append("--pileup-alignment")

        returncode, output = self._execute(parameters, self.output_filenames())

        return returncode

