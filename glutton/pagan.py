from glutton.base import ExternalTool

class Pagan(ExternalTool) :
    def __init__(self) :
        super(Pagan, self).__init__()

    @property
    def version(self) :
        returncode, output = self._execute(["--version"], [], [])

        for line in output.split('\n') :
            if line.startswith('This is PAGAN') :
                v = line.strip().split()[-1]
                return v[:-1]

        raise ExternalToolError('could get version of pagan')

    @property
    def alignment(self) :
        return self.alignment_fname

    def output_filenames(self, outfile) :
        return [ outfile + i for i in ('.codon.fas', '.fas', '') ]

    def run(self, queries_fname, out_fname, alignment_fname, tree_fname=None) :
        parameters = [
                      "--ref-seqfile",  alignment_fname,
                      "--queryfile",    queries_fname,
                      "--outfile",      out_fname,
                      "--fast-placement",
                      "--test-every-terminal-node",
                      "--translate", 
                      "--find-best-orf",
                      "--threads", "1"
                     ]

        if tree_fname :
            parameters.append("--ref-treefile")
            parameters.append(tree_fname)
        else :
            parameters.append("--pileup-alignment")

        self.alignment_fname = out_fname + '.fas'

        returncode, output = self._execute(parameters, self.output_filenames(out_fname))

        #import os, sys
        #print >> sys.stderr, alignment_fname, tree_fname, queries_fname, self.alignment_fname
        #os._exit(0)

        return returncode

if __name__ == '__main__' :
    print Pagan().name, "version is", Pagan().version

