from glutton.base import ExternalTool

class Pagan(ExternalTool) :
    def __init__(self) :
        super(Pagan, self).__init__()

    @property
    def version(self) :
        returncode, output = self._execute(["--version"], [])

        for line in output.split('\n') :
            if line.startswith('This is PAGAN') :
                v = line.strip().split()[-1]
                return v[:-1]

        raise ExternalToolError('could not get version of pagan')

    @property
    def nucleotide_alignment(self) :
        return self.nucleotide_alignment_fname

    @property
    def protein_alignment(self) :
        return self.protein_alignment_fname

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
                      "--threads", "1"
                      #"--find-best-orf",
                      #"--fragments"
                     ]

        if tree_fname :
            parameters.append("--ref-treefile")
            parameters.append(tree_fname)
        # unnecessary (see email 27/11/13 from ari)
        #else :
        #    parameters.append("--pileup-alignment")

        self.protein_alignment_fname = out_fname + '.fas'
        self.nucleotide_alignment_fname = out_fname + '.codon.fas'

        returncode, output = self._execute(parameters, self.output_filenames(out_fname))

        #import os, sys
        #print >> sys.stderr, alignment_fname, tree_fname, queries_fname, self.alignment_fname
        #os._exit(0)

        return returncode

if __name__ == '__main__' :
    print Pagan().name, "version is", Pagan().version

    from sys import argv, stderr, exit

    if len(argv) not in (4, 5) :
        print >> stderr, "Usage: %s qfile ofile afile [tfile]" % argv[0]
        exit(1)

    qfile = argv[1]
    ofile = argv[2]
    afile = argv[3]
    tfile = argv[4] if len(argv) == 5 else None

    p = Pagan()
    p.run(qfile, ofile, afile, tfile)
    
    print "output files are: %s" % ', '.join(p.output_filenames(ofile))
    
