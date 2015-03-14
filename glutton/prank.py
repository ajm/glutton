from glutton.base import ExternalTool, ExternalToolError

class Prank(ExternalTool) :
    def __init__(self) :
        super(Prank, self).__init__()

        self.tree_file = None
        self.alignment_file = None

    @property
    def version(self) :
        returncode, output = self._execute(["-version"], [])

        for line in output.split('\n') :
            if line.startswith('This is PRANK') :
                v = line.strip().split()[-1]
                return v[:-1]

        raise ExternalToolError('could not get version of prank')

    @property
    def tree(self) :
        return self.tree_file

    @property
    def alignment(self) :
        return self.alignment_file

    def output_filenames(self, outfile) :
        return [outfile + x for x in [".best.dnd", ".best.nuc.fas", ".best.pep.fas"]]

    def run(self, d, o) :
        parameters = [
                      "-d=%s" % d, 
                      "-o=%s" % o, 
                      "-showtree",
                      "-translate"
                     ]

        returncode, output = self._execute(parameters, self.output_filenames(o))

        self.tree_file = o + ".best.dnd"
        self.alignment_file = o + ".best.nuc.fas"

        return returncode

if __name__ == '__main__' :
    print Prank().name, "version is", Prank().version

