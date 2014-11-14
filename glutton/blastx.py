from glutton.base import ExternalTool
from glutton.utils import get_log

from sys import exit
import subprocess
import os


class Blastx(ExternalTool) :
    def __init__(self) :
        super(Blastx, self).__init__()

        self._results = []

    @property
    def version(self) :
        returncode, output = self._execute(["-version"], [])

        for line in output.split('\n') :
            if line.startswith('blastx') :
                v = line.strip().split()[-1]
                return v[:-1]

        raise ExternalToolError('could not get version of blastx')
    
    @staticmethod
    def makedb(db_fname, nucleotide) :
        c = ["makeblastdb", "-in", db_fname, "-dbtype", "nucl" if nucleotide else "prot"]

        try :
            get_log().debug(" ".join(c))
            subprocess.check_output(c, stderr=subprocess.STDOUT, close_fds=True)

        except subprocess.CalledProcessError, cpe :
            self.log.error("%s returncode=%d\n%s" % (c[0], cpe.returncode, cpe.output))
            exit(1)

    @property
    def results(self) :
        return self._results

    def run(self, query, database, outfile) :
        parameters = [
            "-query", query,
            "-db", database,
            "-out", outfile,
            "-max_target_seqs", "1",
            "-outfmt", "6"
            ]

        returncode, output = self._execute(parameters, [])

        with open(outfile) as f :
            for line in f :
                line = line.strip()

                if line == "" :
                    continue

                try :
                    contig,gene,identity,length = line.split()[:4]
                    identity = float(identity)
                    length = int(length)
                    self._results.append((contig,gene,identity,length))            

                except ValueError :
                    self.log.warn("bad line returned by blastx (%s)" % line)
                    continue


        return returncode

if __name__ == '__main__' :
     print Blastx().name, "version is", Blastx().version

