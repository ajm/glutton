import subprocess
import os

from glutton.utils import get_log, tmpfile, openmp_num_threads, rm_f


class All_vs_all_search(object) :
    def __init__(self) :
        self.log = get_log()
        self.cleanup_files = []

    def _run(self, c) :
        try :
            self.log.debug(" ".join(c))
            output = subprocess.check_output(c, stderr=subprocess.STDOUT, close_fds=True)

        except subprocess.CalledProcessError, cpe :
            self.log.error("%s returncode=%d\n%s" % (c[0], cpe.returncode, cpe.output))
            os._exit(1)

    def process(self, db, queries, min_identity, min_alignment_length) :
        min_pident = 100 * min_identity

        # creates db + {phr,pin,psq} in same dir as db
        self.log.info("creating blast db...")
        self._run(["makeblastdb", "-in", db, "-dbtype", "prot"])

        outfile = tmpfile()

        print outfile

        self.cleanup_files += [db + i for i in [".phr",".pin",".psq"]]
        #self.cleanup_files.append(outfile)

        # setting threads > 1 makes it crash :-(
        self.log.info("running all vs. all search...")
        self._run(["blastx",
            "-query", queries,
            "-db", db,
            "-out", outfile,
            "-max_target_seqs", "1",
            "-outfmt", "6",
            "-num_threads", str(openmp_num_threads())
            ])

        # read results and output mapping
        # field 3 = pident (as a percentage)
        # field 4 = length
        tmp = {}

        with open(outfile) as f :
            for line in f :
                qseq,gene,ident,leng = line.strip().split()[:4]
                
                print "id=%f alen=%d len(%s) id(%s)" % (float(ident), int(leng), str(min_alignment_length > int(leng)), str(min_pident > float(ident)))

                if min_pident > float(ident) or min_alignment_length > int(leng) :
                    break

                tmp[qseq] = gene #.append((qseq, gene))

        rm_f(self.cleanup_files)
        
        return tmp

    def kill(self) :
        rm_f(self.cleanup_files)

if __name__ == '__main__' :

    from glutton.utils import get_log, glutton_log_defaults

    glutton_log_defaults(get_log())

    ava = All_vs_all_search()
    tmp = ava.process('tc_test.fasta', 'queries_test.fasta', 0.8, 100)
    
    for m in tmp :
        print m, tmp[m]

