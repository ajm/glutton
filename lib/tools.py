import sys
import subprocess
import tempfile
import os

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

class ExonerateHit(object) :
    def __init__(self, query, match) :
        self.query = query
        self.match = match

    def __str__(self) :
        return ' '.join([str(self.query), str(self.match)])

class ExonerateVulgar(object) :
    def __init__(self, seqid, start, stop) :
        self.seqid = seqid
        self.start = start
        self.stop = stop

    def __str__(self) :
        return "%s %d %d" % (self.seqid, self.start, self.stop)

class Exonerate(object) :
    def __init__(self) :
        self.exonerate = None

    def start(self, d) :
        index_file = os.path.join(d, 'exonerate.esi')
        args = ['exonerate-server', 'exonerate.esi',
                '--port', '12887']

        if not (os.path.exists(index_file) and os.path.isfile(index_file)) :
            print >> sys.stderr, "Error: %d does not exist..." % index_file
            sys.exit(-1)

        self.exonerate = subprocess.Popen(args,
                        cwd=d,
                        stderr=subprocess.STDOUT,
                        stdout=subprocess.PIPE)

        # wait until it is listening...
        out = []
        for line in iter(self.exonerate.stdout.readline, '') :
            out.append(line)

            if 'listening on port' in line :
                print >> sys.stderr, "Info: exonerate-server started..."
                break
        else :
            print >> sys.stderr, "Error: exonerate-server did not start properly:\n\n\t%s\n" % '\t'.join(out)
            sys.exit(-1)

        # change the redirect so the pipes buffer does not fill and stall 
        # everything
        self.exonerate.stdout = open('/dev/null', 'w')

    def stop(self) :
        if self.exonerate :
            self.exonerate.terminate()

    def _get_genes(self, fname) :
        tmp = {}
        f = FastqFile(fname)
        f.open()

        for seq in f :
            print seq.id
            tmp[seq.id[1:]] = None

        f.close()
        return tmp

    def query_sequence(self, seq) :
        return self.query_sequences([seq])

    def query_sequences(self, seqlist) :
        f = open(tempfile.mktemp(), 'w')

        for seq in seqlist :
            print >> f, seq.fasta()

        f.close()

        return self.query_file(f.name, dict(zip(map(lambda x : x.id, seqlist), [None]*len(seqlist))), True)

    def query_file(self, fname, seqdict=None, unlink=False) :
        oldargs = ['exonerate', fname,
                'localhost:12887',
                '--bestn', '1',
                '--model', 'affine:local',
                '--showalignment', 'no']

        args = ['exonerate', fname,
                'localhost:12887',
                '--bestn', '1',
                '--showalignment', 'no']

        if seqdict is None :
            seqdict = self._get_genes(fname)

        try :
            out = subprocess.check_output(
                            args,
                            stderr=open('/dev/null', 'w'))

        except subprocess.CalledProcessError, cpe :
            if unlink :
                os.remove(fname)

            return seqdict

        for line in out.split('\n') :
            if not line.startswith('vulgar') :
                continue

            #print line.rstrip()

            # vulgar: aael001690 581 737 + >agap005310 596 752 + 339 M 156 156
            # vulgar: aael009244 667 744 + >agap006707 604 681 + 178 M 77 77
            fields = line.strip().split()
            #seqdict[fields[1]] = fields[5][1:]
            seqdict[fields[1]] = ExonerateHit(
                                    ExonerateVulgar(fields[1],     int(fields[2]), int(fields[3])),
                                    ExonerateVulgar(fields[5][1:], int(fields[6]), int(fields[7]))
                                 )

        if unlink :
            os.remove(fname)

        return seqdict

