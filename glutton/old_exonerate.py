import sys
import subprocess
import tempfile
import os

from glutton.base import ToolBase, Base


class ExonerateError(Exception) :
    pass

class ExonerateServer(Base) :
    def __init__(self, opt, dir, db_name, port=12887) :
        super(ExonerateServer, self).__init__(opt)

        self.dir = dir
        self.db_name = db_name
        self.exonerate = None
        self.port = port

        if not os.path.isfile(os.path.join(self.dir, self.db_name + '.esi')) :
            raise ExonerateError("database called '%s' not found in dir '%s'" % (self.db_name, self.dir))

    def __del__(self) :
        self.stop()

    def _wait_until_listening(self) :
        out = []
        
        for line in iter(self.exonerate.stdout.readline, '') :
            out.append(line)

            if 'listening on port' in line :
                self.info("started...")
                break

        else :
            self.error("did not start properly:\n\n%s\n" % '\n'.join(out))
            sys.exit(1)

    def started(self) :
        return self.exonerate != None

    def start(self) :
        self.info("starting...")

        args = ['exonerate-server', self.db_name + '.esi',
                '--port', str(self.port)]

        subprocess.call(['killall', 'exonerate-server'], 
                        stdout=open('/dev/null', 'w'), 
                        stderr=subprocess.STDOUT)

        self.exonerate = subprocess.Popen(args,
                                    cwd=self.dir,
                                    stderr=subprocess.STDOUT,
                                    stdout=subprocess.PIPE)
        
        self._wait_until_listening()

        # change the redirect so the pipes buffer does not fill and stall everything
        self.exonerate.stdout = open('/dev/null', 'w')

    def stop(self) :
        if self.exonerate :
            self.info("stopping...")
            self.exonerate.terminate()
            self.info("stopped")

    def query(self, fname) :
        args = ['exonerate',
                '--bestn', '10',
                '--model', 'affine:local',
                '--showalignment', 'no',
                fname]

        if self.started() :
            args.append('localhost:' + str(self.port))
        else :
            args.append(self.db_name + '.fa')

        try :
            output = subprocess.check_output(args, 
                                             cwd=self.dir, 
                                             stderr=open('/dev/null', 'w'))

        except subprocess.CalledProcessError, cpe :
            #self.error("exonerate did not run properly, returncode = %d\n\n%s" % (cpe.returncode, cpe.output))
            #sys.exit(1)
            raise ExonerateError("exonerate abnormal termination : returncode = %d" % cpe.returncode)

        hits = []
        high_score = -1

        # field 5 is gene name
        # field 9 is the 'raw score'
        for line in output.split('\n') :
            if line.startswith('vulgar') :
                fields = line.strip().split()
                #return fields[5]

                gene_name = fields[5]
                raw_score = int(fields[9])                

                if high_score == -1 :
                    high_score = raw_score

                if high_score == raw_score :
                    hits.append(gene_name)

            # vulgar: aael001690 581 737 + >agap005310 596 752 + 339 M 156 156
            # vulgar: aael009244 667 744 + >agap006707 604 681 + 178 M 77 77

        if len(hits) == 0 :
            raise ExonerateError("no alignment")

        return hits

class Fastareformat(ToolBase) :
    def __init__(self, opt, fasta_file) :
        super(Fastareformat, self).__init__(opt)
        self.fasta_file = fasta_file

    def run(self) :
        returncode, output = self._execute([self.fasta_file], [])

        if returncode == 0 :
            self._create_file(
                    os.path.basename(self.fasta_file), 
                    output, 
                    os.path.dirname(self.fasta_file)
                    )

        return returncode

class Fasta2esd(ToolBase) :
    def __init__(self, opt, fasta_file, db_file) :
        super(Fasta2esd, self).__init__(opt)
        self.fasta_file = fasta_file
        self.db_file = db_file

    def run(self) :
        returncode, output = self._execute([self.fasta_file, self.db_file], [self.db_file])
        return returncode

class Esd2esi(ToolBase) :
    def __init__(self, opt, db_file, index_file) :
        super(Esd2esi, self).__init__(opt)
        self.db_file = db_file
        self.index_file = index_file

    def run(self) :
        returncode, output = self._execute([self.db_file, self.index_file], [self.index_file])
        return returncode

class ExonerateDBBuilder(Base) :
    def __init__(self, opt, dir, db_name) :
        super(ExonerateDBBuilder, self).__init__(opt)
        self.dir = dir
        self.db_name = db_name

    def filenames(self) :
        return [self.db_name + '.' + ext for ext in ['fa','esi','esd']]

    def build(self, iterable) :
        fasta_file = os.path.join(self.dir, self.db_name + '.fa')
        db_file    = self._change_ext(fasta_file, 'esd')
        index_file = self._change_ext(fasta_file, 'esi')

        self.rm([fasta_file, db_file, index_file])

        self._create_fasta(iterable, fasta_file)

        if Fastareformat(self.opt, fasta_file).run() != 0 :
            return False

        if Fasta2esd(self.opt, fasta_file, db_file).run() != 0 :
            return False

        if Esd2esi(self.opt, db_file, index_file).run() != 0 :
            return False

        return True

    def _create_fasta(self, iterable, fname) :
        f = open(fname, 'w')
        
        genes = set()

        for gf in iterable :
            for k,v in gf.iteritems() :
                if k not in genes :
                    print >> f, ">%s\n%s" % (k,v)
                    genes.add(k)

        f.close()

