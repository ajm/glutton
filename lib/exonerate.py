import sys
import subprocess
import tempfile
import os

from lib.base import ToolBase, Base


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
                self.info("exonerate-server started...")
                break
        else :
            if 'Could not open esd file' in output :
                self.error("exonerate-server did not start properly:" + 
                           "\n\tcould not read esd file, if it was build on another OS run '%s fix -s '%s' -r %d -d '%s' to rebuild esd file'" % 
                           (sys.argv[0], self.opt['species'], self.opt['release'], self.opt['database']))
            else :
                self.error("exonerate-server did not start properly:\n\n%s\n" % '\n'.join(out))

            sys.exit(1)

    def started(self) :
        return self.exonerate != None

    def test(self) :
        test_exonerate = subprocess.Popen(['exonerate-server', self.db_name + '.esi', '--port', str(self.port)],
                                          cwd=self.dir, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        tmp = False

        for line in iter(test_exonerate.stdout.readline, '') :
            if 'listening on port' in line :
                tmp = True
                break

            if 'Could not open esd file' in line :
                tmp = False
                break

        test_exonerate.terminate()
        return tmp

    def start(self) :
        self.info("exonerate-server starting...")

        args = ['exonerate-server', self.db_name + '.esi',
                '--port', str(self.port)]

        self.exonerate = subprocess.Popen(args,
                                    cwd=self.dir,
                                    stderr=subprocess.STDOUT,
                                    stdout=subprocess.PIPE)
        
        if not self._wait_until_listening() :
            sys.exit(1)

        # change the redirect so the pipes buffer does not fill and stall everything
        self.exonerate.stdout = open('/dev/null', 'w')

    def stop(self) :
        if self.exonerate :
            self.info("exonerate-server stopping...")
            self.exonerate.terminate()
            self.info("exonerate-server stopped")

    def query(self, fname) :
        args = ['exonerate', fname,
                'localhost:' + str(self.port),
                '--bestn', '1',
                '--model', 'affine:local',
                '--showalignment', 'no']

        try :
            output = subprocess.check_output(args, stderr=open('/dev/null', 'w'))

        except subprocess.CalledProcessError, cpe :
            self.error("exonerate did not run properly, returncode = %d\n\n%s" % (cpe.returncode, cpe.output))
            sys.exit(1)

        for line in output.split('\n') :
            if line.startswith('vulgar') :
                fields = line.strip().split()
                return fields[5]

            # vulgar: aael001690 581 737 + >agap005310 596 752 + 339 M 156 156
            # vulgar: aael009244 667 744 + >agap006707 604 681 + 178 M 77 77

        raise ExonerateError("no alignment")

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

