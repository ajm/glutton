import sys
import signal
import os
import subprocess
import shutil

from glutton.base import Base
from glutton.ensembl import EnsemblInfo
from glutton.local import LocalInfo
from glutton.transcriptome import Transcriptome
from glutton.queue import WorkQueue
from glutton.scaffolder import Scaffolder


class Subcommand(Base) :
    def __init__(self, options, parameters=[], programs=[], local=True, list_option=True) :
        super(Subcommand, self).__init__(options)

        self.parameters = parameters
        self.programs = programs
        self.local = local

        if list_option and self.opt['list'] :
            info = LocalInfo(self.opt) if self.local else EnsemblInfo(self.opt)
            info.print_species_table()
            sys.exit(0)

        self._check_programs()
        self._check_parameters()
        self._check_species()

    def _program_exists(self, programname) :
        for p in os.environ['PATH'].split(os.pathsep) :
            progpath = os.path.join(p, programname)
            if os.path.isfile(progpath) :
                # there may be another executable with the correct
                # permissions lower down in the path, but the shell
                # would not find it, so just return here...
                return os.access(progpath, os.X_OK)

        return False

    def _check_programs(self) :
        bad = False

        # ensure required programs are installed locally
        for prog in self.programs :
            if not self._program_exists(prog) :
                self.error("'%s' is not installed." % prog)
                bad = True

        if bad :
            sys.exit(1)

    def _check_parameters(self) :
        bad = False

        for p in self.parameters :
            if not self.opt[p] :
                self.error("missing mandatory command line argument '--%s'" % p)
                bad = True

        if bad :
            sys.exit(1)

    def _check_species(self) :
        if 'species' not in self.parameters :
            return

        info = LocalInfo(self.opt) if self.local else EnsemblInfo(self.opt)
        species = self.opt['species']
        release = self.opt['release']
        fill_in = ('species' in self.parameters) and ('release' not in self.parameters)

        if not info.is_valid_species(species) :
            self.error("'%s' is an invalid species" % species)
            sys.exit(1)

        if release is None and fill_in :
            release = info.get_latest_release(species)

            if release == -1 :
                self.error("could not find database for \'%s\'" % species)
                sys.exit(1)

            self.opt['release'] = release
            self.info("release not specified, using release %d" % release)

        else :
            if not info.is_valid_release(species, release) :
                self.error("%d is an invalid release for '%s'" % (release, species))
                sys.exit(1)

    def run(self) :
        return self._run()

    def _run(self) :
        raise NotImplementedError

class BuildCommand(Subcommand) :
    parameters = ['species','database']
    programs = ['prank', 'fastareformat', 'fasta2esd', 'esd2esi']

    def __init__(self, opt) :
        super(BuildCommand, self).__init__(opt, self.parameters, self.programs, local=False)

        self.transcriptome = Transcriptome(opt, WorkQueue(opt, opt['threads']))
        self._init()

    def _init(self) :
        def _cleanup(signal, frame) :
            if not self.transcriptome.is_cancelled() :
                self.transcriptome.stop()
                self.warn("Application will shutdown after current download or user types '^C' again...")
                return

            self.warn("Forced exit from user")
            os._exit(0)

        signal.signal(signal.SIGINT, _cleanup)

    def _run(self) :
        self.transcriptome.build()
        return 0

class ListCommand(Subcommand) :
    def __init__(self, opt) :
        super(ListCommand, self).__init__(opt)

    def _run(self) :
        EnsemblInfo(self.opt).print_species_table()
        LocalInfo(self.opt).print_species_table()
        return 0

class AlignCommand(Subcommand) :
    parameters = ['species', 'input-file', 'output-dir', 'min-length']
    programs = ['pagan', 'exonerate-server', 'exonerate']
    
    def __init__(self, opt) :
        super(AlignCommand, self).__init__(opt, self.parameters, self.programs)

        self.transcriptome = Transcriptome(opt, WorkQueue(opt, opt['threads']), skip_checks=True)
        self._init()

    def _init(self) :
        def _cleanup(signal, frame) :
            self.error("Forced exit from user")
            os._exit(0)

        signal.signal(signal.SIGINT, _cleanup)

    def _run(self) :
        self.transcriptome.align(
                            self.opt['input-file'], 
                            self.opt['output-dir'], 
                            min_length=self.opt['min-length'])
        return 0

class ScaffoldCommand(Subcommand) :
    parameters = ['input-file', 'output-dir', 'min-identity', 'output-file']
    programs = []

    def __init__(self, opt) :
        super(ScaffoldCommand, self).__init__(opt, self.parameters, self.programs)

    def _run(self) :
        Scaffolder(self.opt).scaffold(
                            self.opt['input-file'],
                            self.opt['output-dir'],
                            self.opt['output-file'],
                            self.opt['min-identity'])
        return 0

class FixCommand(Subcommand) :
    parameters = ['species','release']
    programs = ['fastareformat', 'fasta2esd', 'esd2esi']

    def __init__(self, opt) :
        super(FixCommand, self).__init__(opt, self.parameters, self.programs)

        self.transcriptome = Transcriptome(opt, skip_checks=True)

    def _run(self) :
        self.transcriptome.fix()
        return 0

class RmCommand(Subcommand) :
    parameters = ['species', 'release']

    def __init__(self, opt) :
        super(RmCommand, self).__init__(opt, self.parameters)

    def _run(self) :
        self.info("deleting %s:%d ..." % (self.opt['species'], self.opt['release']))

        try :
            shutil.rmtree(os.path.join(self.dbdir, str(self.opt['release']), self.opt['species']))
        
        except OSError, ose :
            self.error(str(ose))
            self.error('halting')
            return 1

        self.info('deleted!')
        return 0

class PackCommand(Subcommand) :
    def __init__(self, opt) :
        super(PackCommand, self).__init__(opt, parameters=['species', 'release'])

        self.transcriptome = Transcriptome(opt, skip_checks=True)

    def _run(self) :
        return self.transcriptome.package()

class UnpackCommand(Subcommand) :
    def __init__(self, opt) :
        super(UnpackCommand, self).__init__(opt, parameters=['input-file'], list_option=False)
        self._init()

    def _init(self) :
        name,ext = os.path.splitext(os.path.basename(self.opt['input-file']))
   
        if ext != '.tgz' :
            self.error("do not know how to deal with '%s' files" % ext)
            sys.exit(1)
   
        try :
            species, release = name.split('_')
            self.opt['species'] = species
            self.opt['release'] = int(release)

        except ValueError, ve:
            self.error("package files must be called 'species_release.tgz'")
            sys.exit(1)

        self.transcriptome = Transcriptome(self.opt)

    def _run(self) :
        return self.transcriptome.unpackage(self.opt['input-file'])         

