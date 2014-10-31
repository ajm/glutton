from sys import exit

from glutton.utils import get_log
from glutton.info import GluttonInformation


class Scaffolder(object) :
    def __init__(self, db, contigs_fname, alignments_dir, scaffolds_fname, identity) :
        self.db                 = db
        self.contigs_fname      = contigs_fname
        self.alignments_dir     = alignments_dir
        self.scaffolds_fname    = scaffolds_fname
        self.identity           = identity

        self.log = get_log()

        self.info = GluttonInformation(db, contigs_fname, alignments_dir)

        if not self.info.alignments_complete() :
            self.log.fatal("alignments are not complete!")
            exit(1)

    def stop(self) :
        pass

    def scaffold(self) :
        pass

