from sys import exit
from glob import glob
from os.path import join
import re

from Bio import SeqIO

from glutton.utils import get_log
from glutton.info import GluttonInformation


class ScaffolderError(Exception) :
    pass

class Scaffolder(object) :
    def __init__(self, db, contigs_fname, alignments_dir, scaffolds_fname, identity) :
        self.db                 = db
        self.contigs_fname      = contigs_fname
        self.alignments_dir     = alignments_dir
        self.scaffolds_fname    = scaffolds_fname
        self.identity           = identity

        self.log = get_log()

        self.info = GluttonInformation(db, contigs_fname, alignments_dir)
        self.scaffold_id = 0

        if not self.info.alignments_complete() :
            self.log.fatal("alignments are not complete!")
            exit(1)

        self.pagan_orfname_regex = re.compile("(.*)\_orf(\d) \[(\-?\d+)\.(\d+)\.(\d+)\]")

    def stop(self) :
        pass

    def _pagan_to_contig_name(self, name) :
        m = self.pagan_orfname_regex.match(name)

        if not m :
            raise ScaffolderError("unexpected pagan name (%s)" % name)

        return m.group(1)

    # we need to do this because an alignment can partially fail
    # e.g.  align two contigs, but only one of them succeeds and the
    #       other is silently(?) dropped
    def read_contigids_in_alignments(self) :
        tmp = set()

        for f in glob(join(self.alignments_dir, "glutton*")) :
            for s in SeqIO.parse(f, 'fasta') :
                if s.description.startswith('query') :
                    try :
                        names.append(self._pagan_to_contig_name(s.description))
                    
                    except ScaffolderError, se :
                        self.log.error(str(se))

            tmp.update(self.info.get_contig_from_query(names))
            
        return tmp

    def fasta_output(self, id, seq, status) :
        return id, status
        #return ">%s %s\n%s" % (id, status, seq)

    def scaffold(self) :
        # previously it was one contig per family, however, they are now performed simultaneously
        #   hence nothing really needs to be stored
        #
        # how do I know which contigs did not get aligned? these need to be output as well 
        #   with a comment
        #
        #       filtered
        #       no gene family assignment
        #       alignment failed
        #
        # ALL
        #   replace query name with original contig name
        #   append as a comment in the fasta file the original gene names (comma separated)
        #
        # SINGLE CONTIG
        #   just output as is
        #
        # NON OVERLAPPING CONTIGS
        #   output combined with gaps or Ns
        #
        # UNAMBIGUOUS OVERLAP
        #   - combine, but indicate somehow...
        #       include all original ids, + give new one?
        #
        # AMBIGUOUS OVERLAP
        #   - do research!
        #   - what do these look like?
        #   - is there a simple way to change the originals
        # 
        
        self.log.info("starting scaffolding")

        fout = open(self.scaffolds_fname, 'w')
        
        # scaffold the results from pagan alignments
        # read what contigs were used
        aligned_contigs = self.read_contigids_in_alignments()
        self.log.info("%d contigs in alignment files" % len(aligned_contigs))


        # output the unscaffolded contigs
        #
        for r in SeqIO.parse(self.contigs_fname, 'fasta') :

            if r.id in aligned_contigs :
                continue

            # unused due to being filtered out
            if not self.info.contig_used(r.id) :
                print >> fout, self.fasta_output(r.id, r.seq, 'filtered')
            
            # unassigned by blast
            elif not self.info.contig_assigned(r.id) :
                print >> fout, self.fasta_output(r.id, r.seq, 'assignment_failed')
            
            # unaligned by pagan
            #elif not self.info.contig_aligned(r.id) :
            elif r.id not in aligned_contigs :
                print >> fout, self.fasta_output(r.id, r.seq, 'alignment_failed')
            
            else :
                self.log.warning("contig cannot be accounted for: %s" % r.id)
                print >> fout, self.fasta_output(r.id, r.seq, 'unaccounted_for')

        fout.close()
        return 0

