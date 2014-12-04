from sys import exit, stderr, stdout
from glob import glob
from os.path import join
from collections import defaultdict, Counter

import re
import operator

from Bio import SeqIO

from glutton.utils import get_log
from glutton.info import GluttonInformation


NOT_DEALT_WITH = 0

class ScaffolderError(Exception) :
    pass

scaffold_counter = -1
scaffold_fmt = "glutton%d"

def scaffold_id() :
    global scaffold_fmt, scaffold_counter

    scaffold_counter += 1

    return scaffold_fmt % scaffold_counter

class Alignment(object) :
    def __init__(self, id, gene_name, start, end, seq, desc='singleton', contigs=None) :
        self.id = id
        self.gene_name = gene_name
        self.start = start
        self.end = end
        self.seq = seq
        self.desc = desc

        if contigs :
            self.contigs = contigs
        else :
            self.contigs = [ self.id ]

        # XXX when being read from the file this should be verified
        #     so i don't need to check then fail
        m = re.match("^(comp\d+\_c\d+)\_seq\d+$", self.id)
        self.gene_id = m.group(1)

    def _ensure_order(self, a2, func) :
        return func(self, a2) if self.start <= a2.start else func(a2, self)

    def extract(self, start, end) :
        return self.seq[start:end]

    # assumes a starts before b
    def _overlaps(self, a, b, consider_isoforms=True) :
        if consider_isoforms and (a.gene_id == b.gene_id) :
            return None

        elif a.end > b.end :
            return (b.start, b.end)

        elif a.end <= b.start :
            return None

        elif a.end > b.start :
            return (b.start, a.end)

        elif a.end == b.start :
            return (b.start, a.end)

        raise ScaffoldError("impossible condition reached") 

    def overlaps(self, a2) :
        return self._ensure_order(a2, self._overlaps)

    def _simple_overlaps(self, a, b) :
        return self._overlaps(a, b, consider_isoforms=False)

    def simple_overlaps(self, a2) :
        return self._ensure_order(a2, self._simple_overlaps)

    def _remove_common_gaps(self, seq_a, seq_b) :
        new_a = ""
        new_b = ""

        for a,b in zip(seq_a, seq_b) :
            if (a,b) == ('-','-') :
                continue

            new_a += a
            new_b += b

        return new_a, new_b

    def _mergeable(self, a, b, consider_isoforms=True) :
        overlap = self._overlaps(a, b, consider_isoforms=consider_isoforms) 

        if not overlap :
            return False

        slice_a, slice_b = self._remove_common_gaps(a.extract(*overlap), b.extract(*overlap))

        for a,b in zip(slice_a, slice_b) :
            if a == '-' or b == '-' :
                return False

        return True

    def mergeable(self, a2) :
        return self._ensure_order(a2, self._mergeable)

    def _simple_mergeable(self, a, b) :
        return self._mergeable(a, b, consider_isoforms=False)

    def simple_mergeable(self, a2) :
        return self._ensure_order(a2, self._simple_mergeable)

    # assumes a starts before b
    def _merge_contigs(self, a, b) :
        # continuous case
        if a.end > b.start :
            return a.seq[:b.start] + b.seq[b.start:b.end] + a.seq[b.end:]
        # discontinuous case
        else :
            return a.seq[:a.end] + ('N' * (b.start - a.end)) + b.seq[b.start:]

    def merge_contigs(self, a2) :
        return self._ensure_order(a2, self._merge_contigs)

    def __add__(self, a2) :
        return Alignment(self.id, 
                        self.gene_name, 
                        min(self.start, a2.start), 
                        max(self.end, a2.end), 
                        self.merge_contigs(a2), 
                        'merged',
                        self.contigs + a2.contigs)

    def __iadd__(self, a2) :
        self.start = min(self.start, a2.start)
        self.end = max(self.end, a2.end)
        self.seq = self.merge_contigs(a2)
        self.contigs += a2.contigs
        self.desc = 'merged'

    def get_desc(self) :
        return self.desc
        #if len(self.contigs) == 1 :
        #    return 'singular'
        #else :
        #    return 'fragmented'

    def format(self, fmt='fasta') :
        assert fmt == 'fasta' # nothing else supported yet

        return ">%s contigs=%s gene=%s desc=%s\n%s" % (\
            scaffold_id(), \
            ','.join(self.contigs), \
            self.gene_name, \
            self.desc, \
            self.seq.replace('-', '')\
          )

    def __str__(self) :
        return str((self.id, self.start, self.end, self.seq))

class Scaffolder(object) :
    def __init__(self, db, contigs_fname, alignments_dir, scaffolds_fname, identity) :
        self.db                 = db
        self.contigs_fname      = contigs_fname
        self.alignments_dir     = alignments_dir
        self.scaffolds_fname    = scaffolds_fname
        self.identity           = identity

        self.log = get_log()

        self.info = GluttonInformation(db, contigs_fname, alignments_dir)
        
        self.scaffold_fmt = "glutton%d"
        self.scaffold_id = 0

#        if not self.info.alignments_complete() :
#            self.log.fatal("alignments are not complete!")
#            exit(1)

        # e.g. query39806_orf1 [2.2.1363]
        #self.pagan_orfname_regex = re.compile("(.*)\_orf(\d+) \[(\-?\d+)\.(\d+)\.(\d+)\]")
        # e.g. query39806_orf1
        self.orfname_regex = re.compile("^(query\d+)\_orf(\d)$")
        # e.g. comp18591_c1_seq4
        self.trinity_regex = re.compile("^(comp\d+\_c\d+)\_seq\d+$") # just group the trinity gene id, seqX is the isoform id

    def stop(self) :
        pass

    def _orf_to_query_name(self, name) :
        m = self.orfname_regex.match(name)

        if not m :
            raise ScaffolderError("unexpected query name (%s)" % name)

        return m.group(1)

#    def read_contigids_in_alignments(self) :
#        tmp = set()
#
#        for f in glob(join(self.alignments_dir, "glutton*")) :
#            names = []
#
#            for s in SeqIO.parse(f, 'fasta') :
#                if s.description.startswith('query') :
#                    try :
#                        names.append(self._orf_to_query_name(s.description))
#                    
#                    except ScaffolderError, se :
#                        self.log.error(str(se))
#
#            tmp.update(self.info.get_contig_from_query(names))
#            
#        return tmp

    # read the alignment file and return a dictionary keyed on the gene name
    #   - there might be multiple orfs for a single query sequence, so keep a track of the best one
    def read_alignment(self, fname) :
        tmp = defaultdict(dict)

        # identity here is the proportion of non-gaps
        def calc_identity(seq1, seq2, start, end) :
            count = 0
            for i,j in zip(seq1[start:end], seq2[start:end]) :
                if (i,j) == ('-','-') :
                    continue
                if (i != '-') and (j != '-') :
                    count += 1

            return count / float(end - start)

        for s in SeqIO.parse(fname, 'fasta') :
            if not s.description.startswith('query') :
                gene_id  = s.description
                gene_seq = str(s.seq)
                continue

            query_id  = self._orf_to_query_name(s.description)
            contig_id = self.info.get_contig_from_query(query_id)
            gene_name = self.db.get_genename_from_geneid(gene_id)

            seq = str(s.seq).replace('N', '-')
            
            contig_start = len(seq) - len(seq.lstrip('-'))
            contig_end   = contig_start + len(seq.strip('-'))

            # if we have seen this before
            if contig_id in tmp[gene_name] :
                new = calc_identity(gene_seq, seq, contig_start, contig_end)
                old = calc_identity(gene_seq, *tmp[gene_name][contig_id])

                if new < old :
                    continue

            # i have seen duplicate sequences before
            dup = False
            for details in tmp[gene_name].values() :
                if details[0] == seq :
                    dup = True
                    break

            if dup :
                continue

            #print gene_name, contig_id
            tmp[gene_name][contig_id] = (seq, contig_start, contig_end)

        # convert from a dict of dicts to a dict of lists
        tmp2 = defaultdict(list)

        #print fname
        for gene in tmp :
            #print "\t", gene
            for contig in tmp[gene] :
                #print "\t\t", contig
                seq,start,end = tmp[gene][contig]
                tmp2[gene].append(Alignment(contig, gene, start, end, seq))

        return tmp2

    def remove_common_gaps(self, seq_a, seq_b) :
        new_a = ""
        new_b = ""

        for a,b in zip(seq_a, seq_b) :
            if (a,b) == ('-','-') :
                continue
            
            new_a += a
            new_b += b

        return new_a, new_b

    def get_alignment_overlaps(self, alignments) :
        tmp = []

        for ind_a in range(0, len(alignments)) :
            for ind_b in range(ind_a + 1, len(alignments)) :
                a = alignments[ind_a]
                b = alignments[ind_b]
                
                overlap = a.simple_overlaps(b)

                if overlap :
                    tmp.append((ind_a, ind_b, overlap))

        return tmp

    # a complex overlap set is defined as one where those overlaps
    # form a non-linear graph
    #
    #   a    b      (simple)
    #   a----b----c (simple)
    #   a----b c----d (simple)
    #
    #   a-------b   (complex)
    #    \
    #     ------c
    #

    def group_alignments(self, alignments) :
        groups = []

        # for each alignment
        for index,align in enumerate(alignments) :
            # for each group
            add_to_group = False
            for group in groups :
                # if 'align' overlaps with a member of that group
                for member in group :
                    if member.simple_overlaps(align) :
                        add_to_group = True
                        break

                # add to the group
                if add_to_group :
                    group.append(align)
                    break

            if not add_to_group :
                groups.append([ align ])

        #print "#alignments=%d #groups=%d (%s)" % (len(alignments), len(groups), ','.join([ str(len(i)) for i in groups ]))

        return groups

    def group_cannot_be_merged_isoforms(self, group) :
        return self.group_cannot_be_merged(group, consider_isoforms=True)

    def group_cannot_be_merged(self, group, consider_isoforms=False) :
        for i in range(0, len(group)) :
            a = group[i]
            for j in range(i+1, len(group)) :
                b = group[j]

                if a.simple_overlaps(b) :
                    if consider_isoforms :
                        if not a.mergeable(b) :
                            return True
                    else :
                        if not a.simple_mergeable(b) :
                            return True

        return False

    def merge_alignments(self, alignments) :
        # perform merges of overlaps first
        complex_group = False
        merged_groups = []


        for group in self.group_alignments(alignments) :
            if self.group_cannot_be_merged(group) :
                
                #self.print_alignments(group)
                
                complex_group = True
                for i in group :
                    i.desc = 'conflict'
                    merged_groups.append(i)

            elif self.group_cannot_be_merged_isoforms(group) :
                #self.print_alignments(group)

                complex_group = True

                # TODO distinguish between groups comprised wholly of isoforms for one gene 
                #      vs a group with other stuff mixed in (isoforms_multiple vs isoforms_ambiguous)

                if len(set([ g.gene_id for g in group ])) == 1 :
                    d = 'single_gene_multiple_isoform'
                else :
                    d = 'multiple_gene_multiple_isoform'

                for i in group :
                    i.desc = d
                    merged_groups.append(i)

            else :
                merged_groups.append( reduce(operator.add, sorted(group, key=lambda x : x.start)) )


        # bail out if there are any complex alignments (non-linear overlaps caused by bad merges or isoforms)
        if complex_group :
            return merged_groups

        # otherwise merge fragments with N's
        return [ reduce(operator.add, sorted(merged_groups, key=lambda x : x.start)) ]

    def print_alignments(self, alignments) :
        if len(alignments) == 0 :
            return

        # many isoforms, same gene
        num_genes = len(set([ '_'.join(i.id.split('_')[:1]) for i in alignments ])) 
        if num_genes == 1 :
            return

        overlaps = self.get_alignment_overlaps(alignments)

        if len(overlaps) == 0 :
            return

        print ""
        print "#alignments=%d #overlaps=%d #genes=%d" % (len(alignments), len(overlaps), num_genes)
        print [ (i.id, i.start, i.end) for i in alignments ]
        print overlaps
        #print ""

        #raise NotImplementedError("")

        # 3. trivial overlaps - no overlap conflicts
        # pairs of contigs overlapping + singletons
        
        for contig_a,contig_b,overlap_indices in overlaps :

            sl = slice(*overlap_indices)
            slice_a = alignments[contig_a].seq[sl]
            slice_b = alignments[contig_b].seq[sl]

            slice_a, slice_b = self.remove_common_gaps(slice_a, slice_b)
            slice_diff = ''.join([ "X" if i != j else " " for i,j in zip(slice_a, slice_b) ])

            if ('-' in slice_a) or ('-' in slice_b) :
                print (contig_a,contig_b,overlap_indices), "!="
            else :
                print (contig_a,contig_b,overlap_indices), "=="

            print slice_a
            print slice_b
            print slice_diff

    def process_alignments(self, fout) :

        for fname in glob(join(self.alignments_dir, 'glutton*.nucleotide')) :
            alignment = self.read_alignment(fname)

            #print fname
            #print alignment

            for gene_name in alignment :
                for a in self.merge_alignments(alignment[gene_name]) :
                    print >> fout, a.format('fasta') 


    def next_scaffold_id(self) :
        self.scaffold_id += 1
        return self.scaffold_fmt % (self.scaffold_id - 1)

    # return a string in the glutton format
    def fasta_output(self, contigs, gene, seq, desc) :
        return ">%s contigs=%s gene=%s desc=%s\n%s" % (self.next_scaffold_id(), ','.join(contigs), gene, desc, seq)

    def scaffold(self) :
        global NOT_DEALT_WITH

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

        with open(self.scaffolds_fname, 'w') as f :
            aligned_contigs = self.process_alignments(f)
            #self.output_unscaffolded_contigs(f, aligned_contigs)

        self.log.info("wrote %s" % self.scaffolds_fname)
        self.log.info("ignored %d alignments (too complex)" % NOT_DEALT_WITH)

        return 0

    def output_unscaffolded_contigs(self, fout, aligned_contigs) :
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

