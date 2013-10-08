import sys
import re
import operator

from os.path import join
from glob import glob
from itertools import izip

from glutton.base import Base
from glutton.filetypes import FastqFile
from glutton.datatypes import Sequence


class MergeError(Exception) :
    pass


class AlignmentRange(object) :
    def __init__(self, start, end, name, seq) :
        self.start = start
        self.end = end
        self.name = name
        self.seq = seq

    # range subsumes 'other'
    def __contains__(self, other) :
        return self.subsumes(other)

    def subsumes(self, other) :
        return (self.start <= other.start) \
                and (self.end >= other.end)

    # ranges overlap
    def __and__(self, other) :
        return self.overlaps(other)

    def overlaps(self, other) :
        return (self.end >= other.start) \
                if (self.start <= other.start) \
                else (other.end >= self.start)

    def __merge_name(self, other) :
        return "%s+%s" % (self.name, other.name)

    def __merge_overlapping(self, other) :
        if other.seq in self.seq :
            return AlignmentRange(self.start, self.end, self.__merge_name(other), self.seq)

        min_overlap = self.end - other.start
        max_overlap = len(other)
        overlap = min_overlap

        seq1 = self.seq
        seq2 = other.seq

        while overlap <= max_overlap :
            if seq1.endswith(seq2[:overlap]) :
                return AlignmentRange(self.start, other.end, self.__merge_name(other), seq1 + seq2[overlap:])
            
            overlap += 1

        
        raise MergeError()       

    def __merge_nonoverlapping(self, other) :
        gap = other.start - self.end
        seq = self.seq + ("N" * gap) + other.seq

        return AlignmentRange(self.start, other.end, self.__merge_name(other), seq)

    # __add__ assumes that self is on the left and other is on the right
    #
    # e.g.  ----left--------------
    #                       ----right----
    # or 
    #       ----left----
    #                       ----right----
    def __add__(self, other) :
        if self.overlaps(other) :
            return self.__merge_overlapping(other)
        else :
            return self.__merge_nonoverlapping(other)

    def __lt__(self, other) :
        return self.end < other.end

    def __len__(self) :
        return self.end - self.start

    def __str__(self) :
        return ">%s\n%s" % (self.name, self.seq)


class Scaffolder(Base) :
    def __init__(self, opt) :
        super(Scaffolder, self).__init__(opt)

        self.pagan_pattern = re.compile("(.*)(\[\-?\d+\ \d+\-\d+\]|\.orf\.\-?\d+\.\d+\.\d+)")
        self.out_file = None

    def __get_protein_alignments(self, d) :
        return [ fn for fn in glob(join(d,'*')) if fn.endswith('.fas') and not fn.endswith('.dna.fas') ]

    def __get_dna_alignments(self, d) :
        return [ fn for fn in glob(join(d,'*')) if fn.endswith('.dna.fas') ]

    def __get_queries(self, fn, allow_empty=False) :
        query = None
        reference = None
        found = 0

        f = FastqFile(fn)
        f.open()

        for seq in f :
            if seq.id.startswith('comp') :
                if not reference :
                    self.error("found query before reference (%s)" % fn)
                    sys.exit(1)

                query = seq
                found += 1

                yield (query, reference)
            else :
                reference = seq

        f.close()

        if found == 0 and not allow_empty:
            self.error("no query in alignment in file (%s)" % fn)
            sys.exit(1)

    def __get_alignment_columns(self, query, reference) :
        alignment_start = len(query.sequence) - len(query.sequence.lstrip('-'))
        alignment_end = alignment_start + len(query.sequence.strip('-'))

        ref_alignment = reference.sequence[alignment_start:alignment_end]
        qry_alignment = query.sequence[alignment_start:alignment_end]

        tmp = [(i,j) for i,j in izip(ref_alignment, qry_alignment) if not ((i == '-') and (j == '-'))]

        id_char = len([ q for r,q in tmp if (r != '-') and (q != '-') ])
        identity = id_char / float(len(ref_alignment))

        # XXX need to correct alignment_end because I want it to be the same as the original 
        #     PRANK alignment for all contigs
        #     alignment_end = alignment_end - total gaps in ref gene not in query contig
        true_alignment_end = alignment_end - len([r for r,q in tmp if (r == '-') and (q != '-')])

        return (alignment_start, true_alignment_end, identity)

    def __merge(self, alignment_ranges) :
        try :
            return reduce(operator.add, sorted(alignment_ranges, key=lambda x : x.start))

        except MergeError, me :
            return None

    def __munge_contigname(self, name) :
        return name.split()[0]

    # needs to be regex to account for ".orf" and "[x" names
    def __munge_paganorfname(self, name) :
        m = self.pagan_pattern.match(name)
        
        if not m :
            self.error("pagan sequence name '%s' not matched by regex" % name)
            sys.exit(1)
        
        return m.group(1)
        #return name.split('.')[0]

    # contig ids actually get mangled by PAGAN (states which orf succeeded, i think)
    # so i need to map it back to the original name at some point
    #
    #   pagan    : >comp0_c0_seq1.orf.2.1.1779 or >comp10099_c0_seq1[-3 2-682]
    #   original : >comp0_c0_seq1 len=395 path=[184:0-140 390:141-394]
    #
    def __read_contigs(self, fn, gene2range) :
        query_ids = set()

        for gene in gene2range :
            for i in gene2range[gene] :
                query_ids.add(i.name)

        query_seqs = {}

        f = FastqFile(fn)
        f.open()

        for seq in f :
            seq_id = self.__munge_contigname(seq.id)
            if seq_id in query_ids :
                query_seqs[seq_id] = seq.sequence

        f.close()

        return query_seqs

    def __extract_alignment_ranges(self, alignment_dir, min_identity) :
        gene2range = {}

        # find out what alignments overlap
        for fn in self.__get_dna_alignments(alignment_dir) :
            for query,reference in self.__get_queries(fn, allow_empty=True) :
                start,stop,identity = self.__get_alignment_columns(query, reference)
                
                if identity < min_identity :
                    continue

                if reference.id not in gene2range :
                    gene2range[reference.id] = []

                gene2range[reference.id].append(\
                        AlignmentRange(start, stop, self.__munge_paganorfname(query.id), ''.join([i for i in query.sequence if i != '-']))
                        )

        return gene2range

    def __output_scaffold(self, contig_filename, output_filename, ignore_these_contigs, extra_contigs) :
        scaf = open(output_filename, "w")

        for ec in extra_contigs :
            print >> scaf, str(ec)

        fq = FastqFile(contig_filename)
        fq.open()

        for seq in fq :
            seq.id = self.__munge_contigname(seq.id)
            if seq.id not in ignore_these_contigs :
                print >> scaf, str(seq)

        fq.close()
        scaf.close()

    def scaffold(self, contig_filename, alignment_dir, output_filename, min_identity) :
        # extract alignment info
        gene2range = self.__extract_alignment_ranges(alignment_dir, min_identity)
        num_contigs = sum([len(i) for i in gene2range.values()])
        self.info("found %d queries aligned to %d genes" % (num_contigs, len(gene2range)))

        # find where contig overlaps are unambiguous and output
        merged_contigs = set()
        merged_versions = []
        num_contigs = 0
        num_super_contigs = 0

        for gene in gene2range :
            # __merge can deal with this easily because it is based on reduce
            # however, I would prefer the stats to not include them
            if len(gene2range[gene]) == 1 :
                continue

            merge = self.__merge(gene2range[gene])
            
            if merge :
                merge.name += (" (aligned to %s)" % gene)

                merged_contigs.update([ar.name for ar in gene2range[gene]])
                merged_versions.append(merge)

                num_contigs += len(gene2range[gene])
                num_super_contigs += 1
            else :
                self.warn("conflict merging contigs aligned to %s" % gene)

        self.info("merged %d contigs into %d super-contigs" % (num_contigs, num_super_contigs))

        self.__output_scaffold(contig_filename, output_filename, merged_contigs, merged_versions)

