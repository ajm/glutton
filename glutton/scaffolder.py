from sys import exit, stderr, stdout
from glob import glob
from os import abort
from os.path import join, getsize
from collections import defaultdict

import re
import operator

from Bio import SeqIO

import pysam

from glutton.utils import get_log, check_dir
from glutton.info import GluttonInformation, GluttonParameters
from glutton.db import GluttonDB
from glutton.assembler_output import AssemblerOutput


#from pycallgraph import PyCallGraph
#from pycallgraph.output import GraphvizOutput



DEBUG = False

class ScaffolderError(Exception) :
    pass

scaffold_counter = -1
scaffold_fmt = "glutton%d"

start_codon = 'ATG'
stop_codons = ('TAG', 'TAA', 'TGA')
protein2codons = { 
        'A' : ('GCT', 'GCC', 'GCA', 'GCG'),
        'R' : ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'N' : ('AAT', 'AAC'),
        'D' : ('GAT', 'GAC'),
        'C' : ('TGT', 'TGC'),
        'Q' : ('CAA', 'CAG'),
        'E' : ('GAA', 'GAG'),
        'G' : ('GGT', 'GGC', 'GGA', 'GGG'),
        'H' : ('CAT', 'CAC'),
        'I' : ('ATT', 'ATC', 'ATA'),
        'L' : ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
        'K' : ('AAA', 'AAG'),
        'M' : ('ATG',),
        'F' : ('TTT', 'TTC'),
        'P' : ('CCT', 'CCC', 'CCA', 'CCG'),
        'S' : ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'T' : ('ACT', 'ACC', 'ACA', 'ACG'),
        'W' : ('TGG',),
        'Y' : ('TAT', 'TAC'),
        'V' : ('GTT', 'GTC', 'GTA', 'GTG'),
        '-' : ('---',),
        'X' : ('NNN',),
        '*' : ('TAG', 'TAA', 'TGA')
    }

codon2protein = {}

for k in protein2codons :
    for v in protein2codons[k] :
        codon2protein[v] = k

# biopython does not support gapped translation
def translate(s) :
    return ''.join([ codon2protein[s[i:i+3]] for i in range(0, len(s), 3) ])

def sequence_limits(s) :
    start = len(s) - len(s.lstrip('-'))
    end = start + len(s.strip('-'))
    return start,end

def scaffold_id() :
    global scaffold_fmt, scaffold_counter

    scaffold_counter += 1

    return scaffold_fmt % scaffold_counter

class Alignment(object) :
    def __init__(self, id, gene_id, gene_name, start, end, seq, label, species, desc='singleton', contigs=None) :
        self.id = id                # contig id (assembler specific)
        self.gene_id = gene_id      # gene id (assembler specific)
        self.gene_name = gene_name  # gene contig is aligned to
        self.start = start          # start pos
        self.end = end              # end pos
        self.seq = seq              # alignment sequence
        self.desc = desc            # scaffold description, e.g. singleton, merged, etc
        self.label = label          # sample label (from cli)
        self.species = species      # species label (from cli)

        if contigs :
            self.contigs = contigs
        else :
            self.contigs = [ self.id ]

        self.start_length = len(self)

        if not len(self) :
            raise ScaffolderError("alignment length zero")

    @property
    def contig_id(self) :
        return str(self.id).split()[0]

    def remove_chars(self, indices) :
        tmp = ""

        for i in range(1, len(indices)) :
            tmp += self.seq[indices[i-1]+1:indices[i]]

        self.seq = tmp
        self.start,self.end = sequence_limits(self.seq)

        if not len(self) :
            raise ScaffolderError("alignment length zero")

        #for i in sorted(indices, reverse=True) :
        #    self.seq = self.seq[:i] + self.seq[i+1:]

    def in_range(self, i) :
        return self.start <= i < self.end

    def __getitem__(self, i) :
        return self.seq[i] if self.start <= i < self.end else None

    def _ensure_order(self, a2, func) :
        return func(self, a2) if self.start <= a2.start else func(a2, self)

    def extract(self, start, end) :
        return self.seq[start:end]

    # assumes a starts before b
    def _overlaps(self, a, b) :
        if a.end > b.end :
            return (b.start, b.end)

        elif a.end <= b.start :
            return None

        elif a.end > b.start :
            return (b.start, a.end)

        elif a.end == b.start :
            return (b.start, a.end)

        raise ScaffolderError("impossible condition reached")

    def overlaps(self, a2) :
        return self._ensure_order(a2, self._overlaps)

    def from_same_file(self, a2) :
        return a2.label == self.label

    def _remove_common_gaps(self, seq_a, seq_b) :
        new_a = ""
        new_b = ""

        for a,b in zip(seq_a, seq_b) :
            if (a,b) == ('-','-') :
                continue

            new_a += a
            new_b += b

        return new_a, new_b

    def _mergeable(self, a, b) :
        overlap = self._overlaps(a, b) 

        if not overlap :
            return False

        slice_a, slice_b = self._remove_common_gaps(a.extract(*overlap), b.extract(*overlap))

        for a,b in zip(slice_a, slice_b) :
            if a == '-' or b == '-' :
                return False

        return True

    def mergeable(self, a2) :
        return self._ensure_order(a2, self._mergeable)

    def isoforms(self, a2) :
        return self.gene_id == a2.gene_id

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
                        self.gene_id,
                        self.gene_name, 
                        min(self.start, a2.start), 
                        max(self.end, a2.end), 
                        self.merge_contigs(a2), 
                        self.label,
                        self.species,
                        'merged',
                        self.contigs + a2.contigs)

    def __iadd__(self, a2) :
        raise ScaffolderError("unsupported")
        self.start = min(self.start, a2.start)
        self.end = max(self.end, a2.end)
        self.seq = self.merge_contigs(a2)
        self.contigs += a2.contigs
        self.desc = 'merged'
    
    def get_desc(self) :
        return self.desc

    def format_contig(self) :
        self.scaffold_id = scaffold_id()
        contig_value = ','.join([ i.split()[0] for i in self.contigs ])

        return ">%s contigs=%s gene=%s desc=%s\n%s" % (\
            self.scaffold_id, \
            contig_value, \
            self.gene_name, \
            self.desc, \
            self.seq.replace('-', '')\
          )

    def trim_at_ATG(self, pos) :
        trim_pos = pos

        for ind in range(0, pos+3, 3)[::-1] :
            codon = self.seq[ind : ind+3]

            if codon == start_codon :
                trim_pos = ind
                break

        self.seq = ('-' * trim_pos) + self.seq[trim_pos:]
        self.start,self.end = sequence_limits(self.seq)

        if not len(self) :
            raise ScaffolderError("alignment length zero")

    def truncate_at_stop_codon(self) :
        self.seq   = self.seq_stop_codon()
        self.start,self.end = sequence_limits(self.seq)
        
        if not len(self) :
            raise ScaffolderError("alignment length zero")

    def seq_stop_codon(self) :
        for i in range(0, len(self.seq), 3) :
            codon = self.seq[i : i+3]
            
            if codon in stop_codons :
                s = self.seq[:i] + ("-" * (len(self.seq) - i))
                return self.munge(s)

        return self.munge(self.seq)

    def munge(self, s) :
        head,sep,tail = s.rpartition('NNN')
        
        if tail.count('-') == len(tail) :
            return head + '---' + tail

        return s

    def format_alignment(self, name) :
        return ">%s species=%s scaffolds=%s id=%.3f cov=%.3f len=%d\n%s" % (name, self.id, ','.join(self.contigs), self.prot_id, self.coverage, len(self), self.seq)

    def non_gap_count(self) :
        return len(self.seq) - self.seq.count('-')

    def __len__(self) :
        return self.end - self.start

    def __str__(self) :
        return str((self.id, self.start, self.end, self.start_length, self.seq))

class Alignment2(Alignment) :
    def __init__(self, id, gene_name, seq, contigs) :
        start,end = sequence_limits(seq)
        super(Alignment2, self).__init__(id, "", gene_name, start, end, seq, "", "", contigs=contigs)

class Scaffolder(object) :
    def __init__(self, top_level_directory, reference_fname, assembler_name, protein_identity, alignment_length, min_gene_coverage, do_not_trim=False, testmode='none') :
        self.alignments_dir     = join(top_level_directory, 'alignments')
        self.output_dir         = join(top_level_directory, 'postprocessing')
        self.protein_identity   = protein_identity
        self.alignment_length   = alignment_length
        self.min_gene_coverage  = min_gene_coverage
        self.trim               = not do_not_trim
        self.testmode           = testmode

        self.scaffold_dir       = join(self.output_dir, 'scaffolds')
        self.genefamily_msa_dir = join(self.output_dir, 'genefamily_msa')
        self.gene_msa_dir       = join(self.output_dir, 'gene_msa')

        check_dir(self.output_dir, create=True)
        check_dir(self.scaffold_dir, create=True)
        check_dir(self.genefamily_msa_dir, create=True)
        check_dir(self.gene_msa_dir, create=True)

        self.log = get_log()

        self.param = GluttonParameters(top_level_directory)
        self.db = GluttonDB(reference_fname)
        self.info = GluttonInformation(self.alignments_dir, self.param, self.db)

        # check reference was the same
        if not self.param.same_reference(self.db) :
            self.log.warn("current reference is %s (checksum=%s), alignments were performed using %s (checksum=%s)" % (reference_fname, self.db.checksum, self.db.filename, self.param['db_checksum']))
            #exit(1);

        # perhaps slightly overambitious to exit, just stick to a warning      
        pending,failures = self.info.num_alignments_not_done()
        if pending != 0 :
            self.log.warn("%d alignments were not run!" % pending)

        self.assembler = AssemblerOutput(assembler_name)

        # e.g. query39806_orf1
        self.orfname_regex = re.compile("^(query\d+)\_orf(\d)$")

    def stop(self) :
        pass

    def _orf_to_query_name(self, name) :
        m = self.orfname_regex.match(name)

        if not m :
            raise ScaffolderError("unexpected query name (%s)" % name)

        return m.group(1)

    def _assembler_gene_name(self, name) :
        m = self.assembler.match(name)

        if not m :
            raise ScaffolderError("unexpected gene name (%s)" % name)

        if m.group(1) :
            return m.group(1)
        else :
            return m.group(2)

    # read the alignment file and return a dictionary keyed on the gene name
    #   - there might be multiple orfs for a single query sequence, so keep a track of the best one
    #@profile
    def read_alignment(self, fname) :
        tmp = defaultdict(dict)
        genes = []

        for s in SeqIO.parse(fname, 'fasta') :
            if not s.description.startswith('query') :
                gene_id     = s.description
                gene_name   = self.db.get_genename_from_geneid(gene_id)
                gene_seq    = str(s.seq)
                gene_prot   = translate(gene_seq)
                gene_start,gene_end = sequence_limits(gene_prot)
                
                genes.append((gene_name, gene_seq))
                continue

            query_id        = self._orf_to_query_name(s.description)
            contig_id,label = self.info.get_contig_from_query(query_id)
            species         = self.param.get_species(label)
            assembler_geneid = self._assembler_gene_name(contig_id)

            seq = str(s.seq).replace('N', '-')
            
            contig_start,contig_end = sequence_limits(seq)

            # pagan bug?
            if (contig_end - contig_start) == 0 :
                continue

            overlap_start = max(contig_start / 3, gene_start)
            overlap_end   = min(contig_end   / 3, gene_end)

            # require all alignments to firmly overlap gene
            if ((overlap_end - overlap_start) * 3) < 100 :
                continue

            ref_identity = self.protein_similarity(gene_prot, 
                                                   translate(seq), 
                                                   overlap_start,
                                                   overlap_end)

            # user defined identity in protein space
            if ref_identity < self.protein_identity :
                continue

            # if we have seen this before
            if (contig_id in tmp[gene_name]) and (ref_identity < tmp[gene_name][contig_id][-1]) : 
                continue

            # first three need to be seq, contig_start, contig_end
            tmp[gene_name][contig_id] = (seq, contig_start, contig_end, label, species, assembler_geneid, s.description, ref_identity)


        # convert from a dict of dicts to a dict of lists
        tmp2 = defaultdict(list)

        self.log.debug("read %s" % fname)
        for gene in tmp :
            self.log.debug("\tgene = %s" % gene)

            for contig in tmp[gene] :
                seq,start,end,label,species,assembler_geneid,queryid,ref_identity = tmp[gene][contig]
                self.log.debug("\t\tquery id = %s (%d,%d)" % (queryid,start,end))
                
                try :
                    tmp2[gene].append(Alignment(contig, assembler_geneid, gene, start, end, seq, label, species))
                
                except ScaffolderError :
                    self.log.debug("empty sequence %s in %s" % (queryid, fname))
                    continue

        return tmp2, genes

    def group_alignments(self, alignments) :
        groups = []

        # for each alignment
        for align in alignments :
            # for each group
            add_to_group = False
            for group in groups :
                # if 'align' overlaps with a member of that group
                # are they were from the same input file (i.e. label)
                for member in group :
                    if member.overlaps(align) and member.from_same_file(align) :
                        add_to_group = True
                        break

                # add to the group
                if add_to_group :
                    group.append(align)
                    break

            if not add_to_group :
                groups.append([ align ])

        return groups

    def group_alignments_by_file(self, alignments) :
        groups = defaultdict(list)

        for a in alignments :
            groups[a.label].append(a)

        return groups

    def group_cannot_be_merged_isoforms(self, group) :
        return self.group_cannot_be_merged(group, consider_isoforms=True)

    def group_cannot_be_merged(self, group, consider_isoforms=False) :
        for i in range(0, len(group)) :
            a = group[i]
            for j in range(i+1, len(group)) :
                b = group[j]
                
                if a.overlaps(b) :
                    if consider_isoforms :
                        if a.isoforms(b) :
                            return True
                    else :
                        if not a.mergeable(b) :
                            return True
        
        return False

    def merge_alignments(self, alignments) :
        global DEBUG

        #if self.testmode != 'none' :
        #    return alignments

        # perform merges of overlaps first
        unmerged_labels = set()
        merged_groups = []

        def label_all(group, label) :
            for i in group :
                i.desc = label

        for group in self.group_alignments(alignments) :
            no_print = False

            # if the assembler thinks these two contigs are isoforms of the
            # same gene then we should not attempt to merge them
            if self.group_cannot_be_merged_isoforms(group) :
                unmerged_labels.add(group[0].label)

                if len(set([ g.gene_id for g in group ])) == 1 :
                    label_all(group, 'single_gene_multiple_isoform')
                else :
                    label_all(group, 'multiple_gene_multiple_isoform')

                merged_groups += group

                no_print = True

            # some contigs mapping to the same region of the gene cannot be
            # merged because the differences between them is too large
            # (defined in Alignment class)
            elif self.group_cannot_be_merged(group) :
                unmerged_labels.add(group[0].label)
                label_all(group, 'conflict')
                merged_groups += group

            # these can be merged trivially
            else :
                merged_groups.append( reduce(operator.add, sorted(group, key=lambda x : x.start)) )


            #if DEBUG and (not no_print) and (len(group) > 1) :
            #    self.print_alignments(group)


        # if there were any conflicts for a given gene in an alignment, then 
        # just output what has been merged, if there were no conflicts for 
        # a gene then concatenate the islands of alignments with N's
        tmp = []
        grouped_by_file = self.group_alignments_by_file(merged_groups)
        for label in grouped_by_file :
            if label in unmerged_labels :
                tmp += grouped_by_file[label]
            else :
                tmp.append( reduce(operator.add, sorted(grouped_by_file[label], key=lambda x : x.start)) )


        return tmp

    def print_alignments(self, alignments) :
        
        f = open('glutton_debug_scaffolds.txt', 'a')

        for a in alignments :
            print >> f, "%s\t%s" % (a.id, a.seq)

        alignment_diff = ""
    
        for column in range(len(alignments[0].seq)) :
            chars = set([ a[column] for a in alignments if a[column] ])
            if ('-' in chars) and (len(chars) > 1) :
                alignment_diff += 'X'
            elif len(chars) < 2 :
                alignment_diff += ' '
            else :
                alignment_diff += 'M'

        print >> f, "difference       \t%s" % alignment_diff
        print >> f, ""

        f.close()

    # XXX now in Alignment class
    def trim_at_ATG(self, seq, pos) :
        trim_pos = pos

        for ind in range(0, pos+3, 3)[::-1] :
            codon = seq[ind : ind+3]

            if codon == start_codon :
                trim_pos = ind
                break

            #if codon in stop_codons :
            #    break

        return ('-' * trim_pos) + seq[trim_pos:]

    def consensus_for_msa(self, reference, alignments, bamfiles) :

        if len(alignments) == 1 :
            a = alignments[0]
            a.trim_at_ATG(reference.start)
            if self.trim :
                a.truncate_at_stop_codon()
            return Alignment2(a.species, a.gene_name, a.seq, [a.contig_id])

        if self.testmode == 'none' :
            return self.consensus_for_msa_glutton(reference, alignments, bamfiles)

        else :
            lengths = [ len(a.seq.replace('-','')) for a in alignments ]
            identities = []
            coverages = []

            gene_prot = translate(reference.seq)

            for a in alignments :
                # identity
                overlap_start = max(a.start, reference.start)
                overlap_end   = min(a.end  , reference.end  )

                ident = self.protein_similarity(gene_prot, translate(a.seq), overlap_start/3, overlap_end/3)
                identities.append(ident)

                # coverage (depth)
                if a.label in bamfiles :
                    try :
                        id = str(a.id).split()[0]
                        coverages.append(bamfiles[a.label].count(id))
                    except :
                        coverages.append(1)

            coverages = [ c / float(l) for c,l in zip(coverages, lengths) ]

            if self.testmode == 'length' :
                top_hit = sorted(zip(lengths, identities, coverages, range(len(lengths))))[-1][-1]
            elif self.testmode == 'identity' :
                top_hit = sorted(zip(identities, lengths, coverages, range(len(lengths))))[-1][-1]
            else :
                top_hit = sorted(zip(coverages, identities, lengths, range(len(lengths))))[-1][-1]

            a = alignments[top_hit]
            a.trim_at_ATG(reference.start)
            if self.trim :
                a.truncate_at_stop_codon()
            #seq = self.trim_at_ATG(a.seq, reference.start)
            return Alignment2(a.species, a.gene_name, a.seq, [a.contig_id])

    #@profile
    def consensus_for_msa_glutton(self, reference, alignments, bamfiles) :

        if len(alignments) == 1 :
            a = alignments[0]
            a.trim_at_ATG(reference.start)
            if self.trim :
                a.truncate_at_stop_codon()
            return Alignment2(a.species, a.gene_name, a.seq, [a.contig_id])

        # this is buggy if within a species there are FAKE and real bam files
        coverage = []

        for a in alignments :
            numerator = 1
            denominator = len(a.seq.replace('-', ''))

            if a.label in bamfiles :
                try :
                    # BWA only uses the fasta id, but we need to store the complete
                    # description line as an id because soapdenovotrans does not provide
                    # the locus information in the first token, but the second
                    #id = str(a.id).split()[0] # moved to a property in Alignment

                    numerator = bamfiles[a.label].count(a.contig_id)

                except ValueError :
                    pass

            coverage.append(numerator / float(denominator))


        s = "-" * len(alignments[0].seq)
        
        for cov,a in sorted(zip(coverage, alignments)) :
            #tmp = ""
            #for c1,c2 in zip(a.seq[a.start:a.end], s[a.start:a.end]) :
            #    tmp += (c1 if c1 not in ('-','N') else c2)

            a.trim_at_ATG(reference.start)
            if self.trim :
                a.truncate_at_stop_codon()
            subseq = a.seq[a.start:a.end]

            if 'N' not in subseq :
                s = s[:a.start] + subseq + s[a.end:]

            else :
                tmp = ""

                for ind,c in enumerate(subseq) :
                    tmp += (c if c != 'N' else s[a.start + ind])

                s = s[:a.start] + tmp + s[a.end:]

        #s = self.trim_at_ATG(s, reference.start)
        return Alignment2(alignments[0].species, alignments[0].gene_name, s, [ a.contig_id for a in alignments ])

    def remove_common_gaps(self, alignment) :
        indices = [-1, len(alignment[0].seq)]
        num_rows = len(alignment)
        
        for index,chars in enumerate(zip(*[ a.seq for a in alignment ])) :
            if (chars.count('-') + chars.count('N')) == num_rows :
            #if chars.count('-') == num_rows :
                indices.append(index)

        indices.sort()

        for a in alignment :
            a.remove_chars(indices)

        return alignment

    def nucleotide_overlap(self, ref, query, start, end) :
        if end <= start :
            return 0

        r = ref[start:end]
        q = query[start:end]

        covered = 0

        for cq,cr in zip(q,r) :
            if (cq,cr) == ('-','-') :
                continue

            if cq == 'N' :
                continue

            covered += 1

        return covered

    def gene_coverage(self, ref, query) :
        total = 0
        covered = 0

        for i in range(ref.start, ref.end) :
            cq = query.seq[i]
            cr = ref.seq[i]

            if (cq,cr) == ('-','-') :
                continue

            if (query.start <= i < query.end) and (cq != 'N') :
                covered += 1

            total += 1

        return covered / float(total)

    def protein_similarity(self, ref, query, start, end) :
        if end <= start :
            return 0.0

        r = ref[start:end]
        q = query[start:end]

        identical = 0
        length = 0

        for cq,cr in zip(q,r) :
            if (cq,cr) == ('-','-') :
                continue

            if cq == 'X' :
                continue

            if cq == cr :
                identical += 1

            length += 1


        return identical / float(length)

    def process_alignments(self, output_files, bam_files) :
        global DEBUG

        counter = -1
        aligned_contigs = defaultdict(set)

        alignment_files = glob(join(self.alignments_dir, 'glutton*.nucleotide'))

        complete_files = 0
        total_files = len(alignment_files)

        stderr.write("\rINFO processed %d / %d alignments " % (complete_files, total_files))
        stderr.flush()

        for fname in alignment_files :
            contigs, genes = self.read_alignment(fname)
            merged_contigs = defaultdict(dict)

            # for each gene, merge the contigs from the same input file
            # and write to output
            for gene_name in contigs :
                for a in contigs[gene_name] :
                    aligned_contigs[a.label].add(a.id)

                    if a.species not in merged_contigs[gene_name] :
                        merged_contigs[gene_name][a.species] = []

                    merged_contigs[gene_name][a.species].append(a)

                for a in self.merge_alignments(contigs[gene_name]) :
                    print >> output_files[a.label], a.format_contig()

            # merge sequences from the same species
            # find stop codon and truncate sequences
            # delete columns with only gaps
            # then write out to a file in self.output_dir
            # in MSA have >species_name contents=gluttonX,gluttonY,gluttonZ
            new_alignment = []
            non_reference_seq = 0

            for gene_name,gene_seq in genes :
                ref = Alignment2(self.db.species, gene_name, gene_seq, [gene_name])
                new_alignment.append(ref)

                #gene_prot = translate(ref.seq)
                ref.prot_id = 1.0
                ref.coverage = 1.0

                for species in merged_contigs[gene_name] :
                    try :
                        tmp = self.consensus_for_msa(ref, merged_contigs[gene_name][species], bam_files)
                        #tmp.truncate_at_stop_codon() # this only needs to be here for the testmodes, otherwise it is redundant

                    except ScaffolderError, se :
                        continue

                    # check length vs alignment_length
                    #if len(tmp) < self.alignment_length :
                    #    continue

                    overlap_start = max(tmp.start, ref.start)
                    overlap_end   = min(tmp.end  , ref.end  )

                    overlap_bases = self.nucleotide_overlap(ref.seq, tmp.seq, overlap_start, overlap_end)

                    if overlap_bases < self.alignment_length :
                        continue

                    coverage = self.gene_coverage(ref, tmp)

                    if coverage < self.min_gene_coverage :
                        continue

                    prot_identity = self.protein_similarity(translate(ref.seq),
                                                            translate(tmp.seq),
                                                            overlap_start / 3,
                                                            overlap_end / 3)

                    if prot_identity < self.protein_identity :
                        continue

                    tmp.prot_id = prot_identity
                    tmp.coverage = coverage

                    new_alignment.append(tmp)
                    non_reference_seq += 1


            if non_reference_seq != 0 :
                counter += 1

                self.write_alignment(join(self.genefamily_msa_dir, "msa%d.fasta" % counter), new_alignment)

                subalignments = defaultdict(list)

                for a in new_alignment :
                    subalignments[a.gene_name].append(a)

                for k,v in subalignments.iteritems() :
                    if len(v) > 1 :
                        self.write_alignment(join(self.gene_msa_dir, "%s.fasta" % k), v)


            complete_files += 1
            stderr.write("\rINFO processed %d / %d alignments " % (complete_files, total_files))
            stderr.flush()

        stderr.write("\rINFO processed %d / %d alignments \n" % (complete_files, total_files))
        stderr.flush()

        self.log.info("created %d multiple sequence alignments" % (counter + 1))

        return aligned_contigs

    def write_alignment(self, fname, alignment) :
        self.remove_common_gaps(alignment)
        with open(fname, 'w') as f :
            for index,a in enumerate(alignment) :
                print >> f, a.format_alignment("seq%d" % (index + 1))

    def scaffold(self) :
        self.log.info("starting scaffolding")

        # open scaffold + BAM files
        output_files = {}
        bam_files = {}

        for label in self.param.get_sample_ids() :
            fname = join(self.scaffold_dir, label + '.fasta')
            self.log.info("creating %s ..." % fname)
            output_files[label] = open(fname, 'w')

            bamfilename = self.param.get_bam(label)

            if bamfilename :
                if bamfilename == 'FAKE' :
                    self.log.warn("bam file missing for sample %s" % label)
                    continue

                bam_files[label] = pysam.AlignmentFile(bamfilename)

#        graphviz = GraphvizOutput()
#        graphviz.output_file = 'scaffolder.png'
#
#        with PyCallGraph(output=graphviz):

        # process alignment files
        aligned_contigs = self.process_alignments(output_files, bam_files)
        
        # append contigs remaining contigs
        self.output_unscaffolded_contigs(output_files, aligned_contigs)

        # close scaffold + BAM files
        for f in output_files.values() + bam_files.values() :
            f.close()

        return 0

    def fasta_output(self, contig, seq, desc) :
        return ">%s contigs=%s desc=%s\n%s" % (scaffold_id(), contig.split()[0], desc, seq)

    def output_unscaffolded_contigs(self, output_files, aligned_contigs) :

        for label in self.param.get_sample_ids() :
            fout = output_files[label]

            for r in SeqIO.parse(self.param.get_contigs(label), 'fasta') :

                contig_name = r.description

                if contig_name in aligned_contigs[label] :
                    continue

                # unused due to being filtered out
                if not self.info.contig_used(contig_name, label) :
                    print >> fout, self.fasta_output(contig_name, r.seq, 'filtered')

                # unassigned by blast
                elif not self.info.contig_assigned(contig_name, label) :
                    print >> fout, self.fasta_output(contig_name, r.seq, 'assignment_failed')

                # unaligned by pagan
                else :
                    print >> fout, self.fasta_output(contig_name, r.seq, 'alignment_failed')

