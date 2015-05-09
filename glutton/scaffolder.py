from sys import exit, stderr, stdout
from glob import glob
from os.path import join
from collections import defaultdict, Counter

import re
import operator

from Bio import SeqIO

import pysam

from glutton.utils import get_log, check_dir
from glutton.info import GluttonInformation


DEBUG = True

class ScaffolderError(Exception) :
    pass

scaffold_counter = -1
scaffold_fmt = "glutton%d"

def scaffold_id() :
    global scaffold_fmt, scaffold_counter

    scaffold_counter += 1

    return scaffold_fmt % scaffold_counter

class Alignment(object) :
    def __init__(self, id, gene_name, start, end, seq, label, species, desc='singleton', contigs=None) :
        self.id = id
        self.gene_name = gene_name
        self.start = start
        self.end = end
        self.seq = seq
        self.desc = desc
        self.label = label
        self.species = species

        if contigs :
            self.contigs = contigs
        else :
            self.contigs = [ self.id ]

        if not self.id or not self.gene_name :
            return

        # XXX when being read from the file this should be verified
        #     so i don't need to check then fail
        m = re.match("^(comp\d+\_c\d+)\_seq\d+$", self.id)
        self.gene_id = m.group(1)

    def remove_chars(self, indices) :
        for i in sorted(indices, reverse=True) :
            self.seq = self.seq[:i] + self.seq[i+1:]

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

        raise ScaffoldError("impossible condition reached") 

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
                        self.gene_name, 
                        min(self.start, a2.start), 
                        max(self.end, a2.end), 
                        self.merge_contigs(a2), 
                        self.label,
                        self.species,
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

    def format_contig(self) :
        self.scaffold_id = scaffold_id()

        return ">%s contigs=%s gene=%s desc=%s\n%s" % (\
            self.scaffold_id, \
            ','.join(self.contigs), \
            self.gene_name, \
            self.desc, \
            self.seq.replace('-', '')\
          )

    def truncate_at_stop_codon(self) :
        self.seq = self.seq_stop_codon()

    def seq_stop_codon(self) :
        for i in range(0, len(self.seq), 3) :
            codon = self.seq[i : i+3]
            
            if codon in ('TAG', 'TAA', 'TGA') :
                s = self.seq[:i] + ("-" * (len(self.seq) - i))
                return s.replace('NNN', '---')

        return self.seq.replace('NNN', '---')

    def format_alignment(self, name) :
        return ">%s species=%s scaffolds=%s\n%s" % (name, self.id, ','.join(self.contigs), self.seq)

    def __len__(self) :
        return self.end - self.start

    def __str__(self) :
        return str((self.id, self.start, self.end, self.seq))

class Alignment2(Alignment) :
    def __init__(self, id, seq, contigs) :
        super(Alignment2, self).__init__(id, "", -1, -1, seq, "", "", contigs=contigs)

class Scaffolder(object) :
    def __init__(self, alignments_dir, db, contig_files, output_dir) :
        self.alignments_dir     = alignments_dir
        self.output_dir         = output_dir

        check_dir(self.output_dir, create=True)

        self.log = get_log()

        self.info = GluttonInformation(alignments_dir, db, contig_files)

        self.db = self.info.get_db()

        # perhaps slightly overambitious to exit, just stick to a warning      
        pending,failures = self.info.num_alignments_not_done()
        if pending != 0 :
            self.log.warn("%d alignments were not run!" % pending)

        # misleading because pagan refuses to align things which are 
        # distantly related anyway...
        #if failures != 0 :
        #    self.log.warn("%d alignments failed!" % failures)

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

    # read the alignment file and return a dictionary keyed on the gene name
    #   - there might be multiple orfs for a single query sequence, so keep a track of the best one
    def read_alignment(self, fname) :
        tmp = defaultdict(dict)
        genes = []

        # identity here is the proportion of non-gaps
        def calc_identity(seq1, seq2, start, end) :

            if start == end :
                return 0.0

            count = 0
            for i,j in zip(seq1[start:end], seq2[start:end]) :
                if (i,j) == ('-','-') :
                    continue
                if (i != '-') and (j != '-') :
                    count += 1

            return count / float(end - start)

        for s in SeqIO.parse(fname, 'fasta') :
            if not s.description.startswith('query') :
                gene_id     = s.description
                gene_name   = self.db.get_genename_from_geneid(gene_id)
                gene_seq    = str(s.seq)

                genes.append((gene_name, gene_seq))
                continue

            query_id        = self._orf_to_query_name(s.description)
            contig_id,label = self.info.get_contig_from_query(query_id)
            species         = self.info.label_to_species(label)

            seq = str(s.seq).replace('N', '-')
            
            contig_start = len(seq) - len(seq.lstrip('-'))
            contig_end   = contig_start + len(seq.strip('-'))

            # if we have seen this before
            if contig_id in tmp[gene_name] :
                new = calc_identity(gene_seq, seq, contig_start, contig_end)
                old = calc_identity(gene_seq, *tmp[gene_name][contig_id][:3])

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
            tmp[gene_name][contig_id] = (seq, contig_start, contig_end, label, species)

        # convert from a dict of dicts to a dict of lists
        tmp2 = defaultdict(list)

        #print fname
        for gene in tmp :
            #print "\t", gene
            for contig in tmp[gene] :
                #print "\t\t", contig
                seq,start,end,label,species = tmp[gene][contig]
                tmp2[gene].append(Alignment(contig, gene, start, end, seq, label, species))

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


            if DEBUG and (not no_print) and (len(group) > 1) :
                self.print_alignments(group)


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

    def consensus_for_msa(self, reference, alignment, bamfiles) :
        if len(alignment) == 1 :
            return Alignment2(alignment[0].species, alignment[0].seq, [alignment[0].scaffold_id])

        # this is buggy if within a species there are FAKE and real bam files
        coverage = []
        for a in alignment :
            numerator = 1
            denominator = len(a) if len(a) > 0 else 1

            if a.label in bamfiles :
                try :
                    numerator = bamfiles[a.label].count(str(a.id))
                except ValueError :
                    pass

            coverage.append(numerator / float(denominator))

        # loop through columns building up consensus sequence
        s = ""
        for column_index,column in enumerate(zip(*[ a.seq for a in alignment ])) :
            counts = Counter()

            for alignment_index,chcov in enumerate(zip(column, coverage)) :
                ch,cov = chcov

                if alignment[alignment_index].in_range(column_index) :
                    counts[ch] += cov

            # only gaps
            if not counts :
                s += '-'
                continue

            # add column from contig with highest coverage
            ch,cov = counts.most_common(1)[0]
            s += ch

        return Alignment2(alignment[0].species, s, [a.scaffold_id for a in alignment])

    def remove_common_gaps(self, alignment) :
        indices = []
        
        for index,chars in enumerate(zip(*[ a.seq for a in alignment ])) :
            if chars.count('-') == len(chars) :
                indices.append(index)

        for a in alignment :
            a.remove_chars(indices)

        return alignment

    def truncate_at_stop_codons(self, alignment) :
        for a in alignment :
            a.truncate_at_stop_codon()

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

                for a in self.merge_alignments(contigs[gene_name]) :
                    print >> output_files[a.label], a.format_contig()

                    if a.species not in merged_contigs[gene_name] :
                        merged_contigs[gene_name][a.species] = []

                    merged_contigs[gene_name][a.species].append(a)


            # merge sequences from the same species
            # find stop codon and truncate sequences
            # delete columns with only gaps
            # then write out to a file in self.output_dir
            # in MSA have >species_name contents=gluttonX,gluttonY,gluttonZ
            new_alignment = []

            for gene_name,gene_seq in genes :
                new_alignment.append(Alignment2(self.db.species, gene_seq, [gene_name]))
                
                for species in merged_contigs[gene_name] :
                    new_alignment.append(self.consensus_for_msa(gene_seq, merged_contigs[gene_name][species], bam_files))
                
            self.truncate_at_stop_codons(new_alignment)
            self.remove_common_gaps(new_alignment)

            counter += 1
            with open(join(self.output_dir, "msa%d.fasta" % counter), 'w') as f :
                for index,a in enumerate(new_alignment) :
                    print >> f, a.format_alignment("seq%d" % (index + 1))

            complete_files += 1
            stderr.write("\rINFO processed %d / %d alignments " % (complete_files, total_files))
            stderr.flush()

        stderr.write("\rINFO processed %d / %d alignments \n" % (complete_files, total_files))
        stderr.flush()

        self.log.info("created %d multiple sequence alignments" % (counter + 1))

        return aligned_contigs

    def scaffold(self) :
        self.log.info("starting scaffolding")

        # open scaffold + BAM files
        output_files = {}
        bam_files = {}

        for label in self.info.get_sample_ids() :
            fname = join(self.output_dir, label + '.fasta')
            self.log.info("creating %s ..." % fname)
            output_files[label] = open(fname, 'w')

            bamfilename = self.info.label_to_bam(label)

            if bamfilename :
                if bamfilename == 'FAKE' :
                    self.log.warn("bam file missing for sample %s" % label)
                    continue

                bam_files[label] = pysam.AlignmentFile(bamfilename)

        # process alignment files
        aligned_contigs = self.process_alignments(output_files, bam_files)
        
        # append contigs remaining contigs
        self.output_unscaffolded_contigs(output_files, aligned_contigs)

        # close scaffold + BAM files
        for f in output_files.values() + bam_files.values() :
            f.close()

        return 0

    def fasta_output(self, contig, seq, desc) :
        return ">%s contigs=%s desc=%s\n%s" % (scaffold_id(), contig, desc, seq)

    def output_unscaffolded_contigs(self, output_files, aligned_contigs) :

        for label in self.info.get_sample_ids() :
            fout = output_files[label]

            for r in SeqIO.parse(self.info.label_to_contigs(label), 'fasta') :

                if r.id in aligned_contigs[label] :
                    continue

                # unused due to being filtered out
                if not self.info.contig_used(r.id, label) :
                    print >> fout, self.fasta_output(r.id, r.seq, 'filtered')

                # unassigned by blast
                elif not self.info.contig_assigned(r.id, label) :
                    print >> fout, self.fasta_output(r.id, r.seq, 'assignment_failed')

                # unaligned by pagan
                else :
                    print >> fout, self.fasta_output(r.id, r.seq, 'alignment_failed')

