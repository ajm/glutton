import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from glutton.utils import get_log

# maybe this should extend Sequence from biopython?
class Gene(object) :

    id_counter = 0

    def __init__(self, name, sequence=None, id=None) :
        self.name = name
        self.sequence = sequence

        if id :
            self.id = id
        else :
            self.id = "gene%d" % Gene.id_counter
            Gene.id_counter += 1

    def __getitem__(self, i) :
        return self.sequence[i]

    def __getslice__(self, i, j) :
        return self.sequence[i:j]

    @property
    def seq(self) :
        return self.sequence

    def format(self, fmt, full=False) :
        if fmt not in ('fasta', 'protein', 'nucleotide') :
            raise NotImplementedError()

        if fmt == 'protein' :
            padding = ''
            r = len(self.sequence) % 3
            if r == 1 :
                padding = 'NN'
            elif r == 2 :
                padding = 'N'

            seq_to_print = Seq(self.sequence + padding, generic_dna).translate() 
        else :
            seq_to_print = self.sequence

        if full :
            return ">%s name=\"%s\"\n%s\n" % (self.id, self.name, seq_to_print)

        return ">%s\n%s\n" % (self.id, seq_to_print)

    # sometimes there are 'N's, remove them so long as they do not cause
    # a frame shift, if they do, replace the remainder with 'G' or something
#    def fix_seq(self, seq) :
#        if 'N' not in seq :
#            return seq
#
#        tmp = str(seq)
#        tmp2 = ""
#        for idx in range(0, len(tmp), 3) :
#            codon = tmp[idx:idx+3]
#
#            if codon == 'NNN' :
#                continue
#
#            tmp2 += codon.replace('N','G')
#
#        return Seq(tmp2)

    def open_reading_frames(self, strand=False) :
        tmp = []
        newid = self.id + "_orf%d"

        for i in range(3) :
            tmp_seq = Gene(self.name, self.sequence[i:], newid % (i + 1))

            # TODO make cli options
            if tmp_seq.max_length_orf() > 100 :
                tmp.append(tmp_seq)

        if strand :
            return tmp

        # ORFs on the other strand
        self.reverse_complement()

        for i in range(3) :
            tmp_seq = Gene(self.name, self.sequence[i:], newid % (i + 4))
            
            # TODO make cli options
            if tmp_seq.max_length_orf() > 100 :
                tmp.append(tmp_seq)

        self.reverse_complement()

        return tmp

    def max_length_orf(self) :
        stop_codons = ('TAA','TAG','TGA')
        codons = [ self.sequence[i:i+3] for i in range(0, len(self.sequence), 3) ]
        indices = [0]
        max_length = 0

        for index,codon in enumerate(codons) :
            if codon in stop_codons :
                indices.append(index)
        
        indices.append(len(codons))

        for i in range(len(indices) - 1) :
            length = indices[i+1] - indices[i]

            if length > max_length :
                max_length = length

        return max_length * 3

    def reverse_complement(self) :
        self.sequence = self.sequence[::-1]

        d = {
            'A':'T',
            'T':'A',
            'G':'C',
            'C':'G',
            'N':'N'
        }

        self.sequence = ''.join([ d[i] for i in self.sequence ])

    def __len__(self) :
        return len(self.sequence)

    def __str__(self) :
        return "%s %s" % (self.id, self.name)

class GeneFamily(list) :

    id_counter = 0

    def __init__(self, genes=[], id=None) :
        list.__init__(self, genes)

        if id :
            self.id = id
        else :
            self.id = "genefamily%d" % GeneFamily.id_counter
            GeneFamily.id_counter += 1

        self.tree = None

    def set_tree(self, tree) :
        self.tree = tree

    def get_tree(self) :
        return self.tree

def seqlen(s) :
    return len(s) if not isinstance(s, list) else sum([ len(i) for i in s ])

def biopy_to_gene(s, id=None) :
    return Gene(s.id, s.seq, id)

def read_alignment_as_genefamily(f, name) :
    tmp = GeneFamily(id=name)
    
    # so long as this allows '-', I don't need AlignIO, I just want the contents
    for s in SeqIO.parse(f, 'fasta') :
        name = ' '.join(s.description.strip().split()[1:])
        tmp.append(Gene(name, s.seq, s.id))
    
    return tmp

# the ensembl downloader should be standalone and not know about
# these objects, this function converts what it returns
#   i.e.: [(name, seq),(name, seq), ... ]
# to gene and genefamily objects
def ensembl_to_glutton(families) :
    tmp = {}
    
    for fam in families :
        gf = GeneFamily([ Gene(g[0], g[1]) for g in fam ])
        tmp[gf.id] = gf
    
    return tmp

#{
#   fam1 :  {
#           gene1 : (name, seq)
#           },
#   fam2 :  {
#           gene2 : (name, seq),
#           gene3 : (name, seq)
#           }
#}
def glutton_to_json(families) :
    tmp = {}

    for famid in families :
        tmp[famid] = {}

        for gene in families[famid] :
            tmp[famid][gene.id] = (gene.name, gene.seq)

    return tmp

def json_to_glutton(families) :
    tmp = {}
    bad_family_count = 0
    bad_gene_count = 0

    for famid in families :
        #tmp[famid] = GeneFamily([ Gene(*families[famid][geneid], id=geneid) for geneid in families[famid] ], id=famid)
        
        fam = []
        bad = False
        for geneid in families[famid] :
            genename,geneseq = families[famid][geneid]

            # if i find this, keep on adding then i can print out the whole family later
            if geneseq == 'Sequenceunavailable' :
                bad = True

            fam.append( Gene(genename, geneseq, id=geneid) )
        
        if not bad :
            tmp[famid] = GeneFamily(fam, id=famid)
        else :
            bad_family_count += 1
            bad_gene_count += len(fam)

    if bad_family_count > 0 :
        get_log().error("%d bad gene families (containing %d genes)" % (bad_family_count, bad_gene_count))

    return tmp

if __name__ == '__main__' :
    from Bio import SeqIO

    s1 = SeqIO.read('test1.fa', 'fasta')
    s2 = SeqIO.read('test2.fa', 'fasta')

    g1 = Gene(s1.id, s1.seq)
    for i in g1.open_reading_frames() :
        print str(i)

    g2 = Gene(s2.id, s2.seq)
    for i in g2.open_reading_frames() :
        print str(i)

