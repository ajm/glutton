from Bio import SeqIO

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
        if fmt != 'fasta' :
            raise NotImplementedError()

        if full :
            return ">%s name=\"%s\"\n%s\n" % (self.id, self.name, self.sequence)

        return ">%s\n%s\n" % (self.id, self.sequence)

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

    for famid in families :
        tmp[famid] = GeneFamily([ Gene(*families[famid][geneid], id=geneid) for geneid in families[famid] ], id=famid)

    return tmp

