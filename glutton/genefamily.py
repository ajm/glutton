
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


# the ensembl downloader should be standalone and not know about
# these objects, this function converts what it returns
#   i.e.: [(name, seq),(name, seq), ... ]
# to gene and genefamily objects
def ensembl_to_glutton_internal(families) :
    tmp = {}
    
    for fam in families :
        gf = GeneFamily([ Gene(g[0], g[1]) for g in fam ])
        tmp[gf.id] = gf
    
    return tmp

