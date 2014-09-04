from cogent.parse.fasta import MinimalFastaParser

class GeneFamily(object) :
    def __init__(self, fname=None) :
        self.genes = {}

        if fname :
            self._load(fname)

    def add(self, name, sequence) :
        self.genes[name] = sequence

    def iteritems(self) :
        return self.genes.iteritems()

    def __getitem__(self, key) :
        return self.genes[key]

    def _load(self, fname) :
        for label,seq in MinimalFastaParser(open(fname)) :
            self.genes[label] = seq

    def __setitem__(self, key, val) :
        self.genes[key] = val

    def _fasta(self, name) :
        return ">%s\n%s" % (name, self.genes[name])

    def __iter__(self) :
        return self.genes.__iter__()

    def __len__(self) :
        return len(self.genes)

    def __str__(self) :
        return "\n".join([self._fasta(name) for name in self.genes])

