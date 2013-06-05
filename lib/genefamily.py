
class GeneFamily(object) :
    def __init__(self) :
        self.genes = {}

    def add(self, name, sequence) :
        self.genes[name] = sequence

    def __getitem__(self, key) :
        return self.genes[key]

    def __setitem__(self, key, val) :
        self.genes[key] = val

    def _fasta(self, name) :
        return ">%s\n%s" % (name, self.genes[name])

    def __str__(self) :
        return "\n".join([self._fasta(name) for name in self.genes])

