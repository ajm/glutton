import re

assembler_regex = {
        'oldtrinity'        : "^(comp\d+\_c\d+)\_seq\d+",                   # e.g. ">comp4_c0_seq1 len=253 path=[3991:0-252]"
        'trinity'           : "^(TRINITY_.*_c\d+_g\d+)_i\d+",              # e.g. ">TRINITY_DN1000_c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]"
        'oases'             : "^(Locus\_\d+)\_",                            # e.g. "Locus_1_Transcript_1/2_Confidence_0.333_Length_2022"
        'soapdenovotrans'   : "^scaffold\d+\ (Locus\_\d+)_\d+\ |^(C\d+)\ ", # e.g. "scaffold14 Locus_83_2 42.9 COMPLEX"
        'transabyss'        : "^(.+)$",                                     # e.g. "R420866 3059 43020273 416906+,...,419962-"      see (1)
        'none'              : "^(.+)$"  
    }

# (1) https://groups.google.com/forum/#!searchin/trans-abyss/isoform/trans-abyss/wKk1xXwySzA/FXu2pZulwfsJ

supported_assemblers = tuple(assembler_regex.keys())

class AssemblerOutput(object) :
    def __init__(self, assembler_name) :
        if assembler_name not in assembler_regex :
            raise KeyError("%s not found" % assembler_name)

        self.name = assembler_name
        self.regex = re.compile(assembler_regex[self.name])

    def match(self, s) :
        return self.regex.match(s)

    def __str__(self) :
        return "%s %s" % (self.__class__.__name__, self.name)

