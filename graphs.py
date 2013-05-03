import sys
import glob
import re

from lib.filetypes import FastqFile
from lib.datatypes import Sequence
from lib.tools import Exonerate, ExonerateHit, ExonerateVulgar

def readall(f) :
    tmp = []
    fq = FastqFile(f)
    fq.open()

    for s in fq :
        s.id = s.id.split()[0] # correct id to only first token
        tmp.append(s)

    fq.close()

    return tmp

def totalseq(contigs) :
    return sum(map(lambda x : len(x), contigs))

def n50(contigs, total) :
    tmp = 0

    for c in sorted(contigs, key=len, reverse=True) :
        tmp += len(c)
        if tmp > (total / 2) :
            return len(c)

    raise "should not be able to get here"

def proportion_aligned(exonerate, contigs) :
    results = exonerate.query_sequences(contigs)

    #for i in results :
    #    print "key = \"%s\", value = \"%s\"" % (i, results[i])

    total_aligned = 0
    total_sequence = 0
    perfect_contigs = 0

    for contig in contigs :
        if results[contig.id] != None :
            qalign = results[contig.id].query
            aligned = qalign.stop - qalign.start
            total_aligned += aligned

            if aligned == len(contig) :
                perfect_contigs += 1

        total_sequence += len(contig)

    return (len(contigs), perfect_contigs, perfect_contigs / float(len(contigs)), 
           total_sequence, total_aligned, total_aligned / float(total_sequence))

def main() :
    if len(sys.argv) != 2 :
        print >> sys.stderr, "Usage: %s <DIR>\n" % sys.argv[0]
        sys.exit(-1)

    directory = sys.argv[1]

    # assembly.SRR535750.m50f0q30.uniq.k41.kt9.mo45.ao45-contigs.fa
    p = re.compile("assembly.SRR535750.m50f0q30.uniq.k(\d+).kt(\d+).mo(\d+).ao(\d+)-contigs.fa")
    delim = "\t"
    read_length = 50
    contig_threshold = 200

    fields = ["k","kmer_threshold","merge_overlap","assemble_overlap",
              "num_contigs","perfect_contigs","proportion_perfect",
              "total_sequence","aligned_sequence","proportion_aligned",
              "n50"]
    print delim.join(fields)

    exonerate = Exonerate()
    exonerate.start(directory)

    for f in glob.glob("*contigs.fa") :
        m = p.match(f)
        if m is None :
            continue

        #d = dict(zip(['k','kt','mo','ao'], map(int, m.groups())))
        contigs = readall(f)
        contigs = filter(lambda x : len(x) >= contig_threshold, contigs)

        # XXX debug
        #contigs = contigs[:10]
        # XXX debug

        #v = map(int, list(m.groups()))
        v = list(m.groups())
        num_contigs, perfect_contigs, prop_perfect, sequence, aligned, prop_aligned = proportion_aligned(exonerate, contigs)
        v.append(num_contigs)
        v.append(perfect_contigs)
        v.append(prop_perfect)
        v.append(sequence)
        v.append(aligned)
        v.append(prop_aligned)
        v.append(n50(contigs, sequence))

        print delim.join(map(str, v))

    exonerate.stop()

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by user..."
        sys.exit(-1)

