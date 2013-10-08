import sys
import glob
import re
import os

from lib.progress import Progress
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

def merge_ranges(ranges) :
    mins = []
    maxs = []

    for start,stop in ranges :
        mins.append(start)
        maxs.append(stop)

    return (min(mins), max(maxs))

def ranges_overlap(r1, r2) :
    if r1[0] <= r2[0] :
        return r1[1] >= r2[0]
    else :
        return r2[1] >= r1[0]

def add_range(ranges, new_range) :
    # find all ranges that overlap, merge this set of ranges
    overlap = [new_range]
    nonoverlap = []

    for i in ranges :
        if ranges_overlap(i, new_range) :
            overlap.append(i)
        else :
            nonoverlap.append(i)

    return nonoverlap + [ merge_ranges(overlap) ]

def add_match_range(match, match_dict) :
    if not match_dict.has_key(match.seqid) :
        match_dict[match.seqid] = [(match.start, match.stop)]
    else :
        match_dict[match.seqid] = add_range(match_dict[match.seqid], (match.start, match.stop))

def sum_ranges(match_dict) :
    return sum(
                map(lambda x : sum(map(lambda y : y[1] - y[0], x)), match_dict.values())
              )

def print_ranges(ranges) :
    for i in ranges :
        if len(ranges[i]) > 5 :
            print i, sorted(ranges[i])

def proportion_aligned(exonerate, contigs, total_transcriptome_length, min_alignment_length) :
    results = exonerate.query_sequences(contigs)

    #for i in results :
    #    print "key = \"%s\", value = \"%s\"" % (i, results[i])

    total_aligned = 0
    total_sequence = 0
    perfect_contigs = 0
    perfect_sequence = 0

    match_ranges = {}

    for contig in contigs :
        total_sequence += len(contig)

        if results[contig.id] != None :
            qalign = results[contig.id].query
            aligned = qalign.stop - qalign.start

            if aligned < min_alignment_length :
                continue

            total_aligned += aligned

            if aligned == len(contig) :
                perfect_contigs += 1
                perfect_sequence += aligned

            # alignment may be on the negative strand
            # but all this code assumes that for range (a,b)
            # a < b
            m = results[contig.id].match
            if m.stop < m.start :
                m.start, m.stop = m.stop, m.start


            add_match_range(m, match_ranges)

    #print print_ranges(match_ranges)
    coverage = sum_ranges(match_ranges) / float(total_transcriptome_length)

    return (len(contigs), perfect_contigs, total_sequence, perfect_sequence, total_aligned, coverage)

def main() :
    if len(sys.argv) != 4 :
        print >> sys.stderr, "Usage: %s <DIR> <min contig length> <min alignment length>\n" % sys.argv[0]
        sys.exit(-1)

    directory = sys.argv[1]
    min_contig_length = int(sys.argv[2])
    min_alignment_length = int(sys.argv[3])

    os.system('killall exonerate-server &> /dev/null')

    # assembly.SRR535750.m50f0q30.uniq.k41.kt9.mo45.ao45-contigs.fa
    #p = re.compile("assembly.SRR535750.m50f0q30.uniq.k(\d+).kt(\d+).mo(\d+).ao(\d+)-contigs.fa")
    p = re.compile("assembly.[^.]+.m\d+f\d+q\d+.uniq.k(\d+).kt(\d+).mo(\d+).ao(\d+)-contigs.fa")
    delim = "\t"

    fields = ["k","kmer_threshold","merge_overlap","assemble_overlap",
              "num_contigs","total_sequence","perfect_contigs","perfect_sequence",
              "aligned_sequence","n50","coverage",
              "prop_perfect_contigs","prop_perfect_sequence","prop_aligned"]

    print delim.join(fields)

    exonerate = Exonerate()
    exonerate.start(directory)

    transcriptome_length = totalseq(readall(os.path.join(directory, 'exonerate.fa')))

    print >> sys.stderr, "transcriptome length = %d" % transcriptome_length

    contig_files = glob.glob("*contigs.fa")

    pb = Progress("Progress:", len(contig_files))
    pb.start()

    for f in contig_files :
        m = p.match(f)
        if m is None :
            continue

        #d = dict(zip(['k','kt','mo','ao'], map(int, m.groups())))
        contigs = readall(f)
        contigs = filter(lambda x : len(x) >= min_contig_length, contigs)

        # XXX debug
        #contigs = contigs[:3]
        # XXX debug

        #v = map(int, list(m.groups()))
        v = list(m.groups())
        num_contigs, perfect_contigs, sequence, perfect_sequence, aligned, coverage = \
                proportion_aligned(exonerate, contigs, transcriptome_length, min_alignment_length)
        v.append(num_contigs)
        v.append(sequence)
        v.append(perfect_contigs)
        v.append(perfect_sequence)
        v.append(aligned)
        v.append(n50(contigs, sequence))
        v.append(coverage)
        v.append(perfect_contigs / float(num_contigs))
        v.append(perfect_sequence / float(sequence))
        v.append(aligned / float(sequence))

        print delim.join(map(str, v))

        pb.increment()

    pb.end()
    exonerate.stop()

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by user..."
        sys.exit(-1)
