import sys
import os
from glob import glob
from subprocess import check_call, CalledProcessError


threads = 10
DEVNULL = None
verbose = True

def run(command) :
    global DEVNULL, verbose

    if verbose :
        print "Running: %s" % command

    try :
        check_call(command.split(), stdout=DEVNULL, stderr=DEVNULL)

    except CalledProcessError, cpe :
        print >> sys.stderr, "Error: %s" % command
        print >> sys.stderr, str(cpe)
        sys.exit(-1)

def run_shell(command) :
    global DEVNULL, verbose

    if verbose :
        print "Running: %s" % command

    try :
        check_call(command, stdout=DEVNULL, stderr=DEVNULL, shell=True)

    except CalledProcessError, cpe :
        print >> sys.stderr, "Error: %s" % command
        print >> sys.stderr, str(cpe)
        sys.exit(-1)

def rm(f) :
    try :
        os.remove(f)
    except :
        pass

def cleanup(files) :
    global verbose

    if verbose :
        print >> sys.stderr, "Deleting temporary files..."

    for f in files :
        rm(f)

def clipext(f) :
    return '.'.join(f.split('.')[:-1])

def sga_init(fastq, m, f, q) :
    global threads

    n = clipext(fastq)
    o = "%s.m%df%dq%d.fq" % (n, m, f, q)
    u = "%s.m%df%dq%d.uniq.fq" % (n, m, f, q)

    run("sga preprocess -p 0 -m %d -f %d -q %d -o %s %s" % (m, f, q, o, fastq))
    run_shell("grep -B1 '^\+$' %s | grep -v '^\+$' | grep -v '^--$' | sort | uniq | python ./seq2fa.py > %s" % (o, u))
    #run("sga index -a ropebwt -t %d --no-reverse %s" % (threads, u))

    return u

# XXX should be filter out everything that cannot be corrected?
def sga_correct(fasta, k, kmer_threshold) :
    global threads

    #run("sga correct -k %d -x %d --discard -t %d -o correct.fa %s" % (k, kmer_threshold, threads, fasta))
    #run("sga index -a ropebwt -t %d correct.fa" % (threads))
    #run("sga filter -x 2 -t %d --homopolymer-check --low-complexity-check correct.fa" % (threads))

    run_shell("cp %s correct.fa" % (fasta))
    run("sga index -a ropebwt -t %d correct.fa" % (threads))
    if kmer_threshold != 0 :
        run("sga filter -k %d -x %d -t %d --homopolymer-check --low-complexity-check correct.fa" % (k, kmer_threshold, threads))
    else :
        run("sga filter --no-kmer-check -t %d --homopolymer-check --low-complexity-check correct.fa" % (threads))

def sga_build(merge_overlap) :
    global threads

    run("sga fm-merge -m %d -t %d -o merged.fa correct.filter.pass.fa" % (merge_overlap, threads))
    run("sga index -d 1000000 -t %d merged.fa" % (threads))
    run("sga rmdup -t %d merged.fa" % (threads))
    run("sga overlap -m %d -t %d merged.rmdup.fa" % (merge_overlap, threads))

# XXX -g and -r
def sga_assemble(fasta, assemble_overlap, out) :
    run("sga assemble -m %d -g 0 -r 10 -o %s merged.rmdup.asqg.gz" % (assemble_overlap, out))

def main() :
    global DEVNULL

    if len(sys.argv) != 2 :
        print >> sys.stderr, "Usage: %s <FASTQ>\n" % sys.argv[0]
        return -1

    DEVNULL = open('/dev/null', 'w')

    f = sga_init(sys.argv[1], 50, 0, 30)

    k_values =  [23,25,27,29,31,33,35,37,39,41]
    kt_values = [0,1,2,3,4,5,6,7,8,9]
    mo_values = [20,25,30,35,40,45]
    ao_values = [20,25,30,35,40,45]

    for kmer in k_values :
        print >> sys.stderr, "Info: k=%d" % kmer
        for kt in kt_values :
            print >> sys.stderr, "Info: k=%d kt=%d" % (kmer, kt)
            sga_correct(f, k=kmer, kmer_threshold=kt)

            for mo in mo_values :
                print >> sys.stderr, "Info: k=%d kt=%d mo=%d" % (kmer, kt, mo)
                sga_build(merge_overlap=mo)

                for ao in ao_values :
                    if ao < mo :
                        continue

                    print >> sys.stderr, "Info: k=%d kt=%d mo=%d ao=%d" % (kmer, kt, mo, ao)
                    outfile = "assembly.%s.k%d.kt%d.mo%d.ao%d" % (clipext(f), kmer, kt, mo, ao)
                    sga_assemble(f, assemble_overlap=ao, out=outfile)

    cleanup(glob("correct.*") + glob("merged.*") + glob(clipext(clipext(f)) + ".*"))

    DEVNULL.close()

    return 0

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        print >> sys.stderr, "Killed by user..."
        sys.exit(-1)

