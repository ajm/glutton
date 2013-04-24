import sys
import glob
import os
import subprocess

from lib.datatypes import Sequence
from lib.filetypes import FastqFile


exonerate = None

def start_exonerate(d) :
    global exonerate

    index_file = os.path.join(d, 'exonerate.esi')
    args = ['exonerate-server', 'exonerate.esi', '--port', '12887']

    if not (os.path.exists(index_file) and os.path.isfile(index_file)) :
        print >> sys.stderr, "Error: %d does not exist..." % index_file
        sys.exit(-1)

    exonerate = subprocess.Popen(args, 
                    cwd=d,
                    stderr=subprocess.STDOUT, 
                    stdout=subprocess.PIPE)

    # wait until it is listening...
    out = []
    for line in iter(exonerate.stdout.readline, '') :
        out.append(line)

        if 'listening on port' in line :
            print >> sys.stderr, "Info: exonerate-server started..."
            break
    else :
        print >> sys.stderr, "Error: exonerate-server did not start properly:\n\n\t%s\n" % '\t'.join(out)
        sys.exit(-1)

    # change the redirect so the pipes buffer does not fill and stall
    # everything
    exonerate.stdout = open('/dev/null', 'w')

def stop_exonerate() :
    global exonerate

    if exonerate :
        exonerate.terminate()

def get_genes(fname) :
    tmp = {}
    f = FastqFile(fname)
    f.open()

    for seq in f :
        tmp[seq.id[1:]] = None

    f.close()
    return tmp

def query_exonerate(fname) :
    args = ['exonerate', fname, 'localhost:12887', '--bestn', '1', '--model', 'affine:local']

    tmp = get_genes(fname)

    try :
        out = subprocess.check_output(
                            args, 
                            stderr=open('/dev/null', 'w'))

    except subprocess.CalledProcessError, cpe:
        #print >> sys.stderr, cpe
        #sys.exit(-1)
        return tmp

    for line in out.split('\n') :
        if not line.startswith('vulgar') :
            continue
        
        # vulgar: aael001690 581 737 + >agap005310 596 752 + 339 M 156 156
        # vulgar: aael009244 667 744 + >agap006707 604 681 + 178 M 77 77
        fields = line.strip().split()
        tmp[fields[1]] = fields[5][1:]

    return tmp

def usage() :
    print >> sys.stderr, "Usage: %s <DIR>" % sys.argv[0]
    print >> sys.stderr, "\twhere DIR is the directory containing all files (eg: ./cache/17/A.gambiae/)\n"

def main() :
    if len(sys.argv) != 2 :
        usage()
        return -1

    directory = sys.argv[1]

    start_exonerate(directory)

    for ofile in glob.glob(os.path.join(directory, 'ortholog_*')) :
        print >> sys.stderr, ofile
        results = query_exonerate(ofile)

        for qid,rid in results.iteritems() :
            print os.path.basename(ofile), qid, rid 

    stop_exonerate()

    return 0

if __name__ == '__main__' :
    try :
        sys.exit(main())

    except KeyboardInterrupt :
        stop_exonerate()

