#!/bin/bash
set -e

SPECIES='D.melanogaster'
RELEASE=72
MINLENGTH=50
THREADS=4

#for i in `ls -d trinity_results/n_dubi* | xargs -I {} basename {}` ; do
for i in `cat dubitatus.txt | grep 'fasta' | cut -d'/' -f1,2 | xargs -I {} basename {}` ; do
    python hoard.py align \
                    -s $SPECIES \
                    -r $RELEASE \
                    -m $MINLENGTH \
                    -i trinity_results/$i/Trinity.fasta \
                    -c $THREADS \
                    -o ${i}_${SPECIES}_${RELEASE} \
                    -t /tmp/fly
done

