#!/bin/bash

MAX_NJOBS=45
SLEEP_TIME=240 # secs

COUNTER=1

for file in *.fasta; do 
    outdir=`basename $file .fasta`
    if [ ! -e $outdir ] && [ `grep -c '>' $file` -lt 1200 ]; then
        let COUNTER++

        echo $outdir
        mkdir $outdir
        cd $outdir
        mv ../$file .
        nextflow -q -bg -c ../nextflow.config run ../nemu-core.nf -with-trace --sequence $file --outdir .
        cd -
    fi
    
    if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
        echo sleep $SLEEP_TIME secs
        sleep $SLEEP_TIME
    fi

done

echo DONE
