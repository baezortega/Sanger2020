#!/bin/bash

MEM=5000
CPUS=1
QUEUE=normal

mkdir -p logs_6.6

NSAMP=`wc -l ../data/original/CrossSpecies_ProjectInfo.txt | cut -f1 -d' '`
NSAMP=$(( NSAMP - 1 ))

for IDX in `seq 1 $NSAMP`; do

    CMD="Rscript 6.6_IndelFiltering.R $IDX"

    bsub -G team78-grp -o logs_6.6/log.%J -q $QUEUE -n $CPUS -R "span[hosts=1] select[mem>=$MEM] rusage[mem=$MEM]" -M $MEM "$CMD"

done
