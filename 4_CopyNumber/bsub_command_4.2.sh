#!/bin/bash

MEM=15000
CPUS=1
QUEUE=long

mkdir -p logs_4.2

NSAMP=`wc -l ../data/original/CrossSpecies_ProjectInfo.txt | cut -f1 -d' '`
NSAMP=$(( NSAMP - 1 ))

for IDX in `seq 1 $NSAMP`; do

    CMD="Rscript 4.2_CN_AlleleSpecific.R $IDX"

    bsub -G team78-grp -o logs_4.2/log.%J -q $QUEUE -n $CPUS -R "span[hosts=1] select[mem>=$MEM] rusage[mem=$MEM]" -M $MEM "$CMD"

done
