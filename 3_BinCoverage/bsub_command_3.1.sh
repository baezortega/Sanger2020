#!/bin/bash

MEM=500
CPUS=1
QUEUE=normal

CMD="Rscript 3.1_BinCoverage.R"

bsub -G team78-grp -o log.%J -q $QUEUE -n $CPUS -R "span[hosts=1] select[mem>=$MEM] rusage[mem=$MEM]" -M $MEM "$CMD"


## TO RESUBMIT FAILED INDIVIDUAL JOBS:

# QUEUE=long; MEM=20000; CPUS=1
# for FILE in `grep -m1 "Exited" logs*/* | cut -f1 -d:`; do
#     CMD=`grep -m1 "Exited" $FILE | cut -f2 -d'<' | cut -f1,2 -d'>'`
#     echo $CMD
#     bsub -o logs_covBed/log.%J -G team78-grp -q $QUEUE -n $CPUS -R "span[hosts=1] select[mem>=$MEM] rusage[mem=$MEM]" -M $MEM "$CMD"
# done

# for FILE in `grep -m1 "Exited" logs*/* | cut -f1 -d:`; do rm $FILE; done


## TO CHECK CORRECT NUMBER OF LINES IN OUTPUT:
# wc -l ../data/processed/BinCoverage/*
