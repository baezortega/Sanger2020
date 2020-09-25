#!/bin/bash

MEM=20000
CPUS=1
QUEUE=long

CMD="Rscript 10.1_Signatures.R"

bsub -G team78-grp -o log.%J -q $QUEUE -n $CPUS -R "span[hosts=1] select[mem>=$MEM] rusage[mem=$MEM]" -M $MEM "$CMD"
