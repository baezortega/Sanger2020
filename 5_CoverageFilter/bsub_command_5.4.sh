#!/bin/bash

MEM=500
CPUS=1
QUEUE=small

CMD="Rscript 5.4_MinCoverage_Union.R"

bsub -G team78-grp -o log.%J -q $QUEUE -n $CPUS -R "span[hosts=1] select[mem>=$MEM] rusage[mem=$MEM]" -M $MEM "$CMD"
