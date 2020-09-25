# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 5.6: MERGE REGIONS ABOVE MAX. COV. PERCENTILE IN EACH SAMPLE AND THOSE IN ITS MATCHED NORMAL


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    COV.PREFIX = "../data/processed/CoverageFilters/CovAbove99Pct_"
)

# Output file paths
OUTPUT = list(
    LOG.DIR = "logs_5.6",
    COV.PREFIX = "../data/processed/CoverageFilters/CovAbove99Pct_NormalUnion_"
)


# Memory, queue and command templates for job submission
MEM = 1000
QUEUE = "small"
BED.CMD = "cat ${SAMPLE} ${NORMAL} | sort -k1,1 -k2,2n | bedtools merge > ${OUTFILE}"
BSUB.CMD = "bsub -G team78-grp -o ${LOG}/log.%J -q ${QUEUE} -n 1 -R \"span[hosts=1] select[mem>=${MEM}] rusage[mem=${MEM}]\" -M ${MEM} \"${CMD}\""


# Create log directory
dir.create(OUTPUT$LOG.DIR, showWarnings=F)


cat("Loading data...\n")
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# For each sample, use bedtools to merge the regions in sample and matched normal
for (i in 1:nrow(sample.info)) {
    sample.id = sample.info$SAMPLE_NAME[i]
    normal.id = sample.info$NORMAL_NAME[i]
    cat("Processing sample ", sample.id, "\n", sep="")
    file.names = paste0(INPUT$COV.PREFIX, c(sample.id, normal.id), ".bed")
    
    if (all(file.exists(file.names))) {
        cmd = gsub("${SAMPLE}", file.names[1],
                   gsub("${NORMAL}", file.names[2],
                        gsub("${OUTFILE}", paste0(OUTPUT$COV.PREFIX, sample.id, ".bed"),
                             BED.CMD, fixed=T), fixed=T), fixed=T)
        
        system(gsub("${QUEUE}", QUEUE,
                    gsub("${MEM}", MEM,
                         gsub("${LOG}", OUTPUT$LOG.DIR,
                              gsub("${CMD}", cmd,
                                   BSUB.CMD, fixed=T), fixed=T), fixed=T), fixed=T))
    }
    else {
        for (j in which(!file.exists(file.names))) {
            cat("WARNING: File", file.names[j], "not found\n")
        }
    }
}

cat("\nDone\n")
