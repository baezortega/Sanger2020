# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 5.1: OBTAIN GENOME COVERAGE HISTOGRAMS FOR EACH SAMPLE


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    BAM.PATHS = "../data/original/Path_SampleBAMs.txt"
)

# Output file paths
OUTPUT = list(
    LOG.DIR = "logs_5.1",
    OUT.DIR = "../data/processed/CoverageHist",
    COV.PREFIX = "CoverageHist_"
)


# Memory, queue and command templates for job submission
MEM = 5000
QUEUE = "normal"
COV.CMD = "bedtools genomecov -ibam ${BAM} | grep genome > ${OUTFILE}"
BSUB.CMD = "bsub -G team78-grp -o ${LOG}/log.%J -q ${QUEUE} -n 1 -R \"span[hosts=1] select[mem>=${MEM}] rusage[mem=${MEM}]\" -M ${MEM} \"${CMD}\""


# Create output and log directories
dir.create(OUTPUT$LOG.DIR, showWarnings=F)
dir.create(OUTPUT$OUT.DIR, showWarnings=F)


cat("Loading data...\n")
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# For each species
for (species in unique(sample.info$SPECIES_NAME)) {
    cat("\nProcessing species:", species, "\n")
    species.idx = sample.info$SPECIES_NAME == species
    ref.name = sample.info$REFERENCE_GENOME[species.idx][1]
    project = sample.info$PROJECT_ID[species.idx][1]
    
    # Use bedtools genomecov to obtain coverage histograms for each sample
    for (sample.id in unique(c(sample.info$SAMPLE_NAME[species.idx],
                               sample.info$NORMAL_NAME[species.idx]))) {
        
        cat("Processing sample ", sample.id, "\n", sep="")
        
        if (species %in% bam.paths$SPECIES) {
            bam.path = bam.paths$PATH[bam.paths$SPECIES == species]
        }
        else {
            bam.path = bam.paths$PATH[bam.paths$SPECIES == "default"]
        }
        
        bam = gsub("${SPECIES}", species,
                   gsub("${REFGENOME}", ref.name, 
                        gsub("${PROJECT}", project,
                             gsub("${SAMPLE}", sample.id,
                                  bam.path, fixed=T), fixed=T), fixed=T), fixed=T)
        
        if (file.exists(bam)) {
            cmd = gsub("${BAM}", bam,
                       gsub("${OUTFILE}",
                            paste0(OUTPUT$OUT.DIR, "/", OUTPUT$COV.PREFIX, sample.id, ".txt"),
                            COV.CMD, fixed=T), fixed=T)
            
            system(gsub("${QUEUE}", QUEUE,
                        gsub("${MEM}", MEM,
                             gsub("${LOG}", OUTPUT$LOG.DIR,
                                  gsub("${CMD}", cmd,
                                       BSUB.CMD, fixed=T), fixed=T), fixed=T), fixed=T))
        }
        else {
            cat("WARNING: File", bam, "not found\n")
        }
    }
}

cat("\nDone\n")
