# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 3.1: CALCULATE BIN COVERAGE IN EACH SAMPLE


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CHROM.LISTS = "../data/original/CrossSpecies_ChromLists.txt",
    BAM.PATHS = "../data/original/Path_SampleBAMs.txt"
)

# Output file paths
OUTPUT = list(
    LOG.DIR = "logs_3.1",
    BINS.DIR = "../data/processed/GenomeBins",
    COV.DIR = "../data/processed/BinCoverage",
    BINS.PREFIX = "GenomeBins_100kb_",
    COV.PREFIX = "BinCov_"
)


# # Species to consider
# SPECIES = c("mouse", "rat", "cow", "dog", "horse", "cat")
# 
# # Minimum contig length (bp)
# MIN.CONTIG.LEN = 11e6

# Genome bin length (bp)
BIN.LEN = 1e5

# Memory, queue and command templates for job submission
MEM = 15000
QUEUE = "long"
COV.CMD = "coverageBed -b ${BINS} -abam ${BAM} > ${OUTFILE}"
BSUB.CMD = "bsub -G team78-grp -o ${LOG}/log.%J -q ${QUEUE} -n 1 -R \"span[hosts=1] select[mem>=${MEM}] rusage[mem=${MEM}]\" -M ${MEM} \"${CMD}\""


# Create output and log directories
dir.create(OUTPUT$LOG.DIR, showWarnings=F)
dir.create(OUTPUT$BINS.DIR, showWarnings=F)
dir.create(OUTPUT$COV.DIR, showWarnings=F)

# Disable scientific notation
options(scipen=999)


cat("Loading data...\n")
bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
chrom.lists = read.table(INPUT$CHROM.LISTS, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# For each species
for (species in unique(chrom.lists$Species)) {
    cat("\nProcessing species:", species, "\n")
    species.idx = sample.info$SPECIES_NAME == species
    ref.name = sample.info$REFERENCE_GENOME[species.idx][1]
    project = sample.info$PROJECT_ID[species.idx][1]
    
    # # Load genome index (fai)
    # genome.fai = read.table(paste0(gsub("${SPECIES}", species,
    #                                     gsub("${REFGENOME}", ref.name, genome.path, fixed=T),
    #                                     fixed=T), ".fai"),
    #                         sep="\t", col.names=c("Chrom", "Length", "x", "y", "z"), as.is=T)
    
    # Create BED file of non-overlapping bins along all chromosomes
    #chrom.lengths = genome.fai[genome.fai$Length > MIN.CONTIG.LEN, c("Chrom", "Length")]
    chrom.lengths = chrom.lists[chrom.lists$Species == species, c("Chromosome", "End")]
    bin.table = NULL
    for (i in 1:nrow(chrom.lengths)) {
        bin.table = rbind(bin.table,
                          cbind(chrom.lengths$Chromosome[i],
                                seq(0, chrom.lengths$End[i] - BIN.LEN, by=BIN.LEN),
                                seq(BIN.LEN, chrom.lengths$End[i], by=BIN.LEN)))
    }
    bins.path = paste0(OUTPUT$BINS.DIR, "/", OUTPUT$BINS.PREFIX, ref.name, ".bed")
    write.table(bin.table, file=bins.path, sep="\t", quote=F, row.names=F, col.names=F)
    
    # Use coverageBed to measure coverage per bin in each sample
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
        
        cmd = gsub("${BINS}", bins.path,
                   gsub("${BAM}", bam,
                        gsub("${OUTFILE}",
                             paste0(OUTPUT$COV.DIR, "/", OUTPUT$COV.PREFIX, sample.id, ".txt"),
                             COV.CMD, fixed=T), fixed=T), fixed=T)
        
        system(gsub("${QUEUE}", QUEUE,
                    gsub("${MEM}", MEM,
                         gsub("${LOG}", OUTPUT$LOG.DIR,
                              gsub("${CMD}", cmd,
                                   BSUB.CMD, fixed=T), fixed=T), fixed=T), fixed=T))
    }
}

cat("\nDone\n")
