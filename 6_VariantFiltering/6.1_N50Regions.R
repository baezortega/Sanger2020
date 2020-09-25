# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 6.1: OBTAIN REGIONS CORRESPONDING TO TRACTS OF 50+ 'N'S FOR EACH GENOME


# Input file paths
INPUT = list(
    GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt"
)

# Output file paths
OUTPUT = list(
    OUT.DIR = "../data/processed/N50Filter",
    BED.PREFIX = "N50plus_1kb_"
)


# Create output directory
dir.create(OUTPUT$OUT.DIR, showWarnings=F)

# Disable scientific notation
options(scipen=999)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(Biostrings))
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n\n")


for (species in unique(sample.info$SPECIES_NAME)) {
    cat("Processing species:", species, "\n")

    # Load reference genome
    ref.name = sample.info$REFERENCE_GENOME[sample.info$SPECIES_NAME == species][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species,
                                   gsub("${REFGENOME}", ref.name, 
                                        genome.path, fixed=T), fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    # Find tracts of 50 consecutive Ns
    match.N50 = vmatchPattern("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", genome)
    
    # Merge regions and extend by 1kb in each direction
    match.N50.merged = NULL
    for (i in 1:length(match.N50)) {
        if (length(match.N50[[i]]) > 0) {
            regions = reduce(match.N50[[i]])
            start(regions) = sapply(start(regions), function(s) max(s - 1000, 1))
            end(regions) = end(regions) + 1000
            regions = reduce(regions)
            
            match.N50.merged = rbind(match.N50.merged,
                                     cbind(chr=names(match.N50)[i],
                                           as.data.frame(regions)[, 1:2]))
        }
    }
    
    # Transform to 0-based and output to BED file
    match.N50.merged$start = match.N50.merged$start - 1
    write.table(match.N50.merged, sep="\t", quote=F, col.names=F, row.names=F,
                file=paste0(OUTPUT$OUT.DIR, "/", OUTPUT$BED.PREFIX, species, ".bed"))
    
    invisible(gc())
}

cat("\nDone\n")
