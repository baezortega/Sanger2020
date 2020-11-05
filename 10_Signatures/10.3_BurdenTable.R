# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 10.3: MAKE TABLE OF BURDEN, SIGNATURE EXPOSURES AND [C>T]pG PER SAMPLE


# Input file paths
INPUT = list(
    GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    SIGS = "../data/processed/Signatures_Definitive.RData",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CALLABLE.PREFIX = "../data/processed/CallableGenome/CallableGenome_"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/Burden_Exposures.RData",
    TABLE = paste0("../output/", Sys.Date(), "_Burden_Exposures_Table.txt")
)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(Biostrings))
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
load(INPUT$SIGS)
rownames(counts.all$sample) = sapply(strsplit(rownames(counts.all$sample), " "), `[`, 1)
rownames(exposures.final$sample$mean) = sapply(strsplit(rownames(exposures.final$sample$mean), " "),
                                               `[`, 1)
cat("Loaded\n")


# # ADD HUMAN SAMPLES WITH EVIDENT COLIBACTIN EXPOSURE TO THE EXCLUSION LIST
# exclude.list = c(exclude.list,
#                  sample.info$SAMPLE_NAME[sample.info$NORMAL_NAME %in%
#                                              c("PD37449b", "PD37512a4", "PD37513a21")])
stopifnot(identical(rownames(counts.all$sample),
                    sample.info$SAMPLE_NAME[!(sample.info$SAMPLE_NAME %in% exclude.list)]))


# Build table of burden and exposures per sample
burden.expos = NULL
for (species in unique(sample.info$SPECIES_NAME)) {
    # Load reference genome
    cat("\nProcessing species:", species, "\n")
    cat("Loading reference genome\n")
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    ref.name = sample.info$REFERENCE_GENOME[species.idx][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species,
                                   gsub("${REFGENOME}", ref.name, genome.path, fixed=T), fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    for (i in which(species.idx)) {
        sample.id = sample.info$SAMPLE_NAME[i]
        normal.id = sample.info$NORMAL_NAME[i]
        if (sample.id %in% exclude.list) {
            cat("WARNING: Sample", sample.id, "is in the exclusion list\n")
            next
        }
        cat("Processing sample", sample.id, "\n")
        
        # Load callable genome regions
        regions = read.table(paste0(INPUT$CALLABLE.PREFIX, sample.id, ".txt"),
                             sep="\t", header=T, as.is=T)
        genome.regions = padAndClip(genome[regions$chrom], IRanges(regions$start, regions$end),
                                    Lpadding.letter=".", Rpadding.letter=".")
        dinuc.freqs = colSums(dinucleotideFrequency(genome.regions))
        
        # Add data to table
        burden.expos = rbind(burden.expos,
                             c("Sample" = sample.id,
                               "Mutations" = sum(counts.all$sample[sample.id, ]),
                               "[C>T]pG" = sum(counts.all$sample[sample.id, c("ACG>ATG", "CCG>CTG",
                                                                              "GCG>GTG", "TCG>TTG")]),
                               "SigA_Exposure" = exposures.final$sample$mean[sample.id, 1],
                               "SigB_Exposure" = exposures.final$sample$mean[sample.id, 2],
                               "SigC_Exposure" = exposures.final$sample$mean[sample.id, 3],
                               "Callable_Genome_Bp" = sum(regions$width),
                               "CpG_Frequency" = as.numeric(dinuc.freqs["CG"])))
    }
}


# Output table
write.table(burden.expos, file=OUTPUT$TABLE, sep="\t", quote=F, row.names=F, col.names=T)
save(burden.expos, file=OUTPUT$DATA)

cat("\nDone\n")
