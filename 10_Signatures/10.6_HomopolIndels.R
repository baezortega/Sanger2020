# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 10.6: PRODUCE COUNTS OF INDELS AT HOMOPOLYMER TRACTS IN EACH SAMPLE


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    VARS.FINAL = "${AC36}/final_indel_calls/${SPECIES}/${SAMPLE}_final_indel_calls.txt"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/HomopolymerIndels_Counts_All.RData",
    COUNTS = paste0("../output/", Sys.Date(), "_HomopolymerIndels_All.txt"),
    PLOT = paste0("../output/", Sys.Date(), "_HomopolymerIndels_PerSpecies.pdf")
)


# Define minimum homopolymer length and homopolymer sequences
MIN.LEN = 5
HOMOPOL = sapply(c("A", "C", "G", "T"), function(x) paste(rep(x, MIN.LEN), collapse=""))


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Initialise count table
indel.counts = sample.names = NULL

for (species in unique(sample.info$SPECIES_NAME)) {
    cat("\nProcessing species:", species, "\n")
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    
    # Load reference genome
    cat("Loading reference genome\n")
    ref.name = sample.info$REFERENCE_GENOME[species.idx][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species,
                                   gsub("${REFGENOME}", ref.name, genome.path, fixed=T), fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    # Process samples
    for (sample.id in sample.info$SAMPLE_NAME[species.idx]) {
        
        # Load filtered indels
        cat("Processing sample", sample.id, "\n")
        var.path = gsub("${AC36}", ac36.path,
                        gsub("${SAMPLE}", sample.id,
                             gsub("${SPECIES}", species,
                                  INPUT$VARS.FINAL, fixed=T), fixed=T), fixed=T)
        
        if (!file.exists(var.path)) {
            cat("WARNING: File", var.path, "not found. Skipping sample.\n")
            next
        }
        vars = read.table(var.path, sep="\t", header=T, check.names=F, as.is=T,
                          colClasses=c("Chrom"="character", "Ref"="character", "Alt"="character"))
        
        # Retrieve 11-bp sequence contexts
        # (context goes between bases -4 and +6 to account for the fact that the
        # variant is called at the position immediately before the insertion/deletion)
        vars$Context = as.character(padAndClip(genome[vars$Chr],
                                               IRanges(vars$Pos - 4, vars$Pos + 6),
                                               Lpadding.letter=".", Rpadding.letter="."))
        stopifnot(identical(substr(vars$Context, 5, 5), substr(vars$Ref, 1, 1)))
        
        # Search each homopolymer sequence in each indel's context
        homopol.idx = sapply(HOMOPOL, function(homopol) {
            grepl(homopol, vars$Context, fixed=T)
        })
        
        counts = colSums(homopol.idx)
        for (i in 1:(ncol(homopol.idx)-1)) {
            for (j in (i+1):ncol(homopol.idx)) {
                counts = c(counts, sum(homopol.idx[, i] & homopol.idx[, j]))
                names(counts)[length(counts)] = paste(names(HOMOPOL)[c(i, j)], collapse="+")
            }
        }
        counts = c("Any" = sum(rowSums(homopol.idx) > 0), counts)
        names(counts) = paste0("Indels_Homopol_", names(counts))
        counts = c("Indels_Total" = nrow(vars), counts)
        
        indel.counts = rbind(indel.counts, counts)
        sample.names = c(sample.names, sample.id)
    }
}
rownames(indel.counts) = sample.names


# Plot indel counts of each type per species
pdf(OUTPUT$PLOT, 10, 6)
par(mar=c(6, 4.5, 4, 1.5), mgp=c(3, 0.8, 0))
for (species in unique(sample.info$SPECIES_NAME)) {
    species.cnt = indel.counts[rownames(indel.counts) %in%
                                   sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species], ]
    homopol.pct = round(median(species.cnt[, 2] / species.cnt[, 1]) * 100, 2)
    boxplot(species.cnt, col="skyblue",
            las=2, names=rep("", ncol(species.cnt)), ylab="Number of indels",
            main=paste0("Indels in ", species, " (", nrow(species.cnt), " samples)\n",
                        homopol.pct, "% homopolymer indels (median)"))
    text(x=1:ncol(species.cnt), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
         labels=str_split_fixed(colnames(species.cnt), "_", 2)[, 2], srt=45, adj=1, xpd=TRUE)
}
invisible(dev.off())


# Output and save count table
cat("\nSaving indel count table\n")
write.table(cbind("Sample"=sample.names, indel.counts),
            file=OUTPUT$COUNTS, sep="\t", quote=F, row.names=F)
save(indel.counts, file=OUTPUT$DATA)

cat("\nDone\n")
