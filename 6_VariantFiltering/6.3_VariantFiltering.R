# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 6.3: APPLY VARIANT FILTERS TO EVERY SAMPLE


# Read sample index from command line
SAMPLE.IDX = as.integer(commandArgs(TRUE)[1])
cat(SAMPLE.IDX, "\n\n")


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    BAM.PATHS = "../data/original/Path_SampleBAMs.txt",
    CHROM.LISTS = "../data/original/CrossSpecies_ChromLists.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    SAMPLE.MATCH = "../data/original/CrossSpecies_SamplesPerIndividual.txt",
    CNT.PREFIX = "../data/processed/AlleleCounts_Indiv/AlleleCounts_",
    N50.BED = "../data/processed/N50Filter/N50plus_1kb_${SPECIES}.bed",
    MIN.COV.BED = "../data/processed/CoverageFilters/CovBelow10_NormalUnion_${SAMPLE}.bed",
    MAX.COV.BED = "../data/processed/CoverageFilters/CovAbove99Pct_NormalUnion_${SAMPLE}.bed",
    CAVEMAN.BED = "${AC36}/${SPECIES}/${SAMPLE}/${SAMPLE}_vs_${NORMAL}.no_analysis_DEDUP.bed",
    ENDS.BED = "${AC36}/callable_genome/1kb_around_ends_of_contigs/${SPECIES}_contig_sizes_1kb_end_contigs.bed",
    VARIANTS = "${AC36}/mathijs_filters_no_shared_var_filter/${SPECIES}/${SAMPLE}/results/${SAMPLE}.txt"
)

# Output file paths
OUTPUT = list(
    RHO.DIR = "../data/processed/BetaBinom_Rho",
    CALLABLE.DIR = "../data/processed/CallableGenome",
    FINAL.DIR = "${AC36}/final_variant_calls/${SPECIES}",
    DISC.DIR = "${AC36}/discarded_variant_calls/${SPECIES}",
    RHO.PREFIX = "BetaBinom_Rho_",
    CALLABLE.PREFIX = "CallableGenome_",
    FINAL.FILE = "${SAMPLE}_final_variant_calls.txt",
    DISC.FILE = "${SAMPLE}_discarded_variant_calls_${FILTER}.txt"
)


# For each sample, we apply the following variant filters
# to the set of variants that have already passed Mathij's filters
# (Filters 1-6 are used to define callable genomes for each sample)
#  1) Discard variants in regions included in the no_analysis.bed file from Caveman
#  2) Discard variants in regions with coverage <10x in sample, matched normal, or both
#  3) Discard variants in regions with coverage >99th percentile in sample, matched normal, or both 
#  4) Discard variants within 1 kb of runs of 50 or more consecutive Ns
#  5) Discard variants within 1 kb of contig ends
#  6) Discard variants outside chromosomal contigs (only for chromosome-level assemblies)
#  7) Discard variants located at <1 kb from the previous or next variant (clustered)
#  8) Discard variants that are likely to be artefacts based on their distribution
#     across samples from the same individual (beta-binomial filter)
#  9) Discard variants with low VAF (VAF below half the median VAF of variants
#     passing the beta-binomial filter)
# 10) Discard variants without ≥1 read on each strand
# 11) Discard variants near indels


# Filter names (for output files)
FILTER.NAMES = c("caveman", "lowcov", "highcov", "N50plus", "contigends",
                 "notchrom", "clustered", "betabinom", "lowVAF", "strandbias", "indels")

# Minimum distance between variants
MIN.DIST = 1000

# Minimum value of rho for beta-binomial filter
MIN.RHO = 0.3

# Maximum distance for near-indel filter
FLANK = 10

# Maximum ratio of number of indel bases to variant bases
INDEL.FACTOR = 3

# Minimum number of supporting reads per strand
MIN.NV.STRAND = 1


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(deepSNV))
suppressPackageStartupMessages(library(GenomicRanges))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
chrom.lists = read.table(INPUT$CHROM.LISTS, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
sample.match = read.table(INPUT$SAMPLE.MATCH, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
stopifnot(!any(duplicated(sample.match$SAMPLE_NAME)))
stopifnot(all(sample.info$SAMPLE_NAME %in% sample.match$SAMPLE_NAME))
rownames(sample.match) = sample.match$SAMPLE_NAME
cat("Loaded\n")


# Create output directories
dir.create(OUTPUT$RHO.DIR, showWarnings=F)
dir.create(OUTPUT$CALLABLE.DIR, showWarnings=F)
for (species in unique(sample.info$SPECIES_NAME)) {
    dir.create(gsub("${AC36}", ac36.path,
                    gsub("${SPECIES}", species, OUTPUT$FINAL.DIR, fixed=T), fixed=T),
               recursive=T, showWarnings=F)
    dir.create(gsub("${AC36}", ac36.path,
                    gsub("${SPECIES}", species, OUTPUT$DISC.DIR, fixed=T), fixed=T),
               recursive=T, showWarnings=F)
}


# Process sample indicated by index
sample.id = sample.info$SAMPLE_NAME[SAMPLE.IDX]
normal.id = sample.info$NORMAL_NAME[SAMPLE.IDX]
species = sample.info$SPECIES_NAME[SAMPLE.IDX]
project = sample.info$PROJECT_ID[SAMPLE.IDX]
ref.name = sample.info$REFERENCE_GENOME[SAMPLE.IDX]
cat("\nProcessing sample", sample.id, "\n\n")

# Find matched samples (from same individual)
idx = match(c("NORMAL_NAME", "SPECIES_NAME"), colnames(sample.match))
matched.samples = as.character(sample.match[sample.id, -idx])
matched.samples = matched.samples[file.exists(paste0(INPUT$CNT.PREFIX,
                                                     sample.id, "_In_", matched.samples, ".txt"))]

# Build input file paths
vars.path = gsub("${AC36}", ac36.path,
                 gsub("${SAMPLE}", sample.id,
                      gsub("${SPECIES}", species,
                           INPUT$VARIANTS, fixed=T), fixed=T), fixed=T)
filt.paths = c(gsub("${AC36}", ac36.path,
                    gsub("${SAMPLE}", sample.id,
                         gsub("${NORMAL}", normal.id,
                              gsub("${SPECIES}", species,
                                   INPUT$CAVEMAN.BED, fixed=T), fixed=T), fixed=T), fixed=T),
               gsub("${SAMPLE}", sample.id, INPUT$MIN.COV.BED, fixed=T),
               gsub("${SAMPLE}", sample.id, INPUT$MAX.COV.BED, fixed=T),
               gsub("${SPECIES}", species, INPUT$N50.BED, fixed=T),
               gsub("${AC36}", ac36.path,
                    gsub("${SPECIES}", species, INPUT$ENDS.BED, fixed=T), fixed=T))
cnt.paths = paste0(INPUT$CNT.PREFIX, sample.id, "_In_", matched.samples, ".txt")
spcs = ifelse(species %in% bam.paths$SPECIES, species, "default")
bam.path = gsub("${SPECIES}", species,
                gsub("${REFGENOME}", ref.name, 
                     gsub("${PROJECT}", project,
                          gsub("${SAMPLE}", sample.id,
                               bam.paths$PATH[bam.paths$SPECIES == spcs],
                               fixed=T), fixed=T), fixed=T), fixed=T)
all.paths = c(vars.path, filt.paths, cnt.paths, bam.path)

if (!all(file.exists(all.paths))) {
    cat("WARNING: File(s)", paste(all.paths[!file.exists(all.paths)], collapse=", "),
        "not found. Skipping sample.\n")
} else {
    
    
    # Read files of variants and filter regions
    variants = read.table(vars.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
    variants = variants[order(variants[, 1], variants[, 2]), ]
    cat(nrow(variants), "variants read from input file\n")
    
    if (nrow(variants) == 0) {
        cat("WARNING: No variants found in input file. Skipping sample.\n")
    } else {
        
        variants.gr = makeGRangesFromDataFrame(variants, seqnames.field="Chr",
                                               start.field="Start", end.field="Start")
        
        filters.gr = lapply(filt.paths, function(fpath) {
            makeGRangesFromDataFrame(read.table(fpath, sep="\t", as.is=T,
                                                col.names=c("chrom", "start", "end")),
                                     starts.in.df.are.0based=T)
        })
        
        
        # (1-5) Identify variants affected by each region-based filter
        discarded.idx = sapply(filters.gr, function(filt.gr) {
            # overlapsAny finds the ranges in query that overlap any of the ranges in subject,
            # and returns a logical vector of length equal to the number of ranges in query
            overlapsAny(variants.gr, filt.gr)
        })
        if (is.null(dim(discarded.idx))) {
            discarded.idx = matrix(discarded.idx, nrow=1)
        }
        
        # Negate index for contig-ends filter, as this defines regions NOT to exclude
        j = match("contigends", FILTER.NAMES)
        discarded.idx[, j] = !discarded.idx[, j]
        
        
        # (6) For chrom-level assemblies, identify variants outside chromosomal contigs
        if (species %in% chrom.lists$Species) {
            not.chrom.idx = !(variants$Chr %in%
                                  chrom.lists$Chromosome[chrom.lists$Species == species])
        } else {
            not.chrom.idx = rep(FALSE, nrow(variants))
        }
        
        
        # (7) Identify clustered variants (below min. distance to the previous/next variant)
        next.chrom.idx = variants[-nrow(variants), 1] == variants[-1, 1]
        next.dist = diff(variants[, 2])
        near.prev.idx = c(FALSE, next.chrom.idx & next.dist < MIN.DIST)
        near.next.idx = c(next.chrom.idx & next.dist < MIN.DIST, FALSE)
        clustered.idx = near.prev.idx | near.next.idx
        
        
        # (8) Apply beta-binomial model to identify artefactual variants based
        # on their distribution across samples from the same individual
        var.nr = var.nv = NULL
        var.str = paste0(variants$Chr, ":", variants$Start)
        for (fpath in cnt.paths) {
            counts = read.table(fpath, header=T, comment.char="", check.names=F, as.is=T)
            stopifnot(identical(paste0(counts$`#CHR`, ":", counts$POS), var.str))
            var.nr = cbind(var.nr, counts$Good_depth)
            var.nv = cbind(var.nv,
                           sapply(1:nrow(variants), function(j) {
                               counts[j, paste0("Count_", variants$Alt[j])]
                           }))
        }
        # Obtain MLE of each variant's "overdispersion" (rho) for a beta-binomial;
        # shared variants with a low value of rho are likely to be artefacts
        # (adapted from Tim Coorens)
        var.nr[var.nr == 0] = 1
        shared.idx = rowSums(var.nv > 0) > 1
        rho.mle = sapply(1:nrow(var.nv), function(j) {
            rho = 10 ^ seq(-6, -0.05, by=0.05)  # rho bounded within 1e-6 and 0.89
            nr = as.numeric(var.nr[j, ])
            nv = as.numeric(var.nv[j, ])
            mu = sum(nv) / sum(nr)
            loglik = sapply(rho, function(r) {
                sum(VGAM::dbetabinom(x=nv, size=nr, prob=mu, rho=r, log=T))
            })
            rho[which.max(loglik)]
        })
        low.rho.idx = shared.idx & rho.mle < MIN.RHO
        
        
        # (9) Identify variants below minimum VAF threshold, defined as
        #     half the median VAF of variants passing the beta-binomial filter
        idx = matched.samples == sample.id
        vaf = var.nv[, idx] / var.nr[, idx]
        stopifnot(!any(is.na(vaf)) | length(vaf) == length(rho.mle))
        if (any(!low.rho.idx)) {
            low.vaf.idx = vaf < median(vaf[!low.rho.idx]) / 2
        } else {
            low.vaf.idx = rep(FALSE, nrow(var.nv))
        }
        
        
        # (10-11) Identify variants with strand bias and variants near indels
        strand.bias.idx = rep(FALSE, nrow(variants))
        indel.idx = rep(FALSE, nrow(variants))
        for (i in 1:nrow(variants)) {
            bases = bam2R(bam.path,
                          variants$Chr[i], variants$Start[i] - FLANK, variants$Start[i] + FLANK,
                          q=30, mask=3844, mq=10)
            # Count number of supporting reads on each strand, and check if any is <1
            nv.strand = bases[FLANK + 1, c(variants$Alt[i], tolower(variants$Alt[i]))]
            if (any(nv.strand < MIN.NV.STRAND)) {
                strand.bias.idx[i] = TRUE
            }
            # Count total support (bases*reads) for an indel around the variant, and check
            # if this is >INDEL.FACTOR times the number of reads supporting the variant
            nv.indels = sum(bases[, c("-", "INS", "DEL", "_", "ins", "del")])
            if (nv.indels > INDEL.FACTOR * sum(nv.strand)) {
                indel.idx[i] = TRUE
            }
        }
        
        
        # Combine filters
        discarded.idx = cbind(discarded.idx, not.chrom.idx, clustered.idx,
                              low.rho.idx, low.vaf.idx, strand.bias.idx, indel.idx)
        stopifnot(ncol(discarded.idx) == length(FILTER.NAMES))
        discarded.any.idx = rowSums(discarded.idx) > 0
        discarded.counts = colSums(discarded.idx)
        discarded.counts.excl = sapply(1:ncol(discarded.idx), function(j) {
            sum(discarded.idx[, j] & rowSums(discarded.idx) == 1)
        })
        
        
        # Write variants affected by each filter to separate files
        for (j in 1:ncol(discarded.idx)) {
            out.path = paste0(gsub("${AC36}", ac36.path,
                                   gsub("${SPECIES}", species, OUTPUT$DISC.DIR, fixed=T), fixed=T),
                              "/", gsub("${SAMPLE}", sample.id,
                                        gsub("${FILTER}", FILTER.NAMES[j],
                                             OUTPUT$DISC.FILE, fixed=T), fixed=T))
            #cat(discarded.counts[j], " variants discarded by filter '", FILTER.NAMES[j], "' written to ", out.path, "\n", sep="")
            cat(discarded.counts[j], " variants discarded by filter '", FILTER.NAMES[j], "' (",
                discarded.counts.excl[j], " exclusively)\n", sep="")
            write.table(variants[discarded.idx[, j], ], file=out.path,
                        quote=F, sep="\t", row.names=F, na="")
        }
        cat(sum(discarded.any.idx), "variants discarded in total\n")
        
        # Write filtered variants to output file
        out.path = paste0(gsub("${AC36}", ac36.path, 
                               gsub("${SPECIES}", species, OUTPUT$FINAL.DIR, fixed=T), fixed=T),
                          "/", gsub("${SAMPLE}", sample.id, OUTPUT$FINAL.FILE, fixed=T))
        #cat(sum(!discarded.any.idx), "filtered variants written to", out.path, "\n")
        cat(sum(!discarded.any.idx), "variants in final variant set\n")
        write.table(variants[!discarded.any.idx, ], file=out.path,
                    quote=F, sep="\t", row.names=F, na="")
        
        
        # Define and output callable genome regions (from filters 1-6)
        # Callable genome is defined by subtracting all 'negative' region
        # filters from the 'positive' contig-ends filter
        idx = match("contigends", FILTER.NAMES)
        callable.genome = filters.gr[[idx]]
        for (i in (1:length(filters.gr))[-idx]) {
            callable.genome = setdiff(callable.genome, filters.gr[[i]])
        }
        # For species with chrom-level assemblies, remove non-chromosomal regions
        if (species %in% chrom.lists$Species) {
            chroms.gr = makeGRangesFromDataFrame(chrom.lists[chrom.lists$Species == species, 1:3])
            callable.genome = intersect(callable.genome, chroms.gr)
        }
        write.table(as.data.frame(callable.genome)[, 1:4],
                    col.names=c("chrom", "start", "end", "width"), sep="\t", row.names=F, quote=F,
                    file=paste0(OUTPUT$CALLABLE.DIR, "/", OUTPUT$CALLABLE.PREFIX, sample.id, ".txt"))
        
        
        # Save rho values, shared variant index and exclusive counts per filter
        save(shared.idx, rho.mle, discarded.counts.excl,
             file=paste0(OUTPUT$RHO.DIR, "/", OUTPUT$RHO.PREFIX, sample.id, ".RData"))
    }
}


cat("\nDone\n")
