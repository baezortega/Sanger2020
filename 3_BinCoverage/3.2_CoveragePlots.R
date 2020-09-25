# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 3.2: PLOT COVERAGE, LOGR AND VAF FOR EACH SAMPLE


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CHROM.LISTS = "../data/original/CrossSpecies_ChromLists.txt",
    COV.PREFIX = "../data/processed/BinCoverage/BinCov_",
    SNPS.PATH = "../data/original/Path_SampleSNPs.txt"
)

# Output file paths
OUTPUT = list(
    #PLOT.PREFIX = paste0("../output/Coverage_LogR_VAF_Plots/", Sys.Date(), "_Coverage_LogR_VAF__")
    PLOT.PREFIX = paste0("../output/", Sys.Date(), "_Coverage_LogR_VAF_")
)


# Read length (bp; for coverage calculation)
READ.LEN = 150

# Create output folder
#dir.create("../output/Coverage_LogR_VAF_Plots", showWarnings=F)

# Disable scientific notation
options(scipen=999)


cat("Loading data and packages...\n")
library(stringr)
suppressPackageStartupMessages(library(GenomicRanges))
snps.path = as.character(as.matrix(read.table(INPUT$SNPS.PATH)))
chrom.lists = read.table(INPUT$CHROM.LISTS, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(length(unique(sample.info$SAMPLE_NAME)) == nrow(sample.info))
cat("Loaded\n")


# For each species
for (species in unique(chrom.lists$Species)) {
    cat("\nProcessing species:", species, "\n")
    
    # Load bin coverage data
    species.idx = sample.info$SPECIES_NAME == species
    samples = unique(c(sample.info$SAMPLE_NAME[species.idx], sample.info$NORMAL_NAME[species.idx]))
    
    cat("Loading bin coverage data\n")
    bin.coverage = lapply(samples, function(sample.id) {
        cov = read.table(paste0(INPUT$COV.PREFIX, sample.id, ".txt"), sep="\t",
                         col.names=c("Chrom", "Start", "End", "Cov", "", "Length", ""), as.is=T)
        cov$Start = cov$Start + 1
        cov$AdjCov = cov$Cov / cov$Length * READ.LEN
        if (species == "cat") {
            cov.order = order(cov$Chrom, cov$Start)
        } else {
            cov.order = order(as.numeric(gsub("X", 100, cov$Chrom)), cov$Start)
        }
        cov[cov.order, c("Chrom", "Start", "End", "Length", "Cov", "AdjCov")]
    })
    names(bin.coverage) = samples
    
    cutoffs = c(match(unique(bin.coverage[[1]]$Chrom), bin.coverage[[1]]$Chrom),
                nrow(bin.coverage[[1]]))
    cols = rep("black", length(unique(bin.coverage[[1]]$Chrom)))
    cols[seq_along(cols) %% 2 == 0] = "darkgrey"
    cols = rep(cols, table(bin.coverage[[1]]$Chrom)[unique(bin.coverage[[1]]$Chrom)])
    
    # Plot coverage, covR and VAF per sample
    #png(paste0(OUTPUT$PLOT.PREFIX, species, "__", sample.id, ".png"), 20, 12, units="in", res=600)
    pdf(paste0(OUTPUT$PLOT.PREFIX, species, ".pdf"), 20, 12.5)
    par(mfrow=c(4, 1), mar=c(2, 5, 1, 2), oma=c(1, 0, 3.5, 0))
    for (sample.id in sample.info$SAMPLE_NAME[species.idx]) {
        
        cat("Processing sample:", sample.id, "\n")
        normal.id = sample.info$NORMAL_NAME[match(sample.id, sample.info$SAMPLE_NAME)]
        
        # If VCF exists, load SNPs; otherwise, skip loading and make blank VAF plots
        sample.snps.path = gsub("${SPECIES}", species,
                                gsub("${SAMPLE}", sample.id,
                                     gsub("${NORMAL}", normal.id,
                                          snps.path, fixed=T), fixed=T), fixed=T)
        
        if (file.exists(sample.snps.path)) {
            cat("   Loading germline SNPs and retrieving VAFs\n")
            snps = read.table(gzfile(gsub("${SPECIES}", species,
                                          gsub("${SAMPLE}", sample.id,
                                               gsub("${NORMAL}", normal.id,
                                                    snps.path, fixed=T), fixed=T), fixed=T)),
                              sep="\t", comment.char="#", as.is=T,
                              col.names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                          "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOUR"))
            
            # Retrieve VAFs and assign SNPs to genome bins
            MIN.COV.N = 15
            bin.cov.ranges = with(bin.coverage[[1]], GRanges(Chrom, IRanges(Start, End)))
            snp.ranges = with(snps, GRanges(CHROM, IRanges(POS, POS)))
            snp.pos = findOverlaps(snp.ranges, bin.cov.ranges, select="first")
            snp.vaf = data.frame(pos = snp.pos,
                                 vaf_T = as.numeric(str_split_fixed(snps$TUMOUR, ":", 10)[, 10]),
                                 vaf_N = as.numeric(str_split_fixed(snps$NORMAL, ":", 10)[, 10]),
                                 cov_N = apply(str_split_fixed(snps$NORMAL, ":", 10)[, 2:9], 1,
                                               function(x) sum(as.numeric(x))))
            snp.vaf = snp.vaf[!is.na(snp.pos) & snp.vaf$cov_N >= MIN.COV.N, ]
            
            # Isolate heterozygous SNPs
            HET.VAF.MIN = 0.35
            HET.VAF.MAX = 0.65
            snp.vaf.het = snp.vaf[snp.vaf$vaf_N > HET.VAF.MIN & snp.vaf$vaf_N < HET.VAF.MAX, ]
            
            # Subsample to first and last SNPs per bin
            idx = match(unique(snp.vaf$pos), snp.vaf$pos)
            idx = sort(c(idx, idx[-1] - 1, rev(idx)[1] + 1))
            snp.vaf = snp.vaf[idx, ]
            idx = match(unique(snp.vaf.het$pos), snp.vaf.het$pos)
            idx = sort(c(idx, idx[-1] - 1, rev(idx)[1] + 1))
            snp.vaf.het = snp.vaf.het[idx, ]
        }
        else {
            cat("   WARNING: Germline SNPs not found. Omitting VAF plots.\n")
            snps = NULL
        }
        
        # Plot coverage
        cat("   Plotting coverage and VAF (if available)\n")
        cn.levels = (0:5) * median(bin.coverage[[sample.id]]$AdjCov) / 2
        plot(bin.coverage[[sample.id]]$AdjCov, pch=16, col=cols, cex=0.5, cex.lab=1.6, cex.axis=1.3,
             xlim=c(500, nrow(bin.coverage[[1]])-500), ylim=c(0, rev(cn.levels)[1]),
             xaxt="n", xlab="", ylab="Coverage",
             panel.first={abline(v=cutoffs, col="grey"); abline(h=cn.levels, lty=3)})
        title(paste0(sample.id, " (", species, ")\n(Matched normal: ", normal.id, ")"),
              cex.main=2, outer=T, line=-0.25)
        mtext(unique(bin.coverage[[sample.id]]$Chrom),
              side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.75)
        mtext(seq_along(cn.levels)-1, side=4, at=cn.levels, las=1, line=0.75)
        
        # Plot coverage ratio
        #covR = bin.coverage[[sample.id]]$AdjCov / bin.coverage[[normal.id]]$AdjCov
        #covR[is.infinite(covR)] = NA
        logR = log2(bin.coverage[[sample.id]]$AdjCov / bin.coverage[[normal.id]]$AdjCov)
        logR[is.infinite(logR)] = NA
        logR = logR - mean(logR, na.rm=T)
        min.logR = floor(min(logR, na.rm=T))
        max.logR = ceiling(max(logR, na.rm=T))
        plot(logR, pch=16, col=cols, cex=0.5, cex.lab=1.6, cex.axis=1.3,
             xlim=c(500, nrow(bin.coverage[[1]])-500), ylim=c(min.logR, max.logR),
             xaxt="n", xlab="", ylab="Sample-Normal logR",
             panel.first={abline(v=cutoffs, col="grey"); abline(h=min.logR:max.logR, lty=3)})
        mtext(unique(bin.coverage[[sample.id]]$Chrom),
              side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.75)
        
        # Plot VAF if available
        if (!is.null(snps)) {
            plot(snp.vaf$pos, snp.vaf$vaf_T, pch=16, col=cols[snp.vaf$pos],
                 cex=0.5, cex.lab=1.6, cex.axis=1.3,
                 xlim=c(500, nrow(bin.coverage[[1]])-500), ylim=c(0, 1),
                 xaxt="n", xlab="", ylab="Germline VAF",
                 panel.first={abline(v=cutoffs, col="grey"); abline(h=c(0, 0.5, 1), lty=3)})
            mtext(unique(bin.coverage[[sample.id]]$Chrom),
                  side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.75)
            
            # Plot VAF for heterozygous SNPs
            plot(snp.vaf.het$pos, snp.vaf.het$vaf_T, pch=16, col=cols[snp.vaf.het$pos],
                 cex=0.5, cex.lab=1.6, cex.axis=1.3,
                 xlim=c(500, nrow(bin.coverage[[1]])-500), ylim=c(0, 1),
                 xaxt="n", xlab="", ylab="Germline VAF (heterozygous SNPs)",
                 panel.first={abline(v=cutoffs, col="grey"); abline(h=c(0, 0.5, 1), lty=3)})
            mtext(unique(bin.coverage[[sample.id]]$Chrom),
                  side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.75)
        }
        else {
            plot(NA, xlim=c(500, nrow(bin.coverage[[1]])-500), ylim=c(0, 1),
                 cex.lab=1.6, cex.axis=1.3, xaxt="n", xlab="", ylab="Germline VAF")
            legend("center", legend="Germline SNPs unavailable", bty="n", cex=1.6)
            plot(NA, xlim=c(500, nrow(bin.coverage[[1]])-500), ylim=c(0, 1), cex.lab=1.6,
                 cex.axis=1.3, xaxt="n", xlab="", ylab="Germline VAF (heterozygous SNPs)")
            legend("center", legend="Germline SNPs unavailable", bty="n", cex=1.6)
        }
    }
    
    invisible(dev.off())
    invisible(gc())
}

cat("\nDone\n")
