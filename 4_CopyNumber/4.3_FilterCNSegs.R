# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 4.3: FILTER COPY NUMBER SEGMENTS 

# For each species, filter copy number segments different from the default
# CN=2 (1-1) which are shared among multiple samples (as we do not expect shared CNVs)


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CN.PREFIX = "../data/processed/CN_AlleleSpecific/CN_AlleleSpecific_"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/CN_AlleleSpecific/CN_AlleleSpecific_All_FilteredSegs.RData",
    CN.DIR = paste0("../output/", Sys.Date(), "_CN_AlleleSpecific_Filtered"),
    CN.A.PREFIX = paste0(Sys.Date(), "_CNSegments_Allele_Filt_"),
    CN.T.PREFIX = paste0(Sys.Date(), "_CNSegments_Total_Filt_"),
    PLOT.PREFIX = paste0(Sys.Date(), "_CNPlots_Filtered_")
)


# Samples to exclude from CN analysis (bad SNP quality)
EXCLUDE.LIST = c("PD36813x15", "PD36813x16")

# Species to consider
SPECIES = c("cat", "cow", "dog", "horse", "human", "mouse", "rabbit", "rat")

# Minimum length for segments of CN≠2 (1-1)
MIN.LEN = 1e6

# Disable scientific notation
options(scipen=999)


# Create output directories
dir.create(OUTPUT$CN.DIR, showWarnings=F)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(GenomicRanges))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Process CN segments for each species
cn.segments.orig = cn.segments.filt = sapply(SPECIES, function(x) NULL, simplify=F)
for (species in SPECIES) {
    
    # Load CN segments in each sample
    cat("\nProcessing species:", species, "\n")
    cat("Loading CN segments\n")
    sample.ids = sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species]
    sample.ids = sample.ids[!(sample.ids %in% EXCLUDE.LIST) &
                                file.exists(paste0(INPUT$CN.PREFIX, sample.ids, ".RData"))]
    cn.segments.orig[[species]] = sapply(sample.ids, function(id) {
        load(paste0(INPUT$CN.PREFIX, id, ".RData"))
        cn.segments
    }, simplify=F)

    # Collect all segments without standard (1-1) CN
    cn.levels.total = unique(unlist(sapply(cn.segments.orig[[species]], function(cn) {
        cn$total$CN[cn$total$CN != 2]
    })))
    cn.levels.allele = unique(unlist(sapply(cn.segments.orig[[species]], function(cn) {
        cn$allele$CN[cn$allele$CN != "1-1"]
    })))
    
    cnvs = data.frame(Chrom="1", Start=1, End=1, stringsAsFactors=F)
    for (cn in cn.segments.orig[[species]]) {
        for (lvl in as.character(cn.levels.total)) {
            cnvs = rbind(cnvs,
                         cn$total[as.character(cn$total$CN) == lvl, c("Chrom","Start","End")])
        }
        for (lvl in cn.levels.allele) {
            cnvs = rbind(cnvs,
                         cn$allele[cn$allele$CN == lvl, c("Chrom","Start","End")])
        }
    }
    cnvs.gr = makeGRangesFromDataFrame(cnvs)
    
    
    # Segment filtering:
    # For each sample, discard segments of CN≠2 (1-1) and presenting
    # >2 overlaps with any segments of CN≠2 (1-1) (including itself).
    # Also discard segments of CN≠2 (1-1) with length < 1 Mb.
    cn.segments.filt[[species]] = cn.segments.orig[[species]]
    for (id in sample.ids) {
        cat("Processing sample", id, "\n")
        cn.sample = cn.segments.filt[[species]][[id]]$total
        keep.total.idx = sapply(1:nrow(cn.sample), function(i) {
            if (cn.sample$CN[i] == 2) {
                TRUE
            } else {
                cn.sample$Length[i] >= MIN.LEN &
                    countOverlaps(makeGRangesFromDataFrame(cn.sample[i, ]), cnvs.gr) < 3
            }
        })
        cn.sample = cn.segments.filt[[species]][[id]]$allele
        keep.allele.idx = sapply(1:nrow(cn.sample), function(i) {
            if (cn.sample$CN[i] == "1-1") {
                TRUE
            } else {
                cn.sample$Length[i] >= MIN.LEN &
                    countOverlaps(makeGRangesFromDataFrame(cn.sample[i, ]), cnvs.gr) < 3
            }
        })
        cn.segments.filt[[species]][[id]]$total = 
            cn.segments.filt[[species]][[id]]$total[keep.total.idx, ]
        cn.segments.filt[[species]][[id]]$allele =
            cn.segments.filt[[species]][[id]]$allele[keep.allele.idx, ]
    }
    
    # # ALTERNATIVE APPROACH: FILTER BASED ON SEGMENTS WITH SAME CN
    # cnvs.total = sapply(as.character(cn.levels.total), function(x) NULL, simplify=F)
    # cnvs.allele = sapply(cn.levels.allele, function(x) NULL, simplify=F)
    # for (cn in cn.segments.filt[[species]]) {
    #     for (lvl in as.character(cn.levels.total)) {
    #         cnvs.total[[lvl]] = rbind(cnvs.total[[lvl]],
    #                                   cn$total[as.character(cn$total$CN) == lvl, ])
    #     }
    #     for (lvl in cn.levels.allele) {
    #         cnvs.allele[[lvl]] = rbind(cnvs.allele[[lvl]],
    #                                    cn$allele[cn$allele$CN == lvl, ])
    #     }
    # }
    # cnvs.total.gr = sapply(cnvs.total, makeGRangesFromDataFrame)
    # cnvs.allele.gr = sapply(cnvs.allele, makeGRangesFromDataFrame)
    # 
    # # Segment filtering:
    # # For each sample, discard segments with CN≠2 (or 1-1) and presenting
    # # >1 overlap with segments of the same CN (including itself)
    # for (id in sample.ids) {
    #     cat("Processing sample", id, "\n")
    #     keep.total.idx = sapply(1:nrow(cn.segments.filt[[species]][[id]]$total), function(i) {
    #         cn = cn.segments.filt[[species]][[id]]$total$CN[i]
    #         if (cn == 2) {
    #             TRUE
    #         } else {
    #             countOverlaps(makeGRangesFromDataFrame(cn.segments.filt[[species]][[id]]$total[i,]),
    #                           cnvs.total.gr[[as.character(cn)]]) < 2
    #         }
    #     })
    #     keep.allele.idx = sapply(1:nrow(cn.segments.filt[[species]][[id]]$allele), function(i) {
    #         cn = cn.segments.filt[[species]][[id]]$allele$CN[i]
    #         if (cn == "1-1") {
    #             TRUE
    #         } else {
    #             countOverlaps(makeGRangesFromDataFrame(cn.segments.filt[[species]][[id]]$allele[i,]),
    #                           cnvs.allele.gr[[cn]]) < 2
    #         }
    #     })
    #     cn.segments.filt[[species]][[id]]$total = 
    #         cn.segments.filt[[species]][[id]]$total[keep.total.idx, ]
    #     cn.segments.filt[[species]][[id]]$allele =
    #         cn.segments.filt[[species]][[id]]$allele[keep.allele.idx, ]
    # }
}


# Plot CN segments
cat("\nPlotting and saving CN segments\n")
for (i in 1:length(cn.segments.filt)) {
    for (j in 1:length(cn.segments.filt[[i]])) {
        sample.id = names(cn.segments.filt[[i]])[j]
        normal.id = sample.info$NORMAL_NAME[match(sample.id, sample.info$SAMPLE_NAME)]
        load(paste0(INPUT$CN.PREFIX, sample.id, ".RData"))
        
        chroms = unique(bin.table$Chrom)
        cutoffs = c(match(chroms, bin.table$Chrom), nrow(bin.table))
        cols = rep("black", length(chroms))
        cols[seq_along(cols) %% 2 == 0] = "darkgrey"
        cols.s2 = rep(cols, table(bin.table$Chrom)[chroms])[set2.idx]
        cols = rep(cols, table(bin.table$Chrom)[chroms])[set1.idx]
        xlm = c(1, nrow(bin.table))
        offset = 0.1; cxlab = 1.5; cxmain = 1.7
        
        png(paste0(OUTPUT$CN.DIR, "/", OUTPUT$PLOT.PREFIX, sample.id, ".png"), 17, 8, "in", res=120)
        par(mfrow=c(3, 1), mar=c(0, 5, 3, 1.25), oma=c(2.5, 0, 3, 0))
        
        # Plot obs/exp coverage ratio
        plot(which(set1.idx), cov.obs / cov.exp,
             col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="Obs/Exp coverage ratio", 
             xpd=NA, xlim=xlm, xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main="Observed/Expected bin coverage ratio", cex.lab=cxlab, cex.main=cxmain)
        points(which(set2.idx), cov.obs.s2 / cov.exp.s2, col=cols.s2, pch=16, cex=0.6)
        title(paste0(sample.id, " (normal ", normal.id, "; ",
                     prettyNum(sum(set1.idx), big.mark=","), " bins with SNPs)"),
              cex.main=2, line=3.5, xpd=NA)
        
        # Plot het BAF
        plot(which(set1.idx), vaf.obs,
             col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="BAF",
             xlim=xlm, ylim=c(0, 1), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main="Heterozygous SNP BAF", cex.lab=cxlab, cex.main=cxmain)
        
        # Plot CN segments
        plot(cn.table$CN_total, type="n", las=2, xlab="", ylab="Copy number",
             xlim=xlm, ylim=c(0, 4), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main="Copy number segments (filtered)", cex.lab=cxlab, cex.main=cxmain)
        segments(x0=cn.segments.filt[[i]][[j]]$total$StartBin,
                 x1=cn.segments.filt[[i]][[j]]$total$EndBin,
                 y0=cn.segments.filt[[i]][[j]]$total$CN, col="forestgreen", lwd=4, lend=2)
        segments(x0=cn.segments.filt[[i]][[j]]$allele$StartBin,
                 x1=cn.segments.filt[[i]][[j]]$allele$EndBin,
                 y0=cn.segments.filt[[i]][[j]]$allele$nMajor + offset, col="red", lwd=4, lend=2)
        segments(x0=cn.segments.filt[[i]][[j]]$allele$StartBin,
                 x1=cn.segments.filt[[i]][[j]]$allele$EndBin,
                 y0=cn.segments.filt[[i]][[j]]$allele$nMinor - offset, col="blue", lwd=4, lend=2)
        
        mtext(chroms, side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.7, font=2)
        invisible(dev.off())
    }
}


# # Plot numbers of segments with CN≠2 per sample before and after filtering
# cnv.count.orig = sapply(cn.segments.orig, function(cn.spc) {
#     sapply(cn.spc, function(cn.smp) {
#         length(table(cn.smp$total$Chrom[cn.smp$total$CN != 2]))
#         #sum(cn.smp$total$CN != 2)  #+ sum(cn.smp$allele$CN != "1-1")
#     })
# })
# cnv.count.filt = sapply(cn.segments.filt, function(cn.spc) {
#     sapply(cn.spc, function(cn.smp) {
#         length(table(cn.smp$total$Chrom[cn.smp$total$CN != 2]))
#         #sum(cn.smp$total$CN != 2)  #+ sum(cn.smp$allele$CN != "1-1")
#     })
# })
# plot(1, type="n", xlim=c(0.5, length(cnv.count.orig)+0.5), ylim=c(0, max(unlist(cnv.count.orig))),
#      xlab="", ylab="Chromosomes with segments of CN ≠ 2", xaxt="n")
# axis(1, 1:length(cnv.count.orig), names(cnv.count.orig))
# for (i in 1:length(cnv.count.orig)) {
#     x0 = rep(i, length(cnv.count.orig[[i]])) - runif(length(cnv.count.orig[[i]]), 0.1, 0.4)
#     x1 = rep(i, length(cnv.count.filt[[i]])) + runif(length(cnv.count.filt[[i]]), 0.1, 0.4)
#     segments(x0=x0, x1=x1,
#              y0=cnv.count.orig[[i]], y1=cnv.count.filt[[i]], col="darkgrey")
#     points(x0, cnv.count.orig[[i]], col="firebrick", pch=16)
#     points(x1, cnv.count.filt[[i]], col="dodgerblue", pch=16)
# }

    
# Output and save CN segments
for (i in 1:length(cn.segments.filt)) {
    for (j in 1:length(cn.segments.filt[[i]])) {
        sample.id = names(cn.segments.filt[[i]])[j]
        write.table(cn.segments.filt[[i]][[j]]$total, sep="\t", quote=F, row.names=F, col.names=T,
                    file=paste0(OUTPUT$CN.DIR, "/", OUTPUT$CN.T.PREFIX, sample.id, ".txt"))
        write.table(cn.segments.filt[[i]][[j]]$allele, sep="\t", quote=F, row.names=F, col.names=T,
                    file=paste0(OUTPUT$CN.DIR, "/", OUTPUT$CN.A.PREFIX, sample.id, ".txt"))
    }
}

save(cn.segments.filt, cn.segments.orig, file=OUTPUT$DATA)

cat("\nDone\n")
