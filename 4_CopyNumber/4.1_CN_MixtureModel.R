# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 4.1: INFER TOTAL COPY NUMBER USING A NORMAL MIXTURE MODEL


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    COV.PREFIX = "../data/processed/BinCoverage/BinCov_"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/CN_PerSpecies.RData",
    PLOT.PREFIX = paste0("../output/", Sys.Date(), "_CopyNumber_")
)


# Species to consider
SPECIES = c("cat", "cow", "dog", "horse", "human", "mouse", "rabbit", "rat")

# Read length (bp; for coverage calculation)
READ.LEN = 150

# Parameters of normal mixture model
PURITY = 0.85
MIN.PROB = 0.9
CN.STATES = 0:8
CN.COVR = PURITY*CN.STATES/2 + (1-PURITY)
SEED = 0xC0FFEE
MIN.SEGMENT = 5


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(GenomicRanges))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(length(unique(sample.info$SAMPLE_NAME)) == nrow(sample.info))
cat("Loaded\n")


cn.per.species = lapply(SPECIES, function(species) {
    #species = "mouse"
    cat("\nProcessing species:", species, "\n")
    
    pdf(paste0(OUTPUT$PLOT.PREFIX, species, ".pdf"), 16, 7)
    par(mar=c(3, 5, 4, 4))
    
    # Load bin coverage data
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    samples = unique(c(sample.info$SAMPLE_NAME[species.idx], sample.info$NORMAL_NAME[species.idx]))
    
    cat("Loading bin coverage data\n")
    bin.coverage = lapply(samples, function(sample.id) {
        cov = read.table(paste0(INPUT$COV.PREFIX, sample.id, ".txt"), sep="\t",
                         col.names=c("Chrom", "Start", "End", "Cov", "", "Length", ""), as.is=T)
        cov$Start = cov$Start + 1
        cov$AdjCov = cov$Cov / cov$Length * READ.LEN
        if (species == "cat") {
            cov.order = order(cov$Chrom, cov$Start)
        }
        else {
            cov.order = order(as.numeric(gsub("X", 100, cov$Chrom)), cov$Start)
        }
        cov[cov.order, c("Chrom", "Start", "End", "Length", "Cov", "AdjCov")]
    })
    names(bin.coverage) = samples
    
    # Call CN in each sample
    sample.cn = lapply(sample.info$SAMPLE_NAME[species.idx], function(sample.id) {
        
        #sample.id = "MD6260ab_lo0004"; normal.id = "MD6260z"
        cat("Processing sample", sample.id, "\n")
        normal.id = sample.info$NORMAL_NAME[match(sample.id, sample.info$SAMPLE_NAME)]
        
        # Select bin coverage data
        bin.idx = (bin.coverage[[1]]$Chrom != "Y") & (bin.coverage[[normal.id]]$AdjCov > 0)
        bin.coords = bin.coverage[[1]][bin.idx, 1:4]
                
        # Calculate normalised coverage ratio
        covR = bin.coverage[[sample.id]]$AdjCov[bin.idx] / bin.coverage[[normal.id]]$AdjCov[bin.idx]
        covR = covR - median(covR, na.rm=T) + 1
        covR[covR < 0] = 0
        if (length(unique(covR)) == 1) {
            covR[1] = covR[1] * 1.1
        }
        #summary(covR); plot(seq_along(covR), covR); plot(density(covR, na.rm=T))
        
        set.seed(SEED)
        mix = normalmixEM(covR, k=length(CN.COVR),
                          maxit=5000, epsilon=1e-10, mean.constr=CN.COVR, arbvar=F)
        
        # Calculate integer CN from mixture model
        bin.cn = apply(mix$posterior, 1, function(x) {
            cn = which(x >= MIN.PROB) - 1
            if (length(cn) == 0)
                NA
            else
                cn
        })
        
        # # Remove CNV segments below minimum length
        # cnv.bins = which((bin.cn != 2) %in% T)
        # if (length(cnv.bins) > 0) {
        #     start.bins = cnv.bins[c(TRUE, diff(cnv.bins) > 1)]
        #     for (i in start.bins) {
        #         idx = i:min(i+MIN.SEGMENT-1, length(bin.cn))
        #         if (!all((bin.cn[idx] == bin.cn[i]) %in% T)) {
        #             bin.cn[idx][bin.cn[idx] == bin.cn[i]] = NA
        #         }
        #     }
        # }
        
        # Infer CN segments above minimum length
        cn.segs.prelim = NULL
        start.idx = 1
        while (start.idx < length(bin.cn)) {
            chrom = bin.coords$Chrom[start.idx]
            chrom.end = rev(which(bin.coords$Chrom == chrom))[1]
            cn = bin.cn[start.idx]
            if (is.na(cn)) {
                start.idx = start.idx + 1
            }
            else {
                # Bins with CN=NA are not counted as ends for segments with CN=2
                if (cn == 2) {
                    end.idx = start.idx + which((bin.cn[-(1:start.idx)] != cn) %in% T) - 1
                }
                else {
                    end.idx = start.idx - 1 + 
                        which(((bin.cn[-(1:start.idx)] != cn) %in% T) |
                                  is.na(bin.cn[-(1:start.idx)]))
                }
                if (length(end.idx) == 0 | end.idx[1] > chrom.end) {
                    end.idx = chrom.end
                }
                else {
                    end.idx = end.idx[1]
                }
                if (end.idx - start.idx + 1 >= MIN.SEGMENT) {
                    cn.segs.prelim = rbind(cn.segs.prelim,
                                           data.frame(StartBin=start.idx, EndBin=end.idx,
                                                      Chrom=chrom,
                                                      Start=bin.coords$Start[start.idx],
                                                      End=bin.coords$End[end.idx],
                                                      Length=bin.coords$End[end.idx] -
                                                          bin.coords$Start[start.idx] + 1,
                                                      CN=cn,
                                                      stringsAsFactors=F))
                }
                start.idx = end.idx + 1
            }
        }
        
        # Merge adjacent segments with the same CN,
        # and separated by a gap below the minimum segment length
        cn.segs = NULL
        i = 1
        while (i <= nrow(cn.segs.prelim)) {
            segment = cn.segs.prelim[i, ]
            j = i + 1
            while (j <= nrow(cn.segs.prelim) &
                   cn.segs.prelim$Chrom[j] == segment$Chrom &
                   cn.segs.prelim$CN[j] == segment$CN &
                   cn.segs.prelim$StartBin[j] - segment$EndBin <= MIN.SEGMENT) {
                
                segment$EndBin = cn.segs.prelim$EndBin[j]
                segment$End = cn.segs.prelim$End[j]
                segment$Length = segment$End - segment$Start + 1
                j = j + 1
                
            }
            cn.segs = rbind(cn.segs, segment)
            i = j
        }
        rownames(cn.segs) = NULL
        
        # Plot and return CN values
        chroms = unique(bin.coords$Chrom)
        cutoffs = c(match(chroms, bin.coords$Chrom), nrow(bin.coords))
        cols = rep("black", length(chroms))
        cols[seq_along(cols) %% 2 == 0] = "darkgrey"
        cols = rep(cols, table(bin.coords$Chrom)[chroms])
        plot(seq_along(covR), covR,
             col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="Normalised read coverage ratio",
             ylim=c(0, min(10, max(max(covR), CN.COVR[max(covR) < CN.COVR][1], CN.COVR[6], na.rm=T))),
             xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main=paste0(sample.id, " (", species, ")\n(Matched normal: ", normal.id, ")"))
        #points(seq_along(covR), CN.COVR[bin.cn + 1], col="red", pch=15, cex=0.45)
        segments(x0=cn.segs$StartBin, x1=cn.segs$EndBin, y0=CN.COVR[cn.segs$CN+1], col="red", lwd=5)
        axis(4, at=CN.COVR, labels=paste0("CN", CN.STATES), col.axis="red", las=2, cex.axis=0.9)
        mtext(chroms, side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.7)
        
        list(CN=cn.segs, bins=cbind(bin.coords, CovR=covR, CN=bin.cn), model=mix)
    })
    
    dev.off()
    
    names(sample.cn) = sample.info$SAMPLE_NAME[species.idx]
    sample.cn
})

names(cn.per.species) = SPECIES


save(cn.per.species, file=OUTPUT$DATA)
cat("\nDone\n")
