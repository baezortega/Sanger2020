# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 5.2: PLOT GENOME COVERAGE HISTOGRAMS FOR EACH SAMPLE


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    COV.PREFIX = "../data/processed/CoverageHist/CoverageHist_"
)

# Output file paths
OUTPUT = list(
    COV.PCT = "../data/processed/Coverage_99Percentiles.RData",
    PDF.DIR = paste0("../output/", Sys.Date(), "_CoverageHist"),
    PDF.PREFIX = paste0(Sys.Date(), "_CoverageHist_")
)


cat("Loading data...\n")
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Create output directory
dir.create(OUTPUT$PDF.DIR, showWarnings=F)


# Initialise object for 99th percentiles
species.list = unique(sample.info$SPECIES_NAME)
coverage.99percentiles = structure(vector(mode="list", length(species.list)), names=species.list)

# For each species
for (species in species.list) {
    cat("\nProcessing species:", species, "\n")
    species.idx = sample.info$SPECIES_NAME == species
    
    # Plot coverage histograms for each sample
    pdf(paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.PREFIX, species, ".pdf"), 10, 6)
    par(mar=c(5, 5, 4, 2))
    
    samples = unique(c(sample.info$SAMPLE_NAME[species.idx],
                       sample.info$NORMAL_NAME[species.idx]))
    
    percentiles = sapply(samples, function(sample.id) {
        
        cat("Processing sample ", sample.id, "\n", sep="")
        if (file.exists(paste0(INPUT$COV.PREFIX, sample.id, ".txt"))) {
            
            cov = read.table(paste0(INPUT$COV.PREFIX, sample.id, ".txt"), sep="\t", as.is=T)
            cum.fraction = sapply(1:nrow(cov), function(i) sum(cov[i:nrow(cov), 5]))
            cov.median = cov[which(cum.fraction < 0.5)[1] - 1, 2]
            cov.99pctile = cov[which(cum.fraction < 0.01)[1] - 1, 2]
            cov.995pctile = cov[which(cum.fraction < 0.005)[1] - 1, 2]
            cov.75.perc = round(sum(cov[cov[, 2] >= 75, 5]) * 100, 2)
            cov.100.perc = round(sum(cov[cov[, 2] >= 100, 5]) * 100, 2)
            
            MAX.COV = 150
            bars = barplot(cov[1:(MAX.COV+1), 5], plot=F)
            barplot(cov[1:(MAX.COV+1), 5],
                    col="dodgerblue4", border="white", xlim=c(0, rev(bars)[1]),
                    xlab="Coverage", ylab="Genome fraction",
                    main=paste0("Coverage distribution in ", sample.id, " (", species, ")"))
            legend("topright",
                   legend=c(paste("Median =", cov.median), paste("99th percentile =", cov.99pctile),
                            paste("99.5th percentile =", cov.995pctile),
                            paste0("Bases with cov. 75+: ", cov.75.perc, "%    "),
                            paste0("Bases with cov. 100+: ", cov.100.perc, "%    ")))
            axis(1, at=bars[seq(1, MAX.COV+1, 25)], labels=cov[seq(1, MAX.COV+1, 25), 2])
            
            # Return 99th percentile
            cov.99pctile
        }
        
        else {
            cat("WARNING: Data not found for sample", sample.id, "\n")
            NA
        }
    })
    
    coverage.99percentiles[[species]] = percentiles
    invisible(dev.off())
}


# Save coverage percentiles
save(coverage.99percentiles, file=OUTPUT$COV.PCT)

cat("\nDone\n")
