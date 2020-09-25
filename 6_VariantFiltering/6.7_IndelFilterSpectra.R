# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 6.7: PRODUCE MUTATIONAL SPECTRA OF INDELS DISCARDED BY EACH FILTER


# Input file paths
INPUT = list(
    #GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    #RHO.PREFIX = "../data/processed/BetaBinom_Rho/BetaBinom_Rho_",
    AC36.PATH = "../data/original/Path_ac36.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CNT.TSV = "../data/processed/AlleleCounts_Indiv_Indels/${SAMPLE}/${NORMAL}_${SAMPLE}_indel_vaf.tsv",
    VARS.ORIG = "../data/processed/AlleleCounts_Indiv_Indels/data/Indels_PASS_${SAMPLE}.txt",
    VARS.DISC = "${AC36}/discarded_indel_calls/${SPECIES}/${SAMPLE}_discarded_indel_calls_${FILTER}.txt",
    VARS.FINAL = "${AC36}/final_indel_calls/${SPECIES}/${SAMPLE}_final_indel_calls.txt"
)

# Output file paths
OUTPUT = list(
    TMP.DIR = "SigProfiler_tmp",
    PDF.DIR.1 = paste0("../output/", Sys.Date(), "_IndelFilters_Spectra_PerSample"),
    PDF.DIR.2 = paste0("../output/", Sys.Date(), "_IndelFilters_Spectra_PerSpecies"),
    PDF.PREFIX = paste0(Sys.Date(), "_IndelFilters_Spectra_"),
    DATA.DIR = "../data/processed/IndelSpectra",
    SPEC.PREFIX = "IndelSpectra_",
    COUNTS = paste0("../output/", Sys.Date(), "_IndelFilters_Counts_All.txt")
    #RHO.DIR = paste0("../output/", Sys.Date(), "_BetaBinomRho_Spectra"),
    #VAF.DIR = paste0("../output/", Sys.Date(), "_BetaBinomRho_VAFHist"),
    #RHO.PREFIX = paste0(Sys.Date(), "_BetaBinomRho_Spectra_"),
    #VAF.PREFIX = paste0(Sys.Date(), "_BetaBinomRho_VAFHist_"),
    #PDF.DIR.3 = paste0("../output/", Sys.Date(), "_IndelFilters_Spectra_PerFilter"),
    #COUNTS.EXCL = paste0("../output/", Sys.Date(), "_IndelFilters_Counts_Exclusive.txt")
)


# Filter names
FILTER.NAMES = c("caveman", "lowcov", "highcov", "N50plus", "contigends",
                 "notchrom", "clustered", "betabinom", "lowVAF", "strandbias")

# SigProfiler reference genome names
REF.NAMES = c("human"="GRCh37", "mouse"="mm10", "rat"="rn6", "dog"="dog")

# Rho thresholds for beta-binomial filter comparisons
#MIN.RHO = seq(0.05, 0.8, 0.05)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(SigProfilerMatrixGeneratorR))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# HUMAN SAMPLES WITH COLIBACTIN EXPOSURE ARE EXCLUDED FROM GLOBAL SPECTRUM
exclude.list = sample.info$SAMPLE_NAME[sample.info$NORMAL_NAME %in%
                                           c("PD37449b", "PD37512a4", "PD37513a21")]


# Create output directories
dir.create(OUTPUT$PDF.DIR.1, showWarnings=F)
dir.create(OUTPUT$PDF.DIR.2, showWarnings=F)
dir.create(OUTPUT$DATA.DIR, showWarnings=F)
#dir.create(OUTPUT$PDF.DIR.3, showWarnings=F)
#dir.create(OUTPUT$RHO.DIR, showWarnings=F)
#dir.create(OUTPUT$VAF.DIR, showWarnings=F)


# Initialise lists of variants per filter and per species
species.list = unique(sample.info$SPECIES_NAME)
#vars.per.species = structure(vector(mode="list", length=length(species.list)), names=species.list)
#vars.per.filter = structure(vector(mode="list", length=length(FILTER.NAMES)), names=FILTER.NAMES)

# Initialise tables of counts per filter and per rho thresold
counts.per.filter.excl = NULL
counts.per.filter = NULL
counts.per.rho = NULL
sample.names = NULL


for (species in species.list) {
    cat("\nProcessing species:", species, "\n")
    species.vars = NULL
    
    for (sample.id in sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species]) {
        cat("Processing sample", sample.id, "\n")
        
        # Create temporary folder for SigProfiler
        dir.create(OUTPUT$TMP.DIR, showWarnings=F)
        
        # Read variant files
        paths = c(gsub("${AC36}", ac36.path,
                       gsub("${SAMPLE}", sample.id,
                            gsub("${SPECIES}", species,
                                 INPUT$VARS.FINAL, fixed=T), fixed=T), fixed=T),
                  gsub("${SAMPLE}", sample.id, 
                       gsub("${SPECIES}", species,
                            INPUT$VARS.ORIG, fixed=T), fixed=T),
                  sapply(FILTER.NAMES, function(name) {
                      gsub("${AC36}", ac36.path,
                           gsub("${SAMPLE}", sample.id,
                                gsub("${SPECIES}", species,
                                     gsub("${FILTER}", name,
                                          INPUT$VARS.DISC, fixed=T), fixed=T), fixed=T), fixed=T)
                  }))
        if (!all(file.exists(paths))) {
            cat("WARNING: File(s)", paste(paths[!file.exists(paths)], collapse=", "),
                "not found. Skipping sample.\n")
            next
        }
        
        vcf.names = paste0(sample.id, "_", LETTERS[1:(length(FILTER.NAMES)+2)],
                           c(":Final", ":Original", paste0(":Discarded_", FILTER.NAMES)))
        spectrum.names = paste(c("Final (filtered) variants",
                                 "Original variants",
                                 paste0("Variants discarded by filter '", FILTER.NAMES, "'")),
                               "in sample", sample.id)
        
        counts = rep(NA, length(paths))
        for (i in 1:length(paths)) {
            if (i == 2) {
                vars = read.table(paths[i], sep="\t", header=F, as.is=T,
                                  col.names=c("Chrom", "Pos", "Ref", "Alt"),
                                  colClasses=c("Chrom"="character", "Ref"="character", "Alt"="character"))
            } else {
                vars = read.table(paths[i], sep="\t", header=T, check.names=F, as.is=T,
                                  colClasses=c("Chrom"="character", "Ref"="character", "Alt"="character"))
            }
            counts[i] = nrow(vars)
            if (nrow(vars) > 0) {
                # Add final variants to species table
                # (except for human colibactin-exposed samples)
                if (!(sample.id %in% exclude.list) & (i == 1)) {
                    species.vars = rbind(species.vars, vars[, 1:4])
                }
                # Output indels to VCF
                if (species %in% names(REF.NAMES)) {
                    write.table(cbind(vars[, 1:2], ".", vars[, 3:4]),
                                file=paste0(OUTPUT$TMP.DIR, "/", vcf.names[i], ".vcf"),
                                sep="\t", quote=F, col.names=F, row.names=F)
                }
            }
        }
        
        # Plot spectra for sample
        if (species %in% names(REF.NAMES)) {
            spc = SigProfilerMatrixGeneratorR(sample.id, REF.NAMES[species], OUTPUT$TMP.DIR,
                                              plot=T, exome=F, chrom_based=F, tsb_stat=F, seqInfo=F)
            # Copy spectra to data and output folders
            file.copy(paste0(OUTPUT$TMP.DIR, "/output/plots/ID_83_plots_", sample.id, ".pdf"),
                      paste0(OUTPUT$PDF.DIR.1, "/", OUTPUT$PDF.PREFIX, sample.id, ".pdf"),
                      overwrite=T)
            file.copy(paste0(OUTPUT$TMP.DIR, "/output/ID/", sample.id, ".ID83.all"),
                      paste0(OUTPUT$DATA.DIR, "/", OUTPUT$SPEC.PREFIX, sample.id, ".txt"),
                      overwrite=T)
            # Delete SigProfiler folder
            unlink(OUTPUT$TMP.DIR, recursive=T)
        }
        
        # Save variant tables and counts
        #vars = variants
        #vars[, 1] = gsub(paste("sample", sample.id), species, vars[, 1])
        #vars.per.species[[species]] = rbind(vars.per.species[[species]], vars)
        #for (filt in FILTER.NAMES) {
        #    vars.per.filter[[filt]] = rbind(vars.per.filter[[filt]], vars[grep(filt, vars[, 1]), ])
        #}
        counts.per.filter = cbind(counts.per.filter, counts)
        sample.names = c(sample.names, sample.id)
        invisible(gc())
        
        # # Produce spectra of variants passing/failing the beta-binomial filter
        # # for a range of rho thresholds
        # #if (sum(grepl("Original", variants[, 1])) > 1 & sum(grepl("Mathij", variants[, 1])) > 1) {
        # if (sum(grepl("Original", variants[, 1])) > 1) {
        #             
        #     # Select original and filtered variants
        #     vars.orig = variants[grepl("Original", variants[, 1]), ]
        #     #vars.shared = variants[grepl("Mathij", variants[, 1]), ]
        #     #vars.shared[, 1] = "Variants discarded by Mathij's shared var. filter (disabled)"
        #     
        #     # Load rho estimates
        #     load(paste0(INPUT$RHO.PREFIX, sample.id, ".RData"))
        #     stopifnot(length(rho.mle) == nrow(vars.orig) & length(rho.mle) == length(shared.idx))
        #     
        #     # For each threshold of rho, plot spectra of passing and failing variants
        #     cairo_pdf(paste0(OUTPUT$RHO.DIR, "/", OUTPUT$RHO.PREFIX, sample.id, ".pdf"),
        #               width=18, height=11, onefile=T)
        #     par(mar=c(3, 3.6, 8, 2), oma=c(2, 0, 2, 0), mfrow=c(2, 1))
        #     #cairo_pdf(paste0(OUTPUT$RHO.DIR, "/", OUTPUT$RHO.PREFIX, sample.id, ".pdf"),
        #     #          width=18, height=13, onefile=T)
        #     #par(mar=c(3, 3.6, 8, 2), oma=c(2, 0, 2, 0), mfrow=c(3, 1))
        # 
        #     counts = NULL
        #     for (min.rho in MIN.RHO) {
        #         low.rho.idx = shared.idx & rho.mle < min.rho
        #         counts = c(counts, sum(low.rho.idx))
        #         vars = rbind(cbind(Sample=paste(sample.id,
        #                                         "\nVariants passing beta-binomial filter with min_rho =",
        #                                         min.rho),
        #                            vars.orig[!low.rho.idx, -1, drop=F]),
        #                      cbind(Sample=paste("Variants failing beta-binomial filter with min_rho =",
        #                                         min.rho),
        #                            vars.orig[low.rho.idx, -1, drop=F]))
        #                      #vars.shared)
        #         spectra = build_catalogues(vars)
        #         plot_spectrum(spectra)
        #         if (nrow(spectra) == 1) {
        #             plot(1, type="n", axes=F, xlab="", ylab="",
        #                  main=paste("No variants failing beta-binomial filter with min_rho =", min.rho))
        #         }
        #     }
        #     invisible(dev.off())
        #     
        #     # Save variant counts
        #     counts.per.rho = cbind(counts.per.rho, counts)
        #     counts.per.filter.excl = cbind(counts.per.filter.excl, discarded.counts.excl)
        #     invisible(gc())
        # }
        # else {
        #     counts.per.rho = cbind(counts.per.rho, rep(0, length(MIN.RHO)))
        #     counts.per.filter.excl = cbind(counts.per.filter.excl, rep(0, length(FILTER.NAMES)))
        # }
        # 
        # # Produce VAF histograms of variants passing/failing the beta-binomial filter
        # # for a range of rho thresholds
        # #if (sum(grepl("Original", variants[, 1])) > 1 & sum(grepl("Mathij", variants[, 1])) > 1) {
        # if (sum(grepl("Original", variants[, 1])) > 1) {
        #     
        #     # Calculate VAFs of variants in current sample
        #     counts = read.table(paste0(INPUT$CNT.PREFIX, sample.id, "_In_", sample.id, ".txt"),
        #                         header=T, comment.char="", check.names=F, as.is=T)
        #     stopifnot(nrow(counts) == nrow(vars.orig))
        #     var.nr = counts$Good_depth
        #     var.nv = sapply(1:nrow(vars.orig), function(j) {
        #         counts[j, paste0("Count_", vars.orig[j, "Alt"])]
        #     })
        #     vaf = var.nv / var.nr
        #     stopifnot(!any(var.nr == 0) & !any(is.na(vaf)) & length(vaf) == length(rho.mle))
        # 
        #     # Plot VAF histograms for each value of rho
        #     cairo_pdf(paste0(OUTPUT$VAF.DIR, "/", OUTPUT$VAF.PREFIX, sample.id, ".pdf"),
        #               width=12, height=8, onefile=T)
        #     par(mfrow=c(2, 1), mar=c(1, 4.5, 4, 0.25), oma=c(3.5, 0, 0.75, 0), mgp=c(2.2, 0.75, 0))
        #     for (min.rho in MIN.RHO) {
        #         low.rho.idx = shared.idx & rho.mle < min.rho
        #         hist(vaf[!low.rho.idx],
        #              breaks=seq(0, 1, 0.01), col="dodgerblue4", border="white", xlab="",
        #              main=paste0(sample.id, "\n", prettyNum(sum(!low.rho.idx), big.mark=","),
        #                          " variants passing beta-binomial filter with min_rho = ", min.rho))
        #         #vaf.dens = density(vaf[!low.rho.idx])
        #         #vaf.dens.max = vaf.dens$x[which.max(vaf.dens$y)]
        #         #lines(vaf.dens$x, vaf.dens$y, col="red", lwd=2)
        #         #abline(v=vaf.dens.max * c(0.5, 1), col="red", lwd=2, lty=2)
        #         abline(v=median(vaf[!low.rho.idx]) * c(0.5, 1), col="red", lty=c(1, 3), lwd=2.5)
        #         hist(vaf[low.rho.idx],
        #              breaks=seq(0, 1, 0.01), col="dodgerblue4", border="white", xlab="VAF", xpd=NA,
        #              main=paste0(prettyNum(sum(low.rho.idx), big.mark=","),
        #                          " variants failing beta-binomial filter with min_rho = ", min.rho))
        #     }
        #     invisible(dev.off())
        #     invisible(gc())
        # }
    }
    
    
    # Plot and save species spectrum
    if (species %in% names(REF.NAMES)) {
        dir.create(OUTPUT$TMP.DIR, showWarnings=F)
        name = toTitleCase(species)
        write.table(cbind(species.vars[, 1:2], ".", species.vars[, 3:4]),
                    file=paste0(OUTPUT$TMP.DIR, "/", name, "_Final.vcf"),
                    sep="\t", quote=F, col.names=F, row.names=F)
        spc = SigProfilerMatrixGeneratorR(name, REF.NAMES[species], OUTPUT$TMP.DIR,
                                          plot=T, exome=F, chrom_based=F, tsb_stat=F, seqInfo=F)
        file.copy(paste0(OUTPUT$TMP.DIR, "/output/plots/ID_83_plots_", name, ".pdf"),
                  paste0(OUTPUT$PDF.DIR.2, "/", OUTPUT$PDF.PREFIX, name, ".pdf"),
                  overwrite=T)
        file.copy(paste0(OUTPUT$TMP.DIR, "/output/ID/", name, ".ID83.all"),
                  paste0(OUTPUT$DATA.DIR, "/", OUTPUT$SPEC.PREFIX, name, ".txt"),
                  overwrite=T)
        unlink(OUTPUT$TMP.DIR, recursive=T)
    }
    
    # Delete SigProfiler folder
    unlink(OUTPUT$TMP.DIR, recursive=T)
}


# # Plot spectra of filtered variants per species and per filter
# for (i in 1:length(vars.per.species)) {
#     if (!is.null(vars.per.species[[i]])) {
#         plot_spectrum(build_catalogues(vars.per.species[[i]]),
#                       paste0(OUTPUT$PDF.DIR.2, "/", OUTPUT$PDF.PREFIX, names(vars.per.species)[i], ".pdf"))
#     }
# }
# for (i in 1:length(vars.per.filter)) {
#     if (nrow(vars.per.filter[[i]]) > 0) {
#         plot_spectrum(build_catalogues(vars.per.filter[[i]]),
#                       paste0(OUTPUT$PDF.DIR.3, "/", OUTPUT$PDF.PREFIX, names(vars.per.filter)[i], ".pdf"))
#     }
# }


# Output tables of counts per filter
colnames(counts.per.filter) =  #colnames(counts.per.filter.excl) = colnames(counts.per.rho) =
    sample.names
rownames(counts.per.filter) = c("Final_variants",
                                "Original_variants",
                                paste0("Discarded_", FILTER.NAMES))
#rownames(counts.per.rho) = paste0("Discarded_betabinom_rho=", MIN.RHO)
counts.all = counts.per.filter  #rbind(counts.per.filter, counts.per.rho)
counts.all = cbind("Set"=rownames(counts.all), counts.all)
#counts.per.filter.excl = cbind("Set"=paste0("Discarded_", FILTER.NAMES), counts.per.filter.excl)

write.table(counts.all, file=OUTPUT$COUNTS, sep="\t", row.names=F, col.names=T, quote=F)
#write.table(counts.per.filter.excl, file=OUTPUT$COUNTS.EXCL, sep="\t", row.names=F, col.names=T, quote=F)


cat("\nDone\n")
