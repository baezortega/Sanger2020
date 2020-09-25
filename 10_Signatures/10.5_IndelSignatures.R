# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 10.5: FIT COSMIC INDEL SIGNATURES TO OBTAIN EXPOSURES PER SAMPLE AND SPECIES


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    SIGS = "../data/original/SigProfiler_ID_signatures.csv",
    SPEC.DIR = "../data/processed/IndelSpectra",
    SPEC.PREFIX = "IndelSpectra_"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/Signatures_Indels.RData",
    PDF.DIR = paste0("../output/", Sys.Date(), "_Signatures_Indels"),
    PDF.SIGS = "Signatures_COSMIC.pdf",
    PDF.CAT.SPECIES = "Catalogues_Species.pdf",
    PDF.CAT.SAMPLE = "Catalogues_Sample.pdf",
    PDF.EXP.SPECIES = "Exposures_Species.pdf",
    PDF.EXP.SAMPLE = "Exposures_Sample.pdf",
    PDF.REC.SPECIES = "Reconstructions_Species.pdf",
    PDF.REC.SAMPLE = "Reconstructions_Sample.pdf",
    TXT.EXP.SPECIES = "Exposures_Species.txt",
    TXT.EXP.SAMPLE = "Exposures_Sample.txt"
)


cat("Loading data and packages...\n")
suppressWarnings(library(sigfit))
suppressPackageStartupMessages(library(tools))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
signatures = read.csv(INPUT$SIGS)
rownames(signatures) = signatures[, 1]
signatures = t(signatures[, -1])
cat("Loaded\n")


# Species to consider
SPECIES = c("human", "mouse", "rat", "dog")

# Colours for the indel spectrum
COLS = c(rep(c("burlywood2", "chocolate2", "darkseagreen3", "forestgreen", "peachpuff",
               "salmon2", "firebrick3", "firebrick4", "slategray2", "lightskyblue3", "steelblue3",
               "dodgerblue4"), each=6),
         "lavender", rep("thistle", 2), rep("mediumpurple", 3), rep("darkorchid4", 5))


# Create output directory
dir.create(OUTPUT$PDF.DIR, showWarnings=F)


# Load sample and species catalogues
counts.all = list(species=NULL, sample=NULL)
sample.names = NULL

for (species in SPECIES) {
    cat("\nProcessing species:", species)
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    for (sample.id in sample.info$SAMPLE_NAME[species.idx]) {
        spec = read.table(paste0(INPUT$SPEC.DIR, "/", INPUT$SPEC.PREFIX, sample.id, ".txt"),
                          sep="\t", header=T, as.is=T)
        counts.all$sample = rbind(counts.all$sample, as.numeric(spec[, 2]))
        sample.names = c(sample.names, sample.id)
    }
    spec = read.table(paste0(INPUT$SPEC.DIR, "/", INPUT$SPEC.PREFIX, toTitleCase(species), ".txt"),
                      sep="\t", header=T, as.is=T)
    counts.all$species = rbind(counts.all$species, as.numeric(spec[, 2]))
}
rownames(counts.all$species) = SPECIES
rownames(counts.all$sample) = sample.names
colnames(counts.all$species) = colnames(counts.all$sample) = colnames(signatures)


# Fit signatures using multinomial model (without opportunities)
ITER = 10000
WARMUP = 5000
CHAINS = 3
SEED = 0xC0FFEE

cat("\n\n\nFITTING SIGNATURES TO SPECIES CATALOGUES",
    "\n----------------------------------------\n\n")
fit.species = fit_signatures(counts.all$species, signatures,
                             iter=ITER, warmup=WARMUP, chains=CHAINS, seed=SEED)

cat("\n\nFITTING SIGNATURES TO SAMPLE CATALOGUES",
    "\n---------------------------------------\n\n")
fit.sample = fit_signatures(counts.all$sample, signatures,
                            iter=ITER, warmup=WARMUP, chains=CHAINS, seed=SEED)


# Retrieve exposures and reconstructions
exposures = list("species" = retrieve_pars(fit.species, "exposures"),
                 "sample" = retrieve_pars(fit.sample, "exposures"))

reconstructions = list("species" = retrieve_pars(fit.species, "reconstructions"),
                       "sample" = retrieve_pars(fit.sample, "reconstructions"))


# Save data
save(signatures, exposures, reconstructions, counts.all, fit.species, fit.sample, file=OUTPUT$DATA)


# Plot signatures, exposures and reconstructions
cat("\n\nPLOTTING SIGNATURES, EXPOSURES AND RECONSTRUCTIONS",
    "\n--------------------------------------------------\n\n")
plot_spectrum(signatures, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.SIGS), colors=COLS)
plot_spectrum(counts.all$species,
              pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.CAT.SPECIES), colors=COLS)
plot_spectrum(counts.all$sample,
              pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.CAT.SAMPLE), colors=COLS)
plot_exposures(fit.species, cex_names=1.9, margin_bottom=12,
               pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.EXP.SPECIES))
plot_exposures(fit.sample, cex_names=0.5, margin_bottom=9,
               pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.EXP.SAMPLE))
plot_reconstruction(fit.species, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.REC.SPECIES),
                    legend_pos="topright", legend_cex=1.5)
plot_reconstruction(fit.sample, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.REC.SAMPLE),
                    legend_pos="topright", legend_cex=1.5)


# Output sample catalogues, opportunities and exposures
write.table(exposures$species$mean, sep="\t", quote=F, col.names=T, row.names=T,
            file=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$TXT.EXP.SPECIES))
write.table(exposures$sample$mean, sep="\t", quote=F, col.names=T, row.names=T,
            file=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$TXT.EXP.SAMPLE))

cat("\nDone\n")
