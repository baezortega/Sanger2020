# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 10.2: FIT MUTATIONAL SIGNATURES TO OBTAIN EXPOSURES PER SPECIES, INDIVIDUAL AND SAMPLE


# Input file paths
INPUT = list(
    DATA = "../data/processed/Signatures_All.RData"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/Signatures_Definitive.RData",
    PDF.DIR = paste0("../output/", Sys.Date(), "_Signatures_Definitive"),
    PDF.EXP.SPECIES = "Exposures_Species_Definitive.pdf",
    PDF.EXP.INDIV = "Exposures_Individual_Definitive.pdf",
    PDF.EXP.SAMPLE = "Exposures_Sample_Definitive.pdf",
    PDF.REC.SPECIES = "Reconstructions_Species_Definitive.pdf",
    PDF.REC.INDIV = "Reconstructions_Individual_Definitive.pdf",
    PDF.REC.SAMPLE = "Reconstructions_Sample_Definitive.pdf",
    PDF.SIGS = "Signatures_Definitive.pdf",
    PDF.SIGS.NORM = "Signatures_Definitive_HumanNormalised.pdf",
    TXT.CAT.SAMPLE = "Catalogues_Sample.txt",
    TXT.OPP.SAMPLE = "Opportunities_Sample.txt",
    TXT.EXP.SAMPLE = "Exposures_Sample_Definitive.txt"
)


cat("Loading data and packages...\n")
suppressWarnings(library(sigfit))
data("cosmic_signatures_v3")
load(INPUT$DATA)
cat("Loaded\n")


# Create output directory
dir.create(OUTPUT$PDF.DIR, showWarnings=F)


## DEFINITIVE SET OF SIGNATURES TO FIT: species, N=3
#signatures.final = signatures.all$species[[3]]

# NEW DEFINITIVE SET OF SIGNATURES:
# To avoid mixing of signatures 1 and 5, the definitive set of signatures is
# obtained by fitting COSMIC SBS1 and extracting 2 signatures with the Fit-Ext model
cat("\n\nEXTRACTING DEFINITIVE SIGNATURES",
    "\n--------------------------------\n\n")
ITER = 15000; WARMUP = 5000; SEED = 0xC0FFEE
fitext.species.3 = fit_extract_signatures(counts.all$species,
                                          convert_signatures(cosmic_signatures_v3[1, ],
                                                             opportunities_from="human-genome"),
                                          num_extra_sigs=2, opportunities=opps.all$species,
                                          iter=ITER, warmup=WARMUP, seed=SEED)
# Retrieve and reorder signatures
sig.idx = c(1, 3, 2)
signatures.final = retrieve_pars(fitext.species.3, "signatures")
for (i in 1:length(signatures.final)) {
    signatures.final[[i]] = signatures.final[[i]][sig.idx, ]
    rownames(signatures.final[[i]]) = rownames(signatures.final[[i]])[sig.idx]
}


# Fit signatures using multinomial model with opportunities
ITER = 10000
WARMUP = 5000
CHAINS = 3
SEED = 0xC0FFEE

cat("\n\nFITTING SIGNATURES TO SPECIES CATALOGUES",
    "\n----------------------------------------\n\n")
fit.species = fit_signatures(counts.all$species, signatures.final,
                             opportunities=opps.all$species,
                             iter=ITER, warmup=WARMUP, chains=CHAINS, seed=SEED)

cat("\n\nFITTING SIGNATURES TO INDIVIDUAL CATALOGUES",
    "\n-------------------------------------------\n\n")
fit.indiv = fit_signatures(counts.all$indiv, signatures.final,
                           opportunities=opps.all$indiv,
                           iter=ITER, warmup=WARMUP, chains=CHAINS, seed=SEED)

cat("\n\nFITTING SIGNATURES TO SAMPLE CATALOGUES",
    "\n---------------------------------------\n\n")
fit.sample = fit_signatures(counts.all$sample, signatures.final,
                            opportunities=opps.all$sample,
                            iter=ITER, warmup=WARMUP, chains=CHAINS, seed=SEED)


# Retrieve exposures and reconstructions
exposures.final = list("species" = retrieve_pars(fit.species, "exposures"),
                       "indiv" = retrieve_pars(fit.indiv, "exposures"),
                       "sample" = retrieve_pars(fit.sample, "exposures"))

reconstructions.final = list("species" = retrieve_pars(fit.species, "reconstructions"),
                             "indiv" = retrieve_pars(fit.indiv, "reconstructions"),
                             "sample" = retrieve_pars(fit.sample, "reconstructions"))


# Save data
save(signatures.final, exposures.final, reconstructions.final, counts.all, opps.all,
     fit.species, fit.indiv, fit.sample, file=OUTPUT$DATA)


# Plot signatures, exposures and reconstructions
cat("\n\nPLOTTING SIGNATURES, EXPOSURES AND RECONSTRUCTIONS",
    "\n--------------------------------------------------\n\n")
rownames(counts.all$indiv) = gsub("Samples with matched normal ", "", rownames(counts.all$indiv))
plot_spectrum(signatures.final, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.SIGS))
plot_spectrum(convert_signatures(signatures.final, opportunities_to="human-genome"),
              pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.SIGS.NORM))
plot_exposures(fit.species, cex_names=1.9, margin_bottom=12,
               pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.EXP.SPECIES))
plot_exposures(counts=counts.all$indiv, exposures=exposures.final$indiv,
               signature_names=rownames(signatures.final$mean),
               cex_names=1.4, margin_bottom=19, pdf_height=12,
               pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.EXP.INDIV))
plot_exposures(fit.sample, cex_names=0.5, margin_bottom=9,
               pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.EXP.SAMPLE))
plot_reconstruction(fit.species, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.REC.SPECIES))
plot_reconstruction(fit.indiv, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.REC.INDIV))
plot_reconstruction(fit.sample, pdf_path=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$PDF.REC.SAMPLE))


# Output sample catalogues, opportunities and exposures
rownames(counts.all$sample) = sapply(strsplit(rownames(counts.all$sample), " "), `[`, 1)
rownames(exposures.final$sample$mean) = sapply(strsplit(rownames(exposures.final$sample$mean), " "),
                                               `[`, 1)
write.table(counts.all$sample, sep="\t", quote=F, col.names=T, row.names=T,
            file=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$TXT.CAT.SAMPLE))
write.table(opps.all$sample, sep="\t", quote=F, col.names=T, row.names=T,
            file=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$TXT.OPP.SAMPLE))
write.table(exposures.final$sample$mean, sep="\t", quote=F, col.names=T, row.names=T,
            file=paste0(OUTPUT$PDF.DIR, "/", OUTPUT$TXT.EXP.SAMPLE))

cat("\nDone\n")
