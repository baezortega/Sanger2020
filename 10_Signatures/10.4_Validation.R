# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 10.4: VALIDATE EXTRACTED SIGNATURES USING MUTATIONALPATTERNS


# Input file paths
INPUT = list(
    DATA = "../data/processed/Signatures_Definitive.RData"
)

# Output file paths
OUTPUT = list(
    PDF = paste0("../output/", Sys.Date(), "_Signatures_Validation.pdf")
)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(NMF))
suppressPackageStartupMessages(library(sigfit))
load(INPUT$DATA)
cat("Loaded\n")


# Retrieve species per individual from catalogue matrix
rownames(counts.all$indiv) = gsub("Samples with matched normal ", "", rownames(counts.all$indiv))
species.per.indiv = gsub(".*\\(", "", gsub("\\)", "", rownames(counts.all$indiv)))
species.list = unique(species.per.indiv)


# Normalise individual catalogues to human-genome opportunities
# counts.indiv.norm = convert_signatures(counts.all$indiv, opportunities_from=opps.all$indiv,
#                                        opportunities_to="human-genome")
counts.indiv.norm = convert_signatures(counts.all$indiv, opportunities_from=opps.all$indiv)
counts.indiv.norm = round(counts.indiv.norm * rowSums(counts.all$indiv))


# Run MutationalPatterns
# Extracting more than 3 signatures yields duplicated versions of signature 1 mixed
# with other components. Extracting 3 signatures yields signatures comparable to
# the sigfit ones (cosine similarities: SBS1=0.999, SBS5=0.956, SBS18=866)

sigs.nmf.3 = MutationalPatterns::extract_signatures(t(counts.indiv.norm)+1e-3, rank=3, nrun=200)
# nmf.est = nmf(t(counts.indiv.norm)+1e-3, rank=2:8, method="brunet", nrun=10, seed=123456)
# plot(nmf.est)
# sigs.nmf.4 = MutationalPatterns::extract_signatures(t(counts.indiv.norm)+1e-3, rank=4, nrun=200)
# sigs.nmf.5 = MutationalPatterns::extract_signatures(t(counts.indiv.norm)+1e-3, rank=5, nrun=200)
# sigs.nmf.6 = MutationalPatterns::extract_signatures(t(counts.indiv.norm)+1e-3, rank=6, nrun=200)
# plot_96_profile(t(convert_signatures(t(sigs.nmf.3$signatures), opportunities_to="human-genome")), cond=T)
# plot_96_profile(t(convert_signatures(t(sigs.nmf.4$signatures), opportunities_to="human-genome")), cond=T)
# plot_96_profile(t(convert_signatures(t(sigs.nmf.5$signatures), opportunities_to="human-genome")), cond=T)
# plot_96_profile(t(convert_signatures(t(sigs.nmf.6$signatures), opportunities_to="human-genome")), cond=T)


# Compare sigfit and MutationalPatterns signatures
#signatures.final = signatures.all$species[[3]]
sigs.sigfit = convert_signatures(signatures.final, opportunities_to="human-genome")
sigs.mutpat = convert_signatures(t(sigs.nmf.3$signatures), opportunities_to="human-genome")
idx = as.numeric(match_signatures(sigs.sigfit, sigs.mutpat))

cairo_pdf(OUTPUT$PDF, width=24, height=15, onefile=T)
par(mar=c(4.5, 7, 6.5, 2), oma=c(1, 0, 1, 0))
for (i in 1:nrow(sigs.mutpat)) {
    cnt = structure(round(sigs.mutpat[i, , drop=F] * 1e4),
                    dimnames=list(paste("MutationalPatterns Signature", i), NULL))
    sgs = structure(rbind(sigs.sigfit[idx[i], ], c(rep(0, 95), 1)),
                     dimnames=list(c(rownames(sigs.sigfit)[idx[i]], "Decoy"), colnames(sigs.sigfit)))
    fit = sigfit::fit_signatures(cnt, sgs)
    plot_reconstruction(fit)
}
invisible(dev.off())

cat("\nDone\n")


## SomaticSignatures
## Gives similar results to MutationalPatterns
###############################################
# suppressPackageStartupMessages(library(SomaticSignatures))
# gof_nmf = assessNumberSignatures(t(counts.indiv.norm), 2:8, nReplicates=5)
# plotNumberSignatures(gof_nmf)
# sigs.nmf = identifySignatures(t(counts.indiv.norm), 3, nmfDecomposition)
# plot_spectrum(convert_signatures(t(signatures(sigs.nmf)), opportunities_to="human-genome"),
#               "~/Desktop/sigs_SomaticSigs.pdf")
# barplot(signatures(sigs.nmf)[,4])

## HDP
## Doesn't work - obtains ≥9 signatures regardless of concentration params
###########################################################################
# suppressPackageStartupMessages(library(hdp))
# # Initialise HDP structure
# # One top grandparent DP drawing from the base distribution (ppindex 0) with its own concentration
# # parameter (cpindex 1); one parent DP per species, drawing from the grandparent distribution
# # (ppindex 1) and sharing a new concentration parameter (cpindex 2), and one DP per individual
# # catalogue, drawing from the parent distribution corresponding to the individual's species.
# hdp.mut = hdp_init(ppindex = c(0,
#                                rep(1, length(species.list)),
#                                1 + match(species.per.indiv, species.list)),  # parental node index
#                    cpindex = c(1,
#                                rep(2, length(species.list)),
#                                2 + match(species.per.indiv, species.list)),  # concentration param index
#                     hh=rep(1, 96),                              # uniform prior over the 96 mut types
#                     alphaa=rep(4, length(species.list) + 2),    # shape hyperparams for each CP
#                     alphab=rep(0.01, length(species.list) + 2))  # rate hyperparams for each CP
# # add data to leaf nodes (one per cancer sample, in row order of mut_count)
# hdp.mut = hdp_setdata(hdp.mut,
#                       dpindex = (length(species.list)+2):numdp(hdp.mut),  # indices of target nodes
#                       counts.indiv.norm)    # mutation counts, rows match up with specified dpindex
# # Run three independent posterior sampling chains
# chlist = vector("list", 3)
# for (i in 1:3) {
#   # Activate DPs, 5 initial components
#   hdp.activated = dp_activate(hdp.mut, 1:numdp(hdp.mut), initcc=5, seed=i*200)
#   chlist[[i]] = hdp_posterior(hdp.activated, burnin=5000, n=100, space=100, cpiter=3, seed=i*1e3)
# }
# hdp.multi.out = hdp_extract_components(hdp_multi_chain(chlist))
# hdp.multi.out
# par(mfrow=c(1,1)); plot_comp_size(hdp.multi.out, bty="L")
# par(mfrow=c(4,3))
# plot_comp_distn(hdp.multi.out,
#                 grouping=as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each=16)),
#                 col=c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70'),
#                 cat_names=substr(sigfit:::mut_types(), 1, 3))
# plot_dp_comp_exposure(hdp.multi.out, main_text="Exposures per individual",
#                       dpindices=(length(species.list)+2):numdp(hdp.mut),
#                       col=RColorBrewer::brewer.pal(12, "Set3"),
#                       incl_nonsig=FALSE, dpnames=rownames(counts.indiv.norm), las=2,
#                       ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
#                       leg.title = 'Signature', oma=c(10,0,0,0))

## EMu
## Doesn't work - obtains 13 signatures; if forced to extract 3, fails to deconvolute
######################################################################################
# write.table(opps.all$indiv, file="~/Desktop/emu/opps.txt", sep="\t", col.names=F, row.names=F)
# write.table(counts.all$indiv, file="~/Desktop/emu/counts.txt", sep="\t", col.names=F, row.names=F)
# # cd ~/Desktop/emu; ./EMu --mut counts.txt --opp opps.txt
# sigs.emu = read.table("~/Desktop/emu/out_3_ml_spectra.txt")
# plot_spectrum(convert_signatures(sigs.emu, opportunities_to="human-genome"), pdf_path="~/Desktop/emu/emu_sigs.pdf")
