# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 10.1: EXTRACT MUTATIONAL SIGNATURES FROM SPECIES, INDIVIDUAL AND SAMPLE CATALOGUES


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CALLABLE.PREFIX = "../data/processed/CallableGenome/CallableGenome_",
    VARS.FINAL = "${AC36}/final_variant_calls/${SPECIES}/${SAMPLE}_final_variant_calls.txt"
)

# Output file paths
OUTPUT = list(
    DATA = "../data/processed/Signatures_All.RData",
    PDF.SPECIES.DIR = paste0("../output/", Sys.Date(), "_SignatureExtraction_Species"),
    PDF.INDIV.DIR = paste0("../output/", Sys.Date(), "_SignatureExtraction_Individual"),
    PDF.SAMPLE.DIR = paste0("../output/", Sys.Date(), "_SignatureExtraction_Sample"),
    NORM.COUNTS.SPECIES = "Species_Catalogues_HumanNormalised.pdf",
    NORM.COUNTS.INDIV = "Individual_Catalogues_HumanNormalised.pdf",
    NORM.COUNTS.SAMPLE = "Sample_Catalogues_HumanNormalised.pdf",
    NORM.SIGS.SPECIES = paste0("Species_Signatures_HumanNormalised_", Sys.Date(), ".pdf"),
    NORM.SIGS.INDIV = paste0("Individual_Signatures_HumanNormalised_", Sys.Date(), ".pdf"),
    NORM.SIGS.SAMPLE = paste0("Sample_Signatures_HumanNormalised_", Sys.Date(), ".pdf"),
    GOF.SPECIES = "Species_GOF.pdf",
    GOF.INDIV = "Individual_GOF.pdf",
    GOF.SAMPLE = "Sample_GOF.pdf"
)


# Function: reverse complement (for string vectors)
rev.comp = function(nucleotide.list) {
    sapply(nucleotide.list, function(nucleotides) {
        paste(
            rev(sapply(strsplit(nucleotides, "")[[1]], function(nuc) {
                if (nuc == "A") "T"
                else if (nuc == "C") "G"
                else if (nuc == "G") "C"
                else if (nuc == "T") "A"
            })),
            collapse="")
    }, USE.NAMES=F)
}

# Function: normalise vector to sum to one
normalise = function(x) {
    x / sum(x)
}


cat("Loading data and packages...\n")
suppressWarnings(library(sigfit))
suppressPackageStartupMessages(library(Biostrings))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# # ADD HUMAN SAMPLES WITH EVIDENT COLIBACTIN EXPOSURE TO THE EXCLUSION LIST
# exclude.list = c(exclude.list,
#                  sample.info$SAMPLE_NAME[sample.info$NORMAL_NAME %in%
#                                              c("PD37449b", "PD37512a4", "PD37513a21")])


# Create output directories
dir.create(OUTPUT$PDF.SPECIES.DIR, showWarnings=F)
dir.create(OUTPUT$PDF.INDIV.DIR, showWarnings=F)
dir.create(OUTPUT$PDF.SAMPLE.DIR, showWarnings=F)


# Initialise objects for variants, counts and opportunities
species.list = unique(sample.info$SPECIES_NAME)
variants = structure(lapply(species.list, function(x) NULL), names=species.list)
counts = opportunities = structure(lapply(species.list, function(x) {
    structure(vector("list", 3), names=c("species", "indiv", "sample"))
}), names=species.list)


# Obtain variants, counts and opportunities for each species
for (species in species.list) {
    # Load reference genome
    cat("\nProcessing species:", species, "\n")
    cat("Loading reference genome\n")
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    ref.name = sample.info$REFERENCE_GENOME[species.idx][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species,
                                   gsub("${REFGENOME}", ref.name, genome.path, fixed=T), fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    # Initialise opportunities matrix
    opportunities[[species]]$sample = matrix(NA, nrow=sum(species.idx), ncol=96,
                                             dimnames=list(sample.info$SAMPLE_NAME[species.idx],
                                                           sigfit:::mut_types()))
    for (i in which(species.idx)) {
        sample.id = sample.info$SAMPLE_NAME[i]
        normal.id = sample.info$NORMAL_NAME[i]
        if (sample.id %in% exclude.list) {
            cat("WARNING: Sample", sample.id, "is in the exclusion list\n")
            next
        }
        
        # Load variants
        cat("Reading variants for sample", sample.id, "\n")
        vars = read.table(gsub("${AC36}", ac36.path,
                               gsub("${SAMPLE}", sample.id,
                                    gsub("${SPECIES}", species,
                                         INPUT$VARS.FINAL, fixed=T), fixed=T), fixed=T),
                          sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
        stopifnot(nrow(vars) > 0)
        
        # Retrieve trinucleotide contexts and build mutation table
        vars$Context = as.character(padAndClip(genome[vars$Chr],
                                               IRanges(vars$Start - 1, vars$Start + 1),
                                               Lpadding.letter=".", Rpadding.letter="."))
        stopifnot(identical(substr(vars$Context, 2, 2), vars$Ref))
        vars$Sample = sample.id
        variants[[species]] = rbind(variants[[species]],
                                    vars[, c("Sample", "Ref", "Alt", "Context")])
        
        # Load and select callable genome regions
        regions = read.table(paste0(INPUT$CALLABLE.PREFIX, sample.id, ".txt"),
                             sep="\t", header=T, as.is=T)
        genome.regions = padAndClip(genome[regions$chrom], IRanges(regions$start, regions$end),
                                    Lpadding.letter=".", Rpadding.letter=".")
        
        # Calculate trinucleotide frequencies and mutational opportunities
        freqs = colSums(trinucleotideFrequency(genome.regions))
        freqs.pyr = sapply(which(substr(names(freqs), 2, 2) %in% c("C", "T")), function(i) {
            rcomp = rev.comp(names(freqs)[i])
            freqs[i] + freqs[rcomp]
        })
        opportunities[[species]]$sample[sample.id, ] =
            normalise(freqs.pyr[substr(sigfit:::mut_types(), 1, 3)])
    }
    
    # Build mutational catalogues per sample, individual and species
    cat("Building mutational catalogues\n")
    indiv.ids = unique(sample.info$NORMAL_NAME[species.idx])
    counts[[species]]$sample = build_catalogues(variants[[species]])
    counts[[species]]$species = colSums(counts[[species]]$sample)
    counts[[species]]$indiv = t(sapply(indiv.ids, function(id) {
        sample.ids = sample.info$SAMPLE_NAME[species.idx & sample.info$NORMAL_NAME == id]
        colSums(counts[[species]]$sample[sample.ids, , drop=F])
    }))
    
    # Calculate mutational opportunities per individual and species
    opportunities[[species]]$species = normalise(colSums(opportunities[[species]]$sample))
    opportunities[[species]]$indiv = t(sapply(indiv.ids, function(id) {
        sample.ids = sample.info$SAMPLE_NAME[species.idx & sample.info$NORMAL_NAME == id]
        normalise(colSums(opportunities[[species]]$sample[sample.ids, , drop=F]))
    }))
    
    rownames(counts[[species]]$indiv) = rownames(opportunities[[species]]$indiv) =
        paste("Samples with matched normal", indiv.ids)
}


# Prepare catalogues and opportunities for extraction (per species/individual/sample)
counts.all = opps.all = list("species"=NULL, "indiv"=NULL, "sample"=NULL)
counts.all$species = t(sapply(counts, `[[`, 1))
opps.all$species = t(sapply(opportunities, `[[`, 1))
for (i in 1:length(counts)) {
    rownames(counts[[i]]$indiv) = paste0(rownames(counts[[i]]$indiv),
                                              " (", names(counts)[i], ")")
    rownames(counts[[i]]$sample) = paste0(rownames(counts[[i]]$sample),
                                          " (", names(counts)[i], ")")
    counts.all$indiv = rbind(counts.all$indiv, counts[[i]]$indiv)
    counts.all$sample = rbind(counts.all$sample, counts[[i]]$sample)
}
for (i in 1:length(opportunities)) {
    opps.all$indiv = rbind(opps.all$indiv, opportunities[[i]]$indiv)
    opps.all$sample = rbind(opps.all$sample, opportunities[[i]]$sample)
}


# Extract signatures using multinomial model with opportunities
NSIGS = 2:7
ITER = 10000
WARMUP = 3000
SEED = 0xC0FFEE

cat("\n\nEXTRACTING SIGNATURES FROM SPECIES CATALOGUES",
    "\n---------------------------------------------\n\n")
pdf(paste0(OUTPUT$PDF.SPECIES.DIR, "/", OUTPUT$GOF.SPECIES), 8, 5)
fit.species = extract_signatures(counts.all$species, nsignatures=NSIGS,
                                 opportunities=opps.all$species,
                                 iter=ITER, warmup=WARMUP, seed=SEED)
invisible(dev.off())

cat("\n\nEXTRACTING SIGNATURES FROM INDIVIDUAL CATALOGUES",
    "\n------------------------------------------------\n\n")
pdf(paste0(OUTPUT$PDF.INDIV.DIR, "/", OUTPUT$GOF.INDIV), 8, 5)
fit.indiv = extract_signatures(counts.all$indiv, nsignatures=NSIGS,
                               opportunities=opps.all$indiv,
                               iter=ITER, warmup=WARMUP, seed=SEED)
invisible(dev.off())

cat("\n\nEXTRACTING SIGNATURES FROM SAMPLE CATALOGUES",
    "\n--------------------------------------------\n\n")
pdf(paste0(OUTPUT$PDF.SAMPLE.DIR, "/", OUTPUT$GOF.SAMPLE), 8, 5)
fit.sample = extract_signatures(counts.all$sample, nsignatures=NSIGS,
                                opportunities=opps.all$sample,
                                iter=ITER, warmup=WARMUP, seed=SEED)
invisible(dev.off())


# Retrieve signatures, exposures and reconstructions
signatures.all = exposures.all = reconstructions.all =
    lapply(c("species"=NA, "indiv"=NA, "sample"=NA), function(x) {
        structure(vector(mode="list", max(NSIGS) + 1),
                  names=c(paste0("N=", 1:max(NSIGS)), "best"))
    })

for (i in NSIGS) {
    signatures.all$species[[i]] = retrieve_pars(fit.species[[i]], "signatures")
    signatures.all$indiv[[i]] = retrieve_pars(fit.indiv[[i]], "signatures")
    signatures.all$sample[[i]] = retrieve_pars(fit.sample[[i]], "signatures")
    exposures.all$species[[i]] = retrieve_pars(fit.species[[i]], "exposures")
    exposures.all$indiv[[i]] = retrieve_pars(fit.indiv[[i]], "exposures")
    exposures.all$sample[[i]] = retrieve_pars(fit.sample[[i]], "exposures")
    reconstructions.all$species[[i]] = retrieve_pars(fit.species[[i]], "reconstructions")
    reconstructions.all$indiv[[i]] = retrieve_pars(fit.indiv[[i]], "reconstructions")
    reconstructions.all$sample[[i]] = retrieve_pars(fit.sample[[i]], "reconstructions")
}


# Save data
save(signatures.all, exposures.all, reconstructions.all, counts.all, opps.all, file=OUTPUT$DATA)


# Plot signatures, exposures and reconstructions
cat("\n\nPLOTTING CATALOGUES, SIGNATURES, EXPOSURES AND RECONSTRUCTIONS",
    "\n--------------------------------------------------------------\n\n")
for (i in NSIGS) {
    plot_all(fit.species[[i]], out_path=paste0(OUTPUT$PDF.SPECIES.DIR, "/", i), prefix="Species")
    plot_all(fit.indiv[[i]], out_path=paste0(OUTPUT$PDF.INDIV.DIR, "/", i), prefix="Individual")
    plot_all(fit.sample[[i]], out_path=paste0(OUTPUT$PDF.SAMPLE.DIR, "/", i), prefix="Sample")
}

# Plot human-genome-normalised signatures and catalogues
for (i in NSIGS) {
    plot_spectrum(convert_signatures(signatures.all$species[[i]], opportunities_to="human-genome"),
                  pdf_path=paste0(OUTPUT$PDF.SPECIES.DIR, "/", i, "/", OUTPUT$NORM.SIGS.SPECIES))
    plot_spectrum(convert_signatures(signatures.all$indiv[[i]], opportunities_to="human-genome"),
                  pdf_path=paste0(OUTPUT$PDF.INDIV.DIR, "/", i, "/", OUTPUT$NORM.SIGS.INDIV))
    plot_spectrum(convert_signatures(signatures.all$sample[[i]], opportunities_to="human-genome"),
                  pdf_path=paste0(OUTPUT$PDF.SAMPLE.DIR, "/", i, "/", OUTPUT$NORM.SIGS.SAMPLE))
}
plot_spectrum(convert_signatures(counts.all$species, opportunities_to="human-genome",
                                 opportunities_from=opps.all$species),
              pdf_path=paste0(OUTPUT$PDF.SPECIES.DIR, "/", OUTPUT$NORM.COUNTS.SPECIES))
plot_spectrum(convert_signatures(counts.all$indiv, opportunities_to="human-genome",
                                 opportunities_from=opps.all$indiv),
              pdf_path=paste0(OUTPUT$PDF.INDIV.DIR, "/", OUTPUT$NORM.COUNTS.INDIV))
plot_spectrum(convert_signatures(counts.all$sample, opportunities_to="human-genome",
                                 opportunities_from=opps.all$sample),
              pdf_path=paste0(OUTPUT$PDF.SAMPLE.DIR, "/", OUTPUT$NORM.COUNTS.SAMPLE))

cat("\nDone\n")
