# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 7: PLOT HEATMAPS OF VARIANT SHARING BETWEEN SAMPLES


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    VARS.STANDARD = "${AC36}/mathijs_filters/${SPECIES}/${SAMPLE}/results/${SAMPLE}.txt",
    VARS.NO.COMP = "${AC36}/mathijs_filters_no_shared_var_filter/${SPECIES}/${SAMPLE}/results/${SAMPLE}.txt",
    VARS.FINAL = "${AC36}/final_variant_calls/${SPECIES}/${SAMPLE}_final_variant_calls.txt"
)

# Output file paths
OUTPUT = list(
    OUT.DIR = paste0("../output/", Sys.Date(), "_VariantSharing"),
    HEATMAP.PREFIX = paste0(Sys.Date(), "_VariantSharing_"),
    SPECTRA.PREFIX = paste0(Sys.Date(), "_SharedVariants_Spectra_")
)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(sigfit))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(colorRamps))
suppressPackageStartupMessages(library(Biostrings))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Create output directory
dir.create(OUTPUT$OUT.DIR, showWarnings=F)


# For each species: find variants shared between each pair of samples,
# for variants processed with/without comparison across samples
for (species in unique(sample.info$SPECIES_NAME)) {
    cat("\nProcessing species:", species, "\n")
    samples = sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species]
    normals = unique(sample.info$NORMAL_NAME[sample.info$SPECIES_NAME == species])
    
    # Load reference genome
    cat("Loading reference genome\n")
    ref.name = sample.info$REFERENCE_GENOME[sample.info$SPECIES_NAME == species][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species,
                                   gsub("${REFGENOME}", ref.name, 
                                        genome.path, fixed=T), fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    # Load variants
    vars.standard = lapply(samples, function(sample.id) {
        cat("Reading standard variants for sample", sample.id, "\n")
        var.path = gsub("${AC36}", ac36.path,
                        gsub("${SAMPLE}", sample.id,
                             gsub("${SPECIES}", species,
                                  INPUT$VARS.STANDARD, fixed=T), fixed=T), fixed=T)
        if (file.exists(var.path)) {
            vars = read.table(var.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="")
            if (nrow(vars) > 0) {
                cbind(vars, "String"=paste0(vars$Chr, ":", vars$Start, vars$Alt))
            }
            else {
                cat("WARNING: File", var.path, "is empty.\n")
                data.frame(String="", stringsAsFactors=F)
            }
        }
        else {
            cat("WARNING: File", var.path, "not found.\n")
            data.frame(String="", stringsAsFactors=F)
        }
    })
    vars.no.comp = lapply(samples, function(sample.id) {
        cat("Reading no-comparison variants for sample", sample.id, "\n")
        var.path = gsub("${AC36}", ac36.path,
                        gsub("${SAMPLE}", sample.id,
                             gsub("${SPECIES}", species,
                                  INPUT$VARS.NO.COMP, fixed=T), fixed=T), fixed=T)
        if (file.exists(var.path)) {
            vars = read.table(var.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="")
            if (nrow(vars) > 0) {
                cbind(vars, "String"=paste0(vars$Chr, ":", vars$Start, vars$Alt))
            }
            else {
                cat("WARNING: File", var.path, "is empty.\n")
                data.frame(String="", stringsAsFactors=F)
            }
        }
        else {
            cat("WARNING: File", var.path, "not found.\n")
            data.frame(String="", stringsAsFactors=F)
        }
    })
    vars.final = lapply(samples, function(sample.id) {
        cat("Reading final filtered variants for sample", sample.id, "\n")
        var.path = gsub("${AC36}", ac36.path,
                        gsub("${SAMPLE}", sample.id,
                             gsub("${SPECIES}", species,
                                  INPUT$VARS.FINAL, fixed=T), fixed=T), fixed=T)
        if (file.exists(var.path)) {
            vars = read.table(var.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="")
            if (nrow(vars) > 0) {
                cbind(vars, "String"=paste0(vars$Chr, ":", vars$Start, vars$Alt))
            }
            else {
                cat("WARNING: File", var.path, "is empty.\n")
                data.frame(String="", stringsAsFactors=F)
            }
        }
        else {
            cat("WARNING: File", var.path, "not found.\n")
            data.frame(String="", stringsAsFactors=F)
        }
    })
    names(vars.standard) = names(vars.no.comp) = names(vars.final) = samples
    
    
    # Count shared variants
    cat("Counting shared variants\n")
    shared.vars.standard = sapply(vars.standard, function(v1) {
        sapply(vars.standard, function(v2) {
            if (v1$String[1] == "" | v2$String[1] == "") {
                -1
            }
            else {
                sum(v1$String %in% v2$String)
            }
        })
    })
    shared.vars.no.comp = sapply(vars.no.comp, function(v1) {
        sapply(vars.no.comp, function(v2) {
            if (v1$String[1] == "" | v2$String[1] == "") {
                -1
            }
            else {
                sum(v1$String %in% v2$String)
            }
        })
    })
    shared.vars.final = sapply(vars.final, function(v1) {
        sapply(vars.final, function(v2) {
            if (v1$String[1] == "" | v2$String[1] == "") {
                -1
            }
            else {
                sum(v1$String %in% v2$String)
            }
        })
    })
    dimnames(shared.vars.standard) = dimnames(shared.vars.no.comp) =
        dimnames(shared.vars.final) = list(samples, samples)
    shared.vars.standard = t(shared.vars.standard[nrow(shared.vars.standard):1, ])
    shared.vars.no.comp = t(shared.vars.no.comp[nrow(shared.vars.no.comp):1, ])
    shared.vars.final = t(shared.vars.final[nrow(shared.vars.final):1, ])
    
    
    # Plot heatmaps
    cat("Plotting heatmaps\n")
    if (length(samples) > 40) {
        pdf(paste0(OUTPUT$OUT.DIR, "/", OUTPUT$HEATMAP.PREFIX, species, ".pdf"), 27, 20)
    } else {
        pdf(paste0(OUTPUT$OUT.DIR, "/", OUTPUT$HEATMAP.PREFIX, species, ".pdf"), 12, 9)
    }
    # NB. Within a noninteractive script, lattice plots are plotted using 'print'
    print(levelplot(shared.vars.final,
                    xlab="", ylab="", scales=list(x=list(rot=45, tck=0), y=list(tck=0)),
                    col.regions=colorRampPalette(c("white", "deepskyblue4")),
                    main=paste("Variant sharing in", species, "\nwith definitive variant filters"),
                    panel=function(x, y, z, ...) {  # to label values in each cell
                        panel.levelplot(x, y, z, ...)
                        panel.text(x, y, ifelse(z == -1, "-", as.character(z)), cex=0.6)
                    }))
    print(levelplot(shared.vars.standard,
                    xlab="", ylab="", scales=list(x=list(rot=45, tck=0), y=list(tck=0)),
                    col.regions=colorRampPalette(c("white", "deepskyblue4")),
                    main=paste("Variant sharing in", species, "\nwith Mathij's shared-variant filter"),
                    panel=function(x, y, z, ...) {  # to label values in each cell
                        panel.levelplot(x, y, z, ...)
                        panel.text(x, y, ifelse(z == -1, "-", as.character(z)), cex=0.6)
                    }))
    print(levelplot(shared.vars.no.comp,
                    xlab="", ylab="", scales=list(x=list(rot=45, tck=0), y=list(tck=0)),
                    col.regions=colorRampPalette(c("white", "deepskyblue4")),
                    main=paste("Variant sharing in", species, "\nwithout Mathij's shared-variant filter"),
                    panel=function(x, y, z, ...) {  # to label values in each cell
                        panel.levelplot(x, y, z, ...)
                        panel.text(x, y, ifelse(z == -1, "-", as.character(z)), cex=0.6)
                    }))
    invisible(dev.off())
    
    
    # Build spectrum of variants shared between samples in each individual
    cat("Plotting spectrum of variants shared between samples in each individual\n")
    shared.vars.indiv = NULL
    
    for (normal in normals) {
        indiv.samples = sample.info$SAMPLE_NAME[sample.info$NORMAL_NAME == normal]
        indiv.samples.idx = which(samples %in% indiv.samples)
        
        for (i in indiv.samples.idx) {
            for (j in indiv.samples.idx) {
                if (j > i) {
                    #v1 = vars.no.comp[[i]]
                    #v2 = vars.no.comp[[j]]
                    v1 = vars.final[[i]]
                    v2 = vars.final[[j]]
                    idx = (v1$String %in% v2$String) & (v1$String != "")
                    if (any(idx)) {
                        # Retrieve trinucleotide contexts
                        context = as.character(padAndClip(genome[v1$Chr[idx]],
                                                          IRanges(v1$Start[idx] - 1, v1$Start[idx] + 1),
                                                          Lpadding.letter=".", Rpadding.letter="."))
                        stopifnot(identical(as.character(substr(context, 2, 2)), v1$Ref[idx]))
                        
                        # Add shared variants to table
                        shared.vars.indiv = rbind(shared.vars.indiv, cbind(SampleID = normal,
                                                                           Ref = v1$Ref[idx],
                                                                           Alt = v1$Alt[idx],
                                                                           Context = context))
                    }
                }
            }
        }
    }
    
    # Plot spectra
    if (!is.null(shared.vars.indiv)) {
        shared.vars.indiv.spectra = build_catalogues(shared.vars.indiv)
        rownames(shared.vars.indiv.spectra) = paste("Final variants shared between samples with normal",
                                                    rownames(shared.vars.indiv.spectra))
        plot_spectrum(shared.vars.indiv.spectra,
                      pdf_path=paste0(OUTPUT$OUT.DIR, "/", OUTPUT$SPECTRA.PREFIX, species, ".pdf"))
    }
    
    invisible(gc())
}
    
cat("\nDone\n")
