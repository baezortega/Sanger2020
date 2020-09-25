# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 9.2: RUN DNDSCV FOR ANNOTATED SPECIES


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    CONTIG.COORDS = "../data/original/Path_ContigCoords.txt",
    REFCDS.PREFIX = "../data/processed/dNdScv_RefCDS/RefCDS_",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    CALLABLE.PREFIX = "../data/processed/CallableGenome/CallableGenome_",
    GNOMAD.LOF = "../data/original/GnomAD.v2.1.1.lof_metrics.by_gene.txt",
    VARS.FINAL = "${AC36}/final_variant_calls/${SPECIES}/${SAMPLE}_final_variant_calls.txt"
)

# Output file paths
OUTPUT = list(
    OUT.DIR = paste0("../output/", Sys.Date(), "_dNdScv"),
    DAT.DIR = "../data/processed/dNdScv_Results",
    OUT.PREFIX = paste0(Sys.Date(), "_dNdScv_"),
    DAT.PREFIX = "dNdScv_Results_",
    COUNTS = paste0(Sys.Date(), "_CodingCounts.txt"),
    PDF = paste0(Sys.Date(), "_dNdS_PerSpecies.pdf")
)


cat("Loading data and packages...\n")
suppressWarnings(library(dndscv))
suppressWarnings(library(RColorBrewer))
suppressPackageStartupMessages(library(GenomicRanges))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
gnomad.lof = read.table(INPUT$GNOMAD.LOF, sep="\t", header=T, as.is=T)
coords.path = as.character(as.matrix(read.table(INPUT$CONTIG.COORDS)))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Create output directories
dir.create(OUTPUT$OUT.DIR, showWarnings=F)
dir.create(OUTPUT$DAT.DIR, showWarnings=F)


# Identify haploinsufficient genes
# (LOEUF values in lower 10% in gnomAD; contained in column 'oe_lof_upper')
loeuf.10.genes = gnomad.lof$gene[(gnomad.lof$oe_lof_upper <=
                                      quantile(gnomad.lof$oe_lof_upper, probs=0.1, na.rm=T)) %in% T]


# Run dNdScv for each species
species.dnds.all = species.dnds.hpli = annot.table = coding.bp = coding.call.bp = NULL
for (species in unique(sample.info$SPECIES_NAME)) {
    cat("\nProcessing species:", species, "\n")
    if (species == "human") {
        refcds = "hg19"
    }
    else {
        refcds = paste0(INPUT$REFCDS.PREFIX, species, ".RData")
        if (!file.exists(refcds)) {
            cat("WARNING: RefCDS not found. Skipping species.\n")
            next
        }
    }
    
    # Build mutation table
    vars.species = NULL
    for (i in which(sample.info$SPECIES_NAME == species)) {
        sample.id = sample.info$SAMPLE_NAME[i]
        if (sample.id %in% exclude.list) {
            cat("INFO: Sample", sample.id, "is in the exclusion list\n")
            next
        }
        
        cat("Reading variants for sample", sample.id, "\n")
        vars.path = gsub("${AC36}", ac36.path,
                         gsub("${SAMPLE}", sample.id,
                              gsub("${SPECIES}", species,
                                   INPUT$VARS.FINAL, fixed=T), fixed=T), fixed=T)
        if (!file.exists(vars.path)) {
            cat("WARNING: file", vars.path, "not found. Skipping sample.\n")
            next
        }
        
        vars = read.table(vars.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
        if (nrow(vars) == 0) {
            cat("WARNING: No variants found in input file. Skipping sample.\n")
            next
        }
        
        cat("Adding", nrow(vars), "variants to table\n")
        vars.species = rbind(vars.species,
                             data.frame(sampleID = sample.id,
                                        normalID = sample.info$NORMAL_NAME[i],
                                        chr = vars$Chr,
                                        pos = vars$Start,
                                        ref = vars$Ref,
                                        alt = vars$Alt,
                                        stringsAsFactors=F))
    }
    
    # Remove variants shared between samples from the same individual
    dup.idx = duplicated(vars.species[, -1])
    vars.species = vars.species[!dup.idx, -2]
    cat("Discarding", sum(dup.idx), "duplicated variants\n")
    cat(nrow(vars.species), "variants after deduplication\n")
    
    # For species with merged assemblies, back-translate variant coordinates
    ref.name = sample.info$REFERENCE_GENOME[match(species, sample.info$SPECIES_NAME)]
    if (grepl("merged", ref.name)) {
        cat("Remapping variant coordinates for genome", ref.name, "\n")
        contig.coords = read.table(gsub("${SPECIES}", species,
                                        gsub("${REFGENOME}", ref.name,
                                             coords.path, fixed=T), fixed=T),
                                   sep="\t", header=F, as.is=T,
                                   col.names=c("contig", "chrom", "start", "end"))
        contig.names = sapply(strsplit(contig.coords$contig, " "), `[`, 1)
        contig.gr = makeGRangesFromDataFrame(contig.coords)
        vars.gr = makeGRangesFromDataFrame(vars.species, start.field="pos", end.field="pos")
        vars.contig.idx = findOverlaps(vars.gr, contig.gr, select="first")
        stopifnot(!any(is.na(vars.contig.idx)))
        vars.species$chr = contig.names[vars.contig.idx]
        vars.species$pos = vars.species$pos - contig.coords$start[vars.contig.idx] + 1
    }
    
    # Identify genes within the callable genome of every sample
    if (species == "human") {
        data("refcds_hg19", package="dndscv")
    } else {
        load(refcds)
    }
    gene.names = gene.chrs = gene.starts = gene.ends = rep(NA, length(RefCDS))
    for (i in 1:length(RefCDS)) {
        gene = RefCDS[[i]]
        gene.names[i] = gene$gene_name
        gene.chrs[i] = gene$chr
        gene.starts[i] = gene$intervals_cds[1, 1]
        gene.ends[i] = gene$intervals_cds[nrow(gene$intervals_cds), 2]
    }
    stopifnot(all(gene.starts < gene.ends))
    gene.gr = makeGRangesFromDataFrame(data.frame(chr = gene.chrs,
                                                  start = gene.starts, end = gene.ends))
    
    callable.idx = rep(TRUE, length(gene.gr))
    for (sample.id in sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species]) {
        if (!(sample.id %in% exclude.list)) {
            if (!file.exists(paste0(INPUT$CALLABLE.PREFIX, sample.id, ".txt"))) {
                cat("WARNING: callable genome file not found for sample", sample.id, "\n")
                next
            }
            call.gen = read.table(paste0(INPUT$CALLABLE.PREFIX, sample.id, ".txt"),
                                  sep="\t", header=T, as.is=T)
            # For  merged assemblies, back-translate variant coordinates
            if (grepl("merged", ref.name)) {
                regions.contig.idx = findOverlaps(makeGRangesFromDataFrame(call.gen),
                                                  contig.gr, select="first")
                stopifnot(!any(is.na(regions.contig.idx)))
                call.gen$chrom = contig.names[regions.contig.idx]
                call.gen$start = call.gen$start - contig.coords$start[regions.contig.idx] + 1
                call.gen$end = call.gen$end - contig.coords$start[regions.contig.idx] + 1
            }
            callable.idx = callable.idx & overlapsAny(gene.gr, makeGRangesFromDataFrame(call.gen))
        }
    }
    hpli.callable.idx = callable.idx & (toupper(gene.names) %in% loeuf.10.genes)
    cat(sum(callable.idx), "total genes in callable genome,", sum(!callable.idx), "genes excluded\n")
    cat(sum(hpli.callable.idx), "haploinsufficient genes in callable genome\n\n")
    
    
    # Run dNdScv on all genes in callable genome
    dnds.out.all = dndscv(vars.species, refdb=refcds, gene_list=gene.names[callable.idx],
                             max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf)
    
    # Run dNdScv on haploinsufficient genes in callable genome
    dnds.out.hpli = dndscv(vars.species, refdb=refcds, gene_list=gene.names[hpli.callable.idx],
                           max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf, outp=1)

    # Output global and genewise dN/dS results
    signif.idx = (dnds.out.all$sel_cv$qallsubs_cv < 0.05 |
        dnds.out.all$sel_cv$qmis_cv < 0.05 | dnds.out.all$sel_cv$qtrunc_cv < 0.05) %in% TRUE
    cat(sum(signif.idx), "genes with significant dN/dS\n")
    out.path = paste0(OUTPUT$OUT.DIR, "/", OUTPUT$OUT.PREFIX, species, ".txt")
    cat("Global dN/dS for ", species,
        "\n\nAll genes in callable genome (n = ", sum(callable.idx), "):\n\n", sep="", file=out.path)
    write.table(dnds.out.all$globaldnds, file=out.path, append=T, sep="\t", quote=F, row.names=F)
    cat("\n\nHaploinsufficient genes in callable genome (n = ", sum(hpli.callable.idx), "):\n\n",
        sep="", file=out.path, append=T)
    write.table(dnds.out.hpli$globaldnds, file=out.path, append=T, sep="\t", quote=F, row.names=F)
    
    cat("\n\n", sum(signif.idx), " genes with significant dN/dS\n\n", sep="", file=out.path, append=T)
    if (any(signif.idx)) {
        write.table(dnds.out.all$sel_cv[signif.idx, ],
                    file=out.path, append=T, sep="\t", quote=F, row.names=F)
    }

    # Collect coding genome length, callable coding genome length,
    # annotated mutation counts, and global dN/dS, and save output
    coding.bp = c(coding.bp, structure(sum(width(gr_genes)), names=species))
    coding.call.bp = c(coding.call.bp, structure(sum(width(gr_genes)[gr_genes@elementMetadata$names %in%
                                                                         gene.names[callable.idx]]),
                                                 names=species))
    annot.table = rbind(annot.table,
                        table(dnds.out.all$annotmuts[, c("sampleID", "impact")])[, c("Essential_Splice", "Missense", "Nonsense", "Synonymous"), drop=F])
    species.dnds.all = c(species.dnds.all, structure(list(dnds.out.all$globaldnds), names=species))
    species.dnds.hpli = c(species.dnds.hpli, structure(list(dnds.out.hpli$globaldnds), names=species))
    save(dnds.out.all, dnds.out.hpli,
         file=paste0(OUTPUT$DAT.DIR, "/", OUTPUT$DAT.PREFIX, species, ".RData"))
}
print(warnings())


# Plot global dN/dS per species
pdf(paste0(OUTPUT$OUT.DIR, "/", OUTPUT$PDF), 15, 6)
par(mar=c(3, 5, 3, 1.75))
cols = c(mis="#377EB8", non="#E41A1C", spl="#4DAF4A", tru="#FF7F00", all="#984EA3", NA, NA)
ylim = c(-4, 4)

# dN/dS in all genes
dnds.mle = sapply(species.dnds.all, function(dnds) structure(dnds$mle, names=rownames(dnds)))
dnds.cilow = sapply(species.dnds.all, function(dnds) structure(dnds$cilow, names=rownames(dnds)))
dnds.cihigh = sapply(species.dnds.all, function(dnds) structure(dnds$cihigh, names=rownames(dnds)))
plot(seq(1, length(species.dnds.all) * (nrow(dnds.mle)+2)),
     log2(as.numeric(rbind(dnds.mle, NA, NA))),
     main="Global dN/dS per species (all genes)",
     col=cols, pch=16, xaxt="n", yaxt="n", ylab="dN/dS", xlab="", cex.lab=1.2, cex=1.1,
     ylim=ylim, xlim=c(2.75, length(species.dnds.all)*(nrow(dnds.mle)+2) - 3.75),
     panel.first=abline(h=0, col="grey"))
segments(x0=seq(1, length(species.dnds.all) * (nrow(dnds.mle)+2)),
         y0=log2(as.numeric(rbind(dnds.cilow, NA, NA))),
         y1=log2(as.numeric(rbind(dnds.cihigh, NA, NA))),
         lwd=2.1, col=cols)
axis(side=2, at=ylim[1]:ylim[2], labels=2^(ylim[1]:ylim[2]), las=1, cex.axis=0.9)
axis(side=1, labels=names(species.dnds.all), tick=F, cex.axis=1.12,
     at=seq((nrow(dnds.mle)+1) / 2, by=nrow(dnds.mle)+2, length=length(species.dnds.all)))
abline(v=seq(nrow(dnds.mle)+1.5, by=nrow(dnds.mle)+2, length=length(species.dnds.all)-1), col="grey")
legend("topleft", legend=rownames(dnds.mle),
       inset=c(0, -0.08), cex=1.1, pch=16, col=cols, horiz=T, bty="n", xpd=T)

# dN/dS in haploinsufficient genes
dnds.mle = sapply(species.dnds.hpli, function(dnds) structure(dnds$mle, names=rownames(dnds)))
dnds.cilow = sapply(species.dnds.hpli, function(dnds) structure(dnds$cilow, names=rownames(dnds)))
dnds.cihigh = sapply(species.dnds.hpli, function(dnds) structure(dnds$cihigh, names=rownames(dnds)))
plot(seq(1, length(species.dnds.hpli) * (nrow(dnds.mle)+2)),
     log2(as.numeric(rbind(dnds.mle, NA, NA))),
     main="Global dN/dS per species (haploinsufficient genes)",
     col=cols, pch=16, xaxt="n", yaxt="n", ylab="dN/dS", xlab="", cex.lab=1.2, cex=1.1,
     ylim=ylim, xlim=c(2.75, length(species.dnds.hpli)*(nrow(dnds.mle)+2) - 3.75),
     panel.first=abline(h=0, col="grey"))
segments(x0=seq(1, length(species.dnds.hpli) * (nrow(dnds.mle)+2)),
         y0=log2(as.numeric(rbind(dnds.cilow, NA, NA))),
         y1=log2(as.numeric(rbind(dnds.cihigh, NA, NA))),
         lwd=2.1, col=cols)
axis(side=2, at=ylim[1]:ylim[2], labels=2^(ylim[1]:ylim[2]), las=1, cex.axis=0.9)
axis(side=1, labels=names(species.dnds.hpli), tick=F, cex.axis=1.12,
     at=seq((nrow(dnds.mle)+1) / 2, by=nrow(dnds.mle)+2, length=length(species.dnds.hpli)))
abline(v=seq(nrow(dnds.mle)+1.5, by=nrow(dnds.mle)+2, length=length(species.dnds.hpli)-1), col="grey")
legend("topleft", legend=rownames(dnds.mle),
       inset=c(0, -0.08), cex=1.1, pch=16, col=cols, horiz=T, bty="n", xpd=T)

invisible(dev.off())


# Output table of annotated mutation counts
species.idx = sample.info$SPECIES_NAME[match(rownames(annot.table), sample.info$SAMPLE_NAME)]
annot.table = cbind("Sample"=rownames(annot.table), as.data.frame(annot.table),
                    "Coding_Bp"=coding.bp[species.idx],
                    "Callable_Coding_Bp"=coding.call.bp[species.idx])
write.table(annot.table, file=paste0(OUTPUT$OUT.DIR, "/", OUTPUT$COUNTS),
            sep="\t", quote=F, row.names=F, col.names=T)

cat("\nDone\n")
