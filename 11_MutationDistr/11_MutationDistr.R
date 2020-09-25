# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 11: PLOT DISTRIBUTIONS OF MUTATIONS ALONG THE GENOME


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    CHROM.LISTS = "../data/original/CrossSpecies_ChromLists.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    VARS.FINAL = "${AC36}/final_variant_calls/${SPECIES}/${SAMPLE}_final_variant_calls.txt",
    INDELS.FINAL = "${AC36}/final_indel_calls/${SPECIES}/${SAMPLE}_final_indel_calls.txt"
)

# Output file paths
OUTPUT = list(
    OUT.DIR = paste0("../output/", Sys.Date(), "_MutationDistr"),
    PDF.PREFIX = paste0(Sys.Date(), "_MutationDistr_")
)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(GenomicRanges))
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
chrom.lists = read.table(INPUT$CHROM.LISTS, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Species to consider
SPECIES = c("cat", "cow", "dog", "horse", "human", "mouse", "rabbit", "rat")

# Mutation types
MUT.TYPES = c("G>T"="C>A", "G>C"="C>G", "G>A"="C>T", "A>T"="T>A",
              "A>G"="T>C", "A>C"="T>G", "Ins"="Ins", "Del"="Del")

# Mutation colours
MUT.COLS = c("deepskyblue", "black", "firebrick2", "gray76",
             "darkolivegreen3", "rosybrown2", "darkorange1", "steelblue4")

# Genome bin length (bp)
BIN.LEN = 1e6

# Function: log10 axis
log10.axis = function(side, at, ...) {
    at.minor = log10(outer(1:9, 10^(min(at):(max(at)-1))))[-1, ]
    lab = sapply(at, function(i) as.expression(bquote(10^ .(i))))
    axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
    axis(side=side, at=at, labels=lab, ...)
}

# Disable scientific notation
options(scipen=999)


# Create output directories
dir.create(OUTPUT$OUT.DIR, showWarnings=F)
#dir.create(OUTPUT$DAT.DIR, showWarnings=F)


for (species in SPECIES) {
    cat("\nProcessing species:", species)
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    ref.name = sample.info$REFERENCE_GENOME[species.idx][1]
    project = sample.info$PROJECT_ID[species.idx][1]
    
    # Define non-overlapping bins along chromosomal contigs
    cat("\nDefining genome bins\n")
    chr.lengths = chrom.lists[chrom.lists$Species == species, c("Chromosome", "End")]
    bin.table = NULL
    for (i in 1:nrow(chr.lengths)) {
        bin.table = rbind(bin.table,
                          cbind(chr.lengths$Chromosome[i],
                                seq(0, chr.lengths$End[i], by=BIN.LEN),
                                c(seq(BIN.LEN, chr.lengths$End[i], by=BIN.LEN), chr.lengths$End[i])))
    }
    bin.table = data.frame("Chr"=bin.table[, 1],
                           "Start"=as.integer(bin.table[, 2]) + 1,
                           "End"=as.integer(bin.table[, 3]),
                           stringsAsFactors=F)
    #if (species != "cat") {
    #    idx = order(as.numeric(gsub("X", 100, gsub("Y", 200, bin.table$Chr))), bin.table$Start)
    #} else {
    #    idx = order(bin.table$Chr, bin.table$Start)
    #}
    #bin.table = bin.table[idx, ]
    
    
    # Process variants for each sample
    vars.species = NULL
    for (sample.id in sample.info$SAMPLE_NAME[species.idx]) {
        cat("Processing sample:", sample.id, "\n")
        normal.id = sample.info$NORMAL_NAME[match(sample.id, sample.info$SAMPLE_NAME)]
        invisible(gc())
        
        # Load SNVs and indels
        paths = c(gsub("${AC36}", ac36.path,
                       gsub("${SAMPLE}", sample.id,
                            gsub("${SPECIES}", species,
                                 INPUT$VARS.FINAL, fixed=T), fixed=T), fixed=T),
                  gsub("${AC36}", ac36.path,
                       gsub("${SAMPLE}", sample.id,
                            gsub("${SPECIES}", species,
                                 INPUT$INDELS.FINAL, fixed=T), fixed=T), fixed=T))
        if (!all(file.exists(paths))) {
            cat("WARNING: File(s)", paste(paths[!file.exists(paths)], collapse=", "),
                "not found. Skipping sample.\n")
            next
        }
                  
        vars = read.table(paths[1], sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
        inds = read.table(paths[2], sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chrom"="character", "Ref"="character", "Alt"="character"))
        if (nrow(vars) == 0) {
            cat("WARNING: No variants found in input file. Skipping sample.\n")
            next
        }
        
        # Obtain variant distances and types for SNVs and indels
        dist.snvs = c(diff(vars$Start), NA)
        dist.snvs[vars$Chr != c(vars$Chr[-1], "")] = NA
        type.snvs = paste0(vars$Ref, ">", vars$Alt)
        for (i in 1:length(MUT.TYPES)) {
            type.snvs[type.snvs == names(MUT.TYPES)[i]] = MUT.TYPES[i]
        }
        if (nrow(inds) > 0) {
            dist.inds = c(diff(inds$Pos), NA)
            dist.inds[inds$Chrom != c(inds$Chrom[-1], "")] = NA
            type.inds = ifelse(nchar(inds$Alt) > nchar(inds$Ref), "Ins", "Del")
        }
        
        # Add variants to table
        cat("Adding", nrow(vars), "SNVs to table\n")
        vars.species = rbind(vars.species,
                             data.frame(SampleID = sample.id,
                                        NormalID = normal.id,
                                        Chr = vars$Chr,
                                        Pos = vars$Start,
                                        Type = factor(type.snvs, levels=MUT.TYPES),
                                        Dist = dist.snvs,
                                        stringsAsFactors=F))
        if (nrow(inds) == 0) {
            cat("No indels found\n")
        } else {
            cat("Adding", nrow(inds), "indels to table\n")
            vars.species = rbind(vars.species,
                                 data.frame(SampleID = sample.id,
                                            NormalID = normal.id,
                                            Chr = inds$Chrom,
                                            Pos = inds$Pos,
                                            Type = factor(type.inds, levels=MUT.TYPES),
                                            Dist = dist.inds,
                                            stringsAsFactors=F))
        }
    }
    
    # Assign variants to genome bins
    cat("Assigning variants to genome bins\n")
    bin.table.gr = with(bin.table, GRanges(Chr, IRanges(Start, End)))
    vars.species.gr = with(vars.species, GRanges(Chr, IRanges(Pos, Pos)))
    vars.species$Bin = findOverlaps(vars.species.gr, bin.table.gr, select="first")
    vars.species$Bin_fct = factor(as.character(vars.species$Bin),
                                  levels=as.character(1:nrow(bin.table)))
    
    # Remove variants shared between samples from the same donor
    dup.idx = duplicated(vars.species[, 2:5])
    vars.species.dedup = vars.species[!dup.idx, ]
    cat("Discarding", sum(dup.idx), "duplicated variants from species table\n")
    cat(nrow(vars.species.dedup), "variants after deduplication\n")
    
    
    # Calculate variant density per bin for each sample and mutation type
    mut.density = structure(vector(mode="list", length=length(unique(vars.species$SampleID)) + 1),
                            names=c(species, unique(vars.species$SampleID)))
    mut.density[[1]] = table(vars.species.dedup[, c("Type", "Bin_fct")])
    for (id in unique(vars.species$SampleID)) {
        mut.density[[id]] = table(vars.species[vars.species$SampleID == id, c("Type", "Bin_fct")])
    }
    
    
    # Plot variant density and distance along the genome
    cairo_pdf(paste0(OUTPUT$OUT.DIR, "/", OUTPUT$PDF.PREFIX, species, ".pdf"), 26, 9, onefile=T)
    par(mfrow=c(2, 1), mar=c(1, 5, 0, 2), oma=c(2, 0, 4, 0))
    cutoffs = c(match(unique(bin.table$Chr), bin.table$Chr), nrow(bin.table))
    
    for (i in 1:length(mut.density)) {
        #barplot(mut.density[[i]], col=MUT.COLS, ...)
        b = barplot(colSums(mut.density[[i]]), col="grey50", border=NA, space=0, cex.lab=1.6,
                    cex.axis=1.3, las=2, xaxs="i", ylim=c(0, max(colSums(mut.density[[i]]))*1.05), 
                    xaxt="n", xlab="", ylab="Mutations per Mb")
        abline(v=b[cutoffs])
        legend("topleft", legend=MUT.TYPES, col=MUT.COLS,
               horiz=T, pch=16, bty="n", pt.cex=1.3, inset=c(0, -0.1), xpd=NA)
        box()
        if (i == 1) {
            idx = rep(T, nrow(vars.species))
        } else {
            idx = vars.species$SampleID == names(mut.density)[i]
        }
        plot(vars.species$Bin[idx], log10(vars.species$Dist[idx]), pch=16, cex=0.5, cex.lab=1.6,
             xaxt="n", yaxt="n", xlab="", xaxs="i", ylab="Inter-mutation distance (bp)",
             col=MUT.COLS[match(vars.species$Type[idx], MUT.TYPES)])
        abline(v=cutoffs)
        mtext(unique(bin.table$Chr), side=1, line=0.75, cex=1.3,
              at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2)
        log10.axis(side=2, at=0:round(max(log10(vars.species$Dist[idx]), na.rm=T)),
                   cex.axis=1.3, las=2)
        title(paste0(ifelse(i == 1, "All mutations in ", "Mutations in "), names(mut.density)[i],
                     " (n = ", prettyNum(sum(idx), big.mark=","), ")"), outer=T, cex.main=1.7)
    }
    invisible(dev.off())
    
}

cat("\nDone\n")
