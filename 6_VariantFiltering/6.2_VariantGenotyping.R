# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 6.2: GENOTYPE THE VARIANTS IN EACH SAMPLE ACROSS THE SAMPLES FROM THE SAME INDIVIDUAL


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    BAM.PATHS = "../data/original/Path_SampleBAMs.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    SAMPLE.MATCH = "../data/original/CrossSpecies_SamplesPerIndividual.txt",
    VARIANTS = "${AC36}/mathijs_filters_no_shared_var_filter/${SPECIES}/${SAMPLE}/results/${SAMPLE}.txt"
)

# Output file paths
OUTPUT = list(
    LOG.DIR = "logs_6.2",
    OUT.DIR = "../data/processed/AlleleCounts_Indiv",
    CNT.PREFIX = "AlleleCounts_",
    POS.PREFIX = "VariantPos_"
)


# Memory, queue and command templates for job submission
MEM = 2000
QUEUE = "normal"
CNT.CMD = "alleleCounter -b ${BAM} -o ${OUTFILE} -m 20 -l ${VARS}"
BSUB.CMD = "bsub -G team78-grp -o ${LOG}/log.%J -q ${QUEUE} -n 1 -R \"span[hosts=1] select[mem>=${MEM}] rusage[mem=${MEM}]\" -M ${MEM} \"${CMD}\""


# Create output and log directories
dir.create(OUTPUT$LOG.DIR, showWarnings=F)
dir.create(OUTPUT$OUT.DIR, showWarnings=F)


cat("Loading data and packages...\n")
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
sample.match = read.table(INPUT$SAMPLE.MATCH, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
stopifnot(!any(duplicated(sample.match$SAMPLE_NAME)))
stopifnot(all(sample.info$SAMPLE_NAME %in% sample.match$SAMPLE_NAME))
rownames(sample.match) = sample.match$SAMPLE_NAME
cat("Loaded\n\n")


for (i in 1:nrow(sample.info)) {
    name = sample.info$SAMPLE_NAME[i]
    species = sample.info$SPECIES_NAME[i]
    ref.name = sample.info$REFERENCE_GENOME[i]
    project = sample.info$PROJECT_ID[i]
    cat("\nProcessing sample", name, "\n")
    
    # Load and output variant positions
    var.path = gsub("${AC36}", ac36.path,
                    gsub("${SAMPLE}", name,
                         gsub("${SPECIES}", species,
                              INPUT$VARIANTS, fixed=T), fixed=T), fixed=T)
    if (!file.exists(var.path)) {
        cat("WARNING: File", var.path, "not found. Skipping sample.\n")
        next
    }
    
    variants = read.table(var.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
    variants = variants[order(variants[, 1], variants[, 2]), ]
    pos.path = paste0(OUTPUT$OUT.DIR, "/", OUTPUT$POS.PREFIX, name, ".txt")
    write.table(variants[, 1:2], file=pos.path, quote=F, sep="\t", row.names=F, col.names=F)
    cat(nrow(variants), "variants read\n")
    
    
    # Run alleleCounter on each matched sample (ie. each sample from same individual)
    # (NB. Blood samples are not used as matched samples of colon samples)
    idx = match(c("NORMAL_NAME", "SPECIES_NAME"), colnames(sample.match))
    samples = as.character(sample.match[name, -idx])
    samples = samples[samples != ""]
    cat("Launching alleleCounter for samples", paste(samples, collapse=", "), "\n")
    
    for (sample.id in samples) {
        if (species %in% bam.paths$SPECIES) {
            bam.path = bam.paths$PATH[bam.paths$SPECIES == species]
        }
        else {
            bam.path = bam.paths$PATH[bam.paths$SPECIES == "default"]
        }
        bam = gsub("${SPECIES}", species,
                   gsub("${REFGENOME}", ref.name, 
                        gsub("${PROJECT}", project,
                             gsub("${SAMPLE}", sample.id,
                                  bam.path, fixed=T), fixed=T), fixed=T), fixed=T)
        if (!file.exists(bam)) {
            cat("WARNING: File", bam, "not found\n")
            next
        }
        
        cmd = gsub("${VARS}", pos.path,
                   gsub("${BAM}", bam,
                        gsub("${OUTFILE}",
                             paste0(OUTPUT$OUT.DIR, "/",
                                    OUTPUT$CNT.PREFIX, name, "_In_", sample.id, ".txt"),
                             CNT.CMD, fixed=T), fixed=T), fixed=T)
        
        system(gsub("${QUEUE}", QUEUE,
                    gsub("${MEM}", MEM,
                         gsub("${LOG}", OUTPUT$LOG.DIR,
                              gsub("${CMD}", cmd,
                                   BSUB.CMD, fixed=T), fixed=T), fixed=T), fixed=T))
    }
}

cat("\nDone\n")
