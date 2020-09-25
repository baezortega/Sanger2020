# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 6.5: GENOTYPE THE INDELS IN EACH SAMPLE ACROSS THE SAMPLES FROM THE SAME INDIVIDUAL


# Input file paths
INPUT = list(
    AC36.PATH = "../data/original/Path_ac36.txt",
    BAM.PATHS = "../data/original/Path_SampleBAMs.txt",
    GENOME.PATH = "../data/original/Path_RefGenomes.txt",
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    SAMPLE.MATCH = "../data/original/CrossSpecies_SamplesPerIndividual.txt",
    VARIANTS = "${AC36}/pindel/output/${SPECIES}/${SAMPLE}/${SAMPLE}_vs_${NORMAL}.flagged.vcf.gz"
)

# Output file paths
OUTPUT = list(
    LOG.DIR = "logs_6.5",
    DATA.DIR = "../data/processed/AlleleCounts_Indiv_Indels/data",
    OUT.DIR = "../data/processed/AlleleCounts_Indiv_Indels",
    BED.PREFIX = "Indels_PASS_"
)


# Memory, queue and command templates for job submission
MEM = 7000
QUEUE = "normal"
VAF.CMD = paste("module purge; module add singularity/3.3.0; set +e;",
                "singularity exec --cleanenv --bind /lustre:/lustre:ro --bind /nfs:/nfs:ro",
                "/software/CASM/singularity/vafcorrect/vafcorrect_5.7.0.sif",
                "cgpVaf.pl -a indel -d ${DATADIR} -o ${OUTDIR} -g ${REF} -b ${BED} -bo 1", # -e .vcf.gz",
                "-nn ${NORMAL} -tn ${SAMPLES}")
BSUB.CMD = "bsub -G team78-grp -o ${LOG}/log.%J -q ${QUEUE} -n 1 -R \"span[hosts=1] select[mem>=${MEM}] rusage[mem=${MEM}]\" -M ${MEM} \"${CMD}\""


# Create output and log directories
dir.create(OUTPUT$LOG.DIR, showWarnings=F)
dir.create(OUTPUT$DATA.DIR, showWarnings=F, recursive=T)


cat("Loading data and packages...\n")
ac36.path = as.character(as.matrix(read.table(INPUT$AC36.PATH)))
bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
sample.match = read.table(INPUT$SAMPLE.MATCH, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
stopifnot(!any(duplicated(sample.match$SAMPLE_NAME)))
stopifnot(all(sample.info$SAMPLE_NAME %in% sample.match$SAMPLE_NAME))
rownames(sample.match) = sample.match$SAMPLE_NAME
cat("Loaded\n\n")


for (i in 1:nrow(sample.info)) {
    name = sample.info$SAMPLE_NAME[i]
    normal = sample.info$NORMAL_NAME[i]
    species = sample.info$SPECIES_NAME[i]
    ref.name = sample.info$REFERENCE_GENOME[i]
    project = sample.info$PROJECT_ID[i]
    cat("\nProcessing sample", name, "\n")
    
    # Build path to reference genome
    genome = gsub("${SPECIES}", species,
                  gsub("${REFGENOME}", ref.name,
                       genome.path, fixed=T), fixed=T)
    
    # Load indels and output 'PASS' indels to BED file
    vars = gsub("${AC36}", ac36.path,
                gsub("${SAMPLE}", name,
                     gsub("${NORMAL}", normal,
                          gsub("${SPECIES}", species,
                               INPUT$VARIANTS, fixed=T), fixed=T), fixed=T), fixed=T)
    if (!file.exists(vars)) {
        cat("WARNING: File", vars, "not found. Skipping sample.\n")
        next
    }
    variants = read.table(gzfile(vars), sep="\t", header=F, as.is=T, col.names=1:11)
    if (nrow(variants) == 0) {
        cat("WARNING: No variants found in input file. Skipping sample.\n")
        next
    } 
    variants = variants[variants[, 7] == "PASS", c(1, 2, 4, 5)]
    variants = variants[order(variants[, 1], variants[, 2]), ]
    bed.path = paste0(OUTPUT$DATA.DIR, "/", OUTPUT$BED.PREFIX, name, ".txt")
    write.table(variants, file=bed.path, quote=F, sep="\t", row.names=F, col.names=F)
    cat(nrow(variants), "'PASS' indels read\n")
    
    # Create symbolic links to BAMs for matched samples (ie. samples from same individual)
    # (NB. Blood samples are not used as matched samples of colon samples)
    if (species %in% bam.paths$SPECIES) {
        bam.path = bam.paths$PATH[bam.paths$SPECIES == species]
    }
    else {
        bam.path = bam.paths$PATH[bam.paths$SPECIES == "default"]
    }
    idx = match(c("NORMAL_NAME", "SPECIES_NAME"), colnames(sample.match))
    samples = as.character(sample.match[name, -idx])
    samples = samples[samples != ""]
    
    for (sample.id in c(samples, normal)) {
        bam = gsub("${SPECIES}", species,
                   gsub("${REFGENOME}", ref.name, 
                        gsub("${PROJECT}", project,
                             gsub("${SAMPLE}", sample.id,
                                  bam.path, fixed=T), fixed=T), fixed=T), fixed=T)
        if (!file.exists(bam)) {
            cat("WARNING: File", bam, "not found\n")
            next
        }
        for (suffix in c(" ", ".bai ", ".bas ")) {
            system(paste0("ln -sf ", bam, suffix, OUTPUT$DATA.DIR, "/", sample.id, ".bam", suffix))
        }
    }
    
    # # Create symbolic links to VCFs for matched samples
    # for (sample.id in samples) {
    #    vars = gsub("${AC36}", ac36.path,
    #                gsub("${SAMPLE}", sample.id,
    #                     gsub("${NORMAL}", normal,
    #                          gsub("${SPECIES}", species,
    #                               INPUT$VARIANTS, fixed=T), fixed=T), fixed=T), fixed=T)
    #     if (!file.exists(vars)) {
    #         cat("WARNING: File", vars, "not found\n"); next
    #     }
    #     system(paste0("ln -sf ", vars, " ", OUTPUT$DATA.DIR, "/", sample.id, ".vcf.gz"))
    # }
    
    # Launch cgpVaf
    cat("Launching cgpVaf for samples", paste(samples, collapse=", "), "\n")
    dir.create(paste0(OUTPUT$OUT.DIR, "/", name), showWarnings=F)
    
    cmd = gsub("${DATADIR}", OUTPUT$DATA.DIR,
               gsub("${OUTDIR}", paste0(OUTPUT$OUT.DIR, "/", name),
                    gsub("${REF}", genome,
                         gsub("${BED}", bed.path,
                              gsub("${NORMAL}", normal,
                                   gsub("${SAMPLES}", paste(samples, collapse=" "), VAF.CMD,
                                        fixed=T), fixed=T), fixed=T), fixed=T), fixed=T), fixed=T)
    system(gsub("${QUEUE}", QUEUE,
                gsub("${MEM}", MEM,
                     gsub("${LOG}", OUTPUT$LOG.DIR,
                          gsub("${CMD}", cmd,
                               BSUB.CMD, fixed=T), fixed=T), fixed=T), fixed=T))
}

cat("\nDone\n")
