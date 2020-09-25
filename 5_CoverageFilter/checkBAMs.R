INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    BAM.PATHS = "../data/original/Path_SampleBAMs.txt"
)

sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))

for (i in 1:nrow(sample.info)) {
    sample.id = sample.info$SAMPLE_NAME[i]
    species = sample.info$SPECIES_NAME[i]
    ref.name = sample.info$REFERENCE_GENOME[i]
    project = sample.info$PROJECT_ID[i]
    #cat(sample.id, "\n", sep="")
    
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
        cat(bam, "not found\n")
    }
}

cat("Done\n")

