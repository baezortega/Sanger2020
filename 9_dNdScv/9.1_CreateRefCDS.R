# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 9.1: GENERATE REFCDS DATABASES (FOR DNDSCV) FOR ANNOTATED SPECIES


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    GENOME.PATH = "../data/original/Path_RefGenomes.txt"
)

# Output file paths
OUTPUT = list(
    OUT.DIR = "../data/processed/dNdScv_RefCDS",
    TABLE.PREFIX = "BioMart_RefTable_",
    REFCDS.PREFIX = "RefCDS_"
)


cat("Loading data and packages...\n")
suppressWarnings(library(dndscv))
suppressPackageStartupMessages(library(biomaRt))
genome.path = as.character(as.matrix(read.table(INPUT$GENOME.PATH)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Create output directory
dir.create(OUTPUT$OUT.DIR, showWarnings=F)


# To list species datasets: listDatasets(useMart("ensembl", host=...))
# To list Ensembl archives: listEnsemblArchives()

# Ensembl v99 has incompatible assembly versions for
# cow (ARS-UCD1.2), horse (EquCab3.0) and human (GRCh38.p13).
# Assemblies are missing for giraffe, lion and ring-tailed lemur
# (https://www.ensembl.org/info/website/archives/assembly.html)

# To use v99: useMart("ensembl", dataset=..., host="http://jan2020.archive.ensembl.org")
# To use GRCh37 (v75): useMart("ensembl", dataset=..., host="http://grch37.ensembl.org")

ensembl.datasets = data.frame(dataset = c(cat = "fcatus_gene_ensembl",
                                          colobus = "capalliatus_gene_ensembl",
                                          cow = "btaurus_gene_ensembl",
                                          dog = "cfamiliaris_gene_ensembl",
                                          ferret = "mpfuro_gene_ensembl",
                                          horse = "ecaballus_gene_ensembl",
                                          #human = "hsapiens_gene_ensembl",
                                          mouse = "mmusculus_gene_ensembl",
                                          naked_mole_rat = "hgfemale_gene_ensembl",
                                          rabbit = "ocuniculus_gene_ensembl",
                                          rat = "rnorvegicus_gene_ensembl",
                                          tiger = "ptaltaica_gene_ensembl"),
                              host = c(cat = "http://jan2020.archive.ensembl.org",     # Felis_catus_9.0
                                       colobus = "http://jan2020.archive.ensembl.org", # Cang.pa_1.0
                                       cow = "http://grch37.ensembl.org",              # UMD3.1
                                       dog = "http://jan2020.archive.ensembl.org",     # CanFam3.1
                                       ferret = "http://jan2020.archive.ensembl.org",  # MusPutFur1.0
                                       horse = "http://grch37.ensembl.org",            # Equ Cab 2
                                       #human = "http://grch37.ensembl.org",           # GRCh37.p13 (not GRCh37d5=GRCh37.p4)
                                       mouse = "http://grch37.ensembl.org",            # GRCm38.p2 (not GRCm38)
                                       naked_mole_rat = "http://jan2020.archive.ensembl.org", # HetGla_female_1.0
                                       rabbit = "http://jan2020.archive.ensembl.org",  # OryCun2.0
                                       rat = "http://jan2020.archive.ensembl.org",     # Rnor_6.0
                                       tiger = "http://jan2020.archive.ensembl.org"),  # PanTig1.0
                              stringsAsFactors=F)


for (species in rownames(ensembl.datasets)) {
    cat("\nProcessing species:", species, "\n")
    
    # Build path to reference genome
    ref.name = gsub("merged_", "",
                    sample.info$REFERENCE_GENOME[match(species, sample.info$SPECIES_NAME)])
    ref.path = gsub("${SPECIES}", species,
                    gsub("${REFGENOME}", ref.name,
                         genome.path, fixed=T), fixed=T)
    
    # Load Ensembl dataset
    ensembl = useMart("ensembl",
                      dataset=ensembl.datasets[species, "dataset"],
                      host=ensembl.datasets[species, "host"])
    
    #View(listAttributes(ensembl))     # attributes define the values to be retrieved
    #View(istFilters(ensembl))         # filters define a restriction on the query
    #filterOptions("biotype", ensembl) # each filter can take a number of values
    
    # Retrieve data from BioMart
    attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_peptide_id",
                   "chromosome_name", "genomic_coding_start", "genomic_coding_end",
                   "cds_start", "cds_end", "cds_length", "strand")
    filters = "transcript_biotype"
    values = "protein_coding"
    col.names = c("gene.id", "gene.name", "cds.id", "chr", "chr.coding.start", "chr.coding.end",
                  "cds.start", "cds.end", "length", "strand")
        
    cds.table = structure(getBM(attributes=attributes, filters=filters, values=values, mart=ensembl),
                          names=col.names)
    
    # Write reference table to file
    cds.path = paste0(OUTPUT$OUT.DIR, "/", OUTPUT$TABLE.PREFIX, species, ".txt")
    write.table(cds.table, file=cds.path, sep="\t", quote=F, row.names=F)
    
    # Build RefCDS
    buildref(cdsfile=cds.path, genomefile=ref.path, excludechrs="MT",
             outfile=paste0(OUTPUT$OUT.DIR, "/", OUTPUT$REFCDS.PREFIX, species, ".RData"))
}

cat("\nDone\n")
