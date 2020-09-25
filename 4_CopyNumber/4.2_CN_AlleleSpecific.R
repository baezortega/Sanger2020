# CROSS-SPECIES ANALYSIS PIPELINE
# Adrian Baez-Ortega, 2020

# STEP 4.2: INFER ALLELE-SPECIFIC COPY NUMBER USING COVERAGE AND BAF


# Read sample index from command line
SAMPLE.IDX = as.integer(commandArgs(TRUE)[1])
cat(SAMPLE.IDX, "\n\n")


# Input file paths
INPUT = list(
    SAMPLE.INFO = "../data/original/CrossSpecies_ProjectInfo.txt",
    EXCLUDE.LIST = "../data/processed/SamplesToExclude.txt",
    COV.PREFIX = "../data/processed/BinCoverage/BinCov_",
    SNPS.PATH = "../data/original/Path_SampleSNPs.txt"
    #BAM.PATHS = "../data/original/Path_SampleBAMs.txt",
    #SAMPLE.MATCH = "../data/original/CrossSpecies_SamplesPerIndividual.txt",
)

# Output file paths
OUTPUT = list(
    DATA.DIR = "../data/processed/CN_AlleleSpecific",
    DATA.PREFIX = "CN_AlleleSpecific_",
    CN.DIR = paste0("../output/", Sys.Date(), "_CN_AlleleSpecific"),
    CN.A.PREFIX = paste0(Sys.Date(), "_CNSegments_Allele_"),
    CN.T.PREFIX = paste0(Sys.Date(), "_CNSegments_Total_"),
    PLOT.PREFIX = paste0(Sys.Date(), "_CNPlots_")
    #HIST.DIR = paste0("../output/", Sys.Date(), "_HetSNPs_BAF_Hist"),
    #HIST.PREFIX = paste0(Sys.Date(), "_HetSNPs_BAF_Hist_"),
)


# Species to consider
SPECIES = c("cat", "cow", "dog", "horse", "human", "mouse", "rabbit", "rat")

# Read length (bp; for coverage calculation)
READ.LEN = 150

# Minimum CN segment length (bins)
MIN.SEGMENT = 5

# Minimum bin coverage
MIN.COV = 10

# Coverage and VAF thresholds for het SNPs
HET.NR = 15
HET.VAF = c(0.4, 0.6)  #c(0.35, 0.65)

# Minimum VAF difference for SNP phasing within a bin
MIN.DIFF = 0.15

# Penalty parameter for CN changes
PENALTY = 0.3

# Disable scientific notation
options(scipen=999)


# Create output directories
#dir.create(OUTPUT$HIST.DIR, showWarnings=F)
dir.create(OUTPUT$DATA.DIR, showWarnings=F)
dir.create(OUTPUT$CN.DIR, showWarnings=F)


cat("Loading data and packages...\n")
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(bbmle))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(deepSNV))
suppressPackageStartupMessages(library(emdbook))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicRanges))
snps.path = as.character(as.matrix(read.table(INPUT$SNPS.PATH)))
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE.LIST)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
#bam.paths = read.table(INPUT$BAM.PATHS, sep="\t", header=T, as.is=T)
#sample.match = read.table(INPUT$SAMPLE.MATCH, sep="\t", header=T, as.is=T)
#stopifnot(!any(duplicated(sample.match$SAMPLE_NAME)))
#stopifnot(all(sample.info$SAMPLE_NAME %in% sample.match$SAMPLE_NAME))
#rownames(sample.match) = sample.match$SAMPLE_NAME
cat("Loaded\n")

#sample.id="MD6260ab_lo0004"; normal.id="MD6260z"; species="mouse"
sample.id = sample.info$SAMPLE_NAME[SAMPLE.IDX]
normal.id = sample.info$NORMAL_NAME[SAMPLE.IDX]
species = sample.info$SPECIES_NAME[SAMPLE.IDX]

#idx = match(c("SAMPLE_NAME", "SPECIES_NAME"), colnames(sample.match))
#matched.samples = as.character(sample.match[sample.id, -idx])
#matched.samples = matched.samples[matched.samples != "" & !(matched.samples %in% exclude.list)]
#spcs = ifelse(species %in% bam.paths$SPECIES, species, "default")
#matched.samples.bam.paths = sapply(matched.samples, function(id) {
#    gsub("${SPECIES}", species,
#         gsub("${REFGENOME}", ref.name, 
#              gsub("${PROJECT}", project,
#                   gsub("${SAMPLE}", sample.id,
#                        bam.paths$PATH[bam.paths$SPECIES == spcs],
#                        fixed=T), fixed=T), fixed=T), fixed=T)
#})

#sample.snps.path="../data/tmp/MD6260ab_lo0004_vs_MD6260z.snps.ids.vcf.gz"
sample.snps.path = gsub("${SPECIES}", species,
                        gsub("${SAMPLE}", sample.id,
                             gsub("${NORMAL}", normal.id,
                                  snps.path, fixed=T), fixed=T), fixed=T)


# Check that species and sample are valid
if (!(species %in% SPECIES)) {
    cat("WARNING:", species, "not included in the species list. Skipping sample.\n")
} else if (sample.id %in% exclude.list) {
    cat("WARNING: Sample", sample.id, "is in the exclusion list. Skipping sample.\n")
} else if (!file.exists(sample.snps.path)) {
    cat("WARNING: File", sample.snps.path, "not found. Skipping sample.\n")
} else {
    
    cat("\nProcessing sample:", sample.id, "\n")
    
    # Load bin coverage for sample and normal
    cat("Loading bin coverage data\n")
    bin.cov = lapply(c(sample.id, normal.id), function(id) {
        cov = read.table(paste0(INPUT$COV.PREFIX, id, ".txt"), sep="\t",
                         col.names=c("Chrom", "Start", "End", "Cov", "", "Length", ""), as.is=T)
        cov$Start = cov$Start + 1
        cov$AdjCov = cov$Cov / cov$Length * READ.LEN
        if (species == "cat") {
            cov.order = order(cov$Chrom, cov$Start)
        } else {
            cov.order = order(as.numeric(gsub("X", 100, cov$Chrom)), cov$Start)
        }
        cov[cov.order, c("Chrom", "Start", "End", "Length", "Cov", "AdjCov")]
    })
    bin.table = bin.cov[[1]][, c("Chrom", "Start", "End")]
    
    
    # Load germline SNPs
    cat("Loading SNPs VCF and retrieving heterozygous SNPs\n")
    snps = read.table(gzfile(sample.snps.path), sep="\t", as.is=T,
                      col.names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                  "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOUR"),
                      colClasses=c("CHROM"="character"))
    
    # Retrieve VAFs and isolate heterozygous germline SNPs
    snps$NR_N = apply(str_split_fixed(snps$NORMAL, ":", 10)[, 2:9], 1,
                      function(x) sum(as.numeric(x)))
    snps$NR_T = apply(str_split_fixed(snps$TUMOUR, ":", 10)[, 2:9], 1,
                      function(x) sum(as.numeric(x)))
    snps$VAF_N = as.numeric(str_split_fixed(snps$NORMAL, ":", 10)[, 10])
    snps$VAF_T = as.numeric(str_split_fixed(snps$TUMOUR, ":", 10)[, 10])
    snps.het = snps[snps$NR_N >= HET.NR &
                    snps$NR_T >= HET.NR &
                    snps$VAF_N >= HET.VAF[1] &
                    snps$VAF_N <= HET.VAF[2], ]
    rm(snps); invisible(gc())
    
    ## Genotype het SNPs across related crypts and filter based on global VAF
    #snps.het.gt = sapply(matched.samples, function(id) {
    #    bam2R(bam.path, het.snps$CHROM[i], het.snps$POS[i], het.snps$POS[i], q=30, mask=3844, mq=10)
    #})
    ## [Sum reads across samples and calculate global BAF]
    ## UNFINISHED...
    
    
    # Assign het SNPs to genome coverage bins
    bin.table.gr = with(bin.table, GRanges(Chrom, IRanges(Start, End)))
    snps.het.gr = with(snps.het, GRanges(CHROM, IRanges(POS, POS)))
    snps.het$BIN = factor(findOverlaps(snps.het.gr, bin.table.gr, select="first"),
                          levels=as.character(1:length(bin.table.gr)))
    snps.het = snps.het[!is.na(snps.het$BIN), ]
    cat("Assigned", nrow(snps.het), "heterozygous SNPs to genome bins\n")
    
    
    ## Plot histograms of SNPs per bin and tumour VAF
    #cairo_pdf(paste0(OUTPUT$HIST.DIR, "/", OUTPUT$HIST.PREFIX, sample.id, ".pdf"), 15, 5.57)
    #par(mfrow=c(1, 2), mar=c(5, 4.75, 6, 0.5))
    #hist(table(snps.het$BIN), breaks=seq(-0.5, max(table(snps.het$BIN))+0.5, 1),
    #     xlim=c(0, 100), col="dodgerblue4", border="white", freq=T,
    #     xlab="SNPs per bin", ylab="Bins", main="\nHeterozygous SNPs per genome bin")
    #mtext(paste0("(Coverage ≥", HET.NR, " and VAF in [", HET.VAF[1], ", ",
    #             HET.VAF[2], "] in matched normal)"), side=3, line=0)
    #mtext(paste("Total:", prettyNum(nrow(snps.het), big.mark=","), "heterozygous SNPs"),
    #      side=3, line=-1.75)
    #mtext(paste(prettyNum(sum(table(snps.het$BIN) == 0)), "/", prettyNum(nrow(bin.table)),
    #            "bins without SNPs"), side=3, line=-3.5)
    #hist(snps.het$VAF_T, breaks=seq(0, 1, 0.01), xlim=c(0, 1), freq=T,
    #     col="dodgerblue4", border="white", xlab="Tumour BAF", ylab="Heterozygous SNPs",
    #     main="\nTumour BAF of heterozygous SNPs")
    #title(paste0(sample.id, " (", species, ")"), outer=T, line=-2.25, cex.main=1.3)
    #invisible(dev.off())
    
    
    # Infer allele-specific CN for bins with variants, total CN for bins without variants
    # (code adapted from Inigo Martincorena, Science 2015 [call_targeted_copy_number.R])
    
    # Define two sets of bins:
    #  (i) bins with coverage ≥15 in sample and matched normal, and ≥1 variant (for allele CN)
    # (ii) bins with cov ≥15 in sample and normal, but with no variants (for total CN only)
    bin.idx = 1:nrow(bin.table) %in% snps.het$BIN
    set1.idx = bin.cov[[1]]$AdjCov >= MIN.COV & bin.cov[[2]]$AdjCov >= MIN.COV & bin.idx
    set2.idx = bin.cov[[1]]$AdjCov >= MIN.COV & bin.cov[[2]]$AdjCov >= MIN.COV & !bin.idx
    cat("Set 1: ", sum(set1.idx), " bins with coverage ≥", MIN.COV, " and ≥1 SNP\n",
        "Set 2: ", sum(set2.idx), " bins with coverage ≥", MIN.COV, " and 0 SNPs\n", sep="")
    
    
    ## a. Beta-binomial model for BAF (assuming no bias in the mean het BAF, ie p=0.5)
    # Calculate variant allele counts from sample BAFs (=VAFs)
    cov.obs = bin.cov[[1]]$AdjCov[set1.idx]            # observed bin coverage (Set 1)
    cov.obs.s2 = round(bin.cov[[1]]$AdjCov[set2.idx])  # observed bin coverage (Set 2)
    vaf.obs = sapply(which(set1.idx), function(bin) {
        idx = snps.het$BIN == bin
        if (sum(idx) > 1) {
            # If the bin has multiple SNPs with very distant VAFs, ensure that
            # all the VAFs are on the same side of 0.5 (to improve CN1/LOH detection)
            #vaf.1st = snps.het$VAF_T[idx][1]
            #for (i in which(idx)[-1]) {
            #    if ((vaf.1st < 0.5 & snps.het$VAF_T[i] > 0.5 & snps.het$VAF_T[i] - vaf.1st > 0.1) |
            #        (vaf.1st > 0.5 & snps.het$VAF_T[i] < 0.5 & vaf.1st - snps.het$VAF_T[i] > 0.1)) {
            #        snps.het$VAF_T[i] = 1 - snps.het$VAF_T[i]
            #    }
            #}
            vaf.1st = snps.het$VAF_T[idx][1]
            chg.idx = rep(FALSE, sum(idx))
            if (vaf.1st < 0.5) {
                chg.idx = snps.het$VAF_T[idx] > 0.5 & snps.het$VAF_T[idx] - vaf.1st > MIN.DIFF
            }
            else if (vaf.1st > 0.5) {
                chg.idx = snps.het$VAF_T[idx] < 0.5 & vaf.1st - snps.het$VAF_T[idx] > MIN.DIFF
            }
            snps.het$VAF_T[idx][chg.idx] = 1 - snps.het$VAF_T[idx][chg.idx]
        }
        median(snps.het$VAF_T[idx])    # median observed bin BAF (Set 1)
    })
    nv.obs = round(vaf.obs * cov.obs)  # median number of supporting reads (Set 1)
    cov.obs = round(cov.obs)
    
    # Fit beta-binomial to supporting and total reads to estimate BAF dispersion
    bb.theta = tryCatch({
        fitdistr(nv.obs, emdbook::dbetabinom, prob=0.5,
                 size=cov.obs, start=list(theta=10), method="L-BFGS-B")$estimate
    }, error = function(e) {
        fitdistr(nv.obs, emdbook::dbetabinom, prob=0.5,
                 size=cov.obs, start=list(theta=10), method="BFGS")$estimate
    })
    cat("bb.theta =", bb.theta, "\n")
    
    
    ## b. Negative-binomial model for the relative coverage
    # Instead of read counts per gene, we use total coverage per 100-kb bin;
    # the expected coverage of a bin is calculated using the median sample coverage
    # and the median ratio between sample and matched normal coverage
    cov.exp = bin.cov[[2]]$AdjCov * median(bin.cov[[1]]$AdjCov / bin.cov[[2]]$AdjCov, na.rm=T)
    cov.exp.s2 = cov.exp[set2.idx]  # Expected coverage (Set 2)
    cov.exp = cov.exp[set1.idx]     # Expected coverage (Set 1)
    
    # Fit negative binomial to total coverage to estimate coverage dispersion
    # ('size' = dispersion parameter = shape par of the gamma mixing distribution)
    nb.size = tryCatch({
        fitdistr(c(cov.obs, cov.obs.s2), dnbinom,
                 mu=c(cov.exp, cov.exp.s2), start=list(size=1), method="L-BFGS-B")$estimate
    }, error = function(e) {
        fitdistr(c(cov.obs, cov.obs.s2), dnbinom,
                 mu=c(cov.exp, cov.exp.s2), start=list(size=1), method="BFGS")$estimate
    })
    cat("nb.size =", nb.size, "\n")
    
    
    ## c. Subfunction: Log-likelihood (joint BAF and coverage)
    loglik.copynumber = function(a1, a2, exp.cov, nA, nB, rho, nb.size, bb.theta, ref.BAF) {
        # Define the beta-binomial implementation to be used
        dbetabinom = emdbook::dbetabinom
        
        # Expected BAF and coverage given the copy number (rho, nA and nB)
        new.baf = (1 - rho + (rho * nA)) / (2 * (1 - rho) + rho * (nA + nB))
        new.baf = new.baf * ref.BAF / (new.baf * ref.BAF + (1 - new.baf) * (1 - ref.BAF))
        new.mu = exp.cov + exp.cov * (nA + nB - 2) * rho / 2
        
        # Likelihoods
        # If a1=0 (no BAF info), return only neg binom likelihood (for total CN);
        # otherwise combine neg binom and beta-binom likelihoods (for allele CN)
        if (a1 == 0) { 
            LL = sum(dnbinom(x=a1+a2, mu=new.mu, size=nb.size, log=T))
        } else {
            ll.baf = sum(dbetabinom(x=a1, size=a1+a2, prob=new.baf, theta=bb.theta, log=T))
            ll.cov = sum(dnbinom(x=a1+a2, mu=new.mu, size=nb.size, log=T))
            LL = ll.baf + ll.cov
        }
        return(-LL)
    }
    
    
    ## d. Subfunction: Optimising nA, nB and rho by exhaustive search of nA and nB
    ##    and numerical optimisation of rho
    optimise.copynumber = function(n.values, a1, a2, exp.cov,
                                   nb.size, bb.theta, penalty.matrix, ref.BAF) {
        LL.mat = array(NA, dim=c(length(n.values), length(n.values)))
        rho.mat = array(NA, dim=c(length(n.values), length(n.values)))
        rho.CIlow.mat = array(NA, dim=c(length(n.values), length(n.values)))
        rho.CIhigh.mat = array(NA, dim=c(length(n.values), length(n.values)))
        for (j in 1:length(n.values)) {
            for (h in 1:length(n.values)) {
                nA = n.values[j]
                nB = n.values[h]
                # Original pars: start=list(rho=0.5), lower=c(rho=0.0001)
                m = mle2(loglik.copynumber, start=list(rho=0.925),
                         data=list(a1=a1, a2=a2, exp.cov=exp.cov, nA=nA, nB=nB,
                                   nb.size=nb.size, bb.theta=bb.theta, ref.BAF=ref.BAF),
                         method="Brent", lower=c(rho=0.85), upper=c(rho=0.9999))
                rho.mat[j, h] = summary(m)@coef[1]
                LL.mat[j, h] = logLik(m)[1]
            }
        }
        
        # Penalise complexity with a prior (obtain posterior probabilities)
        PP.mat = LL.mat + penalty.matrix
        likely.cn = which(exp(PP.mat - max(PP.mat)) >= 0.10, arr.ind=T)
        
        # Report all the likely solutions in order of their posterior prob
        likely.cndf = data.frame(nA=n.values[likely.cn[, 1]], nB=n.values[likely.cn[, 2]],
                                 rhoMLE=rho.mat[likely.cn], LL=LL.mat[likely.cn],
                                 PP=PP.mat[likely.cn], rho2.5=NA, rho97.5=NA)
        ## Confidence interval of rho (by likelihood profile)
        #for (j in 1:nrow(likely.cndf)) {
        #    # Original pars: start=list(rho=0.5), lower=c(rho=0.0001)
        #    m = mle2(loglik.copynumber, start=list(rho=0.925),
        #             data=list(a1=a1, a2=a2, exp.cov=exp.cov, nb.size=nb.size, bb.theta=bb.theta,
        #                       ref.BAF=ref.BAF, nA=likely.cndf[j, 1], nB=likely.cndf[j, 2]),
        #             method="Brent", lower=c(rho=0.85), upper=c(rho=0.9999))
        #    likely.cndf[j, c("rho2.5","rho97.5")] = tryCatch(confint(m), error=function(e) {c(NA,NA)})
        #}
        
        likely.cndf = likely.cndf[order(likely.cndf$PP, decreasing=T), ]
        likely.cndf$rel.prob = exp(likely.cndf$PP - max(likely.cndf$PP))
        return(likely.cndf)
    }
    
    
    ## e. Subfunction: Generating penalty matrix against complex solutions, which
    ##    applies a fixed penalty ('param') for every copy number step away from 1-1;
    ##    values closer to 1/0 imply weaker/stronger penalisation
    get.penalty.matrix = function(n.values, param) {
        n = length(n.values)
        penalty.matrix = array(NA, dim=c(n, n))
        for (j in 1:n) {
            for (h in 1:n) {
                penalty.matrix[j,h] = log(param ^ sum(abs(c(n.values[j], n.values[h]) - 1)))
            }
        }
        ## Alternative: A matrix of the frequency of copy number states from whole-genome TCGA data
        #cntcga = read.table("/nfs/team78pc10/im3/Reference_data/TCGA_CNAs_freqs_diploidSamples.txt")
        #colnames(cntcga) = c("chr","start","end","nMajor","nMinor")    
        #cnafreqs = array(0,dim=c(length(n.values),length(n.values)))
        #for (j in 1:(length(n.values)-1)) {
        #    for (h in j:length(n.values)) {
        #        pos = (cntcga[,4]==n.values[j] & cntcga[,5]==n.values[h]) | (cntcga[,5]==n.values[j] & cntcga[,4]==n.values[h])
        #        cnafreqs[j,h] = sum(cntcga[pos,3]-cntcga[pos,2])
        #    }
        #}
        #cnafreqs = cnafreqs/sum(cnafreqs)
        #for (j in 1:(length(n.values)-1)) {
        #    for (h in j:length(n.values)) {
        #        cnafreqs[h,j] = cnafreqs[j,h]
        #    }
        #}
        #penalty.matrix = log(cnafreqs)
        return(penalty.matrix)
    }
    
    
    ## f. Search for the best copy number for every bin in the current sample
    n.values = 0:4
    p.vec = rep(NA, sum(set1.idx))
    penalty.matrix = get.penalty.matrix(n.values, PENALTY)
    
    ## ORIGINAL VERSION: FIND SIGNIFICANT BINS BASED ON BAF ONLY
    #for (i in 1:sum(set1.idx)) {
    #    # Bin-specific negative binomial (more conservative, downweights effect of coverage)
    #    #nb.size = fitdistr(cov.obs[i], "negative binomial", mu=cov.exp[i], start=list(size=1), method="BFGS")$estimate
    #    # Code to adjust for reference bias in BAF (omitted):
    #    #prob0 = 0.5
    #    #prob1 = vaf.obs[i]
    #    #prob0 = prob0*ref.BAF / (prob0*ref.BAF + (1-prob0)*(1-ref.BAF))
    #    #prob1 = prob1*ref.BAF / (prob1*ref.BAF + (1-prob1)*(1-ref.BAF))
    #    # Obtain p-value using BAF only (LRT)
    #    LL0 = sum(dbetabinom(x=nv.obs[i], size=cov.obs[i], prob=0.5, theta=bb.theta, log=T))
    #    LL1 = sum(dbetabinom(x=nv.obs[i], size=cov.obs[i], prob=vaf.obs[i], theta=bb.theta, log=T))
    #    p.vec[i] = 1 - pchisq(2 * (LL1-LL0), 1)
    #}
    ## Multiple-testing correction
    #qvals = p.adjust(p.vec, method="BH")
    #signif = which(qvals < 0.05)
    #cat(length(signif), "significant CNV candidates found in this sample\n")
    
    ## ORIGINAL VERSION: CALCULATE COPY NUMBER FOR BINS WITH SIGNIFICANT BAF
    ## If there are any significant hits, plot BAF and relative coverage
    #if (length(signif) > 0) {
    #    cn.events = NULL
    #    for (i in signif) {
    #        a1 = nv.obs[i]
    #        a2 = cov.obs[i] - a1
    #        cn = optimise.copynumber(n.values, a1, a2, cov.exp[i], nb.size, bb.theta,
    #                                 penalty.matrix, ref.BAF=0.5)
    #        j = which(set1.idx)[i]
    #        cn = cbind(data.frame(Sample=rep(sample.id, nrow(cn)),
    #                              Bin=rep(paste0(bin.table$Chrom[j], ":",
    #                                             bin.table$Start[j], "-",
    #                                             bin.table$End[j]), nrow(cn)),
    #                              stringsAsFactors=F),
    #                   cn)
    #        cn.events = rbind(cn.events, cn)
    #    }
    #    # Output candidate events
    #    write.table(cn.events, sep="\t", quote=F, row.names=F, col.names=T,
    #                file=paste0(OUTPUT$CN.DIR, "/", OUTPUT$CN.PREFIX, sample.id, ".txt"))
    #} else {
    #    cat("No significant CNV events found in this sample\n",
    #        file=paste0(OUTPUT$CN.DIR, "/", OUTPUT$CN.PREFIX, sample.id, ".txt"))
    #}
    
    ## MODIFIED VERSION: CALCULATE MOST LIKELY COPY NUMBER FOR EACH BIN
    cat("Estimating CN for each bin\n")
    cn.table = data.frame(CN_total = rep(NA, nrow(bin.table)),
                          CN_allele = rep(NA, nrow(bin.table)),
                          nMajor = rep(NA, nrow(bin.table)),
                          nMinor = rep(NA, nrow(bin.table)),
                          nA = rep(NA, nrow(bin.table)),
                          nB = rep(NA, nrow(bin.table)),
                          rho_MLE = rep(NA, nrow(bin.table)),
                          LL = rep(NA, nrow(bin.table)),
                          PP = rep(NA, nrow(bin.table)),
                          rho_2.5 = rep(NA, nrow(bin.table)),
                          rho_97.5 = rep(NA, nrow(bin.table)),
                          rel_prob = rep(NA, nrow(bin.table)))
    
    # Obtain most likely CN state for bins in Sets 1 and 2
    cat("\nSet 1\n")
    for (i in 1:sum(set1.idx)) {
        if (i %% 100 == 0) cat(i, "/", sum(set1.idx), "\n")
        cn = optimise.copynumber(n.values,
                                 a1=nv.obs[i], a2=cov.obs[i]-nv.obs[i], exp.cov=cov.exp[i],
                                 nb.size, bb.theta, penalty.matrix, ref.BAF=0.5)
        n.major = max(cn[1, 1:2])
        n.minor = min(cn[1, 1:2])
        cn.table[which(set1.idx)[i], ] = c(n.major + n.minor,
                                           paste0(n.major, "-", n.minor),
                                           n.major, n.minor, cn[1, ])
    }
    cat("\nSet 2\n")
    for (i in 1:sum(set2.idx)) {
        if (i %% 100 == 0) cat(i, "/", sum(set2.idx), "\n")
        # For a1=0, only total CN is estimated
        cn = optimise.copynumber(n.values,
                                 a1=0, a2=cov.obs.s2[i], exp.cov=cov.exp.s2[i],
                                 nb.size, bb.theta, penalty.matrix, ref.BAF=0.5)
        cn.table[which(set2.idx)[i], ] = c(sum(cn[1, 1:2]), NA, NA, NA, cn[1, ])
    }
    
    
    ## g. Define CN segments above minimum length
    ##    (segmentation is done separately for total CN and allele-specific CN)
    cat("\nDefining CN segments\n")
    cn.segs.prelim = list("total"=NULL, "allele"=NULL)
    for (name in c("total", "allele")) {
        start.idx = 1
        while (start.idx < nrow(cn.table)) {
            chrom = bin.table$Chrom[start.idx]
            chrom.end = rev(which(bin.table$Chrom == chrom))[1]
            colname = paste0("CN_", name)
            cn = cn.table[start.idx, colname]
            if (is.na(cn)) {
                start.idx = start.idx + 1
            }
            else {
                ## Bins with CN=NA are not counted as ends for segments with CN=2
                #if (cn == 2) {
                #    end.idx = start.idx + which(cn.table[-(1:start.idx), colname] != cn) - 1
                #} else {
                end.idx = start.idx - 1 + 
                    which((cn.table[-(1:start.idx), colname] != cn) |
                              is.na(cn.table[-(1:start.idx), colname]))
                #}
                if (length(end.idx) == 0 | end.idx[1] > chrom.end) {
                    end.idx = chrom.end
                }
                else {
                    end.idx = end.idx[1]
                }
                if (end.idx - start.idx + 1 >= MIN.SEGMENT) {
                    if (name == "total") {
                        n.major = n.minor = NA
                    }
                    else {
                        n.major = cn.table$nMajor[start.idx]
                        n.minor = cn.table$nMinor[start.idx]
                    }
                    cn.segs.prelim[[name]] = rbind(cn.segs.prelim[[name]],
                                                   data.frame(StartBin=start.idx, EndBin=end.idx,
                                                              Chrom=chrom,
                                                              Start=bin.table$Start[start.idx],
                                                              End=bin.table$End[end.idx],
                                                              Length=bin.table$End[end.idx] -
                                                                  bin.table$Start[start.idx] + 1,
                                                              CN=cn, nMajor=n.major, nMinor=n.minor,
                                                              stringsAsFactors=F))
                }
                start.idx = end.idx + 1
            }
        }
    }
    
    # Merge adjacent segments with same CN and separated by a gap below the min segment length
    cn.segments = list("total"=NULL, "allele"=NULL)
    for (name in c("total", "allele")) {
        i = 1
        while (i <= nrow(cn.segs.prelim[[name]])) {
            segment = cn.segs.prelim[[name]][i, ]
            j = i + 1
            while (j <= nrow(cn.segs.prelim[[name]]) &
                   cn.segs.prelim[[name]]$Chrom[j] == segment$Chrom &
                   cn.segs.prelim[[name]]$CN[j] == segment$CN &
                   cn.segs.prelim[[name]]$StartBin[j] - segment$EndBin <= MIN.SEGMENT) {
                segment$EndBin = cn.segs.prelim[[name]]$EndBin[j]
                segment$End = cn.segs.prelim[[name]]$End[j]
                segment$Length = segment$End - segment$Start + 1
                j = j + 1
            }
            cn.segments[[name]] = rbind(cn.segments[[name]], segment)
            i = j
        }
        rownames(cn.segments[[name]]) = NULL
    }
    
    
    ## h. Plot coverage, BAF and CN
    cat("Plotting and saving CN segments\n")
    png(paste0(OUTPUT$CN.DIR, "/", OUTPUT$PLOT.PREFIX, sample.id, ".png"), 18, 10, "in", res=120)
    par(mfrow=c(4, 1), mar=c(0, 5, 3, 1.25), oma=c(2.5, 0, 3, 0))
    
    chroms = unique(bin.table$Chrom)
    cutoffs = c(match(chroms, bin.table$Chrom), nrow(bin.table))
    cols = rep("black", length(chroms))
    cols[seq_along(cols) %% 2 == 0] = "darkgrey"
    cols.s2 = rep(cols, table(bin.table$Chrom)[chroms])[set2.idx]
    cols = rep(cols, table(bin.table$Chrom)[chroms])[set1.idx]
    xlm = c(1, nrow(bin.table))
    offset = 0.06; cxlab = 1.5; cxmain = 1.7

    # Plot obs/exp coverage ratio
    plot(which(set1.idx), cov.obs / cov.exp,
         col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="Obs/Exp coverage ratio", 
         xpd=NA, xlim=xlm, xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
         main="Observed/Expected bin coverage ratio", cex.lab=cxlab, cex.main=cxmain)
    points(which(set2.idx), cov.obs.s2 / cov.exp.s2, col=cols.s2, pch=16, cex=0.6)
    title(paste0(sample.id, " (normal ", normal.id, "; ",
                 prettyNum(sum(set1.idx), big.mark=","), " bins with SNPs)"),
          cex.main=2, line=3.5, xpd=NA)
    
    # Plot het BAF
    plot(which(set1.idx), vaf.obs,
         col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="BAF",
         xlim=xlm, ylim=c(0, 1), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
         main="Heterozygous SNP BAF", cex.lab=cxlab, cex.main=cxmain)
    
    # Plot CN per bin
    plot(cn.table$CN_total,
         col="forestgreen", pch=15, cex=0.7, las=2, xlab="", ylab="Copy number",
         xlim=xlm, ylim=c(0, 4), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
         main="Copy number per bin", cex.lab=cxlab, cex.main=cxmain)
    points(cn.table$nMajor + offset, col="red", pch=15, cex=0.7)
    points(cn.table$nMinor - offset, col="blue", pch=15, cex=0.7)
    
    # Plot CN segments
    plot(cn.table$CN_total, type="n", las=2, xlab="", ylab="Copy number",
         xlim=xlm, ylim=c(0, 4), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
         main="Copy number segments (≥5 bins)", cex.lab=cxlab, cex.main=cxmain)
    segments(x0=cn.segments$total$StartBin, x1=cn.segments$total$EndBin,
             y0=cn.segments$total$CN, col="forestgreen", lwd=4, lend=2)
    segments(x0=cn.segments$allele$StartBin, x1=cn.segments$allele$EndBin,
             y0=cn.segments$allele$nMajor + offset, col="red", lwd=4, lend=2)
    segments(x0=cn.segments$allele$StartBin, x1=cn.segments$allele$EndBin,
             y0=cn.segments$allele$nMinor - offset, col="blue", lwd=4, lend=2)
    
    mtext(chroms, side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.7, font=2)
    invisible(dev.off())
    
    
    ## i. Output and save CN segments
    write.table(cn.segments$total, sep="\t", quote=F, row.names=F, col.names=T,
                file=paste0(OUTPUT$CN.DIR, "/", OUTPUT$CN.T.PREFIX, sample.id, ".txt"))
    write.table(cn.segments$allele, sep="\t", quote=F, row.names=F, col.names=T,
                file=paste0(OUTPUT$CN.DIR, "/", OUTPUT$CN.A.PREFIX, sample.id, ".txt"))
    save(bin.table, cn.table, cn.segments,
         set1.idx, set2.idx, cov.obs, cov.exp, vaf.obs, cov.obs.s2, cov.exp.s2,
         file=paste0(OUTPUT$DATA.DIR, "/", OUTPUT$DATA.PREFIX, sample.id, ".RData"))
}

cat("\nDone\n")
