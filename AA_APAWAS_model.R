# module load PLINK/1.9b_5.2 GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0

library(glmnet)
set.seed(1024)

### config
args <- commandArgs(trailingOnly = TRUE)

targetGene <- args[1] ### test: NM_000699.4|AMY2A|chr1|+ / NR_146324.2|MIB2|chr1|+ / NR_146144.1|RNA45SN2|chr21|+
geneid <- paste0(strsplit(targetGene, "\\|")[[1]][1], "_", strsplit(targetGene, "\\|")[[1]][2])

plink_path <- 'plink'

SNPs_in_sumstat <- "/gpfs52/nobackup/sbcs/pingj2/Komen/sumstat/AFR.snps"

workDir <- "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/models/APAWAS/"
setwd(workDir)

tmp_folder <- paste0(workDir, "/tmp/")
exp_folder <- paste0(workDir, "/exp/")
res_folder <- paste0(workDir, "/res/")

covs <- read.table("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/APA_WAS/AA_APAWAS.covs.txt", sep = "\t", header = T, row.names = 1, as.is = T)

genotype_path <- "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/geno/Komen/Komen"

### Some functions ###
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load APA annotation file
cat(' INFO loading APA position annotation ...\n')
apa_info <- read.table("/gpfs52/nobackup/sbcs/pingj2/db/hg38_refseq_extracted_3UTR.bed", header = F, sep = "\t", as.is = T)

# mkdir tmp folder. will be cleaned
cat(' INFO mkdir tmp folders ...\n')
options(warn = -1)
dir.create(tmp_folder)
dir.create(paste0(tmp_folder, '/', geneid))
tmp_folder <- paste0(tmp_folder,'/', geneid)
options(warn = 0)

# Get chr pos for the gene
chr <- as.numeric(sub('^...', '', apa_info$V1[which(apa_info$V4 == targetGene)]))
pos_from <- apa_info$V2[which(apa_info$V4 == targetGene)]
pos_to <- apa_info$V3[which(apa_info$V4 == targetGene)]
pos_from_500K <- max(1, pos_from - 500000)
pos_to_500K <- pos_to + 500000

# Load expression
cat(' INFO loading expression data ...\n')
exptab <- loadRData("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/APA_WAS/Komen.AFR.PDUI_QN.Rdata")
exp.all <- as.data.frame(t(t(exptab[, targetGene])))
rownames(exp.all) <- rownames(exptab)
colnames(exp.all) <- "PDUI"

exp <- t(t(exp.all[which(!is.na(exp.all$PDUI)), ]))
rownames(exp) <- rownames(exp.all)[which(!is.na(exp.all$PDUI))]
colnames(exp) <- "PDUI"

if (shapiro.test(exp)$p < 0.05) {
  exp <- apply(exp, 2, function(x) qnorm(rank(x, ties.method = "r") / (length(x) + 1)) )
}

## Calculate residual for model

exp_cov <- merge(exp, covs, by = 'row.names')
colnames(exp_cov)[1] <- "SUBJID"

expCovModel <- lm(as.formula(paste0(colnames(exp_cov)[2], " ~ ", paste(colnames(exp_cov)[-c(1:2)], collapse = " + "))), data = exp_cov)

exp_residual <- t(t(expCovModel$residuals))
rownames(exp_residual) <- exp_cov$SUBJID
colnames(exp_residual) <- "Residual"

if (shapiro.test(exp_residual)$p < 0.05) {
  exp_residual <- apply(exp_residual, 2, function(x) qnorm(rank(x, ties.method = "r") / (length(x) + 1)) )
}

# Extract genotypes from plink file to dosage file (500Kb)
cat(' INFO generating dosage genotype data ...\n')

dosagecmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_500K, ' --to-bp ', pos_to_500K, ' --recode A --out ', tmp_folder, '/', geneid)
system(dosagecmd, ignore.stdout = T, ignore.stderr = T)
bedcmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_500K, ' --to-bp ', pos_to_500K, ' --make-bed --out ', tmp_folder, '/', geneid)
system(bedcmd, ignore.stdout = T, ignore.stderr = T)

# Load dosage file (500Kb)
dosage_500K <- try(read.table(paste0(tmp_folder, '/', geneid, '.raw'), header = T, stringsAsFactors = F))

if ('try-error' %in% class(dosage_500K)) {
  stop(paste0('no SNP available for ', targetGene))
}

dosage_500K <- dosage_500K[, -c(1, 3:6)]
colnames(dosage_500K) <- c('SUBJID', sapply(colnames(dosage_500K)[-1], function(x) strsplit(x, "[_]")[[1]][1])) #rm the counted allele from rsnum
dosage_500K[, -1] <- round(apply(dosage_500K[, -1], 2, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)), 3) #post imputation imputation. NA replaced by mean

# Load Allele Info (500Kb)
snp_info_500K <- read.table(paste0(tmp_folder, '/', geneid, '.bim'), stringsAsFactors = F)
snp_info_500K$counted_allele <- snp_info_500K$V5
snp_info_500K$ref_allele <- snp_info_500K$V6
snp_info_500K$chr_bp <- paste0(snp_info_500K$V1, '_', snp_info_500K$V4)
colnames(snp_info_500K)[c(2,4)] <- c('rsid', 'bp')
snp_info_500K <- snp_info_500K[, c('rsid', 'chr_bp', 'bp', 'ref_allele', 'counted_allele')]

# TWAS model

#fit single tissue model to get proper window size and a lambda range
cat(' INFO fitting signle tissue prediction model \n')

residual_dosage <- merge(exp_residual, dosage_500K, by.x = 'row.names', by.y = 'SUBJID')

apa_fit <- cv.glmnet(x = as.matrix(residual_dosage[, -c(1:2)]), y = as.matrix(residual_dosage[, 2]), nfolds = 5, keep = T, alpha = 0.5)

fit.df <- data.frame(apa_fit$cvm, apa_fit$lambda, 1:length(apa_fit$cvm)) ##pull info to find best lambda
best.lam <- fit.df[which.min(fit.df[, 1]), ] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
cvm.best <- best.lam[, 1]
lambda.best <- best.lam[, 2]
nrow.best <- best.lam[, 3] ##position of best lambda in cv.glmnet output
ret <- as.data.frame(apa_fit$glmnet.fit$beta[, nrow.best]) # get betas from best lambda
ret[ret == 0.0] <- NA

bestbetas <- as.vector(ret[which(!is.na(ret)), ]) # vector of non-zero betas
names(bestbetas) <- rownames(ret)[which(!is.na(ret))]
bestbetas_snps <- bestbetas[which(names(bestbetas) %in% snp_info_500K$rsid)]

if (length(bestbetas_snps) > 0) {
  res <- cor.test(apa_fit$fit.preval[, nrow.best], residual_dosage[, 2])

  r.test <- res$estimate
  rsq <- r.test ^ 2
  pval <- res$p.value

  cat(' INFO r = ', r.test,', p = ', pval,' for ', targetGene, "\n")

  ### output best shrunken betas for PrediXcan
  bestbetalist <- names(bestbetas_snps)
  bestbetainfo <- snp_info_500K[which(snp_info_500K$rsid %in% bestbetalist), ]
  betatable <- as.matrix(cbind(bestbetainfo, bestbetas_snps))

  extra_res <- t(c(targetGene, strsplit(targetGene, "\\|")[[1]][2], r.test, rsq, length(bestbetas_snps), pval))
  colnames(extra_res) <- c("gene", "genename", "pred.perf.R", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval")
  write.table(extra_res, file = paste0(res_folder, "Extra_Komen_", geneid, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  weight_res <- cbind(targetGene, strsplit(targetGene, "\\|")[[1]][2], betatable[, -3], r.test, rsq, pval, lambda.best)
  colnames(weight_res) <- c("gene", "genename", "rsid", "chr_bp", "refAllele", "effectAllele", "weight", "R", "R2", "P", "Lambda")
  write.table(weight_res, file = paste0(res_folder, "Weight_Komen_", geneid, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  dosage <- dosage_500K[, names(bestbetas_snps)]
  cov <- cov(as.matrix(dosage))
  cov[upper.tri(cov)] <- NA
  cov <- cbind(gene = targetGene, expand.grid(dimnames(cov)), value = as.vector(cov))
  colnames(cov) <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write.table(cov, file = paste0(res_folder, "Cov_Komen_", geneid, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)
} else {
  cat(paste0(' INFO no SNP winner for ', targetGene))
}

#cleaning
cat(' INFO cleaning the tmp folder ... \n')
cmd = paste0('rm -r ', tmp_folder) #will only clean the subfolder under the tmp folder
system(cmd, wait = T)

cat(' INFO done \n')

