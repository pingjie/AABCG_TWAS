# module load PLINK/1.9b_5.2 GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0

library(glmnet)
set.seed(1024)

### config
args <- commandArgs(trailingOnly = TRUE)

targetIsoform <- args[1] ### test: isoform1

plink_path <- 'plink'

SNPs_in_sumstat <- "/gpfs52/nobackup/sbcs/pingj2/Komen/sumstat/AFR.snps"

workDir <- "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/models/spTWAS/"
setwd(workDir)

tmp_folder <- paste0(workDir, "/tmp/")
res_folder <- paste0(workDir, "/res/")

genotype_path <- "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/geno/Komen/Komen"

### Some functions ###
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# mkdir tmp folder. will be cleaned
cat(' INFO mkdir tmp folders ...\n')
options(warn = -1)
dir.create(tmp_folder)
dir.create(paste0(tmp_folder, '/', targetIsoform))
tmp_folder <- paste0(tmp_folder,'/', targetIsoform)
options(warn = 0)

# Load expression
cat(' INFO loading expression data ...\n')
psiTab <- loadRData("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/exp/Komen.AFR.PSI.RMA.PC.lincRNA.Rdata")

psi <- as.data.frame(t(psiTab[targetIsoform, -c(1:5)]))
colnames(psi) <- "Isoform"

## Calculate residual for model
peerFactors <- loadRData("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/exp/Komen.AFR.PEER_factors.rda")

psi_cov <- merge(psi, peerFactors, by = 'row.names')
colnames(psi_cov)[1] <- "SUBJID"

psiCovModel <- lm(as.formula(paste0(colnames(psi_cov)[2], " ~ ", paste(colnames(psi_cov)[-c(1:2)], collapse = " + "))), data = psi_cov)

psi_residual <- t(t(psiCovModel$residuals))
rownames(psi_residual) <- psi_cov$SUBJID
colnames(psi_residual) <- "Residual"

if (shapiro.test(psi_residual)$p < 0.05) {
  psi_residual <- apply(psi_residual, 2, function(x) qnorm(rank(x, ties.method = "r") / (length(x) + 1)) )
}

ensemblID <- psiTab[targetIsoform, "ENSEMBLID"]

# Get chr pos for the gene
dist <- 500000
chr <- as.numeric(sub('^chr', '', psiTab[targetIsoform, "chr"]))
pos_from <- as.numeric(psiTab[targetIsoform, "start"])
pos_to <- as.numeric(psiTab[targetIsoform, "end"])
pos_from_gene <- max(1, pos_from - dist)
pos_to_gene <- pos_to + dist

# Extract genotypes from plink file to dosage file (50Kb)
cat(' INFO generating dosage genotype data ...\n')

dosagecmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_gene, ' --to-bp ', pos_to_gene, ' --recode A --out ', tmp_folder, '/', targetIsoform)
system(dosagecmd, ignore.stdout = T, ignore.stderr = T)
bedcmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_gene, ' --to-bp ', pos_to_gene, ' --make-bed --out ', tmp_folder, '/', targetIsoform)
system(bedcmd, ignore.stdout = T, ignore.stderr = T)

# Load dosage file (50Kb)
dosage_dist <- try(read.table(paste0(tmp_folder, '/', strsplit(targetIsoform, "@")[[1]][1], '.raw'), header = T, stringsAsFactors = F))

if ('try-error' %in% class(dosage_dist)) {
  stop(paste0('no SNP available for ', targetIsoform))
}

dosage_dist <- dosage_dist[, -c(1, 3:6)]
colnames(dosage_dist) <- c('SUBJID', sapply(colnames(dosage_dist)[-1], function(x) strsplit(x, "[_]")[[1]][1])) #rm the counted allele from rsnum
dosage_dist[, -1] <- round(apply(dosage_dist[, -1], 2, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)), 3) #post imputation imputation. NA replaced by mean

# Load Allele Info (500Kb)
snp_info_dist <- read.table(paste0(tmp_folder, '/', targetIsoform, '.bim'), stringsAsFactors = F)
snp_info_dist$counted_allele <- snp_info_dist$V5
snp_info_dist$ref_allele <- snp_info_dist$V6
snp_info_dist$chr_bp <- paste0(snp_info_dist$V1, '_', snp_info_dist$V4)
colnames(snp_info_dist)[c(2,4)] <- c('rsid', 'bp')
snp_info_dist <- snp_info_dist[, c('rsid', 'chr_bp', 'bp', 'ref_allele', 'counted_allele')]

# spTWAS model

#fit single tissue model to get proper window size and a lambda range
cat(' INFO fitting signle tissue prediction model \n')

residual_dosage <- merge(psi_residual, dosage_dist, by.x = 'row.names', by.y = 'SUBJID')

psi_fit <- cv.glmnet(x = as.matrix(residual_dosage[, -c(1:2)]), y = as.matrix(residual_dosage[, 2]), nfolds = 5, keep = T, alpha = 0.5)

fit.df <- data.frame(psi_fit$cvm, psi_fit$lambda, 1:length(psi_fit$cvm)) ##pull info to find best lambda
best.lam <- fit.df[which.min(fit.df[, 1]), ] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
cvm.best <- best.lam[, 1]
lambda.best <- best.lam[, 2]
nrow.best <- best.lam[, 3] ##position of best lambda in cv.glmnet output
ret <- as.data.frame(psi_fit$glmnet.fit$beta[, nrow.best]) # get betas from best lambda
ret[ret == 0.0] <- NA

bestbetas <- as.vector(ret[which(!is.na(ret)), ]) # vector of non-zero betas
names(bestbetas) <- rownames(ret)[which(!is.na(ret))]
bestbetas_snps <- bestbetas[which(names(bestbetas) %in% snp_info_dist$rsid)]

if (length(bestbetas_snps) > 0) {
  res <- cor.test(psi_fit$fit.preval[, nrow.best], residual_dosage[, 2])

  r.test <- res$estimate
  rsq <- r.test ^ 2
  pval <- res$p.value

  cat(' INFO r = ', r.test,', p = ', pval,' for ', targetIsoform, "\n")

  ### output best shrunken betas for spTWAS
  bestbetalist <- names(bestbetas_snps)
  bestbetainfo <- snp_info_dist[which(snp_info_dist$rsid %in% bestbetalist), ]
  betatable <- as.matrix(cbind(bestbetainfo, bestbetas_snps))

  extra_res <- t(c(targetIsoform, ensemblID, r.test, rsq, length(bestbetas_snps), pval))
  colnames(extra_res) <- c("psiName", "genename", "pred.perf.R", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval")
  write.table(extra_res, file = paste0(res_folder, "Extra_Komen_", targetIsoform, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  weight_res <- cbind(targetIsoform, ensemblID, betatable[, -3], r.test, rsq, pval, lambda.best)
  colnames(weight_res) <- c("psiName", "genename", "rsid", "chr_bp", "refAllele", "effectAllele", "weight", "R", "R2", "P", "Lambda")
  write.table(weight_res, file = paste0(res_folder, "Weight_Komen_", targetIsoform, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  dosage <- dosage_dist[, names(bestbetas_snps)]
  cov <- cov(as.matrix(dosage))
  cov[upper.tri(cov)] <- NA
  cov <- cbind(gene = targetIsoform, expand.grid(dimnames(cov)), value = as.vector(cov))
  colnames(cov) <- c('psiName', 'RSID1', 'RSID2', 'VALUE')
  write.table(cov, file = paste0(res_folder, "Cov_Komen_", targetIsoform, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)
} else {
  cat(paste0(' INFO no SNP winner for ', targetIsoform))
}

#cleaning
cat(' INFO cleaning the tmp folder ... \n')
cmd = paste0('rm -r ', tmp_folder) #will only clean the subfolder under the tmp folder
system(cmd, wait = T)

cat(' INFO done \n')

