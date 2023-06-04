# module load PLINK/1.9b_5.2 GCC/10.2.0 OpenMPI/4.0.5 R/4.0.5

library(glmnet)
set.seed(1024)

### config
args <- commandArgs(trailingOnly = TRUE)

gene <- args[1] ### test: ENSG00000007237.18

plink_path <- 'plink'

SNPs_in_sumstat <- "/gpfs52/nobackup/sbcs/pingj2/sumstat/breast/AFR/AFR.variants"

workDir <- "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/models"
setwd(workDir)

tmp_folder <- paste0(workDir, "/tmp/")
res_folder <- paste0(workDir, "/res/")

genotype_path <- "/gpfs52/nobackup/sbcs/pingj2/Komen/geno/Komen/Komen.AFR"

### Some functions ###
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load gencode annotation file
cat(' INFO loading gencode annotation ...\n')
gencode_path <- "/gpfs52/nobackup/sbcs/pingj2/db/gencode.v26.GRCh38.genes.txt"
gencode <- read.table(gencode_path, header = T, stringsAsFactors = F)

# mkdir tmp folder. will be cleaned
cat(' INFO mkdir tmp folders ...\n')
options(warn = -1)
dir.create(tmp_folder)
dir.create(res_folder)

dir.create(paste0(tmp_folder, '/', gene))
tmp_folder <- paste0(tmp_folder,'/', gene)
options(warn = 0)

# Get chr pos for the gene
dist <- 500000
chr <- as.numeric(sub('^...', '', gencode$chr[which(gencode$geneid_full == gene)]))
pos_from <- gencode$left[which(gencode$geneid_full == gene)]
pos_to <- gencode$right[which(gencode$geneid_full == gene)]
pos_from_gene <- max(1, pos_from - dist)
pos_to_gene <- pos_to + dist

# Load expression
cat(' INFO loading expression data ...\n')
exptab <- loadRData("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.V26.PEER_residuals.rda")
exp.all <- as.data.frame(t(t(exptab[, gene])))
colnames(exp.all) <- "exp"

# Extract genotypes from plink file to dosage file
cat(' INFO generating dosage genotype data ...\n')

dosagecmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_gene, ' --to-bp ', pos_to_gene, ' --recode A --out ', tmp_folder, '/', gene)
system(dosagecmd, ignore.stdout = T, ignore.stderr = T)
bedcmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, ' --chr ', chr, ' --from-bp ', pos_from_gene, ' --to-bp ', pos_to_gene, ' --make-bed --out ', tmp_folder, '/', gene)
system(bedcmd, ignore.stdout = T, ignore.stderr = T)

# Load dosage file
dosage_from_gene <- try(read.table(paste0(tmp_folder, '/', gene, '.raw'), header = T, stringsAsFactors = F))

if ('try-error' %in% class(dosage_from_gene)) {
  stop(paste0('no SNP available for ', gene))
}

dosage_from_gene <- dosage_from_gene[, -c(1, 3:6)]
names(dosage_from_gene) <- c('SUBJID', gsub('.{2}$', '', sub('.', '', names(dosage_from_gene)[-1])))
dosage_from_gene[, -1] <- round(apply(dosage_from_gene[, -1], 2, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)), 3) #post imputation imputation. NA replaced by mean

# Load Allele Info
snp_info_from_gene <- read.table(paste0(tmp_folder, '/', gene, '.bim'), stringsAsFactors = F)
snp_info_from_gene$counted_allele <- snp_info_from_gene$V5
snp_info_from_gene$ref_allele <- snp_info_from_gene$V6
snp_info_from_gene$chr_bp <- paste0(snp_info_from_gene$V1, '_', snp_info_from_gene$V4)
colnames(snp_info_from_gene)[c(2,4)] <- c('rsid', 'bp')
snp_info_from_gene <- snp_info_from_gene[, c('rsid', 'chr_bp', 'bp', 'ref_allele', 'counted_allele')]

# PrediXcan model

#fit single tissue model to get proper window size and a lambda range
cat(' INFO fitting signle tissue prediction model \n')

exp_dosage <- merge(exp.all, dosage_from_gene, by.x = "row.names", by.y = 'SUBJID')

fit <- cv.glmnet(x = as.matrix(exp_dosage[, -c(1:2)]), y = as.matrix(exp_dosage[, 2]), nfolds = 5, keep = T, alpha = 0.5)

#fit <- cv.glmnet(x = as.matrix(exp_dosage[, -c(1:2)]), y = as.matrix(exp_dosage[, 2]), nfolds = 5, keep = T, alpha = 0.5, nlambda = 50, pmax = 200)

fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm)) ##pull info to find best lambda
best.lam <- fit.df[which.min(fit.df[, 1]), ] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
cvm.best <- best.lam[, 1]
lambda.best <- best.lam[, 2]
nrow.best <- best.lam[, 3] ##position of best lambda in cv.glmnet output
ret <- as.data.frame(fit$glmnet.fit$beta[, nrow.best]) # get betas from best lambda
ret[ret == 0.0] <- NA

bestbetas <- as.vector(ret[which(!is.na(ret)), ]) # vector of non-zero betas
names(bestbetas) <- rownames(ret)[which(!is.na(ret))]
bestbetas_snps <- bestbetas[which(names(bestbetas) %in% snp_info_from_gene$rsid)]

if (length(bestbetas_snps) > 0) {
  res <- cor.test(fit$fit.preval[, nrow.best], exp_dosage[, 2])

  r.test <- res$estimate
  rsq <- r.test ^ 2
  pval <- res$p.value

  cat(' INFO r = ', r.test,', p = ', pval,' for ', gene, "\n")

  ### output best shrunken betas for PrediXcan
  bestbetalist <- names(bestbetas_snps)
  bestbetainfo <- snp_info_from_gene[which(snp_info_from_gene$rsid %in% bestbetalist), ]
  betatable <- as.matrix(cbind(bestbetainfo, bestbetas_snps))

  extra_res <- t(c(gene, gencode$genename[which(gencode$geneid_full == gene)], r.test, rsq, length(bestbetas_snps), pval))
  colnames(extra_res) <- c("gene", "genename", "pred.perf.R", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval")
  write.table(extra_res, file = paste0(res_folder, "Extra_AFR_Komen_V26_", gene, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  weight_res <- cbind(gene, gencode$genename[which(gencode$geneid_full == gene)], betatable[, -3], r.test, rsq, pval, lambda.best)
  colnames(weight_res) <- c("gene", "genename", "rsid", "chr_bp", "refAllele", "effectAllele", "weight", "R", "R2", "P", "Lambda")
  write.table(weight_res, file = paste0(res_folder, "Weight_AFR_Komen_V26_", gene, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

  dosage <- dosage_from_gene[, names(bestbetas_snps)]
  cov <- cov(as.matrix(dosage))
  cov[upper.tri(cov)] <- NA
  cov <- cbind(gene = gene, expand.grid(dimnames(cov)), value = as.vector(cov))
  colnames(cov) <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write.table(cov, file = paste0(res_folder, "Cov_AFR_Komen_V26_", gene, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)
} else {
  cat(paste0(' INFO no SNP winner for ', gene))
}

#cleaning
cat(' INFO cleaning the tmp folder ... \n')
cmd = paste0('rm -r ', tmp_folder) #will only clean the subfolder under the tmp folder
system(cmd, wait = T)

cat(' INFO done \n')

