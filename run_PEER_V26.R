library(preprocessCore)

expTab <- read.table("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.RNAseQC.TPM.gencode_v26.txt", sep = "\t", header = T, row.names = 1, as.is = T)

gencode <- read.table("/gpfs52/nobackup/sbcs/pingj2/db/gencode.v26.GRCh38.genes.txt", sep = "\t", header = T, as.is = T)
gencode_auto <- subset(gencode, !(chr %in% c("chrX", "chrY", "chrM")))

### Remove genes with no expression data and genes on chrX and chrY (n=10,424)
expTab_rowSums <- rowSums(expTab)
expTab <- expTab[expTab_rowSums > 0, ]
expTab <- expTab[which(rownames(expTab) %in% gencode_auto$geneid_full), ]

expTab_log_QN <- normalize.quantiles(as.matrix(log2(expTab + 1)))
rownames(expTab_log_QN) <- rownames(expTab)
colnames(expTab_log_QN) <- colnames(expTab)

expTab_log_QN_Inverse <- apply(expTab_log_QN, 1, function(x) qnorm(rank(x, ties.method = "r") / (length(x) + 1)) )

exp <- as.data.frame(cbind(rownames(expTab_log_QN_Inverse), expTab_log_QN_Inverse))
colnames(exp)[1] <- "ID"

write.table(exp, file = "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.TPM.gencode_v26.noZero.autosome.log2.QN.Inversed.txt", sep = "\t", row.names = F, quote = F)


## singularity exec /gpfs23/scratch/sbcs/pingj2/r-peer_1.3--r341h470a237_1.sif R

set.seed(1024)
library(peer)

expTab <- read.table("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.TPM.gencode_v26.noZero.autosome.log2.QN.Inversed.txt", sep = "\t", header = T, row.names = 1, as.is = T)
covTab <- read.table("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/clin/Komen.AFR.cov.txt", sep = "\t", header = T, row.names = 1, as.is = T)
covTab <- covTab[rownames(expTab), ]

model <- PEER()
PEER_setPhenoMean(model, as.matrix(expTab))
PEER_setCovariates(model, as.matrix(covTab))

peerN <- 15

# the num here need to be decided per no. of subjects in the tissue of interest, see https://www.gtexportal.org/home/documentationPage
# for eQTL: 15 factors for N<150, 30 factors for 150<= N <250, 45 factors for 250<= N <350, and 60 factors for N>=350
# as a result of optimizing for the number of eGenes discovered
if ( nrow(PEER_getPhenoMean(model)) < 150 ) {
    peerN <- 15
} else if ( nrow(PEER_getPhenoMean(model)) < 250 & nrow(PEER_getPhenoMean(model)) >= 150 ) {
    peerN <- 30
} else if ( nrow(PEER_getPhenoMean(model)) < 350 & nrow(PEER_getPhenoMean(model)) >= 250 ) {
    peerN <- 45
} else if ( nrow(PEER_getPhenoMean(model)) >= 350 ) {
    peerN <- 60
}

PEER_setNk(model, peerN)
PEER_update(model)  # 100+ iterations; take 50+ hours

factors <- PEER_getX(model)
rownames(factors) <- rownames(expTab)
colnames(factors) <- c(colnames(covTab), paste0("PEER_", 1:peerN))

residuals <- PEER_getResiduals(model)  # This is the Residuals after adjusting for PEER factors and PC variables, which are to be used for following analyses

save(factors, file = "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.V26.PEER_factors.rda")
write.table(factors, file = "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.V26.PEER_factors.txt", sep = "\t", row.names = T, quote = F)

pdf("/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.V26.PEER.diagnostics.pdf")
PEER_plotModel(model)
dev.off()

rownames(residuals) <- rownames(expTab)
colnames(residuals) <- colnames(expTab)
save(residuals, file = "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.V26.PEER_residuals.rda")
write.table(residuals, file = "/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/V26/exp/Komen.AFR.V26.PEER_residuals.txt", sep = "\t", row.names = T, quote = F)
