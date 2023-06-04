#!/bin/sh

ml load GCC/8.2.0 OpenMPI/3.1.4 pandas/0.24.2 R/3.6.0

dbDir=/nobackup/sbcs/pingj2/AA_TWAS/models/R0.1/500K
resDir=/nobackup/sbcs/pingj2/AA_TWAS/res/R0.1/spTWAS

twasDB=${dbDir}/Komen.AFR.spTWAS.db
twasCov=${dbDir}/Komen.AFR.spTWAS.Covs.txt

### Overall
/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${twasDB} \
  --covariance ${twasCov} \
  --gwas_folder /nobackup/sbcs/pingj2/AA_TWAS/sumstat \
  --gwas_file_pattern "AFR.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${resDir}/Komen.AFR.spTWAS.Overall.results.csv 2>&1 | tee ${resDir}/Komen.AFR.spTWAS.Overall.results.log

### ER Positive
/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${twasDB} \
  --covariance ${twasCov} \
  --gwas_folder /nobackup/sbcs/pingj2/AA_TWAS/sumstat/ER_pos \
  --gwas_file_pattern "AFR_ER_pos.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${resDir}/Komen.AFR.spTWAS.ER_Pos.results.csv 2>&1 | tee ${resDir}/Komen.AFR.spTWAS.ER_Pos.results.log

### ER Negative
/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${twasDB} \
  --covariance ${twasCov} \
  --gwas_folder /nobackup/sbcs/pingj2/AA_TWAS/sumstat/ER_neg \
  --gwas_file_pattern "AFR_ER_neg.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${resDir}/Komen.AFR.spTWAS.ER_Neg.results.csv 2>&1 | tee ${resDir}/Komen.AFR.spTWAS.ER_Neg.results.log

### TNBC
/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${twasDB} \
  --covariance ${twasCov} \
  --gwas_folder /nobackup/sbcs/pingj2/AA_TWAS/sumstat/TNBC \
  --gwas_file_pattern "AFR_TNBC.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${resDir}/Komen.AFR.spTWAS.TNBC.results.csv 2>&1 | tee ${resDir}/Komen.AFR.spTWAS.TNBC.results.log

Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/Komen.AFR.spTWAS.Overall.results.csv ${resDir}/Komen.AFR.spTWAS.Overall.results.tsv
Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/Komen.AFR.spTWAS.ER_Pos.results.csv ${resDir}/Komen.AFR.spTWAS.ER_Pos.results.tsv
Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/Komen.AFR.spTWAS.ER_Neg.results.csv ${resDir}/Komen.AFR.spTWAS.ER_Neg.results.tsv
Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/Komen.AFR.spTWAS.TNBC.results.csv ${resDir}/Komen.AFR.spTWAS.TNBC.results.tsv
