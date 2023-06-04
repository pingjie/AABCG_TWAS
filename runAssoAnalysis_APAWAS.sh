#!/bin/sh

ml load GCC/8.2.0 OpenMPI/3.1.4 pandas/0.24.2 R/3.6.0

gwasDir=/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/sumstat

dbDir=/nobackup/sbcs/pingj2/AA_TWAS/models/R0.1/500K
outDir=/gpfs52/nobackup/sbcs/pingj2/AA_TWAS/res/R0.1/APAWAS

dbFile=${dbDir}/Komen.AFR.APAWAS.db
covFile=${dbDir}/Komen.AFR.APAWAS.Covs.txt

# Overall

/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${dbFile} \
  --covariance ${covFile} \
  --gwas_folder ${gwasDir} \
  --gwas_file_pattern "AFR.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${outDir}/Komen.AFR.APAWAS.Overall.results.csv 2>&1 | tee ${outDir}/Komen.AFR.APAWAS.Overall.results.log

Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${outDir}/Komen.AFR.APAWAS.Overall.results.csv ${outDir}/Komen.AFR.APAWAS.Overall.results.tsv

# ER Positive

/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${dbFile} \
  --covariance ${covFile} \
  --gwas_folder ${gwasDir}/ER_pos \
  --gwas_file_pattern "AFR_ER_pos.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${outDir}/Komen.AFR.APAWAS.ER_Pos.results.csv 2>&1 | tee ${outDir}/Komen.AFR.APAWAS.ER_Pos.results.log

Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${outDir}/Komen.AFR.APAWAS.ER_Pos.results.csv ${outDir}/Komen.AFR.APAWAS.ER_Pos.results.tsv

# ER Negative

/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${dbFile} \
  --covariance ${covFile} \
  --gwas_folder ${gwasDir}/ER_neg \
  --gwas_file_pattern "AFR_ER_neg.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${outDir}/Komen.AFR.APAWAS.ER_Neg.results.csv 2>&1 | tee ${outDir}/Komen.AFR.APAWAS.ER_Neg.results.log

Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${outDir}/Komen.AFR.APAWAS.ER_Neg.results.csv ${outDir}/Komen.AFR.APAWAS.ER_Neg.results.tsv

# TNBC

/nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
  --model_db_path ${dbFile} \
  --covariance ${covFile} \
  --gwas_folder ${gwasDir}/TNBC \
  --gwas_file_pattern "AFR_TNBC.SNPs.sumstat.withMimicID.txt" \
  --snp_column mimicID \
  --effect_allele_column TEST \
  --non_effect_allele_column OTHER \
  --beta_column BETA \
  --pvalue_column P \
  --output_file ${outDir}/Komen.AFR.APAWAS.TNBC.results.csv 2>&1 | tee ${outDir}/Komen.AFR.APAWAS.TNBC.results.log

Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${outDir}/Komen.AFR.APAWAS.TNBC.results.csv ${outDir}/Komen.AFR.APAWAS.TNBC.results.tsv
