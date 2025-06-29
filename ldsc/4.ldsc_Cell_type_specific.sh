#!/bin/bash
"C:\ProgramData\miniconda3\envs\ldsc/python.exe" "C:/Program Files/R/R-4.3.3/library/bin/ldsc/ldsc.py" \
   --h2-cts "../08.heritability_ldsc/result/LCIS.sumstats.gz"  \
   --ref-ld-chr  "../03.reference/ldsc/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."  \
   --w-ld-chr  "../03.reference/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."  \
   --ref-ld-chr-cts  "GTEx.ldcts"  \
   --out "C:/Users/Administrator/Desktop/GWAS_POST/7.21/12.S-LDSC_tissue_specific/result/LCIS_GTEx_tissue_type_specific"
