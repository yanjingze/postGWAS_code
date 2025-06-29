#!/bin/bash
"C:\ProgramData\miniconda3\envs\ldsc/python.exe" "C:/Program Files/R/R-4.3.3/library/bin/ldsc/ldsc.py" \
   --rg "../08.heritability_ldsc/result/BC.sumstats.gz","../08.heritability_ldsc/result/breast_adenosis.sumstats.gz","../08.heritability_ldsc/result/breast_fibrocystic_disease.sumstats.gz","../08.heritability_ldsc/result/intraductal_breast_neoplasm.sumstats.gz","../08.heritability_ldsc/result/mastitis.sumstats.gz","../08.heritability_ldsc/result/breast_cysts.sumstats.gz","../08.heritability_ldsc/result/LCIS.sumstats.gz"  \
   --ref-ld-chr  "../03.reference/ldsc/eur_w_ld_chr/"  \
   --w-ld-chr  "../03.reference/ldsc/eur_w_ld_chr/"  \
   --out "C:/Users/Administrator/Desktop/GWAS_POST/7.21/09.genetic_correlation_ldsc/result/BC_ldsc_rg"
