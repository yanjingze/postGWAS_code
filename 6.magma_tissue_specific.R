
setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\magma_tissue_specific")

rm(list=ls())
gc()
getwd()
.libPaths()

library(data.table)

index <- data.table::fread("../export_index/gwasfile_info2.txt")

trait="finngen_R11_C3_BREAST_CANCER_QC_hg19"

phenotype <- index[!index$name%in%trait,]$phenotype
alldata <- data.frame()
for (phen in phenotype) {
  map_gene <- read.table(paste0("../FUMA/",phen,"/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"),header = T)
  map_gene$Trait_pair <- paste0(index[index$name%in%trait,]$phenotype,"-",phen)
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
  map_gene$Pbon=p.adjust(p = map_gene$P, method = "fdr")
  # alldata$Pfdr=p.adjust(p = alldata$P, method = "fdr")
  alldata <- rbind(alldata,map_gene)
  alldata <- alldata[order(alldata$P),]
}

data.table::fwrite(alldata,
                   file = "magma_tissue_specific.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")

alldata_sig <- subset(alldata,P < 0.05)

data.table::fwrite(alldata_sig,
                   file = "magma_tissue_specific_sig.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")
