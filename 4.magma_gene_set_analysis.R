setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\magma_gene_set_analysis")

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
  map_gene <- read.table(paste0("../FUMA/",phen,"/magma.gsa.out"),header = T)
  map_gene$Trait_pair <- paste0(index[index$name%in%trait,]$phenotype,"-",phen)
  map_gene$Pbon=p.adjust(p = map_gene$P, method = "bonferroni")
  alldata <- rbind(alldata,map_gene)
  alldata <- alldata[order(alldata$P),]
}

data.table::fwrite(alldata,
                   file = "magma_gene_set_analysis.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")

alldata_sig <- subset(alldata,Pbon < 0.05)

if (nrow(alldata_sig)<10) {
  all_trait_alldata_five <- data.frame()
  Trait_pair <- unique(alldata$Trait_pair)
  for (pair in Trait_pair) {
    trait_alldata <- alldata[alldata$Trait_pair==pair,]
    trait_alldata_five <- trait_alldata[c(1:5),]
    all_trait_alldata_five <- rbind(trait_alldata_five,all_trait_alldata_five)
  }
  alldata_sig <- all_trait_alldata_five
}

data.table::fwrite(alldata_sig,
                   file = "magma_gene_set_analysis_sig.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")