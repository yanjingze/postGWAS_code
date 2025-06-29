setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\predata_circular_diagram")

rm(list=ls())
gc()
getwd()
.libPaths()

if (!requireNamespace("NetWeaver", quietly = TRUE))install.packages("NetWeaver")
library(NetWeaver)

index <- data.table::fread("../export_index/gwasfile_info2.txt")

build=37

trait="finngen_R11_C3_BREAST_CANCER_QC_hg19"

coloc_result <- data.table::fread("../colocalization/coloc_result.txt")

data(ucsc.hg19.cytoband)
data(ucsc.hg38.cytoband)
if (build==37) {
  cytoband <- force(ucsc.hg19.cytoband)
}else if(build==38) {
  cytoband <- force(ucsc.hg38.cytoband)
}

all_Region <- c()
for (i in 1:nrow(coloc_result)) {
  Region <- as.character(subset(cytoband,Chr==paste0("chr",coloc_result[i,]$chr)& Start <=coloc_result[i,]$pos & End >=coloc_result[i,]$pos )$Band)
  all_Region <- c(all_Region,Region)
}

coloc_result$Region <- paste0(coloc_result$chr,all_Region)

coloc_result$Nearest_genes <- "NA"

phenotype <- index[!index$name%in%trait,]$phenotype
allmerge <- data.frame()
for (phen in phenotype) {
  map_gene <- data.table::fread(paste0("../FUMA/",phen,"/genes.txt"))
  map_gene$GenomicLocus <- as.numeric(map_gene$GenomicLocus)
  coloc_result_phen <- subset(coloc_result,phen2==phen)
  merge <- merge(coloc_result_phen,map_gene,by="GenomicLocus")
  allmerge <- rbind(allmerge,merge)
  allmerge$Nearest_genes <- allmerge$symbol
}


alldata <- data.frame()
for (phen in phenotype) {
  map_gene <- data.table::fread(paste0("../FUMA/",phen,"/genes.txt"))
  coloc_result_phen <- subset(coloc_result,phen2==phen)
  # merge_coloc_result_phen <- merge(coloc_result_phen,map_gene, by ="GenomicLocus",all.x =T)
  for (i in coloc_result_phen$GenomicLocus) {
    # if (! i%in% map_gene$GenomicLocus ) next
    coloc_result_phen[coloc_result_phen$GenomicLocus==i,]$Nearest_genes <- paste0(map_gene[map_gene$GenomicLocus==i,]$symbol,collapse = " ")
  }
  alldata <- rbind(coloc_result_phen,alldata)
}

data.table::fwrite(allmerge,
                   file = "plot_circular_diagram.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")

data.table::fwrite(alldata,
                   file = "Shared_pleiotropic_loci_idendified by_PLACO.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")