
setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\PLACO_FUMA_SMR")

rm(list=ls())
gc()
getwd()
.libPaths()

if (!requireNamespace("biomaRt", quietly = TRUE))BiocManager::install("biomaRt")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("GenomicRanges", quietly = TRUE))BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE))BiocManager::install("IRanges")
library(biomaRt)
library(org.Hs.eg.db)
library(GenomicRanges)
library(IRanges)

index <- data.table::fread("../export_index/gwasfile_info2.txt")

trait="finngen_R11_C3_BREAST_CANCER_QC_hg19"
inputFile <- index[index$name!=trait]$phenotype

alltissue=c("Whole_Blood","Breast_Mammary_Tissue","eQTLGen")

all_PLACO_FUMA_SMR <- data.frame()
for (phen in inputFile) {
  FUMA_path=paste0("../FUMA/",phen)

  GenomicRiskLoci=data.table::fread(paste0(FUMA_path,"/GenomicRiskLoci.txt"))

  geneManhattan=data.table::fread(paste0(FUMA_path,"/magma.genes.out"))
  geneManhattan$Trait_pair=paste0(index[index$name%in%trait,]$phenotype,"-",phen)
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
  geneManhattan$Pbon=p.adjust(p = geneManhattan$P, method = "bonferroni")
  geneManhattan_sig <- subset(geneManhattan,Pbon<0.05)

  data.table::fwrite(geneManhattan_sig,
                     file =paste0(phen,"_magma_gene_based_bon_sig.txt"),
                     sep = "\t",
                     na = "NA",
                     quote = FALSE)

  if (nrow(geneManhattan_sig)==0)next

  PLACO_data <- data.table::fread(paste0("../PLACO_pleiotropic/",phen,"_placo.txt"))
  geneManhattan_sig$sig.placo <- "NA"

  for (i in 1:nrow(geneManhattan_sig)) {
    Chr=geneManhattan_sig[i,]$CHR
    START=geneManhattan_sig[i,]$START
    STOP=geneManhattan_sig[i,]$STOP
    sub_PLACO_data <- subset(PLACO_data,CHR==Chr& BP>START & BP < STOP)
    geneManhattan_sig[i,]$sig.placo <- ifelse(any(PLACO_data$p.placo < 5e-8),1,0)
  }

  eqtl_data <- data.table::fread(paste0(FUMA_path,"/eqtl.txt"))
  eqtl_data$db_tissue <- paste0(eqtl_data$db,"_",eqtl_data$tissue)
  new_col_names <- paste0("eQTL_",unique(eqtl_data$db_tissue))

  new_data <- c(rep.int(0, nrow(geneManhattan_sig)))

  df <- data.frame(matrix(0,ncol = length(new_col_names), nrow = nrow(geneManhattan_sig)))
  colnames(df) <- new_col_names
  geneManhattan_sig <- cbind(geneManhattan_sig,df) %>% as.data.frame()
  
  eqtl_data_Overlaps <- eqtl_data[eqtl_data$gene %in% geneManhattan_sig$GENE,]
  
  for (i in 1:nrow(geneManhattan_sig)) {
    for (col in new_col_names) {
      if (geneManhattan_sig[i,]$GENE %in% eqtl_data_Overlaps$gene) {
        status2 <- 1
      }else{
        status2 <-0
      }
      
      j=  which(col ==colnames(geneManhattan_sig))
      geneManhattan_sig[i,j] <- status2
    }
  }
  
  new_data_SMR <- c(rep.int(0, nrow(geneManhattan_sig)))

  df <- data.frame(matrix(0,ncol = length(alltissue), nrow = nrow(geneManhattan_sig)))
  colnames(df) <- paste0("SMR_",alltissue)
  geneManhattan_sig <- cbind(geneManhattan_sig,df) %>% as.data.frame()
  for (i in 1:nrow(geneManhattan_sig)) {
    geneManhattan_sig[i,]$GENE
    
    for (tissue in alltissue) {
      SMR_file <- paste0("../SMR/result/",phen,"_",tissue,".smr")
      if (!file.exists(SMR_file))next
      SMR_data <- data.table::fread(SMR_file)
      SMR_data_sig <- subset(SMR_data,p_SMR< 0.05 & p_HEIDI > 0.05)
      if (geneManhattan_sig[i,]$GENE %in% SMR_data_sig$probeID) {
        status3 <- 1
      }else{
        status3 <-0
      }
      
      k =  which(paste0("SMR_",tissue) ==colnames(geneManhattan_sig))
      geneManhattan_sig[i,k] <- status3
    }
  }

  if (!"eQTL_eQTLGen_eQTLGen_trans_eQTLs" %in%names(geneManhattan_sig)) {
    geneManhattan_sig$eQTL_eQTLGen_eQTLGen_trans_eQTLs <- 0
  }
  
  all_PLACO_FUMA_SMR <- rbind(all_PLACO_FUMA_SMR,geneManhattan_sig)
}  

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- all_PLACO_FUMA_SMR$SYMBOL
gene_cytobands <-   getBM(
  attributes = c("hgnc_symbol", "band"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

gene_cytobands_filtered <- gene_cytobands[apply(gene_cytobands, 1, function(row) all(row != "")), ]

print(gene_cytobands_filtered)

all_PLACO_FUMA_SMR_merge <- merge(all_PLACO_FUMA_SMR,gene_cytobands_filtered,by.x = "SYMBOL", by.y = "hgnc_symbol", all.x = T)

data.table::fwrite(all_PLACO_FUMA_SMR_merge,
                   file ="all_PLACO_FUMA_SMR.txt",
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)