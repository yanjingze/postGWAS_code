
setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\plot_FUMA")


rm(list=ls())
gc()
getwd()
.libPaths()

if (!requireNamespace("qqman", quietly = TRUE))install.packages("qqman")
# if (!requireNamespace("topr", quietly = TRUE))install.packages("topr")
library(ggplot2)
library(qqman)

width_LocusZoom = 8
height_LocusZoom = 6

width=14
height=6

sig.thres = 5e-8

build=37

index <- data.table::fread("../export_index/gwasfile_info2.txt")

trait="finngen_R11_C3_BREAST_CANCER_QC_hg19"
inputFile <- index[index$name!=trait]$phenotype

GenomicRiskLoci_path="./GenomicRiskLoci/"
if (!dir.exists(GenomicRiskLoci_path)) {
  dir.create(GenomicRiskLoci_path)
}

for (phen in inputFile) {

  FUMA_path=paste0("../FUMA/",phen)
  
  GenomicRiskLoci=data.table::fread(paste0(FUMA_path,"/GenomicRiskLoci.txt"))
  message("\n当前第",which(phen == inputFile),"个GWAS:",phen)

  for (snp in GenomicRiskLoci$rsID) {
    i=which(snp == GenomicRiskLoci$rsID)
    message("当前第",i,"个GenomicRiskLoci,LeadSNPs：",GenomicRiskLoci[i,]$rsID)
    gwasfile <- paste0("../PLACO_pleiotropic/",phen,"_placo.txt")
    gwas_data=data.table::fread(gwasfile)
    gwas_data=transform(gwas_data,c=-qnorm(p.placo/2))

    gwas_data$z=ifelse(gwas_data$beta>0,gwas_data$c,-gwas_data$c)
    gwas_data=gwas_data[,c("SNP","CHR","BP","z")]
    colnames(gwas_data)=c("marker","chr","pos","z")
    head(gwas_data)

    filter <- GenomicRiskLoci[i,]
    chrom <- filter$chr
    endbp <- filter$end
    startbp <- filter$start
    markers <- subset(gwas_data,(chr==chrom & pos <= endbp & pos >= startbp))
    head(markers)  

    corr= ieugwasr::ld_matrix_local(markers$marker, 
                                    bfile="../03.reference/1000Gv3/EUR", 
                                    plink_bin=friendlypostGWASsoft::plink_binary(), 
                                    with_alleles = F)
    
    final_SNP=intersect(markers$marker, colnames(corr))
    SNP=snp
    if (!is.null(SNP)) {
      if (!SNP %in% final_SNP) {
        warning(snp," 不在重新构建的markers, corr数据集中，LocusZoom图不包含该位点")
        SNP=NULL
      }
    }

    plot <- gassocplot2::assoc_plot(markers[markers$marker%in%final_SNP,],
                                    build=build, 
                                    corr,top.marker=SNP,
                                    # label = add_SNP_legend,
                                    title = "LocusZoom plots of GWAS top SNPs",
                                    legend=TRUE,
                                    sig.thres=sig.thres)

    ggplot2::ggsave(paste0(GenomicRiskLoci_path,phen,"_",snp,"_locusZoom.pdf"), 
                    plot = plot,
                    width = width_LocusZoom, 
                    height = height_LocusZoom)
    
  }
}


geneManhattan_path="./geneManhattan/"
if (!dir.exists(geneManhattan_path)) {
  dir.create(geneManhattan_path)
}

for (phen in inputFile) {

  FUMA_path=paste0("../FUMA/",phen)

  GenomicRiskLoci=data.table::fread(paste0(FUMA_path,"/GenomicRiskLoci.txt"))
  geneManhattan=data.table::fread(paste0(FUMA_path,"/magma.genes.out"))
  geneManhattan=subset(geneManhattan,CHR!="X")
  geneManhattan$CHR=as.numeric(geneManhattan$CHR)
  
  geneManhattan$Pfdr=p.adjust(p = geneManhattan$P, method = "fdr")
  geneManhattan$Pbon=p.adjust(p = geneManhattan$P, method = "bonferroni")
  
  data.table::fwrite(geneManhattan,
                     file =paste0(geneManhattan_path,phen,"_magma_gene_based_fdr.txt"),
                     sep = "\t",
                     na = "NA",
                     quote = FALSE)

  pdf(file = paste0(geneManhattan_path,phen,"_geneManhattan.pdf"), width = width, height = height)
  manhattan(geneManhattan,
            chr = "CHR",
            bp = "START",
            p = "P",
            col = rainbow(11),
            snp = "SYMBOL",
            main = "",
            ylab = '-log10(P-value)',
            ylim = c(0, ifelse(10*ceiling(max(-log10(na.omit(geneManhattan$P)))/10)>300,300,10*ceiling(max(-log10(na.omit(geneManhattan$P)))/10))), #设定Y轴最高上限为300
            genomewideline = -log10(0.05 / nrow(geneManhattan)), 
            suggestiveline = -log10(0.05),
            annotatePval = 0.05 / nrow(geneManhattan),
            annotateTop = FALSE,
            logp = TRUE)
  dev.off()
}


snpManhattan_path="./snpManhattan/"
if (!dir.exists(snpManhattan_path)) {
  dir.create(snpManhattan_path)
}

for (phen in inputFile) {
  message("\n当前第",which(phen == inputFile),"个GWAS:",phen)
  snpfile <- paste0("../PLACO_pleiotropic/",phen,"_placo.txt")
  snpManhattan=data.table::fread(snpfile)
  snpManhattan=subset(snpManhattan,CHR!="X")
  snpManhattan$CHR=as.numeric(snpManhattan$CHR)

  png(filename =paste0(snpManhattan_path,phen,"_snpManhattan.png"), width = width*100, height = height*100, units = "px")
  manhattan(snpManhattan,
            chr = "CHR",
            bp = "BP",
            p = "p.placo",
            col = rainbow(11),
            snp = "SNP",
            main = "",
            ylab = '-log10(P-value)',
            ylim = c(0, ifelse(10*ceiling(max(-log10(na.omit(snpManhattan$pval)))/10)>300,300,10*ceiling(max(-log10(na.omit(snpManhattan$pval)))/10))), #设定Y轴最高上限为300
            genomewideline = -log10(sig.thres), 
            suggestiveline = -log10(0.05),
            # annotatePval = sig.thres,
            # annotateTop = FALSE,
            logp = TRUE)
  dev.off()
}
