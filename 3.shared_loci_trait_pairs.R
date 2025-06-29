
rm(list=ls())
gc()
getwd()
.libPaths()

library(dplyr)

data <- data.table::fread("../predata_circular_diagram/Shared_pleiotropic_loci_idendified by_PLACO.txt")
data$trait_pairs <- paste0(data$phen1,"-",data$phen2)
data$Locus_boundary <- paste0(data$chr,":",data$start,"-",data$end)

data <- distinct(data, Region, trait_pairs, .keep_all = TRUE)

dup_loc <- unique(data$Region)
shared_loci <- data.frame()
for (i in dup_loc) {
  line <- data.frame(Number=which(i==dup_loc),
                     Locus_boundary=paste0(data[data$Region==i,]$Locus_boundary,collapse = ", "),
                     Region=unique(data[data$Region==i,]$Region),
                     Trait_pair=paste0(data[data$Region==i,]$trait_pairs,collapse = ", "))
  shared_loci <- rbind(shared_loci,line)
}

data.table::fwrite(shared_loci,
                   file = "shared_loci_trait_pairs.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")
