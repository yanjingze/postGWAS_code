
setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\plot_S-LDSC_tissue_specific")

rm(list=ls())
gc()
getwd()
.libPaths()

library(ggplot2)
library(reshape)

index <- data.table::fread("../export_index/gwasfile_info2.txt")

width = 10
height = 12

all_ldsc_data <- data.frame()
all_plot_data <- list()
for (phen in index$phenotype) {
  ldsc_data <- data.table::fread(paste0("../S-LDSC_tissue_specific/",phen,"_GTEx_tissue_type_specific.cell_type_results_fdr.txt"))
  ldsc_data$Trait <- phen
  plot_data <- ldsc_data[,c("Name","Coefficient_P_value")]
  colnames(plot_data)[2]= phen
  all_plot_data[[phen]] <- plot_data
  all_ldsc_data <- rbind(all_ldsc_data,ldsc_data)
}

merged_plot_data <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), all_plot_data)

data.table::fwrite(all_ldsc_data,
                   file = "all_S_ldsc_data.txt",
                   sep = "\t")

data.table::fwrite(merged_plot_data,
                   file = "all_S_ldsc_data_Pval.txt",
                   sep = "\t")

plotdata = melt(merged_plot_data, id.var = c("Name"))
plotdata$variable = factor(plotdata$variable,levels = index$phenotype)
merged_plot_data$Name = factor(merged_plot_data$Name,levels = merged_plot_data$Name)

label_text = merged_plot_data$Name
theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x  = element_text(margin = margin(l = 3),size = 13,colour="black",angle = 45, hjust = 1),
  axis.text.y  = element_text(margin = margin(l = 3),size = 13,colour="black"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 13),
  legend.key.size=unit(0.7,'cm')
)

p <- ggplot(plotdata, aes(y = Name, x = variable, fill = value)) +
  geom_tile(color = "grey", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "#87CEFA", limits = c(0, 1), na.value = "white") +
  theme_change +
  scale_y_discrete(limits = rev(unique(label_text)), position = "left") +
  guides(fill = guide_legend(title = "Coefficient_P_value"))

print(p)

ggsave("plot_S-LDSC_tissue_specific.pdf",p,width = width,height = height)