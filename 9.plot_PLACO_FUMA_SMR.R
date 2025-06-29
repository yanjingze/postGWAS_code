setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\plot_PLACO_FUMA_SMR")

rm(list=ls())
gc()
getwd()
.libPaths()

if (!requireNamespace("reshape", quietly = TRUE))install.packages("reshape")
library(dplyr)
library(reshape)
library(ggplot2)

width = 10
height = 12

cc = data.table::fread("../PLACO_FUMA_SMR/all_PLACO_FUMA_SMR.txt")
dd = subset(cc,select=-c(GENE,CHR,START,STOP,NSNPS,NPARAM,N,ZSTAT,P,Trait_pair,Pbon,band))

dd <- dd %>%
  group_by(SYMBOL) %>%
  summarize(across(everything(), sum, na.rm = TRUE))

colnames(dd)[2] = c("Loci_signal")

ee = data.frame(dd)

if (F) {
c0 = dd$Loci_signal
c0[which(c0>1)] = 1
c1 = apply(ee%>% select(starts_with("eQTL_")),1,sum)
c1[which(c1>1)] = 1
c2 = apply(ee%>% select(starts_with("SMR_")),1,sum)
c2[which(c2>1)] = 1
c3 = c0 + c1 + c2
ind = which(c3 == 3)
ee = ee[ind,]
}

ff = melt(ee, id.var = c("SYMBOL"))
ff$variable = factor(gsub("\\.","/",ff$variable),levels = colnames(dd)[-c(1)])
ee$SYMBOL = factor(ee$SYMBOL,levels = ee$SYMBOL)

label_text = ee$SYMBOL
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

## output the graphics
p = ggplot(ff, aes(y = SYMBOL, x = variable, fill = value)) +
  geom_tile(color = "grey", size = 0.5) +
  scale_fill_gradient(low = "white", high = "#87CEFA", na.value = "white", breaks = seq(0, max(ff$value), by = 1)) +
  theme_change +
  scale_y_discrete(limits = rev(unique(label_text)), position = "left") +
  guides(fill = guide_legend(title = "")) #+theme(legend.position="top")

print(p)

ggsave("plot_PLACO_FUMA_SMR.pdf",p,width = width,height = height)