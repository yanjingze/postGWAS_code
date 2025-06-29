
setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\plot_magma_tissue_specific")

rm(list=ls())
gc()
getwd()
.libPaths()

library(ggplot2)
library(RColorBrewer)

width = 10
height = 14

MAGMA_data <- data.table::fread("../magma_tissue_specific/magma_tissue_specific_sig.txt")[,c("Trait_pair","FULL_NAME","BETA","SE","P")]
colnames(MAGMA_data) = c("exp","Tissues","BETA","SE","P")

dd = MAGMA_data
dd$Tissues = gsub("_"," ",dd$Tissues)
dd$exp=factor(dd$exp,levels = unique(dd$exp), 
              labels = unique(dd$exp))
dd$Tissues = factor(dd$Tissues)
dd$logP=(-1)*log10(dd$P)
p = ggplot(data=dd,
           aes(x = Tissues,y = logP,fill=logP))+
  geom_bar(stat="identity",color="grey",fill="#277DA1")+
  xlab('Tissues')+ ylab("-log10(P value)")+
  facet_grid(exp~.,scales = "free", space = "free") +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "blue") +   #-log10(0.05)显著性阈值
  geom_hline(yintercept = 3.033424, linetype = "dashed", color = "red") +   #-log10(0.05/54)显著性阈值
  #scale_fill_gradient(low = "lightblue", high = "darkblue") + 
  #scale_y_continuous(limits = c(0,1.5))+
  theme_bw() + 
  theme(legend.position = "None") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size=15),
        strip.text.x = element_text(margin = margin(b = 7),size = 15),#x label
        axis.text.x = element_text(margin = margin(l = 3),size = 15),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +  
  coord_flip()+ guides(shape = guide_legend(order = 0))
p

ggsave("plot_magma_tissue_specific.pdf",p,width = width,height = height)
