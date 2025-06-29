
setwd("C:\Users\Administrator\Desktop\GWAS_POST\7.21\plot_circular_diagram")

rm(list=ls())
gc()
getwd()
.libPaths()

if (!requireNamespace("networkD3", quietly = TRUE))install.packages("networkD3")
if (!requireNamespace("webshot2", quietly = TRUE))install.packages("webshot2")
if (!requireNamespace("htmlwidgets", quietly = TRUE))install.packages("htmlwidgets")
library(webshot2)
library(networkD3)
library(htmlwidgets)

index <- data.table::fread("../export_index/gwasfile_info2.txt")

trait="finngen_R11_C3_BREAST_CANCER_QC_hg19"

data <- data.table::fread("../predata_circular_diagram/plot_circular_diagram.txt")

df <- data.frame(
  Level1=index[index$name%in%trait,]$phenotype,
  Level2 = data$phen2,
  Level3 = data$Region,
  Level4 = data$symbol,
  Level5 = data$PP.H4.abf
)

create_hierarchy <- function(df, levels) {
  if (length(levels) == 1) {
    return(unique(df[[levels]]))
  }
  
  lapply(unique(df[[levels[1]]]), function(x) {
    sub_df <- df[df[[levels[1]]] == x, ]
    list(name = x, children = create_hierarchy(sub_df, levels[-1]))
  })
}

hierarchy <- create_hierarchy(df, colnames(df))

print(hierarchy)
p=radialNetwork(List = hierarchy[[1]], fontSize = 12, opacity = 0.9, margin=0)
p

saveWidget(p, file="circus.html")

# webshot("circus.html", "circus.pdf")
webshot("circus.html", 
        "circus.pdf",
        vwidth = 1488,
        vheight = 1116)

# if (!requireNamespace("pagedown", quietly = TRUE))install.packages("pagedown")
# library(pagedown)
# pagedown::chrome_print("circus.html", output = "circus.pdf")