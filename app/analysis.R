# ### analysis.R ###

library(trackViewer)
library(ggplot2)

File_location <- "."

vstNormalizedCounts_Macrophage <- read.table(paste0(File_location,"/","vstNormalizedCounts_Macrophage.txt"),header = T)
UCSC.hg19.genes<- read.table(paste0(File_location,"/","Gene_Symbols.txt"),header = F)

`%notin%` <- Negate(`%in%`)
