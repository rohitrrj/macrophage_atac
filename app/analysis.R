# ### analysis.R ###

library(Gviz)
library(rtracklayer)
library(trackViewer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(session)
# library(limma)
# library(car)
# library(githubinstall)
library(FusionExpressionPlot)
library(ggplot2)
# library(DESeq2)

# HC_1<-"http://web.stanford.edu/~jadhav/Tracks/HC_1.bw"
# HC_2<-"http://web.stanford.edu/~jadhav/Tracks/HC_2.bw"
# HC_3<-"http://web.stanford.edu/~jadhav/Tracks/HC_3.bw"
# HC_4<-"http://web.stanford.edu/~jadhav/Tracks/HC_4.bw"
# CAD_1<-"http://web.stanford.edu/~jadhav/Tracks/CAD_1.bw"
# CAD_2<-"http://web.stanford.edu/~jadhav/Tracks/CAD_2.bw"
# CAD_3<-"http://web.stanford.edu/~jadhav/Tracks/CAD_3.bw"
# CAD_4<-"http://web.stanford.edu/~jadhav/Tracks/CAD_4.bw"

# HC_1<-"/mnt/Projects/Tuantuan/Tracks/HC_1.bw"
# HC_2<-"/mnt/Projects/Tuantuan/Tracks/HC_2.bw"
# HC_3<-"/mnt/Projects/Tuantuan/Tracks/HC_3.bw"
# HC_4<-"/mnt/Projects/Tuantuan/Tracks/HC_4.bw"
# CAD_1<-"/mnt/Projects/Tuantuan/Tracks/CAD_1.bw"
# CAD_2<-"/mnt/Projects/Tuantuan/Tracks/CAD_2.bw"
# CAD_3<-"/mnt/Projects/Tuantuan/Tracks/CAD_3.bw"
# CAD_4<-"/mnt/Projects/Tuantuan/Tracks/CAD_4.bw"

HC_1<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_1.bw"
HC_2<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_2.bw"
HC_3<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_3.bw"
HC_4<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_4.bw"
CAD_1<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_1.bw"
CAD_2<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_2.bw"
CAD_3<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_3.bw"
CAD_4<-"https://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_4.bw"

# fitqwr2.cqn<-readRDS("./Macrophage_fitqwr2cqn.Rds")
# tryCatch(
#   restore.session("./Tuantuan_Run2_Differential.Rda"),
#   error=function(e) e
# )

comparison <- "HCvsCAD"
File_location <- "."
ymax <- 8
xloc <- 6 # px
yloc <- 6 # py
# tops<-topTable(fitqwr2.cqn, coef = comparison, number = Inf,sort.by = "none")
# write.table(x = tops,file = "./HCvsCAD.txt",row.names = T,col.names = T,quote = F,sep = '\t')
tops <- read.table(paste0(File_location,"/",comparison,".txt"),header = T)
tops$Significance<-"NA"
tops[tops$adj.P.Val <= 0.05,"Significance"]<-"< 0.05"
tops[tops$adj.P.Val > 0.05,"Significance"]<-"> 0.05"
tops_sig<-subset(tops,adj.P.Val<0.05)
tops_sig_up<-subset(tops_sig,logFC>0)
tops_sig_down<-subset(tops_sig,logFC < 0)
tops_not_sig<-subset(tops,adj.P.Val>=0.05)
PvsU_up<-tops_sig_up
PvsU_down<-tops_sig_down
PvsU_not_sig<-tops_not_sig
tops_sig_mod<-tops_sig_up
tops_sig_mod<-rbind(tops_sig_mod,tops_sig_down)


