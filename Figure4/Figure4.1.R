rm(list = ls())
library(limma)
library(ggrepel)
library(ggthemes)
library(tidyverse)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 


FE<-read.table("SENESCENCE.txt",header=T,sep="\t")
FEgene<-FE[,1]

dat<-read.csv("merge matrix new.csv",header = T,sep = ",",row.names = 1)
group_list <- read.table("clinicaldata.txt", sep = "\t", header = T)
range(dat)
ccc=intersect(group_list$Accession,colnames(dat))
dat=dat[,ccc]
aaa<-intersect(rownames(dat),FEgene)
express.norm <-dat[aaa,]


x<-read.csv("result for DEGs.csv",header = T,row.names = 1)
x<-x[aaa,]
write.csv(x, "result for DEGs merge.csv", quote = F)

p.cut<-0.05
logFC.cut<-0.15


volcano<-x
volcano$type[(volcano$P.Value > p.cut|volcano$P.Value=="NA")|(volcano$logFC < logFC.cut)& volcano$logFC > -logFC.cut] <- "none significant"
volcano$type[volcano$P.Value <= p.cut & volcano$logFC >= logFC.cut] <- "up-regulated"
volcano$type[volcano$P.Value <= p.cut & volcano$logFC <= -logFC.cut] <- "down-regulated"  

p = ggplot(volcano,aes(logFC,-1*log10(P.Value),color=type))       
p + geom_point()

x_lim <- max(volcano$logFC,-x$logFC) 
gg=p + geom_point( aes(size = abs(logFC)),alpha = 0.4) + xlim(-x_lim,x_lim) +labs(x="log2(FC)",y="-log10(P.Value)")+
  scale_color_manual(values =c("blue","grey","red"))+
  geom_hline(aes(yintercept=-1*log10(p.cut)),colour="black", linetype="dashed") +
  geom_vline(xintercept=c(-logFC.cut,logFC.cut),colour="black", linetype="dashed")
print(gg) 

pdf(file="merge volcano.PDF", width=6, height=6)
gg
dev.off()



volcano1<-dplyr::filter(volcano,P.Value < 0.05)

genename<-read.table("DEGs.txt",header = T)
genename=genename$genename
express1<-as.data.frame(express.norm)
express1$ID<-rownames(express1)
phonedata<-dplyr::filter(express1,ID %in% genename)
phonedata<-dplyr::select(phonedata,-ID)
group_list<-group_list[order(group_list$condition,decreasing = T),]

group_list[,2]

annotation_col<-data.frame(group_list[,2])   
rownames(annotation_col)<-group_list[,1]
colnames(annotation_col)<-"Group"
phonedata<-select(phonedata,rownames(annotation_col))

pheatmap(phonedata,scale = "row",cluster_row=T,cluster_col=F, annotation_col = annotation_col,show_colnames = F
         ,border_color = F,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
pdf(file="merge heatmap.PDF", width=8, height=6)
pheatmap(phonedata,scale = "row",cluster_row=T,cluster_col=F, annotation_col = annotation_col,show_colnames = F
         ,border_color = F,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()



