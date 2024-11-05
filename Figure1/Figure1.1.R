rm(list = ls())

library(limma)
library(sva)
library(tidyverse)
library(limma)
library(ggrepel)
library(ggthemes)
library(tidyverse)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 



express <- read.csv("GSE56815_matrix.csv", row.names = 1)
group_list <- read.table("GSE56815_sample.txt", sep = "\t", header = T)
express<-dplyr::select(express,group_list$Accession)


# 
# express$ID=rownames(express)
# library(data.table)
# library(org.Hs.eg.db)
# library(clusterProfiler)
# library(biomaRt)
# genename=rownames(express)
# gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="ENTREZID", columns=c("SYMBOL"))
# colnames(gene_map)[1]="ID"
# express=inner_join(gene_map,express,by="ID")
# express<-distinct(express,SYMBOL,.keep_all = T)
# express=na.omit(express)
# rownames(express)=express$SYMBOL
# express=express[,-1]
# express=express[,-1]


# express$ID=NA
# for (i in 1:nrow(express)) {
#   express$ID[i]<-unlist(str_split_fixed(rownames(express)[i],"//",3))[2]
# }
# express$ID<-str_replace_all(express$ID," ","")
# express<-distinct(express,ID,.keep_all = T)
# rownames(express)<-express$ID
# express<-dplyr::select(express,-ID)



express<-express[which(unlist(str_split_fixed(rownames(express),"///",2))[,2]==""),]
express=na.omit(express)
range(express)
# if yes no need to log transfer, if above this range, have to do log transfer.
# express <- log2(express+1)
# range(express)
express=rbind(geneNames=colnames(express), express)
write.table(express, file="GSE56815.txt", sep="\t", quote=F, col.names=F)

rm(list = ls())



files=c("GSE56814.txt", "GSE56815.txt")       


geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)


allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))

  
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] != "TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}


outTab=ComBat(allTab, batchType, par.prior=TRUE)

library(FactoMineR)
library(factoextra)
ddb.pca <- PCA(t(allTab), graph = FALSE)

pheno<-data.frame(ID=colnames(allTab))
pheno$cancer<-pheno$ID
pheno[1:73,2]<-"GSE56814"
pheno[74:153,2]<-"GSE56815"



pdf(file = "Before removing batch effects.pdf",height=8,width=8)
fviz_pca_ind(ddb.pca,
             geom.ind = "point", 
             pointsize =2, 
             pointshape = 21, 
             fill.ind = pheno$cancer, 
             palette = "lacent", 
             addEllipses = TRUE, 
             legend.title = "Groups", 
             title="") +
  theme_bw() + 
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )
dev.off()



library(FactoMineR)
library(factoextra)
ddb.pca <- PCA(t(outTab), graph = FALSE)
pdf(file = "After removing batch effects.pdf",height=8,width=8)
fviz_pca_ind(ddb.pca,
             geom.ind = "point", 
             pointsize =2, 
             pointshape = 21, 
             fill.ind = pheno$cancer, 
             palette = "lacent", 
             addEllipses = TRUE, 
             legend.title = "Groups", 
             title="") +
  theme_bw() + 
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )

dev.off()



clin=read.table("clinicaldata.txt",header = T)
clin=dplyr::filter(clin,condition == "OP")
clin=clin$Accession
colnames(outTab)=str_replace_all(colnames(outTab),"GSE56814_","")
colnames(outTab)=str_replace_all(colnames(outTab),"GSE56815_","")



outTab1=outTab[,clin]
outTab1=as.data.frame(outTab1)
outTab=as.data.frame(outTab)

save(outTab,file = "merge_all.RDATA")
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge_all.txt", sep="\t", quote=F, col.names=F)

save(outTab1,file = "merge_OP.RDATA")
outTab1=rbind(geneNames=colnames(outTab1), outTab1)
write.table(outTab1, file="merge_OP.txt", sep="\t", quote=F, col.names=F)

