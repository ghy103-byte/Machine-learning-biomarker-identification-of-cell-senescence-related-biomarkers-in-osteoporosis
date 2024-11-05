rm(list = ls())
df <- read.table("mergeROC.txt",head=T,sep="\t",check.names = F)
head(df)
df<-df[,-1]


library("pROC")

mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2","grey","green","red")


pdf("ROC.pdf",height=6,width=6)
auc.out <- c()

x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, 
              ci=TRUE, 
              main="",
              
              col=mycol[2],
              lwd=2, 
              legacy.axes=T)

ci.lower <- round(as.numeric(x$ci[1]),3) 
ci.upper <- round(as.numeric(x$ci[3]),3) 

auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)



for (i in 3:ncol(df)){
  x <- plot.roc(df[,1],df[,i],
                add=T, 
                smooth=F,
                ci=TRUE,
                col=mycol[i],
                lwd=2,
                legacy.axes=T)
  
  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)
  
  auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
}



auc.out <- as.data.frame(auc.out)
colnames(auc.out) <- c("Name","AUC","AUC CI")




legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:length(df)],
       lwd = 2,
       bty="n")
dev.off()



