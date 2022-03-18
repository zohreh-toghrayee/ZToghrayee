rm(list=ls())
library(dplyr)
##install.packages("stringr", repos='http://cran.us.r-project.org')
library("stringr")
library("dplyr")
library("tidyr")
library("lme4")
library(lmerTest)#
library("dplyr")
library("tidyr")
library("lme4")
library(lmerTest)
library(data.table)
#######################END OF MUTATIONAL MATRIX
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
sps1 = as.data.frame(fread("../Processed_Data/15.genesforspsplotshrna.csv"))

sps2<-readRDS("../Processed_Data/15.genesforspsplotdrive.RDS")
head(sps1)
dim(sps1)
dim(sps2)
head(sps2)
#################Adjust p-values
unique(sps1$GENESYMBOLS)
genes<-c("TP53","BRAF","NRAS","KRAS","PIK3CA")
spsf<-sps1[which(sps1$GENESYMBOLS %in% genes),]
#for(i in 1:length(genes)){
  corrfile<-merge(spsf,sps2,by="SEQ")
  dim(corrfile)
  unique(corrfile$GENESYMBOLS)
#}
  corrfile$logFC<-log((corrfile$SAMPLE_COUNT+0.5)/(corrfile$PLASMID_COUNT+1))
    aglogfc<-data.frame(aggregate(corrfile$logFC,by=list(corrfile$SEQ,corrfile$POOL),mean))
  agsps<-data.frame(aggregate(corrfile$SPS,by=list(corrfile$SEQ),mean))
  head(aglogfc)
  head(agsps)
  colnames(aglogfc)<-c("SEQ","pool","logFC")
  colnames(agsps)<-c("SEQ","SPS")
 filef<- merge(aglogfc,agsps,by="SEQ")
 head(filef)
 colnames(filef)<-c("SEQ","pool","logFC","SPS")
 myfile<-merge(filef,sps1,by="SEQ")
 head(myfile)
 dim(myfile)
  myfile$pool2<-factor(myfile$pool)
  table(myfile$pool2[which(myfile$GENESYMBOLS=="KRAS")])
  ggplot(myfile, aes(x=logFC, y=SPS.x,color=pool )) + 
    geom_point(shape=18)+
    geom_smooth(method=lm,  linetype="dashed",
                color="darkred", fill="red")+
    facet_wrap(~GENESYMBOLS)+ggtitle("Boxplot of Mean-based Ranks of Top Genes in Pancancers")
  corrgenes<-rep(0,length(genes))
  for(i in 1:genes[i]){
    genefile<-myfile[which(myfile$GENESYMBOLS %in% genes[i]),]
   corrgenes[i]<- cor.test(genefile$logFC.x,genefile$SPS.x)
  }
  
  
  library(GGally)
mycorfile<-myfile[,c("logFC","SPS.x","pool","GENESYMBOLS")]  
mycorfilegene<-mycorfile[which(mycorfile$GENESYMBOLS %in% "BRAF"),]
ggpairs(mycorfilegene, columns = 1:2, aes(color = pool, alpha = 0.5),
upper = list(continuous = wrap("smooth", size = 2.5)))+
  
  
  ggtitle("Plot of correlation of logFC and SPS of shRNAs in three POOLs in BRAF ")

mycorfilegene<-mycorfile[which(mycorfile$GENESYMBOLS %in% "BRAF"),]
ggpairs(mycorfilegene, columns = 1:2, aes(color = pool, alpha = 0.5),
        upper = list(combo = "facetdensity"))+
          ggtitle("Plot of correlation of logFC and SPS of shRNAs in three POOLs in BRAF ")