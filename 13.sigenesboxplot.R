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
cancerd1 = as.data.frame(fread("../Processed_Data/4.finaldata.csv"))
siggenes<-as.data.frame(fread("../OutPut/finalsiggenes.csv"))
annotation = read.csv("../Input_Data/Annotation/TableS2.csv", stringsAsFactors = F)
cancers2<-strsplit(annotation$PATHOLOGIST_ANNOTATION,"[:]")
annotation[,"cancer"]<-unlist(lapply(cancers2,function(x) {tolower(x[1])}))
cancers<-unique(annotation$cancer)
colnames(annotation)<-c("CLEANNAME","PRIMARY_SITE","PATHOLOGY_ANNOTATION","inCCLE","SOURCE","PRIDE","cancer")
cancerd<-merge(cancerd1,annotation,by="CLEANNAME")
wild_matix<-readRDS("wild_matix.R")
cancer_mutations3<-readRDS("cancer_mutations.R")

#################

cancerd$muts<-0
finalcancer2<-list()

cancers<-unique(cancerd$cancer)
print(cancers)
for(j in 1:length(cancers)){
cancer<-cancerd[which(cancerd$cancer %in% cancers[j]),]
genes<-unique(siggenes$Sigenes[which(siggenes$Ctype %in% cancers[j])])
for(i in 1:length(genes)){
cancerf<-cancer[which(cancer$genename %in% genes[i]),]
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes[i]),]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes[i])]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is mutated
  cancerf$muts[which(cancerf$CLEANNAME %in% wildtype )]<-0
  cancerf$muts[which(cancerf$CLEANNAME %in% mutant )]<-1
finalcancer2[[i]]<-cancerf
}
}

##########append files of genes
cancerdf<-do.call(rbind, finalcancers2)
cancer6topg<-cancerdf[which(cancerdf$genename %in% c("KRAS","NRAS","PIK3CA","BRAF","CTNNB1","TP53",
                                                     "DHX9","TOP2A","PA2G4")),]
cancer6topg$mutation<-factor(cancer6topg$muts,labels=c("wildtype","mutant")) 
ggplot(cancer6topg,aes(x=mutation,y=rank,fill=mutation))+geom_boxplot(alpha=0.3)+
  facet_wrap(~genename)+ggtitle("Boxplot of Quantile1-based Ranks of Top Genes in Pancancers")
 

 write.csv(cancerdf,"Siggenesdata.csv")
