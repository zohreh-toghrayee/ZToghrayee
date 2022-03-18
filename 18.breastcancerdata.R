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
library(doParallel)
library(data.table)
myCluster <-makeCluster(24, type = "PSOCK")
registerDoParallel(myCluster)

##########################
annotation = read.csv("TableS2.csv", stringsAsFactors = F)
cancers2<-strsplit(annotation$PATHOLOGIST_ANNOTATION,"[:]")
annotation[,"cancer"]<-unlist(lapply(cancers2,function(x) {tolower(x[1])}))
cancers<-unique(annotation$cancer)
colnames(annotation)<-c("CLEANNAME","PRIMARY_SITE","PATHOLOGY_ANNOTATION","inCCLE","SOURCE","PRIDE","cancer")
##print(unique(annotation$CLEANNAME))
###############################
##print(colnames(annotation))

#######################
cancerd1 = as.data.frame(fread("4.finaldata.csv"))
#print(colnames(cancerd1))
###print(unique(cancerd$genename))
cancerd<-merge(cancerd1,annotation,by="CLEANNAME")
print(colnames(cancerd))


###############################MUTATIONAL MATRIX

wild_matix<-readRDS("wild_matix.R")
cancer_mutations3<-readRDS("cancer_mutations.R")
shRNAdata=as.data.frame(fread("shRNAdata.csv"))

###############well-known methods 
#####Define cancer specific
#
cancers<-c("leukemia","eye","soft_tissue","gastric","cns","kidney","thyroid","skin"                     
  ,"salivary_gland","ovary","lung","bone","endometrium","bladder","upper_aerodigestive_tract","pnet","breast","pancreas","colorectal",
  "prostate","undefined","liver","lymphoma","oesophagus","cervix","biliary_tract")  
library(glmmTMB)
finalcancer2<-list()
cancerd1<-cancerd[which(cancerd$cancer=="breast" ),]
mygenes2<-c("PAX5","RASGRP2","USP19","ROS1","HIPK3","NPFFR2","KRAS","BAZ2B","LYST","TINF2","NAV2","SLX4")
 
  cancerd<-cancerd1[which(cancerd1$genename %in% mygenes2),]
  
for(i in 1:length(mygenes2)){
  cancer<-cancerd[which(cancerd$genename %in% mygenes2[i]),]

  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes2[i])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  cancer$muts<-0
  cancer$muts[which(cancer$CLEANNAME %in% wildtype)]<-0
  cancer$muts[which(cancer$CLEANNAME %in% mutant)]<-1
  finalcancer2[[i]]<-cancer
}
##########append files of genes
cancerdf<-do.call(rbind, finalcancer2)

write.csv(cancerdf,"18.BreastCancerData.csv")
