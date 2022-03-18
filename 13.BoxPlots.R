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


##########33
#
##################33practice


################final performance
###step1: append all results of all 26 cancer types
folder <- "../OutPut/EachCancerResults/"       
# path to folder that holds multiple .csv files 
setwd<-"F:/THESIS-CODE/DriveZT/OutPut/EachCancerResults/"
getwd()
file_list <- list.files(path=folder, pattern="*.csv")  
# create list of all .csv files in folder
####detect significant genes in each cancer type and set in in one file
siggenes<-list()
allcancers<-list()
cancers<-rep(0,length(file_list))
for(i in 1:length(file_list)){
  siggenes1<- as.data.frame(fread(file_list[[i]]))
  #head(siggenes1)
  siggenes[[i]]<- siggenes1$V1[which(siggenes1$PvalueNbModel<0.05)]
  cancers[i]<-unlist(strsplit(file_list,".csv")[[i]][2])
  allcancers[[i]]<-data.frame(Ctype=rep(cancers[i],length(siggenes[[i]])),Sigenes= siggenes[[i]])
}

cancerdd<-do.call(rbind, allcancers)

####step2: call merged cancer data and annotation(did in server) and mutation matrices
cancerdf<-as.data.frame(fread("../Processed_Data/cancerdf.csv"))
head(cancerdf)
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
mydata<-as.data.frame(fread("../Processed_Data/Siggenesdata.csv"))
head(mydata)
genes<-unique(mydata$genename)
cancers<-unique(cancerdf$cancer)
finalgenes<-list()
finaleachcancer<-list()
for(j in 2:length(cancers)){
  #j=1
  cancer<-cancerdf[which(cancerdf$cancer %in% cancers[j]),]
  genes<-unique(cancerdd$Sigenes[which(cancerdd$Ctype %in% cancers[j])])
  finalgenes<-list()
  for(i in 1:length(genes)){
    #i=1
    if(length(which(cancer$genename %in% genes[i]))>0){
    cancerf<-cancer[which(cancer$genename %in% genes[i]),]
    cancerf$muts<-0
    wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes[i]),]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is not mutated
    mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes[i])]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is mutated
    cancerf$muts[which(cancerf$CLEANNAME %in% wildtype )]<-0
    cancerf$muts[which(cancerf$CLEANNAME %in% mutant )]<-1
    #print(head(cancerf))
    finalgenes[[i]]<-cancerf
    }
  }
  finaleachcancer[[j]]<-do.call(rbind,finalgenes)
}

finalfile<-do.call(rbind,finaleachcancer)
head(finalfile)
table(finalfile$muts)
finalfile<-finalfile[,-1]
#############step 3 plot
pdf("allplots.pdf", onefile = TRUE)
p <- list()
library(ggplot2)
for(i in 2:length(cancers)){
  if(length(which(finalfile$cancer %in% cancers[i]))>0){
  cancerf<-finalfile[which(finalfile$cancer %in% cancers[i]),]
 
  
  cancerf$muts<-factor(cancerf$muts,labels=c("wildtype","mutant")) 
  p[[i]]<-ggplot(cancerf,aes(x=muts,y=rank,fill=muts))+geom_boxplot(alpha=0.3)+
    facet_wrap(~genename)+ggtitle(cancers[i])
}
}
dev.off()
p[[2]]
######
library(gridExtra)
pdf("plots.pdf", onefile = TRUE)
do.call(grid.arrange, p)
dev.off()
pp = p[1:4]
do.call(marrangeGrob, c(p, list(nrow = 4, ncol = 1))) 
pp<-list(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]])
pdf("plots.pdf", onefile = TRUE)
for (i in seq(length(p))) {
  if(length(p[[i]])>0){
  
  do.call(grid.arrange, p[[i]]) 
}
}
dev.off()
getwd()
