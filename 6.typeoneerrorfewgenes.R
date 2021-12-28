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
cancerd = read.csv("../Processed_Data/finaldata.csv")
mygenes2<-unique(cancerd$genename)

wild_matix<-readRDS("../wild_matix.R")
cancer_mutations3<-readRDS("../cancer_mutations.R")
################################################3
vvm<-rep(0,length(mygenes2))
names(vvm)<-mygenes2
n<-1000
library("TMB")
library(glmmTMB)
i=which(mygenes2=="KRAS")
for(i in 1:length(mygenes2)){
  cancer<-cancerd[which(cancerd$genename %in% mygenes2[i]),]
  cancelf<-unique(cancer$CLEANNAME)
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)=="KRAS")]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  others<-intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==2)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  celinesf<-cancelf[-which(cancelf %in% others)]
  cancer2<-cancer[which(cancer$CLEANNAME %in% celinesf),]
  cancer2$muts<-0
  cancer2$muts[which(cancer2$CLEANNAME %in% wildtype)]<-0
  cancer2$muts[which(cancer2$CLEANNAME %in% mutant)]<-1
  
  if(length(table(cancer2$muts))==2){
   
    manyp2<-rep(0,n)
    nm<-mutant
    for(j in 1:n){
      cancer2$mutationgroup<-0
      am<-sample(celinesf,length(nm))
      aw<-celinesf[-which(celinesf %in% am)]
      cancer2[which(cancer2$CLEANNAME %in% am),"mutationgroup"]<-1
      cancer2[which(cancer2$CLEANNAME %in% aw),"mutationgroup"]<-0
      as<-summary(glmer(rank ~mutationgroup+(1|CLEANNAME),data=cancer2,family="poisson"))
      manyp2[j]<- pnorm(as.matrix(as$coefficients)[2,3])
      print(j)
    }
    vvm[i]<-length(which(manyp2<0.05))/n
  }
  print(i)
}

