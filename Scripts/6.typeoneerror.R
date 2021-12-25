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
onep = read.csv("../Input_Data/Processed_Data/2.datasetonepoolf.csv")
twop = read.csv("../Input_Data/Processed_Data/2.datasettwopoolf.csv")
onep <-data.frame(onep)
twop<-data.frame(twop)


onepool<-onep[,c("CLEANNAME","genename","logFC")]
colnames(onepool)<-c("CLEANNAME","genename","NewlogFC")
twppool<-twop[,c("CLEANNAME","genename","NewlogFC")]
finaldrive<-rbind(onepool,twppool)
finaldrive<-finaldrive[-which(finaldrive$NewlogFC==1000),]


denoisedlogfc<-data.frame(finaldrive)
cancerd1=denoisedlogfc %>%  group_by(CLEANNAME) %>%  mutate(rank = min_rank(NewlogFC))
cancerd<-data.frame( cancerd1)

###############################MUTATIONAL MATRIX

wild_matix<-readRDS("../Input_Data/Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Input_Data/Processed_Data/cancer_mutations.R")
shRNAdata=read.csv("../Input_Data/Processed_Data/shRNAdata.csv")
################################################

################################################3
mygenes2<-which(wgenes %in% c("TP53","KRAS","PIK3CA","BRAF","NRAS","HRAS","CTNNB1","CDK4","APC","AXIN1"))
vvm<-rep(0,length(mygenes2))
names(vvm)<-wgenes
n<-1000
for(i in 1:length(mygenes2)){
  cancer<-cancerd[which(cancerd$genename %in% mygenes2[i]),]
  cancelf<-unique(cancer$CLEANNAME)
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes2[i])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  others<-intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==2)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  celinesf<-cancelf[-which(cancelf %in% others)]
  cancer2<-cancer[which(cancer$CLEANNAME %in% celinesf),]
  cancer2$muts<-0
  cancer2$muts[which(cancer2$CLEANNAME %in% wildtype)]<-0
  cancer2$muts[which(cancer2$CLEANNAME %in% mutant)]<-1

  if(length(table(cancer$muts))==2){
  cellines<-unique(cancer$CLEANNAME)
  manyp2<-rep(0,n)
  nm<-which(cancer$muts==0)
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

write.csv(vvm,"../Output/typeoneerror.csv")