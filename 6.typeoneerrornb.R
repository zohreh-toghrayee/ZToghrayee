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
onep = read.csv("2.datasetonepoolf.csv")
twop = read.csv("2.datasettwopoolf.csv")
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

wild_matix<-readRDS("wild_matix.R")
cancer_mutations3<-readRDS("cancer_mutations.R")

################################################

################################################3
  mygenes2<- c("TP53","KRAS","PIK3CA","BRAF","NRAS","HRAS","CTNNB1","CDK4","APC","AXIN1")
  vvm<-rep(0,length(mygenes2))
  names(vvm)<-mygenes2
  n<-200
library(glmmTMB)
i=which(mygenes2=="KRAS")
##for(i in 1:length(mygenes2)){
  cancer<-cancerd[which(cancerd$genename %in% mygenes2[i]),]
  celinesf<-unique(cancer$CLEANNAME)
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes2[i])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  #others<-intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[i]),]==2)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  #celinesf<-cancelf[-which(cancelf %in% others)]
  #cancer2<-cancer[which(cancer$CLEANNAME %in% c(wildtype,mutant)),]
  cancer$muts<-0
  cancer$muts[which(cancer$CLEANNAME %in% wildtype)]<-0
  cancer$muts[which(cancer$CLEANNAME %in% mutant)]<-1
  
  if(length(table(cancer$muts))==2){
   
    manyp2<-rep(0,n)
    nm<-mutant
    for(j in 1:n){
      cancer$mutationgroup<-0
      am<-sample(celinesf,length(nm))
      aw<-celinesf[-which(celinesf %in% am)]
      cancer[which(cancer$CLEANNAME %in% am),"mutationgroup"]<-1
      cancer[which(cancer$CLEANNAME %in% aw),"mutationgroup"]<-0
      mnb <-summary(glmer.nb(rank ~mutationgroup+(1|CLEANNAME),data=cancer, verbose=TRUE))
     manyp2[j]<- pnorm(as.matrix(mnb$coefficients)[2,3])
   print(j)
    }
    vvm[i]<-length(which(manyp2<0.05))/n
  }
  #print(i)
#}


write.csv(vvm,"6.typeoneerror.csv")