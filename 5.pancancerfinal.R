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

#######################END OF MUTATIONAL MATRIX
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
onep = read.csv("../Processed_Data/2.datasetonepoolf.csv")
twop = read.csv("../Processed_Data/2.datasettwopoolf.csv")
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

wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
shRNAdata=read.csv("../Processed_Data/shRNAdata.csv")

#####
###############well-known methods 
ATARIS<-readRDS("../Input_Data/DRIVE/DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))

DEMETER2<-read.csv("../Input_Data/DRIVE/Demeter2_DRIVE_gene_dep_scores.csv")###demeter results
tmpd = strsplit(colnames(DEMETER2), "[_]")
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
DEMETER2[,1] = word(DEMETER2[,1], 1)
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))

RSA<-readRDS("../Input_Data/DRIVE/DRIVE_RSA_data.RDS")###RSA results
tmpd = strsplit(colnames(RSA), "[_]")
colnames(RSA) = unlist(lapply(tmpd, function(x){tolower(x[1])}))


#library(TMB)
library(glmmTMB)
############# data2 here is drive data after preprocesing

mygenes2<-unique(cancerd$genename)
twp2<-matrix(0,nrow=length(mygenes2),ncol=6)
rownames(twp2)<-mygenes2
colnames(twp2)<-c("N.wildtype","N.Mutant","ATARiS","DEMETER","RSA","Poisson")
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

twp2[i,1]<-length(wildtype)
twp2[i,2]<-length(mutant)
cellines<-unique(cancer$CLEANNAME)
if(length(wildtype) >1 & length(mutant)>1){
####################
######ataris/demeter/rsa
 ataris=matrix(0,nrow=length(colnames(ATARIS) ),ncol=2)
  if(length(which(rownames(ATARIS)==mygenes2[i])!=0)){
  rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS)==mygenes2[i]),which(colnames(ATARIS) %in% cancelf)])
  ataris[,1]= t(ATARIS[which(rownames(ATARIS)==mygenes2[i]),which(colnames(ATARIS) %in% cancelf)])
  ataris[which(rownames(ataris) %in% mutant),2]<-1
  ataris[which(rownames(ataris) %in% wildtype),2]<-0
if(length(table(ataris[,2]))>1){
dfr1<-summary(lm(ataris[,1]~ataris[,2]))
df<-lm(ataris[,1]~ataris[,2])
andf1<-anova(df)
    dff1<-andf1$Df[2]
  twp2[i,3]<-pt(as.matrix(dfr1$coefficients)[2,3],dff1)

} 
  }

demeter2=matrix(0,nrow=length(colnames(DEMETER2)),ncol=2)
    if(length(which(DEMETER2[,1]==mygenes2[i]))!=0){
      
    rownames(demeter2)=colnames(DEMETER2[which(DEMETER2[,1]==mygenes2[i]),which(colnames(DEMETER2) %in% cancelf)])
    demeter2[,2]<-0
    demeter2[,1]= t(DEMETER2[which(DEMETER2[,1]==mygenes2[i]),which(colnames(DEMETER2) %in% cancelf)])
    demeter2[which(rownames(demeter2) %in% mutant),2]<-1
    demeter2[which(rownames(demeter2) %in% wildtype),2]<-0
if(length(table(demeter2[,2]))>1){
dfr2<-summary(lm(demeter2[,1]~demeter2[,2]))
df2<-lm(demeter2[,1]~demeter2[,2])
andf2<-anova(df2)
dff2<-andf2$Df[2]

  twp2[i,4]<-pt(as.matrix(dfr2$coefficients)[2,3],dff2)
   
    }
}

rsa=matrix(0,nrow=length(colnames(RSA)),ncol=2)
    if(length(which(rownames(RSA)==mygenes2[i]))!=0){
    rownames(rsa)=colnames(RSA[which(rownames(RSA)==mygenes2[i]),which(colnames(RSA) %in% cancelf)])
    rsa[,2]<-0
    rsa[,1]= t(RSA[which(rownames(RSA)==mygenes2[i]),which(colnames(RSA) %in% cancelf)])
    rsa[which(rownames(rsa) %in% mutant),2]<-1
    rsa[which(rownames(rsa) %in% wildtype),2]<-0 
 if(length(table(rsa[,2]))>1){
dfr3<-summary(lm(rsa[,1]~ rsa[,2]))   
df3<-lm(rsa[,1]~ rsa[,2])
andf3<-anova(df3)
dff3<-andf3$Df[2]

  twp2[i,5]<-pt(as.matrix(dfr3$coefficients)[2,3],dff3)

    }
    }

####
  
  
   as<-summary(glmer(cancer$rank ~cancer$muts+(1|cancer$CLEANNAME),family=poisson))
    twp2[i,6]<-pnorm(as.matrix(as$coefficients)[2,3])
  }

print(i)
}



twp2<-twp2[which(twp2[,1]>0),]
twp2<-twp2[order(twp2[,1]),]

write.csv(twp2,"../Input_Data/Output/PancancerResults.csv")
