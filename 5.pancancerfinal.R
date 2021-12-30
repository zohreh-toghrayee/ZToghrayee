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
###Required packages for "foreach" package
library(doParallel)
library(data.table)
myCluster <-makeCluster(24, type = "PSOCK")
registerDoParallel(myCluster)
#######################END OF MUTATIONAL MATRIX
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
onep = as.data.frame(fread("../Processed_Data/2.datasetonepoolf.csv"))
twop = as.data.frame(fread("../Processed_Data/2.datasettwopoolf.csv"))

onepool<-onep[,c("CLEANNAME","genename","logFC")]
colnames(onepool)<-c("CLEANNAME","genename","NewlogFC")
twppool<-twop[,c("CLEANNAME","genename","NewlogFC")]
finaldrive<-rbind(onepool,twppool)
finaldrive<-finaldrive[-which(finaldrive$NewlogFC==1000),]


denoisedlogfc<-data.frame(finaldrive)
cancerd1=denoisedlogfc %>%  group_by(CLEANNAME) %>%  mutate(rank = min_rank(NewlogFC))
cancerd<-data.frame(cancerd1)

###############################MUTATIONAL MATRIX

wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
shRNAdata=as.data.frame(fread("../Processed_Data/shRNAdata.csv"))


#library(TMB)
library(glmmTMB)
############# data2 here is drive data after preprocesing

mygenes2<-unique(cancerd$genename)
mygenes3<-mygenes2
twp= foreach(gene=mygenes3, .combine="cbind", .packages=c("glmmTMB","lme4")) %dopar%  {
  twp2<-data.frame()  
  cancer<-cancerd[which(cancerd$genename %in% gene),]
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==gene),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==gene)]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  twp2[1,1]<-length(wildtype)
  twp2[2,1]<-length(mutant)
  cellines<-unique(cancer$CLEANNAME)
  if(length(wildtype) >1 & length(mutant)>1){
    ####
    cancer$muts<-0
    cancer$muts[which(cancer$CLEANNAME %in% wildtype)]<-0
    cancer$muts[which(cancer$CLEANNAME %in% mutant)]<-1
    
    ###mnb <-summary(glmer.nb(rank ~muts+(1|CLEANNAME),data=cancer,control=glmerControl(optimizer="nloptwrap",tolPwrss=1e-3,optCtrl = list(maxfun = 100000))))##not converged
    mnb <-summary(glmmTMB(rank ~muts+(1|CLEANNAME),data=cancer,family="nbinom2"))
    cv<-mnb$coefficients
    twp2[3,1]<- pnorm(as.matrix(cv$cond)[2,3])
    
  }
  
  twp2
}


colnames(twp)<-mygenes3
rownames(twp)<-c("#WildType","#Mutant","P-valueNbModel")
write.csv(twp,"../Output/PancancerResults.csv")

stopCluster(myCluster)