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
cancerd = as.data.frame(fread("../Processed_Data/4.finaldata.csv"))
print(colnames(cancerd))
###print(unique(cancerd$genename))

###############################MUTATIONAL MATRIX

wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")



################################################3
 # mygenes2<- c("TP53","KRAS","PIK3CA","BRAF","NRAS","HRAS","CTNNB1","CDK4","APC","AXIN1",
##"ZNF423","BAZ2A","MUC4","ABCA13","COPS4","AKNA","MYO15A")
 
  mygenes2<-c("TP53","KRAS","PIK3CA")
  n<-50
library(glmmTMB)
myresults<-data.frame()

for(j in 1:length(mygenes2)){
cancer<-cancerd[which(cancerd$genename %in% mygenes2[j]),]
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes2[j]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes2[j])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  cancer$muts<-0
  cancer$muts[which(cancer$CLEANNAME %in% wildtype)]<-0
  cancer$muts[which(cancer$CLEANNAME %in% mutant)]<-1
   
  celinesf<-unique(cancer$CLEANNAME)
  nm<-mutant
  num<-c(1:50)
twp= foreach(i=num, .combine="rbind", .packages=c("glmmTMB","lme4")) %dopar%  {
  twp2<-data.frame() 
if(length(table(cancer$muts))==2){  
      cancer$mutationgroup<-0
      am<-sample(celinesf,length(nm))
      aw<-celinesf[-which(celinesf %in% am)]
      cancer[which(cancer$CLEANNAME %in% am),"mutationgroup"]<-1
      cancer[which(cancer$CLEANNAME %in% aw),"mutationgroup"]<-0
      mnb <-summary(glmmTMB(rank ~mutationgroup+(1|CLEANNAME),data=cancer,family="nbinom2"))
     cv<-mnb$coefficients
     twp2[1,1]<- pnorm(as.matrix(cv$cond)[2,3])
  
    }else{
twp2[1,1]<-99
}

twp2
}
myresults[,j]<-twp
print(j)
}

print("finish")
write.csv(myresults,"twp.csv",row.names=FALSE)


row.names(myresults)<-c(1:50)
colnames(myresults)<-mygenes2
write.csv(myresults,"TypeOneErrorResults.csv")
stopCluster(myCluster)
