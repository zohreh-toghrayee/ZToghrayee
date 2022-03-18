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

#######################################################################3
cancerd = as.data.frame(fread("../Processed_Data/afewgenesdataf.csv"))
head(cancerd)
cancerd<-cancerd[,-c(1,2)]
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
###################
library(glmmTMB)
head(cancerd)
cancer<-cancerd[which(cancerd$genename=="TP53"),]
  #mnb <-summary(glmer.nb(rank ~muts+(1|CLEANNAME),data=cancer,control=glmerControl(optimizer="nloptwrap",tolPwrss=1e-3,optCtrl = list(maxfun = 100000))))
  ##mnb <-summary(glmmTMB(rank ~muts+(1|CLEANNAME),data=cancer,family=Gamma))##error
  mnb1 <-summary(glmmTMB(rank ~muts+(1|CLEANNAME),data=cancer,family=nbinom2))
  cv<-mnb1$coefficients
   pnorm(as.matrix(cv$cond)[2,3])
  
  
  mnb2 <-summary(glmmTMB(rank ~muts+(1|CLEANNAME),data=cancer,family=nbinom1))
  cv<-mnb2$coefficients
  pnorm(as.matrix(cv$cond)[2,3])
  