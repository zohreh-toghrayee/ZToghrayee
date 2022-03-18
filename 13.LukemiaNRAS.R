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




####step2: call merged cancer data and annotation(did in server) and mutation matrices
cancerdf<-as.data.frame(fread("../Processed_Data/LukemiaData.csv"))
head(cancerdf)
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")


mydata<-cancerdf[which(cancerdf$genename %in% "NRAS"),]  
wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)=="NRAS"),]==1)),unique(mydata$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)=="NRAS")]==1)),unique(mydata$CLEANNAME))###cellines where the ith gene is mutated
md<-mydata[which(mydata$CLEANNAME %in% mutant),]
head(md)
