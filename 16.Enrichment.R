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
#######################END OF MUTATIONAL MATRIX
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
breast<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvbreast.csv"))
n<-length(unique(breast$PvalueNbModel))
p<-breast$PvalueNbModel
breast$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
head(breast)
colnames(breast)<-c("gene","NWildType","NMutant","PvalueNbModel","adjPvalue")
breast[which(breast$adjPvalue<0.05),]
