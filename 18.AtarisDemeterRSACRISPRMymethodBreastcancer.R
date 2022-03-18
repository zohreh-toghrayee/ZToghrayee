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
breastataris<-as.data.frame(read.csv("../Processed_Data/BreastATARIS.csv"))
head(breastataris)
dim(breastataris)

breastdemeter<-as.data.frame(read.csv("../Processed_Data/BreastDEMETER.csv"))
head(breastdemeter)
dim(breastdemeter)

breastrsa<-as.data.frame(read.csv("../Processed_Data/BreastRSA.csv"))
head(breastrsa)
dim(breastrsa)
breastcrispr<-as.data.frame(read.csv("../Processed_Data/BreastCRISPR.csv"))
head(breastcrispr)
dim(breastcrispr)
breastmymethod<-as.data.frame(read.csv("../OutPut/EachCancerResults/eachcancerResults.csvbreast.csv"))
head(breastmymethod)
dim(breastmymethod)
  n<-length(unique(breastataris[,2]))
  p<-breastataris[,2]
  breastataris[,3]<-p.adjust(p, method = "fdr", n = length(p))
  
  n<-length(unique(breastdemeter[,2]))
  p<-breastdemeter[,2]
  breastdemeter[,3]<-p.adjust(p, method = "fdr", n = length(p))
  n<-length(unique(breastrsa[,2]))
  p<-breastrsa[,2]
  breastrsa[,3]<-p.adjust(p, method = "fdr", n = length(p))
  n<-length(unique(breastcrispr[,2]))
  p<-breastcrispr[,2]
  breastcrispr[,3]<-p.adjust(p, method = "fdr", n = length(p))
  n<-length(unique(breastmymethod$PvalueNbModel))
  p<-breastmymethod$PvalueNbModel
  breastmymethod$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
  
  colnames(breastataris)<-c("gene","Pvalue","adjPvalue")
  colnames(breastdemeter)<-c("gene","Pvalue","adjPvalue")
  colnames(breastrsa)<-c("gene","Pvalue","adjPvalue")
  colnames(breastcrispr)<-c("gene","Pvalue","adjPvalue")
  colnames(breastmymethod)<-c("gene","Nwildtype","Nmutant","Pvalue","adjPvalue")
  ####join my method
file1<-left_join(breastmymethod,breastataris,by="gene")
file2<-left_join(file1,breastdemeter,by="gene")
file3<-left_join(file2,breastrsa,by="gene")
allmethodsbreast<-left_join(file3,breastcrispr,by="gene")
head(allmethodsbreast)
write.csv(allmethodsbreast,"../OutPut/allmethods.csv")

panr<-as.data.frame(read.csv("../OutPut/PancancerResults.csv"))
head(panr)
n<-length(unique(panr$PvalueNbModel))
p<-panr$PvalueNbModel
panr$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
panrf<-panr[which(panr$adjPvalue<0.05),]
write.csv(panrf,"../OutPut/PancancerResultsSignificantGenes.csv")
