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

write.csv(cancerd,"../Input_Data/Processed_Data/KRASrankeddata.csv")
#################raw data
data2= readRDS("../Input_Data/Processed_Data/DRIVEdata.RDS")
cancer_data<-data.frame(data2)
shRNAdata=read.csv("../Input_Data/Processed_Data/shRNAdata.csv")

  shrna<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% "KRAS"),"SEQ" ]
  cancer_dataf<-cancerd[which(cancerd$SEQ %in% shrna),]
  shRNAGeneMap2<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% "KRAS"), ]
  cancer= left_join( cancer_dataf, shRNAGeneMap2, by="SEQ")
cancer$logFC=cancer$SAMPLE_COUNT/cancer$PLASMID_COUNT
write.csv(cancer,"../Input_Data/Processed_Data/KRASrawdata.csv")

