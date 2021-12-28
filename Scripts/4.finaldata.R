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


write.csv(cancerd,"4.finaldata.csv")
