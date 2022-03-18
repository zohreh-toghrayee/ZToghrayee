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

######################################################################
library(corrplot)
library(RColorBrewer)
corm<-matrix(0,nrow=11,ncol=6)
rownames(corm)<-c("PAX5","RASGRP2","USP19","ROS1","HIPK3","NPFFR2",
                  "KRAS","BAZ2B","LYST","TINF2","NAV2")
colnames(corm)<-c("OurMethod","ATARIS","DEMETER","RSA","CRISPR","Scientific Literature")
corm[,1]<-1
corm[,2:5]<-0.5
corm[,6]<-1
M <-corm
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))