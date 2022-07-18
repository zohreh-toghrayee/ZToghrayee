rm(list=ls())
library(dplyr)
##install.packages("stringr", repos='http://cran.us.r-project.org')
library("stringr")
library("dplyr")
library("tidyr")
library("lme4")
library(lmerTest)#
library(sva)

#data(bladderdata)
library(pamr)
library(limma)
library(data.table)
library(reshape2)
##step 1: call two preprocessed data(shRNA nad Drive)
data2=readRDS("../Input_Data/DRIVE/DRIVECountData.RDS")
cancer_data<-data.frame(data2)
shRNAdata=read.csv("../Input_Data/DRIVE/ShrnaGeneMap.RDS")
shRNAdata=read.csv("../Input_Data/DRIVE/shRNAdataf.csv")

########################################

head(cancer_data)
colnames(cancer_data)
mygenes<-unique(shRNAdata$GENESYMBOLS)
myshran<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% c("KRAS","VMP1","SLC47A1","HRH3")),]
shr<-unique(myshran$SEQ)
mydrive<-cancer_data[which(cancer_data$SEQ %in% shr),]
myfile<-left_join(mydrive,myshran,by="SEQ")
head(myfile)

myfile$pool3<-factor(myfile$POOL)
myfile$pool2<-as.numeric(myfile$pool3)
myfile1<-myfile[which(myfile$GENESYMBOLS %in% "KRAS"),]
batch<-aggregate(myfile1$pool2,by=list(myfile1$SHRNA_ID),mean)
batchc<-table(batch$x)
length(which(table(batchc)==1))
mygenes<-unique(myfile$GENESYMBOLS)
onepool<-read.csv("../Processed_Data/2.datasetonepoolf.csv")
twopool<-read.csv("../Processed_Data/2.datasettwopoolf.csv")

genes1<-unique(onepool$genename)
genes2<-unique(twopool$genename)
mygenes1<-c(genes1,genes2)

mygenes<-mygenes1
numpool<-rep(0,length(mygenes))
for(i in 1:length(mygenes)){
  myshran<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% mygenes1[i]),]
  shr<-unique(myshran$SEQ)
  mydrive<-cancer_data[which(cancer_data$SEQ %in% shr),]
  myfile<-left_join(mydrive,myshran,by="SEQ")
  myfile$pool3<-factor(myfile$POOL)
  myfile$pool2<-as.numeric(myfile$pool3)
  batch<-aggregate(myfile$pool2,by=list(myfile$SHRNA_ID),mean)
  batchc<-table(batch$x)
  #length(batchc)
  
    if((length(batchc)>2 & length(which(batchc==1))<=1 )| (length(batchc)==2 & length(which(batchc==1))==0)){
      numpool[i]<-mygenes[i]
    }else{
      numpool[i]<-0
      
    }
  print(i)
  }
  
write.csv(numpool,"../Processed_Data/NumberOfPools.csv")

