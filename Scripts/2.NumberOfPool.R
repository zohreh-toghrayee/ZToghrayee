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
data2=readRDS("../Input_Data/PreprocessingData/DRIVEdata.RDS")
cancer_data<-data.frame(data2)
shRNAdata=read.csv("../Input_Data/PreprocessingData/shRNAdata.csv")

##########################################################determine genes
fewgenes<-unique(shRNAdata$GENESYMBOLS)
myfinaltwopoolgenes<-rep(0,length(fewgenes))

for(k in 1:length(fewgenes)){
  shrna<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% fewgenes[k]),"SEQ" ]
  cancer_data2<-cancer_data[which(cancer_data$SEQ %in% shrna),]
  shRNAdata2<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% fewgenes[k]), ]
  cancer_dataf= left_join(cancer_data2, shRNAdata2, by="SEQ")
    ##############making cancer data
   batch <- aggregate(cancer$spscat,by=list(cancer$SHRNA_ID),mean)
  batchc<-rep((batch$x),1)###sps
  if(length(which(table(batchc)==1))==0){
    myfinalonepoolgenes[k]<-fewgenes1[k]
  }else{
 myfinalonepoolgenes[k]<-0
}
print("k1")
print(k)
}
names(myfinaltwopoolgenes)<-fewgenes1
#####################################################original file for get denoised logFC
write.csv(myfinaltwopoolgenes,"../Input_Data/Processed_Data/NumberOfPools.csv")
