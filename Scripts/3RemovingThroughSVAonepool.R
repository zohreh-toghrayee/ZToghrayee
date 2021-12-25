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
##step 1: Call data including shRNA /DRIVE/NumberOfPool
data2= readRDS("../Input_Data/Processed_Data/DRIVEdata.RDS")
cancer_data<-data.frame(data2)
shRNAdata=read.csv("../Input_Data/Processed_Data/shRNAdata.csv")

myfinalgenes<-read.csv("../Input_Data/Processed_Data/NumberOfPools.csv")
##########################################################determine genes
mygenes<-data.frame(myfinalgenes)
##head(mygenes)
colnames(mygenes)<-c("id","genes")
geneonepool<-mygenes$id[which(mygenes$genes=="0")]
fewgenes<-intersect(unique(shRNAdata$GENESYMBOLS),geneonepool)

##########################################for one pool

vv<-matrix(1000,nrow=1,ncol=6)
colnames(vv)<-c("CLEANNAME","sampleshrna","SAMPLE_COUNT","SHRNA_ID","mutationgroup","logFC")
finalcancer2<-list()

for(k in 1:length(fewgenes)){
   shrna<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% fewgenes[k]),"SEQ" ]
  cancer_dataf<-cancer_data[which(cancer_data$SEQ %in% shrna),]
  shRNAGeneMap2<-shRNAdata[which(shRNAdata$GENESYMBOLS %in% fewgenes[k]), ]
  cancer= left_join(cancer_dataf, shRNAGeneMap2, by="SEQ")
  
 cancermm1<-cancer[,c("SHRNA_ID","CLEANNAME","PLASMID_COUNT")]
cancermm2<-cancer[,c("SHRNA_ID","CLEANNAME","SAMPLE_COUNT")]

cancer$POOL2<-factor(cancer$POOL)
cancer$pool3<-as.numeric(cancer$POOL2)

###########################adjust for plasmid

mybatchwide1<-reshape(data=cancermm1,
                        idvar="CLEANNAME",
                        v.name="PLASMID_COUNT",
                        timevar="SHRNA_ID",
                        direction="wide")
 mybatchwide1<-mybatchwide1[which(complete.cases(mybatchwide1)), ]
 ##mybatchwide1[is.na(mybatchwide1)]<-0
  mybatchwide12<-as.matrix(mybatchwide1[,-1])
  mybatchwide12<-mybatchwide12+0.5
  ss<-as.vector(summary(cancer$SPS))
  cancer$spscat<-0
  cancer$spscat[which(cancer$SPS <=ss[2])]<-0
  cancer$spscat[which(cancer$SPS > ss[2] & cancer$SPS <=ss[5])]<-1
  cancer$spscat[which(cancer$SPS >ss[5])]<-2
  
  ####for sample count
  mybatchwide<-reshape(data=cancermm2,
                       idvar="CLEANNAME",
                       v.name="SAMPLE_COUNT",
                       timevar="SHRNA_ID",
                       direction="wide")
 mybatchwide<-mybatchwide[which(complete.cases(mybatchwide)), ]
 ##mybatchwide[is.na(mybatchwide)]<-0
  mybatchwide2<-as.matrix(mybatchwide[,-1])
  mybatchwide2<-mybatchwide2+0.5
  batch <- aggregate(cancer$spscat,by=list(cancer$SHRNA_ID),mean)
  #batchrc<-matrix(rep((batchr$x),1),nrow=dim(mybatchwide2)[2],ncol=1)
  batchc<-rep((batch$x),1)
if(length(which(table(batchc)==1))==0 & length(table(batchc))>1){
  modcombat = model.matrix(~ 1,data=cancermm2)
  adjusted_counts2 <- ComBat_seq(mybatchwide2, batch=batchc)
  adjusted_counts<-(adjusted_counts2/(mybatchwide12))
  newdata<-data.frame(cbind(adjusted_counts,mybatchwide[,1]))
  

newdata2<-melt(data=newdata,
        id.vars=tail(names(newdata),1),
        variable.name = "SAMPLE_COUNT",
        value.name = "COUNT"
        )


newdata2$SAMPLE_COUNT<-as.character(newdata2$SAMPLE_COUNT)
tmp2 = strsplit(newdata2$SAMPLE_COUNT, "[.]")
newdata2$SAMPLE_COUNT2 = unlist(lapply(tmp2, function(x){tolower(x[2])}))
colnames(newdata2)<-c("CLEANNAME","sampleshrna","SAMPLE_COUNT","SHRNA_ID")
newdata2$mutationgroup<-0
newdata2$logFC<-log(as.numeric(newdata2$SAMPLE_COUNT))

finalcancer2[[k]]<-newdata2
}else{
    
    finalcancer2[[k]]<-vv
  }
print(k)

}


#####################################add name of genes to the each file
finalcancers3<-list()
for (k in  1:length(fewgenes)){
  
  # Create the first data if no data exist yet
  finalcancers1<-data.frame(finalcancer2[[k]])
  finalcancers1$genename<-fewgenes[k]
  finalcancers3[[k]]<-finalcancers1[,c("CLEANNAME","SAMPLE_COUNT","genename","logFC")]
 
  
  print(k)
}##########append files of genes
dataset2<-do.call(rbind, finalcancers3)
######################################################FOR RANKING
write.csv(dataset2,"../Input_Data/Processed_Data/2.datasetonepoolf.csv")


