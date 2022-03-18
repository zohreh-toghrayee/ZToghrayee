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
crispr<-as.data.frame(fread("../Input_Data/CRISPR/CRISPR_gene_effect.csv"))
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")

dim(crispr)
crispr[1:5,1:5]
rownames(crispr)
tmp = strsplit(colnames(crispr),"[()]")
colnames(crispr) = unlist(lapply(tmp, function(x){(x[2])}))
colnames(crispr)[1]<-"DepMap_ID"
which(colnames(crispr)=="10235")
colnames(crispr)[10000:10020]
sampleinfo =as.data.frame(fread("../Input_Data/Mutations/sample_info.csv"))

intersect(unique(sampleinfo$DepMap_ID),unique(crispr$DepMap_ID))
cellinescrispr<-sampleinfo[,1:4]###depmapid and cell lines names
crisprfinal<-merge(crispr,cellinescrispr,by="DepMap_ID")
dim(crisprfinal)
colnames(crisprfinal)[17388:17390]
#################3
crisprfinal[1:5,17388:17390]
#################
breast<-as.data.frame(fread("../Processed_Data/14.breastgenedata.csv"))
dim(breast)
unique(breast$CLEANNAME)

############
genes1<-c("PAX5","RASGRP2","USP19","ROS1","HIPK3","NPFFR2","KRAS","BAZ2B","LYST","TINF2","NAV2","SLX4")
###ENTREZ ID of genes
genes<-c("5079","10235","10869","6098","10114","10886","3845","29994","1130","26277","89797","84464")
tya<-rep(0,length(genes))
tyt<-rep(0,length(genes))
names(tya)<-genes1
names(tya)<-genes1

for(i in 1:length(genes)){
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes1[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes1[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
  cellinesf<-c(wildtype,mutant)
  llln<-length(which(tolower(unique(crisprfinal$stripped_cell_line_name)) %in% cellinesf))
  criss=matrix(0,nrow=llln,ncol=2)
 ccle<-crisprfinal$stripped_cell_line_name[which(tolower(unique(crisprfinal$stripped_cell_line_name)) %in% cellinesf)]
  colnames(criss)<-c("score","group")
  criss[,2]<-0
  rownames(criss)<-tolower(ccle)
  criss[,1]= crisprfinal[which(tolower(unique(crisprfinal$stripped_cell_line_name)) %in% cellinesf),which(colnames(crisprfinal)==genes[i])]
  dim(criss)
  criss[which(rownames(criss) %in% wildtype),2]<-0
  criss[which(rownames(criss) %in% mutant),2]<-1
  
  df<-summary(lm(criss[,1]~criss[,2]))
  df1<-lm(criss[,1]~criss[,2])
  andf<-anova(df1)
  dff<-andf$Df[2]
  # wa[k,j]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
  tya[i]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
  if(length(which(criss[,2]==1))>1){
    tt<-t.test(criss[which(criss[,2]==0),1],criss[which(criss[,2]==1),1])
    tt<-t.test(criss[,1]~criss[,2])
    
    tyt[i]<-pt(tt$statistic,tt$parameter, lower.tail = TRUE)
  }
  
}

write.csv(tya,"../OutPut/CRISPRoverlappResults.csv")
