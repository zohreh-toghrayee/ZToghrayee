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
entrez<-as.data.frame(fread("../Input_Data/entrezz.csv"))
head(entrez)
class(entrez)
entrez[,1]
entrez2 = strsplit(entrez[,1], "[\033]")
entrez[,1] = unlist(lapply(entrez2, function(x){x[1]}))
entrez[,2] = unlist(lapply(entrez2, function(x){x[2]}))
colnames(entrez)<-c("gene","EntrezId")
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
crisprfinal[1:4,1:4]

colnames(crisprfinal)[17388:17390]
#################3
crisprfinal[1:5,17388:17390]
#################
breastg<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvbreast.csv"))
breastg$V1
breast<-as.data.frame(fread("../Processed_Data/18.BreastCancerData.csv"))
dim(breast)
unique(breast$CLEANNAME)
head(breast)
genes1<-unique(breast$genename)
############
genes1<-c("PAX5","RASGRP2","USP19","ROS1","HIPK3","NPFFR2","KRAS","BAZ2B","LYST","TINF2","NAV2","SLX4")
###ENTREZ ID of genes
genes2<-unique(colnames(crisprfinal))
genes2<-c("5079","10235","10869","6098","10114","10886","3845","29994","1130","26277","89797","84464")


rownames(tya)<-genes1
crisprg<-crisprfinal[,which(colnames(crisprfinal) %in% entrez$EntrezId)]
for(i in 1:length(colnames(crisprg))){
  colnames(crisprg)[i]<-entrez$gene[which(entrez$EntrezId %in% colnames(crisprg)[i] )]
  
}
crisprg[1:4,1:4]
dim(crisprg)
stripped_cell_line_name<-crisprfinal$stripped_cell_line_name
crisprfinal2<-cbind(crisprg,stripped_cell_line_name)
genes<-colnames(crisprg)
tyac<-rep(0,length(genes))
for(i in 1:length(genes)){
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes1[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
  cellinesf<-c(wildtype,mutant)
  llln<-length(which(tolower(unique(crisprfinal2$stripped_cell_line_name)) %in% cellinesf))
  criss=matrix(0,nrow=llln,ncol=2)
 ccle<-crisprfinal2$stripped_cell_line_name[which(tolower(unique(crisprfinal2$stripped_cell_line_name)) %in% cellinesf)]
  colnames(criss)<-c("score","group")
  criss[,2]<-0
  rownames(criss)<-tolower(ccle)
  criss[,1]= crisprfinal2[which(tolower(unique(crisprfinal2$stripped_cell_line_name)) %in% cellinesf),which(colnames(crisprfinal2)==genes[i])]
  dim(criss)
  criss[which(rownames(criss) %in% wildtype),2]<-0
  criss[which(rownames(criss) %in% mutant),2]<-1
  
  df<-summary(lm(criss[,1]~criss[,2]))
  df1<-lm(criss[,1]~criss[,2])
  andf<-anova(df1)
  dff<-andf$Df[2]
  # wa[k,j]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
  tya[i,4]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
  
  
}

#########
ATARIS<-readRDS("../Input_Data/DRIVE/DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))

DEMETER2<-read.csv("../Input_Data/DRIVE/Demeter2_DRIVE_gene_dep_scores.csv")###demeter results
tmpd = strsplit(colnames(DEMETER2), "[_]")
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
DEMETER2[,1] = word(DEMETER2[,1], 1)
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))

RSA<-readRDS("../Input_Data/DRIVE/DRIVE_RSA_data.RDS")###RSA results
tmpd = strsplit(colnames(RSA), "[_]")
colnames(RSA) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
genes2<-rownames(ATARIS)[which(rownames(ATARIS) %in% genes1)]
tyaa<-rep(0,length(genes2))
for(i in 1:length(genes2)){
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes2[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes2[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
  cellines<-c(wildtype,mutant)
  ataris=matrix(0,nrow=length(which(colnames(ATARIS) %in% cellines )),ncol=2)
  
  rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS)==genes2[i]),which(colnames(ATARIS) %in% cellines)])
  ataris[,1]= t(ATARIS[which(rownames(ATARIS)==genes2[i]),which(colnames(ATARIS) %in% cellines)])
  ataris[which(rownames(ataris) %in% mutant),2]<-1
  ataris[which(rownames(ataris) %in% wildtype),2]<-0
  if(length(table(ataris[,2]))>1){
    dfr1<-summary(lm(ataris[,1]~ataris[,2]))
    df<-lm(ataris[,1]~ataris[,2])
    andf1<-anova(df)
    dff1<-andf1$Df[2]
    tyaa[i]<-pt(as.matrix(dfr1$coefficients)[2,3],dff1)
    
  } else{
    tyaa[i]<-99
  }

}

genesd<-DEMETER2[which(DEMETER2[,1] %in% genes1),1]
  for(i in 1:length(genesd)){
    wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genesd[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
    mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genesd[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
    cellines<-c(wildtype,mutant)

    demeter2=matrix(0,nrow=length(which(colnames(DEMETER2) %in% cellines)),ncol=2)
    
  rownames(demeter2)=colnames(DEMETER2[wchih(DEMETER2[,1]==genesd[i]),which(colnames(DEMETER2) %in% cellines)])
  demeter2[,2]<-0
  demeter2[,1]= t(DEMETER2[which(DEMETER2[,1]==genesd[i]),which(colnames(DEMETER2) %in% cellines)])
  demeter2[which(rownames(demeter2) %in% mutant),2]<-1
  demeter2[which(rownames(demeter2) %in% wildtype),2]<-0
  if(length(table(demeter2[,2]))>1){
    dfr2<-summary(lm(demeter2[,1]~demeter2[,2]))
    df2<-lm(demeter2[,1]~demeter2[,2])
    andf2<-anova(df2)
    dff2<-andf2$Df[2]
    
    tya[i,2]<-pt(as.matrix(dfr2$coefficients)[2,3],dff2)
    
  }else{
    tya[i,2]<-99
  }
}
  genesr<-rownames(RSA)[which(rownames(RSA) %in% genes1)]
  tyaa<-rep(0,length(genes2))
  
  rsa=matrix(0,nrow=length(which(colnames(RSA) %in% cellines)),ncol=2)
  for(i in 1:length(genesr)){
    wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genesr[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
    mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genesr[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
    cellines<-c(wildtype,mutant)
    
  rownames(rsa)=colnames(RSA[which(rownames(RSA)==genesr[i]),which(colnames(RSA) %in% cellines)])
  rsa[,2]<-0
  rsa[,1]= t(RSA[which(rownames(RSA)==genesr[i]),which(colnames(RSA) %in% cellines)])
  rsa[which(rownames(rsa) %in% mutant),2]<-1
  rsa[which(rownames(rsa) %in% wildtype),2]<-0 
  if(length(table(rsa[,2]))>1){
    dfr3<-summary(lm(rsa[,1]~ rsa[,2]))   
    df3<-lm(rsa[,1]~ rsa[,2])
    andf3<-anova(df3)
    dff3<-andf3$Df[2]
    
    tya[i,3]<-pt(as.matrix(dfr3$coefficients)[2,3],dff3)
    
  }else{
    tya[i,3]<-99
  }

}
i=1

######
names(tyaa)<-genes2
genes2
genes1

write.csv(tya,"../OutPut/breastothermethodResults1.csv")
write.csv(tyaa,"../OutPut/breastothermethodResults2.csv")


