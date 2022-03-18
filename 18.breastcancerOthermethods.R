rm(list=ls())
library(dplyr)
library("stringr")
library("dplyr")
library("tidyr")
library("lme4")
library(lmerTest)#
library("dplyr")
library("tidyr")
library("lme4")
library(lmerTest)
library(doParallel)
library(data.table)

##########################bind variable including cancer name to denoised DRIVE data
annotation = read.csv("TableS2.csv", stringsAsFactors = F)
cancers2<-strsplit(annotation$PATHOLOGIST_ANNOTATION,"[:]")
annotation[,"cancer"]<-unlist(lapply(cancers2,function(x) {tolower(x[1])}))
cancers<-unique(annotation$cancer)
colnames(annotation)<-c("CLEANNAME","PRIMARY_SITE","PATHOLOGY_ANNOTATION","inCCLE","SOURCE","PRIDE","cancer")
cancerd1 = as.data.frame(fread("4.finaldata.csv"))
cancerd<-merge(cancerd1,annotation,by="CLEANNAME")
###############################MUTATIONAL MATRIX
wild_matix<-readRDS("wild_matix.R")
cancer_mutations3<-readRDS("cancer_mutations.R")

cancers<-c("leukemia","eye","soft_tissue","gastric","cns","kidney","thyroid","skin"                     
  ,"salivary_gland","ovary","lung","bone","endometrium","bladder","upper_aerodigestive_tract","pnet","breast","pancreas","colorectal",
  "prostate","undefined","liver","lymphoma","oesophagus","cervix","biliary_tract")  


breast<-cancerd[which(cancerd$cancer=="breast" ),]
##############Since colnames of CRISPR data are ENTREZID of egenes:
entrez<-as.data.frame(fread("entrezz.csv"))
entrez2 = strsplit(entrez[,1], "[\033]")
entrez[,1] = unlist(lapply(entrez2, function(x){x[1]}))
entrez[,2] = unlist(lapply(entrez2, function(x){x[2]}))
colnames(entrez)<-c("gene","EntrezId")
###
crispr<-as.data.frame(fread("CRISPR_gene_effect.csv"))
tmp = strsplit(colnames(crispr),"[()]")
colnames(crispr) = unlist(lapply(tmp, function(x){(x[2])}))
colnames(crispr)[1]<-"DepMap_ID"
sampleinfo =as.data.frame(fread("sample_info.csv"))

cellinescrispr<-sampleinfo[,1:4]###depmapid and cell lines names
crisprfinal<-merge(crispr,cellinescrispr,by="DepMap_ID")
#################
breastg<-as.data.frame(fread("eachcancerResults.csvbreast.csv"))
genes1<-breastg$V1

###ENTREZ ID of genes: replace entrezid in crispr file by genes 

crisprg<-crisprfinal[,which(colnames(crisprfinal) %in% entrez$EntrezId)]
for(i in 1:length(colnames(crisprg))){
  colnames(crisprg)[i]<-entrez$gene[which(entrez$EntrezId %in% colnames(crisprg)[i] )]
  
}
stripped_cell_line_name<-crisprfinal$stripped_cell_line_name
crisprfinal2<-cbind(crisprg,stripped_cell_line_name)
genes<-colnames(crisprg)
tyac<-rep(0,length(genes))
names(tyac)<-genes
for(i in 1:length(genes)){
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
  cellinesf<-c(wildtype,mutant)
  llln<-length(which(tolower(unique(crisprfinal2$stripped_cell_line_name)) %in% cellinesf))
  criss=matrix(0,nrow=llln,ncol=2)
 ccle<-crisprfinal2$stripped_cell_line_name[which(tolower(unique(crisprfinal2$stripped_cell_line_name)) %in% cellinesf)]
  colnames(criss)<-c("score","group")
   rownames(criss)<-tolower(ccle)
  criss[,1]= crisprfinal2[which(tolower(unique(crisprfinal2$stripped_cell_line_name)) %in% cellinesf),which(colnames(crisprfinal2)==genes[i])]
   criss[which(rownames(criss) %in% wildtype),2]<-0
  criss[which(rownames(criss) %in% mutant),2]<-1
  
  df<-summary(lm(criss[,1]~criss[,2]))
  df1<-lm(criss[,1]~criss[,2])
  andf<-anova(df1)
  dff<-andf$Df[2]
  tyac[i]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
  
  }
write.csv(tyac,"BreastCRISPR.csv")

##########
ATARIS<-readRDS("DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))

DEMETER2<-read.csv("Demeter2_DRIVE_gene_dep_scores.csv")###demeter results
tmpd = strsplit(colnames(DEMETER2), "[_]")
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
DEMETER2[,1] = word(DEMETER2[,1], 1)
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))

RSA<-readRDS("DRIVE_RSA_data.RDS")###RSA results
tmpd = strsplit(colnames(RSA), "[_]")
colnames(RSA) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
###########ATARIS
genes2<-rownames(ATARIS)[which(rownames(ATARIS) %in% genes1)]
tyaa<-rep(0,length(genes2))
names(tyaa)<-genes2
for(i in 1:length(genes2)){
print(i)
print(genes2[i])
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
print(dfr1)
    df<-lm(ataris[,1]~ataris[,2])
    andf1<-anova(df)
    dff1<-andf1$Df[2]
if(!is.na(dfr1$coefficients)){
    tyaa[i]<-pt(as.matrix(dfr1$coefficients)[2,3],dff1)
    }
  } else{
    tyaa[i]<-99
  }

}
write.csv(tyaa,"BreastATARIS.csv")

#############DEMETER
genesd<-DEMETER2[which(DEMETER2[,1] %in% genes1),1]
tyad<-rep(0,length(genesd))
names(tyad)<-genesd
  for(i in 1:length(genesd)){
print(i)
    wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genesd[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
    mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genesd[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
    cellines<-c(wildtype,mutant)
    demeter2=matrix(0,nrow=length(which(colnames(DEMETER2) %in% cellines)),ncol=2)
    
  rownames(demeter2)=colnames(DEMETER2[which(DEMETER2[,1]==genesd[i]),which(colnames(DEMETER2) %in% cellines)])
  demeter2[,2]<-0
  demeter2[,1]= t(DEMETER2[which(DEMETER2[,1]==genesd[i]),which(colnames(DEMETER2) %in% cellines)])
  demeter2[which(rownames(demeter2) %in% mutant),2]<-1
  demeter2[which(rownames(demeter2) %in% wildtype),2]<-0
  if(length(table(demeter2[,2]))>1){
    dfr2<-summary(lm(demeter2[,1]~demeter2[,2]))
    df2<-lm(demeter2[,1]~demeter2[,2])
    andf2<-anova(df2)
    dff2<-andf2$Df[2]
    
    tyad[i]<-pt(as.matrix(dfr2$coefficients)[2,3],dff2)
    
  }else{
    tyad[i]<-99
  }
}

write.csv(tyad,"BreastDEMETER.csv")

######RSA
  genesr<-rownames(RSA)[which(rownames(RSA) %in% genes1)]
  tyar<-rep(0,length(genesr))
  names(tyar)<-genesr
  for(i in 1:length(genesr)){
print(i)
    wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genesr[i]),]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is not mutated
    mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genesr[i])]==1)),unique(breast$CLEANNAME))###cellines where the ith gene is mutated
    cellines<-c(wildtype,mutant)
      rsa=matrix(0,nrow=length(which(colnames(RSA) %in% cellines)),ncol=2)

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
    
    tyar[i]<-pt(as.matrix(dfr3$coefficients)[2,3],dff3)
    
  }else{
    tyar[i]<-99
  }

}


write.csv(tyar,"BreastRSA.csv")


