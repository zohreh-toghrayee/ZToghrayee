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
panr<-as.data.frame(read.csv("../OutPut/PancancerResults.csv"))
head(panr)
n<-length(unique(panr$PvalueNbModel))
p<-panr$PvalueNbModel
panr$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
mygenes<-panr$X
g1<-panr$X[which(panr$adjPvalue< 0.05)]
g2<-sample(panr$X[which(panr$adjPvalue> 0.05)],50)
expr<-as.data.frame(read.csv("../OutPut/15.expresults.csv"))
head(expr)
expr$group<-0
expr$group[which(expr$X %in% g1)]<-1
summary(expr$V2[which(expr$group==0)])
summary(expr$V2[which(expr$group==1)])
summary(lm(expr$V2~expr$group))
write.csv(expr,"../OutPut/15.expresultsf.csv")

cancerd = as.data.frame(fread("../Processed_Data/FewGenesData.csv"))
head(cancerd)
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
###################
mygenes<-unique(cancerd$genename)
cancerd$muts<-0
cancerg<-list()
for(i in 1:length(mygenes)){
  cancer<-cancerd[which(cancerd$genename %in% mygenes[i]),]
wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes[i])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated

cancer$muts[which(cancer$CLEANNAME %in% wildtype & cancer$genename==mygenes[i])]<-0
cancer$muts[which(cancer$CLEANNAME %in% mutant & cancer$genename==mygenes[i])]<-1
cancerg[[i]]<-cancer
}
dataset1<-do.call(rbind, cancerg)
write.csv(dataset1,"../Processed_Data/afewgenesdataf.csv")
head(dataset1)
dim(dataset1)
dim(cancerd)
table(dataset1$muts)
dataset1$mutsF<-factor(dataset1$muts,levels=c(0,1),labels=c("Wildtype","mutant"))
head(dataset1)
table(dataset1$mutsF)
dk<-dataset1[which(dataset1$genename %in% "KRAS"),]
table(dk$muts)
########
ccled<-as.data.frame(read.csv("../Input_Data/CCLE-Expresssion/CCLE_expression.csv"))
head(ccled)
ccled[1:5,1:5]
colnames(ccled)###cell lines

dim(ccled)###cols are genes and rows are cell lines
sampleinfo = read.csv("E:\\Git\\genesilencing_project\\org_data\\CCLE\\sample_info.csv")
colnames(ccled)[colnames(ccled) == "X"] <- "DepMap_ID"
#intersect(unique(sampleinfo$DepMap_ID),unique(ccled$DepMap_ID))
cellinesccled<-sampleinfo[,1:4]
ccledfinal<-merge(ccled,cellinesccled,by="DepMap_ID")
ccledfinal$cellines<-tolower(ccledfinal$stripped_cell_line_name)
##ccledfinal$stripped_cell_line_name
dim(ccledfinal)
unique(colnames(ccledfinal))
tmp = strsplit(colnames(ccledfinal), "[..]")
colnames(ccledfinal) = unlist(lapply(tmp, function(x){(x[1])}))
presults<-matrix(0,nrow=length(mygenes),ncol=2)


mygenes2<-colnames(ccledfinal[,which(colnames(ccledfinal) %in% mygenes)])
presults<-data.frame()
i=which(mygenes2 %in% "TOP2A")
for(i in 1:length(mygenes2)){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% mygenes2[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-dataset1[which(dataset1$genename %in% mygenes2[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),summary))
  head(myfileag)
  ww<-myfile[which(myfile$muts==0),"CLEANNAME"]
  mm<-myfile[which(myfile$muts==1),"CLEANNAME"]
  myfileag$Group.1
  myfileagf<-data.frame(celline=myfileag$Group.1  ,q1=myfileag$x[,3])
  colnames(myfileagf)<-c("stripped_cell_line_name","Q1")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
  
  anadata$expscore<-as.numeric(anadata$expscore)
  anadata$Q1<-as.numeric(anadata$Q1)
  anadata$muts[which(anadata$stripped_cell_line_name %in% ww)]<-0
  anadata$muts[which(anadata$stripped_cell_line_name %in% mm)]<-1
  summary(anadata$expscore[which(anadata$muts==0)])
  summary(anadata$expscore[which(anadata$muts==1)])
  table(anadata$muts)
  if(sd(anadata$expscore)>0 & sd(anadata$meanrank)>0){
    CORTEST<-cor.test(anadata$expscore,anadata$Q1)
    presults[i,1]<-CORTEST$p.value
    presults[i,2]<-CORTEST$estimate
  }
  
  print(i)
}
rownames(presults)<-mygenes2

#########3
oncogenes<-as.data.frame(read.csv("../Input_Data/Oncogenes.csv"))
suppresors<-as.data.frame(read.csv("../Input_Data/suppresor.csv"))
head(oncogenes)
oncog<-mygenes[which(mygenes %in% oncogenes$Gene)]
supp<-mygenes[which(mygenes %in% suppresors$Gene)]
presults<-rep(0,length(oncog))
for(i in 1:length(oncog)){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% oncog[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-dataset1[which(dataset1$genename %in% oncog[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),median))
  head(myfileag)
  ww<-myfile[which(myfile$muts==0),"CLEANNAME"]
  mm<-myfile[which(myfile$muts==1),"CLEANNAME"]
  myfileag$Group.1
  myfileagf<-data.frame(celline=myfileag$Group.1  ,med=myfileag$x)
  colnames(myfileagf)<-c("stripped_cell_line_name","Med")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
  
  anadata$expscore<-as.numeric(anadata$expscore)
  anadata$med<-as.numeric(anadata$Med)
  anadata$muts[which(anadata$stripped_cell_line_name %in% ww)]<-0
  anadata$muts[which(anadata$stripped_cell_line_name %in% mm)]<-1
  #summary(anadata$expscore[which(anadata$muts==0)])
  #summary(anadata$expscore[which(anadata$muts==1)])
  #table(anadata$muts)
    CORTEST<-cor.test(anadata$expscore,anadata$med)
    
    presults[i]<-CORTEST$estimate
  
  
  print(i)
}
presultsS<-rep(0,length(supp))
for(i in 1:length(supp)){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% supp[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-dataset1[which(dataset1$genename %in% supp[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),median))
  head(myfileag)
  ww<-myfile[which(myfile$muts==0),"CLEANNAME"]
  mm<-myfile[which(myfile$muts==1),"CLEANNAME"]
  myfileag$Group.1
  myfileagf<-data.frame(celline=myfileag$Group.1  ,med=myfileag$x)
  colnames(myfileagf)<-c("stripped_cell_line_name","Med")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
  
  anadata$expscore<-as.numeric(anadata$expscore)
  anadata$med<-as.numeric(anadata$Med)
  anadata$muts[which(anadata$stripped_cell_line_name %in% ww)]<-0
  anadata$muts[which(anadata$stripped_cell_line_name %in% mm)]<-1
  #summary(anadata$expscore[which(anadata$muts==0)])
  #summary(anadata$expscore[which(anadata$muts==1)])
  #table(anadata$muts)
  if(sd(anadata$expscore)>0 & sd(anadata$meanrank)>0){
    CORTEST<-cor.test(anadata$expscore,anadata$med)
    if(!is.na(CORTEST$estimate)){
    presultsS[i]<-CORTEST$estimate
  }
  
  print(i)
}

names(presultsS)<-supp
  
names(presults)<-oncog
