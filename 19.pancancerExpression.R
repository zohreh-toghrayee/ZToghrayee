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
panr<-as.data.frame(read.csv("PancancerResults.csv"))
n<-length(unique(panr$PvalueNbModel))
p<-panr$PvalueNbModel
panr$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
mygenes<-panr$X

cancerd = as.data.frame(fread("4.finaldata.csv"))

########
ccled<-as.data.frame(read.csv("CCLE_expression.csv"))
sampleinfo = read.csv("sample_info.csv")
colnames(ccled)[colnames(ccled) == "X"] <- "DepMap_ID"
cellinesccled<-sampleinfo[,1:4]
ccledfinal<-merge(ccled,cellinesccled,by="DepMap_ID")
ccledfinal$cellines<-tolower(ccledfinal$stripped_cell_line_name)

tmp = strsplit(colnames(ccledfinal), "[..]")
colnames(ccledfinal) = unlist(lapply(tmp, function(x){(x[1])}))
presults<-matrix(0,nrow=length(mygenes),ncol=2)
#########3
oncogenes<-as.data.frame(read.csv("Oncogenes.csv"))
suppresors<-as.data.frame(read.csv("suppresor.csv"))
oncog<-mygenes[which(mygenes %in% oncogenes$Gene)]
supp<-mygenes[which(mygenes %in% suppresors$Gene)]
oncosupp<-cancerd[which(cancerd$genename %in% c(oncog,supp)),]
write.csv(oncosupp,"19.oncosuppcancerdata.csv")
presults<-rep(0,length(oncog))
presultsg<-rep("0",length(oncog))

for(i in 1:length(oncog)){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% oncog[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-cancerd[which(cancerd$genename %in% oncog[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),median))
   myfileagf<-data.frame(celline=myfileag$Group.1  ,med=myfileag$x)
  colnames(myfileagf)<-c("stripped_cell_line_name","Med")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
   anadata$expscore<-as.numeric(anadata$expscore)
  anadata$med<-as.numeric(anadata$Med)
   
    CORTEST<-cor.test(anadata$expscore,anadata$med)
    presults[i]<-CORTEST$estimate
 presultsg[i]<-oncog[i]
  print(i)
}
presultsS<-rep(0,length(supp))
presultsSg<-rep("0",length(supp))

for(i in 1:length(supp)){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% supp[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-cancerd[which(cancerd$genename %in% supp[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),median))
  myfileagf<-data.frame(celline=myfileag$Group.1  ,med=myfileag$x)
  colnames(myfileagf)<-c("stripped_cell_line_name","Med")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
  
  anadata$expscore<-as.numeric(anadata$expscore)
  anadata$med<-as.numeric(anadata$Med)
  
    CORTEST<-cor.test(anadata$expscore,anadata$med)
    
    presultsS[i]<-CORTEST$estimate
   presultsSg[i]<-supp[i]
  print(i)
}


write.csv(presults,"19.pancancerOncoExp.csv")
write.csv(presultsg,"19.pancancerOncoExpg.csv")

write.csv(presultsS,"19.pancancerSuppExp.csv")
write.csv(presultsSg,"19.pancancerSuppExpg.csv")

