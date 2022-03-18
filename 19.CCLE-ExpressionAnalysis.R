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
oncosupp<-as.data.frame(read.csv("../processed_Data/19.oncosuppcancerdata.csv"))
head(oncosupp)
mygenes<-unique(oncosupp$genename)
ccled<-as.data.frame(read.csv("../Input_Data/CCLE-Expresssion/CCLE_expression.csv"))
##head(ccled)
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
oncogenes<-as.data.frame(read.csv("../Input_Data/Oncogenes.csv"))
suppresors<-as.data.frame(read.csv("../Input_Data/suppresor.csv"))
head(oncogenes)
oncog<-mygenes[which(mygenes %in% oncogenes$Gene)]
supp<-mygenes[which(mygenes %in% suppresors$Gene)]
presults<-rep(0,length(oncog))
for(i in 1:length(oncog)){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% oncog[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-oncosupp[which(oncosupp$genename %in% oncog[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),median))
  head(myfileag)
  
  myfileagf<-data.frame(celline=myfileag$Group.1  ,med=myfileag$x)
  colnames(myfileagf)<-c("stripped_cell_line_name","Med")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
  
  anadata$expscore<-as.numeric(anadata$expscore)
  anadata$med<-as.numeric(anadata$Med)
    CORTEST<-cor.test(anadata$expscore,anadata$med)
    
    presults[i]<-CORTEST$estimate
  
  
  print(i)
}
names(presults)<-oncog
presultsS<-rep(0,length(supp))
for(i in 1:length(supp)){
  if(length(which(colnames(ccledfinal) %in% supp[i]))>0){
  expdata<-as.data.frame(cbind(ccledfinal[,which(colnames(ccledfinal) %in% supp[i])],ccledfinal$cellines))
  colnames(expdata)<-c("expscore","stripped_cell_line_name")
  myfile<-oncosupp[which(oncosupp$genename %in% supp[i]),]
  myfileag<-as.data.frame(aggregate(myfile$rank,by=list(myfile$CLEANNAME),median))
  head(myfileag)
  myfileagf<-data.frame(celline=myfileag$Group.1  ,med=myfileag$x)
  colnames(myfileagf)<-c("stripped_cell_line_name","Med")
  anadata<-merge(expdata,myfileagf,by="stripped_cell_line_name")
  
  anadata$expscore<-as.numeric(anadata$expscore)
  anadata$med<-as.numeric(anadata$Med)
    CORTEST<-cor.test(anadata$expscore,anadata$med)
    
    presultsS[i]<-CORTEST$estimate
  }
  
  print(i)
}

names(presultsS)<-supp
  length(which(supp<0))
names(presults)<-oncog
length(which(oncog<0))
write.csv(presults,"../OutPut/oncogCorrRankExp.csv")
write.csv(presultsS,"../OutPut/suppCorrRankExp.csv")

panr<-as.data.frame(read.csv("../OutPut/PancancerResults.csv"))
head(panr)
n<-length(unique(panr$PvalueNbModel))
p<-panr$PvalueNbModel
panr$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
mygenes<-panr$X
#####relationship between oncogenes & TS and only significant genes
g1<-panr$X[which(panr$adjPvalue< 0.05)]
sigg<-panr[which(panr$X %in% g1),c("X","adjPvalue")]
presultsonco<-data.frame(names(presults),presults)
colnames(presultsonco)<-c("gene","corr")
head(presultsonco)
colnames(sigg)<-c("gene","adjpval")
oncofile<-left_join(sigg,presultsonco,by="gene")
oncofilef<-oncofile[complete.cases(oncofile$corr),]

presultssupp<-data.frame(names(presultsS),presultsS)
colnames(presultssupp)<-c("gene","corr")
head(presultssupp)
colnames(sigg)<-c("gene","adjpval")
suppfile<-left_join(sigg,presultssupp,by="gene")
suppfilef<-suppfile[complete.cases(suppfile$corr),]

#####relationship between oncogenes & TS and ALL genes
g1<-panr$X
sigg<-panr[which(panr$X %in% g1),c("X","adjPvalue")]
presultsonco<-data.frame(names(presults),presults)
colnames(presultsonco)<-c("gene","correlation")
head(presultsonco)
colnames(sigg)<-c("gene","adjpval")
oncofile<-left_join(sigg,presultsonco,by="gene")
oncofilef<-oncofile[complete.cases(oncofile$correlation),]
oncofilef[which(oncofilef$adjpval<0.05),"sig"]<-1
oncofilef[which(oncofilef$adjpval>0.05),"sig"]<-0
oncofilef$sigg<-factor(oncofilef$sig,labels=c("NonSignificant","Significant"))
oncofilef$genes<-factor(oncofilef$gene)
library(ggplot2)
ggplot() + 
  geom_point(data = oncofilef, aes(x = genes, y = correlation, color = sigg)
             ,size=4) +
  scale_color_manual(values = c("NonSignificant" = "black", "Significant" = "red"))+
  theme(axis.text=element_text(size=6))+
  ggtitle("Correlation between mean rank and expression of Oncogenes in Cell lines in Pancancer")  


presultssupp<-data.frame(names(presultsS),presultsS)
g1<-panr$X
sigg<-panr[which(panr$X %in% g1),c("X","adjPvalue")]
presultssupp<-data.frame(names(presultsS),presultsS)
colnames(presultssupp)<-c("gene","correlation")
head(presultssupp)
colnames(sigg)<-c("gene","adjpval")
suppfile<-left_join(sigg,presultssupp,by="gene")
suppfilef<-suppfile[complete.cases(suppfile$correlation),]
suppfilef[which(suppfilef$adjpval<0.05),"sig"]<-1
suppfilef[which(suppfilef$adjpval>0.05),"sig"]<-0
suppfilef$sigg<-factor(suppfilef$sig,labels=c("NonSignificant","Significant"))
suppfilef$genes<-factor(suppfilef$gene)
library(ggplot2)
ggplot() + 
  geom_point(data = suppfilef, aes(x = genes, y = correlation, color = sigg)
             ,size=4) +
  scale_color_manual(values = c("NonSignificant" = "black", "Significant" = "red"))+
  theme(axis.text=element_text(size=6))+
  ggtitle("Correlation between mean rank and expression of Tummor Suppressors in Cell lines in Pancancer")  

cor.test(oncofilef$adjpval,oncofilef$correlation)
cor.test(suppfilef$adjpval,suppfilef$correlation)
plot(density(oncofilef$correlation),x=oncofilef$gene)

######This is for paper
colnames(oncofile)
oncofile$group<-1
suppfile$group<-0
oncosupp<-rbind(oncofile,suppfile)
oncosuppf<-oncosupp[complete.cases(oncosupp$correlation),]
dim(oncosuppf)
cor.test(-log(oncosupp$adjpval),oncosupp$correlation)
oncosuppf$groupg<-factor(oncosuppf$group,labels=c("Tummor Suppressor","Oncogenes"))


p<-ggplot(oncosuppf, aes(x=correlation, color=groupg)) +
  geom_density()+ggtitle("Density plot of Correlation between ranks and expression of Oncogenes and Tumor suppresors in cell lines in Pancancer ")
p+theme(title=element_text(size=6))+
  geom_label(
    label="P-value=0.002013", 
    x=-0.25,
    y=3.5,
    label.padding = unit(0.35, "lines"), # Rectangle size around label
    label.size = 0.25,
    color = "black",
    fill="#69b3a2"
  )
wilcox.test(oncosuppf$correlation~oncosuppf$group)




