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
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
cancerd1 = as.data.frame(fread("../Processed_Data/4.finaldata.csv"))
siggenes<-as.data.frame(fread("../OutPut/finalsiggenes.csv"))
siggenesP<-as.data.frame(fread("../OutPut/PancancerResults.csv"))
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")

#################Adjust p-values
head(siggenesP)
n<-length(unique(siggenesP$PvalueNbModel))
p<-siggenesP$PvalueNbModel
siggenesP$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
#########
cancerd1$muts<-0
finalcancer2<-list()
genes<-siggenesP$X[which(siggenesP$adjPvalue<0.05)]

for(i in 1:length(genes)){
cancerf<-cancerd1[which(cancerd1$genename %in% genes[i]),]
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes[i]),]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes[i])]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is mutated
  cancerf$muts[which(cancerf$CLEANNAME %in% wildtype )]<-0
  cancerf$muts[which(cancerf$CLEANNAME %in% mutant )]<-1
finalcancer2[[i]]<-cancerf
}


##########append files of genes
cancerdf<-do.call(rbind, finalcancer2)

  write.csv(cancerdf,"../Processed_Data/PancancerSigFile.csv")
  ####saving significant genes in pancancer with afjusted p-values
  siggenesAdjP<-siggenesP[which(siggenesP$adjPvalue<0.05),]
  write.csv(  siggenesAdjP,"../OutPut/PancancerAdjPvalueSigGenes.csv")
  
  #########################Waterfall Plot
  cancer<-cancerdf[which(cancerdf$genename %in% "PIK3CA"),]
  wildtype<-cancer$CLEANNAME[which(cancer$muts==0)]
  mutant<-cancer$CLEANNAME[which(cancer$muts==1)]
  
  
  cancer2<-aggregate(cancer$rank,by=list(cancer$CLEANNAME),mean)
  cancer2<-data.frame(cancer2)
  cancer2$mutation<-0
  colnames(cancer2)<-c("Cell lines","Meanrank","mutation")
  cancer2$mutation[which(cancer2$`Cell lines` %in% wildtype)]<-0
  cancer2$mutation[which(cancer2$`Cell lines` %in% mutant)]<-1
  cancer2$mutationn<-factor(cancer2$mutation,levels=c(0,1),labels=c("Wildtype","Mutant"))
  ##aggregate(cancer2$Meanrank,by=list(cancer2$mutationn),median)
  mm<-median(cancer2$Meanrank)
  cancer2$rank<-cancer2$Meanrank-mm
  cancerf<-cancer2[order(cancer2$rank),]
  col <- ifelse(cancerf$mutationn == "Wildtype", 
                "steelblue", # if dose = 80 mg, then the color will be steel blue
                "red") # if dose != 80 mg (i.e. 150 mg here), then the color will be cadet blue
  
  MCC<- barplot(cancerf$rank, 
                col=col, 
                border=col, 
                space=0.5, 
                cex.lab=1,  
                ##ylim=c(-50,50),
                main = "Waterfall plot for free-batch-effects- ranks of PIK3CA in Pancancer", 
                ylab="free-batch-effect ranks",
                cex.axis=1,
                
                legend.text= c( "Wildtype", "Mutant"),
                args.legend=list(title="Ranked logFC", fill=c("steelblue", "red"), border=NA, cex=0.75,
                                x= "top"))
  
  #
  
  ###########Density of P-values
  library(ggplot2)
  ggplot(siggenesAdjP, aes(x=PvalueNbModel)) + geom_histogram(binwidth=0.0000095)+
    labs(title="Histogram plot of P-values in Pancancer")+
    theme_classic()
  
  #######################Box Plots of first Quantile of Ranks
  finalcancersq<-list()
  for(i in 1:length(genes)){
    cancer<-cancerdf[which(cancerdf$genename %in% genes[i]),]
    ss<-summary(cancer$rank)
    finalcancersq[[i]]<-cancer[which(cancer$rank<=ss[2]),] 
    
  }
  
  cancersqq<-do.call(rbind,finalcancersq)
  cancer6topg<-cancersqq[which(cancersqq$genename %in% c("KRAS","NRAS","BRAF","PIK3CA",
                                                         "CTNNB1","DHX9","TOP2A","CDK4")),]
  cancer6topg$mutation<-factor(cancer6topg$muts,labels=c("wildtype","mutant")) 
  ggplot(cancer6topg,aes(x=mutation,y=rank,fill=mutation))+geom_boxplot(alpha=0.3)+
    facet_wrap(~genename)+ggtitle("Boxplot of Quantile1-based Ranks of Top Genes in Pancancers")
  
  ##########
  cancerf2<-cancerf[which(cancerf$genename %in% c("KRAS","NRAS","BRAF","PIK3CA",
                                                         "CTNNB1","DHX9","TOP2A","CDK4")),]
 ## cancerf2$mutation<-factor(cancerf2$muts,labels=c("wildtype","mutant"))
  cancermeanbox<-list()
   for(i in 1:length(genes)){
   
  cancer<-cancerdf[which(cancerdf$genename %in% genes[i]),]
  wildtype<-cancer$CLEANNAME[which(cancer$muts==0)]
  mutant<-cancer$CLEANNAME[which(cancer$muts==1)]
  cancer2<-aggregate(cancer$rank,by=list(cancer$CLEANNAME),mean)
  cancer2<-data.frame(cancer2)
  cancer2$mutation<-0
  colnames(cancer2)<-c("Cell lines","Meanrank","mutation")
  cancer2$mutation[which(cancer2$`Cell lines` %in% wildtype)]<-0
  cancer2$mutation[which(cancer2$`Cell lines` %in% mutant)]<-1
  cancer2$mutationn<-factor(cancer2$mutation,levels=c(0,1),labels=c("Wildtype","Mutant"))
  ##aggregate(cancer2$Meanrank,by=list(cancer2$mutationn),median)
  mm<-median(cancer2$Meanrank)
  cancer2$rank<-cancer2$Meanrank-mm
  cancer2$genename<-genes[i]
  cancermeanbox[[i]]<-cancer2
  }
  
  cancermeanboxf<-do.call(rbind,cancermeanbox)
  head(cancermeanboxf)
  cancermeanboxf2<-cancermeanboxf[which(cancermeanboxf$genename %in% c("KRAS","NRAS",
                             "BRAF","PIK3CA","TP53","CTNNB1","DHX9","TOP2A","PA2G4")),]
  ggplot(cancermeanboxf2,aes(x=mutationn,y=rank,fill=mutationn))+geom_boxplot(alpha=0.3)+
    facet_wrap(~genename)+ggtitle("Boxplot of Mean-based Ranks of Top Genes in Pancancers")