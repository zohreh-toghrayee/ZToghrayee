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
library(glmmTMB)
#######################END OF MUTATIONAL MATRIX
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
KRASk =as.data.frame(fread("../Processed_Data/KRASdata.csv"))
cancer<-data.frame(KRASk)
cellinesf<-unique(cancer$CLEANNAME)
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")

##############significance of KRAS in real data
wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)=="KRAS"),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)=="KRAS")]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
cancer$muts<-0
cancer$muts[which(cancer$CLEANNAME %in% wildtype )]<-0
cancer$muts[which(cancer$CLEANNAME %in% mutant )]<-1
mnb <-summary(glmmTMB(rank ~muts+(1|CLEANNAME),data=cancer,family="nbinom2"))
cv<-mnb$coefficients
pnorm(as.matrix(cv$cond)[2,3])
###########################prepare for power analysis
####files will be used in power section

cellinesf<-unique(cancer$CLEANNAME)


ATARIS<-readRDS("../Input_Data/DRIVE/DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))
ataris=matrix(0,nrow=length(which(colnames(ATARIS) %in% cellinesf)),ncol=4)

colnames(ataris)<-c("score","group","score2","samescale")
ataris[,3]<-0
rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS)=="KRAS"),which(colnames(ATARIS) %in% cellinesf)])
ataris[,1]= t(ATARIS[which(rownames(ATARIS)=="KRAS"),which(colnames(ATARIS) %in% cellinesf)])
dim(ataris)

###makin same scale atarius to the ranks 
g<-dim(ataris)[1]
srank<-summary(cancer$rank)
atarisor<-ataris[order(ataris[,1]),]
dim(atarisor)
atarisor[1,3]<-srank[1]
atarisor[g,3]<-srank[6]
dif<-(atarisor[g,1]-atarisor[1,1])
dh<-dim(atarisor)[1]
rot<-(srank[6]-srank[1])
for(i in 2:dh){
  atarisor[i,3]<- (rot*(atarisor[i,1]-atarisor[1,1])/dif)+2
}

#####for a specified effect size and different sample size
b<-100###number of simulations
ncm<-c(4,6,16,25,35,55)
ncw<-c(5,12,18,30,40,60)
resf<-matrix(0,nrow=b,ncol=length(ncm))
wa<-matrix(0,nrow=b,ncol=length(ncm)) 
raw<-matrix(0,nrow=b,ncol=length(ncm)) #resf<-rep(0,b)
d<-0.5
effects<-d*sqrt(var(cancer$rank))
for(j in 1:length(ncm)){
  for(k in 1:b){
###cancer
    cancer$mutshufle<-0
    cancer$rankv<-0
    WW<-sample(cellinesf,ncw[j])
    MM<-sample(cellinesf[-which(cellinesf %in% WW)],ncm[j])
    cancer[which(cancer$CLEANNAME %in% WW),"mutshufle"]<-0
    cancer[which(cancer$CLEANNAME %in% MM),"mutshufle"]<-1
    cancer[which(cancer$CLEANNAME %in% WW),"rankv"]<-cancer$rank[which(cancer$CLEANNAME %in% WW)]+round(effects)
    cancer[which(cancer$CLEANNAME %in% MM),"rankv"]<-cancer$rank[which(cancer$CLEANNAME %in% MM)]
        
     mnb <-summary(glmmTMB(rankv ~mutshufle+(1|CLEANNAME),data=cancer,family="nbinom2"))
     cv<-mnb$coefficients
     resf[k,j]<- pnorm(as.matrix(cv$cond)[2,3])

##ataris
    ataris2<-atarisor[which(rownames(atarisor) %in% cellinesf),]
    ataris2[which(rownames(ataris2) %in% WW),2]<-0
    ataris2[which(rownames(ataris2) %in% MM),2]<-1
    ataris2[which(rownames(ataris2) %in% WW),4]<-ataris2[which(rownames(ataris2) %in% WW),3]+round(effects)
    ataris2[which(rownames(ataris2) %in% MM),4]<-ataris2[which(rownames(ataris2) %in% MM),3]
    df<-summary(lm(ataris2[,4]~ataris2[,2]))
    df1<-lm(ataris2[,4]~ataris2[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]
    wa[k,j]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
###raw
    print(k)
    
  }
  print(j)
}

####################plotting power of two methods
#############these are the results file
wa<-as.data.frame(read.csv("../Processed_Data/ATARISpower.csv"))
resf<-as.data.frame(read.csv("../Processed_Data/MyMethodpower.csv"))
head(wa)
head(resf)
wa<-wa[,-1]
resf<-resf[,-1]
powerplotmymethod<-rep(0,length(ncm))
powerataris<-rep(0,length(ncm))

for(i in 1:length(ncm)){
  powerplotmymethod[i]<-length(which(resf[,i]<0.05))
  powerataris[i]<-length(which(wa[,i]<0.05))
  
  
}
step<-ncw+ncm
m<-2####the number of methods to compare
KRASpowerplot<-matrix(0,nrow=length(ncm),ncol=m+1)
KRASpowerplot[,1]<-step
KRASpowerplot[,2]<-powerplotmymethod/100
KRASpowerplot[,3]<-powerataris/100

colnames(KRASpowerplot)<-c("samplesize","Mymethod","ATARIS")
KRASpowerplot2<-data.frame(rep(KRASpowerplot[,1],2),c(KRASpowerplot[,2],KRASpowerplot[,3]),rep(1:2,each=6))
colnames(KRASpowerplot2)<-c("samplesize","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS"))
library(ggplot2)
p<-ggplot(KRASpowerplot2, aes(x=samplesize, y=power, group=Method)) +
  geom_line(aes(color=Method))

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS, and raw data in Pancancer based on simulation(ES=0.5)")

#############
#################################################for different Effect Sizes
zarib<-seq(0.001,1,by=0.05)
resf<-matrix(0,nrow=b,ncol=length(zarib))
wa<-matrix(0,nrow=b,ncol=length(zarib))
raw<-matrix(0,nrow=b,ncol=length(zarib)) 


for(j in 1:length(zarib)){
####cancer data
    effects<-zarib[j]*sqrt(var(cancer$rank))
    effects2<-zarib[j]*sqrt(var(krasd$meanlogfc))
    for(k in 1:b){
    cancer$mutshufle<-0
    cancer$rankv<-0
    WW<-sample(cellinesf,length(wildtype))
    MM<-sample(cellinesf[-which(cellinesf %in% WW)],length(mutant))
    cancer2[which(cancer2$CLEANNAME %in% WW),"mutshufle"]<-0
    cancer2[which(cancer2$CLEANNAME %in% MM),"mutshufle"]<-1
    cancer2[which(cancer2$CLEANNAME %in% WW),"rankv"]<-cancer2$rank[which(cancer2$CLEANNAME %in% WW)]+round(effects)
    cancer2[which(cancer2$CLEANNAME %in% MM),"rankv"]<-cancer2$rank[which(cancer2$CLEANNAME %in% MM)]
    msh2<-summary(glmer(rankv ~mutshufle+(1|CLEANNAME),data=cancer2,family="poisson"))
    resf[k,j]<-pnorm(as.matrix(msh2$coefficients)[2,3])
####ataris data
    ataris2<-atarisor[which(rownames(atarisor) %in% cellinesf),]
    ataris2[which(rownames(ataris2) %in% WW),2]<-0
    ataris2[which(rownames(ataris2) %in% MM),2]<-1
    ataris2[which(rownames(ataris2) %in% WW),4]<-ataris2[which(rownames(ataris2) %in% WW),3]+round(effects)
    ataris2[which(rownames(ataris2) %in% MM),4]<-ataris2[which(rownames(ataris2) %in% MM),3]
    df<-summary(lm(ataris2[,4]~ataris2[,2]))
    df1<-lm(ataris2[,4]~ataris2[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]
    wa[k,j]<-pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
#####raw data
    krasd2<-krasd[which(krasd$celline %in% cellinesf),]
    krasd2[which(krasd2$celline %in% WW),"muts"]<-0
    krasd2[which(krasd2$celline %in% MM),"muts"]<-1
    krasd2[which(krasd2$celline %in% WW),"news"]<-krasd2[which(krasd2$celline %in% WW),"meanlogfc"]+effects2
    krasd2[which(krasd2$celline %in% MM),"news"]<-krasd2[which(krasd2$celline %in% MM),"meanlogfc"]
    dfr<-summary(lm(krasd2[,"news"]~krasd2[,"muts"]))
    dfr1<-lm(krasd2[,"news"]~krasd2[,"muts"])
    andf1<-anova(dfr1)
    dff1<-andf1$Df[2]
    raw[k,j]<-pt(as.matrix(dfr$coefficients)[2,3],dff1)
    
    print(k)
    
  }
  print(j)
}



#####plotting power plot for different three methods
powerplotmymethod<-rep(0,length(zarib))
powerataris<-rep(0,length(zarib))
powerraw<-rep(0,length(zarib))

for(i in 1:length(zarib)){
  powerplotmymethod[i]<-length(which(resf[,i]<0.05))
  powerataris[i]<-length(which(wa[,i]<0.05))
  powerraw[i]<-length(which(raw[,i]<0.05))
  
}
step<-ncw+ncm
m<-3####the number of methods to compare
KRASpowerplot<-matrix(0,nrow=length(zarib),ncol=4)
KRASpowerplot[,1]<-zarib
KRASpowerplot[,2]<-powerplotmymethod/100
KRASpowerplot[,3]<-powerataris/100
KRASpowerplot[,4]<-powerraw/100

colnames(KRASpowerplot)<-c("ES","Mymethod","ATARIS","powerraw")
KRASpowerplot2<-data.frame(rep(KRASpowerplot[,1],3),c(KRASpowerplot[,2],KRASpowerplot[,3],KRASpowerplot[,4]),rep(1:3,each=20))
colnames(KRASpowerplot2)<-c("ES","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS","raw"))
write.csv(KRASpowerplot2,"F:/THESIS-CODE/Finalresults/finalonesidetest/powerKRASallmethodsfordifferentzaribs.csv")
library(ggplot2)
p<-ggplot(KRASpowerplot2, aes(x=ES, y=power, group=Method)) +
  geom_line(aes(color=Method))

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS, and raw data in Pancancer based on simulations for different ES")
pdf("../Output/my_plot.pdf")

# Creating a plot
plot(p+ggtitle("Power plot of our method,ATARIS, and raw data in Pancancer based on simulations for different ES"))

# Closing the graphical device
dev.off() 
#################