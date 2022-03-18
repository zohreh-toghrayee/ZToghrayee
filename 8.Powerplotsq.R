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
KRASk =as.data.frame(fread("../OutPut/8.powermSIM2Sq2.csv"))
dim(KRASk)

q<-2
ncm<-seq(from=2,to=24,by=2)
ncm<-seq(from=2,to=12,by=1)

ncw<-q*ncm


b<-500
dim(KRASk)
head(KRASk)
KRASk[,1]<-rep(ncm,each=b)

####################plotting power of two methods
#############these are the results file
powerplotmymethod<-rep(0,length(ncm))
powerataris<-rep(0,length(ncm))
powerdemeter<-rep(0,length(ncm))
powerrsa<-rep(0,length(ncm))

  for(j in 1:length(ncm)){
  powerplotmymethod[j]<-length(which(KRASk[which(KRASk[,1]==ncm[j]),2]<0.05))
  powerataris[j]<-length(which(KRASk[which(KRASk[,1]==ncm[j]),3]<0.05))
  powerdemeter[j]<-length(which(KRASk[which(KRASk[,1]==ncm[j]),4]<0.05))
  powerrsa[j]<-length(which(KRASk[which(KRASk[,1]==ncm[j]),5]<0.05))

  }

step<-ncm
m<-4####the number of methods to compare

KRASpowerplot<-matrix(0,nrow=length(ncm),ncol=m+1)
KRASpowerplot[,1]<-step
KRASpowerplot[,2]<-powerplotmymethod/b
KRASpowerplot[,3]<-powerataris/b
KRASpowerplot[,4]<-powerdemeter/b
KRASpowerplot[,5]<-powerrsa/b

colnames(KRASpowerplot)<-c("Mutantsize","Mymethod","ATARIS","DEMETER","RSA")
KRASpowerplot2<-data.frame(samplesize=rep(KRASpowerplot[,1],m),power=c(KRASpowerplot[,2],KRASpowerplot[,3],
                           KRASpowerplot[,4],KRASpowerplot[,5]),method=rep(1:m,each=length(ncm)))
colnames(KRASpowerplot2)<-c("Mutantsize","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS","DEMETER","RSA"))
library(ggplot2)
p<-ggplot(KRASpowerplot2, aes(x=Mutantsize, y=power, group=Method)) +
  geom_line(aes(color=Method))

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS,DEMETER,RSA in Pancancer based on simulation(q=1)")

#############
KRASpowerplot2f<-KRASpowerplot2[order(KRASpowerplot2$Mutantsize),]
##KRASpowerplot2f<-KRASpowerplot2f[-1,]
head(KRASpowerplot2f)
library(ggplot2)
library(stats)
library(tibble)
lo1<-loess(power[which(KRASpowerplot2f$method==1)]~Mutantsize[which(KRASpowerplot2f$method==1)],data=KRASpowerplot2f,span=0.50)
lo2<-loess(power[which(KRASpowerplot2f$method==2)]~Mutantsize[which(KRASpowerplot2f$method==2)],data=KRASpowerplot2f,span=0.50)
lo3<-loess(power[which(KRASpowerplot2f$method==3)]~Mutantsize[which(KRASpowerplot2f$method==3)],data=KRASpowerplot2f,span=0.50)
lo4<-loess(power[which(KRASpowerplot2f$method==4)]~Mutantsize[which(KRASpowerplot2f$method==4)],data=KRASpowerplot2f,span=0.50)

smoothed101 <- predict(lo1)
smoothed102 <- predict(lo2)
smoothed103 <- predict(lo3)
smoothed104 <- predict(lo4)

#
plot(KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==1)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==1)], type="l", main="Power of different methods in q=4,500 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==2)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==2)], type="l", main="Power of different methods in q=4,500 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==3)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==3)], type="l", main="Power of different methods in q=4,500 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==4)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==4)], type="l", main="Power of different methods in q=4,500 simulations", xlab="Sample-size", ylab="Power")

lines(smoothed101, x=KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==1)], col="red")
lines(smoothed102, x=KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==2)], col="blue")
lines(smoothed103, x=KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==3)], col="green")
lines(smoothed104, x=KRASpowerplot2f$Mutantsize[which(KRASpowerplot2f$method==4)], col="black")
legend(30, 0.5, legend=c("Mymethod", "ATARIS","DEMETER","RSA"),
       col=c("red", "blue","green","black"), lty=1:2, cex=1)

seq(from=2,to=24,by=2)

#########################
mm<-seq(from=2,to=24,by=2)
length(mm)
KRASk =as.data.frame(fread("../OutPut/8.powermSIM2Sq51.csv"))
dim(KRASk)
ncm<-seq(from=1,to=12,by=1)
ncm<-seq(from=2,to=24,by=2)
q<-5
step<-ncm+q*ncm
KRASk[,1]<-rep(step,each=b)

head(KRASk)

b<-500
mydata1<-data.frame()
for(j in 1:length(step)){
  mydata1[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata1[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata1[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata1[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata1[j,5]<-step[j]
  
  }


mydata2<-data.frame()
for(j in 1:length(step)){
  mydata2[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata2[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata2[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata2[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata2[j,5]<-step[j]
  
}

mydata3<-data.frame()
for(j in 1:length(step)){
  mydata3[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata3[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata3[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata3[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata3[j,5]<-step[j]
  
}
mydata4<-data.frame()
for(j in 1:length(step)){
  mydata4[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata4[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata4[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata4[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata4[j,5]<-step[j]
  
}
mydata5<-data.frame()
for(j in 1:length(step)){
  mydata5[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata5[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata5[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata5[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata5[j,5]<-step[j]
  
}
mydata6<-data.frame()
for(j in 1:length(step)){
  mydata6[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata6[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata6[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata6[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata6[j,5]<-step[j]
  
}
mydata61<-data.frame()
for(j in 1:length(step)){
  mydata61[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata61[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata61[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata61[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata61[j,5]<-step[j]
  
}
mydata51<-data.frame()
for(j in 1:length(step)){
  mydata51[j,1]<-length(which(KRASk[which(KRASk[,1]==step[j]),2]<0.05))/b
  mydata51[j,2]<-length(which(KRASk[which(KRASk[,1]==step[j]),3]<0.05))/b
  mydata51[j,3]<-length(which(KRASk[which(KRASk[,1]==step[j]),4]<0.05))/b
  mydata51[j,4]<-length(which(KRASk[which(KRASk[,1]==step[j]),5]<0.05))/b
  mydata51[j,5]<-step[j]
  
}

mydata<-list(mydata1,mydata2,mydata3,mydata4,mydata5,mydata6,mydata61,mydata51)
final<-do.call("rbind",mydata)
dim(final)
table(final$V5)
write.csv(final,"../OutPut/mydatapower.csv")
mym<-as.data.frame(aggregate(final$V1,by=list(final$V5),mean))
am<-as.data.frame(aggregate(final$V2,by=list(final$V5),mean))
dm<-as.data.frame(aggregate(final$V3,by=list(final$V5),mean))
rm<-as.data.frame(aggregate(final$V4,by=list(final$V5),mean))

final2<-list(mym,am,dm,rm)
mydataf<-do.call("cbind",final2)
KRASpowerplot<-mydataf[,c(1,2,4,6,8)]
colnames(KRASpowerplot)<-c("samplesize","Mymethod","ATARIS","DEMETER","RSA")
ncm2<-unique(KRASpowerplot$samplesize)
KRASpowerplot2<-data.frame(samplesize=rep(KRASpowerplot[,1],m),power=c(KRASpowerplot[,2],KRASpowerplot[,3],
                                                                       KRASpowerplot[,4],KRASpowerplot[,5]),method=rep(1:m,each=length(ncm2)))
colnames(KRASpowerplot2)<-c("samplesize","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS","DEMETER","RSA"))
library(ggplot2)
KRASpowerplot2fg<-KRASpowerplot2[which(KRASpowerplot2$samplesize<50),]
p<-ggplot(KRASpowerplot2fg, aes(x=samplesize, y=power, group=Method)) +
  geom_line(aes(color=Method))

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS,DEMETER,RSA in Pancancer based on total simulation")

library(ggplot2)
library(stats)
library(tibble)
KRASpowerplot2f<-KRASpowerplot2fg[order(KRASpowerplot2fg$samplesize),]

lo1<-loess(power[which(KRASpowerplot2f$method==1)]~samplesize[which(KRASpowerplot2f$method==1)],data=KRASpowerplot2f,span=0.50)
lo2<-loess(power[which(KRASpowerplot2f$method==2)]~samplesize[which(KRASpowerplot2f$method==2)],data=KRASpowerplot2f,span=0.50)
lo3<-loess(power[which(KRASpowerplot2f$method==3)]~samplesize[which(KRASpowerplot2f$method==3)],data=KRASpowerplot2f,span=0.50)
lo4<-loess(power[which(KRASpowerplot2f$method==4)]~samplesize[which(KRASpowerplot2f$method==4)],data=KRASpowerplot2f,span=0.50)

smoothed101 <- predict(lo1)
smoothed102 <- predict(lo2)
smoothed103 <- predict(lo3)
smoothed104 <- predict(lo4)

#
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==1)], type="l", main="Power of different methods in 500 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==2)], type="l", main="Power of different methods in 500 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==3)], type="l", main="Power of different methods in 500 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==4)], type="l", main="Power of different methods in 500 simulations", xlab="Sample-size", ylab="Power")

lines(smoothed101, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)], col="red")
lines(smoothed102, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)], col="blue")
lines(smoothed103, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)], col="green")
lines(smoothed104, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)], col="black")
legend(30, 0.5, legend=c("Mymethod", "ATARIS","DEMETER","RSA"),
       col=c("red", "blue","green","black"), lty=1:2, cex=1)

#sam<-read.csv("../OutPut/SIMTSample2.csv")
#ncm<-as.numeric(sam$nr_tot)
#prop<-sam$nr_mut/sam$nr_wt
dim(powerm)
head(powerm)
mean(powerm$V1)
powerm<-powerm[,-1]
powerplotmymethod<-rep(0,length(ncm))

for(i in 1:length(ncm)){
  powerplotmymethod[i]<-100*(length(which(powerm[,i]<0.05))/100)
  
  
  
}

step<-ncm+ncw
m<-1####the number of methods to compare
KRASpowerplot<-matrix(0,nrow=length(ncm),ncol=m+1)
KRASpowerplot[,1]<-step
KRASpowerplot[,2]<-powerplotmymethod


colnames(KRASpowerplot)<-c("samplesize","Mymethod")
KRASpowerplot2<-data.frame(samplesize=rep(KRASpowerplot[,1],m),power=c(KRASpowerplot[,2]),method=rep(1:m,each=length(ncm)))
colnames(KRASpowerplot2)<-c("samplesize","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS","DEMETER","RSA"))
library(ggplot2)
p<-ggplot(KRASpowerplot2, aes(x=samplesize, y=power)) +
  geom_line()

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS,DEMETER,RSA in Pancancer based on simulation(ES=0.5)")



#######################line bar plot
head(final)
library(ggplot2)

colnames(final)<-c("mymethod","ATARiS","DEMETER","RSA","Samplesize")
final2<-final[which(final$Samplesize<50),]
sd<-sd(final2$mymethod)
p<- ggplot(final2, aes(x=Samplesize, y=mymethod)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mymethod-sd, ymax=mymethod+sd), width=.2,
                position=position_dodge(0.08))
print(p)

#################################################for different Effect Sizes
zarib<-seq(0.1,1.5,by=0.1)
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

