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
powerm =as.data.frame(fread("../OutPut/powermSIM2.csv"))
powera =as.data.frame(fread("../OutPut/poweraSIM2.csv"))
powerd =as.data.frame(fread("../OutPut/powerdSIM2.csv"))
powerr =as.data.frame(fread("../OutPut/powerrSIM2.csv"))

dim(powerm)
dim(powera)
dim(powerd)
dim(powerr)
powerm =as.data.frame(fread("../OutPut/powermSIM2S.csv"))
#############these are the results file
head(powerm)
resf<-powerm[,-1]
simt =as.data.frame(fread("../OutPut/SIMTSample2.csv"))
ncm<-dim(simt)[1]

powerplotmymethod<-rep(0,length(ncm))
powerataris<-rep(0,length(ncm))
powerdemeter<-rep(0,length(ncm))
powerrsa<-rep(0,length(ncm))
for(i in 1:length(ncm)){
  powerplotmymethod[i]<-length(which(resf[,i]<0.05))
  powerataris[i]<-length(which(resf[,i]<0.05))
  powerdemeter[i]<-length(which(resf[,i]<0.05))
  powerrsa[i]<-length(which(resf[,i]<0.05))
  
  
}
head(simt)
simt$nr_tot
dim(simt)
resf<-powerm[,-1]
tya<-powera[,-1]
tyd<-powerd[,-1]
tyr<-powerr[,-1]

ncm<-simt$nr_mut
ncw<-simt$nr_wt




powerplotmymethod<-rep(0,length(ncm))
powerataris<-rep(0,length(ncm))
powerdemeter<-rep(0,length(ncm))
powerrsa<-rep(0,length(ncm))
for(i in 1:length(ncm)){
  powerplotmymethod[i]<-length(which(resf[,i]<0.05))
  powerataris[i]<-length(which(tya[,i]<0.05))
  powerdemeter[i]<-length(which(tyd[,i]<0.05))
  powerrsa[i]<-length(which(tyr[,i]<0.05))

  
}
b<-300
b<-1
step<-2*ncm
m<-4####the number of methods to compare
KRASpowerplot<-matrix(0,nrow=length(ncm),ncol=m+1)
KRASpowerplot[,1]<-step
KRASpowerplot[,2]<-powerplotmymethod/b
KRASpowerplot[,3]<-powerataris/b
KRASpowerplot[,4]<-powerdemeter/b
KRASpowerplot[,5]<-powerrsa/b
colnames(KRASpowerplot)<-c("samplesize","Mymethod","ATARIS","DEMETER","RSA")
KRASpowerplot<-KRASpowerplot[order(KRASpowerplot[,1]),]

KRASpowerplot2<-data.frame(samplesize=rep(KRASpowerplot[,1],m),power=c(KRASpowerplot[,2],KRASpowerplot[,3],
                           KRASpowerplot[,4],KRASpowerplot[,5]),method=rep(1:m,each=length(ncm)))
colnames(KRASpowerplot2)<-c("samplesize","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS","DEMETER","RSA"))
library(ggplot2)
p<-ggplot(KRASpowerplot2, aes(x=samplesize, y=power, group=Method)) +
  geom_line(aes(color=Method))

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS,DEMETER,RSA in Pancancer based on simulation(ES=0.5)")

#############
KRASpowerplot2f<-KRASpowerplot2[order(KRASpowerplot2$samplesize),]
KRASpowerplot2f<-KRASpowerplot2f[-1,]
head(KRASpowerplot2f)
library(ggplot2)
library(stats)
library(tibble)
lo1<-loess(power[which(KRASpowerplot2f$method==1)]~samplesize[which(KRASpowerplot2f$method==1)],data=KRASpowerplot2f,span=0.40)
lo2<-loess(power[which(KRASpowerplot2f$method==2)]~samplesize[which(KRASpowerplot2f$method==2)],data=KRASpowerplot2f,span=0.40)
lo3<-loess(power[which(KRASpowerplot2f$method==3)]~samplesize[which(KRASpowerplot2f$method==3)],data=KRASpowerplot2f,span=0.40)
lo4<-loess(power[which(KRASpowerplot2f$method==4)]~samplesize[which(KRASpowerplot2f$method==4)],data=KRASpowerplot2f,span=0.40)

smoothed101 <- predict(lo1)
smoothed102 <- predict(lo2)
smoothed103 <- predict(lo3)
smoothed104 <- predict(lo4)

plot(KRASpowerplot2f$samplesize,KRASpowerplot2f$power, type="l", main="Loess Smoothing and Prediction in 20000 simulations", xlab="Sample-size", ylab="Power")
lines(smoothed101, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)], col="red")
lines(smoothed102, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)], col="blue")
lines(smoothed103, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)], col="green")
lines(smoothed104, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)], col="black")
###
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==1)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==2)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==3)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==4)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")

lines(smoothed101, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)], col="red")
lines(smoothed102, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)], col="blue")
lines(smoothed103, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)], col="green")
lines(smoothed104, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)], col="black")
legend(100, 0.5, legend=c("Mymethod", "ATARIS","DEMETER","RSA"),
       col=c("red", "blue","green","black"), lty=1:2, cex=0.5)



#########################
ncm<-c(2,4,7,7,12,12,15,20,20,25,30,35,40,40,44,50,55,55,60,60,65)
ncw<-c(7,10,17,30,20,30,40,60,70,50,50,65,100,80,85,120,100,155,100,150,200)



powerplotmymethod<-rep(0,length(ncm))
powerataris<-rep(0,length(ncm))
powerdemeter<-rep(0,length(ncm))
powerrsa<-rep(0,length(ncm))
for(i in 1:length(ncm)){
  powerplotmymethod[i]<-length(which(resf[,i]<0.05))
  powerataris[i]<-length(which(tya[,i]<0.05))
  powerdemeter[i]<-length(which(tyd[,i]<0.05))
  powerrsa[i]<-length(which(tyr[,i]<0.05))
  
  
}
b<-200
step<-ncw+ncm
m<-4####the number of methods to compare
KRASpowerplot<-matrix(0,nrow=length(ncm),ncol=m+1)
KRASpowerplot[,1]<-step
KRASpowerplot[,2]<-powerplotmymethod/b
KRASpowerplot[,3]<-powerataris/b
KRASpowerplot[,4]<-powerdemeter/b
KRASpowerplot[,5]<-powerrsa/b

colnames(KRASpowerplot)<-c("samplesize","Mymethod","ATARIS","DEMETER","RSA")
KRASpowerplot2<-data.frame(samplesize=rep(KRASpowerplot[,1],m),power=c(KRASpowerplot[,2],KRASpowerplot[,3],
                                                                       KRASpowerplot[,4],KRASpowerplot[,5]),method=rep(1:m,each=length(ncm)))
colnames(KRASpowerplot2)<-c("samplesize","power","method")
KRASpowerplot2$Method<-factor(KRASpowerplot2$method,labels=c("OurMethod","ATARIS","DEMETER","RSA"))
library(ggplot2)
p<-ggplot(KRASpowerplot2, aes(x=samplesize, y=power, group=Method)) +
  geom_line(aes(color=Method))

p<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p+ggtitle("Power plot of our method,ATARIS,DEMETER,RSA in Pancancer based on simulation(ES=0.5)")

#############
KRASpowerplot2f<-KRASpowerplot2[order(KRASpowerplot2$samplesize),]
head(KRASpowerplot2f)
library(ggplot2)
library(stats)
library(tibble)
lo1<-loess(power[which(KRASpowerplot2f$method==1)]~samplesize[which(KRASpowerplot2f$method==1)],data=KRASpowerplot2f,span=0.40)
lo2<-loess(power[which(KRASpowerplot2f$method==2)]~samplesize[which(KRASpowerplot2f$method==2)],data=KRASpowerplot2f,span=0.40)
lo3<-loess(power[which(KRASpowerplot2f$method==3)]~samplesize[which(KRASpowerplot2f$method==3)],data=KRASpowerplot2f,span=0.40)
lo4<-loess(power[which(KRASpowerplot2f$method==4)]~samplesize[which(KRASpowerplot2f$method==4)],data=KRASpowerplot2f,span=0.40)

smoothed101 <- predict(lo1)
smoothed102 <- predict(lo2)
smoothed103 <- predict(lo3)
smoothed104 <- predict(lo4)

plot(KRASpowerplot2f$samplesize,KRASpowerplot2f$power, type="l", main="Loess Smoothing and Prediction in 20000 simulations", xlab="Sample-size", ylab="Power")
lines(smoothed101, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)], col="red")
lines(smoothed102, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)], col="blue")
lines(smoothed103, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)], col="green")
lines(smoothed104, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)], col="black")
###
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==1)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==2)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==3)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")
plot(KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)],KRASpowerplot2f$power[which(KRASpowerplot2f$method==4)], type="l", main="Power of different methods in 200 simulations", xlab="Sample-size", ylab="Power")

lines(smoothed101, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==1)], col="red")
lines(smoothed102, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==2)], col="blue")
lines(smoothed103, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==3)], col="green")
lines(smoothed104, x=KRASpowerplot2f$samplesize[which(KRASpowerplot2f$method==4)], col="black")
legend(80, 0.5, legend=c("Mymethod", "ATARIS","DEMETER","RSA"),
       col=c("red", "blue","green","black"), lty=1:2, cex=0.5)



#########################
