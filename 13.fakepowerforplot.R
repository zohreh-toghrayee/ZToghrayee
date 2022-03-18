rm(list=ls())
set.seed(10)

nrRepeat = 10000
genes = c("KRAS", "TP53", "NRAS")
dat1 = data.frame(matrix(NA, nrRepeat , 4))
dat2 = data.frame(matrix(NA, nrRepeat , 4))
dat3 = data.frame(matrix(NA, nrRepeat , 4))


nr_genes_wt = c(180, 200, 190)
nr_genes_mut = c(72, 100, 50)
names(nr_genes_wt) = names(nr_genes_mut ) = genes

k = 1
colnames(dat1) = c("gene", "nr_wt", "nr_mut", "nr_tot")
colnames(dat2) = c("gene", "nr_wt", "nr_mut", "nr_tot")
colnames(dat3) = c("gene", "nr_wt", "nr_mut", "nr_tot")

##for(i in 1:3) {
  for(k in 1:nrRepeat) {
    nr_wt = sample(2:nr_genes_wt[1], 1)
    nr_mut = sample(2:min(nr_genes_wt[1]+1, nr_genes_mut[1]), 1)
    dat1[k, ] = c(genes[1], nr_wt, nr_mut, nr_wt + nr_mut)
    
  }
for(k in 1:nrRepeat) {
  nr_wt = sample(2:nr_genes_wt[2], 1)
  nr_mut = sample(2:min(nr_genes_wt[2]+1, nr_genes_mut[2]), 1)
  dat2[k, ] = c(genes[2], nr_wt, nr_mut, nr_wt + nr_mut)
  
}
for(k in 1:nrRepeat) {
  nr_wt = sample(2:nr_genes_wt[3], 1)
  nr_mut = sample(2:min(nr_genes_wt[3]+1, nr_genes_mut[3]), 1)
  dat3[k, ] = c(genes[3], nr_wt, nr_mut, nr_wt + nr_mut)
  
}
#}
datt<-rbind(dat1,dat2)
myfile<-rbind(datt,dat3)
head(myfile)
myfile$nr_tot<-as.numeric(myfile$nr_tot)
myfile$nr_wt<-as.numeric(myfile$nr_wt)
myfile$nr_mut<-as.numeric(myfile$nr_mut)
summary(myfile$nr_tot)
# fake p-values for type I error
myfile$pvalue = runif(nrow(myfile))
head(myfile)

# plotting
##s for type I error
plot(myfile$nr_tot,  myfile$pvalue)

# plotting
KRAS<-myfile[which(myfile$gene %in% "NRAS"),]
length(which(KRAS$pvalue<0.05))
dim(KRAS)
106/2000
plot(KRAS$nr_tot,KRAS$pvalue)


rm(list=ls())
set.seed(10)

nrRepeat = 10000
genes = c("KRAS", "TP53", "NRAS")
dat = data.frame(matrix(NA, nrRepeat * length(genes), 4))


nr_genes_wt = c(180, 200, 190)
nr_genes_mut = c(72, 100, 50)
names(nr_genes_wt) = names(nr_genes_mut ) = genes

k = 1
colnames(dat) = c("gene", "nr_wt", "nr_mut", "nr_tot")
for(gene in genes) {
  for(i in 1:nrRepeat) {
    nr_wt = sample(2:nr_genes_wt[gene], 1)
    nr_mut = sample(2:min(nr_genes_wt[gene]+1, nr_genes_mut[gene]), 1)
    dat[k, ] = c(gene, nr_wt, nr_mut, nr_wt + nr_mut)
    k = k + 1
  }
  
}
head(dat)
summary(dat$nr_tot)
dat[1500:1505,]
# fake p-values for type I error
dat$pvalue = runif(nrow(dat))
head(dat)

dat$nr_tot<-as.numeric(dat$nr_tot)

library(ggplot2)
mygene<-dat[which(dat$gene %in% "KRAS"),]
lu<-length(unique(mygene$nr_tot))
su<-unique(mygene$nr_tot)
tu<-as.data.frame(table(mygene$nr_tot))
head(tu)
tu$Freq[which(tu$Var1==su[i])]
five<-data.frame(samplesize=su,typeone=rep(0,lu))
for(i in 1:lu){
  five[i,1]<-su[i]
  five[i,2]<-(length(which(mygene$pvalue[which(mygene$nr_tot==su[i])]<0.05))/tu$Freq[which(tu$Var1==su[i])])
}


##five$samplesize<-factor(five$samplesize)

ggplot(five,aes(x=samplesize,y=typeone))+geom_line()
mygene<-dat
lu<-length(unique(mygene$nr_tot))
su<-unique(mygene$nr_tot)
tu<-as.data.frame(table(mygene$nr_tot))
head(tu)
tu$Freq[which(tu$Var1==su[i])]
five<-data.frame(samplesize=su,typeone=rep(0,lu))
for(i in 1:lu){
  five[i,1]<-su[i]
  five[i,2]<-(length(which(mygene$pvalue[which(mygene$nr_tot==su[i])]<0.05))/tu$Freq[which(tu$Var1==su[i])])
}

##five$samplesize<-factor(five$samplesize)
ggplot(five,aes(x=samplesize,y=typeone))+geom_boxplot()
mean(five$typeone)
five2<-tibble(five)
five2<-five2[order(five2$samplesize),]
library(stats)
library(tibble)
lo1<-loess(typeone~samplesize,data=five2,span=0.1)
lo2<-loess(typeone~samplesize,data=five2,span=0.4)
lo3<-loess(typeone~samplesize,data=five2,span=0.6)

smoothed101 <- predict(lo1)
smoothed102 <- predict(lo2)
smoothed103 <- predict(lo3)

plot(five2$samplesize,five2$typeone, type="l", main="Loess Smoothing and Prediction", xlab="Date", ylab="Unemployment (Median)")
lines(smoothed101, x=five2$samplesize, col="red")
lines(smoothed102, x=five2$samplesize, col="blue")
lines(smoothed103, x=five2$samplesize, col="green")

method = c("loess", "model.frame")
x<-rnorm(40)
y=2*x+rnorm(40)
lo<-loess(y~x,span = 0.5)
smoothed10 <- predict(lo)
plot(y,x, type="l", main="Loess Smoothing and Prediction", xlab="Date", ylab="Unemployment (Median)")
lines(smoothed10, x=x, col="red")

calcSSE <- function(x){
  loessMod <- try(loess(uempmed ~ index, data=economics, span=x), silent=T)
  res <- try(loessMod$residuals, silent=T)
  if(class(res)!="try-error"){
    if((sum(res, na.rm=T) > 0)){
      sse <- sum(res^2)  
    }
  }else{
    sse <- 99999
  }
  return(sse)
}

# Run optim to find span that gives min SSE, starting at 0.5
optim(par=c(0.5), calcSSE, method="SANN")
#> $par