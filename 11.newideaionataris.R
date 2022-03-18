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

ATARIS<-readRDS("../Input_Data/DRIVE/DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))
###rank data in each cell line
head(ATARIS)
rownames(ATARIS)##genes
colnames(ATARIS)###Celllines
atarisrank<-data.frame(ATARIS)
rownames(atarisrank)<-rownames(ATARIS)
colnames(atarisrank)<-colnames(ATARIS)
for(i in 1:length(rownames(ATARIS))){
  atarisrank[i,]<-rank(ATARIS[i,])
}
atarisrank[1:3,1:4]
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
gene<-"KRAS"
wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==gene),]==1)),unique(colnames(ATARIS)))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==gene)]==1)),unique(colnames(ATARIS)))###cellines where the ith gene is mutated
# others<-intersect(names(which(wild_matix[which(rownames(wild_matix)==g
cancelf<-unique(colnames(ATARIS))
mygenes2<-rownames(ATARIS)
dim(ATARIS)
ataris=matrix(0,nrow=length(colnames(ATARIS) ),ncol=2)
#if(length(which(rownames(ATARIS)==mygenes2[i])!=0)){
  rownames(ataris)=colnames(atarisrank[which(rownames(atarisrank)==gene),which(colnames(atarisrank) %in% cancelf)])
  ataris[,1]= t(atarisrank[which(rownames(atarisrank)==gene),which(colnames(atarisrank) %in% cancelf)])
  ataris[which(rownames(ataris) %in% mutant),2]<-1
  ataris[which(rownames(ataris) %in% wildtype),2]<-0
  #if(length(table(ataris[,2]))>1){
    dfr1<-summary(lm(ataris[,1]~ataris[,2]))
    df<-lm(ataris[,1]~ataris[,2])
    andf1<-anova(df)
    dff1<-andf1$Df[2]
    pt(as.matrix(dfr1$coefficients)[2,3],dff1)
    
#  } 
#}
ataris<-data.frame(ataris)
colnames(ataris)<-c("rank","muts")
agata<-data.frame(aggregate(ataris$rank,by=list(ataris$muts),sum))
agata2<-data.frame(aggregate(ataris$rank,by=list(ataris$muts),mean))
#agata3<-data.frame(aggregate(ataris$rank,by=list(ataris$muts),sd))
lamda<-mean(ataris$rank)
ataris$pois<--log(ppois(ataris$rank,lamda))
ataris$gam<--log(pgamma(ataris$rank,lamda))
ataris<-ataris[which(is.finite(ataris$gam)),]
dfr1<-summary(lm(ataris$gam~ataris$muts))
df<-lm(ataris$gam~ataris$muts)
df<-lm(ataris$pois~ataris$muts)
andf1<-anova(df)
andf1$`Pr(>F)`[1]
andf1<-anova(df)
dff1<-andf1$Df[2]
pt(as.matrix(dfr1$coefficients)[2,3],dff1)
###############for gamma
atarisg<-ataris[-57,]##inf
df<-lm(atarisg$gam~atarisg$muts)
andf1<-anova(df)
andf1$`Pr(>F)`
dff1<-andf1$Df[2]
#pt(as.matrix(dfr1$coefficients)[2,4],dff1)
#as.matrix(andf1$coefficients)[2,5]
####
ataris$rank2<-round(ataris$rank)
ataris$celline<-rownames(ataris)
library(glmmTMB)
mnb<-summary(glmmTMB(rank2 ~ muts + (1|celline), data=ataris, family="poisson"))
cv<-mnb$coefficients
 pnorm(as.matrix(cv$cond)[2,3])
 ataris$muts2<-factor(ataris$muts)
# mnb<-summary(glmmTMB(pois ~ muts + (1|celline), data=ataris, family=gaussian))
 cv<-mnb$coefficients
 pnorm(as.matrix(cv$cond)[2,3])
 mnb<-lmer(pois ~ muts2 , data=ataris)
 cv<-mnb$coefficients
 pnorm(as.matrix(cv$cond)[2,3])
table(ataris$muts)
table(ataris$celline)
####
colnames(agata)<-c("group","ranks")

size.vec <- 2:4 
set.seed(250) 
dat <- rnbinom(3, size = size.vec, prob = 0.2) 
dat 
enbinom(dat, size = size.vec, method = "mvue") 

n = 500  # Example 2: simulated data
x = runif(n)
y1 = rnbinom(n, mu=exp(3+x), size=exp(1)) # k is size
y2 = rnbinom(n, mu=exp(2-x), size=exp(0))
library("VGAM")
fit = vglm(cbind(y1,y2) ~ x, negbinomial, tra=TRUE) # multivariate response
coef(fit, matrix=TRUE)
#######################################33
##################################333
gene<-"KRAS"
if(length(which(rownames(ATARIS)==mygenes2[i])!=0)){
  rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS)==gene),which(colnames(ATARIS) %in% cancelf)])
  ataris[,1]= t(ATARIS[which(rownames(ATARIS)==gene),which(colnames(ATARIS) %in% cancelf)])
  ataris[which(rownames(ataris) %in% mutant),2]<-1
  ataris[which(rownames(ataris) %in% wildtype),2]<-0
  if(length(table(ataris[,2]))>1){
    dfr1<-summary(lm(ataris[,1]~ataris[,2]))
    df<-lm(ataris[,1]~ataris[,2])
    andf1<-anova(df)
    dff1<-andf1$Df[2]
   pt(as.matrix(dfr1$coefficients)[2,3],dff1)
    
  } 
}

aggregate(ataris[])