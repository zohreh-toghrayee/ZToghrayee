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
library(glmmTMB)
myCluster <-makeCluster(20, type = "PSOCK")
registerDoParallel(myCluster)
#######################END OF MUTATIONAL MATRIX
cancerdd = as.data.frame(read.csv("4.finaldata.csv"))
wild_matix<-readRDS("wild_matix.R")
cancer_mutations3<-readRDS("cancer_mutations.R")

 KRASk<-cancerdd[which(cancerdd$genename %in% c("KRAS","NRAS","BRAF","TP53","UGT8","TRPV3","FUS")), ]
 cancerd<-data.frame(KRASk)

 

##KRASk =as.data.frame(fread("KRASdata.csv"))
cancer<-cancerd[which(cancerd$genename %in% "KRAS"),]
wildtype2=intersect(names(which(wild_matix[which(rownames(wild_matix)=="KRAS"),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)=="KRAS")]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
cellinesf<-unique(cancer$CLEANNAME)
wildtype<-cellinesf[-which(cellinesf %in% mutant)]
nw<-length(wildtype)
nm<-length(mutant)
q<-round(nw/nm)
ncm<-c(2:12)
ncw<-q*ncm

ATARIS<-readRDS("DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))
ataris=matrix(0,nrow=length(which(colnames(ATARIS) %in% cellinesf)),ncol=4)
colnames(ataris)<-c("score","group","score2","samescale")
ataris[,3]<-0
rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS)=="KRAS"),which(colnames(ATARIS) %in% cellinesf)])
ataris[,1]= t(ATARIS[which(rownames(ATARIS)=="KRAS"),which(colnames(ATARIS) %in% cellinesf)])
dim(ataris)

#####################Deemetr
DEMETER2<-read.csv("Demeter2_DRIVE_gene_dep_scores.csv")###demeter results
tmpd = strsplit(colnames(DEMETER2), "[_]")
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
DEMETER2[,1] = word(DEMETER2[,1], 1)
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
##unique(DEMETER2$x)
demeter=matrix(0,nrow=length(which(colnames(DEMETER2) %in% cellinesf)),ncol=4)

colnames(demeter)<-c("score","group","score2","samescale")
demeter[,3]<-0
rownames(demeter)=colnames(DEMETER2[which(rownames(DEMETER2)=="KRAS"),which(colnames(DEMETER2) %in% cellinesf)])
demeter[,1]= t(DEMETER2[which(DEMETER2$x=="KRAS"),which(colnames(DEMETER2) %in% cellinesf)])
dim(demeter)
####################333

RSA<-readRDS("DRIVE_RSA_data.RDS")###RSA results
tmpd = strsplit(colnames(RSA), "[_]")
colnames(RSA) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
rsa=matrix(0,nrow=length(which(colnames(RSA) %in% cellinesf)),ncol=4)

colnames(rsa)<-c("score","group","score2","samescale")
#demeter[,3]<-0
rownames(rsa)=colnames(RSA[which(rownames(RSA)=="KRAS"),which(colnames(RSA) %in% cellinesf)])
rsa[,1]= t(RSA[which(rownames(RSA)=="KRAS"),which(colnames(RSA) %in% cellinesf)])
dim(rsa)
d<-0.5
effects<-0
##############33
ataris=matrix(0,nrow=length(which(colnames(ATARIS) %in% cellinesf)),ncol=4)
  colnames(ataris)<-c("score","group","score2","samescale")
  ataris[,3]<-0
  rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS) %in% "KRAS"),which(colnames(ATARIS) %in%  cellinesf)])
  ataris[,1]= t(ATARIS[which(rownames(ATARIS) %in% "KRAS"),which(colnames(ATARIS) %in%  cellinesf)])
g<-dim(ataris)[1]
srank<-summary(cancer$rank)
atarisor<-ataris[order(ataris[,1]),]
##dim(atarisor)
atarisor[1,3]<-srank[1]
atarisor[g,3]<-srank[6]
dif<-(atarisor[g,1]-atarisor[1,1])
dh<-dim(atarisor)[1]
rot<-(srank[6]-srank[1])
for(i in 2:dh){
  atarisor[i,3]<- (rot*(atarisor[i,1]-atarisor[1,1])/dif)+2
}

demeter=matrix(0,nrow=length(which(colnames(DEMETER2) %in% cellinesf)),ncol=4)
  
  colnames(demeter)<-c("score","group","score2","samescale")
  demeter[,3]<-0
  rownames(demeter)=colnames(DEMETER2[which(rownames(DEMETER2) %in% "KRAS"),which(colnames(DEMETER2) %in% cellinesf)])
  demeter[,1]= t(DEMETER2[which(DEMETER2$x %in% "KRAS"),which(colnames(DEMETER2) %in% cellinesf)])
  ##dim(demeter)
g<-dim(demeter)[1]
srank<-summary(cancer$rank)
demeteror<-demeter[order(demeter[,1]),]
##dim(demeteror)
demeteror[1,3]<-srank[1]
demeteror[g,3]<-srank[6]
dif<-(demeteror[g,1]-demeteror[1,1])
dh<-dim(demeteror)[1]
rot<-(srank[6]-srank[1])
for(i in 2:dh){
  demeteror[i,3]<- (rot*(demeteror[i,1]-demeteror[1,1])/dif)+2
}
rsa=matrix(0,nrow=length(which(colnames(RSA) %in% cellinesf)),ncol=4)
  colnames(rsa)<-c("score","group","score2","samescale")
  
  rownames(rsa)=colnames(RSA[which(rownames(RSA) %in% "KRAS"),which(colnames(RSA) %in% cellinesf)])
  rsa[,1]= t(RSA[which(rownames(RSA) %in% "KRAS"),which(colnames(RSA) %in% cellinesf)])
  dim(rsa)
g<-dim(rsa)[1]
srank<-summary(cancer$rank)
rsaor<-rsa[order(rsa[,1]),]
#dim(rsaor)
rsaor[1,3]<-srank[1]
rsaor[g,3]<-srank[6]
dif<-(rsaor[g,1]-rsaor[1,1])
dh<-dim(rsaor)[1]
rot<-(srank[6]-srank[1])
for(i in 2:dh){
  rsaor[i,3]<- (rot*(rsaor[i,1]-rsaor[1,1])/dif)+2
}


###############33
mymethod<-list()
nww<-length(wildtype)
nmm<-length(mutant)
q<-round(nww/nmm)
q<-1
##ncm<-seq(from=2,to=24,by=2)
ncm<-c(2:12)
ncw<-q*ncm
num<-1:500##number of simulation
     for(j in 1:length(ncm)){
  nw<-ncw[j]
  nm<-ncm[j]
    twp<-data.frame()
    twp= foreach(i=num, .combine="rbind", .packages=c("glmmTMB","lme4")) %dopar%  {
    twp2<-data.frame()
       
    cancer$mutshufle<-0
    cancer$rankv<-0
    WW<-sample(wildtype,nw)
    MM<-sample(mutant,nm)
    cancerg<-cancer[which(cancer$CLEANNAME %in% c(WW,MM)),]
    cancerg[which(cancerg$CLEANNAME %in% WW),"mutshufle"]<-0
    cancerg[which(cancerg$CLEANNAME %in% MM),"mutshufle"]<-1
    cancerg[which(cancerg$CLEANNAME %in% WW),"rankv"]<-cancerg$rank[which(cancerg$CLEANNAME %in% WW)]+round(effects)
    cancerg[which(cancerg$CLEANNAME %in% MM),"rankv"]<-cancerg$rank[which(cancerg$CLEANNAME %in% MM)]
      
     mnb <-summary(glmmTMB(rankv ~mutshufle+(1|CLEANNAME),data=cancerg,family="nbinom2"))
     cv<-mnb$coefficients
     twp2[1,1]<- pnorm(as.matrix(cv$cond)[2,3])

     ataris2<-atarisor[which(rownames(atarisor) %in% c(WW,MM)),]
    ataris2[which(rownames(ataris2) %in% WW),2]<-0
    ataris2[which(rownames(ataris2) %in% MM),2]<-1
    atarisf<-ataris2[which(rownames(ataris2) %in% c(WW,MM)),]
    atarisf[which(rownames(atarisf) %in% WW),4]<-atarisf[which(rownames(atarisf) %in% WW),3]+round(effects)
    atarisf[which(rownames(atarisf) %in% MM),4]<-atarisf[which(rownames(atarisf) %in% MM),3]
    df<-summary(lm(atarisf[,4]~atarisf[,2]))
    df1<-lm(atarisf[,4]~atarisf[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]
     twp2[1,2]<- pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
demetrer2<-demeteror[which(rownames(demeteror) %in% c(WW,MM)),]
    demetrer2[which(rownames(demetrer2) %in% WW),2]<-0
    demetrer2[which(rownames(demetrer2) %in% MM),2]<-1
    demetrerf<-demetrer2[which(rownames(demetrer2) %in% c(WW,MM)),]
    demetrerf[which(rownames(demetrerf) %in% WW),4]<-demetrerf[which(rownames(demetrerf) %in% WW),3]+round(effects)
    demetrerf[which(rownames(demetrerf) %in% MM),4]<-demetrerf[which(rownames(demetrerf) %in% MM),3]
    df<-summary(lm(demetrerf[,4]~demetrerf[,2]))
    df1<-lm(demetrer2[,4]~demetrer2[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]
     twp2[1,3]<- pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)

rsa2<-rsaor[which(rownames(rsaor) %in% c(WW,MM)),]
    rsa2[which(rownames(rsa2) %in% WW),2]<-0
    rsa2[which(rownames(rsa2) %in% MM),2]<-1
    rsaf<-rsa2[which(rownames(rsa2) %in% c(WW,MM)),]
    rsaf[which(rownames(rsaf) %in% WW),4]<-rsaf[which(rownames(rsaf) %in% WW),3]+round(effects)
    rsaf[which(rownames(rsaf) %in% MM),4]<-rsaf[which(rownames(rsaf) %in% MM),3]
    df<-summary(lm(rsaf[,4]~rsaf[,2]))
    df1<-lm(rsaf[,4]~rsaf[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]     
    twp2[1,4]<- pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)

twp2
}

print(j)
mymethod[[j]]<-twp
}

 powerm<-as.data.frame(do.call("rbind",mymethod))



write.csv(powerm,"8.powermSIM2Sq11.csv")
print("mymethod")

##########################
##saveRDS(mymethod, file="powermymethod.RData")
stopCluster(myCluster)