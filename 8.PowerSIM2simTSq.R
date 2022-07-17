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
ATARIS<-readRDS("DRIVE_ATARiS_data.RDS")###ATARIS results
tmp = strsplit(colnames(ATARIS), "[_]")
colnames(ATARIS) = unlist(lapply(tmp, function(x){tolower(x[1])}))

DEMETER2<-read.csv("Demeter2_DRIVE_gene_dep_scores.csv")###demeter results
tmpd = strsplit(colnames(DEMETER2), "[_]")
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
DEMETER2[,1] = word(DEMETER2[,1], 1)
colnames(DEMETER2) = unlist(lapply(tmpd, function(x){tolower(x[1])}))

RSA<-readRDS("DRIVE_RSA_data.RDS")###RSA results
tmpd = strsplit(colnames(RSA), "[_]")
colnames(RSA) = unlist(lapply(tmpd, function(x){tolower(x[1])}))
 genedata<-cancerdd[which(cancerdd$genename %in% c("KRAS","NRAS","CTNNB1")), ]
 cancerd<-data.frame(genedata)

 gene<-c("KRAS","NRAS","CTNNB1")
num<-1:1000##number of simulation
for( i in 1:length(gene)){
cancer<-cancerd[which(cancerd$genename %in% gene[i]),]
wildtype2=intersect(names(which(wild_matix[which(rownames(wild_matix)==gene[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==gene[i])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
cellinesf<-unique(cancer$CLEANNAME)
wildtype<-cellinesf[-which(cellinesf %in% mutant)]

ataris=matrix(0,nrow=length(which(colnames(ATARIS) %in% cellinesf)),ncol=4)
colnames(ataris)<-c("score","group")
rownames(ataris)=colnames(ATARIS[which(rownames(ATARIS)==gene[i]),which(colnames(ATARIS) %in% cellinesf)])
ataris[,1]= t(ATARIS[which(rownames(ATARIS)==gene[i]),which(colnames(ATARIS) %in% cellinesf)])

demeter=matrix(0,nrow=length(which(colnames(DEMETER2) %in% cellinesf)),ncol=2)
colnames(demeter)<-c("score","group")
rownames(demeter)=colnames(DEMETER2[which(rownames(DEMETER2)==gene[i]),which(colnames(DEMETER2) %in% cellinesf)])
demeter[,1]= t(DEMETER2[which(DEMETER2$x==gene[i]),which(colnames(DEMETER2) %in% cellinesf)])

rsa=matrix(0,nrow=length(which(colnames(RSA) %in% cellinesf)),ncol=2)
colnames(rsa)<-c("score","group")
rownames(rsa)=colnames(RSA[which(rownames(RSA)==gene[i]),which(colnames(RSA) %in% cellinesf)])
rsa[,1]= t(RSA[which(rownames(RSA)==gene[i]),which(colnames(RSA) %in% cellinesf)])

######33
mymethod<-list()
nww<-length(wildtype)
nmm<-length(mutant)

ncm<-c(2:12)
coef<-c(1:5)
for(q in length(coef)){
ncw<-coef[q]*ncm

  for(j in 1:length(ncm)){
  nw<-ncw[j]
  nm<-ncm[j]
    twp<-data.frame()
    twp= foreach(k=num, .combine="rbind", .packages=c("glmmTMB","lme4")) %dopar%  {
    twp2<-data.frame()
       
    cancer$mutshufle<-0
    cancer$rankv<-0
    WW<-sample(wildtype,nw)
    MM<-sample(mutant,nm)
    cancerg<-cancer[which(cancer$CLEANNAME %in% c(WW,MM)),]
    cancerg[which(cancerg$CLEANNAME %in% WW),"mutshufle"]<-0
    cancerg[which(cancerg$CLEANNAME %in% MM),"mutshufle"]<-1
          
    mnb <-summary(glmmTMB(rank ~mutshufle+(1|CLEANNAME),data=cancerg,family="nbinom2"))
    cv<-mnb$coefficients
    twp2[1,1]<- pnorm(as.matrix(cv$cond)[2,3])

     ataris2<-ataris[which(rownames(ataris) %in% c(WW,MM)),]
    ataris2[which(rownames(ataris2) %in% WW),2]<-0
    ataris2[which(rownames(ataris2) %in% MM),2]<-1
    df<-summary(lm(ataris2[,1]~ataris2[,2]))
    df1<-lm(ataris2[,2]~ataris2[,1])
    andf<-anova(df1)
    dff<-andf$Df[2]
    twp2[1,2]<- pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)
 
    demetrer2<-demeter[which(rownames(demeter) %in% c(WW,MM)),]
    demetrer2[which(rownames(demetrer2) %in% WW),2]<-0
    demetrer2[which(rownames(demetrer2) %in% MM),2]<-1
    df<-summary(lm(demetrer2[,1]~demetrer2[,2]))
    df1<-lm(demetrer2[,1]~demetrer2[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]
    twp2[1,3]<- pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)

    rsa2<-rsa[which(rownames(rsa) %in% c(WW,MM)),]
    rsa2[which(rownames(rsa2) %in% WW),2]<-0
    rsa2[which(rownames(rsa2) %in% MM),2]<-1
    df<-summary(lm(rsa2[,1]~rsa2[,2]))
    df1<-lm(rsa2[,4]~rsa2[,2])
    andf<-anova(df1)
    dff<-andf$Df[2]     
    twp2[1,4]<- pt(as.matrix(df$coefficients)[2,3],dff, lower.tail = TRUE)

twp2
}

print(j)
mymethod[[j]]<-twp
}

 powerm<-as.data.frame(do.call("rbind",mymethod))##all for one gene and one q (all 11 ncm each 1000 times)

write.csv(powerm,paste0("8.powermSIM2Sq",q,genes[i],".csv").csv")
}
}


##########################

stopCluster(myCluster)