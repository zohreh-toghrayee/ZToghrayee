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

 
#################3
nrRepeat = 1000
genes = c("KRAS","NRAS","BRAF","TP53","UGT8","TRPV3","FUS")
dat = data.frame(matrix(NA, nrRepeat * length(genes), 4))
nr_genes_wt = rep(0,length(genes))
nr_genes_mut = rep(0,length(genes))
names(nr_genes_wt) = names(nr_genes_mut ) = genes

colnames(dat) = c("gene", "nr_wt", "nr_mut", "nr_tot")
for(j in 1:length(genes)) {
cancerf<-cancerd[which(cancerd$genename==genes[j]),]
wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==genes[j]),]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==genes[j])]==1)),unique(cancerf$CLEANNAME))###cellines where the ith gene is mutated
nr_genes_wt[j] = length(wildtype)
nr_genes_mut[j]= length(mutant)
}


dat = data.frame(matrix(NA, nrRepeat , 4))
colnames(dat) = c("gene", "nr_wt", "nr_mut", "nr_tot")


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

dat$nr_wt<-as.numeric(dat$nr_wt)
dat$nr_mut<-as.numeric(dat$nr_mut)
dat$nr_tot<-as.numeric(dat$nr_tot)
SIMT<-data.frame(dat)
###############################

write.csv(SIMT,"SIMTTypeone.csv")
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 

cellinesf<-unique(cancerd$CLEANNAME)

##########################
n<-dim(SIMT)[1]
mymethod<-list()

##for(j in 1:length(genes)){
  tableg<-SIMT[which(SIMT$gene %in% "BRAF"),]
  cancer<-cancerd[which(cancerd$genename %in% "BRAF"),]
  wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)=="BRAF"),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
  mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)=="BRAF")]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated
  cancer$muts<-0
  cancer$muts[which(cancer$CLEANNAME %in% wildtype)]<-0
  cancer$muts[which(cancer$CLEANNAME %in% mutant)]<-1
   
  celinesf<-unique(cancer$CLEANNAME)
  num<-c(1:dim(tableg)[1])
twp= foreach(i=num, .combine="rbind", .packages=c("glmmTMB","lme4")) %dopar%  {

  twp2<-data.frame() 

      cancer$mutationgroup<-0
nm<-tableg$nr_mut[i]
nw<-tableg$nr_wt[i]
      am<-sample(celinesf,nm)
      aw<-sample(celinesf[-which(celinesf %in% am)],nw)
      cancer[which(cancer$CLEANNAME %in% am),"mutationgroup"]<-1
      cancer[which(cancer$CLEANNAME %in% aw),"mutationgroup"]<-0
      mnb <-summary(glmmTMB(rank ~mutationgroup+(1|CLEANNAME),data=cancer,family="nbinom2"))
     cv<-mnb$coefficients
     twp2[1,1]<- pnorm(as.matrix(cv$cond)[2,3])
  
   
twp2
}
mymethod[[1]]<-twp
###}

typeone<-matrix(0,nrow=nrRepeat,ncol=2)
##for(j in 1:length(genes)){
typeone[,2]<-mymethod[[1]]$V1
##}
write.csv(typeone,"6.typeoneSIM1000BRAF.csv")







stopCluster(myCluster)
