rm(list=ls())
library(dplyr)
##install.packages("stringr", repos='http://cran.us.r-project.org')
library("stringr")
library("dplyr")
library("tidyr")
library("lme4")
library(lmerTest)#
library(sva)
library(pamr)
library(limma)
library(data.table)
library(reshape2)
##step 1: create shRNA data based on common genes with copy number data
shRNAGeneMap = readRDS("../Input_Data/DRIVE/ShrnaGeneMap.RDS")
data2<-data.frame(shRNAGeneMap)

copyn<-read.csv("../Input_Data/Mutations/data_CNA.csv")
copyn<-data.frame(copyn)
tmp = strsplit(colnames(copyn),"[_]")
colnames(copyn) = unlist(lapply(tmp, function(x){tolower(x[1])}))

genes<-intersect(data2$GENESYMBOLS ,copyn$hugo)

shRNAdata<-data2[which(data2$GENESYMBOLS %in% genes),]

##################step2: integrating shRNA data with targetscan data(SPS feature of seed of shRNA reagents)
data = read.csv("../Input_Data/TARGETSCAN/spsthermo.csv")
spsdata<-data.frame(data)
head(data)
colnames(spsdata)<-c("Seed2","8mer","TA","TA_percentile_rank","	6mer SPS","7mer-m8 SPS"
                     ,"SPS","Mean-SPS-percentile-rank",	"miRNA","Seed") 
shRNAdata[,4]=substr(shRNAdata$SEQ, 14, 20)
sl<-length(spsdata[,1])
dd<-rep(0,sl)
seedg<-unique(shRNAdata[,4])

for(i in 1:sl){
  sg<-spsdata$Seed[i]
  dd[i]<-length(which(spsdata$Seed==sg))
  shRNAdata[which(shRNAdata[,4]==sg),5]<-spsdata$TA[which(spsdata$Seed==sg)]
  shRNAdata[which(shRNAdata[,4]==sg),6]<-spsdata$SPS[which(spsdata$Seed==sg)]
  
}
colnames(shRNAdata)<-c("SEQ","GENESYMBOLS","ENTREZ_ID","seed","TA","SPS")

##############################step3: DRIVE data based on new shRNA data
dataD = readRDS("../Input_Data/DRIVE/DriveCountData.RDS")
reagents<-shRNAdata[,"SEQ"]
cancer_data<-dataD[which(dataD$SEQ %in% reagents),]
data2<-data.frame(cancer_data)
####1) remove NA
cancer_data2 = data2[which(!is.na(data2$PLASMID_COUNT)), ]

saveRDS(cancer_data2,"../Input_Data/Processed_Data/DRIVEdata.RDS")
write.csv(shRNAdata,"../Input_Data/Processed_Data/shRNAdata.csv")
###########################################################################################
###############################################step4: creating mutational matrix

annotation = read.csv("../Input_Data/Annotation/TableS2.csv", stringsAsFactors = F)####has information on cancer specific and its cell lines
cancers2<-strsplit(annotation$PATHOLOGIST_ANNOTATION,"[:]")
annotation[,"cancer"]<-unlist(lapply(cancers2,function(x) {tolower(x[1])}))
colnames(annotation)<-c("CLEANNAME","PRIMARY_SITE","PATHOLOGIST_ANNOTATION","in.CCLE","Source","RRID","cancer") 

##
mutationtotal = read.csv("../Input_Data/Mutations/CCLE_mutations.csv")####has information on mutational status of each gene
mutationtotal<-data.frame(mutationtotal)
sampleinfo = read.csv("../Input_Data/Mutations/sample_info.csv")#####has information on cell lines
mutData<-merge(mutationtotal,sampleinfo,by="DepMap_ID")

mutData$Cellline<-tolower(mutData$stripped_cell_line_name)
#colnames(mutationtotal);#colnames(sampleinfo);#dim(mutationtotal);#dim(mutData);#length(unique(mutData$Hugo_Symbol##head(mutData);#unique(mutData$stripped_cell_line_name)

# extracting cell line name according to the annotation file (Table S2)
#tmp = strsplit(mutData$Tumor_Sample_Barcode, "[_]")
#mutData$Cellline = unlist(lapply(tmp, function(x){tolower(x[1])}))

# select only cell lines in the DRIVE project
mutData = mutData[which(mutData$Cellline %in% annotation$CLEANNAME), ]

mydata = mutData[, c("Hugo_Symbol","Variant_Classification","Cellline")]

cls<-unique(mydata$Cellline)
################################################
cancer_cellines = annotation[, "CLEANNAME"]

##step3 : create missense mutation matrix
select_submatrix <- function(data, celllines, threshold = 1) {
  nData = data[celllines, ]
  nData[, which(apply(nData, 2, sum) >= threshold)]
}




###########################################################

mutation = names(sort(table(mydata$Variant_Classification)))#####different types of mutations

mydata$vv<-0
mydata[which(mydata$Variant_Classification %in% c("Start_Codon_Ins","Stop_Codon_Ins","Start_Codon_Del",     
                                                  "Stop_Codon_Del","In_Frame_Ins","In_Frame_Del","Splice_Site"
                                                  ,"Frame_Shift_Ins","Nonsense_Mutation","Frame_Shift_Del"
                                                  )),"vv"]<-1
mydata[which(mydata$Variant_Classification %in% "Missense_Mutation"),"vv"]<-2

mydata[which(is.na(mydata$vv)),"vv"]<-0
mydata2<-mydata[,c("Hugo_Symbol","Cellline","vv")]
#############################################
mutation_matix=list()
vv<-sort(unique(mydata$vv))
for(i in 1:length(vv)){
  tmpMutData <-mydata2[mydata2$vv==vv[i],]
  mutation_matix[[i]]=as.matrix(xtabs(~ tmpMutData$Cellline + tmpMutData$Hugo_Symbol))
  dimnames(mutation_matix[[i]]) = list(dimnames(mutation_matix[[i]])[[1]], dimnames(mutation_matix[[i]])[[2]])
}
select_submatrix <- function(data, celllines, threshold = 1) {
  nData = data[celllines, ]
  nData[, which(apply(nData, 2, sum) >= threshold)]
}

# prepare missense data
mutation_matix[[3]][mutation_matix[[3]]>=1]=1
missenceData=mutation_matix[[3]]

cancer_cellines3 = intersect(annotation[, "CLEANNAME"], 
                            rownames(mutation_matix[[3]]))


# build sub matrix for the mutation data
data = missenceData
celllines = cancer_cellines3
threshold = 1

cancer_mutations3 = select_submatrix(missenceData, cancer_cellines3, threshold = 1)##########this is used for identify mutant cell lines

####the two other mutation matrices are used to make the wildtype matrix
mutation_matix[[2]][mutation_matix[[2]]>=1]=1
otherData=mutation_matix[[2]]

cancer_cellines2 = intersect(annotation[, "CLEANNAME"], 
                             rownames(mutation_matix[[2]]))
data = otherData
celllines = cancer_cellines2
threshold = 1

cancer_mutations2 = select_submatrix(otherData, cancer_cellines2, threshold = 1)

mutation_matix[[1]][mutation_matix[[1]]>=1]=1
zeroData=mutation_matix[[1]]

cancer_cellines1 = intersect(annotation[, "CLEANNAME"], 
                             rownames(mutation_matix[[1]]))
data = zeroData
celllines = cancer_cellines1
threshold = 1

cancer_mutations0 = select_submatrix(zeroData, cancer_cellines1, threshold = 1)

######################Copy number data
mygenes<-sort(genes)###genes are intersect shRNA genes and copy number genes
mycelline<-sort(unique(annotation$CLEANNAME))

copyn<-read.csv("../Input_Data/Mutations/data_CNA.csv")
copyn<-data.frame(copyn)
tmp = strsplit(colnames(copyn),"[_]")
colnames(copyn) = unlist(lapply(tmp, function(x){tolower(x[1])}))

copynf<-copyn[which(copyn$hugo %in% genes),-c(1,2)]
rownames(copynf)<-copynf$hugo
copynmat<-t(copynf)
copynmat<-as.matrix(copynmat)
##########################
wild_matix=matrix(0,nrow=length(mygenes),ncol=length(mycelline))
rownames(wild_matix)<-mygenes
colnames(wild_matix)<-mycelline
cancell<-mycelline
for(i in 1:length(mygenes)){
  cg<-copynmat[,which(colnames(copynmat)==mygenes[i])]
  int12<-intersect(names(which(cancer_mutations2[,which(colnames(cancer_mutations2)==mygenes[i])]==0)),
                     names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes[i])]==0)))
  int123<-intersect(int12,names(which(cancer_mutations0[,which(colnames(cancer_mutations0)==mygenes[i])]==0)))
  cell2<-intersect(int123,names(which(abs(cg)!=2)))
  cel01<-sort(c(names(which(cancer_mutations2[,which(colnames(cancer_mutations2)==mygenes[i])]!=0)),
            names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes[i])]!=0)),
            names(which(cancer_mutations0[,which(colnames(cancer_mutations0)==mygenes[i])]!=0)),
            names(which(abs(cg)==2))))
dif1<-cancel[-which(cancell %in% cel01)]
ds<-sort(c(names(which(cancer_mutations2[,which(colnames(cancer_mutations2)==mygenes[i])]!=0)),
      names(which(cancer_mutations0[,which(colnames(cancer_mutations0)==mygenes[i])]!=0)),
      names(which(abs(cg)==2))))
dif2<-setdiff(cancell,dif1)
  wild_matix[which(rownames(wild_matix)==mygenes[i]),which(colnames(wild_matix) %in% dif1)]<-1
  wild_matix[which(rownames(wild_matix)==mygenes[i]),which(colnames(wild_matix) %in% dif2)]<-0
  wild_matix[which(rownames(wild_matix)==mygenes[i]),which(colnames(wild_matix) %in% ds)]<-2
  
}
#############

saveRDS(wild_matix,"../Input_Data/Processed_Data/wild_matix.R")
saveRDS(cancer_mutations3,"../Input_Data/Processed_Data/cancer_mutations.R")


