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

#######################################################################3
cancerd = as.data.frame(fread("../Processed_Data/FewGenesData.csv"))
head(cancerd)
wild_matix<-readRDS("../Processed_Data/wild_matix.R")
cancer_mutations3<-readRDS("../Processed_Data/cancer_mutations.R")
###################
mygenes<-unique(cancerd$genename)
cancerd$muts<-0
cancerg<-list()
for(i in 1:length(mygenes)){
  cancer<-cancerd[which(cancerd$genename %in% mygenes[i]),]
wildtype=intersect(names(which(wild_matix[which(rownames(wild_matix)==mygenes[i]),]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is not mutated
mutant=intersect(names(which(cancer_mutations3[,which(colnames(cancer_mutations3)==mygenes[i])]==1)),unique(cancer$CLEANNAME))###cellines where the ith gene is mutated

cancer$muts[which(cancer$CLEANNAME %in% wildtype & cancer$genename==mygenes[i])]<-0
cancer$muts[which(cancer$CLEANNAME %in% mutant & cancer$genename==mygenes[i])]<-1
cancerg[[i]]<-cancer
}
dataset1<-do.call(rbind, cancerg)
write.csv(dataset1,"../Processed_Data/afewgenesdataf.csv")
head(dataset1)
dim(dataset1)
dim(cancerd)

dataset1$mutsF<-factor(dataset1$muts,levels=c(0,1),labels=c("Wildtype","mutant"))
head(dataset1)
table(dataset1$mutsF)
library(ggplot2)
p<-ggplot(dataset1,aes(x=mutsF,y=rank,fill=mutsF))+geom_boxplot(alpha=0.3)+
  facet_wrap(~genename)

p2<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
             axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
             legend.title = element_text(size = 8),
             legend.text = element_text(size = 8))

p2+ggtitle("Boxplot of ranks in Pancancer for a few significant genes")


#################
library(ggplot2)
  
ng<-length(mygenes)
p3<-alist()
for (i in seq(1,ng)){ 
  cancerT<-dataset1[which(dataset1$genename==mygenes[i]),-1]
  cancerT$mutsF<-factor(cancerT$muts,levels=c(0,1),labels=c("Wildtype","mutant"))
    p<-ggplot(cancerT,aes(x=mutsF,y=rank),fill=mutsF)+geom_boxplot(alpha=0.3)+
      scale_color_manual(values=c("#999999", "#E69F00"))
    
  p2<-p + theme(plot.title = element_text(size=8),axis.text.x=element_text(size=8),
                axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 8))
    
  
 p3[[i]]<- p2+ggtitle("Boxplot of ranks in Pancancer for a few significant genes")
  
  
} 


library(gridExtra)
somePDFPath = "F:\\THESIS-CODE\\myboxplot.pdf"
pdf(file=somePDFPath,onefile = TRUE)

for (i in seq(length(p3))) {
  do.call("grid.arrange", p3[[i]])  
}
dev.off()

################################################
dataset1<-as.data.frame(fread("../Processed_Data/afewgenesdataf.csv"))
######waterplot
library(tidyverse)
library(dplyr)
library(knitr)
mygenes<-sort(unique(dataset1$genename))
i=1
cancerp<-dataset1[which(dataset1$genename %in% mygenes[i]),]

aggregate(cancerp$rank,by=list(cancerp$muts),mean)
wwtp53<-t.test(cancerp$rank[cancerp$muts==0],cancerp$rank[cancerp$muts==1])
format(wwtp53$p.value,scientific=TRUE)

cancer2<-aggregate(cancerp$rank,by=list(cancerp$CLEANNAME),mean)
cancer2<-data.frame(cancer2)
cancer2$mutation<-0
colnames(cancer2)<-c("Cell lines","Meanrank","mutation")
wildtype<-unique(cancerp$CLEANNAME[which(cancerp$muts==0 & cancerp$genename==mygenes[i])])
mutant<-unique(cancerp$CLEANNAME[which(cancerp$muts==1 & cancerp$genename==mygenes[i])])

cancer2$mutation[which(cancer2$`Cell lines` %in% wildtype)]<-0
cancer2$mutation[which(cancer2$`Cell lines` %in% mutant)]<-1
cancer2$mutationn<-factor(cancer2$mutation,levels=c(0,1),labels=c("Wildtype","Mutant"))
aggregate(cancer2$Meanrank,by=list(cancer2$mutationn),median)
mm<-round(median(cancer2$Meanrank))
cancer2$rank<-cancer2$Meanrank-mm
cancerf<-cancer2[order(cancer2$rank),]
##
library(tidyverse)
library(dplyr)
library(knitr)   
##
col <- ifelse(cancerf$mutationn == "Wildtype", 
              "steelblue", # if dose = 80 mg, then the color will be steel blue
              "red") # if dose != 80 mg (i.e. 150 mg here), then the color will be cadet blue

p1<-barplot(cancerf$rank, 
              col=col, 
              border=col, 
              space=0.5, 
              ##ylim=c(-50,50),
              main = "Waterfall plot for BRAF", 
              ylab="Ranked free-batch-effect logFC",
              cex.axis=1.5, 
              legend.text= c( "Wildtype", "Mutant"),
              args.legend=list(title="Ranked logFC", fill=c("steelblue", "red"), border=NA, cex=0.1))


####for cancer specific we needs to select significant genes in 
##each of cancers from final data to one data and joint to annotation
###to separate cancer type as a variable and like above plot for each cancer plot for all selected 
#########significant genes and like bellow save in a pdf


###############################example1


somePDFPath = "F:\\THESIS-CODE\\some.pdf"
pdf(file=somePDFPath)  

for (i in seq(5,10))   
{   
  par(mfrow = c(2,1))
  VAR1=rnorm(i)  
  VAR2=rnorm(i)  
  plot(VAR1,VAR2)   
} 
dev.off() 


##########################SOME EXAMPLE
genes<-c("KRAS","TP53","PIK3CA")
generate.PDF <- function(data) {    
  pdf("F:/THESIS-CODE/myplot.pdf", width=8.5, height=5,onefile=T)
  plot1 <- plot(x,y)
  plot2 <- plot(y,x)
  plot3 <- plot(x,y*2)
  print(plot1)
  print(plot2)
  print(plot3)
  dev.off()
}

generate.PDF(data)

library(grid)
install.packages("graphics")
r <- raster(ncol=3, nrow=3)
r[] <- 1:ncell(r)
as.raster(r)
# }
