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
#######################END OF MUTATIONAL MATRIX
###############THIS IS DENOISE RANKED DATA AFTER SVA PACKAGE 
pvaladj<-function(results){
  n<-length(unique(results$PvalueNbModel))
  p<-results$PvalueNbModel
  results$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
  
  colnames(results)<-c("gene","NWildType","NMutant","PvalueNbModel","adjPvalue")
  results[which(results$adjPvalue<0.05),] 
}
breast<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvbreast.csv"))
breastf<-pvaladj(breast)
unique(breastf$gene)
write.csv(breastf,"../OutPut/BreastCancerAdjPvalue.csv")
##biliary_tract<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvbiliary_tract.csv"))
##biliary_tractf<-pvaladj(biliary_tract)

bladder<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvbladder.csv"))
bladderf<-pvaladj(bladder)

bone<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvbone.csv"))
bonef<-pvaladj(bone)

#cervix<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvcervix.csv"))
#cervixf<-pvaladj(cervix)

cns<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvcns.csv"))
cnsf<-pvaladj(cns)

colorectal<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvcolorectal.csv"))
colorectalf<-pvaladj(colorectal)

endometrium<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvendometrium.csv"))
endometriumf<-pvaladj(endometrium)
endometriumf$gene

#eye<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csveye.csv"))
#eyef<-pvaladj(eye)

gastric<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvgastric.csv"))
gastricf<-pvaladj(gastric)

kidney<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvkidney.csv"))
kidneyf<-pvaladj(kidney)

leukemia<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvleukemia.csv"))
leukemiaf<-pvaladj(leukemia)

liver<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvliver.csv"))
liverf<-pvaladj(liver)

lung<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvlung.csv"))
lungf<-pvaladj(lung)

lymphoma<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvlymphoma.csv"))
lymphomaf<-pvaladj(lymphoma)

oesophagus<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvoesophagus.csv"))
oesophagusf<-pvaladj(oesophagus)

ovary<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvovary.csv"))
ovaryf<-pvaladj(ovary)

pancreas<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvpancreas.csv"))
pancreasf<-pvaladj(pancreas)

##pnet<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvpnet.csv"))
##pnetf<-pvaladj(pnet)

#prostate<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvprostate.csv"))
#prostatef<-pvaladj(prostate)

#salivary_gland<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvsalivary_gland.csv"))
#salivary_glandf<-pvaladj(salivary_gland)

skin<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvskin.csv"))
skinf<-pvaladj(skin)

soft_tissue<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvsoft_tissue.csv"))
soft_tissuef<-pvaladj(soft_tissue)

##thyroid<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvthyroid.csv"))
##thyroidf<-pvaladj(thyroid)

##undefined<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvundefined.csv"))
##undefinedf<-pvaladj(undefined)

upper_aerodigestive_tract<-as.data.frame(fread("../OutPut/EachCancerResults/eachcancerResults.csvupper_aerodigestive_tract.csv"))
upper_aerodigestive_tractf<-pvaladj(upper_aerodigestive_tract)

topgenes<-c("KRAS","NRAS","TP53","PIK3CA","CTNNB1")
cancers<-list(breastf,bladderf,bonef,cnsf,colorectalf,endometriumf,gastricf,kidneyf,leukemiaf,
     liverf,lungf,lymphomaf,oesophagusf,ovaryf,pancreasf,skinf,soft_tissuef,upper_aerodigestive_tractf)
cancersN<-c("Breast","Bladder","Bone","Cns","Colorectalf","Endometrium","Gastric",
               "Kidney","Leukemia",
              "Liver","Lung","Lymphoma","Oesophagus","Ovary","Pancreas","Skin",
              "Soft_tissue","Upper_aerodigestive_tractf")

cancertopgenes<-matrix(0,nrow=length(cancersN),ncol=length(topgenes))
rownames(cancertopgenes)<-cancersN
colnames(cancertopgenes)<-topgenes
#names(cancertopgenes)<-cancersN
  for(i in 1:length(cancersN)){
    for(j in 1:length(topgenes)){
      index<-which(cancers[[i]]$gene %in% topgenes[j])
      if(length(index)>0){
      cancertopgenes[i,which(colnames(cancertopgenes) %in% topgenes[j])]<-cancers[[i]][index,]$adjPvalue
    }
    }
  }
myfilefinal<-data.frame(cancertopgenes)
myfilefinal$cancers<-rownames(cancertopgenes)
library("gplots")
heatmap.2(cancertopgenes, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
 library("pheatmap")
  pheatmap(cancertopgenes, cutree_rows = 4)

  ##############common genes in cancers
  topgenes2<-unique(c(breastf$gene,bladderf$gene,bonef$gene,cnsf$gene,colorectalf$gene
               ,endometriumf$gene,gastricf$gene,kidneyf$gene,leukemiaf$gene,
                liverf$gene,lungf$gene,lymphomaf$gene,oesophagusf$gene,ovaryf$gene
               ,pancreasf$gene,skinf$gene,soft_tissuef$gene,upper_aerodigestive_tractf$gene))
  
  #names(cancertopgenes)<-cancersN
  topgenes22<-c(breastf$gene,bladderf$gene,bonef$gene,cnsf$gene,colorectalf$gene
                      ,endometriumf$gene,gastricf$gene,kidneyf$gene,leukemiaf$gene,
                      liverf$gene,lungf$gene,lymphomaf$gene,oesophagusf$gene,ovaryf$gene
                      ,pancreasf$gene,skinf$gene,soft_tissuef$gene,upper_aerodigestive_tractf$gene)
 ttg<- data.frame(table(topgenes22))
 genest<-ttg$topgenes22[which(ttg$Freq>1)]##these genes are more than one cancer
 #genest[1:12]
 # ttg$Freq[(which(ttg$topgenes2 %in% "NRAS"))]
  cancertopgenes<-matrix(0,nrow=length(cancersN),ncol=length(genest))
  rownames(cancertopgenes)<-cancersN
  colnames(cancertopgenes)<-genest
  for(i in 1:length(cancersN)){
    for(j in 1:length(genest)){
      index<-which(cancers[[i]]$gene %in% genest[j])
      if(length(index)>0 &   ttg$Freq[(which(ttg$topgenes22 %in% genest[j]))]){
        cancertopgenes[i,which(colnames(cancertopgenes) %in% genest[j])]<-cancers[[i]][index,]$adjPvalue
      }
    }
  }
  library("pheatmap")
  cancertopgeneslog<-log(cancertopgenes)
  cancertopgeneslog[which(!is.finite(cancertopgeneslog))] <- 0
  pheatmap(cancertopgeneslog, cutree_rows = 4,main="Heatmap of log(P-value) of Top genes in cancers")
  
  ##library(ComplexHeatmap)
  
 ## Heatmap(cancertopgeneslog, name="TopGenes",row_names_gp = gpar(fontsize = 7))
  
  #####################Fisher Score of p-values
  cancers<-list(breast,bladder,bone,cns,colorectal,endometrium,gastric,kidney,leukemia,
                liver,lung,lymphoma,oesophagus,ovary,pancreas,skin,soft_tissue,upper_aerodigestive_tract)

  pvaladj2<-function(results){
    n<-length(unique(results$PvalueNbModel))
    p<-results$PvalueNbModel
    results$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
    
    colnames(results)<-c("gene","NWildType","NMutant","PvalueNbModel","adjPvalue")
     return(results)
  }
  mygenes<-c("KRAS","NRAS","PIK3CA","TP53","BRAF","HRAS","APC","AXIN1","CTNNB1",
             "TOP2A","RRAS2","CDK4","JHDM1D","PA2G4","RPS6KA4","FLT3")
    FI<-matrix(0,nrow=length(cancers),ncol=length(mygenes))
  for(i in 1:length(cancers)){
    cancers[[i]]<-pvaladj2(cancers[[i]])
    for(j in 1:length(mygenes)){
      if(length(which(cancers[[i]]$gene %in% mygenes[j]))>0){
        
      
     FI[i,j] <-2*(-log(cancers[[i]]$adjPvalue [which(cancers[[i]]$gene %in% mygenes[j])]))
      
    }
    }
  }
  
   FIscores<-colSums(FI)
   names(FIscores)<-mygenes
   FIscores<- FIscores[order( -FIscores)]
   FIscores2<-data.frame(FIscores)
   FIscores2$gene<- names(FIscores)
   
   colnames(FIscores2)<-c("FisherScore","Gene")
   FIscores2<-FIscores2[order(-FIscores2$FisherScore),]
   FIscores2$genes<-factor(FIscores2$Gene)
    library(ggplot2)
   ggplot(FIscores2,aes(x=Gene,y=FisherScore,order=-FisherScore))+geom_point()

   ggplot(FIscores2) +
     geom_line(aes(x = Gene, y = FisherScore))  
   library(dplyr)
   ggplot(FIscores2 %>%
            arrange(FisherScore),
          aes(REORDER(FisherScore),x = genes, y = FisherScore)) +
     geom_point()
   
   bb<-barplot(FIscores2$FisherScore,names.arg=FIscores2$genes
   ,main="Plot of Fisher Scores of Top Genes in All Cancers")###this the plot of ineterst
  ## text(bb,FIscores2$FisherScore,round(FIscores2$FisherScore,0),cex=1)
   
   p<-ggplot(data=FIscores2, aes(x=reorder(genes,-FisherScore), y=FisherScore)) +
     geom_bar(stat="identity", fill="steelblue")+
     theme_minimal()
   p+labs(x="Genes")+ggtitle("Barplot of Fisher Scores of Top Genes in all Cancers")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
   
   ##########
   library(ggplot2)
   library(dplyr)
   
   # The dataset is provided in the gapminder library
   library(gapminder)
   data <- gapminder %>% filter(year=="2007") %>% dplyr::select(-year)
   
   # Most basic bubble plot
   FIscores2 %>%
     arrange(desc(FisherScore)) %>%
     mutate(gene = factor(Gene, Gene)) %>%
     ggplot(aes(x=Gene, y=FisherScore, size = FisherScore)) +
     geom_point(alpha=0.5) +
     scale_size(range = c(.1, 24), name="Population (M)")
   ##########Bubble plot
   siggenesP<-as.data.frame(fread("../OutPut/PancancerResults.csv"))
   n<-length(unique(siggenesP$PvalueNbModel))
   p<-siggenesP$PvalueNbModel
   siggenesP$adjPvalue<-p.adjust(p, method = "fdr", n = length(p))
   mygenes<-siggenesP$X[which(siggenesP$adjPvalue<0.05)]
   
   
   cancertopgenes<-matrix(0,nrow=length(cancersN),ncol=length(mygenes))
   rownames(cancertopgenes)<-cancersN
   colnames(cancertopgenes)<-mygenes
   #names(cancertopgenes)<-cancersN
   for(i in 1:length(cancersN)){
     for(j in 1:length(mygenes)){
       index<-which(cancers[[i]]$gene %in% mygenes[j])
       if(length(index)>0){
         cancertopgenes[i,which(colnames(cancertopgenes) %in% mygenes[j])]<-(-2*log(cancers[[i]][index,]$adjPvalue))
       }
     }
   }
   myfilefinal<-data.frame(cancertopgenes)
   myfilefinal$cancers<-rownames(cancertopgenes)
   library("gplots")
   
   
   #########
   heatmap.2(cancertopgenes, scale = "none", col = bluered(100), 
             trace = "none", density.info = "none")
   library("pheatmap")
   pheatmap(cancertopgenes, cutree_rows = 4,fontsize_row = 8, fontsize_col = 7,
            main = "Heatmap of -2log(P-values) of top genes in Pancancer in cancer types")
   
   