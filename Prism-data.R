

################33prism data
sampleinfo = read.csv("E:\\Git\\genesilencing_project\\org_data\\CCLE\\sample_info.csv")
prism = read.csv("E:/Git/genesilencing_project/org_data/PRISM/prismdata.csv", stringsAsFactors = F)
prism<-data.frame(prism)
sampleinfo<-data.frame(sampleinfo)
head(sampleinfo)
dim(prism)
prism[1:4,1:4]

colnames(prism)[colnames(prism) == "X"] <- "DepMap_ID"
colnames(prism)[1]
colnames(sampleinfo)[1]
intersect(unique(sampleinfo$DepMap_ID),unique(prism$DepMap_ID))
cellinesprism<-sampleinfo[,1:4]
prismfinal<-merge(prism,cellinesprism,by="DepMap_ID")
prismfinal[1:4,1:4]
dim(prismfinal)
unique(colnames(prismfinal))
wprism<-sampleinfo[which(tolower(sampleinfo$stripped_cell_line_name) %in% wildtype),"DepMap_ID"]
mprism<-sampleinfo[which(tolower(sampleinfo$stripped_cell_line_name) %in% mutant),"DepMap_ID"]
vd<-prism[which(prism$DepMap_ID %in% mprism),-1]
sum(vd[1:2])
t.test(prism[which(prism$DepMap_ID %in% wprism),-1],prism[which(prism$DepMap_ID %in% mprism),-1])
#