

library(gprofiler2)
genes<-c("PAX5","RASGRP2","USP19","ROS1","HIPK3","NPFFR2","KRAS","BAZ2B"
,"LYST","TINF2","NAV2")
gostres <- gost(query = genes, organism = "hsapiens")

# The result is a named list where "result" is a data.frame with the enrichment analysis results
# and "meta" containing a named list with all the metadata for the query.
head(gostres$result)
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p
