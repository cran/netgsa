## ----setup, include=FALSE------------------------------------------------
rm(list=ls())
knitr::opts_chunk$set(echo = TRUE)

## ----package, eval=FALSE-------------------------------------------------
#  library(glassoFast)
#  library(graphite)
#  library(igraph)
#  library(msigdbr)
#  library(netgsa)
#  library(Rgraphviz)

## ----data----------------------------------------------------------------
data(breastcancer2012, package = "netgsa")
ls()

## ----rownames------------------------------------------------------------
head(rownames(x))

## ----ER status-----------------------------------------------------------
table(group)

## ----edgelist------------------------------------------------------------
head(edgelist)

## ----packages, include=FALSE---------------------------------------------
library(glassoFast)
library(glmnet)
library(graphite)
library(igraph)
library(msigdbr) 
library(netgsa)

## ----pathways------------------------------------------------------------
paths <- pathways('hsapiens','kegg')
paths[[1]]
head(nodes(paths[[1]]))

## ----preparePathways-----------------------------------------------------
pathwayList <- preparePathways('kegg')
head(pathwayList[[1]])

## ----MSigDB--------------------------------------------------------------
pathwayList <- preparePathways('MSigDB')
head(pathwayList[[1]])

## ----our pathways--------------------------------------------------------
pathways[1:2]

## ----csv-----------------------------------------------------------------
write.csv(edgelist,file='edgelist.txt',row.names = FALSE)
out <- prepareAdjacencyMatrix(x, group, pathways, FALSE, 'edgelist.txt', NULL)

## ----prepareAdjacencyMatrix----------------------------------------------
genenames <- unique(c(pathways[[24]], pathways[[52]]))
genenames <- intersect(genenames, rownames(x))
p <- length(genenames)
p
sx <- x[match(genenames, rownames(x)),]
sout <- prepareAdjacencyMatrix(sx, group, pathways, FALSE, 'edgelist.txt', NULL)

## ----B-------------------------------------------------------------------
# pathway indicator matrix
dim(sout$B)

## ----estimate matrices---------------------------------------------------
ncond <- length(unique(group))
Amat <- vector("list",ncond)
sx <- sx[match(colnames(sout$B), rownames(sx)),]

for (k in 1:ncond){
  data_c <- sx[,(group==k)]
  # select the tuning parameter
  fitBIC <- bic.netEst.undir(data_c,one=sout$Adj,
                             lambda=seq(1,10)*sqrt(log(p)/ncol(data_c)),eta=0.1)
  # refit the network
  fit <- netEst.undir(data_c,one=sout$Adj,
                      lambda=which.min(fitBIC$BIC)*sqrt(log(p)/ncol(data_c)),eta=0.1)
  Amat[[k]] <- fit$Adj
}

## ----netgsa 1------------------------------------------------------------
test1 <- NetGSA(Amat, sx, group, pathways = sout$B, lklMethod = 'REHE')
head(test1$results)

## ----netgsa 2------------------------------------------------------------
sout <- prepareAdjacencyMatrix(sx, group, pathways, FALSE, 'edgelist.txt', NULL, estimate_network=TRUE, lambda_c = 9, eta=0.1)
test2 <- NetGSA(sout$Amat, sx, group, pathways = sout$B, lklMethod = 'REHE')
head(test2$results)

## ----netgsa DAG----------------------------------------------------------
# e.g. the "Adrenergic signaling in cardiomyocytes" pathway from KEGG is a DAG.
print(is_dag(g))

genenames <- V(g)$name
p <- length(genenames)

# reorder the variables and get the adjacency matrix
reOrder <- topo_sort(g,"in")
Adj <- as.matrix(get.adjacency(g))
Adj <- Adj[reOrder,reOrder]

B <- matrix(rep(1,p),nrow=1)
rownames(B) <- "Adrenergic signaling in cardiomyocytes"
colnames(B) <- rownames(Adj)
gx <- x[match(rownames(Adj), rownames(x)),]

Amat <- vector("list", 2)
for (k in 1:2){
  data_c <- gx[,which(group==k)]
  Amat[[k]] <- netEst.dir(data_c, one = Adj)$Adj
}
test <- NetGSA(Amat, gx, group, pathways = B, lklMethod = 'REHE')


