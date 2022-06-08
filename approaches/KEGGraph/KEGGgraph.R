### R code from vignette source 'KEGGgraph.Rnw'

###################################################
### code chunk number 1: lib
###################################################
library(KEGGgraph)


###################################################
### code chunk number 2: remoteRetrieval (eval = FALSE)
###################################################
## tmp <- tempfile()
## pName <- "p53 signaling pathway"
## data(KEGGPATHNAME2ID)
## pId <- mget(pName, KEGGPATHNAME2ID)[[1]]
## retrieveKGML(pId, organism="cel", destfile=tmp, method="wget", quiet=TRUE)


###################################################
### code chunk number 3: library
###################################################
mapkKGML <- system.file("extdata/hsa04010.xml",
                     package="KEGGgraph")


###################################################
### code chunk number 4: parsemapk
###################################################
mapkG <- parseKGML2Graph(mapkKGML,expandGenes=TRUE)
mapkG


###################################################
### code chunk number 5: parsemapk2
###################################################
mapkpathway <- parseKGML(mapkKGML)
mapkpathway
mapkG2 <- KEGGpathway2Graph(mapkpathway, expandGenes=TRUE)
mapkG2


###################################################
### code chunk number 6: nodeandedge
###################################################
mapkNodes <- nodes(mapkG)
nodes(mapkG)[1:3]
mapkEdges <- edges(mapkG)
edges(mapkG)[1]


###################################################
### code chunk number 7: keggnodedata
###################################################
mapkGnodedata <- getKEGGnodeData(mapkG)
mapkGnodedata[[2]]


###################################################
### code chunk number 8: keggnodedataalt (eval = FALSE)
###################################################
## getKEGGnodeData(mapkG, 'hsa:5924')


###################################################
### code chunk number 9: keggedgedata
###################################################
mapkGedgedata <- getKEGGedgeData(mapkG)
mapkGedgedata[[4]]


###################################################
### code chunk number 10: keggedgedataalt (eval = FALSE)
###################################################
## getKEGGedgeData(mapkG,'hsa:627~hsa:4915')


###################################################
### code chunk number 11: inout
###################################################
mapkGoutdegrees <- sapply(edges(mapkG), length)
mapkGindegrees <- sapply(inEdges(mapkG), length)
topouts <- sort(mapkGoutdegrees, decreasing=T)
topins <- sort(mapkGindegrees, decreasing=T)
topouts[1:3]
topins[1:3]


###################################################
### code chunk number 12: subsetprepare
###################################################
library(Rgraphviz)
set.seed(123)
randomNodes <- sample(nodes(mapkG), 25)
mapkGsub <- subGraph(randomNodes, mapkG)
mapkGsub


###################################################
### code chunk number 13: makeattr
###################################################
makeAttr <- function(graph, default, valNodeList) {
  tmp <- nodes(graph)
  x <- rep(default, length(tmp)); names(x) <- tmp
  
  if(!missing(valNodeList)) {
    stopifnot(is.list(valNodeList))
    allnodes <- unlist(valNodeList)
    stopifnot(all(allnodes %in% tmp))
    for(i in seq(valNodeList)) {
      x[valNodeList[[i]]] <- names(valNodeList)[i]
    }
  }
  return(x)
}


###################################################
### code chunk number 14: subsetplot
###################################################
outs <- sapply(edges(mapkGsub), length) > 0
ins <- sapply(inEdges(mapkGsub), length) > 0
ios <- outs | ins

## translate the KEGG IDs into Gene Symbol
if(require(org.Hs.eg.db)) {
  ioGeneID <- translateKEGGID2GeneID(names(ios))
  nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
} else {
  nodesNames <- names(ios)
}
names(nodesNames) <- names(ios)

nAttrs <- list();
nAttrs$fillcolor <- makeAttr(mapkGsub, "lightgrey", list(orange=names(ios)[ios]))
nAttrs$label <- nodesNames
plot(mapkGsub, "neato", nodeAttrs=nAttrs,
     attrs=list(node=list(fillcolor="lightgreen",
                  width="0.75", shape="ellipse"), 
       edge=list(arrowsize="0.7")))


###################################################
### code chunk number 15: mergedemo
###################################################
wntKGML <- system.file("extdata/hsa04310.xml",package="KEGGgraph")
wntG <- parseKGML2Graph(wntKGML)
graphs <- list(mapk=mapkG, wnt=wntG)
merged <- mergeGraphs(graphs)
merged


###################################################
### code chunk number 16: bcc (eval = FALSE)
###################################################
## library(RBGL)
## bcc <- brandes.betweenness.centrality(mapkG)
## rbccs <- bcc$relative.betweenness.centrality.vertices[1L,]
## toprbccs <- sort(rbccs,decreasing=TRUE)[1:4]
## toprbccs


###################################################
### code chunk number 17: bccplot (eval = FALSE)
###################################################
## toprbccName <- names(toprbccs)
## toprin <- sapply(toprbccName, function(x) inEdges(mapkG)[x])
## toprout <- sapply(toprbccName, function(x) edges(mapkG)[x])
## toprSubnodes <- unique(unname(c(unlist(toprin), unlist(toprout), toprbccName)))
## toprSub <- subGraph(toprSubnodes, mapkG)
## 
## nAttrs <- list()
## tops <- c("MAPK3K1","GRB2","MAP2K2","MAP2K1")
## topLabels <- lapply(toprbccName, function(x) x); names(topLabels) <- tops
## nAttrs$label <- makeAttr(toprSub, "", topLabels)
## nAttrs$fillcolor <- makeAttr(toprSub, "lightblue", list(orange=toprbccName))
## nAttrs$width <- makeAttr(toprSub,"",list("0.8"=toprbccName))
## 
## plot(toprSub, "twopi", nodeAttrs=nAttrs, attrs=list(graph=list(start=2)))


###################################################
### code chunk number 18: help (eval = FALSE)
###################################################
## help(package=KEGGgraph)


###################################################
### code chunk number 19: reactionexample
###################################################
mapfile <-  system.file("extdata/map00260.xml",package="KEGGgraph")
map <- parseKGML(mapfile)
map
reactions <- getReactions(map)
reactions[[1]]


###################################################
### code chunk number 20: cnexample
###################################################
chemicalGraph <- KEGGpathway2reactionGraph(map)

outDegrees <- sapply(edges(chemicalGraph), length)
maxout <- names(sort(outDegrees,decreasing=TRUE))[1:3]

nAttrs <- list()
maxoutlabel <- as.list(maxout); names(maxoutlabel) <- maxout
nAttrs$label <- makeAttr(chemicalGraph, "", maxoutlabel)
nAttrs$fillcolor <- makeAttr(chemicalGraph, "lightblue", list(orange=maxout))
nAttrs$width <- makeAttr(chemicalGraph,"0.8", list("1.8"=maxout))
plot(chemicalGraph, nodeAttrs=nAttrs)


###################################################
### code chunk number 21: mapk14expand (eval = FALSE)
###################################################
## mapkGembed <- parseKGMLexpandMaps(mapkKGML)


###################################################
### code chunk number 22: subgraphbynode
###################################################
mapkGall <- parseKGML2Graph(mapkKGML,genesOnly=FALSE)
mapkGall
mapkGsub <- subGraphByNodeType(mapkGall, "gene")
mapkGsub


###################################################
### code chunk number 23: biomart (eval = FALSE)
###################################################
## toprbccKEGGID <- names(toprbccs)
## toprbccKEGGID
## toprbccGeneID <- translateKEGGID2GeneID(toprbccKEGGID)
## toprbccGeneID


###################################################
### code chunk number 24: orgHuman (eval = FALSE)
###################################################
## if(require(org.Hs.eg.db)) {
##   tnodes <- nodes(toprSub)
##   tgeneids <- translateKEGGID2GeneID(tnodes)
##   tgenesymbols <- sapply(mget(tgeneids, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
##   toprSubSymbol <- toprSub
##   nodes(toprSubSymbol) <- tgenesymbols
##   plot(toprSubSymbol, "neato",attrs=list(node=list(font=5, fillcolor="lightblue")))
## }


###################################################
### code chunk number 25: biomart2 (eval = FALSE)
###################################################
## library(biomaRt)
## hsapiens <- useMart("ensembl","hsapiens_gene_ensembl" )
## filters <- listFilters(hsapiens)
## getBM(attributes=c("entrezgene","hgnc_symbol"), 
##       filters="entrezgene", 
##       values=toprbccGeneID, mart=hsapiens)


###################################################
### code chunk number 26: sessionInfo
###################################################
sessionInfo()


