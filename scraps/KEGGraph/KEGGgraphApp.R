### R code from vignette source 'KEGGgraphApp.Rnw'

###################################################
### code chunk number 1: loadtoy
###################################################
library(KEGGgraph)
toyKGML <- system.file("extdata/kgml-ed-toy.xml", package="KEGGgraph")
toyGraph <- parseKGML2Graph(toyKGML, genesOnly=FALSE)
toyGraph
nodes(toyGraph)


###################################################
### code chunk number 2: vistoy
###################################################
library(Rgraphviz)
nodeInfo <- getKEGGnodeData(toyGraph)
nodeType <- sapply(nodeInfo, getType)
makeNodeRenderAttrs <- function(g, label=nodes(g),
                                shape="ellipse", fill="#e0e0e0",...) {
  rv <- list(label=label, shape=shape, fill=fill, ...)
  nA <- nodeRenderInfo(g)
  for(i in seq(along=rv)) {
    if (length(rv[[i]]) == 1) {
      rv[[i]] <- rep(rv[[i]], numNodes(g))
    } else {
      if (length(rv[[i]]) != numNodes(g))
        stop("Attribute vector must have as many elements as 'g' has nodes.")
    }
    names(rv[[i]]) <- nodes(g)
    nA[[ names(rv)[[i]] ]] <- rv[[i]]
  }
  nodeRenderInfo(g) <- nA
  return(g)
}
toyDrawn <- plotKEGGgraph(toyGraph)
toyDrawnRefine <- makeNodeRenderAttrs(toyDrawn, fill=ifelse(nodeType=="gene", "lightblue", "orange"),
                                      shape=ifelse(nodeType=="gene", "ellipse","rectangle"))


###################################################
### code chunk number 3: vistoyReal
###################################################
renderGraph(toyDrawnRefine)


###################################################
### code chunk number 4: spiaLoad
###################################################
if(require(SPIA)) {
  data(colorectalcancer,package="SPIA")
} else {
  load(system.file("extdata/colorectalcancerSPIA.RData", package="KEGGgraph"))
}
head(top)


###################################################
### code chunk number 5: spiaTrans
###################################################
library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
top$ENTREZ <- unlist(as.list(x[top$ID]))
top <- top[!is.na(top$ENTREZ),]
top <- top[!duplicated(top$ENTREZ),]
tg1 <- top[top$adj.P.Val < 0.05,]
DE_Colorectal <- tg1$logFC
names(DE_Colorectal) <- as.vector(tg1$ENTREZ)
ALL_Colorectal <- top$ENTREZ


###################################################
### code chunk number 6: remoteRetreval (eval = FALSE)
###################################################
## tmp <- "hsa05210.xml"
## retrieveKGML(pathwayid="05210", organism="hsa", destfile=tmp)


###################################################
### code chunk number 7: spiaPer
###################################################
colFile <- system.file("extdata/hsa05210.xml", 
                       package="KEGGgraph")
g <- parseKGML2Graph(colFile)
deKID <- translateGeneID2KEGGID(names(DE_Colorectal))
allKID <- translateGeneID2KEGGID(ALL_Colorectal)
isDiffExp <- nodes(g) %in% deKID
sprintf("%2.2f%% genes differentially-expressed", mean(isDiffExp)*100)


###################################################
### code chunk number 8: spiaVis
###################################################
library(RColorBrewer)
library(org.Hs.eg.db)
library(RBGL)
library(grid)
ar <- 20
cols <- rev(colorRampPalette(brewer.pal(6, "RdBu"))(ar))
logfcs <- DE_Colorectal[match(nodes(g), deKID)]
names(logfcs) <- nodes(g)
logfcs[is.na(logfcs)] <- 0
incol <- round((logfcs+2)*5); incol[incol>ar] <- ar
undetected <- !nodes(g) %in% allKID
logcol <- cols[incol]; logcol[logfcs==0] <- "darkgrey"; logcol[undetected] <- "yellow"
names(logcol) <- names(logfcs)
nA <- makeNodeAttrs(g, fillcolor=logcol, label="", width=10, height=1.2)
par(mar=c(3,5,0,5), mgp=c(0,0,0))
layout(mat=matrix(c(rep(1,8),2), ncol=1, byrow=TRUE))
plot(g, "dot", nodeAttrs=nA)
image(as.matrix(seq(1,ar)), col=cols, yaxt="n", xaxt="n")
mtext("down-regulation", side=1,  at=0, line=1)
mtext("up-regulation", side=1,  at=1, line=1)


###################################################
### code chunk number 9: spiaSingleNode
###################################################
gDeg <- degree(g)
gIsSingle <- gDeg[[1]] + gDeg[[2]] == 0
options(digits=3)
gGeneID <- translateKEGGID2GeneID(nodes(g))
gSymbol <-  sapply(gGeneID, function(x) mget(x, org.Hs.egSYMBOL, ifnotfound=NA)[[1]])
isUp <- logfcs > 0
isDown <- logfcs < 0
singleUp <- isUp & gIsSingle
singleDown <- isDown & gIsSingle


###################################################
### code chunk number 10: spiaSub
###################################################
ups <- nodes(g)[logfcs > 0]
upNs <- unique(unlist(neighborhood(g, ups, return.self=TRUE)))
upSub <- subKEGGgraph(upNs, g)
upNeighbor <- nodes(upSub)[sapply(neighborhood(upSub, nodes(upSub)), length)>0]
upNeighbor <- setdiff(upNeighbor, nodes(g)[undetected])
upSub <- subKEGGgraph(upNeighbor, upSub)
upSubGID <- translateKEGGID2GeneID(nodes(upSub))
upSymbol <- gSymbol[upSubGID]
upnA <- makeNodeAttrs(upSub, fillcolor=logcol[nodes(upSub)], label=upSymbol, fixedsize=TRUE, width=10, height=10, font=20)
plot(upSub, "dot", nodeAttrs=upnA)


