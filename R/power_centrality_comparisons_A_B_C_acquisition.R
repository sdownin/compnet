setwd("C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2")
img_dir <-"C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2/img"
#.libPaths('C:/Users/sdowning/Documents/R/win-library/3.2')
library(igraph)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

##=============================================================
##              Helper Functions
##-------------------------------------------------------------
## wrapper for igraph power centrality function
power <- function(graph, beta){
  if (class(graph) != 'igraph') 
    stop("graph must be an igraph object.")
  return(igraph::power_centrality(graph = graph, 
                                  nodes = V(graph), 
                                  exponent = beta, 
                                  rescale = F,
                                  loops = F))
}

## map a vector of power centralities to values between [min,max] node sizes for plotting
map <- function(x, min=5, max=23) {
  xmin <- x[which.min(x)]
  xmax <- x[which.max(x)]
  maxdiff <- abs(xmax - xmin)
  scalediff <- max - min
  return(sapply(x, function(i) (scalediff*(i - xmin)/maxdiff) + min ))
}

## plot graph
plotGraph <- function(graph, ...) {
  .par <- par()
  par(mar=c(.1,.1,.1,.1))
  igraph::plot.igraph(
    x=graph,
    vertex.label.cex=.7,
    vertex.size=9,
    edge.wdith=1.5,
    ...
  )
  par(.par)
}

##=============================================================
##              Make Competition Network
##-------------------------------------------------------------

## Make graph by edges
g1 <- igraph::graph.data.frame(data.frame(
  from=c('B','B','B','C','C','C','4','4','4', '8','8','8',   '13','13','13', 'A','A','A'), 
  to  =c('1',   '2',   '3',   '4',   '13',  '8',   '5','6','7', '9','10','11', '14','15','12', '16',   '17',   '18')      
), directed=F)

gbc <- igraph::graph.data.frame(data.frame(
  from=c('B', 'B','B','B','C','C','C','4','4','4', '8','8','8',   '13','13','13', 'A','A','A'), 
  to  =c('C', '1',   '2',   '3',   '4',   '13',  '8',   '5','6','7', '9','10','11', '14','15','12', '16',   '17',   '18')      
), directed=F)

gba <- igraph::graph.data.frame(data.frame(
  from=c('B', 'B','B','B','C','C','C','4','4','4', '8','8','8',   '13','13','13', 'A','A','A'), 
  to  =c('A', '1',   '2',   '3',   '4',   '13',  '8',   '5','6','7', '9','10','11', '14','15','12', '16',   '17',   '18')      
), directed=F)


vertex.attr.comb <- list(name=function(x)paste(x, collapse="|"))

mapping1 <- as.integer(V(g1))
mapping1[V(g1)$name=='C'] <- which(V(g1)$name == 'B')
gb.c <- igraph::contract(g1, mapping = mapping1, vertex.attr.comb = vertex.attr.comb)

mapping2 <- as.integer(V(g1))
mapping2[V(g1)$name=='A'] <- which(V(g1)$name == 'B')
gb.a <- igraph::contract(g1, mapping = mapping2, vertex.attr.comb = vertex.attr.comb)


##=============================================================
##            Power Centrality
##-------------------------------------------------------------

betas <- c('-0.4'=-0.4,'-0.3'=-0.3,'-0.2'=-0.2,'-0.1'=-0.1,
           '0.1'=0.1,'0.2'=0.2,'0.3'=0.3,'0.4'=0.4)

df.pow <- as.data.frame(sapply(betas, function(x)power(g1, x)))
print(df.pow[which(rownames(df.pow) %in% c('B','C','A')), c('-0.1','0.1')])

df.pow <- as.data.frame(sapply(betas, function(x)power(gb.c, x)))
print(df.pow[which(rownames(df.pow) %in% c('B|C','C','A')), c('-0.1','0.1')])

df.pow <- as.data.frame(sapply(betas, function(x)power(gb.a, x)))
print(df.pow[which(rownames(df.pow) %in% c('B|A','C','A')), c('-0.1','0.1')])



##=============================================================
##            Power Centrality
##-------------------------------------------------------------
## compute power centralities for range of beta values
# betas <- c('-3'=-3,'-2'=-2,'-1.5'=-1.5,'-0.5'=-0.5,'0.5'=0.5,'1.5'=1.5,'2'=2,'3'=3)
# betas <- c('-0.2'=-0.2,'0.2'=0.2)
betas <- c('-0.4'=-0.4,'-0.3'=-0.3,'-0.2'=-0.2,'-0.1'=-0.1,
           '0.1'=0.1,'0.2'=0.2,'0.3'=0.3,'0.4'=0.4)

df.pow <- as.data.frame(sapply(betas, function(x)power(g1, x)))
print(df.pow[which(rownames(df.pow) %in% c('B','C','A')), c('-0.1','0.1')])

df.pow.bf <- as.data.frame(sapply(betas, function(x)power(gbf, x)))
print(df.pow.bf[which(rownames(df.pow.bf) %in% c('B','C','A')), c('-0.1','0.1')])

df.pow.bl <- as.data.frame(sapply(betas, function(x)power(gbl, x)))
print(df.pow.bl[which(rownames(df.pow.bl) %in% c('B','C','A')), c('-0.1','0.1')])


##=============================================================
##            Visualize
##-------------------------------------------------------------
gx <- gbl
df.pow <- as.data.frame(sapply(betas, function(x)power(gx, x)))
print(df.pow[which(rownames(df.pow) %in% c('B','C','A')), c('-0.1','0.1')])
par(mfrow=c(1,1))
m <- ncol(df.pow)
n <- nrow(df.pow)
SEED <- 333
font <- 17
legcex <- 1
nodeNames <- c('B','C','A')
layoutFunction <- layout.kamada.kawai
#
beta <- '-0.4'
plot.igraph(gx, layout=layoutFunction,
            vertex.size = map(df.pow[,beta]),
            vertex.label.cex = map(df.pow[,beta],min = .8,max = 1.6),
            main=sprintf('Power Centrality (beta = %s)',beta),
            vertex.color=sapply(rownames(df.pow),function(name)ifelse(name %in% nodeNames, 'pink', 'gray')),
            vertex.label.color=sapply(rownames(df.pow),function(name)ifelse(name %in% nodeNames, 'darkred','black')),
            vertex.label.font=font)


### comparison
par(mfrow=c(1,2), mar=c(0,0,3,0))
nodes <- c('C',as.character(4:15))
gx2 <- igraph::induced.subgraph(gx, V(gx)[V(gx)$name %in% nodes])
df.pow2 <- df.pow[row.names(df.pow) %in% nodes,]
beta <- '0.4'
plot.igraph(gx2, layout=layoutFunction,
            vertex.size = map(df.pow2[,beta]),
            vertex.label.cex = map(df.pow2[,beta],min = .8,max = 1.6),
            main=sprintf('Cooperative Power Centrality\n(beta = %s)',beta),
            vertex.color=sapply(rownames(df.pow2),function(name)ifelse(name %in% nodeNames, 'pink', 'gray')),
            vertex.label.color=sapply(rownames(df.pow2),function(name)ifelse(name %in% nodeNames, 'darkred','black')),
            vertex.label.font=font)
beta <- '-0.4'
plot.igraph(gx2, layout=layoutFunction,
            vertex.size = map(df.pow2[,beta]),
            vertex.label.cex = map(df.pow2[,beta],min = .8,max = 1.6),
            main=sprintf('Adversarial Power Centrality\n(beta = %s)',beta),
            vertex.color=sapply(rownames(df.pow2),function(name)ifelse(name %in% nodeNames, 'pink', 'gray')),
            vertex.label.color=sapply(rownames(df.pow2),function(name)ifelse(name %in% nodeNames, 'darkred','black')),
            vertex.label.font=font)
