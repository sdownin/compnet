library(igraph)
library(Matrix)

## save plots
dirname <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\acquisitions"
setwd(dirname)

## cache default params
.par = par()

##
#
##
plot2 <- function(gx, layout=layout.fruchterman.reingold, vertex.size=15, focal.firm=NA, fam='sans', edge.curved=F, seed=11111, ...)
{
  if ('type' %in% igraph::list.vertex.attributes(gx)) {
    vcolors <- sapply(V(gx)$type, function(x)ifelse(x, "SkyBlue2", "gray"))
    lcolors <-  sapply(V(gx)$type, function(x)ifelse(x, "darkblue", "black"))
    vshapes <- sapply(1:vcount(gx),function(x)ifelse(V(gx)$type[x], "circle", "square"))
    isBipartite <- length(unique(V(gx)$type)) > 1
  } else {
    vcolors <- rep("SkyBlue2", vcount(gx))
    lcolors <-  rep("darkblue", vcount(gx))
    vshapes <- rep("circle", vcount(gx))
    isBipartite <- FALSE
  }
  vAttrs = names(igraph::vertex.attributes(gx))
  fonts <- rep(1, vcount(gx))
  framecols <- rep('black', vcount(gx))
  framewidths <- rep(1, vcount(gx)) 
  if(!is.na(focal.firm)) {
    vcolors[V(gx)$name==focal.firm] <- 'darkblue'
    lcolors[V(gx)$name==focal.firm] <- 'white'
  }
  if(!isBipartite) {
    adjmat <- as_adjacency_matrix(gx, attr = 'weight', sparse = F)
    ffidx <- which(V(gx)$name==focal.firm)
    mmcidx <- unname(which(adjmat[ , ffidx] > 1))
    framecols[mmcidx] <- 'darkred'
    lcolors[mmcidx] <- 'darkred'
    framewidths[mmcidx] <- 5
    fonts[mmcidx] <- 4
  }
  set.seed(seed)
  plot(gx, 
       layout = layout, 
       layout.par = list(), 
       labels = NULL, 
       label.color = lcolors, 
       label.font = NULL, 
       label.degree = -pi/4, 
       label.dist = 0, 
       vertex.label=sapply(1:vcount(gx), function(x) ifelse("name" %in% vAttrs, V(gx)$name[x], x)),
       vertex.color = vcolors, 
       vertex.shape = vshapes,
       vertex.size = vertex.size, 
       vertex.frame.color=framecols, 
       vertex.frame.width=framewidths, 
       vertex.label.family=fam,  # Font family of the label (e.g."Times", "Helvetica")
       vertex.label.font=fonts,  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
       vertex.label.color=lcolors,
       edge.color = "darkgrey", 
       edge.width = E(gx)$weight^2,
       edge.labels = NA, 
       edge.lty=1, 
       margin=0,
       loop.angle=0, 
       axes = FALSE, 
       xlab = "", 
       ylab = "",
       xlim=c(-1,1), 
       ylim=c(-1,1), 
       edge.curved=edge.curved,
       ...)
}


##
# Bipartite Graph Acquisition
##
biAcq <- function(gi, acquirer, target, decay=-0.2, project=T, verbose=T)
{
  if (project) {
    gi <- bipartite.projection(gi, remove.type = F)$proj2
    V(gi)$type <- unlist(V(gi)$type)
    V(gi)$name <- unlist(V(gi)$name)
  }
  
  tdf <- as_data_frame(gi, what='vertices')
  tdf$before <- power_centrality(gi, exponent = decay)
  tdf$after <- NA
  
  vnamemap <- names(V(gi))
  vmap <- as.integer(V(gi))
  
  revord <- which(vnamemap==target) < which(vnamemap==acquirer)
  comb.func <- ifelse(revord, 'last', 'first')
  
  vmap[vnamemap==target] <- vmap[vnamemap==acquirer]
  
  vertex.attr.comb <- list(type=ifelse(revord, 'last', 'first'),
                           name=ifelse(revord, 'last', 'first'))
  
  gi.2 <- igraph::contract.vertices(gi, vmap, vertex.attr.comb = vertex.attr.comb)
  gi.2 <- igraph::simplify(gi.2, remove.multiple = T, remove.loops = T, edge.attr.comb = list(weight='sum'))
  gi.2 <- igraph::induced.subgraph(gi.2, V(gi.2)[igraph::degree(gi.2)>0])
  
  tdf$after[tdf$name!=target] <- power_centrality(gi.2, exponent = decay)
  tdf$delta <- tdf$after - tdf$before
  
  if (verbose)
    print(tdf)
  
  return(list(df=tdf, g=gi.2))
}

##
#
##
mapTo <- function(x,   #vector of degrees
                  minmax=c(9,20),  # vector of length 2: min and max
                  log=F           # logicial for log transform
) { 
  if (any(minmax < 0)) stop ("Negative output range is not allowed.\nPlease assign minmax argument as vector of 2 non-negative values (x>=0) and rerun function.")
  n <- length(x)
  dist <- max(x) - min(x)  #scalar
  #output
  M <- max(minmax)
  m <- min(minmax)
  range <- M - m   # interval to be mapped to
  scale <- (x-min(x)) / dist 
  
  if(log) {
    if(all(x>=1 | x==0)) {
      lx <- log(x)
      maxlx <- max(lx[lx<Inf & lx>-Inf])
      scale <- lx / maxlx
      y <- m + range*scale
      y[is.na(y)|is.nan(y)|y==-Inf|y==Inf] <- m
    } else {
      #augment x proportions to > 1 yielding values suitable for log transform
      scalepos <- (scale+1)/min(scale+1)
      #account for log transform while still proportional to elements of x 
      # by normalizing with the log of the maximum
      scale <- log(scalepos) / log(max(scalepos))
      y <- m + range*scale
    }
  } else {
    y <- m + range*scale
  }
  return(y)
}

centPow <- function(gx, decay=-0.1)
{
  return(power_centrality(gx, exponent = decay))
}

df.pow <- function(gx, betas=c(-.3,-.2,-.1,-.01,0))
{
  df <- data.frame(name=V(gx)$name)
  for (beta in betas) df[ , as.character(beta)] <- centPow(gx, beta)
  return(df)
}

##
# Combine two bipartite networks
#   - keeps original names of vertex and edge properties  (unlike igraph "+" operator:  g3 <- g1 + g2)
##
bipartiteCombine <- function(gx1, gx2) {
  .vt <- unique(rbind(as_data_frame(gx1,'vertices'),as_data_frame(gx2,'vertices')))
  nz <- names(.vt)
  if ((! 'name' %in% nz) & (! 'type' %in% nz)) 
    stop('graphs must have name and type attributes.')
  idx.name <- which(nz=='name')
  idx.type <- which(nz=='type')
  .vt <- .vt[ ,c(idx.name,idx.type)] ## rearrange "name" column first
  .el <- rbind(as_data_frame(gx1,'edges'),as_data_frame(gx2,'edges'))
  gx <- graph.data.frame(d = .el, directed = F, vertices = .vt)
  return(gx)
}

##
# Gets vertex indices of firms in dyads that have multi-market contact 
##
which.mmc <- function(g, focal, keep.focal=F, proj.max=T) {
  if (class(g) != 'igraph') stop('g must be an igraph object')
  if (!igraph::is.weighted(g)) E(g)$weight <- 1
  if (igraph::is.bipartite(g))
    return(which.mmc.bipartite(g, focal, keep.focal, proj.max))
  ## NOT BIPARTITE
  if (! focal %in% 1:vcount(g)) 
    stop('focal firm index must be in vertices of g')
  adjm <- igraph::as_adjacency_matrix(g, attr = 'weight', sparse = F)
  vids <- unname(which(adjm[focal,] > 1))
  if (keep.focal) {
    return(sort(unique(c(vids, focal))))
  } else {
    return(sort(unique(vids)))
  }
}
##
# Gets BIPARTITE graph vertex indices of firms in dyads that have multi-market contact 
##
which.mmc.bipartite <- function(g, focal, keep.focal=F, proj.max=T) {
  g2.l <- bipartite.projection(g, multiplicity = T, remove.type = F)
  vcs <- sapply(g2.l,vcount)
  idx.proj <- ifelse(proj.max, which.max(vcs), which.min(vcs))
  g2 <- g2.l[[idx.proj]]
  if (! focal %in% 1:vcount(g2)) 
    stop('focal firm index must be in vertices of g')
  adjm2 <- igraph::as_adjacency_matrix(g2, attr = 'weight', sparse = F)
  proj.vids <- unname(which(adjm2[focal,] > 1))
  if (keep.focal) {
    proj.vids <- sort(c(proj.vids, focal))
  }
  proj.names <- V(g2)$name[proj.vids]
  ##
  vids.f <- which(V(g)$name %in% proj.names)
  vids.m <- c()
  for (v in vids.f) {
    vids.m <- unique(c(vids.m, as.integer(igraph::neighbors(g, v))))
  }
  vids.focal <- which(V(g)$name==as.character(focal))
  if (keep.focal) {
    return(sort(c(vids.f, vids.m, vids.focal)))
  } else {
    return(sort(c(vids.f, vids.m)))
  }
}

##
# Creates MMC subgraph
#   - subsets to firms with MMC relations to another firm
#   - removes non-MMC edges (weight <= 1)
##
mmcSubgraph <- function(g, focal) {
  vids <- which.mmc(g, focal, keep.focal=T)
  ## MMC VERTEX SUBGRAPH
  g.sub <- igraph::induced.subgraph(g, vids = vids)
  ## DROP NON-MMC EDGES
  edges <- which(E(g.sub)$weight <= 1)
  g.sub <- igraph::delete.edges(g.sub, edges)
  return(g.sub)
}

##
# Gets adjacency matrix -- for either Bipartite or unipartite graphs 
#   - unipartite, just return adjmat
#   - bipartite, return adjmat for the mode with more (proj.max=T) or fewer (proj.max=F) nodes
##
getAdjacencyMatrix <- function(g, focal, proj.max=T) {
  if (igraph::is.bipartite(g)) {
    g2.l <- bipartite.projection(g, multiplicity = T, remove.type = F)
    vcs <- sapply(g2.l,vcount)
    idx.proj <- ifelse(proj.max, which.max(vcs), which.min(vcs))
    g2 <- g2.l[[idx.proj]]
    adjm <- igraph::as_adjacency_matrix(g2, attr = 'weight', sparse = F)
  } else {
    adjm <- igraph::as_adjacency_matrix(g, attr = 'weight', sparse = F)
  }
  return(adjm)
}


mmcSum <- function(g, focal, proj.max=T) {
  adjm <- getAdjacencyMatrix(g, focal, proj.max)
  vids.mmc <- which(adjm[focal,] > 1)
  return(sum(adjm[focal,vids.mmc]))
}


mmcCount <- function(g, focal, proj.max=T) {
  adjm <- getAdjacencyMatrix(g, focal, proj.max)
  return(length(which(adjm[focal,] > 1)))
}

##-----------------------------------------------------------------------------------






##==================================
## FIRM-MARKET GRAPH
##----------------------------------
# ## EXAMPLE OF ALL 4 QUADRANTS
# n1 <- 4
# n2 <- 12
# focal.firm <- as.character(4)
# set.seed(1133241)  #1133241
# gx=sample_bipartite(n1,n2,'gnp',.62)
# ## SPARSE 1
# n1 <- 5
# n2 <- 12
# focal.firm <- as.character(4)
# set.seed(11111)
# gx=sample_bipartite(n1,n2,'gnp',.6)
##
# ## DENSE 2
# n1 <- 4
# n2 <- 12
# focal.firm <- as.character(4)
# set.seed(1133241)
# gx=sample_bipartite(n1,n2,'gnp',.70)
##
##----------------------------------

## Main Cluster
# n1 <- 4   ## markets
# n2 <- 12  ## firms
c1 <- list(m = 4, f = 12)
## Cluster 2
c2 <- list(m = 2, f = 6)


focal.firm <- '4'

## CREATE RANDOM BIPARTITE FIRM_MARKET
set.seed(1133241)  #1133241
gx1 <- sample_bipartite(c1$m,c1$f,'gnp',.62)
V(gx1)$name <- c(LETTERS[1:c1$m], 1:c1$f)
E(gx1)$weight <- 1

set.seed(11341)  #1133241
gx2 <- sample_bipartite(c2$m,c2$f,'gnp',.72)
V(gx2)$name <- c(LETTERS[(c1$m+1):(c1$m+c2$m)], (c1$f+1):(c1$f+c2$f))
E(gx2)$weight <- 1

gx <- bipartiteCombine(gx1, gx2)
# .vt <- unique(rbind(as_data_frame(gx1,'vertices'),as_data_frame(gx2,'vertices')))
# nz <- names(.vt)
# idx.name <- which(nz=='name')
# idx.type <- which(nz=='type')
# .vt <- .vt[ ,c(idx.name,idx.type)] ## rearrange "name" column first
# .el <- rbind(as_data_frame(gx1,'edges'),as_data_frame(gx2,'edges'))
# gx <- graph.data.frame(d = .el, directed = F, vertices = .vt)


## BIMODAL FIRM_MARKET PLOT
vshapes <- sapply(V(gx)$type,function(x)ifelse(x,'circle','square'))
par(mar=c(.1,.1,.1,.1), mfrow=c(1,2))
plot2(gx, 
      layout=layout.bipartite,
      vertex.shape=vshapes,
      vertex.size=18,
      focal.firm=focal.firm)
plot2(gx, 
      layout=layout.kamada.kawai,
      vertex.shape=vshapes,
      vertex.size=18,
      focal.firm=focal.firm, edge.curved = F)
par(mfrow=c(1,1))

## UNIMODAL FIRM_FIRM
gx.ff <- bipartite.projection(gx, remove.type = F)$proj2
V(gx.ff)$type <- unlist(V(gx.ff)$type)

## UNIMODAL FIRM_FIRM ADJACENCY MATRIX
adjmat <- as_adjacency_matrix(gx.ff, attr = 'weight', sparse = F)
print(adjmat)
## SUM MMC
ffidx <- which(V(gx.ff)$name==focal.firm)
mmcidx <- which(adjmat[, ffidx] > 1)
print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
## PLOT FIRM_FIRM NETWORk
vshapes <- sapply(V(gx.ff)$type,function(x)ifelse(x,'circle','square'))
## Save plot of bipartite --> firm-firm competition networ
# png(sprintf("firm_market_firm_firm_2side_N%s_M%s.png",n2,n1), height = 4.5, width = 8.5, units = 'in', res = 250)
par(mar=c(.1,.1,.1,.1), mfrow=c(1,2))
plot2(gx, 
      layout=layout.bipartite,
      vertex.shape=vshapes,
      vertex.size=18,
      focal.firm=focal.firm)
plot2(gx.ff, 
      layout=layout.fruchterman.reingold,
      vertex.shape=vshapes,
      vertex.size= 18, ##1.1*mapTo(centPow(gx.ff, beta = -.01)),
      focal.firm=focal.firm
)
# dev.off()


##==================================
##
## ACQUISITION LOOP -- ALL OTHER FIRMS
##
##
##
##
##
##
##----------------------------------
acq.df <- data.frame(
  name=unlist(V(gx.ff)$name), 
  mmc.sum=NA, 
  mmc.degree=NA, 
  mmc.clust=NA,
  mmc.max.clique=NA,
  mmc.centralization=NA,
  mmc.closeness=NA,
  mmc.eigen=NA,
  mmc.eccen=NA,
  mmc.subgraph=NA,
  exposure=NA, 
  constraint=NA,
  diff.mmc.dyads=NA, 
  diff.mmc.sum=NA, 
  diff.exposure=NA, 
  ##
  stringsAsFactors = F
)
acq.df.0 <- as.list(acq.df[1,])
## base adjacency matrix
adjmat <- getAdjacencyMatrix(gx.ff, attr = 'weight', sparse = F)
## base mmc measures
ffidx <- which(V(gx.ff)$name==focal.firm)
mmcidx <- which(adjmat[, ffidx] > 1)
mmc.dyads0 <- length(adjmat[mmcidx, ffidx])
mmc.sum0 <- sum(adjmat[mmcidx, ffidx])

acq.df.0$mmc.sum <- mmcSum(gx.ff, 4)
acq.df.0$mmc.degree <-  mmcCount(gx.ff, 4)


for (i in 1:vcount(gx.ff)) {
  target.firm <- as.character(i)
  
  if (target.firm != focal.firm) {
    ## NODE COLLAPSE BIPARTITE GRAPH
    tmp <- biAcq(gx, focal.firm, target.firm, project = T)
    gx2.ff <- tmp$g
    exposure.df <- tmp$df
    V(gx2.ff)$type <- unlist(V(gx2.ff)$type)
    V(gx2.ff)$name <- unlist(V(gx2.ff)$name)
    
    ## ADJACENCY
    adjmat <- as_adjacency_matrix(gx2.ff, attr = 'weight', sparse = F)
    ## SUM MMC
    ffidx <- which(V(gx2.ff)$name==focal.firm)
    mmcidx <- which(adjmat[, ffidx] > 1)
    mmc.dyads <- length(adjmat[mmcidx, ffidx])
    mmc.sum <- sum(adjmat[mmcidx, ffidx])
    ## Exposure centrality
    exposure <- exposure.df[which(exposure.df$name==focal.firm), 'after']
    exposure.diff <- exposure.df[which(exposure.df$name==focal.firm), 'delta']
    ## dataframe
    idx <- which(as.character(acq.df$name)==target.firm)
    acq.df$mmc.dyads[idx] <- mmc.dyads
    acq.df$mmc.sum[idx] <- mmc.sum
    acq.df$exposure[idx] <- exposure
    acq.df$diff.mmc.dyads[idx] <- mmc.dyads - mmc.dyads0
    acq.df$diff.mmc.sum[idx] <- mmc.sum - mmc.sum0
    acq.df$diff.exposure[idx] <- exposure.diff
    ## PLOT
    vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
    pngfile <- sprintf("%s\\firm_firm_mmc_acquisition%s_1.png",dirname, target.firm)
    png(pngfile, width = 5, height = 5, units = 'in', res = 250)
    par(mar=c(.1,.1,.1,.1))
    plot2(gx2.ff, 
          layout=layout.fruchterman.reingold,
          vertex.shape=vshapes,
          vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
          focal.firm=focal.firm
    )
    dev.off()
  }
}
print(acq.df)
csvfilename <- sprintf("%s\\acquisition_mmc_synergies_structure_position_compare_M%s_N%s.csv", dirname, n1, n2)
write.csv(acq.df, file = csvfilename)



##=========================
## EXAMPLE HIGH MARKETS
##------------------------

n1 <- 2
n2 <- 12
focal.firm <- as.character(4)

## CREATE RANDOM BIPARTITE FIRM_MARKET
set.seed(1133241)  #1133241
gx=sample_bipartite(n1,n2,'gnp',.72)
V(gx)$name <- c(LETTERS[1:n1], 1:n2)
E(gx)$weight <- 1

## BIMODAL FIRM_MARKET PLOT
vshapes <- sapply(V(gx)$type,function(x)ifelse(x,'circle','square'))
par(mar=c(.1,.1,.1,.1), mfrow=c(1,2))
plot2(gx, 
      layout=layout.bipartite,
      vertex.shape=vshapes,
      vertex.size=18,
      focal.firm=focal.firm)
plot2(gx, 
      layout=layout.kamada.kawai,
      vertex.shape=vshapes,
      vertex.size=18,
      focal.firm=focal.firm, edge.curved = F)
par(mfrow=c(1,1))

## UNIMODAL FIRM_FIRM
gx.ff <- bipartite.projection(gx, remove.type = F)$proj2
V(gx.ff)$type <- unlist(V(gx.ff)$type)
## UNIMODAL FIRM_FIRM ADJACENCY MATRIX
adjmat <- as_adjacency_matrix(gx.ff, attr = 'weight', sparse = F)
print(adjmat)
## SUM MMC
ffidx <- which(V(gx.ff)$name==focal.firm)
mmcidx <- which(adjmat[, ffidx] > 1)
print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
## PLOT FIRM_FIRM NETWORk
vshapes <- sapply(V(gx.ff)$type,function(x)ifelse(x,'circle','square'))
## Save plot of bipartite --> firm-firm competition networ
png(sprintf("firm_market_firm_firm_2side_N%s_M%s.png",n2,n1), height = 4.5, width = 8.5, units = 'in', res = 250)
par(mar=c(.1,.1,.1,.1), mfrow=c(1,2))
plot2(gx, 
      layout=layout.bipartite,
      vertex.shape=vshapes,
      vertex.size=18,
      focal.firm=focal.firm)
plot2(gx.ff, 
      layout=layout.fruchterman.reingold,
      vertex.shape=vshapes,
      vertex.size= 18, ##1.1*mapTo(centPow(gx.ff, beta = -.01)),
      focal.firm=focal.firm
)
dev.off()



# 
# ##==================================
# ## ACQUISITION 3
# ##----------------------------------
# target.firm <- as.character(3)
# ## ACQUISITION UNIMODAL FIRM_FIRM
# gx2.ff <- biAcq(gx, focal.firm, target.firm, project = T)$g
# V(gx2.ff)$type <- unlist(V(gx2.ff)$type)
# V(gx2.ff)$name <- unlist(V(gx2.ff)$name)
# ## ADJACENCY
# adjmat <- as_adjacency_matrix(gx2.ff, attr = 'weight', sparse = F)
# print(adjmat)
# ## SUM MMC
# ffidx <- which(V(gx2.ff)$name==focal.firm)
# mmcidx <- which(adjmat[, ffidx] > 1)
# print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
# ## PLOT
# vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
# plot2(gx2.ff, 
#       layout=layout.fruchterman.reingold,
#       vertex.shape=vshapes,
#       vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
#       focal.firm=focal.firm
# )
# 
# ##==================================
# ## ACQUISITION 6
# ##----------------------------------
# target.firm <- as.character(6)
# ## ACQUISITION UNIMODAL FIRM_FIRM
# gx2.ff <- biAcq(gx, focal.firm, target.firm, project = T)$g
# V(gx2.ff)$type <- unlist(V(gx2.ff)$type)
# V(gx2.ff)$name <- unlist(V(gx2.ff)$name)
# ## ADJACENCY
# adjmat <- as_adjacency_matrix(gx2.ff, attr = 'weight', sparse = F)
# print(adjmat)
# ## SUM MMC
# ffidx <- which(V(gx2.ff)$name==focal.firm)
# mmcidx <- which(adjmat[, ffidx] > 1)
# print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
# ## PLOT
# vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
# plot2(gx2.ff, 
#       layout=layout.fruchterman.reingold,
#       vertex.shape=vshapes,
#       vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
#       focal.firm=focal.firm
# )
# 
# ##==================================
# ## ACQUISITION 2
# ##----------------------------------
# target.firm <- as.character(2)
# ## ACQUISITION UNIMODAL FIRM_FIRM
# gx2.ff <- biAcq(gx, focal.firm, target.firm, project = T)$g
# V(gx2.ff)$type <- unlist(V(gx2.ff)$type)
# V(gx2.ff)$name <- unlist(V(gx2.ff)$name)
# ## ADJACENCY
# adjmat <- as_adjacency_matrix(gx2.ff, attr = 'weight', sparse = F)
# print(adjmat)
# ## SUM MMC
# ffidx <- which(V(gx2.ff)$name==focal.firm)
# mmcidx <- which(adjmat[, ffidx] > 1)
# print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
# ## PLOT
# vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
# plot2(gx2.ff, 
#       layout=layout.fruchterman.reingold,
#       vertex.shape=vshapes,
#       vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
#       focal.firm=focal.firm
# )
# 
# ##==================================
# ## ACQUISITION 5
# ##----------------------------------
# target.firm <- as.character(5)
# ## ACQUISITION UNIMODAL FIRM_FIRM
# gx2.ff <- biAcq(gx, focal.firm, target.firm, project = T)$g
# V(gx2.ff)$type <- unlist(V(gx2.ff)$type)
# V(gx2.ff)$name <- unlist(V(gx2.ff)$name)
# ## ADJACENCY
# adjmat <- as_adjacency_matrix(gx2.ff, attr = 'weight', sparse = F)
# print(adjmat)
# ## SUM MMC
# ffidx <- which(V(gx2.ff)$name==focal.firm)
# mmcidx <- which(adjmat[, ffidx] > 1)
# print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
# ## PLOT
# vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
# plot2(gx2.ff, 
#       layout=layout.fruchterman.reingold,
#       vertex.shape=vshapes,
#       vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
#       focal.firm=focal.firm
# )
# 
# ##==================================
# ## ACQUISITION 10
# ##----------------------------------
# target.firm <- as.character(10)
# ## ACQUISITION UNIMODAL FIRM_FIRM
# gx2.ff <- biAcq(gx, focal.firm, target.firm, project = T)$g
# V(gx2.ff)$type <- unlist(V(gx2.ff)$type)
# V(gx2.ff)$name <- unlist(V(gx2.ff)$name)
# ## ADJACENCY
# adjmat <- as_adjacency_matrix(gx2.ff, attr = 'weight', sparse = F)
# print(adjmat)
# ## SUM MMC
# ffidx <- which(V(gx2.ff)$name==focal.firm)
# mmcidx <- which(adjmat[, ffidx] > 1)
# print(sprintf("FOCAL FIRM %s SUM OF MMC: %s", focal.firm, sum(adjmat[mmcidx, ffidx])))
# ## PLOT
# vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
# plot2(gx2.ff, 
#       layout=layout.fruchterman.reingold,
#       vertex.shape=vshapes,
#       vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
#       focal.firm=focal.firm
# )





##----------- END ------------------

as_adjacency_matrix(bipartite.projection(gx2)$proj2, attr = 'weight', sparse = F)

power_centrality(gx, exponent = -0.2, nodes = V(gx)$type)

biAcq(gi, '1', '2', T)



plot2(gx,vertex.shape=vshapes, layout=layout.kamada.kawai)
# plot(gx,vertex.shape=vshapes, layout=layout.fruchterman.reingold)

E(gx)$weight <- 1
gx.bp <- bipartite.projection(gx)
plot(gx.bp$proj1, edge.width=E(gx.bp$proj1)*.3, vertex.shape='square')
plot(gx.bp$proj2, edge.width=E(gx.bp$proj2)*.025, vertex.shape='circle', 
     layout=layout.fruchterman.reingold)

as_data_frame(gx, what='vertices')

biAcq(gi, '1', '2', T)


##--------------------------------------------------------------------------------------

gx=sample_bipartite(5,5,'gnp',.5)
vshapes <- sapply(V(gx)$type,function(x)ifelse(x,'circle','square'))
plot(gx, 
     layout=layout.bipartite,
     vertex.shape=vshapes)

as_data_frame(gx, what='vertices')


## BIPARTITE EDGE DATAFRAME
df <- data.frame(
  market = c('A','A','A','A', 'B','B','B', 'C','C','C'),
  firm =   c(1,  2,  3,  4,   3,  4,  5,   4,  5,  6)
)

## SPARSE INCIDENCE MATRIX
R <- spMatrix(nrow=length(unique(df$firm)),
              ncol=length(unique(df$market)),
              i = as.numeric(factor(df$firm)),
              j = as.numeric(factor(df$market)),
              x = rep(1, length(as.numeric(df$firm))) )
row.names(R) <- levels(factor(df$firm))
colnames(R) <- levels(factor(df$market))
R

## FIRM_FIRM MATRIX
Rrow <- tcrossprod(R)

## t(R) 
## MODE1::MARKETS
## MODE2: FIRMS
gi <- graph.incidence(t(R))
vshapes <- sapply(V(gi)$type, function(x)ifelse(x,'circle','square'))

plot(gi,  vertex.shape=vshapes)

plot(gi, 
     layout=layout.bipartite,
     vertex.shape=vshapes,
     vertex.size=power_centrality(gi, exponent = -0.2)*30
)

df <- as_data_frame(gi, what='vertices')
for (i in 1:6) {
  for (j in 1:6) {
    if (i != j) {
      df[ ,paste0(i,j)] <- biAcq(gi, as.character(i), as.character(j))$delta
    }
  }
}

## WHICH MIN
apply(df[df$type==T,3:ncol(df)], 2, which.min)

biAcq(gi, '1', '2', T)

## FIRM-FIRM ADJACENCY
ga <- graph.adjacency(tcrossprod(R), diag = F, mode = 'undirected')
E(ga)$weight <- 1
ga <- simplify(ga, remove.multiple = T, remove.loops = T, edge.attr.comb = list(weight='sum'))
set.seed(2)
plot(ga, edge.width=E(ga)$weight^2)




