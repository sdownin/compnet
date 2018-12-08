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
  vAttrs <- igraph::list.vertex.attributes(gx) 
  if ('type' %in% vAttrs) {
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
       edge.width = 1 + 2 * (E(gx)$weight-1),
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


# ##
# # Bipartite Graph Acquisition -- PREVIOUS VERSION
# ##
# biAcq.prev <- function(gi, acquirer, target, decay=-0.2, project=T, verbose=T)
# {
#   if (project) {
#     gi.l <- bipartite.projection(gi, multiplicity = T, remove.type = F)
#     vcs <- sapply(gi.l,vcount)
#     gi <- gi.l[[ which.max(vcs) ]]
#     V(gi)$type <- unlist(V(gi)$type)
#     V(gi)$name <- unlist(V(gi)$name)
#   }
#   
#   tdf <- as_data_frame(gi, what='vertices')
#   tdf$before <- power_centrality(gi, exponent = decay)
#   tdf$after <- NA
#   
#   vnamemap <- names(V(gi))
#   vmap <- as.integer(V(gi))
#   
#   revord <- which(vnamemap==target) < which(vnamemap==acquirer)
#   comb.func <- ifelse(revord, 'last', 'first')
#   
#   vmap[vnamemap==target] <- vmap[vnamemap==acquirer]
#   
#   vertex.attr.comb <- list(type=ifelse(revord, 'last', 'first'),
#                            name=ifelse(revord, 'last', 'first'))
#   
#   gi.2 <- igraph::contract.vertices(gi, vmap, vertex.attr.comb = vertex.attr.comb)
#   gi.2 <- igraph::simplify(gi.2, remove.multiple = T, remove.loops = T, edge.attr.comb = list(weight='sum'))
#   gi.2 <- igraph::induced.subgraph(gi.2, V(gi.2)[igraph::degree(gi.2)>0])
#   
#   tdf$after[tdf$name!=target] <- power_centrality(gi.2, exponent = decay)
#   tdf$delta <- tdf$after - tdf$before
#   
#   if (verbose)
#     print(tdf)
#   
#   return(list(df=tdf, g=gi.2))
# }

##
# Bipartite Graph Acquisition
##
biAcq <- function(gi, acquirer.name, target.name, project=F, verbose=T)
{
  if ( ! 'name' %in% igraph::list.vertex.attributes(gi)) 
    stop('gi must have name attribute.')
  
  is.bi <- is.bipartite.safe(gi)
  
  if (project & is.bi) {
    gi.l <- bipartite.projection(gi, multiplicity = T, remove.type = F)
    vcs <- sapply(gi.l,vcount)
    gi <- gi.l[[ which.max(vcs) ]]
    V(gi)$type <- unlist(V(gi)$type)
    V(gi)$name <- unlist(V(gi)$name)
    is.bi <- is.bipartite.safe(gi)
  }
  
  acquirer <- which(V(gi)$name==acquirer.name)
  target   <- which(V(gi)$name==target.name)
  
  if (length(acquirer)==0 | length(target)==0) {
    stop(sprintf('has acquirer=%s; has target=%s', length(acquirer)>0, length(target)>0))
  }
  
  vnamemap <- names(V(gi))
  vmap <- as.integer(V(gi))
  
  revord <- which(vnamemap==target.name) < which(vnamemap==acquirer.name)
  comb.func <- ifelse(revord, 'last', 'first')
  
  vmap[vnamemap==target.name] <- vmap[vnamemap==acquirer.name]
  
  vertex.attr.comb <- list(type=ifelse(revord, 'last', 'first'),
                           name=ifelse(revord, 'last', 'first'))
  
  if (is.bi) {
    edge.attr.comb <- list(weight=ifelse(revord, 'last', 'first'))
  } else {
    edge.attr.comb <- list(weight='sum')
  }
  
  gi.2 <- igraph::contract.vertices(gi, vmap, vertex.attr.comb = vertex.attr.comb)
  gi.2 <- igraph::simplify(gi.2, remove.multiple = T, remove.loops = T, edge.attr.comb = edge.attr.comb)
  gi.2 <- igraph::induced.subgraph(gi.2, V(gi.2)[igraph::degree(gi.2)>0])
  
  return(gi.2)
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
# checks if graph is actually bipartite by `type` attribute
#   - if only one `type` then not functionally bipartite
#   - if more than one `type` then is functionally biparite
##
is.bipartite.safe <- function(g)
{
  if (!igraph::is.bipartite(g) | ! 'type' %in% igraph::list.vertex.attributes(g)) 
    return(FALSE)
  return(length(unique(V(g)$type)) > 1)
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
  idx.rest <- which( ! nz %in% c('name','type'))
  .vt <- .vt[ ,c(idx.name,idx.type,idx.rest)] ## rearrange "name" column first
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
  if (is.bipartite.safe(g))
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
mmcSubgraph <- function(g, focal=NA) {
  is.bi <- is.bipartite.safe(g)
  
  if (!is.na(focal)) {
    ## MMC VERTEX SUBGRAPH if `focal` is set
    vids <- which.mmc(g, focal, keep.focal=T)
    g.sub <- igraph::induced.subgraph(g, vids = vids)
  } else {
    g.sub <- g
  }
  
  ## DROP NON-MMC EDGES
  if (is.bi) {
    g.sub.proj <- getGraphProjection(g.sub)
    eids <- which( E(g.sub.proj)$weight > 1 )
    edf <- as.data.frame(igraph::get.edges(g.sub.proj, eids))
    ##
    ## TODO: FIX FOR BIPARTITE NETWORKS
    ##       IF CANNOT REMOVE NONMMC-BIPARTITE EDGES SEPARATE FROM NODES
    #        THEN ONLY REMOVE NONMMC-NODES
    ##
    uvids <- unique(c(edf[,1],edf[,2]))
    ##
    bi.edf <- edf
    bi.eids <- c()
    for (i in 1:nrow(edf)) {
      name1 <- V(g.sub.proj)$name[ edf[i,1] ]
      name2 <- V(g.sub.proj)$name[ edf[i,2] ]
      bi.edf[i,1] <- which(V(g.sub)$name == name1)
      bi.edf[i,2] <- which(V(g.sub)$name == name2)
      bi.eids <- c(bi.eids, igraph::get.edge.ids(g.sub, c(bi.edf[i,1], bi.edf[i,2])))
    }
  } else {
    edges <- which(E(g.sub)$weight <= 1)
  }
  
  if (length(edges) > 0) {
    g.sub <- igraph::delete.edges(g.sub, edges)
  }
  
  return(g.sub)
}

##
# Gets adjacency matrix -- for either Bipartite or unipartite graphs 
#   - unipartite, just return adjmat
#   - bipartite, return adjmat for the mode with more (proj.max=T) or fewer (proj.max=F) nodes
##
getAdjacencyMatrix <- function(g, proj.max=T) {
  if (is.bipartite.safe(g)) {
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

##
#
##
getGraphProjection <- function(g, max.proj=T, remove.type=T)
{
  if (!is.bipartite.safe(g))
    return(g)
  g2.l <- bipartite.projection(g, multiplicity = T, remove.type = remove.type)
  vcs <- sapply(g2.l,vcount)
  which.proj <- ifelse(max.proj, which.max, which.min)
  g2 <- g2.l[[ which.proj(vcs) ]]
  return(g2)
}

## 
# Get vertex IDs 
#  - if bipartite, return vids of maximal projection (largest size mode)
##
getMaxProjNames <- function(gx)
{
  if (! 'name' %in% igraph::list.vertex.attributes(gx))
    stop('gx must have vertex name attribute')
  g2 <- getGraphProjection(gx)
  return(V(g2)$name)
}

# mmcSum <- function(g, focal, proj.max=T) {
#   adjm <- getAdjacencyMatrix(g, focal, proj.max)
#   vids.mmc <- which(adjm[focal,] > 1)
#   return(sum(adjm[focal,vids.mmc]))
# }
# 
# 
# mmcCount <- function(g, focal, proj.max=T) {
#   adjm <- getAdjacencyMatrix(g, focal, proj.max)
#   return(length(which(adjm[focal,] > 1)))
# }

getMmcEdgeSum <- function(g, name.remove=T) {
  if ( ! 'weight' %in% igraph::list.edge.attributes(g)) 
    stop('g must have edge weights')
  adj <- getAdjacencyMatrix(g, T)
  sums <- apply(adj, 1, function(x)  sum(x[x>1]) )
  if (name.remove)
    sums <- unname(sums)
  return(sums)
}

getMmcEdgeCount <- function(g, name.remove=T) {
  if ( ! 'weight' %in% igraph::list.edge.attributes(g)) 
    stop('g must have edge weights')
  adj <- getAdjacencyMatrix(g, T)
  counts <- apply(adj, 1, function(x) length(x[x>1]) )
  if (name.remove)
    counts <- unname(counts)
  return(counts)
}

getMaxMmcCliqueSize <- function(g, vid, min=3)
{
  cls <- igraph::cliques(g, min = min)
  if (length(cls)==0)
    return(0)
  idx <- which(sapply(cls,function(x) vid %in% x))
  return(max(sapply(cls[idx], length)))
}



##
#
##
# getMmcTargetDataframe <- function(gx.m, vid.a, is.ego=FALSE)
# {
#   if (is.ego) {
#     gx.m <- igraph::make_ego_graph(gx.m, 1, vid.a)[[1]]  
#     vid.a <- which(as.character(V(gx.m)$name) == as.character(vid.a))
#   }
#   return(data.frame(
#     name=unlist(V(gx.m)$name), 
#     ##
#     sum=mmcEdgeSum(gx.m)[vid.a], ## sum of mmc
#     degree=mmcEdgeCount(gx.m)[vid.a],  ## number of mmc competitiors 
#     ##
#     clust=igraph::transitivity(gx.m, type = 'global'),
#     closeness=unname(igraph::closeness(gx.m, vid.a)),
#     eigen=unname(igraph::eigen_centrality(gx.m)$vector[vid.a]),
#     pow.n1=unname(igraph::power_centrality(gx.m, vid.a, exponent = -0.1)),
#     pow.n3=unname(igraph::power_centrality(gx.m, vid.a, exponent = -0.3)),
#     eccen=unname(igraph::eccentricity(gx.m, vid.a)),
#     ##
#     central.clos=igraph::centr_clo(gx.m)$centralization / igraph::centr_clo_tmax(gx.m),
#     central.eign=igraph::centr_eigen(gx.m)$centralization / igraph::centr_eigen_tmax(gx.m),
#     central.betw=igraph::centr_betw(gx.m)$centralization / igraph::centr_betw_tmax(gx.m),
#     central.degr=igraph::centr_degree(gx.m)$centralization / igraph::centr_degree_tmax(gx.m),
#     ##
#     subgraph=unname(igraph::subgraph.centrality(gx.m)[vid.a]),
#     density=igraph::graph.density(gx.m), 
#     constraint=unname(igraph::constraint(gx.m, vid.a)),
#     max.clique=maxCliqueSize(gx.m, vid.a),
#     ##
#     stringsAsFactors = F
#   ))
# }

##
#
##
getMmcDf <- function(gx.m, vert.name, ego.order=NA, proj.uni=FALSE)
{
  vid.a <- which(V(gx.m)$name==vert.name)
  if (length(vid.a)==0) 
    stop('vert.name not in graph gx.m')
  if (proj.uni) {
    gx.m <- getGraphProjection(gx.m)
    vid.a <- which(V(gx.m)$name==vert.name)
  }
  if (!is.na(ego.order) & ego.order >= 1) {
    ord <- ifelse(is.bipartite.safe(gx.m), 2*ego.order, 1*ego.order)  ## bipartite twice distance
    gx.m <- igraph::make_ego_graph(gx.m, ord, vid.a)[[1]]  
    vid.a <- which(V(gx.m)$name == vert.name)
  }
  is.bi <- is.bipartite.safe(gx.m)
  df <- data.frame(
    name=unlist(V(gx.m)$name[vid.a]), 
    is.bi=is.bi,
    v=vcount(gx.m),
    e=ecount(gx.m),
    ##
    sum=ifelse(is.bi, NA, getMmcEdgeSum(gx.m)),
    degree=ifelse(is.bi, NA, getMmcEdgeCount(gx.m)),
    max.clique=ifelse(is.bi, NA, getMaxMmcCliqueSize(gx.m, vid.a)),
    ##
    clust=igraph::transitivity(gx.m, type = 'global'),
    closeness=unname(igraph::closeness(gx.m, vid.a)),
    eigen=unname(igraph::eigen_centrality(gx.m)$vector[vid.a]),
    pow.n1=unname(igraph::power_centrality(gx.m, vid.a, exponent = -0.1)),
    pow.n3=unname(igraph::power_centrality(gx.m, vid.a, exponent = -0.3)),
    eccen=unname(igraph::eccentricity(gx.m, vid.a)),
    ##
    central.clos=igraph::centr_clo(gx.m)$centralization / igraph::centr_clo_tmax(gx.m),
    central.eign=igraph::centr_eigen(gx.m)$centralization / igraph::centr_eigen_tmax(gx.m),
    central.betw=igraph::centr_betw(gx.m)$centralization / igraph::centr_betw_tmax(gx.m),
    central.degr=igraph::centr_degree(gx.m)$centralization / igraph::centr_degree_tmax(gx.m),
    ##
    subgraph=unname(igraph::subgraph.centrality(gx.m)[vid.a]),
    density=igraph::graph.density(gx.m), 
    constraint=unname(igraph::constraint(gx.m, vid.a)),
    ##
    stringsAsFactors = F
  )
  return(df)
}

# ##
# #
# ##
# getMmcAcquirerDf <- function(gx.m, vert.name, is.ego=FALSE)
# {
#   vid.a <- which(V(gx.m)$name==vert.name)
#   if (is.ego) {
#     ord <- ifelse(is.bipartite.safe(gx.m), 2, 1)  ## 
#     gx.m <- igraph::make_ego_graph(gx.m, ord, vid.a)[[1]]  
#     vid.a <- which(as.character(V(gx.m)$name) == as.character(vid.a))
#   }
#   return(data.frame(
#     name=unlist(V(gx.m)$name[vid.a]), 
#     ## 
#     sum=mmcEdgeSum(gx.m)[vid.a], ## sum of mmc
#     degree=mmcEdgeCount(gx.m)[vid.a],  ## number of mmc competitiors
#     max.clique=maxCliqueSize(gx.m, vid.a),
#     ##
#     clust=igraph::transitivity(gx.m, type = 'global'),
#     closeness=unname(igraph::closeness(gx.m, vid.a)),
#     eigen=unname(igraph::eigen_centrality(gx.m)$vector[vid.a]),
#     pow.n1=unname(igraph::power_centrality(gx.m, vid.a, exponent = -0.1)),
#     pow.n3=unname(igraph::power_centrality(gx.m, vid.a, exponent = -0.3)),
#     eccen=unname(igraph::eccentricity(gx.m, vid.a)),
#     ##
#     central.clos=igraph::centr_clo(gx.m)$centralization / igraph::centr_clo_tmax(gx.m),
#     central.eign=igraph::centr_eigen(gx.m)$centralization / igraph::centr_eigen_tmax(gx.m),
#     central.betw=igraph::centr_betw(gx.m)$centralization / igraph::centr_betw_tmax(gx.m),
#     central.degr=igraph::centr_degree(gx.m)$centralization / igraph::centr_degree_tmax(gx.m),
#     ##
#     subgraph=unname(igraph::subgraph.centrality(gx.m)[vid.a]),
#     density=igraph::graph.density(gx.m), 
#     constraint=unname(igraph::constraint(gx.m, vid.a)),
#     ##
#     stringsAsFactors = F
#   ))
# }

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
gx_1 <- sample_bipartite(c1$m,c1$f,'gnp',.62)
V(gx_1)$name <- c(LETTERS[1:c1$m], 1:c1$f)
E(gx_1)$weight <- 1

set.seed(11341)  #1133241
gx_2 <- sample_bipartite(c2$m,c2$f,'gnp',.72)
V(gx_2)$name <- c(LETTERS[(c1$m+1):(c1$m+c2$m)], (c1$f+1):(c1$f+c2$f))
E(gx_2)$weight <- 1

## COMBINE
gx <- bipartiteCombine(gx_1, gx_2)

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

vid.a <- 4
vid.ts <- c(15,16)
par(mfrow=c(2,3), mar=c(.1,.1,1.5,.1))
for (vid.t in vid.ts) 
{
  plot2(gx, main="Pre-Acquisition")
  plot2(biAcq(gx, vid.a, vid.t), main=sprintf("%s==>%s",vid.a,vid.t))
  plot2(mmcSubgraph(biAcq(gx, vid.a, vid.t), vid.a), 
        main=sprintf("%s==>%s MMC Subgraph",vid.a,vid.t))
}


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

## focal firm
focal.firm <- 4
focal.name <- as.character(focal.firm)

## vert names for either uni or bipartite
gnames <- getMaxProjNames(gx)

df0.a   <- getMmcDf(gx, focal.name, ego.order=NA)
df0.a.e <- getMmcDf(gx, focal.name, ego.order=1)

df.a.diff <- data.frame()
df.a.e.diff <- data.frame()

meta.attrs <- c('name','targ','is.bi','v','e')
mmc.attrs <- names(df0.a)[which( ! names(df0.a) %in% meta.attrs)]

for (i in gnames) {
  
  target.firm <- as.numeric(i)
  target.name <- as.character(i)
  
  if (target.firm != focal.firm) {
    
    ## NODE COLLAPSE BIPARTITE GRAPH
    gx2 <- biAcq(gx, focal.name, target.name, project = F)
    # plot2(gx2)
    
    V(gx2)$type <- unlist(V(gx2)$type)
    V(gx2)$name <- unlist(V(gx2)$name)
    
    dfi.a   <- getMmcDf(gx2, focal.name, ego.order=NA)
    dfi.a.e <- getMmcDf(gx2, focal.name, ego.order=1)
    
    ## make diff df
    dfi.a.diff   <- dfi.a
    dfi.a.e.diff <- dfi.a.e
    
    ## set diff values
    dfi.a.diff[,mmc.attrs]   <- dfi.a[,mmc.attrs] - df0.a[,mmc.attrs]
    dfi.a.e.diff[,mmc.attrs] <- dfi.a.e[,mmc.attrs] - df0.a.e[,mmc.attrs]
    
    ## add target
    dfi.a.diff$targ   <- target.name
    dfi.a.e.diff$targ <- target.name
    
    ## append
    df.a.diff <- rbind(df.a.diff, dfi.a.diff)
    df.a.e.diff <- rbind(df.a.e.diff, dfi.a.e.diff)
    
    ## dataframe
    # idx <- which(as.character(acq.df$name)==target.firm)
    
    
    # ## PLOT
    # vshapes <- sapply(V(gx2.ff)$type,function(x)ifelse(x,'circle','square'))
    # pngfile <- sprintf("%s\\firm_firm_mmc_acquisition%s_1.png",dirname, target.firm)
    # png(pngfile, width = 5, height = 5, units = 'in', res = 250)
    #   par(mar=c(.1,.1,.1,.1))
    #   plot2(gx2.ff, 
    #         layout=layout.fruchterman.reingold,
    #         vertex.shape=vshapes,
    #         vertex.size= 18, ##1.1*mapTo(centPow(gx2.ff, beta = -.1)),
    #         focal.firm=focal.firm
    #   )
    # dev.off()
  }
  
}

plot2(gx)

sepattrs <- c()
for (attr in mmc.attrs) {
  if (length(unique(df.a.diff[,attr] < 0)) > 1)
    sepattrs <- c(sepattrs, attr)
}
cat(sprintf('all separating attrs +/-:\n  %s\n\n', paste(sepattrs, collapse = ", ")))
View(df.a.diff[,c(meta.attrs,sepattrs)])

sepattrs <- c()
for (attr in mmc.attrs) {
  if (length(unique(df.a.e.diff[,attr] < 0)) > 1)
    sepattrs <- c(sepattrs, attr)
}
cat(sprintf('EGO separating attrs +/-:\n  %s\n\n', paste(sepattrs, collapse = ", ")))
View(df.a.e.diff[,c(meta.attrs,sepattrs)])


print(df.a.diff)
csvfilename <- sprintf("%s\\acquisition_mmc_synergies_structure_position_compare_c1M%s_c1N%s_c2M%s_c2N%s.csv", dirname, c1$m, c1$f, c2$m, c2$f)
write.csv(df.a.diff, file = csvfilename)

print(df.a.e.diff)
csvfilename <- sprintf("%s\\acquisition_mmc_synergies_structure_position_compare_EGO_c1M%s_c1N%s_c2M%s_c2N%s.csv", dirname, c1$m, c1$f, c2$m, c2$f)
write.csv(df.a.e.diff, file = csvfilename)







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




