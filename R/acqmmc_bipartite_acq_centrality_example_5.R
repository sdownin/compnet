
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
  
  vertex.attr.comb <- list(type=function(x)ifelse(revord, x[length(x)], x[1]),
                           name=function(x)ifelse(revord, x[length(x)], x[1]))
  
  if (is.bi) {
    edge.attr.comb <- list(weight=function(x)ifelse(revord, x[length(x)], x[1]))
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
#
##
getNotMmcBipartiteEdgeIds <- function(g.sub)
{
  # adjm <- getAdjacencyMatrix(g.sub, proj.max = T)
  firmnames <- getMaxProjNames(g.sub)
  vids <- which(V(g.sub)$name %in% firmnames)
  # mids <- which( ! V(g.sub)$name %in% firmnames)
  bic <- igraph::bibcoupling(g.sub)  ## shared neighbors (## possibly large matrix, need to subject to only firm vids)
  
  ##------------------------------
  ##  NOT MMC eids
  ##------------------------------
  ## 1. F-F non-MMC dyads
  ffno <- which(bic == 1, arr.ind=T) ## 2-col matrix of (row,col) id tuples for non-mmc elements 
  ffno <- ffno[which(ffno[,1] %in% vids & ffno[,2] %in% vids), ]
  
  ## 2. F-M-F No-MMC 3-paths: cache as (F,M),(M,F) tuples
  fm.no <- c()  ## Firm-Market-Firm No-MMC paths: saved as (F1-M, M-F2, ...)
  urow1 <- unique(ffno[,1])
  cat(' fetching Non-MMC bipartite dyads...\n')
  for (i in 1:length(urow1)) {
    if (i %% 50 == 0) cat(sprintf('  %s (%.2f%s)\n',i,100*i/length(urow1),'%'))
    r1i.r2js <- ffno[which(ffno[,1] == urow1[i]),2]
    xi.paths <- igraph::all_shortest_paths(g.sub, urow1[i], r1i.r2js)$res
    ls <- sapply(xi.paths, length)
    idx <- which(ls==3)
    for (j in idx) {
      x <- as.integer(xi.paths[[j]])
      fm.no <- c(fm.no, c(x[1],x[2], x[2],x[3]))
    }
  }
  cat(' done.\n')
  ## 3. save bipartite F-M edge IDs for non-MMC dyads
  eid.no <- unique(igraph::get.edge.ids(g.sub, vp = fm.no, directed = F, error = F, multi = F))
  
  ##------------------------------
  ##  MMC eids  (filter out)
  ##-------------------------------
  ## 4. MMC firm-firm dyads
  ff.mmc <- which(bic > 1, arr.ind=T) ## 
  ff.mmc <- ff.mmc[which(ff.mmc[,1] %in% vids & ff.mmc[,2] %in% vids), ]
  
  ## 5. F-M-F MMC 3-paths
  fm.mmc <- integer()  ##  cache MMC tuples (F,M),(M,F),(...)
  urow1 <- unique(ff.mmc[,1])
  cat(' fetching MMC bipartite dyads...\n')
  for (i in 1:length(urow1)) {
    if (i %% 50 == 0) cat(sprintf('  %s (%.2f%s)\n',i,100*i/length(urow1),'%'))
    r1i.r2js <- ff.mmc[which(ff.mmc[,1] == urow1[i]),2]
    xi.paths <- igraph::all_shortest_paths(g.sub, urow1[i], r1i.r2js)$res
    ls <- sapply(xi.paths, length)
    idx <- which(ls==3)
    for (j in idx) {
      x <- as.integer(xi.paths[[j]])
      fm.mmc <- c(fm.mmc, c(x[1],x[2], x[2],x[3]))
    }
  }
  cat(' done.\n')
  ## 6. save bipartite F-M edge IDs for MMC dyads
  eid.mmc <- unique(igraph::get.edge.ids(g.sub, vp = fm.mmc, directed = F, error = F, multi = F))
  
  ##--------------------------------
  ## Check is Non-MMC and is NOT MMC
  ##--------------------------------
  ## 7. filter only the Non-MMC F-M dyads that are not included in any MMC F-M dyads (which comprise MMC F-F dyads)
  eids <- eid.no[ which( ! eid.no %in% eid.mmc) ]
  
  return(eids)
}

##
# Creates MMC subgraph
#   - subsets to firms with MMC relations to another firm
#   - removes non-MMC edges (weight <= 1)
##
# ## filter to ego-mmc network (?)
# focal.firm <- which(V(g)$name == focal.name)
# if (length(focal.firm)>0) {
#   ## MMC VERTEX SUBGRAPH if `focal` is set
#   vids <- which.mmc(g, focal.firm, keep.focal=T)
#   g.sub <- igraph::induced.subgraph(g, vids = vids)
# } else {
#   g.sub <- g
# }
##
mmcSubgraph <- function(g, remove.isolates=F) 
{
  is.bi <- is.bipartite.safe(g)

  g.sub <- g
  
  ## DROP NON-MMC EDGES
  if (is.bi) { 
    eids <- getNotMmcBipartiteEdgeIds(g.sub) 
  } else { 
    which(E(g.sub)$weight <= 1) 
  }
  
  if (length(eids) > 0) {
    g.sub <- igraph::delete.edges(g.sub, eids)
  }
  
  if (remove.isolates) {
    g.sub <- igraph::induced.subgraph(g.sub, vids = which(igraph::degree(g.sub)>0) )
  }
  
  return(g.sub)
}

# mmcSubgraphDEBUG <- function(g, focal.name=NA, remove.isolates=F) 
# {
#   is.bi <- is.bipartite.safe(g)
#   focal.firm <- which(V(g)$name == focal.name)
#   
#   # if (length(focal.firm)>0) {
#   #   ## MMC VERTEX SUBGRAPH if `focal` is set
#   #   vids <- which.mmc(g, focal.firm, keep.focal=T)
#   #   g.sub <- igraph::induced.subgraph(g, vids = vids)
#   # } else {
#   #   g.sub <- g
#   # }
#   g.sub <- g
#   
#   ## DROP NON-MMC EDGES
#   if (is.bi) {
#     # adjm <- getAdjacencyMatrix(g.sub, proj.max = T)
#     firmnames <- getMaxProjNames(g.sub)
#     vids <- which(V(g.sub)$name %in% firmnames)
#     # mids <- which( ! V(g.sub)$name %in% firmnames)
#     bic <- igraph::bibcoupling(g.sub)  ## shared neighbors (## possibly large matrix, need to subject to only firm vids)
#     ## MMC firm-firm dyads
#     mmc <- which(bic > 1, arr.ind=T) ## 
#     mmc <- mmc[which(mmc[,1] %in% vids & mmc[,2] %in% vids), ]
#     ## non-MMC firm-firm dyads
#     nmmc <- which(bic == 1, arr.ind=T) ## 2-col matrix of (row,col) id tuples for non-mmc elements 
#     nmmc <- nmmc[which(nmmc[,1] %in% vids & nmmc[,2] %in% vids), ]
#     # ###
#     # bnmids <- sort(unique(c(nmmc[,1],nmmc[,2])))
#     # sapply(bnmids,function(i){
#     #   ls <- sapply(igraph::all_shortest_paths(g.sub, i, vids[ ! vids %in% i])$res, length)
#     #   return(length(ls[ls==3]))
#     # })
#     # ###
#     edge.l <- list()
#     for (i in 1:nrow(nmmc)) {
#       ps <- igraph::all_shortest_paths(g.sub, nmmc[i,1], nmmc[i,2])$res
#       if (length(ps)==1) {
#         pv <- as.integer(ps[[1]])
#         chk1 <- length(which( (mmc[,1]==pv[1] & mmc[,2]==pv[2]) | (mmc[,1]==pv[2] & mmc[,2]==pv[1]) )) == 0
#         chk2 <- length(which( (mmc[,1]==pv[2] & mmc[,2]==pv[3]) | (mmc[,1]==pv[3] & mmc[,2]==pv[2]) )) == 0
#         if (chk1) edge.l <- c(edge.l, list(c(idx[1],idx[2])))
#         if (chk2) edge.l <- c(edge.l, list(c(idx[2],idx[3])))
#       }
#     }
#     eids <- sort(unique(igraph::get.edge.ids(g.sub, vp = edges, directed = F))) ## edge ids of non-mmc firm
#   } else {
#     eids <- which(E(g.sub)$weight <= 1)
#   }
#   
#   if (length(eids) > 0) {
#     g.sub <- igraph::delete.edges(g.sub, eids)
#   }
#   
#   if (remove.isolates) {
#     g.sub <- igraph::induced.subgraph(g.sub, vids = which(igraph::degree(g.sub)==0))
#   }
#   
#   return(g.sub)
# }

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
  bps <- igraph::bipartite.projection.size(gx)
  types <- unique(V(gx)$type)
  idx <- ifelse(bps$vcount1 > bps$vcount2, 1, 2)
  return(V(gx)$name[ which(V(gx)$type == types[idx]) ])
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
    stop(sprintf('vert.name `%s` not in graph gx.m',vert.name))
  
  if (proj.uni) {
    gx.m <- getGraphProjection(gx.m)
    vid.a <- which(V(gx.m)$name==vert.name)
  } else {
    .proj.gx.m <- getGraphProjection(gx.m)
    .proj.vid.a <- which(V(.proj.gx.m)$name==vert.name)
  }
  
  if (!is.na(ego.order) & ego.order >= 1) {
    ord <- ifelse(is.bipartite.safe(gx.m), 2*ego.order, 1*ego.order)  ## bipartite twice distance
    gx.m <- igraph::make_ego_graph(gx.m, ord, vid.a)[[1]]  
    vid.a <- which(V(gx.m)$name == vert.name)
    ##
    ord <- ifelse(is.bipartite.safe(.proj.gx.m), 2*ego.order, 1*ego.order)  ## bipartite twice distance
    .proj.gx.m <- igraph::make_ego_graph(.proj.gx.m, ord, .proj.vid.a)[[1]]  
    .proj.vid.a <- which(V(.proj.gx.m)$name == vert.name)
  }
  
  is.bi <- is.bipartite.safe(gx.m)
  
  df <- data.frame(
    name=unlist(V(gx.m)$name[vid.a]), 
    is.bi=is.bi,
    v=vcount(gx.m),
    e=ecount(gx.m),
    ##
    sum=ifelse(is.bi, getMmcEdgeSum(.proj.gx.m), getMmcEdgeSum(gx.m)),
    degree=ifelse(is.bi, getMmcEdgeCount(.proj.gx.m), getMmcEdgeCount(gx.m)),
    max.clique=ifelse(is.bi, getMaxMmcCliqueSize(.proj.gx.m, .proj.vid.a), getMaxMmcCliqueSize(gx.m, vid.a)),
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

# ## Main Cluster
# # n1 <- 4   ## markets
# # n2 <- 12  ## firms
c1 <- list(m = 4, f = 12)
## Cluster 2
c2 <- list(m = 2, f = 6)

# c1 <- list(m = 10, f = 1200)
# ## Cluster 2
# c2 <- list(m = 5, f = 600)

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

df.a <- data.frame()
df.a.e <- data.frame()

df0.t   <- getMmcDf(gx, focal.name, ego.order=NA)
df0.t.e <- getMmcDf(gx, focal.name, ego.order=1)

df.t.diff <- data.frame()
df.t.e.diff <- data.frame()

meta.attrs <- c('name','targ','is.bi','v','e')
mmc.attrs <- names(df0.t)[which( ! names(df0.t) %in% meta.attrs)]

for (i in gnames) {
  
  ## Acquirer MMC Metrics
  df.a <- rbind(df.a, getMmcDf(gx, as.character(i)) )
  df.a.e <- rbind(df.a.e, getMmcDf(gx, as.character(i), ego.order=1) )
  
  ## Target Synergy MMC Metrics
  target.firm <- as.numeric(i)
  target.name <- as.character(i)
  
  if (target.firm != focal.firm) {
    
    ## NODE COLLAPSE BIPARTITE GRAPH
    gx2 <- biAcq(gx, focal.name, target.name, project = F)
    
    cat(sprintf('%s(%s)-->%s(%s)\n',focal.name, focal.firm, target.name, target.firm))

    ## MMC Subgraph
    gx2.sub <- mmcSubgraph(gx2, remove.isolates=T)
    plot2(gx2.sub)
        
    ## Get MMC metrics
    dfi.t   <- getMmcDf(gx2.sub, focal.name)
    dfi.t.e <- getMmcDf(gx2.sub, focal.name, ego.order=1)
    
    ## make diff df
    dfi.t.diff   <- dfi.t
    dfi.t.e.diff <- dfi.t.e
    
    ## set diff values
    dfi.t.diff[,mmc.attrs]   <- dfi.t[,mmc.attrs] - df0.t[,mmc.attrs]
    dfi.t.e.diff[,mmc.attrs] <- dfi.t.e[,mmc.attrs] - df0.t.e[,mmc.attrs]
    
    ## add target
    dfi.t.diff$targ   <- target.name
    dfi.t.e.diff$targ <- target.name
    
    ## append
    df.t.diff <- rbind(df.t.diff, dfi.t.diff)
    df.t.e.diff <- rbind(df.t.e.diff, dfi.t.e.diff)
    
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

## SAVE DATAFRAMES
print(df.a)
csvfilename <- sprintf("%s\\acquisition_acquirer_mmc_compare_c1M%s_c1N%s_c2M%s_c2N%s.csv", dirname, c1$m, c1$f, c2$m, c2$f)
write.csv(df.a, file = csvfilename)

print(df.a.e)
csvfilename <- sprintf("%s\\acquisition_acquirer_mmc_compare_EGO_c1M%s_c1N%s_c2M%s_c2N%s.csv", dirname, c1$m, c1$f, c2$m, c2$f)
write.csv(df.a.e, file = csvfilename)

print(df.t.diff)
csvfilename <- sprintf("%s\\acquisition_mmc_synergies_structure_position_compare_c1M%s_c1N%s_c2M%s_c2N%s.csv", dirname, c1$m, c1$f, c2$m, c2$f)
write.csv(df.t.diff, file = csvfilename)

print(df.t.e.diff)
csvfilename <- sprintf("%s\\acquisition_mmc_synergies_structure_position_compare_EGO_c1M%s_c1N%s_c2M%s_c2N%s.csv", dirname, c1$m, c1$f, c2$m, c2$f)
write.csv(df.t.e.diff, file = csvfilename)




## ACQUIRER
par(mfrow=c(3,3), mar=c(1,1,2.5,1))
for (attr in mmc.attrs) {
  if (is.numeric(df.a[,attr]))
    hist(df.a[,attr], col='gray', main=attr)
}

## TARGET SYNERGY
par(mfrow=c(3,3), mar=c(1,1,2.5,1))
for (attr in mmc.attrs) {
  x <- df.t.diff[,attr]
  if (is.numeric(x) & length(unique(x)) > 1)
    hist(df.t.diff[,attr], col='gray', main=attr)
}

plot2(gx)

sepattrs <- c()
for (attr in mmc.attrs) {
  if (length(unique(df.t.diff[,attr] < 0)) > 1)
    sepattrs <- c(sepattrs, attr)
}
cat(sprintf('all separating attrs +/-:\n  %s\n\n', paste(sepattrs, collapse = ", ")))
View(df.t.diff[,c(meta.attrs,sepattrs)])

sepattrs <- c()
for (attr in mmc.attrs) {
  if (length(unique(df.t.e.diff[,attr] < 0)) > 1)
    sepattrs <- c(sepattrs, attr)
}
cat(sprintf('EGO separating attrs +/-:\n  %s\n\n', paste(sepattrs, collapse = ", ")))
View(df.t.e.diff[,c(meta.attrs,sepattrs)])








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




