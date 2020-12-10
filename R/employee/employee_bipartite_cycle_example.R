
library(igraph)
library(Matrix)
library(intergraph)


## save plots
dirname <- "C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\envelopment"
setwd(dirname)

## cache default params
.par = par()


##
#
##
plot3 <- function(gx, layout=layout.fruchterman.reingold, vertex.size=15, focal.firm=NA, fam='sans', edge.curved=F, seed=11111, ...)
{
  vAttrs <- igraph::list.vertex.attributes(gx) 
  if ('type' %in% vAttrs) {
    vcolors <- sapply(V(gx)$type, function(x)ifelse(x, "gray","SkyBlue2"))
    lcolors <-  sapply(V(gx)$type, function(x)ifelse(x, "black","darkblue"))
    vshapes <- sapply(1:vcount(gx),function(x)ifelse(V(gx)$type[x], "square", "circle"))
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

##
#  FUNCTION TO TEST CYCLE INTERPRETATION BY MARKETS
##
testBipartiteCycles <- function(mat)
{
  N <- nrow(mat)
  M <- ncol(mat)
  rownames(mat) <- 1:N
  colnames(mat) <- LETTERS[1:M]
  b1 <- graph_from_incidence_matrix(mat, directed = F, add.names = 'name')
  g1 <- bipartite.projection(b1, multiplicity = T)$proj1
  ## MARKET INCIDENCE DISTRIBUTION (diversification)
  degcnt <- plyr::count(igraph::degree(b1)[1:N]) 
  lo <- round(100 * sum(degcnt$freq[which(degcnt$x == 1)]) / N)
  mi <- round(100 * sum(degcnt$freq[which(degcnt$x %in% 2:3)]) / N)
  hi <- round(100 * sum(degcnt$freq[which(degcnt$x > 3)]) / N)
  cy <- sna::kcycle.census(asNetwork(g1), 6)$cycle.count
  mod <- modularity(g1, multilevel.community(g1)$membership)
  ## PLOTTING
  plotname <- sprintf("bipart_cycle_niche_width_N%s_M%s_mod%s_nar%s_mid%s_br%s.png",
                      N,M, round(100*mod), lo,mi,hi)
  png(plotname, height = 5.5, width = 9, units = 'in', res = 200)
    par(mar=c(.1,.1,.1,.1), mfrow=c(1,3), oma=c(0, 0, 2, 0))
    image(t(mat))
    plot3(b1)
    plot3(g1)
    mtext(sprintf("[Narrow:%s|Mid:%s|Broad:%s] (mod=%.2f) 3c=%.1f; 4c=%.1f; 5c=%.1f; 6c=%.1f", 
                  lo,mi,hi, mod, cy[2,1]/N, cy[3,1]/N, cy[4,1]/N, cy[5,1]/N), 
          outer = TRUE, cex = 1.5)
  dev.off()
}


## Disconnected components
## NARROW:High | MIDW:Low | BROAD:None
mat <- matrix(c(1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,0,
                0,1,0,0,1,
                0,0,0,1,0),
              byrow=T,nrow=6)
testBipartiteCycles(mat)

## Star-like hub-n-spoke
## NARROW:High | MIDW:None | BROAD:Low
mat <- matrix(c(1,0,0,0,0,
                0,1,0,0,0,
                0,0,1,0,0,
                0,0,0,1,0,
                1,1,1,1,1,
                0,0,0,0,1),
              byrow=T,nrow=6)
testBipartiteCycles(mat)

## 
## NARROW:Low | MIDW:High | BROAD:None
mat <- matrix(c(1,1,0,0,0,
                0,1,1,0,0,
                0,0,1,1,0,
                1,0,0,1,1,
                0,0,0,1,0,
                0,0,0,1,1),
              byrow=T,nrow=6)
testBipartiteCycles(mat)

## NARROW:Low | MIDW:High | BROAD:Low
mat <- matrix(c(1,1,0,0,0,
                0,1,1,0,0,
                0,0,1,1,0,
                1,0,0,1,1,
                0,0,0,1,0,
                1,1,1,1,1),
              byrow=T,nrow=6)
testBipartiteCycles(mat)

## NARROW:Low | MIDW:Low | BROAD:High
mat <- matrix(c(1,0,0,0,0,
                0,1,1,1,1,
                0,0,1,1,0,
                1,1,1,0,1,
                1,1,1,1,1,
                1,1,0,1,1),
              byrow=T,nrow=6)
testBipartiteCycles(mat)

## NARROW:None | MIDW:Low | BROAD:High
mat <- matrix(c(1,0,1,1,1,
                0,1,1,1,1,
                0,0,1,1,0,
                1,1,1,0,1,
                1,1,1,1,1,
                1,1,0,1,1),
              byrow=T,nrow=6)
testBipartiteCycles(mat)


## Core-periphery
## NARROW:Mid | MIDW:None | BROAD:Mid
mat <- matrix(c(1,1,1,1,1,
                0,1,1,1,1,
                1,1,1,1,0,
                0,1,0,0,0,
                0,0,0,0,1,
                0,0,0,1,0),
              byrow=T,nrow=6)
testBipartiteCycles(mat)

## Small Worlds
## NARROW:Low | MIDW:High | BROAD:None
mat <- matrix(c(1,1,1,0,0,
                0,1,0,0,0,
                1,1,1,0,0,
                0,1,0,1,0,
                0,0,0,1,1,
                0,0,0,1,1),
              byrow=T,nrow=6)
testBipartiteCycles(mat)









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




