##------------------------------------
##
##
##
##-------------------------------------
work_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2"
setwd(work_dir)

library(parallel)
library(btergm)
library(igraph)
library(ggplot2)
library(intergraph)
library(strinr)
library(psych)
# library(doParallel)
# library(jtools)
# library(interactions)

# registerDoParallel(cores = 4)
# system.time({
#    foreach(i=1:4000) %dopar% {
#       sqrt(i)
#    }
# })
# 
# system.time({
#    for (i in 1:4000) sqrt(i)
# })

##----BOOTSTRAP REPLICATATIONS-----
R <- 200
ncpus <- detectCores()
parallel <- "multicore"
##---------------------------------

## Sample Time periods
periods <- list(2008:2009,2010:2011,2012:2013,2014:2015)
startYr <- min(unlist(periods))
endYr <- max(unlist(periods))

##---------------------------------
## Data load 
##---------------------------------
empnetsfile <- sprintf('employee_competition_networks_lists_pd%s_yr%s_%s-%s.rds',
                       length(periods),endYr-startYr+1,startYr,endYr)
# saveRDS(list(geli=geli, gcli=gcli), file = empnetsfile)
.nets <- readRDS(empnetsfile)
geli <- .nets$geli
gcli <- .nets$gcli


npds <- length(geli)

ttriple <- lapply(geli, function(net) as.matrix(net %n% 'ttriple') )
cycle3 <- lapply(geli, function(net) as.matrix(net %n% 'cycle3') )
status <- lapply(geli, function(net){
   x <- net %v% 'eigen_centrality'
   return(outer(x,x,'+'))   ## change statistics are defined as sum{V_i, V_j}
})
#
samefmclust <- lapply(geli, function(net){
   x <- net %v% 'cluster_between'
   return(outer(x,x,'==') * 1)
})
tofrexecdiff <- lapply(geli, function(net) as.matrix(net %n% 'to_from_exec_role_diff_avg') )
agediff <- lapply(geli, function(net) as.matrix(net %n% 'age_diff') * -1 )  ## multiply by -1 for forward fiff [to]-[from]
mmc <- lapply(gcli, function(net) as.matrix(net %n% 'mmc') )
cossim <- lapply(gcli, function(net) as.matrix(net %n% 'cat_cos_sim') )
pmuncommon <- lapply(gcli, function(net) { # product market uncommonality dummy
   D <- as.matrix(igraph::distances(asIgraph(net)))
   PMU <- D
   diag(PMU) <- 0
   PMU[D == 1] <- 0
   PMU[D > 1] <- 1
   return(PMU)
})
##---------- A. INTERACTIONS: COMPNET TIE-------------------------------
a.pmu.status             <- lapply(1:npds, function(i)  pmuncommon[[i]] * status[[i]] )
a.pmu.cossim             <- lapply(1:npds, function(i)  pmuncommon[[i]] * cossim[[i]] )
# a.pmu.samefmclust        <- lapply(1:npds, function(i)  pmuncommon[[i]] * samefmclust[[i]] )
a.pmu.cycle3             <- lapply(1:npds, function(i)  pmuncommon[[i]] * cycle3[[i]] )
a.pmu.ttriple            <- lapply(1:npds, function(i)  pmuncommon[[i]] * ttriple[[i]] )



##-----------------------------------------------
## Model Load
##-----------------------------------------------
structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- readRDS(file=sprintf('%s.rds',filebase))
gc()



# [1] "edges"                     
# [2] "edgecov.cycle3[[i]]"       
# [3] "edgecov.ttriple[[i]]"      
# [4] "edgecov.memory[[i]]"       
# [5] "edgecov.delrecip[[i]]"     
# [6] "edgecov.tofrexecdiff[[i]]" 
# [7] "nodecov.age"               
# [8] "nodecov.ipo_status"        
# [9] "edgecov.status[[i]]"       
# [10] "edgecov.agediff[[i]]"      
# [11] "edgecov.mmc[[i]]"          
# [12] "nodematch.state_code"      
# [13] "edgecov.cossim[[i]]"       
# [14] "edgecov.pmuncommon[[i]]"   
# [15] "edgecov.a.pmu.status[[i]]" 
# [16] "edgecov.a.pmu.cossim[[i]]" 
# [17] "edgecov.a.pmu.cycle3[[i]]" 
# [18] "edgecov.a.pmu.ttriple[[i]]"
##-----------------------------------------------
## Names mapping
##-----------------------------------------------
## COEF NAME MAP
coef.name.map <- list(
  `edges`='Edges (network constant)',
  # EMPLOYEE (HUMAN CAPITAL)
  `edgecov.tofrexecdiff[[i]]`='Executive Role Promotion',
  `edgecov.to.fr.exec.diff[[i]]`='Executive Role Promotion',
  `edgecov.frtenure[[i]]`='From Job Tenure (years)',
  `edgecov.frexec[[i]]`='From Executive Role',
  `edgecov.toexec[[i]]`='To Executive Role',
  `edgecov.female[[i]]`='Female Ratio',
  # FIRM
  `nodecov.age`='Firms Age',
  `nodecov.ipo_status`='Public Firm (1=Yes)',
  `edgecov.status[[i]]`='Status',
  # DYAD
  `nodematch.ipo_status`='Ownership Status Homophily',
  `nodematch.state_code`='HQ Region Homophily',
  `edgecov.agediff[[i]]`='Age Difference',
  `edgecov.mmc[[i]]`='Firm Branch MMC',
  # DYAD - PRODUCT MARKET
  `edgecov.invdist[[i]]`='Product Market Closeness (inverse distance)',
  # Dyad - FACTOR MARKET
  `edgecov.cossim[[i]]`='Category Similarity',
  # STRUCTURE
  `ctriple`='Cyclic Triads', `edgecov.cycle3[[i]]`='Cyclic Triads',
  `ttriple`='Transitive Triads', `edgecov.ttriple[[i]]`='Transitive Triads',
  # TEMPORAL
  `edgecov.delrecip[[i]]`='Delayed Reciprocity',
  `edgecov.memory[[i]]`='Temporal Stability',
  #
  `edgecov.cossim[[i]]` = 'Category Similarity',   
  `edgecov.pmuncommon[[i]]`= 'PMU',
  `edgecov.a.pmu.status[[i]]`= 'PMU * Status',
  `edgecov.a.pmu.cossim[[i]]`='PMU * Category Similarity',
  `edgecov.a.pmu.cycle3[[i]]`= 'PMU * Cyclic Triads',
  `edgecov.a.pmu.ttriple[[i]]`='PMU * Transitive Triads',
  # DV
  `Y` = 'Personnel Mobility',
  `DV` = 'Personnel Mobility'
)







##-----------------------------------------------
## Correlatons and Descriptives
##-----------------------------------------------
X <- m5a@effects
X <- X[,which(names(X) != 'edges')]
df <- cbind(Y=m5a@response, X)
names <- names(df)
for (i in 1:length(names)) {
   names[i] <- coef.name.map[[ names[i] ]]
}
names(df) <- names

xdes <- psych::describe(df)
xcor <- cor(df)

write.csv(xdes, file = sprintf('descriptives_%s.csv',filebase))
write.csv(xcor, file = sprintf('correlations_%s.csv',filebase))





##-----------------------------------------------
## GOF
##-----------------------------------------------
nsim <- 30

## GOF-----------------------------------
gf.dsp <- btergm::gof(m5a, statistics = c(dsp), nsim=nsim)
saveRDS(gf.dsp, file = sprintf('gof_dsp_%s.rds',filebase))
write.csv(gf.dsp[[1]]$stats, file = sprintf('gof_dsp_%s.csv',filebase))
png(filename = sprintf('gof_dsp_%s.png',filebase), width = 6.5, height = 5.5, units = 'in', res = 200)
   plot(gf.dsp)
dev.off()
rm(gf.dsp)
gc()

gf.esp <- btergm::gof(m5a, statistics = c(esp), nsim=nsim)
saveRDS(gf.esp, file = sprintf('gof_esp_%s.rds',filebase))
write.csv(gf.esp[[1]]$stats, file = sprintf('gof_esp_%s.csv',filebase))
png(filename = sprintf('gof_esp_%s.png',filebase), width = 6.5, height = 5.5, units = 'in', res = 200)
   plot(gf.esp)
dev.off()
rm(gf.esp)
gc()

gf.odeg <- btergm::gof(m5a, statistics = c(odeg), nsim=nsim)
saveRDS(gf.odeg, file = sprintf('gof_odeg_%s.rds',filebase))
write.csv(gf.odeg[[1]]$stats, file = sprintf('gof_odeg_%s.csv',filebase))
png(filename = sprintf('gof_odeg_%s.png',filebase), width = 6.5, height = 5.5, units = 'in', res = 200)
   plot(gf.odeg)
dev.off()
rm(gf.odeg)
gc()

gf.ideg <- btergm::gof(m5a, statistics = c(ideg), nsim=nsim)
saveRDS(gf.ideg, file = sprintf('gof_ideg_%s.rds',filebase))
write.csv(gf.ideg[[1]]$stats, file = sprintf('gof_ideg_%s.csv',filebase))
png(filename = sprintf('gof_ideg_%s.png',filebase), width = 6.5, height = 5.5, units = 'in', res = 200)
   plot(gf.ideg)
dev.off()
rm(gf.ideg)
gc()

gf.geodesic <- btergm::gof(m5a, statistics = c(geodesic), nsim=nsim)
saveRDS(gf.geodesic, file = sprintf('gof_geodesic_%s.rds',filebase))
write.csv(gf.geodesic[[1]]$stats, file = sprintf('gof_geodesic_%s.csv',filebase))
png(filename = sprintf('gof_geodesic_%s.png',filebase), width = 6.5, height = 5.5, units = 'in', res = 200)
   plot(gf.geodesic)
dev.off()
rm(gf.geodesic)
gc()

rm(m5a)
gc()
## -------------------------------------------------------





##--------------------------------------
##    INTERPRET
##--------------------------------------

structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
# fit <- readRDS(file=sprintf('%s.rds',filebase))
fit <- m5a
rm(m5a)
gc()



## FUNCTIONS FOR NAs REMOVED
.med <- function(x){return(median(x, na.rm=TRUE))}
.min <- function(x){return(min(x, na.rm=TRUE))}
.max <- function(x){return(max(x, na.rm=TRUE))}
.avg <- function(x){return(mean(x, na.rm=TRUE))}
.std <- function(x){return(sd(x, na.rm=TRUE))}
.qtl <- function(x, probs){return(quantile(x, probs=probs, na.rm=TRUE))}
.iqr <- function(x){return(IQR(x, na.rm=TRUE))}



# ## TERGM RESULT
# results_file <- if (d==2){
#    file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s_d%s.rds',firm_i,nPeriod,R,m_x,d))
# } else {
#    file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
# } 
# fits <- readRDS(results_file)
# fit <- fits[[firm_i]][[m_x]]


## create igraph lists
# gs <- lapply(nets, asIgraph)


## visualize effects bootstrap normality
effects <- names(fit@effects)
par(mfrow=c(3,3))
for (effect in effects) {
   qqnorm(fit@boot$t[,effect], main=effect)
}
## visualize effects bootstrap distributions
par(mfrow=c(3,3))
for (effect in effects) {
   x <- fit@boot$t[,effect]
   dr <- diff(range(x))
   xl1 <- min(0-(.1*dr), min(x))
   xl2 <- max(.1*dr, max(x))
   ci <- quantile(x, c(.025,.975))
   hist(x, main=effect, breaks=13, col=rgb(.2,.2,.2,.2), xlim=c(xl1,xl2))
   abline(v=ci, col='red',lwd=2)
   abline(v=0, col='black',lwd=1, lty=2)
   segments(x0=ci[1],x1=ci[2],y0=0,y1=0,col='red',lwd=2)
}

# library()
# par(mfrow=c(1,1))
# ci95 <- c(.025,.5,.975)
# coef.mat <- t(apply(fits@boot$t,2,function(x)quantile(x,ci95)))
# matplot(x=df.coef, y=1:nrow(coef.mat))


##===============================================
## interpreation data frame i-->j (for all t)
##-----------------------------------------------

name_i <- 'microsoft'
d <- 2

nets <- geli
g <- asIgraph(nets[[length(nets)]])
vcnt <- vcount(g)
time.steps <- fit@time.steps
firm.names <-  V(g)$vertex.names
v.focal <- as.integer( V(g)[V(g)$vertex.names==name_i] )

## indices of firms that are high degree to compute their probability
mindeg <- 5
j.sub.deg.t1 <- which(igraph::degree(asIgraph(nets[[1]]), mode='all') >= mindeg)
j.sub.deg.t2 <- which(igraph::degree(asIgraph(nets[[2]]), mode='all') >= mindeg)
j.sub.deg.t3 <- which(igraph::degree(asIgraph(nets[[3]]), mode='all') >= mindeg)
j.sub.deg.t4 <- which(igraph::degree(asIgraph(nets[[4]]), mode='all') >= mindeg)
j.sub.deg <- unique(c(j.sub.deg.t1, j.sub.deg.t2, j.sub.deg.t3, j.sub.deg.t4))
j.sub.deg <- sort(j.sub.deg)
j.sub.deg <- j.sub.deg[which(j.sub.deg != v.focal)]
length(j.sub.deg)

## distance matrices
dl.in <- list()
dl.out <- list()
dl.all <- list()
for (t in 1:fit@time.steps) {
   dl.in[[t]]  <- igraph::distances(asIgraph(nets[[t+1]]),v = v.focal, mode = 'in')[1, ]
   dl.out[[t]] <- igraph::distances(asIgraph(nets[[t+1]]),v = v.focal, mode = 'out')[1, ]
   dl.all[[t]] <- igraph::distances(asIgraph(nets[[t+1]]),v = v.focal, mode = 'all')[1, ]
}


## make data frame of predicted probabilities
# for (t in 1:time.steps)
#   idf[,paste0('t',t)]<- NA
# ## time columns
# tcolnames <- names(idf)[grep(names(idf), pattern = 't\\d{1,3}', perl = T)]

## competitor probability micro-interpretation  file
interpfilebase <- sprintf('interpret_%s_%s-%s_R%s_%s',
                          name_i, startYr, endYr, R, structvar)

## main loop: competitor probability micro-interpretation
# idf <- data.frame()
# for (t in 1:fit@time.steps) {
#    cat(sprintf('t=%s, time.step=%s, period=%s\n', t, t+1, names(geli)[t+1]))
#    probs <- btergm::interpret(fit, target=fit@data$networks[[t]], t=t, type='tie', i=v.focal, j=1:10)
# }
# for (j in j.sub.deg) {  #1:vcount(g)
#    cat(sprintf('%s:  %s (%.3f%s)\n',j,V(g)$vertex.names[j], 100*which(j %in% j.sub.deg)/length(j.sub.deg),'%'))
#    j.name <- V(g)$vertex.names[j]
#    if (j == v.focal) {
#       probs <- rep(0, time.steps)
#    } else {
#       probs <- btergm::interpret(fit, type='tie', i=v.focal, j=j)
#    }
#    for (t in 1:length(probs)) {
#       tmp <- data.frame(i=v.focal, j=j, i.name=name_i, j.name=j.name, t=t+1, pd=names(geli)[t+1],
#                         d.in=dl.in[[t]][j], d.out=dl.out[[t]][j], dl.all=dvec.all[[t]][j], 
#                         p=probs[t])
#       idf <- rbind(idf, tmp)
#    }
#    if (j %% 10 == 0)  write.csv(x = idf, file = sprintf('%s.csv', interpfilebase))
# }
# write.csv(x = idf, file = sprintf('%s.csv', interpfilebase))


## main loop: competitor probability micro-interpretation

# # idf <- data.frame()
# j.sub.deg <- 1:vcount(g)
# for (j in j.sub.deg) {  #1:vcount(g)
#    j.name <- V(g)$vertex.names[j]
#    if (j %in% idf$j) {
#       cat(sprintf('skipping j=%s, %s in idf\n',j, j.name))
#       next
#    }
#    
#    cat(sprintf('%s:  %s (%.3f%s)\n',j,V(g)$vertex.names[j], 100*which(j.sub.deg %in% j)/length(j.sub.deg),'%'))
#    
#    t <- fit@time.steps
#       
#    if (j == v.focal) {
#       prob <- 0
#    } else {
#       prob <- btergm::interpret(fit, type='tie', i=v.focal, j=j, t=t)
#    }
# 
#    tmp <- data.frame(i=v.focal, j=j, i.name=name_i, j.name=j.name, t=t+1, pd=names(fit@data$networks)[t],
#                      d.in=dl.in[[t]][j], d.out=dl.out[[t]][j], d.all=dl.all[[t]][j],
#                      p=prob)
#    idf <- rbind(idf, tmp)
#    
#    if (nrow(idf) %% 5 == 0)  write.csv(x = idf, file = sprintf('%s.csv', interpfilebase), row.names = F)
# }
# write.csv(x = idf, file = sprintf('%s.csv', interpfilebase), row.names = F)

# idf <- data.frame()
# j.sub.deg <- 1:vcount(g)

for (i in 1:vcount(g)) {
   
   for (j in 1:vcount(g)) {  #1:vcount(g)
      
      i.name <- V(g)$vertex.names[i]
      j.name <- V(g)$vertex.names[j]
      if (length(which(idf$j == j & idf$i == i)) > 0) {
         cat(sprintf('skipping i=%s %s, j=%s %s, in idf\n', i,i.name, j, j.name))
         next
      }
      
      cat(sprintf('%s->%s: %s -> %s (%.3f%s : %.3f%s)\n',i,j,i.name,j.name, 100*which(1:vcount(g) %in% i)/vcount(g),'%', 100*which(1:vcount(g) %in% j)/vcount(g),'%'))
      
      t <- fit@time.steps
      
      if (j == i) {
         prob <- 0
      } else {
         prob <- btergm::interpret(fit, type='tie', i=i, j=j, t=t)
      }
      
      tmp <- data.frame(i=i, j=j, i.name=i.name, j.name=j.name, t=t+1, pd=names(fit@data$networks)[t],
                        d.in=dl.in[[t]][j], d.out=dl.out[[t]][j], d.all=dl.all[[t]][j],
                        p=prob)
      idf <- rbind(idf, tmp)
      
      if (nrow(idf) %% 5 == 0)  write.csv(x = idf, file = sprintf('%s.csv', interpfilebase), row.names = F)
   }
   
   write.csv(x = idf, file = sprintf('%s.csv', interpfilebase), row.names = F)
   
}
write.csv(x = idf, file = sprintf('%s.csv', interpfilebase), row.names = F)



idf <- read.csv( sprintf('%s.csv', interpfilebase))
idf <- idf[order(idf$j),]

# j.sub.deg <- 1:4
# t <- fit@time.steps
# ge <- asIgraph(fit@data$networks$`2014|2015`)
# vedf <- as_data_frame(ge, 'verteices')
# 
# registerDoParallel(cores = 4)
# system.time({
#    
#    foreach(j = j.sub.deg) %dopar% {
#       
#       cat(sprintf('%s:  %s (%.3f%s)\n',j,V(g)$vertex.names[j], 100*which(j.sub.deg %in% j)/length(j.sub.deg),'%'))
#       j.name <- V(g)$vertex.names[j]
#       
#       
#       if (j == v.focal) {
#          prob <- 0
#       } else {
#          prob <- btergm::interpret(fit, type='tie', i=v.focal, j=j, t=t)
#       }
#       
#       tmp <- data.frame(i=v.focal, j=j, i.name=name_i, j.name=j.name, t=t+1, pd=names(fit@data$networks)[t],
#                         d.in=dl.in[[t]][j], d.out=dl.out[[t]][j], d.all=dl.all[[t]][j],
#                         p=prob)
#       # idf <- rbind(idf, tmp)
#       
#       # if (nrow(idf) %% 5 == 0)  write.csv(x = idf, file = sprintf('%s.csv', interpfilebase))
#       
#    }
#    
# })




qmap <- function(x, probs=seq(0,1, 1/3)) {
   probs <- sort(probs)
   qs <- quantile(x, probs)
   y <- c()
   for (i in 1:length(x)) {
      for (j in 2:length(qs)) {
         # cat(sprintf('i=%s, j=%s\n',i,j))
         if (x[i] <= qs[j]) {
            y[i] <- j-1
            break
         }
      }
   }
   return(y)
}

##-------------------------
##  PLOT NETWORK COLOR CLUSTER
##_------------------------
seed <- 5215215
t <- 4
par(mar=c(.1,.1,.1,.1))
.gx <- asIgraph(geli[[t]])
gxx <- induced_subgraph(.gx, which(igraph::degree(.gx)>0))
gcc <- induced.subgraph(asIgraph(gcli[[t]]), which(gcli[[t]] %v% 'vertex.names' %in% V(gxx)$vertex.names))
gxx.v.focal <- which(V(gxx)$vertex.names == 'microsoft')
#
V(gxx)$vertex.color <- NA
# # # LEAVE DISCONNECTED COMPONENT GRAY - COLOR CLUSTERS
# decomp <- decompose(gxx,mode = 'weak',min.vertices = 3)
# lcc <- decomp[[which.max(sapply(decomp,vcount))]]
# lccidx <- which(V(gxx)$vertex.names %in% V(lcc)$vertex.names)
# lcccomm <- igraph::multilevel.community(as.undirected(lcc))$membership
# #
# lcccomm <- igraph::edge.betweenness.community(lcc)$membership
# nclust <- length(unique(lcccomm))
# V(gxx)$cluster <- 0
# V(gxx)$cluster[lccidx] <- lcccomm
# V(gxx)$vertex.color[lccidx] <- rainbow(nclust,alpha = .6)[lcccomm]
vid.color <- which( V(gxx)$vertex.names %in% idf$j.name )
heatcols <-  heat.colors(3, rev=T, alpha = .8)
qs <- qmap(idf$p[which(idf$j.name %in% V(gxx)$vertex.names[vid.color])], seq(0,1, 1/3))
V(gxx)$vertex.color[vid.color] <- heatcols[qs]
#
V(gxx)$vertex.shape <- 'circle'
V(gxx)$vertex.shape[gxx.v.focal] <- 'square'
V(gxx)$vertex.color[gxx.v.focal] <- 'blue'
#
idx.comp <- igraph::neighbors(gcc, v = which(V(gxx)$vertex.names==name_i))
idx.comp.in.fmnet <- which(V(gxx)$vertex.names %in% V(gxx)$vertex.names[idx.comp])
V(gxx)$vertex.shape[idx.comp.in.fmnet] <- 'square'
#
V(gxx)$vertex.size <- 3
V(gxx)$vertex.size[c(gxx.v.focal,idx.comp.in.fmnet)] <- 5
#
set.seed(seed)
plot.igraph(gxx, layout=layout.fruchterman.reingold(gxx, grid="nogrid"),
            vertex.color=V(gxx)$vertex.color,
            vertex.shape=V(gxx)$vertex.shape,
            vertex.size=V(gxx)$vertex.size,
            vertex.label=NA,
            vertex.label.cex = 0,
            edge.width=.2,
            edge.arrow.size=.2)
legend('topright', title='FM Threats:',
       legend=c('FM (>4 ties)','Low','Mid','High','PM Rivals','Low','Mid','High'), 
       col=c(NA,heatcols), pch=c(rep(16,4),rep(15,4)))
###----------------------------------------------------------
gcc.v.focal <- which(V(gcc)$vertex.names==name_i)
V(gcc)$vertex.shape <- 'circle'
V(gcc)$vertex.shape[gcc.v.focal] <- 'square'
V(gcc)$vertex.shape[idx.comp] <- 'square'
V(gcc)$vertex.color <- NA
V(gcc)$vertex.color[gcc.v.focal] <- 'blue'
V(gcc)$vertex.color[idx.comp] <- 'darkred'
V(gcc)$vertex.size <- 3
V(gcc)$vertex.size[gcc.v.focal] <- 5
V(gcc)$vertex.size[idx.comp.in.fmnet] <- 4
#
gccsub <- induced.subgraph(gcc, which(igraph::degree(gcc)>0))
set.seed(seed+3)
plot.igraph(gccsub, layout=layout.fruchterman.reingold(gccsub, grid="nogrid"),
            vertex.color=V(gccsub)$vertex.color,
            vertex.shape=V(gccsub)$vertex.shape,
            vertex.size=V(gccsub)$vertex.size,
            vertex.label=NA,
            vertex.label.cex = 0,
            edge.width=.2,
            edge.arrow.size=.2, 
            set.seed=seed)
##--------------------------------------





# ##-------------------------
# ##  PLOT NETWORK COLOR CLUSTER
# ##_------------------------
# par(mar=c(.1,.1,.1,.1))
# .gx <- asIgraph(geli[[2]])
# gxx <- induced_subgraph(.gx, which(igraph::degree(.gx)>0))
# #
# V(gxx)$vertex.color <- 'gray'
# decomp <- decompose(gxx,mode = 'weak',min.vertices = 3)
# lcc <- decomp[[which.max(sapply(decomp,vcount))]]
# lccidx <- which(V(gxx)$vertex.names %in% V(lcc)$vertex.names)
# # lcccomm <- igraph::multilevel.community(as.undirected(lcc))$membership
# lcccomm <- igraph::edge.betweenness.community(lcc)$membership
# nclust <- length(unique(lcccomm))
# V(gxx)$cluster <- 0
# V(gxx)$cluster[lccidx] <- lcccomm
# V(gxx)$vertex.color[lccidx] <- rainbow(nclust,alpha = .6)[lcccomm]
# #
# plot.igraph(gxx, layout=layout.fruchterman.reingold(gxx, grid="nogrid"),
#             vertex.color=V(gxx)$vertex.color,
#             vertex.size=2.3,
#             vertex.label=NA,
#             vertex.label.cex = 0,
#             edge.width=.5,
#             edge.arrow.size=.3)
# ##--------------------------------------





# glmfit <- glm(fit@response ~ . -1, data=data.frame(fit@effects), weights = fit@weights, family="binomial")
#
# data(api)
# dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw,
#           data = apistrat, fpc = ~fpc)
# glmfit <- svyglm(api00 ~ avg.ed * growth, design = dstrat)

# data <- cbind(data.frame(response=fit@response),fit@effects)
# interact_plot(glmfit, pred='edgecov.pmuncommon[[i]]', modx='edgecov.a.pmu.status[[i]]', 
#               data = data)


y <- fit@response
X <- fit@effects
X.med <- apply(X,2,function(x) median(x, na.rm = T))


prednm <- 'edgecov.pmuncommon[[i]]'



# q.pred <- quantile(X[,prednm], probs = c(.16, .5, .84))
pred.1 <- 1
pred.0 <- 0

modnames <- c('Status', 'Similarity', 'Closure')
mods <- c('edgecov.status[[i]]', 'edgecov.cossim[[i]]', 'edgecov.cycle3[[i]]')
modints <- c('edgecov.a.pmu.status[[i]]', 'edgecov.a.pmu.cossim[[i]]', 'edgecov.a.pmu.cycle3[[i]]')
pred <- c(0,1)
qs <- c(.16, .5, .84)
qnames <- c('Low', 'Mid', 'High')
prdf <- data.frame()

for (k in 1:3) {  # moderator hypotheses
   
   if (modnames[k] == 'Closure') {
      q.mod <- c(0,1,2)
   } else {
      q.mod <- quantile(X[,mods[k]], probs = c(.16, .5, .84))
   }
      
   for (j in 1:3) { # quantiles low,mid,high   
         
      for (i in 1:2) { # predictor (PMU=0,1)

         X.med.q <- X.med 
         X.med.q[ which(names(X.med.q)==prednm) ] <- pred[i]
         X.med.q[ which(names(X.med.q)==mods[k]) ] <- q.mod[j]
         X.med.q[ which(names(X.med.q)==modints[k]) ] <- q.mod[j] * pred[i]
         lp <- t(X.med.q) %*% cbind(fit@coef)
         p.hat <- 1/(1+exp(-lp))
         #
         tmp <- data.frame(X=pred[i], q=qs[j], level=qnames[j], 
                           Mname=modnames[k], M=q.mod[j], 
                           X.M=q.mod[j] * pred[i], 
                           lp=lp, p=p.hat)
         prdf <- rbind(prdf, tmp)
      }
      
   }
   
}


ggplot(aes(x=X, y=p,  colour=level, fill=level), data=prdf) + geom_line() + geom_point() +
   facet_wrap(.~Mname) + theme_bw()




mindeg <- 3
nets <- geli
j.sub.deg.t1 <- which(igraph::degree(asIgraph(nets[[1]]), mode='all') >= mindeg)
j.sub.deg.t2 <- which(igraph::degree(asIgraph(nets[[2]]), mode='all') >= mindeg)
j.sub.deg.t3 <- which(igraph::degree(asIgraph(nets[[3]]), mode='all') >= mindeg)
j.sub.deg.t4 <- which(igraph::degree(asIgraph(nets[[4]]), mode='all') >= mindeg)
j.sub.deg <- unique(c(j.sub.deg.t1, j.sub.deg.t2, j.sub.deg.t3, j.sub.deg.t4))
j.sub.deg <- sort(j.sub.deg)
j.sub.deg <- j.sub.deg[which(j.sub.deg != v.focal)]
length(j.sub.deg)


# geli ~ edges + edgecov(cycle3) + edgecov(ttriple) + memory(type = "stability", lag = 1) + delrecip(mutuality = TRUE, lag = 1) + edgecov(tofrexecdiff) + 
#    nodecov("age") + nodecov("ipo_status") + edgecov(status) + 
#    edgecov(agediff) + edgecov(mmc) + nodematch("state_code", diff = F) + edgecov(cossim) + edgecov(pmuncommon) + edgecov(a.pmu.status) + 
#    edgecov(a.pmu.cossim) + edgecov(a.pmu.cycle3) + edgecov(a.pmu.ttriple)

j.sub.deg <- sort(unique(c(j.sub.deg, v.focal)))
t <- 4  ## period network
X <- data.frame(
   1,
   c(cycle3[[t]][j.sub.deg,j.sub.deg]),
   c(ttriple[[t]][j.sub.deg,j.sub.deg]),
   c(ttriple[[t]][j.sub.deg,j.sub.deg]),
)




y <- fit@response
X <- fit@effects[which(fit@weights < 2),]
X.med <- apply(X,2,function(x) median(x, na.rm = T))


prednm <- 'edgecov.pmuncommon[[i]]'



# q.pred <- quantile(X[,prednm], probs = c(.16, .5, .84))
pred.1 <- 1
pred.0 <- 0


modnames <- c('Status', 'Similarity', 'Closure')
mods <- c('edgecov.status[[i]]', 'edgecov.cossim[[i]]', 'edgecov.cycle3[[i]]')
modints <- c('edgecov.a.pmu.status[[i]]', 'edgecov.a.pmu.cossim[[i]]', 'edgecov.a.pmu.cycle3[[i]]')
pred <- c(0,1)
qs <- c(.16, .84)
qnames <- c('Low', 'High')
prdf <- data.frame()

for (k in 1:3) {  # moderator hypotheses
   
   if (modnames[k] == 'Closure') {
      q.mod <- c(0,1)
   } else {
      q.mod <- quantile(X[,mods[k]], probs = qs)
   }
   
   for (j in 1:2) { # quantiles low,mid,high   
      
      for (i in 1:2) { # predictor (PMU=0,1)
         
         X.med.q <- X.med 
         X.med.q[ which(names(X.med.q)==prednm) ] <- pred[i]
         X.med.q[ which(names(X.med.q)==mods[k]) ] <- q.mod[j]
         X.med.q[ which(names(X.med.q)==modints[k]) ] <- q.mod[j] * pred[i]
         lp <- t(X.med.q) %*% cbind(fit@coef)
         p.hat <- 1/(1+exp(-lp))
         #
         tmp <- data.frame(X=pred[i], level=qnames[j], 
                           Mname=modnames[k], M=q.mod[j], 
                           X.M=q.mod[j] * pred[i], 
                           lp=lp, p=p.hat)
         prdf <- rbind(prdf, tmp)
      }
      
   }
   
}


ggplot(aes(x=X, y=p,  colour=level, fill=level), data=prdf) + geom_line() + geom_point() +
   facet_wrap(.~Mname) + theme_bw()






# ## main loop: competitor probability micro-interpretation
# if (file.exists(sprintf('%s.csv', interpfilebase))) {
#   ## READ IN PREDICTED PROBS
#   idf <- read.csv(sprintf('%s.csv', interpfilebase))
# } else {
#   idf <- data.frame()
#   for (j in j.sub.deg) {  #1:vcount(g)
#     cat(sprintf('%s:  %s (%.3f%s)\n',j,V(g)$vertex.names[j], 100*which(j %in% j.sub.deg)/length(j.sub.deg),'%'))
#     j.name <- V(g)$vertex.names[j]
#     if (j == v.focal) {
#       probs <- rep(0, time.steps)
#     } else {
#       probs <- btergm::interpret(fit, type='tie', i=v.focal, j=j)
#     }
#     for (t in 1:length(probs)) {
#       tmp <- data.frame(i=v.focal, j=j, i.name=name_i, j.name=j.name, t=t+1, pd=names(geli)[t+1],
#                         d.in=dl.in[[t]][j], d.out=dl.out[[t]][j], dl.all=dvec.all[[t]][j],
#                         p=probs[t])
#       idf <- rbind(idf, tmp)
#     }
#     if (j %% 5 == 0)  write.csv(x = idf, file = sprintf('%s.csv', interpfilebase))
#   }
#   write.csv(x = idf, file = sprintf('%s.csv', interpfilebase))
# }

###
## READ IN Interpreataion Predicted Probabilities
##
idf <- read.csv(sprintf('%s.csv', interpfilebase))



## add covariates
idf$genidx_multilevel <- NA
idf$njobs_multilevel <- NA
idf$cent_pow_n0_4 <- NA
idf$absdiff_pow_n0_4 <- NA
idf$year <- NA
for (row in 1:nrow(idf)) {
   t <- idf$t[row]
   i <- idf$i[row]
   j <- idf$j[row]
   # idf$d[row] <- igraph::distances(gs[[t+1]], v = i, to = j)[1,1]
   idf$year[row] <- 2006 + as.integer(t)
   idf$genidx_multilevel[row] <- (nets[[t+1]] %v% 'genidx_multilevel')[j]
   idf$njobs_multilevel[row] <- (nets[[t+1]] %v% 'njobs_multilevel')[j]
   idf$cent_pow_n0_4[row] <- (nets[[t+1]] %v% 'cent_pow_n0_4')[j]
   idf$absdiff_pow_n0_4[row] <- abs((nets[[t+1]] %v% 'cent_pow_n0_4')[j] - (nets[[t+1]] %v% 'cent_pow_n0_4')[i])
   if(row %% 300 == 0) cat(sprintf('row %s\n',row))
}

## remove probabilities and covariates (set NA) for Inf distance (firm's not in component that pd)
idf$p[idf$d==Inf] <- NA
idf$genidx_multilevel[idf$d==Inf] <- NA
idf$njobs_multilevel[idf$d==Inf] <- NA
idf$cent_pow_n0_4[idf$d==Inf] <- NA
idf$absdiff_pow_n0_4[idf$d==Inf] <- NA

## distance category
idf$d_cat <- as.factor(sapply(idf$d,function(x){
   xstr <- as.character(x)
   ifelse(xstr %in% c('1','2','3','4'), xstr, ifelse(x==Inf, NA, '5+'))
}))
idf$i <- as.factor(idf$i)
idf$year <- as.factor(idf$year)

## distance category
idf$d_cat <- as.factor(sapply(idf$d,function(x){
   xstr <- as.character(x)
   ifelse(xstr %in% c('1','2','3','4'), xstr, ifelse(x==Inf, NA, '5+'))
}))
idf$i <- as.factor(idf$i)
idf$year <- as.factor(idf$year)

## aware.all.cutoff
aware.frac <- .28
xp <- idf$p[idf$year=='2016' & !is.na(idf$p)]
aware.all.cutoff <- .qtl(xp, 1 - aware.frac)

## mutate group-period factors
idf2 <-
   idf[idf$d != 0 & idf$d !='0', ] %>%
   group_by(d) %>%
   mutate(
      outlier = p > .med(p) + .iqr(p) * 2, 
      low.otl = p < .med(p) - .iqr(p) * 2,
      aware.med = p > .qtl(p, .5),
      aware.cutoff = p > aware.all.cutoff
   ) %>% 
   ungroup

## manual colors
colors2rw <- c('black','red')
colors2 <- c( "#333333", '#aaaaaa')
colors4 <- c( '#222222', '#aaaaaa', '#555555', '#888888')
colors5 <- c( '#111111', '#333333','#555555', '#777777','#999999')







##---------------- 
structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                 ##---- STRUCTURAL ----------------------------
              edgecov(cycle3) + edgecov(ttriple) +
                 ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                 ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                 ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") +  edgecov(status) +
                 ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) + edgecov(mmc) +  nodematch("state_code", diff = F) + ## GEOGRAPHIC
                 ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
              edgecov(cossim) +
                 ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(pmuncommon) + 
                 ##---- INTERACTIONS -----------------
              edgecov(a.pmu.status) + 
                 edgecov(a.pmu.cossim) + 
                 edgecov(a.pmu.cycle3) + edgecov(a.pmu.ttriple) 
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()























# categories
cat <- m5a@data$networks$`2014|2015` %v% 'category_list'

cz <- c(unname(unlist(sapply(cat, function(x) {
   strsplit(x, split = '[|]')[[1]]
}))))

cnt <- plyr::count(cz)
cnt <- cnt[order(cnt$freq, decreasing=T), ]










