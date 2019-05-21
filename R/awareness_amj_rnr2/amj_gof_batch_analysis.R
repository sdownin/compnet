##
##
##  Goodness of Fit 
## 
##  Batch Simulations List Analysis
##
##

library(btergm)
library(MASS)
library(plyr)
library(dplyr)


work_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2'
data_dir <- file.path(work_dir, 'firm_nets_rnr2\\firmctrl')
results_dir <- file.path(work_dir,'amj_rnr2_results\\firmctrl')

##====================================
## Functions
##------------------------------------


##
# cbind list of matrices intro one matrix
##
lcbind <- function(l){
  x <- l[[1]]
  for (i in 2:length(l)) {
    x <- cbind(x, l[[i]])
  }
  return(x)
}

###
## Goodness of Fit for batch computed simluations
## NOTE:  Nets and covarites used to create simluations 
##        must be cached in memory
##
## @param  l [list]   The simulated 
## @param  st [character] name of function
## @param  alpha [numeric] The alpha level for obs-vs-sim group means t test
## @param  na.rm [bool]  Flag to return GOF DF without statistic values with NaN|NA pvals
## @param  plot [bool]  Flag to create plot
## @param  plotname [character] The plot file name (or full path)
##
## @return df [dataframe]
###
gof.batch <- function(l, st,  alpha=0.05, na.rm=TRUE,
                      plot=TRUE, plotname='gof_batch.png')
{
  
  stats <- list(deg=deg, dsp=dsp, esp=esp, geodesic=geodesic)
  stnames <- c(
    dsp='Dyad-wise Shared Partners',
    esp='Edge-wise Shared Partners',
    deg='Degree',
    geodesic='Geodesic Distance'
  )  
  
  ## AUXILIARY STATISTIC
  stfunc <- stats[[st]]
  
  cat(sprintf('GOF statistic: %s\n simulations\n', st))
  
  ## combine all networks over simulations & batches to one matrix 
  simmat <- lcbind(l[[st]][[1]])
  for (i in 2:length(l[[st]])) {
    simmat <- cbind(simmat, lcbind(l[[st]][[i]]))
  }
  
  cat(' observations\n')
  
  ## observed networks stats
  obs <- sapply(nets[2:length(nets)], function(net) stfunc(as.matrix(net[,])))
  
  cat(' summary data.frame\n')
  
  ## main GOF DF
  df <- data.frame(
    `o.min`=apply(obs, 1, min),
    `o.q1`=apply(obs, 1, function(x)quantile(x, probs = .25)),
    `o.mean`=apply(obs, 1, mean),
    `o.med`=apply(obs, 1, median),
    `o.q3`=apply(obs, 1,  function(x)quantile(x, probs = .75)),
    `o.max`=apply(obs, 1, max),
    #
    `s.min`=apply(simmat, 1, min),
    `s.q1`=apply(simmat, 1, function(x)quantile(x, probs = .25)),
    `s.mean`=apply(simmat, 1, mean),
    `s.med`=apply(simmat, 1, median),
    `s.q3`=apply(simmat, 1,  function(x)quantile(x, probs = .75)),
    `s.max`=apply(simmat, 1, max)
  )
  
  cat(' pvals\n')
  
  ## test gorup means
  df$`Pr(>z)` <- NA
  df$`.` <- NA
  for (i in 1:nrow(df)) {
    p <- t.test(obs[i, ], simmat[i, ])$p.value
    df$`Pr(>z)`[i] <- p
    df$`.`[i] <- ifelse(p < alpha, '* ', '  ')
  }
  
  df2 <- df[which(!is.nan(df$`Pr(>z)`) & !is.na(df$`Pr(>z)`)),]
  
  if (plot) {  
    cat(' plotting\n')
    
    #lower and upper outlier cutoffs (whiskerys of box plot)
    dfp2 <- df2
    
    ## replace values with relative frequences
    o.cols <- grep('^o\\.', names(dfp2)) 
    s.cols <- grep('^s\\.', names(dfp2)) 
    dfp2[,o.cols] <- dfp2[,o.cols] / max(dfp2$o.max)
    dfp2[,s.cols] <- dfp2[,s.cols] / max(dfp2$s.max)
    
    # Outlier range (whiskers in box plot)
    dfp2$s.ol <- apply(dfp2[,c('s.q1','s.q3')], 1, function(x) {
      return(x[1] - (1.5 * (x[2] - x[2])))
    })
    dfp2$s.ou <- apply(dfp2[,c('s.q1','s.q3')], 1, function(x) {
      return(x[2] + (1.5 * (x[2] - x[2])))
    })
    
    cols <- which(names(dfp2) %in% c('s.q1','s.med','s.q3'))
    dfp2t <- t(dfp2[,cols])
    
    ## save plot
    png(filename = plotname, height = 5, width = 6, units = 'in', res=200)
      par(mar=c(4.5,4.1,1.8,.8))
      boxplot(dfp2t, ylab='Frequency', xlab=stnames[st], main=stnames[st])
      lines(dfp2$s.mean, lty=2)
      lines(dfp2$s.ou, lty=3)
      lines(dfp2$s.ol, lty=3)
      lines(dfp2$o.mean, lty=1, lwd=2)
    dev.off()
  }
  
  cat(sprintf('%s simulations done.\n\n', ncol(simmat)))
  
  if (na.rm==TRUE) {
    return(df2)
  } else {
    return(df)
  }
  
}











# ##============================================================
# #
# #  Modified GOF from btergm for batch processed network statistics
# #
# ##------------------------------------------------------------
# 
## SETTINGS
firm_i <- 'qualtrics'
d <- 3
R <- 2000
m_x <- 'm4'
nsim <- 10
nPeriod <- 11


## NETWORKS LIST
data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
nets <- readRDS(data_file)
# 
# ## make MMC nets list
# mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
# cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
# cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
# # cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
# cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
# # centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
# # centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
# shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor'))
# shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))
# 
# ## TERGM RESULT
# results_file <- if (d==2){
#   file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s_d%s.rds',firm_i,nPeriod,R,m_x,d))
# } else {
#   file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
# } 
# fits <- readRDS(results_file)
# fit <- fits[[firm_i]][[m_x]]

##--------------------------------------------------
##  Load Batch GOF Stats List
##--------------------------------------------------
gffile <- 'gof_batch_sim_stats_qualtrics_pd11_R2000_m4.rds'
## load
l <-  readRDS(file.path(results_dir, gffile))
## remove last incomplete batch
nbatch <- length(l[[1]])
for (i in 1:length(l)) {
  l[[i]] <- l[[i]][1:(nbatch-1)]
}

##--------------------------------------------------
## Run GOF for all stratics 
## output plots
##--------------------------------------------------
gfl <- list()
for (st in c('deg','dsp','esp','geodesic')) {
  plotfilename <- sprintf('gof_batch_plot_%s_pd%s_R%s_%s_%s.png', firm_i, nPeriod, R, m_x, st)
  plotfile <- file.path(results_dir, plotfilename)
  gfl[[st]] <- gof.batch(l=l, st=st, plotname = plotfile)
}


gfl$deg
gfl$dsp
gfl$esp
gfl$geodesic



## END





















filename <- sprintf('gof_%s_pd%s_R%s_%s_nsim%s_all.rds', firm_i, nPeriod, R, m_x, nsim)
gf <- readRDS(file.path(results_dir, filename))

summary(gf)




pattern <- sprintf('gof_%s_pd%s_R%s_%s_nsim%s_all_[0-9]+\\.rds', firm_i, nPeriod, R, m_x, nsim)
files <- dir(results_dir, pattern = pattern)

for (filename in files) {
  
  num <- gsub('(gof_qualtrics_pd11_R2000_m4_nsim10_all_|\\.rds)','',filename)
  cat(sprintf('num = %s', num))
  
  gf <- readRDS(file.path(results_dir, filename))
  
  print(gf)
  
  plotfile <- sprintf('gof_PLOT_%s_pd%s_R%s_%s_nsim%s_all_%s.png', firm_i, nPeriod, R, m_x, nsim, num)
  png(file.path(results_dir,plotfile), width = 8, height = 7, units = 'in', res= 200)
  plot(gf)
  dev.off()
}

###
### Best results 20190325020152
###

gf


###---------------------------------------------------
##  Mahalanobis Distance Test
###---------------------------------------------------
gffile <- 'gof_batch_sim_stats_qualtrics_pd11_R2000_m4.rds'

l <-  readRDS(file.path(results_dir, gffile))

nbatch <- length(l[[1]])
## remove last incomplete batch
for (i in 1:length(l)) {
  l[[i]] <- l[[i]][1:(nbatch-1)]
}
summary(l$deg)
## simulation 1
summary(l$deg[[1]])
## simulation i period t  MATRIX
l$deg[[1]][[1]]

stats <- list(deg=deg, dsp=dsp, esp=esp, geodesic=geodesic)

##limit periods
pdsub <- 1:length(l$deg[[1]])


## AUXILIARY STATISTIC
st <- 'geodesic'
stfunc <- stats[[st]]

## observed
obs.pd <- sapply(nets[2:length(nets)], function(net) stfunc(as.matrix(net[,])))
if (length(pdsub)==1) {
  obs <- obs.pd[,pdsub]
} else {
  obs <- rowMeans(obs.pd[,pdsub])
}

## reduce each simulation:  Average of periods 
lsim <- lapply(l[[st]], function(sim.pds) Reduce('+', sim.pds[pdsub])/length(pdsub))

## combine simulations over batches to one matrix by column binding
simmat <- lsim[[1]]
for (i in 2:length(lsim)) {
  simmat <- cbind(simmat, lsim[[i]])
}
dim(simmat)

##  Mu
mu <- rowMeans(simmat)

## fix length mismatch (if obs has Inf term)
if (length(obs) > length(mu)) {
  obs <- obs[1:(length(obs)-1)]
}
# cat(sprintf('length obs %s, mu %s', length(obs), length(mu)))



## subset aux stat range
# idx <- 1:min(length(obs), 80)
# idx <- 1:length(obs)
idx0 <- which(obs > 0 & rownames(simmat) != 'Inf')
idx <- seq(min(idx0), max(idx0), 1)

##
sigma <- cov( t(simmat[idx,idx]) )
sigma.inv <- ginv(sigma)
d.obs <- mahalanobis(obs[idx], mu[idx], sigma.inv, inverted = TRUE)

## for each simulation
test <- list(bin=c(), p=c())
for (i in 1:ncol(simmat)) {
  sim <- as.numeric( simmat[idx,i] )
  d.sim <- mahalanobis(sim[idx], mu[idx], sigma.inv, inverted = TRUE)
  bin <- as.integer( d.sim >= d.obs )
  test$d.sim[i] <- d.sim
  test$bin[i] <- bin
  test$p[i] <- sum(test$bin) / length(test$bin)      
}

cat(sprintf('%s: nsim=%s; obs MMD = %0.1f, sim MMD avg %.1f (%.1f,%.1f), p = %0.3f;\n  idx %s', 
            st, ncol(simmat), d.obs, mean(test$d.sim), min(test$d.sim), max(test$d.sim),
              test$p[length(test$p)], paste(rownames(simmat)[idx], collapse = ',')))

plot(test$p, type='l', ylim=c(0,1.01*max(test$p)), main=st); abline(h=.05)

