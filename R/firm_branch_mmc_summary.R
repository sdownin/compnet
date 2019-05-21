
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(igraph)
library(intergraph)
library(texreg)
library(plyr)

## DIRECTORIES
data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')
net_dir <- file.path(work_dir,'firm_nets_rnr2','firmctrl')
sup_data_dir <- file.path(work_dir,'amj_rnr2_sup_data')  ## supplmental data dir

## set woring dir
setwd(work_dir)

## load awareness functions
aaf <- source(file.path(version_dir,'amj_awareness_functions.R'))$value
cb  <- source(file.path(version_dir,'amj_cb_data_prep.R'))$value           ## cb : CrunchBase
g.full <- source(file.path(version_dir,'amj_make_full_graph.R'))$value        ## g.full : full competition graph


br <- cb$co_br[which(cb$co_br$company_name_unique %in% V(g.full)$name), ]

firms <- V(g.full)$name[which(V(g.full)$name %in% br$company_name_unique)]

fidx <- which(cb$co$company_name_unique %in% firms)

head(cb$fund)
head(cb$co)
head(cb$co_rou)
head(cb$inv_rou)
head(cb$inv[which(cb$inv$investor_type=='corporate_venture_capital'),])

plyr::count(cb$inv$investor_type)

vc <- cb$inv[which(cb$inv$investor_type=='venture_capital'),]
vccnt <- plyr::count(vc$investor_name_unique)
vccnt <- vccnt[order(vccnt$freq, decreasing = T), ]
vccnt

ir <- cb$inv_rou[which(cb$inv_rou$investor_uuid %in% vc$uuid),]
plyr::count(ir$investor_uuid)

results_dir <- 'C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2/amj_rnr2_results'


# ('clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')

cat(sprintf('\n\n checking %s\n', name_i))
file_i <- file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',name_i,nPeriod,R,mx))



patch.todo <- c('qualtrics')

for (name_i in patch.todo) {
  nPeriod <- 11
  mx <- 'm4_7cycle'
  R <- 2000
  cat(sprintf('\n\n checking %s\n', name_i))
  file_i <- file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',name_i,nPeriod,R,mx))
  fits <- readRDS(file_i)
  fit <- fits[[name_i]][[mx]]
  rm(fits); gc()
  screenreg(fit, digits = 4, ci.force = T, ci.force.level = .5)
}



## CHECK GOF
nets <- readRDS('')














##-----------------------------------------
library("network")
set.seed(5)

networks <- list()
for(i in 1:10){            # create 10 random networks with 10 actors
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0           # loops are excluded
  nw <- network(mat)       # create network object
  networks[[i]] <- nw      # add network to the list
}

covariates <- list()
for (i in 1:10) {          # create 10 matrices as covariate
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  covariates[[i]] <- mat   # add matrix to the list
}

fit <- btergm(networks ~ edges + istar(2) +
                edgecov(covariates), R = 100)

summary(fit)               # show estimation results

# }
# NOT RUN {
# The same example using MCMC MLE:

fit2 <- mtergm(networks ~ edges + istar(2) + 
                 edgecov(covariates))

summary(fit2)


gf <- gof(fit, statistics=c(deg,esp,dsp,geodesic), nsim=10)

summary(gf)

saveRDS(gf, file = 'test_save_gof_rds.rds')
gf <- readRDS('test_save_gof_rds.rds')


 
gfc <- readRDS(file.path(results_dir,'firmctrl','gof_qualtrics_pd11_R2000_m4_nsim5_geodesic.rds'))















##--------------------------
## DEGREE
##--------------------------
## OBSERVED MAX
x <- c(
  386,149,86,73,40,29,18,17,12,10,5,10,6,6,7,6,
  2,4,3,2,2,1,3,2,2,1,1,1,1,2,1,1,1,1,1,1,1,2,1,
  1,1,1,1,1,1,1)
xdf <- data.frame(i=0:(length(x)-1),x=x,cumu=cumsum(x))
xdf$p <- xdf$cumu / sum(xdf$x)

## OBSERVED MEDIAN
x <- c(
  193, 89.5, 63, 60.5, 29.5, 24.5, 15.5, 
  8, 7.5, 6.5, 6.5, 2.5, 5.0, 3.0, 3.5, 
  2.5, 1.5, 3.5,
  1, 1.5, 1, 1.5, 0, .5, 0, 1, .5)
xdf <- data.frame(i=0:(length(x)-1),x=x,cumu=cumsum(x))
xdf$p <- xdf$cumu / sum(xdf$x)




