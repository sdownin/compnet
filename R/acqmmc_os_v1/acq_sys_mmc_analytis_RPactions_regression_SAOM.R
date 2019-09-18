##--------------------------------------------------------------
##
##  MMC & ACQUISITIONS 
##
##--------------------------------------------------------------
# .libPaths('C:/Users/steph/Documents/R/win-library/3.2')
library(igraph)
library(intergraph)
library(pglm)
library(lme4)
library(texreg)
library(readxl)
library(haven)
library(stringr)
library(stringdist)

## DIRECTORIES
data_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
version_dir <- file.path(work_dir,'R','acqmmc_os_v1')
net_dir <- file.path(work_dir,'acqmmc_os_v1','data')
result_dir <- file.path(work_dir,'acqmmc_os_v1','data')
sup_data_dir <- file.path(work_dir,'acqmmc_os_v1','sup_data')  ## supplmental data dir

## set woring dir
setwd(work_dir)


###
## Complexity measure of actions 
##  @see Yu, Subramaniam, & Cannella, 2009
###
complexity <- function(x, scale=FALSE) {
  sumx <- sum(x)
  sq.prop <- sapply(x, function(k) (k/sumx)^2 )
  sum.sq.prop <- sum(sq.prop)
  y <- ifelse(sum.sq.prop==0, 0, 1/sum.sq.prop)
  if (scale)
    return(y / length(x))
  return(y)
}



###
## PGLM FUNCTION FOR TEXREG TABLE
###
extract.pglm <- function (model, include.nobs = TRUE, include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$estimate)
  coefficients <- s$estimate[, 1]
  standard.errors <- s$estimate[, 2]
  significance <- s$estimate[, 4]
  loglik.value <- s$loglik
  n <- nrow(model$model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, loglik.value)
    gof.names <- c(gof.names, "Log-Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                     se = standard.errors, pvalues = significance, gof.names = gof.names, 
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)
}

## In order for this code to work, you should also register the function so that it handles the pglm maxLik objects by default when extract is called:
setMethod("extract", signature = className("maxLik", "maxLik"), 
          definition = extract.pglm)
# ## set firms to create networks (focal firm or replication study focal firms)

# head(cnt, 20)
# # Top Acquirers:             x   freq
# # 545                   google   88
# # 599                      ibm   58
# # 815                microsoft   47
# # 1496                   yahoo   46
# # 77                       aol   43
# # 944                   oracle   39
# # 64                    amazon   34
# # 88                     apple   32
# # 453                 facebook   31
# # 403                     ebay   28
# # 259                    cisco   27
# # 1109              salesforce   26
# # 1362                 twitter   24
# # 578          hewlett-packard   21
# # 414                      emc   19
# # 640                    intel   18
# # 556                  groupon   17
# # 421  endurance-international   16
# # 1250                symantec   15
# # 1347             tripadvisor   15






# firms.todo <- c('ibm','qualtrics','medallia','first-mile-geo','verint')
firm_i <- 'microsoft'
# firms.todo <- c('microsoft')
# for (firm_i in firms.todo)
# {
  # firm_i <- firms.todo[1]
  ## -- settings --
  # firm_i <- firm_i
  name_i <- firm_i
  firm_i_ego <- firm_i
  d <- 3
  yrpd <- 1
  startYr <- 2006
  endYr <- 2016            ## dropping first for memory term; actual dates 2007-2016
  ## --------------  
  periods <- seq(startYr,endYr,yrpd)
  
  
  
  ## READ IN DATA FILE
  
  y1 <- 2006
  y2 <- 2016
  
  ## load DATA
  firmfile <- sprintf('acq_sys_mmc_panel_RPactions_regression_df_ego-%s_d%s_%s-%s.csv', 
                      firm_i_ego, d, y1, y2)
  dfl <- read.csv(file.path(result_dir,firmfile), stringsAsFactors = F)
  # dfl$name <- as.character(dfl$name)
  

  ys <- y1:y2
  nl <- list()
  bl <- list()
  for (t in 1:10) 
  {
    netpdfile <- sprintf('acq_sys_mmc_panel_RPactions_NETLIST_ego-%s_d%s_%s-%s_t%s.rds', 
                         firm_i_ego, d, y1, y2, t)
    nl[[t]] <- readRDS(file.path(result_dir, netpdfile))
    # gt <- nl[[t]]$gt
    # #
    # membership <- V(gt)$com_multilevel
    # markets <- unique(membership)
    # markets <- markets[order(markets)]
    # adj <- igraph::as_adjacency_matrix(gt, sparse = F)
    # df.ms <- ldply(1:nrow(adj), function(i){
    #   i.nbr <- unname(which( adj[i, ] == 1 ))  ## rivals (neighbors in the network) ## i.nbr <- as.integer(igraph::neighbors(g, V(g)[i]))  
    #   i.ms <- unique(membership[i.nbr])  ## markets of firm i (the markets of the rivals of firm i)
    #   i.ms.row <- rep(0, length(markets)) ## dummy row for [Firm x Market] matrix
    #   i.ms.row[ i.ms ] <- 1 ## assign [Firm x Market] matrix row value to 1
    #   names(i.ms.row) <- sapply(1:length(markets),function(x)paste0('m',x))
    #   return(i.ms.row)
    # })
    # ## convert df to matrix
    # m.ms <- as.matrix(df.ms)
    ## bipartite firm-market
    # gb <- igraph::graph_from_incidence_matrix(as.matrix(df.ms), directed = F)
    # bl[[t]] <- 
  }
  names(nl) <- ys[(length(ys)-10+1):length(ys)]
  npds <- length(nl)
  
  firmnamesall <- unique(unlist(sapply(nl,function(x)V(x$gmmcsub)$vertex.names)))
  firmnames0 <- V(nl$`2007`$gmmcsub)$vertex.names
  adj0 <- as_adj(nl$`2007`$gmmcsub, sparse = F)
  n0 <- length(firmnames0)
  # firmnames <- unique(unlist(sapply(nl,function(x)V(x$gmmcsub)$vertex.names)))
  
  # ## which names in all years
  # fnidxl <- sapply(firmnamesall,function(x){
  #     all(sapply(nl, function(z)x %in% V(z$gmmcsub)$vertex.names))
  #   })
  # fnidx <- which(fnidxl)
  # firmnamesub <- firmnamesall[fnidx] ## firms in all yeras of gmmcsub MMC network
  # nsub <- length(firmnamesub)
  ## which firm made at least 1 acquisitions
  .filtername <- which(
    dfl$cent_deg_mmc>0 &
    dfl$ipo>0
  )
  firmnamesub <- unique(as.character(dfl$name[.filtername]))
  nsub <- length(firmnamesub)
  #*********************88888888888
  
  ##------------------------------------
  ## NETWORK DV ARRAY
  ##------------------------------------
  mmcarr <- array(0, dim=c(nsub, nsub, npds))
  for (t in 1:length(nl)) {
    x <- nl[[t]]
    # gmmc.e <- igraph::delete.edges(x$gmmc, which(E(x$gmmc)$weight < 2))
    # ## create subgraph removing isolates
    # gmmcsub <- igraph::induced.subgraph(gmmc.e, which(igraph::degree(gmmc.e) > 0))
    # vids <- which(
    #   V(x$gmmc)$vertex.names %in% firmnamesub & 
    #     V(x$gmmc)$vertex.names %in% V(gmmcsub)$vertex.names
    # )
    vids <- which(V(x$gmmc)$vertex.names %in% firmnamesub)
    mmcarr[,,t] <- as_adj(induced.subgraph(x$gmmc,vids), sparse = F)
  }
  # for (t in 1:length(nl)) {
  #   x <- nl[[t]]
  #   gmmcadj <- matrix(0, nrow=vcount(x$gmmc), ncol=vcount(x$gmmc))
  #   vids.gmmcsub <- which(V(x$gmmcsub)$vertex.names %in% firmnamesub)
  #   ## period t gmmcsub adjmat inserted into whole gmmc adjmat
  #   gmmcadj[vids.gmmcsub,vids.gmmcsub] <- as_adj(induced.subgraph(x$gmmcsub,vids.gmmcsub), sparse = F)
  #   ## subset all years firms from gmmc adjmat into allyears array
  #   vids.gmmc <- which(V(x$gmmc)$vertex.names %in% firmnamesub)
  #   mmcarr[,,t] <- gmmcadj[vids.gmmc,vids.gmmc]
  # }
  ## 2-year periods
  mmcarr2 <- array(0, dim=c(nsub, nsub, npds/2))
  for (t in 1:(npds/2)) {
    a <- (t-1)*2 + 1
    b <- (t-1)*2 + 2
    mab <- mmcarr[,,a] + mmcarr[,,b]
    mab[mab>0] <- 1
    mmcarr2[,,t] <- mab
  }
  
  ## COVARIATES (Acquisitions)
  arrAcq <- array(0,dim=c(nsub,npds))
  # arrW <- array(0,dim=c(nsub,npds))
  for (t in 1:length(nl)) {
    idx.acq <- which(dfl$year==ys[t] & dfl$name %in% firmnamesub)
    arrAcq[,t] <- round(log(1+dfl$rp_Acquisitions[idx.acq]))
  }
  
  ##-------------------
  ## SAOM
  ##------------------
  depMMC <- sienaDependent(mmcarr, type="oneMode", nodeSet=c('FIRMS'), 
                           sparse = FALSE, allowOnly = FALSE)
  # TREATMENT
  depAcq <- sienaDependent(arrAcq, type="behavior", nodeSet=c('FIRMS'), 
                             sparse = FALSE, allowOnly = FALSE)

  # Controls
  # moTrtLinearArray <- array(moTrtLinearWaves, dim=c(nrows, nyrs))
  # moCovTrtLinear <- varCovar(moTrtLinearArray, nodeSet="FIRMS")
  # #
  # activWavesArray <- array(activWaves, dim=c(nrows, nyrs))
  # 
  
  # NODES
  firms <- firmnamesub
  nrows <- length(firms)
  FIRMS <- sienaNodeSet(nrows, 'FIRMS', firms)
  
  ##
  sysDat <- sienaDataCreate(depMMC, depAcq,    ##dyCovTrt, dyCovTrtLinear, 
                               nodeSets=list(FIRMS) ) #,smoke1, alcohol
  #
  pharmaEffects <- getEffects(pharmaDat)
  #
  effectsDocumentation(pharmaEffects)
  # print01Report(pharmaDat, modelname="pharma-test-1")
  
  # pharmaEffects <- includeEffects(pharmaEffects, linear, name="depActiv")
  pharmaEffects <- includeEffects(pharmaEffects, outdeg, name="depActiv", interaction1 = 'depMMC')
  pharmaEffects <- includeEffects(pharmaEffects, outRate, name="depActiv", interaction1 = 'depMMC', type='rate')
  # pharmaEffects <- includeEffects(pharmaEffects, isolate, name="depActiv", interaction1 = 'depMMC')
  # pharmaEffects <- includeEffects(pharmaEffects, popAlt, name="depActiv", interaction1 = 'depMMC')
  pharmaEffects <- includeEffects(pharmaEffects, density, name="depMMC")
  # pharmaEffects <- includeEffects(pharmaEffects, inPop, name="depMMC")
  # pharmaEffects <- includeEffects(pharmaEffects, outAct, name="depMMC")
  pharmaEffects <- includeEffects(pharmaEffects, cycle4, name="depMMC")
  pharmaEffects <- includeEffects(pharmaEffects, sameXCycle4, name="depMMC", interaction1 = 'depActiv')
  #
  pharmaEffects <- includeTimeDummy(pharmaEffects, outdeg, name="depActiv", interaction1 = 'depMMC', timeDummy = '5,6,7,8')
  pharmaEffects <- includeTimeDummy(pharmaEffects, sameXCycle4, name="depMMC", interaction1 = 'depActiv', timeDummy = '5,6,7,8')
  
  #
  pharmaModel <- sienaAlgorithmCreate(projname='pharma-mmc-test-proj-1')
  
  pharmaResults <- siena07(pharmaModel, data=pharmaDat, effects=pharmaEffects, 
                           batch = T,   returnDeps = T, ## necessary for GOF
                           useCluster = T, nbrNodes = ncores, clusterType = 'PSOCK')
  
  screenreg(pharmaResults, digits = 4)
  
  
  
  # mmcl <- as.array(c(lapply(nl,function(x){
  #   vids <- which(V(x$gmmcsub)$vertex.names %in% firmnames0)
  #   as_adj(induced.subgraph(x$gmmcsub,vids), sparse = F)
  # })), 
  # dim=c(n0, n0, npds))
  
  
  
  
  
  
  
  
  
  
  ##=======================================
  ##  TEST PANEL DATA RERESSION
  ##---------------------------------------
  dfl$name <- as.character(dfl$name)
  yrmin <- min(unique(dfl$year))
  ipo.name <- dfl$name
  allyr.name <- dfl$name
  gte1acq.name <- dfl$name
  ## Do NOT filter Public
  #
  ##----------------
  ## filter public
  ##----------------
  ## get public firms to subset data
  ipo.name <- dfl$name[which(dfl$ipo==1)]
  # ipo.yrmin.name <- V(gmmc)$vertex.names[ which(V(gmmc)$ipo_status == 1) ]
  ## get firms by any(type==TRUE)
  names.type.true <- unique(dfl$name[dfl$type])
  ##----------------
  # ##----------------
  # ## filter Firms w/at least 1 acquisition
  # ##----------------
  # gte1acq.name <- unique(dfl$name[which(dfl$rp_Acquisitions > 0)])
  #
  # ##----------------
  # ## filter firms present in ALL Years
  # ##----------------
  # allyr.name <- unique(dfl$name[which(as.character(dfl$name) %in% firmnamesub)])
  #
  ## Regression data subset
  dfreg <- dfl[which(dfl$name %in% names.type.true & 
                       dfl$name %in% ipo.name & 
                       dfl$name %in% allyr.name &
                       dfl$name %in% gte1acq.name), ]
  # dfreg <- dfl[which(dfl$name %in% names.type.true), ]
  
  # ## compute feedback measure
  # dfreg$feedback1 <- NA
  # dfreg$feedback2 <- NA
  # dfreg$feedback3 <- NA
  # dfreg$feedback4 <- NA
  # dfreg$feedback5 <- NA
  # for (t in unique(dfreg$year)) {
  #   if (t == min(unique(dfreg$year))) {
  #     dfreg$feedback1[which(dfreg$year==t)] <- NA
  #     dfreg$feedback2[which(dfreg$year==t)] <- NA
  #     dfreg$feedback3[which(dfreg$year==t)] <- NA
  #     dfreg$feedback4[which(dfreg$year==t)] <- NA
  #   } else {
  #     dfreg$feedback1[which(dfreg$year==t)] <- dfreg$smmc1[which(dfreg$year==t)] - dfreg$smmc1[which(dfreg$year==(t-1))] 
  #     dfreg$feedback2[which(dfreg$year==t)] <- dfreg$smmc2[which(dfreg$year==t)] - dfreg$smmc2[which(dfreg$year==(t-1))] 
  #     dfreg$feedback3[which(dfreg$year==t)] <- dfreg$smmc3[which(dfreg$year==t)] - dfreg$smmc3[which(dfreg$year==(t-1))] 
  #     dfreg$feedback4[which(dfreg$year==t)] <- dfreg$smmc4[which(dfreg$year==t)] - dfreg$smmc4[which(dfreg$year==(t-1))] 
  #   }
  # }
  
  ##factors
  dfreg$i <- as.factor(dfreg$i)
  dfreg$year <- as.factor(dfreg$year)
  dfreg$y <- dfreg$y.cur + dfreg$y.new
  
  ## Competitive Pressure aggregations
  dfreg$wdeg_rp_all <- (dfreg$wdeg_rp_Acquisitions + dfreg$wdeg_rp_Capacity + 
                          dfreg$wdeg_rp_Legal + dfreg$wdeg_rp_Market_expansions + 
                          dfreg$wdeg_rp_Marketing + dfreg$wdeg_rp_New_product + 
                          dfreg$wdeg_rp_Pricing + dfreg$wdeg_rp_Strategic_alliances)
  ## acquisitions|cohort (1); non_acquisitions (7); relations (2); invariant (5)
  dfreg$wdeg_rp_NON_acquisitions <- dfreg$wdeg_rp_all - dfreg$wdeg_rp_Acquisitions
  dfreg$wdeg_rp_net_relations <- dfreg$wdeg_rp_Market_expansions + dfreg$wdeg_rp_New_product
  dfreg$wdeg_rp_net_restruct <- dfreg$wdeg_rp_Acquisitions + dfreg$wdeg_rp_net_relations
  dfreg$wdeg_rp_net_invariant <- dfreg$wdeg_rp_NON_acquisitions - dfreg$wdeg_rp_net_relations
  
  ## complexity
  xcols <- c( "wdeg_rp_Acquisitions", "wdeg_rp_Capacity",           
    "wdeg_rp_Legal", "wdeg_rp_Market_expansions",  "wdeg_rp_Marketing",          
    "wdeg_rp_New_product", "wdeg_rp_Pricing", "wdeg_rp_Strategic_alliances")
  dfreg$wdeg_rp_complexity <- apply(dfreg[,xcols],1,complexity)
  
  ## MUTATE CONTROLS
  dfreg$dum.crisis <- 0
  dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
  dfreg$acq_sum_1_sc <- scale(dfreg$acq_sum_1, center = T, scale = T)
  dfreg$employee_na_age_sc <- scale(dfreg$employee_na_age, center = T, scale = T)
  dfreg$sales_na_0_mn_sc <- scale(dfreg$sales_na_0_mn, center = T, scale = T)
  dfreg$cent_deg_all_sc <- scale(dfreg$cent_deg_all, center = T, scale = T)
  
  # ## remove NAs in first year due to feedback lag
  # dfreg <- dfreg[!is.na(dfreg$feedback1),]
  
  # ## check 2-year periods
  # dfreg$yeven <- sapply(as.character(dfreg$year),function(x)as.numeric(x)%%2==0)
  # dfreg$year2 <- as.factor(sapply(as.character(dfreg$year),function(y){
  #     y <- as.numeric(y)
  #     if (y <= 2008) return(1)
  #     if (y <= 2010) return(2)
  #     if (y <= 2012) return(3)
  #     if (y <= 2014) return(4)
  #     if (y <= 2016) return(5)
  #   }))
  # dfreg0 <- dfreg[dfreg$yeven,]
  # dfreg1 <- dfreg[!dfreg$yeven,]
  # #
  # dfreg1$rp_Acquisitions <- dfreg1$rp_Acquisitions + dfreg0$rp_Acquisitions
  # dfreg1$rp_all <- dfreg1$rp_all + dfreg0$rp_all
  
  dfreg$rp_Acquisitions_bin <- as.integer(dfreg$rp_Acquisitions > 0)
  
  
  # ## filter observations with enough MMC rivals to be relevant ???
  dfreg <- dfreg[which(dfreg$cent_deg_mmc > 0),]
  # 
  
  ##===================================
  ## NEW RAVENPACK Acquisition DV - ACTIONS MODERATOR
  ##-----------------------------------
  pow <- 'smmc4n'
  ##-----------------------------------
  
  # dfreg$smmc4nz <- as.numeric(scale(dfreg$smmc4n))
  mract0 <- pglm(rp_Acquisitions ~ dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + 
                   # log(1 + cent_deg_all) +
                   log(1+cent_deg_mmc) +
                   I(log(1+cent_deg_mmc)^2),
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  # summary(mract0)
  mract1 <- pglm(rp_Acquisitions ~  dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + 
                   # log(1 + cent_deg_all) +
                   # log(1+wdeg_rp_net_invariant) + 
                   # log(1+wdeg_rp_net_relations) + 
                   log(1+wdeg_rp_NON_acquisitions) +
                   log(1+wdeg_rp_Acquisitions) + 
                   log(1+cent_deg_mmc) +
                   log(1+cent_deg_mmc):log(1+wdeg_rp_Acquisitions) + 
                   I(log(1+cent_deg_mmc)^2) +
                   I(log(1+cent_deg_mmc)^2):log(1+wdeg_rp_Acquisitions),
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  # summary(mract1)
  screenreg(list(mract0,mract1),single.row = T)
  
  
  
  
  mract0 <- pglm(rp_Acquisitions ~  dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + 
                   log(1 + cent_deg_all) +
                   # log(1+wdeg_rp_NON_acquisitions) +
                   # log(1+wdeg_rp_net_invariant) + 
                   # log(1+wdeg_rp_net_relations) + 
                   # log(1+wdeg_rp_Acquisitions) + 
                   smmc35n +
                   I(smmc35n^2) ,
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  summary(mract0)
  
  mract1 <- pglm(rp_Acquisitions ~  dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + 
                   log(1 + cent_deg_all) +
                   log(1+wdeg_rp_NON_acquisitions) +
                   # log(1+wdeg_rp_net_invariant) + 
                   # log(1+wdeg_rp_net_relations) + 
                   log(1+wdeg_rp_Acquisitions) + 
                   smmc35n +
                   smmc35n:log(1+wdeg_rp_Acquisitions) + 
                   smmc35n:log(1+wdeg_rp_net_relations) + 
                   smmc35n:log(1+wdeg_rp_net_invariant) + 
                   I(smmc35n^2) +
                   I(smmc35n^2):log(1+wdeg_rp_Acquisitions) + 
                   I(smmc35n^2):log(1+wdeg_rp_net_relations) + 
                   I(smmc35n^2):log(1+wdeg_rp_net_invariant),
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  summary(mract1)
  screenreg(list(mract0,mract1),single.row = T)
  
  
  mract1 <- pglm(rp_Acquisitions ~  dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + 
                   # log(1 + cent_deg_all) +
                   log(1+wdeg_rp_net_relations) +
                   log(1+wdeg_rp_net_invariant) +
                   log(1+wdeg_rp_Acquisitions) +
                   log(1+cent_deg_mmc) +
                   log(1+cent_deg_mmc):log(1+wdeg_rp_net_relations) +
                   log(1+cent_deg_mmc):log(1+wdeg_rp_net_invariant) +
                   log(1+cent_deg_mmc):log(1+wdeg_rp_Acquisitions) +
                   I(log(1+cent_deg_mmc)^2) +
                   I(log(1+cent_deg_mmc)^2):log(1+wdeg_rp_net_relations) +
                   I(log(1+cent_deg_mmc)^2):log(1+wdeg_rp_net_invariant) +
                   I(log(1+cent_deg_mmc)^2):log(1+wdeg_rp_Acquisitions),
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  
  
  
  
  
  mract2 <- pglm(rp_Acquisitions ~  dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                   wdeg_rp_Acquisitions +
                   wdeg_rp_NON_acquisitions +
                   smmc4n +
                   smmc4n:wdeg_rp_Acquisitions + 
                   smmc4n:wdeg_rp_NON_acquisitions + 
                   I(smmc4n^2) +
                   I(smmc4n^2):wdeg_rp_Acquisitions + 
                   I(smmc4n^2):wdeg_rp_NON_acquisitions,
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  summary(mract2)
  
  
  mract2b <- pglm(rp_Acquisitions ~  dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                    I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                    I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                    wdeg_rp_net_restruct +
                    wdeg_rp_net_invariant +
                    smmc4n +
                    smmc4n:wdeg_rp_net_restruct + 
                    smmc4n:wdeg_rp_net_invariant + 
                    I(smmc4n^2) +
                    I(smmc4n^2):wdeg_rp_net_restruct + 
                    I(smmc4n^2):wdeg_rp_net_invariant,
                  data=dfreg, family = poisson, 
                  model = 'within', effect = 'twoways',
                  R = 100, method='nr',
                  index=c('i','year'))
  summary(mract2b)
  
  mract3 <- pglm(rp_Acquisitions ~ dum.crisis + acq_cnt_5 + rp_NON_acquisitions +
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                   wdeg_rp_Acquisitions + 
                   wdeg_rp_net_relations +
                   wdeg_rp_net_invariant +
                   smmc4n +
                   smmc4n:wdeg_rp_Acquisitions + 
                   smmc4n:wdeg_rp_net_relations +
                   smmc4n:wdeg_rp_net_invariant +
                   I(smmc4n^2) +
                   I(smmc4n^2):wdeg_rp_Acquisitions +
                   I(smmc4n^2):wdeg_rp_net_relations +
                   I(smmc4n^2):wdeg_rp_net_invariant,
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
  summary(mract3)  
  
  
  
  #
  screenreg(list(mract0, mract1, mract2, mract2b, mract3), digits = 3, single.row = T)
  
  #
  htmlfile <- sprintf('acq_sys_mmc_panel_RPactions_regression_RESULTS_table_ego-%s_d%s_nrow%s_%s.html', 
                      firm_i_ego, d, nrow(dfreg), pow)
  htmlreg(list(mract0, mract1, mract2, mract2b, mract3), 
          file = file.path(result_dir,htmlfile), digits = 3)
  
    
# }


















