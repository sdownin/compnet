##--------------------------------------------------------------
##
##  MMC & ACQUISITIONS 
##
##--------------------------------------------------------------
# .libPaths('C:/Users/steph/Documents/R/win-library/3.2')
library(igraph)
library(network)
library(intergraph)
library(pglm)
library(lme4)
library(texreg)
library(stringdist)
library(plyr)
library(stringr)
library(lubridate)

## DIRECTORIES
work_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2"
nets_dir <- file.path(work_dir,"firm_nets_rnr2","firmctrl")
nets_repl_dir <- file.path(work_dir, 'firm_nets_rnr2_repl')
version_dir <- file.path(work_dir,'R','acqmmc_os_v1')
data_dir <- file.path(work_dir,'acqmmc_os_v1','data')
result_dir <- file.path(work_dir,'acqmmc_os_v1','data')
rvn_dir <- 'C:/Data/RavenPack'

setwd(work_dir)


## LOAD DATA AND DEPENDENCIES
acf    <- source(file.path(version_dir,'acq_compnet_functions.R'))$value    ## acf: awareness functions
cb     <- source(file.path(version_dir,'acq_cb_data_prep.R'))$value           ## cb : CrunchBase
# sdc    <- source(file.path(version_dir,'acq_sdc_coop.R'))$value               ## sdc: Thompson SDC
# si     <- source(file.path(version_dir,'acq_firm_size_controls.R'))$value     ## si : size controls from mergent intellect
# ih     <- source(file.path(version_dir,'acq_institutional_holdings.R'))$value ## ih : institutional holdings
# g.full <- source(file.path(version_dir,'acq_make_full_graph.R'))$value        ## g.full : full competition graph

# ## set firms to create networks (focal firm or replication study focal firms)



##============================
##   ACTION CATEGORIES
##----------------------------
# dtafile <- 'rp40_action_categories_2000_2018.dta'
# csvfile <- 'rp40_action_categories_2000_2018.csv'
csvfile <- 'ravenpack_resolved_actions_2008-2016.csv'     ## %%CHANGED%%
##convert
rvn <- read.csv(file.path(rvn_dir, csvfile), stringsAsFactors = F)

# ## all data 7GB
# dtafile <- 'rp40_2000_2016.dta'
# rvn.all <- read.dta13(file.path(rvn_dir, dtafile), select.rows=10)



##===============================
## set analysis params
##-------------------------------
firms.todo <- c('cisco')
firm_i <- firms.todo[1]
firm_i_ego <- firm_i
d <- 3
ncpus <- 4
parallel <- "multicore"
nPeriods <- 11  ## 5
startYr <- 2006
endYr <- 2017            ## dropping first for memory term; actual dates 2007-2016
##-------------------------------



##----------------------------
## LOAD   
##   from list file or collection of period-files to be combined into a list
##----------------------------
data_file <- file.path(nets_dir,sprintf('%s_d%s.rds',firm_i,d))
if (file.exists(data_file)) {
  ## 1. load networks list file
  nets <- readRDS(data_file)
} else {
  ## 2. Combine nets files for larger networks
  nets <- list()
  for (t in startYr:endYr) {
    data_file <- file.path(nets_repl_dir,sprintf('%s_d%s_y%s.rds',firm_i,d,t-1))
    nets[[length(nets)+1]] <- readRDS(data_file)
  }
  names(nets) <- as.integer(startYr:endYr) 
}

#filter list periods
if (nPeriods < length(nets))   
  nets <- nets[(length(nets)-nPeriods+1):length(nets)] 

#periods
periods <- as.integer(names(nets)) -1




##================================
##
##    FUNCTIONS
##
##-------------------------------

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

##
# Gets the count of actions of given type for given firms in RavenPack data
# @param [dataframe] rvnsub   RavenPack subset by period (year)
# @param [string[n]] names    CrunchBase company_name_unique vector for network vertices
# @param [dataframe] cbrpmap  CB-to-RP name (entity ID) mapping
# @param [string] var_rp      Name of action type in RavenPack dataframe 
# @return [integer[n]]        The action counts vector
##
getRpActionCount <- function(rvnsub, names, cbrpmap, var_rp)
{
  n <- length(names)
  x <- rep(0,n)
  
  if (length(rvnsub)==0 | nrow(rvnsub)==0)
    return(x)
  
  ## network vertex indices of firms in RavenPack data
  vids <- which(names %in% cbrpmap$company_name_unique)
  for (i in 1:length(vids)) 
  {
    actids <- which(
      rvnsub$rp_entity_id==cbrpmap$rp_entity_id[i] & 
      rvnsub$action_category==var_rp
    )
    v <- vids[i]  ## the v'th firm in network is the i'th RP entity
    x[v] <- length(actids) ## count of actions is the length of the action ids vector
  }
  
  return(x)
}

##
# Gets each vertex's sum of neighbors' values of given attribute (attr)
# @param [igraph] gt    Igraph object (for period t)
# @param [integer] vid  Vertex index of focal firm to fetch neighbors' attributes
# @param [string] attr  Name of vertex attribute
# @return [integer[m]]  Vector of neighbors attribute values (for m neighbors)
##
getNeighborsAttrSum <- function(gt, attr) {
  adjmat <- igraph::get.adjacency(gt)
  attr <- igraph::get.vertex.attribute(gt, attr)
  return(as.numeric(adjmat %*% attr))
}
# getNeighborsAttr <- function(gt, vid, attr) {
#   nvids <- as.integer(igraph::neighbors(gt, vid))
#   return(igraph::get.vertex.attribute(gt, attr, nvids))
# }


##==========================================
##
## Ravepack 
##
##------------------------------------------
## Load RavenPack
g0 <- asIgraph(nets[[1]])
firms.all.df <- data.frame(
  company_name_unique=as.character(V(g0)$vertex.names),
  company_name=as.character(V(g0)$company_name),
  stringsAsFactors = F
)
## assign unique name to company_name if missing
idxm <- which(is.na(firms.all.df$company_name) | str_to_lower(firms.all.df$company_name)=='na')
firms.all.df$company_name[idxm] <- firms.all.df$company_name_unique[idxm]

fdf <- igraph::as_data_frame(g0, what='vertices')
# View(head(fdf, 30))

## firms that have sales > 0 in ANY period
sidx <- sapply(nets, function(net) which(net %v% 'sales_na_0' > 0))
sidx <- sort(unique(unlist(sapply(nets, function(net) which(net %v% 'sales_na_0' > 0)))))
firms <- firms.all.df[sidx, ]
# View(firms)



##------------------------------------------
## RavenPack to CrunchBase Name Map
##------------------------------------------
# RAVENPACK NAMES
rpn <- unique(rvn[,c('rp_entity_id','entity_name')])
rpn$name.l <- str_to_lower(rpn$entity_name)

# COHORT OF FIRMS
cbrpmap <- data.frame(
  company_name_unique = firms$company_name_unique,
  company_name = firms$company_name,
  idx1 = NA, ## mode
   idx1pct = NA, ## confidence (%)
   rp_entity_id_1 = NA, ## top match ID
   entity_name_1 = NA, ##  top match name
  idx2 = NA, ## 
   idx2pct = NA, ## 
   rp_entity_id_2 = NA, ## 
   entity_name_2 = NA, ## 
  idx3 = NA, ## 
   idx3pct = NA, ## 
   rp_entity_id_3 = NA, ## 
   entity_name_3 = NA, ## 
  stringsAsFactors = F
)

## checking for string in ravenpack firm names
# rvn[grep('duke',rvn$entity_name,ignore.case = T,perl = T), c('entity_name','rp_entity_id') ]
rpn$entity_name[grep('mapsense',rpn$entity_name, ignore.case = T, perl=T)]

## string distance matrix for CrunchBase-to-RavenPack Name Mapping
for (i in 1:nrow(firms)) {
  firm_i <- firms$company_name[i]
  
  if (grepl('\\(now',firm_i,ignore.case=T,perl=T)) {
    ## use new name if "old_name (now new_name)"
    firm_i <- gsub('(.+\\(now\\s+|\\))', '', firm_i)
  } else if (grepl('[|]',firm_i,ignore.case=T,perl=T)) {
    ## use full name if "abbrev | full name"
    .parts <- str_trim(unlist(strsplit(firm_i,'[|]')), side='both')
    firm_i <- .parts[length(.parts)]
  }
  
  if (i %% 10 == 0) cat(sprintf('i=%03s (%03.2f%s) %s \n',i,100*i/nrow(firms),'%',firms$company_name_unique[i]))
  ##
  ## FIRST DO grep() string search, see if excatly match 1
  ##  - else if multiple matches, test distance on matches
  ##  - else test distance on all unmateched RavenPack firms
  ##
  
  grepidx <- grep(sprintf('^%s',firm_i),rpn$entity_name,ignore.case=T,perl=T)
  if (length(grepidx)==1) {
    
    cbrpmap$idx1[i] <- grepidx
    cbrpmap$idx1pct[i] <- 1
    cbrpmap$rp_entity_id_1[i] <- rpn$rp_entity_id[grepidx]
    cbrpmap$entity_name_1[i] <- rpn$entity_name[grepidx]
    
  } else {
    
    rpn_i <- if (length(grepidx) > 1) { rpn[grepidx, ] } else { rpn }
    
    idf <- data.frame(
      name=rpn_i$entity_name,
      name.l=rpn_i$name.l,
      jw=stringdist(rpn_i$name.l, firm_i, method = 'jw'),
      lcs=stringdist(rpn_i$name.l, firm_i, method = 'lcs'),
      cosine=stringdist(rpn_i$name.l, firm_i, method = 'cosine'),
      jaccard=stringdist(rpn_i$name.l, firm_i, method = 'jaccard'),
      qgram=stringdist(rpn_i$name.l, firm_i, method = 'qgram'),
      osa=stringdist(rpn_i$name.l, firm_i, method = 'osa')#,
      # lv=stringdist(rpn_i$name.l, firm_i, method = 'lv'),
      # dl=stringdist(rpn_i$name.l, firm_i, method = 'dl'),
    )
    idx <- apply(idf[,3:ncol(idf)], 2, which.min)
    idxcnt <- plyr::count(idx)
    idxcnt <- merge(idxcnt, data.frame(metric=names(idx),x=idx), by='x')
    
    if (all(idxcnt$freq == 1)) {
      ## keep original order from above, if all metrics only = 1
      idxcnt <- idxcnt[sapply(idx,function(x)which(idxcnt$x==x)), ]
    }  else if (max(idxcnt$freq) <= 2) {
      ## keep 'jw' first then descending order by metric freq
      idxcnt <- idxcnt[order(idxcnt$freq, decreasing = T),]
      jwrow <- which(idxcnt$metric=='jw')
      idxcnt <- rbind(idxcnt[jwrow,],idxcnt[-jwrow,])
    }  else {
      ## else order by descending metric freq
      idxcnt <- idxcnt[order(idxcnt$freq, decreasing = T),]
    }
    
    if (nrow(idxcnt)==0 | length(idxcnt)==0) {
      cat(sprintf('missing strdist for i=%03s, `%s\n`',i,firms$company_name_unique[i]))
      next
    }
    
    ## top match
    idx1 <- idxcnt$x[1]
    cbrpmap$idx1[i] <- idx1
    cbrpmap$idx1pct[i] <- idxcnt$freq[1] / length(idxcnt$freq)
    cbrpmap$rp_entity_id_1[i] <- rpn_i$rp_entity_id[idx1]
    cbrpmap$entity_name_1[i] <- rpn_i$entity_name[idx1]
    ## secondary / tertiary matches
    if (nrow(idxcnt) >= 2) {
      idx2 <- idxcnt$x[2]
      cbrpmap$idx2[i] <- idx2
      cbrpmap$idx2pct[i] <- idxcnt$freq[2] / length(idxcnt$freq)
      cbrpmap$rp_entity_id_2[i] <- rpn_i$rp_entity_id[idx2]
      cbrpmap$entity_name_2[i] <- rpn_i$entity_name[idx2]
    }
    if (nrow(idxcnt) >= 3) {
      idx3 <- idxcnt$x[3]
      cbrpmap$idx3[i] <- idx3
      cbrpmap$idx3pct[i] <- idxcnt$freq[3] / length(idxcnt$freq)
      cbrpmap$rp_entity_id_3[i] <- rpn_i$rp_entity_id[idx3]
      cbrpmap$entity_name_3[i] <- rpn_i$entity_name[idx3]
    }
    
  }
  
}

## save mapping dataframe
cbrpmapfile <- sprintf('cb_ravenpack_name_map_ego-%s_d%s.csv',firm_i_ego,d)
write.csv(cbrpmap, file=file.path(result_dir,cbrpmapfile))

# ## Create dfreg dataframe
# ## Filter out firms that don't have RavenPack ID
# dfreg <- cbrpmap[!is.na(cbrpmap$rp_entity_id),]

## Fill in values to new firm other manually checked firm
ckfirm <- 'microsoft'
ckmatch <- read.csv(file.path(result_dir,sprintf('cb_ravenpack_name_map_ego-%s_d%s.csv',ckfirm,d)), stringsAsFactors = F)
cbrpmap <- merge(cbrpmap, ckmatch[,c('OK','company_name_unique')], by ='company_name_unique', all.x=T, all.y=F)
colnms <- names(cbrpmap)
cbrpmap <- cbrpmap[,c(which(colnms=='OK'),which(colnms!='OK'))]
write.csv(cbrpmap, file=file.path(result_dir,cbrpmapfile))
#










# firms.todo <- c('ibm','google')
firms.todo <- c('cisco','google') ## c('apple','amazon','facebook','google')
for (firm_i_ego in firms.todo)
{
  
  
  
  ##--------------------------------------------------------
  ##
  ## AFTER Creating CB-RP Matches List
  ##
  ##--------------------------------------------------------
  firm_i <- firm_i_ego
  name_i <- firm_i
  #
  cbrpmapfile <- sprintf('cb_ravenpack_name_map_ego-%s_d%s.csv',firm_i_ego,d)
  ## LOAD CB-Ravenpack approximate matches
  cbrpmatch <- read.csv(file=file.path(data_dir,cbrpmapfile), 
                        stringsAsFactors = F, na.strings = c("'","-","'-"))
  cbrpmatch$name <- cbrpmatch$company_name_unique
  ## Filter subset of firms with match
  cbrpmapsub <- cbrpmatch[!is.na(cbrpmatch$OK),]
  ## Create CB-RP Map
  cbrpmap <- data.frame(stringsAsFactors = F)
  for (i in 1:nrow(cbrpmapsub)) 
  {
    if(cbrpmapsub$OK[i] == 1) {
      rp_entity_id <- cbrpmapsub$rp_entity_id_1[i]
      rp_entity_name <- cbrpmapsub$entity_name_1[i]
    } else if (cbrpmapsub$OK[i] == 2) {
      rp_entity_id <- cbrpmapsub$rp_entity_id_2[i]
      rp_entity_name <- cbrpmapsub$entity_name_2[i]
    } else if (cbrpmapsub$OK[i] == 3) {
      rp_entity_id <- cbrpmapsub$rp_entity_id_3[i]
      rp_entity_name <- cbrpmapsub$entity_name_3[i]
    } else {
      next
    }
    cat(sprintf('i=%s, OK=%s, rp_ent_id=%s\n',
                i, cbrpmapsub$OK[i], rp_entity_id))
    cbrpmap <- rbind(cbrpmap, data.frame(
      company_name_unique=cbrpmapsub$name[i],
      rp_entity_id=rp_entity_id, 
      rp_entity_name=rp_entity_name,
      stringsAsFactors = F
    )
    )
  }
  
  head(cbrpmap)
  
  ##----------end Ravenpack name map----------
  
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
  
  
  
  
  ##----------------------------------------
  ##  Compute Ravenpack ACTIONS per Firm-Year
  ##----------------------------------------
  # ## init DVs
  # dfreg$rp_Acquisitions <- 0
  # dfreg$rp_Capacity <- 0
  # dfreg$rp_Legal <- 0
  # dfreg$rp_Market_expansions <- 0
  # dfreg$rp_Marketing <- 0
  # dfreg$rp_New_product <- 0
  # dfreg$rp_Pricing <- 0
  # dfreg$rp_Strategic_alliances <- 0
  # dfreg$rp_NON_acquisitions <- 0
  # dfreg$rp_all <- 0
  # dfreg$rp_complexity <- 0
  ##  tmp mapping df
  .rptmp <- data.frame(
    rp = c('Acquisitions', 'Capacity_related_actions',
           'Legal_actions', 'Market_expansion', 'Marketing_actions',
           'New_Product_Action', 'pricing_actions',
           'strategic_alliance'),
    dfreg= c('rp_Acquisitions','rp_Capacity','rp_Legal',
             'rp_Market_expansions','rp_Marketing',
             'rp_New_product','rp_Pricing',
             'rp_Strategic_alliances'),
    val = rep(0,8),
    stringsAsFactors = F
  )
  
  
  
  
  ##--------------------------------------------------------
  ##--------------------------------------------------------
  ##------------Network Time Period List--------------------
  ##--------------------------------------------------------
  ##--------------------------------------------------------
  # load('acq_sys_mmc_analysis_test_1.Rdata')
  nl <- list(gt=list(),gmmc=list(),gmmcsub=list())
  dfl <- data.frame()
  
  ## look through period 1:(T-1)
  ##   use node-collapsed network structure from
  ##   year 1 (lagged year) to predict year 2 acquisitions (1st non-lag year)
  ##   year T-1 to predict year T aquisitions (last year in sample)
  for (t in 1:(length(nets)-1)) 
  {
    cat(sprintf('t=%s\n',t))
    ## year 1 end of year collapsed network to predict year 2 aquisitions
    ## drop year 1 lag (so year 2 is year 1)
    net <- nets[[t]]
    
    ##=================================================================
    ## MMC & Acquisitions
    ## -- CRUNCHBASE ACQUISITIONS DATA
    ##-----------------------------------------------------------------
    gt <- asIgraph(net)
    V(gt)$com_multilevel <- igraph::multilevel.community(gt)$membership
    
    ## Create [Firm x Market] incidence matrix
    ## markets of each firm (which markets they are in based on NC membership of their rivals)
    cat('  computing firm-market matrix\n')
    membership <- V(gt)$com_multilevel
    markets <- unique(membership)
    markets <- markets[order(markets)]
    adj <- igraph::as_adjacency_matrix(gt, sparse = F)
    df.ms <- ldply(1:nrow(adj), function(i){
      i.nbr <- unname(which( adj[i, ] == 1 ))  ## rivals (neighbors in the network) ## i.nbr <- as.integer(igraph::neighbors(g, V(g)[i]))  
      i.ms <- unique(membership[i.nbr])  ## markets of firm i (the markets of the rivals of firm i)
      i.ms.row <- rep(0, length(markets)) ## dummy row for [Firm x Market] matrix
      i.ms.row[ i.ms ] <- 1 ## assign [Firm x Market] matrix row value to 1
      names(i.ms.row) <- sapply(1:length(markets),function(x)paste0('m',x))
      return(i.ms.row)
    })
    ## convert df to matrix
    m.ms <- as.matrix(df.ms)
    ## bipartite firm-market
    gb <- igraph::graph_from_incidence_matrix(m.ms, directed = F)
    ## MMC weighted unipartite projections
    gmmc <- igraph::bipartite.projection(gb, multiplicity = T)$proj1
    E(gmmc)$invweight <- 1/E(gmmc)$weight
    ##--------------------
    ##  SET COVARIATES
    ##--------------------
    cat(sprintf('  computing controls\n'))
    V(gmmc)$vertex.names <- V(gt)$vertex.names
    V(gmmc)$ipo_status <- V(gt)$ipo_status
    V(gmmc)$employee_na_age <- V(gt)$employee_na_age
    V(gmmc)$sales_na_0_mn <- V(gt)$sales_na_0_mn
    V(gmmc)$cent_deg_all <- igraph::degree(gt)
    ## CRUNCHBASE ACQUISTIONS
    a10<- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t]-10) & cb$co_acq$acquired_year < periods[t]), ]
    a5 <- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t]-5) & cb$co_acq$acquired_year < periods[t]), ]
    a2 <- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t]-2) & cb$co_acq$acquired_year < periods[t]), ]
    a1 <- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t]-1) & cb$co_acq$acquired_year < periods[t]), ]
    #
    V(gmmc)$acq_cnt_10 <- 0
    V(gmmc)$acq_cnt_5 <- 0
    V(gmmc)$acq_cnt_1 <- 0
    V(gmmc)$acq_sum_1 <- 0
    V(gmmc)$acq_mean_1 <- 0
    V(gmmc)$acq_median_1 <- 0
    V(gmmc)$acq_sum_2 <- 0
    V(gmmc)$acq_mean_2 <- 0
    V(gmmc)$acq_median_2 <- 0
    
    for (v in 1:vcount(gmmc)) {
      firm_v <- V(gmmc)$vertex.names[v]
      ## ACQUISITION EXPERIENCE (count of acquisitions in last 5 years)
      V(gmmc)$acq_cnt_10[v] <- length(which(a10$acquirer_name_unique==firm_v))
      V(gmmc)$acq_cnt_5[v] <- length(which(a5$acquirer_name_unique==firm_v))
      V(gmmc)$acq_cnt_1[v] <- length(which(a1$acquirer_name_unique==firm_v))
      ## LACK OF SLACK RESOURCES (Value of Acquisitions in previous year)
      vals1 <- a1$price_usd[which(a1$acquirer_name_unique==firm_v)]
      vals2 <- a2$price_usd[which(a2$acquirer_name_unique==firm_v)]
      V(gmmc)$acq_sum_1[v] <- sum(vals1, na.rm = T)
      V(gmmc)$acq_mean_1[v] <-  mean(vals1, na.rm = T)
      V(gmmc)$acq_median_1[v] <-  median(vals1, na.rm = T)
      V(gmmc)$acq_sum_2[v] <-  sum(vals2, na.rm = T)
      V(gmmc)$acq_mean_2[v] <-  mean(vals2, na.rm = T)
      V(gmmc)$acq_median_2[v] <-  median(vals2, na.rm = T)
    }
    # ##-----------------
    # ## COMPUTE Competitive Pressure
    # ##-----------------
    # pres1 <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmc, exp= 0.1), error = function(e)e)
    # pres2 <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmc, exp= 0.2), error = function(e)e)
    # pres3 <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmc, exp= 0.3), error = function(e)e)
    # pres4 <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmc, exp= 0.4), error = function(e)e)
    # pres1n <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmc, exp=-0.1), error = function(e)e)
    # pres2n <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmc, exp=-0.2), error = function(e)e)
    # pres3n <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmc, exp=-0.3), error = function(e)e)
    # pres4n <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmc, exp=-0.4), error = function(e)e)
    # preseig <- igraph::eigen_centrality(gmmc)$vector
    # presconstr <- igraph::constraint(gmmc)
    # presclose <- igraph::closeness(gmmc, normalized = T)
    # presbetw  <- igraph::betweenness(gmmc, normalized = T)
    # # presalphaw  <- igraph::alpha.centrality(gmmc, weights = 'weight', sparse = T)
    # # presalphainvw  <- igraph::alpha.centrality(gmmc, weights = 'invweight', sparse = T)
    # presauthw <- igraph::authority.score(gmmc)$vector
    # presdivers <- igraph::diversity(gmmc)
    ##-----------------
    ## COMPUTE SYSTEM MMC VECTORS
    ##-----------------
    cat('  computing system MMC measure\n')
    ## remove edges weight < 2
    gmmc.e <- igraph::delete.edges(gmmc, which(E(gmmc)$weight < 2))
    ## create subgraph removing isolates
    gmmcsub <- igraph::induced.subgraph(gmmc.e, which(igraph::degree(gmmc.e) > 0))
    E(gmmcsub)$invweight <- 1/E(gmmcsub)$weight
    mmcnames <- V(gmmcsub)$vertex.names
    ## Centrality -- SYSTEMS MMC for only firms with MMC eges
    pc0.1 <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmcsub, exp= 0.1), error = function(e)e)
    # pc0.2 <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmcsub, exp= 0.2), error = function(e)e)
    # pc0.3 <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmcsub, exp= 0.3), error = function(e)e)
    pc0.4 <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmcsub, exp= 0.4), error = function(e)e)
    #
    pcn0.1  <- tryCatch(tmpn0.1 <- igraph::power_centrality(gmmcsub, exp=-0.1), error = function(e)e)
    pcn0.2  <- tryCatch(tmpn0.2 <- igraph::power_centrality(gmmcsub, exp=-0.2), error = function(e)e)
    pcn0.3  <- tryCatch(tmpn0.3 <- igraph::power_centrality(gmmcsub, exp=-0.3), error = function(e)e)
    pcn0.35 <- tryCatch(tmpn0.35<- igraph::power_centrality(gmmcsub, exp=-0.35), error = function(e)e)
    pcn0.4  <- tryCatch(tmpn0.4 <- igraph::power_centrality(gmmcsub, exp=-0.4), error = function(e)e)
    pcn0.45 <- tryCatch(tmpn0.45<- igraph::power_centrality(gmmcsub, exp=-0.45), error = function(e)e)
    pcn0.5  <- tryCatch(tmpn0.5 <- igraph::power_centrality(gmmcsub, exp=-0.5), error = function(e)e)
    ## init vector of system MMC for ALL singlemarket and multimarket firms
    smmc1n <- rep(0, vcount(gmmc))
    smmc2n <- rep(0, vcount(gmmc))
    smmc3n <- rep(0, vcount(gmmc))
    smmc35n <- rep(0, vcount(gmmc))
    smmc4n <- rep(0, vcount(gmmc))
    smmc45n <- rep(0, vcount(gmmc))
    smmc5n <- rep(0, vcount(gmmc))
    #
    smmc1 <- rep(0, vcount(gmmc))
    # smmc2 <- rep(0, vcount(gmmc))
    # smmc3 <- rep(0, vcount(gmmc))
    smmc4 <- rep(0, vcount(gmmc))
    smmceig <- rep(0, vcount(gmmc))
    smmcconstr <- rep(0, vcount(gmmc))
    smmcclose <- rep(0, vcount(gmmc))
    smmcbetw <- rep(0, vcount(gmmc))
    # smmcalphaw<- rep(0, vcount(gmmc))
    # smmcalphainvw<- rep(0, vcount(gmmc))
    smmcauthw <- rep(0, vcount(gmmc))
    smmcdivers <- rep(0, vcount(gmmc))
    ## indices of mmc firm subset within full graph
    idx.sub.in.mmc <- unname(sapply(mmcnames, function(x) which(V(gmmc)$vertex.names == x)))
    ## assign subset MMC firms' system MMC into vector for all firms
    smmc1n[idx.sub.in.mmc]  <- pcn0.1
    smmc2n[idx.sub.in.mmc]  <- pcn0.2
    smmc3n[idx.sub.in.mmc]  <- pcn0.3
    smmc35n[idx.sub.in.mmc] <- pcn0.35
    smmc4n[idx.sub.in.mmc]  <- pcn0.4
    smmc45n[idx.sub.in.mmc] <- pcn0.45
    smmc5n[idx.sub.in.mmc]  <- pcn0.5
    #
    smmc1[idx.sub.in.mmc] <- pc0.1
    # smmc2[idx.sub.in.mmc] <- pc0.2
    # smmc3[idx.sub.in.mmc] <- pc0.3
    smmc4[idx.sub.in.mmc] <- pc0.4
    #
    smmceig[idx.sub.in.mmc] <- igraph::eigen_centrality(gmmcsub)$vector
    smmcconstr[idx.sub.in.mmc] <- igraph::constraint(gmmcsub)
    #
    smmcclose[idx.sub.in.mmc] <- igraph::closeness(gmmcsub, normalized = T)
    smmcbetw[idx.sub.in.mmc]  <- igraph::betweenness(gmmcsub, normalized = T)
    # smmcalphaw[idx.sub.in.mmc]  <- igraph::alpha.centrality(gmmcsub, weights = 'weight', sparse = T)
    # smmcalphainvw[idx.sub.in.mmc]  <- igraph::alpha.centrality(gmmcsub, weights = 'invweight', sparse = T)
    #
    smmcauthw[idx.sub.in.mmc] <- igraph::authority.score(gmmcsub)$vector
    smmcdivers[idx.sub.in.mmc] <- igraph::diversity(gmmcsub)
    ##-----------------
    ## Number of competitors
    ##-----------------
    cent_deg_mmc        <- rep(0, vcount(gmmc))
    cent_deg_mmc_weight <- rep(0, vcount(gmmc))
    # 
    cent_deg_mmc[idx.sub.in.mmc]        <- igraph::degree(gmmcsub)        ## number of MMC rivals
    cent_deg_all <- V(gmmc)$cent_deg_all
    cent_deg_single <- cent_deg_all - cent_deg_mmc
    #
    cent_deg_mmc_weight[idx.sub.in.mmc] <- sapply(1:vcount(gmmcsub), function(v){
      eids <- unname(unlist(igraph::incident_edges(gmmcsub, v)))
      return(sum(E(gmmcsub)$weight[eids]))
    })
    
    ##-----------------
    ## init period dataframe
    ##------------------
    tdf <- data.frame(name=V(gmmc)$vertex.names, 
                      i=as.integer(V(gmmc)), 
                      year=periods[t], 
                      type=as.integer(V(gmmc)) %in% idx.sub.in.mmc,
                      ipo = V(gmmc)$ipo_status,
                      employee_na_age=V(gmmc)$employee_na_age,
                      sales_na_0_mn=V(gmmc)$sales_na_0_mn,
                      cent_deg_mmc=cent_deg_mmc,
                      cent_deg_mmc_weight=cent_deg_mmc_weight,
                      cent_deg_single=cent_deg_single,
                      cent_deg_all= cent_deg_all,
                      cent_deg_mmc_diff= cent_deg_mmc - cent_deg_single,
                      cent_deg_mmc_weight_diff= cent_deg_mmc_weight - cent_deg_single,
                      acq_cnt_10=V(gmmc)$acq_cnt_10,
                      acq_cnt_5=V(gmmc)$acq_cnt_5,
                      acq_cnt_1=V(gmmc)$acq_cnt_1,
                      acq_sum_1=V(gmmc)$acq_sum_1,
                      acq_mean_1=V(gmmc)$acq_mean_1,
                      acq_median_1=V(gmmc)$acq_median_1,
                      acq_sum_2=V(gmmc)$acq_sum_2,
                      acq_mean_2=V(gmmc)$acq_mean_2,
                      acq_median_2=V(gmmc)$acq_median_2,
                      # pres1=pres1, pres2=pres2, pres3=pres3, pres4=pres4,
                      smmc1=smmc1, 
                      # smmc2=smmc2, smmc3=smmc3, 
                      smmc4=smmc4,
                      # pres1n=pres1n, pres2n=pres2n, pres3n=pres3n, pres4n=pres4n,
                      smmc1n=smmc1n, smmc2n=smmc2n, smmc3n=smmc3n, smmc35n=smmc35n,
                      smmc4n=smmc4n, smmc45n=smmc45n, smmc5n=smmc5n,
                      # smmc5n=smmc5n,
                      # preseig=preseig,  presconstr=presconstr,
                      smmceig=smmceig, smmcconstr=smmcconstr,
                      smmcclose=smmcclose, smmcbetw=smmcbetw,
                      # smmcalphaw=smmcalphaw, smmcalphainvw=smmcalphainvw,
                      smmcdivers=smmcdivers,
                      # presclose=presclose, presbetw=presbetw,
                      # presalphaw=presalphaw, presalphainvw=presalphainvw,
                      # presdivers=presdivers,
                      y.new=NA, y.cur=NA)
    
    ##-----------------
    ## CRUNCHBASE -  DV Acquisitions
    ##-----------------
    cat('  computing DV acquisition counts\n')
    al <- lapply(V(gmmc)$vertex.names, function(.)list(cur=0, new=0))
    names(al) <- V(gmmc)$vertex.names
    ## subset CB acquisitions
    t1 <- sprintf('%d-01-01',periods[t]) ## inclusive start date 'YYYY-MM-DD'
    t2 <- sprintf('%d-12-31',periods[t]) ## inclusive end date 'YYYY-MM-DD'
    ## period years indicate [start, end) -- start inclusive; end exclusive ; acquisitions witin period
    acqs.pd <- cb$co_acq[cb$co_acq$acquired_on >= t1 & cb$co_acq$acquired_on <= t2, ]
    ## acquisitions in period
    for (i in 1:nrow(acqs.pd)) {
      name.t <- acqs.pd$acquiree_name_unique[i]
      name.a <- acqs.pd$acquirer_name_unique[i]
      ## vertex index
      idx.t <- which(V(gmmc)$vertex.names == name.t)
      idx.a <- which(V(gmmc)$vertex.names == name.a)
      if (length(idx.a)==0) {
        if (i %% 100 == 0) cat(sprintf(' skipping [i=%04s] %s-->%s\n',i,name.a,name.t))
        next  ## skip if acquirer isn't in network
      }
      ## indices of firms in markets of acquiring firm
      js <- unname(which(m.ms[idx.a,] > 0))
      ## indices of competitors of acquiring firm
      idx.comp <- which(apply(m.ms, 1, function(x) sum(x[js]) > 0))
      ## increment 
      target.in.net <- length(idx.t)>0  ## target exists in network
      is.cur.market.acq <- ifelse(target.in.net, idx.t %in% idx.comp, FALSE)  ## target in competitors
      if (is.cur.market.acq) {
        al[[name.a]]$cur <- al[[name.a]]$cur + 1
      } else {
        al[[name.a]]$new <- al[[name.a]]$new + 1
      }
    }
    ## ASSIGN DVs
    tdf$y.new <- sapply(al, function(x)x$new)
    tdf$y.cur <- sapply(al, function(x)x$cur)
    
    
    ##-----------------
    ## RavenPack Actions  (alternative Acquisitions DV; 2nd )
    ##-----------------
    cat('  computing RavenPack action counts\n')
    rvnsub <- rvn[rvn$year==periods[t], ] ## ravepack data subset
    n <- vcount(gmmc)
    
    ## loop each action type to set
    V(gmmc)$rp_all <- 0
    for (j in 1:nrow(.rptmp)) {
      val <- getRpActionCount(rvnsub, V(gmmc)$vertex.names, cbrpmap, .rptmp$rp[j])
      ## increment action count to total intensity "rp_all"
      V(gmmc)$rp_all <- V(gmmc)$rp_all + val
      ## set action count to new attributes
      gmmc <- igraph::set.vertex.attribute(gmmc, name=.rptmp$dfreg[j], value=val)
      ##---------------------------------------------------------------
      ## PRESSURE 1: vertex-weighted degree (count of rivals actions)
      ##---------------------------------------------------------------
      wdeg_val <- getNeighborsAttrSum(gmmc, .rptmp$dfreg[j])
      wdeg_name <- sprintf('wdeg_%s',.rptmp$dfreg[j])
      gmmc <- igraph::set.vertex.attribute(gmmc, name=wdeg_name, value=wdeg_val)
      ## TODO: add weighted eigen centrality (using vertex dyad sums as edge weights)
      ##        as alternative measure PRESSURE 2
      ##---------------------------------------------------------------
      ##  Add attributes to Regression Dataframe
      ##---------------------------------------------------------------
      tdf[, .rptmp$dfreg[j] ] <- val
      tdf[, wdeg_name ] <- wdeg_val
      tdf[, 'rp_all'] <- V(gmmc)$rp_all
    }
    
    ## add aggregate action measures
    tdf$rp_NON_acquisitions <- tdf$rp_all - tdf$rp_Acquisitions
    tdf$rp_net_relations <- tdf$rp_New_product + tdf$rp_Market_expansions ## market expansions + product market entires
    tdf$rp_net_invariant <- (tdf$rp_Capacity + tdf$rp_Legal + 
                               tdf$rp_Marketing + tdf$rp_Pricing +
                               tdf$rp_Strategic_alliances) ## excluding acquisitions, new_product, markte_expansions
    
    
    ##-------------------------
    ## append dataframe
    dfl <- rbind(dfl, tdf)   
    
    ##------------------------
    ## save network lists
    glistfile <- sprintf('acq_sys_mmc_panel_RPactions_NETLIST_ego-%s_d%s_%s-%s_t%s.rds', 
                         firm_i_ego, d, periods[1], periods[length(periods)], t)
    saveRDS(list(gt=gt, gmmc=gmmc, gmmcsub=gmmcsub), file=file.path(data_dir,glistfile))
    
    
  } ## end firm_i period loop
  
  
  ## SAVE DATA
  firmfile <- sprintf('acq_sys_mmc_panel_RPactions_regression_df_ego-%s_d%s_%s-%s.csv', 
                      firm_i_ego, d, periods[1], periods[length(periods)])
  write.csv(dfl, file=file.path(data_dir,firmfile))
  
  
  
} ## end firms.todo loop

  
  
  
  






  
  
  
  
  
## READ IN DATA FILE
firmfile <- sprintf('acq_sys_mmc_panel_regression_df_col%s_%s_%s-%s.csv',
                    50, name_i, periods[1], periods[length(periods)])
dfl <- read.csv(file.path(result_dir,firmfile))








##=======================================
##  TEST PANEL DATA RERESSION
##---------------------------------------
dfl$name <- as.character(dfl$name)
yrmin <- min(unique(dfl$year))
## get public firms to subset data
ipo.yrmin.name <- dfl$name[which(dfl$ipo==1)]
# ipo.yrmin.name <- V(gmmc)$vertex.names[ which(V(gmmc)$ipo_status == 1) ]
## get firms by any(type==TRUE)
names.type.true <- unique(dfl$name[dfl$type])
## Regression data subset
dfreg <- dfl[which(dfl$name %in% names.type.true & dfl$name %in% ipo.yrmin.name), ]
# dfreg <- dfl[which(dfl$name %in% names.type.true), ]

## compute feedback measure
dfreg$feedback1 <- NA
dfreg$feedback2 <- NA
dfreg$feedback3 <- NA
dfreg$feedback4 <- NA
dfreg$feedback5 <- NA
for (t in unique(dfreg$year)) {
  if (t == min(unique(dfreg$year))) {
    dfreg$feedback1[which(dfreg$year==t)] <- NA
    dfreg$feedback2[which(dfreg$year==t)] <- NA
    dfreg$feedback3[which(dfreg$year==t)] <- NA
    dfreg$feedback4[which(dfreg$year==t)] <- NA
  } else {
    dfreg$feedback1[which(dfreg$year==t)] <- dfreg$smmc1[which(dfreg$year==t)] - dfreg$smmc1[which(dfreg$year==(t-1))] 
    dfreg$feedback2[which(dfreg$year==t)] <- dfreg$smmc2[which(dfreg$year==t)] - dfreg$smmc2[which(dfreg$year==(t-1))] 
    dfreg$feedback3[which(dfreg$year==t)] <- dfreg$smmc3[which(dfreg$year==t)] - dfreg$smmc3[which(dfreg$year==(t-1))] 
    dfreg$feedback4[which(dfreg$year==t)] <- dfreg$smmc4[which(dfreg$year==t)] - dfreg$smmc4[which(dfreg$year==(t-1))] 
  }
}

##factors
dfreg$i <- as.factor(dfreg$i)
dfreg$year <- as.factor(dfreg$year)
dfreg$y <- dfreg$y.cur + dfreg$y.new



##=======================================
##  GLMER
##---------------------------------------
coef.map <- list(dum.crisis='Financial Crisis Dummy',
                 `I(acq_cnt_5 > 0)TRUE`='Acquisition Experience (5 yr binary)',
                 `I(acq_sum_1/1e+09)`='Prior Yr Acquisition Sum ($ Mn.)',
                 `I(employee_na_age/1e+06)`='Employees',
                 `I(sales_na_0_mn/1e+06)`='Sales ($ Mn.)',
                 `log(1 + cent_deg_all)`='Ln Competitors',
                 `smmc1n`='System MMC',
                 `I(smmc1n^2)`='System MMC Squared',
                 `pres1n`='Competitive Pressure',
                 `lag(feedback1, 0:0)`='System Feedback',
                 `I(smmc1n^2):pres1n`='System MMC * Pressure',
                 `I(smmc1n^2):lag(feedback1, 0:0)`='System MMC * Feedback',
                 acq_sum_1_sc = 'Prior Yr Acquisition Sum ($ Mn.)',
                 employee_na_age_sc = 'Employees',
                 sales_na_0_mn_sc = 'Sales ($ Mn.)',
                 cent_deg_all_sc = 'Competitors'
)

dfreg$dum.crisis <- 0
dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
dfreg$acq_sum_1_sc <- scale(dfreg$acq_sum_1, center = T, scale = T)
dfreg$employee_na_age_sc <- scale(dfreg$employee_na_age, center = T, scale = T)
dfreg$sales_na_0_mn_sc <- scale(dfreg$sales_na_0_mn, center = T, scale = T)
dfreg$cent_deg_all_sc <- scale(dfreg$cent_deg_all, center = T, scale = T)

##-----------------------------------
## NEGBIN
##-----------------------------------
lmn0 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc
                 ,
                 data=dfreg, 
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=1e5)))

lmn1 <- glmer.nb(y ~ (1 | i) + (1 | year) + 
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc +
                   smmc1n + I(smmc1n^2),
                 data=dfreg, 
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=3e5)))
##***********convergence warning: 0.00135935***********
lmn2 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc +
                   smmc1n + I(smmc1n^2) + 
                   pres1n +
                   I(smmc1n^2):pres1n,
                 data=dfreg, 
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=4e5)))
##***********convergence warning: 0.0611716***********
lmn3 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc +
                   smmc1n + I(smmc1n^2) +
                   lag(feedback1, 0:0) + 
                   I(smmc1n^2):lag(feedback1, 0:0),
                 data=dfreg, 
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=4e5)))

lmnAll <- glmer.nb(y ~ (1 | i) + (1 | year) + 
                     dum.crisis + I(acq_cnt_5>0) + 
                     acq_sum_1_sc  + employee_na_age_sc + 
                     sales_na_0_mn_sc + cent_deg_all_sc +
                     smmc1n + I(smmc1n^2) + 
                     pres1n +
                     lag(feedback1, 0:0) + # I(presconstr^2) +
                     I(smmc1n^2):pres1n +
                     I(smmc1n^2):lag(feedback1, 0:0),
                   data=dfreg, 
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=4e5)))

lmn.list <- list(lmn0,lmn1,lmn2,lmn3,lmnAll)
screenreg(lmn.list, digits = 3, 
          custom.coef.map = coef.map)

saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_glmernb_list.rds'))
texreg::htmlreg(lmn.list,
                file.path(result_dir, 'acqmmc_glmernb_list.html'),
                digits = 3, custom.coef.map = coef.map)




##-----------------------------------
## POISSON
##-----------------------------------
lm0 <- glmer(y ~ (1 | i) + (1 | year) +  
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + cent_deg_all_sc
             ,
             data=dfreg, family = poisson,
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=1e5)))

lm1 <- glmer(y ~ (1 | i) + (1 | year) + 
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2),
             data=dfreg, family = poisson,
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=2e5)))

lm2 <- glmer(y ~ (1 | i) + (1 | year) +  
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + log(1+cent_deg_all) +
               smmc1n + I(smmc1n^2) + 
               pres1n +
               I(smmc1n^2):pres1n,
             data=dfreg, family = poisson,
             control=glmerControl(optimizer="bobyqa",
                                  boundary.tol = 1e-4,
                                  optCtrl=list(maxfun=3e5)))

lmpg2 <- pglm(y  ~  
                dum.crisis + I(acq_cnt_5>0) + 
                acq_sum_1_sc  + employee_na_age_sc + 
                sales_na_0_mn_sc + log(1+cent_deg_all) +
                smmc1n + I(smmc1n^2) + 
                pres1n +
                I(smmc1n^2):pres1n,
              data=dfreg, family = poisson, 
              model = 'random', effect = 'twoways',
              R = 100, method='nr',
              index=c('i','year')); summary(lmpg2)

lm3 <- glmer(y ~ (1 | i) + (1 | year) +  
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2) +
               lag(feedback1, 0:0) + 
               I(smmc1n^2):lag(feedback1, 0:0),
             data=dfreg, family = poisson,
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=3e5)))

lmAll <- glmer(y ~ (1 | i) + (1 | year) + 
                 dum.crisis + I(acq_cnt_5>0) + 
                 acq_sum_1_sc  + employee_na_age_sc + 
                 sales_na_0_mn_sc + cent_deg_all_sc +
                 smmc1n + I(smmc1n^2) + 
                 pres1n +
                 lag(feedback1, 0:0) + # I(presconstr^2) +
                 I(smmc1n^2):pres1n +
                 I(smmc1n^2):lag(feedback1, 0:0),
               data=dfreg, family = poisson,
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=3e5)))

lmn.list <- list(lm0,lm1,lm2,lm3,lmAll)
screenreg(lmn.list, digits = 3, 
          custom.coef.map = coef.map)

saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_glmernb_poisson_list.rds'))
texreg::htmlreg(lmn.list,
                file.path(result_dir, 'acqmmc_glmernb_poisson_list.html'),
                digits = 3, custom.coef.map = coef.map)





# pglm(formula, data, subset, na.action,
#      effect = c("individual", "time", "twoways"),
#      model = c("random", "pooling", "within", "between"),
#      family, other = NULL, index = NULL, start = NULL, R = 20,  ...) 

# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#            dfreg,
#            family = negbin, effect = 'twoways', model = "within", 
#            print.level = 3, method = "nr", ## 'bfgs'
#            index = c('i', 'year'), R=20)
# 
# m2p <- pglm(y.new ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "pooling", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4p <- pglm(y.new ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "pooling",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# summary(m2w); summary(m2p); summary(m4w); summary(m4p)
# 
# 
# 
# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2b <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "between", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4b <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "between",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2w); summary(m2b); summary(m4w); summary(m4b)
# 
# 
# 
# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2b <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "between", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4b <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "between",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2w); summary(m2b); summary(m4w); summary(m4b)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# m2x <- pglm(y.new ~ smmc2 + feedback2 + smmc2:feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "random", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2xc <- pglm(y.cur ~ smmc2 + I(smmc2^2) + feedback2 + I(smmc2^2) :feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "random", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2x); summary(m2xc)
# 
# 
# 
# 
# m2x <- pglm(y.new ~ smmc2 + feedback2 + smmc2:feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2xc <- pglm(y.cur ~ smmc2 + I(smmc2^2) + feedback2 + I(smmc2^2) :feedback2,
#              dfreg,
#              family = negbin, effect = 'twoways', model = "within", 
#              print.level = 3, method = "nr", ## 'bfgs'
#              index = c('i', 'year'), R=20)
# summary(m2x); summary(m2xc)
# 
# 
# 
# m4x <- pglm(y.new ~ smmc4 + feedback4 + smmc4:feedback4,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4xc <- pglm(y.cur ~ smmc4 + I(smmc4^2) + feedback4 + I(smmc4^2) :feedback4,
#              dfreg,
#              family = negbin, effect = 'twoways', model = "within", 
#              print.level = 3, method = "nr", ## 'bfgs'
#              index = c('i', 'year'), R=20)
# summary(m4x); summary(m4xc)
# 
# 
# 
# 
# 
# m2yl <- pglm(y ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#             dfreg,
#             family = poisson, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20); summary(m2yl)
# 
# m4yl <- pglm(y ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = poisson, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20); summary(m4yl)

####################################################################
m4y <- pglm(y ~ smmc4 + lag(feedback4, 0) + smmc4:lag(feedback4, 0),
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4y)
####################################################################


m4y2 <- pglm(y ~ smmc4 + I(smmc4^2) + lag(feedback4, 0) + 
               smmc4:lag(feedback4, 0) + I(smmc4^2):lag(feedback4, 0),
             dfreg,
             family = poisson, effect = 'twoways', model = "within", 
             print.level = 3, method = "nr", ## 'bfgs'
             index = c('i', 'year'), R=20); summary(m4y2)

m4yf <- pglm(y ~ smmc4 +  lag(feedback4, 0:1) + lag(I(feedback4^2), 0:1) +
               smmc4:lag(I(feedback4^2), 0:1),
             dfreg,
             family = poisson, effect = 'twoways', model = "within", 
             print.level = 3, method = "nr", ## 'bfgs'
             index = c('i', 'year'), R=20); summary(m4yf)






############## CONTROLS ##################################
m4a <- pglm(y~ smmc1 + pres1 + smmc1:pres1 + 
              acq_cnt_5 +  cent_deg_all,
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4a)
############## //CONTROLS ##################################

m4a <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              # cent_deg_all +
              smmc1n + # I(smmc1n^2) +
              pres1n + # I(presconstr^2) +
              smmc1n:pres1n              # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4a)



m4b <- pglm(y ~  acq_cnt_5 + I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n + # I(smmc1n^2) +
              pres1n + 
              lag(feedback1, 0:1) + # I(presconstr^2) +
              smmc1n:lag(feedback1, 0:1)
            # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4b)


m4c <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n +  I(smmc1n^2) +
              pres1n +  I(pres1n^2) +
              smmc1n:pres1n +
              I(smmc1n^2):I(pres1n^2) 
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4c)




m4d <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n +  I(smmc1n^2) +
              lag(feedback1, 0:1) + 
              I(smmc1n^2):lag(feedback1, 0:1)
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4d)


dfreg$smmc1n.pres1n <- dfreg$smmc1n-dfreg$pres1n
m4e <- pglm(y.new ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + 
              pres1n + # I(presconstr^2) +
              smmc1n:pres1n # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4e)

m5e <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) +
              pres1n + 
              I(smmc1n^2):pres1n 
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m5e)

m6e <- pglm(y.new ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + 
              lag(feedback1, 0:1) + # I(presconstr^2) +
              smmc1n:lag(feedback1, 0:1) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m6e)

dfreg$dum.crisis <- 0
dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
m7e <- pglm(y ~  dum.crisis + I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) + 
              pres1n +
              # lag(feedback1, 0:0) + # I(presconstr^2) +
              I(smmc1n^2):pres1n #+
            # I(smmc1n^2):lag(feedback1, 0:0) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m7e)

m8e <- pglm(y ~  dum.crisis + I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) + 
              pres1n +
              lag(feedback1, 0:0) + # I(presconstr^2) +
              I(smmc1n^2):pres1n +
              I(smmc1n^2):lag(feedback1, 0:0) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m8e)


mx <- m4e
y <- mx$model[,1]
X <- mx$model[,-1]
X$`I(acq_cnt_5 > 0)` <- as.numeric(X$`I(acq_cnt_5 > 0)`)
X$`smmc1.pres1n` <- X$smmc1n * X$pres1n
beta <- mx$estimate

m4e <- as.matrix(X) %*% beta

y <- mx$model$y

##========================================
##  LINEAR
##----------------------------------------

qs.smmc1n <- quantile(X$smmc1n,  c(1/3,.5,2/3))
smmc1n.f <- unname(sapply(X$smmc1n, function(x){
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))
#
qs.pres1n <- quantile(X$pres1n, c(1/3,.5,2/3))
pres1n.f <- unname(sapply(X$pres1n, function(x){
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))

xdf <- data.frame(y=y, smmc1n.f=smmc1n.f, 
                  pres1n.f=pres1n.f,
                  year=as.numeric(as.character(dfreg$year)))
# xdf <- xdf[which(xdf$year < 2012),]

interaction.plot(xdf$smmc1n.f, xdf$pres1n.f, xdf$y,
                 fun = function(x)mean(x, na.rm=T),
                 xlab='System MMC',
                 ylab='Avg. Acquisitions',
                 trace.label = 'Pressure',
                 type = 'b', pch=c(2,16,8))

##========================================
##  Curvilinear
##----------------------------------------

qs.smmc1n <- quantile(X$smmc1n,  c(1/3,.5,2/3))
smmc1n.f <- unname(sapply(X$smmc1n, function(x){
  if(x >= qs.smmc1n[3]) return('3High')
  if(x >= qs.smmc1n[1]) return('2Mid')
  return('1Low')
}))
#
qs.pres1n <- quantile(X$pres1n, c(1/3,.5,2/3))
pres1n.f <- unname(sapply(X$pres1n, function(x){
  # if(x >= qs.pres1n[3]) return('3High')
  # if(x >= qs.pres1n[1]) return('2Mid')
  # return('1Low')
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))

interaction.plot(smmc1n.f, pres1n.f, y,
                 fun = function(x)mean(x, na.rm=T),
                 xlab='System MMC',
                 ylab='Avg. Acquisitions',
                 trace.label = 'Pressure',
                 type = 'b', pch=c(2,16,8))






# ##--------------------------------------------
# ## Patch category similarity
# ##--------------------------------------------
# patch.todo <- c('qualtrics',
#                 'clarabridge','confirmit','medallia','snap-surveys-ltd')
# for (name_i in patch.todo) {
#   cat(sprintf(' patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %n% 'cat_cos_sim' <- acf$.cov.categoryCosineSimilarity(nets[[t]])
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }

##--------------------------------------------
## Patch  sizes variables
##--------------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
for (name_i in patch.todo) {
  cat(sprintf('\n\n patching %s\n', name_i))
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  end <- as.numeric(names(nets)[t])
  start <- end-1
  for (t in 1:length(nets)) {
    nets[[t]] <- acf$setCovariates(nets[[t]], start, end,
                                   covlist=c('employee','sales'), size=si)
  }
  saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
}

sapply(lapply(nets,function(net)net%v%'employee_na_age'), summary)
sapply(lapply(nets,function(net)net%v%'employee_na_catage'), summary)
sapply(lapply(nets,function(net)net%v%'sales_na_0_mn'), summary)


###############
# for (name_i in patch.todo) {
#   cat(sprintf('\n\n patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %v% 'employee_na_0_log' <- log(1 + (nets[[t]] %v% 'employee_na_0'))
#     nets[[t]] %v% 'employee_na_cat_log' <- log(1 + (nets[[t]] %v% 'employee_na_cat'))
#     nets[[t]] %v% 'employee_na_age_log' <- log(1 + (nets[[t]] %v% 'employee_na_age'))
#     nets[[t]] %v% 'employee_na_catage_log' <- log(1 + (nets[[t]] %v% 'employee_na_catage'))
#     nets[[t]] %v% 'employee_na_0_std' <- c(scale(nets[[t]] %v% 'employee_na_0'))
#     nets[[t]] %v% 'employee_na_cat_std' <- c(scale(nets[[t]] %v% 'employee_na_cat'))
#     nets[[t]] %v% 'employee_na_age_std' <-  c(scale(nets[[t]] %v% 'employee_na_age'))
#     nets[[t]] %v% 'employee_na_catage_std' <-  c(scale(nets[[t]] %v% 'employee_na_catage'))
#     #
#     nets[[t]] %v% 'sales_na_0_log' <- log(1 + (nets[[t]] %v% 'sales_na_0'))
#     nets[[t]] %v% 'sales_na_cat_log' <- log(1 + (nets[[t]] %v% 'sales_na_cat'))
#     nets[[t]] %v% 'sales_na_age_log' <- log(1 + (nets[[t]] %v% 'sales_na_age'))
#     nets[[t]] %v% 'sales_na_catage_log' <- log(1 + (nets[[t]] %v% 'sales_na_catage'))
#     nets[[t]] %v% 'sales_na_0_std' <- c(scale(nets[[t]] %v% 'sales_na_0'))
#     nets[[t]] %v% 'sales_na_cat_std' <- c(scale(nets[[t]] %v% 'sales_na_cat'))
#     nets[[t]] %v% 'sales_na_age_std' <-  c(scale(nets[[t]] %v% 'sales_na_age'))
#     nets[[t]] %v% 'sales_na_catage_std' <-  c(scale(nets[[t]] %v% 'sales_na_catage'))
# 
#     nets[[t]] %v% 'sales_na_0_mn' <- (nets[[t]] %v% 'sales_na_0') / 1e6
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }






##--------------------------------------------
## Patch Centrality
##--------------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
for (name_i in patch.todo) {
  cat(sprintf('\n\n patching %s\n', name_i))
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  for (t in 1:length(nets)) {
    nets[[t]] <- acf$.cov.centrality(nets[[t]])
    year <- as.numeric(names(nets)[t])
    nets[[t]] %n% 'shared_investor_nd' <- acf$.cov.sharedInvestor(nets[[t]], ih, cb$co_rou, cb$inv_rou, cb$inv, year, off.diagonal.blocks=F)
  }
  saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
}


##-----------------------------------------
## Check encounters distribution per year
##-----------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
ckl <- list()
for (name_i in patch.todo) {
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  ckl[[name_i]] <- sapply(nets, function(net){
    g <- asIgraph(net)
    g.sub <- igraph::induced.subgraph(g,which(igraph::degree(g)>0))
    deg <- igraph::degree(g.sub)
    return(c(e=ecount(g.sub),v=vcount(g.sub),
             deg_avg=round(mean(deg),1), deg_med=median(deg),
             deg_min=min(deg), deg_max=max(deg)))
  })
  # ec <- sapply(nets, function(net) {
  #     m <- net[,]
  #     return(sum(m[lower.tri(m)])) 
  #   })
  # vc <- sapply(nets, function(net) nrow(net[,]))
  # ckdf <- rbind(ckdf, data.frame(
  #     name=name_i,
  #     yrs=length(nets),
  #     v_sum=sum(vc), 
  #     v_avg=mean(vc),
  #     v_min=min(vc),
  #     v_max=max(vc),
  #     e_sum=sum(ec),
  #     e_avg=mean(ec),
  #     e_min=min(ec),
  #     e_max=max(ec)
  #   ))
}

ckrm <- lapply(ckl, rowMeans)
ckrs <- lapply(ckl, rowSums)
ckrmdf <- t(round(as.data.frame(ckrm),1))
ckrsdf <- t(round(as.data.frame(ckrs),1))

colSums(ckrmdf)

colSums(ckrsdf)

##-----------------------------------------------
## TEST BTERGM
##-----------------------------------------------
# name_i <- 'qualtrics'
# nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
# acf$covSummaryPlot(nets, name_i, net_dir)
library(btergm)
name_i <- 'qualtrics'
d <- 3
nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))

mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor'))
shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))

####################### DEFINE MODELS ###################################

m4 <-   nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed=T) +
  nodematch("ipo_status", diff = F) +
  nodematch("state_code", diff = F) +
  nodecov("age") + absdiff("age") +
  nodecov("employee_na_age") + nodecov("sales_na_0") +
  edgecov(cossim) + edgecov(centjoin) + edgecov(centratio) + edgecov(shcomp) + edgecov(shinv) +
  edgecov(mmc) +
  ##edgecov(cpa) +
  ##edgecov(cpc) +
  ##edgecov(cpp) +
  memory(type = "stability", lag = 1) +
  timecov(transform = function(t) t) +
  nodecov("genidx_multilevel") +
  nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") +
  cycle(3) + cycle(4) + cycle(5)




# ##--------------------------------------------
# ## Patch category similarity 
# ##--------------------------------------------
# patch.todo <- c('qualtrics',
#                 'clarabridge','confirmit','medallia','snap-surveys-ltd')
# for (name_i in patch.todo) {
#   cat(sprintf(' patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %n% 'cat_cos_sim' <- acf$.cov.categoryCosineSimilarity(nets[[t]])
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }


# ### Network size by scope of awareness
# name_i <- 'qualtrics'
# ds <- 1:10
# names(ds) <- as.character(ds)
# gstat <- sapply(ds, function(d){
#   g.base <- g.full  
#   g.d.sub <- igraph::make_ego_graph(g.base, d, V(g.base)[V(g.base)$name==name_i], 'all')[[1]]
#   return(c(v=vcount(g.d.sub),e=ecount(g.d.sub),d=igraph::graph.density(g.d.sub)))
# })
# 
# plot(gstat[1,], log='', type='b')
# plot(c(NA,diff(gstat[1,])), log='', type='b')





## get focal firm distances to specific firms over time
net.yrs <- 2007:2016
gs <- list()
for (yr in net.yrs) {
  gs[[as.character(yr)]] <- asIgraph(readRDS(file.path('firm_nets_rnr',sprintf('facebook_d3_y%d.rds',yr))))
}


g = gs[[length(gs)]]

co.idx <- which(
  grepl('content and publishing',V(g.full)$category_group_list) 
  & !(V(g.full)$employee_count %in% c('NA','-','1-10','11-50','51-100'))
)

####################################################################


##
## run main network period creation loop
##
for (i in 1:length(firms.todo)) {

  name_i <- firms.todo[i]
  cat(sprintf('\n\n------------ %s -------------\n\n',name_i))
  periods <- seq(startYr,endYr,yrpd)
  company.name <- 'company_name_unique'
  g.base <- g.full  
  
  ## focal firm ego network sample
  g.d.sub <- igraph::make_ego_graph(graph = g.base, nodes = V(g.base)[V(g.base)$name==name_i], order = d, mode = 'all')[[1]]
  
  ## convert to network object
  net.d.sub <- asNetwork(g.d.sub)
  net <- net.d.sub
  net %n% 'ego' <- name_i
  
  ##------------------------------------------------------
  ##-------preprocess parent-subsidiary relationships-----
  ##----------node collapse like acquisitions-------------
  ##------------------------------------------------------
  cat(' collapsing parent-subsidiary relations...')
  ## load in manually checked parent-subsidiary relations
  # dfpar.xl <- file.path(sup_data_dir, 'private_firm_financials_all_firm_networks_delete_nodes.xlsx')
  # dfpar <- read_excel(dfpar.xl,  na = c('','-',"'-"))
  dfpar.csv <- file.path(sup_data_dir, 'private_firm_financials_all_firm_networks_delete_nodes.csv')
  dfpar <- read.csv(dfpar.csv, na.strings = c('','-',"'-","--"), stringsAsFactors = F)
  dfpar <- dfpar[!is.na(dfpar$parent), ]
  dfpar$parent <- str_to_lower(dfpar$parent)
  ## merge in parent uuid
  .parent.uuid <- cb$co[,c('company_name_unique','company_uuid')]
  names(.parent.uuid) <- c('parent_name_unique', 'parent_uuid') 
  dfpar <- merge(dfpar[,c('parent','firm')], .parent.uuid, by.x='parent', by.y='parent_name_unique', all.x=T, all.y=F)
  ## merge in firm uuid
  .firm.uuid <- cb$co[,c('company_name_unique','company_uuid')]
  dfpar <- merge(dfpar, .firm.uuid, by.x='firm', by.y='company_name_unique', all.x=T, all.y=F)
  ## CrunchBase parent relationships to collapse
  cb.par <- cb$co_parent[,c('company_name_unique','parent_name_unique','parent_org_uuid','org_uuid')]
  names(cb.par) <- c('firm','parent','parent_uuid','company_uuid')
  ## combine CrunchBase and manual parent-child mappings
  par.chi <- rbind(dfpar, cb.par)
  ## filter both parent, child in full graph
  gfuuid <- V(g.d.sub)$company_uuid
  par.chi <- par.chi[which(!is.na(par.chi$parent_uuid) & par.chi$parent_uuid %in% gfuuid 
                           & !is.na(par.chi$company_uuid) & par.chi$company_uuid %in% gfuuid), ]
  ## merge in founded_on date of parent
  par.chi <- merge(par.chi, cb$co[,c('company_name_unique','founded_on')], by.x='parent', by.y='company_name_unique', all.x=T, all.y=F)
  ## quasi-node collapse subsidiaries to parent nodes
  par.chi.nc <- data.frame(
    acquirer_uuid=par.chi$parent_uuid,
    acquiree_uuid=par.chi$company_uuid,
    acquired_on=par.chi$founded_on,   ## for parent-subsidiary mapping, just use parent company founded_on date
    stringsAsFactors = F
  )
  g.d.sub <- acf$nodeCollapseGraph(g.d.sub, par.chi.nc, remove.isolates=T, verbose = T)
  
  cat('done.\n')
  
  ##_-----------------------------------------------------
  ##-------process pre-start-year acquisitions------------
  ##------------------------------------------------------
  acqs.pd <- cb$co_acq[cb$co_acq$acquired_on <= sprintf('%d-12-31',startYr-1), ]
  g.d.sub <- acf$nodeCollapseGraph(g.d.sub, acqs.pd, remove.isolates=T, verbose = T)
  net.d.sub <- asNetwork(g.d.sub)
  cat(sprintf('v = %d, e = %d\n',vcount(g.d.sub),ecount(g.d.sub)))
  
  # ## subset to firms with employees count > 10
  # idx.employee <- which( !(V(g.d.sub)$employee_count %in% c('NA','-','1-10')) )
  # g.d.sub <- igraph::induced.subgraph(g.d.sub, vids = V(g.d.sub)[idx.employee])
  # cat(sprintf('filtered >10 employee count: v = %d, e = %d\n',vcount(g.d.sub),ecount(g.d.sub)))
  
  
  ##------------Network Time Period List--------------------
  # load('acq_sys_mmc_analysis_test_1.Rdata')
  nl <- list()
  dfl <- data.frame()
  
  for (t in 2:length(periods)) 
  {
    ## period dates
    cat(sprintf('\n\nmaking period %s-%s:\n', periods[t-1],periods[t]))
    t1 <- sprintf('%d-01-01',periods[t-1]) ## inclusive start date 'YYYY-MM-DD'
    t2 <- sprintf('%d-12-31',periods[t-1]) ## inclusive end date 'YYYY-MM-DD'

    ## period years indicate [start, end) -- start inclusive; end exclusive 
    ## acquisitions witin period
    acqs.pd <- cb$co_acq[cb$co_acq$acquired_on >= t1 & cb$co_acq$acquired_on <= t2, ]
    
    ## don't process node collapse acquisitions until after computing MMC and DVs
    ## 1. MOVED TO END OF LOOP
    
    ## 2. Subset Period Network
    nl[[t]] <- acf$makePdNetwork(asNetwork(g.d.sub), periods[t-1], periods[t], isolates.remove = F) 
    
    ## 3. Set Covariates for updated Period Network
    nl[[t]] <- acf$setCovariates(nl[[t]], periods[t-1], periods[t],
                                 covlist = c('age','ipo_status','employee','sales'),
                                 acq=cb$co_acq, br=cb$co_br, ipo=cb$co_ipo, 
                                 rou=cb$co_rou, inv_rou=cb$inv_rou, inv=cb$inv,
                                 coop=sdc, ih=ih, size=si)
    
    ##=================================================================
    ## MMC & Acquisitions
    ##-----------------------------------------------------------------
    gt <- asIgraph(nl[[t]])
    V(gt)$com_multilevel <- igraph::multilevel.community(gt)$membership

    ## Create [Firm x Market] incidence matrix
    ## markets of each firm (which markets they are in based on NC membership of their rivals)
    cat('computing firm-market matrix\n')
    membership <- V(gt)$com_multilevel
    markets <- unique(membership)
    markets <- markets[order(markets)]
    adj <- igraph::as_adjacency_matrix(gt, sparse = F)
    df.ms <- ldply(1:nrow(adj), function(i){
      i.nbr <- unname(which( adj[i, ] == 1 ))  ## rivals (neighbors in the network) ## i.nbr <- as.integer(igraph::neighbors(g, V(g)[i]))  
      i.ms <- unique(membership[i.nbr])  ## markets of firm i (the markets of the rivals of firm i)
      i.ms.row <- rep(0, length(markets)) ## dummy row for [Firm x Market] matrix
      i.ms.row[ i.ms ] <- 1 ## assign [Firm x Market] matrix row value to 1
      names(i.ms.row) <- sapply(1:length(markets),function(x)paste0('m',x))
      return(i.ms.row)
    })
    ## convert df to matrix
    m.ms <- as.matrix(df.ms)
    ## bipartite firm-market
    gb <- igraph::graph_from_incidence_matrix(m.ms, directed = F)
    ## MMC weighted unipartite projections
    gmmc <- igraph::bipartite.projection(gb, multiplicity = T)$proj1
    E(gmmc)$invweight <- 1/E(gmmc)$weight
    ##--------------------
    ##  SET COVARIATES
    ##--------------------
    cat(sprintf('computing controls\n'))
    V(gmmc)$vertex.names <- V(gt)$vertex.names
    V(gmmc)$ipo_status <- V(gt)$ipo_status
    V(gmmc)$employee_na_age <- V(gt)$employee_na_age
    V(gmmc)$sales_na_0_mn <- V(gt)$sales_na_0_mn
    V(gmmc)$cent_deg_all <- igraph::degree(gt)
    #
    a10<- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t-1]-10) & cb$co_acq$acquired_year < periods[t-1]), ]
    a5 <- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t-1]-5) & cb$co_acq$acquired_year < periods[t-1]), ]
    a2 <- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t-1]-2) & cb$co_acq$acquired_year < periods[t-1]), ]
    a1 <- cb$co_acq[which(cb$co_acq$acquired_year >= (periods[t-1]-1) & cb$co_acq$acquired_year < periods[t-1]), ]
    #
    V(gmmc)$acq_cnt_10 <- 0
    V(gmmc)$acq_cnt_5 <- 0
    V(gmmc)$acq_cnt_1 <- 0
    V(gmmc)$acq_sum_1 <- 0
    V(gmmc)$acq_mean_1 <- 0
    V(gmmc)$acq_median_1 <- 0
    V(gmmc)$acq_sum_2 <- 0
    V(gmmc)$acq_mean_2 <- 0
    V(gmmc)$acq_median_2 <- 0
    for (v in 1:vcount(gmmc)) {
      firm_v <- V(gmmc)$vertex.names[v]
      ## ACQUISITION EXPERIENCE (count of acquisitions in last 5 years)
      V(gmmc)$acq_cnt_10[v] <- length(which(a10$acquirer_name_unique==firm_v))
      V(gmmc)$acq_cnt_5[v] <- length(which(a5$acquirer_name_unique==firm_v))
      V(gmmc)$acq_cnt_1[v] <- length(which(a1$acquirer_name_unique==firm_v))
      ## LACK OF SLACK RESOURCES (Value of Acquisitions in previous year)
      vals1 <- a1$price_usd[which(a1$acquirer_name_unique==firm_v)]
      vals2 <- a2$price_usd[which(a2$acquirer_name_unique==firm_v)]
      V(gmmc)$acq_sum_1[v] <- sum(vals1, na.rm = T)
      V(gmmc)$acq_mean_1[v] <-  mean(vals1, na.rm = T)
      V(gmmc)$acq_median_1[v] <-  median(vals1, na.rm = T)
      V(gmmc)$acq_sum_2[v] <-  sum(vals2, na.rm = T)
      V(gmmc)$acq_mean_2[v] <-  mean(vals2, na.rm = T)
      V(gmmc)$acq_median_2[v] <-  median(vals2, na.rm = T)
    }
    
    ##-----------------
    ## COMPUTE Competitive Pressure
    ##-----------------
    pres1 <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmc, exp= 0.1), error = function(e)e)
    pres2 <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmc, exp= 0.2), error = function(e)e)
    pres3 <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmc, exp= 0.3), error = function(e)e)
    pres4 <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmc, exp= 0.4), error = function(e)e)
    pres1n <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmc, exp=-0.1), error = function(e)e)
    pres2n <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmc, exp=-0.2), error = function(e)e)
    pres3n <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmc, exp=-0.3), error = function(e)e)
    pres4n <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmc, exp=-0.4), error = function(e)e)
    preseig <- igraph::eigen_centrality(gmmc)$vector
    presconstr <- igraph::constraint(gmmc)
    presclose <- igraph::closeness(gmmc, normalized = T)
    presbetw  <- igraph::betweenness(gmmc, normalized = T)
    # presalphaw  <- igraph::alpha.centrality(gmmc, weights = 'weight', sparse = T)
    # presalphainvw  <- igraph::alpha.centrality(gmmc, weights = 'invweight', sparse = T)
    presauthw <- igraph::authority.score(gmmc)$vector
    presdivers <- igraph::diversity(gmmc)
    ##-----------------
    ## COMPUTE SYSTEM MMC VECTORS
    ##-----------------
    cat('computing system MMC measure\n')
    ## remove edges weight < 2
    gmmc.e <- igraph::delete.edges(gmmc, which(E(gmmc)$weight < 2))
    ## create subgraph removing isolates
    gmmcsub <- igraph::induced.subgraph(gmmc.e, which(igraph::degree(gmmc.e) > 0))
    E(gmmcsub)$invweight <- 1/E(gmmcsub)$weight
    mmcnames <- V(gmmcsub)$vertex.names
    ## Centrality -- SYSTEMS MMC for only firms with MMC eges
    pc0.1 <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmcsub, exp= 0.1), error = function(e)e)
    pc0.2 <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmcsub, exp= 0.2), error = function(e)e)
    pc0.3 <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmcsub, exp= 0.3), error = function(e)e)
    pc0.4 <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmcsub, exp= 0.4), error = function(e)e)
    #
    pcn0.1 <- tryCatch(tmpn0.1<- igraph::power_centrality(gmmcsub, exp=-0.1), error = function(e)e)
    pcn0.2 <- tryCatch(tmpn0.2<- igraph::power_centrality(gmmcsub, exp=-0.2), error = function(e)e)
    pcn0.3 <- tryCatch(tmpn0.3<- igraph::power_centrality(gmmcsub, exp=-0.3), error = function(e)e)
    pcn0.4 <- tryCatch(tmpn0.4<- igraph::power_centrality(gmmcsub, exp=-0.4), error = function(e)e)
    ## init vector of system MMC for ALL singlemarket and multimarket firms
    smmc1n <- rep(0, vcount(gmmc))
    smmc2n <- rep(0, vcount(gmmc))
    smmc3n <- rep(0, vcount(gmmc))
    smmc4n <- rep(0, vcount(gmmc))
    smmc1 <- rep(0, vcount(gmmc))
    smmc2 <- rep(0, vcount(gmmc))
    smmc3 <- rep(0, vcount(gmmc))
    smmc4 <- rep(0, vcount(gmmc))
    smmceig <- rep(0, vcount(gmmc))
    smmcconstr <- rep(0, vcount(gmmc))
    smmcclose <- rep(0, vcount(gmmc))
    smmcbetw <- rep(0, vcount(gmmc))
    # smmcalphaw<- rep(0, vcount(gmmc))
    # smmcalphainvw<- rep(0, vcount(gmmc))
    smmcauthw <- rep(0, vcount(gmmc))
    smmcdivers <- rep(0, vcount(gmmc))
    ## indices of mmc firm subset within full graph
    idx.sub.in.mmc <- unname(sapply(mmcnames, function(x) which(V(gmmc)$vertex.names == x)))
    ## assign subset MMC firms' system MMC into vector for all firms
    smmc1n[idx.sub.in.mmc] <- pcn0.1
    smmc2n[idx.sub.in.mmc] <- pcn0.2
    smmc3n[idx.sub.in.mmc] <- pcn0.3
    smmc4n[idx.sub.in.mmc] <- pcn0.4
    #
    smmc1[idx.sub.in.mmc] <- pc0.1
    smmc2[idx.sub.in.mmc] <- pc0.2
    smmc3[idx.sub.in.mmc] <- pc0.3
    smmc4[idx.sub.in.mmc] <- pc0.4
    #
    smmceig[idx.sub.in.mmc] <- igraph::eigen_centrality(gmmcsub)$vector
    smmcconstr[idx.sub.in.mmc] <- igraph::constraint(gmmcsub)
    #
    smmcclose[idx.sub.in.mmc] <- igraph::closeness(gmmcsub, normalized = T)
    smmcbetw[idx.sub.in.mmc]  <- igraph::betweenness(gmmcsub, normalized = T)
    # smmcalphaw[idx.sub.in.mmc]  <- igraph::alpha.centrality(gmmcsub, weights = 'weight', sparse = T)
    # smmcalphainvw[idx.sub.in.mmc]  <- igraph::alpha.centrality(gmmcsub, weights = 'invweight', sparse = T)
    #
    smmcauthw[idx.sub.in.mmc] <- igraph::authority.score(gmmcsub)$vector
    smmcdivers[idx.sub.in.mmc] <- igraph::diversity(gmmcsub)
    ##-----------------
    ## Number of competitors
    ##-----------------
    cent_deg_mmc        <- rep(0, vcount(gmmc))
    cent_deg_mmc_weight <- rep(0, vcount(gmmc))
    # 
    cent_deg_mmc[idx.sub.in.mmc]        <- igraph::degree(gmmcsub)        ## number of MMC rivals
    cent_deg_all <- V(gmmc)$cent_deg_all
    cent_deg_single <- cent_deg_all - cent_deg_mmc
    #
    cent_deg_mmc_weight[idx.sub.in.mmc] <- sapply(1:vcount(gmmcsub), function(v){
        eids <- unname(unlist(igraph::incident_edges(gmmcsub, v)))
        return(sum(E(gmmcsub)$weight[eids]))
      })
      
    ##-----------------
    ## init period dataframe
    ##------------------
    tdf <- data.frame(name=V(gmmc)$vertex.names, i=as.integer(V(gmmc)), 
                      year=periods[t-1], type=as.integer(V(gmmc)) %in% idx.sub.in.mmc,
                      ipo = V(gmmc)$ipo_status,
                      employee_na_age=V(gmmc)$employee_na_age,
                      sales_na_0_mn=V(gmmc)$sales_na_0_mn,
                      cent_deg_mmc=cent_deg_mmc,
                      cent_deg_mmc_weight=cent_deg_mmc_weight,
                      cent_deg_single=cent_deg_single,
                      cent_deg_all= cent_deg_all,
                      cent_deg_mmc_diff= cent_deg_mmc - cent_deg_single,
                      cent_deg_mmc_weight_diff= cent_deg_mmc_weight - cent_deg_single,
                      acq_cnt_10=V(gmmc)$acq_cnt_10,
                      acq_cnt_5=V(gmmc)$acq_cnt_5,
                      acq_cnt_1=V(gmmc)$acq_cnt_1,
                      acq_sum_1=V(gmmc)$acq_sum_1,
                      acq_mean_1=V(gmmc)$acq_mean_1,
                      acq_median_1=V(gmmc)$acq_median_1,
                      acq_sum_2=V(gmmc)$acq_sum_2,
                      acq_mean_2=V(gmmc)$acq_mean_2,
                      acq_median_2=V(gmmc)$acq_median_2,
                      pres1=pres1, pres2=pres2, pres3=pres3, pres4=pres4,
                      smmc1=smmc1, smmc2=smmc2, smmc3=smmc3, smmc4=smmc4,
                      pres1n=pres1n, pres2n=pres2n, pres3n=pres3n, pres4n=pres4n,
                      smmc1n=smmc1n, smmc2n=smmc2n, smmc3n=smmc3n, smmc4n=smmc4n,
                      preseig=preseig,  presconstr=presconstr,
                      smmceig=smmceig, smmcconstr=smmcconstr,
                      smmcclose=smmcclose, smmcbetw=smmcbetw,
                      # smmcalphaw=smmcalphaw, smmcalphainvw=smmcalphainvw,
                      smmcdivers=smmcdivers,
                      presclose=presclose, presbetw=presbetw,
                      # presalphaw=presalphaw, presalphainvw=presalphainvw,
                      presdivers=presdivers,
                      y.new=NA, y.cur=NA)
    ##-----------------
    ## COMPUTE DV acquisition counts
    ##-----------------
    cat('computing DV acquisition counts\n')
    al <- lapply(V(gmmc)$vertex.names, function(.)list(cur=0, new=0))
    names(al) <- V(gmmc)$vertex.names
    ## acquisitions in period
    for (i in 1:nrow(acqs.pd)) {
      name.t <- acqs.pd$acquiree_name_unique[i]
      name.a <- acqs.pd$acquirer_name_unique[i]
      ## vertex index
      idx.t <- which(V(gmmc)$vertex.names == name.t)
      idx.a <- which(V(gmmc)$vertex.names == name.a)
      if (length(idx.a)==0) {
        if (i %% 100 == 0) cat(sprintf(' skipping [i=%04s] %s-->%s\n',i,name.a,name.t))
        next  ## skip if acquirer isn't in network
      }
      ## indices of firms in markets of acquiring firm
      js <- unname(which(m.ms[idx.a,] > 0))
      ## indices of competitors of acquiring firm
      idx.comp <- which(apply(m.ms, 1, function(x) sum(x[js]) > 0))
      ## increment 
      target.in.net <- length(idx.t)>0  ## target exists in network
      is.cur.market.acq <- ifelse(target.in.net, idx.t %in% idx.comp, FALSE)  ## target in competitors
      if (is.cur.market.acq) {
        al[[name.a]]$cur <- al[[name.a]]$cur + 1
      } else {
        al[[name.a]]$new <- al[[name.a]]$new + 1
      }
    }
    ## ASSIGN DVs
    tdf$y.new <- sapply(al, function(x)x$new)
    tdf$y.cur <- sapply(al, function(x)x$cur)
    
    ## append dataframe
    dfl <- rbind(dfl, tdf)   
    
    ##+==============================================
    ## 1. Node Collapse acquisitions within period
    ##-----------------------------------------------
    g.d.sub <- acf$nodeCollapseGraph(g.d.sub, acqs.pd, verbose = T)

  }
  
  ## SAVE DATA
  firmfile <- sprintf('acq_sys_mmc_panel_regression_df_col%s_%s_%s-%s.csv', 
                      ncol(dfl), name_i, periods[1], periods[length(periods)])
  write.csv(dfl, file=file.path(result_dir,firmfile))
  
}



## READ IN DATA FILE
firmfile <- sprintf('acq_sys_mmc_panel_regression_df_col%s_%s_%s-%s.csv',
                    50, name_i, periods[1], periods[length(periods)])
dfl <- read.csv(file.path(result_dir,firmfile))








##=======================================
##  TEST PANEL DATA RERESSION
##---------------------------------------
dfl$name <- as.character(dfl$name)
yrmin <- min(unique(dfl$year))
## get public firms to subset data
ipo.yrmin.name <- dfl$name[which(dfl$ipo==1)]
# ipo.yrmin.name <- V(gmmc)$vertex.names[ which(V(gmmc)$ipo_status == 1) ]
## get firms by any(type==TRUE)
names.type.true <- unique(dfl$name[dfl$type])
## Regression data subset
dfreg <- dfl[which(dfl$name %in% names.type.true & dfl$name %in% ipo.yrmin.name), ]
# dfreg <- dfl[which(dfl$name %in% names.type.true), ]

## compute feedback measure
dfreg$feedback1 <- NA
dfreg$feedback2 <- NA
dfreg$feedback3 <- NA
dfreg$feedback4 <- NA
dfreg$feedback5 <- NA
for (t in unique(dfreg$year)) {
  if (t == min(unique(dfreg$year))) {
    dfreg$feedback1[which(dfreg$year==t)] <- NA
    dfreg$feedback2[which(dfreg$year==t)] <- NA
    dfreg$feedback3[which(dfreg$year==t)] <- NA
    dfreg$feedback4[which(dfreg$year==t)] <- NA
  } else {
    dfreg$feedback1[which(dfreg$year==t)] <- dfreg$smmc1[which(dfreg$year==t)] - dfreg$smmc1[which(dfreg$year==(t-1))] 
    dfreg$feedback2[which(dfreg$year==t)] <- dfreg$smmc2[which(dfreg$year==t)] - dfreg$smmc2[which(dfreg$year==(t-1))] 
    dfreg$feedback3[which(dfreg$year==t)] <- dfreg$smmc3[which(dfreg$year==t)] - dfreg$smmc3[which(dfreg$year==(t-1))] 
    dfreg$feedback4[which(dfreg$year==t)] <- dfreg$smmc4[which(dfreg$year==t)] - dfreg$smmc4[which(dfreg$year==(t-1))] 
  }
}
  
##factors
dfreg$i <- as.factor(dfreg$i)
dfreg$year <- as.factor(dfreg$year)
dfreg$y <- dfreg$y.cur + dfreg$y.new



##=======================================
##  GLMER
##---------------------------------------
coef.map <- list(dum.crisis='Financial Crisis Dummy',
                 `I(acq_cnt_5 > 0)TRUE`='Acquisition Experience (5 yr binary)',
                 `I(acq_sum_1/1e+09)`='Prior Yr Acquisition Sum ($ Mn.)',
                 `I(employee_na_age/1e+06)`='Employees',
                 `I(sales_na_0_mn/1e+06)`='Sales ($ Mn.)',
                 `log(1 + cent_deg_all)`='Ln Competitors',
                 `smmc1n`='System MMC',
                 `I(smmc1n^2)`='System MMC Squared',
                 `pres1n`='Competitive Pressure',
                 `lag(feedback1, 0:0)`='System Feedback',
                 `I(smmc1n^2):pres1n`='System MMC * Pressure',
                 `I(smmc1n^2):lag(feedback1, 0:0)`='System MMC * Feedback',
                 acq_sum_1_sc = 'Prior Yr Acquisition Sum ($ Mn.)',
                 employee_na_age_sc = 'Employees',
                 sales_na_0_mn_sc = 'Sales ($ Mn.)',
                 cent_deg_all_sc = 'Competitors'
                 )

dfreg$dum.crisis <- 0
dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
dfreg$acq_sum_1_sc <- scale(dfreg$acq_sum_1, center = T, scale = T)
dfreg$employee_na_age_sc <- scale(dfreg$employee_na_age, center = T, scale = T)
dfreg$sales_na_0_mn_sc <- scale(dfreg$sales_na_0_mn, center = T, scale = T)
dfreg$cent_deg_all_sc <- scale(dfreg$cent_deg_all, center = T, scale = T)

##-----------------------------------
## NEGBIN
##-----------------------------------
lmn0 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + cent_deg_all_sc
               ,
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=1e5)))

lmn1 <- glmer.nb(y ~ (1 | i) + (1 | year) + 
                dum.crisis + I(acq_cnt_5>0) + 
                  acq_sum_1_sc  + employee_na_age_sc + 
                  sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2),
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=3e5)))
##***********convergence warning: 0.00135935***********
lmn2 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                dum.crisis + I(acq_cnt_5>0) + 
                  acq_sum_1_sc  + employee_na_age_sc + 
                  sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2) + 
               pres1n +
               I(smmc1n^2):pres1n,
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=4e5)))
##***********convergence warning: 0.0611716***********
lmn3 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                dum.crisis + I(acq_cnt_5>0) + 
                  acq_sum_1_sc  + employee_na_age_sc + 
                  sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2) +
               lag(feedback1, 0:0) + 
               I(smmc1n^2):lag(feedback1, 0:0),
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=4e5)))

lmnAll <- glmer.nb(y ~ (1 | i) + (1 | year) + 
                  dum.crisis + I(acq_cnt_5>0) + 
                    acq_sum_1_sc  + employee_na_age_sc + 
                    sales_na_0_mn_sc + cent_deg_all_sc +
                 smmc1n + I(smmc1n^2) + 
                 pres1n +
                 lag(feedback1, 0:0) + # I(presconstr^2) +
                 I(smmc1n^2):pres1n +
                 I(smmc1n^2):lag(feedback1, 0:0),
               data=dfreg, 
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=4e5)))

lmn.list <- list(lmn0,lmn1,lmn2,lmn3,lmnAll)
screenreg(lmn.list, digits = 3, 
          custom.coef.map = coef.map)

saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_glmernb_list.rds'))
texreg::htmlreg(lmn.list,
                file.path(result_dir, 'acqmmc_glmernb_list.html'),
                digits = 3, custom.coef.map = coef.map)




##-----------------------------------
## POISSON
##-----------------------------------
lm0 <- glmer(y ~ (1 | i) + (1 | year) +  
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc
                 ,
                 data=dfreg, family = poisson,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=1e5)))

lm1 <- glmer(y ~ (1 | i) + (1 | year) + 
                  dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc +
                   smmc1n + I(smmc1n^2),
                 data=dfreg, family = poisson,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

lm2 <- glmer(y ~ (1 | i) + (1 | year) +  
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + log(1+cent_deg_all) +
                   smmc1n + I(smmc1n^2) + 
                   pres1n +
                   I(smmc1n^2):pres1n,
                 data=dfreg, family = poisson,
                 control=glmerControl(optimizer="bobyqa",
                                      boundary.tol = 1e-4,
                                      optCtrl=list(maxfun=3e5)))

lmpg2 <- pglm(y  ~  
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + log(1+cent_deg_all) +
               smmc1n + I(smmc1n^2) + 
               pres1n +
               I(smmc1n^2):pres1n,
             data=dfreg, family = poisson, 
             model = 'random', effect = 'twoways',
             R = 100, method='nr',
             index=c('i','year')); summary(lmpg2)

lm3 <- glmer(y ~ (1 | i) + (1 | year) +  
                   dum.crisis + I(acq_cnt_5>0) + 
                   acq_sum_1_sc  + employee_na_age_sc + 
                   sales_na_0_mn_sc + cent_deg_all_sc +
                   smmc1n + I(smmc1n^2) +
                   lag(feedback1, 0:0) + 
                   I(smmc1n^2):lag(feedback1, 0:0),
                 data=dfreg, family = poisson,
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=3e5)))

lmAll <- glmer(y ~ (1 | i) + (1 | year) + 
                     dum.crisis + I(acq_cnt_5>0) + 
                     acq_sum_1_sc  + employee_na_age_sc + 
                     sales_na_0_mn_sc + cent_deg_all_sc +
                     smmc1n + I(smmc1n^2) + 
                     pres1n +
                     lag(feedback1, 0:0) + # I(presconstr^2) +
                     I(smmc1n^2):pres1n +
                     I(smmc1n^2):lag(feedback1, 0:0),
                   data=dfreg, family = poisson,
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=3e5)))

lmn.list <- list(lm0,lm1,lm2,lm3,lmAll)
screenreg(lmn.list, digits = 3, 
          custom.coef.map = coef.map)

saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_glmernb_poisson_list.rds'))
texreg::htmlreg(lmn.list,
                file.path(result_dir, 'acqmmc_glmernb_poisson_list.html'),
                digits = 3, custom.coef.map = coef.map)





# pglm(formula, data, subset, na.action,
#      effect = c("individual", "time", "twoways"),
#      model = c("random", "pooling", "within", "between"),
#      family, other = NULL, index = NULL, start = NULL, R = 20,  ...) 

# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#            dfreg,
#            family = negbin, effect = 'twoways', model = "within", 
#            print.level = 3, method = "nr", ## 'bfgs'
#            index = c('i', 'year'), R=20)
# 
# m2p <- pglm(y.new ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "pooling", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4p <- pglm(y.new ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "pooling",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# summary(m2w); summary(m2p); summary(m4w); summary(m4p)
# 
# 
# 
# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2b <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "between", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4b <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "between",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2w); summary(m2b); summary(m4w); summary(m4b)
# 
# 
# 
# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2b <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "between", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4b <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "between",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2w); summary(m2b); summary(m4w); summary(m4b)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# m2x <- pglm(y.new ~ smmc2 + feedback2 + smmc2:feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "random", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2xc <- pglm(y.cur ~ smmc2 + I(smmc2^2) + feedback2 + I(smmc2^2) :feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "random", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2x); summary(m2xc)
# 
# 
# 
# 
# m2x <- pglm(y.new ~ smmc2 + feedback2 + smmc2:feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2xc <- pglm(y.cur ~ smmc2 + I(smmc2^2) + feedback2 + I(smmc2^2) :feedback2,
#              dfreg,
#              family = negbin, effect = 'twoways', model = "within", 
#              print.level = 3, method = "nr", ## 'bfgs'
#              index = c('i', 'year'), R=20)
# summary(m2x); summary(m2xc)
# 
# 
# 
# m4x <- pglm(y.new ~ smmc4 + feedback4 + smmc4:feedback4,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4xc <- pglm(y.cur ~ smmc4 + I(smmc4^2) + feedback4 + I(smmc4^2) :feedback4,
#              dfreg,
#              family = negbin, effect = 'twoways', model = "within", 
#              print.level = 3, method = "nr", ## 'bfgs'
#              index = c('i', 'year'), R=20)
# summary(m4x); summary(m4xc)
# 
# 
# 
# 
# 
# m2yl <- pglm(y ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#             dfreg,
#             family = poisson, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20); summary(m2yl)
# 
# m4yl <- pglm(y ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = poisson, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20); summary(m4yl)

####################################################################
m4y <- pglm(y ~ smmc4 + lag(feedback4, 0) + smmc4:lag(feedback4, 0),
             dfreg,
             family = poisson, effect = 'twoways', model = "within", 
             print.level = 3, method = "nr", ## 'bfgs'
             index = c('i', 'year'), R=20); summary(m4y)
####################################################################


m4y2 <- pglm(y ~ smmc4 + I(smmc4^2) + lag(feedback4, 0) + 
               smmc4:lag(feedback4, 0) + I(smmc4^2):lag(feedback4, 0),
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4y2)

m4yf <- pglm(y ~ smmc4 +  lag(feedback4, 0:1) + lag(I(feedback4^2), 0:1) +
               smmc4:lag(I(feedback4^2), 0:1),
             dfreg,
             family = poisson, effect = 'twoways', model = "within", 
             print.level = 3, method = "nr", ## 'bfgs'
             index = c('i', 'year'), R=20); summary(m4yf)






############## CONTROLS ##################################
m4a <- pglm(y~ smmc1 + pres1 + smmc1:pres1 + 
              acq_cnt_5 +  cent_deg_all,
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4a)
############## //CONTROLS ##################################

m4a <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              # cent_deg_all +
              smmc1n + # I(smmc1n^2) +
              pres1n + # I(presconstr^2) +
              smmc1n:pres1n              # I(smmc1n^2):I(presconstr^2) +
              ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4a)



m4b <- pglm(y ~  acq_cnt_5 + I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n + # I(smmc1n^2) +
              pres1n + 
              lag(feedback1, 0:1) + # I(presconstr^2) +
              smmc1n:lag(feedback1, 0:1)
            # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4b)


m4c <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n +  I(smmc1n^2) +
              pres1n +  I(pres1n^2) +
              smmc1n:pres1n +
              I(smmc1n^2):I(pres1n^2) 
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4c)




m4d <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n +  I(smmc1n^2) +
              lag(feedback1, 0:1) + 
              I(smmc1n^2):lag(feedback1, 0:1)
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4d)


dfreg$smmc1n.pres1n <- dfreg$smmc1n-dfreg$pres1n
m4e <- pglm(y.new ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + 
              pres1n + # I(presconstr^2) +
              smmc1n:pres1n # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4e)

m5e <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) +
              pres1n + 
              I(smmc1n^2):pres1n 
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m5e)

m6e <- pglm(y.new ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + 
              lag(feedback1, 0:1) + # I(presconstr^2) +
              smmc1n:lag(feedback1, 0:1) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m6e)

dfreg$dum.crisis <- 0
dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
m7e <- pglm(y ~  dum.crisis + I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) + 
              pres1n +
              # lag(feedback1, 0:0) + # I(presconstr^2) +
              I(smmc1n^2):pres1n #+
              # I(smmc1n^2):lag(feedback1, 0:0) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m7e)

m8e <- pglm(y ~  dum.crisis + I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) + 
              pres1n +
              lag(feedback1, 0:0) + # I(presconstr^2) +
              I(smmc1n^2):pres1n +
              I(smmc1n^2):lag(feedback1, 0:0) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m8e)


mx <- m4e
y <- mx$model[,1]
X <- mx$model[,-1]
X$`I(acq_cnt_5 > 0)` <- as.numeric(X$`I(acq_cnt_5 > 0)`)
X$`smmc1.pres1n` <- X$smmc1n * X$pres1n
beta <- mx$estimate

m4e <- as.matrix(X) %*% beta

y <- mx$model$y

##========================================
##  LINEAR
##----------------------------------------

qs.smmc1n <- quantile(X$smmc1n,  c(1/3,.5,2/3))
smmc1n.f <- unname(sapply(X$smmc1n, function(x){
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))
#
qs.pres1n <- quantile(X$pres1n, c(1/3,.5,2/3))
pres1n.f <- unname(sapply(X$pres1n, function(x){
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))

xdf <- data.frame(y=y, smmc1n.f=smmc1n.f, 
                  pres1n.f=pres1n.f,
                  year=as.numeric(as.character(dfreg$year)))
# xdf <- xdf[which(xdf$year < 2012),]

interaction.plot(xdf$smmc1n.f, xdf$pres1n.f, xdf$y,
                 fun = function(x)mean(x, na.rm=T),
                 xlab='System MMC',
                 ylab='Avg. Acquisitions',
                 trace.label = 'Pressure',
                 type = 'b', pch=c(2,16,8))

##========================================
##  Curvilinear
##----------------------------------------

qs.smmc1n <- quantile(X$smmc1n,  c(1/3,.5,2/3))
smmc1n.f <- unname(sapply(X$smmc1n, function(x){
  if(x >= qs.smmc1n[3]) return('3High')
  if(x >= qs.smmc1n[1]) return('2Mid')
  return('1Low')
  }))
#
qs.pres1n <- quantile(X$pres1n, c(1/3,.5,2/3))
pres1n.f <- unname(sapply(X$pres1n, function(x){
    # if(x >= qs.pres1n[3]) return('3High')
    # if(x >= qs.pres1n[1]) return('2Mid')
    # return('1Low')
    if(x >= qs.pres1n[2]) return('2High')
    return('1Low')
  }))

interaction.plot(smmc1n.f, pres1n.f, y,
                 fun = function(x)mean(x, na.rm=T),
                 xlab='System MMC',
                 ylab='Avg. Acquisitions',
                 trace.label = 'Pressure',
                 type = 'b', pch=c(2,16,8))






# ##--------------------------------------------
# ## Patch category similarity
# ##--------------------------------------------
# patch.todo <- c('qualtrics',
#                 'clarabridge','confirmit','medallia','snap-surveys-ltd')
# for (name_i in patch.todo) {
#   cat(sprintf(' patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %n% 'cat_cos_sim' <- acf$.cov.categoryCosineSimilarity(nets[[t]])
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }

##--------------------------------------------
## Patch  sizes variables
##--------------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
for (name_i in patch.todo) {
  cat(sprintf('\n\n patching %s\n', name_i))
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  end <- as.numeric(names(nets)[t])
  start <- end-1
  for (t in 1:length(nets)) {
    nets[[t]] <- acf$setCovariates(nets[[t]], start, end,
                                  covlist=c('employee','sales'), size=si)
  }
  saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
}

sapply(lapply(nets,function(net)net%v%'employee_na_age'), summary)
sapply(lapply(nets,function(net)net%v%'employee_na_catage'), summary)
sapply(lapply(nets,function(net)net%v%'sales_na_0_mn'), summary)


###############
# for (name_i in patch.todo) {
#   cat(sprintf('\n\n patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %v% 'employee_na_0_log' <- log(1 + (nets[[t]] %v% 'employee_na_0'))
#     nets[[t]] %v% 'employee_na_cat_log' <- log(1 + (nets[[t]] %v% 'employee_na_cat'))
#     nets[[t]] %v% 'employee_na_age_log' <- log(1 + (nets[[t]] %v% 'employee_na_age'))
#     nets[[t]] %v% 'employee_na_catage_log' <- log(1 + (nets[[t]] %v% 'employee_na_catage'))
#     nets[[t]] %v% 'employee_na_0_std' <- c(scale(nets[[t]] %v% 'employee_na_0'))
#     nets[[t]] %v% 'employee_na_cat_std' <- c(scale(nets[[t]] %v% 'employee_na_cat'))
#     nets[[t]] %v% 'employee_na_age_std' <-  c(scale(nets[[t]] %v% 'employee_na_age'))
#     nets[[t]] %v% 'employee_na_catage_std' <-  c(scale(nets[[t]] %v% 'employee_na_catage'))
#     #
#     nets[[t]] %v% 'sales_na_0_log' <- log(1 + (nets[[t]] %v% 'sales_na_0'))
#     nets[[t]] %v% 'sales_na_cat_log' <- log(1 + (nets[[t]] %v% 'sales_na_cat'))
#     nets[[t]] %v% 'sales_na_age_log' <- log(1 + (nets[[t]] %v% 'sales_na_age'))
#     nets[[t]] %v% 'sales_na_catage_log' <- log(1 + (nets[[t]] %v% 'sales_na_catage'))
#     nets[[t]] %v% 'sales_na_0_std' <- c(scale(nets[[t]] %v% 'sales_na_0'))
#     nets[[t]] %v% 'sales_na_cat_std' <- c(scale(nets[[t]] %v% 'sales_na_cat'))
#     nets[[t]] %v% 'sales_na_age_std' <-  c(scale(nets[[t]] %v% 'sales_na_age'))
#     nets[[t]] %v% 'sales_na_catage_std' <-  c(scale(nets[[t]] %v% 'sales_na_catage'))
# 
#     nets[[t]] %v% 'sales_na_0_mn' <- (nets[[t]] %v% 'sales_na_0') / 1e6
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }






##--------------------------------------------
## Patch Centrality
##--------------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
for (name_i in patch.todo) {
  cat(sprintf('\n\n patching %s\n', name_i))
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  for (t in 1:length(nets)) {
    nets[[t]] <- acf$.cov.centrality(nets[[t]])
    year <- as.numeric(names(nets)[t])
    nets[[t]] %n% 'shared_investor_nd' <- acf$.cov.sharedInvestor(nets[[t]], ih, cb$co_rou, cb$inv_rou, cb$inv, year, off.diagonal.blocks=F)
  }
  saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
}


##-----------------------------------------
## Check encounters distribution per year
##-----------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
ckl <- list()
for (name_i in patch.todo) {
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  ckl[[name_i]] <- sapply(nets, function(net){
      g <- asIgraph(net)
      g.sub <- igraph::induced.subgraph(g,which(igraph::degree(g)>0))
      deg <- igraph::degree(g.sub)
      return(c(e=ecount(g.sub),v=vcount(g.sub),
               deg_avg=round(mean(deg),1), deg_med=median(deg),
               deg_min=min(deg), deg_max=max(deg)))
    })
  # ec <- sapply(nets, function(net) {
  #     m <- net[,]
  #     return(sum(m[lower.tri(m)])) 
  #   })
  # vc <- sapply(nets, function(net) nrow(net[,]))
  # ckdf <- rbind(ckdf, data.frame(
  #     name=name_i,
  #     yrs=length(nets),
  #     v_sum=sum(vc), 
  #     v_avg=mean(vc),
  #     v_min=min(vc),
  #     v_max=max(vc),
  #     e_sum=sum(ec),
  #     e_avg=mean(ec),
  #     e_min=min(ec),
  #     e_max=max(ec)
  #   ))
}

ckrm <- lapply(ckl, rowMeans)
ckrs <- lapply(ckl, rowSums)
ckrmdf <- t(round(as.data.frame(ckrm),1))
ckrsdf <- t(round(as.data.frame(ckrs),1))

colSums(ckrmdf)

colSums(ckrsdf)

##-----------------------------------------------
## TEST BTERGM
##-----------------------------------------------
# name_i <- 'qualtrics'
# nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
# acf$covSummaryPlot(nets, name_i, net_dir)
library(btergm)
name_i <- 'qualtrics'
d <- 3
nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))

mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor'))
shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))

####################### DEFINE MODELS ###################################

m4 <-   nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed=T) +
  nodematch("ipo_status", diff = F) +
  nodematch("state_code", diff = F) +
  nodecov("age") + absdiff("age") +
  nodecov("employee_na_age") + nodecov("sales_na_0") +
  edgecov(cossim) + edgecov(centjoin) + edgecov(centratio) + edgecov(shcomp) + edgecov(shinv) +
  edgecov(mmc) +
  ##edgecov(cpa) +
  ##edgecov(cpc) +
  ##edgecov(cpp) +
  memory(type = "stability", lag = 1) +
  timecov(transform = function(t) t) +
  nodecov("genidx_multilevel") +
  nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") +
  cycle(3) + cycle(4) + cycle(5)




# ##--------------------------------------------
# ## Patch category similarity 
# ##--------------------------------------------
# patch.todo <- c('qualtrics',
#                 'clarabridge','confirmit','medallia','snap-surveys-ltd')
# for (name_i in patch.todo) {
#   cat(sprintf(' patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %n% 'cat_cos_sim' <- acf$.cov.categoryCosineSimilarity(nets[[t]])
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }


# ### Network size by scope of awareness
# name_i <- 'qualtrics'
# ds <- 1:10
# names(ds) <- as.character(ds)
# gstat <- sapply(ds, function(d){
#   g.base <- g.full  
#   g.d.sub <- igraph::make_ego_graph(g.base, d, V(g.base)[V(g.base)$name==name_i], 'all')[[1]]
#   return(c(v=vcount(g.d.sub),e=ecount(g.d.sub),d=igraph::graph.density(g.d.sub)))
# })
# 
# plot(gstat[1,], log='', type='b')
# plot(c(NA,diff(gstat[1,])), log='', type='b')





## get focal firm distances to specific firms over time
net.yrs <- 2007:2016
gs <- list()
for (yr in net.yrs) {
  gs[[as.character(yr)]] <- asIgraph(readRDS(file.path('firm_nets_rnr',sprintf('facebook_d3_y%d.rds',yr))))
}


g = gs[[length(gs)]]

co.idx <- which(
  grepl('content and publishing',V(g.full)$category_group_list) 
  & !(V(g.full)$employee_count %in% c('NA','-','1-10','11-50','51-100'))
)
co.names <- c(V(g.full)$name[co.idx],'amazon','aws')
ds = sapply(gs[2:10],function(g)igraph::distances(g, v = which(V(g)$vertex.names=='facebook'),
                                            to = V(g)[V(g)$vertex.names %in% co.names]))
# ds2 = t(igraph::distances(gs[[11]], 
#                         v = which(V(gs[[11]])$vertex.names=='facebook'),
#                         to = which(V(gs[[11]])$vertex.names %in% co.names)))
# ds = cbind(ds1,ds2)
# colnames(ds) <- names(gs)
rownames(ds) <- V(gs[[2]])$vertex.names[V(gs[[2]])$vertex.names %in% co.names]

row.names(ds) = V(g)$vertex.names[V(g)$vertex.names %in% co.names]
print(ds)

write.csv(ds, "facebook_news_dist_2008-2016_d3.csv")

igraph::distances(g.full, V(g.full)[V(g.full)$name=='facebook'],V(g.full)[V(g.full)$name=='hearstcorporation'])

## get cb category coverage
cgl = unique(c(str_split(V(g)$category_group_list,'[|]',simplify = T)))
cl = unique(c(str_split(V(g)$category_list,'[|]',simplify = T)))


# load('netrisk_dynamic_firm_nets.RData')

# # ADD CONSTRAINT NODE PROPERTY
# for (t in 1:length(nets)) {
#   g.tmp <- getIgraphFromNet(nets[[t]])
#   if (vcount(g.tmp)>0 & ecount(g.tmp)>0) {
#     cons <-  igraph::constraint(g.tmp)
#     cons[is.nan(cons) | is.na(cons)] <- 0 ### ???
#     nets[[t]] %v% 'constraint' <- cons
#   }
# }

#------------------------------------------------------
#              Predictors Diagnostics
# #------------------------------------------------------
# ## Plot density
# n <- ceiling(sqrt(length(firm.nets)))
# m <- ifelse(n*(n-1) >= length(firm.nets), n-1, n)
# par(mfrow=c(m,n), mar=c(2.5,2.5,2,1))
# for (firm_i in names(firm.nets)) {
#   nets <- firm.nets[[firm_i]]
#   plot(as.numeric(names(nets))-1,   
#        sapply(nets,function(net) {
#          sum(net[lower.tri(net)])/(nrow(net[,])*nrow(net[,]-1)/2) 
#        }), 
#        ylab='density', xlab='year', type='b', main=firm_i)
# }
# 
# ## Net Risk
# par(mfrow=c(3,3), mar=c(2.5,2.5,2,1))
# for (firm_i in names(firm.nets)) {
#   nets <- firm.nets[[firm_i]]
#   sapply(seq_along(nets), function(j) {
#     hist(nets[[j]] %v% 'net_risk', breaks=25, main=sprintf('%s %s',firm_i,names(nets)[j]))
#   })
# }

#------------------------------------------------------
