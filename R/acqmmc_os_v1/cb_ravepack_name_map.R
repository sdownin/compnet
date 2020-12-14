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

## LOAD DATA AND DEPENDENCIES
acf    <- source(file.path(version_dir,'acq_compnet_functions.R'))$value    ## acf: awareness functions
cb     <- source(file.path(version_dir,'acq_cb_data_prep.R'))$value           ## cb : CrunchBase
sdc    <- source(file.path(version_dir,'acq_sdc_coop.R'))$value               ## sdc: Thompson SDC
si     <- source(file.path(version_dir,'acq_firm_size_controls.R'))$value     ## si : size controls from mergent intellect
ih     <- source(file.path(version_dir,'acq_institutional_holdings.R'))$value ## ih : institutional holdings
g.full <- source(file.path(version_dir,'acq_make_full_graph.R'))$value        ## g.full : full competition graph

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

## -- settings --
firms.todo <- c('ibm')
d <- 3
yrpd <- 1
startYr <- 2006
endYr <- 2017            ## dropping first for memory term; actual dates 2007-2016
## --------------  

name_i <- firms.todo[1]

#############################################################
## READ IN DATA FILE
firmfile <- sprintf('acq_sys_mmc_panel_regression_df_col%s_%s_%s-%s.csv',
                    50, name_i, startYr, endYr)
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

##factors
dfreg$i <- as.factor(dfreg$i)
dfreg$year <- as.factor(dfreg$year)
dfreg$y <- dfreg$y.cur + dfreg$y.new


##========================================
## RavenPack NAME MAPPING BEGIN 
##----------------------------------------
firms <- sort(unique(dfreg$name))
# RAVENPACK NAMES
rpn <- unique(rvn[,c('rp_entity_id','entity_name')])
rpn$name.l <- str_to_lower(rpn$entity_name)
# COHORT OF FIRMS
cbrpmap <- data.frame(
  name = unique(as.character(dfreg$name)),
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
  stringsAsFactors = F)

## string distance matrix for CrunchBase-to-RavenPack Name Mapping
for (i in 1:length(firms)) {
  firm_i <- gsub('\\W+', ' ', firms[i])
  if (i %% 10 == 0) cat(sprintf('i=%03s (%03.2f%s) %s \n',i,100*i/length(firms),'%',firms[i]))
  ##
  ## FIRST DO grep() string search, see if excatly match 1
  ##  - else if multiple matches, test distance on matches
  ##  - else test distance on all unmateched RavenPack firms
  ##
  idf <- data.frame(
    name=rpn$entity_name,
    name.l=rpn$name.l,
    lcs=stringdist(rpn$name.l, firm_i, method = 'lcs'),
    cosine=stringdist(rpn$name.l, firm_i, method = 'cosine'),
    jaccard=stringdist(rpn$name.l, firm_i, method = 'jaccard'),
    jw=stringdist(rpn$name.l, firm_i, method = 'jw'),
    qgram=stringdist(rpn$name.l, firm_i, method = 'qgram'),
    osa=stringdist(rpn$name.l, firm_i, method = 'osa')#,
    # lv=stringdist(rpn$name.l, firm_i, method = 'lv'),
    # dl=stringdist(rpn$name.l, firm_i, method = 'dl'),
  )
  idx <- apply(idf[,3:ncol(idf)], 2, which.min)
  idxcnt <- plyr::count(idx)
  idxcnt <- idxcnt[order(idxcnt$freq, decreasing = T),]
  if (nrow(idxcnt)==0 | length(idxcnt)==0) {
    cat(sprintf('missing strdist for i=%03s, `%s\n`',i,firms[i]))
    next
  }
  ## top match
  idx1 <- idxcnt$x[1]
  cbrpmap$idx1[i] <- idx1
  cbrpmap$idx1pct[i] <- idxcnt$freq[1] / sum(idxcnt$freq)
  cbrpmap$rp_entity_id_1[i] <- rpn$rp_entity_id[idx1]
  cbrpmap$entity_name_1[i] <- rpn$entity_name[idx1]
  ## secondary / tertiary matches
  if (nrow(idxcnt) >= 2) {
    idx2 <- idxcnt$x[2]
    cbrpmap$idx2[i] <- idx2
    cbrpmap$idx2pct[i] <- idxcnt$freq[2] / sum(idxcnt$freq)
    cbrpmap$rp_entity_id_2[i] <- rpn$rp_entity_id[idx2]
    cbrpmap$entity_name_2[i] <- rpn$entity_name[idx2]
  }
  if (nrow(idxcnt) >= 3) {
    idx3 <- idxcnt$x[3]
    cbrpmap$idx3[i] <- idx3
    cbrpmap$idx3pct[i] <- idxcnt$freq[3] / sum(idxcnt$freq)
    cbrpmap$rp_entity_id_3[i] <- rpn$rp_entity_id[idx3]
    cbrpmap$entity_name_3[i] <- rpn$entity_name[idx3]
  }
}

## save mapping dataframe
cbrpmapfile <- sprintf('crunchbase_ravenpack_name_mapping_n%s.csv', length(firms))
write.csv(cbrpmap, file=file.path(result_dir,cbrpmapfile))


