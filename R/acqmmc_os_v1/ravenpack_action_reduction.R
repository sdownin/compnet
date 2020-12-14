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
rvn_dir <- 'C:/Data/RavenPack'

setwd(work_dir)


## LOAD DATA AND DEPENDENCIES
# acf    <- source(file.path(version_dir,'acq_compnet_functions.R'))$value    ## acf: awareness functions
# cb     <- source(file.path(version_dir,'acq_cb_data_prep.R'))$value           ## cb : CrunchBase
# sdc    <- source(file.path(version_dir,'acq_sdc_coop.R'))$value               ## sdc: Thompson SDC
# si     <- source(file.path(version_dir,'acq_firm_size_controls.R'))$value     ## si : size controls from mergent intellect
# ih     <- source(file.path(version_dir,'acq_institutional_holdings.R'))$value ## ih : institutional holdings
# g.full <- source(file.path(version_dir,'acq_make_full_graph.R'))$value        ## g.full : full competition graph

# ## set firms to create networks (focal firm or replication study focal firms)


# t0 <- "2008-01-01 23:26:17.536" 
# t1 <- "2008-01-01 23:36:34.123"
# 
# as.POSIXct(t1, tz = 'UTC') - as.POSIXct(t0, tz = 'UTC')
# 
# strptime(t0, format = '%Y-%m-%d %H:%m:%s', tz = 'UTC')

##============================
##   ACTION CATEGORIES
##----------------------------
# dtafile <- 'rp40_action_categories_2000_2018.dta'
csvfile <- 'rp40_action_categories_2000_2018.csv'
##convert
rvn.all <- read.csv(file.path(rvn_dir, csvfile), stringsAsFactors = F)


## keep only 2006-2016
y1 <- 2006
y2 <- 2016
rvn <- rvn.all[which(rvn.all$year>=y1 & rvn.all$year<= y2), ]
## Change timestamp to POSIXct and add order ID
rvn <- rvn[order(rvn$timestamp_utc),]
rvn$timestamp_utc <- as.POSIXct(rvn$timestamp_utc, tz = 'UTC')
rvn$id <- 1:nrow(rvn)
# View(rbind(head(rvn),tail(rvn)))


## UTC time difference in *minutes*
## ens_elapsed in *milliseconds*
## so multimple by 60 for seconds, then by 1000 for milliseconds 
## >  as.numeric(rvn$timestamp_utc[4] - rvn$timestamp_utc[2]) * 60 * 1000
## >  1516464


## drop columns with all NA
col.drop <- c('evaluation_method','maturity')
rvn <- rvn[,which( ! names(rvn) %in% col.drop)]

# ## Reorder columns
# colnms <- names(rvn)
# colfront <- c('id','timestamp_utc','rp_entity_id','entity_type','entity_name',
#               'action_category','category','group','type','sub_type','property', 
#               'ens','ens_similarity_gap','ens_key','ens_elapsed',
#               'g_ens','g_ens_similarity_gap','g_ens_key','g_ens_elapsed',
#               'event_similarity_key',  'news_type','source','rp_story_id','rp_story_event_index','rp_story_event_count',
#               'product_key','company','isin', 'year','yr2','month')
# idx <- 1:ncol(rvn)
# idx.front <- which(colnms %in% colfront)
# idx.other <- which( ! colnms %in% colfront)
# rvn <- rvn[,c(idx.front, idx.other)]
# View(rbind(head(rvn,15),tail(rvn,15)))

         

### reduce RavenPack news items to unique actions (combine stories about same event) by "ens_key" (or "g_ens_key")
ra <- plyr::ddply(rvn, .(entity_name,ens_key,g_ens_key), summarize, .progress = 'text',
                  # Timestamp of EARLIEST report of event
                  timestamp_utc= min(timestamp_utc),
                  year= min(year),
                  ens_key= paste(unique(ens_key),collapse = '|'),
                  g_ens_key= paste(unique(g_ens_key),collapse = '|'),
                  # FIRM
                  entity_name= paste(unique(entity_name), collapse = '|'),
                  entity_type= paste(unique(entity_type), collapse = '|'),
                  rp_entity_id= paste(unique(rp_entity_id), collapse = '|'),
                  company= paste(unique(company), collapse = '|'),
                  product_key= paste(unique(product_key), collapse = '|'),
                  # EVENT
                  action_category= paste(unique(action_category), collapse = '|'),
                  category= paste(unique(category), collapse = '|'),
                  group= paste(unique(group), collapse = '|'),
                  type= paste(unique(type), collapse = '|'),
                  sub_type= paste(unique(sub_type), collapse = '|'),
                  property= paste(unique(property), collapse = '|'),
                  country_code= paste(unique(country_code), collapse = '|'),
                  relevance= paste(unique(relevance), collapse = '|'),
                  topic= paste(unique(topic), collapse = '|'),
                  # METADATA TO COMPARE STORIES OF SAME EVENT
                  ens_elapsed= paste(ens_elapsed,collapse = '|'),
                  g_ens_elapsed= paste(g_ens_elapsed,collapse = '|'),
                  ens= paste(ens,collapse = '|'),
                  ens_similarity_gap= paste(ens_similarity_gap,collapse = '|'),
                  g_ens= paste(g_ens,collapse = '|'),
                  g_ens_similarity_gap= paste(g_ens_similarity_gap,collapse = '|'),
                  event_similarity_key= paste(event_similarity_key,collapse = '|'),
                  news_type= paste(news_type,collapse = '|'),
                  source= paste(source,collapse = '|'),
                  # RP 
                  rp_story_id= paste(unique(rp_story_id),collapse = '|'),
                  rp_story_event_index= paste(unique(rp_story_event_index),collapse = '|'),
                  rp_story_event_count= paste(rp_story_event_count,collapse = '|'),
                  rp_position_id= paste(rp_position_id,collapse = '|'),
                  position_name= paste(position_name,collapse = '|'),
                  #  OTHER FIELDS (UNKNOWN)
                  isin= paste(unique(isin),collapse = '|'),
                  yr2= paste(unique(yr2), collapse = '|'),
                  ess= paste(unique(ess), collapse = '|'),
                  aes= paste(unique(aes), collapse = '|'),
                  aev= paste(unique(aev), collapse = '|'),
                  css= paste(unique(css), collapse = '|'),
                  nip= paste(unique(nip), collapse = '|'),
                  peq= paste(unique(peq), collapse = '|'),
                  bee= paste(unique(bee), collapse = '|'),
                  bmq= paste(unique(bmq), collapse = '|'),
                  bam= paste(unique(bam), collapse = '|'),
                  bca= paste(unique(bca), collapse = '|'),
                  ber= paste(unique(ber), collapse = '|'),
                  anl_chg= paste(unique(anl_chg), collapse = '|'),
                  mcq= paste(unique(mcq), collapse = '|'),
                  X_merge= paste(unique(X_merge), collapse = '|')
                  )
ra <- ra[order(ra$timestamp_utc),]
print(dim(ra))
View(ra)

## Save Ravenpack
rafile <- sprintf('ravenpack_resolved_actions_%s-%s.csv', y1,y2)
write.csv(ra, file = file.path(rvn_dir,rafile), row.names = F)


# ### Check for combined entities 
# ra[grep('[|]',ra$entity_name,ignore.case = T,perl = T), c('entity_name')]




# s1 <- "Zynga Inc.TENCEN|T HOLDINGS LTD."
idx.n2 <- grep('[|]',ra$entity_name)
View(ra[idx.n2,])


# ## all data 7GB
# dtafile <- 'rp40_2000_2016.dta'
# rvn.all <- read.dta13(file.path(rvn_dir, dtafile), select.rows=10)

## summarize action types
cnt <- plyr::ddply(ra, .(action_category, category), summarize, cnt=length(timestamp_utc))
cnt <- cnt[order(cnt$action_category),]
cnt$pct <- 100 * cnt$cnt / sum(cnt$cnt)
head(cnt, 20)


length(unique(rvn$entity_name))


































