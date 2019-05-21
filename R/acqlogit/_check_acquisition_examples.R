##########################################################################################
#
# MMC Network & Acquisitions
#
# Create acquisition logit model covariates arrays
#
# @author   Stephen Downing <sdowning.bm02g@nctu.edu.tw>
#
# @export [list] cb
#
#
# ## update founded_on,closed_on dates  - Jin-Su's Email 2018-04-23
# ## C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\founded_on_date_edit
# co.date <- cb$readCsv('founded_on_date_edit/missing_companies_20180330.csv')
#
##########################################################################################

setwd("C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2")
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(parallel)
library(network, quietly = T)
library(texreg, quietly = T)
library(igraph, quietly = T)
library(plyr, quietly = T)
library(lattice, quietly = T)
library(latticeExtra, quietly = T)
library(ggplot2, quietly = T)
library(reshape2)

data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/"

script_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\R\\awareness_amj_rnr2'

## LOAD Scripts and Data
acf <- source(file.path(script_dir,'amj_awareness_functions.R'))$value ## FUNCTIONS 
cb  <- source(file.path(script_dir,'amj_cb_data_prep.R'))$value      ## DATA 

is.missing <- function(x)
{
  if(is.null(x)) 
    return(TRUE)
  return(is.na(x) | is.nan(x) | x == '')
}

##=======================================
##  PREP DATA and LOAD GRAPH
##---------------------------------------
## comptetition network
g.full <- read.graph('g_full.graphml', format='graphml')

## add comp net vertex IDs to acquisitions dataframe
co_acq <- cb$co_acq
gdf <- data.frame(acquirer_vid=as.integer(V(g.full)), 
                  acquirer_name_unique=V(g.full)$name,
                  acquirer_net_uuid=V(g.full)$company_uuid)
co_acq <- merge(co_acq, gdf, by='acquirer_name_unique', all.x = T, all.y = F)
gdf <- data.frame(acquiree_vid=as.integer(V(g.full)), 
                  acquiree_name_unique=V(g.full)$name,
                  acquiree_net_uuid=V(g.full)$company_uuid)
co_acq <- merge(co_acq, gdf, by='acquiree_name_unique', all.x=T, all.y = F)

## SORT CO_ACQ BY acquisition date
co_acq <- co_acq[order(co_acq$acquired_on, decreasing = F), ]


##==========================================
## CHECK EGO FIRM NETWORK PROPERTIES 
##------------------------------------------
# # 
# # ## filter acquisitions in recent period (if necessary)
# # co_acq_d <- co_acq[which(co_acq$acquired_on >= '1988-01-01'), ]
# # 
# # ## acquisition filtered by comp net firms union (either aquirer or acquired)
# # ## 14208
# # co_acq_g_u <- co_acq[which(co_acq$acquirer_name_unique %in% V(g.full)$name
# #                            | co_acq$acquiree_name_unique %in% V(g.full)$name), ]
# 
# ## acquisition filtered by comp net firms intersection (both aquirer and acquired)
# ## 3442
# co_acq_g_i <- co_acq[which(co_acq$acquirer_name_unique %in% V(g.full)$name
#                            & co_acq$acquiree_name_unique %in% V(g.full)$name), ]
# ## acquisition filtered by comp net firms intersection, recent date (>=2010)
# ## 2786
# co_acq_g_i_d <- co_acq_g_i[which(co_acq_g_i$acquired_on >= '1988-01-01'), ]
# 
# ## check top filtered acquirers
# cnt <- plyr::count(co_acq_g_i_d$acquirer_name_unique)
# cnt <- cnt[order(cnt$freq, decreasing = T),]
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
# 
# ## Check Acquisitions distribution by network distance included
# df.ego <- data.frame()
# name_i <- 'ibm'
# for (d in 1:6) {
#   g.ego <- igraph::make_ego_graph(graph = g.full, 
#                                   nodes = V(g.full)[V(g.full)$name==name_i], 
#                                   order = d, mode = 'all')[[1]]
#   mem <- igraph::multilevel.community(g.ego)$membership
#   mem.cnt <- plyr::count(mem)
#   mem.cnt <- mem.cnt[order(mem.cnt$freq, decreasing = T), ]
#   dim(mem.cnt)
#   ##  
#   acq.src <- co_acq_d[which(co_acq_d$acquirer_name_unique %in% V(g.ego)$name), ]
#   acq.src.trg <- co_acq_d[which(co_acq_d$acquirer_name_unique %in% V(g.ego)$name
#                                 & co_acq_d$acquiree_name_unique %in% V(g.ego)$name), ]
#   ##
#   df.tmp <- data.frame(d=d,v=vcount(g.ego),e=ecount(g.ego),
#                        acq.src=nrow(acq.src),acq.src.trg=nrow(acq.src.trg),
#                        r.in.out=round(nrow(acq.src.trg)/(nrow(acq.src)+nrow(acq.src.trg)),3),
#                        first.comp=min(E(g.ego)$relation_began_on))
#   df.ego <- rbind(df.ego, df.tmp)
# }; df.ego 


##--------------------------------------------------------------
##
##            CREATE FIRM NETWORK PERIOD LISTS  
##
##--------------------------------------------------------------

## SETTINGS
name_i <- 'google'
d <- 3
years <- 2006:2017

## get date periods and ego network based on settings
times <- sapply(years, function(x)paste0(x,'-01-01'))
start <- times[1]
end <- times[length(times)]

## EGO NETWORK
g.ego <- igraph::make_ego_graph(graph = g.full,
                                nodes = V(g.full)[V(g.full)$name==name_i],
                                order = d, mode = 'all')[[1]]

## NETWORKS IN TIMEFRAME TO PROCESS NODE COLLAPSE AND POCESS COVARIATES
g.pd      <- acf$makePdGraph(g.ego, start, end, isolates.remove=TRUE)   ## ego network
g.full.pd <- acf$makePdGraph(g.full, start, end, isolates.remove=TRUE)  ## full network

V(g.full.pd)$com <- igraph::multilevel.community(g.full.pd)$membership
cm <- plyr::count(V(g.full.pd)$com)
cm <- cm[order(cm$freq, decreasing = T), ]
cm

V(g.pd)$com <- igraph::multilevel.community(g.pd)$membership
# V(g.pd)$com <- igraph::multilevel.community(g.pd)$membership
cme <- plyr::count(V(g.pd)$com)
cme <- cme[order(cme$freq, decreasing = T), ]
cme

fml <- lapply(1:vcount(g.pd), function(v){
  nvids <- neighbors(g.pd, v)
  return(V(g.pd)$com[nvids])
})

fm.i <- c()
fm.j <- c()
for (i in 1:length(fml)) {
  fm.i <- c(fm.i, rep(i, length(fml[[i]])))
  fm.j <- c(fm.j, fml[[i]])
}
fmspmat <- sparseMatrix(i=fm.i, j=fm.j, x = 1,
                        dims = c(length(unique(fm.i)),length(unique(fm.j))))
fmmat <- as.matrix(fmspmat)
fmmat[fmmat > 1] <- 1 ## ??????

mc <- lapply(sort(unique(V(g.pd)$com)), function(x){
  vids <- which(V(g.pd)$com == x)
  return(unlist(strsplit(V(g.pd)$category_list[vids], '[|]')))
})

mcdf <- data.frame()
for (i in 1:length(mc)) {
  x <- mc[[i]]
  .df <- plyr::count(mc[[i]])
  .df <- .df[order(.df$freq, decreasing = T), ]
  .df$i <- i
  .df$name <- paste(.df$x[1:5], collapse = "|")
  mcdf <- rbind(mcdf, .df)
}
print(unique(mcdf$name))

# [1]  378        484         2554        3174        4909
# [1] "amazon"    "apple"     "facebook"  "google"    "microsoft"
firms <- c('amazon', 'apple', 'facebook', 'google', 'microsoft')
acq <- cb$co_acq[which(cb$co_acq$acquirer_name_unique %in% firms & cb$co_acq$acquired_year==2011), ]
.acqcats <- cb$co[,c('company_name_unique','category_list','category_group_list')]
names(.acqcats) <- c('company_name_unique','acquiree_category_list','acquiree_category_group_list')
acq <- merge(acq, .acqcats, by.x='acquiree_name_unique',by.y='company_name_unique',all.x=T,all.y = F)
acq$acquiree_commnity <- sapply(acq$acquiree_name_unique, function(x){
    if (x %in% V(g.pd)$name)  return(V(g.pd)$com[which(V(g.pd)$name %in% x)])
    return(NA)
  })
View(acq)
write.csv(acq, 'acq_5warships_2011_acquisitions.csv', row.names = F)

firmidx <- which(V(g.pd)$name %in% firms)

fmsub <- fmmat[firmidx, ]

gb <- graph_from_incidence_matrix(fmsub)

plot(gb, 
     vertex.label=firms,
     vertex.shape=sapply(V(gb)$type, function(x)ifelse(x, 'square','circle')),
     vertex.color=sapply(V(gb)$type, function(x)ifelse(x, 'gray','steelblue')),
     vertex.size=sapply(V(gb)$type, function(x)ifelse(x, 9, 15))
     )

g1 <- bipartite.projection(gb)$proj1
plot(g1, layout=layout.circle,
     edge.width=.2*E(g1)$weight,
     edge.label=E(g1)$weight,
     vertex.label=firms)



####












## CHECK NETWORK PERIOD SIZES
sapply(2:length(times), function(i){gi=acf$makePdGraph(g.ego, times[i-1], times[i], TRUE); return(c(e=ecount(gi),v=vcount(gi)))})

## PROCESS NODE COLLAPSE OF ACQUISITIONS BEFORE START OF TIMEFRAME
acqs.init <- co_acq[co_acq$acquired_on < start, ]
g.pd <- acf$nodeCollapseGraph(g.pd, acqs.init)  #remove.isolates ?
g.full.pd <- acf$nodeCollapseGraph(g.full.pd, acqs.init, verbose = TRUE)  ## remove.isolates

## SAVE INITIALIZED EGO NETWORK AND GLOBAL NETWORK
igraph::write.graph(g.pd, file=sprintf('g_%s_NCINIT_%s_%s',name_i,start,end), format = 'graphml')
igraph::write.graph(g.full.pd, file=sprintf('g_full_NCINIT_%s_%s',start,end), format = 'graphml')

## LOAD IN EGO NETWORK AND GLOBAL NETWORK
g.pd <- igraph::read.graph(sprintf('g_%s_NCINIT_%s_%s',name_i,start,end), format='graphml')
g.full.pd <- igraph::read.graph(sprintf('g_full_NCINIT_%s_%s',start,end), format='graphml')

## Full timeframe Clusters
V(g.pd)$nc <- as.integer(igraph::multilevel.community(g.pd)$membership)
V(g.full.pd)$nc <- as.integer(igraph::multilevel.community(g.full.pd)$membership)

## keep original timeframe graph
g.pd.orig <- g.pd
g.full.pd.orig <- g.full.pd



## YEAR PERIODS: DEFINE NICHE CLUSTERS
l <- list()
df.mmc <- data.frame()
df.rem <- data.frame()
df.reg <- data.frame()
lidx <- 0  ## acquisition list index
timeval <- timeval.last <- 0



## GET ALL ACQ EVENT VERTICES 
## keep only acquisitions with acquirer in ego network and target in global competition network
acq.src <- co_acq[ co_acq$acquirer_name_unique %in% V(g.pd)$name & co_acq$acquiree_name_unique %in% V(g.full.pd)$name, ]
acq.src.allpd <- acq.src[ acq.src$acquired_on >= start & acq.src$acquired_on < end , ]
acq.src.allpd <- acq.src.allpd[order(acq.src.allpd$acquired_on, decreasing = F), ]


##===============================
##
##  MAIN LOOP: COMPUTE COVARIATES
##
##-------------------------------

## ACQUISITION EVENTS:  UPDATE MMC & DYNAMIC EFFs
for (j in 1:nrow(acq.src.allpd)) {
  
  date_j <- acq.src.allpd$acquired_on[j]
  ## g.pd            d2 updated each acquisition
  ## g.pd.orig       d2 original
  ## g.full.pd.orig  global network within timeframe start, end
  
  cat(sprintf('\n\nstart %s end %s : acquisition %s\n\n',start,end,j))
  if ( !(acq.src.allpd$acquiree_name_unique[j] %in% V(g.full.pd.orig)$name) ) 
    next
  if ( !(acq.src.allpd$acquirer_name_unique[j] %in% V(g.pd)$name) ) 
    next
  ## skip if acquirer is not public
  isPublicAcq <- (acq.src.allpd$acquirer_name_unique[j] %in% cb$co_ipo$company_name_unique 
                  & cb$co_ipo$went_public_on[cb$co_ipo$company_name_unique==acq.src.allpd$acquirer_name_unique[j]] <= acq.src.allpd$acquired_on[j])
  if (length(isPublicAcq)==0)
    next
  if ( ! isPublicAcq)
    next

  lidx <- lidx + 1
  l[[lidx]] <- list()
  
  ## Update MMC after acquisition
  l[[lidx]]$mmc <- acf$getFmMmc(g.pd, as.integer(V(g.pd)$nc))
  
  ## SUM FM MMC over markets  ??????
  V(g.pd)$fm.mmc.sum <- rowSums(l[[lidx]]$mmc)
  V(g.pd)$num.mkts <- apply(l[[lidx]]$mmc, MARGIN=1, FUN=function(x){
    return(length(x[x>0]))
  })
  
  ## GET  DATAFRAME VARS
  
  ## Acquirer d2 original org.vid
  xi.orig.vid <- V(g.pd.orig)$orig.vid[which(acq.src.allpd$acquirer_name_unique[j] == V(g.pd.orig)$name)]
  ## target d2 original org.vid
  xj.orig.vid <- V(g.pd.orig)$orig.vid[which(acq.src.allpd$acquiree_name_unique[j] == V(g.pd.orig)$name)]
  xj.orig.vid <- ifelse(length(xj.orig.vid) > 1, xj.orig.vid, NA)
  ## acquirer  d2 t=j id
  xi <- which(V(g.pd)$name==acq.src.allpd$acquirer_name_unique[j])
  ## target  d2 t=j id
  xj <- which(V(g.pd)$name==acq.src.allpd$acquirer_name_unique[j])
  # acquirer id in original graph (at start of period)
  xi.orig <- as.integer(V(g.pd.orig)[V(g.pd.orig)$name==acq.src.allpd$acquirer_name_unique[j]])
  xi.nc <- as.integer(V(g.pd.orig)$nc[xi.orig]) ## original nc for the period
  #
  xi.mmc.sum <-  V(g.pd)$fm.mmc.sum[xi]
  xi.num.mkts <-  V(g.pd)$num.mkts[xi]
  ##
  xj <- as.integer(V(g.pd)[V(g.pd)$name==acq.src.allpd$acquiree_name_unique[j]])
  xj.orig <- ifelse( !is.na(xj.orig.vid), as.integer(V(g.pd.orig)[V(g.pd.orig)$orig.vid==xj.orig.vid]), NA)
  xj.orig <- ifelse(length(xj.orig) > 1, xj.orig, NA)
  xj.nc <- ifelse(length(xj)==0,NA,  V(g.pd.orig)$nc[xj.orig] )  ## original nc for the period
  
  ##--------------------------------------
  ## TARGET ALTERNATIVES SET
  ##--------------------------------------
  cat('target set ...')
  ## target alternative set vids
  targ.id <- which(V(g.full.pd)$name == acq.src.allpd$acquiree_name_unique[j])
  targ.vids.d2 <- igraph::neighborhood(graph = g.full.pd, order = 2, nodes = targ.id)[[1]]
  targ.vids.d2 <- targ.vids.d2[which( !(names(targ.vids.d2) %in% V(g.full.pd)$name[targ.id]))]
  targ.vids.d1 <- igraph::neighborhood(graph = g.full.pd, order = 1, nodes = targ.id)[[1]]
  targ.vids.d1 <- targ.vids.d1[which( !(names(targ.vids.d1) %in% V(g.full.pd)$name[targ.id]))]
  ## Target alternatives dataframe
  df.targ.alt <- cb$co[which(cb$co$company_name_unique %in% c(names(targ.vids.d1),names(targ.vids.d2),V(g.full.pd)$name[targ.id])), ]
  df.targ.alt$d <- sapply(df.targ.alt$company_name_unique, function(x){ return(
      ifelse(x ==  V(g.full.pd)$name[targ.id], 0, 
             ifelse(x %in% names(targ.vids.d1), 1,   2))
    )})
  ## ipo status
  df.targ.alt$is.public <- sapply(1:nrow(df.targ.alt), function(x){
    isNotOperating <- df.targ.alt$status[x] != 'operating'
    ipo.date <- cb$co_ipo$went_public_on[which(cb$co_ipo$company_name_unique == df.targ.alt$company_name_unique[x])]
    if (length(ipo.date)<1) 
      return(0)
    return(ifelse( isNotOperating & ipo.date <= date_j, 1, 0))
  })
  ## target had IPO
  df.targ <- df.targ.alt[which(df.targ.alt$company_name_unique == V(g.full.pd)$name[targ.id]), ]
  
  if (nrow(df.targ) == 0)
    next
  
  tmp <- df.targ.alt[df.targ.alt$company_name_unique %in% names(targ.vids.d2), ]
  ## select based on ownership status
  tmp <- if (df.targ$is.public == 1) {
      tmp[tmp$is.public == 1, ]
    } else {
      tmp[tmp$is.public == 0, ] 
    }
  # tmp.alt <- tmp[sample(1:nrow(tmp),size = min(9,nrow(tmp)),replace = F), ]
  tmp.alt <- tmp  ## keep all alternative targets
  ## combine target and alternatives for target set
  df.targ.alt <- rbind(tmp.alt, df.targ)
  ## add MMC
  df.targ.alt$fm.mmc.sum <- sapply(df.targ.alt$company_name_unique, function(name){
    ifelse(name %in% V(g.pd)$name, V(g.pd)$fm.mmc.sum[which(V(g.pd)$name == name)] , NA)
  })
  df.targ.alt$num.mkts <- sapply(df.targ.alt$company_name_unique, function(name){
    ifelse(name %in% V(g.pd)$name, V(g.pd)$num.mkts[which(V(g.pd)$name == name)] , NA)
  })
  cat('done.\n')
  ##--------------------------------------
  
  ##--------------------------------------
  ## ACQUIRER ALTERNATIVES SET
  ##--------------------------------------
  cat('acquirer set ...')
  ## acquirer alternative set vids
  acq.id <- which(V(g.pd)$name == acq.src.allpd$acquirer_name_unique[j])
  acq.vids.d2 <- igraph::neighborhood(graph = g.pd, order = 2, nodes = acq.id)[[1]]
  acq.vids.d2 <- acq.vids.d2[which( !(names(acq.vids.d2) %in% V(g.pd)$name[acq.id]))]
  acq.vids.d1 <- igraph::neighborhood(graph = g.pd, order = 1, nodes = acq.id)[[1]]
  acq.vids.d1 <- acq.vids.d2[which( !(names(acq.vids.d1) %in% V(g.pd)$name[acq.id]))]
  ## acquirer alternatives dataframe
  # length(acq.vids.d1)
  df.acq.alt <- cb$co[which(cb$co$company_name_unique %in% c(names(acq.vids.d1),names(acq.vids.d2),V(g.pd)$name[acq.id])), ]
  df.acq.alt$d <- sapply(df.acq.alt$company_name_unique, function(x){ return(
    ifelse(x ==  V(g.pd)$name[acq.id], 0, 
           ifelse(x %in% names(acq.vids.d1), 1,   2))
  )})
  ## ipo status
  df.acq.alt$is.public <- sapply(1:nrow(df.acq.alt), function(x){
    isNotOperating <- df.acq.alt$status[x] != 'operating'
    ipo.date <- cb$co_ipo$went_public_on[which(cb$co_ipo$company_name_unique == df.acq.alt$company_name_unique[x])]
    if (length(ipo.date)<1) 
      return(0)
    return(ifelse( isNotOperating & ipo.date <= date_j, 1, 0))
  })
  ## target had IPO
  df.acq <- df.acq.alt[which(df.acq.alt$company_name_unique == V(g.pd)$name[acq.id]), ]
  
  if (nrow(df.acq) == 0)
    next
  
  tmp <- df.acq.alt[df.acq.alt$company_name_unique %in% names(acq.vids.d2), ]
  ## select based on ownership status
  if (df.acq$is.public == 1) {
    tmp <- tmp[tmp$is.public == 1, ]
  } else {
    tmp <- tmp[tmp$is.public == 0, ] 
  }
  # tmp.acq.alt <- tmp[sample(1:nrow(tmp),size = min(9,nrow(tmp)),replace = F), ]
  tmp.acq.alt <- tmp ## keep all
  ## combine target and alternatives for target set
  df.acq.alt <- rbind(tmp.acq.alt, df.acq)
  ## set MMC
  df.acq.alt$fm.mmc.sum <- sapply(df.acq.alt$company_name_unique, function(name){
      ifelse(name %in% V(g.pd)$name, as.numeric(V(g.pd)$fm.mmc.sum[which(V(g.pd)$name == name)]) , NA)
    })
  df.acq.alt$num.mkts <- sapply(df.acq.alt$company_name_unique, function(name){
      ifelse(name %in% V(g.pd)$name, as.numeric(V(g.pd)$num.mkts[which(V(g.pd)$name == name)]) , NA)
    })
  cat('done.\n')
  # cat('acquirer set ...')
  # ## acquirer alternative set vids
  # acq.id <- which(V(g.full.pd)$name == acq.src.allpd$acquirer_name_unique[j])
  # acq.vids.d2 <- igraph::neighborhood(graph = g.full.pd, order = 2, nodes = acq.id)[[1]]
  # acq.vids.d2 <- acq.vids.d2[which( !(names(acq.vids.d2) %in% V(g.full.pd)$name[acq.id]))]
  # acq.vids.d1 <- igraph::neighborhood(graph = g.full.pd, order = 1, nodes = acq.id)[[1]]
  # acq.vids.d1 <- acq.vids.d2[which( !(names(acq.vids.d1) %in% V(g.full.pd)$name[acq.id]))]
  # ## acquirer alternatives dataframe
  # # length(acq.vids.d1)
  # df.acq.alt <- co[which(co$company_name_unique %in% c(names(acq.vids.d1),names(acq.vids.d2),V(g.full.pd)$name[acq.id])), ]
  # df.acq.alt$d <- sapply(df.acq.alt$company_name_unique, function(x){ return(
  #   ifelse(x ==  V(g.full.pd)$name[acq.id], 0, 
  #          ifelse(x %in% names(acq.vids.d1), 1,   2))
  # )})
  # ## ipo status
  # df.acq.alt$is.public <- sapply(1:nrow(df.acq.alt), function(x){
  #   isNotOperating <- df.acq.alt$status[x] != 'operating'
  #   ipo.date <- co_ipo$went_public_on[which(co_ipo$company_name_unique == df.acq.alt$company_name_unique[x])]
  #   if (length(ipo.date)<1) 
  #     return(0)
  #   return(ifelse( isNotOperating & ipo.date <= date_j, 1, 0))
  # })
  # ## target had IPO
  # df.acq <- df.acq.alt[which(df.acq.alt$company_name_unique == V(g.full.pd)$name[acq.id]), ]
  # tmp <- df.acq.alt[df.acq.alt$company_name_unique %in% names(acq.vids.d2), ]
  # if (df.acq$is.public == 1) {
  #   tmp <- tmp[tmp$is.public == 1, ]
  # } else {
  #   tmp <- tmp[tmp$is.public == 0, ] 
  # }
  # tmp.acq.alt <- tmp[sample(1:nrow(tmp),size = min(9,nrow(tmp)),replace = F), ]
  # ## combine target and alternatives for target set
  # df.acq.alt <- rbind(tmp.acq.alt, df.acq)
  # cat('done.\n')
  ##--------------------------------------

  ##--------------------------------------
  ##  NETWORK COVARIATES
  ##--------------------------------------
  cat('computing network covariates...')
  df.acq.alt$set <- 'acquirer'
  df.acq.alt$event <- sapply(df.acq.alt$d, function(d)ifelse(as.integer(d)==0, 1, 0))
  df.targ.alt$set <- 'target'
  df.targ.alt$event <- sapply(df.targ.alt$d, function(d)ifelse(as.integer(d)==0, 1, 0))
  df.alt <- rbind(df.acq.alt, df.targ.alt)
  df.alt$t <- j ## acquisition index
  df.alt <- df.alt[order(which(V(g.full.pd.orig)$name %in% df.alt$company_name_unique )), ] ## confirm ascencing order
  
  # ## Create Diff Graph (removed|acquired nodes are represented as isolates)
  vids <- which( V(g.full.pd)$name %in% df.alt$company_name_unique )
  vids.orig <- which( V(g.full.pd.orig)$name %in% df.alt$company_name_unique )
  vids.orig.rm <- vids[which( !(vids.orig %in% vids))]
  mapping <- V(g.full.pd.orig)[which(V(g.full.pd.orig)$orig.vid %in% V(g.full.pd)$orig.vid) ]
  g.diff <- igraph::contract.vertices(g.full.pd, mapping = mapping)
  V(g.diff)$name <- V(g.full.pd.orig)$name
  vids.diff <- as.integer( V(g.diff)[which( V(g.diff)$name %in% df.alt$company_name_unique )] )
  
  ## global covars
  tmp.cov <- data.frame(
    company_name_unique = unlist(V(g.diff)$name[vids.diff]),
    pow.n3 = unname(unlist(igraph::power_centrality(g.diff, nodes = vids.diff, exponent = -0.3))),
    pow.n2 = unname(unlist(igraph::power_centrality(g.diff, nodes = vids.diff, exponent = -0.2))),
    pow.n1 = unname(unlist(igraph::power_centrality(g.diff, nodes = vids.diff, exponent = -0.1))),
    deg = unname(igraph::degree(g.diff, v = vids.diff)),
    constraint = unname(unlist(igraph::constraint(g.diff, nodes = vids.diff)))
  )
  df.alt <- merge(df.alt, tmp.cov, by = 'company_name_unique', all.x = T, all.y = F)
  
  df.alt$acq.experience <- unlist(sapply(1:nrow(df.alt), function(x){ return(
    nrow(acq.src.allpd[which(acq.src.allpd$acquirer_name_unique == df.alt$company_name_unique[x]
                             & acq.src.allpd$acquired_on <= date_j), ]) / j ## scale experience to num observed acquisitions
  )}))
  ## local covars in pd graph
  df.alt$fm.mmc.sum <- sapply(df.alt$company_name_unique, function(name){
    ifelse(name %in% V(g.pd)$name, as.numeric(V(g.pd)$fm.mmc.sum[which(V(g.pd)$name == name)]), NA)
  })
  df.alt$num.mkts   <-  sapply(df.alt$company_name_unique, function(name){
    ifelse(name %in% V(g.pd)$name, as.numeric(V(g.pd)$num.mkts[which(V(g.pd)$name == name)]), NA)
  })
  
  cat('done.\n')
  ##---------------------------------------------
  
  l[[lidx]]$df.alt <- df.alt
  
  ##---------------------------------------------
  cat('appending dyadic regression dataframe...')
  
  for (k in 1:nrow(df.alt[df.alt$set=='acquirer', ])) {
    ix <- which( df.alt$company_name_unique == df.alt[df.alt$set=='acquirer', ][k, ]$company_name_unique )
    
    if (length(df.alt$event[ix])==0)
      next
    
    for (r in 1:nrow(df.alt[df.alt$set=='target', ])) {
      jx <- which( df.alt$company_name_unique == df.alt[df.alt$set=='target', ][r, ]$company_name_unique )
      
      ## skip pairing if neither acquirer nor target were in this actual event 
      if (length(df.alt$event[jx])==0)
        next
      if ( !as.integer(df.alt$event[ix]) & !as.integer(df.alt$event[jx]) )
        next
      
      if (df.alt$company_name_unique[ix] != df.alt$company_name_unique[jx]) {
        # cat(sprintf('ix %s jx %s\n',ix,jx))
        
        ij.dist <- igraph::distances(g.full.pd, 
                                     v = V(g.full.pd)[which(V(g.full.pd)$name == df.alt$company_name_unique[ix])], 
                                     to =V(g.full.pd)[which(V(g.full.pd)$name == df.alt$company_name_unique[jx])] )
        df.tmp.dyad <- data.frame(
          # event features 
          y = ifelse(as.integer(df.alt$event[ix]) & as.integer(df.alt$event[jx]), 1, 0),
          t = j,
          date = date_j,
          i = df.alt$company_name_unique[ix], 
          j = df.alt$company_name_unique[jx], 
          # acquirer covars 
          i.pow.n1 = df.alt$pow.n1[ix],
          i.pow.n2 = df.alt$pow.n2[ix],
          i.pow.n3 = df.alt$pow.n3[ix],
          i.deg = df.alt$deg[ix],
          i.fm.mmc.sum = ifelse(is.missing(df.alt$fm.mmc.sum[ix]), NA, df.alt$fm.mmc.sum[ix]),
          i.num.mkts = ifelse(is.missing(df.alt$num.mkts[ix]), NA, df.alt$num.mkts[ix]),
          i.constraint = df.alt$constraint[ix],
          i.acq.experience = df.alt$acq.experience[ix],
          # target covars 
          j.pow.n1 = df.alt$pow.n1[jx],
          j.pow.n2 = df.alt$pow.n2[jx],
          j.pow.n3 = df.alt$pow.n3[jx],
          j.deg = df.alt$deg[jx],
          j.fm.mmc.sum = ifelse(is.missing(df.alt$fm.mmc.sum[jx]), NA, df.alt$fm.mmc.sum[jx]),
          j.num.mkts = ifelse(is.missing(df.alt$num.mkts[jx]), NA, df.alt$num.mkts[jx]),
          j.constraint = df.alt$constraint[jx],
          j.acq.experience = df.alt$acq.experience[jx],
          # dyadic covars 
          ij.same.region = ifelse(df.alt$region[ix] == df.alt$region[jx], 1, 0),
          ij.same.state = ifelse(df.alt$state_code[ix] == df.alt$state_code[jx], 1, 0),
          ij.same.country = ifelse(df.alt$country_code[ix] == df.alt$country_code[jx], 1, 0),
          ij.same.employee.range = ifelse(df.alt$employee_count[ix] == df.alt$employee_count[jx], 1, 0),
          ij.dist = ifelse( class(ij.dist)=='matrix' & nrow(ij.dist)>0 & ncol(ij.dist)>0, ij.dist[1,1], Inf),
          ij.diff.pow.n1 = as.numeric(df.alt$pow.n1[ix]) - as.numeric(df.alt$pow.n1[jx]),
          ij.diff.pow.n2 = as.numeric(df.alt$pow.n2[ix]) - as.numeric(df.alt$pow.n2[jx]),
          ij.diff.pow.n3 = as.numeric(df.alt$pow.n3[ix]) - as.numeric(df.alt$pow.n3[jx]),
          ij.diff.deg = as.numeric(df.alt$deg[ix]) - as.numeric(df.alt$deg[jx]),
          ij.diff.fm.mmc.sum = ifelse(any(is.missing(df.alt$fm.mmc.sum[ix]),is.missing(df.alt$fm.mmc.sum[jx])), NA, as.numeric(df.alt$fm.mmc.sum[ix]) - as.numeric(df.alt$fm.mmc.sum[jx])),
          ij.diff.num.mkts = ifelse(any(is.missing(df.alt$num.mkts[ix]), is.missing(df.alt$num.mkts[jx])), NA, as.numeric(df.alt$num.mkts[ix]) - as.numeric(df.alt$num.mkts[jx])),
          ij.diff.constraint = as.numeric(df.alt$constraint[ix]) - as.numeric(df.alt$constraint[jx]),
          ij.diff.acq.experience = as.numeric(df.alt$acq.experience[ix]) - as.numeric(df.alt$acq.experience[jx])
        )
        df.reg <- rbind(df.reg, df.tmp.dyad)
      }
    }
  }
  
  cat('done.\n')
  ##---------------------------------------

  
  ##--------------------------------------------------------------------------    
  ## NODE COLLAPSE update network
  g.pd <- acf$nodeCollapseGraph(g.pd, acq.src.allpd[j,])
  g.full.pd <- acf$nodeCollapseGraph(g.full.pd, acq.src.allpd[j,])
  
  ## save incrementally
  if (lidx %% 10 == 0) {   
    saveRDS(list(l=l,df.reg=df.reg), file = sprintf("acqlogit_compnet_covs_list_%s.rds",name_i))
  }
  
  gc()
}

## final save
saveRDS(list(l=l,df.reg=df.reg), file = sprintf("acqlogit_compnet_covs_list_%s.rds",name_i))



tmp <- readRDS(sprintf("acqlogit_compnet_covs_list_%s.rds",name_i))
l <- tmp$l
df.reg <- tmp$df.reg

## save regression dataframe to seprate table file
write.csv(df.reg, file = sprintf("acqlogit_compnet_covs_df_%s.csv",name_i), row.names = F)




