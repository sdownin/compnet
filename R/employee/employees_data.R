#############################################################################################
#
#  Competition Networks and Employee Mobility
#
#  Data frame preparation
#
#############################################################################################
work_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2"
setwd(work_dir)

library(plyr)
library(reshape2)
library(lattice); library(latticeExtra)
library(ggplot2)
library(igraph)
library(stringr)
library(xergm)
library(btergm)
library(intergraph)
library(texreg)
library(parallel)

# source(file.path(work_dir,'R','employee','employee_cb_data_prep.R'))
eaf <- source(file.path(work_dir,'R','employee','employee_awareness_functions.R'))$value

sup_data_dir <- file.path(work_dir,'employee_sup_data')  ## supplmental data dir
# cb list filename
cb.list.file <- file.path(sup_data_dir,'cb_list.rds')
#
cb <- readRDS(cb.list.file)

data_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
par.default <- par()
lattice::trellis.par.set(strip.background=list(col="lightgrey"))


##
## Network stability measure (Jaccard index)
##
netJaccard <- function(m1, m2)
{
  ## num. maintained
  n11 <- sum(m1 * m2) ## element-wise multiplication
  ##  num. dropped or added
  n10.n01 <- sum( (m2 - m1) != 0 )
  ## jaccard index
  return( n11 / (n11 + n10.n01) )
}



##---------------------------------
## Compnet graph
##---------------------------------
## full compention network graph
gcfile <- file.path('employee','compnet_full.graphml')
gc.full <- igraph::read.graph(gcfile, format='graphml')
#---------------------------------------------------------------------------------

# cb$csv$jobs      <- 'jobs.csv'
# cb$csv$ppl       <- 'people.csv'
# cb$csv$ppl_desc  <- 'people_descriptions.csv'
# job.all       <- cb$readCsv(file.path(data_dir, cb$csv$jobs))
# ppl       <- cb$readCsv(file.path(data_dir, cb$csv$ppl))
# ppl_desc  <- cb$readCsv(file.path(data_dir, cb$csv$ppl_desc))

job.all  <- cb$job
ppl      <- cb$ppl
ppl_desc <- cb$ppl_desc


# job.all$title1 <- NA
# job.all$title2 <- NA
# job.all[,c('title1','title2')] <- unname(t(sapply(job.all$title, function(x){
#   parts <- str_trim(stringr::str_split(x, '[&/|]|(and)',simplify=T))
#   parts <- parts[which(parts != '')]
#   j1 <- (parts[1])
#   j2 <- ifelse(length(parts) > 1, parts[2], NA)
#   return(c(title1=j1, title2=j2))
# })))
# 
# job.all$started_on_year <- year(ymd(job.all$started_on))
# 
# ##update TRUE/FALSE bool values in job.all boolean fields
# parseBool <- function(x) {
#   x[which(x=='t')] <- TRUE
#   x[which(x=='f')] <- FALSE
#   x[which(x=='')] <- NA
#   x <- as.logical(x)
#   return(x)
# }
# job.all$is_current <- parseBool(job.all$is_current)
# job.all$job_role <- parseBool(job.all$job_role)
# job.all$executive_role <- parseBool(job.all$executive_role)
# job.all$advisory_role <- parseBool(job.all$advisory_role)
# 
# ##-----------------------------------------------------------
# ## Add company founded_on date for job started_on date
# ## if job attributes:  job_role=t, executive_role=t
# ##-----------------------------------------------------------
# job.all <- merge(job.all, cb$co[,c('company_uuid','employee_count','founded_on','closed_on')], 
#                  by.x='org_uuid',by.y='company_uuid',all.x=T,all.y=F)
# job.all$founded_on <- ymd(job.all$founded_on)
# job.all$closed_on <- ymd(job.all$closed_on)
# job.all$started_on <- ymd(job.all$started_on)
# job.all$ended_on <- ymd(job.all$ended_on)
# # jobs that are executive/job role for company not known to be >1000 employees
# .idx <- which(
#   is.na(job.all$started_on) 
#   & job.all$job_role==TRUE 
#   & job.all$executive_role==FALSE
#   & job.all$employee_count %in% c('1-10','11-50','51-100','51-200','101-250','201-500','251-500','501-1000',NA,'-')
# )
# job.all$started_on[.idx] <- job.all$founded_on[.idx]
# 
# ##-----------------------------------------------------------
# ## add tenure to jobs with end date
# ##-----------------------------------------------------------
# .idx <- which(!is.na(job.all$started_on) & !is.na(job.all$ended_on))
# job.all$tenure <- NA
# job.all$tenure[.idx] <- as.numeric(job.all$ended_on[.idx] - job.all$started_on[.idx]) / 365.2422 ## days per year
# # add truncated tenure to current jobs based on export date
# trunc_date <- ymd('2016-10-24')
# .idx <- which(job.all$is_current & is.na(job.all$ended_on))
# job.all$tenure[.idx] <- as.numeric(trunc_date - job.all$started_on[.idx]) / 365.2422 ## days per year
# job.all$tenure_trunc <- FALSE
# job.all$tenure_trunc[.idx] <- TRUE  ## current job with no ended_on date (as of trunc 2016-10-24)
# ##remove erroneous job data with negative tenure or tenure over 100 yrs:  876 rows (0.000853427 of cases)
# .idx <- which(
#   (job.all$tenure >= 0 & job.all$tenure < 100) | 
#     is.na(job.all$tenure)
# )
# ## CREATE job dataframe
# job <- job.all[.idx, ]
# dim(job)
# 
# ## only keep jobs at organizations that are companies (not universities, nonprofits, etc.)
# job <- job[which(job$org_uuid %in% cb$co$company_uuid),]
# # ## keep only jobs with a tenure
# job <- job[!is.na(job$tenure),]
# dim(job)
# 
# ## check distribution
# ggplot(aes(x=job$tenure, fill=tenure_trunc, colour=tenure_trunc), data=job) + 
#   geom_histogram(alpha=.5, bins=55) + scale_x_log10() + theme_bw() +
#   xlab('Years')
# 
# ##---------------------------------------------
# ## filter only from-->to pairs of jobs
# ##---------------------------------------------
# ## filter jobs by persons with >1 job
# percnt <- plyr::count(job$person_uuid)
# percnt <- percnt[which(percnt$freq > 1), ]
# dim(percnt)
# persub <- percnt$x
# job <- job[which(job$person_uuid %in% persub),]
# ## order
# job <- job[order(job$started_on, decreasing = F),]

## loop to create job edgelist -----------
# eljob <- data.frame()
# for (i in 1:length(persub)) {  #length(persub)
#   peruuid <- persub[i]
#   if (i %% 100 == 0) cat(sprintf('\ni=%s %s',i,peruuid))
#   ji <- job[which(job$person_uuid == peruuid), ]
#   if (nrow(ji) < 2)
#     next
#   for (k2 in 2:nrow(ji)) {
#     k1 <- k2-1
#     job2 <- ji[k2,]   ## to df
#     names(job2) <- sapply(names(job2),function(x)sprintf('to_%s',x), USE.NAMES = F)
#     job1 <- ji[k1,]   ## from df
#     names(job1) <- sapply(names(job1),function(x)sprintf('from_%s',x), USE.NAMES = F)
#     jobmv <- cbind(job1,job2)
#     jobmv$days_between <- as.numeric(job2$to_started_on - job1$from_ended_on)
#     ## append to job edgelist
#     eljob <- rbind(jobmv, eljob)
#   }
# }
# eljob$person_uuid <- eljob$from_person_uuid
# eljob$to_person_uuid <- NULL
# eljob$from_person_uuid <- NULL
# eljobnms <- names(eljob)
# .id1 <- which(eljobnms == 'from_org_uuid')
# .id2 <- which(eljobnms == 'to_org_uuid')
# .id3 <- which(eljobnms == 'person_uuid')
# .id4 <- which(eljobnms == 'from_ended_on')
# .id5 <- which(eljobnms == 'to_started_on')
# .id6 <- which(eljobnms == 'days_between')
# .idz <- 1:ncol(eljob)
# eljob <- eljob[ , c(.id1,.id2,.id3,.id4,.id5,.id6,.idz[-c(.id1,.id2,.id3,.id4,.id5,.id6)])]
## SAVE
# saveRDS(eljob, file='edgelist_job.rds')
# write.csv(eljob, file = 'edgelist_job.csv', row.names = F)






## READ SAVED edge list data frame -------------
eljob <- readRDS('edgelist_job.rds')
dim(eljob)
## Remove NAs
eljob <- eljob[which(eljob$to_org_uuid != 'NA' & eljob$from_org_uuid != 'NA'), ]

# eljob <- read.csv(edgelist_job.csv, header = T)

# eljob$id <- apply(eljob[,1:2],1,function(x)paste(x,collapse = '|'))

## Employee attributes (gender ratios to be averaged when simplifying multiple edges) ------------------
pplgender <- ppl[,c('uuid','gender')]
pplgender$male_avg <- 0
pplgender$female_avg <- 0
pplgender$male_avg[which(pplgender$gender=='Male')] <- 1
pplgender$female_avg[which(pplgender$gender=='Female')] <- 1
pplgender$gender <- NULL
eljob <- merge(eljob, pplgender, by.x='person_uuid',by.y='uuid',all.x=T,all.y=F)
# eljob <- eljob[,c(which(names(eljob)=='from_org_uuid'),which(names(eljob)=='to_org_uuid'),
#                   which( ! names(eljob) %in% c('from_org_uuid','to_org_uuid')))]

## Dyadic job role measure
eljob$to_from_exec_role_diff_avg <- eljob$to_executive_role - eljob$from_executive_role
eljob$from_exec_role_avg <- as.numeric( eljob$from_executive_role )
eljob$to_exec_role_avg <- as.numeric( eljob$to_executive_role )

##==============================================
## CREATE Organization vertex dataframe
##----------------------------------------------


##==============================================
## GRAPH of employees
##----------------------------------------------
## Reset all dates to characters -- igraph doesn't handle dates format
strdatecols <- c('to_founded_on','from_founded_on','from_started_on','to_started_on',
                 'from_ended_on','to_ended_on')
for (col in strdatecols) {
  eljob[,col] <- as.character(eljob[,col])
}

## from_company_name_unique and to_company_name_unique
.tcnu <- data.frame(to_company_name_unique=cb$co$company_name_unique, to_org_uuid=cb$co$company_uuid, stringsAsFactors = F)
.fcnu <- data.frame(from_company_name_unique=cb$co$company_name_unique, from_org_uuid=cb$co$company_uuid, stringsAsFactors = F)
#
eljob <- merge(x=eljob, y=.tcnu, by='to_org_uuid', all.x=T, all.y=F)
eljob <- merge(x=eljob, y=.fcnu, by='from_org_uuid', all.x=T, all.y=F)
## set to->from company_name_unique as first two columns
eljobnms <- names(eljob)
colidx <- c(which(eljobnms=='from_company_name_unique'), which(eljobnms=='to_company_name_unique'), 
            which( ! eljobnms %in% c('from_company_name_unique', 'to_company_name_unique')))
eljob <- eljob[,colidx]

## Filter only edges where to_started_on_year >= 2007
eljob <- eljob[which(eljob$to_started_on_year >= 2007), ]

## filter only edges with vertices in full compnet
eljob <- eljob[which(eljob$from_org_uuid %in% V(gc.full)$company_uuid & eljob$to_org_uuid %in% V(gc.full)$company_uuid),]

## which uuids in el.from el.to and compnet
ev1 <- unique(eljob$from_org_uuid)
ev2 <- unique(eljob$to_org_uuid)
vv <- V(gc.full)$company_uuid
vall <- intersect(intersect(ev1, ev2), vv)
##
eljob <- eljob[which(eljob$from_org_uuid %in% vall & eljob$to_org_uuid %in% vall),]
##  still too many vids in 'vall' need to double check remove extras here
euuid <- unique(c(eljob$from_org_uuid,eljob$to_org_uuid))

##------------------------------
## CREATE vertex dataframe for employee graph 
##------------------------------
vdf <- as_data_frame(gc.full, 'vertices')
vdf <- vdf[which(vdf$company_uuid %in% vall & vdf$company_uuid %in% euuid),]
# vdf <- vdf[,c(which(names(vdf)=='company_uuid'), which(names(vdf)!='company_uuid'))]
names(vdf)[c(which(names(vdf)=='name'))] <- c('company_name_unique')
# 
# # # # ## add vertices from edgelist (missing in vdf) into vdf
# ej <- unique(c(eljob$from_org_uuid,eljob$to_org_uuid))
# vj <- unique(vdf$company_uuid)
# ej[which( ! ej %in% vj)]
# vj[which( ! vj %in% ej) ]

##
ge <- graph.data.frame(eljob, directed = T, vertices = vdf)
## remove loops
E(ge)$weight <- 1
ge <- igraph::simplify(ge, remove.multiple = FALSE, remove.loops = T)
##-----------------------------------------------------
## remove edges where delay between leaving job and staring job is < 0 days or > 180 days
##  -- filters out 6096/16355   (37.27%) ?????? ***********
#------------------------------------------------------
deleids <- which(E(ge)$days_between < 0 | E(ge)$days_between > 180 )
ge <- igraph::delete.edges(ge, deleids)
# ## remove isolates (after removing loops and taking compnet subgraph)
# ge <- induced_subgraph(ge, which(igraph::degree(ge)>0))

# ## DECOMPOSE INTO COMPONENTS
# ge.comp <- decompose(ge, mode = 'weak', min.vertices = 3)
# ge <- ge.comp[[which.max(sapply(ge.comp,vcount))]]

## create compnet sugraph for nodes in empnet
gc <- induced.subgraph(gc.full, which(V(gc.full)$company_uuid %in% V(ge)$company_uuid))
# gc <- induced.subgraph(gc, which(V(gc)$company_uuid != 'NA'))

## Check compnet and empnet have same node order
all(V(ge)$name==V(gc)$name)
all(V(ge)$company_name_unique==V(gc)$company_name_unique)

## edge data frame
# ge.el <- igraph::as_data_frame(ge, 'edges')
# dim(ge.el[ge.el$days_between>=0,])





# ## COMMUNITY FOR COLORATION
# comm <- igraph::multilevel.community(as.undirected(ge))
# ncomm <- length(unique(comm$membership))
# # comm <- igraph::edge.betweenness.community(ge)
# V(ge)$cluster <- comm$membership
# 
# ## PLOT
# # png('employee_network_v1.png', height = 10, width = 10, units = 'in', res = 200)
# #   par(mar=c(.1,.1,.1,.1))
# #   plot.igraph(ge, layout=layout.fruchterman.reingold(ge, grid="nogrid"),
# #               vertex.color=rainbow(ncomm,alpha = .6)[V(ge)$cluster],
# #               vertex.size=2,
# #               vertex.label=NA,
# #               vertex.label.cex = 0,
# #               edge.width=.5,
# #               edge.arrow.size=.3)
# # dev.off()

### -------------------------------
##  START WITH ROOT FIRM EGO NETWORKS 
##   -- employee net --> filter compnet
##---------------------------------
## d-th order ego network from root firm
# rootid <- cb$co$company_uuid[which(cb$co$company_name_unique=='microsoft')]
gesub <- igraph::make_ego_graph(ge, order = 2, nodes = which(V(ge)$name=='microsoft'), mode = 'all')[[1]]
gcsub <- induced.subgraph(gc, which(V(gc)$name %in% V(gesub)$name))

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
gfuuid <- V(gcsub)$company_uuid
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
gcsub <- eaf$nodeCollapseGraph(gcsub, par.chi.nc, remove.isolates=T, verbose = T)
cat('done.\n')

##_-----------------------------------------------------
##-------process pre-start-year acquisitions------------
##------------------------------------------------------
## PERIODS
periods <- list(2008:2009,2010:2011,2012:2013,2014:2015)
# periods <- list(2010,2011,2012,2013,2014,2015)
startYr <- min(unlist(periods))
endYr <- max(unlist(periods))
## preprocess acquisitions before first period
acqs.pd <- cb$co_acq[cb$co_acq$acquired_on <= sprintf('%d-12-31',startYr-1), ]
gcsub <- eaf$nodeCollapseGraph(gcsub, acqs.pd, remove.isolates=T, verbose = T)
netsub <- asNetwork(gcsub)
### remove vertices not yet founded or already closed
vids <- which(V(gcsub)$founded_year <= endYr & !(V(gcsub)$closed_on < startYr))
gcsub <- induced.subgraph(gcsub, vids)
## and filter employee network to only vertices from firms in compnet (those that were founded in time and not )
gesub <- induced.subgraph(gesub, which(V(gesub)$company_uuid %in% V(gcsub)$company_uuid))
## summarize network sizes
cat(sprintf('empnet sub:  v = %d, e = %d\ncompnet sub: v = %d, e = %d\n',
            vcount(gesub),ecount(gesub),vcount(gcsub),ecount(gcsub)))

# llop  EMPLOYEE Period Network
geli <- list()
gcli <- list()
for (years in periods) {
  
  cat(sprintf('\n\n\n-------------------%s----------------------\n\n',paste(years,collapse = '-')))
  
  pdstr <- paste(years,collapse='|')
  
  ##-----------------------------------------------------
  ## EMPLOYEE NETWORK------------------------------------
  ##-----------------------------------------------------
  # print(years)
  eidsdel <- which( ! E(gesub)$to_started_on_year %in% years )
  gesubt <- igraph::delete_edges(gesub, edges = eidsdel)
  #
  edge.attr.comb = list(weight='sum',
                        from_exec_role_avg=function(x)mean(x,na.rm=T),
                        to_exec_role_avg=function(x)mean(x,na.rm=T),
                        to_from_exec_role_diff_avg=function(x)mean(x,na.rm=T),
                        male_avg=function(x)mean(x,na.rm=T),
                        female_avg=function(x)mean(x,na.rm=T),
                        from_tenure=function(x)mean(x,na.rm=T), 
                        days_between=function(x)mean(x,na.rm=T), 
                        'ignore')
  gesubt <- simplify(gesubt, remove.multiple = T, edge.attr.comb = edge.attr.comb)
  
  ## save network attrs from edge attrs (for TERGM edgecov predictors)
  enet <- asNetwork(gesubt)
  # #
  enet %n% 'from_exec_role_avg' <- as_adjacency_matrix(gesubt, attr='from_exec_role_avg', sparse = F)
  enet %n% 'to_exec_role_avg' <- as_adjacency_matrix(gesubt, attr='to_exec_role_avg', sparse = F)
  enet %n% 'to_from_exec_role_diff_avg' <- as_adjacency_matrix(gesubt, attr='to_from_exec_role_diff_avg', sparse = F)
  # enet %n% 'male_avg' <- as_adjacency_matrix(gesubt, attr='male_avg', sparse = F)
  # enet %n% 'female_avg' <- as_adjacency_matrix(gesubt, attr='female_avg', sparse = F)
  enet %n% 'from_tenure' <- as_adjacency_matrix(gesubt, attr='from_tenure', sparse = F)
  enet %n% 'days_between' <- as_adjacency_matrix(gesubt, attr='days_between', sparse = F)
  # enet %n% 'weight' <- as_adjacency_matrix(gesubt, attr='weight', sparse = F)
  ### AGE
  enet %v% 'age' <- (max(years)+1) - as.numeric(V(gesubt)$founded_year)
  enet %n% 'age_diff' <- outer( enet %v% 'age', enet %v% 'age', '-' )  ## this is backward [from]-[to] ## later need to multiple by -1 to get forward [to]-[from]


  cl <- edge.betweenness.community(gesubt, weights = NULL)
  enet %v% 'cluster_between' <- cl$membership
  enet %v% 'between' <- igraph::betweenness(gesubt, weights = NULL)
  enet %v% 'degree_out' <- igraph::degree(gesubt, mode = 'out')
  enet %v% 'degree_in' <- igraph::degree(gesubt, mode = 'in')
  enet %v% 'degree_total' <- igraph::degree(gesubt, mode = 'total')
  enet %v% 'eigen_centrality' <- eigen_centrality(gesubt, directed = T, weights = NULL, scale = T)$vector

  n <- length(enet$val)
  x <- ergm::ergmMPLE(formula = enet ~ triangle + transitive + twopath +
                        dgwesp(0, fixed=T, type="OTP") + dgwesp(0, fixed=T, type="OSP") ,
                      output = 'array', control=control.ergm(MPLE.max.dyad.types=n^2) )$predictor
  x1 <- x[,,1]
  x1[is.na(x1)] <- 0
  x2 <- x[,,2]
  x2[is.na(x2)] <- 0
  x3 <- x[,,3]
  x3[is.na(x3)] <- 0
  x4 <- x[,,4]
  x4[is.na(x4)] <- 0
  x5 <- x[,,5]
  x5[is.na(x5)] <- 0

  enet %n% 'triangle'   <- x1
  enet %n% 'transitive' <- x2
  enet %n% 'twopath'    <- x3
  enet %n% 'dgwesp.otp' <- x4
  enet %n% 'dgwesp.osp' <- x5

  ##-----------------------------------------------------
  ## COMPETITION NETWORK --------------------------------
  ##-----------------------------------------------------
  ## compnet functions use end year EXCLUSIVE; so end=2017 actually stops at 2016
  ##  --> adjust employee network end year +1 = compnet end year
  start <- min(years)
  end <- max(years) + 1
  ##
  gcsubt <- induced.subgraph(gcsub, which(V(gcsub)$company_uuid %in% V(gesubt)$company_uuid))
  #
  pd <- paste(years,collapse='|')
  t1 <- sprintf('%s-01-01', start)
  t2 <- sprintf('%s-12-31', end-1)
  ## 1. Node Collapse acquisitions within period
  acqs.pd <- cb$co_acq[cb$co_acq$acquired_on >= t1 & cb$co_acq$acquired_on <= t2, ]
  gcsubt <- eaf$nodeCollapseGraph(gcsubt, acqs.pd, verbose = T)
  ## 2. Subset Period Network
  cnet <- eaf$makePdNetwork(asNetwork(gcsubt), start, end, isolates.remove = F)
  ## 3. Set Covariates for updated Period Network
  cnet <- eaf$setCovariates(cnet, start, end,
                            covlist=c(
                                    'mmc',
                                      'dist','invdist','ipo_status',
                                      # 'constraint','similarity','centrality',
                                      # 'generalist','coop','employee','sales',
                                      'cat_cos_sim'
                                      # 'shared_competitor','shared_investor'
                            ),
                            acq=cb$co_acq, br=cb$co_br, ipo=cb$co_ipo,
                            rou=cb$co_rou, inv_rou=cb$inv_rou, inv=cb$inv,
                            coop=sdc, ih=ih, size=si)

  cl <- edge.betweenness.community(gcsubt, weights = NULL)$membership
  cnet %v% 'cluster_between' <- cl
  cnet %n% 'same_cluster_between' <- as.matrix(outer(cl, cl, '==') * 1)
  cnet %v% 'eigen_centrality' <- eigen_centrality(gcsubt, directed = F, weights = NULL, scale = T)$vector

  cl <- multilevel.community(gcsubt)$membership
  cnet %v% 'cluster_multilevel' <- cl
  cnet %n% 'same_cluster_multilevel' <- as.matrix(outer(cl, cl, '==') * 1)

  ## add compnet vertex attrs to empnet
  enet %v% 'ipo_status' <- cnet %v% 'ipo_status'

  ## UPDATE
  gcli[[pdstr]] <- cnet
  geli[[pdstr]] <- enet
  
} 
summary(geli)

## SAVE ------------------
empnetsfile <- sprintf('employee_competition_networks_lists_pd%s_yr%s_%s-%s.rds',
                       length(periods),endYr-startYr+1,startYr,endYr)
# saveRDS(list(geli=geli, gcli=gcli), file = empnetsfile)
.nets <- readRDS(empnetsfile)
geli <- .nets$geli
gcli <- .nets$gcli

# for (i2 in 2:length(geli)) {
#   i1 <- i2-1    # class(net[,])
#   cat(sprintf('\n\n%s-%s\n',names(geli)[i1],names(geli)[i2]))
#   print( netJaccard(geli[[i1]][,], geli[[i2]][,]) )
# }

for (i in 1:length(geli)) {
  network::delete.network.attribute(geli[[i]], attrname = 'dgwesp.osp')
}

##-------------PATCH ------------------
## add eigenvector centrality
## load saved data lists 
##-------------------------------------
.nets <- readRDS(empnetsfile)
geli <- .nets$geli
gcli <- .nets$gcli
#
for (i in 1:length(geli)) {
  cat(sprintf('\n\n----------%s: %s------------\n',i,names(geli)[i]))
  enet <- geli[[i]]
  cnet <- gcli[[i]]
  # ge <- asIgraph(enet)
  # gc <- asIgraph(cnet)
  # # E(ge)$weight <- NULL
  # # E(gc)$weight <- NULL
  # ## eigenvector centrality
  # enet %v% 'eigen_centrality' <- eigen_centrality(ge, directed = T, weights = NULL, scale = T)$vector
  # cnet %v% 'eigen_centrality' <- eigen_centrality(gc, directed = F, weights = NULL, scale = T)$vector
  # #
  # # compnet cluster
  # cl <- multilevel.community(gc)$membership
  # cnet %v% 'cluster_multilevel' <- cl
  # cnet %n% 'same_cluster_multilevel' <- as.matrix(outer(cl, cl, '==') * 1)
  
  clust <- enet %v% 'cluster_between'
  xlocal <- outer(clust, clust, '==') * 1
  
  ## Add change statistic predictors
  n <- length(enet$val)
  x <- ergm::ergmMPLE(formula = enet ~ triangle + ttriple + cycle(3) + cycle(4) + 
                        nearsimmelian + localtriangle(xlocal) + 
                        dgwesp(0, fixed=T, type="OTP") + 
                        dgwesp(0, fixed=T, type="OSP") + 
                        dgwnsp(0, fixed=T, type="OTP") + 
                        dgwnsp(0, fixed=T, type="OSP") 
                      ,output = 'array',
                      control=control.ergm(MPLE.max.dyad.types=n^2) )$predictor
  x1 <- x[,,1]
  x1[is.na(x1)] <- 0
  x2 <- x[,,2]
  x2[is.na(x2)] <- 0
  x3 <- x[,,3]
  x3[is.na(x3)] <- 0
  x4 <- x[,,4]
  x4[is.na(x4)] <- 0
  x5 <- x[,,5]
  x5[is.na(x5)] <- 0
  x6 <- x[,,6]
  x6[is.na(x6)] <- 0
  x7 <- x[,,7]
  x7[is.na(x7)] <- 0
  x8 <- x[,,8]
  x8[is.na(x8)] <- 0
  x9 <- x[,,9]
  x9[is.na(x9)] <- 0
  x10 <- x[,,10]
  x10[is.na(x10)] <- 0
  
  enet %n% 'triangle'     <- x1
  enet %n% 'ttriple'   <- x2
  enet %n% 'cycle3' <- x3
  enet %n% 'cycle4' <- x4
  enet %n% 'nearsimmelian' <- x5
  enet %n% 'localtriangle' <- x6
  # enet %n% 'dgwesp.otp' <- x7
  # enet %n% 'dgwesp.osp' <- x8
  # enet %n% 'dgwnsp.otp' <- x9
  # enet %n% 'dgwnsp.osp' <- x10

  ## ##UPDATE list
  geli[[i]] <- enet
  gcli[[i]] <- cnet
}

# enet <- geli[[2]]
# pmu <- pmuncommon[[2]]
# status.t <- a.pmu.status[[2]]
# cossim.t <- a.pmu.cossim[[2]]
# cdf1 <- data.frame(
#   y=c(enet[,]),
#   status = c(status.t),
#   cossim = c( cossim.t ),
#   triangle= c( pmu * enet %n% 'triangle'),
#   ttriple=  c( pmu * enet %n% 'ttriple'),
#   cycle3= c( pmu * enet %n% 'cycle3'),
#   cycle4= c( pmu * enet %n% 'cycle4'),
#   nearsimmelian= c( pmu * enet %n% 'nearsimmelian'),
#   localtriangle= c( pmu * enet %n% 'localtriangle'),
#   `dgwesp.otp`= c( pmu * enet %n% 'dgwesp.otp'),
#   `dgwesp.osp`= c( pmu * enet %n% 'dgwesp.osp'),
#   `dgwnsp.otp`= c( pmu * enet %n% 'dgwnsp.otp'),
#   `dgwnsp.osp`= c( pmu * enet %n% 'dgwnsp.osp')
# )
# enet <- geli[[3]]
# pmu <- pmuncommon[[3]]
# status.t <- a.pmu.status[[3]]
# cossim.t <- a.pmu.cossim[[3]]
# cdf2 <- data.frame(
#   y=c(enet[,]),
#   status = c(status.t),
#   cossim = c( cossim.t ),
#   triangle= c( pmu * enet %n% 'triangle'),
#   ttriple=  c( pmu * enet %n% 'ttriple'),
#   cycle3= c( pmu * enet %n% 'cycle3'),
#   cycle4= c( pmu * enet %n% 'cycle4'),
#   nearsimmelian= c( pmu * enet %n% 'nearsimmelian'),
#   localtriangle= c( pmu * enet %n% 'localtriangle'),
#   `dgwesp.otp`= c( pmu * enet %n% 'dgwesp.otp'),
#   `dgwesp.osp`= c( pmu * enet %n% 'dgwesp.osp'),
#   `dgwnsp.otp`= c( pmu * enet %n% 'dgwnsp.otp'),
#   `dgwnsp.osp`= c( pmu * enet %n% 'dgwnsp.osp')
# )
# cdf <- rbind(cdf1,cdf2)
# round(cor(cdf), 3)

## SAVE PATCHED 
saveRDS(list(geli=geli, gcli=gcli), file = empnetsfile)

# ## PATCH FACTOR MARKET DIRECTED NETWORK BETWEENNESS CLUSTER
# for (i in 1:length(geli)) {
#   cat(sprintf('\n\n----------%s: %s------------\n',i,names(geli)[i]))
#   enet <- geli[[i]]
#   ge <- asIgraph(enet)
#   # cl <- edge.betweenness.community(ge, weights = NULL)
#   # enet %v% 'cluster_between' <- cl$membership
#   # enet %v% 'between' <- igraph::betweenness(ge, weights = NULL)
#   enet %v% 'degree_out' <- igraph::degree(ge, mode = 'out')
#   enet %v% 'degree_in' <- igraph::degree(ge, mode = 'in')
#   enet %v% 'degree_total' <- igraph::degree(ge, mode = 'total')
#   geli[[i]] <- enet
#   # # #
#   # cnet <- gcli[[i]]
#   # cl <- edge.betweenness.community(asIgraph(cnet), weights = NULL)$membership
#   # cnet %v% 'cluster_between' <- cl
#   # cnet %n% 'same_cluster_between' <- as.matrix(outer(cl, cl, '==') * 1)
#   # gcli[[i]] <- cnet
# }
# 
# ## SAVE PATCHED 
# saveRDS(list(geli=geli, gcli=gcli), file = empnetsfile)

##-------------------------
##  PLOT NETWORK COLOR CLUSTER
##_------------------------
par(mar=c(.1,.1,.1,.1))
.gx <- asIgraph(geli[[2]])
gxx <- induced_subgraph(.gx, which(igraph::degree(.gx)>0))
#
V(gxx)$vertex.color <- 'gray'
decomp <- decompose(gxx,mode = 'weak',min.vertices = 3)
lcc <- decomp[[which.max(sapply(decomp,vcount))]]
lccidx <- which(V(gxx)$vertex.names %in% V(lcc)$vertex.names)
# lcccomm <- igraph::multilevel.community(as.undirected(lcc))$membership
lcccomm <- igraph::edge.betweenness.community(lcc)$membership
nclust <- length(unique(lcccomm))
V(gxx)$cluster <- 0
V(gxx)$cluster[lccidx] <- lcccomm
V(gxx)$vertex.color[lccidx] <- rainbow(nclust,alpha = .6)[lcccomm]
#
plot.igraph(gxx, layout=layout.fruchterman.reingold(gxx, grid="nogrid"),
            vertex.color=V(gxx)$vertex.color,
            vertex.size=2.3,
            vertex.label=NA,
            vertex.label.cex = 0,
            edge.width=.5,
            edge.arrow.size=.3)
##--------------------------------------


##----------------------------------------
## FIT TERGM MODEL
##----------------------------------------
## EMPLOYEE NET ----------------------------
# dgwnsp.otp <- lapply(geli, function(net) as.matrix(net %n% 'dgwnsp_otp') )
# dgwnsp.itp <- lapply(geli, function(net) as.matrix(net %n% 'dgwnsp_itp') )
# dgwnsp.osp <- lapply(geli, function(net) as.matrix(net %n% 'dgwnsp_osp') )
# dgwnsp.isp <- lapply(geli, function(net) as.matrix(net %n% 'dgwnsp_isp') )
# intransitive <- lapply(geli, function(net) as.matrix(net %n% 'intransitive') )
# ctriple <- lapply(geli, function(net) as.matrix(net %n% 'ctriple') )
ttriple <- lapply(geli, function(net) as.matrix(net %n% 'ttriple') )
cycle3 <- lapply(geli, function(net) as.matrix(net %n% 'cycle3') )
# cycle4 <- lapply(geli, function(net) as.matrix(net %n% 'cycle4') )
# localtriangle <- lapply(geli, function(net) as.matrix(net %n% 'localtriangle') )
# triangle   <- lapply(geli, function(net) as.matrix(net %n% 'triangle') )
# transitive <- lapply(geli, function(net) as.matrix(net %n% 'transitive') )
# twopath    <- lapply(geli, function(net) as.matrix(net %n% 'twopath') )
# dgwnsp.otp <- lapply(geli, function(net) as.matrix(net %n% 'dgwnsp.otp') )
# dgwesp.osp <- lapply(geli, function(net) as.matrix(net %n% 'dgwesp.osp') )

# eigendiff <- lapply(geli, function(net) {
#   x <- net %v% 'eigen_centrality'
#   return( outer(x,x,'-') * -1 )   ## multiply by -1 for forward diff [to]-[from] 
# })
status <- lapply(geli, function(net){
  x <- net %v% 'eigen_centrality'
  return(outer(x,x,'+'))   ## change statistics are defined as sum{V_i, V_j}
})
#
samefmclust <- lapply(geli, function(net){
  x <- net %v% 'cluster_between'
  return(outer(x,x,'==') * 1)
})
# relposition <- lapply(geli, function(net){
#   x <- net %v% 'between'
#   return(outer(x,x,'-') * -1)  ## multiply by -1 for forward diff [to]-[from];  [central]-[peripheral] > 0
# })
#
# female <- lapply(geli, function(net) as.matrix(net %n% 'female_avg') )
# frexec <- lapply(geli, function(net) as.matrix(net %n% 'from_exec_role_avg') )
# toexec <- lapply(geli, function(net) as.matrix(net %n% 'to_exec_role_avg') )
tofrexecdiff <- lapply(geli, function(net) as.matrix(net %n% 'to_from_exec_role_diff_avg') )
# frtenure <- lapply(geli, function(net) as.matrix(net %n% 'from_tenure') )
# frtenure.sq <- lapply(geli, function(net) as.matrix(net %n% 'from_tenure')^2 )
agediff <- lapply(geli, function(net) as.matrix(net %n% 'age_diff') * -1 )  ## multiply by -1 for forward fiff [to]-[from]
# empweight <- lapply(geli, function(net) as.matrix(net %n% 'weight') )
## COMPETITION NET ----------------------------
mmc <- lapply(gcli, function(net) as.matrix(net %n% 'mmc') )
cossim <- lapply(gcli, function(net) as.matrix(net %n% 'cat_cos_sim') )
# invdist <- lapply(gcli, function(net) {
#   D <- 1 / as.matrix(igraph::distances(asIgraph(net)))
#   diag(D) <- 1
#   return(D) 
#   })
#
# entrepopp <- lapply(gcli, function(net){
#   x <- net %v% 'ipo_status'
#   return(outer(x,x,'-'))  ## don't need to multiple by -1; already backwards diff [from]-[to] == same as == -1 * forward diff [to]-[from]
# })
#
diffcl <- lapply(gcli, function(net) {
  x <- net %n% 'same_cluster_multilevel'
  return( 1 - x )  ## take opposite of same cluster membership
})
#
pmuncommon <- lapply(gcli, function(net) { # product market uncommonality dummy
  D <- as.matrix(igraph::distances(asIgraph(net)))
  PMU <- D
  diag(PMU) <- 0
  PMU[D == 1] <- 0
  PMU[D > 1] <- 1
  return(PMU)
})
# pmr <- lapply(gcli, function(net) {   # product market rival dummy
#   D <- as.matrix(igraph::distances(asIgraph(net)))
#   PMR <- D
#   diag(PMR) <- 0
#   PMR[D != 1] <- 0
#   return(PMR)
# })
# #  # Single matrix of 'same neighborhood' if have one same cateogry_group
### SLOW ------------------
# cgl <- gcli[[1]] %v% 'category_group_list'
# samecatgrp.mat <- outer(cgl, cgl, Vectorize(function(a,b){
#   x <- c(str_split(a,'[|]',simplify = T))
#   y <- c(str_split(b,'[|]',simplify = T))
#   return( as.integer(length(intersect(x,y)) > 0) )
# }))
##------------------------
npds <- length(geli)
# samecatgrp <- lapply(1:npds, function(i)  samecatgrp.mat )
# diffcatgrp <- lapply(1:npds, function(i)  1-samecatgrp.mat )
##---------- A. INTERACTIONS: COMPNET TIE-------------------------------
a.pmu.status             <- lapply(1:npds, function(i)  pmuncommon[[i]] * status[[i]] )
# a.pmu.samefmclust        <- lapply(1:npds, function(i)  pmuncommon[[i]] * samefmclust[[i]] )
a.pmu.cossim             <- lapply(1:npds, function(i)  pmuncommon[[i]] * cossim[[i]] )
# a.pmu.relposition        <- lapply(1:npds, function(i)  pmuncommon[[i]] * relposition[[i]] )
# a.pmu.triangle           <- lapply(1:npds, function(i)  pmuncommon[[i]] * triangle[[i]] )
# a.pmu.transitive         <- lapply(1:npds, function(i)  pmuncommon[[i]] * transitive[[i]] )
a.pmu.samefmclust        <- lapply(1:npds, function(i)  pmuncommon[[i]] * samefmclust[[i]] )
a.pmu.cycle3             <- lapply(1:npds, function(i)  pmuncommon[[i]] * cycle3[[i]] )
# a.pmu.cycle4             <- lapply(1:npds, function(i)  pmuncommon[[i]] * cycle4[[i]] )
# a.pmu.localtriangle      <- lapply(1:npds, function(i)  pmuncommon[[i]] * localtriangle[[i]] )
# a.pmu.dgwesp.otp         <- lapply(1:npds, function(i)  pmuncommon[[i]] * dgwesp.otp[[i]] )
# a.pmu.dgwesp.osp         <- lapply(1:npds, function(i)  pmuncommon[[i]] * dgwesp.osp[[i]] )
a.pmu.ttriple            <- lapply(1:npds, function(i)  pmuncommon[[i]] * ttriple[[i]] )
# a.pmu.intransitive       <- lapply(1:npds, function(i)  pmuncommon[[i]] * intransitive[[i]] )
# a.pmu.dgwnsp.otp         <- lapply(1:npds, function(i)  pmuncommon[[i]] * dgwnsp.otp[[i]] )
# a.pmu.dgwnsp.osp         <- lapply(1:npds, function(i)  pmuncommon[[i]] * dgwnsp.osp[[i]] )
##---------- B. INTERACTIONS: COMPNET DIFFERENT CLUSTER---------------------------
# b.diffcl.status          <- lapply(1:npds, function(i)  diffcl[[i]] * status[[i]] )
# b.diffcl.samefmclust     <- lapply(1:npds, function(i)  diffcl[[i]] * samefmclust[[i]] )
# b.diffcl.triangle        <- lapply(1:npds, function(i)  diffcl[[i]] * triangle[[i]] )
##---------- C. INTERACTIONS: DIFFERENT CB CATEGORY GROUP (none in common)------------------
# c.diffcatgrp.status      <- lapply(1:npds, function(i)  diffcatgrp[[i]] * status[[i]] )
# c.diffcatgrp.samefmclust <- lapply(1:npds, function(i)  diffcatgrp[[i]] * samefmclust[[i]] )
# c.diffcatgrp.triangle    <- lapply(1:npds, function(i)  diffcatgrp[[i]] * triangle[[i]] )

# vt <- matrix(rep(v,3),byrow=T, nrow=length(v))
# w * vt

##----BOOTSTRAP REPLICATATIONS-----
R <- 200
ncpus <- detectCores()
parallel <- "multicore"
##---------------------------------

#GOF-----------------------------------
# structvar <- 'cycle3_ttriple'
# filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
# m5a <- readRDS(sprintf('%s.rds',filebase))
# 
# gf.dsp <- btergm::gof(m5a, statistics = c(dsp), nsim=100, parallel = parallel, ncpus = ncpus)
# saveRDS(gf.dsp, file = sprintf('gof_dsp_%s.rds',filebase))
# rm(gf.dsp)
# gc()
# 
# gf.esp <- btergm::gof(m5a, statistics = c(esp), nsim=100, parallel = parallel, ncpus = ncpus)
# saveRDS(gf.dsp, file = sprintf('gof_esp_%s.rds',filebase))
# rm(gf.esp)
# gc()
# 
# gf.odeg <- btergm::gof(m5a, statistics = c(odeg), nsim=100, parallel = parallel, ncpus = ncpus)
# saveRDS(gf.dsp, file = sprintf('gof_odeg_%s.rds',filebase))
# rm(gf.odeg)
# gc()
# 
# gf.ideg <- btergm::gof(m5a, statistics = c(ideg), nsim=100, parallel = parallel, ncpus = ncpus)
# saveRDS(gf.dsp, file = sprintf('gof_ideg_%s.rds',filebase))
# rm(gf.ideg)
# gc()
# 
# gf.geodesic <- btergm::gof(m5a, statistics = c(geodesic), nsim=100, parallel = parallel, ncpus = ncpus)
# saveRDS(gf.dsp, file = sprintf('gof_geodesic_%s.rds',filebase))
# rm(gf.geodesic)
# gc()
# 
# rm(m5a)
# gc()
#-------------------------------------------------------

##----------------
## 0. baseline
##---------------- 
structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m0a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") +
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) + edgecov(mmc) +  nodematch("state_code", diff = F) + ## GEOGRAPHIC
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(pmuncommon) 
                ##---- INTERACTIONS -----------------
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()


##----------------
## H1. status
##---------------- 
structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m1a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") +  edgecov(status) +
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) + edgecov(mmc) +  nodematch("state_code", diff = F) + ## GEOGRAPHIC
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(pmuncommon) + 
                ##---- INTERACTIONS -----------------
              edgecov(a.pmu.status) 
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()


##----------------
## H2. COSSIM
##---------------- 
structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m2a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") +  
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) + edgecov(mmc) +  nodematch("state_code", diff = F) + ## GEOGRAPHIC
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
              edgecov(cossim) +
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(pmuncommon) + 
                ##---- INTERACTIONS -----------------
                edgecov(a.pmu.cossim) 
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()


##----------------
## H3 TRIADS
##---------------- 
structvar <- 'cycle3_ttriple'
filebase <- sprintf('employee_net_btergm_m3a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
              edgecov(cycle3) + edgecov(ttriple) +
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") + 
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) + edgecov(mmc) +  nodematch("state_code", diff = F) + ## GEOGRAPHIC
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(pmuncommon) + 
                ##---- INTERACTIONS -----------------
                edgecov(a.pmu.cycle3) + edgecov(a.pmu.ttriple) 
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()

##----------------
## A#. Product Market Uncommonality - cycle3 + cycle4
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

###############################################################################################

##----------------
## A#. Product Market Uncommonality - cycle3 + cycle4
##---------------- 
structvar <- 'cycle34'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
              edgecov(cycle3) + edgecov(cycle4) +
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
                edgecov(a.pmu.cycle3) + edgecov(a.pmu.cycle4) 
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()

##----------------
## A#. Product Market Uncommonality - cycle3 + cycle4
##---------------- 
structvar <- 'cycle34'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                 ##---- STRUCTURAL ----------------------------
               edgecov(cycle3) + edgecov(cycle4) +
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
                 edgecov(a.pmu.cycle3) + edgecov(a.pmu.cycle4) 
               , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()

##----------------
## A#. Product Market Uncommonality - localtriangle
##---------------- 
structvar <- 'localtriangle'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
              edgecov(localtriangle) +
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
                edgecov(a.pmu.localtriangle) 
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()

##----------------
## A#. Product Market Uncommonality - samefmclust
##---------------- 
structvar <- 'samefmclust'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
              edgecov(samefmclust) +
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
                edgecov(a.pmu.samefmclust)   
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a, file=sprintf('%s.rds',filebase))
rm(m5a)
gc()



##----------------
## A2. Product Market Uncommonality - TWOPATHS
##---------------- 
structvar <- 'intransitive'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a1 <- btergm(geli ~ edges + 
                 ##---- STRUCTURAL ----------------------------
               edgecov(intransitive) +
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
                 edgecov(a.pmu.intransitive)   
               , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a1),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a1), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a1, file=sprintf('%s.rds',filebase))
rm(m5a1)
gc()

##----------------
## A2. Product Market Uncommonality - TWOPATHS
##---------------- 
structvar <- 'twopath'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a2 <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
              edgecov(twopath) +
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
              edgecov(a.pmu.twopath)   
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a2),  digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a2), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a2, file=sprintf('%s.rds',filebase))
rm(m5a2)
gc()

##----------------
## A3. Product Market Uncommonality - GWESP.OTP
##---------------- 
structvar <- 'gwesp_otp'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a3 <- btergm(geli ~ edges + 
                 ##---- STRUCTURAL ----------------------------
               edgecov(dgwesp.otp) +
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
                 edgecov(a.pmu.dgwesp.otp)   
               , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a3), digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a3), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a3, file=sprintf('%s.rds',filebase))
rm(m5a3)
gc()


##----------------
## A4. Product Market Uncommonality - TRANSITIVE 
##---------------- 
structvar <- 'transitive'
filebase <- sprintf('employee_net_btergm_m5a_%s_y%s-%s_R%s',structvar,startYr,endYr,R)
m5a4 <- btergm(geli ~ edges + 
                 ##---- STRUCTURAL ----------------------------
               edgecov(transitive) +
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
                 edgecov(a.pmu.transitive)   
               , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a4), digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
htmlreg(list(m5a4), file=sprintf('%s.html',filebase), digits = 3, single.row = T)
saveRDS(m5a4, file=sprintf('%s.rds',filebase))
rm(m5a4)
gc()



# ## COEF NAME MAP
# coef.name.map <- list(
#   `edges`='Edges',
#   # EMPLOYEE (HUMAN CAPITAL)
#   `edgecov.tofrexecdiff[[i]]`='Executive Role Mobility (-1,0,1)',
#   `edgecov.to.fr.exec.diff[[i]]`='Executive Role Mobility (-1,0,1)',
#   `edgecov.frtenure[[i]]`='From Job Tenure (years)',
#   `edgecov.frexec[[i]]`='From Executive Role',
#   `edgecov.toexec[[i]]`='To Executive Role',
#   `edgecov.female[[i]]`='Female Ratio',
#   # FIRM
#   `nodecov.age`='Firms Age',
#   `nodecov.ipo_status`='Ownership Status (1=IPO)',
#   # DYAD
#   `nodematch.ipo_status`='Ownership Status Homophily',
#   `nodematch.state_code`='HQ Region Homophily',
#   `edgecov.agediff[[i]]`='Age Difference',
#   `edgecov.mmc[[i]]`='Firm Branch MMC',
#   # DYAD - PRODUCT MARKET
#   `edgecov.invdist[[i]]`='Product Market Closeness (inverse distance)',
#   # Dyad - FACTOR MARKET
#   `edgecov.cossim[[i]]`='Category Similarity',
#   # STRUCTURE
#   `ctriple`='Cyclic Triads',
#   `ttriple`='Transitive Triads',
#   # TEMPORAL
#   `edgecov.delrecip[[i]]`='Delayed Reciprocity',
#   `edgecov.memory[[i]]`='Temporal Stability'
# )


##----------------
## A. Product Market Uncommonality
##---------------- 
m5a <- btergm(geli ~ edges + 
                ##---- STRUCTURAL ----------------------------
              # odegree +   ## total degree
                edgecov(triangle) +
              # edgecov(transitive) +
                # edgecov(ctriple) +  edgecov(ttriple) + #localtriangle(catgroupsame) +  #mutual
                # edgecov(dgwnsp.otp) + edgecov(dgwnsp.osp) + 
                # edgecov(intransitive) +
                # gwodegree(0, fixed=T) + #gwodegree(0, fixed=T) + #gwesp(0, fixed=T) +
                # dgwesp(decay=0, fixed=T, cutoff=20, type="OTP") #+  ##transitive shared partners
                # dgwesp(decay=0, fixed=T, cutoff=30, type="ITP") +  ##cyclical shared partners
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                # timecov(transform = function(t) t) +
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                # edgecov(frexec) +
                # edgecov(frexec) +  edgecov(toexec) +  ###edgecov(female) +
                # edgecov(frtenure) + #edgecov(frtenure.sq) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") + 
                # nodecov("degree_out") +  #nodecov("degree_in") +
              edgecov(status) +
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) +
                nodematch("state_code", diff = F) + ## GEOGRAPHIC
                edgecov(mmc) +
                # edgecov(samefmclust) +  # edgecov(eigendiff) + 
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
              edgecov(cossim) +
                # edgecov(samecatgrp) + 
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(pmuncommon) + 
                #edgecov(diffcl) + 
                #edgecov(diffcatgrp) +
                ##---- INTERACTIONS -----------------
              edgecov(a.pmu.status) + 
              # edgecov(a.pmu.samefmclust) + 
              edgecov(a.pmu.cossim) + 
              edgecov(a.pmu.triangle) # +
              # edgecov(a.pmu.transitive) # +
              # edgecov(a.pmu.ctriple) + 
              # edgecov(a.pmu.ttriple) 
              # edgecov(a.pmu.intransitive)    
              , R=R, parallel = parallel, ncpus = ncpus)
screenreg(list(m5a), digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
saveRDS(m5a, file=sprintf('employee_net_btergm_m5a_triangle_y%s-%s_R%s.rds',startYr,endYr,R))
rm(m5a)
gc()



##----------------
## B. DIFFERENT PRODUCT MARKET CLUSTER
##---------------- 
m5b <- btergm(geli ~ 
                ##---- STRUCTURAL ----------------------------
              edges + 
                edgecov(triangle) + # edgecov(ttriple) + #localtriangle(catgroupsame) +  #mutual
                # gwodegree(0, fixed=T) + #gwodegree(0, fixed=T) + #gwesp(0, fixed=T) +
                # dgwesp(decay=0, fixed=T, cutoff=20, type="OTP") #+  ##transitive shared partners
                # dgwesp(decay=0, fixed=T, cutoff=30, type="ITP") +  ##cyclical shared partners
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                # timecov(transform = function(t) t) +
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                # edgecov(frexec) +
                # edgecov(frexec) +  edgecov(toexec) +  ###edgecov(female) +
                # edgecov(frtenure) + #edgecov(frtenure.sq) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") + 
                edgecov(status) +
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) +
                nodematch("state_code", diff = F) + ## GEOGRAPHIC
                edgecov(mmc) +
                edgecov(samefmclust) +  # edgecov(eigendiff) + 
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
              edgecov(cossim) + 
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              # edgecov(pmuncommon) + 
                edgecov(diffcl) +
                #edgecov(diffcatgrp) +
                ##---- INTERACTIONS -----------------
              edgecov(b.diffcl.status) + 
                edgecov(b.diffcl.samefmclust) + 
                edgecov(b.diffcl.triangle) 
              , R=100)
screenreg(list(m5b), digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
saveRDS(m5b, file='employee_net_btergm_m5b.rds')
rm(m5b)
gc()



##----------------
## C. DIFFERENT CATEGORY GROUP (none in common)
##---------------- 
m5c <- btergm(geli ~ 
                ##---- STRUCTURAL ----------------------------
              edges + 
                edgecov(triangle) + # edgecov(ttriple) + #localtriangle(catgroupsame) +  #mutual
                # gwodegree(0, fixed=T) + #gwodegree(0, fixed=T) + #gwesp(0, fixed=T) +
                # dgwesp(decay=0, fixed=T, cutoff=20, type="OTP") #+  ##transitive shared partners
                # dgwesp(decay=0, fixed=T, cutoff=30, type="ITP") +  ##cyclical shared partners
                ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
                # timecov(transform = function(t) t) +
                ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
                # edgecov(frexec) +
                # edgecov(frexec) +  edgecov(toexec) +  ###edgecov(female) +
                # edgecov(frtenure) + #edgecov(frtenure.sq) +
                ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") + 
                edgecov(status) +
                ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) +
                nodematch("state_code", diff = F) + ## GEOGRAPHIC
                edgecov(mmc) +
                edgecov(samefmclust) +  # edgecov(eigendiff) + 
                ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
              edgecov(cossim) + 
                ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              # edgecov(pmuncommon) + 
                #edgecov(diffcl) + 
                edgecov(diffcatgrp) +
                ##---- INTERACTIONS -----------------
              edgecov(c.diffcatgrp.status) + 
                edgecov(c.diffcatgrp.samefmclust) + 
                edgecov(c.diffcatgrp.triangle) 
              , R=100)
screenreg(list(m5c), digits = 3, single.row = T)  ##,  custom.coef.map = coef.name.map
saveRDS(m5c, file='employee_net_btergm_m5c.rds')
rm(m5c)
gc()



# female 
# frexec 
# toexec 
# tofrexecdiff
# frtenure 
# agediff 
# mmc 
# cossim 
# invdist 
# pmuncommon 

mod9 <- btergm(geli ~ 
              ##---- STRUCTURAL ----------------------------
              edges + ctriple + ttriple + #localtriangle(catgroupsame) +  #mutual 
              # gwodegree(0, fixed=T) + #gwodegree(0, fixed=T) + #gwesp(0, fixed=T) +
              # dgwesp(decay=0, fixed=T, cutoff=20, type="OTP") #+  ##transitive shared partners
              # dgwesp(decay=0, fixed=T, cutoff=30, type="ITP") +  ##cyclical shared partners
              ##---- TEMPORAL -----------------
              memory(type = "stability",lag=1) + delrecip(mutuality = TRUE,lag=1) + 
              # timecov(transform = function(t) t) +
              ##---- INDIVIDUAL / HUMAN CAPITAL -----------------
              edgecov(tofrexecdiff) +
              # edgecov(frexec) +
              # edgecov(frexec) +  edgecov(toexec) +  
              edgecov(female) +
              # edgecov(frtenure) + #edgecov(frtenure.sq) +
              ##---- FIRM -----------------
              nodecov("age") + nodecov("ipo_status") + 
              ##---- FIRM DYAD - CONTROL -----------------
              edgecov(agediff) +
              nodematch("ipo_status", diff = F) + 
              nodematch("state_code", diff = F) + ## GEOGRAPHIC
              edgecov(mmc) +
              ##---- FIRM DYAD - PRODUCT MARKET COMP -----------------
              edgecov(invdist) +
              ##---- FIRM DYAD - RESOURCE RELATEDNESS -----------------
              edgecov(cossim)
               , R=100)
screenreg(list(mod9), digits = 3, single.row = T)  # custom.coef.map = coef.name.map
saveRDS(mod9, file='employee_net_btergm_mod9.rds')



## roles c('director', '(vp|svp|evp)', '(founder|co-founder)', '(manager|managing)' )
## 
## 

## job.all
.titles <- na.omit(c(job.all$title1,job.all$title2))
ja.ti <- plyr::count(.titles)
ja.ti <- ja.ti[order(ja.ti$freq, decreasing = T), ]
dim(ja.ti)
View(head(ja.ti,100))

##========================================================
## Employee Distribution
##--------------------------------------------------------



## started on year
job.all$started_on_year <- as.numeric(substr(job.all$started_on, 1, 4))
## job year count
jycnt <- plyr::count(job.all$started_on_year)




## SUBSET COMPANY EGO NET
name_i = 'microsoft'
d <- 2
g.sub <- igraph::make_ego_graph(g.full, order=d, nodes=V(g.full)[V(g.full)$name==name_i])[[1]]

co.uuid <- V(g.sub)$company_uuid

job <- job.all[job.all$org_uuid %in% co.uuid, ]
## filter job years 2010-2015
years <- 2009:2015
job <- job[job$started_on_year %in% years, ]

pcnt <- plyr::count(job$person_uuid)
pcnt <- pcnt[order(pcnt$freq, decreasing=T),]

ocnt <- plyr::count(job$org_uuid)
ocnt <- ocnt[order(ocnt$freq, decreasing=T),]

cat(sprintf('people: cnt %d, frac>1 %.3f, min %d, mid %d, mean %.3f, max %d',
            nrow(pcnt), length(pcnt$freq[pcnt$freq>1])/nrow(pcnt), 
            min(pcnt$freq), median(pcnt$freq), mean(pcnt$freq), max(pcnt$freq)))
cat(sprintf('companies: cnt %d, frac>1 %.3f, min %d, mid %d, mean %.3f, max %d',
            nrow(ocnt), length(ocnt$freq[ocnt$freq>1])/nrow(ocnt), 
            min(ocnt$freq), median(ocnt$freq), mean(ocnt$freq), max(ocnt$freq)))


## make association matrix for 





## Jobs
job.f <- merge(co[, c('company_name_unique','company_uuid')],job,
               by.x='company_uuid',by.y='org_uuid',all=F)
jb <- merge(ppl[, c('uuid','person_name_unique')], job.f,
            by.x='uuid',by.y='person_uuid',all=F)

pcnt <- plyr::count(jb$person_name_unique)

ocnt <- plyr::count(jb$company_uuid)















