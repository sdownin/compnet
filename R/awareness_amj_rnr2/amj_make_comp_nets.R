##--------------------------------------------------------------
##
##  AMJ 2018 SPECIAL ISSUE 
##  CREATE FIRM COMPETITION NETWORK PERIODS AS LIST OBJECTS
##
##--------------------------------------------------------------
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(igraph)
library(intergraph)

## DIRECTORIES
data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')
net_dir <- file.path(work_dir,'firm_nets_rnr2')
sup_data_dir <- file.path(work_dir,'amj_rnr2_sup_data')  ## supplmental data dir

## set woring dir
setwd(work_dir)

aaf <- source(file.path(version_dir,'amj_awareness_functions.R'))$value
cb  <- source(file.path(version_dir,'amj_cb_data_prep.R'))$value           ## cb : CrunchBase
sdc <- source(file.path(version_dir,'amj_sdc_coop.R'))$value               ## sdc: Thompson SDC
si  <- source(file.path(version_dir,'amj_firm_size_controls.R'))$value     ## si : size controls from mergent intellect
ih  <- source(file.path(version_dir,'amj_institutional_holdings.R'))$value ## ih : institutional holdings
g.full <- source(file.path(version_dir,'amj_make_full_graph.R'))$value     ## g.full

# ## set firms to create networks (focal firm or replication study focal firms)
# firms.todo <- c('qualtrics','cloudcherry',
#                 'abroad101','checkmarket','clarabridge',
#                 'confirmit','customergauge','cx-index','empathica',
#                 'feedback-lite','first-mile-geo','getfeedback',
#                 'inqwise','leaderamp', 'medallia','myfeelback',
#                 'promoter-io','satmetrix',
#                 'snap-surveys-ltd','super-simple-survey','survata',
#                 'surveygizmo','surveymonkey',
#                 'surveyrock','typeform','userate','verint','voice-polls')
     
## set firms to create networks (focal firm or replication study focal firms)
firms.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')


# firms.todo <- c('facebook')
# firms.todo <- c('cnnmoney','fox-business-network','bloomberg',
#                 'hearstcorporation','newscorporation')

## -- settings --
d <- 3
yrpd <- 1
startYr <- 2005
endYr <- 2017            ## dropping first for memory term; actual dates 2007-2016
lg.cutoff <- 1100        ## large network size cutoff to save periods seprately 
force.overwrite <- FALSE ## if network files in directory should be overwritten
## --------------  


##
## run main network period creation loop
##
# for (i in 1:length(firms.todo)) {
for (i in 1:1) {

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
  dfpar <- read_excel(file.path(sup_data_dir, 'private_firm_financials_all_firm_networks_delete_nodes.xlsx'),
                      sheet = 1,  na = c('','-',"'-"))
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
  g.d.sub <- aaf$nodeCollapseGraph(g.d.sub, par.chi.nc, remove.isolates=T, verbose = T)
  
  cat('done.\n')
  
  ##_-----------------------------------------------------
  ##-------process pre-start-year acquisitions------------
  ##------------------------------------------------------
  acqs.pd <- cb$co_acq[cb$co_acq$acquired_on <= sprintf('%d-12-31',startYr-1), ]
  g.d.sub <- aaf$nodeCollapseGraph(g.d.sub, acqs.pd, remove.isolates=T, verbose = T)
  net.d.sub <- asNetwork(g.d.sub)
  cat(sprintf('v = %d, e = %d\n',vcount(g.d.sub),ecount(g.d.sub)))
  
  # ## subset to firms with employees count > 10
  # idx.employee <- which( !(V(g.d.sub)$employee_count %in% c('NA','-','1-10')) )
  # g.d.sub <- igraph::induced.subgraph(g.d.sub, vids = V(g.d.sub)[idx.employee])
  # cat(sprintf('filtered >10 employee count: v = %d, e = %d\n',vcount(g.d.sub),ecount(g.d.sub)))
  
  ##------------Network Time Period List--------------------
  nl <- list()
  
  for (t in 2:length(periods)) 
  {
    ## period dates
    cat(sprintf('\nmaking period %s-%s:\n', periods[t-1],periods[t]))
    t1 <- sprintf('%d-01-01',periods[t-1]) ## inclusive start date 'YYYY-MM-DD'
    t2 <- sprintf('%d-12-31',periods[t-1]) ## inclusive end date 'YYYY-MM-DD'
    
    ## check if period network file exists (skip if not force overwrite)
    file.rds <- file.path(net_dir,sprintf('%s_d%d_y%s.rds',name_i,d,periods[t-1]))
    if (!force.overwrite & file.exists(file.rds)) {
      cat(sprintf('file exists: %s\nskipping.\n', file.rds))
      next
    }
    
    ## period years indicate [start, end) -- start inclusive; end exclusive 
    
    ## 1. Node Collapse acquisitions within period
    acqs.pd <- cb$co_acq[cb$co_acq$acquired_on >= t1 & cb$co_acq$acquired_on <= t2, ]
    g.d.sub <- aaf$nodeCollapseGraph(g.d.sub, acqs.pd, verbose = T)
    
    ## 2. Subset Period Network
    nl[[t]] <- aaf$makePdNetwork(asNetwork(g.d.sub), periods[t-1], periods[t], isolates.remove = F) 
    
    ## 3. Set Covariates for updated Period Network
    nl[[t]] <- aaf$setCovariates(nl[[t]], periods[t-1], periods[t],
                                 acq=cb$co_acq,br=cb$co_br,rou=cb$co_rou,ipo=cb$co_ipo,
                                 coop=coop, ih=ih, size=size)
    
    ## save each period if large network (would exceed memory as full list of time periods)
    if (vcount(g.d.sub) >= lg.cutoff) {
      saveRDS(nl[[t]], file = file.rds)
      nv <- length(nl[[t]]$val)
      names(nv)[1] <- as.character(periods[t-1])
      write.csv(nv, file = file.path(net_dir, sprintf('%s_d%s.csv',name_i,d)),append = TRUE)
      nl[[t]] <- NULL ## remove from memory
    }
    
  }
  
  ##---------Small Networks: clean period list and save whole -----------
  if (vcount(g.d.sub) < lg.cutoff) 
  {
    ## ----drop null and skipped periods----
    nl.bak <- nl
    nl <- nl[which(sapply(nl, length)>0)]
    
    if (length(nl) > 1) {
      names(nl) <- periods[2:length(periods)]
    }
    
    # ## ---------- add LAGS ----------------
    # if (length(nl) > 1) {
    #   for (t in 2:length(nl)) { 
    #     nl[[t]] %n% 'DV_lag' <- nl[[t-1]][,]
    #   }
    # }  
    
    ##--------------- GET TERGM NETS LIST -----------
    ## only nets with edges > 0
    if (length(nl) > 1) {
      nets.all <- nl[2:length(nl)]
    } else {
      nets.all <- nl
    }
    nets <- nets.all[ which(sapply(nets.all, aaf$getNetEcount) > 0) ]
    ## record network sizes
    write.csv(sapply(nets,function(x)length(x$val)), file = file.path(net_dir, sprintf('%s_d%s.csv',name_i,d)))
    
    #-------------------------------------------------
    
    ## CAREFUL TO OVERWRITE 
    saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
    
  }
  
}










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
