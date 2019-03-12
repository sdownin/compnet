##--------------------------------------------------------------
##
##        Make Full Global Competition Netowrk
##
##--------------------------------------------------------------
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(igraph)
library(readxl)


.make.full.graph <- function(x=NA)
{
  
  
  ## DIRECTORIES
  data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
  work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
  img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
  version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')
  net_dir <- file.path(work_dir,'firm_nets_rnr2')
  sup_data_dir <- file.path(work_dir,'amj_rnr2_sup_data')  ## supplmental data dir
  
  ## set woring dir
  setwd(work_dir)
  
  source(file.path(version_dir,'amj_awareness_functions.R'))
  source(file.path(version_dir,'amj_cb_data_prep.R'))           ## cb : CrunchBase
  source(file.path(version_dir,'amj_sdc_coop.R'))               ## sdc: Thompson SDC
  source(file.path(version_dir,'amj_firm_size_controls.R'))     ## mi : mergent intellect
  source(file.path(version_dir,'amj_institutional_holdings.R')) ## ih : institutional holdings

  # graph filename
  g.full.file <- file.path(net_dir,'g_full.graphml')
  
  ## load full graph, else make full graph if not exists in working directory
  if (file.exists(g.full.file)) 
  {
    
    cat('\nloading full graph...')
    g.full <- read.graph(g.full.file, format='graphml')
    cat('done.')
    return(g.full)
    
  } else {
    
    cat('\nmaking full graph...')
    
    max.year <- 2016
    
    ## delete edges at or later than this date (the year after max.year)
    exclude.date <- sprintf('%d-01-01', max.year+1)
    
    ## make graph
    g.full <- aaf$makeGraph(comp = cb$co_comp, vertdf = cb$co)
    
    ## cut out confirmed dates >= 2016
    g.full <- igraph::induced.subgraph(g.full, vids=V(g.full)[which(V(g.full)$founded_year <= max.year
                                                                    | is.na(V(g.full)$founded_year)
                                                                    | V(g.full)$founded_year=='' ) ] )
    g.full <- igraph::delete.edges(g.full, E(g.full)[which(E(g.full)$relation_created_at >= exclude.date)])
    
    ## SIMPLIFY
    g.full <- igraph::simplify(g.full, remove.loops=T,remove.multiple=T,
                               edge.attr.comb = list(weight='sum',
                                                     relation_began_on='max',
                                                     relation_ended_on='min'))
    
    ## DELETE NODES MANUALLY CHECKED
    dfck <- read_excel(file.path(sup_data_dir, 'private_firm_financials_all_firm_networks_delete_nodes.xlsx'),
                       sheet = 1,  na = c('','-',"'-"))
    dfck <- dfck[which(dfck$MI_note=='DELETE'), ]
    g.full <- igraph::induced.subgraph(g.full, which( ! V(g.full)$name %in% dfck$firm))
    
    ##------------------------------------------------------
    ##-------preprocess parent-subsidiary relationships-----
    ##----------node collapse like acquisitions-------------
    ##------------------------------------------------------
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
    gfuuid <- V(g.full)$company_uuid
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
    g.full <- aaf$nodeCollapseGraph(g.full, par.chi.nc, remove.isolates=T, verbose = T)
    
    cat('done.\n')
    
    ##============================================
    ## ADD MANUAL UPDATES (EDGES | NODES)
    ##--------------------------------------------
    
    cat('\nmanually updating competitive relations in full graph...')
    
    ## add Qualtrics - Medallia competitive relation
    vid1 <- which(V(g.full)$company_uuid=='2f6ed0df-e019-f0ad-10bc-d7eee4710103')  ## qualtrics
    vid2 <- which(V(g.full)$company_uuid=='405c6579-fce0-ff76-6870-aa0236bafde7')  ## medallia
    edgeAttrs <- list(weight=1, relation_began_on=max('2002-01-01'), relation_ended_on=NA)
    g.full <- igraph::add.edges(g.full, c(vid1,vid2), attr = edgeAttrs)
    g.full <- igraph::simplify(g.full, remove.loops=T,remove.multiple=T,
                               edge.attr.comb = list(weight='sum',
                                                     relation_began_on='max',
                                                     relation_ended_on='min'))
    V(g.full)$weight <- 1
  
    cat('done.\n')
    
    return(g.full)  
  
  }
  
}


g.full <- .make.full.graph()

## save graph file
g.full.file <- file.path(net_dir,'g_full.graphml')
igraph::write.graph(graph = g.full, file=g.full.file, format = 'graphml')



