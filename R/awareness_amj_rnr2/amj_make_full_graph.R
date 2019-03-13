##--------------------------------------------------------------
##
##        Make Full Global Competition Netowrk
##
##--------------------------------------------------------------
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(igraph)
library(readxl)


##
#
##
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

  # graph filename
  g.full.file <- file.path(net_dir,'g_full.graphml')
  
  ## load full graph, else make full graph if not exists in working directory
  if (file.exists(g.full.file)) 
  {
    
    cat('\nloading full graph...')
    g.full <- read.graph(g.full.file, format='graphml')
    cat('done.')
    return(g.full)
    
  } 
    

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
  
  ## save graph file
  g.full.file <- file.path(net_dir,'g_full.graphml')
  igraph::write.graph(graph = g.full, file=g.full.file, format = 'graphml')
  
  return(g.full)  
  
  
}

##
# Export
##
.make.full.graph()





