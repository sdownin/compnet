setwd("C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2")
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(parallel)
library(statnet, quietly = T)
library(network, quietly = T)
library(xergm, quietly = T)  ## includes rem, tnam, GERGM
library(texreg, quietly = T)
library(igraph, quietly = T)
library(plyr, quietly = T)
library(lattice, quietly = T)
library(latticeExtra, quietly = T)
library(ggplot2, quietly = T)
library(reshape2)

data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/"
## FUNCTIONS
source(file.path(getwd(),'R','comp_net_functions.R'))
source(file.path(getwd(),'R','netrisk_functions.R'))
## DATA
source(file.path(getwd(),'R','cb_data_prep.R'))


## comptetition network
## 37828
g.full <- read.graph('g_full.graphml', format='graphml')

## All acquisitions
## 29805
# co_acq

## acquisition filtered by comp net firms union (either aquirer or acquired)
## 14208
co_acq_g_u <- co_acq[which(co_acq$acquirer_name_unique %in% V(g.full)$name
                          | co_acq$acquiree_name_unique %in% V(g.full)$name), ]
## acquisition filtered by comp net firms intersection (both aquirer and acquired)
## 3442
co_acq_g_i <- co_acq[which(co_acq$acquirer_name_unique %in% V(g.full)$name
                         & co_acq$acquiree_name_unique %in% V(g.full)$name), ]
## acquisition filtered by comp net firms intersection, recent date (>=2010)
## 2786
co_acq_g_i_d <- co_acq_g_i[which(co_acq_g_i$acquired_on >= '2010-01-01'), ]

## check top filtered acquirers
cnt <- plyr::count(co_acq_g_i_d$acquirer_name_unique)
cnt <- cnt[order(cnt$freq, decreasing = T),]
head(cnt, 20)
# Top Acquirers:             x   freq
# 496                   google   73
# 544                      ibm   46
# 414                 facebook   32
# 847                   oracle   30
# 1358                   yahoo   29
# 80                     apple   28
# 740                microsoft   26
# 1002              salesforce   26
# 1229                 twitter   23
# 58                    amazon   20
# 240                    cisco   18
# 368                     ebay   17
# 505                  groupon   17
# 581                    intel   16
# 68                       aol   15
# 385  endurance-international   15
# 945                  rakuten   14
# 243           citrix-systems   13
# 323                     dell   12
# 379                      emc   12

## 


# Top Acquirers:         x  freq
# 5343              Google  204
# 2656               Cisco  197
# 8070           Microsoft  186
# 6040                 IBM  170
# 13985             Yahoo!  120
# 9234  Oracle Corporation  113
# 5779     Hewlett-Packard   99
# 941                Apple   82
# 6449               Intel   82
# 9435     Parker Hannifin   74
# 4093                 EMC   73
# 890                  AOL   68
# 712               Amazon   66
# 48                    3M   65
# 4528            Facebook   64
# 11968           Symantec   60
# 3923                eBay   59
# 1338               Avnet   54
# 12908            Twitter   48
# 1596  Berkshire Hathaway   46

## Competition Network (d=2) Size:
# apple      1101 2187
# google     2100 4598
# IBM         828 1710 
# Microsoft  1873 3925
# cisco       220  406
# oracle     1248 2456

## FIND FOCAL FIRM 
name_i <- 'cisco'
d <- 2

g.ego <- igraph::make_ego_graph(graph = g.full, 
                                nodes = V(g.full)[V(g.full)$name==name_i], 
                                order = d, mode = 'all')[[1]]
g.ego

par(mar=c(.1,.1,.1,.1))
deg=log(1+igraph::degree(g.ego))
plot(g.ego, vertex.size=deg*4, vertex.label.cex=deg*.3)

##--------------------------------------------------------------
##--------------------------------------------------------------
##--------- CREATE FIRM NETWORK PERIOD LISTS  ------------------
##--------------------------------------------------------------
##--------------------------------------------------------------

# ## cache original
# firm.nets.orig <- firm.nets



## run main network period creation loop
for (i in 1:length(firms.todo)) {
  ## -- settings --
  d <- 2
  yrpd <- 1
  startYr <- 2007
  endYr <- 2017
  ## --------------
  name_i <- firms.todo[i]
  cat(sprintf('\n---------%s----------\n',name_i))
  periods <- seq(startYr,endYr,yrpd)
  company.name <- 'company_name_unique'
  verbose <- TRUE
  #
  #g.base <- igraph::make_ego_graph(g.full,order=k,nodes=V(g.full)[V(g.full)$name=='surveymonkey'])[[1]]
  g.base <- g.full
  g.k.sub <- igraph::make_ego_graph(graph = g.base, nodes = V(g.full)[V(g.full)$name==name_i], 
                                    order = d, mode = 'all')[[1]]
  net.k.sub <- getNetFromIgraph(g.k.sub)
  net.k.sub %n% 'ego' <- name_i
  net <- net.k.sub
  #----------------Network List-------------------
  nl <- list()
  for (t in 2:length(periods)) {
    cat(sprintf('\nmaking period %s-%s:\n', periods[t-1],periods[t]))
    tmp.net <- makePdNetwork(net.k.sub,
                             start=periods[t-1], end=periods[t])
    nl[[t]] <- setCovariates(tmp.net, periods[t-1], periods[t],
                             netRiskCommunityAlgo='multilevel.community',
                             downweight.env.risk=FALSE,
                             acq=co_acq,br=co_br,rou=co_rou,ipo=co_ipo)
  }
  nl.bak <- nl
  nl <- nl[which(sapply(nl, length)>0)]
  names(nl) <- periods[2:length(periods)]
  ## ---------- add LAGS ----------------
  for (t in 2:length(nl)) {
    nl[[t]] %n% 'DV_lag' <- nl[[t-1]][,]
    nl[[t]] %n% 'dist_lag' <- nl[[t-1]] %n% 'dist'
    ##-------------------------------------------
    # nl[[t]] %v% 'net_risk_lag' <- nl[[t-1]] %v% 'net_risk'
    # nl[[t]] %v% 'cent_deg_lag' <- nl[[t-1]] %v% 'cent_deg'
    # nl[[t]] %v% 'genidx_multilevel_lag' <- nl[[t-1]] %v% 'genidx_multilevel'
    # nl[[t]] %v% 'cent_pow_n0_5_lag' <- nl[[t-1]] %v% 'cent_pow_n0_5'
    # nl[[t]] %v% 'cent_pow_n0_1_lag' <- nl[[t-1]] %v% 'cent_pow_n0_1'
    # nl[[t]] %v% 'cent_eig_lag' <- nl[[t-1]] %v% 'cent_eig'
    # nl[[t]] %v% 'cent_deg_lag' <- nl[[t-1]] %v% 'cent_deg'
    # g.tmp <- getIgraphFromNet(nl[[t]])
    # if (vcount(g.tmp)>0 & ecount(g.tmp)>0) {
    #   nl[[t]] %v% 'constraint' <- igraph::constraint(g.tmp)
    # }
  }
  ##--------------- GET TERGM NETS LIST -----------
  ## only nets with edges > 0
  nets.all <- nl[2:length(nl)]
  nets <- nets.all[ which(sapply(nets.all, getNetEcount) > 0) ]
  #-------------------------------------------------
  
  # ## SAVE variable in image
  # # firm.nl <- list()
  # firm.nets[[net_group]][[name_i]] <- nets
  
  # ## CAREFUL TO OVERWRITE
  # file.name <- sprintf('tergm_firm_nets_1yr_6pd_v4_%s.rds',net_group)
  # saveRDS(firm.nets, file=file.name)
  ## CAREFUL TO OVERWRITE
  saveRDS(nets, file=sprintf('firm_nets_cem/%s_d%s.rds', name_i, d))
  gc()
}



##--------------------------------------------------------------
##--------------------------------------------------------------
##--------- CREATE FIRM NETWORK PERIOD LISTS  ------------------
##--------- SAVE EACH PD NET SEPARATELY       ------------------
##--------------------------------------------------------------
##--------------------------------------------------------------
firms.todo <- 'qualtrics'

## -- settings --
d <- 2   ## large graph > 1000 nodes, save years separately for one firm
yrpd <- 1
startYr <- 2010
endYr <- 2017

## --------------
i <- 1
name_i <- firms.todo[i]
cat(sprintf('\n---------%s----------\n',name_i))
periods <- seq(startYr,endYr,yrpd)
company.name <- 'company_name_unique'
verbose <- TRUE
#
#g.base <- igraph::make_ego_graph(g.full,order=k,nodes=V(g.full)[V(g.full)$name=='surveymonkey'])[[1]]
g.base <- g.full
g.k.sub <- igraph::make_ego_graph(graph = g.base, nodes = V(g.full)[V(g.full)$name==name_i], 
                                  order = d, mode = 'all')[[1]]
net.k.sub <- getNetFromIgraph(g.k.sub)
net.k.sub %n% 'ego' <- name_i
net <- net.k.sub
#----------------Network List-------------------
for (t in 2:length(periods)) {
  nl <- list()
  pd <- as.character(periods[t])
  cat(sprintf('\nmaking period %s-%s:\n', periods[t-1],periods[t]))
  tmp.net <- makePdNetwork(net.k.sub,
                           start=periods[t-1], end=periods[t])
  cat(sprintf("\n nodes : %s \n",nrow(tmp.net[,])))
  nl[[pd]] <- setCovariates(tmp.net, periods[t-1], periods[t],
                            netRiskCommunityAlgo='multilevel.community',
                            downweight.env.risk=FALSE,
                            acq=co_acq,br=co_br,rou=co_rou,ipo=co_ipo)
  
  ## CAREFUL TO OVERWRITE
  saveRDS(nl, file=sprintf('firm_nets_cem/%s_d%s_pd%s.rds', name_i, d, periods[t]))
  rm(nl)
  gc()
}
################################################################
################################################################
