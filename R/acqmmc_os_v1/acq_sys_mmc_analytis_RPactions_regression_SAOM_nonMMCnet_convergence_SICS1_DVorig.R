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
library(parallel)
library(RSiena)
library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(lubridate)

## DIRECTORIES
data_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
cs_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2/compustat"
rvn_dir <- 'C:/Data/RavenPack'
version_dir <- file.path(work_dir,'R','acqmmc_os_v1')
net_dir <- file.path(work_dir,'acqmmc_os_v1','data')
result_dir <- file.path(work_dir,'acqmmc_os_v1','data')
sup_data_dir <- file.path(work_dir,'acqmmc_os_v1','sup_data')  ## supplmental data dir

## set woring dir
setwd(work_dir)


## LOAD DATA AND DEPENDENCIES
cb     <- source(file.path(version_dir,'acq_cb_data_prep.R'))$value 
acf    <- source(file.path(version_dir,'acq_compnet_functions.R'))$value 


###
## Compute HHI 
###
hhi <- function(x, na.rm=FALSE){
  if (na.rm) {
    x <- x[ which(!is.nan(x) & !is.na(x)) ]
  }
  z <- 100 * x / sum(x)
  return(sum(z^2))
}

##
# Plot from array of adjacency matrices
##
arrPlot <- function(arr, t=1, min.deg=0, ...)
{
  top <- ifelse('main' %in% names(list(...)), 2, .1)
  par(mar=c(.1,.1,top,.1))
  gx=graph.adjacency(arr[,,t])
  plot(induced.subgraph(gx,which(igraph::degree(gx)>=min.deg)), 
       edge.arrow.width=.2, edge.arrow.size=.1, 
       vertex.label.cex=.5, vertex.size=5, 
       ...)
}

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
#  Check Siena Model Converged
#   - max overall t ratio < 0.25
#   - all t-ratio < 0.1
#  @return bool
##
checkSienaConv <- function(res, t.lim=0.1,max.lim=0.25) {
  tmax <- res$tconv.max
  ts <- res$tconv
  ck.tmax <- tmax < max.lim
  ck.ts <- all(abs(ts) < t.lim)
  ck <- all(ck.tmax, ck.ts)
  cat(sprintf('\nCONVERGED:  %s\n  tconv.max (%.3f) < 0.25  %s\n  all t (max %.3f) < 0.10  %s\n\n',
              ck, tmax, ck.tmax, max(ts), ck.ts))
  return(ck)
}


###
## CrunchBase category Cosine Similarity
##
cbCatCosSim <- function(df) 
{
  all_cat_vec <- sort(unique(unlist(strsplit(df$category_list, '[|]'))))
  
  ## CATEGORY-FIRM MATRIX  [M,N] : M cateories, N firms
  firm_cats <- df$category_list
  cfm <- unname(sapply(firm_cats, function(x){
    as.integer(sapply(all_cat_vec,function(cat)grepl(cat,x)))
  }))
  ## FIRM-CATEGORY MATRIX [N,M] :  N firms, M categories
  fcm <- t(cfm)
  
  ## COMPUTE Cosine Similarity:
  ## =  u.v / |u||v|
  ## Firm-Firm COVARIANCE MATRIX [[u1.v1],[u2.v2], ...]
  X <- fcm %*% t(fcm)
  ## Firm Norms  |u|;  and |v| = |u|
  u <- apply(fcm, MARGIN = 1, FUN = function(x) sqrt(sum(x^2)) )
  ## NORM-NORM MATRIX   [ [|u1||v1|, |u1||v2|, ...], [|u2||v1|, |u2||v2|, ...], ...]
  UV <- outer(u, u, '*')
  ## COSINE SIMILARITY MATRIX equals elementwise division of 
  ##  covariance matrix by the norms matrix
  ##  = X / UV
  sim <- X / UV
  
  ## replace NULL, NA, NaN, diags with zero
  sim[is.null(sim)] <- 0
  sim[is.nan(sim)] <- 0
  sim[is.na(sim)] <- 0
  diag(sim) <- 0
  
  ## RETURN
  return(sim)
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

##===============================================
## Load Ravenpack Data
##-----------------------------------------------
csvfile <- 'rp40_action_categories_2000_2018.csv'
##convert
rvn <- read.csv(file.path(rvn_dir, csvfile), stringsAsFactors = F)


##===============================================
##  Financial controls
##-----------------------------------------------
# Size -- Ln employees, Firm age
# Financials -- quick ratio (slack)
# Performance -- lagged ROA
# Org. learning -- Acquisition experience 
# Market overlap -- lagged proportion of sales in overlapped segments (Compustat segments)

## compustat - CrunchBase names for patching missing gvkey 
## patch names with missing gvkey from crunchbase cb$co
nkpatch <- data.frame(
  company_name_unique=c('blackberry', 'broadsoft','cnova', 'compuware',
                        'csc', 'digital-river','fujitsu','geeknet',
                        'lg-display-co', 'looksmart', 'naspers', 'opera-software',
                        'siemens','sony','sybase','syniverse-technologies',
                        'telstra','velti','violin-memory','xerox'),
  conm=c('BLACKBERRY LTD','BROADSOFT INC',  'CNOVA NV', 'COMPUWARE CORP',
         'CSC HOLDINGS LLC', 'DIGITAL RIVER INC','FUJITSU LTD', 'GEEKNET INC',
         'LG DISPLAY CO LTD', 'LOOKSMART GROUP INC', 'NASPERS LTD', 'OPERA LTD -ADR',
         'SIEMENS AG', 'SONY CORP', 'SYBASE INC','SYNIVERSE HOLDINGS',
         'TELSTRA CORP LTD', 'VELTI PLC', 'VIOLIN MEMORY INC','XEROX CORP'
  ), stringsAsFactors=F)

## Compustat SEGMENTS
css<- read.csv(file.path(cs_dir,'segments.csv'), nrows=Inf, stringsAsFactors = F)
css <- css[which(css$stype=='BUSSEG'),]
css.cols <- c('gvkey','conm','sales',
  "datadate", "curcds", "isosrc", "naicsh", 
  "srcs",  "upds", "NAICSS1", "NAICSS2", "SICS1",   
  "SICS2","geotp", "snms", "soptp1","soptp2",   
  "tic", "cusip", "cik", "sic",      
  "naics",'sales_dc','sales_fn','salexg','salexg_dc')
css2 <- css[,css.cols]
css2$year <- year(ymd(css2$datadate))
css2 <- unique(css2)  ## remove duplicates from diff sourcedate (only using datadate)
View(css2[sample(1:1e5,200,F),])

## Compustat segments name gvkey 
cssnm <- unique(css2[,c('conm','gvkey','tic','sic','cusip')])

## segment totals
css2$snms2 <- sapply(css2$snms,function(x) str_to_lower(trimws(x)) )
SEGT <- ddply(css2,.(SICS1,year),summarize,.progress = 'text',
              sum=sum(sales, na.rm = T),
              sd=sd(sales, na.rm = T),
              min=min(sales, na.rm = T),
              avg=mean(sales, na.rm = T),
              max=max(sales, na.rm = T),
              hhi=hhi(sales, na.rm = T),
              n=length(which(!is.na(sales))),
              snms=paste(unique(snms),collapse = "|"),
              sic=paste(unique(c(SICS1,SICS2,sic)),collapse = "|"),
              conm=paste(unique(conm),collapse = "|"))
# # View(segt)
# css2[grepl('alphabet',css2$conm,T,T) & css2$year==2014, ]
# segt[grepl('google',segt$snms2,T,T) & segt$year==2013, ]


## Compustat annual fundamentals
csa <- read.csv(file.path(cs_dir,'fundamentals-annual.csv'), stringsAsFactors = F)
## SELECT COLUMNS FROM COMPUSTAT
csa.cols <- c('conm','conml','gvkey','datadate','fyear','indfmt','consol','popsrc','tic','cusip',
          'act', ## total assets  (ln for size proxy)
          'che', ## cash and short term investments (scale by total assets for cash holdings)
          'emp', ## employees (ln employee size proxy) 
          'ebitda', ## ebidta (scale by total assets for ROA performance proxy)
          'rect', ##  receivables total
          'lct', ## Current liabilities total
          'act',  ## current assets
          'invt'  ## inventories total
          
)
csa <- csa[,csa.cols] ## replace in memory
csa2 <- csa
## coalese quick assets (else current assets)
csa2$quick.a <- apply(csa2[,c('act','invt','che','rect')],1,function(x){
  a <- x[1] - x[2]
  b <- x[3] + x[4]
  return(ifelse(!is.na(a), a, b))
})
csa2$quick <- csa2$quick.a / csa2$lct
csa2$roa <- csa2$ebitda / csa2$act

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
  ## year names off by 1 (2006 - 2015) <-- +1
  dfl$year <- dfl$year + 1
  # dfl$name <- as.character(dfl$name)
  

  ys <- y1:y2
  nl <- list()
  bl <- list()
  for (t in 1:10) 
  {
    netpdfile <- sprintf('acq_sys_mmc_panel_RPactions_NETLIST_ego-%s_d%s_%s-%s_t%s.rds', 
                         firm_i_ego, d, y1, y2, t)
    nl[[t]] <- readRDS(file.path(result_dir, netpdfile))
  }
  names(nl) <- ys[(length(ys)-10+1):length(ys)]
  npds <- length(nl)
  
  firmnamesall <- unique(c(unlist(sapply(nl,function(x)V(x$gt)$vertex.names))))
  firmnames0 <- V(nl$`2007`$gt)$vertex.names
  adj0 <- as_adj(nl$`2007`$gt, sparse = T)
  n0 <- length(firmnames0)
  
  ## which firm made at least 1 acquisitions
  firm.acq.1 <- unique(dfl$name[which(dfl$rp_Acquisitions > 0)])
  .filtername <- unique(dfl$name[which(
    dfl$cent_deg_mmc>0 
    & dfl$ipo>0 
    #& (dfl$name %in% firm.acq.1)  ## filter by who made at least one acquisition ?
  )])
  ## filter names in vertex order with MMC rivals (is IPO firm)
  firmnamesub <- V(nl$`2007`$gt)$vertex.names[which(V(nl$`2007`$gt)$vertex.names %in% .filtername)]
  nsub <- length(firmnamesub)
  
  ##----------------------------
  ## NAMEKEY for compustat controls
  ##_---------------------------
  ## MAP names to COMPUSTAT 
  namekey <- data.frame(name=firmnamesub, stringsAsFactors = F)
  key.cols <- c('company_name','company_name_unique','cusip','tic','gvkey','conm','stock_symbol','stock_exchange_symbol')
  namekey <- merge(namekey, cb$co[,key.cols], by.x='name', by.y='company_name_unique', all.x=T, all.y=F)
  ## print missing firms to patch
  for (i in which(is.na(namekey$gvkey))) {
    idx <- grep(namekey$company_name[i],cssnm$conm,T,T)
    if (length(idx)) {
      cat(sprintf('\n\nidx=%s %s (%s)\n', paste(idx,collapse = '|'), namekey$company_name[i], namekey$name[i]))
      print(cssnm[idx,])
    }
  }
  nkpatch2 <- merge(nkpatch,cssnm[which(cssnm$conm%in%nkpatch$conm),],by='conm',all.x=T,all.y=F)
  for (i in which(is.na(namekey$gvkey))) {
    cat(sprintf('%s\n',i))
    .nkid <- which( namekey$name[i] %in% nkpatch2$company_name_unique )
    if (length(.nkid) & length(nkpatch2$cusip[.nkid])) {
      namekey$cusip[i] <- nkpatch2$cusip[.nkid]
      namekey$gvkey[i] <-  nkpatch2$gvkey[.nkid]
      namekey$tic[i] <-  nkpatch2$tic[.nkid]  # namekey$sic[i] <-  nkpatch2$sic[.nkid]
    }
  }
  #SYNIVERSE HOLDINGS INC-OLD	158473	SVR.1	7373
  # write.csv(namekey, 'acq_sys_mmc_cs_name_key_patch.csv')
  #*********************88888888888
  
  ##=============================================
  ## Select Years
  ##---------------------------------------------
  ## Aggregate DV Acquisitions by period (2 years, 3 years)
  agpdlist <- list(c('2010'),c('2011'),c('2012'),c('2013'),
                   c('2014'),c('2015'),c('2016'))
  npds <- length(agpdlist)
  ## number of years in each period to aggregate
  pdyrs <- length(agpdlist[[1]])
  ## network waves (snapshots)
  netwavepds <- as.character(c(sapply(agpdlist, function(x)as.numeric(x[1]))))
  nwaves <- length(netwavepds)  ## number of network snapshots 
  
  
  ## FIRM-lEVEL PREDICTORS
  dfla <- data.frame()
  for (t in 1:npds) {
    t_yrs <- agpdlist[[t]]
    yrs <- as.integer(t_yrs)
    cat(sprintf('t_yrs %s\n',paste(t_yrs, collapse = ',')))
    for (i in 1:length(firmnamesub)) {
      gvi <- namekey$gvkey[which(namekey$name==firmnamesub[i])]
      ##-----------------------------------------------------
      ## compustat segments for market growth
      ## compute wieghted average market growth for firm by its segments sales
      sgti <- css2[which(css2$gvkey==gvi & css2$year %in% yrs),]
      sgti1 <- css2[which(css2$gvkey==gvi & css2$year %in% (yrs-1)),] ## lag 1 year
      mkgrowth <- 1 ## default no change in markets ratio 1/1
      if (nrow(sgti) > 0) {  ## if segments data update market growth avg
        mkgvec <- c() ## get growth * proportion for each of firm's markets
        for (j in 1:nrow(sgti)) {
          sj <-  sgti$SICS1[j] ## segment
          rj <- sum( SEGT$sum[which(SEGT$SICS1==sj & SEGT$year %in% yrs)] ) # total sales revenue  in sgement
          rj1<- sum( SEGT$sum[which(SEGT$SICS1==sj & SEGT$year %in% (yrs-1))] ) # total lagged sales revenue in sgement
          gj <- ifelse( (is.na(rj)|is.na(rj1)|rj1==0),  1,  rj/rj1 ) # growth if exists
          pj <- sgti$sales[j] / sum(sgti$sales, na.rm = T) ## proportion of portfolio
          mkgvec <- c(mkgvec, pj * gj ) # weighted market growth 
        }
        mkgrowth <- sum(mkgvec, na.rm=T)
      } 
      ##compustat controls ---------------
      csti <- csa2[which(csa2$gvkey==gvi & csa2$fyear %in% yrs),]
      csti <- csti[order(csti$fyear, decreasing = F), ]
      cs.exists <- nrow(csti) > 0
      ## lagged 1 year
      csti1 <-  csa2[which(csa2$gvkey==gvi & csa2$fyear %in% (yrs-1)),]
      ## computed values from previous script  -------------------
      xti <- dfl[which(dfl$year %in% t_yrs & dfl$name==firmnamesub[i]), ]
      xti <- xti[order(xti$year, decreasing = F),]
      ## previous year
      yn1 <- as.integer(min(t_yrs)) - 1
      xtin1 <- dfl[which(dfl$year == yn1 & dfl$name==firmnamesub[i]), ]
      ## lag 2 years
      yn2 <- as.integer(min(t_yrs)) - 2
      xtin2 <- dfl[which(dfl$year == yn2 & dfl$name==firmnamesub[i]), ]
      ## lag 3 years
      yn3 <- as.integer(min(t_yrs)) - 3
      xtin3 <- dfl[which(dfl$year == yn3 & dfl$name==firmnamesub[i]), ]
      # pd tmp df to rbind
      .tmp <- data.frame(
        i=i, name=firmnamesub[i],
        pd=t, y1=t_yrs[1], y2=t_yrs[length(t_yrs)],
        type=xti$type[1], 
        # weighted average market growth
        cs_mkgrowth=mkgrowth,
        ## computstat controls
        cs_employee= ifelse(cs.exists, csti$emp[1], 0),
        cs_roa= ifelse(cs.exists, csti$roa[1], 0),
        cs_roa_1= ifelse(cs.exists, csti1$roa[1], 0),
        cs_quick= ifelse(cs.exists, csti$quick[1], 0),
        cs_quick_1= ifelse(cs.exists, csti1$quick[1], 0),
        ## State (size, position) in peirod t is ending value from range in period (t-1)
        employee_na_age= xti$employee_na_age[1],
        sales_na_0_mn= xti$sales_na_0_mn[1],
        acq_cnt_5= xti$acq_cnt_5[1],
        acq_sum_1= xti$acq_sum_1[1],
        cent_deg_mmc= xti$cent_deg_mmc[1],
        ## DV Acquisitions
        rp_Acquisitions   =mean(xti$rp_Acquisitions, na.rm=T),
        ## mean(t,t-1) of rolling 3-yr period for DVs and covars of BEHVAIOR (action counts)
        rp_Acquisitions_ma=mean(c(xti$rp_Acquisitions,
                               xtin1$rp_Acquisitions,
                               xtin2$rp_Acquisitions), na.rm=T),
        rp_NON_acquisitions=mean(xti$rp_NON_acquisitions, na.rm=T),
        #
        rp_net_restruct =sum(c(xti$rp_Acquisitions, 
                               xti$rp_New_product), na.rm=T),
        rp_net_invariant=sum(c(xti$rp_Capacity, 
                               xti$rp_Legal, 
                               xti$rp_Market_expansions,
                               xti$rp_Marketing,
                               xti$rp_Pricing,
                               xti$rp_Strategic_alliances), na.rm=T),
        ## sum(t-1,t-2,t-3) previous 3yr acquisitions == Acq Experience
        acq_exp_3=sum(c(xtin1$rp_Acquisitions,
                        xtin2$rp_Acquisitions,
                        xtin3$rp_Acquisitions), na.rm=T),
        ## competitive pressure
        wdeg_rp_Acquisitions=mean(xti$wdeg_rp_Acquisitions, na.rm=T),
        wdeg_rp_all=sum(c(xti$wdeg_rp_Acquisitions, xti$wdeg_rp_Capacity,
                          xti$wdeg_rp_Legal, xti$wdeg_rp_Market_expansions,
                          xti$wdeg_rp_Marketing,  xti$wdeg_rp_New_product,
                          xti$wdeg_rp_Pricing,  xti$wdeg_rp_Strategic_alliances), na.rm=T)
      )
      .tmp$wdeg_rp_NON_acquisitions <- .tmp$wdeg_rp_all - .tmp$wdeg_rp_Acquisitions
      #
      dfla <- rbind(dfla,.tmp)
      #
    }
  }
  
  # igraph::coreness()
  # igraph::eccentricity()
  # igraph::triangles()
  # igraph::
  
  ##------------------------------------
  ## NETWORK DV ARRAY
  ##------------------------------------
  mmcarr <- array(0, dim=c(nsub, nsub, nwaves))
  comparr <- array(0, dim=c(nsub, nsub, nwaves))
  arrSmmc <- array(0,dim=c(nsub,nwaves))
  arrSmmcSq <- array(0,dim=c(nsub,nwaves))
  arrAge <- array(0,dim=c(nsub))
  for (t in 1:nwaves) {
    wave <- netwavepds[t]
    cat(sprintf('wave %s\n',wave))
    t_nl <- which(names(nl) == wave)
    x <- nl[[t_nl]]
    ## adjmat from compnet
    compadj <- as_adj(x$gt)
    ## weighted adjmat from MMC network
    mmcadj <- as_adj(x$gmmc,attr = 'weight')
    ## MMC values ported to compnet adjmat (x * 0 --> 0; x * 1 --> x)
    compadj.mmcadj <- as.matrix(compadj * mmcadj)
    compadj.w.mmc <- compadj.mmcadj
    ## indices of MMC edges
    idx.compadj.w.mmc.gt1 <- which(compadj.w.mmc > 1)
    ## set all edges to 0 and assign only MMC edges to 1
    compadj.w.mmc[!is.na(compadj.w.mmc)] <- 0
    compadj.w.mmc[idx.compadj.w.mmc.gt1] <- 1
    #
    vids <- which(V(x$gt)$vertex.names %in% firmnamesub)
    ## Compnet subgraph for firms cohort (public acquirers)
    gtsub <- induced.subgraph(x$gt,vids)
    comparr[,,t] <- as_adj(gtsub, sparse = F)
    ## MMCnet subgraph for firms cohort
    gtsubmmc <- induced.subgraph(graph.adjacency(compadj.w.mmc),vids)
    mmcarr[,,t] <- as_adj(gtsubmmc, sparse = F)
    # ## mmc spaces weighted edges of MMC subset network (no extra spaces)
    # gtsubmmc.w <- induced.subgraph(graph.adjacency(compadj.mmcadj, weighted = T),vids)
    ##-----AGE-------------------------
    # AGE
    arrAge <- V(x$gt)$age[vids]
    ##---SMMC-CENTRALITY-----------------
    arrSmmc[,t] <- log( 1 + igraph::degree(gtsubmmc) )
    # arrSmmc[,t] <- log(1 + igraph::authority_score(gtsubmmc.w)$vector * 100)
    # arrSmmc[,t] <- log( 1 + igraph::power_centrality(gtsubmmc, exponent = 0)*100 )
    # arrSmmc[,t] <- log(  igraph::subgraph.centrality(gtsubmmc, diag = F) )
    # arrSmmc[,t] <- igraph::coreness(gtsubmmc)
    # arrSmmc[,t] <- log (1 + 100 / igraph::eccentricity(gtsubmmc) )
    # trans <- igraph::transitivity(gtsubmmc, type = 'local') * 100
    # trans[is.nan(trans) | is.na(trans)] <- 0
    # arrSmmc[,t] <- trans
    # trans <- igraph::transitivity(gtsubmmc.w, type = 'barrat')
    # trans[is.nan(trans) | is.na(trans)] <- 0
    # arrSmmc[,t] <- trans * 100
    # arrSmmc[,t] <-  log(1 + igraph::betweenness(gtsubmmc))
    # arrSmmc[,t] <-  log( 1 + igraph::eigen_centrality(gtsubmmc)$vector * 100)
    # arrSmmc[,t] <- log(1 + authority_score(gtsubmmc)$vector*100)
    # arrSmmc[,t] <-  constraint(gtsubmmc)*100
    # const <- constraint(gtsubmmc) * 100
    # const[is.nan(const) | is.na(const)] <- 0
    # arrSmmc[,t] <- const
    # eig <- eigen_centrality(gtsubmmc.w)
    # arrSmmc[,t] <- bonpow(gtsubmmc.w, exponent = 1/max(eig$value))
    arrSmmcSq[,t] <- arrSmmc[,t]^2
  }
  
  # #########################################
  ## Save / Load Workspace Image for firm_i
  ##----------------------------------------
  workspace_file <- sprintf(sprintf('sys_mmc_workspace_pre-saom_%s_%s-%s.RData',
                                    firm_i_ego,netwavepds[1],netwavepds[length(netwavepds)]))
  # save.image(file = workspace_file)
  load(file = workspace_file)
  #########################################
  
  
  ##---------------------------------------
  ## FIRM COVARIATES (Acquisitions)
  arrNetR <- array(0,dim=c(nsub,npds))
  arrNetI <- array(0,dim=c(nsub,npds))
  arrNetRw <- array(0,dim=c(nsub,npds))
  arrNetIw <- array(0,dim=c(nsub,npds))
  arrMktGro <- array(0,dim=c(nsub,npds))
  arrAcq <- array(0,dim=c(nsub,npds))
  arrDegMmc <- array(0,dim=c(nsub,npds))
  arrDegMmcSq <- array(0,dim=c(nsub,npds))
  arrEmploy <- array(0,dim=c(nsub,npds))
  arrSales <- array(0,dim=c(nsub,npds))
  arrSlack <- array(0,dim=c(nsub,npds))
  arrAcqExper <- array(0,dim=c(nsub,npds))
  arrDumCris <- array(0,dim=c(nsub,npds)) 
  arrAcqSum <- array(0,dim=c(nsub,npds)) 
  arrNonAcqAct <- array(0,dim=c(nsub,npds)) 
  arrWdegAcq <- array(0,dim=c(nsub,npds))
  arrWdegAll <- array(0,dim=c(nsub,npds))
  arrSmmc_WdegAcq <- array(0,dim=c(nsub,npds))
  arrSmmcSq_WdegAcq <- array(0,dim=c(nsub,npds))
  arrSmmc_WdegAll <- array(0,dim=c(nsub,npds))
  arrSmmcSq_WdegAll <- array(0,dim=c(nsub,npds))
  # arrW <- array(0,dim=c(nsub,npds))
  for (t in 1:npds) {
    for (i in 1:length(firmnamesub)) {
      .namei <- firmnamesub[i]
      idx <- which(dfla$pd==t & dfla$name==.namei)
      # behavior cutoffs
      cut.r <- Inf
      cut.i <- Inf
      min.r <- 1
      min.i <- 1
      logXf <- log2
      ## RESTRUCT--------
      xr <- dfla$rp_net_restruct[idx] 
      xr <- ceiling( logXf( 1 + xr ) )
      arrNetR[i,t] <-  ifelse( xr > cut.r, cut.r, ifelse(xr < min.r, min.r, xr))
      #
      xrw <- dfla$rp_net_restruct[idx] / dfla$cs_mkgrowth[idx]
      xrw <- ceiling( logXf( 1 + xrw ) )
      arrNetRw[i,t] <-  ifelse( xrw > cut.r, cut.r, ifelse(xrw < min.r, min.r, xrw))
      ## INVARIANT -------
      xi <- dfla$rp_net_invariant[idx] 
      xi <- ceiling( logXf( 1 + xi ) )
      arrNetI[i,t] <-  ifelse( xi > cut.i, cut.i, ifelse(xi < min.i, min.i, xi))
      #
      xiw <- dfla$rp_net_invariant[idx] / dfla$cs_mkgrowth[idx]
      xiw <- ceiling( logXf( 1 + xiw ) )
      arrNetIw[i,t] <-  ifelse( xiw > cut.i, cut.i, ifelse(xiw < min.i, min.i, xiw))
      ## market growth
      arrMktGro[i,t] <- dfla$cs_mkgrowth[idx]
      ## Acquisitions
      x <- dfla$rp_Acquisitions_ma[idx]  ## Acq MOVING AVERAGE
      # x <- dfla$rp_Acquisitions[idx]   ## Acq
      z <- ceiling( logXf( 1 + x ) )
      thresh <- 4
      arrAcq[i,t] <-  ifelse(z > thresh, thresh, z)
      #
      w <- dfla$cent_deg_mmc[idx]
      arrDegMmc[i,t] <- log(1+w)
      arrDegMmcSq[i,t] <- log(1+w)^2
      #
      arrWdegAcq[i,t] <- log(1 + dfla$wdeg_rp_Acquisitions[idx])
      arrWdegAll[i,t] <- log(1 + dfla$wdeg_rp_all[idx])
      #
      arrEmploy[i,t]    <- dfla$cs_employee[idx]
      arrSales[i,t]     <- dfla$cs_roa_1[idx]
      arrSlack[i,t]     <- ifelse(is.na(dfla$cs_quick[idx]), 0, dfla$cs_quick[idx])
      # arrAcqExper[i,t]  <- log( 1 + dfla$acq_cnt_5[idx] )  # log(1 + dfl$acq_cnt_5[idx])
      arrAcqExper[i,t]  <- log(1 + dfla$acq_exp_3[idx])  # log(1 + dfl$acq_cnt_5[idx])
      arrAcqSum[i,t]    <- dfla$acq_sum_1[idx]/1e+09
      arrNonAcqAct[i,t] <- log(1 + dfla$rp_NON_acquisitions[idx])
      #
      arrSmmc_WdegAcq[i,t]   <- arrWdegAcq[i,t] * arrSmmc[i,t]
      arrSmmcSq_WdegAcq[i,t] <- arrWdegAcq[i,t] * arrSmmcSq[i,t]
      arrSmmc_WdegAll[i,t]   <- arrWdegAll[i,t] * arrSmmc[i,t]
      arrSmmcSq_WdegAll[i,t] <- arrWdegAll[i,t] * arrSmmcSq[i,t]
    }
  }
  
  ## FILTER YEARS
  yridx <- 2:(length(netwavepds))#c(1,3,5,7) #1:length(netwavepds) # c(1,4,7) 
  ## FILTER FIRMS (nm=has actions; dex=has MMC relations)
  # nmsale <- which(firmnamesub %in% unique(dfla$name[!is.na(dfla$cs_roa_1)]))
  idxnm <- which(firmnamesub %in% unique(dfla$name[which(dfla$rp_net_restruct>0 | dfla$rp_net_invariant>0)]))
  idxdeg <- c(unlist(sapply(yridx, function(t)which(rowSums(mmcarr[,,t])>0))))
  
  # nmactall <- unique(dfla$name[which(dfla$rp_net_restruct>0)])
  ##
  firmsubidx <- intersect(idxnm,idxdeg) ## ## FILTER ACQUISITIONS > 0
  
  # firmsubidx <- 1:length(firmnamesub)  ## KEEP ALL
  ##
  firmnamesub2 <- firmnamesub[firmsubidx]
  nsub2 <- length(firmnamesub2)
  # nlyrs <- as.numeric(names(nl))
  # yrsubidx <- (length(nlyrs)-4):(length(nlyrs)-1)
  # yrsub <- nlyrs[yrsubidx]
  #
  mmcarr2 <- mmcarr[firmsubidx,firmsubidx, yridx ]
  comparr2 <- comparr[firmsubidx,firmsubidx, yridx ]
  #
  arrNetR2 <- arrNetR[firmsubidx, yridx ]
  arrNetI2 <- arrNetI[firmsubidx, yridx ]
  arrNetRw2 <- arrNetRw[firmsubidx, yridx ]
  arrNetIw2 <- arrNetIw[firmsubidx, yridx ]
  arrMktGro2 <- arrMktGro[firmsubidx, yridx ]
  arrAcq2 <- arrAcq[firmsubidx, yridx ]
  ##-----------------------------
  # ##*****MMC subset of compnet (ONLY CrunchBase relations filtered if has MMC)******
  arrSmmc2 <- arrSmmc[firmsubidx, yridx ]
  arrSmmcSq2 <- arrSmmcSq[firmsubidx, yridx ]
  # ## All MMC spaces count (not just CrunchBase relations) Degree Centrality
  # arrSmmc2 <- arrDegMmc[firmsubidx,  ]
  # arrSmmcSq2 <- arrDegMmcSq[firmsubidx,  ]
  ##------------------------------
  #
  arrWdegAcq2 <- arrWdegAcq[firmsubidx, yridx ]
  arrWdegAll2 <- arrWdegAll[firmsubidx, yridx ]
  arrSmmc_WdegAcq2   <- arrSmmc_WdegAcq[firmsubidx, yridx ]
  arrSmmcSq_WdegAcq2 <- arrSmmcSq_WdegAcq[firmsubidx, yridx ]
  arrSmmc_WdegAll2   <- arrSmmc_WdegAll[firmsubidx, yridx ]
  arrSmmcSq_WdegAll2 <- arrSmmcSq_WdegAll[firmsubidx, yridx ]
  #
  arrEmploy2 <- arrEmploy[firmsubidx, yridx ]
  arrSales2 <- arrSales[firmsubidx, yridx ]
  arrSlack2 <- arrSlack[firmsubidx, yridx ]
  arrAcqExper2 <- arrAcqExper[firmsubidx, yridx ]
  arrAcqSum2 <- arrAcqSum[firmsubidx, yridx ]
  arrNonAcqAct2 <- arrNonAcqAct[firmsubidx, yridx ]
  arrAge2 <- arrAge[ firmsubidx ]
  
  # par(mfrow=c(2,2), mar=c(4.5,2,.3,2))
  par(mfrow=c(1,1))
  behdf <-              within(plyr::count(c(arrNetR2[])),{beh<-'Restructuring'; MarketGrowth<-F })
  behdf <- rbind(behdf, within(plyr::count(c(arrNetI2[])),{beh<-'Invariant'; MarketGrowth<-F }))
  behdf <- rbind(behdf, within(plyr::count(c(arrNetRw2[])),{beh<-'Restructuring'; MarketGrowth<-T }))
  behdf <- rbind(behdf, within(plyr::count(c(arrNetIw2[])),{beh<-'Invariant'; MarketGrowth<-T }))
  matplot(t(arrNetRw2[c(30,45),]), type='b')
  arrPlot(mmcarr2,5,0)
  ggplot(behdf, aes(x=x, y=freq, fill=MarketGrowth)) + 
    geom_bar(stat="identity", width=.5, position = "dodge") +
    facet_wrap(.~beh) + ggtitle('Competitive Aggressiveness') + 
    xlab(NULL) + ylab('Frequency') +
    theme_bw() + theme(legend.position='bottom')
  ##----------------------------
  
  # ## DYAD FIXED COVARIATES
  arrMktOvr <- array(1, dim=c(nsub2, nsub2))
  simdf <- cb$co[which(cb$co$company_name_unique %in% firmnamesub2), c('company_name_unique','category_list')]
  arrMktOvr2 <- as.array(cbCatCosSim(simdf))
  
  ##-------------------
  ## SAOM INIT VARS
  ##------------------
  depMMC <- sienaDependent(mmcarr2, type="oneMode", nodeSet=c('FIRMS'), 
                           sparse = FALSE, allowOnly = FALSE)
  depComp <- sienaDependent(comparr2, type="oneMode", nodeSet=c('FIRMS'), 
                           sparse = FALSE, allowOnly = FALSE)
  # TREATMENT
  depAcq <- sienaDependent(arrAcq2, type="behavior", nodeSet=c('FIRMS'), 
                             sparse = FALSE, allowOnly = FALSE)
  depNetR <- sienaDependent(arrNetRw2, type="behavior", nodeSet=c('FIRMS'), 
                              sparse = FALSE, allowOnly = FALSE)
  depNetI <- sienaDependent(arrNetIw2, type="behavior", nodeSet=c('FIRMS'), 
                              sparse = FALSE, allowOnly = FALSE)
  depNonAcq <- sienaDependent(arrNonAcqAct2, type="behavior", nodeSet=c('FIRMS'), 
                           sparse = FALSE, allowOnly = FALSE)
  # #PREDICTOR
  covSmmc <- varCovar(arrSmmc2, nodeSet="FIRMS")
  covSmmcSq <- varCovar(arrSmmcSq2, nodeSet="FIRMS")
  covWdegAcq <- varCovar(arrWdegAcq2, nodeSet="FIRMS")
  covWdegAll <- varCovar(arrWdegAll2, nodeSet="FIRMS")
  # DV scaled
  covNetRw <- varCovar(arrNetRw2, nodeSet="FIRMS")
  covNetIw <- varCovar(arrNetIw2, nodeSet="FIRMS")
  covMktGro <- varCovar(arrMktGro2, nodeSet="FIRMS")
  # #CONTROLS
  covEmploy <- varCovar(arrEmploy2, nodeSet="FIRMS")
  covSales <- varCovar(arrSales2, nodeSet="FIRMS")
  covAcqExper <- varCovar(arrAcqExper2, nodeSet="FIRMS")
  covAcqSum <- varCovar(arrAcqSum2, nodeSet="FIRMS")
  covSlack <- varCovar(arrSlack2, nodeSet="FIRMS")
  # covNonAcqAct <- varCovar(arrNonAcqAct2, nodeSet="FIRMS")
  covSmmc_covWdegAcq <- varCovar(arrSmmc_WdegAcq2, nodeSet="FIRMS")
  covSmmcSq_covWdegAcq <- varCovar(arrSmmcSq_WdegAcq2, nodeSet="FIRMS")
  covSmmc_covWdegAll <- varCovar(arrSmmc_WdegAll2, nodeSet="FIRMS")
  covSmmcSq_covWdegAll <- varCovar(arrSmmcSq_WdegAll2, nodeSet="FIRMS")
  ## Constant covar
  covAge <- coCovar(arrAge2, nodeSet = "FIRMS")
  ## Constant Dyadic covar
  covCatSim <- coDyadCovar(arrMktOvr2, nodeSets = c("FIRMS","FIRMS"))
  #
  # NODES
  firms <- firmnamesub2
  nrows <- length(firms)
  FIRMS <- sienaNodeSet(nrows, 'FIRMS', firms)
  ## create system data object
  # sysEffects1 <- includeEffects(sysEffects1, inPopSqrt, name="depActiv", interaction1 = 'depMMC')
  #
  # pharmaEffects <- includeTimeDummy(pharmaEffects, outdeg, name="depActiv", interaction1 = 'depMMC', timeDummy = '5,6,7,8')

  
  ##--------------------------------------
  ##
  ##
  ##   H1ab / H2
  ##
  ##
  ##--------------------------------------
  sysDat <- sienaDataCreate(depNetR, depNetI, depMMC,  #depNonAcq,
                            covEmploy, covAcqExper, covAge, covSales,
                            covCatSim, covSlack,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  sysEff <- includeEffects(sysEff, gwesp,
                           name="depMMC")
  sysEff <- includeEffects(sysEff, outRateLog,
                           name="depMMC", type='rate')
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, X,
                           name="depMMC", interaction1 = 'covCatSim')
  # sysEff <- includeEffects(sysEff, totDist2,
  #                          name="depMMC", type='creation', interaction1 = 'covNetRw')
  sysEff <- includeEffects(sysEff, totDist2, 
                              name="depMMC", type='creation', interaction1 = 'depNetR')
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, isolate,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSlack')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depNetI", interaction1 = 'covMktGro')
  #
  sysEff <- includeEffects(sysEff, isolate,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covAcqExper')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covSlack')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depNetR", interaction1 = 'covMktGro')
  ## BEHAVIOR
  # control rivals' pressure
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetR", interaction1 = 'depMMC')
  # H1a H2a (degree effect on behavior level)
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetI", interaction1 = 'depMMC')
  # H1b H2b  (degree effect on rate of behavior change)
  sysEff <- includeEffects(sysEff, outRate,
                           name="depNetR", type='rate', interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, outRate,
                           name="depNetI", type='rate', interaction1 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-DVorig-0', 
                                 firstg = 0.05,  ## default: 0.2
                                 n2start=120,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResAMC0 <- siena07(sysMod, data=sysDat, effects=sysEff, 
                       # prevAns = sysResAMC0,
                       batch = T,   returnDeps = T, ## necessary for GOF
                       useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResAMC0)
  # 
  print(sysResAMC0); screenreg(sysResAMC0, digits = 4,single.row = T)
  gfAMC0.br = RSiena::sienaGOF(sysResAMC0, BehaviorDistribution,
                              varName="depNetR"); plot(gfAMC0.br)
  gfAMC0.bi = RSiena::sienaGOF(sysResAMC0, BehaviorDistribution,
                               varName="depNetI"); plot(gfAMC0.bi)
  gfAMC0.od = RSiena::sienaGOF(sysResAMC0, OutdegreeDistribution,
                               varName="depMMC"); plot(gfAMC0.od)
  gfAMC0.tc = RSiena::sienaGOF(sysResAMC0, TriadCensus,
                               varName="depMMC"); plot(gfAMC0.tc)

  print01report(sysResAMC0)
  siena.table(sysResAMC0, 'acq_sys_mmc_sysResAMC0_H1ab_H2_converg_DVorig.html', type = 'html', sig = T, d = 3, vertLine = T)
  siena.table(sysResAMC0, 'acq_sys_mmc_sysResAMC0_H1ab_H2_converg_DVorig.tex', type = 'tex', sig = T, d = 3, vertLine = T)
  saveRDS(list(res=sysResAMC0, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_H1ab_H2_converg_DVorig.rds')
  
  ###  
###  
###  
### 
  ###  
  ###  
  ###  
  
  
  
  ##--------------------------------------
  ##
  ##
  ##   CONTROL
  ##
  ##
  ##--------------------------------------
  sysDat <- sienaDataCreate(depNetR, depNetI, depMMC,  #depNonAcq,
                            covEmploy, covAcqExper, covAge, covSales,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  sysEff <- includeEffects(sysEff, gwesp,
                           name="depMMC")
  # sysEff <- includeEffects(sysEff, totDist2,
  #                          name="depMMC", type='creation', interaction1 = 'depNetR')
  # sysEff <- includeEffects(sysEff, altX,
  #                          name="depMMC", interaction1 = 'depNetR')
  # sysEff <- includeEffects(sysEff, altX,
  #                          name="depMMC", interaction1 = 'depNetR')
  # sysEff <- includeEffects(sysEff, totDist2,
  #                          name="depMMC", interaction1 = 'depNetI')
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  #
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covAcqExper')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  ## BEHAVIOR
  # sysEff <- includeEffects(sysEff, outdeg,
  #                          name="depNetR", interaction1 = 'depMMC')
  # sysEff <- includeEffects(sysEff, outdeg,
  #                          name="depNetI", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, avAlt,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, avAlt,
                           name="depNetR", interaction1 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-AMC-0c', 
                                 firstg = 0.05,  ## default: 0.2
                                 n2start=100,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResAMC0c <- siena07(sysMod, data=sysDat, effects=sysEff, 
                        # prevAns = sysResAMC0c,
                        batch = T,   returnDeps = T, ## necessary for GOF
                        useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResAMC0c)
  # 
  print(sysResAMC0c); screenreg(sysResAMC0c, digits = 4,single.row = T)
  gfAMC0.br = RSiena::sienaGOF(sysResAMC0c, BehaviorDistribution,
                               varName="depNetR"); plot(gfAMC0.br)
  gfAMC0.bi = RSiena::sienaGOF(sysResAMC0c, BehaviorDistribution,
                               varName="depNetI"); plot(gfAMC0.bi)
  gfAMC0.od = RSiena::sienaGOF(sysResAMC0c, OutdegreeDistribution,
                               varName="depMMC"); plot(gfAMC0.od)
  gfAMC0.tc = RSiena::sienaGOF(sysResAMC0c, TriadCensus,
                               varName="depMMC"); plot(gfAMC0.tc)
  
  
  saveRDS(list(res=sysResAMC0, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_H1ab_H2_support.rds')
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##-------------------------------------
  ## SAOM *** BEHAVIOR *** H1 
  ##-------------------------------------
  sysDat <- sienaDataCreate(depAcq, #depNonAcq,
                            covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                            covEmploy, covSales, covAcqExper,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covAcqExper')
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmc')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2t-H1b', 
                                 # behModelType = c(depAcq=2),
                                 firstg = 0.1,  ## default: 0.2
                                 n2start=80,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  
  sysResH1b <- siena07(sysMod, data=sysDat, effects=sysEff, 
                       # prevAns = sysResH1b,
                      batch = T,   returnDeps = T, ## necessary for GOF
                      useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH1b)
  # 
  summary(sysResH1b); screenreg(sysResH1b, digits = 4,single.row = T)
  gfh1b.be = RSiena::sienaGOF(sysResH1b, BehaviorDistribution,
                             varName="depAcq"); plot(gfh1b.be)
  
  ##-------------------------------------
  ## SAOM *** BEHAVIOR *** H2
  ##-------------------------------------
  sysDat <- sienaDataCreate(depAcq, #depNonAcq,
                            covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                            covEmploy, covSales, covAcqExper,
                            covWdegAcq, covSmmcSq_covWdegAcq,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covAcqExper')
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmc')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covWdegAcq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2t-H2b', 
                                 # behModelType = c(depAcq=2),
                                 firstg = 0.1,  ## default: 0.2
                                 n2start=120,    ## default: 2.52*(p+7)
                                 nsub = 5,       ## default: 4
                                 seed=135, maxlike=F)
  
  sysResH2b <- siena07(sysMod, data=sysDat, effects=sysEff, 
                       # prevAns = sysResH2b,
                       batch = T,   returnDeps = T, ## necessary for GOF
                       useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH2b)
  # 
  summary(sysResH2b); screenreg(sysResH2b, digits=4,single.row=T)
  gfh2b.be = RSiena::sienaGOF(sysResH2b, BehaviorDistribution,
                              varName="depAcq"); plot(gfh2b.be)
  
  
  ##-------------------------------------
  ## SAOM *** BEHAVIOR + STRUCT *** H2
  ##-------------------------------------
  sysDat <- sienaDataCreate(depAcq, depMMC,
                            covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                            covEmploy, covSales, covAcqExper,
                            covWdegAcq, covSmmcSq_covWdegAcq,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covAcqExper')
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmc')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covWdegAcq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2t-H2s', 
                                 # modelType = c(depMMC=2),                                  
                                 # behModelType = c(depAcq=2),
                                 firstg = 0.1,  ## default: 0.2
                                 n2start=110,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  
  sysResH2s <- siena07(sysMod, data=sysDat, effects=sysEff, 
                       # prevAns = sysResH2b,
                       batch = T,   returnDeps = T, ## necessary for GOF
                       useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH2s)
  # 
  summary(sysResH2s); screenreg(sysResH2s, digits = 4,single.row = T)
  gfh2s.be = RSiena::sienaGOF(sysResH2s, BehaviorDistribution,
                              varName="depAcq"); plot(gfh2s.be)
  
  
  
  
  ##-------------------------------------
  ## SAOM *** COEVOLUTION *** H2
  ##-------------------------------------
  sysDat <- sienaDataCreate(depAcq, depMMC,
                            covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                            covEmploy, covSales, covAcqExper,
                            covWdegAcq, 
                            covSmmc_covWdegAcq,  covSmmcSq_covWdegAcq,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  ## NETWORK CONTROL
  sysEff <- includeEffects(sysEff, density, gwesp, #density
                            name="depMMC")
  ## NETWORK SELECTION
  # sysEff <- includeEffects(sysEff, RateX, #density
  #                          name="depMMC",interaction1 = 'depAcq')
  sysEff <- includeEffects(sysEff, altX,
                            name="depMMC",interaction1 = 'depAcq')
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covAcqExper')
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmc')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covWdegAcq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmc_covWdegAcq')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  # ## BEHAVIOR INFLUENCE
  sysEff <- includeEffects(sysEff, avSimEgoX,
                           name="depAcq", interaction1 = 'covSmmcSq', interaction2 = 'depMMC')
  
  sysEff <- includeEffects(sysEff, avAlt,
                          name="depAcq", interaction1="depMMC")
  sysEff <- includeInteraction(sysEff, quad, avAlt,
                              name="depAcq", interaction1=c("","depMMC"))
  
  
  # sysEff <- includeEffects(sysEff, totSimEgoX,
  #                                name="depAcq", interaction1 = 'covSmmcSq', interaction2 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2t-H2c', 
                                 # modelType = c(depMMC=2), 
                                 # behModelType = c(depAcq=2),
                                 firstg = 0.1,  ## default: 0.2
                                 n2start=120,    ## default: 2.52*(p+7)
                                 nsub = 5,       ## default: 4
                                 seed=135, maxlike=F)
  
  sysResH2c <- siena07(sysMod, data=sysDat, effects=sysEff, #prevAns = sysResH2c,
                       #prevAns = sysResH2c,
                       batch = T,   returnDeps = T, ## necessary for GOF
                       useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH2c)
  #
  screenreg(sysResH2c,single.row = T,digits=4)
  summary(sysResH2c) 
  screenreg(sysResH2c, digits = 4,single.row = T)
  
  
  cnt = plyr::count(c(arrAcq2))
  cnt$pct = round(100 * cnt$freq / sum(cnt$freq), 1)
  print(cnt)
  gf2c.be = RSiena::sienaGOF(sysResH2c, BehaviorDistribution,
                          varName="depAcq"); plot(gf2c.be)
  gf2c.od = RSiena::sienaGOF(sysResH2c, OutdegreeDistribution,
                          varName="depMMC"); plot(gf2c.od)
  gf2c.tr = RSiena::sienaGOF(sysResH2c, TriadCensus,
                          varName="depMMC"); plot(gf2c.tr)
  #
  

  
  # ## COMINED REGRESSION TABLE
  sysResList <- list(sysResH1b,sysResH2b,sysResH2s,sysResH2c)
  screenreg(sysResList, single.row = T, digits = 3)
  
  
  ## SIMULATE
  sim_model <- sienaAlgorithmCreate( projname = "sim_model", cond = FALSE,
                                     useStdInits = FALSE, nsub = 0 , simOnly = TRUE)
  sim_ans <- siena07( sim_model, data = sysDat, effects = sysEff )
  
  
  
  
  
  
  
  
  
  
  
  ##-------------------------------------
  ##
  ##
  ## SAOM *** COEVOLUTION *** ALTERNATIVE MODEL
  ##
  ##
  ##-------------------------------------
  sysDat <- sienaDataCreate(depAcq, depMMC,
                            covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                            covEmploy, covSales, covAcqExper,
                            covWdegAcq, covSmmc_covWdegAcq, covSmmcSq_covWdegAcq,
                            covWdegAll, covSmmc_covWdegAll, covSmmcSq_covWdegAll,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  ## NETWORK CONTROL
  sysEff <- includeEffects(sysEff, density, gwesp, #density
                           name="depMMC")
  ## NETWORK SELECTION
  # sysEff <- includeEffects(sysEff, RateX, #density
  #                          name="depMMC",interaction1 = 'depAcq')
  sysEff <- includeEffects(sysEff, altX,
                           name="depMMC",interaction1 = 'depAcq')
  ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covAcqExper')
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmc')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depAcq", interaction1 = 'covSmmcSq')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depAcq", interaction1 = 'covWdegAcq')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depAcq", interaction1 = 'covSmmc_covWdegAcq')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  # sysEff <- includeEffects(sysEff, effFrom,
  #                          name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  # ## BEHAVIOR INFLUENCE
  sysEff <- includeEffects(sysEff, avSimEgoX, #degAbsContrX,
                           name="depAcq", interaction1 = 'covSmmcSq', interaction2 = 'depMMC')
  
  sysEff <- includeTimeDummy(sysEff, effFrom, timeDummy="3,4", 
                             name="depAcq", interaction1 = 'covSmmcSq')
  # 
  # sysEff <- includeEffects(sysEff, avAlt,
  #                          name="depAcq", interaction1="depMMC")
  # sysEff <- includeInteraction(sysEff, quad, avAlt,
  #                              name="depAcq", interaction1=c("","depMMC"))
  
  
  # sysEff <- includeEffects(sysEff, totSimEgoX,
  #                                name="depAcq", interaction1 = 'covSmmcSq', interaction2 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2t-H2c', 
                                 firstg = 0.1,  ## default: 0.2
                                 n2start=120,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  
  sysResH2c <- siena07(sysMod, data=sysDat, effects=sysEff, #prevAns = sysResH2c,
                       prevAns = sysResH2c,
                       batch = T,   returnDeps = T, ## necessary for GOF
                       useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH2c)
  #
  screenreg(sysResH2c,single.row = T,digits=4)
  summary(sysResH2c) 
  screenreg(sysResH2c, digits = 4,single.row = T)
  
  
  cnt = plyr::count(c(arrAcq2))
  cnt$pct = round(100 * cnt$freq / sum(cnt$freq), 1)
  print(cnt)
  gf2c.be = RSiena::sienaGOF(sysResH2c, BehaviorDistribution,
                             varName="depAcq"); plot(gf2c.be)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##-------------------------------------
  ## SAOM  *****TESTING TESTING TESTING****
  ##-------------------------------------
  sysDat2t <- sienaDataCreate(depMMC, depAcq, #depNonAcq,
                              covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                              covEmploy, covSales, covAcqExper,
                              covWdegAcq, covSmmcSq_covWdegAcq,
                              nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEffects2t <- getEffects(sysDat2t)
  effectsDocumentation(sysEffects2t)
  ## NETWORK  COMP
  # sysEffects2t <- includeEffects(sysEffects2t, density, # transTies, gwesp,
  #                                name="depComp")
  ## NETWORK  MMC
  sysEffects2t <- includeEffects(sysEffects2t, gwesp, #density
                                 name="depMMC")
  sysEffects2t <- includeEffects(sysEffects2t, egoX, #density
                                 name="depMMC", interaction1 = 'depAcq')
  # sysEffects2t <- includeEffects(sysEffects2t, sameXTransTrip, #density
  #                                name="depMMC", interaction1 = 'depAcq')
  # sysEffects2t <- includeEffects(sysEffects2t, density, # transTies, gwesp, 
  #                                name="depComp")
  # sysEffects2t <- includeEffects(sysEffects2t, egoX, 
  #                                name="depMMC", type='creation', interaction1 = 'covSmmc')
  ## CONTROL BEHVAIOR
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covEmploy')
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covSales')
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covAcqExper')
  ## BEHAVIOR
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covSmmc')
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covSmmcSq')
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covWdegAcq')
  sysEffects2t <- includeEffects(sysEffects2t, effFrom,
                                 name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  # sysEffects2t <- includeEffects(sysEffects2t, avSim, # avSimAltX, totAlt
  #                                name="depAcq", interaction1 = 'depComp')
  # sysEffects2t <- includeEffects(sysEffects2t, totSimEgoX, # avSimAltX,
  #                                name="depAcq", interaction1= 'covSmmcSq', interaction2= 'depMMC')
  # sysEffects2t <- includeEffects(sysEffects2t, totAltEgoX,
  #                                name="depAcq", interaction1 = 'covSmmcSq', interaction2 = 'depMMC')
  # sysEffects2t <- includeInteraction(sysEffects2t, totAlt,
  #                                    name="depAcq", interaction1="depMMC")
  # sysEffects2t <- includeInteraction(sysEffects2t, quad, totAlt,
  #                                   name="depAcq", interaction1=c("","depMMC"))
  ##
  sysModel2t <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2t')
  
  sysResults2t <- siena07(sysModel2t, data=sysDat2t, effects=sysEffects2t, 
                          batch = T,   returnDeps = T, ## necessary for GOF
                          useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  
  # 
  summary(sysResults2t)
  screenreg(sysResults2t, digits = 4,single.row = T)
  
  #
  sysSAOMfile <- sprintf('acq_sys_mmc_SAOM_FIT_LIST_mod2t_nonMMCnet_ego-%s_d%s_%s.rds', 
                         firm_i_ego, d, paste(netwavepds,collapse = '-'))
  saveRDS(list(results=sysResults2t, model=sysModel2,
               data=sysDat2,  effects=sysEffects2,
               FIRMS=FIRMS, depMMC=depMMC, depAcq=depMMC), 
          file = file.path(result_dir, sysSAOMfile))
  
  cnt = plyr::count(c(arrAcq2))
  cnt$pct = round(100 * cnt$freq / sum(cnt$freq), 1)
  print(cnt)
  gf2b = RSiena::sienaGOF(sysResults2t, BehaviorDistribution,
                          varName="depAcq"); plot(gf2b)
  gf2o = RSiena::sienaGOF(sysResults2t, OutdegreeDistribution,
                          varName="depMMC"); plot(gf2o)
  gf2c = RSiena::sienaGOF(sysResults2t, TriadCensus,
                          varName="depMMC"); plot(gf2c)
  #
  

 
  ##-------------------------------------
  ## SAOM MODEL 1 - H1
  ##-------------------------------------
  sysDat1 <- sienaDataCreate(depMMC, depAcq, 
                            covWdegAcq, covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                            covEmploy, covSales, covAcqExper, 
                            covAcqSum, covNonAcqAct, 
                            covSmmc_covWdegAcq, covSmmcSq_covWdegAcq,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEffects1 <- getEffects(sysDat1)
  effectsDocumentation(sysEffects1)
  ## NETWORK 
  sysEffects1 <- includeEffects(sysEffects1, density, #transTies, gwesp, 
                                name="depMMC")
  # sysEffects1 <- includeEffects(sysEffects1, absDiffX, 
  #                               name="depMMC", interaction1 = 'depAcq')
  # sysEffects1 <- includeEffects(sysEffects1, simX, 
  #                               name="depMMC", interaction1 = 'depAcq')
  ## BEHAVIOR
  # sysEffects1 <- includeEffects(sysEffects1, linear,  
  #                               name="depAcq")
  sysEffects1 <- includeEffects(sysEffects1, effFrom,
                                name="depAcq", interaction1 = 'covSmmc')
  sysEffects1 <- includeEffects(sysEffects1, effFrom,
                                name="depAcq", interaction1 = 'covSmmcSq')
  # (controls to BEHAVIOR)
  sysEffects1 <- includeEffects(sysEffects1, effFrom,
                                name="depAcq", interaction1 = 'covEmploy')
  sysEffects1 <- includeEffects(sysEffects1, effFrom,
                                name="depAcq", interaction1 = 'covSales')
  sysEffects1 <- includeEffects(sysEffects1, effFrom,
                                name="depAcq", interaction1 = 'covAcqExper')
  # sysEffects1 <- includeEffects(sysEffects1, effFrom,
  #                               name="depAcq", interaction1 = 'covAcqSum')
  # sysEffects1 <- includeEffects(sysEffects1, effFrom,
  #                               name="depAcq", interaction1 = 'covNonAcqAct')
  
  ## model
  sysModel1 <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-1')
  ## estimate
  sysResults1 <- siena07(sysModel1, data=sysDat1, effects=sysEffects1, 
                           batch = T,   returnDeps = T, ## necessary for GOF
                           useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  
  # screenreg(sysResults1, digits = 4)
  summary(sysResults1)
  # save results
  sysSAOMfile <- sprintf('acq_sys_mmc_SAOM_FIT_LIST_mod1_nonMMCnet_ego-%s_d%s_%s-%s.rds', 
                       firm_i_ego, d, yrsub[1], yrsub[length(yrsub)])
  saveRDS(list(results=sysResults1, model=sysModel1,
               data=sysDat1,  effects=sysEffects1,
               FIRMS=FIRMS, depMMC=depMMC, depAcq=depMMC), 
          file = file.path(result_dir, sysSAOMfile))
  

  gf1b = RSiena::sienaGOF(sysResults1, BehaviorDistribution,
                         varName="depAcq")
  gf1o = RSiena::sienaGOF(sysResults1, OutdegreeDistribution,
                          varName="depMMC")
  gf1c = RSiena::sienaGOF(sysResults1, TriadCensus,
                          varName="depMMC")
  ##-------------------------------------
  ## SAOM MODEL 2 - H2
  ##-------------------------------------
  sysDat2 <- sienaDataCreate(depMMC, depAcq, 
                             covSmmc, covSmmcSq,    ##dyCovTrt, dyCovTrtLinear, 
                             covEmploy, covSales, covAcqExper, 
                             covAcqSum, covNonAcqAct, covWdegAcq,
                             covSmmc_covWdegAcq,
                             covSmmcSq_covWdegAcq ,
                             nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEffects2 <- getEffects(sysDat2)
  effectsDocumentation(sysEffects2)
  # print01Report(pharmaDat, modelname="pharma-test-1")
  ## NETWORK 
  sysEffects2 <- includeEffects(sysEffects2, density, # transTies, gwesp, 
                                name="depMMC")
  # sysEffects1 <- includeEffects(sysEffects1, absDiffX, 
  #                               name="depMMC", interaction1 = 'depAcq')
  # sysEffects1 <- includeEffects(sysEffects1, simX, 
  #                               name="depMMC", interaction1 = 'depAcq')
  ## BEHAVIOR
  # sysEffects1 <- includeEffects(sysEffects1, linear,  
  #                               name="depAcq")
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covSmmc')
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covSmmcSq')
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covWdegAcq')
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covSmmc_covWdegAcq')
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covSmmcSq_covWdegAcq')
  # (controls to BEHAVIOR)
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covEmploy')
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covSales')
  sysEffects2 <- includeEffects(sysEffects2, effFrom,
                                name="depAcq", interaction1 = 'covAcqExper')
  # sysEffects1 <- includeEffects(sysEffects1, effFrom,
  #                               name="depAcq", interaction1 = 'covAcqSum')
  # sysEffects1 <- includeEffects(sysEffects1, effFrom,
  #                               name="depAcq", interaction1 = 'covNonAcqAct')
  #
  ##
  sysModel2 <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-2')
  
  sysResults2 <- siena07(sysModel2, data=sysDat2, effects=sysEffects2, 
                         batch = T,   returnDeps = T, ## necessary for GOF
                         useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  
  # screenreg(sysResults1, digits = 4)
  summary(sysResults2)
  
  sysSAOMfile <- sprintf('acq_sys_mmc_SAOM_FIT_LIST_mod2_nonMMCnet_ego-%s_d%s_%s-%s.rds', 
                         firm_i_ego, d, yrsub[1], yrsub[length(yrsub)])
  saveRDS(list(results=sysResults2, model=sysModel2,
               data=sysDat2,  effects=sysEffects2,
               FIRMS=FIRMS, depMMC=depMMC, depAcq=depMMC), 
          file = file.path(result_dir, sysSAOMfile))
  
  
  
  
  gf2b = RSiena::sienaGOF(sysResults2, BehaviorDistribution,
                          varName="depAcq")
  gf2o = RSiena::sienaGOF(sysResults2, OutdegreeDistribution,
                          varName="depMMC")
  gf2c = RSiena::sienaGOF(sysResults2, TriadCensus,
                          varName="depMMC")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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















###Compustat Segments market overlap -- Shows no overlap between GOOGLE & MICROSOFT ?
# for (t in 1:npds) {
#   cat(sprintf('\n t=%s  %s\n  i = ', t,netwavepds[t]))
#   for (i in 1:length(firmnamesub2)) {
#     if (i %% 20 == 0) cat(sprintf('%s ', i))
#     for (j in 1:length(firmnamesub2)) {
#       if (i < j) ## only do upper triangle
#         next
#       wave <- netwavepds[t]
#       waveyr <- as.integer(wave)
#       gvi <- namekey$gvkey[which(namekey$name==firmnamesub2[i])]
#       gvj <- namekey$gvkey[which(namekey$name==firmnamesub2[j])]
#       sgi <- css2[which(css2$gvkey==gvi & css2$year==waveyr), ]
#       sgj <- css2[which(css2$gvkey==gvj & css2$year==waveyr), ]
#       # sgi$SICS1[is.na(sgi$SICS1)] <- sgi$sic[1]
#       # sgj$SICS1[is.na(sgj$SICS1)] <- sgj$sic[1]
#       sgi2 <- ddply(sgi, .(SICS1),summarize,sum=sum(sales, na.rm=T))
#       sgj2 <- ddply(sgj, .(SICS1),summarize,sum=sum(sales, na.rm=T))
#       seg.all <- unique(c(sgi2$SICS1, sgj2$SICS1))
#       s.tot <- sum(c(sgi2$sum, sgj2$sum), na.rm = T)
#       s.ovr <- 0
#       for (seg in seg.all) {
#         if (seg %in% sgi2$SICS1 & seg %in% sgj2$SICS1) {
#           s.s.i <- sgi$sales[which(sgi$SICS1==seg)]
#           s.s.j <- sgj$sales[which(sgj$SICS1==seg)]
#           s.ovr <- s.ovr + sum(c(s.s.i, s.s.j), na.rm=T)
#         }
#       }
#       arrMktOvr[i,j,t] <- ifelse(s.tot==0, 1, s.ovr / s.tot )
#     }
#   }
#   ## symmetric covariates 
#   ## replace lower triangle with upper trangle
#   mat <- arrMktOvr[,,t]
#   tmat <- t(mat)
#   mat[lower.tri(mat, diag = F)] <- tmat[lower.tri(tmat, diag = F)]
#   arrMktOvr[,,t] <- mat
# }

