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
library(RColorBrewer)


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
# Plot from adjacency matrix
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
  ts <- abs(res$tconv)
  ck.tmax <- tmax < max.lim
  ck.ts <- all(ts < t.lim)
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

##
#
##
pngPlot <- function(x, filename, height = 4, width = 6.5, units = 'in', res = 300)
{
  png(filename,height = height, width = width, units = units, res =res )
  plot(x)
  dev.off()
}

##
#
##
getPvalStars <- function(p) {
  if (is.na(p)|is.nan(p))return('   ')
  if(p < 0.001)return('***')
  if(p < 0.01) return('** ')
  if(p < 0.05) return('*  ')
  return('   ')
}

##
#
##
saomResDf <- function(res, digits=3) 
{
  df <- data.frame(
    DV=res$effects$name,
    Effect=res$effects$effectName,
    Type=res$effects$type,
    Est=round(res$theta, digits = digits),
    se=round(res$se, digits = digits),
    t=round(res$theta/res$se, digits = digits),
    p=round(pt(abs(res$theta/res$se), df=Inf, lower.tail = F) * 2, digits = digits),
    stringsAsFactors = F
  )
  idx.rate <- grep('rate',df$Effect,T,T)
  df$t[idx.rate] <- NA
  df$p[idx.rate] <- NA
  return(df)
}

##
#
##
cbindDfList <- function(dfList) {
  df <- dfList[[1]]
  if (length(dfList) > 1) {
    for (i in 2:length(dfList)) {
      df <- cbind(df, dfList[[i]])
    }
  }
  return(df)
}


##
#  Create SAOM regression comparison Table
##
saomTable <- function(resList, file=NA, nameMap=NA, digits=3, drop.p.col=TRUE, drop.dv.col=TRUE)
{
  if (class(resList) != 'list') {
    resList <- list(resList)
  }
  behNames <- c()
  netNames <- c()
  for (res in resList) {
    behNames <- c(behNames, names(res$f$Data1$Behaviors))
    netNames <- c(netNames, names(res$f$Data1$nets))
  }
  dvNames <- unique(c(behNames, netNames))
  
  nameList <- list()
  for (res in resList) {
    resDvdf <- data.frame()
    for (dv in dvNames) {
      .effDf <- as.data.frame(res$effects)
      effTypeNames <- .effDf[grep(dv,res$effects$name,T,T),c('effectName','type')]
      nameList[[dv]] <- unique(rbind(resDvdf, effTypeNames))
    }
  }
  
  dfl <- list()
  for (res in resList) {
    for (dv in dvNames) {
      if (dv %in% res$effects$name) {
        .df <- saomResDf(res, digits=digits)
        dfl[[dv]][[ length(dfl[[dv]]) + 1 ]] <- .df[which(.df$DV==dv),]
      }
    }
  }
  
  mod.cols <- c('Est','se','p')
  
  tdf <- data.frame(stringsAsFactors = F)
  for (dv in names(nameList)) {
    hasRowname <- FALSE
    for (rowi in 1:nrow(nameList[[dv]])) { ## effect row
      eff <- nameList[[dv]][rowi,]
      effRow <- list()
      for (modDf in dfl[[dv]] ) { ## model dataframe in DV group
        effId <- which(modDf$Effect == eff$effectName & modDf$Type == eff$type)
        if (length(effId) > 0) {
          effRow[[length(effRow)+1]] <- modDf[effId,mod.cols]
        } else {
          .nadf <- data.frame()
          for (col in mod.cols) .nadf[1,col] <- NA
          effRow[[length(effRow)+1]] <- .nadf
        }
      }
      effRowDf <- cbind(data.frame(DV=dv,Effect=eff$effectName, Type=eff$type, stringsAsFactors = F), cbindDfList(effRow))
      if (!hasRowname) {
        effRowDf <- rbind(effRowDf, effRowDf)
        effRowDf[1,]  <- c('', sprintf('Dynamics: %s', dv), rep(NA, ncol(effRowDf)-2))
        hasRowname <- TRUE
      }
      tdf <- rbind(tdf, effRowDf)
    }
  }
  
  # move rate rows to end
  .tmp.rate.row <- tdf[1,]
  .tmp.rate.row$Effect <- 'Rate Parameters'
  rate.idx <- which(tdf$Type=='rate')
  tdf <- rbind(tdf[-rate.idx,], .tmp.rate.row, tdf[rate.idx, ])
  
  obs <- c()
  ns <- c()
  conv <- c()
  convt <- c()
  iter <- c()
  for (res in resList) {
    obs <-c(obs, res$observations)
    ns <- c(ns, attributes(res$f$Data1$nets[[1]][[1]][[1]])$nActors)
    conv <- c(conv, res$tconv.max)
    convt<- c(convt, max(abs(res$tconv)))
    iter <- c(iter, res$n)
  }  
  
  # est idx
  idx.est <- which(names(tdf)%in% 'Est')
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Time Periods'
  tdf[nrow(tdf), idx.est] <-  obs
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Num. Firms'
  tdf[nrow(tdf), idx.est] <-  ns
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Largest converg. t ratio'
  tdf[nrow(tdf), idx.est] <-  round(convt, digits = digits)
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Overall max. converg. ratio'
  tdf[nrow(tdf), idx.est] <-  round(conv, digits = digits)
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Iterations'
  tdf[nrow(tdf), idx.est] <-  iter
  
  idx.se  <- grep('^se\\.{0,1}',names(tdf),T,T)
  idx.p <- grep('^p\\.{0,}', names(tdf),T,T)
  for (i in 1:length(idx.se)) {
    sei <- idx.se[i]
    pi <- idx.p[i]
    tdf[,sei] <- apply(tdf[,c(sei,pi)],1,function(x) {
      se <- x[1]
      p <- x[2]
      spfstr <- sprintf('(%s%s.%sf)%s','%',digits+2,digits,getPvalStars(p))
      ifelse(is.na(se)|se=='NA','',sprintf(spfstr,as.numeric(se)))
    })
  }
  
  ## add Type to name (not eval or rate)
  idx.nonrate <- which(tdf$Type %in% c('endow','creation'))
  tdf$Effect[idx.nonrate] <- apply(tdf[idx.nonrate,c('Effect','Type')],1,function(x){
    sprintf('%s: %s',x[2],x[1])
  })
  tdf <- tdf[,which(names(tdf) != 'Type')]
  
  ## name mapping for effects
  if (!any(is.na(nameMap))) 
  {
    ord <- c()
    for (eff in names(nameMap)) {
      idx.nm <- which(tdf$Effect == eff)
      ord <- c(ord, idx.nm)
      if (length(idx.nm) > 0) {
        tdf$Effect[idx.nm] <- nameMap[[eff]]
      }
    }
    tdf <- rbind(tdf[ord,], tdf[-ord,])
  }
  
  if (drop.dv.col) {
    tdf <- tdf[,-1]
  }
  if (drop.p.col) {
    idx.p.col <- grep('^p[\\.\\d]{0,}',names(tdf),T,T)
    tdf <- tdf[, -idx.p.col]
  }
  
  if (!is.na(file)) {
    write.csv(tdf, file = file, na = "", row.names = F)
  }
  
  return(tdf)
  
}


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



rvc <- plyr::count(rvn$entity_name[which(rvn$action_category=='Acquisitions' & rvn$year>=2010)])
rvc <- rvc[order(rvc$freq, decreasing = T),]
View(rvc)
head(rvc,15)
# x freq
# 5957                              Moody's Corp.  278
# 4475                                  IBM Corp.  265
# 766                   Arthur J. Gallagher & Co.  251
# 6537                               Oracle Corp.  251
# 848                                   AT&T Inc.  250
# 6885                      Pinetree Capital Ltd.  206
# 2682                         Deutsche Boerse AG  203
# 2080                         Cisco Systems Inc.  199
# 5819                            Microsoft Corp.  177
# 1189                    Berkshire Hathaway Inc.  175
# 9203 Valeant Pharmaceuticals International Inc.  175
# 3757                       General Electric Co.  170
# 6845                                Pfizer Inc.  168
# 7873                                  Shire PLC  157
# 3840                       Glencore Xstrata PLC  151

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
firm_i <- 'microsoft' #'facebook' #'ibm'  # 'microsoft'  'amazon' 'apple' 'facebook' 'ibm'
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
        ## DV product
        rp_Product = sum(c(xti$rp_New_product), na.rm=T),
        ## DV Acquisitions
        rp_Acquisitions = sum(c(xti$rp_Acquisitions), na.rm=T),
        ##        
        rp_NON_acquisitions=mean(xti$rp_NON_acquisitions, na.rm=T),
        # RESTRUCTURING VS INVARIANT
        rp_net_restruct =sum(c(xti$rp_Acquisitions, 
                               xti$rp_New_product), na.rm=T),
        rp_net_invariant=sum(c(xti$rp_Capacity, 
                               xti$rp_Legal, 
                               xti$rp_Market_expansions,
                               xti$rp_Marketing,
                               xti$rp_Pricing,
                               xti$rp_Strategic_alliances), na.rm=T),
        # STRATEGIC VS INVARIANT
        rp_strat =sum(c(xti$rp_Acquisitions, 
                        xti$rp_New_product,
                        xti$rp_Capacity,
                        xti$rp_Market_expansions,
                        xti$rp_Strategic_alliances), na.rm=T),
        rp_tacti =sum(c(xti$rp_Legal, 
                        xti$rp_Marketing,
                        xti$rp_Pricing), na.rm=T),
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
  
  
  
  
  
  
  ##########################################################
  ## Save / Load Workspace Image for firm_i
  ##--------------------------------------------------------
  # workspace_file <- sprintf(sprintf('sys_mmc_workspace_pre-saom_%s_%s-%s.RData',
  #                                   firm_i_ego,netwavepds[1],netwavepds[length(netwavepds)]))
  #  # save.image(file = workspace_file)
  # load(file = workspace_file)
  ##########################################################
  
  
  
  
  ##------------------
  ## FIRM COVARIATES (Acquisitions)
  arrNetR <- array(0,dim=c(nsub,npds))
  arrNetI <- array(0,dim=c(nsub,npds))
  arrNetRw <- array(0,dim=c(nsub,npds))
  arrNetIw <- array(0,dim=c(nsub,npds))
  arrAcq <- array(0,dim=c(nsub,npds))
  #
  arrAcqw <- array(0,dim=c(nsub,npds))
  arrProd <- array(0,dim=c(nsub,npds))
  arrProdw <- array(0,dim=c(nsub,npds))
  arrStr <- array(0,dim=c(nsub,npds))
  arrTac <- array(0,dim=c(nsub,npds))
  arrStrw <- array(0,dim=c(nsub,npds))
  arrTacw <- array(0,dim=c(nsub,npds))
  #
  arrMktGro <- array(0,dim=c(nsub,npds))
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
      cut0 <-  Inf
      cut.r <- Inf
      cut.i <- Inf
      min0 <-  1
      min.r <- 1
      min.i <- 1
      logXf <- function(x)log(x, base= exp(1) )
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
      ## Acquisitions --------
      xa <- dfla$rp_Acquisitions[idx] 
      xa <- ceiling( logXf( 1 + xa ) )
      arrAcq[i,t] <-  ifelse( xa > cut0, cut0, ifelse(xa < min0, min0, xa))
      #
      xaw <- dfla$rp_Acquisitions[idx] / dfla$cs_mkgrowth[idx]
      xaw <- ceiling( logXf( 1 + xaw ) )
      arrAcqw[i,t] <-  ifelse( xaw > cut0, cut0, ifelse(xaw < min0, min0, xaw))
      ## Product  --------
      xp <- dfla$rp_Product[idx] 
      xp <- ceiling( logXf( 1 + xp ) )
      arrProd[i,t] <-  ifelse( xp > cut0, cut0, ifelse(xp < min0, min0, xp))
      #
      xpw <- dfla$rp_Product[idx] / dfla$cs_mkgrowth[idx]
      xpw <- ceiling( logXf( 1 + xpw ) )
      arrProdw[i,t] <-  ifelse( xpw > cut0, cut0, ifelse(xpw < min0, min0, xpw))
      ## STRATEGIC******************************
      xs <- dfla$rp_strat[idx] 
      xs <- ceiling( logXf( 1 + xs ) )
      arrStr[i,t] <-  ifelse( xs > cut0, cut0, ifelse(xs < min0, min0, xs))
      #
      xsw <- dfla$rp_strat[idx] / dfla$cs_mkgrowth[idx]
      xsw <- ceiling( logXf( 1 + xsw ) )
      arrStrw[i,t] <-  ifelse( xsw > cut0, cut0, ifelse(xsw < min0, min0, xsw))
      ## TACTICAL***  --------
      xt <- dfla$rp_tacti[idx] 
      xt <- ceiling( logXf( 1 + xt ) )
      arrTac[i,t] <-  ifelse( xt > cut0, cut0, ifelse(xt < min0, min0, xt))
      #
      xtw <- dfla$rp_tacti[idx] / dfla$cs_mkgrowth[idx]
      xtw <- ceiling( logXf( 1 + xtw ) )
      arrTacw[i,t] <-  ifelse( xtw > cut0, cut0, ifelse(xtw < min0, min0, xtw))
      ##-----------------------------------
      ##-----------------------------------
      ## market growth
      arrMktGro[i,t] <- dfla$cs_mkgrowth[idx]
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
  yridx <- 1:(length(netwavepds))#c(1,3,5,7) #1:length(netwavepds) # c(1,4,7) 
  ## FILTER FIRMS 
  # nmsale <- which(firmnamesub %in% unique(dfla$name[!is.na(dfla$cs_roa_1)]))
  ## FILTER FIRMS that have at least one action in RavenPack
  idxnm <- which(firmnamesub %in% unique(dfla$name[which(dfla$rp_net_restruct>0 | dfla$rp_net_invariant>0)]))
  idxdeg <- c(unlist(sapply(yridx, function(t)which(rowSums(mmcarr[,,t])>0))))
  
  # nmactall <- unique(dfla$name[which(dfla$rp_net_restruct>0)])
  ##
  firmsubidx <- intersect(idxnm,idxdeg) ## ## FILTER ACQUISITIONS > 0
  
  # firmsubidx <- 1:length(firmnamesub)  ## KEEP ALL
  ##
  ##**********************************************************************************
  # firmnamesub2 <- firmnamesub[firmsubidx]  ###  SUBSET FIRMS BY RAVENPACK
  #-----
  firmnamesub2 <- firmnamesub                ###  NOT SUBSET FIRMS BY RAVENPACK
  firmsubidx <- 1:length(firmnamesub)
  ##**********************************************************************************
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
  arrAcq2 <- arrAcq[firmsubidx, yridx ]
  arrAcqw2 <- arrAcqw[firmsubidx, yridx ]
  arrProd2 <- arrProd[firmsubidx, yridx ]
  arrProdw2 <- arrProdw[firmsubidx, yridx ]
  arrStr2 <- arrStr[firmsubidx, yridx ]
  arrStrw2 <- arrStrw[firmsubidx, yridx ]
  arrTac2 <- arrTac[firmsubidx, yridx ]
  arrTacw2 <- arrTacw[firmsubidx, yridx ]
  #
  arrMktGro2 <- arrMktGro[firmsubidx, yridx ]
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
  
  # # par(mfrow=c(2,2), mar=c(4.5,2,.3,2))
  # ####
  # idf <-  plyr::count(c(arrNetIw2[]))
  # idf$pct <-  round(100 * idf$freq / sum(idf$freq), 1)
  # print(idf)
  # rdf <- plyr::count(c(arrNetRw2[]))
  # rdf$pct <- round(100 * rdf$freq / sum(rdf$freq), 1)
  # print(rdf)
  ######
  par(mfrow=c(1,1))
  behdf <-              within(plyr::count(c(arrNetR2[])),{beh<-'Restructuring'; MarketGrowth<-F })
  behdf <- rbind(behdf, within(plyr::count(c(arrNetI2[])),{beh<-'Invariant'; MarketGrowth<-F }))
  behdf <- rbind(behdf, within(plyr::count(c(arrNetRw2[])),{beh<-'Restructuring'; MarketGrowth<-T }))
  behdf <- rbind(behdf, within(plyr::count(c(arrNetIw2[])),{beh<-'Invariant'; MarketGrowth<-T }))
  matplot(t(arrNetRw2[c(30,45),]), type='b')
  for (t in yridx) arrPlot(mmcarr2,t,0, main=sprintf('Year %s',netwavepds[t]))
  ## PLOT BEHAVIOR (all years) by Market Growth - FACET BY AGGRESSIVENESS (RESTRUCT vs INVARIANT
  ggplot(behdf, aes(x=x, y=freq, fill=MarketGrowth)) + 
    geom_bar(stat="identity", width=.6, position = "dodge") +
    facet_wrap(.~beh) + xlab('Competitive Aggressiveness') + ylab('Frequency') +
    theme_bw() + theme(legend.position='top') + scale_fill_manual(values = c('darkgrey','steelblue'))
  ## PLOT BEHAVIOR BY YEAR - FACET BY AGGRESSIVENESS (RESTRUCT vs INVARIANT)
  behyrdf <- within(plyr::count(c(arrNetRw2[,1])),{beh<-'Restructuring'; year<-netwavepds[yridx[1]] })
  for (tt in 2:length(yridx)) 
    behyrdf <- rbind(behyrdf, within(plyr::count(c(arrNetRw2[,tt])),{beh<-'Restructuring'; year<-netwavepds[yridx[tt]] }) )
  behyrdf <- rbind(behyrdf, within(plyr::count(c(arrNetIw2[,1])),{beh<-'Invariant'; year<-netwavepds[yridx[1]] }) )
  for (tt in 2:length(yridx)) 
    behyrdf <- rbind(behyrdf, within(plyr::count(c(arrNetIw2[,tt])),{beh<-'Invariant'; year<-netwavepds[yridx[tt]] }) )
  ggplot(behyrdf, aes(x=x, y=freq, fill=year)) + 
    geom_bar(stat="identity", width=.8, position = "dodge") +
    facet_wrap(.~beh) + xlab('Competitive Aggressiveness') + ylab('Frequency') +
    theme_bw() + theme(legend.position='top') + scale_fill_manual(values=brewer.pal(9,"PuBu")[(9-length(yridx)+1):9])
  ##------------------------------
 
  ## NETWORK STABILITY CHECK > 0.3 Jaccard Index
  cat('\n JACCARD INDEX: ')
  for (t in 2:dim(mmcarr2)[3]) {
    cat(sprintf('  t[%s-%s] = %.3f\n', t-1,t,netJaccard(mmcarr2[,,t],mmcarr2[,,t-1])))
  }; cat('\n')
  ##------------------------------
  
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
  depNetR <- sienaDependent(arrNetRw2, type="behavior", nodeSet=c('FIRMS'), 
                              sparse = FALSE, allowOnly = FALSE)
  depNetI <- sienaDependent(arrNetIw2, type="behavior", nodeSet=c('FIRMS'), 
                              sparse = FALSE, allowOnly = FALSE)
  depAcq <- sienaDependent(arrAcqw2, type="behavior", nodeSet=c('FIRMS'), 
                             sparse = FALSE, allowOnly = FALSE)  
  depProd <- sienaDependent(arrProdw2, type="behavior", nodeSet=c('FIRMS'), 
                           sparse = FALSE, allowOnly = FALSE)  
  depStr <- sienaDependent(arrStrw2, type="behavior", nodeSet=c('FIRMS'), 
                            sparse = FALSE, allowOnly = FALSE)  
  depTac <- sienaDependent(arrTacw2, type="behavior", nodeSet=c('FIRMS'), 
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

  
  ##=======================================================================
  ##
  ##
  ##   H1 / H2 /     FULL
  ##
  ##
  ##-----------------------------------------------------------------------
  sysDat <- sienaDataCreate(depMMC, depNetI, depNetR,  #depNonAcq, depNetR, depNetI
                            covEmploy, covAcqExper, covAge, covSales,
                            covCatSim, covSlack,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  # sysEff <- includeEffects(sysEff, RateX,
  #                          name="depMMC", type='rate', interaction1 = 'depNetR')
  #
  sysEff <- includeEffects(sysEff, outInAss, nbrDist2, transTies, #gwesp, #inPop, #inPopSqrt, 
                           name="depMMC")
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, X,
                           name="depMMC", interaction1 = 'covCatSim')
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'depNetR')
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='eval', include = TRUE)
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='endow')
   ## CONTROL BEHVAIOR -----------------------------------------
  sysEff <- includeEffects(sysEff, isolate,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSlack')
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
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, behDenseTriads,# behDenseTriads,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, behDenseTriads,# behDenseTriads,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname=sprintf('sys-mmc-acq-proj-%s-converg-VDcut-behDenseTriad_restruct_1-Inf_H1ab_H2-0_noGWESP',firm_i_ego), 
                                 firstg = 0.12, #0.07,  ## default: 0.2
                                 n2start=240,    ## default: 2.52*(p+7)
                                 nsub = 5,       ## default: 4
                                 seed=133, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResAMC0 <- siena07(sysMod, data=sysDat, effects=sysEff, 
                       #prevAns = sysResAMC0,
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
  
  .prefix <- 'SAOM_1-Inf_H1ab_H2_noGWESP_GOF'
  png(sprintf('%s_%s.png',.prefix,'behavior_restruct'), height = 4, width = 6.5, units = 'in', res = 300)
    plot(gfAMC0.br); dev.off();
  png(sprintf('%s_%s.png',.prefix,'behavior_invariant'), height = 4, width = 6.5, units = 'in', res = 300)
    plot(gfAMC0.bi); dev.off();
  png(sprintf('%s_%s.png',.prefix,'net_degree'), height = 4, width = 6.5, units = 'in', res = 300)
    plot(gfAMC0.od); dev.off();
  png(sprintf('%s_%s.png',.prefix,'net_triad_census'), height = 4, width = 6.5, units = 'in', res = 300)
    plot(gfAMC0.tc); dev.off();
  saveRDS(list(res=sysResAMC0, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_converg_DVcut_behDenseTriad_restruct_1-Inf_H1ab_H2_noGWESP.rds')
  
  x <- readRDS('acq_sys_mmc_sysResAMC0_converg_DVcut_behDenseTriad_restruct_1-Inf_H1ab_H2_noGWESP.rds')
    
  # ans <- sysResAMC0
  ans <- x$res
  
  mat <- matrix(0,1,ans$pp)
  id <- 8:9
  mat[1,id] <- 1
  
  Wald.RSiena(mat, ans)
  
  th1 <- ans$theta[id][1]
  th2 <- ans$theta[id][2]
  se12 <- sqrt(mat %*% ans$covtheta %*% t(mat))
  cat(sprintf('linear combination:\n %.3f + %.3f = %.3f \n / %.3f\n = %.3f',
              th1, th2, th1+th2, se12, (th1+th2)/se12))
  
  sum(ans$theta[id]) / sqrt(mat %*% ans$covtheta %*% t(mat))
  
  
  
  ## GET MMC CHAINS
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-converg-VDcut-behDenseTriad_restruct_1-Inf_H1ab_H2-0_noGWESP_CHAINS')
  sysResAMC0Chains <- siena07(sysMod, data=sysDat, effects=sysEff, 
                        batch = T,   returnDeps = T, returnChains=T, ## necessary for GOF
                        useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  res$f$basicEffects$depMMC
  z = sysResAMC0Chains

  
  ###   sysResAMC0 = readRDS('acq_sys_mmc_sysResAMC0_converg_DVcut_behDenseTriad_restruct_1-6_H1ab_H2ab_H3ab.rds')
###  
###  
### 
  ###  
  ###  
  ###  
  
  
  ##=======================================================================
  ##
  ##
  ##   H1  only   Dense Triads
  ##
  ##
  ##-----------------------------------------------------------------------
  sysDat <- sienaDataCreate(depMMC, depNetI, depNetR,  #depNonAcq, depNetR, depNetI
                            covEmploy, covAcqExper, covAge, covSales,
                            covCatSim, covSlack,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  sysEff <- includeEffects(sysEff, outInAss, nbrDist2, transTies, #gwesp, #inPop, #inPopSqrt, 
                           name="depMMC")
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, X,
                           name="depMMC", interaction1 = 'covCatSim')
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'depNetR')
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='eval', include = FALSE)#***
  ## CONTROL BEHVAIOR -----------------------------------------
  sysEff <- includeEffects(sysEff, isolate,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSlack')
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
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, behDenseTriads,# behDenseTriads,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, behDenseTriads,# behDenseTriads,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-converg-VDcut-H1ab_noGWESP', 
                                 firstg = 0.13,  ## default: 0.2
                                 n2start=220,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=133, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResH1ab <- siena07(sysMod, data=sysDat, effects=sysEff, 
                        # prevAns = sysResAMC0,
                        batch = T,   returnDeps = T, ## necessary for GOF
                        useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH1ab)
  # 
  print(sysResH1ab); screenreg(sysResH1ab, digits = 4,single.row = T)
  gfH1.br = RSiena::sienaGOF(sysResH1ab, BehaviorDistribution,
                               varName="depNetR"); plot(gfH1.br)
  gfH1.bi = RSiena::sienaGOF(sysResH1ab, BehaviorDistribution,
                               varName="depNetI"); plot(gfH1.bi)
  gfH1.od = RSiena::sienaGOF(sysResH1ab, OutdegreeDistribution,
                               varName="depMMC"); plot(gfH1.od)
  gfH1.tc = RSiena::sienaGOF(sysResH1ab, TriadCensus,
                               varName="depMMC"); plot(gfH1.tc)
  
  saveRDS(list(res=sysResH1ab, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_converg_DVcut_H1ab_noGWESP.rds')
  
  
  
  ##=======================================================================
  ##
  ##
  ##   H2 only   PURPOSEFUL CREATION OF MMC DENSE TRIADS
  ##
  ##
  ##-----------------------------------------------------------------------
  sysDat <- sienaDataCreate(depMMC, depNetI, depNetR,  #depNonAcq, depNetR, depNetI
                            covEmploy, covAcqExper, covAge, covSales,
                            covCatSim, covSlack,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  sysEff <- includeEffects(sysEff, outInAss, nbrDist2, transTies, #gwesp, #inPop, #inPopSqrt, 
                           name="depMMC")
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, X,
                           name="depMMC", interaction1 = 'covCatSim')
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'depNetR')
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='eval', include = TRUE)
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='endow')
  ## CONTROL BEHVAIOR -----------------------------------------
  sysEff <- includeEffects(sysEff, isolate,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSlack')
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
  ## BEHAVIOR
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-converg-VDcut-H2_noGWESP', 
                                 firstg = 0.1,  ## default: 0.2
                                 n2start=280,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=137, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResH2 <- siena07(sysMod, data=sysDat, effects=sysEff, 
                        # prevAns = sysResAMC0,
                        batch = T,   returnDeps = T, ## necessary for GOF
                        useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResH2)
  # 
  print(sysResH2); screenreg(sysResH2, digits = 4,single.row = T)
  gfH2.br = RSiena::sienaGOF(sysResH2, BehaviorDistribution,
                             varName="depNetR"); plot(gfH2.br)
  gfH2.bi = RSiena::sienaGOF(sysResH2, BehaviorDistribution,
                             varName="depNetI"); plot(gfH2.bi)
  gfH2.od = RSiena::sienaGOF(sysResH2, OutdegreeDistribution,
                             varName="depMMC"); plot(gfH2.od)
  gfH2.tc = RSiena::sienaGOF(sysResH2, TriadCensus,
                             varName="depMMC"); plot(gfH2.tc)
  
  saveRDS(list(res=sysResH2, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_converg_DVcut_H2_noGWESP.rds')
  
  
  
  ##=======================================================================
  ##
  ##
  ##   CONTROL  COEVOLUTION CONTROs
  ##
  ##
  ##-----------------------------------------------------------------------
  sysDat <- sienaDataCreate(depNetR, depNetI, depMMC,  #depNonAcq,
                            # covNetRw, covMktGro,
                            covEmploy, covAcqExper, covAge, covSales,
                            covCatSim, covSlack,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  sysEff <- includeEffects(sysEff, outInAss, nbrDist2, transTies, #gwesp, #inPop, #inPopSqrt, 
                           name="depMMC")
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, X,
                           name="depMMC", interaction1 = 'covCatSim')
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'depNetR')
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='eval', include = FALSE)
    ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, isolate,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSlack')
  #
  ## CONTROL BEHVAIOR
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
  #
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetI", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, totAlt,
                           name="depNetR", interaction1 = 'depMMC')
  #
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetR", interaction1 = 'depMMC')
  sysEff <- includeEffects(sysEff, outdeg,
                           name="depNetI", interaction1 = 'depMMC')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-converg-VDcut-ctrl2_noGWESP', 
                                 firstg = 0.13,  ## default: 0.2
                                 n2start=230,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResAMC0ctrl2 <- siena07(sysMod, data=sysDat, effects=sysEff, 
                            # prevAns = sysResAMC0,
                            batch = T,   returnDeps = T, ## necessary for GOF
                            useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResAMC0ctrl2)
  # 
  print(sysResAMC0ctrl2); screenreg(sysResAMC0ctrl2, digits = 4,single.row = T)
  gfAMC0ctrl2.br = RSiena::sienaGOF(sysResAMC0ctrl2, BehaviorDistribution,
                                   varName="depNetR"); plot(gfAMC0ctrl2.br)
  gfAMC0ctrl2.bi = RSiena::sienaGOF(sysResAMC0ctrl2, BehaviorDistribution,
                                   varName="depNetI"); plot(gfAMC0ctrl2.bi)
  gfAMC0ctrl2.od = RSiena::sienaGOF(sysResAMC0ctrl2, OutdegreeDistribution,
                                   varName="depMMC"); plot(gfAMC0ctrl2.od)
  gfAMC0ctrl2.tc = RSiena::sienaGOF(sysResAMC0ctrl2, TriadCensus,
                                   varName="depMMC"); plot(gfAMC0ctrl2.tc)
  
  saveRDS(list(res=sysResAMC0ctrl2, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_converg_DVcut_ctrl2_noGWESP.rds')
  
  # meta0 <- siena08(sysResAMC0,
  #                  sysResAMC0ctrl, 
  #                  projname = 'sysMeta')
  # 
  
  
  
  
  
  ##=======================================================================
  ##
  ##
  ##   CONTROL
  ##
  ##
  ##-----------------------------------------------------------------------
  sysDat <- sienaDataCreate(depNetR, depNetI, depMMC,  #depNonAcq,
                            # covNetRw, covMktGro,
                            covEmploy, covAcqExper, covAge, covSales,
                            covCatSim, covSlack,
                            nodeSets=list(FIRMS) ) #,smoke1, alcohol
  ##
  sysEff <- getEffects(sysDat)
  # effectsDocumentation(sysEff)
  sysEff <- includeEffects(sysEff, outInAss, nbrDist2, transTies, #gwesp, #inPop, #inPopSqrt, 
                           name="depMMC")
  sysEff <- includeEffects(sysEff, egoX,
                           name="depMMC", interaction1 = 'covAge')
  sysEff <- includeEffects(sysEff, X,
                           name="depMMC", interaction1 = 'covCatSim')
  sysEff <- includeEffects(sysEff, transTriads,
                           name="depMMC", type='eval', include = FALSE)
    ## CONTROL BEHVAIOR
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetI", interaction1 = 'covSlack')
  #
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covAcqExper')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covEmploy')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covSales')
  sysEff <- includeEffects(sysEff, effFrom,
                           name="depNetR", interaction1 = 'covSlack')
  #
  sysMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-converg-VDcut-ctrl_noGWESP',
                                 firstg = 0.13,  ## default: 0.2
                                 n2start=180,    ## default: 2.52*(p+7)
                                 nsub = 4,       ## default: 4
                                 seed=135, maxlike=F)
  # ##***  save.image('acq_sys_mmc_SAOM_AMC.rda')  ##***
  sysResAMC0ctrl <- siena07(sysMod, data=sysDat, effects=sysEff,
                            # prevAns = sysResAMC0,
                            batch = T,   returnDeps = T, ## necessary for GOF
                            useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  conv <- checkSienaConv(sysResAMC0ctrl)
  #
  print(sysResAMC0ctrl); screenreg(sysResAMC0ctrl, digits = 4,single.row = T)
  gfAMC0ctrl.br = RSiena::sienaGOF(sysResAMC0ctrl, BehaviorDistribution,
                                   varName="depNetR"); plot(gfAMC0ctrl.br)
  gfAMC0ctrl.bi = RSiena::sienaGOF(sysResAMC0ctrl, BehaviorDistribution,
                                   varName="depNetI"); plot(gfAMC0ctrl.bi)
  gfAMC0ctrl.od = RSiena::sienaGOF(sysResAMC0ctrl, OutdegreeDistribution,
                                   varName="depMMC"); plot(gfAMC0ctrl.od)
  gfAMC0ctrl.tc = RSiena::sienaGOF(sysResAMC0ctrl, TriadCensus,
                                   varName="depMMC"); plot(gfAMC0ctrl.tc)

  saveRDS(list(res=sysResAMC0ctrl, mod=sysMod, dat=sysDat, eff=sysEff),
          file = 'acq_sys_mmc_sysResAMC0_converg_DVcut_ctrl_noGWESP.rds')





  
  
  
  
  
  ##================================================================
  ##
  ## print SAOM Regression table compariosn 
  ##
  ##---------------------------------------------------------------- 
  sysNameMap <- list(
    `Dynamics: depNetI`='Invariant Aggressiveness Dynamics',
    #
    `depNetI linear shape`='Linear Shape',
    `depNetI quadratic shape`='Quadratic Shape',
    `depNetI: effect from covEmploy`='Employees (Thousands)',
    `depNetI: effect from covSales`='ROA Lag',
    `depNetI: effect from covSlack`='Quick Ratio',
    `depNetI isolate`='Network Isolate Effect',
    `depNetI degree`='Competitor Count',
    `depNetI total alter`='Competitors\' Invariant Aggressiveness',
    `depNetI dense triads`='System MMC',
    #---------------------------------------------
    `Dynamics: depNetR`='Restructuring Aggressiveness Dynamics',
    #
    `depNetR linear shape`='Linear Shape',
    `depNetR quadratic shape`='Quadratic Shape',
    `depNetR: effect from covEmploy`='Employees (Thousands)',
    `depNetR: effect from covAcqExper`='Acquisition Experience',
    `depNetR: effect from covSales`='ROA Lag',
    `depNetR: effect from covSlack`='Quick Ratio',
    `depNetR isolate`='Network Isolate Effect',
    `depNetR degree`='Competitor Count',
    `depNetR total alter`='Competitors\' Restructuring Aggressiveness',
    `depNetR dense triads`='System MMC',
    #---------------------------------------------
    `Dynamics: depMMC`='Network Dynamics',
    #
    `covAge ego`='Firm Age',
    `covCatSim`='Category Similarity',
    `degree (density)`='Network Density',
    `degree^(1/2) assortativity`='Competitive Assortativity',
    `transitive ties`='Transitive Ties (direct + indirect MMC rivals)',
    `number of actor pairs at dist 2`='Indirect Competitors',
    `GWESP (69)`='GWESC',
    `sqrt degree of alter`='Competitor\'s Sqrt. of Competitors',
    `depNetR ego`='Focal Firm\'s Restructuring Aggressiveness',
    `depNetR tot alter at dist 2 (1)`='Indirect Competitors\' Pressure',
    `transitive triads`='System MMC Creation',
    `creation: transitive triads`='System MMC Creation',
    `endow: transitive triads`='System MMC Maintenance',
    #---------------------------------------------
    `Rate Parameters`='Rate Parameters',
    #
    `rate depNetI (period 1)`='Invariant Aggr. Rate 2010-2011',
    `rate depNetI (period 2)`='Invariant Aggr. Rate 2011-2012',
    `rate depNetI (period 3)`='Invariant Aggr. Rate 2012-2013',
    `rate depNetI (period 4)`='Invariant Aggr. Rate 2013-2014',
    `rate depNetI (period 5)`='Invariant Aggr. Rate 2014-2015',
    `rate depNetI (period 6)`='Invariant Aggr. Rate 2015-2016',
    #
    `rate depNetR (period 1)`='Restructuring Aggr. Rate 2010-2011',
    `rate depNetR (period 2)`='Restructuring Aggr. Rate 2011-2012',
    `rate depNetR (period 3)`='Restructuring Aggr. Rate 2012-2013',
    `rate depNetR (period 4)`='Restructuring Aggr. Rate 2013-2014',
    `rate depNetR (period 5)`='Restructuring Aggr. Rate 2014-2015',
    `rate depNetR (period 6)`='Restructuring Aggr. Rate 2015-2016',
    #
    `constant depMMC rate (period 1)`='Network Rate 2010-2011',
    `constant depMMC rate (period 2)`='Network Rate 2011-2012',
    `constant depMMC rate (period 3)`='Network Rate 2012-2013',
    `constant depMMC rate (period 4)`='Network Rate 2013-2014',
    `constant depMMC rate (period 5)`='Network Rate 2014-2015',
    `constant depMMC rate (period 6)`='Network Rate 2015-2016'
  )
  
  sysResList <- list(sysResAMC0ctrl, sysResAMC0ctrl2, 
                     sysResH1ab, sysResH2, sysResAMC0)
  saveRDS(sysResList, file="sys_MMC_SAOM_res_list_microsoft_DVcut_1-Inf_2010-2016_noGWESP.rds")
  tab <- saomTable(sysResList, digits = 3, nameMap = sysNameMap,
            file="sys_MMC_SAOM_res_table_microsoft_DVcut_1-Inf_2010-2016_noGWESP.csv")
  
  ## multiparameter Wald tests
  mc0 <- c('depNetR isolate', 'depNetR degree', 'depNetR total alter',
           'depNetI isolate', 'depNetI degree', 'depNetI total alter', 'depNetR ego' )
  mc1 <- c('depNetR dense triads', 'depNetI dense triads')
  # mc2 <- c('depNetR tot alter at dist 2')
  mc2 <- c()
  # TEST ALL COEVOLUTION
  ## check param IDs
  sort(sapply(mc0,function(x)grep(x,sysResAMC0ctrl2$effects$effectName,T,T)))
  sort(sapply(unique(c(mc0,mc1)),function(x)grep(x,sysResH1ab$effects$effectName,T,T)))
  sort(sapply(unique(c(mc0,mc2)),function(x)grep(x,sysResH2$effects$effectName,T,T)))
  sort(sapply(unique(c(mc0,mc1,mc2)),function(x)grep(x,sysResAMC0$effects$effectName,T,T)))
  # TEST ALL COEVOLUTION
  mtf <- ldply(list(
   Multipar.RSiena(sysResAMC0ctrl2, 13,22,23,24,37,38,39),
   Multipar.RSiena(sysResH1ab, 13,22,23,24,25,37,38,39,40),
   Multipar.RSiena(sysResH2, 15,24,25,26,38,39,40),
   Multipar.RSiena(sysResAMC0, 15,24,25,26,27,39,40,41,42) 
  ), data.frame)
  mtf$pvalue <- as.numeric( round(mtf$pvalue, 4) )
  print(mtf)

  # Multipar.RSiena(sysResAMC0ctrl2, 21,22,23,36,37,38),
  # Multipar.RSiena(sysResH1ab, 21,22,23,24,36,37,38,39),
  # Multipar.RSiena(sysResH2, 22,23,24,36,37,38),
  # Multipar.RSiena(sysResAMC0, 22,23,24,25,37,38,39,40) 
  
  ## HYPOTHESIZED ONLY
  sort(sapply(unique(c(mc1)),function(x)grep(x,sysResH1ab$effects$effectName,T,T)))
  # sort(sapply(unique(c(mc2)),function(x)grep(x,sysResH2$effects$effectName,T,T)))
  sort(sapply(unique(c(mc1,mc2)),function(x)grep(x,sysResAMC0$effects$effectName,T,T)))
  ##
  mth <- ldply(list(
    Multipar.RSiena(sysResH1ab, 25,40),
    # Multipar.RSiena(sysResH2, 11),
    Multipar.RSiena(sysResAMC0, 27,42)
  ), data.frame)
  mth$pvalue <- as.numeric( round(mth$pvalue, 4) )
  print(mth)

  
  
  
  
  
  
  
  
  
  
  
  ##------------------------
  ## LOAD data from cold start
  ##------------------------
  sysResList <- readRDS("sys_MMC_SAOM_res_list_microsoft_DVcut_1-Inf_2010-2016_noGWESP.rds")
  sysResAMC0ctrl <- sysResList[[1]]
  sysResAMC0ctrl2 <- sysResList[[2]]
  sysResH1ab <- sysResList[[3]]
  sysResH2 <- sysResList[[4]]
  sysResAMC0 <- sysResList[[5]]
  
  ## descriptives 
  res <- sysResAMC0
  beh <- 'depNetR'
  getBeh <- function(res, beh)
  {
    return(res$f$Data1$behavs[[beh]]$beh[,])
  }
  
  describeBeh <- function(res, beh)
  {
    dv <- getBeh(res, beh)
    return(data.frame(
      nodes=nrow(dv),
      periods=ncol(dv),
      nobs=nrow(dv)*ncol(dv),
      mean=mean(c(dv)),
      sd=sd(c(dv)),
      min=min(c(dv)),
      median=median(c(dv)),
      max=max(c(dv))))
  }
  
  describeBeh(res, 'depNetR')
  describeBeh(res, 'depNetI')
  
  
  library(sna)
  
  
  ## GET DESCRIPTIVES
  mmcarr2  
  
  smmc <- sapply(1:dim(mmcarr2)[3], function(t){
      kcyc <- sna::kcycle.census(mmcarr2[,,t], maxlen = 3, 
                     mode = 'graph')  #cycle.comembership = 'bylength'
      return(kcyc$cycle.count[2,-1])  ## skip "agg" first column
  })

  
  
  ## Firm networks variables
  deg <- sapply(1:dim(mmcarr2)[3],function(t)rowSums(mmcarr2[,,t]))
  mean(c(deg))
  sd(c(deg))
  
  ## competitors' INVARIANT aggressiveness  
  invBeh <- getBeh(sysResAMC0, 'depNetI')
  invAggr <- sapply(1:dim(mmcarr2)[3], function(t){
    (mmcarr2[,,t] %*% invBeh[,t]) ## diagonals == 0; so only competitors' aggressiveness is summed
  })
  
  ## competitors' RESTRUCTURING aggressiveness  
  resBeh <- getBeh(sysResAMC0, 'depNetR')
  resAggr <- sapply(1:dim(mmcarr2)[3], function(t){
    (mmcarr2[,,t] %*% resBeh[,t]) ## diagonals == 0; so only competitors' aggressiveness is summed
  })
  
  ## Firm variale - Covariates
  arrEmploy2
  arrSales2
  arrAcqExper2
  arrSlack2
  # arrAcqSum2

  ## Firm static
  arrAge2
  ## dyad static
  arrMktOvr2
  
  
  # ## assortativity
  # assort <- sapply(1:dim(mmcarr2)[3], function(t){
  #   outer(sqrt(deg[,t]),sqrt(deg[,t]), '*') ## diagonals == 0; so only competitors' aggressiveness is summed
  # })
  
  
  
  vdf <- data.frame(
    smmc = c(smmc),
    deg = c(deg),
    invAggr = c(invAggr), 
    resAggr = c(resAggr),
    invBeh = c(invBeh),
    resBeh = c(resBeh),
    employ = c(arrEmploy2),
    sales = c(arrSales2),
    exper = c(arrAcqExper2),
    slack = c(arrSlack2),
    age = c(arrAge2)
  )
  vdf[is.na(vdf)] <- 0
  
  edf <- data.frame(
    mmc = c(mmcarr2),
    catsim = c(arrMktOvr2)
  )
  edf[is.na(edf)] <- 0
  
  vdf
  edf[ row(edf) != col(edf) ]
  
  library(psych)
  cor.vdf <- cor(vdf)
  descr.vdf <- describe(vdf)
  
  write.csv(cor.vdf, file = 'smmc_DESCRIPT_vertex_COR.csv', row.names = T)
  write.csv(descr.vdf, file = 'smmc_DESCRIPT_vertex_SUMMARY.csv', row.names = F)
  
  
  critical.r <- function( n, alpha = .05 ) {
    df <- n - 2
    critical.t <- qt( alpha/2, df, lower.tail = F )
    critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
    return( critical.r )
  }
  critical.r( 518 )
  
  


  Firm Age
  Category Similarity
  Network Density
  Indirect Competitors
  Competitive Assortativity
  Transitive MMC Rivals (direct + indirect MMC)
  Focal Firms Restructuring Aggressiveness
System MMC Creation
System MMC Maintenance-Creation Difference

  
  
  
  
  
  
  
  ##=======================================================
  ##
  ##
  ##    PLOT  NETWORKS
  ##
  ##
  ##-------------------------------------------------------
  firmnamesub2
  yr1 <- netwavepds[yridx[1]]
  gt1 <- nl[[ yr1 ]]$gt
  yr2 <- netwavepds[yridx[length(yridx)]]
  gt2 <- nl[[ yr2 ]]$gt
  
  
  
  
  ##
  plotCohortNetColorMmc <- function(gs, firms, 
                                    filename = NA,
                                    fileprefix = 'sys_mmc_cohort_net',
                                    filetime = TRUE,
                                    min.degree = 1,
                                    coord=NA, seed=1111,  ...) 
  {
    gs <- igraph::induced.subgraph(gs, vids = which(igraph::degree(gs)>=min.degree) )
    #
    idx.firms <- which(V(gs)$vertex.names %in% firms)
    # ## color ------------------------------
    # V(gs)$color <- rgb(.1,.1,.9,.05)  ## light translucent blue
    # V(gs)$color[idx.firms] <- 'red'
    # V(gs)$color.frame <- rgb(.1,.1,.9,.45) 
    # V(gs)$color.frame[idx.firms] <- 'red'
    # ## color CLUSTER ----------------------
    V(gs)$cluster <- igraph::multilevel.community(gs)$membership
    ncl <- length(unique(V(gs)$cluster))
    rbvertcolors <- rainbow(ncl, s = 1, v = 1, alpha = .35)
    rbframecolors <- rainbow(ncl, s = 1, v = 1, alpha = .45)
    V(gs)$color <- rbvertcolors[ V(gs)$cluster  ]
    V(gs)$color[idx.firms] <- 'black'
    V(gs)$color.frame <- rbframecolors[ V(gs)$cluster  ]
    V(gs)$color.frame[idx.firms] <- 'black'
    
    V(gs)$size <- 2
    V(gs)$size[idx.firms] <- 3
    ##label
    vertex.label <- ''
    ##
    file.parts <- c(fileprefix, filename)
    if (filetime)
      file.parts <- c(file.parts, as.integer(Sys.time()))
    file.full <- sprintf('%s.png', paste(file.parts[!is.na(file.parts)], collapse = '_'))
    
    if (all(is.na(coord))) {
      set.seed(seed)
      minC <- rep(-Inf, vcount(gs))
      maxC <- rep(Inf, vcount(gs))
      minC[1] <- maxC[1] <- 0
      coord <- layout_with_fr(gs, minx=minC, maxx=maxC,
                              miny=minC, maxy=maxC,
                              grid='nogrid') ## defaults to grid layout for >= 1000 nodes
    }
    # coord <- layout_with_kk(gs)
    
    png(filename = file.full, width = 8, height = 8, res = 400, units = 'in')
      par(mar=c(.1,.1,.1,.1))
      set.seed(seed)
      plot(gs, layout=coord
          , vertex.size=V(gs)$size
          , vertex.color=V(gs)$color
          , vertex.frame.color = V(gs)$color.frame
          , vertex.label=vertex.label
          , vertex.label.cex=.01
          , vertex.label.color=NA
          , vertex.label.font = 2
          , vertex.label.family = 'sans'
          , vertex.shape = 'circle'
          , edge.arrow.width = .1
          , edge.arrow.size = .1
          , edge.width = .3
          , edge.color = 'darkgray'
          , rescale = TRUE
          , ...)
      legend('topright', legend=filename, bty = "n")
    dev.off()
    return(coord)
  }
  
  
  
  # coord1 <- plotCohortNetColorMmc(gt1, firmnamesub2, yr1)
  coord1 <- plotCohortNetColorMmc(gt1, firmnamesub2, yr1, coord=coord1)
  
  # coord2 <- plotCohortNetColorMmc(gt2, firmnamesub2, yr2)
  coord2 <- plotCohortNetColorMmc(gt2, firmnamesub2, yr2, coord=coord2)
  
  

  
  
  
  plotMmcNetColorBeh <- function(mat, beh,
                                  filename = NA,
                                  fileprefix = 'sys_mmc_net_beh',
                                  legend.title = 'Aggressiv.',
                                  filetime = TRUE,
                                  min.degree = 1,
                                  coord=NA, seed=111,  ...) 
  {
    gs <- igraph::graph.adjacency(mat, diag = F)
    #
    # ## color ------------------------------
    # V(gs)$color <- rgb(.1,.1,.9,.05)  ## light translucent blue
    # V(gs)$color[idx.firms] <- 'red'
    # V(gs)$color.frame <- rgb(.1,.1,.9,.45) 
    # V(gs)$color.frame[idx.firms] <- 'red'
    # ## color CLUSTER ----------------------
    ncl <- length(unique(beh))
    rbvertcolors <- heat.colors(ncl, alpha = .75, rev = T)
    V(gs)$color <- rbvertcolors[ beh  ]
    V(gs)$color.frame <- 'black'
    
    V(gs)$size <- 9

    ##label
    vertex.label <- ''
    ##
    file.parts <- c(fileprefix, filename)
    if (filetime)
      file.parts <- c(file.parts, as.integer(Sys.time()))
    file.full <- sprintf('%s.png', paste(file.parts[!is.na(file.parts)], collapse = '_'))
    
    png(filename = file.full, width = 8, height = 8, res = 400, units = 'in')
    par(mar=c(.1,.1,.1,.1))
    set.seed(seed)
    plot(gs, layout=layout.fruchterman.reingold
         , vertex.size=V(gs)$size
         , vertex.color=V(gs)$color
         , vertex.frame.color = V(gs)$color.frame
         , vertex.label=vertex.label
         , vertex.label.cex=.01
         , vertex.label.color=NA
         , vertex.label.font = 2
         , vertex.label.family = 'sans'
         , vertex.shape = 'circle'
         , edge.arrow.width = .1
         , edge.arrow.size = .1
         , edge.width = 2
         , edge.color = 'darkgray'
         , rescale = TRUE
         , ...)
    legend('topright', title = legend.title, legend=sort(unique(beh)), 
           col=rbvertcolors, pch=rep(19, ncl), cex=1.6, pt.cex = 2.5)
    dev.off()
  }
  
  ## RESTRUCTURING AGGRESSIVENSS
  yrsubidx.1 <- 1
  yridx.1 <- yridx[yrsubidx.1]
  mat1 <- mmcarr[,,yrsubidx.1]
  beh1 <- arrNetRw[,yrsubidx.1]
    
  yrsubidx.2 <- length(yridx)
  yridx.2 <- yridx[yrsubidx.2]
  mat2 <- mmcarr[,,yrsubidx.2]
  beh2 <- arrNetRw[,yrsubidx.2]
  
  plotMmcNetColorBeh(mat1, beh1, netwavepds[yridx.1], legend.title = 'Restruct.', seed=1357)
  plotMmcNetColorBeh(mat2, beh2, netwavepds[yridx.2], legend.title = 'Restruct.', seed=13579)
    
  ## INVARIANT AGGRESSIVENSS
  yrsubidx.1 <- 1
  yridx.1 <- yridx[yrsubidx.1]
  mat1 <- mmcarr2[,,yrsubidx.1]
  beh1 <- arrNetIw2[,yrsubidx.1]
  
  yrsubidx.2 <- length(yridx)
  yridx.2 <- yridx[yrsubidx.2]
  mat2 <- mmcarr2[,,yrsubidx.2]
  beh2 <- arrNetIw2[,yrsubidx.2]
  
  plotMmcNetColorBeh(mat1, beh1, netwavepds[yridx.1], legend.title = 'Invariant', seed=1357)
  plotMmcNetColorBeh(mat2, beh2, netwavepds[yridx.2], legend.title = 'Invariant', seed=13579)
  
  
  
 ######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  
  # 
  # 
  # ## NETWORK AUTOCORRELATION
  # MoranGeary <- function(i, data, sims, wave, groupName, varName, levls=1:2){
  #   #unloadNamespace("igraph") # to avoid package clashes
  #   require(sna)
  #   require(network)
  # 
  #   x <- as.sociomatrix(networkExtraction(i, data, sims, wave, groupName, varName[1]))
  #   z <- behaviorExtraction(i,data,sims,wave,groupName,varName[2])
  #   n <- length(z)
  #   z.ave <- mean(z,na.rm=TRUE)
  #   numerator <- n*sum(x*outer(z-z.ave,z-z.ave),na.rm=TRUE)
  #   denominator <- sum(x,na.rm=TRUE)*sum((z-z.ave)^2,na.rm=TRUE)
  #   res <- numerator/denominator
  #   numerator <- (n-1)*sum(x*(outer(z,z,FUN='-')^2),na.rm=TRUE)
  #   denominator <- 2*sum(x,na.rm=TRUE)*sum((z-z.ave)^2,na.rm=TRUE)
  #   res[2] <- numerator/denominator
  #   names(res) <- c("Moran","Geary")
  #   return(res)
  # }
  # 
  # sysResAMC0$a
  # 
  # sims[[i]][[groupName]][[varName[1]]][[period]][, 1]
  # , 
  # sims[[i]][[groupName]][[varName]][[period]][, 2], x = sims[[i]][[groupName]][[varName]][[period]][,3], dims = dimsOfDepVar[1:2]) 
  # 
  #   
  # simMod <- sienaAlgorithmCreate(projname='sys-mmc-acq-test-proj-converg-VDcut-sims-0', simOnly = T, 
  #                                cond = FALSE, useStdInits = FALSE, nsub = 0 ,
  #                                seed=135)
  # 
  # sims  <- siena07(simMod, data=sysDat, effects=sysEff, 
  #                 useCluster = T, nbrNodes = detectCores(), clusterType = 'PSOCK')
  # 
  # i <- 1
  # wave <- 1
  # groupName <- 'Data1'
  # varName <- c('depMMC','depNetR')
  # data <- sysDat
  # mc <- MoranGeary(i, sysResAMC0$f, sims, wave, groupName, varName)
  # 
  #
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  #######################################################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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

