##########################################################################################
#
# AMJ 2018 SPECIAL ISSUE COMPETITION NETWORKS AND AWARENESS
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
#.libPaths('C:/Users/sdowning/Documents/R/win-library/3.2')
library(plyr, quietly = T)
library(magrittr, quietly = T)
library(reshape2, quietly = T)
library(lubridate, quietly = T)
library(stringr, quietly = T)


##
#
##
.main.cbdp <- function(x=NA)
{
  
  ## DIRECTORIES
  data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
  work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
  cs_dir <- file.path(work_dir, 'compustat')   ### just using same compustat directly for all wrds data
  version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')
  sup_data_dir <- file.path(work_dir,'amj_rnr2_sup_data')  ## supplmental data dir
  
  ## set woring dir
  setwd(work_dir)
  
  # cb list filename
  cb.list.file <- file.path(sup_data_dir,'cb_list.rds')
  
  ## load full graph, else make full graph if not exists in working directory
  if (file.exists(cb.list.file)) 
  {
    cat('\nloading CrunchBase `cb` list object...')
    cb <- readRDS(cb.list.file)
    cat('done.')
    return(cb)
    
  } 
  
  
  
  ## Data list object
  cb <- list()
  
  
  
  ##=================================
  ##  FUNCTIONS
  ##----------------------------------
  
  ##
  # Return vector of indices of  x which are in y
  # @param [vector] x
  # @param [vector] y
  # @return [vector] length is length(intersect(x,y))
  ##
  cb$match <- function (x,y)
  {
    return(which(x %in% y))
  }
  
  ##
  # Returns vector of values from dataframe df 
  #  by checking for non-NA values each row
  #  in columns (old.col, new.col)
  #  - if neither is non-NA, return default value provided (defaults to NA)
  #  - if 1 is non-NA, return the non-NA value
  #  - if both are non-NA, return the value of column `new.col`
  # @param [dataframe] df
  # @param [character] old.col  The old column name (to be replaced by new.col)
  # @param [character] new.col  The new column name (to replace old.col)
  # @return  [vector]
  ##
  cb$parseNonNa <- function(df, old.col, new.col, default=NA)
  {
    return(unname(
      apply(df[ , c(old.col, new.col)], MARGIN = 1, FUN = function(x){
        idx <- which(!is.na(x))
        switch(as.character(length(idx)),
               '0'=default,## `default` if neither is non-NA
               '1'=x[idx], ## the non-NA value if only 1 is non-NA
               '2'=x[2]    ## the replacement column if both are non-NA
        )
      })
    ))
  }
  
  ##
  # Checks if the argument `x` is a falsy (null, <nan>, <na>, 'NA', or "" empty string)
  # @param [any] x  The element to check
  # @return [bool]  
  ##
  cb$falsy <- function(x)
  {
    if (is.null(x)) 
      return(TRUE)
    return (is.nan(x) | is.na(x) | x=='NA' | x=='')
  }
  
  ##
  # Returns vector of competitive relation beginning date strings (or NA) by the following logic:
  #  - take the latest founding date when both have valid founding dates
  #  - else NA
  # @param [dataframe] df  The dataframe to apply the relationBeganOn logic to
  # @return [character[]]  The vector of beginning date strings
  ##
  cb$relationBeganOn <- function(df)
  {
    cols <- c('company_founded_on','competitor_founded_on')
    for(col in cols){
      if (!(col %in% names(df))) stop(sprintf('dataframe missing column `%s`.',col))
    }
    return(unname(apply(X = df[,cols],MARGIN = 1,FUN = function(x){
      max(x)
    })))
  }
  
  ##
  # Returns vector of competitive relation ending date strings (or NA) by the following logic:
  #  - take the earliest ending date (closed, acquired) if any (company,competitor) are not NA
  #  - else NA
  # @param [dataframe] df  The dataframe to apply the relationEndedOn logic to
  # @return [character[]]  The vector of ending date strings
  ##
  cb$relationEndedOn <- function(df)
  {
    cols <- c('company_closed_on','company_acquired_on','competitor_closed_on','competitor_acquired_on')
    for(col in cols){
      if (!(col %in% names(df))) stop(sprintf('dataframe missing column `%s`.',col))
    }
    return(unname(apply(X = df[,cols],MARGIN = 1,FUN = function(x){
      ifelse(any(!cb$falsy(x)), min(x[!cb$falsy(x)]), NA) 
    })))
  }
  
  ##
  # Loads and returns a dataframe from a CSV file
  # @param [character] filename         The file name (absolute path or relative to current working directory)
  # @param [boolean] header             A flag to include header row
  # @param [character] quote            The quotation character
  # @param [character] na.strings       The character strings to represent <NA> values
  # @param [character] stringsAsFactors A flag to convert character strings to factor type
  # @param [character] fill             A flag to fill missing columns per row
  # @return [dataframe] filename
  ##
  cb$readCsv <- function(filename, header=T, quote='"', na.strings=c('NA',''), stringsAsFactors=F, fill=T)
  {
    return(read.csv(file = filename, header=header, quote=quote, na.strings = na.strings, 
                    stringsAsFactors = stringsAsFactors, fill=fill))
  }
  
  
  
  ##=================================
  ##  Sources
  ##----------------------------------
  cb$csv <- list()
  cb$csv$co          <- 'organizations.csv'
  cb$csv$co_comp     <- 'competitors.csv'
  cb$csv$co_cust     <- 'customers.csv'
  cb$csv$co_parent   <- 'org_parents.csv'
  cb$csv$co_prod     <- 'products.csv'
  cb$csv$co_acq      <- 'acquisitions.csv'
  cb$csv$co_br       <- 'branches.csv'
  cb$csv$co_rou      <- 'funding_rounds.csv'  ## company--funding_round
  cb$csv$co_ipo      <- 'ipos.csv'
  cb$csv$fund        <- 'funds.csv'
  cb$csv$inv         <- 'investors.csv'
  cb$csv$inv_rou     <- 'investments.csv'  ## investor--funding_round
  cb$csv$inv_part    <- 'investment_partners.csv'
  # cb$csv$ev        <- 'events.csv'
  # cb$csv$ev_rel    <- 'event_relationships.csv'
  # cb$csv$categ     <- 'category_groups.csv'
  # cb$csv$jobs      <- 'jobs.csv'
  # cb$csv$ppl       <- 'people.csv'
  # cb$csv$ppl_desc  <- 'people_descriptions.csv'
  
  ## COMPUSTAT SEGMENTS FILE
  cs.segments.file <- file.path(cs_dir,'segments-historical.csv')
  
  ##=================================
  ##  Data import
  ##----------------------------------
  cat('\nloading dataframes...')
  
  co        <- cb$readCsv(file.path(data_dir, cb$csv$co))
  co_comp   <- cb$readCsv(file.path(data_dir, cb$csv$co_comp)) 
  co_cust   <- cb$readCsv(file.path(data_dir, cb$csv$co_cust)) 
  co_parent <- cb$readCsv(file.path(data_dir, cb$csv$co_parent))
  co_prod   <- cb$readCsv(file.path(data_dir, cb$csv$co_prod)) 
  co_acq    <- cb$readCsv(file.path(data_dir, cb$csv$co_acq)) 
  co_br     <- cb$readCsv(file.path(data_dir, cb$csv$co_br)) 
  co_rou    <- cb$readCsv(file.path(data_dir, cb$csv$co_rou)) 
  co_ipo    <- cb$readCsv(file.path(data_dir, cb$csv$co_ipo)) 
  fund      <- cb$readCsv(file.path(data_dir, cb$csv$fund)) 
  inv       <- cb$readCsv(file.path(data_dir, cb$csv$inv)) 
  inv_rou   <- cb$readCsv(file.path(data_dir, cb$csv$inv_rou)) 
  inv_part  <- cb$readCsv(file.path(data_dir, cb$csv$inv_part))  
  # ev        <- cb$readCsv(file.path(data_dir, cb$csv$ev))  
  # ev_rel    <- cb$readCsv(file.path(data_dir, cb$csv$ev_rel))  
  # categ     <- cb$readCsv(file.path(data_dir, cb$csv$categ))  
  # job       <- cb$readCsv(file.path(data_dir, cb$csv$jobs))  
  # ppl       <- cb$readCsv(file.path(data_dir, cb$csv$ppl))  
  # ppl_desc  <- cb$readCsv(file.path(data_dir, cb$csv$ppl_desc))  
  
  
  cat('done.\n')
  cat('cleaning data...')
  
  
  ##==========================
  ## COMPANIES 
  ##--------------------------
  ## update founded_on,closed_on dates  - Jin-Su's Email 2018-04-23
  co.date <- cb$readCsv('founded_on_date_edit/missing_companies_20180330.csv')
  names(co.date) <- c('com_uuid','com_name','com_founded_on_UPDATE','com_closed_on_UPDATE','com_status','note')
  ## fix date formatting  '12/13/2011' ==> '2011-12-31'
  co.date$com_founded_on_UPDATE <- sapply(co.date$com_founded_on_UPDATE, function(x){
    return(ifelse(cb$falsy(x), NA,  as.character(mdy(x))))
  })
  co.date$com_closed_on_UPDATE <- sapply(co.date$com_closed_on_UPDATE, function(x){
    return(ifelse(cb$falsy(x), NA,  as.character(mdy(x))))
  })
  
  ## merge udpatese into company df
  co <- merge(x=co,y=co.date,by.x='company_uuid',by.y='com_uuid',all.x=T,all.y=F)
  ## APPLY UPDATE OLD_COLUMN TO NEW_COLUMN  
  ## FASTER THAN FOR-LOOP WHEN DATAFRAME IS VERY LARGE
  co$company_name <- cb$parseNonNa(co, 'company_name','com_name')
  co$status_update <- cb$parseNonNa(co, 'status','com_status')
  co$founded_on <- cb$parseNonNa(co, 'founded_on','com_founded_on_UPDATE')
  co$closed_on <- cb$parseNonNa(co, 'closed_on','com_closed_on_UPDATE')
  ## DROP REPLACEMENT COLUMNS
  co$com_name <- NULL
  co$com_status <- NULL
  co$com_founded_on_UPDATE <- NULL
  co$com_closed_on_UPDATE <- NULL
  co$note <- NULL
  
  # ##  DROP COMPANY row with mising name
  co <- co[which(co$company_name_unique!="" & !is.na(co$company_name_unique)), ]
  
  ## ONLY COMPANIES (no universities, etc)
  co <- co[which(co$primary_role == 'company'), ]
  
  ## DROP DUPLICATES
  co.cnt <- plyr::count(co$company_name_unique)
  co.cnt <- co.cnt[order(co.cnt$freq, decreasing = T), ]
  # co.cnt[co.cnt$freq > 1, ]
  co.dups <- as.character(co.cnt$x[which(co.cnt$freq > 1)])
  if(length(co.dups) > 0) {
    co <- co[which( !(co$company_name_unique %in% co.dups) ), ]
  }
  
  ## ADD YEAR variables
  co$founded_year <- as.numeric(stringr::str_sub( co$founded_on,1,4))
  co$closed_year <-  as.numeric(stringr::str_sub( co$closed_on,1,4))
  co$acquired_year <-  as.numeric(stringr::str_sub( co$acquired_on,1,4))
  
  ## ADD STOCK SYMBOL FROM IPO TABLE
  ipo_dedup <- ddply(co_ipo, .(company_name_unique), summarize,
                     stock_symbol=min(stock_symbol),
                     stock_exchange_symbol=min(stock_exchange_symbol),
                     went_public_on=min(went_public_on))
  co <- merge(co, ipo_dedup, by.x='company_name_unique',by.y='company_name_unique',all.x=T,all.y=F)
  
  ## ADD KEYS & CODES
  ## - COMPUSTAT GVKEY
  ## - SIC CODES
  ## - CUSIP (for SDC relations)
  csseg <- read.csv(cs.segments.file, header = T, stringsAsFactors = F)
  csply <- ddply(csseg, .(tic), summarize, .progress = 'text',
                 sic=paste(unique(c(SICS1,SICS2,sic)),collapse='|'),
                 naics=paste(unique(c(NAICSS1,NAICSS2,naics)),collapse='|'),
                 cusip=paste(unique(c(cusip)),collapse='|'),
                 gvkey=paste(unique(c(gvkey)),collapse='|'),
                 conm=paste(unique(c(conm)),collapse='|'))
  
  ## CHECK WHICH CB FIRM IPO NOT HAVE SYMBOLS IN COMPUSTAT?
  # idx.cb.ipo <- which(!is.na(cb$co$stock_symbol))
  # cb.ipo.no.cs <- sapply(cb$co$stock_symbol[idx.cb.ipo], function(x) (!is.na(x) & (! x %in% csply$tic)) )
  # idx.cb.ipo.no.cs <- idx.cb.ipo[cb.ipo.no.cs]
  # ##  name                                    symbol    exchange
  # ##  google	                                GOOG	    NASDAQ
  # ##  alibaba-health-information-technology	  241	      SEHK
  # ##  corporate-executive-board-ceb	          EXBD	    NASDAQ
  # ##  servicepower	                          SVR	      LON ## CS MISSING
  # ##  software-ag	                            SOW	      ETR ## CS MISSING
  # ##  wynyard-group	                          WYN	      NZX ## CS MISSING
  df.stockfix <- data.frame(
    company_name_unique=c('google','alibaba-health-information-technology',
                          'corporate-executive-board-ceb',
                          'servicepower','software-ag','wynyard-group'),
    cb_stock_symbol=c('GOOG','241','EXBD','SVR','SOW','WYN'),
    tic=c('GOOGL','BABA','CEB',NA,NA,NA),
    stringsAsFactors = F
  )
  
  ## MANUALLY FIX STOCK SYMBOLS FROM CRUNCHBASE THAT DON'T MATCH COMPUSTAT
  co$tic <- co$stock_symbol
  for (fix_i in 1:nrow(df.stockfix)) {
    idx_fix_i <- which(co$company_name_unique == df.stockfix$company_name_unique[fix_i])
    co$tic[idx_fix_i] <- df.stockfix$tic[fix_i]
  }
  
  ## MERGE IN KEYS & CODES (gvkey, cusip, sic, naics) BY TIC SYMBOL
  co <- merge(co, csply, by.x='tic', by.y='tic', all.x=T, all.y=F)
  
  
  ##==========================
  ## ACQUISITIONS
  ##--------------------------
  ## make unique acquired_on dates for each company
  ##  - as min of acquired_on dates when multiples (acquired multiple times) 
  cat('reshaping acquisitions dataframe...')
  co_acq_acquired <- ddply(co_acq, 'acquiree_name_unique', summarize, .progress='text',
                           acquiree_uuid=paste(unique(acquiree_uuid),collapse="|"),
                           acquired_on=min(acquired_on),
                           acquired_on_concat=paste(acquired_on,collapse="|"))
  co_acq_acquired <- co_acq_acquired[which(!cb$falsy(co_acq_acquired$acquiree_name_unique) & !cb$falsy(co_acq_acquired$acquired_on)),]
  ##
  co_acq$company_name_unique <- co_acq$acquirer_name_unique
  co_acq$acquired_year <- as.numeric(stringr::str_sub(co_acq$acquired_on,1,4))
  
  
  ##==========================
  ## COMPETITIVE RELATIONS
  ##--------------------------
  ## add dates to compute relation_began_on, relation_ended_on
  ## remove acquired_on|closed_on|founed_on dates from python pre-processing script
  ## if these exist:
  drop.cols <- c('company_closed_on','competitor_closed_on','company_founded_on','competitor_founded_on',
                 'relation_ended_on','max_founded_on','relation_began_on',
                 'acquiree_name_unique_concat','acquired_on_concat','acquired_on')
  for (col in drop.cols){
    if (col %in% names(co_comp)) {
      co_comp[,col] <- NULL
    }
  }
  ## rename entity_uuid to company_uuid
  names(co_comp)[which(names(co_comp)=='entity_uuid')] <- 'company_uuid'
  ##  - COMPANY FOUNDED | CLOSED
  tmp <- co[,c('company_uuid','founded_on','closed_on')]
  names(tmp) <- c('company_uuid','company_founded_on','company_closed_on')
  co_comp <- merge(x=co_comp,y=tmp, by.x='company_uuid',by.y='company_uuid', all.x=T,all.y=F)
  ##  - COMPANY ACQUIRED  ## using co_acq_acquired for unique acquiree dates
  tmp <- co_acq_acquired[,c('acquiree_uuid','acquired_on')]
  names(tmp) <- c('acquiree_uuid','company_acquired_on')
  co_comp <- merge(x=co_comp,y=tmp, by.x='company_uuid',by.y='acquiree_uuid', all.x=T,all.y=F)
  ##  - COMPETITOR FOUNDED | CLOSED
  tmp <- co[,c('company_uuid','founded_on','closed_on')]
  names(tmp) <- c('competitor_uuid','competitor_founded_on','competitor_closed_on')
  co_comp <- merge(x=co_comp,y=tmp, by.x='company_uuid',by.y='competitor_uuid', all.x=T,all.y=F)
  ##  - COMPETITOR ACQUIRED  ## using co_acq_acquired for unique acquiree dates
  tmp <- co_acq_acquired[,c('acquiree_uuid','acquired_on')]
  names(tmp) <- c('acquiree_uuid','competitor_acquired_on')
  co_comp <- merge(x=co_comp,y=tmp, by.x='company_uuid',by.y='acquiree_uuid', all.x=T,all.y=F)
  
  ##
  ## keep only intersection of competitive relations by the following conditions:
  ##
  ## 1. company or competitor names is missing
  idx.name.exists <- which(co_comp$company_name_unique != '' & !is.na(co_comp$company_name_unique)
                           & co_comp$competitor_name_unique != '' & !is.na(co_comp$competitor_name_unique) )
  ## 2. only relations where both parties are in the company table 
  idx.name.in.co <- which(co_comp$company_uuid %in% co$company_uuid
                          | co_comp$competitor_uuid %in% co$company_uuid)
  ## 3. only relations where founded_on date is not missing or "no exist"
  idx.date.exists <- which( !cb$falsy(co_comp$company_founded_on) 
                            & co_comp$company_founded_on != "no exist"
                            & !cb$falsy(co_comp$competitor_founded_on) 
                            & co_comp$competitor_founded_on != "no exist" )
  idx.all <- intersect( intersect(idx.name.exists,idx.name.in.co), idx.date.exists )
  co_comp <- co_comp[idx.all, ]
  ## remove previous relation_began_on|relation_ended_on
  co_comp$relation_began_year <- NULL ##  as.numeric(stringr::str_sub( co_comp$relation_began_on,1,4))
  co_comp$relation_ended_year <- NULL ##  as.numeric(stringr::str_sub( co_comp$relation_ended_on,1,4))
  
  ## set relation_began_on|relation_ended_on by new rules use new founed_on dates
  
  ## RELATION BEGAN ON
  co_comp$relation_began_on <- cb$relationBeganOn(co_comp)
  ## RELATION ENDED ON
  co_comp$relation_ended_on <- cb$relationEndedOn(co_comp)
  
  
  
  ##==========================
  ## MULTI_MARKET_CONTACT CODE
  ##--------------------------
  co_br$mmc_code <- apply(co_br[,c('country_code3','region_code2')],1,function(x){
    paste(x, collapse="_")
  })
  
  
  ##==========================
  ## FUNDING ROUND
  ##--------------------------
  co_rou$funded_year <- as.numeric(stringr::str_sub(co_rou$announced_on,1,4))
  
  
  ##==========================
  ## BRANCHES
  ##--------------------------
  co_br$created_year <- as.numeric(stringr::str_sub(co_br$created_at,1,4))
  
  
  ##==========================
  ## IPO
  ##--------------------------
  co_ipo$went_public_year <- as.numeric(stringr::str_sub(co_ipo$went_public_on,1,4))
  
  
  ##==========================
  ## CUSTOMER
  ##--------------------------
  co_cust$created_year <- as.numeric(stringr::str_sub(co_cust$created_at,1,4))
  
  
  ##==========================
  ## Unique Branches
  ##--------------------------
  co_br <- unique(co_br)
  
  
  cat('done.\n')
  
  
  
  ##==========================
  ##
  ##  CLEAN ENVIRONMENT
  ##
  ##--------------------------
  cat('clearing environment...')
  ## dataframes
  cb$co        <- co
  cb$co_comp   <- co_comp
  cb$co_cust   <- co_cust
  cb$co_parent <- co_parent
  cb$co_prod   <- co_prod
  cb$co_acq    <- co_acq
  cb$co_br     <- co_br
  cb$co_rou    <- co_rou
  cb$co_ipo    <- co_ipo
  cb$fund      <- fund
  cb$inv       <- inv
  cb$inv_rou   <- inv_rou 
  cb$inv_part  <- inv_part
  # cb$ev        <- ev
  # cb$ev_rel    <- ev_rel
  # cb$categ     <- categ
  # cb$job       <- job
  # cb$ppl       <- ppl
  # cb$ppl_desc  <- ppl_desc
  
  

  ## SAVE
  ## load full graph, else make full graph if not exists in working directory
  saveRDS(cb, file = cb.list.file)
  
  
  ## RETURN
  return(cb)
}

##
# CRUNCHBASE DATA LIST OBJECT
# @export [list] cb
##
.main.cbdp()


















