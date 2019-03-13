##--------------------------------------------------------------
##
##    INSTITUTIONAL HOLDINGS 
##     - for public companies computing dyadic shared investors
##     - - (private companies use CrunchBase investment rounds and investors)
##
##--------------------------------------------------------------

.main.institutional.holdings <- function(x=NA)
{
  
  ## DIRECTORIES
  data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
  work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
  img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
  version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')
  wrds_dir <- file.path(work_dir, 'compustat')   ### just using same compustat directly for all wrds data
  
  ## set woring dir
  setwd(work_dir)
  
  # ## LOAD AND PREP CRUNCHBASE DATA IF NOT IN MEMORY
  # if(!('cb' %in% ls())) 
  #   source(file.path(version_dir,'amj_cb_data_prep.R'))
  # 
  # ## CACHE ENVIRONMENT to keep when clearing tmp objects added here
  # ## excluding directory variables ending in `_dir`
  # .ls <- ls()[grep('(?<!_dir)$',ls(),perl = T)]
  
  
  cat('loading institutional holdings data...')
  
  
  ## instituational holdings filename
  ih.combined.file <- file.path(wrds_dir, 'wrds_institutional_holdings_2007-2017.csv')
  
  if (file.exists(ih.combined.file)) {
    
    ih <- read.csv(ih.combined.file, stringsAsFactors = F)
    
  } else {
    
    ## all instituational holdings records
    ih.all <- data.frame()
    
    for (yr in 2007:2017) {
      cat(sprintf(' loading institutional holdings for year %s\n', yr))
      ## LOAD Instituation Holdings Data
      ih.file <- file.path(wrds_dir, sprintf('wrds_institutional_holdings_%s.csv',yr))
      ih <- read.csv(ih.file, header = T, stringsAsFactors = F)
      ##
      ih <- ih[which(ih$ticker %in% cb$co$tic), ]
      ih$cusip_6 <- substr(ih$cusip, 1, 6)
      ih$year <- as.integer(substr(ih$rdate, 1, 4))
      ih$month <- as.integer(substr(ih$rdate, 5, 6))
      ## subset 4th quarter /  12th month
      ih <- ih[which(ih$month == 12), ]
      ## drop columns 
      cols.drop <- c('shrout1','shrout2','prdate')
      ih <- ih[, ! names(ih) %in% cols.drop]
      ## append to all Instituation Holdings records
      ih.all <- rbind(ih.all, ih)
    }
    ih <- ih.all
    
    ## SAVE COMBINED File
    write.csv(ih, file = ih.combined.file, row.names = F)
    
  }
  
  return(ih)
  
}

##------------------------------------
## EXPORT `mi` object with firm size controls 
##  from mergent intellect and compustat combined
##-------------------------------------
.main.institutional.holdings()
