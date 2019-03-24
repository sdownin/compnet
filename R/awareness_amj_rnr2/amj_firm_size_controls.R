##--------------------------------------------------------------
##
##        FIRM SIZE CONTROLS DATA
##         - Employees
##         - Sales
##      
##         : private & public companies (Year>=2010) -- Mergent Intellect
##         : public companies (Year<2010) -- Compustat Annual Fundamentals
##
##--------------------------------------------------------------

.main.size.ctrl <- function(x=NA)
{
  
  ## DIRECTORIES
  data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
  work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
  img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
  cs_dir <- file.path(work_dir, 'compustat')
  version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')
  sup_data_dir <- file.path(work_dir,'amj_rnr2_sup_data')  ## supplmental data dir
  
  ## set woring dir
  setwd(work_dir)
  
  # ## LOAD AND PREP CRUNCHBASE DATA IF NOT IN MEMORY
  # if(!('cb' %in% ls())) 
  #   source(file.path(version_dir,'amj_cb_data_prep.R'))
  # 
  # ## CACHE ENVIRONMENT to keep when clearing tmp objects added here
  # ## excluding directory variables ending in `_dir`
  # .ls <- ls()[grep('(?<!_dir)$',ls(),perl = T)]
  
  
  cat('loading firm size controls...')
  
  
  ## R script to create mergent intellect file if not exists
  mi.file.R.script <- 'private_firm_ctrl_mergent_intellect_combine_excel.R'
  ## mergent intellect file 
  
  mi.files <- c(
    'awareness_focal_firms_mi_size_controls.csv',
    'awareness_focal_firms_mi_size_controls_o19.csv'
  )
  
  ## init combined dataframe for all files in mi.files
  miall <- data.frame()
  
  ## loop each file in mi.files to combine all mergent intellect private firm controls
  for (mi.file in mi.files)
  {
      mi.file.path <- file.path(sup_data_dir, mi.file)
    
    if (!file.exists(mi.file.path))
      stop(sprintf('Mergent Intellect firm data file `%s` does not exist. First run script to create it: %s',
                   mi.file, mi.file.R.script))
    
    ### load Mergent Intellect data file
    mi <- read.csv(mi.file.path, na.strings = c('NA','','-','--'), stringsAsFactors = F)
    
    ## Add Tic symbol
    mi <- merge(mi, cb$co[,c('company_name_unique','tic')], by.x='firm',by.y='company_name_unique',all.x=T,all.y=F)
    
    ## LOAD COMPUSTAT ANNUAL FUNDAMENTALS
    cs.funda.file <- file.path(cs_dir, 'segments-historical.csv')  ## 'segments.csv'; 'fundamentals-annual-RNR2.csv'
    cs <- read.csv(cs.funda.file, header = T, stringsAsFactors = F)
    
    ## set year var
    cs$data_year <- as.integer(sapply(cs$datadate,function(x)substr(x,1,4)))
    # cs$src_year <- as.integer(sapply(cs$srcdate,function(x)substr(x,1,4)))
    
    # ## merge in stock tic symbol
    # cs <- merge(cs, cb$co[,c('gvkey','tic')], by.x='gvkey', by.y='gvkey',all.x=T, all.y=F)
    
    ### public firms
    tics <- unique(mi$tic[!is.na(mi$tic)])
    
    ## loop to fill in public firm data for Year<2010
    for (tic in tics) {
      # cs_tic <- cs[which(cs$tic==tic & cs$stype=='BUSSEG' & cs$sid==1), ]
      cs_tic <- cs[which(cs$tic==tic), ]
      years <- mi$year[which(mi$tic==tic)]
      
      for (year in years) {
        
        if (year %in% cs_tic$data_year) {
          
          idx.mi <- which(mi$tic==tic & mi$year==year)
          
          cs_yr <- cs_tic[which(cs_tic$data_year==year), ]
          sids <- unique(cs_yr$sid)
          
          .sales <- 0
          .emps <- 0
          
          for (sid in sids) {
            cs_yr_sid <- cs_yr[which(cs_yr$sid==sid), ]
            ## ADD SALES if exists
            if (!all(is.na(cs_yr_sid$sales))) {
              sid_sales <- max(cs_yr_sid$sales, na.rm = T)
              if (!is.na(sid_sales)) {
                .sales <- .sales + sid_sales *1e6  ## COMPUSTAT LISTS SEGMENT SALES IN MILLIONS
              }
            }
            ## ADD EMPLOYEES
            if (!all(is.na(cs_yr_sid$emps))) {
              sid_emps <- max(cs_yr_sid$emps, na.rm = T)
              if (!is.na(sid_emps)) {
                .emps <- .emps + sid_emps  ## COMPUSTAT LISTS SEGMENT EMPLOYEES IN MILLIONS
              }
            }
          }
          
          ## UPDATE VALUES IF EXIST AND POSITIVE 
          if (.sales > 0) mi$sales[idx.mi] <- .sales
          if (.emps > 0)  mi$employee_all[idx.mi] <- .emps
          
        }
        
      }
      
    }
    
    ## 2018 data not complete -- drop rows here
    mi <- mi[which(mi$year < 2018), ]
    
    ## combine mergent intellect dataframes per file
    miall <- rbind(miall, mi)
  }
  
  
  ##RETURN combined dataframe
  return(miall)
}


  
##------------------------------------
## EXPORT `mi` object with firm size controls 
##  from mergent intellect and compustat combined
##-------------------------------------
.main.size.ctrl()


