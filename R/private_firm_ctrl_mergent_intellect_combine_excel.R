####################################################
##
##  Compnet Awareness -- Private Firm Data
##
##   - Combine excel files extracted from PDFs
##   - PDF File Source: Mergent Intellect
##
####################################################

library(stringr)
library(readxl)
library(tibble)
library(uuid)
library(dplyr)

## dir names
proj_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
mi_dir <- file.path(proj_dir,"private_firm_controls","mergentintellect")
pub_dir <- file.path(mi_dir,'public_excel')
pri_dir <- file.path(mi_dir,'private_excel')

## set working dir
setwd(proj_dir)

## firm name to data file mapping
mapdf.file <- file.path(proj_dir,"private_firm_controls","private_firm_financials_six_focal_firm_networks_v2.xlsx")
mapdf <- read_excel(mapdf.file, sheet = 1, na=c("-",""), skip = 0)


##============================================
##
## PUBLIC FIRMS DATA
##
##--------------------------------------------

##============================================
## Check number of sheets in excel files
##--------------------------------------------
mapdf$num_sheets <- NA
idx  <- which(!is.na(mapdf$mergent_intellect))
for (i in idx) {
  mapdf$mergent_intellect[i]
  ## file absolute path
  dir <- ifelse( !is.na(mapdf$is_pub[i]) & mapdf$is_pub[i]==1, pub_dir, pri_dir)
  file_i <- file.path(dir,sprintf("%s.xlsx",mapdf$mergent_intellect[i]))
  ## get number of sheets if file exists
  if (file.exists(file_i)) {
    sheets <- excel_sheets(file_i)
    sheets <- sheets[!grepl("^Sheet",sheets)]
    mapdf$num_sheets[i] <- length(sheets)
  }
  cat(sprintf('%s %s\n',which(i==idx),mapdf$mergent_intellect[i]))
}

##============================================
## Helper Functions
##--------------------------------------------

## Simple dataframe init function
initDf <- function(firm=NA, years=2006:2018) {
  return(
    data.frame(year=years, firm=firm, employee_all=NA, employee_site=NA, sales=NA)
  )
}

## Check if dataframe has pattern in any row,col
greplDf <- function(df, pattern='Sales') {
  x <- paste(unlist(df), collapse="|")
  return(grepl(pattern, x, ignore.case = T, perl = T))
}

##============================================
## Initialize data list
##--------------------------------------------

## set tmp data_dir for this loop
.data_dir <- pri_dir

## Excel files
files <- dir(.data_dir, pattern = "\\.xlsx{0,1}$")

# ## lines to skip for Medtrack header atop data table
# skip.lines <- 1

## fields to save for all files (regardless of sheet count)
base.attrs <- c('D-U-N-S Number','Subsidiary Status','Primary SIC Code','Primary NAICS Code',
                'Latitude','Longitude','PhysicalAddress','Zip code',
                'Year of Founding','Company Type','Manufacturer',
                'Employees (All Sites)','Employees (This Site)','Employees Total (Year 1)',
                'Sales')
## Fixed effects data attrs as exported table column titles (excluding sheet 1)
fixeff.attrs <- c('Sales volume','Employees This Site','Employees All Sites')

fixeff.map <- function(x) {
  switch(x,
   `Employees This Site`='employee_site',
   `Employees (This Site)`='employee_site',
   `Employees All Sites`='employee_all',
   `Employees (All Sites)`='employee_all',
   `Sales volume`='sales',
   `Sales`='sales'
  )
}


## init list of combined dataframes
l <- list()

##============================================
## Load and combine data files
##--------------------------------------------

## load files loop
for (file in files) {
  cat(sprintf("loading file: %s\n", file))
  
  file_abr <- str_replace_all(file,'\\.xlsx{0,1}$','')
  firm <- mapdf$firm[which(mapdf$mergent_intellect == file_abr)]
  if (length(firm)==0) cat(sprintf('missing firm for file:  %s\n',file_abr))
  
  ## get sheets (excluding unedited sheets with default name "Sheet__")
  sheets <- excel_sheets(file.path(.data_dir, file))
  sheets <- sheets[!grepl("^Sheet",sheets)]
  num_sheets <- length(sheets)
  
  ## absolute path of data file
  file_path <- file.path(.data_dir, file)
  
  ## loop sheets in file
  for (i in 1:length(sheets)) 
  {
    cat(sprintf('  sheet %s %s\n',i,sheets[i]))
    if (i == 1) 
    {
      ## col_names = FALSE
      df <- read_excel(file_path, sheet = sheets[i], na=c("-",""), col_names = F)
      
      ## FIRST SHEET:  init df; set firm base attrs
      l[[firm]] <- list(df = initDf(firm))
      for (attr in base.attrs) {
        val <- unname(unlist(df[which(df[,1] == attr),2]))
        l[[firm]][[attr]] <- ifelse(length(val)>0, val, NA)
      }
      
      ## set last year(s) df values from 1st sheet base attrs
      n <- nrow(l[[firm]]$df)
      l[[firm]]$df$employee_all[n] <- l[[firm]]$`Employees (All Sites)`
      l[[firm]]$df$employee_site[n] <- l[[firm]]$`Employees (This Site)`
      l[[firm]]$df$employee_all[n-1] <- l[[firm]]$`Employees Total (Year 1)`
      l[[firm]]$df$sales[n] <- l[[firm]]$Sales
      
    } else {
      
      ## after sheet 1, identify if contains data, else skip
      ## col_names = TRUE
      df <- read_excel(file_path, sheet = sheets[i], na=c("-",""), col_names = T)
      
      ## check any fix effect table exists on sheet
      if (!any(sapply(fixeff.attrs, function(pattern) greplDf(df,pattern) )))
        next
      
      ## number of years of data (number of rows of data within df wrongly containing multiple tables)
      num_years <- 0
      num_tables <- 0
      if ('Year' %in% names(df)) {
        idx.yr <- which(df$Year %in% as.character(l[[firm]]$df$year))
        num_years <- length(unique(df$Year[idx.yr]))
        num_tables <- max(plyr::count(df$Year[idx.yr])$freq)
      }
      
      ## separate tables (if num_tables >= 1)
      if (num_tables > 0) {
        for (j in 1:num_tables) 
        {
          
          ## j'th table (subset of rows) on current sheet of df
          dfj <- read_excel(file_path, sheet = sheets[i], na=c("-",""), 
                            skip=(j-1)*(num_years+1), n_max = j*num_years, col_names = T)
          
          ## fill in firm dataframe fixed effects for the columns in current table
          for (eff in fixeff.attrs) {
            if (eff %in% names(dfj)) {
              col.eff <- fixeff.map(eff)
              row.yrs <- which( l[[firm]]$df$year %in% dfj$Year )
              l[[firm]]$df[row.yrs,col.eff] <- dfj[[eff]]
            }
          }
          
        }
      }

      
      df <- read_excel(file_path, sheet = sheets[i], na=c("-",""), n_max = num_years, col_names = T)
      
    }
    
  } ##/end sheets loop in file
  
} ##/end files loop

print(summary(l))





