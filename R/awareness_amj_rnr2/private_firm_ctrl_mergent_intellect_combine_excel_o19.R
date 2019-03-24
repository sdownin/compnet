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
library(reshape2)

## dir names
proj_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
mi_dir <- file.path(proj_dir,"private_firm_controls","Mergent_intellect_o19")
pub_dir <- file.path(mi_dir,'public_excel')
pri_dir <- file.path(mi_dir,'private_excel')
data_dirs <- c(pub_dir, pri_dir)
version_dir <- file.path(proj_dir,'R','awareness_amj_rnr2')
sup_data_dir <- file.path(proj_dir,'amj_rnr2_sup_data')  ## supplmental data dir

## set working dir
setwd(proj_dir)


##============================================
## Load firm-filename mapping df 
##  - *reload after any changes to mapdf file
##--------------------------------------------
mapdf.file <- file.path(sup_data_dir,"private_firm_financials_19other_focal_firm_networks_JS.xlsx")
mapdf <- read_excel(mapdf.file, sheet = 1, na=c("-",""), skip = 0)


##============================================
##
## PRIVATE FIRMS DATA
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

## Convert numeric string to numeric by removing commas and dollar signs
##   example:  "$99,123,123.00" --> 99123123.00
str2num <- function(str) {
  if (is.null(str)) return(NA)
  if (is.na(str) | str=='') return(NA)
  x <- str_replace_all(str,'[^\\w\\.]','')
  if (!is.na(x) & x != '') return(as.numeric(x))
  return(NA)
}

##============================================
## Initialize data list
##--------------------------------------------

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


##============================================
## Precheck for All files
##--------------------------------------------

for (.data_dir in data_dirs)
{
  ## Excel files
  files <- dir(.data_dir, pattern = "\\.xlsx{0,1}$")
  ## load files loop
  cat(sprintf("checking files in dir: %s\n", .data_dir))
  for (file in files) {
    file_abr <- str_replace_all(file,'\\.xlsx{0,1}$','')
    firm <- mapdf$firm[which(mapdf$mergent_intellect == file_abr)]
    ## check for firm file
    if (length(firm)==0) {
      cat(sprintf(' - missing firm file:  %s\n',file_abr))
    } else if (length(firm) > 1) {
      cat(sprintf(' - - multiple firms `%s` for file:  %s\n',paste(firm,collapse = "|"),file_abr))
    }
  }
}; cat('done.\n')


##============================================
## Load and combine data files
##--------------------------------------------

## init list of combined dataframes
l <- list()

## loop directories in project
for (.data_dir in data_dirs)
{
  
  ## Excel files
  files <- dir(.data_dir, pattern = "\\.xlsx{0,1}$")
  
  ## load files loop
  for (file in files) 
  {
    cat(sprintf("loading file: %s\n", file))
    
    file_abr <- str_replace_all(file,'\\.xlsx{0,1}$','')
    firm <- mapdf$firm[which(mapdf$mergent_intellect == file_abr)]
    
    ## check firm for file
    stars <- '\n\n**************************************************************\n\n'
    if (length(firm)==0) {
      cat(sprintf('%s - missing firm for file:  %s%s',stars,file_abr,stars))
      next
    } else if (length(firm) > 1) {
      cat(sprintf('%s - - multiple firms `%s` for file:  %s%s',stars,paste(firm,collapse = "|"),file_abr,stars))
      next
    }
    
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
        
        ## set crunchbase is_pub attribute separately from mergent intellect's attribute
        is_pub <-  mapdf$is_pub[which(mapdf$mergent_intellect == file_abr)]
        l[[firm]]$is_pub <- ifelse(is_pub==1 | is_pub=='1', 1, 0)
        
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
        
        ## Skip if no fixed effect table exists on sheet
        if (!any(sapply(fixeff.attrs, function(pattern) greplDf(df,pattern) )))
          next
        
        ##==============================================
        ##  Main data tables extraction logic
        ##----------------------------------------------
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
            dfj <- read_excel(file_path, sheet = sheets[i], na=c("-",""), col_names = T, 
                              skip=(j-1)*(num_years+1), n_max = j*num_years)
            
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
        ##----------------------------------------------
        ##  /end main logic
        ##----------------------------------------------
        
      } ## /end if-else sheet[i]
      
      
    } ##/end sheets loop in file
    
  } ##/end files loop in directory

} ##/end directories loop in project

# warnings()
# print(head(summary(l))); print(length(l))
# 
# ## save data list as serialized file
# saveRDS(l, file = 'awareness_six_focal_firms_mergent_intellect_private_firm_controls.rds')


##============================================
## Combine firm dataframes and melt to longform
##--------------------------------------------
ldf <- data.frame()
for (k in 1:length(l)) {
  ldf <- rbind(ldf, l[[k]]$df)
}
ldf <- as_tibble(ldf)

##============================================
## Convert USD strings --> numeric column format
##--------------------------------------------
cols.fix <- c('employee_all','employee_site','sales')
for (col in cols.fix) {
  ldf[,col] <- sapply(ldf[[col]], str2num)  ## tibble column to vector : tbl[[col]] --> c(..., ...)
}

## save dataframe as CSV
out.file <- file.path(sup_data_dir,'awareness_focal_firms_mi_size_controls_o19.csv')
write.csv(ldf, file = out.file, row.names = F)





