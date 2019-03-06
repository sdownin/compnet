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
library(dbplyr)

## dir names
proj_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
mi_dir <- file.path(proj_dir,"private_firm_controls","mergentintellect")
pub_dir <- file.path(mi_dir,'public_excel')
pri_dir <- file.path(mi_dir,'private_excel')

## set working dir
setwd(proj_dir)

##============================================
##
## PUBLIC FIRMS DATA
##
##--------------------------------------------

##============================================
## Initialize data list
##--------------------------------------------

## Excel files
files <- dir(pub_dir, pattern = "\\.xlsx{0,1}$")

# ## lines to skip for Medtrack header atop data table
# skip.lines <- 1

## init list of combined dataframes
l <- list()
for (sheet in sheets) {
  l[[sheet]] <- data.frame()
}


##============================================
## Load and combine data files
##--------------------------------------------

## load files loop
for (file in files) {
  cat(sprintf("loading file: %s\n", file))
  
  
  ## get sheets (excluding unedited sheets with default name "Sheet__")
  sheets <- excel_sheets(file.path(medtrack_product_dir, file))
  sheets <- sheets[!grepl("^Sheet",sheets)]
  
  ## extract category from filename 
  ## following last underscore "_"; preceding file extension ".xls|x"
  parts <- str_split(file,"[_]")[[1]]
  category <- str_split(parts[length(parts)],"[.]")[[1]][1]
  
  ## loops sheets in workbook
  for (sheet in sheets) {
    cat(sprintf("  sheet %s\n", sheet))
    
    ## absolute path of data file
    file_full_path <- file.path(medtrack_product_dir, file)
    
    ## already deleted header from Product Synopsis sheet 
    ##  but not from other sheets in workbook;
    ##  skip 11 lines of header material in other sheets
    skip.lines <- ifelse(sheet == "Product Synopsis", 0, 11)
    
    ## load data
    df <- read_excel(file_full_path, sheet = sheet, na="--", skip = skip.lines)
    
    ## clean column names
    names(df) <- str_to_lower(str_replace_all(names(df),"[\\s\\/]+","_"))
    
    ## add category
    df$therapeutic_class <- category
    
    ## append rows to combined dataframe
    l[[sheet]] <- rbind(l[[sheet]], df)
  }
  
}

print(summary(l))