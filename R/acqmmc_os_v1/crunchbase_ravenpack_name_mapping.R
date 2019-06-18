##
##
##  CREATE CRUNCHBASE - TO - RAVENPACK MAPPING
##  
##    - INPUTE: FUZZY NAME MATCHING TABLE (with manual edits)
##    - OUPUT:  FINAL MAPPING (of manually checked names)
##
library(readxl)


## DIRECTORIES
data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
version_dir <- file.path(work_dir,'R','acqmmc_os_v1')
result_dir <- file.path(work_dir,'acqmmc_os_v1','data')

## set woring dir
setwd(work_dir)

##============================
##   ACTION CATEGORIES
##----------------------------
xlsxfile <- 'crunchbase_ravenpack_name_mapping_n187_EDIT.xlsx'

##convert
df.in <- read_xlsx(file.path(result_dir, xlsxfile), na = c('NA','-', ' ','',' - '))

map <- df.in[,c('name')]
map$rp_entity_id <- NA
map$rp_entity_name <- NA
names(map)[which(names(map)=='name')] <- 'cb_company_name_unique'

for (i in 1:nrow(df.in)) {
  USE <- df.in$USE[i]
  if (is.na(USE) | USE=='\'-') {
    cat(sprintf('skipping %s  %s\n', i,df.in$name[i]))
    next
  }
  USE <- as.numeric(USE)
  if (USE == 1) {
    map$rp_entity_id[i] <- df.in$rp_entity_id_1[i]
    map$rp_entity_name[i] <- df.in$entity_name_1[i]
  } else if (USE == 2) {
    map$rp_entity_id[i] <- df.in$rp_entity_id_2[i]
    map$rp_entity_name[i] <- df.in$entity_name_2[i]    
  } else if (USE == 3) {
    map$rp_entity_id[i] <- df.in$rp_entity_id_3[i]
    map$rp_entity_name[i] <- df.in$entity_name_3[i]
  } else if (USE == 4) {
    map$rp_entity_id[i] <- df.in$rp_entity_id_4[i]
    map$rp_entity_name[i] <- df.in$entity_name_4[i]
  }
}


### REMEMBER rp_entity_id & entity_name have multiples separated by pipes "|"
csvmapfile <- 'crunchbase_ravenpack_name_mapping_n187_FINAL.csv'

write.csv(map, file.path(result_dir, csvmapfile), row.names = F)





