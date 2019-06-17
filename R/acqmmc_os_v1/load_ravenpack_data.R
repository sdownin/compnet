
library(haven)
library(rio)

data_dir <- 'D:\\RavenPack'

##============================
##   ACTION CATEGORIES
##----------------------------
dtafile <- 'rp40_action_categories_2000_2018.dta'
csvfile <- 'rp40_action_categories_2000_2018.csv'
# xlsxfile <- 'rp40_action_categories_2000_2018.xlsx'


##convert
convert(file.path(data_dir, dtafile),
        file.path(data_dir, csvfile))

# df <- haven::read_dta(file.path(data_dir,dtafile))

##============================
##   RAVENPACK FULL DATA
##----------------------------
dtafile <- 'rp40_2000_2016.dta'
csvfile <- 'rp40_2000_2016.csv'
xlsxfile <- 'rp40_2000_2016.xlsx'

##convert
convert(file.path(data_dir, dtafile),
        file.path(data_dir,xlsxfile))

# df <- haven::read_dta(file.path(data_dir,dtafile))