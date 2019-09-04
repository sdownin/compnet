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

## DIRECTORIES
data_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/steph/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
version_dir <- file.path(work_dir,'R','acqmmc_os_v1')
net_dir <- file.path(work_dir,'acqmmc_os_v1','data')
result_dir <- file.path(work_dir,'acqmmc_os_v1','data')
sup_data_dir <- file.path(work_dir,'acqmmc_os_v1','sup_data')  ## supplmental data dir

## set woring dir
setwd(work_dir)

rvn_dir <- 'C:\\Data\\RavenPack'

##============================
##   ACTION CATEGORIES
##----------------------------
dtafile <- 'rp40_action_categories_2000_2018.dta'

##convert
rvn <- read_dta(file.path(rvn_dir, dtafile))


###
## Complexity measure of actions 
##  @see Yu, Subramaniam, & Cannella, 2009
###
complexity <- function(x, scale=FALSE) {
  sumx <- sum(x)
  sq.prop <- sapply(x, function(k) (k/sumx)^2 )
  y <- 1/sum(sq.prop)
  if (scale)
    return(y / length(x))
  return(y)
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

## -- settings --
firms.todo <- c('ibm')
d <- 3
yrpd <- 1
startYr <- 2006
endYr <- 2017            ## dropping first for memory term; actual dates 2007-2016
## --------------  
periods <- seq(startYr,endYr,yrpd)
name_i <- 'ibm'


## READ IN DATA FILE
firmfile <- sprintf('acq_sys_mmc_panel_regression_df_col%s_%s_%s-%s.csv',
                    50, name_i, periods[1], periods[length(periods)])
dfl <- read.csv(file.path(result_dir,firmfile))








##=======================================
##  TEST PANEL DATA RERESSION
##---------------------------------------
dfl$name <- as.character(dfl$name)
yrmin <- min(unique(dfl$year))
## get public firms to subset data
ipo.yrmin.name <- dfl$name[which(dfl$ipo==1)]
# ipo.yrmin.name <- V(gmmc)$vertex.names[ which(V(gmmc)$ipo_status == 1) ]
## get firms by any(type==TRUE)
names.type.true <- unique(dfl$name[dfl$type])
## Regression data subset
dfreg <- dfl[which(dfl$name %in% names.type.true & dfl$name %in% ipo.yrmin.name), ]
# dfreg <- dfl[which(dfl$name %in% names.type.true), ]

## compute feedback measure
dfreg$feedback1 <- NA
dfreg$feedback2 <- NA
dfreg$feedback3 <- NA
dfreg$feedback4 <- NA
dfreg$feedback5 <- NA
for (t in unique(dfreg$year)) {
  if (t == min(unique(dfreg$year))) {
    dfreg$feedback1[which(dfreg$year==t)] <- NA
    dfreg$feedback2[which(dfreg$year==t)] <- NA
    dfreg$feedback3[which(dfreg$year==t)] <- NA
    dfreg$feedback4[which(dfreg$year==t)] <- NA
  } else {
    dfreg$feedback1[which(dfreg$year==t)] <- dfreg$smmc1[which(dfreg$year==t)] - dfreg$smmc1[which(dfreg$year==(t-1))] 
    dfreg$feedback2[which(dfreg$year==t)] <- dfreg$smmc2[which(dfreg$year==t)] - dfreg$smmc2[which(dfreg$year==(t-1))] 
    dfreg$feedback3[which(dfreg$year==t)] <- dfreg$smmc3[which(dfreg$year==t)] - dfreg$smmc3[which(dfreg$year==(t-1))] 
    dfreg$feedback4[which(dfreg$year==t)] <- dfreg$smmc4[which(dfreg$year==t)] - dfreg$smmc4[which(dfreg$year==(t-1))] 
  }
}
  
##factors
dfreg$i <- as.factor(dfreg$i)
dfreg$year <- as.factor(dfreg$year)
dfreg$y <- dfreg$y.cur + dfreg$y.new


##========================================
## RavenPack NAME MAPPING BEGIN 
##----------------------------------------
firms <- sort(unique(dfreg$name))
# RAVENPACK NAMES
rpn <- unique(rvn[,c('rp_entity_id','entity_name')])
rpn$name.l <- str_to_lower(rpn$entity_name)
# COHORT OF FIRMS
cbrpmap <- data.frame(
  name = unique(as.character(dfreg$name)),
  idx1 = NA, ## mode
  idx1pct = NA, ## confidence (%)
  rp_entity_id_1 = NA, ## top match ID
  entity_name_1 = NA, ##  top match name
  idx2 = NA, ## 
  idx2pct = NA, ## 
  rp_entity_id_2 = NA, ## 
  entity_name_2 = NA, ## 
  idx3 = NA, ## 
  idx3pct = NA, ## 
  rp_entity_id_3 = NA, ## 
  entity_name_3 = NA, ## 
stringsAsFactors = F)

## string distance matrix for CrunchBase-to-RavenPack Name Mapping
for (i in 1:length(firms)) {
  firm_i <- gsub('\\W+', ' ', firms[i])
  if (i %% 10 == 0) cat(sprintf('i=%03s (%03.2f%s) %s \n',i,100*i/length(firms),'%',firms[i]))
  ##
  ## FIRST DO grep() string search, see if excatly match 1
  ##  - else if multiple matches, test distance on matches
  ##  - else test distance on all unmateched RavenPack firms
  ##
  idf <- data.frame(
    name=rpn$entity_name,
    name.l=rpn$name.l,
    lcs=stringdist(rpn$name.l, firm_i, method = 'lcs'),
    cosine=stringdist(rpn$name.l, firm_i, method = 'cosine'),
    jaccard=stringdist(rpn$name.l, firm_i, method = 'jaccard'),
    jw=stringdist(rpn$name.l, firm_i, method = 'jw'),
    qgram=stringdist(rpn$name.l, firm_i, method = 'qgram'),
    osa=stringdist(rpn$name.l, firm_i, method = 'osa')#,
    # lv=stringdist(rpn$name.l, firm_i, method = 'lv'),
    # dl=stringdist(rpn$name.l, firm_i, method = 'dl'),
  )
  idx <- apply(idf[,3:ncol(idf)], 2, which.min)
  idxcnt <- plyr::count(idx)
  idxcnt <- idxcnt[order(idxcnt$freq, decreasing = T),]
  if (nrow(idxcnt)==0 | length(idxcnt)==0) {
    cat(sprintf('missing strdist for i=%03s, `%s\n`',i,firms[i]))
    next
  }
  ## top match
  idx1 <- idxcnt$x[1]
  cbrpmap$idx1[i] <- idx1
  cbrpmap$idx1pct[i] <- idxcnt$freq[1] / sum(idxcnt$freq)
  cbrpmap$rp_entity_id_1[i] <- rpn$rp_entity_id[idx1]
  cbrpmap$entity_name_1[i] <- rpn$entity_name[idx1]
  ## secondary / tertiary matches
  if (nrow(idxcnt) >= 2) {
    idx2 <- idxcnt$x[2]
    cbrpmap$idx2[i] <- idx2
    cbrpmap$idx2pct[i] <- idxcnt$freq[2] / sum(idxcnt$freq)
    cbrpmap$rp_entity_id_2[i] <- rpn$rp_entity_id[idx2]
    cbrpmap$entity_name_2[i] <- rpn$entity_name[idx2]
  }
  if (nrow(idxcnt) >= 3) {
    idx3 <- idxcnt$x[3]
    cbrpmap$idx3[i] <- idx3
    cbrpmap$idx3pct[i] <- idxcnt$freq[3] / sum(idxcnt$freq)
    cbrpmap$rp_entity_id_3[i] <- rpn$rp_entity_id[idx3]
    cbrpmap$entity_name_3[i] <- rpn$entity_name[idx3]
  }
}

## save mapping dataframe
cbrpmapfile <- sprintf('crunchbase_ravenpack_name_mapping_n%s.csv', length(firms))
write.csv(cbrpmap, file=file.path(result_dir,cbrpmapfile))


##----------------------------------------
## LOAD FINALIZED RAVENPACK NAME MAPPING
##----------------------------------------
cbrpmapfinalfile <- sprintf('crunchbase_ravenpack_name_mapping_n%s_FINAL.csv', length(firms))
cbrpmap <- read.csv(file.path(result_dir,cbrpmapfinalfile),
                    stringsAsFactors = F, na.strings = c('', ' ', 'NA', '-'))

##----------------------------------------
##  MERGE Ravenpack IDs into CrunchBase Data
##----------------------------------------
dfreg <- merge(dfreg, cbrpmap, 
               by.x='name',by.y='cb_company_name_unique',all.x=T,all.y=F)
## Filter out firms that don't have RavenPack ID
dfreg <- dfreg[!is.na(dfreg$rp_entity_id),]

##----------------------------------------
##  Compute Ravenpack ACTIONS per Firm-Year
##----------------------------------------
## init DVs
dfreg$rp_Acquisitions <- 0
dfreg$rp_Capacity <- 0
dfreg$rp_Legal <- 0
dfreg$rp_Market_expansions <- 0
dfreg$rp_Marketing <- 0
dfreg$rp_New_product <- 0
dfreg$rp_Pricing <- 0
dfreg$rp_Strategic_alliances <- 0
##
dfreg$rp_NON_acquisitions <- 0
dfreg$rp_all <- 0
dfreg$rp_complexity <- 0
##
# .rptmp <- data.frame(
#   rp=Acquisitions='rp_Acquisitions',
#   Capacity_related_actions='rp_Capacity',
#   Legal_actions='rp_Legal',
#   Market_expansion='rp_Market_expansions',
#   Marketing_actions='rp_Marketing',
#   New_Product_Action='rp_New_product',
#   pricing_actions='rp_Pricing',
#   strategic_alliance='rp_Strategic_alliances',
# )
.rptmp <- data.frame(
  rp = c('Acquisitions', 'Capacity_related_actions', 
         'Legal_actions', 'Market_expansion', 'Marketing_actions',
         'New_Product_Action', 'pricing_actions',
         'strategic_alliance'),
  dfreg= c('rp_Acquisitions','rp_Capacity','rp_Legal',
           'rp_Market_expansions','rp_Marketing',
           'rp_New_product','rp_Pricing',
           'rp_Strategic_alliances'),
  val = rep(0,8), 
  stringsAsFactors = F
)
## COMPUTE ACTION DVS FOR EACH FIRM_YEAR OBS
for (i in 1:nrow(dfreg)) {
  x <- dfreg[i,]
  year <- as.numeric(as.character(x$year))
  idx <- which(rvn$year==year & rvn$rp_entity_id==x$rp_entity_id)
  sub <- rvn[idx,]
  if (length(sub)>0 & nrow(sub)>0) {
    cnt <- plyr::count(sub$action_category)
    cnt <- merge(cnt,.rptmp[,c('rp','val')],by.x = 'x',by.y='rp',all=T)
    cnt$freq <- apply(cnt[,c('freq','val')],1,function(x)max(x,na.rm=T))
    cnt$val <- NULL
    ## Assign simple count values by action type
    for (k in 1:nrow(.rptmp)) {
      rp_action <- .rptmp$rp[k]
      dfreg_action <- .rptmp$dfreg[k]
      dfreg[i,dfreg_action] <- cnt$freq[which(cnt$x==rp_action)]
    }
    ## Compute covariates of action distributions
    dfreg$rp_complexity[i] <- complexity(cnt$freq)
    dfreg$rp_all[i] <- sum(cnt$freq, na.rm = T)
    dfreg$rp_NON_acquisitions[i] <- dfreg$rp_all[i] - dfreg$rp_Acquisitions[i]
  }
}

###
### ****************** HERE ******************
###



##========================================
## COEF NAME MAPPING
##----------------------------------------
coef.map <- list(
  `(Intercept)` = '(Intercept)',
  sigma = '(Random effect SD)',
  acq_sum_1_sc = 'Prior Yr Acquisition Sum ($ Mn.)',
  employee_na_age_sc = 'Employees',
  sales_na_0_mn_sc = 'Sales ($ Mn.)',
  cent_deg_all_sc = 'Competitors',
  dum.crisis='Financial Crisis Dummy',
 `I(acq_cnt_5 > 0)TRUE`='Acquisition Experience (5 yr binary)',
 `log(1 + acq_cnt_5)`='Acquisition Experience (ln 5 yr count + 1)',
 `I(acq_sum_1/1e+09)`='Prior Yr Acquisition Sum ($ Bn.)',
 `I(acq_sum_1/1e+06)`='Prior Yr Acquisition Sum ($ Mn.)',
 `I(employee_na_age/1000)`='Employees (1000)',
 `I(sales_na_0_mn/1e+06)`='Sales ($ Mn.)',
 `log(1 + cent_deg_all)`='Ln Competitors',
 `cent_deg_all`='Competitors',
 `smmc1n`='System MMC',
 `I(smmc1n^2)`='System MMC Squared [H1]',
 `pres1n`='Competitive Pressure',
 `feedback1`='System Feedback',
 `smmc1n:pres1n`='System MMC * Pressure',
 `smmc1n:feedback1`='System MMC * Feedback',
 `lag(feedback1, 0:0)`='System Feedback',
 `lag(feedback1, 0:0)0`='System Feedback',
 `lag(feedback1, 0:1)0`='System Feedback',
 `lag(feedback1, 0:1)1`='System Feedback (Lag 1)',
 `I(smmc1n^2):pres1n`='System MMC Squared * Pressure [H2]',
 `pres1n:I(smmc1n^2)`='System MMC Squared * Pressure [H2]',
 `I(smmc1n^2):lag(feedback1, 0:0)`='System MMC Squared * Feedback [H3]',
 `lag(feedback1, 0:0):I(smmc1n^2)`='System MMC Squared * Feedback [H3]',
 `I(smmc1n^2):feedback1`='System MMC Squared * Feedback [H3]',
 `feedback1:I(smmc1n^2)`='System MMC Squared * Feedback [H3]',
 `poly(smmc1n, 2)1` = 'System MMC',
 `poly(smmc1n, 2)2` = 'System MMC Squared [H1]',
 `poly(smmc1n, 2)1:pres1n` = 'System MMC * Pressure',
 `poly(smmc1n, 2)2:pres1n` = 'System MMC Squared * Pressure [H2]',
 `poly(smmc1n, 2)1:feedback1` = 'System MMC * Feedback',
 `poly(smmc1n, 2)2:feedback1` = 'System MMC Squared * Feedback [H3]'
 )

## MUTATE CONTROLS
dfreg$dum.crisis <- 0
dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
dfreg$acq_sum_1_sc <- scale(dfreg$acq_sum_1, center = T, scale = T)
dfreg$employee_na_age_sc <- scale(dfreg$employee_na_age, center = T, scale = T)
dfreg$sales_na_0_mn_sc <- scale(dfreg$sales_na_0_mn, center = T, scale = T)
dfreg$cent_deg_all_sc <- scale(dfreg$cent_deg_all, center = T, scale = T)

## remove NAs in first year due to feedback lag
dfreg <- dfreg[!is.na(dfreg$feedback1),]



##===================================
## RAVENPACK COMPARE Acquisition DV - Actions Moderator
##-----------------------------------
mract <- pglm(rp_Acquisitions ~  dum.crisis + I(acq_cnt_5 > 0) + 
                I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                lag(rp_NON_acquisitions,1) +
                smmc1n +  # pres1n +  feedback1 + 
                #smmc1n:pres1n +  smmc1n:feedback1 + 
                smmc1n:lag(rp_NON_acquisitions,1) + 
                I(smmc1n^2) +
                I(smmc1n^2):lag(rp_NON_acquisitions,1),
              data=dfreg, family = poisson, 
              model = 'within', effect = 'twoways',
              R = 100, method='nr',
              index=c('i','year'))
summary(mract)

mractl1 <- pglm(rp_Acquisitions ~  dum.crisis + I(acq_cnt_5 > 0) + 
                I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                lag(rp_NON_acquisitions,0:1) +
                smmc1n +  # pres1n +  feedback1 + 
                #smmc1n:pres1n +  smmc1n:feedback1 + 
                smmc1n:lag(rp_NON_acquisitions,0:1) + 
                I(smmc1n^2) +
                I(smmc1n^2):lag(rp_NON_acquisitions,0:1),
              data=dfreg, family = poisson, 
              model = 'within', effect = 'twoways',
              R = 100, method='nr',
              index=c('i','year'))
summary(mractl1)

mractl1 <- pglm(rp_Acquisitions ~  dum.crisis + I(acq_cnt_5 > 0) + 
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                   smmc1n +
                   I(smmc1n^2),
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
summary(mractl1)

mractall <- pglm(rp_Acquisitions ~  dum.crisis + I(acq_cnt_5 > 0) + 
                  I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                  I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                  lag(rp_NON_acquisitions,1) +
                  smmc1n +  # pres1n +  feedback1 + 
                  #smmc1n:pres1n +  smmc1n:feedback1 + 
                  smmc1n:lag(rp_NON_acquisitions,1) + 
                  I(smmc1n^2) +
                  I(smmc1n^2):lag(rp_NON_acquisitions,1),
                data=dfreg, family = poisson, 
                model = 'within', effect = 'twoways',
                R = 100, method='nr',
                index=c('i','year'))
summary(mractall)


mractall <- pglm(rp_Acquisitions ~  dum.crisis + I(acq_cnt_5 > 0) + 
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                   lag(rp_NON_acquisitions,1) +
                   smmc1n +  # pres1n +  feedback1 + 
                   #smmc1n:pres1n +  smmc1n:feedback1 + 
                   smmc1n:lag(rp_NON_acquisitions,1) + 
                   I(smmc1n^2) +
                   I(smmc1n^2):lag(rp_NON_acquisitions,1),
                 data=dfreg, family = poisson, 
                 model = 'within', effect = 'twoways',
                 R = 100, method='nr',
                 index=c('i','year'))
summary(mractall)






# mrcomp.list <- list(mrall,mrna, mracq,mrc,mrl,mrme,mrmkt,mrnp,mrp,mrsa)
# model.names <- c('All','Non-Acq',   'Acquisition',
#                  'Capacity','Legal','Mkt.Expan.','Marketing','Product',
#                  'Pricing','Alliance')
# # screenreg(lmn.list, digits = 3)
# screenreg(mrcomp.list, digits = 3, 
#           custom.coef.map = coef.map,
#           custom.model.names = model.names)
# 
# saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_pglm_poisson_list_ravenpack_H1_compare_DV.rds'))
# texreg::htmlreg(mrcomp.list,
#                 file.path(result_dir, 'acqmmc_pglm_poisson_list_ravenpack_H1_compare_DV.html'),
#                 digits = 3, custom.coef.map = coef.map,
#                 custom.model.names = model.names)








