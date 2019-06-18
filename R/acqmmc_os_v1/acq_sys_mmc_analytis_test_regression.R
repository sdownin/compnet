##--------------------------------------------------------------
##
##  MMC & ACQUISITIONS 
##
##--------------------------------------------------------------
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
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
data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/crunchbase_export_20161024"
work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"
version_dir <- file.path(work_dir,'R','acqmmc_os_v1')
net_dir <- file.path(work_dir,'acqmmc_os_v1','data')
result_dir <- file.path(work_dir,'acqmmc_os_v1','data')
sup_data_dir <- file.path(work_dir,'acqmmc_os_v1','sup_data')  ## supplmental data dir

## set woring dir
setwd(work_dir)

rvn_dir <- 'D:\\RavenPack'

##============================
##   ACTION CATEGORIES
##----------------------------
dtafile <- 'rp40_action_categories_2000_2018.dta'

##convert
rvn <- read_dta(file.path(rvn_dir, dtafile))




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
##  Compute Ravenpack ACTIONS per Firm-Year
##----------------------------------------


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
 `I(smmc1n^2)`='System MMC Squared',
 `pres1n`='Competitive Pressure',
 `feedback1`='System Feedback',
 `smmc1n:pres1n`='System MMC * Pressure',
 `smmc1n:feedback1`='System MMC * Feedback',
 `lag(feedback1, 0:0)`='System Feedback',
 `lag(feedback1, 0:0)0`='System Feedback',
 `lag(feedback1, 0:1)0`='System Feedback',
 `lag(feedback1, 0:1)1`='System Feedback (Lag 1)',
 `I(smmc1n^2):pres1n`='System MMC Squared * Pressure',
 `pres1n:I(smmc1n^2)`='System MMC Squared * Pressure',
 `I(smmc1n^2):lag(feedback1, 0:0)`='System MMC Squared * Feedback',
 `lag(feedback1, 0:0):I(smmc1n^2)`='System MMC Squared * Feedback',
 `I(smmc1n^2):feedback1`='System MMC Squared * Feedback',
 `feedback1:I(smmc1n^2)`='System MMC Squared * Feedback',
 `poly(smmc1n, 2)1` = 'System MMC',
 `poly(smmc1n, 2)2` = 'System MMC Squared',
 `poly(smmc1n, 2)1:feedback1` = 'System MMC * Feedback',
 `poly(smmc1n, 2)2:feedback1` = 'System MMC Squared * Feedback',
 `poly(smmc1n, 2)1:pres1n` = 'System MMC * Pressure',
 `poly(smmc1n, 2)2:pres1n` = 'System MMC Squared * Pressure'
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

##-----------------------------------
## POISSON - PGLM
##-----------------------------------
lm0 <- pglm(y ~ dum.crisis + I(acq_cnt_5 > 0) + 
               I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
               I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) 
             ,
             data=dfreg, family = poisson, 
            model = 'within', effect = 'twoways',
            R = 100, method='nr',
            index=c('i','year'))

lm1 <- pglm(y ~ dum.crisis + I(acq_cnt_5 > 0) + 
              I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
              I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
              smmc1n + I(smmc1n^2),
            data=dfreg, family = poisson, 
            model = 'within', effect = 'twoways',
            R = 100, method='nr',
            index=c('i','year'))

lm2 <- pglm(y ~  dum.crisis + I(acq_cnt_5 > 0) + 
              I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
              I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
              smmc1n + 
              pres1n +
              smmc1n:pres1n + 
              I(smmc1n^2) +
              I(smmc1n^2):pres1n,
            data=dfreg, family = poisson, 
            model = 'within', effect = 'twoways',
            R = 100, method='nr',
            index=c('i','year'))

lm3 <- pglm(y ~  dum.crisis + I(acq_cnt_5 > 0) + 
              I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
              I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
              smmc1n +
              feedback1 + 
              smmc1n:feedback1 + 
              I(smmc1n^2) +
              I(smmc1n^2):feedback1,
            data=dfreg, family = poisson, 
            model = 'within', effect = 'twoways',
            R = 100, method='nr',
            index=c('i','year'))

lmAll <- pglm(y ~  dum.crisis + I(acq_cnt_5 > 0) + 
                I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                smmc1n +
                pres1n +               
                feedback1 + 
                smmc1n:pres1n +                 
                smmc1n:feedback1 + 
                I(smmc1n^2) +
                I(smmc1n^2):pres1n +
                I(smmc1n^2):feedback1,
              data=dfreg, family = poisson, 
              model = 'within', effect = 'twoways',
              R = 100, method='nr',
              index=c('i','year'))


lmn.list <- list(lm0,lm1,lm2,lm3,lmAll)
# screenreg(lmn.list, digits = 3)
screenreg(lmn.list, digits = 3, 
          custom.coef.map = coef.map
          )

saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_pglm_poisson_list.rds'))
texreg::htmlreg(lmn.list,
                file.path(result_dir, 'acqmmc_pglm_poisson_list.html'),
                digits = 3, custom.coef.map = coef.map)




##=================================
##  GLMM ML / BOOT
##---------------------------------
gmAll <- glmmML(y ~  dum.crisis + I(acq_cnt_5 > 0) + 
                  I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                  I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                  poly(smmc1n, 2) + 
                  pres1n +
                  feedback1 + # I(presconstr^2) +
                  poly(smmc1n, 2):pres1n +
                  poly(smmc1n, 2):feedback1,
                data=dfreg, family = poisson, 
                cluster= dfreg$i, 
                boot=100)
gmbAll <- glmmboot(y ~  dum.crisis + I(acq_cnt_5 > 0) + 
                   I(acq_sum_1/1e9)  + I(employee_na_age/1e3) + 
                   I(sales_na_0_mn/1e6) + log(1 + cent_deg_all) +
                   poly(smmc1n, 2) + 
                   pres1n +
                   feedback1 + # I(presconstr^2) +
                   poly(smmc1n, 2):pres1n +
                   poly(smmc1n, 2):feedback1,
                 data=dfreg, family = poisson, 
                 cluster= dfreg$i, 
                 boot=500)
round(summary(gmbAll))

##=======================================
##  GLMER
##---------------------------------------

##-----------------------------------
## NEGBIN
##-----------------------------------
lmn0 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
               dum.crisis + I(acq_cnt_5>0) + 
               acq_sum_1_sc  + employee_na_age_sc + 
               sales_na_0_mn_sc + cent_deg_all_sc
               ,
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=1e5)))

lmn1 <- glmer.nb(y ~ (1 | i) + (1 | year) + 
                dum.crisis + I(acq_cnt_5>0) + 
                  acq_sum_1_sc  + employee_na_age_sc + 
                  sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2),
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=3e5)))
##***********convergence warning: 0.00135935***********
lmn2 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                dum.crisis + I(acq_cnt_5>0) + 
                  acq_sum_1_sc  + employee_na_age_sc + 
                  sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2) + 
               pres1n +
               I(smmc1n^2):pres1n,
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=4e5)))
##***********convergence warning: 0.0611716***********
lmn3 <- glmer.nb(y ~ (1 | i) + (1 | year) +  
                dum.crisis + I(acq_cnt_5>0) + 
                  acq_sum_1_sc  + employee_na_age_sc + 
                  sales_na_0_mn_sc + cent_deg_all_sc +
               smmc1n + I(smmc1n^2) +
               lag(feedback1, 0:0) + 
               I(smmc1n^2):lag(feedback1, 0:0),
             data=dfreg, 
             control=glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=4e5)))

lmnAll <- glmer.nb(y ~ (1 | i) + (1 | year) + 
                  dum.crisis + I(acq_cnt_5>0) + 
                    acq_sum_1_sc  + employee_na_age_sc + 
                    sales_na_0_mn_sc + cent_deg_all_sc +
                 smmc1n + I(smmc1n^2) + 
                 pres1n +
                 lag(feedback1, 0:0) + # I(presconstr^2) +
                 I(smmc1n^2):pres1n +
                 I(smmc1n^2):lag(feedback1, 0:0),
               data=dfreg, 
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=4e5)))

lmn.list <- list(lmn0,lmn1,lmn2,lmn3,lmnAll)
screenreg(lmn.list, digits = 3, 
          custom.coef.map = coef.map)

saveRDS(lmn.list, file = file.path(result_dir, 'acqmmc_glmernb_list.rds'))
texreg::htmlreg(lmn.list,
                file.path(result_dir, 'acqmmc_glmernb_list.html'),
                digits = 3, custom.coef.map = coef.map)






# pglm(formula, data, subset, na.action,
#      effect = c("individual", "time", "twoways"),
#      model = c("random", "pooling", "within", "between"),
#      family, other = NULL, index = NULL, start = NULL, R = 20,  ...) 

# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#            dfreg,
#            family = negbin, effect = 'twoways', model = "within", 
#            print.level = 3, method = "nr", ## 'bfgs'
#            index = c('i', 'year'), R=20)
# 
# m2p <- pglm(y.new ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "pooling", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4p <- pglm(y.new ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = negbin, effect = 'twoways', model = "pooling",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# summary(m2w); summary(m2p); summary(m4w); summary(m4p)
# 
# 
# 
# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2b <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "between", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4b <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'time', model = "between",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2w); summary(m2b); summary(m4w); summary(m4b)
# 
# 
# 
# m2w <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2b <- pglm(y.new ~ smmc2 + lag(feedback2, 1) + smmc2:lag(feedback2, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "between", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4w <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# m4b <- pglm(y.new ~ smmc4 + lag(feedback4, 1) + smmc4:lag(feedback4, 1),
#             dfreg,
#             family = negbin, effect = 'individual', model = "between",
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2w); summary(m2b); summary(m4w); summary(m4b)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# m2x <- pglm(y.new ~ smmc2 + feedback2 + smmc2:feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "random", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2xc <- pglm(y.cur ~ smmc2 + I(smmc2^2) + feedback2 + I(smmc2^2) :feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "random", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# summary(m2x); summary(m2xc)
# 
# 
# 
# 
# m2x <- pglm(y.new ~ smmc2 + feedback2 + smmc2:feedback2,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m2xc <- pglm(y.cur ~ smmc2 + I(smmc2^2) + feedback2 + I(smmc2^2) :feedback2,
#              dfreg,
#              family = negbin, effect = 'twoways', model = "within", 
#              print.level = 3, method = "nr", ## 'bfgs'
#              index = c('i', 'year'), R=20)
# summary(m2x); summary(m2xc)
# 
# 
# 
# m4x <- pglm(y.new ~ smmc4 + feedback4 + smmc4:feedback4,
#             dfreg,
#             family = negbin, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20)
# 
# m4xc <- pglm(y.cur ~ smmc4 + I(smmc4^2) + feedback4 + I(smmc4^2) :feedback4,
#              dfreg,
#              family = negbin, effect = 'twoways', model = "within", 
#              print.level = 3, method = "nr", ## 'bfgs'
#              index = c('i', 'year'), R=20)
# summary(m4x); summary(m4xc)
# 
# 
# 
# 
# 
# m2yl <- pglm(y ~ smmc2 + lag(feedback2, 0:1) + smmc2:lag(feedback2, 0:1),
#             dfreg,
#             family = poisson, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20); summary(m2yl)
# 
# m4yl <- pglm(y ~ smmc4 + lag(feedback4, 0:1) + smmc4:lag(feedback4, 0:1),
#             dfreg,
#             family = poisson, effect = 'twoways', model = "within", 
#             print.level = 3, method = "nr", ## 'bfgs'
#             index = c('i', 'year'), R=20); summary(m4yl)

####################################################################
m4y <- pglm(y ~ smmc4 + lag(feedback4, 0) + smmc4:lag(feedback4, 0),
             dfreg,
             family = poisson, effect = 'twoways', model = "within", 
             print.level = 3, method = "nr", ## 'bfgs'
             index = c('i', 'year'), R=20); summary(m4y)
####################################################################


m4y2 <- pglm(y ~ smmc4 + I(smmc4^2) + lag(feedback4, 0) + 
               smmc4:lag(feedback4, 0) + I(smmc4^2):lag(feedback4, 0),
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4y2)

m4yf <- pglm(y ~ smmc4 +  lag(feedback4, 0:1) + lag(I(feedback4^2), 0:1) +
               smmc4:lag(I(feedback4^2), 0:1),
             dfreg,
             family = poisson, effect = 'twoways', model = "within", 
             print.level = 3, method = "nr", ## 'bfgs'
             index = c('i', 'year'), R=20); summary(m4yf)






############## CONTROLS ##################################
m4a <- pglm(y~ smmc1 + pres1 + smmc1:pres1 + 
              acq_cnt_5 +  cent_deg_all,
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4a)
############## //CONTROLS ##################################

m4a <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              # cent_deg_all +
              smmc1n + # I(smmc1n^2) +
              pres1n + # I(presconstr^2) +
              smmc1n:pres1n              # I(smmc1n^2):I(presconstr^2) +
              ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4a)



m4b <- pglm(y ~  acq_cnt_5 + I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n + # I(smmc1n^2) +
              pres1n + 
              lag(feedback1, 0:1) + # I(presconstr^2) +
              smmc1n:lag(feedback1, 0:1)
            # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg,
            family = poisson, effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4b)


m4c <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n +  I(smmc1n^2) +
              pres1n +  I(pres1n^2) +
              smmc1n:pres1n +
              I(smmc1n^2):I(pres1n^2) 
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4c)




m4d <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              cent_deg_all +
              smmc1n +  I(smmc1n^2) +
              lag(feedback1, 0:1) + 
              I(smmc1n^2):lag(feedback1, 0:1)
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "within", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4d)


dfreg$smmc1n.pres1n <- dfreg$smmc1n-dfreg$pres1n
m4e <- pglm(y.new ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + 
              pres1n + # I(presconstr^2) +
              smmc1n:pres1n # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m4e)

m5e <- pglm(y ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) +
              pres1n + 
              I(smmc1n^2):pres1n 
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m5e)

m6e <- pglm(y.new ~  I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + 
              lag(feedback1, 0:1) + # I(presconstr^2) +
              smmc1n:lag(feedback1, 0:1) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m6e)

dfreg$dum.crisis <- 0
dfreg$dum.crisis[which(dfreg$year %in% as.character(2009:2017))] <- 1
m7e <- pglm(y ~  dum.crisis + I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) + 
              pres1n +
              # lag(feedback1, 0:0) + # I(presconstr^2) +
              I(smmc1n^2):pres1n #+
              # I(smmc1n^2):lag(feedback1, 0:0) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m7e)

m8e <- pglm(y ~  dum.crisis + I(acq_cnt_5>0) + 
              I(acq_sum_1/1e9)  + 
              I(employee_na_age/1e6) + 
              I(sales_na_0_mn/1e6) +
              log(1+cent_deg_all) +
              smmc1n + I(smmc1n^2) + 
              pres1n +
              lag(feedback1, 0:0) + # I(presconstr^2) +
              I(smmc1n^2):pres1n +
              I(smmc1n^2):lag(feedback1, 0:0) # I(smmc1n^2):I(presconstr^2) +
            ,
            dfreg, family = poisson, 
            effect = 'twoways', model = "random", 
            print.level = 3, method = "nr", ## 'bfgs'
            index = c('i', 'year'), R=20); summary(m8e)


mx <- m4e
y <- mx$model[,1]
X <- mx$model[,-1]
X$`I(acq_cnt_5 > 0)` <- as.numeric(X$`I(acq_cnt_5 > 0)`)
X$`smmc1.pres1n` <- X$smmc1n * X$pres1n
beta <- mx$estimate

m4e <- as.matrix(X) %*% beta

y <- mx$model$y

##========================================
##  LINEAR
##----------------------------------------

qs.smmc1n <- quantile(X$smmc1n,  c(1/3,.5,2/3))
smmc1n.f <- unname(sapply(X$smmc1n, function(x){
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))
#
qs.pres1n <- quantile(X$pres1n, c(1/3,.5,2/3))
pres1n.f <- unname(sapply(X$pres1n, function(x){
  if(x >= qs.pres1n[2]) return('2High')
  return('1Low')
}))

xdf <- data.frame(y=y, smmc1n.f=smmc1n.f, 
                  pres1n.f=pres1n.f,
                  year=as.numeric(as.character(dfreg$year)))
# xdf <- xdf[which(xdf$year < 2012),]

interaction.plot(xdf$smmc1n.f, xdf$pres1n.f, xdf$y,
                 fun = function(x)mean(x, na.rm=T),
                 xlab='System MMC',
                 ylab='Avg. Acquisitions',
                 trace.label = 'Pressure',
                 type = 'b', pch=c(2,16,8))

##========================================
##  Curvilinear
##----------------------------------------

qs.smmc1n <- quantile(X$smmc1n,  c(1/3,.5,2/3))
smmc1n.f <- unname(sapply(X$smmc1n, function(x){
  if(x >= qs.smmc1n[3]) return('3High')
  if(x >= qs.smmc1n[1]) return('2Mid')
  return('1Low')
  }))
#
qs.pres1n <- quantile(X$pres1n, c(1/3,.5,2/3))
pres1n.f <- unname(sapply(X$pres1n, function(x){
    # if(x >= qs.pres1n[3]) return('3High')
    # if(x >= qs.pres1n[1]) return('2Mid')
    # return('1Low')
    if(x >= qs.pres1n[2]) return('2High')
    return('1Low')
  }))

interaction.plot(smmc1n.f, pres1n.f, y,
                 fun = function(x)mean(x, na.rm=T),
                 xlab='System MMC',
                 ylab='Avg. Acquisitions',
                 trace.label = 'Pressure',
                 type = 'b', pch=c(2,16,8))






# ##--------------------------------------------
# ## Patch category similarity
# ##--------------------------------------------
# patch.todo <- c('qualtrics',
#                 'clarabridge','confirmit','medallia','snap-surveys-ltd')
# for (name_i in patch.todo) {
#   cat(sprintf(' patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %n% 'cat_cos_sim' <- acf$.cov.categoryCosineSimilarity(nets[[t]])
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }

##--------------------------------------------
## Patch  sizes variables
##--------------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
for (name_i in patch.todo) {
  cat(sprintf('\n\n patching %s\n', name_i))
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  end <- as.numeric(names(nets)[t])
  start <- end-1
  for (t in 1:length(nets)) {
    nets[[t]] <- acf$setCovariates(nets[[t]], start, end,
                                  covlist=c('employee','sales'), size=si)
  }
  saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
}

sapply(lapply(nets,function(net)net%v%'employee_na_age'), summary)
sapply(lapply(nets,function(net)net%v%'employee_na_catage'), summary)
sapply(lapply(nets,function(net)net%v%'sales_na_0_mn'), summary)


###############
# for (name_i in patch.todo) {
#   cat(sprintf('\n\n patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %v% 'employee_na_0_log' <- log(1 + (nets[[t]] %v% 'employee_na_0'))
#     nets[[t]] %v% 'employee_na_cat_log' <- log(1 + (nets[[t]] %v% 'employee_na_cat'))
#     nets[[t]] %v% 'employee_na_age_log' <- log(1 + (nets[[t]] %v% 'employee_na_age'))
#     nets[[t]] %v% 'employee_na_catage_log' <- log(1 + (nets[[t]] %v% 'employee_na_catage'))
#     nets[[t]] %v% 'employee_na_0_std' <- c(scale(nets[[t]] %v% 'employee_na_0'))
#     nets[[t]] %v% 'employee_na_cat_std' <- c(scale(nets[[t]] %v% 'employee_na_cat'))
#     nets[[t]] %v% 'employee_na_age_std' <-  c(scale(nets[[t]] %v% 'employee_na_age'))
#     nets[[t]] %v% 'employee_na_catage_std' <-  c(scale(nets[[t]] %v% 'employee_na_catage'))
#     #
#     nets[[t]] %v% 'sales_na_0_log' <- log(1 + (nets[[t]] %v% 'sales_na_0'))
#     nets[[t]] %v% 'sales_na_cat_log' <- log(1 + (nets[[t]] %v% 'sales_na_cat'))
#     nets[[t]] %v% 'sales_na_age_log' <- log(1 + (nets[[t]] %v% 'sales_na_age'))
#     nets[[t]] %v% 'sales_na_catage_log' <- log(1 + (nets[[t]] %v% 'sales_na_catage'))
#     nets[[t]] %v% 'sales_na_0_std' <- c(scale(nets[[t]] %v% 'sales_na_0'))
#     nets[[t]] %v% 'sales_na_cat_std' <- c(scale(nets[[t]] %v% 'sales_na_cat'))
#     nets[[t]] %v% 'sales_na_age_std' <-  c(scale(nets[[t]] %v% 'sales_na_age'))
#     nets[[t]] %v% 'sales_na_catage_std' <-  c(scale(nets[[t]] %v% 'sales_na_catage'))
# 
#     nets[[t]] %v% 'sales_na_0_mn' <- (nets[[t]] %v% 'sales_na_0') / 1e6
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }






##--------------------------------------------
## Patch Centrality
##--------------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
for (name_i in patch.todo) {
  cat(sprintf('\n\n patching %s\n', name_i))
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  for (t in 1:length(nets)) {
    nets[[t]] <- acf$.cov.centrality(nets[[t]])
    year <- as.numeric(names(nets)[t])
    nets[[t]] %n% 'shared_investor_nd' <- acf$.cov.sharedInvestor(nets[[t]], ih, cb$co_rou, cb$inv_rou, cb$inv, year, off.diagonal.blocks=F)
  }
  saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
}


##-----------------------------------------
## Check encounters distribution per year
##-----------------------------------------
patch.todo <- c('qualtrics',
                'clarabridge','confirmit','medallia','snap-surveys-ltd','getfeedback')
ckl <- list()
for (name_i in patch.todo) {
  nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
  ckl[[name_i]] <- sapply(nets, function(net){
      g <- asIgraph(net)
      g.sub <- igraph::induced.subgraph(g,which(igraph::degree(g)>0))
      deg <- igraph::degree(g.sub)
      return(c(e=ecount(g.sub),v=vcount(g.sub),
               deg_avg=round(mean(deg),1), deg_med=median(deg),
               deg_min=min(deg), deg_max=max(deg)))
    })
  # ec <- sapply(nets, function(net) {
  #     m <- net[,]
  #     return(sum(m[lower.tri(m)])) 
  #   })
  # vc <- sapply(nets, function(net) nrow(net[,]))
  # ckdf <- rbind(ckdf, data.frame(
  #     name=name_i,
  #     yrs=length(nets),
  #     v_sum=sum(vc), 
  #     v_avg=mean(vc),
  #     v_min=min(vc),
  #     v_max=max(vc),
  #     e_sum=sum(ec),
  #     e_avg=mean(ec),
  #     e_min=min(ec),
  #     e_max=max(ec)
  #   ))
}

ckrm <- lapply(ckl, rowMeans)
ckrs <- lapply(ckl, rowSums)
ckrmdf <- t(round(as.data.frame(ckrm),1))
ckrsdf <- t(round(as.data.frame(ckrs),1))

colSums(ckrmdf)

colSums(ckrsdf)

##-----------------------------------------------
## TEST BTERGM
##-----------------------------------------------
# name_i <- 'qualtrics'
# nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
# acf$covSummaryPlot(nets, name_i, net_dir)
library(btergm)
name_i <- 'qualtrics'
d <- 3
nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))

mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor'))
shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))

####################### DEFINE MODELS ###################################

m4 <-   nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed=T) +
  nodematch("ipo_status", diff = F) +
  nodematch("state_code", diff = F) +
  nodecov("age") + absdiff("age") +
  nodecov("employee_na_age") + nodecov("sales_na_0") +
  edgecov(cossim) + edgecov(centjoin) + edgecov(centratio) + edgecov(shcomp) + edgecov(shinv) +
  edgecov(mmc) +
  ##edgecov(cpa) +
  ##edgecov(cpc) +
  ##edgecov(cpp) +
  memory(type = "stability", lag = 1) +
  timecov(transform = function(t) t) +
  nodecov("genidx_multilevel") +
  nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") +
  cycle(3) + cycle(4) + cycle(5)




# ##--------------------------------------------
# ## Patch category similarity 
# ##--------------------------------------------
# patch.todo <- c('qualtrics',
#                 'clarabridge','confirmit','medallia','snap-surveys-ltd')
# for (name_i in patch.todo) {
#   cat(sprintf(' patching %s\n', name_i))
#   nets <- readRDS(file.path(net_dir, sprintf('%s_d%d.rds',name_i,d)))
#   for (t in 1:length(nets)) {
#     nets[[t]] %n% 'cat_cos_sim' <- acf$.cov.categoryCosineSimilarity(nets[[t]])
#   }
#   saveRDS(nets, file = file.path(net_dir, sprintf('%s_d%d_patched.rds',name_i,d)))
# }


# ### Network size by scope of awareness
# name_i <- 'qualtrics'
# ds <- 1:10
# names(ds) <- as.character(ds)
# gstat <- sapply(ds, function(d){
#   g.base <- g.full  
#   g.d.sub <- igraph::make_ego_graph(g.base, d, V(g.base)[V(g.base)$name==name_i], 'all')[[1]]
#   return(c(v=vcount(g.d.sub),e=ecount(g.d.sub),d=igraph::graph.density(g.d.sub)))
# })
# 
# plot(gstat[1,], log='', type='b')
# plot(c(NA,diff(gstat[1,])), log='', type='b')





## get focal firm distances to specific firms over time
net.yrs <- 2007:2016
gs <- list()
for (yr in net.yrs) {
  gs[[as.character(yr)]] <- asIgraph(readRDS(file.path('firm_nets_rnr',sprintf('facebook_d3_y%d.rds',yr))))
}


g = gs[[length(gs)]]

co.idx <- which(
  grepl('content and publishing',V(g.full)$category_group_list) 
  & !(V(g.full)$employee_count %in% c('NA','-','1-10','11-50','51-100'))
)
co.names <- c(V(g.full)$name[co.idx],'amazon','aws')
ds = sapply(gs[2:10],function(g)igraph::distances(g, v = which(V(g)$vertex.names=='facebook'),
                                            to = V(g)[V(g)$vertex.names %in% co.names]))
# ds2 = t(igraph::distances(gs[[11]], 
#                         v = which(V(gs[[11]])$vertex.names=='facebook'),
#                         to = which(V(gs[[11]])$vertex.names %in% co.names)))
# ds = cbind(ds1,ds2)
# colnames(ds) <- names(gs)
rownames(ds) <- V(gs[[2]])$vertex.names[V(gs[[2]])$vertex.names %in% co.names]

row.names(ds) = V(g)$vertex.names[V(g)$vertex.names %in% co.names]
print(ds)

write.csv(ds, "facebook_news_dist_2008-2016_d3.csv")

igraph::distances(g.full, V(g.full)[V(g.full)$name=='facebook'],V(g.full)[V(g.full)$name=='hearstcorporation'])

## get cb category coverage
cgl = unique(c(str_split(V(g)$category_group_list,'[|]',simplify = T)))
cl = unique(c(str_split(V(g)$category_list,'[|]',simplify = T)))


# load('netrisk_dynamic_firm_nets.RData')

# # ADD CONSTRAINT NODE PROPERTY
# for (t in 1:length(nets)) {
#   g.tmp <- getIgraphFromNet(nets[[t]])
#   if (vcount(g.tmp)>0 & ecount(g.tmp)>0) {
#     cons <-  igraph::constraint(g.tmp)
#     cons[is.nan(cons) | is.na(cons)] <- 0 ### ???
#     nets[[t]] %v% 'constraint' <- cons
#   }
# }

#------------------------------------------------------
#              Predictors Diagnostics
# #------------------------------------------------------
# ## Plot density
# n <- ceiling(sqrt(length(firm.nets)))
# m <- ifelse(n*(n-1) >= length(firm.nets), n-1, n)
# par(mfrow=c(m,n), mar=c(2.5,2.5,2,1))
# for (firm_i in names(firm.nets)) {
#   nets <- firm.nets[[firm_i]]
#   plot(as.numeric(names(nets))-1,   
#        sapply(nets,function(net) {
#          sum(net[lower.tri(net)])/(nrow(net[,])*nrow(net[,]-1)/2) 
#        }), 
#        ylab='density', xlab='year', type='b', main=firm_i)
# }
# 
# ## Net Risk
# par(mfrow=c(3,3), mar=c(2.5,2.5,2,1))
# for (firm_i in names(firm.nets)) {
#   nets <- firm.nets[[firm_i]]
#   sapply(seq_along(nets), function(j) {
#     hist(nets[[j]] %v% 'net_risk', breaks=25, main=sprintf('%s %s',firm_i,names(nets)[j]))
#   })
# }

#------------------------------------------------------
