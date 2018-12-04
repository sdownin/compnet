##--------------------------------
## NETWORK MMC Acquisition selection rEGRESSION
#
## @see http://www.elff.eu/software/#mclogit
# The package 'mclogit' implements the estimation of mixed conditional logit models via the PQL method as used in my article published in Electoral Studies. It is published on CRAN. Development occurs on GitHub, where both releases and the development tree can be found.
# 
# The probability that individual i chooses alternative j from choice set Si is
# ??ij=expa(??ij)???k???Siexpa(??ik)
# with
# ??ij=??1x1ij+???+??qxqij+Uij
# where xhij are values of independent variables, ??h are parameters (coefficients), and Uij are random effects with a normal distribution.
# 
# The package allows to specify that random effects are equal for all individuals within clusters Cc, that is Ui1j=Ui2j for i1???Cc and i2???Cc, where such clusters also may be nested in a "multi-level" manner.
# 
# The "dependent variable" yij may be a "dummy variable" that is equal to 1 if individual i has chosen alternative j and 0 if s/he has chosen another alternative. For example, if all individuals i face the same set of five alternatives, then five values of the dependent variable would correspond to each individual with only one of the values being equal to one and the other four values being equal to zero. (This is sometimes called that the data are in "stacked" format.)
# 
# Also, if "covariate classes" of individuals are formed that share the same values of the independent variables and are members of the same cluster, and i indicates such a covariate class, yij may be the count of individuals from covariate class i that have chosen alternative j.
##-------------------------------

setwd("C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2")
# .libPaths('C:/Users/T430/Documents/R/win-library/3.2')
library(texreg, quietly = T)
library(memisc, quietly = T)
library(plyr, quietly = T)
library(lattice, quietly = T)
library(latticeExtra, quietly = T)
library(ggplot2, quietly = T)
library(reshape2)
library(mlogit)
library(mclogit)
library(lme4)
library(lmerTest)
library(lubridate)

data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/crunchbase/"
# syn_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/acqlogit_tmp_20181123"
acq_data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2/acqlogit_data"

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##============================================
## load data 
##   **SLOW TO LOAD**
##--------------------------------------------
## 1. CrunchBase data
cb  <- source(file.path(getwd(),'R','acqlogit','cb_data_prep.R'))$value      ## DATA 


##=================================
## SETTINGS:
## focal firm's competition network
##---------------------------------
name_i <- 'ibm'
d <- 2

## LOAD DATALIST
# data.in <- readRDS(sprintf("acqlogit_data/acqlogit_compnet_processed_acquisitions_synergies_list_%s_d%s.rds",name_i,d))
# data.in <- readRDS(sprintf("%s/acqlogit_compnet_processed_acquisitions_synergies_list_%s_d%s.rds",acq_data_dir,name_i,d))
data.in <- readRDS(sprintf("%s/acqlogit_compnet_processed_acquisitions_synergies_ABSORB_list_%s_d%s.rds",acq_data_dir,name_i,d))
l <- data.in$l
l.cov <- data.in$l.cov

df.reg <- ldply(l.cov, function(x) {  as.data.frame(x) })
# df.reg <- ldply(l.cov, function(x) {
#   names(x)[which(names(x)=="")] <- c('j.fund.v.cnt','j.fund.v.amt','j.fund.cnt','j.fund.amt')
#   as.data.frame(x)
# })
print(length(unique(df.reg$uuid)))
# View(df.reg)


## COPY (or subset) REGRESSION DATAFRAME
df.sub <- df.reg

## ADD MISSING COVARS
df.sub$ij.age.diff <- df.sub$i.age - df.sub$j.age

## FIX INFINITE DISTANCES
dists.nonInf <- df.sub$ij.dist[df.sub$ij.dist < Inf]
df.sub$ij.dist[df.sub$ij.dist == Inf] <- round(sd(dists.nonInf)) + max(dists.nonInf) ## 1 SD above max

## scale closeness
df.sub$ij.syn.closeness2 <- df.sub$ij.syn.closeness * 1e4

## scale acquisition experience ratio
df.sub$i.acqs.100 <- df.sub$i.acqs * 100
df.sub$j.acqs.100 <- df.sub$j.acqs * 100

## COUNTERFACTUAL ACQUISITION PAIRING SIMILARITY
df.sub$ij.sim <- NA
df.sub$ij.cossim <- NA
for (i in 1:nrow(df.sub)) {
  firm.i <- as.character(df.sub$i[i])
  firm.j <- as.character(df.sub$j[i])
  cats.i <- str_split(cb$co$category_list[cb$co$company_name_unique == firm.i], '[|]')[[1]]
  cg.i <- str_split(cb$co$category_group_list[cb$co$company_name_unique == firm.i], '[|]')[[1]]
  c.i <- unique(c(cats.i,cg.i))
  cats.j <- str_split(cb$co$category_list[cb$co$company_name_unique == firm.j], '[|]')[[1]]
  cg.j <- str_split(cb$co$category_group_list[cb$co$company_name_unique == firm.j], '[|]')[[1]]
  c.j <- unique(c(cats.j,cg.j))
  c.all <- unique(c(c.i,c.j))
  v.i <- as.integer(c.all %in% c.i)
  v.j <- as.integer(c.all %in% c.j)
  df.sub$ij.sim[i] <- (v.i %*% v.j)[1,1]
  df.sub$ij.cossim[i] <- (v.i %*% v.j)[1,1] / (sqrt(sum(v.i^2)) * sqrt(sum(v.j^2)))
  if (i %% 50 == 0) cat(sprintf('%s\n',i))
}
write.csv(df.sub, file = sprintf("%s/acqlogit_compnet_processed_acquisitions_synergies_df_COSSIM_%s_d%s.rds",tmp_dir,name_i,d))
## CACHE SUBSETTED DATAFRAME AFTER ADDING SIMILARITIES
df.sub.orig <- df.sub


# ## PATCH GLOBAL DEGREE OF TARGET
# if (length(unique(df.sub$t))==length(l)) {
#   df.sub$j.deg.full <- NA
#   for (i in 1:length(l)) {
#     cat(sprintf('patching acquisition %s (%.2f%s)\n',i,100*i/length(l),'%'))
#     tmp <- l[[i]]$df.targ.alt
#     t <- unique(df.sub$t)[i]
#     for (row in 1:nrow(tmp)) {
#       df.sub$j.deg.full[which(df.sub$t==t & df.sub$j==tmp$company_name_unique[row])] <- tmp$deg.full[row]
#     }
#   }
# } else {
#   print("cannot patch: length of list l != number of acquisition in df.sub")
# }
# ## PATCH GLOBAL DEGREE OF ACQUIRER
# if (length(unique(df.sub$t))==length(l)) {
#   df.sub$i.deg.full <- NA
#   for (i in 1:length(l)) {
#     cat(sprintf('patching acquisition %s (%.2f%s)\n',i,100*i/length(l),'%'))
#     tmp <- l[[i]]$df.acq.alt
#     t <- unique(df.sub$t)[i]
#     for (row in 1:nrow(tmp)) {
#       df.sub$i.deg.full[which(df.sub$t==t & df.sub$i==tmp$company_name_unique[row])] <- tmp$deg.full[row]
#     }
#   }
# } else {
#   print("cannot patch: length of list l != number of acquisition in df.sub")
# }

# ## SUBSET PAIRINGS WITH USABLE COMPUSTAT CONTROLS
check.cols <- c('i.ln_asset','i.cash','i.roa','i.m2b',
                'j.age','i.age','i.deg.full','j.deg.full','ij.same.region',
                'i.fm.mmc.sum','i.acqs','j.acqs','ij.cossim','ij.dist',
                'j.constraint','i.constraint',
                'ij.syn.constraint','ij.syn.closeness',
                'i.div')
for (col in check.cols) {
  df.sub <- df.sub[which(!is.na(df.sub[,col]) & !is.nan(df.sub[,col])), ]
}

##==================================================
##  ORGSCI 2018-11-08
##--------------------------------------------------

mg0 <- mclogit(
  cbind(y,t) ~  
    i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + 
    # I(j.age^2) + 
    i.age +
    i.deg.full + j.deg.full +
    ij.same.region +
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim +
    ij.dist +
    j.constraint +
    i.constraint + 
    i.pow.n3 +
    ij.syn.constraint +
    # ij.syn.degree +
    # ij.syn.closeness2 +
    ij.syn.pow.n3
  ,
  data = df.sub)
summary(mg0)

mg1 <- mclogit(
  cbind(y,t) ~  
    i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age +  i.age + 
    # I(j.age^2) +
    i.deg.full + j.deg.full +
    ij.same.region +  
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim + 
    ij.dist +
    j.constraint +
    i.constraint + i.pow.n3 +
    ij.syn.constraint +
    # ij.syn.degree +
    # ij.syn.closeness2 +
    ij.syn.pow.n3 +
    I(i.fm.mmc.sum^2)
  ,
  data = df.sub)
summary(mg1)

mg2 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + 
    # I(j.age^2) + 
    i.age +
    i.deg.full + j.deg.full +
    ij.same.region +  
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim + 
    ij.dist +
    j.constraint +
    i.constraint + i.pow.n3 +
    ij.syn.constraint +
    # ij.syn.degree +
    # ij.syn.closeness2 +
    ij.syn.pow.n3 + 
    I(i.fm.mmc.sum^2) + 
    I(i.fm.mmc.sum^2):ij.syn.constraint
  ,
  data = df.sub)
summary(mg2)


# 
# mg3 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
#     j.age + I(j.age^2) +  i.age +
#     i.deg.full + j.deg.full +
#     ij.same.region +  
#     i.fm.mmc.sum + i.acqs + j.acqs +
#     ij.cossim + 
#     ij.dist +
#     j.constraint + i.constraint + i.pow.n3 +
#     ij.syn.constraint +
#     # ij.syn.degree +
#     ij.syn.closeness2 + 
#     ij.syn.pow.n1 + 
#     I(i.fm.mmc.sum^2) + 
#     I(i.fm.mmc.sum^2):ij.syn.degree
#   ,
#   data = df.sub)
# summary(mg3)

mg4 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + 
    # I(j.age^2) +  
    i.age +
    i.deg.full + j.deg.full +
    ij.same.region +  
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim + 
    ij.dist +
    j.constraint +
    i.constraint + i.pow.n3 +
    ij.syn.constraint +
    # ij.syn.degree +
    # ij.syn.closeness2 +
    ij.syn.pow.n3 + 
    I(i.fm.mmc.sum^2) + 
    I(i.fm.mmc.sum^2):ij.syn.pow.n3
  ,
  data = df.sub)
summary(mg4)

mg5 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + 
    # I(j.age^2) +  
    i.age +
    i.deg.full + j.deg.full +
    ij.same.region +  
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim + 
    ij.dist +
    j.constraint +
    i.constraint + i.pow.n3 +
    ij.syn.constraint +
    # ij.syn.degree +
    # ij.syn.closeness2 +
    ij.syn.pow.n3 + 
    I(i.fm.mmc.sum^2) + 
    I(i.fm.mmc.sum^2):ij.syn.constraint + 
    # I(i.fm.mmc.sum^2):ij.syn.closeness2 
    # I(i.fm.mmc.sum^2):ij.syn.degree + 
    I(i.fm.mmc.sum^2):ij.syn.pow.n3
  ,
  data = df.sub)
summary(mg5)

mg6 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + 
    # I(j.age^2) +  
    i.age +
    i.deg.full + j.deg.full +
    ij.same.region +
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim +
    ij.dist +
    j.constraint +
    i.constraint + i.pow.n3 +
    ij.syn.constraint +
    # ij.syn.degree +
    ij.syn.closeness2 +
    ij.syn.pow.n3 +
    I(i.fm.mmc.sum^2) +
    I(i.fm.mmc.sum^2):ij.syn.closeness2
  ,
  data = df.sub)
summary(mg6)

tabmost <- mtable(mg0,mg1,mg2,mg4,mg5,mg6)
print(tabmost)
## SAVE FILE
memisc::write.mtable(tabmost, 
        file = sprintf("acqlogit_data/acqlogit_reg_table_ABSORB_%s_d%s_OrgSci20181111.tsv",name_i,d))


## Chisq Deviance Test
anova(mg0,mg1, test='Chisq')
anova(mg0,mg2, test='Chisq')
anova(mg0,mg3, test='Chisq')
anova(mg0,mg4, test='Chisq')
anova(mg0,mg5, test='Chisq')
anova(mg0,mg6, test='Chisq')

## Correlations and descriptives stats
library(psych)

X <- mg4$model
X <- cbind(y=X[,1][,1],t=X[,1][,2],X[,-1])  ## unbind the cbind(y,t)

# X$react.rival.2_mmc2 <- X$react.rival.2 * X$`I(i.fm.mmc.sum^2)`
# X$ij.syn.closeness_mmc2 <- X$ij.syn.closeness * X$`I(i.fm.mmc.sum^2)`
X$ij.syn.constraint_mmc2 <- X$ij.syn.constraint * X$`I(i.fm.mmc.sum^2)`
X$ij.syn.degree_mmc2 <- X$ij.syn.degree * X$`I(i.fm.mmc.sum^2)`

##----------------------------
## Rename
##----------------------------
namemap <- list(
  y = 'y',
  t = 't',
  i.ln_asset = 'Ln Assets (Acquirer)',            
  i.cash = 'Cash Holdings (Acquirer)',          
  i.roa = 'ROA (Acquirer)',          
  i.ln_emp = 'Ln Employees (Acquirer)',  
  i.m2b = 'Market-to-Book (Acquirer)', 
  i.div = 'Diversification (Acquirer)', 
  j.age = 'Age (Target)',             
  i.age = 'Age (Acquirer)', 
  `I(j.age^2)` = 'Age^2 (Target)',
  ij.age.diff = 'Age Diff.',         
  i.deg = 'Competitors (Acquirer)', ## from focal firm's ego network
  i.deg.full = 'Competitors (Acquirer)', ## from global network
  j.deg = 'Competitors (Target)',   ## from focal firm's ego network
  j.deg.full = 'Competitors (Target)',  ## from global network
  ij.same.region = 'Geographic Region Homophily' ,    
  i.fm.mmc.sum = 'FM-MMC',        
  i.acqs = 'Acquisition Experience (Acquirer)',
  j.acqs = 'Acquisition Experience (Target)',
  ij.cossim = 'Product Dissimilarity',  
  ij.dist = 'Competitive Distance',        
  react.rival.2 = 'Rivals Acquisitions (Acquirer)',  
  j.constraint = 'Constraint (Target)',
  i.constraint = 'Constraint (Acquirer)', 
  ij.syn.constraint = 'Structural Synergy',
  ij.syn.degree = 'Positional Synergy (degree)',
  ij.syn.closeness = 'Positional Synergy (closeness)',        
  `I(i.fm.mmc.sum^2)` = 'FM-MMC^2 (Acquirer)',       
  react.rival.2_mmc2 = 'FM-MMC^2 * Rivals Acquisitions',      
  ij.syn.closeness_mmc2 = 'FM-MMC^2 * Positional Synergy (closeness)',  
  ij.syn.degree_mmc2 = 'FM-MMC^2 * Positional Synergy (degree)',  
  ij.syn.constraint_mmc2 = 'FM-MMC^2 * Structural Synergy'
)

## Descriptives
# png('acqlogit_descriptives_PLOT_313_ibm_MOST_GRANT.png', width = 8, height = 8, units = 'in', res = 200)
par(mfrow=c(3,3), mar=c(4,4,2.5,1))
for (var in names(X)) {
  if(var %in% c('y','t')) next
  hist(X[,var],col='gray',breaks=17, main=namemap[[var]], xlab='')
  abline(v=quantile(X[,var],c(.025,.5,.975)),col='darkred',lty=c(2,1,2),lwd=c(1,2,1))
}
# dev.off()

## REPLACE NAMES IN EFFECTS DF 
for (j in 1:ncol(X)) {
  names(X)[j] <- namemap[[names(X)[j]]]
}

xcor <- psych::corr.test(X, adjust = 'bonferroni')
write.csv(as.data.frame(xcor$r), 
          file = sprintf('acqlogit_data/acqlogit_correlations_table_%s_d%s_OrgSci20181107.csv',name_i,d))
write.csv(as.data.frame(xcor$p), 
          file = sprintf('acqlogit_data/acqlogit_corr_pvals_table_%s_d%s_OrgSci20181107.csv',name_i,d))

## DESCRIPTIVES
xdescr <- psych::describe(X[, ! names(X) %in% c('y','t')])
write.csv(as.data.frame(xdescr), 
          file = sprintf('acqlogit_data/acqlogit_descriptives_table_%s_d%s_OrgSci20181107.csv',name_i,d))




mgz <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + I(j.age^2) +  i.age +
    i.deg.full + j.deg.full +
    ij.same.region +
    i.fm.mmc.sum + i.acqs.100 + j.acqs.100 +
    ij.cossim +
    ij.dist +
    j.constraint +
    i.constraint + i.pow.n3 +
    i.rival.acq.3 + 
    ij.syn.constraint +
    # ij.syn.degree +
    # ij.syn.closeness2 +
    ij.syn.pow.n3 +
    I(i.fm.mmc.sum^2) +
    I(i.fm.mmc.sum^2):ij.syn.constraint
  ,
  data = df.sub)
summary(mgz)

# mg3 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
#     j.age + I(j.age^2) +  i.age +
#     i.deg.full + j.deg.full +
#     ij.same.region +  
#     i.fm.mmc.sum + i.acqs + j.acqs +
#     ij.cossim + 
#     ij.dist +
#     j.constraint + i.constraint +
#     ij.syn.constraint +
#     ij.syn.degree +
#     I(i.fm.mmc.sum^2) + 
#     I(i.fm.mmc.sum^2):ij.syn.degree
#   ,
#   data = df.sub)
# summary(mg3)

mg4 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.m2b + i.div +
    j.age + I(j.age^2) +  i.age +
    i.deg.full + j.deg.full +
    ij.same.region +  
    i.fm.mmc.sum + i.acqs + j.acqs +
    ij.cossim + 
    ij.dist +
    j.constraint + i.constraint +
    ij.syn.constraint +
    ij.syn.degree +
    I(i.fm.mmc.sum^2) + 
    I(i.fm.mmc.sum^2):ij.syn.constraint + 
    I(i.fm.mmc.sum^2):ij.syn.degree  
  ,
  data = df.sub)
summary(mg4)


mg3 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp + i.m2b +
    j.age + i.age + 
    i.deg + j.deg + 
    ij.same.region +  
    i.fm.mmc.sum + I(i.acqs * 100) +
    ij.cossim + 
    ij.dist + 
    j.constraint + i.constraint +
    react.rival.2 +
    ij.syn.constraint + 
    ij.syn.closeness + 
    I(i.fm.mmc.sum^2) + 
    I(i.fm.mmc.sum^2):ij.syn.constraint,
  data = df.sub)
summary(mg3)

mg4 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp + i.m2b +
    j.age + i.age + 
    i.deg + j.deg + 
    ij.same.region +  
    i.fm.mmc.sum + I(i.acqs * 100) +
    ij.cossim + 
    ij.dist + 
    j.constraint + i.constraint +
    react.rival.2 +
    ij.syn.constraint + 
    ij.syn.closeness + 
    I(i.fm.mmc.sum^2) + 
    I(i.fm.mmc.sum^2):ij.syn.closeness,
  data = df.sub)
summary(mg4)

mg5 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp + i.m2b +
    j.age + i.age + 
    i.deg + j.deg + 
    ij.same.region +  
    i.fm.mmc.sum + I(i.acqs * 100) +
    ij.cossim + 
    ij.dist + 
    react.rival.2 +
    j.constraint + i.constraint +
    ij.syn.constraint + 
    ij.syn.closeness + 
    I(i.fm.mmc.sum^2) + 
    # ij.syn.degree + 
    I(i.fm.mmc.sum^2):react.rival.2 +
    I(i.fm.mmc.sum^2):ij.syn.constraint + 
    I(i.fm.mmc.sum^2):ij.syn.closeness,
  # I(i.fm.mmc.sum^2):ij.syn.degree 
  data = df.sub)
summary(mg5)



##=======================================
##
## PLOT U-SHAPE MMC INTERACTIONS
##
##---------------------------------------
# ms <- c('-3 SD','-2 SD','-1 SD','Median','+1 SD','+2 SD','+3 SD')
# ms <- c('Min','-1 SD','Median','+1 SD','Max')
ms <- c('L99','1Q','Median','3Q','U99')
# ms <- c('L95','1Q','Median','3Q','U95')
# ms <- c('Min','1Q','Median','3Q','Max')
lwds <- c(3,3,3,1.5,3,3,3)
labs <- list(
  list(
    var='react.rival.2',
    ylab='Ln Acquisition Likelihood',
    xlab='Standardized FM-MMC',
    leg.title='Rival Acquisitions'
  ),
  list(
    var='ij.syn.constraint',
    ylab='Ln Acquisition Likelihood',
    xlab='Standardized FM-MMC',
    leg.title='Structural Synergy'
  )
)
fit <- mg5

png(sprintf('mmc_curve_interaction_%s_d%s_%s_mg9.png',name_i,d,length(ms)),width = 10,height = 4,units = 'in',res = 200)
par(mfrow=c(1,2), mar=c(4.2,4.2,.5,1.5))
for (i in 1:length(labs))
{
  var <- labs[[i]]$var
  mmc <- 'i.fm.mmc.sum'
  mmc.z <- 'mmc.z'
  fit.data <-  find_data(fit, parent.frame())
  
  ## MMC Z SCORE
  fit.data[,mmc.z] <- scale(fit.data[,mmc], center = T, scale = T)
  
  ## CREATE MEDIAN DATAFRAME (1 row)
  df.med <- fit.data
  for (col in names(df.med)) {
    if (is.numeric(df.med[[col]]) & col!=mmc.z) 
      df.med[,col] <- median(df.med[,col],na.rm = T)
  }
  df.med <- droplevels.data.frame(df.med)
  
  ## PREP PREDICTION DATAFRAME
  z <- seq(-3,3.3,length.out = 50)
  df.med <- df.med[1:length(z),]
  df.med[,mmc.z] <- z
  df.med[,mmc] <- z
  ## replace predicted data column
  df.comb <- df.med[mmc.z]  ## data.frame
  
  for (m in ms) 
  {
    ## VAR SETTING
    df.med[,var] <- if(m == '-1 SD'){
      median(fit.data[[var]]) -   sd(fit.data[[var]])
    }else if (m == '-2 SD'){
      median(fit.data[[var]]) - 2*sd(fit.data[[var]])
    }else if (m == '-3 SD'){
      median(fit.data[[var]]) - 3*sd(fit.data[[var]])
    }else if (m == '+1 SD'){
      median(fit.data[[var]]) +   sd(fit.data[[var]])
    }else if (m == '+2 SD'){
      median(fit.data[[var]]) + 2*sd(fit.data[[var]])
    }else if (m == '+3 SD'){
      median(fit.data[[var]]) + 3*sd(fit.data[[var]])
    }else if (m == '1Q'){
      quantile(fit.data[[var]], .25)
    }else if (m == '3Q'){
      quantile(fit.data[[var]], .75)
    }else if (m == 'L95'){  ## 2-sided
      quantile(fit.data[[var]], .025)
    }else if (m == 'U95'){  ## 2-sided
      quantile(fit.data[[var]], .975)
    }else if (m == 'L99'){  ## 2-sided
      quantile(fit.data[[var]], .005)
    }else if (m == 'U99'){  ## 2-sided
      quantile(fit.data[[var]], .995)
    }else if (m == 'Min'){
      min(fit.data[[var]]) 
    }else if (m == 'Max'){
      max(fit.data[[var]])
    }else{
      median(fit.data[[var]])
    }
    ## predictions
    pred <- prediction(fit, data = df.med)
    df.comb[,m] <- pred$fitted
  }
  
  ## PLOT
  pch <- c(15,2,16)
  matplot(df.comb[,1], df.comb[,-1],
          ylab=labs[[i]]$ylab,
          xlab=labs[[i]]$xlab,
          type='l', log='y', lwd=lwds,lty=1:length(ms),col=1:length(ms))
  # abline(v=c(-1.29,1.29),lty=4,col='gray',pch=pch)
  legend('bottomright',title=labs[[i]]$leg.title,
         legend=ms,lty=1:length(ms),col=1:length(ms), 
         lwd=lwds, cex=.8)
}
dev.off()

##---------------------------------------------------
##
##                end curviliean plots 
##
##---------------------------------------------------











###------------- ORG SCIENCE --------------------------------------




mc2 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim),
  data = df.sub)
summary(mc2)

mc3 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    I(i.fm.mmc.sum^2),
  data = df.sub)
summary(mc3)

mc3a <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    ij.dist,
  data = df.sub)
summary(mc3a)

mc3b <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    j.constraint,
  data = df.sub)
summary(mc3b)

mc3c <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    ij.syn.constraint,
  data = df.sub)
summary(mc3c)

mc3d <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    ij.syn.degree,
  data = df.sub)
summary(mc3d)

mc3e <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    ij.syn.closeness,
  data = df.sub)
summary(mc3e)

mc3z <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    I(100 * i.ij.cossim) + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint + 
    ij.syn.constraint + 
    ij.syn.degree + 
    ij.syn.closeness,
  data = df.sub)
summary(mc3z)

tabp <- mtable(mc2,mc3,mc3a,mc3b,mc3c,mc3d,mc3e,mc3z)
print(tabp)
## SAVE FILE
memisc::write.mtable(tabp, 
    file = sprintf("acqlogit_reg_table_313_primary_%s_d%s.tsv",name_i,d))

##=================================
## ORGSCI: FM-MMC INTERACTIONS
##--------------------------------
mc4 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree ,
  data = df.sub)
summary(mc4)

mc4a <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    I(i.fm.mmc.sum^2):ij.dist,
  data = df.sub)
summary(mc4a)

mc4b <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    I(i.fm.mmc.sum^2):j.constraint,
  data = df.sub)
summary(mc4b)

mc4c <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    I(i.fm.mmc.sum^2):ij.syn.closeness,
  data = df.sub)
summary(mc4c)

mc4d <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    I(i.fm.mmc.sum^2):ij.syn.constraint,
  data = df.sub)
summary(mc4d)

mc4e <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    I(i.fm.mmc.sum^2):ij.syn.degree,
  data = df.sub)
summary(mc4e)

mc4f <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    I(i.fm.mmc.sum^2):i.ij.cossim,
  data = df.sub)
summary(mc4f)

mc4g <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree + 
    I(i.fm.mmc.sum^2):ij.dist +
    I(i.fm.mmc.sum^2):j.constraint + 
    I(i.fm.mmc.sum^2):ij.syn.closeness + 
    I(i.fm.mmc.sum^2):ij.syn.constraint + 
    I(i.fm.mmc.sum^2):ij.syn.degree + 
    I(i.fm.mmc.sum^2):i.ij.cossim,
  data = df.sub)
summary(mc4g)


tabh1i <- mtable(mc4,mc4a,mc4b,mc4c,mc4d,mc4e,mc4f,mc4g)
print(tabh1i)
## SAVE FILE
memisc::write.mtable(tabh1i, 
file = sprintf("acqlogit_reg_table_313_h1interactions_%s_d%s.tsv",name_i,d))



##=================================
## ORGSCI: DISTANCE INTERACTIONS
##--------------------------------
mc5 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree ,
  data = df.sub)
summary(mc5)

mc5a <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.dist:j.constraint,
  data = df.sub)
summary(mc5a)

mc5b <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.dist:ij.syn.closeness,
  data = df.sub)
summary(mc5b)

mc5c <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.dist:ij.syn.constraint,
  data = df.sub)
summary(mc5c)

mc5d <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.dist:ij.syn.degree,
  data = df.sub)
summary(mc5d)


mc5e <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.dist:i.ij.cossim,
  data = df.sub)
summary(mc5e)

mc5f <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.dist:j.constraint + 
    ij.dist:ij.syn.closeness + 
    ij.dist:ij.syn.constraint + 
    ij.dist:ij.syn.degree + 
    ij.dist:i.ij.cossim,
  data = df.sub)
summary(mc5f)

tabdi <- mtable(mc5,mc5a,mc5b,mc5c,mc5d,mc5e,mc5f)
print(tabdi)
## SAVE FILE
memisc::write.mtable(tabdi, 
                     file = sprintf("acqlogit_reg_table_distance_interactions_%s_d%s.tsv",name_i,d))



##=================================
## ORGSCI: SYNERGEY DEGREE INTERACTIONS
##--------------------------------
mc6 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree ,
  data = df.sub)
summary(mc6)

mc6a <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:ij.dist,
  data = df.sub)
summary(mc6a)

mc6b <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:j.constraint,
  data = df.sub)
summary(mc6b)

mc6c <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:i.ij.cossim,
  data = df.sub)
summary(mc6c)

mc6d <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:ij.syn.closeness,
  data = df.sub)
summary(mc6d)

mc6e <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:ij.syn.constraint,
  data = df.sub)
summary(mc6e)

mc6z <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:j.constraint,
  data = df.sub)
summary(mc6z)

mc6f <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.degree:ij.dist +
    ij.syn.degree:j.constraint +
    ij.syn.degree:i.ij.cossim +
    ij.syn.degree:ij.syn.closeness +
    ij.syn.degree:ij.syn.constraint,
  data = df.sub)
summary(mc6f)

tabsdi <- mtable(mc6,mc6a,mc6b,mc6c,mc6d,mc6e,mc6f)
print(tabsdi)
## SAVE FILE
memisc::write.mtable(tabsdi, 
                     file = sprintf("acqlogit_reg_table_313_synergy_degree_interactions_%s_d%s.tsv",name_i,d)))


##=================================
## ORGSCI: SYNERGEY CLOSENESS INTERACTIONS
##--------------------------------
mc7 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree ,
  data = df.sub)
summary(mc7)

mc7a <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.closeness:ij.dist,
  data = df.sub)
summary(mc7a)

mc7b <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.closeness:j.constraint,
  data = df.sub)
summary(mc7b)

mc7c <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.closeness:i.ij.cossim,
  data = df.sub)
summary(mc7c)

mc7d <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.closeness:ij.syn.degree,
  data = df.sub)
summary(mc7d)

mc7e <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.closeness:ij.syn.constraint,
  data = df.sub)
summary(mc7e)

mc7f <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + ij.diff.deg + 
    ij.same.region +  ij.diff.deg   +
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    I(i.fm.mmc.sum^2) + 
    ij.dist + 
    j.constraint +
    ij.syn.closeness + 
    ij.syn.constraint + 
    ij.syn.degree +
    ij.syn.closeness:ij.dist +
    ij.syn.closeness:j.constraint +
    ij.syn.closeness:i.ij.cossim +
    ij.syn.closeness:ij.syn.degree +
    ij.syn.closeness:ij.syn.constraint,
  data = df.sub)
summary(mc7f)

tabsci <- mtable(mc7,mc7a,mc7b,mc7c,mc7d,mc7e,mc7f)
print(tabsci)
## SAVE FILE
memisc::write.mtable(tabsci, 
                     file = sprintf("acqlogit_reg_table_313_synergy_closeness_interactions_%s_d%s.tsv",name_i,d))















mg9 <- mclogit(
  cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
    j.age + i.age + ij.age.diff + 
    i.deg + j.deg + 
    ij.same.region +  
    i.fm.mmc.sum + I(i.acqs * 100) +
    i.ij.cossim + 
    ij.dist + 
    react.rival.2 +
    j.constraint + i.constraint +
    ij.syn.constraint + 
    ij.syn.closeness + 
    I(i.fm.mmc.sum^2) + 
    # ij.syn.degree +
    I(i.fm.mmc.sum^2):react.rival.2 +
    I(i.fm.mmc.sum^2):ij.syn.constraint + 
    I(i.fm.mmc.sum^2):i.ij.cossim,
  # I(i.fm.mmc.sum^2):ij.syn.degree 
  data = df.sub)
summary(mg9)








# 
# 
# 
# ##=======================================
# ##
# ## **QUANTILES**
# ## PLOT U-SHAPE MMC INTERACTIONS
# ##
# ##---------------------------------------
# # ms <- c('-3 SD','-2 SD','-1 SD','Median','+1 SD','+2 SD','+3 SD')
# # ms <- c('Min','-1 SD','Median','+1 SD','Max')
# ms <- c('L99','1Q','Median','3Q','U99')
# # ms <- c('Min','1Q','Median','3Q','Max')
# lwds <- c(2,2,2,1,2,2,2)
# labs <- list(
#   list(
#     var='react.rival.2',
#     ylab='Ln Acquisition Likelihood',
#     xlab='Standardized FM-MMC',
#     leg.title='Rival Acquisitions'
#   ),
#   list(
#     var='ij.syn.constraint',
#     ylab='Ln Acquisition Likelihood',
#     xlab='Standardized FM-MMC',
#     leg.title='Structural Synergy'
#   )
# )
# fit <- mg9
# 
# png(sprintf('mmc_curve_interaction_%s_sds%s_mg9.png',length(ms),'ibm'),width = 10,height = 4,units = 'in',res = 200)
# par(mfrow=c(1,2), mar=c(4.2,4.2,.5,1.5))
# for (i in 1:length(labs))
# {
#   var <- labs[[i]]$var
#   mmc <- 'i.fm.mmc.sum'
#   mmc.z <- 'mmc.z'
#   mmc.q <- 'mmc.q'
#   fit.data <-  find_data(fit, parent.frame())
#   
#   ## MMC Z SCORE
#   fit.data[,mmc.z] <- scale(fit.data[,mmc], center = T, scale = T)
#   ## Quantiles
#   fit.data[,mmc.q] <- fit.data[,mmc]
#   
#   ## CREATE MEDIAN DATAFRAME (1 row)
#   df.med <- fit.data
#   for (col in names(df.med)) {
#     if (is.numeric(df.med[[col]]) & col!=mmc.z & col!=mmc.q) 
#       df.med[,col] <- median(df.med[,col],na.rm = T)
#   }
#   df.med <- droplevels.data.frame(df.med)
#   
#   ## PREP PREDICTION DATAFRAME
#   z <- seq(-2.5,3.2,length.out = 50)
#   q <- seq(0,1,length.out = 50)
#   df.med <- df.med[1:length(z),]
#   df.med[,mmc.z] <- z
#   df.med[,mmc.q] <- quantile(fit.data[,mmc], q)
#   ##---------------
#   df.med[,mmc] <- q 
#   ##---------------
#   ## replace predicted data column
#   df.comb <- df.med[mmc.z]  ## data.frame
#   
#   for (m in ms) 
#   {
#     ## VAR SETTING
#     df.med[,var] <- if(m == '-1 SD'){
#       median(fit.data[[var]]) -   sd(fit.data[[var]])
#     }else if (m == '-2 SD'){
#       median(fit.data[[var]]) - 2*sd(fit.data[[var]])
#     }else if (m == '-3 SD'){
#       median(fit.data[[var]]) - 3*sd(fit.data[[var]])
#     }else if (m == '+1 SD'){
#       median(fit.data[[var]]) +   sd(fit.data[[var]])
#     }else if (m == '+2 SD'){
#       median(fit.data[[var]]) + 2*sd(fit.data[[var]])
#     }else if (m == '+3 SD'){
#       median(fit.data[[var]]) + 3*sd(fit.data[[var]])
#     }else if (m == '1Q'){
#       quantile(fit.data[[var]], .25)
#     }else if (m == '3Q'){
#       quantile(fit.data[[var]], .75)
#     }else if (m == 'L95'){  ## 2-sided
#       quantile(fit.data[[var]], .025)
#     }else if (m == 'U95'){  ## 2-sided
#       quantile(fit.data[[var]], .975)
#     }else if (m == 'L99'){  ## 2-sided
#       quantile(fit.data[[var]], .005)
#     }else if (m == 'U99'){  ## 2-sided
#       quantile(fit.data[[var]], .995)
#     }else if (m == 'Min'){
#       min(fit.data[[var]]) 
#     }else if (m == 'Max'){
#       max(fit.data[[var]])
#     }else{
#       median(fit.data[[var]])
#     }
#     ## predictions
#     pred <- prediction(fit, data = df.med)
#     df.comb[,m] <- pred$fitted
#   }
#   
#   ## PLOT
#   pch <- c(15,2,16)
#   matplot(df.comb[,1], df.comb[,-1],
#           ylab=labs[[i]]$ylab,
#           xlab=labs[[i]]$xlab,
#           type='l', log='y', lwd=lwds,lty=1:length(ms),col=1:length(ms))
#   # abline(v=c(-1.29,1.29),lty=4,col='gray',pch=pch)
#   legend('bottomright',title=labs[[i]]$leg.title,
#          legend=ms,lty=1:length(ms),col=1:length(ms), 
#          lwd=lwds, cex=.8)
# }
# dev.off()
# 
# ##---------------------------------------------------
# ##
# ##                end curviliean plots 
# ##
# ##---------------------------------------------------
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
# 
# 
# 
# 
# 
# 
# mr2 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age + ij.age.diff + 
#     i.deg + j.deg + ij.diff.deg + 
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.ij.cossim + 
#     I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     j.constraint +
#     ij.syn.closeness + 
#     ij.syn.constraint + 
#     ij.syn.degree +
#     react.rival.2 + 
#     I(i.fm.mmc.sum^2):react.rival.2,
#   data = df.sub)
# summary(mr2)
# 
# 
# mr3 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age + ij.age.diff + 
#     i.deg + j.deg + ij.diff.deg + 
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.ij.cossim + 
#     I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     j.constraint +
#     ij.syn.closeness + 
#     ij.syn.constraint + 
#     ij.syn.degree +
#     react.rival.3 + 
#     I(i.fm.mmc.sum^2):react.rival.3,
#   data = df.sub)
# summary(mr3)
# 
# mr4 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age + ij.age.diff + 
#     i.deg + j.deg + ij.diff.deg + 
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.ij.cossim + 
#     I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     j.constraint +
#     ij.syn.closeness + 
#     ij.syn.constraint + 
#     ij.syn.degree +
#     react.rival.4 + 
#     I(i.fm.mmc.sum^2):react.rival.4,
#   data = df.sub)
# summary(mr4)
# 
# mr5 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age + ij.age.diff + 
#     i.deg + j.deg + ij.diff.deg + 
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.ij.cossim + 
#     I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     j.constraint +
#     ij.syn.closeness + 
#     ij.syn.constraint + 
#     ij.syn.degree +
#     react.rival.5 + 
#     I(i.fm.mmc.sum^2):react.rival.5,
#   data = df.sub)
# summary(mr5)
# 
# mr10 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age + ij.age.diff + 
#     i.deg + j.deg + ij.diff.deg + 
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.ij.cossim + 
#     I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     j.constraint +
#     ij.syn.closeness + 
#     ij.syn.constraint + 
#     ij.syn.degree +
#     react.rival.10 + 
#     I(i.fm.mmc.sum^2):react.rival.10,
#   data = df.sub)
# summary(mr10)
# 
# mtable(mr1,mr2,mr3)
# 
# mc4g <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age + ij.age.diff + 
#     i.deg + j.deg +  
#     ij.same.region +  
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.ij.cossim + 
#     I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     j.constraint +
#     ij.syn.closeness + 
#     ij.syn.constraint + 
#     ij.syn.degree + 
#     I(i.fm.mmc.sum^2):ij.dist +
#     I(i.fm.mmc.sum^2):j.constraint + 
#     I(i.fm.mmc.sum^2):ij.syn.closeness + 
#     I(i.fm.mmc.sum^2):ij.syn.constraint + 
#     I(i.fm.mmc.sum^2):ij.syn.degree + 
#     I(i.fm.mmc.sum^2):i.ij.cossim,
#   data = df.sub)
# summary(mc4g)
# 
# 
# 
# prediction(mc4g, data = find_data(mc4g, parent.frame()), calculate_se = TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# mc4 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist,
#   data = df.sub)
# summary(mc4)
# 
# 
# mc5 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist ,
#   data = df.sub)
# summary(mc5)
# 
# mc5b <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist + 
#     I(i.fm.mmc.sum^2):j.constraint,
#   data = df.sub)
# summary(mc5b)
# 
# mc5c <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist + 
#     I(i.fm.mmc.sum^2):ij.dist,
#   data = df.sub)
# summary(mc5c)
# 
# mc5d <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist + 
#     I(i.fm.mmc.sum^2):ij.dist,
#   data = df.sub)
# summary(mc5d)
# 
# mc6 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist + ij.syn.closeness ,
#   data = df.sub)
# summary(mc6)
# 
# mc6a <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist + ij.syn.closeness +
#     I(i.fm.mmc.sum^2):ij.syn.closeness,
#   data = df.sub)
# summary(mc6a)
# 
# mc6b <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist + ij.syn.degree +
#     I(i.fm.mmc.sum^2):ij.syn.degree,
#   data = df.sub)
# summary(mc6b)
# 
# mc6c <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist + ij.syn.constraint +
#     I(i.fm.mmc.sum^2):ij.syn.constraint,
#   data = df.sub)
# summary(mc6c)
# 
# mc6d <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     j.age + i.age +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist  +
#     I(i.fm.mmc.sum^2):ij.dist,
#   data = df.sub)
# summary(mc6d)
# 
# ## SUMMARY TABLE
# mtable(mc0,mc1,mc2,mc3,mc4,mc6,mc6a,mc6b,mc6c,mc6d)
# 
# ## SAVE FILE
# mtab <- memisc::mtable(mc0,mc1,mc2,mc3,mc4,mc5,mc5b,mc5c)
# memisc::write.mtable(mtab, file = "acqlogit_reg_table_ibm.tsv")
# 
# write.csv(df.sub, file="acqlogit_reg_data_sample_ibm.csv", row.names = F)
# 
# 
# 
# 
# 
# 
# ###------------- SMS Proposal --------------------------------------
# 
# mc0 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp + 
#     
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 +  I(100 * i.ij.cossim) ,
#   data = df.sub)
# summary(mc0)
# 
# mc1 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2),
#   data = df.sub)
# summary(mc1)
# 
# mc2 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     j.constraint,
#   data = df.sub)
# summary(mc2)
# 
# mc3 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     ij.dist,
#   data = df.sub)
# summary(mc3)
# 
# mc4 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist,
#   data = df.sub)
# summary(mc4)
# 
# 
# mc5 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist + 
#     I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.sub)
# summary(mc5)
# 
# mc5b <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist + 
#     I(i.fm.mmc.sum^2):j.constraint,
#   data = df.sub)
# summary(mc5b)
# 
# mc5c <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist + 
#     I(i.fm.mmc.sum^2):ij.dist,
#   data = df.sub)
# summary(mc5c)
# 
# mc5d <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) +
#     ij.dist + 
#     I(i.fm.mmc.sum^2):ij.dist,
#   data = df.sub)
# summary(mc5d)
# 
# mc6 <- mclogit(
#   cbind(y,t) ~  i.ln_asset + i.cash + i.roa + i.ln_emp +
#     ij.same.region +  ij.diff.deg   +
#     i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + I(100 * i.ij.cossim) + 
#     I(i.fm.mmc.sum^2) + 
#     j.constraint + 
#     ij.dist + 
#     I(ij.dist^2),
#   data = df.sub)
# summary(mc6)
# 
# 
# 
# ## SUMMARY TABLE
# mtable(mc0,mc1,mc2,mc3,mc4,mc5,mc5b,mc5c)
# 
# ## SAVE FILE
# mtab <- memisc::mtable(mc0,mc1,mc2,mc3,mc4,mc5,mc5b,mc5c)
# memisc::write.mtable(mtab, file = "acqlogit_reg_table_ibm.tsv")
# 
# write.csv(df.sub, file="acqlogit_reg_data_sample_ibm.csv", row.names = F)
# 
# 
# 
# mtab <- memisc::mtable(mc0,mc1,mc2,mc3,mc4, signif.symbols = "")
# memisc::write.mtable(mtab, file = "sms_acq_logit_reg_table_ibm.txt", signif.symbols = "")
# 
# # correlations
# round(cor(mc4$model),2)
# 
# (smc0 <- mclogit::getSummary.mclogit(mc0, alpha=.05) )
# (smc1 <- mclogit::getSummary.mclogit(mc1, alpha=.05) )
# (smc2 <- mclogit::getSummary.mclogit(mc2, alpha=.05) )
# (smc3 <- mclogit::getSummary.mclogit(mc3, alpha=.05) )
# (smc4 <- mclogit::getSummary.mclogit(mc4, alpha=.05) )
# 
# smc4$coef
# 
# round(mc4$covmat, 3)
# 
# df.model <- cbind(predict(mc4),mc4$model)
# 
# 
# 
# 
# ## SANITY CHECK
# mc0 <- mclogit(
#   cbind(y,t) ~ i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.sub)
# summary(mc0)
# 
# ## add all controls to check
# mc1 <- mclogit(
#   cbind(y,t) ~ ij.same.state + ij.same.region + ij.same.country + ij.same.employee.range + ij.diff.deg + 
#     i.num.mkts + i.deg  + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.sub)
# summary(mc1)
# 
# ## preferable model (only useful controls)
# mc2 <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     i.deg  +  
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.sub)
# summary(mc2)
# 
# ## period effect  acquirer
# mc2r <- mclogit(
#   cbind(y,t) ~ ij.same.region +  ij.same.employee.range + ij.diff.deg +
#     i.deg  +  
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   random = ~ 1 | t,
#   control = mclogit.control(epsilon=1e-08, maxit=2000, trace=TRUE),
#   data = df.sub)
# summary(mc2r)
# 
# mc3r <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     i.deg  + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   random = ~ 1 | t,
#   control = mclogit.control(epsilon=1e-08, maxit=2000, trace=TRUE),
#   data = df.sub)
# summary(mc3r)
# 
# ## NESTED i in time t
# mc4r <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     i.deg  + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   random = ~ 1 | t / i,
#   control = mclogit.control(epsilon=1e-08, maxit=2000, trace=TRUE),
#   data = df.sub)
# summary(mc4r)
# 
# mc4 <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     i.deg  + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + 
#     i.pow.n2 ,
#   control = mclogit.control(epsilon=1e-07, maxit=200, trace=TRUE),
#   data = df.sub)
# summary(mc4)
# 
# 
# 
# ## preferable model (only useful controls)
# cols <- c('i.deg','i.fm.mmc.sum','ij.dist','i.pow.n2')
# df.sub.std <- df.sub
# df.sub.std[,cols] <- scale(df.sub.std[,cols], center = T, scale = T)
# 
# 
# mc2 <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     log(1 + i.deg) + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.sub.std)
# summary(mc2)
# 
# 
# 
# 
# 
# 
# 
# ## >= 10 constraint , power 
# ## >= 9 constraint no, power no
# ## >= 8 constraint no, power no
# ## >= 7 constraint no, power no
# ## >= 6 constraint no, power   no
# ## >= 5 constraint no, power  no
# ## >= 4 constraint no, power yes
# ## >= 3 constraint yes, power yes
# ## >= 2 constraint yes, power no
# mc3 <- mclogit(
#   cbind(y,t) ~  ij.same.region + ij.same.employee.range + ij.diff.deg   +
#     i.deg  +  i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + j.pow.n1 + ij.dist  +
#     i.constraint + ij.diff.constraint +
#     I(i.fm.mmc.sum^2) +
#     i.fm.mmc.sum:i.pow.n2 +
#     i.fm.mmc.sum:i.constraint +
#     i.fm.mmc.sum:ij.dist  , 
#   data = df.sub)
# summary(mc3)
# 
# mc2 <- mclogit(
#   cbind(y,t) ~  ij.same.region + ij.same.employee.range + ij.diff.deg   +
#     i.deg  +  i.fm.mmc.sum + I(i.acqs * 100) +
#     i.pow.n2 + j.pow.n1 + ij.dist  +
#     i.constraint + 
#     I(i.fm.mmc.sum^2) +
#     i.fm.mmc.sum:i.pow.n2 +
#     i.fm.mmc.sum:i.constraint +
#     i.fm.mmc.sum:j.constraint +
#     i.fm.mmc.sum:ij.dist  , 
#   data = df.sub)
# summary(mc2)
# 
# 
# memisc::mtable(mc1,mc2,mc3)
# 
# mtab <- memisc::mtable(mc0,mc1,mc2,mc3,mc4, signif.symbols = "")
# memisc::write.mtable(mtab, file = "sms_acq_logit_reg_table_ibm.txt", signif.symbols = "")
# 
# # correlations
# round(cor(mc4$model),2)
# 
# (smc0 <- mclogit::getSummary.mclogit(mc0, alpha=.05) )
# (smc1 <- mclogit::getSummary.mclogit(mc1, alpha=.05) )
# (smc2 <- mclogit::getSummary.mclogit(mc2, alpha=.05) )
# (smc3 <- mclogit::getSummary.mclogit(mc3, alpha=.05) )
# (smc4 <- mclogit::getSummary.mclogit(mc4, alpha=.05) )
# 
# smc4$coef
# 
# round(mc4$covmat, 3)
# 
# df.model <- cbind(predict(mc4),mc4$model)
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
# 
# 
# 
# 
# 
# 
# ## ----------------------------------------
# ##  Logit Model 
# ##------------------------------------------
# mc2 <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     log(1 + i.deg) + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.sub)
# summary(mc2)
# 
# mc2.i <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     log(1 + i.deg) + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   random = ~ 1 | i,
#   data = df.sub)
# summary(mc2.i)
# 
# ## ----------------------------------------
# ##  MIXED EFFECTS LOGISTIC
# ##------------------------------------------
# glm1 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 +
#     (1 | i) + (1 | t), 
#   verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=1e4)),
#   data=df.sub
# )
# summary(glm1)
# 
# 
# ##--------------------------------------
# ## only acquirer (focus on target selection)
# ##--------------------------------------
# acq.only <- df.sub[df.sub$y==1,c('t','i')]
# acq.only <-droplevels.data.frame(acq.only)
# #
# df.acq <- data.frame()
# df.cnt <- data.frame()
# for (i in 1:nrow(acq.only)) {
#   t <- as.numeric(acq.only$t[i])
#   firm_i <- as.character(acq.only$i[i])
#   tmp <- df.sub[as.numeric(df.sub$t) == t, ]
#   tmp <- tmp[as.character(tmp$i) == firm_i, ]
#   cat(sprintf("t %s : nrow %s\n",t,nrow(tmp)))
#   df.cnt <- rbind(df.cnt, data.frame(t=t,n=nrow(tmp)))
#   df.acq <- rbind(df.acq, tmp)
# }
# ## keep only nrow >= 3
# t.keep <- df.cnt$t[which(df.cnt$n >= 3)]
# df.acq <- df.acq[which(df.acq$t %in% t.keep), ]
# 
# acq1 <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     ij.dist + ij.diff.pow.n1 +
#     I(i.fm.mmc.sum^2):ij.dist ,
#   data = df.acq)
# summary(acq1)
# 
# 
# glm.acq1 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     ij.dist + 
#     I(i.fm.mmc.sum^2):ij.dist  +
#     (1  | i) , 
#   verbose = 9, family=binomial(link='logit'),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.acq
# )
# summary(glm.acq1)
# 
# 
# ##--------------------------------------
# ## only target (focus on MMC motivating acquirer)
# ##--------------------------------------
# targ.only <- df.sub[df.sub$y==1,c('t','j')]
# targ.only <-droplevels.data.frame(targ.only)
# #
# df.targ <- data.frame()
# df.cnt <- data.frame()
# for (i in 1:nrow(targ.only)) {
#   t <- as.numeric(targ.only$t[i])
#   j <- as.character(targ.only$j[i])
#   tmp <- df.sub[as.numeric(df.sub$t) == t, ]
#   tmp <- tmp[as.character(tmp$j) == j, ]
#   cat(sprintf("t %s : nrow %s\n",t,nrow(tmp)))
#   df.cnt <- rbind(df.cnt, data.frame(t=t,n=nrow(tmp)))
#   df.targ <- rbind(df.targ, tmp)
# }
# ## keep only nrow >= 3
# t.keep <- df.cnt$t[which(df.cnt$n >= 4)]
# df.targ <- df.targ[which(df.targ$t %in% t.keep), ]
# 
# write.csv(df.targ, file="sms_reg_df_targ_set_ibm.csv",row.names = F)
# 
# ##-----------------------------------------------
# ## LOGIT Acquirer set only
# ##-----------------------------------------------
# 
# 
# mc2t <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     log(1 + i.deg) + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   data = df.targ)
# summary(mc2t)
# 
# mc2tr.t <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     log(1 + i.deg) + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   random =  ~ 1 | t,
#   data = df.targ)
# summary(mc2tr.t)
# 
# mc2tr.i <- mclogit(
#   cbind(y,t) ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     log(1 + i.deg) + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2,
#   random =  ~ 1 | i,
#   data = df.targ)
# summary(mc2tr.i)
# 
# 
# ##--------------------------------------------
# ##
# ## SMS MIXED LOGIT REG TABLE
# ##
# ##--------------------------------------------
# m0 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     ij.dist +
#     i.pow.n2 +  
#     ( 1 | i ) + ( 1 | t ), 
#   # verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.targ
# )
# summary(m0)
# 
# m1 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist +
#     i.pow.n2 +  
#     ( 1 | i ) + ( 1 | t ), 
#   # verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.targ
# )
# summary(m1)
# 
# m2 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 +
#     ( 1 | i ) + ( 1 | t ), 
#   # verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.targ
# )
# summary(m2)
# 
# m3 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 + 
#     ( 1 | i ) + ( 1 | t ), 
#   # verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.targ
# )
# summary(m3)
# 
# m4 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 +
#     ( 1 | i ) + ( 1 | t ), 
#   # verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.targ
# )
# summary(m4)
# 
# 
# mlist <- list(m0=m0,m1=m1,m2=m2,m3=m3,m4=m4)
# screenreg(mlist, digits = 3, ci.force = T, ci.force.level = .95)
# screenreg(mlist, digits=3)
# htmlreg(mlist, file = "sms_acq_glmm_results_ibm.html", ci.force = T, ci.force.level = .95, digits = 3, star.symbol = '')
# 
# df.targ2 <- cbind(predict(m4),(1/(1+exp(-predict(m4)))),df.targ)
# write.csv(df.targ2, file="acqlogit_covs_pred_df_ibm.csv",row.names = F)
# 
# ##--------------end sms-----------------------
# 
# 
# ## NESTED -- NO
# glm.targ1 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 +
#     (1  | t/i ), 
#   verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=1e3)),
#   data=df.targ
# )
# summary(glm.targ1)
# 
# ## CROSSED -- YES
# glm.targ2 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 +
#     ( 1 | i ) + ( 1 | t ), 
#   verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=1e4)),
#   data=df.targ
# )
# summary(glm.targ2)
# 
# 
# glm.targ2 <- glmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg  +
#     log(1 + i.deg) +  i.fm.mmc.sum +
#     I(i.fm.mmc.sum^2) +
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 +
#     (1  | i/t ), 
#   verbose = 9, 
#   family=binomial(link = "logit"),
#   control=glmerControl(optimizer="bobyqa",
#                        optCtrl=list(maxfun=2e5)),
#   data=df.targ
# )
# summary(glm.targ2)
# 
# ##----------------------------------------
# ## Mixed LOGISTIC REGRESSION
# ##---------------------------------------
# cols <- c('i.deg','i.fm.mmc.sum','ij.dist','i.pow.n2')
# df.sub.std <- df.sub
# df.sub.std[,cols] <- scale(df.sub.std[,cols], center = T, scale = T)
# 
# 
# lm1 <- lmer(
#   y ~ ij.same.region + ij.same.employee.range + ij.diff.deg +
#     i.deg  + 
#     i.fm.mmc.sum + I(i.fm.mmc.sum^2) + 
#     ij.dist + I(i.fm.mmc.sum^2):ij.dist + 
#     i.pow.n2 +  I(i.fm.mmc.sum^2):i.pow.n2 + 
#     (1 | j) + (1 | i), 
#   verbose = 9, family=binomial,
#   data=df.sub.std
# )
# summary(lm1)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##--------------------------------------------------------------
# ## Examples
# ##--------------------------------------------------------------
# 
# 
# ##################################################
# ### MCLOGIT MIXED CONDITIONAL LOGIT
# data(Transport)
# 
# mc1 <- mclogit(
#   cbind(resp,suburb) ~ distance+cost,
#   data=Transport
# )
# summary(mc1)
# 
# data(electors)
# 
# mc2 <- mclogit(
#   cbind(Freq,interaction(time,class)) ~ econ.left/class+welfare/class+auth/class,
#   random=~1|party.time,
#   data=within(electors,party.time<-interaction(party,time)))
# summary(mc2)
# 
# 
# 
# ######################################################
# ## MCLOGIT  
# 
# data("Fishing", package = "mlogit")
# Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
# 
# ## a pure "conditional" model
# 
# summary(mlogit(mode ~ price + catch, data = Fish))
# 
# ## a pure "multinomial model"
# 
# summary(mlogit(mode ~ 0 | income, data = Fish))
# 
# ## which can also be estimated using multinom (package nnet)
# 
# library("nnet")
# summary(multinom(mode ~ income, data = Fishing))
# 
# ## a "mixed" model
# 
# m <- mlogit(mode ~ price+ catch | income, data = Fish)
# summary(m)
# 
# ## same model with charter as the reference level
# 
# m <- mlogit(mode ~ price+ catch | income, data = Fish, reflevel = "charter")
# 
# ## same model with a subset of alternatives : charter, pier, beach
# 
# m <- mlogit(mode ~ price+ catch | income, data = Fish,
#             alt.subset = c("charter", "pier", "beach"))
# 
# ## model on unbalanced data i.e. for some observations, some
# ## alternatives are missing
# 
# # a data.frame in wide format with two missing prices
# Fishing2 <- Fishing
# Fishing2[1, "price.pier"] <- Fishing2[3, "price.beach"] <- NA
# mlogit(mode~price+catch|income, Fishing2, shape="wide", choice="mode", varying = 2:9)
# 
# # a data.frame in long format with three missing lines
# data("TravelMode", package = "AER")
# Tr2 <- TravelMode[-c(2, 7, 9),]
# mlogit(choice~wait+gcost|income+size, Tr2, shape = "long",
#        chid.var = "individual", alt.var="mode", choice = "choice")
# 
# ## An heteroscedastic logit model
# 
# data("TravelMode", package = "AER")
# hl <- mlogit(choice ~ wait + travel + vcost, TravelMode,
#              shape = "long", chid.var = "individual", alt.var = "mode",
#              method = "bfgs", heterosc = TRUE, tol = 10)
# 
# ## A nested logit model
# 
# TravelMode$avincome <- with(TravelMode, income * (mode == "air"))
# TravelMode$time <- with(TravelMode, travel + wait)/60
# TravelMode$timeair <- with(TravelMode, time * I(mode == "air"))
# TravelMode$income <- with(TravelMode, income / 10)
# 
# # Hensher and Greene (2002), table 1 p.8-9 model 5
# TravelMode$incomeother <- with(TravelMode, ifelse(mode %in% c('air', 'car'), income, 0))
# nl <- mlogit(choice~gcost+wait+incomeother, TravelMode,
#              shape='long', alt.var='mode',
#              nests=list(public=c('train', 'bus'), other=c('car','air')))
# 
# # same with a comon nest elasticity (model 1)
# nl2 <- update(nl, un.nest.el = TRUE)
# 
# ## a probit model
# ## Not run: 
# pr <- mlogit(choice ~ wait + travel + vcost, TravelMode,
#              shape = "long", chid.var = "individual", alt.var = "mode",
#              probit = TRUE)
# 
# ## End(Not run)
# 
# 
# ## a mixed logit model
# rpl <- mlogit(mode ~ price+ catch | income, Fishing, varying = 2:9,
#               shape = 'wide', rpar = c(price= 'n', catch = 'n'),
#               correlation = TRUE, halton = NA,
#               R = 10, tol = 10, print.level = 0)
# summary(rpl)
# rpar(rpl)
# cor.mlogit(rpl)
# cov.mlogit(rpl)
# rpar(rpl, "catch")
# summary(rpar(rpl, "catch"))
# 
# ## End(Not run)
# 
# # # a ranked ordered model
# # data("Game", package = "mlogit")
# # g <- mlogit(ch~own|hours, Game, choice='ch', varying = 1:12,
# #             ranked=TRUE, shape="wide", reflevel="PC")
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
