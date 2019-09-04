cat('\n\n');timestamp();cat('\n')
library(btergm)
library(parallel)
library(texreg)

## DIRECTORIES
data_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/awareness-cues-amj"
work_dir <- data_dir

##===============================
## set analysis params
##-------------------------------
firm_i <- 'qualtrics'
d <- 2
ncpus <- 4
parallel <- "multicore"
nPeriods <- 11  ## 5
##-------------------------------

## load data
data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
nets <- readRDS(data_file)

#cache lists
if (nPeriods < length(nets))   
  nets <- nets[(length(nets)-nPeriods+1):length(nets)] 

cat("\n------------ estimating TERGM for:",firm_i,'--------------\n')
cat(sprintf("Using %s cores\n", detectCores()))

## make MMC nets list
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
  nodecov("employee_na_age") +
  nodecov("sales_na_0_mn") +
  edgecov(cossim) +
  edgecov(shinv) +   
  edgecov(mmc) + 
  memory(type = "stability", lag = 1) + 
  timecov(transform = function(t) t) +
  nodecov("genidx_multilevel") + 
  nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + 
  cycle(3) + cycle(4) + cycle(5) 

################################ end models#######################

##
# DEFINE MODEL and MODEL NAME TO COMPUTE
## 
m_x <- 'm4'
##

##==============================================
## BTERGM
##----------------------------------------------

# SET RESAMPLES
##
R <- 100

## RUN TERGM
fit4 <- btergm(get(m_x), R=R, parallel = parallel, ncpus = ncpus)

## SAVE SERIALIZED
fits.file <- sprintf('fit_%s_pd%s_R%s_%s_d%s.rds', firm_i, nPeriods, R, m_x, d)
saveRDS(fit4, file=file.path(data_dir,fits.file))
## SAVE FORMATTED REGRESSION TABLE
html.file <- sprintf('%s_tergm_results_pd%s_R%s_%s_d%s.html',  firm_i, nPeriods, R, m_x, d)
htmlreg(fit4, digits = 2, file=file.path(data_dir,html.file))


##================================================
## Goodness of Fit
##------------------------------------------------
gof4 <- gof(fit4, statistics=c(dsp, esp, deg, geodesic), nsim=10)
print(gof4)

## plot
plot(gof4)


##================================================
## COMPARE ESTIMATION ALGORITHM
##------------------------------------------------

set.seed(1111)
## estimate the TERGM with bootstrapped PMLE
fit4mc <- mtergm(get(m_x), ctrl=control.ergm(seed = 1111))

## SAVE SERIALIZED
fits.file <- sprintf('fit_%s_pd%s_R%s_%s_d%s.rds', firm_i, nPeriods, R, m_x, d)
saveRDS(fit4mc, file=file.path(data_dir,fits.file))
## SAVE FORMATTED REGRESSION TABLE
html.file <- sprintf('%s_tergm_results_pd%s_R%s_%s_d%s.html',  firm_i, nPeriods, R, m_x, d)
htmlreg(fit4mc, digits = 2, file=file.path(data_dir,html.file))
cat('finished successfully.')


## COMPARISON TABLE 
## Cache model fits list
fits <- list(PMLE=fit4, MCMCMLE=fit4mc)
## Echo model comparison table to screen
screenreg(fits, digits = 3)

## SAVE SERIALIZED
fits.file <- sprintf('fit_compare_%s_pd%s_R%s_%s_d%s.rds', firm_i, nPeriods, R, m_x, d)
saveRDS(fit4mc, file=file.path(data_dir,fits.file))
## SAVE FORMATTED REGRESSION TABLE
html.file <- sprintf('%s_compare_tergm_results_pd%s_R%s_%s_d%s.html',  firm_i, nPeriods, R, m_x, d)
htmlreg(fit4mc, digits = 2, file=file.path(data_dir,html.file))
cat('finished successfully.')


