
library(btergm)


work_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2'
data_dir <- file.path(work_dir, 'firm_nets_rnr2\\firmctrl')
results_dir <- file.path(work_dir,'amj_rnr2_results\\firmctrl')


##---------------------------------------------
##  D3 - Main Results
##---------------------------------------------
## SETTINGS
firm_i <- 'qualtrics'
d <- 3
R <- 2000
m_x <- 'm4'
nPeriod <- 11
# Number of Simulations
nsim <- 2

## NETWORKS LIST
data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
nets <- readRDS(data_file)

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

## TERGM RESULT
results_file <- if (d==2){
  file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s_d%s.rds',firm_i,nPeriod,R,m_x,d))
} else {
  file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
} 
fits <- readRDS(results_file)
fit <- fits[[firm_i]][[m_x]]

gf <- gof(fit, nsim=nsim, statistics=c(deg, dsp, esp, geodesic), verbose=T)  ## rocpr
gf.file <- sprintf('gof_%s_pd%s_R%s_%s_nsim%s_WINLOCAL.rds',firm_i,nPeriod,R,m_x,nsim)
saveRDS(gf, file = file.path(results_dir,gf.file))

plot(gf)

sim



##---------------------------------------------
##  D2
##---------------------------------------------
## SETTINGS
firm_i <- 'qualtrics'
d <- 2
R <- 2000
m_x <- 'm4'
nPeriod <- 11
nsim <- 10


## NETWORKS LIST
data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
nets <- readRDS(data_file)

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

## TERGM RESULT
results_file <- if (d==2){
    file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s_d%s.rds',firm_i,nPeriod,R,m_x,d))
  } else {
    file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
  } 
fits <- readRDS(results_file)
rm(fits)
fit <- fits[[firm_i]][[m_x]]

gf <- gof(fit, nsim=nsim, statistics=c(deg, dsp, esp, geodesic))  ## rocpr

plot(gf)
plot(gf1)


##---------------- GOF -------------------------------
## READ IN GOT FILES
nsim <- 10
R <- 2000
d <- 2
m_x <- 'm4'
firm_i <- 'qualtrics'
nPeriod <- 11 ## before deducline lag 1

fits <- readRDS(file.path(results_dir,'fit_qualtrics_pd11_R2000_m4_d2.rds'))
fit <- fits[[firm_i]][[m_x]]
rm(fits)
gf <- gof(fit, nsim=10, statistics=c(deg, dsp, esp, geodesic))  ## rocpr

gof_file <- sprintf('gof_%s_pd%s_R%s_%s_nsim%s_ALL_WINLOCAL.rds',firm_i,nPeriod,R,mx,nsim)
saveRDS(gf, file=file.path(results_dir,gof_file))



## READ IN GOT FILES
nsim <- 50
R <- 2000
m_x <- 'm4'
firm_i <- 'qualtrics'
nPeriod <- 11 ## before deducline lag 1

items <- c('deg','esp','dsp','geodesic')
gfl <- list()
for (item in items) {
  filestr <- sprintf('gof_%s_pd%s_R%s_%s_nsim%s_%s.rds', firm_i, nPeriod, R, m_x, nsim, item)
  gfl[[item]] <- readRDS(file.path(results_dir, filestr))
}
gfl$deg


plot(gfl$geodesic)


##----------------------- CORRELATIONS
library(psych)

X <- fits$m4@effects
X <- X[X$nodecov.age < 110 & X$nodecov.age >= 0, ]

(desc <- psych::describe(X))
write.csv(desc, file = 'qualtrics_m4_describe.csv', row.names = T)

(cr <- cor(X))
write.csv(cr, file = 'qualtrics_m4_correlations.csv', row.names = T)


(ct <- psych::corr.test(X))
write.csv(ct$p, file = 'qualtrics_m4_corr_test_p.csv', row.names = T)

