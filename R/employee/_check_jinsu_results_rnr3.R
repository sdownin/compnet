
library(igraph)
library(intergraph)
library(btergm)
library(texreg)
library(parallel)

tutorial.nets <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\awareness-cues-amj\\qualtrics_d2.rds'
tutorial.fit <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\awareness-cues-amj\\fit_compare_qualtrics_pd11_R100_m4_d2.rds'
nets <- readRDS(tutorial.file)
fit1 <- readRDS(tutorial.fit)

fit1@coef

##

tutorial.fit2 <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\awareness-cues-amj\\fit_qualtrics_pd11_R100_m4_d2.rds'
fit2 <- readRDS(tutorial.fit)

fit2@coef




firm_i <- 'qualtrics'
d <- 2
ncpus <- 4
parallel <- "multicore"
nPeriods <- 11  ## 5

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
# SET RESAMPLES
##
R <- 30

## RUN TERGM
fitd2 <- btergm(get(m_x), R=R, parallel = parallel, ncpus = ncpus)
fitd2data <- btergm(get(m_x), R=8, parallel = parallel, ncpus = ncpus,
                  returndata = T)




##
## Jin-Su's results
nets3file <- 'C:\\Users\\T430\\Downloads\\jinsu_qualtrics_d3.rds'
nets <- readRDS(nets3file)
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
## change statistics covariates
fitd3data <- btergm(get(m_x), R=8, parallel = parallel, ncpus = ncpus,
                    returndata = T)




## Jin-Su's results
rep.fit3 <- 'C:\\Users\\T430\\Downloads\\jinsu_qualtrics_d2.rds'
fit3 <- readRDS(rep.fit3)

screenreg(list(fit1,fit2,fit3), 
          ci.force = T, ci.force.level =.95, single.row = T)



fit0 <- btergm(nets ~ edges + gwesp(0, fixed = T) + 
                 gwdegree(0, fixed = T) + cycle(3) + 
                 cycle(4) + cycle(5), 
               R=30, ncpus = detectCores(), )


## Jin-Su's results
rep.nets4 <- 'C:\\Users\\T430\\Downloads\\jinsu_qualtrics_d3.rds'
nets <- readRDS(rep.nets4)



fitd2file <- 'C:\\Users\\T430\\Downloads\\jinsu_fit_qualtrics_pd11_R2000_m4_d2.rds'
fitd2 <- readRDS(fitd2file)
fitd2 <- fitd2$qualtrics$m4

### Jin-Su's replication on 2019-07-16
fitd3file <- 'C:\\Users\\T430\\Downloads\\jinsu_fit_qualtrics_pd11_R2000_m4.rds'
fitd3 <- readRDS(fitd3file)
fitd3 <- fitd3$qualtrics$m4


#### check original fitted model from Venus server R2 submission
m4venusfile <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\employee_rnr2_results\\firmctrl\\fit_qualtrics_pd11_R2000_m4.rds'
fm4v <- readRDS(m4venusfile)
fm4v <- fm4v$qualtrics$m4

##
sapply(fm4v@data$networks, function(x)sum(x[,])/2)
sapply(fitd3@data$networks, function(x)sum(x[,])/2)
##
sapply(fm4v@data$cossim, function(x)sum(x[,])/2)
sapply(fitd3@data$cossim, function(x)sum(x[,])/2)
##XXXXXXX
sapply(fm4v@data$shinv, function(x)sum(x[,])/2)
sapply(fitd3@data$shinv, function(x)sum(x[,])/2)
##
sapply(fm4v@data$mmc, function(x)sum(x[,])/2)
sapply(fitd3@data$mmc, function(x)sum(x[,])/2)
##
sapply(fm4v@data$memory, function(x)sum(x[,])/2)
sapply(fitd3@data$memory, function(x)sum(x[,])/2)
##
sapply(fm4v@data$timecov1, function(x)sum(x[,])/2)
sapply(fitd3@data$timecov1, function(x)sum(x[,])/2)
##
sapply(fm4v@data$offsmat, function(x)sum(x[,])/2)
sapply(fitd3@data$offsmat, function(x)sum(x[,])/2)

# ####
# nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + 
#   nodematch("ipo_status", diff = F) + 
#   nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + 
#   nodecov("employee_na_age") + nodecov("sales_na_0_mn") + 
#   edgecov(cossim) + edgecov(shinv) + 
#   edgecov(mmc) + memory(type = "stability", lag = 1) + 
#   timecov(transform = function(t) t) + 
#   nodecov("genidx_multilevel") + 
#   nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + 
#   cycle(3) + cycle(4) + cycle(5)
##
sapply(fm4v@data$networks, function(net)summary(net %v% 'age'))
sapply(fitd3@data$networks, function(net)summary(net %v% 'age'))
##
sapply(fm4v@data$networks, function(net)summary(net %v% 'genidx_multilevel'))
sapply(fitd3@data$networks, function(net)summary(net %v% 'genidx_multilevel'))
##
sapply(fm4v@data$networks, function(net)summary(net %v% 'cent_pow_n0_4'))
sapply(fitd3@data$networks, function(net)summary(net %v% 'cent_pow_n0_4'))
##
sapply(fm4v@data$networks, function(net)summary(net %v% 'sales_na_0_mn'))
sapply(fitd3@data$networks, function(net)summary(net %v% 'sales_na_0_mn'))
##
sapply(fm4v@data$networks, function(net)summary(net %v% 'employee_na_age'))
sapply(fitd3@data$networks, function(net)summary(net %v% 'employee_na_age'))


##
shv <- lapply(fm4v@data$networks, function(net)as.matrix(net %n% 'shared_investor_nd'))
shj <- lapply(fitd3@data$networks, function(net)as.matrix(net %n% 'shared_investor_nd'))



###---------------------------------------
nl <- fm4v$qualtrics$m4@data$networks
rm(fm4v)
net <- nl$`2008`

dlpath <- 'C:\\Users\\T430\\Downloads'
cb2 <- readRDS(file.path(dlpath, 'cb_test.rds'))
ih2 <- readRDS(file.path(dlpath, 'ih_test.rds'))
nl2 <- readRDS(file.path(dlpath, 'nl_test.rds'))

dim(ih)
dim(ih2)

dim(cb$inv)
dim(cb2$inv)

dim(cb$co_rou)
dim(cb2$co_rou)

dim(cb$inv_rou)
dim(cb2$inv_rou)

dim(cb$co)
dim(cb2$co)


length(nl2)
net2 <- nl2$`2008`

any(c(net %n% 'shared_investor_nd' != net2 %n% 'shared_investor_nd'))

## main:  Venus R2 vs. jin-su whole net list
all(c(fm4v@data$networks$`2008` %n% 'shared_investor_nd' == fitd3@data$networks$`2008` %n% 'shared_investor_nd'))

## check jin-su 1period check vs. Venus R2 
all(c(nl2$`2008` %n% 'shared_investor_nd' == fm4v@data$networks$`2008` %n% 'shared_investor_nd'))

## check jin-su 1period check vs. jin-su whole net list
all(c(nl2$`2008` %n% 'shared_investor_nd' == fitd3@data$networks$`2008` %n% 'shared_investor_nd'))






##-------------------------------------------------------------
dlpath <- 'C:\\Users\\T430\\Downloads'
nl2f <- readRDS(file.path(dlpath, 'nl_full.rds'))

### jin-su full nets 2 vs. Venus R2
all(c(nl2f$`2008` %n% 'shared_investor_nd' == fitd3@data$networks$`2008` %n% 'shared_investor_nd'))


all(c(nl2f$`2008` %n% 'shared_investor_nd' == nl2$`2008` %n% 'shared_investor_nd'))

## partial computation of networks list is correct
## completing the networks list makes them wrong as of >=2008
## what is changed at the end of the make comp nets script?

# ####
# nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + 
#   nodematch("ipo_status", diff = F) + 
#   nodematch("state_code", diff = F) + nodecov("age") + absdiff("age") + 
#   nodecov("employee_na_age") + nodecov("sales_na_0_mn") + 
#   edgecov(cossim) + edgecov(shinv) + 
#   edgecov(mmc) + memory(type = "stability", lag = 1) + 
#   timecov(transform = function(t) t) + 
#   nodecov("genidx_multilevel") + 
#   nodecov("cent_pow_n0_4") + absdiff("cent_pow_n0_4") + 
#   cycle(3) + cycle(4) + cycle(5)
##

##---------
## cossim different in 2008, 2007, 2006
## shinv only different in 2008
## BUT DV nets same in all years
##---------

##XXX
which(c(nl2f$`2008` %n% 'shared_investor_nd' != nl2$`2008` %n% 'shared_investor_nd'))
##XXX
which(c(nl2f$`2008` %n% 'cat_cos_sim' != nl2$`2008` %n% 'cat_cos_sim'))

which(c(nl2f$`2008` %n% 'mmc' != nl2$`2008` %n% 'mmc'))

which(c(nl2f$`2008` %v% 'genidx_multilevel' != nl2$`2008` %v% 'genidx_multilevel'))
which(c(nl2f$`2008` %v% 'ipo_status' != nl2$`2008` %v% 'ipo_status'))
which(c(nl2f$`2008` %v% 'age' != nl2$`2008` %v% 'age'))
which(c(nl2f$`2008` %v% 'cent_pow_n0_4' != nl2$`2008` %v% 'cent_pow_n0_4'))
which(c(nl2f$`2008` %v% 'sales_na_0_mn' != nl2$`2008` %v% 'sales_na_0_mn'))
which(c(nl2f$`2008` %v% 'employee_na_age' != nl2$`2008` %v% 'employee_na_age'))
which(c(nl2f$`2008` %v% 'state_code' != nl2$`2008` %v% 'state_code'))

## DV nets compare all same
all(nl2f$`2006`[,] == nl2$`2006`[,])
all(nl2f$`2007`[,] == nl2$`2007`[,])
all(nl2f$`2008`[,] == nl2$`2008`[,])



all(nl2f$`2008` %v% 'category_list' == nl2$`2008` %v% 'category_list')

a <- aaf$.cov.categoryCosineSimilarity(nl2f$`2008`)
b <- aaf$.cov.categoryCosineSimilarity(nl2$`2008`)
assertthat::are_equal(a,b)

assertthat::are_equal(nl2f$`2008` %v% 'category_list', nl2$`2008` %v% 'category_list')


part <- nl2$`2008` %v% 'category_list'
full <- nl2f$`2008` %v% 'category_list'
for (i in 1:length(part)) {
  if (any(is.na(c(part[i],full[i])))) {
    cat(sprintf('i=%s\n [part] %s\n [full] %s\n\n',i,part[i],full[i]))
  } else if (part[i] != full[i]) {
    cat(sprintf('i=%s\n [part] %s\n [full] %s\n\n',i,part[i],full[i]))
  }
}



assertthat::are_equal(nl2f$`2008` %n% 'shared_investor_nd', nl2$`2008` %n% 'shared_investor_nd')
which(nl2f$`2008` %n% 'shared_investor_nd' != nl2$`2008` %n% 'shared_investor_nd')



####------------------------------------
### new stephen's laptop vs. Venus R2
slaptop.quald3 <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\firm_nets_rnr2_repl\\qualtrics_d3.rds'
slnets <- readRDS(slaptop.quald3)

assertthat::are_equal(slnets$`2008` %n% 'shared_investor_nd',
                      fm4v@data$networks$`2008` %n% 'shared_investor_nd')

assertthat::are_equal(slnets$`2016` %n% 'cat_cos_sim',
                      fm4v@data$networks$`2016` %n% 'cat_cos_sim')



###----------------
## D2 networks check
### new stephen's laptop vs. Venus R2
r2.quald2 <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\employee_rnr2_results\\firmctrl\\fit_qualtrics_pd11_R2000_m4_d2.rds'
r2fitd2 <- readRDS(r2.quald2)
r2fitd2 <- r2fitd2$qualtrics$m4

slaptop.quald2 <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\firm_nets_rnr2_repl\\qualtrics_d2.rds'
slnets2 <- readRDS(slaptop.quald2)

assertthat::are_equal(slnets2$`2008` %n% 'shared_investor_nd',
                      r2fitd2@data$networks$`2008` %n% 'shared_investor_nd')

assertthat::are_equal(slnets2$`2008` %n% 'cat_cos_sim',
                      r2fitd2@data$networks$`2008` %n% 'cat_cos_sim')





###====================================
##  Replicate Results
##-------------------------------------
.repl <- function() {
  net_repl_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\firm_nets_rnr2_repl'
  res_repl_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2\\employee_rnr2_results_repl'
  for (d in 2:3)
  {
    cat(sprintf('\n\nd = %s\n\n',d))
    ## load
    nets <- readRDS(file.path(net_repl_dir, sprintf('qualtrics_d%s.rds',d)))
    len <- length(nets)
    nets <- nets[(len-10):len]
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
    ## make MMC nets list
    mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
    cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
    cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
    cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
    cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
    # centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
    # centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
    # shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor')) 
    shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))
    ## change statistics covariates
    fitdx <- btergm(m4, R=40, parallel='multicore', ncpus=detectCores())
    ## check observations
    cat(sprintf(' [nobs  NO ] %s\n [nrows YES] %s\n',fitdx@nobs, nrow(fitdx@effects)))
    ## echo results table 
    screenreg(fitdx, ci.force = T, ci.force.level = .95)
    ## save fitted output
    saveRDS(fitdx, file=file.path(res_repl_dir,sprintf('fit_qualtrics_R40_d%s.rds',d)))
  }
}
.repl()


library(psych)
X <- fitdx@effects[,-1]  ## exclude edges intercept
descr <- describe(X, na.rm = T)
rcor <- psych::corr.test(X)

saveRDS(list(descr=descr, rcor=rcor),
        file=file.path(res_repl_dir, 'cor_summary_qualtrics_d3.rds'))
##
View(cbind(rownames(descr),round(cbind(descr$mean,descr$sd),5)))
