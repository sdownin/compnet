
library(btergm)
library(psych)
library(igraph)
library(intergraph)
library(boot)


work_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2'
data_dir <- file.path(work_dir, 'firm_nets_rnr2\\firmctrl')
results_dir <- file.path(work_dir,'amj_rnr2_results\\firmctrl')
version_dir <- file.path(work_dir,'R','awareness_amj_rnr2')

aaf    <- source(file.path(version_dir,'amj_awareness_functions.R'))$value    ## aaf: awareness functions
cb     <- source(file.path(version_dir,'amj_cb_data_prep.R'))$value   
sdc    <- source(file.path(version_dir,'amj_sdc_coop.R'))$value  

##---------------------------------------------
##  D3 - Main Results
##---------------------------------------------
## SETTINGS
firm_i <- 'qualtrics'
d <- 3
R <- 2000
m_x <- 'm4'
nPeriod <- 11


## NETWORKS LIST
data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
nets <- readRDS(data_file)

## make MMC nets list
mmc <- lapply(nets, function(net) as.matrix(net %n% 'mmc'))
# cpc <- lapply(nets, function(net) as.matrix(net %n% 'coop'))
cpp <- lapply(nets, function(net) as.matrix(net %n% 'coop_past'))
cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
cossim <- lapply(nets, function(net) as.matrix(net %n% 'cat_cos_sim'))
# centjoin <- lapply(nets, function(net) as.matrix(net %n% 'joint_cent_pow_n0_4'))
# centratio <- lapply(nets, function(net) as.matrix(net %n% 'cent_ratio_pow_n0_4'))
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


##----------------------------------
## DEGREE OF SEPARATION RANGE
##----------------------------------
firm_i <- 'qualtrics'
sapply(nets, function(net){
  g <- asIgraph(net)
  d <- c(igraph::distances(g))
  return(quantile(d[which(d>0 & d<Inf)], c(0, .025, .25, .5, .75, .975, 1)))
})


## effects matrix (excluding edges)
X <- fit@effects[, -1]


###---------------------------------
## COREELATIONS
##----------------------------------
ct <- psych::corr.test(X, use = 'complete', method = 'pearson')
## write corrlations
corfile <- sprintf('correlation_%s_pd%s_R%s_%s_d%s.csv',firm_i,nPeriod,R,m_x,d)
write.csv(round(ct$r,2), file=file.path(results_dir,corfile))
## write pvals
cor.p.file <- sprintf('correlation_Pvals_%s_pd%s_R%s_%s_d%s.csv',firm_i,nPeriod,R,m_x,d)
write.csv(round(t(ct$p),2), file=file.path(results_dir,cor.p.file))

##-----------------------------------
## Summary Statistics
##-----------------------------------
de <- psych::describe(X)
## save descriptions
defile <- sprintf('descriptives_%s_pd%s_R%s_%s_d%s.csv',firm_i,nPeriod,R,m_x,d)
write.csv(round(de,2), file=file.path(results_dir,defile))


###=======================================
##
##  Number of Firms & Encounters 
##  
##----------------------------------------

firms <- c(
  'checkmarket',
  'clarabridge',				
  'cloudcherry',						
  'confirmit',				
  'customergauge',				
  'cx-index',					
  'empathica',					
  'feedback-lite',				
  'first-mile-geo',				
  'getfeedback',				
  'inqwise'	,				
  'leaderamp',				
  'medallia'	,				
  'qualtrics'	,				
  'myfeelback'	,			
  'promoter-io'	,			
  'satmetrix'	,						
  'snap-surveys-ltd'	,			
  'super-simple-survey'	,			
  'survata'	,						
  'surveyrock'	,						
  'typeform'	,			
  'userate'	,					
  'verint',				
  'voice-polls'	
)
d <- 3

ics <- c() ## indirect competitors & focal firm names vectors
edf <- data.frame()
el <- list()
for (firm_i in firms) {
  if (length(edf) > 1 & firm_i %in% edf$name)
    next
  cat(sprintf(' %s\n',firm_i))
  ## Encountrs 
  data_file <- file.path(data_dir,sprintf('%s_d%s.rds',firm_i,d))
  nets <- readRDS(data_file)
  n <- dim(nets[[1]][,])[1]
  es <- sapply(nets, function(net) {
      g <- asIgraph(net)
      return(c(e=ecount(g), v=length(which(igraph::degree(g)>0))))
    })
  ics <- c(ics,  nets[[5]] %v% 'vertex.names')
  cumu <- length(unique(ics))
  cats <- unlist(strsplit(nets[[5]] %v% 'category_list', '[|]'))
  rm(nets)
  el[[firm_i]] <- list(
    n=n, 
    e=es, 
    cumu=cumu,
    cats=cats
  )
  ##
  # edf <- rbind(edf, data.frame(name=firm_i, n=n, e=sum(es), cumu=cumu))
}
# edf <- rbind(edf, data.frame(name='TOTAL', n=sum(edf$n),e=sum(edf$e), cumu=edf$cumu[nrow(edf)]))
edf

## save encoutners count
ecrdsfile <- sprintf('firm_nets_encounters_sizes_pd%s_R%s_d%s.rds',nPeriod,R,d)
saveRDS(edf, file=file.path(data_dir, ecrdsfile))

es <- c()
vs <- c()
for (i in 1:length(el)) {
  es <- c(es, unname(el[[i]]$e[1,]))
  vs <- c(vs, unname(el[[i]]$e[2,]))
}

cat(sprintf('e=%s, v=%s, E[v]=%.2f, e/v=%.2f, med/med=%.2f',
            sum(es),sum(vs), mean(vs), sum(es)/sum(vs), median(es)/median(vs)))

## save encoutners count
ecfile <- sprintf('firm_nets_encounters_sizes_pd%s_R%s_d%s.csv',nPeriod,R,d)
write.csv(edf, file=file.path(data_dir, ecfile))



##---------------------------------------------
##
##---------------------------------------------
emp <- c(sapply(nets, function(net){net %v% 'employee_na_age'}))
summary(emp)
sal <- c(sapply(nets, function(net){net %v% 'sales_na_0_mn'}))
summary(sal)


## TESTING
m1 <- matrix(c(0,1,1,1, 0,0,1,0, 0,0,0,0, 0,0,0,0), nrow=4)
rank <- matrix(1:16, nrow=4)
val <- c(1.1,2.22,3.333,4.4444)
n1 <- network(m1, directed = F)
n1 %e% 'rank' <- rank
n1 %v% 'val' <- val

x1 <- ergmMPLE(n1 ~ edges + cycle(3) + nodecov('val') + absdiff('val') + edgecov(rank), output = 'array')

##---------------------------------------------
## Qualtrics number of (nonisolate) competitors
##---------------------------------------------
# comps <- sapply(nets, function(net){
#     g <- asIgraph(net)
#     return(length(which(igraph::degree(g)>0)))
#   })
# > comps
# 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 
# 119  146  192  223  274  308  370  410  440  477  477 
# > sum(comps)
# [1] 3436



##=============================================
##
##  COMPLEMENTS VS SUBSTITUTES
##
##---------------------------------------------
classfile <- 'CLASSIFICATION_qualtrics_d3_2016_comp_dist_table_JS_20180214.csv'
cdf <- read.csv(file.path(data_dir, classfile), stringsAsFactors = F)
cdf$Substitutes <- as.numeric(cdf$Substitutes)
cdf$Unrelated <- as.numeric(cdf$Unrelated)
cdf <- cdf[,1:8]
cdf$both <- apply(cdf[,c('Complements','Substitutes')],1,function(x) as.integer(sum(x)>1))
cdf$distance <- as.factor(cdf$distance)
cdf2 <- ddply(cdf, .(distance), summarize,
              Complements=sum(Complements, na.rm=T),
              Substitutes=sum(Substitutes, na.rm=T),
              Unrelated=sum(Unrelated, na.rm=T),
              both=sum(both, na.rm=T))
ds <- plyr::count(cdf$distance[cdf$distance!='0'])
cdf2p <- cdf2[-1,]
cdf2p[2,2:5] <- cdf2p[2,2:5] / ds$freq[2]
cdf2p[3,2:5] <- cdf2p[3,2:5] / ds$freq[3]
cdf2p[4,2:5] <- cdf2p[4,2:5] / ds$freq[4]
cdf2p    

library(reshape2)
library(ggplot2)
cdfl <- melt(cdf2, id.vars = 'distance',
             measure.vars = c('Complements','Substitutes','Unrelated'), 
             variable.name = 'Type')
cdfl$Separation <- as.integer(gsub('[+]','',cdfl$distance))
cdfl <- cdfl[cdfl$Separation>0,]
ggplot(aes(x=distance,y=value, colour=Type), data=cdfl) + 
  geom_boxplot()

##=============================================
##
##  RESULTS INTERPRET EFFECT SIZE
##
##---------------------------------------------
## Inverse logit
invlogit <- function(x) 1/(1+exp(-x))

## interpret effect size
effectsize <- function(X, theta, h, ps=c(.16,.5,.84), x.med=NA) {
  cat(sprintf('%s\n',h))
  if (all(is.na(x.med)))
    x.med <- apply(X, 2, median)
  
  q <- quantile(X[,h], ps)
  xdf <- data.frame(q=q, a=NA, p=NA)
  
  idx <- which(names(X) == h)
  xmlist <- list()
  for (i in 1:length(q)) {
    xmlist[[i]]      <- x.med
    xmlist[[i]][idx] <- q[i]
    xdf$a[i] <-  (xmlist[[i]] %*% theta)[1]
  }
  
  xdf$p <- sapply(xdf$a,invlogit)
  
  xdf$d.pct <- sapply(xdf$p, function(x) 100*(x - xdf$p[1])/xdf$p[1] )
  
  return(xdf)
}

ps <- c(.16,.5,.84)

### H1
m_x <- 'm1'
results_file <- file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
fits <- readRDS(results_file)
X <- fits[[firm_i]][[m_x]]@effects
theta <- fits[[firm_i]][[m_x]]@coef
x.med <- apply(X, 2, median)
(dh1 <- effectsize(X, theta, 'nodecov.genidx_multilevel', ps=ps, x.med=x.med))

### H2
m_x <- 'm2'
results_file <- file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
fits <- readRDS(results_file)
X <- fits[[firm_i]][[m_x]]@effects
theta <- fits[[firm_i]][[m_x]]@coef
x.med <- apply(X, 2, median)
(dh2 <- effectsize(X, theta, 'absdiff.cent_pow_n0_4', ps=ps, x.med=x.med))

### H3
m_x <- 'm3'
results_file <- file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
fits <- readRDS(results_file)
X <- fits[[firm_i]][[m_x]]@effects
theta <- fits[[firm_i]][[m_x]]@coef
x.med <- apply(X, 2, median)
ps <- seq(.55,.95,.05)
(dh31 <- effectsize(X, theta, 'cycle3', ps=ps, x.med=x.med))
(dh32 <- effectsize(X, theta, 'cycle4', ps=ps, x.med=x.med))
(dh33 <- effectsize(X, theta, 'cycle5', ps=ps, x.med=x.med))


dh1
dh2
dh31
dh32
dh33


## 
ps=c(.025,.16,.5,.84,.975)


## Qunatile (16, 84) // Robust alterative instead of +/- 1SD
## H1

hs <- c(
  'nodecov.genidx_multilevel',
  'absdiff.cent_pow_n0_4',
  'cycle3',
  'cycle4',
  'cycle5'
)


## Qunatile (16, 84) // Robust alterative instead of +/- 1SD
## H1

# hs <- c(
#   'nodecov.genidx_multilevel',
#   'absdiff.cent_pow_n0_4',
#   'cycle3',
#   'cycle4',
#   'cycle5'
# )
# hl <- list()
# for (h in hs)
#   hl[[h]] <- effsize(X, theta, h, x.med=x.med)
# 
# hl
# 
# effsize(X, theta, 'cycle3', ps=c(.025,.16,.5,.84,.975), x.med=x.med)

# > summary(X$nodecov.genidx_multilevel)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.         min for active firms
# 0.00000 0.08347 0.17584 0.28459 0.38333 2.90892      0.01149425
# 
# > summary(X$absdiff.cent_pow_n0_4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3044  0.7394  1.0976  1.4954 14.7126 
# 
# > summary(X$cycle3)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000  0.00000  0.00000  0.09188  0.00000 24.00000 
# 
# > summary(X$cycle4)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0000   0.0000   0.0000   0.7418   0.0000 148.0000 
# 
# > summary(X$cycle5)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    0.000    0.000    6.659    2.000 1329.000

## Effects
th1 <-  2.14 ##H1 (model 1)
th2 <-  0.32 ## H2 (model 2)
th3 <-  0.61 ## H3a (model 3)
th4 <-  0.22 ## H3b-4cycle (model 3)
th5 <- -0.02 ## H3b-5cycle (model 3)
## ODDS
sprintf('H1 +odds 100*(exp(%.2f)-1) =  %.1f%s', th1,  100 * (exp(th1) - 1) ,'%')
sprintf('H2 +odds 100*(exp(%.2f)-1) = %.1f%s', th2,  100 * (exp(th2) - 1) ,'%')
sprintf('H3a +odds 100*(exp(%.2f)-1) = %.1f%s', th3,  100 * (exp(th3) - 1) ,'%')
sprintf('H3b1 +odds 100*(exp(%.2f)-1) = %.1f%s', th4,  100 * (exp(th4) - 1) ,'%')
sprintf('H3b2 +odds 100*(exp(%.2f)-1) = %.1f%s', th5,  100 * (exp(th5) - 1) ,'%')



##---------------------------------------------
##  GOF COMPARISON (good vs bad) D2
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
shcomp <- lapply(nets, function(net) as.matrix(net %n% 'shared_competitor'))
shinv <- lapply(nets, function(net) as.matrix(net %n% 'shared_investor_nd'))

# ## TERGM RESULT
# results_file <- if (d==2){
#     file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s_d%s.rds',firm_i,nPeriod,R,m_x,d))
#   } else {
#     file.path(results_dir, sprintf('fit_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x))
#   } 
# fits <- readRDS(results_file)
# rm(fits)
# fit <- fits[[firm_i]][[m_x]]

## Model 4
form1 <- nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + 
  nodematch("ipo_status", diff = F) + nodematch("state_code",  diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + 
  nodecov("sales_na_0_mn") + edgecov(cossim) + edgecov(shinv) + 
  edgecov(mmc) + memory(type = "stability", lag = 1) + timecov(transform = function(t) t) + 
  nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + 
  absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
## No Network controls
form2 <- form2 <- nets ~ edges + 
  nodematch("ipo_status", diff = F) + nodematch("state_code",  diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + 
  nodecov("sales_na_0_mn") + edgecov(cossim) + edgecov(shinv) + 
  edgecov(mmc) + memory(type = "stability", lag = 1) +
  nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + 
  absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
## No Dyad Controls
form3 <- nets ~ edges + gwesp(0, fixed = T) + gwdegree(0, fixed = T) + 
  nodecov("age") +
  nodecov("employee_na_age") + 
  nodecov("sales_na_0_mn") + 
  memory(type = "stability", lag = 1) + timecov(transform = function(t) t) + 
  nodecov("genidx_multilevel") + nodecov("cent_pow_n0_4") + 
  absdiff("cent_pow_n0_4") + cycle(3) + cycle(4) + cycle(5)
## No structure
form4 <- nets ~ edges + 
  nodematch("ipo_status", diff = F) + nodematch("state_code",  diff = F) + nodecov("age") + absdiff("age") + nodecov("employee_na_age") + 
  nodecov("sales_na_0_mn") + edgecov(cossim) + edgecov(shinv) + 
  edgecov(mmc) 

# Run goodness of fit
fit1 <- btergm(form1, R = 100, parallel = 'multicore', ncpus = detectCores())
fit2 <- btergm(form2, R = 100, parallel = 'multicore', ncpus = detectCores())
fit3 <- btergm(form3, R = 100, parallel = 'multicore', ncpus = detectCores())
fit4 <- btergm(form4, R = 100, parallel = 'multicore', ncpus = detectCores())
#
gf1 <- gof(fit1, nsim=20, statistics=c(deg, dsp, esp, geodesic))  ## rocpr
gf2 <- gof(fit2, nsim=20, statistics=c(deg, dsp, esp, geodesic))  ## rocpr
gf3 <- gof(fit3, nsim=20, statistics=c(deg, dsp, esp, geodesic))  ## rocpr
gf4 <- gof(fit4, nsim=20, statistics=c(deg, dsp, esp, geodesic))  ## rocpr
#
gf1file <- sprintf('gof_compare_1_GOOD_%s_R%s_pd%s_d%s.png', firm_i, R, nPeriod, d)
gf2file <- sprintf('gof_compare_2_NET_%s_R%s_pd%s_d%s.png', firm_i, R, nPeriod, d)
gf3file <- sprintf('gof_compare_3_DYAD_%s_R%s_pd%s_d%s.png', firm_i, R, nPeriod, d)
gf4file <- sprintf('gof_compare_4_STRUCT_%s_R%s_pd%s_d%s.png', firm_i, R, nPeriod, d)
#
png(file.path(results_dir, gf1file), height = 8, width = 8, units = 'in', res=200)
  plot(gf1)
dev.off()
png(file.path(results_dir, gf2file), height = 8, width = 8, units = 'in', res=200)
  plot(gf2)
dev.off()
png(file.path(results_dir, gf3file), height = 8, width = 8, units = 'in', res=200)
  plot(gf3)
dev.off()

png(file.path(results_dir, gf4file), height = 8, width = 8, units = 'in', res=200)
  plot(gf4)
dev.off()

gf123file <- sprintf('gof_compare_1234_%s_R%s_pd%s_d%s.rds', firm_i, R, nPeriod, d)
saveRDS(list(gf1,gf2,gf3,gf4), file = file.path(results_dir, gf123file))

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

