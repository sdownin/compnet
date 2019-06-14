##----------------------------------------------------------
##
##   AMJ TERGM INTERPRET
##
##----------------------------------------------------------

work_dir <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/compnet2"
data_dir <- file.path(work_dir, 'firm_nets_rnr2','firmctrl')
results_dir <- file.path(work_dir,'amj_rnr2_results','firmctrl')
img_dir  <- "C:/Users/T430/Google Drive/PhD/Dissertation/competition networks/envelopment/img"


## set woring dir
setwd(work_dir)

## load awareness functions
aaf <- source(file.path(getwd(),'R','awareness_amj_rnr2','amj_awareness_functions.R'))$value

## -----------Model Results Settings-----
name_i <- 'qualtrics'
firm_i <- name_i
d <- 3
R <- 2000
nPeriod <- 11
m_x <- 'm4'
##----------------------------------------

## FUNCTIONS FOR NAs REMOVED
.med <- function(x){return(median(x, na.rm=TRUE))}
.min <- function(x){return(min(x, na.rm=TRUE))}
.max <- function(x){return(max(x, na.rm=TRUE))}
.avg <- function(x){return(mean(x, na.rm=TRUE))}
.std <- function(x){return(sd(x, na.rm=TRUE))}
.qtl <- function(x, probs){return(quantile(x, probs=probs, na.rm=TRUE))}
.iqr <- function(x){return(IQR(x, na.rm=TRUE))}


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


## create igraph lists
# gs <- lapply(nets, asIgraph)


## visualize effects bootstrap normality
effects <- names(fit@effects)
par(mfrow=c(3,3))
for (effect in effects) {
  qqnorm(fit@boot$t[,effect], main=effect)
}
## visualize effects bootstrap distributions
par(mfrow=c(3,3))
for (effect in effects) {
  x <- fit@boot$t[,effect]
  dr <- diff(range(x))
  xl1 <- min(0-(.1*dr), min(x))
  xl2 <- max(.1*dr, max(x))
  ci <- quantile(x, c(.025,.975))
  hist(x, main=effect, breaks=13, col=rgb(.2,.2,.2,.2), xlim=c(xl1,xl2))
  abline(v=ci, col='red',lwd=2)
  abline(v=0, col='black',lwd=1, lty=2)
  segments(x0=ci[1],x1=ci[2],y0=0,y1=0,col='red',lwd=2)
}

# library()
# par(mfrow=c(1,1))
# ci95 <- c(.025,.5,.975)
# coef.mat <- t(apply(fits@boot$t,2,function(x)quantile(x,ci95)))
# matplot(x=df.coef, y=1:nrow(coef.mat))

# 
# ##===============================================
# ## interpreation data frame i-->j (for all t)
# ##-----------------------------------------------
# g <- asIgraph(nets[[length(nets)]])
# vcnt <- vcount(g)
# time.steps <- fits@time.steps
# firm.names <-  V(g)$vertex.names
# v.focal <- as.integer( V(g)[V(g)$vertex.names==name_i] )
# 
# ## make data frame of predicted probabilities
# # for (t in 1:time.steps)
# #   idf[,paste0('t',t)]<- NA
# # ## time columns
# # tcolnames <- names(idf)[grep(names(idf), pattern = 't\\d{1,3}', perl = T)]
# 
# ## competitor probability micro-interpretation  file
# interp.df.file <- sprintf('%s/interpret_%s_pd%s_d%s_R%s_%s.csv', 
#                           results_dir, name_i, nPeriods, d, R, m_x)
# 
# ## main loop: competitor probability micro-interpretation
# if (file.exists(interp.df.file)) {
#   ## READ IN PREDICTED PROBS
#   idf <- read.csv(interp.df.file)
# } else {
#   idf <- data.frame()
#   for (j in 1:vcount(g)) {
#     cat(sprintf('%s:  %s\n',j,V(g)$vertex.names[j]))
#     j.name <- V(g)$vertex.names[j]
#     if (j == v.focal) {
#       probs <- rep(0, time.steps)
#     } else {
#       probs <- btergm::interpret(fits, type='tie', i=v.focal, j=j)
#     }
#     for (t in 1:length(probs)) {
#       d <- igraph::distances(asIgraph(nets[[t+1]]), v = v.focal, to = j)[1,1]
#       tmp <- data.frame(i=v.focal, j=j, i.name=name_i, j.name=j.name, t=t, d=d, p=probs[t])
#       idf <- rbind(idf, tmp)
#     }
#     if (j %% 10 == 0)  write.csv(x = idf, file = interp.df.file)
#   }
#   write.csv(x = idf, file = interp.df.file)
# }

###
## READ IN Interpreataion Predicted Probabilities
##
idf.file <- sprintf('interpret_%s_pd%s_R%s_%s.csv', firm_i, nPeriod, R, m_x)
idf <- read.csv(file.path(results_dir,idf.file))


## add covariates
idf$genidx_multilevel <- NA
idf$njobs_multilevel <- NA
idf$cent_pow_n0_4 <- NA
idf$absdiff_pow_n0_4 <- NA
idf$year <- NA
for (row in 1:nrow(idf)) {
  t <- idf$t[row]
  i <- idf$i[row]
  j <- idf$j[row]
  # idf$d[row] <- igraph::distances(gs[[t+1]], v = i, to = j)[1,1]
  idf$year[row] <- 2006 + as.integer(t)
  idf$genidx_multilevel[row] <- (nets[[t+1]] %v% 'genidx_multilevel')[j]
  idf$njobs_multilevel[row] <- (nets[[t+1]] %v% 'njobs_multilevel')[j]
  idf$cent_pow_n0_4[row] <- (nets[[t+1]] %v% 'cent_pow_n0_4')[j]
  idf$absdiff_pow_n0_4[row] <- abs((nets[[t+1]] %v% 'cent_pow_n0_4')[j] - (nets[[t+1]] %v% 'cent_pow_n0_4')[i])
  if(row %% 300 == 0) cat(sprintf('row %s\n',row))
}

## remove probabilities and covariates (set NA) for Inf distance (firm's not in component that pd)
idf$p[idf$d==Inf] <- NA
idf$genidx_multilevel[idf$d==Inf] <- NA
idf$njobs_multilevel[idf$d==Inf] <- NA
idf$cent_pow_n0_4[idf$d==Inf] <- NA
idf$absdiff_pow_n0_4[idf$d==Inf] <- NA

## distance category
idf$d_cat <- as.factor(sapply(idf$d,function(x){
  xstr <- as.character(x)
  ifelse(xstr %in% c('1','2','3','4'), xstr, ifelse(x==Inf, NA, '5+'))
}))
idf$i <- as.factor(idf$i)
idf$year <- as.factor(idf$year)

## distance category
idf$d_cat <- as.factor(sapply(idf$d,function(x){
  xstr <- as.character(x)
  ifelse(xstr %in% c('1','2','3','4'), xstr, ifelse(x==Inf, NA, '5+'))
}))
idf$i <- as.factor(idf$i)
idf$year <- as.factor(idf$year)

## aware.all.cutoff
aware.frac <- .28
xp <- idf$p[idf$year=='2016' & !is.na(idf$p)]
aware.all.cutoff <- .qtl(xp, 1 - aware.frac)

## mutate group-period factors
idf2 <-
  idf[idf$d != 0 & idf$d !='0', ] %>%
  group_by(d) %>%
  mutate(
    outlier = p > .med(p) + .iqr(p) * 2, 
    low.otl = p < .med(p) - .iqr(p) * 2,
    aware.med = p > .qtl(p, .5),
    aware.cutoff = p > aware.all.cutoff
  ) %>% 
  ungroup

## manual colors
colors2rw <- c('black','red')
colors2 <- c( "#333333", '#aaaaaa')
colors4 <- c( '#222222', '#aaaaaa', '#555555', '#888888')
colors5 <- c( '#111111', '#333333','#555555', '#777777','#999999')


## Data for hypotheses interaction plots
# idf3 <- idf2[idf2$year %in% c('2007','2008','2009','2010','2011','2012','2013','2014','2015','2016'), ]

##====================================================
##
## Hostility Profile (awareness set) Scatter Facet plot
##  - EARLY -
##
##----------------------------------------------------
## set years
yrs <- as.character(2007:2016)
for (yr in yrs) {
  ## aware.all.cutoff
  aware.frac <- ifelse(yr < '2013', .25, .35)
  xp <- idf$p[idf$year==yr & !is.na(idf$p)]
  aware.all.cutoff <- .qtl(xp, 1 - aware.frac)
  ## mutate group-period factors
  idf2 <-
    idf[idf$d != 0 & idf$d !='0', ] %>%
    group_by(d) %>%
    mutate(
      outlier = p > .med(p) + .iqr(p) * 2, 
      low.otl = p < .med(p) - .iqr(p) * 2,
      aware.med = p > .qtl(p, .5),
      aware.cutoff = p > aware.all.cutoff
    ) %>% 
    ungroup
  ##
  idf3 <- idf2[idf2$year %in% yr, ]
  # idf3 <- idf3[ idf3$d != Inf | idf3$year != '2016', ]
  set.seed(269951)
  dodge.width <- 0
  ##
  ggplot(idf3) + aes(x = d_cat, y = p, color=aware.cutoff, pch=aware.cutoff) +
    geom_point(position=position_jitterdodge(dodge.width=dodge.width), lwd=2.5,
               data=function(x)dplyr::filter_(x, ~ !low.otl)) +
    scale_color_manual(values=colors2rw) +
    scale_shape_manual(values=c(1,16)) +
    facet_grid(year ~ .) +  # facet_wrap(~ year) +
    xlab("Degree of Separation from Focal Firm") + 
    scale_y_log10("Conditional Probability of Competitive Encounter", limits=c(2e-5, 40)) +
    geom_hline(yintercept = aware.all.cutoff, lty=3) +
    theme_classic() +  theme(legend.position='none')
  fig3.scatter.name <- sprintf('interpret_RNR2_%s_Fig3_scatter_cut_classic_%.4f_y%s.png', name_i, aware.all.cutoff, yr)
  ggsave(filename = fig3.scatter.name, width = 5.5, height = 5.5, units = 'in', dpi = 250)
}



## 2 year factor plot ----------------------
## set years
yrs <- c('2010','2016')
##
idf3 <- idf2[idf2$year %in% yrs, ]
idf3 <- idf3[ idf3$d != Inf | idf3$year != '2016', ]
set.seed(269951)
dodge.width <- 0
##
ggplot(idf3) + aes(x = d_cat, y = p, color=aware.cutoff, pch=aware.cutoff) +
  geom_point(position=position_jitterdodge(dodge.width=dodge.width), lwd=2.5,
             data=function(x)dplyr::filter_(x, ~ !low.otl)) +
  scale_color_manual(values=colors2rw) +
  scale_shape_manual(values=c(1,16)) +
  facet_grid(year ~ .) +  # facet_wrap(~ year) +
  xlab("Degree of Separation from Focal Firm") + 
  scale_y_log10("Conditional Probability of Competitive Encounter", limits=c(2e-5, 40)) +
  geom_hline(yintercept = aware.all.cutoff, lty=3) +
  theme_classic() +  theme(legend.position='none')
fig3.scatter.name <- sprintf('interpret_RNR2_%s_Fig3_scatter_cut_classic_%.4f_y%s-%s.png', name_i, aware.all.cutoff, yrs[1],yrs[2])
ggsave(filename = fig3.scatter.name, width = 4.5, height = 9, units = 'in', dpi = 250)




##====================================================
##
##     AERIAL MAPS
##
##-----------------------------------------------------
# quantile.cutoff <- .8
# cutoff <- quantile(idf$p[idf$year==2016], quantile.cutoff, na.rm = T)
seeds <- 11120:11130
plot.years <- c(2007:2009,2015:2016)
for (seed in seeds) {
  for (year in plot.years) {
    yr <- as.character(year)
    #
    aware.frac <- ifelse(yr < '2013', .25, .35)
    xp <- idf$p[idf$year==yr & !is.na(idf$p)]
    aware.all.cutoff <- .qtl(xp, 1 - aware.frac)
    ## mutate group-period factors
    idf2 <-
      idf[idf$d != 0 & idf$d !='0', ] %>%
      group_by(d) %>%
      mutate(
        outlier = p > .med(p) + .iqr(p) * 2, 
        low.otl = p < .med(p) - .iqr(p) * 2,
        aware.med = p > .qtl(p, .5),
        aware.cutoff = p > aware.all.cutoff
      ) %>% 
      ungroup
    ##
    idf3 <- idf2[idf2$year %in% yr, ]
    #
    gLatest <- asIgraph(nets[[ which(names(nets)==as.character(year+1)) ]])
    #
    idf.sub <- idf[idf$year==as.character(year), ]
    idf.sub[is.na(idf.sub$p), 'p'] <- 0
    #
    ## COMPUTE CUTOFF FRACTION FROM aware.all.cutoff SPECIFIED ABOVE
    quantile.cutoff <- length(idf.sub$p[idf.sub$p >= aware.all.cutoff] ) / length(idf.sub$p)
    #
    aerial.map.name <- sprintf('interpret_RNR2_%s_Fig3_aerial_map_cut%.4f_quant%s_y%s_sd%s.png', name_i, aware.all.cutoff, round(100*quantile.cutoff), year, seed)
    png(aerial.map.name, height=7.5, width=7.5, units='in', res = 250)
    par(mar=c(.1,.1,.1,.1), mfrow=c(1,1))
    V(gLatest)$prob <- idf.sub$p
    aaf$plotCompNetColPredict(gs = gLatest,
                              focal.firm = name_i,
                              cutoff=aware.all.cutoff,
                              probAttrName='prob',
                              layout.algo = layout.fruchterman.reingold,
                              seed=seed)
    dev.off()
  }
}

###--------------






##====================================================
##
##     CORRELATIONS CRITICAL VALUES
##
##-----------------------------------------------------
library(psych)

y = fits$qualtrics$m4@response
X = fits$qualtrics$m4@effects
rm(fits); gc()
cr = psych::corr.test(X, alpha = 0.05, adjust = 'none')

n <- 815351
r <- seq(-.000005,.000005,.00000005)
tstat <- r * sqrt(n-2)/sqrt(1-r^2)
se <- sqrt( (1-r^2) / (n-2))
p <- pt(tstat/se, df=n-2, lower.tail = F)
plot(x=r, y=p, type='o', ylab='probability', xlab='correlation'); abline(h=0.05)

colnames(cr$p) <- 1:ncol(cr$p)
colnames(cr$t) <- 1:ncol(cr$t)
colnames(cr$r) <- 1:ncol(cr$r)

rr <- cr$r[5,2]
tt <- cr$t[5,2]
pp <- cr$p[5,2]
se <- cr$se[5,2]
  
r <- cr$r[5,2]

## CRITICAL VALUES OF CORRELATIONS
tstat <- r * sqrt(n-2)/sqrt(1-r^2)
pval <- pt(tstat, df=n-2, lower.tail = F) * 2

rcrit <- function(t,n)t/sqrt(n-2+t^2)
alpha <- 0.05
tcrit <- qt(alpha/2, df=n-2, lower.tail = F)
print(sprintf('|r| > %.5f significant at p<0.05.',rcrit(tcrit, n)))

abs(cr$r) > rcrit(tcrit, n)

## test all equal comparing corelation vs pvals to their critical values
(abs(cr$r) > rcrit(tcrit, n))  == (cr$p < alpha)



cr$r[10,17]
cr$r[ 4,14]
cr$r[ 9,14]
cr$r[13,14]

