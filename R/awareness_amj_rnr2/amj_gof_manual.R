##
##
##  GOODNESS OF FIT
##  BATCH RUNS
##
##

library(btergm)
library(MASS)


work_dir <- 'C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\competition networks\\compnet2'
data_dir <- file.path(work_dir, 'firm_nets_rnr2\\firmctrl')
results_dir <- file.path(work_dir,'amj_rnr2_results\\firmctrl')

listsum <- function(l){
  x <- as.matrix(l[[1]])
  for (i in 2:length(l))
    x <- x + as.matrix(l[[i]])
  return(x)
}


# function which reduces a statistic x nsim matrix and computes summary stats
# input: two matrices of a certain type of statistics (simulated and observed)
# goal: get rid of empty rows at the end, e.g., where dsp(37) or so is usually
# not observed; return value: object containing the summary statistics
reduce.matrix <- function(sim, obs) {
  
  numsim <- ncol(as.matrix(sim))
  numobs <- ncol(as.matrix(obs))
  xlabels <- rownames(obs)
  
  # if geodist statistic: put aside the last 'Inf' row
  if (is.null(rownames(sim)) || rownames(sim)[nrow(sim)] == "Inf") {
    geo <- TRUE
    inf.sim <- sim[nrow(sim), ]  # put aside this row for now and reuse later
    sim <- sim[-nrow(sim), ]
    if (class(obs) == "matrix") {
      inf.obs <- obs[nrow(obs), ]
      obs <- matrix(obs[-nrow(obs), ], ncol = numobs)
    } else {
      inf.obs <- obs[length(obs)]
      obs <- matrix(obs[-length(obs)], ncol = numobs)
    }
  } else {
    geo <- FALSE
  }
  
  # find first empty row for simulation matrix
  sim.rs <- rowSums(sim)
  sim.remove <- length(sim.rs)  # at which row index can removal start?
  for (i in (length(sim.rs) - 1):1) {
    if (sim.rs[i] == 0 && sim.remove == (i + 1)) {
      sim.remove <- i  # remember which is the first empty row
    }
  }
  
  if (class(obs) != "matrix") {  # one network is compared
    obs.remove <- length(obs)
    for (i in (length(obs) - 1):1) {
      if (obs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  } else {  # several networks are compared
    obs.rs <- rowSums(obs)
    obs.remove <- length(obs.rs)
    for (i in (length(obs.rs) - 1):1) {
      if (obs.rs[i] == 0 && obs.remove == (i + 1)) {
        obs.remove <- i
      }
    }
  }
  rem <- max(c(obs.remove, sim.remove), na.rm = TRUE)  # which one is longer?
  
  # remove unnecessary rows
  if (class(obs) != "matrix") {  # get rid of empty observations or rows of obs
    obs <- matrix(obs[-(rem:length(obs))], ncol = numobs)
  } else {
    obs <- matrix(obs[-(rem:nrow(obs)), ], ncol = numobs)
  }
  sim <- matrix(sim[-(rem:nrow(sim)), ], ncol = numsim)  # same for sim stats
  
  if (nrow(obs) < rem) {
    for (i in (nrow(obs) + 1):rem) {
      obs <- rbind(obs, rep(0, ncol(obs)))
    }
  }
  if (nrow(sim) < rem) {
    for (i in (nrow(sim) + 1):rem) {
      sim <- rbind(sim, rep(0, ncol(sim)))
    }
  }
  
  # for geodist, add Inf row again
  if (geo == TRUE) {
    sim <- rbind(sim, inf.sim)
    obs <- rbind(obs, inf.obs)
    rownames(sim) <- c(1:(nrow(sim) - 1), "Inf")
  }
  
  # create final object which will contain the raw simulations + the comparison
  reducedobject <- list()
  reducedobject$sim <- as.data.frame(sim)
  
  rownames(sim) <- NULL
  rownames(obs) <- NULL
  
  # compute means, p values, etc. and put them in a data frame
  x <- matrix()
  if (ncol(obs) == 1 || ncol(sim) == 1) {  # compute all the summary statistics
    x.obs <- obs
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    zval <- (x.mean - x.obs) / sd(x.mean)
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    x.pval <- pval
    x <- data.frame(x.obs, x.mean, x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs", "sim: mean", "median", "min", "max", "Pr(>z)")
  } else {  # for several target networks, compute distribution
    x.obs.mean <- apply(obs, 1, mean)
    x.obs.min <- apply(obs, 1, min)
    x.obs.max <- apply(obs, 1, max)
    x.obs.median <- apply(obs, 1, median)
    x.mean <- apply(sim, 1, mean)
    x.min <- apply(sim, 1, min)
    x.max <- apply(sim, 1, max)
    x.median <- apply(sim, 1, median)
    x.pval <- numeric()
    for (i in 1:nrow(sim)) {
      tryCatch(
        expr = {
          x.pval[i] <- t.test(obs[i, ], sim[i, ])$p.value  # compare group means
        }, 
        error = function(e) {
          x.pval[i] <- 1  # if both are 0, a t.test cannot be computed...
        }, 
        finally = {}
      )
      
      if (is.nan(x.pval[i])) {  # geodist contains "Inf"
        x.pval[i] <- 1
      }
    }
    x <- data.frame(x.obs.mean, x.obs.median, x.obs.min, x.obs.max, x.mean, 
                    x.median, x.min, x.max, x.pval)
    colnames(x) <- c("obs: mean", "median", "min", "max", 
                     "sim: mean", "median", "min", "max", "Pr(>z)")
  }
  if (geo == TRUE) {
    rownames(x) <- c(xlabels[c(1:(nrow(x) - 1))], "Inf")
  } else {
    rownames(x) <- xlabels[1:nrow(x)]
  }
  rownames(reducedobject$sim) <- rownames(x)
  
  reducedobject$comparison <- as.data.frame(x)
  
  alpha <- 0.05
  reducedobject$comparison$. <- sapply(reducedobject$comparison$`Pr(>z)`, function(x){
    return(ifelse(x < alpha, '* ', ''))
  })
  
  return(reducedobject)
}









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
# cpa <- lapply(nets, function(net) as.matrix(net %n% 'coop') + as.matrix(net %n% 'coop_past') )
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




# SAVE BTERGM GOF
# gf <- gof(fit, nsim=2, statistics=c(deg,dsp,esp,geodesic))
# gf.file.path <- file.path(results_dir, sprintf('gof_winlocal_%s_pd%s_R%s_%s_d%s.rds',firm_i,nPeriod,R,m_x,d))
# saveRDS(gf, file = gf.file.path)
# plot(gf)


















##========================================================================
###--------------------------------------
##
##  BATCH AUXILIARY STATISTICS
##
##---------------------------------------
## USER SETTINGS
nsim.total <- 2000
nsim.batch <- 2  ## per period each batch

# gof_%s_pd%s_R%s_%s_nsim%s_WINLOCAL.rds',firm_i,nPeriod,R,m_x,nsim
gof_batch_filename <- sprintf('gof_batch_sim_stats_%s_pd%s_R%s_%s.rds',firm_i,nPeriod,R,m_x)
gof_batch_file <- file.path(results_dir, gof_batch_filename)

## load or init gof list
if (file.exists(gof_batch_file)) {
  l.gof <- readRDS(gof_batch_file)
} else {
  ## init if not exists, else load from file to start
  l.gof <- list()
}

## set run stats (remaining batches to run to reach total)
stats <- list(deg=deg,dsp=dsp,esp=esp,geodesic=geodesic)
len <- length(nets)-1
periods <- 1:len
#
num.batches <- nsim.total / nsim.batch

# MAIN LOOP
for (batch in 1:num.batches)
{
  
  cat(sprintf('\nbatch %s; nsim = %s:  t = ',batch,nsim.batch))
  stat.mats <- list()
  
  ## simulation loop for each time period network
  for (t in periods) {
    cat(sprintf(' %s ',t))
    sims <- simulate(fit, nsim=nsim.batch, index=t, verbose = F, statsonly=FALSE)
    #
    for (stat_z_name in names(stats)) {
      stat_z <- stats[[stat_z_name]]
      stat.mats[[stat_z_name]][[t]] <- sapply(sims, function(sim) stat_z(as.matrix(sim)))
    }
  }
  
  
  for (stat_z_name in names(stats)) {
    # x <- Reduce('+',  stat.mats[[stat_z_name]])
    x <- stat.mats[[stat_z_name]]
    
    if (length(l.gof[[stat_z_name]])<1) {
      l.gof[[stat_z_name]] <- list()
      l.gof[[stat_z_name]][[1]] <-  x
    } else {
      .len <- length(l.gof[[stat_z_name]])
      l.gof[[stat_z_name]][[(.len+1)]] <- x
    }
  }

  ## Save progress
  saveRDS(l.gof, file = gof_batch_file)
  
}


##===================== END =================================================






















###-----------------------------------
##
##  CHECK GOF RESULTS SAVED
##  MAHALANOBIS DISTANCE
##
##------------------------------------
prev_gof_batch_file <- 'gof_batch_sim_stats_qualtrics_pd11_R2000_m4 - Copy.rds'
l.gof <- readRDS(file.path(results_dir, prev_gof_batch_file))


## stats
stats <- list(deg=deg,dsp=dsp,esp=esp,geodesic=geodesic)

## test list
test <- list(
  deg=list(bin=c(), p=c()),
  dsp=list(bin=c(), p=c()),
  esp=list(bin=c(), p=c()),
  degree=list(bin=c(), p=c())
)

## l.gof[[stat]][[batch]][[period]] = simulations.matrix

## for each statistic
for (stat_z_name in names(stats)) {
  stat_z <- stats[[stat_z_name]]u

  
  
  mu <- rowMeans(l.gof[[stat_z_name]])
  
  ## observed distance
  obs.sub <- sapply(nets[2:length(nets)], function(net) stat_z(as.matrix(net)))
  obs <- rowSums(obs.sub)
  
  ## fix length mismatch
  if (length(obs) > length(mu)) {
    obs <- obs[1:(length(obs)-1)]
  }
  # cat(sprintf('length obs %s, mu %s', length(obs), length(mu)))
  
  sigma <- cov( t(l.gof[[stat_z_name]]) )
  sigma.inv <- ginv(sigma)
  d.obs <- mahalanobis(obs, mu, sigma.inv, inverted = TRUE)
  
  ## for each simulation
  for (i in 1:ncol(l.gof[[stat_z_name]])) {
    sim <- as.numeric( l.gof[[stat_z_name]][,i] )
    d.sim <- mahalanobis(sim, mu, sigma.inv, inverted = TRUE)
    #
    bin <- as.integer( d.sim >= d.obs )
    test[[stat_z_name]]$bin <- c(test[[stat_z_name]]$bin, bin)
    #
    p <- sum(test[[stat_z_name]]$bin) / length(test[[stat_z_name]]$bin) 
    test[[stat_z_name]]$P <-   c(test[[stat_z_name]]$P , p)   
  }
  
  
}
































##---------------------------------------
##  GOODNESS OF FIT
##---------------------------------------
## USER SETTINGS
nsim.total <- 100
nsim.batch <- 10  ## per period each batch
target <- nets[2:length(nets)]  ## skip lagged year & first 3 years

# gof_%s_pd%s_R%s_%s_nsim%s_WINLOCAL.rds',firm_i,nPeriod,R,m_x,nsim
gof_batch_filename <- sprintf('gof_batch_%s_pd%s_R%s_%s_simtot%s.rds',firm_i,nPeriod,R,m_x, nsim.total)
gof_batch_file <- file.path(results_dir, gof_batch_filename)

#
len <- length(nets)-1
obs.indices <- 1:len
stats <- list(deg=deg,dsp=dsp,esp=esp,geodesic=geodesic)

## load or init gof list
if (file.exists(gof_batch_file)) {
  l.gof <- readRDS(gof_batch_file)
} else {
  
  ## init if not exists, else load from file to start
  l.gof <- list(
    obs=list(),
    obs.sum=list(),
    sim=list(),
    sim.sum=list(),
    sim.mean=list(),
    test=list(),
    p=list(),
    nsim.batches=c(),  ## batch size record
    nsim.cumulative=0, ## completed progress (sum of batches)
    nsim.total=0       ## total to run (set by user)
  )
  
  for (stat_z_name in names(stats)) {
    l.gof$sim[[stat_z_name]] <- list()
    l.gof$sim.sum[[stat_z_name]] <- list()
    l.gof$sim.mean[[stat_z_name]] <- list()
    #
    for (index in obs.indices) {
      l.gof$obs[[stat_z_name]][[index]] <- stat_z(as.matrix(target[[index]]))  ### index + 1
    }
    l.gof$obs.sum[[stat_z_name]] = Reduce('+', l.gof$obs[[stat_z_name]])
  }
  
}

## set run stats (remaining batches to run to reach total)
l.gof$nsim.total <- nsim.total
sim.todo <- l.gof$nsim.total - l.gof$nsim.cumulative
num.batches <- round(sim.todo / (nsim.batch*len))


# MAIN LOOP
for (batch in 1:num.batches)
{
  cat(sprintf('\nbatch %s; nsim = %s:  t = ',batch,nsim.batch))
  
  ## simulation loop for each time period network
  ## minimize memory by only caching one period of network simulations before overwriting
  for (index in obs.indices) {  ## index of observed network (skipping first lag)
    cat(sprintf(' %s ',index))
    sims <- simulate(fit, nsim=nsim.batch, index=index, verbose = F, statsonly=FALSE)
    
    ## compute statistics
    for (stat_z_name in names(stats)) {
      stat_z <- stats[[stat_z_name]]
      cat(stat_z_name)
      sim.z.mat <- sapply(sims, function(sim) stat_z(as.matrix(sim)))
      
      if (length(l.gof$sim[[stat_z_name]]) < index) {
        l.gof$sim[[stat_z_name]][[index]] <-  sim.z.mat 
      } else {
        l.gof$sim[[stat_z_name]][[index]] <- cbind( l.gof$sim[[stat_z_name]][[index]], sim.z.mat ) 
      }
    }
    zdim <- dim(l.gof$sim[[stat_z_name]][[index]] )
    cat(sprintf('dim (%s,%s)',dim(zdim)[1],dim(zdim)[2]))
  }
  
  ## for each statistic
  for (stat_z_name in names(stats)) {
    l.gof$sim.sum[[stat_z_name]] <-  Reduce('+', l.gof$sim[[stat_z_name]])
    l.gof$sim.mean[[stat_z_name]] <- rowMeans(l.gof$sim.sum[[stat_z_name]])
    
    ## indices to keep in test
    zlen <- length(l.gof$sim.mean[[stat_z_name]])
    zidx <- 1:round(zlen * test.prop, 0)
    
    ## observed distance
    obs <- as.numeric( l.gof$obs.sum[[stat_z_name]] )
    mu <- rowMeans(l.gof$sim.sum[[stat_z_name]])
    
    ## fix length mismatch
    if (length(obs) > length(mu)) {
      obs <- obs[1:(length(obs)-1)]
    }
    # cat(sprintf('length obs %s, mu %s', length(obs), length(mu)))
    
    sigma <- cov( t(l.gof$sim.sum[[stat_z_name]]) )
    sigma.inv <- ginv(sigma)
    d.obs <- mahalanobis(obs, mu, sigma.inv, inverted = TRUE)
    
    ## for each simulation
    for (i in 1:ncol(l.gof$sim.sum[[stat_z_name]])) {
      sim <- as.numeric( l.gof$sim.sum[[stat_z_name]][,i] )
      d.sim <- mahalanobis(sim, mu, sigma.inv, inverted = TRUE)
      test <- as.integer( d.sim >= d.obs )
      l.gof$test[[stat_z_name]] <- c(l.gof$test[[stat_z_name]], test)
      l.gof$p[[stat_z_name]] <- sum(l.gof$test[[stat_z_name]]) / length(l.gof$test[[stat_z_name]])      
    }
    
    
  }
  
  print(l.gof$test)
  print(l.gof$p)
  
  ## Save progress
  saveRDS(l.gof, file = gof_batch_file)
}



























##---------------------------------------
##  BATCH GOODNESS OF FIT FROM btergm package
##---------------------------------------
## USER SETTINGS
nsim.total <- 2000
nsim.batch <- 2  ## per period each batch
target <- fit@data$networks[4:11]  ## skip lagged year & first 3 years

# gof_%s_pd%s_R%s_%s_nsim%s_WINLOCAL.rds',firm_i,nPeriod,R,m_x,nsim
gof_batch_filename <- sprintf('gof_batch_%s_pd%s_R%s_%s_simtot%s.rds',firm_i,nPeriod,R,m_x, nsim.total)
gof_batch_file <- file.path(results_dir, gof_batch_filename)

## load or init gof list
if (file.exists(gof_batch_file)) {
  l.gof <- readRDS(gof_batch_file)
} else {
  ## init if not exists, else load from file to start
  l.gof <- list(
    data=list(
      sim=list(), 
      obs=list()
    ),
    summary=list(
      sim=list(),
      table=list()
    ),
    nsim.batches=c(),  ## batch size record
    nsim.cumulative=0, ## completed progress (sum of batches)
    nsim.total=0       ## total to run (set by user)
  )
}

## set run stats (remaining batches to run to reach total)
stats <- list(deg=deg,dsp=dsp,esp=esp,geodesic=geodesic)
len <- length(nets)-1
periods <- 1:len
l.gof$nsim.total <- nsim.total
sim.todo <- l.gof$nsim.total - l.gof$nsim.cumulative
num.batches <- round(sim.todo / (nsim.batch*len))

# MAIN LOOP
for (batch in 1:num.batches)
{
  sims <- list()
  cat(sprintf('\nbatch %s; nsim = %s:  t = ',batch,nsim.batch))
  
  ## simulation loop for each time period network
  for (t in periods) {
    cat(sprintf(' %s ',t))
    .tmp <- simulate(fit, nsim=nsim.batch, index=t, verbose = F, statsonly=FALSE)
    sims <- c(as.list(sims), as.list(.tmp))
  }
  
  ## compute statistics
  for (stat_z_name in names(stats)) {
    stat_z <- stats[[stat_z_name]]
    sim.z.mat <- sapply(sims, function(sim) stat_z(as.matrix(sim)))
    obs.z.mat <- sapply(target, function(net) stat_z(as.matrix(net)))
    ##
    l.gof$data$sim[[stat_z_name]] <- cbind(l.gof$data$sim[[stat_z_name]], sim.z.mat)
    l.gof$data$obs[[stat_z_name]] <- cbind(l.gof$data$obs[[stat_z_name]], obs.z.mat)
    ##
    reduced <-  reduce.matrix(
        l.gof$data$sim[[stat_z_name]], 
        l.gof$data$obs[[stat_z_name]]
      ) 
    l.gof$summary$sim[[stat_z_name]] <- reduced$sim
    l.gof$summary$table[[stat_z_name]] <- reduced$comparison
  }
  
  ## Update progress
  l.gof$nsim.batches <- c(l.gof$nsim.batches, nsim.batch*len)
  l.gof$nsim.cumulative <- sum(l.gof$nsim.batches)
  
  ## Save progress
  saveRDS(l.gof, file = gof_batch_file)
}
















##-----------------------------------
## BATCH GOF from btergm package
##-----------------------------------








##
stat <- 'deg'
obs.mean <- l.gof$summary$table[[stat]]$`obs: mean`
sim.mean <- l.gof$summary$table[[stat]]$`sim: mean`
stat.len <- length(sim.mean)
stat.data <- t(l.gof$summary$sim[[stat]])
sigma <- cov(stat.data)
stat.idx <- 1:(stat.len-10)
md2 <- mahalanobis(x = sim.mean[stat.idx], 
                  center = obs.mean[stat.idx], 
                  cov = sigma[stat.idx,stat.idx])
md <- sqrt(md2)

pchisq(md, df=length(stat.idx)-1, lower.tail = FALSE) ## lower.tail=FALSE   P[X>x]

##  Mahalanobis Distance = 16.384; p =.173)







##-----------------------------------------------
##
##-----------------------------------------------

l <- btergm::tergmprepare(formula = getformula(fit), offset = fit@offset, verbose = T)
offset <- object@offset
form <- as.formula(l$form, env = environment())

target <- nets

# check and rearrange target network(s)
if (is.null(target)) {
  if (verbose == TRUE) {
    message(paste("\nNo 'target' network(s) provided. Using networks on the",
                  "left-hand side of the model formula as observed networks.\n"))
  }
  target <- l$networks
} else if (class(target) == "network" || class(target) == "matrix") {
  target <- list(target)
  if (verbose == TRUE) {
    message("\nOne observed ('target') network was provided.\n")
  }
} else if (class(target) == "list") {
  if (verbose == TRUE) {
    message(paste("\n", length(target), "observed ('target') networks were",
                  "provided.\n"))
  }
} else {
  stop("'target' must be a network, matrix, or list of matrices or networks.")
}

# extract coefficients from object
if (offset == TRUE) {
  coefs <- c(coef(object), -Inf)  # -Inf for offset matrix
} else {
  coefs <- coef(object)
}

# adjust formula at each step, and simulate networks
sim <- list()
degen <- list()
for (index in 1:l$time.steps) {
  i <- index  # index 'i' is used in formula construction!
  # simulations for statnet-style and rocpr GOF
  if (verbose == TRUE) {
    if ("btergm" %in% class(object) || "mtergm" %in% class(object)) {
      f.i <- gsub("\\[\\[i\\]\\]", paste0("[[", index, "]]"), 
                  paste(deparse(form), collapse = ""))
      f.i <- gsub("\\s+", " ", f.i)
      if ("btergm" %in% class(object)) {
        f.i <- gsub("^networks", l$lhs.original, f.i)
      }
    } else if ("ergm" %in% class(object)) {
      f.i <- paste(deparse(formula), collapse = "")
      f.i <- gsub("\\s+", " ", f.i)
    } else {
      stop(paste("Unknown object type:", class(object)))
    }
    message(paste("Simulating", nsim, 
                  "networks from the following formula:\n", f.i, "\n"))
  }
  tryCatch(
    {
      sim[[index]] <- simulate(form,
                               nsim = nsim,
                               coef = coefs,
                               constraints = ~ .,
                               control = control.simulate.formula(MCMC.interval = MCMC.interval,
                                                                  MCMC.burnin = MCMC.burnin))
    }, error = function(cond) {
      sim[[index]] <- NULL
      warning("There was a problem with the simulation at t = ", index, 
              ". Sometimes this may be due to node attributes being ",
              "incomplete at some time steps. For example, not all levels ", 
              "of a nodefactor term may be present at all time steps. ", 
              "Original error message: ", cond)
    }
  )
  
}
