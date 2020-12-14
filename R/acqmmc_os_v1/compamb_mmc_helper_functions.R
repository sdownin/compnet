####################################
##
## Competitive Ambidexterity
## SAOM analysis helper functions
##
##
####################################

###
## Compute HHI 
###
hhi <- function(x, na.rm=FALSE){
  if (na.rm) {
    x <- x[ which(!is.nan(x) & !is.na(x)) ]
  }
  z <- 100 * x / sum(x)
  return(sum(z^2))
}

##
# Plot from adjacency matrix
##
arrPlot <- function(arr, t=1, min.deg=0, ...)
{
  top <- ifelse('main' %in% names(list(...)), 2, .1)
  par(mar=c(.1,.1,top,.1))
  gx=graph.adjacency(arr[,,t])
  plot(induced.subgraph(gx,which(igraph::degree(gx)>=min.deg)), 
       edge.arrow.width=.2, edge.arrow.size=.1, 
       vertex.label.cex=.5, vertex.size=5, 
       ...)
}

###
## Complexity measure of actions 
##  @see Yu, Subramaniam, & Cannella, 2009
###
complexity <- function(x, scale=FALSE) {
  sumx <- sum(x)
  sq.prop <- sapply(x, function(k) (k/sumx)^2 )
  sum.sq.prop <- sum(sq.prop)
  y <- ifelse(sum.sq.prop==0, 0, 1/sum.sq.prop)
  if (scale)
    return(y / length(x))
  return(y)
}

##
#  Check Siena Model Converged
#   - max overall t ratio < 0.25
#   - all t-ratio < 0.1
#  @return bool
##
checkSienaConv <- function(res, t.lim=0.1,max.lim=0.25) {
  tmax <- res$tconv.max
  ts <- abs(res$tconv)
  ck.tmax <- tmax < max.lim
  ck.ts <- all(ts < t.lim)
  ck <- all(ck.tmax, ck.ts)
  cat(sprintf('\nCONVERGED:  %s\n  tconv.max (%.3f) < 0.25  %s\n  all t (max %.3f) < 0.10  %s\n\n',
              ck, tmax, ck.tmax, max(ts), ck.ts))
  return(ck)
}


###
## CrunchBase category Cosine Similarity
##
cbCatCosSim <- function(df) 
{
  all_cat_vec <- sort(unique(unlist(strsplit(df$category_list, '[|]'))))
  
  ## CATEGORY-FIRM MATRIX  [M,N] : M cateories, N firms
  firm_cats <- df$category_list
  cfm <- unname(sapply(firm_cats, function(x){
    as.integer(sapply(all_cat_vec,function(cat)grepl(cat,x)))
  }))
  ## FIRM-CATEGORY MATRIX [N,M] :  N firms, M categories
  fcm <- t(cfm)
  
  ## COMPUTE Cosine Similarity:
  ## =  u.v / |u||v|
  ## Firm-Firm COVARIANCE MATRIX [[u1.v1],[u2.v2], ...]
  X <- fcm %*% t(fcm)
  ## Firm Norms  |u|;  and |v| = |u|
  u <- apply(fcm, MARGIN = 1, FUN = function(x) sqrt(sum(x^2)) )
  ## NORM-NORM MATRIX   [ [|u1||v1|, |u1||v2|, ...], [|u2||v1|, |u2||v2|, ...], ...]
  UV <- outer(u, u, '*')
  ## COSINE SIMILARITY MATRIX equals elementwise division of 
  ##  covariance matrix by the norms matrix
  ##  = X / UV
  sim <- X / UV
  
  ## replace NULL, NA, NaN, diags with zero
  sim[is.null(sim)] <- 0
  sim[is.nan(sim)] <- 0
  sim[is.na(sim)] <- 0
  diag(sim) <- 0
  
  ## RETURN
  return(sim)
}

##
#
##
pngPlot <- function(x, filename, height = 4, width = 6.5, units = 'in', res = 300)
{
  png(filename,height = height, width = width, units = units, res =res )
  plot(x)
  dev.off()
}

##
#
##
getPvalStars <- function(p) {
  if (is.na(p)|is.nan(p))return('   ')
  if(p < 0.001)return('***')
  if(p < 0.01) return('** ')
  if(p < 0.05) return('*  ')
  return('   ')
}

##
#
##
saomResDf <- function(res, digits=3) 
{
  df <- data.frame(
    DV=res$effects$name,
    Effect=res$effects$effectName,
    Type=res$effects$type,
    Est=round(res$theta, digits = digits),
    se=round(res$se, digits = digits),
    t=round(res$theta/res$se, digits = digits),
    p=round(pt(abs(res$theta/res$se), df=Inf, lower.tail = F) * 2, digits = digits),
    stringsAsFactors = F
  )
  idx.rate <- grep('rate',df$Effect,T,T)
  df$t[idx.rate] <- NA
  df$p[idx.rate] <- NA
  return(df)
}

##
#
##
cbindDfList <- function(dfList) {
  df <- dfList[[1]]
  if (length(dfList) > 1) {
    for (i in 2:length(dfList)) {
      df <- cbind(df, dfList[[i]])
    }
  }
  return(df)
}


##
#  Create SAOM regression comparison Table
##
saomTable <- function(resList, file=NA, nameMap=NA, digits=3, drop.p.col=TRUE, drop.dv.col=TRUE)
{
  if (class(resList) != 'list') {
    resList <- list(resList)
  }
  behNames <- c()
  netNames <- c()
  for (res in resList) {
    behNames <- c(behNames, names(res$f$Data1$Behaviors))
    netNames <- c(netNames, names(res$f$Data1$nets))
  }
  dvNames <- unique(c(behNames, netNames))
  
  nameList <- list()
  for (res in resList) {
    resDvdf <- data.frame()
    for (dv in dvNames) {
      .effDf <- as.data.frame(res$effects)
      effTypeNames <- .effDf[grep(dv,res$effects$name,T,T),c('effectName','type')]
      nameList[[dv]] <- unique(rbind(resDvdf, effTypeNames))
    }
  }
  
  dfl <- list()
  for (res in resList) {
    for (dv in dvNames) {
      if (dv %in% res$effects$name) {
        .df <- saomResDf(res, digits=digits)
        dfl[[dv]][[ length(dfl[[dv]]) + 1 ]] <- .df[which(.df$DV==dv),]
      }
    }
  }
  
  mod.cols <- c('Est','se','p')
  
  tdf <- data.frame(stringsAsFactors = F)
  for (dv in names(nameList)) {
    hasRowname <- FALSE
    for (rowi in 1:nrow(nameList[[dv]])) { ## effect row
      eff <- nameList[[dv]][rowi,]
      effRow <- list()
      for (modDf in dfl[[dv]] ) { ## model dataframe in DV group
        effId <- which(modDf$Effect == eff$effectName & modDf$Type == eff$type)
        if (length(effId) > 0) {
          effRow[[length(effRow)+1]] <- modDf[effId,mod.cols]
        } else {
          .nadf <- data.frame()
          for (col in mod.cols) .nadf[1,col] <- NA
          effRow[[length(effRow)+1]] <- .nadf
        }
      }
      effRowDf <- cbind(data.frame(DV=dv,Effect=eff$effectName, Type=eff$type, stringsAsFactors = F), cbindDfList(effRow))
      if (!hasRowname) {
        effRowDf <- rbind(effRowDf, effRowDf)
        effRowDf[1,]  <- c('', sprintf('Dynamics: %s', dv), rep(NA, ncol(effRowDf)-2))
        hasRowname <- TRUE
      }
      tdf <- rbind(tdf, effRowDf)
    }
  }
  
  # move rate rows to end
  .tmp.rate.row <- tdf[1,]
  .tmp.rate.row$Effect <- 'Rate Parameters'
  rate.idx <- which(tdf$Type=='rate')
  tdf <- rbind(tdf[-rate.idx,], .tmp.rate.row, tdf[rate.idx, ])
  
  obs <- c()
  ns <- c()
  conv <- c()
  convt <- c()
  iter <- c()
  for (res in resList) {
    obs <-c(obs, res$observations)
    ns <- c(ns, attributes(res$f$Data1$nets[[1]][[1]][[1]])$nActors)
    conv <- c(conv, res$tconv.max)
    convt<- c(convt, max(abs(res$tconv)))
    iter <- c(iter, res$n)
  }  
  
  # est idx
  idx.est <- which(names(tdf)%in% 'Est')
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Time Periods'
  tdf[nrow(tdf), idx.est] <-  obs
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Num. Firms'
  tdf[nrow(tdf), idx.est] <-  ns
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Largest converg. t ratio'
  tdf[nrow(tdf), idx.est] <-  round(convt, digits = digits)
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Overall max. converg. ratio'
  tdf[nrow(tdf), idx.est] <-  round(conv, digits = digits)
  #
  tdf[nrow(tdf)+1, ] <- NA
  tdf[nrow(tdf), 'Effect'] <- 'Iterations'
  tdf[nrow(tdf), idx.est] <-  iter
  
  idx.se  <- grep('^se\\.{0,1}',names(tdf),T,T)
  idx.p <- grep('^p\\.{0,}', names(tdf),T,T)
  for (i in 1:length(idx.se)) {
    sei <- idx.se[i]
    pi <- idx.p[i]
    tdf[,sei] <- apply(tdf[,c(sei,pi)],1,function(x) {
      se <- x[1]
      p <- x[2]
      spfstr <- sprintf('(%s%s.%sf)%s','%',digits+2,digits,getPvalStars(p))
      ifelse(is.na(se)|se=='NA','',sprintf(spfstr,as.numeric(se)))
    })
  }
  
  ## add Type to name (not eval or rate)
  idx.nonrate <- which(tdf$Type %in% c('endow','creation'))
  tdf$Effect[idx.nonrate] <- apply(tdf[idx.nonrate,c('Effect','Type')],1,function(x){
    sprintf('%s: %s',x[2],x[1])
  })
  tdf <- tdf[,which(names(tdf) != 'Type')]
  
  ## name mapping for effects
  if (!any(is.na(nameMap))) 
  {
    ord <- c()
    for (eff in names(nameMap)) {
      idx.nm <- which(tdf$Effect == eff)
      ord <- c(ord, idx.nm)
      if (length(idx.nm) > 0) {
        tdf$Effect[idx.nm] <- nameMap[[eff]]
      }
    }
    tdf <- rbind(tdf[ord,], tdf[-ord,])
  }
  
  if (drop.dv.col) {
    tdf <- tdf[,-1]
  }
  if (drop.p.col) {
    idx.p.col <- grep('^p[\\.\\d]{0,}',names(tdf),T,T)
    tdf <- tdf[, -idx.p.col]
  }
  
  if (!is.na(file)) {
    write.csv(tdf, file = file, na = "", row.names = F)
  }
  
  return(tdf)
  
}


##
## Network stability measure (Jaccard index)
##
netJaccard <- function(m1, m2)
{
  ## num. maintained
  n11 <- sum(m1 * m2) ## element-wise multiplication
  ##  num. dropped or added
  n10.n01 <- sum( (m2 - m1) != 0 )
  ## jaccard index
  return( n11 / (n11 + n10.n01) )
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





##
# Gets the count of actions of given type for given firms in RavenPack data
# @param [dataframe] rvnsub   RavenPack subset by period (year)
# @param [string[n]] names    CrunchBase company_name_unique vector for network vertices
# @param [dataframe] cbrpmap  CB-to-RP name (entity ID) mapping
# @param [string] var_rp      Name of action type in RavenPack dataframe 
# @return [integer[n]]        The action counts vector
##
getRpActionCount <- function(rvnsub, names, cbrpmap, var_rp)
{
  n <- length(names)
  x <- rep(0,n)
  
  if (length(rvnsub)==0 | nrow(rvnsub)==0)
    return(x)
  
  ## network vertex indices of firms in RavenPack data
  vids <- which(names %in% cbrpmap$company_name_unique)
  for (i in 1:length(vids)) 
  {
    actids <- which(
      rvnsub$rp_entity_id==cbrpmap$rp_entity_id[i] & 
        rvnsub$action_category==var_rp
    )
    v <- vids[i]  ## the v'th firm in network is the i'th RP entity
    x[v] <- length(actids) ## count of actions is the length of the action ids vector
  }
  
  return(x)
}

##
# Gets each vertex's sum of neighbors' values of given attribute (attr)
# @param [igraph] gt    Igraph object (for period t)
# @param [integer] vid  Vertex index of focal firm to fetch neighbors' attributes
# @param [string] attr  Name of vertex attribute
# @return [integer[m]]  Vector of neighbors attribute values (for m neighbors)
##
getNeighborsAttrSum <- function(gt, attr) {
  adjmat <- igraph::get.adjacency(gt)
  attr <- igraph::get.vertex.attribute(gt, attr)
  return(as.numeric(adjmat %*% attr))
}
# getNeighborsAttr <- function(gt, vid, attr) {
#   nvids <- as.integer(igraph::neighbors(gt, vid))
#   return(igraph::get.vertex.attribute(gt, attr, nvids))
# }


