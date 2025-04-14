aldex.sampler <- function(Y, X, nsample, scale=NULL, scale.args, return.pars) {
  N <- ncol(Y)
  D <- nrow(Y)

  
  ## dirichlet sample
  logComp <- array(NA, c(D, N, nsample))
  for (i in 1:N) {
    logComp[,i,] <- log2(rDirichlet(nsample, Y[,i]+0.5))
  }

  ## sample from scale model
  if (is.null(scale)){
    stop("You probably want a scale model :)")
  } else if (is.matrix(scale)) {
    stopifnot(dim(scale)==c(N, nsample))
    logScale <- scale
  } else if (is.function(scale)) {
    req.args <- formalArgs(scale)
    scale.args <- scale.args[intersect(names(scale.args), req.args)]
    if ("logComp" %in% req.args) scale.args$logComp <- logComp
    if ("X" %in% req.args) scale.args$X <- X
    if ("Y" %in% req.args) scale.args$Y <- Y
    logScale <- do.call(scale, scale.args)
    ## some error checking to make sure whats returned has right dimensions
  } 

  out <- list()
  ## compute scaled abundances (W)
  out$logW <- sweep(logComp, c(2,3), logScale, FUN=`+`)

  if ("logScale" %in% return.pars) out$logScale <- logScale
  if ("logComp" %in% return.pars) out$logComp <- logComp

  return(out)
 }

##' Calculate p-values adjusting for changes in sign as described by Nixon
##' et al. (2024) in Beyond Normalization: Incorporating Scale Uncertainty
##' in Microbiome and Gene Expression Analysis (internal only)
##'
##' @param p.lower A P x D x S matrix for P covariates, D taxa/genes, and S
##'   monte carlo samples representing the lower tail p.values
##' @param p.upper A P x D x S matrix for P covariates, D taxa/genes, and S
##'   monte carlo samples representing the upper tail p.values
##' @param p.adjust.method An adjutment method for p.adjust
##' @return A list with P x D matrices with the non-adjusted and adjusted
##'   p-values.
##' @author Justin Silverman, Kyle McGovern
aldex.pvals <- function(p.lower, p.upper, p.adjust.method) {
  #### p-value calculations, accounting for sign changes ####
  p.lower.local <- array(pmin(1, 2 * p.lower), dim = dim(p.lower))
  p.upper.local <- array(pmin(1, 2 * p.upper), dim = dim(p.upper))

  p.lower.adj <- apply(p.lower.local, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.upper.adj <- apply(p.upper.local, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })

  p.lower.mean <- apply(p.lower.local, c(2,1), mean)
  p.upper.mean <- apply(p.upper.local, c(2,1), mean)
  rm(p.lower.local, p.upper.local)


  p.res <- c()
  for(col_i in 1:ncol(p.lower.mean)) {
    tmp_mat <- cbind(p.lower.mean[,col_i],
                     p.upper.mean[,col_i])
    p.res <- rbind(p.res, apply(tmp_mat, 1, min))
  }
  rm(p.lower.mean, p.upper.mean)

  p.lower.mean.adj <- apply(p.lower.adj, c(1,2), mean)
  p.upper.mean.adj <- apply(p.upper.adj, c(1,2), mean)
  p.adj.res <- c()
  for(col_i in 1:ncol(p.lower.mean.adj)) {
    tmp_mat <- cbind(p.lower.mean.adj[,col_i],
                     p.upper.mean.adj[,col_i])
    p.adj.res <- rbind(p.adj.res, apply(tmp_mat, 1, min))
  }
  rm(p.lower.mean.adj, p.upper.mean.adj)
  return(list(p.res=p.res, p.adj.res=p.adj.res))
 }
