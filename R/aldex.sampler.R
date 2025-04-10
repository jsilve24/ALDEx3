##' Get chunks of estimates of the absolute abundances, composition, and scale
##' given a scale model (internal only).
##'
##' @param Y an (D x N) matrix of sequence count, N is number of samples, D
##'   is number of taxa or genes
##' @param X A model matrix of dimension P x N (P is number of linear model
##'   covariates).
##' @param chunk.size The number of monte carlo samples of the absolute
##'   abudnances, composition, and scale to return.
##' @param scale the scale model, can be a function or an N x nsample matrix.
##' @param scale.args A list of parameters to pass to the scale model
##' @param return.pars Can include logWperp and logWpara to return the scale
##'   and composition samples, respetively
##' @return A list with the absolute abundance estimates logW, and optionally
##'   also logWperp and logWpara (see return.pars)
##' @author Justin Silverman, Kyle McGovern
aldex.chunk <- function(Y, X, chunk.size, scale,
                        scale.args, return.pars) {
  N <- ncol(Y)
  D <- nrow(Y)

  ## dirichlet sample
  logWpara <- array(NA, c(D, N, chunk.size))
  for (i in 1:N) {
    logWpara[,i,] <- log2(rDirichlet(chunk.size, Y[,i]+0.5))
  }

  ## sample from scale model
  if (is.null(scale)){
    stop("You probably want a scale model :)")
  } else if (is.matrix(scale)) {
    stopifnot(dim(scale)==c(N, chunk.size))
    logWperp <- scale
  } else if (is.function(scale)) {
    req.args <- formalArgs(scale)
    scale.args <- scale.args[intersect(names(scale.args), req.args)]
    if ("logWpara" %in% req.args) scale.args$logWpara <- logWpara
    if ("X" %in% req.args) scale.args$X <- X
    if ("Y" %in% req.args) scale.args$Y <- Y
    logWperp <- do.call(scale, scale.args)
    ## some error checking to make sure whats returned has right dimensions
  }

  ## compute scaled abundances (W)
  logW <- sweep(logWpara, c(2,3), logWperp, FUN=`+`)

  ## memory management
  if (!("logWperp" %in% return.pars)) rm(logWperp)
  if (!("logWpara" %in% return.pars)) rm(logWpara)

  ## Return Results
  out <- list()
  out$logW <- logW
  if ("logWperp" %in% return.pars) {
    out$logWperp <- logWperp
  }
  if ("logWpara" %in% return.pars) {
    out$logWpara <- logWpara
  }
  return(out)
}

##' Get chunks of estimates of the absolute abundances, composition, and scale
##' given a scale model (internal only).
##'
##' @param Y an (D x N) matrix of sequence count, N is number of samples, D
##'   is number of taxa or genes
##' @param X A model matrix of dimension P x N (P is number of linear model
##'   covariates).
##' @param nsample The total number of monte carlo samples to draw for aldex
##'   analysis.
##' @param streamsize memory footprint (approximate) at which to
##'   use streaming.
##' @return A vector of values to use as the chunk.size in aldex.chunk
##' @author Justin Silverman, Kyle McGovern
aldex.getchunksizes <- function(Y, X, nsample, streamsize) {
  N <- ncol(Y)
  D <- nrow(Y)

  ## calculate streaming threshold
  stream <- FALSE
  if (N*D*nsample*8/100000 > streamsize) stream <- TRUE

  nsample.remaining <- nsample
  iter <- 1
  if (stream) {
    nsample.local <- floor(streamsize*100000/(N*D*8))
    if (nsample.local < 1) stop("streamsize too small")
  } else {
    nsample.local <- nsample
  }

  chunk_sizes <- c()
  while (nsample.remaining > 0) {
    if(nsample.remaining<nsample.local) {
      chunk_sizes <- c(chunk_sizes, nsample.remaining)
    } else {
      chunk_sizes <- c(chunk_sizes, nsample.local)
    }
    nsample.remaining <- nsample.remaining - nsample.local
  }
  return(chunk_sizes)
}
