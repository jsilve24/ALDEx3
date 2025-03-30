##' ALDEx3 Linear Models
##'
##' 
##' @title ALDEx3 Linear Modles
##' @param counts an (D x N) matrix of sequence count, N is number of samples, D
##'   is number of taxa or genes
##' @param X either a formula (in which case DATA must be non-null) or a model
##'   matrix of dimension P x N (P is number of linear model covariates). If a
##'   formula is passed it should not include the target Y, e.g., should simply
##'   be "~condition-1" (note the lack of the left hand side).
##' @param data a data frame for use with formula, must have N rows
##' @param nsample number of monte carlo replicates
##' @param GAMMA the scale model, can be a function or an N x nsample matrix. If
##'   a function is passed, it should take one argument (pi) which is a N x D x
##'   nsample array of log-transformed relative abundances. That function must
##'   in turn output an N x nsample matrix of scale samples on the log2 scale.
##' @param streamsize (default 8000) memory footprint (approximate) at which to
##'   use streaming. This should be thought of as the number of Mb for each
##'   streaming chunk. If D*N*nsample*8/1000000 is less than streamsize then no
##'   streaming will be performed. Note, to conserve memory, samples from the
##'   Dirichlet and scale models will not be returned if streaming is used.
##'   Streaming can be turned off by setting streamsize=Inf.
##' @param return.pars what results should be returned, see return section
##'   below.
##' @param p.adjust.method (default BH) The method for multiple hypothesis test
##'   correction. See `p.adjust` for all available methods.
##' @param robust.se (default FALSE) should White's Heteroscedasticity-robust
##'   Standard Errors be used rather then the traditional approach. (slower)
##' @return a list with elements controled by parameter resturn.pars. Options
##'   include:
##'   - X: P x N covariate matrix 
##'   - estimate: (P x D x nsample array) of linear model estimates
##'   - std.error: (P x D x nsample array) of standard error for the estimates
##'   - p.val: (P x D) matrix, unadjusted p-value for two-sided t-test
##'   - p.val.adj: (P x D) matrix, p-value for two-sided t-test adjusted
##'      according to `p.adj.method`
##'   - logWperp: (N x S) matrix, samples of the log scale from the scale model
##'   - logWpara: (D x N x S) array, samples of the log composition from
##'      the multinomial-Dirichlet
##'   - streaming: boolean, detnote if streaming was used. 
##'  Note, logWperp and logWpara are not returned if streaming is active
##' @export
##' @author Justin Silverman
aldex <- function(Y, X, data=NULL, nsample=2000,  GAMMA=NULL,
                  streamsize=8000,
                  return.pars=c("X", "estimate", "std.error", "p.val", "p.val.adj"),
                  return.samples=FALSE, p.adjust.method="BH", robust.se = FALSE) {
  N <- ncol(Y)
  D <- nrow(Y)

  ## calculate streaming threshold
  stream <- FALSE
  if (N*D*nsample*8/1000000 > streamsize) stream <- TRUE

  ## compute model matrix
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    X <- t(model.matrix(X, data))
    ## some error checking to make sure whats returned has right dimensions
  } 

  ## perform streaming 
  out <- list()
  nsample.remaining <- nsample
  iter <- 1
  if (stream) {
    nsample.local <- floor(streamsize*1000000/(N*D*8))
    if (nsample.local < 1) stop("streamsize too small")
  } else {
    nsample.local <- nsample
  }
  while (nsample.remaining > 0) {
    nsample.remaining <- nsample.remaining - nsample.local
    out[[iter]] <- aldex.lm.internal(Y, X, nsample.local, GAMMA, stream, robust.se)
    iter <- iter+1
  }
  ## combine output of the different streams
  out <- combine.streams(out)

  # p-value calculations, accounting for sign changes
  p.lower <- array(pmin(1, 2 * out$p.lower), dim = dim(out$p.lower))
  p.upper <- array(pmin(1, 2 * out$p.upper), dim = dim(out$p.upper))
  p.lower.adj <- apply(p.lower, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.upper.adj <- apply(p.upper, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.lower.mean <- apply(p.lower, c(2,1), mean)
  p.upper.mean <- apply(p.upper, c(2,1), mean)
  p.res <- c()
  for(col_i in 1:ncol(p.lower.mean)) {
    tmp_mat <- cbind(p.lower.mean[,col_i],
                     p.upper.mean[,col_i])
    p.res <- rbind(p.res, apply(tmp_mat, 1, min))
  }
  p.lower.mean.adj <- apply(p.lower.adj, c(1,2), mean)
  p.upper.mean.adj <- apply(p.upper.adj, c(1,2), mean)
  p.adj.res <- c()
  for(col_i in 1:ncol(p.lower.mean.adj)) {
    tmp_mat <- cbind(p.lower.mean.adj[,col_i],
                     p.upper.mean.adj[,col_i])
    p.adj.res <- rbind(p.adj.res, apply(tmp_mat, 1, min))
  }

  res <- list()
  res$streaming <- stream
  if ("estimate" %in% return.pars) res$estimate <- out$estimate
  if ("std.error" %in% return.pars) res$std.error <- out$std.error
  if ("p.val" %in% return.pars) res$p.val <- p.res
  if ("p.val.adj" %in% return.pars) res$p.val.adj <- p.adj.res
  if ("logWperp" %in% return.pars) res$logWperp <- out$logWperp
  if ("logWpara" %in% return.pars) res$logWpara <- out$logWpara
  return(res)
  ## TODO write a good "summary" function and wrap this all in S3 class
}


aldex.lm.internal <- function(Y, X, nsample, GAMMA=NULL, stream, robust.se=FALSE) {
  N <- ncol(Y)
  D <- nrow(Y)
  
  ## dirichlet sample
  logWpara <- array(NA, c(D, N, nsample))
  for (i in 1:N) {
    logWpara[,i,] <- log2(rDirichlet(nsample, Y[,i]+0.5))
  }

  ## sample from scale model
  if (is.null(GAMMA)){
    stop("You probably want a scale model :)")
  } else if (is.function(GAMMA)) {
    logWperp <- GAMMA(X, Y, logWpara)
    ## some error checking to make sure whats returned has right dimensions
  } 

  ## compute scaled abundances (W)
  logW <- sweep(logWpara, c(2,3), logWperp, FUN=`+`)

  ## fit linear model
  out <- fflm(aperm(logW, c(2,1,3)),t(X), robust.se) # TODO change fflm so it
                                          # gives correct output dimensions and
                                          # takes correct inputs dimensions --
                                          # without needing aperm
  if (!stream){
    out$logWpara <- logWpara
    out$logWperp <- logWperp
  }

  return(out)
 }
