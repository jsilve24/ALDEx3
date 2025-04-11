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
##' @param scale the scale model, can be a function or an N x nsample matrix.
##'   The API for writing your own scale models is documented below in examples.
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
##' @param test (default t.HC3), "t", t test is performed for each covariate
##'   (fast); "t.HC0" Heteroskedasticsity-Robust Standard Errors used (HC0;
##'   White's; slower); "t.HC3" (default) Heteroskedasticsity-Robust Standard
##'   Errors used (HC3; unlike HC0, this includes a leverage adjustment and is
##'   better for small sample sizes or when there are data with high leverage;
##'   slowest). To learn more about these, loko at Long and Ervin (2000) Using
##'   Heteroscedasticity Consistent Standard Errors in the Linear Regression
##'   Model, The American Statistician.
##' @param ... parameters to be passed to the scale model (if a function is
##'   provided)
##' @return a list with elements controled by parameter resturn.pars. Options
##'   include: - X: P x N covariate matrix - estimate: (P x D x nsample) array
##'   of linear model estimates - std.error: (P x D x nsample) array of standard
##'   error for the estimates - p.val: (P x D) matrix, unadjusted p-value for
##'   two-sided t-test - p.val.adj: (P x D) matrix, p-value for two-sided t-test
##'   adjusted according to `p.adj.method` - logWperp: (N x S) matrix, samples
##'   of the log scale from the scale model - logWpara: (D x N x S) array,
##'   samples of the log composition from the multinomial-Dirichlet - streaming:
##'   boolean, detnote if streaming was used. Note, logWperp and logWpara are
##'   not returned if streaming is active
##' @examples
##' \dontrun{
##' Y <- matrix(1:110, 10, 11)
##' condition <- c(rep(0, 5), rep(1, 6))
##' data <- data.frame(condition=condition)
##' ## demonstrate formula interface and passing optional argument (gamma) to
##' ## the scale model (clr)
##' res <- aldex(Y, ~condition, data, nsample=2000 scale=clr, gamma=0.5)
##' 
##' ## demonstrating how to write a custom scale model, I will write a model that
##' ## generalizes total sum scaling (where we assume no change between 
##' ## conditions)
##' ## Functions can include parameters X (model matrix), Y, and logWpara.
##' ## If included in the function definition, those parameters will be passed
##' ## dynamically when aldex is running. Other optional parameters (gamma)
##' ## can be passed as additional arguments to the aldex function
##' tss <- function(X, logWpara, gamma=0.5) {
##'   P <- nrow(X)
##'   nsample <- dim(logWpara)[3]
##'   Lambdaperp <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
##'   logWperp <- t(X)%*% Lambdaperp
##'   return(logWperp)
##' }
##' res <- aldex(Y, ~condition, data, nsample=2000, scale=tss, gamma=0.75)
##' }
##' @export
##' @author Justin Silverman
aldex <- function(Y, X, data=NULL, nsample=2000,  scale=NULL,
                  streamsize=8000,
                  return.pars=c("X", "estimate", "std.error", "p.val", "p.val.adj",
                                "logWpara", "logWperp"),
                  return.samples=FALSE, p.adjust.method="BH", test="t.HC3", ...) {
  scale.args <- list(...)

  N <- ncol(Y)
  D <- nrow(Y)

  ## calculate streaming threshold
  stream <- FALSE
  if (N*D*nsample*8/100000 > streamsize) stream <- TRUE

  ## compute model matrix
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    formula <- X
    X <- t(model.matrix(X, data))
    ## some error checking to make sure whats returned has right dimensions
  } else {
    formula <- NULL
  }

  ## perform streaming 
  out <- list()
  nsample.remaining <- nsample
  iter <- 1
  if (stream) {
    nsample.local <- floor(streamsize*100000/(N*D*8))
    if (nsample.local < 1) stop("streamsize too small")
  } else {
    nsample.local <- nsample
  }
  while (nsample.remaining > 0) {
    ## update sample sizes
    sample.size <- min(nsample.local, nsample.remaining)
    nsample.remaining <- nsample.remaining - sample.size

    ## stream
    out[[iter]] <- aldex.internal(Y, X, sample.size, scale, stream, test,
                                  scale.args, return.pars)
    iter <- iter+1
  }
  ## combine output of the different streams
  out <- combine.streams(out)


  #### p-value calculations, accounting for sign changes ####
  p.lower <- array(pmin(1, 2 * out$p.lower), dim = dim(out$p.lower))
  out$p.lower <- NULL # free up memory
  p.upper <- array(pmin(1, 2 * out$p.upper), dim = dim(out$p.upper))
  out$p.upper <- NULL # free up memory

  p.lower.adj <- apply(p.lower, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })
  p.upper.adj <- apply(p.upper, c(1,3), function(item) {
    p.adjust(item, method=p.adjust.method)
  })

  p.lower.mean <- apply(p.lower, c(2,1), mean)
  p.upper.mean <- apply(p.upper, c(2,1), mean)
  rm(p.lower, p.upper)


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
  #### END p value computation

  res <- list()
  res$streaming <- stream
  res$X <- X
  if ("estimate" %in% return.pars) res$estimate <- out$estimate
  out$estimate <- NULL
  if ("std.error" %in% return.pars) res$std.error <- out$std.error
  out$std.error <- NULL
  if ("p.val" %in% return.pars) res$p.val <- p.res
  rm(p.res)
  if ("p.val.adj" %in% return.pars) res$p.val.adj <- p.adj.res
  rm(p.adj.res)
  if ("logWperp" %in% return.pars) res$logWperp <- out$logWperp
  out$logWperp <- NULL
  if ("logWpara" %in% return.pars) res$logWpara <- out$logWpara
  out$logWpara <- NULL
  if (!is.null(data)) res$data <- data
  if (!is.null(formula)) res$formula <- formula
  return(res)
  ## TODO write a good "summary" function and wrap this all in S3 class
}


aldex.internal <- function(Y, X, nsample, scale=NULL, stream,
                           test, scale.args, return.pars) {
  N <- ncol(Y)
  D <- nrow(Y)

  
  ## dirichlet sample
  logWpara <- array(NA, c(D, N, nsample))
  for (i in 1:N) {
    logWpara[,i,] <- log2(rDirichlet(nsample, Y[,i]+0.5))
  }

  ## sample from scale model
  if (is.null(scale)){
    stop("You probably want a scale model :)")
  } else if (is.matrix(scale)) {
    stopifnot(dim(scale)==c(N, nsample))
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

  ## fit linear model
  out <- fflm(aperm(logW, c(2,1,3)),t(X), test) 

  if (!stream & ("logWperp" %in% return.pars)){
    out$logWperp <- logWperp
  }
  if (!stream & ("logWpara" %in% return.pars)){
    out$logWpara <- logWpara
  }
  return(out)
 }
