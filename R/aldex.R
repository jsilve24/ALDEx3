##' ALDEx3 Linear Models
##'
##' 
##' @title ALDEx3 Linear Modles
##' @param Y an (D x N) matrix of sequence count, N is number of samples, D
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

  ## compute model matrix
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    formula <- X
    X <- t(model.matrix(X, data))
    ## some error checking to make sure whats returned has right dimensions
  } else {
    formula <- NULL
  }

  out <- list()
  iter <- 1
  chunk.sizes <- aldex.getchunksizes(Y, X, nsample, streamsize)
  for(chunk.size in chunk.sizes) {
    chunk <- aldex.chunk(Y, X, chunk.size, scale, scale.args,
                         return.pars)
    chunk <- c(chunk, fflm(aperm(chunk$logW, c(2,1,3)), t(X), test))
    out[[iter]] <- chunk
    iter <- iter + 1
  }
  out <- combine.streams(out)
  pval.l <- aldex.pvals(out$p.lower, out$p.upper,
                        p.adjust.method)

  res <- list()
  res$streaming <- length(chunk.sizes)>1
  res$X <- X
  if ("estimate" %in% return.pars) res$estimate <- out$estimate
  out$estimate <- NULL
  if ("std.error" %in% return.pars) res$std.error <- out$std.error
  out$std.error <- NULL
  if ("p.val" %in% return.pars) res$p.val <- pval.l$p.res
  if ("p.val.adj" %in% return.pars) res$p.val.adj <- pval.l$p.adj.res
  rm(pval.l)
  if ("logWperp" %in% return.pars) res$logWperp <- out$logWperp
  out$logWperp <- NULL
  if ("logWpara" %in% return.pars) res$logWpara <- out$logWpara
  out$logWpara <- NULL
  if (!is.null(data)) res$data <- data
  if (!is.null(formula)) res$formula <- formula
  return(res)
  ## TODO write a good "summary" function and wrap this all in S3 class
}

