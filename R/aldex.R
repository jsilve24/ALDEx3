##' ALDEx3 Linear Models
##'
##' 
##' @title ALDEx3 Linear Modles
##' @param counts an (D x N) matrix of sequence count, N is number of samples, D
##'   is number of taxa or genes
##' @param X either a formula (in which case DATA must be non-null) or a model
##'   matrix of dimension P x N (P is number of linear model covariates). If a
##'   formula is passed it should not include the target Y, e.g., should simply
##'   be "~condition-1" (note the lack of the left hand side). If using lme4,
##'   this should be the formula including random effects.
##' @param data a data frame for use with formula, must have N rows
##' @param method (default lm) The regression method; "lm": linear regression,
##'   "lme4": linear mixed effects regression with lme4; "nlme": linear mixed
##'    effects models with nlme, REQUIRES "random" argument representing random
##'    effects be passed into `aldex` function, can also pass "correlation" as
##'    argument (see nlme documentation for how to use "correlation" argument).
##' @param nsample number of monte carlo replicates
##' @param scale the scale model, can be a function or an N x nsample matrix
##'   (samples should be given on log2-scale; e.g., samples should be of log of
##'   system scale). The API for writing your own scale models is documented
##'   below in examples.
##' @param streamsize (default 8000) memory footprint (approximate) at which to
##'   use streaming. This should be thought of as the number of Mb for each
##'   streaming chunk. If D*N*nsample*8/1000000 is less than streamsize then no
##'   streaming will be performed. Note, to conserve memory, samples from the
##'   Dirichlet and scale models will not be returned if streaming is used.
##'   Streaming can be turned off by setting streamsize=Inf.
##' @param n.cores (default detectCores()-1) If method is 'lme4', use this many
##'   cores for running mixed effects models in parallel.
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
##'   provided), may also be random or correlation arguments for nlme.
##' @return a list with elements controled by parameter resturn.pars. Options
##'   include: - X: P x N covariate matrix - estimate: (P x D x nsample) array
##'   of linear model estimates - std.error: (P x D x nsample) array of standard
##'   error for the estimates - p.val: (P x D) matrix, unadjusted p-value for
##'   two-sided t-test - p.val.adj: (P x D) matrix, p-value for two-sided t-test
##'   adjusted according to `p.adj.method` - logScale: (N x S) matrix, samples
##'   of the log scale from the scale model - logComp: (D x N x S) array,
##'   samples of the log composition from the multinomial-Dirichlet - streaming:
##'   boolean, detnote if streaming was used. - random.effects (Pr x N x S): if
##'   using mixed effects models, return all Pr random effects. Note, logScale
##'   and logComp are not returned if streaming is active.
##' @examples
##' \dontrun{
##' Y <- matrix(1:110, 10, 11)
##' condition <- c(rep(0, 5), rep(1, 6))
##' data <- data.frame(condition=condition)
##' ## demonstrate formula interface and passing optional argument (gamma) to
##' ## the scale model (clr)
##' res <- aldex(Y, ~condition, data, nsample=2000, scale=clr.sm, gamma=0.5)
##' 
##' ## demonstrating how to write a custom scale model, I will write a model
##' ## that generalizes total sum scaling (where we assume no change between 
##' ## conditions)
##' ## Functions can include parameters X (model matrix), Y, and logComp.
##' ## If included in the function definition, those parameters will be passed
##' ## dynamically when aldex is running. Other optional parameters (gamma)
##' ## can be passed as additional arguments to the aldex function
##' tss <- function(X, logComp, gamma=0.5) {
##'   P <- nrow(X)
##'   nsample <- dim(logComp)[3]
##'   LambdaScale <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
##'   logScale <- t(X)%*% LambdaScale
##'   return(logScale)
##' }
##' res <- aldex(Y, ~condition, data, nsample=2000, scale=tss, gamma=0.75)
##' }
##' @importFrom lme4 nobars
##' @importFrom parallel detectCores
##' @export
##' @author Justin Silverman, Kyle McGovern
aldex <- function(Y, X, data=NULL, method="lm", nsample=2000,  scale=NULL,
                  streamsize=8000, n.cores=detectCores()-1,
                  return.pars=c("X", "estimate", "std.error", "p.val",
                                "p.val.adj", "logComp", "logScale"),
                  return.samples=FALSE, p.adjust.method="BH", test="t.HC3",
                  ...) {
  scale.args <- list(...)

  mem.args <- list()
  # Elements to move
  to_remove <- c("random", "correlation")
  for(to_move in to_remove) {
    # Check and move
    if (to_move %in% names(scale.args)) {
      mem.args[[to_move]] <- scale.args[[to_move]]
      scale.args[[to_move]] <- NULL
    }
  }

  N <- ncol(Y)
  D <- nrow(Y)

  ## Checks if scale is matrix
  if(is.matrix(scale)) {
    if(any(dim(scale)!=c(N, nsample))) {
      stop("When scale is a matrix it must be dim c(N, nsample)")
    }
    if (!all(is.finite(scale))) {
      stop("Some elements of scale matrix are infinite or NA")
    }
  }

  ## Checks for mixed effects modeling
  if ((method=="lme4")|(method=="nlme")) {
    if(is.null(data)) stop("data should not be null if method=\"lme4\"")
    if(!inherits(X, "formula")) stop("X should be a mixed effects ",
                                     "formula if method=\"lme4\" ",
                                     "or a fixed effects formula ",
                                     "if method=\"nlme\"")
  }

  ## compute model matrix
  if (inherits(X, "formula")) {
    if (is.null(data)) stop("data should not be null if X is a formula")
    formula <- X
    ## nobars removes random effects
    X <- t(model.matrix(nobars(X), data))
    ## some error checking to make sure whats returned has right dimensions
  } else {
    formula <- NULL
  }

  ## perform streaming 
  out <- list()
  nsample.remaining <- nsample
  iter <- 1
  stream <- ifelse(N*D*nsample*8/100000 > streamsize, TRUE, FALSE)
  if (stream) {
    nsample.local <- floor(streamsize*100000/(N*D*8))
    if (nsample.local < 1) stop("streamsize too small")
    ## override if logComp and logScale were requested 
    return.pars <- return.pars[!(return.pars %in% c("logComp", "logScale"))]
  } else {
    nsample.local <- nsample
  }
  samples.processed <- 0
  while (nsample.remaining > 0) {
    ## update sample sizes
    sample.size <- min(nsample.local, nsample.remaining)
    nsample.remaining <- nsample.remaining - sample.size
    if(is.matrix(scale)) {
      ## Subset scale matrix
      matrix.start <- nsample-nsample.remaining-sample.size+1
      matrix.end <- nsample-nsample.remaining
      out[[iter]] <- aldex.sampler(Y, X, sample.size,
                                   scale[,matrix.start:matrix.end,drop=FALSE],
                                   scale.args, return.pars)
    } else {
      out[[iter]] <- aldex.sampler(Y, X, sample.size, scale, scale.args,
                                   return.pars)
    }
    ## fit linear model
    if (method=="lm") {
      res <- fflm(aperm(out[[iter]]$logW, c(2,1,3)), t(X), test)
    } else if ((method=="lme4")|method=="nlme") {
      res <- sr.mem(aperm(out[[iter]]$logW, c(2,1,3)), formula,
                    data, n.cores, method, mem.args)
    }
    out[[iter]]$logW <- NULL # don't duplicate info in logComp and logScale 
    ## while ugly, the following loop should avoid shallow copy of logW and
    ## logComp currently in out[[iter]], that could be memory intensive. 
    for (name in names(res)) { out[[iter]][[name]] <- res[[name]] }
    iter <- iter+1
  }
  ## combine output of the different streams
  out <- combine.streams(out)

  ## p-value calculations, accounting for sign changes
  pval.l <- aldex.pvals(out$p.lower, out$p.upper, p.adjust.method)

  ## collect names to make output nicer
  if (is.null(colnames(Y))) {
    sample.names <- paste0("sample_", 1:N)
  } else {
    sample.names <- colnames(Y)
  }
  if (is.null(rownames(Y))) {
    entity.names <- paste0("entity_", 1:D)
  } else {
    entity.names <- rownames(Y)
  }
  if (is.null(rownames(X))) {
    parameter.names <- paste0("parameter_", 1:nrow(X))
  } else {
    parameter.names <- rownames(X)
  }

  res <- list()
  res$streaming <- stream
  res$X <- X
  if ("estimate" %in% return.pars) {
    res$estimate <- out$estimate
    rownames(res$estimate) <- parameter.names
    colnames(res$estimate) <- entity.names
  }
  out$estimate <- NULL
  if ("std.error" %in% return.pars) {
    res$std.error <- out$std.error
    rownames(res$std.error) <- parameter.names
    colnames(res$std.error) <- entity.names
  }
  out$std.error <- NULL
  if ("p.val" %in% return.pars) {
    res$p.val <- pval.l$p.res
    rownames(res$p.val) <- parameter.names
    colnames(res$p.val) <- entity.names
  }
  pval.l$p.res <- NULL
  if ("p.val.adj" %in% return.pars) {
    res$p.val.adj <- pval.l$p.adj.res
    rownames(res$p.val.adj) <- parameter.names
    colnames(res$p.val.adj) <- entity.names
  }
  pval.l$p.adj.res <- NULL
  if ("logScale" %in% return.pars) {
    res$logScale <- out$logScale
    rownames(res$logScale) <- sample.names
  }
  out$logScale <- NULL
  if ("logComp" %in% return.pars) {
    res$logComp <- out$logComp
    rownames(res$logComp) <- entity.names
    colnames(res$logComp) <- sample.names
  }
  out$logComp <- NULL
  if ("random.effects" %in% return.pars) {
    res$random.effects <- out$random.eff
    colnames(res$random.effects) <- entity.names
  }
  out$random.eff <- NULL

  if (!is.null(data)) res$data <- data
  if (!is.null(formula)) res$formula <- formula
  attr(res, "class") <- "aldex"
  return(res)
}


