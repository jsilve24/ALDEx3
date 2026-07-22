##' ALDEx3 Linear Models
##'
##' Sequence counts determine the relative composition of each sample, but not
##' its total abundance or scale. The \code{scale} argument tells ALDEx3 what is
##' known about that missing quantity and how uncertain that knowledge is.
##'
##' Choose a scale model based on the information available:
##' \itemize{
##'   \item Use \code{\link{sample.sm}} for external measurements such as flow
##'     cytometry, qPCR, or spike-ins for individual samples.
##'   \item Use \code{\link{coefficient.sm}} when prior knowledge concerns a
##'     group difference, treatment effect, time trend, or another model term.
##'   \item Use \code{\link{clr.sm}} when CLR normalization is a reasonable
##'     center but its implied scale differences are uncertain.
##'   \item Use \code{\link{tss.sm}} when the starting assumption is no
##'     systematic scale difference between groups.
##'   \item Write a custom scale-model function when the built-in models do not
##'     match the information available for the study. The "Writing a custom
##'     scale model" example below describes the required function interface
##'     and \code{N x nsample} output.
##' }
##'
##' @title ALDEx3 Linear Models
##' @param Y an (D x N) matrix of sequence count, N is number of samples, D
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
##'    argument (see nlme documentation for how to use "correlation" argument);
##'    or "blmm": an ALDEx3-specific approximate mixed-effects engine that uses
##'    one batched profiled mixed-model anchor fit per feature, draw-specific
##'    local covariance updates, and exact conditional fixed-effect solves.
##'    The approximation is only in the variance-component optimisation step.
##'    See the mixed-effects vignette for details, validation guidance, and the
##'    runtime comparison with exact \code{lme4}.
##' @param nsample number of monte carlo replicates
##' @param scale A scale-model function, such as \code{clr.sm}, \code{tss.sm},
##'   \code{sample.sm}, or \code{coefficient.sm}; alternatively, an
##'   \code{N x nsample} matrix of log2-scale draws. Additional arguments for a
##'   scale-model function are passed through \code{...}. See "Writing a custom
##'   scale model" in the Examples section for the required interface.
##' @param streamsize (default 8000) memory footprint (approximate) at which to
##'   use streaming. This should be thought of as the number of Mb for each
##'   streaming chunk. If D*N*nsample*8/1000000 is less than streamsize then no
##'   streaming will be performed. Note, to conserve memory, samples from the
##'   Dirichlet and scale models will not be returned if streaming is used.
##'   Streaming can be turned off by setting streamsize=Inf.
##' @param n.cores (default detectCores()-1) If method is 'lme4' or 'blmm',
##'   use this many cores for feature-level parallelism across taxa/features.
##'   For \code{method = "blmm"}, parallelism is across features, not across
##'   Monte Carlo draws within a feature.
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
##' @param onesided (default: FALSE) if sided return p-values for two-sided
##'   test. Otherwise if "lower" or "upper" return one-sided test corresponding
##'   to test that estimate is negative or positive respectively.
##' @param ... Additional scale-model arguments, such as \code{gamma},
##'   \code{s.mu}, or \code{c.mu}. For \code{method = "nlme"}, this may also
##'   include \code{random} or \code{correlation}.
##' @return a list with elements controled by parameter resturn.pars. Options
##'   include: - X: P x N covariate matrix - estimate: (P x D x nsample) array
##'   of linear model estimates - std.error: (P x D x nsample) array of standard
##'   error for the estimates - p.val: (P x D) matrix, unadjusted p-value for
##'   two-sided t-test - p.val.adj: (P x D) matrix, p-value for two-sided t-test
##'   adjusted according to `p.adj.method` - logScale: (N x S) matrix, samples
##'   of the log scale from the scale model - logComp: (D x N x S) array,
##'   samples of the log composition from the multinomial-Dirichlet - streaming:
##'   boolean, detnote if streaming was used. - random.effects (Pr x N x S): if
##'   using mixed effects models, return all Pr random-effects covariance terms
##'   and the residual variance. For \code{method = "blmm"}, correlated
##'   random-effect covariance terms are returned when present, matching the
##'   \code{lme4} naming convention as closely as possible. Note, logScale and
##'   logComp are not returned if streaming is active.
##' @examples
##' \donttest{
##' Y <- matrix(1:110, 10, 11)
##' condition <- c(rep(0, 5), rep(1, 6))
##' data <- data.frame(condition=condition)
##' ## Use CLR as the center and allow coefficient-scale uncertainty of SD 0.5.
##' res <- aldex(Y, ~condition, data, nsample=2000, scale=clr.sm, gamma=0.5)
##' 
##' ## Writing a custom scale model
##' ## A custom model must return one log2-scale draw per sample and Monte
##' ## Carlo replicate. Here, external_scale supplies a known value per sample.
##' known_scale <- function(X, logComp, external_scale) {
##'   nsample <- dim(logComp)[3]
##'   matrix(rep(external_scale, nsample), nrow=ncol(X), ncol=nsample)
##' }
##' res <- aldex(Y, ~condition, data, nsample=2000, scale=known_scale,
##'              external_scale=rep(0, ncol(Y)))
##' }
##' @importFrom reformulas nobars
##' @importFrom parallel detectCores
##' @importFrom methods formalArgs
##' @import stats
##' @export
##' @author Justin Silverman, Kyle McGovern
aldex <- function(Y, X, data=NULL, method="lm", nsample=2000,  scale=NULL,
                  streamsize=8000, n.cores=detectCores()-1,
                  return.pars=c("X", "estimate", "std.error", "p.val",
                                "p.val.adj", "logComp", "logScale"),
                  p.adjust.method="BH", test="t.HC3",
                  onesided=FALSE,
                  ...) {
  scale.args <- list(...)
  scale.args <- .upgrade_deprecated_scale_args(scale, scale.args)

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
  if ((method=="lme4")|(method=="nlme")|(method=="blmm")) {
    if(is.null(data)) stop("data should not be null if method is a mixed-effects engine")
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
    } else if (method=="blmm") {
      res <- blmm(aperm(out[[iter]]$logW, c(2,1,3)), formula,
                  data, n.cores)
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
  pval.l <- aldex.pvals(out$p.lower, out$p.upper, p.adjust.method,onesided)

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
