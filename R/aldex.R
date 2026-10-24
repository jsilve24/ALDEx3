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
##' @param Y an (D x N) matrix of sequence count, N is number of samples, D is
##'   number of taxa or genes
##' @param X either a formula (in which case DATA must be non-null) or a model
##'   matrix of dimension P x N (P is number of linear model covariates). If a
##'   formula is passed it should not include the target Y, e.g., should simply
##'   be "~condition-1" (note the lack of the left hand side). If using lme4,
##'   this should be the formula including random effects.
##' @param data a data frame for use with formula, must have N rows
##' @param method (default lm) The regression method; "lm": linear regression,
##'   "lme4": linear mixed effects regression with lme4; "nlme": linear mixed
##'   effects models with nlme, REQUIRES "random" argument representing random
##'   effects be passed into `aldex` function, can also pass "correlation" as
##'   argument (see nlme documentation for how to use "correlation" argument);
##'   or "blmm": an ALDEx3-specific approximate mixed-effects engine. For each
##'   feature, BLMM estimates covariance parameters jointly across Monte Carlo
##'   draws, then makes a local adjustment for each draw. This shared fit is
##'   called the anchor. Fixed effects and standard errors are calculated
##'   using each draw's adjusted covariance parameters; p-values use the
##'   approximation selected by \code{test}. If the shared fit fails or its
##'   curvature cannot be used, BLMM fits every draw for that feature separately
##'   with \code{lmerTest} (which uses \code{lme4} for fitting). If only a
##'   draw's local update or inference fails validation, only that draw is
##'   refitted. A warning reports these refits. BLMM reuses matrix calculations
##'   and automatically chooses sparse or dense computation; no extra argument
##'   is needed. See the mixed-effects vignette for details and timing examples.
##' @param nsample number of monte carlo replicates
##' @param scale A scale-model function, such as \code{clr.sm}, \code{tss.sm},
##'   \code{sample.sm}, or \code{coefficient.sm}; alternatively, an \code{N x
##'   nsample} matrix of log2-scale draws. Additional arguments for a
##'   scale-model function are passed through \code{...}. See
##'   "Writing a custom scale model" in the Examples section for the required
##'   interface.
##' @param streamsize (default 8000) working-size threshold for streaming. The
##'   existing heuristic uses ten times the raw response-array size in decimal
##'   MB: streaming occurs when \code{D*N*nsample*8/100000 > streamsize}. Thus
##'   the default corresponds to about 800 MB of raw response draws, not an 8 GB
##'   array. Each chunk contains at most \code{floor(streamsize*100000/(D*N*8))}
##'   draws. This is a size heuristic, not a hard process-memory limit. To
##'   conserve memory, Dirichlet and scale samples are not returned when
##'   streaming is used. Set \code{streamsize=Inf} to disable streaming. BLMM
##'   fits a separate anchor per chunk, so changing this threshold can change
##'   its approximate results.
##' @param n.cores (default detectCores()-1) Number of workers for mixed models.
##'   The fully fitted engines distribute feature-draw fits across workers. For
##'   \code{method = "blmm"}, parallelism is across features, not across Monte
##'   Carlo draws within a feature.
##' @param return.pars what results should be returned, see return section
##'   below.
##' @param p.adjust.method (default BH) The method for multiple hypothesis test
##'   correction. See `p.adjust` for all available methods.
##' @param test inference selector. Defaults depend on `method`: `"t.HC3"` for
##'   `"lm"`, `"t.satterthwaite"` for `"blmm"`, and native `"t"` inference for
##'   `"lme4"` and `"nlme"`. For `lm`, choose `"t"`, `"t.HC0"`, or `"t.HC3"`;
##'   for `blmm`, choose `"t.satterthwaite"` or `"t.residual"`; for `lme4` and
##'   `nlme`, only `"t"` is supported.
##' @details For \code{method = "lm"}, \code{"t"} uses ordinary residual
##'   variance, \code{"t.HC0"} uses the White sandwich covariance, and
##'   \code{"t.HC3"} applies its leverage correction. For \code{"blmm"},
##'   \code{"t.satterthwaite"} uses the shared-anchor approximation and
##'   \code{"t.residual"} uses residual degrees of freedom. Draws that fail
##'   BLMM's checks are refitted with \code{lmerTest} and use its Satterthwaite
##'   inference. This also applies when \code{test = "t.residual"}: draws that
##'   pass the checks use residual df, but refitted draws use Satterthwaite df.
##'   The \code{lme4}
##'   engine uses per-draw lmerTest Satterthwaite inference internally;
##'   \code{nlme} uses its native denominator df. Their public selector is
##'   \code{"t"}; neither accepts the two BLMM selectors. Names are
##'   case-sensitive.
##'
##'   Each contrast row tests \eqn{H_0: c^T\beta=0}. Its SE uses the full
##'   fitted coefficient covariance, including off-diagonal entries. This
##'   argument specifies linear hypotheses, not R's factor-coding contrasts.
##'   For \code{nlme}, a contrast may involve only coefficients with the same
##'   native denominator df; unsupported combinations produce an error.
##'
##'   Each coefficient and contrast forms a separate adjustment family across
##'   features. For two-sided tests, directional probabilities are doubled and
##'   capped at one, then adjusted within each draw and direction. ALDEx3
##'   averages over draws and selects the smaller directional aggregate.
##'   Combining discovery lists across coefficients or contrasts carries no
##'   joint FDR guarantee. The supplied scale model and original design are
##'   unchanged. Streaming may change BLMM results because each chunk has its
##'   own shared anchor; contrasts use the same fits within each chunk.
##' @param contrasts optional numeric K x P matrix with one fixed-effect
##'   contrast per row. Its uniquely named columns must be the complete
##'   coefficient set and are reordered to the fitted design. Contrast rows
##'   have unique nonempty names, are finite and nonzero, and need not sum to
##'   zero or be independent. Contrasts require a full-rank fixed-effect design
##'   with positive residual degrees of freedom. They do not change composition
##'   or scale sampling or adjustment families. Coefficient results are
##'   unchanged unless checks on an added contrast cause a draw to be refitted.
##' @param onesided (default: FALSE) if sided return p-values for two-sided
##'   test. Otherwise if "lower" or "upper" return one-sided test corresponding
##'   to test that estimate is negative or positive respectively.
##' @param ... Additional scale-model arguments, such as \code{gamma},
##'   \code{s.mu}, or \code{c.mu}. For \code{method = "nlme"}, this may also
##'   include \code{random} or \code{correlation}.
##' @return An \code{aldex} list. \code{X} (P x N), \code{streaming},
##'   \code{method}, and the resolved \code{test} are always recorded. Here P is
##'   the coefficient count, D the feature count, N the observation count, and
##'   S the Monte Carlo draw count. \code{return.pars} controls:
##'   \describe{
##'     \item{estimate, std.error, df}{P x D x S arrays of fitted coefficients,
##'       within-fit standard errors, and degrees of freedom.}
##'     \item{p.val, p.val.adj}{P x D matrices of final raw and adjusted
##'       probabilities, respecting \code{onesided} and \code{p.adjust.method}.}
##'     \item{logComp, logScale}{D x N x S log2-composition draws and N x S
##'       log2-scale draws. These are not retained during streaming.}
##'     \item{random.effects}{Pr x D x S arrays of Pr random-effect variance
##'       and covariance terms, including residual variance, for mixed models.}
##'   }
##'   When contrasts are requested, \code{contrasts$matrix} always records their
##'   reordered matrix. \code{contrasts$estimate}, \code{contrasts$std.error},
##'   and \code{contrasts$df} have dimensions K x D x S; \code{contrasts$p.val}
##'   and \code{contrasts$p.val.adj} have dimensions K x D. K is the number of
##'   contrast rows. The same \code{return.pars} fields control coefficient and
##'   contrast results together. Names and singleton dimensions are retained.
##'   Use \code{summary(fit, component = "contrasts")} for stored contrast
##'   summaries; the default summary continues to show coefficients.
##' @examples
##' set.seed(15)
##' Y <- matrix(rpois(3 * 12, 50), 3, 12)
##' metadata <- data.frame(group = rep(c(0, 1), each = 6))
##' X <- model.matrix(~ group, metadata)
##' C <- matrix(c(0, 1), 1, dimnames = list("group difference", colnames(X)))
##' fit <- aldex(Y, ~ group, metadata, nsample = 8, scale = clr.sm,
##'              gamma = 0.2, contrasts = C)
##' summary(fit, component = "contrasts")
##'
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
                  return.pars=c("X", "estimate", "std.error", "df", "p.val",
                                "p.val.adj", "logComp", "logScale"),
                  p.adjust.method="BH", test=NULL,
                  onesided=FALSE, contrasts=NULL,
                  ...) {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% c("lm", "blmm", "lme4", "nlme"))) {
    stop("method must be one of 'lm', 'blmm', 'lme4', or 'nlme'")
  }
  if (is.null(test)) {
    test <- switch(method, lm="t.HC3", blmm="t.satterthwaite",
                   lme4="t", nlme="t")
  }
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("test must be a single character inference selector")
  }
  allowed.tests <- switch(method, lm=c("t", "t.HC0", "t.HC3"),
                          blmm=c("t.satterthwaite", "t.residual"),
                          lme4="t", nlme="t")
  if (!(test %in% allowed.tests)) {
    stop("test '", test, "' is not supported for method '", method,
         "'; allowed selectors are: ", paste(allowed.tests, collapse=", "))
  }
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

  if (!is.matrix(X) || !is.numeric(X) || ncol(X) != N || !all(is.finite(X))) {
    stop("X must be a finite numeric P x N model matrix")
  }
  if (is.null(rownames(X)) || anyNA(rownames(X)) || any(rownames(X) == "")) {
    rownames(X) <- paste0("parameter_", seq_len(nrow(X)))
  }
  if (anyDuplicated(rownames(X))) {
    stop("X must have unique, nonempty coefficient names")
  }
  contrasts <- .aldex_validate_contrasts(contrasts, t(X))
  if (!is.null(contrasts) && qr(t(X))$rank < nrow(X)) {
    stop("contrasts require a full-rank fixed-effect design; rank-deficient ",
         "designs are unsupported")
  }
  if (!is.null(contrasts) && N <= nrow(X)) {
    stop("contrasts require positive residual degrees of freedom (N > P)")
  }
  if (method == "lm" && qr(t(X))$rank < nrow(X)) {
    stop("rank-deficient designs are not supported by method 'lm'")
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
      res <- fflm(aperm(out[[iter]]$logW, c(2,1,3)), t(X), test, contrasts)
    } else if ((method=="lme4")|method=="nlme") {
      res <- sr.mem(aperm(out[[iter]]$logW, c(2,1,3)), formula,
                    data, n.cores, method, mem.args, contrasts=contrasts)
    } else if (method=="blmm") {
      res <- blmm(aperm(out[[iter]]$logW, c(2,1,3)), formula,
                  data, n.cores,
                  inference=if (test == "t.satterthwaite") "anchor" else "residual",
                  contrasts=contrasts)
    }
    out[[iter]]$logW <- NULL # don't duplicate info in logComp and logScale 
    ## while ugly, the following loop should avoid shallow copy of logW and
    ## logComp currently in out[[iter]], that could be memory intensive. 
    for (name in names(res)) { out[[iter]][[name]] <- res[[name]] }
    iter <- iter+1
  }
  ## Combine streams without simplifying singleton contrast dimensions.
  out <- combine.streams(out)

  ## p-value calculations, accounting for sign changes
  pval.l <- aldex.pvals(out$p.lower, out$p.upper, p.adjust.method,onesided)
  contrast.pval <- if (!is.null(contrasts)) {
    aldex.pvals(out$contrast.p.lower, out$contrast.p.upper,
                p.adjust.method, onesided)
  } else NULL

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
  res$method <- method
  res$test <- test
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
  if ("df" %in% return.pars) {
    res$df <- out$df
    rownames(res$df) <- parameter.names
    colnames(res$df) <- entity.names
  }
  out$df <- NULL
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
  if (!is.null(contrasts)) {
    cres <- list(matrix=contrasts)
    if ("estimate" %in% return.pars) {
      cres$estimate <- out$contrast.estimate
      rownames(cres$estimate) <- rownames(contrasts)
      colnames(cres$estimate) <- entity.names
    }
    if ("std.error" %in% return.pars) {
      cres$std.error <- out$contrast.std.error
      rownames(cres$std.error) <- rownames(contrasts)
      colnames(cres$std.error) <- entity.names
    }
    if ("df" %in% return.pars) {
      cres$df <- out$contrast.df
      rownames(cres$df) <- rownames(contrasts)
      colnames(cres$df) <- entity.names
    }
    if ("p.val" %in% return.pars) {
      cres$p.val <- contrast.pval$p.res
      rownames(cres$p.val) <- rownames(contrasts)
      colnames(cres$p.val) <- entity.names
    }
    if ("p.val.adj" %in% return.pars) {
      cres$p.val.adj <- contrast.pval$p.adj.res
      rownames(cres$p.val.adj) <- rownames(contrasts)
      colnames(cres$p.val.adj) <- entity.names
    }
    res$contrasts <- cres
  }
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
