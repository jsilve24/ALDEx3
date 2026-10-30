##' @title ALDEx3 Linear Models
##' @description Fit linear or linear mixed models to sequence count data
##' while accounting for uncertainty in both composition and scale (total
##' abundance).
##'
##' For t-based inference, for each of \code{nsample} Monte Carlo draws, \code{aldex} samples each
##' sample's composition from a Dirichlet-multinomial posterior, samples its
##' log2 scale from the scale model, and adds the two to obtain log2
##' abundances. It then fits the chosen model to every feature, adjusts the
##' p-values across features within the draw, and averages estimates, standard
##' errors and adjusted p-values over draws.
##'
##' @param Y a D x N matrix of sequence counts, where D is the number of taxa or
##'   genes and N the number of samples.
##' @param X a one-sided model formula (no response), used with \code{data},
##'   or for \code{method = "lm"} a P x N numeric model matrix (P is the number
##'   of coefficients). For \code{"blmm"} and \code{"lme4"}, the formula
##'   includes random effects in \code{lme4} syntax, e.g.
##'   \code{~ treat + (1 | subject)}. For \code{"nlme"}, give the fixed-effects
##'   formula here and the random effects through \code{random} in \code{...}.
##' @param data a data frame with N rows, required when \code{X} is a formula.
##' @param method (default \code{"lm"}) the regression engine:
##'   \describe{
##'     \item{\code{"lm"}}{fast ordinary linear regression.}
##'     \item{\code{"blmm"}}{ALDEx3's approximate linear mixed model, with
##'       \code{lmerTest} refits for draws that fail its checks (see Details).}
##'     \item{\code{"lme4"}}{a full \code{lmerTest} (\code{lme4}) mixed model fit
##'       to every draw.}
##'     \item{\code{"nlme"}}{an \code{nlme} mixed model; requires \code{random}
##'       and optionally accepts \code{correlation} (see \code{\link[nlme]{lme}}).}
##'   }
##' @param nsample number of Monte Carlo draws (at least two for bootstrap
##'   Scale Simulation Random Variables, or SSRVs).
##' @param bootstrap.strata Optional length-N nonmissing vector for independent
##'   subject resampling within strata in test="ssrv.bootstrap". By default,
##'   all designs use global paired row resampling. Strata affect the bootstrap
##'   branch only and do not provide cluster inference.
##' @param bootstrap.max.attempts Maximum attempts per full-rank bootstrap
##'   design (default 100); exhaustion stops with rejection diagnostics.
##' @param scale A scale-model function, such as \code{clr.sm}, \code{tss.sm},
##'   \code{sample.sm}, or \code{coefficient.sm}; alternatively, an \code{N x
##'   nsample} matrix of log2-scale draws. Additional arguments for a
##'   scale-model function are passed through \code{...}. With
##'   \code{test="ssrv.bootstrap"}, also accepts named coefficient specifications from
##'   \code{\link{compose.sm}} or a custom coefficient-scale callback. See
##'   "Writing a custom scale model" in the Examples section for the required
##'   interface.
##' @param streamsize (default 8000) size threshold for processing draws in
##'   chunks. Streaming occurs when \code{D*N*nsample*8/100000 > streamsize}, so
##'   the default corresponds to about 800 MB of response draws. Each chunk
##'   contains at most \code{floor(streamsize*100000/(D*N*8))} draws. This is a
##'   heuristic, not a hard memory limit. Composition and scale draws are not
##'   returned when streaming. Set \code{streamsize=Inf} to disable streaming.
##'   BLMM fits a separate anchor per chunk, so changing this threshold can
##'   change its results slightly.
##' @param n.cores (default detectCores()-1) Number of workers for mixed models.
##'   \code{"lme4"} and \code{"nlme"} distribute feature-draw fits across
##'   workers; \code{"blmm"} distributes features.
##' @param return.pars what results should be returned, see return section
##'   below.
##' @param p.adjust.method (default BH) The method for multiple hypothesis test
##'   correction. See `p.adjust` for all available methods.
##' @param test inference method. Defaults depend on `method`: `"t.HC3"` for
##'   `"lm"`, `"t.satterthwaite.adaptive"` for `"blmm"`, and `"t"` for
##'   `"lme4"` and `"nlme"`. For `lm`, choose `"t"`, `"t.HC0"`, `"t.HC3"`, or `"ssrv.bootstrap"`;
##'   for `blmm`, choose `"t.satterthwaite.adaptive"`, `"t.satterthwaite"`,
##'   or the deprecated `"t.residual"`; for `lme4` and `nlme`, only `"t"` is
##'   supported. Names are case-sensitive.
##' @details
##' \strong{Bootstrap SSRVs (1.8.0).} Select \code{test="ssrv.bootstrap"}
##'   with \code{method="lm"} for independent-subject fixed-effects inference.
##'   Each draw samples original-subject compositions. The bootstrap branch
##'   resamples subjects and design rows together; sparse Laplace uses the
##'   original subjects. Each fits ordinary least squares for all features and adds a single
##'   coefficient-scale shift across taxa. Hybrid selects a complete final-effect
##'   sequence per sparse coefficient by median feature 95% interval width,
##'   choosing Laplace on exact ties. When at least one sparse specification is
##'   present and all use Laplace, no subject or feature bootstrap is drawn;
##'   Gaussian coefficients also use the original subjects. Gaussian-only, callback, and existing
##'   scale models use subject bootstrap. No residual t errors or degrees of
##'   freedom are computed. Use \code{\link{compose.sm}} for sparse/Gaussian mixtures;
##'   every coefficient including the intercept must be explicitly assigned.
##'   Composable and custom coefficient models use original log2 proportions,
##'   without automatic CLR centering. Gaussian shifts describe total-scale
##'   coefficients with the same meaning as \code{coefficient.sm}.
##'   Existing sample-scale functions and matrices are projected through the
##'   same resampled fitting operator, preserving sample identity.
##'   Default \code{return.pars=c("posterior","estimate")} retains final
##'   absolute coefficient draws plus means and pointwise equal-tailed 95%
##'   intervals and raw and adjusted sign probabilities. \code{return.pars="posterior"} drops draws after computing
##'   exact quantiles. Summaries are stored in \code{posterior} (mean, lower,
##'   upper, level, interval) and \code{contrasts$posterior}. Probability
##'   matrices are stored as \code{p.val} and \code{p.val.adj} for each.
##'   They can also be requested individually via \code{return.pars}. The draw count
##'   and rank-rejection/scale-selection diagnostics are recorded. Sample arrays
##'   are processed one draw at a time; coefficient working memory remains
##'   proportional to P*D*nsample. RNG order is independent of streaming and
##'   retention. Inference-dependent plotting/effect/p-adjust helpers are
##'   unavailable. See \code{vignette("ALDEx3-bootstrap-SSRV")} for assumptions,
##'   custom joint coefficient draws, and examples.
##' \strong{Scale models.} Sequence counts determine the relative composition of each sample, but not
##' its total abundance or scale. The \code{scale} argument tells ALDEx3 what is
##' known about that missing quantity and how uncertain that knowledge is.
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
##'   \strong{Linear models.} For \code{method = "lm"}, \code{"t"} uses the
##'   ordinary residual variance, \code{"t.HC0"} the White sandwich covariance,
##'   and \code{"t.HC3"} its leverage-corrected version.
##'
##'   \strong{BLMM.} For each feature, BLMM fits one set of covariance
##'   parameters shared by all draws (the anchor). With the default
##'   \code{"t.satterthwaite.adaptive"}, each draw then refines its own
##'   covariance fit from the anchor with at most eight steps, and its
##'   estimates, standard errors and Satterthwaite degrees of freedom are
##'   computed at that fit. A draw is accepted only if it passes numerical
##'   checks of convergence, curvature and stability of its inference; draws
##'   that fail are refit with \code{lmerTest}, and a warning reports how many.
##'   Fits with all covariance parameters at zero, or with whole zero columns of
##'   the covariance Cholesky factor, get local inference without a refit. The
##'   checks are numerical; they do not guarantee calibrated p-values.
##'   \code{"t.satterthwaite"} is the cheaper shared-anchor approximation used
##'   by default in ALDEx3 1.5.0-1.6.2. The deprecated \code{"t.residual"} uses
##'   residual degrees of freedom for accepted draws and warns once per call.
##'   If \code{lmerTest} cannot produce valid inference for a refit draw's
##'   contrast, that draw's p-values are \code{NA} with a warning; errors during
##'   fitting stop the analysis. BLMM chooses dense or sparse computation
##'   automatically. See \code{vignette("ALDEx3-mixed-effects")} for details.
##'
##'   \strong{Other mixed models.} \code{"lme4"} uses \code{lmerTest}
##'   Satterthwaite inference for every draw; \code{"nlme"} uses its native
##'   denominator degrees of freedom.
##'
##'   \strong{Contrasts.} Each contrast row tests \eqn{H_0: c^T\beta=0}, using
##'   the full fitted coefficient covariance for t-based inference. Bootstrap
##'   SSRV scores count the signs of each within-draw joint sum of coefficients,
##'   rather than combining coefficient scores. These are linear hypotheses, not
##'   R's factor-coding contrasts. For \code{nlme}, a contrast may involve only
##'   coefficients with the same denominator degrees of freedom.
##'
##'   \strong{Multiple testing for t-based inference.} Each coefficient and contrast is a separate
##'   adjustment family across features. Within each draw, p-values are
##'   adjusted with \code{p.adjust.method} and then averaged over draws; for
##'   two-sided tests this is done for each direction and the smaller result
##'   is reported. Use \code{\link{aldex.p.adjust}} on the fitted object to
##'   adjust several coefficients or contrasts as one family without refitting.
##'
##'   \strong{Multiple testing for bootstrap SSRVs.} Each coefficient and contrast is adjusted separately across all
##'   D supplied features after computing the raw probability. For final draws
##'   theta, the raw two-sided score is min(1, 2 * min(mean(theta <= 0),
##'   mean(theta >= 0))). Zeros count in both tails; all-zero draws give 1.
##'   A nonfinite draw makes the score unavailable. Scores use S draws without
##'   smoothing; without zero ties, resolution is 2/S. A zero score means no
##'   contrary draw was observed, not certainty. Adjusted scores do not give
##'   the probability of a true null or guarantee frequentist FDR control.
##'   \code{aldex.p.adjust} does not support bootstrap SSRV fits.
##' @param contrasts optional numeric K x P matrix with one fixed-effect
##'   contrast per row. Its uniquely named columns must be the complete
##'   coefficient set and are reordered to the fitted design. Contrast rows
##'   have unique nonempty names, are finite and nonzero, and need not sum to
##'   zero or be independent. Contrasts require a full-rank fixed-effect design
##'   with positive residual degrees of freedom for existing t-based inference.
##'   The bootstrap path requires full rank but does not calculate residual
##'   degrees of freedom. Contrasts do not change composition
##'   or scale sampling or adjustment families. For BLMM, checks on an added
##'   contrast can cause further refinement or a refit, which can also change
##'   coefficient results slightly.
##' @param onesided (default \code{FALSE}) \code{FALSE} gives two-sided
##'   p-values; \code{"lower"} or \code{"upper"} gives one-sided p-values for
##'   the alternative that the estimate is negative or positive, respectively.
##'   Bootstrap SSRVs support only \code{FALSE}.
##' @param ... Additional scale-model arguments, such as \code{gamma},
##'   \code{s.mu}, or \code{c.mu}. For \code{method = "nlme"}, this may also
##'   include \code{random} or \code{correlation}.
##' @return An \code{aldex} list. \code{X} (P x N), \code{streaming},
##'   \code{method}, and the resolved \code{test} are always recorded.
##'   All inference records \code{p.adjust.method} and \code{onesided}.
##'   Bootstrap SSRVs also record \code{nsample},
##'   \code{scale.basis}, and \code{diagnostics}; their posterior outputs and
##'   defaults are described in Details. Unbounded/undefined Laplace approximations
##'   warn and record their conditional variance, including infinity. Affected
##'   posterior results are unavailable (NA), not replaced by bootstrap. Summary
##'   adds a status column for such fits; see the bootstrap SSRV vignette.
##'   Here P is
##'   the coefficient count, D the feature count, N the observation count, and
##'   S the Monte Carlo draw count. \code{return.pars} controls:
##'   \describe{
##'     \item{estimate}{A P x D x S array of coefficient draws.}
##'     \item{std.error, df}{For t-based inference, P x D x S arrays of
##'       within-fit standard errors and degrees of freedom.}
##'     \item{posterior}{For bootstrap SSRVs, P x D matrices of posterior means
##'       and equal-tailed 95% limits (mean, lower, upper), plus level and
##'       interval metadata. Requesting this also retains p.val and p.val.adj.}
##'     \item{p.val, p.val.adj}{P x D matrices of final raw and adjusted
##'       probabilities, respecting \code{onesided} and \code{p.adjust.method}.}
##'     \item{logComp, logScale}{D x N x S log2-composition draws and N x S
##'       log2-scale draws for t-based inference. These are not retained during
##'       streaming or supported for bootstrap SSRVs.}
##'     \item{random.effects}{Pr x D x S arrays of Pr random-effect variance
##'       and covariance terms, including residual variance, for mixed models.}
##'   }
##'   When contrasts are requested, \code{contrasts$matrix} always records their
##'   reordered matrix. \code{contrasts$estimate} has dimensions K x D x S,
##'   as do \code{contrasts$std.error} and \code{contrasts$df} for t-based
##'   inference. Bootstrap \code{contrasts$posterior} contains K x D summary
##'   matrices. \code{contrasts$p.val}
##'   and \code{contrasts$p.val.adj} have dimensions K x D. K is the number of
##'   contrast rows. The same \code{return.pars} fields control coefficient and
##'   contrast results together. Names and singleton dimensions are retained.
##'   Use \code{summary(fit, component = "contrasts")} for stored contrast
##'   summaries; the default summary includes coefficients and contrasts.
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
                  bootstrap.strata=NULL, bootstrap.max.attempts=100L,
                  ...) {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% c("lm", "blmm", "lme4", "nlme"))) {
    stop("method must be one of 'lm', 'blmm', 'lme4', or 'nlme'")
  }
  if (is.null(test)) {
    test <- switch(method, lm="t.HC3", blmm="t.satterthwaite.adaptive",
                   lme4="t", nlme="t")
  }
  if (!is.character(test) || length(test) != 1L || is.na(test)) {
    stop("test must be a single character inference selector")
  }
  if (test == "ssrv.bootstrap") {
    if (method != "lm") stop("test = 'ssrv.bootstrap' supports only method = 'lm'")
    if (inherits(X, "formula") && any(all.names(X) %in% c("|", "||")))
      stop("random-effects terms are not supported by test = 'ssrv.bootstrap'")
    args <- list(...)
    if (any(names(args) %in% c("random", "correlation", "weights", "cluster", "groups")))
      stop("random/correlation/cluster arguments are unsupported for ssrv.bootstrap")
    if (!identical(onesided, FALSE))
      stop("one-sided inference is unavailable for ssrv.bootstrap")
    p.adjust.method <- match.arg(p.adjust.method, stats::p.adjust.methods)
    if (missing(return.pars)) return.pars <- c("posterior", "estimate")
    return(.aldex_ssrv_bootstrap(Y, X, data, nsample, scale, streamsize,
                                return.pars, contrasts, bootstrap.strata,
                                bootstrap.max.attempts, p.adjust.method,
                                .upgrade_deprecated_scale_args(scale, args)))
  }
  if (inherits(scale, c("aldex_coefficient_scale", "aldex_sparse_scale", "aldex_gaussian_scale")))
    stop("composable coefficient scales require test = 'ssrv.bootstrap'")
  allowed.tests <- switch(method, lm=c("t", "t.HC0", "t.HC3"),
                          blmm=c("t.satterthwaite", "t.satterthwaite.adaptive", "t.residual"),
                          lme4="t", nlme="t")
  if (!(test %in% allowed.tests)) {
    stop("test '", test, "' is not supported for method '", method,
         "'; allowed selectors are: ", paste(allowed.tests, collapse=", "))
  }
  if (method == "blmm" && test == "t.residual") {
    warning("test = 't.residual' is deprecated; use 't.satterthwaite.adaptive' ",
            "or explicitly choose 't.satterthwaite' for shared-anchor inference",
            call.=FALSE)
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
                  inference=switch(test, t.satterthwaite="anchor",
                                   t.satterthwaite.adaptive="adaptive",
                                   t.residual="residual"),
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
  parameter.names <- rownames(X)
  label <- function(x, row.names) {
    rownames(x) <- row.names
    colnames(x) <- entity.names
    x
  }
  fitted.fields <- c("estimate", "std.error", "df")
  pval.fields <- c(p.val = "p.res", p.val.adj = "p.adj.res")

  res <- list()
  res$streaming <- stream
  res$X <- X
  res$method <- method
  res$test <- test
  res$p.adjust.method <- p.adjust.method
  res$onesided <- onesided
  ## Drop each array from out/pval.l once copied to keep peak memory down.
  for (field in fitted.fields) {
    if (field %in% return.pars) res[[field]] <- label(out[[field]], parameter.names)
    out[[field]] <- NULL
  }
  for (field in names(pval.fields)) {
    if (field %in% return.pars)
      res[[field]] <- label(pval.l[[pval.fields[[field]]]], parameter.names)
    pval.l[[pval.fields[[field]]]] <- NULL
  }
  if (!is.null(contrasts)) {
    cres <- list(matrix=contrasts)
    for (field in intersect(fitted.fields, return.pars))
      cres[[field]] <- label(out[[paste0("contrast.", field)]], rownames(contrasts))
    for (field in intersect(names(pval.fields), return.pars))
      cres[[field]] <- label(contrast.pval[[pval.fields[[field]]]], rownames(contrasts))
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
