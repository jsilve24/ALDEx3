##' Sparse scale model for regression coefficients
##'
##' Estimate a total-abundance coefficient by assuming that the dominant peak
##' of absolute feature effects is at zero. For each draw, the estimated scale
##' coefficient is added to the corresponding log2-proportion coefficient
##' of every feature. Use this assumption only when a dominant cluster of
##' features with little or no absolute change is biologically plausible.
##'
##' Combine this specification with other coefficient models using
##' \code{compose.sm}, and supply the result as \code{scale} to
##' \code{aldex(test="ssrv.bootstrap")}. Assign every exact model-matrix
##' coefficient name, including the intercept.
##' @details The default \code{"hybrid"} compares bootstrap and Laplace
##' candidates using the median across features of their 95% final-effect
##' interval widths, and selects the complete sequence with the larger median;
##' exact ties select Laplace. Bootstrap resamples subjects and features;
##' Laplace fits the original subjects and uses inverse negative log-density
##' curvature for the scale variance, rather than a LASSO prior. When all
##' sparse specifications use Laplace, accompanying Gaussian coefficients also
##' use the original subjects, without subject or feature resampling.
##'
##' Gaussian kernel density estimation uses \code{bw.nrd0} by default,
##' re-estimating bandwidth for feature resamples, without a bandwidth floor.
##' Local maxima are refined and exact density ties choose the leftmost mode.
##' Zero curvature gives infinite Laplace variance; an unbounded or undefined
##' approximation warns and makes affected results unavailable, without an
##' automatic bootstrap substitution.
##' @param coefficients Exact model-matrix coefficient names.
##' @param uncertainty \code{"hybrid"} (default), \code{"bootstrap"}, or
##'   \code{"laplace"}; see Details.
##' @param bandwidth Optional finite positive Gaussian kernel SD in log2
##'   coefficient units. \code{NULL} estimates bandwidth separately for each
##'   profile. Constant profiles require an explicit positive bandwidth.
##' @return A sparse coefficient scale specification for \code{compose.sm}.
##' @export
sparse.sm <- function(coefficients, uncertainty=c("hybrid", "bootstrap", "laplace"),
                      bandwidth=NULL) {
  .ssrv_names(coefficients)
  uncertainty <- match.arg(uncertainty)
  if (!is.null(bandwidth) && (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
                            !is.finite(bandwidth) || bandwidth <= 0))
    stop("bandwidth must be finite and positive")
  structure(list(coefficients=coefficients, uncertainty=uncertainty,
                 bandwidth=bandwidth), class="aldex_sparse_scale")
}

##' Gaussian prior for total-abundance coefficients
##'
##' Specify a joint Gaussian prior on regression coefficients of log2 total
##' abundance. Each sampled coefficient is added to the corresponding
##' log2-proportion coefficient for every feature. For example, a treatment
##' coefficient with mean -1 represents an expected halving of total abundance.
##' A zero intercept coefficient fixes the reference total abundance at one in
##' the chosen units. Predictor recoding changes the prior unless its parameters
##' are transformed consistently.
##'
##' Use exact model-matrix names and combine specifications with \code{compose.sm}.
##' Coefficients in one call can be correlated through \code{covariance};
##' separate calls are independent. Covariance rows and columns follow the
##' order of \code{coefficients}.
##' @param coefficients Exact model-matrix coefficient names in covariance order.
##' @param mean Scalar or vector of prior means in log2-total coefficient units.
##' @param sd Scalar or vector of nonnegative prior standard deviations in the
##'   same units, supplied instead of \code{covariance}.
##' @param covariance Joint covariance matrix, supplied instead of \code{sd}.
##' @return A Gaussian coefficient scale specification for \code{compose.sm}.
##' @export
gaussian.sm <- function(coefficients, mean=0, sd=NULL, covariance=NULL) {
  .ssrv_names(coefficients)
  n <- length(coefficients)
  if (length(mean) == 1L) mean <- rep(mean, n)
  mean <- .validate_scale_mean(mean, n, "mean")
  if (!is.null(sd) && length(sd) == 1L) sd <- rep(sd, n)
  u <- .resolve_scale_uncertainty(sd, covariance, n=n,
                                 sd.name="sd", covariance.name="covariance")
  structure(list(coefficients=coefficients, mean=mean, sd=u$sd,
                 covariance=u$covariance), class="aldex_gaussian_scale")
}

##' Combine scale models for regression coefficients
##'
##' Assign a sparse model or Gaussian prior to each model-matrix coefficient,
##' including the intercept, exactly once. Pass the result to \code{aldex} as
##' \code{scale} with \code{test="ssrv.bootstrap"}. The models describe
##' log2 total-abundance coefficients added to fitted log2-proportion
##' coefficients; there is no automatic CLR centering.
##' @param ... Specifications made by \code{sparse.sm} or \code{gaussian.sm}.
##' @return A coefficient-scale specification for \code{aldex(scale=...)}.
##' @examples
##' sc <- compose.sm(gaussian.sm("(Intercept)", sd=0),
##'                  sparse.sm("treatment"), gaussian.sm("age", sd=0.1))
##' @export
compose.sm <- function(...) {
  blocks <- list(...)
  if (!length(blocks) || !all(vapply(blocks, function(b)
      inherits(b, "aldex_sparse_scale") || inherits(b, "aldex_gaussian_scale"), logical(1))))
    stop("compose.sm requires sparse.sm or gaussian.sm blocks")
  nm <- unlist(lapply(blocks, `[[`, "coefficients"), use.names=FALSE)
  if (anyDuplicated(nm)) stop("duplicate or overlapping coefficient assignments")
  structure(blocks, class="aldex_coefficient_scale")
}

.ssrv_names <- function(x) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop("coefficients must be unique nonempty exact coefficient names")
}

## Centering and bandwidth scaling here are algebraic coordinate changes only.
.ssrv_kde <- function(a, bandwidth=NULL) {
  if (!length(a) || any(!is.finite(a))) stop("KDE requires finite coefficients")
  origin <- a[1L]
  centered <- a-origin
  if (is.null(bandwidth) && (length(a)<2L || max(centered)==min(centered)))
    stop("KDE bandwidth cannot be estimated from a constant feature profile; supply an explicit positive bandwidth")
  h <- if (is.null(bandwidth)) bw.nrd0(centered) else bandwidth
  if (length(h)!=1L || !is.finite(h) || h<=0) stop("invalid KDE bandwidth")
  if (max(centered)==min(centered))
    return(list(mode=origin,variance=h^2,bandwidth=h))
  z <- centered/h
  if (any(!is.finite(z))) stop("KDE coordinates overflow at the requested bandwidth")
  weights <- function(x) {
    u <- -0.5*(x-z)^2
    w <- exp(u-max(u))
    w/sum(w)
  }
  logdensity <- function(x) {
    u <- -0.5*(x-z)^2
    m <- max(u)
    if (!is.finite(m)) return(-Inf)
    m+log(mean(exp(u-m)))
  }
  score <- function(x) sum(weights(x)*(z-x))
  # Target spacing is h/8 in original units. Kernel centers are also included.
  grid.points <- max(33L,ceiling(diff(range(z))*8)+1L)
  if(!is.finite(grid.points) || grid.points>.Machine$integer.max)
    stop("KDE grid is not representable at the requested bandwidth; no coarsening was applied")
  grid <- sort(unique(c(z,seq(min(z),max(z),length.out=grid.points))))
  val <- vapply(grid,logdensity,numeric(1))
  peaks <- which(val>=c(-Inf,utils::head(val,-1L)) &
                 val>=c(utils::tail(val,-1L),-Inf))
  candidates <- vapply(peaks,function(i) {
    lo <- if(i==1L) grid[i]-1 else grid[i-1L]
    hi <- if(i==length(grid)) grid[i]+1 else grid[i+1L]
    # Refine a stationary maximum using the analytic score where bracketed.
    # This also resolves the symmetric flat-mode fixture without a curvature cutoff.
    if (score(lo)>=0 && score(hi)<=0)
      uniroot(score,c(lo,hi),tol=.Machine$double.eps^0.75)$root
    else optimize(logdensity,c(lo,hi),maximum=TRUE,
                  tol=.Machine$double.eps^0.75)$maximum
  },numeric(1))
  density <- vapply(candidates,logdensity,numeric(1))
  best <- min(candidates[density==max(density)])
  w <- weights(best)
  gap <- 1-(sum(w*(z-best)^2)-sum(w*(z-best))^2)
  # No positive-curvature cutoff, variance cap, or replacement distribution.
  variance <- if (is.na(gap) || gap<0) NA_real_ else if(gap==0) Inf else h^2/gap
  list(mode=origin+h*best,variance=variance,bandwidth=h)
}

.ssrv_interval_width <- function(effects) {
  if (any(!is.finite(effects))) return(NA_real_)
  median(apply(effects, 1L, function(x) diff(quantile(x, c(0.025,0.975), names=FALSE))))
}

.ssrv_select <- function(laplace, bootstrap, mode) {
  if (mode!="hybrid") return(mode)
  vl <- .ssrv_interval_width(laplace); vb <- .ssrv_interval_width(bootstrap)
  if(is.na(vl) || is.na(vb)) return(NA_character_)
  if(vl>=vb) "laplace" else "bootstrap"
}
