##' Sparse coefficient scale model
##'
##' Assign a shared scalar scale shift to named log2-proportion coefficients. Gaussian KDE
##' modes use bw.nrd0 without a bandwidth floor, refined local
##' maxima, and leftmost ties. Feature bootstrap re-estimates bandwidth.
##' Laplace uses inverse negative log-density curvature, not a LASSO prior.
##' Hybrid selects one complete candidate sequence by its all-draw variance
##' (bootstrap on exact ties), without substituting a candidate on numerical failure.
##' The sparse assumption places the dominant absolute-effect density peak at
##' zero; the conditional model estimates total-scale coefficients from the
##' compositional feature profile. Zero log-density curvature gives infinite
##' Laplace variance. Unbounded or undefined approximations are warned about;
##' affected posterior results are unavailable rather than replaced by bootstrap.
##' @param coefficients Exact model-matrix coefficient names.
##' @param uncertainty One of hybrid, bootstrap, or laplace.
##' @param bandwidth Optional finite positive Gaussian kernel SD in log2
##'   coefficient units. NULL estimates bandwidth separately for each profile.
##'   Constant profiles require an explicit positive bandwidth; no smoothing
##'   variance is assigned automatically.
##' @return A coefficient scale block for compose.sm.
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

##' Gaussian coefficient scale block
##'
##' A joint Gaussian prior on regression coefficients of log2 absolute total
##' abundance. Each draw is added to the corresponding log2-proportion coefficient
##' for every feature. Separate blocks are independent; covariance within a block
##' is preserved. Assign the intercept explicitly, for example with mean=0, sd=0.
##' Zero shifts give the log-proportion coefficients without correction: a zero
##' intercept scale shift fixes reference total abundance at one, not its
##' geometric mean. These are total-scale coefficients, as in coefficient.sm.
##' Recoding predictors changes these priors unless transformed consistently.
##' @param coefficients Exact model-matrix coefficient names in block order.
##' @param mean Scalar or vector of expected scale shifts.
##' @param sd Scalar or vector of nonnegative standard deviations.
##' @param covariance Joint covariance matrix, instead of sd.
##' @return A coefficient scale block for compose.sm.
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

##' Compose named coefficient scale blocks
##'
##' Every coefficient, including the intercept, must be assigned exactly once.
##' Composable models use log2 proportions as baseline; no sample-scale model
##' is added. Only test="ssrv.bootstrap" supports these specifications.
##' @param ... Blocks made by sparse.sm or gaussian.sm.
##' @return A composable coefficient-scale specification for aldex(scale=...).
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

.ssrv_select <- function(laplace, bootstrap, mode) {
  if (mode!="hybrid") return(mode)
  # An unavailable sequence has no sample variance to compare. Never substitute
  # the other candidate or invent an empirical variance for an unsampled Gaussian.
  if (any(!is.finite(laplace)) || any(!is.finite(bootstrap))) return(NA_character_)
  vl <- var(laplace); vb <- var(bootstrap)
  if(is.na(vl) || is.na(vb)) return(NA_character_)
  if(vl>vb) "laplace" else "bootstrap"
}
