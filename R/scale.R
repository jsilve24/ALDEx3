## scale models

.validate_scale_mean <- function(x, n, name) {
  if (is.null(x)) {
    stop(name, " cannot be NULL")
  }
  if (!is.numeric(x) || length(x) != n) {
    stop(name, " should be a numeric vector of length ", n)
  }
  if (!all(is.finite(x))) {
    stop(name, " should contain only finite values")
  }
  x
}

.validate_scale_sd <- function(x, n, name) {
  if (!is.numeric(x) || length(x) != n) {
    stop(name, " should be a numeric vector of length ", n)
  }
  if (!all(is.finite(x))) {
    stop(name, " should contain only finite values")
  }
  if (any(x < 0)) {
    stop(name, " should contain only non-negative values")
  }
  x
}

.validate_scale_variance <- function(x, n, name) {
  if (!is.numeric(x) || length(x) != n) {
    stop(name, " should be a numeric vector of length ", n)
  }
  if (!all(is.finite(x))) {
    stop(name, " should contain only finite values")
  }
  if (any(x < 0)) {
    stop(name, " should contain only non-negative values")
  }
  x
}

.validate_scale_covariance <- function(x, n, name) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(n, n))) {
    stop(name, " should be a numeric ", n, "x", n, " covariance matrix")
  }
  if (!all(is.finite(x))) {
    stop(name, " should contain only finite values")
  }

  magnitude <- max(1, max(abs(x)))
  if (max(abs(x - t(x))) > sqrt(.Machine$double.eps) * magnitude) {
    stop(name, " should be symmetric")
  }
  x <- (x + t(x)) / 2
  values <- eigen(x, symmetric=TRUE, only.values=TRUE)$values
  if (min(values) < -1e-8 * max(1, max(abs(values)))) {
    stop(name, " should be positive semidefinite")
  }
  x
}

.draw_scale_mvn <- function(nsample, mu, covariance) {
  draws <- MASS::mvrnorm(nsample, mu, covariance)
  t(matrix(draws, nrow=nsample, ncol=length(mu)))
}

.resolve_scale_uncertainty <- function(sd, covariance,
                                       legacy.variance=NULL,
                                       legacy.covariance=NULL,
                                       n, sd.name, covariance.name,
                                       legacy.variance.name=NULL,
                                       legacy.covariance.name=NULL,
                                       warn=TRUE) {
  supplied <- c(!is.null(sd), !is.null(covariance),
                !is.null(legacy.variance), !is.null(legacy.covariance))
  if (sum(supplied) != 1) {
    stop("Exactly one of ", sd.name, ", ", covariance.name,
         if (!is.null(legacy.variance.name)) {
           paste0(", ", legacy.variance.name)
         } else "",
         if (!is.null(legacy.covariance.name)) {
           paste0(", or ", legacy.covariance.name)
         } else "",
         " should be provided")
  }

  if (!is.null(legacy.variance)) {
    legacy.variance <- .validate_scale_variance(
      legacy.variance, n, legacy.variance.name
    )
    if (warn) {
      warning(legacy.variance.name, " is deprecated; use ", sd.name,
              " = sqrt(", legacy.variance.name, ") instead",
              call.=FALSE)
    }
    sd <- sqrt(legacy.variance)
  } else if (!is.null(legacy.covariance)) {
    legacy.covariance <- .validate_scale_covariance(
      legacy.covariance, n, legacy.covariance.name
    )
    if (warn) {
      warning(legacy.covariance.name, " is deprecated; use ",
              covariance.name, " instead", call.=FALSE)
    }
    covariance <- legacy.covariance
  }

  if (!is.null(sd)) {
    sd <- .validate_scale_sd(sd, n, sd.name)
  } else {
    covariance <- .validate_scale_covariance(
      covariance, n, covariance.name
    )
  }
  list(sd=sd, covariance=covariance)
}

.upgrade_deprecated_scale_args <- function(scale, args) {
  if (identical(scale, sample.sm)) {
    supplied <- c(!is.null(args$s.sd), !is.null(args$s.cov),
                  !is.null(args$s.var), !is.null(args$s.cor))
    if (sum(supplied) > 1) {
      stop("Exactly one of s.sd, s.cov, s.var, or s.cor should be provided")
    }
    if (!is.null(args$s.var)) {
      args$s.var <- .validate_scale_variance(
        args$s.var, length(args$s.var), "s.var"
      )
      warning("s.var is deprecated; use s.sd = sqrt(s.var) instead",
              call.=FALSE)
      args$s.sd <- sqrt(args$s.var)
      args$s.var <- NULL
    } else if (!is.null(args$s.cor)) {
      warning("s.cor is deprecated; use s.cov instead", call.=FALSE)
      args$s.cov <- args$s.cor
      args$s.cor <- NULL
    }
  } else if (identical(scale, coefficient.sm)) {
    supplied <- c(!is.null(args$c.sd), !is.null(args$c.cov),
                  !is.null(args$c.cor))
    if (sum(supplied) > 1) {
      stop("Exactly one of c.sd, c.cov, or c.cor should be provided")
    }
    if (!is.null(args$c.cor)) {
      warning("c.cor is deprecated; use c.cov instead", call.=FALSE)
      args$c.cov <- args$c.cor
      args$c.cor <- NULL
    }
  }
  args
}

center.sm <- function(logComp) {
  g <- apply(logComp, c(2,3), FUN=`mean`)
  return(sweep(logComp, c(2,3), g, FUN=`-`))
 }

##' CLR-based scale model
##'
##' Use this model when centered log-ratio (CLR) normalization is a reasonable
##' starting point, but you do not want to assume that it recovers the sample
##' scales exactly. For example, in a case-control study, CLR may provide a
##' plausible center while \code{gamma} allows the true difference in total
##' abundance between cases and controls to depart from that center.
##'
##' Set \code{gamma = 0} to use the CLR-implied scales without additional
##' uncertainty. A value such as \code{gamma = 0.5} allows each scale-model
##' coefficient to vary with SD 0.5 on the log2 scale. Larger values express
##' less confidence in the CLR assumption.
##'
##' The uncertainty applies to model coefficients, not independently to every
##' sample. Samples with the same covariate values therefore receive the same
##' random shift in a Monte Carlo draw. This differs from \code{\link{sample.sm}},
##' which can draw a separate scale for each sample.
##'
##' @param X Model matrix passed automatically by \code{aldex()}. Rows are model
##'   coefficients and columns are samples.
##'
##' @param logComp Monte Carlo log-compositions passed automatically by
##'   \code{aldex()}, with dimensions \code{features x samples x nsample}.
##'
##' @param gamma Non-negative scalar. Standard deviation of the Gaussian random
##'   shift that relaxes the CLR assumption about scale. \code{gamma = 0}
##'   yields the fixed CLR assumption. Larger values express less confidence in
##'   the CLR-implied differences between sample scales.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo samples
##'   of the log-scale for each sample (rows) and each Monte Carlo draw (columns).
##'
##' @references Nixon G, Gloor GB, Silverman JD (2025).
##'   "Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization".
##'   Genome Biology. \doi{10.1186/s13059-025-03609-3}
##' 
##' @author Justin Silverman
##' @examples
##' # Allow moderate uncertainty around CLR-implied scale differences.
##' Y <- matrix(seq_len(60), nrow = 10)
##' condition <- factor(c("control", "control", "control",
##'                       "case", "case", "case"))
##' fit <- aldex(Y, ~ condition, data.frame(condition), nsample = 20,
##'              scale = clr.sm, gamma = 0.5)
##' @export
clr.sm <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  gamma <- .validate_scale_sd(gamma, 1, "gamma")
  logScale <- -colMeans(logComp, dims=1)

  tmp <- P*nsample
  LambdaScale <- matrix(rnorm(tmp,0,gamma), P, nsample)
  logScale <- logScale + t(X)%*% LambdaScale
  return(logScale)
}

##' Per-sample scale model using external information
##'
##' Use this model when each sample has its own external scale measurement, such
##' as a flow-cytometry cell count, qPCR concentration, or spike-in estimate.
##' Put the log2 measurement for each sample in \code{s.mu}. Use \code{s.sd} to
##' describe the measurement uncertainty for each sample.
##'
##' For example, if four samples have estimated total cell counts of
##' \code{1e8}, \code{2e8}, \code{8e7}, and \code{1.5e8}, set
##' \code{s.mu = log2(c(1e8, 2e8, 8e7, 1.5e8))}. If each measurement has an SD
##' of 0.5 on the log2 scale, set \code{s.sd = rep(0.5, 4)}.
##'
##' Use \code{s.sd} when measurement errors are independent. Use \code{s.cov}
##' instead when errors are correlated, for
##' example because several samples share a calibration standard. The diagonal
##' of \code{s.cov} contains variances, and its off-diagonal entries contain
##' covariances. Thus independent SDs of 0.5 can also be written as
##' \code{s.cov = diag(0.5^2, 4)}.
##'
##' Supply exactly one of \code{s.sd} or \code{s.cov}. With \code{s.sd}, each
##' sample receives its own independent scale draw. Samples in the same group
##' do not automatically share uncertainty as they do in coefficient-based
##' models.
##'
##' @param X Model matrix passed automatically by \code{aldex()}. This function
##'   uses its number of columns to determine the number of samples.
##'
##' @param logComp Monte Carlo log-compositions passed automatically by
##'   \code{aldex()}. This function uses the third dimension to determine the
##'   number of draws.
##'
##' @param s.mu Expected log2 scale for each sample, in the same order as the
##'   columns of the count matrix. Must have length \code{N}.
##'
##' @param s.var Deprecated numeric vector of per-sample log2-scale variances.
##'   Use \code{s.sd = sqrt(s.var)} instead. Retained until ALDEx3 2.0.0 for
##'   backward compatibility.
##'
##' @param s.cor Deprecated \code{N x N} covariance matrix. Despite its name,
##'   this argument has always represented covariance rather than correlation.
##'   Use \code{s.cov} instead. Retained until ALDEx3 2.0.0.
##'
##' @param s.sd SD of the log2 scale for each sample. Must have length \code{N}.
##'   Use this when measurement errors are independent across samples.
##'
##' @param s.cov An \code{N x N} covariance matrix for the samples' log2 scales.
##'   Use this instead of \code{s.sd} when measurement errors are correlated.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo
##'   draws of the log2 scale for each sample (rows) across \code{nsample} draws
##'   (columns).
##'
##' @importFrom MASS mvrnorm
##' @author Kyle McGovern
##' @examples
##' # External total-cell estimates for four samples.
##' cell_count <- c(1e8, 2e8, 8e7, 1.5e8)
##' X <- matrix(1, 1, 4)
##' logComp <- array(0, c(2, 4, 10))
##' logScale <- sample.sm(
##'   X, logComp, s.mu = log2(cell_count), s.sd = rep(0.5, 4)
##' )
##' @export
sample.sm <- function(X, logComp, s.mu=NULL, s.var=NULL, s.cor=NULL,
                      s.sd=NULL, s.cov=NULL) {
  N <- ncol(X)
  nsample <- dim(logComp)[3]
  s.mu <- .validate_scale_mean(s.mu, N, "s.mu")
  uncertainty <- .resolve_scale_uncertainty(
    sd=s.sd, covariance=s.cov,
    legacy.variance=s.var, legacy.covariance=s.cor,
    n=N, sd.name="s.sd", covariance.name="s.cov",
    legacy.variance.name="s.var", legacy.covariance.name="s.cor"
  )

  if (!is.null(uncertainty$sd)) {
    logScale <- matrix(
      rnorm(N * nsample, rep(s.mu, nsample),
            rep(uncertainty$sd, nsample)),
      N, nsample
    )
  } else {
    logScale <- .draw_scale_mvn(nsample, s.mu, uncertainty$covariance)
  }

  return(logScale)
}


##' Scale model using prior information about groups or treatments
##'
##' Use this model when prior information describes a treatment, group, time, or
##' other model coefficient rather than an individual sample. Samples with the
##' same covariate values share the same scale shift in each Monte Carlo draw.
##'
##' Suppose an effect plot from a control-versus-treatment experiment is
##' centered near -8, even though prior knowledge suggests that most features
##' should not change. One way to represent that knowledge is to shift the
##' treatment samples upward by 8 log2 units relative to the controls. With
##' control as the reference level, use \code{c.mu = c(0, 8)}: the first value
##' leaves the common intercept unchanged, and the second adds 8 to treatment.
##' Setting \code{c.sd = c(0, 0.5)} expresses an SD of 0.5 around that treatment
##' shift. If the plot were centered near +8, use -8 instead. Adding 8 to both
##' groups would not recenter the plot because it would leave their difference
##' unchanged.
##'
##' More generally, \code{c.mu} gives the expected log2 scale change for each
##' model coefficient and \code{c.sd} gives its SD. Use \code{c.cov} when the
##' coefficient uncertainties are correlated. For each Monte Carlo draw,
##' ALDEx3 samples a coefficient vector and applies it to the model matrix
##' \code{X}; mathematically, \eqn{b^{(m)} \sim N(c.mu, Sigma)} and the sample
##' scales are \eqn{b^{(m)T} X}.
##'
##' Exactly one of \code{c.sd} or \code{c.cov} must be provided. The deprecated
##' \code{c.cor} argument may be used in place of \code{c.cov} during the
##' compatibility period.
##'
##' @param X Model matrix passed automatically by \code{aldex()}. Rows are model
##'   coefficients and columns are samples.
##'
##' @param logComp Monte Carlo log-compositions passed automatically by
##'   \code{aldex()}. This function uses the third dimension to determine the
##'   number of draws.
##'
##' @param c.mu Expected log2 scale change for each model coefficient, in the
##'   same order as the rows of \code{X}. Must have length \code{P}.
##'
##' @param c.cor Deprecated \code{P x P} covariance matrix for the fixed-effect
##'   coefficients. Use \code{c.cov} instead. Retained until ALDEx3 2.0.0.
##'
##' @param c.sd SD of the log2 scale change for each model coefficient. Must
##'   have length \code{P}; coefficients are treated as independent.
##'
##' @param c.cov A \code{P x P} covariance matrix for the model coefficients'
##'   log2 scale changes. Use this instead of \code{c.sd} when coefficient
##'   uncertainties are correlated.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo
##'   draws of the log2 scale for each sample (rows) across \code{nsample} draws
##'   (columns).
##'
##' @importFrom MASS mvrnorm
##' @author Kyle McGovern
##' @examples
##' # Three control samples followed by three treatment samples.
##' condition <- factor(
##'   c("control", "control", "control", "treatment", "treatment", "treatment"),
##'   levels = c("control", "treatment")
##' )
##' X <- t(model.matrix(~ condition))
##' logComp <- array(0, c(2, 6, 10))
##' logScale <- coefficient.sm(
##'   # Shift treatment upward by 8, with SD 0.5 around that shift.
##'   X, logComp, c.mu = c(0, 8), c.sd = c(0, 0.5)
##' )
##' @export
coefficient.sm <- function(X, logComp, c.mu=NULL, c.cor=NULL,
                           c.sd=NULL, c.cov=NULL) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]

  c.mu <- .validate_scale_mean(c.mu, P, "c.mu")
  uncertainty <- .resolve_scale_uncertainty(
    sd=c.sd, covariance=c.cov, legacy.covariance=c.cor,
    n=P, sd.name="c.sd", covariance.name="c.cov",
    legacy.covariance.name="c.cor"
  )

  if (!is.null(uncertainty$sd)) {
    coefficients <- matrix(
      rnorm(P * nsample, rep(c.mu, nsample), rep(uncertainty$sd, nsample)),
      P, nsample
    )
  } else {
    coefficients <- matrix(NA_real_, P, nsample)
    for (i in seq_len(nsample)) {
      coefficients[, i] <- MASS::mvrnorm(
        1, c.mu, uncertainty$covariance
      )
    }
  }
  logScale <- t(X) %*% coefficients
  return(logScale)
}

##' TSS-centered scale model
##'
##' Total-sum scaling (TSS) starts from the assumption that total scale does not
##' change systematically between experimental groups. Use this model when that
##' is a reasonable starting point but you want to express uncertainty about it.
##' For example, if a treatment is expected to redistribute features without
##' changing total biomass, TSS centers the treatment-versus-control scale
##' difference at zero.
##'
##' Set \code{gamma = 0} to enforce no scale difference. Setting
##' \code{gamma = 0.5} keeps the expected difference at zero but allows each
##' scale-model coefficient to vary with SD 0.5 on the log2 scale. Larger values
##' express less confidence that total scale is unchanged.
##'
##' ALDEx3 applies this uncertainty to model coefficients: for each draw,
##' \eqn{b^{(m)} \sim N(0, \gamma^2 I)}, and the sample scales are
##' \eqn{b^{(m)T} X}. Samples with the same covariate values therefore share a
##' random shift. This differs from \code{\link{sample.sm}}, which can draw
##' independent uncertainty for every sample.
##'
##' @param X Model matrix passed automatically by \code{aldex()}. Rows are model
##'   coefficients and columns are samples.
##'
##' @param logComp Monte Carlo log-compositions passed automatically by
##'   \code{aldex()}. This function uses the third dimension to determine the
##'   number of draws.
##'
##' @param gamma Non-negative scalar. Standard deviation of the Gaussian random
##'   shift applied to the scale-model coefficients (in log2 space).
##'   \code{gamma = 0} implies no scale uncertainty (all draws are centered at
##'   zero effect); larger values allow greater departures from the TSS-centered
##'   assumption.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo
##'   draws of the log2 scale for each sample (rows) across \code{nsample} draws
##'   (columns).
##'
##' @references Nixon G, Gloor GB, Silverman JD (2025).
##'   "Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization".
##'   Genome Biology. \doi{10.1186/s13059-025-03609-3}
##'
##' @author Justin Silverman
##' @examples
##' # Center the group-scale difference at zero, with moderate uncertainty.
##' Y <- matrix(seq_len(60), nrow = 10)
##' condition <- factor(c("control", "control", "control",
##'                       "treatment", "treatment", "treatment"))
##' fit <- aldex(Y, ~ condition, data.frame(condition), nsample = 20,
##'              scale = tss.sm, gamma = 0.5)
##' @export
tss.sm <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  gamma <- .validate_scale_sd(gamma, 1, "gamma")
  LambdaScale <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
  logScale <- t(X)%*% LambdaScale
  return(logScale)
}
