## scale models

center.sm <- function(logComp) {
  g <- apply(logComp, c(2,3), FUN=`mean`)
  return(sweep(logComp, c(2,3), g, FUN=`-`))
 }

##' Default CLR-based scale model (with optional scale uncertainty)
##'
##' Implements the default scale model described in Nixon et al. (Beyond
##' Normalizations / scale-uncertainty framework). This model generalizes the
##' centered log-ratio (CLR) normalization by treating the (log) scale as a
##' latent random variable and allowing additive uncertainty around the
##' CLR-implied scale differences via a Gaussian term with standard deviation
##' \code{gamma}.
##'
##' In the limit \code{gamma = 0}, this reduces to the CLR assumption (no scale
##' uncertainty beyond the CLR-implied scale). Larger \code{gamma} values
##' represent increasing uncertainty about the CLR-implied scale differences.
##'
##' @param X A numeric design matrix used to model scale variation across
##'   samples. This is the covariate/design matrix passed internally by
##'   \code{aldex()} to the scale model. Rows correspond to regression
##'   coefficients (e.g., intercept and covariates after contrasts/encoding) and
##'   columns correspond to samples. If the analysis includes only an intercept
##'   (no covariates), \code{X} is typically a 1 x N matrix of ones. (This
##'   parameter is automatically passed by aldex)
##'
##' @param logComp A numeric array of Monte Carlo log-compositions with
##'   dimensions \code{features x samples x nsample}. This is produced
##'   internally by ALDEx3 from Dirichlet-multinomial Monte Carlo sampling and
##'   log-ratio representation. ##'   (This parameter is automatically passed by aldex)
##'
##' @param gamma Non-negative scalar. Standard deviation of the Gaussian
##'   perturbation that relaxes the CLR assumption about scale. \code{gamma = 0}
##'   yields the pure CLR assumption; recommended default values in the
##'   scale-uncertainty literature are often around \code{0.5}, but appropriate
##'   values depend on how strongly you trust the CLR scale assumption in the
##'   current study.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo samples
##'   of the log-scale for each sample (rows) and each Monte Carlo draw (columns).
##'
##' @references Nixon, Gloor, and Silverman. Incorporating scale uncertainty in
##'   microbiome and gene expression analysis as an extension of normalization.
##'   Genome Biology (2025).
##'
##' @author Justin Silverman
##' @export
clr.sm <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  logScale <- -colMeans(logComp, dims=1)

  tmp <- P*nsample
  LambdaScale <- matrix(rnorm(tmp,0,gamma), P, nsample)
  logScale <- logScale + t(X)%*% LambdaScale
  return(logScale)
}

##' Sample-specific scale model with user-specified mean and variance/covariance
##'
##' Draws Monte Carlo samples of the log2 scale for each sample using user-supplied
##' moments. This scale model is useful when external measurements (e.g., qPCR,
##' flow cytometry, spike-ins) provide information about absolute scale, or when
##' you want to encode prior information about scale on a per-sample basis.
##'
##' Exactly one of \code{s.var} or \code{s.cor} must be provided:
##' \itemize{
##'   \item \code{s.var}: independent per-sample log2-scale variance (diagonal covariance)
##'   \item \code{s.cor}: full \code{N x N} log2-scale covariance matrix
##' }
##'
##' The returned matrix has \code{N} rows (samples) and \code{nsample} columns
##' (Monte Carlo draws), consistent with the ALDEx3 scale-model interface.
##'
##' @param X A numeric design matrix passed internally by \code{aldex()} to the
##'   scale model. Columns correspond to samples (\code{N = ncol(X)}). This scale
##'   model does not use \code{X} directly, but \code{N} is inferred from it.
##'   (Automatically supplied by \code{aldex()}.)
##'
##' @param logComp A numeric array of Monte Carlo log-compositions with dimensions
##'   \code{features x samples x nsample}. This scale model uses \code{nsample}
##'   to determine the number of Monte Carlo draws, but does not otherwise use
##'   \code{logComp}. (Automatically supplied by \code{aldex()}.)
##'
##' @param s.mu Numeric vector of length \code{N} giving the mean of the log2
##'   scale for each sample. Must not be \code{NULL}.
##'
##' @param s.var Numeric vector of length \code{N} giving the marginal variance of
##'   the log2 scale for each sample. Use this when assuming samples' log2 scales
##'   are independent. Must be \code{NULL} if \code{s.cor} is provided.
##'
##' @param s.cor Numeric \code{N x N} covariance matrix for the log2 scale across
##'   samples. Use this when encoding correlations between samples' log2 scales.
##'   Must be \code{NULL} if \code{s.var} is provided.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo
##'   draws of the log2 scale for each sample (rows) across \code{nsample} draws
##'   (columns).
##'
##' @importFrom MASS mvrnorm
##' @author Kyle McGovern
##' @export
sample.sm <- function(X, logComp, s.mu=NULL, s.var=NULL, s.cor=NULL) {
  N <- ncol(X)
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  if(is.null(s.mu)) {
    stop("s.mu cannot be NULL")
  }
  if(length(s.mu)!=N) {
    stop("s.mu should have same length as ncol(X)")
  }

  if((!is.null(s.var))&is.null(s.cor)) {
    if(length(s.var)!=N) {
      stop("s.var should have same length as ncol(X)")
    }
    logScale <- replicate(nsample, rnorm(N, s.mu, sqrt(s.var)))
  } else if((!is.null(s.cor))&is.null(s.var)) {
    if(any(dim(s.cor) != c(N, N))) {
      stop("s.cor should be an NxN matrix where N=ncol(X)")
    }
    logScale <- t(mvrnorm(nsample, s.mu, s.cor))
  } else {
    stop("One of s.var and s.cor should be NULL")
  }

  return(logScale)
}


##' Coefficient-based scale model with user-specified prior on fixed effects
##'
##' Draws Monte Carlo samples of the log2 scale by sampling fixed-effect
##' coefficients from a multivariate normal distribution and mapping them
##' through the design matrix \code{X}. This scale model is useful when you
##' want to encode prior information about how covariates (e.g., treatment,
##' batch, time) affect scale, rather than specifying scale moments directly
##' per sample.
##'
##' Specifically, for each Monte Carlo draw \eqn{b^{(m)} \sim N(c.mu, c.cor)},
##' the per-sample log2 scale is computed as \eqn{b^{(m)T} X}, producing an
##' \code{N x nsample} matrix of log2-scale draws.
##'
##' For example, with an intercept and a treatment indicator where treatment is
##' expected to increase log2 scale by ~1 on average, one might use
##' \code{c.mu = c(0, 1)} and \code{c.cor = diag(c(0.25, 0.25))} (i.e., SD 0.5
##' for each coefficient, independent).
##'
##' @param X A numeric design matrix passed internally by \code{aldex()} to the
##'   scale model. Rows correspond to fixed-effect coefficients/covariates
##'   (\code{P = nrow(X)}) and columns correspond to samples
##'   (\code{N = ncol(X)}). (Automatically supplied by \code{aldex()}.)
##'
##' @param logComp A numeric array of Monte Carlo log-compositions with
##'   dimensions \code{features x samples x nsample}. This scale model uses
##'   \code{nsample} (the number of Monte Carlo draws) but does not otherwise use
##'   \code{logComp}. (Automatically supplied by \code{aldex()}.)
##'
##' @param c.mu Numeric vector of length \code{P} giving the mean of the fixed
##'   effect coefficients in log2-scale space. Must not be \code{NULL}.
##'
##' @param c.cor Numeric \code{P x P} covariance matrix for the fixed effect
##'   coefficients in log2-scale space. Must not be \code{NULL}.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo
##'   draws of the log2 scale for each sample (rows) across \code{nsample} draws
##'   (columns).
##'
##' @importFrom MASS mvrnorm
##' @author Kyle McGovern
##' @export
coefficient.sm <- function(X, logComp, c.mu=NULL, c.cor=NULL) {
  N <- ncol(X)
  P <- nrow(X)
  nsample <- dim(logComp)[3]

  if(is.null(c.mu)) {
    stop("c.mu cannot be NULL")
  } else if(length(c.mu)!=P) {
    stop("c.mu should have length of P=nrow(X)")
  }
  if(is.null(c.cor)) {
    stop("c.cor cannot be NULL")
  } else if(any(dim(c.cor)!=c(P, P))) {
    stop("c.cor should be a PxP matrix where P=nrow(X)")
  }

  logScale <- matrix(NA, N, nsample)
  for(i in 1:nsample) {
    logScale[,i] <- mvrnorm(1, c.mu, c.cor)%*%X
  }
  return(logScale)
}

##' TSS-centered scale model (with optional scale uncertainty)
##'
##' Implements a total-sum-scaling (TSS)-centered variant of the default
##' scale-uncertainty model described in Nixon et al. (Beyond Normalizations /
##' scale-uncertainty framework). Unlike \code{\link{clr.sm}}, which is centered
##' on the CLR-implied scale, this model is centered on the TSS assumption that
##' there is no systematic change in scale across.
##'
##' Scale uncertainty is introduced via an additive Gaussian perturbation on the
##' (log2) fixed effects. For each Monte Carlo draw, a coefficient vector is
##' sampled as \eqn{b^{(m)} \sim N(0, \gamma^2 I)}, and the per-sample log2 scale
##' is computed as \eqn{b^{(m)T} X}. Larger values of \code{gamma} correspond to
##' weaker confidence in the TSS-centered assumption (more allowed scale
##' variation); \code{gamma = 0} yields no scale variation beyond the model
##' center.
##'
##' Note: \code{logComp} is included to match the ALDEx3 scale-model interface
##' and to determine \code{nsample}, but it is not otherwise used by this model.
##'
##' @param X A numeric design matrix passed internally by \code{aldex()} to the
##'   scale model. Rows correspond to fixed-effect coefficients/covariates
##'   (\code{P = nrow(X)}) and columns correspond to samples
##'   (\code{N = ncol(X)}). (Automatically supplied by \code{aldex()}.)
##'
##' @param logComp A numeric array of Monte Carlo log-compositions with
##'   dimensions \code{features x samples x nsample}. This scale model uses
##'   \code{nsample} (the number of Monte Carlo draws) but does not otherwise use
##'   \code{logComp}. (Automatically supplied by \code{aldex()}.)
##'
##' @param gamma Non-negative scalar. Standard deviation of the Gaussian
##'   perturbation applied to the scale-model coefficients (in log2 space).
##'   \code{gamma = 0} implies no scale uncertainty (all draws are centered at
##'   zero effect); larger values allow greater departures from the TSS-centered
##'   assumption.
##'
##' @return A numeric matrix of dimension \code{N x nsample} giving Monte Carlo
##'   draws of the log2 scale for each sample (rows) across \code{nsample} draws
##'   (columns).
##'
##' @references Nixon, Gloor, and Silverman. Incorporating scale uncertainty in
##'   microbiome and gene expression analysis as an extension of normalization.
##'   Genome Biology (2025).
##'
##' @author Justin Silverman
##' @export
tss.sm <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  LambdaScale <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
  logScale <- t(X)%*% LambdaScale
  return(logScale)
}
