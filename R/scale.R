## scale models

center.sm <- function(logComp) {
  g <- apply(logComp, c(2,3), FUN=`mean`)
  return(sweep(logComp, c(2,3), g, FUN=`-`))
 }

##' Default scale model discussed in Nixon et al. Beyond Normalizations
##'
##' @param gamma variance in scale, documented in Nixon et al. (TODO document
##'   here as well)
##' @return N x nsample matrix
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

##' Scale model where the user may specify the mean and variance of the
##' log2 of the scale for each individual sample.
##'
##' Use this scale model if you want to specify the mean and (co)variance
##' of the log2 scale. This scale model is particularly useful if you have
##' paired qPCR or flow cytometry measurements of scale. Broadly useful if
##' you want to specify the scale for each individual sample.
##'
##' @param s.mu A vector of length N (# samples; ncol(X))
##'   representing the log2 mean of each sample's scale (cannot be NULL)
##' @param s.var A vector of length N (# samples; ncol(X))
##'   representing the variance of each sample in log2 space (cannot be
##'   NULL if s.cor is NULL)
##' @param s.cor A NxN vector where N is # of samples (ncol(X)) representing
##'   the covariance matrix in log2 space (cannot be NULL if s.var is NULL)
##' @return N x nsample matrix
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

##' Scale model where the user may specify the mean and variance for the
##' log2 fixed effects of the linear (mixed) model.
##'
##' For instance, if you have an intercept and a treatment where the
##' treatment is expected to increase the log2 scale, then we could
##' have c.mu=c(0, 1) and c.cor=diag(0.25, 2).
##'
##' @param s.mu A vector of length P (# covariates; nrow(X))
##'   representing the log2 mean of each sample's scale (cannot be NULL)
##' @param s.cor A PxP vector where P is # of covariates (nrow(X))
##'   represeting the covariance matrix in log2 space
##' @return N x nsample matrix
##' @importFrom MASS mvrnorm
##' @author Kyle McGovern
##' @export
coef.sm <- function(X, logComp, c.mu=NULL, c.cor=NULL) {
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

##' Scale model centered on TSS normalization and using log-normal to account
##' for uncertainty
##' 
##' This scale model is identical to the default scale model in Nixon et al.
##' Beyond Normalizations but is centered on the TSS normalization rather than
##' the CLR normalization. In practice this equates to an assumption that there
##' is no change in scale between conditions. The uncertainty in that assumption
##' is controled by the parameter gamma which is the standard deviation of a
##' normal distribution. 
##' 
##' Further documentation will be provided in future versions.  
##'
##' @param gamma variance in scale, documented in Nixon et al. (TODO document
##'   here as well)
##' @return N x nsample matrix
##' @author Justin Silverman
##' @export
tss.sm <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  LambdaScale <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
  logScale <- t(X)%*% LambdaScale
  return(logScale)
}
