## scale models

center <- function(logComp) {
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
clr <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  logScale <- -colMeans(logComp, dims=1)

  tmp <- P*nsample
  LambdaScale <- matrix(rnorm(tmp,0,gamma), P, nsample)
  logScale <- logScale + t(X)%*% LambdaScale
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
tss <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  LambdaScale <- matrix(rnorm(P*nsample,0,gamma), P, nsample)
  logScale <- t(X)%*% LambdaScale
  return(logScale)
}
