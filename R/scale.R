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
clr <- function(X, logComp, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logComp)[3]
  logScale <- -colMeans(logComp, dims=1)

  tmp <- P*nsample
  LambdaScale <- matrix(rnorm(tmp,0,gamma), P, nsample)
  logScale <- logScale + t(X)%*% LambdaScale
  return(logScale)
}
