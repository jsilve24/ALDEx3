## scale models

center <- function(logWpara) {
  g <- apply(logWpara, c(2,3), FUN=`mean`)
  return(sweep(logWpara, c(2,3), g, FUN=`-`))
 }

##' Default scale model discussed in Nixon et al. Beyond Normalizations
##'
##' @param gamma variance in scale, documented in Nixon et al. (TODO document
##'   here as well)
##' @return N x nsample matrix
##' @author Justin Silverman
clr <- function(X, logWpara, gamma=0.5) {
  P <- nrow(X)
  nsample <- dim(logWpara)[3]
  logWperp <- -colMeans(logWpara, dims=1)

  tmp <- P*nsample
  Lambdaperp <- matrix(rnorm(tmp,0,gamma), P, nsample)
  logWperp <- logWperp + t(X)%*%Lambdaperp
  return(logWperp)
}
