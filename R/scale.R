## scale models

center <- function(logWpara) {
  g <- apply(logWpara, c(2,3), FUN=`mean`)
  return(sweep(logWpara, c(2,3), g, FUN=`-`))
 }

##' Default scale model discussed in Nixon et al. Beyond Normalizations
##'
##' @param gamma TODO fix documentation
##' @return N x nsample matrix
##' @author Justin Silverman
default <- function(X, Y, logWpara) {
  gamma <- 0.5 # TODO provide functionality so this is no longer hard coded
  P <- nrow(X)
  nsample <- dim(logWpara)[3]
  logWperp <- -apply(logWpara, c(2,3), FUN=`mean`)
  tmp <- P*nsample
  Lambdaperp <- matrix(rnorm(tmp,0,gamma), P, nsample)
  logWperp <- logWperp + t(X)%*% Lambdaperp
  return(logWperp)
}
