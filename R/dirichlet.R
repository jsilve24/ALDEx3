##' Function for sampling Dirichlet random variables
##'
##' @param n number of samples
##' @param Alpha DxN matrix of Dirichlet parameters
##' @return an D x n matrix of samples 
##' @author Justin Silverman
rDirichletMat <- function(n, Alpha) {
  D <- nrow(Alpha)
  N <- ncol(Alpha)
  logWpara <- array(rep(Alpha, times=n), dim=c(D, N, n))
  logWpara[] <- rgamma(N*D*n, logWpara[], 1)
  s <- colSums(logWpara, dims=1)
  sweep(logWpara, c(2,3), s, FUN=`/`)
 }
