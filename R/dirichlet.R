

##' Function for sampling Dirichlet random variables
##'
##' @param n number of samples
##' @param alpha D-vector of Dirichlet parameters
##' @return an D x n matrix of samples 
##' @author Justin Silverman
rDirichlet <- function(n, alpha) {
  D <- length(alpha)
  out <- matrix(NA, D, n)
  for (d in 1:D) {
    out[d,] <- rgamma(n, alpha[d],1)
  }
  return(miniclo(out))
 }
