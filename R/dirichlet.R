##' Function for sampling Dirichlet random variables (base-2 normalized via log-sum-exp for stability)
##'
##' @param n number of samples
##' @param alpha D-vector of Dirichlet parameters
##' @return a D x n matrix of samples 
##' @author Justin Silverman
rLogDirichlet <- function(n, alpha) {
  D <- length(alpha)
  out <- matrix(NA, D, n)
  for (d in 1:D) {
    out[d, ] <- log(rgamma(n, alpha[d], 1))
  }
  
  # Use colLogSumExps (natural log) and convert to log base 2
  log2_norm <- matrixStats::colLogSumExps(out) 

  res <- sweep(out, 2, log2_norm, FUN = "-")/ log(2)
  
  return(res)
}
