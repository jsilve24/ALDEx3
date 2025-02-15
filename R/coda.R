##' Closure operation
##'
##' @param X should be a D x N matrix, closes along the D dimension 
##' @return D x N matrix with columns summing to 1
##' @author Justin Silverman
miniclo <- function(X) {
   tot <- colSums(X)
   return(sweep(X, 2, tot, FUN="/"))
 }

