## These are functions that provide additional calculations/summaries after aldex computation

##' Function to compute cohensd on the results provided by the aldex function
##'
##' WARNING: this function is experimental and requires users read the
##' documetation fully.
##' @title Cohen's D 
##' @param m the output of a call to `aldex`
##' @param var if aldex was called with X being a pre-computed model matrix,
##'   this var should be an integer corresponding to a binary covariate
##'   indicateing the desired effect size to calculate (an effect size between
##'   two groups indicated by the binary covariate). For example, if the third
##'   covariate in the model is an indicator denoting health (0) and disease (1)
##'   then set `var=3`. In contrast, if X was a formula (in which case the
##'   `data` argument should have been specified) then `var` can be set to
##'   the unquoted name of the binary condition variable (e.g.,
##'   `var=condition`).
##' @return A (D x nsample)-matrix of Cohen's D statistics for the variable of
##'   interest
##' @author Justin Silverman
cohensd <- function(m, var) {
  expr <- substitute(var)
  if (is.numeric(expr)) {
    ## don't need to do anything, var is already numeric
  } else { # assume unquoted string
    var <- deparse(expr)
    var <- which(rownames(m$X)==var)
  }
  diff.mean <- m$estimate[var,,]
   D <- dim(m$estimate)[2] # number of taxa
   S <- dim(m$estimate)[3] # number of samples
   x <- m$X[var,]
   x0idx <- which(x==0)
   x1idx <- which(x==1)
   n0 <- length(x0idx) # number in group 0
   n1 <- length(x1idx) # number in group 1
   ## we need estimate of scaled abundande
   if (!all(c("logWpara", "logWperp") %in% names(m))) {
     stop("m must contain logWpara and logWperp samples")
   } else {
     logW <- sweep(m$logWpara, c(2,3), m$logWperp, FUN=`+`)
   }
   cohensd <- matrix(NA, D, S)
   for (d in 1:D) { # for each taxa
     for (s in 1:S) { # for each sample
       var0 <- var(logW[d,x0idx,s])
       var1 <- var(logW[d,x1idx,s])
       cohensd[d,s] <- diff.mean[d,s]/sqrt(((n0-1)*var0+(n1-1)*var1)/(n0+n1-2))
     }
   }
  return(cohensd)
 }
