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
    # already numeric
  } else {
    var <- deparse(expr)
    var <- which(rownames(m$X) == var)
  }

  diff.mean <- m$estimate[var,,]  # D x S
  D <- dim(m$estimate)[2]
  S <- dim(m$estimate)[3]

  x <- m$X[var,]
  x0idx <- which(x == 0)
  x1idx <- which(x == 1)
  n0 <- length(x0idx)
  n1 <- length(x1idx)

  if (!all(c("logWpara", "logWperp") %in% names(m))) {
    stop("m must contain logWpara and logWperp samples")
  }

  logW <- sweep(m$logWpara, c(2,3), m$logWperp, FUN = `+`)  # D x N x S

  # subset over group 0 and group 1
  logW0 <- logW[, x0idx, , drop = FALSE]  # D x n0 x S
  logW1 <- logW[, x1idx, , drop = FALSE]  # D x n1 x S

  # We need to compute var0[d, s] = var(logW[d, x0idx, s])
  # We'll do this by reshaping into matrices D*S x n0 and D*S x n1

  reshape_and_compute_var <- function(logWg, n) {
    D <- dim(logWg)[1]
    n_samples <- dim(logWg)[2]
    S <- dim(logWg)[3]
    # D x n x S -> (D*S) x n
    mat <- matrix(aperm(logWg, c(1,3,2)), nrow = D*S, ncol = n)
    # Compute row variances (sample variance with Bessel's correction)
    row_vars <- matrixStats::rowVars(mat)
    matrix(row_vars, nrow = D, ncol = S)
  }

  var0 <- reshape_and_compute_var(logW0, n0)
  var1 <- reshape_and_compute_var(logW1, n1)

  # Pooled variance
  pooled_var <- ((n0 - 1) * var0 + (n1 - 1) * var1) / (n0 + n1 - 2)
  cohensd <- diff.mean / sqrt(pooled_var)

  return(cohensd)
}
