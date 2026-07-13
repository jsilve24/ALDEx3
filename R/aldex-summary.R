## These are functions that provide additional calculations/summaries after aldex computation

##' Function to compute cohensd on the results provided by the aldex function
##'
##' WARNING: this function is experimental and requires users read the
##' documentation fully.
##' @title Cohen's D 
##' @param m the output of a call to `aldex`
##' @param var if aldex was called with X being a pre-computed model matrix,
##'   this var should be an integer corresponding to a binary covariate
##'   indicating the desired effect size to calculate (an effect size between
##'   two groups indicated by the binary covariate). For example, if the third
##'   covariate in the model is an indicator denoting health (0) and disease (1)
##'   then set `var=3`. In contrast, if X was a formula (in which case the
##'   `data` argument should have been specified) then `var` can be set to
##'   the unquoted or quoted name of the binary condition variable (e.g.,
##'   `var=condition` or `var="condition"`).
##' @import matrixStats stats
##' @return A (D x nsample)-matrix of Cohen's d statistics for the variable of
##'   interest.
##' @author Justin Silverman
cohensd <- function(m, var) {
  idx <- .aldex_resolve_contrast(m, substitute(var))
  pooled <- .aldex_pooled_var(m, idx)
  diff.mean <- m$estimate[idx,,]  # D x S
  cohens.d <- diff.mean / sqrt(pooled$pooled.var)

  return(cohens.d)
}

##' ALDEx2-inspired effect diagnostics for an ALDEx3 binary contrast
##'
##' `aldex.effect` is a helper for effect-size diagnostics that Greg Gloor and
##' ALDEx2 users commonly inspect when evaluating pairwise effects. It is not a
##' replacement for `summary.aldex`, which reports model estimates, standard
##' errors and adjusted p-values. This function summarizes a single binary
##' ALDEx3 model contrast using Cohen's-d-based Monte Carlo diagnostics.
##'
##' For a selected contrast row `k`, feature `d`, Monte Carlo draw `s`, and
##' binary group indicators `x = object$X[k, ]`, let `B[d, s]` be the fitted
##' ALDEx3 coefficient stored in `object$estimate[k, d, s]`. Let
##' `logW[d, i, s] = object$logComp[d, i, s] + object$logScale[i, s]` be the
##' reconstructed log abundance for sample `i`. `aldex.effect` requires
##' `logComp` and `logScale`, so it is unavailable when `aldex` streams large
##' jobs and omits those arrays. For each feature and Monte Carlo draw, the
##' pooled within-group variance is
##' `((n0 - 1) * var(logW[d, x == 0, s]) + (n1 - 1) *
##' var(logW[d, x == 1, s])) / (n0 + n1 - 2)`, where `n0` and `n1` are the
##' group sizes.
##'
##' The returned columns are calculated as follows:
##' * `estimate`: the mean of `B[d, s]` over Monte Carlo draws.
##' * `pooled.SD`: the mean of the pooled within-group standard deviation,
##'   `sqrt(pooled variance)`, over Monte Carlo draws.
##' * `cohens.d`: the mean of `B[d, s] / pooled.SD[d, s]` over Monte Carlo
##'   draws. This is Cohen's-d-based and is not the same as the historical
##'   robust ALDEx2 effect statistic based on between/within dispersion.
##' * `overlap`: an ALDEx2-style directional uncertainty diagnostic. It counts
##'   the Monte Carlo Cohen's d draws below and above zero, smooths those two
##'   counts with a 0.5 pseudocount, converts them with the ALDEx/Aitchison
##'   mean, and returns the smaller smoothed sign probability. Values near 0
##'   indicate that nearly all Monte Carlo effect draws have the same sign.
##'   Values near 0.5 indicate uncertainty about direction. `overlap` is not a
##'   p-value, a confidence interval probability, or a literal area of overlap
##'   between two density curves.
##'
##' @title ALDEx3 Effect Diagnostics
##' @param object An object returned by `aldex`.
##' @param contrast The exact model coefficient or unambiguous binary data
##'   column to summarize.
##' @return A data.frame with columns `parameter`, `entity`, `estimate`,
##'   `pooled.SD`, `cohens.d`, and `overlap`.
##' @author Greg Gloor, Justin Silverman
##' @export
aldex.effect <- function(object, contrast) {
  idx <- .aldex_resolve_contrast(object, substitute(contrast), contrast)
  pooled <- .aldex_pooled_var(object, idx)
  diff.mean <- object$estimate[idx,,]  # D x S
  cohens.d <- diff.mean / sqrt(pooled$pooled.var)

  res <- data.frame(
    parameter = rep(rownames(object$X)[idx], dim(object$estimate)[2]),
    entity = colnames(object$estimate),
    estimate = rowMeans(diff.mean),
    pooled.SD = rowMeans(sqrt(pooled$pooled.var)),
    cohens.d = rowMeans(cohens.d),
    overlap = .aldex_effect_overlap(cohens.d),
    row.names = NULL
  )
  return(res)
}

## Shared exact contrast resolver for effect-size helpers and plotting.
.aldex_resolve_contrast <- function(object, expr, value) {
  req(object, c("X"))
  rn <- rownames(object$X)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(object$X)))

  resolve_name <- function(nm) {
    if (length(nm) != 1 || is.na(nm) || nm == "") {
      stop("contrast must identify exactly one binary coefficient")
    }
    exact <- which(rn == nm)
    if (length(exact) == 1) {
      return(exact)
    } else if (length(exact) > 1) {
      stop("contrast matches more than one row of object$X: ", nm)
    } else if (!is.null(object$data) && nm %in% colnames(object$data)) {
      candidates <- which(startsWith(rn, nm))
      if (length(candidates) != 1) {
        stop("contrast must map to exactly one binary coefficient: ", nm)
      }
      return(candidates)
    }
    return(NULL)
  }

  idx <- NULL
  if (is.numeric(expr)) {
    if (length(expr) != 1 || is.na(expr) || expr < 1 || expr > nrow(object$X)) {
      stop("contrast index must identify exactly one row of object$X")
    }
    idx <- as.integer(expr)
  } else {
    if (is.character(expr)) {
      nm <- expr
    } else {
      nm <- deparse(expr)
    }
    idx <- resolve_name(nm)
    if (is.null(idx) && !missing(value)) {
      val <- tryCatch(value, error=function(e) NULL)
      if (is.numeric(val) && length(val) == 1) {
        if (is.na(val) || val < 1 || val > nrow(object$X)) {
          stop("contrast index must identify exactly one row of object$X")
        }
        idx <- as.integer(val)
      } else if (is.character(val)) {
        idx <- resolve_name(val)
      }
    }
    if (is.null(idx)) {
      stop("contrast not found in object$X: ", nm)
    }
  }

  x <- object$X[idx,]
  if (!all(x %in% c(0, 1))) {
    stop("contrast must identify a binary 0/1 coefficient")
  }
  if (sum(x == 0) < 2 || sum(x == 1) < 2) {
    stop("contrast must have at least two samples in each group")
  }
  return(idx)
}

.aldex_pooled_var <- function(m, idx) {
  if (!all(c("logComp", "logScale") %in% names(m))) {
    stop("m must contain logComp and logScale samples\ntry reducing the nsample parameter")
  }

  x <- m$X[idx,]
  x0idx <- which(x == 0)
  x1idx <- which(x == 1)
  n0 <- length(x0idx)
  n1 <- length(x1idx)
  logW <- sweep(m$logComp, c(2,3), m$logScale, FUN = `+`)  # D x N x S

  logW0 <- logW[, x0idx, , drop = FALSE]  # D x n0 x S
  logW1 <- logW[, x1idx, , drop = FALSE]  # D x n1 x S

  reshape_and_compute_var <- function(logWg, n) {
    D <- dim(logWg)[1]
    S <- dim(logWg)[3]
    # D x n x S -> (D*S) x n
    mat <- matrix(aperm(logWg, c(1,3,2)), nrow = D*S, ncol = n)
    row_vars <- matrixStats::rowVars(mat)
    matrix(row_vars, nrow = D, ncol = S)
  }

  var0 <- reshape_and_compute_var(logW0, n0)
  var1 <- reshape_and_compute_var(logW1, n1)
  pooled_var <- ((n0 - 1) * var0 + (n1 - 1) * var1) / (n0 + n1 - 2)
  return(list(pooled.var = pooled_var, x0idx = x0idx, x1idx = x1idx))
}

.aldex_effect_overlap <- function(effect.draws) {
  apply(effect.draws, 1, function(row) {
    if (all(is.na(row))) warning("NAs in effect")
    row[is.na(row)] <- 0
    min(aitchison.mean(c(sum(row < 0), sum(row > 0)) + 0.5))
  })
}

# Private function from ALDEx2, via Andrew Fernandes.
aitchison.mean <- function( n, log=FALSE ) {

    # Input is a vector of non-negative integer counts.
    # Output is a probability vector of expected frequencies.
    # If log-frequencies are requested, the uninformative subspace is removed.

    n <- round( as.vector( n, mode="numeric" ) )
    if ( any( n < 0 ) ) stop("counts cannot be negative")

    a <- n + 0.5
    sa <- sum(a)

    log.p <- digamma(a) - digamma(sa)
    log.p <- log.p - mean(log.p)

    if ( log ) return(log.p)

    p <- exp( log.p - max(log.p) )
    p <- p / sum(p)
    return(p)
}

## summary <- function(object, ...) {
##   UseMethod("summary", object)
## }

##' Summarize an ALDEx3 result object
##'
##' Provides a summary of the adjusted p-values, estimates, and standard errors
##' from an ALDEx3 result object.
##'
##' This method extracts adjusted p-values from `object$p.val.adj`, along with
##' posterior estimates and standard errors averaged across Monte Carlo samples.
##' The result is returned as a long-format data.frame suitable for downstream
##' analysis or visualization.
##' 
##' @title Summary Method for ALDEx3 Objects
##' @param object An object of class \code{aldex}
##' @param ignore.intercept (default=TRUE), ignore intercept when creating
##'   summary table
##' @param ... Additional arguments (currently ignored).
##' @return A \code{data.frame} with columns \code{parameter}, \code{entity},
##'   \code{p.val.adjusted}, \code{estimate}, and \code{std.error}.
##' @export
##' @author Justin Silverman
summary.aldex <- function(object, ignore.intercept=TRUE, ...) {
  req(object, c("p.val.adj", "estimate", "std.error"))
  res <- array_to_df(object$p.val.adj)
  colnames(res) <- c("parameter", "entity", "p.val.adj")
  res$estimate <- array_to_df(rowMeans(object$estimate, dims=2))$value
  res$std.error <- array_to_df(rowMeans(object$std.error, dims=2))$value
  if (ignore.intercept) {
    res <- res[res$parameter != "(Intercept)",]
  }
  rownames(res) <- NULL
  return(res[,c("parameter", "entity", "estimate", "std.error", "p.val.adj")])
}

array_to_df <- function(arr) {
  dimn <- dimnames(arr)
  idx_grid <- expand.grid(dimn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx_grid$value <- as.vector(arr)
  return(idx_grid)
}
