aldex.sampler <- function(Y, X, nsample, scale=NULL, scale.args, return.pars) {
  N <- ncol(Y)
  D <- nrow(Y)

  
  ## dirichlet sample
  logComp <- array(NA, c(D, N, nsample))
  for (i in 1:N) {
    logComp[,i,] <- rLogDirichlet(nsample, Y[,i]+0.5)
  }

  ## sample from scale model
  if (is.null(scale)){
    stop("You probably want a scale model :)")
  } else if (is.matrix(scale)) {
    stopifnot(dim(scale)==c(N, nsample))
    logScale <- scale
  } else if (is.function(scale)) {
    req.args <- formalArgs(scale)
    scale.args <- scale.args[intersect(names(scale.args), req.args)]
    if ("logComp" %in% req.args) scale.args$logComp <- logComp
    if ("X" %in% req.args) scale.args$X <- X
    if ("Y" %in% req.args) scale.args$Y <- Y
    logScale <- do.call(scale, scale.args)
    ## some error checking to make sure whats returned has right dimensions
  } 

  out <- list()
  ## compute scaled abundances (W)
  out$logW <- sweep(logComp, c(2,3), logScale, FUN=`+`)

  if ("logScale" %in% return.pars) out$logScale <- logScale
  if ("logComp" %in% return.pars) out$logComp <- logComp

  return(out)
 }

##' Calculate p-values adjusting for changes in sign as described by Nixon
##' et al. (2024) in Beyond Normalization: Incorporating Scale Uncertainty
##' in Microbiome and Gene Expression Analysis (internal only)
##'
##' This internal function operates on arrays of directional probabilities.
##' To adjust an already fitted object using coefficient or contrast names,
##' use \code{\link{aldex.p.adjust}} instead.
##'
##' @param p.lower A P x D x S array of lower-tail probabilities, with P
##'   coefficient or contrast rows, D taxa/genes, and S Monte Carlo draws.
##' @param p.upper A P x D x S array of upper-tail probabilities, with the
##'   same row, feature, and draw order as \code{p.lower}.
##' @param p.adjust.method An adjutment method for p.adjust
##' @param onesided (default: FALSE) if sided return p-values for two-sided
##'   test. Otherwise if "lower" or "upper" return one-sided test corresponding
##'   to test that estimate is negative or positive respectively.
##' @param groups \code{NULL} (the default) adjusts each row separately across
##'   its D features. Otherwise, supply a list of integer vectors indexing the
##'   first dimension of the probability arrays. Every row index from 1 to P
##'   must appear exactly once across the list; missing or repeated indices
##'   are errors. For example, with three rows, \code{list(c(1L, 2L), 3L)}
##'   pools all \code{2 * D} hypotheses from rows 1 and 2 into one family and adjusts
##'   row 3's D hypotheses separately. Each family is adjusted within every
##'   draw and direction, before Monte Carlo aggregation. This internal
##'   argument accepts row indices, not names or \code{"contrast:"} prefixes.
##'   Unlike \code{aldex.p.adjust}, it does not add unlisted singleton rows.
##' @return A list with P x D matrices with the non-adjusted and adjusted
##'   p-values.
##' @seealso \code{\link{aldex.p.adjust}} for the public interface on fitted
##'   objects, which accepts names and automatically handles unlisted rows.
##' @examples
##' # Three hypothesis rows, two features, and two Monte Carlo draws.
##' lower <- array(c(0.01, 0.04, 0.20, 0.03, 0.30, 0.60,
##'                  0.02, 0.08, 0.25, 0.05, 0.35, 0.70), c(3, 2, 2))
##' upper <- 1 - lower
##' # Rows 1 and 2 share a family of four hypotheses; row 3 has two.
##' grouped <- ALDEx3:::aldex.pvals(lower, upper, "BH",
##'                                groups = list(c(1L, 2L), 3L))
##' grouped$p.adj.res
##' # All six hypotheses in one family, or the default separate families.
##' all_rows <- ALDEx3:::aldex.pvals(lower, upper, "BH", groups = list(1:3))
##' separate <- ALDEx3:::aldex.pvals(lower, upper, "BH")
##' @references Nixon G, Gloor GB, Silverman JD (2025).
##'   "Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization".
##'   Genome Biology. \doi{10.1186/s13059-025-03609-3}
##' @author Justin Silverman, Kyle McGovern
##' @keywords internal
aldex.pvals <- function(p.lower, p.upper, p.adjust.method, onesided=FALSE,
                        groups=NULL) {
  #### p-value calculations, accounting for sign changes ####
  if (onesided==FALSE) {
    p.lower.local <- array(pmin(1, 2 * p.lower), dim = dim(p.lower))
    p.upper.local <- array(pmin(1, 2 * p.upper), dim = dim(p.upper))
  } else {
    p.lower.local <- p.lower
    p.upper.local <- p.upper
  }

  if (is.null(groups)) {
    p.lower.adj <- apply(p.lower.local, c(1,3), function(item) {
      p.adjust(item, method=p.adjust.method)
    })
    p.upper.adj <- apply(p.upper.local, c(1,3), function(item) {
      p.adjust(item, method=p.adjust.method)
    })

    # apply() drops the feature axis when p.adjust() returns a single value.
    # Restore D x P x S so the subsequent means always average over draws.
    if (dim(p.lower)[2L] == 1L) {
      axes <- c(2L, 1L, 3L)
      p.lower.adj <- array(p.lower.adj, dim(p.lower)[axes], dimnames(p.lower)[axes])
      p.upper.adj <- array(p.upper.adj, dim(p.upper)[axes], dimnames(p.upper)[axes])
    }
  } else {
    groups <- .aldex_validate_p_adjust_groups(groups, dim(p.lower)[1L])
    p.lower.adj <- .aldex_p_adjust_groups(p.lower.local, groups,
                                           p.adjust.method)
    p.upper.adj <- .aldex_p_adjust_groups(p.upper.local, groups,
                                           p.adjust.method)
  }

  p.lower.mean <- apply(p.lower.local, c(2,1), mean)
  p.upper.mean <- apply(p.upper.local, c(2,1), mean)
  rm(p.lower.local, p.upper.local)

  p.res <- c()
  for(col_i in 1:ncol(p.lower.mean)) {
    tmp_mat <- cbind(p.lower.mean[,col_i],
                     p.upper.mean[,col_i])
    p.res <- rbind(p.res, apply(tmp_mat, 1, min))
  }
  ## rm(p.lower.mean, p.upper.mean)

  p.lower.mean.adj <- apply(p.lower.adj, c(1,2), mean)
  p.upper.mean.adj <- apply(p.upper.adj, c(1,2), mean)
  p.adj.res <- c()
  for(col_i in 1:ncol(p.lower.mean.adj)) {
    tmp_mat <- cbind(p.lower.mean.adj[,col_i],
                     p.upper.mean.adj[,col_i])
    p.adj.res <- rbind(p.adj.res, apply(tmp_mat, 1, min))
  }
  ## rm(p.lower.mean.adj, p.upper.mean.adj)
  if (onesided==FALSE) {
    return(list(p.res=p.res, p.adj.res=p.adj.res))
  } else if (onesided=="lower") {
    return(list(p.res=t(p.lower.mean), p.adj.res=t(p.lower.mean.adj)))
  } else if (onesided =="upper") {
    return(list(p.res=t(p.upper.mean), p.adj.res=t(p.upper.mean.adj)))
  } else {
    stop("not a recognized option for argument `onesided`")
  }
 }

.aldex_validate_p_adjust_groups <- function(groups, n.rows) {
  if (!is.list(groups) || !length(groups)) {
    stop("`groups` must be a nonempty list of row-index vectors")
  }

  groups <- lapply(groups, function(group) {
    if (!is.numeric(group) || !is.null(dim(group)) || !length(group) ||
        anyNA(group) || any(!is.finite(group)) ||
        any(group < 1 | group > .Machine$integer.max | group != floor(group))) {
      stop("each `groups` entry must be a nonempty vector of finite integer row indices")
    }
    as.integer(group)
  })
  indices <- unlist(groups, use.names=FALSE)
  if (length(indices) != n.rows || !identical(sort(indices), seq_len(n.rows))) {
    stop("`groups` must partition all p-value array rows exactly once")
  }
  groups
}

.aldex_p_adjust_groups <- function(p.values, groups, method) {
  axes <- c(2L, 1L, 3L)
  adjusted <- array(NA_real_, dim(p.values)[axes], dimnames(p.values)[axes])
  for (draw in seq_len(dim(p.values)[3L])) {
    for (group in groups) {
      values <- matrix(p.adjust(
        as.vector(p.values[group, , draw, drop=FALSE]), method=method
      ), nrow=length(group), ncol=dim(p.values)[2L])
      adjusted[, group, draw] <- t(values)
    }
  }
  adjusted
}
