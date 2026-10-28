# A zero Cholesky column is invariant to changing the sign of every entry
# in that column. Only remove complete, exactly zero columns: a zero diagonal
# with nonzero entries below it does not have this symmetry at the fitted point.
blmm_zero_columns <- function(theta, metadata) {
  dropped <- integer(0)
  for (block in metadata$blocks) {
    column <- (block$lower - 1L) %/% block$size + 1L
    for (j in seq_len(block$size)) {
      indices <- block$theta[column == j]
      if (all(theta[indices] == 0)) dropped <- c(dropped, indices)
    }
  }
  dropped <- sort(unique(dropped))
  list(drop = dropped, keep = setdiff(seq_along(theta), dropped))
}

# Verify that the omitted coordinates decouple from the retained raw
# information. No regularisation of the Hessian is permitted for inference.
blmm_partial_information <- function(theta, lower, metadata, H, gradient) {
  n <- length(theta)
  if (!n || length(lower) != n || any(!is.finite(theta)) || anyNA(lower) ||
      any(lower == Inf) || any(theta < lower) || length(gradient) != n ||
      any(!is.finite(gradient)) || !is.matrix(H) ||
      !identical(dim(H), c(n, n)) || any(!is.finite(H)))
    stop("blmm: invalid partial boundary derivatives or coordinates")
  if (max(abs(H - t(H))) > 64 * .Machine$double.eps * max(1, max(abs(H))))
    stop("blmm: partial boundary Hessian is not symmetric")
  H <- (H + t(H)) / 2
  indices <- blmm_zero_columns(theta, metadata)
  dropped <- indices$drop
  keep <- indices$keep
  if (!length(dropped) || !length(keep))
    stop("blmm: partial boundary requires zero and retained coordinates")
  active <- is.finite(lower) & theta <= lower +
    1e-10 * (1 + abs(replace(lower, !is.finite(lower), 0)))
  if (any(active[keep])) stop("blmm: unsupported active partial boundary")
  if (any(gradient[dropped] != 0) || any(H[dropped, keep, drop = FALSE] != 0))
    stop("blmm: partial boundary derivatives do not decouple")
  if (min(eigen(H[dropped, dropped, drop = FALSE], symmetric = TRUE,
                only.values = TRUE)$values) < 0)
    stop("blmm: partial boundary has negative dropped curvature")
  retained <- blmm_validate_H(H[keep, keep, drop = FALSE])
  list(indices = keep, H = retained, chol = chol(retained), dropped = dropped)
}
