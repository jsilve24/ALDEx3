## Contrast validation shared by ALDEx3 model engines.

#' Validate fixed-effect contrasts
#'
#' @param contrasts numeric matrix with one contrast per row and named columns.
#' @param X fixed-effect design matrix with observations in rows and named
#'   coefficient columns.
#' @return A numeric contrast matrix whose columns have been reordered to match
#'   `X`.
#' @noRd
.aldex_validate_contrasts <- function(contrasts, X) {
  if (is.null(contrasts)) return(NULL)
  if (!is.matrix(X) || !is.numeric(X) || !all(is.finite(X))) {
    stop("X must be a finite numeric design matrix when validating contrasts")
  }
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(colnames(X) == "") ||
      anyDuplicated(colnames(X))) {
    stop("X must have unique, nonempty coefficient names for contrasts")
  }
  if (!is.matrix(contrasts) || !is.numeric(contrasts)) {
    stop("contrasts must be a numeric matrix")
  }
  if (nrow(contrasts) < 1L || ncol(contrasts) < 1L) {
    stop("contrasts must have at least one row and one column")
  }
  if (!all(is.finite(contrasts))) stop("contrasts must contain only finite values")
  if (is.null(rownames(contrasts)) || anyNA(rownames(contrasts)) ||
      any(rownames(contrasts) == "") || anyDuplicated(rownames(contrasts))) {
    stop("contrasts must have unique, nonempty row names")
  }
  if (is.null(colnames(contrasts)) || anyNA(colnames(contrasts)) ||
      any(colnames(contrasts) == "") || anyDuplicated(colnames(contrasts))) {
    stop("contrasts must have unique, nonempty column names")
  }
  if (!setequal(colnames(contrasts), colnames(X)) ||
      ncol(contrasts) != ncol(X)) {
    stop("contrast columns must name the complete coefficient set in X")
  }
  contrasts <- contrasts[, colnames(X), drop = FALSE]
  if (any(rowSums(abs(contrasts)) == 0)) {
    stop("contrast rows must be nonzero")
  }

  ## A contrast c is estimable precisely when it lies in row(X).  Checking this
  ## before the full-rank engine guard gives a useful error for aliased terms.
  q <- qr(t(X))
  if (q$rank < ncol(X)) {
    Q <- qr.Q(q, complete = FALSE)[, seq_len(q$rank), drop = FALSE]
    projected <- Q %*% (crossprod(Q, t(contrasts)))
    residual <- t(contrasts) - projected
    tolerance <- 100 * .Machine$double.eps *
      pmax(1, sqrt(colSums(t(contrasts)^2)))
    bad <- which(sqrt(colSums(residual^2)) > tolerance)
    if (length(bad)) {
      stop("nonestimable contrast(s) for rank-deficient X: ",
           paste(rownames(contrasts)[bad], collapse = ", "))
    }
  }
  contrasts
}
