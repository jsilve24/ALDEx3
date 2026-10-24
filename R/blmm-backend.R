# Computational backends for BLMM. Both use the same cached cross-products;
# only the random-effect system's factorization and solves differ.

# Internal controls, deliberately independent of the inference selector.
blmm_backend_control <- function(backend = getOption("ALDEx3.blmm.backend", "auto"),
                                 memory = getOption("ALDEx3.blmm.memory", Inf)) {
  backend <- match.arg(backend, c("auto", "dense", "sparse"))
  if (!is.numeric(memory) || length(memory) != 1L || is.na(memory) || memory <= 0)
    stop("blmm: memory budget must be a positive number of bytes")
  list(backend = backend, memory = memory)
}

# A structural factor built with positive entries contains every possible fill
# location, including locations whose fitted covariance contribution is zero.
blmm_symbolic_factor <- function(A_pattern) {
  Matrix::Cholesky(A_pattern, LDL = FALSE, super = FALSE)
}

blmm_design <- function(reTrms, X, backend = getOption("ALDEx3.blmm.backend", "auto"),
                        draws = 1L, workers = 1L, features = 1L) {
  control <- blmm_backend_control(backend)
  Zt <- methods::as(reTrms$Zt, "dgCMatrix")
  T <- methods::as(reTrms$Lambdat, "dgCMatrix")
  G <- methods::as(Matrix::tcrossprod(Zt), "generalMatrix")
  symbolic <- NULL
  factor_nnz <- factor_work <- NA_real_
  if (control$backend != "dense") {
    Z_pattern <- Zt
    Z_pattern@x[] <- 1
    T_pattern <- T
    T_pattern@x[] <- 1
    A_pattern <- Matrix::forceSymmetric(
      Matrix::tcrossprod(T_pattern %*% Z_pattern) + Matrix::Diagonal(nrow(T)))
    # A lower bound rejects already-infeasible symbolic setup. Actual fill and
    # AD working estimates are checked below before allocating a dense backend.
    if (12 * (length(G@x) + length(A_pattern@x)) > control$memory)
      stop("blmm: sparse design exceeds the configured working-memory estimate")
    symbolic <- blmm_symbolic_factor(A_pattern)
    factor_nnz <- sum(symbolic@colcount)
    factor_work <- sum(as.double(symbolic@colcount)^2)
  }
  selection <- blmm_select_backend(
    control, nrow(X), ncol(X), nrow(T), length(reTrms$theta), draws, workers,
    factor_nnz, factor_work, features)
  structure(list(
    Zt = Zt, Lambdat = T, Lind = reTrms$Lind,
    G = G, ZX = as.matrix(Zt %*% X), XX = crossprod(X),
    symbolic = if (selection$backend == "sparse") symbolic else NULL,
    backend = selection$backend, selection = selection
  ), class = "blmm_design")
}

# Conservative operation-count proxy calibrated on five declared small/medium
# training designs; held-out designs are used only for validation. The margin
# avoids switching for predicted savings below a millisecond. This is not a
# universal timing guarantee; callers can force either computational backend.
blmm_select_backend <- function(control, n, p, q, n_theta, draws, workers,
                                factor_nnz, factor_work, features = 1L) {
  # Conservative working estimates include tape operations, response buffers,
  # results and PSOCK input copies. They are not an operating-system RSS cap.
  dense_work <- q^3 / 3 + 2 * q^2 * (p + draws)
  sparse_work <- factor_work + 2 * factor_nnz * (p + draws)
  input_copies <- if (.Platform$OS.type == "windows") workers + 1 else 1
  shared_bytes <- 8 * features * draws * (input_copies * n + 5 * p + n_theta + 1)
  response_bytes <- 8 * (4 * q * draws + n * draws + 4 * p * draws)
  dense_bytes <- shared_bytes + workers * (8 * 6 * q^2 + response_bytes + 64 * dense_work)
  sparse_bytes <- shared_bytes + workers * (48 * factor_nnz + response_bytes + 64 * sparse_work)
  predicted_saving <- 0.013 * dense_work / 1e6 -
    0.20 * sparse_work / 1e6 - 0.0014
  if (control$backend == "auto") {
    dense_ok <- dense_bytes <= control$memory
    sparse_ok <- sparse_bytes <= control$memory
    if (!dense_ok && !sparse_ok)
      stop("blmm: neither backend fits the configured working-memory estimate")
    selected <- if (sparse_ok && (!dense_ok || predicted_saving > 0.001)) "sparse" else "dense"
    reason <- if (!dense_ok || !sparse_ok) "working-memory feasibility" else
      "calibrated factor-work estimate with setup margin"
  } else {
    selected <- control$backend
    reason <- "internal override"
  }
  required <- if (selected == "dense") dense_bytes else sparse_bytes
  if (required > control$memory)
    stop("blmm: requested backend exceeds the configured working-memory estimate")
  list(backend = selected, requested = control$backend,
       reason = reason, selector_version = 1L,
       n = n, p = p, q = q, n_theta = n_theta, draws = draws, workers = workers,
       features = features, dense_work = dense_work, sparse_work = sparse_work,
       predicted_sparse_saving = predicted_saving, selection_margin = 0.001,
       factor_nnz = factor_nnz, factor_work = factor_work,
       dense_bytes = dense_bytes, sparse_bytes = sparse_bytes,
       memory_budget = control$memory)
}

blmm_response_cache <- function(design, X, Y) {
  Y <- as.matrix(Y)
  list(ZY = as.matrix(design$Zt %*% Y), XY = crossprod(X, Y),
       yy = colSums(Y^2))
}

blmm_factor_solve <- function(factor, rhs, backend) {
  rhs <- as.matrix(rhs)
  if (backend == "sparse") {
    as.matrix(Matrix::solve(factor, rhs, system = "A"))
  } else {
    backsolve(factor, forwardsolve(t(factor), rhs))
  }
}

# A strict validation run exposes backend failures instead of exact rescue.
# Ordinary runs retain the existing statistical rescue boundaries unchanged.
blmm_strict <- function() isTRUE(getOption("ALDEx3.blmm.strict", FALSE))
