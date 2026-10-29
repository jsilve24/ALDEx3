# Regression tests for the explicit BLMM dense and sparse computational paths.
# These tests compare the observable profiled objective and conditional solves;
# factor permutations are backend details and are intentionally not compared.

blmm_sparse_ad_fixture <- function(seed = 1201L, S = 3L) {
  set.seed(seed)
  groups <- factor(rep(seq_len(10L), each = 5L))
  x <- rep(seq(-1, 1, length.out = 5L), 10L)
  data <- data.frame(g = groups, x = x)
  a <- rnorm(nlevels(groups), sd = 0.7)
  b <- rnorm(nlevels(groups), sd = 0.25)
  mean_y <- 0.5 + 0.8 * x + a[groups] + b[groups] * x
  Y <- vapply(seq_len(S), function(s) {
    mean_y + rnorm(length(mean_y), sd = 0.03 + 0.005 * s)
  }, numeric(length(mean_y)))
  list(
    formula = ~ x + (1 + x | g), data = data,
    Y = Y, logW = array(Y, c(nrow(Y), 1L, ncol(Y)))
  )
}

blmm_sparse_structure_fixture <- function(design = c("nested", "crossed"),
                                          seed = 1202L, S = 2L) {
  design <- match.arg(design)
  set.seed(seed)
  if (design == "nested") {
    g <- factor(rep(seq_len(12L), each = 4L))
    block <- factor(rep(rep(seq_len(6L), each = 2L), each = 4L))
    x <- rep(seq(-1, 1, length.out = 4L), 12L)
    data <- data.frame(g = g, block = block, x = x)
    form <- ~ x + (1 | block) + (1 | g)
    mean_y <- 0.4 + 0.7 * x + rnorm(nlevels(block), sd = 0.45)[block] +
      rnorm(nlevels(g), sd = 0.35)[g]
  } else {
    g <- factor(rep(seq_len(6L), each = 6L))
    h <- factor(rep(seq_len(6L), times = 6L))
    x <- rep(seq(-1, 1, length.out = 6L), 6L)
    data <- data.frame(g = g, h = h, x = x)
    form <- ~ x + (1 | g) + (1 | h)
    mean_y <- 0.4 + 0.7 * x + rnorm(nlevels(g), sd = 0.45)[g] +
      rnorm(nlevels(h), sd = 0.35)[h]
  }
  Y <- vapply(seq_len(S), function(s) {
    mean_y + rnorm(length(mean_y), sd = 0.025 + 0.005 * s)
  }, numeric(length(mean_y)))
  list(formula = form, data = data, Y = Y,
       logW = array(Y, c(nrow(Y), 1L, ncol(Y))))
}

blmm_sparse_contrast <- function(formula, data) {
  X <- model.matrix(reformulas::nobars(formula), data)
  C <- matrix(0, nrow = 2L, ncol = ncol(X),
              dimnames = list(c("x", "intercept_plus_x"), colnames(X)))
  C["x", "x"] <- 1
  C["intercept_plus_x", c("(Intercept)", "x")] <- c(1, 0.5)
  C
}

blmm_sparse_ad_objects <- function(z, backend) {
  p <- ALDEx3:::blmm_parse_formula(
    z$formula, z$data, backend = backend, draws = ncol(z$Y)
  )
  obj <- ALDEx3:::blmm_make_adfun(
    p$X, z$Y, p$basis, p$is_log, p$phi_init
  )
  list(parsed = p, obj = obj)
}

test_that("BLMM dense and sparse AD paths agree at fixed parameters", {
  z <- blmm_sparse_ad_fixture()
  dense <- blmm_sparse_ad_objects(z, "dense")
  sparse <- blmm_sparse_ad_objects(z, "sparse")
  expect_identical(dense$parsed$X, sparse$parsed$X)
  expect_identical(dense$parsed$phi_init, sparse$parsed$phi_init)
  expect_identical(dense$parsed$is_log, sparse$parsed$is_log)
  expect_true(inherits(sparse$parsed$basis, "blmm_design"))
  expect_identical(sparse$parsed$basis$backend, "sparse")

  correlation <- which(!is.finite(dense$parsed$lower))
  expect_length(correlation, 1L)
  for (value in c(-0.2, 0, 0.2)) {
    phi <- dense$parsed$phi_init
    phi[correlation] <- value
    fn_dense <- dense$obj$fn(phi)
    fn_sparse <- sparse$obj$fn(phi)
    expect_equal(fn_sparse, fn_dense, tolerance = 1e-8,
                 info = as.character(value))
    expect_equal(sparse$obj$gr(phi), dense$obj$gr(phi),
                 tolerance = 1e-6, info = paste("gradient", value))
    expect_equal(sparse$obj$he(phi), dense$obj$he(phi),
                 tolerance = 1e-6, info = paste("Hessian", value))
    rd <- dense$obj$report()
    rs <- sparse$obj$report()
    expect_equal(rs$PWRSS_s, rd$PWRSS_s, tolerance = 1e-8,
                 info = paste("PWRSS", value))
    expect_equal(rs$LX, rd$LX, tolerance = 1e-8,
                 info = paste("LX", value))
    expect_equal(rs$ldL2, rd$ldL2, tolerance = 1e-8,
                 info = paste("logdet", value))
  }
})

test_that("BLMM AD derivatives agree with central differences without numDeriv",
          {
  z <- blmm_sparse_ad_fixture(S = 2L)
  p <- ALDEx3:::blmm_parse_formula(z$formula, z$data, backend = "sparse",
                                    draws = ncol(z$Y))
  obj <- ALDEx3:::blmm_make_adfun(p$X, z$Y, p$basis, p$is_log, p$phi_init)
  phi <- p$phi_init
  h <- 2e-4
  gradient_fd <- numeric(length(phi))
  hessian_fd <- matrix(NA_real_, length(phi), length(phi))
  for (j in seq_along(phi)) {
    plus <- phi
    minus <- phi
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    gradient_fd[j] <- (obj$fn(plus) - obj$fn(minus)) / (2 * h)
    hessian_fd[, j] <- (obj$gr(plus) - obj$gr(minus)) / (2 * h)
  }
  expect_equal(as.numeric(obj$gr(phi)), gradient_fd, tolerance = 1e-4)
  expect_equal(obj$he(phi), hessian_fd, tolerance = 1e-4)

  obj$fn(phi)
  report <- obj$report()
  expected <- blmm_anchor_objective_from_report(
    report, nrow(p$X), ncol(p$X)
  )
  expect_equal(obj$fn(phi), expected, tolerance = 1e-8)
})

test_that("anchor sparse and dense paths agree for slope, nested, and crossed designs",
          {
  old <- options(ALDEx3.blmm.backend = "dense",
                 ALDEx3.blmm.strict = TRUE)
  on.exit(options(old), add = TRUE)
  for (design in c("slope", "nested", "crossed")) {
    z <- if (design == "slope") blmm_sparse_ad_fixture(S = 2L) else
      blmm_sparse_structure_fixture(design, S = 2L)
    C <- blmm_sparse_contrast(z$formula, z$data)
    fits <- vector("list", 2L)
    names(fits) <- c("dense", "sparse")
    for (backend in names(fits)) {
      options(ALDEx3.blmm.backend = backend)
      fits[[backend]] <- ALDEx3:::blmm(
        z$logW, z$formula, z$data, n.cores = 1L,
        inference = "anchor", contrasts = C
      )
    }
    for (field in names(fits$dense))
      expect_equal(fits$sparse[[field]], fits$dense[[field]], tolerance = 1e-5,
                   info = paste(design, field))
  }
})

# Only the anchor path performs R-level sparse factor solves.
test_that("strict sparse backend failures are exposed instead of rescued", {
  z <- blmm_sparse_ad_fixture(S = 2L)
  old <- options(ALDEx3.blmm.backend = "sparse", ALDEx3.blmm.strict = TRUE)
  on.exit(options(old), add = TRUE)
  testthat::local_mocked_bindings(
    blmm_factor_solve = function(...) stop("injected sparse backend failure"),
    .package = "ALDEx3"
  )
  expect_error(
    ALDEx3:::blmm(z$logW, z$formula, z$data, n.cores = 1L, inference = "anchor"),
    "injected sparse backend failure"
  )
})

test_that("dense design preparation does not require sparse symbolic analysis", {
  z <- blmm_sparse_ad_fixture(S = 2L)
  old <- options(ALDEx3.blmm.backend = "dense")
  on.exit(options(old), add = TRUE)
  testthat::local_mocked_bindings(
    blmm_symbolic_factor = function(...) stop("dense path called sparse symbolic factor"),
    .package = "ALDEx3"
  )
  parsed <- ALDEx3:::blmm_parse_formula(
    z$formula, z$data, backend = "dense", draws = ncol(z$Y)
  )
  expect_identical(parsed$basis$backend, "dense")
  expect_null(parsed$basis$symbolic)
})

test_that("strict forked backend failures rethrow the original condition", {
  skip_if(.Platform$OS.type == "windows", "POSIX fork behavior is required")
  z <- inference_fixture(seed = 731, groups = 12, design = "between",
                         noise = .01, S = 2L, effect = .6)
  logW <- abind::abind(z$logW, z$logW, along = 2L)
  old <- options(ALDEx3.blmm.backend = "sparse", ALDEx3.blmm.strict = TRUE)
  on.exit(options(old), add = TRUE)
  testthat::local_mocked_bindings(
    blmm_make_adfun = function(...) stop("injected fork backend failure"),
    .package = "ALDEx3"
  )

  warns <- character()
  err <- tryCatch(
    withCallingHandlers(
      ALDEx3:::blmm(logW, z$formula, z$data, n.cores = 2L),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "injected fork backend failure")
  if (length(warns) > 0L)
    expect_true(any(grepl("fork|core|worker|scheduled", warns, ignore.case = TRUE)))
})

test_that("both backends preserve whole-feature rescue identities", {
  whole_z <- inference_fixture(seed = 731, groups = 12, design = "between",
                               noise = .01, S = 2, effect = .6)
  old <- options(ALDEx3.blmm.backend = "dense", ALDEx3.blmm.strict = FALSE)
  on.exit(options(old), add = TRUE)

  whole <- whole_z$logW
  attr(whole, "blmm.force_fail_features") <- 1L
  exact <- ALDEx3:::blmm_exact_feature(whole_z$Y, whole_z$formula, whole_z$data)
  results <- vector("list", 2L)
  names(results) <- c("dense", "sparse")
  for (backend in names(results)) {
    options(ALDEx3.blmm.backend = backend)
    warns <- character()
    results[[backend]] <- withCallingHandlers(
      ALDEx3:::blmm(whole, whole_z$formula, whole_z$data, n.cores = 1L),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_length(warns, 1L)
    expect_match(warns[[1L]], "refit with lmerTest")
    for (field in names(exact))
      expect_equal(results[[backend]][[field]][, 1, , drop = FALSE],
                   exact[[field]], tolerance = 1e-10,
                   info = paste("whole", backend, field))
  }
})
