test_that("BLMM selection uses factor work, overrides, and memory feasibility", {
  choose <- function(q, nnz, work, draws = 30L, backend = "auto", memory = Inf,
                     workers = 1L) {
    ALDEx3:::blmm_select_backend(
      ALDEx3:::blmm_backend_control(backend, memory), n = 5 * q, p = 2L,
      q = q, n_theta = 1L, draws = draws, workers = workers,
      factor_nnz = nnz, factor_work = work)
  }
  small <- choose(8, 8, 8)
  large <- choose(100, 100, 100)
  filled <- choose(100, 5050, sum(seq_len(100)^2))
  expect_identical(small$backend, "dense")
  expect_identical(large$backend, "sparse")
  expect_identical(filled$backend, "dense")
  expect_identical(choose(100, 100, 100, backend = "dense")$backend, "dense")
  expect_identical(choose(8, 8, 8, backend = "sparse")$backend, "sparse")
  budget <- mean(c(small$dense_bytes, small$sparse_bytes))
  constrained <- choose(8, 8, 8, memory = budget)
  expect_identical(constrained$backend, "sparse")
  expect_identical(constrained$reason, "working-memory feasibility")
  expect_error(choose(8, 8, 8, memory = budget, backend = "dense"), "memory")
  expect_error(choose(8, 8, 8, memory = 1), "neither backend")
  expect_gt(choose(100, 100, 100, workers = 2L)$sparse_bytes, large$sparse_bytes)
  expect_error(ALDEx3:::blmm_backend_control(memory = NA_real_), "positive")
})

test_that("BLMM records selection through its opt-in diagnostic callback", {
  z <- inference_fixture(S = 1L, noise = 0)
  selected <- NULL
  old <- options(ALDEx3.blmm.backend = "auto", ALDEx3.blmm.strict = TRUE,
                 ALDEx3.blmm.diagnostics = function(x) selected <<- x)
  on.exit(options(old), add = TRUE)
  a <- ALDEx3:::blmm(z$logW, z$formula, z$data)
  expect_identical(selected$requested, "auto")
  expect_identical(selected$backend, "dense")
  expect_identical(selected$selector_version, 1L)
  expect_true(is.finite(selected$factor_nnz))
  options(ALDEx3.blmm.backend = "dense", ALDEx3.blmm.diagnostics = NULL)
  b <- ALDEx3:::blmm(z$logW, z$formula, z$data)
  expect_equal(a, b, tolerance = 1e-8)
})

test_that("sparse computation preserves public aggregation and streaming chunks", {
  z <- inference_fixture(S = 4L, noise = .002)
  N <- nrow(z$Y)
  saved <- aperm(z$logW, c(2, 1, 3))
  cursor <- 0L
  testthat::local_mocked_bindings(aldex.sampler = function(Y, X, nsample, ...) {
    take <- cursor + seq_len(nsample)
    cursor <<- cursor + nsample
    list(logW = saved[, , take, drop = FALSE])
  }, .package = "ALDEx3")
  old <- options(ALDEx3.blmm.strict = TRUE, ALDEx3.blmm.backend = "dense")
  on.exit(options(old), add = TRUE)
  C <- matrix(c(1, -.5), 1, dimnames = list("contrast", c("(Intercept)", "x")))
  for (streamsize in c(8000, N * 8 * 2 / 100000)) {
    run <- function(backend) {
      options(ALDEx3.blmm.backend = backend)
      cursor <<- 0L
      ALDEx3::aldex(matrix(10, 1, N), z$formula, z$data,
                    method = "blmm", n.cores = 1L, nsample = 4L,
                    scale = matrix(0, N, 4), contrasts = C, streamsize = streamsize,
                    return.pars = c("estimate", "std.error", "df", "p.val", "p.val.adj"))
    }
    a <- run("dense")
    b <- run("sparse")
    expect_equal(b, a, tolerance = 1e-5)
    expect_identical(b$streaming, streamsize < 8000)
  }
})
