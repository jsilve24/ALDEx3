test_that("p-value aggregation preserves singleton axes and BH ordering", {
  for (P in c(1L, 2L)) for (D in c(1L, 3L)) for (S in c(1L, 8L)) {
    lo <- array(rep(c(.01, .7, .13, .95, .04, .35, .8), length.out=P*D*S),
                c(P, D, S))
    up <- 1 - lo
    for (side in list(FALSE, "lower", "upper")) {
      expected.raw <- expected.adj <- matrix(NA_real_, P, D)
      # Explicit oracle: adjust across features within each coefficient/draw,
      # then average each direction, then select the requested direction.
      for (j in seq_len(P)) {
        lower <- matrix(lo[j, , ], D, S)
        upper <- matrix(up[j, , ], D, S)
        if (identical(side, FALSE)) {
          lower <- pmin(2 * lower, 1)
          upper <- pmin(2 * upper, 1)
        }
        adjusted.lower <- adjusted.upper <- matrix(NA_real_, D, S)
        for (s in seq_len(S)) {
          adjusted.lower[, s] <- p.adjust(lower[, s], "BH")
          adjusted.upper[, s] <- p.adjust(upper[, s], "BH")
        }
        select <- function(lower, upper) {
          if (identical(side, FALSE)) pmin(rowMeans(lower), rowMeans(upper))
          else if (side == "lower") rowMeans(lower) else rowMeans(upper)
        }
        expected.raw[j, ] <- select(lower, upper)
        expected.adj[j, ] <- select(adjusted.lower, adjusted.upper)
      }
      actual <- ALDEx3:::aldex.pvals(lo, up, "BH", side)
      expect_equal(actual$p.res, expected.raw)
      expect_equal(actual$p.adj.res, expected.adj)
      if (D == 1L) expect_equal(actual$p.adj.res, actual$p.res)
    }
  }
})

test_that("aldex returns named single-feature p-values with and without streaming", {
  counts <- matrix(rep(100L, 12), 1, dimnames=list("gene", paste0("s", 1:12)))
  dat <- data.frame(x=rep(c(-.5, .5), each=6))
  logscale <- matrix(3 + .4 * dat$x + rep(c(-.3, .1, .2), 4), 12, 8)
  logscale <- logscale + outer(seq_len(12), seq_len(8), function(i, s) .05 * sin(i+s))
  run <- function(streamsize) aldex(counts, ~x, data=dat, scale=logscale,
    nsample=8, method="lm", test="t", n.cores=1, streamsize=streamsize,
    return.pars=c("p.val", "p.val.adj"))
  full <- run(Inf)
  streamed <- run(.003) # Three draws per batch, with a final two-draw batch.
  expect_false(full$streaming)
  expect_true(streamed$streaming)
  expect_identical(dim(full$p.val.adj), c(2L, 1L))
  expect_identical(dimnames(full$p.val.adj), list(c("(Intercept)", "x"), "gene"))
  expect_equal(full$p.val.adj, full$p.val)
  expect_equal(streamed$p.val.adj, full$p.val.adj)
})
