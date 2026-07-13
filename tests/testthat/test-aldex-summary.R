test_that("cohensd runs", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  X <- formula(~condition)
  data <- data.frame(condition=condition)
  nsample <- 50
  foo <- aldex(Y, X, data=data, nsample=nsample, scale=clr.sm)
  a <- cohensd(foo, condition)
  b <- cohensd(foo, 2)
  expect_equal(a, b)
  expect_equal(dim(a), c(nrow(Y), nsample))
  expect_true(!any(is.na(a)))

  rm(condition)
  expect_equal(cohensd(foo, condition), b)
})


test_that("cohensd correct", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  X <- formula(~condition)
  data <- data.frame(condition=condition)
  nsample <- 50
  set.seed(322217)
  foo <- aldex(Y, X, data, nsample=nsample, scale=clr.sm)
  a <- cohensd(foo, "condition")
  s <- 5
  d <- 2
  logW <- sweep(foo$logComp, c(2,3), foo$logScale, FUN=`+`)[d,,s]
  x0 <- logW[condition==0]
  x1 <- logW[condition==1]
  mean0 <- mean(x0)
  mean1 <- mean(x1)
  n0 <- sum(condition==0)
  n1 <- sum(condition==1)
  denom <- sqrt(((n0-1)*var(x0) + (n1-1)*var(x1))/(n0+n1-2))
  expect_equal(unname(a[d,s]), (mean1-mean0)/denom)
})

test_that("cohensd uses exact coefficient matching", {
  Y <- matrix(1:120, 10, 12)
  cond <- c(rep(0, 6), rep(1, 6))
  cond2 <- rep(c(0, 1), 6)
  data <- data.frame(cond=cond, cond2=cond2)
  foo <- aldex(Y, ~cond + cond2, data=data, nsample=20, scale=clr.sm)

  expect_equal(cohensd(foo, "cond"), cohensd(foo, 2))
  expect_equal(cohensd(foo, "cond2"), cohensd(foo, 3))

  bad <- foo
  rownames(bad$X) <- c("(Intercept)", "condA", "condB")
  bad$X[3,] <- bad$X[2,]
  expect_error(cohensd(bad, "cond"), "exactly one binary coefficient")
})

test_that("aldex.effect reports documented effect diagnostics", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  data <- data.frame(condition=condition)
  foo <- aldex(Y, ~condition, data=data, nsample=20, scale=clr.sm)

  eff <- aldex.effect(foo, "condition")
  expect_equal(colnames(eff), c("parameter", "entity", "estimate",
                                "pooled.SD", "cohens.d", "overlap"))
  expect_equal(nrow(eff), nrow(Y))
  expect_true(all(eff$overlap >= 0 & eff$overlap <= 0.5))
  expect_equal(eff$cohens.d, unname(rowMeans(cohensd(foo, "condition"))))
})

test_that("summary.aldex works", {
  Y <- matrix(1:110, 10, 11)
  data <- data.frame(disease=c(rep(0, 5), rep(1, 6)))
  nsample <- 50
  foo <- aldex(Y, ~disease, data=data, nsample=nsample, scale=clr.sm)
  class(foo) <- "aldex"
  bar <- summary(foo)
  expect_equal(colnames(bar), c("parameter", "entity",  "estimate",
                                "std.error", "p.val.adj"))
})

test_that("aldex.plot handles all plot types and edge cases", {
  Y <- matrix(1:110, 10, 11)
  condition <- c(rep(0, 5), rep(1, 6))
  data <- data.frame(condition=condition)
  foo <- aldex(Y, ~condition, data=data, nsample=20, scale=clr.sm)

  grDevices::pdf(file=tempfile(fileext=".pdf"))
  on.exit(grDevices::dev.off(), add=TRUE)

  expect_silent(aldex.plot(foo, contrast="condition", plot="volcano"))
  expect_silent(aldex.plot(foo, contrast="condition", plot="effect"))
  expect_silent(aldex.plot(foo, contrast="condition", plot="MA"))
  expect_silent(aldex.plot(foo, contrast="condition", plot="water",
                           threshold=0, min.diff=Inf))

  streamed <- foo
  streamed$logComp <- NULL
  streamed$logScale <- NULL
  expect_error(aldex.plot(streamed, contrast="condition", plot="MA"),
               "logComp/logScale")
})
