make_p_adjust_fixture <- function() {
  coefficients <- array(c(
    -3.2, -.4,  2.6,  .2,
     3.1,  .5, -2.5, -.3,
    -1.1, -3.0,  .3,  2.7
  ), dim=c(2, 2, 3),
  dimnames=list(c("alpha", "beta"), c("taxon_a", "taxon_b"), NULL))
  # An NA must stay NA through both tail aggregation and p.adjust().
  coefficients["alpha", "taxon_b", 3] <- NA_real_
  contrasts <- array(c(
    -2.7, .7,
     2.8, -.8,
    -.2, 3.2
  ), dim=c(1, 2, 3),
  dimnames=list("beta", c("taxon_a", "taxon_b"), NULL))

  structure(list(
    estimate=coefficients,
    std.error=array(1, dim(coefficients), dimnames=dimnames(coefficients)),
    df=array(9, dim(coefficients), dimnames=dimnames(coefficients)),
    contrasts=list(
      estimate=contrasts,
      std.error=array(1, dim(contrasts), dimnames=dimnames(contrasts)),
      df=array(9, dim(contrasts), dimnames=dimnames(contrasts))
    )
  ), class="aldex")
}

# This deliberately does not call aldex.pvals(): it is the per-draw,
# per-direction, per-family definition against which aldex.p.adjust() is tested.
p_adjust_group_oracle <- function(object, families, method, onesided) {
  pieces <- list(object[c("estimate", "std.error", "df")])
  if (!is.null(object$contrasts)) {
    pieces[[2]] <- object$contrasts[c("estimate", "std.error", "df")]
  }
  estimate <- abind::abind(lapply(pieces, `[[`, "estimate"), along=1)
  std.error <- abind::abind(lapply(pieces, `[[`, "std.error"), along=1)
  df <- abind::abind(lapply(pieces, `[[`, "df"), along=1)
  lower <- array(stats::pt(estimate / std.error, df, lower.tail=TRUE),
                 dim(estimate), dimnames(estimate))
  upper <- array(stats::pt(estimate / std.error, df, lower.tail=FALSE),
                 dim(estimate), dimnames(estimate))
  if (identical(onesided, FALSE)) {
    lower <- array(pmin(1, 2 * lower), dim(lower), dimnames(lower))
    upper <- array(pmin(1, 2 * upper), dim(upper), dimnames(upper))
  }

  adjust <- function(values) {
    out <- array(NA_real_, dim(values), dimnames(values))
    for (s in seq_len(dim(values)[3])) for (family in families) {
      out[family, , s] <- stats::p.adjust(
        as.vector(values[family, , s, drop=FALSE]), method=method
      )
    }
    out
  }
  adjusted.lower <- adjust(lower)
  adjusted.upper <- adjust(upper)
  collect <- function(lower, upper) {
    low <- apply(lower, c(1, 2), mean)
    high <- apply(upper, c(1, 2), mean)
    if (identical(onesided, FALSE)) pmin(low, high)
    else if (identical(onesided, "lower")) low
    else high
  }
  list(raw=collect(lower, upper), adjusted=collect(adjusted.lower, adjusted.upper))
}

test_that("aldex.p.adjust adjusts combined named families before draw aggregation", {
  fit <- make_p_adjust_fixture()
  # Row 1 is coefficient alpha, row 2 coefficient beta, row 3 contrast beta.
  # beta is deliberately left unlisted and must form its own family.
  families <- list(c(1L, 3L), 2L)
  expected.groups <- list(c("coefficient:alpha", "contrast:beta"),
                          "coefficient:beta")

  for (method in c("BH", "holm")) for (side in list(FALSE, "lower", "upper")) {
    actual <- aldex.p.adjust(fit, groups=c("alpha", "contrast:beta"),
                             method=method, onesided=side)
    expected <- p_adjust_group_oracle(fit, families, method, side)
    expect_equal(unname(actual$p.val), unname(expected$raw[1:2, ]))
    expect_equal(unname(actual$contrasts$p.val), unname(expected$raw[3, , drop=FALSE]))
    expect_equal(unname(actual$p.val.adj), unname(expected$adjusted[1:2, ]))
    expect_equal(unname(actual$contrasts$p.val.adj),
                 unname(expected$adjusted[3, , drop=FALSE]))
    expect_identical(actual$p.adjust.groups, expected.groups)
    expect_identical(actual$p.adjust.method, method)
    expect_identical(actual$onesided, side)
  }

  # The ranks swap among draws, so adjusting final averaged values is wrong.
  expected <- p_adjust_group_oracle(fit, families, "BH", FALSE)
  adjusted.after_mean <- expected$raw
  for (family in families) {
    adjusted.after_mean[family, ] <- matrix(
      stats::p.adjust(as.vector(expected$raw[family, ]), "BH"),
      nrow=length(family)
    )
  }
  expect_false(isTRUE(all.equal(expected$adjusted, adjusted.after_mean)))
  expect_true(is.na(actual$p.val["alpha", "taxon_b"]))
  expect_true(is.na(actual$p.val.adj["alpha", "taxon_b"]))
})

test_that("aldex.p.adjust resolves families, defaults, and invalid references", {
  fit <- make_p_adjust_fixture()
  explicit <- aldex.p.adjust(fit,
                              groups=list("coefficient:alpha", "contrast:beta"),
                              method="holm", onesided="upper")
  expect_identical(explicit$p.adjust.groups,
                   list("coefficient:alpha", "contrast:beta", "coefficient:beta"))

  fallback <- aldex.p.adjust(fit)
  fallback.expected <- p_adjust_group_oracle(fit, list(1L, 2L, 3L), "BH", FALSE)
  expect_equal(unname(fallback$p.val.adj), unname(fallback.expected$adjusted[1:2, ]))
  expect_identical(fallback$p.adjust.method, "BH")
  expect_identical(fallback$onesided, FALSE)

  stored <- fit
  stored$p.adjust.method <- "holm"
  stored$onesided <- "upper"
  expected <- p_adjust_group_oracle(stored, list(1L, 2L, 3L), "holm", "upper")
  defaults <- aldex.p.adjust(stored)
  expect_equal(unname(defaults$p.val.adj), unname(expected$adjusted[1:2, ]))
  expect_equal(unname(defaults$contrasts$p.val.adj),
               unname(expected$adjusted[3, , drop=FALSE]))
  expect_identical(defaults$p.adjust.groups,
                   list("coefficient:alpha", "coefficient:beta", "contrast:beta"))

  expect_error(aldex.p.adjust(fit, groups="beta"), "ambiguous")
  expect_error(aldex.p.adjust(fit, groups="missing"), "unknown")
  expect_error(aldex.p.adjust(fit, groups=c("alpha", "alpha")), "groups")
  expect_error(aldex.p.adjust(fit, groups=character()), "group")
  expect_error(aldex.p.adjust(fit, groups=list(NA_character_)), "group")
  expect_error(aldex.p.adjust(fit, groups=list(character())), "group")
})

test_that("aldex.p.adjust validates stored inferential arrays and singleton dimensions", {
  fit <- make_p_adjust_fixture()
  missing.coefficient.df <- fit
  missing.coefficient.df$df <- NULL
  expect_error(aldex.p.adjust(missing.coefficient.df), "df")
  missing.contrast.se <- fit
  missing.contrast.se$contrasts$std.error <- NULL
  expect_error(aldex.p.adjust(missing.contrast.se), "std.error")

  single <- structure(list(
    estimate=array(-2, c(1, 1, 2), dimnames=list("x", "taxon", NULL)),
    std.error=array(1, c(1, 1, 2), dimnames=list("x", "taxon", NULL)),
    df=array(7, c(1, 1, 2), dimnames=list("x", "taxon", NULL))
  ), class="aldex")
  adjusted <- aldex.p.adjust(single, groups="x")
  expect_identical(dim(adjusted$p.val), c(1L, 1L))
  expect_equal(adjusted$p.val, adjusted$p.val.adj)
  expect_identical(adjusted$p.adjust.groups, list("coefficient:x"))

  extreme.upper <- structure(list(
    estimate=array(1e10, c(1, 1, 1), dimnames=list("x", "taxon", NULL)),
    std.error=array(1, c(1, 1, 1), dimnames=list("x", "taxon", NULL)),
    df=array(10, c(1, 1, 1), dimnames=list("x", "taxon", NULL))
  ), class="aldex")
  upper <- aldex.p.adjust(extreme.upper, onesided="upper")
  expect_gt(upper$p.val[1, 1], 0)
  expect_equal(upper$p.val[1, 1], stats::pt(1e10, 10, lower.tail=FALSE))
})

test_that("aldex.p.adjust reuses fitted and streamed draws without fitting or sampling", {
  counts <- rbind(
    c(20, 19, 22, 90, 95, 88),
    c(90, 94, 87, 21, 20, 24),
    c(50, 48, 51, 49, 52, 50)
  )
  rownames(counts) <- paste0("taxon_", 1:3)
  data <- data.frame(x=rep(c(-.5, .5), each=3))
  identities <- diag(2)
  dimnames(identities) <- list(c("(Intercept)", "x"), c("(Intercept)", "x"))
  common <- list(Y=counts, X=~x, data=data, method="lm", test="t",
                 nsample=3, scale=matrix(0, 6, 3), n.cores=1,
                 contrasts=identities,
                 return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  set.seed(902)
  fitted <- do.call(aldex, c(common, list(streamsize=Inf)))
  set.seed(903)
  streamed <- do.call(aldex, c(common, list(streamsize=.002)))
  expect_identical(fitted$p.adjust.method, "BH")
  expect_identical(fitted$onesided, FALSE)
  expect_true(streamed$streaming)
  expect_false(any(c("logComp", "logScale") %in% names(streamed)))

  testthat::local_mocked_bindings(
    aldex.sampler=function(...) stop("sampling must not occur"),
    fflm=function(...) stop("fitting must not occur"), .package="ALDEx3"
  )
  unchanged <- aldex.p.adjust(fitted)
  expect_equal(unchanged$p.val, fitted$p.val)
  expect_equal(unchanged$p.val.adj, fitted$p.val.adj)
  grouped <- aldex.p.adjust(fitted,
                             groups=c("coefficient:x", "contrast:x"))
  for (field in c("estimate", "std.error", "df", "X", "method", "test")) {
    expect_identical(grouped[[field]], fitted[[field]], info=field)
  }
  expect_identical(grouped$contrasts$estimate, fitted$contrasts$estimate)
  expect_identical(grouped$contrasts$std.error, fitted$contrasts$std.error)
  expect_identical(grouped$contrasts$df, fitted$contrasts$df)
  expect_equal(grouped$p.val, fitted$p.val)
  expect_equal(summary(grouped, component="contrasts")$p.val.adj,
               as.vector(grouped$contrasts$p.val.adj))

  reconstructed <- aldex.p.adjust(streamed, groups="coefficient:x")
  expect_identical(reconstructed$estimate, streamed$estimate)
  expect_identical(reconstructed$contrasts$estimate, streamed$contrasts$estimate)
  expect_equal(reconstructed$p.val, streamed$p.val)
  expect_equal(reconstructed$p.val.adj, streamed$p.val.adj)
})

test_that("aldex.p.adjust reconstructs stored mixed-model coefficient and contrast tails", {
  data <- data.frame(
    g=factor(rep(seq_len(8), each=4)),
    x=rep(c(-1, -.25, .25, 1), 8)
  )
  subject <- rep(seq(-.35, .35, length.out=8), each=4)
  counts <- rbind(
    round(70 + 15 * data$x + 11 * subject + seq_len(nrow(data)) %% 5),
    round(90 - 12 * data$x + 7 * subject + (2 * seq_len(nrow(data))) %% 7)
  )
  rownames(counts) <- c("taxon_a", "taxon_b")
  contrast <- matrix(c(-1, 1), 1,
                     dimnames=list("x_minus_intercept", c("(Intercept)", "x")))
  common <- list(Y=counts, data=data, n.cores=1, nsample=2,
                 scale=matrix(0, nrow(data), 2), contrasts=contrast,
                 return.pars=c("estimate", "std.error", "df", "p.val", "p.val.adj"))
  specifications <- list(
    list(method="lme4", X=~x + (1 | g)),
    list(method="nlme", X=~x, random=~1 | g),
    list(method="blmm", X=~x + (1 | g), test="t.satterthwaite"),
    list(method="blmm", X=~x + (1 | g), test="t.residual")
  )
  for (specification in specifications) {
    set.seed(940 + match(specification$method, c("lme4", "nlme", "blmm")))
    if (identical(specification$test, "t.residual")) {
      expect_warning(fit <- do.call(aldex, c(common, specification)),
                     "t[.]residual.*deprecated")
    } else {
      fit <- do.call(aldex, c(common, specification))
    }
    reconstructed <- aldex.p.adjust(fit)
    expect_equal(reconstructed$p.val, fit$p.val, tolerance=1e-12,
                 info=paste(specification$method, specification$test))
    expect_equal(reconstructed$p.val.adj, fit$p.val.adj, tolerance=1e-12,
                 info=paste(specification$method, specification$test))
    expect_equal(reconstructed$contrasts$p.val, fit$contrasts$p.val,
                 tolerance=1e-12, info=paste(specification$method, specification$test))
    expect_equal(reconstructed$contrasts$p.val.adj, fit$contrasts$p.val.adj,
                 tolerance=1e-12, info=paste(specification$method, specification$test))
  }
})
