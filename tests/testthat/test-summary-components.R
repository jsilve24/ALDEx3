summary_component_fixture <- function(contrasts = TRUE) {
  values <- function(parameters, offset, probabilities) {
    list(
      estimate = array(seq_len(8) + offset, c(2, 2, 2),
                       list(parameters, c("taxon_a", "taxon_b"), NULL)),
      std.error = array((seq_len(8) + offset) / 10, c(2, 2, 2),
                        list(parameters, c("taxon_a", "taxon_b"), NULL)),
      p.val.adj = matrix(probabilities, 2, 2,
                         dimnames = list(parameters, c("taxon_a", "taxon_b")))
    )
  }
  fit <- values(c("(Intercept)", "slope"), 0, c(.11, .22, .33, .44))
  if (contrasts) {
    # A requested contrast may have the same label as a coefficient intercept.
    fit$contrasts <- values(c("(Intercept)", "requested"), 10,
                            c(.55, .66, .77, .88))
  }
  structure(fit, class = "aldex")
}

test_that("the default combines coefficients and requested contrasts without recalculation", {
  fit <- summary_component_fixture()
  before <- serialize(fit, NULL)
  expected <- data.frame(
    parameter = c("slope", "slope", "(Intercept)", "requested",
                  "(Intercept)", "requested"),
    entity = c("taxon_a", "taxon_b", "taxon_a", "taxon_a", "taxon_b", "taxon_b"),
    estimate = c(4, 6, 13, 14, 15, 16),
    std.error = c(.4, .6, 1.3, 1.4, 1.5, 1.6),
    p.val.adj = c(.22, .44, .55, .66, .77, .88),
    type = c(rep("coefficient", 2), rep("contrast", 4)))

  expect_equal(summary(fit), expected)
  expect_identical(summary(fit), summary(fit, component = "all"))
  expect_identical(serialize(fit, NULL), before)
})

test_that("explicit summary filters preserve their schemas and numerical results", {
  fit <- summary_component_fixture()
  combined <- summary(fit)
  columns <- c("parameter", "entity", "estimate", "std.error", "p.val.adj")
  for (component in c("coefficients", "contrasts")) {
    type <- if (component == "coefficients") "coefficient" else "contrast"
    expected <- combined[combined$type == type, columns]
    rownames(expected) <- NULL
    expect_identical(summary(fit, component = component), expected)
  }
  expect_error(summary(fit, component = "invalid"), "arg.*should be one of")
})

test_that("without requested contrasts the default preserves the coefficient-only schema", {
  fit <- summary_component_fixture(contrasts = FALSE)
  for (ignore in c(TRUE, FALSE)) {
    expected <- summary(fit, ignore.intercept = ignore, component = "coefficients")
    expect_identical(summary(fit, ignore.intercept = ignore), expected)
    expect_named(expected, c("parameter", "entity", "estimate", "std.error", "p.val.adj"))
  }
  expect_error(summary(fit, component = "contrasts"), "does not contain requested contrasts")
})

test_that("intercept filtering affects coefficients only, including empty coefficient tables", {
  fit <- summary_component_fixture()
  shown <- summary(fit, ignore.intercept = FALSE)
  hidden <- summary(fit)
  expect_equal(shown$parameter[shown$type == "coefficient"],
               rep(c("(Intercept)", "slope"), 2))
  expect_equal(hidden$parameter[hidden$type == "coefficient"], rep("slope", 2))
  expect_identical(summary(fit, component = "contrasts", ignore.intercept = FALSE),
                   summary(fit, component = "contrasts", ignore.intercept = TRUE))

  fit$estimate <- fit$estimate[1, , , drop = FALSE]
  fit$std.error <- fit$std.error[1, , , drop = FALSE]
  fit$p.val.adj <- fit$p.val.adj[1, , drop = FALSE]
  only_contrasts <- summary(fit)
  expect_identical(only_contrasts$type, rep("contrast", 4))
  expect_identical(only_contrasts$parameter, rep(c("(Intercept)", "requested"), 2))
  expect_equal(nrow(summary(fit, component = "coefficients")), 0L)
})
