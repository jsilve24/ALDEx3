plot_contrast_collision_fixture <- function() {
  parameters <- c("(Intercept)", "condition")
  entities <- c("taxon_a", "taxon_b")
  estimate <- array(c(0, 1, 0, -2), c(2, 2, 1),
                    list(parameters, entities, "draw"))
  std.error <- array(1, dim(estimate), dimnames=dimnames(estimate))
  p.val.adj <- matrix(c(.9, .01, .9, .2), 2, 2,
                      dimnames=list(parameters, entities))

  contrast.estimate <- array(c(100, -100), c(1, 2, 1),
                             list("condition", entities, "draw"))
  contrast.std.error <- array(1, dim(contrast.estimate),
                              dimnames=dimnames(contrast.estimate))
  contrast.p.val.adj <- matrix(c(.001, .001), 1, 2,
                               dimnames=list("condition", entities))

  structure(list(
    X = rbind(`(Intercept)`=rep(1, 4), condition=c(0, 0, 1, 1)),
    data = data.frame(condition=c(0, 0, 1, 1)),
    estimate = estimate,
    std.error = std.error,
    p.val.adj = p.val.adj,
    contrasts = list(
      matrix = matrix(c(0, 1), 1, 2,
                      dimnames=list("condition", parameters)),
      estimate = contrast.estimate,
      std.error = contrast.std.error,
      p.val.adj = contrast.p.val.adj
    )
  ), class="aldex")
}

test_that("aldex.plot uses coefficient results when a contrast shares its name", {
  fit <- plot_contrast_collision_fixture()
  initial <- NULL
  plotted <- NULL
  # Capture graphics in a local function environment, independent of R's
  # cached S3 plot methods or the location of the plot generic.
  plotter <- aldex.plot
  environment(plotter) <- list2env(list(
    plot = function(x, y, ...) {
      initial <<- list(x=x, y=y)
    },
    abline = function(...) NULL,
    mtext = function(...) NULL,
    points = function(x, y, ...) {
      plotted <<- list(x=x, y=y, args=list(...))
    },
    aldex.effect = function(...) {
      data.frame(parameter="condition", entity=c("taxon_a", "taxon_b"),
                 estimate=c(1, -2), pooled.SD=c(1, 2),
                 cohens.d=c(1, -1), overlap=c(0, 0))
    }
  ), parent=environment(aldex.plot))

  expect_silent(plotter(fit, contrast="condition", plot="volcano"))

  ## The base volcano plot has one point per coefficient feature, and only
  ## the coefficient's first feature is significant.  The stored contrast
  ## also has significant values, so selecting summary() without a component
  ## filter would add both its estimates and its points here.
  expect_length(initial$x, 2)
  expect_equal(initial$x, c(1, -2))
  expect_equal(plotted$x, 1)
  expect_equal(plotted$y, -log10(.01))

  initial <- NULL
  plotted <- NULL
  expect_silent(plotter(fit, contrast="condition", plot="effect"))
  expect_length(initial$x, 2)
  expect_equal(plotted$x, 1)
  expect_equal(plotted$y, 1)
})
