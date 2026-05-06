skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Prior Plot Data Semantics
# ============================================================================ #
#
# PURPOSE:
#   Semantic ggplot-layer tests for plot.prior outputs. These complement the
#   visual-profile vdiffr snapshots in test-priors-plot.R.
#
# TAGS: @evaluation, @unit, @priors, @plots
# ============================================================================ #

prior_plot_layer_data <- function(plot) {
  ggplot2::ggplot_build(plot)$data
}

test_that("continuous prior ggplot data matches analytic density and labels", {
  x <- seq(-1, 1, length.out = 5)
  p <- prior("normal", list(mean = 0, sd = 1))

  g <- plot(
    p,
    plot_type = "ggplot",
    x_seq = x,
    xlab = "Effect",
    ylab = "Mass density",
    main = "Normal prior"
  )
  layers <- prior_plot_layer_data(g)

  expect_s3_class(g, "ggplot")
  expect_length(layers, 1)
  expect_equal(layers[[1]]$x, x)
  expect_equal(layers[[1]]$y, stats::dnorm(x), tolerance = 1e-12)
  expect_equal(g$labels$title, "Normal prior")
  expect_equal(g$scales$get_scales("x")$name, "Effect")
  expect_equal(g$scales$get_scales("y")$name, "Mass density")
})

test_that("continuous prior ggplot density preserves selected x grid and bounded mass", {
  x <- seq(-4, 4, length.out = 401)
  p <- prior("normal", list(mean = 0, sd = 1))

  layer <- prior_plot_layer_data(plot(p, plot_type = "ggplot", x_seq = x))[[1]]
  dx <- diff(layer$x)
  trapezoid_mass <- sum(dx * (head(layer$y, -1) + tail(layer$y, -1)) / 2)

  expect_equal(layer$x, x)
  expect_true(all(dx > 0))
  expect_equal(trapezoid_mass, diff(stats::pnorm(c(-4, 4))), tolerance = 1e-4)
  expect_true(all(layer$y >= 0))
  expect_equal(layer$y[which.max(layer$y)], stats::dnorm(0), tolerance = 1e-12)
})

test_that("transformed continuous prior ggplot data applies the Jacobian", {
  original_x <- seq(-1, 1, length.out = 5)
  transformed_x <- exp(original_x)

  g <- plot(
    prior("normal", list(mean = 0, sd = 1)),
    plot_type = "ggplot",
    x_seq = original_x,
    transformation = "exp"
  )
  layer <- prior_plot_layer_data(g)[[1]]

  expect_equal(layer$x, transformed_x)
  expect_equal(layer$y, stats::dnorm(original_x) / transformed_x, tolerance = 1e-12)
  expect_true(all(layer$x > 0))
})

test_that("discrete and point prior ggplot layers preserve exact probability mass", {
  discrete_layer <- prior_plot_layer_data(
    plot(prior("bernoulli", list(.33)), plot_type = "ggplot")
  )[[1]]

  expect_equal(discrete_layer$x, c(0, 1))
  expect_equal(discrete_layer$y, c(.67, .33), tolerance = 1e-12)
  expect_probability_mass(discrete_layer$y)

  point_layer <- prior_plot_layer_data(
    plot(prior("point", list(1.5)), plot_type = "ggplot")
  )[[1]]

  expect_equal(point_layer$x, 1.5)
  expect_equal(point_layer$xend, 1.5)
  expect_equal(point_layer$y, 0)
  expect_equal(point_layer$yend, 1)
})

test_that("point prior ggplot layer honors explicit probability scaling", {
  point_layer <- prior_plot_layer_data(
    plot(prior("point", list(location = -1)), plot_type = "ggplot", scale_y2 = 3)
  )[[1]]

  expect_equal(point_layer$x, -1)
  expect_equal(point_layer$xend, -1)
  expect_equal(point_layer$y, 0)
  expect_equal(point_layer$yend, 3)
  expect_probability_mass(point_layer$yend / 3)
})

test_that("spike-and-slab prior ggplot layers separate slab density and spike mass", {
  x <- seq(-2, 2, length.out = 7)
  p <- prior_spike_and_slab(prior("normal", list(mean = 0, sd = 1)))

  g <- plot(p, plot_type = "ggplot", x_seq = x)
  layers <- prior_plot_layer_data(g)

  expect_length(layers, 2)
  expect_equal(layers[[1]]$x, x)
  expect_equal(layers[[1]]$y, .5 * stats::dnorm(x), tolerance = 1e-12)
  expect_equal(layers[[2]]$x, 0)
  expect_equal(layers[[2]]$xend, 0)
  expect_equal(layers[[2]]$y, 0)
  expect_equal(layers[[2]]$yend / attr(g, "scale_y2"), .5, tolerance = 1e-12)
})

test_that("weightfunction prior ggplot layers preserve step locations and intervals", {
  p <- prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .25)))
  layers <- prior_plot_layer_data(plot(p, plot_type = "ggplot"))

  expect_length(layers, 2)

  expected_x <- c(0, .05, .05, 1)
  expected_y <- c(1, 1, .25, .25)
  expected_lower <- expected_y
  expected_upper <- expected_y

  expect_equal(layers[[1]]$x, c(expected_x, rev(expected_x)))
  expect_equal(layers[[1]]$y, c(expected_lower, rev(expected_upper)), tolerance = 1e-12)
  expect_equal(layers[[2]]$x, expected_x)
  expect_equal(layers[[2]]$y, expected_y, tolerance = 1e-12)
})

test_that("factor prior ggplot data uses factor-contrast density scale", {
  p <- prior_factor("normal", list(mean = 0, sd = 2), contrast = "treatment")
  p$parameters$K <- 3
  x <- seq(-4, 4, length.out = 9)

  layer <- prior_plot_layer_data(plot(p, plot_type = "ggplot", x_seq = x))[[1]]

  expect_equal(layer$x, x)
  expect_equal(layer$y, stats::dnorm(x, mean = 0, sd = 2), tolerance = 1e-12)
})

test_that("geom_prior contributes exact density and point-mass layers", {
  x <- seq(-1, 1, length.out = 5)

  density_plot <- ggplot2::ggplot() +
    geom_prior(prior("normal", list(mean = 0, sd = 1)), x_seq = x)
  density_layer <- ggplot2::ggplot_build(density_plot)$data[[1]]

  expect_equal(density_layer$x, x)
  expect_equal(density_layer$y, stats::dnorm(x), tolerance = 1e-12)

  point_plot <- ggplot2::ggplot() +
    geom_prior(prior("point", list(location = 0.5)), scale_y2 = 2)
  point_layer <- ggplot2::ggplot_build(point_plot)$data[[1]]

  expect_equal(point_layer$x, 0.5)
  expect_equal(point_layer$xend, 0.5)
  expect_equal(point_layer$y, 0)
  expect_equal(point_layer$yend, 1)
})

test_that("prior plot data rejects invalid plotting options before rendering", {
  p <- prior("normal", list(mean = 0, sd = 1))

  expect_error(plot(p, plot_type = "grid"), "'plot_type'")
  expect_error(plot(p, plot_type = "ggplot", show_figures = 1.5), "'show_figures'")
  expect_error(
    plot(p, plot_type = "ggplot", transformation = "not-a-transformation"),
    "function or character string"
  )
  expect_error(ggplot2::ggplot() + geom_prior(p, scale_y2 = -1), "'scale_y2'")
})
