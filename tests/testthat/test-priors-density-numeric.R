skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Prior Density Numeric Regression Tests
# ============================================================================ #

test_that("density transformation Jacobians use absolute inverse derivatives", {
  expect_equal(
    BayesTools:::.density.prior_transformation_y(
      x = c(-2, -4),
      y = c(10, 20),
      transformation = "lin",
      transformation_arguments = list(a = 0, b = -2)
    ),
    c(5, 10)
  )

  transformed_x <- exp(log(c(2, 4)) * 2)
  expect_equal(
    BayesTools:::.density.prior_transformation_y(
      x = transformed_x,
      y = c(10, 20),
      transformation = "exp_lin",
      transformation_arguments = list(a = 0, b = 2)
    ),
    c(10, 20) * c(2, 4) / (2 * transformed_x)
  )

  transformed_x_dec <- exp(-log(c(2, 4)))
  expect_equal(
    BayesTools:::.density.prior_transformation_y(
      x = transformed_x_dec,
      y = c(10, 20),
      transformation = "exp_lin",
      transformation_arguments = list(a = 0, b = -1)
    ),
    c(10, 20) * c(2, 4) / transformed_x_dec
  )

  exp_transformation <- list(fun = exp, inv = log, jac = exp)
  expect_equal(
    BayesTools:::.density.prior_transformation_y(
      x = c(2, 4),
      y = c(10, 20),
      transformation = exp_transformation
    ),
    c(10, 20) / c(2, 4)
  )
})


test_that("sampled factor prior density resolves stored factor levels", {
  p <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  attr(p, "levels") <- 2

  expect_silent(d <- density(p, force_samples = TRUE, n_samples = 20, n_points = 20))
  expect_s3_class(d, "density.prior.orthonormal")
})


test_that("density defaults and forced sampling keep x and y aligned", {
  discrete <- density(prior("bernoulli", list(.3)), truncate_end = FALSE)
  expect_equal(discrete$x, c(0, 1))
  expect_equal(discrete$y, c(.7, .3))

  set.seed(1)
  sampled <- density(
    prior("normal", list(0, 1)),
    x_seq = seq(-1, 1, length.out = 5),
    n_points = 25,
    force_samples = TRUE,
    n_samples = 200
  )
  expect_length(sampled$x, 25)
  expect_length(sampled$y, 25)

  set.seed(2)
  sampled_discrete <- density(
    prior("bernoulli", list(.3)),
    force_samples = TRUE,
    n_samples = 200,
    truncate_end = FALSE
  )
  expect_equal(sampled_discrete$x, c(0, 1))
  expect_equal(sum(sampled_discrete$y), 1)
})


test_that("density rejects invalid range and transformation shapes", {
  p <- prior("normal", list(0, 1))

  expect_error(
    density(p, x_range = c(1, -1)),
    "lower range limit"
  )
  expect_error(
    density(
      p,
      x_seq = 0,
      transformation = list(fun = exp, inv = log, jac = exp, extra = identity)
    ),
    "must have length '3'"
  )
})


test_that("exp transformed range settings are interpreted on output scale", {
  transformed_range <- exp(c(-1, 1))

  d <- density(
    prior("normal", list(0, 1)),
    x_range = transformed_range,
    n_points = 5,
    transformation = "exp",
    transformation_settings = TRUE,
    truncate_end = FALSE
  )

  expected_x <- seq(transformed_range[1], transformed_range[2], length.out = 5)

  expect_equal(d$x, expected_x)
  expect_equal(attr(d, "x_range"), transformed_range)
  expect_equal(d$y, stats::dnorm(log(expected_x)) / expected_x, tolerance = 1e-12)
  expect_true(all(d$x > 0))
})


test_that("exp transformed simple prior preserves continuous mass", {
  d <- density(
    prior("uniform", list(0, 1)),
    x_range = c(0, 1),
    n_points = 1001,
    transformation = "exp",
    truncate_end = FALSE
  )

  transformed_mass <- sum(diff(d$x) * (head(d$y, -1) + tail(d$y, -1)) / 2)

  expect_equal(attr(d, "x_range"), c(1, exp(1)))
  expect_equal(d$y, 1 / d$x, tolerance = 1e-12)
  expect_equal(transformed_mass, 1, tolerance = 1e-6)
})


test_that("spike-and-slab density scales slab density and spike mass", {
  p <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(mean = 1, sd = 2)),
    prior_inclusion = prior("point", list(location = .25))
  )
  x_seq <- c(-1, 1, 3)

  d <- density(
    p,
    x_seq = x_seq,
    force_samples = FALSE,
    truncate_end = FALSE
  )

  expect_s3_class(d, "density.prior.spike_and_slab")
  expect_equal(d$variable$x, x_seq)
  expect_equal(d$variable$y, .25 * stats::dnorm(x_seq, mean = 1, sd = 2))
  expect_equal(d$inclusion$x, 0)
  expect_equal(d$inclusion$y, .75)
  expect_equal(attr(d$variable, "y_range"), c(0, max(d$variable$y)))
  expect_equal(attr(d$inclusion, "y_range"), c(0, .75))
  expect_equal(attr(d, "y_range_variable"), attr(d$variable, "y_range"))
  expect_equal(attr(d, "y_range_inclusion"), attr(d$inclusion, "y_range"))
})


test_that("transformed spike-and-slab density reports transformed top-level range", {
  d <- density(
    prior_spike_and_slab(prior("normal", list(0, 1))),
    transformation = "exp"
  )

  component_range <- range(c(
    attr(d$variable, "x_range"),
    attr(d$inclusion, "x_range")
  ))
  expect_equal(attr(d, "x_range"), component_range)
  expect_true(attr(d, "x_range")[1] > 0)
})
