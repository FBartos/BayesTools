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
