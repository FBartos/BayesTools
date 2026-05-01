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
})


test_that("sampled factor prior density resolves stored factor levels", {
  p <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  attr(p, "levels") <- 2

  expect_silent(d <- density(p, force_samples = TRUE, n_samples = 20, n_points = 20))
  expect_s3_class(d, "density.prior.orthonormal")
})
