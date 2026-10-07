skip_if_not_test_profile("unit")

test_that("narrow huge finite pieces retain physical geometry and kernel scale", {

  bounds <- c(1e300, 1e300 * (1 + 1e-12))
  slope <- prior("normal", list(0, 1e-300))
  attr(slope, "multiply_by") <- "scale"
  law <- .prior_linear_combination_density(list(slope = slope,
    scale = prior("uniform", as.list(bounds))), c(slope = 1), n_grid = 8192)
  reference <- stats::integrate(function(probability){
    scale <- bounds[1L] + diff(bounds) * probability
    stats::dnorm(1, sd = 1e-300 * scale)
  }, 0, 1, rel.tol = 1e-12)$value
  result <- prior_density_ordinate(law, 1)
  expect_identical(result$behavior, "regular")
  expect_lt(abs(exp(result$log_density) / reference - 1), 1e-4)
  expect_true(all(result$provenance$integration$piece_evaluations <= 8192L))
})

test_that("finite inverse region loss refuses while genuine support events remain exact", {

  law <- .prior_linear_combination_density(list(theta = prior("gamma", list(.001, 1))),
    c(theta = 1e30), n_grid = 512)
  region <- .hypothesis_prior_region(quote(theta < 1e-300), "theta")
  expect_error(.prior_linear_density_region_probability(law, region),
    class = "BayesTools_numerical_condition")
  side <- hypothesis_parse("theta < 1e-300")$statements[[1L]]$left
  expect_error(.hypothesis_prior_density_prob(law, side, "theta"),
    class = "BayesTools_hypothesis_region")
  empty <- .hypothesis_prior_region(quote(theta < -1e-300), "theta")
  expect_equal(as.numeric(.prior_linear_density_region_probability(law, empty)), 0)
  ordinary <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))),
    c(theta = 2), n_grid = 512)
  expect_equal(as.numeric(.prior_linear_density_region_probability(ordinary,
    .hypothesis_prior_region(quote(theta < 1), "theta"))), stats::pnorm(.5))
})
