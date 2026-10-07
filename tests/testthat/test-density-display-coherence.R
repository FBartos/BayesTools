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
