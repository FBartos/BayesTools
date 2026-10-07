skip_if_not_test_profile("unit")

test_that("affine point coordinates avoid collapsed translations", {

  arguments <- list(posterior = c(1, 2, 3), prior = prior("normal", list(0, 1)),
                    parameter = "theta", density_method = "normal",
                    logBF = TRUE, columns = "all")
  baseline <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 0")))
  shifted <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta + 1e100 = 1e100")))
  expect_equal(as.numeric(shifted$BF), as.numeric(baseline$BF), tolerance = 1e-12)
  shifted_diagnostics <- attr(shifted, "numerical_diagnostics")[[1L]]
  expect_equal(shifted_diagnostics$evaluation_value, 0)
  expect_equal(shifted_diagnostics$requested_value, 1e100)
  scaled <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "2 * theta = 100")))
  direct <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta = 50")))
  expect_equal(as.numeric(scaled$BF), as.numeric(direct$BF), tolerance = 1e-12)
  scaled_diagnostics <- attr(scaled, "numerical_diagnostics")[[1L]]
  direct_diagnostics <- attr(direct, "numerical_diagnostics")[[1L]]
  expect_equal(scaled_diagnostics$log_prior_height,
    direct_diagnostics$log_prior_height - log(2), tolerance = 1e-12)
  expect_equal(scaled_diagnostics$log_posterior_height,
    direct_diagnostics$log_posterior_height - log(2), tolerance = 1e-12)
  divided <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta / 49 = 1")))
  diagnostics <- attr(divided, "numerical_diagnostics")[[1L]]
  expect_equal(diagnostics$log_prior_height, stats::dnorm(49, log = TRUE) + log(49), tolerance = 1e-12)
  expect_equal(diagnostics$log_posterior_height, stats::dnorm(49, 2, 1, log = TRUE) + log(49), tolerance = 1e-12)
  expect_equal(diagnostics$evaluation_value, 49)
  expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta - theta = 0"))),
               class = "BayesTools_point_mass_at_null")
  expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = "theta * 1e-200 * 1e-200 = 0"))),
               class = "BayesTools_hypothesis_numerical_unavailable")
  for(expression in c("theta * (1e-200 * 1e-200) = 0", "theta * exp(-1000) = 0")){
    expect_error(do.call(hypothesis_BF, c(arguments, list(hypothesis = expression))),
      class = "BayesTools_hypothesis_numerical_unavailable")
  }
  constant_function <- do.call(hypothesis_BF, c(arguments, list(hypothesis = "sqrt(4) * theta = 100")))
  expect_equal(as.numeric(constant_function$BF), as.numeric(direct$BF), tolerance = 1e-12)
})

test_that("constant divisors preserve strict and inclusive atom boundaries", {

  draws <- data.frame(theta = c(48, 49, 50))
  forms <- c("theta / 49 >= 1", "theta / 49 > 1", "theta / -49 >= -1", "theta / -49 > -1")
  expected <- list(c(FALSE, TRUE, TRUE), c(FALSE, FALSE, TRUE),
                   c(TRUE, TRUE, FALSE), c(TRUE, FALSE, FALSE))
  observed <- lapply(forms, function(text){
    side <- hypothesis_parse(text)$statements[[1L]]$left
    .hypothesis_draw_region_indicator(side, draws)
  })
  expect_identical(observed, expected)
  expect_error(.hypothesis_linear_coefficients(parse(text = "theta * 1e-200 * 1e-200")[[1L]],
    "theta", draws), class = "BayesTools_hypothesis_numerical_unavailable")
})
