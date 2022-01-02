context("Prior distribution tool functions")

test_that("Prior handling works", {


  # direct checks for positive/negative
  expect_error(BayesTools:::.check_parameter_positive(0, "par"), "The 'par' must be positive.")
  expect_null(BayesTools:::.check_parameter_positive(0, "par", TRUE))
  expect_error(BayesTools:::.check_parameter_positive(-.01, "par", TRUE), "The 'par' must be non-negative.")
  expect_error(BayesTools:::.check_parameter_negative(0, "par"), "The 'par' must be negative.")
  expect_null(BayesTools:::.check_parameter_negative(0, "par", TRUE))
  expect_error(BayesTools:::.check_parameter_negative(.01, "par", TRUE), "The 'par' must be non-positive.")

  # check different ordering of names withing lists
  expect_equal(prior("normal", list(0, 1)), prior("normal", list(0, 1), list(-Inf, Inf)))
  expect_equal(prior("normal", list(0, 1), list(0, Inf)), prior("normal", list(0, 1), list("lower" = 0)))
  expect_equal(prior("normal", list(0, 1), list(-Inf, 1)), prior("normal", list(0, 1), list("upper" = 1)))
  expect_equal(prior("normal", list(0, 1), list("upper" = Inf, lower = -Inf)), prior("normal", list(0, 1), list(-Inf, Inf)))
  expect_equal(prior("normal", list(0, 1)), prior("normal", list("sd" = 1, "mean" = 0)))

  # check errors messages
  expect_error(prior("normal", list(0)), "normal prior distribution requires 2 parameters.")
  expect_error(prior("normal", list(c(0, 0), 1)), "The 'mean' must be a numeric vector of length 1.")
  expect_error(prior("normal", list("a", 1)), "The 'mean' must be a numeric vector of length 1.")
  expect_error(prior("normal", list(0, "location" = 1)), "Parameters 'location' are not supported for a normal distribution.")
  expect_error(prior("normal", list(0, 1), list(Inf, -Inf)), "The lower truncation point must be lower than the upper truncation points.")
  expect_error(prior("lognormal", list(0, 1), list(-5, Inf)), "Lower truncation point must be larger or equal to 0.")
  expect_error(prior("beta", list(0, 1), list(0, 2)), "Upper truncation point must be smaller or equal to 1.")
  expect_error(prior("normal", list(0, -1)), "The 'sd' must be positive.")
  expect_error(prior_weightfunction("one-sided", list(c(1), c(1, 1))), "Parameter 'steps' must be higher than 0 and lower than 1.")
  expect_error(prior_weightfunction("one-sided", list(c(.05, 0.1), c(1, 1))), "The parameter alpha needs to have one more argument then there are steps.")
  expect_error(prior_weightfunction("one-sided", list(c(.05, 0.01), c(1, 1))), "Parameters 'steps' must be monotonically increasing.")
  expect_error(prior_weightfunction("two-sided.fixed", list(c(1), c(1, 1))), "Parameter 'steps' must be higher than 0 and lower than 1.")
  expect_error(prior_weightfunction("two-sided.fixed", list(c(.05, 0.1), c(1, 1))), "The parameter omega needs to have one more argument then there are steps.")
  expect_error(prior_weightfunction("two-sided.fixed", list(c(.05, 0.01), c(1, 1))), "Parameters 'steps' must be monotonically increasing.")
  expect_error(prior_weightfunction("one-sided", list(c(1), c(1, 1), c(1, 1))), "Parameter 'steps' must be higher than 0 and lower than 1.")
  expect_error(prior_weightfunction("one-sided", list(c(0.05, 0.10), c(1, 1), c(1, 1))), "The parameter alpha1 needs to have one more argument then there are steps <= .5.")
  expect_error(prior_weightfunction("one-sided", list(c(0.05, 0.55, 0.65), c(1, 1), c(1, 1))), "The parameter alpha2 needs to have one more argument then there are steps > .5.")
  expect_error(prior_weightfunction("one-sided", list(c(.05, 0.55, .40), c(1, 1), c(1, 1))), "Parameters 'steps' must be monotonically increasing.")

})
