skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Prior Distribution Tool Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for prior handling utilities, parameter checks, and prior type
#   detection functions.
#
# DEPENDENCIES:
#   - None (pure R)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# TAGS: @evaluation, @priors, @tools
# ============================================================================ #

test_that("Prior handling works", {


  # direct checks for positive/negative
  expect_error(BayesTools:::.check_parameter_positive(0, "par"), "The 'par' must be positive.")
  expect_null(BayesTools:::.check_parameter_positive(0, "par", TRUE))
  expect_error(BayesTools:::.check_parameter_positive(-.01, "par", TRUE), "The 'par' must be non-negative.")
  expect_error(BayesTools:::.check_parameter_negative(0, "par"), "The 'par' must be negative.")
  expect_null(BayesTools:::.check_parameter_negative(0, "par", TRUE))
  expect_error(BayesTools:::.check_parameter_negative(.01, "par", TRUE), "The 'par' must be non-positive.")

  expect_error(BayesTools:::.check_parameter_range(-.01, "par", lower = 0, upper = 1, include_bounds = TRUE), "The 'par' must be higher than 0 and lower than 1.")
  expect_error(BayesTools:::.check_parameter_range(3, "par", lower = 1, upper = 3, include_bounds = FALSE), "The 'par' must be higher or equal to than 1 and lower or equal to than 3.")


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
  expect_error(prior_weightfunction("one-sided", c(1), wf_cumulative(c(1, 1))), "'steps' must be higher than 0 and lower than 1.")
  expect_error(prior_weightfunction("one-sided", c(.05, 0.1), wf_cumulative(c(1, 1))), "length.*3")
  expect_error(prior_weightfunction("one-sided", c(.05, 0.01), wf_cumulative(c(1, 1, 1))), "'steps' must be monotonically increasing.")
  expect_error(prior_weightfunction("two-sided", c(1), wf_fixed(c(1, 1))), "'steps' must be higher than 0 and lower than 1.")
  expect_error(prior_weightfunction("two-sided", c(.05, 0.1), wf_fixed(c(1, 1))), "length.*3")
  expect_error(prior_weightfunction("two-sided", c(.05, 0.01), wf_fixed(c(1, 1, 1))), "'steps' must be monotonically increasing.")
  expect_error(prior_weightfunction("one-sided", c(1), wf_independent(prior("beta", list(1, 1)))), "'steps' must be higher than 0 and lower than 1.")
  expect_error(wf_independent(prior("normal", list(0, 1))), "non-negative support")
  expect_silent(wf_independent(prior("normal", list(0, 1)), "log_omega"))
  expect_error(prior_weightfunction("one-sided", c(.05, 0.55, .40), wf_independent(prior("beta", list(1, 1)))), "'steps' must be monotonically increasing.")

})


test_that("Truncated simple prior distribution functions respect support", {

  p <- prior("normal", list(0, 1), list(0, 2))

  expect_equal(cdf(p, c(-1, 0, 2, 3)), c(0, 0, 1, 1), tolerance = 1e-12)
  expect_equal(ccdf(p, c(-1, 0, 2, 3)), c(1, 1, 0, 0), tolerance = 1e-12)

  probabilities <- c(.1, .5, .9)
  expect_equal(cdf(p, quant(p, probabilities)), probabilities, tolerance = 1e-10)
})


test_that("point priors validate explicit truncation against their actual support", {
  expect_silent(prior("point", list(1), truncation = list(0, 2)))
  expect_error(
    prior("point", list(1), truncation = list(2, 3)),
    "must contain the point location"
  )

  expect_silent(prior("mpoint", list(1, 2), truncation = list(0, 2)))
  expect_error(
    prior("mpoint", list(1, 2), truncation = list(2, 3)),
    "must contain the point location"
  )
})


test_that("vector prior K validation allows deferred dimensions but rejects invalid values", {
  expect_silent(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"))
  expect_error(prior("mnormal", list(0, 1, 2.5)), "integer")
  expect_error(prior("mt", list(0, 1, 3, -2)), "positive integer")
})


test_that("mixture coercion preserves component prior weights", {
  simple <- prior_mixture(list(
    prior_none(prior_weights = 2),
    prior("normal", list(0, 1), prior_weights = 7)
  ))
  expect_equal(attr(simple, "prior_weights"), c(2, 7))
  expect_equal(vapply(simple, function(x) x$prior_weights, numeric(1)), c(2, 7))

  factor <- prior_mixture(list(
    prior("point", list(0), prior_weights = 9),
    prior_factor("mnormal", list(0, 1), contrast = "orthonormal", prior_weights = 1)
  ))
  expect_equal(attr(factor, "prior_weights"), c(9, 1))
  expect_equal(vapply(factor, function(x) x$prior_weights, numeric(1)), c(9, 1))
})


test_that("prior constructors reject invalid model weights", {
  expect_error(
    prior("normal", list(0, 1), prior_weights = 0),
    "higher than 0"
  )
  expect_error(
    prior("normal", list(0, 1), prior_weights = Inf),
    "must be finite"
  )
  expect_error(
    prior_none(prior_weights = NA_real_),
    "cannot contain NA/NaN"
  )
  expect_error(
    prior_weightfunction(
      side = "one-sided",
      steps = c(.05),
      weights = wf_cumulative(c(1, 1)),
      prior_weights = NaN
    ),
    "cannot contain NA/NaN"
  )
  expect_error(
    prior_spike_and_slab(
      prior_parameter = prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(.5)),
      prior_weights = Inf
    ),
    "must be finite"
  )
})


test_that("raw vector-prior marginal generics return numeric arrays", {
  p_mnormal <- prior("mnormal", list(0, 1, 2))
  expect_equal(quant(p_mnormal, .5), matrix(0, nrow = 1, ncol = 2))
  expect_equal(mcdf(p_mnormal, 0), matrix(.5, nrow = 1, ncol = 2))
  expect_equal(mccdf(p_mnormal, 0), matrix(.5, nrow = 1, ncol = 2))
  expect_equal(mlpdf(p_mnormal, 0), matrix(stats::dnorm(0, log = TRUE), nrow = 1, ncol = 2))
  expect_equal(mean(p_mnormal), c(0, 0))
  expect_equal(var(p_mnormal), c(1, 1))
  expect_equal(sd(p_mnormal), c(1, 1))

  p_mpoint <- prior("mpoint", list(2, 3))
  expect_equal(mquant(p_mpoint, .5), matrix(2, nrow = 1, ncol = 3))
  expect_equal(mcdf(p_mpoint, 2), matrix(1, nrow = 1, ncol = 3))
})


test_that("is.prior functions work correctly", {

  # Create priors for testing
  p_normal <- prior("normal", list(0, 1))
  p_point  <- prior("point", list(0))
  p_discrete <- prior("bernoulli", list(0.5))
  p_vector <- prior("mnormal", list(0, 1, 3))
  p_pet    <- prior_PET("normal", list(0, 1))
  p_peese  <- prior_PEESE("normal", list(0, 1))
  p_wf     <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))
  p_factor_t <- prior_factor("normal", contrast = "treatment", list(0, 1))
  p_factor_o <- prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
  p_factor_m <- prior_factor("mnormal", contrast = "meandif", list(0, 1))
  p_factor_i <- prior_factor("beta", contrast = "independent", list(1, 1))
  p_spike_slab <- prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("spike", list(0.5)))
  p_none <- prior_none()

  # Test is.prior
  expect_true(is.prior(p_normal))
  expect_true(is.prior(p_point))
  expect_false(is.prior("not a prior"))
  expect_false(is.prior(list(a = 1)))

  # Test is.prior.simple
  expect_true(is.prior.simple(p_normal))
  expect_true(is.prior.simple(p_point))
  expect_false(is.prior.simple(p_wf))
  expect_false(is.prior.simple(p_vector))

  # Test is.prior.point
  expect_true(is.prior.point(p_point))
  expect_false(is.prior.point(p_normal))

  # Test is.prior.none
  expect_true(is.prior.none(p_none))
  expect_false(is.prior.none(p_normal))

  # Test is.prior.discrete
  expect_true(is.prior.discrete(p_discrete))
  expect_false(is.prior.discrete(p_normal))

  # Test is.prior.vector
  expect_true(is.prior.vector(p_vector))
  expect_false(is.prior.vector(p_normal))

  # Test is.prior.PET
  expect_true(is.prior.PET(p_pet))
  expect_false(is.prior.PET(p_normal))

  # Test is.prior.PEESE
  expect_true(is.prior.PEESE(p_peese))
  expect_false(is.prior.PEESE(p_normal))

  # Test is.prior.weightfunction
  expect_true(is.prior.weightfunction(p_wf))
  expect_false(is.prior.weightfunction(p_normal))

  # Test is.prior.factor
  expect_true(is.prior.factor(p_factor_t))
  expect_true(is.prior.factor(p_factor_o))
  expect_true(is.prior.factor(p_factor_m))
  expect_true(is.prior.factor(p_factor_i))
  expect_false(is.prior.factor(p_normal))

  # Test is.prior.treatment
  expect_true(is.prior.treatment(p_factor_t))
  expect_false(is.prior.treatment(p_factor_o))

  # Test is.prior.orthonormal
  expect_true(is.prior.orthonormal(p_factor_o))
  expect_false(is.prior.orthonormal(p_factor_t))

  # Test is.prior.meandif
  expect_true(is.prior.meandif(p_factor_m))
  expect_false(is.prior.meandif(p_factor_o))

  # Test is.prior.independent
  expect_true(is.prior.independent(p_factor_i))
  expect_false(is.prior.independent(p_factor_t))

  # Test is.prior.spike_and_slab
  expect_true(is.prior.spike_and_slab(p_spike_slab))
  expect_false(is.prior.spike_and_slab(p_normal))

})


test_that(".check_prior works correctly", {

  p_normal <- prior("normal", list(0, 1))

  # Valid prior should pass
  expect_null(BayesTools:::.check_prior(p_normal))

  # Non-prior should fail
  expect_error(BayesTools:::.check_prior("not a prior"), "must be a valid prior object")
  expect_error(BayesTools:::.check_prior(list(a = 1)), "must be a valid prior object")

})


test_that(".check_prior_list works correctly", {

  p_normal <- prior("normal", list(0, 1))
  p_point  <- prior("point", list(0))
  p_wf     <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))

  # Valid prior list should pass
  expect_null(BayesTools:::.check_prior_list(list(p_normal, p_point)))

  # Empty list with allow_NULL
  expect_null(BayesTools:::.check_prior_list(NULL, allow_NULL = TRUE))
  expect_null(BayesTools:::.check_prior_list(list(), allow_NULL = TRUE))

  # Non-list should fail
  expect_error(BayesTools:::.check_prior_list("not a list"), "must be a list")

  # List with non-prior should fail
  expect_error(BayesTools:::.check_prior_list(list("not a prior")), "must be a prior distribution")

  # Disallowing specific types
  expect_error(
    BayesTools:::.check_prior_list(list(p_point), allow_prior.point = FALSE),
    "must not contain point priors"
  )
  expect_error(
    BayesTools:::.check_prior_list(list(p_wf), allow_prior.weightfunction = FALSE),
    "must not contain weightfunction priors"
  )

})


test_that("prior-like objects missing classes or fields are rejected at public boundaries", {
  p_normal <- prior("normal", list(0, 1))

  classless_prior <- p_normal
  class(classless_prior) <- "list"
  expect_error(
    prior_mixture(list(classless_prior, p_normal)),
    "must be a valid prior object"
  )
  expect_error(
    prior_spike_and_slab(classless_prior),
    "must be a prior distribution"
  )

  classless_weights <- unclass(wf_cumulative(c(1, 1)))
  expect_error(
    prior_weightfunction("one-sided", c(.05), classless_weights),
    "created by wf_cumulative\\(\\), wf_fixed\\(\\), or wf_independent\\(\\)"
  )

  missing_omega <- wf_fixed(c(1, 1))
  missing_omega$omega <- NULL
  expect_error(
    prior_weightfunction("one-sided", c(.05), missing_omega),
    "omega.*cannot be NULL"
  )

  unsupported_type <- wf_fixed(c(1, 1))
  unsupported_type$type <- "unsupported"
  expect_error(
    prior_weightfunction("one-sided", c(.05), unsupported_type),
    "Unsupported weightfunction weight prior type"
  )
})


test_that(".check_and_name_parameters works correctly", {

  # Valid parameters
  params <- list(0, 1)
  result <- BayesTools:::.check_and_name_parameters(params, c("mean", "sd"), "normal")
  expect_equal(names(result), c("mean", "sd"))

  # Named parameters in different order
  params2 <- list(sd = 1, mean = 0)
  result2 <- BayesTools:::.check_and_name_parameters(params2, c("mean", "sd"), "normal")
  expect_equal(result2$mean, 0)
  expect_equal(result2$sd, 1)

  # Wrong number of parameters
  expect_error(
    BayesTools:::.check_and_name_parameters(list(0), c("mean", "sd"), "normal"),
    "requires 2 parameters"
  )

  # Invalid parameter names
  expect_error(
    BayesTools:::.check_and_name_parameters(list(location = 0, sd = 1), c("mean", "sd"), "normal"),
    "Parameters 'location' are not supported"
  )

})


test_that(".check_and_set_truncation works correctly", {

  # Default truncation
  result <- BayesTools:::.check_and_set_truncation(list())
  expect_equal(result$lower, -Inf)
  expect_equal(result$upper, Inf)

  # Named truncation
  result2 <- BayesTools:::.check_and_set_truncation(list(lower = 0))
  expect_equal(result2$lower, 0)
  expect_equal(result2$upper, Inf)

  result3 <- BayesTools:::.check_and_set_truncation(list(upper = 1))
  expect_equal(result3$lower, -Inf)
  expect_equal(result3$upper, 1)

  # Positional truncation
  result4 <- BayesTools:::.check_and_set_truncation(list(0, 1))
  expect_equal(result4$lower, 0)
  expect_equal(result4$upper, 1)

  # Single positional (becomes lower)
  result5 <- BayesTools:::.check_and_set_truncation(list(0))
  expect_equal(result5$lower, 0)

  # Distribution-specific defaults
  result6 <- BayesTools:::.check_and_set_truncation(list(), lower = 0)
  expect_equal(result6$lower, 0)

  # Error conditions
  expect_error(
    BayesTools:::.check_and_set_truncation(list(1, 2, 3)),
    "More than two truncation points"
  )

  expect_error(
    BayesTools:::.check_and_set_truncation(list(bad_name = 0)),
    "must be named 'lower' and 'upper'"
  )

  expect_error(
    BayesTools:::.check_and_set_truncation(list(1, 0)),
    "lower truncation point must be lower"
  )

  expect_error(
    BayesTools:::.check_and_set_truncation(list(-1, Inf), lower = 0),
    "Lower truncation point must be larger or equal to 0"
  )

  expect_error(
    BayesTools:::.check_and_set_truncation(list(-Inf, 2), upper = 1),
    "Upper truncation point must be smaller or equal to 1"
  )

})


test_that(".check_parameter works correctly", {

  # Valid parameters
  expect_null(BayesTools:::.check_parameter(1, "param"))
  expect_null(BayesTools:::.check_parameter(c(1, 2, 3), "param", length = 3))
  expect_null(BayesTools:::.check_parameter(c(1, 2), "param", length = 0))

  # Expressions should pass through
  expect_null(BayesTools:::.check_parameter(expression(x + 1), "param"))

  # Invalid parameters
  expect_error(
    BayesTools:::.check_parameter("a", "param"),
    "must be a numeric vector"
  )

  expect_error(
    BayesTools:::.check_parameter(c(1, 2), "param", length = 3),
    "must be a numeric vector of length 3"
  )

})


test_that(".check_parameter_dimensions works correctly", {

  # Valid dimensions
  expect_null(BayesTools:::.check_parameter_dimensions(3, "K"))
  expect_null(BayesTools:::.check_parameter_dimensions(NA, "K", allow_NA = TRUE))

  # Expressions should pass through
  expect_null(BayesTools:::.check_parameter_dimensions(expression(K), "K"))

  # Invalid dimensions
  expect_error(
    BayesTools:::.check_parameter_dimensions(NA, "K", allow_NA = FALSE),
    "must be defined"
  )

  # Note: The function has some implementation quirks with vector input,
  # so we just test that invalid inputs throw some error
  expect_error(BayesTools:::.check_parameter_dimensions(c(1, 2), "K"))
  expect_error(BayesTools:::.check_parameter_dimensions(1.5, "K"))

})


test_that(".get_prior_factor_levels works correctly", {

  # Treatment contrast - levels - 1
  p_treatment <- prior_factor("normal", contrast = "treatment", list(0, 1))
  attr(p_treatment, "levels") <- 3
  expect_equal(BayesTools:::.get_prior_factor_levels(p_treatment), 2)

  # Independent contrast - all levels
  p_independent <- prior_factor("beta", contrast = "independent", list(1, 1))
  attr(p_independent, "levels") <- 3
  expect_equal(BayesTools:::.get_prior_factor_levels(p_independent), 3)

  # Orthonormal contrast - levels - 1
  p_orthonormal <- prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
  attr(p_orthonormal, "levels") <- 4
  expect_equal(BayesTools:::.get_prior_factor_levels(p_orthonormal), 3)

  # Meandif contrast - levels - 1
  p_meandif <- prior_factor("mnormal", contrast = "meandif", list(0, 1))
  attr(p_meandif, "levels") <- 3
  expect_equal(BayesTools:::.get_prior_factor_levels(p_meandif), 2)

})


test_that(".prior_clean_input_name works correctly", {

  expect_equal(BayesTools:::.prior_clean_input_name("Normal"), "normal")
  expect_equal(BayesTools:::.prior_clean_input_name("Log-Normal"), "lognormal")
  expect_equal(BayesTools:::.prior_clean_input_name("Student_t"), "studentt")
  expect_equal(BayesTools:::.prior_clean_input_name("one.sided"), "onesided")
  expect_equal(BayesTools:::.prior_clean_input_name("two-sided"), "twosided")
  expect_equal(BayesTools:::.prior_clean_input_name(" Sp ace "), "space")

})
