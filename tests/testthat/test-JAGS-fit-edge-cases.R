# ============================================================================ #
# TEST FILE: JAGS Fit Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for JAGS fitting functions including input validation,
#   error handling, and boundary conditions.
#
# DEPENDENCIES:
#   - rjags: For JAGS model syntax generation and testing
#   - common-functions.R: REFERENCE_DIR, test_reference_text, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"): For all tests
#
# MODELS/FIXTURES:
#   - Some tests use pre-fitted models from test-00-model-fits.R
#
# TAGS: @edge-cases, @JAGS, @input-validation
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-fit-edge-cases")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: Input validation tests
# ============================================================================ #
test_that("JAGS_add_priors input validation works", {

  # Empty prior_list returns original syntax
  expect_equal(JAGS_add_priors("model{}", list()), "model{}")

  # prior_list must be a list of priors
  expect_error(JAGS_add_priors("model{}", list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_add_priors("model{}", prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that("JAGS_get_inits input validation works", {

  # Empty prior_list returns empty list
  expect_equal(JAGS_get_inits(list(), chains = 2, seed = 1), list())

  # Input validation
  expect_error(JAGS_get_inits(list(x = 1), chains = 2, seed = 1), "'prior_list' must be a list of priors.")
  expect_error(JAGS_get_inits(prior("normal", list(0, 1)), chains = 2, seed = 1), "'prior_list' must be a list of priors.")

})


test_that("JAGS_to_monitor input validation works", {

  # Empty prior_list returns empty string
  expect_equal(JAGS_to_monitor(list()), "")

  # Input validation
  expect_error(JAGS_to_monitor(list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_to_monitor(prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that(".check_JAGS_syntax validates syntax correctly", {

  # Test with valid syntax
  expect_silent(JAGS_add_priors("model{}", list(mu = prior("normal", list(0, 1)))))

  # Test with missing "model" keyword
  expect_error(
    JAGS_add_priors("invalid{}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing opening brace
  expect_error(
    JAGS_add_priors("model}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing closing brace
  expect_error(
    JAGS_add_priors("model{", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with non-character input
  expect_error(
    JAGS_add_priors(123, list(mu = prior("normal", list(0, 1)))),
    "must be a character"
  )

})


test_that("JAGS_extend error handling", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Test error when fit is not a BayesTools_fit
  expect_error(
    JAGS_extend(list(), autofit_control = list()),
    "'fit' must be a 'BayesTools_fit'"
  )

})


# ============================================================================ #
# SECTION 2: Convergence edge cases
# ============================================================================ #
test_that("JAGS_check_convergence handles single chain (R-hat warning)", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 1,  # Single chain - R-hat cannot be computed
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent.jags = TRUE
  ))

  # Should warn about single chain R-hat
  expect_warning(
    JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = 1.05),
    "Only one chain was run"
  )

})


test_that("JAGS_check_convergence handles ESS and error checks", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 2,
    adapt = 50,
    burnin = 50,
    sample = 50,  # Small sample for testing convergence failures
    silent.jags = TRUE
  ))

  # Test with very strict ESS requirement (should fail)
  result_ess <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = 10000, max_error = NULL, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_ess)
  expect_true(!is.null(attr(result_ess, "errors")))

  # Test with very strict error requirement
  result_err <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = 0.00001, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_err)

  # Test with very strict SD error requirement
  result_sd <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = 0.00001, fail_fast = FALSE)
  expect_false(result_sd)

})


# ============================================================================ #
# SECTION 3: JAGS_fit with is_JASP mode
# ============================================================================ #
test_that("JAGS_fit works with is_JASP mode", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Simple model for testing is_JASP mode
  set.seed(1)
  data <- list(
    y = rnorm(20, 0.5, 1),
    N = 20
  )

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1), list(0, Inf))
  )

  model_syntax <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu, 1/pow(sigma, 2))
    }
  }"

  # Mock JASP progress bar functions (they should be skipped if not available)
  # The is_JASP mode should work but simply skip progress bars if functions don't exist
  fit_jasp <- capture.output(tryCatch({
    suppressWarnings(JAGS_fit(
      model_syntax = model_syntax,
      data = data,
      prior_list = prior_list,
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      seed = 1,
      silent = TRUE,
      is_JASP = TRUE,
      is_JASP_prefix = "Test"
    ))
  }, error = function(e) {
    # If JASP functions don't exist, this should still produce a fit
    # or fail gracefully
    if (grepl("JASP", e$message)) {
      skip("JASP progress bar functions not available")
    }
    stop(e)
  }))

  test_reference_text(paste0(fit_jasp, collapse = ","), "fit_jasp.txt")

})
