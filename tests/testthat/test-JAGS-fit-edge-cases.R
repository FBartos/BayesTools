# ============================================================================ #
# TEST FILE: JAGS Fit Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case and comprehensive tests for JAGS fitting functions including
#   JAGS_add_priors, JAGS_fit, and related utilities.
#
# DEPENDENCIES:
#   - rjags: For JAGS model syntax generation and testing
#   - common-functions.R: REFERENCE_DIR, test_reference_text, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"): For all tests
#   - skip_if_no_fits(): For tests using pre-fitted models
#
# MODELS/FIXTURES:
#   - Some tests use pre-fitted models from test-00-model-fits.R
#
# TAGS: @evaluation, @edge-cases, @JAGS
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-fit-edge-cases")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: JAGS_add_priors tests
# ============================================================================ #
test_that("JAGS_add_priors handles various prior types", {

  skip_if_not_installed("rjags")

  # Test with simple priors
  syntax_simple <- "model{}"
  priors_simple <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )

  result_simple <- JAGS_add_priors(syntax_simple, priors_simple)
  test_reference_text(result_simple, "JAGS_add_priors_simple.txt")

  # Test with truncated priors
  priors_truncated <- list(
    mu = prior("normal", list(0, 1), list(0, Inf))
  )

  result_truncated <- JAGS_add_priors(syntax_simple, priors_truncated)
  test_reference_text(result_truncated, "JAGS_add_priors_truncated.txt")

  # Test with point prior
  priors_point <- list(
    mu = prior("point", list(0))
  )

  result_point <- JAGS_add_priors(syntax_simple, priors_point)
  test_reference_text(result_point, "JAGS_add_priors_point.txt")

  # Test with factor priors
  priors_factor <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_factor[[1]], "levels") <- 3

  result_factor <- JAGS_add_priors(syntax_simple, priors_factor)
  test_reference_text(result_factor, "JAGS_add_priors_factor.txt")

  # Test with weightfunction priors
  priors_wf <- list(
    omega = prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))
  )

  result_wf <- JAGS_add_priors(syntax_simple, priors_wf)
  test_reference_text(result_wf, "JAGS_add_priors_weightfunction.txt")

})


# ============================================================================ #
# SECTION 2: JAGS_get_inits tests
# ============================================================================ #
test_that("JAGS_get_inits handles various prior types", {

  skip_if_not_installed("rjags")

  # Test with simple priors
  priors_simple <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )

  inits1 <- JAGS_get_inits(priors_simple, chains = 2, seed = 1)
  expect_equal(length(inits1), 2)
  expect_true("mu" %in% names(inits1[[1]]))
  expect_true("sigma" %in% names(inits1[[1]]))

  # Same seed should give same results
  inits2 <- JAGS_get_inits(priors_simple, chains = 2, seed = 1)
  expect_equal(inits1, inits2)

  # Different seeds should give different results
  inits3 <- JAGS_get_inits(priors_simple, chains = 2, seed = 123)
  expect_false(isTRUE(all.equal(inits1, inits3)))

  # Test with truncated priors
  priors_truncated <- list(
    mu = prior("normal", list(0, 1), list(0, Inf))
  )

  inits_truncated <- JAGS_get_inits(priors_truncated, chains = 2, seed = 1)
  expect_true(all(sapply(inits_truncated, function(i) i$mu >= 0)))

  # Test with point prior
  priors_point <- list(
    mu = prior("point", list(5))
  )

  inits_point <- JAGS_get_inits(priors_point, chains = 2, seed = 1)
  # Point priors should not generate inits (they're fixed)
  expect_true(!("mu" %in% names(inits_point[[1]])) || all(sapply(inits_point, function(i) i$mu == 5)))

  # Test with factor priors
  priors_factor <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_factor[[1]], "levels") <- 3

  inits_factor <- JAGS_get_inits(priors_factor, chains = 2, seed = 1)
  expect_true("p1" %in% names(inits_factor[[1]]))

})


# ============================================================================ #
# SECTION 3: JAGS_check_convergence tests
# ============================================================================ #
test_that("JAGS_check_convergence works with fitted models", {

  skip_if_not_installed("rjags")
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  prior_list <- attr(fit_simple, "prior_list")

  # Test convergence check with prior_list
  convergence <- JAGS_check_convergence(fit_simple, prior_list = prior_list)
  expect_true(is.logical(convergence) || is.list(convergence))

  # Test with NULL prior_list
  convergence_null <- JAGS_check_convergence(fit_simple, prior_list = NULL)
  expect_true(is.logical(convergence_null) || is.list(convergence_null))

})


# ============================================================================ #
# SECTION 4: JAGS_to_monitor tests
# ============================================================================ #
test_that("JAGS_to_monitor generates correct monitor strings", {

  skip_if_not_installed("rjags")

  # Test with simple priors
  priors_simple <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )

  monitor <- JAGS_to_monitor(priors_simple)
  test_reference_text(paste(sort(monitor), collapse = ","), "JAGS_to_monitor_simple.txt")

  # Test with point prior
  priors_with_point <- list(
    mu = prior("normal", list(0, 1)),
    fixed = prior("point", list(0))
  )

  monitor_point <- JAGS_to_monitor(priors_with_point)
  test_reference_text(paste(sort(monitor), collapse = ", "), "JAGS_to_monitor_point.txt")

  # Test with factor priors
  priors_factor <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_factor[[1]], "levels") <- 3

  monitor_factor <- JAGS_to_monitor(priors_factor)
  test_reference_text(paste(sort(monitor_factor), collapse = ","), "JAGS_to_monitor_factor.txt")

})


# ============================================================================ #
# SECTION 5: JAGS_fit attribute preservation
# ============================================================================ #
test_that("JAGS_fit preserves attributes", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # Check that prior_list attribute is preserved
  prior_list <- attr(fit_simple, "prior_list")
  expect_true(!is.null(prior_list))
  expect_true(is.list(prior_list))

  # Check class
  expect_true(inherits(fit_simple, "BayesTools_fit") || inherits(fit_simple, "runjags"))

})


# ============================================================================ #
# SECTION 6: runjags_estimates_table tests (diagnostics via summary-tables)
# ============================================================================ #
test_that("runjags_estimates_table works with fitted models", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # Test basic estimates table
  estimates_table <- runjags_estimates_table(fit_simple)
  test_reference_table(estimates_table, "runjags_estimates_simple.txt")

  # Test without specific parameters
  estimates_table_param <- runjags_estimates_table(fit_simple, remove_parameters = "m")
  test_reference_table(estimates_table_param, "runjags_estimates_param_m.txt")

})


# ============================================================================ #
# SECTION 7: JAGS_extend tests
# ============================================================================ #
test_that("JAGS_extend works correctly", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # Test extending a fitted model
  fit_extended <- JAGS_extend(
    fit_simple,
    autofit_control = list(
      max_Rhat = 1.05,
      min_ESS = 100,
      max_error = 0.01,
      max_SD_error = 0.05,
      max_time = list(time = 1, unit = "mins"),
      sample_extend = 100,
      restarts = 2,
      max_extend = 2
    ),
    silent = TRUE,
    seed = 1
  )

  # Test extending a fitted model
  fit_extended2 <- JAGS_extend(
    fit_simple,
    autofit_control = list(
      max_Rhat = 1.05,
      min_ESS = 100,
      max_error = 0.01,
      max_SD_error = 0.05,
      max_time = list(time = 1, unit = "mins"),
      sample_extend = 100,
      restarts = 2,
      max_extend = 2
    ),
    parallel = TRUE,
    cores = 2,
    silent = TRUE,
    seed = 1
  )

  # Check that the extended fit is still a BayesTools_fit

  expect_true(inherits(fit_extended, "BayesTools_fit"))
  expect_true(inherits(fit_extended, "runjags"))
  expect_true(inherits(fit_extended2, "BayesTools_fit"))
  expect_true(inherits(fit_extended2, "runjags"))

  # Check that attributes are preserved
  expect_true(!is.null(attr(fit_extended, "prior_list")))
  expect_true(!is.null(attr(fit_extended, "model_syntax")))
  expect_true(!is.null(attr(fit_extended2, "prior_list")))
  expect_true(!is.null(attr(fit_extended2, "model_syntax")))

  # Check that the extended fit has more samples
  original_samples  <- nrow(suppressWarnings(coda::as.mcmc(fit_simple)))
  extended_samples  <- nrow(suppressWarnings(coda::as.mcmc(fit_extended)))
  extended_samples2 <- nrow(suppressWarnings(coda::as.mcmc(fit_extended2)))
  expect_true(extended_samples  >= original_samples)
  expect_true(extended_samples2 >= original_samples)

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
# SECTION 8: .check_JAGS_syntax error handling
# ============================================================================ #
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


# ============================================================================ #
# SECTION 9: JAGS_fit with is_JASP mode
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


# ============================================================================ #
# SECTION 10: .JAGS_prior.mixture with PEESE prior
# ============================================================================ #
test_that("JAGS_add_priors handles mixture with PEESE prior", {

  skip_if_not_installed("rjags")

  # Create a bias mixture with PEESE prior
  bias_mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PEESE("normal", list(0, 1), prior_weights = 1)
  ))

  priors_peese <- list(
    bias = bias_mixture
  )

  result_peese <- JAGS_add_priors("model{}", priors_peese)
  test_reference_text(result_peese, "JAGS_add_priors_peese_mixture.txt")

})

test_that("JAGS_add_priors handles mixture with PET prior", {

  skip_if_not_installed("rjags")

  # Create a bias mixture with PET prior
  bias_mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))

  priors_pet <- list(
    bias = bias_mixture
  )

  result_pet <- JAGS_add_priors("model{}", priors_pet)
  test_reference_text(result_pet, "JAGS_add_priors_pet_mixture.txt")

})


# ============================================================================ #
# SECTION 11: Additional coverage tests for uncovered code paths
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


test_that("JAGS_check_and_list_autofit_settings validates all parameters", {

  # Valid settings
  valid_settings <- list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 1, unit = "mins"),
    sample_extend = 100,
    restarts = 3,
    max_extend = 10
  )
  expect_silent(JAGS_check_and_list_autofit_settings(valid_settings))

  # max_time without names - should auto-assign
  unnamed_time <- list(
    max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05,
    max_time = list(1, "mins"), sample_extend = 100
  )
  expect_silent(JAGS_check_and_list_autofit_settings(unnamed_time))

})


test_that("JAGS_add_priors handles spike_and_slab priors", {

  skip_if_not_installed("rjags")

  priors_sas <- list(
    mu = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    )
  )

  result <- JAGS_add_priors("model{}", priors_sas)
  expect_true(grepl("mu_variable", result))
  expect_true(grepl("mu_inclusion", result))
  expect_true(grepl("mu_indicator", result))

  # Test inits
  inits <- JAGS_get_inits(priors_sas, chains = 2, seed = 1)
  expect_true("mu_variable" %in% names(inits[[1]]) || "mu_inclusion" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_sas)
  expect_true("mu_indicator" %in% monitor)

})


test_that("JAGS_add_priors handles standard prior_mixture (non-bias)", {

  skip_if_not_installed("rjags")

  # Standard mixture (not bias mixture)
  mix <- prior_mixture(list(
    prior("normal", list(0, 0.5)),
    prior("normal", list(0, 1))
  ), is_null = c(TRUE, FALSE))

  priors_mix <- list(mu = mix)

  result <- JAGS_add_priors("model{}", priors_mix)
  expect_true(grepl("mu_indicator", result))
  expect_true(grepl("mu_component_1", result))
  expect_true(grepl("mu_component_2", result))

  # Test inits
  inits <- JAGS_get_inits(priors_mix, chains = 2, seed = 1)
  expect_true("mu_indicator" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_mix)
  expect_true("mu_indicator" %in% monitor)
  expect_true("mu" %in% monitor)

})


test_that("JAGS handles invgamma prior", {

  skip_if_not_installed("rjags")

  priors_inv <- list(tau = prior("invgamma", list(3, 2)))

  # Test syntax
  result <- JAGS_add_priors("model{}", priors_inv)
  expect_true(grepl("inv_tau", result))
  expect_true(grepl("dgamma", result))

  # Test inits
  inits <- JAGS_get_inits(priors_inv, chains = 2, seed = 1)
  expect_true("inv_tau" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_inv)
  expect_true("tau" %in% monitor)

})


test_that("JAGS handles weightfunction one.sided with alpha1/alpha2", {

  skip_if_not_installed("rjags")

  # One-sided with steps crossing 0.5 uses alpha1/alpha2 parametrization
  priors_wf2 <- list(omega = prior_weightfunction("one.sided", list(c(0.05, 0.60), c(1, 1), c(1, 1))))

  # Test syntax
  result <- JAGS_add_priors("model{}", priors_wf2)
  expect_true(grepl("eta1", result))
  expect_true(grepl("eta2", result))

  # Test inits
  inits <- JAGS_get_inits(priors_wf2, chains = 2, seed = 1)
  expect_true("eta1" %in% names(inits[[1]]))
  expect_true("eta2" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_wf2)
  expect_true("eta1" %in% monitor)
  expect_true("eta2" %in% monitor)

})


test_that("JAGS handles weightfunction fixed prior", {

  skip_if_not_installed("rjags")

  priors_wf_fixed <- list(omega = prior_weightfunction("one.sided.fixed", list(steps = c(0.05), omega = c(1, 0.5))))

  # Test syntax - fixed weightfunction has no eta parameters to sample
  result <- JAGS_add_priors("model{}", priors_wf_fixed)
  expect_true(grepl("omega", result))

  # Test inits - fixed weightfunction should return empty inits for eta
  inits <- JAGS_get_inits(priors_wf_fixed, chains = 2, seed = 1)
  # Should not have eta since it's fixed
  expect_true(!("eta" %in% names(inits[[1]])))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_wf_fixed)
  expect_true("omega" %in% monitor)

})


test_that("JAGS handles factor treatment/independent priors", {

  skip_if_not_installed("rjags")

  # Treatment contrast
  prior_treat <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior_treat, "levels") <- 3

  priors_treat <- list(fac = prior_treat)
  result_treat <- JAGS_add_priors("model{}", priors_treat)
  expect_true(grepl("fac\\[i\\]", result_treat))

  # Independent contrast
  prior_indep <- prior_factor("gamma", list(2, 1), contrast = "independent")
  attr(prior_indep, "levels") <- 2

  priors_indep <- list(fac = prior_indep)
  result_indep <- JAGS_add_priors("model{}", priors_indep)
  expect_true(grepl("dgamma", result_indep))

})


test_that("JAGS handles vector mt prior", {

  skip_if_not_installed("rjags")

  prior_mt <- prior("mt", list(location = 0, scale = 1, df = 5, K = 2))
  priors_mt <- list(p = prior_mt)

  # Test syntax
  result <- JAGS_add_priors("model{}", priors_mt)
  expect_true(grepl("prior_par_s_p", result))
  expect_true(grepl("prior_par_z_p", result))

  # Test inits
  inits <- JAGS_get_inits(priors_mt, chains = 2, seed = 1)
  expect_true("prior_par_s_p" %in% names(inits[[1]]))
  expect_true("prior_par_z_p" %in% names(inits[[1]]))

})


test_that("JAGS handles bias mixture with weightfunction", {

  skip_if_not_installed("rjags")

  bias_mix_wf <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one.sided", list(c(0.05), c(1, 1)), prior_weights = 1)
  ))

  priors_bias_wf <- list(bias = bias_mix_wf)

  result <- JAGS_add_priors("model{}", priors_bias_wf)
  expect_true(grepl("bias_indicator", result))
  expect_true(grepl("omega", result))
  expect_true(grepl("eta", result))

  # Test inits
  inits <- JAGS_get_inits(priors_bias_wf, chains = 2, seed = 1)
  expect_true("bias_indicator" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_bias_wf)
  expect_true("bias_indicator" %in% monitor)
  expect_true("omega" %in% monitor)

})
