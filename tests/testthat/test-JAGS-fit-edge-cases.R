context("JAGS fit edge cases and comprehensive tests")

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

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # Test basic estimates table
  estimates_table <- runjags_estimates_table(fit_simple)
  test_reference_table(estimates_table, "runjags_estimates_simple.txt")

  # Test without specific parameters
  estimates_table_param <- runjags_estimates_table(fit_simple, remove_parameters = "m")
  test_reference_table(estimates_table_param, "runjags_estimates_param_m.txt")

})
