# ============================================================================ #
# TEST FILE: JAGS Formula Standardization
# ============================================================================ #
#
# PURPOSE:
#   Tests for automatic standardization of continuous predictors in JAGS_formula
#
# DEPENDENCIES:
#   - common-functions.R: Test helpers
#
# SKIP CONDITIONS:
#   - None (pure R tests, no JAGS fitting required)
#
# TAGS: @formula, @standardization, @fast
# ============================================================================ #

# Load common test helpers
source(testthat::test_path("common-functions.R"))

test_that("JAGS_formula accepts and validates formula_scale parameter", {
  
  # Setup test data
  set.seed(1)
  df <- data.frame(
    y      = rnorm(60),
    x_cont = rnorm(60),
    x_fac  = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  
  # Test 1: formula_scale = NULL (no standardization)
  prior_list <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x_cont"    = prior("normal", list(0, 1)),
    "x_fac"     = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  
  result <- JAGS_formula(
    formula       = ~ x_cont + x_fac,
    parameter     = "mu",
    data          = df,
    prior_list    = prior_list,
    formula_scale = NULL
  )
  
  expect_false("formula_scale" %in% names(result))
  
  # Test 2: formula_scale with standardization
  result_scaled <- JAGS_formula(
    formula       = ~ x_cont + x_fac,
    parameter     = "mu",
    data          = df,
    prior_list    = prior_list,
    formula_scale = list(x_cont = TRUE)
  )
  
  expect_true("formula_scale" %in% names(result_scaled))
  expect_true("mu_x_cont" %in% names(result_scaled$formula_scale))
  expect_equal(names(result_scaled$formula_scale$mu_x_cont), c("mean", "sd"))
  
  # Test 3: Check that scale info is correct
  original_mean <- mean(df$x_cont)
  original_sd   <- sd(df$x_cont)
  
  expect_equal(result_scaled$formula_scale$mu_x_cont$mean, original_mean)
  expect_equal(result_scaled$formula_scale$mu_x_cont$sd, original_sd)
  
  # Test 4: formula_scale with FALSE should not standardize
  result_not_scaled <- JAGS_formula(
    formula       = ~ x_cont + x_fac,
    parameter     = "mu",
    data          = df,
    prior_list    = prior_list,
    formula_scale = list(x_cont = FALSE)
  )
  
  expect_false("formula_scale" %in% names(result_not_scaled))
  
})

test_that("JAGS_formula standardization preserves data correctly", {
  
  set.seed(2)
  df <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 3),
    x2 = rnorm(50, mean = -5, sd = 2)
  )
  
  prior_list <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x1"        = prior("normal", list(0, 1)),
    "x2"        = prior("normal", list(0, 1))
  )
  
  # Standardize both predictors
  result <- JAGS_formula(
    formula       = ~ x1 + x2,
    parameter     = "beta",
    data          = df,
    prior_list    = prior_list,
    formula_scale = list(x1 = TRUE, x2 = TRUE)
  )
  
  # Check that both predictors are standardized
  expect_length(result$formula_scale, 2)
  expect_true("beta_x1" %in% names(result$formula_scale))
  expect_true("beta_x2" %in% names(result$formula_scale))
  
  # Verify scale parameters
  expect_equal(result$formula_scale$beta_x1$mean, 10, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x1$sd, 3, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$mean, -5, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$sd, 2, tolerance = 0.5)
})

test_that("transform_scale_samples transforms coefficients correctly", {
  
  # Create mock posterior samples
  set.seed(3)
  n_samples <- 100
  
  # Simulated standardized coefficients
  posterior <- matrix(
    c(
      rnorm(n_samples, mean = 0.5, sd = 0.1),  # mu_intercept
      rnorm(n_samples, mean = 0.3, sd = 0.05), # mu_x_cont (standardized)
      rnorm(n_samples, mean = 0.2, sd = 0.05)  # mu_x_fac (not standardized)
    ),
    nrow = n_samples,
    ncol = 3
  )
  colnames(posterior) <- c("mu_intercept", "mu_x_cont", "mu_x_fac")
  
  # Scale information (x_cont was standardized with mean=5, sd=2)
  formula_scale <- list(
    mu_x_cont = list(mean = 5, sd = 2)
  )
  
  # Transform back to original scale
  posterior_original <- transform_scale_samples(posterior, formula_scale)
  
  # Check that x_cont coefficient is rescaled (divided by sd)
  expect_equal(posterior_original[, "mu_x_cont"], posterior[, "mu_x_cont"] / 2)
  
  # Check that x_fac is unchanged (not in formula_scale)
  expect_equal(posterior_original[, "mu_x_fac"], posterior[, "mu_x_fac"])
  
  # Check that intercept is adjusted
  # intercept_original = intercept_std - beta_original * mean
  # where beta_original = beta_std / sd (already done above)
  expected_intercept <- posterior[, "mu_intercept"] - 
    (posterior[, "mu_x_cont"] / 2 * 5)
  expect_equal(posterior_original[, "mu_intercept"], expected_intercept)
})

test_that("Manual and automatic scaling produce equivalent results", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  
  # Load pre-fitted models
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  
  # Check that automatic scaling has formula_scale attribute
  expect_true(!is.null(attr(fit_auto, "formula_scale")))
  expect_true("mu_x_cont1" %in% names(attr(fit_auto, "formula_scale")))
  expect_true("mu_x_cont2" %in% names(attr(fit_auto, "formula_scale")))
  
  # Check that manual scaling has the scale info stored
  expect_true(!is.null(attr(fit_manual, "manual_scale")))
  
  # Compare scaling parameters
  # The automatic and manual scaling should have stored the same mean/sd
  manual_scale <- attr(fit_manual, "manual_scale")
  auto_scale <- attr(fit_auto, "formula_scale")
  
  expect_equal(manual_scale$mu_x_cont1$mean, auto_scale$mu_x_cont1$mean, tolerance = 1e-10)
  expect_equal(manual_scale$mu_x_cont1$sd, auto_scale$mu_x_cont1$sd, tolerance = 1e-10)
  expect_equal(manual_scale$mu_x_cont2$mean, auto_scale$mu_x_cont2$mean, tolerance = 1e-10)
  expect_equal(manual_scale$mu_x_cont2$sd, auto_scale$mu_x_cont2$sd, tolerance = 1e-10)
  
  # Extract posterior samples
  posterior_manual <- as.matrix(fit_manual$mcmc[[1]])
  posterior_auto <- as.matrix(fit_auto$mcmc[[1]])
  
  # The raw posterior samples should be very similar (both are on scaled space)
  # since both models were fit with the same seed and same scaled data
  
  # Compare means of main effects
  mean_manual_x1 <- mean(posterior_manual[, "mu_x_cont1"])
  mean_auto_x1 <- mean(posterior_auto[, "mu_x_cont1"])
  
  mean_manual_x2 <- mean(posterior_manual[, "mu_x_cont2"])
  mean_auto_x2 <- mean(posterior_auto[, "mu_x_cont2"])
  
  mean_manual_interaction <- mean(posterior_manual[, "mu_x_cont1__xXx__x_cont2"])
  mean_auto_interaction <- mean(posterior_auto[, "mu_x_cont1__xXx__x_cont2"])
  
  # These should be very close since both use scaled data
  expect_equal(mean_manual_x1, mean_auto_x1, tolerance = 0.05)
  expect_equal(mean_manual_x2, mean_auto_x2, tolerance = 0.05)
  expect_equal(mean_manual_interaction, mean_auto_interaction, tolerance = 0.05)
  
  # Compare standard deviations
  sd_manual_x1 <- sd(posterior_manual[, "mu_x_cont1"])
  sd_auto_x1 <- sd(posterior_auto[, "mu_x_cont1"])
  
  sd_manual_x2 <- sd(posterior_manual[, "mu_x_cont2"])
  sd_auto_x2 <- sd(posterior_auto[, "mu_x_cont2"])
  
  sd_manual_interaction <- sd(posterior_manual[, "mu_x_cont1__xXx__x_cont2"])
  sd_auto_interaction <- sd(posterior_auto[, "mu_x_cont1__xXx__x_cont2"])
  
  expect_equal(sd_manual_x1, sd_auto_x1, tolerance = 0.05)
  expect_equal(sd_manual_x2, sd_auto_x2, tolerance = 0.05)
  expect_equal(sd_manual_interaction, sd_auto_interaction, tolerance = 0.05)
  
  # Compare intercepts (these should also be similar)
  mean_manual_int <- mean(posterior_manual[, "mu_intercept"])
  mean_auto_int <- mean(posterior_auto[, "mu_intercept"])
  
  expect_equal(mean_manual_int, mean_auto_int, tolerance = 0.05)
})

test_that("Downstream functions work with scaled models", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  
  # Load pre-fitted models
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  
  # Test that summary functions work
  expect_no_error(summary(fit_manual))
  expect_no_error(summary(fit_auto))
  
  # Test JAGS_estimates_table (if available)
  if(requireNamespace("BayesTools", quietly = TRUE)){
    expect_no_error(JAGS_estimates_table(fit_manual))
    expect_no_error(JAGS_estimates_table(fit_auto))
  }
})

test_that("Visual comparison of manual vs automatic scaling", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  skip_if_not_installed("vdiffr")
  
  # Load pre-fitted models
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  
  # Extract posteriors (both are on scaled space)
  posterior_manual <- as.matrix(fit_manual$mcmc[[1]])
  posterior_auto <- as.matrix(fit_auto$mcmc[[1]])
  
  # Visual test: Compare density plots
  # For x_cont1 coefficient (scaled space)
  vdiffr::expect_doppelganger(
    "scaling-comparison-x_cont1",
    {
      par(mfrow = c(1, 2))
      plot(density(posterior_manual[, "mu_x_cont1"]), 
           main = "Manual Scaling", xlab = "mu_x_cont1 (scaled)")
      plot(density(posterior_auto[, "mu_x_cont1"]), 
           main = "Automatic Scaling", xlab = "mu_x_cont1 (scaled)")
    }
  )
  
  # For x_cont2 coefficient (scaled space)
  vdiffr::expect_doppelganger(
    "scaling-comparison-x_cont2",
    {
      par(mfrow = c(1, 2))
      plot(density(posterior_manual[, "mu_x_cont2"]), 
           main = "Manual Scaling", xlab = "mu_x_cont2 (scaled)")
      plot(density(posterior_auto[, "mu_x_cont2"]), 
           main = "Automatic Scaling", xlab = "mu_x_cont2 (scaled)")
    }
  )
  
  # For interaction term (scaled space)
  vdiffr::expect_doppelganger(
    "scaling-comparison-interaction",
    {
      par(mfrow = c(1, 2))
      plot(density(posterior_manual[, "mu_x_cont1__xXx__x_cont2"]), 
           main = "Manual Scaling", xlab = "mu_x_cont1:x_cont2 (scaled)")
      plot(density(posterior_auto[, "mu_x_cont1__xXx__x_cont2"]), 
           main = "Automatic Scaling", xlab = "mu_x_cont1:x_cont2 (scaled)")
    }
  )
  
  # For intercept (scaled space)
  vdiffr::expect_doppelganger(
    "scaling-comparison-intercept",
    {
      par(mfrow = c(1, 2))
      plot(density(posterior_manual[, "mu_intercept"]), 
           main = "Manual Scaling", xlab = "mu_intercept (scaled)")
      plot(density(posterior_auto[, "mu_intercept"]), 
           main = "Automatic Scaling", xlab = "mu_intercept (scaled)")
    }
  )
})
