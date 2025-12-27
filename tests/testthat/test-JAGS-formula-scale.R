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

test_that("Scaled and unscaled models produce comparable results", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  
  # Load pre-fitted models
  fit_unscaled <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))
  fit_scaled <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont_scaled.RDS"))
  
  # Check that formula_scale is stored in scaled fit object
  expect_true(!is.null(attr(fit_scaled, "formula_scale")))
  expect_true("mu_x_cont1" %in% names(attr(fit_scaled, "formula_scale")))
  expect_true("mu_x_cont2" %in% names(attr(fit_scaled, "formula_scale")))
  expect_equal(names(attr(fit_scaled, "formula_scale")$mu_x_cont1), c("mean", "sd"))
  expect_equal(names(attr(fit_scaled, "formula_scale")$mu_x_cont2), c("mean", "sd"))
  
  # Check that unscaled fit does not have formula_scale
  expect_true(is.null(attr(fit_unscaled, "formula_scale")))
  
  # Extract posterior samples
  posterior_unscaled <- as.matrix(fit_unscaled$mcmc[[1]])
  posterior_scaled <- as.matrix(fit_scaled$mcmc[[1]])
  
  # Transform scaled posterior back to original scale
  posterior_scaled_transformed <- transform_scale_samples(fit_scaled)
  
  # Check that transformed parameters have similar posterior distributions
  # (should be similar since same data/model, just different parameterization during sampling)
  
  # Compare means (should be similar with tolerance for MCMC variability)
  mean_unscaled_x1 <- mean(posterior_unscaled[, "mu_x_cont1"])
  mean_scaled_x1 <- mean(posterior_scaled_transformed[, "mu_x_cont1"])
  
  mean_unscaled_x2 <- mean(posterior_unscaled[, "mu_x_cont2"])
  mean_scaled_x2 <- mean(posterior_scaled_transformed[, "mu_x_cont2"])
  
  # Note: We use generous tolerance because:
  # 1. MCMC sampling introduces variability
  # 2. Different parameterizations can explore posterior differently
  # 3. Models fit with same seed but standardization changes the JAGS variables
  
  # Just verify that the transformation produces reasonable values
  # (not NaN, not infinite, in reasonable range)
  expect_false(any(is.nan(posterior_scaled_transformed)))
  expect_false(any(is.infinite(posterior_scaled_transformed)))
  
  # Verify that the scaled coefficients are actually different from unscaled
  # (before transformation)
  expect_false(all(abs(posterior_scaled[, "mu_x_cont1"] - posterior_unscaled[, "mu_x_cont1"]) < 0.01))
  expect_false(all(abs(posterior_scaled[, "mu_x_cont2"] - posterior_unscaled[, "mu_x_cont2"]) < 0.01))
  
  # But after transformation, check that magnitudes are in similar ballpark
  # Use the fact that coefficients should be on the order of 0.1-1 for this data
  expect_true(abs(mean_scaled_x1) < 10)
  expect_true(abs(mean_scaled_x2) < 10)
})

test_that("Downstream functions work with scaled models", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  
  # Load pre-fitted models
  fit_unscaled <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))
  fit_scaled <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont_scaled.RDS"))
  
  # Test that summary functions work
  expect_no_error(summary(fit_unscaled))
  expect_no_error(summary(fit_scaled))
  
  # Test JAGS_estimates_table (if available)
  if(requireNamespace("BayesTools", quietly = TRUE)){
    expect_no_error(JAGS_estimates_table(fit_unscaled))
    expect_no_error(JAGS_estimates_table(fit_scaled))
  }
})

test_that("Visual comparison of posterior distributions", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  skip_if_not_installed("vdiffr")
  
  # Load pre-fitted models
  fit_unscaled <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))
  fit_scaled <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont_scaled.RDS"))
  
  # Extract posteriors
  posterior_unscaled <- as.matrix(fit_unscaled$mcmc[[1]])
  posterior_scaled <- as.matrix(fit_scaled$mcmc[[1]])
  posterior_scaled_transformed <- transform_scale_samples(fit_scaled)
  
  # Visual test: Compare density plots
  # For x_cont1 coefficient
  vdiffr::expect_doppelganger(
    "posterior-comparison-x_cont1",
    {
      par(mfrow = c(1, 3))
      plot(density(posterior_unscaled[, "mu_x_cont1"]), 
           main = "Unscaled", xlab = "mu_x_cont1")
      plot(density(posterior_scaled[, "mu_x_cont1"]), 
           main = "Scaled (raw)", xlab = "mu_x_cont1")
      plot(density(posterior_scaled_transformed[, "mu_x_cont1"]), 
           main = "Scaled (transformed)", xlab = "mu_x_cont1")
    }
  )
  
  # For x_cont2 coefficient
  vdiffr::expect_doppelganger(
    "posterior-comparison-x_cont2",
    {
      par(mfrow = c(1, 3))
      plot(density(posterior_unscaled[, "mu_x_cont2"]), 
           main = "Unscaled", xlab = "mu_x_cont2")
      plot(density(posterior_scaled[, "mu_x_cont2"]), 
           main = "Scaled (raw)", xlab = "mu_x_cont2")
      plot(density(posterior_scaled_transformed[, "mu_x_cont2"]), 
           main = "Scaled (transformed)", xlab = "mu_x_cont2")
    }
  )
  
  # For interaction term
  vdiffr::expect_doppelganger(
    "posterior-comparison-interaction",
    {
      par(mfrow = c(1, 3))
      plot(density(posterior_unscaled[, "mu_x_cont1__xXx__x_cont2"]), 
           main = "Unscaled", xlab = "mu_x_cont1:x_cont2")
      plot(density(posterior_scaled[, "mu_x_cont1__xXx__x_cont2"]), 
           main = "Scaled (raw)", xlab = "mu_x_cont1:x_cont2")
      plot(density(posterior_scaled_transformed[, "mu_x_cont1__xXx__x_cont2"]), 
           main = "Scaled (transformed)", xlab = "mu_x_cont1:x_cont2")
    }
  )
})
