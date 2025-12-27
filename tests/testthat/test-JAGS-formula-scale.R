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
    x_cont = rnorm(60, mean = 3, sd = 5),
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
  expect_equal(unname(result$data$mu_data_x_cont), as.numeric(df$x_cont))

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
  expect_equal(unname(result_scaled$data$mu_data_x_cont), as.numeric(scale(df$x_cont)))

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
  expect_equal(unname(result_not_scaled$data$mu_data_x_cont), as.numeric(df$x_cont))

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


  # Standardize one predictors
  result <- JAGS_formula(
    formula       = ~ x1 + x2,
    parameter     = "beta",
    data          = df,
    prior_list    = prior_list,
    formula_scale = list(x1 = FALSE, x2 = TRUE)
  )

  # Check that both predictors are standardized
  expect_length(result$formula_scale, 1)
  expect_true(!"beta_x1" %in% names(result$formula_scale))
  expect_true("beta_x2" %in% names(result$formula_scale))

  # Verify scale parameters
  expect_equal(result$formula_scale$beta_x2$mean, -5, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$sd, 2, tolerance = 0.5)
  expect_equal(unname(result$data$beta_data_x1), as.numeric(df$x1))
  expect_equal(unname(result$data$beta_data_x2), as.numeric(scale(df$x2)))
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

  skip_if_no_fits()

  # Load pre-fitted models
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto   <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))

  # Check that automatic scaling has formula_scale attribute
  expect_true(!is.null(attr(fit_auto, "formula_scale")))
  expect_true("mu_x_cont1" %in% names(attr(fit_auto, "formula_scale")))
  expect_true("mu_x_cont2" %in% names(attr(fit_auto, "formula_scale")))

  # Check that manual scaling has the scale info stored
  expect_true(!is.null(attr(fit_manual, "manual_scale")))

  # Compare scaling parameters
  # The automatic and manual scaling should have stored the same mean/sd
  manual_scale <- attr(fit_manual, "manual_scale")
  auto_scale   <- attr(fit_auto, "formula_scale")

  expect_equal(manual_scale$mu_x_cont1$mean, auto_scale$mu_x_cont1$mean, tolerance = 1e-10)
  expect_equal(manual_scale$mu_x_cont1$sd,   auto_scale$mu_x_cont1$sd, tolerance = 1e-10)
  expect_equal(manual_scale$mu_x_cont2$mean, auto_scale$mu_x_cont2$mean, tolerance = 1e-10)
  expect_equal(manual_scale$mu_x_cont2$sd,   auto_scale$mu_x_cont2$sd, tolerance = 1e-10)

  # Extract posterior samples
  posterior_manual <- as.matrix(fit_manual$mcmc[[1]])
  posterior_auto   <- as.matrix(fit_auto$mcmc[[1]])

  # The raw posterior samples should be very similar (both are on scaled space)
  # since both models were fit with the same seed and same scaled data

  # Compare means of main effects
  mean_manual_x1 <- mean(posterior_manual[, "mu_x_cont1"])
  mean_auto_x1   <- mean(posterior_auto[, "mu_x_cont1"])

  mean_manual_x2 <- mean(posterior_manual[, "mu_x_cont2"])
  mean_auto_x2   <- mean(posterior_auto[, "mu_x_cont2"])

  mean_manual_interaction <- mean(posterior_manual[, "mu_x_cont1__xXx__x_cont2"])
  mean_auto_interaction   <- mean(posterior_auto[, "mu_x_cont1__xXx__x_cont2"])

  # These should be very close since both use scaled data
  expect_equal(mean_manual_x1, mean_auto_x1)
  expect_equal(mean_manual_x2, mean_auto_x2)
  expect_equal(mean_manual_interaction, mean_auto_interaction)

  # Compare standard deviations
  sd_manual_x1 <- sd(posterior_manual[, "mu_x_cont1"])
  sd_auto_x1   <- sd(posterior_auto[, "mu_x_cont1"])

  sd_manual_x2 <- sd(posterior_manual[, "mu_x_cont2"])
  sd_auto_x2   <- sd(posterior_auto[, "mu_x_cont2"])

  sd_manual_interaction <- sd(posterior_manual[, "mu_x_cont1__xXx__x_cont2"])
  sd_auto_interaction   <- sd(posterior_auto[, "mu_x_cont1__xXx__x_cont2"])

  expect_equal(sd_manual_x1, sd_auto_x1)
  expect_equal(sd_manual_x2, sd_auto_x2)
  expect_equal(sd_manual_interaction, sd_auto_interaction)

  # Compare intercepts (these should also be similar)
  mean_manual_int <- mean(posterior_manual[, "mu_intercept"])
  mean_auto_int   <- mean(posterior_auto[, "mu_intercept"])

  expect_equal(mean_manual_int, mean_auto_int, tolerance = 0.05)
})

test_that("Downstream functions work with scaled models", {

  skip_if_no_fits()

  # Load pre-fitted models
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto   <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))

  expect_equal(JAGS_estimates_table(fit_manual), JAGS_estimates_table(fit_auto))
})

test_that("JAGS_evaluate_formula applies scaling correctly", {
  skip_if_not_installed("rjags")
  skip_if_no_fits()
  
  # Load pre-fitted models with scaling
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  
  # Create new data with same scale as original (unscaled)
  set.seed(3)
  new_data <- data.frame(
    x_cont1 = rnorm(10, mean = 1000, sd = 500),
    x_cont2 = rnorm(10, mean = 0.5, sd = 0.1)
  )
  
  # Get prior lists from fit attributes
  prior_list_auto <- attr(fit_auto, "prior_list")
  prior_list_manual <- attr(fit_manual, "prior_list")
  
  # For manual scaling, we need to manually scale the new data
  manual_scale <- attr(fit_manual, "manual_scale")
  new_data_manual <- new_data
  new_data_manual$x_cont1 <- (new_data$x_cont1 - manual_scale$mu_x_cont1$mean) / manual_scale$mu_x_cont1$sd
  new_data_manual$x_cont2 <- (new_data$x_cont2 - manual_scale$mu_x_cont2$mean) / manual_scale$mu_x_cont2$sd
  
  # For automatic scaling, JAGS_evaluate_formula should apply scaling automatically
  # (using the formula_scale attribute from fit_auto)
  
  # Evaluate formula on new data
  pred_manual <- JAGS_evaluate_formula(
    fit = fit_manual,
    formula = ~ x_cont1 * x_cont2,
    parameter = "mu",
    data = new_data_manual,
    prior_list = prior_list_manual
  )
  
  pred_auto <- JAGS_evaluate_formula(
    fit = fit_auto,
    formula = ~ x_cont1 * x_cont2,
    parameter = "mu",
    data = new_data,  # Note: passing unscaled data
    prior_list = prior_list_auto
  )
  
  # The predictions should be very similar (both models use same scaled data internally)
  expect_equal(dim(pred_manual), dim(pred_auto))
  expect_equal(dim(pred_manual), c(10, 1000))  # 10 data points, 1000 posterior samples
  
  # Compare mean predictions across all data points
  mean_pred_manual <- rowMeans(pred_manual)
  mean_pred_auto <- rowMeans(pred_auto)
  
  # Should be very close since both use same scaling
  expect_equal(mean_pred_manual, mean_pred_auto, tolerance = 0.1)
  
  # Also check that without scaling, predictions would be different
  # (this verifies that scaling is actually being applied)
  pred_auto_no_scale <- JAGS_evaluate_formula(
    fit = fit_manual,  # Use manual fit which doesn't have formula_scale attribute
    formula = ~ x_cont1 * x_cont2,
    parameter = "mu",
    data = new_data,  # Unscaled data
    prior_list = prior_list_manual
  )
  
  mean_pred_no_scale <- rowMeans(pred_auto_no_scale)
  
  # These should be very different from the correctly scaled predictions
  expect_true(max(abs(mean_pred_no_scale - mean_pred_auto)) > 1)
})

