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

  ### Standardize both predictors
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
  expect_equal(result$formula_scale$beta_x1$sd,    3, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$mean, -5, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$sd,    2, tolerance = 0.5)

  ### Standardize both predictors (lazily)
  result <- JAGS_formula(
    formula       = ~ x1 + x2,
    parameter     = "beta",
    data          = df,
    prior_list    = prior_list,
    formula_scale = TRUE
  )

  # Check that both predictors are standardized
  expect_length(result$formula_scale, 2)
  expect_true("beta_x1" %in% names(result$formula_scale))
  expect_true("beta_x2" %in% names(result$formula_scale))

  # Verify scale parameters
  expect_equal(result$formula_scale$beta_x1$mean, 10, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x1$sd,    3, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$mean, -5, tolerance = 0.5)
  expect_equal(result$formula_scale$beta_x2$sd,    2, tolerance = 0.5)

  ### Standardize one predictors
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
  expected_intercept <- posterior[, "mu_intercept"] - (posterior[, "mu_x_cont"] / 2 * 5)
  expect_equal(posterior_original[, "mu_intercept"], expected_intercept)
})

test_that("transform_scale_samples handles interaction terms correctly", {

  # Create mock posterior samples with interaction
  set.seed(4)
  n_samples <- 100

  # Simulated standardized coefficients
  posterior <- matrix(
    c(
      rnorm(n_samples, mean = 1.0, sd = 0.1),   # mu_intercept
      rnorm(n_samples, mean = 0.3, sd = 0.05),  # mu_x1 (standardized)
      rnorm(n_samples, mean = 0.2, sd = 0.05),  # mu_x2 (standardized)
      rnorm(n_samples, mean = 0.1, sd = 0.02)   # mu_x1__xXx__x2 (interaction)
    ),
    nrow = n_samples,
    ncol = 4
  )
  colnames(posterior) <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")

  # Scale information
  formula_scale <- list(
    mu_x1 = list(mean = 5, sd = 2),
    mu_x2 = list(mean = 10, sd = 4)
  )

  # Transform back to original scale
  posterior_original <- transform_scale_samples(posterior, formula_scale)

  # The interaction coefficient should be divided by (sd_x1 * sd_x2) = 2 * 4 = 8
  expect_equal(
    posterior_original[, "mu_x1__xXx__x2"],
    posterior[, "mu_x1__xXx__x2"] / (2 * 4),
    tolerance = 1e-10
  )

  # The main effect x1 should be: beta_x1_orig = beta_x1_z/sd_x1 - beta_int_orig * mean_x2
  beta_int_orig <- posterior[, "mu_x1__xXx__x2"] / 8
  beta_x1_z_div_sd <- posterior[, "mu_x1"] / 2
  expected_beta_x1 <- beta_x1_z_div_sd - beta_int_orig * 10
  expect_equal(posterior_original[, "mu_x1"], expected_beta_x1, tolerance = 1e-10)

  # The main effect x2 should be: beta_x2_orig = beta_x2_z/sd_x2 - beta_int_orig * mean_x1
  beta_x2_z_div_sd <- posterior[, "mu_x2"] / 4
  expected_beta_x2 <- beta_x2_z_div_sd - beta_int_orig * 5
  expect_equal(posterior_original[, "mu_x2"], expected_beta_x2, tolerance = 1e-10)

  # The intercept should be:
  # alpha_orig = alpha_z - (beta_x1_z/sd_x1)*mean_x1 - (beta_x2_z/sd_x2)*mean_x2 + beta_int_orig*mean_x1*mean_x2
  # Note: uses beta_z/sd (intermediate values), not beta_orig (interaction-adjusted values)
  expected_intercept <- posterior[, "mu_intercept"] -
    beta_x1_z_div_sd * 5 - beta_x2_z_div_sd * 10 + beta_int_orig * 5 * 10
  expect_equal(posterior_original[, "mu_intercept"], expected_intercept, tolerance = 1e-10)
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

test_that("Marginal likelihoods match for manual and automatic scaling", {

  skip_if_no_fits()

  # Load pre-fitted marginal likelihoods
  marglik_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled_marglik.RDS"))
  marglik_auto   <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled_marglik.RDS"))

  # The log marginal likelihoods should be very similar
  # (both models use same scaled data internally)
  expect_equal(marglik_manual$logml, marglik_auto$logml, tolerance = 0.1)
})

test_that("JAGS_evaluate_formula applies scaling correctly", {

  skip_if_no_fits()

  # Load pre-fitted models with scaling
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))
  fit_auto   <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))

  # Create new data with same scale as original (unscaled)
  set.seed(3)
  new_data <- data.frame(
    x_cont1 = rnorm(10, mean = 1000, sd = 1000),
    x_cont2 = rnorm(10, mean = 0.5,  sd = 0.1)
  )

  # Get prior lists from fit attributes
  prior_list_auto   <- attr(fit_auto, "prior_list")
  prior_list_manual <- attr(fit_manual, "prior_list")

  # For manual scaling, we need to manually scale the new data
  manual_scale <- attr(fit_manual, "manual_scale")
  new_data_manual         <- new_data
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

  # The predictions should be very similar
  # (both models use same scaled data internally, and seed)
  expect_equal(pred_manual, pred_auto)

  # Also check that without scaling, predictions would be different
  # (this verifies that scaling is actually being applied)
  pred_auto_no_scale <- JAGS_evaluate_formula(
    fit = fit_manual,  # Use manual fit which doesn't have formula_scale attribute
    formula = ~ x_cont1 * x_cont2,
    parameter = "mu",
    data = new_data,  # Unscaled data
    prior_list = prior_list_manual
  )

  # These should be very different from the correctly scaled predictions
  expect_true(any(rowMeans(pred_manual) - rowMeans(pred_auto_no_scale) > 1))
})

test_that("runjags_estimates_table with transform_scaled unscales coefficients", {
  # TODO: something is wrong here with the intercept handling
  skip_if_no_fits()

  # Load pre-fitted model with automatic scaling
  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))

  # Get formula_scale attribute
  formula_scale <- attr(fit_auto, "formula_scale")
  expect_true(!is.null(formula_scale))

  # Get estimates without unscaling
  estimates_scaled <- JAGS_estimates_table(fit_auto, transform_scaled = FALSE)

  # Get estimates with unscaling
  estimates_unscaled <- JAGS_estimates_table(fit_auto, transform_scaled = TRUE)

  # The scaled coefficient for x_cont1 should be divided by sd
  # to get the unscaled coefficient
  sd_x_cont1   <- formula_scale$mu_x_cont1$sd
  sd_x_cont2   <- formula_scale$mu_x_cont2$sd
  mean_x_cont1 <- formula_scale$mu_x_cont1$mean
  mean_x_cont2 <- formula_scale$mu_x_cont2$mean

  # Check that the interaction term is correctly unscaled (divided by product of SDs)
  scaled_coef_int   <- estimates_scaled["(mu) x_cont1:x_cont2", "Mean"]
  unscaled_coef_int <- estimates_unscaled["(mu) x_cont1:x_cont2", "Mean"]
  expect_equal(unscaled_coef_int, scaled_coef_int / (sd_x_cont1 * sd_x_cont2), tolerance = 1e-10)

  # The main effects are adjusted for interaction contributions
  # beta_x1_orig = beta_x1_z/sd_x1 - beta_int_orig * mean_x2
  scaled_coef_x1 <- estimates_scaled["(mu) x_cont1", "Mean"]
  expected_x1    <- scaled_coef_x1 / sd_x_cont1 - unscaled_coef_int * mean_x_cont2
  expect_equal(estimates_unscaled["(mu) x_cont1", "Mean"], expected_x1, tolerance = 1e-10)

  # beta_x2_orig = beta_x2_z/sd_x2 - beta_int_orig * mean_x1
  scaled_coef_x2 <- estimates_scaled["(mu) x_cont2", "Mean"]
  expected_x2    <- scaled_coef_x2 / sd_x_cont2 - unscaled_coef_int * mean_x_cont1
  expect_equal(estimates_unscaled["(mu) x_cont2", "Mean"], expected_x2, tolerance = 1e-10)

  # The intercept should be adjusted
  # alpha_orig = alpha_z - beta_x1_orig*mean_x1 - beta_x2_orig*mean_x2 - beta_int_orig*mean_x1*mean_x2
  scaled_intercept   <- estimates_scaled["(mu) intercept", "Mean"]
  expected_intercept <- scaled_intercept - expected_x1 * mean_x_cont1 - expected_x2 * mean_x_cont2 -
    unscaled_coef_int * mean_x_cont1 * mean_x_cont2
  expect_equal(estimates_unscaled["(mu) intercept", "Mean"], expected_intercept, tolerance = 1e-10)
})

test_that("runjags_estimates_table transform_scaled with return_samples works", {

  skip_if_no_fits()

  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  formula_scale <- attr(fit_auto, "formula_scale")

  # Get samples without unscaling
  samples_scaled <- JAGS_estimates_table(fit_auto, transform_scaled = FALSE, return_samples = TRUE)

  # Get samples with unscaling
  samples_unscaled <- JAGS_estimates_table(fit_auto, transform_scaled = TRUE, return_samples = TRUE)

  # Check that x_cont1 samples are correctly unscaled
  sd_x_cont1 <- formula_scale$mu_x_cont1$sd
  expect_equal(
    samples_unscaled[, "(mu) x_cont1"],
    samples_scaled[, "(mu) x_cont1"] / sd_x_cont1,
    tolerance = 1e-10
  )
})

test_that("ensemble_estimates_table with transform_scaled unscales coefficients", {

  skip_if_no_fits()
  skip_if_not_installed("bridgesampling")

  # Load pre-fitted models
  fit_auto    <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  marglik_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled_marglik.RDS"))

  formula_scale <- attr(fit_auto, "formula_scale")

  # Create a simple model list for mix_posteriors
  model_list <- list(
    list(
      fit           = fit_auto,
      marglik       = marglik_auto,
      prior_weights = 1
    )
  )

  # Get mixed posteriors
  mixed_posteriors <- mix_posteriors(
    model_list  = model_list,
    parameters  = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
    is_null_list = list(
      mu_intercept = 1,
      mu_x_cont1   = 1,
      mu_x_cont2   = 1
    ),
    seed = 1
  )

  # Get estimates without unscaling
  estimates_scaled <- ensemble_estimates_table(
    samples          = mixed_posteriors,
    parameters       = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
    transform_scaled = FALSE
  )

  # Get estimates with unscaling
  estimates_unscaled <- ensemble_estimates_table(
    samples          = mixed_posteriors,
    parameters       = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
    transform_scaled = TRUE,
    formula_scale    = formula_scale
  )

  # The scaled coefficient for x_cont1 should be divided by sd
  sd_x_cont1 <- formula_scale$mu_x_cont1$sd

  scaled_coef_x1   <- estimates_scaled["mu_x_cont1", "Mean"]
  unscaled_coef_x1 <- estimates_unscaled["mu_x_cont1", "Mean"]

  expect_equal(unscaled_coef_x1, scaled_coef_x1 / sd_x_cont1, tolerance = 1e-10)
})

test_that("transform_scaled = FALSE is the default behavior", {

  skip_if_no_fits()

  fit_auto <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))

  # Default behavior should be no unscaling
  estimates_default <- JAGS_estimates_table(fit_auto)
  estimates_false   <- JAGS_estimates_table(fit_auto, transform_scaled = FALSE)

  expect_equal(estimates_default, estimates_false)
})

test_that("transform_scaled has no effect when formula_scale is NULL", {

  skip_if_no_fits()

  # Load model without automatic scaling (manual scaling doesn't have formula_scale attr)
  fit_manual <- readRDS(file.path(temp_fits_dir, "fit_formula_manual_scaled.RDS"))

  # transform_scaled = TRUE should have no effect when formula_scale is NULL
  estimates_false <- JAGS_estimates_table(fit_manual, transform_scaled = FALSE)
  estimates_true  <- JAGS_estimates_table(fit_manual, transform_scaled = TRUE)

  expect_equal(estimates_false, estimates_true)
})


