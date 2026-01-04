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
  # Use nested structure keyed by parameter name
  formula_scale <- list(
    mu = list(
      mu_x_cont = list(mean = 5, sd = 2)
    )
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

  # Scale information - use nested structure keyed by parameter name
  formula_scale <- list(
    mu = list(
      mu_x1 = list(mean = 5, sd = 2),
      mu_x2 = list(mean = 10, sd = 4)
    )
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

  # Check that automatic scaling has formula_scale attribute with nested structure
  expect_true(!is.null(attr(fit_auto, "formula_scale")))
  expect_true("mu" %in% names(attr(fit_auto, "formula_scale")))
  expect_true("mu_x_cont1" %in% names(attr(fit_auto, "formula_scale")$mu))
  expect_true("mu_x_cont2" %in% names(attr(fit_auto, "formula_scale")$mu))

  # Check that manual scaling has the scale info stored
  expect_true(!is.null(attr(fit_manual, "manual_scale")))

  # Compare scaling parameters
  # The automatic and manual scaling should have stored the same mean/sd
  manual_scale <- attr(fit_manual, "manual_scale")
  auto_scale   <- attr(fit_auto, "formula_scale")$mu

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
  marglik_manual <- readRDS(file.path(temp_marglik_dir, "fit_formula_manual_scaled.RDS"))
  marglik_auto   <- readRDS(file.path(temp_marglik_dir, "fit_formula_auto_scaled.RDS"))

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
  # to get the unscaled coefficient (nested structure: formula_scale$mu$...)
  sd_x_cont1   <- formula_scale$mu$mu_x_cont1$sd
  sd_x_cont2   <- formula_scale$mu$mu_x_cont2$sd
  mean_x_cont1 <- formula_scale$mu$mu_x_cont1$mean
  mean_x_cont2 <- formula_scale$mu$mu_x_cont2$mean

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

  # For models with interactions, the transformation is more complex (nested structure)
  sd_x_cont1   <- formula_scale$mu$mu_x_cont1$sd
  sd_x_cont2   <- formula_scale$mu$mu_x_cont2$sd
  mean_x_cont1 <- formula_scale$mu$mu_x_cont1$mean
  mean_x_cont2 <- formula_scale$mu$mu_x_cont2$mean

  # First, compute the unscaled interaction coefficient
  unscaled_int <- samples_scaled[, "(mu) x_cont1:x_cont2"] / (sd_x_cont1 * sd_x_cont2)

  # Check that x_cont1 samples are correctly unscaled (with interaction adjustment)
  # beta_x1_orig = beta_x1_z/sd_x1 - beta_int_orig * mean_x2
  expected_x1 <- samples_scaled[, "(mu) x_cont1"] / sd_x_cont1 - unscaled_int * mean_x_cont2
  expect_equal(
    samples_unscaled[, "(mu) x_cont1"],
    expected_x1,
    tolerance = 1e-10
  )
})

test_that("ensemble_estimates_table with transform_scaled unscales coefficients", {

  skip_if_no_fits()
  skip_if_not_installed("bridgesampling")

  # Load pre-fitted models
  fit_auto     <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  marglik_auto <- readRDS(file.path(temp_marglik_dir, "fit_formula_auto_scaled.RDS"))

  formula_scale <- attr(fit_auto, "formula_scale")

  # Create a simple model list for mix_posteriors
  model_list <- list(
    list(
      fit           = fit_auto,
      marglik       = marglik_auto,
      prior_weights = 1
    )
  )

  # Get mixed posteriors - include interaction term for proper unscaling
  mixed_posteriors <- mix_posteriors(
    model_list   = model_list,
    parameters   = c("mu_intercept", "mu_x_cont1", "mu_x_cont2", "mu_x_cont1__xXx__x_cont2"),
    is_null_list = list(
      mu_intercept               = 1,
      mu_x_cont1                 = 1,
      mu_x_cont2                 = 1,
      "mu_x_cont1__xXx__x_cont2" = 1
    ),
    seed = 1
  )

  # Get estimates without unscaling
  estimates_scaled <- ensemble_estimates_table(
    samples          = mixed_posteriors,
    parameters       = c("mu_intercept", "mu_x_cont1", "mu_x_cont2", "mu_x_cont1__xXx__x_cont2"),
    transform_scaled = FALSE
  )

  # Get estimates with unscaling
  estimates_unscaled <- ensemble_estimates_table(
    samples          = mixed_posteriors,
    parameters       = c("mu_intercept", "mu_x_cont1", "mu_x_cont2", "mu_x_cont1__xXx__x_cont2"),
    transform_scaled = TRUE,
    formula_scale    = formula_scale
  )

  # For models with interactions, the transformation is more complex (nested structure)
  sd_x_cont1   <- formula_scale$mu$mu_x_cont1$sd
  sd_x_cont2   <- formula_scale$mu$mu_x_cont2$sd
  mean_x_cont1 <- formula_scale$mu$mu_x_cont1$mean
  mean_x_cont2 <- formula_scale$mu$mu_x_cont2$mean

  # Check that the interaction term is correctly unscaled (divided by product of SDs)
  scaled_coef_int   <- estimates_scaled["(mu) x_cont1:x_cont2", "Mean"]
  unscaled_coef_int <- estimates_unscaled["(mu) x_cont1:x_cont2", "Mean"]
  expect_equal(unscaled_coef_int, scaled_coef_int / (sd_x_cont1 * sd_x_cont2), tolerance = 1e-10)

  # The main effects are adjusted for interaction contributions
  # beta_x1_orig = beta_x1_z/sd_x1 - beta_int_orig * mean_x2
  scaled_coef_x1 <- estimates_scaled["mu_x_cont1", "Mean"]
  expected_x1    <- scaled_coef_x1 / sd_x_cont1 - unscaled_coef_int * mean_x_cont2
  expect_equal(estimates_unscaled["mu_x_cont1", "Mean"], expected_x1, tolerance = 1e-10)

  # beta_x2_orig = beta_x2_z/sd_x2 - beta_int_orig * mean_x1
  scaled_coef_x2 <- estimates_scaled["mu_x_cont2", "Mean"]
  expected_x2    <- scaled_coef_x2 / sd_x_cont2 - unscaled_coef_int * mean_x_cont1
  expect_equal(estimates_unscaled["mu_x_cont2", "Mean"], expected_x2, tolerance = 1e-10)
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


# ============================================================================ #
# DUAL PARAMETER REGRESSION WITH LOG(INTERCEPT) TESTS
# ============================================================================ #

test_that("Dual parameter model with log(intercept) has correct formula_scale structure", {

  skip_if_no_fits()

  # Load pre-fitted dual parameter regression model
  fit_dual <- readRDS(file.path(temp_fits_dir, "fit_dual_param_regression.RDS"))

  # Check that formula_scale attribute exists
  formula_scale <- attr(fit_dual, "formula_scale")
  expect_true(!is.null(formula_scale))


  # Check that both parameters have scaling info
  expect_true("mu" %in% names(formula_scale))
  expect_true("log_sigma" %in% names(formula_scale))
  
  # Check nested structure
  expect_true("mu_x_mu" %in% names(formula_scale$mu))
  expect_true("log_sigma_x_sigma" %in% names(formula_scale$log_sigma))

  # Verify scale info structure for mu parameter
  expect_equal(names(formula_scale$mu$mu_x_mu), c("mean", "sd"))
  expect_true(is.numeric(formula_scale$mu$mu_x_mu$mean))
  expect_true(is.numeric(formula_scale$mu$mu_x_mu$sd))

  # Verify scale info structure for log_sigma parameter
  expect_equal(names(formula_scale$log_sigma$log_sigma_x_sigma), c("mean", "sd"))
  expect_true(is.numeric(formula_scale$log_sigma$log_sigma_x_sigma$mean))
  expect_true(is.numeric(formula_scale$log_sigma$log_sigma_x_sigma$sd))
  
  # Verify log_intercept attribute is stored correctly
  # mu should NOT have log_intercept (or be FALSE)
  expect_false(isTRUE(attr(formula_scale$mu, "log_intercept")))
  # log_sigma SHOULD have log_intercept = TRUE
  expect_true(isTRUE(attr(formula_scale$log_sigma, "log_intercept")))

  # Verify the model has expected parameters
  param_names <- colnames(fit_dual$mcmc[[1]])
  expect_true("mu_intercept" %in% param_names)
  expect_true("mu_x_mu" %in% param_names)
  expect_true("log_sigma_intercept" %in% param_names)
  expect_true("log_sigma_x_sigma" %in% param_names)
})

test_that("transform_scale_samples works with dual parameter model", {

  skip_if_no_fits()

  # Load pre-fitted dual parameter regression model
  fit_dual <- readRDS(file.path(temp_fits_dir, "fit_dual_param_regression.RDS"))
  formula_scale <- attr(fit_dual, "formula_scale")

  # Extract posterior samples
  posterior <- as.matrix(fit_dual$mcmc[[1]])

  # Transform to original scale
  posterior_transformed <- transform_scale_samples(posterior, formula_scale)

  # Get scale parameters (nested structure)
  mu_scale        <- formula_scale$mu$mu_x_mu
  log_sigma_scale <- formula_scale$log_sigma$log_sigma_x_sigma

  # Check mu_x_mu coefficient is correctly unscaled (divided by sd)
  expected_mu_x_mu <- posterior[, "mu_x_mu"] / mu_scale$sd
  expect_equal(posterior_transformed[, "mu_x_mu"], expected_mu_x_mu, tolerance = 1e-10)

  # Check log_sigma_x_sigma coefficient is correctly unscaled (divided by sd)
  expected_log_sigma_x_sigma <- posterior[, "log_sigma_x_sigma"] / log_sigma_scale$sd
  expect_equal(posterior_transformed[, "log_sigma_x_sigma"], expected_log_sigma_x_sigma, tolerance = 1e-10)

  # Check mu intercept is adjusted: intercept_orig = intercept_z - beta_orig * mean
  expected_mu_intercept <- posterior[, "mu_intercept"] - expected_mu_x_mu * mu_scale$mean
  expect_equal(posterior_transformed[, "mu_intercept"], expected_mu_intercept, tolerance = 1e-10)

  # Check log_sigma intercept is adjusted with multiplicative transformation (due to log(intercept)):
  # intercept_orig = intercept_z * exp(-beta_orig * mean)
  expected_log_sigma_intercept <- posterior[, "log_sigma_intercept"] * exp(-expected_log_sigma_x_sigma * log_sigma_scale$mean)
  expect_equal(posterior_transformed[, "log_sigma_intercept"], expected_log_sigma_intercept, tolerance = 1e-10)
})

test_that("JAGS_estimates_table with transform_scaled works for dual parameter model", {

  skip_if_no_fits()

  # Load pre-fitted dual parameter regression model
  fit_dual <- readRDS(file.path(temp_fits_dir, "fit_dual_param_regression.RDS"))
  formula_scale <- attr(fit_dual, "formula_scale")

  # Get estimates without unscaling
  estimates_scaled <- JAGS_estimates_table(fit_dual, transform_scaled = FALSE)

  # Get estimates with unscaling
  estimates_unscaled <- JAGS_estimates_table(fit_dual, transform_scaled = TRUE)

  # Get scale parameters (nested structure)
  mu_sd             <- formula_scale$mu$mu_x_mu$sd
  mu_mean           <- formula_scale$mu$mu_x_mu$mean
  log_sigma_sd      <- formula_scale$log_sigma$log_sigma_x_sigma$sd
  log_sigma_mean    <- formula_scale$log_sigma$log_sigma_x_sigma$mean

  # Check mu_x_mu coefficient is correctly unscaled
  scaled_mu_coef   <- estimates_scaled["(mu) x_mu", "Mean"]
  unscaled_mu_coef <- estimates_unscaled["(mu) x_mu", "Mean"]
  expect_equal(unscaled_mu_coef, scaled_mu_coef / mu_sd, tolerance = 1e-10)

  # Check log_sigma_x_sigma coefficient is correctly unscaled
  scaled_log_sigma_coef   <- estimates_scaled["(log_sigma) x_sigma", "Mean"]
  unscaled_log_sigma_coef <- estimates_unscaled["(log_sigma) x_sigma", "Mean"]
  expect_equal(unscaled_log_sigma_coef, scaled_log_sigma_coef / log_sigma_sd, tolerance = 1e-10)

  # Check mu intercept is correctly adjusted
  scaled_mu_int   <- estimates_scaled["(mu) intercept", "Mean"]
  expected_mu_int <- scaled_mu_int - unscaled_mu_coef * mu_mean
  expect_equal(estimates_unscaled["(mu) intercept", "Mean"], expected_mu_int, tolerance = 1e-10)

  # Check log_sigma intercept is correctly adjusted with multiplicative transformation
  # Due to log(intercept): intercept_orig = intercept_z * exp(-beta_orig * mean)
  # For means, we can't use the simple relationship because E[X * exp(Y)] != E[X] * exp(E[Y])
  # Instead, verify that the unscaled intercept is close to the true value (0.5)
  # and that it differs from the scaled intercept (which would be biased)
  # Note: with transform_scaled=TRUE, the intercept is renamed to exp(intercept)
  unscaled_log_sigma_int <- estimates_unscaled["(log_sigma) exp(intercept)", "Mean"]
  
  # The unscaled intercept should be reasonably close to the true value of 0.5
  expect_true(abs(unscaled_log_sigma_int - 0.5) < 0.15)
  
  # The scaled intercept should NOT be close to 0.5 (it's on the wrong scale)
  scaled_log_sigma_int <- estimates_scaled["(log_sigma) intercept", "Mean"]
  expect_true(abs(scaled_log_sigma_int - 0.5) > abs(unscaled_log_sigma_int - 0.5))
})

test_that("JAGS_evaluate_formula applies scaling correctly for dual parameter model", {

  skip_if_no_fits()

  # Load pre-fitted dual parameter regression model
  fit_dual <- readRDS(file.path(temp_fits_dir, "fit_dual_param_regression.RDS"))
  formula_scale <- attr(fit_dual, "formula_scale")
  prior_list <- attr(fit_dual, "prior_list")

  # Create new data (on original unscaled scale)
  set.seed(123)
  new_data_mu <- data.frame(x_mu = rnorm(5, mean = 5, sd = 2))
  new_data_sigma <- data.frame(x_sigma = rnorm(5, mean = 3, sd = 1.5))

  # Evaluate mu formula (standard intercept)
  pred_mu <- JAGS_evaluate_formula(
    fit       = fit_dual,
    formula   = ~ x_mu,
    parameter = "mu",
    data      = new_data_mu,
    prior_list = prior_list
  )

  # Evaluate log_sigma formula (log intercept)
  formula_log_sigma <- ~ x_sigma
  attr(formula_log_sigma, "log(intercept)") <- TRUE

  pred_log_sigma <- JAGS_evaluate_formula(
    fit       = fit_dual,
    formula   = formula_log_sigma,
    parameter = "log_sigma",
    data      = new_data_sigma,
    prior_list = prior_list
  )

  # Basic sanity checks
  expect_equal(nrow(pred_mu), 5)
  expect_equal(nrow(pred_log_sigma), 5)

  # The predictions should be matrices with n_samples columns
  expect_true(ncol(pred_mu) > 1)
  expect_true(ncol(pred_log_sigma) > 1)

  # Verify manually: predictions should match manual calculation
  posterior <- as.matrix(coda::as.mcmc.list(fit_dual))
  mu_scale <- formula_scale$mu$mu_x_mu
  log_sigma_scale <- formula_scale$log_sigma$log_sigma_x_sigma

  # Scale the new data as the function should do internally
  x_mu_scaled <- (new_data_mu$x_mu - mu_scale$mean) / mu_scale$sd
  x_sigma_scaled <- (new_data_sigma$x_sigma - log_sigma_scale$mean) / log_sigma_scale$sd

  # For first observation in new_data_mu
  # mu[i] = intercept + x_mu * x_mu_scaled[i]
  expected_mu_1 <- posterior[, "mu_intercept"] + posterior[, "mu_x_mu"] * x_mu_scaled[1]
  expect_equal(pred_mu[1, ], expected_mu_1, tolerance = 1e-10)

  # For first observation in new_data_sigma (with log intercept)
  # log_sigma[i] = log(intercept) + x_sigma * x_sigma_scaled[i]
  expected_log_sigma_1 <- log(posterior[, "log_sigma_intercept"]) + posterior[, "log_sigma_x_sigma"] * x_sigma_scaled[1]
  expect_equal(pred_log_sigma[1, ], expected_log_sigma_1, tolerance = 1e-10)
})


# ============================================================================ #
# LM-BASED VALIDATION TESTS
# ============================================================================ #
#
# These tests validate the unscaling transformation by comparing against lm():
# 1. Fit lm() with scaled predictors -> extract coefficients
# 2. Transform coefficients using transform_scale_samples()
# 3. Compare against lm() with unscaled predictors
#
# This approach validates both the implementation AND the derivation.
# ============================================================================ #

# Helper: Create formula_scale from data frame and variable names
# Creates nested structure matching JAGS_fit output: list(mu = list(mu_x1 = list(mean, sd)))
.make_formula_scale <- function(df, var_names, prefix = "mu") {
  param_scale <- list()
  for (var in var_names) {
    param_name <- paste0(prefix, "_", var)
    param_scale[[param_name]] <- list(
      mean = mean(df[[var]]),
      sd   = sd(df[[var]])
    )
  }
  # Return nested structure keyed by parameter name
  result <- list()
  result[[prefix]] <- param_scale
  result
}

# Helper: Convert lm coefficients to posterior matrix format (repeated rows)
# Uses the same naming convention as JAGS (__xXx__ for interactions)
.lm_coefs_to_posterior <- function(coefs, prefix = "mu", n_rep = 10) {
  # Convert names: "(Intercept)" -> "mu_intercept", "x1:x2" -> "mu_x1__xXx__x2"
  new_names <- names(coefs)
  new_names <- gsub("\\(Intercept\\)", "intercept", new_names)
  new_names <- gsub(":", "__xXx__", new_names)
  new_names <- paste0(prefix, "_", new_names)

  # Remove scale() wrapper from names if present
  new_names <- gsub("scale\\(([^)]+)\\)", "\\1", new_names)

  posterior <- matrix(rep(coefs, each = n_rep), nrow = n_rep, ncol = length(coefs))
  colnames(posterior) <- new_names
  posterior
}

# Helper to reorder lm coefficients to match posterior column order
.reorder_lm_coefs <- function(coef_unscaled, posterior_transformed) {
  # Build mapping from posterior names to lm names
  posterior_names <- colnames(posterior_transformed)
  lm_names <- sapply(posterior_names, function(nm) {
    # Remove mu_ prefix
    stripped <- sub("^mu_", "", nm)
    if (stripped == "intercept") return("(Intercept)")
    # Replace __xXx__ with :
    gsub("__xXx__", ":", stripped)
  })
  coef_unscaled[lm_names]
}


test_that("lm validation: simple standardization (one predictor)", {

  set.seed(42)
  df <- data.frame(
    x1 = rnorm(500, mean = 10, sd = 3),
    y  = rnorm(500)
  )
  df$y <- 5 + 2 * scale(df$x1) + rnorm(500, 0, 0.5)

  # Fit with scaled predictor
  fit_scaled <- lm(y ~ scale(x1), data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictor (ground truth)
  fit_unscaled <- lm(y ~ x1, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform scaled coefficients
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, "x1")
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: multiple predictors (no interaction)", {

  set.seed(43)
  df <- data.frame(
    x1 = rnorm(500, mean = 3, sd = 5),
    x2 = rnorm(500, mean = -10, sd = 2)
  )
  df$y <- 2 - 0.5 * scale(df$x1) + 1.5 * scale(df$x2) + rnorm(500, 0, 0.3)

  # Fit with scaled predictors
  fit_scaled <- lm(y ~ scale(x1) + scale(x2), data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 + x2, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: two-way interaction (both scaled)", {

  set.seed(44)
  df <- data.frame(
    x1 = rnorm(500, mean = 5, sd = 2),
    x2 = rnorm(500, mean = -3, sd = 4)
  )
  df$y <- 3 + 0.8 * scale(df$x1) - 0.5 * scale(df$x2) +
          0.3 * scale(df$x1) * scale(df$x2) + rnorm(500, 0, 0.5)

  # Fit with scaled predictors
  fit_scaled <- lm(y ~ scale(x1) * scale(x2), data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare all coefficients
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: two-way interaction (partial scaling)", {

  set.seed(45)
  df <- data.frame(
    x1 = rnorm(500, mean = 8, sd = 3),
    x2 = rnorm(500, mean = -2, sd = 5)
  )
  # Only x1 is scaled
  df$y <- 1 + 0.6 * scale(df$x1) - 0.4 * df$x2 +
          0.25 * scale(df$x1) * df$x2 + rnorm(500, 0, 0.4)

  # Fit with partial scaling (only x1 scaled)
  fit_scaled <- lm(y ~ scale(x1) * x2, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform - only x1 is in formula_scale
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, "x1")  # Only x1 scaled
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: three-way interaction (all scaled)", {

  set.seed(46)
  df <- data.frame(
    x1 = rnorm(500, mean = 3, sd = 2),
    x2 = rnorm(500, mean = -5, sd = 3),
    x3 = rnorm(500, mean = 10, sd = 4)
  )
  df$y <- 2 +
          0.5 * scale(df$x1) - 0.3 * scale(df$x2) + 0.4 * scale(df$x3) +
          0.2 * scale(df$x1) * scale(df$x2) +
          0.15 * scale(df$x1) * scale(df$x3) +
          0.1 * scale(df$x2) * scale(df$x3) +
          0.08 * scale(df$x1) * scale(df$x2) * scale(df$x3) +
          rnorm(500, 0, 0.3)

  # Fit with scaled predictors
  fit_scaled <- lm(y ~ scale(x1) * scale(x2) * scale(x3), data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2 * x3, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2", "x3"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare all coefficients
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: three-way interaction (partial scaling)", {

  set.seed(47)
  df <- data.frame(
    x1 = rnorm(500, mean = 4, sd = 2),
    x2 = rnorm(500, mean = -3, sd = 3),
    x3 = rnorm(500, mean = 7, sd = 1)  # This one not scaled
  )
  # x1 and x2 scaled, x3 not scaled
  df$y <- 1 +
          0.4 * scale(df$x1) - 0.2 * scale(df$x2) + 0.3 * df$x3 +
          0.15 * scale(df$x1) * scale(df$x2) +
          0.12 * scale(df$x1) * df$x3 +
          0.08 * scale(df$x2) * df$x3 +
          0.05 * scale(df$x1) * scale(df$x2) * df$x3 +
          rnorm(500, 0, 0.2)

  # Fit with partial scaling
  fit_scaled <- lm(y ~ scale(x1) * scale(x2) * x3, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2 * x3, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform - only x1 and x2 are scaled
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare all coefficients
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: four-way interaction", {

  set.seed(48)
  df <- data.frame(
    x1 = rnorm(1000, mean = 2, sd = 1),
    x2 = rnorm(1000, mean = -4, sd = 2),
    x3 = rnorm(1000, mean = 6, sd = 3),
    x4 = rnorm(1000, mean = -1, sd = 0.5)
  )
  # Complex model with 4-way interaction
  df$y <- 3 +
          0.3 * scale(df$x1) - 0.2 * scale(df$x2) +
          0.4 * scale(df$x3) - 0.1 * scale(df$x4) +
          0.05 * scale(df$x1) * scale(df$x2) * scale(df$x3) * scale(df$x4) +
          rnorm(1000, 0, 0.5)

  # Fit with scaled predictors
  fit_scaled <- lm(y ~ scale(x1) * scale(x2) * scale(x3) * scale(x4), data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2 * x3 * x4, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2", "x3", "x4"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare all coefficients
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: five-way interaction (warning test)", {

  set.seed(49)
  df <- data.frame(
    x1 = rnorm(2000, mean = 1, sd = 0.5),
    x2 = rnorm(2000, mean = -2, sd = 1),
    x3 = rnorm(2000, mean = 3, sd = 1.5),
    x4 = rnorm(2000, mean = -1, sd = 0.3),
    x5 = rnorm(2000, mean = 4, sd = 2)
  )
  df$y <- 2 +
          0.2 * scale(df$x1) - 0.1 * scale(df$x2) +
          0.3 * scale(df$x3) - 0.15 * scale(df$x4) + 0.1 * scale(df$x5) +
          0.02 * scale(df$x1) * scale(df$x2) * scale(df$x3) * scale(df$x4) * scale(df$x5) +
          rnorm(2000, 0, 0.5)

  # Fit with scaled predictors
  fit_scaled <- lm(y ~ scale(x1) * scale(x2) * scale(x3) * scale(x4) * scale(x5), data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2 * x3 * x4 * x5, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform - expect warning about 5+ way interaction
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2", "x3", "x4", "x5"))

  expect_warning(
    posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale),
    "5-way or higher interactions"
  )

  # Should still produce correct results despite warning
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: complex model from user example", {

  # This is the exact example pattern from the user's request
  set.seed(1)
  df_orig <- data.frame(
    x1 = rnorm(1000, mean =   3, sd = 5),
    x2 = rnorm(1000, mean = -10, sd = 80),
    x3 = rnorm(1000, mean = -20, sd = 0.07),
    x4 = rnorm(1000, mean =  50, sd = 30),
    x5 = rnorm(1000, mean =  20, sd = 0.2)
  )

  # DGP with specific structure
  df_orig$y <- with(
    df_orig,
    5 - 0.1 * scale(x1) + 0.2 * scale(x2) + 0.3 * scale(x1) * scale(x2) -
    0.25 * scale(x3) * scale(x4) * scale(x5) + 0.40 * scale(x3) * scale(x4) +
    rnorm(1000, 0, 1)
  )

  # Fit the model with scaled predictors (matching DGP)
  fit_scaled <- lm(y ~ scale(x1) * scale(x2) + scale(x3) * scale(x4) * scale(x5), data = df_orig)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2 + x3 * x4 * x5, data = df_orig)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df_orig, c("x1", "x2", "x3", "x4", "x5"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare all coefficients
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: factor + scaled continuous interaction", {

  set.seed(50)
  df <- data.frame(
    x1 = rnorm(500, mean = 5, sd = 3),
    f1 = factor(sample(letters[1:2], 500, TRUE))
  )
  df$y <- 2 + 0.5 * scale(df$x1) +
          ifelse(df$f1 == "b", 0.3, 0) +
          ifelse(df$f1 == "b", 0.2, 0) * scale(df$x1) +
          rnorm(500, 0, 0.4)

  # Fit with scaled predictor
  fit_scaled <- lm(y ~ scale(x1) * f1, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictor (ground truth)
  fit_unscaled <- lm(y ~ x1 * f1, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform - only x1 is scaled (f1 is factor, not scaled)
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, "x1")
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: factor + unscaled continuous interaction", {

  set.seed(51)
  df <- data.frame(
    x1 = rnorm(500, mean = 8, sd = 2),  # Will NOT be scaled
    x2 = rnorm(500, mean = -3, sd = 4), # Will be scaled
    f1 = factor(sample(letters[1:2], 500, TRUE))
  )
  df$y <- 1 + 0.3 * df$x1 + 0.4 * scale(df$x2) +
          ifelse(df$f1 == "b", 0.5, 0) +
          ifelse(df$f1 == "b", 0.1, 0) * df$x1 +
          ifelse(df$f1 == "b", 0.15, 0) * scale(df$x2) +
          rnorm(500, 0, 0.3)

  # Fit with partial scaling (only x2 scaled)
  fit_scaled <- lm(y ~ x1 * f1 + scale(x2) * f1, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * f1 + x2 * f1, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform - only x2 is scaled
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, "x2")
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: multi-level factor with scaled continuous", {

  set.seed(52)
  df <- data.frame(
    x1 = rnorm(600, mean = 3, sd = 5),
    f1 = factor(sample(letters[1:3], 600, TRUE))
  )
  df$y <- 2 + 0.6 * scale(df$x1) +
          ifelse(df$f1 == "b", 0.4, ifelse(df$f1 == "c", -0.3, 0)) +
          ifelse(df$f1 == "b", 0.2, ifelse(df$f1 == "c", 0.1, 0)) * scale(df$x1) +
          rnorm(600, 0, 0.5)

  # Fit with scaled predictor
  fit_scaled <- lm(y ~ scale(x1) * f1, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictor (ground truth)
  fit_unscaled <- lm(y ~ x1 * f1, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, "x1")
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: two factors with scaled continuous interaction", {

  set.seed(53)
  df <- data.frame(
    x1 = rnorm(800, mean = 10, sd = 4),
    f1 = factor(sample(letters[1:2], 800, TRUE)),
    f2 = factor(sample(letters[1:3], 800, TRUE))
  )
  # Complex model with factor-factor and factor-continuous interactions
  df$y <- 3 + 0.5 * scale(df$x1) + rnorm(800, 0, 0.6)

  # Fit with scaled predictor - full three-way interaction
  fit_scaled <- lm(y ~ scale(x1) * f1 * f2, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictor (ground truth)
  fit_unscaled <- lm(y ~ x1 * f1 * f2, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, "x1")
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: complex model with factors and mixed scaling", {

  # Comprehensive test with the user's data structure
  set.seed(1)
  df <- data.frame(
    x1 = rnorm(1000, mean =   3, sd = 5),
    x2 = rnorm(1000, mean = -10, sd = 80),
    x3 = rnorm(1000, mean = -20, sd = 0.07),
    x4 = rnorm(1000, mean =  50, sd = 30),
    x5 = rnorm(1000, mean =  20, sd = 0.2),
    f1 = factor(sample(letters[1:2], 1000, TRUE)),
    f2 = factor(sample(letters[1:3], 1000, TRUE))
  )

  # Model with scaled continuous, unscaled continuous, and factors
  # x1, x2, x3 are scaled; x4, x5 are NOT scaled
  df$y <- 5 +
          0.3 * scale(df$x1) - 0.2 * scale(df$x2) + 0.1 * scale(df$x3) +
          0.15 * df$x4 - 0.1 * df$x5 +
          ifelse(df$f1 == "b", 0.4, 0) +
          0.2 * scale(df$x1) * ifelse(df$f1 == "b", 1, 0) +
          0.1 * df$x4 * ifelse(df$f1 == "b", 1, 0) +
          rnorm(1000, 0, 1)

  # Fit with partial scaling
  fit_scaled <- lm(y ~ scale(x1) * f1 + scale(x2) + scale(x3) + x4 * f1 + x5, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * f1 + x2 + x3 + x4 * f1 + x5, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform - only x1, x2, x3 are scaled
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2", "x3"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})


test_that("lm validation: factor interactions with multiple scaled continuous", {

  set.seed(54)
  df <- data.frame(
    x1 = rnorm(800, mean = 5, sd = 3),
    x2 = rnorm(800, mean = -2, sd = 6),
    f1 = factor(sample(letters[1:2], 800, TRUE))
  )
  # Continuous-continuous and continuous-factor interactions
  df$y <- 2 +
          0.4 * scale(df$x1) - 0.3 * scale(df$x2) +
          0.25 * scale(df$x1) * scale(df$x2) +
          ifelse(df$f1 == "b", 0.5, 0) +
          0.15 * scale(df$x1) * ifelse(df$f1 == "b", 1, 0) +
          0.1 * scale(df$x2) * ifelse(df$f1 == "b", 1, 0) +
          rnorm(800, 0, 0.5)

  # Fit with scaled predictors - three-way interaction x1 * x2 * f1
  fit_scaled <- lm(y ~ scale(x1) * scale(x2) * f1, data = df)
  coef_scaled <- coef(fit_scaled)

  # Fit with unscaled predictors (ground truth)
  fit_unscaled <- lm(y ~ x1 * x2 * f1, data = df)
  coef_unscaled <- coef(fit_unscaled)

  # Transform
  posterior_scaled <- .lm_coefs_to_posterior(coef_scaled)
  formula_scale <- .make_formula_scale(df, c("x1", "x2"))
  posterior_transformed <- transform_scale_samples(posterior_scaled, formula_scale)

  # Compare
  expect_equal(
    unname(posterior_transformed[1, ]),
    unname(.reorder_lm_coefs(coef_unscaled, posterior_transformed)),
    tolerance = 1e-10
  )
})

