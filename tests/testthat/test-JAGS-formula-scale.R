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

test_that("JAGS_fit propagates formula_scale correctly", {
  skip_if_not_installed("rjags")
  
  # Simple test data
  set.seed(4)
  df <- data.frame(
    y = rnorm(30),
    x = rnorm(30, mean = 10, sd = 3)
  )
  
  # Simple model
  model_syntax <- "
    model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], tau)
      }
      tau ~ dgamma(0.01, 0.01)
    }
  "
  
  data_list <- list(y = df$y)
  
  prior_list <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x"         = prior("normal", list(0, 1))
  )
  
  formula_list <- list(mu = ~ x)
  formula_data_list <- list(mu = df)
  formula_prior_list <- list(mu = prior_list)
  formula_scale_list <- list(mu = list(x = TRUE))
  
  # Fit model (with minimal iterations for speed)
  fit <- JAGS_fit(
    model_syntax       = model_syntax,
    data               = data_list,
    formula_list       = formula_list,
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    chains             = 2,
    adapt              = 100,
    burnin             = 100,
    sample             = 200,
    silent             = TRUE
  )
  
  # Check that formula_scale is stored in fit object
  expect_true(!is.null(attr(fit, "formula_scale")))
  expect_true("mu_x" %in% names(attr(fit, "formula_scale")))
  expect_equal(names(attr(fit, "formula_scale")$mu_x), c("mean", "sd"))
  
  # Verify values
  expect_equal(attr(fit, "formula_scale")$mu_x$mean, 10, tolerance = 0.5)
  expect_equal(attr(fit, "formula_scale")$mu_x$sd, 3, tolerance = 0.5)
})
