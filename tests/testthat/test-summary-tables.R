# ============================================================================ #
# TEST FILE: Summary Tables
# ============================================================================ #
#
# PURPOSE:
#   Tests for summary table functions including ensemble_estimates_table,
#   ensemble_inference_table, ensemble_summary_table, ensemble_diagnostics_table,
#   model_summary_table, and print methods.
#
# DEPENDENCIES:
#   - rjags, bridgesampling: For tests using pre-fitted models
#   - common-functions.R: temp_fits_dir, test_reference_table
#
# SKIP CONDITIONS:
#   - skip_if_no_fits(): Pre-fitted models required for most tests
#   - skip_if_not_installed("rjags"), skip_if_not_installed("bridgesampling")
#
# MODELS/FIXTURES:
#   - fit_summary*, fit_simple_normal, fit_simple_spike, fit_orthonormal_*
#
# TAGS: @evaluation, @summary-tables
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "summary-tables")
source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: ensemble_estimates_table tests
# ============================================================================ #
test_that("ensemble_estimates_table handles matrix posteriors", {

  skip_if_no_fits()
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Load fits with margliks for creating mixed posteriors
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE), "omega" = c(TRUE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  # Test basic table creation
  estimates_table <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = c("m", "omega")
  )

  test_reference_table(estimates_table, "ensemble_estimates_basic.txt")

  # Test with custom probs
  estimates_table_probs <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = c("m", "omega"),
    probs = c(0.10, 0.50, 0.90)
  )

  test_reference_table(estimates_table_probs, "ensemble_estimates_custom_probs.txt")

})


test_that("ensemble_estimates_table handles transform_factors", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load orthonormal models with marginal likelihoods
  fit_orthonormal_0 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0.RDS"))
  marglik_orthonormal_0 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_0.RDS"))

  fit_orthonormal_1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik_orthonormal_1 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_1.RDS"))

  models <- list(
    list(fit = fit_orthonormal_0, marglik = marglik_orthonormal_0, prior_weights = 1),
    list(fit = fit_orthonormal_1, marglik = marglik_orthonormal_1, prior_weights = 1)
  )

  # Get factor parameter names from the model
  prior_list <- attr(fit_orthonormal_1, "prior_list")
  factor_params <- names(prior_list)[sapply(prior_list, is.prior.factor)]

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = factor_params,
    is_null_list = setNames(list(c(TRUE, FALSE)), factor_params),
    seed = 1,
    n_samples = 1000
  )

  # Test with transform_factors = TRUE
  estimates_table_transform <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = factor_params,
    transform_factors = TRUE
  )

  test_reference_table(estimates_table_transform, "ensemble_estimates_transform_factors.txt")

})


test_that("ensemble_estimates_table handles formula posteriors", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Use orthonormal models (have formulas and marginal likelihoods)
  fit_formula <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0.RDS"))
  marglik_formula <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_0.RDS"))

  fit_formula2 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik_formula2 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_1.RDS"))

  models <- list(
    list(fit = fit_formula, marglik = marglik_formula, prior_weights = 1),
    list(fit = fit_formula2, marglik = marglik_formula2, prior_weights = 1)
  )

  prior_list <- attr(fit_formula, "prior_list")
  params <- names(prior_list)[!sapply(prior_list, is.null)]

  is_null_list <- setNames(
    lapply(params, function(p) c(FALSE, FALSE)),
    params
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = params,
    is_null_list = is_null_list,
    seed = 1,
    n_samples = 1000
  )

  # Test with formula_prefix = TRUE
  estimates_prefix_true <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = params,
    formula_prefix = TRUE
  )

  # Test with formula_prefix = FALSE
  estimates_prefix_false <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = params,
    formula_prefix = FALSE
  )

  test_reference_table(estimates_prefix_true, "ensemble_estimates_formula_prefix_true.txt")
  test_reference_table(estimates_prefix_false, "ensemble_estimates_formula_prefix_false.txt")

})


# ============================================================================ #
# SECTION 2: ensemble_inference_table tests
# ============================================================================ #
test_that("ensemble_inference_table handles multiple parameters", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1)
  )

  inference <- ensemble_inference(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE), "omega" = c(TRUE, FALSE))
  )

  # Basic table
  inference_table <- ensemble_inference_table(inference, names(inference))
  test_reference_table(inference_table, "ensemble_inference_basic.txt")

  # With logBF
  inference_table_log <- ensemble_inference_table(inference, names(inference), logBF = TRUE)
  test_reference_table(inference_table_log, "ensemble_inference_logBF.txt")

  # With BF01
  inference_table_bf01 <- ensemble_inference_table(inference, names(inference), BF01 = TRUE)
  test_reference_table(inference_table_bf01, "ensemble_inference_BF01.txt")

  # With both
  inference_table_both <- ensemble_inference_table(inference, names(inference), logBF = TRUE, BF01 = TRUE)
  test_reference_table(inference_table_both, "ensemble_inference_both.txt")

})


# ============================================================================ #
# SECTION 3: ensemble_summary_table and ensemble_diagnostics_table tests
# ============================================================================ #
test_that("ensemble_summary_table handles different model configurations", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Use models with and without spike-at-zero to test remove_spike_0
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  models <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_normal)),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_spike))
  )
  models <- models_inference(models)

  # Test summary table
  summary_table <- ensemble_summary_table(models, c("m", "s"))
  test_reference_table(summary_table, "ensemble_summary_basic.txt")

  # Test with short_name
  summary_table_short <- ensemble_summary_table(models, c("m", "s"), short_name = TRUE)
  test_reference_table(summary_table_short, "ensemble_summary_short_name.txt")

  # Test with logBF and BF01
  summary_table_bf <- ensemble_summary_table(models, c("m", "s"), logBF = TRUE, BF01 = TRUE)
  test_reference_table(summary_table_bf, "ensemble_summary_bf_options.txt")

  # Test with remove_spike_0
  summary_table_no_spike <- ensemble_summary_table(models, c("m", "s"), remove_spike_0 = FALSE)
  test_reference_table(summary_table_no_spike, "ensemble_summary_no_spike.txt")

})


test_that("ensemble_summary_table handles parameters as list", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  models <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_normal)),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_spike))
  )
  models <- models_inference(models)

  # Test with parameters supplied as a list
  pars <- list("m" = "m", "renamed 2" = "s")
  summary_table_list <- ensemble_summary_table(models, pars)
  test_reference_table(summary_table_list, "ensemble_summary_params_list.txt")

})


test_that("ensemble_diagnostics_table handles different configurations", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Use models with and without spike-at-zero to test remove_spike_0
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  models <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_normal)),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_spike))
  )
  models <- models_inference(models)

  # Test diagnostics table
  diagnostics_table <- ensemble_diagnostics_table(models, c("m", "s"))
  test_reference_table(diagnostics_table, "ensemble_diagnostics_basic.txt")

  # Test with short_name
  diagnostics_short <- ensemble_diagnostics_table(models, c("m", "s"), short_name = TRUE)
  test_reference_table(diagnostics_short, "ensemble_diagnostics_short_name.txt")

  # Test with remove_spike_0
  diagnostics_no_spike <- ensemble_diagnostics_table(models, c("m", "s"), remove_spike_0 = FALSE)
  test_reference_table(diagnostics_no_spike, "ensemble_diagnostics_no_spike.txt")

})


# ============================================================================ #
# SECTION 4: marginal_estimates_table tests
# ============================================================================ #
test_that("marginal_estimates_table handles various inputs", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Create sample data for marginal inference testing
  set.seed(1)
  samples <- list(
    mu = rnorm(1000, 0, 1)
  )

  inference <- list(
    mu = structure(list(
      BF = 2.5,
      prior_probs = c(0.5, 0.5),
      post_probs = c(0.4, 0.6)
    ), class = c("list", "marginal_inference"))
  )

  attr(inference$mu, "is_null") <- c(TRUE, FALSE)
  attr(inference$mu, "prior_list") <- list(
    prior("spike", list(0)),
    prior("normal", list(0, 1))
  )

  marginal_table <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = "mu"
  )

  test_reference_table(marginal_table, "marginal_estimates_basic.txt")

  # With logBF
  marginal_table_log <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = "mu",
    logBF = TRUE
  )
  test_reference_table(marginal_table_log, "marginal_estimates_logBF.txt")

  # With BF01
  marginal_table_bf01 <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = "mu",
    BF01 = TRUE
  )
  test_reference_table(marginal_table_bf01, "marginal_estimates_BF01.txt")

})


# ============================================================================ #
# SECTION 5: model_summary_table tests
# ============================================================================ #
test_that("model_summary_table handles various configurations", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Use model with spike-at-zero to test remove_spike_0
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  model <- list(
    fit = fit_simple_spike,
    marglik = marglik_simple_spike,
    prior_weights = 1,
    fit_summary = runjags_estimates_table(fit_simple_spike)
  )
  model_list <- list(model)
  model_list <- models_inference(model_list)

  # Basic model summary
  summary_table <- model_summary_table(model_list[[1]])
  test_reference_table(summary_table, "model_summary_basic.txt")

  # With short_name
  summary_short <- model_summary_table(model_list[[1]], short_name = TRUE)
  test_reference_table(summary_short, "model_summary_short_name.txt")

  # With remove_spike_0 (should remove 'm' which has spike at zero)
  summary_no_spike <- model_summary_table(model_list[[1]], remove_spike_0 = TRUE)
  test_reference_table(summary_no_spike, "model_summary_no_spike.txt")

})


# ============================================================================ #
# SECTION 6: update.BayesTools_table tests
# ============================================================================ #
test_that("update.BayesTools_table works correctly", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary0)),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary1))
  )
  models <- models_inference(models)

  inference <- ensemble_inference(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE), "omega" = c(TRUE, FALSE))
  )

  # Create inference table
  inference_table <- ensemble_inference_table(inference, names(inference))

  # Update with new title
  updated_table <- update(inference_table, title = "Updated Title")
  test_reference_table(updated_table, "update_table_new_title.txt")

  # Update with footnotes
  updated_footnotes <- update(inference_table, footnotes = "This is a footnote")
  test_reference_table(updated_footnotes, "update_table_footnotes.txt")

  # Update with warnings
  updated_warnings <- update(inference_table, warnings = "This is a warning")
  test_reference_table(updated_warnings, "update_table_warnings.txt")

  # Update with logBF
  updated_logbf <- update(inference_table, logBF = TRUE)
  test_reference_table(updated_logbf, "update_table_logBF.txt")

  # Update with BF01
  updated_bf01 <- update(inference_table, BF01 = TRUE)
  test_reference_table(updated_bf01, "update_table_BF01.txt")

})
