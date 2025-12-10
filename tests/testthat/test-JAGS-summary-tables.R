context("JAGS summary tables functions")

REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-summary-tables")
source(testthat::test_path("common-functions.R"))

# ============================================================================ #
# SECTION 1: Test Empty Tables
# ============================================================================ #
test_that("Empty summary tables work correctly", {

  runjags_summary_empty <- runjags_estimates_empty_table()

  expect_equivalent(nrow(runjags_summary_empty), 0)

  # Test that empty tables have correct structure
  expect_s3_class(runjags_summary_empty, "BayesTools_table")

  test_reference_table(runjags_summary_empty, "empty_runjags_estimates.txt", "Empty runjags_estimates table mismatch")
})

# ============================================================================ #
# SECTION 2: Test Advanced Features (Transformations, Formula Handling, etc.)
# ============================================================================ #
test_that("Summary table advanced features work correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Use fit_formula_interaction_cont for testing advanced features
  # This model has continuous interactions and formulas
  fit_complex <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))

  # Test 1: Parameter transformations
  runjags_summary_transform <- runjags_estimates_table(
    fit_complex,
    transformations = list("mu_intercept" = list(fun = exp))
  )

  # Test 2: Formula handling with prefix
  runjags_summary_prefix_true <- runjags_estimates_table(fit_complex, formula_prefix = TRUE)
  runjags_summary_prefix_false <- runjags_estimates_table(fit_complex, formula_prefix = FALSE)

  # Test 3: Conditional vs unconditional
  runjags_summary_conditional <- runjags_estimates_table(fit_complex, conditional = TRUE)
  runjags_summary_unconditional <- runjags_estimates_table(fit_complex, conditional = FALSE)

  # Test 4: Factor transformations (use fit_factor_treatment for this)
  fit_factor <- readRDS(file.path(temp_fits_dir, "fit_factor_treatment.RDS"))

  runjags_summary_factor <- runjags_estimates_table(fit_factor)

  # Test 5: Use fit with spike and slab
  fit_spike <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_factor.RDS"))

  runjags_summary_spike <- runjags_estimates_table(fit_spike)
  runjags_inference_spike <- runjags_inference_table(fit_spike)

  # Test 6: Orthonormal contrast transformations to differences from the mean
  fit_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS"))

  runjags_summary_orthonormal <- suppressMessages(runjags_estimates_table(fit_orthonormal, transform_factors = TRUE))

  # Test 7: Custom transformations with transform_factors = FALSE
  # Use a model with factor parameters for transformation testing
  runjags_summary_custom_transform <- suppressMessages(runjags_estimates_table(
    fit_factor,
    transform_factors = FALSE,
    transformations = list("mu_x_fac3t[2]" = list(fun = exp))
  ))

  # Test 8: Conditional with remove_inclusion
  runjags_summary_remove_inclusion <- suppressMessages(runjags_estimates_table(
    fit_spike,
    transform_factors = TRUE,
    conditional = TRUE,
    remove_inclusion = TRUE
  ))

  # Test 9: Custom probs parameter
  runjags_summary_custom_probs <- runjags_estimates_table(fit_complex)

  # Test basic properties
  expect_s3_class(runjags_summary_transform, "BayesTools_table")
  expect_s3_class(runjags_summary_prefix_true, "BayesTools_table")
  expect_s3_class(runjags_summary_prefix_false, "BayesTools_table")
  expect_s3_class(runjags_summary_conditional, "BayesTools_table")
  expect_s3_class(runjags_summary_unconditional, "BayesTools_table")
  expect_s3_class(runjags_summary_factor, "BayesTools_table")
  expect_s3_class(runjags_summary_spike, "BayesTools_table")
  expect_s3_class(runjags_inference_spike, "BayesTools_table")
  expect_s3_class(runjags_summary_orthonormal, "BayesTools_table")
  expect_s3_class(runjags_summary_custom_transform, "BayesTools_table")
  expect_s3_class(runjags_summary_remove_inclusion, "BayesTools_table")
  expect_s3_class(runjags_summary_custom_probs, "BayesTools_table")

  # Test that row names differ with different formula_prefix settings
  expect_false(identical(rownames(runjags_summary_prefix_true),
                        rownames(runjags_summary_prefix_false)))

  # Test that remove_inclusion reduces the number of rows
  expect_true(nrow(runjags_summary_remove_inclusion) <= nrow(runjags_summary_spike))

  test_reference_table(runjags_summary_transform, "advanced_transform.txt", "Transform table mismatch")
  test_reference_table(runjags_summary_prefix_true, "advanced_formula_prefix_true.txt", "Formula prefix true table mismatch")
  test_reference_table(runjags_summary_prefix_false, "advanced_formula_prefix_false.txt", "Formula prefix false table mismatch")
  test_reference_table(runjags_summary_conditional, "advanced_conditional.txt", "Conditional table mismatch")
  test_reference_table(runjags_summary_unconditional, "advanced_unconditional.txt", "Unconditional table mismatch")
  test_reference_table(runjags_summary_factor, "advanced_factor_treatment.txt", "Factor treatment table mismatch")
  test_reference_table(runjags_summary_spike, "advanced_spike_slab_estimates.txt", "Spike slab estimates table mismatch")
  test_reference_table(runjags_inference_spike, "advanced_spike_slab_inference.txt", "Spike slab inference table mismatch")
  test_reference_table(runjags_summary_orthonormal, "advanced_orthonormal_transform.txt", "Orthonormal transform table mismatch")
  test_reference_table(runjags_summary_custom_transform, "advanced_custom_transform.txt", "Custom transform table mismatch")
  test_reference_table(runjags_summary_remove_inclusion, "advanced_remove_inclusion.txt", "Remove inclusion table mismatch")
  test_reference_table(runjags_summary_custom_probs, "advanced_custom_probs.txt", "Custom probs table mismatch")

})

# ============================================================================ #
# SECTION 3: Test Summary Tables for All Saved Models
# ============================================================================ #
test_that("Summary tables for all saved models", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  # Load model registry to get list of all fitted models
  registry_file <- file.path(temp_fits_dir, "model_registry.RDS")

  model_registry <- readRDS(registry_file)
  model_names <- model_registry$model_name

  print_dir <- testthat::test_path("..", "results", "print")

  for (model_name in model_names) {
    fit_file <- file.path(temp_fits_dir, paste0(model_name, ".RDS"))
    marglik_file <- file.path(temp_fits_dir, paste0(model_name, "_marglik.RDS"))

    fit <- readRDS(fit_file)
    has_marglik <- file.exists(marglik_file)

    if (has_marglik) {
      marglik <- readRDS(marglik_file)
    }

    # Process model summary table
    if (has_marglik) {
      model_list <- list(
        list(fit = fit, marglik = marglik, prior_weights = 1,
             fit_summary = runjags_estimates_table(fit))
      )
      model_list <- models_inference(model_list)
      model_summary <- model_summary_table(model_list[[1]])
      test_reference_table(model_summary, paste0(model_name, "_model_summary.txt"),
                       paste0("Model summary mismatch for ", model_name))
    }

    # Process runjags estimates table
    runjags_summary <- runjags_estimates_table(fit)
    test_reference_table(runjags_summary, paste0(model_name, "_runjags_estimates.txt"),
                     paste0("Runjags estimates mismatch for ", model_name))

  }
})
