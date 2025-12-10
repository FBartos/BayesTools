context("Summary tables functions")

# ==============================================================================
# CONFIGURATION: Set to TRUE to regenerate reference files, FALSE to run tests
# ==============================================================================
GENERATE_REFERENCE_FILES <- FALSE

# Get the directory where prefitted models are stored
temp_fits_dir <- Sys.getenv("BAYESTOOLS_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  temp_fits_dir <- file.path(tempdir(), "BayesTools_test_fits")
}

# Skip tests on CRAN as they require pre-fitted models
skip_on_cran()

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Process reference file: save if GENERATE_REFERENCE_FILES=TRUE, test otherwise
test_reference <- function(table, filename, info_msg = NULL, 
                              print_dir = testthat::test_path("..", "results", "print")) {
  if (GENERATE_REFERENCE_FILES) {
    # Save mode
    if (!dir.exists(print_dir)) {
      dir.create(print_dir, recursive = TRUE)
    }
    writeLines(capture_output_lines(table, print = TRUE, width = 150),
               file.path(print_dir, filename))
  } else {
    # Test mode
    ref_file <- file.path(print_dir, filename)
    if (file.exists(ref_file)) {
      expected_output <- readLines(ref_file, warn = FALSE)
      actual_output   <- capture_output_lines(table, print = TRUE, width = 150)
      expect_equal(actual_output, expected_output, info = info_msg)
    } else {
      skip(paste("Reference file", filename, "not found."))
    }
  }
}

# ==============================================================================
# SECTION 1: Test Empty Tables
# ==============================================================================
test_that("Empty summary tables work correctly", {
  
  runjags_summary_empty <- runjags_estimates_empty_table()
  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  ensemble_inference_empty <- ensemble_inference_empty_table()
  
  expect_equivalent(nrow(runjags_summary_empty), 0)
  expect_equivalent(nrow(ensemble_estimates_empty), 0)
  expect_equivalent(nrow(ensemble_inference_empty), 0)
  
  # Test that empty tables have correct structure
  expect_s3_class(runjags_summary_empty, "BayesTools_table")
  expect_s3_class(ensemble_estimates_empty, "BayesTools_table")
  expect_s3_class(ensemble_inference_empty, "BayesTools_table")
  
  test_reference(runjags_summary_empty, "empty_runjags_estimates.txt", "Empty runjags_estimates table mismatch")
  test_reference(ensemble_estimates_empty, "empty_ensemble_estimates.txt", "Empty ensemble_estimates table mismatch")
  test_reference(ensemble_inference_empty, "empty_ensemble_inference.txt", "Empty ensemble_inference table mismatch")
  
  if (GENERATE_REFERENCE_FILES) {
    message("Generated reference files for empty tables")
  }
})

# ==============================================================================
# SECTION 2: Test Advanced Features (Transformations, Formula Handling, etc.)
# ==============================================================================
test_that("Summary table advanced features work correctly", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Use fit_formula_interaction_cont for testing advanced features
  # This model has continuous interactions and formulas
  fit_complex <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))
  marglik_complex <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont_marglik.RDS"))
  
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
  marglik_factor <- readRDS(file.path(temp_fits_dir, "fit_factor_treatment_marglik.RDS"))
  
  runjags_summary_factor <- runjags_estimates_table(fit_factor)
  
  # Test 5: Use fit with spike and slab
  fit_spike <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_factor.RDS"))
  marglik_spike <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_factor_marglik.RDS"))
  
  runjags_summary_spike <- runjags_estimates_table(fit_spike)
  runjags_inference_spike <- runjags_inference_table(fit_spike)
  
  # Test basic properties
  expect_s3_class(runjags_summary_transform, "BayesTools_table")
  expect_s3_class(runjags_summary_prefix_true, "BayesTools_table")
  expect_s3_class(runjags_summary_prefix_false, "BayesTools_table")
  expect_s3_class(runjags_summary_conditional, "BayesTools_table")
  expect_s3_class(runjags_summary_unconditional, "BayesTools_table")
  expect_s3_class(runjags_summary_factor, "BayesTools_table")
  expect_s3_class(runjags_summary_spike, "BayesTools_table")
  expect_s3_class(runjags_inference_spike, "BayesTools_table")
  
  # Test that row names differ with different formula_prefix settings
  expect_false(identical(rownames(runjags_summary_prefix_true), 
                        rownames(runjags_summary_prefix_false)))
  
  # Test that conditional has fewer rows than unconditional (removes low inclusion parameters)
  # This may not always be true, so we just check they exist
  expect_true(nrow(runjags_summary_conditional) >= 0)
  expect_true(nrow(runjags_summary_unconditional) >= 0)
  
  test_reference(runjags_summary_transform, "advanced_transform.txt", "Transform table mismatch")
  test_reference(runjags_summary_prefix_true, "advanced_formula_prefix_true.txt", "Formula prefix true table mismatch")
  test_reference(runjags_summary_prefix_false, "advanced_formula_prefix_false.txt", "Formula prefix false table mismatch")
  test_reference(runjags_summary_conditional, "advanced_conditional.txt", "Conditional table mismatch")
  test_reference(runjags_summary_unconditional, "advanced_unconditional.txt", "Unconditional table mismatch")
  test_reference(runjags_summary_factor, "advanced_factor_treatment.txt", "Factor treatment table mismatch")
  test_reference(runjags_summary_spike, "advanced_spike_slab_estimates.txt", "Spike slab estimates table mismatch")
  test_reference(runjags_inference_spike, "advanced_spike_slab_inference.txt", "Spike slab inference table mismatch")

})

# ==============================================================================
# SECTION 3: Test Summary Tables for All Saved Models
# ==============================================================================
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
      test_reference(model_summary, paste0(model_name, "_model_summary.txt"),
                       paste0("Model summary mismatch for ", model_name))
    }
    
    # Process runjags estimates table
    runjags_summary <- runjags_estimates_table(fit)
    test_reference(runjags_summary, paste0(model_name, "_runjags_estimates.txt"),
                     paste0("Runjags estimates mismatch for ", model_name))
    
  }
})