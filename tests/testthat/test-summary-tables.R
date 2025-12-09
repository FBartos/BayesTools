context("Summary tables functions")

# Get the directory where prefitted models are stored
temp_fits_dir <- Sys.getenv("BAYESTOOLS_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  temp_fits_dir <- file.path(tempdir(), "BayesTools_test_fits")
}

# Skip tests on CRAN as they require pre-fitted models
skip_on_cran()

test_that("Summary table printing works", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Check if reference files exist
  # test_path() constructs paths relative to tests/testthat/, so we use ".." to reach tests/results/print
  print_dir <- testthat::test_path("..", "results", "print")
  if (!dir.exists(print_dir)) {
    skip("Print reference directory not found. Set if(FALSE) to if(TRUE) in generation block and run to generate.")
  }
  
  # ============================================================================
  # Test 1-6: Basic summary tables (weightfunction priors)
  # ============================================================================
  fit0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))
  
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))
  
  fit2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))
  
  models_wf <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
  )
  models_wf <- models_inference(models_wf)
  inference_wf <- ensemble_inference(model_list = models_wf, parameters = c("m", "omega"), 
                                      is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors_wf <- mix_posteriors(model_list = models_wf, parameters = c("m", "omega"), 
                                          is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  
  # ============================================================================
  # Test 7-11: Simple priors
  # ============================================================================
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))
  
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))
  
  models_simple <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_simple_normal)),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_simple_spike))
  )
  models_simple <- models_inference(models_simple)
  inference_simple <- ensemble_inference(model_list = models_simple, parameters = c("m", "s"), 
                                          is_null_list = list("m" = 0, "s" = 1), conditional = FALSE)
  
  # ============================================================================
  # Test 12-13: Transformations
  # ============================================================================
  runjags_summary_transform <- runjags_estimates_table(fit_simple_normal, 
                                                         transformations = list("m" = list(fun = exp)))
  
  # ============================================================================
  # Test 14-16: Empty tables
  # ============================================================================
  runjags_summary_empty <- runjags_estimates_empty_table()
  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  ensemble_inference_empty <- ensemble_inference_empty_table()
  
  # ============================================================================
  # Create list of all tables to test
  # ============================================================================
  fits <- list(
    # 1-6: Weightfunction models (basic test)
    model_summary_table(models_wf[[2]]),
    runjags_estimates_table(fit1),
    ensemble_estimates_table(mixed_posteriors_wf, parameters = c("m", "omega"), probs = c(.025, 0.95)),
    ensemble_inference_table(inference_wf, names(inference_wf)),
    ensemble_summary_table(models_wf, c("m", "omega")),
    ensemble_diagnostics_table(models_wf, c("m", "omega")),
    
    # 7-11: Simple models
    model_summary_table(models_simple[[1]]),
    runjags_estimates_table(fit_simple_normal),
    ensemble_inference_table(inference_simple, names(inference_simple)),
    ensemble_summary_table(models_simple, c("m", "s")),
    ensemble_diagnostics_table(models_simple, c("m", "s")),
    
    # 12-13: Transformations
    runjags_summary_transform,
    runjags_estimates_table(fit2, transformations = list("m" = list(fun = exp))),
    
    # 14-16: Empty tables
    runjags_summary_empty,
    ensemble_estimates_empty,
    ensemble_inference_empty
  )
  
  # Compare printed output with saved reference files
  for(i in seq_along(fits)){
    ref_file <- file.path(print_dir, paste0(i, ".txt"))
    if (!file.exists(ref_file)) {
      skip(paste0("Reference file ", i, ".txt not found. Set if(FALSE) to if(TRUE) in generation block and run to generate."))
    }
    
    # Use readLines for simpler and more robust file reading
    expected_output <- readLines(ref_file, warn = FALSE)
    actual_output <- capture_output_lines(fits[[i]], print = TRUE, width = 150)
    
    expect_equal(actual_output, expected_output, 
                 info = paste0("Print output mismatch for table ", i))
  }
})

# ==============================================================================
# REFERENCE OUTPUT GENERATION
# ==============================================================================
# To regenerate reference print output files, set the condition below to TRUE
# and run this test file. This will overwrite existing reference files.
if (FALSE) {
  
  # Ensure models are fitted first
  if (!dir.exists(temp_fits_dir) || length(list.files(temp_fits_dir, pattern = "\\.RDS$")) == 0) {
    message("Pre-fitted models not found. Running test-00-model-fits.R...")
    testthat::test_file(testthat::test_path("test-00-model-fits.R"))
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # ============================================================================
  # Generate 1-6: Basic summary tables (weightfunction priors)
  # ============================================================================
  fit0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))
  
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))
  
  fit2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))
  
  models_wf <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
  )
  models_wf <- models_inference(models_wf)
  
  inference_wf <- ensemble_inference(model_list = models_wf, parameters = c("m", "omega"), 
                                      is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors_wf <- mix_posteriors(model_list = models_wf, parameters = c("m", "omega"), 
                                          is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  
  # ============================================================================
  # Generate 7-11: Simple priors
  # ============================================================================
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))
  
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))
  
  models_simple <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_simple_normal)),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_simple_spike))
  )
  models_simple <- models_inference(models_simple)
  inference_simple <- ensemble_inference(model_list = models_simple, parameters = c("m", "s"), 
                                          is_null_list = list("m" = 0, "s" = 1), conditional = FALSE)
  
  # ============================================================================
  # Generate 12-13: Transformations
  # ============================================================================
  runjags_summary_transform <- runjags_estimates_table(fit_simple_normal, 
                                                         transformations = list("m" = list(fun = exp)))
  
  # ============================================================================
  # Generate 14-16: Empty tables
  # ============================================================================
  runjags_summary_empty <- runjags_estimates_empty_table()
  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  ensemble_inference_empty <- ensemble_inference_empty_table()
  
  # ============================================================================
  # Create list of all tables to generate
  # ============================================================================
  fits <- list(
    # 1-6: Weightfunction models (basic test)
    model_summary_table(models_wf[[2]]),
    runjags_estimates_table(fit1),
    ensemble_estimates_table(mixed_posteriors_wf, parameters = c("m", "omega"), probs = c(.025, 0.95)),
    ensemble_inference_table(inference_wf, names(inference_wf)),
    ensemble_summary_table(models_wf, c("m", "omega")),
    ensemble_diagnostics_table(models_wf, c("m", "omega")),
    
    # 7-11: Simple models
    model_summary_table(models_simple[[1]]),
    runjags_estimates_table(fit_simple_normal),
    ensemble_inference_table(inference_simple, names(inference_simple)),
    ensemble_summary_table(models_simple, c("m", "s")),
    ensemble_diagnostics_table(models_simple, c("m", "s")),
    
    # 12-13: Transformations
    runjags_summary_transform,
    runjags_estimates_table(fit2, transformations = list("m" = list(fun = exp))),
    
    # 14-16: Empty tables
    runjags_summary_empty,
    ensemble_estimates_empty,
    ensemble_inference_empty
  )
  
  # Create print directory if it doesn't exist
  print_dir <- testthat::test_path("..", "results", "print")
  if (!dir.exists(print_dir)) {
    dir.create(print_dir, recursive = TRUE)
    message("Created directory: ", print_dir)
  }
  
  # Generate print files
  for(i in seq_along(fits)){
    output_lines <- capture_output_lines(fits[[i]], print = TRUE, width = 150)
    output_file <- file.path(print_dir, paste0(i, ".txt"))
    writeLines(output_lines, output_file)
    message("Generated reference file: ", output_file)
  }
  
  message("\nReference print files generated successfully in ", normalizePath(print_dir))
  message("Total files created: ", length(fits))
}
