context("Summary tables functions")

# This file tests summary table functions using pre-fitted models from test-00-model-fits.R
# The pre-fitted models are loaded from the temporary directory set by test-00-model-fits.R

skip_on_cran()

# Set UPDATE_OUTPUT to TRUE to regenerate the print output files
# This flag should be set to TRUE when you want to regenerate reference files
# after making changes to summary table functions or print methods
UPDATE_OUTPUT <- FALSE

# Get the directory where pre-fitted models are stored
temp_fits_dir <- Sys.getenv("BAYESTOOLS_TEST_FITS_DIR")
if(temp_fits_dir == "" || !dir.exists(temp_fits_dir)){
  skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
}

# Load model registry to know which models are available
model_registry_file <- file.path(temp_fits_dir, "model_registry.RDS")
if(!file.exists(model_registry_file)){
  skip("Model registry not found. Run test-00-model-fits.R first.")
}
model_registry <- readRDS(model_registry_file)

# Helper function to load a pre-fitted model
load_fit <- function(model_name) {
  # Check if model is in registry
  if(!model_name %in% model_registry$model_name) {
    stop(paste("Model", model_name, "not found in registry. Available models:",
               paste(head(model_registry$model_name, 10), collapse=", ")))
  }
  
  fit_file <- file.path(temp_fits_dir, paste0(model_name, ".RDS"))
  if(!file.exists(fit_file)){
    stop(paste("Pre-fitted model file", model_name, "not found at", fit_file,
               "\nPlease run test-00-model-fits.R first to generate required models."))
  }
  fit <- readRDS(fit_file)
  # Check if it's a valid runjags fit (not an error)
  if(!inherits(fit, "runjags")){
    stop(paste("Model", model_name, "is not a valid runjags fit.",
               "\nThis may indicate an error during model fitting in test-00-model-fits.R.",
               "\nPlease check that test file for issues."))
  }
  fit
}

# =============================================================================
# TEST 1: Basic summary table functions with formula models
# =============================================================================
test_that("Summary tables basic functions work", {
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load formula-based models
  fit1 <- load_fit("fit_formula_simple")
  fit2 <- load_fit("fit_formula_treatment")
  fit3 <- load_fit("fit_formula_orthonormal")
  
  # Create dummy marginal likelihoods for model comparison
  marglik1 <- list(logml = -88.22)
  marglik2 <- list(logml = -88.66)
  marglik3 <- list(logml = -89.89)
  class(marglik1) <- class(marglik2) <- class(marglik3) <- "bridge"
  
  # Create model list
  models <- list(
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2)),
    list(fit = fit3, marglik = marglik3, prior_weights = 1, fit_summary = runjags_estimates_table(fit3))
  )
  models <- models_inference(models)
  
  # Test model summary
  model_summary <- model_summary_table(models[[1]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  "))
  expect_true(nrow(model_summary) >= 5)
  
  # Test runjags summary structure
  runjags_summary <- models[[1]]$fit_summary
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_true(nrow(runjags_summary) > 0)
  
  # Test empty tables
  runjags_summary_empty <- runjags_estimates_empty_table()
  expect_equivalent(nrow(runjags_summary_empty), 0)
  expect_equal(colnames(runjags_summary_empty), colnames(runjags_summary))
  
  model_summary_empty <- model_summary_empty_table()
  expect_equivalent(nrow(model_summary_empty), 5)
})

# =============================================================================
# TEST 2: Ensemble summary functions
# =============================================================================
test_that("Ensemble summary functions work", {
  
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load models
  fit1 <- load_fit("fit_formula_simple")
  fit2 <- load_fit("fit_formula_treatment") 
  fit3 <- load_fit("fit_formula_orthonormal")
  
  # Create dummy marginal likelihoods
  marglik1 <- list(logml = -88.22)
  marglik2 <- list(logml = -88.66)
  marglik3 <- list(logml = -89.89)
  class(marglik1) <- class(marglik2) <- class(marglik3) <- "bridge"
  
  # Create model list
  models <- list(
    list(fit = fit1, marglik = marglik1, fit_summary = runjags_estimates_table(fit1), prior_weights = 1),
    list(fit = fit2, marglik = marglik2, fit_summary = runjags_estimates_table(fit2), prior_weights = 1),
    list(fit = fit3, marglik = marglik3, fit_summary = runjags_estimates_table(fit3), prior_weights = 1)
  )
  models <- models_inference(models)
  
  # Compute ensemble inference
  inference <- ensemble_inference(
    model_list = models,
    parameters = c("mu_x_cont1"),
    is_null_list = list("mu_x_cont1" = c(TRUE, FALSE, TRUE)),
    conditional = FALSE
  )
  
  # Mix posteriors
  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = c("mu_x_cont1"),
    is_null_list = list("mu_x_cont1" = c(TRUE, FALSE, TRUE)),
    seed = 1, n_samples = 1000
  )
  
  # Test ensemble estimates
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("mu_x_cont1"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025", "0.95"))
  expect_true(nrow(estimates_table) > 0)
  
  # Test ensemble inference
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_true(nrow(inference_table) > 0)
  
  # Test ensemble summary
  summary_table <- ensemble_summary_table(models, c("mu_x_cont1"))
  expect_true("(mu) x_cont1" %in% colnames(summary_table))
  
  # Test ensemble diagnostics
  diagnostics_table <- ensemble_diagnostics_table(models, c("mu_x_cont1"))
  expect_true("(mu) x_cont1" %in% colnames(diagnostics_table))
  
  # Test empty tables
  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  expect_equivalent(nrow(ensemble_estimates_empty), 0)
  expect_equal(colnames(ensemble_estimates_empty), colnames(estimates_table))
  
  ensemble_inference_empty <- ensemble_inference_empty_table()
  expect_equivalent(nrow(ensemble_inference_empty), 0)
  expect_equal(colnames(ensemble_inference_empty), colnames(inference_table))
})

# =============================================================================
# TEST 3: Formula prefix and transformations
# =============================================================================
test_that("Formula prefix removal and transformations work", {
  
  skip_on_os(c("mac", "linux", "solaris"))
  
  fit <- load_fit("fit_formula_treatment")
  
  # Test with formula prefix
  runjags_summary1 <- runjags_estimates_table(fit, formula_prefix = TRUE)
  expect_true(any(grepl("^\\(mu\\)", rownames(runjags_summary1))))
  
  # Test without formula prefix
  runjags_summary2 <- runjags_estimates_table(fit, formula_prefix = FALSE)
  expect_true(all(!grepl("^\\(mu\\)", rownames(runjags_summary2))))
  
  # Test numeric transformation
  fit_simple <- load_fit("fit_simple_normal")
  summary_transformed <- runjags_estimates_table(
    fit_simple,
    transformations = list("m" = list(fun = exp))
  )
  expect_true(nrow(summary_transformed) > 0)
  
  # Test factor transformation
  fit_orth <- load_fit("fit_formula_orthonormal")
  summary_factor <- suppressMessages(
    runjags_estimates_table(fit_orth, transform_factors = TRUE)
  )
  expect_true(any(grepl("\\[dif:", rownames(summary_factor))))
})

# =============================================================================
# TEST 4: Interaction models
# =============================================================================
test_that("Summary tables work with interaction models", {
  
  skip_on_os(c("mac", "linux", "solaris"))
  
  # Load interaction models
  fit_cont <- load_fit("fit_formula_interaction_cont")
  fit_mix <- load_fit("fit_formula_interaction_mix")
  fit_fac <- load_fit("fit_formula_interaction_fac")
  
  # Test estimates for interactions
  est_cont <- runjags_estimates_table(fit_cont)
  expect_true(any(grepl(":", rownames(est_cont))))
  
  est_mix <- runjags_estimates_table(fit_mix)
  expect_true(any(grepl(":", rownames(est_mix))))
  
  est_fac <- runjags_estimates_table(fit_fac)
  expect_true(any(grepl(":", rownames(est_fac))))
})

# =============================================================================
# TEST 5: Random effects models
# =============================================================================
test_that("Summary tables work with random effects models", {
  
  skip_on_os(c("mac", "linux", "solaris"))
  
  # Load random effects models
  fit_ri <- load_fit("fit_random_intercept")
  fit_rs <- load_fit("fit_random_slope")
  fit_rfs <- load_fit("fit_random_factor_slope")
  
  # Test estimates for random effects
  est_ri <- runjags_estimates_table(fit_ri)
  expect_true(any(grepl("\\|", rownames(est_ri))))
  
  est_rs <- runjags_estimates_table(fit_rs)
  expect_true(any(grepl("\\|", rownames(est_rs))))
  
  est_rfs <- runjags_estimates_table(fit_rfs)
  expect_true(any(grepl("\\|", rownames(est_rfs))))
})

# =============================================================================
# TEST 6: Joint complex models
# =============================================================================
test_that("Summary tables work with complex joint models", {
  
  skip_on_os(c("mac", "linux", "solaris"))
  
  # Load joint complex model
  fit <- load_fit("fit_joint_complex")
  
  # Test estimates
  est <- runjags_estimates_table(fit)
  expect_true(nrow(est) > 0)
  expect_true(any(grepl("inclusion", rownames(est)))) # mixture/spike-and-slab components
})

# =============================================================================
# TEST 7: Advanced JAGS features
# =============================================================================
test_that("Summary tables work with advanced JAGS features", {
  
  # Load models with advanced features
  fit_add_param <- load_fit("fit_add_parameters")
  fit_autofit <- load_fit("fit_autofit_error")
  fit_parallel <- load_fit("fit_parallel")
  
  # Test add_parameters
  est_add <- runjags_estimates_table(fit_add_param)
  expect_true("g" %in% rownames(est_add)) # Additional parameter
  
  # Test autofit
  est_autofit <- runjags_estimates_table(fit_autofit)
  expect_true(nrow(est_autofit) > 0)
  summary_autofit <- summary(fit_autofit)
  expect_true(summary_autofit[1,"MCerr"] < 0.05) # Should meet convergence criterion
  
  # Test parallel
  est_parallel <- runjags_estimates_table(fit_parallel)
  expect_true(nrow(est_parallel) > 0)
  expect_equal(length(fit_parallel$mcmc), 2) # 2 chains
})

# =============================================================================
# TEST 8: Column addition and removal
# =============================================================================
test_that("Column addition and removal functions work", {
  
  fit <- load_fit("fit_simple_normal")
  runjags_summary <- runjags_estimates_table(fit)
  
  # Test adding columns - error cases
  expect_error(
    add_column(runjags_summary, column_title = "New Title",
               column_values = c(0.2, 0.3, 0.4, 0.5)),
    "The 'column_values' must be a vector of the same length as has the table rows."
  )
  
  expect_error(
    add_column(data.frame(a = 1:3, b = c("A", "B", "C")),
               column_title = "New Title", column_values = c(0.2, 0.3, 0.4)),
    "The 'table' must be of class 'BayesTools_table'."
  )
  
  # Test adding columns - success
  new_table <- add_column(runjags_summary, column_title = "Test",
                         column_values = rep(0.5, nrow(runjags_summary)))
  expect_true("Test" %in% colnames(new_table))
  expect_equal(ncol(new_table), ncol(runjags_summary) + 1)
  
  # Test removing columns - error cases
  expect_error(
    remove_column(runjags_summary, column_position = 100),
    "The 'column_position' must be equal or lower than"
  )
  
  # Test removing columns - success
  removed_table <- remove_column(runjags_summary, column_position = 1)
  expect_equal(ncol(removed_table), ncol(runjags_summary) - 1)
})

# =============================================================================
# TEST 9: Print output comparison with saved files
# =============================================================================
test_that("Print output matches saved files", {
  
  # Create list of fits to test
  fits <- list(
    fit_simple_normal = load_fit("fit_simple_normal"),
    fit_formula_simple = load_fit("fit_formula_simple"),
    fit_formula_treatment = load_fit("fit_formula_treatment"),
    fit_formula_orthonormal = load_fit("fit_formula_orthonormal"),
    fit_formula_interaction_cont = load_fit("fit_formula_interaction_cont"),
    fit_random_intercept = load_fit("fit_random_intercept")
  )
  
  # Check if print output files exist
  # Using relative path from tests/testthat/ directory
  print_dir <- file.path("../results/print")
  if(!dir.exists(print_dir)) {
    skip("Print output directory not found. Set UPDATE_OUTPUT=TRUE to generate outputs.")
  }
  
  # Compare each fit's print output with saved file
  for(i in seq_along(fits)) {
    print_file <- file.path(print_dir, paste0(i, ".txt"))
    if(file.exists(print_file)) {
      saved_output <- readLines(print_file)
      
      current_output <- capture_output_lines(fits[[i]], print = TRUE, width = 150)
      
      expect_equal(current_output, saved_output,
                   info = paste("Model", names(fits)[i], "print output mismatch"))
    }
  }
})

# =============================================================================
# UTILITY: Generate print output files (UPDATE_OUTPUT flag at top of file)
# =============================================================================
# The UPDATE_OUTPUT flag at the top of this file controls whether to regenerate
# the print output reference files. Set it to TRUE to regenerate outputs.

if(UPDATE_OUTPUT) {
  
  test_that("Generate print output files", {
    
    # Create print directory if it doesn't exist
    # Using relative path from tests/testthat/ directory
    print_dir <- file.path("../results/print")
    dir.create(print_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Select models to generate print outputs for
    fits <- list(
      load_fit("fit_simple_normal"),
      load_fit("fit_formula_simple"),
      load_fit("fit_formula_treatment"),
      load_fit("fit_formula_orthonormal"),
      load_fit("fit_formula_interaction_cont"),
      load_fit("fit_random_intercept")
    )
    
    # Generate print files
    for(i in seq_along(fits)) {
      output_lines <- capture_output_lines(fits[[i]], print = TRUE, width = 150)
      writeLines(output_lines, con = file.path(print_dir, paste0(i, ".txt")))
    }
    
    message("Print output files generated successfully in ", print_dir)
  })
}
