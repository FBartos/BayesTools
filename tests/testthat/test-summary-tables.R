context("Summary tables functions")

# Get the directory where prefitted models are stored
temp_fits_dir <- Sys.getenv("BAYESTOOLS_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  temp_fits_dir <- file.path(tempdir(), "BayesTools_test_fits")
}

# Skip tests on CRAN as they require pre-fitted models
skip_on_cran()

test_that("Summary tables work with simple priors", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load pre-fitted models from test-00-model-fits.R
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))
  
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))
  
  # Create models list for model averaging
  models <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_simple_normal)),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_simple_spike))
  )
  models <- models_inference(models)
  
  # Test model_summary_table
  model_summary <- model_summary_table(models[[1]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  "))
  expect_s3_class(model_summary, "BayesTools_table")
  
  # Test runjags_estimates_table
  runjags_summary <- runjags_estimates_table(fit_simple_normal)
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary), c("m", "s"))
  expect_s3_class(runjags_summary, "BayesTools_table")
  
  # Test ensemble_inference
  inference <- ensemble_inference(model_list = models, parameters = c("m", "s"), 
                                   is_null_list = list("m" = 0, "s" = 1), conditional = FALSE)
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_s3_class(inference_table, "BayesTools_table")
  
  # Test ensemble_summary_table
  summary_table <- ensemble_summary_table(models, c("m", "s"))
  expect_equal(colnames(summary_table), c("Model", "m", "s", "prior_prob", "marglik", "post_prob", "inclusion_BF"))
  expect_s3_class(summary_table, "BayesTools_table")
  
  # Test ensemble_diagnostics_table
  diagnostics_table <- ensemble_diagnostics_table(models, c("m", "s"))
  expect_true("max_MCMC_error" %in% colnames(diagnostics_table))
  expect_s3_class(diagnostics_table, "BayesTools_table")
})

test_that("Summary tables work with weightfunction priors", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load pre-fitted models for summary tables (from SECTION 1B)
  fit0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))
  
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))
  
  fit2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))
  
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
  )
  models <- models_inference(models)
  
  # Test with weightfunction parameters
  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), 
                                   is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), 
                                      is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  
  # Test ensemble_estimates_table
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("m", "omega"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_true("m" %in% rownames(estimates_table))
  expect_s3_class(estimates_table, "BayesTools_table")
})

test_that("Summary tables work with vector priors", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load vector prior models
  fit_vector_mnormal <- readRDS(file.path(temp_fits_dir, "fit_vector_mnormal.RDS"))
  fit_vector_mcauchy <- readRDS(file.path(temp_fits_dir, "fit_vector_mcauchy.RDS"))
  
  # Test runjags_estimates_table with multivariate parameters
  runjags_summary_mn <- runjags_estimates_table(fit_vector_mnormal)
  expect_true(any(grepl("p1\\[", rownames(runjags_summary_mn))))
  expect_s3_class(runjags_summary_mn, "BayesTools_table")
  
  runjags_summary_mc <- runjags_estimates_table(fit_vector_mcauchy)
  expect_true(any(grepl("p1\\[", rownames(runjags_summary_mc))))
  expect_s3_class(runjags_summary_mc, "BayesTools_table")
})

test_that("Summary tables work with factor priors", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load factor prior models
  fit_factor_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS"))
  fit_factor_treatment <- readRDS(file.path(temp_fits_dir, "fit_factor_treatment.RDS"))
  fit_factor_meandif <- readRDS(file.path(temp_fits_dir, "fit_factor_meandif.RDS"))
  
  # Test runjags_estimates_table with factors
  runjags_summary_orth <- runjags_estimates_table(fit_factor_orthonormal)
  expect_s3_class(runjags_summary_orth, "BayesTools_table")
  
  runjags_summary_treat <- runjags_estimates_table(fit_factor_treatment)
  expect_s3_class(runjags_summary_treat, "BayesTools_table")
  
  runjags_summary_md <- runjags_estimates_table(fit_factor_meandif)
  expect_s3_class(runjags_summary_md, "BayesTools_table")
})

test_that("Summary tables work with mixture priors", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load mixture prior models
  fit_mixture_simple <- readRDS(file.path(temp_fits_dir, "fit_mixture_simple.RDS"))
  fit_mixture_components <- readRDS(file.path(temp_fits_dir, "fit_mixture_components.RDS"))
  
  # Test runjags_estimates_table with mixture priors
  runjags_summary_mix <- runjags_estimates_table(fit_mixture_simple)
  expect_s3_class(runjags_summary_mix, "BayesTools_table")
  expect_true(any(grepl("inclusion", rownames(runjags_summary_mix))))
  
  runjags_summary_comp <- runjags_estimates_table(fit_mixture_components)
  expect_s3_class(runjags_summary_comp, "BayesTools_table")
})

test_that("Summary tables work with spike-and-slab priors", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load spike-and-slab models
  fit_spike_slab_simple <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_simple.RDS"))
  
  # Test runjags_estimates_table with spike-and-slab
  runjags_summary_ss <- runjags_estimates_table(fit_spike_slab_simple)
  expect_s3_class(runjags_summary_ss, "BayesTools_table")
  expect_true(any(grepl("inclusion", rownames(runjags_summary_ss))))
  
  # Test runjags_inference_table
  inference_ss <- runjags_inference_table(fit_spike_slab_simple)
  expect_s3_class(inference_ss, "BayesTools_table")
  expect_equal(colnames(inference_ss), c("prior_prob", "post_prob", "inclusion_BF"))
})

test_that("Summary tables work with formula-based models", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load formula-based models
  fit_formula_simple <- readRDS(file.path(temp_fits_dir, "fit_formula_simple.RDS"))
  marglik_formula_simple <- readRDS(file.path(temp_fits_dir, "fit_formula_simple_marglik.RDS"))
  
  fit_formula_treatment <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  marglik_formula_treatment <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment_marglik.RDS"))
  
  # Test runjags_estimates_table with formulas
  runjags_summary_formula <- runjags_estimates_table(fit_formula_simple)
  expect_s3_class(runjags_summary_formula, "BayesTools_table")
  expect_true(any(grepl("\\(mu\\)", rownames(runjags_summary_formula))))
  
  # Test with formula_prefix = FALSE
  runjags_summary_no_prefix <- runjags_estimates_table(fit_formula_simple, formula_prefix = FALSE)
  expect_s3_class(runjags_summary_no_prefix, "BayesTools_table")
  expect_false(any(grepl("\\(mu\\)", rownames(runjags_summary_no_prefix))))
  
  # Test model averaging with formula models
  models <- list(
    list(fit = fit_formula_simple, marglik = marglik_formula_simple, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_formula_simple)),
    list(fit = fit_formula_treatment, marglik = marglik_formula_treatment, prior_weights = 1, 
         fit_summary = runjags_estimates_table(fit_formula_treatment))
  )
  models <- models_inference(models)
  
  summary_table <- ensemble_summary_table(models, c("mu_intercept", "mu_x_cont1", "sigma"))
  expect_s3_class(summary_table, "BayesTools_table")
})

test_that("Summary tables work with interaction models", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load interaction models
  fit_formula_interaction_cont <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))
  fit_formula_interaction_mix <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix.RDS"))
  
  # Test runjags_estimates_table with interactions
  runjags_summary_int <- runjags_estimates_table(fit_formula_interaction_cont)
  expect_s3_class(runjags_summary_int, "BayesTools_table")
  expect_true(any(grepl(":", rownames(runjags_summary_int))))
  
  runjags_summary_mix_int <- runjags_estimates_table(fit_formula_interaction_mix)
  expect_s3_class(runjags_summary_mix_int, "BayesTools_table")
})

test_that("Summary tables work with random effects models", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load random effects models
  fit_random_intercept <- readRDS(file.path(temp_fits_dir, "fit_random_intercept.RDS"))
  fit_random_slope <- readRDS(file.path(temp_fits_dir, "fit_random_slope.RDS"))
  
  # Test runjags_estimates_table with random effects
  runjags_summary_re_int <- runjags_estimates_table(fit_random_intercept)
  expect_s3_class(runjags_summary_re_int, "BayesTools_table")
  expect_true(any(grepl("\\|", rownames(runjags_summary_re_int))))
  
  runjags_summary_re_slope <- runjags_estimates_table(fit_random_slope)
  expect_s3_class(runjags_summary_re_slope, "BayesTools_table")
})

test_that("Summary tables work with expression priors", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load expression prior models
  fit_expression_simple <- readRDS(file.path(temp_fits_dir, "fit_expression_simple.RDS"))
  fit_expression_spike_slab <- readRDS(file.path(temp_fits_dir, "fit_expression_spike_slab.RDS"))
  
  # Test runjags_estimates_table with expression priors
  runjags_summary_expr <- runjags_estimates_table(fit_expression_simple)
  expect_s3_class(runjags_summary_expr, "BayesTools_table")
  
  runjags_summary_expr_ss <- runjags_estimates_table(fit_expression_spike_slab)
  expect_s3_class(runjags_summary_expr_ss, "BayesTools_table")
})

test_that("Summary tables work with advanced JAGS features", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load advanced feature models
  fit_add_parameters <- readRDS(file.path(temp_fits_dir, "fit_add_parameters.RDS"))
  fit_autofit_error <- readRDS(file.path(temp_fits_dir, "fit_autofit_error.RDS"))
  fit_parallel <- readRDS(file.path(temp_fits_dir, "fit_parallel.RDS"))
  
  # Test with add_parameters
  runjags_summary_add <- runjags_estimates_table(fit_add_parameters)
  expect_s3_class(runjags_summary_add, "BayesTools_table")
  expect_true("g" %in% rownames(runjags_summary_add))
  
  # Test with autofit
  runjags_summary_autofit <- runjags_estimates_table(fit_autofit_error)
  expect_s3_class(runjags_summary_autofit, "BayesTools_table")
  
  # Test with parallel
  runjags_summary_parallel <- runjags_estimates_table(fit_parallel)
  expect_s3_class(runjags_summary_parallel, "BayesTools_table")
})

test_that("Summary tables transformations work", {
  
  skip_if_not_installed("rjags")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load model for transformation testing
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  
  # Test transformations parameter
  runjags_summary_transform <- runjags_estimates_table(fit_simple_normal, 
                                                        transformations = list("m" = list(fun = exp)))
  expect_s3_class(runjags_summary_transform, "BayesTools_table")
  
  # Load factor model for factor transformation
  fit_factor_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS"))
  
  # Test transform_factors parameter
  runjags_summary_tf <- suppressMessages(runjags_estimates_table(fit_factor_orthonormal, 
                                                                   transform_factors = TRUE))
  expect_s3_class(runjags_summary_tf, "BayesTools_table")
})

test_that("Empty summary tables work", {
  
  # Test all empty table functions
  runjags_summary_empty <- runjags_estimates_empty_table()
  expect_s3_class(runjags_summary_empty, "BayesTools_table")
  expect_equivalent(nrow(runjags_summary_empty), 0)
  
  runjags_inference_empty <- runjags_inference_empty_table()
  expect_s3_class(runjags_inference_empty, "BayesTools_table")
  expect_equivalent(nrow(runjags_inference_empty), 0)
  
  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  expect_s3_class(ensemble_estimates_empty, "BayesTools_table")
  expect_equivalent(nrow(ensemble_estimates_empty), 0)
  
  ensemble_inference_empty <- ensemble_inference_empty_table()
  expect_s3_class(ensemble_inference_empty, "BayesTools_table")
  expect_equivalent(nrow(ensemble_inference_empty), 0)
  
  ensemble_summary_empty <- ensemble_summary_empty_table()
  expect_s3_class(ensemble_summary_empty, "BayesTools_table")
  expect_equivalent(nrow(ensemble_summary_empty), 0)
  
  ensemble_diagnostics_empty <- ensemble_diagnostics_empty_table()
  expect_s3_class(ensemble_diagnostics_empty, "BayesTools_table")
  expect_equivalent(nrow(ensemble_diagnostics_empty), 0)
  
  model_summary_empty <- model_summary_empty_table()
  expect_s3_class(model_summary_empty, "BayesTools_table")
})

test_that("Summary table printing works", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  if (!dir.exists(temp_fits_dir)) {
    skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
  }
  
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  
  # Load models for print testing
  fit0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))
  
  fit1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))
  
  fit2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))
  
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
  )
  models <- models_inference(models)
  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), 
                                   is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), 
                                      is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  
  # Create list of table objects to test printing
  fits <- list(
    model_summary_table(models[[2]]),
    runjags_estimates_table(fit1),
    ensemble_estimates_table(mixed_posteriors, parameters = c("m", "omega"), probs = c(.025, 0.95)),
    ensemble_inference_table(inference, names(inference)),
    ensemble_summary_table(models, c("m", "omega")),
    ensemble_diagnostics_table(models, c("m", "omega"))
  )
  
  # Check if reference files exist
  # test_path() constructs paths relative to tests/testthat/, so we use ".." to reach tests/results/print
  print_dir <- testthat::test_path("..", "results", "print")
  if (!dir.exists(print_dir)) {
    skip("Print reference directory not found. Run tests/generate_print_references.R to generate.")
  }
  
  # Compare printed output with saved reference files
  for(i in 1:length(fits)){
    ref_file <- file.path(print_dir, paste0(i, ".txt"))
    if (!file.exists(ref_file)) {
      skip(paste0("Reference file ", i, ".txt not found. Run tests/generate_print_references.R to generate."))
    }
    
    # Use readLines for simpler and more robust file reading
    expected_output <- readLines(ref_file, warn = FALSE)
    actual_output <- capture_output_lines(fits[[i]], print = TRUE, width = 150)
    
    expect_equal(actual_output, expected_output, 
                 info = paste0("Print output mismatch for table ", i))
  }
})


