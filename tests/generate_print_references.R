#!/usr/bin/env Rscript
# Script to generate reference print output files for test-summary-tables.R
# Run this script whenever the print format changes or to initialize reference files

# Load required packages
library(BayesTools)
library(testthat)

# Silence JAGS output
runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

# Get the directory where prefitted models are stored
temp_fits_dir <- Sys.getenv("BAYESTOOLS_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  temp_fits_dir <- file.path(tempdir(), "BayesTools_test_fits")
  Sys.setenv(BAYESTOOLS_TEST_FITS_DIR = temp_fits_dir)
}

# If models don't exist, run test-00-model-fits.R first
if (!dir.exists(temp_fits_dir) || length(list.files(temp_fits_dir, pattern = "\\.RDS$")) == 0) {
  message("Pre-fitted models not found. Running test-00-model-fits.R...")
  # Look for test file in testthat subdirectory
  test_file_path <- "testthat/test-00-model-fits.R"
  if (file.exists(test_file_path)) {
    test_file(test_file_path)
  } else {
    stop("Could not find test-00-model-fits.R. Please run it manually first or run this script from the tests/ directory.")
  }
}

# Load models for print testing
fit0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
marglik0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))

fit1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
marglik1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))

fit2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
marglik2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))

# Create models list
models <- list(
  list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
  list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
  list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
)
models <- models_inference(models)

# Generate inference and posteriors
inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), 
                                 is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), 
                                    is_null_list = list("m" = 0, "omega" = 1), seed = 1)

# Create list of table objects to print
fits <- list(
  model_summary_table(models[[2]]),
  runjags_estimates_table(fit1),
  ensemble_estimates_table(mixed_posteriors, parameters = c("m", "omega"), probs = c(.025, 0.95)),
  ensemble_inference_table(inference, names(inference)),
  ensemble_summary_table(models, c("m", "omega")),
  ensemble_diagnostics_table(models, c("m", "omega"))
)

# Create print directory if it doesn't exist
print_dir <- "results/print"
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
