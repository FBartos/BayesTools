---
applyTo: "**/tests/testthat/*.R"
description: Guidelines for organizing and maintaining tests in BayesTools, including model fitting, model averaging, and summary table tests. Ensures consistency and avoids duplication.
---

# BayesTools Test Organization Guidelines

## Overview

Tests in BayesTools follow a structured organization where model fitting is centralized in `test-00-model-fits.R` and consumed by other test files. This approach ensures consistency, avoids duplication, and speeds up test execution.

## Key Principles

### 1. Single Source of Truth for Model Fitting

**All model fitting and marginal likelihood computation must be done in `test-00-model-fits.R`.**

- `test-00-model-fits.R` is the **only** file that should:
  - Fit JAGS models using `JAGS_fit()`
  - Compute marginal likelihoods using `JAGS_bridgesampling()`
  - Save fitted models as RDS files
  - Save marginal likelihoods as separate RDS files

- Other test files should:
  - **Only load** pre-fitted models using `readRDS()`
  - **Only load** pre-computed marginal likelihoods using `readRDS()`
  - Test the functionality they are designed for (e.g., model averaging, plotting, etc.)

### 2. Avoid Duplication

**Before adding a new model to `test-00-model-fits.R`, check if a similar model already exists.**

Models are duplicates if they have the same model structure, prior types, and data structure. Use one model per prior type.

### 3. Model Naming Convention

Use pattern: `fit_{category}_{descriptor}` (e.g., `fit_simple_normal`, `fit_formula_treatment`)

**Model Registry**: `test-00-model-fits.R` maintains a registry of all fitted models in `model_registry.RDS`. Other test files should load this registry to discover available models rather than hardcoding model names:

```r
registry_file <- file.path(temp_fits_dir, "model_registry.RDS")
model_registry <- readRDS(registry_file)
model_names <- model_registry$model_name
```

### 4. Saving and Loading Models

```r
# In test-00-model-fits.R: Save with save_fit() helper
result <- save_fit(fit_model_name, "fit_model_name",
                   marglik = marglik_model_name,  # If available
                   simple_priors = TRUE, note = "Description")

# In other test files: Load with readRDS()
fit_model_name <- readRDS(file.path(temp_fits_dir, "fit_model_name.RDS"))
marglik_model_name <- readRDS(file.path(temp_fits_dir, "fit_model_name_marglik.RDS"))
```

**Note**: Marginal likelihoods are only computed for models with actual data (not spike-and-slab or mixture priors).

### 5. Helper Functions for Reference Files

**Use a helper function to reduce repetition** when saving/testing reference files:

```r
# Define at the top of test files with reference outputs
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
```

### 6. Test File Organization

#### test-00-model-fits.R
- **Purpose**: Fit all models and compute marginal likelihoods
- **Run order**: First (prefix `00-`)
- **Outputs**: RDS files in `tempdir()` via `BAYESTOOLS_TEST_FITS_DIR` env var
- **Registry**: Maintains `model_registry` with model metadata

#### Tests Using JAGS Models
All tests that use JAGS models (e.g., `test-model-averaging.R`, `test-JAGS-*.R`, `test-summary-tables.R`) must:
- Load pre-fitted models from `temp_fits_dir` using `readRDS()`
- **Never** fit models directly (only `test-00-model-fits.R` fits models)
- Check model availability with `if (!dir.exists(temp_fits_dir))` and skip appropriately
- Use `skip_if_not_installed("rjags")` and `skip_if_not_installed("bridgesampling")`

**For tests with reference files** (e.g., `test-summary-tables.R`, visual regression tests):
- **Configuration**: `GENERATE_REFERENCE_FILES` flag (FALSE = test, TRUE = generate)
  - **IMPORTANT**: **Never** modify this flag. Only the package maintainer changes this flag when intentionally updating reference files after format changes.
  - Default value is `FALSE` (testing mode)
  - Changing to `TRUE` regenerates all reference files (tables, figures, etc.) and should only be done by the maintainer
- **Outputs**: Reference files (`.txt`, `.svg`, `.png`, etc.) stored in `tests/results/` subdirectories

## Maintenance Checklist

**Adding a new model:**
- [ ] Check for duplicates in `test-00-model-fits.R`
- [ ] Add model to `test-00-model-fits.R` with `save_fit()` and appropriate metadata
**Using pre-fitted models:**
- [ ] Load with `readRDS()`, never fit models outside `test-00-model-fits.R`
- [ ] Add skip conditions for missing models/packages
- [ ] Check marginal likelihood file existence before loading

**Updating summary table tests (MAINTAINER ONLY):**
- [ ] Set `GENERATE_REFERENCE_FILES <- TRUE` in `test-summary-tables.R`
- [ ] Run tests to generate reference files
- [ ] Review diffs carefully before committing
- [ ] Reset flag to `FALSE`
- **Note**: Contributors/agents should **never** modify `GENERATE_REFERENCE_FILES`CE_FILES <- TRUE` in `test-summary-tables.R`
- [ ] Run tests to generate reference files
- [ ] Review diffs carefully before committing
- [ ] Reset flag to `FALSE`

## Quick Examples

### Adding and Using a Model

```r
# 1. In test-00-model-fits.R
fit_new     <- JAGS_fit(model_syntax, data, priors, ...)
marglik_new <- JAGS_bridgesampling(fit_new, log_posterior, data, priors)
result      <- save_fit(fit_new, "fit_new", marglik = marglik_new, note = "Description")

# 2. In any test file using JAGS models
fit_new      <- readRDS(file.path(temp_fits_dir, "fit_new.RDS"))
marglik_file <- file.path(temp_fits_dir, "fit_new_marglik.RDS")
if (file.exists(marglik_file)) {
  marglik_new <- readRDS(marglik_file)
}

# 3. Add to test-summary-tables.R model_names vector
model_names <- c(..., "fit_new")
```
## Common Pitfalls

❌ Fitting models outside `test-00-model-fits.R`  
❌ Creating duplicate models with different parameters  
❌ **Modifying `GENERATE_REFERENCE_FILES` flag** (maintainer only)

✅ Always load pre-fitted models with `readRDS()`  
✅ Use one model per prior type  
✅ Leave `GENERATE_REFERENCE_FILES <- FALSE` unchanged
✅ Set `GENERATE_REFERENCE_FILES <- TRUE` when updating formats
## Troubleshooting

- **"Pre-fitted models not available"**: Run `devtools::test(filter = "00-model-fits")`
- **Summary table mismatch**: Contact maintainer; **do not** modify `GENERATE_REFERENCE_FILES`
- **Marginal likelihood not found**: Check model has data and isn't spike-and-slab/mixture
- **Marginal likelihood not found**: Check model has data and isn't spike-and-slab/mixture