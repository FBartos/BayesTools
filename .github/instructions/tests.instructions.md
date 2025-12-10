---
applyTo: "**/tests/testthat/*.R"
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

### 2. STRICTLY Avoid Duplication

**Before adding a new model to `test-00-model-fits.R`, you MUST exhaustively check if an existing model can be used.**

- **Do not create a new model just to test a specific function** (e.g., a plot or summary). Use an existing model that has the necessary components (e.g., if you need a model with a factor prior, use `fit_factor_independent` or `fit_formula_interaction_fac`).
- **Models are duplicates** if they have the same model structure, prior types, and data structure.
- **Reuse Strategy**:
    1. Read `test-00-model-fits.R` to see available models.
    2. Identify a model that has the features you need (e.g., "I need a model with a spike-and-slab prior").
    3. Use that model in your test.
    4. **Only** if no such model exists, add a new one to `test-00-model-fits.R`.

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

# Define at the top of test files with reference outputs
# Load common test helpers that define test_reference_table() and test_reference_text()
REFERENCE_DIR <<- testthat::test_path("..", "results", "print")
source(testthat::test_path("common-functions.R"))

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

## AI Agent Protocol

When asked to write or refactor tests:

1.  **Scan `test-00-model-fits.R` FIRST.** Understand the inventory of available models.
2.  **Map requirements to existing models.** If the user needs a test for "diagnostic plots for factor priors", find an existing model with factor priors (e.g., `fit_formula_interaction_fac`).
3.  **Refuse to create new models** unless the test requires a specific mathematical structure not present in the entire suite.
4.  **Never** add a model to `test-00-model-fits.R` without explicitly explaining why none of the existing 15+ models suffice.

## Maintenance Checklist

**Adding a new model:**
- [ ] Check for duplicates in `test-00-model-fits.R`
- [ ] Add model to `test-00-model-fits.R` with `save_fit()` and appropriate metadata
**Using pre-fitted models:**
- [ ] Load with `readRDS()`, never fit models outside `test-00-model-fits.R`
- [ ] Add skip conditions for missing models/packages
- [ ] Check marginal likelihood file existence before loading

**Updating summary table tests (MAINTAINER ONLY):**
- [ ] Set `GENERATE_REFERENCE_FILES <- TRUE`
- [ ] Run tests to generate reference files
- [ ] Review diffs carefully before committing
- [ ] Reset flag to `FALSE`
- **Note**: Contributors/agents should **never** modify `GENERATE_REFERENCE_FILES <- TRUE`
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
