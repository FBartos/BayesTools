# BayesTools Test Organization Guide

## Overview

This document describes the organization and maintenance of tests in the BayesTools package, specifically regarding model fitting and model averaging tests.

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

Models are considered duplicates if they have:
- The same model structure (same JAGS syntax)
- The same type of priors (even with different parameter values)
- The same data structure

For example:
- ✅ GOOD: One model with normal prior `prior("normal", list(0, 1))` and another with spike prior `prior("spike", list(0))`
- ❌ BAD: Two models both with normal priors but different parameters `list(0, 1)` vs `list(0, 2)`

### 3. Model Naming Convention

Models in `test-00-model-fits.R` should follow this naming pattern:
- `fit_{category}_{descriptor}`

Examples:
- `fit_simple_normal` - Simple model with normal priors
- `fit_simple_spike` - Simple model with spike prior
- `fit_formula_simple` - Formula-based model with simple regression
- `fit_formula_treatment` - Formula-based model with treatment factors

### 4. Saving Models with Marginal Likelihoods

When saving a model that has a marginal likelihood:

```r
# Fit the model
fit_model_name <- JAGS_fit(...)

# Compute marginal likelihood (only if model has data)
log_posterior <- function(parameters, data) { ... }
marglik_model_name <- JAGS_bridgesampling(fit_model_name, 
                                          log_posterior = log_posterior,
                                          data = data,
                                          prior_list = priors)

# Save both
result <- save_fit(fit_model_name, "fit_model_name",
                   marglik = marglik_model_name,  # Include marglik here
                   simple_priors = TRUE,
                   note = "Description of the model")
model_registry[["fit_model_name"]] <<- result$registry_entry
fit_model_name <- result$fit
```

This will save:
- `fit_model_name.RDS` - The fitted model
- `fit_model_name_marglik.RDS` - The marginal likelihood

### 5. Loading Models in Other Tests

In `test-model-averaging.R` and other test files:

```r
# Load pre-fitted model
fit_model_name <- readRDS(file.path(temp_fits_dir, "fit_model_name.RDS"))

# Load pre-computed marginal likelihood (if available)
marglik_model_name <- readRDS(file.path(temp_fits_dir, "fit_model_name_marglik.RDS"))

# Use in tests
models <- list(
  list(fit = fit_model_name, marglik = marglik_model_name, prior_weights = 1)
)
```

### 6. When Marginal Likelihoods Are Available

Not all models have marginal likelihoods. Marginal likelihoods are only computed for models that:
- Have actual data (not empty `data = list()`)
- Have a proper log posterior function
- Are **NOT** spike-and-slab or mixture prior models (bridgesampling is not implemented for these)

### 7. Test File Organization

#### test-00-model-fits.R
- **Purpose**: Fit all models and compute marginal likelihoods
- **Run order**: Runs first (prefix `00-`)
- **Dependencies**: None
- **Outputs**: RDS files in temp directory

#### test-model-averaging.R
- **Purpose**: Test model averaging functionality
- **Dependencies**: Requires `test-00-model-fits.R` to run first
- **Inputs**: Loads pre-fitted models and marginal likelihoods
- **Sections**:
  1. Basic function tests (no JAGS required)
  2. JAGS model averaging with pre-fitted models
  3. Print tests with text file comparison
  4. Output generation script (UPDATE_OUTPUT flag)

## Maintenance Checklist

When adding a new test that requires model fitting:

- [ ] Check if a similar model already exists in `test-00-model-fits.R`
- [ ] If not, add the model to `test-00-model-fits.R` (not to the test file itself)
- [ ] Compute and save marginal likelihood if the model has data
- [ ] Update the model registry
- [ ] Add file existence check at the end of the test block
- [ ] In your test file, load the pre-fitted model using `readRDS()`
- [ ] Load the pre-computed marginal likelihood if available
- [ ] Document the model purpose in the `note` parameter

## Example: Adding a New Model for Model Averaging

### Step 1: Add to test-00-model-fits.R

```r
# In the appropriate section of test-00-model-fits.R

# Model: Description
priors_new_model <- list(
  param1 = prior("normal", list(0, 1)),
  param2 = prior("gamma", list(2, 1))
)

model_syntax_new <- "model { ... }"
data_new <- list(y = ..., N = ...)

fit_new_model <- JAGS_fit(model_syntax_new, data_new, priors_new_model,
                          chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)

# Compute marginal likelihood
log_posterior_new <- function(parameters, data) {
  # Define log likelihood
}
marglik_new_model <- JAGS_bridgesampling(fit_new_model,
                                         log_posterior = log_posterior_new,
                                         data = data_new,
                                         prior_list = priors_new_model)

result <- save_fit(fit_new_model, "fit_new_model",
                   marglik = marglik_new_model,
                   simple_priors = TRUE,
                   note = "Description")
model_registry[["fit_new_model"]] <<- result$registry_entry
fit_new_model <- result$fit

# Add to file existence check
expect_true(file.exists(file.path(temp_fits_dir, "fit_new_model.RDS")))
```

### Step 2: Use in test-model-averaging.R

```r
test_that("Test new model averaging", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Load pre-fitted model and marginal likelihood
  fit_new_model <- readRDS(file.path(temp_fits_dir, "fit_new_model.RDS"))
  marglik_new_model <- readRDS(file.path(temp_fits_dir, "fit_new_model_marglik.RDS"))
  
  # Load another model for comparison
  fit_other_model <- readRDS(file.path(temp_fits_dir, "fit_other_model.RDS"))
  marglik_other_model <- readRDS(file.path(temp_fits_dir, "fit_other_model_marglik.RDS"))
  
  # Create model list
  models <- list(
    list(fit = fit_new_model, marglik = marglik_new_model, prior_weights = 1),
    list(fit = fit_other_model, marglik = marglik_other_model, prior_weights = 1)
  )
  
  # Test ensemble inference
  inference <- ensemble_inference(model_list = models, ...)
  
  # Tests
  expect_true(is.list(inference))
  # ... more tests
})
```

## Benefits of This Organization

1. **No duplication**: Each model is fitted exactly once
2. **Faster tests**: Pre-computed marginal likelihoods save time
3. **Easier maintenance**: Model definitions in one place
4. **Cleaner code**: Test files focus on testing, not setup
5. **Consistency**: All models use the same data and parameters across tests
