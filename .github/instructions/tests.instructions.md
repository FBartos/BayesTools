````instructions
---
applyTo: "**/tests/testthat/*.R"
---

# BayesTools Test Guidelines

## Overview

- Model fitting is centralized in `test-00-model-fits.R`; other tests load cached models
- **testthat Edition 3** - do not use `context()` calls
- Tests use `common-functions.R` for shared helpers

## Test Caching (TDD Workflow)

Model fitting is slow. The caching system lets you run the full suite once and reuse fits.

### Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `BAYESTOOLS_TEST_FILES_DIR` | Cache directory location | `../temp/BayesTools_test_files` |
| `BAYESTOOLS_TEST_SKIP_REFIT` | Skip fitting if cache exists | TRUE |

### Recommended TDD Workflow

```r
# 1. Run full suite once to verify current code and populate cache (if missing) 
devtools::test()

# 2. Iterate on your feature (uses cached fits)
devtools::test(filter = "your-feature")

# 3. Final verification (disable cache if fit / marglik code or its dependencies changed)
clean_cached_fits()
devtools::test()
```

### When to Clear Cache

Clear with `clean_cached_fits()` when you modify:
- `JAGS_fit()` or `JAGS_bridgesampling()` logic (or any of its dependencies)
- Model definitions in `test-00-model-fits.R`

## Key Rules

### Model Fitting

- **Only `test-00-model-fits.R`** fits models and computes marginal likelihoods
- Other tests load with `readRDS(file.path(temp_fits_dir, "model_name.RDS"))`
- Check `model_registry.RDS` for available models before creating new ones

### File Naming

| Pattern | Purpose |
|---------|---------|
| `test-{feature}.R` | Main tests |
| `test-{feature}-input.R` | Input validation |
| `test-{feature}-coverage.R` | Edge cases |

### Skip Conditions

| Condition | When to Use |
|-----------|-------------|
| `skip_if_no_fits()` | Test loads pre-fitted models |
| `skip_if_not_installed("rjags")` | Test requires JAGS |
| `skip_if_not_installed("bridgesampling")` | Test computes marginal likelihoods |
| `skip_if_not_installed("vdiffr")` | Visual regression tests |

### Helper Functions (common-functions.R)

```r
source(testthat::test_path("common-functions.R"))

# Prior testing
test_prior(prior)
test_weightfunction(prior)
test_orthonormal(prior)

# Reference file testing
test_reference_table(table, filename)
test_reference_text(text, filename)

# Skip/cache helpers
skip_if_no_fits()
skip_refit_if_cached(name)
clean_cached_fits()
```

## AI Agent Protocol

1. **Scan `test-00-model-fits.R` first** - understand available models
2. **Reuse existing models** - don't create duplicates
3. **Never fit models** outside `test-00-model-fits.R`
4. **Never modify** `GENERATE_REFERENCE_FILES` flag (maintainer only)

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Pre-fitted models not available" | Run `devtools::test(filter = "00-model-fits")` |
| Stale cache causing failures | `clean_cached_fits()` then rerun |
| Tests pass locally, fail on CI | Clear cache, run full suite |

````
