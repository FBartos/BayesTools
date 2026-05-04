---
description: 'R unit tests: testing standards and Copilot guidance for interaction with unit tests.'
applyTo: "**/tests/testthat/*.R"
---

# BayesTools Test Guidelines

## Overview

- Model fitting is centralized in `test-00-model-fits.R`; other tests load cached models
- **testthat Edition 3** - do not use `context()` calls
- Tests use `common-functions.R` for shared helpers
- Tests are split into package-critical, fixture/cache, pure visual, cached-JAGS visual, and model-fit profiles

## Test Profiles

Use `tools/test-profile.R` instead of calling `devtools::test()` directly for profiled test runs. The runner sets the right environment variables and passes a testthat file filter, which prevents skipped visual files from pruning local vdiffr snapshots.

| Profile | Command | Purpose |
|---------|---------|---------|
| `unit` | `Rscript tools/test-profile.R unit` | Fast package-critical tests; no real JAGS fits or visual regression expectations |
| `fixture` | `Rscript tools/test-profile.R fixture` | Tests that load cached model fits, tables, and reference files |
| `visual` | `Rscript tools/test-profile.R visual` | Pure vdiffr plot tests; sets `VDIFFR_RUN_TESTS=true` and does not require cached JAGS fits |
| `visual-fixture` | `Rscript tools/test-profile.R visual-fixture` | vdiffr plot tests that load cached JAGS fits; sets `VDIFFR_RUN_TESTS=true` and `NOT_CRAN=true` |
| `fit` | `Rscript tools/test-profile.R fit` | Slow model fitting and marginal-likelihood tests; sets `NOT_CRAN=true` and refreshes cached fits by default |
| `all` | `Rscript tools/test-profile.R all` | Full developer verification across every profile |

Supported aliases include `fixtures`, `plots`, `snapshot`, `snapshots`, `jags-visual`, `jags_visual`, `cached-visual`, `cached_visual`, `visual_fixture`, `visual-fixtures`, `visual_fixtures`, `heavy`, and `slow`.

## Test Caching

Model fitting is slow. The caching system lets you run the full suite once and reuse fits.

### Environment Variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `BAYESTOOLS_TEST_FILES_DIR` | Cache directory location | `file.path(tempdir(), "BayesTools_test_files")` |
| `BAYESTOOLS_TEST_SKIP_REFIT` | Skip fitting if a validated cache exists | TRUE in helpers; `tools/test-profile.R fit` sets FALSE when unset |
| `BAYESTOOLS_TEST_PROFILE` | Active test profile when not using `tools/test-profile.R` | `all` in helpers, `unit` in runner |
| `VDIFFR_RUN_TESTS` | Enables vdiffr visual regression checks | unset |
| `NOT_CRAN` | Enables slow/heavy checks | unset |

### Recommended TDD Workflow

```r
# 1. Run package-critical tests during normal development
system("Rscript tools/test-profile.R unit")

# 2. If your change touches cached-fit consumers, run fixture tests
system("Rscript tools/test-profile.R fixture")

# 3. If your change touches plots or snapshots, run visual tests
system("Rscript tools/test-profile.R visual")

# 4. If your change touches cached JAGS plot output, run cached visual tests
system("Rscript tools/test-profile.R visual-fixture")

# 5. Final verification when fit / marglik code or its dependencies changed.
#    The fit profile refreshes cached fits by default.
system("Rscript tools/test-profile.R fit")

# 6. Full local verification before a large release candidate
system("Rscript tools/test-profile.R all")
```

### When to Clear Cache

`Rscript tools/test-profile.R fit` refreshes cached fits by default. Clear with `clean_cached_fits()` when you want to remove old artifacts before fitting, or set `BAYESTOOLS_TEST_SKIP_REFIT=TRUE` only when you intentionally want to reuse a previously validated fit cache.

Clear cache when you modify:
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
| `skip_if_not_test_profile("unit")` | Fast package-critical test files |
| `skip_if_not_test_profile("fixture")` | Tests that require cached fit/reference fixtures |
| `skip_if_not_test_profile("visual")` | Visual-only test files |
| `skip_if_not_test_profile("visual-fixture")` | Cached-JAGS visual test files |
| `skip_if_not_test_profile("fit")` | Slow fitting/marginal-likelihood tests |
| `skip_if_no_fits()` | Test loads pre-fitted models |
| `skip_if_no_fit_cache()` | Alias for tests that load cached fits |
| `skip_if_not_visual_tests()` | vdiffr expectations inside mixed unit/visual files |
| `skip_if_not_visual_fixture_tests()` | vdiffr expectations that require cached JAGS fits |
| `skip_if_not_heavy_tests()` | Heavy tests inside mixed files |
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
skip_if_not_test_profile("unit")
skip_if_not_visual_tests()
skip_if_not_visual_fixture_tests()
skip_if_not_heavy_tests()
skip_if_no_fits()
skip_if_no_fit_cache()
skip_refit_if_cached(name)
mark_refit_cache_complete(name)
clean_cached_fits()
```

## AI Agent Protocol

1. **Scan `test-00-model-fits.R` first** - understand available models
2. **Reuse existing models** - don't create duplicates
3. **Never fit models** outside `test-00-model-fits.R`
4. **Never modify** `GENERATE_REFERENCE_FILES` flag (maintainer only)
5. **Use `tools/test-profile.R`** for profile runs; do not rely on plain `devtools::test()` for non-all profiles

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Pre-fitted models not available" | Run `Rscript tools/test-profile.R fit` |
| Stale cache causing failures | Run `Rscript tools/test-profile.R fit`; use `clean_cached_fits()` first if you need a fully empty cache |
| Tests pass locally, fail on CI | Clear cache, run `Rscript tools/test-profile.R all` |
| vdiffr snapshots disappear during non-visual runs | Use `Rscript tools/test-profile.R unit` or `fixture`; the runner applies the required file filter |
