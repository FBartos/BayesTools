# BayesTools Copilot Instructions

BayesTools is an R package for Bayesian analyses, JAGS model automation, and Bayesian model averaging.

## Architecture & Key Components

- **Priors (`R/priors.R`)**: The core `prior` S3 class.
  - Constructed via `prior(distribution, parameters, ...)`.
  - Must implement S3 methods: `print`, `plot`, `rng` (random generation), `quant` (quantile), `cdf`, `pdf`.
  - See `R/priors.R` for the factory function and `R/distributions-*.R` for specific implementations.
- **JAGS Integration (`R/JAGS-fit.R`, `R/JAGS-formula.R`)**:
  - `JAGS_fit` wraps `runjags::run.jags`.
  - `JAGS_formula` handles formula parsing, data preparation, and prior assignment for formula-based models.
  - Automates model syntax generation (`JAGS_add_priors`), data preparation, and initialization.
  - Supports `autofit` for convergence.
- **Model Averaging (`R/model-averaging.R`)**:
  - `compute_inference` and `ensemble_inference` calculate posterior probabilities and Bayes factors.
  - Relies on `bridgesampling` for marginal likelihoods.
- **Input Validation (`R/tools.R`)**:
  - Centralized validation functions used throughout the package.

## Coding Conventions

- **Input Validation**: **STRICTLY** use the project's internal check functions from `R/tools.R` for all user inputs.
  - `check_bool(x, name)` 
  - `check_char(x, name)` 
  - `check_int(x, name)` 
  - `check_real(x, name)` 
  - `check_list(x, name)` 
  - Do not use `stopifnot` or manual `if` checks for standard type validation.
- **Naming**: Use `snake_case` for all function and variable names.
- **Documentation**: Full `roxygen2` documentation for all exported functions.
- **Error Handling**: Use `stop()`, `warning()`, `message()` with clear, descriptive messages.

## Development Workflow

### Build & Test
- **Install**: `devtools::install()` (Timeout: 5m+)
- **Package-critical tests**: `Rscript tools/test-profile.R unit`
- **Cached-fit tests**: `Rscript tools/test-profile.R fixture`
- **Pure visual tests**: `Rscript tools/test-profile.R visual`
- **Cached-JAGS visual tests**: `Rscript tools/test-profile.R visual-fixture`
- **Fit/marglik tests**: `Rscript tools/test-profile.R fit` refreshes cached fits by default; set `BAYESTOOLS_TEST_SKIP_REFIT=TRUE` only to reuse a validated cache.
- **Full developer verification**: `Rscript tools/test-profile.R all`
- **Targeted Test**: `devtools::test(filter = 'priors')` only for ad hoc local iteration after choosing the relevant profile
- **Check**: `rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'), error_on = 'never')` 
- **Docs**: `devtools::document()` 

### Testing Strategy (`tests/testthat/`)
- **Priors**: Use the `test_prior(prior)` helper in `tests/testthat/test-priors.R` to validate new priors (checks rng, pdf, cdf, quant consistency).
- **Visual Regression**: Use `vdiffr::expect_doppelganger` for plot tests.
- **Profiles**: Add `skip_if_not_test_profile()` at the top of every `test-*.R` file. Use `visual-fixture` for plot tests that load cached JAGS fits.
- **Performance**: Use `skip_on_cran()` for computationally intensive tests (JAGS fitting).
- **JAGS**: Ensure `JAGS >= 4.3.0` is installed.

## Common Tasks

### Adding a New Prior
1.  Add distribution name to `prior` function in `R/priors.R`.
2.  Implement `rng`, `pdf`, `cdf`, `quant` methods in `R/distributions-*.R`.
3.  Add `plot.prior` support if needed.
4.  Add unit test in `tests/testthat/test-priors.R` using `test_prior()`.

### Modifying JAGS Fitting
1.  Edit `R/JAGS-fit.R` for general fitting logic or `R/JAGS-formula.R` for formula handling.
2.  Verify package-critical code with `Rscript tools/test-profile.R unit`.
3.  Rebuild and verify fit fixtures with `Rscript tools/test-profile.R fit`.
4.  Verify cached-fit consumers with `Rscript tools/test-profile.R fixture` and `Rscript tools/test-profile.R visual-fixture`.

## Feature Addition Steps
1. Create unit tests verifying the desired behavior.
2. Implement the feature.
3. Verify the tests pass.
4. Update documentation as needed.
5. Update NEWS.md with a summary of changes.
