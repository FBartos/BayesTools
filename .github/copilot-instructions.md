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
- **Test**: `devtools::test()` (Timeout: 5m+)
- **Targeted Test**: `devtools::test(filter = 'priors')` (Recommended for dev)
- **Check**: `rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'), error_on = 'never')` 
- **Docs**: `devtools::document()` 

### Testing Strategy (`tests/testthat/`)
- **Priors**: Use the `test_prior(prior)` helper in `tests/testthat/test-priors.R` to validate new priors (checks rng, pdf, cdf, quant consistency).
- **Visual Regression**: Use `vdiffr::expect_doppelganger` for plot tests.
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
2.  Verify with `devtools::test(filter = 'JAGS')`.
3.  Ensure backward compatibility with `runjags` objects. 

## Feature Addition Steps
1. Create unit tests verifying the desired behavior.
2. Implement the feature.
3. Verify the tests pass.
4. Update documentation as needed.
5. Update NEWS.md with a summary of changes.