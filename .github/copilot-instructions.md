# BayesTools

BayesTools is an R package that provides tools for conducting Bayesian analyses and Bayesian model averaging. The package contains functions for creating prior distribution objects, automating JAGS model syntax generation, mixing posterior samples, and generating summary tables for Bayesian inference.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

Bootstrap, build, and test the repository:
- Install system dependencies: `sudo apt update && sudo apt install -y jags r-cran-devtools r-cran-testthat r-cran-rjags`
- Install the package: `R -e "devtools::install()"` -- takes 2 minutes. NEVER CANCEL. Set timeout to 5+ minutes.
- Run tests: `R -e "devtools::test()"` -- takes 15 minutes. NEVER CANCEL. Set timeout to 30+ minutes.
- Run R CMD check: `R -e "rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran', '--ignore-vignettes'), build_args = c('--no-build-vignettes'), error_on = 'never', check_dir = 'check')"` -- takes 1 minute. Set timeout to 5+ minutes.
- Update documentation: `R -e "devtools::document()"` -- takes 30 seconds. Set timeout to 2+ minutes.

## Validation

- **CRITICAL VALIDATION**: Always test basic package functionality after making changes by running: `R -e "library(BayesTools); p <- prior('normal', list(0, 1)); print(p); plot(p)"`
- **NEVER CANCEL** any R CMD operations - builds and tests may take 15+ minutes but must complete
- Always run `devtools::document()` after modifying roxygen2 comments or function signatures
- Always run `rcmdcheck::rcmdcheck()` before finalizing changes to ensure CRAN compliance
- Test the specific functionality you modified - BayesTools has three main areas:
  1. **Prior distributions**: Test with `prior()` function and plotting
  2. **JAGS integration**: Test with `JAGS_fit()` for model fitting (requires longer running time)
  3. **Model averaging**: Test with posterior mixing functions

## CRAN Compliance and Package Standards

- **STRICT REQUIREMENT**: All code must follow CRAN guidelines
- Use snake_case for function and variable names throughout the codebase
- All exported functions require complete roxygen2 documentation with @param, @return, @examples
- Place R scripts in R/ directory, tests in tests/testthat/, vignettes in vignettes/
- Use message(), warning(), and stop() for user feedback instead of hardcoded messages
- Maintain API compatibility with downstream packages that re-export BayesTools functions
- Type-check all user-facing output since functions are often re-exported

## System Requirements

- **JAGS >= 4.3.0** is required and must be installed system-wide before running tests
- R version 4.0+ with key dependencies: rjags, runjags, bridgesampling, coda, mvtnorm, extraDistr
- For development: devtools, testthat, vdiffr, rcmdcheck
- **Note**: Some tests are skipped on different OS platforms due to numerical precision differences

## Common Tasks

### Build and Install Package
```r
# From repository root
devtools::install()
```
**Timing**: 2 minutes. **NEVER CANCEL** - set timeout to 5+ minutes.

### Run Complete Test Suite
```r
devtools::test()
```
**Timing**: 15 minutes for full suite. **NEVER CANCEL** - set timeout to 30+ minutes.
Many tests are skipped on CRAN to reduce build time. Visual regression tests (vdiffr) may fail in different environments - this is expected.

### R CMD Check (CRAN Validation)
```r
rcmdcheck::rcmdcheck(
  args = c("--no-manual", "--as-cran", "--ignore-vignettes"),
  build_args = c("--no-build-vignettes"),
  error_on = "never",
  check_dir = "check"
)
```
**Timing**: 1 minute. **NEVER CANCEL** - set timeout to 5+ minutes.

### Update Documentation
```r
devtools::document()
```
**Timing**: 30 seconds. Always run after modifying roxygen2 comments.

## Key Repository Structure

### Repository root:
```
.
├── DESCRIPTION          # Package metadata and dependencies
├── NAMESPACE           # Exported functions (auto-generated)
├── README.Rmd          # Source for README.md (edit this, not README.md)
├── README.md           # Generated from README.md
├── BayesTools.Rproj    # RStudio project file
├── R/                  # All R source code
├── man/                # Documentation (auto-generated from roxygen2)
├── tests/              # Test suite
├── vignettes/          # Package vignettes
├── data/               # Package datasets
├── inst/               # Additional package files
└── .github/            # GitHub Actions workflows
```

### Key R Files in R/:
- `priors.R` - Core prior distribution functions
- `JAGS-fit.R` - JAGS model fitting automation
- `JAGS-formula.R` - Formula interface for JAGS
- `JAGS-marglik.R` - Marginal likelihood computation
- `model-averaging.R` - Bayesian model averaging functions
- `summary-tables.R` - Result summary and formatting
- `tools.R` - Input validation and utility functions

### Tests Structure in tests/testthat/:
- `test-*.R` files contain comprehensive test cases
- Uses testthat framework with vdiffr for visual regression testing
- Tests integrate with JAGS for realistic Bayesian analysis scenarios
- Many tests are conditional (`skip_on_cran()`, `skip_if_not_installed()`)

## GitHub Actions CI/CD

The repository uses several workflows:
- **R-CMD-check.yaml**: Runs across Windows, Ubuntu, and macOS with full CRAN checks
- **R-test.yaml**: Focused unit testing on Windows
- **test-coverage.yaml**: Code coverage reporting
- **pkgdown.yaml**: Documentation site generation

All workflows handle JAGS installation automatically for their respective platforms.

## Development Workflow

1. **Make changes**: Edit R files in R/ directory
2. **Update docs**: Run `devtools::document()` if you modified roxygen2 comments
3. **Install locally**: Run `devtools::install()` to test changes
4. **Validate basic functionality**: Test core features work as expected
5. **Run tests**: Run `devtools::test()` for comprehensive validation
6. **CRAN check**: Run `rcmdcheck::rcmdcheck()` to ensure CRAN compliance
7. **Commit**: Ensure .gitignore excludes build artifacts (check/, Rplots.pdf, etc.)

## Performance Considerations

- Functions may be called repeatedly in downstream packages - optimize accordingly
- JAGS model fitting can be computationally intensive - tests reflect realistic usage
- Use efficient data structures for prior objects and posterior samples
- Consider memory usage for large posterior sample arrays

## Troubleshooting

- **JAGS not found**: Ensure JAGS >= 4.3.0 is installed system-wide
- **Test failures**: Visual tests (vdiffr) may fail on different systems - this is often acceptable
- **Package won't install**: Check that all dependencies are available, especially rjags
- **Long test times**: This is normal - Bayesian computation is inherently time-intensive
- **Memory issues**: Some tests work with large posterior samples - ensure adequate RAM

## Example Validation Commands

After making changes, validate with these commands:

```r
# Test basic functionality
library(BayesTools)
p <- prior("normal", list(0, 1))
print(p)
plot(p)

# Test JAGS integration (if JAGS-related changes)
data_test <- list(x = rnorm(10), N = 10)
priors_test <- list(mu = prior("normal", list(0, 1)))
model_syntax <- "model{ for(i in 1:N){ x[i] ~ dnorm(mu, 1) } }"
# Note: Full JAGS test requires more setup and time
```

Always ensure changes maintain backward compatibility and follow the existing code patterns in the package.
