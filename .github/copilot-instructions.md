This is an R package that must strictly follow CRAN guidelines for code, documentation, and package structure.

Always use snake_case for function and variable names throughout the codebase.

Place all R scripts in the R/ directory, tests in tests/testthat/, and vignettes in vignettes/.

All exported functions require roxygen2 documentation with clear parameter descriptions, return values, and examples.

Use testthat for all tests and ensure tests cover edge cases, error handling, and input validation.

When deprecating functions or arguments, use lifecycle badges and provide deprecation warnings with clear migration paths.

Optimize for performance, especially for functions that may be called repeatedly in downstream packages.

Use message(), warning(), and stop() for user feedback instead of hardcoded messages.

Ensure all examples and vignettes are reproducible and do not rely on random seeds unless explicitly set.

Avoid breaking changes to the API; maintain compatibility with downstream packages that re-export functions.

All user-facing output from functions like prior() must be type-checked since these are often re-exported in downstream packages.

BayesTools automates interaction with JAGS and posterior samples from Stan for Bayesian analysis workflows.

The package provides interfaces for creating prior distributions, transforming R formulas into JAGS code, and generating posterior summaries.

Avoid using non-CRAN dependencies to maintain CRAN compliance.

Follow test-driven development practices and ensure comprehensive test coverage for both user-facing and developer-facing functions.
