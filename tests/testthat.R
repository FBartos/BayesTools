library(testthat)
library(BayesTools)

reporter <- NULL
if (identical(Sys.getenv("AGENT"), "1") &&
    exists("LlmReporter", envir = asNamespace("testthat"), inherits = FALSE)) {
  reporter <- testthat::LlmReporter$new()
}

# CRAN runs only the smoke test; the profiles cover everything else.
on_cran <- get("on_cran", envir = asNamespace("testthat"), inherits = FALSE)
filter  <- NULL
if (on_cran()) {
  Sys.setenv(BAYESTOOLS_TEST_PROFILE = "cran")
  filter <- "cran-smoke"
}

if (is.null(reporter)) {
  test_check("BayesTools", filter = filter)
} else {
  test_check("BayesTools", filter = filter, reporter = reporter)
}
