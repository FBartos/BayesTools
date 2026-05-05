library(testthat)
library(BayesTools)

reporter <- NULL
if (identical(Sys.getenv("AGENT"), "1") &&
    exists("LlmReporter", envir = asNamespace("testthat"), inherits = FALSE)) {
  reporter <- testthat::LlmReporter$new()
}

if (is.null(reporter)) {
  test_check("BayesTools")
} else {
  test_check("BayesTools", reporter = reporter)
}
