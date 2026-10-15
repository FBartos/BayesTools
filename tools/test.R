# The fit cache is the configured BAYESTOOLS_TEST_FILES_DIR, or a directory of
# the session's temporary directory (removed with the session).
cache_dir <- Sys.getenv("BAYESTOOLS_TEST_FILES_DIR")
if (!nzchar(cache_dir)) {
  cache_dir <- file.path(tempdir(), "BayesTools_test_files")
}
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# The test helpers read the fit cache location from BAYESTOOLS_TEST_FILES_DIR;
# without it they fall back to tempdir() and this directory stays unused.
Sys.setenv(
  AGENT                       = "1",
  BAYESTOOLS_TEST_FILES_DIR   = normalizePath(cache_dir, winslash = "/", mustWork = TRUE),
  BAYESTOOLS_TEST_PROFILE     = "all",
  BAYESTOOLS_TEST_SKIP_REFIT  = "false",
  NOT_CRAN                    = "true",
  VDIFFR_RUN_TESTS            = "true"
)

devtools::test(
  stop_on_failure = TRUE,
  reporter = testthat::LlmReporter$new()
)
