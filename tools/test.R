cache_dir <- file.path(
  "C:/R-Packages/temp",
  paste0("BayesTools_test_files-", format(Sys.time(), "%Y%m%d-%H%M%S"))
)
dir.create(cache_dir, recursive = TRUE)

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
