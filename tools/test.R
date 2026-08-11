cache_dir <- file.path(
  "C:/R-Packages/temp",
  paste0("BayesTools_test_files-", format(Sys.time(), "%Y%m%d-%H%M%S"))
)
dir.create(cache_dir, recursive = TRUE)

Sys.setenv(
  AGENT                       = "1",
  BAYESTOOLS_TEST_PROFILE     = "all",
  BAYESTOOLS_TEST_FILES_DIR   = normalizePath(cache_dir, winslash = "/"),
  BAYESTOOLS_TEST_SKIP_REFIT  = "false",
  NOT_CRAN                    = "true",
  VDIFFR_RUN_TESTS            = "true"
)

devtools::test(
  stop_on_failure = TRUE,
  reporter = testthat::LlmReporter$new()
)
