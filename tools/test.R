setwd("C:/R-Packages/BayesTools")
clean_cached_fits()
test_files_dir <- file.path(tempdir(), "BayesTools_test_files")
dir.create(test_files_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  BAYESTOOLS_TEST_PROFILE = "all",
  BAYESTOOLS_TEST_SKIP_REFIT = "false",
  BAYESTOOLS_TEST_FILES_DIR = normalizePath(test_files_dir, winslash = "/", mustWork = TRUE),
  NOT_CRAN = "true",
  VDIFFR_RUN_TESTS = "true"
)

Sys.unsetenv("AGENT")

devtools::test()
