clean_cached_fits()
Sys.unsetenv("AGENT")
Sys.setenv(
  BAYESTOOLS_TEST_PROFILE    = "all",
  BAYESTOOLS_TEST_SKIP_REFIT = "true",
  NOT_CRAN                   = "true",
  VDIFFR_RUN_TESTS           = "true"
)

devtools::test()

source(file.path("tests", "testthat", "helper-00-reference-table-review.R"))
review_test_snapshots()
