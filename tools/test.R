args <- commandArgs(trailingOnly = TRUE)
profile <- if(length(args) > 0L && nzchar(args[[1L]])) args[[1L]] else "all"

Sys.setenv(
  BAYESTOOLS_TEST_PROFILE = profile,
  BAYESTOOLS_TEST_SKIP_REFIT = "false"
)

source(file.path("tools", "test-profile.R"))
