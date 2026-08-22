if (interactive()) {

  library(devtools)
  library(testthat)
  library(vdiffr)

  .bayestools_profile_root <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
  source(
    file.path(.bayestools_profile_root, ".dev", "test-tests.R"),
    local = .GlobalEnv
  )
  rm(.bayestools_profile_root, envir = .GlobalEnv)
}
