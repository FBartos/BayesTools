skip_if_not_test_profile("unit")

test_that("live fitting is confined to the centralized fit file", {

  test_directory <- testthat::test_path()
  test_files <- list.files(
    test_directory,
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
  centralized_file <- "test-00-model-fits.R"

  declares_fit_profile <- vapply(test_files, function(path){
    lines <- readLines(path, warn = FALSE, n = 20L)
    any(grepl(
      "^\\s*skip_if_not_test_profile\\([\"']fit[\"']\\)",
      lines
    ))
  }, logical(1))
  expect_setequal(
    basename(test_files[declares_fit_profile]),
    centralized_file
  )

  backend_patterns <- c(
    "rjags::jags\\.model\\s*\\(",
    "rjags::coda\\.samples\\s*\\(",
    "runjags::run\\.jags\\s*\\("
  )
  noncentralized <- test_files[basename(test_files) != centralized_file]
  violations <- unlist(lapply(noncentralized, function(path){
    lines <- readLines(path, warn = FALSE)
    matched <- unique(unlist(lapply(backend_patterns, function(pattern){
      grep(pattern, lines, value = TRUE)
    })))
    if(length(matched) == 0L){
      return(character())
    }
    paste0(basename(path), ": ", trimws(matched))
  }), use.names = FALSE)

  expect_identical(violations, character())
})
