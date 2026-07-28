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

test_that("vignettes use the BayesTools marginal-likelihood contract", {

  vignette_directory <- testthat::test_path("..", "..", "vignettes")
  skip_if_not(
    dir.exists(vignette_directory),
    paste0(
      "Repository vignette sources are not available in this ",
      "installed-package test context."
    )
  )

  vignette_files <- list.files(
    vignette_directory,
    pattern = "\\.Rmd$",
    full.names = TRUE
  )
  violations <- unlist(lapply(vignette_files, function(path){
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    matches <- grep(
      "bridgesampling::bf\\s*\\(",
      lines,
      value = TRUE
    )
    if(length(matches) == 0L){
      return(character())
    }
    paste0(basename(path), ": ", trimws(matches))
  }), use.names = FALSE)

  expect_identical(violations, character())
})

test_that("CI does not bypass vignette construction or checking", {

  workflow_file <- testthat::test_path(
    "..",
    "..",
    ".github",
    "workflows",
    "R-CMD-check.yaml"
  )
  skip_if_not(
    file.exists(workflow_file),
    "Repository workflow sources are not available in this installed-package context."
  )

  workflow <- readLines(
    workflow_file,
    warn = FALSE,
    encoding = "UTF-8"
  )
  bypasses <- workflow[
    grepl("--no-build-vignettes|--ignore-vignettes", workflow)
  ]

  expect_identical(bypasses, character())
  expect_true(any(grepl(
    'pattern = "\\\\.Rmd$"',
    workflow,
    fixed = TRUE
  )))
})

test_that("vignettes never fit models during ordinary rendering", {

  vignette_directory <- testthat::test_path("..", "..", "vignettes")
  skip_if_not(
    dir.exists(vignette_directory),
    paste0(
      "Repository vignette sources are not available in this ",
      "installed-package test context."
    )
  )

  vignette_files <- list.files(
    vignette_directory,
    pattern = "\\.Rmd$",
    full.names = TRUE
  )
  violations <- character()
  for(path in vignette_files){
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    chunk_header <- NULL
    for(line_number in seq_along(lines)){
      line <- lines[[line_number]]
      if(grepl("^```\\{r(?:[ ,}]|$)", line)){
        chunk_header <- line
        next
      }
      if(identical(line, "```")){
        chunk_header <- NULL
        next
      }
      if(is.null(chunk_header)){
        next
      }
      fits_model <- grepl(
        "(JAGS_fit|JAGS_bridgesampling|rstanarm::stan_lmer)\\s*\\(",
        line
      )
      if(fits_model && !grepl("eval\\s*=\\s*FALSE", chunk_header)){
        violations <- c(
          violations,
          paste0(
            basename(path), ":", line_number, ": ",
            trimws(line)
          )
        )
      }
    }
  }

  cache_files <- file.path(
    vignette_directory,
    paste0(tools::file_path_sans_ext(basename(vignette_files)), ".RDS")
  )
  missing_caches <- basename(cache_files[!file.exists(cache_files)])

  expect_identical(violations, character())
  expect_identical(missing_caches, character())
})
