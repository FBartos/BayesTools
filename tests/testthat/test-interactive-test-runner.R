skip_if_not_test_profile("unit")


test_that("quiet LLM reporter counts skips without printing each reason", {

  output_path <- tempfile("bayestools-quiet-reporter-", fileext = ".txt")
  on.exit(unlink(output_path), add = TRUE)
  reporter <- bayestools_quiet_llm_reporter(file = output_path)
  skip <- structure(
    list(message = "routine profile skip"),
    class = c("expectation_skip", "expectation", "condition")
  )
  success <- structure(
    list(),
    class = c("expectation_success", "expectation", "condition")
  )
  failure <- tryCatch(
    testthat::expectation("failure", "deliberate failure"),
    error = identity
  )

  reporter$add_result("profiles", "quiet output", skip)
  reporter$add_result("profiles", "quiet output", success)
  reporter$add_result("profiles", "quiet output", failure)
  reporter$end_reporter()
  output <- readLines(output_path, warn = FALSE)

  expect_equal(reporter$n_skip, 1L)
  expect_equal(reporter$n_ok, 1L)
  expect_equal(reporter$n_fail, 1L)
  expect_true(any(grepl("deliberate failure", output, fixed = TRUE)))
  expect_true(any(grepl("SKIP 1", output, fixed = TRUE)))
  expect_false(any(grepl("routine profile skip", output, fixed = TRUE)))
})


test_that("interactive runner dispatches comprehensive and filtered profiles", {

  runner_env <- new.env(parent = globalenv())
  source(
    testthat::test_path("..", "..", ".dev", "test-tests.R"),
    local = runner_env
  )
  calls   <- list()
  cleaned <- 0L
  assign(
    ".run_bayestools_test_profile",
    function(profile, filter = NULL) {

      calls[[length(calls) + 1L]] <<- list(
        profile  = profile,
        filter   = filter,
        reporter = Sys.getenv("BAYESTOOLS_TEST_REPORTER"),
        agent     = Sys.getenv("AGENT", unset = NA_character_),
        quiet     = Sys.getenv("BAYESTOOLS_TEST_QUIET_SKIPS")
      )
      invisible(TRUE)
    },
    envir = runner_env
  )
  assign(
    ".bayestools_clean_test_cache",
    function(...) {
      cleaned <<- cleaned + 1L
      invisible(TRUE)
    },
    envir = runner_env
  )

  expected_formals <- c(
    "filter", "reporter", "refit", "update", "update_timings",
    "regenerate", "load_package", "stop_on_failure", "root"
  )
  expect_identical(names(formals(runner_env$test_tests)), expected_formals)
  expect_identical(formals(runner_env$test_tests)[["reporter"]], "progress")

  runner_env$test_tests(
    filter          = "JAGS-parameter-catalog",
    refit           = TRUE,
    load_package    = FALSE,
    stop_on_failure = TRUE,
    root            = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    c("fit", "filter")
  )
  expect_identical(calls[[2L]][["filter"]], "JAGS-parameter-catalog")
  expect_identical(calls[[2L]][["reporter"]], "progress")
  expect_true(is.na(calls[[2L]][["agent"]]))
  expect_identical(calls[[2L]][["quiet"]], "FALSE")
  expect_identical(cleaned, 1L)

  calls <- list()
  runner_env$test_tests(
    filter       = "interpret",
    load_package = FALSE,
    root         = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    "filter"
  )
  expect_identical(cleaned, 1L)

  calls <- list()
  runner_env$test_tests(
    filter       = "interpret",
    reporter     = "llm",
    load_package = FALSE,
    root         = testthat::test_path()
  )
  expect_identical(calls[[1L]][["reporter"]], "llm")
  expect_identical(calls[[1L]][["agent"]], "1")
  expect_identical(calls[[1L]][["quiet"]], "TRUE")

  calls <- list()
  runner_env$test_tests(
    load_package = FALSE,
    root         = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    "all"
  )
  expect_identical(cleaned, 1L)

  calls <- list()
  runner_env$test_tests(
    regenerate   = TRUE,
    load_package = FALSE,
    root         = testthat::test_path()
  )
  expect_identical(
    vapply(calls, `[[`, character(1), "profile"),
    "all"
  )
  expect_identical(cleaned, 2L)

  expect_error(
    runner_env$test_tests(
      update_timings = TRUE,
      load_package   = FALSE,
      root           = testthat::test_path()
    ),
    "tests have no timing baselines"
  )
})
