skip_if_not_test_profile("unit")

test_that("JAGS fit settings reject missing and non-finite controls", {
  valid <- list(
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    thin = 1,
    autofit = FALSE,
    parallel = FALSE,
    cores = 1,
    silent = TRUE,
    seed = 1
  )
  invalid <- list(
    chains = NA_real_,
    adapt = Inf,
    burnin = NA_real_,
    sample = Inf,
    thin = NA_real_,
    autofit = NA,
    parallel = NA,
    cores = Inf,
    silent = NA,
    seed = NaN
  )

  expect_silent(do.call(JAGS_check_and_list_fit_settings, valid))
  for(name in names(invalid)){
    settings <- valid
    settings[[name]] <- invalid[[name]]
    expect_error(
      do.call(JAGS_check_and_list_fit_settings, settings),
      paste0("'", name, "'")
    )
  }
})

test_that("JAGS autofit settings reject missing and non-finite controls", {
  valid <- list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10,
    check_indicators = FALSE,
    monitor = c("mu", "sigma"),
    allow_not_assessable = FALSE
  )
  invalid <- list(
    max_Rhat = Inf,
    min_ESS = Inf,
    max_error = NaN,
    max_SD_error = NA_real_,
    max_time = list(time = Inf, unit = "mins"),
    sample_extend = Inf,
    restarts = NA_real_,
    max_extend = Inf,
    check_indicators = NA,
    monitor = NA_character_,
    allow_not_assessable = NA
  )

  expect_silent(JAGS_check_and_list_autofit_settings(valid))
  for(name in names(invalid)){
    settings <- valid
    settings[[name]] <- invalid[[name]]
    expect_error(
      JAGS_check_and_list_autofit_settings(settings),
      paste0("'", if(name == "max_time") "max_time:time" else name, "'")
    )
  }

  invalid_unit <- valid
  invalid_unit$max_time$unit <- NA_character_
  expect_error(
    JAGS_check_and_list_autofit_settings(invalid_unit),
    "'max_time:unit'"
  )

  empty_monitor <- valid
  empty_monitor$monitor <- character()
  expect_error(
    JAGS_check_and_list_autofit_settings(empty_monitor),
    "'monitor' argument must select at least one parameter"
  )
})

test_that("JAGS_extend validates runtime controls before extension", {
  fit <- structure(list(), class = "BayesTools_fit")
  invalid <- list(
    parallel = NA,
    cores = Inf,
    silent = NA
  )

  for(name in names(invalid)){
    arguments <- list(fit = fit)
    arguments[[name]] <- invalid[[name]]
    expect_error(
      do.call(JAGS_extend, arguments),
      paste0("'", name, "'")
    )
  }
})

.jags_extend_test_fit <- function(){

  fit <- structure(
    list(),
    class = c("runjags", "BayesTools_fit")
  )
  attr(fit, "prior_list") <- list()
  attr(fit, "model_syntax") <- "model{}"
  attr(fit, "required_packages") <- character()
  attr(fit, "jags_modules") <- character()
  attr(fit, "add_parameters") <- character()
  attr(fit, "parameter_registry") <- build_test_parameter_registry(character())
  fit
}

.jags_extend_test_control <- function(max_time = list(time = 60, unit = "secs")){
  list(
    max_Rhat = NULL,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    max_time = max_time,
    sample_extend = 1,
    restarts = 1,
    max_extend = 1,
    check_indicators = FALSE
  )
}

test_that("JAGS_extend preserves the last valid fit after a backend error", {

  skip_if_not_installed("runjags")
  fit <- .jags_extend_test_fit()
  testthat::local_mocked_bindings(
    extend.jags = function(...){
      stop("backend exploded")
    },
    .package = "runjags"
  )

  expect_false("seed" %in% names(formals(JAGS_extend)))
  expect_warning(
    result <- JAGS_extend(
      fit,
      autofit_control = .jags_extend_test_control()
    ),
    "returning the last valid fit.*backend exploded"
  )
  expect_identical(result, fit)
})

test_that("JAGS_extend resets its time budget for every call", {

  skip_if_not_installed("runjags")
  fit <- .jags_extend_test_fit()
  extension_calls <- 0L
  clock_calls <- 0L
  clock_values <- as.POSIXct(
    c(
      "2026-01-01 00:00:00",
      "2026-01-01 00:00:01",
      "2026-01-01 01:00:00",
      "2026-01-01 01:00:01"
    ),
    tz = "UTC"
  )
  testthat::local_mocked_bindings(
    extend.jags = function(runjags.object, ...){
      extension_calls <<- extension_calls + 1L
      runjags.object
    },
    .package = "runjags"
  )
  testthat::local_mocked_bindings(
    .bt_jags_extend_time = function(){
      clock_calls <<- clock_calls + 1L
      clock_values[[clock_calls]]
    },
    JAGS_check_convergence = function(...) TRUE,
    .package = "BayesTools"
  )

  control <- .jags_extend_test_control(
    max_time = list(time = 2, unit = "secs")
  )
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_equal(extension_calls, 2L)
  expect_equal(clock_calls, 4L)
})

test_that("JAGS_extend forwards explicit convergence monitor policy", {

  skip_if_not_installed("runjags")
  fit <- .jags_extend_test_fit()
  convergence_arguments <- NULL
  testthat::local_mocked_bindings(
    extend.jags = function(runjags.object, ...){
      runjags.object
    },
    .package = "runjags"
  )
  testthat::local_mocked_bindings(
    JAGS_check_convergence = function(...){
      convergence_arguments <<- list(...)
      TRUE
    },
    .package = "BayesTools"
  )

  control <- .jags_extend_test_control()
  control$monitor <- "mu"
  control$allow_not_assessable <- TRUE
  expect_silent(JAGS_extend(fit, autofit_control = control))
  expect_identical(convergence_arguments$monitor, "mu")
  expect_true(convergence_arguments$allow_not_assessable)
})
