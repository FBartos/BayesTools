skip_if_not_test_profile("unit")

test_that("parallel JAGS checks the actual workers' package builds", {

  native <- getLoadedDLLs()[["BayesTools"]][["path"]]
  expect_identical(
    unname(.JAGS_package_builds("BayesTools")$BayesTools$dll),
    unname(tools::md5sum(native))
  )
  expected <- .JAGS_package_builds("stats")
  actual <- expected
  testthat::local_mocked_bindings(
    clusterCall = function(cl, fun, packages) list(actual, actual),
    .package = "parallel"
  )
  expect_identical(unname(.JAGS_require_packages("stats", cl = list())), TRUE)
  message <- paste0(
    "Parallel JAGS fitting is unavailable with mismatched package versions, R code, or native builds: 'stats'. ",
    "Install the parent-session builds into a library and set 'R_LIBS_USER' ",
    "to that library before starting R and its workers."
  )
  actual$stats$version <- "0.0.0"
  expect_error(.JAGS_require_packages("stats", cl = list()), message, fixed = TRUE)
  actual <- expected
  actual$stats$dll[] <- "different-build"
  expect_error(.JAGS_require_packages("stats", cl = list()), message, fixed = TRUE)
  actual <- expected
  actual$stats$r_code <- "different-code"
  expect_error(.JAGS_require_packages("stats", cl = list()), message, fixed = TRUE)
  expect_error(.JAGS_require_packages("stats", cl = list(),
    operation = "Parallel zplot density computation"),
    sub("Parallel JAGS fitting", "Parallel zplot density computation", message, fixed = TRUE),
    fixed = TRUE)
  actual$stats <- NULL
  expect_error(
    .JAGS_require_packages("stats", cl = list()),
    "Required packages are not available: 'stats'.", fixed = TRUE
  )
})

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
  attr(fit, "parameter_map") <- .bt_build_parameter_map(character())
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

test_that("JAGS_extend rejects stale fitted metadata before backend work", {

  stale_contract <- .bt_attach_fit_contract(.jags_extend_test_fit())
  attr(stale_contract, "fit_contract")$formula_design_version <- 2L
  package_calls <- 0L
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...){
      package_calls <<- package_calls + 1L
    },
    .package = "BayesTools"
  )

  expect_error(
    JAGS_extend(
      stale_contract,
      autofit_control = .jags_extend_test_control()
    ),
    "missing or unsupported 'formula_design' metadata",
    fixed = TRUE
  )
  expect_identical(package_calls, 0L)

  formula_result <- JAGS_formula(
    ~ 1,
    "mu",
    data.frame(row = 1:2),
    list(intercept = prior("normal", list(0, 1)))
  )
  stale_design <- formula_result$formula_design
  stale_design$schema_version <- 2L
  actual_mismatch <- .jags_extend_test_fit()
  attr(actual_mismatch, "formula_design") <- list(mu = stale_design)
  actual_mismatch <- .bt_attach_fit_contract(actual_mismatch)

  expect_error(
    JAGS_extend(
      actual_mismatch,
      autofit_control = .jags_extend_test_control()
    ),
    "cannot replay this fitted formula",
    fixed = TRUE
  )
  expect_identical(package_calls, 0L)
})

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
  expected_fit <- fit
  attr(expected_fit, "warnings") <-
    "The model extension failed; returning the last valid fit. Backend error: backend exploded"
  expect_identical(result, expected_fit)
  expect_identical(
    attr(result, "warnings"),
    "The model extension failed; returning the last valid fit. Backend error: backend exploded"
  )
})

test_that("JAGS_fit autofit preserves the last valid fit after a backend error", {

  skip_if_not_installed("runjags")
  initial_fit <- structure(
    list(mcmc = list(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "mu")))),
    class = "runjags"
  )
  testthat::local_mocked_bindings(
    run.jags = function(...) initial_fit,
    extend.jags = function(...) stop("backend exploded"),
    add.summary = function(x, ...) x,
    .package = "runjags"
  )
  testthat::local_mocked_bindings(
    JAGS_check_convergence = function(...) FALSE,
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )

  expect_warning(
    result <- JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      autofit = TRUE,
      autofit_control = list(
        max_Rhat = NULL,
        min_ESS = NULL,
        max_error = NULL,
        max_SD_error = NULL,
        max_time = list(time = 60, unit = "secs"),
        sample_extend = 1,
        restarts = 1,
        max_extend = 1,
        check_indicators = FALSE
      ),
      silent = TRUE,
      seed = 1
    ),
    "returning the last valid fit.*backend exploded"
  )
  expect_s3_class(result, "BayesTools_fit")
  expect_false(inherits(result, "error"))
  expect_identical(result$mcmc, initial_fit$mcmc)
  expect_identical(
    attr(result, "warnings"),
    "The model extension failed; returning the last valid fit. Backend error: backend exploded"
  )
})

test_that("JAGS build checks fingerprint loaded R definitions", {

  original <- .JAGS_package_builds("stats")$stats
  sd <- stats::sd
  body(sd) <- quote(stop("modified implementation"))
  testthat::local_mocked_bindings(sd = sd, .package = "stats")
  changed <- .JAGS_package_builds("stats")$stats
  expect_identical(changed$version, original$version)
  expect_identical(changed$dll, original$dll)
  expect_false(identical(changed$r_code, original$r_code))
})

test_that("JAGS_fit reports and records backend errors recovered by a restart", {

  skip_if_not_installed("runjags")
  backend_calls <- 0L
  recovered_fit <- structure(
    list(mcmc = list(matrix(0, nrow = 2, ncol = 1,
                            dimnames = list(NULL, "mu")))),
    class = "runjags"
  )
  testthat::local_mocked_bindings(
    run.jags = function(...){
      backend_calls <<- backend_calls + 1L
      if(backend_calls == 1L){
        stop("backend exploded")
      }
      recovered_fit
    },
    .package = "runjags"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .package = "BayesTools"
  )

  fit_once <- function(silent){

    JAGS_fit(
      model_syntax = "model{ mu ~ dnorm(0, 1) }",
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      autofit_control = list(
        max_Rhat = NULL,
        min_ESS = NULL,
        max_error = NULL,
        max_SD_error = NULL,
        max_time = list(time = 60, unit = "secs"),
        sample_extend = 1,
        restarts = 2,
        max_extend = 1,
        check_indicators = FALSE
      ),
      silent = silent,
      seed = 1
    )
  }
  expected_message <-
    "JAGS fitting attempt 1 failed and was restarted: backend exploded."

  for(silent in c(FALSE, TRUE)){
    backend_calls <- 0L
    emitted <- character()
    result <- withCallingHandlers(fit_once(silent), warning = function(w){
      expect_identical(backend_calls, 1L)
      expect_null(conditionCall(w))
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

    expect_identical(emitted, if(silent) character() else expected_message)
    expect_identical(backend_calls, 2L)
    expect_s3_class(result, "BayesTools_fit")
    expect_identical(attr(result, "warnings"), expected_message)
  }
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


test_that("loaded numerical constants and immutable attributes enter build parity", {

  original <- .JAGS_package_builds("BayesTools")$BayesTools
  nested <- structure(list(orders = list(c(3L, 7L), c(15L, 31L))),
    source = list(kind = "quadrature", version = 1L))
  testthat::local_mocked_bindings(.bt_parameter_map_version = nested, .package = "BayesTools")
  changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_identical(changed$version, original$version)
  expect_identical(changed$dll, original$dll)
  expect_false(identical(changed$r_code, original$r_code))
  nested$orders[[2L]][2L] <- 63L
  testthat::local_mocked_bindings(.bt_parameter_map_version = nested, .package = "BayesTools")
  value_changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_false(identical(value_changed$r_code, changed$r_code))
  attr(nested, "source")$version <- 2L
  testthat::local_mocked_bindings(.bt_parameter_map_version = nested, .package = "BayesTools")
  attribute_changed <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_false(identical(attribute_changed$r_code, value_changed$r_code))

  cache <- new.env(parent = emptyenv())
  cache$value <- 1
  mutable <- structure(list(orders = c(3L, 7L)), metadata = list(cache = cache))
  testthat::local_mocked_bindings(.bt_parameter_map_version = mutable, .package = "BayesTools")
  first <- .JAGS_package_builds("BayesTools")$BayesTools
  cache$value <- 2
  second <- .JAGS_package_builds("BayesTools")$BayesTools
  expect_identical(second, first)
})


test_that("post-fit cleanup errors identify the actual parallel operation", {

  testthat::local_mocked_bindings(
    stopCluster = function(cl) stop("socket is unavailable"),
    .package = "parallel"
  )
  expect_warning(.JAGS_finish_runtime_setup(
    function(context) stop("finish must not run"), chains = 2L, cl = list(),
    operation = "Parallel zplot density computation"),
    paste0("Parallel zplot density computation worker cleanup failed: socket is unavailable. ",
      "The runtime finish callback was not run."), fixed = TRUE)
})
