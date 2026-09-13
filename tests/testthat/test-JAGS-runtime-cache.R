skip_if_not_test_profile("unit")

test_that("runtime caches assign each prior process shard once and restore empty pools", {

  expect_error(.JAGS_validate_runtime_cache(list()),
    "'runtime_cache' must be NULL or a function accepting 'context' and 'state' arguments.", fixed = TRUE)
  cluster <- as.list(1:3)
  tasks <- list()
  calls <- list()
  callback <- function(context, state = NULL){

    calls[[length(calls) + 1L]] <<- list(context = context, state = state)
    if(context$phase == "capture") list(process = context$process_id) else NULL
  }
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, callback){
      expect_identical(environment(fun), baseenv())
      tasks <<- x
      lapply(x, fun, callback = callback)
    },
    .package = "parallel"
  )
  states <- lapply(1:4, function(index) list(identity = "fixture", payload = as.raw(index)))
  set.seed(12)
  seed <- .Random.seed
  .JAGS_run_runtime_cache(callback, "restore", 10L, cluster, states)
  expect_identical(lapply(tasks, `[[`, "state"), list(states[c(1L, 4L)], states[2L], states[3L]))
  expect_identical(vapply(tasks, function(task) task$context$process_chains, integer(1L)), c(4L, 3L, 3L))
  expect_equal(sum(lengths(lapply(tasks, `[[`, "state"))), 4L)
  .JAGS_run_runtime_cache(callback, "restore", 10L, cluster, states[1L])
  expect_identical(lapply(tasks, `[[`, "state"), list(states[1L], list(), list()))
  calls <- list()
  .JAGS_run_runtime_cache(callback, "restore", 10L, cluster)
  expect_length(calls, 3L)
  expect_true(all(vapply(calls, function(call) identical(call$state, list()), logical(1L))))
  calls <- list()
  .JAGS_run_runtime_cache(callback, "restore", 10L, state = states)
  expect_identical(calls[[1L]]$state, states)
  expect_identical(calls[[1L]]$context$role, "local")
  expect_warning(.JAGS_run_runtime_cache(callback, "restore", 10L, state = "invalid"),
    "Runtime cache restore is unavailable because the retained state is not a list.", fixed = TRUE)
  expect_identical(tail(calls, 1L)[[1L]]$state, list())
  captured <- .JAGS_run_runtime_cache(callback, "capture", 10L, cluster)
  expect_identical(captured, lapply(1:3, function(index) list(process = index)))
  expect_identical(.Random.seed, seed)
})


test_that("runtime cache warnings and capture failures are returned to the caller", {

  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, callback) lapply(x, fun, callback = callback),
    .package = "parallel"
  )
  callback <- function(context, state = NULL){

    if(context$process_id == 1L) warning("Only part of the optional cache was available.", call. = FALSE)
    if(context$process_id == 2L) stop("Snapshot export failed.", call. = FALSE)
    as.raw(context$process_id)
  }
  messages <- character()
  value <- withCallingHandlers(.JAGS_run_runtime_cache(callback, "capture", 3L, as.list(1:3)),
    warning = function(condition){
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
  expect_identical(value, list(as.raw(1L), NULL, as.raw(3L)))
  expect_identical(messages, c(
    "Runtime cache capture in worker 1: Only part of the optional cache was available.",
    "Runtime cache capture in worker 2: Snapshot export failed."))
  expect_warning(value <- .JAGS_run_runtime_cache(function(context, state = NULL) {
    stop("Snapshot export failed.", call. = FALSE)
  }, "capture", 1L), "Runtime cache capture in the local process: Snapshot export failed.", fixed = TRUE)
  expect_null(value)
})


test_that("saved fit caches restore without entering backend extension payloads", {

  skip_if_not_installed("runjags")
  basic_fit <- structure(list(end.state = rep("", 2L)), class = c("runjags", "BayesTools_fit"))
  attr(basic_fit, "parameter_map") <- .bt_build_parameter_map(character())
  callback <- function(context, state = NULL){

    if(context$phase == "capture") list(process = context$process_id, identity = "fixture") else NULL
  }
  environment(callback) <- baseenv()
  cache_calls <- list()
  events <- character()
  original_cache <- .JAGS_run_runtime_cache
  check_failure <- FALSE
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .bt_attach_fit_contract = function(fit, ...) fit,
    JAGS_check_convergence = function(...){
      if(check_failure) stop("Unexpected convergence failure.")
      TRUE
    },
    .JAGS_run_runtime_cache = function(runtime_cache, phase, chains, cl = NULL, state = NULL){
      events <<- c(events, phase)
      cache_calls[[length(cache_calls) + 1L]] <<- list(phase = phase, state = state)
      original_cache(runtime_cache, phase, chains, cl, state)
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    run.jags = function(...) basic_fit,
    extend.jags = function(runjags.object, ...){
      expect_null(attr(runjags.object, "runtime_state", exact = TRUE))
      expect_null(attr(runjags.object, "runtime_cache", exact = TRUE))
      runjags.object
    },
    add.summary = function(fit, ...) fit,
    .package = "runjags"
  )
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores) as.list(seq_len(cores)),
    stopCluster = function(cl) events <<- c(events, "stop"),
    clusterApply = function(cl, x, fun, ...) lapply(x, fun, ...),
    .package = "parallel"
  )
  control <- list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL, max_time = NULL, sample_extend = 1,
    restarts = 1, max_extend = 1, check_indicators = FALSE)
  fit <- JAGS_fit("model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
    prior_list = list(mu = prior("normal", list(0, 1))), chains = 2,
    adapt = 50, burnin = 50, sample = 100, seed = 1, silent = TRUE,
    runtime_cache = callback)
  expect_identical(events, c("restore", "capture"))
  expect_null(cache_calls[[1L]]$state)
  expect_identical(attr(fit, "runtime_state"), list(list(process = 1L, identity = "fixture")))
  saved <- unserialize(serialize(fit, NULL))
  old_state <- attr(saved, "runtime_state", exact = TRUE)
  events <- character()
  extended <- JAGS_extend(saved, autofit_control = control, parallel = TRUE, cores = 2L)
  expect_identical(events, c("restore", "capture", "stop"))
  expect_identical(tail(cache_calls, 2L)[[1L]]$state, old_state)
  expect_identical(attr(saved, "runtime_state", exact = TRUE), old_state)
  expect_identical(attr(extended, "runtime_state"), list(
    list(process = 1L, identity = "fixture"), list(process = 2L, identity = "fixture")))
  no_capture <- function(context, state = NULL) NULL
  environment(no_capture) <- baseenv()
  events <- character()
  cold <- JAGS_extend(extended, autofit_control = control, runtime_cache = no_capture)
  expect_identical(events, c("restore", "capture"))
  expect_identical(tail(cache_calls, 2L)[[1L]]$state, attr(extended, "runtime_state"))
  expect_null(attr(cold, "runtime_state", exact = TRUE))
  expect_false(is.null(attr(extended, "runtime_state", exact = TRUE)))
  failing <- function(context, state = NULL){
    if(context$phase == "capture") stop("Snapshot export failed.", call. = FALSE)
    NULL
  }
  environment(failing) <- baseenv()
  expect_warning(valid <- JAGS_extend(extended, autofit_control = control,
    runtime_cache = failing), "Runtime cache capture in the local process: Snapshot export failed.", fixed = TRUE)
  expect_s3_class(valid, "BayesTools_fit")
  expect_null(attr(valid, "runtime_state", exact = TRUE))
  events <- character()
  check_failure <- TRUE
  expect_error(JAGS_extend(extended, autofit_control = control, parallel = TRUE, cores = 2L),
    "Unexpected convergence failure.", fixed = TRUE)
  expect_identical(events, c("restore", "stop"))
})
