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

  basic_fit <- contract_test_backend_fit(2L)
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
    .bt_check_convergence = function(...){
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
    .bt_rjags_fit = function(...) basic_fit,
    .bt_rjags_extend = function(fit, ...){
      expect_null(attr(fit, "runtime_state", exact = TRUE))
      expect_null(attr(fit, "runtime_cache", exact = TRUE))
      fit
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores) as.list(seq_len(cores)),
    clusterCall = function(cl, fun, ...) lapply(cl, function(worker) fun(...)),
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
  expect_identical(events, c("restore", "capture", "stop"))
  fail_extend <- TRUE
  testthat::local_mocked_bindings(
    .bt_rjags_extend = function(...) stop("JAGS model compilation failed"),
    .package = "BayesTools"
  )
  events <- character()
  expect_warning(
    failed <- JAGS_extend(extended, autofit_control = control, parallel = TRUE, cores = 2L),
    "The model extension failed; returning the last valid fit.",
    fixed = TRUE
  )
  expect_identical(events, c("restore", "stop"))
  expect_identical(attr(failed, "runtime_state", exact = TRUE),
    attr(extended, "runtime_state", exact = TRUE))

  # The initial successful fit also needs a snapshot before an extension
  # fails; capturing afterward could save the failed backend's partial state.
  check_failure <- FALSE
  testthat::local_mocked_bindings(
    .bt_check_convergence = function(...) FALSE,
    .package = "BayesTools"
  )
  events <- character()
  expect_warning(
    failed_first <- JAGS_fit(
      "model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
      prior_list = list(mu = prior("normal", list(0, 1))), chains = 2,
      adapt = 50, burnin = 50, sample = 100, seed = 1, silent = TRUE,
      autofit = TRUE, autofit_control = control, runtime_cache = callback
    ),
    "The model extension failed; returning the last valid fit.",
    fixed = TRUE
  )
  expect_identical(events, c("restore", "capture"))
  expect_identical(attr(failed_first, "runtime_state", exact = TRUE),
    list(list(process = 1L, identity = "fixture")))
})

test_that("runtime cache capture warnings are recorded in the fit warnings", {

  basic_fit <- contract_test_backend_fit(2L)
  worker_alive <- TRUE
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .bt_attach_fit_contract = function(fit, ...) fit,
    .bt_check_convergence = function(...) TRUE,
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...) basic_fit,
    .bt_rjags_extend = function(fit, ...) fit,
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores) as.list(seq_len(cores)),
    clusterCall = function(cl, fun, ...) lapply(cl, function(worker) fun(...)),
    stopCluster = function(cl) invisible(NULL),
    clusterApply = function(cl, x, fun, ...){
      if(!worker_alive && identical(x[[1L]]$context$phase, "capture")){
        stop("error reading from connection", call. = FALSE)
      }
      lapply(x, fun, ...)
    },
    .package = "parallel"
  )
  failing <- function(context, state = NULL){
    if(context$phase == "capture") stop("Snapshot export failed.", call. = FALSE)
    NULL
  }
  environment(failing) <- baseenv()
  control <- list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL, max_time = NULL, sample_extend = 1,
    restarts = 1, max_extend = 1, check_indicators = FALSE)

  local_message <- "Runtime cache capture in the local process: Snapshot export failed."
  expect_warning(
    fit <- JAGS_fit("model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
      prior_list = list(mu = prior("normal", list(0, 1))), chains = 2,
      adapt = 50, burnin = 50, sample = 100, seed = 1, silent = TRUE,
      runtime_cache = failing),
    local_message, fixed = TRUE
  )
  expect_null(attr(fit, "runtime_state", exact = TRUE))
  expect_identical(attr(fit, "warnings"), local_message)

  expect_warning(
    extended <- JAGS_extend(fit, autofit_control = control, runtime_cache = failing),
    local_message, fixed = TRUE
  )
  expect_identical(attr(extended, "warnings"), local_message)

  # A dead worker drops the whole capture; the warning stays with the fit.
  keeping <- function(context, state = NULL){
    if(context$phase == "capture") list(process = context$process_id) else NULL
  }
  environment(keeping) <- baseenv()
  worker_alive <- FALSE
  dead_message <- "Runtime cache capture could not be completed: error reading from connection"
  expect_warning(
    parallel_fit <- JAGS_extend(fit, autofit_control = control,
      parallel = TRUE, cores = 2L, runtime_cache = keeping),
    dead_message, fixed = TRUE
  )
  expect_null(attr(parallel_fit, "runtime_state", exact = TRUE))
  expect_identical(attr(parallel_fit, "warnings"), c(local_message, dead_message))
})

# A mocked backend for runtime reports. 'state' scripts the backend failures
# (the numbers of the failing run.jags and extend.jags calls and their
# messages, and the number of unconverged checks) and records the calls and the
# objects passed to the backend. A test mocks the backend once (the mocked
# fit-contract helpers would leave a second backend fit without a contract) and
# scripts each scenario on the same state.
.report_test_script <- function(state, fail_runs = integer(), fail_extensions = integer(),
                                run_error = "retry probe",
                                extend_error = "JAGS model compilation failed",
                                unconverged = 0L, reset_reports = TRUE){

  state$fail_runs <- fail_runs
  state$fail_extensions <- fail_extensions
  state$run_error <- run_error
  state$extend_error <- extend_error
  state$unconverged <- unconverged
  state$runs <- state$extensions <- state$checks <- state$cluster_calls <- 0L
  state$objects <- list()
  # read-and-clear state of each process that every backend call leaves
  state$pending <- c(0L, 0L)
  if(reset_reports){
    state$reports <- 0L
    state$contexts <- list()
  }
  invisible(state)
}

.report_test_state <- function(...){

  .report_test_script(new.env(parent = emptyenv()), ...)
}

.report_test_backend <- function(state, env = parent.frame()){

  basic_fit <- contract_test_backend_fit(2L)
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .bt_attach_fit_contract = function(fit, ...) fit,
    .bt_check_convergence = function(...){
      state$checks <- state$checks + 1L
      state$checks > state$unconverged
    },
    .JASP_progress_bar_start = function(...) invisible(NULL),
    .JASP_progress_bar_tick = function(...) invisible(NULL),
    .package = "BayesTools", .env = env
  )
  testthat::local_mocked_bindings(
    .bt_rjags_fit = function(...){
      state$runs <- state$runs + 1L
      state$pending <- state$pending + 1L
      if(state$runs %in% state$fail_runs) stop(state$run_error, call. = FALSE)
      basic_fit
    },
    .bt_rjags_extend = function(fit, ...){
      state$extensions <- state$extensions + 1L
      state$pending <- state$pending + 1L
      state$objects <- c(state$objects, list(fit))
      if(state$extensions %in% state$fail_extensions) stop(state$extend_error, call. = FALSE)
      fit
    },
    .package = "BayesTools", .env = env
  )
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores, ...) as.list(seq_len(cores)),
    clusterCall = function(cl, fun, ...) lapply(cl, function(worker) fun(...)),
    stopCluster = function(cl) invisible(NULL),
    clusterApply = function(cl, x, fun, callback){
      expect_identical(environment(fun), baseenv())
      state$cluster_calls <- state$cluster_calls + 1L
      lapply(x, fun, callback = callback)
    },
    .package = "parallel", .env = env
  )
  invisible(state)
}

# A report callback returning the process label and a call counter.
.report_test_callback <- function(state){

  force(state)
  function(context){

    state$reports <- state$reports + 1L
    state$contexts <- c(state$contexts, list(context))
    list(process = paste(context$role, context$process_id), call = state$reports)
  }
}

.report_test_entry <- function(process, call) list(process = process, call = call)

# A report callback that reads and clears the state that the backend calls
# since its previous call left in its process.
.report_test_reader <- function(state){

  force(state)
  function(context){

    value <- state$pending[[context$process_id]]
    state$pending[[context$process_id]] <- 0L
    value
  }
}

.report_test_fit <- function(..., chains = 2){

  JAGS_fit("model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
    prior_list = list(mu = prior("normal", list(0, 1))), chains = chains,
    adapt = 50, burnin = 50, sample = 100, seed = 1, silent = TRUE, ...)
}

.report_test_control <- function(restarts = 1, max_extend = 1){

  list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL, max_time = NULL, sample_extend = 1,
    restarts = restarts, max_extend = max_extend, check_indicators = FALSE)
}


test_that("runtime reports collect one element per process of every kept run", {

  entry <- .report_test_entry
  state <- .report_test_backend(.report_test_state())
  report <- .report_test_callback(state)

  # serial fit: one run, the calling process
  fit <- .report_test_fit(runtime_report = report)
  expect_identical(attr(fit, "runtime_report_values"), list(list(entry("local 1", 1L))))
  expect_identical(attr(fit, "runtime_report"), report)
  expect_identical(state$contexts[[1L]], list(phase = "report", role = "local",
    chains = 2L, processes = 1L, process_id = 1L, process_chains = 2L, parallel = FALSE))
  expect_identical(state$cluster_calls, 0L)

  # parallel fit: one element per worker, in worker order
  .report_test_script(state)
  fit <- .report_test_fit(runtime_report = report, chains = 3, parallel = TRUE, cores = 2)
  expect_identical(attr(fit, "runtime_report_values"),
    list(list(entry("worker 1", 1L), entry("worker 2", 2L))))
  expect_identical(vapply(state$contexts, `[[`, integer(1L), "process_chains"), c(2L, 1L))
  expect_identical(state$cluster_calls, 1L)

  # autofit: the initial run and each extension, serial and parallel
  for(parallel in c(FALSE, TRUE)){
    .report_test_script(state, unconverged = 2L)
    fit <- .report_test_fit(runtime_report = report, autofit = TRUE,
      autofit_control = .report_test_control(max_extend = 5), parallel = parallel, cores = 2)
    expect_identical(state$extensions, 2L)
    runs <- attr(fit, "runtime_report_values")
    expect_length(runs, 3L)
    if(parallel){
      expect_identical(runs[[3L]], list(entry("worker 1", 5L), entry("worker 2", 6L)))
    }else{
      expect_identical(runs, lapply(1:3, function(call) list(entry("local 1", call))))
    }
  }

  # a callback that returns NULL keeps no report but is stored and called
  .report_test_script(state, unconverged = 1L)
  silent_report <- function(context){
    state$reports <- state$reports + 1L
    NULL
  }
  fit <- .report_test_fit(runtime_report = silent_report, autofit = TRUE,
    autofit_control = .report_test_control(), parallel = TRUE, cores = 2)
  expect_identical(state$reports, 4L)
  expect_false("runtime_report_values" %in% names(attributes(fit)))
  expect_null(attr(fit, "runtime_report_values"))
  expect_identical(attr(fit, "runtime_report"), silent_report)
})

test_that("runtime reports drop failed attempts and failed extensions but still call the callback", {

  entry <- .report_test_entry
  state <- .report_test_backend(.report_test_state())
  report <- .report_test_callback(state)

  # a backend failure followed by a restart: the failed attempt's value is
  # absent, and the counter shows that the callback ran after it
  for(parallel in c(FALSE, TRUE)){
    .report_test_script(state, fail_runs = 1L)
    fit <- .report_test_fit(runtime_report = report, parallel = parallel, cores = 2,
      autofit_control = .report_test_control(restarts = 2))
    expect_identical(state$runs, 2L)
    expected <- if(parallel){
      list(list(entry("worker 1", 3L), entry("worker 2", 4L)))
    }else{
      list(list(entry("local 1", 2L)))
    }
    expect_identical(attr(fit, "runtime_report_values"), expected)
    expect_identical(attr(fit, "warnings"),
      "JAGS fitting attempt 1 failed and was restarted: retry probe.")
  }

  # every attempt failed: an error fit keeps the callback but carries no report
  .report_test_script(state, fail_runs = 1:2)
  failed <- .report_test_fit(runtime_report = report,
    autofit_control = .report_test_control(restarts = 2))
  expect_s3_class(failed, "error")
  expect_identical(state$reports, 2L)
  expect_false("runtime_report_values" %in% names(attributes(failed)))
  expect_null(attr(failed, "runtime_report_values"))
  expect_identical(attr(failed, "runtime_report"), report)

  # the callback reads the state that each backend attempt leaves (one
  # increment per run.jags or extend.jags call, read and cleared): every kept
  # run reads exactly one attempt, so the callback also ran after the failed
  # initial attempt and the failed extension, whose reads were dropped
  reader <- .report_test_reader(state)
  for(parallel in c(FALSE, TRUE)){
    .report_test_script(state, fail_runs = 1L, fail_extensions = 2L, unconverged = 10L)
    expect_warning(
      read <- .report_test_fit(runtime_report = reader, autofit = TRUE, parallel = parallel,
        cores = 2, autofit_control = .report_test_control(restarts = 2, max_extend = 5)),
      "The model extension failed; returning the last valid fit.", fixed = TRUE)
    expect_identical(c(state$runs, state$extensions), c(2L, 2L))
    expect_identical(attr(read, "runtime_report_values"),
      rep(list(rep(list(1L), if(parallel) 2L else 1L)), 2L))
    expect_identical(state$pending, c(0L, if(parallel) 0L else 4L))
    .report_test_script(state, fail_extensions = 2L, unconverged = 10L)
    expect_warning(
      read_extended <- JAGS_extend(read, autofit_control = .report_test_control(max_extend = 3),
        parallel = parallel, cores = 2),
      "The model extension failed; returning the last valid fit.", fixed = TRUE)
    expect_identical(attr(read_extended, "runtime_report_values"),
      list(rep(list(1L), if(parallel) 2L else 1L)))
    expect_identical(state$pending, c(0L, if(parallel) 0L else 2L))
  }

  # a failed automatic extension: the runs before it are kept, its own is not
  .report_test_script(state, fail_extensions = 2L, unconverged = 10L)
  expect_warning(
    fit <- .report_test_fit(runtime_report = report, autofit = TRUE,
      autofit_control = .report_test_control(max_extend = 5)),
    "The model extension failed; returning the last valid fit.", fixed = TRUE)
  expect_identical(state$reports, 3L)
  expect_identical(attr(fit, "runtime_report_values"),
    list(list(entry("local 1", 1L)), list(entry("local 1", 2L))))

  # the same in JAGS_extend(): a call without a successful extension returns
  # no report (the input's report is not carried over); a later failure keeps
  # the runs before it
  .report_test_script(state, fail_extensions = 1L, unconverged = 10L, reset_reports = FALSE)
  expect_warning(
    none <- JAGS_extend(fit, autofit_control = .report_test_control(max_extend = 3)),
    "The model extension failed; returning the last valid fit.", fixed = TRUE)
  expect_identical(state$reports, 4L)
  expect_false("runtime_report_values" %in% names(attributes(none)))
  expect_null(attr(none, "runtime_report_values"))
  expect_identical(attr(none, "runtime_report"), report)
  .report_test_script(state, fail_extensions = 2L, unconverged = 10L, reset_reports = FALSE)
  expect_warning(
    partial <- JAGS_extend(fit, autofit_control = .report_test_control(max_extend = 3)),
    "The model extension failed; returning the last valid fit.", fixed = TRUE)
  expect_identical(state$reports, 6L)
  expect_identical(attr(partial, "runtime_report_values"), list(list(entry("local 1", 5L))))
  # neither report attribute travels in the object passed to the backend
  expect_length(state$objects, 2L)
  for(object in state$objects){
    expect_null(attr(object, "runtime_report_values", exact = TRUE))
    expect_null(attr(object, "runtime_report", exact = TRUE))
  }
})

test_that("runtime reports are not requested from workers after a connection failure", {

  entry <- .report_test_entry
  connection <- "error reading from connection"
  state <- .report_test_backend(.report_test_state())
  report <- .report_test_callback(state)

  # the initial run: no report, no retry, an error fit without a report
  .report_test_script(state, fail_runs = 1L, run_error = connection)
  failed <- .report_test_fit(runtime_report = report, parallel = TRUE, cores = 2,
    autofit_control = .report_test_control(restarts = 3))
  expect_s3_class(failed, "BayesTools_JAGS_worker_connection_error")
  expect_identical(c(state$runs, state$cluster_calls, state$reports), c(1L, 0L, 0L))
  expect_false("runtime_report_values" %in% names(attributes(failed)))
  expect_identical(attr(failed, "runtime_report"), report)

  # an automatic extension and a JAGS_extend() extension: the kept runs remain
  .report_test_script(state, fail_extensions = 1L, extend_error = connection, unconverged = 10L)
  expect_warning(
    fit <- .report_test_fit(runtime_report = report, autofit = TRUE, parallel = TRUE,
      cores = 2, autofit_control = .report_test_control(max_extend = 3)),
    "Parallel JAGS worker communication failed", fixed = TRUE)
  expect_identical(c(state$extensions, state$cluster_calls), c(1L, 1L))
  expect_identical(attr(fit, "runtime_report_values"),
    list(list(entry("worker 1", 1L), entry("worker 2", 2L))))
  .report_test_script(state, fail_extensions = 1L, extend_error = connection, unconverged = 10L)
  expect_warning(
    extended <- JAGS_extend(fit, autofit_control = .report_test_control(),
      parallel = TRUE, cores = 2),
    "Parallel JAGS worker communication failed", fixed = TRUE)
  expect_identical(c(state$extensions, state$cluster_calls, state$reports), c(1L, 0L, 0L))
  expect_false("runtime_report_values" %in% names(attributes(extended)))
})

test_that("JAGS_extend replays the stored report callback, accepts a new one, and NULL disables it", {

  entry <- .report_test_entry
  state <- .report_test_backend(.report_test_state(unconverged = 10L))
  report <- .report_test_callback(state)
  fit <- .report_test_fit(runtime_report = report)

  # the stored callback, with the current topology; only this call's runs
  extended <- JAGS_extend(fit, autofit_control = .report_test_control(max_extend = 2),
    parallel = TRUE, cores = 2)
  expect_identical(attr(extended, "runtime_report"), report)
  expect_identical(attr(extended, "runtime_report_values"), list(
    list(entry("worker 1", 2L), entry("worker 2", 3L)),
    list(entry("worker 1", 4L), entry("worker 2", 5L))))
  expect_identical(attr(fit, "runtime_report_values"), list(list(entry("local 1", 1L))))
  expect_length(state$objects, 2L)
  for(object in state$objects){
    expect_null(attr(object, "runtime_report_values", exact = TRUE))
    expect_null(attr(object, "runtime_report", exact = TRUE))
  }

  # a saved and reloaded fit replays its stored callback
  saved <- unserialize(serialize(fit, NULL))
  reloaded <- JAGS_extend(saved, autofit_control = .report_test_control())
  expect_identical(lengths(attr(reloaded, "runtime_report_values")), 1L)
  expect_identical(attr(reloaded, "runtime_report_values")[[1L]][[1L]]$process, "local 1")

  # a new callback replaces the stored one
  other <- function(context) "other"
  replaced <- JAGS_extend(fit, autofit_control = .report_test_control(), runtime_report = other)
  expect_identical(attr(replaced, "runtime_report_values"), list(list("other")))
  expect_identical(attr(replaced, "runtime_report"), other)

  # an explicit NULL disables it and removes the input's report
  reports <- state$reports
  disabled <- JAGS_extend(fit, autofit_control = .report_test_control(),
    parallel = TRUE, cores = 2, runtime_report = NULL)
  expect_identical(state$reports, reports)
  expect_false(any(c("runtime_report_values", "runtime_report") %in% names(attributes(disabled))))
  expect_error(JAGS_extend(fit, autofit_control = .report_test_control(), runtime_report = "report"),
    "'runtime_report' must be NULL or a function accepting one context argument.", fixed = TRUE)
})

test_that("runtime report warnings and errors are recorded for kept runs and never lose the fit", {

  collect <- function(expr){
    messages <- character()
    value <- withCallingHandlers(expr, warning = function(condition){
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    list(value = value, messages = messages)
  }
  state <- .report_test_backend(.report_test_state())

  # warnings accumulate over the kept runs; a dropped attempt's warning is
  # signalled but not recorded
  .report_test_script(state, fail_runs = 1L, unconverged = 1L)
  warning_report <- function(context){
    state$reports <- state$reports + 1L
    warning("Counter read ", state$reports, ".", call. = FALSE)
    state$reports
  }
  run <- collect(.report_test_fit(runtime_report = warning_report, autofit = TRUE,
    autofit_control = .report_test_control(restarts = 2)))
  expect_identical(run$messages, paste0("Runtime report in the local process: Counter read ", 1:3, "."))
  expect_identical(attr(run$value, "runtime_report_values"), list(list(2L), list(3L)))
  expect_identical(attr(run$value, "warnings"), c(
    "JAGS fitting attempt 1 failed and was restarted: retry probe.",
    paste0("Runtime report in the local process: Counter read ", 2:3, ".")))

  # an error leaves NULL for its process and is recorded once
  .report_test_script(state, unconverged = 1L)
  error_report <- function(context){
    if(context$process_id == 2L || !context$parallel) stop("Counter read failed.", call. = FALSE)
    context$process_id
  }
  run <- collect(.report_test_fit(runtime_report = error_report, autofit = TRUE,
    autofit_control = .report_test_control()))
  expect_s3_class(run$value, "BayesTools_fit")
  expect_false("runtime_report_values" %in% names(attributes(run$value)))
  expect_identical(run$messages, rep("Runtime report in the local process: Counter read failed.", 2L))
  expect_identical(attr(run$value, "warnings"), "Runtime report in the local process: Counter read failed.")
  run <- collect(.report_test_fit(runtime_report = error_report, parallel = TRUE, cores = 2))
  expect_identical(attr(run$value, "runtime_report_values"), list(list(1L, NULL)))
  expect_identical(attr(run$value, "warnings"), "Runtime report in worker 2: Counter read failed.")

  # a failed transport of the report: NULL per process and one warning
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, callback) stop("error writing to connection", call. = FALSE),
    .package = "parallel"
  )
  run <- collect(.report_test_fit(runtime_report = error_report, parallel = TRUE, cores = 2))
  expect_s3_class(run$value, "BayesTools_fit")
  expect_false("runtime_report_values" %in% names(attributes(run$value)))
  expect_identical(run$messages, "Runtime report could not be completed: error writing to connection")
  expect_identical(attr(run$value, "warnings"), run$messages)
})

test_that("fits without a runtime report are unchanged and never contact the workers for one", {

  state <- .report_test_backend(.report_test_state(unconverged = 1L))
  # the attributes of the base commit (58f0495) for this mocked fit and its
  # extension, with and without an explicit NULL
  base_attributes <- c("names", "class", "parameter_map", "draw_geometry", "fit_contract",
    "prior_list", "model_syntax", "add_parameters", "required_packages", "jags_modules",
    "adaptation")
  for(explicit in c(FALSE, TRUE)){
    .report_test_script(state, unconverged = 1L)
    fit <- if(explicit){
      .report_test_fit(autofit = TRUE, autofit_control = .report_test_control(),
        parallel = TRUE, cores = 2, runtime_report = NULL)
    }else{
      .report_test_fit(autofit = TRUE, autofit_control = .report_test_control(),
        parallel = TRUE, cores = 2)
    }
    expect_identical(names(attributes(fit)), base_attributes)
    extended <- JAGS_extend(fit, autofit_control = .report_test_control(), parallel = TRUE, cores = 2)
    expect_identical(names(attributes(extended)), base_attributes)
    expect_identical(c(state$extensions, state$cluster_calls), c(2L, 0L))
  }
})

test_that("JASP ticks follow completed work on both chain groups", {

  completed <- c(NA_real_, NA_real_)
  observed <- list()
  worker <- 0L
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, ...){
      lapply(seq_along(x), function(i){
        worker <<- i
        fun(x[[i]], ...)
      })
    }, .package = "parallel"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_group = function(inits, ...){
      expect_length(inits, if(worker == 1L) 2L else 1L)
      completed[worker] <<- 0
    },
    .bt_rjags_batch = function(sample, ...) completed[worker] <<- completed[worker] + sample,
    .bt_rjags_merge_workers = function(...) list(),
    .JASP_progress_bar_start = function(...) NULL,
    .JASP_progress_bar_tick = function() observed[[length(observed) + 1L]] <<- completed,
    .package = "BayesTools"
  )
  .bt_rjags_run("", list(), list(1L, 2L, 3L), "x", 50L, 50L, 103L, 3L,
    list(1L, 2L), TRUE, jasp_prefix = "")
  expect_identical(observed, lapply(c(0, 21, 42, 63, 83, 103), rep, 2L))
})

test_that("JASP progress preserves one runtime report per fit or autofit extension", {

  state <- .report_test_backend(.report_test_state())
  report <- .report_test_callback(state)
  for(parallel in c(FALSE, TRUE)){
    .report_test_script(state, unconverged = 1L)
    fit <- .report_test_fit(runtime_report = report, parallel = parallel, cores = 2,
      is_JASP = TRUE, autofit = TRUE, autofit_control = .report_test_control())
    expect_identical(c(state$runs, state$extensions), c(1L, 1L))
    runs <- attr(fit, "runtime_report_values")
    if(parallel){
      expect_identical(runs, list(
        list(.report_test_entry("worker 1", 1L), .report_test_entry("worker 2", 2L)),
        list(.report_test_entry("worker 1", 3L), .report_test_entry("worker 2", 4L))
      ))
    }else{
      expect_identical(runs, list(
        list(.report_test_entry("local 1", 1L)), list(.report_test_entry("local 1", 2L))
      ))
    }
  }

  .report_test_script(state, fail_runs = 1L)
  failed <- .report_test_fit(runtime_report = report, is_JASP = TRUE,
    autofit_control = .report_test_control())
  expect_s3_class(failed, "error")
  expect_identical(state$reports, 1L)
  expect_null(attr(failed, "runtime_report_values", exact = TRUE))
  expect_identical(attr(failed, "runtime_report"), report)
})
