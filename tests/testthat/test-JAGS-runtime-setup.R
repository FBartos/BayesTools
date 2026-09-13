skip_if_not_test_profile("unit")

test_that("worker output paths are validated without creating logs", {

  expect_null(.JAGS_validate_worker_output(NULL))
  for(value in list("", NA_character_, c("first.log", "second.log"), 1, list("log"))){
    expect_error(.JAGS_validate_worker_output(value), "'worker_output'", fixed = TRUE)
  }
  missing_parent <- file.path(tempfile("missing-worker-log-directory-"), "worker.log")
  expect_error(.JAGS_validate_worker_output(missing_parent), "'worker_output'", fixed = TRUE)
  expect_error(.JAGS_validate_worker_output(tempdir()), "'worker_output'", fixed = TRUE)
  expect_false(dir.exists(dirname(missing_parent)))
  log_path <- tempfile("worker output ", fileext = ".log")
  normalized <- file.path(normalizePath(dirname(log_path), winslash = "/"), basename(log_path))
  expect_identical(.JAGS_validate_worker_output(log_path), normalized)
  expect_false(file.exists(log_path))

  cluster <- structure(list(1L, 2L), class = "runtime_test_cluster")
  arguments <- list()
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores, ...){
      expect_identical(cores, 2L)
      arguments[[length(arguments) + 1L]] <<- list(...)
      cluster
    },
    .package = "parallel"
  )
  expect_identical(.JAGS_make_cluster(2L), cluster)
  expect_identical(.JAGS_make_cluster(2L, worker_output = normalized), cluster)
  expect_identical(arguments, list(list(), list(outfile = normalized)))
  expect_false(file.exists(log_path))
})

test_that("runtime topology reaches every process and cleanup follows worker exit", {

  message <- "'runtime_setup' must be NULL or a function accepting one context argument."
  expect_error(.JAGS_validate_runtime_setup(list()), message, fixed = TRUE)
  expect_invisible(.JAGS_run_runtime_setup(NULL, chains = 10L))
  expect_error(.JAGS_run_runtime_setup(function(context) stop("setup failed"), 10L),
    "setup failed", fixed = TRUE)

  contexts <- list()
  events <- character()
  setup <- function(context){

    contexts[[length(contexts) + 1L]] <<- context
    events <<- c(events, paste(context$phase, context$role))
    invisible(NULL)
  }
  cluster <- structure(as.list(seq_len(3L)), class = "runtime_test_cluster")
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, setup){
      expect_identical(cl, cluster)
      expect_identical(environment(fun), baseenv())
      lapply(x, fun, setup = setup)
    },
    stopCluster = function(cl) events <<- c(events, "stop"),
    .package = "parallel"
  )
  set.seed(5)
  original_rng <- .Random.seed
  .JAGS_run_runtime_setup(setup, chains = 10L)
  expect_identical(contexts[[1L]], list(phase = "start", role = "local",
    chains = 10L, processes = 1L, process_id = 1L, process_chains = 10L,
    parallel = FALSE))
  .JAGS_run_runtime_setup(setup, chains = 10L, cl = cluster)
  expect_identical(contexts[[2L]], list(phase = "start", role = "coordinator",
    chains = 10L, processes = 3L, process_id = 0L, process_chains = 0L,
    parallel = TRUE))
  expect_identical(vapply(contexts[3:5], `[[`, integer(1), "process_id"), 1:3)
  expect_identical(vapply(contexts[3:5], `[[`, integer(1), "process_chains"),
    c(4L, 3L, 3L))
  expect_identical(sum(vapply(contexts[3:5], `[[`, integer(1), "process_chains")), 10L)
  .JAGS_finish_runtime_setup(setup, chains = 10L, cl = cluster)
  expect_identical(tail(events, 2L), c("stop", "finish coordinator"))
  expect_identical(contexts[[6L]]$phase, "finish")
  expect_identical(.Random.seed, original_rng)

  testthat::local_mocked_bindings(
    stopCluster = function(cl) stop("worker socket is unavailable"),
    .package = "parallel"
  )
  expect_warning(.JAGS_finish_runtime_setup(setup, chains = 10L, cl = cluster),
    paste0("Parallel JAGS worker cleanup failed: worker socket is unavailable. ",
      "The runtime finish callback was not run."), fixed = TRUE)
  expect_length(contexts, 6L)
})

test_that("backend transport classification preserves the original condition", {

  backend <- function(...) list(...)
  arguments <- list(value = 3, monitor = c("mu", "tau"))
  expect_identical(.JAGS_run_backend(backend, arguments, parallel = TRUE), arguments)

  for(message in c("error reading from connection", "error writing to connection",
      "invalid connection", "connection is not open",
      "one node produced an error: error writing to connection")){
    original <- simpleError(message, call = quote(unserialize(node$con)))
    failing_backend <- function(...) stop(original)
    condition <- .JAGS_run_backend(failing_backend, list(), parallel = TRUE)
    expect_s3_class(condition, "BayesTools_JAGS_worker_connection_error")
    expect_identical(condition$parent, original)
    expect_match(conditionMessage(condition), message, fixed = TRUE)
    expect_match(conditionMessage(condition), "Automatic retries were stopped.", fixed = TRUE)

    # A local backend may fail on an unrelated file connection. Without a
    # worker topology, preserve its original condition and retry policy.
    expect_identical(.JAGS_run_backend(failing_backend, list(), parallel = FALSE), original)
  }
  original <- simpleError("Node inconsistent with parents")
  expect_identical(.JAGS_run_backend(function(...) stop(original), list(),
    parallel = TRUE), original)
})

test_that("parallel transport failures stop initialization retries immediately", {

  skip_if_not_installed("runjags")
  calls <- stops <- 0L
  original <- simpleError("error writing to connection",
    call = quote(serialize(data, node$con)))
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores){
      structure(as.list(seq_len(cores)), class = "runtime_test_cluster")
    },
    stopCluster = function(cl) stops <<- stops + 1L,
    .package = "parallel"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .bt_attach_fit_contract = function(fit, ...) fit,
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    run.jags = function(...){
      calls <<- calls + 1L
      stop(original)
    },
    .package = "runjags"
  )
  for(silent in c(TRUE, FALSE)){
    calls <- stops <- 0L
    messages <- character()
    fit <- withCallingHandlers(JAGS_fit(
      model_syntax = "model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 2, adapt = 50, burnin = 50, sample = 100,
      autofit_control = list(restarts = 10, sample_extend = 100), parallel = TRUE,
      cores = 2, silent = silent, seed = 1
    ), warning = function(condition){
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    expect_identical(calls, 1L)
    expect_identical(stops, 1L)
    expect_s3_class(fit, "BayesTools_JAGS_worker_connection_error")
    expect_identical(fit$parent, original)
    expect_match(conditionMessage(fit), conditionMessage(original), fixed = TRUE)
    expect_match(conditionMessage(fit), "Automatic retries were stopped.", fixed = TRUE)
    expect_length(messages, if(silent) 0L else 1L)
    expect_false(any(grepl("was restarted", messages, fixed = TRUE)))
    if(!silent){
      expect_match(messages[[1L]], conditionMessage(fit), fixed = TRUE)
    }
  }
})

test_that("cleanup reaches surviving nodes after the first worker disconnects", {

  broken <- rawConnection(raw(), "r+")
  healthy <- rawConnection(raw(), "r+")
  on.exit(try(close(broken), silent = TRUE), add = TRUE)
  on.exit(try(close(healthy), silent = TRUE), add = TRUE)
  cluster <- structure(list(list(id = 1L, con = broken), list(id = 2L, con = healthy)),
    class = "runtime_test_cluster")
  contacted <- integer()
  finished <- FALSE
  testthat::local_mocked_bindings(
    stopCluster = function(cl){
      expect_identical(class(cl), class(cluster))
      if(length(cl) > 1L) stop("first worker disconnected")
      contacted <<- c(contacted, cl[[1L]]$id)
      if(cl[[1L]]$id == 1L) stop("error writing to connection")
      close(cl[[1L]]$con)
    },
    .package = "parallel"
  )
  expect_warning(.JAGS_finish_runtime_setup(function(context){
    finished <<- TRUE
  }, chains = 2L, cl = cluster),
    paste0("Parallel JAGS worker cleanup failed: first worker disconnected. ",
      "The runtime finish callback was not run."), fixed = TRUE)
  expect_identical(contacted, 1:2)
  expect_error(isOpen(broken), "invalid connection", fixed = TRUE)
  expect_error(isOpen(healthy), "invalid connection", fixed = TRUE)
  expect_false(finished)
})

test_that("fit retries and extensions use the actual current worker topology", {

  skip_if_not_installed("runjags")
  basic_fit <- structure(list(end.state = rep("", 4L)),
    class = c("runjags", "BayesTools_fit"))
  attr(basic_fit, "parameter_map") <- .bt_build_parameter_map(character())
  attr(basic_fit, "prior_list") <- list()
  attr(basic_fit, "model_syntax") <- "model{}"
  attr(basic_fit, "required_packages") <- character()
  attr(basic_fit, "jags_modules") <- character()
  attr(basic_fit, "add_parameters") <- character()
  control <- list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL, max_time = NULL, sample_extend = 1,
    restarts = 2, max_extend = 1, check_indicators = FALSE)
  events <- character()
  contexts <- list()
  backend_calls <- convergence_calls <- 0L
  extension_arguments <- fit_arguments <- list()
  cluster_arguments <- list()
  worker_output <- tempfile("runtime-worker-output-", fileext = ".log")
  normalized_output <- file.path(normalizePath(dirname(worker_output), winslash = "/"),
    basename(worker_output))
  setup_failure <- FALSE
  setup <- function(context){

    contexts[[length(contexts) + 1L]] <<- context
    events <<- c(events, paste(context$phase, context$role))
    if(setup_failure && context$role == "worker") stop("worker setup failed")
    invisible(NULL)
  }
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores, ...){
      cluster_arguments[[length(cluster_arguments) + 1L]] <<- list(...)
      structure(as.list(seq_len(cores)), class = "runtime_test_cluster")
    },
    stopCluster = function(cl) events <<- c(events, "stop"),
    clusterApply = function(cl, x, fun, setup) lapply(x, fun, setup = setup),
    .package = "parallel"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) events <<- c(events, "packages"),
    .JAGS_load_modules = function(...) events <<- c(events, "modules"),
    .bt_attach_parameter_map = function(fit, ...) fit,
    .bt_attach_draw_geometry = function(fit, ...) fit,
    .bt_attach_fit_contract = function(fit, ...) fit,
    JAGS_check_convergence = function(...){
      convergence_calls <<- convergence_calls + 1L
      convergence_calls > 1L
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    run.jags = function(...){
      events <<- c(events, "fit")
      fit_arguments[[length(fit_arguments) + 1L]] <<- list(...)
      backend_calls <<- backend_calls + 1L
      if(backend_calls == 1L) stop("retry probe")
      basic_fit
    },
    extend.jags = function(...){
      events <<- c(events, "extend")
      extension_arguments[[length(extension_arguments) + 1L]] <<- list(...)
      basic_fit
    },
    add.summary = function(fit, ...) fit,
    .package = "runjags"
  )

  for(parallel in c(FALSE, TRUE)){
    events <- character()
    contexts <- list()
    backend_calls <- convergence_calls <- 0L
    extension_arguments <- fit_arguments <- list()
    cluster_arguments <- list()
    fit <- JAGS_fit(
      model_syntax = "model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 4, adapt = 50, burnin = 50, sample = 100,
      autofit = TRUE, autofit_control = control, parallel = parallel,
      cores = 3, silent = TRUE, seed = 1, runtime_setup = setup,
      worker_output = worker_output
    )
    start_events <- if(parallel){
      c("start coordinator", rep("start worker", 3L))
    }else "start local"
    finish_events <- if(parallel) c("stop", "finish coordinator") else "finish local"
    expect_identical(events, c("packages", "modules", start_events,
      "fit", "fit", "extend", finish_events))
    expect_identical(attr(fit, "runtime_setup", exact = TRUE), setup)
    expect_null(attr(fit, "worker_output", exact = TRUE))
    expect_false(file.exists(worker_output))
    if(parallel){
      expect_identical(cluster_arguments, list(list(outfile = normalized_output)))
    }else{
      expect_length(cluster_arguments, 0L)
    }
    expect_identical(extension_arguments[[1L]]$method, if(parallel) "rjparallel" else "rjags")
    expect_identical(extension_arguments[[1L]]$cl, fit_arguments[[1L]]$cl)
    expect_identical(extension_arguments[[1L]]$n.sims, if(parallel) 3L else NULL)
    if(parallel){
      expect_identical(vapply(contexts[2:4], `[[`, integer(1), "process_chains"),
        c(2L, 1L, 1L))
    }

    # Changing worker count must override the backend's stored n.sims.
    events <- character()
    contexts <- list()
    extended <- JAGS_extend(fit, autofit_control = control,
      parallel = TRUE, cores = 2, silent = TRUE, worker_output = worker_output)
    expect_identical(events, c("packages", "modules", "start coordinator",
      "start worker", "start worker", "extend", "stop", "finish coordinator"))
    expect_identical(vapply(contexts[2:3], `[[`, integer(1), "process_chains"), c(2L, 2L))
    expect_identical(tail(extension_arguments, 1L)[[1L]]$n.sims, 2L)
    expect_identical(attr(extended, "runtime_setup", exact = TRUE), setup)
    expect_identical(tail(cluster_arguments, 1L)[[1L]], list(outfile = normalized_output))
    expect_null(attr(extended, "worker_output", exact = TRUE))

    events <- character()
    contexts <- list()
    capped <- JAGS_extend(extended, autofit_control = control,
      parallel = TRUE, cores = 8, silent = TRUE)
    # Logging is an explicit setting for this invocation, not replayed from fit.
    expect_identical(tail(cluster_arguments, 1L)[[1L]], list())
    expect_identical(tail(extension_arguments, 1L)[[1L]]$n.sims, 4L)
    expect_identical(vapply(contexts[2:5], `[[`, integer(1), "process_chains"), rep(1L, 4L))

    events <- character()
    local <- JAGS_extend(capped, autofit_control = control,
      parallel = FALSE, silent = TRUE)
    expect_identical(events, c("packages", "modules", "start local", "extend", "finish local"))
    expect_identical(tail(extension_arguments, 1L)[[1L]]$method, "rjags")
    expect_null(tail(extension_arguments, 1L)[[1L]]$n.sims)

    events <- character()
    disabled <- JAGS_extend(local, autofit_control = control,
      parallel = TRUE, cores = 2, silent = TRUE, runtime_setup = NULL)
    expect_identical(events, c("packages", "modules", "extend", "stop"))
    expect_null(attr(disabled, "runtime_setup", exact = TRUE))
  }

  events <- character()
  setup_failure <- TRUE
  expect_error(JAGS_extend(fit, autofit_control = control,
    parallel = TRUE, cores = 2, silent = TRUE), "worker setup failed", fixed = TRUE)
  expect_identical(events, c("packages", "modules", "start coordinator",
    "start worker", "stop", "finish coordinator"))

  # A failed parallel extension retains usable draws, while reporting the
  # worker transport failure instead of suggesting another initialization.
  setup_failure <- FALSE
  extension_calls <- 0L
  testthat::local_mocked_bindings(
    extend.jags = function(...){
      extension_calls <<- extension_calls + 1L
      stop("error reading from connection")
    },
    .package = "runjags"
  )
  expect_warning(retained <- JAGS_extend(fit, autofit_control = control,
    parallel = TRUE, cores = 2, silent = TRUE),
    "returning the last valid fit.*Parallel JAGS worker communication failed.*error reading from connection")
  expect_identical(extension_calls, 1L)
  expect_identical(retained$end.state, fit$end.state)
  expect_s3_class(retained, "runjags")
})
