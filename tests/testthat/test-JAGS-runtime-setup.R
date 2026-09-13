skip_if_not_test_profile("unit")

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
  setup_failure <- FALSE
  setup <- function(context){

    contexts[[length(contexts) + 1L]] <<- context
    events <<- c(events, paste(context$phase, context$role))
    if(setup_failure && context$role == "worker") stop("worker setup failed")
    invisible(NULL)
  }
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores){
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
    fit <- JAGS_fit(
      model_syntax = "model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
      prior_list = list(mu = prior("normal", list(0, 1))),
      chains = 4, adapt = 50, burnin = 50, sample = 100,
      autofit = TRUE, autofit_control = control, parallel = parallel,
      cores = 3, silent = TRUE, seed = 1, runtime_setup = setup
    )
    start_events <- if(parallel){
      c("start coordinator", rep("start worker", 3L))
    }else "start local"
    finish_events <- if(parallel) c("stop", "finish coordinator") else "finish local"
    expect_identical(events, c("packages", "modules", start_events,
      "fit", "fit", "extend", finish_events))
    expect_identical(attr(fit, "runtime_setup", exact = TRUE), setup)
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
      parallel = TRUE, cores = 2, silent = TRUE)
    expect_identical(events, c("packages", "modules", "start coordinator",
      "start worker", "start worker", "extend", "stop", "finish coordinator"))
    expect_identical(vapply(contexts[2:3], `[[`, integer(1), "process_chains"), c(2L, 2L))
    expect_identical(tail(extension_arguments, 1L)[[1L]]$n.sims, 2L)
    expect_identical(attr(extended, "runtime_setup", exact = TRUE), setup)

    events <- character()
    contexts <- list()
    capped <- JAGS_extend(extended, autofit_control = control,
      parallel = TRUE, cores = 8, silent = TRUE)
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
})
