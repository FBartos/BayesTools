skip_if_not_test_profile("unit")

.mock_bridge_sampler <- function(...){

  arguments <- list(...)
  repetitions <- arguments[["repetitions"]]
  if(is.null(repetitions)){
    repetitions <- 1L
  }
  maxiter <- arguments[["maxiter"]]
  niter <- if(maxiter <= 1L){
    c(2L, rep.int(1L, repetitions - 1L))
  }else{
    rep.int(7L, repetitions)
  }

  structure(
    list(
      logml = seq_len(repetitions),
      niter = niter
    ),
    class = "bridge_list"
  )
}

test_that("JAGS_bridgesampling preserves repeated bridge estimates", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  prior_list <- list(mu = prior("normal", list(0, 1)))

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = prior_list,
    repetitions = 2,
    maxiter = 1000
  )

  expect_s3_class(result, "bridge_list")
  expect_length(result[["logml"]], 2L)
  expect_length(result[["niter"]], 2L)
  expect_true(all(is.finite(result[["logml"]])))
  expect_true(all(result[["niter"]] <= 1000L))
  expect_null(attr(result, "warning"))
})

test_that("JAGS_bridgesampling checks repeated iteration limits collectively", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  prior_list <- list(mu = prior("normal", list(0, 1)))

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = prior_list,
    repetitions = 2,
    maxiter = 1
  )

  expect_s3_class(result, "bridge_list")
  expect_length(result[["niter"]], 2L)
  expect_true(any(result[["niter"]] > 1L))
  expect_identical(
    attr(result, "warning"),
    paste(
      "Marginal likelihood could not be estimated within the maximum number",
      "of iterations and might be more variable than usual."
    )
  )
})

test_that("JAGS_bridgesampling validates the log-posterior callback eagerly", {

  expect_error(
    JAGS_bridgesampling(
      fit = NULL,
      log_posterior = 1
    ),
    "'log_posterior' must be a function.",
    fixed = TRUE
  )
})

test_that("global log-posterior callbacks survive a PSOCK round trip", {

  cluster <- tryCatch(
    parallel::makePSOCKcluster(
      1L,
      setup_timeout = 10,
      timeout = 30
    ),
    error = function(e) e
  )
  if(inherits(cluster, "error")){
    skip(paste(
      "A PSOCK worker could not be started:",
      conditionMessage(cluster)
    ))
  }
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      arguments <- list(...)
      worker_value <- parallel::clusterCall(
        cluster,
        function(compiled_callback){
          user_callback <- get(
            "log_posterior",
            envir = environment(compiled_callback),
            inherits = FALSE
          )
          user_callback(
            parameters = list(mu = 1),
            data = list(offset = 2)
          )
        },
        arguments[["log_posterior"]]
      )
      structure(
        list(logml = unlist(worker_value), niter = 1L),
        class = "bridge"
      )
    },
    .package = "bridgesampling"
  )

  global_names <- paste0(
    ".BayesTools_psock_",
    Sys.getpid(),
    c("_callback", "_posterior", "_priors")
  )
  callback_name <- global_names[[1L]]
  posterior_name <- global_names[[2L]]
  prior_name <- global_names[[3L]]
  assign(
    callback_name,
    function(parameters, data){
      parameters[["mu"]] + data[["offset"]]
    },
    envir = .GlobalEnv
  )
  assign(
    posterior_name,
    coda::as.mcmc(matrix(
      seq_len(20),
      ncol = 1,
      dimnames = list(NULL, "mu")
    )),
    envir = .GlobalEnv
  )
  assign(
    prior_name,
    list(mu = prior("normal", list(0, 1))),
    envir = .GlobalEnv
  )
  on.exit(
    rm(list = global_names, envir = .GlobalEnv),
    add = TRUE
  )

  bridge_call <- substitute(
    BayesTools::JAGS_bridgesampling(
      fit = FIT,
      log_posterior = CALLBACK,
      data = list(offset = 2),
      prior_list = PRIORS,
      maxiter = 1000
    ),
    list(
      FIT = as.name(posterior_name),
      CALLBACK = as.name(callback_name),
      PRIORS = as.name(prior_name)
    )
  )
  result <- eval(bridge_call, envir = .GlobalEnv)

  expect_false(inherits(result, "error"))
  expect_equal(result[["logml"]], 3)
  expect_identical(result[["niter"]], 1L)
})
