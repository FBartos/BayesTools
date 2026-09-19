skip_if_not_test_profile("unit")

test_that("bridge callbacks reject undefined priors and retain zero density outside support", {

  posterior <- coda::mcmc(cbind(mu = seq(-.5, .5, length.out = 20), auxiliary = 1:20))
  proposal <- NA_real_
  callback_value <- NULL
  random_prior <- NULL
  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      arguments <- list(...)
      samples <- arguments$samples[1L, ]
      samples[["mu"]] <- proposal
      if(!is.null(random_prior)){
        arguments$bridge_formula_random_prior_evaluator <- random_prior
      }
      callback_names <- c(
        "data", "bridge_prior_evaluator", "bridge_formula_prior_evaluator",
        "bridge_formula_random_prior_evaluator", "bridge_formula_parameter_evaluator",
        "add_parameters", "fixed_random_latent", "bridge_context", "bridge_context_evaluator"
      )
      callback_value <<- do.call(
        arguments$log_posterior,
        c(list(samples.row = samples), arguments[callback_names])
      )
      structure(list(logml = 0, niter = 1L, mcse_logml = 0, method = "normal"), class = "bridge")
    },
    .package = "bridgesampling"
  )
  for(with_auxiliary in c(FALSE, TRUE)){
    arguments <- list(
      fit = posterior,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1), list(-1, 1)))
    )
    if(with_auxiliary){
      arguments$add_parameters <- "auxiliary"
      arguments$add_bounds <- list(lb = c(auxiliary = -Inf), ub = c(auxiliary = Inf))
    }
    for(invalid in c(NA_real_, NaN)){
      proposal <- invalid
      expect_error(
        do.call(JAGS_bridgesampling, arguments),
        "Bridge log prior evaluated to NA or NaN. Check the prior specification and monitored posterior samples.",
        fixed = TRUE
      )
    }
    proposal <- 2
    do.call(JAGS_bridgesampling, arguments)
    expect_identical(callback_value, -Inf)
  }

  formula_output <- JAGS_formula(
    ~ 1 + diag(1 | id), "mu", data.frame(id = factor(c("a", "b"))),
    list(intercept = prior("point", list(0))),
    prior_random = prior_random(id = random_block(sd = prior("point", list(1))))
  )
  testthat::local_mocked_bindings(
    .bt_JAGS_bridge_compile_random_effect_lkj_prior = function(...) function(samples) NaN,
    .package = "BayesTools"
  )
  compiled_random <- BayesTools:::.bt_JAGS_bridge_compile_formula_random_prior_evaluator(
    list(mu = formula_output$formula_design)
  )
  random_term <- formula_output$formula_design$random_effects[[1L]]
  latent_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term, random_term$n_groups, random_term$n_columns
  ))
  latent_samples <- stats::setNames(rep(0, length(latent_names)), latent_names)
  random_prior <- list(
    log_prior = function(samples) compiled_random$log_prior(latent_samples),
    uses_posterior_row = FALSE
  )
  proposal <- .2
  expect_error(
    do.call(JAGS_bridgesampling, arguments),
    "Bridge log prior evaluated to NA or NaN. Check the prior specification and monitored posterior samples.",
    fixed = TRUE
  )
})

test_that("bridge chain metadata preserves runjags conversion errors", {

  testthat::local_mocked_bindings(
    as.mcmc.list = function(...) stop("Malformed runjags chains.", call. = FALSE),
    .package = "coda"
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_chain_metadata(
      structure(list(), class = "runjags"), matrix(1, nrow = 2)
    ),
    "Malformed runjags chains.", fixed = TRUE
  )
})

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
      niter = niter,
      mcse_logml = rep.int(0.05, repetitions),
      method = "normal"
    ),
    class = "bridge_list"
  )
}

test_that("JAGS_bridgesampling aggregates repeated bridge estimates explicitly", {

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

  expect_s3_class(result, "BayesTools_marglik")
  expect_identical(result[["scale"]], "natural_log")
  expect_equal(result[["logml"]], 1.5)
  expect_identical(result[["aggregation"]][["rule"]], "median_finite_logml")
  expect_identical(result[["aggregation"]][["n_included"]], 2L)
  expect_identical(result[["aggregation"]][["n_failed"]], 0L)
  expect_identical(result[["repetitions"]][["logml"]], c(1, 2))
  expect_identical(result[["repetitions"]][["niter"]], c(7, 7))
  expect_identical(result[["repetitions"]][["mcse"]], c(0.05, 0.05))
  expect_true(all(result[["repetitions"]][["success"]]))
  expect_true(all(result[["repetitions"]][["within_maxiter"]]))
  expect_identical(result[["repetitions"]][["method"]], c("normal", "normal"))
  expect_identical(result[["repetitions"]][["n_chains"]], c(1L, 1L))
  expect_identical(result[["repetitions"]][["n_draws"]], c(20L, 20L))
  expect_s3_class(result[["diagnostics"]][["upstream"]], "bridge_list")
})

test_that("JAGS_bridgesampling forwards bridge controls and fitted-chain neff", {

  seen <- new.env(parent = emptyenv())
  bridge_sampler <- function(...){
    arguments <- list(...)
    seen$controls <- arguments[c(
      "cores", "repetitions", "method", "maxiter", "silent", "use_neff"
    )]
    .mock_bridge_sampler(...)
  }
  testthat::local_mocked_bindings(
    bridge_sampler = bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))

  JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = list(mu = prior("normal", list(0, 1))),
    cores       = 3L,
    repetitions = 4L,
    method      = "warp3",
    maxiter     = 2500L,
    silent      = FALSE
  )

  expect_identical(
    seen$controls,
    list(
      cores       = 3L,
      repetitions = 4L,
      method      = "warp3",
      maxiter     = 2500L,
      silent      = FALSE,
      use_neff    = TRUE
    )
  )
})

test_that("JAGS_bridgesampling passes use_neff as the logical flag upstream defines", {

  # bridgesampling::bridge_sampler() branches on `if (use_neff)` and computes
  # the ESS itself. Anything but a logical scalar is either silently ignored
  # (length 1) or aborts the sampler with "the condition has length > 1".
  seen <- new.env(parent = emptyenv())
  bridge_sampler <- function(...){
    arguments <- list(...)
    seen$use_neff <- arguments[["use_neff"]]
    .mock_bridge_sampler(...)
  }
  testthat::local_mocked_bindings(
    bridge_sampler = bridge_sampler,
    .package = "bridgesampling"
  )

  chain1 <- coda::mcmc(matrix(
    c(rnorm(20, 0, 1), rnorm(20, 0.1, 1)),
    ncol = 2,
    dimnames = list(NULL, c("mu", "sigma"))
  ))
  chain2 <- coda::mcmc(matrix(
    c(rnorm(20, 0, 1), rnorm(20, 0.1, 1)),
    ncol = 2,
    dimnames = list(NULL, c("mu", "sigma"))
  ))
  posterior <- coda::mcmc.list(chain1, chain2)

  JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = list(
      mu = prior("normal", list(0, 1)),
      sigma = prior("normal", list(0, 1))
    )
  )

  expect_true(is.logical(seen$use_neff))
  expect_length(seen$use_neff, 1L)
  expect_false(is.na(seen$use_neff))
})


test_that("JAGS_bridgesampling drives the real upstream sampler on multiple chains", {

  # Deliberately NOT mocked: every argument this wrapper forwards has to be
  # one bridgesampling::bridge_sampler() actually accepts. A mocked sampler
  # cannot catch a forwarded value that upstream rejects.
  set.seed(1)
  make_chain <- function() coda::mcmc(matrix(
    c(rnorm(400), rnorm(400)),
    ncol = 2,
    dimnames = list(NULL, c("mu", "sigma"))
  ))
  posterior <- coda::mcmc.list(make_chain(), make_chain(), make_chain())

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = list(
      mu = prior("normal", list(0, 1)),
      sigma = prior("normal", list(0, 1))
    ),
    seed = 1
  )

  expect_s3_class(result, "BayesTools_marglik")
  expect_true(is.finite(result[["logml"]]))
  expect_identical(result[["diagnostics"]][["chains"]][["count"]], 3L)
})

test_that("JAGS_bridgesampling seeds bridge proposal draws explicitly", {

  seen <- new.env(parent = emptyenv())
  seen$draws <- numeric()
  bridge_sampler <- function(...){
    seen$draws <- c(seen$draws, stats::runif(1))
    .mock_bridge_sampler(...)
  }
  testthat::local_mocked_bindings(
    bridge_sampler = bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  call_bridge <- function(seed){
    JAGS_bridgesampling(
      fit = posterior,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1))),
      seed = seed
    )
  }

  call_bridge(11)
  invisible(stats::runif(10))
  call_bridge(11)

  expect_identical(seen$draws[[1L]], seen$draws[[2L]])
  expect_error(
    call_bridge(1.5),
    "The 'seed' argument must be an integer vector.",
    fixed = TRUE
  )
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

  expect_s3_class(result, "BayesTools_marglik")
  expect_true(any(!result[["repetitions"]][["within_maxiter"]]))
  expect_false(result[["repetitions"]][["success"]][[1L]])
  expect_identical(
    result[["repetitions"]][["warning"]][[1L]],
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

test_that("JAGS_bridgesampling bypasses formula and context replay for ordinary models", {

  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      arguments <- list(...)
      callback <- arguments[["log_posterior"]]
      required <- setdiff(names(formals(callback)), c("samples.row", "..."))
      callback_arguments <- arguments[required]
      callback_arguments[["samples.row"]] <- c(mu = 0.25, tau = 0.50)
      value <- do.call(callback, callback_arguments)
      structure(
        list(
          logml = value,
          niter = 1L,
          mcse_logml = 0.01,
          method = "normal"
        ),
        class = "bridge"
      )
    },
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    rep(c(0.25, 0.50), each = 20L),
    ncol = 2L,
    dimnames = list(NULL, c("mu", "tau"))
  ))
  priors <- list(
    mu  = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1), truncation = list(0, Inf))
  )

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data){
      seen$parameters <- parameters
      data$offset + parameters$mu - parameters$tau
    },
    data = list(offset = -2),
    prior_list = priors,
    bridge_context = FALSE
  )

  expect_identical(seen$parameters, list(mu = 0.25, tau = 0.50))
  expected <- lpdf(priors$mu, 0.25) + lpdf(priors$tau, 0.50) - 2.25
  expect_equal(result$logml, expected, tolerance = 0)
})

test_that("JAGS_bridgesampling evaluates fixed scalar and vector models exactly", {

  posterior <- coda::as.mcmc(matrix(
    rep(c(0.25, 2, 2), each = 20),
    nrow = 20,
    dimnames = list(NULL, c("mu", "beta[1]", "beta[2]"))
  ))
  prior_list <- list(
    mu   = prior("point", list(0.25)),
    beta = prior("mpoint", list(2, 2))
  )

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data, bridge_context){
      expect_s3_class(
        bridge_context,
        "BayesTools_bridge_nodes_context"
      )
      expect_named(bridge_context, "nodes")
      expect_equal(bridge_context$nodes[["mu"]], parameters[["mu"]])
      expect_equal(
        unname(bridge_context$nodes[c("beta[1]", "beta[2]")]),
        parameters[["beta"]]
      )
      data[["offset"]] + parameters[["mu"]] + sum(parameters[["beta"]])
    },
    data = list(offset = -10),
    prior_list = prior_list,
    bridge_context = "nodes",
    packages = "BayesTools"
  )

  expect_s3_class(result, "BayesTools_marglik")
  expect_equal(result[["logml"]], -5.75)
  expect_identical(
    result[["aggregation"]][["rule"]],
    "exact_zero_dimensional"
  )
  expect_identical(result[["aggregation"]][["n_repetitions"]], 0L)
  expect_null(result[["diagnostics"]][["upstream"]])
  expect_identical(result[["diagnostics"]][["chains"]][["count"]], 1L)
  expect_identical(result[["diagnostics"]][["chains"]][["draws"]], 20L)
})

test_that("BayesTools fits require the parameter map for bridge replay", {

  fit <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  class(fit) <- c("BayesTools_fit", class(fit))

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      prior_list = list(mu = prior("normal", list(0, 1)))
    ),
    "Refit the model with the current BayesTools version.",
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
        list(
          logml = unlist(worker_value),
          niter = 1L,
          mcse_logml = 0.01,
          method = "normal"
        ),
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
  expect_s3_class(result, "BayesTools_marglik")
  expect_equal(result[["logml"]], 3)
  expect_identical(result[["repetitions"]][["niter"]], 1)
  expect_identical(result[["diagnostics"]][["chains"]][["count"]], 1L)
})

test_that("JAGS_bridgesampling aborts on non-finite repetitions by default", {

  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      structure(
        list(
          logml = c(-10, NA_real_, -12),
          niter = c(4L, 1000L, 5L),
          mcse_logml = c(0.1, NA_real_, 0.2),
          method = "warp3"
        ),
        class = "bridge_list"
      )
    },
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))

  expect_error(
    JAGS_bridgesampling(
      fit = posterior,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1)))
    ),
    class = "BayesTools_marglik_repetition_failure"
  )
})

test_that("explicit drop policy records non-finite repetitions", {

  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      warning("upstream diagnostic")
      structure(
        list(
          logml = c(-10, Inf, -14),
          niter = c(4L, 1000L, 5L),
          mcse_logml = c(0.1, NA_real_, 0.2),
          method = "warp3"
        ),
        class = "bridge_list"
      )
    },
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))

  expect_warning(
    result <- JAGS_bridgesampling(
      fit = posterior,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1))),
      nonfinite = "drop"
    ),
    "Dropped non-finite bridge-sampling repetition\\(s\\) 2"
  )

  expect_equal(result[["logml"]], -12)
  expect_identical(result[["aggregation"]][["n_included"]], 2L)
  expect_identical(result[["aggregation"]][["n_failed"]], 1L)
  expect_false(result[["repetitions"]][["finite"]][[2L]])
  expect_match(
    result[["repetitions"]][["error"]][[2L]],
    "non-finite log marginal likelihood",
    fixed = TRUE
  )
  expect_identical(
    result[["diagnostics"]][["upstream_warnings"]],
    "upstream diagnostic"
  )
})

test_that("manual marginal-likelihood objects declare their scale", {

  result <- bridgesampling_object(-Inf)

  expect_s3_class(result, "BayesTools_marglik")
  expect_identical(result[["scale"]], "natural_log")
  expect_identical(result[["aggregation"]][["rule"]], "supplied_scalar")
  expect_equal(nrow(result[["repetitions"]]), 0L)
  expect_error(bridgesampling_object(Inf), "may only be infinite when it is -Inf")
  expect_error(bridgesampling_object(c(1, 2)), "one numeric")
})
