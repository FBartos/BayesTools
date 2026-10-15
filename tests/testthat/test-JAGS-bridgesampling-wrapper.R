skip_if_not_test_profile("unit")

test_that("bridge callbacks reject undefined priors and retain zero density outside support", {

  posterior <- coda::mcmc(cbind(mu = seq(-.5, .5, length.out = 20),
    auxiliary = c(seq(2, 20, 2), seq(1, 19, 2))))
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

test_that("JAGS_bridgesampling forwards bridge controls and the use_neff default", {

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

  arguments <- list(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = list(
      mu = prior("normal", list(0, 1)),
      sigma = prior("normal", list(0, 1))
    ),
    seed = 1
  )
  result <- do.call(JAGS_bridgesampling, arguments)

  # The target (two standard-normal priors, zero likelihood) is normalized, so
  # its log marginal likelihood is exactly 0: the estimate has to agree with 0
  # within its reported Monte Carlo error.
  expect_s3_class(result, "BayesTools_marglik")
  expect_true(is.finite(result[["logml"]]))
  expect_true(is.finite(result[["repetitions"]][["mcse"]]))
  expect_lt(abs(result[["logml"]]), 4 * result[["repetitions"]][["mcse"]])
  expect_identical(result[["diagnostics"]][["chains"]][["count"]], 3L)

  explicit <- do.call(JAGS_bridgesampling, c(arguments, list(use_neff = FALSE)))
  expect_s3_class(explicit, "BayesTools_marglik")
  expect_true(is.finite(explicit[["logml"]]))
  expect_true(is.finite(explicit[["repetitions"]][["mcse"]]))
  expect_lt(abs(explicit[["logml"]]), 4 * explicit[["repetitions"]][["mcse"]])
  expect_identical(explicit[["diagnostics"]][["chains"]][["count"]], 3L)
})

test_that("JAGS_bridgesampling restores the caller's random-number state on the sampler route", {

  # The sampler is mocked to consume the stream; the seed argument must make
  # its draws reproducible and leave .Random.seed and RNGkind() untouched.
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

  set.seed(77)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  call_bridge(5)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
  call_bridge(5)
  expect_identical(.Random.seed, seed_before)
  expect_identical(seen$draws[[1L]], seen$draws[[2L]])
})

test_that("print.BayesTools_marglik reports a successful sampler result", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  output <- capture.output(
    result <- print(JAGS_bridgesampling(
      fit = coda::as.mcmc(matrix(seq_len(20), ncol = 1, dimnames = list(NULL, "mu"))),
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1))),
      repetitions = 2L
    ))
  )
  expect_identical(output[[1L]], "BayesTools marginal-likelihood result")
  expect_identical(output[[2L]], paste0("  logml (natural log): ", format(result[["logml"]])))
  expect_identical(output[[3L]], "  aggregation: median_finite_logml")
  expect_identical(output[[4L]], "  repetitions: 2 included, 0 failed")
  expect_length(output, 4L)
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

test_that("a bridge target without any declared parameter is refused, not evaluated exactly", {

  constant <- function(parameters, data) -2.5
  raw <- coda::as.mcmc(matrix(0, nrow = 20, ncol = 1, dimnames = list(NULL, "mu")))
  fit <- raw
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(mu = prior("point", list(0)))
  fit <- attach_test_parameter_map(fit)

  # raw draws with the prior list omitted, and an empty prior list
  expect_error(JAGS_bridgesampling(raw, constant), class = "BayesTools_marglik_input")
  expect_error(JAGS_bridgesampling(raw, constant, prior_list = list()), class = "BayesTools_marglik_input")
  expect_error(JAGS_bridgesampling(fit, constant, prior_list = list()), class = "BayesTools_marglik_input")
  declared <- tryCatch(JAGS_bridgesampling(raw, constant), error = function(e) conditionMessage(e))
  expect_match(declared, "without any declared parameter", fixed = TRUE)
  expect_match(declared, "'prior_list'", fixed = TRUE)

  # declared point coordinates keep the exact zero-dimensional route
  for(result in list(
    JAGS_bridgesampling(raw, constant, prior_list = list(mu = prior("point", list(0)))),
    JAGS_bridgesampling(fit, constant, prior_list = list(mu = prior("point", list(0)))),
    JAGS_bridgesampling(fit, constant)
  )){
    expect_identical(result[["aggregation"]][["rule"]], "exact_zero_dimensional")
    expect_equal(result[["logml"]], -2.5, tolerance = 0)
  }
})

test_that("supplied non-formula priors must match the priors stored in the fit", {

  set.seed(16)
  draws <- cbind(mu = abs(stats::rnorm(400)), sigma = stats::rgamma(400, 2, 2))
  fit <- coda::as.mcmc(draws)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(
    mu    = prior("normal", list(0, 1), list(0, Inf)),
    sigma = prior("gamma", list(2, 2))
  )
  fit <- attach_test_parameter_map(fit)
  callback <- function(parameters, data) 0

  # the stored priors, omitted or supplied identically, give one target
  stored <- JAGS_bridgesampling(fit, callback, seed = 1)
  identical_priors <- JAGS_bridgesampling(fit, callback, seed = 1, prior_list = list(
    mu    = prior("normal", list(0, 1), list(0, Inf)),
    sigma = prior("gamma", list(2, 2))
  ))
  expect_identical(identical_priors[["logml"]], stored[["logml"]])
  expect_true(is.finite(stored[["logml"]]))

  # a supplied prior that differs from the stored one targets another model
  truncation <- expect_error(
    JAGS_bridgesampling(fit, callback, seed = 1, prior_list = list(
      mu = prior("normal", list(0, 1)), sigma = prior("gamma", list(2, 2))
    )),
    class = "BayesTools_refit_required"
  )
  expect_match(conditionMessage(truncation), "differ from the fitted priors: mu.", fixed = TRUE)
  parameters <- expect_error(
    JAGS_bridgesampling(fit, callback, seed = 1, prior_list = list(
      mu = prior("normal", list(0, 2), list(0, Inf)), sigma = prior("gamma", list(2, 1))
    )),
    class = "BayesTools_refit_required"
  )
  expect_match(conditionMessage(parameters), "differ from the fitted priors: mu, sigma.", fixed = TRUE)

  # a deliberate subset of the stored priors stays allowed
  subset <- JAGS_bridgesampling(fit, callback, seed = 1, prior_list = list(
    mu = prior("normal", list(0, 1), list(0, Inf))
  ))
  expect_true(is.finite(subset[["logml"]]))
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
    fixed = TRUE,
    class = "BayesTools_refit_required"
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

test_that("JAGS_bridgesampling reports the original error of a failed fit", {

  skip_if_not_installed("runjags")
  testthat::local_mocked_bindings(
    run.jags = function(...) stop("Compilation error on line 1: unknown variable tau"),
    .package = "runjags"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .package = "BayesTools"
  )
  failed <- JAGS_fit(
    model_syntax = "model{ mu ~ dnorm(0, tau) }",
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
      restarts = 1,
      max_extend = 1,
      check_indicators = FALSE
    ),
    silent = TRUE,
    seed = 1
  )
  expect_s3_class(failed, "error")
  expect_s3_class(failed, "BayesTools_fit")

  expect_error(
    JAGS_bridgesampling(
      fit = failed,
      log_posterior = function(parameters, data) 0,
      data = list()
    ),
    paste0(
      "Bridge sampling is unavailable because the model fit failed: ",
      "Compilation error on line 1: unknown variable tau."
    ),
    fixed = TRUE
  )
})

test_that("JAGS_bridgesampling rejects deterministic bridge coordinates", {

  # Review scenario: a row-shaped source tau[i] <- s * tau_data[i] supplied as
  # tau[1..n] bridge coordinates. The draws span only the direction of s,
  # so the bridge target is improper and must not produce an estimate.
  sampler_calls <- 0L
  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      sampler_calls <<- sampler_calls + 1L
      .mock_bridge_sampler(...)
    },
    .package = "bridgesampling"
  )
  set.seed(81)
  n_draws <- 400L
  s <- stats::rgamma(n_draws, 2, 2)
  tau_data <- c(0.6, 1.1, 1.4, 0.9)
  tau_names <- paste0("tau[", seq_along(tau_data), "]")
  deterministic <- cbind(s = s, outer(s, tau_data))
  colnames(deterministic) <- c("s", tau_names)
  call_bridge <- function(posterior){
    JAGS_bridgesampling(
      fit = coda::as.mcmc(posterior),
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(s = prior("gamma", list(2, 2))),
      add_parameters = tau_names,
      add_bounds = list(
        lb = stats::setNames(rep(0, length(tau_names)), tau_names),
        ub = stats::setNames(rep(Inf, length(tau_names)), tau_names)
      )
    )
  }

  expect_error(
    call_bridge(deterministic),
    paste0(
      "The bridge-sampling target was rejected by diagnostics: the posterior ",
      "draws of the bridge coordinates are rank-deficient (rank 1 of 5 varying ",
      "coordinates; linearly dependent coordinate(s): "
    ),
    fixed = TRUE
  )
  expect_identical(sampler_calls, 0L)

  # Free stochastic nodes span the full coordinate space.
  stochastic <- deterministic
  stochastic[, tau_names] <- matrix(
    stats::rgamma(n_draws * length(tau_names), 2, 2),
    ncol = length(tau_names)
  )
  expect_s3_class(call_bridge(stochastic), "BayesTools_marglik")
  expect_identical(sampler_calls, 1L)

  # Fewer draws than coordinates cannot span the coordinate space; the remedy
  # is more draws, not a different source route.
  expect_error(
    call_bridge(stochastic[1:4, , drop = FALSE]),
    paste0(
      "The bridge-sampling target was rejected by diagnostics: 4 posterior ",
      "draws span only rank 3 of 5 varying bridge coordinates. Increase the ",
      "number of posterior draws, for example with a larger 'sample' in ",
      "JAGS_fit() or with JAGS_extend()."
    ),
    fixed = TRUE
  )
  expect_identical(sampler_calls, 1L)
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

test_that("non-finite repetitions name their cause and keep the muffled warnings (F020)", {

  mock_upstream <- function(logml, niter){
    force(logml)
    force(niter)
    function(...){
      warning("upstream diagnostic")
      structure(
        list(
          logml = logml,
          niter = niter,
          mcse_logml = rep(0.1, length(logml)),
          method = "warp3"
        ),
        class = "bridge_list"
      )
    }
  }
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  bridge <- function(...){
    JAGS_bridgesampling(
      fit = posterior,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1))),
      ...
    )
  }
  catch <- function(expr) tryCatch(expr, error = function(e) e)

  # iteration limit exceeded in one of several repetitions
  testthat::local_mocked_bindings(
    bridge_sampler = mock_upstream(c(-10, NA_real_, -12), c(4L, 1001L, 5L)),
    .package = "bridgesampling"
  )
  error <- catch(bridge(maxiter = 1000))
  expect_s3_class(error, "BayesTools_marglik_repetition_failure")
  expect_match(
    conditionMessage(error),
    "repetition(s) 2 (the iteration reached its limit 'maxiter' = 1000)",
    fixed = TRUE
  )
  expect_match(conditionMessage(error), "Upstream warning(s): upstream diagnostic.", fixed = TRUE)
  expect_match(conditionMessage(error), "Increase 'maxiter'", fixed = TRUE)
  expect_match(conditionMessage(error), "set 'nonfinite = \"drop\"' only when", fixed = TRUE)
  expect_identical(error$upstream_warnings, "upstream diagnostic")
  expect_identical(error$repetitions$repetition, 2L)

  # a single non-finite repetition: no iteration limit, and "drop" is no remedy
  testthat::local_mocked_bindings(
    bridge_sampler = mock_upstream(NA_real_, 5L),
    .package = "bridgesampling"
  )
  error <- catch(bridge())
  expect_s3_class(error, "BayesTools_marglik_repetition_failure")
  expect_match(
    conditionMessage(error),
    "repetition(s) 1 (the estimate was non-finite after 5 iteration(s))",
    fixed = TRUE
  )
  expect_match(conditionMessage(error), "Inspect the 'log_posterior' callback", fixed = TRUE)
  expect_false(grepl("drop", conditionMessage(error), fixed = TRUE))
  expect_false(grepl("Increase 'maxiter'", conditionMessage(error), fixed = TRUE))

  # every repetition dropped: the same condition family
  testthat::local_mocked_bindings(
    bridge_sampler = mock_upstream(c(NA_real_, Inf), c(3L, 4L)),
    .package = "bridgesampling"
  )
  error <- catch(bridge(nonfinite = "drop"))
  expect_s3_class(error, "BayesTools_marglik_repetition_failure")
  expect_match(conditionMessage(error), "did not return any finite", fixed = TRUE)
  expect_identical(error$repetitions$repetition, 1:2)
  expect_identical(error$upstream_warnings, "upstream diagnostic")
  expect_false(grepl("drop", sub("^.*finite", "", conditionMessage(error)), fixed = TRUE))
})

test_that("callback condition classes survive the bridge-sampling context (F020)", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  callback_condition <- structure(
    list(message = "callback refused", call = NULL, detail = "kept"),
    class = c("custom_callback_error", "error", "condition")
  )
  error <- tryCatch(
    JAGS_bridgesampling(
      fit = posterior,
      log_posterior = function(parameters, data) stop(callback_condition),
      data = NULL,
      prior_list = list(theta = prior("point", list(location = 0))),
      seed = 1
    ),
    error = function(e) e
  )
  expect_s3_class(error, "custom_callback_error")
  expect_identical(conditionMessage(error), "Bridge sampling failed: callback refused")
  expect_identical(error$detail, "kept")
  expect_identical(error$upstream_warnings, character())
})

test_that("an NA of the callback at a posterior draw is refused (F018)", {

  # No mocking: the counts come from the warnings of the real upstream sampler.
  # The exact posterior of mu | y under mu ~ N(0, 1) and unit-variance data.
  withr::local_seed(9)
  y <- c(0.3, 1.2, -0.4, 0.8, 0.5)
  draws <- stats::rnorm(4000, sum(y) / 6, sqrt(1 / 6))
  fit   <- coda::mcmc(matrix(draws, ncol = 1L, dimnames = list(NULL, "mu")))
  priors <- list(mu = prior("normal", list(0, 1)))
  loglik <- function(parameters, data) sum(stats::dnorm(data$y, parameters$mu, 1, log = TRUE))
  bridge <- function(callback){
    JAGS_bridgesampling(fit, callback, data = list(y = y), prior_list = priors, seed = 1)
  }

  # the valid callback is the reference: no NA, no condition
  expect_no_warning(reference <- bridge(loglik))
  expect_true(reference$success)

  # NA above mu = 0.7 hits the posterior draws the sampler evaluates (the
  # second half of the supplied draws)
  evaluated <- draws[2001:4000]
  error <- tryCatch(
    bridge(function(parameters, data) if(parameters$mu > 0.7) NA_real_ else loglik(parameters, data)),
    error = function(e) e
  )
  expect_s3_class(error, "error")
  expect_s3_class(error, "BayesTools_marglik_callback_na")
  expect_s3_class(error, "BayesTools_marglik_input")
  expect_identical(error$n_posterior_na, as.numeric(sum(evaluated > 0.7)))
  expect_match(
    conditionMessage(error),
    paste0("returned NA at ", sum(evaluated > 0.7), " of 2000 posterior draws"),
    fixed = TRUE
  )
  expect_match(conditionMessage(error), "-Inf for a state of zero density", fixed = TRUE)
  expect_gt(length(error$upstream_warnings), 0L)
})

test_that("an NA of the callback at a proposal draw is evaluated as zero density and reported (F018)", {

  # NA only beyond the range of the posterior draws is met by proposal draws
  # alone (uniform draws have lighter tails than the normal proposal): the
  # estimate is finite, the NA is evaluated as zero density and its share is
  # reported
  withr::local_seed(9)
  bounded <- coda::mcmc(matrix(stats::runif(4000, -1, 1), ncol = 1L, dimnames = list(NULL, "mu")))
  tails   <- function(parameters, data) if(abs(parameters$mu) > 1.5) NA_real_ else 0
  warned  <- NULL
  result  <- withCallingHandlers(
    JAGS_bridgesampling(bounded, tails, data = list(), prior_list = list(mu = prior("normal", list(0, 1))),
                        seed = 1, repetitions = 2L),
    warning = function(w){
      warned <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(warned, "BayesTools_marglik_callback_na")
  expect_s3_class(warned, "BayesTools_marglik_input")
  expect_s3_class(warned, "warning")
  expect_gt(warned$n_proposal_na, 0)
  expect_identical(warned$n_proposal, 4000)
  expect_match(
    conditionMessage(warned),
    paste0("returned NA at ", warned$n_proposal_na, " of 4000 bridge-sampling proposal draws (",
           format(100 * warned$n_proposal_na / 4000, digits = 2), "%)"),
    fixed = TRUE
  )
  expect_true(result$success)
  expect_true(is.finite(result$logml))
  expect_true(any(grepl("on the proposal draws produced NA", result$diagnostics$upstream_warnings, fixed = TRUE)))
})

test_that("the NA counts of the callback are summed over the repetitions of the sampler (F018)", {

  signal <- BayesTools:::.bt_JAGS_bridge_signal_callback_na
  proposal <- "5 evaluation(s) of log_prob() on the proposal draws produced NA and have been replaced by -Inf."
  other    <- c("3 of the 2000 log_prob() evaluations on the proposal draws produced -Inf/Inf.", "upstream diagnostic")

  expect_true(signal(other, 2000, 3))
  expect_warning(
    signal(c(other, proposal, proposal, proposal), n_posterior = 2000, repetitions = 3),
    "returned NA at 15 of 6000 bridge-sampling proposal draws (0.25%)",
    fixed = TRUE,
    class = "BayesTools_marglik_callback_na"
  )
  expect_error(
    signal(c(proposal, "12 evaluation(s) of log_prob() on the posterior draws produced NA and have been replaced by -Inf."), 2000, 1),
    "returned NA at 12 of 2000 posterior draws",
    fixed = TRUE,
    class = "BayesTools_marglik_input"
  )
})

test_that("an NA of the callback with method = 'warp3' is outside the NA policy (F018)", {

  # The roxygen states the scope of the policy: bridgesampling stops on an NA
  # of the warp3 route with an upstream error, before it words any warning.
  withr::local_seed(9)
  y <- c(0.3, 1.2, -0.4, 0.8, 0.5)
  draws <- stats::rnorm(4000, sum(y) / 6, sqrt(1 / 6))
  fit   <- coda::mcmc(matrix(draws, ncol = 1L, dimnames = list(NULL, "mu")))
  priors <- list(mu = prior("normal", list(0, 1)))
  loglik <- function(parameters, data) sum(stats::dnorm(data$y, parameters$mu, 1, log = TRUE))
  bridge <- function(callback, method){
    JAGS_bridgesampling(fit, callback, data = list(y = y), prior_list = priors, seed = 1, method = method)
  }
  with_na <- function(parameters, data) if(parameters$mu > 0.7) NA_real_ else loglik(parameters, data)

  expect_no_warning(valid <- bridge(loglik, "warp3"))
  expect_true(valid$success)
  expect_s3_class(tryCatch(bridge(with_na, "normal"), error = function(e) e), "BayesTools_marglik_callback_na")

  error <- tryCatch(bridge(with_na, "warp3"), error = function(e) e)
  expect_s3_class(error, "error")
  expect_false(inherits(error, "BayesTools_marglik_callback_na"))
  expect_match(conditionMessage(error), "Bridge sampling failed: ", fixed = TRUE)
  expect_match(conditionMessage(error), "NAs are not allowed in subscripted assignments", fixed = TRUE)

  # the warnings that bridgesampling words for the warp3 route are not counted
  signal <- BayesTools:::.bt_JAGS_bridge_signal_callback_na
  expect_true(signal(
    "12 evaluation(s) of log_prob() on the warp-transformed posterior draws produced NA and have been replaced by -Inf.",
    2000, 1
  ))
})

test_that("the inclusion-gate refusal of bridge sampling lists at most three indicators", {

  indicator_design <- function(...){
    structure(
      list(
        parameter = "mu",
        random_allocations = list(list(
          inclusion = lapply(c(...), function(name) list(indicator_name = name))
        ))
      ),
      class = "BayesTools_formula_design"
    )
  }
  expect_true(BayesTools:::.bt_JAGS_bridge_check_no_allocation_inclusion(list(mu = indicator_design())))
  expect_error(
    BayesTools:::.bt_JAGS_bridge_check_no_allocation_inclusion(list(mu = indicator_design("inclusion_g"))),
    "Affected indicator node(s): 'inclusion_g'.",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_check_no_allocation_inclusion(
      list(mu = indicator_design("g1", "g2", "g3", "g4"))
    ),
    "Affected indicator node(s): 'g1', 'g2', 'g3', ....",
    fixed = TRUE
  )
})

test_that("empty formula inputs are no formula input for a fit without formulas (F074)", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  bridge <- function(...){
    JAGS_bridgesampling(
      fit = posterior,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(mu = prior("normal", list(0, 1))),
      seed = 1,
      ...
    )
  }
  reference <- bridge()
  for(argument in c("formula_list", "formula_data_list", "formula_prior_list",
                    "formula_scale_list", "formula_random_prior_list",
                    "formula_random_effects_compile_list")){
    result <- do.call(bridge, structure(list(list()), names = argument))
    expect_identical(result[["logml"]], reference[["logml"]])
  }

  # a supplied formula input still cannot replace missing fitted metadata
  error <- tryCatch(
    bridge(formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1))))),
    error = function(e) e
  )
  expect_s3_class(error, "BayesTools_refit_required")
  expect_match(conditionMessage(error), "fitted formula-design metadata are missing", fixed = TRUE)
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

test_that("P-F8 bridge wrapper preserves explicit use_neff and its default", {

  seen <- list()
  testthat::local_mocked_bindings(
    bridge_sampler = function(use_neff, ...){
      arguments <- c(list(use_neff = use_neff), list(...))
      seen[[length(seen) + 1L]] <<- arguments[names(arguments) == "use_neff"]
      structure(list(logml = 0, niter = 7L, mcse_logml = .05,
                     method = "normal"), class = "bridge")
    },
    .package = "bridgesampling"
  )
  posterior <- coda::mcmc(matrix(seq_len(20), ncol = 1L,
                                dimnames = list(NULL, "mu")))
  arguments <- list(fit = posterior, log_posterior = function(parameters, data) 0,
                    data = list(), prior_list = list(mu = prior("normal", list(0, 1))))
  default <- do.call(JAGS_bridgesampling, arguments)
  expect_identical(seen, list(list(use_neff = TRUE)))
  explicit <- do.call(JAGS_bridgesampling, c(arguments, list(use_neff = FALSE)))
  expect_identical(seen, list(list(use_neff = TRUE), list(use_neff = FALSE)))
  expect_s3_class(default, "BayesTools_marglik")
  expect_s3_class(explicit, "BayesTools_marglik")
})

test_that("P-F8 deferred callback arguments record warnings without escape", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::mcmc(matrix(seq_len(20), ncol = 1L,
                                dimnames = list(NULL, "mu")))
  expect_no_warning(result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data, extra) extra,
    data = list(),
    prior_list = list(mu = prior("normal", list(0, 1))),
    extra = { warning("extra warning"); 0 }
  ))
  expect_identical(result$diagnostics$upstream_warnings, "extra warning")
  expect_true(result$repetitions$success)
})

test_that("P-F8 deferred callback argument errors retain bridge context", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::mcmc(matrix(seq_len(20), ncol = 1L,
                                dimnames = list(NULL, "mu")))
  expect_error(JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data, extra) extra,
    data = list(),
    prior_list = list(mu = prior("normal", list(0, 1))),
    extra = stop("extra unavailable")
  ), "Bridge sampling failed: extra unavailable", fixed = TRUE)
})

test_that("P-F8 exact point route ignores upstream use_neff controls", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  arguments <- list(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = NULL,
    prior_list = list(theta = prior("point", list(location = 0))),
    seed = 1
  )
  default <- do.call(JAGS_bridgesampling, arguments)
  expect_identical(default$logml, 0)
  expect_identical(default$aggregation$rule, "exact_zero_dimensional")
  for(use_neff in list(TRUE, FALSE, NULL)){
    explicit <- do.call(JAGS_bridgesampling, c(arguments, list(use_neff = use_neff)))
    expect_identical(explicit$logml, default$logml)
    expect_identical(explicit$aggregation, default$aggregation)
  }
})

test_that("P-F8 exact point route forwards genuine likelihood arguments", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data, extra) extra,
    data = NULL,
    prior_list = list(theta = prior("point", list(location = 0))),
    use_neff = FALSE,
    packages = "stats",
    extra = -3,
    seed = 1
  )
  expect_identical(result$logml, -3)
  expect_identical(result$aggregation$rule, "exact_zero_dimensional")
})

test_that("P-F8 exact point route ignores all named upstream controls", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  arguments <- list(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = NULL,
    prior_list = list(theta = prior("point", list(location = 0))),
    seed = 1
  )
  controls <- list(packages = character(), use_neff = FALSE, varlist = character(),
                   envir = .GlobalEnv, rcppFile = NULL, param_types = character(),
                   verbose = FALSE)
  for(control in c("varlist", "envir", "rcppFile", "param_types", "verbose")){
    explicit <- do.call(JAGS_bridgesampling, c(arguments, controls[control]))
    expect_identical(explicit$logml, 0)
  }
  combined <- do.call(JAGS_bridgesampling, c(arguments, controls))
  expect_identical(combined$logml, 0)
  expect_identical(combined$aggregation$rule, "exact_zero_dimensional")

  arguments$log_posterior <- function(parameters, data, foo) foo
  genuine <- do.call(JAGS_bridgesampling, c(arguments, controls, list(foo = -3)))
  expect_identical(genuine$logml, -3)
})

test_that("P-F8 exact point callback warnings are recorded without escape", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  call_bridge <- function(log_posterior = function(parameters, data, foo) foo, ...){
    JAGS_bridgesampling(
      fit = posterior, log_posterior = log_posterior, data = NULL,
      prior_list = list(theta = prior("point", list(location = 0))),
      seed = 1, ...
    )
  }
  default <- call_bridge(foo = -3)
  set.seed(55)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  expect_no_warning(deferred <- call_bridge(
    foo = { stats::runif(1); warning("55 deferred exact warning"); -3 }
  ))
  expect_identical(deferred$diagnostics$upstream_warnings, "55 deferred exact warning")
  expect_no_warning(body <- call_bridge(function(parameters, data, foo){
    stats::runif(1)
    warning("55 exact body warning")
    warning("55 exact body warning")
    foo
  }, foo = -3))
  expect_identical(body$diagnostics$upstream_warnings, "55 exact body warning")
  for(result in list(deferred, body)){
    expect_identical(result$logml, -3)
    expect_identical(result$aggregation, default$aggregation)
    expect_identical(result$repetitions, default$repetitions)
    expect_identical(result$diagnostics$chains, default$diagnostics$chains)
  }
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
})

test_that("P-F8 exact point callback errors retain bridge context and RNG state", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  call_bridge <- function(log_posterior = function(parameters, data) 0, ...){
    JAGS_bridgesampling(
      fit = posterior, log_posterior = log_posterior, data = NULL,
      prior_list = list(theta = prior("point", list(location = 0))),
      seed = 1, ...
    )
  }
  set.seed(55)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  expect_error(call_bridge(function(parameters, data, foo) foo,
    foo = { stats::runif(1); stop("55 deferred exact error") }),
    "Bridge sampling failed: 55 deferred exact error", fixed = TRUE)
  expect_error(call_bridge(function(parameters, data){
    stats::runif(1)
    stop("55 exact body error")
  }), "Bridge sampling failed: 55 exact body error", fixed = TRUE)
  expect_error(call_bridge(varlist = stop("55 reserved control evaluation")),
    "Bridge sampling failed: 55 reserved control evaluation", fixed = TRUE)
  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)
})

test_that("P-F8 exact point result validation retains its own error context", {

  posterior <- coda::mcmc(matrix(0, nrow = 100L, ncol = 1L,
                                dimnames = list(NULL, "theta")))
  call_bridge <- function(value){
    JAGS_bridgesampling(
      fit = posterior, log_posterior = function(parameters, data) value,
      data = NULL, prior_list = list(theta = prior("point", list(location = 0)))
    )
  }
  na_error <- tryCatch(call_bridge(NA_real_), error = identity)
  expect_identical(conditionMessage(na_error), paste0(
    "The exact zero-dimensional log marginal likelihood evaluated to NA. ",
    "Check that the log-posterior callback returns a numeric value."
  ))
  inf_error <- tryCatch(call_bridge(Inf), error = identity)
  expect_identical(conditionMessage(inf_error), paste0(
    "'logml' must be one numeric natural-log marginal likelihood (NA for a ",
    "failed computation) and may only be infinite when it is -Inf."
  ))
})
