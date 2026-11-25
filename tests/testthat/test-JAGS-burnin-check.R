skip_if_not_test_profile("unit")

# The readiness check between burn-in and sampling: its statistic, its
# selection, its extension loop and its integration in JAGS_fit(), with the
# backend mocked.

.burnin_rhat <- function(chains){

  values <- unlist(chains, use.names = FALSE)
  chain_id <- rep(seq_along(chains), lengths(chains))
  BayesTools:::.bt_burnin_split_Rhat(values, chain_id, length(chains))
}

.burnin_unsplit_point <- function(chains){

  coda::gelman.diag(coda::mcmc.list(lapply(chains, coda::mcmc)),
    autoburnin = FALSE, multivariate = FALSE)$psrf[1L, 1L]
}

test_that("the split-chain point estimate separates agreeing, shifted and drifting chains", {

  set.seed(11)
  agreeing <- replicate(3L, stats::rnorm(500), simplify = FALSE)
  expect_lt(.burnin_rhat(agreeing), 1.1)

  shifted <- agreeing
  shifted[[3L]] <- shifted[[3L]] + 3
  expect_gt(.burnin_rhat(shifted), 1.5)

  # a common drift of 6 standard deviations over the block: the chains agree
  # with each other, so only the split halves reveal it
  drift <- lapply(agreeing, function(chain) chain + seq(-3, 3, length.out = 500))
  expect_lt(.burnin_unsplit_point(drift), 1.1)
  expect_gt(.burnin_rhat(drift), 1.5)

  # the point estimate of coda over the 2 x chains half-chains
  halves <- unlist(lapply(shifted, function(chain) list(chain[1:250], chain[251:500])), recursive = FALSE)
  expect_equal(.burnin_rhat(shifted), .burnin_unsplit_point(halves), tolerance = 1e-12)
})

test_that("degenerate blocks are decided before coda", {

  # every half-chain constant, but at different values
  expect_identical(.burnin_rhat(list(rep(0, 100), rep(1, 100), rep(0, 100))), Inf)
  expect_identical(.burnin_rhat(list(c(rep(0, 50), rep(1, 50)), rep(0, 100))), Inf)

  # half-chains that vary with equal means (identical state counts in every
  # half): coda gives NaN, the chains agree
  occupancy <- c(rep(0, 40), rep(1, 10))
  equal <- list(c(occupancy, rev(occupancy)), c(rev(occupancy), occupancy), c(occupancy, occupancy))
  expect_true(is.nan(.burnin_unsplit_point(unlist(lapply(equal, function(chain){
    list(chain[1:50], chain[51:100])
  }), recursive = FALSE))))
  expect_identical(.burnin_rhat(equal), 1)

  # nothing to assess: constant everywhere
  expect_identical(.burnin_rhat(list(rep(2, 100), rep(2, 100))), NA_real_)

  # a value that is not finite in some half-chains only: the chains do not
  # agree; in every half-chain: not assessed
  expect_identical(.burnin_rhat(list(c(stats::rnorm(99), NA), stats::rnorm(100))), Inf)
  expect_identical(.burnin_rhat(list(c(stats::rnorm(99), Inf), stats::rnorm(100))), Inf)
  expect_identical(.burnin_rhat(list(c(-Inf, stats::rnorm(99)), rep(NaN, 100))), Inf)
  expect_identical(.burnin_rhat(list(rep(Inf, 100), rep(Inf, 100))), NA_real_)
  expect_identical(.burnin_rhat(list(rep(c(1, NA), 50), rep(c(NaN, 2), 50))), NA_real_)
})

test_that("an odd block leaves the middle draw out of both halves", {

  set.seed(3)
  chains <- list(stats::rnorm(101), stats::rnorm(101) + 1)
  halves <- unlist(lapply(chains, function(chain) list(chain[1:50], chain[52:101])), recursive = FALSE)
  expect_equal(.burnin_rhat(chains), .burnin_unsplit_point(halves), tolerance = 1e-12)

  # the degenerate cases are decided on the halves: a block that varies only
  # at its dropped middle draw is constant, and so is one that is not finite
  # only there
  constant <- rep(0, 101)
  expect_identical(.burnin_rhat(list(replace(constant, 51, 1), constant)), NA_real_)
  expect_identical(.burnin_rhat(list(replace(constant, 51, NA), constant)), NA_real_)
  expect_identical(.burnin_rhat(list(replace(constant, 51, 1), constant + 1)), Inf)
})

# A block in the form the backend returns it, with the coordinates that
# JAGS_fit() builds for its columns.
.burnin_block <- function(chains, prior_list){

  block <- coda::mcmc.list(lapply(chains, coda::mcmc))
  columns <- colnames(chains[[1L]])
  list(
    block = block,
    coordinates = BayesTools:::.bt_build_parameter_coordinates(
      columns = columns,
      monitor_names = unique(sub("\\[.*$", "", columns)),
      prior_list = prior_list
    )
  )
}

test_that("a block is assessed on the selection of the convergence check", {

  set.seed(5)
  priors <- list(
    mu    = prior("normal", list(0, 1)),
    theta = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1)))
  )
  chain <- function(indicator){
    cbind(mu = stats::rnorm(200), theta = stats::rnorm(200),
      theta_indicator = indicator, theta_inclusion = stats::runif(200))
  }
  # the indicator disagrees: every chain constant at its own state
  input <- .burnin_block(list(chain(0), chain(0), chain(1)), priors)
  assess <- function(...){
    BayesTools:::.bt_burnin_block_Rhat(input$block, input$coordinates, priors, ...)
  }

  default <- assess(monitor = NULL, check_indicators = FALSE)
  expect_setequal(names(default), c("mu", "theta"))
  expect_true(all(default < 1.5))

  indicators <- assess(monitor = NULL, check_indicators = TRUE)
  expect_setequal(names(indicators), c("mu", "theta", "theta_indicator"))
  expect_identical(unname(indicators[["theta_indicator"]]), Inf)

  expect_identical(names(assess(monitor = "mu", check_indicators = FALSE)), "mu")
  expect_setequal(names(assess(monitor = "mu", check_indicators = TRUE)), c("mu", "theta_indicator"))
  # an auxiliary coordinate is assessed when requested
  expect_identical(names(assess(monitor = "theta_inclusion", check_indicators = FALSE)), "theta_inclusion")
  expect_error(
    assess(monitor = "mu[2]", check_indicators = FALSE),
    class = "BayesTools_parameter_not_found"
  )

  # categorical indicators are assessed per observed state
  states <- list(rep(1:3, length.out = 200), rep(c(1, 2, 3, 3), length.out = 200), rep(1:3, length.out = 200))
  input <- .burnin_block(lapply(states, function(state){
    cbind(mu = stats::rnorm(200), theta = stats::rnorm(200), theta_indicator = state - 1,
      theta_inclusion = stats::runif(200))
  }), priors)
  occupancies <- assess(monitor = "theta_indicator", check_indicators = FALSE)
  expect_length(occupancies, 3L)
  expect_true(all(grepl("^theta_indicator \\(state [0-2]\\)$", names(occupancies))))

  # nothing assessable
  input <- .burnin_block(replicate(3L, cbind(mu = rep(1, 200)), simplify = FALSE), priors["mu"])
  expect_length(assess(monitor = NULL, check_indicators = FALSE), 0L)
})

test_that("a requested indicator state is resolved whatever states the block visits", {

  set.seed(6)
  priors <- list(
    mu    = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1)), prior("normal", list(2, 1)))),
    theta = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1)))
  )
  chain <- function(states){
    cbind(mu = stats::rnorm(200), mu_indicator = rep(states, length.out = 200), theta = stats::rnorm(200),
      theta_indicator = rep(0:1, length.out = 200), theta_inclusion = stats::runif(200))
  }
  # the block visits states 1 and 2 only: its occupancy target carries the
  # indicator's own name
  chains <- list(chain(1:2), chain(c(2, 1, 1)), chain(c(1, 2, 2, 2)))
  input <- .burnin_block(chains, priors)
  assess <- function(monitor){
    BayesTools:::.bt_burnin_block_Rhat(input$block, input$coordinates, priors, monitor, FALSE)
  }
  occupancy_rhat <- function(state){
    unname(.burnin_rhat(lapply(chains, function(chain) as.numeric(chain[, "mu_indicator"] == state))))
  }
  expect_identical(names(assess("mu_indicator")), "mu_indicator")

  # a declared state that the block has not visited has a constant occupancy:
  # not assessable, not an error
  expect_length(assess("mu_indicator (state 3)"), 0L)
  expect_identical(assess("mu_indicator (state 2)"), c("mu_indicator (state 2)" = occupancy_rhat(2)))
  several <- assess(c("mu_indicator (state 1)", "mu_indicator (state 3)", "mu"))
  expect_setequal(names(several), c("mu", "mu_indicator (state 1)"))
  expect_identical(several[["mu_indicator (state 1)"]], occupancy_rhat(1))
  expect_identical(several[["mu"]], unname(.burnin_rhat(lapply(chains, function(chain) chain[, "mu"]))))

  # a selector that names no categorical indicator or a state it cannot take
  for(monitor in c("mu_indicator (state 4)", "mu_indicator (state 0)", "mu_indicator (state 1.5)",
                   "nu_indicator (state 1)", "theta_indicator (state 1)")){
    expect_error(assess(monitor), class = "BayesTools_parameter_not_found", info = monitor)
  }

  # a block that visits every state keeps the targets of the convergence check
  chains <- list(chain(1:3), chain(c(3, 1, 2)), chain(c(1, 2, 3, 3)))
  input <- .burnin_block(chains, priors)
  expect_identical(assess("mu_indicator (state 3)"), c("mu_indicator (state 3)" = occupancy_rhat(3)))
  expect_setequal(names(assess("mu_indicator")), paste0("mu_indicator (state ", 1:3, ")"))
})

# The extension loop with scripted assessments.
.burnin_loop <- function(values, autofit = TRUE, max_extensions = 10, time_exceeded = function() FALSE,
                         extension_room = Inf){

  calls <- new.env()
  calls$assess <- 0L
  calls$blocks <- integer()
  record <- BayesTools:::.bt_burnin_readiness(
    first = "first block",
    run_block = function(sweeps){
      calls$blocks <- c(calls$blocks, sweeps)
      paste("block", length(calls$blocks))
    },
    assess = function(draws){
      calls$assess <- calls$assess + 1L
      calls$draws <- c(calls$draws, draws)
      values[[min(calls$assess, length(values))]]
    },
    burnin_max_Rhat = 1.5,
    block = 400L,
    burnin_extend = 500,
    max_extensions = max_extensions,
    autofit = autofit,
    time_exceeded = time_exceeded,
    extension_room = extension_room
  )
  list(record = record, calls = calls)
}

test_that("the extension loop stops when ready, at the cap, at max_time and without autofit", {

  # ready at once
  run <- .burnin_loop(list(c(a = 1.2, b = 1.1)))
  expect_identical(run$record$status, "ready")
  expect_identical(run$record$stop_reason, "ready")
  expect_identical(run$record$extension, 0)
  expect_length(run$calls$blocks, 0L)
  expect_identical(run$record$blocks,
    data.frame(block = 0L, sweeps = 400L, Rhat = 1.2, parameter = "a", stringsAsFactors = FALSE))
  expect_identical(nrow(run$record$failing), 0L)

  # ready after two extensions, each block checked alone
  run <- .burnin_loop(list(c(a = 3, b = 2), c(a = 1.2, b = 1.7), c(a = 1.1, b = 1.4)))
  expect_identical(run$record$status, "ready")
  expect_identical(run$record$extension, 1000)
  expect_identical(run$calls$blocks, c(500, 500))
  expect_identical(run$calls$draws, c("first block", "block 1", "block 2"))
  expect_identical(run$record$blocks$Rhat, c(3, 1.7, 1.4))
  expect_identical(run$record$blocks$parameter, c("a", "b", "b"))
  expect_identical(run$record$blocks$sweeps, c(400L, 500L, 500L))

  # the cap
  run <- .burnin_loop(list(c(a = 3, b = 2, c = 1.2)), max_extensions = 2)
  expect_identical(run$record$status, "not_ready")
  expect_identical(run$record$stop_reason, "max_extend")
  expect_identical(run$record$extension, 1000)
  expect_identical(run$record$failing,
    data.frame(parameter = c("a", "b"), Rhat = c(3, 2), stringsAsFactors = FALSE))
  expect_identical(.burnin_loop(list(c(a = 3)), max_extensions = 0)$record$stop_reason, "max_extend")

  # max_time, read before each extension block is started
  reads <- 0L
  run <- .burnin_loop(list(c(a = 3)), time_exceeded = function(){
    reads <<- reads + 1L
    reads > 1L
  })
  expect_identical(run$record$stop_reason, "max_time")
  expect_identical(run$record$extension, 500)
  expect_identical(reads, 2L)

  # the integer range of the iteration numbers, checked before each block and
  # before the clock
  run <- .burnin_loop(list(c(a = 3)), extension_room = 1000, time_exceeded = function() FALSE)
  expect_identical(run$record[c("status", "stop_reason", "extension")],
    list(status = "not_ready", stop_reason = "iteration_limit", extension = 1000))
  run <- .burnin_loop(list(c(a = 3)), extension_room = 999)
  expect_identical(run$record[c("stop_reason", "extension")], list(stop_reason = "iteration_limit", extension = 500))
  run <- .burnin_loop(list(c(a = 3)), extension_room = 499, time_exceeded = function() stop("clock read"))
  expect_identical(run$record$stop_reason, "iteration_limit")
  expect_identical(run$calls$blocks, integer())
  expect_identical(.burnin_loop(list(c(a = 3)), max_extensions = 1, extension_room = 500)$record$stop_reason,
    "max_extend")

  # without autofit: one check, no extension and no clock
  run <- .burnin_loop(list(c(a = 3)), autofit = FALSE, time_exceeded = function() stop("clock read"))
  expect_identical(run$record$status, "not_ready")
  expect_identical(run$record$stop_reason, "no_autofit")
  expect_length(run$calls$blocks, 0L)

  # nothing assessable: no extension and no claim of readiness
  run <- .burnin_loop(list(numeric()))
  expect_identical(run$record$status, "not_assessable")
  expect_identical(run$record$stop_reason, "not_assessable")
  expect_length(run$calls$blocks, 0L)
  expect_identical(run$record$blocks$Rhat, NA_real_)
  expect_null(BayesTools:::.bt_burnin_warning(run$record))
})

test_that("the warning names the quantities with the largest values and the remedy", {

  values <- c(p1 = Inf, p2 = 3.14159, p3 = 2.5, p4 = 2, p5 = 1.8, p6 = 1.6, p7 = 1.2)
  record <- .burnin_loop(list(values), autofit = FALSE)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("The chains had not reached agreement at the end of burn-in: the split-chain R-hat of the ",
      "last 400 burn-in iterations was Inf for 'p1', 3.14 for 'p2', 2.50 for 'p3', 2.00 for 'p4', ",
      "1.80 for 'p5' and 1 more (above 'burnin_max_Rhat' of 1.5). Increase 'burnin', or fit with ",
      "'autofit = TRUE'.")
  )
  record <- .burnin_loop(list(c(mu = 2)), max_extensions = 1)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("The chains had not reached agreement at the end of burn-in after 1 extension of 500 ",
      "iterations: the split-chain R-hat of the last 500 burn-in iterations was 2.00 for 'mu' (above ",
      "'burnin_max_Rhat' of 1.5). Increase 'burnin' or 'burnin_max_extend'.")
  )
  record <- .burnin_loop(list(c(mu = 2)), time_exceeded = function() TRUE)$record
  expect_match(BayesTools:::.bt_burnin_warning(record), "Increase 'burnin' or 'max_time'.", fixed = TRUE)
  record <- .burnin_loop(list(c(mu = 2)), extension_room = 700)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("The chains had not reached agreement at the end of burn-in after 1 extension of 500 ",
      "iterations: the split-chain R-hat of the last 500 burn-in iterations was 2.00 for 'mu' (above ",
      "'burnin_max_Rhat' of 1.5). A further extension would number the draws beyond the largest ",
      "iteration number (2147483647); reduce 'adapt' or the sampling iterations ('sample' times 'thin') ",
      "to leave room for it.")
  )
  expect_null(BayesTools:::.bt_burnin_warning(.burnin_loop(list(c(mu = 1.2)))$record))
  expect_null(BayesTools:::.bt_burnin_warning(NULL))
})

test_that("the readiness settings default when omitted and validate when given", {

  base <- list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10)
  settings <- function(...) BayesTools:::.bt_burnin_check_settings(
    JAGS_check_and_list_autofit_settings(utils::modifyList(base, list(...), keep.null = TRUE)))

  expect_identical(settings(), list(burnin_max_Rhat = 1.5, burnin_extend = 500, burnin_max_extend = 10))
  # the validator does not add the elements, so that an omission stays one
  expect_false(any(grepl("^burnin_", names(JAGS_check_and_list_autofit_settings(base)))))
  expect_null(settings(burnin_max_Rhat = NULL)$burnin_max_Rhat)
  expect_identical(settings(burnin_max_Rhat = 1.2, burnin_extend = 100, burnin_max_extend = 0),
    list(burnin_max_Rhat = 1.2, burnin_extend = 100, burnin_max_extend = 0))
  expect_identical(
    JAGS_check_and_list_autofit_settings(c(base, list(burnin_max_Rhat = NULL)), skip_sample_extend = TRUE)$burnin_max_Rhat,
    NULL
  )

  invalid <- list(
    list(burnin_max_Rhat = 0.9, "'burnin_max_Rhat' must be equal or higher than 1."),
    list(burnin_max_Rhat = Inf, "The 'burnin_max_Rhat' argument must contain only finite values."),
    list(burnin_max_Rhat = NA_real_, "The 'burnin_max_Rhat' argument cannot contain NA/NaN values."),
    list(burnin_extend = 99, "'burnin_extend' must be equal or higher than 100."),
    list(burnin_extend = NULL, "The 'burnin_extend' argument cannot be NULL."),
    list(burnin_extend = 250.5, "The 'burnin_extend' argument must be an integer vector."),
    list(burnin_max_extend = -1, "'burnin_max_extend' must be equal or higher than 0."),
    list(burnin_max_extend = Inf, "The 'burnin_max_extend' argument must contain only finite values."),
    list(burnin_max_extend = NULL, "The 'burnin_max_extend' argument cannot be NULL.")
  )
  for(case in invalid){
    control <- c(base, case[1L])
    expect_error(JAGS_check_and_list_autofit_settings(control), case[[2L]], fixed = TRUE, info = names(case)[[1L]])
  }
})

# JAGS_fit() with the backend mocked below the readiness check: the chain
# groups return the scripted burn-in blocks and the sampled draws.
.readiness_backend <- function(blocks, chains = 3L, fail_batches = 0L, env = parent.frame()){

  state <- new.env()
  state$groups <- list()
  state$sweeps <- numeric()
  state$batches <- 0L
  state$block_index <- 0L
  next_block <- function(){
    state$block_index <- state$block_index + 1L
    blocks[[min(state$block_index, length(blocks))]]
  }
  sampled <- function(sample, thin){
    list(
      mcmc = coda::mcmc.list(lapply(seq_len(chains), function(chain){
        coda::mcmc(matrix(stats::rnorm(sample), ncol = 1L, dimnames = list(NULL, "mu")), thin = thin)
      })),
      end_state = lapply(seq_len(chains), function(chain){
        list(.RNG.name = "base::Wichmann-Hill", .RNG.state = c(chain, chain + 1L, chain + 2L))
      }),
      samplers = data.frame(Index.No = integer(), Sampler = character(), Node = character()),
      live_model = "compiled"
    )
  }
  batch <- function(sample, thin){
    state$batches <- state$batches + 1L
    if(state$batches <= fail_batches) stop("sampling exploded", call. = FALSE)
    sampled(sample, thin)
  }
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(...) invisible(NULL),
    .JAGS_load_modules = function(...) invisible(NULL),
    .bt_rjags_group = function(inits, model, data, monitor, adapt, burnin, sample, thin, silent,
                               keep_live, live_model = NULL, probability_plan = NULL,
                               prepare = FALSE, block = 0L, session){
      state$groups[[length(state$groups) + 1L]] <- list(chains = length(inits), burnin = burnin,
        block = block, prepare = prepare)
      if(!prepare) return(batch(sample, thin))
      session$model <- "compiled"
      if(block > 0L) next_block() else NULL
    },
    .bt_rjags_burnin_block = function(sweeps, monitor, progress, session){
      state$sweeps <- c(state$sweeps, sweeps)
      next_block()
    },
    .bt_rjags_batch = function(sample, monitor, thin, ...) batch(sample, thin),
    .package = "BayesTools", .env = env
  )
  state
}

.readiness_chains <- function(shift = 0, chains = 3L, n = 500L){

  coda::mcmc.list(lapply(seq_len(chains), function(chain){
    coda::mcmc(matrix(stats::rnorm(n) + if(chain == chains) shift else 0, ncol = 1L,
      dimnames = list(NULL, "mu")))
  }))
}

.readiness_fit <- function(autofit = FALSE, burnin = 600, chains = 3L, thin = 1L, adapt = 100, ...){

  control <- utils::modifyList(list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL, max_time = list(time = 60, unit = "mins"), sample_extend = 100,
    restarts = 2, max_extend = 1), list(...), keep.null = TRUE)
  JAGS_fit("model{ x ~ dnorm(mu, 1) }", data = list(x = 0),
    prior_list = list(mu = prior("normal", list(0, 1))), chains = chains, adapt = adapt,
    burnin = burnin, sample = 200, thin = thin, seed = 1, silent = TRUE, autofit = autofit,
    autofit_control = control)
}

test_that("the check runs on the last block of burn-in when it can and the fit records it", {

  set.seed(1)
  state <- .readiness_backend(list(.readiness_chains()))
  fit <- .readiness_fit()
  expect_identical(state$groups[[1L]][c("chains", "burnin", "block", "prepare")],
    list(chains = 3L, burnin = 600, block = 500L, prepare = TRUE))
  record <- attr(fit, "burnin_check")
  expect_identical(record$status, "ready")
  expect_identical(record$block, 500L)
  expect_identical(fit$burnin, 600)
  expect_identical(stats::start(fit$mcmc[[1L]]), 701)
  expect_null(attr(fit, "warnings"))
  expect_null(fit$live_model)

  # a burn-in shorter than a block is checked whole
  state <- .readiness_backend(list(.readiness_chains(n = 150L)))
  fit <- .readiness_fit(burnin = 150)
  expect_identical(state$groups[[1L]]$block, 150L)
  expect_identical(attr(fit, "burnin_check")$block, 150L)

  # no check: burn-in below 100 sweeps, one chain, or switched off; the fit
  # then runs as before, in one call of the chain group
  for(case in list(list(burnin = 99), list(chains = 1L), list(burnin_max_Rhat = NULL))){
    state <- .readiness_backend(list(), chains = if(is.null(case$chains)) 3L else case$chains)
    fit <- do.call(.readiness_fit, case)
    expect_identical(state$groups[[1L]]$block, 0L)
    expect_false(state$groups[[1L]]$prepare)
    expect_null(attr(fit, "burnin_check"))
  }
})

test_that("a model without samplers is not checked", {

  # the chain group returns no block when the compiled model has no sampler
  state <- .readiness_backend(list(NULL))
  fit <- .readiness_fit(autofit = TRUE)
  expect_identical(state$groups[[1L]][c("block", "prepare")], list(block = 500L, prepare = TRUE))
  expect_length(state$sweeps, 0L)
  expect_null(attr(fit, "burnin_check"))
  expect_identical(fit$burnin, 600)
  expect_identical(stats::start(fit$mcmc[[1L]]), 701)

  # the processes of a parallel fit compile the same graph and must agree
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, ...) lapply(x, fun, ...),
    clusterCall = function(cl, fun, ...) lapply(cl, function(worker) fun(...)),
    .package = "parallel"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_group = function(inits, ...){
      if(length(inits) == 2L) .readiness_chains(chains = 2L, n = 200L) else NULL
    },
    .package = "BayesTools"
  )
  readiness <- list(block = 200L, burnin_max_Rhat = 1.5, burnin_extend = 500, max_extensions = 0,
    autofit = FALSE, time_exceeded = function() FALSE, assess = function(draws) c(mu = 1))
  expect_error(BayesTools:::.bt_rjags_run("", list(), list(1L, 2L, 3L), "mu", 100L, 300L, 100L, 1L,
    cl = list(1L, 2L), silent = TRUE, readiness = readiness),
    "Internal error: the JAGS processes disagree on whether the model has samplers.", fixed = TRUE)
})

test_that("chains that disagree without autofit sample after one warning", {

  set.seed(2)
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  expect_warning(
    fit <- .readiness_fit(),
    "The chains had not reached agreement at the end of burn-in: the split-chain R-hat",
    class = "BayesTools_burnin_not_converged"
  )
  expect_length(state$sweeps, 0L)
  record <- attr(fit, "burnin_check")
  expect_identical(record[c("status", "stop_reason", "extension")],
    list(status = "not_ready", stop_reason = "no_autofit", extension = 0))
  expect_identical(fit$burnin, 600)
  expect_identical(attr(fit, "warnings"), BayesTools:::.bt_burnin_warning(record))
  expect_true(inherits(
    tryCatch(BayesTools:::.bt_burnin_warn(record), warning = identity),
    "BayesTools_fit_diagnostic"
  ))
})

test_that("autofit extends burn-in in blocks and the iteration numbers count them", {

  set.seed(3)
  state <- .readiness_backend(list(.readiness_chains(shift = 10), .readiness_chains(shift = 10),
    .readiness_chains()))
  fit <- .readiness_fit(autofit = TRUE, thin = 2L, burnin_extend = 300)
  expect_identical(state$groups[[1L]]$block, 300L)
  expect_identical(state$sweeps, c(300, 300))
  record <- attr(fit, "burnin_check")
  expect_identical(record$status, "ready")
  expect_identical(record$extension, 600)
  expect_identical(record$blocks$sweeps, c(300L, 300L, 300L))
  expect_identical(fit$burnin, 1200)
  expect_identical(stats::start(fit$mcmc[[1L]]), 100 + 1200 + 1)
  expect_identical(coda::thin(fit$mcmc[[1L]]), 2)
  geometry <- JAGS_draw_geometry(fit)
  expect_identical(geometry$chains$start, rep(1301L, 3L))
  expect_identical(geometry$chains$end, rep(1301L + 199L * 2L, 3L))
  expect_null(attr(fit, "warnings"))

  # the cap: the warning names the extensions
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  expect_warning(
    fit <- .readiness_fit(autofit = TRUE, burnin_max_extend = 2),
    "after 2 extensions of 500 iterations",
    class = "BayesTools_burnin_not_converged"
  )
  expect_identical(attr(fit, "burnin_check")$stop_reason, "max_extend")
  expect_identical(fit$burnin, 1600)
  expect_length(attr(fit, "warnings"), 1L)

  # max_time on the fit's clock, read before each block
  clock <- 0
  testthat::local_mocked_bindings(
    .bt_autofit_clock = function(){
      clock <<- clock + 60
      as.POSIXct(clock, origin = "1970-01-01")
    },
    .package = "BayesTools"
  )
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  expect_warning(
    fit <- .readiness_fit(autofit = TRUE, max_time = list(time = 100, unit = "secs")),
    "Increase 'burnin' or 'max_time'.",
    fixed = TRUE
  )
  expect_identical(attr(fit, "burnin_check")$stop_reason, "max_time")
  expect_identical(state$sweeps, 500)
})

test_that("no extension numbers the draws beyond the integer range", {

  limit <- .Machine$integer.max

  # burn-in: one block of 500 fits before the last draw (adapt + burnin +
  # extension + 200) would pass the limit, a second does not
  set.seed(7)
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  expect_warning(
    fit <- .readiness_fit(autofit = TRUE, adapt = limit - 1500, burnin = 600),
    "A further extension would number the draws beyond the largest iteration number",
    class = "BayesTools_burnin_not_converged"
  )
  record <- attr(fit, "burnin_check")
  expect_identical(record[c("stop_reason", "extension")], list(stop_reason = "iteration_limit", extension = 500))
  expect_identical(state$sweeps, 500)
  expect_identical(fit$burnin, 1100)
  expect_identical(JAGS_draw_geometry(fit)$chains$end, rep(limit - 200L, 3L))

  # sampling: the first automatic extension reaches the limit exactly, the
  # second would pass it and is not started; the fit keeps its draws
  message <- paste0("The automatic model fitting was terminated due to the iteration limit: an extension of ",
    "200 draws per chain would number the draws beyond the largest iteration number (2147483647). ",
    "Refit with a smaller 'adapt', 'burnin', or 'sample' times 'thin' to leave room for extensions.")
  state <- .readiness_backend(list(.readiness_chains()))
  fit <- .readiness_fit(autofit = TRUE, adapt = limit - 1000, burnin = 600, min_ESS = 1e6,
    sample_extend = 200, max_extend = 5)
  expect_identical(state$batches, 2L)
  expect_identical(NROW(fit$mcmc[[1L]]), 400L)
  expect_identical(JAGS_draw_geometry(fit)$chains$end, rep(limit, 3L))
  expect_identical(attr(fit, "warnings"), message)

  # an "auto" extension of at least 1000 draws does not fit at all
  state <- .readiness_backend(list(.readiness_chains()))
  fit <- .readiness_fit(autofit = TRUE, adapt = limit - 1000, burnin = 600, min_ESS = 1e6,
    sample_extend = "auto", max_extend = 5, max_time = NULL)
  expect_identical(state$batches, 1L)
  expect_identical(NROW(fit$mcmc[[1L]]), 200L)
  expect_match(attr(fit, "warnings"), "an extension of 1000 draws per chain", fixed = TRUE)

  # JAGS_extend(): an extension of 300 draws would pass the limit (the last
  # draw is at limit - 200), so it is not started and the draws are kept
  expect_warning(
    extended <- JAGS_extend(fit, autofit_control = list(max_Rhat = NULL, min_ESS = 1e6, max_error = NULL,
      max_SD_error = NULL, max_time = NULL, sample_extend = 300, max_extend = 5), silent = FALSE),
    sub("200 draws", "300 draws", message, fixed = TRUE), fixed = TRUE
  )
  expect_identical(state$batches, 1L)
  expect_identical(extended$mcmc, fit$mcmc)
  expect_true(sub("200 draws", "300 draws", message, fixed = TRUE) %in% attr(extended, "warnings"))
})

test_that("a restart carries no state or warning of the failed attempt", {

  set.seed(4)
  state <- .readiness_backend(list(.readiness_chains(shift = 10), .readiness_chains()), fail_batches = 1L)
  expect_no_warning(fit <- .readiness_fit(restarts = 2))
  expect_length(state$groups, 2L)
  expect_identical(attr(fit, "burnin_check")$status, "ready")
  expect_identical(attr(fit, "warnings"),
    "JAGS fitting attempt 1 failed and was restarted: sampling exploded.")
})

test_that("the convergence monitor is validated also without autofit", {

  state <- .readiness_backend(list(.readiness_chains()))
  expect_error(.readiness_fit(monitor = "nu"), class = "BayesTools_parameter_not_found")
  expect_length(state$groups, 0L)
  # without the check, a fit with autofit = FALSE does not read the monitor
  expect_s3_class(.readiness_fit(monitor = "nu", burnin_max_Rhat = NULL), "BayesTools_fit")

  # an indexed element is resolved on the block, and its failure is not
  # retried as a sampler failure
  state <- .readiness_backend(list(.readiness_chains()))
  expect_error(.readiness_fit(monitor = "mu[2]", restarts = 5), class = "BayesTools_parameter_not_found")
  expect_length(state$groups, 1L)
})

test_that("parallel blocks return to the chain order and workers are released after a failure", {

  chains <- lapply(1:5, function(chain) coda::mcmc(matrix(chain, nrow = 2L, dimnames = list(NULL, "mu"))))
  groups <- list(c(1L, 3L, 5L), c(2L, 4L))
  merged <- BayesTools:::.bt_rjags_merge_blocks(
    list(coda::mcmc.list(chains[c(1L, 3L, 5L)]), coda::mcmc.list(chains[c(2L, 4L)])), groups)
  expect_identical(vapply(merged, function(chain) chain[1L, 1L], numeric(1)), as.numeric(1:5))

  released <- 0L
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, ...) lapply(x, fun, ...),
    clusterCall = function(cl, fun, ...){
      released <<- released + 1L
      lapply(cl, function(worker) fun(...))
    },
    .package = "parallel"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_group = function(inits, ...) .readiness_chains(chains = length(inits), n = 200L),
    .bt_rjags_batch = function(...) stop("worker failed"),
    .package = "BayesTools"
  )
  readiness <- list(block = 200L, burnin_max_Rhat = 1.5, burnin_extend = 500, max_extensions = 0,
    autofit = FALSE, time_exceeded = function() FALSE, assess = function(draws){
      expect_length(draws, 3L)
      c(mu = 1)
    })
  expect_error(BayesTools:::.bt_rjags_run("", list(), list(1L, 2L, 3L), "mu", 100L, 300L, 100L, 1L,
    cl = list(1L, 2L), silent = TRUE, readiness = readiness), "worker failed")
  expect_identical(released, 1L)

  # a lost connection is not followed by a release call
  testthat::local_mocked_bindings(
    .bt_rjags_batch = function(...) stop("error reading from connection"),
    .package = "BayesTools"
  )
  expect_error(BayesTools:::.bt_rjags_run("", list(), list(1L, 2L, 3L), "mu", 100L, 300L, 100L, 1L,
    cl = list(1L, 2L), silent = TRUE, readiness = readiness), "connection")
  expect_identical(released, 1L)
})

test_that("the JASP warm-up bar ticks after the check and its extensions", {

  events <- character()
  testthat::local_mocked_bindings(
    .bt_rjags_group = function(inits, ..., block, session){
      events <<- c(events, paste("warm-up block", block))
      .readiness_chains(shift = 10, n = block)
    },
    .bt_rjags_burnin_block = function(sweeps, ...){
      events <<- c(events, paste("extension", sweeps))
      .readiness_chains(n = sweeps)
    },
    .bt_rjags_batch = function(sample, ...){
      events <<- c(events, paste("batch", sample))
      list(mcmc = .readiness_chains(n = sample), end_state = list(), samplers = data.frame())
    },
    .JASP_progress_bar_start = function(n, label) events <<- c(events, paste("start", label)),
    .JASP_progress_bar_tick = function() events <<- c(events, "tick"),
    .package = "BayesTools"
  )
  readiness <- list(block = 200L, burnin_max_Rhat = 1.5, burnin_extend = 300, max_extensions = 5, extension_room = Inf,
    autofit = TRUE, time_exceeded = function() FALSE, assess = function(draws){
      c(mu = .burnin_rhat(lapply(draws, as.numeric)))
    })
  out <- BayesTools:::.bt_rjags_run("", list(), list(1L, 2L, 3L), "mu", 100L, 300L, 10L, 1L,
    cl = NULL, silent = TRUE, jasp_prefix = "", readiness = readiness)
  expect_identical(events[1:5], c("start Adapting and burnin the model", "warm-up block 200",
    "extension 300", "tick", "start Sampling the model"))
  expect_identical(out$burnin_check$extension, 300)
})
