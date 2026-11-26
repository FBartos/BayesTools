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
                         extension_room = Inf, thin = 1L, burnin_max_Rhat = 1.5, burnin_extend = 500){

  calls <- new.env()
  calls$assess <- 0L
  calls$blocks <- integer()
  record <- BayesTools:::.bt_burnin_readiness(
    first = "first block",
    run_block = function(draws){
      calls$blocks <- c(calls$blocks, draws)
      paste("block", length(calls$blocks))
    },
    assess = function(draws){
      calls$assess <- calls$assess + 1L
      calls$draws <- c(calls$draws, draws)
      values[[min(calls$assess, length(values))]]
    },
    burnin_max_Rhat = burnin_max_Rhat,
    block = 400L,
    burnin_extend = burnin_extend,
    thin = thin,
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

test_that("with thinning the loop runs blocks of draws and records iterations", {

  # the first block holds 400 draws, every third iteration, and each
  # extension 500 draws: 1500 iterations
  run <- .burnin_loop(list(c(a = 3), c(a = 3), c(a = 1.2)), thin = 3L)
  expect_identical(run$calls$blocks, c(500, 500))
  expect_identical(run$record[c("block", "burnin_extend", "thin", "extension", "status")],
    list(block = 1200L, burnin_extend = 500, thin = 3L, extension = 3000, status = "ready"))
  expect_identical(run$record$blocks$sweeps, c(1200L, 1500L, 1500L))

  # the room for iteration numbers counts iterations: two blocks of 1500 fit
  # in 3000, a third does not; none fits in 1499
  run <- .burnin_loop(list(c(a = 3)), thin = 3L, extension_room = 3000)
  expect_identical(run$record[c("stop_reason", "extension")], list(stop_reason = "iteration_limit", extension = 3000))
  expect_identical(run$calls$blocks, c(500, 500))
  run <- .burnin_loop(list(c(a = 3)), thin = 3L, extension_room = 1499)
  expect_identical(run$record[c("stop_reason", "extension")], list(stop_reason = "iteration_limit", extension = 0))
  expect_length(run$calls$blocks, 0L)

  # 'extension' keeps the type of 'burnin_extend', also at the default
  # thin = 1 (a double)
  expect_identical(.burnin_loop(list(c(a = 3), c(a = 1.2)), burnin_extend = 500L, thin = 1)$record$extension, 500L)
  expect_identical(.burnin_loop(list(c(a = 1.2)), burnin_extend = 500L, thin = 1)$record$extension, 0L)
  expect_identical(.burnin_loop(list(c(a = 3), c(a = 1.2)), burnin_extend = 500L, thin = 3)$record$extension, 1500L)
  expect_identical(.burnin_loop(list(c(a = 3), c(a = 1.2)), burnin_extend = 500, thin = 1L)$record$extension, 500)
})

test_that("the warning names the quantities with the largest values and the remedy", {

  values <- c(p1 = Inf, p2 = 3.14159, p3 = 2.5, p4 = 2, p5 = 1.8, p6 = 1.6, p7 = 1.2)
  record <- .burnin_loop(list(values), autofit = FALSE)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("Chains did not agree at the end of burn-in: R-hat Inf for 'p1', 3.14 for 'p2', 2.50 for 'p3', ",
      "2.00 for 'p4', 1.80 for 'p5' and 1 more (limit 1.50). Increase 'adapt' and 'burnin'; if the chains ",
      "mix slowly, increase 'thin'. Or fit with 'autofit = TRUE'.")
  )
  # a value just above the threshold does not print as the threshold
  record <- .burnin_loop(list(c(mu = 1.50001, tau = 1.5049, nu = 1.5051)), autofit = FALSE)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("Chains did not agree at the end of burn-in: R-hat 1.51 for 'nu', 1.505 for 'tau', 1.50001 for ",
      "'mu' (limit 1.50). Increase 'adapt' and 'burnin'; if the chains mix slowly, increase 'thin'. Or fit ",
      "with 'autofit = TRUE'.")
  )
  # stopped by 'burnin_max_extend' or by 'max_time'
  record <- .burnin_loop(list(c(mu = 2)), max_extensions = 1)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("Chains did not agree at the end of burn-in after 1 extension: R-hat 2.00 for 'mu' (limit 1.50). ",
      "Increase 'adapt', 'burnin' or 'thin'.")
  )
  record <- .burnin_loop(list(c(mu = 2.31)), thin = 4L)$record
  expect_identical(record$stop_reason, "max_extend")
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("Chains did not agree at the end of burn-in after 10 extensions: R-hat 2.31 for 'mu' (limit 1.50). ",
      "Increase 'adapt', 'burnin' or 'thin'.")
  )
  reads <- 0L
  record <- .burnin_loop(list(c(mu = 2)), time_exceeded = function(){
    reads <<- reads + 1L
    reads > 2L
  })$record
  expect_identical(record$stop_reason, "max_time")
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("Chains did not agree at the end of burn-in after 2 extensions: R-hat 2.00 for 'mu' (limit 1.50). ",
      "Increase 'adapt', 'burnin' or 'thin'.")
  )
  record <- .burnin_loop(list(c(mu = 2)), time_exceeded = function() TRUE)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    "Chains did not agree at the end of burn-in: R-hat 2.00 for 'mu' (limit 1.50). Increase 'adapt', 'burnin' or 'thin'."
  )
  # stopped by the iteration-number limit
  record <- .burnin_loop(list(c(mu = 2)), extension_room = 700)$record
  expect_identical(
    BayesTools:::.bt_burnin_warning(record),
    paste0("Chains did not agree at the end of burn-in after 1 extension: R-hat 2.00 for 'mu' (limit 1.50). ",
      "A further extension would number the draws beyond the largest iteration number (2147483647); reduce ",
      "'adapt' or the sampling iterations ('sample' times 'thin') to leave room for it.")
  )
  expect_null(BayesTools:::.bt_burnin_warning(.burnin_loop(list(c(mu = 1.2)))$record))
  expect_null(BayesTools:::.bt_burnin_warning(NULL))
})

test_that("the printed limit never reads above the threshold and a listed value always reads above it", {

  limit <- BayesTools:::.bt_burnin_format_limit
  # two decimals when that reads back as the threshold, otherwise rounded
  # down to three; 1.7 + 0.9 / 6 is the double below 1.85
  automatic <- 1.7 + 0.9 / 2:12
  expect_identical(vapply(automatic, limit, character(1)),
    c("2.15", "2.00", "1.925", "1.88", "1.849", "1.828", "1.812", "1.80", "1.79", "1.781", "1.775"))
  expect_identical(limit(1.5), "1.50")
  expect_identical(limit(1.2345), "1.234")
  # a given value is never rounded up: just below 1.501 prints 1.500
  # (1.5009999999999999 is the double of 1.501)
  expect_identical(limit(1.501 - 1e-13), "1.500")
  expect_identical(limit(1.5009999999999999), "1.501")

  warning_at <- function(threshold, values){
    BayesTools:::.bt_burnin_warning(
      .burnin_loop(list(values), autofit = FALSE, burnin_max_Rhat = threshold)$record)
  }
  warning_for <- function(chains, values) warning_at(1.7 + 0.9 / chains, values)
  printed <- function(threshold, value){
    message <- warning_at(threshold, c(mu = value))
    parts <- regmatches(message, regexec("R-hat ([^ ]+) for 'mu' \\(limit ([^)]+)\\)", message))[[1L]]
    c(value = as.numeric(parts[[2L]]), limit = as.numeric(parts[[3L]]))
  }
  # for the automatic values of 2 to 12 chains and for given values, each
  # value above the threshold (the next doubles above it, relative steps of
  # 1e-15 to 0.1, and steps of 0.0005 up to 0.02) reads above the printed
  # limit, which does not read above the threshold
  next_double <- function(x, steps) x + steps * 2^(floor(log2(x)) - 52)
  for(threshold in c(automatic, 1.5, 1.05, 2.345678, 10, 1.5009999999999999, 1, 1.501 - 1e-13)){
    grid <- c(next_double(threshold, 1:3), threshold * (1 + 10^-(15:1)), threshold + seq(0.0005, 0.02, by = 0.0005))
    expect_true(all(grid > threshold))
    shown <- vapply(grid, function(value) printed(threshold, value), numeric(2))
    expect_true(all(shown["limit", ] <= threshold), info = sprintf("%.17g", threshold))
    expect_true(all(shown["value", ] > shown["limit", ]), info = sprintf("%.17g", threshold))
  }

  remedy <- " Increase 'adapt' and 'burnin'; if the chains mix slowly, increase 'thin'. Or fit with 'autofit = TRUE'."
  # six chains: R-hat in [1.85, 1.855) reads above the limit 1.849
  expect_identical(warning_for(6L, c(mu = 1.854, tau = 1.852)),
    paste0("Chains did not agree at the end of burn-in: R-hat 1.85 for 'mu', 1.85 for 'tau' (limit 1.849).", remedy))
  # where fifteen decimals are not enough: the shortest representation that
  # reads back as the value
  expect_identical(warning_at(2, c(mu = 2 + 2 * .Machine$double.eps)),
    paste0("Chains did not agree at the end of burn-in: R-hat 2.0000000000000004 for 'mu' (limit 2.00).", remedy))
  # seven chains: 1.82857...; a value just above it prints with two decimals
  expect_identical(warning_for(7L, c(mu = 1.8286, tau = 1.82858)),
    paste0("Chains did not agree at the end of burn-in: R-hat 1.83 for 'mu', 1.83 for 'tau' (limit 1.828).", remedy))
  # eight chains: 1.8125, printed 1.812; 1.8126 needs three decimals to exceed it
  expect_identical(warning_for(8L, c(mu = 1.8126, tau = 1.81251)),
    paste0("Chains did not agree at the end of burn-in: R-hat 1.813 for 'mu', 1.813 for 'tau' (limit 1.812).",
      remedy))
  # four chains: 1.925
  expect_identical(warning_for(4L, c(mu = 1.9251)),
    paste0("Chains did not agree at the end of burn-in: R-hat 1.93 for 'mu' (limit 1.925).", remedy))
})

test_that("the readiness settings default when omitted and validate when given", {

  base <- list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10)
  settings <- function(...) BayesTools:::.bt_burnin_check_settings(
    JAGS_check_and_list_autofit_settings(utils::modifyList(base, list(...), keep.null = TRUE)))

  expect_identical(settings(), list(burnin_max_Rhat = "auto", burnin_extend = 500, burnin_max_extend = 10))
  # the validator does not add the elements, so that an omission stays one
  expect_false(any(grepl("^burnin_", names(JAGS_check_and_list_autofit_settings(base)))))
  expect_null(settings(burnin_max_Rhat = NULL)$burnin_max_Rhat)
  expect_identical(settings(burnin_max_Rhat = 1.2, burnin_extend = 100, burnin_max_extend = 0),
    list(burnin_max_Rhat = 1.2, burnin_extend = 100, burnin_max_extend = 0))
  expect_identical(settings(burnin_max_Rhat = "auto")$burnin_max_Rhat, "auto")
  # a named string is stored as the plain string
  expect_identical(settings(burnin_max_Rhat = c(value = "auto"))$burnin_max_Rhat, "auto")
  expect_identical(
    JAGS_check_and_list_autofit_settings(c(base, list(burnin_max_Rhat = NULL)), skip_sample_extend = TRUE)$burnin_max_Rhat,
    NULL
  )

  # "auto" resolves with the number of chains; a number or NULL stays
  threshold <- BayesTools:::.bt_burnin_threshold
  expect_identical(threshold("auto", 2), 1.7 + 0.9 / 2)
  expect_equal(vapply(c(2, 3, 4, 7, 8), function(chains) threshold("auto", chains), numeric(1)),
    c(2.15, 2, 1.925, 1.7 + 0.9 / 7, 1.8125), tolerance = 1e-15)
  expect_identical(threshold(1.3, 4), 1.3)
  expect_null(threshold(NULL, 4))

  invalid <- list(
    list(burnin_max_Rhat = 0.9, "'burnin_max_Rhat' must be equal or higher than 1."),
    list(burnin_max_Rhat = Inf, "The 'burnin_max_Rhat' argument must contain only finite values."),
    list(burnin_max_Rhat = NA_real_, "The 'burnin_max_Rhat' argument cannot contain NA/NaN values."),
    list(burnin_max_Rhat = "Auto", "The 'Auto' values are not recognized by the 'burnin_max_Rhat' argument."),
    list(burnin_max_Rhat = "", "The '' values are not recognized by the 'burnin_max_Rhat' argument."),
    list(burnin_max_Rhat = c("auto", "auto"), "The 'burnin_max_Rhat' argument must have length '1'."),
    list(burnin_max_Rhat = NA_character_, "The 'burnin_max_Rhat' argument cannot contain NA/NaN values."),
    list(burnin_max_Rhat = TRUE, "The 'burnin_max_Rhat' argument must be a numeric vector."),
    # empty values are not NULL, which turns the check off
    list(burnin_max_Rhat = numeric(), "The 'burnin_max_Rhat' argument must have length '1'."),
    list(burnin_max_Rhat = list(), "The 'burnin_max_Rhat' argument must have length '1'."),
    list(burnin_max_Rhat = character(), "The 'burnin_max_Rhat' argument must have length '1'."),
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
  state$events <- character()
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
    state$events <- c(state$events, "sampling")
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
        block = block, prepare = prepare, thin = thin)
      if(!prepare) return(batch(sample, thin))
      session$model <- "compiled"
      if(block > 0L) state$events <- c(state$events, "block")
      if(block > 0L) next_block() else NULL
    },
    .bt_rjags_burnin_block = function(draws, thin, monitor, progress, session){
      state$sweeps <- c(state$sweeps, draws)
      state$block_thin <- c(state$block_thin, thin)
      state$events <- c(state$events, "block")
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

test_that("JAGS_fit() resolves the automatic threshold with the fit's number of chains", {

  # four chains: 1.7 + 0.9 / 4
  set.seed(1)
  state <- .readiness_backend(list(.readiness_chains(chains = 4L)), chains = 4L)
  fit <- .readiness_fit(chains = 4L)
  expect_identical(state$groups[[1L]]$chains, 4L)
  expect_identical(attr(fit, "burnin_check")$burnin_max_Rhat, 1.925)

  # a number given in 'autofit_control' reaches the record
  state <- .readiness_backend(list(.readiness_chains(chains = 4L)), chains = 4L)
  fit <- .readiness_fit(chains = 4L, burnin_max_Rhat = 1.5)
  expect_identical(attr(fit, "burnin_check")$burnin_max_Rhat, 1.5)
})

test_that("the block holds draws at the fit's thinning and the record counts its iterations", {

  # burnin = 2000: min(500, 2000 %/% thin) draws over the last draws * thin
  # iterations of burn-in
  for(case in list(list(thin = 1L, draws = 500L, block = 500L),
                   list(thin = 2L, draws = 500L, block = 1000L),
                   list(thin = 4L, draws = 500L, block = 2000L),
                   list(thin = 5L, draws = 400L, block = 2000L))){
    info <- paste("thin", case$thin)
    set.seed(1)
    state <- .readiness_backend(list(.readiness_chains(n = case$draws)))
    fit <- .readiness_fit(burnin = 2000, thin = case$thin)
    expect_identical(state$groups[[1L]][c("burnin", "block", "prepare", "thin")],
      list(burnin = 2000, block = case$draws, prepare = TRUE, thin = case$thin), info = info)
    record <- attr(fit, "burnin_check")
    expect_identical(record[c("block", "thin", "extension", "status")],
      list(block = case$block, thin = case$thin, extension = 0, status = "ready"), info = info)
    expect_identical(record$blocks$sweeps, case$block, info = info)
    expect_identical(fit$burnin, 2000, info = info)
  }
  # thin = 30: 66 draws are too few, the check does not run
  state <- .readiness_backend(list())
  fit <- .readiness_fit(burnin = 2000, thin = 30L)
  expect_identical(state$groups[[1L]][c("block", "prepare")], list(block = 0L, prepare = FALSE))
  expect_null(attr(fit, "burnin_check"))

  # thin = 101 above burnin = 100: no draw at all; the check does not run
  # and the fit equals the fit with the check off
  set.seed(1)
  state <- .readiness_backend(list())
  fit <- .readiness_fit(burnin = 100, thin = 101L)
  expect_identical(state$groups[[1L]][c("block", "prepare")], list(block = 0L, prepare = FALSE))
  expect_null(attr(fit, "burnin_check"))
  set.seed(1)
  state <- .readiness_backend(list())
  off <- .readiness_fit(burnin = 100, thin = 101L, burnin_max_Rhat = NULL)
  # the parameter map's cache identifier is unique to each fit
  without_id <- function(fit){
    attr(attr(fit, "parameter_map"), "runtime_cache_id") <- NULL
    fit
  }
  expect_identical(without_id(fit), without_id(off))

  # a burn-in that is not a multiple of 'thin': 1000 %/% 3 = 333 draws over
  # the last 999 iterations
  state <- .readiness_backend(list(.readiness_chains(n = 333L)))
  fit <- .readiness_fit(burnin = 1000, thin = 3L)
  expect_identical(state$groups[[1L]]$block, 333L)
  expect_identical(attr(fit, "burnin_check")[c("block", "thin")], list(block = 999L, thin = 3L))
})

test_that("a chain group runs the unmonitored part of burn-in, then the block at the fit's thinning", {

  skip_if_not_installed("rjags")
  run_group <- function(block, thin){
    session <- new.env(parent = emptyenv())
    inits <- lapply(1:2, function(chain) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = chain))
    draws <- BayesTools:::.bt_rjags_group(inits, "model{ mu ~ dnorm(0, 1)\n y ~ dnorm(mu, 1) }",
      data = list(y = 0.5), monitor = "mu", adapt = 100L, burnin = 1000L, sample = 10L, thin = thin,
      silent = TRUE, keep_live = FALSE, prepare = TRUE, block = block, session = session)
    list(draws = draws, iteration = session$model$iter(), state = session$model$state(internal = TRUE))
  }
  # 1000 %/% 3 = 333 draws over the last 999 iterations, after one
  # unmonitored iteration: the chains end where the unmonitored burn-in ends
  off <- run_group(0L, 3L)
  on  <- run_group(333L, 3L)
  expect_null(off$draws)
  expect_identical(vapply(on$draws, nrow, integer(1)), c(333L, 333L))
  expect_identical(coda::thin(on$draws), 3)
  expect_identical(on$iteration, off$iteration)
  expect_identical(on$state, off$state)
})

test_that("the block holds the iterations that JAGS records at the fit's thinning", {

  skip_if_not_installed("rjags")
  # the burn-in block and one extension block of a seeded run at thin = 3,
  # and the same blocks of an unthinned run: JAGS records the first
  # iteration of each window of 'thin' iterations, as in sampling
  run_blocks <- function(block, thin){
    session <- new.env(parent = emptyenv())
    inits <- lapply(1:2, function(chain) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = chain))
    first <- BayesTools:::.bt_rjags_group(inits, "model{ mu ~ dnorm(0, 1)\n y ~ dnorm(mu, 1) }",
      data = list(y = 0.5), monitor = "mu", adapt = 100L, burnin = 1000L, sample = 10L, thin = thin,
      silent = TRUE, keep_live = FALSE, prepare = TRUE, block = block, session = session)
    extension <- BayesTools:::.bt_rjags_burnin_block(600L %/% thin, thin, "mu", "none", session)
    lapply(list(first = first, extension = extension), function(draws){
      lapply(draws, function(chain) as.numeric(chain[, "mu"]))
    })
  }
  unthinned <- run_blocks(1000L, 1L)
  thinned   <- run_blocks(333L, 3L)
  # after 1000 - 333 * 3 = 1 unmonitored iteration, the block holds burn-in
  # iterations 2, 5, ..., 998; the extension its iterations 1, 4, ..., 598
  expect_identical(lengths(thinned$first), c(333L, 333L))
  expect_identical(thinned$first, lapply(unthinned$first, `[`, 2L + 3L * 0:332))
  expect_identical(lengths(thinned$extension), c(200L, 200L))
  expect_identical(thinned$extension, lapply(unthinned$extension, `[`, 1L + 3L * 0:199))
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
    "Chains did not agree at the end of burn-in: R-hat",
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

test_that("the readiness warning follows the kept attempt's sampling and precedes automatic extensions", {

  warned <- function(expr, state){
    withCallingHandlers(expr, BayesTools_burnin_not_converged = function(w){
      state$events <- c(state$events, "warning")
      invokeRestart("muffleWarning")
    })
  }
  set.seed(8)
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  warned(.readiness_fit(), state)
  expect_identical(state$events, c("block", "sampling", "warning"))

  # with autofit: the extensions of burn-in, the sampling, the warning, then
  # the automatic extension after sampling
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  warned(.readiness_fit(autofit = TRUE, burnin_max_extend = 1, min_ESS = 1e6, max_extend = 1), state)
  expect_identical(state$events, c("block", "block", "sampling", "warning", "sampling"))

  # a failed attempt that is restarted warns nothing
  state <- .readiness_backend(list(.readiness_chains(shift = 10)), fail_batches = 1L)
  warned(.readiness_fit(restarts = 2), state)
  expect_identical(state$events, c("block", "sampling", "block", "sampling", "warning"))
})

test_that("automatic extensions and JAGS_extend() never run the check", {

  set.seed(9)
  state <- .readiness_backend(list(.readiness_chains()))
  assessed <- 0L
  block_Rhat <- BayesTools:::.bt_burnin_block_Rhat
  testthat::local_mocked_bindings(
    .bt_burnin_block_Rhat = function(...){
      assessed <<- assessed + 1L
      block_Rhat(...)
    },
    .package = "BayesTools"
  )
  # the convergence check after sampling fails: two automatic extensions
  fit <- .readiness_fit(autofit = TRUE, min_ESS = 1e6, max_extend = 2)
  expect_identical(assessed, 1L)
  expect_identical(state$events, c("block", "sampling", "sampling", "sampling"))
  expect_identical(vapply(state$groups, function(group) group$block, integer(1)), c(500L, 0L, 0L))
  expect_identical(attr(fit, "burnin_check")$status, "ready")

  extended <- JAGS_extend(fit, autofit_control = list(max_Rhat = NULL, min_ESS = 1e6, max_error = NULL,
    max_SD_error = NULL, max_time = NULL, sample_extend = 100, max_extend = 2))
  expect_identical(assessed, 1L)
  expect_identical(state$events, c("block", rep("sampling", 5L)))
  expect_identical(lapply(state$groups[2:5], `[`, c("block", "burnin", "prepare")),
    rep(list(list(block = 0L, burnin = 0L, prepare = FALSE)), 4L))
  expect_identical(NROW(extended$mcmc[[1L]]), 600L)
  expect_identical(attr(extended, "burnin_check"), attr(fit, "burnin_check"))
})

test_that("autofit extends burn-in in blocks and the iteration numbers count them", {

  # thin = 2: blocks of 300 draws over 600 iterations each
  set.seed(3)
  state <- .readiness_backend(list(.readiness_chains(shift = 10, n = 300L), .readiness_chains(shift = 10, n = 300L),
    .readiness_chains(n = 300L)))
  fit <- .readiness_fit(autofit = TRUE, thin = 2L, burnin_extend = 300)
  expect_identical(state$groups[[1L]]$block, 300L)
  expect_identical(state$sweeps, c(300, 300))
  expect_identical(state$block_thin, c(2L, 2L))
  record <- attr(fit, "burnin_check")
  expect_identical(record$status, "ready")
  expect_identical(record[c("block", "burnin_extend", "thin", "extension")],
    list(block = 600L, burnin_extend = 300, thin = 2L, extension = 1200))
  expect_identical(record$blocks$sweeps, c(600L, 600L, 600L))
  expect_identical(fit$burnin, 1800)
  expect_identical(stats::start(fit$mcmc[[1L]]), 100 + 1800 + 1)
  expect_identical(coda::thin(fit$mcmc[[1L]]), 2)
  geometry <- JAGS_draw_geometry(fit)
  expect_identical(geometry$chains$start, rep(1901L, 3L))
  expect_identical(geometry$chains$end, rep(1901L + 199L * 2L, 3L))
  expect_null(attr(fit, "warnings"))

  # 'fit$burnin' keeps the type it was given
  state <- .readiness_backend(list(.readiness_chains(shift = 10), .readiness_chains()))
  fit <- .readiness_fit(autofit = TRUE, burnin = 600L)
  expect_identical(fit$burnin, 1100L)
  state <- .readiness_backend(list(.readiness_chains()))
  expect_identical(.readiness_fit(burnin = 600L)$burnin, 600L)

  # the record's 'extension' keeps the type of 'burnin_extend' at the
  # default thin = 1
  state <- .readiness_backend(list(.readiness_chains(shift = 10, n = 300L), .readiness_chains(n = 300L)))
  fit <- .readiness_fit(autofit = TRUE, burnin_extend = 300L, thin = 1)
  expect_identical(attr(fit, "burnin_check")[c("block", "extension")], list(block = 300L, extension = 300L))

  # the cap: the warning names the extensions
  state <- .readiness_backend(list(.readiness_chains(shift = 10)))
  expect_warning(
    fit <- .readiness_fit(autofit = TRUE, burnin_max_extend = 2),
    "after 2 extensions: R-hat",
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
    "after 1 extension: R-hat",
    fixed = TRUE
  )
  expect_identical(attr(fit, "burnin_check")$stop_reason, "max_time")
  expect_identical(state$sweeps, 500)
})

test_that("no extension numbers the draws beyond the integer range", {

  limit <- .Machine$integer.max

  # a burn-in beyond the iteration numbers is refused before any backend work
  # (its block of 2^30 draws at thin = 2 would overflow an integer product)
  state <- .readiness_backend(list())
  expect_error(.readiness_fit(burnin = 2147483648, thin = 2L, burnin_extend = 1073741824L),
    "The 'burnin' must be equal or lower than 2147483647.", fixed = TRUE)
  expect_length(state$groups, 0L)

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

  # with thin = 2 a block of 500 draws runs 1000 iterations: the room after
  # the last draw (adapt + 600 + 1 + 199 * 2) is 2000 iterations, so two
  # blocks fit and the last draw reaches the limit exactly
  state <- .readiness_backend(list(.readiness_chains(shift = 10, n = 300L), .readiness_chains(shift = 10)))
  expect_warning(
    fit <- .readiness_fit(autofit = TRUE, adapt = limit - 2999, burnin = 600, thin = 2L),
    "A further extension would number the draws beyond the largest iteration number",
    class = "BayesTools_burnin_not_converged"
  )
  record <- attr(fit, "burnin_check")
  expect_identical(record[c("stop_reason", "extension")], list(stop_reason = "iteration_limit", extension = 2000))
  expect_identical(record$blocks$sweeps, c(600L, 1000L, 1000L))
  expect_identical(state$sweeps, c(500, 500))
  expect_identical(fit$burnin, 2600)
  expect_identical(JAGS_draw_geometry(fit)$chains$end, rep(limit, 3L))

  # sampling: the first automatic extension reaches the limit exactly, the
  # second would pass it and is not started; the fit keeps its draws
  message <- paste0("The automatic model fitting was terminated due to the iteration limit: an extension of ",
    "200 draws per chain would number the draws beyond the largest iteration number (2147483647). ",
    "Use a smaller 'adapt', 'burnin', or 'sample' times 'thin' to leave room for extensions.")
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
  state <- .readiness_backend(list(.readiness_chains(shift = 10, chains = 4L), .readiness_chains(chains = 4L)),
    chains = 4L, fail_batches = 1L)
  # each attempt checks with the threshold of the fit's four chains
  thresholds <- numeric()
  readiness <- BayesTools:::.bt_burnin_readiness
  testthat::local_mocked_bindings(
    .bt_burnin_readiness = function(..., burnin_max_Rhat){
      thresholds <<- c(thresholds, burnin_max_Rhat)
      readiness(..., burnin_max_Rhat = burnin_max_Rhat)
    },
    .package = "BayesTools"
  )
  expect_no_warning(fit <- .readiness_fit(restarts = 2, chains = 4L))
  expect_length(state$groups, 2L)
  expect_identical(thresholds, c(1.925, 1.925))
  expect_identical(attr(fit, "burnin_check")$burnin_max_Rhat, 1.925)
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

  # an error that a worker returns is not a transport failure, whatever its
  # wording: it keeps its class and the models are released
  testthat::local_mocked_bindings(
    .bt_rjags_batch = function(...) stop("error reading from connection"),
    .package = "BayesTools"
  )
  failed <- tryCatch(BayesTools:::.bt_rjags_run("", list(), list(1L, 2L, 3L), "mu", 100L, 300L, 100L, 1L,
    cl = list(1L, 2L), silent = TRUE, readiness = readiness), error = identity)
  expect_false(inherits(failed, "BayesTools_JAGS_worker_connection_error"))
  expect_identical(conditionMessage(failed), "error reading from connection")
  expect_identical(released, 2L)
})

test_that("a staged parallel run releases the workers' models on every exit while they can be reached", {

  live <- BayesTools:::.bt_rjags_worker_session
  on.exit(live$model <- NULL, add = TRUE)
  released <- calls <- 0L
  failing_call <- interrupted_call <- 0L
  cancel <- FALSE
  interrupt <- function(...){
    signalCondition(structure(class = c("interrupt", "condition"), list(message = "interrupted", call = NULL)))
  }
  testthat::local_mocked_bindings(
    clusterApply = function(cl, x, fun, ...){
      calls <<- calls + 1L
      if(calls == failing_call) stop("error reading from connection")
      if(calls == interrupted_call) interrupt()
      lapply(x, fun, ...)
    },
    clusterCall = function(cl, fun, ...){
      released <<- released + 1L
      lapply(cl, function(worker) fun(...))
    },
    .package = "parallel"
  )
  testthat::local_mocked_bindings(
    .bt_rjags_group = function(inits, ...){
      live$model <- "compiled"
      .readiness_chains(chains = length(inits), n = 200L)
    },
    .bt_rjags_batch = function(...){
      list(mcmc = .readiness_chains(chains = 1L, n = 10L), end_state = list(), samplers = data.frame())
    },
    .JASP_progress_bar_start = function(...) NULL,
    .JASP_progress_bar_tick = function() if(cancel) stop("JASP cancelled", call. = FALSE),
    .package = "BayesTools"
  )
  run <- function(assess, jasp_prefix = NULL){
    live$model <- NULL
    released <<- calls <<- 0L
    readiness <- list(block = 200L, burnin_max_Rhat = 1.5, burnin_extend = 500, max_extensions = 0,
      autofit = FALSE, time_exceeded = function() FALSE, assess = assess)
    tryCatch(withRestarts(
      BayesTools:::.bt_rjags_run("", list(), list(1L, 2L, 3L), "mu", 100L, 300L, 100L, 1L,
        cl = list(1L, 2L), silent = TRUE, jasp_prefix = jasp_prefix, readiness = readiness),
      BayesTools_jasp_progress_error = function(e, ran) stop(e)),
      error = identity, interrupt = identity)
  }

  # an interrupt on the coordinator between worker calls
  result <- run(interrupt)
  expect_s3_class(result, "interrupt")
  expect_identical(released, 1L)
  expect_null(live$model)

  # a JASP cancellation after the check
  cancel <- TRUE
  result <- run(function(draws) c(mu = 1), jasp_prefix = "")
  cancel <- FALSE
  expect_identical(conditionMessage(result), "JASP cancelled")
  expect_identical(released, 1L)
  expect_null(live$model)

  # a coordinator error that names a parameter 'connection[2]' keeps its class
  selector <- errorCondition(
    "The requested convergence monitor 'connection[2]' is not available in the fitted model.",
    class = "BayesTools_parameter_not_found")
  result <- run(function(draws) stop(selector))
  expect_identical(result, selector)
  expect_identical(BayesTools:::.JAGS_run_backend(function() stop(result), list()), selector)
  expect_identical(released, 1L)
  expect_null(live$model)

  # the worker call itself fails (the transport), or is interrupted while in
  # flight (its reply is pending): the workers are not contacted again; the
  # owner of the cluster stops it
  failing_call <- 2L
  result <- run(function(draws) c(mu = 1))
  expect_s3_class(result, "BayesTools_JAGS_worker_connection_error")
  expect_identical(conditionMessage(result$parent), "error reading from connection")
  expect_identical(released, 0L)
  failing_call <- 0L
  interrupted_call <- 2L
  result <- run(function(draws) c(mu = 1))
  expect_s3_class(result, "interrupt")
  expect_identical(released, 0L)
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
