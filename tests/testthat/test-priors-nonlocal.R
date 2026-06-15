skip_if_not_test_profile("unit")

moment_reference_lpdf <- function(x, location, tau, order){
  delta <- x - location
  out <- rep(-Inf, length(delta))
  keep <- is.finite(delta) & delta != 0
  log_double_factorial <- lgamma(2 * order + 1) - order * log(2) - lgamma(order + 1)
  out[keep] <- 2 * order * log(abs(delta[keep])) +
    stats::dnorm(delta[keep], 0, sqrt(tau), log = TRUE) -
    order * log(tau) -
    log_double_factorial
  out
}

invmoment_reference_lpdf <- function(x, location, tau, order, df){
  delta <- x - location
  out <- rep(-Inf, length(delta))
  keep <- is.finite(delta) & delta != 0
  out[keep] <- log(order) +
    (df / 2) * log(tau) -
    lgamma(df / (2 * order)) -
    (df + 1) * log(abs(delta[keep])) -
    (tau / delta[keep]^2)^order
  out
}

test_that("moment prior constructor stores canonical parameters", {
  p_mode <- prior("moment", list(mode = -.5, location = .25))
  expect_equal(p_mode$distribution, "moment")
  expect_equal(p_mode$parameters$mode, .5)
  expect_equal(p_mode$parameters$tau, .125)
  expect_equal(p_mode$parameters$order, 1L)
  expect_equal(p_mode$parameters$location, .25)

  p_tau <- prior("pmom", list(tau = .125, location = .25))
  expect_equal(p_tau$parameters, p_mode$parameters)

  p_default_order_location <- prior("moment", list(.5))
  expect_equal(p_default_order_location$parameters$location, 0)
  expect_equal(p_default_order_location$parameters$order, 1L)
  expect_equal(p_default_order_location$parameters$mode, .5)
  expect_equal(p_default_order_location$parameters$tau, .125)

  p_tau_default_order_location <- prior("moment", list(tau = .125))
  expect_equal(p_tau_default_order_location$parameters, p_default_order_location$parameters)

  p_positional_alias <- prior("pmom", list(.125))
  expect_equal(p_positional_alias$parameters$mode, .125)
  expect_equal(p_positional_alias$parameters$tau, .125^2 / 2)

  p_near_order <- prior("moment", list(mode = .5, order = 2 - 1e-9))
  expect_equal(p_near_order$parameters$order, 2L)
  expect_equal(p_near_order$parameters$tau, .5^2 / 4)
})

test_that("inverse-moment prior constructor stores canonical parameters", {
  p_mode <- prior("invmoment", list(mode = -.5, df = 3, location = .25))
  expect_equal(p_mode$distribution, "invmoment")
  expect_equal(p_mode$parameters$mode, .5)
  expect_equal(p_mode$parameters$tau, .5)
  expect_equal(p_mode$parameters$order, 1L)
  expect_equal(p_mode$parameters$df, 3)
  expect_equal(p_mode$parameters$location, .25)

  p_tau <- prior("pimom", list(tau = .5, nu = 3, location = .25))
  expect_equal(p_tau$parameters, p_mode$parameters)

  p_default_order_location <- prior("invmoment", list(.5, 3))
  expect_equal(p_default_order_location$distribution, "invmoment")
  expect_equal(p_default_order_location$parameters$location, 0)
  expect_equal(p_default_order_location$parameters$order, 1L)
  expect_equal(p_default_order_location$parameters$mode, .5)
  expect_equal(p_default_order_location$parameters$tau, .5)

  p_tau_default_order_location <- prior("pimom", list(tau = .5, nu = 3))
  expect_equal(p_tau_default_order_location$parameters, p_default_order_location$parameters)

  p_near_order <- prior("invmoment", list(mode = .5, order = 2 - 1e-9, df = 3))
  expect_equal(p_near_order$parameters$order, 2L)
  expect_equal(p_near_order$parameters$tau, .25)
})

test_that("nonlocal priors reject invalid parameterizations", {
  expect_error(prior("moment", list(mode = .5, tau = .1)), "exactly one")
  expect_error(prior("moment", list(order = 1)), "exactly one")
  expect_error(prior("moment", list(mode = 0)), "mode")
  expect_error(prior("moment", list(.5, 1)), "single 'mode'")
  expect_error(prior("moment", list(mode = .5, order = 1.5)), "order")
  expect_error(prior("moment", list(mode = .5, location = Inf)), "location")
  expect_error(prior("moment", list(m = .5)), "not supported")
  expect_error(prior("moment", list(mode = .Machine$double.xmin)), "mode")
  expect_error(prior("moment", list(mode = .Machine$double.xmax)), "mode")
  expect_error(prior("moment", list(tau = .Machine$double.xmax)), "tau")
  expect_error(
    prior("moment", list(mode = .5, order = .Machine$integer.max + 1)),
    "order"
  )

  expect_error(prior("invmoment", list(mode = .5)), "df")
  expect_error(prior("invmoment", list(mode = .5, order = 1)), "df")
  expect_error(prior("invmoment", list(.5, 1, 3)), "requires 'mode' and 'df'")
  expect_error(prior("invmoment", list(mode = .5, df = 0)), "df")
  expect_error(prior("invmoment", list(tau = 0, df = 3)), "tau")
  expect_error(prior("invmoment", list(mode = .5, df = 3, nu = 3)), "only one")
  expect_error(prior("invmoment", list(mode = .5, d = 3)), "not supported")
  expect_error(prior("invmoment", list(mode = .Machine$double.xmin, df = 3)), "mode")
  expect_error(prior("invmoment", list(mode = .Machine$double.xmax, df = 3)), "mode")
  expect_error(
    prior("invmoment", list(mode = .5, order = .Machine$integer.max + 1, df = 3)),
    "order"
  )
})

test_that("moment prior density, distribution, and quantiles match reference identities", {
  p <- prior("moment", list(mode = .5, location = .25))
  x <- c(-.75, -.25, .25, .75, 1.25)

  expect_equal(lpdf(p, x), moment_reference_lpdf(x, .25, .125, 1), tolerance = 1e-12)
  expect_equal(pdf(p, .25), 0)
  expect_equal(lpdf(p, .25), -Inf)
  expect_equal(pdf(p, .25 - .5), pdf(p, .25 + .5), tolerance = 1e-12)

  probs <- c(.05, .25, .5, .75, .95)
  q <- quant(p, probs)
  expect_equal(q[3], .25)
  expect_equal(cdf(p, q), probs, tolerance = 1e-12)
  expect_equal(ccdf(p, q), 1 - probs, tolerance = 1e-12)
  expect_true(q[1] < q[2] && q[2] < q[3] && q[3] < q[4] && q[4] < q[5])

  probs_no_median <- c(.05, .25, .75, .95)
  q_no_median <- quant(p, probs_no_median)
  expect_equal(
    BayesTools:::.pmoment_prior(q_no_median, .25, .125, 1, lower.tail = FALSE),
    1 - probs_no_median,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.pmoment_prior(q_no_median, .25, .125, 1, log.p = TRUE),
    log(probs_no_median),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.pmoment_prior(q_no_median, .25, .125, 1, lower.tail = FALSE, log.p = TRUE),
    log1p(-probs_no_median),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.qmoment_prior(log(probs_no_median), .25, .125, 1, log.p = TRUE),
    q_no_median,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.qmoment_prior(log1p(-probs_no_median), .25, .125, 1, lower.tail = FALSE, log.p = TRUE),
    q_no_median,
    tolerance = 1e-12
  )

  far_q <- 10
  expect_gt(ccdf(p, far_q), 0)
  expect_true(is.finite(BayesTools:::.qmoment_prior(1e-170, .25, .125, 1, lower.tail = FALSE)))
  expect_true(is.finite(BayesTools:::.qmoment_prior(log(1e-300), .25, .125, 1, lower.tail = FALSE, log.p = TRUE)))

  expect_equal(
    stats::integrate(function(z) pdf(p, z), lower = -Inf, upper = Inf)$value,
    1,
    tolerance = 1e-8
  )
})

test_that("inverse-moment prior density, distribution, and quantiles match reference identities", {
  p <- prior("invmoment", list(mode = .5, df = 3, location = .25))
  x <- c(-1.25, -.25, .25, .75, 1.75)

  expect_equal(lpdf(p, x), invmoment_reference_lpdf(x, .25, .5, 1, 3), tolerance = 1e-12)
  expect_equal(pdf(p, .25), 0)
  expect_equal(lpdf(p, .25), -Inf)
  expect_equal(pdf(p, .25 - .5), pdf(p, .25 + .5), tolerance = 1e-12)

  probs <- c(.05, .25, .5, .75, .95)
  q <- quant(p, probs)
  expect_equal(q[3], .25)
  expect_equal(cdf(p, q), probs, tolerance = 1e-12)
  expect_equal(ccdf(p, q), 1 - probs, tolerance = 1e-12)
  expect_true(q[1] < q[2] && q[2] < q[3] && q[3] < q[4] && q[4] < q[5])

  probs_no_median <- c(.05, .25, .75, .95)
  q_no_median <- quant(p, probs_no_median)
  expect_equal(
    BayesTools:::.pinvmoment_prior(q_no_median, .25, .5, 1, 3, lower.tail = FALSE),
    1 - probs_no_median,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.pinvmoment_prior(q_no_median, .25, .5, 1, 3, log.p = TRUE),
    log(probs_no_median),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.pinvmoment_prior(q_no_median, .25, .5, 1, 3, lower.tail = FALSE, log.p = TRUE),
    log1p(-probs_no_median),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.qinvmoment_prior(log(probs_no_median), .25, .5, 1, 3, log.p = TRUE),
    q_no_median,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.qinvmoment_prior(log1p(-probs_no_median), .25, .5, 1, 3, lower.tail = FALSE, log.p = TRUE),
    q_no_median,
    tolerance = 1e-12
  )

  far_q <- 100
  expect_gt(ccdf(p, far_q), 0)
  expect_true(is.finite(BayesTools:::.qinvmoment_prior(1e-12, .25, .5, 1, 3, lower.tail = FALSE)))
  expect_true(is.finite(BayesTools:::.qinvmoment_prior(log(1e-20), .25, .5, 1, 3, lower.tail = FALSE, log.p = TRUE)))

  expect_equal(
    stats::integrate(function(z) pdf(p, z), lower = -Inf, upper = .25)$value +
      stats::integrate(function(z) pdf(p, z), lower = .25, upper = Inf)$value,
    1,
    tolerance = 1e-8
  )
})

test_that("nonlocal priors report analytic moments and finite random draws", {
  p_moment <- prior("moment", list(mode = .5, location = .25))
  p_invmoment <- prior("invmoment", list(mode = .5, df = 3, location = .25))

  expect_equal(mean(p_moment), .25)
  expect_equal(var(p_moment), .125 * 3)
  expect_equal(mean(p_invmoment), .25)
  expect_equal(var(p_invmoment), 1, tolerance = 1e-12)
  expect_true(is.nan(mean(prior("invmoment", list(mode = .5, df = 1)))))
  expect_true(is.nan(var(prior("invmoment", list(mode = .5, df = 2)))))

  set.seed(1)
  draws <- rng(p_moment, 100)
  expect_length(draws, 100)
  expect_true(all(is.finite(draws)))
  expect_false(any(draws == p_moment$parameters$location))

  set.seed(2)
  draws <- rng(p_invmoment, 100)
  expect_length(draws, 100)
  expect_true(all(is.finite(draws)))
  expect_false(any(draws == p_invmoment$parameters$location))
})

test_that("nonlocal priors normalize truncation", {
  priors <- list(
    prior("moment", list(mode = .5, location = .25), truncation = list(lower = .25, upper = Inf)),
    prior("invmoment", list(mode = .5, df = 3, location = .25), truncation = list(lower = .25, upper = Inf))
  )
  probs <- c(.1, .5, .9)

  for(p in priors){
    expect_equal(cdf(p, .25), 0, tolerance = 1e-12)
    expect_equal(ccdf(p, .25), 1, tolerance = 1e-12)

    q <- quant(p, probs)
    expect_true(all(q >= .25))
    expect_equal(cdf(p, q), probs, tolerance = 1e-12)

    expect_equal(
      stats::integrate(function(z) pdf(p, z), lower = .25, upper = Inf)$value,
      1,
      tolerance = 1e-8
    )

    set.seed(3)
    draws <- rng(p, 100)
    expect_true(all(is.finite(draws)))
    expect_true(all(draws >= .25))
  }

  tail_priors <- list(
    prior("moment", list(mode = .5), truncation = list(lower = 10, upper = Inf)),
    prior("invmoment", list(mode = .5, df = 3), truncation = list(lower = 100, upper = Inf))
  )

  for(p in tail_priors){
    expect_gt(BayesTools:::.prior_C(p), 0)

    q <- quant(p, probs)
    expect_true(all(is.finite(q)))
    expect_true(all(q >= p$truncation$lower))
    expect_equal(cdf(p, q), probs, tolerance = 1e-10)
    expect_equal(ccdf(p, q), 1 - probs, tolerance = 1e-10)
    expect_true(is.finite(pdf(p, q[2])))

    set.seed(4)
    draws <- rng(p, 100)
    expect_true(all(is.finite(draws)))
    expect_true(all(draws >= p$truncation$lower))
  }
})

test_that("nonlocal priors print user-facing parameters", {
  p_moment_default <- prior("moment", list(mode = .5))
  p_invmoment_default <- prior("invmoment", list(mode = .5, df = 3))
  p_moment_order <- prior("moment", list(mode = .5, order = 2))
  p_invmoment_order <- prior("invmoment", list(mode = .5, order = 2, df = 3))
  p_moment <- prior("moment", list(mode = .5, location = .25))
  p_invmoment <- prior("invmoment", list(mode = .5, df = 3, location = .25))

  expect_equal(utils::capture.output(print(p_moment_default)), "Moment(0.5)")
  expect_equal(utils::capture.output(print(p_moment_default, short_name = TRUE)), "MOM(0.5)")
  expect_equal(
    utils::capture.output(print(p_moment_default, parameter_names = TRUE)),
    "Moment(mode = 0.5)"
  )

  expect_equal(utils::capture.output(print(p_moment_order)), "Moment(mode = 0.5, order = 2)")
  expect_equal(
    utils::capture.output(print(p_moment_order, parameter_names = TRUE)),
    "Moment(mode = 0.5, order = 2)"
  )

  expect_equal(utils::capture.output(print(p_moment)), "Moment(mode = 0.5, location = 0.25)")
  expect_equal(utils::capture.output(print(p_moment, short_name = TRUE)), "MOM(mode = 0.5, location = 0.25)")
  expect_equal(
    utils::capture.output(print(p_moment, parameter_names = TRUE)),
    "Moment(mode = 0.5, location = 0.25)"
  )
  expect_false(grepl("tau", print(p_moment, parameter_names = TRUE, silent = TRUE), fixed = TRUE))

  expect_equal(utils::capture.output(print(p_invmoment_default)), "InvMoment(0.5, 3)")
  expect_equal(utils::capture.output(print(p_invmoment_default, short_name = TRUE)), "iMOM(0.5, 3)")
  expect_equal(
    utils::capture.output(print(p_invmoment_default, parameter_names = TRUE)),
    "InvMoment(mode = 0.5, df = 3)"
  )

  expect_equal(utils::capture.output(print(p_invmoment_order)), "InvMoment(mode = 0.5, order = 2, df = 3)")
  expect_equal(
    utils::capture.output(print(p_invmoment_order, parameter_names = TRUE)),
    "InvMoment(mode = 0.5, order = 2, df = 3)"
  )

  expect_equal(utils::capture.output(print(p_invmoment)), "InvMoment(mode = 0.5, df = 3, location = 0.25)")
  expect_equal(utils::capture.output(print(p_invmoment, short_name = TRUE)), "iMOM(mode = 0.5, df = 3, location = 0.25)")
  expect_equal(
    utils::capture.output(print(p_invmoment, parameter_names = TRUE)),
    "InvMoment(mode = 0.5, df = 3, location = 0.25)"
  )
  expect_false(grepl("tau", print(p_invmoment, parameter_names = TRUE, silent = TRUE), fixed = TRUE))
})

test_that("nonlocal prior density ranges are finite and nonnegative", {
  priors <- list(
    prior("moment", list(mode = .5, location = .25)),
    prior("invmoment", list(mode = .5, df = 3, location = .25))
  )

  for(p in priors){
    p_range <- range(p)
    expect_true(all(is.finite(p_range)))
    expect_true(p_range[1] < p_range[2])

    p_density <- density(p, n_points = 100)
    expect_true(all(is.finite(p_density$x)))
    expect_true(all(is.finite(p_density$y)))
    expect_true(all(p_density$y >= 0))
  }
})

test_that("nonlocal priors generate direct JAGS syntax and standard metadata", {
  p_moment <- prior("moment", list(mode = .5))
  p_invmoment <- prior("invmoment", list(mode = .5, df = 3))

  syntax_moment <- JAGS_add_priors("model{}", list(theta = p_moment))
  expect_match(syntax_moment, "theta ~ dbt_moment(0,0.125,1)", fixed = TRUE)

  syntax_invmoment <- JAGS_add_priors("model{}", list(theta = p_invmoment))
  expect_match(syntax_invmoment, "theta ~ dbt_invmoment(0,0.5,1,3)", fixed = TRUE)

  syntax_truncated <- JAGS_add_priors(
    "model{}",
    list(theta = prior("moment", list(mode = .5), truncation = list(lower = 0, upper = Inf)))
  )
  expect_match(syntax_truncated, "theta ~ dbt_moment(0,0.125,1)T(0,)", fixed = TRUE)

  expect_equal(JAGS_to_monitor(list(theta = p_moment)), "theta")
  inits <- JAGS_get_inits(list(theta = p_moment), chains = 1, seed = 1)[[1]]
  expect_named(inits, c("theta", ".RNG.seed", ".RNG.name"), ignore.order = TRUE)
  expect_true(is.finite(inits$theta))

  expect_equal(
    JAGS_marglik_priors(c(theta = .5), list(theta = p_moment)),
    lpdf(p_moment, .5),
    tolerance = 1e-12
  )
})

test_that("nonlocal priors request the BayesTools JAGS module recursively", {
  p_moment <- prior("moment", list(mode = .5))
  p_inclusion <- prior("beta", list(1, 1))
  p_spike_and_slab <- prior_spike_and_slab(p_moment, p_inclusion)
  p_mixture <- prior_mixture(list(prior("normal", list(0, 1)), p_moment))

  expect_true(BayesTools:::.JAGS_prior_list_uses_nonlocal(list(theta = p_moment)))
  expect_true(BayesTools:::.JAGS_prior_list_uses_nonlocal(list(theta = p_spike_and_slab)))
  expect_true(BayesTools:::.JAGS_prior_list_uses_nonlocal(list(theta = p_mixture)))
  expect_false(BayesTools:::.JAGS_prior_list_uses_nonlocal(list(theta = prior("normal", list(0, 1)))))
})
