# TEST FILE: Weightfunction prior redesign
# ============================================================================ #
#
# PURPOSE:
#   Contract tests for the unified prior_weightfunction() API, independent
#   weight priors, one/two-sided mapping, JAGS code generation, and bridge
#   sampling parameter extraction.
#
# TAGS: @priors, @jags, @weightfunctions
# ============================================================================ #

test_that("prior_weightfunction stores canonical geometry and weight priors", {

  wf <- prior_weightfunction(
    side = "one-sided",
    steps = c(.025, .05),
    weights = wf_cumulative(c(1, 2, 3)),
    prior_weights = 2
  )

  expect_true(is.prior.weightfunction(wf))
  expect_equal(wf$distribution, "weightfunction")
  expect_equal(wf$side, "one-sided")
  expect_equal(wf$steps, c(.025, .05))
  expect_equal(wf$bins$lower, c(0, .025, .05))
  expect_equal(wf$bins$upper, c(.025, .05, 1))
  expect_equal(wf$bins$reference, c(TRUE, FALSE, FALSE))
  expect_equal(wf$weights$type, "cumulative")
  expect_equal(wf$weights$alpha, c(1, 2, 3))
  expect_equal(wf$prior_weights, 2)
})

test_that("weightfunction constructors validate independent scales", {

  expect_silent(wf_independent(prior("beta", list(1, 1))))
  expect_silent(wf_independent(prior("gamma", list(2, 1))))
  expect_error(
    wf_independent(prior("normal", list(0, 1))),
    "non-negative support"
  )

  expect_silent(wf_independent(
    prior("normal", list(0, 1)),
    scale = "log_omega"
  ))

  expect_error(
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(.9, .5))),
    "reference-bin"
  )
})

test_that("weightfunctions_mapping expands two-sided priors onto one-sided cuts", {

  one_sided <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_cumulative(c(1, 1, 1))
  )
  two_sided <- prior_weightfunction(
    "two-sided",
    c(.05),
    wf_fixed(c(1, .5))
  )

  expect_equal(
    weightfunctions_mapping(list(one_sided, two_sided), cuts_only = TRUE, one_sided = TRUE),
    c(0, .025, .05, .975, 1)
  )
  expect_equal(
    weightfunctions_mapping(list(one_sided, two_sided), one_sided = TRUE),
    list(c(1L, 2L, 3L, 3L), c(1L, 2L, 2L, 1L))
  )
})

test_that("JAGS generation uses component-local omega for bias mixtures", {

  fixed <- prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .5)))
  independent <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_independent(prior("beta", list(2, 3)))
  )
  bias <- prior_mixture(list(prior_none(), fixed, independent))

  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_match(syntax, "omega_component_1\\[1\\] <- 1")
  expect_match(syntax, "omega_local_component_2\\[2\\] <- 0.5")
  expect_match(syntax, "omega_local_component_3\\[2\\] ~ dbeta\\(2,3\\)")
  expect_match(syntax, "omega\\[1\\] <- omega_component_1\\[1\\] \\* equals\\(bias_indicator, 1\\)")
  expect_false(grepl("eta2omega", syntax, fixed = TRUE))
})

test_that("JAGS bridge helpers use natural latent weight parameters", {

  cumulative <- prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)))
  independent <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_independent(prior("beta", list(2, 3)))
  )
  log_independent <- prior_weightfunction(
    "one-sided",
    c(.05),
    wf_independent(
      prior("normal", list(0, 1)),
      "log_omega"
    )
  )

  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(cumulative)), paste0("eta[", 1:3, "]"))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(independent)), paste0("omega[", 2:3, "]"))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(log_independent)), "log_omega[2]")
  expect_equal(attr(.JAGS_bridgesampling_posterior_info.weightfunction(log_independent), "ub"), c("log_omega[2]" = Inf))

  samples <- c("eta[1]" = 1, "eta[2]" = 2, "eta[3]" = 3, "omega[2]" = 1.2, "omega[3]" = .1, "log_omega[2]" = .5)
  expect_equal(JAGS_marglik_parameters(samples, list(omega = cumulative))$omega, c(1, 5/6, 1/2))
  expect_equal(JAGS_marglik_parameters(samples, list(omega = independent))$omega, c(1, 1.2, .1))
  expect_equal(JAGS_marglik_parameters(samples, list(omega = log_independent))$omega, c(1, exp(.5)))
  expect_equal(
    JAGS_marglik_priors(samples, list(omega = log_independent)),
    mlpdf(log_independent$weights$prior, .5),
    tolerance = 1e-12
  )
})

test_that("JAGS syntax and fitting allow independent omega weights above one", {

  skip_if_not_installed("rjags")

  omega_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega")
  )

  omega_syntax <- JAGS_add_priors("model{}", list(omega = omega_prior))
  log_syntax   <- JAGS_add_priors("model{}", list(omega = log_prior))

  expect_match(omega_syntax, "omega\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(log_syntax, "log_omega\\[2\\] ~ dnorm")
  expect_match(log_syntax, "omega\\[2\\] <- exp\\(log_omega\\[2\\]\\)")

  omega_fit <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(omega = omega_prior),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 300,
    seed       = 11
  ))
  log_fit <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(omega = log_prior),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 300,
    seed       = 12
  ))

  omega_samples <- as.matrix(.fit_to_posterior(omega_fit))
  log_samples   <- as.matrix(.fit_to_posterior(log_fit))

  expect_true("omega[2]" %in% colnames(omega_samples))
  expect_true("omega[2]" %in% colnames(log_samples))
  expect_true("log_omega[2]" %in% colnames(log_samples))
  expect_gt(mean(omega_samples[, "omega[2]"] > 1), .90)
  expect_gt(mean(log_samples[, "omega[2]"] > 1), .90)
  expect_equal(
    unname(log_samples[, "omega[2]"]),
    unname(exp(log_samples[, "log_omega[2]"])),
    tolerance = 1e-8
  )
})

test_that("point(1) weightfunction null components are handled explicitly", {

  wf <- prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1)), prior_weights = 3)
  point_null <- prior("point", list(1), prior_weights = 1)

  mixed <- .mix_priors.weightfunction(
    list(point_null, wf),
    parameter = "omega",
    seed = 1,
    n_samples = 40
  )

  expect_equal(colnames(mixed), c("omega[0,0.05]", "omega[0.05,1]"))
  expect_equal(unname(mixed[attr(mixed, "models_ind") == 1, ]), matrix(1, nrow = 10, ncol = 2))
  expect_error(
    .mix_priors.weightfunction(
      list(prior("point", list(.5), prior_weights = 1), wf),
      parameter = "omega",
      n_samples = 40
    ),
    "point\\(1\\)/none null priors"
  )
})

test_that("omega diagnostics reject bias mixtures without weightfunctions", {

  fit <- structure(list(), class = c("runjags", "BayesTools_fit"))
  attr(fit, "prior_list") <- list(
    bias = prior_mixture(list(
      prior_PET("normal", list(0, 1), prior_weights = 1),
      prior_PEESE("normal", list(0, 1), prior_weights = 1)
    ))
  )

  expect_error(
    JAGS_diagnostics_density(fit, parameter = "omega"),
    "at least one weightfunction component"
  )
})
