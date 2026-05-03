skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS Marginal Likelihood Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for JAGS marginal likelihood computation functions.
#   Uses simple models where the log marginal likelihood is known to be 0
#   (for prior samples, the marginal likelihood for any proper prior is 1).
#
# DEPENDENCIES:
#   - rjags: For JAGS model fitting
#   - bridgesampling: For marginal likelihood computation
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags")
#   - Note: Creates fresh models, does not need pre-fitted models
#
# MODELS/FIXTURES:
#   - Creates models with known analytical marginal likelihoods for validation
#
# TAGS: @evaluation, @JAGS, @marginal-likelihood
# ============================================================================ #

# Load common test helpers
source(testthat::test_path("common-functions.R"))

test_that("direct composed bias priors support bridge-sampling helpers", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))
  phacking  <- prior_phacking(form = "linear", alpha = prior("beta", list(2, 3)))
  bias      <- prior_bias(selection = selection, phacking = phacking)

  posterior <- matrix(
    c(
      1.5, 2.5, .2,
      1.1, 2.1, .4
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(posterior) <- c("eta[1]", "eta[2]", "alpha")

  prepared <- JAGS_bridgesampling_posterior(posterior, list(bias = bias))
  expect_equal(colnames(prepared), c("eta[1]", "eta[2]", "alpha"))
  expect_equal(attr(prepared, "lb"), c("eta[1]" = 0, "eta[2]" = 0, "alpha" = 0))
  expect_equal(attr(prepared, "ub"), c("eta[1]" = Inf, "eta[2]" = Inf, "alpha" = 1))

  samples <- posterior[1, ]
  expected_prior_density <-
    sum(stats::dgamma(samples[c("eta[1]", "eta[2]")], shape = c(1, 2), rate = 1, log = TRUE)) +
    stats::dbeta(samples[["alpha"]], 2, 3, log = TRUE)
  expect_equal(
    JAGS_marglik_priors(samples, list(bias = bias)),
    expected_prior_density,
    tolerance = 1e-12
  )

  parameters <- JAGS_marglik_parameters(samples, list(bias = bias))
  constants <- phack_backend_constants(phacking$form, phacking$source, phacking$destination, target = phacking$target)
  expect_equal(parameters$omega, c(1, samples[["eta[2]"]] / sum(samples[c("eta[1]", "eta[2]")])))
  expect_equal(parameters$alpha, samples[["alpha"]])
  expect_equal(
    parameters$pi_null,
    samples[["alpha"]] * constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    parameters$beta_null,
    samples[["alpha"]] * constants$beta_null_per_alpha,
    tolerance = 1e-12
  )
})

test_that("bias mixtures fail explicitly in bridge-sampling helpers", {

  bias <- prior_mixture(list(
    prior_none(),
    prior_phacking(form = "linear")
  ))
  samples <- c("bias_indicator" = 1, "alpha" = .2)
  posterior <- matrix(samples, nrow = 1)

  expect_error(
    JAGS_bridgesampling_posterior(posterior, list(bias = bias)),
    "bias mixture priors"
  )
  expect_error(
    JAGS_marglik_priors(samples, list(bias = bias)),
    "bias mixture priors"
  )
  expect_error(
    JAGS_marglik_parameters(samples, list(bias = bias)),
    "bias mixture priors"
  )
})

skip_refit_if_cached("JAGS-marglik")

# This file tests the JAGS marginal likelihood computation functions
# It uses simple models where the log marginal likelihood is known to be 0
# (for prior samples, the marginal likelihood for any proper prior is 1, log(1) = 0)
# More complex consistency tests (e.g., including formulas etc part of `test-00-model-fits.R`)

test_that("JAGS model functions work (simple)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("normal", list(0, 1), list(1, Inf)),
    p3  = prior("lognormal", list(0, .5)),
    p4  = prior("t", list(0, .5, 5)),
    p5  = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    p6  = prior("gamma", list(2, 1)),
    p7  = prior("invgamma", list(3, 2), list(1, 3)),
    p8  = prior("exp", list(1.5)),
    p9  = prior("beta", list(3, 2)),
    p10 = prior("uniform", list(1, 5)),
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
    #p12 = prior("bernoulli", list(0.75)) discrete priors are not supported with bridgesampling
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

# skip the rest as it takes too long
skip_on_cran()

test_that("JAGS model functions work (vector)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("mnormal", list(mean = 0, sd = 1, K = 3),),
    p2  = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    p3  = prior("mt",      list(location = 2, scale = 0.5, df = 5, K = 2))
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 5*1e-2) # the mCauchy is a bit more variable
  }

})

test_that("JAGS model functions work (factor)", {

  skip_if_not_installed("rjags")
  all_priors   <- list(
    p1  = prior_factor("mnorm", list(mean = 0, sd = 1),    contrast = "orthonormal"),
    p2  = prior_factor("beta",  list(alpha = 1, beta = 1), contrast = "treatment"),
    p3  = prior_factor("beta",  list(alpha = 2, beta = 2), contrast = "treatment"),
    p4  = prior_factor("gamma",   list(shape = 2, rate = 3), contrast = "independent"),
    p5  = prior_factor("uniform", list(a = -0.5, b = 1.5),   contrast = "independent"),
    p6  = prior_factor("mnorm", list(mean = 0, sd = 1),     contrast = "meandif")
  )

  # add levels
  attr(all_priors[[1]], "levels") <- 3
  attr(all_priors[[2]], "levels") <- 2
  attr(all_priors[[3]], "levels") <- 3
  attr(all_priors[[4]], "levels") <- 1
  attr(all_priors[[5]], "levels") <- 3
  attr(all_priors[[6]], "levels") <- 3
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (spike and slab)", {
  skip("Marginal likelihood computation for spike and slab priors is not implemented.")
  skip_if_not_installed("rjags")
  all_priors   <- list(
    p1  = prior_spike_and_slab(prior("normal",   list(0, 1)), prior_inclusion = prior("beta", list(1, 1))),
    p2  = prior_spike_and_slab(prior("gamma",    list(3, 4)), prior_inclusion = prior("beta", list(5, 1))),
    p3  = prior_spike_and_slab(prior("invgamma", list(4, 5)), prior_inclusion = prior("point", list(.3)))
  )

  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-3)
  }

})

test_that("JAGS model functions work (weightfunctions)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1))),
    prior_weightfunction("one-sided", c(.05, 0.10), wf_cumulative(c(1, 2, 3))),
    prior_weightfunction("one-sided", c(.05, 0.60), wf_independent(prior("beta", list(1, 1)))),
    prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)))
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (spikes)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1    = prior("spike", list(1)),
    p2.2  = prior_factor("spike", list(location = 2), contrast = "treatment"),
    p3.2  = prior_factor("spike", list(location = 3), contrast = "independent"),
    p4.2  = prior_factor("spike", list(location = 0), contrast = "orthonormal"),
    p5.2  = prior_factor("spike", list(location = 0), contrast = "meandif"),
    p2.5  = prior_factor("spike", list(location = 2), contrast = "treatment"),
    p3.5  = prior_factor("spike", list(location = 3), contrast = "independent"),
    p4.5  = prior_factor("spike", list(location = 0), contrast = "orthonormal"),
    p5.5  = prior_factor("spike", list(location = 0), contrast = "meandif")
  )
  attr(all_priors$p2.2, "levels") <- 2
  attr(all_priors$p3.2, "levels") <- 2
  attr(all_priors$p4.2, "levels") <- 2
  attr(all_priors$p5.2, "levels") <- 2
  attr(all_priors$p2.5, "levels") <- 2
  attr(all_priors$p3.5, "levels") <- 2
  attr(all_priors$p4.5, "levels") <- 2
  attr(all_priors$p5.5, "levels") <- 2
  nuisance_prior <- list(sigma = prior("normal", list(0, 1)))
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- c(all_priors[i], nuisance_prior)
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("bridge sampling object function works",{

  marglik0 <- bridgesampling_object()
  marglik1 <- bridgesampling_object(1)

  expect_equal(marglik0$logml, -Inf)
  expect_equal(marglik1$logml, 1)
  expect_s3_class(marglik0, "bridge")

})

test_that("JAGS marglik with formula works", {

  # Test marginal likelihood computation with formula interface
  # Uses intercept-only formula with various priors
  # When sampling from prior and computing marglik, the result should be ~0 (log(1))

  skip_if_not_installed("rjags")

  # Simple data for the formula
  set.seed(1)
  df_test <- data.frame(x = rnorm(10))
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Create formula prior list with intercept only
  prior_list <- list(
    "intercept" = prior("gamma",  list(2, 2)),
    "x"         = prior("normal", list(0, 1))
  )

  # Process formula to get JAGS syntax
  formula_result <- JAGS_formula(~ 1 + x, parameter = "mu", data = df_test, prior_list = prior_list)

  # Build JAGS model with formula priors
  model_syntax <- JAGS_add_priors("model{}", formula_result$prior_list)
  monitor      <- JAGS_to_monitor(formula_result$prior_list)
  inits        <- JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1)

  # Sample from prior using JAGS
  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")

  # Compute marginal likelihood using formula interface
  marglik <- JAGS_bridgesampling(
    fit                = samples,
    log_posterior      = log_posterior,
    data               = list(),
    prior_list         = NULL,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = list(mu = df_test),
    formula_prior_list = list(mu = prior_list)
  )

  expect_equal(marglik$logml, 0, tolerance = 1e-3)
})

test_that("JAGS marglik with exp(intercept) formula works", {

  # Test marginal likelihood computation with formula interface
  # Uses intercept-only formula with various priors
  # When sampling from prior and computing marglik, the result should be ~0 (log(1))

  skip_if_not_installed("rjags")

  # Simple data for the formula
  set.seed(1)
  df_test <- data.frame(x = rnorm(10))
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Create formula prior list with intercept only
  prior_list <- list(
    "intercept" = prior("gamma",  list(2, 2)),
    "x"         = prior("normal", list(0, 1))
  )

  # Process formula to get JAGS syntax
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- TRUE
  formula_result <- JAGS_formula(formula, parameter = "mu", data = df_test, prior_list = prior_list)
  expect_equal(formula_result$formula_syntax, "for(i in 1:N_mu){\n  mu[i] = log(mu_intercept) + mu_x * mu_data_x[i]\n}\n")

  # Build JAGS model with formula priors
  model_syntax <- JAGS_add_priors("model{}", formula_result$prior_list)
  monitor      <- JAGS_to_monitor(formula_result$prior_list)
  inits        <- JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1)

  # Sample from prior using JAGS
  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")

  # Compute marginal likelihood using formula interface
  marglik <- JAGS_bridgesampling(
    fit                = samples,
    log_posterior      = log_posterior,
    data               = list(),
    prior_list         = NULL,
    formula_list       = list(mu = formula),
    formula_data_list  = list(mu = df_test),
    formula_prior_list = list(mu = prior_list)
  )

  expect_equal(marglik$logml, 0, tolerance = 1e-3)
})

test_that("JAGS bridgesampling infers formula scaling metadata from fits", {

  df_test <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x"         = prior("normal", list(0, 1))
  )

  scaled_formula <- JAGS_formula(
    ~ 1 + x,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  fit <- matrix(c(0, 1), ncol = 2)
  attr(fit, "formula_scale") <- list(mu = scaled_formula$formula_scale)

  inferred_scale <- BayesTools:::.JAGS_formula_scale_list_from_fit(fit, "mu")
  expect_equal(inferred_scale, list(mu = list(x = TRUE)))

  bridge_formula <- JAGS_formula(
    ~ 1 + x,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    formula_scale = inferred_scale$mu
  )

  expect_equal(bridge_formula$data$mu_data_x, scaled_formula$data$mu_data_x)
})

test_that("JAGS formula marglik reconstructs inverse-gamma terms on original scale", {

  samples <- c(
    "inv_mu_intercept" = 2,
    "inv_mu_x"         = 4
  )
  formula_data_list <- list(
    mu = list(
      N_mu      = 2,
      mu_data_x = c(10, 20)
    )
  )
  formula_prior_list <- list(
    mu = list(
      mu_intercept = prior("invgamma", list(2, 1)),
      mu_x         = prior("invgamma", list(2, 1))
    )
  )

  parameters <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(parameters$mu, c(0.5 + 0.25 * 10, 0.5 + 0.25 * 20))

  formula_log_intercept <- ~ 1 + x
  attr(formula_log_intercept, "log(intercept)") <- TRUE
  parameters_log <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = formula_log_intercept),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(parameters_log$mu, c(log(0.5) + 0.25 * 10, log(0.5) + 0.25 * 20))
})


test_that("JAGS marglik reconstructs indexed factor inverse-gamma auxiliaries", {
  theta_prior <- prior_factor("invgamma", list(2, 1), contrast = "independent")
  theta_prior$parameters$K <- 2

  parameters <- BayesTools:::.JAGS_marglik_parameters.factor(
    samples = c("inv_theta[1]" = 2, "inv_theta[2]" = 4),
    prior = theta_prior,
    parameter_name = "theta"
  )

  expect_equal(parameters$theta, c(0.5, 0.25))
})


test_that("JAGS formula marglik preserves predictor names containing _data", {
  samples <- c(
    "mu_intercept" = 1,
    "mu_x_data"   = 2
  )
  formula_data_list <- list(
    mu = list(
      N_mu           = 2,
      mu_data_x_data = c(10, 20)
    )
  )
  formula_prior_list <- list(
    mu = list(
      mu_intercept = prior("normal", list(0, 1)),
      mu_x_data    = prior("normal", list(0, 1))
    )
  )

  parameters <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = ~ 1 + x_data),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(parameters$mu, c(1 + 2 * 10, 1 + 2 * 20))
})


# Targeted tests for uncovered code paths in JAGS-marglik.R

test_that("JAGS_bridgesampling_posterior input validation works", {

  posterior <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(posterior) <- c("mu", "sigma", "x")

  # Input validation errors

  expect_error(JAGS_bridgesampling_posterior(data.frame(x = 1), prior_list = NULL), "'posterior' must be a matrix")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = "x"), "'prior_list' must be a list.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = 1), "'add_parameters' must be a character")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = "x"), "'add_bounds' must be a list")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(a = 1)), "'add_bounds' must contain lower and upper bounds")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = c("x", "y"), add_bounds = list(lb = 0, ub = 1)), "lb' and 'ub' must have the same lenght")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(lb = "a", ub = "b")), "lb' and 'ub' must be numeric")

  # Unsupported prior types
  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list = list(p1 = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1))))),
    "spike and slab"
  )
  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list = list(p1 = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1))), is_null = c(TRUE, FALSE)))),
    "prior mixture"
  )

  # Missing parameters
  posterior_small <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(posterior_small) <- c("a", "b")
  expect_error(JAGS_bridgesampling_posterior(posterior_small, prior_list = list(x = prior("normal", list(0, 1)))), "'posterior' does not contain all")

  # Successful case with add_parameters
  result <- JAGS_bridgesampling_posterior(posterior, prior_list = list(mu = prior("normal", list(0, 1))), add_parameters = "x", add_bounds = list(lb = -Inf, ub = Inf))
  expect_true(is.matrix(result))
  expect_true("x" %in% colnames(result))

})

test_that("JAGS_marglik_priors input validation and edge cases work", {

  # Empty prior_list returns empty list

  expect_equal(JAGS_marglik_priors(list(), prior_list = list()), list())

  # Input validation
  expect_error(JAGS_marglik_priors(list(), prior_list = "x"), "'prior_list' must be a list.")
  expect_error(JAGS_marglik_priors(list(), prior_list = prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")
  expect_error(JAGS_marglik_priors(list(), prior_list = list(x = 1)), "'prior_list' must be a list of priors.")

})

test_that("JAGS_marglik_parameters input validation and edge cases work", {

  # Test: empty prior_list returns empty list
  result <- JAGS_marglik_parameters(list(), prior_list = list())
  expect_equal(result, list())

  # Test: prior_list must be a list
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = "not_a_list"),
    "'prior_list' must be a list."
  )

  # Test: prior_list must be a list of priors (single prior passed)
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = prior("normal", list(0, 1))),
    "'prior_list' must be a list of priors."
  )

  # Test: prior_list must be a list of priors (non-prior elements)
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = list(x = 1)),
    "'prior_list' must be a list of priors."
  )

})

test_that(".fit_to_posterior handles different input types", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)

  # mcmc.list (rjags::coda.samples)
  samples_mcmc_list <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 100, quiet = TRUE, progress.bar = "none")
  marglik <- JAGS_bridgesampling(samples_mcmc_list, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik, "bridge")

  # mcmc (coda::as.mcmc)
  samples_mcmc <- coda::as.mcmc(samples_mcmc_list[[1]])
  marglik_mcmc <- JAGS_bridgesampling(samples_mcmc, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_mcmc, "bridge")

  # Error for unsupported input
  expect_error(JAGS_bridgesampling("bad_input", prior_list = prior_list, data = list(), log_posterior = log_posterior), "not implemented")

})

test_that(".fit_to_posterior handles jags.samples output", {

  skip_if_not_installed("rjags")

  # Scalar parameter
  prior_list <- list(mu = prior("normal", list(0, 1)), sigma = prior("gamma", list(1, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples_jags <- rjags::jags.samples(model = model, variable.names = monitor, n.iter = 100, progress.bar = "none")
  marglik_jags <- JAGS_bridgesampling(samples_jags, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_jags, "bridge")

})

test_that(".fit_to_posterior handles vector parameters in jags.samples", {

  skip_if_not_installed("rjags")

  # Vector parameter (K > 1)
  prior_list <- list(p = prior("mnormal", list(mean = 0, sd = 1, K = 3)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples_jags <- rjags::jags.samples(model = model, variable.names = monitor, n.iter = 100, progress.bar = "none")
  marglik_jags <- JAGS_bridgesampling(samples_jags, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_jags, "bridge")

})

test_that("JAGS_bridgesampling handles runjags output", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  log_posterior <- STANDARD_LOG_POSTERIOR
  old_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.runjags = old_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.runjags = TRUE)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 2,
    adapt = 100,
    burnin = 100,
    sample = 500,
    silent.jags = TRUE,
    modules = "glm"
  ))

  marglik <- JAGS_bridgesampling(fit, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik, "bridge")
  expect_equal(marglik$logml, 0, tolerance = 0.1)

})
