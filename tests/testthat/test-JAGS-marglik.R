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

make_bridge_random_fixture <- function(formula = ~ 1 + x + us(1 + x | id),
                                       data = NULL,
                                       prior_list = NULL,
                                       prior_random_list = NULL){

  if(is.null(data)){
    data <- data.frame(
      x = c(-1, 0, 1, 2),
      id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
    )
  }
  if(is.null(prior_list)){
    prior_list <- list(intercept = prior("normal", list(0, 1)))
    if("x" %in% all.vars(formula)){
      prior_list$x <- prior("normal", list(0, 1))
    }
  }
  if(is.null(prior_random_list)){
    prior_random_list <- prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1, include_correlation = FALSE)
      )
    )
  }

  list(
    data = data,
    formula = formula,
    prior_list = prior_list,
    prior_random_list = prior_random_list,
    result = JAGS_formula(
      formula = formula,
      parameter = "mu",
      data = data,
      prior_list = prior_list,
      prior_random = prior_random_list
    )
  )
}

expect_formula_random_prior_only_bridge <- function(formula, data, prior_list,
                                                    prior_random_list,
                                                    n_iter = 8000,
                                                    tolerance = 0.08,
                                                    seed = 1,
                                                    maxiter = 2000){

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(seed)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = seed),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = n_iter,
    quiet = TRUE,
    progress.bar = "none"
  )
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = maxiter
  )

  expect_s3_class(marglik, "bridge")
  expect_equal(marglik$logml, 0, tolerance = tolerance)

  invisible(list(
    marglik = marglik,
    formula_result = formula_result,
    samples = samples
  ))
}

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

test_that("JAGS marginal-likelihood helpers reject spike-and-slab priors explicitly", {

  prior_list <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    )
  )
  posterior <- matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "theta"))
  samples <- c(theta = 0, theta_inclusion = .5)

  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list),
    "spike and slab priors is not implemented"
  )
  expect_error(
    JAGS_marglik_priors(samples, prior_list),
    "prior mixture priors is not implemented"
  )
  expect_error(
    JAGS_marglik_parameters(samples, prior_list),
    "prior mixture priors is not implemented"
  )
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

test_that("JAGS bridgesampling posterior supports add-only parameters", {

  posterior <- matrix(
    c(
      -1, 0.2,
       0, 0.5,
       1, 0.8
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(posterior) <- c("x", "prob")

  info <- BayesTools:::.JAGS_bridgesampling_posterior_info(NULL)
  expect_type(info, "character")
  expect_equal(length(info), 0L)
  expect_equal(attr(info, "lb"), numeric())
  expect_equal(attr(info, "ub"), numeric())

  result_null <- JAGS_bridgesampling_posterior(
    posterior = posterior,
    prior_list = NULL,
    add_parameters = c("x", "prob"),
    add_bounds = list(
      lb = c(x = -Inf, prob = 0),
      ub = c(x = Inf, prob = 1)
    )
  )
  expect_equal(colnames(result_null), c("x", "prob"))
  expect_equal(attr(result_null, "lb"), c(x = -Inf, prob = 0))
  expect_equal(attr(result_null, "ub"), c(x = Inf, prob = 1))

  result_empty <- JAGS_bridgesampling_posterior(
    posterior = posterior,
    prior_list = list(),
    add_parameters = "x",
    add_bounds = list(lb = -Inf, ub = Inf)
  )
  expect_equal(colnames(result_empty), "x")
  expect_equal(attr(result_empty, "lb"), c(x = -Inf))
  expect_equal(attr(result_empty, "ub"), c(x = Inf))
})

test_that("JAGS bridgesampling validates rebuilt formula random design metadata", {

  fixture <- make_bridge_random_fixture()
  fitted <- list(mu = fixture$result$formula_design)
  rebuilt <- list(mu = fixture$result$formula_design)

  expect_silent(BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, rebuilt))

  changed <- rebuilt
  changed$mu$random_effects[[1]]$structure <- "diag"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "covariance structure",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$structure <- NULL
  changed$mu$random_effects[[1]]$covariance <- "us"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "missing canonical 'random_term\\$structure'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$homogeneous_sd <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$correlation <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "missing canonical 'random_term\\$correlation'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$block_name <- "other"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "block names",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$group_levels <- rev(changed$mu$random_effects[[1]]$group_levels)
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "group levels",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$column_names[1] <- "changed"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "column names",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$n_columns <- changed$mu$random_effects[[1]]$n_columns + 1L
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "group or column counts",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$model_matrix <- changed$mu$random_effects[[1]]$model_matrix[-1, , drop = FALSE]
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "model matrix shape or columns",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$group_map <- rev(changed$mu$random_effects[[1]]$group_map)
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "group map",
    fixed = TRUE
  )

  structured_fixture <- make_bridge_random_fixture(
    formula = ~ 1 + ar1(idx | id),
    data = data.frame(
      idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
      id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  structured_fitted <- list(mu = structured_fixture$result$formula_design)
  structured_changed <- structured_fitted
  structured_changed$mu$random_effects[[1]]$structured_index$label <- "changed"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      structured_fitted,
      structured_changed
    ),
    "structured random-effect index metadata",
    fixed = TRUE
  )

  allocation_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )
  allocation_formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  allocation_prior_list <- list(intercept = prior("normal", list(0, 1)))
  allocation_fixture <- make_bridge_random_fixture(
    formula = allocation_formula,
    data = allocation_data,
    prior_list = allocation_prior_list,
    prior_random_list = prior_random(
      allocation = random_variance_allocation(
        sd = prior("gamma", list(2, 2)),
        allocation = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  independent_fixture <- make_bridge_random_fixture(
    formula = allocation_formula,
    data = allocation_data,
    prior_list = allocation_prior_list,
    prior_random_list = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2))),
      drug = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      list(mu = allocation_fixture$result$formula_design),
      list(mu = independent_fixture$result$formula_design)
    ),
    "scale/allocation",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling rejects fitted/rebuilt random design mismatches", {

  fixture <- make_bridge_random_fixture()
  fit <- coda::mcmc(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "dummy")))
  attr(fit, "formula_design") <- list(mu = fixture$result$formula_design)

  mismatch_data <- fixture$data
  mismatch_data$id <- factor(as.character(mismatch_data$id), levels = c("b", "a"))

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = NULL,
      formula_list = list(mu = fixture$formula),
      formula_data_list = list(mu = mismatch_data),
      formula_prior_list = list(mu = fixture$prior_list),
      formula_random_prior_list = list(mu = fixture$prior_random_list),
      maxiter = 10
    ),
    "group levels",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling supports formula random effects through prior_random", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  if(!isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE))){
    skip("BayesTools JAGS module is not available")
  }

  df_test <- data.frame(
    x = c(-0.5, 0.5),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 + x + us(1 + x | id)
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior_lkj(eta = 1, include_correlation = FALSE)
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_true("mu__xREx__id_xRE_CORx_lkj_u[1]" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 5000,
    quiet = TRUE,
    progress.bar = "none"
  )

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000
  )

  expect_s3_class(marglik, "bridge")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling supports continuous-time CAR formula random effects", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 + car(0 + time | id)
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      rho = prior("normal", list(0, 0.5))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_equal(formula_result$formula_design$random_effects[[1]]$structure, "car")
  expect_equal(formula_result$formula_design$random_effects[[1]]$correlation$bounds, c(lower = 0, upper = 1))
  expect_true("mu__xREx__id_rho" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000
  )

  expect_s3_class(marglik, "bridge")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling supports Dirichlet variance-allocation random effects", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    allocation = random_variance_allocation(
      sd = prior("gamma", list(2, 2)),
      allocation = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_true("mu__xRE_ALLOCx_allocation_weight" %in% names(formula_result$prior_list))
  expect_true("mu__xREx__study_xRE_Zx" %in% formula_result$add_parameters)
  expect_true("mu__xREx__drug_xRE_Zx" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000
  )

  expect_s3_class(marglik, "bridge")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling gives unit marglik for prior-only random-effect settings", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  continuous_data <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"))
  )
  factor_data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 3), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 3), levels = c("g1", "g2", "g3"))
  )

  diag_heterogeneous <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + diag(1 + x | id, hom = FALSE),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        terms = list(x = prior("gamma", list(3, 1)))
      )
    ),
    seed = 10
  )
  expect_false(diag_heterogeneous$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    diag_heterogeneous$formula_result$formula_design$random_effects[[1]]$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_x")
  )

  id_homogeneous <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + id(1 + x | id),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        monitor = random_monitor(
          latent = TRUE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    ),
    seed = 11
  )
  expect_true(id_homogeneous$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    unique(id_homogeneous$formula_result$formula_design$random_effects[[1]]$sd_parameter_names),
    "mu__xREx__id_sd"
  )

  lkj_syntax <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + us(1 + x | id),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(
          eta = 2,
          backend = "syntax",
          include_correlation = FALSE
        )
      )
    ),
    seed = 12
  )
  expect_equal(lkj_syntax$formula_result$jags_modules, character())
  expect_true(any(grepl("_xRE_CORx_lkj_u", lkj_syntax$formula_result$add_parameters, fixed = TRUE)))

  cs_fisher_z <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + cs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = prior("normal", list(0, 0.5))
      )
    ),
    seed = 13
  )
  expect_true(cs_fisher_z$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    cs_fisher_z$formula_result$formula_design$random_effects[[1]]$correlation$rho_scale,
    "fisher_z"
  )
  expect_true("mu__xREx__id_rho_z" %in% names(cs_fisher_z$formula_result$prior_list))

  hcs_logit <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + hcs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5)),
          rho_scale = "logit"
        )
      )
    ),
    seed = 14
  )
  expect_false(hcs_logit$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_true("mu__xREx__id_rho_logit" %in% names(hcs_logit$formula_result$prior_list))
  expect_equal(
    hcs_logit$formula_result$formula_design$random_effects[[1]]$correlation$rho_scale,
    "logit"
  )

  ar1_fisher_z <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + ar1(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = prior("normal", list(0, 0.5))
      )
    ),
    seed = 15
  )
  expect_equal(ar1_fisher_z$formula_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_true("mu__xREx__id_rho_z" %in% names(ar1_fisher_z$formula_result$prior_list))

  har_raw <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + har(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5), truncation = list(lower = -1, upper = 1)),
          rho_scale = "rho"
        )
      )
    ),
    seed = 16
  )
  expect_equal(har_raw$formula_result$formula_design$random_effects[[1]]$structure, "har")
  expect_false(har_raw$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_true("mu__xREx__id_rho" %in% names(har_raw$formula_result$prior_list))
  expect_false("mu__xREx__id_rho_z" %in% names(har_raw$formula_result$prior_list))

  nested_allocation <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    data = data.frame(
      study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
      paper = factor(c("p1", "p2", "p1", "p2"), levels = c("p1", "p2")),
      drug = factor(c("d1", "d1", "d2", "d2"), levels = c("d1", "d2"))
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        allocation = prior("dirichlet", list(alpha = c(3, 2)))
      )
    ),
    n_iter = 10000,
    seed = 17
  )
  expect_true("mu__xRE_ALLOCx_total_re_weight" %in% names(nested_allocation$formula_result$prior_list))
  expect_true("mu__xRE_ALLOCx_nested_split_weight" %in% names(nested_allocation$formula_result$prior_list))
  expect_equal(
    length(nested_allocation$formula_result$formula_design$random_effects[[1]]$allocation$factors),
    2L
  )

  sd_leaf_allocation <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + hcs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      allocation = random_variance_allocation(
        name = "leaf_alloc",
        terms = "id",
        components = "sd",
        scale = "mean_variance",
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = c(2, 3, 4)))
      ),
      id = random_block(rho = prior("normal", list(0, 0.5)))
    ),
    n_iter = 10000,
    seed = 18
  )
  expect_true("mu__xRE_ALLOCx_leaf_alloc_weight" %in% names(sd_leaf_allocation$formula_result$prior_list))
  expect_equal(
    sd_leaf_allocation$formula_result$formula_design$random_effects[[1]]$allocation$components,
    "sd"
  )
  expect_equal(
    sd_leaf_allocation$formula_result$formula_design$random_effects[[1]]$allocation$scale,
    "mean_variance"
  )
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
  attr(theta_prior, "levels") <- 2

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
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = c("x", "y"), add_bounds = list(lb = 0, ub = 1)), "'lb' and 'ub' must have the same length")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(lb = "a", ub = "b")), "'lb' and 'ub' must be numeric")
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = stats::setNames(-Inf, "wrong"), ub = stats::setNames(Inf, "x"))
    ),
    "names must match",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = 1, ub = 0)
    ),
    "smaller than upper",
    fixed = TRUE
  )

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

  # Empty prior_list contributes zero log prior density

  expect_equal(JAGS_marglik_priors(list(), prior_list = list()), 0)
  expect_equal(JAGS_marglik_priors_formula(list(), formula_prior_list = list(mu = list())), 0)

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
