skip_if_not_test_profile("unit")

test_that("compiled simple-prior densities preserve the canonical calculation", {

  priors <- list(
    prior("normal", list(0, 1)),
    prior("normal", list(0, 1), list(0, 2)),
    prior("gamma", list(2, 1), list(.2, 3)),
    prior("point", list(.25))
  )
  values <- list(
    c(-2, 0, 2, NA_real_),
    c(-1, 0, .5, 2, 3, NA_real_),
    c(0, .2, 1, 3, 4, NA_real_),
    c(0, .25, 1, NA_real_)
  )

  for(i in seq_along(priors)){
    canonical <- BayesTools:::.prior_simple_lpdf(priors[[i]], values[[i]])
    compiled  <- BayesTools:::.prior_simple_lpdf_evaluator(priors[[i]])
    expect_identical(compiled(values[[i]]), canonical)
  }
})

test_that("compiled factor densities preserve product priors and boundaries", {

  values <- c(.25, .5, 1, 1.5)
  samples <- stats::setNames(values, paste0("theta[", seq_along(values), "]"))
  for(contrast in c("independent", "treatment")){
    specifications <- list(
      prior_factor("normal", list(0, 1), list(0, 2), contrast = contrast),
      prior_factor("invgamma", list(2, 1), contrast = contrast),
      prior_factor("point", list(.25), contrast = contrast)
    )
    expected <- c(
      sum(stats::dnorm(values, log = TRUE)) -
        length(values) * log(stats::pnorm(2) - stats::pnorm(0)),
      sum(stats::dgamma(1 / values, shape = 2, rate = 1, log = TRUE) -
            2 * log(values)),
      0
    )
    for(i in seq_along(specifications)){
      prior <- specifications[[i]]
      attr(prior, "levels") <- length(values) + as.integer(contrast == "treatment")
      evaluator <- BayesTools:::.bt_JAGS_bridge_compile_factor_evaluator(prior, "theta")
      expect_equal(evaluator$log_prior(samples), expected[[i]], tolerance = 1e-12)
      expect_equal(evaluator$log_prior(as.list(samples)), expected[[i]], tolerance = 1e-12)
      expect_equal(evaluator$parameters(samples)$theta,
                   if(i == 3L) rep(.25, length(values)) else values)
      if(i < 3L){
        invalid <- samples
        invalid[[2L]] <- -1
        expect_equal(evaluator$log_prior(invalid), -Inf)
      }
      if(i == 1L){
        expect_error(
          evaluator$log_prior(samples[-1L]),
          "'samples' does not contain all monitored factor prior parameters.",
          fixed = TRUE
        )
      }
    }
  }
})

test_that("compiled bridge prior evaluators match public marglik helpers", {

  theta_prior <- prior_factor("invgamma", list(2, 1), list(.1, 2), contrast = "independent")
  attr(theta_prior, "levels") <- 2

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))
  phacking <- prior_phacking(form = "linear", alpha = prior("beta", list(2, 3)))
  bias <- prior_bias(selection, phacking)

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    none  = prior_none(),
    sigma = prior("invgamma", list(3, 2), list(1, 3)),
    point = prior("point", list(.25)),
    v     = prior("mnormal", list(0, 1, 2)),
    z     = prior("mpoint", list(0, 2)),
    w     = prior("dirichlet", list(alpha = c(2, 3))),
    theta = theta_prior,
    PET   = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(2, 1)),
    bias  = bias
  )

  samples <- c(
    mu = .2,
    sigma = 2,
    "v[1]" = -.1,
    "v[2]" = .3,
    "prior_par_eta_w[1]" = 1.2,
    "prior_par_eta_w[2]" = 2.4,
    "theta[1]" = .5,
    "theta[2]" = 1 / 3,
    PET = .4,
    PEESE = .6,
    "omega[2]" = .625,
    alpha = .4
  )

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)

  expect_equal(
    compiled$log_prior(samples),
    JAGS_marglik_priors(samples, prior_list),
    tolerance = 1e-12
  )
  expect_equal(
    compiled$parameters(samples),
    JAGS_marglik_parameters(samples, prior_list),
    tolerance = 1e-12
  )

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma test.
  legacy_samples <- samples[!names(samples) %in% c("sigma", "theta[1]", "theta[2]")]
  legacy_samples <- c(
    legacy_samples,
    inv_sigma = .5,
    "inv_theta[1]" = 2,
    "inv_theta[2]" = 3
  )
  expect_equal(
    compiled$parameters(legacy_samples),
    compiled$parameters(samples),
    tolerance = 1e-12
  )
})

test_that("compiled bridge prior evaluators preserve positive support behavior", {

  prior_list <- list(
    sigma = prior("invgamma", list(3, 2), list(1, 3)),
    w     = prior("dirichlet", list(alpha = c(2, 3)))
  )
  samples <- c(
    sigma = 2,
    "prior_par_eta_w[1]" = 1.2,
    "prior_par_eta_w[2]" = 2.4
  )
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)

  bad_invgamma <- samples
  bad_invgamma[["sigma"]] <- 0
  expect_equal(compiled$log_prior(bad_invgamma), -Inf)
  invgamma_error <- tryCatch(
    compiled$parameters(bad_invgamma),
    error = function(e) e
  )
  expect_s3_class(invgamma_error, "BayesTools_marglik_out_of_support")

  bad_dirichlet <- samples
  bad_dirichlet[["prior_par_eta_w[1]"]] <- 0
  expect_equal(compiled$log_prior(bad_dirichlet), -Inf)
  dirichlet_error <- tryCatch(
    compiled$parameters(bad_dirichlet),
    error = function(e) e
  )
  expect_s3_class(dirichlet_error, "BayesTools_marglik_out_of_support")

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma test.
  legacy_invgamma <- samples[!names(samples) %in% "sigma"]
  legacy_invgamma[["inv_sigma"]] <- .5
  expect_equal(
    compiled$parameters(legacy_invgamma),
    compiled$parameters(samples),
    tolerance = 1e-12
  )
})

test_that("row prior evaluation preserves joint prior boundaries", {

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    sigma = prior("invgamma", list(3, 2), list(1, 3)),
    w     = prior("dirichlet", list(alpha = c(2, 3)))
  )
  samples <- cbind(
    mu = c(-.2, .4, .8),
    sigma = c(1.2, 2.1, 1.6),
    "prior_par_eta_w[1]" = c(1.2, 1.8, 2.1),
    "prior_par_eta_w[2]" = c(2.4, 1.1, 3.2)
  )
  expected <- vapply(seq_len(nrow(samples)), function(i){
    JAGS_marglik_priors(samples[i, ], prior_list)
  }, numeric(1))

  expect_equal(JAGS_marglik_priors_rows(samples, prior_list), expected)
  evaluator <- JAGS_marglik_priors_rows_evaluator(prior_list)
  expect_equal(evaluator(samples), expected)
  expect_equal(evaluator(samples[c(3L, 1L), , drop = FALSE]), expected[c(3L, 1L)])
  expect_equal(
    JAGS_marglik_priors_rows(samples[1L, ], prior_list),
    expected[1L]
  )

  invalid <- samples
  invalid[2L, "prior_par_eta_w[1]"] <- 0
  invalid[3L, "sigma"] <- 0
  invalid_density <- JAGS_marglik_priors_rows(invalid, prior_list)
  expect_true(is.finite(invalid_density[1L]))
  expect_equal(invalid_density[2L], -Inf)
  expect_equal(invalid_density[3L], -Inf)
})

test_that("row prior evaluation uses PET and PEESE backend coordinates", {

  prior_list <- list(
    location   = prior("normal", list(0, 1)),
    bias_pet   = prior_PET("normal", list(0, 2), list(0, 3)),
    bias_peese = prior_PEESE("gamma", list(2, 1)),
    allocation = prior("dirichlet", list(alpha = c(2, 3)))
  )
  samples <- cbind(
    location = c(-.2, .4, .8, 0, -.1),
    PET = c(0, .3, 3, 3.1, 1),
    PEESE = c(.2, 1, 2, 1, 0),
    "prior_par_eta_allocation[1]" = c(1.2, 1.8, 2.1, 1, 1),
    "prior_par_eta_allocation[2]" = c(2.4, 1.1, 3.2, 1, 1)
  )
  expected <- stats::dnorm(samples[, "location"], log = TRUE) +
    stats::dnorm(samples[, "PET"], sd = 2, log = TRUE) -
    log(stats::pnorm(3 / 2) - .5) +
    stats::dgamma(samples[, "PEESE"], shape = 2, rate = 1, log = TRUE) +
    stats::dgamma(samples[, "prior_par_eta_allocation[1]"],
                  shape = 2, rate = 1, log = TRUE) +
    stats::dgamma(samples[, "prior_par_eta_allocation[2]"],
                  shape = 3, rate = 1, log = TRUE)
  expected[samples[, "PET"] > 3] <- -Inf

  evaluator <- JAGS_marglik_priors_rows_evaluator(prior_list)
  expect_equal(evaluator(samples), expected, tolerance = 1e-12)
  expect_equal(evaluator(as.data.frame(samples)), expected, tolerance = 1e-12)
  expect_equal(evaluator(samples[3L, ]), expected[3L], tolerance = 1e-12)
})

test_that("row prior evaluation retains an exact scalar fallback", {

  prior_list <- list(
    sigma = prior("invgamma", list(3, 2), list(1, 3)),
    v     = prior("mnormal", list(0, 1, 2))
  )
  samples <- cbind(
    sigma = c(1.2, 2.1),
    "v[1]" = c(-.1, .3),
    "v[2]" = c(.4, -.2)
  )
  expected <- vapply(seq_len(nrow(samples)), function(i){
    JAGS_marglik_priors(samples[i, ], prior_list)
  }, numeric(1))

  expect_equal(JAGS_marglik_priors_rows(samples, prior_list), expected)
  expect_equal(
    JAGS_marglik_priors_rows(samples[, FALSE, drop = FALSE], list()),
    numeric(nrow(samples))
  )

  theta_prior <- prior_factor(
    "normal",
    parameters = list(0, 1),
    contrast   = "independent"
  )
  attr(theta_prior, "levels") <- 2L
  factor_samples <- cbind(
    "theta[1]" = c(-.2, .4),
    "theta[2]" = c(.3, -.1)
  )
  factor_list <- list(theta = theta_prior)
  factor_expected <- vapply(seq_len(nrow(factor_samples)), function(i){
    JAGS_marglik_priors(factor_samples[i, ], factor_list)
  }, numeric(1))
  expect_equal(
    JAGS_marglik_priors_rows(factor_samples, factor_list),
    factor_expected,
    tolerance = 0
  )
})

test_that("compiled formula prior evaluator matches public formula density helper", {

  formula_prior_list <- list(
    mu = list(
      mu_intercept = prior("normal", list(0, 1)),
      mu_x         = prior("invgamma", list(3, 2), list(1, 3))
    ),
    sigma = list(
      sigma_intercept = prior("gamma", list(2, 1))
    )
  )
  samples <- c(
    mu_intercept = .2,
    mu_x = 2,
    sigma_intercept = .8
  )

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(formula_prior_list)

  expect_equal(
    compiled$log_prior(samples),
    JAGS_marglik_priors_formula(samples, formula_prior_list),
    tolerance = 1e-12
  )
  expect_equal(
    compiled$parameters(samples),
    c(
      JAGS_marglik_parameters(samples, formula_prior_list$mu),
      JAGS_marglik_parameters(samples, formula_prior_list$sigma)
    ),
    tolerance = 1e-12
  )
})

test_that("shared ordered allocations contribute one prior density across scopes", {

  formula_data <- data.frame(
    f = ordered(c("low", "mid", "high"), levels = c("low", "mid", "high"))
  )
  make_prior_list <- function(parameter, id){
    JAGS_formula(
      formula = ~ f,
      parameter = parameter,
      data = formula_data,
      prior_list = list(
        intercept = prior("point", list(0)),
        f = prior_ordered(
          total = prior("normal", list(0, 1)),
          allocation = prior("dirichlet", list(alpha = c(2, 3))),
          id = id
        )
      )
    )$prior_list
  }

  formula_prior_list <- list(
    mu = make_prior_list("mu", "shared"),
    tau = make_prior_list("tau", "shared")
  )
  samples <- c(
    mu_f_ordered_total = 1,
    tau_f_ordered_total = 2,
    "prior_par_eta_ordered_alloc_shared_f[1]" = 1,
    "prior_par_eta_ordered_alloc_shared_f[2]" = 3
  )
  total_log_prior <- sum(stats::dnorm(c(1, 2), log = TRUE))
  allocation_log_prior <- sum(stats::dgamma(
    c(1, 3),
    shape = c(2, 3),
    rate = 1,
    log = TRUE
  ))
  expected <- total_log_prior + allocation_log_prior

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list
  )
  expect_equal(
    JAGS_marglik_priors_formula(samples, formula_prior_list),
    expected,
    tolerance = 1e-12
  )
  expect_equal(compiled$log_prior(samples), expected, tolerance = 1e-12)
  expect_length(compiled$allocation_keys, 1L)

  model_evaluators <-
    BayesTools:::.bt_JAGS_bridge_compile_model_prior_evaluators(
      prior_list = formula_prior_list$mu,
      formula_prior_list = list(tau = formula_prior_list$tau)
    )
  expect_equal(
    model_evaluators$prior$log_prior(samples),
    stats::dnorm(1, log = TRUE) + allocation_log_prior,
    tolerance = 1e-12
  )
  expect_equal(
    model_evaluators$formula$log_prior(samples),
    stats::dnorm(2, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    model_evaluators$prior$log_prior(samples) +
      model_evaluators$formula$log_prior(samples),
    expected,
    tolerance = 1e-12
  )
})

test_that("distinct ordered allocations retain independent prior densities", {

  formula_data <- data.frame(
    f = ordered(c("low", "mid", "high"), levels = c("low", "mid", "high"))
  )
  make_prior_list <- function(parameter, id){
    JAGS_formula(
      formula = ~ f,
      parameter = parameter,
      data = formula_data,
      prior_list = list(
        intercept = prior("point", list(0)),
        f = prior_ordered(
          total = prior("normal", list(0, 1)),
          allocation = prior("dirichlet", list(alpha = c(2, 3))),
          id = id
        )
      )
    )$prior_list
  }

  formula_prior_list <- list(
    mu = make_prior_list("mu", "mu_shape"),
    tau = make_prior_list("tau", "tau_shape")
  )
  samples <- c(
    mu_f_ordered_total = 1,
    tau_f_ordered_total = 2,
    "prior_par_eta_ordered_alloc_mu_shape_f[1]" = 1,
    "prior_par_eta_ordered_alloc_mu_shape_f[2]" = 3,
    "prior_par_eta_ordered_alloc_tau_shape_f[1]" = 2,
    "prior_par_eta_ordered_alloc_tau_shape_f[2]" = 4
  )
  mu_allocation_log_prior <- sum(stats::dgamma(
    c(1, 3),
    shape = c(2, 3),
    rate = 1,
    log = TRUE
  ))
  tau_allocation_log_prior <- sum(stats::dgamma(
    c(2, 4),
    shape = c(2, 3),
    rate = 1,
    log = TRUE
  ))
  expected <- sum(stats::dnorm(c(1, 2), log = TRUE)) +
    mu_allocation_log_prior +
    tau_allocation_log_prior

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list
  )
  expect_equal(
    JAGS_marglik_priors_formula(samples, formula_prior_list),
    expected,
    tolerance = 1e-12
  )
  expect_equal(compiled$log_prior(samples), expected, tolerance = 1e-12)
  expect_length(compiled$allocation_keys, 2L)

  model_evaluators <-
    BayesTools:::.bt_JAGS_bridge_compile_model_prior_evaluators(
      prior_list = formula_prior_list$mu,
      formula_prior_list = list(tau = formula_prior_list$tau)
    )
  expect_equal(
    model_evaluators$prior$log_prior(samples) +
      model_evaluators$formula$log_prior(samples),
    expected,
    tolerance = 1e-12
  )
  expect_length(model_evaluators$formula$allocation_keys, 2L)
})

test_that("bridge callback dispatcher exposes context only when requested", {

  context <- structure(
    list(state = c(theta = .1)),
    class = c("BayesTools_bridge_context", "list")
  )
  parameters <- list(theta = .1)

  expect_equal(
    BayesTools:::.bt_JAGS_bridge_call_log_posterior(
      log_posterior = function(parameters, data){
        expect_null(attr(parameters, "bridge_context", exact = TRUE))
        data$value
      },
      parameters = parameters,
      data = list(value = 1),
      context = context,
      bridge_context = FALSE
    ),
    1
  )

  expect_equal(
    BayesTools:::.bt_JAGS_bridge_call_log_posterior(
      log_posterior = function(parameters, data, bridge_context){
        expect_s3_class(bridge_context, "BayesTools_bridge_context")
        data$value
      },
      parameters = parameters,
      data = list(value = 2),
      context = context,
      bridge_context = TRUE
    ),
    2
  )

  nodes_context <- structure(
    list(nodes = c(theta = .1)),
    class = c(
      "BayesTools_bridge_nodes_context",
      "BayesTools_bridge_context",
      "list"
    )
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_call_log_posterior(
      log_posterior = function(parameters, data, bridge_context){
        expect_s3_class(
          bridge_context,
          "BayesTools_bridge_nodes_context"
        )
        data$value
      },
      parameters = parameters,
      data = list(value = 3),
      context = nodes_context,
      bridge_context = "nodes"
    ),
    3
  )

  expect_identical(
    BayesTools:::.bt_JAGS_bridge_context_mode(FALSE),
    "none"
  )
  expect_identical(
    BayesTools:::.bt_JAGS_bridge_context_mode(TRUE),
    "full"
  )
  expect_identical(
    BayesTools:::.bt_JAGS_bridge_context_mode("nodes"),
    "nodes"
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_context_mode("metadata"),
    "must be FALSE, TRUE"
  )
})

test_that("compiled random-effect prior evaluator matches public helper", {

  formula_data <- data.frame(
    id = factor(c("a", "b", "a", "c"))
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 1)),
        monitor = random_monitor(latent = TRUE)
      )
    )
  )
  random_term <- formula_output$formula_design$random_effects[[1]]
  z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  samples <- stats::setNames(c(.1, -.2, .3), z_names)
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_random_prior_evaluator(
    formula_design_list = list(mu = formula_output$formula_design)
  )

  expect_true(compiled$uses_posterior_row)
  expect_equal(
    compiled$log_prior(samples),
    bayestools_reference_random_effect_log_prior(
      samples = samples,
      random_term = random_term
    ),
    tolerance = 1e-12
  )
  expect_error(
    compiled$log_prior(samples[-1]),
    "standardized latent random effects"
  )
})

test_that("known group covariance random-effect prior uses MVN kernel density", {

  formula_data <- data.frame(
    id = factor(c("a", "b", "a", "c"), levels = c("a", "b", "c"))
  )
  K <- matrix(
    c(2, .4, .2,
      .4, 3, .5,
      .2, .5, 4),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  formula_output <- JAGS_formula(
    formula = random_effects,
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 1)),
        monitor = random_monitor(latent = TRUE)
      )
    )
  )
  random_term <- formula_output$formula_design$random_effects[[1]]
  z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  z <- c(.1, -.2, .3)
  samples <- stats::setNames(z, z_names)
  chol_K <- chol(K[random_term$group_levels, random_term$group_levels])
  expected <- -0.5 * (
    length(z) * log(2 * pi) +
      2 * sum(log(diag(chol_K))) +
      as.numeric(crossprod(z, chol2inv(chol_K) %*% z))
  )
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_random_prior_evaluator(
    formula_design_list = list(mu = formula_output$formula_design)
  )

  expect_equal(compiled$log_prior(samples), expected, tolerance = 1e-12)
  expect_equal(
    bayestools_reference_random_effect_log_prior(
      samples = samples,
      random_term = random_term
    ),
    expected,
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    expected,
    sum(stats::dnorm(z, mean = 0, sd = 1, log = TRUE)),
    tolerance = 1e-8
  )))

  changed_K <- K
  changed_K["a", "b"] <- changed_K["b", "a"] <- .1
  changed_random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(changed_K, scale = "none")
  )
  changed_output <- JAGS_formula(
    formula = changed_random_effects,
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 1)))
    )
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = formula_output$formula_design,
      rebuilt = changed_output$formula_design
    ),
    "group covariance metadata differ",
    fixed = TRUE
  )

  marginalized_output <- JAGS_formula(
    formula = random_effects,
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 1)))
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  changed_marginalized_output <- JAGS_formula(
    formula = changed_random_effects,
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 1)))
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = marginalized_output$formula_design,
      rebuilt = changed_marginalized_output$formula_design
    ),
    "group covariance metadata differ",
    fixed = TRUE
  )
})

test_that("compiled formula parameter evaluator matches design reconstruction", {

  formula_data <- data.frame(
    x = c(-1, 0, 2, 3),
    g = factor(c("a", "b", "a", "c"))
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 + x + g,
    parameter = "mu",
    data = formula_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("invgamma", list(3, 2), list(1, 3)),
      g         = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  formula_prior <- formula_output$prior_list
  attr(formula_prior$mu_x, "multiply_by") <- "x_scale"
  attr(formula_prior$mu_g, "multiply_by") <- 2
  formula_prior_list <- list(mu = formula_prior)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(
    mu_intercept = .2,
    mu_x = 2,
    "mu_g[1]" = .3,
    "mu_g[2]" = -.4
  )
  prior_list_parameters <- list(x_scale = 1.5)

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )

  expect_equal(
    compiled$parameters(samples, prior_list_parameters),
    JAGS_marglik_parameters_formula(
      samples = samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = prior_list_parameters,
      formula_design_list = formula_design_list,
      model_data = list()
    ),
    tolerance = 1e-12
  )
})

test_that("compiled random formula plans retain analytic contributions", {

  data <- data.frame(x = c(-1, 0, 2, 3), group = factor(c("a", "b", "a", "b")))
  for(structure in c("id", "us")){
    output <- JAGS_formula(
      formula = if(structure == "id") ~ id(1 + x | group) else ~ 1 + (1 | group),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("point", list(0))),
      prior_random = prior_random(group = random_block(sd = prior("point", list(2))))
    )
    term <- output$formula_design$random_effects[[1L]]
    latent <- if(structure == "id") c(.5, 2, -1, .25) else c(.5, 2)
    samples <- stats::setNames(latent, as.vector(BayesTools:::.bt_random_effect_latent_names(
      term, n_groups = 2L, n_columns = term$n_columns)))
    plan <- BayesTools:::.bt_JAGS_bridge_compile_random_value_plan(term, output$prior_list)
    expected <- if(structure == "id") c(3, 4, -3, 5.5) else c(1, 4, 1, 4)
    expect_false(is.null(plan))
    expect_equal(BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples, term, output$prior_list, value_plan = plan), expected)
    evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
      formula_list = list(mu = output$formula),
      formula_data_list = list(mu = output$data),
      formula_prior_list = list(mu = output$prior_list),
      formula_design_list = list(mu = output$formula_design),
      model_data = list())
    expect_equal(evaluator$parameters(samples, list())$mu, expected)
    expect_error(
      evaluator$parameters(samples[-1L], list()),
      "Bridge samples are missing standardized latent random effects for block 'group'.",
      fixed = TRUE)
  }
})

test_that("ordered formula parameters are reconstructed from bridge coordinates", {

  formula_data <- data.frame(
    f = ordered(
      c("low", "mid", "high", "low"),
      levels = c("low", "mid", "high")
    )
  )
  cases <- list(
    fixed_allocation = list(
      total = prior("normal", list(0, 1)),
      allocation = c(.25, .75),
      samples = c(mu_f_ordered_total = 2),
      coefficients = c(.5, 1.5),
      formula_value = c(0, .5, 2, 0)
    ),
    dirichlet_allocation = list(
      total = prior("normal", list(0, 1)),
      allocation = prior("dirichlet", list(alpha = c(2, 3))),
      samples = c(
        mu_f_ordered_total = 2,
        "prior_par_eta_mu_f_ordered_alloc_f_1[1]" = 1,
        "prior_par_eta_mu_f_ordered_alloc_f_1[2]" = 3
      ),
      coefficients = c(.5, 1.5),
      formula_value = c(0, .5, 2, 0)
    ),
    point_total = list(
      total = prior("point", list(4)),
      allocation = c(.25, .75),
      samples = numeric(),
      coefficients = c(1, 3),
      formula_value = c(0, 1, 4, 0)
    )
  )

  for(case_name in names(cases)){
    case <- cases[[case_name]]
    formula_output <- JAGS_formula(
      formula = ~ f,
      parameter = "mu",
      data = formula_data,
      prior_list = list(
        intercept = prior("point", list(0)),
        f = prior_ordered(
          total = case$total,
          allocation = case$allocation
        )
      )
    )
    formula_prior_list <- list(mu = formula_output$prior_list)
    formula_list <- list(mu = formula_output$formula)
    formula_data_list <- list(mu = formula_output$data)
    formula_design_list <- list(mu = formula_output$formula_design)
    coefficient_names <- BayesTools:::.JAGS_prior_factor_names(
      "mu_f",
      formula_output$prior_list$mu_f
    )

    expect_false(
      any(coefficient_names %in% names(case$samples)),
      info = case_name
    )

    public_prior_parameters <- JAGS_marglik_parameters(
      case$samples,
      formula_output$prior_list
    )
    compiled_prior_parameters <-
      BayesTools:::.bt_JAGS_bridge_compile_prior_list_evaluator(
        formula_output$prior_list
      )$parameters(case$samples)

    expect_equal(
      public_prior_parameters$mu_f,
      case$coefficients,
      tolerance = 1e-12,
      info = case_name
    )
    expect_equal(
      compiled_prior_parameters$mu_f,
      case$coefficients,
      tolerance = 1e-12,
      info = case_name
    )

    public_formula_parameters <- JAGS_marglik_parameters_formula(
      samples = case$samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list,
      model_data = list()
    )
    compiled_formula_parameters <-
      BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
        formula_list = formula_list,
        formula_data_list = formula_data_list,
        formula_prior_list = formula_prior_list,
        formula_design_list = formula_design_list,
        model_data = list()
      )$parameters(case$samples, list())

    expect_equal(
      public_formula_parameters$mu,
      case$formula_value,
      tolerance = 1e-12,
      info = case_name
    )
    expect_equal(
      compiled_formula_parameters$mu,
      case$formula_value,
      tolerance = 1e-12,
      info = case_name
    )
  }
})

test_that("formula parameter evaluators reject missing named multiply_by parameters", {

  formula_data <- data.frame(x = c(-1, 0, 2))
  formula_output <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = formula_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )
  )
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(mu_intercept = .2, mu_x = -.4)

  for(prior_name in c("mu_intercept", "mu_x")){
    formula_prior_list <- list(mu = formula_output$prior_list)
    multiplier_name <- paste0(prior_name, "_scale")
    attr(
      formula_prior_list$mu[[prior_name]],
      "multiply_by"
    ) <- multiplier_name

    compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      formula_design_list = formula_design_list,
      model_data = list()
    )
    expected_error <- paste0(
      "Formula prior 'multiply_by' parameter '",
      multiplier_name,
      "' is missing from 'prior_list_parameters'."
    )

    expect_error(
      compiled$parameters(samples, list()),
      expected_error,
      fixed = TRUE
    )
    expect_error(
      JAGS_marglik_parameters_formula(
        samples = samples,
        formula_list = formula_list,
        formula_data_list = formula_data_list,
        formula_prior_list = formula_prior_list,
        prior_list_parameters = list(),
        formula_design_list = formula_design_list,
        model_data = list()
      ),
      expected_error,
      fixed = TRUE
    )
    expect_error(
      JAGS_marglik_parameters_formula(
        samples = samples,
        formula_list = formula_list,
        formula_data_list = formula_data_list,
        formula_prior_list = formula_prior_list,
        prior_list_parameters = list()
      ),
      expected_error,
      fixed = TRUE
    )
  }
})

test_that("compiled formula parameter evaluator preserves log-intercept reconstruction", {

  formula_obj <- ~ 1 + x
  attr(formula_obj, "log(intercept)") <- TRUE
  formula_data <- data.frame(x = c(0, 1, 2))
  formula_output <- JAGS_formula(
    formula = formula_obj,
    parameter = "mu",
    data = formula_data,
    prior_list = list(
      intercept = prior("gamma", list(2, 1)),
      x         = prior("normal", list(0, 1))
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(mu_intercept = 1.25, mu_x = .4)

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )

  expect_equal(
    compiled$parameters(samples, list()),
    JAGS_marglik_parameters_formula(
      samples = samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list,
      model_data = list()
    ),
    tolerance = 1e-12
  )
})

test_that("formula replay rejects invalid log-intercept prior support", {
  formula_obj <- ~ 1
  attr(formula_obj, "log(intercept)") <- TRUE
  formula_prior_list <- list(
    mu = list(mu_intercept = prior("point", list(0)))
  )

  expect_error(
    JAGS_marglik_parameters_formula(
      samples = c(mu_intercept = 0),
      formula_list = list(mu = formula_obj),
      formula_data_list = list(mu = list(N_mu = 1)),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    "must have strictly positive support",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
      formula_list = list(mu = formula_obj),
      formula_data_list = list(mu = list(N_mu = 1)),
      formula_prior_list = formula_prior_list,
      formula_design_list = NULL,
      model_data = list()
    ),
    "must have strictly positive support",
    fixed = TRUE
  )
})

test_that("formula reconstruction rejects unknown priors instead of returning zero", {
  formula_output <- JAGS_formula(
    ~ x,
    "mu",
    data.frame(x = c(1, 2)),
    list(
      intercept = prior("point", list(1)),
      x = prior("normal", list(0, 1))
    )
  )
  unsupported_prior <- structure(
    list(distribution = "unsupported", prior_weights = 1),
    class = c("prior", "prior.unsupported")
  )
  formula_output$prior_list$mu_x <- unsupported_prior
  formula_prior_list <- list(mu = formula_output$prior_list)

  expect_error(
    JAGS_marglik_parameters_formula(
      samples = c(mu_x = 2),
      formula_list = list(mu = formula_output$formula),
      formula_data_list = list(mu = formula_output$data),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = list(mu = formula_output$formula_design)
    ),
    "Unsupported formula reconstruction prior for 'mu_x'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
      formula_list = list(mu = formula_output$formula),
      formula_data_list = list(mu = formula_output$data),
      formula_prior_list = formula_prior_list,
      formula_design_list = list(mu = formula_output$formula_design),
      model_data = list()
    ),
    "Unsupported formula reconstruction prior for 'mu_x'",
    fixed = TRUE
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples = c(mu_x = 2),
      formula_list = list(mu = formula_output$formula),
      formula_data_list = list(mu = formula_output$data),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    "Unsupported formula reconstruction prior for 'mu_x'",
    fixed = TRUE
  )
})

test_that("formula reconstruction rejects unsupported fixed-formula calls", {
  invalid_formula <- ~ I(x^2)
  formula_prior_list <- list(
    mu = list(mu_intercept = prior("point", list(0)))
  )

  expect_error(
    JAGS_marglik_parameters_formula(
      samples = numeric(),
      formula_list = list(mu = invalid_formula),
      formula_data_list = list(mu = list(N_mu = 2)),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    "Unsupported fixed-formula call 'I(x^2)'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
      formula_list = list(mu = invalid_formula),
      formula_data_list = list(mu = list(N_mu = 2)),
      formula_prior_list = formula_prior_list,
      formula_design_list = NULL,
      model_data = list()
    ),
    "Unsupported fixed-formula call 'I(x^2)'",
    fixed = TRUE
  )
})

test_that("compiled formula parameter evaluator matches legacy fallback reconstruction", {

  samples <- c(
    mu_intercept = 1,
    mu_x_data = 2
  )
  formula_data_list <- list(
    mu = list(
      N_mu = 2,
      mu_data_x_data = c(10, 20)
    )
  )
  formula_prior_list <- list(
    mu = list(
      mu_intercept = prior("normal", list(0, 1)),
      mu_x_data    = prior("normal", list(0, 1))
    )
  )
  formula_list <- list(mu = ~ 1 + x_data)

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = NULL,
    model_data = list()
  )

  expect_equal(
    compiled$parameters(samples, list()),
    JAGS_marglik_parameters_formula(
      samples = samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    tolerance = 1e-12
  )
})

test_that("compiled formula parameter evaluator still requires zero-multiplied samples", {

  formula_data <- data.frame(x = c(0, 1, 2))
  formula_output <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = formula_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  attr(formula_prior_list$mu$mu_x, "multiply_by") <- 0
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  missing_term_samples <- c(mu_intercept = .2)

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )

  expect_error(
    compiled$parameters(missing_term_samples, list()),
    "formula prior parameters"
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples = missing_term_samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list,
      model_data = list()
    ),
    "formula prior parameters"
  )
})

test_that("compiled formula source parameter overlay matches public helper", {

  samples <- c(theta = .1, phi = .2)
  prior_list_parameters <- list(theta = 1, tau = 2)
  formula_parameters <- list(mu = c(1, 2, 3), sigma = c(4, 5, 6))
  source_base <- BayesTools:::.bt_JAGS_bridge_formula_source_base_parameters(
    samples = samples,
    prior_list_parameters = prior_list_parameters
  )

  expect_identical(
    BayesTools:::.bt_JAGS_bridge_formula_source_parameters(
      source_base = source_base,
      formula_parameters = formula_parameters
    ),
    BayesTools:::.bt_JAGS_marglik_parameter_source_parameters(
      samples = samples,
      prior_list_parameters = prior_list_parameters,
      formula_parameters = formula_parameters
    )
  )
})

test_that("row-indexed source parameter shortcut reuses complete bridge lists", {

  posterior <- matrix(
    c(.1, .2),
    nrow = 1,
    dimnames = list(NULL, c("theta", "phi"))
  )
  complete_parameters <- list(theta = 1, phi = 2, mu = c(3, 4))
  partial_parameters <- list(theta = 1, mu = c(3, 4))

  expect_identical(
    BayesTools:::.bt_JAGS_marglik_row_indexed_source_parameters(
      posterior = posterior,
      draw = 1,
      parameters = complete_parameters
    ),
    complete_parameters
  )
  expect_identical(
    BayesTools:::.bt_JAGS_marglik_row_indexed_source_parameters(
      posterior = posterior,
      draw = 1,
      parameters = partial_parameters
    ),
    BayesTools:::.bt_parameter_source_draw_parameters(
      posterior = posterior,
      draw = 1,
      parameters = partial_parameters
    )
  )
})

test_that("compiled bridge allocation plans reuse cached Dirichlet draws", {

  formula_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    drug = factor(c("a", "b", "a", "b", "a", "b"))
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(
    "mu_intercept" = 10,
    "mu__xRE_ALLOCx_allocation__allocation_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__study_xRE_Zx[3,1]" = 0.3,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  cached_samples <- BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE)
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )

  expect_equal(
    compiled$parameters(cached_samples, list()),
    JAGS_marglik_parameters_formula(
      samples = samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list,
      model_data = list()
    ),
    tolerance = 1e-12
  )
  expect_equal(
    compiled$parameters(cached_samples, list())$mu,
    10 +
      c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
      sqrt(3) * c(1, 2, 1, 2, 1, 2),
    tolerance = 1e-12
  )

  posterior_row <- BayesTools:::.bt_JAGS_marglik_random_effect_posterior_row(cached_samples)
  cache <- attr(
    posterior_row,
    "BayesTools_random_effect_dirichlet_draw_cache",
    exact = TRUE
  )
  cache_key <- BayesTools:::.bt_random_effect_dirichlet_cache_key(
    "mu__xRE_ALLOCx_allocation__weight",
    2,
    "eta"
  )
  expect_true(is.environment(cache))
  expect_true(exists(cache_key, envir = cache, inherits = FALSE))

  prepared_design <- BayesTools:::.bt_JAGS_bridge_prepare_random_effect_allocation_design(
    formula_output$formula_design
  )
  factor_plan <- attr(
    prepared_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors,
    "BayesTools_random_effect_allocation_factor_plan",
    exact = TRUE
  )
  expect_length(factor_plan, 1L)
  expect_equal(factor_plan[[1]]$weight_name, "mu__xRE_ALLOCx_allocation__weight")

  bad_samples <- samples
  bad_samples[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]"]] <- 0
  bad_error <- tryCatch(
    compiled$parameters(
      BayesTools:::.bt_JAGS_bridge_cache_posterior_row(bad_samples, TRUE),
      list()
    ),
    error = function(e) e
  )
  expect_s3_class(bad_error, "BayesTools_marglik_out_of_support")
  expect_match(
    conditionMessage(bad_error),
    "Dirichlet allocation auxiliary samples"
  )
})

test_that("compiled bridge allocation cache preserves row-indexed source reconstruction", {

  formula_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    tau_factor = c(0.5, 0.75, 1, 1.25)
  )
  tau_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      data$tau_factor[seq_len(n_rows)]
    }
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(tau_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(
    "mu_intercept" = 10,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )
  parameters <- compiled$parameters(
    BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE),
    list()
  )
  expected_random <- formula_data$tau_factor * (
    c(0.1, 0.1, 0.2, 0.2) * sqrt(1 / 4) +
      c(1, 2, 1, 2) * sqrt(3 / 4)
  )

  expect_equal(
    parameters,
    JAGS_marglik_parameters_formula(
      samples = samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list,
      model_data = list()
    ),
    tolerance = 1e-12
  )
  expect_equal(parameters$mu, 10 + expected_random, tolerance = 1e-12)
})

test_that("compiled row sources receive natural formula-prior parameters", {

  formula_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    tau_factor = c(0.5, 0.75, 1, 1.25)
  )
  weight_name <- "mu__xRE_ALLOCx_allocation__weight"
  eta_name <- paste0("prior_par_eta_", weight_name, "[1]")
  tau_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters[[weight_name]][1] *
        parameters[[eta_name]] *
        data$tau_factor[seq_len(n_rows)]
    }
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(tau_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  samples <- c(
    "mu_intercept" = 10,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 6,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  samples[[weight_name]] <- 99
  formula_prior_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list
  )
  formula_prior_parameters <- formula_prior_evaluator$parameters(samples)
  expect_equal(
    formula_prior_parameters[[weight_name]],
    c(0.25, 0.75),
    tolerance = 1e-12
  )

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = list(mu = formula_output$formula),
    formula_data_list = list(mu = formula_data),
    formula_prior_list = formula_prior_list,
    formula_design_list = list(mu = formula_output$formula_design),
    model_data = list()
  )
  parameters <- compiled$parameters(
    BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE),
    list(),
    formula_prior_parameters
  )
  source_values <- 0.25 * 2 * formula_data$tau_factor
  expected_random <- source_values * (
    c(0.1, 0.1, 0.2, 0.2) * sqrt(0.25) +
      c(1, 2, 1, 2) * sqrt(0.75)
  )

  expect_equal(parameters$mu, 10 + expected_random, tolerance = 1e-12)
})

test_that("row-source callbacks reject sampled random formula dependencies", {

  formula_data <- data.frame(
    id = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  dependent_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      exp(parameters$log_sigma)
    }
  )
  mu_output <- JAGS_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd_source = random_sd_source(dependent_source))
    )
  )
  log_sigma_output <- JAGS_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    parameter = "log_sigma",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 1)))
    )
  )
  formula_list <- list(
    mu = mu_output$formula,
    log_sigma = log_sigma_output$formula
  )
  formula_data_list <- list(
    mu = mu_output$data,
    log_sigma = log_sigma_output$data
  )
  formula_prior_list <- list(
    mu = mu_output$prior_list,
    log_sigma = log_sigma_output$prior_list
  )
  formula_design_list <- list(
    mu = mu_output$formula_design,
    log_sigma = log_sigma_output$formula_design
  )
  samples <- c(
    "mu_intercept" = 0,
    "mu__xREx__id_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__id_xRE_Zx[2,1]" = 0.2,
    "log_sigma_intercept" = 0,
    "log_sigma__xREx__id_intercept" = 1,
    "log_sigma__xREx__id_xRE_Zx[1,1]" = 0.3,
    "log_sigma__xREx__id_xRE_Zx[2,1]" = 0.4
  )
  dependency_error <- "source 'tau\\[row\\]'.*formula parameter 'log_sigma'"

  expect_error(
    JAGS_marglik_parameters_formula(
      samples = samples,
      formula_list = formula_list,
      formula_data_list = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list,
      model_data = list()
    ),
    dependency_error,
    perl = TRUE
  )

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )
  expect_error(
    compiled$parameters(
      BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE),
      list()
    ),
    dependency_error,
    perl = TRUE
  )
})

test_that("bridge posterior-row cache preserves random-effect helper fallback shape", {

  samples <- c("z[1,1]" = .1, "z[1,2]" = .2)
  expected <- matrix(
    unname(samples),
    nrow = 1L,
    dimnames = list(NULL, names(samples))
  )

  expect_identical(
    BayesTools:::.bt_JAGS_marglik_random_effect_posterior_row(samples),
    expected
  )

  cached <- BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE)
  cached_row <- BayesTools:::.bt_JAGS_marglik_random_effect_posterior_row(cached)
  cache <- attr(
    cached_row,
    "BayesTools_random_effect_dirichlet_draw_cache",
    exact = TRUE
  )
  attr(cached_row, "BayesTools_random_effect_dirichlet_draw_cache") <- NULL
  expect_identical(
    cached_row,
    expected
  )
  expect_null(
    attr(
      BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, FALSE),
      "BayesTools_marglik_posterior_row",
      exact = TRUE
    )
  )
  expect_true(is.environment(cache))
})

test_that("bridge context exposes resolved formula allocation nodes", {

  formula_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    drug = factor(c("a", "b", "a", "b", "a", "b"))
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(
    "mu_intercept" = 10,
    "mu__xRE_ALLOCx_allocation__allocation_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__study_xRE_Zx[3,1]" = 0.3,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  samples <- BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE)
  formula_prior_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list
  )
  formula_parameter_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )
  prior_parameters <- list()
  formula_prior_parameters <- formula_prior_evaluator$parameters(samples)
  formula_parameters <- formula_parameter_evaluator$parameters(
    samples,
    prior_parameters
  )

  context <- BayesTools:::.bt_JAGS_bridge_context(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )

  expect_s3_class(context, "BayesTools_bridge_context")
  expect_true(is.matrix(context$state_matrix))
  expect_equal(
    colnames(context$state_matrix),
    names(context$state)
  )
  expect_true(is.environment(
    attr(
      context$state_matrix,
      "BayesTools_random_effect_dirichlet_draw_cache",
      exact = TRUE
    )
  ))
  expect_equal(
    context$formula_prior_parameters[["mu__xRE_ALLOCx_allocation__weight"]],
    c(1, 3) / 4,
    tolerance = 1e-12
  )
  expect_equal(
    context$nodes[["mu__xRE_ALLOCx_allocation__weight[1]"]],
    1 / 4,
    tolerance = 1e-12
  )
  expect_equal(
    context$nodes[["mu__xRE_ALLOCx_allocation__weight[2]"]],
    3 / 4,
    tolerance = 1e-12
  )
  expect_true(
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" %in%
      names(context$state)
  )
  expect_false(
    "mu__xRE_ALLOCx_allocation__weight[1]" %in%
      names(context$state)
  )
  weight_info <- context$node_info[
    context$node_info$name == "mu__xRE_ALLOCx_allocation__weight[1]",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(weight_info), 1L)
  expect_equal(weight_info$owner, "formula")
  expect_equal(weight_info$role, "random_allocation")
  expect_true(is.na(weight_info$block_name))
  expect_equal(
    context$random$mu$study$block_name,
    "study"
  )
  expect_equal(
    context$random$mu$study$dimensions$n_groups,
    3L
  )
  expect_equal(
    unname(context$random$mu$study$allocation$weights[["mu__xRE_ALLOCx_allocation__weight"]]),
    c(1, 3) / 4,
    tolerance = 1e-12
  )
  expect_equal(
    context$random$mu$study$scale$type,
    "column"
  )
  expect_equal(
    unname(context$random$mu$study$scale$column_sd),
    1,
    tolerance = 1e-12
  )
  expect_equal(
    unname(context$random$mu$drug$scale$column_sd),
    sqrt(3),
    tolerance = 1e-12
  )
  expect_equal(
    context$random$mu$study$correlation$matrix,
    matrix(1, 1L, 1L)
  )
  expect_equal(
    context$formula_parameters$mu,
    formula_parameters$mu,
    tolerance = 1e-12
  )

  context_arguments <- list(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters
  )
  full_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_context_evaluator(
    mode = TRUE,
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )
  compiled_full <- do.call(full_evaluator$context, context_arguments)
  expect_equal(compiled_full, context, tolerance = 0)

  nodes_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_context_evaluator(
    mode = "nodes",
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )
  nodes_context <- do.call(nodes_evaluator$context, context_arguments)
  expect_s3_class(nodes_context, "BayesTools_bridge_nodes_context")
  expect_named(nodes_context, "nodes")
  expect_identical(nodes_context$nodes, context$nodes)

  selected_names <- c(
    formula_output$formula_design$random_effects[[2L]]$sd_parameter_names,
    formula_output$formula_design$random_effects[[1L]]$sd_parameter_names
  )
  selected_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_context_evaluator(
    mode = "nodes",
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list(),
    node_names = selected_names
  )
  selected_context <- do.call(selected_evaluator$context, context_arguments)
  expect_identical(
    selected_context$nodes,
    context$nodes[selected_names]
  )

  missing_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_context_evaluator(
    mode = "nodes",
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list(),
    node_names = "absent_node"
  )
  expect_error(
    do.call(missing_evaluator$context, context_arguments),
    "Requested bridge context node(s) are unavailable: absent_node",
    fixed = TRUE
  )
})

test_that("marginal bridge context skips an explicitly empty node selection", {

  random_evaluator <- list(
    nodes = function(...) stop("random nodes should not be evaluated")
  )
  marginal_evaluator <- list(
    covariance = function(...) list(mu = list(representation = "factor_state"))
  )

  context <- BayesTools:::.bt_JAGS_bridge_marginal_context(
    samples = c(theta = 0),
    prior_parameters = list(),
    formula_prior_parameters = list(),
    formula_parameters = list(),
    add_parameters = NULL,
    random_context_evaluator = random_evaluator,
    marginal_random_evaluator = marginal_evaluator,
    node_names = character(),
    random_only = TRUE
  )

  expect_s3_class(context, "BayesTools_bridge_marginal_context")
  expect_identical(context$nodes, numeric())
  expect_identical(
    context$marginalized_random$mu$representation,
    "factor_state"
  )
})

test_that("bridge context exposes marginalized random blocks without latent draws", {

  formula_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"))
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | estimate, name = "estimate", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "total_re",
        terms = c(study = "study", estimate = "estimate"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    ),
    random_effects_compile = random_effects_compile(
      marginalized = "estimate"
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_output$data)
  formula_design_list <- list(mu = formula_output$formula_design)
  random_terms <- formula_output$formula_design$random_effects
  study_term <- random_terms[[1L]]
  estimate_term <- random_terms[[2L]]
  z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = study_term,
    n_groups = study_term$n_groups,
    n_columns = study_term$n_columns
  ))
  weight_name <- "mu__xRE_ALLOCx_total_re__weight"
  samples <- c(
    "mu_intercept" = 10,
    "mu__xRE_ALLOCx_total_re__allocation_sd" = 2,
    stats::setNames(
      c(1, 3),
      paste0(
        BayesTools:::.JAGS_prior_dirichlet_eta_name(weight_name),
        "[",
        1:2,
        "]"
      )
    ),
    stats::setNames(c(.1, .2), z_names)
  )
  samples <- BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE)
  formula_prior_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list
  )
  formula_parameter_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )
  prior_parameters <- list()
  formula_prior_parameters <- formula_prior_evaluator$parameters(samples)
  formula_parameters <- formula_parameter_evaluator$parameters(
    samples,
    prior_parameters
  )

  context <- BayesTools:::.bt_JAGS_bridge_context(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )

  expect_equal(context$random$mu$study$compile_mode, "sampled")
  expect_equal(context$random$mu$estimate$compile_mode, "marginalized")
  expect_true(is.matrix(context$random$mu$study$latent))
  expect_null(context$random$mu$estimate$latent)
  expect_false(any(grepl(
    "mu__xREx__estimate_xRE_Zx",
    names(context$state),
    fixed = TRUE
  )))
  expect_equal(
    unname(context$random$mu$estimate$allocation$weights[[weight_name]]),
    c(1, 3) / 4,
    tolerance = 1e-12
  )
  expect_equal(
    unname(context$random$mu$estimate$scale$column_sd),
    sqrt(3),
    tolerance = 1e-12
  )
  expect_equal(
    context$random$mu$estimate$covariance,
    matrix(3, 1L, 1L),
    tolerance = 1e-12
  )
  expect_true(estimate_term$sd_parameter_names %in% names(context$nodes))
  expect_equal(
    context$nodes[[estimate_term$sd_parameter_names]],
    sqrt(3),
    tolerance = 1e-12
  )

  nodes_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_context_evaluator(
    mode = "nodes",
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )
  nodes_context <- nodes_evaluator$context(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters
  )
  expect_identical(nodes_context$nodes, context$nodes)
})

test_that("compiled row sources retain per-draw validation", {

  source <- parameter_source(
    "tau", shape = "row",
    values = function(parameters, data, n_rows){

      parameters$value
    }
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("s1", "s1"))),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(sd_source = random_sd_source(source))
    )
  )
  evaluator <- BayesTools:::.bt_JAGS_bridge_compile_row_source_evaluator(
    formula_output$formula_design$random_effects[[1L]],
    n_rows = 2L,
    context = "Bridge context"
  )
  posterior <- matrix(numeric(), nrow = 1L, ncol = 0L)
  expect_equal(
    evaluator(posterior, parameters = list(value = c(0, 2))),
    matrix(c(0, 2), nrow = 1L, dimnames = list(NULL, c("tau[1]", "tau[2]")))
  )
  expect_error(
    evaluator(posterior, parameters = list(value = "invalid")),
    "Bridge context for source 'tau[row]' must return a numeric vector.",
    fixed = TRUE
  )
  expect_error(
    evaluator(posterior, parameters = list(value = 1)),
    "Bridge context for source 'tau[row]' must return a numeric vector of length 2.",
    fixed = TRUE
  )
  expect_error(
    evaluator(posterior, parameters = list(value = c(1, NA_real_))),
    "Bridge context for source 'tau[row]' returned missing values.",
    fixed = TRUE
  )
})


test_that("bridge context exposes row-indexed external SD source nodes", {

  formula_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    tau_factor = c(0.5, 0.75, 1, 1.25)
  )
  tau_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      data$tau_factor[seq_len(n_rows)]
    }
  )
  formula_output <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(tau_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_list <- list(mu = formula_output$formula)
  formula_data_list <- list(mu = formula_data)
  formula_design_list <- list(mu = formula_output$formula_design)
  samples <- c(
    "mu_intercept" = 10,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  samples <- BayesTools:::.bt_JAGS_bridge_cache_posterior_row(samples, TRUE)
  formula_prior_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list
  )
  formula_parameter_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_formula_parameter_evaluator(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_design_list = formula_design_list,
    model_data = list()
  )
  prior_parameters <- list()
  formula_prior_parameters <- formula_prior_evaluator$parameters(samples)
  formula_parameters <- formula_parameter_evaluator$parameters(
    samples,
    prior_parameters
  )

  context <- BayesTools:::.bt_JAGS_bridge_context(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )

  expect_equal(
    unname(context$random$mu$study$scale$row_sd_source),
    formula_data$tau_factor,
    tolerance = 1e-12
  )
  expect_equal(
    context$random$mu$study$scale$type,
    "row_indexed"
  )
  expect_equal(
    context$nodes[names(context$random$mu$study$scale$row_sd_source)],
    context$random$mu$study$scale$row_sd_source,
    tolerance = 1e-12
  )
  row_sd_info <- context$node_info[
    context$node_info$name == names(context$random$mu$study$scale$row_sd_source)[1L],
    ,
    drop = FALSE
  ]
  expect_equal(row_sd_info$role, "random_effect")
  expect_true(is.na(row_sd_info$block_name))
  expect_equal(
    unname(context$random$mu$study$allocation$weights[["mu__xRE_ALLOCx_allocation__weight"]]),
    c(1, 3) / 4,
    tolerance = 1e-12
  )
  expect_equal(
    context$random$mu$study$allocation$block_multiplier,
    sqrt(1 / 4),
    tolerance = 1e-12
  )
  expect_equal(
    context$random$mu$drug$allocation$block_multiplier,
    sqrt(3 / 4),
    tolerance = 1e-12
  )

  nodes_evaluator <- BayesTools:::.bt_JAGS_bridge_compile_context_evaluator(
    mode = "nodes",
    add_parameters = NULL,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = list()
  )
  nodes_context <- nodes_evaluator$context(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters
  )
  expect_identical(nodes_context$nodes, context$nodes)
})

test_that("compiled bridge prior evaluator rejects unsupported mixtures at setup", {

  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_prior_list_evaluator(
      list(theta = prior_mixture(list(
        prior("normal", list(0, 1)),
        prior("normal", list(1, 1))
      )))
    ),
    "mixture"
  )
})


test_that("bridge posterior resolves only complete singleton coordinate families", {

  raw <- matrix(c(-.3, .4, 1.1), ncol = 1L, dimnames = list(NULL, "z"))
  for(coordinate in c("z", "z[1]", "z[1,1]")){
    value <- JAGS_bridgesampling_posterior(
      posterior = raw, prior_list = list(), add_parameters = coordinate,
      add_bounds = list(lb = stats::setNames(-Inf, coordinate),
                        ub = stats::setNames(Inf, coordinate))
    )
    expect_identical(colnames(value), coordinate)
    expect_identical(unname(value[, 1L]), unname(raw[, 1L]))
    expect_identical(names(attr(value, "lb")), coordinate)
  }
  expect_identical(colnames(raw), "z")
  indexed <- raw
  colnames(indexed) <- "z[1,1]"
  value <- JAGS_bridgesampling_posterior(
    posterior = indexed, prior_list = list(), add_parameters = "z",
    add_bounds = list(lb = c(z = -Inf), ub = c(z = Inf))
  )
  expect_identical(colnames(value), "z")
  expect_identical(unname(value[, 1L]), unname(raw[, 1L]))
  # Neither a partly monitored array nor a competing indexed family is scalar.
  for(posterior in list(cbind(raw, "z[2,1]" = 0),
                       stats::setNames(data.frame(z = raw[, 1L]), "z[2,1]"))){
    expect_error(JAGS_bridgesampling_posterior(
      posterior = as.matrix(posterior), prior_list = list(),
      add_parameters = "z[1,1]",
      add_bounds = list(lb = c("z[1,1]" = -Inf), ub = c("z[1,1]" = Inf))
    ), "'posterior' does not contain all of the parameters corresponding to the 'prior_list' and the 'add_parameter' argument.",
    fixed = TRUE)
  }
})
