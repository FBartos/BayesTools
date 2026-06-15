skip_if_not_test_profile("unit")

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
    inv_sigma = .5,
    "v[1]" = -.1,
    "v[2]" = .3,
    "prior_par_eta_w[1]" = 1.2,
    "prior_par_eta_w[2]" = 2.4,
    "inv_theta[1]" = 2,
    "inv_theta[2]" = 3,
    PET = .4,
    PEESE = .6,
    "eta[1]" = 1.5,
    "eta[2]" = 2.5,
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
})

test_that("compiled bridge prior evaluators preserve auxiliary support behavior", {

  prior_list <- list(
    sigma = prior("invgamma", list(3, 2), list(1, 3)),
    w     = prior("dirichlet", list(alpha = c(2, 3)))
  )
  samples <- c(
    inv_sigma = .5,
    "prior_par_eta_w[1]" = 1.2,
    "prior_par_eta_w[2]" = 2.4
  )
  compiled <- BayesTools:::.bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)

  bad_invgamma <- samples
  bad_invgamma[["inv_sigma"]] <- 0
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
    inv_mu_x = .5,
    sigma_intercept = .8
  )

  compiled <- BayesTools:::.bt_JAGS_bridge_compile_formula_prior_evaluator(formula_prior_list)

  expect_equal(
    compiled$log_prior(samples),
    JAGS_marglik_priors_formula(samples, formula_prior_list),
    tolerance = 1e-12
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
    BayesTools:::.bt_JAGS_marglik_random_effect_prior(
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
    inv_mu_x = .5,
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
      allocation = random_variance_allocation(
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
    "mu__xRE_ALLOCx_allocation__total_sd" = 2,
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
      allocation = random_variance_allocation(
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
