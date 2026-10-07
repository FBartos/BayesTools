skip_if_not_test_profile("unit")

.natural_formula_case <- function(multiplier = "mu_intercept",
                                   intercept = prior("point", list(2)),
                                   ordinary = list(), model_data = NULL,
                                   random = FALSE, sd_source = "tau2"){
  data <- data.frame(x = c(-2, -1, 1, 2), id = factor(c("a", "a", "b", "b")))
  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- multiplier
  priors <- list(intercept = intercept, x = slope)
  formula <- if(random) ~ x + diag(1 | id) else ~ x
  prior_random <- if(random) prior_random(id = random_block(sd_source = random_sd_source(sd_source))) else NULL
  compiled <- JAGS_formula(formula, "mu", data, priors, prior_random = prior_random)
  row <- c(mu_x = .5)
  if(!is.prior.point(intercept)) row <- c(mu_intercept = 2, row)
  if(random){
    term <- compiled$formula_design$random_effects[[1L]]
    row <- c(row, setNames(c(1, -1), as.vector(.bt_random_effect_latent_names(term, term$n_groups, term$n_columns))))
  }
  draws <- matrix(rep(row, each = 24), 24, dimnames = list(NULL, names(row)))
  # Vary every independent coordinate separately to preserve the bridge rank.
  for(column in seq_len(ncol(draws))) draws[, column] <- draws[, column] + (seq_len(24) / 100)^column
  attr(draws, "prior_list") <- ordinary
  carrier <- JAGS_formula_draws(draws, formula, "mu", data, priors,
    prior_random = prior_random, model_data = model_data)
  fit <- coda::mcmc.list(coda::mcmc(as.matrix(carrier)))
  class(fit) <- c("BayesTools_fit", class(fit))
  for(name in c("prior_list", "formula_design", "formula_scale")) attr(fit, name) <- attr(carrier, name)
  fit <- .bt_attach_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)
  list(data = data, priors = priors, formula = formula, row = row, carrier = carrier,
    fit = fit, design = attr(carrier, "formula_design")$mu,
    formula_priors = attr(carrier, "formula_design")$mu$prior_list, ordinary = ordinary)
}

test_that("fixed reconstruction receives formula, extra, and retained scalar multipliers", {
  for(multiplier in c("mu_intercept", "extra", "sigma")){
    case <- .natural_formula_case(multiplier, model_data = if(multiplier == "sigma") list(sigma = 3) else NULL)
    row <- case$row
    if(multiplier == "extra") row <- c(row, extra = 4)
    expected <- 2 + switch(multiplier, mu_intercept = 2, extra = 4, sigma = 3) * .5 * case$data$x
    inputs <- list(formula_list = list(mu = case$formula), formula_data_list = list(mu = case$data),
      formula_prior_list = list(mu = case$formula_priors), formula_design_list = list(mu = case$design),
      model_data = list(sigma = 999))
    public <- do.call(JAGS_marglik_parameters_formula, c(list(samples = row, prior_list_parameters = list()), inputs))
    compiled <- do.call(.bt_JAGS_bridge_compile_formula_parameter_evaluator, inputs)
    expect_equal(public$mu, expected, tolerance = 0)
    expect_equal(compiled$parameters(row, list())$mu, expected, tolerance = 0)
  }
  case <- .natural_formula_case()
  expected <- 2 + 2 * outer(case$data$x, as.matrix(case$carrier)[, "mu_x"])
  expect_equal(unname(JAGS_evaluate_formula(case$carrier, parameter = "mu")), expected)
  case <- .natural_formula_case(intercept = prior("invgamma", list(2, 1)))
  row <- c(mu_intercept = 2, inv_mu_intercept = .5, mu_x = .25)
  decoded <- .bt_JAGS_bridge_compile_formula_prior_evaluator(list(mu = case$formula_priors))$parameters(row)
  expect_equal(decoded$mu_intercept, 2)
  compiled <- .bt_JAGS_bridge_compile_formula_parameter_evaluator(list(mu = case$formula), list(mu = case$data),
    list(mu = case$formula_priors), list(mu = case$design), list())
  expect_equal(compiled$parameters(row, list(), decoded)$mu, 2 + .5 * case$data$x)
  aliased <- row
  aliased[["mu_intercept"]] <- .5
  expect_equal(compiled$parameters(aliased, list(), decoded)$mu, 2 + .5 * case$data$x)
  expect_error(compiled$parameters(row, list(), list(mu_intercept = NA_real_)), class = "BayesTools_formula_transform_unavailable")
})

test_that("sampled scalar SD reconstruction forwards natural point state and cached rows", {
  case <- .natural_formula_case(random = TRUE, ordinary = list(tau2 = prior("point", list(.4))))
  expected <- 2 + case$data$x + c(.4, .4, -.4, -.4)
  inputs <- list(formula_list = list(mu = case$formula), formula_data_list = list(mu = case$data),
    formula_prior_list = list(mu = case$formula_priors), formula_design_list = list(mu = case$design), model_data = list())
  compiled <- do.call(.bt_JAGS_bridge_compile_formula_parameter_evaluator, inputs)
  cached <- .bt_JAGS_bridge_cache_posterior_row(case$row, TRUE)
  for(row in list(case$row, cached, c(case$row, tau2 = 99))){
    expect_equal(do.call(JAGS_marglik_parameters_formula, c(list(samples = row, prior_list_parameters = list(tau2 = .4)), inputs))$mu, expected)
    expect_equal(compiled$parameters(row, list(tau2 = .4))$mu, expected)
  }
  term <- case$design$random_effects[[1L]]
  sd <- .bt_JAGS_bridge_compile_random_sd_evaluator(term, case$formula_priors)
  expect_equal(sd$values(cached, list(tau2 = .4)), .4)
  expect_equal(.bt_JAGS_marglik_random_effect_sd_values(cached, term, case$formula_priors, list(tau2 = .4)), .4)
  for(value in list(NA_real_, c(.4, .5), "0.4")){
    expect_error(sd$values(c(case$row, tau2 = 99), list(tau2 = value)), class = "BayesTools_formula_transform_unavailable")
  }
  expect_error(compiled$parameters(case$row, list(tau2 = .7)), class = "BayesTools_formula_transform_unavailable")
  data_case <- .natural_formula_case(random = TRUE, sd_source = "sigma", model_data = list(sigma = 3))
  evaluator <- .bt_JAGS_bridge_compile_formula_parameter_evaluator(list(mu = data_case$formula),
    list(mu = data_case$data), list(mu = data_case$formula_priors), list(mu = data_case$design), list(sigma = 99))
  expect_equal(evaluator$parameters(c(data_case$row, sigma = 99), list())$mu,
    2 + data_case$data$x + c(3, 3, -3, -3))
})

test_that("public bridge callbacks receive natural state in every context", {
  callback_count <- 0L
  testthat::local_mocked_bindings(bridge_sampler = function(...){
    arguments <- list(...)
    names <- c("data", "bridge_prior_evaluator", "bridge_formula_prior_evaluator",
      "bridge_formula_random_prior_evaluator", "bridge_formula_parameter_evaluator", "add_parameters",
      "fixed_random_latent", "bridge_context", "bridge_context_evaluator")
    logml <- do.call(arguments$log_posterior,
      c(list(samples.row = arguments$samples[1L, ]), arguments[names]))
    callback_count <<- callback_count + 1L
    structure(list(logml = logml, niter = 1L, mcse_logml = 0, method = "normal"), class = "bridge")
  }, .package = "bridgesampling")
  for(random in c(FALSE, TRUE)){
    case <- .natural_formula_case(random = random,
      ordinary = if(random) list(tau2 = prior("point", list(.4))) else list())
    row <- as.matrix(case$carrier)[1, ]
    expected <- 2 + 2 * row[["mu_x"]] * case$data$x
    if(random){
      latent_names <- as.vector(.bt_random_effect_latent_names(case$design$random_effects[[1]], 2, 1))
      expected <- expected + .4 * row[latent_names][c(1, 1, 2, 2)]
    }
    for(mode in c("none", "full", "nodes", "marginal")){
      value <- JAGS_bridgesampling(case$fit, data = list(y = rep(0, 4)), bridge_context = mode,
        log_posterior = function(parameters, data, ...){
          expect_equal(unname(parameters$mu), unname(expected))
          0
        })
      expect_s3_class(value, "BayesTools_marglik")
    }
  }
  expect_identical(callback_count, 8L)
})

test_that("structured scalar SD state agrees with the declared formula-owned control", {
  data <- data.frame(x = c(-1, 0, 1, -1, 0, 1),
    id = factor(rep(c("a", "b"), each = 3)), t = factor(rep(1:3, 2)))
  formulas <- list(~ 1 + id(1 + x | id), ~ 1 + diag(1 + x | id), ~ 1 + ar1(t | id))
  for(i in seq_along(formulas)){
    formula <- formulas[[i]]
    source_block <- if(i == 3L) random_block(sd_source = random_sd_source("tau2"), cor = prior("point", list(.25))) else random_block(sd_source = random_sd_source("tau2"))
    control_block <- if(i == 3L) random_block(sd = prior("point", list(.4)), cor = prior("point", list(.25))) else random_block(sd = prior("point", list(.4)))
    external <- JAGS_formula(formula, "mu", data, list(intercept = prior("point", list(2))),
      prior_random = prior_random(id = source_block))
    control <- JAGS_formula(formula, "mu", data, list(intercept = prior("point", list(2))),
      prior_random = prior_random(id = control_block))
    term <- external$formula_design$random_effects[[1L]]
    latent_names <- if(inherits(term$latent_layout, "BayesTools_random_effect_structured_local_layout")){
      term$latent_layout$node_names
    }else as.vector(.bt_random_effect_latent_names(term, term$n_groups, term$n_columns))
    row <- setNames(seq_along(latent_names) / 10, latent_names)
    evaluate <- function(compiled, parameters){
      inputs <- list(formula_list = list(mu = formula), formula_data_list = list(mu = data),
        formula_prior_list = list(mu = compiled$prior_list), formula_design_list = list(mu = compiled$formula_design), model_data = list())
      public <- do.call(JAGS_marglik_parameters_formula, c(list(samples = row, prior_list_parameters = parameters), inputs))$mu
      plan <- do.call(.bt_JAGS_bridge_compile_formula_parameter_evaluator, inputs)
      expect_equal(plan$parameters(row, parameters)$mu, public, tolerance = 1e-14)
      public
    }
    expect_equal(evaluate(external, list(tau2 = .4)), evaluate(control, list()), tolerance = 1e-14)
  }
})

test_that("marginal batch and compact SD routes respect their explicit natural state", {
  case <- .natural_formula_case(random = TRUE, ordinary = list(tau2 = prior("point", list(.4))))
  block <- case$design$random_effects[[1L]]$block_name
  request <- .bt_JAGS_bridge_marginal_random_spec(list(mu = case$design),
    list(mu = list(blocks = block, row_blocks = list(1:2, 3:4), factor_state = TRUE)), "marginal")
  evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(list(mu = case$design), request,
    list(mu = case$data), list(mu = case$formula_priors), list(), posterior_names = names(case$row))
  posterior <- matrix(rep(case$row, 2), nrow = 2, byrow = TRUE, dimnames = list(NULL, names(case$row)))
  expect_equal(evaluator$coefficient_scales(posterior, "mu", block, parameters = list(tau2 = .4)), matrix(.4, 2, 1))
  components <- evaluator$factor_components(posterior, parameters = list(tau2 = .4))
  expect_equal(components$mu[[block]]$coefficient_scale, matrix(.4, 2, 1))
  expect_false(is.null(evaluator$factor_states(posterior, parameters = list(tau2 = .4))))
  for(value in list(NA_real_, "0.4", c(.4, .5, .6))){
    expect_error(evaluator$coefficient_scales(posterior, "mu", block, parameters = list(tau2 = value)), class = "BayesTools_formula_transform_unavailable")
  }
})
