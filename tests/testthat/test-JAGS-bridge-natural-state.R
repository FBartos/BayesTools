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
