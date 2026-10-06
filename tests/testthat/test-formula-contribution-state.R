skip_if_not_test_profile("unit")

.formula_state_test_fit <- function(multiplier = "sigma", mean = 20, sd = 10,
                                    intercept = prior("point", list(5)),
                                    slope = prior("point", list(2)),
                                    extra_priors = list(sigma = prior("normal", list(0, 1))),
                                    values = cbind(mu_intercept = c(5, 5), mu_x = c(2, 2), sigma = c(2, 4)),
                                    log_intercept = FALSE){

  attr(slope, "multiply_by") <- multiplier
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- log_intercept
  compiled <- JAGS_formula(formula, "mu", data.frame(x = mean + c(-sd, 0, sd)),
    list(intercept = intercept, x = slope), formula_scale = TRUE)
  .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(values)),
    c(compiled$prior_list, extra_priors), list(mu = compiled$formula_design),
    list(mu = compiled$formula_scale))
}

test_that("original coefficients read immutable fitted multiplier states", {

  fit <- .formula_state_test_fit()
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  expect_identical(descriptor$schema_version, 3L)
  expect_true(all(is.na(descriptor$matrix["mu_intercept", ])))
  expect_identical(descriptor$matrix["mu_x", ], c(mu_intercept = 0, mu_x = .1))
  expect_identical(descriptor$basis_matrix["mu_intercept", ], c(mu_intercept = 1, mu_x = -2))
  expect_identical(descriptor$prior_recipes$mu_intercept$type, "contribution_affine")
  expect_identical(descriptor$prior_recipes$mu_x$type, "raw_affine")
  corrected <- transform_scale_samples(fit)
  expect_identical(unname(corrected[, "mu_intercept"]), c(-3, -11))
  expect_identical(unname(corrected[, "mu_x"]), rep(.2, 2))
  x <- c(0, 20, 37)
  fitted <- outer(c(2, 4) * 2, (x - 20) / 10) + 5
  original <- outer(corrected[, "mu_x"] * c(2, 4), x) + corrected[, "mu_intercept"]
  expect_equal(original, fitted, tolerance = 1e-14)
  expect_identical(JAGS_formula_coefficient_transform_schema()$field, names(descriptor))
})

test_that("same multiplier cancellation at zero and selected scope are exact", {

  fit <- .formula_state_test_fit(values = cbind(mu_intercept = 5, mu_x = 2, sigma = c(0, 2)))
  corrected <- transform_scale_samples(fit)
  expect_identical(unname(corrected[, "mu_x"]), rep(.2, 2))
  expect_identical(unname(corrected[, "mu_intercept"]), c(5, -3))
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  supplied <- cbind(mu_x = c(2, 2))
  selected <- .bt_apply_formula_coefficient_transform(supplied, descriptor, "mu_x")
  expect_identical(unname(selected[, "mu_x"]), rep(.2, 2))
  expect_error(.bt_apply_formula_coefficient_transform(supplied, descriptor, "mu_intercept"),
    class = "BayesTools_formula_transform_unavailable")
})

test_that("logged intercepts and numeric overrides retain fitted declarations", {

  fit <- .formula_state_test_fit(intercept = prior("point", list(exp(5))), log_intercept = TRUE,
    values = cbind(mu_intercept = c(exp(5), exp(5)), mu_x = c(2, 2), sigma = c(2, 4)))
  transformed <- transform_scale_samples(fit)
  expect_equal(unname(transformed[, "mu_intercept"]), exp(c(-3, -11)), tolerance = 1e-14)
  mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  row <- marginal_posterior(mixed, "mu_intercept", ~ x, list(x = 20))$intercept
  expect_equal(as.numeric(row), rep(5, 2), tolerance = 1e-14)
  override <- list(mu = list(mu_x = list(mean = 0, sd = 2)))
  altered <- transform_scale_samples(.formula_state_test_fit(), formula_scale = override)
  expect_identical(unname(altered[, "mu_intercept"]), c(5, 5))
  expect_identical(unname(altered[, "mu_x"]), c(1, 1))
  missing <- .formula_state_test_fit(values = cbind(mu_intercept = c(5, 5), mu_x = c(2, 2)))
  expect_error(transform_scale_samples(missing), class = "BayesTools_formula_transform_unavailable")
  descriptor <- JAGS_formula_coefficient_transform(.formula_state_test_fit(), "mu")
  for(value in c(NA_real_, Inf)){
    bad <- cbind(mu_intercept = c(5, 5), mu_x = c(2, 2), sigma = c(2, value))
    expect_error(.bt_apply_formula_coefficient_transform(bad, descriptor, "mu_intercept"),
      class = "BayesTools_formula_transform_unavailable")
  }
})

test_that("distinct zero denominators refuse only requested numeric targets", {

  priors <- list(intercept = prior("point", list(5)), x = prior("point", list(2)),
    z = prior("point", list(3)), `x:z` = prior("point", list(4)))
  attr(priors$x, "multiply_by") <- "u"
  attr(priors$z, "multiply_by") <- attr(priors$`x:z`, "multiply_by") <- "v"
  compiled <- JAGS_formula(~ x * z, "mu", data.frame(x = c(10, 20, 30), z = c(8, 10, 12)),
    priors, formula_scale = TRUE)
  ordinary <- list(u = prior_spike_and_slab(prior("normal", list(0, 1))),
    v = prior("normal", list(0, 1)))
  sources <- cbind(mu_intercept = c(5, 5), mu_x = c(2, 2), mu_z = c(3, 3),
    mu_x__xXx__z = c(4, 4), u = c(0, 2), v = c(1, 1), u_indicator = c(0, 1))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(sources)),
    c(compiled$prior_list, ordinary), list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
  descriptor <- JAGS_formula_coefficient_transform(fit, "mu")
  condition <- tryCatch(.bt_apply_formula_coefficient_transform(sources, descriptor, "mu_x"), error = identity)
  expect_s3_class(condition, "BayesTools_formula_transform_unavailable")
  expect_false(inherits(condition, "BayesTools_formula_measure_unavailable"))
  expect_identical(condition$reason, "zero_multiplier_denominator")
  sibling <- as_mixed_posteriors(fit, "mu_z", transform_scaled = TRUE)$mu_z
  expect_identical(as.numeric(sibling), c(-2.5, -2.5))
  expect_identical(.posterior_atoms_get(sibling)$mass, 1)
  sources[, "u"] <- c(2, 4)
  sources[, "u_indicator"] <- 1
  expect_equal(as.numeric(.bt_apply_formula_coefficient_transform(sources, descriptor, "mu_x")[, "mu_x"]),
    c(-.8, -.3), tolerance = 1e-14)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(sources)),
    c(compiled$prior_list, ordinary), list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
  selected <- as_mixed_posteriors(fit, c("mu_x", "mu_z"), transform_scaled = TRUE)
  expect_equal(as.numeric(selected$mu_x), c(-.8, -.3), tolerance = 1e-14)
  expect_identical(posterior_metadata(selected$mu_x, "measure_unavailable")$measure,
    c("prior_density", "atoms", "support"))
  expect_error(JAGS_formula_prior_density(fit, "mu", "mu_x"), class = "BayesTools_formula_prior_density_unavailable")
  expect_error(Savage_Dickey_BF(marginal_posterior(selected, "mu_x", use_formula = FALSE), 0, silent = TRUE),
    class = "BayesTools_formula_measure_unavailable")
  expect_equal(as.numeric(marginal_posterior(selected, "mu_x", use_formula = FALSE)), c(-.8, -.3), tolerance = 1e-14)
  expect_identical(JAGS_formula_prior_density(fit, "mu", "mu_z")$points$x, -2.5)
})

test_that("declared zero contributions do not require irrelevant states", {

  fit <- .formula_state_test_fit(slope = prior("point", list(0)),
    values = cbind(mu_intercept = c(5, 5), mu_x = c(0, 0)))
  original <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  expect_identical(as.numeric(original$mu_intercept), c(5, 5))
  expect_identical(as.numeric(original$mu_x), c(0, 0))
  row <- marginal_posterior(original, "mu_intercept", ~ x, list(x = 1), prior_samples = TRUE)$intercept
  expect_identical(as.numeric(row), c(5, 5))
  expect_identical(.posterior_atoms_get(row)$mass, 1)
  fit <- .formula_state_test_fit(multiplier = 0, slope = prior("normal", list(0, 1)),
    values = cbind(mu_intercept = c(5, 5), mu_x = c(2, 3)), extra_priors = list())
  transformed <- transform_scale_samples(fit)
  expect_equal(unname(transformed[, "mu_x"]), c(.2, .3), tolerance = 1e-14)
  expect_identical(unname(transformed[, "mu_intercept"]), c(5, 5))
})

test_that("current draws-only data constants finalize both carriers", {

  slope <- prior("point", list(2))
  attr(slope, "multiply_by") <- "sigma"
  draws <- JAGS_formula_draws(cbind(mu_intercept = c(5, 5), mu_x = c(2, 2)), ~ x, "mu",
    data.frame(x = c(10, 20, 30)), list(intercept = prior("point", list(5)), x = slope),
    formula_scale = TRUE, model_data = list(sigma = 3))
  design <- attr(draws, "formula_design")$mu
  scale <- attr(draws, "formula_scale")$mu
  expect_identical(attr(scale, "unscale_design"), attr(design$formula_scale, "unscale_design"))
  expect_identical(attr(scale, "unscale_design")$state_constants, c(sigma = 3))
  expect_identical(attr(scale, "unscale_design")$owner_scope, "fit")
  expect_true(attr(scale, "unscale_design")$complete)
  descriptor <- .bt_formula_coefficient_transform(c("mu_intercept", "mu_x"), scale, "mu")
  expect_identical(unname(.bt_apply_formula_coefficient_transform(as.matrix(draws), descriptor)[, "mu_intercept"]), c(-7, -7))
})

test_that("fresh multiplier primitives are retained without changing RNG", {

  fit <- .formula_state_test_fit()
  seed <- .caller_rng_state()
  first <- transform_prior_samples(fit, n_samples = 2000, seed = 711)
  expect_identical(.caller_rng_state(), seed)
  attr(fit, "original_fit") <- list(posterior = matrix(1e100, 2, 2))
  second <- transform_prior_samples(fit, n_samples = 2000, seed = 711)
  expect_identical(first, second)
  expect_equal(first[, "mu_intercept"], 5 - 4 * first[, "sigma"], tolerance = 1e-14)
  expect_identical(first[, "mu_x"], rep(.2, 2000))
  plain <- list(theta = prior("normal", list(0, 1)), sigma = prior("gamma", list(2, 1)))
  expected <- .generate_prior_sample_matrix(plain, 20, c("theta", "sigma"), seed = 71)
  actual <- .generate_transformed_prior_samples(plain, c("theta", "sigma"), 20, seed = 71)
  expect_identical(actual, expected)
})

test_that("multiple prefixes use one original fitted snapshot", {

  slope <- prior("point", list(2))
  attr(slope, "multiply_by") <- "nu_z"
  mu <- JAGS_formula(~ x, "mu", data.frame(x = c(10, 20, 30)),
    list(intercept = prior("point", list(5)), x = slope), formula_scale = TRUE)
  nu <- JAGS_formula(~ z, "nu", data.frame(z = c(0, 2, 4)),
    list(intercept = prior("normal", list(0, 1)), z = prior("normal", list(0, 1))), formula_scale = TRUE)
  source <- cbind(mu_intercept = c(5, 5), mu_x = c(2, 2), nu_intercept = c(7, 7), nu_z = c(2, 4))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(source)),
    c(mu$prior_list, nu$prior_list), list(mu = mu$formula_design, nu = nu$formula_design),
    list(mu = mu$formula_scale, nu = nu$formula_scale))
  scales <- attr(fit, "formula_scale")
  first <- .apply_unscale_transform(source, scales)
  second <- .apply_unscale_transform(source, scales[c("nu", "mu")])
  expect_identical(first, second)
  expect_identical(unname(first[, "mu_intercept"]), c(-3, -11))
  expect_identical(unname(first[, "nu_z"]), c(1, 2))
  expect_identical(source[, "nu_z"], c(2, 4))
})

.formula_state_runjags_transport <- function(fit){

  out <- structure(list(mcmc = fit, sample = nrow(as.matrix(fit)), monitor = colnames(as.matrix(fit)),
    summary.pars = list(mutate = NULL)), class = c("runjags", "BayesTools_fit", "list"))
  for(field in c("prior_list", "formula_design", "formula_scale", "parameter_map")) attr(out, field) <- attr(fit, field, exact = TRUE)
  .bt_attach_fit_contract(.bt_attach_draw_geometry(out))
}

test_that("old contexts and incomplete fitted owners require refitting", {

  fit <- .formula_state_test_fit()
  context <- .prior_density_context(attr(fit, "prior_list"), c("mu_intercept", "mu_x"), attr(fit, "formula_scale"))
  old <- context
  old$schema_version <- NULL
  expect_error(.prior_density_from_context(old, c(mu_x = 1)), class = "BayesTools_refit_required")
  old <- context
  old$transforms$mu$descriptor$schema_version <- 2L
  expect_error(.prior_density_from_context(old, c(mu_x = 1)), class = "BayesTools_refit_required")
  design <- attr(fit, "formula_design")
  attr(design$mu$formula_scale, "unscale_design")$complete <- FALSE
  attr(fit, "formula_design") <- design
  expect_error(JAGS_formula_coefficient_transform(fit, "mu"), class = "BayesTools_refit_required")
})

test_that("declared point siblings reduce mixed recipes without pseudo priors", {

  fit <- .formula_state_test_fit()
  law <- JAGS_formula_prior_density(fit, "mu", weights = c(mu_intercept = 1, mu_x = 20))
  expect_equal(prior_density_ordinate(law, 9)$log_density, stats::dnorm(9, 9, 4, log = TRUE), tolerance = 1e-12)
  context <- .prior_density_context(attr(fit, "prior_list"), c("mu_intercept", "mu_x"), attr(fit, "formula_scale"))
  shared <- .prior_density_from_context(context, c(mu_intercept = 1, mu_x = 20))
  expect_equal(prior_density_ordinate(shared, 9)$log_density, stats::dnorm(9, 9, 4, log = TRUE), tolerance = 1e-12)
})
