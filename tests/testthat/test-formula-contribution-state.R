skip_if_not_test_profile("unit")

test_that("complete declarations certify absent terms with a valid empty source state", {
  compiled <- JAGS_formula(~ 1, "mu", data.frame(x = c(1, 2, 3)),
    list(intercept = prior("normal", list(0, 1))))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(cbind(mu_intercept = c(1, 2, 3)))),
    compiled$prior_list, list(mu = compiled$formula_design))
  state <- .bt_formula_state_new(fit, as.matrix(fit), "mu_g")
  expect_identical(dim(state$values), c(3L, 0L))
  expect_identical(names(state$models[[1L]]$declaration_priors), "mu_intercept")
  expect_null(.bt_formula_state_validate(state))
  malformed <- state
  malformed$models[[1L]]$required <- "missing"
  expect_identical(.bt_formula_state_validate(malformed), "required own-model formula states must be finite and present")
  broken <- fit
  design <- attr(broken, "formula_design")
  design$mu$prior_list <- list()
  attr(broken, "formula_design") <- design
  expect_error(.bt_formula_state_new(broken, as.matrix(broken), "mu_g"), class = "BayesTools_formula_transform_unavailable")
})

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

test_that("fitted and original formula values and measures use fitted rows", {

  for(geometry in list(c(20, 10), c(10, 3), c(3, 7), c(0, 10))){
    fit <- .formula_state_test_fit(mean = geometry[[1L]], sd = geometry[[2L]])
    raw <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"))
    original <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
    for(original_scale in c(FALSE, TRUE)){
      view <- if(original_scale) original else raw
      row <- marginal_posterior(view, "mu_intercept", formula = ~ 1 + x,
        at = list(x = if(original_scale) geometry[[1L]] else 0), prior_samples = TRUE)$intercept
      expect_identical(as.numeric(row), c(5, 5))
      expect_identical(.posterior_atoms_get(row)$mass, 1)
      expect_identical(unname(.posterior_atoms_get(row)$locations), matrix(5))
      expect_identical(posterior_metadata(row, "linear_weight_space"), "formula_contribution")
      expect_identical(posterior_metadata(row, "linear_weights")[1L, ], c(mu_intercept = 1, mu_x = 0))
    }
  }
  fit <- .formula_state_test_fit(mean = 0)
  children <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  children <- children[c("mu_intercept", "mu_x")]
  expect_false("sigma" %in% names(children))
  rm(fit)
  row <- marginal_posterior(children, "mu_intercept", formula = ~ 1 + x,
    at = list(x = 1))$intercept
  expect_equal(as.numeric(row), c(5.4, 5.8), tolerance = 1e-14)
})

test_that("public at values preserve fitted SD and original data units", {

  fit <- .formula_state_test_fit(multiplier = NULL, extra_priors = list(),
    values = cbind(mu_intercept = c(5, 5), mu_x = c(2, 2)))
  fitted <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"))
  original <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  for(at in c(0, 1, 20, 30)){
    row <- marginal_posterior(fitted, "mu_intercept", ~ x, list(x = at))$intercept
    expect_identical(as.numeric(row), rep(5 + 2 * at, 2L))
  }
  for(at in c(0, 1)){
    fitted_row <- marginal_posterior(fitted, "mu_intercept", ~ x, list(x = at), prior_samples = TRUE)$intercept
    original_row <- marginal_posterior(original, "mu_intercept", ~ x,
      list(x = 20 + 10 * at), prior_samples = TRUE)$intercept
    expect_identical(as.numeric(fitted_row), as.numeric(original_row))
    expect_identical(posterior_metadata(fitted_row, "linear_weights"), posterior_metadata(original_row, "linear_weights"))
    expect_identical(.posterior_support_get(fitted_row), .posterior_support_get(original_row))
    expect_identical(.posterior_atoms_get(fitted_row), .posterior_atoms_get(original_row))
  }
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

test_that("compiled signs and products determine contribution support", {

  fit <- .formula_state_test_fit(multiplier = -3,
    slope = prior("uniform", list(0, 1)), values = cbind(mu_intercept = c(5, 5), mu_x = c(.2, .7)),
    extra_priors = list())
  for(original in c(FALSE, TRUE)){
    mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = original)
    row <- marginal_posterior(mixed, "mu_intercept", formula = ~ x,
      at = list(x = if(original) 30 else 1), prior_samples = TRUE)$intercept
    expect_identical(.posterior_support_get(row)$bounds, c(2, 5))
    expect_equal(as.numeric(row), c(4.4, 2.9), tolerance = 1e-14)
  }
  fit <- .formula_state_test_fit(slope = prior("beta", list(2, 3)),
    extra_priors = list(sigma = prior("normal", list(0, 1), truncation = list(lower = 0))),
    values = cbind(mu_intercept = c(5, 5), mu_x = c(.2, .7), sigma = c(2, 4)))
  mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  expect_identical(.posterior_support_get(mixed$mu_intercept)$bounds, c(-Inf, 5))
  expect_identical(.posterior_support_get(mixed$mu_x)$bounds, c(0, .1))
})

test_that("posterior atoms use actual joint gate frequencies", {

  sigma <- prior_mixture(list(prior("point", list(0)),
    prior("normal", list(0, 1), truncation = list(lower = 0), prior_weights = 3)))
  fit <- .formula_state_test_fit(extra_priors = list(sigma = sigma),
    values = cbind(mu_intercept = rep(5, 4), mu_x = rep(2, 4),
      sigma = c(0, 0, 0, 2), sigma_indicator = c(1, 1, 1, 2)))
  mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"))
  row <- marginal_posterior(mixed, "mu_intercept", formula = ~ x,
    at = list(x = 1), prior_samples = TRUE)$intercept
  expect_identical(as.numeric(row), c(5, 5, 5, 9))
  expect_identical(.posterior_atoms_get(row)$mass, .75)
  expect_identical(unname(.posterior_atoms_get(row)$locations), matrix(5))
  law <- posterior_metadata(row, "prior_density")
  expect_equal(law$points$p[law$points$x == 5], .25, tolerance = 1e-15)
  components <- .posterior_components_get(row)
  expect_identical(components$index, c(1L, 1L, 1L, 2L))
  expect_identical(components$supports[[1L]]$bounds, c(5, 5))
  expect_identical(components$supports[[2L]]$bounds, c(5, Inf))
})

test_that("catalog ordinary mixture components preserve estimator parity", {

  for(prior in list(prior_mixture(list(prior("uniform", list(0, 1)),
    prior("uniform", list(10, 11))), is_null = c(FALSE, FALSE)),
    prior_spike_and_slab(prior("uniform", list(0, 1)), prior("point", list(.5))))){
    spike <- is.prior.spike_and_slab(prior)
    values <- if(spike) c(rep(0, 500), seq(.0005, .9995, length.out = 500)) else
      c(seq(.0005, .9995, length.out = 500), seq(10.0005, 10.9995, length.out = 500))
    indicator <- if(spike) rep(0:1, each = 500) else rep(1:2, each = 500)
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(cbind(theta = values,
      theta_indicator = indicator))), list(theta = prior))
    selection <- parameter_catalog_resolve(parameter_catalog(fit), "theta")
    semantic <- parameter_mixed_posterior(fit, selection)
    ordinary <- marginal_posterior(as_mixed_posteriors(fit, "theta"), "theta",
      prior_samples = TRUE, use_formula = FALSE)
    expect_identical(as.numeric(semantic), as.numeric(ordinary))
    expect_identical(.posterior_components_get(semantic), .posterior_components_get(ordinary))
    expect_equal(.posterior_atoms_get(semantic)$mass, .posterior_atoms_get(ordinary)$mass, tolerance = 0)
    expect_identical(.bt_meta_condition(semantic, "averaged"), .bt_meta_condition(ordinary, "averaged"))
    expect_equal(as.numeric(Savage_Dickey_BF(semantic, 1, silent = TRUE)),
      as.numeric(Savage_Dickey_BF(ordinary, 1, silent = TRUE)), tolerance = 1e-14)
    expect_identical(colnames(as.matrix(parameter_draws(fit, selection))), "theta")
    expect_identical(.bt_meta_get(semantic, "draw_index"), seq_along(values))
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

test_that("authoritative scalar coefficient measures follow requested transformations", {

  fit <- .formula_state_test_fit()
  mixed <- as_mixed_posteriors(fit, "mu_x", transform_scaled = TRUE)
  result <- marginal_posterior(mixed, "mu_x", use_formula = FALSE, prior_samples = TRUE,
    transformation = "lin", transformation_arguments = list(a = 1, b = -2))
  expect_equal(as.numeric(result), c(.6, .6), tolerance = 1e-14)
  expect_equal(as.numeric(.posterior_atoms_get(result)$locations), .6, tolerance = 1e-14)
  expect_identical(.posterior_atoms_get(result)$mass, 1)
  expect_equal(posterior_metadata(result, "prior_density")$points$x, .6, tolerance = 1e-14)
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

test_that("positive posterior models without selected rows keep their atoms", {

  fit1 <- .formula_state_runjags_transport(.formula_state_test_fit())
  fit2 <- .formula_state_runjags_transport(.formula_state_test_fit(extra_priors = list(sigma = prior("point", list(0))),
    values = cbind(mu_intercept = c(5, 5), mu_x = c(2, 2), sigma = c(0, 0))))
  model_list <- list(list(fit = fit1, marglik = bridgesampling_object(0), prior_weights = .95),
    list(fit = fit2, marglik = bridgesampling_object(0), prior_weights = .05))
  mixed <- mix_posteriors(model_list, c("mu_intercept", "mu_x"),
    list(mu_intercept = c(FALSE, FALSE), mu_x = c(FALSE, FALSE)), seed = 73, n_samples = 2)
  state <- .bt_meta_get(mixed, "formula_state")
  expect_equal(state$posterior_model_probabilities, c(.95, .05), tolerance = 1e-15)
  expect_identical(state$model, c(1L, 1L))
  row <- marginal_posterior(mixed, "mu_intercept", ~ x, list(x = 1), prior_samples = TRUE)$intercept
  expect_identical(unname(.posterior_atoms_get(row)$locations[1L, 1L]), 5)
  expect_equal(.posterior_atoms_get(row)$mass, .05, tolerance = 1e-15)
  expect_length(.posterior_components_get(row)$supports, 2L)
})

test_that("heterogeneous model contributions retain foreign NA holes and leaf laws", {

  fit1 <- .formula_state_runjags_transport(.formula_state_test_fit())
  fit2 <- .formula_state_runjags_transport(.formula_state_test_fit(multiplier = -3, extra_priors = list(),
    values = cbind(mu_intercept = c(5, 5), mu_x = c(2, 2))))
  mixed <- mix_posteriors(list(list(fit = fit1, marglik = bridgesampling_object(0), prior_weights = .6),
    list(fit = fit2, marglik = bridgesampling_object(0), prior_weights = .4)),
    c("mu_intercept", "mu_x"), list(c(FALSE, FALSE), c(FALSE, FALSE)), seed = 74, n_samples = 10)
  state <- .bt_meta_get(mixed, "formula_state")
  expect_true(all(is.na(state$values[state$model == 2L, "sigma"])))
  row <- marginal_posterior(mixed, "mu_intercept", ~ x, list(x = 1), prior_samples = TRUE)$intercept
  expect_true(all(as.numeric(row)[state$model == 2L] == -1))
  expect_equal(.posterior_atoms_get(row)$mass, .4, tolerance = 1e-15)
  law <- posterior_metadata(row, "prior_density")
  region <- list(intervals = cbind(-Inf, 7), indicator = function(x) x <= 7)
  probability <- .prior_linear_density_region_probability(law, region)
  expect_equal(as.numeric(probability), .6 * stats::pnorm(1) + .4, tolerance = 1e-12)
})

test_that("conditioned gate evidence uses original eligible rows", {

  sigma <- prior_mixture(list(prior("point", list(0)),
    prior("point", list(2)), prior("normal", list(0, 1), truncation = list(lower = 0))),
    is_null = c(TRUE, FALSE, FALSE))
  fit <- .formula_state_test_fit(extra_priors = list(sigma = sigma),
    values = cbind(mu_intercept = rep(5, 6), mu_x = rep(2, 6),
      sigma = c(0, 0, 2, 2, 2, 4), sigma_indicator = c(1, 1, 2, 2, 2, 3)))
  mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "sigma"), conditional = "sigma")
  selected <- mixed[c("mu_intercept", "mu_x")]
  state <- .bt_meta_get(selected, "formula_state")
  expect_identical(state$draw_index, 3:6)
  expect_identical(state$models[[1L]]$eligible_n, 4L)
  row <- marginal_posterior(selected, "mu_intercept", ~ x, list(x = 1), prior_samples = TRUE)$intercept
  expect_identical(unname(.posterior_atoms_get(row)$locations[1L, 1L]), 9)
  expect_identical(.posterior_atoms_get(row)$mass, .75)
  expect_identical(posterior_metadata(row, "prior_density")$points$p, .5)
})

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

test_that("ordered projections distinguish raw and compiled contribution spaces", {

  data <- expand.grid(f = ordered(c("early", "middle", "late"),
    levels = c("early", "middle", "late")), x = c(10, 20, 30))
  ordered <- prior_ordered(prior("point", list(1)))
  attr(ordered, "multiply_by") <- 2
  slope <- prior("point", list(2))
  attr(slope, "multiply_by") <- "sigma"
  compiled <- JAGS_formula(~ f + x, "mu", data,
    list(intercept = prior("point", list(5)), f = ordered, x = slope), formula_scale = TRUE)
  prior <- compiled$prior_list$mu_f
  set.seed(177)
  draws <- .prior_ordered_draws(prior, 2L)
  coefficients <- draws$coefficients
  colnames(coefficients) <- .JAGS_prior_factor_names("mu_f", prior)
  sources <- .prior_ordered_total_samples(draws, "mu_f")
  spec <- .bt_ordered_spec("mu_f", prior)
  for(record in .prior_ordered_dirichlet_records(prior)){
    gamma <- draws$allocation_samples[[record$key]] * seq_len(2L)
    colnames(gamma) <- spec$allocations[[record$key]]$gamma_coordinates
    sources <- cbind(sources, gamma)
  }
  values <- cbind(mu_intercept = c(5, 5), mu_x = c(2, 2), sigma = c(3, 3), coefficients, sources)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(values)),
    c(compiled$prior_list, list(sigma = prior("point", list(3)))), list(mu = compiled$formula_design),
    list(mu = compiled$formula_scale))
  raw <- as_mixed_posteriors(fit, c("mu_intercept", "mu_f", "mu_x", "sigma"))
  weights <- c(mu_intercept = 0, `mu_f[1]` = 1, `mu_f[2]` = 1, mu_x = 0)
  coefficient <- .bt_ordered_formula_projections(raw, matrix(weights, 1L,
    dimnames = list(NULL, names(weights))), weight_space = "coefficient")
  expect_equal(coefficient[[1L]]$values, c(1, 1), tolerance = 1e-14)
  weights <- c(mu_intercept = 1, `mu_f[1]` = 1, `mu_f[2]` = 1, mu_x = 1)
  contribution <- .bt_ordered_formula_projections(raw, matrix(weights, 1L,
    dimnames = list(NULL, names(weights))), weight_space = "formula_contribution")
  expect_equal(contribution[[1L]]$values, c(13, 13), tolerance = 1e-14)
  at_x <- compiled$formula_scale$mu_x$mean + compiled$formula_scale$mu_x$sd
  for(original in c(FALSE, TRUE)){
    mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_f", "mu_x"), transform_scaled = original)
    row <- marginal_posterior(mixed, "mu_intercept", ~ f + x,
      list(f = "late", x = if(original) at_x else 1), prior_samples = TRUE)$intercept
    expect_equal(as.numeric(row), c(13, 13), tolerance = 1e-14)
  }
})
test_that("ordered contribution replay declines lost finite multiplier folds", {

  data <- expand.grid(f = ordered(c("early", "middle", "late"), levels = c("early", "middle", "late")), x = c(10, 20, 30))
  for(multiplier in list(1e-200, "sigma", 0)){
    ordered <- prior_ordered(prior("point", list(1e200)))
    attr(ordered, "multiply_by") <- multiplier
    compiled <- JAGS_formula(~ f + x, "mu", data,
      list(intercept = prior("point", list(0)), f = ordered, x = prior("point", list(0))), formula_scale = TRUE)
    prior <- compiled$prior_list$mu_f
    set.seed(177)
    draws <- .prior_ordered_draws(prior, 2L)
    coefficients <- draws$coefficients
    colnames(coefficients) <- .JAGS_prior_factor_names("mu_f", prior)
    sources <- .prior_ordered_total_samples(draws, "mu_f")
    spec <- .bt_ordered_spec("mu_f", prior)
    for(record in .prior_ordered_dirichlet_records(prior)){
      gamma <- draws$allocation_samples[[record$key]] * seq_len(2L)
      colnames(gamma) <- spec$allocations[[record$key]]$gamma_coordinates
      sources <- cbind(sources, gamma)
    }
    values <- cbind(mu_intercept = c(0, 0), mu_x = c(0, 0), sigma = rep(1e-200, 2L), coefficients, sources)
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(values)),
      c(compiled$prior_list, list(sigma = prior("point", list(1e-200)))), list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
    mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_f", "mu_x", "sigma"))
    weights <- c(mu_intercept = 0, `mu_f[1]` = 1e-200, `mu_f[2]` = 1e-200, mu_x = 0)
    actual <- sum(coefficients[1L, ] * if(identical(multiplier, 0)) 0 else 1e-200) * 1e-200
    if(identical(multiplier, 0)) expect_identical(actual, 0) else expect_equal(actual / 1e-200, 1, tolerance = 1e-14)
    replay <- .bt_ordered_formula_projections(mixed, matrix(weights, 1L, dimnames = list(NULL, names(weights))), weight_space = "formula_contribution")
    if(identical(multiplier, 0)){
      expect_identical(replay[[1L]]$values, c(0, 0))
    }else{
      expect_null(replay)
    }
  }
})
