skip_if_not_test_profile("unit")

test_that("unrepresentable numeric and point folds preserve numerical scale provenance", {

  for(scale in c(1e-200, 1e200)){
    for(named in c(FALSE, TRUE)){
      prior <- prior("normal", list(0, 1 / scale))
      attr(prior, "multiply_by") <- if(named) "s" else scale
      priors <- list(x = prior)
      if(named) priors$s <- prior("point", list(scale))
      weights <- c(x = scale)
      route <- .prior_density_route_linear(priors, weights, NULL, 4096L)
      expect_identical(route$type, "unknown")
      expect_null(route$recipe)
      expect_true(is.character(route$reason) && any(grepl("Numerical scale", route$reason, fixed = TRUE)))
      expect_identical(route$provenance$kind, "numerical_scale_unavailable")
      context <- .prior_density_build_context(priors, names(priors))
      context_route <- .prior_density_route_context(context, weights, NULL, NULL, NULL)
      expect_identical(context_route$type, "unknown")
      expect_null(context_route$recipe)
      expect_error(.prior_density_from_context(context, weights), class = "BayesTools_numerical_condition")
      if(!named) expect_error(.prior_linear_split_multiply_groups(priors, weights),
        class = "BayesTools_numerical_unavailable")
    }
  }
})

test_that("formula numerical fold failures keep values and refuse only unavailable measures", {

  for(scale in c(1e-200, 1e200)) for(named in c(FALSE, TRUE)){
    x_prior <- prior("normal", list(0, 1 / scale))
    attr(x_prior, "multiply_by") <- if(named) "s" else scale
    compiled <- JAGS_formula(~ 0 + x, "mu", data.frame(x = c(-1, 0, 1)), list(x = x_prior))
    draws <- cbind(mu_x = c(-1, 0, 1) / scale, mu_intercept = 0)
    priors <- compiled$prior_list
    if(named){
      priors$s <- prior("point", list(scale))
      draws <- cbind(draws, s = scale)
    }
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), priors,
      list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
    mixed <- as_mixed_posteriors(fit, names(compiled$prior_list), transform_scaled = FALSE)
    leaf <- marginal_posterior(mixed, "mu_intercept", formula = ~ 0 + x,
      at = list(x = scale), prior_samples = TRUE)$intercept
    expect_identical(as.numeric(leaf), c(-1, 0, 1) * scale)
    unavailable <- posterior_metadata(leaf, "measure_unavailable")
    expect_identical(sort(unavailable$measure), sort(c("prior_density", "atoms", "support")))
    expect_true(all(grepl("Numerical scale", unavailable$reason, fixed = TRUE)))
    expect_null(posterior_metadata(leaf, "prior_density"))
    expect_null(posterior_metadata(leaf, "atoms"))
    expect_null(posterior_metadata(leaf, "support"))
    expect_error(.bt_formula_measure_check(leaf, "prior_density", "intercept"),
      class = "BayesTools_formula_measure_unavailable")
    sibling <- marginal_posterior(mixed, "mu_intercept", formula = ~ 0 + x,
      at = list(x = 0), prior_samples = TRUE)$intercept
    expect_identical(as.numeric(sibling), rep(0, 3))
    expect_null(posterior_metadata(sibling, "measure_unavailable"))
    expect_identical(prior_density_ordinate(posterior_metadata(sibling, "prior_density"), 0)$point_mass, 1)
  }
})

test_that("ordinary and genuine zero multiply folds retain declared measures", {

  for(named in c(FALSE, TRUE)){
    for(scale in c(0, 2)){
      x <- prior("normal", list(0, 1))
      attr(x, "multiply_by") <- if(named) "s" else scale
      priors <- list(x = x)
      if(named) priors$s <- prior("point", list(scale))
      route <- .prior_density_route_linear(priors, c(x = 3), NULL, 4096L)
      ordinate <- .prior_density_route_ordinate(route, 0)
      expect_identical(ordinate$point_mass, if(scale == 0) 1 else 0)
      if(scale != 0) expect_equal(ordinate$log_density, dnorm(0, 0, 6, log = TRUE), tolerance = 1e-14)
    }
  }
})
test_that("independent additive formula grids retain their numerical laws", {

  data <- data.frame(x1 = c(-1, 0, 1, -1, 0, 1), x2 = c(1, -1, 0, 0, 1, -1))
  compiled <- JAGS_formula(~ x1 + x2, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x1 = prior("t", list(0, 1, 5)),
      x2 = prior("t", list(0, 1, 7))))
  draws <- cbind(mu_intercept = seq(-1, 1, length.out = 40),
    mu_x1 = seq(-2, 2, length.out = 40), mu_x2 = cos(1:40) / 3)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), compiled$prior_list,
    list(mu = compiled$formula_design), list(mu = compiled$formula_scale))
  mixed <- as_mixed_posteriors(fit, names(compiled$prior_list), transform_scaled = FALSE)
  levels <- marginal_posterior(mixed, "mu_x1", formula = ~ x1 + x2,
    at = list(x2 = 1), prior_samples = TRUE)
  for(i in seq_along(levels)){
    x1 <- c(-1, 0, 1)[i]
    leaf <- levels[[i]]
    expect_equal(as.numeric(leaf), draws[, "mu_intercept"] + x1 * draws[, "mu_x1"] + draws[, "mu_x2"])
    law <- posterior_metadata(leaf, "prior_density")
    expect_s3_class(law, "prior_linear_density")
    expect_true(posterior_atoms_free(leaf))
    support <- posterior_metadata(leaf, "support")
    expect_true(support$exact)
    expect_identical(unname(support$bounds), c(-Inf, Inf))
    if(x1 != 0){
      ordinate <- prior_density_ordinate(law, 0)
      expect_identical(ordinate$behavior, "unknown")
      expect_false(ordinate$exact)
      expect_error(hypothesis_BF(leaf, hypothesis = "theta = 0", parameter = "theta"),
        class = "BayesTools_inexact_ordinate")
      route <- .prior_density_route_linear(compiled$prior_list,
        c(mu_intercept = 1, mu_x1 = x1, mu_x2 = 1), NULL, law$n_grid)
      reference <- .prior_density_route_recipe_grid(route$recipe)
      expect_equal(law$density, reference$density)
      expect_equal(law$point, reference$point)
    }
  }
  recipe <- .prior_density_route_recipe(list(x = prior("uniform", list(1, 2)),
    y = prior("t", list(0, 1, 5)), z = prior("point", list(3))), c(x = -2, y = 1, z = 1), NULL, 4096)
  expect_identical(.prior_density_route_additive_measure(recipe)$type, "atom_free")
  for(transform in c("log", "unknown")){
    unsupported <- recipe
    unsupported$source_transforms <- c(x = transform)
    expect_null(.prior_density_route_additive_measure(unsupported))
  }
  unsupported <- recipe
  attr(unsupported$prior_list$x, "multiply_by") <- 2
  expect_null(.prior_density_route_additive_measure(unsupported))
  unsupported <- recipe
  unsupported$prior_list$z <- prior("point", list(expression(x)))
  expect_null(.prior_density_route_additive_measure(unsupported))
  unsupported <- .prior_density_route_recipe(list(x = prior("mnormal", list(0, 1, 2)),
    y = recipe$prior_list$y), c("x[1]" = 1, "x[2]" = 1, y = 1), NULL, 4096)
  expect_null(.prior_density_route_additive_measure(unsupported))
})
