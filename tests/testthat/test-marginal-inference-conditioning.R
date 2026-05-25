skip_if_not_test_profile("unit")

.mock_marginal_fit <- function(posterior, prior_list) {
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  fit
}

test_that("as_mixed_posteriors applies AND and OR conditioning exactly", {

  prior_list <- list(
    mu_a = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(0.5))
    ),
    mu_b = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(0.5))
    )
  )
  posterior <- cbind(
    mu_a = c(0, 10, 0, 30, 31),
    mu_b = c(0, 0, 20, 40, 41),
    mu_a_indicator = c(0, 1, 0, 1, 1),
    mu_b_indicator = c(0, 0, 1, 1, 1)
  )
  fit <- .mock_marginal_fit(posterior, prior_list)

  and_samples <- as_mixed_posteriors(
    fit,
    parameters = c("mu_a", "mu_b"),
    conditional = c("mu_a", "mu_b"),
    conditional_rule = "AND"
  )

  expect_equal(as.numeric(and_samples$mu_a), c(30, 31))
  expect_equal(as.numeric(and_samples$mu_b), c(40, 41))
  expect_equal(attr(and_samples$mu_a, "models_ind"), c(1, 1))
  expect_equal(attr(and_samples$mu_b, "models_ind"), c(1, 1))

  or_samples <- as_mixed_posteriors(
    fit,
    parameters = c("mu_a", "mu_b"),
    conditional = c("mu_a", "mu_b"),
    conditional_rule = "OR"
  )

  expect_equal(as.numeric(or_samples$mu_a), c(10, 0, 30, 31))
  expect_equal(as.numeric(or_samples$mu_b), c(0, 20, 40, 41))
  expect_equal(attr(or_samples$mu_a, "models_ind"), c(1, 0, 1, 1))
  expect_equal(attr(or_samples$mu_b, "models_ind"), c(0, 1, 1, 1))
})

test_that("as_mixed_posteriors propagates named upstream posterior densities", {

  prior_list <- list(theta = prior("normal", list(0, 1)))
  posterior <- cbind(theta = seq(-1, 1, length.out = 21))
  fit <- .mock_marginal_fit(posterior, prior_list)

  stored_x <- seq(-2, 2, length.out = 41)
  stored_y <- stats::dnorm(stored_x, mean = .25, sd = .8)
  attr(fit, "posterior_density") <- list(
    theta = list(
      parameter = "theta",
      x         = stored_x,
      y         = stored_y,
      method    = "iwmde"
    )
  )

  mixed <- as_mixed_posteriors(fit, parameters = "theta")
  marginal <- marginal_posterior(
    samples       = mixed,
    parameter     = "theta",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 128
  )
  plot_data <- BayesTools:::.plot_data_samples.simple(
    samples                  = mixed,
    parameter                = "theta",
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  expect_equal(attr(mixed$theta, "posterior_density")$method, "iwmde")
  expect_equal(attr(marginal, "posterior_density")$x, stored_x)
  expect_equal(plot_data$density$x, stored_x)
  expect_equal(plot_data$density$y, stored_y)
})

test_that("as_mixed_posteriors does not reuse stale conditional densities", {

  prior_list <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(0.5))
    ),
    phi = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(0.5))
    )
  )
  posterior <- cbind(
    theta           = c(0, 1, 0, 2, 3),
    phi             = c(0, 0, 4, 5, 6),
    theta_indicator = c(0, 1, 0, 1, 1),
    phi_indicator   = c(0, 0, 1, 1, 1)
  )
  stored_density <- list(
    parameter = "theta",
    x         = seq(-3, 3, length.out = 61),
    y         = rep(1, 61),
    method    = "iwmde"
  )
  fit <- .mock_marginal_fit(posterior, prior_list)
  attr(fit, "posterior_density") <- list(theta = stored_density)

  stale <- as_mixed_posteriors(
    fit,
    parameters       = c("theta", "phi"),
    conditional      = c("theta", "phi"),
    conditional_rule = "OR"
  )
  expect_null(attr(stale$theta, "posterior_density"))

  stored_density$conditional <- c("phi", "theta")
  stored_density$conditional_rule <- "OR"
  attr(fit, "posterior_density") <- list(theta = stored_density)
  matched <- as_mixed_posteriors(
    fit,
    parameters       = c("theta", "phi"),
    conditional      = c("theta", "phi"),
    conditional_rule = "OR"
  )
  mismatched <- as_mixed_posteriors(
    fit,
    parameters       = c("theta", "phi"),
    conditional      = c("theta", "phi"),
    conditional_rule = "AND"
  )

  expect_equal(attr(matched$theta, "posterior_density")$method, "iwmde")
  expect_null(attr(mismatched$theta, "posterior_density"))
})

test_that("formula marginals do not inherit raw coefficient densities", {

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x         = prior("normal", list(0, 1))
  )
  attr(prior_list$mu_intercept, "parameter") <- "mu"
  attr(prior_list$mu_x, "parameter") <- "mu"

  fit <- .mock_marginal_fit(
    cbind(
      mu_intercept = seq(-.5, .5, length.out = 11),
      mu_x         = seq(1, 2, length.out = 11)
    ),
    prior_list
  )
  attr(fit, "posterior_density") <- list(
    mu_x = list(
      parameter = "mu_x",
      x         = seq(-2, 2, length.out = 41),
      y         = rep(1, 41),
      method    = "iwmde"
    )
  )

  samples <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x"))
  marginal <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_x",
    formula       = ~ x,
    prior_samples = TRUE,
    n_samples     = 128
  )

  expect_true(all(vapply(marginal, function(level) {
    is.null(attr(level, "posterior_density"))
  }, logical(1))))
})

test_that("factor marginals attach only level-matched densities", {

  formula_result <- JAGS_formula(
    formula    = ~ fac,
    parameter  = "mu",
    data       = data.frame(fac = factor(c("A", "B", "C"), levels = c("A", "B", "C"))),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      fac       = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  fit <- .mock_marginal_fit(
    cbind(
      "mu_fac[1]" = seq(-1, 0, length.out = 11),
      "mu_fac[2]" = seq(1, 2, length.out = 11)
    ),
    formula_result[["prior_list"]]
  )
  stored_x <- seq(-3, 3, length.out = 31)
  attr(fit, "posterior_density") <- list(
    "mu_fac[1]" = list(
      parameter = "mu_fac[1]",
      x         = stored_x,
      y         = rep(100, length(stored_x)),
      method    = "raw-coefficient"
    ),
    B = list(
      parameter = "B",
      x         = stored_x,
      y         = stats::dnorm(stored_x, mean = -0.5, sd = .8),
      method    = "iwmde"
    ),
    C = list(
      parameter = "C",
      x         = stored_x,
      y         = stats::dnorm(stored_x, mean = 1.5, sd = .8),
      method    = "iwmde"
    )
  )

  samples <- as_mixed_posteriors(fit, parameters = "mu_fac")
  marginal <- marginal_posterior(
    samples     = samples,
    parameter   = "mu_fac",
    use_formula = FALSE
  )

  expect_null(attr(marginal$A, "posterior_density"))
  expect_equal(attr(marginal$B, "posterior_density")$method, "iwmde")
  expect_equal(attr(marginal$C, "posterior_density")$method, "iwmde")
  expect_false(identical(attr(marginal$B, "posterior_density")$method, "raw-coefficient"))
})

test_that("as_marginal_inference rejects precomputed marginal-inference BFs", {

  prior_list <- list(theta = prior("normal", list(0, 1)))
  fit <- .mock_marginal_fit(
    cbind(theta = seq(-3, 3, length.out = 301)),
    prior_list
  )
  stored_x <- seq(-4, 4, length.out = 401)
  stored_y <- stats::dnorm(stored_x, mean = .75, sd = .9)
  attr(fit, "posterior_density") <- list(
    theta = list(
      parameter = "theta",
      x         = stored_x,
      y         = stored_y,
      method    = "iwmde"
    )
  )

  expect_error(
    as_marginal_inference(
      model                = fit,
      marginal_parameters = "theta",
      parameters          = "theta",
      conditional_list    = list(theta = NULL),
      conditional_rule    = "AND",
      formula             = NULL,
      n_samples           = 256,
      silent              = TRUE,
      density_method      = "precomputed"
    ),
    "not supported"
  )
  inference <- as_marginal_inference(
    model                = fit,
    marginal_parameters = "theta",
    parameters          = "theta",
    conditional_list    = list(theta = NULL),
    conditional_rule    = "AND",
    formula             = NULL,
    n_samples           = 256,
    silent              = TRUE,
    density_method      = "KDE"
  )

  expect_equal(attr(inference, "density_method"), "KDE")
})

test_that("as_marginal_inference does not consume raw stored posterior ordinates", {

  prior_list <- list(theta = prior("normal", list(0, 1)))
  fit <- .mock_marginal_fit(
    cbind(theta = seq(-3, 3, length.out = 301)),
    prior_list
  )
  attr(fit, "posterior_ordinate") <- list(
    theta = list(
      parameter   = "theta",
      value       = .25,
      ordinate    = .5,
      method      = "qCMDE",
      diagnostics = list(relative_mcse = .2)
    )
  )

  expect_error(
    as_marginal_inference(
      model                = fit,
      marginal_parameters = "theta",
      parameters          = "theta",
      conditional_list    = list(theta = NULL),
      conditional_rule    = "AND",
      formula             = NULL,
      null_hypothesis     = .25,
      n_samples           = 256,
      silent              = TRUE,
      density_method      = "precomputed"
    ),
    "not supported"
  )
})

test_that("marginal_posterior evaluates formula terms exactly on default grid", {

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x = prior("normal", list(0, 1))
  )
  attr(prior_list$mu_intercept, "parameter") <- "mu"
  attr(prior_list$mu_x, "parameter") <- "mu"

  fit <- .mock_marginal_fit(
    cbind(mu_intercept = c(10, 20), mu_x = c(1, 2)),
    prior_list
  )
  samples <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x"))
  marginal <- marginal_posterior(
    samples = samples,
    parameter = "mu_x",
    formula = ~ x,
    prior_samples = FALSE
  )

  expect_equal(names(marginal), c("-1SD", "0SD", "1SD"))
  expect_equal(as.numeric(marginal[["-1SD"]]), c(9, 18))
  expect_equal(as.numeric(marginal[["0SD"]]), c(10, 20))
  expect_equal(as.numeric(marginal[["1SD"]]), c(11, 22))
})

test_that("conditional spike-and-slab prior densities use the slab", {

  prior_list <- list(
    mu_x = prior_spike_and_slab(
      prior("normal", list(mean = 1, sd = 0.2)),
      prior_inclusion = prior("point", list(location = 0.5))
    )
  )

  posterior <- cbind(
    mu_x           = c(rep(0, 50), seq(0.5, 1.5, length.out = 50)),
    mu_x_indicator = c(rep(0, 50), rep(1, 50))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list

  samples <- as_mixed_posteriors(
    model       = fit,
    parameters  = "mu_x",
    conditional = "mu_x"
  )
  marginal <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_x",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 128
  )

  expect_equal(.prior_linear_density_point_mass(attr(marginal, "prior_density"), 0), 0)
})


test_that("marginal inference conditions formula levels by active weights", {

  prior_list <- list(
    mu_intercept = prior_spike_and_slab(
      prior("normal", list(mean = 1, sd = 0.2)),
      prior_inclusion = prior("point", list(location = 0.5))
    ),
    mu_x = prior_spike_and_slab(
      prior("normal", list(mean = 0.5, sd = 0.2)),
      prior_inclusion = prior("point", list(location = 0.5))
    )
  )
  attr(prior_list[["mu_intercept"]], "parameter") <- "mu"
  attr(prior_list[["mu_x"]], "parameter") <- "mu"

  indicators <- expand.grid(
    mu_intercept_indicator = c(0, 1),
    mu_x_indicator         = c(0, 1)
  )
  indicators <- indicators[rep(seq_len(nrow(indicators)), each = 50), ]
  posterior <- cbind(
    mu_intercept = ifelse(
      indicators[["mu_intercept_indicator"]] == 1,
      seq(0.75, 1.25, length.out = nrow(indicators)),
      0
    ),
    mu_x = ifelse(
      indicators[["mu_x_indicator"]] == 1,
      seq(0.25, 0.75, length.out = nrow(indicators)),
      0
    ),
    indicators
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list

  inference <- as_marginal_inference(
    model                = fit,
    marginal_parameters = c("mu_intercept", "mu_x"),
    parameters          = c("mu_intercept", "mu_x"),
    conditional_list    = list(
      mu_intercept = c("mu_intercept", "mu_x"),
      mu_x         = c("mu_intercept", "mu_x")
    ),
    conditional_rule    = "OR",
    formula             = ~ x,
    n_samples           = 128,
    silent              = TRUE
  )

  zero_level <- inference[["conditional"]][["mu_x"]][["0SD"]]
  intercept <- inference[["conditional"]][["mu_intercept"]][["intercept"]]

  expect_equal(attr(zero_level, "effective_conditional"), "mu_intercept")
  expect_equal(mean(as.numeric(zero_level) == 0), 0)
  expect_equal(.prior_linear_density_point_mass(attr(zero_level, "prior_density"), 0), 0)
  expect_equal(attr(intercept, "effective_conditional"), "mu_intercept")
  expect_equal(mean(as.numeric(intercept) == 0), 0)
  expect_equal(.prior_linear_density_point_mass(attr(intercept, "prior_density"), 0), 0)

  inference_mu_only <- as_marginal_inference(
    model                = fit,
    marginal_parameters = "mu_x",
    parameters          = c("mu_intercept", "mu_x"),
    conditional_list    = list(mu_x = "mu_x"),
    conditional_rule    = "OR",
    formula             = ~ x,
    n_samples           = 128,
    silent              = TRUE
  )

  zero_level_mu_only <- inference_mu_only[["conditional"]][["mu_x"]][["0SD"]]

  expect_equal(attr(zero_level_mu_only, "effective_conditional"), character())
  expect_true(is.numeric(zero_level_mu_only))
  expect_gt(length(zero_level_mu_only), 0)
})


test_that("marginal inference conditions treatment factor levels by active weights", {

  formula_result <- JAGS_formula(
    formula    = ~ fac,
    parameter  = "mu",
    data       = data.frame(
      fac = factor(c("A", "B", "C"), levels = c("A", "B", "C"))
    ),
    prior_list = list(
      intercept = prior_spike_and_slab(
        prior("normal", list(mean = 1, sd = 0.2)),
        prior_inclusion = prior("point", list(location = 0.5))
      ),
      fac = prior_spike_and_slab(
        prior_factor("normal", list(mean = 0.5, sd = 0.2), contrast = "treatment"),
        prior_inclusion = prior("point", list(location = 0.5))
      )
    )
  )

  indicators <- expand.grid(
    mu_intercept_indicator = c(0, 1),
    mu_fac_indicator       = c(0, 1)
  )
  indicators <- indicators[rep(seq_len(nrow(indicators)), each = 50), ]
  posterior <- cbind(
    mu_intercept = ifelse(
      indicators[["mu_intercept_indicator"]] == 1,
      seq(0.75, 1.25, length.out = nrow(indicators)),
      0
    ),
    "mu_fac[1]" = ifelse(
      indicators[["mu_fac_indicator"]] == 1,
      seq(0.25, 0.75, length.out = nrow(indicators)),
      0
    ),
    "mu_fac[2]" = ifelse(
      indicators[["mu_fac_indicator"]] == 1,
      seq(0.50, 1.00, length.out = nrow(indicators)),
      0
    ),
    indicators
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result[["prior_list"]]

  inference <- as_marginal_inference(
    model                = fit,
    marginal_parameters = "mu_fac",
    parameters          = c("mu_intercept", "mu_fac"),
    conditional_list    = list(mu_fac = c("mu_intercept", "mu_fac")),
    conditional_rule    = "OR",
    formula             = ~ fac,
    n_samples           = 128,
    silent              = TRUE
  )

  factor_levels <- inference[["conditional"]][["mu_fac"]]
  averaged_levels <- inference[["averaged"]][["mu_fac"]]

  expect_equal(attr(factor_levels[["A"]], "effective_conditional"), "mu_intercept")
  expect_equal(attr(factor_levels[["B"]], "effective_conditional"), c("mu_intercept", "mu_fac"))
  expect_equal(attr(factor_levels[["C"]], "effective_conditional"), c("mu_intercept", "mu_fac"))
  expect_equal(mean(as.numeric(factor_levels[["A"]]) == 0), 0)
  expect_equal(mean(as.numeric(factor_levels[["B"]]) == 0), 0)
  expect_equal(mean(as.numeric(factor_levels[["C"]]) == 0), 0)
  expect_equal(.prior_linear_density_point_mass(attr(factor_levels[["A"]], "prior_density"), 0), 0)
  expect_equal(.prior_linear_density_point_mass(attr(factor_levels[["B"]], "prior_density"), 0), 0)
  expect_equal(.prior_linear_density_point_mass(attr(factor_levels[["C"]], "prior_density"), 0), 0)
  expect_false(length(factor_levels[["A"]]) == length(factor_levels[["B"]]))

  averaged_plot_data <- .plot_data_marginal_samples(
    samples                  = inference[["averaged"]],
    parameter                = "mu_fac",
    prior                    = TRUE,
    n_points                 = 32,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE
  )
  averaged_points <- averaged_plot_data[vapply(averaged_plot_data, inherits, logical(1), what = "density.prior.point")]
  averaged_density <- averaged_plot_data[vapply(averaged_plot_data, inherits, logical(1), what = "density.prior.simple")]

  expect_equal(length(averaged_points), 3L)
  expect_equal(length(averaged_density), 3L)
  for(level in names(averaged_levels)){
    level_points <- averaged_points[vapply(averaged_points, function(x) identical(attr(x, "level_name"), level), logical(1))]
    level_density <- averaged_density[vapply(averaged_density, function(x) identical(attr(x, "level_name"), level), logical(1))]

    expect_equal(length(level_points), 1L)
    expect_equal(length(level_density), 1L)
    expect_equal(level_points[[1]][["x"]], 0)
    expect_equal(level_points[[1]][["y"]], mean(as.numeric(averaged_levels[[level]]) == 0))
    expect_false(any(level_density[[1]][["samples"]] == 0))
  }

  expect_warning(
    marginal_estimates_table(
      samples    = inference[["conditional"]],
      inference  = inference[["inference"]],
      parameters = "mu_fac"
    ),
    NA
  )
  expect_warning(
    .plot_data_marginal_samples(
      samples                 = inference[["conditional"]],
      parameter               = "mu_fac",
      prior                   = TRUE,
      n_points                = 32,
      transformation          = NULL,
      transformation_arguments = NULL,
      transformation_settings = FALSE
    ),
    NA
  )

  inference_factor_only <- as_marginal_inference(
    model                = fit,
    marginal_parameters = "mu_fac",
    parameters          = c("mu_intercept", "mu_fac"),
    conditional_list    = list(mu_fac = "mu_fac"),
    conditional_rule    = "OR",
    formula             = ~ fac,
    n_samples           = 128,
    silent              = TRUE
  )

  factor_only_levels <- inference_factor_only[["conditional"]][["mu_fac"]]

  expect_equal(attr(factor_only_levels[["A"]], "effective_conditional"), character())
  expect_equal(attr(factor_only_levels[["B"]], "effective_conditional"), "mu_fac")
  expect_equal(attr(factor_only_levels[["C"]], "effective_conditional"), "mu_fac")
})
