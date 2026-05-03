skip_if_not_test_profile("unit")

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
