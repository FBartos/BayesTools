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
