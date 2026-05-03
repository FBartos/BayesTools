skip_if_not_test_profile("unit")

test_that("as_mixed_posteriors handles treatment factor-continuous interaction coefficients", {

  df <- data.frame(
    alloc = factor(
      rep(c("alternate", "random", "systematic"), each = 4),
      levels = c("alternate", "random", "systematic")
    ),
    year = seq(1960, 1971, length.out = 12)
  )
  formula_result <- JAGS_formula(
    formula       = ~ alloc * year,
    parameter     = "mu",
    data          = df,
    formula_scale = list(year = TRUE),
    prior_list    = list(
      intercept    = prior("normal", list(0, 1)),
      alloc        = prior_factor("normal", list(0, 1), contrast = "treatment"),
      year         = prior("normal", list(0, 1)),
      "alloc:year" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_alloc__xXx__year

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_alloc__xXx__year[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list")    <- formula_result$prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  samples <- as_mixed_posteriors(
    model            = fit,
    parameters       = "mu_alloc__xXx__year",
    transform_scaled = TRUE,
    n_prior_samples  = 64
  )

  expect_equal(
    colnames(samples$mu_alloc__xXx__year),
    c(
      "mu_alloc[random]__xXx__year",
      "mu_alloc[systematic]__xXx__year"
    )
  )
  expect_equal(attr(samples$mu_alloc__xXx__year, "factor_design"), attr(interaction_prior, "factor_design"))

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples                  = samples,
    parameter                = "mu_alloc__xXx__year",
    n_points                 = 32,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE
  )
  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.factor")]

  expect_equal(length(density_entries), 2L)
  expect_equal(
    unname(vapply(density_entries, attr, character(1), which = "level_name")),
    colnames(samples$mu_alloc__xXx__year)
  )

  prior_plot_data <- BayesTools:::.plot_data_prior_factor_density_transformed(
    prior_density_context = attr(samples, "prior_density_context"),
    samples               = samples,
    parameter             = "mu_alloc__xXx__year",
    prior_list            = attr(samples$mu_alloc__xXx__year, "prior_list"),
    n_points              = 32
  )

  expect_equal(length(prior_plot_data), 2L)
  expect_equal(
    unname(vapply(prior_plot_data, attr, character(1), which = "level_name")),
    colnames(samples$mu_alloc__xXx__year)
  )
})


test_that("marginal_posterior handles treatment factor-continuous interaction coefficients", {

  df <- data.frame(
    alloc = factor(
      rep(c("alternate", "random", "systematic"), each = 4),
      levels = c("alternate", "random", "systematic")
    ),
    year = seq(1960, 1971, length.out = 12)
  )
  formula_result <- JAGS_formula(
    formula    = ~ alloc * year,
    parameter  = "mu",
    data       = df,
    prior_list = list(
      intercept    = prior("normal", list(0, 1)),
      alloc        = prior_factor("normal", list(0, 1), contrast = "treatment"),
      year         = prior("normal", list(0, 1)),
      "alloc:year" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_alloc__xXx__year

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_alloc__xXx__year[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list

  samples <- as_mixed_posteriors(
    model      = fit,
    parameters = "mu_alloc__xXx__year"
  )
  marginal <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_alloc__xXx__year",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 32
  )
  expected <- posterior %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(names(marginal), c("alternate", "random", "systematic"))
  expect_equal(as.numeric(marginal$alternate), as.numeric(expected[, 1]))
  expect_equal(as.numeric(marginal$random), as.numeric(expected[, 2]))
  expect_equal(as.numeric(marginal$systematic), as.numeric(expected[, 3]))
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(marginal$alternate, "prior_density"), 0),
    1
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(marginal$random, "prior_density"), 0),
    0
  )
})
