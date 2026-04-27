# ============================================================================ #
# TEST FILE: Model Averaging Plots Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for plot functions including input validation and
#   error handling for invalid prior configurations.
#
# DEPENDENCIES:
#   - None (pure R testing)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @edge-cases, @plots, @input-validation
# ============================================================================ #


# ============================================================================ #
# SECTION 1: plot_prior_list input validation
# ============================================================================ #
test_that("plot_prior_list rejects non-list input", {

  expect_error(
    plot_prior_list(prior("normal", list(0, 1))),
    "must be a list of priors"
  )

})


test_that("plot_prior_list rejects PET-PEESE without prior_list_mu", {

  pet_list <- list(
    p1 = prior_PET("normal", list(0, 1))
  )
  expect_error(
    plot_prior_list(pet_list),
    "prior_list_mu"
  )

})


test_that("plot_prior_list rejects prior_list_mu when not needed", {

  simple_list <- list(
    p1 = prior("normal", list(0, 1))
  )
  expect_error(
    plot_prior_list(simple_list, prior_list_mu = list(prior("spike", list(0)))),
    "prior_list_mu"
  )

})

.scaled_multi_factor_interaction_samples <- function(n_posterior = 20, n_prior = 128){

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3")),
    x = seq_len(4)
  )
  formula_result <- JAGS_formula(
    formula = ~ x * a * b,
    parameter = "mu",
    data = df,
    formula_scale = list(x = TRUE),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      a = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "x:a" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x:b" = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x:a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  prior_list <- formula_result$prior_list
  posterior_columns <- unlist(lapply(names(prior_list), function(parameter) {
    if(is.prior.factor(prior_list[[parameter]])){
      BayesTools:::.JAGS_prior_factor_names(parameter, prior_list[[parameter]])
    }else{
      parameter
    }
  }), use.names = FALSE)

  set.seed(515)
  posterior <- matrix(rnorm(n_posterior * length(posterior_columns)), nrow = n_posterior)
  colnames(posterior) <- posterior_columns
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  samples <- as_mixed_posteriors(
    fit,
    parameters = "mu_x__xXx__a__xXx__b",
    transform_scaled = TRUE,
    n_prior_samples = n_prior
  )

  list(samples = samples, prior_list = prior_list)
}

.scaled_default_treatment_interaction_samples <- function(n_posterior = 20, n_prior = 128){

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
      "__default_continuous" = prior("normal", list(0, 1)),
      "__default_factor"     = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  prior_list <- formula_result$prior_list
  posterior_columns <- unlist(lapply(names(prior_list), function(parameter) {
    if(is.prior.factor(prior_list[[parameter]])){
      BayesTools:::.JAGS_prior_factor_names(parameter, prior_list[[parameter]])
    }else{
      parameter
    }
  }), use.names = FALSE)

  set.seed(616)
  posterior <- matrix(rnorm(n_posterior * length(posterior_columns)), nrow = n_posterior)
  colnames(posterior) <- posterior_columns
  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  samples <- as_mixed_posteriors(
    fit,
    parameters = "mu_alloc",
    transform_scaled = TRUE,
    n_prior_samples = n_prior
  )

  list(samples = samples, prior_list = prior_list)
}

test_that("transformed factor prior plot data uses indexed prior density columns", {

  fixture <- .scaled_multi_factor_interaction_samples()
  samples <- fixture$samples
  prior_list <- fixture$prior_list

  plot_data <- BayesTools:::.plot_data_prior_factor_density_transformed(
    prior_density_context = attr(samples, "prior_density_context"),
    samples = samples,
    parameter = "mu_x__xXx__a__xXx__b",
    prior_list = attr(samples$mu_x__xXx__a__xXx__b, "prior_list"),
    n_points = 32
  )

  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.factor")]
  expect_equal(length(density_entries), 6L)
  expect_equal(
    unname(vapply(density_entries, attr, character(1), which = "level_name")),
    paste0(
      "mu_x__xXx__a[dif: ",
      rep(c("a1", "a2"), times = 3),
      "]__xXx__b[dif: ",
      rep(c("b1", "b2", "b3"), each = 2),
      "]"
    )
  )
})

test_that("transformed treatment factor prior plot data omits reference coefficient", {

  fixture <- .scaled_default_treatment_interaction_samples()
  samples <- fixture$samples

  factor_weights <- BayesTools:::.prior_factor_level_weight_matrix(
    sample_metadata = samples$mu_alloc,
    parameter       = "mu_alloc",
    samples         = samples
  )
  expect_equal(rownames(factor_weights), colnames(samples$mu_alloc))
  expect_equal(colnames(factor_weights), c("mu_alloc[1]", "mu_alloc[2]"))
  expect_false(any(rowSums(abs(factor_weights)) == 0))

  plot_data <- BayesTools:::.plot_data_prior_factor_density_transformed(
    prior_density_context = attr(samples, "prior_density_context"),
    samples               = samples,
    parameter             = "mu_alloc",
    prior_list            = attr(samples$mu_alloc, "prior_list"),
    n_points              = 32
  )

  expect_equal(length(plot_data), 2L)
  expect_false(any(vapply(plot_data, inherits, logical(1), what = "density.prior.point")))
  expect_equal(
    unname(vapply(plot_data, attr, character(1), which = "level_name")),
    colnames(samples$mu_alloc)
  )
})

test_that("plot_posterior handles transformed scaled factor-interaction priors", {

  skip_if_not_installed("ggplot2")

  fixture <- .scaled_multi_factor_interaction_samples(n_posterior = 24, n_prior = 128)

  plot <- plot_posterior(
    samples = fixture$samples,
    parameter = "mu_x__xXx__a__xXx__b",
    plot_type = "ggplot",
    prior = TRUE,
    n_points = 32
  )

  expect_s3_class(plot, "ggplot")
  expect_gt(length(plot$layers), 1)
})
