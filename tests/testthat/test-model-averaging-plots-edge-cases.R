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

test_that("PET-PEESE prior plot data uses deterministic linear-combination summaries", {

  prior_list <- list(
    PET = prior_PET("normal", list(0, 1), truncation = list(-Inf, Inf))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )

  plot_data <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = c(0, 0.5, 1),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 3,
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )

  expect_null(plot_data$samples)
  expect_equal(plot_data$y, c(0, 0, 0), tolerance = 1e-8)
  expect_equal(plot_data$y_lCI, stats::qnorm(.025) * c(0, 0.5, 1), tolerance = 0.02)
  expect_equal(plot_data$y_uCI, stats::qnorm(.975) * c(0, 0.5, 1), tolerance = 0.02)
})

test_that("PET-PEESE prior plot data uses CDF quantiles for half-Cauchy PET", {

  x_seq <- c(0, 0.25, 0.5, 1)
  prior_list <- list(
    PET = prior_PET("cauchy", list(0, 1))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )

  positive <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu,
    effect_direction         = "positive"
  )

  expected_positive <- function(p) x_seq * tan(pi * p / 2)
  expect_null(positive$samples)
  expect_equal(positive$y_lCI, expected_positive(.025), tolerance = 1e-4)
  expect_equal(positive$y,     expected_positive(.500), tolerance = 1e-4)
  expect_equal(positive$y_uCI, expected_positive(.975), tolerance = 1e-4)

  negative <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu,
    effect_direction         = "negative"
  )

  expect_equal(negative$y_lCI, -positive$y_uCI, tolerance = 1e-4)
  expect_equal(negative$y,     -positive$y,     tolerance = 1e-4)
  expect_equal(negative$y_uCI, -positive$y_lCI, tolerance = 1e-4)
})

test_that("PET-PEESE prior plot data uses CDF quantiles for PET and PEESE mixtures", {

  x_seq <- c(0, 0.25, 0.5, 1)
  prior_list <- list(
    PET   = prior_PET("cauchy", list(0, 1)),
    PEESE = prior_PEESE("cauchy", list(0, 5))
  )
  prior_list_mu <- list(
    mu1 = prior("spike", list(0)),
    mu2 = prior("spike", list(0))
  )

  half_cauchy_mixture_quantile <- function(se, p){
    if(se == 0){
      return(0)
    }

    a <- se
    b <- 5 * se^2
    if(isTRUE(all.equal(p, .5))){
      return(sqrt(a * b))
    }

    tangent <- tan(pi * p)
    (-(a + b) + sign(tangent) * sqrt((a + b)^2 + 4 * tangent^2 * a * b)) /
      (2 * tangent)
  }

  plot_data <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )

  expect_null(plot_data$samples)
  expect_equal(plot_data$y_lCI, vapply(x_seq, half_cauchy_mixture_quantile, numeric(1), p = .025), tolerance = 5e-4)
  expect_equal(plot_data$y,     vapply(x_seq, half_cauchy_mixture_quantile, numeric(1), p = .500), tolerance = 5e-4)
  expect_equal(plot_data$y_uCI, vapply(x_seq, half_cauchy_mixture_quantile, numeric(1), p = .975), tolerance = 5e-4)

  dense <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = seq(0, 1, length.out = 101),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 101,
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )
  expect_lt(max(abs(diff(diff(dense$y_uCI)))), 0.02)
})

test_that("PET-PEESE prior plot data transforms ordered CDF quantiles", {

  prior_list <- list(
    PET = prior_PET("normal", list(0, 1), truncation = list(-Inf, Inf))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )
  x_seq <- c(0, 0.5, 1)

  exp_plot <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = "exp",
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )
  expect_equal(exp_plot$y,     exp(rep(0, length(x_seq))), tolerance = 1e-8)
  expect_equal(exp_plot$y_lCI, exp(stats::qnorm(.025) * x_seq), tolerance = 1e-4)
  expect_equal(exp_plot$y_uCI, exp(stats::qnorm(.975) * x_seq), tolerance = 1e-4)

  lin_plot <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = "lin",
    transformation_arguments = list(a = 2, b = -3),
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )
  expect_equal(lin_plot$y, rep(2, length(x_seq)), tolerance = 1e-8)
  expect_equal(lin_plot$y_lCI, 2 - 3 * stats::qnorm(.975) * x_seq, tolerance = 1e-4)
  expect_equal(lin_plot$y_uCI, 2 - 3 * stats::qnorm(.025) * x_seq, tolerance = 1e-4)
})

test_that("PET-PEESE prior plot data falls back to samples for custom transformations", {

  set.seed(1)
  prior_list <- list(
    PET = prior_PET("normal", list(0, 1), truncation = list(-Inf, Inf))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )
  custom_transformation <- list(
    fun = function(x) x^2,
    inv = sqrt,
    jac = function(x) 1 / (2 * sqrt(x))
  )

  plot_data <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = c(0, 0.5, 1),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 3,
    n_samples                = 100,
    transformation           = custom_transformation,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )

  expect_false(is.null(plot_data$samples))
  expect_equal(dim(plot_data$samples), c(100, 3))
  expect_true(all(plot_data$y_lCI <= plot_data$y))
  expect_true(all(plot_data$y <= plot_data$y_uCI))
})

test_that("PET-PEESE posterior plot data does not recycle coefficient rows", {

  make_samples <- function(x, prior) {
    attr(x, "prior_list") <- list(prior)
    attr(x, "models_ind") <- rep(1, length(x))
    class(x) <- c("mixed_posteriors", "mixed_posteriors.simple")
    x
  }

  samples <- list(
    mu    = make_samples(c(0, 1), prior("normal", list(0, 1))),
    PET   = make_samples(c(2, 3), prior_PET("normal", list(0, 1))),
    PEESE = make_samples(c(0, 0), prior("point", list(0)))
  )
  class(samples) <- c("mixed_posteriors", "list")

  plot_data <- BayesTools:::.plot_data_samples.PETPEESE(
    samples                  = samples,
    x_seq                    = c(0, 1),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 2,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE
  )

  expect_equal(dim(plot_data$samples), c(2, 2))
  expect_equal(plot_data$y, c(0.5, 3))
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

  expect_silent({
    plot <- plot_posterior(
      samples = fixture$samples,
      parameter = "mu_x__xXx__a__xXx__b",
      plot_type = "ggplot",
      prior = TRUE,
      n_points = 32
    )
  })

  expect_s3_class(plot, "ggplot")
  expect_gt(length(plot$layers), 1)
})
