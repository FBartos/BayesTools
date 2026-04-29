# TEST FILE: Weightfunction analytical plot priors
# ============================================================================ #
#
# PURPOSE:
#   Numeric tests for analytical weightfunction prior plot data used by
#   plot_prior_list() and plot_posterior(..., prior = TRUE).
#
# TAGS: @evaluation, @plots, @weightfunctions
# ============================================================================ #

test_that("mixed weightfunction plot data uses analytical cumulative-Dirichlet marginals", {

  wf_prior <- prior_weightfunction("one.sided", list(c(.05), c(2, 4)), prior_weights = 3)
  no_bias  <- prior_none(prior_weights = 1)

  plot_data_1 <- .plot_data_prior_list.weightfunction(
    list(no_bias, wf_prior),
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 50,
    n_samples = 1
  )
  plot_data_2 <- .plot_data_prior_list.weightfunction(
    list(no_bias, wf_prior),
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 50,
    n_samples = 10000
  )

  expect_null(plot_data_1$samples)
  expect_equal(plot_data_1$y, plot_data_2$y)
  expect_equal(plot_data_1$y_lCI, plot_data_2$y_lCI)
  expect_equal(plot_data_1$y_uCI, plot_data_2$y_uCI)

  expect_equal(plot_data_1$x, c(0, .05, .05, 1))
  expect_equal(plot_data_1$y, c(1, 1, .75, .75), tolerance = 1e-8)
  expect_equal(plot_data_1$y_lCI, c(1, 1, stats::qbeta(.025 / .75, 4, 2), stats::qbeta(.025 / .75, 4, 2)), tolerance = 1e-7)
  expect_equal(plot_data_1$y_uCI, c(1, 1, 1, 1), tolerance = 1e-8)
})

test_that("individual omega prior plot data is analytical for mixtures and points", {

  wf_prior <- prior_weightfunction("one.sided", list(c(.05), c(2, 4)), prior_weights = 3)
  no_bias  <- prior_none(prior_weights = 1)

  plot_data <- .plot_data_prior_list.weightparameter(
    list(no_bias, wf_prior),
    parameter = "omega[0.05,1]",
    n_points = 5,
    n_samples = 1
  )

  expect_named(plot_data, c("density", "points1"))
  expect_null(plot_data$density$samples)
  expect_equal(plot_data$density$x, seq(0, 1, length.out = 5))
  expect_equal(plot_data$density$y, .75 * c(0, stats::dbeta(c(.25, .5, .75), 4, 2), 0), tolerance = 1e-8)
  expect_equal(plot_data$points1$x, 1)
  expect_equal(plot_data$points1$y, .25)
})

test_that("mixed omega prior components preserve total probability mass", {

  prior_list <- list(
    prior_none(prior_weights = 1),
    prior_weightfunction("two.sided", list(alpha = c(1, 1), steps = c(.05)), prior_weights = 1/3),
    prior_weightfunction("one.sided", list(alpha = c(1, 1, 1), steps = c(.025, .05)), prior_weights = 1/3),
    prior_PET("normal", list(0, 1), prior_weights = 1/3)
  )

  context <- .weightfunction_prior_list_context(prior_list)
  expected_point_mass <- c(1, 2/3, 2/3, 5/6)
  names(expected_point_mass) <- context$omega_names

  for(parameter in context$omega_names){
    plot_data <- .plot_data_prior_list.weightparameter(
      prior_list,
      parameter = parameter,
      n_points = 1001,
      n_samples = 1
    )

    density_mass <- if("density" %in% names(plot_data)){
      sum(diff(plot_data$density$x) * (head(plot_data$density$y, -1) + tail(plot_data$density$y, -1)) / 2)
    }else{
      0
    }
    point_mass <- sum(vapply(plot_data, function(component){
      if(inherits(component, "density.prior.point")) component$y else 0
    }, numeric(1)))

    expect_equal(point_mass, unname(expected_point_mass[parameter]), tolerance = 1e-8)
    expect_equal(density_mass + point_mass, 1, tolerance = 1e-6)
  }
})

test_that("posterior omega parameter densities are scaled by all samples", {

  prior_list <- list(
    prior_none(prior_weights = 3),
    prior_weightfunction("one.sided", list(c(.05), c(1, 1)), prior_weights = 1)
  )

  samples <- matrix(
    c(rep(1, 3000), seq(.001, .999, length.out = 1000)),
    ncol = 1
  )
  colnames(samples) <- "omega[0.05,1]"
  attr(samples, "prior_list") <- prior_list
  attr(samples, "models_ind") <- c(rep(1, 3000), rep(2, 1000))

  plot_data <- .plot_data_samples.weightparameter(
    list(omega = samples),
    parameter = "omega[0.05,1]",
    n_points = 512
  )

  density_mass <- sum(diff(plot_data$density$x) * (head(plot_data$density$y, -1) + tail(plot_data$density$y, -1)) / 2)

  expect_equal(plot_data$points1$y, .75)
  expect_gt(density_mass, .20)
  expect_lt(density_mass, .30)
})

test_that("conditional bias posteriors zero null bias prior weights", {

  bias_prior <- prior_mixture(
    list(
      prior_none(prior_weights = 1),
      prior_weightfunction("one.sided", list(c(.05), c(1, 1)), prior_weights = 1),
      prior_PET("normal", list(0, 1), prior_weights = 1)
    ),
    is_null = c(TRUE, FALSE, FALSE)
  )

  model <- matrix(
    c(
      1, 1.0, 0.0, 0.0,
      2, 0.2, 0.5, 0.0,
      3, 0.3, 1.0, 0.7
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(model) <- c("bias_indicator", "omega[2]", "omega[1]", "PET")
  class(model) <- c("matrix", "BayesTools_fit")
  attr(model, "prior_list") <- list(bias = bias_prior)

  mixed <- as_mixed_posteriors(model, parameters = "bias", conditional = "bias")
  conditioned_prior_weights <- sapply(attr(mixed$bias, "prior_list"), function(prior) prior$prior_weights)

  expect_equal(conditioned_prior_weights, c(0, 1, 1))
  expect_equal(attr(mixed$bias, "models_ind"), c(2, 3))
})

test_that("non-monotone one-sided weightfunction marginals use product-beta CDFs", {

  wf_prior <- prior_weightfunction("one.sided", list(c(.05, .60), c(1, 1), c(1, 1)))
  context  <- .weightfunction_prior_list_context(list(wf_prior))

  expect_equal(context$omega_names, c("omega[0,0.05]", "omega[0.05,0.6]", "omega[0.6,1]"))

  components <- lapply(seq_along(context$omega_names), function(i){
    .weightfunction_prior_marginal_components(context, i)[[1]]
  })

  expect_equal(vapply(components, `[[`, character(1), "type"), c("point", "beta", "one_minus_product_beta"))
  expect_equal(vapply(components, .weightfunction_component_mean, numeric(1)), c(1, .5, .75), tolerance = 1e-8)

  product_cdf <- 1 - .5 + .5 * log(.5)
  expect_equal(.weightfunction_component_cdf(components[[3]], .5), product_cdf, tolerance = 1e-6)

  plot_data_1 <- .plot_data_prior_list.weightfunction(
    list(wf_prior),
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 20,
    n_samples = 1
  )
  plot_data_2 <- .plot_data_prior_list.weightfunction(
    list(wf_prior),
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 20,
    n_samples = 10000
  )

  expect_null(plot_data_1$samples)
  expect_equal(plot_data_1$y, plot_data_2$y)
  expect_equal(plot_data_1$y_lCI, plot_data_2$y_lCI)
  expect_equal(plot_data_1$y_uCI, plot_data_2$y_uCI)
})
