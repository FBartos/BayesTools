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

  wf_prior <- prior_weightfunction("one-sided", c(.05), wf_cumulative(c(2, 4)), prior_weights = 3)
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

  wf_prior <- prior_weightfunction("one-sided", c(.05), wf_cumulative(c(2, 4)), prior_weights = 3)
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
    prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)), prior_weights = 1/3),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 1, 1)), prior_weights = 1/3),
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
    prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1)), prior_weights = 1)
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
      prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1)), prior_weights = 1),
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
  colnames(model) <- c("bias_indicator", "omega[1]", "omega[2]", "PET")
  class(model) <- c("matrix", "BayesTools_fit")
  attr(model, "prior_list") <- list(bias = bias_prior)

  mixed <- as_mixed_posteriors(model, parameters = "bias", conditional = "bias")
  conditioned_prior_weights <- sapply(attr(mixed$bias, "prior_list"), function(prior) prior$prior_weights)

  expect_equal(conditioned_prior_weights, c(0, 1, 1))
  expect_equal(attr(mixed$bias, "models_ind"), c(2, 3))
})

test_that("independent and log-independent weightfunction marginals are analytical", {

  beta_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("beta", list(2, 3))),
    prior_weights = 1
  )
  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(0, 1), truncation = list(lower = -Inf, upper = 0)), "log_omega"),
    prior_weights = 1
  )

  beta_component <- .weightfunction_prior_marginal_components(
    .weightfunction_prior_list_context(list(beta_prior)), 2
  )[[1]]
  log_component <- .weightfunction_prior_marginal_components(
    .weightfunction_prior_list_context(list(log_prior)), 2
  )[[1]]

  expect_equal(beta_component$type, "prior")
  expect_equal(.weightfunction_component_pdf(beta_component, .5), stats::dbeta(.5, 2, 3), tolerance = 1e-8)
  expect_equal(.weightfunction_component_cdf(beta_component, .5), stats::pbeta(.5, 2, 3), tolerance = 1e-8)

  expect_equal(log_component$type, "prior")
  expect_equal(
    .weightfunction_component_pdf(log_component, .5),
    pdf(log_component$prior, log(.5)) / .5,
    tolerance = 1e-8
  )
  expect_equal(
    .weightfunction_component_cdf(log_component, .5),
    mcdf(log_component$prior, log(.5)),
    tolerance = 1e-8
  )
})

test_that("analytical weightfunction plotting supports independent weights above one", {

  omega_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3))),
    prior_weights = 1
  )
  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(mean = log(2), sd = .2)), "log_omega"),
    prior_weights = 1
  )

  omega_component <- .weightfunction_prior_marginal_components(
    .weightfunction_prior_list_context(list(omega_prior)), 2
  )[[1]]
  log_component <- .weightfunction_prior_marginal_components(
    .weightfunction_prior_list_context(list(log_prior)), 2
  )[[1]]

  expect_equal(.weightfunction_component_pdf(omega_component, 1.5), stats::dgamma(1.5, shape = 9, rate = 3), tolerance = 1e-8)
  expect_equal(.weightfunction_component_cdf(omega_component, 1.5), stats::pgamma(1.5, shape = 9, rate = 3), tolerance = 1e-8)
  expect_equal(.weightfunction_component_pdf(log_component, 2), stats::dnorm(log(2), mean = log(2), sd = .2) / 2, tolerance = 1e-8)
  expect_equal(.weightfunction_component_cdf(log_component, 2), stats::pnorm(log(2), mean = log(2), sd = .2), tolerance = 1e-8)
  expect_equal(mlpdf(log_prior, 2)[,2], stats::dnorm(log(2), mean = log(2), sd = .2, log = TRUE) - log(2), tolerance = 1e-8)
  expect_equal(mcdf(log_prior, 2)[,2], .5, tolerance = 1e-8)
  expect_gt(mquant(log_prior, .75)[,2], 1)
  expect_gt(range(log_prior)[2], 1)

  omega_parameter_data <- .plot_data_prior_list.weightparameter(
    list(omega_prior),
    parameter = "omega[0.05,1]",
    n_points = 100,
    n_samples = 1
  )
  log_parameter_data <- .plot_data_prior_list.weightparameter(
    list(log_prior),
    parameter = "omega[0.05,1]",
    n_points = 100,
    n_samples = 1
  )
  log_weightfunction_data <- .plot_data_prior_list.weightfunction(
    list(log_prior),
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 100,
    n_samples = 1
  )

  expect_gt(max(omega_parameter_data$density$x), 1)
  expect_gt(max(log_parameter_data$density$x), 1)
  expect_gt(attr(omega_parameter_data$density, "x_range")[2], 1)
  expect_gt(attr(log_parameter_data$density, "x_range")[2], 1)
  expect_gt(log_weightfunction_data$y[3], 1)
  expect_gt(log_weightfunction_data$y_uCI[3], 1)
  expect_gt(attr(log_weightfunction_data, "y_range")[2], 1)

  log_density <- density(log_prior, individual = TRUE, n_points = 100)
  log_weightfunction_density <- density(log_prior, n_points = 100)

  expect_gt(max(log_density[[2]]$x), 1)
  expect_gt(attr(log_density[[2]], "x_range")[2], 1)
  expect_gt(attr(log_weightfunction_density, "y_range")[2], 1)
  expect_s3_class(plot(log_prior, plot_type = "ggplot"), "ggplot")
})

test_that("posterior weightfunction plotting ranges include omega samples above one", {

  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"),
    prior_weights = 1
  )

  omega_samples <- matrix(
    c(rep(1, 200), seq(1.05, 2.25, length.out = 200)),
    ncol = 2
  )
  colnames(omega_samples) <- c("omega[0,0.05]", "omega[0.05,1]")
  attr(omega_samples, "prior_list") <- list(log_prior)
  attr(omega_samples, "models_ind") <- rep(1, nrow(omega_samples))

  parameter_data <- .plot_data_samples.weightparameter(
    list(omega = omega_samples),
    parameter = "omega[0.05,1]",
    n_points = 100
  )
  weightfunction_data <- .plot_data_samples.weightfunction(
    list(omega = omega_samples),
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 100
  )

  expect_gt(max(parameter_data$density$x), 1)
  expect_gt(attr(parameter_data$density, "x_range")[2], 1)
  expect_gt(weightfunction_data$y[3], 1)
  expect_gt(weightfunction_data$y_uCI[3], 1)
  expect_gt(attr(weightfunction_data, "y_range")[2], 1)
})

test_that("analytical plotting handles heterogeneous weightfunction mixtures", {

  prior_list <- list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 2),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 3),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 4),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 5),
    prior_PET("normal", list(0, 1), prior_weights = 6),
    prior_PEESE("normal", list(0, 1), prior_weights = 7)
  )

  context <- .weightfunction_prior_list_context(prior_list)
  expect_equal(context$omega_names, c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.1]", "omega[0.1,0.975]", "omega[0.975,1]"))

  parameter_ind <- match("omega[0.05,0.1]", context$omega_names)
  components <- .weightfunction_prior_marginal_components(context, parameter_ind)
  weights <- vapply(components, `[[`, numeric(1), "weight")

  expected_cdf_1_5 <- sum(weights * c(
    1,
    stats::pbeta(1.5, 3, 3),
    stats::pgamma(1.5, shape = 9, rate = 3),
    stats::pnorm(log(1.5), mean = log(1.5), sd = .15),
    1,
    1,
    1
  ))
  expected_pdf_1_5 <- sum(weights * c(
    0,
    0,
    stats::dgamma(1.5, shape = 9, rate = 3),
    stats::dnorm(log(1.5), mean = log(1.5), sd = .15) / 1.5,
    0,
    0,
    0
  ))

  expect_equal(.weightfunction_mixture_cdf(components, 1.5), expected_cdf_1_5, tolerance = 1e-8)
  expect_equal(
    sum(vapply(components, function(component) component$weight * .weightfunction_component_pdf(component, 1.5), numeric(1))),
    expected_pdf_1_5,
    tolerance = 1e-8
  )
  expect_gt(.weightfunction_mixture_mean(components), 1)
  expect_gt(.weightfunction_mixture_quantile(components, .975), 1)

  parameter_plot_data <- .plot_data_prior_list.weightparameter(
    prior_list,
    parameter = "omega[0.05,0.1]",
    n_points = 400,
    n_samples = 1
  )
  point_mass <- sum(vapply(parameter_plot_data, function(component){
    if(inherits(component, "density.prior.point")) component$y else 0
  }, numeric(1)))

  expect_null(parameter_plot_data$density$samples)
  expect_gt(max(parameter_plot_data$density$x), 1)
  expect_equal(point_mass, (1 + 5 + 6 + 7) / 28, tolerance = 1e-8)

  weightfunction_plot_data <- .plot_data_prior_list.weightfunction(
    prior_list,
    x_seq = NULL,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = 400,
    n_samples = 1
  )

  expect_null(weightfunction_plot_data$samples)
  expect_gt(max(weightfunction_plot_data$y), 1)
  expect_gt(max(weightfunction_plot_data$y_uCI), 1)
  expect_gt(attr(weightfunction_plot_data, "y_range")[2], 1)

  model <- matrix(
    c(
      1, 1.00, 1.00, 1.00, 1.00, 1.00, 0, 0,
      2, 1.00, 0.70, 0.50, 0.50, 0.50, 0, 0,
      3, 1.00, 1.00, 2.00, 2.50, 2.50, 0, 0,
      4, 1.00, 1.80, 1.80, 1.80, 1.80, 0, 0,
      5, 1.00, 0.40, 0.40, 0.40, 1.00, 0, 0,
      6, 1.00, 1.00, 1.00, 1.00, 1.00, 0, 0,
      7, 1.00, 1.00, 1.00, 1.00, 1.00, 0, 0
    ),
    nrow = 7,
    byrow = TRUE
  )
  colnames(model) <- c("bias_indicator", paste0("omega[", 1:5, "]"), "PET", "PEESE")
  class(model) <- c("matrix", "BayesTools_fit")
  attr(model, "prior_list") <- list(bias = prior_mixture(prior_list))

  mixed <- as_mixed_posteriors(model, parameters = "bias", conditional = "omega")
  posterior_plot_data <- .plot_data_samples.weightparameter(
    list(bias = mixed$bias),
    parameter = "omega[0.05,0.1]",
    n_points = 100
  )

  expect_true(all(attr(mixed$bias, "models_ind") %in% 2:5))
  expect_gt(max(posterior_plot_data$density$x), 1)
  expect_s3_class(plot_posterior(mixed, "omega", prior = TRUE, plot_type = "ggplot"), "ggplot")
})
