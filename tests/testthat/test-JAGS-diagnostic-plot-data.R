skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Diagnostic Plot Data
# ============================================================================ #
#
# PURPOSE:
#   Unit-profile semantic tests for non-rendered diagnostic plot data. These
#   assert chain/iteration alignment, selected parameter columns, transformations,
#   and prepared density/trace/autocorrelation data without cached JAGS fits.
#
# TAGS: @unit, @diagnostics, @plot-data
# ============================================================================ #

.mock_diagnostics_fit <- function() {
  chain_1 <- cbind(
    theta = c(-1.0, -0.4, -0.1, 0.3, 0.7, 1.1),
    "beta[1]" = c(1, 2, 3, 4, 5, 6),
    "beta[2]" = c(6, 5, 4, 3, 2, 1)
  )
  chain_2 <- cbind(
    theta = c(-0.8, -0.2, 0.0, 0.5, 0.9, 1.3),
    "beta[1]" = c(2, 3, 4, 5, 6, 7),
    "beta[2]" = c(7, 6, 5, 4, 3, 2)
  )

  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(
    theta = prior("normal", list(0, 1)),
    beta = prior("mnormal", list(mean = 0, sd = 1, K = 2))
  )
  fit
}

diagnostic_plot_layer_data <- function(plot) {
  ggplot2::ggplot_build(plot)$data
}

test_that("diagnostic plot data preserves chain, iteration, and transformations", {
  fit <- .mock_diagnostics_fit()
  prior_list <- attr(fit, "prior_list")

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "theta",
    prior_list = prior_list,
    transformations = list(theta = list(fun = function(x, shift) x + shift, arg = list(shift = 2))),
    transform_factors = FALSE
  )

  expect_equal(dim(plot_data), c(12, 1))
  expect_equal(colnames(plot_data), "theta")
  expect_equal(as.numeric(plot_data[1:6, "theta"]), c(-1.0, -0.4, -0.1, 0.3, 0.7, 1.1) + 2)
  expect_equal(as.numeric(plot_data[7:12, "theta"]), c(-0.8, -0.2, 0.0, 0.5, 0.9, 1.3) + 2)
  expect_equal(attr(plot_data, "chain"), rep(1:2, each = 6))
  expect_equal(attr(plot_data, "iter"), rep(1:6, times = 2))
  expect_equal(attr(plot_data, "parameter"), "theta")
  expect_identical(attr(plot_data, "prior"), prior_list$theta)

  trace_data <- .diagnostics_plot_data_trace(plot_data, n_points = 10, ylim = NULL)
  expect_equal(names(trace_data), "theta")
  expect_equal(attr(trace_data, "chains"), 2)
  expect_equal(attr(trace_data, "parameter_name"), "theta")
  expect_s3_class(trace_data$theta[[1]], "BayesTools_chain")
  expect_equal(trace_data$theta[[1]]$x, 1:6)
  expect_equal(trace_data$theta[[1]]$y, as.numeric(plot_data[1:6, "theta"]))
  expect_equal(trace_data$theta[[2]]$x, 1:6)
  expect_equal(trace_data$theta[[2]]$y, as.numeric(plot_data[7:12, "theta"]))

  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 32, xlim = c(0, 4))
  expect_equal(attr(density_data, "chains"), 2)
  expect_equal(attr(density_data$theta, "x_range"), c(0, 4))
  expect_s3_class(density_data$theta[[1]], "density")
  expect_equal(length(density_data$theta[[1]]$x), 32)
  expect_true(all(is.finite(density_data$theta[[1]]$y)))
  expect_true(all(density_data$theta[[1]]$y >= 0))

  autocorrelation_data <- .diagnostics_plot_data_autocorrelation(plot_data, n_points = 10, lags = 3)
  expect_equal(attr(autocorrelation_data, "chains"), 2)
  expect_s3_class(autocorrelation_data$theta[[1]], "BayesTools_autocorrelation")
  expect_equal(autocorrelation_data$theta[[1]]$x, 0:3)
  expect_equal(autocorrelation_data$theta[[1]]$y[[1]], 1)
  expect_equal(attr(autocorrelation_data$theta[[1]], "x_range"), c(0, 3))
})

test_that("diagnostic ggplot geoms expose exact density trace and autocorrelation data", {
  fit <- .mock_diagnostics_fit()
  prior_list <- attr(fit, "prior_list")

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "theta",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )

  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 32, xlim = c(-2, 2))
  density_plot <- .ggplot.prior_empty("simple", list(xlab = "theta", ylab = "Density")) +
    .geom_diagnostics.density(density_data$theta[[1]], col = "red", lwd = 2, lty = 2)
  density_layer <- diagnostic_plot_layer_data(density_plot)[[1]]

  expect_equal(density_layer$x, density_data$theta[[1]]$x)
  expect_equal(density_layer$y, density_data$theta[[1]]$y)
  expect_equal(unique(density_layer$colour), "red")
  expect_equal(unique(density_layer$linewidth), 2)
  expect_equal(unique(density_layer$linetype), 2)
  expect_equal(density_plot$scales$get_scales("x")$name, "theta")
  expect_equal(density_plot$scales$get_scales("y")$name, "Density")

  trace_data <- .diagnostics_plot_data_trace(plot_data, n_points = 10, ylim = c(-2, 2))
  trace_plot <- .ggplot.prior_empty("simple", list(xlab = "Iteration", ylab = "theta")) +
    .geom_diagnostics.trace(trace_data$theta[[2]], col = "blue", lwd = 3)
  trace_layer <- diagnostic_plot_layer_data(trace_plot)[[1]]

  expect_equal(trace_layer$x, 1:6)
  expect_equal(trace_layer$y, as.numeric(plot_data[attr(plot_data, "chain") == 2, "theta"]))
  expect_equal(unique(trace_layer$colour), "blue")
  expect_equal(unique(trace_layer$linewidth), 3)

  autocorrelation_data <- .diagnostics_plot_data_autocorrelation(plot_data, n_points = 10, lags = 3)
  autocorrelation_plot <- .ggplot.prior_empty("simple", list(xlab = "Lag", ylab = "Autocorrelation(theta)")) +
    .geom_diagnostics.autocorrelation(autocorrelation_data$theta[[1]], col = "black")
  autocorrelation_layer <- diagnostic_plot_layer_data(autocorrelation_plot)[[1]]

  expect_equal(autocorrelation_layer$x, 0:3)
  expect_equal(autocorrelation_layer$y, autocorrelation_data$theta[[1]]$y)
  expect_equal(unique(autocorrelation_layer$fill), "black")
  expect_equal(autocorrelation_plot$scales$get_scales("x")$name, "Lag")
  expect_equal(autocorrelation_plot$scales$get_scales("y")$name, "Autocorrelation(theta)")
})

test_that("diagnostic plot data selects vector prior columns in model order", {
  fit <- .mock_diagnostics_fit()
  prior_list <- attr(fit, "prior_list")

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "beta",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )

  expect_equal(dim(plot_data), c(12, 2))
  expect_equal(colnames(plot_data), c("beta[1]", "beta[2]"))
  expect_equal(as.numeric(plot_data[1:6, "beta[1]"]), 1:6)
  expect_equal(as.numeric(plot_data[1:6, "beta[2]"]), 6:1)
  expect_equal(as.numeric(plot_data[7:12, "beta[1]"]), 2:7)
  expect_equal(as.numeric(plot_data[7:12, "beta[2]"]), 7:2)
  expect_equal(attr(plot_data, "chain"), rep(1:2, each = 6))
  expect_equal(attr(plot_data, "parameter"), "beta")
  expect_identical(attr(plot_data, "prior"), prior_list$beta)

  trace_data <- .diagnostics_plot_data_trace(plot_data, n_points = 10, ylim = c(0, 8))
  expect_equal(names(trace_data), c("beta[1]", "beta[2]"))
  expect_equal(attr(trace_data$`beta[1]`, "y_range"), c(0, 8))
  expect_equal(trace_data$`beta[1]`[[1]]$y, 1:6)
  expect_equal(trace_data$`beta[2]`[[1]]$y, 6:1)
})

test_that("diagnostic plot data selects and names weightfunction omega bins", {
  chain_1 <- cbind(
    "omega[1]" = rep(1, 6),
    "omega[2]" = seq(.2, .7, length.out = 6),
    "omega[3]" = seq(.8, .3, length.out = 6)
  )
  chain_2 <- cbind(
    "omega[1]" = rep(1, 6),
    "omega[2]" = seq(.3, .8, length.out = 6),
    "omega[3]" = seq(.7, .2, length.out = 6)
  )
  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  prior_list <- list(
    omega = prior_weightfunction("one-sided", c(.025, .05), wf_fixed(c(1, .7, .3)))
  )

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "omega",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )

  expect_equal(dim(plot_data), c(12, 2))
  expect_equal(colnames(plot_data), c("omega[0.025,0.05]", "omega[0.05,1]"))
  expect_equal(as.numeric(plot_data[1:6, 1]), chain_1[, "omega[2]"])
  expect_equal(as.numeric(plot_data[1:6, 2]), chain_1[, "omega[3]"])
  expect_equal(as.numeric(plot_data[7:12, 1]), chain_2[, "omega[2]"])
  expect_equal(as.numeric(plot_data[7:12, 2]), chain_2[, "omega[3]"])
  expect_equal(attr(plot_data, "chain"), rep(1:2, each = 6))
  expect_equal(attr(plot_data, "iter"), rep(1:6, times = 2))
  expect_equal(attr(plot_data, "parameter"), "omega")
  expect_identical(attr(plot_data, "prior"), prior_list$omega)
})

test_that("diagnostic plot data masks spike-and-slab spike rows", {
  chain_1 <- cbind(
    theta = c(-2, -1, 0, 1, 2, 3, 4, 5),
    theta_indicator = c(0, 1, 1, 0, 1, 1, 1, 1)
  )
  chain_2 <- cbind(
    theta = c(6, 7, 8, 9, 10, 11, 12, 13),
    theta_indicator = c(1, 1, 0, 1, 1, 1, 1, 0)
  )
  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  slab_prior <- prior("normal", list(mean = 0, sd = 1))
  prior_list <- list(
    theta = prior_spike_and_slab(
      prior_parameter = slab_prior,
      prior_inclusion = prior("point", list(location = .5))
    )
  )

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "theta",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )

  expected_theta <- c(chain_1[, "theta"], chain_2[, "theta"])
  expected_indicator <- c(chain_1[, "theta_indicator"], chain_2[, "theta_indicator"])
  expected_theta[expected_indicator == 0] <- NA

  expect_equal(as.numeric(plot_data[, "theta"]), expected_theta)
  expect_equal(attr(plot_data, "chain"), rep(1:2, each = 8))
  expect_equal(attr(plot_data, "iter"), rep(1:8, times = 2))
  expect_equal(attr(plot_data, "parameter"), "theta")
  expect_false(inherits(attr(plot_data, "prior"), "prior.spike_and_slab"))
  expect_equal(attr(plot_data, "prior")$distribution, slab_prior$distribution)
  expect_equal(attr(plot_data, "prior")$parameters, slab_prior$parameters)
  expect_equal(attr(plot_data, "prior")$truncation, slab_prior$truncation)
})

test_that("diagnostic plot data rejects malformed transformations and insufficient slab mass", {
  fit <- .mock_diagnostics_fit()
  prior_list <- attr(fit, "prior_list")

  expect_error(
    .diagnostics_plot_data(
      fit = fit,
      parameter = "theta",
      prior_list = prior_list,
      transformations = list(theta = list(fun = "not-a-function")),
      transform_factors = FALSE
    ),
    "'transformations'"
  )

  chain_1 <- cbind(theta = 1:6, theta_indicator = c(1, 0, 0, 0, 0, 0))
  chain_2 <- cbind(theta = 7:12, theta_indicator = c(0, 1, 0, 0, 0, 0))
  sparse_fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(sparse_fit) <- c("BayesTools_fit", class(sparse_fit))
  sparse_prior_list <- list(
    theta = prior_spike_and_slab(prior("normal", list(mean = 0, sd = 1)))
  )

  expect_error(
    .diagnostics_plot_data(
      fit = sparse_fit,
      parameter = "theta",
      prior_list = sparse_prior_list,
      transformations = NULL,
      transform_factors = FALSE
    ),
    "enough samples under the slab"
  )
})
