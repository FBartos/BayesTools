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

.diagnostic_density_mass <- function(density) {
  sum(diff(density$x) * (head(density$y, -1) + tail(density$y, -1)) / 2)
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
  expect_true(isTRUE(attr(plot_data, "density_support_transformed")))

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

test_that("sampled diagnostic plots reject degenerate inputs clearly", {

  constant <- matrix(
    1,
    nrow = 4,
    ncol = 1,
    dimnames = list(NULL, "theta")
  )
  attr(constant, "chain") <- rep(1:2, each = 2)
  attr(constant, "prior") <- prior("normal", list(0, 1))

  expect_error(
    .diagnostics_plot_data_density(
      constant,
      n_points = 32,
      xlim = NULL
    ),
    "Density diagnostics.*not assessable.*constant"
  )
  expect_error(
    .diagnostics_plot_data_autocorrelation(
      constant,
      n_points = 10,
      lags = 1
    ),
    "Autocorrelation diagnostics.*not assessable.*constant"
  )

  too_short <- constant[1:2, , drop = FALSE]
  attr(too_short, "chain") <- 1:2
  attr(too_short, "prior") <- attr(constant, "prior")
  expect_error(
    .diagnostics_plot_data_density(
      too_short,
      n_points = 32,
      xlim = NULL
    ),
    "require at least two finite posterior samples"
  )

  empty <- matrix(numeric(), nrow = 0L, ncol = 0L)
  attr(empty, "chain") <- integer()
  expect_error(
    .diagnostics_plot_data_autocorrelation(
      empty,
      n_points = 10,
      lags = 1
    ),
    "require at least one parameter"
  )
})

test_that("diagnostic density plot data reflects bounded prior support", {
  chain_1 <- cbind(theta = seq(.005, .995, length.out = 400))
  chain_2 <- cbind(theta = rev(chain_1[, "theta"]))

  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  prior_list <- list(theta = prior("uniform", list(0, 1)))

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "theta",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )
  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 512, xlim = c(0, 1))

  for(chain_density in density_data$theta){
    expect_equal(range(chain_density$x), c(0, 1))
    expect_true(all(is.finite(chain_density$y)))
    expect_true(all(chain_density$y >= 0))
    expect_equal(.diagnostic_density_mass(chain_density), 1, tolerance = .03)
    expect_gt(chain_density$y[1], .75)
    expect_gt(chain_density$y[length(chain_density$y)], .75)
    expect_true(isTRUE(attr(chain_density, "boundary_reflection")))
  }
})

test_that("diagnostic density plot data reflects simplex coordinate support", {
  w1 <- seq(.005, .995, length.out = 400)
  chain_1 <- cbind("w[1]" = w1, "w[2]" = 1 - w1)
  chain_2 <- cbind("w[1]" = rev(w1), "w[2]" = rev(1 - w1))

  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  prior_list <- list(w = prior("dirichlet", list(alpha = c(1, 1))))

  expect_equal(
    BayesTools:::.diagnostics_prior_bounds(prior_list$w, "w"),
    list(lower = 0, upper = 1)
  )

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "w",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )
  expect_equal(colnames(plot_data), c("w[1]", "w[2]"))

  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 256, xlim = c(0, 1))

  for(parameter_density in density_data){
    for(chain_density in parameter_density){
      expect_equal(range(chain_density$x), c(0, 1))
      expect_true(all(is.finite(chain_density$y)))
      expect_true(all(chain_density$y >= 0))
      expect_true(isTRUE(attr(chain_density, "boundary_reflection")))
    }
  }
})

test_that("diagnostic density resolves heterogeneous composite-bias support by column", {
  selection <- prior_weightfunction(
    side = "one-sided",
    steps = .05,
    weights = wf_fixed(c(1, 1.5))
  )
  bias <- prior_bias(
    selection = selection,
    phacking = prior_phacking(report_scale = "alpha")
  )
  omega <- seq(1.05, 1.45, length.out = 400)
  alpha <- seq(.005, .995, length.out = 400)
  chain_1 <- cbind(
    "omega[1]" = 1,
    "omega[2]" = omega,
    alpha = alpha
  )
  chain_2 <- cbind(
    "omega[1]" = 1,
    "omega[2]" = rev(omega),
    alpha = rev(alpha)
  )

  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  prior_list <- list(bias = bias)

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "bias",
    prior_list = prior_list,
    transformations = NULL,
    transform_factors = FALSE
  )
  expect_equal(colnames(plot_data), c("omega[0.05,1]", "alpha"))

  density_data <- .diagnostics_plot_data_density(
    plot_data,
    n_points = 256,
    xlim = NULL
  )

  for(chain_density in density_data[["omega[0.05,1]"]]){
    expect_gt(max(chain_density$x), 1)
    expect_gt(max(chain_density$y), 0)
    expect_true(isTRUE(attr(chain_density, "boundary_reflection")))
  }
  for(chain_density in density_data$alpha){
    expect_true(all(chain_density$x >= 0 & chain_density$x <= 1))
    expect_gt(max(chain_density$y), 0)
    expect_true(isTRUE(attr(chain_density, "boundary_reflection")))
  }
})

test_that("diagnostic density does not reflect after custom transformations", {
  chain_1 <- cbind(theta = seq(.005, .995, length.out = 100))
  chain_2 <- cbind(theta = rev(chain_1[, "theta"]))

  fit <- coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  class(fit) <- c("BayesTools_fit", class(fit))
  prior_list <- list(theta = prior("uniform", list(0, 1)))

  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "theta",
    prior_list = prior_list,
    transformations = list(theta = list(fun = function(x, shift) x + shift, arg = list(shift = 2))),
    transform_factors = FALSE
  )
  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 128, xlim = c(2, 3))

  expect_false(isTRUE(attr(density_data$theta[[1]], "boundary_reflection")))
  expect_equal(range(density_data$theta[[1]]$x), c(2, 3))
  expect_true(max(density_data$theta[[1]]$y) > 0)
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

.sparse_acf_fit_for_test <- function(chains = 4L){

  values <- rep(0, 200L)
  values[c(1L, 30L, 90L)] <- c(1, 2, 4)
  indicator <- as.integer(seq_len(200L) %in% c(1L, 30L, 90L))
  fit <- list(mcmc = coda::mcmc.list(lapply(seq_len(chains), function(i){
    coda::mcmc(cbind(theta = values, theta_indicator = indicator))
  })), summary.pars = list(mutate = NULL))
  class(fit) <- c("BayesTools_fit", "runjags")
  attr(fit, "prior_list") <- list(theta = prior_spike_and_slab(
    prior("normal", list(0, 1)), prior_inclusion = prior("point", list(.5))))
  attach_test_parameter_map(fit)
}

test_that("R116 D10 sparse ACF preserves original lags and explicit pair diagnostics", {

  fit <- .sparse_acf_fit_for_test()
  input <- .diagnostics_plot_data(fit, "theta", attr(fit, "prior_list"), NULL, FALSE)
  warnings <- list()
  data <- withCallingHandlers(.diagnostics_plot_data_autocorrelation(input, 128L, 30L), warning = function(w){
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  values <- rep(NA_real_, 200L)
  values[c(1L, 30L, 90L)] <- c(1, 2, 4)
  reference <- as.numeric(stats::acf(values, lag.max = 30L, plot = FALSE, na.action = stats::na.pass)$acf)
  pairs <- vapply(0:30, function(lag){
    sum(is.finite(values[seq_len(200L - lag)]) & is.finite(values[lag + seq_len(200L - lag)]))
  }, integer(1))
  expect_length(warnings, 4L)
  expect_equal(pairs[c(1L, 30L)], c(3L, 1L))
  expect_identical(which(is.finite(reference)) - 1L, c(0L, 29L))
  for(i in seq_len(4L)){
    chain <- data$theta[[i]]
    expect_identical(names(chain), c("x", "y"))
    expect_identical(chain$x, 0:30)
    expect_identical(chain$y, reference)
    expect_identical(attr(chain, "pair_count"), pairs)
    expect_identical(attr(chain, "unavailability_reason"), ifelse(is.finite(reference), NA_character_, "no_finite_pairs"))
    expect_equal(attr(chain, "y_range"), c(0, 1))
    if(length(warnings) >= i){
      w <- warnings[[i]]
      unavailable <- (0:30)[!is.finite(reference)]
      expect_s3_class(w, "BayesTools_autocorrelation_unavailable")
      expect_s3_class(w, "BayesTools_plot_condition")
      expect_null(conditionCall(w))
      expect_identical(w$parameter, "theta")
      expect_identical(w$chain, i)
      expect_identical(w$unavailable_lags, unavailable)
      expect_identical(w$pair_count, pairs[!is.finite(reference)])
      expect_identical(w$reason, rep("no_finite_pairs", 29L))
      expect_identical(conditionMessage(w), paste0("Autocorrelation diagnostics for 'theta' in chain ", i,
        " are unavailable at lag(s) ", paste(unavailable, collapse = ", "),
        ": the retained chain is too short, has insufficient finite observation pairs, or produced a nonfinite autocorrelation. These lag values were set to NA; original iteration spacing was retained."))
    }
  }
  warnings <- list()
  plot <- withCallingHandlers(JAGS_diagnostics(fit, "theta", type = "autocorrelation", plot_type = "ggplot", lags = 30), warning = function(w){
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 4L)
  built <- ggplot2::ggplot_build(plot)
  for(layer in built$data){
    expect_identical(as.integer(layer$x), c(0L, 29L))
    expect_equal(layer$y, reference[c(1L, 30L)], tolerance = 0)
    expect_true(all(is.finite(layer$y)))
  }
  expect_error(.diagnostics_plot_data(.sparse_acf_fit_for_test(2L), "theta", attr(fit, "prior_list"), NULL, FALSE), "The parameter with a spike and slab prior did not result in enough samples under the slab for producing a diagnostic figure.", fixed = TRUE)
})

test_that("R116 D10 physical short chains pad requested lags without recomputation", {

  values <- c(-1, -.4, -.1, .3, .7, 1.1)
  input <- matrix(rep(values, 2L), ncol = 1L, dimnames = list(NULL, "theta"))
  attr(input, "chain") <- rep(1:2, each = 6L)
  attr(input, "parameter") <- "theta"
  warnings <- list()
  data <- withCallingHandlers(.diagnostics_plot_data_autocorrelation(input, 128L, 30L), warning = function(w){
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 2L)
  chain <- data$theta[[1L]]
  expect_identical(chain$x, 0:30)
  expect_equal(chain$y[1:6], as.numeric(stats::acf(values, lag.max = 30L, plot = FALSE)$acf), tolerance = 0)
  expect_true(all(is.na(chain$y[7:31])))
  expect_identical(attr(chain, "pair_count"), c(6:1, rep(0L, 25L)))
  expect_identical(attr(chain, "unavailability_reason"), c(rep(NA_character_, 6L), rep("chain_too_short", 25L)))
  expect_no_warning(.diagnostics_plot_data_autocorrelation(input, 128L, 3L))
  plot <- ggplot2::ggplot() + .geom_diagnostics.autocorrelation(chain)
  expect_identical(as.integer(ggplot2::ggplot_build(plot)$data[[1L]]$x), 0:5)
})
