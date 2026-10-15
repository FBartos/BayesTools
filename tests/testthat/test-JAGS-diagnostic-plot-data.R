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
    "Density diagnostics for 'theta' in chain 1 are unavailable because the posterior samples are constant. Use type = \"trace\" to inspect the chain.",
    fixed = TRUE,
    class = "BayesTools_plot_condition"
  )
  expect_error(
    .diagnostics_plot_data_autocorrelation(
      constant,
      n_points = 10,
      lags = 1
    ),
    "Autocorrelation diagnostics for 'theta' in chain 1 are unavailable because the posterior samples are constant. Use type = \"trace\" to inspect the chain.",
    fixed = TRUE,
    class = "BayesTools_plot_condition"
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
    expect_equal(layer$xmax - layer$xmin, rep(.9, 2L), tolerance = 1e-12)
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

# The chains of a diagnostics test fit: the columns of `draws` for the first
# chain and the same draws in reverse order for the second one.
.diagnostics_two_chain_fit_for_test <- function(draws){

  fit <- coda::mcmc.list(
    coda::mcmc(draws),
    coda::mcmc(draws[rev(seq_len(nrow(draws))), , drop = FALSE])
  )
  class(fit) <- c("BayesTools_fit", class(fit))
  fit
}

# The mass of a diagnostic density between the bounds, and its largest ordinate
# outside of them (zero for a density that respects the support).
.diagnostic_density_support_check <- function(density, bounds){

  inside <- density$x >= bounds[1] & density$x <= bounds[2]
  list(
    mass    = .diagnostic_density_mass(list(x = density$x[inside], y = density$y[inside])),
    outside = if(any(!inside)) max(density$y[!inside]) else 0
  )
}

test_that("diagnostics of the components of a single bias prior use the bounds of the prior", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  phacking  <- prior_phacking(form = "linear")
  single    <- prior_bias(selection, phacking)
  mixture   <- prior_mixture(list(prior_none(), single))

  pi_null_bounds <- .diagnostics_prior_bounds(mixture, "pi_null")
  expect_equal(pi_null_bounds, list(lower = 0, upper = 0.0836822), tolerance = 1e-6)

  grid  <- seq(.005, .995, length.out = 300)
  fit   <- .diagnostics_two_chain_fit_for_test(cbind(
    "omega[1]" = 1,
    "omega[2]" = grid,
    alpha      = grid,
    pi_null    = pi_null_bounds$upper * grid
  ))
  xlim  <- list(
    omega   = c(-.25, 1.25),
    alpha   = c(-.25, 1.25),
    pi_null = pi_null_bounds$upper * c(-.25, 1.25)
  )
  names <- c(omega = "omega[0.025,1]", alpha = "alpha", pi_null = "pi_null")
  upper <- c(omega = 1, alpha = 1, pi_null = pi_null_bounds$upper)

  for(parameter in c("omega", "alpha", "pi_null")){

    single_data  <- .diagnostics_plot_data(fit, parameter, list(bias = single),  NULL, FALSE)
    mixture_data <- .diagnostics_plot_data(fit, parameter, list(bias = mixture), NULL, FALSE)

    # the component keeps its own column, named by its parameter
    expect_identical(colnames(single_data), names[[parameter]])
    expect_identical(as.numeric(single_data), as.numeric(mixture_data))
    expect_identical(attr(single_data, "prior"), single)

    # the same bounds as the component inside a bias mixture
    bounds_single  <- .diagnostics_density_bounds(attr(single_data, "prior"), parameter, colnames(single_data))
    bounds_mixture <- .diagnostics_density_bounds(attr(mixture_data, "prior"), parameter, colnames(mixture_data))
    expect_equal(bounds_single, bounds_mixture)
    expect_equal(unlist(bounds_single), c(lower = 0, upper = upper[[parameter]]))

    # the density is the one of the mixture and stays on the support
    density_single  <- .diagnostics_plot_data_density(single_data,  n_points = 256, xlim = xlim[[parameter]])
    density_mixture <- .diagnostics_plot_data_density(mixture_data, n_points = 256, xlim = xlim[[parameter]])
    expect_equal(density_single, density_mixture)
    expect_identical(names(density_single), names[[parameter]])
    for(chain_density in density_single[[1]]){
      check <- .diagnostic_density_support_check(chain_density, c(0, upper[[parameter]]))
      expect_true(isTRUE(attr(chain_density, "boundary_reflection")))
      expect_equal(check$mass, 1, tolerance = .03)
      expect_identical(check$outside, 0)
    }
  }

  # the composite parameter reports the omega bin and the null mass
  composite_data <- .diagnostics_plot_data(fit, "bias", list(bias = single), NULL, FALSE)
  expect_identical(colnames(composite_data), c("omega[0.025,1]", "pi_null"))
})

test_that("diagnostics of a single PET, PEESE, weightfunction or p-hacking bias prior use its bounds", {

  # `component` names the parameter of the same column in the bias mixture
  expect_bounds <- function(fit, parameter, component, bias, bounds, xlim){

    plot_data <- .diagnostics_plot_data(fit, parameter, list(bias = bias), NULL, FALSE)
    mixture   <- .diagnostics_plot_data(fit, component, list(bias = prior_mixture(list(prior_none(), bias))), NULL, FALSE)
    expect_identical(attr(plot_data, "prior"), bias)
    expect_identical(as.numeric(plot_data), as.numeric(mixture))
    expect_equal(
      unlist(.diagnostics_density_bounds(attr(plot_data, "prior"), attr(plot_data, "parameter"), colnames(plot_data)[1])),
      c(lower = bounds[1], upper = bounds[2])
    )
    density <- .diagnostics_plot_data_density(plot_data, n_points = 256, xlim = xlim)
    for(chain_density in density[[1]]){
      check <- .diagnostic_density_support_check(chain_density, bounds)
      expect_true(isTRUE(attr(chain_density, "boundary_reflection")))
      expect_identical(check$outside, 0)
      expect_gt(check$mass, .9)
    }
    invisible(plot_data)
  }

  grid <- seq(.005, 2.995, length.out = 300)

  # PET and PEESE: the truncation of the prior; 'bias' names the PET column
  for(bias in list(prior_PET("normal", list(0, 1), truncation = list(0, Inf)),
                   prior_PEESE("normal", list(0, 1), truncation = list(0, Inf)))){
    column <- if(is.prior.PET(bias)) "PET" else "PEESE"
    fit    <- .diagnostics_two_chain_fit_for_test(matrix(grid, ncol = 1L, dimnames = list(NULL, column)))
    for(parameter in c(column, "bias")){
      plot_data <- expect_bounds(fit, parameter, column, bias, c(0, Inf), c(-1, 4))
      expect_identical(colnames(plot_data), column)
    }
  }

  # a standalone weightfunction: its weights
  weightfunction <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  fit <- .diagnostics_two_chain_fit_for_test(cbind("omega[1]" = 1, "omega[2]" = seq(.005, .995, length.out = 300)))
  for(parameter in c("omega", "bias")){
    plot_data <- expect_bounds(fit, parameter, "omega", weightfunction, c(0, 1), c(-.25, 1.25))
    expect_identical(colnames(plot_data), "omega[0.025,1]")
  }

  # a standalone p-hacking prior: the share of rejections and the null mass
  phacking <- prior_phacking(form = "linear")
  upper    <- .diagnostics_prior_bounds(phacking, "pi_null")$upper
  fit <- .diagnostics_two_chain_fit_for_test(cbind(
    alpha   = seq(.005, .995, length.out = 300),
    pi_null = upper * seq(.005, .995, length.out = 300)
  ))
  expect_bounds(fit, "alpha",   "alpha",   phacking, c(0, 1),     c(-.25, 1.25))
  expect_bounds(fit, "pi_null", "pi_null", phacking, c(0, upper), upper * c(-.25, 1.25))
})

test_that("a bias mixture is drawn as the components of its branches (F221)", {

  bias <- prior_mixture(list(
    prior_none(),
    prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1))),
    prior_PET("normal", list(0, 1), list(0, Inf))
  ), is_null = c(TRUE, FALSE, FALSE))
  grid <- seq(.01, .99, length.out = 120)
  draws <- cbind(mu = seq(-1, 1, length.out = 120), bias_indicator = rep(1:3, 40L),
                 "omega[1]" = 1, "omega[2]" = grid, PET = grid)
  fit <- structure(list(
    mcmc = coda::mcmc.list(coda::mcmc(draws[1:60, ]), coda::mcmc(draws[61:120, ])),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)), bias = bias)
  fit <- attach_test_parameter_map(fit)

  for(type in c("density", "trace", "autocorrelation")){
    plots <- JAGS_diagnostics(fit, "bias", type, plot_type = "ggplot", n_points = 32)
    expect_length(plots, 2L)
    for(plot in plots){
      expect_s3_class(plot, "ggplot")
    }
  }
})

test_that("the trace axis shows the retained MCMC iteration of thinned chains (F223)", {

  chain_1 <- coda::mcmc(cbind(theta = seq(-1, 1, length.out = 6)), start = 11, thin = 2)
  chain_2 <- coda::mcmc(cbind(theta = seq(1, -1, length.out = 6)), start = 11, thin = 2)
  fit <- coda::mcmc.list(chain_1, chain_2)
  class(fit) <- c("BayesTools_fit", class(fit))
  plot_data <- .diagnostics_plot_data(
    fit = fit,
    parameter = "theta",
    prior_list = list(theta = prior("normal", list(0, 1))),
    transformations = NULL,
    transform_factors = FALSE
  )
  expect_identical(attr(plot_data, "iter"), rep(seq(11, by = 2, length.out = 6), 2L))

  trace_data <- .diagnostics_plot_data_trace(plot_data, n_points = 10, ylim = NULL)
  expect_identical(trace_data$theta[[1L]]$x, seq(11, by = 2, length.out = 6))
  expect_identical(attr(trace_data$theta, "x_range"), c(11, 21))
})

test_that("PET and PEESE diagnostics of a bias prior without that branch are refused clearly (F216)", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  bias      <- prior_bias(selection, prior_phacking(form = "linear"))
  draws <- cbind(mu = seq(-1, 1, length.out = 60), "omega[1]" = 1,
                 "omega[2]" = seq(.01, .99, length.out = 60))
  fit <- structure(list(
    mcmc = coda::mcmc.list(coda::mcmc(draws)),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)), bias = bias)
  fit <- attach_test_parameter_map(fit)

  for(parameter in c("PET", "PEESE")){
    expect_error(
      JAGS_diagnostics(fit, parameter, "density", plot_type = "ggplot"),
      paste0("'", parameter, "' diagnostics require at least one ", parameter, " component in the 'bias' prior."),
      fixed = TRUE
    )
  }
})

test_that("JAGS_diagnostics() draws the components of a single bias prior", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  bias      <- prior_bias(selection, prior_phacking(form = "linear"))
  upper     <- .diagnostics_prior_bounds(bias, "pi_null")$upper

  grid  <- seq(.005, .995, length.out = 120)
  draws <- cbind(
    mu         = seq(-1, 1, length.out = 120),
    "omega[1]" = 1,
    "omega[2]" = grid,
    alpha      = rev(grid),
    phack_kind = 1,
    pi_null    = upper * grid
  )
  fit <- structure(list(
    mcmc = coda::mcmc.list(coda::mcmc(draws)),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)), bias = bias)
  fit <- attach_test_parameter_map(fit)

  for(parameter in c("omega", "alpha", "pi_null", "bias")){
    for(type in c("density", "trace", "autocorrelation")){
      plots <- JAGS_diagnostics(fit, parameter, type, plot_type = "ggplot", n_points = 64)
      if(inherits(plots, "ggplot")) plots <- list(plots)
      expect_length(plots, if(parameter == "bias") 2L else 1L)
      for(plot in plots){
        expect_s3_class(plot, "ggplot")
        expect_true(all(is.finite(ggplot2::ggplot_build(plot)$data[[1]]$y)))
      }
    }
  }

  # the density of alpha lives on [0, 1]
  alpha_density <- ggplot2::ggplot_build(
    JAGS_diagnostics(fit, "alpha", "density", plot_type = "ggplot", xlim = c(-.25, 1.25), n_points = 256)
  )$data[[1]]
  outside <- alpha_density$x < 0 | alpha_density$x > 1
  expect_true(any(outside))
  expect_identical(alpha_density$y[outside], rep(0, sum(outside)))
  expect_gt(max(alpha_density$y), .5)
})

test_that("diagnostic plot data masks the spike rows of every coefficient of a factor spike-and-slab", {

  info <- JAGS_formula(
    ~ f, "mu", data.frame(f = factor(c("a", "b", "c"))),
    list(
      intercept = prior("point", list(0)),
      f = prior_spike_and_slab(prior_factor("normal", list(0, 1), contrast = "treatment"))
    )
  )
  indicator <- rep(c(0, 1, 1), 6)
  draws <- cbind(
    mu_f_indicator = indicator,
    "mu_f[1]"      = indicator * seq_along(indicator),
    "mu_f[2]"      = indicator * (seq_along(indicator) + 1)
  )
  fit <- structure(list(
    mcmc = coda::mcmc.list(coda::mcmc(draws)),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list")     <- info$prior_list
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  fit <- attach_test_parameter_map(fit)

  plot_data <- .diagnostics_plot_data(fit, "mu_f", info$prior_list, NULL, FALSE)

  # six of the 18 draws are excluded: the same rows are NA in each coefficient
  excluded <- which(indicator == 0)
  expect_identical(colnames(plot_data), c("mu_f[b]", "mu_f[c]"))
  expect_identical(sum(is.na(plot_data)), 12L)
  expect_identical(which(is.na(plot_data[, "mu_f[b]"])), excluded)
  expect_identical(which(is.na(plot_data[, "mu_f[c]"])), excluded)
  expect_identical(as.numeric(plot_data[-excluded, "mu_f[b]"]), draws[-excluded, "mu_f[1]"])
  expect_identical(as.numeric(plot_data[-excluded, "mu_f[c]"]), draws[-excluded, "mu_f[2]"])
  expect_equal(attr(plot_data, "iter"), seq_len(18L))

  # the diagnostics describe the draws under the slab only
  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = NULL)
  for(coefficient in c("mu_f[b]", "mu_f[c]")){
    included <- as.numeric(plot_data[-excluded, coefficient])
    expect_equal(attr(density_data[[coefficient]], "x_range"), range(included))
    expect_equal(density_data[[coefficient]][[1]]$bw, stats::bw.nrd0(included))
  }
  trace_data <- .diagnostics_plot_data_trace(plot_data, n_points = 64, ylim = NULL)
  expect_identical(which(is.na(trace_data[["mu_f[b]"]][[1]]$y)), excluded)
})

# A fit of the three-level factor `f` of the formula parameter `mu` with the
# prior `f_prior`; `coefficients` holds the draws of the two coefficient
# columns.
.diagnostics_factor_fit_for_test <- function(f_prior, coefficients){

  info <- JAGS_formula(
    ~ f, "mu", data.frame(f = factor(c("a", "b", "c"))),
    list(intercept = prior("point", list(0)), f = f_prior)
  )
  draws <- cbind(
    mu_f_indicator = rep(1:2, length.out = nrow(coefficients)),
    "mu_f[1]"      = coefficients[, 1],
    "mu_f[2]"      = coefficients[, 2]
  )
  fit <- structure(list(
    mcmc = coda::mcmc.list(coda::mcmc(draws)),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list")     <- info$prior_list
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  attach_test_parameter_map(fit)
}

test_that("density diagnostics of a factor mixture use the support of the whole mixture", {

  # the atom of the first component and the draws of the second one
  positive <- cbind(rep(c(0, .1, 0, .3, 0, .5, 0, .7, 0, .9), 2), rep(c(0, .2, 0, .4, 0, .6, 0, .8, 0, 1), 2))
  mixed    <- cbind(rep(c(0, .1, 0, -.3, 0, .5, 0, -.7, 0, .9), 2), rep(c(0, -.2, 0, .4, 0, -.6, 0, .8, 0, -1), 2))
  maxima   <- function(density) vapply(density, function(coefficient) max(coefficient[[1]]$y), numeric(1))

  # a point at 0 first and a normal second: the support is the one of the normal
  mixture <- prior_mixture(list(
    prior("point", list(0)),
    prior_factor("normal", list(0, 1), contrast = "treatment")
  ))
  fit       <- .diagnostics_factor_fit_for_test(mixture, positive)
  prior_list <- attr(fit, "prior_list")
  plot_data <- .diagnostics_plot_data(fit, "mu_f", prior_list, NULL, FALSE)
  expect_identical(colnames(plot_data), c("mu_f[b]", "mu_f[c]"))
  expect_identical(attr(plot_data, "prior"), prior_list$mu_f)
  expect_equal(.diagnostics_prior_bounds(attr(plot_data, "prior"), "mu_f"), list(lower = -Inf, upper = Inf))
  density <- .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = c(.05, 1.1))
  expect_true(all(maxima(density) > 0))
  for(coefficient in density){
    expect_false(isTRUE(attr(coefficient[[1]], "boundary_reflection")))
  }

  # through the public function, as a ggplot per coefficient
  plots <- JAGS_diagnostics(fit, "mu_f", "density", plot_type = "ggplot", xlim = c(.05, 1.1), n_points = 64)
  expect_length(plots, 2L)
  for(plot in plots){
    expect_gt(max(ggplot2::ggplot_build(plot)$data[[1]]$y), 0)
  }

  # a point outside of the truncation of the second component is no support
  mixture <- prior_mixture(list(
    prior("point", list(5)),
    prior_factor("normal", list(.5, 1), list(0, 1), contrast = "treatment")
  ))
  expect_equal(.diagnostics_prior_bounds(mixture, "mu_f"), list(lower = 0, upper = 1))

  # two continuous components: the union of their supports
  mixture <- prior_mixture(list(
    prior_factor("normal", list(0, 1), list(0, Inf), contrast = "treatment"),
    prior_factor("normal", list(0, 1), list(-Inf, 0), contrast = "treatment")
  ))
  expect_equal(.diagnostics_prior_bounds(mixture, "mu_f"), list(lower = -Inf, upper = Inf))
  fit        <- .diagnostics_factor_fit_for_test(mixture, mixed)
  plot_data  <- .diagnostics_plot_data(fit, "mu_f", attr(fit, "prior_list"), NULL, FALSE)
  density    <- .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = c(-1.2, 1.2))
  expect_true(all(maxima(density) > 0))
  for(coefficient in density){
    negative <- coefficient[[1]]$x < 0
    expect_true(any(negative))
    expect_gt(max(coefficient[[1]]$y[negative]), 0)
  }

  # a point and an orthonormal (or meandif) component: the real line
  for(contrast in c("orthonormal", "meandif")){
    mixture <- prior_mixture(list(
      prior("point", list(0)),
      prior_factor("mnormal", list(0, 1), contrast = contrast)
    ))
    fit       <- .diagnostics_factor_fit_for_test(mixture, positive)
    plot_data <- .diagnostics_plot_data(fit, "mu_f", attr(fit, "prior_list"), NULL, FALSE)
    expect_equal(.diagnostics_prior_bounds(attr(plot_data, "prior"), "mu_f"), list(lower = -Inf, upper = Inf))
    density   <- .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = c(.05, 1.1))
    expect_true(all(maxima(density) > 0))
  }

  # a mixture of atoms keeps the range of the atoms
  expect_equal(
    .diagnostics_prior_bounds(prior_mixture(list(prior("point", list(0)), prior("point", list(2)))), "theta"),
    list(lower = 0, upper = 2)
  )
  # a null component of a publication bias mixture is an atom as well
  expect_equal(
    .diagnostics_prior_bounds(prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1), truncation = list(1, Inf)))), "PET"),
    list(lower = 1, upper = Inf)
  )

  # a single factor prior is its own support
  single <- prior_factor("normal", list(0, 1), list(0, Inf), contrast = "treatment")
  fit        <- .diagnostics_factor_fit_for_test(single, positive)
  prior_list <- attr(fit, "prior_list")
  plot_data  <- .diagnostics_plot_data(fit, "mu_f", prior_list, NULL, FALSE)
  expect_identical(attr(plot_data, "prior"), prior_list$mu_f)
  expect_equal(.diagnostics_prior_bounds(attr(plot_data, "prior"), "mu_f"), list(lower = 0, upper = Inf))
  density <- .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = c(0, 1.1))
  for(coefficient in density){
    expect_true(isTRUE(attr(coefficient[[1]], "boundary_reflection")))
    expect_gt(max(coefficient[[1]]$y), 0)
  }
})

# A fit of the bias prior `bias` with the monitored columns of its omega bins,
# alpha, kind and null mass (and `mu`, as in a model fit).
.diagnostics_bias_fit_for_test <- function(bias, draws){

  fit <- structure(list(
    mcmc = coda::mcmc.list(coda::mcmc(draws), coda::mcmc(draws[rev(seq_len(nrow(draws))), , drop = FALSE])),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)), bias = bias)
  attach_test_parameter_map(fit)
}

test_that("density and autocorrelation diagnostics skip the weight bins a prior fixes", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  bias      <- prior_bias(selection, prior_phacking(form = "linear"))
  upper     <- .diagnostics_prior_bounds(bias, "pi_null")$upper
  grid      <- seq(.005, .995, length.out = 120)
  draws     <- cbind(
    mu         = seq(-1, 1, length.out = 120),
    "omega[1]" = 1,
    "omega[2]" = .5,
    alpha      = grid,
    phack_kind = 1,
    pi_null    = upper * grid
  )
  fit <- .diagnostics_bias_fit_for_test(bias, draws)

  # the composite parameter has a declared constant (omega[0.025,1]) and the null mass
  plot_data <- .diagnostics_plot_data(fit, "bias", attr(fit, "prior_list"), NULL, FALSE)
  expect_identical(colnames(plot_data), c("omega[0.025,1]", "pi_null"))
  expect_identical(attr(plot_data, "structural"), c(TRUE, FALSE))

  # density and autocorrelation return the null mass only; the trace both
  density_data <- .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = NULL)
  expect_identical(names(density_data), "pi_null")
  expect_identical(attr(density_data, "parameter_name"), "pi_null")
  expect_length(attr(density_data, "label_parts"), 1L)
  expect_identical(.bt_label(attr(density_data, "label_parts"), style = "selector"), "pi_null")

  autocorrelation_data <- .diagnostics_plot_data_autocorrelation(plot_data, n_points = 64, lags = 5)
  expect_identical(names(autocorrelation_data), "pi_null")
  expect_identical(autocorrelation_data$pi_null[[1]]$x, 0:5)

  trace_data <- .diagnostics_plot_data_trace(plot_data, n_points = 64, ylim = NULL)
  expect_identical(names(trace_data), c("omega[0.025,1]", "pi_null"))

  # the same through the public function: one plot of the null mass, two traces
  for(type in c("density", "autocorrelation")){
    plot <- JAGS_diagnostics(fit, "bias", type, plot_type = "ggplot", n_points = 64)
    expect_s3_class(plot, "ggplot")
    expect_true(all(is.finite(ggplot2::ggplot_build(plot)$data[[1]]$y)))
  }
  plots <- JAGS_diagnostics(fit, "bias", "trace", plot_type = "ggplot")
  expect_length(plots, 2L)

  # the omega bin on its own has nothing to show, as for a spike prior
  for(type in c("density", "autocorrelation")){
    expect_message(
      result <- JAGS_diagnostics(fit, "omega", type, plot_type = "ggplot"),
      "No diagnostic plots are produced for a spike prior distribution",
      fixed = TRUE
    )
    expect_null(result)
  }
  expect_s3_class(JAGS_diagnostics(fit, "omega", "trace", plot_type = "ggplot"), "ggplot")
  expect_s3_class(JAGS_diagnostics(fit, "alpha", "density", plot_type = "ggplot"), "ggplot")
})

test_that("diagnostics refuse constant draws that no prior declares constant", {

  # the bins of cumulative weights are sampled: constant draws are no structure
  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  bias      <- prior_bias(selection, prior_phacking(form = "linear"))
  upper     <- .diagnostics_prior_bounds(bias, "pi_null")$upper
  grid      <- seq(.005, .995, length.out = 120)
  draws     <- cbind(
    mu         = seq(-1, 1, length.out = 120),
    "omega[1]" = 1,
    "omega[2]" = .5,
    alpha      = grid,
    phack_kind = 1,
    pi_null    = upper * grid
  )
  fit       <- .diagnostics_bias_fit_for_test(bias, draws)
  plot_data <- .diagnostics_plot_data(fit, "bias", attr(fit, "prior_list"), NULL, FALSE)
  expect_identical(attr(plot_data, "structural"), c(FALSE, FALSE))
  expect_error(
    .diagnostics_plot_data_density(plot_data, n_points = 64, xlim = NULL),
    "Density diagnostics for 'omega\\[0.025,1\\]'.*unavailable.*constant"
  )
  expect_error(
    .diagnostics_plot_data_autocorrelation(plot_data, n_points = 64, lags = 5),
    "Autocorrelation diagnostics for 'omega\\[0.025,1\\]'.*unavailable.*constant"
  )

  # a bin the prior fixes is diagnosed when its draws vary (a corrupted sample)
  fixed_bias <- prior_bias(
    prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5))),
    prior_phacking(form = "linear")
  )
  draws[, "omega[2]"] <- seq(.1, .9, length.out = 120)
  fit       <- .diagnostics_bias_fit_for_test(fixed_bias, draws)
  plot_data <- .diagnostics_plot_data(fit, "bias", attr(fit, "prior_list"), NULL, FALSE)
  expect_identical(attr(plot_data, "structural"), c(FALSE, FALSE))
  expect_identical(
    names(.diagnostics_plot_data_density(plot_data, n_points = 64, xlim = NULL)),
    c("omega[0.025,1]", "pi_null")
  )
})

test_that("a weightfunction with fixed weights is diagnosed like a spike prior", {

  fit <- structure(list(
    mcmc = coda::mcmc.list(
      coda::mcmc(cbind("omega[1]" = rep(1, 40), "omega[2]" = rep(.5, 40))),
      coda::mcmc(cbind("omega[1]" = rep(1, 40), "omega[2]" = rep(.5, 40)))
    ),
    summary.pars = list(mutate = NULL)
  ), class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(omega = prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5))))
  fit <- attach_test_parameter_map(fit)

  for(type in c("density", "autocorrelation")){
    expect_message(
      result <- JAGS_diagnostics(fit, "omega", type, plot_type = "ggplot"),
      "No diagnostic plots are produced for a spike prior distribution",
      fixed = TRUE
    )
    expect_null(result)
  }
  expect_s3_class(JAGS_diagnostics(fit, "omega", "trace", plot_type = "ggplot"), "ggplot")
})

test_that("a diagnostic transformation applies to the PET and PEESE column and the omega bins it plots", {

  inc <- function(x) x + 10
  raw_values <- function(fit, parameter, prior_list){
    as.numeric(.diagnostics_plot_data(fit, parameter, prior_list, NULL, FALSE))
  }

  # a bias prior plotted as PET or PEESE: the column, as under its column name
  grid <- seq(.005, 2.995, length.out = 120)
  for(bias in list(prior_PET("normal", list(0, 1), truncation = list(0, Inf)),
                   prior_PEESE("normal", list(0, 1), truncation = list(0, Inf)))){
    column <- if(is.prior.PET(bias)) "PET" else "PEESE"
    fit    <- .diagnostics_two_chain_fit_for_test(matrix(grid, ncol = 1L, dimnames = list(NULL, column)))
    raw    <- raw_values(fit, "bias", list(bias = bias))
    for(parameter in c(column, "bias")){
      plot_data <- .diagnostics_plot_data(fit, parameter, list(bias = bias),
        stats::setNames(list(list(fun = inc)), parameter), FALSE)
      expect_identical(as.numeric(plot_data), raw + 10)
      expect_true(isTRUE(attr(plot_data, "density_support_transformed")))
    }
  }

  # a weightfunction: every omega bin it plots
  weightfunction <- prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 1, 1)))
  draws <- cbind("omega[1]" = 1, "omega[2]" = seq(.2, .7, length.out = 100), "omega[3]" = seq(.8, .3, length.out = 100))
  fit   <- .diagnostics_two_chain_fit_for_test(draws)
  raw   <- .diagnostics_plot_data(fit, "omega", list(omega = weightfunction), NULL, FALSE)
  expect_identical(ncol(raw), 2L)
  for(case in list(list("omega", list(omega = weightfunction)), list("bias", list(bias = weightfunction)))){
    plot_data <- .diagnostics_plot_data(fit, case[[1L]], case[[2L]],
      stats::setNames(list(list(fun = inc)), case[[1L]]), FALSE)
    expect_identical(as.numeric(plot_data), as.numeric(raw) + 10)
  }

  # the omega bins of a bias prior with a selection component; a requested
  # name that is a column keeps its column
  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  bias      <- prior_bias(selection, prior_phacking(form = "linear"))
  upper     <- .diagnostics_prior_bounds(bias, "pi_null")$upper
  g         <- seq(.005, .995, length.out = 100)
  fit       <- .diagnostics_two_chain_fit_for_test(cbind("omega[1]" = 1, "omega[2]" = g, alpha = g, pi_null = upper * g))
  for(parameter in c("omega", "alpha")){
    transformed <- .diagnostics_plot_data(fit, parameter, list(bias = bias),
      stats::setNames(list(list(fun = inc)), parameter), FALSE)
    expect_identical(as.numeric(transformed), raw_values(fit, parameter, list(bias = bias)) + 10)
  }

  # through JAGS_diagnostics(): the plotted values follow the transformation
  draws <- cbind(mu = seq(-1, 1, length.out = 120), PET = seq(.005, 2.995, length.out = 120))
  bias  <- prior_PET("normal", list(0, 1), truncation = list(0, Inf))
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)), summary.pars = list(mutate = NULL)),
    class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)), bias = bias)
  fit <- attach_test_parameter_map(fit)
  plain   <- JAGS_diagnostics(fit, "bias", "trace", plot_type = "ggplot")
  shifted <- JAGS_diagnostics(fit, "bias", "trace", transformations = list(bias = list(fun = inc)), plot_type = "ggplot")
  expect_equal(ggplot2::ggplot_build(shifted)$data[[1]]$y, ggplot2::ggplot_build(plain)$data[[1]]$y + 10)
})
