skip_if_not_test_profile("unit")

.mock_public_diagnostics_fit <- function(n = 6, theta = NULL){

  if(is.null(theta)){
    theta <- seq(-1, 1, length.out = n)
  }
  fit <- list(
    mcmc = coda::mcmc.list(
      coda::mcmc(cbind(theta = theta)),
      coda::mcmc(cbind(theta = rev(theta)))
    ),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- c("BayesTools_fit", "runjags")
  attr(fit, "prior_list") <- list(
    theta = prior("normal", list(0, 1))
  )

  return(attach_test_parameter_map(fit))
}


test_that("the released JAGS_diagnostics_* wrappers select the diagnostic type", {

  fit <- .mock_public_diagnostics_fit()
  exports <- getNamespaceExports("BayesTools")
  expect_true(all(c(
    "JAGS_diagnostics",
    "JAGS_diagnostics_density",
    "JAGS_diagnostics_trace",
    "JAGS_diagnostics_autocorrelation"
  ) %in% exports))

  layer_data <- function(plot){
    lapply(ggplot2::ggplot_build(plot)$data, function(layer) layer[, intersect(colnames(layer), c("x", "y")), drop = FALSE])
  }

  for(type in c("density", "trace", "autocorrelation")){
    wrapper <- get(paste0("JAGS_diagnostics_", type), envir = asNamespace("BayesTools"))
    extra   <- if(type == "autocorrelation") list(lags = 3) else list()
    expect_equal(
      layer_data(do.call(wrapper, c(list(fit, parameter = "theta", plot_type = "ggplot"), extra))),
      layer_data(do.call(JAGS_diagnostics, c(list(fit, parameter = "theta", type = type, plot_type = "ggplot"), extra)))
    )
  }

  # the argument removed from the wrappers is rejected as in JAGS_diagnostics()
  expect_error(
    JAGS_diagnostics_density(fit, parameter = "theta", transform_orthonormal = TRUE),
    "The 'transform_orthonormal' argument was removed; use 'transform_factors = TRUE'.",
    fixed = TRUE
  )
})


test_that("JAGS autocorrelation diagnostics validate lags", {

  fit <- .mock_public_diagnostics_fit()
  invalid_lags <- list(
    c(1, 2),
    -1,
    1.5,
    NA_real_,
    NaN,
    Inf,
    -Inf
  )

  for(lags in invalid_lags){
    expect_error(
      JAGS_diagnostics(type = "autocorrelation", 
        fit,
        parameter = "theta",
        plot_type = "ggplot",
        lags = lags
      ),
      "'lags'"
    )
  }
})


test_that("JAGS autocorrelation diagnostics display negative correlations", {

  theta <- rep(c(-1, 1), 5)
  fit <- .mock_public_diagnostics_fit(theta = theta)
  expected_min <- min(stats::acf(theta, lag.max = 4, plot = FALSE)$acf)

  plot <- JAGS_diagnostics(type = "autocorrelation", 
    fit,
    parameter = "theta",
    plot_type = "ggplot",
    lags = 4
  )
  built_plot <- ggplot2::ggplot_build(plot)
  expect_lt(min(built_plot$data[[1]]$y), 0)
  expect_lte(built_plot$layout$panel_scales_y[[1]]$limits[[1]], expected_min)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  JAGS_diagnostics(type = "autocorrelation", 
    fit,
    parameter = "theta",
    plot_type = "base",
    lags = 4
  )
  expect_lte(graphics::par("usr")[[3]], expected_min)
})


test_that("JAGS autocorrelation diagnostics preserve valid scalar lag limits", {

  fit <- .mock_public_diagnostics_fit()

  plot <- JAGS_diagnostics(type = "autocorrelation", 
    fit,
    parameter = "theta",
    plot_type = "ggplot",
    lags = 3
  )
  plot_layers <- ggplot2::ggplot_build(plot)$data
  for(layer in plot_layers){
    expect_equal(layer$x, 0:3)
  }

  zero_lag_plot <- JAGS_diagnostics(type = "autocorrelation", 
    fit,
    parameter = "theta",
    plot_type = "ggplot",
    lags = 0
  )
  zero_lag_layers <- ggplot2::ggplot_build(zero_lag_plot)$data
  for(layer in zero_lag_layers){
    expect_equal(layer$x, 0)
  }
})


test_that("JAGS autocorrelation diagnostics retain requested unavailable lags", {

  fit <- .mock_public_diagnostics_fit(n = 6)
  warnings <- list()
  plot <- withCallingHandlers(JAGS_diagnostics(type = "autocorrelation",
    fit, parameter = "theta", plot_type = "ggplot", lags = 30), warning = function(w){
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 2L)
  expect_s3_class(warnings[[1L]], "BayesTools_autocorrelation_unavailable")
  expect_identical(warnings[[1L]]$unavailable_lags, 6:30)
  raw <- withCallingHandlers(.diagnostics_plot_data_autocorrelation(
    .diagnostics_plot_data(fit, "theta", attr(fit, "prior_list"), NULL, FALSE), 128L, 30L),
    warning = function(w) invokeRestart("muffleWarning"))
  expect_identical(raw$theta[[1L]]$x, 0:30)
  expect_true(all(is.na(raw$theta[[1L]]$y[7:31])))
  plot_layers <- ggplot2::ggplot_build(plot)$data

  expect_true(all(vapply(plot_layers, function(layer) nrow(layer) == 6L, logical(1))))
  for(layer in plot_layers){
    expect_equal(layer$x, 0:5)
  }
})


.indicator_trace_test_fit <- function(prior, states){

  chains <- lapply(states, function(state){
    coda::mcmc(cbind(theta = seq_along(state), theta_indicator = state),
               start = 101, thin = 2)
  })
  fit <- structure(list(mcmc = coda::mcmc.list(chains), summary.pars = list(mutate = NULL)),
                   class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(theta = prior)
  attach_test_parameter_map(fit)
}

test_that("mixture trace statistics preserve separate chains, batches and named states", {

  mixture <- prior_mixture(list(
    narrow = prior("normal", list(0, 1)), spike = prior("point", list(0)),
    wide = prior("normal", list(0, 5))
  ), is_null = c(FALSE, TRUE, FALSE))
  fit <- .indicator_trace_test_fit(mixture, list(
    c(1, 2, 3, 3, 2, 3, 2, 2, 1, 1, 1, 1),
    c(2, 2, 2, 2, 2, 1, 2, 3, 1, 2, 1, 2)
  ))

  inclusion <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot",
                                      statistic = "inclusion", batch_size = 4)
  layers <- ggplot2::ggplot_build(inclusion)$data
  expect_length(layers, 2L)
  expect_equal(layers[[1]]$x, c(104, 112, 120))
  expect_equal(layers[[2]]$x, layers[[1]]$x)
  expect_equal(layers[[1]]$y, c(.75, .25, 1))
  expect_equal(layers[[2]]$y, c(0, .5, .5))

  switching <- JAGS_diagnostics(fit, "theta", "trace", plot_type = "ggplot",
                                statistic = "switching", batch_size = 4)
  layers <- ggplot2::ggplot_build(switching)$data
  expect_equal(layers[[1]]$y, c(2/3, 2/3, 0))
  expect_equal(layers[[2]]$y, c(0, 1, 1))

  plots <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot",
                                  statistic = "component", batch_size = 4)
  expect_named(plots, names(mixture))
  expected <- list(c(.25, 0, 1), c(.25, .75, 0), c(.5, .25, 0))
  for (i in seq_along(plots)) {
    expect_equal(ggplot2::ggplot_build(plots[[i]])$data[[1]]$y, expected[[i]])
  }
})

test_that("binary traces use the active state and preserve the final partial batch", {

  state <- c(1, 0, 0, 0, 1, 1, 1, 0, 1)
  fit <- .indicator_trace_test_fit(prior_spike_and_slab(prior("normal", list(0, 1))),
                                   list(state))
  for (statistic in c("inclusion", "component")) {
    plot <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot", statistic = statistic)
    layer <- ggplot2::ggplot_build(plot)$data[[1]]
    expect_equal(layer$x, c(103, 109, 115))
    expect_equal(layer$y, c(1/3, 2/3, 2/3))
  }
  plot <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot",
                                 statistic = "inclusion", batch_size = 4)
  layer <- ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(layer$x, c(104, 112, 117))
  expect_equal(layer$y, c(.25, .75, 1))

  switching <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot",
                                      statistic = "switching", batch_size = 4)
  expect_equal(ggplot2::ggplot_build(switching)$data[[1]]$y, c(1/3, 1/3, NA_real_))
  single <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot",
                                   statistic = "switching", batch_size = 9)
  expect_s3_class(single$layers[[1]]$geom, "GeomPoint")
  expect_equal(ggplot2::ggplot_build(single)$data[[1]]$y, .5)
})

test_that("indicator traces require sampled states and meaningful group metadata", {

  fit <- .mock_public_diagnostics_fit()
  expect_error(JAGS_diagnostics_trace(fit, "theta", statistic = "inclusion"),
               class = "BayesTools_plot_condition")
  expect_error(JAGS_diagnostics(fit, "theta", "density", statistic = "component"))

  mixture <- prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(0, 5))))
  fit <- .indicator_trace_test_fit(mixture, list(c(2, 2, 1, 1)))
  expect_error(JAGS_diagnostics_trace(fit, "theta", statistic = "inclusion"),
               class = "BayesTools_plot_condition")
  plots <- JAGS_diagnostics_trace(fit, "theta", plot_type = "ggplot", statistic = "component")
  expect_length(plots, 2L)
  expect_error(JAGS_diagnostics_trace(fit, "theta", statistic = "switching", batch_size = 1))
  bad_fit <- .indicator_trace_test_fit(mixture, list(c(2, 2, 3, 1)))
  expect_error(JAGS_diagnostics_trace(bad_fit, "theta", statistic = "component"),
               class = "BayesTools_plot_condition")
})

test_that("each weight-function prior has its own component trace even with no visits", {

  branch <- function(cut) prior_weightfunction("one-sided", cut, wf_cumulative(c(1, 1)))
  priors <- prior_mixture(list(prior_none(), branch(.025), branch(.05)),
                          is_null = c(TRUE, FALSE, FALSE))
  draws <- cbind(mu = 1:6, bias_indicator = c(1, 2, 2, 1, 1, 1),
                 "omega[1]" = 1, "omega[2]" = .8, "omega[3]" = .6)
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)),
                        summary.pars = list(mutate = NULL)),
                   class = c("BayesTools_fit", "runjags"))
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)), bias = priors)
  fit <- attach_test_parameter_map(fit)

  plots <- JAGS_diagnostics_trace(fit, "bias", plot_type = "ggplot",
                                  statistic = "component", batch_size = 3)
  expect_length(plots, 3L)
  expect_length(unique(names(plots)), 3L)
  expect_equal(ggplot2::ggplot_build(plots[[2]])$data[[1]]$y, c(2/3, 0))
  expect_equal(ggplot2::ggplot_build(plots[[3]])$data[[1]]$y, c(0, 0))
})
