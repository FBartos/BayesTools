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

  return(fit)
}


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
      JAGS_diagnostics_autocorrelation(
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

  plot <- JAGS_diagnostics_autocorrelation(
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
  JAGS_diagnostics_autocorrelation(
    fit,
    parameter = "theta",
    plot_type = "base",
    lags = 4
  )
  expect_lte(graphics::par("usr")[[3]], expected_min)
})


test_that("JAGS autocorrelation diagnostics preserve valid scalar lag limits", {

  fit <- .mock_public_diagnostics_fit()

  plot <- JAGS_diagnostics_autocorrelation(
    fit,
    parameter = "theta",
    plot_type = "ggplot",
    lags = 3
  )
  plot_layers <- ggplot2::ggplot_build(plot)$data
  for(layer in plot_layers){
    expect_equal(layer$x, 0:3)
  }

  zero_lag_plot <- JAGS_diagnostics_autocorrelation(
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


test_that("JAGS autocorrelation diagnostics cap output at available lags", {

  fit <- .mock_public_diagnostics_fit(n = 6)
  plot <- JAGS_diagnostics_autocorrelation(
    fit,
    parameter = "theta",
    plot_type = "ggplot",
    lags = 30
  )
  plot_layers <- ggplot2::ggplot_build(plot)$data

  expect_true(all(vapply(plot_layers, function(layer) nrow(layer) == 6L, logical(1))))
  for(layer in plot_layers){
    expect_equal(layer$x, 0:5)
  }
})
