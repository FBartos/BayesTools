skip_if_not_test_profile("unit")

.strict_bounded_plot_transform <- function(){

  list(
    fun = tanh,
    inv = function(x){
      if(any(!is.finite(x)) || any(x < -1 | x > 1)){
        stop("inverse called outside correlation support")
      }
      atanh(x)
    },
    jac = function(x) 1 - tanh(x)^2,
    output_support = c(-1, 1)
  )
}

test_that("bounded transformed grids do not invert out-of-support display coordinates", {

  transformation <- .strict_bounded_plot_transform()
  grid <- c(-1.5, -1, -.5, 0, .5, 1, 1.5)
  prior <- prior("normal", list(0, .4))
  density <- density(prior, x_seq = grid, transformation = transformation,
                     transformation_settings = TRUE)

  expect_equal(density$x, c(-.5, 0, .5))
  expect_equal(density$y, stats::dnorm(atanh(density$x), sd = .4) / (1 - density$x^2))
  expect_identical(attr(density, "x_range"), c(-1.5, 1.5))
  expect_error(transformation$inv(1.5), "outside correlation support", fixed = TRUE)

  builtin <- density(prior, x_seq = grid, transformation = "tanh",
                     transformation_settings = TRUE)
  expect_equal(builtin$x, density$x)
  expect_equal(builtin$y, density$y)
  expect_error(density(prior, x_range = c(2, 3), transformation = transformation,
                       transformation_settings = TRUE), "does not contain values", fixed = TRUE)
})

test_that("linear prior density overlays retain bounded image and wider display range", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(0, .4))),
    weights = c(theta = 1), n_grid = 512
  )
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    prior_density, n_points = 301, x_range = c(-1.5, 1.5),
    transformation = .strict_bounded_plot_transform(), transformation_settings = TRUE
  )$density
  expect_true(all(plot_data$x > -1 & plot_data$x < 1))
  expect_identical(attr(plot_data, "x_range"), c(-1.5, 1.5))
  expect_equal(plot_data$y, stats::dnorm(atanh(plot_data$x), sd = .4) / (1 - plot_data$x^2),
               tolerance = 1e-10)
})

test_that("closed transformed boundaries preserve continuous heights and atom locations", {

  transformation <- list(
    fun = function(x) x + 1,
    inv = function(x){
      if(any(x < 1 | x > 2)) stop("outside support")
      x - 1
    },
    jac = function(x) rep(1, length(x)),
    output_support = c(1, 2)
  )
  priors <- list(
    prior("uniform", list(0, 1), prior_weights = .5),
    prior("point", list(0), prior_weights = .2),
    prior("point", list(1), prior_weights = .3)
  )
  density <- density(priors[[1]], x_seq = c(0, 1, 1.5, 2, 3),
                     transformation = transformation, transformation_settings = TRUE,
                     truncate_end = FALSE)
  expect_equal(density$x, c(1, 1.5, 2))
  expect_equal(density$y, rep(1, 3))
  plot_data <- BayesTools:::.plot_data_prior_list.simple(
    priors, x_seq = NULL, x_range = c(0, 3), x_range_quant = NULL,
    n_points = 301, n_samples = 1000, force_samples = FALSE, individual = FALSE,
    transformation = transformation, transformation_arguments = NULL,
    transformation_settings = TRUE
  )
  expect_true(all(plot_data$density$x >= 1 & plot_data$density$x <= 2))
  expect_equal(c(plot_data$points1$x, plot_data$points2$x), c(1, 2))
  expect_equal(c(plot_data$points1$y, plot_data$points2$y), c(.2, .3))
  expect_identical(attr(plot_data$density, "x_range"), c(0, 3))

  linear_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior_mixture(priors)), weights = c(theta = 1), n_grid = 512
  )
  linear_plot <- BayesTools:::.prior_linear_density_to_plot_data(
    linear_density, n_points = 301, x_range = c(0, 3),
    transformation = transformation, transformation_settings = TRUE
  )
  expect_equal(c(linear_plot$points1$x, linear_plot$points2$x), c(1, 2))
  expect_equal(c(linear_plot$points1$y, linear_plot$points2$y), c(.2, .3))
})

test_that("bounded transformations preserve display limits in base and ggplot", {

  prior <- prior("normal", list(0, .4))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  graphics::par(xaxs = "i")
  plot(prior, xlim = c(-1.5, 1.5), transformation = .strict_bounded_plot_transform(),
       transformation_settings = TRUE)
  expect_equal(graphics::par("usr")[1:2], c(-1.5, 1.5))

  skip_if_not_installed("ggplot2")
  plot <- plot(prior, plot_type = "ggplot", xlim = c(-1.5, 1.5),
               transformation = .strict_bounded_plot_transform(),
               transformation_settings = TRUE)
  expect_equal(plot$scales$get_scales("x")$limits, c(-1.5, 1.5))
  built <- ggplot2::ggplot_build(plot)
  expect_true(all(built$data[[1]]$x >= -1 & built$data[[1]]$x <= 1))
})

test_that("density transformations use analytic exp_lin limits and omit saturated knots", {

  half_normal <- prior("normal", list(0, 1), truncation = list(0, Inf))
  identity_map <- density(half_normal, transformation = "exp_lin",
                          transformation_arguments = list(a = 1, b = 1))
  expect_true(all(is.finite(identity_map$y)))
  expect_true(all(is.finite(attr(identity_map, "y_range"))))
  # The source knot at zero maps to the density f(0) / exp(a) of exp(a) x.
  expect_equal(max(identity_map$y[identity_map$x == 0]),
               2 * stats::dnorm(0) / exp(1))
  root_map <- density(half_normal, transformation = "exp_lin",
                      transformation_arguments = list(a = 1, b = .5))
  expect_true(all(is.finite(root_map$y)))
  expect_true(all(root_map$y[root_map$x == 0] == 0))

  cauchy <- prior("cauchy", list(0, .707))
  conditions <- list()
  tanh_map <- withCallingHandlers(density(cauchy, transformation = "tanh"),
    BayesTools_prior_curve_unavailable = function(condition){
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    })
  expect_length(conditions, 1L)
  source <- density(cauchy)
  indices <- attr(tanh_map, "source_indices", exact = TRUE)
  expect_identical(tanh_map$x, tanh(source$x[indices]))
  expect_identical(conditions[[1L]]$unresolved_values, source$x[conditions[[1L]]$indices])
  expect_true(any(diff(tanh_map$x) == 0))
  expect_true(all(is.finite(tanh_map$y)))
  expect_true(all(abs(tanh_map$x) <= 1))
  interior <- abs(tanh_map$x) < .99
  expect_equal(
    tanh_map$y[interior],
    stats::dcauchy(atanh(tanh_map$x[interior]), 0, .707) / (1 - tanh_map$x[interior]^2),
    tolerance = 1e-10
  )

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  expect_warning(plot(cauchy, transformation = "tanh"), class = "BayesTools_prior_curve_unavailable")
  skip_if_not_installed("ggplot2")
  expect_warning(ggplot <- plot(cauchy, transformation = "tanh", plot_type = "ggplot"),
    class = "BayesTools_prior_curve_unavailable")
  expect_s3_class(ggplot, "ggplot")
  expect_error(.density.prior_transformation_grid(c(20, 21), c(1, 1), "tanh"),
    class = "BayesTools_transformation")
  custom <- list(fun = tanh, inv = atanh, jac = function(x) 1 - tanh(x)^2, output_support = c(-1, 1))
  expect_error(.density.prior_transformation_grid(c(0, 17, 17.001), c(1, 1, 1), custom),
    class = "BayesTools_nonmonotone_transformation")
})

test_that("transformation output support metadata is validated", {

  transformation <- .strict_bounded_plot_transform()
  for(support in list(c(1, -1), c(0, 0), c(NA_real_, 1), 1, "bounded")){
    transformation$output_support <- support
    expect_error(density(prior("normal", list(0, 1)), transformation = transformation),
                 "output_support", fixed = TRUE)
  }
})

test_that("exp_lin keeps the curve of a root map when the source density is infinite at zero", {

  roots <- list(
    gamma = list(prior = prior("gamma", list(.5, 1)), upper = Inf,
                 density = function(y) 2 * y * stats::dgamma(y^2, .5, 1)),
    beta  = list(prior = prior("beta", list(.5, .5)), upper = 1,
                 density = function(y) 2 * y * stats::dbeta(y^2, .5, .5)))
  for(root in roots){
    curves <- lapply(c(TRUE, FALSE), function(truncate_end){
      density(root$prior, transformation = "exp_lin",
              transformation_arguments = list(a = 0, b = .5), truncate_end = truncate_end)
    })
    for(curve in curves){
      # finite below the source's own upper boundary, where a Beta density stays infinite as untransformed
      below <- curve$x < root$upper
      expect_true(all(is.finite(curve$y[below])))
      inside <- curve$x > 0 & below
      expect_equal(curve$y[inside], root$density(curve$x[inside]), tolerance = 1e-8)
    }
    # the default closes the curve at the source zero with a zero height, as for finite boundaries
    expect_true(all(curves[[1L]]$y[curves[[1L]]$x == 0] == 0))
    expect_false(any(curves[[2L]]$x == 0))
    expect_identical(curves[[1L]]$x[curves[[1L]]$x > 0 & curves[[1L]]$x < root$upper],
                     curves[[2L]]$x[curves[[2L]]$x > 0 & curves[[2L]]$x < root$upper])
    expect_identical(curves[[1L]]$y[curves[[1L]]$x > 0 & curves[[1L]]$x < root$upper],
                     curves[[2L]]$y[curves[[2L]]$x > 0 & curves[[2L]]$x < root$upper])
  }
  gamma_curve <- density(roots$gamma$prior, transformation = "exp_lin",
                         transformation_arguments = list(a = 0, b = .5))
  expect_true(all(is.finite(gamma_curve$y)))
  expect_true(all(is.finite(attr(gamma_curve, "y_range"))))
  heights <- stats::approx(gamma_curve$x, gamma_curve$y, xout = c(.1, .5, 1))$y
  expect_equal(heights, roots$gamma$density(c(.1, .5, 1)), tolerance = 5e-3)
})

# ---------------------------------------------------------------------------- #
# KDE routes of the marginal and posterior plots use the checked transformation grid
# ---------------------------------------------------------------------------- #

.kde_marginal_samples <- function(x, parameter){

  class(x) <- c("marginal_posterior.simple", "marginal_posterior", class(x))
  out <- list(BayesTools:::.bt_meta_set(x, "atoms", posterior_atom_attribute()))
  names(out) <- parameter
  out
}
.kde_mixed_posteriors <- function(columns, prior_list){

  fit <- structure(coda::mcmc(columns), class = c("mcmc", "BayesTools_fit"))
  attr(fit, "prior_list") <- prior_list
  fit <- BayesTools:::.bt_attach_parameter_map(fit, monitor_names = NULL)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  BayesTools:::.bt_attach_fit_contract(fit)
}
# evaluates 'expr', returning its value and the warnings it signalled
.kde_collect_warnings <- function(expr){

  warnings <- list()
  value <- withCallingHandlers(expr, warning = function(condition){
    warnings[[length(warnings) + 1L]] <<- condition
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
.kde_has_warning <- function(collected, pattern = NULL, class = NULL){

  any(vapply(collected$warnings, function(condition){
    (is.null(pattern) || grepl(pattern, conditionMessage(condition), fixed = TRUE)) &&
      (is.null(class) || inherits(condition, class))
  }, logical(1)))
}
# renders both backends; the limits are taken from the retained, finite curve
.kde_expect_finite_plots <- function(draw, warning_pattern = NULL, warning_class = NULL){

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  base <- .kde_collect_warnings(draw("base"))
  expect_true(all(is.finite(graphics::par("usr"))))
  if(!is.null(warning_pattern) || !is.null(warning_class)){
    expect_true(.kde_has_warning(base, warning_pattern, warning_class))
  }

  gg <- .kde_collect_warnings(draw("ggplot"))
  if(!is.null(warning_pattern) || !is.null(warning_class)){
    expect_true(.kde_has_warning(gg, warning_pattern, warning_class))
  }
  built <- ggplot2::ggplot_build(gg$value)
  expect_true(all(is.finite(unlist(lapply(built$data, function(layer) c(layer$x, layer$y))))))
  expect_true(all(is.finite(built$layout$panel_params[[1]]$x.range)))
  expect_true(all(is.finite(built$layout$panel_params[[1]]$y.range)))
  invisible(NULL)
}

test_that("marginal KDE plots keep the finite part of the transformed curve", {

  skip_if_not_installed("ggplot2")
  set.seed(1)
  samples <- .kde_marginal_samples(stats::rgamma(400, 3, 3), "tau")
  log_map  <- list(fun = log, inv = exp, jac = function(x) 1 / x, output_support = c(-Inf, Inf))
  sqrt_map <- list(fun = sqrt, inv = function(x) x^2, jac = function(x) 1 / (2 * sqrt(x)), output_support = c(0, Inf))
  cases <- list(
    exp_lin_identity = list("exp_lin", NULL),
    exp_lin_power    = list("exp_lin", list(a = 0, b = 2)),
    log              = list(log_map, NULL),
    sqrt             = list(sqrt_map, NULL)
  )

  for(case in cases){
    # the KDE tail below zero has no image: the out-of-support points are dropped
    data <- .kde_collect_warnings(
      .plot_data_marginal_samples(samples, "tau", FALSE, 200, case[[1]], case[[2]], FALSE))
    density <- data$value$density1
    expect_true(.kde_has_warning(data, "NaNs produced"))
    expect_true(all(is.finite(density$x)) && all(is.finite(density$y)))
    expect_identical(attr(density, "x_range"), range(density$x))
    expect_identical(attr(density, "y_range"), c(0, max(density$y)))
    expect_lt(length(density$x), 200)

    .kde_expect_finite_plots(function(type){
      plot_marginal(samples, "tau", plot_type = type, transformation = case[[1]],
                    transformation_arguments = case[[2]])
    }, warning_pattern = "NaNs produced")
  }
})

test_that("marginal KDE plots keep the curve of the exp and lin transformations", {

  set.seed(1)
  tau <- stats::rgamma(400, 3, 3)
  samples <- .kde_marginal_samples(tau, "tau")
  kde <- stats::density(tau, n = 200)

  exp_data <- .plot_data_marginal_samples(samples, "tau", FALSE, 200, "exp", NULL, FALSE)$density1
  expect_equal(exp_data$x, exp(kde$x))
  expect_equal(exp_data$y, kde$y / exp(kde$x))

  lin_data <- .plot_data_marginal_samples(samples, "tau", FALSE, 200, "lin", list(a = 1, b = 2), FALSE)$density1
  expect_equal(lin_data$x, 1 + 2 * kde$x)
  expect_equal(lin_data$y, kde$y / 2)
})

test_that("marginal KDE plots of saturated tanh images drop the knots without an ordinate", {

  skip_if_not_installed("ggplot2")
  samples <- .kde_marginal_samples(seq(-30, 30, length.out = 201), "theta")

  collected <- .kde_collect_warnings(.plot_data_marginal_samples(samples, "theta", FALSE, 512, "tanh", NULL, FALSE))
  expect_true(.kde_has_warning(collected, class = "BayesTools_prior_curve_unavailable"))
  data <- collected$value$density1
  expect_true(all(is.finite(data$x)) && all(is.finite(data$y)))
  expect_true(all(abs(data$x) < 1))
  expect_identical(attr(data, "x_range"), range(data$x))
  expect_identical(attr(data, "y_range"), c(0, max(data$y)))
  kde <- stats::density(seq(-30, 30, length.out = 201), n = 512)
  index <- match(data$x, tanh(kde$x))
  expect_false(anyNA(index))
  interior <- abs(data$x) < .9
  expect_equal(data$y[interior], kde$y[index[interior]] / (1 - data$x[interior]^2), tolerance = 1e-10)

  .kde_expect_finite_plots(function(type){
    plot_marginal(samples, "theta", plot_type = type, transformation = "tanh")
  }, warning_class = "BayesTools_prior_curve_unavailable")
})

test_that("posterior KDE plots of saturated tanh images drop the knots without an ordinate", {

  skip_if_not_installed("ggplot2")
  theta <- as_mixed_posteriors(
    .kde_mixed_posteriors(cbind(theta = seq(0, 25, length.out = 100)),
                          list(theta = prior("normal", list(0, 10)))),
    "theta")

  data <- .kde_collect_warnings(.plot_data_samples.simple(theta, "theta", 256, "tanh", NULL, FALSE))
  expect_true(.kde_has_warning(data, class = "BayesTools_prior_curve_unavailable"))
  density <- data$value$density
  expect_true(all(is.finite(density$x)) && all(is.finite(density$y)))
  expect_lt(length(density$x), 256)
  expect_identical(attr(density, "x_range"), range(density$x))
  expect_identical(attr(density, "y_range"), c(0, max(density$y)))

  .kde_expect_finite_plots(function(type){
    plot_posterior(theta, "theta", plot_type = type, transformation = "tanh", n_points = 256)
  }, warning_class = "BayesTools_prior_curve_unavailable")

  # a transformation without saturation keeps every knot
  lin <- .plot_data_samples.simple(theta, "theta", 256, "lin", list(a = 1, b = 2), FALSE)$density
  raw <- .plot_data_samples.simple(theta, "theta", 256, NULL, NULL, FALSE)$density
  expect_equal(lin$x, 1 + 2 * raw$x)
  expect_equal(lin$y, raw$y / 2)
})

test_that("posterior plots of a very wide prior under exp omit only the unavailable prior curve", {

  skip_if_not_installed("ggplot2")
  set.seed(1)
  theta <- as_mixed_posteriors(
    .kde_mixed_posteriors(cbind(theta = stats::rnorm(200, .5)),
                          list(theta = prior("normal", list(0, 1000)))),
    "theta")
  draw <- function(type, prior, ...){
    plot_posterior(theta, "theta", plot_type = type, prior = prior, n_points = 64,
                   transformation = "exp", transformation_settings = TRUE, ...)
  }
  layer_data <- function(plot) ggplot2::ggplot_build(plot)$data

  # the 0.995 prior quantile has no finite image under exp: the posterior is
  # drawn alone, as with prior = FALSE
  collected <- .kde_collect_warnings(draw("ggplot", TRUE))
  curve_warnings <- Filter(function(condition) inherits(condition, "BayesTools_prior_curve_unavailable"),
                           collected$warnings)
  expect_length(curve_warnings, 1L)
  expect_s3_class(curve_warnings[[1L]], "BayesTools_plot_condition")
  expect_match(conditionMessage(curve_warnings[[1L]]), "'theta'", fixed = TRUE)
  expect_match(conditionMessage(curve_warnings[[1L]]),
               "automatic transformed prior display range is unavailable", fixed = TRUE)
  expect_match(conditionMessage(curve_warnings[[1L]]), "'xlim'", fixed = TRUE)
  expect_s3_class(collected$value, "ggplot")
  posterior_alone <- layer_data(draw("ggplot", FALSE))
  expect_length(posterior_alone, 1L)
  expect_identical(layer_data(collected$value), posterior_alone)

  .kde_expect_finite_plots(function(type) draw(type, TRUE), warning_class = "BayesTools_prior_curve_unavailable")

  # explicit display limits keep the prior curve
  with_limits <- .kde_collect_warnings(draw("ggplot", TRUE, xlim = c(.05, 20)))
  expect_false(.kde_has_warning(with_limits, class = "BayesTools_prior_curve_unavailable"))
  expect_length(layer_data(with_limits$value), 2L)
})

test_that("independent factor posterior KDE plots of saturated tanh images drop the knots without an ordinate", {

  skip_if_not_installed("ggplot2")
  formula <- JAGS_formula(~ 0 + fac, "mu", data.frame(fac = factor(c("A", "B", "C"))),
    prior_list = list(fac = prior_factor("normal", list(0, 1), contrast = "independent")))
  columns <- BayesTools:::.JAGS_prior_factor_names("mu_fac", formula$prior_list$mu_fac)
  draws <- matrix(rep(seq(0, 25, length.out = 100), 3), 100, dimnames = list(NULL, columns))
  draws[, 2] <- draws[, 2] * .9
  draws[, 3] <- draws[, 3] * 1.1
  factor_samples <- as_mixed_posteriors(.kde_mixed_posteriors(draws, formula$prior_list), "mu_fac")

  data <- .kde_collect_warnings(.plot_data_samples.factor(factor_samples, "mu_fac", 256, "tanh", NULL, FALSE))
  expect_true(.kde_has_warning(data, class = "BayesTools_prior_curve_unavailable"))
  curves <- data$value[grepl("^density", names(data$value))]
  expect_length(curves, 3)
  for(curve in curves){
    expect_true(all(is.finite(curve$x)) && all(is.finite(curve$y)))
    expect_lt(length(curve$x), 256)
    expect_identical(attr(curve, "x_range"), range(curve$x))
  }

  .kde_expect_finite_plots(function(type){
    plot_posterior(factor_samples, "mu_fac", plot_type = type, transformation = "tanh", n_points = 256)
  }, warning_class = "BayesTools_prior_curve_unavailable")
})
