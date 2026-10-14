skip_if_not_test_profile("unit")

test_that("orthonormal point mixtures retain every continuous transformed knot", {
  normal <- prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  point <- prior_factor("spike", list(0), contrast = "orthonormal")
  attr(normal, "levels") <- attr(point, "levels") <- 3L
  priors <- list(normal = normal, point = point)
  data <- suppressMessages(.plot_data_prior_list.simple(priors, x_seq = NULL,
    x_range = c(.01, 5), x_range_quant = NULL, n_points = 500L, n_samples = 10000L,
    force_samples = FALSE, individual = FALSE, transformation = "exp",
    transformation_arguments = NULL, transformation_settings = TRUE))
  expect_length(data$density$x, 500L)
  expect_length(data$density$y, 500L)
  # The declared orthonormal row is a projection of I - J/3; its variance
  # under independent unit-variance coefficients is the row's squared norm.
  contrast_row <- contr.orthonormal(3L)[1L, ]
  variance <- sum(contrast_row^2)
  expect_equal(variance, 2/3, tolerance = 1e-15)
  expect_equal(as.numeric(data$density$y), .5 * stats::dnorm(log(data$density$x),
    mean = 0, sd = sqrt(variance)) / data$density$x, tolerance = 1e-13)
  expect_identical(as.numeric(data$points1$x), 1)
  expect_identical(as.numeric(data$points1$y), .5)
  figure <- suppressMessages(plot_prior_list(priors, plot_type = "ggplot", transformation = "exp",
    transformation_settings = TRUE, xlim = c(.01, 5)))
  expect_s3_class(figure, "ggplot")
  expect_true(any(vapply(ggplot2::ggplot_build(figure)$data, nrow, integer(1)) == 500L))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(suppressMessages(plot_prior_list(priors, transformation = "exp",
    transformation_settings = TRUE, xlim = c(.01, 5))))
  ordinary <- .plot_data_prior_list.simple(list(prior("normal", list(0, 1)), prior("point", list(0))),
    x_seq = NULL, x_range = c(-3, 3), x_range_quant = NULL, n_points = 31L,
    n_samples = 10000L, force_samples = FALSE, individual = FALSE, transformation = NULL,
    transformation_arguments = NULL, transformation_settings = FALSE)
  expect_length(ordinary$density$x, 31L)
  expect_equal(as.numeric(ordinary$density$y), .5 * stats::dnorm(ordinary$density$x), tolerance = 1e-14)
})

test_that("spike plotting preserves distinct nearby atoms and exact duplicates", {

  for(locations in list(c(0, 1e-10), c(1e6, 1e6 + .005), c(1, 1))){
    priors <- lapply(locations, function(location) prior("spike", list(location)))
    samples <- rep(locations, each = 2)
    attr(samples, "prior_list") <- priors
    samples <- .bt_meta_set(samples, "atoms", posterior_atom_attribute(
      data.frame(x = locations, mass = c(.5, .5))
    ))
    plot_data <- BayesTools:::.plot_data_samples.simple(
      list(theta = samples), "theta", 64, NULL, NULL, FALSE
    )
    points <- plot_data[vapply(plot_data, inherits, logical(1), "density.prior.point")]
    expect_equal(sum(vapply(points, `[[`, numeric(1), "y")), 1)
    expect_equal(vapply(points, `[[`, numeric(1), "x"), unique(locations), ignore_attr = TRUE)
    expect_equal(
      vapply(points, `[[`, numeric(1), "y"),
      rep(1 / length(unique(locations)), length(unique(locations))),
      ignore_attr = TRUE
    )
  }

  plot_data <- list(x = c(0, 0, 1), y = c(.2, .3, .5))
  layer <- BayesTools:::.geom_prior.point(plot_data)
  expect_equal(layer$data$x, plot_data$x)
  expect_equal(layer$data$yend, plot_data$y)
})

test_that("weight prior plotting preserves adjacent representable atom locations", {

  locations <- c(1, 1 + .Machine$double.eps, 1, .5)
  components <- lapply(locations, function(location){
    list(type = "point", location = location, weight = 1 / 4)
  })
  plotted <- BayesTools:::.plot_data_prior_weightparameter_components(
    components, parameter = "omega", n_points = 10
  )
  expect_length(plotted, 3L)
  expect_equal(unname(vapply(plotted, function(x) x$y, numeric(1))), c(1 / 4, 1 / 2, 1 / 4))
  expect_equal(unname(vapply(plotted, function(x) x$x, numeric(1))), sort(unique(locations)), tolerance = 0)
})

test_that("PET-PEESE plot functions retain truncated-normal tail probabilities", {

  p <- prior_PET("normal", list(0, 1), truncation = list(40, 41))
  functions <- BayesTools:::.petpeese_prior_simple_functions(p)
  probabilities <- c(.1, .5, .9)
  quantiles <- functions$quant(probabilities)
  expect_true(all(is.finite(quantiles)))
  expect_equal(functions$cdf(quantiles), probabilities, tolerance = 1e-10)
  expect_equal(functions$ccdf(quantiles), 1 - probabilities, tolerance = 1e-10)
  expected_pdf <- exp(stats::dnorm(quantiles, log = TRUE) -
    stats::pnorm(40, lower.tail = FALSE, log.p = TRUE))
  # The omitted upper tail has relative mass below exp(-40).
  expect_equal(functions$pdf(quantiles), expected_pdf, tolerance = 1e-12)
})

# ============================================================================ #
# TEST FILE: Prior Plot Data Semantics
# ============================================================================ #
#
# PURPOSE:
#   Semantic ggplot-layer tests for plot.prior outputs. These complement the
#   visual-profile vdiffr snapshots in test-priors-plot.R.
#
# TAGS: @evaluation, @unit, @priors, @plots
# ============================================================================ #

prior_plot_layer_data <- function(plot) {
  ggplot2::ggplot_build(plot)$data
}

test_that("continuous prior ggplot data matches analytic density and labels", {
  x <- seq(-1, 1, length.out = 5)
  p <- prior("normal", list(mean = 0, sd = 1))

  g <- plot(
    p,
    plot_type = "ggplot",
    x_seq = x,
    xlab = "Effect",
    ylab = "Mass density",
    main = "Normal prior"
  )
  layers <- prior_plot_layer_data(g)

  expect_s3_class(g, "ggplot")
  expect_length(layers, 1)
  expect_equal(layers[[1]]$x, x)
  expect_equal(layers[[1]]$y, stats::dnorm(x), tolerance = 1e-12)
  expect_equal(g$labels$title, "Normal prior")
  expect_equal(g$scales$get_scales("x")$name, "Effect")
  expect_equal(g$scales$get_scales("y")$name, "Mass density")
})

test_that("continuous prior ggplot density preserves selected x grid and bounded mass", {
  x <- seq(-4, 4, length.out = 401)
  p <- prior("normal", list(mean = 0, sd = 1))

  layer <- prior_plot_layer_data(plot(p, plot_type = "ggplot", x_seq = x))[[1]]
  dx <- diff(layer$x)
  trapezoid_mass <- sum(dx * (head(layer$y, -1) + tail(layer$y, -1)) / 2)

  expect_equal(layer$x, x)
  expect_true(all(dx > 0))
  expect_equal(trapezoid_mass, diff(stats::pnorm(c(-4, 4))), tolerance = 1e-4)
  expect_true(all(layer$y >= 0))
  expect_equal(layer$y[which.max(layer$y)], stats::dnorm(0), tolerance = 1e-12)
})

test_that("transformed continuous prior ggplot data applies the Jacobian", {
  original_x <- seq(-1, 1, length.out = 5)
  transformed_x <- exp(original_x)

  g <- plot(
    prior("normal", list(mean = 0, sd = 1)),
    plot_type = "ggplot",
    x_seq = original_x,
    transformation = "exp"
  )
  layer <- prior_plot_layer_data(g)[[1]]

  expect_equal(layer$x, transformed_x)
  expect_equal(layer$y, stats::dnorm(original_x) / transformed_x, tolerance = 1e-12)
  expect_true(all(layer$x > 0))
})

test_that("discrete and point prior ggplot layers preserve exact probability mass", {
  discrete_layer <- prior_plot_layer_data(
    plot(prior("bernoulli", list(.33)), plot_type = "ggplot")
  )[[1]]

  expect_equal(discrete_layer$x, c(0, 1))
  expect_equal(discrete_layer$y, c(.67, .33), tolerance = 1e-12)
  expect_probability_mass(discrete_layer$y)

  point_layer <- prior_plot_layer_data(
    plot(prior("point", list(1.5)), plot_type = "ggplot")
  )[[1]]

  expect_equal(point_layer$x, 1.5)
  expect_equal(point_layer$xend, 1.5)
  expect_equal(point_layer$y, 0)
  expect_equal(point_layer$yend, 1)
})

test_that("point prior ggplot layer honors explicit probability scaling", {
  point_layer <- prior_plot_layer_data(
    plot(prior("point", list(location = -1)), plot_type = "ggplot", scale_y2 = 3)
  )[[1]]

  expect_equal(point_layer$x, -1)
  expect_equal(point_layer$xend, -1)
  expect_equal(point_layer$y, 0)
  expect_equal(point_layer$yend, 3)
  expect_probability_mass(point_layer$yend / 3)
})

test_that("spike-and-slab prior ggplot layers separate slab density and spike mass", {
  x <- seq(-2, 2, length.out = 7)
  p <- prior_spike_and_slab(prior("normal", list(mean = 0, sd = 1)))

  g <- plot(p, plot_type = "ggplot", x_seq = x)
  layers <- prior_plot_layer_data(g)

  expect_length(layers, 2)
  expect_equal(layers[[1]]$x, x)
  expect_equal(layers[[1]]$y, .5 * stats::dnorm(x), tolerance = 1e-12)
  expect_equal(layers[[2]]$x, 0)
  expect_equal(layers[[2]]$xend, 0)
  expect_equal(layers[[2]]$y, 0)
  expect_equal(layers[[2]]$yend / attr(g, "scale_y2"), .5, tolerance = 1e-12)
})

test_that("weightfunction prior ggplot layers preserve step locations and intervals", {
  p <- prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .25)))
  layers <- prior_plot_layer_data(plot(p, plot_type = "ggplot"))

  expect_length(layers, 2)

  expected_x <- c(0, .05, .05, 1)
  expected_y <- c(1, 1, .25, .25)
  expected_lower <- expected_y
  expected_upper <- expected_y

  expect_equal(layers[[1]]$x, c(expected_x, rev(expected_x)))
  expect_equal(layers[[1]]$y, c(expected_lower, rev(expected_upper)), tolerance = 1e-12)
  expect_equal(layers[[2]]$x, expected_x)
  expect_equal(layers[[2]]$y, expected_y, tolerance = 1e-12)
})

test_that("factor prior ggplot data uses factor-contrast density scale", {
  p <- prior_factor("normal", list(mean = 0, sd = 2), contrast = "treatment")
  p$parameters$K <- 3
  x <- seq(-4, 4, length.out = 9)

  layer <- prior_plot_layer_data(plot(p, plot_type = "ggplot", x_seq = x))[[1]]

  expect_equal(layer$x, x)
  expect_equal(layer$y, stats::dnorm(x, mean = 0, sd = 2), tolerance = 1e-12)
})

test_that("geom_prior contributes exact density and point-mass layers", {
  x <- seq(-1, 1, length.out = 5)

  density_plot <- ggplot2::ggplot() +
    geom_prior(prior("normal", list(mean = 0, sd = 1)), x_seq = x)
  density_layer <- ggplot2::ggplot_build(density_plot)$data[[1]]

  expect_equal(density_layer$x, x)
  expect_equal(density_layer$y, stats::dnorm(x), tolerance = 1e-12)

  point_plot <- ggplot2::ggplot() +
    geom_prior(prior("point", list(location = 0.5)), scale_y2 = 2)
  point_layer <- ggplot2::ggplot_build(point_plot)$data[[1]]

  expect_equal(point_layer$x, 0.5)
  expect_equal(point_layer$xend, 0.5)
  expect_equal(point_layer$y, 0)
  expect_equal(point_layer$yend, 2)
})

test_that("prior plot data rejects invalid plotting options before rendering", {
  p <- prior("normal", list(mean = 0, sd = 1))

  expect_error(plot(p, plot_type = "grid"), "'plot_type'")
  expect_error(plot(p, plot_type = "ggplot", show_figures = 1.5), "'show_figures'")
  expect_error(
    plot(p, plot_type = "ggplot", transformation = "not-a-transformation"),
    "function or character string"
  )
  expect_error(ggplot2::ggplot() + geom_prior(p, scale_y2 = -1), "'scale_y2'")
})


test_that("base probability axes fit panels without moving subsequent overlays", {

  grDevices::pdf(NULL, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2))
  original_mar <- graphics::par("mar")
  priors <- list(
    prior("normal", list(mean = 0, sd = 1)),
    prior("point", list(location = 0))
  )

  for (panel in seq_len(2L)) {
    plot_prior_list(priors)
    expect_equal(graphics::par("mfg")[1:2], c(1L, panel))
    expect_equal(graphics::par("mar"), original_mar)

    # The title is three lines beyond the plot's right edge, inside its panel.
    right_edge <- graphics::grconvertX(graphics::par("usr")[[2L]], "user", "nfc")
    right_space <- (1 - right_edge) * graphics::par("fin")[[1L]]
    expect_gt(right_space, 3 * graphics::par("csi") * graphics::par("mex"))

    overlay_before <- graphics::grconvertX(c(-1, 0, 1), "user", "ndc")
    lines_prior_list(priors, col = "red", lty = 2)
    expect_equal(graphics::grconvertX(c(-1, 0, 1), "user", "ndc"), overlay_before)
  }

  plot(priors[[1L]])
  right_edge <- graphics::grconvertX(graphics::par("usr")[[2L]], "user", "nfc")
  expect_equal(
    (1 - right_edge) * graphics::par("fin")[[1L]],
    original_mar[[4L]] * graphics::par("csi") * graphics::par("mex")
  )

  custom_mar <- original_mar
  custom_mar[[4L]] <- 6.1
  graphics::par(mar = custom_mar)
  plot_prior_list(priors)
  expect_equal(graphics::par("mar"), custom_mar)
})

test_that("weightfunction rescale keeps two-cut coordinates aligned", {

  two_cut <- list(
    x     = c(0, 1),
    y     = c(1, 1),
    y_lCI = c(1, 1),
    y_uCI = c(1, 1)
  )
  x_at <- BayesTools:::.weightfunction_plot_x_at(two_cut, TRUE)
  expect_equal(x_at, c(0, 1))
  expect_length(x_at, length(two_cut$y))

  three_cut <- list(
    x = c(0, .05, .05, 1),
    y = c(1, 1, .25, .25)
  )
  x_at3 <- BayesTools:::.weightfunction_plot_x_at(three_cut, TRUE)
  expect_equal(x_at3, c(0, .5, .5, 1))
  expect_length(x_at3, length(three_cut$y))
})

test_that("individual weightfunction and PET overlays use the selected component", {

  weight <- prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .25)))
  selected <- density(weight, individual = TRUE)[[1]]
  weight_plot <- ggplot2::ggplot() +
    geom_prior(weight, individual = TRUE, show_parameter = 1)
  weight_layer <- ggplot2::ggplot_build(weight_plot)$data[[1]]

  expect_equal(unique(weight_layer$x), unique(selected$x[selected$y != 0]))
  expect_equal(weight_layer$yend, selected$y[selected$y != 0])

  pet <- prior_PET("normal", list(0, 1))
  x <- seq(0, 2, length.out = 5)
  pet_data <- density(pet, individual = TRUE, x_seq = x)
  pet_plot <- ggplot2::ggplot() + geom_prior(pet, individual = TRUE, x_seq = x)
  pet_layer <- ggplot2::ggplot_build(pet_plot)$data[[1]]

  expect_equal(pet_layer$x, pet_data$x)
  expect_equal(pet_layer$y, pet_data$y, tolerance = 1e-12)
})

test_that("individual weight-function index k is the k-th omega of the summary tables", {

  # Cumulative weights: omega_k = sum_{j >= k} theta_j with theta ~
  # Dirichlet(alpha), so by Dirichlet aggregation omega_1 = 1 and
  # omega_k ~ Beta(sum(alpha[k:J]), sum(alpha[1:(k - 1)])) for k >= 2.
  cases <- list(
    list(prior = prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3))),
         alpha = c(1, 2, 3)),
    list(prior = prior_weightfunction("two-sided", c(.05, .10), wf_cumulative(c(2, 1, 4))),
         alpha = c(2, 1, 4))
  )

  drawn <- list()
  record <- function(plot_data, ...){
    drawn[[length(drawn) + 1L]] <<- plot_data
    invisible(NULL)
  }
  testthat::local_mocked_bindings(
    .plot.prior.simple  = function(x, plot_type, plot_data, ...) record(plot_data),
    .plot.prior.point   = function(x, plot_type, plot_data, ...) record(plot_data),
    .lines.prior.simple = function(plot_data, ...) record(plot_data),
    .lines.prior.point  = function(plot_data, ...) record(plot_data),
    .package = "BayesTools"
  )
  expect_omega <- function(x, y, k, alpha, info){
    if(k == 1L){
      expect_equal(unique(x), 1, info = info)
      return(invisible())
    }
    inside <- x > 0 & x < 1
    expect_gt(sum(inside), 10)
    expect_equal(
      y[inside],
      stats::dbeta(x[inside], sum(alpha[k:length(alpha)]), sum(alpha[seq_len(k - 1L)])),
      tolerance = 1e-10,
      info = info
    )
  }

  for(case in cases){
    p <- case$prior
    J <- nrow(p$bins)
    # The k-th omega printed by the summary tables (one-sided p-value scale,
    # as in RoBMA) maps to the prior's k-th weight, and its interval is the
    # prior's k-th bin (two-sided bins are twice the one-sided p-values).
    table_cuts <- weightfunctions_mapping(list(p), cuts_only = TRUE, one_sided = TRUE)
    expect_equal(weightfunctions_mapping(list(p), one_sided = TRUE)[[1]][seq_len(J)], seq_len(J))

    for(k in seq_len(J)){
      info <- paste(p$side, "k =", k)
      table_interval <- table_cuts[c(k, k + 1L)]
      if(p$side == "two-sided"){
        table_interval <- pmin(2 * table_interval, 1)
      }
      expect_equal(c(p$bins$lower[k], p$bins$upper[k]), table_interval, info = info)

      drawn <- list()
      plot(p, individual = TRUE, show_figures = k)
      lines(p, individual = TRUE, show_parameter = k)
      expect_length(drawn, 2L)
      for(plot_data in drawn){
        expect_equal(attr(plot_data, "steps"), table_interval, info = info)
        expect_omega(plot_data$x, plot_data$y, k, case$alpha, info)
      }

      layer <- ggplot2::ggplot_build(
        ggplot2::ggplot() + geom_prior(p, individual = TRUE, show_parameter = k)
      )$data[[1]]
      expect_omega(layer$x, layer$y, k, case$alpha, info)
    }

    # The default omits the reference weight fixed at 1.
    drawn <- list()
    plot(p, individual = TRUE)
    expect_equal(
      lapply(drawn, attr, "steps"),
      lapply(2:J, function(k) c(p$bins$lower[k], p$bins$upper[k]))
    )

    expect_error(
      plot(p, individual = TRUE, show_figures = J + 1L),
      paste0("'show_figures' must be between -", J, " and ", J,
             ", excluding 0, for a weight function with ", J, " publication weights."),
      fixed = TRUE
    )
    expect_error(
      lines(p, individual = TRUE, show_parameter = J + 1L),
      paste0("'show_parameter' must be between 1 and ", J,
             " for a weight function with ", J, " publication weights."),
      fixed = TRUE
    )
  }
})

test_that("individual posterior omega figures follow the summary-table columns", {

  priors <- list(
    prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1))),
    .set_prior_model_weight(prior_weightfunction("one-sided", c(.01), wf_cumulative(c(1, 1))), 0),
    prior_none()
  )
  omega_cuts <- weightfunctions_mapping(priors, cuts_only = TRUE)
  models_ind <- rep(c(1, 3), c(60, 40))
  omega <- matrix(1, nrow = length(models_ind), ncol = 4)
  omega[models_ind == 1, 3] <- seq(.2, .9, length.out = 60)
  colnames(omega) <- .weightfunction_omega_names(omega_cuts)
  omega <- .model_probability_component_set(omega, models_ind, priors, c(.6, 0, .4))
  attr(omega, "prior_list") <- priors
  omega <- .bt_meta_set(omega, "atoms", .posterior_atoms_from_priors(
    priors,
    c(.6, 0, .4),
    n_columns = ncol(omega),
    column_names = colnames(omega),
    null_location = 1
  ))
  class(omega) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")
  samples <- list(omega = omega)
  class(samples) <- c("mixed_posteriors", "list")

  for(k in seq_len(ncol(omega))){
    figure <- plot_posterior(samples, "omega", plot_type = "ggplot", individual = TRUE, show_figures = k)
    expect_identical(names(figure), colnames(omega)[k])
  }
})

test_that("plot_prior_list individual PET-PEESE uses the parameter density", {

  pet <- prior_PET("normal", list(0, 1))
  x <- seq(0, 2, length.out = 5)
  pet_data <- density(pet, individual = TRUE, x_seq = x)
  g <- plot_prior_list(
    list(pet),
    individual = TRUE,
    plot_type  = "ggplot",
    x_seq      = x
  )
  layer <- ggplot2::ggplot_build(g)$data[[1]]

  expect_s3_class(g, "ggplot")
  expect_equal(layer$x, pet_data$x)
  expect_equal(layer$y, pet_data$y, tolerance = 1e-12)
})

test_that("individual PET-PEESE prior lists refuse PET and PEESE together", {

  pet   <- prior_PET("normal", list(0, 1))
  peese <- prior_PEESE("normal", list(0, 1))
  none  <- prior_none()
  message <- "'individual = TRUE' is unavailable for prior lists that contain both PET and PEESE priors"

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for(priors in list(list(pet, peese), list(pet, peese, none))){
    expect_error(plot_prior_list(priors, individual = TRUE), message, fixed = TRUE)
    expect_error(plot_prior_list(priors, individual = TRUE, plot_type = "ggplot"), message, fixed = TRUE)
    expect_error(lines_prior_list(priors, individual = TRUE), message, fixed = TRUE)
    expect_error(geom_prior_list(priors, individual = TRUE), message, fixed = TRUE)
  }
})

test_that("individual PET and PEESE prior lists agree across the plot, lines and geom functions", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  mu <- prior("normal", list(0, 1))
  x  <- seq(0, 2, length.out = 5)

  for(coefficient in c("PET", "PEESE")){
    bias <- if(coefficient == "PET") prior_PET("normal", list(0, 1)) else prior_PEESE("normal", list(0, 1))
    expected <- density(bias, individual = TRUE, x_seq = x)

    plot_layer <- ggplot2::ggplot_build(
      plot_prior_list(list(bias), individual = TRUE, plot_type = "ggplot", x_seq = x)
    )$data[[1]]
    geom_layer <- ggplot2::ggplot_build(
      ggplot2::ggplot() + geom_prior_list(list(bias), individual = TRUE, x_seq = x)
    )$data[[1]]
    drawn <- NULL
    testthat::local_mocked_bindings(
      .lines.prior.simple = function(plot_data, ...) drawn <<- plot_data
    )
    plot(0, xlim = c(0, 2), ylim = c(0, 1))
    lines_prior_list(list(bias), individual = TRUE, x_seq = x)

    expect_equal(plot_layer$x, expected$x)
    expect_equal(plot_layer$y, expected$y, tolerance = 1e-12)
    expect_equal(geom_layer$x, plot_layer$x)
    expect_equal(geom_layer$y, plot_layer$y)
    expect_equal(drawn$x, expected$x)
    expect_equal(drawn$y, expected$y, tolerance = 1e-12)

    # the mu priors are used only by the full regression bands
    for(call in list(
      function() plot_prior_list(list(bias), individual = TRUE, prior_list_mu = list(mu)),
      function() lines_prior_list(list(bias), individual = TRUE, prior_list_mu = list(mu)),
      function() geom_prior_list(list(bias), individual = TRUE, prior_list_mu = list(mu))
    )){
      expect_error(call(), "'prior_list_mu' is required only for PET-PEESE plots.", fixed = TRUE)
    }
    expect_error(lines_prior_list(list(bias)), "prior_list_mu")
    expect_error(geom_prior_list(list(bias)), "prior_list_mu")
  }
})

test_that("individual PET-PEESE lists with prior_none() models draw the coefficient density and a spike at zero", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  mu   <- prior("normal", list(0, 1))
  none <- prior_none(prior_weights = 2)

  # the plot data each function draws
  captured <- list()
  testthat::local_mocked_bindings(
    .plot_prior_list.both = function(plot_data, ...) captured$plot <<- plot_data,
    .lines.prior.simple   = function(plot_data, ...) captured$lines_density <<- plot_data,
    .lines.prior.point    = function(plot_data, ...) captured$lines_point <<- plot_data
  )

  for(bias in list(prior_PET("normal", list(0, 1), prior_weights = 1),
                   prior_PEESE("normal", list(0, 1), prior_weights = 1))){
    priors <- list(bias, none)

    # the coefficient density with its model probability and an atom at zero
    # with the probability of the models without the coefficient
    plot_prior_list(priors, individual = TRUE, xlim = c(0, 2), n_points = 21)
    expect_named(captured$plot, c("density", "points1"))
    expect_s3_class(captured$plot$points1, "density.prior.point")
    expect_equal(captured$plot$points1$x, 0)
    expect_equal(captured$plot$points1$y, 2 / 3)
    single <- density(bias, individual = TRUE, x_seq = captured$plot$density$x[-1], truncate_end = FALSE)
    expect_equal(captured$plot$density$y[-1], single$y / 3, tolerance = 1e-12)

    lines_prior_list(priors, individual = TRUE, xlim = c(0, 2), n_points = 21)
    expect_equal(captured$lines_density$x, captured$plot$density$x)
    expect_equal(captured$lines_density$y, captured$plot$density$y)
    expect_equal(captured$lines_point$x, 0)
    expect_equal(captured$lines_point$y, 2 / 3)

    geom_gg <- ggplot2::ggplot() + geom_prior_list(priors, individual = TRUE, xlim = c(0, 2))
    expect_s3_class(geom_gg, "ggplot")
    geom_layers <- unname(vapply(geom_gg$layers, function(layer) class(layer$geom)[1], character(1)))
    expect_equal(sort(geom_layers), c("GeomLine", "GeomSegment"))
    expect_equal(geom_gg$layers[[which(geom_layers == "GeomSegment")]]$data$x, 0)
  }

  # the full regression bands keep requiring and using the mu priors
  expect_error(plot_prior_list(list(prior_PET("normal", list(0, 1)), none)), "prior_list_mu")
  expect_error(lines_prior_list(list(prior_PET("normal", list(0, 1)), none)), "prior_list_mu")
  expect_error(geom_prior_list(list(prior_PET("normal", list(0, 1)), none)), "prior_list_mu")
  bands <- plot_prior_list(list(prior_PET("normal", list(0, 1)), none), prior_list_mu = list(mu, mu), plot_type = "ggplot")
  expect_s3_class(bands, "ggplot")
})

test_that("Bernoulli prior-list components are atoms at zero and one", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  bernoulli <- prior("bernoulli", list(.3))

  # an atom-only list renders as probability masses, as plot() of the prior does
  plot <- plot_prior_list(list(bernoulli), plot_type = "ggplot")
  geoms <- unname(vapply(plot$layers, function(layer) class(layer$geom)[1], character(1)))
  segments <- do.call(rbind, lapply(plot$layers, function(layer) layer$data))
  segments <- segments[order(segments$x), ]
  expect_equal(geoms, c("GeomSegment", "GeomSegment"))
  expect_equal(segments$x, c(0, 1))
  expect_equal(segments$yend, c(.7, .3))
  expect_equal(plot$scales$get_scales("y")$name, "Probability")
  direct <- plot(bernoulli, plot_type = "ggplot")
  expect_equal(direct$scales$get_scales("y")$name, "Probability")

  # the weighted masses of several Bernoulli components sum by location
  data <- .plot_data_prior_list.simple(
    list(bernoulli, prior("bernoulli", list(.5), prior_weights = 3)),
    NULL, NULL, NULL, 101, 10000, FALSE, FALSE, NULL, NULL, FALSE)
  expect_named(data, c("points1", "points2"))
  masses <- stats::setNames(c(data$points1$y, data$points2$y), c(data$points1$x, data$points2$x))
  expect_equal(masses[["0"]], (.7 + 3 * .5) / 4)
  expect_equal(masses[["1"]], (.3 + 3 * .5) / 4)

  # the masses are the normalized pmf of a truncated Bernoulli
  truncated <- prior("bernoulli", list(.3), truncation = list(.5, 1))
  data <- .plot_data_prior_list.simple(list(truncated), NULL, NULL, NULL, 101, 10000, FALSE, FALSE, NULL, NULL, FALSE)
  expect_named(data, "points1")
  expect_equal(data$points1$x, 1)
  expect_equal(data$points1$y, 1)
  data <- .plot_data_prior_list.simple(list(truncated, bernoulli), NULL, NULL, NULL, 101, 10000, FALSE, FALSE, NULL, NULL, FALSE)
  masses <- stats::setNames(c(data$points1$y, data$points2$y), c(data$points1$x, data$points2$x))
  expect_equal(masses[["0"]], .7 / 2)
  expect_equal(masses[["1"]], (1 + .3) / 2)

  # a Bernoulli and a continuous prior follow the mixed-measure rule
  mixed <- .plot_data_prior_list.simple(
    list(bernoulli, prior("normal", list(0, 1)), prior("spike", list(1))),
    NULL, NULL, NULL, 101, 10000, FALSE, FALSE, NULL, NULL, FALSE)
  expect_named(mixed, c("density", "points1", "points2"))
  masses <- stats::setNames(c(mixed$points1$y, mixed$points2$y), c(mixed$points1$x, mixed$points2$x))
  expect_equal(masses[["0"]], .7 / 3)
  expect_equal(masses[["1"]], (.3 + 1) / 3)
  expect_equal(mixed$density$y, stats::dnorm(mixed$density$x) / 3)

  # the overlays draw the same atoms
  captured <- list()
  testthat::local_mocked_bindings(
    .lines.prior.point = function(plot_data, ...) captured$lines <<- c(captured$lines, list(plot_data))
  )
  plot(0, xlim = c(0, 1), ylim = c(0, 1))
  lines_prior_list(list(bernoulli))
  expect_equal(vapply(captured$lines, function(point) point$x, numeric(1)), c(0, 1))
  expect_equal(vapply(captured$lines, function(point) point$y, numeric(1)), c(.7, .3))
  overlay <- ggplot2::ggplot() + geom_prior_list(list(bernoulli))
  overlay_geoms <- unname(vapply(overlay$layers, function(layer) class(layer$geom)[1], character(1)))
  expect_equal(overlay_geoms, c("GeomSegment", "GeomSegment"))
})

test_that("prior-list plots refuse ordered priors and their mixtures clearly", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  general <- "Prior-list plots do not support ordered priors (or mixtures of ordered priors)"
  single  <- paste0(general, "; plot a single bound ordered prior with plot().")
  levels  <- c("lo", "mid", "hi")
  bound   <- prior_factor_levels(prior_ordered(prior("normal", list(0, 1))), levels)
  other   <- prior_factor_levels(prior_ordered(prior("normal", list(0, 2))), levels)
  unbound <- prior_ordered(prior("normal", list(0, 1)))

  # a single ordered prior points to plot()
  for(ordered in list(bound, unbound)){
    expect_error(plot_prior_list(list(ordered)), single, fixed = TRUE)
    expect_error(plot_prior_list(list(ordered), plot_type = "ggplot"), single, fixed = TRUE)
    expect_error(lines_prior_list(list(ordered)), single, fixed = TRUE)
    expect_error(geom_prior_list(list(ordered)), single, fixed = TRUE)
  }

  # plot() of a mixture of ordered priors goes through the list plot: no pointer back to it
  mixture <- prior_mixture(list(bound, other))
  expect_error(plot(mixture), paste0(general, "."), fixed = TRUE)
  expect_error(plot(mixture, plot_type = "ggplot"), paste0(general, "."), fixed = TRUE)

  # an ordered prior beside other priors
  expect_error(plot_prior_list(list(prior("normal", list(0, 1)), bound)), paste0(general, "."), fixed = TRUE)
  expect_error(plot_prior_list(list(bound, other)), paste0(general, "."), fixed = TRUE)

  # a plain list of ordinary priors is unaffected
  expect_no_error(plot_prior_list(list(prior("normal", list(0, 1)), prior("normal", list(1, 1))), plot_type = "ggplot"))
})

test_that("plot_prior_list rejects individual weightfunction lists", {

  weight <- prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .25)))
  expect_error(
    plot_prior_list(list(weight), individual = TRUE),
    "individual = TRUE"
  )
})

test_that("lines.prior reuses the active probability scale", {

  point_scale <- NULL
  testthat::local_mocked_bindings(
    .plot_scale_y2_state_current = function(){
      list(scale_y2 = 9, ylim2 = c(0, 1), usr = c(0, 1, 0, 9))
    },
    .lines.prior.point = function(plot_data, scale_y2 = 1, ...){

      point_scale <<- scale_y2
      return(invisible(NULL))
    },
    .package = "BayesTools"
  )

  lines(prior("point", list(location = 0)))
  expect_equal(point_scale, 9)
})

test_that("spike-and-slab prior plots use the inclusion probability", {

  x <- seq(-2, 2, length.out = 9)
  p <- prior_spike_and_slab(
    prior("normal", list(mean = 0, sd = .5)),
    prior_inclusion = prior("spike", list(.3))
  )
  expected_slab <- .3 * stats::dnorm(x, 0, .5)

  # ggplot plot(): slab .3 * N(0, .5) (peak .239) and spike mass .7
  g <- plot(p, plot_type = "ggplot", x_seq = x)
  layers <- prior_plot_layer_data(g)
  expect_equal(layers[[1]]$x, x)
  expect_equal(layers[[1]]$y, expected_slab, tolerance = 1e-12)
  expect_equal(layers[[2]]$x, 0)
  expect_equal(layers[[2]]$yend / attr(g, "scale_y2"), .7, tolerance = 1e-12)

  # geom_prior() on its own probability scale
  geom_layers <- prior_plot_layer_data(ggplot2::ggplot() + geom_prior(p, x_seq = x, scale_y2 = 2))
  expect_equal(geom_layers[[1]]$y, expected_slab, tolerance = 1e-12)
  expect_equal(geom_layers[[2]]$yend / 2, .7, tolerance = 1e-12)

  # base plot() and lines() agree
  drawn <- list()
  testthat::local_mocked_bindings(
    .lines.prior.simple = function(plot_data, ...){
      drawn$slab <<- c(drawn$slab, max(plot_data$y))
      invisible(NULL)
    },
    .lines.prior.point = function(plot_data, scale_y2 = 1, ...){
      drawn$spike <<- c(drawn$spike, plot_data$y)
      invisible(NULL)
    },
    .package = "BayesTools"
  )
  device_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(device_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(p, x_seq = x)
  lines(p, x_seq = x)
  expect_equal(drawn$slab, rep(.3 * stats::dnorm(0, 0, .5), 2), tolerance = 1e-12)
  expect_equal(drawn$spike, rep(.7, 2), tolerance = 1e-12)
})

test_that("geom_prior spike-and-slab xlim includes the spike at zero", {

  p <- prior_spike_and_slab(
    prior("normal", list(mean = 3, sd = .2), truncation = list(2, 4))
  )
  g <- ggplot2::ggplot() + geom_prior(p)
  layers <- ggplot2::ggplot_build(g)$data

  expect_true(any(layers[[2]]$x == 0))
  expect_lte(min(layers[[1]]$x), 0)
})
