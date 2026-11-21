skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

.month_level_fit <- function(month_prior, shift = 0){

  months <- c("May", "Jun", "Jul")
  data <- data.frame(Month = factor(rep(months, 2L), levels = months))
  info <- JAGS_formula(~Month, "mu", data, list(intercept = prior("normal", list(0, 1)),
    Month = month_prior))
  columns <- .JAGS_prior_factor_names("mu_Month", info$prior_list$mu_Month)
  x <- seq(-1, 1, length.out = 120L)
  coefficients <- cbind(3 * x + shift, cos(x * 6) + 2 * shift)
  if(is.prior.point(month_prior)) coefficients[] <- month_prior$parameters$location
  colnames(coefficients) <- columns
  draws <- cbind(mu_intercept = x, coefficients)
  if(is.prior.spike_and_slab(month_prior) || is.prior.mixture(month_prior)){
    included <- seq_along(x) %% 4L != 0L
    draws[!included, columns] <- 0
    indicator <- if(is.prior.spike_and_slab(month_prior)) as.integer(included) else ifelse(included, 1L, 2L)
    draws <- cbind(draws, mu_Month_indicator = indicator)
  }
  if(is.prior.regularized_term(info$prior_list$mu_Month)){
    draws <- cbind(draws, mu_Month_variance_1 = 1 + x^2)
  }
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)), sample = nrow(draws)),
    class = c("runjags", "BayesTools_fit", "list"))
  attr(fit, "prior_list") <- info$prior_list
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  attach_test_parameter_map(fit)
}

test_that("a selected centered level retains its curve, prior, atoms and colour", {
  for(contrast in c("meandif", "orthonormal")){
    slab <- prior_factor("mnormal", list(0, .5), contrast = contrast)
    fit <- .month_level_fit(prior_spike_and_slab(slab, prior("point", list(.65))))
    samples <- as_mixed_posteriors(fit, "mu_Month", n_prior_samples = 128L)
    expect_identical(posterior_factor_columns(samples, "mu_Month", c("Jul", "May")),
      c("mu_Month[dif: Jul]", "mu_Month[dif: May]"))
    draw <- function(...){
      plot_posterior(samples, "mu_Month", plot_type = "ggplot", prior = TRUE,
        n_points = 64L, xlim = c(-6, 6), ylim = c(0, 2), ylim2 = c(0, 1),
        legend_labels = c("Spring", "June", "July"), ...)
    }
    full <- draw()
    july <- draw(level = "Jul")
    expect_s3_class(july, "ggplot")
    layers <- ggplot2::ggplot_build(full)$data
    selected <- ggplot2::ggplot_build(july)$data
    # One common slab prior, one prior spike and one posterior spike are shared
    # across levels. The full plot's last layer contains all three posterior curves.
    expect_length(layers, 4L)
    expect_length(selected, 4L)
    expect_equal(selected[[1L]][c("x", "y")], layers[[1L]][c("x", "y")])
    for(i in 2:3){
      expected_mass <- c(.35, .25)[i - 1L]
      expect_equal(unique(selected[[i]]$x), 0)
      expect_equal(unique(selected[[i]]$yend) / july$bt_scale_y2_state$scale_y2, expected_mass)
      expect_equal(unique(layers[[i]]$yend) / full$bt_scale_y2_state$scale_y2, expected_mass)
    }
    colour_scale <- full$scales$get_scales("colour")
    july_colour <- unname(colour_scale$palette(3)[["July"]])
    july_rows <- layers[[4L]]$colour == july_colour
    expect_equal(selected[[4L]]$x, layers[[4L]]$x[july_rows])
    expect_equal(selected[[4L]]$y, layers[[4L]]$y[july_rows])
    expect_equal(unique(selected[[4L]]$colour), july_colour)
    expect_equal(unique(selected[[3L]]$colour), july_colour)
    subset_plot <- ggplot2::ggplot_build(draw(level = c("Jul", "Jun")))$data
    expect_equal(unique(subset_plot[[3L]]$colour), july_colour)
    reordered <- ggplot2::ggplot_build(draw(level = c("Jul", "Jun", "May")))$data
    expect_equal(unique(reordered[[3L]]$colour), july_colour)
    expect_identical(july$scales$get_scales("colour")$breaks, "July")
    # An overlay is built against the receiving ggplot's probability scale.
    overlay <- plot_posterior(samples, "mu_Month", level = "Jul", add = TRUE,
      prior = TRUE, plot_type = "ggplot", n_points = 64L, lty = 2)
    combined <- july + overlay
    overlay_layers <- utils::tail(ggplot2::ggplot_build(combined)$data, 4L)
    expect_identical(combined$bt_scale_y2_state, july$bt_scale_y2_state)
    expect_equal(unique(overlay_layers[[2L]]$yend) / july$bt_scale_y2_state$scale_y2, .35)
    expect_equal(unique(overlay_layers[[3L]]$yend) / july$bt_scale_y2_state$scale_y2, .25)

    figures <- draw(individual = TRUE)
    expect_named(figures, c("May", "Jun", "Jul"))
    expect_true(all(vapply(figures, inherits, logical(1), "ggplot")))
    expect_equal(ggplot2::ggplot_build(figures[["Jul"]])$data, selected)
    expect_named(draw(individual = TRUE, level = "Jul"), "Jul")
    expect_named(draw(individual = TRUE, level = c("Jul", "May")), c("Jul", "May"))
    expect_error(draw(level = "July"), "May.*Jun.*Jul")
  }
})

test_that("a transformed conditioned level retains the matching prior and posterior", {
  slab <- prior_factor("mnormal", list(0, .5), contrast = "meandif")
  fit <- .month_level_fit(prior_spike_and_slab(slab))
  conditional <- as_mixed_posteriors(fit, "mu_Month", conditional = "mu_Month", n_prior_samples = 128L)
  draw <- function(...){
    suppressMessages(plot_posterior(conditional, "mu_Month", transformation = "exp", prior = TRUE,
      plot_type = "ggplot", n_points = 64L, ...))
  }
  full <- ggplot2::ggplot_build(draw())$data
  july <- ggplot2::ggplot_build(draw(level = "Jul"))$data
  expect_length(july, 2L)
  # Conditioning on inclusion removes the atom. Each prior has a separate layer.
  expect_equal(july[[1L]]$x, full[[3L]]$x)
  expect_equal(july[[1L]]$y, full[[3L]]$y)
  july_rows <- full[[4L]]$colour == full[[3L]]$colour[1L]
  expect_equal(july[[2L]]$x, full[[4L]]$x[july_rows])
  expect_equal(july[[2L]]$y, full[[4L]]$y[july_rows])
  expect_true(all(july[[2L]]$x > 0))
})

test_that("base individual plots follow level order and overlays keep the active axes", {
  slab <- prior_factor("mnormal", list(0, .5), contrast = "meandif")
  p <- prior_spike_and_slab(slab)
  samples <- as_mixed_posteriors(.month_level_fit(p), "mu_Month", n_prior_samples = 128L)
  other <- as_mixed_posteriors(.month_level_fit(p, shift = 1), "mu_Month", n_prior_samples = 128L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  legends <- list()
  mass_tips <- numeric()
  spike_colours <- character()
  original_new <- graphics::plot.new
  pages <- 0L
  testthat::local_mocked_bindings(
    plot.new = function(){
      pages <<- pages + 1L
      original_new()
    },
    legend = function(x, legend, ...){
      legends[[length(legends) + 1L]] <<- list(labels = legend, col = list(...)$col)
    },
    arrows = function(x0, y0, x1 = x0, y1 = y0, ...){
      mass_tips <<- c(mass_tips, y1)
      spike_colours <<- c(spike_colours, list(...)$col)
    }, .package = "graphics")
  plot_posterior(samples, "mu_Month", individual = TRUE, n_points = 64L)
  expect_identical(pages, 3L)
  expect_error(plot_posterior(samples, "mu_Month", individual = TRUE, add = TRUE),
    "Separate factor figures require 'add = FALSE'", fixed = TRUE)
  expect_identical(pages, 3L)
  expect_equal(unlist(lapply(legends, "[[", "labels"), use.names = FALSE), c("May", "Jun", "Jul"))
  expect_equal(unname(spike_colours), unlist(lapply(legends, "[[", "col"), use.names = FALSE))
  plot_posterior(samples, "mu_Month", level = "Jul", prior = TRUE,
    n_points = 64L, ylim = c(0, 3), ylim2 = c(0, 1))
  expect_identical(pages, 4L)
  state <- .plot_scale_y2_state_current()
  usr <- graphics::par("usr")
  mass_tips <- numeric()
  plot_posterior(other, "mu_Month", level = "Jul", add = TRUE, lty = 2, n_points = 64L)
  expect_identical(pages, 4L)
  expect_identical(graphics::par("usr"), usr)
  expect_identical(.plot_scale_y2_state_current(), state)
  expect_equal(unique(mass_tips) / state$scale_y2, .25)
  expect_equal(legends[[5L]]$col, legends[[3L]]$col)
})

test_that("mixed factor overlays scale prior and posterior spikes to the receiving axes", {
  slab <- prior_factor("mnormal", list(0, .5), contrast = "meandif")
  spike_slab <- prior_spike_and_slab(slab, prior("point", list(.65)))
  mixed <- as_mixed_posteriors(.month_level_fit(spike_slab), "mu_Month", n_prior_samples = 128L)
  continuous <- as_mixed_posteriors(.month_level_fit(slab), "mu_Month", n_prior_samples = 128L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  tips <- numeric()
  testthat::local_mocked_bindings(
    arrows = function(x0, y0, x1 = x0, y1 = y0, ...) tips <<- c(tips, y1),
    .package = "graphics")
  for(level in list(NULL, "Jul")){
    draw <- function(samples, ...){
      plot_posterior(samples, "mu_Month", n_points = 64L, level = level, ...)
    }
    draw(mixed, prior = TRUE)
    standalone_scale <- .plot_scale_y2_state_current()$scale_y2
    expected <- c(.35, .25) * standalone_scale
    expect_equal(tips, expected)
    draw(continuous)
    expect_null(.plot_scale_y2_state_current())
    tips <- numeric()
    draw(mixed, prior = TRUE, add = TRUE)
    expect_equal(tips, expected)
    tips <- numeric()

    full <- draw(mixed, prior = TRUE, plot_type = "ggplot")
    first <- draw(continuous, plot_type = "ggplot")
    expect_null(first$bt_scale_y2_state)
    combined <- first + draw(mixed, prior = TRUE, add = TRUE, plot_type = "ggplot")
    layers <- ggplot2::ggplot_build(combined)$data
    spikes <- Filter(function(layer) "yend" %in% names(layer), layers)
    expect_equal(unname(vapply(spikes, function(spike) unique(spike$yend), numeric(1))),
      c(.35, .25) * full$bt_scale_y2_state$scale_y2)
  }
  # On a pure-point receiver, the primary y axis is already probability.
  point <- list(prior("point", list(0)))
  mixed_prior <- list(prior("normal", list(0, 1)), point[[1L]])
  for(receiver in list(
    function(type) plot(point[[1L]], xlim = c(-3, 3), plot_type = type),
    function(type) plot_prior_list(point, xlim = c(-3, 3), plot_type = type)
  )){
    receiver("base")
    expect_identical(.plot_scale_y2_state_current()$scale_y2, 1)
    tips <- numeric()
    plot_posterior(mixed, "mu_Month", level = "Jul", prior = TRUE, add = TRUE, n_points = 64L)
    plot_prior_list(mixed_prior, add = TRUE, n_points = 64L)
    lines_prior_list(mixed_prior, n_points = 64L)
    expect_equal(tips, c(.35, .25, .5, .5))
    combined <- receiver("ggplot") +
      plot_posterior(mixed, "mu_Month", level = "Jul", prior = TRUE, add = TRUE, plot_type = "ggplot", n_points = 64L) +
      plot_prior_list(mixed_prior, add = TRUE, plot_type = "ggplot", n_points = 64L) +
      geom_prior_list(mixed_prior, n_points = 64L)
    expect_identical(combined$bt_scale_y2_state$scale_y2, 1)
    spikes <- Filter(function(layer) "yend" %in% names(layer), ggplot2::ggplot_build(combined)$data)
    expect_equal(unname(vapply(spikes, function(spike) unique(spike$yend), numeric(1))),
      c(1, .35, .25, .5, .5))
  }
})

test_that("mixture, regularized, model-averaged and ordered effects use their fitted level labels", {
  slab <- prior_factor("mnormal", list(0, .5), contrast = "meandif")
  spike <- prior_factor("point", list(0), contrast = "meandif")
  mixture <- prior_mixture(list(slab, spike), is_null = c(FALSE, TRUE))
  regularized <- prior_regularized("lasso", .5, contrast = "meandif")
  check_july <- function(samples, parameter, prior = TRUE){
    figures <- plot_posterior(samples, parameter, plot_type = "ggplot", prior = prior,
      individual = TRUE, level = "Jul", n_points = 64L)
    expect_named(figures, "Jul")
    expect_s3_class(figures[["Jul"]], "ggplot")
    full <- plot_posterior(samples, parameter, plot_type = "ggplot", prior = prior, n_points = 64L)
    full_curve <- utils::tail(ggplot2::ggplot_build(full)$data, 1L)[[1L]]
    july_curve <- utils::tail(ggplot2::ggplot_build(figures[["Jul"]])$data, 1L)[[1L]]
    scale <- full$scales$get_scales("colour")
    july_colour <- scale$palette(length(scale$breaks))[["Jul"]]
    july_rows <- full_curve$colour == july_colour
    expect_equal(july_curve$x, full_curve$x[july_rows])
    expect_equal(july_curve$y, full_curve$y[july_rows])
  }
  for(p in list(mixture, regularized)){
    check_july(as_mixed_posteriors(.month_level_fit(p), "mu_Month", n_prior_samples = 128L), "mu_Month")
  }
  models <- lapply(c(0, 1), function(shift){
    list(fit = .month_level_fit(slab, shift), marglik = bridgesampling_object(0), prior_weights = 1)
  })
  mixed <- mix_posteriors(models, "mu_Month", list(c(FALSE, FALSE)), seed = 90, n_samples = 120L)
  check_july(mixed, "mu_Month", prior = FALSE)
  ordered <- ordered_plot_test_fixture(levels = c("May", "Jun", "Jul"))$samples
  check_july(ordered, "mu_f")
  expect_named(plot_posterior(ordered, "mu_f", plot_type = "ggplot", individual = TRUE,
    n_points = 64L), c("Jun", "Jul"))
  expect_error(plot_posterior(ordered, "mu_f", level = "May"), "Jun.*Jul")
})

test_that("a shared point mass identifies every selected level in the ggplot legend", {
  fit <- .month_level_fit(prior_factor("point", list(0), contrast = "treatment"))
  samples <- as_mixed_posteriors(fit, "mu_Month")
  for(levels in list(c("Jun", "Jul"), "Jul")){
    plot <- plot_posterior(samples, "mu_Month", plot_type = "ggplot", level = levels)
    built <- ggplot2::ggplot_build(plot)
    scale <- built$plot$scales$get_scales("colour")
    expect_identical(as.character(scale$get_breaks()), levels)
    expect_length(built$data, 1L)
    expect_equal(built$data[[1L]]$colour, scale$map(levels[[1L]]))
    expect_equal(built$data[[1L]]$yend, 1)
  }
})

test_that("model-averaged spikes shared by all levels use one fitted colour", {
  priors <- list(prior_factor("normal", list(0, 1), contrast = "treatment"),
    prior_factor("point", list(0), contrast = "treatment"),
    prior_factor("point", list(.5), contrast = "treatment"))
  models <- lapply(priors, function(p){
    list(fit = .month_level_fit(p), marglik = bridgesampling_object(0), prior_weights = 1)
  })
  samples <- mix_posteriors(models, "mu_Month", list(c(FALSE, TRUE, TRUE)), seed = 90, n_samples = 300L)
  for(levels in list(c("Jun", "Jul"), c("Jul", "Jun"))){
    plot <- plot_posterior(samples, "mu_Month", level = levels, n_points = 32L, plot_type = "ggplot")
    built <- ggplot2::ggplot_build(plot)
    spikes <- Filter(function(layer) "yend" %in% names(layer), built$data)
    expect_equal(vapply(spikes, "[[", numeric(1), "x"), c(0, .5))
    expected <- built$plot$scales$get_scales("colour")$map(levels[[1L]])
    expect_equal(vapply(spikes, "[[", character(1), "colour"), rep(expected, 2L))
  }
})

test_that("selected levels keep their fitted colours beyond the first eight levels", {
  samples <- ordered_plot_test_fixture(levels = month.abb)$samples
  draw <- function(...){
    plot_posterior(samples, "mu_f", plot_type = "ggplot", n_points = 32L, ...)
  }
  full <- draw()$scales$get_scales("colour")
  palette <- full$palette(length(full$breaks))
  expect_length(unique(palette), 11L)
  for(colors in list(NULL, c("red", "blue"), grDevices::rainbow(12L))){
    combined <- draw(col = colors)$scales$get_scales("colour")
    december <- draw(col = colors, level = "Dec")
    expected <- combined$palette(length(combined$breaks))[["Dec"]]
    expect_identical(unique(ggplot2::ggplot_build(december)$data[[1L]]$colour), expected)
  }
})

test_that("factor-valued random-slope SDs select the fitted treatment levels", {
  months <- c("May", "Jun", "Jul")
  data <- data.frame(f = factor(rep(months, 4L), levels = months),
    g = factor(rep(LETTERS[1:4], each = 3L)))
  spec <- JAGS_formula(~1 + diag(0 + f | g), "mu", data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(g = random_block(sd = prior("gamma", list(2, 2)))))
  term <- spec$formula_design$random_effects[[1L]]
  draws <- cbind(seq(0, 1, length.out = 100), seq(.1, .5, length.out = 100),
    seq(.8, 1.6, length.out = 100))
  colnames(draws) <- c("mu_intercept", term$sd_parameter_names)
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(draws)), sample = nrow(draws)),
    class = c("runjags", "BayesTools_fit", "list"))
  attr(fit, "prior_list") <- spec$prior_list
  attr(fit, "formula_design") <- list(mu = spec$formula_design)
  parameter <- setdiff(names(spec$prior_list), "mu_intercept")
  samples <- as_mixed_posteriors(attach_test_parameter_map(fit), parameter)
  draw <- function(...){
    plot_posterior(samples, parameter, plot_type = "ggplot", n_points = 64L, ...)
  }
  full_plot <- draw()
  full <- ggplot2::ggplot_build(full_plot)$data[[1L]]
  selected <- ggplot2::ggplot_build(draw(level = "Jul"))$data[[1L]]
  expect_named(draw(individual = TRUE), months[-1L])
  expect_error(draw(level = "May"), "Jun.*Jul")
  july_colour <- full_plot$scales$get_scales("colour")$palette(2)[["Jul"]]
  expect_equal(unique(selected$colour), july_colour)
  rows <- full$colour == july_colour
  expect_equal(selected$x, full$x[rows])
  expect_equal(selected$y, full$y[rows])
})

test_that("factor marginal plots select the original cells and their stored curves", {
  fit <- .month_level_fit(prior_factor("mnormal", list(0, .5), contrast = "orthonormal"))
  samples <- as_mixed_posteriors(fit, "mu_Month", n_prior_samples = 128L)
  marginal <- marginal_posterior(samples, "mu_Month", use_formula = FALSE,
    prior_samples = TRUE, n_samples = 128L)
  draw <- function(...){
    plot_marginal(list(months = marginal), "months", plot_type = "ggplot", prior = TRUE,
      n_points = 64L, legend_labels = c("Spring", "June", "July"), ...)
  }
  full <- ggplot2::ggplot_build(draw())$data
  figures <- draw(individual = TRUE, level = c("Jul", "May"))
  expect_named(figures, c("Jul", "May"))
  july <- ggplot2::ggplot_build(figures[["Jul"]])$data
  expect_length(july, 2L)
  expect_equal(july[[1L]]$x, full[[3L]]$x)
  expect_equal(july[[1L]]$y, full[[3L]]$y)
  rows <- full[[4L]]$colour == july[[2L]]$colour[[1L]]
  expect_equal(july[[2L]]$x, full[[4L]]$x[rows])
  expect_equal(july[[2L]]$y, full[[4L]]$y[rows])
  expect_identical(figures[["Jul"]]$scales$get_scales("colour")$breaks, "July")
  expect_error(draw(level = "July"), "May.*Jun.*Jul")
  expect_error(draw(individual = TRUE, add = TRUE), "Separate factor figures require 'add = FALSE'", fixed = TRUE)

  # Only July has a stored curve. Other levels retain their own KDEs.
  grid <- seq(-4, 4, length.out = 37)
  stored <- posterior_density_attribute(grid, stats::dnorm(grid), "iwmde", "precomputed", parameter = "Jul")
  marginal[["Jul"]] <- .bt_meta_set(marginal[["Jul"]], "posterior_density", stored)
  selected <- ggplot2::ggplot_build(draw(level = "Jul", density_method = "precomputed"))$data
  expect_equal(selected[[2L]]$x, grid)
  expect_equal(selected[[2L]]$y, stats::dnorm(grid))
  expect_error(draw(level = "May", density_method = "precomputed"), "Precomputed posterior density is unavailable")

  ordered <- ordered_plot_test_fixture(levels = c("May", "Jun", "Jul"))
  om <- marginal_posterior(ordered$samples, "mu_f", use_formula = FALSE)
  expect_named(plot_marginal(list(months = om), "months", plot_type = "ggplot",
    individual = TRUE, n_points = 64L), c("Jun", "Jul"))
  expect_error(plot_marginal(list(months = om), "months", level = "May"), "Jun.*Jul")
})
