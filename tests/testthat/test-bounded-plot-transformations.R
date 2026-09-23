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
  tanh_map <- density(cauchy, transformation = "tanh")
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
  expect_no_error(plot(cauchy, transformation = "tanh"))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(cauchy, transformation = "tanh", plot_type = "ggplot"), "ggplot")
})

test_that("transformation output support metadata is validated", {

  transformation <- .strict_bounded_plot_transform()
  for(support in list(c(1, -1), c(0, 0), c(NA_real_, 1), 1, "bounded")){
    transformation$output_support <- support
    expect_error(density(prior("normal", list(0, 1)), transformation = transformation),
                 "output_support", fixed = TRUE)
  }
})
