skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

test_that("ordered points-only totals retain exact scalar measures and sampled metadata", {
  cases <- list(
    positive = list(location = 2, probability = 1, alpha = c(2, 2, .5)),
    slab = list(location = 2, probability = .5, alpha = c(2, 2, .5)),
    finite = list(location = 2, probability = .5, alpha = c(2, 2, 2)),
    negative = list(location = -2, probability = 1, alpha = c(2, 2, .5)),
    negative_slab = list(location = -2, probability = .5, alpha = c(2, 2, .5)),
    fixed = list(location = 2, probability = 1, allocation = c(.2, .3, .5)),
    fixed_zero = list(location = 2, probability = .5, allocation = c(0, .3, .7)))
  transformations <- list(identity = list(name = NULL, arguments = NULL),
    exp = list(name = "exp", arguments = NULL),
    tanh = list(name = "tanh", arguments = NULL),
    lin = list(name = "lin", arguments = list(a = 3, b = -2)))
  for(case in cases){
    total <- prior("point", list(case$location))
    if(case$probability < 1) total <- prior_spike_and_slab(total, prior("point", list(case$probability)))
    allocation <- if(is.null(case$alpha)) case$allocation else prior("dirichlet", list(case$alpha))
    p <- ordered_plot_test_fixture(total, allocation)$prior
    for(transformation in transformations){
      forward <- function(x){
        switch(if(is.null(transformation$name)) "identity" else transformation$name,
          identity = x, exp = exp(x), tanh = tanh(x), lin = 3 - 2 * x)
      }
      inverse <- function(x){
        switch(if(is.null(transformation$name)) "identity" else transformation$name,
          identity = x, exp = log(x), tanh = atanh(x), lin = (x - 3) / -2)
      }
      jacobian <- function(x){
        switch(if(is.null(transformation$name)) "identity" else transformation$name,
          identity = rep(1, length(x)), exp = 1 / x, tanh = 1 / (1 - x^2), lin = rep(.5, length(x)))
      }
      args <- list(x = p, n_points = 64L, n_samples = 128L,
        transformation = transformation$name, transformation_arguments = transformation$arguments)
      set.seed(600)
      analytic <- do.call(density, args)
      set.seed(600)
      forced <- do.call(density, c(args, list(force_samples = TRUE)))
      set.seed(600)
      samples <- rng(p, 128L, transform_factor_samples = TRUE)
      expect_identical(attr(analytic, "method"), "analytic_mixed_measure")
      for(i in seq_len(4L)){
        m <- i - 1L
        fractional <- !is.null(case$alpha) && m > 0L && m < 3L
        share <- if(m == 0L) 0 else if(is.null(case$alpha)) sum(case$allocation[seq_len(m)]) else 1
        expected_atoms <- if(fractional){
          if(case$probability == 1) data.frame(location = numeric(), mass = numeric()) else
            data.frame(location = forward(0), mass = 1 - case$probability)
        }else if(share == 0){
          data.frame(location = forward(0), mass = 1)
        }else if(case$probability == 1){
          data.frame(location = forward(case$location * share), mass = 1)
        }else{
          data.frame(location = forward(c(0, case$location * share)), mass = c(1 - case$probability, case$probability))
        }
        atoms <- analytic[[i]]$atoms
        expect_equal(unname(as.matrix(atoms[order(atoms$location), , drop = FALSE])),
          unname(as.matrix(expected_atoms[order(expected_atoms$location), , drop = FALSE])), tolerance = 2e-12)
        expect_null(analytic[[i]]$samples)
        expect_identical(forced[[i]]$samples, forward(samples[, i]))
        forced[[i]]["samples"] <- list(NULL)
        expect_identical(forced[[i]], analytic[[i]])
        expect_equal(analytic[[i]]$diagnostics$continuous_mass, if(fractional) case$probability else 0)
        if(fractional){
          a <- sum(case$alpha[seq_len(m)])
          b <- sum(case$alpha[(m + 1L):3L])
          curve <- analytic[[i]]$continuous
          source <- inverse(curve$x)
          expected <- case$probability * stats::dbeta(source / case$location, a, b) /
            abs(case$location) * jacobian(curve$x)
          expect_equal_each(curve$density, expected, tolerance = 2e-12)
          expect_equal(case$probability * diff(stats::pbeta(c(0, 1), a, b)),
            analytic[[i]]$diagnostics$continuous_mass, tolerance = 0)
          expect_equal(analytic[[i]]$diagnostics$continuous_integral,
            sum(diff(curve$x) * (head(curve$density, -1L) + tail(curve$density, -1L)) / 2), tolerance = 0)
          expect_true(all(is.finite(curve$x) & is.finite(curve$density)))
          if(is.null(transformation$name) && b >= 1){
            expect_equal(curve$density[c(1L, nrow(curve))],
              case$probability * stats::dbeta(if(case$location > 0) c(0, 1) else c(1, 0), a, b) /
                abs(case$location), tolerance = 2e-12)
          }
        }else{
          expect_null(analytic[[i]]$continuous)
        }
      }
    }
  }
})

test_that("ordered direct plotting omits certified curves before product quadrature", {
  fixture <- ordered_plot_test_fixture(prior("normal", list(0, 1)), prior("dirichlet", list(c(.5, .25, .25))))
  p <- fixture$prior
  original <- p
  set.seed(600)
  expect_error(density(p, x_range = c(-3, 3), n_points = 64L, n_samples = 128L),
    "Ordered-prior product density quadrature failed")
  calls <- 0L
  local_mocked_bindings(.density.prior.ordered_dirichlet_product = function(...){
    calls <<- calls + 1L
    stop("An omitted product curve reached quadrature.")
  })
  input <- .plot_ordered_prior_density_input(p)
  expect_identical(attr(input, "ordered_plot_skip", exact = TRUE), c(FALSE, TRUE, TRUE, FALSE))
  displayed <- density(input, x_range = c(-3, 3), n_points = 64L, n_samples = 128L)
  expect_identical(unname(vapply(displayed, attr, integer(1), "component")), seq_len(4L))
  expect_identical(unname(vapply(displayed, attr, character(1), "component_name")), names(displayed))
  expect_true(all(vapply(displayed[2:3], inherits, logical(1), "density.prior.display_empty")))
  expect_null(attr(displayed[[2L]], "x_range", exact = TRUE))
  expect_null(attr(displayed[[3L]], "y_range", exact = TRUE))
  full <- density(prior("normal", list(0, 1)), x_range = c(-3, 3), n_points = 64L)
  expect_identical(displayed[[4L]]$x, full$x)
  expect_identical(displayed[[4L]]$y, full$y)
  expect_identical(displayed[[4L]]$diagnostics, full$diagnostics)
  expect_equal(attr(displayed, "y_range", exact = TRUE), c(0, 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(600)
  expect_no_warning(plot(p, plot_type = "base", xlim = c(-3, 3), n_points = 64L, n_samples = 128L))
  set.seed(600)
  plots <- plot(p, plot_type = "ggplot", xlim = c(-3, 3), n_points = 64L, n_samples = 128L)
  expect_length(plots, 4L)
  expect_null(plots[[2L]])
  expect_null(plots[[3L]])
  expect_true(all(vapply(plots[c(1L, 4L)], inherits, logical(1), "ggplot")))
  expect_no_warning(lines(p, show_parameter = 2L, xlim = c(-3, 3), n_points = 64L, n_samples = 128L))
  expect_length(geom_prior(p, show_parameter = 3L, xlim = c(-3, 3), n_points = 64L, n_samples = 128L), 0L)
  expect_error(plot(p, show_figures = 2L, plot_type = "ggplot", xlim = c(-3, 3),
    n_points = 64L, n_samples = 128L), class = "BayesTools_ordered_prior_display_empty")
  expect_identical(calls, 0L)
  expect_identical(p, original)
  expect_null(attr(p, "ordered_plot_skip", exact = TRUE))
  expect_error(density(p, x_range = c(-3, 3), n_points = 64L, n_samples = 128L),
    "An omitted product curve reached quadrature", fixed = TRUE)
  expect_identical(calls, 1L)
  fixed <- ordered_plot_test_fixture(prior("normal", list(0, 1)), c(.2, .3, .5))$prior
  attr(fixed, "ordered_plot_skip") <- c(FALSE, TRUE, FALSE, FALSE)
  scaled_total <- .density.prior.ordered_scaled_total
  scales_seen <- numeric()
  local_mocked_bindings(.density.prior.ordered_scaled_total = function(...){
    args <- list(...)
    scales_seen <<- c(scales_seen, args$scale)
    scaled_total(...)
  })
  fixed_display <- density(fixed, x_range = c(-3, 3), n_points = 64L)
  expect_s3_class(fixed_display[[2L]], "density.prior.display_empty")
  expect_identical(scales_seen, c(0, .5, 1))
})

test_that("ordered plot omission preserves sampling and mixed atoms in every density path", {
  p <- ordered_plot_test_fixture(prior("normal", list(0, 1)), prior("dirichlet", list(c(.5, .25, .25))))$prior
  set.seed(600)
  original <- density(p, force_samples = TRUE, x_range = c(-3, 3), n_points = 64L, n_samples = 128L)
  original_rng <- .Random.seed
  kde <- .density_kde_boundary
  calls <- 0L
  local_mocked_bindings(.density_kde_boundary = function(...){calls <<- calls + 1L; kde(...)})
  set.seed(600)
  omitted <- density(.plot_ordered_prior_density_input(p), force_samples = TRUE,
    x_range = c(-3, 3), n_points = 64L, n_samples = 128L)
  expect_identical(.Random.seed, original_rng)
  expect_identical(omitted[c(1L, 4L)], original[c(1L, 4L)])
  expect_identical(calls, 1L)
  mixed <- ordered_plot_test_fixture(prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.5))),
    prior("dirichlet", list(c(.5, .25, .25))))$prior
  route_mixed <- .density.prior.ordered_route_mixed
  weights_seen <- list()
  local_mocked_bindings(.density.prior.ordered_route_mixed = function(...){
    args <- list(...)
    weights_seen[[length(weights_seen) + 1L]] <<- args$weights
    route_mixed(...)
  })
  for(force in c(FALSE, TRUE)){
    set.seed(600)
    displayed <- .plot_data_ordered_prior_display(mixed,
      density(.plot_ordered_prior_density_input(mixed), force_samples = force,
        x_range = c(-3, 3), n_points = 64L, n_samples = 128L))
    for(i in 2:3){
      expect_equal(displayed[[i]]$atoms, data.frame(location = 0, mass = .5))
      expect_null(displayed[[i]]$continuous)
      expect_false(inherits(displayed[[i]], "density.prior.display_empty"))
      expect_identical(attr(displayed[[i]], "component"), i)
    }
  }
  expect_length(weights_seen, 4L)
  expect_true(all(vapply(weights_seen, function(w) all(w == 0) || all(w == 1), logical(1))))
})

test_that("ordered private plot masks validate declared levels and retain unsuppressed controls", {
  p <- prior_ordered(prior("normal", list(0, 1)), allocation = prior("dirichlet", list(c(.5, .25, .25))))
  attr(p, "levels") <- 4L
  input <- .plot_ordered_prior_density_input(p)
  expect_identical(attr(input, "ordered_plot_skip", exact = TRUE), c(FALSE, TRUE, TRUE, FALSE))
  expect_null(attr(p, "ordered_metadata", exact = TRUE))
  expect_null(attr(p, "ordered_plot_skip", exact = TRUE))
  for(mask in list(c(FALSE, TRUE), c(FALSE, NA, TRUE, FALSE), c(0, 1, 1, 0))){
    attr(input, "ordered_plot_skip") <- mask
    expect_error(density(input, x_range = c(-3, 3), n_points = 64L),
      "The private ordered plot mask must contain one non-missing logical value per declared level.", fixed = TRUE)
  }
  controls <- list(
    intrinsic = ordered_plot_test_fixture(prior("gamma", list(.5, 1)))$prior,
    fixed = ordered_plot_test_fixture(prior("normal", list(0, 1)), c(.2, .3, .5))$prior,
    finite = ordered_plot_test_fixture(prior("point", list(2)), prior("dirichlet", list(c(2, 2, 2))))$prior)
  for(control in controls){
    expect_null(attr(.plot_ordered_prior_density_input(control), "ordered_plot_skip", exact = TRUE))
  }
  ordinary <- prior("dirichlet", list(c(.5, .25, .25)))
  expect_identical(.plot_ordered_prior_density_input(ordinary), ordinary)
  local_mocked_bindings(.prior_density_route_linear = function(...) list(type = "unknown"))
  expect_null(attr(.plot_ordered_prior_density_input(p), "ordered_plot_skip", exact = TRUE))
})

test_that("ordered complex points-only totals retain their existing sampled fallback", {
  for(two_ordered in c(FALSE, TRUE)){
    data <- expand.grid(f = ordered(c("early", "middle", "late", "last"),
      levels = c("early", "middle", "late", "last")),
      g = if(two_ordered) ordered(c("low", "mid", "high"), levels = c("low", "mid", "high")) else factor(c("a", "b", "c")))
    build <- function(total){
      JAGS_formula(~f*g, "mu", data, list(intercept = prior("point", list(0)),
        f = prior_ordered(prior("point", list(2))),
        g = if(two_ordered) prior_ordered(prior("point", list(2))) else
          prior_factor("normal", list(0, 1), contrast = "treatment"),
        "f:g" = prior_ordered(total)))$prior_list$mu_f__xXx__g
    }
    p <- build(prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))))
    for(force in c(FALSE, TRUE)){
      set.seed(600)
      sampled <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = force)
      expect_null(attr(sampled, "method", exact = TRUE))
      expect_identical(unname(vapply(sampled, function(component) length(component$samples), integer(1))),
        rep(128L, length(sampled)))
    }
    mixed <- build(prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.5))))
    expect_error(density(mixed, x_range = c(-3, 3), n_points = 64L, n_samples = 128L),
      paste0("Mixed-measure ordered densities currently require one ordered term and one scalar total. ",
        "Split the interaction into explicitly named terms before requesting its density."), fixed = TRUE)
  }
})
