skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

test_that("ordered direct displays retain both weighted non-reference priors", {
  fixture <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))
  p <- fixture$prior
  x <- c(-.7, -.2, -.05, 0, .05, .2, .7)
  original <- density(p, x_seq = x, n_points = length(x))
  physical <- vapply(x[x != 0], function(value){
    integral <- stats::integrate(function(share){
      stats::dnorm(value / share, 0, .5) / share
    }, 0, 1, rel.tol = 1e-10, abs.tol = 1e-12)
    expect_lt(integral$abs.error, 1e-10)
    integral$value
  }, numeric(1))
  expect_equal_each(original[[2L]]$y[x != 0], physical, tolerance = 2e-10)
  expect_identical(original[[2L]]$y[x == 0], Inf)
  expect_equal_each(original[[3L]]$y, stats::dnorm(x, 0, .5), tolerance = 1e-14)
  expect_equal(stats::integrate(function(share) .25 * share^2, 0, 1)$value, .25 / 3)
  columns <- .JAGS_prior_factor_names("mu_f", p)
  context <- posterior_metadata(fixture$samples, "prior_context")
  partial <- .prior_density_from_context(context, stats::setNames(c(1, 0), columns))
  complete <- .prior_density_from_context(context, stats::setNames(c(1, 1), columns))
  expect_identical(prior_density_ordinate(partial, 0)$behavior, "infinite")
  expect_identical(prior_density_ordinate(partial, 0)$point_mass, 0)
  expect_true(prior_density_ordinate(partial, 0)$exact)
  expect_identical(prior_density_ordinate(complete, 0)$behavior, "regular")
  displayed <- .plot_data_ordered_prior_display(p, original)
  expect_identical(displayed[[2L]]$x, x[x != 0])
  expect_identical(displayed[[2L]]$y, original[[2L]]$y[x != 0])
  expect_identical(displayed[[3L]], original[[3L]])
  expect_identical(original, density(p, x_seq = x, n_points = length(x)))
  plots <- plot(p, plot_type = "ggplot", x_seq = x, n_points = length(x))
  expect_length(plots, 3L)
  expect_null(plots[[1L]])
  expect_true(all(vapply(plots[2:3], inherits, logical(1), "ggplot")))
  expect_equal(ggplot2::ggplot_build(plots[[2L]])$data[[1L]]$y, physical, tolerance = 2e-10)
  expect_s3_class(plot(p, show_figures = 2L, plot_type = "ggplot", x_seq = x), "ggplot")
  expect_s3_class(plot(p, show_figures = 1L, plot_type = "ggplot", x_seq = x), "ggplot")
  selected <- plot(p, show_figures = -1L, plot_type = "ggplot", x_seq = x)
  expect_length(selected, 3L)
  expect_null(selected[[1L]])
  expect_true(all(vapply(selected[2:3], inherits, logical(1), "ggplot")))
  expect_length(geom_prior(p, x_seq = x), 2L)
  expect_length(geom_prior(p, show_parameter = 2L, x_seq = x), 1L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(p, x_seq = x))
  expect_no_warning(lines(p, x_seq = x))
})

test_that("visible ordered priors preserve genuine atoms and honest density failures", {
  for(total in list(prior("gamma", list(.5, 1)), prior("point", list(2)),
    prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))))){
    p <- ordered_plot_test_fixture(total, prior("dirichlet", list(c(2, 2, .5))))$prior
    original <- density(p, n_points = 64L)
    displayed <- .plot_data_ordered_prior_display(p, original)
    expect_true(any(displayed[[3L]]$y > 0))
    if(!is.null(original[[3L]]$atoms)){
      expect_identical(displayed[[3L]]$atoms, original[[3L]]$atoms)
      expect_identical(displayed[[3L]]$continuous, original[[3L]]$continuous)
    }
    expect_s3_class(plot(p, show_figures = 3L, plot_type = "ggplot", n_points = 64L), "ggplot")
  }
  p <- ordered_plot_test_fixture(prior("normal", list(0, 1)),
    prior("dirichlet", list(c(.5, .25, .25))))$prior
  expect_error(density(p, x_range = c(-3, 3), n_points = 64L),
    "Ordered-prior product density quadrature failed")
  expect_error(plot(p, xlim = c(-3, 3), n_points = 64L),
    "Ordered-prior product density quadrature failed")
  set.seed(600)
  original <- density(p, force_samples = TRUE, n_points = 64L, n_samples = 128L)
  original_rng <- .Random.seed
  set.seed(600)
  displayed <- .plot_data_ordered_prior_display(p,
    density(p, force_samples = TRUE, n_points = 64L, n_samples = 128L))
  expect_identical(displayed, original)
  expect_identical(.Random.seed, original_rng)
})

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

test_that("ordered points-only family guard certifies literal mixtures and preserves discrete sampling", {
  total <- prior_mixture(list(prior("point", list(0), prior_weights = .3),
    prior("point", list(2), prior_weights = .7)))
  for(allocation in list(c(.2, .3, .5), prior("dirichlet", list(c(2, 2, .5))))){
    p <- ordered_plot_test_fixture(total, allocation)$prior
    set.seed(600)
    initial_rng <- .Random.seed
    analytic <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L)
    expect_identical(.Random.seed, initial_rng)
    expect_identical(attr(analytic, "method"), "analytic_mixed_measure")
    set.seed(600)
    forced <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = TRUE)
    forced_rng <- .Random.seed
    set.seed(600)
    samples <- rng(p, 128L, transform_factor_samples = TRUE)
    expect_identical(.Random.seed, forced_rng)
    random <- is.prior(allocation)
    for(i in seq_len(4L)){
      partial <- random && i %in% 2:3
      share <- if(i == 1L) 0 else if(random) 1 else sum(allocation[seq_len(i - 1L)])
      expected_atoms <- if(partial) data.frame(location = 0, mass = .3) else if(share == 0)
        data.frame(location = 0, mass = 1) else data.frame(location = c(0, 2 * share), mass = c(.3, .7))
      expect_equal(analytic[[i]]$atoms, expected_atoms, tolerance = 2e-12)
      expect_equal(sum(analytic[[i]]$atoms$mass) + analytic[[i]]$diagnostics$continuous_mass, 1, tolerance = 2e-12)
      expect_null(analytic[[i]]$samples)
      expect_identical(forced[[i]]$samples, samples[, i])
      forced[[i]]["samples"] <- list(NULL)
      expect_identical(forced[[i]], analytic[[i]])
      if(partial){
        curve <- analytic[[i]]$continuous
        shapes <- if(i == 2L) c(2, 2.5) else c(4, .5)
        expect_equal_each(curve$density, .7 * stats::dbeta(curve$x / 2, shapes[1L], shapes[2L]) / 2, tolerance = 2e-12)
      }
    }
  }
  mixed_producer <- .density.prior.ordered_mixed
  unsupported_totals <- list(prior("bernoulli", list(.5)),
    prior("bernoulli", list(.5), truncation = list(lower = .5)),
    prior_mixture(list(prior("bernoulli", list(.5), prior_weights = .7),
      prior("point", list(2), prior_weights = .3))))
  for(total in unsupported_totals){
    for(allocation in list(c(.2, .3, .5), prior("dirichlet", list(c(2, 2, .5))))){
      p <- ordered_plot_test_fixture(total, allocation)$prior
      for(force in c(FALSE, TRUE)){
        set.seed(600)
        actual <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = force)
        actual_rng <- .Random.seed
        expect_null(attr(actual, "method", exact = TRUE))
        expect_identical(unname(vapply(actual, function(component) length(component$samples), integer(1))), rep(128L, 4L))
        local_mocked_bindings(.density.prior.ordered_mixed = function(...) NULL)
        set.seed(600)
        fallback <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = force)
        expect_identical(actual, fallback)
        expect_identical(.Random.seed, actual_rng)
        local_mocked_bindings(.density.prior.ordered_mixed = mixed_producer)
      }
    }
  }
})

test_that("ordered mixed-measure labels describe the actual displayed measure", {
  allocation <- prior("dirichlet", list(c(2, 2, 2)))
  point <- ordered_plot_test_fixture(prior("point", list(2)), allocation)$prior
  slab <- ordered_plot_test_fixture(prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))), allocation)$prior
  expect_identical(plot(point, show_figures = 4L, plot_type = "ggplot", n_points = 64L)$scales$get_scales("y")$name, "Probability")
  expect_identical(plot(point, show_figures = 2L, plot_type = "ggplot", n_points = 64L)$scales$get_scales("y")$name, "Density")
  expect_identical(plot(slab, show_figures = 2L, plot_type = "ggplot", n_points = 64L)$scales$get_scales("y")$name, "Density / probability mass")
  expect_identical(plot(point, show_figures = 4L, plot_type = "ggplot", n_points = 64L, ylab = "Custom")$scales$get_scales("y")$name, "Custom")
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
