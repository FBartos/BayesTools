skip_if_not_test_profile("unit")

test_that("horseshoe ordinate evaluation defers the coefficient density grid", {
  group <- prior_regularized("horseshoe", 0.35, slab_scale = 1)
  priors <- JAGS_formula(~ wt, "mu", mtcars,
    list(intercept = prior("point", list(0)), pool = group))$prior_list
  linear <- .prior_linear_combination_density(priors, c(mu_wt = 1), n_grid = 32L)
  expect_true(.prior_linear_density_is_deferred(linear))
  expect_identical(prior_density_ordinate(linear, 0)$behavior, "infinite")
  expect_equal(prior_density_ordinate(linear, 0.4)$log_density, lpdf(priors$mu_wt, 0.4))
  positive <- list(intervals = matrix(c(0, Inf), 1L), indicator = function(x) x > 0)
  expect_equal(as.numeric(.prior_linear_density_region_probability(linear, positive)), 0.5)
  interval <- list(intervals = matrix(c(-0.6, -0.2), 1L),
                   indicator = function(x) x > -0.6 & x < -0.2)
  expect_equal(as.numeric(.prior_linear_density_region_probability(linear, interval)),
               cdf(priors$mu_wt, -0.2) - cdf(priors$mu_wt, -0.6))
  # The deferred structural curve retains the plotting range of the eager grid.
  eager <- .prior_linear_combination_density(priors, c(mu_wt = 1), n_grid = 32L, .defer = FALSE)
  expect_identical(.prior_linear_density_to_plot_data(linear), .prior_linear_density_to_plot_data(eager))

  # The structural pole also exists when the global prior's tails make the
  # display quantiles numerically unavailable.
  heavy <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)),
    pool = prior_regularized("horseshoe", prior("invgamma", list(.25, .2)))))$prior_list
  deferred <- .prior_linear_combination_density(heavy, c(mu_wt = 1), n_grid = 32L)
  expect_identical(prior_density_ordinate(deferred, 0)$behavior, "infinite")
  expect_equal(as.numeric(.prior_linear_density_region_probability(deferred, positive)), 0.5)
})
test_that("automatic shrinkage plot limits do not build the coefficient grid", {
  pool <- prior_regularized("horseshoe", prior("cauchy", list(0, .2), list(0, Inf)))
  priors <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(2)), sh = pool))$prior_list
  weights <- c(mu_intercept = 1, mu_wt = -.5)
  tail_prob <- .002
  context <- .prior_density_context(priors, names(priors), tail_prob = tail_prob)
  testthat::local_mocked_bindings(
    .prior_linear_combination_grid = function(...) stop("A scalar plot tried to materialize its full prior grid."))
  law <- .prior_density_from_context(context, weights)
  limits <- sort(2 - .5 * mquant(priors$mu_wt, c(tail_prob, 1 - tail_prob)))
  affine <- .prior_linear_combination_density(priors, weights, tail_prob = tail_prob,
    output_transformation = "lin", output_transformation_arguments = list(a = -1, b = -2))
  curve <- .prior_linear_density_to_plot_data(law, n_points = 64)$density
  expect_length(curve$x, 64L)
  expect_equal(range(curve$x), limits)
  expect_equal(curve$y, mpdf(priors$mu_wt, (curve$x - 2) / -.5) / .5)
  log_point <- .prior_linear_combination_density(priors, weights, tail_prob = tail_prob,
    source_transforms = c(mu_intercept = "log", mu_wt = NA_character_))
  expect_equal(range(.prior_linear_density_to_plot_data(log_point, n_points = 20)$density$x),
               (limits - 2) + log(2))
  expect_equal(range(.prior_linear_density_to_plot_data(affine, n_points = 40)$density$x),
               sort(-1 - 2 * limits))
  displayed <- .prior_linear_density_to_plot_data(law, n_points = 40, x_range = c(2, 3),
    transformation = "lin", transformation_arguments = list(a = 3, b = -1), transformation_settings = TRUE)
  expect_equal(range(displayed$density$x), c(2, 3))
  gg <- ggplot2::ggplot_build(plot(law, plot_type = "ggplot", n_points = 64))
  expect_equal(range(gg$data[[1L]]$x), limits)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(law, n_points = 64))
  expect_invisible(plot(law, xlim = c(-.5, .5), n_points = 64))
})

test_that("regularization declarations distinguish fixed widths and learned scales", {
  expect_true(is.prior.simple(prior_regularized("lasso", 0.4)))
  expect_equal(pdf(prior_regularized("ridge", 0.4), 1), dnorm(1, sd = 0.4))
  hs <- prior_regularized("horseshoe", 0.5, terms = "wt", slab_scale = 2)
  expect_true(is.prior.regularized(hs))
  expect_match(print(hs, silent = TRUE), "Regularized horseshoe.*slab df = 4.*wt")
  expect_error(prior_regularized("horseshoe", prior("normal", list(0, 1))), "positive support")
  expect_error(prior_regularized("ridge", 1, terms = "intercept"), "non-intercept")
  expect_error(prior_regularized("ridge", Inf), "finite")
  expect_error(prior_regularized("lasso", 1, slab_scale = 1), "only to the horseshoe")
  expect_error(prior_regularized("horseshoe", 1, slab_df = 7), "Specify 'slab_scale'")
  expect_error(rng(hs, 3), "named 'prior_list'")
})

test_that("regularized region probabilities integrate small intervals through the pole", {
  hc <- prior("cauchy", list(0, .2), list(0, Inf))
  law <- function(group){
    priors <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)),
      shrink = group))$prior_list
    .prior_linear_combination_density(priors, c(mu_wt = 1), n_grid = 32L)
  }
  log_mass <- function(d, lo, hi){
    region <- list(intervals = matrix(c(lo, hi), 1L), indicator = function(x) lo < x & x < hi)
    attr(.prior_linear_density_region_probability(d, region), "log_mass")
  }
  expect_log_mass <- function(actual, expected, absolute_tolerance){
    expect_lt(abs(actual - expected), absolute_tolerance)
  }
  # With a fixed global scale s and half-Cauchy local, the density is
  # exp(u)*E1(u)/(s*sqrt(2*pi^3)), u=x²/(2s²). At |x|<=2e-12 the omitted
  # terms O(u*log(u)) are below 1e-20. A ridge with a half-Cauchy width
  # has exactly the same one-coefficient law.
  central <- function(x, s){
    x * (2 - 0.5772156649015329 - 2 * log(x) + log(2 * s^2)) / (s * sqrt(2 * pi^3))
  }
  a <- 1e-12
  b <- 2e-12
  for(s in c(.35, .2)){
    d <- law(if(s == .35) prior_regularized("horseshoe", s) else prior_regularized("ridge", hc))
    expect_log_mass(log_mass(d, a, b), log(central(b, s) - central(a, s)), 1e-8)
    expect_equal(log_mass(d, -b, -a), log_mass(d, a, b))
    expect_log_mass(log_mass(d, 0, a), log(central(a, s)), 1e-6)
    expect_log_mass(log_mass(d, -a, b), log(central(a, s) + central(b, s)), 1e-6)
    expect_equal(log_mass(d, 0, Inf), -log(2))
  }

  # Independent reference integrals of the Normal mixture on log width.
  # Two independent half-Cauchys have log-product density 2v/(pi² sinh v);
  # for the finite slab, integrate also over chi-square(4) = 4/c².
  learned <- law(prior_regularized("horseshoe", hc))
  expect_log_mass(log_mass(learned, 1e-8, 2e-8), log(1.1066614270590066e-6), 1e-8)
  slab <- law(prior_regularized("horseshoe", .35, slab_scale = 1, slab_df = 4))
  expect_log_mass(log_mass(slab, 1e-8, 2e-8), log(1.2464104824557918e-7), 1e-8)
  expect_log_mass(log_mass(slab, 1e-3, 1.1e-3), log(4.3550463195e-4), 1e-7)
  # The independently computed mass from zero must compose consistently
  # with a finite interval on the same side of the integrable pole.
  for(d in list(learned, slab)){
    expect_log_mass(log_mass(d, 0, b), .model_probability_log_sum(c(log_mass(d, 0, a), log_mass(d, a, b))), 1e-6)
  }
})

test_that("learned scale marginal ordinates agree with independent integration", {
  bound <- function(family, scale, ...){
    p <- JAGS_formula(~ wt, "mu", mtcars, list(
      intercept = prior("point", list(0)), sh = prior_regularized(family, scale, ...)))
    p$prior_list[["mu_wt"]]
  }
  global <- prior("gamma", list(2, 1))
  ridge <- bound("ridge", global)
  lasso <- bound("lasso", global)
  for(value in c(0, 0.3, 2)){
    normal <- integrate(function(t) dnorm(value, sd = t) * dgamma(t, 2), 0, Inf)
    laplace <- integrate(function(t) exp(-value/t) / (2*t) * dgamma(t, 2), 0, Inf)
    expect_equal(pdf(ridge, value), normal$value, tolerance = 1e-6)
    expect_equal(pdf(lasso, value), laplace$value, tolerance = 1e-6)
  }
  expect_equal(quant(ridge, c(0, 0.5, 1)), c(-Inf, 0, Inf))
  expect_equal(cdf(lasso, quant(lasso, .82)), .82, tolerance = 1e-6)
  expect_equal(var(ridge), 6)
  expect_equal(var(lasso), 12)
  expect_equal(prior_density_ordinate(bound("ridge", prior("exp", list(1))), 0)$behavior,
               "infinite")
  hs <- bound("horseshoe", 0.35)
  expect_equal(pdf(hs, 0), Inf)
  expect_true(is.nan(mean(hs)))
  expect_true(is.nan(var(hs)))
  expected <- integrate(function(t) dnorm(1, sd = 0.35 * t) * 2 * dcauchy(t), 0, Inf)
  expect_equal(pdf(hs, 1), expected$value, tolerance = 1e-6)
  finite <- bound("horseshoe", 0.35, slab_scale = 1)
  expect_equal(prior_density_ordinate(finite, 0)$behavior, "infinite")
  # Integrate the local scale and the inverse-chi-square slab independently.
  expected <- integrate(function(u) vapply(u, function(v){
    integrate(function(t) pnorm(-0.8 * sqrt(1 / (0.35*t)^2 + v/4)) * 2 * dcauchy(t),
              0, Inf)$value * dchisq(v, 4)
  }, numeric(1)), 0, Inf)
  expect_equal(ccdf(finite, 0.8), expected$value, tolerance = 1e-5)
})

test_that("random-global finite-slab horseshoe agrees with the direct scale-mixture integral", {
  global <- prior("cauchy", list(0, 0.2), list(0, Inf))
  hs <- JAGS_formula(~ wt, "mu", mtcars, list(
    intercept = prior("point", list(0)),
    shrink = prior_regularized("horseshoe", global, slab_scale = 2, slab_df = 4)
  ))$prior_list$mu_wt

  # For independent standard half-Cauchys the density of the log product
  # z is 2*z/(pi^2*sinh(z)). Mix the Normal over this coordinate and the
  # independent chi-square(4) coordinate of c^2 = 16 / chi-square(4).
  reference <- function(x, density){
    integrand <- function(z){
      mass <- ifelse(z == 0, 2 / pi^2, 2 * abs(z) / (pi^2 * sinh(abs(z))))
      vapply(seq_along(z), function(i){
        if(mass[i] == 0) return(0)
        stats::integrate(function(v){
          sd <- 1 / sqrt(25 * exp(-2 * z[i]) + v / 16)
          conditional <- if(density) stats::dnorm(x, sd = sd) else stats::pnorm(-x, sd = sd)
          conditional * stats::dchisq(v, 4)
        }, 0, Inf, rel.tol = 1e-8)$value * mass[i]
      }, numeric(1))
    }
    cuts <- sort(unique(c(-Inf, log(x / .2), 0, Inf)))
    sum(vapply(seq_len(length(cuts) - 1L), function(i){
      stats::integrate(integrand, cuts[i], cuts[i + 1L], rel.tol = 1e-7)$value
    }, numeric(1)))
  }
  x <- c(1e-10, 0.001, 1.4)
  expect_equal(pdf(hs, x), vapply(x, reference, numeric(1), density = TRUE), tolerance = 2e-6)
  expect_equal(ccdf(hs, x), vapply(x, reference, numeric(1), density = FALSE), tolerance = 2e-6)
  # Integrating the two half-Cauchys gives c² E[tau/(tau+c)].
  variance <- stats::integrate(function(v){
    c2 <- 16 / v
    a <- sqrt(c2) / 0.2
    c2 * (1 + 2 * a * log(a) / pi) / (1 + a^2) * stats::dchisq(v, 4)
  }, 0, Inf, rel.tol = 1e-10)$value
  expect_equal(var(hs), variance, tolerance = 1e-6)
  expect_equal(ccdf(hs, quant(hs, 0.99)), .01, tolerance = 1e-6)
})

test_that("ridge and lasso resolve the peak under a learned half-Cauchy width", {
  # In u=log(tau/0.2) coordinates the Gaussian or Laplace kernel's 1/tau
  # cancels the Jacobian. Integrate on both sides of the kernel transition.
  x <- 1e-10
  cuts <- c(-Inf, log(x / .2), 0, Inf)
  for(family in c("ridge", "lasso")){
    model <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)),
      sh = prior_regularized(family, prior("cauchy", list(0, .2), list(0, Inf)))))
    integrand <- function(u){
      z <- exp(log(x / .2) - u)
      kernel <- if(family == "ridge") stats::dnorm(z) else exp(-z) / 2
      kernel * (2 / (.2 * pi)) * stats::plogis(-2 * u)
    }
    exact <- sum(vapply(seq_len(3), function(i){
      stats::integrate(integrand, cuts[i], cuts[i + 1L], rel.tol = 1e-10)$value
    }, numeric(1)))
    expect_equal(pdf(model$prior_list$mu_wt, x), exact, tolerance = 1e-7)
  }
})

test_that("regularized prior plots produce finite curves on their automatic ranges", {
  hc <- prior("cauchy", list(0, .2), list(0, Inf))
  groups <- list(prior_regularized("ridge", hc),
    prior_regularized("ridge", prior("gamma", list(2, 4))), prior_regularized("lasso", hc),
    prior_regularized("horseshoe", 0.35), prior_regularized("horseshoe", hc),
    prior_regularized("horseshoe", hc, slab_scale = 2))
  for(group in groups){
    model <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)), sh = group))
    curve <- ggplot2::ggplot_build(plot(model$prior_list$mu_wt, plot_type = "ggplot", n_points = 64))$data[[1L]]
    expect_true(all(is.finite(curve$x)) && sum(is.finite(curve$y) & curve$y > 0) >= 60)
  }
  # A far-tail bracket used to fail even for an ordinary fitted prior
  # range, preventing the structural ordinate from being returned.
  model <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)),
    sh = prior_regularized("horseshoe", hc)))
  hs <- model$prior_list$mu_wt
  expect_equal(ccdf(hs, quant(hs, 1 - 1e-5)), 1e-5, tolerance = 1e-6)
  expect_identical(prior_density_ordinate(hs, 0)$behavior, "infinite")
})
test_that("multi-column shrinkage terms distinguish joint and marginal methods", {
  group <- prior_regularized("ridge", prior("gamma", list(2, 1)))
  data <- data.frame(g = factor(rep(letters[1:3], 2)))
  term <- JAGS_formula(~ g, "mu", data, list(intercept = prior("point", list(0)),
                                            shrinkage = group))$prior_list$mu_g
  expect_error(pdf(term, c(0, 1)), class = "BayesTools_regularized_joint_unavailable")
  expect_equal(mpdf(term, 0), dnorm(0))
  expect_equal(dim(rng(term, 5, transform_factor_samples = FALSE)), c(5L, 2L))
  levels <- rng(term, 5)
  expect_equal(dim(levels), c(5L, 3L))
  expect_true(all(levels[, 1L] == 0))
})
