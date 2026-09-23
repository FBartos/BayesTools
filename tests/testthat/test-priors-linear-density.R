skip_if_not_test_profile("unit")

test_that("bounded conditional-normal prior ordinates retain their numerical error", {

  priors <- list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1)),
                 s = prior("cauchy", list(0, 1), list(0, 5)))
  attr(priors$b, "multiply_by") <- "s"
  density <- .prior_linear_combination_density(priors, c(a = 1, b = -1), n_grid = 1024)
  expected <- 5 / (sqrt(26) * sqrt(2 * pi) * atan(5))
  ordinate <- prior_density_ordinate(density, 0)
  expect_identical(ordinate$behavior, "regular")
  expect_identical(ordinate$method, "conditional_normal_mixture")
  expect_true(ordinate$exact)
  expect_false(ordinate$provenance$integration$exact)
  expect_equal(exp(ordinate$log_density), expected, tolerance = 1e-9)
  integration <- ordinate$provenance$integration
  expect_lte(integration$absolute_error, integration$error_bound)
  expect_lte(integration$evaluations, 1024)
  expect_equal(as.numeric(.prior_linear_density_height(density, 0)), expected, tolerance = 1e-9)

  # atan(s) is uniform: an independent change-of-variable integral avoids
  # the multiplier-density routine and its original integration coordinates.
  for(value in c(-2, .75, 4)){
    reference <- stats::integrate(function(angle){
      stats::dnorm(value, sd = 1 / cos(angle)) / atan(5)
    }, 0, atan(5), rel.tol = 1e-12)$value
    result <- prior_density_ordinate(density, value)
    expect_lt(abs(exp(result$log_density) - reference),
              .prior_linear_density_refinement_tolerance()$absolute +
                .prior_linear_density_refinement_tolerance()$relative * reference)
  }

  shifted_priors <- list(a = prior("normal", list(.7, .5)),
                         b = prior("normal", list(.3, 1.2)),
                         s = prior("uniform", list(-2, 3)))
  attr(shifted_priors$b, "multiply_by") <- "s"
  shifted <- .prior_linear_combination_density(shifted_priors, c(a = 1, b = 1), n_grid = 512)
  shifted_reference <- stats::integrate(function(probability){
    s <- -2 + 5 * probability
    stats::dnorm(1.1, mean = .7 + .3 * s, sd = sqrt(.5^2 + 1.2^2 * s^2))
  }, 0, 1, rel.tol = 1e-12)$value
  expect_equal(as.numeric(.prior_linear_density_height(shifted, 1.1)),
               shifted_reference, tolerance = 1e-7)

  pure_product <- .prior_linear_combination_density(priors, c(b = 1), n_grid = 128)
  expect_identical(prior_density_ordinate(pure_product, 0)$behavior, "infinite")
  expect_identical(prior_density_ordinate(pure_product, 0)$method, "unsupported_provenance")
  split <- .prior_linear_split_multiply_groups(priors, c(a = 1, b = 1))
  split$additive_weights <- c(b = 1)
  expect_null(.prior_conditional_normal_spec(priors, split, c(a = NA_character_, b = NA_character_)))
})

test_that("unbounded conditional-normal ordinates retain exact support and counted budgets", {

  multipliers <- list(prior("normal", list(0, 1)),
                       prior("cauchy", list(0, 1), list(0, Inf)),
                       prior("cauchy", list(0, 1), list(-Inf, 0)))
  # s = sinh(t) reduces the normal-multiplier center to the K0 integral:
  # https://dlmf.nist.gov/10.32.E9
  expected <- c(besselK(.25, 0, expon.scaled = TRUE) / (2 * pi),
                 rep(sqrt(2) / pi^(3/2), 2))
  original_integrate <- stats::integrate
  counted <- 0L
  testthat::local_mocked_bindings(
    integrate = function(f, ...){
      original_integrate(function(x){
        counted <<- counted + length(x)
        f(x)
      }, ...)
    },
    .package = "stats"
  )
  for(i in seq_along(multipliers)){
    priors <- list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1)),
                   s = multipliers[[i]])
    attr(priors$b, "multiply_by") <- "s"
    density <- .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 512)
    counted <- 0L
    ordinate <- prior_density_ordinate(density, 0)
    expect_identical(ordinate$behavior, "regular")
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), expected[[i]], tolerance = 1e-7)
    expect_identical(ordinate$provenance$integration$evaluations, counted)
    expect_lte(counted, 512)
    expect_equal(as.numeric(.prior_linear_density_height(density, 0)), expected[[i]], tolerance = 1e-7)
  }
})

test_that("legacy product flags cannot establish an infinite density", {

  priors <- list(a = prior("uniform", list(-1, 1)), b = prior("normal", list(0, 1)),
                 s = prior("normal", list(0, 1)))
  attr(priors$b, "multiply_by") <- "s"
  unsupported <- .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 128)
  expect_true(0 %in% attr(unsupported, "singular_density_points"))
  expect_identical(prior_density_ordinate(unsupported, 0)$behavior, "unknown")
  expect_error(.prior_linear_density_height(unsupported, 0),
               "The prior density at the flagged product ordinate is unavailable from supported deterministic provenance.",
               fixed = TRUE)
  pure <- .prior_linear_combination_density(priors, c(b = 1), n_grid = 128)
  expect_identical(prior_density_ordinate(pure, 0)$behavior, "infinite")
  expect_identical(.prior_linear_density_height(pure, 0), Inf)
  attr(pure, "singular_density_points") <- NULL
  expect_identical(.prior_linear_density_height(pure, 0), Inf)
})

test_that("conditional-normal ordinates preserve model and inclusion mixture weights", {

  priors <- list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1)),
                 s = prior("cauchy", list(0, 1), list(0, 5)))
  attr(priors$b, "multiply_by") <- "s"
  alternative_height <- 5 / (sqrt(26) * sqrt(2 * pi) * atan(5))
  second <- priors
  second$b <- prior("point", list(0))
  attr(second$b, "multiply_by") <- "s"
  models <- lapply(names(priors), function(parameter){
    list(.set_prior_model_weight(priors[[parameter]], 3),
         .set_prior_model_weight(second[[parameter]], 1))
  })
  names(models) <- names(priors)
  context <- .prior_density_model_mixture_context(models, names(priors), n_grid = 512)
  density <- .prior_density_from_context(context, c(a = 1, b = -1))
  expected <- .75 * alternative_height + .25 * stats::dnorm(0)
  expect_equal(as.numeric(.prior_linear_density_height(density, 0)), expected, tolerance = 1e-9)
  expect_true(prior_density_ordinate(density, 0)$exact)

  priors$b <- prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.3)))
  attr(priors$b, "multiply_by") <- "s"
  context <- .prior_density_build_context(priors, names(priors), conditional = "b", n_grid = 512)
  density <- .prior_density_from_context(context, c(a = 1, b = -1))
  expect_equal(as.numeric(.prior_linear_density_height(density, 0)), alternative_height, tolerance = 1e-9)
  context <- .prior_density_build_context(priors, names(priors), n_grid = 512)
  density <- .prior_density_from_context(context, c(a = 1, b = -1))
  expect_equal(as.numeric(.prior_linear_density_height(density, 0)),
               .3 * alternative_height + .7 * stats::dnorm(0), tolerance = 1e-9)

  rows <- rbind(c(a = 1, b = -1), c(a = 1, b = 0))
  row_density <- .prior_density_from_context_rows(context, rows)
  expected <- .5 * (.3 * alternative_height + .7 * stats::dnorm(0)) + .5 * stats::dnorm(0)
  expect_equal(as.numeric(.prior_linear_density_height(row_density, 0)), expected, tolerance = 1e-9)
  expect_false(prior_density_ordinate(row_density, 0)$provenance$integration$exact)
})

test_that("conditional-normal classification survives exhausted numerical budgets", {

  priors <- list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1)),
                 s = prior("beta", list(.1, .1)))
  attr(priors$b, "multiply_by") <- "s"
  spec <- .prior_conditional_normal_spec(
    priors, .prior_linear_split_multiply_groups(priors, c(a = 1, b = 1)),
    c(a = NA_character_, b = NA_character_)
  )
  result <- .prior_conditional_normal_ordinate(spec, 0, n_grid = 21)
  expect_identical(result$behavior, "regular")
  expect_true(result$exact)
  expect_true(is.na(result$log_density))
  expect_false(result$provenance$integration$converged)
})

test_that("failed conditional-normal quadrature cannot fall back to grid heights", {

  priors <- list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1)),
                 s = prior("cauchy", list(0, 1), list(0, 5)))
  attr(priors$b, "multiply_by") <- "s"
  density <- .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 128)
  integration_reply <- list(value = .2, abs.error = .1, subdivisions = 1L,
                            message = "maximum number of subdivisions reached")
  testthat::local_mocked_bindings(
    integrate = function(...) integration_reply,
    .package = "stats"
  )
  testthat::local_mocked_bindings(
    .prior_linear_density_grid_height = function(...) stop("Grid fallback is forbidden."),
    .package = "BayesTools"
  )
  result <- prior_density_ordinate(density, 0)
  expect_identical(result$behavior, "regular")
  expect_true(result$exact)
  expect_true(is.na(result$log_density))
  expect_false(result$provenance$integration$converged)
  expect_error(.prior_linear_density_height(density, 0),
               "Conditional-normal prior density was rejected by diagnostics: integration reported",
               fixed = TRUE)
  integration_reply <- list(value = 0, abs.error = 0, subdivisions = 1L, message = "OK")
  zero <- prior_density_ordinate(density, 0)
  expect_identical(zero$behavior, "regular")
  expect_true(zero$exact)
  expect_true(is.na(zero$log_density))
  expect_error(.prior_linear_density_height(density, 0),
               "zero ordinate for a structurally positive density", fixed = TRUE)
})

test_that("conditional-normal mixtures preflight expansion before numerical evaluation", {

  slab <- prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.5)))
  attr(slab, "multiply_by") <- "s"
  names_b <- paste0("b", seq_len(25))
  priors <- c(list(a = prior("normal", list(0, 1)),
                   s = prior("cauchy", list(0, 1), list(0, 5))),
               stats::setNames(rep(list(slab), 25), names_b))
  weights <- c(a = 1, stats::setNames(rep(1, 25), names_b))
  calls <- 0L
  testthat::local_mocked_bindings(
    .prior_conditional_normal_ordinate = function(...) {
      calls <<- calls + 1L
      stop("Numerical expansion must not start.")
    }
  )
  result <- .prior_density_ordinate_linear_base(priors, weights, NULL, 0, n_grid = 512)
  expect_identical(result$behavior, "unknown")
  expect_identical(calls, 0L)
  priors$s <- prior("point", list(1))
  expect_identical(.prior_density_ordinate_linear_base(priors, weights, NULL, 0, n_grid = 512)$behavior, "unknown")
  expect_identical(calls, 0L)

  # Nested stored metadata are counted recursively even though constructors
  # currently reject nesting. Three positive leaves cannot fit 42 evaluations.
  nested <- prior_mixture(list(prior("point", list(0)), prior("normal", list(0, 1))))
  nested[[2L]] <- prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(0, 2))))
  attr(nested, "multiply_by") <- "s"
  priors <- list(a = prior("normal", list(0, 1)), b = nested,
                 s = prior("cauchy", list(0, 1), list(0, 5)))
  expect_identical(.prior_density_ordinate_linear_base(priors, c(a = 1, b = 1), NULL, 0, n_grid = 42)$behavior, "unknown")
  expect_identical(calls, 0L)
})

test_that("feasible conditional-normal leaves share one evaluation budget", {

  budget <- .prior_linear_density_default_grid()
  slab <- prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.5)))
  attr(slab, "multiply_by") <- "s"
  priors <- list(a = prior("normal", list(0, 1)),
                 s = prior("cauchy", list(0, 1), list(0, 5)),
                 b1 = slab, b2 = slab, b3 = slab)
  original <- .prior_conditional_normal_ordinate
  budgets <- evaluations <- numeric()
  testthat::local_mocked_bindings(
    .prior_conditional_normal_ordinate = function(spec, value, n_grid){
      result <- original(spec, value, n_grid)
      budgets <<- c(budgets, n_grid)
      evaluations <<- c(evaluations, result$provenance$integration$evaluations)
      result
    }
  )
  result <- .prior_density_ordinate_linear_base(priors, c(a = 1, b1 = 1, b2 = 1, b3 = 1), NULL, 0, n_grid = budget)
  expected <- sum(vapply(0:3, function(k){
    stats::dbinom(k, 3, .5) * stats::integrate(function(angle){
      stats::dnorm(0, sd = sqrt(1 + k * tan(angle)^2)) / atan(5)
    }, 0, atan(5), rel.tol = 1e-12)$value
  }, numeric(1)))
  expect_equal(exp(result$log_density), expected, tolerance = 1e-7)
  expect_length(budgets, 7L)
  expect_lte(sum(budgets), budget)
  expect_lte(sum(evaluations), budget)
  expect_lte(.prior_density_ordinate_integration(result$provenance)$budget, budget)
})

test_that("FFT removed-mass diagnostics have probability units", {

  set.seed(135)
  a <- stats::runif(64)
  b <- stats::runif(64)
  a[1:40] <- 0
  b[25:64] <- 0
  a <- a / sum(a)
  b <- b / sum(b)
  distribution <- function(y, mass, scale){
    list(
      density = list(x = (seq_along(y) - 1) * scale,
                     y = y / scale, mass = mass),
      points = data.frame(x = 0, p = 1 - mass), n_grid = length(y)
    )
  }
  clipping <- function(scale, mass_a = 1, mass_b = 1){
    result <- .prior_linear_density_convolve(
      distribution(a, mass_a, scale), distribution(b, mass_b, scale),
      dx = scale
    )
    attr(result, "fft_clipping")[[1L]]
  }
  original <- clipping(1)
  rescaled <- clipping(16)
  mixed <- clipping(16, .5, .25)
  expect_gt(original$clipped_value_count, 0L)
  expect_identical(rescaled$clipped_value_count, original$clipped_value_count)
  expect_equal(rescaled$clipped_negative_mass / original$clipped_negative_mass, 1)
  expect_equal(mixed$clipped_negative_mass / rescaled$clipped_negative_mass, .125)
  expect_equal(mixed$continuous_component_mass, .125)
})

test_that("underflow cannot turn a continuous prior into zero density", {

  priors <- list(a = prior("normal", list(mean = 0, sd = 1e200)),
                 b = prior("normal", list(mean = 0, sd = 1e200)))
  # The convolution is a proper normal with a representable positive height;
  # its FFT product underflows before multiplication by the grid spacing.
  expect_gt(stats::dnorm(0, sd = sqrt(2) * 1e200), 0)
  expect_error(
    .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 256),
    paste0("Continuous prior density is unavailable because its numerical grid ",
           "has no finite positive mass. Rescale the modeled quantity and its ",
           "prior parameters before evaluating this density."),
    fixed = TRUE
  )
})

test_that("unrepresentable continuous grids cannot become atoms or non-finite densities", {

  normal <- prior("normal", list(mean = 1e20, sd = 1))
  expect_identical(prior_density_ordinate(normal, 1e20)$behavior, "regular")
  expect_error(
    .prior_linear_combination_density(list(a = normal), c(a = 1), n_grid = 256),
    paste0("Continuous prior density is unavailable because its numerical range ",
           "collapses to one representable value. Center or rescale the modeled ",
           "quantity and its prior parameters before evaluating this density."),
    fixed = TRUE
  )
  point <- .prior_linear_combination_density(
    list(a = prior("spike", list(location = 1e20))), c(a = 1), n_grid = 256
  )
  expect_equal(point$points, data.frame(x = 1e20, p = 1))

  narrow <- list(a = prior("normal", list(mean = 0, sd = 1e-200)),
                 b = prior("normal", list(mean = 0, sd = 1e-200)))
  expect_true(is.finite(stats::dnorm(0, sd = sqrt(2) * 1e-200)))
  expect_error(
    .prior_linear_combination_density(narrow, c(a = 1, b = 1), n_grid = 256),
    paste0("FFT prior convolution is unavailable because its numerical values ",
           "are non-finite. Rescale the modeled quantity and its prior parameters ",
           "before evaluating this density."),
    fixed = TRUE
  )
})

test_that("boundary-singular prior densities keep exact edge-cell masses", {

  # Shapes below one give an integrable infinite density at a support bound.
  # The bound knot carries the exact CDF mass of its grid cell. The relative
  # tolerance covers the grid renormalisation, which absorbs the O(sqrt(dx))
  # midpoint error of the neighbouring singular cells (about 1.2e-3 for
  # gamma(0.5, 1) on the default grid).
  singular <- list(
    list(prior = prior("beta", list(0.5, 0.5)), lower = TRUE, upper = TRUE),
    list(prior = prior("gamma", list(0.5, 1)), lower = TRUE, upper = FALSE),
    list(prior = prior("beta", list(2, 0.7)), lower = FALSE, upper = TRUE),
    list(prior = prior("beta", list(0.5, 1.5)), lower = TRUE, upper = FALSE)
  )
  for(case in singular){
    density <- .prior_linear_combination_density(list(p = case$prior), c(p = 1))
    x  <- density$density$x
    y  <- density$density$y
    dx <- x[2] - x[1]
    expect_true(all(is.finite(y)))
    if(case$lower){
      expect_equal(x[1], 0)
      expect_equal(y[1] * dx, cdf(case$prior, dx / 2), tolerance = 5e-3)
    }
    if(case$upper){
      expect_equal(x[length(x)], 1)
      expect_equal(y[length(y)] * dx, ccdf(case$prior, 1 - dx / 2), tolerance = 5e-3)
    }
    expect_equal(
      .prior_linear_density_height(density, .3),
      pdf(case$prior, .3)
    )

    sum_density <- .prior_linear_combination_density(
      list(p = case$prior, q = prior("normal", list(0, 1))),
      c(p = 1, q = 1)
    )
    expect_true(all(is.finite(sum_density$density$y)))
    context <- .prior_density_build_context(list(p = case$prior), "p")
    expect_s3_class(.prior_density_from_context(context, c(p = 1)), "prior_linear_density")
  }

  fit <- coda::mcmc(cbind(p = stats::qgamma(stats::ppoints(64), .5, 1)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(p = prior("gamma", list(.5, 1)))
  fit <- .bt_attach_parameter_map(fit, monitor_names = "p")
  marginal <- marginal_posterior(
    as_mixed_posteriors(fit, "p"), "p", use_formula = FALSE,
    prior_samples = TRUE, n_samples = 128
  )
  expect_equal(
    .prior_linear_density_height(attr(marginal, "prior_density"), .5),
    stats::dgamma(.5, .5, 1)
  )
})

test_that("boundary-singular densities keep the bound knot when the grid ends one rounding error short", {

  # seq(0, 1, by = 1 / (n - 1)) stops at 1 - 2^-53 for n = 500 and 3000. The
  # knot next to the singular bound then kept a finite but huge ordinate
  # (about 3e7 for beta(.5, .5)) that took almost all of the grid mass after
  # normalisation (grid heights off by -96% to -100%); it now sits on the bound
  # and carries its cell mass. Reference: the prior density; 3e-3 covers the
  # grid renormalisation of the O(sqrt(dx)) neighbouring-cell error (observed
  # at most 1.4e-3 at n = 500).
  cases <- list(
    list(prior = prior("beta", list(.5, .5)), n_grid = 500L),
    list(prior = prior("beta", list(2, .7)),  n_grid = 3000L)
  )
  for(case in cases){
    expect_lt(max(seq(0, 1, by = 1 / (case$n_grid - 1L))), 1)
    density <- .prior_linear_combination_density(
      list(p = case$prior), c(p = 1), n_grid = case$n_grid
    )
    expect_true(all(is.finite(density$density$y)))
    for(value in c(.1, .5, .9)){
      expect_equal(.prior_linear_density_grid_height(density, value),
                   pdf(case$prior, value), tolerance = 3e-3)
    }
  }
})

test_that("exp_lin output transformations keep analytic limits at a zero source knot", {

  # y = 2 x^b: at the source knot x = 0 the density is f(0) / 2 for b = 1, zero
  # for 0 < b < 1, and the knot is dropped for b > 1 (singular) and b < 0
  # (image at infinity). Grid heights are compared with the analytic density
  # of the transformed half-normal / lognormal; 1e-3 covers the source-grid
  # normalisation (captured mass and Riemann versus trapezoid rule).
  half_normal <- prior("normal", list(0, 1), truncation = list(0, Inf))
  lognormal   <- prior("lognormal", list(0, 1))
  transformed <- function(source, b){
    .prior_linear_combination_density(
      list(mu = source), c(mu = 1), output_transformation = "exp_lin",
      output_transformation_arguments = list(a = log(2), b = b)
    )
  }
  height <- function(density, value){
    .prior_linear_density_grid_height(density, value)
  }
  exact <- function(density_x, value, b){
    source <- (value / 2)^(1 / b)
    density_x(source) / abs(2 * b * source^(b - 1))
  }
  half_normal_pdf <- function(x) 2 * stats::dnorm(x)

  identity_map <- transformed(half_normal, 1)
  expect_equal(identity_map$density$x[1], 0)
  expect_equal(height(identity_map, 0), stats::dnorm(0), tolerance = 1e-3)
  expect_equal(height(identity_map, 1), exact(half_normal_pdf, 1, 1), tolerance = 1e-3)

  root_map <- transformed(half_normal, .5)
  expect_equal(root_map$density$x[1], 0)
  expect_identical(root_map$density$y[1], 0)
  expect_equal(height(root_map, 1), exact(half_normal_pdf, 1, .5), tolerance = 1e-3)

  square_map <- transformed(half_normal, 2)
  expect_gt(square_map$density$x[1], 0)
  expect_equal(height(square_map, 1), exact(half_normal_pdf, 1, 2), tolerance = 1e-3)

  inverse_map <- transformed(half_normal, -1)
  expect_true(all(is.finite(inverse_map$density$x)))
  expect_equal(height(inverse_map, 1), exact(half_normal_pdf, 1, -1), tolerance = 1e-3)

  lognormal_map <- transformed(lognormal, 1)
  expect_identical(lognormal_map$density$y[1], 0)
  expect_equal(height(lognormal_map, 1), stats::dlnorm(.5) / 2, tolerance = 1e-3)

  # Saturating transformations of heavy-tailed priors omit the knots whose
  # transformed ordinate is not representable instead of failing.
  cauchy <- list(mu = prior("cauchy", list(0, .707)))
  tanh_map <- .prior_linear_combination_density(cauchy, c(mu = 1),
                                               output_transformation = "tanh")
  expect_true(all(abs(tanh_map$density$x) <= 1))
  expect_true(all(is.finite(tanh_map$density$y)))
  exp_map <- .prior_linear_combination_density(cauchy, c(mu = 1),
                                              output_transformation = "exp")
  expect_true(all(is.finite(exp_map$density$x) & exp_map$density$x >= 0))
  expect_true(all(is.finite(exp_map$density$y)))
})

test_that("nonlinear output transformations keep dense knots of wide priors", {

  # Heights against the analytic lognormal / tanh-normal densities. The
  # source grid omits 1e-4 per tail and is renormalised, a relative bias of
  # about 2e-4; 1e-3 bounds it with margin.
  exp_map <- .prior_linear_combination_density(
    list(mu = prior("normal", list(0, 8))), c(mu = 1), output_transformation = "exp"
  )
  expect_equal(.prior_linear_density_grid_height(exp_map, 1),
               stats::dlnorm(1, 0, 8), tolerance = 1e-3)
  expect_equal(.prior_linear_density_grid_height(exp_map, 3),
               stats::dlnorm(3, 0, 8), tolerance = 1e-3)

  tanh_map <- .prior_linear_combination_density(
    list(mu = prior("normal", list(0, 4))), c(mu = 1), output_transformation = "tanh"
  )
  expect_equal(.prior_linear_density_grid_height(tanh_map, 0),
               stats::dnorm(0, 0, 4), tolerance = 1e-3)
  expect_equal(.prior_linear_density_grid_height(tanh_map, .5),
               stats::dnorm(atanh(.5), 0, 4) / (1 - .5^2), tolerance = 1e-3)
  expect_true(all(diff(tanh_map$density$x) > 0))
})

test_that("heavy-tailed combinations resolve their narrowest source and mixtures stay bounded", {

  # Normal-type combinations keep the default grid (the narrowest source
  # spans more than 64 of its knots per robust SD).
  normal_sum <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1))),
    c(a = 1, b = 1)
  )
  expect_equal(diff(normal_sum$density$x[1:2]),
               diff(range(normal_sum$density$x)) / (length(normal_sum$density$x) - 1))
  expect_equal(diff(normal_sum$density$x[1:2]),
               2 * 2 * stats::qnorm(1 - 1e-4) / (4096 - 1), tolerance = 1e-10)

  # A Cauchy slab sets a range of about +/-2250; the N(0, .2) source is
  # resolved with at least 64 knots per SD (previously dx was about 1.1 and
  # the height 0.315). Reference: 0.5 phi(0; .2) + 0.5 (N(0, .2) * Cauchy)(0)
  # by quadrature; 1e-3 covers the omitted 1e-4 tails.
  slab <- prior_spike_and_slab(prior("cauchy", list(0, .707)),
                               prior_inclusion = prior("spike", list(.5)))
  reference <- .5 * stats::dnorm(0, 0, .2) + .5 * stats::integrate(
    function(t) stats::dnorm(t, 0, .2) * stats::dcauchy(t, 0, .707),
    -Inf, Inf, rel.tol = 1e-12
  )$value
  heavy <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, .2)), b = slab), c(a = 1, b = 1)
  )
  expect_lte(diff(heavy$density$x[1:2]), .2 / 64)
  expect_equal(.prior_linear_density_grid_height(heavy, 0), reference, tolerance = 1e-3)

  # Refinement halves the source spacing and omits less tail probability, by
  # 10 when polynomial tails would more than double the range. t3 + N(0, .1)
  # converges to its quadrature reference; Cauchy combinations exceed the grid
  # limit before converging and stop loudly instead of reporting a height
  # biased by their omitted tail mass.
  student <- .prior_linear_combination_density(
    list(a = prior("t", list(0, 1, 3)), b = prior("normal", list(0, .1))), c(a = 1, b = 1)
  )
  student_height <- .prior_linear_density_height(student, 0)
  expect_true(isTRUE(attr(student_height, "adaptive_evaluation")$converged))
  expect_equal(
    as.numeric(student_height),
    stats::integrate(function(t) stats::dt(t, 3) * stats::dnorm(-t, 0, .1),
                     -Inf, Inf, rel.tol = 1e-12)$value,
    tolerance = 1e-4
  )
  expect_error(
    .prior_linear_density_height(heavy, 0),
    "Adaptive prior-density evaluation did not converge within the documented grid-refinement error criterion.",
    fixed = TRUE
  )

  # Density jumps (half-normal, truncated normal and uniform components)
  # converge once the spacing strictly halves; references by quadrature.
  half_normal <- prior("normal", list(0, 1), truncation = list(0, Inf))
  jumps <- list(
    list(priors = list(a = half_normal, b = half_normal), value = .5,
         reference = stats::integrate(function(t) 2 * stats::dnorm(t) * 2 * stats::dnorm(.5 - t),
                                      0, .5, rel.tol = 1e-12)$value),
    list(priors = list(a = half_normal, b = prior("normal", list(0, 1))), value = -1,
         reference = stats::integrate(function(t) 2 * stats::dnorm(t) * stats::dnorm(-1 - t),
                                      0, Inf, rel.tol = 1e-12)$value),
    list(priors = list(a = prior("uniform", list(0, 1)), b = prior("normal", list(0, .3))), value = -.2,
         reference = stats::integrate(function(t) stats::dnorm(-.2 - t, 0, .3),
                                      0, 1, rel.tol = 1e-12)$value)
  )
  for(jump in jumps){
    density <- .prior_linear_combination_density(jump$priors, c(a = 1, b = 1))
    height  <- .prior_linear_density_height(density, jump$value)
    expect_true(isTRUE(attr(height, "adaptive_evaluation")$converged))
    expect_equal(as.numeric(height), jump$reference, tolerance = 1e-4)
  }

  # Mixing a narrow and a Cauchy model needs about 5e7 knots at the finest
  # spacing (previously 2.4 GB); it now stops before allocating.
  context <- .prior_density_build_context(
    list(mu = list(prior("normal", list(0, .05), prior_weights = 1),
                   prior("cauchy", list(0, .707), prior_weights = 1))),
    "mu"
  )
  expect_error(
    .prior_density_from_context(context, c(mu = 1)),
    "Mixed prior density is unavailable because its components have incompatible scales",
    fixed = TRUE
  )
})

test_that("mixture grids beyond the limit end adaptive refinement as non-convergence", {

  # The mixing guard is a grid-limit condition.
  empty <- data.frame(x = numeric(), p = numeric())
  narrow <- list(density = list(x = seq(-1e-3, 1e-3, length.out = 101), y = rep(500, 101), mass = 1),
                 points = empty, n_grid = 101L)
  wide <- list(density = list(x = seq(-50, 50, length.out = 101), y = rep(.01, 101), mass = 1),
               points = empty, n_grid = 101L)
  expect_error(
    .prior_linear_density_mix(list(narrow, wide), c(.5, .5), dx = 2e-5),
    class = "BayesTools_prior_grid_limit"
  )

  # The initial model mixture fits the grid (about 4.8e5 knots). After halving
  # the spacing each model stays below the limit, but their union at the finest
  # model spacing needs about 3.3e6 knots: refinement ends and the height is
  # reported as not converged, instead of failing with the mixing error meant
  # for incompatible scales in the requested density itself.
  context <- .prior_density_build_context(
    list(a = list(prior("normal", list(0, .01), prior_weights = 1),
                  prior("t", list(0, 1, 3), prior_weights = 1)),
         b = list(prior("gamma", list(3, 2), prior_weights = 1),
                  prior("normal", list(0, 1), prior_weights = 1))),
    c("a", "b")
  )
  density <- .prior_density_from_context(context, c(a = 1, b = 1))
  expect_lt(length(density$density$x), .prior_linear_density_max_grid())
  expect_error(
    .prior_linear_density_height(density, .3),
    "Adaptive prior-density evaluation did not converge within the documented grid-refinement error criterion.",
    fixed = TRUE
  )
})

test_that("linear group ranges accept omitted source transformations", {

  group <- list(prior = prior("normal", list(0, 1)), weights = c(mu = 1), indices = 1L)
  expect_equal(
    BayesTools:::.prior_linear_group_range(group, tail_prob = .001),
    stats::qnorm(c(.001, .999))
  )
})

test_that("adaptive ordinates stop when a halved grid spacing exceeds the grid limit", {

  density <- structure(
    list(density = list(x = c(-1, 1), y = c(1, 1), mass = 1), points = NULL),
    class = "prior_linear_density"
  )
  attr(density, "adaptive_evaluation") <- list(
    kind = "linear_combination",
    arguments = list(prior_list = list(a = prior("normal", list(0, 1)),
                                       b = prior("uniform", list(0, 1))),
                     weights = c(a = 1, b = 1), n_grid = 4096, tail_prob = 1e-4)
  )
  attr(density, "grid_resolution") <- c(spacing = 1e-6, n_grid = 2097152)
  requested <- NULL
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(prior_list, weights, n_grid, tail_prob,
                                                 grid_spacing, .record_evaluation){
      requested <<- c(tail_prob = tail_prob, grid_spacing = grid_spacing)
      stop(BayesTools:::.prior_linear_density_grid_limit_error(2, grid_spacing))
    },
    .package = "BayesTools"
  )
  expect_null(BayesTools:::.prior_linear_density_refinement(density))
  expect_equal(requested, c(tail_prob = 1e-7, grid_spacing = 5e-7))
  expect_error(
    BayesTools:::.prior_linear_density_height(density, 0),
    "Adaptive prior-density evaluation did not converge within the documented grid-refinement error criterion.",
    fixed = TRUE
  )
  side <- hypothesis_parse("theta < 0")$statements[[1L]]$left
  expect_error(
    BayesTools:::.hypothesis_prior_density_prob(density, side, "theta"),
    "Adaptive prior-probability evaluation did not converge within the documented grid-refinement error criterion.",
    fixed = TRUE
  )
})

test_that("linear density normalize warns only for large mass deviation", {

  quiet <- BayesTools:::.prior_linear_density_normalize(list(
    density = list(x = 0:1, y = c(1, 1), mass = 1 + 1e-8),
    points = NULL
  ), warn = TRUE)
  expect_equal(quiet$density$mass, 1, tolerance = 1e-12)

  # Intermediate partial masses stay quiet by default.
  expect_silent(
    BayesTools:::.prior_linear_density_normalize(list(
      density = list(x = 0:1, y = c(1, 1), mass = 0.5),
      points = NULL
    ))
  )

  expect_warning(
    noisy <- BayesTools:::.prior_linear_density_normalize(list(
      density = list(x = 0:1, y = c(1, 1), mass = 2.5),
      points = data.frame(x = 0.5, p = 0.5)
    ), warn = TRUE),
    "renormalizing to 1",
    fixed = TRUE
  )
  expect_equal(noisy$density$mass + sum(noisy$points$p), 1, tolerance = 1e-12)

  expect_error(
    BayesTools:::.prior_linear_density_normalize(list(
      density = list(x = 0:1, y = c(0, 0), mass = 0),
      points = NULL
    )),
    "zero total mass",
    fixed = TRUE
  )
})

test_that("linear prior density matches analytic normal sums", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    weights = c(x = 2, y = -0.5),
    n_grid  = 1024
  )

  x <- density$density$x
  y <- density$density$y
  expected <- stats::dnorm(x, mean = -0.5, sd = sqrt(2^2 + 1^2))
  expected <- expected / (sum(expected) * (x[2] - x[1]))

  expect_equal(
    sum(abs(y - expected)) * (x[2] - x[1]),
    0,
    tolerance = 0.03
  )
  expect_equal(density$density$mass, 1, tolerance = 1e-8)
})

test_that("linear prior ordinates adapt across center and omitted tails", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(0, 1))
    ),
    weights   = c(x = 1, y = 1),
    n_grid    = 512,
    tail_prob = 1e-3
  )

  # A normal sum has an exact structural ordinate, which is returned directly
  # instead of refining the grid.
  skewed <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("gamma", list(3, 2))
    ),
    weights   = c(x = 1, y = 1),
    n_grid    = 512,
    tail_prob = 1e-3
  )
  skewed_reference <- function(value){
    stats::integrate(function(t) stats::dnorm(value - t) * stats::dgamma(t, 3, 2),
                     0, Inf, rel.tol = 1e-12)$value
  }

  evaluate_density <- BayesTools:::.prior_linear_combination_density
  refinement_calls <- 0L
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(...) {

      refinement_calls <<- refinement_calls + 1L
      evaluate_density(...)
    },
    .package = "BayesTools"
  )

  exact_center <- BayesTools:::.prior_linear_density_height(density, 0)
  exact_tail   <- BayesTools:::.prior_linear_density_height(density, 8)
  expect_identical(refinement_calls, 0L)
  expect_equal(exact_center, stats::dnorm(0, sd = sqrt(2)), tolerance = 1e-12)
  expect_equal(exact_tail, stats::dnorm(8, sd = sqrt(2)), tolerance = 1e-12)

  # A normal plus gamma sum has no structural ordinate and is refined; the
  # reference is the convolution integral.
  center <- BayesTools:::.prior_linear_density_height(skewed, 0)
  expect_identical(refinement_calls, 2L)
  # Each refinement halves the source spacing (1024 knots over [-3.09, 8.69]
  # initially) and omits 1000 times less tail probability.
  expect_equal(
    attr(center, "adaptive_evaluation")[c("n_grid", "tail_prob", "refinements")],
    list(n_grid = 16384L, tail_prob = 1e-9, refinements = 2L)
  )
  refinement_calls <- 0L
  tail <- BayesTools:::.prior_linear_density_height(skewed, 9)
  expect_gt(9, max(skewed$density$x))
  expect_identical(refinement_calls, 3L)
  expect_lt(
    abs(as.numeric(center) / skewed_reference(0) - 1),
    1e-4
  )
  expect_lt(
    abs(as.numeric(tail) / skewed_reference(9) - 1),
    1e-4
  )
  expect_true(isTRUE(attr(center, "adaptive_evaluation")$converged))
  expect_true(isTRUE(attr(tail, "adaptive_evaluation")$converged))
  expect_gt(attr(tail, "adaptive_evaluation")$refinements, 0)

  refinement_calls <- 0L
  side <- hypothesis_parse("theta < 0")$statements[[1L]]$left
  probability <- BayesTools:::.hypothesis_prior_density_prob(
    density, side, "theta"
  )
  expect_lt(abs(probability / .5 - 1), 1e-4)
  expect_identical(refinement_calls, 1L)

  diagnostics <- attr(density, "numerical_diagnostics", exact = TRUE)
  expect_equal(diagnostics$tail_probability_per_source, 1e-3)
  expect_equal(
    diagnostics$intended_captured_probability_per_continuous_source,
    .998
  )
  expect_true(is.list(diagnostics$grid_normalization))
  expect_true(is.list(diagnostics$fft_clipping))
})

test_that("prior heights use exact ordinates at density jumps and zero outside support", {

  # Model-averaged prior with a component truncated at the null: the grid
  # cannot converge across the jump, the exact mixture ordinate is the
  # right-limit 0.5 phi(0) + 0.5 phi(0; .5, 1) / (1 - Phi(0; .5, 1)).
  prior_list <- list(mu = list(
    prior("normal", list(0, 1), prior_weights = 1),
    prior("normal", list(.5, 1), truncation = list(0, Inf), prior_weights = 1)
  ))
  context <- .prior_density_build_context(prior_list, "mu")
  mixture <- .prior_density_from_context(context, c(mu = 1))
  jump <- .5 * stats::dnorm(0) + .5 * stats::dnorm(0, .5, 1) / stats::pnorm(0, .5, 1, lower.tail = FALSE)
  expect_equal(.prior_linear_density_height(mixture, 0), jump, tolerance = 1e-12)
  expect_equal(exp(prior_density_ordinate(mixture, 0)$log_density), jump, tolerance = 1e-12)
  expect_equal(.prior_linear_density_height(mixture, -.3), .5 * stats::dnorm(-.3),
               tolerance = 1e-12)

  fit <- coda::mcmc(cbind(mu = stats::qnorm(stats::ppoints(64), .2, .5)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)))
  fit <- .bt_attach_parameter_map(fit, monitor_names = "mu")
  posterior <- marginal_posterior(as_mixed_posteriors(fit, "mu"), "mu",
                                  use_formula = FALSE, prior_samples = TRUE,
                                  n_samples = 64)
  attr(posterior, "prior_density") <- mixture
  bf <- Savage_Dickey_BF(posterior, null_hypothesis = 0, silent = TRUE)
  expect_true(is.finite(bf) && bf > 0)

  # Outside the exactly known support the continuous height is zero, with or
  # without an exact structural ordinate.
  half_normal <- prior("normal", list(0, 1), truncation = list(0, Inf))
  shifted <- .prior_linear_combination_density(
    list(p = prior("point", list(.5)), a = half_normal), c(p = 1, a = 1)
  )
  expect_identical(.prior_linear_density_height(shifted, .4), 0)
  convolved <- .prior_linear_combination_density(
    list(a = half_normal, b = half_normal), c(a = 1, b = 1)
  )
  expect_identical(prior_density_ordinate(convolved, -1)$behavior, "unknown")
  expect_identical(.prior_linear_density_height(convolved, -1), 0)
  transformed <- .prior_linear_combination_density(
    list(a = half_normal, b = half_normal), c(a = 1, b = 1),
    output_transformation = "exp"
  )
  expect_identical(.prior_linear_density_height(transformed, .5), 0)
  expect_identical(
    .prior_linear_density_support_hull(attr(transformed, "adaptive_evaluation")),
    c(1, Inf)
  )

  # An ordered term has no exactly known support here, so no zero is implied.
  ordered <- prior_ordered(prior("normal", list(0, 1)))
  attr(ordered, "levels") <- 3
  ordered <- .prior_ordered_default_bound(ordered, "mu_f")
  expect_null(.prior_linear_combination_support_hull(
    list(mu_intercept = half_normal, mu_f = ordered),
    c(mu_intercept = 1, "mu_f[1]" = 1)
  ))
})

test_that("linear prior density handles multiply_by products and point mass", {

  priors <- list(
    beta = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE)
    ),
    sigma = prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(priors$beta, "multiply_by") <- "sigma"

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = c(beta = 1),
    n_grid     = 512
  )

  expect_equal(BayesTools:::.prior_linear_density_point_mass(density, 0), 0.5, tolerance = 1e-8)
  expect_gt(BayesTools:::.prior_linear_density_height(density, 0), 0)
})

test_that("a multiply_by scale with its own weight is rejected as dependent", {

  # x * sigma + sigma shares sigma between both terms (density at 1: 0.280 by
  # Monte Carlo); convolving them as independent gave 0.383.
  x <- prior("normal", list(0, 1))
  attr(x, "multiply_by") <- "sigma"
  priors <- list(x = x, sigma = prior("normal", list(0, 1), truncation = list(0, Inf)))
  expect_error(
    .prior_linear_combination_density(priors, c(x = 1, sigma = 1)),
    paste0("The prior density of this linear combination is unavailable because 'sigma' ",
           "enters it both as the 'multiply_by' scale of other coefficients and with its ",
           "own weight, which makes the terms dependent. Evaluate the terms separately."),
    fixed = TRUE
  )
  context <- .prior_density_build_context(priors, c("x", "sigma"))
  expect_error(.prior_density_from_context(context, c(x = 1, sigma = 1)),
               "enters it both as the 'multiply_by' scale", fixed = TRUE)

  # Each term alone keeps its density.
  expect_s3_class(.prior_linear_combination_density(priors, c(x = 1), n_grid = 256),
                  "prior_linear_density")
  expect_s3_class(.prior_linear_combination_density(priors, c(sigma = 1), n_grid = 256),
                  "prior_linear_density")
})

test_that("product densities count each mixed-measure component once", {

  make_dist <- function(inclusion){
    p <- prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior("spike", list(inclusion))
    )
    BayesTools:::.prior_linear_group_distribution(
      group = list(
        prior = p,
        weights = c(x = 1),
        indices = 1L
      ),
      dx = .02,
      tail_prob = 1e-4,
      source_transforms = c(x = NA_character_),
      n_grid = 256
    )
  }

  lhs <- make_dist(.5)
  rhs <- make_dist(.75)
  product <- BayesTools:::.prior_linear_density_product(
    lhs,
    rhs,
    n_grid = 256
  )

  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(product, 0),
    1 - .5 * .75,
    tolerance = 1e-12
  )
  expect_equal(product$density$mass, .5 * .75, tolerance = 1e-12)

  scaled <- BayesTools:::.prior_linear_density_product(
    lhs,
    BayesTools:::.prior_linear_density_point(.25),
    n_grid = 256
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(scaled, 0),
    .5,
    tolerance = 1e-12
  )
  expect_equal(scaled$density$mass, .5, tolerance = 1e-12)
})

test_that("linear prior density treats zero-weight combinations as point priors", {

  priors <- list(beta = prior("normal", list(0, 1)))

  empty_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = numeric(),
    n_grid     = 128
  )
  expect_equal(BayesTools:::.prior_linear_density_point_mass(empty_density, 0), 1)
  expect_null(empty_density$density)

  zero_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = c(beta = 0),
    n_grid     = 128
  )
  expect_equal(BayesTools:::.prior_linear_density_point_mass(zero_density, 0), 1)
  expect_null(zero_density$density)

  context <- BayesTools:::.prior_density_context(
    prior_list   = priors,
    column_names = "beta",
    n_grid       = 128
  )
  context_density <- BayesTools:::.prior_density_from_context(
    context,
    weights = c(beta = 0)
  )
  expect_equal(BayesTools:::.prior_linear_density_point_mass(context_density, 0), 1)
  expect_null(context_density$density)
})

test_that("linear prior density preserves every finite nonzero coefficient", {

  priors <- list(beta = prior("normal", list(0, 1)))
  tiny_weight <- .Machine$double.eps

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = c(beta = tiny_weight),
    n_grid     = 256
  )

  expect_false(is.null(density$density))
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(density, 0),
    0
  )
  expect_equal(
    range(density$density$x) / tiny_weight,
    stats::qnorm(c(1e-4, 1 - 1e-4)),
    tolerance = 1e-8
  )

  support <- BayesTools:::.posterior_support_from_prior_list_weights(
    priors,
    c(beta = tiny_weight)
  )
  expect_equal(support$bounds, c(-Inf, Inf))
})


test_that("linear prior density rejects non-finite coefficient weights", {

  priors <- list(
    x = prior("normal", list(0, 1)),
    y = prior("normal", list(0, 1))
  )

  for(bad_weight in c(NA_real_, NaN)){
    expect_error(
      BayesTools:::.prior_linear_combination_density(
        prior_list = priors,
        weights    = c(x = 1, y = bad_weight),
        n_grid     = 128
      ),
      "The 'weights' argument cannot contain NA/NaN values.",
      fixed = TRUE
    )
  }

  for(bad_weight in c(Inf, -Inf)){
    expect_error(
      BayesTools:::.prior_linear_combination_density(
        prior_list = priors,
        weights    = c(x = 1, y = bad_weight),
        n_grid     = 128
      ),
      "The 'weights' argument must contain only finite values.",
      fixed = TRUE
    )
  }
})


test_that("linear prior density honors named scalar source transforms", {
  p <- prior("lognormal", list(0, 1))

  expect_equal(
    BayesTools:::.prior_linear_scalar_range(
      p,
      weight = 1,
      tail_prob = .025,
      source_transform = c(beta = "log")
    ),
    stats::qnorm(c(.025, .975)),
    tolerance = 1e-8
  )
})


test_that("row-wise prior densities mix row predictions, not averaged weights", {
  context <- BayesTools:::.prior_density_context(
    prior_list   = list(beta = prior("normal", list(0, 1))),
    column_names = "beta",
    n_grid       = 1024
  )

  row_density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(c(1, 2), ncol = 1, dimnames = list(NULL, "beta"))
  )
  averaged_density <- BayesTools:::.prior_density_from_context(
    context,
    weights = c(beta = 1.5)
  )

  density_second_moment <- function(d){
    x <- d$density$x
    y <- d$density$y
    dx <- x[2] - x[1]
    sum(x^2 * y) * dx + sum(d$points$x^2 * d$points$p)
  }

  expect_equal(density_second_moment(row_density), mean(c(1^2, 2^2)), tolerance = .08)
  expect_equal(density_second_moment(averaged_density), 1.5^2, tolerance = .08)
})

test_that("row-wise prior densities transform the row mixture once", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(mu_intercept = prior("normal", list(0, 1)),
                        mu_x         = prior("normal", list(0, 1))),
    column_names = c("mu_intercept", "mu_x"),
    n_grid       = 4096
  )
  weights <- rbind(c(mu_intercept = 1, mu_x = .5), c(mu_intercept = 1, mu_x = 2))
  sds <- sqrt(1 + weights[, "mu_x"]^2)

  # Rows are N(0, 1 + w^2) on the linear-predictor scale, so the transformed
  # mixture is an equal mixture of lognormal (exp) or tanh-normal densities.
  # The grid keeps the size of an untransformed row mixture (previously 25.7M
  # knots for two rows); 1e-3 covers the omitted 1e-4 source tails.
  exp_density <- BayesTools:::.prior_density_from_context_rows(
    context, weights, output_transformation = "exp"
  )
  expect_lte(length(exp_density$density$x), 2 * 4096)
  for(value in c(.5, 1, 2)){
    expect_equal(BayesTools:::.prior_linear_density_grid_height(exp_density, value),
                 mean(stats::dlnorm(value, 0, sds)), tolerance = 1e-3)
  }
  tanh_density <- BayesTools:::.prior_density_from_context_rows(
    context, weights, output_transformation = "tanh"
  )
  for(value in c(-.5, 0, .5)){
    expect_equal(BayesTools:::.prior_linear_density_grid_height(tanh_density, value),
                 mean(stats::dnorm(atanh(value), 0, sds)) / (1 - value^2),
                 tolerance = 1e-3)
  }

  symmetric <- rbind(c(mu_intercept = 1, mu_x = .5), c(mu_intercept = 1, mu_x = -.5))
  symmetric_density <- BayesTools:::.prior_density_from_context_rows(
    context, symmetric, output_transformation = "exp"
  )
  expect_lte(length(symmetric_density$density$x), 2 * 4096)
})

test_that("row-wise prior densities preserve bitwise-distinct design rows", {

  context <- BayesTools:::.prior_density_context(
    prior_list = list(beta = prior("point", list(location = 1))),
    column_names = "beta",
    n_grid = 64
  )
  nearby <- 1 + .Machine$double.eps
  density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(
      c(1, nearby),
      ncol = 1,
      dimnames = list(NULL, "beta")
    )
  )

  expect_identical(density$points$x, c(1, nearby))
  expect_equal(density$points$p, c(.5, .5))
})

test_that("conditional log-intercept prior densities mix the conditioned models", {

  formula_scale <- list(mu = list(mu_x = list(mean = 5, sd = 2)))
  attr(formula_scale$mu, "log_intercept") <- TRUE
  slab <- prior("normal", list(0, 1))
  prior_list <- list(
    mu_intercept = prior("lognormal", list(0, .5)),
    mu_x = prior_spike_and_slab(slab, prior_inclusion = prior("point", list(.5)))
  )
  columns <- c("mu_intercept", "mu_x")

  conditional <- .generate_transformed_prior_densities(
    prior_list, columns, formula_scale, conditional = "mu_x"
  )
  # Reference: the unconditional density of the prior list filtered to the
  # alternative (slab) component, computed through the unconditional bypass.
  filtered <- .generate_transformed_prior_densities(
    list(mu_intercept = prior_list$mu_intercept, mu_x = slab), columns, formula_scale
  )
  unconditional <- .generate_transformed_prior_densities(prior_list, columns, formula_scale)
  for(value in c(.5, 1, 1.5, 3)){
    expect_equal(
      .prior_linear_density_grid_height(conditional$mu_intercept, value),
      .prior_linear_density_grid_height(filtered$mu_intercept, value),
      tolerance = 1e-8
    )
    expect_equal(
      as.numeric(.prior_linear_density_height(conditional$mu_intercept, value)),
      as.numeric(.prior_linear_density_height(filtered$mu_intercept, value)),
      tolerance = 1e-8
    )
  }
  expect_false(isTRUE(all.equal(
    .prior_linear_density_grid_height(conditional$mu_intercept, 1),
    .prior_linear_density_grid_height(unconditional$mu_intercept, 1)
  )))

  # Public route: conditional unscaling of a log(intercept) fit.
  posterior <- cbind(
    mu_intercept   = stats::qlnorm(stats::ppoints(20), 0, .5),
    mu_x           = c(rep(0, 8), seq(.1, 1, length.out = 12)),
    mu_x_indicator = c(rep(0L, 8), rep(1L, 12))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- formula_scale
  fit <- .bt_attach_parameter_map(fit, monitor_names = colnames(posterior))
  mixed <- as_mixed_posteriors(fit, columns, conditional = "mu_x", transform_scaled = TRUE)
  expect_equal(
    as.numeric(.prior_linear_density_height(attr(mixed, "prior_densities")$mu_intercept, 1)),
    as.numeric(.prior_linear_density_height(filtered$mu_intercept, 1)),
    tolerance = 1e-8
  )
})

test_that("plot_transformed_prior exposes transformed prior plotting as a public wrapper", {

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(1, .25), prior_weights = 3)
      ),
      is_null = c(TRUE, FALSE)
    )
  )
  formula_scale <- list(mu = list(mu_x = list(mean = 5, sd = 2)))
  attr(formula_scale$mu, "log_intercept") <- FALSE

  plot <- plot_transformed_prior(
    prior_list    = prior_list,
    column_names  = c("mu_intercept", "mu_x"),
    formula_scale = formula_scale,
    parameter     = "mu_x",
    n_points      = 128,
    plot_type     = "ggplot",
    par_name      = "x"
  )

  expect_s3_class(plot, "ggplot")
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1))))
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomSegment"), logical(1))))

  device_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(device_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot_transformed_prior(
    prior_list    = prior_list,
    column_names  = c("mu_intercept", "mu_x"),
    formula_scale = formula_scale,
    parameter     = "mu_x",
    n_points      = 64,
    plot_type     = "base"
  ))
})

test_that("plot_transformed_prior returns NULL only for untransformed identity parameters", {

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x         = prior("normal", list(0, 1))
  )

  expect_null(plot_transformed_prior(
    prior_list   = prior_list,
    column_names = c("mu_intercept", "mu_x"),
    parameter    = "mu_x",
    plot_type    = "ggplot"
  ))

  transformed_plot <- plot_transformed_prior(
    prior_list     = prior_list,
    column_names   = c("mu_intercept", "mu_x"),
    parameter      = "mu_x",
    transformation = "exp",
    x_range        = c(-2, 2),
    n_points       = 64,
    plot_type      = "ggplot"
  )

  expect_s3_class(transformed_plot, "ggplot")
  expect_error(
    plot_transformed_prior(
      prior_list   = prior_list,
      column_names = c("mu_intercept", "mu_x"),
      parameter    = "mu_z"
    ),
    "Parameter 'mu_z' was not found in 'column_names'.",
    fixed = TRUE
  )
})

test_that("plot_transformed_prior draws raw coefficients without multiply_by", {

  # the monitored slope is its own N(0, 1) prior; multiply_by = "sigma" scales
  # only the linear predictor
  x_prior <- prior("normal", list(0, 1))
  attr(x_prior, "multiply_by") <- "sigma"
  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x         = x_prior,
    sigma        = prior("lognormal", list(0, .5))
  )
  formula_scale <- list(mu = list(mu_x = list(mean = 3, sd = 2)))
  attr(formula_scale$mu, "log_intercept") <- FALSE
  columns <- c("mu_intercept", "mu_x", "sigma")

  # original-scale raw slope b / s ~ N(0, 1 / s); raw intercept
  # b0 - (m / s) b ~ N(0, sqrt(1 + (m / s)^2))
  expected <- list(
    mu_x         = function(x) stats::dnorm(x, 0, 1 / 2),
    mu_intercept = function(x) stats::dnorm(x, 0, sqrt(1 + (3 / 2)^2))
  )
  densities <- .generate_transformed_prior_densities(prior_list, columns, formula_scale)
  for(parameter in names(expected)){
    for(value in c(0, .3, 1)){
      expect_equal(
        as.numeric(.prior_linear_density_height(densities[[parameter]], value)),
        expected[[parameter]](value),
        tolerance = 1e-6
      )
    }
    # plotted grid values; 1e-3 covers the display interpolation (~1e-4)
    plot <- plot_transformed_prior(
      prior_list, columns, formula_scale, parameter,
      n_points = 101, plot_type = "ggplot"
    )
    line <- ggplot2::ggplot_build(plot)$data[[1]]
    expect_lt(max(abs(line$y - expected[[parameter]](line$x))), 1e-3)
  }
})
