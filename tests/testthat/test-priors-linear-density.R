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

test_that("linear group ranges accept omitted source transformations", {

  group <- list(prior = prior("normal", list(0, 1)), weights = c(mu = 1), indices = 1L)
  expect_equal(
    BayesTools:::.prior_linear_group_range(group, tail_prob = .001),
    stats::qnorm(c(.001, .999))
  )
})

test_that("adaptive ordinates cannot converge by repeating the capped grid", {

  density <- structure(
    list(density = list(x = c(-1, 1), y = c(1, 1), mass = 1), points = NULL),
    class = "prior_linear_density"
  )
  attr(density, "adaptive_evaluation") <- list(
    kind = "linear_combination", arguments = list(n_grid = 32768, tail_prob = 1e-12)
  )
  expect_null(BayesTools:::.prior_linear_density_refinement(density))
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

  attr(density, "adaptive_evaluation")$arguments$n_grid <- 16384L
  # The refined grid must describe a different distribution: region
  # probabilities are ratios of grid integrals, so a pure rescaling of the
  # ordinates would legitimately converge.
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(n_grid, tail_prob, .record_evaluation){
      density$density$y <- c(1, 3)
      density
    },
    .package = "BayesTools"
  )
  expect_error(
    BayesTools:::.prior_linear_density_height(density, 0),
    "did not converge", fixed = TRUE
  )
  expect_error(
    BayesTools:::.hypothesis_prior_density_prob(density, side, "theta"),
    "did not converge", fixed = TRUE
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

  evaluate_density <- BayesTools:::.prior_linear_combination_density
  refinement_calls <- 0L
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(...) {

      refinement_calls <<- refinement_calls + 1L
      evaluate_density(...)
    },
    .package = "BayesTools"
  )

  center <- BayesTools:::.prior_linear_density_height(density, 0)
  expect_identical(refinement_calls, 2L)
  expect_equal(
    attr(center, "adaptive_evaluation")[c("n_grid", "tail_prob", "refinements")],
    list(n_grid = 8192L, tail_prob = 1e-9, refinements = 2L)
  )
  refinement_calls <- 0L
  tail <- BayesTools:::.prior_linear_density_height(density, 8)
  expect_identical(refinement_calls, 4L)
  expect_lt(
    abs(as.numeric(center) / stats::dnorm(0, sd = sqrt(2)) - 1),
    1e-4
  )
  expect_lt(
    abs(as.numeric(tail) / stats::dnorm(8, sd = sqrt(2)) - 1),
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
