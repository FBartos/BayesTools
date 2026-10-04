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

  # without the additive normal term the product is a pure scale mixture; the
  # half-Cauchy multiplier's positive density at zero makes it infinite there
  pure_product <- .prior_linear_combination_density(priors, c(b = 1), n_grid = 128)
  expect_identical(prior_density_ordinate(pure_product, 0)$behavior, "infinite")
  expect_identical(prior_density_ordinate(pure_product, 0)$method, "conditional_normal_mixture")
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
    # each piece of the split integral has the full budget
    integration <- ordinate$provenance$integration
    expect_identical(integration$evaluations, counted)
    expect_identical(sum(integration$piece_evaluations), counted)
    expect_true(all(integration$piece_evaluations <= 512))
    expect_equal(as.numeric(.prior_linear_density_height(density, 0)), expected[[i]], tolerance = 1e-7)
  }
})

test_that("legacy product flags cannot establish an infinite density", {

  # A non-normal additive term with a product has no structural route; its
  # capped product grid is not used for heights either. Products are not
  # flagged from their grids, and stale flags are ignored: b * s with b, s ~
  # N(0, 1) has the density K0(|x|) / pi, infinite only at 0.
  priors <- list(a = prior("uniform", list(-1, 1)), b = prior("normal", list(0, 1)),
                 s = prior("normal", list(0, 1)))
  attr(priors$b, "multiply_by") <- "s"
  unsupported <- .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 128)
  expect_null(attr(unsupported, "singular_density_points"))
  expect_identical(prior_density_ordinate(unsupported, 0)$behavior, "unknown")
  expect_error(.prior_linear_density_height(unsupported, 0),
               paste0("The prior density of this linear combination is unavailable: its ",
                      "'multiply_by' product has no structural density route"),
               fixed = TRUE)
  pure <- .prior_linear_combination_density(priors, c(b = 1), n_grid = 128)
  expect_identical(prior_density_ordinate(pure, 0)$behavior, "infinite")
  expect_identical(.prior_linear_density_height(pure, 0), Inf)
  attr(pure, "singular_density_points") <- c(0, 1)
  expect_identical(.prior_linear_density_height(pure, 0), Inf)
  expect_equal(as.numeric(.prior_linear_density_height(pure, 1)), besselK(1, 0) / pi,
               tolerance = 1e-8)
})

test_that("exp_lin images of a zero source have their structural boundary limit", {

  # Y = exp(a) X^b with X's density ~ C x^(b - 1) at 0 has the density
  # C / (b exp(a)) at 0: X ~ Beta(1/2, 1) (C = 1 / B(1/2, 1) = 1/2) with
  # sqrt(4 X), and X ~ gamma(1/2, 2) (C = sqrt(2) / Gamma(1/2)) with
  # 3 sqrt(X); references in closed form, also for the plotted curve.
  cases <- list(
    list(prior = prior("beta", list(.5, 1)), a = log(4) / 2, limit = .5 / (.5 * 2),
         density = function(y) stats::dbeta(y^2 / 4, .5, 1) * y / 2),
    list(prior = prior("gamma", list(.5, 2)), a = log(3),
         limit = sqrt(2) / gamma(.5) / (.5 * 3),
         density = function(y) stats::dgamma((y / 3)^2, .5, 2) * 2 * y / 9)
  )
  for(case in cases){
    density <- .prior_linear_combination_density(
      list(x = case$prior), c(x = 1),
      output_transformation = "exp_lin",
      output_transformation_arguments = list(a = case$a, b = .5)
    )
    ordinate <- prior_density_ordinate(density, 0)
    expect_identical(ordinate$behavior, "regular")
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), case$limit, tolerance = 1e-14)
    expect_equal(as.numeric(.prior_linear_density_height(density, 0)), case$limit,
                 tolerance = 1e-14)
    curve <- .prior_linear_density_to_plot_data(density, n_points = 11, x_range = c(0, 1))$density
    # the curve leaves the support bound 0 by the edge (0, 0) -> (0, limit)
    expect_identical(curve$x[1:2], c(0, 0))
    expect_equal(curve$y, c(0, case$limit, case$density(curve$x[-(1:2)])), tolerance = 1e-12)
  }
})

test_that("products without an additive normal term are pure scale mixtures", {

  # b * s with b ~ N(0, 1): f(x) = int phi(x / s) / s f_s(s) ds, and at the
  # offset 0 f(0) = phi(0) E[1 / s], finite exactly when f_s vanishes at zero
  # (closed-form inverse moments: gamma rate / (shape - 1), lognormal
  # exp(-mu + sigma^2 / 2), inverse gamma shape / scale). References:
  # integrate() at rel.tol 1e-12 with stats:: densities, split around the
  # scale peak |x|. The capped product grid gave -20.3% at 0 (gamma), +8.2% at
  # .1 (half-normal) and stopped for the lognormal and inverse-gamma scales.
  split_reference <- function(f, points){
    sum(vapply(seq_len(length(points) - 1L), function(i){
      stats::integrate(f, points[i], points[i + 1L], rel.tol = 1e-12,
                       subdivisions = 5000L)$value
    }, numeric(1)))
  }
  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "s"
  cases <- list(
    list(s = prior("gamma", list(2, 2)), f = function(s) stats::dgamma(s, 2, 2),
         upper = Inf, at_zero = stats::dnorm(0) * 2),
    list(s = prior("lognormal", list(0, 1)), f = function(s) stats::dlnorm(s),
         upper = Inf, at_zero = stats::dnorm(0) * exp(.5)),
    list(s = prior("invgamma", list(1, .15)), f = function(s) stats::dgamma(1 / s, 1, .15) / s^2,
         upper = Inf, at_zero = stats::dnorm(0) / .15),
    list(s = prior("normal", list(0, 1), list(0, Inf)), f = function(s) 2 * stats::dnorm(s),
         upper = Inf, at_zero = Inf),
    list(s = prior("uniform", list(0, 1)), f = function(s) stats::dunif(s),
         upper = 1, at_zero = Inf)
  )
  for(case in cases){
    density <- .prior_linear_combination_density(list(b = slope, s = case$s), c(b = 1))
    zero <- prior_density_ordinate(density, 0)
    expect_true(zero$exact)
    expect_identical(zero$method, "conditional_normal_mixture")
    expect_identical(zero$behavior, if(is.finite(case$at_zero)) "regular" else "infinite")
    expect_equal(as.numeric(.prior_linear_density_height(density, 0)), case$at_zero, tolerance = 1e-12)
    for(value in c(-.1, .05, 1, 3)){
      points <- sort(unique(c(0, abs(value) * c(.1, 1, 10), 20)))
      points <- c(points[points < case$upper], case$upper)
      reference <- split_reference(function(s) stats::dnorm(value / s) / s * case$f(s), points)
      ordinate <- prior_density_ordinate(density, value)
      expect_identical(ordinate$behavior, "regular")
      expect_true(ordinate$provenance$integration$converged)
      expect_equal(as.numeric(.prior_linear_density_height(density, value)), reference,
                   tolerance = 1e-8)
    }
  }

  # regions: P(b s > c) = int P(b > c / s) f_s(s) ds; the grid gave +1.3%
  # and -15.5% for the central interval
  density <- .prior_linear_combination_density(list(b = slope, s = prior("gamma", list(2, 2))), c(b = 1))
  upper_tail <- function(c0){
    split_reference(function(s) stats::pnorm(c0 / s, lower.tail = FALSE) * stats::dgamma(s, 2, 2),
                    c(0, .1, 1, 3, Inf))
  }
  probability <- function(hypothesis){
    .hypothesis_prior_density_prob(density, hypothesis_parse(hypothesis)$statements[[1L]]$left, "theta")
  }
  for(c0 in c(.05, .5, 3)){
    expect_equal(probability(paste("theta >", c0)), upper_tail(c0), tolerance = 1e-8)
  }
  expect_equal(probability("theta > -0.05 & theta < 0.05"), 1 - 2 * upper_tail(.05), tolerance = 1e-8)

  # a deterministic additive part and a nonzero multiplied mean: at the offset
  # .2, f = phi(b_m / b_s) / b_s E[1 / s]
  shifted_slope <- prior("normal", list(.4, .3))
  attr(shifted_slope, "multiply_by") <- "s"
  shifted <- .prior_linear_combination_density(
    list(p = prior("point", list(.2)), b = shifted_slope, s = prior("gamma", list(3, 2))),
    c(p = 1, b = 1)
  )
  expect_equal(as.numeric(.prior_linear_density_height(shifted, .2)),
               stats::dnorm(.4 / .3) / .3 * 2 / (3 - 1), tolerance = 1e-12)
  for(value in c(.5, 1)){
    expect_equal(
      as.numeric(.prior_linear_density_height(shifted, value)),
      split_reference(function(s) stats::dnorm(value, .2 + .4 * s, .3 * s) * stats::dgamma(s, 3, 2),
                      c(0, .1, .5, 1, 3, Inf)),
      tolerance = 1e-8
    )
  }

  # a point multiplier or a deterministic multiplied part is an affine term
  fixed_scale <- .prior_linear_combination_density(list(b = slope, s = prior("point", list(2))), c(b = 1))
  expect_identical(prior_density_ordinate(fixed_scale, .5)$method, "scalar_affine")
  expect_equal(as.numeric(.prior_linear_density_height(fixed_scale, .5)), stats::dnorm(.5, 0, 2),
               tolerance = 1e-14)
  fixed_slope <- prior("point", list(.5))
  attr(fixed_slope, "multiply_by") <- "s"
  scaled_multiplier <- .prior_linear_combination_density(
    list(b = fixed_slope, s = prior("gamma", list(2, 2))), c(b = 1)
  )
  expect_identical(prior_density_ordinate(scaled_multiplier, .3)$method, "scalar_affine")
  expect_equal(as.numeric(.prior_linear_density_height(scaled_multiplier, .3)),
               2 * stats::dgamma(.6, 2, 2), tolerance = 1e-14)

  # a non-normal additive term with a product has no structural route and no
  # grid height or grid probability
  general <- .prior_linear_combination_density(
    list(a = prior("t", list(0, 1, 5)), b = slope, s = prior("lognormal", list(0, 1))),
    c(a = 1, b = 1)
  )
  expect_identical(prior_density_ordinate(general, 1)$behavior, "unknown")
  message <- paste0("The prior density of this linear combination is unavailable: its ",
                    "'multiply_by' product has no structural density route")
  expect_error(.prior_linear_density_height(general, 1), message, fixed = TRUE)
  expect_error(
    .hypothesis_prior_density_prob(general, hypothesis_parse("theta > 0")$statements[[1L]]$left, "theta"),
    message, fixed = TRUE
  )
})

test_that("products of a non-normal term and a scale are scale-mixture ordinates", {

  # b * s with b ~ t(0, 1, 5) and s ~ gamma(3, 2): f(x) = int f_t(x / s) /
  # s f_s(s) ds and f(0) = f_t(0) E[1 / s]. References: integrate() at
  # rel.tol 1e-12 with stats:: densities.
  slope <- prior("t", list(0, 1, 5))
  attr(slope, "multiply_by") <- "s"
  density <- .prior_linear_combination_density(list(b = slope, s = prior("gamma", list(3, 2))), c(b = 1))
  zero <- prior_density_ordinate(density, 0)
  expect_identical(zero$method, "scale_mixture")
  expect_equal(exp(zero$log_density), stats::dt(0, 5) * 2 / (3 - 1), tolerance = 1e-12)
  for(value in c(-.7, .2, 4)){
    reference <- sum(vapply(list(c(0, .1), c(.1, 1), c(1, 5), c(5, Inf)), function(piece){
      stats::integrate(function(s) stats::dt(value / s, 5) / s * stats::dgamma(s, 3, 2),
                       piece[1L], piece[2L], rel.tol = 1e-12)$value
    }, numeric(1)))
    expect_equal(as.numeric(.prior_linear_density_height(density, value)), reference, tolerance = 1e-8)
  }
  side <- hypothesis_parse("theta > 0.5")$statements[[1L]]$left
  expect_equal(
    .hypothesis_prior_density_prob(density, side, "theta"),
    stats::integrate(function(s) stats::pt(.5 / s, 5, lower.tail = FALSE) * stats::dgamma(s, 3, 2),
                     0, Inf, rel.tol = 1e-12)$value,
    tolerance = 1e-8
  )

  # a positive factor and a positive scale: zero density below the offset, and
  # at the offset f_s(0) E[1 / b] when only the scale is positive at zero
  positive <- prior("lognormal", list(0, .5))
  attr(positive, "multiply_by") <- "s"
  density <- .prior_linear_combination_density(list(b = positive, s = prior("exp", list(1))), c(b = 1))
  expect_identical(prior_density_ordinate(density, -.1)$behavior, "zero")
  expect_equal(as.numeric(.prior_linear_density_height(density, 0)), exp(.5^2 / 2), tolerance = 1e-12)
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
  # a rejected quadrature is not exact (exact = TRUE never comes with a
  # missing regular log density)
  expect_false(result$exact)
  expect_true(is.na(result$log_density))
  expect_false(result$provenance$integration$converged)
  expect_match(result$reason, "The prior-density quadrature was rejected by its diagnostics: integration reported",
               fixed = TRUE)
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
  expect_false(result$exact)
  expect_true(is.na(result$log_density))
  expect_false(result$provenance$integration$converged)
  expect_match(result$reason, "maximum number of subdivisions reached", fixed = TRUE)
  expect_error(.prior_linear_density_height(density, 0),
               "Conditional-normal prior density was rejected by diagnostics: integration reported",
               fixed = TRUE)
  integration_reply <- list(value = 0, abs.error = 0, subdivisions = 1L, message = "OK")
  zero <- prior_density_ordinate(density, 0)
  expect_identical(zero$behavior, "regular")
  expect_false(zero$exact)
  expect_true(is.na(zero$log_density))
  expect_match(zero$reason, "zero ordinate for a structurally positive density", fixed = TRUE)
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
  # currently reject nesting. Three positive leaves exceed the leaf cap of
  # 42 / 15 = 2 initial quadrature rules.
  nested <- prior_mixture(list(prior("point", list(0)), prior("normal", list(0, 1))))
  nested[[2L]] <- prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(0, 2))))
  attr(nested, "multiply_by") <- "s"
  priors <- list(a = prior("normal", list(0, 1)), b = nested,
                 s = prior("cauchy", list(0, 1), list(0, 5)))
  expect_identical(.prior_density_ordinate_linear_base(priors, c(a = 1, b = 1), NULL, 0, n_grid = 42)$behavior, "unknown")
  expect_identical(calls, 0L)
})

test_that("feasible conditional-normal leaves each receive the full evaluation budget", {

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
  expect_true(all(budgets == budget))
  expect_true(all(evaluations <= budget))
  integration <- .prior_density_ordinate_integration(result$provenance)
  expect_identical(integration$budget, 7 * budget)
  expect_identical(integration$evaluations, sum(evaluations))
  expect_true(integration$converged)
})

test_that("conditional-normal quadrature accepts a result converged at the budget", {

  # One QUADPACK piece of the split integral: b * s with s ~ half-normal over
  # the whole support (0, Inf). At 0 the integral is exp(z) K0(z) / (2 pi)
  # with z = 1 / 4 (s = sinh(t) gives the K0 integral,
  # https://dlmf.nist.gov/10.32.E9). The one-sided 15-point rule converges
  # after three intervals, i.e. 15 * (2 * 3 - 1) = 75 evaluations.
  multiplier <- prior("normal", list(0, 1), list(0, Inf))
  integrand <- function(s){
    exp(stats::dnorm(0, 0, sqrt(1 + s^2), log = TRUE) + lpdf(multiplier, s))
  }
  expected <- besselK(.25, 0, expon.scaled = TRUE) / (2 * pi)
  tolerance <- .prior_linear_density_refinement_tolerance()

  exact_budget <- .prior_conditional_normal_piece(
    integrand, 0, Inf, n_grid = 75,
    relative = tolerance$relative, absolute = tolerance$quadrature_floor
  )
  expect_identical(exact_budget$message, "OK")
  expect_identical(exact_budget$evaluations, 75L)
  expect_equal(exact_budget$value, expected, tolerance = 1e-8)

  # one evaluation short of the third interval: the budget, not QUADPACK's
  # interval limit, stops the quadrature
  short_budget <- .prior_conditional_normal_piece(
    integrand, 0, Inf, n_grid = 74,
    relative = tolerance$relative, absolute = tolerance$quadrature_floor
  )
  expect_identical(short_budget$message, "the integration evaluation budget was exhausted")
  expect_lte(short_budget$evaluations, 74L)

  # the ordinate splits the same integral into pieces, each with the budget
  priors <- list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1)),
                 s = multiplier)
  attr(priors$b, "multiply_by") <- "s"
  spec <- .prior_conditional_normal_spec(
    priors, .prior_linear_split_multiply_groups(priors, c(a = 1, b = 1)),
    c(a = NA_character_, b = NA_character_)
  )
  ordinate <- .prior_conditional_normal_ordinate(spec, 0, n_grid = 75)
  integration <- ordinate$provenance$integration
  expect_true(integration$converged)
  expect_true(all(integration$piece_evaluations <= 75L))
  expect_identical(length(integration$piece_evaluations), length(integration$breakpoints) - 1L)
  expect_equal(exp(ordinate$log_density), expected, tolerance = 1e-8)
})

test_that("conditional-normal quadratures split scale-disparate integrals at breakpoints", {

  # Split-integral references: integrate() at rel.tol 1e-12 over hand-chosen
  # pieces around the Gaussian peak and the other term's mass (beta(.5, .5):
  # u = sin(t)^2 removes the bound singularities; normal x half-normal
  # multiplier at 0: exp(z) K0(z) / (2 pi sigma) with z = 1 / (4 sigma^2),
  # https://dlmf.nist.gov/10.32.E9). One QUADPACK call over the whole support
  # missed the mass of all but the beta cases.
  split_reference <- function(f, points){
    sum(vapply(seq_len(length(points) - 1L), function(i){
      stats::integrate(f, points[i], points[i + 1L], rel.tol = 1e-12,
                       subdivisions = 5000L)$value
    }, numeric(1)))
  }
  density_of <- function(priors){
    .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 4096)
  }
  height <- function(density, value){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$method, "conditional_normal_mixture")
    expect_true(ordinate$provenance$integration$converged)
    as.numeric(.prior_linear_density_height(density, value))
  }
  cases <- list(
    list(priors = list(a = prior("normal", list(0, .001)), b = prior("gamma", list(3, 2))), value = 1,
         reference = split_reference(function(u) stats::dnorm(1 - u, 0, .001) * stats::dgamma(u, 3, 2),
                                     c(0, .98, 1, 1.02, Inf))),
    list(priors = list(a = prior("normal", list(0, .1)), b = prior("cauchy", list(0, 1))), value = 5,
         reference = split_reference(function(u) stats::dnorm(5 - u, 0, .1) * stats::dcauchy(u),
                                     c(-Inf, 0, 3, 5, 7, Inf))),
    list(priors = list(a = prior("normal", list(0, 1)), b = prior("lognormal", list(3, .01))), value = 20,
         reference = split_reference(function(u) stats::dnorm(20 - u) * stats::dlnorm(u, 3, .01),
                                     c(0, exp(3) - 1, exp(3), exp(3) + 1, Inf))),
    list(priors = list(a = prior("normal", list(0, .3)), b = prior("beta", list(.5, .5))), value = 0,
         reference = 2 / pi * split_reference(function(t) stats::dnorm(0 - sin(t)^2, 0, .3),
                                              c(0, pi / 4, pi / 2))),
    list(priors = list(a = prior("normal", list(0, .01)), b = prior("beta", list(.5, .5))), value = .999,
         reference = 2 / pi * split_reference(function(t) stats::dnorm(.999 - sin(t)^2, 0, .01),
                                              c(0, asin(sqrt(.9)), asin(sqrt(.999)), pi / 2))),
    list(priors = list(a = prior("normal", list(0, 10)), b = prior("gamma", list(1e4, 1e4))), value = 1,
         reference = split_reference(function(u) stats::dnorm(1 - u, 0, 10) * stats::dgamma(u, 1e4, 1e4),
                                     c(0, .9, 1, 1.1, Inf)))
  )
  for(case in cases){
    expect_equal(height(density_of(case$priors), case$value), case$reference, tolerance = 1e-8)
  }
  # the references agree with the review's values (5 significant digits)
  expect_equal(cases[[1L]]$reference, .54134, tolerance = 1e-4)
  expect_equal(cases[[2L]]$reference, .012256, tolerance = 1e-4)
  expect_equal(cases[[3L]]$reference, .38973, tolerance = 1e-4)

  # normal x half-normal multiplier with a very narrow and a very wide scale
  for(sigma in c(1e-3, 1e3)){
    slope <- prior("normal", list(0, 1))
    attr(slope, "multiply_by") <- "s"
    priors <- list(a = prior("normal", list(0, 1)), b = slope,
                   s = prior("normal", list(0, sigma), list(0, Inf)))
    z <- 1 / (4 * sigma^2)
    density <- density_of(priors)
    expect_equal(height(density, 0), besselK(z, 0, expon.scaled = TRUE) / (2 * pi * sigma),
                 tolerance = 1e-8)
    expect_equal(
      height(density, 1.5),
      split_reference(function(m) stats::dnorm(1.5, 0, sqrt(1 + m^2)) * 2 * stats::dnorm(m, 0, sigma),
                      c(0, sigma * c(.1, 1, 3, 10), Inf)),
      tolerance = 1e-8
    )
  }
})

test_that("multiplier quadratures resolve the location peak of a narrow multiplied normal", {

  # a + b * s with b ~ N(b_m, b_s) and b_s <= |b_m| / 2: in s the integrand
  # peaks at s* = (v - a_m) / b_m with local SD sqrt(a_s^2 + b_s^2 s*^2) / |b_m|.
  # References: integrate() at rel.tol 1e-12 with stats:: densities over pieces
  # around that peak; they agree with the review's values to the digits given.
  split_reference <- function(f, points){
    sum(vapply(seq_len(length(points) - 1L), function(i){
      stats::integrate(f, points[i], points[i + 1L], rel.tol = 1e-12,
                       subdivisions = 5000L)$value
    }, numeric(1)))
  }
  density_of <- function(priors){
    .prior_linear_combination_density(priors, c(a = 1, b = 1), n_grid = 4096)
  }
  height <- function(density, value){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$method, "conditional_normal_mixture")
    expect_true(ordinate$provenance$integration$converged)
    as.numeric(.prior_linear_density_height(density, value))
  }
  scale_mixture <- function(a, b, s){
    attr(b, "multiply_by") <- "s"
    list(a = a, b = b, s = s)
  }

  # gamma(3, 2) multiplier: without a breakpoint at the peak every piece missed
  # it (4.9e-115 instead of .1226)
  priors <- scale_mixture(prior("normal", list(0, .001)), prior("normal", list(3, .001)),
                          prior("gamma", list(3, 2)))
  width <- sqrt(.001^2 + (.001 * .5)^2) / 3
  reference <- split_reference(
    function(s) stats::dnorm(1.5, 3 * s, sqrt(.001^2 + (.001 * s)^2)) * stats::dgamma(s, 3, 2),
    c(0, .5 + c(-30, -10, -3, -1, 0, 1, 3, 10, 30) * width, Inf)
  )
  expect_equal(reference, .1226265, tolerance = 1e-6)
  expect_equal(height(density_of(priors), 1.5), reference, tolerance = 1e-8)

  # inverse-gamma(3, 2) multiplier at the image of its median: the median
  # breakpoint alone found half of the peak (.14695 instead of .29388)
  priors <- scale_mixture(prior("normal", list(.2, 1e-4)), prior("normal", list(3, 1e-4)),
                          prior("invgamma", list(shape = 3, scale = 2)))
  median <- 2 / stats::qgamma(.5, 3)
  value <- .2 + 3 * median
  width <- sqrt(1e-8 + (1e-4 * median)^2) / 3
  reference <- split_reference(
    function(s) stats::dnorm(value, .2 + 3 * s, sqrt(1e-8 + (1e-4 * s)^2)) * stats::dgamma(1 / s, 3, 2) / s^2,
    c(0, median + c(-30, -10, -3, -1, 0, 1, 3, 10, 30) * width, Inf)
  )
  expect_equal(reference, .29388, tolerance = 1e-4)
  expect_equal(height(density_of(priors), value), reference, tolerance = 1e-8)

  # b_s much larger than |b_m| (33 times): the window around s* is not a peak,
  # and breakpoints there made the pieces miss mass elsewhere (-0.7% at -1,
  # -2.3e-4 at -.5)
  priors <- scale_mixture(prior("normal", list(.2, 1e-4)), prior("normal", list(3, 100)),
                          prior("normal", list(0, 1e-3), list(0, Inf)))
  density <- density_of(priors)
  for(value in c(-1, -.5)){
    reference <- split_reference(
      function(s) stats::dnorm(value, .2 + 3 * s, sqrt(1e-8 + (100 * s)^2)) * 2 * stats::dnorm(s, 0, 1e-3),
      c(0, 1e-3 * c(1e-3, .01, .1, .5, 1, 2, 3, 5, 10), Inf)
    )
    expect_equal(height(density, value), reference, tolerance = 1e-8)
  }

  # just beyond the earlier guard b_s <= |b_m| / 10 the integrand is still a
  # narrow peak at s* (the standardized distance grows to |b_m| / b_s ~ 8-10
  # away from it); without peak breakpoints every piece missed it (1e-21 to
  # 3e-18 instead of .036-.113, reported as converged). a ~ N(.2, 1e-4),
  # b ~ N(3, b_s), value .2 (s* = 0). References: 30-digit mpmath tanh-sinh
  # quadrature over dense breakpoints (densities written from their formulas),
  # agreeing with an R QUADPACK reference to 1e-8.
  cases <- list(
    list(s = prior("t", list(1, .5, 3)), b_s = 3 * (.1 * (1 + 1e-12)), reference = .045470734073157164),
    list(s = prior("normal", list(1, .5)), b_s = 3 * .11, reference = .036446362420311275),
    list(s = prior("uniform", list(-1, 2)), b_s = 3 * .12, reference = .11278578719472659)
  )
  for(case in cases){
    priors <- scale_mixture(prior("normal", list(.2, 1e-4)), prior("normal", list(3, case$b_s)), case$s)
    expect_equal(height(density_of(priors), .2), case$reference, tolerance = 1e-8)
  }
})

test_that("conditional-normal breakpoints keep their distance from bounds with infinite density", {

  # A breakpoint much closer to such a bound than the next breakpoint made
  # QUADPACK integrate its piece as if the singularity were at the breakpoint,
  # counting the mass next to the bound twice (2-8% too high); strongly
  # singular betas stopped with non-finite function values. References remove
  # the singularities: for gamma(shape, rate), u = x^(1 / shape) gives
  # g(u) du = rate^shape / Gamma(shape + 1) exp(-rate u) dx; for beta(a, b),
  # u = x^(1 / a) below 1 / 2 and 1 - u = y^(1 / b) above it; integrate() at
  # rel.tol 1e-12 with stats:: densities. They agree with the review's
  # parabolic-cylinder and mpmath values to the digits given.
  split_reference <- function(f, points){
    sum(vapply(seq_len(length(points) - 1L), function(i){
      stats::integrate(f, points[i], points[i + 1L], rel.tol = 1e-12,
                       subdivisions = 5000L)$value
    }, numeric(1)))
  }
  gamma_reference <- function(shape, rate, conditional, points){
    rate^shape / gamma(shape + 1) * split_reference(function(x){
      u <- x^(1 / shape)
      conditional(u) * exp(-rate * u)
    }, c(0, points^shape, Inf))
  }
  beta_reference <- function(a, b, conditional, points){
    edge <- c(0, 1e-6, 1e-3, .01, .1, .3, 1) * .5
    points <- points[points > 0 & points < 1]
    lower <- sort(unique(c(edge, points[points < .5])))
    upper <- sort(unique(c(edge, 1 - points[points > .5])))
    lower <- split_reference(function(x){
      u <- x^(1 / a)
      conditional(u) * (1 - u)^(b - 1)
    }, lower^a) / (a * beta(a, b))
    upper <- split_reference(function(y){
      distance <- y^(1 / b)
      conditional(1 - distance) * (1 - distance)^(a - 1)
    }, upper^b) / (b * beta(a, b))
    lower + upper
  }
  density <- function(priors, weights){
    .prior_linear_combination_density(priors, weights, n_grid = 4096)
  }
  height <- function(density, value){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$method, "conditional_normal_mixture")
    expect_true(ordinate$provenance$integration$converged)
    as.numeric(.prior_linear_density_height(density, value))
  }

  # N(.3, 1) + gamma(.15, 1) at .3
  reference <- gamma_reference(.15, 1, function(u) stats::dnorm(-u), c(1e-3, .1, 1, 3, 10))
  expect_equal(reference, .3817084456, tolerance = 1e-9)
  convolution <- density(list(a = prior("normal", list(.3, 1)), b = prior("gamma", list(.15, 1))),
                         c(a = 1, b = 1))
  expect_equal(height(convolution, .3), reference, tolerance = 1e-8)

  # N(.3, .001) + gamma(.15, .1) at .3 - .001
  reference <- gamma_reference(.15, .1, function(u) stats::dnorm(-.001 - u, 0, .001),
                               c(1e-4, 1e-3, 3e-3, .01, .03, .1))
  expect_equal(reference, 58.15411, tolerance = 1e-6)
  convolution <- density(list(a = prior("normal", list(.3, .001)), b = prior("gamma", list(.15, .1))),
                         c(a = 1, b = 1))
  expect_equal(height(convolution, .299), reference, tolerance = 1e-8)

  # N(.3, .01) + .001 gamma(.1, .1) at .31: the Gaussian-peak point
  # (v - m) / w - s / w is a rounding residue of about 1e-14 above the bound
  reference <- gamma_reference(.1, .1, function(u) stats::dnorm(.01 - .001 * u, 0, .01),
                               c(1e-3, 1, 5, 10, 20, 40, 110, 300))
  expect_equal(reference, 25.6337125, tolerance = 1e-8)
  convolution <- density(list(a = prior("normal", list(.3, .01)), b = prior("gamma", list(.1, .1))),
                         c(a = 1, b = .001))
  expect_equal(height(convolution, .31), reference, tolerance = 1e-8)

  # scale mixture N(0, 1) + N(0, 1) * s with s ~ gamma(.15, 1) at 0
  reference <- gamma_reference(.15, 1, function(s) stats::dnorm(0, 0, sqrt(1 + s^2)),
                               c(1e-3, .1, 1, 3, 10, 30))
  expect_equal(reference, .3859880, tolerance = 1e-6)
  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "s"
  mixture <- density(list(a = prior("normal", list(0, 1)), b = slope, s = prior("gamma", list(.15, 1))),
                     c(a = 1, b = 1))
  expect_equal(height(mixture, 0), reference, tolerance = 1e-8)

  # beta(.1, .1) and beta(.3, .3) plus N(.3, s), which stopped on non-finite
  # function values; digits are lost in 1 - u next to the upper bound (errors
  # up to about 1e-7), so the values are checked at the documented criterion
  tolerance <- .prior_linear_density_refinement_tolerance()
  for(shape in c(.1, .3)){
    for(sd in c(.001, .1, 1)){
      convolution <- density(list(a = prior("normal", list(.3, sd)), b = prior("beta", list(shape, shape))),
                             c(a = 1, b = 1))
      for(offset in c(0, .5, 1)){
        reference <- beta_reference(shape, shape, function(u) stats::dnorm(offset - u, 0, sd),
                                    offset + c(-10, -3, -1, 0, 1, 3, 10) * sd)
        expect_lt(abs(height(convolution, .3 + offset) - reference),
                  tolerance$relative * reference)
      }
    }
  }

  # A narrow Gaussian peak next to such a bound: its breakpoints (at least one
  # local SD from the bound) are kept by the spacing rule. Dropping them left
  # the peak inside the piece from the bound, which missed it: 6e-66 to 7e-45
  # instead of 20-30 reported as converged, or a zero-ordinate stop.
  # N(.3, sd) + gamma(shape, 1) at .3 + k * sd, as in the review's probe.
  # References: 30-digit mpmath tanh-sinh quadrature over dense breakpoints
  # (gamma density written from its formula); with sd = 1e-11 the value's
  # rounding error is 5e-6 of sd, so that case is checked at the criterion.
  cases <- list(
    list(shape = .8, sd = 1e-8, k = 0, reference = 19.963935864753457),
    list(shape = .8, sd = 1e-8, k = .3, reference = 23.795313625905310),
    list(shape = .8, sd = 1e-8, k = 1, reference = 29.555247422099176),
    list(shape = .8, sd = 1e-8, k = 3, reference = 27.884821065620117),
    list(shape = .3, sd = 1e-8, k = 0, reference = 183208.23150700132)
  )
  # the cases share their densities: one per shape and SD, built once (the
  # numerical grid of a convolution with a narrow Gaussian peak is the costly
  # part, and the ordinates read the same density at every k)
  convolutions <- list()
  for(case in cases){
    key <- paste(case$shape, case$sd)
    if(is.null(convolutions[[key]])){
      convolutions[[key]] <- density(
        list(a = prior("normal", list(.3, case$sd)), b = prior("gamma", list(case$shape, 1))),
        c(a = 1, b = 1)
      )
    }
    expect_equal(height(convolutions[[key]], .3 + case$k * case$sd), case$reference, tolerance = 1e-8)
  }
  convolution <- density(list(a = prior("normal", list(.3, 1e-11)), b = prior("gamma", list(.5, 1))),
                         c(a = 1, b = 1))
  reference <- 153441.73487990270
  expect_lt(abs(height(convolution, .3) - reference), tolerance$relative * reference)
})

test_that("conditional-normal breakpoints merge near-coincident points but never a Gaussian peak", {

  height <- function(priors, weights, value){
    density <- .prior_linear_combination_density(priors, weights, n_grid = 4096)
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$method, "conditional_normal_mixture")
    expect_true(ordinate$provenance$integration$converged)
    as.numeric(.prior_linear_density_height(density, value))
  }

  # A Gaussian-peak point a rounding error from a support bound left a piece a
  # few ulps wide, and QUADPACK stopped on roundoff. The value is the image of
  # the bound plus one Gaussian SD, as in the review's grid; the references
  # (integrate() at rel.tol 1e-12 with stats:: densities over the whole
  # support) agree with the review's to 1e-12.
  value <- .3 + .001 * 3 + 1
  reference <- stats::integrate(function(u){
    stats::dnorm(value - .3 - .001 * u) * stats::dgamma(u, 2, 1)
  }, 1, 3, rel.tol = 1e-12)$value / diff(stats::pgamma(c(1, 3), 2, 1))
  expect_equal(reference, .24169258799295, tolerance = 1e-12)
  priors <- list(a = prior("normal", list(.3, 1)),
                 b = prior("gamma", list(2, 1), list(lower = 1, upper = 3)))
  expect_equal(height(priors, c(a = 1, b = .001), value), reference, tolerance = 1e-8)

  value <- .3 + .001 * 1 - .01
  reference <- stats::integrate(function(u){
    stats::dnorm(value - .3 - .001 * u, 0, .01) * stats::dbeta(u, 2, 5)
  }, 0, 1, rel.tol = 1e-12)$value
  expect_equal(reference, 25.9220102711423, tolerance = 1e-12)
  priors <- list(a = prior("normal", list(.3, .01)), b = prior("beta", list(2, 5)))
  expect_equal(height(priors, c(a = 1, b = .001), value), reference, tolerance = 1e-8)

  # A Gaussian peak narrower than the merge width 1e-9 * max(1, |u|) keeps its
  # breakpoints: merging them left a peak shoulder at the end of a wide piece
  # (-15.9%, reported as converged). With SD s / |w| = 1e-10 the ordinate is
  # g(1) / |w| + O((s / w)^2) for the gamma(3, 2) density g.
  priors <- list(a = prior("normal", list(0, 1e-10)), b = prior("gamma", list(3, 2)))
  expect_equal(height(priors, c(a = 1, b = 1), 1), stats::dgamma(1, 3, 2), tolerance = 1e-6)
  priors <- list(a = prior("normal", list(0, 1e-7)), b = prior("gamma", list(3, 2)))
  expect_equal(height(priors, c(a = 1, b = 1e3), 1e3), stats::dgamma(1, 3, 2) / 1e3, tolerance = 1e-6)
})

test_that("lockstep mixture grids add the components' absolute changes", {

  # Two grid components whose refinements change in opposite directions: the
  # signed change of the mixture height is 0 at every refinement, the
  # weighted absolute changes are .01, .0099, .000099.
  sequences <- list(
    a = c(1, 1.01, 1.0001, 1.000001, 1.00000001),
    b = c(1, .99, .9999, .999999, .99999999)
  )
  grid <- function(name, level){
    structure(list(name = name, level = level, density = list(x = c(-10, 10))),
              class = "prior_linear_density")
  }
  testthat::local_mocked_bindings(
    .prior_linear_density_grid_height = function(x, value) sequences[[x$name]][x$level + 1L],
    .prior_linear_density_refinement = function(x) grid(x$name, x$level + 1L)
  )
  height <- .prior_linear_mixture_lockstep_height(list(
    list(weight = .5, method = "grid", density = grid("a", 0L)),
    list(weight = .5, method = "grid", density = grid("b", 0L))
  ), 0)
  evaluation <- attr(height, "adaptive_evaluation")
  expect_identical(evaluation$refinements, 3L)
  expect_equal(evaluation$absolute_change, .5 * (1.0001 - 1.000001) + .5 * (.999999 - .9999))
  expect_equal(as.numeric(height), 1)
})

test_that("grid refinement accepts small heights only within the relative criterion", {

  # Refinements of a far-tail height whose changes (1e-15) are far below an
  # absolute floor but 10% of the height: never converged. Relative changes
  # below 1e-4 converge, and an exactly zero height is kept.
  grid <- function(name, level){
    structure(list(name = name, level = level, density = list(x = c(-10, 10))),
              class = "prior_linear_density")
  }
  sequences <- list(
    tail     = 1e-14 * c(1, 1.1, 1.2, 1.3, 1.4),
    relative = 1e-14 * c(1, 1 + 2e-5, 1 + 3e-5, 1 + 3e-5, 1 + 3e-5),
    zero     = rep(0, 5)
  )
  testthat::local_mocked_bindings(
    .prior_linear_density_refinement = function(x) grid(x$name, x$level + 1L)
  )
  refine <- function(name){
    .prior_linear_density_refine_grids(
      list(grid(name, 0L)), 1,
      evaluate = function(density) sequences[[density$name]][density$level + 1L]
    )
  }
  expect_false(refine("tail")$converged)
  relative <- refine("relative")
  expect_true(relative$converged)
  expect_identical(relative$refinements, 1L)
  expect_equal(relative$error_bound, 1e-4 * 1e-14 * (1 + 2e-5))
  zero <- refine("zero")
  expect_true(zero$converged)
  expect_identical(zero$total, 0)
})

test_that("row and leaf conditional-normal ordinates each receive the full evaluation budget", {

  # Normal intercept + x * slope * sigma with sigma ~ half-normal. At 0 row r
  # has the closed-form ordinate exp(z) K0(z) / (2 pi |x_r|), z = 1 / (4 x_r^2)
  # (s = sinh(t) gives the K0 integral, https://dlmf.nist.gov/10.32.E9), and
  # the pooled rows are their average. The Gauss-Kronrod rules resolve these
  # smooth integrands to about 1e-10 relative; the 1e-4 acceptance gate is
  # QUADPACK's conservative error estimate.
  sigma <- prior("normal", list(0, 1), list(0, Inf))
  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "sigma"
  priors <- list(mu_intercept = prior("normal", list(0, 1)), mu_x = slope, sigma = sigma)
  row_ordinate <- function(x){
    besselK(1 / (4 * x^2), 0, expon.scaled = TRUE) / (2 * pi * abs(x))
  }
  # a change-of-variable reference away from 0 (sigma = tan(angle))
  row_reference <- function(x, value){
    stats::integrate(function(angle){
      s <- tan(angle)
      stats::dnorm(value, sd = sqrt(1 + x^2 * s^2)) * 2 * stats::dnorm(s) / cos(angle)^2
    }, 0, pi / 2, rel.tol = 1e-12)$value
  }
  row_ordinate_record <- function(context, rows, value){
    .prior_density_ordinate_from_adaptive(
      list(kind = "density_context_rows",
           arguments = list(context = context, weights = rows)),
      value
    )
  }
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

  # 150 distinct covariate rows with the default marginal_posterior() budget
  budget <- 10000L
  x <- stats::qnorm(seq(.005, .995, length.out = 150L))
  rows <- cbind(mu_intercept = 1, mu_x = x, sigma = 0)
  context <- .prior_density_build_context(priors, names(priors), n_grid = budget)
  result <- row_ordinate_record(context, rows, 0)
  expect_identical(result$provenance$unique_rows, 150L)
  expect_equal(exp(result$log_density), mean(row_ordinate(x)), tolerance = 1e-8)
  expect_length(budgets, 150L)
  expect_true(all(budgets == budget))
  expect_true(all(evaluations <= budget))
  row_integration <- lapply(result$provenance$row_classifications, `[[`, "integration")
  expect_true(all(vapply(row_integration, `[[`, logical(1), "converged")))
  expect_true(all(vapply(row_integration, `[[`, numeric(1), "budget") == budget))
  integration <- .prior_density_ordinate_integration(result$provenance)
  expect_true(integration$converged)
  expect_identical(integration$budget, 150 * budget)
  expect_identical(integration$evaluations, sum(evaluations))

  shifted <- row_ordinate_record(context, rows, .3)
  expect_equal(exp(shifted$log_density),
               mean(vapply(x, row_reference, numeric(1), value = .3)),
               tolerance = 1e-8)

  # a spike-and-slab slope: each of the 60 rows has an exact spike leaf and a
  # slab leaf whose quadrature gets the full budget
  spike_slope <- prior_spike_and_slab(prior("normal", list(0, 1)))
  attr(spike_slope, "multiply_by") <- "sigma"
  spike_priors <- priors
  spike_priors$mu_x <- spike_slope
  x_60 <- stats::qnorm(seq(.01, .99, length.out = 60L))
  rows_60 <- cbind(mu_intercept = 1, mu_x = x_60, sigma = 0)
  spike_context <- .prior_density_build_context(spike_priors, names(spike_priors), n_grid = budget)
  budgets <- evaluations <- numeric()
  spike_result <- row_ordinate_record(spike_context, rows_60, 0)
  expect_equal(exp(spike_result$log_density),
               mean(.5 * row_ordinate(x_60) + .5 * stats::dnorm(0)),
               tolerance = 1e-8)
  expect_length(budgets, 60L)
  expect_true(all(budgets == budget))
  expect_true(.prior_density_ordinate_integration(spike_result$provenance)$converged)

  # the public density height no longer rejects the pooled rows
  density <- .prior_density_from_context_rows(
    .prior_density_build_context(priors, names(priors), n_grid = 512L),
    rows
  )
  height <- .prior_linear_density_height(density, 0)
  expect_equal(as.numeric(height), mean(row_ordinate(x)), tolerance = 1e-8)
  expect_true(attr(height, "numerical_diagnostics")$converged)
})

test_that("mixture prior ordinates evaluate every component exactly at a density jump", {

  # N(0.5, 1)T(0, Inf) jumps at 0; its ordinate there is the one-sided limit
  # inside the support. N(0, 1) + N(0.5, 1)T(0, Inf) has the closed form below
  # (complete the square in the convolution integral).
  truncated <- prior("normal", list(.5, 1), list(0, Inf))
  f_truncated <- function(v) ifelse(v >= 0, stats::dnorm(v, .5) / stats::pnorm(.5), 0)
  f_sum <- function(v){
    stats::dnorm(v, .5, sqrt(2)) * stats::pnorm((v + .5) / sqrt(2)) / stats::pnorm(.5)
  }
  expect_height <- function(density, reference){
    for(value in c(0, -.5, .05)){
      ordinate <- prior_density_ordinate(density, value)
      expect_identical(ordinate$behavior, "regular")
      expect_true(ordinate$exact)
      expect_equal(as.numeric(.prior_linear_density_height(density, value)),
                   reference(value), tolerance = 1e-10)
    }
  }
  mock_fit <- function(samples, prior_list){
    fit <- structure(
      list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples),
           summary.pars = list(mutate = NULL), monitor = colnames(samples)),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(fit, "prior_list") <- prior_list
    fit <- attach_test_parameter_map(fit)
    fit
  }
  set.seed(71)
  n <- 200L
  data <- data.frame(t = factor(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi")))

  # Model mixture (mix_posteriors): the reference level is N(0, 1) or the
  # truncated intercept; the 'mid' level adds a treatment level that is
  # truncated in the first model and fixed at 0 in the second.
  first <- JAGS_formula(~ 1 + t, "mu", data = data, prior_list = list(
    intercept = prior("normal", list(0, 1)),
    t = prior_factor("normal", list(.5, 1), list(0, Inf), contrast = "treatment")
  ))$prior_list
  second <- JAGS_formula(~ 1 + t, "mu", data = data, prior_list = list(
    intercept = truncated,
    t = prior_factor("point", list(location = 0), contrast = "treatment")
  ))$prior_list
  models <- list(
    list(fit = mock_fit(cbind(mu_intercept = stats::rnorm(n), "mu_t[1]" = rng(truncated, n),
                              "mu_t[2]" = rng(truncated, n)), first),
         marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = mock_fit(cbind(mu_intercept = rng(truncated, n), "mu_t[1]" = 0,
                              "mu_t[2]" = 0), second),
         marglik = bridgesampling_object(0), prior_weights = 1)
  )
  mixed <- mix_posteriors(models, parameters = c("mu_intercept", "mu_t"),
                          is_null_list = list(mu_intercept = c(FALSE, FALSE),
                                              mu_t = c(FALSE, FALSE)),
                          seed = 1, n_samples = n)
  levels <- marginal_posterior(mixed, "mu_t", formula = ~ 1 + t, prior_samples = TRUE)
  reference_level <- .bt_meta_get(levels[["lo"]], "prior_density")
  expect_equal(as.numeric(.prior_linear_density_height(reference_level, 0)),
               .5 * stats::dnorm(0) + .5 * stats::dnorm(0, .5) / stats::pnorm(.5),
               tolerance = 1e-10)
  expect_height(reference_level, function(v) .5 * stats::dnorm(v) + .5 * f_truncated(v))
  # the already exact per-model ordinates are unchanged
  expect_identical(
    vapply(prior_density_ordinate(reference_level, 0)$provenance$components,
           function(component) component$provenance$kind, character(1)),
    c("scalar_affine", "scalar_affine")
  )
  two_term <- .bt_meta_get(levels[["mid"]], "prior_density")
  expect_height(two_term, function(v) .5 * f_sum(v) + .5 * f_truncated(v))
  # the normal + truncated normal component is the closed-form convolution
  closed_form <- prior_density_ordinate(two_term, 0)$provenance$components[[1L]]$provenance
  expect_identical(closed_form$kind, "truncated_normal_convolution")
  expect_null(closed_form$integration)

  # Single fit (as_mixed_posteriors): mixture intercept and a spike-or-normal
  # slope; each indicator combination is its own component.
  single <- JAGS_formula(~ x, "mu", data = data.frame(x = c(-1, 0, 1)), prior_list = list(
    intercept = prior_mixture(list(prior("normal", list(0, 1)), truncated),
                              is_null = c(FALSE, FALSE)),
    x = prior_mixture(list(prior("spike", list(0)), prior("normal", list(0, 1))),
                      is_null = c(TRUE, FALSE))
  ))
  intercept_indicator <- sample(1:2, n, TRUE)
  slope_indicator <- sample(1:2, n, TRUE)
  fit <- coda::mcmc(cbind(
    mu_intercept = ifelse(intercept_indicator == 1L, stats::rnorm(n), rng(truncated, n)),
    mu_x = ifelse(slope_indicator == 1L, 0, stats::rnorm(n)),
    mu_intercept_indicator = intercept_indicator,
    mu_x_indicator = slope_indicator
  ))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- single$prior_list
  fit <- attach_test_parameter_map(fit)
  single_mixed <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x"))
  single_levels <- marginal_posterior(single_mixed, "mu_x", formula = ~ x, prior_samples = TRUE)
  expect_height(.bt_meta_get(single_levels[["1SD"]], "prior_density"), function(v){
    .25 * (stats::dnorm(v) + stats::dnorm(v, sd = sqrt(2)) + f_truncated(v) + f_sum(v))
  })
  expect_height(.bt_meta_get(single_levels[["0SD"]], "prior_density"), function(v){
    .5 * stats::dnorm(v) + .5 * f_truncated(v)
  })
  expect_identical(
    prior_density_ordinate(.bt_meta_get(single_levels[["0SD"]], "prior_density"), 0)$method,
    "scalar_affine"
  )

  # Conditional mixture: conditioning on the slope's alternative leaves the
  # N(0, 1) + N(0, 1) and truncated + N(0, 1) components.
  context <- .prior_density_build_context(single$prior_list, c("mu_intercept", "mu_x"),
                                          conditional = "mu_x", n_grid = 4096)
  expect_s3_class(context, "prior_density_conditional_context")
  expect_height(.prior_density_from_context(context, c(mu_intercept = 1, mu_x = 1)),
                function(v) .5 * stats::dnorm(v, sd = sqrt(2)) + .5 * f_sum(v))

  # Original-scale intercept under transform_scaled: b0 - (m / s) * b1 with the
  # intercept mixture above and a spike-and-slab standardized slope. Reference:
  # the four component densities, the bounded + normal one by a 1-D
  # convolution integral at rel.tol 1e-10.
  scaled <- JAGS_formula(~ x, "mu", data = data.frame(x = c(1, 2, 3.5, 4, 6, 8.5)), prior_list = list(
    intercept = prior_mixture(list(prior("normal", list(0, 1)), truncated),
                              is_null = c(FALSE, FALSE)),
    x = prior_spike_and_slab(prior("normal", list(0, 1)))
  ), formula_scale = list(x = TRUE))
  ratio <- scaled$formula_scale[["mu_x"]]$mean / scaled$formula_scale[["mu_x"]]$sd
  scaled_indicator <- sample(1:2, n, TRUE)
  slope_inclusion <- stats::rbinom(n, 1, .5)
  posterior <- cbind(
    mu_intercept = ifelse(scaled_indicator == 1L, stats::rnorm(n), rng(truncated, n)),
    mu_x = slope_inclusion * stats::rnorm(n),
    mu_intercept_indicator = scaled_indicator,
    mu_x_indicator = slope_inclusion
  )
  scaled_fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
         monitor = colnames(posterior), sample = n),
    class = c("runjags", "BayesTools_fit")
  )
  attr(scaled_fit, "prior_list") <- scaled$prior_list
  attr(scaled_fit, "formula_design") <- list(mu = scaled$formula_design)
  attr(scaled_fit, "formula_scale") <- list(mu = scaled$formula_scale)
  scaled_fit <- attach_test_parameter_map(scaled_fit)
  scaled_fit <- .bt_attach_fit_contract(.bt_attach_draw_geometry(.bt_attach_parameter_map(scaled_fit)))
  scaled_mixed <- as_mixed_posteriors(scaled_fit, c("mu_intercept", "mu_x"),
                                      transform_scaled = TRUE, n_prior_samples = 2000)
  original_intercept <- marginal_posterior(scaled_mixed, "mu_intercept", use_formula = FALSE,
                                           prior_samples = TRUE)
  expect_equal(as.numeric(original_intercept),
               posterior[, "mu_intercept"] - ratio * posterior[, "mu_x"])
  expect_height(.bt_meta_get(original_intercept, "prior_density"), function(v){
    .25 * stats::dnorm(v) + .25 * stats::dnorm(v, sd = sqrt(1 + ratio^2)) +
      .25 * f_truncated(v) + .25 * stats::integrate(function(t){
        f_truncated(t) * stats::dnorm(v - t, sd = ratio)
      }, 0, Inf, rel.tol = 1e-10)$value
  })

  # Outside mixtures the same convolution uses the same closed form.
  plain <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = truncated), c(a = 1, b = 1)
  )
  expect_identical(prior_density_ordinate(plain, .05)$method, "truncated_normal_convolution")
  expect_height(plain, f_sum)
})

test_that("mixture components without exact ordinates are refined on their own grids", {

  # An intercept mixture (N(0, 1) | T = N(0.5, 1)T(0, Inf)) and a treatment
  # mixture (0 | T for both coefficients). The densities of T, of T + T,
  # f2(v) = phi((v - 1) / sqrt(2)) (2 Phi(v / sqrt(2)) - 1) / (sqrt(2) Phi(0.5)^2),
  # and of N + T are closed forms; three-term sums by quadrature over f2.
  truncated <- prior("normal", list(.5, 1), list(0, Inf))
  f_truncated <- function(v) ifelse(v >= 0, stats::dnorm(v, .5) / stats::pnorm(.5), 0)
  f_sum <- function(v){
    stats::dnorm(v, .5, sqrt(2)) * stats::pnorm((v + .5) / sqrt(2)) / stats::pnorm(.5)
  }
  f_double <- function(v){
    ifelse(v > 0, stats::dnorm(v, 1, sqrt(2)) * (2 * stats::pnorm(v / sqrt(2)) - 1) /
             stats::pnorm(.5)^2, 0)
  }
  f_normal_double <- function(v){
    stats::integrate(function(u) stats::dnorm(v - u) * f_double(u), 0, Inf, rel.tol = 1e-12)$value
  }
  f_triple <- function(v){
    if(v <= 0) return(0)
    stats::integrate(function(u) f_truncated(u) * f_double(v - u), 0, v, rel.tol = 1e-12)$value
  }
  data <- data.frame(t = factor(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi")))
  priors <- JAGS_formula(~ 1 + t, "mu", data = data, prior_list = list(
    intercept = prior_mixture(list(prior("normal", list(0, 1)), truncated),
                              is_null = c(FALSE, FALSE)),
    t = prior_mixture(list(prior_factor("point", list(location = 0), contrast = "treatment"),
                           prior_factor("normal", list(.5, 1), list(0, Inf), contrast = "treatment")),
                      is_null = c(TRUE, FALSE))
  ))$prior_list
  context <- .prior_density_build_context(priors, c("mu_intercept", "mu_t[1]", "mu_t[2]"),
                                          n_grid = 10000)

  # the level intercept + t[mid] has components N, N + T (Gaussian
  # convolution), T (jump at 0) and T + T (two-term convolution): every
  # component is structural, also at the jump and next to the T + T kink
  level <- .prior_density_from_context(context, c(mu_intercept = 1, "mu_t[1]" = 1, "mu_t[2]" = 0))
  for(value in c(-.05, -.01, 0, .05, 1)){
    height <- .prior_linear_density_height(level, value)
    expect_null(attr(height, "adaptive_evaluation"))
    expect_equal(as.numeric(height),
                 .25 * (stats::dnorm(value) + f_sum(value) + f_truncated(value) + f_double(value)),
                 tolerance = 1e-10)
  }

  # intercept + t[mid] + t[hi] has components N, N + T + T, T and T + T + T;
  # the three-term sums have no structural ordinate and are evaluated on
  # their own grids, refined in lockstep with the documented criterion on the
  # mixture height
  density <- .prior_density_from_context(context, c(mu_intercept = 1, "mu_t[1]" = 1, "mu_t[2]" = 1))
  reference <- function(v){
    .25 * (stats::dnorm(v) + f_normal_double(v) + f_truncated(v) + f_triple(v))
  }
  for(value in c(-.05, 0, .05, 1)){
    height <- .prior_linear_density_height(density, value)
    evaluation <- attr(height, "adaptive_evaluation")
    expect_true(evaluation$converged)
    # T + T + T is exactly zero below its support
    grids <- if(value < 0) 1L else 2L
    expect_identical(evaluation$components, grids)
    components <- attr(height, "component_heights")
    expect_identical(vapply(components, `[[`, character(1), "method"),
                     c("exact", "grid", "exact", if(value < 0) "exact" else "grid"))
    expect_equal(vapply(components, `[[`, numeric(1), "weight"), rep(.25, 4))
    expect_lte(abs(as.numeric(height) - reference(value)), evaluation$error_bound)
  }
})

test_that("model-averaged mean-difference coordinates are exact mixtures of their components", {

  # A RoBMA-style product-space factor prior: the mean-difference slab
  # mNormal(0, .35) and the spike at 0, each with probability 1/2. Every
  # coordinate of the slab is N(0, .35), so the combination with weights d is
  # N(0, .35 ||d||) in the slab and the atom at 0 in the spike: its continuous
  # density is .5 * dnorm(x, 0, .35 ||d||), also for a single coordinate (a
  # level whose mean-difference design row has one nonzero entry), which is
  # not a scalar term.
  data <- data.frame(g = factor(c("a", "b", "c"), levels = c("a", "b", "c")))
  priors <- JAGS_formula(~ 1 + g, "mu", data = data, prior_list = list(
    intercept = prior("normal", list(0, 1)),
    g = prior_mixture(list(prior_factor("mnormal", list(0, .35), contrast = "meandif"),
                           prior_factor("spike", list(0), contrast = "meandif")),
                      is_null = c(FALSE, TRUE))
  ))$prior_list
  context <- .prior_density_build_context(priors, c("mu_intercept", "mu_g[1]", "mu_g[2]"))
  design <- contr.meandif(3)
  weight_sets <- c(
    list(c("mu_g[1]" = 1), c("mu_g[2]" = 1), c("mu_g[2]" = -2)),
    lapply(1:3, function(level) c("mu_g[1]" = design[level, 1], "mu_g[2]" = design[level, 2]))
  )
  for(weights in weight_sets){
    density <- .prior_density_from_context(context, weights)
    sd <- .35 * sqrt(sum(weights^2))
    for(value in c(-.4, .1, 1)){
      ordinate <- prior_density_ordinate(density, value)
      expect_true(ordinate$exact)
      expect_identical(ordinate$behavior, "regular")
      expect_identical(ordinate$method, "finite_mixture")
      expect_equal(exp(ordinate$log_density), .5 * stats::dnorm(value, 0, sd), tolerance = 1e-14)
      expect_equal(as.numeric(.prior_linear_density_height(density, value)),
                   .5 * stats::dnorm(value, 0, sd), tolerance = 1e-14)
    }
    zero <- prior_density_ordinate(density, 0)
    expect_identical(zero$behavior, "point_mass")
    expect_equal(zero$point_mass, .5)
    expect_identical(zero$provenance$continuous_behavior, "regular")
    status <- prior_ordinate_status(density, c(0, .1))
    expect_identical(status$eligible, c(FALSE, TRUE))
    expect_identical(status$condition[[1L]], "BayesTools_point_mass_at_null")
    probability <- .hypothesis_prior_density_prob(
      density, hypothesis_parse("theta > 0.1")$statements[[1L]]$left, "theta"
    )
    expect_equal(as.numeric(probability), .5 * stats::pnorm(.1, 0, sd, lower.tail = FALSE),
                 tolerance = 1e-14)
  }
})

test_that("linear combinations of multivariate t priors are univariate t terms", {

  # A multivariate t vector prior X = mu 1 + z / sqrt(w), z ~ N(0, s^2 I),
  # w ~ Gamma(nu / 2, nu / 2) (the scale-matrix parameterization of the JAGS
  # emitter, checked below) has a'X ~ t(mu sum(a), s ||a||, nu). References:
  # the closed-form t density and distribution function (extraDistr), the
  # Cauchy and Gaussian-convolution integrals by integrate() at rel.tol 1e-12,
  # and a Monte Carlo sample of 1e6 draws of X built as the emitter builds it,
  # compared by bins and regions within 4 binomial standard errors. The bin
  # probability of an ordinate is Simpson's rule over the bin (width .1; its
  # error is below 1e-8, against standard errors above 5e-5). Before, every
  # such combination was a general convolution: 'unknown', with grid values
  # (region probabilities without a structural route).
  syntax <- .JAGS_prior.vector(prior("mt", list(location = 2, scale = .5, df = 3, K = 3)), "p")
  expect_match(syntax, "prior_par_s_p ~ dgamma(1.5, 1.5)", fixed = TRUE)
  expect_match(syntax, "prior_par2_p[i,i] <- 4", fixed = TRUE)
  expect_match(syntax, "p[i] <- prior_par_z_p[i]/sqrt(prior_par_s_p) + 2", fixed = TRUE)

  draw_mt <- function(n, K, location, scale, df){
    z <- matrix(stats::rnorm(n * K, 0, scale), nrow = n, ncol = K)
    location + z / sqrt(stats::rgamma(n, shape = df / 2, rate = df / 2))
  }
  region_probability <- function(density, lower, upper){
    .prior_linear_density_region_probability(density, list(
      intervals = .prior_region_intervals(lower, upper),
      indicator = function(values) values > lower & values < upper
    ))
  }
  continuous_density <- function(density, value){
    exp(prior_density_ordinate(density, value)$log_density)
  }
  integrate_pieces <- function(f, breaks){
    sum(vapply(seq_len(length(breaks) - 1L), function(i){
      stats::integrate(f, breaks[i], breaks[i + 1L], rel.tol = 1e-12)$value
    }, numeric(1)))
  }
  expect_monte_carlo <- function(density, draws, bins, regions, atom = 0){
    n <- length(draws)
    for(center in bins){
      edges <- center + c(-.05, 0, .05)
      heights <- vapply(edges, continuous_density, numeric(1), density = density)
      expected <- (heights[1L] + 4 * heights[2L] + heights[3L]) * .1 / 6
      observed <- mean(draws > edges[1L] & draws < edges[3L])
      expect_lte(abs(observed - expected), 4 * sqrt(expected * (1 - expected) / n))
    }
    for(lower in regions){
      expected <- as.numeric(region_probability(density, lower, Inf))
      observed <- mean(draws > lower)
      expect_lte(abs(observed - expected), 4 * sqrt(expected * (1 - expected) / n))
    }
  }

  data <- data.frame(g = factor(rep(letters[1:4], 3)), h = factor(rep(letters[1:3], each = 4)))
  columns <- c("mu_intercept", paste0("mu_g[", 1:3, "]"), paste0("mu_h[", 1:2, "]"))
  for(contrast in c("meandif", "orthonormal")){
    priors <- JAGS_formula(~ 1 + g + h, "mu", data = data, prior_list = list(
      intercept = prior("normal", list(0, 1)),
      g = prior_factor("mt", list(location = 0, scale = .5, df = 3), contrast = contrast),
      h = prior_factor("mcauchy", list(location = 0, scale = .25), contrast = contrast)
    ))$prior_list
    context <- .prior_density_build_context(priors, columns)
    design <- if(contrast == "meandif") contr.meandif(4) else contr.orthonormal(4)
    level <- function(i) stats::setNames(design[i, ], paste0("mu_g[", 1:3, "]"))
    weight_sets <- list(c("mu_g[2]" = 1), level(1), level(4), level(2) - level(3),
                        c("mu_g[1]" = .3, "mu_g[2]" = -.7, "mu_g[3]" = 2))
    # every check of the contrast is made; the failures are collected and
    # asserted once per contrast
    problems <- expectation_problems({
      for(weights in weight_sets){
        density <- .prior_density_from_context(context, weights)
        scale <- .5 * sqrt(sum(weights^2))
        for(value in c(-1.3, 0, .2, 25)){
          ordinate <- prior_density_ordinate(density, value)
          expect_true(ordinate$exact)
          expect_identical(ordinate$behavior, "regular")
          expect_identical(ordinate$method, "scalar_affine")
          expect_equal(ordinate$log_density,
                       extraDistr::dlst(value, df = 3, mu = 0, sigma = scale, log = TRUE),
                       tolerance = 1e-13)
        }
        record <- ordinate$provenance$multivariate_t[[1L]]
        expect_identical(record$parameter, "mu_g")
        expect_equal(record$weights, weights[weights != 0])
        expect_equal(record$t, c(location = 0, scale = scale, df = 3), tolerance = 1e-15)
        expect_equal(as.numeric(.prior_linear_density_height(density, .2)),
                     extraDistr::dlst(.2, df = 3, mu = 0, sigma = scale), tolerance = 1e-13)
        probability <- region_probability(density, -.4, .9)
        expect_identical(attr(probability, "numerical_diagnostics")$method, "exact")
        expect_equal(as.numeric(probability),
                     diff(extraDistr::plst(c(-.4, .9), df = 3, mu = 0, sigma = scale)),
                     tolerance = 1e-13)
        curve <- .prior_linear_density_to_plot_data(density, n_points = 51, x_range = c(-3, 3))$density
        expect_equal(curve$y, extraDistr::dlst(curve$x, df = 3, mu = 0, sigma = scale),
                     tolerance = 1e-13)
        expect_identical(prior_ordinate_status(density, .2)$eligible, TRUE)
      }

      # a Cauchy level (mcauchy, one degree of freedom) plus the normal
      # intercept is a Gaussian convolution, and a mt level plus a Cauchy level
      # a two-term convolution
      h_level <- stats::setNames(if(contrast == "meandif") contr.meandif(3)[1, ] else contr.orthonormal(3)[1, ],
                                 paste0("mu_h[", 1:2, "]"))
      h_scale <- .25 * sqrt(sum(h_level^2))
      combinations <- list(
        list(weights = c(mu_intercept = 1, h_level), method = "conditional_normal_mixture",
             reference = function(value){
               integrate_pieces(function(x) stats::dcauchy(x, 0, h_scale) * stats::dnorm(value - x),
                                c(-Inf, sort(c(0, value)), Inf))
             }),
        list(weights = c(level(1), h_level), method = "convolution",
             reference = function(value){
               integrate_pieces(function(x){
                 stats::dcauchy(x, 0, h_scale) *
                   extraDistr::dlst(value - x, df = 3, mu = 0, sigma = .5 * sqrt(sum(level(1)^2)))
               }, c(-Inf, sort(c(0, value)), Inf))
             })
      )
      for(combination in combinations){
        density <- .prior_density_from_context(context, combination$weights)
        for(value in c(-.8, .3, 4)){
          ordinate <- prior_density_ordinate(density, value)
          expect_true(ordinate$exact)
          expect_identical(ordinate$method, combination$method)
          expect_equal(exp(ordinate$log_density), combination$reference(value), tolerance = 1e-6)
        }
      }
    })
    expect_identical(problems, character(), info = contrast)
  }

  # a non-factor vector prior with a nonzero location: the location is
  # mu sum(a); two Cauchy vector terms sum to one Cauchy term; the zero
  # combination is the point at 0
  vector_priors <- list(p = prior("mt", list(location = 2, scale = .5, df = 5, K = 3)),
                        q = prior("mcauchy", list(location = -1, scale = .5, K = 2)),
                        r = prior("mcauchy", list(location = 1, scale = .25, K = 2)))
  weights <- c("p[1]" = .5, "p[2]" = 1, "p[3]" = -2)
  density <- .prior_linear_combination_density(vector_priors, weights)
  expect_equal(prior_density_ordinate(density, 1)$log_density,
               extraDistr::dlst(1, df = 5, mu = -1, sigma = .5 * sqrt(5.25), log = TRUE),
               tolerance = 1e-13)
  cauchy <- .prior_linear_combination_density(vector_priors, c("q[1]" = 1, "q[2]" = 1, "r[2]" = -2))
  ordinate <- prior_density_ordinate(cauchy, .5)
  expect_true(ordinate$exact)
  expect_identical(ordinate$method, "scalar_affine")
  expect_equal(exp(ordinate$log_density), stats::dcauchy(.5, -4, .5 * sqrt(2) + .5), tolerance = 1e-13)
  expect_length(ordinate$provenance$multivariate_t, 2L)
  zero <- prior_density_ordinate(.prior_linear_combination_density(vector_priors, c("p[1]" = 0)), 0)
  expect_identical(zero$behavior, "point_mass")
  expect_identical(zero$point_mass, 1)

  # a truncation set after construction (vector priors do not support one) is
  # not dropped: the group keeps the general route, as a truncated mnormal
  # keeps it (the untruncated t was exact at -1, outside the truncation)
  truncated <- vector_priors["p"]
  truncated$p$truncation <- list(lower = 0, upper = Inf)
  ordinate <- prior_density_ordinate(.prior_linear_combination_density(truncated, weights), -1)
  expect_false(ordinate$exact)
  expect_null(ordinate$provenance$multivariate_t)

  # a 'multiply_by' scale: the product of the t combination and the scale
  scaled <- vector_priors["p"]
  attr(scaled$p, "multiply_by") <- "s"
  scaled$s <- prior("lognormal", list(0, .5))
  product <- .prior_linear_combination_density(scaled, weights)
  for(value in c(-3, .5, 6)){
    ordinate <- prior_density_ordinate(product, value)
    expect_true(ordinate$exact)
    expect_identical(ordinate$method, "scale_mixture")
    reference <- integrate_pieces(function(s){
      extraDistr::dlst(value / s, df = 5, mu = -1, sigma = .5 * sqrt(5.25)) / s * stats::dlnorm(s, 0, .5)
    }, sort(unique(c(0, 1, abs(value), Inf))))
    expect_equal(exp(ordinate$log_density), reference, tolerance = 1e-6)
  }

  # Monte Carlo of the multivariate t: a level of a mt factor prior, the
  # normal intercept plus a contrast, and the mixture rule (a model-averaged
  # mixture with a spike, and a spike-and-slab with inclusion Beta(2, 3))
  set.seed(20260926)
  n <- 1e6
  x <- draw_mt(n, 3, 0, .5, 3)
  priors <- JAGS_formula(~ 1 + g, "mu", data = data, prior_list = list(
    intercept = prior("normal", list(0, 1)),
    g = prior_factor("mt", list(location = 0, scale = .5, df = 3), contrast = "meandif")
  ))$prior_list
  context <- .prior_density_build_context(priors, columns[1:4])
  design <- contr.meandif(4)
  expect_monte_carlo(.prior_density_from_context(context, stats::setNames(design[3, ], columns[2:4])),
                     draws = drop(x %*% design[3, ]), bins = c(-1.2, -.3, 0, .45, 1.6),
                     regions = c(-.9, .1, 2))
  contrast_weights <- design[1, ] - design[4, ]
  expect_monte_carlo(.prior_density_from_context(context, stats::setNames(c(1, contrast_weights), columns[1:4])),
                     draws = drop(x %*% contrast_weights) + stats::rnorm(n),
                     bins = c(-2.5, -.6, .1, 1.3), regions = c(-1.5, .4, 3))

  for(kind in c("mixture", "spike_and_slab")){
    slab <- prior_factor("mt", list(location = 0, scale = .5, df = 3), contrast = "orthonormal")
    g <- if(kind == "mixture"){
      prior_mixture(list(slab, prior_factor("spike", list(0), contrast = "orthonormal")),
                    is_null = c(FALSE, TRUE))
    }else{
      prior_spike_and_slab(slab, prior_inclusion = prior("beta", list(2, 3)))
    }
    inclusion <- if(kind == "mixture") .5 else .4
    priors <- JAGS_formula(~ 1 + g, "mu", data = data, prior_list = list(
      intercept = prior("normal", list(0, 1)), g = g
    ))$prior_list
    context <- .prior_density_build_context(priors, columns[1:4])
    weights <- contr.orthonormal(4)[2, ]
    density <- .prior_density_from_context(context, stats::setNames(weights, columns[2:4]))
    scale <- .5 * sqrt(sum(weights^2))
    ordinate <- prior_density_ordinate(density, .3)
    expect_true(ordinate$exact)
    expect_identical(ordinate$method, "finite_mixture")
    expect_equal(exp(ordinate$log_density),
                 inclusion * extraDistr::dlst(.3, df = 3, mu = 0, sigma = scale), tolerance = 1e-13)
    zero <- prior_density_ordinate(density, 0)
    expect_identical(zero$behavior, "point_mass")
    expect_equal(zero$point_mass, 1 - inclusion, tolerance = 1e-15)
    expect_identical(zero$provenance$continuous_behavior, "regular")
    draws <- drop(x %*% weights) * (stats::runif(n) < inclusion)
    expect_monte_carlo(density, draws, bins = c(-1.1, -.2, .35, 1.5), regions = c(-.7, .15, 1.8))
  }
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
  fit <- attach_test_parameter_map(fit)
  fit <- .bt_attach_parameter_map(fit, monitor_names = "p")
  marginal <- marginal_posterior(
    as_mixed_posteriors(fit, "p"), "p", use_formula = FALSE,
    prior_samples = TRUE, n_samples = 128
  )
  expect_equal(
    .prior_linear_density_height(.bt_meta_get(marginal, "prior_density"), .5),
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
  # 10 when polynomial tails would more than double the range.
  # t3 + t30(0, .1) + t30(0, .1) converges to its nested-quadrature reference;
  # Cauchy combinations exceed the grid limit before converging and stop
  # loudly instead of reporting a height biased by their omitted tail mass.
  # (A normal term plus one other term is a Gaussian convolution and two
  # non-normal terms a two-term convolution, both evaluated by quadrature, so
  # these grid vehicles have three terms.) A spike-and-slab mixture is
  # evaluated per component (the slab component is a Gaussian convolution by
  # quadrature), so it matches the reference above.
  student <- .prior_linear_combination_density(
    list(a = prior("t", list(0, 1, 3)), b = prior("t", list(0, .1, 30)),
         c = prior("t", list(0, .1, 30))),
    c(a = 1, b = 1, c = 1)
  )
  student_height <- .prior_linear_density_height(student, 0)
  expect_true(isTRUE(attr(student_height, "adaptive_evaluation")$converged))
  narrow <- function(v) stats::dt(v / .1, 30) / .1
  narrow_sum <- function(u){
    vapply(u, function(value){
      stats::integrate(function(s) narrow(s) * narrow(value - s), -Inf, Inf, rel.tol = 1e-12)$value
    }, numeric(1))
  }
  expect_equal(
    as.numeric(student_height),
    stats::integrate(function(t) stats::dt(t, 3) * narrow_sum(-t), -Inf, Inf, rel.tol = 1e-10)$value,
    tolerance = 1e-4
  )
  cauchy_sum <- .prior_linear_combination_density(
    list(a = prior("t", list(0, .2, 30)), b = prior("cauchy", list(0, .707)),
         c = prior("t", list(0, .2, 30))),
    c(a = 1, b = 1, c = 1)
  )
  expect_error(
    .prior_linear_density_height(cauchy_sum, 0),
    "Adaptive prior-density evaluation did not converge within the documented grid-refinement error criterion.",
    fixed = TRUE
  )
  expect_equal(as.numeric(.prior_linear_density_height(heavy, 0)), reference, tolerance = 1e-8)

  # Density jumps and kinks (half-normal and uniform terms) converge once the
  # spacing strictly halves; references by quadrature over the closed-form
  # density of two half-normals, f2(x) = 4 phi(x / sqrt(2)) (2 Phi(x / sqrt(2))
  # - 1) / sqrt(2) on x > 0, and of two uniforms (triangular).
  half_normal <- prior("normal", list(0, 1), truncation = list(0, Inf))
  half_normal_sum <- function(x){
    ifelse(x > 0, 4 * stats::dnorm(x / sqrt(2)) * (2 * stats::pnorm(x / sqrt(2)) - 1) / sqrt(2), 0)
  }
  triangular <- function(x) ifelse(x < 0 | x > 2, 0, ifelse(x < 1, x, 2 - x))
  uniform_sum <- function(t) triangular(t) * stats::dt((-.2 - t) / .3, 30) / .3
  jumps <- list(
    list(priors = list(a = half_normal, b = half_normal, c = half_normal), value = 1.5,
         reference = stats::integrate(function(t) 2 * stats::dnorm(t) * half_normal_sum(1.5 - t),
                                      0, 1.5, rel.tol = 1e-12)$value),
    list(priors = list(a = half_normal, b = half_normal, c = prior("t", list(0, 1, 30))), value = -1,
         reference = stats::integrate(function(t) half_normal_sum(t) * stats::dt(-1 - t, 30),
                                      0, Inf, rel.tol = 1e-12)$value),
    list(priors = list(a = prior("uniform", list(0, 1)), b = prior("uniform", list(0, 1)),
                       c = prior("t", list(0, .3, 30))), value = -.2,
         reference = stats::integrate(uniform_sum, 0, 1, rel.tol = 1e-12)$value +
           stats::integrate(uniform_sum, 1, 2, rel.tol = 1e-12)$value)
  )
  for(jump in jumps){
    density <- .prior_linear_combination_density(jump$priors, c(a = 1, b = 1, c = 1))
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
  # model spacing needs about 3.3e6 knots: refinement of the mixture grid (used
  # for prior probabilities without a structural representation) ends and is
  # reported as not converged, instead of failing with the mixing error meant
  # for incompatible scales in the requested density itself. The height is the
  # weighted sum of the models' own ordinates (both Gaussian convolutions) and
  # never refines the mixture grid; reference by quadrature. So is the region
  # probability (reference: 30-digit mpmath, see
  # test-hypothesis-prior-region-grid.R).
  context <- .prior_density_build_context(
    list(a = list(prior("normal", list(0, .01), prior_weights = 1),
                  prior("t", list(0, 1, 3), prior_weights = 1)),
         b = list(prior("gamma", list(3, 2), prior_weights = 1),
                  prior("normal", list(0, 1), prior_weights = 1))),
    c("a", "b")
  )
  density <- .prior_density_from_context(context, c(a = 1, b = 1))
  expect_lt(length(density$density$x), .prior_linear_density_max_grid())
  side <- hypothesis_parse("theta < .3")$statements[[1L]]$left
  expect_error(
    .hypothesis_prior_density_grid_prob(
      density, .hypothesis_simple_parameter_comparison(side, "theta"),
      .hypothesis_side_expression(side), "theta"
    ),
    "Adaptive prior-probability evaluation did not converge within the documented grid-refinement error criterion.",
    fixed = TRUE
  )
  expect_lt(abs(.hypothesis_prior_density_prob(density, side, "theta") -
                  .2993480255493350221), 1e-8)
  expect_equal(
    as.numeric(.prior_linear_density_height(density, .3)),
    .5 * stats::integrate(function(t) stats::dnorm(.3 - t, 0, .01) * stats::dgamma(t, 3, 2),
                          0, Inf, rel.tol = 1e-12)$value +
      .5 * stats::integrate(function(t) stats::dt(t, 3) * stats::dnorm(.3 - t),
                            -Inf, Inf, rel.tol = 1e-12)$value,
    tolerance = 1e-8
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
  # (three terms: a normal plus one other term and two non-normal terms have
  # structural ordinates and would not reach the grid)
  attr(density, "adaptive_evaluation") <- list(
    kind = "linear_combination",
    arguments = list(prior_list = list(a = prior("t", list(0, 1, 30)),
                                       b = prior("uniform", list(0, 1)),
                                       c = prior("uniform", list(0, 1))),
                     weights = c(a = 1, b = 1, c = 1), n_grid = 4096, tail_prob = 1e-4)
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
      x = prior("t", list(0, 1, 30)),
      y = prior("gamma", list(3, 2)),
      z = prior("gamma", list(2, 2))
    ),
    weights   = c(x = 1, y = 1, z = 1),
    n_grid    = 512,
    tail_prob = 1e-3
  )
  # the two gamma terms sum to gamma(5, 2)
  skewed_reference <- function(value){
    stats::integrate(function(t) stats::dt(value - t, 30) * stats::dgamma(t, 5, 2),
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

  # A t30 plus two gamma terms has no structural ordinate and is refined (a
  # normal plus one other term would be a Gaussian convolution, two non-normal
  # terms a two-term convolution); the reference is the convolution integral.
  center <- BayesTools:::.prior_linear_density_height(skewed, 0)
  expect_identical(refinement_calls, 2L)
  # Each refinement halves the source spacing (2046 knots over [-3.39, 13.6]
  # initially) and omits 1000 times less tail probability.
  expect_equal(
    attr(center, "adaptive_evaluation")[c("n_grid", "tail_prob", "refinements")],
    list(n_grid = 32768L, tail_prob = 1e-9, refinements = 2L)
  )
  refinement_calls <- 0L
  tail <- BayesTools:::.prior_linear_density_height(skewed, 14)
  expect_gt(14, max(skewed$density$x))
  expect_identical(refinement_calls, 4L)
  expect_lt(
    abs(as.numeric(center) / skewed_reference(0) - 1),
    1e-4
  )
  expect_lt(
    abs(as.numeric(tail) / skewed_reference(14) - 1),
    1e-4
  )
  expect_true(isTRUE(attr(center, "adaptive_evaluation")$converged))
  expect_true(isTRUE(attr(tail, "adaptive_evaluation")$converged))
  expect_gt(attr(tail, "adaptive_evaluation")$refinements, 0)

  # The normal sum also has an exact region probability; its grid evaluation
  # (used for combinations without a structural representation) refines once.
  refinement_calls <- 0L
  side <- hypothesis_parse("theta < 0")$statements[[1L]]$left
  probability <- BayesTools:::.hypothesis_prior_density_prob(
    density, side, "theta"
  )
  expect_equal(probability, .5, tolerance = 1e-15)
  expect_identical(refinement_calls, 0L)
  probability <- BayesTools:::.hypothesis_prior_density_grid_prob(
    density, BayesTools:::.hypothesis_simple_parameter_comparison(side, "theta"),
    BayesTools:::.hypothesis_side_expression(side), "theta"
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
  fit <- attach_test_parameter_map(fit)
  fit <- .bt_attach_parameter_map(fit, monitor_names = "mu")
  posterior <- marginal_posterior(as_mixed_posteriors(fit, "mu"), "mu",
                                  use_formula = FALSE, prior_samples = TRUE,
                                  n_samples = 64)
  posterior <- .bt_meta_set(posterior, "prior_density", mixture)
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
  expect_identical(prior_density_ordinate(convolved, -1)$behavior, "zero")
  expect_identical(.prior_linear_density_height(convolved, -1), 0)
  convolved <- .prior_linear_combination_density(
    list(a = half_normal, b = half_normal, c = half_normal), c(a = 1, b = 1, c = 1)
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

  formula_scale <- list(mu = formula_scale_for_test(~ x, list(x = list(mean = 5, sd = 2))))
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
  fit <- attach_test_parameter_map(fit)
  fit <- .bt_attach_parameter_map(fit, monitor_names = colnames(posterior))
  mixed <- as_mixed_posteriors(fit, columns, conditional = "mu_x", transform_scaled = TRUE)
  expect_equal(
    as.numeric(.prior_linear_density_height(.bt_meta_get(mixed, "prior_densities")$mu_intercept, 1)),
    as.numeric(.prior_linear_density_height(filtered$mu_intercept, 1)),
    tolerance = 1e-8
  )
})

test_that("marginal posteriors of log-intercept unscaled intercepts use the log-scale combination", {

  # log(intercept) scaling: the unscaled intercept is
  # Y = exp(log(b0) - r b1) = b0 exp(-r b1), r = mean(x) / sd(x), with the
  # positive fitted intercept b0 ~ N(0, .5)T(0, Inf) and the standardized
  # slope b1 ~ N(0, 1). Its density is f(y) = int f_b0(y e^(r t)) e^(r t)
  # phi(t) dt; its log is linear in (log(b0), b1), so the prior density is the
  # exp of the log-scale combination, its support (0, Inf), and it is the
  # scale product of b0 and the lognormal exp(-r b1).
  formula <- ~ x
  attr(formula, "log(intercept)") <- TRUE
  scaled <- JAGS_formula(formula, "mu", data = data.frame(x = c(1, 2, 3.5, 4, 6, 8.5)), prior_list = list(
    intercept = prior("normal", list(0, .5), list(0, Inf)),
    x = prior("normal", list(0, 1))
  ), formula_scale = list(x = TRUE))
  ratio <- scaled$formula_scale[["mu_x"]]$mean / scaled$formula_scale[["mu_x"]]$sd
  set.seed(3)
  n <- 50L
  posterior <- cbind(mu_intercept = abs(stats::rnorm(n, .3, .1)), mu_x = stats::rnorm(n, .1, .2))
  scaled_fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
         monitor = colnames(posterior), sample = n),
    class = c("runjags", "BayesTools_fit")
  )
  attr(scaled_fit, "prior_list") <- scaled$prior_list
  attr(scaled_fit, "formula_design") <- list(mu = scaled$formula_design)
  attr(scaled_fit, "formula_scale") <- list(mu = scaled$formula_scale)
  scaled_fit <- attach_test_parameter_map(scaled_fit)
  scaled_fit <- .bt_attach_fit_contract(.bt_attach_draw_geometry(.bt_attach_parameter_map(scaled_fit)))
  scaled_mixed <- as_mixed_posteriors(scaled_fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)

  reference <- function(y){
    stats::integrate(function(t){
      exp(log(2) + stats::dnorm(y * exp(ratio * t), 0, .5, log = TRUE) + ratio * t +
            stats::dnorm(t, log = TRUE))
    }, -Inf, Inf, rel.tol = 1e-12)$value
  }
  for(prior_samples in c(FALSE, TRUE)){
    intercept <- marginal_posterior(scaled_mixed, "mu_intercept", use_formula = FALSE,
                                    prior_samples = prior_samples)
    expect_equal(as.numeric(intercept),
                 posterior[, "mu_intercept"] * exp(-ratio * posterior[, "mu_x"]))
    support <- .bt_meta_get(intercept, "support")
    expect_identical(support$bounds, c(0, Inf))
    expect_true(support$exact)
  }
  density <- .bt_meta_get(intercept, "prior_density")
  expect_identical(
    .bt_meta_get(intercept, "linear_weights")[c("mu_intercept", "mu_x")],
    c(mu_intercept = 1, mu_x = 0)
  )
  # the scale product has exact ordinates and heights (point hypotheses are
  # eligible) and quadrature region probabilities
  for(value in c(.05, .5, 2)){
    ordinate <- prior_density_ordinate(density, value)
    expect_true(ordinate$exact)
    expect_identical(ordinate$method, "scale_mixture")
    expect_equal(exp(ordinate$log_density), reference(value), tolerance = 1e-10)
    expect_equal(as.numeric(.prior_linear_density_height(density, value)),
                 reference(value), tolerance = 1e-10)
  }
  expect_true(prior_ordinate_status(density, .5)$eligible)
  probability <- .hypothesis_prior_density_prob(
    density, hypothesis_parse("theta > 0.5")$statements[[1L]]$left, "theta"
  )
  expect_equal(as.numeric(probability), stats::integrate(function(t){
    2 * stats::pnorm(.5 * exp(ratio * t), 0, .5, lower.tail = FALSE) * stats::dnorm(t)
  }, -Inf, Inf, rel.tol = 1e-12)$value, tolerance = 1e-10)

  # the marginal of the linear predictor at x = 0 (log(Y)) on the original
  # scale: the log image of the scale product, with density f(e^z) e^z
  predictor <- marginal_posterior(scaled_mixed, "mu_intercept", formula = formula,
                                  prior_samples = TRUE)
  expect_equal(as.numeric(predictor[["intercept"]]), log(as.numeric(intercept)))
  predictor_density <- .bt_meta_get(predictor[["intercept"]], "prior_density")
  expect_true(prior_ordinate_status(predictor_density, -1)$eligible)
  for(value in c(-3, -1, .5)){
    ordinate <- prior_density_ordinate(predictor_density, value)
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), reference(exp(value)) * exp(value),
                 tolerance = 1e-10)
    expect_equal(as.numeric(.prior_linear_density_height(predictor_density, value)),
                 reference(exp(value)) * exp(value), tolerance = 1e-10)
  }
  predictor_probability <- .hypothesis_prior_density_prob(
    predictor_density, hypothesis_parse("theta > -0.7")$statements[[1L]]$left, "theta"
  )
  expect_equal(as.numeric(predictor_probability), stats::integrate(function(t){
    2 * stats::pnorm(exp(-.7) * exp(ratio * t), 0, .5, lower.tail = FALSE) * stats::dnorm(t)
  }, -Inf, Inf, rel.tol = 1e-12)$value, tolerance = 1e-10)
  # its support (the real line) and atoms (none) are declared from structure,
  # so Savage-Dickey runs: the exact prior ordinate (reference integral) over
  # the Gaussian kernel sum of the draws at the null (no support bound)
  level <- predictor[["intercept"]]
  expect_identical(.bt_meta_get(level, "support")$bounds, c(-Inf, Inf))
  expect_true(posterior_atoms_free(level))
  bf <- hypothesis_BF(level, hypothesis = "x = -1", parameter = "x")
  draws <- as.numeric(level)
  expect_equal(attr(bf, "raw_BF")[[1L]],
               reference(exp(-1)) * exp(-1) /
                 mean(stats::dnorm(-1, draws, stats::bw.nrd0(draws))),
               tolerance = 1e-9)
  # its exp (the predictor on the original scale) is the scale product again
  exp_predictor <- marginal_posterior(scaled_mixed, "mu_intercept", formula = formula,
                                      prior_samples = TRUE, transformation = "exp")
  exp_density <- .bt_meta_get(exp_predictor[["intercept"]], "prior_density")
  for(value in c(.05, .5, 2)){
    ordinate <- prior_density_ordinate(exp_density, value)
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), reference(value), tolerance = 1e-10)
  }
})

test_that("formula marginals of log-intercept combinations derive supports, components and atoms from structure", {

  # log(intercept) scaling with an intercept mixture of the point 1.5 and
  # N(1, .5)T(.5, 2), and a spike-and-slab standardized slope: on the
  # original (transform_scaled) scale the linear predictor at x is
  # log(b0) + k(x) b1, linear in (log(b0), b1). Independent derivation from
  # the priors: the component (point, spike) is the point log(1.5), the
  # component (truncated, spike) the interval [log(.5), log(2)], and the
  # components with the slab (k(x) != 0 at the three levels) the real line;
  # the atom log(1.5) has the draw frequency of (point, spike) as its mass.
  formula <- ~ x
  attr(formula, "log(intercept)") <- TRUE
  truncated <- prior("normal", list(1, .5), list(.5, 2))
  scaled <- JAGS_formula(formula, "mu", data = data.frame(x = c(1, 2, 3.5, 4, 6, 8.5)), prior_list = list(
    intercept = prior_mixture(list(prior("point", list(1.5)), truncated), is_null = c(FALSE, FALSE)),
    x = prior_spike_and_slab(prior("normal", list(0, 1)))
  ), formula_scale = list(x = TRUE))
  set.seed(5)
  n <- 200L
  intercept_indicator <- sample(1:2, n, TRUE)
  slope_inclusion <- stats::rbinom(n, 1, .5)
  posterior <- cbind(
    mu_intercept = ifelse(intercept_indicator == 1L, 1.5, rng(truncated, n)),
    mu_x = slope_inclusion * stats::rnorm(n),
    mu_intercept_indicator = intercept_indicator,
    mu_x_indicator = slope_inclusion
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
         monitor = colnames(posterior), sample = n),
    class = c("runjags", "BayesTools_fit")
  )
  attr(fit, "prior_list") <- scaled$prior_list
  attr(fit, "formula_design") <- list(mu = scaled$formula_design)
  attr(fit, "formula_scale") <- list(mu = scaled$formula_scale)
  fit <- attach_test_parameter_map(fit)
  fit <- .bt_attach_fit_contract(.bt_attach_draw_geometry(.bt_attach_parameter_map(fit)))
  mixed <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  levels <- marginal_posterior(mixed, "mu_x", formula = formula, prior_samples = TRUE)

  # the slope's spike-and-slab components are c("alternative", "null")
  slab_first <- attr(scaled$prior_list$mu_x, "components")
  expect_identical(slab_first, c("alternative", "null"))
  expected_keys <- cbind(intercept_indicator, ifelse(slope_inclusion == 1L, 1L, 2L))
  expected_support <- function(key){
    if(key[[2L]] == 1L){
      return(list(bounds = c(-Inf, Inf), points = numeric(), type = "interval"))
    }
    if(key[[1L]] == 1L){
      return(list(bounds = rep(log(1.5), 2L), points = log(1.5), type = "points"))
    }
    list(bounds = log(c(.5, 2)), points = numeric(), type = "interval")
  }
  atom_mass <- mean(intercept_indicator == 1L & slope_inclusion == 0L)
  for(level in names(levels)){
    marginal <- levels[[level]]
    support <- .bt_meta_get(marginal, "support")
    expect_identical(support$bounds, c(-Inf, Inf))
    expect_equal(support$points, log(1.5))
    expect_identical(support$type, "mixed")
    expect_true(support$exact)

    atoms <- .bt_meta_get(marginal, "atoms")
    expect_true(atoms$declared)
    expect_equal(as.numeric(atoms$locations), log(1.5))
    expect_equal(atoms$mass, atom_mass)

    components <- .bt_meta_get(marginal, "components")
    keys <- unname(components$keys[components$index, , drop = FALSE])
    expect_equal(keys, unname(expected_keys), ignore_attr = TRUE)
    for(key_i in seq_len(nrow(components$keys))){
      expected <- expected_support(components$keys[key_i, ])
      observed <- components$supports[[key_i]]
      expect_equal(observed$bounds, expected$bounds)
      expect_equal(observed$points, expected$points)
      expect_identical(observed$type, expected$type)
      expect_true(observed$exact)
    }
  }
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
  formula_scale <- list(mu = formula_scale_for_test(~ x, list(x = list(mean = 5, sd = 2))))
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
  formula_scale <- list(mu = formula_scale_for_test(~ x, list(x = list(mean = 3, sd = 2))))
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

test_that("plotted prior densities evaluate the exact route at every plotted value", {

  # References: closed forms (normal sum, skew-normal sum of a normal and a
  # half-normal, lognormal image) and integrate() at rel.tol 1e-12. The
  # interpolated grid was 8% low at 5 in the normal-sum tail (its sources are
  # truncated at 1e-4 per tail), 6.8% low at the peak of a conditional-normal
  # curve, and 43% off next to a row-mixture jump.
  plotted <- function(density, x_range, n_points = 101){
    .prior_linear_density_to_plot_data(density, n_points = n_points, x_range = x_range)$density
  }
  normal_sum <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1))), c(a = 1, b = 1)
  )
  curve <- plotted(normal_sum, c(-5, 5))
  expect_equal(curve$y, stats::dnorm(curve$x, 0, sqrt(2)), tolerance = 1e-14)

  half_normal <- prior("normal", list(0, 1), list(0, Inf))
  skewed <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, .1)), b = half_normal), c(a = 1, b = 1)
  )
  curve <- plotted(skewed, c(-.5, 3))
  expect_equal(curve$y, 2 * stats::dnorm(curve$x, 0, sqrt(1.01)) * stats::pnorm(10 * curve$x / sqrt(1.01)),
               tolerance = 1e-10)

  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "s"
  mixture <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = slope, s = prior("lognormal", list(0, 1))), c(a = 1, b = 1)
  )
  curve <- plotted(mixture, c(-4, 4), n_points = 21)
  reference <- vapply(curve$x, function(value){
    sum(vapply(list(c(0, 1), c(1, Inf)), function(piece){
      stats::integrate(function(s) stats::dnorm(value, 0, sqrt(1 + s^2)) * stats::dlnorm(s),
                       piece[1L], piece[2L], rel.tol = 1e-12)$value
    }, numeric(1)))
  }, numeric(1))
  expect_equal(curve$y, reference, tolerance = 1e-8)

  context <- .prior_density_context(list(a = half_normal, b = prior("normal", list(0, 1))), c("a", "b"))
  rows <- .prior_density_from_context_rows(context, rbind(c(a = 1, b = 0), c(a = 1, b = 1)))
  curve <- plotted(rows, c(-2, 3))
  expect_equal(curve$y, .5 * ifelse(curve$x >= 0, 2 * stats::dnorm(curve$x), 0) +
                 .5 * stats::dnorm(curve$x, 0, sqrt(2)) * 2 * stats::pnorm(curve$x / sqrt(2)),
               tolerance = 1e-12)

  transformed <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = prior("normal", list(0, 1))), c(a = 1, b = 1),
    output_transformation = "exp"
  )
  curve <- plotted(transformed, c(.01, 5))
  expect_equal(curve$y, stats::dlnorm(curve$x, 0, sqrt(2)), tolerance = 1e-12)
})

test_that("plotted quadrature routes share one batched quadrature on a bounded display grid", {

  # A route with quadrature leaves is plotted on at most 200 equally spaced
  # values plus the values it must include, and all regular values share one
  # batched quadrature: no per-value QUADPACK piece is integrated (the
  # per-value quadrature at each of 1000 values took 1-3.5 s per overlay).
  # References: integrate() at rel.tol 1e-12 and closed forms; the batched
  # quadrature's relative tolerance is 1e-8.
  pieces <- 0L
  piece <- .prior_conditional_normal_piece
  local_mocked_bindings(.prior_conditional_normal_piece = function(...){
    pieces <<- pieces + 1L
    piece(...)
  })

  # the transform-scaled intercept of a t-slab coefficient: N(0, 1) + 2.5 T,
  # T ~ t(0, .5, 3), a Gaussian convolution with its peak at the offset 0
  convolution <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = prior("t", list(0, .5, 3))), c(a = 1, b = 2.5)
  )
  curve <- .prior_linear_density_to_plot_data(convolution, x_range = c(-8, 8))$density
  expect_lte(length(curve$x), 200L)
  expect_true(0 %in% curve$x)
  expect_identical(curve$x[which.max(curve$y)], 0)
  reference <- vapply(curve$x, function(value){
    stats::integrate(function(u) stats::dnorm(value - 2.5 * u) * stats::dt(u / .5, 3) / .5,
                     -Inf, Inf, rel.tol = 1e-12)$value
  }, numeric(1))
  expect_equal(curve$y, reference, tolerance = 1e-8)

  # a mixture with an atom at 1, a jump at 1.5 (a t(0, .5, 3) term truncated
  # to [.5, Inf) next to the spike) and a Gaussian convolution
  a <- prior_mixture(list(prior("spike", list(1), prior_weights = 1),
                          prior("normal", list(1, .3), prior_weights = 1)), is_null = c(TRUE, FALSE))
  b <- prior_mixture(list(prior("spike", list(0), prior_weights = 1),
                          prior("t", list(0, .5, 3), list(.5, Inf), prior_weights = 1)),
                     is_null = c(TRUE, FALSE))
  mixture <- .prior_linear_combination_density(list(a = a, b = b), c(a = 1, b = 1))
  plot_data <- .prior_linear_density_to_plot_data(mixture, x_range = c(-1, 4))
  curve <- plot_data$density
  expect_lte(length(curve$x), 205L)
  jump_delta <- 1e-6 * 5
  expect_true(all(c(1, 1.5, 1.5 - jump_delta) %in% curve$x))
  # 1.5 is a lower bound whose plotted value is the limit inside the support,
  # so a value just inside it would repeat that value
  expect_false((1.5 + jump_delta) %in% curve$x)
  expect_equal(plot_data$points1$x, 1)
  expect_equal(plot_data$points1$y, .25)
  truncated_t <- function(u){
    ifelse(u >= .5, stats::dt(u / .5, 3) / .5 / stats::pt(1, 3, lower.tail = FALSE), 0)
  }
  reference <- vapply(curve$x, function(value){
    convolution <- stats::integrate(function(u) stats::dnorm(value - u, 1, .3) * truncated_t(u),
                                    .5, Inf, rel.tol = 1e-12)$value
    .25 * truncated_t(value - 1) + .25 * stats::dnorm(value, 1, .3) + .25 * convolution
  }, numeric(1))
  expect_equal(curve$y, reference, tolerance = 1e-8)

  # a peak off every structural point: N(0, 1) + G, G ~ gamma(3, 1); the vertex
  # of the parabola through the largest plotted value and its neighbours is
  # plotted, so the plotted maximum is the density's maximum (optimize() at
  # tol 1e-10 over the reference) to 1e-7
  skewed <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = prior("gamma", list(3, 1))), c(a = 1, b = 1)
  )
  curve <- .prior_linear_density_to_plot_data(skewed, x_range = c(-5, 15))$density
  skewed_reference <- function(value){
    stats::integrate(function(u) stats::dnorm(value - u) * stats::dgamma(u, 3, 1),
                     0, Inf, rel.tol = 1e-12)$value
  }
  maximum <- stats::optimize(skewed_reference, c(0, 5), maximum = TRUE, tol = 1e-10)
  expect_equal(max(curve$y), maximum$objective, tolerance = 1e-7)
  expect_identical(pieces, 0L)
})

test_that("plotted densities draw a support bound with one value outside it", {

  # The plotted value at a finite support bound is the density's limit inside
  # the support, so the display grid adds the value 1e-6 of the plotted range
  # outside the bound (below a lower bound, above an upper bound) and not the
  # one inside it, which repeated the value at the bound (the first two points
  # of the exp(affine) scale-intercept prior coincided in the plot). The
  # inside value stays where the density at the bound is unavailable,
  # and both sides where one term's support ends and another's
  # starts. References: the scale-product offset f(0+) E[1 / W] in closed
  # form, and truncated-t densities with integrate() at rel.tol 1e-12.
  z0 <- -0.98946052195266332
  intercept <- .prior_linear_combination_density(
    list(b0 = prior("normal", list(0, 1 / sqrt(8)), list(0, Inf)), b1 = prior("normal", list(0, .5))),
    c(b0 = 1, b1 = z0), source_transforms = c(b0 = "log", b1 = NA), output_transformation = "exp"
  )
  curve <- .prior_linear_density_to_plot_data(intercept, x_range = c(0, 2))$density
  # the bound is repeated once, with density 0 (the edge of density.prior())
  expect_identical(curve$x[1:2], c(0, 0))
  expect_identical(curve$y[1L], 0)
  expect_equal(curve$y[2L], 2 * stats::dnorm(0, 0, 1 / sqrt(8)) * exp((.5 * z0)^2 / 2), tolerance = 1e-12)
  expect_gt(min(diff(curve$x[-1L])), 1e-3)
  expect_lte(length(curve$x), 201L)

  # an upper bound: the negated product of a gamma(2, 1) term and a
  # lognormal scale has support (-Inf, 0]
  term <- prior("gamma", list(2, 1))
  attr(term, "multiply_by") <- "s"
  negated <- .prior_linear_combination_density(
    list(b = term, s = prior("lognormal", list(0, .5))), c(b = 1),
    output_transformation = "lin", output_transformation_arguments = list(a = 0, b = -1)
  )
  curve <- .prior_linear_density_to_plot_data(negated, x_range = c(-5, 1))$density
  delta <- 1e-6 * 6
  expect_true(all(c(0, delta) %in% curve$x))
  expect_false(-delta %in% curve$x)
  expect_equal(curve$y[curve$x == delta], 0)

  # A singular bound stays off the curve without an artificial edge to zero
  # or a nearby point whose arbitrary height would dominate the vertical axis.
  term <- prior("gamma", list(.5, 1))
  attr(term, "multiply_by") <- "s"
  singular <- .prior_linear_combination_density(
    list(b = term, s = prior("lognormal", list(0, .5))), c(b = 1)
  )
  curve <- .prior_linear_density_to_plot_data(singular, x_range = c(0, 5))$density
  expect_identical(prior_density_ordinate(singular, 0)$behavior, "infinite")
  expect_gt(curve$x[1L], 1e-3)
  expect_false((1e-6 * 5) %in% curve$x)
  expect_true(all(is.finite(curve$y)))
  expect_equal(curve$y[1L], exp(prior_density_ordinate(singular, curve$x[1L])$log_density), tolerance = 1e-8)

  # one term's support ends at 1 and another's starts there: both sides
  a <- prior_mixture(list(prior("spike", list(0)), prior("normal", list(0, .3))), is_null = c(TRUE, FALSE))
  b <- prior_mixture(list(prior("t", list(0, .5, 3), list(1, Inf)), prior("t", list(0, .5, 3), list(-Inf, 1))),
                     is_null = c(FALSE, FALSE))
  split <- .prior_linear_combination_density(list(a = a, b = b), c(a = 1, b = 1))
  curve <- .prior_linear_density_to_plot_data(split, x_range = c(-2, 4))$density
  delta <- 1e-6 * 6
  at <- c(1 - delta, 1, 1 + delta)
  expect_true(all(at %in% curve$x))
  upper_part <- function(u) ifelse(u >= 1, stats::dt(u / .5, 3) / .5 / stats::pt(2, 3, lower.tail = FALSE), 0)
  lower_part <- function(u) ifelse(u <= 1, stats::dt(u / .5, 3) / .5 / stats::pt(2, 3), 0)
  convolution <- function(value, part, lower, upper){
    stats::integrate(function(u) stats::dnorm(value - u, 0, .3) * part(u), lower, upper, rel.tol = 1e-12)$value
  }
  reference <- vapply(at, function(value){
    .25 * (upper_part(value) + lower_part(value)) +
      .25 * convolution(value, upper_part, 1, Inf) + .25 * convolution(value, lower_part, -Inf, 1)
  }, numeric(1))
  expect_equal(curve$y[match(at, curve$x)], reference, tolerance = 1e-8)
})

test_that("allocated variance prior plots retain regular densities without singular-boundary spikes", {

  # Kearon's group variance has the form Z = 2 U T^2, where U is uniform
  # and T is half-normal. Integrating f_T(t) / (2 t^2) from sqrt(z / 2)
  # to infinity gives this closed form, independently of the route quadrature.
  scale <- .86
  reference <- function(z){
    lower <- sqrt(z / 2)
    stats::dnorm(lower, sd = scale) / lower -
      stats::pnorm(lower / scale, lower.tail = FALSE) / scale^2
  }
  for(mass in c(1, .75)){
    density <- .prior_allocation_product_density(
      prior("normal", list(0, scale), list(0, Inf)),
      multipliers = list(.prior_allocation_point(0, 1 - mass),
                         .prior_allocation_share(1, 1, 2, mass)),
      n_grid = 256, square = TRUE
    )
    before <- density
    plot_data <- .prior_linear_density_to_plot_data(density, n_points = 1000, x_range = c(0, 2))
    curve <- plot_data$density
    expect_gt(min(curve$x), 1e-3)
    expect_false(2e-6 %in% curve$x)
    expect_true(all(is.finite(curve$y)))
    expect_equal(curve$y, mass * reference(curve$x), tolerance = 1e-8)
    expect_identical(density, before)
    expect_identical(prior_density_ordinate(density, 0)$behavior,
                     if(mass == 1) "infinite" else "point_mass")
    if(mass < 1){
      expect_identical(plot_data$points1$x, 0)
      expect_identical(plot_data$points1$y, 1 - mass)
    }else{
      expect_null(plot_data$points1)
    }
  }

  # A decreasing transformation maps the singular lower bound to an upper
  # bound. Its regular curve still has the same densities and declared atom.
  reflected <- .prior_linear_density_to_plot_data(
    density, n_points = 1000, x_range = c(0, 2),
    transformation = "lin", transformation_arguments = list(a = 0, b = -1)
  )
  expect_true(all(reflected$density$x < 0))
  expect_equal(reflected$density$y, .75 * reference(-reflected$density$x), tolerance = 1e-8)
  expect_identical(reflected$points1$x, 0)
  expect_identical(reflected$points1$y, .25)
})

# Plotted curves of the prior densities of the support-edge tests.
.edge_curve <- function(density, ...){
  .prior_linear_density_to_plot_data(density, n_points = 101, ...)$density
}

# The first two points (lower) and the last two points (upper) of a curve: a
# bound repeated with density 0 on the outer side, then its density.
expect_support_edge <- function(curve, lower = NULL, upper = NULL){
  n <- length(curve$x)
  if(!is.null(lower)){
    expect_identical(curve$x[1:2], c(lower$x, lower$x))
    expect_identical(curve$y[1L], 0)
    expect_equal(curve$y[2L], lower$y, tolerance = 1e-12)
  }
  if(!is.null(upper)){
    expect_equal(curve$x[(n - 1L):n], c(upper$x, upper$x), tolerance = 1e-14)
    expect_equal(curve$y[n - 1L], upper$y, tolerance = 1e-12)
    expect_identical(curve$y[n], 0)
  }
}

test_that("plotted prior curves drop to zero at the lower bound of a half-normal prior", {

  # density.prior() repeats a truncation bound with density 0 (truncate_end):
  # a density that jumps from 0 to a positive value at a finite support bound
  # is drawn with a vertical edge, (lower, 0) then (lower, f), and (upper, f)
  # then (upper, 0) at an upper bound. Reference: the closed-form density (the
  # closed-form route is exact, tolerance 1e-12).
  scale <- .86
  half_normal <- .prior_linear_combination_density(
    list(x = prior("normal", list(0, scale), list(0, Inf))), c(x = 1), n_grid = 512
  )
  curve <- .edge_curve(half_normal)
  expect_support_edge(curve, lower = list(x = 0, y = 2 * stats::dnorm(0, 0, scale)))
  expect_length(curve$x, 102L)
  expect_false(is.unsorted(curve$x))
  expect_equal(curve$y[-1L], 2 * stats::dnorm(curve$x[-1L], 0, scale), tolerance = 1e-12)

  # the same prior reached through the posterior-sample metadata, the path of
  # the prior curves of the RoBMA plots
  theta <- .bt_meta_update(
    structure(abs(stats::qnorm(seq(.01, .99, length.out = 64))),
              class = c("mixed_posteriors.simple", "mixed_posteriors"),
              prior_list = list(prior("normal", list(0, scale), list(0, Inf)))),
    prior_density = half_normal
  )
  attached <- .plot_data_attached_prior_density(list(theta = theta), "theta", n_points = 101)$density
  expect_identical(attached, curve)
})

test_that("plotted prior curves drop to zero at both bounds of a uniform prior and one of its root", {

  # Uniform(0, 2) directly and as 2 B for B ~ Beta(1, 1) under a linear map (the
  # shape of the variance allocation multiplier of a heterogeneity prior): the
  # density .5 jumps at 0 and 2. Its square root sqrt(2 B), through exp_lin, has
  # the density t on (0, sqrt(2)): 0 at 0 (no zero value is added there) and
  # sqrt(2) at the upper bound. References: closed forms, tolerance 1e-12.
  uniform <- .prior_linear_combination_density(
    list(x = prior("uniform", list(0, 2))), c(x = 1), n_grid = 512
  )
  both <- list(lower = list(x = 0, y = .5), upper = list(x = 2, y = .5))
  expect_support_edge(.edge_curve(uniform), lower = both$lower, upper = both$upper)
  beta_unit <- prior("beta", list(1, 1))
  doubled <- .prior_linear_combination_density(
    list(source = beta_unit), c(source = 1), n_grid = 512,
    output_transformation = "lin", output_transformation_arguments = list(a = 0, b = 2)
  )
  curve <- .edge_curve(doubled)
  expect_support_edge(curve, lower = both$lower, upper = both$upper)
  expect_equal(curve$y[2:(length(curve$y) - 1L)], rep(.5, length(curve$y) - 2L), tolerance = 1e-12)

  root <- .prior_linear_combination_density(
    list(source = beta_unit), c(source = 1), n_grid = 512,
    output_transformation = "exp_lin", output_transformation_arguments = list(a = log(2) / 2, b = .5)
  )
  curve <- .edge_curve(root)
  n <- length(curve$x)
  expect_support_edge(curve, upper = list(x = sqrt(2), y = sqrt(2)))
  expect_identical(curve$x[1L], 0)
  expect_identical(curve$y[1L], 0)
  expect_gt(curve$x[2L], 0)
  expect_equal(curve$y[-n], curve$x[-n], tolerance = 1e-12)
})

test_that("plotted prior curve edges follow the plotted range and the transformation", {

  # A bound outside the plotted range has no zero value; a bound strictly
  # inside it has the whole edge, below which the density is 0 (references:
  # closed forms, tolerance 1e-12).
  scale <- .86
  half_normal <- .prior_linear_combination_density(
    list(x = prior("normal", list(0, scale), list(0, Inf))), c(x = 1), n_grid = 512
  )
  uniform <- .prior_linear_combination_density(
    list(x = prior("uniform", list(0, 2))), c(x = 1), n_grid = 512
  )
  curve <- .edge_curve(half_normal, x_range = c(.5, 3))
  expect_identical(curve$x[1L], .5)
  expect_equal(curve$y[1L], 2 * stats::dnorm(.5, 0, scale), tolerance = 1e-12)
  expect_false(anyDuplicated(curve$x) > 0L)
  curve <- .edge_curve(uniform, x_range = c(.5, 1.5))
  expect_false(anyDuplicated(curve$x) > 0L)
  expect_true(all(curve$y == .5))
  curve <- .edge_curve(half_normal, x_range = c(-1, 3))
  at_zero <- which(curve$x == 0)
  expect_length(at_zero, 2L)
  expect_identical(at_zero[2L], at_zero[1L] + 1L)
  expect_identical(curve$y[at_zero[1L]], 0)
  expect_equal(curve$y[at_zero[2L]], 2 * stats::dnorm(0, 0, scale), tolerance = 1e-12)
  expect_true(all(curve$y[curve$x < 0] == 0))
  expect_equal(curve$y[curve$x > 0], 2 * stats::dnorm(curve$x[curve$x > 0], 0, scale), tolerance = 1e-12)
  curve <- .edge_curve(uniform, x_range = c(-1, 3))
  expect_identical(sum(curve$x == 0), 2L)
  expect_identical(sum(curve$x == 2), 2L)
  inner <- curve$x > 0 & curve$x < 2
  expect_equal(curve$y[inner], rep(.5, sum(inner)), tolerance = 1e-12)
  expect_true(all(curve$y[curve$x < 0 | curve$x > 2] == 0))

  # a transformation maps the edge by the rules of density.prior(): the
  # logarithm sends the bound 0 to -Inf, which is dropped with its zero value,
  # and the bound 2 to log(2), where the density of log(T) is 2 f(2) = 1
  log_map <- list(fun = log, inv = exp, jac = function(x) 1 / x)
  curve <- .edge_curve(half_normal, transformation = log_map)
  expect_true(all(is.finite(curve$x)))
  expect_false(anyDuplicated(curve$x) > 0L)
  curve <- .edge_curve(uniform, transformation = log_map)
  expect_true(all(is.finite(curve$x)))
  n <- length(curve$x)
  expect_equal(curve$x[(n - 1L):n], rep(log(2), 2L), tolerance = 1e-14)
  expect_equal(curve$y[n - 1L], 1, tolerance = 1e-12)
  expect_identical(curve$y[n], 0)
  expect_false(anyDuplicated(curve$x[-((n - 1L):n)]) > 0L)

  # a plotted range on the transformed scale: 1 + 3 T has the support [1, 7],
  # whose lower bound is inside the range [.1, 4] and upper bound outside it
  curve <- .edge_curve(uniform, x_range = c(.1, 4), transformation = "lin",
                       transformation_arguments = list(a = 1, b = 3), transformation_settings = TRUE)
  at_edge <- which(curve$x == 1)
  expect_length(at_edge, 2L)
  expect_identical(at_edge[2L], at_edge[1L] + 1L)
  expect_identical(curve$y[at_edge[1L]], 0)
  expect_equal(curve$y[at_edge[2L]], .5 / 3, tolerance = 1e-12)
  expect_true(all(curve$y[curve$x < 1] == 0))
  expect_equal(curve$y[curve$x > 1], rep(.5 / 3, sum(curve$x > 1)), tolerance = 1e-12)
  expect_identical(attr(curve, "x_range"), c(.1, 4))
  expect_false(anyDuplicated(curve$x[curve$x > 1]) > 0L)

  # a decreasing map keeps the zero value on the outer side of each bound: 1 - 3 T
  # has the support [-5, 1]; below the lower bound of T (plotted at 1) the
  # plotted values continue to the right, above its upper bound (plotted at -5)
  # to the left
  curve <- .edge_curve(uniform, x_range = c(-6, 2), transformation = "lin",
                       transformation_arguments = list(a = 1, b = -3), transformation_settings = TRUE)
  at_edge <- which(curve$x == 1)
  expect_length(at_edge, 2L)
  expect_equal(curve$y[at_edge], c(.5 / 3, 0), tolerance = 1e-12)
  expect_identical(at_edge[2L], at_edge[1L] + 1L)
  at_edge <- which(curve$x == -5)
  expect_length(at_edge, 2L)
  expect_equal(curve$y[at_edge], c(0, .5 / 3), tolerance = 1e-12)
  expect_identical(at_edge[2L], at_edge[1L] + 1L)
  inner <- curve$x > -5 & curve$x < 1
  expect_equal(curve$y[inner], rep(.5 / 3, sum(inner)), tolerance = 1e-12)
  curve <- .edge_curve(uniform, transformation = "lin", transformation_arguments = list(a = 0, b = -1))
  n <- length(curve$x)
  expect_equal(curve$x[1:2], c(0, 0))
  expect_equal(curve$y[1:2], c(0, .5), tolerance = 1e-12)
  expect_equal(curve$x[(n - 1L):n], c(-2, -2), tolerance = 1e-14)
  expect_equal(curve$y[(n - 1L):n], c(.5, 0), tolerance = 1e-12)

  # exp_lin maps the source bound 0 to 0 with the analytic limits of
  # density.prior(): exp(1) T keeps the edge (0, 0), (0, f(0) / exp(1)), and
  # sqrt(T) has the density 0 at 0. Reference: density.prior() of the same
  # prior on the same grid (both evaluate the closed form, tolerance 1e-12).
  for(arguments in list(list(a = 1, b = 1), list(a = 0, b = .5))){
    curve <- .edge_curve(half_normal, x_range = c(0, 3), transformation = "exp_lin",
                         transformation_arguments = arguments)
    reference <- density(prior("normal", list(0, scale), list(0, Inf)), x_range = c(0, 3),
                         n_points = 101, transformation = "exp_lin",
                         transformation_arguments = arguments)
    expect_identical(curve$x[1:2], c(0, 0))
    expect_identical(curve$y[1L], 0)
    expect_equal(curve$x, reference$x, tolerance = 1e-12)
    expect_equal(curve$y, reference$y, tolerance = 1e-12)
  }
  expect_equal(curve$y[2L], 0)
})

test_that("plotted curves keep the edge at an output-transformed bound that rounds outside its source support", {

  # The inverse image of a mapped bound can round one step outside the source
  # support: (0.1 + 0.3 - 0.1) / 0.3 is 1.0000000000000002, so the route density
  # of 0.1 + 0.3 U is 0 at its upper bound 0.4, as it is for sqrt(5 B) at
  # sqrt(5). The plotted value at the bound is the one-sided limit inside the
  # support that prior_density_ordinate() returns (it matches the endpoint), and
  # the edge (bound, f), (bound, 0) follows. References: the closed forms (1 / b
  # for a + b U, 2 y / c for sqrt(c B)) and the ordinate, tolerance 1e-12.
  u01 <- prior("uniform", list(0, 1))
  b11 <- prior("beta", list(1, 1))
  ordinate <- function(density, value){
    exp(prior_density_ordinate(density, value)$log_density)
  }
  shifted <- function(a, b){
    .prior_linear_combination_density(
      list(s = u01), c(s = 1), n_grid = 64,
      output_transformation = "lin", output_transformation_arguments = list(a = a, b = b)
    )
  }
  root <- function(c){
    .prior_linear_combination_density(
      list(s = b11), c(s = 1), n_grid = 64,
      output_transformation = "exp_lin", output_transformation_arguments = list(a = log(c) / 2, b = .5)
    )
  }
  # both edges of a + b U: bounds a and a + b with the value 1 / |b|, equal to the
  # ordinate there
  expect_shifted_edges <- function(a, b, ...){
    density <- shifted(a, b)
    curve <- .edge_curve(density, ...)
    n <- length(curve$x)
    bounds <- sort(c(a, a + b))
    expect_equal(curve$x[1:2], rep(bounds[1L], 2L), tolerance = 1e-14)
    expect_equal(curve$x[(n - 1L):n], rep(bounds[2L], 2L), tolerance = 1e-14)
    expect_identical(curve$y[1L], 0)
    expect_identical(curve$y[n], 0)
    expect_equal(curve$y[2L], 1 / abs(b), tolerance = 1e-12)
    expect_equal(curve$y[n - 1L], 1 / abs(b), tolerance = 1e-12)
    expect_equal(curve$y[2L], ordinate(density, curve$x[2L]), tolerance = 1e-12)
    expect_equal(curve$y[n - 1L], ordinate(density, curve$x[n]), tolerance = 1e-12)
    expect_false(is.unsorted(curve$x))
  }

  # 0.1 + 0.3 U
  expect_gt((.1 + .3 - .1) / .3, 1)
  expect_shifted_edges(.1, .3)
  curve <- .edge_curve(shifted(.1, .3))
  expect_equal(curve$y[-c(1L, length(curve$y))], rep(1 / .3, length(curve$y) - 2L), tolerance = 1e-12)

  # sqrt(5 B): the density 2 y / 5 on (0, sqrt(5)) is 0 at the lower bound (no
  # edge there) and 2 / sqrt(5) at the upper bound
  curve <- .edge_curve(root(5))
  n <- length(curve$x)
  expect_equal(curve$x[(n - 1L):n], rep(sqrt(5), 2L), tolerance = 1e-14)
  expect_equal(curve$y[n - 1L], 2 / sqrt(5), tolerance = 1e-12)
  expect_equal(curve$y[n - 1L], ordinate(root(5), sqrt(5)), tolerance = 1e-12)
  expect_identical(curve$y[n], 0)
  expect_identical(curve$x[1L], 0)
  expect_identical(curve$y[1L], 0)
  expect_gt(curve$x[2L], 0)
  expect_equal(curve$y[-n], 2 * curve$x[-n] / 5, tolerance = 1e-12)

  # sqrt(c B) for c = 2, ..., 30 (c = 5 rounded outside, among others)
  for(c in 2:30){
    density <- root(c)
    curve <- .edge_curve(density)
    n <- length(curve$x)
    expect_equal(curve$x[(n - 1L):n], rep(sqrt(c), 2L), tolerance = 1e-14)
    expect_equal(curve$y[n - 1L], 2 / sqrt(c), tolerance = 1e-12)
    expect_equal(curve$y[n - 1L], ordinate(density, curve$x[n]), tolerance = 1e-12)
    expect_identical(curve$y[n], 0)
    expect_identical(curve$y[1L], 0)
    expect_gt(curve$x[2L], 0)
  }

  # the 30 cases a + b U (a in 0, .1, .2, .7, -1.3; b in .1, .3, .7, 1.3, 2.9, 3.7) and a
  # decreasing map a - |b| U: both bounds have the edge
  for(a in c(0, .1, .2, .7, -1.3)){
    for(b in c(.1, .3, .7, 1.3, 2.9, 3.7)){
      expect_shifted_edges(a, b)
    }
  }
  for(a in c(.7, 1.3, 2.9)){
    for(b in c(-.1, -.3, -.7, -1.3, -2.9)){
      expect_shifted_edges(a, b)
    }
  }

  # a bound strictly inside the plotted range has the whole edge with zeros beyond
  # it, and a bound outside the range none
  curve <- .edge_curve(shifted(.1, .3), x_range = c(0, 1))
  at_bound <- which(curve$x == .4)
  expect_length(at_bound, 2L)
  expect_equal(curve$y[at_bound], c(1 / .3, 0), tolerance = 1e-12)
  expect_identical(at_bound[2L], at_bound[1L] + 1L)
  expect_true(all(curve$y[curve$x > .4] == 0))
  curve <- .edge_curve(shifted(.1, .3), x_range = c(.2, .35))
  expect_false(anyDuplicated(curve$x) > 0L)
  expect_equal(curve$y, rep(1 / .3, length(curve$y)), tolerance = 1e-12)

  # a plotted transformation maps the edge: 1 + 2 (0.1 + 0.3 U) has the bound 1.8
  # and the density 1 / .6
  curve <- .edge_curve(shifted(.1, .3), transformation = "lin",
                       transformation_arguments = list(a = 1, b = 2))
  n <- length(curve$x)
  expect_equal(curve$x[(n - 1L):n], rep(1.8, 2L), tolerance = 1e-14)
  expect_equal(curve$y[(n - 1L):n], c(1 / .6, 0), tolerance = 1e-12)

  # only the plotted value at a bound changes: the route density used by grids is
  # the route's own (0 at the rounded bound), and the interior plotted values are
  # the route's densities
  route <- .prior_density_route_from_adaptive(
    attr(shifted(.1, .3), "adaptive_evaluation", exact = TRUE)
  )
  expect_identical(.prior_density_route_density(route, .1 + .3), 0)
  curve <- .edge_curve(shifted(.1, .3))
  inner <- seq(3L, length(curve$x) - 2L)
  expect_identical(curve$y[inner], .prior_density_route_density(route, curve$x[inner]))
})

test_that("plotted curves keep the edge at a bound of a scalar affine prior that rounds outside its source support", {

  # w U + o for U ~ U(0, 1) and a point term o (the scalar route with an
  # offset): the inverse image (o + w - o) / w of the upper bound rounds just
  # above 1 for (w, o) = (.1, .2), (.3, .1), (.3, .7), where prior_density_ordinate()
  # was 0 (outside the support) and so was the edge. The ordinate there is the
  # one-sided limit 1 / |w|, so both bounds of every combination have the edge
  # (bound, 0), (bound, 1 / |w|). References: the closed form 1 / |w| and the
  # ordinate, tolerance 1e-12.
  affine <- function(w, o){
    .prior_linear_combination_density(
      list(s = prior("uniform", list(0, 1)), p = prior("point", list(o))),
      c(s = w, p = 1), n_grid = 64
    )
  }
  expect_affine_edges <- function(w, o, ...){
    density <- affine(w, o)
    label <- sprintf("w = %g, o = %g", w, o)
    expect_identical(
      .prior_density_route_from_adaptive(attr(density, "adaptive_evaluation", exact = TRUE))$type,
      "scalar", info = label
    )
    curve <- .edge_curve(density, ...)
    n <- length(curve$x)
    bounds <- sort(c(o, o + w))
    expect_equal(curve$x[1:2], rep(bounds[1L], 2L), tolerance = 1e-14, info = label)
    expect_equal(curve$x[(n - 1L):n], rep(bounds[2L], 2L), tolerance = 1e-14, info = label)
    expect_identical(curve$y[1L], 0, info = label)
    expect_identical(curve$y[n], 0, info = label)
    expect_equal(curve$y[2L], 1 / abs(w), tolerance = 1e-12, info = label)
    expect_equal(curve$y[n - 1L], 1 / abs(w), tolerance = 1e-12, info = label)
    expect_equal(
      curve$y[n - 1L], exp(prior_density_ordinate(density, curve$x[n])$log_density),
      tolerance = 1e-12, info = label
    )
    expect_false(is.unsorted(curve$x), info = label)
    expect_equal(curve$y[-c(1L, n)], rep(1 / abs(w), n - 2L), tolerance = 1e-12, info = label)
  }

  # the rounded combinations and the 12 combinations of the probe, w in .1, .3,
  # .7, 1.3 and o in .1, .2, .7, and the decreasing maps
  for(combination in list(c(.1, .2), c(.3, .1), c(.3, .7))){
    expect_gt((combination[2L] + combination[1L] - combination[2L]) / combination[1L], 1)
  }
  for(w in c(.1, .3, .7, 1.3, -.1, -.3, -.7, -1.3)){
    for(o in c(.1, .2, .7, -.35)){
      expect_affine_edges(w, o)
    }
  }

  # a bound strictly inside the plotted range keeps the whole edge, with zeros
  # beyond it
  curve <- .edge_curve(affine(.1, .2), x_range = c(0, 1))
  upper <- .2 + .1
  at_bound <- which(abs(curve$x - upper) < 1e-14)
  expect_length(at_bound, 2L)
  expect_equal(curve$y[at_bound], c(10, 0), tolerance = 1e-12)
  expect_identical(at_bound[2L], at_bound[1L] + 1L)
  expect_true(all(curve$y[curve$x > upper + 1e-14] == 0))
})

test_that("plotted quadrature curves drop to zero at a support bound and keep interior jumps", {

  # A quadrature route (an exponential(1) term scaled by a lognormal(0, .5)
  # factor, a scale product) with a positive density at its bound 0: the
  # density at 0 is f_X(0) E[1 / W] = exp(.5^2 / 2), and the bound is repeated
  # with density 0 as for the closed forms (it was drawn without it); the other
  # values are integrate() references at rel.tol 1e-12, split around the peak
  # of the scale density, within the batched quadrature's 1e-8.
  term <- prior("exp", list(1))
  attr(term, "multiply_by") <- "s"
  scaled <- .prior_linear_combination_density(
    list(b = term, s = prior("lognormal", list(0, .5))), c(b = 1)
  )
  curve <- .edge_curve(scaled, x_range = c(0, 6))
  expect_identical(curve$x[1:2], c(0, 0))
  expect_identical(curve$y[1L], 0)
  expect_equal(curve$y[2L], exp(.5^2 / 2), tolerance = 1e-8)
  expect_lte(length(curve$x), 205L)
  inner <- curve$x > 0
  reference <- vapply(curve$x[inner], function(value){
    sum(vapply(list(c(0, .5), c(.5, 1), c(1, 2), c(2, Inf)), function(piece){
      stats::integrate(function(w) exp(-value / w) / w * stats::dlnorm(w, 0, .5),
                       piece[1L], piece[2L], rel.tol = 1e-12)$value
    }, numeric(1)))
  }, numeric(1))
  expect_equal(curve$y[inner], reference, tolerance = 1e-8)

  # an interior jump: the bound 1.5 of the truncated t term is inside the
  # support of the normal term (positive on both sides), so no value with
  # density 0 is added there (the value just outside the bound still is)
  a <- prior_mixture(list(prior("spike", list(1), prior_weights = 1),
                          prior("normal", list(1, .3), prior_weights = 1)), is_null = c(TRUE, FALSE))
  b <- prior_mixture(list(prior("spike", list(0), prior_weights = 1),
                          prior("t", list(0, .5, 3), list(.5, Inf), prior_weights = 1)),
                     is_null = c(TRUE, FALSE))
  jump <- .prior_linear_combination_density(list(a = a, b = b), c(a = 1, b = 1))
  curve <- .edge_curve(jump, x_range = c(-1, 4))
  expect_false(anyDuplicated(curve$x) > 0L)
  expect_true(all(curve$y > 0))
})

test_that("plotted curves interpolated from a numerical grid drop to zero at an exact support bound", {

  # Three gamma(.3, 1) terms have no structural route (three non-normal terms):
  # the curve is interpolated from the numerical grid. The support hull of its
  # provenance is [0, Inf), so the curve leaves the bound 0 by the edge
  # (0, 0) -> (0, the grid's value at 0); there is no edge for a plotted range
  # that starts above the bound. The grid is the reference of its own curve.
  terms <- lapply(1:3, function(i) prior("gamma", list(.3, 1)))
  names(terms) <- c("a", "b", "c")
  gamma_sum <- .prior_linear_combination_density(terms, c(a = 1, b = 1, c = 1), n_grid = 512)
  expect_identical(
    .prior_density_route_from_adaptive(attr(gamma_sum, "adaptive_evaluation", exact = TRUE))$type,
    "unknown"
  )
  expect_identical(gamma_sum$density$x[1L], 0)
  curve <- .edge_curve(gamma_sum)
  expect_identical(curve$x[1:2], c(0, 0))
  expect_identical(curve$y[1L], 0)
  expect_equal(curve$y[2L], gamma_sum$density$y[1L] * gamma_sum$density$mass, tolerance = 1e-12)
  expect_false(is.unsorted(curve$x))
  curve <- .edge_curve(gamma_sum, x_range = c(.5, 3))
  expect_false(anyDuplicated(curve$x) > 0L)
  expect_true(all(curve$y > 0))
})

test_that("batched densities compute the single-value breakpoints and integrals for all values at once", {

  # The batched plans compute the breakpoints of all values together; they
  # are identical to the single-value breakpoints of the ordinates, including
  # peak windows, scale points, the thinning next to singular bounds (Beta
  # shares at 0 and 1) and the merging of close points.
  multipliers <- list(
    prior("normal", list(0, 1), list(0, Inf)), prior("normal", list(1, 2), list(-1, 3)),
    prior("beta", list(.5, 1.5)), prior("beta", list(1.5, .5)), prior("beta", list(.3, .3)),
    prior("gamma", list(.5, 1)), prior("invgamma", list(1, .15)), prior("lognormal", list(0, 1)),
    prior("cauchy", list(0, 1), list(0, 5)), prior("uniform", list(-2, 3))
  )
  values <- c(0, 1e-12, -1e-9, -3.7, -.41, .2, .2, .7 + 1e-15, 2.9, 55, -1e7)
  for(multiplier in multipliers){
    bounds <- unlist(multiplier$truncation[c("lower", "upper")], use.names = FALSE)
    setup <- .prior_conditional_normal_breakpoint_setup(multiplier, bounds)
    for(config in list(c(0, 0, 0, 1), c(.7, 0, 1, .3), c(0, 1, -2.5, .3), c(.7, 1, 0, 4))){
      spec <- list(additive_mean = config[1L], additive_sd = config[2L],
                   product_mean = config[3L], product_sd = config[4L],
                   multiplier = multiplier, bounds = bounds)
      expect_identical(
        .prior_conditional_normal_breakpoints_values(spec, values, setup),
        lapply(values, function(value){
          .prior_conditional_normal_breakpoints(spec, value, setup = setup)
        })
      )
    }
  }
  for(map in list(NULL, list(type = "sqrt", scale = 2.5))){
    for(share in list(prior("beta", list(1.5, .5)), prior("beta", list(2, .7)))){
      spec <- .prior_scale_product_spec(
        offset = .2, scale = -1.7, factor = prior("t", list(0, 1, 3)),
        multiplier = share, sources = list(), map = map
      )
      setup <- .prior_scale_product_breakpoint_setup(spec)
      distances <- (values - spec$offset) / spec$scale
      expect_identical(
        .prior_scale_product_breakpoints_values(spec, distances, setup),
        lapply(distances, .prior_scale_product_breakpoints, spec = spec, setup = setup)
      )
    }
  }

  # The batched quadrature computes the nodes and the value-independent terms
  # of the split integrand once per distinct interval: every value's integral
  # is identical to its integral computed alone, and to the integral of the
  # unsplit integrand.
  x <- seq(-4, 4, length.out = 61)
  spec <- list(additive_mean = 0, additive_sd = 0, product_mean = 0, product_sd = 1,
               multiplier = prior("normal", list(0, 1), list(0, Inf)), bounds = c(0, Inf))
  plans <- list(
    .prior_conditional_normal_density_plan(spec, x),
    .prior_scale_product_density_plan(.prior_scale_product_spec(
      offset = 0, scale = 1, factor = prior("t", list(0, 1, 3)),
      multiplier = prior("beta", list(2, 2)), sources = list()
    ), x)
  )
  for(plan in plans){
    expect_true(plan$batch)
    batched <- .prior_density_quadrature_batch(plan$integrand, plan$breakpoints)
    expect_false(anyNA(batched))
    alone <- vapply(seq_along(plan$breakpoints), function(j){
      .prior_density_quadrature_batch(
        list(shared = plan$integrand$shared,
             value  = function(shared, index) plan$integrand$value(shared, rep(j, length(index)))),
        plan$breakpoints[j]
      )
    }, numeric(1))
    expect_identical(batched, alone)
    unsplit <- function(nodes, index) plan$integrand$value(plan$integrand$shared(nodes), index)
    expect_identical(.prior_density_quadrature_batch(unsplit, plan$breakpoints), batched)
  }
})

test_that("route-product display grids batch the values of leaves with a singular integrand", {

  # X = A + B * S with A, B ~ N(0, 1) and a share S ~ Beta(1.5, .5), whose
  # density is infinite at 1 (an ordered level's cumulative share): the
  # batched bisection cannot reach 1e-8 next to the singularity, so plotted
  # curves take the per-value ordinates, while the display grid of a route
  # product batches all values with the ordinates' acceptance criterion
  # (relative error estimate at most 1e-4). Reference: s = 1 - w^2 removes the
  # singularity, f(x) = int_0^1 2 phi(x; 0, sqrt(1 + s^2)) s^(1/2) / B(1.5, .5) dw,
  # integrate() at rel.tol 1e-12.
  spec <- list(additive_mean = 0, additive_sd = 1, product_mean = 0, product_sd = 1,
               multiplier = prior("beta", list(1.5, .5)), bounds = c(0, 1), sources = list())
  route <- list(type = "conditional_normal", spec = spec, n_grid = 4096L)
  reference <- function(value){
    stats::integrate(function(w){
      s <- 1 - w^2
      2 * stats::dnorm(value, 0, sqrt(1 + s^2)) * sqrt(s) / beta(1.5, .5)
    }, 0, 1, rel.tol = 1e-12)$value
  }
  x <- seq(-4, 4, length.out = 41)
  expected <- vapply(x, reference, numeric(1))
  pieces <- 0L
  piece <- .prior_conditional_normal_piece
  local_mocked_bindings(.prior_conditional_normal_piece = function(...){
    pieces <<- pieces + 1L
    piece(...)
  })

  batched <- .prior_density_route_density(route, x, batch_singular = TRUE)
  expect_identical(pieces, 0L)
  expect_lte(max(abs(batched / expected - 1)), 1e-4)
  ordinates <- .prior_density_route_density(route, x)
  expect_gt(pieces, 0L)
  expect_lte(max(abs(ordinates / expected - 1)), 1e-8)

  pieces <- 0L
  grid <- .prior_linear_density_route_product(
    route, range = c(-4, 4), points = .prior_linear_density_empty_points(), n_grid = 41L
  )
  expect_identical(pieces, 0L)
  # the coalesced grid (knots rebuilt from the spacing, the curve normalized to
  # its continuous mass) interpolates the batched values
  expect_equal(grid$density$x, x, tolerance = 1e-12)
  expect_equal(grid$density$y / sum(grid$density$y), batched / sum(batched), tolerance = 1e-9)
})

test_that("route-product display grids take per-value ordinates for leaves with a strong singularity", {

  # Next to a density that is infinite at a finite bound with exponent below
  # 0.1 (e.g. a Beta share with shape < 0.1 there) the batched bisection's
  # error estimate is not reliable: values accepted at an estimate of 1e-4
  # were up to 4.3e-4 off. Such leaves are not batched, so their display-grid
  # values are the per-value ordinates. References: integrals over the Beta
  # term with s = u^(1 / a) on [0, 1/2] and s = 1 - v^(1 / b) on [1/2, 1],
  # which remove both bound singularities, integrate() at rel.tol 1e-13; the
  # bound checked is the ordinates' acceptance criterion 1e-4.
  beta_integral <- function(kernel, a, b){
    lower <- stats::integrate(function(u){
      s <- u^(1 / a)
      kernel(s) * (1 - s)^(b - 1) / (a * beta(a, b))
    }, 0, .5^a, rel.tol = 1e-13, subdivisions = 2000L)$value
    upper <- stats::integrate(function(v){
      s <- 1 - v^(1 / b)
      kernel(s) * s^(a - 1) / (b * beta(a, b))
    }, 0, .5^b, rel.tol = 1e-13, subdivisions = 2000L)$value
    lower + upper
  }
  check_leaf <- function(route, x, kernel, a, b){
    plan <- if(route$type == "scale_product"){
      .prior_scale_product_density_plan(route$spec, x, singular = TRUE)
    }else{
      .prior_conditional_normal_density_plan(route$spec, x, singular = TRUE)
    }
    expect_false(plan$batch)
    expect_null(plan$integrand)
    grid <- .prior_density_route_density(route, x, batch_singular = TRUE)
    expect_identical(grid, .prior_density_route_density(route, x))
    expected <- vapply(x, function(value) beta_integral(function(s) kernel(s, value), a, b),
                       numeric(1))
    expect_lte(max(abs(grid / expected - 1)), 1e-4)
  }

  # an ordered level with a half-normal total and Dirichlet(0.05) allocations
  # of four coefficients: the total times its Beta(0.05, 0.15) share, on the
  # 1024-value display grid of the product range (its 4th and 6th values were
  # 3.0e-4 and 4.3e-4 off)
  df <- data.frame(y = seq_len(10), f = ordered(rep(letters[1:5], 2), levels = letters[1:5]))
  formula_info <- JAGS_formula(y ~ f, "mu", data = df, prior_list = list(
    intercept = prior("normal", list(0, 1)),
    f = prior_ordered(prior("normal", list(0, 1), list(0, Inf)),
                      allocation = prior("dirichlet", list(alpha = rep(.05, 4))))
  ))
  level <- .prior_linear_combination_density(
    list(mu_f = formula_info$prior_list[["mu_f"]]),
    stats::setNames(c(1, 0, 0, 0), paste0("mu_f[", 1:4, "]"))
  )
  route <- .prior_density_route_from_adaptive(attr(level, "adaptive_evaluation", exact = TRUE))
  expect_identical(route$type, "scale_product")
  expect_equal(unlist(route$spec$multiplier$parameters), c(alpha = .05, beta = .15))
  range <- .prior_linear_density_range(level)
  x <- seq(range[1L], range[2L], length.out = 1024L)[c(2L, 4L, 6L, 50L, 400L)]
  check_leaf(route, x, function(s, value) 2 * stats::dnorm(value / s) / s, .05, .15)

  # a conditional-normal leaf A + B S, A ~ N(0, .3), B ~ N(2, 1),
  # S ~ Beta(0.05, 0.5) (1.4e-4 and 1.6e-4 off at -1.75 and 4 / 3), and a
  # scale product L S with L ~ gamma(2, 2) and S ~ Beta(0.3, 0.05), strong at
  # its upper bound (2.8e-4 off at 0.002)
  normal_leaf <- list(type = "conditional_normal", n_grid = 1024L, spec = list(
    additive_mean = 0, additive_sd = .3, product_mean = 2, product_sd = 1,
    multiplier = prior("beta", list(.05, .5)), bounds = c(0, 1), sources = list()
  ))
  check_leaf(normal_leaf, seq(-4, 4, length.out = 97)[c(1, 27, 28, 45, 49, 53, 65, 66, 97)],
             function(s, value) stats::dnorm(value, 2 * s, sqrt(.3^2 + s^2)), .05, .5)
  product_leaf <- list(type = "scale_product", n_grid = 1024L, spec = .prior_scale_product_spec(
    offset = 0, scale = 1, factor = prior("gamma", list(2, 2)),
    multiplier = prior("beta", list(.3, .05)), sources = list()
  ))
  check_leaf(product_leaf, c(.002, .01, .05, .2, 1, 2.5),
             function(s, value) stats::dgamma(value / s, 2, 2) / s, .3, .05)

  # a singular leaf whose exponents are at least 0.1 keeps the batch
  for(shape in list(c(.1, .3), c(.5, 1.5))){
    batched <- .prior_scale_product_density_plan(.prior_scale_product_spec(
      offset = 0, scale = 1, factor = prior("normal", list(0, 1), list(0, Inf)),
      multiplier = prior("beta", as.list(shape)), sources = list(),
      map = list(type = "sqrt", scale = 2)
    ), c(.1, 1), singular = TRUE)
    expect_false(batched$batch)
    expect_false(is.null(batched$integrand))
  }
  expect_false(.prior_density_strong_singularity(prior("beta", list(.1, .3))))
  expect_true(.prior_density_strong_singularity(prior("beta", list(.3, .05))))
  expect_true(.prior_density_strong_singularity(prior("gamma", list(.05, 1))))
  expect_false(.prior_density_strong_singularity(prior("gamma", list(.5, 1))))
  expect_false(.prior_density_strong_singularity(prior("normal", list(0, 1), list(0, Inf))))
})

test_that("row mixtures without atoms build their numerical grid only for grid consumers", {

  # Normal intercept + x * slope * sigma (sigma half-normal) over 12 distinct
  # design rows: every row's product component has an exact display grid, so
  # the mixture's grid is deferred until a grid field is read. The deferred
  # density equals the eagerly built one (the same builder) once built.
  sigma <- prior("normal", list(0, 1), list(0, Inf))
  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "sigma"
  priors <- list(mu_intercept = prior("normal", list(0, 1)), mu_x = slope, sigma = sigma)
  context <- .prior_density_build_context(priors, names(priors), n_grid = 256L)
  rows <- cbind(mu_intercept = 1, mu_x = stats::qnorm(seq(.05, .95, length.out = 12L)), sigma = 0)
  builds <- 0L
  build <- .prior_density_rows_grid
  local_mocked_bindings(.prior_density_rows_grid = function(...){
    builds <<- builds + 1L
    build(...)
  })
  eager_of <- function(rows){
    local_mocked_bindings(.prior_density_rows_atom_free = function(...) FALSE)
    .prior_density_from_context_rows(context, rows)
  }
  eager <- eager_of(rows)
  expect_false(inherits(eager$density, "prior_linear_density_deferred"))
  builds <- 0L

  deferred <- .prior_density_from_context_rows(context, rows)
  expect_s3_class(deferred$density, "prior_linear_density_deferred")
  expect_identical(attr(deferred, "adaptive_evaluation"), attr(eager, "adaptive_evaluation"))
  expect_identical(deferred$points, eager$points)
  expect_identical(deferred$n_grid, eager$n_grid)
  expect_identical(deferred$density$mass, eager$density$mass)

  # ordinates, exact heights and region probabilities use the route only
  for(value in c(-1, 0, .4)){
    expect_identical(prior_density_ordinate(deferred, value), prior_density_ordinate(eager, value))
    expect_identical(.prior_linear_density_height(deferred, value),
                     .prior_linear_density_height(eager, value))
  }
  region <- list(intervals = matrix(c(-.5, 1), 1L), indicator = function(x) x > -.5 & x < 1)
  expect_identical(
    .prior_density_route_region(.prior_density_route_from_adaptive(attr(deferred, "adaptive_evaluation")), region),
    .prior_density_route_region(.prior_density_route_from_adaptive(attr(eager, "adaptive_evaluation")), region)
  )
  # plots over a given range evaluate the route at the plotted values
  expect_identical(.prior_linear_density_to_plot_data(deferred, x_range = c(-3, 3)),
                   .prior_linear_density_to_plot_data(eager, x_range = c(-3, 3)))
  expect_identical(builds, 0L)

  # the first grid read builds the grid once; fields and the grid-reading
  # consumers then equal the eager density's
  expect_identical(deferred$density$x, eager$density$x)
  expect_identical(deferred$density[["y"]], eager$density$y)
  expect_identical(names(deferred$density), names(eager$density))
  expect_identical(.prior_linear_density_materialize(deferred), eager)
  expect_identical(.prior_linear_density_to_plot_data(deferred),
                   .prior_linear_density_to_plot_data(eager))
  expect_identical(.prior_linear_density_grid_height(deferred, .3),
                   .prior_linear_density_grid_height(eager, .3))
  expect_identical(builds, 1L)

  # rows with a point mass (a row without a weight is the point at zero, a
  # spike-and-slab slope has atoms) build their grid eagerly
  zero_row <- rbind(rows[1:2, ], c(mu_intercept = 0, mu_x = 0, sigma = 0))
  expect_false(inherits(.prior_density_from_context_rows(context, zero_row)$density,
                        "prior_linear_density_deferred"))
  spike_priors <- priors
  spike_priors$mu_x <- prior_spike_and_slab(prior("normal", list(0, 1)))
  attr(spike_priors$mu_x, "multiply_by") <- "sigma"
  spike_context <- .prior_density_build_context(spike_priors, names(spike_priors), n_grid = 256L)
  expect_false(inherits(.prior_density_from_context_rows(spike_context, rows[1:2, ])$density,
                        "prior_linear_density_deferred"))
})

# A small fitted object whose formula parameter is standardized (x), with one
# more parameter (sigma) outside the formula.
transformed_density_test_fit <- function(x_prior = prior("normal", list(0, 1))){

  scaled <- JAGS_formula(
    ~ x, "mu",
    data = data.frame(x = c(1, 2, 3.5, 4, 6, 8.5)),
    prior_list = list(intercept = prior("normal", list(0, 1)), x = x_prior),
    formula_scale = list(x = TRUE)
  )
  set.seed(4)
  n <- 40L
  posterior <- cbind(
    mu_intercept = stats::rnorm(n),
    mu_x         = stats::rnorm(n, .2, .2),
    sigma        = abs(stats::rnorm(n, 1, .1))
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
         monitor = colnames(posterior), sample = n),
    class = c("runjags", "BayesTools_fit")
  )
  attr(fit, "prior_list") <- c(scaled$prior_list, list(sigma = prior("gamma", list(2, 2))))
  attr(fit, "formula_design") <- list(mu = scaled$formula_design)
  attr(fit, "formula_scale") <- list(mu = scaled$formula_scale)
  fit <- attach_test_parameter_map(fit)
  .bt_attach_fit_contract(.bt_attach_draw_geometry(.bt_attach_parameter_map(fit)))
}

# The values of a density, without the record of the context it was built from
# (which holds the evaluation arguments, not the density).
transformed_density_payload <- function(density){

  values <- attributes(density)
  values[["adaptive_evaluation"]] <- NULL
  out <- unclass(density)
  attributes(out) <- values
  out
}

# Reference: .generate_transformed_prior_densities() without a memo, which is the
# pre-memo computation of every prior column of the list.
test_that("transformed prior densities are computed once per fit and only for the requested parameters", {

  fit <- transformed_density_test_fit()
  original <- BayesTools:::.prior_linear_combination_density
  calls <- 0L
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(...){
      calls <<- calls + 1L
      original(...)
    },
    .package = "BayesTools"
  )
  request <- function(parameters, n_prior_samples = 512L, fit. = fit){
    as_mixed_posteriors(fit., parameters, transform_scaled = TRUE, n_prior_samples = n_prior_samples)
  }
  densities <- function(samples) .bt_meta_get(samples, "prior_densities")

  # only the prior columns of the requested parameters are built
  calls <- 0L
  intercept <- request("mu_intercept")
  expect_identical(names(densities(intercept)), "mu_intercept")
  first_cost <- calls
  expect_gt(first_cost, 0L)

  # the same request again, and a subset of it, build nothing
  calls <- 0L
  again <- request("mu_intercept")
  expect_identical(calls, 0L)
  expect_identical(densities(again)$mu_intercept, densities(intercept)$mu_intercept)

  # a wider request builds only the columns that are new
  calls <- 0L
  both <- request(c("mu_intercept", "mu_x"))
  expect_identical(names(densities(both)), c("mu_intercept", "mu_x"))
  expect_gt(calls, 0L)
  expect_lt(calls, 2L * first_cost)
  expect_identical(densities(both)$mu_intercept, densities(intercept)$mu_intercept)
  calls <- 0L
  request(c("mu_x", "mu_intercept"))
  expect_identical(calls, 0L)

  # every column equals the unmemoized computation from the same inputs
  prior_list <- attr(fit, "prior_list")
  column_names <- colnames(as.matrix(fit$mcmc[[1L]]))
  reference <- .generate_transformed_prior_densities(
    prior_list, column_names, attr(fit, "formula_scale"), n_grid = 512L
  )
  all_columns <- request(c("mu_intercept", "mu_x", "sigma"))
  expect_identical(names(densities(all_columns)), c("mu_intercept", "mu_x", "sigma"))
  for(column in names(reference)){
    expect_identical(
      transformed_density_payload(densities(all_columns)[[column]]),
      transformed_density_payload(reference[[column]]),
      info = column
    )
  }

  # any change of an input is recomputed and never answered from the memo
  calls <- 0L
  request("mu_intercept", n_prior_samples = 1024L)
  expect_gt(calls, 0L)
  # (the memo keeps two input sets: use the original one so that the changed
  # prior below displaces the grid of 1024 instead)
  calls <- 0L
  request(c("mu_intercept", "mu_x"))
  expect_identical(calls, 0L)
  changed <- fit
  attr(changed, "prior_list")$mu_x <- prior("normal", list(0, 2))
  calls <- 0L
  shifted <- request(c("mu_intercept", "mu_x"), fit. = changed)
  expect_gt(calls, 0L)
  shifted_reference <- .generate_transformed_prior_densities(
    attr(changed, "prior_list"), column_names, attr(fit, "formula_scale"), n_grid = 512L
  )
  for(column in c("mu_intercept", "mu_x")){
    expect_identical(
      transformed_density_payload(densities(shifted)[[column]]),
      transformed_density_payload(shifted_reference[[column]]),
      info = column
    )
  }
  expect_false(identical(
    transformed_density_payload(densities(shifted)$mu_x),
    transformed_density_payload(densities(both)$mu_x)
  ))
  # ... while the original inputs are still answered from their own entry
  calls <- 0L
  request(c("mu_intercept", "mu_x"))
  expect_identical(calls, 0L)

  # what was derived from a fitted map is dropped when the map's tables are
  # replaced, so a fit with edited metadata never reuses the old densities
  edited <- fit
  edited_map <- attr(edited, "parameter_map")
  edited_map$quantities$display_label[[1L]] <- "relabelled"
  attr(edited, "parameter_map") <- edited_map
  request(c("mu_intercept", "mu_x"))
  calls <- 0L
  request(c("mu_intercept", "mu_x"))
  expect_identical(calls, 0L)
  request(c("mu_intercept", "mu_x"), fit. = edited)
  expect_gt(calls, 0L)

  # the memo is bounded and does not consume random numbers
  set.seed(9)
  seed <- .Random.seed
  for(n_grid in seq(600L, by = 100L, length.out = 8L)){
    request("mu_x", n_prior_samples = n_grid)
  }
  expect_identical(.Random.seed, seed)
  memo <- BayesTools:::.bt_prior_density_memo(fit)
  expect_identical(length(memo$entries), BayesTools:::.bt_prior_density_memo_limit())

  # a fit without a parameter map has nothing to memoize against and computes
  expect_null(BayesTools:::.bt_prior_density_memo(structure(list(), class = "BayesTools_fit")))
  calls <- 0L
  .generate_transformed_prior_densities(prior_list, column_names, attr(fit, "formula_scale"), n_grid = 512L)
  .generate_transformed_prior_densities(prior_list, column_names, attr(fit, "formula_scale"), n_grid = 512L)
  expect_gt(calls, first_cost)
})

# Reference: the unmemoized computation of the same inputs (no memo), which is
# what a recomputation after eviction has to reproduce bit for bit.
test_that("the transformed prior density memo keeps the two most recently used input sets of a fit", {

  fit <- transformed_density_test_fit()
  original <- BayesTools:::.prior_linear_combination_density
  calls <- 0L
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(...){
      calls <<- calls + 1L
      original(...)
    },
    .package = "BayesTools"
  )
  columns <- c("mu_intercept", "mu_x")
  payloads <- function(densities){
    stats::setNames(lapply(columns, function(column){
      transformed_density_payload(densities[[column]])
    }), columns)
  }
  # the density builds and the density payloads of a request with one grid
  built <- function(n_prior_samples){
    calls <<- 0L
    samples <- as_mixed_posteriors(fit, columns, transform_scaled = TRUE, n_prior_samples = n_prior_samples)
    list(calls = calls, payload = payloads(.bt_meta_get(samples, "prior_densities")))
  }
  prior_list <- attr(fit, "prior_list")
  column_names <- colnames(as.matrix(fit$mcmc[[1L]]))
  reference <- function(n_grid){
    payloads(.generate_transformed_prior_densities(
      prior_list, column_names, attr(fit, "formula_scale"), n_grid = n_grid
    ))
  }

  expect_identical(BayesTools:::.bt_prior_density_memo_limit(), 2L)
  memo <- BayesTools:::.bt_prior_density_memo(fit)
  set.seed(11)
  seed <- .Random.seed

  # three distinct input sets (grids A, B, C); a hit moves its set to the front
  a1 <- built(512L)
  b1 <- built(768L)
  expect_gt(a1$calls, 0L)
  expect_gt(b1$calls, 0L)
  expect_length(memo$entries, 2L)
  expect_identical(built(768L)$calls, 0L)
  expect_identical(built(512L)$calls, 0L)

  # the third set evicts the least recently used one: B, although A is older,
  # because A was used last
  c1 <- built(1024L)
  expect_gt(c1$calls, 0L)
  expect_length(memo$entries, 2L)
  expect_identical(built(512L)$calls, 0L)
  expect_identical(built(1024L)$calls, 0L)

  # a re-request after eviction builds the same densities again, and evicts the
  # set used least recently (now A)
  b2 <- built(768L)
  expect_identical(b2$calls, b1$calls)
  expect_identical(b2$payload, b1$payload)
  expect_length(memo$entries, 2L)
  expect_identical(built(1024L)$calls, 0L)
  a2 <- built(512L)
  expect_identical(a2$calls, a1$calls)
  expect_identical(a2$payload, a1$payload)
  expect_length(memo$entries, 2L)

  # with and without eviction the values are those of the unmemoized computation
  expect_identical(a1$payload, reference(512L))
  expect_identical(b1$payload, reference(768L))
  expect_identical(c1$payload, reference(1024L))
  expect_identical(a2$payload, reference(512L))
  expect_identical(b2$payload, reference(768L))
  expect_identical(.Random.seed, seed)
})
