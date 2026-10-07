skip_if_not_test_profile("unit")

test_that("ordinary truncated densities preserve available natural tail range results", {

  normal <- prior("normal", list(0, 1))
  truncated <- prior("normal", list(0, 1), list(lower = 0, upper = Inf))
  expect_warning(expect_identical(pdf(normal, 40), 0), NA)
  expect_warning(expect_identical(pdf(truncated, 40), 0), NA)
  expect_equal(lpdf(truncated, 40), stats::dnorm(40, log = TRUE) + log(2),
               tolerance = 2e-13)
})

test_that("inverse-gamma coordinates retain scale and subnormal root precision", {

  scale <- 1e-310
  p <- prior("invgamma", list(shape = 3, scale = scale))
  expected <- -1 - log(2) - log(scale)
  expect_equal(lpdf(p, scale), expected, tolerance = 3e-13)
  ordinate <- prior_density_ordinate(p, scale)
  expect_true(ordinate$exact)
  expect_equal(ordinate$log_density, expected, tolerance = 3e-13)
  expect_equal(cdf(p, scale), stats::pgamma(1, 3, lower.tail = FALSE), tolerance = 2e-15)
  expect_equal(ccdf(p, scale), stats::pgamma(1, 3), tolerance = 2e-15)
  values <- quant(p, c(.1, .5, .9))
  expect_true(all(is.finite(values) & values > 0))
  expect_equal(values / scale, 1 / stats::qgamma(c(.1, .5, .9), 3, lower.tail = FALSE),
               tolerance = 3e-12)
  expect_warning(pdf(p, scale), class = "BayesTools_numerical_range_limit")
  # Exact shape-1 identity. qgamma's positive subnormal root rounds to 2*minsub;
  # using its rounded logarithm would move the final ordinary value by 22%.
  expect_equal(.qinvgamma_prior(-744, 1, scale, lower.tail = FALSE, log.p = TRUE),
               exp(log(scale) + 744), tolerance = 2e-13)
  expect_equal(.pinvgamma_prior(1, 1, scale, lower.tail = FALSE, log.p = TRUE),
               log(scale), tolerance = 2e-13)
  expect_equal(.pinvgamma_prior(1, 1, scale, lower.tail = FALSE) / scale, 1,
               tolerance = 2e-13)
})

test_that("tiny Gamma shapes use available log tails and declare the remaining limit", {

  a <- .Machine$double.xmin * .Machine$double.eps
  expect_identical(.pinvgamma_prior(2, a, 1), a)
  # Two direct 460/520-digit references at exact binary64 a and radial r=.5/1.
  expect_equal(.pinvgamma_prior(c(2, 1), a, 1, log.p = TRUE),
               c(-745.02029479342605, -745.95700388038331), tolerance = 3e-13)
  expect_equal(.dinvgamma_prior(.5, a, 1, log = TRUE),
               log(a) - 2 + log(2), tolerance = 3e-13)
  warning <- NULL
  values <- withCallingHandlers(.pinvgamma_prior(c(.5, NA_real_, NaN), a, 1, log.p = TRUE),
    BayesTools_numerical_unavailable = function(condition){
      warning <<- condition
      invokeRestart("muffleWarning")
    })
  expect_true(is.nan(values[1L]))
  expect_identical(values[2:3], c(NA_real_, NaN))
  expect_s3_class(warning, "BayesTools_numerical_condition")
  expect_identical(warning$indices, 1L)
  expect_null(warning$call)
  expect_identical(warning$family, "invgamma")
  expect_identical(warning$requested_scale, "log")
  # Analytic upper-tail bound establishes this correctly rounded natural zero.
  expect_warning(expect_identical(.pinvgamma_prior(.5, a, 1), 0), NA)
  # Domain-invalid low-level values remain distinct from valid unavailability.
  expect_warning(expect_true(is.nan(.pinvgamma_prior(1, -1, 1))), NA)
  expect_identical(.dinvgamma_prior(1, -1, 1, log = TRUE), -Inf)
})

test_that("nonlocal scale coordinates and high-order quantiles remain finite", {

  expect_warning(expect_identical(.qmoment_prior(.5, 0, .125, 1), 0), NA)
  expect_warning(expect_identical(.qinvmoment_prior(.5, 0, 1, 1000, 3), 0), NA)
  p <- prior("invmoment", list(tau = 1e308, order = 1, df = 1))
  expect_equal(cdf(p, 2e154), .739750061093476738, tolerance = 2e-14)
  expect_equal(ccdf(p, 2e154), .260249938906523262, tolerance = 2e-14)
  expected <- log(.5) - .25 - lgamma(.5) - log(2e154)
  expect_equal(lpdf(p, 2e154), expected, tolerance = 3e-13)
  high <- prior("invmoment", list(tau = 1, order = 1000, df = 3))
  q <- quant(high, c(.1, .5, .9))
  expect_true(all(is.finite(q)))
  expect_identical(q[2L], 0)
  expect_equal(q[1L], -q[3L], tolerance = 2e-13)
  expect_equal(cdf(high, q[c(1L, 3L)]), c(.1, .9), tolerance = 3e-12)
  moment <- prior("moment", list(tau = 1e308, order = 2))
  expect_true(is.finite(moment$parameters$mode))
  expect_true(is.finite(quant(moment, .9)))
  largest <- prior("moment", list(tau = .Machine$double.xmax))
  expect_true(is.finite(largest$parameters$mode))
  expect_equal(log(largest$parameters$mode),
               (log(2) + log(.Machine$double.xmax)) / 2, tolerance = 2e-15)
})

test_that("inverse-moment final magnitudes retain original probability and parameter precision", {

  orders <- c(100, 1000, .Machine$integer.max, 1e308)
  # Independent 460/520-digit exact-binary Gamma references under the certified
  # leading identity. The order-100 control retains its normal qgamma root.
  distances <- c(5.014348381131724187, 5.001442219288464544,
                 5.000000000671967197, 5)
  for(i in seq_along(orders)){
    for(lower_tail in c(TRUE, FALSE)){
      expected <- if(lower_tail) c(-distances[i], distances[i]) else c(distances[i], -distances[i])
      expect_equal(.qinvmoment_prior(c(.1, .9), 0, 1, orders[i], 1, lower.tail = lower_tail),
                   expected, tolerance = 3e-13)
      expect_equal(.qinvmoment_prior(log(c(.1, .9)), 0, 1, orders[i], 1,
                                    lower.tail = lower_tail, log.p = TRUE),
                   expected, tolerance = 3e-13)
    }
  }
  expect_equal(.qinvmoment_prior(.1, .25, 4, 1e308, 1), -9.75, tolerance = 3e-13)
  expect_equal(.qinvmoment_prior(.9, 2, 1e308, 1e308, 1), 5e154, tolerance = 3e-13)
  expect_equal(.qinvmoment_prior(.1, 0, 1e-310, 1e308, 1) / sqrt(1e-310), -5,
               tolerance = 3e-13)
  expect_identical(.qinvmoment_prior(c(0, .5, 1), .25, 1, 1e308, 1), c(-Inf, .25, Inf))
  expect_identical(.qinvmoment_prior(c(-Inf, -log(2), 0), .25, 1, 1e308, 1, log.p = TRUE),
                   c(-Inf, .25, Inf))
  expect_warning(overflow <- .qinvmoment_prior(.1, 0, 1, 1e308, .001),
                 class = "BayesTools_numerical_range_limit")
  expect_identical(overflow, -Inf)
  expect_error(.prior_numerical_finite(overflow, "invmoment"), class = "BayesTools_prior_rng_unavailable")
  expect_warning(unavailable <- .qinvmoment_prior(.1, 0, 1, 1e308,
    .Machine$double.xmin * .Machine$double.eps), class = "BayesTools_numerical_unavailable")
  expect_true(is.nan(unavailable))
  expect_warning(collapsed <- .qinvmoment_prior(.1, 1e308, 1, 1e308, 1),
                 class = "BayesTools_numerical_unavailable")
  expect_true(is.nan(collapsed))
  # Native acceptance of large finite integer orders does not change the public cap.
  expect_error(prior("invmoment", list(tau = 1, order = 1e308, df = 1)), "'order'", fixed = TRUE)
})

test_that("far-tail and central truncations retain normalization and conditional quantiles", {

  gamma <- prior("gamma", list(shape = 1, rate = 1), list(lower = 800, upper = Inf))
  expect_equal(lpdf(gamma, 801), -1, tolerance = 2e-13)
  expect_equal(pdf(gamma, 801), exp(-1), tolerance = 2e-13)
  expect_equal(cdf(gamma, 801), -expm1(-1), tolerance = 2e-13)
  expect_equal(ccdf(gamma, 801), exp(-1), tolerance = 2e-13)
  expect_equal(quant(gamma, c(0, .5, 1)), c(800, 800 + log(2), Inf), tolerance = 2e-13)
  moment <- prior("moment", list(mode = .5), list(lower = 14, upper = Inf))
  log_mass <- -log(2) + stats::pchisq(14^2 / .125, 3, lower.tail = FALSE, log.p = TRUE)
  reference <- 2 * log(14.001) + stats::dnorm(14.001, sd = sqrt(.125), log = TRUE) - log(.125)
  expect_equal(lpdf(moment, 14.001), reference - log_mass, tolerance = 4e-12)
  expect_equal(quant(moment, .5), sqrt(.125 * stats::qchisq(log_mass, 3, lower.tail = FALSE, log.p = TRUE)),
               tolerance = 2e-13)
  priors <- list(
    prior("moment", list(tau = .125), list(lower = -1e-8, upper = 1e-8)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = -.01, upper = .01))
  )
  for(p in priors){
    q <- quant(p, c(0, .25, .5, .75, 1))
    expect_true(all(is.finite(q)))
    expect_identical(q[c(1L, 3L, 5L)], c(p$truncation$lower, 0, p$truncation$upper))
    expect_equal(q[2L], -q[4L], tolerance = 3e-12)
    expect_equal(cdf(p, q[c(2L, 4L)]), c(.25, .75), tolerance = 3e-11)
    expect_equal(cdf(p, c(NA_real_, NaN)), c(NA_real_, NaN))
  }
  shifted <- prior("moment", list(tau = .125, location = .25), list(lower = -.75, upper = 1.25))
  expect_identical(quant(shifted, .5), .25)
  near_symmetric <- prior("moment", list(tau = .125, location = 1e-300), list(lower = -1, upper = 1))
  expect_warning(unavailable <- quant(near_symmetric, .5), class = "BayesTools_numerical_unavailable")
  expect_true(is.nan(unavailable))
})

test_that("unresolved normalization refuses once before evaluator and bridge work", {

  p <- prior("invgamma", list(shape = .Machine$double.xmin * .Machine$double.eps, scale = 1),
             list(lower = .4, upper = .5))
  expect_error(.prior_simple_lpdf_evaluator(p), class = "BayesTools_numerical_unavailable")
  expect_error(lpdf(p, .45), class = "BayesTools_numerical_unavailable")
  expect_error(JAGS_marglik_priors(c(theta = .45), list(theta = p)),
               class = "BayesTools_numerical_unavailable")
  ordinate <- prior_density_ordinate(p, .45)
  expect_identical(ordinate$behavior, "regular")
  expect_false(ordinate$exact)
  expect_true(is.na(ordinate$log_density))
  expect_match(ordinate$reason, "truncation mass", fixed = TRUE)
  expect_error(.prior_nonlocal_log_interval_mass(prior("moment", list(mode = .5)), c(0, 1), 2),
               "equal lengths", fixed = TRUE)
})

test_that("ordinary Gamma subnormal truncation refuses unreliable backend normalization", {

  shape <- .Machine$double.xmin * .Machine$double.eps
  p <- prior("gamma", list(shape = shape, rate = 1), list(lower = 1, upper = 2))
  expect_identical(p$parameters$shape, shape)
  expect_error(.prior_simple_lpdf_evaluator(p), class = "BayesTools_numerical_unavailable")
  expect_error(pdf(p, c(1.5, 2)), class = "BayesTools_numerical_unavailable")
  expect_error(lpdf(p, c(1.5, 2)), class = "BayesTools_numerical_unavailable")
  expect_error(cdf(p, 1.5), class = "BayesTools_numerical_unavailable")
  expect_error(quant(p, .5), class = "BayesTools_numerical_unavailable")
  expect_error(rng(p, 4), class = "BayesTools_prior_rng_unavailable")
  expect_error(.generate_prior_sample_matrix(list(theta = p), 4), class = "BayesTools_prior_rng_unavailable")
  expect_error(JAGS_get_inits(list(theta = p), chains = 1, seed = 631),
               class = "BayesTools_numerical_unavailable")
  expect_error(JAGS_marglik_priors(c(theta = 1.5), list(theta = p)),
               class = "BayesTools_numerical_unavailable")
  ordinate <- prior_density_ordinate(p, 1.5)
  expect_identical(ordinate$behavior, "regular")
  expect_false(ordinate$exact)
  expect_true(is.na(ordinate$log_density))
  expect_identical(quant(p, c(0, 1)), c(1, 2))
  expect_identical(cdf(p, c(1, 2)), c(0, 1))
  whole <- prior("gamma", list(shape = shape, rate = 1))
  expect_identical(.prior_simple_log_C(whole), 0)
  expect_identical(.prior_C(whole), 1)
  # Independent 50/80-digit E1 limit references. At these adjacent normal
  # shapes the finite-shape correction is far below binary64 resolution.
  for(a in c(.Machine$double.xmin, 2 * .Machine$double.xmin)){
    supported <- prior("gamma", list(shape = a, rate = 1), list(lower = 1, upper = 2))
    expect_equal(pdf(supported, c(1.5, 2)), c(.8725390239209256, .3969162523528322),
                 tolerance = 2e-12)
    expect_equal(lpdf(supported, c(1.5, 2)), c(-.13634789934882265, -.9240299718006036),
                 tolerance = 2e-12)
    expect_equal(cdf(supported, 1.5), .7001522459316267, tolerance = 2e-12)
  }
})

test_that("nonlocal missing pairs and unrepresentable narrow quantiles are explicit", {

  p <- prior("moment", list(tau = .125))
  out <- .prior_nonlocal_log_interval_mass(p, c(NA_real_, NaN, 0), c(NaN, NA_real_, 0))
  expect_identical(out, c(NA_real_, NaN, -Inf))
  narrow <- prior("moment", list(tau = .125), list(lower = 1, upper = 1 + .Machine$double.eps))
  condition <- NULL
  value <- tryCatch(withCallingHandlers(quant(narrow, .5),
    BayesTools_numerical_unavailable = function(warning){
      condition <<- warning
      if(inherits(warning, "warning")) invokeRestart("muffleWarning")
    }), BayesTools_numerical_unavailable = function(error){
      condition <<- error
      NaN
    })
  expect_true(is.nan(value))
  expect_s3_class(condition, "BayesTools_numerical_unavailable")
  expect_identical(quant(narrow, c(0, 1)), c(1, 1 + .Machine$double.eps))
})

test_that("uniform budgets remain unchanged and finite-required consumers are strict", {

  set.seed(604)
  ig <- rng(prior("invgamma", list(shape = 3, scale = 1e-310)), 8)
  after_ig <- .Random.seed
  set.seed(604)
  stats::runif(8)
  expect_identical(.Random.seed, after_ig)
  expect_true(all(is.finite(ig) & ig > 0))
  p <- prior("invmoment", list(tau = 1, order = 1000, df = 3))
  set.seed(605)
  draws <- rng(p, 8)
  after <- .Random.seed
  set.seed(605)
  uniforms <- matrix(stats::runif(16), nrow = 2)
  expect_identical(.Random.seed, after)
  expect_identical(sign(draws), ifelse(uniforms[2L, ] < .5, -1, 1))
  # Check the declared signed Gamma law for both the ordinary-root and
  # certified leading-log-root branches.
  expect_equal(.pinvmoment_prior(abs(draws), 0, 1, 1000, 3, lower.tail = FALSE) * 2,
               uniforms[1L, ], tolerance = 4e-12)
  set.seed(608)
  extreme <- .rinvmoment_prior(8, 0, 1, 1e308, 1)
  extreme_after <- .Random.seed
  set.seed(608)
  extreme_uniforms <- matrix(stats::runif(16), nrow = 2)
  expect_identical(.Random.seed, extreme_after)
  # At this exact accepted native order, the certified correction is far below
  # binary64 resolution; the original size draw supplies the final magnitude.
  expect_equal(extreme, ifelse(extreme_uniforms[2L, ] < .5, -1, 1) / extreme_uniforms[1L, ],
               tolerance = 3e-13)
  truncated <- prior("moment", list(mode = .5), list(lower = 14, upper = Inf))
  set.seed(606)
  bounded <- rng(truncated, 8)
  bounded_after <- .Random.seed
  set.seed(606)
  u <- stats::runif(8)
  expect_identical(.Random.seed, bounded_after)
  expect_equal(bounded, quant(truncated, u), tolerance = 0)
  range_prior <- prior("invmoment", list(tau = 1, df = 1e-310))
  set.seed(607)
  expect_warning(raw <- rng(range_prior, 8), class = "BayesTools_numerical_range_limit")
  range_after <- .Random.seed
  set.seed(607)
  stats::runif(16)
  expect_identical(.Random.seed, range_after)
  expect_true(all(is.infinite(raw)))
  expect_error(.prior_numerical_finite(raw, "invmoment"), class = "BayesTools_prior_rng_unavailable")
  expect_warning(expect_error(
    .generate_prior_sample_matrix(list(theta = range_prior), 8, seed = 607),
    class = "BayesTools_prior_rng_unavailable"), class = "BayesTools_numerical_range_limit")
  expect_warning(expect_error(
    JAGS_get_inits(list(theta = range_prior), chains = 1, seed = 607),
    class = "BayesTools_numerical_unavailable"), class = "BayesTools_numerical_range_limit")
  expect_warning(expect_error(
    density(range_prior, x_range = c(-10, 10), n_points = 16, n_samples = 8, force_samples = TRUE),
    class = "BayesTools_prior_rng_unavailable"), class = "BayesTools_numerical_range_limit")
})

test_that("bounded R sampling retains exact inverse-CDF sign boundaries", {

  requested <- integer()
  priors <- list(
    prior("moment", list(mode = .5, location = .25), list(lower = -.75, upper = 1.25)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = -.01, upper = .01)))
  testthat::with_mocked_bindings({
    for(p in priors){
      expect_warning(draws <- rng(p, 3), NA)
      expect_identical(draws, rep(p$parameters$location, 3))
    }
    near_symmetric <- prior("moment", list(tau = .125, location = 1e-300),
                            list(lower = -1, upper = 1))
    expect_error(rng(near_symmetric, 3), class = "BayesTools_prior_rng_unavailable")
  }, runif = function(n){
    requested <<- c(requested, n)
    rep(.5, n)
  }, .package = "stats")
  expect_identical(requested, c(3, 3, 3))
})

test_that("nonlocal moment integrals retain values without raw PDF tail warnings", {

  priors <- list(
    prior("moment", list(mode = .5), list(lower = -Inf, upper = 0)),
    prior("invmoment", list(mode = .5, df = 6), list(lower = -Inf, upper = 0)))
  # Genuine installed pre-correction values, with the same integration settings.
  expected <- list(c(-.564189583547719, .238096858059049),
                   c(-.621742035370958, .225691917118892))
  for(i in seq_along(priors)){
    expect_warning(values <- c(mean(priors[[i]]), sd(priors[[i]])), NA)
    expect_equal(values, expected[[i]], tolerance = 1e-12)
  }
  unresolved <- prior("invmoment", list(tau = 1, order = 1000,
    df = .Machine$double.xmin * .Machine$double.eps), list(lower = -1, upper = 1))
  expect_error(mean(unresolved), class = "BayesTools_numerical_unavailable")
  expect_error(sd(unresolved), class = "BayesTools_numerical_unavailable")
})

test_that("inverse moments retain finite logs and refuse unresolved offset arithmetic", {

  lognormal <- prior("lognormal", list(meanlog = 0, sdlog = 40))
  moment <- .prior_density_inverse_moment(lognormal, 1024)
  expect_identical(moment$value, Inf)
  expect_identical(moment$log_value, 800)
  expect_true(moment$available)
  spec <- list(additive_mean = 0, additive_sd = 0, product_mean = 0, product_sd = 1,
               multiplier = lognormal, bounds = c(0, Inf), sources = list())
  ordinate <- .prior_conditional_normal_offset_ordinate(spec, 0, 1024)
  expect_identical(ordinate$behavior, "regular")
  expect_true(ordinate$exact)
  expect_equal(ordinate$log_density, 799.0810614667953, tolerance = 2e-13)
  huge <- .prior_density_inverse_moment(prior("lognormal", list(0, 1.5e154)), 1024)
  expect_true(is.finite(huge$log_value))
  expect_identical(huge$log_value, (1.5e154 * .5) * 1.5e154)
  expect_false(huge$available)
  spec$multiplier <- prior("lognormal", list(0, 1.5e154))
  unavailable <- .prior_conditional_normal_offset_ordinate(spec, 0, 1024)
  expect_identical(unavailable$behavior, "regular")
  expect_false(unavailable$exact)
  expect_true(is.na(unavailable$log_density))
  mapped <- .prior_scale_product_inverse_moment(list(multiplier = prior("beta", list(2^60, 2^60)),
                                                   map = list(scale = 1)), 1024)
  expect_false(mapped$available)
  expect_true(is.na(mapped$log_value))
  expect_match(mapped$reason, "shape shift", fixed = TRUE)
  moderate <- .prior_scale_product_inverse_moment(list(multiplier = prior("beta", list(2, 3)),
                                                     map = list(scale = 2)), 1024)
  expect_true(moderate$available)
  expect_equal(moderate$log_value, lbeta(1.5, 3) - lbeta(2, 3) - log(2) / 2, tolerance = 2e-15)
})

test_that("public inference refuses unavailable inverse-moment ordinates", {

  slope <- prior("normal", list(0, 1))
  attr(slope, "multiply_by") <- "s"
  context <- .prior_density_context(
    list(b = slope, s = prior("lognormal", list(0, 1.5e154))), "b", n_grid = 1024)
  density <- .prior_density_from_context_rows(
    context, matrix(c(1, 2), ncol = 1, dimnames = list(NULL, "b")))
  ordinate <- prior_density_ordinate(density, 0)
  expect_identical(ordinate$behavior, "regular")
  expect_false(ordinate$exact)
  expect_true(is.na(ordinate$log_density))
  status <- prior_ordinate_status(density, 0)
  expect_false(status$eligible)
  expect_identical(status$condition, "BayesTools_inexact_ordinate")
  posterior <- structure(seq(-1, 1, length.out = 200),
    class = c("marginal_posterior.simple", "marginal_posterior", "numeric"))
  posterior <- .bt_meta_set(posterior, "prior_density", density)
  posterior <- .bt_meta_set(posterior, "atoms", posterior_atom_attribute())
  expect_error(hypothesis_BF(posterior = posterior, hypothesis = "theta = 0", parameter = "theta"),
               class = "BayesTools_inexact_ordinate")
})
test_that("nonlocal JAGS starts require usable normalized log density without redrawing", {

  requested <- integer()
  priors <- list(prior("moment", list(mode = .5, location = .25), list(lower = -.75, upper = 1.25)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = -.01, upper = .01)))
  testthat::with_mocked_bindings({
    for(p in priors){
      expect_identical(rng(p, 1L), p$parameters$location)
      expect_identical(quant(p, .5), p$parameters$location)
      condition <- tryCatch(.JAGS_init.simple(p, "theta"), error = identity)
      expect_s3_class(condition, "BayesTools_numerical_unavailable")
      expect_identical(condition$operation, "initialization")
      expect_identical(condition$indices, 1L)
      expect_identical(condition$values, p$parameters$location)
    }
  }, runif = function(n){requested <<- c(requested, n); rep(.5, n)}, .package = "stats")
  expect_identical(requested, c(1L, 1, 1L, 1))
  for(p in priors){
    set.seed(616)
    expected <- rng(p, 1L)
    seed <- .Random.seed
    set.seed(616)
    observed <- .JAGS_init.simple(p, "theta")$theta
    expect_identical(observed, expected)
    expect_identical(.Random.seed, seed)
    expect_true(is.finite(lpdf(p, observed)))
  }
})
test_that("nonlocal finite starts retain range and unresolved normalized-density refusals", {

  priors <- list(prior("moment", list(tau = 1)),
    prior("invmoment", list(tau = 1, order = 1000, df = .Machine$double.xmin * .Machine$double.eps), list(lower = -1, upper = 1)))
  values <- c(1e200, .5)
  for(i in seq_along(priors)){
    testthat::with_mocked_bindings({
      condition <- tryCatch(.JAGS_init.simple(priors[[i]], "theta"), error = identity)
      available <- .Call("BayesTools_native_range_environment", PACKAGE = "BayesTools")$available
      expect_s3_class(condition, if(i == 1L && available) "BayesTools_numerical_range_limit" else "BayesTools_numerical_unavailable")
      expect_identical(condition$operation, "initialization")
      expect_identical(condition$values, values[i])
      expect_identical(condition$indices, 1L)
      if(i == 1L && available){
        expect_identical(condition$log_density, -Inf)
        expect_s3_class(condition$parent, "BayesTools_numerical_range_limit")
      }
    }, rng = function(...) values[i])
  }
})
