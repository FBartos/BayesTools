skip_if_not_test_profile("unit")

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
