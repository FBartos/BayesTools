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
