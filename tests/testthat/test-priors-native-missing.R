skip_if_not_test_profile("unit")

test_that("R133 native scalar methods preserve NA and NaN separately", {
  specifications <- list(
    list("moment", list(tau = 1, order = 1), list(lower = .2, upper = 3)),
    list("invmoment", list(tau = 1, order = 1, df = 3), list(lower = .2, upper = 3)),
    list("invgamma", list(shape = 3, scale = 2), list(lower = .2, upper = 3))
  )
  observations <- c(NA_real_, NaN, 1, -Inf, Inf)
  probabilities <- c(NA_real_, NaN, .5, 0, 1)
  for(specification in specifications){
    for(truncation in list(NULL, specification[[3L]])){
      p <- if(is.null(truncation)) prior(specification[[1L]], specification[[2L]]) else
        prior(specification[[1L]], specification[[2L]], truncation = truncation)
      outputs <- c(lapply(list(pdf, lpdf, mpdf, mlpdf, cdf, ccdf, mcdf, mccdf),
                          function(method) method(p, observations)),
                   lapply(list(quant, mquant), function(method) method(p, probabilities)))
      expect_true(all(vapply(outputs, function(value) identical(value[1:2], c(NA_real_, NaN)), logical(1))))
      expect_true(all(vapply(outputs, function(value) all(!is.na(value[3:5])), logical(1))))
      expect_equal(cdf(p, observations[3:5]) + ccdf(p, observations[3:5]), rep(1, 3), tolerance = 1e-12)
    }
  }
})

test_that("R133 native log flags retain missing values and finite reference identities", {
  missing <- c(NA_real_, NaN)
  for(family in c("moment", "invmoment", "invgamma")){
    parameters <- switch(family, moment = list(0, 1, 1), invmoment = list(0, 1, 1, 3), invgamma = list(3, 2))
    for(log_flag in c(FALSE, TRUE)){
      density <- get(paste0(".d", family, "_prior"), asNamespace("BayesTools"))
      expect_identical(do.call(density, c(list(missing), parameters, list(log = log_flag))), missing)
      for(lower_tail in c(FALSE, TRUE)){
        for(operation in c("p", "q")){
          method <- get(paste0(".", operation, family, "_prior"), asNamespace("BayesTools"))
          expect_identical(do.call(method, c(list(missing), parameters,
            list(lower.tail = lower_tail, log.p = log_flag))), missing)
        }
      }
    }
  }
  expect_equal(pdf(prior("moment", list(tau = 1)), 1), stats::dnorm(1), tolerance = 1e-12)
  expect_equal(cdf(prior("moment", list(tau = 1)), 1), .5 + .5 * stats::pchisq(1, 3), tolerance = 1e-12)
  expect_equal(lpdf(prior("invmoment", list(tau = 1, df = 3)), 1), -lgamma(1.5) - 1, tolerance = 1e-12)
  expect_equal(cdf(prior("invmoment", list(tau = 1, df = 3)), 1), .5 + .5 * stats::pgamma(1, 1.5, lower.tail = FALSE), tolerance = 1e-12)
  p <- prior("invgamma", list(shape = 3, scale = 2))
  expect_equal(lpdf(p, 1), 3 * log(2) - lgamma(3) - 2, tolerance = 1e-12)
  expect_equal(cdf(p, 1), stats::pgamma(2, 3, lower.tail = FALSE), tolerance = 1e-12)
})

test_that("R133 generic truncated Normal CDF preserves both missing values", {
  p <- prior("normal", list(0, 1), truncation = list(-1, 2))
  q <- c(NA_real_, NaN, .5, -Inf, Inf)
  reference <- (stats::pnorm(.5) - stats::pnorm(-1)) / (stats::pnorm(2) - stats::pnorm(-1))
  for(method in list(cdf, mcdf, ccdf, mccdf)) expect_identical(method(p, q)[1:2], c(NA_real_, NaN))
  for(method in list(cdf, mcdf)) expect_equal(method(p, q)[3:5], c(reference, 0, 1), tolerance = 1e-12)
  for(method in list(ccdf, mccdf)) expect_equal(method(p, q)[3:5], c(1 - reference, 1, 0), tolerance = 1e-12)
})
