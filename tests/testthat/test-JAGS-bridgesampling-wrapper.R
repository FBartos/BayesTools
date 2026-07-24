skip_if_not_test_profile("unit")

test_that("JAGS_bridgesampling preserves repeated bridge estimates", {

  set.seed(11)
  posterior <- coda::as.mcmc(matrix(
    stats::rnorm(1000),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  prior_list <- list(mu = prior("normal", list(0, 1)))

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = prior_list,
    repetitions = 2,
    maxiter = 1000
  )

  expect_s3_class(result, "bridge_list")
  expect_length(result[["logml"]], 2L)
  expect_length(result[["niter"]], 2L)
  expect_true(all(is.finite(result[["logml"]])))
  expect_true(all(result[["niter"]] <= 1000L))
  expect_null(attr(result, "warning"))
})

test_that("JAGS_bridgesampling checks repeated iteration limits collectively", {

  set.seed(11)
  posterior <- coda::as.mcmc(matrix(
    stats::rnorm(1000),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))
  prior_list <- list(mu = prior("normal", list(0, 1)))

  result <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = function(parameters, data) 0,
    data = list(),
    prior_list = prior_list,
    repetitions = 2,
    maxiter = 1
  )

  expect_s3_class(result, "bridge_list")
  expect_length(result[["niter"]], 2L)
  expect_true(any(result[["niter"]] > 1L))
  expect_identical(
    attr(result, "warning"),
    paste(
      "Marginal likelihood could not be estimated within the maximum number",
      "of iterations and might be more variable than usual."
    )
  )
})
