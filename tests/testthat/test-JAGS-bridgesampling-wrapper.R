skip_if_not_test_profile("unit")

.mock_bridge_sampler <- function(...){

  arguments <- list(...)
  repetitions <- arguments[["repetitions"]]
  if(is.null(repetitions)){
    repetitions <- 1L
  }
  maxiter <- arguments[["maxiter"]]
  niter <- if(maxiter <= 1L){
    c(2L, rep.int(1L, repetitions - 1L))
  }else{
    rep.int(7L, repetitions)
  }

  structure(
    list(
      logml = seq_len(repetitions),
      niter = niter
    ),
    class = "bridge_list"
  )
}

test_that("JAGS_bridgesampling preserves repeated bridge estimates", {

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
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

  testthat::local_mocked_bindings(
    bridge_sampler = .mock_bridge_sampler,
    .package = "bridgesampling"
  )
  posterior <- coda::as.mcmc(matrix(
    seq_len(20),
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
