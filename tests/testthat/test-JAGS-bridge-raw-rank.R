skip_if_not_test_profile("unit")

test_that("R133 raw affine dependence refuses bounded bridge coordinates before sampling", {
  x2 <- seq(1, 4, length.out = 200)
  dependent <- cbind(x1 = 2 * x2 + 1, x2 = x2)
  sampler_calls <- 0L
  testthat::local_mocked_bindings(bridge_sampler = function(...){
    sampler_calls <<- sampler_calls + 1L
    structure(list(logml = -1, niter = 1L, method = "normal", mcse_logml = .1), class = "bridge")
  }, .package = "bridgesampling")
  for(bounds in list(c(0, 10), c(0, Inf), c(-Inf, 10), c(-Inf, Inf))){
    posterior <- dependent
    attr(posterior, "lb") <- c(x1 = bounds[1], x2 = bounds[1])
    attr(posterior, "ub") <- c(x1 = bounds[2], x2 = bounds[2])
    expect_error(.bt_JAGS_bridge_check_coordinate_rank(posterior), "rank-deficient (rank 1 of 2", fixed = TRUE)
    expect_error(JAGS_bridgesampling(coda::as.mcmc(dependent),
      log_posterior = function(parameters, data) 0, data = list(),
      prior_list = list(x1 = prior("normal", list(0, 1), truncation = list(bounds[1], bounds[2]))),
      add_parameters = "x2", add_bounds = list(lb = c(x2 = bounds[1]), ub = c(x2 = bounds[2]))),
      "rank-deficient (rank 1 of 2", fixed = TRUE)
  }
  expect_identical(sampler_calls, 0L)
  full_rank <- dependent
  full_rank[, "x2"] <- dependent[c(seq(2, 200, 2), seq(1, 199, 2)), "x2"]
  expect_s3_class(JAGS_bridgesampling(coda::as.mcmc(full_rank),
    log_posterior = function(parameters, data) 0, data = list(),
    prior_list = list(x1 = prior("normal", list(0, 1), truncation = list(0, 10))),
    add_parameters = "x2", add_bounds = list(lb = c(x2 = 0), ub = c(x2 = 10))), "BayesTools_marglik")
  expect_identical(sampler_calls, 1L)
})

test_that("R133 bridge rank checks preserve transformed-only and ordinary controls", {
  t <- seq(1, 2, length.out = 200)
  transformed_dependence <- cbind(x1 = exp(t), x2 = exp(2 * t))
  attr(transformed_dependence, "lb") <- c(x1 = 0, x2 = 0)
  attr(transformed_dependence, "ub") <- c(x1 = Inf, x2 = Inf)
  expect_error(.bt_JAGS_bridge_check_coordinate_rank(transformed_dependence), "rank-deficient (rank 1 of 2", fixed = TRUE)
  full_rank <- cbind(x1 = t, x2 = sin(t))
  expect_true(.bt_JAGS_bridge_check_coordinate_rank(full_rank))
  expect_true(.bt_JAGS_bridge_check_coordinate_rank(cbind(x1 = t, constant = 1)))
  expect_error(.bt_JAGS_bridge_check_coordinate_rank(full_rank[1:2, ]), "2 posterior draws span only rank 1 of 2", fixed = TRUE)
  expect_error(.bt_JAGS_bridge_check_coordinate_rank(cbind(x1 = t, x2 = 2 * t)), "rank-deficient (rank 1 of 2", fixed = TRUE)
})
