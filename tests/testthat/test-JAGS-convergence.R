skip_if_not_test_profile("unit")

.mock_convergence_fit <- function(chain_1, chain_2){

  fit <- list(
    mcmc = coda::mcmc.list(
      coda::mcmc(chain_1),
      coda::mcmc(chain_2)
    ),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"

  return(fit)
}


test_that("JAGS_check_convergence preserves internal indicator and inclusion substrings", {

  parameters <- c(
    "theta_indicator_effect",
    "theta_indicator_effect[1]",
    "theta_inclusion_effect",
    "theta_inclusion_effect[1]"
  )
  for(parameter in parameters){
    chain_1 <- matrix(
      rep(-5, 100),
      ncol = 1,
      dimnames = list(NULL, parameter)
    )
    chain_2 <- matrix(
      rep(5, 100),
      ncol = 1,
      dimnames = list(NULL, parameter)
    )
    prior_list <- list(prior("normal", list(0, 1)))
    names(prior_list) <- sub("\\[.*$", "", parameter)

    convergence <- JAGS_check_convergence(
      .mock_convergence_fit(chain_1, chain_2),
      prior_list = prior_list,
      max_Rhat = 1.05,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL
    )

    expect_false(convergence, info = parameter)
    expect_match(attr(convergence, "errors"), "R-hat", info = parameter)
  }
})


test_that("JAGS_check_convergence recognizes scalar and indexed auxiliary suffixes", {

  set.seed(1)
  chain_1 <- cbind(
    mu = stats::rnorm(100),
    mu_indicator = 0,
    "mu_indicator[1]" = 0,
    mu_inclusion = 0,
    "mu_inclusion[1]" = 0
  )
  chain_2 <- cbind(
    mu = stats::rnorm(100),
    mu_indicator = 1,
    "mu_indicator[1]" = 1,
    mu_inclusion = 1,
    "mu_inclusion[1]" = 1
  )
  fit <- .mock_convergence_fit(chain_1, chain_2)
  prior_list <- list(mu = prior("normal", list(0, 1)))

  expect_true(JAGS_check_convergence(
    fit,
    prior_list = prior_list,
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    check_indicators = FALSE
  ))

  convergence <- JAGS_check_convergence(
    fit,
    prior_list = prior_list,
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    check_indicators = TRUE
  )

  expect_false(convergence)
  expect_match(attr(convergence, "errors"), "R-hat")
})

test_that("sampled constants are not assessable but point priors are structural", {

  constant_chain <- matrix(
    rep(1, 20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )
  sampled <- JAGS_check_convergence(
    .mock_convergence_fit(constant_chain, constant_chain),
    prior_list = list(mu = prior("normal", list(0, 1))),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL
  )
  expect_false(sampled)
  expect_match(attr(sampled, "errors"), "R-hat.*not assessable")
  sampled_diagnostics <- attr(sampled, "diagnostics")
  expect_s3_class(
    sampled_diagnostics,
    "BayesTools_convergence_diagnostics"
  )
  expect_equal(sampled_diagnostics$state, "not_assessable")

  structural <- JAGS_check_convergence(
    .mock_convergence_fit(constant_chain, constant_chain),
    prior_list = list(mu = prior("point", list(1))),
    max_Rhat = 1.05,
    min_ESS = 10,
    max_error = 0.01,
    max_SD_error = 0.05
  )
  expect_true(structural)
  expect_equal(
    attr(structural, "diagnostics")$state,
    "structural_constant"
  )
})

test_that("explicit convergence monitors distinguish omitted parameters", {

  set.seed(42)
  chain_1 <- cbind(mu = rnorm(100), tau = rep(1, 100))
  chain_2 <- cbind(mu = rnorm(100), tau = rep(1, 100))
  fit <- .mock_convergence_fit(chain_1, chain_2)
  priors <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1))
  )

  selected <- JAGS_check_convergence(
    fit,
    prior_list = priors,
    max_Rhat = NULL,
    min_ESS = 1,
    max_error = NULL,
    max_SD_error = NULL,
    monitor = "mu"
  )
  expect_true(selected)
  diagnostics <- attr(selected, "diagnostics")
  expect_equal(
    setNames(diagnostics$state, diagnostics$parameter),
    c(mu = "assessable", tau = "not_requested")
  )

  all_parameters <- JAGS_check_convergence(
    fit,
    prior_list = priors,
    max_Rhat = NULL,
    min_ESS = 1,
    max_error = NULL,
    max_SD_error = NULL
  )
  expect_false(all_parameters)
  expect_equal(
    attr(all_parameters, "diagnostics")$state,
    c("assessable", "not_assessable")
  )
  expect_error(
    JAGS_check_convergence(
      fit,
      prior_list = priors,
      monitor = "misspelled"
    ),
    "requested convergence monitor 'misspelled' is not available"
  )
})

test_that("a base convergence monitor selects every indexed element", {

  chain_1 <- cbind(
    mu = 1:20,
    "beta[1]" = 1:20,
    "beta[2]" = 20:1
  )
  chain_2 <- cbind(
    mu = 2:21,
    "beta[1]" = 2:21,
    "beta[2]" = 21:2
  )
  result <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2),
    prior_list = list(
      mu = prior("normal", list(0, 1)),
      beta = prior("mnormal", list(mean = 0, sd = 1, K = 2))
    ),
    max_Rhat = NULL,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    monitor = "beta"
  )

  expect_true(result)
  diagnostics <- attr(result, "diagnostics")
  expect_equal(
    diagnostics$parameter[diagnostics$state == "assessable"],
    c("beta[1]", "beta[2]")
  )
  expect_equal(
    diagnostics$state[diagnostics$parameter == "mu"],
    "not_requested"
  )
})

test_that("empty convergence selections do not claim convergence", {

  set.seed(43)
  chain_1 <- cbind(mu = rnorm(20))
  chain_2 <- cbind(mu = rnorm(20))
  result <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2),
    prior_list = list(mu = prior("normal", list(0, 1))),
    monitor = character()
  )

  expect_identical(as.vector(result), logical())
  expect_s3_class(
    attr(result, "diagnostics"),
    "BayesTools_convergence_diagnostics"
  )
  expect_equal(nrow(attr(result, "diagnostics")), 0L)
})

test_that("empty available convergence parameters are vacuously TRUE", {

  set.seed(44)
  chain_1 <- cbind(mu = rnorm(20))
  chain_2 <- cbind(mu = rnorm(20))
  # Drop every sample column via add_parameters removal with no structural priors.
  result <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2),
    prior_list = list(),
    add_parameters = "mu",
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL
  )

  expect_true(result)
  expect_equal(nrow(attr(result, "diagnostics")), 0L)
})

test_that("not-assessable diagnostics require an explicit opt-in to ignore", {

  constant_chain <- matrix(
    rep(1, 20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )
  result <- JAGS_check_convergence(
    .mock_convergence_fit(constant_chain, constant_chain),
    prior_list = list(mu = prior("normal", list(0, 1))),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    allow_not_assessable = TRUE
  )

  expect_true(result)
  expect_equal(
    attr(result, "diagnostics")$state,
    "not_assessable"
  )
  expect_null(attr(result, "errors"))
})

test_that("one chain cannot satisfy an enabled R-hat criterion", {

  set.seed(44)
  chain <- cbind(mu = rnorm(20))
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(chain)),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"
  result <- JAGS_check_convergence(
    fit,
    prior_list = list(mu = prior("normal", list(0, 1))),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL
  )

  expect_false(result)
  expect_match(attr(result, "errors"), "R-hat.*not assessable")
  expect_equal(attr(result, "diagnostics")$state, "not_assessable")
})
