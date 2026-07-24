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
  mu <- stats::rnorm(100)
  chain_1 <- cbind(
    mu = mu,
    mu_indicator = 0,
    "mu_indicator[1]" = 0,
    mu_inclusion = 0,
    "mu_inclusion[1]" = 0
  )
  chain_2 <- cbind(
    mu = mu,
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
