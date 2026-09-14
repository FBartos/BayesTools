skip_if_not_test_profile("unit")

test_that("cumulative selection initializations stay strictly positive", {

  cumulative <- prior_weightfunction(
    side = "one-sided",
    steps = .05,
    weights = wf_cumulative(c(1e-300, 1e-300))
  )

  expect_warning(
    direct_inits <- JAGS_get_inits(
      prior_list = list(bias = cumulative),
      chains = 1,
      seed = 1
    ),
    "deterministic, order-one rescaling"
  )
  for(chain_inits in direct_inits){
    expect_true(is.finite(chain_inits$omega_ratio))
    expect_true(chain_inits$omega_ratio > 0)
    expect_true(chain_inits$omega_ratio < 1)
  }

  mixture <- prior_mixture(list(
    prior_none(),
    cumulative
  ))
  expect_warning(
    mixture_inits <- JAGS_get_inits(
      prior_list = list(bias = mixture),
      chains = 1,
      seed = 1
    ),
    "deterministic, order-one rescaling"
  )
  for(chain_inits in mixture_inits){
    expect_true(is.finite(chain_inits$omega_ratio_component_2))
    expect_true(chain_inits$omega_ratio_component_2 > 0)
    expect_true(chain_inits$omega_ratio_component_2 < 1)
  }

  expect_warning(
    component_inits <- BayesTools:::.JAGS_init.weightfunction(
      cumulative,
      component_id = 2
    ),
    "deterministic, order-one rescaling"
  )
  expect_true(is.finite(component_inits$omega_ratio_component_2))
  expect_true(component_inits$omega_ratio_component_2 > 0)
  expect_true(component_inits$omega_ratio_component_2 < 1)
})

test_that("independent weight functions emit JAGS inits for free bins", {

  omega <- prior_weightfunction(
    side = "one-sided",
    steps = c(.025, .05),
    weights = wf_independent(prior("gamma", list(2, 1)))
  )
  log_omega <- prior_weightfunction(
    side = "one-sided",
    steps = .05,
    weights = wf_independent(prior("normal", list(0, 1)), scale = "log_omega")
  )

  set.seed(1)
  omega_inits <- JAGS_get_inits(list(bias = omega), chains = 1, seed = 1)[[1]]
  expect_equal(length(omega_inits$omega), 3L)
  expect_true(is.na(omega_inits$omega[[1L]]))
  expect_true(all(is.finite(omega_inits$omega[-1L])))
  expect_true(all(omega_inits$omega[-1L] > 0))

  log_inits <- JAGS_get_inits(list(bias = log_omega), chains = 1, seed = 1)[[1]]
  expect_equal(length(log_inits$log_omega), 2L)
  expect_true(is.na(log_inits$log_omega[[1L]]))
  expect_true(is.finite(log_inits$log_omega[[2L]]))

  mixture <- prior_mixture(list(prior_none(), omega))
  mixture_inits <- JAGS_get_inits(list(bias = mixture), chains = 1, seed = 1)[[1]]
  expect_equal(length(mixture_inits$omega_component_2), 3L)
  expect_true(is.na(mixture_inits$omega_component_2[[1L]]))
  expect_true(all(is.finite(mixture_inits$omega_component_2[-1L])))
})

test_that("Gamma initialization fallback preserves representable median proportions", {

  set.seed(1)
  expect_warning(
    initialization <- BayesTools:::.JAGS_positive_gamma_initialization(
      rep(.001, 4),
      "test latent Gamma variables"
    ),
    "deterministic, order-one rescaling"
  )
  expect_identical(initialization, rep(1, 4))

  set.seed(1)
  expect_warning(
    expect_error(
      BayesTools:::.JAGS_positive_gamma_initialization(
        c(1e-323, 2e-323),
        "test latent Gamma variables"
      ),
      "medians are not representable"
    ),
    "deterministic, order-one rescaling"
  )
})
