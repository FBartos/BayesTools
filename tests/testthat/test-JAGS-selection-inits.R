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
    expect_true(all(is.finite(chain_inits$eta)))
    expect_true(all(chain_inits$eta > 0))
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
    expect_true(all(is.finite(chain_inits$eta_component_2)))
    expect_true(all(chain_inits$eta_component_2 > 0))
  }

  expect_warning(
    component_inits <- BayesTools:::.JAGS_init.weightfunction(
      cumulative,
      component_id = 2
    ),
    "deterministic, order-one rescaling"
  )
  expect_true(all(is.finite(component_inits$eta_component_2)))
  expect_true(all(component_inits$eta_component_2 > 0))
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
