skip_if_not_test_profile("unit")

test_that("cumulative selection initializations stay strictly positive", {

  cumulative <- prior_weightfunction(
    side = "one-sided",
    steps = .05,
    weights = wf_cumulative(c(1e-300, 1e-300))
  )

  direct_inits <- JAGS_get_inits(
    prior_list = list(bias = cumulative),
    chains = 2,
    seed = 1
  )
  for(chain_inits in direct_inits){
    expect_true(all(is.finite(chain_inits$eta)))
    expect_true(all(chain_inits$eta > 0))
  }

  mixture <- prior_mixture(list(
    prior_none(),
    cumulative
  ))
  mixture_inits <- JAGS_get_inits(
    prior_list = list(bias = mixture),
    chains = 2,
    seed = 1
  )
  for(chain_inits in mixture_inits){
    expect_true(all(is.finite(chain_inits$eta_component_2)))
    expect_true(all(chain_inits$eta_component_2 > 0))
  }

  component_inits <- BayesTools:::.JAGS_init.weightfunction(
    cumulative,
    component_id = 2
  )
  expect_true(all(is.finite(component_inits$eta_component_2)))
  expect_true(all(component_inits$eta_component_2 > 0))
})
