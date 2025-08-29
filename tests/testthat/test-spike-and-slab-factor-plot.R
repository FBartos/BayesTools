# Test for plot_posterior with spike-and-slab factor priors
# This test addresses the issue: "plot_posterior does not show differences from grand mean for spike-and-slab factor priors"

test_that("plot_posterior handles spike-and-slab factor priors correctly", {
  
  skip_on_cran()
  
  set.seed(1)
  data_formula <- data.frame(
    x_fac = factor(rep(c("A", "B", "C"), 10), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(30, ifelse(data_formula$x_fac == "B", -0.2, 0.4)),
    N = 30
  )
  
  formula_list <- list(mu = ~ x_fac - 1)
  formula_data_list <- list(mu = data_formula)
  prior_list <- list(sigma = prior("lognormal", list(0, 1)))
  model_syntax <- "model{\nfor(i in 1:N){\n  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n}\n}"
  
  # Test regular factor prior (baseline)
  formula_prior_list_normal <- list(
    mu = list("x_fac" = prior_factor("mnormal", contrast = "meandif", list(0, 1)))
  )
  
  fit_normal <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list, silent = TRUE,
    formula_list = formula_list, formula_data_list = formula_data_list, 
    formula_prior_list = formula_prior_list_normal,
    n_iter = 100, n_burnin = 50, n_chains = 2
  )
  
  # Test spike-and-slab factor prior (the issue)
  formula_prior_list_spike <- list(
    mu = list("x_fac" = prior_spike_and_slab(prior_factor("mnormal", contrast = "meandif", list(0, 1))))
  )
  
  fit_spike <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list, silent = TRUE,
    formula_list = formula_list, formula_data_list = formula_data_list, 
    formula_prior_list = formula_prior_list_spike,
    n_iter = 100, n_burnin = 50, n_chains = 2
  )
  
  this_parameter <- JAGS_parameter_names(parameters = "x_fac", formula_parameter = "mu")
  this_posterior_normal <- as_mixed_posteriors(model = fit_normal, parameters = this_parameter)
  this_posterior_spike <- as_mixed_posteriors(model = fit_spike, parameters = this_parameter)
  
  # Verify that both prior types are correctly detected as factors
  prior_list_normal <- attr(this_posterior_normal[[this_parameter]], "prior_list")
  prior_list_spike <- attr(this_posterior_spike[[this_parameter]], "prior_list")
  
  prior_list_normal_simplified <- BayesTools:::.simplify_prior_list(prior_list_normal)
  prior_list_spike_simplified <- BayesTools:::.simplify_prior_list(prior_list_spike)
  
  # Test the factor detection logic (this was the bug)
  has_factor_normal <- any(sapply(prior_list_normal_simplified, function(p) {
    if(is.prior.factor(p)){
      TRUE
    } else if(is.prior.spike_and_slab(p) && is.prior.factor(p$variable)){
      TRUE
    } else {
      FALSE
    }
  }))
  
  has_factor_spike <- any(sapply(prior_list_spike_simplified, function(p) {
    if(is.prior.factor(p)){
      TRUE
    } else if(is.prior.spike_and_slab(p) && is.prior.factor(p$variable)){
      TRUE
    } else {
      FALSE
    }
  }))
  
  expect_true(has_factor_normal, "Normal factor prior should be detected as factor")
  expect_true(has_factor_spike, "Spike-and-slab factor prior should be detected as factor")
  
  # Test that plot_posterior works for both types (this was failing before the fix)
  expect_silent(plot_posterior(this_posterior_normal, parameter = this_parameter, prior = TRUE))
  expect_silent(plot_posterior(this_posterior_spike, parameter = this_parameter, prior = TRUE))
  
  # Additional verification: the spike-and-slab should have the correct structure
  spike_prior <- prior_list_spike[[1]]
  expect_true(is.prior.spike_and_slab(spike_prior), "First prior should be spike-and-slab")
  expect_true(is.prior.factor(spike_prior$variable), "Variable part should be a factor prior")
  expect_true(is.prior.meandif(spike_prior$variable), "Variable part should be meandif contrast")
})