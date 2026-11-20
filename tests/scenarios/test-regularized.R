if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized")
# testthat::test_file("tests/scenarios/test-regularized.R")

### Description
# Which predictors of fertility retain an association after shrinkage? Compare
# four priors on the five slopes using data from 47 Swiss provinces. Fertility
# and the predictors are standardized, and the residual SD is fixed at 1.
# For each prior, BayesTools and RStan estimate the same regression. The table
# and figure show each engine's posterior for the slopes and their shared scale.

testthat::test_that("Shrinkage of Swiss fertility predictors", {

  dat <- as.data.frame(scale(datasets::swiss))
  slopes <- ~ Agriculture + Examination + Education + Catholic + Infant.Mortality
  X <- model.matrix(slopes, dat)[, -1]

  likelihood <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu[i], 1)
    }
  }"
  jags_data <- list(N = nrow(dat), y = dat$Fertility)
  stan_data <- list(N = nrow(dat), K = ncol(X), X = X, y = dat$Fertility)
  intercept <- prior("normal", list(0, 1))
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))

  # The intercept has its own prior. Each of the five slopes shares the global
  # scale; under a horseshoe each also has its own local scale.
  horseshoe_data <- c(stan_data, list(J = ncol(X), unit = seq_len(ncol(X))))

  # Match coefficients by predictor name. RStan's beta vector follows the
  # column order in X. Both fits report their shared global scale as well.
  parameters <- c("Intercept", colnames(X), "Global scale")
  jags_columns <- setNames(c("mu_intercept", paste0("mu_", colnames(X)), "mu_shrinkage_global"), parameters)
  stan_columns <- setNames(c("alpha", sprintf("beta[%d]", seq_len(ncol(X))), "tau"), parameters)
  posterior_draws <- function(jags, stan) {

    list(BayesTools = scenario_draws(coda::as.mcmc.list(jags), jags_columns),
         RStan      = scenario_draws(rstan::As.mcmc.list(stan), stan_columns))
  }

  ### Ridge ----
  # Each slope has a Normal prior with the shared scale as its SD.
  ridge <- scenario_fit("jags_ridge", {
    JAGS_fit(likelihood, jags_data, prior_list = list(),
             formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
             formula_prior_list = list(mu = list(
               intercept = intercept, shrinkage = prior_regularized("ridge", scale = global))),
             chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
  })
  ridge_stan <- scenario_fit("stan_ridge", {
    rstan::sampling(scenario_stan_model("ridge.stan"), data = stan_data,
                    chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
                    control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
  }, cache_version = scenario_stan_version("ridge.stan"))
  scenario_comparison("comparison-ridge", parameters, main = "Swiss fertility: ridge",
                      posteriors = posterior_draws(ridge, ridge_stan), show_table = TRUE)

  ### Lasso ----
  # Each slope has a Laplace prior with the shared scale.
  lasso <- scenario_fit("jags_lasso", {
    JAGS_fit(likelihood, jags_data, prior_list = list(),
             formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
             formula_prior_list = list(mu = list(
               intercept = intercept, shrinkage = prior_regularized("lasso", scale = global))),
             chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
  })
  lasso_stan <- scenario_fit("stan_lasso", {
    rstan::sampling(scenario_stan_model("lasso.stan"), data = stan_data,
                    chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
                    control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
  }, cache_version = scenario_stan_version("lasso.stan"))
  scenario_comparison("comparison-lasso", parameters, main = "Swiss fertility: lasso",
                      posteriors = posterior_draws(lasso, lasso_stan), show_table = TRUE)

  ### Horseshoe ----
  # The shared global and slope-specific local scales determine each slope's SD.
  horseshoe <- scenario_fit("jags_horseshoe", {
    JAGS_fit(likelihood, jags_data, prior_list = list(),
             formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
             formula_prior_list = list(mu = list(
               intercept = intercept, shrinkage = prior_regularized("horseshoe", scale = global))),
             chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
  })
  horseshoe_stan <- scenario_fit("stan_horseshoe", {
    rstan::sampling(scenario_stan_model("horseshoe.stan"), data = horseshoe_data,
                    chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
                    control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
  }, cache_version = scenario_stan_version("horseshoe.stan"))
  scenario_comparison("comparison-horseshoe", parameters, main = "Swiss fertility: horseshoe",
                      posteriors = posterior_draws(horseshoe, horseshoe_stan), show_table = TRUE)

  ### Horseshoe with a finite slab ----
  # A shared Student-t slab (scale 2, 4 degrees of freedom) limits large effects.
  finite_slab <- scenario_fit("jags_regularized_horseshoe", {
    JAGS_fit(likelihood, jags_data, prior_list = list(),
             formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
             formula_prior_list = list(mu = list(
               intercept = intercept, shrinkage = prior_regularized("horseshoe", scale = global,
                                                                    slab_scale = 2, slab_df = 4))),
             chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
  })
  finite_slab_stan <- scenario_fit("stan_regularized_horseshoe", {
    rstan::sampling(scenario_stan_model("finite-slab-horseshoe.stan"), data = horseshoe_data,
                    chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
                    control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
  }, cache_version = scenario_stan_version("finite-slab-horseshoe.stan"))
  scenario_comparison("comparison-regularized_horseshoe", parameters,
                      main = "Swiss fertility: horseshoe with a finite slab",
                      posteriors = posterior_draws(finite_slab, finite_slab_stan), show_table = TRUE)
})
