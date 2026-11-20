if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized")
source(file.path(.scenario_config()$root, "helpers", "shrinkage-plot.R"), local = TRUE)
source(file.path(.scenario_config()$root, "helpers", "regularized-bridge.R"), local = TRUE)
# testthat::test_file("tests/scenarios/test-regularized.R")

### Description
# How do BayesTools and RStan shrink the five fertility slopes in 47 Swiss
# provinces? Compare the slope distributions and the evidence for each prior
# in a standardized regression with known residual SD 1.

testthat::test_that("Shrinkage of Swiss fertility predictors", {

  dat <- as.data.frame(scale(datasets::swiss))
  slopes <- ~ Agriculture + Examination + Education + Catholic + Infant.Mortality
  X <- model.matrix(slopes, dat)[, -1]
  y <- dat$Fertility

  ### The same regression in BayesTools and RStan ----
  likelihood <- "model {
    for (i in 1:N) { y[i] ~ dnorm(mu[i], 1) }
  }"
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))
  stan_data <- list(N = nrow(X), K = ncol(X), X = X, y = y)
  # Each horseshoe slope has its own local scale.
  horseshoe_data <- c(stan_data, list(J = ncol(X), unit = seq_len(ncol(X))))

  fit_pair <- function(name, shrinkage, file, data = stan_data) {
    jags <- scenario_fit(paste0("jags_", name), {
      JAGS_fit(likelihood, data = list(y = y, N = length(y)), prior_list = list(),
        formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
        formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)), shrinkage = shrinkage)),
        chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
    }, cache_version = regularized_fit_version(likelihood, deparse(slopes), dat, y, shrinkage))
    stan <- scenario_fit(paste0("stan_", name), {
      rstan::sampling(scenario_stan_model(file), data = data,
        chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011, refresh = 0,
        control = list(adapt_delta = 0.995, max_treedepth = 14))
    }, cache_version = regularized_fit_version(file, scenario_stan_version(file), data))
    list(name = name, jags = jags, stan = stan, stan_data = data)
  }

  ### What the tables and figures compare ----
  slope_priors <- setNames(paste0("mu_", colnames(X)), colnames(X))
  parameters <- c("Intercept", colnames(X), "Global scale")
  jags_columns <- setNames(c("mu_intercept", slope_priors, "mu_shrinkage_global"), parameters)
  stan_columns <- setNames(c("alpha", sprintf("beta[%d]", seq_len(ncol(X))), "tau"), parameters)
  show_fertility <- function(pair, prior_label) {
    draws <- list(BayesTools = scenario_draws(pair$jags$mcmc, jags_columns),
                  RStan = scenario_draws(rstan::As.mcmc.list(pair$stan), stan_columns))
    title <- paste("Swiss fertility:", prior_label)
    scenario_comparison(paste0("comparison-", pair$name), parameters, posteriors = draws,
                        main = title, show_table = TRUE)
    scenario_plot(paste0("distributions-", pair$name),
                  scenario_shrinkage_plot(pair$jags, draws, slope_priors, title))
  }

  ### Four priors on the same five slopes ----
  ridge <- fit_pair("ridge", prior_regularized("ridge", global), "ridge.stan")
  show_fertility(ridge, "ridge")

  lasso <- fit_pair("lasso", prior_regularized("lasso", global), "lasso.stan")
  show_fertility(lasso, "lasso")

  horseshoe <- fit_pair("horseshoe", prior_regularized("horseshoe", global), "horseshoe.stan", horseshoe_data)
  show_fertility(horseshoe, "horseshoe")

  finite_slab <- fit_pair("regularized_horseshoe", prior_regularized("horseshoe", global, slab_scale = 2, slab_df = 4),
                          "finite-slab-horseshoe.stan", horseshoe_data)
  show_fertility(finite_slab, "horseshoe with a finite slab")

  ### Marginal likelihoods and log Bayes factors relative to ridge ----
  # The Examination local scale mixes slowly in the finite-slab fit. A fresh
  # set of chains can move its log evidence by several hundredths.
  models <- list(Ridge = ridge, Lasso = lasso, Horseshoe = horseshoe, "Finite-slab horseshoe" = finite_slab)
  bridges <- lapply(models, function(pair) regularized_bridge(pair$name, pair$jags, pair$stan, pair$stan_data))
  scenario_text("model-comparison", regularized_evidence_table(bridges))
})
