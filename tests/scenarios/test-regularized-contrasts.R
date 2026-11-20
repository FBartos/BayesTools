if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized-contrasts")
source(file.path(.scenario_config()$root, "helpers", "shrinkage-plot.R"), local = TRUE)
source(file.path(.scenario_config()$root, "helpers", "regularized-bridge.R"), local = TRUE)
# testthat::test_file("tests/scenarios/test-regularized-contrasts.R")

### Description
# Do BayesTools and RStan agree on shrinkage of monthly differences in New
# York ozone? Each month effect is a deviation from the equally weighted mean
# of May through September at average wind. Compare the four priors under
# mean-difference and orthonormal coding, with standardized ozone and wind and
# known residual SD 1.

testthat::test_that("Monthly ozone effects under centered shrinkage priors", {

  ozone <- na.omit(airquality[c("Ozone", "Month", "Wind")])
  ozone$Month <- factor(ozone$Month, labels = month.abb[5:9])
  ozone$Ozone <- as.numeric(scale(ozone$Ozone))
  ozone$Wind  <- as.numeric(scale(ozone$Wind))
  months <- levels(ozone$Month)
  y <- ozone$Ozone

  ### Contrasts for the independent RStan design ----
  helmert <- contr.helmert(months)
  orthonormal <- sweep(helmert, 2, sqrt(colSums(helmert^2)), "/")
  meandif <- orthonormal / sqrt(1 - 1 / length(months))

  ### The same regression in BayesTools and RStan ----
  likelihood <- "model {
    for (i in 1:N) { y[i] ~ dnorm(mu[i], 1) }
  }"
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))

  fit_pair <- function(name, shrinkage, file, data) {
    jags <- scenario_fit(paste0("jags_", name), {
      JAGS_fit(likelihood, data = list(y = y, N = length(y)), prior_list = list(),
        formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
        formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)), sh = shrinkage)),
        chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    }, cache_version = regularized_fit_version(likelihood, ozone, y, shrinkage))
    stan <- scenario_fit(paste0("stan_", name), {
      rstan::sampling(scenario_stan_model(file), data = data,
        chains = 3, cores = 3, iter = 6500, warmup = 1500, seed = 6283, refresh = 0,
        control = list(adapt_delta = 0.995, max_treedepth = 14))
    }, cache_version = regularized_fit_version(file, scenario_stan_version(file), data))
    list(name = name, jags = jags, stan = stan, stan_data = data)
  }

  ### Compare actual month deviations and the wind slope ----
  effects <- c(months, "Wind")
  jags_columns <- setNames(c(paste0("(mu) Month[dif: ", months, "]"), "(mu) Wind"), effects)
  effect_selectors <- setNames(c(paste0("mu_Month[", months, "]"), "mu_Wind"), effects)
  show_ozone <- function(pair, basis, title) {
    jags <- JAGS_estimates_table(pair$jags, transform_factors = TRUE, return_samples = TRUE)
    beta <- rstan::extract(pair$stan, pars = "beta")$beta
    # The first four Stan coefficients are contrasts. Multiply by the
    # contrast row for each month to put them on the BayesTools level scale.
    month_effects <- beta[, 1:4] %*% t(basis)
    colnames(month_effects) <- months
    draws <- list(BayesTools = scenario_draws(jags, jags_columns),
                  RStan = cbind(month_effects, Wind = beta[, 5]))
    scenario_comparison(paste0("comparison-", pair$name), effects, posteriors = draws,
                        main = title, show_table = TRUE)
    scenario_plot(paste0("distributions-", pair$name),
                  scenario_shrinkage_plot(pair$jags, draws, effect_selectors, title))
  }

  ### Fit, plot, and bridge one prior under one coding ----
  bases <- list(meandif = meandif, orthonormal = orthonormal)
  coding_labels <- c(meandif = "mean-difference", orthonormal = "orthonormal")
  prior_labels <- c(ridge = "Ridge", lasso = "Group lasso", horseshoe = "Horseshoe",
                    regularized_horseshoe = "Finite-slab horseshoe")
  compare_prior <- function(contrast, family, stan_file, slab_scale = NULL) {
    basis <- bases[[contrast]]
    X <- cbind(basis[as.integer(ozone$Month), ], Wind = ozone$Wind)
    data <- list(N = nrow(X), K = ncol(X), X = X, y = y)
    if (family == "horseshoe") {
      # One local scale for the month contrasts, another for wind.
      data <- c(data, list(J = 2L, unit = c(1L, 1L, 1L, 1L, 2L)))
    }
    # A finite slab uses the default four degrees of freedom in both engines.
    shrinkage <- prior_regularized(family, global, contrast = contrast, slab_scale = slab_scale)
    name <- if (is.null(slab_scale)) family else "regularized_horseshoe"
    pair <- fit_pair(paste(contrast, name, sep = "_"), shrinkage, stan_file, data)
    show_ozone(pair, basis, paste(prior_labels[[name]], "with", coding_labels[[contrast]], "month contrasts"))
    regularized_bridge(pair$name, pair$jags, pair$stan, pair$stan_data)
  }

  ### Mean-difference contrasts ----
  md <- list()
  md$Ridge     <- compare_prior("meandif", "ridge",     "ridge.stan")
  md$Lasso     <- compare_prior("meandif", "lasso",     "grouped-lasso.stan")
  md$Horseshoe <- compare_prior("meandif", "horseshoe", "horseshoe.stan")
  md[["Finite-slab horseshoe"]] <- compare_prior("meandif", "horseshoe", "finite-slab-horseshoe.stan", slab_scale = 2)
  scenario_text("model-comparison-meandif", regularized_evidence_table(md))

  ### Orthonormal contrasts ----
  ortho <- list()
  ortho$Ridge     <- compare_prior("orthonormal", "ridge",     "ridge.stan")
  ortho$Lasso     <- compare_prior("orthonormal", "lasso",     "grouped-lasso.stan")
  ortho$Horseshoe <- compare_prior("orthonormal", "horseshoe", "horseshoe.stan")
  ortho[["Finite-slab horseshoe"]] <- compare_prior("orthonormal", "horseshoe", "finite-slab-horseshoe.stan", slab_scale = 2)
  scenario_text("model-comparison-orthonormal", regularized_evidence_table(ortho))
})
