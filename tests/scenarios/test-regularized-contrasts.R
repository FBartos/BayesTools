if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized-contrasts")
# testthat::test_file("tests/scenarios/test-regularized-contrasts.R")

### Description
# Do the daily New York ozone readings differ by month after accounting for
# wind? Compare May through September at average wind. Each month effect is
# a deviation from the unweighted mean of the five months. Both ozone and
# wind are standardized, and the residual SD is fixed at one.
#
# BayesTools and RStan use the same prior on the coded coefficients. With
# either contrast the four month coefficients shrink together. Wind has its
# own shrinkage unit and shares the global scale with Month.

testthat::test_that("Monthly ozone effects under centered shrinkage priors", {

  ozone <- na.omit(airquality[c("Ozone", "Month", "Wind")])
  ozone$Ozone <- as.numeric(scale(ozone$Ozone))
  ozone$Wind  <- as.numeric(scale(ozone$Wind))
  ozone$Month <- factor(ozone$Month, labels = month.abb[5:9])
  months <- levels(ozone$Month)

  ### BayesTools regression ----
  likelihood <- "model {
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], 1)
    }
  }"
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))
  fit_ozone <- function(shrinkage) {
    JAGS_fit(
      likelihood, data = list(y = ozone$Ozone, N = nrow(ozone)), prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)), sh = shrinkage)),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429
    )
  }

  ### RStan regression ----
  # Build the same two design matrices independently from base R Helmert
  # contrasts. The mean-difference rows have unit length, whereas the
  # orthonormal columns have unit length.
  helmert <- contr.helmert(months)
  orthonormal <- sweep(helmert, 2, sqrt(colSums(helmert^2)), "/")
  meandif <- orthonormal / sqrt(1 - 1 / length(months))

  stan_data <- function(basis) {
    X <- cbind(basis[as.integer(ozone$Month), ], Wind = ozone$Wind)
    list(N = nrow(X), K = ncol(X), X = X, y = ozone$Ozone)
  }
  data_orthonormal <- stan_data(orthonormal)
  data_meandif     <- stan_data(meandif)
  # One horseshoe local scale for all four month contrasts, another for wind.
  local_scales <- list(J = 2L, unit = c(1L, 1L, 1L, 1L, 2L))

  fit_stan <- function(file, data) {
    rstan::sampling(
      scenario_stan_model(file), data = data, chains = 3, cores = 3,
      iter = 6500, warmup = 1500, seed = 6283, refresh = 0,
      control = list(adapt_delta = 0.995, max_treedepth = 14)
    )
  }

  ### Outcomes to compare ----
  # BayesTools can report the five month deviations directly. In RStan,
  # multiply its four month coefficients by the contrast row for each month.
  # The last coefficient is the change for a one-SD increase in wind.
  effects <- c(months, "Wind")
  ozone_posteriors <- function(jags, stan, basis) {
    jags_levels <- JAGS_estimates_table(jags, transform_factors = TRUE, return_samples = TRUE)
    jags_columns <- setNames(c(paste0("(mu) Month[dif: ", months, "]"), "(mu) Wind"), effects)
    beta <- rstan::extract(stan, pars = "beta")$beta
    stan_levels <- beta[, 1:4] %*% t(basis)
    colnames(stan_levels) <- months
    list(BayesTools = scenario_draws(jags_levels, jags_columns),
         RStan = cbind(stan_levels, Wind = beta[, 5]))
  }

  ### Mean-difference contrasts ----

  jags_ridge_meandif <- scenario_fit("jags_meandif_ridge", {
    fit_ozone(prior_regularized("ridge", global, contrast = "meandif"))
  })
  stan_ridge_meandif <- scenario_fit("stan_meandif_ridge", {
    fit_stan("ridge.stan", data_meandif)
  }, cache_version = scenario_stan_version("ridge.stan"))
  scenario_comparison("comparison-meandif_ridge", effects,
    posteriors = ozone_posteriors(jags_ridge_meandif, stan_ridge_meandif, meandif),
    main = "Ridge with mean-difference month contrasts", show_table = TRUE)

  jags_lasso_meandif <- scenario_fit("jags_meandif_lasso", {
    fit_ozone(prior_regularized("lasso", global, contrast = "meandif"))
  })
  stan_lasso_meandif <- scenario_fit("stan_meandif_lasso", {
    fit_stan("grouped-lasso.stan", data_meandif)
  }, cache_version = scenario_stan_version("grouped-lasso.stan"))
  scenario_comparison("comparison-meandif_lasso", effects,
    posteriors = ozone_posteriors(jags_lasso_meandif, stan_lasso_meandif, meandif),
    main = "Group lasso with mean-difference month contrasts", show_table = TRUE)

  jags_horseshoe_meandif <- scenario_fit("jags_meandif_horseshoe", {
    fit_ozone(prior_regularized("horseshoe", global, contrast = "meandif"))
  })
  stan_horseshoe_meandif <- scenario_fit("stan_meandif_horseshoe", {
    fit_stan("horseshoe.stan", c(data_meandif, local_scales))
  }, cache_version = scenario_stan_version("horseshoe.stan"))
  scenario_comparison("comparison-meandif_horseshoe", effects,
    posteriors = ozone_posteriors(jags_horseshoe_meandif, stan_horseshoe_meandif, meandif),
    main = "Horseshoe with mean-difference month contrasts", show_table = TRUE)

  jags_regularized_horseshoe_meandif <- scenario_fit("jags_meandif_regularized_horseshoe", {
    fit_ozone(prior_regularized("horseshoe", global, contrast = "meandif", slab_scale = 2, slab_df = 4))
  })
  stan_regularized_horseshoe_meandif <- scenario_fit("stan_meandif_regularized_horseshoe", {
    fit_stan("finite-slab-horseshoe.stan", c(data_meandif, local_scales))
  }, cache_version = scenario_stan_version("finite-slab-horseshoe.stan"))
  scenario_comparison("comparison-meandif_regularized_horseshoe", effects,
    posteriors = ozone_posteriors(jags_regularized_horseshoe_meandif, stan_regularized_horseshoe_meandif, meandif),
    main = "Finite-slab horseshoe with mean-difference month contrasts", show_table = TRUE)

  ### Orthonormal contrasts ----

  jags_ridge_orthonormal <- scenario_fit("jags_orthonormal_ridge", {
    fit_ozone(prior_regularized("ridge", global, contrast = "orthonormal"))
  })
  stan_ridge_orthonormal <- scenario_fit("stan_orthonormal_ridge", {
    fit_stan("ridge.stan", data_orthonormal)
  }, cache_version = scenario_stan_version("ridge.stan"))
  scenario_comparison("comparison-orthonormal_ridge", effects,
    posteriors = ozone_posteriors(jags_ridge_orthonormal, stan_ridge_orthonormal, orthonormal),
    main = "Ridge with orthonormal month contrasts", show_table = TRUE)

  jags_lasso_orthonormal <- scenario_fit("jags_orthonormal_lasso", {
    fit_ozone(prior_regularized("lasso", global, contrast = "orthonormal"))
  })
  stan_lasso_orthonormal <- scenario_fit("stan_orthonormal_lasso", {
    fit_stan("grouped-lasso.stan", data_orthonormal)
  }, cache_version = scenario_stan_version("grouped-lasso.stan"))
  scenario_comparison("comparison-orthonormal_lasso", effects,
    posteriors = ozone_posteriors(jags_lasso_orthonormal, stan_lasso_orthonormal, orthonormal),
    main = "Group lasso with orthonormal month contrasts", show_table = TRUE)

  jags_horseshoe_orthonormal <- scenario_fit("jags_orthonormal_horseshoe", {
    fit_ozone(prior_regularized("horseshoe", global, contrast = "orthonormal"))
  })
  stan_horseshoe_orthonormal <- scenario_fit("stan_orthonormal_horseshoe", {
    fit_stan("horseshoe.stan", c(data_orthonormal, local_scales))
  }, cache_version = scenario_stan_version("horseshoe.stan"))
  scenario_comparison("comparison-orthonormal_horseshoe", effects,
    posteriors = ozone_posteriors(jags_horseshoe_orthonormal, stan_horseshoe_orthonormal, orthonormal),
    main = "Horseshoe with orthonormal month contrasts", show_table = TRUE)

  jags_regularized_horseshoe_orthonormal <- scenario_fit("jags_orthonormal_regularized_horseshoe", {
    fit_ozone(prior_regularized("horseshoe", global, contrast = "orthonormal", slab_scale = 2, slab_df = 4))
  })
  stan_regularized_horseshoe_orthonormal <- scenario_fit("stan_orthonormal_regularized_horseshoe", {
    fit_stan("finite-slab-horseshoe.stan", c(data_orthonormal, local_scales))
  }, cache_version = scenario_stan_version("finite-slab-horseshoe.stan"))
  scenario_comparison("comparison-orthonormal_regularized_horseshoe", effects,
    posteriors = ozone_posteriors(jags_regularized_horseshoe_orthonormal, stan_regularized_horseshoe_orthonormal, orthonormal),
    main = "Finite-slab horseshoe with orthonormal month contrasts", show_table = TRUE)
})
