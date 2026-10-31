if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("pastes")
# testthat::test_file("tests/scenarios/test-pastes.R")

### Description
# nested random intercepts for batches and casks within batches of the lme4
# Pastes data (the nested example of the RandomEffects vignette): a variance
# allocation prior and independent half-normal SD priors, compared with the
# lme4 maximum-likelihood fit, including the batch and cask effects and the
# fitted values; the batch block is mean-centered

testthat::test_that("Pastes nested random intercepts", {

  data("Pastes", package = "lme4")
  Pastes$strength_z <- as.numeric(scale(Pastes$strength))

  likelihood <- "model{
    for(i in 1:N){
      strength_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  data_list     <- list(strength_z = Pastes$strength_z, N = nrow(Pastes))
  unit_prior    <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  formula_list  <- list(mu = ~ 1 + (1 | batch/cask))

  ### model fits ----
  fit_lme4 <- lme4::lmer(strength_z ~ 1 + (1 | batch/cask), data = Pastes, REML = FALSE)

  fit_allocation <- scenario_fit("fit_allocation", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = formula_list, formula_data_list = list(mu = Pastes),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      # the intercept and the ten batch effects trade off; the mean-centered
      # batch block samples them around the intercept
      formula_random_prior_list = list(mu = prior_random(
        batch      = random_block(parameterization = "mean_centered"),
        allocation = random_variance_allocation(
          name    = "nested_re",
          terms   = c(batch = "batch", cask = "cask_batch"),
          sd      = unit_sd_prior,
          weights = prior("dirichlet", parameters = list(alpha = c(1, 1)))
        )
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 4000, seed = 41
    )
  })
  fit_independent <- scenario_fit("fit_independent", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = formula_list, formula_data_list = list(mu = Pastes),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(
        batch      = random_block(sd = unit_sd_prior, parameterization = "mean_centered"),
        cask_batch = random_block(sd = unit_sd_prior)
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 4000, seed = 42
    )
  })

  ### model summaries ----
  fit_lme4
  scenario_text("estimates-allocation",       JAGS_estimates_table(fit_allocation))
  scenario_text("estimates-allocation-full",  JAGS_estimates_table(fit_allocation, random_effects_summary = "full"))
  scenario_text("estimates-independent",      JAGS_estimates_table(fit_independent))
  scenario_text("estimates-independent-full", JAGS_estimates_table(fit_independent, random_effects_summary = "full"))

  ### lme4 comparison ----
  # the allocation quantities of lme4 and of the independent-prior fit are
  # derived from their two SDs
  samples_allocation  <- JAGS_estimates_table(fit_allocation,  random_effects_summary = "full", return_samples = TRUE)
  samples_independent <- JAGS_estimates_table(fit_independent, random_effects_summary = "full", return_samples = TRUE)
  allocation_summaries <- function(sd_batch, sd_cask) {
    variance_total <- sd_batch^2 + sd_cask^2
    cbind(sd_total = sqrt(variance_total), var_prop_batch = sd_batch^2 / variance_total)
  }

  scenario_text("lme4-comparison", {
    lme4_sd <- vapply(lme4::VarCorr(fit_lme4), attr, numeric(1), "stddev")
    lme4_values <- c(
      intercept      = lme4::fixef(fit_lme4)[["(Intercept)"]],
      sigma          = stats::sigma(fit_lme4),
      sd_batch       = lme4_sd[["batch"]],
      sd_cask        = lme4_sd[["cask:batch"]],
      allocation_summaries(lme4_sd[["batch"]], lme4_sd[["cask:batch"]])[1L, ]
    )
    posterior_summary <- function(draws) {
      c(estimate = mean(draws), median = stats::median(draws),
        lCI = unname(stats::quantile(draws, 0.025)), uCI = unname(stats::quantile(draws, 0.975)))
    }
    allocation_draws <- samples_allocation[, c(
      "(mu) intercept", "sigma", "(mu) batch: sd(intercept)", "(mu) cask_batch: sd(intercept)",
      "(mu) nested_re: sd_total", "(mu) nested_re: var_prop(batch)"
    )]
    independent_draws <- cbind(
      samples_independent[, c("(mu) intercept", "sigma", "(mu) batch: sd(intercept)", "(mu) cask_batch: sd(intercept)")],
      allocation_summaries(samples_independent[, "(mu) batch: sd(intercept)"], samples_independent[, "(mu) cask_batch: sd(intercept)"])
    )
    comparison <- rbind(
      data.frame(fit = "lme4", parameter = names(lme4_values), estimate = unname(lme4_values),
                 median = NA_real_, lCI = NA_real_, uCI = NA_real_),
      data.frame(fit = "allocation", parameter = names(lme4_values),
                 t(apply(allocation_draws, 2, posterior_summary)), row.names = NULL),
      data.frame(fit = "independent", parameter = names(lme4_values),
                 t(apply(independent_draws, 2, posterior_summary)), row.names = NULL)
    )
    comparison[, 3:6] <- round(comparison[, 3:6], 3)
    comparison <- comparison[order(match(comparison$parameter, names(lme4_values))), ]
    rownames(comparison) <- NULL
    comparison
  })

  ### random effects and fitted values ----
  # the effect of a level is its conditional prediction with only its block
  # minus the fixed part; casks are identified by cask and batch (lme4 names
  # them "cask:batch")
  groups <- list(
    batch      = as.character(Pastes$batch),
    cask_batch = paste(Pastes$cask, Pastes$batch, sep = ":")
  )
  level_effects <- function(fit, block) {
    prediction <- JAGS_predict_formula(fit, "mu", formula_target = "conditional", blocks = block, components = TRUE)
    effects    <- rowMeans(prediction$value - prediction$components$fixed)
    tapply(effects, groups[[block]], mean)
  }
  ranef_lme4 <- lme4::ranef(fit_lme4)
  lme4_effects <- list(
    batch      = stats::setNames(ranef_lme4[["batch"]][["(Intercept)"]],      rownames(ranef_lme4[["batch"]])),
    cask_batch = stats::setNames(ranef_lme4[["cask:batch"]][["(Intercept)"]], rownames(ranef_lme4[["cask:batch"]]))
  )

  # the batch effects are less shrunk than lme4's: with ten batches the
  # posterior SD of batch (about 0.45-0.5) exceeds the ML estimate (0.338; REML
  # 0.398), and the allocation prior moves further variance to batch. The
  # differences therefore grow linearly with the lme4 effect without a shift;
  # averaging the closed-form conditional batch effects over the posterior
  # draws reproduces the BayesTools effects within Monte Carlo error
  scenario_plot("ranef-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (block in c("batch", "cask_batch")) {
      reference <- lme4_effects[[block]]
      for (fit_name in c("allocation", "independent")) {
        fit      <- if (fit_name == "allocation") fit_allocation else fit_independent
        estimate <- level_effects(fit, block)
        scenario_agreement_plot(reference, unname(estimate[names(reference)]), paste0(block, " effects (", fit_name, ")"))
      }
    }
  })

  scenario_plot("fitted-comparison", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    reference <- stats::fitted(fit_lme4)
    for (fit_name in c("allocation", "independent")) {
      fit <- if (fit_name == "allocation") fit_allocation else fit_independent
      prediction <- JAGS_predict_formula(fit, "mu", formula_target = "conditional")
      scenario_agreement_plot(reference, rowMeans(prediction$value), paste0("Fitted values (", fit_name, ")"))
    }
  })

  ### prior and posterior distributions ----
  posterior_of <- function(fit, name) {
    selection <- parameter_catalog_resolve(parameter_catalog(fit), name)
    stats::setNames(list(parameter_mixed_posterior(fit, selection)), name)
  }
  plot_quantity <- function(fit, name) {
    plot_posterior(posterior_of(fit, name), name, prior = TRUE, main = name)
  }

  scenario_plot("posterior-allocation", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    plot_quantity(fit_allocation, "(mu) nested_re: sd_total")
    plot_quantity(fit_allocation, "(mu) nested_re: var_prop(batch)")
    plot_quantity(fit_allocation, "(mu) batch: sd(intercept)")
    plot_quantity(fit_allocation, "(mu) cask_batch: sd(intercept)")
  })

  scenario_plot("posterior-independent", {
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    plot_quantity(fit_independent, "(mu) batch: sd(intercept)")
    plot_quantity(fit_independent, "(mu) cask_batch: sd(intercept)")
    plot_quantity(fit_independent, "sigma")
  })
})
