if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("cbpp")
# testthat::test_file("tests/scenarios/test-cbpp.R")

### Description
# a binomial GLMM for the lme4 cbpp data (new cases of contagious bovine
# pleuropneumonia among the animals of 15 herds in four periods): random herd
# intercepts with a period effect (treatment contrasts), compared with glmer
# (Laplace and adaptive Gauss-Hermite quadrature) and with brms under matched
# priors; Bayes factors for the period effect, the herd SD (bridge sampling in
# BayesTools and brms, spike and slab, and Savage-Dickey), and
# observation-level overdispersion (bridge sampling in BayesTools and brms);
# fitted, herd-effect, and new-herd predictions on the probability scale

testthat::test_that("cbpp binomial random herd intercepts", {

  data("cbpp", package = "lme4")
  cbpp$obs <- factor(seq_len(nrow(cbpp)))

  likelihood <- "model{
    for(i in 1:N){
      incidence[i] ~ dbin(ilogit(mu[i]), size[i])
    }
  }
  "
  data_list <- list(incidence = cbpp$incidence, size = cbpp$size, N = nrow(cbpp))
  log_posterior <- function(parameters, data) {
    sum(stats::dbinom(data[["incidence"]], data[["size"]], stats::plogis(parameters[["mu"]]), log = TRUE))
  }

  ### glmer fits ----
  glmer_herd       <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  glmer_herd_agq   <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial, nAGQ = 25)
  glmer_overdispersion <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd) + (1 | obs), data = cbpp, family = binomial)

  ### BayesTools fits ----
  # priors on the logit scale: intercept ~ N(0, 2.5), one N(0, 1) per
  # treatment-contrast period coefficient, and half-normal(0, 1) SDs; each
  # block used in a bridge-sampling Bayes factor also returns its marginal
  # likelihood
  herd <- scenario_fit("fit_herd", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ period + (1 | herd)), formula_data_list = list(mu = cbpp),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 2.5)),
        period    = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
      )),
      formula_random_prior_list = list(mu = prior_random(
        herd = random_block(sd = prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf)))
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  no_period <- scenario_fit("fit_no_period", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ 1 + (1 | herd)), formula_data_list = list(mu = cbpp),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(mean = 0, sd = 2.5)))),
      formula_random_prior_list = list(mu = prior_random(
        herd = random_block(sd = prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf)))
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  no_herd <- scenario_fit("fit_no_herd", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ period), formula_data_list = list(mu = cbpp),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 2.5)),
        period    = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  fit_herd_spike <- scenario_fit("fit_herd_spike", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ period + (1 | herd)), formula_data_list = list(mu = cbpp),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 2.5)),
        period    = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
      )),
      formula_random_prior_list = list(mu = prior_random(herd = random_block(
        sd = prior_spike_and_slab(prior_parameter = prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf)))
      ))),
      chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
    )
  })
  overdispersion <- scenario_fit("fit_overdispersion", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ period + (1 | herd) + (1 | obs)), formula_data_list = list(mu = cbpp),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 2.5)),
        period    = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
      )),
      formula_random_prior_list = list(mu = prior_random(
        sd = prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  fit_herd           <- herd[["fit"]]
  fit_overdispersion <- overdispersion[["fit"]]

  ### brms fits ----
  # the same priors: `0 + Intercept` puts N(0, 2.5) on the intercept, the
  # period coefficients keep treatment contrasts with N(0, 1) priors, and brms
  # bounds the SD priors at zero; each block returns the fit and its
  # bridge-sampling marginal likelihood
  brms_herd <- scenario_fit("brms_herd", {
    fit <- brms::brm(
      incidence | trials(size) ~ 0 + Intercept + period + (1 | herd), data = cbpp, family = binomial(),
      prior = c(brms::prior(normal(0, 2.5), class = b, coef = Intercept), brms::prior(normal(0, 1), class = b),
                brms::prior(normal(0, 1), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_no_period <- scenario_fit("brms_no_period", {
    fit <- brms::brm(
      incidence | trials(size) ~ 0 + Intercept + (1 | herd), data = cbpp, family = binomial(),
      prior = c(brms::prior(normal(0, 2.5), class = b, coef = Intercept), brms::prior(normal(0, 1), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_no_herd <- scenario_fit("brms_no_herd", {
    fit <- brms::brm(
      incidence | trials(size) ~ 0 + Intercept + period, data = cbpp, family = binomial(),
      prior = c(brms::prior(normal(0, 2.5), class = b, coef = Intercept), brms::prior(normal(0, 1), class = b)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_overdispersion <- scenario_fit("brms_overdispersion", {
    fit <- brms::brm(
      incidence | trials(size) ~ 0 + Intercept + period + (1 | herd) + (1 | obs), data = cbpp, family = binomial(),
      prior = c(brms::prior(normal(0, 2.5), class = b, coef = Intercept), brms::prior(normal(0, 1), class = b),
                brms::prior(normal(0, 1), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99)
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })

  ### model summaries ----
  scenario_text("estimates-herd",           JAGS_estimates_table(fit_herd))
  scenario_text("estimates-overdispersion", JAGS_estimates_table(fit_overdispersion))

  ### glmer and brms comparisons ----
  # glmer reports ML estimates without the priors' shrinkage, so the period
  # effects and the herd SD differ from the posterior means; brms shares the
  # priors and agrees with BayesTools within Monte Carlo error
  glmer_values <- function(model) {
    sds <- vapply(lme4::VarCorr(model), attr, numeric(1), "stddev")
    c(stats::setNames(lme4::fixef(model), c("intercept", "period[2]", "period[3]", "period[4]")),
      "sd(herd)" = sds[["herd"]], "sd(obs)" = if ("obs" %in% names(sds)) sds[["obs"]] else NA_real_)
  }
  bayes_names <- c(
    intercept = "(mu) intercept", "period[2]" = "(mu) period[2]", "period[3]" = "(mu) period[3]",
    "period[4]" = "(mu) period[4]"
  )
  brms_names <- c(
    intercept = "b_Intercept", "period[2]" = "b_period2", "period[3]" = "b_period3", "period[4]" = "b_period4",
    "sd(herd)" = "sd_herd__Intercept"
  )

  scenario_text("comparison-herd", {
    scenario_comparison_table(
      c(names(bayes_names), "sd(herd)"),
      references = list("glmer Laplace" = glmer_values(glmer_herd), "glmer AGQ" = glmer_values(glmer_herd_agq)),
      posteriors = list(
        BayesTools = scenario_draws(JAGS_estimates_table(fit_herd, return_samples = TRUE),
                                    c(bayes_names, "sd(herd)" = "(mu) sd(intercept)")),
        brms       = scenario_draws(as.matrix(brms_herd[["fit"]]), brms_names)
      )
    )
  })
  scenario_text("comparison-overdispersion", {
    scenario_comparison_table(
      c(names(bayes_names), "sd(herd)", "sd(obs)"),
      references = list("glmer Laplace" = glmer_values(glmer_overdispersion)),
      posteriors = list(
        BayesTools = scenario_draws(JAGS_estimates_table(fit_overdispersion, return_samples = TRUE),
                                    c(bayes_names, "sd(herd)" = "(mu) herd: sd(intercept)", "sd(obs)" = "(mu) obs: sd(intercept)")),
        brms       = scenario_draws(as.matrix(brms_overdispersion[["fit"]]), c(brms_names, "sd(obs)" = "sd_obs__Intercept"))
      )
    )
  })

  ### Bayes factors ----
  # log BF10 (positive values favour H1); bridge sampling compares separately
  # fitted models in BayesTools and in brms, the spike-and-slab prior and
  # Savage-Dickey obtain the herd test from a single fit
  scenario_text("bayes-factor-period", {
    scenario_bf_table(
      H1 = "period + (1 | herd)",
      H0 = "1 + (1 | herd)",
      log_BF10 = c(
        "bridge sampling"        = scenario_log_bf_bridge(herd[["marglik"]], no_period[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_herd[["marglik"]], brms_no_period[["marglik"]])
      )
    )
  })
  scenario_text("bayes-factor-herd", {
    scenario_bf_table(
      H1 = "period + (1 | herd)",
      H0 = "period",
      log_BF10 = c(
        "bridge sampling"        = scenario_log_bf_bridge(herd[["marglik"]], no_herd[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_herd[["marglik"]], brms_no_herd[["marglik"]]),
        "spike and slab"         = log(as.numeric(JAGS_inference_table(fit_herd_spike)[["inclusion_BF"]])),
        "Savage-Dickey"          = log(as.numeric(hypothesis_BF(fit_herd, "sd(intercept) = 0", seed = 1)[["BF"]]))
      )
    )
  })
  scenario_text("bayes-factor-overdispersion", {
    scenario_bf_table(
      H1 = "period + (1 | herd) + (1 | obs)",
      H0 = "period + (1 | herd)",
      log_BF10 = c(
        "bridge sampling"        = scenario_log_bf_bridge(overdispersion[["marglik"]], herd[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_overdispersion[["marglik"]], brms_herd[["marglik"]])
      )
    )
  })
  scenario_text("inference-herd-spike", JAGS_inference_table(fit_herd_spike))

  ### predictions for a new herd ----
  # at a herd effect of 0 (glmer re.form = NA, the BayesTools "zero" policy)
  # and averaged over a new herd's effect (the "sample" policy and brms
  # sample_new_levels = "gaussian"), which pulls the probabilities toward 0.5
  scenario_text("predictions-new-herd", {
    new_herd <- data.frame(herd = factor("new"), period = factor(1:4, levels = levels(cbpp$period)), size = 1)
    predict_new <- function(new_levels, seed = NULL) {
      stats::plogis(JAGS_predict_formula(fit_herd, "mu", data = new_herd, formula_target = "conditional",
                                         new_levels = new_levels, seed = seed)$value)
    }
    scenario_prediction_table(new_herd[, "period", drop = FALSE], list(
      "glmer"             = stats::predict(glmer_herd, newdata = new_herd, re.form = NA, type = "response"),
      "BayesTools zero"   = predict_new("zero"),
      "BayesTools sample" = predict_new("sample", seed = 1),
      "brms sample"       = t(brms::posterior_epred(brms_herd[["fit"]], newdata = new_herd,
                                                    allow_new_levels = TRUE, sample_new_levels = "gaussian"))
    ))
  })

  ### fitted probabilities and herd effects ----
  # fitted probabilities are posterior means of the conditional probability;
  # a herd effect is the conditional linear predictor minus its fixed part
  fitted_probability <- rowMeans(stats::plogis(JAGS_predict_formula(fit_herd, "mu", formula_target = "conditional")$value))
  herd_prediction    <- JAGS_predict_formula(fit_herd, "mu", formula_target = "conditional", components = TRUE)
  herd_effects       <- tapply(rowMeans(herd_prediction$value - herd_prediction$components$fixed), cbpp$herd, mean)
  scenario_plot("fitted-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    brms_probability <- colMeans(brms::posterior_epred(brms_herd[["fit"]])) / cbpp$size
    scenario_agreement_plot(brms_probability, fitted_probability, "Fitted probabilities", reference_label = "brms")
    scenario_agreement_plot(stats::fitted(glmer_herd), fitted_probability, "Fitted probabilities", reference_label = "glmer")
    scenario_agreement_plot(brms::ranef(brms_herd[["fit"]])[["herd"]][levels(cbpp$herd), "Estimate", "Intercept"],
                            unname(herd_effects[levels(cbpp$herd)]), "Herd effects", reference_label = "brms")
    scenario_agreement_plot(lme4::ranef(glmer_herd)[["herd"]][levels(cbpp$herd), "(Intercept)"],
                            unname(herd_effects[levels(cbpp$herd)]), "Herd effects", reference_label = "glmer")
  })

  ### prior and posterior distributions ----
  scenario_plot("posterior-sd", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    plot_quantity <- function(fit, name, main) {
      selection <- parameter_catalog_resolve(parameter_catalog(fit), name)
      plot_posterior(stats::setNames(list(parameter_mixed_posterior(fit, selection)), name), name, prior = TRUE, main = main)
    }
    plot_quantity(fit_herd_spike, "(mu) sd(intercept)", "sd(herd), spike and slab")
    plot_quantity(fit_overdispersion, "(mu) obs: sd(intercept)", "sd(obs)")
  })
})
