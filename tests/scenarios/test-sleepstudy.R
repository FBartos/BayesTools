if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("sleepstudy")
# testthat::test_file("tests/scenarios/test-sleepstudy.R")

### Description
# random intercepts and slopes and a heterogeneous compound-symmetry covariance
# for the lme4 sleepstudy data (the first two examples of the RandomEffects
# vignette): estimates compared with lme4 on the fitted and the original Days
# scale, Bayes factors for the random slope (bridge sampling and spike and
# slab), for the intercept-slope correlation (bridge sampling), and for the
# HCS correlation equal to 0.9 (bridge sampling, Savage-Dickey, and a
# mixture prior)

testthat::test_that("sleepstudy random intercepts, slopes, and HCS", {

  data("sleepstudy", package = "lme4")
  sleepstudy$Subject    <- factor(sleepstudy$Subject)
  sleepstudy$Reaction_z <- as.numeric(scale(sleepstudy$Reaction))
  sleepstudy$Days_s     <- as.numeric(scale(sleepstudy$Days))
  sleepstudy$Days_f     <- factor(sleepstudy$Days)

  likelihood <- "model{
    for(i in 1:N){
      Reaction_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  data_list     <- list(Reaction_z = sleepstudy$Reaction_z, N = nrow(sleepstudy))
  unit_prior    <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["Reaction_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  ### lme4 fits ----
  lme4_intercept   <- lme4::lmer(Reaction_z ~ Days_s + (1 | Subject),           data = sleepstudy, REML = FALSE)
  lme4_slope       <- lme4::lmer(Reaction_z ~ Days_s + (1 + Days_s | Subject),  data = sleepstudy, REML = FALSE)
  lme4_slope_days  <- lme4::lmer(Reaction_z ~ Days + (1 + Days | Subject),      data = sleepstudy, REML = FALSE)
  lme4_independent <- lme4::lmer(Reaction_z ~ Days_s + (1 + Days_s || Subject), data = sleepstudy, REML = FALSE)
  lme4_hcs         <- lme4::lmer(Reaction_z ~ Days_s + cs(0 + Days_f | Subject), data = sleepstudy, REML = FALSE)

  ### BayesTools fits ----
  # each fit block returns the fit and, where a Bayes factor needs it, its
  # bridge-sampling marginal likelihood
  intercept <- scenario_fit("fit_intercept", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + (1 | Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(sd = unit_sd_prior, parameterization = "mean_centered"))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 1000, sample = 4000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(
      fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1))
  })
  slope <- scenario_fit("fit_slope", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + (1 + Days | Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(sd = unit_sd_prior, cor = prior_lkj(eta = 1)))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 1000, sample = 4000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(
      fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1))
  })
  independent <- scenario_fit("fit_independent", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + (1 + Days || Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(sd = unit_sd_prior))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 1000, sample = 4000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(
      fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1))
  })
  fit_slope_spike <- scenario_fit("fit_slope_spike", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + (1 + Days | Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(
        cor   = prior_lkj(eta = 1),
        terms = list(intercept = unit_sd_prior, Days = prior_spike_and_slab(prior_parameter = unit_sd_prior))
      ))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
  })
  hcs <- scenario_fit("fit_hcs", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + hcs(Days_f | Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(sd = unit_sd_prior))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 2000, sample = 5000, seed = 31
    )
    list(fit = fit, marglik = JAGS_bridgesampling(
      fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1))
  })
  # the null model of the HCS correlation test: the correlation fixed at 0.9
  hcs_null <- scenario_fit("fit_hcs_null", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + hcs(Days_f | Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(
        sd         = unit_sd_prior,
        covariance = random_covariance(cor = prior("spike", list(location = 0.9)), cor_scale = "cor")
      ))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 2000, sample = 5000, seed = 31
    )
    list(fit = fit, marglik = JAGS_bridgesampling(
      fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1))
  })
  fit_hcs_test <- scenario_fit("fit_hcs_test", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Days + hcs(Days_f | Subject)), formula_data_list = list(mu = sleepstudy),
      formula_prior_list = list(mu = list(intercept = unit_prior, Days = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(Subject = random_block(
        sd         = unit_sd_prior,
        covariance = random_covariance(
          cor       = prior_mixture(
            prior_list = list(
              prior("spike",   list(location = 0.9)),
              prior("uniform", list(a = -1/9, b = 1))
            ),
            is_null = c(TRUE, FALSE)
          ),
          cor_scale = "cor"
        )
      ))),
      formula_scale_list = list(mu = list(Days = TRUE)),
      chains = 4, adapt = 1000, burnin = 2000, sample = 5000, seed = 31
    )
  })
  fit_intercept   <- intercept[["fit"]]
  fit_slope       <- slope[["fit"]]
  fit_independent <- independent[["fit"]]
  fit_hcs         <- hcs[["fit"]]

  ### model summaries ----
  scenario_text("estimates-intercept",   JAGS_estimates_table(fit_intercept))
  scenario_text("estimates-slope",       JAGS_estimates_table(fit_slope))
  scenario_text("estimates-independent", JAGS_estimates_table(fit_independent))
  scenario_text("estimates-slope-spike", JAGS_estimates_table(fit_slope_spike))
  scenario_text("estimates-hcs",         JAGS_estimates_table(fit_hcs))
  scenario_text("estimates-hcs-test",    JAGS_estimates_table(fit_hcs_test))

  ### lme4 comparisons ----
  posterior_summary <- function(draws) {
    c(estimate = mean(draws), median = stats::median(draws),
      lCI = unname(stats::quantile(draws, 0.025)), uCI = unname(stats::quantile(draws, 0.975)))
  }
  # one row per parameter and source: the lme4 estimate, then the BayesTools
  # posterior summaries; `bayes_names` maps the parameters to table columns
  compare_with_lme4 <- function(lme4_values, fit, bayes_names, transform_scaled = FALSE) {
    samples <- JAGS_estimates_table(fit, transform_scaled = transform_scaled, return_samples = TRUE)
    bayes   <- t(apply(samples[, bayes_names, drop = FALSE], 2, posterior_summary))
    out <- rbind(
      data.frame(source = "lme4", parameter = names(bayes_names), estimate = unname(lme4_values[names(bayes_names)]),
                 median = NA_real_, lCI = NA_real_, uCI = NA_real_),
      data.frame(source = "BayesTools", parameter = names(bayes_names), bayes, row.names = NULL)
    )
    out[, 3:6] <- round(out[, 3:6], 3)
    out <- out[order(match(out$parameter, names(bayes_names))), ]
    rownames(out) <- NULL
    out
  }
  lme4_coefficients <- function(model, days = "Days_s") {
    vc  <- lme4::VarCorr(model)
    sds <- unlist(lapply(vc, attr, "stddev"))
    cor <- attr(vc[["Subject"]], "correlation")
    c(
      intercept       = lme4::fixef(model)[["(Intercept)"]],
      Days            = lme4::fixef(model)[[days]],
      sigma           = stats::sigma(model),
      "sd(intercept)" = unname(sds[grepl("Intercept", names(sds))]),
      "sd(Days)"      = if (any(grepl("Days", names(sds)))) unname(sds[grepl("Days", names(sds))]) else NA_real_,
      "cor(intercept,Days)" = if (ncol(cor) == 2L) cor[1L, 2L] else NA_real_
    )
  }
  slope_names <- c(
    intercept = "(mu) intercept", Days = "(mu) Days", sigma = "sigma",
    "sd(intercept)" = "(mu) sd(intercept)", "sd(Days)" = "(mu) sd(Days)",
    "cor(intercept,Days)" = "(mu) cor(intercept,Days)"
  )

  scenario_text("lme4-comparison-intercept",
    compare_with_lme4(lme4_coefficients(lme4_intercept), fit_intercept, slope_names[1:4]))
  scenario_text("lme4-comparison-slope",
    compare_with_lme4(lme4_coefficients(lme4_slope), fit_slope, slope_names))
  scenario_text("lme4-comparison-slope-original-scale",
    compare_with_lme4(lme4_coefficients(lme4_slope_days, days = "Days"), fit_slope, slope_names, transform_scaled = TRUE))
  scenario_text("lme4-comparison-independent",
    compare_with_lme4(lme4_coefficients(lme4_independent), fit_independent, slope_names[1:5]))

  scenario_text("lme4-comparison-hcs", {
    lme4_vc  <- as.matrix(lme4::VarCorr(lme4_hcs)[["Subject"]])
    day_sds  <- stats::setNames(sqrt(diag(lme4_vc)), paste0("sd(Days_f[", levels(sleepstudy$Days_f), "])"))
    lme4_values <- c(
      intercept = lme4::fixef(lme4_hcs)[["(Intercept)"]], Days = lme4::fixef(lme4_hcs)[["Days_s"]],
      sigma = stats::sigma(lme4_hcs), day_sds, cor = stats::cov2cor(lme4_vc)[2L, 1L]
    )
    hcs_names <- c(
      intercept = "(mu) intercept", Days = "(mu) Days", sigma = "sigma",
      stats::setNames(paste0("(mu) ", names(day_sds)), names(day_sds)), cor = "(mu) cor"
    )
    compare_with_lme4(lme4_values, fit_hcs, hcs_names)
  })

  ### Bayes factors ----
  # each test compares H1 against H0 as log BF10 (positive values favour H1):
  # bridge sampling compares two separately fitted models, the other methods
  # obtain the same Bayes factor from a single fit
  log_bf10_bridge <- function(marglik_1, marglik_0) {
    bridgesampling::bf(marglik_1, marglik_0, log = TRUE)[["bf"]]
  }
  bf_comparison <- function(H1, H0, log_BF10, note = "") {
    data.frame(H1 = H1, H0 = H0, method = names(log_BF10),
               log_BF10 = round(unname(log_BF10), 3), note = note)
  }

  # random slope: the spike-and-slab prior on sd(Days) includes or excludes
  # the slope within one fit; with every draw including it, the inclusion BF
  # is only a bound set by the number of draws
  scenario_text("bayes-factor-random-slope", {
    slope_spike <- JAGS_inference_table(fit_slope_spike)[["inclusion_BF"]]
    spike_bound <- identical(attr(slope_spike, "bound_operator"), ">")
    bf_comparison(
      H1 = "(1 + Days | Subject)",
      H0 = "(1 | Subject)",
      log_BF10 = c(
        "bridge sampling" = log_bf10_bridge(slope[["marglik"]], intercept[["marglik"]]),
        "spike and slab"  = log(as.numeric(slope_spike))
      ),
      note = c("", if (spike_bound) "lower bound: all draws include the slope" else "")
    )
  })

  # intercept-slope correlation: LKJ(1) correlation against independent
  # intercepts and slopes
  scenario_text("bayes-factor-correlation", {
    bf_comparison(
      H1 = "(1 + Days | Subject), cor ~ LKJ(1)",
      H0 = "(1 + Days || Subject)",
      log_BF10 = c(
        "bridge sampling" = log_bf10_bridge(slope[["marglik"]], independent[["marglik"]])
      )
    )
  })

  # HCS correlation: the default uniform(-1/9, 1) prior against cor = 0.9;
  # Savage-Dickey compares the prior and posterior densities at 0.9 of the
  # free fit, the mixture prior samples both hypotheses in one fit
  scenario_text("bayes-factor-hcs-correlation", {
    hcs_cor <- parameter_mixed_posterior(
      fit_hcs, parameter_catalog_resolve(parameter_catalog(fit_hcs), "(mu) cor")
    )
    savage_dickey <- hypothesis_BF(hcs_cor, hypothesis = "`(mu) cor` = 0.9", BF01 = TRUE, seed = 1)
    mixture <- JAGS_inference_table(fit_hcs_test)[["inclusion_BF"]]
    bf_comparison(
      H1 = "cor ~ uniform(-1/9, 1)",
      H0 = "cor = 0.9",
      log_BF10 = c(
        "bridge sampling" = log_bf10_bridge(hcs[["marglik"]], hcs_null[["marglik"]]),
        "Savage-Dickey"   = -log(as.numeric(savage_dickey[["BF"]])),
        "mixture prior"   = log(as.numeric(mixture))
      )
    )
  })
  scenario_text("inference-slope-spike", JAGS_inference_table(fit_slope_spike))
  scenario_text("inference-hcs-test",    JAGS_inference_table(fit_hcs_test))

  ### random effects and fitted values ----
  # the random-effect contribution of each observation: the conditional
  # prediction minus its fixed part (lme4: fitted minus population-level).
  # Each subject has one observation per day, so its HCS day effects separate
  # from the residual only through the covariance structure; those panels
  # scatter more around lme4 than the intercept and slope models
  random_contribution <- function(fit) {
    prediction <- JAGS_predict_formula(fit, "mu", formula_target = "conditional", components = TRUE)
    rowMeans(prediction$value - prediction$components$fixed)
  }
  lme4_contribution <- function(model) {
    stats::fitted(model) - stats::predict(model, re.form = NA)
  }

  scenario_plot("ranef-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(lme4_contribution(lme4_intercept),   random_contribution(fit_intercept),   "Random intercepts")
    scenario_agreement_plot(lme4_contribution(lme4_slope),       random_contribution(fit_slope),       "Random intercepts and slopes")
    scenario_agreement_plot(lme4_contribution(lme4_independent), random_contribution(fit_independent), "Independent intercepts and slopes")
    scenario_agreement_plot(lme4_contribution(lme4_hcs),         random_contribution(fit_hcs),         "HCS day effects")
  })

  scenario_plot("fitted-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    fitted_mean <- function(fit) rowMeans(JAGS_predict_formula(fit, "mu", formula_target = "conditional")$value)
    scenario_agreement_plot(stats::fitted(lme4_intercept),   fitted_mean(fit_intercept),   "Random intercepts")
    scenario_agreement_plot(stats::fitted(lme4_slope),       fitted_mean(fit_slope),       "Random intercepts and slopes")
    scenario_agreement_plot(stats::fitted(lme4_independent), fitted_mean(fit_independent), "Independent intercepts and slopes")
    scenario_agreement_plot(stats::fitted(lme4_hcs),         fitted_mean(fit_hcs),         "HCS")
  })

  ### prior and posterior distributions ----
  posterior_of <- function(fit, name) {
    selection <- parameter_catalog_resolve(parameter_catalog(fit), name)
    stats::setNames(list(parameter_mixed_posterior(fit, selection)), name)
  }
  plot_quantity <- function(fit, name, prior = TRUE, main = name) {
    plot_posterior(posterior_of(fit, name), name, prior = prior, main = main)
  }

  # the random-effect summaries are on the original Days scale; there the
  # intercept SD and the correlation combine the fitted intercept and slope
  # (Days is not centered at 0) and have no deterministic prior density
  scenario_plot("posterior-slope", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    plot_quantity(fit_slope, "(mu) sd(intercept)", prior = FALSE)
    plot_quantity(fit_slope, "(mu) cor(intercept,Days)", prior = FALSE)
    plot_quantity(fit_slope, "(mu) sd(Days)")
    plot_quantity(fit_slope_spike, "(mu) sd(Days)", main = "(mu) sd(Days), spike and slab")
  })

  scenario_plot("posterior-hcs", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    plot_quantity(fit_hcs, "(mu) sd(Days_f[0])")
    plot_quantity(fit_hcs, "(mu) sd(Days_f[9])")
    plot_quantity(fit_hcs, "(mu) cor")
    plot_quantity(fit_hcs_test, "(mu) cor", main = "(mu) cor, spike at 0.9")
  })
})
