if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("sleepstudy")
# testthat::test_file("tests/scenarios/test-sleepstudy.R")

### Description
# random intercepts and slopes and a heterogeneous compound-symmetry covariance
# for the lme4 sleepstudy data (the first two examples of the RandomEffects
# vignette): estimates compared with lme4 and with brms under matched priors on
# the fitted and the original Days scale, Bayes factors for the random slope
# (bridge sampling in BayesTools and brms, and spike and slab), for the
# intercept-slope correlation (bridge sampling in BayesTools and brms), and for
# the HCS correlation equal to 0.9 (bridge sampling, Savage-Dickey, and a
# mixture prior), and fitted, population-level, and new-subject predictions
# compared with lme4 and brms

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

  ### brms fits ----
  # the intercept, slope, and independent models with the BayesTools priors:
  # normal(0, 1) on the intercept of `0 + Intercept` and on the standardized
  # Days_s coefficient, and normal(0, 1) on the SDs and sigma, which brms
  # bounds at zero (half-normal); each block returns the fit and its
  # bridge-sampling marginal likelihood
  brms_intercept <- scenario_fit("brms_intercept", {
    fit <- brms::brm(
      Reaction_z ~ 0 + Intercept + Days_s + (1 | Subject), data = sleepstudy,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sd),
                brms::prior(normal(0, 1), class = sigma)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_slope <- scenario_fit("brms_slope", {
    fit <- brms::brm(
      Reaction_z ~ 0 + Intercept + Days_s + (1 + Days_s | Subject), data = sleepstudy,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sd),
                brms::prior(normal(0, 1), class = sigma), brms::prior(lkj(1), class = cor)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_independent <- scenario_fit("brms_independent", {
    fit <- brms::brm(
      Reaction_z ~ 0 + Intercept + Days_s + (1 + Days_s || Subject), data = sleepstudy,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sd),
                brms::prior(normal(0, 1), class = sigma)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })

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

  ### lme4 and brms comparisons ----
  # one row per parameter and source: the lme4 estimate, then the BayesTools
  # and brms posterior summaries; `bayes_names` maps the parameters to
  # BayesTools table columns
  compare_with_lme4 <- function(lme4_values, fit, bayes_names, transform_scaled = FALSE, brms_draws = NULL) {
    samples <- JAGS_estimates_table(fit, transform_scaled = transform_scaled, return_samples = TRUE)
    scenario_comparison_table(
      names(bayes_names), references = list(lme4 = lme4_values),
      posteriors = c(list(BayesTools = scenario_draws(samples, bayes_names)),
                     if (!is.null(brms_draws)) list(brms = brms_draws))
    )
  }
  brms_names <- c(
    intercept = "b_Intercept", Days = "b_Days_s", sigma = "sigma",
    "sd(intercept)" = "sd_Subject__Intercept", "sd(Days)" = "sd_Subject__Days_s",
    "cor(intercept,Days)" = "cor_Subject__Intercept__Days_s"
  )
  brms_draws <- function(fit, parameters) {
    scenario_draws(as.matrix(fit), brms_names[parameters])
  }
  # brms draws on the original Days scale: with Days_s = (Days - m) / s, the
  # coefficients transform linearly and the random intercepts and slopes by
  # A = [1, -m / s; 0, 1 / s], so their covariance becomes A Sigma A'
  brms_original_scale <- function(fit) {
    draws <- brms_draws(fit, names(brms_names))
    m <- mean(sleepstudy$Days)
    s <- stats::sd(sleepstudy$Days)
    var_intercept <- draws[, "sd(intercept)"]^2
    var_days      <- draws[, "sd(Days)"]^2
    covariance    <- draws[, "cor(intercept,Days)"] * draws[, "sd(intercept)"] * draws[, "sd(Days)"]
    sd_intercept  <- sqrt(var_intercept - 2 * m / s * covariance + (m / s)^2 * var_days)
    sd_days       <- draws[, "sd(Days)"] / s
    cbind(
      intercept = draws[, "intercept"] - draws[, "Days"] * m / s, Days = draws[, "Days"] / s,
      sigma = draws[, "sigma"], "sd(intercept)" = sd_intercept, "sd(Days)" = sd_days,
      "cor(intercept,Days)" = (covariance - m / s * var_days) / s / (sd_intercept * sd_days)
    )
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
    compare_with_lme4(lme4_coefficients(lme4_intercept), fit_intercept, slope_names[1:4],
                      brms_draws = brms_draws(brms_intercept[["fit"]], names(slope_names)[1:4])))
  scenario_text("lme4-comparison-slope",
    compare_with_lme4(lme4_coefficients(lme4_slope), fit_slope, slope_names,
                      brms_draws = brms_draws(brms_slope[["fit"]], names(slope_names))))
  scenario_text("lme4-comparison-slope-original-scale",
    compare_with_lme4(lme4_coefficients(lme4_slope_days, days = "Days"), fit_slope, slope_names, transform_scaled = TRUE,
                      brms_draws = brms_original_scale(brms_slope[["fit"]])))
  scenario_text("lme4-comparison-independent",
    compare_with_lme4(lme4_coefficients(lme4_independent), fit_independent, slope_names[1:5],
                      brms_draws = brms_draws(brms_independent[["fit"]], names(slope_names)[1:5])))

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
  # bridge sampling compares two separately fitted models (in BayesTools and
  # in brms), the other methods obtain the same Bayes factor from a single fit

  # random slope: the spike-and-slab prior on sd(Days) includes or excludes
  # the slope within one fit; with every draw including it, the inclusion BF
  # is only a bound set by the number of draws
  scenario_text("bayes-factor-random-slope", {
    slope_spike <- JAGS_inference_table(fit_slope_spike)[["inclusion_BF"]]
    spike_bound <- identical(attr(slope_spike, "bound_operator"), ">")
    scenario_bf_table(
      H1 = "(1 + Days | Subject)",
      H0 = "(1 | Subject)",
      log_BF10 = c(
        "bridge sampling"        = scenario_log_bf_bridge(slope[["marglik"]], intercept[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_slope[["marglik"]], brms_intercept[["marglik"]]),
        "spike and slab"         = log(as.numeric(slope_spike))
      ),
      note = c("", "", if (spike_bound) "lower bound: all draws include the slope" else "")
    )
  })

  # intercept-slope correlation: LKJ(1) correlation against independent
  # intercepts and slopes
  scenario_text("bayes-factor-correlation", {
    scenario_bf_table(
      H1 = "(1 + Days | Subject), cor ~ LKJ(1)",
      H0 = "(1 + Days || Subject)",
      log_BF10 = c(
        "bridge sampling"        = scenario_log_bf_bridge(slope[["marglik"]], independent[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_slope[["marglik"]], brms_independent[["marglik"]])
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
    scenario_bf_table(
      H1 = "cor ~ uniform(-1/9, 1)",
      H0 = "cor = 0.9",
      log_BF10 = c(
        "bridge sampling" = scenario_log_bf_bridge(hcs[["marglik"]], hcs_null[["marglik"]]),
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

  # brms shares the priors, so the posterior means of the conditional and the
  # population-level (fixed-part) predictions agree within Monte Carlo error
  scenario_plot("brms-fitted-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    fitted_mean <- function(fit, target) rowMeans(JAGS_predict_formula(fit, "mu", formula_target = target)$value)
    brms_mean   <- function(fit, re_formula = NULL) colMeans(brms::posterior_epred(fit, re_formula = re_formula))
    scenario_agreement_plot(brms_mean(brms_intercept[["fit"]]),   fitted_mean(fit_intercept, "conditional"),
                            "Random intercepts", reference_label = "brms")
    scenario_agreement_plot(brms_mean(brms_slope[["fit"]]),       fitted_mean(fit_slope, "conditional"),
                            "Random intercepts and slopes", reference_label = "brms")
    scenario_agreement_plot(brms_mean(brms_independent[["fit"]]), fitted_mean(fit_independent, "conditional"),
                            "Independent intercepts and slopes", reference_label = "brms")
    scenario_agreement_plot(brms_mean(brms_slope[["fit"]], re_formula = NA), fitted_mean(fit_slope, "fixed"),
                            "Population level (slopes)", reference_label = "brms")
  })

  ### predictions for a new subject ----
  # lme4 predicts a new subject at the population level, which the "zero"
  # policy reproduces; the "sample" policy draws the new subject's intercept
  # and slope from the fitted covariance, as brms sample_new_levels = "gaussian"
  scenario_text("predictions-new-subject", {
    new_subject <- data.frame(Subject = factor("new"), Days = 0:9)
    new_subject$Days_s <- (new_subject$Days - mean(sleepstudy$Days)) / stats::sd(sleepstudy$Days)
    predict_new <- function(new_levels, seed = NULL) {
      JAGS_predict_formula(fit_slope, "mu", data = new_subject, formula_target = "conditional",
                           new_levels = new_levels, seed = seed)$value
    }
    scenario_prediction_table(new_subject[, "Days", drop = FALSE], list(
      "lme4"              = stats::predict(lme4_slope, newdata = new_subject, allow.new.levels = TRUE),
      "BayesTools zero"   = predict_new("zero"),
      "BayesTools sample" = predict_new("sample", seed = 1),
      "brms sample"       = t(brms::posterior_epred(brms_slope[["fit"]], newdata = new_subject,
                                                    allow_new_levels = TRUE, sample_new_levels = "gaussian"))
    ))
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
