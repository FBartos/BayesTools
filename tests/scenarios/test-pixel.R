if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("pixel")
# testthat::test_file("tests/scenarios/test-pixel.R")

### Description
# autoregressive dog-by-day random effects for the nlme Pixel data (lymph-node
# pixel intensity of 10 dogs on both sides over unequally spaced days 0-21,
# with each dog observed on some of the days): a quadratic time trend with
# correlated random intercepts and slopes per dog, AR(1) day effects over the
# ordered days, and continuous-time AR(1) (CAR) day effects over the actual
# day distances, each with a dog-by-side intercept. Compared with lme4 (random
# slopes and AR(1); lme4 2.0 has no CAR) and with brms for the random slopes
# under matched priors; bridge-sampling marginal likelihoods compare the
# structures, and predictions for dog-days without observations use the
# Gaussian conditional of the unobserved day effects (lme4: their BLUPs).
# HAR is not included: its day-specific SDs are not identified here (the lme4
# fit is singular)

testthat::test_that("Pixel AR(1) and CAR day effects", {

  data("Pixel", package = "nlme")
  Pixel <- as.data.frame(Pixel)
  Pixel$Dog     <- factor(as.character(Pixel$Dog), levels = as.character(1:10))
  Pixel$pixel_z <- as.numeric(scale(Pixel$pixel))
  Pixel$day_s   <- as.numeric(scale(Pixel$day))
  Pixel$day_s2  <- Pixel$day_s^2
  Pixel$day_f   <- factor(Pixel$day)

  likelihood <- "model{
    for(i in 1:N){
      pixel_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  data_list <- list(pixel_z = Pixel$pixel_z, N = nrow(Pixel))
  log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["pixel_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  # priors: N(0, 1) fixed coefficients, half-normal(0, 1) sigma and SDs,
  # LKJ(1), and the default uniform(-1, 1) AR(1) and uniform(0, 1) CAR
  # correlations
  unit_prior    <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  fixed_priors  <- list(mu = list(intercept = unit_prior, day_s = unit_prior, day_s2 = unit_prior))
  marglik <- function(fit) {
    JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1)
  }

  ### lme4 fits ----
  # the structured term comes last: lme4 2.0.1 sorts random terms by their
  # number of levels without moving the structure tags, so an earlier AR(1)
  # term would be fitted as unstructured; the check below guards that
  lme4_slope <- lme4::lmer(pixel_z ~ day_s + day_s2 + (1 + day_s | Dog) + (1 | Dog:Side), data = Pixel, REML = FALSE)
  lme4_ar1   <- lme4::lmer(pixel_z ~ day_s + day_s2 + (1 | Dog:Side) + ar1(0 + day_f | Dog), data = Pixel, REML = FALSE)
  structure_classes <- vapply(lme4::getReCovs(lme4_ar1), function(x) class(x)[[1L]], character(1))
  testthat::expect_identical(unname(structure_classes), c("Covariance.us", "Covariance.ar1"))

  ### BayesTools fits ----
  # long thinned chains: the intercept and the dog effects trade off
  slope <- scenario_fit("fit_slope", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ day_s + day_s2 + (1 + day_s | Dog) + (1 | Dog:Side)), formula_data_list = list(mu = Pixel),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(sd = unit_sd_prior, cor = prior_lkj(eta = 1))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  ar1 <- scenario_fit("fit_ar1", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ day_s + day_s2 + ar1(day_f | Dog) + (1 | Dog:Side)), formula_data_list = list(mu = Pixel),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(sd = unit_sd_prior)),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  car <- scenario_fit("fit_car", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ day_s + day_s2 + car(day | Dog) + (1 | Dog:Side)), formula_data_list = list(mu = Pixel),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(sd = unit_sd_prior)),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  fit_slope <- slope[["fit"]]
  fit_ar1   <- ar1[["fit"]]
  fit_car   <- car[["fit"]]

  ### brms fit ----
  brms_slope <- scenario_fit("brms_slope", {
    fit <- brms::brm(
      pixel_z ~ 0 + Intercept + day_s + day_s2 + (1 + day_s | Dog) + (1 | Dog:Side), data = Pixel,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sd),
                brms::prior(normal(0, 1), class = sigma), brms::prior(lkj(1), class = cor)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99)
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })

  ### model summaries ----
  scenario_text("estimates-slope", JAGS_estimates_table(fit_slope))
  scenario_text("estimates-ar1",   JAGS_estimates_table(fit_ar1))
  scenario_text("estimates-car",   JAGS_estimates_table(fit_car))

  ### lme4 and brms comparisons ----
  lme4_fixed <- function(model) {
    fixed <- lme4::fixef(model)
    c(intercept = fixed[["(Intercept)"]], day_s = fixed[["day_s"]], day_s2 = fixed[["day_s2"]],
      sigma = stats::sigma(model),
      "Dog_Side: sd(intercept)" = attr(lme4::VarCorr(model)[["Dog:Side"]], "stddev")[["(Intercept)"]])
  }
  bayes_columns <- function(parameters) {
    stats::setNames(ifelse(parameters == "sigma", "sigma", paste0("(mu) ", parameters)), parameters)
  }
  slope_dog    <- lme4::VarCorr(lme4_slope)[["Dog"]]
  slope_values <- c(
    lme4_fixed(lme4_slope),
    "Dog: sd(intercept)" = attr(slope_dog, "stddev")[["(Intercept)"]], "Dog: sd(day_s)" = attr(slope_dog, "stddev")[["day_s"]],
    "Dog: cor(intercept,day_s)" = attr(slope_dog, "correlation")[2L, 1L]
  )
  scenario_comparison(
    "comparison-slope", names(slope_values), main = "Random intercepts and slopes per dog",
    references = list(lme4 = slope_values),
    posteriors = list(
      BayesTools = scenario_draws(JAGS_estimates_table(fit_slope, return_samples = TRUE), bayes_columns(names(slope_values))),
      brms       = scenario_draws(as.matrix(brms_slope[["fit"]]), c(
        intercept = "b_Intercept", day_s = "b_day_s", day_s2 = "b_day_s2", sigma = "sigma",
        "Dog_Side: sd(intercept)" = "sd_Dog:Side__Intercept", "Dog: sd(intercept)" = "sd_Dog__Intercept",
        "Dog: sd(day_s)" = "sd_Dog__day_s", "Dog: cor(intercept,day_s)" = "cor_Dog__Intercept__day_s"))
    )
  )
  # the CAR rows describe a different correlation function (cor^|day distance|
  # against cor^|index distance|), shown next to the AR(1) fits
  ar1_dog    <- as.matrix(lme4::VarCorr(lme4_ar1)[["Dog"]])
  ar1_values <- c(lme4_fixed(lme4_ar1), "Dog: sd" = sqrt(ar1_dog[1L, 1L]), "Dog: cor" = stats::cov2cor(ar1_dog)[2L, 1L])
  scenario_comparison(
    "comparison-ar1", names(ar1_values), main = "AR(1) and CAR day effects",
    references = list(lme4 = ar1_values),
    posteriors = list(
      BayesTools       = scenario_draws(JAGS_estimates_table(fit_ar1, return_samples = TRUE), bayes_columns(names(ar1_values))),
      "BayesTools CAR" = scenario_draws(JAGS_estimates_table(fit_car, return_samples = TRUE), bayes_columns(names(ar1_values)))
    )
  )

  ### exact marginal likelihoods ----
  # each structure's marginal likelihood reduces to a low-dimensional integral
  # over sigma, the SDs, and the correlation (LKJ(1) on a 2 x 2 matrix is
  # uniform(-1, 1); AR1 uniform(-1, 1); CAR uniform(0, 1)), which importance
  # sampling evaluates with an error far below that of bridge sampling
  exact_logml <- scenario_fit("exact_marglik", {
    same_dog      <- outer(Pixel$Dog, Pixel$Dog, "==") * 1
    same_dog_side <- same_dog * outer(Pixel$Side, Pixel$Side, "==")
    residual      <- diag(nrow(Pixel))
    index_lag     <- abs(outer(as.integer(Pixel$day_f), as.integer(Pixel$day_f), "-"))
    day_lag       <- abs(outer(Pixel$day, Pixel$day, "-"))
    # theta: sigma, the dog SDs, the dog-by-side SD, and the correlation
    covariances <- list(
      slope = function(t) {
        same_dog * (t[2]^2 + t[5] * t[2] * t[3] * outer(Pixel$day_s, Pixel$day_s, "+") + t[3]^2 * outer(Pixel$day_s, Pixel$day_s)) +
          t[4]^2 * same_dog_side + t[1]^2 * residual
      },
      ar1 = function(t) same_dog * t[2]^2 * t[4]^index_lag + t[3]^2 * same_dog_side + t[1]^2 * residual,
      car = function(t) same_dog * t[2]^2 * t[4]^day_lag   + t[3]^2 * same_dog_side + t[1]^2 * residual
    )
    columns <- list(
      slope = c("sigma", "(mu) Dog: sd(intercept)", "(mu) Dog: sd(day_s)", "(mu) Dog_Side: sd(intercept)",
                "(mu) Dog: cor(intercept,day_s)"),
      ar1   = c("sigma", "(mu) Dog: sd", "(mu) Dog_Side: sd(intercept)", "(mu) Dog: cor"),
      car   = c("sigma", "(mu) Dog: sd", "(mu) Dog_Side: sd(intercept)", "(mu) Dog: cor")
    )
    bounds <- list(
      slope = list("(mu) Dog: cor(intercept,day_s)" = c(-1, 1)),
      ar1   = list("(mu) Dog: cor" = c(-1, 1)),
      car   = list("(mu) Dog: cor" = c(0, 1))
    )
    fits <- list(slope = fit_slope, ar1 = fit_ar1, car = fit_car)
    set.seed(1)
    vapply(names(covariances), function(structure) {
      draws <- JAGS_estimates_table(fits[[structure]], return_samples = TRUE)[, columns[[structure]]]
      scenario_gaussian_logml(
        Pixel$pixel_z, stats::model.matrix(~ day_s + day_s2, Pixel), covariances[[structure]], draws,
        correlation_bounds = bounds[[structure]]
      )[["logml"]]
    }, numeric(1))
  })

  ### structure comparison ----
  scenario_text("structure-comparison", {
    out <- data.frame(
      structure        = c("(1 + day_s | Dog)", "ar1(day_f | Dog)", "car(day | Dog)"),
      lme4_logLik      = c(as.numeric(stats::logLik(lme4_slope)), as.numeric(stats::logLik(lme4_ar1)), NA_real_),
      exact_logml      = unname(exact_logml[c("slope", "ar1", "car")]),
      BayesTools_logml = c(slope[["marglik"]][["logml"]], ar1[["marglik"]][["logml"]], car[["marglik"]][["logml"]]),
      brms_logml       = c(brms_slope[["marglik"]][["logml"]], NA_real_, NA_real_)
    )
    out$BayesTools_log_BF_vs_slope <- out$BayesTools_logml - out$BayesTools_logml[[1L]]
    out[, -1] <- round(out[, -1], 3)
    out
  })

  ### fitted values and unobserved dog-days ----
  # dog-days without observations (left side): the conditional predictions
  # include the day effects drawn from their Gaussian conditional given the
  # dog's observed days; lme4 predicts with the BLUPs of the same effects
  unobserved <- expand.grid(Dog = levels(Pixel$Dog), day = sort(unique(Pixel$day)))
  unobserved <- unobserved[!paste(unobserved$Dog, unobserved$day) %in% paste(Pixel$Dog, Pixel$day), ]
  unobserved$Side   <- factor("L", levels = levels(Pixel$Side))
  unobserved$Dog    <- factor(unobserved$Dog, levels = levels(Pixel$Dog))
  unobserved$day_f  <- factor(unobserved$day, levels = levels(Pixel$day_f))
  unobserved$day_s  <- (unobserved$day - mean(Pixel$day)) / stats::sd(Pixel$day)
  unobserved$day_s2 <- unobserved$day_s^2
  predicted_mean <- function(fit, data = NULL) {
    rowMeans(JAGS_predict_formula(fit, "mu", data = data, formula_target = "conditional", seed = 1)$value)
  }
  scenario_plot("fitted-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(stats::fitted(lme4_slope), predicted_mean(fit_slope), "Fitted (random slopes)")
    scenario_agreement_plot(colMeans(brms::posterior_epred(brms_slope[["fit"]])), predicted_mean(fit_slope),
                            "Fitted (random slopes)", reference_label = "brms")
    scenario_agreement_plot(stats::fitted(lme4_ar1), predicted_mean(fit_ar1), "Fitted (AR1)")
    scenario_agreement_plot(stats::predict(lme4_ar1, newdata = unobserved), predicted_mean(fit_ar1, unobserved),
                            "Unobserved dog-days (AR1)")
  })
  scenario_plot("car-comparison", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(predicted_mean(fit_ar1), predicted_mean(fit_car), "Fitted",
                            reference_label = "AR1", estimate_label = "CAR")
    scenario_agreement_plot(predicted_mean(fit_ar1, unobserved), predicted_mean(fit_car, unobserved),
                            "Unobserved dog-days", reference_label = "AR1", estimate_label = "CAR")
  })

  ### prior and posterior distributions ----
  scenario_plot("posterior-cor", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    fits <- list("AR1 cor" = fit_ar1, "CAR cor" = fit_car)
    for (main in names(fits)) {
      selection <- parameter_catalog_resolve(parameter_catalog(fits[[main]]), "(mu) Dog: cor")
      plot_posterior(stats::setNames(list(parameter_mixed_posterior(fits[[main]], selection)), "(mu) Dog: cor"),
                     "(mu) Dog: cor", prior = TRUE, main = main)
    }
  })
})
