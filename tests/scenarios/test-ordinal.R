if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("ordinal")
# testthat::test_file("tests/scenarios/test-ordinal.R")

### Description
# ordered-factor (monotone) predictors: the nlme Gun data (rounds fired by
# nine teams of slight, average, or heavy physique, with two loading methods
# and random team intercepts) and the lme4 cake data (breaking angle at six
# ordered baking temperatures, three recipes, random recipe-by-replicate
# intercepts). The ordered prior puts a N(0, 1) prior on the total effect of
# the highest level and a flat Dirichlet on the shares of the increments, so
# a level's effect is the total times the cumulative share; brms monotonic
# effects mo() with b ~ N(0, 1 / D) and simo ~ Dirichlet(1) are the same
# model. Compared with brms and lme4, against unordered (mean-difference)
# and linear codings, with Bayes factors from bridge sampling, a
# spike-and-slab prior on the ordered total, and, for Gun, exact marginal
# likelihoods (the fixed and random effects integrated analytically, the SDs
# and the increment share importance-sampled)

testthat::test_that("Ordered-factor predictors for the Gun and cake data", {

  ### Gun data ----
  data("Gun", package = "nlme")
  Gun <- as.data.frame(Gun)
  Gun$rounds_z <- as.numeric(scale(Gun$rounds))
  Gun$Physique <- factor(as.character(Gun$Physique), levels = c("Slight", "Average", "Heavy"))
  Gun$Team     <- factor(as.character(Gun$Team))
  method_contrast   <- contr.meandif(levels(Gun$Method))
  physique_contrast <- contr.meandif(levels(Gun$Physique))
  Gun$method_c  <- method_contrast[as.integer(Gun$Method), 1L]
  Gun$physique_1 <- physique_contrast[as.integer(Gun$Physique), 1L]
  Gun$physique_2 <- physique_contrast[as.integer(Gun$Physique), 2L]
  Gun$Physique_mo <- factor(Gun$Physique, ordered = TRUE)

  gun_likelihood <- "model{
    for(i in 1:N){
      rounds_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  gun_data <- list(rounds_z = Gun$rounds_z, N = nrow(Gun))
  gun_log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["rounds_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  unit_prior    <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  meandif_prior <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
  gun_marglik <- function(fit) {
    JAGS_bridgesampling(fit, log_posterior = gun_log_posterior, data = gun_data, prior_list = list(sigma = unit_sd_prior), seed = 1)
  }

  gun <- list()
  gun[["ordered"]] <- scenario_fit("fit_gun_ordered", {
    fit <- JAGS_fit(gun_likelihood, gun_data, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ Method + Physique + (1 | Team)), formula_data_list = list(mu = Gun),
                    formula_prior_list = list(mu = list(intercept = unit_prior, Method = meandif_prior, Physique = prior_ordered(unit_prior))),
                    formula_random_prior_list = list(mu = prior_random(Team = random_block(sd = unit_sd_prior))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 10000, seed = 1)
    list(fit = fit, marglik = gun_marglik(fit))
  })
  gun[["unordered"]] <- scenario_fit("fit_gun_unordered", {
    fit <- JAGS_fit(gun_likelihood, gun_data, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ Method + Physique + (1 | Team)), formula_data_list = list(mu = Gun),
                    formula_prior_list = list(mu = list(intercept = unit_prior, Method = meandif_prior, Physique = meandif_prior)),
                    formula_random_prior_list = list(mu = prior_random(Team = random_block(sd = unit_sd_prior))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 10000, seed = 1)
    list(fit = fit, marglik = gun_marglik(fit))
  })
  gun[["null"]] <- scenario_fit("fit_gun_null", {
    fit <- JAGS_fit(gun_likelihood, gun_data, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ Method + (1 | Team)), formula_data_list = list(mu = Gun),
                    formula_prior_list = list(mu = list(intercept = unit_prior, Method = meandif_prior)),
                    formula_random_prior_list = list(mu = prior_random(Team = random_block(sd = unit_sd_prior))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 10000, seed = 1)
    list(fit = fit, marglik = gun_marglik(fit))
  })
  # spike-and-slab on the ordered total: one fit includes or excludes physique
  fit_gun_spike <- scenario_fit("fit_gun_spike", {
    JAGS_fit(gun_likelihood, gun_data, prior_list = list(sigma = unit_sd_prior),
             formula_list = list(mu = ~ Method + Physique + (1 | Team)), formula_data_list = list(mu = Gun),
             formula_prior_list = list(mu = list(intercept = unit_prior, Method = meandif_prior,
                                                 Physique = prior_ordered(prior_spike_and_slab(unit_prior)))),
             formula_random_prior_list = list(mu = prior_random(Team = random_block(sd = unit_sd_prior))),
             chains = 4, adapt = 1000, burnin = 2000, sample = 20000, seed = 1)
  })

  # brms: mo() with D = 2 increments, so b ~ N(0, 1/2) gives the N(0, 1)
  # total; the unordered model uses the mean-difference columns
  gun_brms_priors <- c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sd),
                       brms::prior(normal(0, 1), class = sigma))
  brms_marglik <- function(fit) {
    set.seed(1)
    brms::bridge_sampler(fit, silent = TRUE)
  }
  gun_brms <- list()
  gun_brms[["ordered"]] <- scenario_fit("brms_gun_ordered", {
    fit <- brms::brm(rounds_z ~ 0 + Intercept + method_c + mo(Physique_mo) + (1 | Team), data = Gun,
                     prior = c(gun_brms_priors, brms::prior(normal(0, 0.5), class = b, coef = moPhysique_mo),
                               brms::prior(dirichlet(1), class = simo, coef = moPhysique_mo1)),
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
                     control = list(adapt_delta = 0.99))
    list(fit = fit, marglik = brms_marglik(fit))
  })
  gun_brms[["unordered"]] <- scenario_fit("brms_gun_unordered", {
    fit <- brms::brm(rounds_z ~ 0 + Intercept + method_c + physique_1 + physique_2 + (1 | Team), data = Gun,
                     prior = gun_brms_priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
                     control = list(adapt_delta = 0.99))
    list(fit = fit, marglik = brms_marglik(fit))
  })
  gun_brms[["null"]] <- scenario_fit("brms_gun_null", {
    fit <- brms::brm(rounds_z ~ 0 + Intercept + method_c + (1 | Team), data = Gun,
                     prior = gun_brms_priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
                     control = list(adapt_delta = 0.99))
    list(fit = fit, marglik = brms_marglik(fit))
  })
  gun_lme4 <- lme4::lmer(rounds_z ~ Method + Physique + (1 | Team), data = Gun, REML = FALSE)

  # exact marginal likelihoods: with the fixed coefficients N(0, 1) and the
  # team intercepts integrated out, y ~ N(0, X X' + tau^2 Z Z' + sigma^2 I);
  # for the ordered model the physique column is 0, s, 1 for the three levels
  # with the share s ~ uniform(0, 1)
  gun_exact <- scenario_fit("exact_gun", {
    same_team <- outer(Gun$Team, Gun$Team, "==") * 1
    residual  <- diag(nrow(Gun))
    base      <- cbind(1, Gun$method_c)
    level     <- as.integer(Gun$Physique)
    sd_columns <- c("sigma", "(mu) sd(intercept)")
    draws <- function(name, columns = sd_columns) JAGS_estimates_table(gun[[name]][["fit"]], return_samples = TRUE)[, columns]
    set.seed(1)
    c(
      ordered = scenario_gaussian_logml(
        Gun$rounds_z, matrix(0, nrow(Gun), 0L),
        function(t) {
          X <- cbind(base, c(0, t[3], 1)[level])
          tcrossprod(X) + t[2]^2 * same_team + t[1]^2 * residual
        },
        draws("ordered", c(sd_columns, "(mu) Physique_ordered_allocation[Average]")),
        correlation_bounds = list("(mu) Physique_ordered_allocation[Average]" = c(0, 1))
      )[["logml"]],
      unordered = scenario_gaussian_logml(
        Gun$rounds_z, cbind(base, Gun$physique_1, Gun$physique_2),
        function(t) t[2]^2 * same_team + t[1]^2 * residual, draws("unordered")
      )[["logml"]],
      null = scenario_gaussian_logml(
        Gun$rounds_z, base, function(t) t[2]^2 * same_team + t[1]^2 * residual, draws("null")
      )[["logml"]]
    )
  })

  scenario_text("estimates-gun-ordered", JAGS_estimates_table(gun[["ordered"]][["fit"]]))
  scenario_text("estimates-gun-ordered-levels", JAGS_estimates_table(gun[["ordered"]][["fit"]], transform_factors = TRUE))
  scenario_text("inference-gun-spike", JAGS_inference_table(fit_gun_spike))

  # physique effects relative to slight: lme4 treatment coefficients, the
  # ordered levels, the unordered level differences, and brms mo() (b D times
  # the cumulative simplex)
  scenario_text("comparison-gun-physique", {
    ordered   <- JAGS_estimates_table(gun[["ordered"]][["fit"]], transform_factors = TRUE, return_samples = TRUE)
    unordered <- JAGS_estimates_table(gun[["unordered"]][["fit"]], transform_factors = TRUE, return_samples = TRUE)
    unordered_levels <- unordered[, c("(mu) Physique[dif: Average]", "(mu) Physique[dif: Heavy]")] -
      unordered[, "(mu) Physique[dif: Slight]"]
    brms_draws <- as.matrix(gun_brms[["ordered"]][["fit"]])
    brms_levels <- cbind(2 * brms_draws[, "bsp_moPhysique_mo"] * brms_draws[, "simo_moPhysique_mo1[1]"],
                         2 * brms_draws[, "bsp_moPhysique_mo"])
    parameters <- c("Average - Slight", "Heavy - Slight", "total share of Average")
    scenario_comparison_table(
      parameters,
      references = list(lme4 = c("Average - Slight" = lme4::fixef(gun_lme4)[["PhysiqueAverage"]],
                                 "Heavy - Slight" = lme4::fixef(gun_lme4)[["PhysiqueHeavy"]])),
      posteriors = list(
        ordered = cbind(scenario_draws(ordered, c("Average - Slight" = "(mu) Physique[Average]", "Heavy - Slight" = "(mu) Physique[Heavy]")),
                        "total share of Average" = JAGS_estimates_table(gun[["ordered"]][["fit"]], return_samples = TRUE)[, "(mu) Physique_ordered_allocation[Average]"]),
        "brms mo()" = cbind("Average - Slight" = brms_levels[, 1L], "Heavy - Slight" = brms_levels[, 2L],
                            "total share of Average" = brms_draws[, "simo_moPhysique_mo1[1]"]),
        unordered = stats::setNames(as.data.frame(unordered_levels), parameters[1:2])
      )
    )
  })

  scenario_text("bayes-factor-gun", {
    log_bf <- function(fits, h1, h0) scenario_log_bf_bridge(fits[[h1]][["marglik"]], fits[[h0]][["marglik"]])
    rbind(
      scenario_bf_table("ordered physique", "no physique effect", c(
        "exact" = gun_exact[["ordered"]] - gun_exact[["null"]],
        "bridge sampling" = log_bf(gun, "ordered", "null"),
        "bridge sampling (brms)" = log_bf(gun_brms, "ordered", "null"),
        "spike and slab" = log(as.numeric(JAGS_inference_table(fit_gun_spike)[["inclusion_BF"]])))),
      scenario_bf_table("unordered physique", "no physique effect", c(
        "exact" = gun_exact[["unordered"]] - gun_exact[["null"]],
        "bridge sampling" = log_bf(gun, "unordered", "null"),
        "bridge sampling (brms)" = log_bf(gun_brms, "unordered", "null"))),
      scenario_bf_table("ordered physique", "unordered physique", c(
        "exact" = gun_exact[["ordered"]] - gun_exact[["unordered"]],
        "bridge sampling" = log_bf(gun, "ordered", "unordered"),
        "bridge sampling (brms)" = log_bf(gun_brms, "ordered", "unordered")))
    )
  })
  # hypotheses on the ordered levels: Heavy is the total, so its Savage-Dickey
  # BF against 0 is the BF of the ordered model against the null (exact
  # exp(-0.454) = 0.635); Average = share x Heavy with the share in (0, 1), so
  # Average < Heavy is the complement of Heavy < 0 and its BF is the inverse.
  # Point nulls on an intermediate level are refused: a Beta(1, 1) share gives
  # an infinite prior ordinate.
  scenario_text("hypotheses-gun-ordered", {
    hypothesis_BF(gun[["ordered"]][["fit"]], c("Physique[Heavy] = 0", "Physique[Heavy] < 0", "Physique[Average] < Physique[Heavy]"), seed = 1)
  })
  scenario_text("predictions-gun", {
    grid <- expand.grid(Method = levels(Gun$Method), Physique = levels(Gun$Physique))
    grid$Method   <- factor(grid$Method, levels = levels(Gun$Method))
    grid$Physique <- factor(grid$Physique, levels = levels(Gun$Physique))
    grid$Physique_mo <- factor(grid$Physique, levels = levels(Gun$Physique), ordered = TRUE)
    grid$method_c <- method_contrast[as.integer(grid$Method), 1L]
    scenario_prediction_table(grid[, c("Method", "Physique")], list(
      "lme4"       = stats::predict(gun_lme4, newdata = grid, re.form = NA),
      "BayesTools" = JAGS_predict_formula(gun[["ordered"]][["fit"]], "mu", data = grid[, c("Method", "Physique")], formula_target = "fixed")$value,
      "brms mo()"  = t(brms::posterior_epred(gun_brms[["ordered"]][["fit"]], newdata = grid, re_formula = NA))
    ))
  })

  ### cake data ----
  data("cake", package = "lme4")
  cake$angle_z    <- as.numeric(scale(cake$angle))
  cake$block      <- interaction(cake$recipe, cake$replicate, drop = TRUE)
  cake$temp_s     <- (cake$temp - mean(cake$temp)) / stats::sd(cake$temp)
  cake$temperature_f <- factor(as.character(cake$temperature), levels = levels(cake$temperature))
  recipe_contrast      <- contr.meandif(levels(cake$recipe))
  temperature_contrast <- contr.meandif(levels(cake$temperature_f))
  cake$recipe_1 <- recipe_contrast[as.integer(cake$recipe), 1L]
  cake$recipe_2 <- recipe_contrast[as.integer(cake$recipe), 2L]
  for (j in seq_len(ncol(temperature_contrast))) {
    cake[[paste0("temperature_", j)]] <- temperature_contrast[as.integer(cake$temperature_f), j]
  }

  cake_likelihood <- "model{
    for(i in 1:N){
      angle_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  cake_data <- list(angle_z = cake$angle_z, N = nrow(cake))
  cake_log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["angle_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  cake_marglik <- function(fit) {
    JAGS_bridgesampling(fit, log_posterior = cake_log_posterior, data = cake_data, prior_list = list(sigma = unit_sd_prior), seed = 1)
  }
  cake_fits <- list()
  cake_fits[["ordered"]] <- scenario_fit("fit_cake_ordered", {
    fit <- JAGS_fit(cake_likelihood, cake_data, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ recipe + temperature_f + (1 | block)), formula_data_list = list(mu = cake),
                    formula_prior_list = list(mu = list(intercept = unit_prior, recipe = meandif_prior, temperature_f = prior_ordered(unit_prior))),
                    formula_random_prior_list = list(mu = prior_random(block = random_block(sd = unit_sd_prior, parameterization = "mean_centered"))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 10000, seed = 1)
    list(fit = fit, marglik = cake_marglik(fit))
  })
  cake_fits[["unordered"]] <- scenario_fit("fit_cake_unordered", {
    fit <- JAGS_fit(cake_likelihood, cake_data, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ recipe + temperature_f + (1 | block)), formula_data_list = list(mu = cake),
                    formula_prior_list = list(mu = list(intercept = unit_prior, recipe = meandif_prior, temperature_f = meandif_prior)),
                    formula_random_prior_list = list(mu = prior_random(block = random_block(sd = unit_sd_prior, parameterization = "mean_centered"))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 10000, seed = 1)
    list(fit = fit, marglik = cake_marglik(fit))
  })
  cake_fits[["linear"]] <- scenario_fit("fit_cake_linear", {
    fit <- JAGS_fit(cake_likelihood, cake_data, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ recipe + temp + (1 | block)), formula_data_list = list(mu = cake),
                    formula_prior_list = list(mu = list(intercept = unit_prior, recipe = meandif_prior, temp = unit_prior)),
                    formula_scale_list = list(mu = list(temp = TRUE)),
                    formula_random_prior_list = list(mu = prior_random(block = random_block(sd = unit_sd_prior, parameterization = "mean_centered"))),
                    chains = 4, adapt = 1000, burnin = 2000, sample = 10000, seed = 1)
    list(fit = fit, marglik = cake_marglik(fit))
  })

  # brms: mo() with D = 5 increments, so b ~ N(0, 1/5)
  cake_brms <- list()
  cake_brms[["ordered"]] <- scenario_fit("brms_cake_ordered", {
    fit <- brms::brm(angle_z ~ 0 + Intercept + recipe_1 + recipe_2 + mo(temperature) + (1 | block), data = cake,
                     prior = c(gun_brms_priors, brms::prior(normal(0, 0.2), class = b, coef = motemperature),
                               brms::prior(dirichlet(1), class = simo, coef = motemperature1)),
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  cake_brms[["unordered"]] <- scenario_fit("brms_cake_unordered", {
    fit <- brms::brm(angle_z ~ 0 + Intercept + recipe_1 + recipe_2 + temperature_1 + temperature_2 + temperature_3 +
                       temperature_4 + temperature_5 + (1 | block), data = cake,
                     prior = gun_brms_priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  cake_brms[["linear"]] <- scenario_fit("brms_cake_linear", {
    fit <- brms::brm(angle_z ~ 0 + Intercept + recipe_1 + recipe_2 + temp_s + (1 | block), data = cake,
                     prior = gun_brms_priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  cake_lme4 <- lme4::lmer(angle_z ~ recipe + temperature_f + (1 | block), data = cake, REML = FALSE,
                          contrasts = list(temperature_f = "contr.treatment"))

  scenario_text("estimates-cake-ordered-levels", JAGS_estimates_table(cake_fits[["ordered"]][["fit"]], transform_factors = TRUE))

  # temperature effects relative to 175 degrees
  scenario_text("comparison-cake-temperature", {
    levels_t   <- levels(cake$temperature_f)[-1L]
    parameters <- paste0(levels_t, " - 175")
    ordered    <- JAGS_estimates_table(cake_fits[["ordered"]][["fit"]], transform_factors = TRUE, return_samples = TRUE)
    brms_draws <- as.matrix(cake_brms[["ordered"]][["fit"]])
    simplex    <- brms_draws[, paste0("simo_motemperature1[", 1:5, "]")]
    brms_levels <- 5 * brms_draws[, "bsp_motemperature"] * t(apply(simplex, 1, cumsum))
    colnames(brms_levels) <- parameters
    lme4_levels <- lme4::fixef(cake_lme4)[paste0("temperature_f", levels_t)]
    scenario_comparison_table(
      parameters,
      references = list(lme4 = stats::setNames(lme4_levels, parameters)),
      posteriors = list(
        ordered     = scenario_draws(ordered, stats::setNames(paste0("(mu) temperature_f[", levels_t, "]"), parameters)),
        "brms mo()" = brms_levels
      )
    )
  })
  scenario_text("bayes-factor-cake", {
    log_bf <- function(fits, h1, h0) scenario_log_bf_bridge(fits[[h1]][["marglik"]], fits[[h0]][["marglik"]])
    rbind(
      scenario_bf_table("ordered temperature", "linear temperature", c(
        "bridge sampling" = log_bf(cake_fits, "ordered", "linear"),
        "bridge sampling (brms)" = log_bf(cake_brms, "ordered", "linear"))),
      scenario_bf_table("ordered temperature", "unordered temperature", c(
        "bridge sampling" = log_bf(cake_fits, "ordered", "unordered"),
        "bridge sampling (brms)" = log_bf(cake_brms, "ordered", "unordered")))
    )
  })

  scenario_plot("level-effects", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    for (case in list(list(fit = gun[["ordered"]][["fit"]], term = "Physique", levels = levels(Gun$Physique), main = "Gun: physique"),
                      list(fit = cake_fits[["ordered"]][["fit"]], term = "temperature_f", levels = levels(cake$temperature_f),
                           main = "cake: temperature"))) {
      samples <- JAGS_estimates_table(case$fit, transform_factors = TRUE, return_samples = TRUE)
      columns <- grep(paste0("^\\(mu\\) ", case$term, "\\["), colnames(samples), value = TRUE)
      draws   <- cbind(0, samples[, columns, drop = FALSE])
      summary <- apply(draws, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
      graphics::plot(seq_len(ncol(draws)), summary[2L, ], ylim = range(summary), pch = 19, xaxt = "n",
                     xlab = "level", ylab = "effect relative to the first level", main = case$main)
      graphics::axis(1, at = seq_along(case$levels), labels = case$levels)
      graphics::segments(seq_len(ncol(draws)), summary[1L, ], seq_len(ncol(draws)), summary[3L, ])
      graphics::abline(h = 0, lty = 2, col = "grey50")
    }
  })
})
