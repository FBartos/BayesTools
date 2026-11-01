if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("linear")
# testthat::test_file("tests/scenarios/test-linear.R")

### Description
# fixed-effect linear models for the ToothGrowth data (tooth length by
# supplement and vitamin C dose): a mean-difference supplement factor with a
# standardized continuous dose (formula_scale_list) and their interaction, and
# orthonormal supplement and three-level dose factors with their interaction.
# With N(0, 1) priors on the intercept and every contrast coefficient and a
# half-normal(0, 1) residual SD, the posterior given sigma is normal and the
# marginal likelihood is a one-dimensional integral over sigma, so posterior
# means and SDs, predictions, marginal likelihoods, Savage-Dickey ratios, and
# posterior model probabilities have exact references. Compared with lm (ML),
# brms under the same priors, bridge sampling, Savage-Dickey (hypothesis_BF),
# a spike-and-slab prior on the interaction, region hypotheses, and model
# averaging over five nested models

testthat::test_that("ToothGrowth fixed-effect linear models", {

  tg <- ToothGrowth
  tg$len_z  <- as.numeric(scale(tg$len))
  tg$dose_f <- factor(tg$dose)
  # the BayesTools contrast codings as numeric columns for lm, brms, and the
  # exact reference: meandif supplement, standardized dose, orthonormal dose
  meandif     <- contr.meandif(levels(tg$supp))
  orthonormal <- contr.orthonormal(levels(tg$dose_f))
  tg$supp_c   <- meandif[as.integer(tg$supp), 1L]
  tg$dose_s   <- (tg$dose - mean(tg$dose)) / stats::sd(tg$dose)
  tg$dose_1   <- orthonormal[as.integer(tg$dose_f), 1L]
  tg$dose_2   <- orthonormal[as.integer(tg$dose_f), 2L]
  tg$supp_o   <- contr.orthonormal(levels(tg$supp))[as.integer(tg$supp), 1L]

  likelihood <- "model{
    for(i in 1:N){
      len_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  data_list <- list(len_z = tg$len_z, N = nrow(tg))
  log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["len_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  unit_prior    <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  meandif_prior <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
  ortho_prior   <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  marglik <- function(fit) {
    JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1)
  }

  ### exact reference ----
  # y | sigma ~ N(0, X X' + sigma^2 I) with beta ~ N(0, I); given sigma, beta
  # is normal with covariance V = (X'X / sigma^2 + I)^-1 and mean
  # V X'y / sigma^2. Everything is integrated over sigma on a 0.0005 grid up
  # to 4 (the posterior SD of sigma is about 0.05)
  designs <- list(
    null = cbind(intercept = rep(1, nrow(tg))),
    supp = cbind(intercept = 1, supp = tg$supp_c),
    dose = cbind(intercept = 1, dose = tg$dose_s),
    add  = cbind(intercept = 1, supp = tg$supp_c, dose = tg$dose_s),
    full = cbind(intercept = 1, supp = tg$supp_c, dose = tg$dose_s, "supp:dose" = tg$supp_c * tg$dose_s),
    factor_add  = cbind(intercept = 1, supp = tg$supp_o, "dose_f{1}" = tg$dose_1, "dose_f{2}" = tg$dose_2),
    factor_full = cbind(intercept = 1, supp = tg$supp_o, "dose_f{1}" = tg$dose_1, "dose_f{2}" = tg$dose_2,
                        "supp:dose_f{1}" = tg$supp_o * tg$dose_1, "supp:dose_f{2}" = tg$supp_o * tg$dose_2)
  )
  exact_linear <- function(X, y = tg$len_z) {
    sigma  <- seq(0.0005, 4, by = 0.0005)
    gram   <- eigen(crossprod(X), symmetric = TRUE)
    outer  <- eigen(tcrossprod(X), symmetric = TRUE)
    proj   <- drop(crossprod(outer$vectors, y))
    lambda <- pmax(outer$values, 0)
    log_weight <- vapply(sigma, function(s) {
      -0.5 * (length(y) * log(2 * pi) + sum(log(lambda + s^2)) + sum(proj^2 / (lambda + s^2)))
    }, numeric(1)) + log(2) + stats::dnorm(sigma, log = TRUE)
    scale  <- max(log_weight)
    weight <- exp(log_weight - scale)
    total  <- sum(weight) * 0.0005
    weight <- weight * 0.0005 / total
    Xy     <- drop(crossprod(X, y))
    moments <- lapply(sigma, function(s) {
      V <- gram$vectors %*% (t(gram$vectors) / (gram$values / s^2 + 1))
      list(mean = drop(V %*% Xy) / s^2, cov = V)
    })
    beta_mean <- Reduce(`+`, Map(function(m, w) w * m$mean, moments, weight))
    beta_cov  <- Reduce(`+`, Map(function(m, w) w * (m$cov + tcrossprod(m$mean)), moments, weight)) - tcrossprod(beta_mean)
    density_at_zero <- vapply(seq_len(ncol(X)), function(j) {
      sum(weight * vapply(moments, function(m) stats::dnorm(0, m$mean[j], sqrt(m$cov[j, j])), numeric(1)))
    }, numeric(1))
    names(beta_mean) <- names(density_at_zero) <- colnames(X)
    list(logml = scale + log(total), mean = c(beta_mean, sigma = sum(weight * sigma)),
         sd = c(sqrt(diag(beta_cov)), sigma = sqrt(sum(weight * sigma^2) - sum(weight * sigma)^2)),
         cov = beta_cov, density_at_zero = density_at_zero)
  }
  exact <- scenario_fit("exact_reference", lapply(designs, exact_linear))

  ### BayesTools fits ----
  # five nested models with the mean-difference supplement and the scaled
  # dose, the full model with a spike-and-slab interaction prior, and the
  # orthonormal factor models; each nested block returns its fit and
  # bridge-sampling marginal likelihood
  fits <- list()
  fits[["null"]] <- scenario_fit("fit_null", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ 1), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })
  fits[["supp"]] <- scenario_fit("fit_supp", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ supp), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior, supp = meandif_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })
  fits[["dose"]] <- scenario_fit("fit_dose", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ dose), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior, dose = unit_prior)),
                    formula_scale_list = list(mu = list(dose = TRUE)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })
  fits[["add"]] <- scenario_fit("fit_add", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ supp + dose), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior, supp = meandif_prior, dose = unit_prior)),
                    formula_scale_list = list(mu = list(dose = TRUE)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })
  fits[["full"]] <- scenario_fit("fit_full", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ supp * dose), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior, supp = meandif_prior, dose = unit_prior,
                                                        "supp:dose" = meandif_prior)),
                    formula_scale_list = list(mu = list(dose = TRUE)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })
  fit_spike <- scenario_fit("fit_full_spike", {
    JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
             formula_list = list(mu = ~ supp * dose), formula_data_list = list(mu = tg),
             formula_prior_list = list(mu = list(intercept = unit_prior, supp = meandif_prior, dose = unit_prior,
                                                 "supp:dose" = prior_spike_and_slab(meandif_prior))),
             formula_scale_list = list(mu = list(dose = TRUE)),
             chains = 4, adapt = 1000, burnin = 1000, sample = 20000, seed = 1)
  })
  fits[["factor_add"]] <- scenario_fit("fit_factor_add", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ supp + dose_f), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior, supp = ortho_prior, dose_f = ortho_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })
  fits[["factor_full"]] <- scenario_fit("fit_factor_full", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(sigma = unit_sd_prior),
                    formula_list = list(mu = ~ supp * dose_f), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = unit_prior, supp = ortho_prior, dose_f = ortho_prior,
                                                        "supp:dose_f" = ortho_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit))
  })

  ### lm and brms fits ----
  # brms uses the same numeric contrast columns, `0 + Intercept`, N(0, 1) on
  # every coefficient, and a half-normal(0, 1) sigma
  lm_full        <- stats::lm(len_z ~ supp_c * dose_s, data = tg)
  lm_full_raw    <- stats::lm(len_z ~ supp_c * dose, data = tg)
  lm_factor_full <- stats::lm(len_z ~ supp_o * (dose_1 + dose_2), data = tg)
  brms_priors    <- c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sigma))
  brms_marglik   <- function(fit) {
    set.seed(1)
    brms::bridge_sampler(fit, silent = TRUE)
  }
  brms_fits <- list()
  brms_fits[["dose"]] <- scenario_fit("brms_dose", {
    fit <- brms::brm(len_z ~ 0 + Intercept + dose_s, data = tg, prior = brms_priors,
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["add"]] <- scenario_fit("brms_add", {
    fit <- brms::brm(len_z ~ 0 + Intercept + supp_c + dose_s, data = tg, prior = brms_priors,
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["full"]] <- scenario_fit("brms_full", {
    fit <- brms::brm(len_z ~ 0 + Intercept + supp_c * dose_s, data = tg, prior = brms_priors,
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["factor_add"]] <- scenario_fit("brms_factor_add", {
    fit <- brms::brm(len_z ~ 0 + Intercept + supp_o + dose_1 + dose_2, data = tg, prior = brms_priors,
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["factor_full"]] <- scenario_fit("brms_factor_full", {
    fit <- brms::brm(len_z ~ 0 + Intercept + supp_o * (dose_1 + dose_2), data = tg, prior = brms_priors,
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })

  ### model summaries ----
  scenario_text("estimates-full",        JAGS_estimates_table(fits[["full"]][["fit"]]))
  scenario_text("estimates-factor-full", JAGS_estimates_table(fits[["factor_full"]][["fit"]]))
  scenario_text("inference-full-spike",  JAGS_inference_table(fit_spike))

  ### exact, lm, and brms comparisons ----
  full_columns <- c(intercept = "(mu) intercept", supp = "(mu) supp{1}", dose = "(mu) dose",
                    "supp:dose" = "(mu) supp:dose{1}", sigma = "sigma")
  scenario_text("comparison-full", {
    scenario_comparison_table(
      names(full_columns),
      references = list(
        exact = exact[["full"]][["mean"]],
        lm    = c(stats::setNames(stats::coef(lm_full), names(full_columns)[1:4]), sigma = stats::sigma(lm_full))
      ),
      posteriors = list(
        BayesTools = scenario_draws(JAGS_estimates_table(fits[["full"]][["fit"]], return_samples = TRUE), full_columns),
        brms       = scenario_draws(as.matrix(brms_fits[["full"]][["fit"]]), c(
          intercept = "b_Intercept", supp = "b_supp_c", dose = "b_dose_s", "supp:dose" = "b_supp_c:dose_s", sigma = "sigma"))
      )
    )
  })
  # exact posterior SDs: their agreement checks the posterior spread
  scenario_text("posterior-sd-full", {
    samples <- scenario_draws(JAGS_estimates_table(fits[["full"]][["fit"]], return_samples = TRUE), full_columns)
    data.frame(parameter = names(full_columns), exact = round(exact[["full"]][["sd"]], 4),
               BayesTools = round(apply(samples, 2, stats::sd), 4), row.names = NULL)
  })
  # original dose units: transform_scaled undoes the standardization; the
  # exact means transform linearly with the dose mean m and SD s
  scenario_text("comparison-full-original-scale", {
    m <- mean(tg$dose)
    s <- stats::sd(tg$dose)
    b <- exact[["full"]][["mean"]]
    exact_original <- c(intercept = b[["intercept"]] - b[["dose"]] * m / s, supp = b[["supp"]] - b[["supp:dose"]] * m / s,
                        dose = b[["dose"]] / s, "supp:dose" = b[["supp:dose"]] / s)
    scenario_comparison_table(
      names(exact_original),
      references = list(exact = exact_original,
                        lm = stats::setNames(stats::coef(lm_full_raw), names(exact_original))),
      posteriors = list(BayesTools = scenario_draws(
        JAGS_estimates_table(fits[["full"]][["fit"]], transform_scaled = TRUE, return_samples = TRUE), full_columns[1:4]))
    )
  })
  factor_columns <- c(intercept = "(mu) intercept", supp = "(mu) supp{1}", "dose_f{1}" = "(mu) dose_f{1}",
                      "dose_f{2}" = "(mu) dose_f{2}", "supp:dose_f{1}" = "(mu) supp:dose_f{1}",
                      "supp:dose_f{2}" = "(mu) supp:dose_f{2}", sigma = "sigma")
  scenario_text("comparison-factor-full", {
    scenario_comparison_table(
      names(factor_columns),
      references = list(
        exact = exact[["factor_full"]][["mean"]],
        lm    = c(stats::setNames(stats::coef(lm_factor_full), names(factor_columns)[1:6]), sigma = stats::sigma(lm_factor_full))
      ),
      posteriors = list(
        BayesTools = scenario_draws(JAGS_estimates_table(fits[["factor_full"]][["fit"]], return_samples = TRUE), factor_columns),
        brms       = scenario_draws(as.matrix(brms_fits[["factor_full"]][["fit"]]), c(
          intercept = "b_Intercept", supp = "b_supp_o", "dose_f{1}" = "b_dose_1", "dose_f{2}" = "b_dose_2",
          "supp:dose_f{1}" = "b_supp_o:dose_1", "supp:dose_f{2}" = "b_supp_o:dose_2", sigma = "sigma"))
      )
    )
  })

  ### Bayes factors ----
  # log BF10: exact, bridge sampling in BayesTools and brms, Savage-Dickey
  # (hypothesis_BF; exact: the prior over the posterior density at 0), and the
  # spike-and-slab inclusion BF from a single fit
  exact_log_bf <- function(h1, h0) exact[[h1]][["logml"]] - exact[[h0]][["logml"]]
  scenario_text("bayes-factor-interaction", {
    savage_dickey <- hypothesis_BF(fits[["full"]][["fit"]], "supp:dose = 0", seed = 1)
    scenario_bf_table(
      H1 = "supp * dose",
      H0 = "supp + dose",
      log_BF10 = c(
        "exact"                  = exact_log_bf("full", "add"),
        "bridge sampling"        = scenario_log_bf_bridge(fits[["full"]][["marglik"]], fits[["add"]][["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_fits[["full"]][["marglik"]], brms_fits[["add"]][["marglik"]]),
        "Savage-Dickey (exact)"  = log(stats::dnorm(0) / exact[["full"]][["density_at_zero"]][["supp:dose"]]),
        "Savage-Dickey"          = log(as.numeric(savage_dickey[["BF"]][1L])),
        "spike and slab"         = log(as.numeric(JAGS_inference_table(fit_spike)[["inclusion_BF"]]))
      )
    )
  })
  scenario_text("bayes-factor-supplement", {
    scenario_bf_table(
      H1 = "supp + dose",
      H0 = "dose",
      log_BF10 = c(
        "exact"                  = exact_log_bf("add", "dose"),
        "bridge sampling"        = scenario_log_bf_bridge(fits[["add"]][["marglik"]], fits[["dose"]][["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_fits[["add"]][["marglik"]], brms_fits[["dose"]][["marglik"]])
      )
    )
  })
  scenario_text("bayes-factor-factor-interaction", {
    scenario_bf_table(
      H1 = "supp * dose_f",
      H0 = "supp + dose_f",
      log_BF10 = c(
        "exact"                  = exact_log_bf("factor_full", "factor_add"),
        "bridge sampling"        = scenario_log_bf_bridge(fits[["factor_full"]][["marglik"]], fits[["factor_add"]][["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_fits[["factor_full"]][["marglik"]], brms_fits[["factor_add"]][["marglik"]])
      )
    )
  })
  # a region hypothesis per supplement level: the VC supplement is below the
  # mean length difference (its meandif cell is negative); the exact
  # posterior probability is a normal mixture over sigma
  scenario_text("hypothesis-supplement-region", {
    hypothesis_BF(fits[["full"]][["fit"]], "supp < 0", seed = 1)
  })

  ### model averaging ----
  # posterior model probabilities of the five nested models (equal prior
  # probabilities) from the bridge-sampling marginal likelihoods and exactly
  scenario_text("model-probabilities", {
    nested <- c("null", "supp", "dose", "add", "full")
    models <- models_inference(lapply(nested, function(name) {
      list(model = fits[[name]][["fit"]], marglik = fits[[name]][["marglik"]], prior_weights = 1)
    }))
    exact_logml <- vapply(nested, function(name) exact[[name]][["logml"]], numeric(1))
    data.frame(
      model      = nested,
      exact      = round(exp(exact_logml - max(exact_logml)) / sum(exp(exact_logml - max(exact_logml))), 4),
      BayesTools = round(vapply(models, function(model) model[["inference"]][["post_prob"]], numeric(1)), 4),
      exact_logml      = round(exact_logml, 3),
      BayesTools_logml = round(vapply(nested, function(name) fits[[name]][["marglik"]][["logml"]], numeric(1)), 3),
      row.names = NULL
    )
  })

  ### predictions ----
  # population means on a supplement-by-dose grid: the exact posterior mean
  # and SD of x'beta, lm, and brms
  scenario_text("predictions-grid", {
    grid <- expand.grid(supp = factor(levels(tg$supp), levels = levels(tg$supp)), dose = c(0.5, 1, 1.5, 2))
    grid$supp_c <- meandif[as.integer(grid$supp), 1L]
    grid$dose_s <- (grid$dose - mean(tg$dose)) / stats::sd(tg$dose)
    X_grid <- cbind(1, grid$supp_c, grid$dose_s, grid$supp_c * grid$dose_s)
    scenario_prediction_table(grid[, c("supp", "dose")], list(
      "exact"         = drop(X_grid %*% exact[["full"]][["mean"]][1:4]),
      "exact sd"      = sqrt(rowSums((X_grid %*% exact[["full"]][["cov"]]) * X_grid)),
      "lm"            = stats::predict(lm_full, newdata = grid),
      "BayesTools"    = JAGS_predict_formula(fits[["full"]][["fit"]], "mu", data = grid[, c("supp", "dose")], formula_target = "fixed")$value,
      "brms"          = t(brms::posterior_epred(brms_fits[["full"]][["fit"]], newdata = grid))
    ))
  })

  scenario_plot("fitted-comparison", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    fitted_mean <- function(fit) rowMeans(JAGS_predict_formula(fit, "mu", formula_target = "fixed")$value)
    scenario_agreement_plot(drop(designs[["full"]] %*% exact[["full"]][["mean"]][1:4]), fitted_mean(fits[["full"]][["fit"]]),
                            "supp * dose", reference_label = "exact")
    scenario_agreement_plot(drop(designs[["factor_full"]] %*% exact[["factor_full"]][["mean"]][1:6]),
                            fitted_mean(fits[["factor_full"]][["fit"]]), "supp * dose_f", reference_label = "exact")
  })
})
