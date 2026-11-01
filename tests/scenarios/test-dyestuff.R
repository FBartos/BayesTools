if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("dyestuff")
# testthat::test_file("tests/scenarios/test-dyestuff.R")

### Description
# one-way random intercepts for the lme4 Dyestuff data (a clearly positive
# batch SD) and Dyestuff2 data (the lme4 ML batch SD is 0, a singular fit),
# checked against an exact reference: with the intercept integrated
# analytically, the marginal likelihoods, the posterior means, and the
# posterior distribution of the batch SD tau are two-dimensional integrals over
# (sigma, tau) evaluated on a grid. Compared with BayesTools bridge sampling,
# the spike-and-slab inclusion BF, and Savage-Dickey BFs for tau = 0 under the
# noncentered and the mean-centered parameterization, with brms under matched
# priors, and with lme4. The mean-centered parameterization samples tau near 0
# sparsely, which inflates its Savage-Dickey BF; the diagnostics table shows
# which MCMC diagnostics reveal it

testthat::test_that("Dyestuff one-way random intercepts against an exact reference", {

  datasets <- list()
  for (name in c("Dyestuff", "Dyestuff2")) {
    env <- new.env()
    utils::data(list = name, package = "lme4", envir = env)
    datasets[[name]] <- env[[name]]
    datasets[[name]]$Yield_z <- as.numeric(scale(datasets[[name]]$Yield))
  }

  likelihood <- "model{
    for(i in 1:N){
      Yield_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["Yield_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  data_list <- function(data) list(Yield_z = data$Yield_z, N = nrow(data))

  ### exact reference ----
  # priors: intercept ~ N(0, 1), sigma and tau ~ half-normal(0, 1). Given
  # (sigma, tau), y ~ N(0, V + 1 1') with V = sigma^2 I + tau^2 Z Z', and the
  # intercept posterior mean is 1'V^-1 y / (1 + 1'V^-1 1). The (sigma, tau)
  # integrals use trapezoid weights on a 0.01 grid up to 4, where the priors on
  # the standardized scale leave negligible mass; the integrand vanishes at
  # sigma = 0
  exact_reference <- function(data) {
    y  <- data$Yield_z
    n  <- length(y)
    ZZ <- tcrossprod(stats::model.matrix(~ 0 + Batch, data))
    grid_sigma <- seq(0.01, 4, by = 0.01)
    grid_tau   <- seq(0, 4, by = 0.01)
    weights_sigma <- c(rep(0.01, length(grid_sigma) - 1L), 0.005)
    weights_tau   <- c(0.005, rep(0.01, length(grid_tau) - 2L), 0.005)
    log_prior <- function(x) log(2) + stats::dnorm(x, log = TRUE)
    evaluate <- function(sigma, tau) {
      R      <- chol(diag(sigma^2, n) + tau^2 * ZZ)
      solved <- backsolve(R, forwardsolve(t(R), cbind(y, 1)))
      y_V_y  <- sum(y * solved[, 1])
      one_V_y   <- sum(solved[, 1])
      one_V_one <- sum(solved[, 2])
      log_det <- 2 * sum(log(diag(R))) + log1p(one_V_one)
      c(log_lik = -0.5 * (n * log(2 * pi) + log_det + y_V_y - one_V_y^2 / (1 + one_V_one)),
        mu      = one_V_y / (1 + one_V_one))
    }
    values  <- array(NA_real_, c(2L, length(grid_sigma), length(grid_tau)))
    for (j in seq_along(grid_tau)) {
      values[, , j] <- vapply(grid_sigma, evaluate, numeric(2), tau = grid_tau[j])
    }
    log_lik  <- values[1L, , ]
    log_post <- log_lik + outer(log_prior(grid_sigma), log_prior(grid_tau), "+")
    scale    <- max(log_post)
    mass     <- exp(log_post - scale) * outer(weights_sigma, weights_tau)
    total    <- sum(mass)
    tau_density <- colSums(exp(log_post - scale) * weights_sigma) / total
    tau_cdf     <- c(0, cumsum((tau_density[-1L] + tau_density[-length(tau_density)]) / 2 * 0.01))
    log_null    <- log_lik[, 1L] + log_prior(grid_sigma)
    list(
      logml_h1 = scale + log(total),
      logml_h0 = max(log_null) + log(sum(exp(log_null - max(log_null)) * weights_sigma)),
      mean = c(intercept = sum(mass * values[2L, , ]) / total, sigma = sum(mass * grid_sigma) / total,
               "sd(Batch)" = sum(t(mass) * grid_tau) / total),
      tau = grid_tau, tau_density = tau_density, tau_cdf = tau_cdf
    )
  }
  exact <- list(
    Dyestuff  = scenario_fit("exact_dyestuff",  exact_reference(datasets[["Dyestuff"]])),
    Dyestuff2 = scenario_fit("exact_dyestuff2", exact_reference(datasets[["Dyestuff2"]]))
  )

  ### BayesTools fits ----
  # the random-intercept model (noncentered and mean-centered), the model
  # without batch effects, and a spike-and-slab prior on the batch SD
  fit_random <- function(name, data, parameterization) {
    scenario_fit(name, {
      unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
      fit <- JAGS_fit(
        model_syntax = likelihood, data = data_list(data), prior_list = list(sigma = unit_sd_prior),
        formula_list = list(mu = ~ 1 + (1 | Batch)), formula_data_list = list(mu = data),
        formula_prior_list = list(mu = list(intercept = prior("normal", list(mean = 0, sd = 1)))),
        formula_random_prior_list = list(mu = prior_random(Batch = random_block(
          sd = unit_sd_prior, parameterization = parameterization))),
        chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
      )
      list(fit = fit, marglik = JAGS_bridgesampling(
        fit, log_posterior = log_posterior, data = data_list(data), prior_list = list(sigma = unit_sd_prior), seed = 1))
    })
  }
  fit_null <- function(name, data) {
    scenario_fit(name, {
      unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
      fit <- JAGS_fit(
        model_syntax = likelihood, data = data_list(data), prior_list = list(sigma = unit_sd_prior),
        formula_list = list(mu = ~ 1), formula_data_list = list(mu = data),
        formula_prior_list = list(mu = list(intercept = prior("normal", list(mean = 0, sd = 1)))),
        chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
      )
      list(fit = fit, marglik = JAGS_bridgesampling(
        fit, log_posterior = log_posterior, data = data_list(data), prior_list = list(sigma = unit_sd_prior), seed = 1))
    })
  }
  fit_spike <- function(name, data) {
    scenario_fit(name, {
      unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
      JAGS_fit(
        model_syntax = likelihood, data = data_list(data), prior_list = list(sigma = unit_sd_prior),
        formula_list = list(mu = ~ 1 + (1 | Batch)), formula_data_list = list(mu = data),
        formula_prior_list = list(mu = list(intercept = prior("normal", list(mean = 0, sd = 1)))),
        formula_random_prior_list = list(mu = prior_random(Batch = random_block(
          sd = prior_spike_and_slab(prior_parameter = unit_sd_prior)))),
        chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
      )
    })
  }
  noncentered <- list(
    Dyestuff  = fit_random("fit_dyestuff_noncentered",  datasets[["Dyestuff"]],  "noncentered"),
    Dyestuff2 = fit_random("fit_dyestuff2_noncentered", datasets[["Dyestuff2"]], "noncentered")
  )
  mean_centered <- list(
    Dyestuff  = fit_random("fit_dyestuff_mean_centered",  datasets[["Dyestuff"]],  "mean_centered"),
    Dyestuff2 = fit_random("fit_dyestuff2_mean_centered", datasets[["Dyestuff2"]], "mean_centered")
  )
  null <- list(
    Dyestuff  = fit_null("fit_dyestuff_null",  datasets[["Dyestuff"]]),
    Dyestuff2 = fit_null("fit_dyestuff2_null", datasets[["Dyestuff2"]])
  )
  spike <- list(
    Dyestuff  = fit_spike("fit_dyestuff_spike",  datasets[["Dyestuff"]]),
    Dyestuff2 = fit_spike("fit_dyestuff2_spike", datasets[["Dyestuff2"]])
  )

  ### brms and lme4 fits ----
  # brms with the same priors (it bounds the SD and sigma priors at zero), with
  # and without batch effects; lme4 maximum likelihood
  fit_brms <- function(name, data, random) {
    scenario_fit(name, {
      formula <- if (random) Yield_z ~ 0 + Intercept + (1 | Batch) else Yield_z ~ 0 + Intercept
      priors  <- c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sigma),
                   if (random) brms::prior(normal(0, 1), class = sd))
      fit <- brms::brm(
        formula, data = data, prior = priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
        chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
        control = list(adapt_delta = 0.99)
      )
      set.seed(1)
      list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
    })
  }
  brms_random <- list(
    Dyestuff  = fit_brms("brms_dyestuff",  datasets[["Dyestuff"]],  random = TRUE),
    Dyestuff2 = fit_brms("brms_dyestuff2", datasets[["Dyestuff2"]], random = TRUE)
  )
  brms_null <- list(
    Dyestuff  = fit_brms("brms_dyestuff_null",  datasets[["Dyestuff"]],  random = FALSE),
    Dyestuff2 = fit_brms("brms_dyestuff2_null", datasets[["Dyestuff2"]], random = FALSE)
  )
  lme4_fits <- lapply(datasets, function(data) {
    suppressMessages(lme4::lmer(Yield_z ~ 1 + (1 | Batch), data = data, REML = FALSE))
  })

  ### model summaries ----
  scenario_text("estimates-dyestuff",  JAGS_estimates_table(noncentered[["Dyestuff"]][["fit"]]))
  scenario_text("estimates-dyestuff2", JAGS_estimates_table(noncentered[["Dyestuff2"]][["fit"]]))

  ### exact, lme4, and brms comparisons ----
  # the exact rows are posterior means; lme4 reports ML estimates (Dyestuff2:
  # a singular fit with batch SD 0)
  comparison_references <- function(name) {
    lme4_fit <- lme4_fits[[name]]
    list(
      exact = exact[[name]][["mean"]],
      lme4  = c(intercept = lme4::fixef(lme4_fit)[["(Intercept)"]], sigma = stats::sigma(lme4_fit),
                "sd(Batch)" = attr(lme4::VarCorr(lme4_fit)[["Batch"]], "stddev")[["(Intercept)"]])
    )
  }
  comparison_posteriors <- function(name) {
    list(
      BayesTools = scenario_draws(JAGS_estimates_table(noncentered[[name]][["fit"]], return_samples = TRUE),
                                  c(intercept = "(mu) intercept", sigma = "sigma", "sd(Batch)" = "(mu) sd(intercept)")),
      brms       = scenario_draws(as.matrix(brms_random[[name]][["fit"]]),
                                  c(intercept = "b_Intercept", sigma = "sigma", "sd(Batch)" = "sd_Batch__Intercept"))
    )
  }
  scenario_comparison(
    "comparison-dyestuff", c("intercept", "sigma", "sd(Batch)"), main = "Dyestuff",
    references = comparison_references("Dyestuff"), posteriors = comparison_posteriors("Dyestuff")
  )
  scenario_comparison(
    "comparison-dyestuff2", c("intercept", "sigma", "sd(Batch)"), main = "Dyestuff2",
    references = comparison_references("Dyestuff2"), posteriors = comparison_posteriors("Dyestuff2")
  )

  ### Bayes factors ----
  # the batch SD tau against no batch effects as log BF10: the exact reference,
  # bridge sampling in BayesTools and brms, the spike-and-slab inclusion BF,
  # and Savage-Dickey (the half-normal prior density at 0 against the
  # posterior density of tau at 0) under both parameterizations
  savage_dickey <- function(fit) {
    log(as.numeric(hypothesis_BF(fit, "sd(intercept) = 0", seed = 1)[["BF"]]))
  }
  compare_bayes_factors <- function(name) {
    scenario_bf_table(
      H1 = "(1 | Batch)",
      H0 = "no batch effects",
      log_BF10 = c(
        "exact"                        = exact[[name]][["logml_h1"]] - exact[[name]][["logml_h0"]],
        "bridge sampling"              = scenario_log_bf_bridge(noncentered[[name]][["marglik"]], null[[name]][["marglik"]]),
        "bridge sampling (brms)"       = scenario_log_bf_bridge(brms_random[[name]][["marglik"]], brms_null[[name]][["marglik"]]),
        "spike and slab"               = log(as.numeric(JAGS_inference_table(spike[[name]])[["inclusion_BF"]])),
        "Savage-Dickey (noncentered)"  = savage_dickey(noncentered[[name]][["fit"]]),
        "Savage-Dickey (mean-centered)" = savage_dickey(mean_centered[[name]][["fit"]])
      )
    )
  }
  scenario_text("bayes-factor-dyestuff",  compare_bayes_factors("Dyestuff"))
  scenario_text("bayes-factor-dyestuff2", compare_bayes_factors("Dyestuff2"))
  scenario_text("inference-spike-dyestuff",  JAGS_inference_table(spike[["Dyestuff"]]))
  scenario_text("inference-spike-dyestuff2", JAGS_inference_table(spike[["Dyestuff2"]]))

  ### MCMC diagnostics near tau = 0 ----
  # the table's ESS and R-hat of tau describe its bulk and look fine under both
  # parameterizations; the Savage-Dickey BF depends on the draws near 0, which
  # the ESS of the indicator tau < 0.05 and its spread across chains describe
  # (pooled draws are chain-major). The exact probability is the reference
  tau_draws <- function(fit) {
    tau <- JAGS_estimates_table(fit, return_samples = TRUE)[, "(mu) sd(intercept)"]
    matrix(tau, ncol = 4L)
  }
  scenario_text("diagnostics-tau", {
    rows <- list()
    for (name in names(datasets)) {
      for (parameterization in c("noncentered", "mean-centered")) {
        fit   <- if (parameterization == "noncentered") noncentered[[name]][["fit"]] else mean_centered[[name]][["fit"]]
        table <- JAGS_estimates_table(fit)
        tau   <- tau_draws(fit)
        rows[[length(rows) + 1L]] <- data.frame(
          data = name, parameterization = parameterization,
          ESS = round(table["(mu) sd(intercept)", "ESS"]), R_hat = round(table["(mu) sd(intercept)", "R_hat"], 3),
          ESS_tail_log_tau = round(posterior::ess_tail(log(tau))),
          ESS_tau_below_0.05 = round(posterior::ess_basic(1 * (tau < 0.05))),
          P_tau_below_0.05 = round(mean(tau < 0.05), 4),
          chain_min = round(min(colMeans(tau < 0.05)), 4), chain_max = round(max(colMeans(tau < 0.05)), 4),
          exact = round(stats::approx(exact[[name]][["tau"]], exact[[name]][["tau_cdf"]], 0.05)$y, 4)
        )
      }
    }
    do.call(rbind, rows)
  })

  ### posterior distribution of tau ----
  # top: the exact marginal posterior density of tau with the noncentered
  # draws; bottom: the posterior CDF near 0, where the mean-centered draws fall
  # short of the exact probability
  scenario_plot("posterior-tau", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (name in names(datasets)) {
      reference <- exact[[name]]
      graphics::hist(tau_draws(noncentered[[name]][["fit"]]), breaks = seq(0, 6, by = 0.05), freq = FALSE,
                     xlim = c(0, 2.5), col = "grey85", border = "white", main = name, xlab = "sd(Batch)")
      graphics::lines(reference[["tau"]], reference[["tau_density"]], lwd = 2)
    }
    for (name in names(datasets)) {
      reference <- exact[[name]]
      near_zero <- reference[["tau"]] <= 0.3
      graphics::plot(reference[["tau"]][near_zero], reference[["tau_cdf"]][near_zero], type = "l", lwd = 2,
                     main = paste0(name, ": CDF near 0"), xlab = "sd(Batch)", ylab = "P(sd(Batch) < x)")
      x <- seq(0, 0.3, by = 0.005)
      graphics::lines(x, stats::ecdf(tau_draws(noncentered[[name]][["fit"]]))(x),   col = "steelblue", lty = 2)
      graphics::lines(x, stats::ecdf(tau_draws(mean_centered[[name]][["fit"]]))(x), col = "firebrick", lty = 3)
      graphics::legend("topleft", c("exact", "noncentered", "mean-centered"), bty = "n",
                       col = c("black", "steelblue", "firebrick"), lty = 1:3, lwd = c(2, 1, 1))
    }
  })

  ### spike-and-slab posterior ----
  scenario_plot("posterior-spike", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    for (name in names(datasets)) {
      fit       <- spike[[name]]
      selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(intercept)")
      plot_posterior(stats::setNames(list(parameter_mixed_posterior(fit, selection)), "(mu) sd(intercept)"),
                     "(mu) sd(intercept)", prior = TRUE, main = name)
    }
  })
})
