if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("grouseticks")
# testthat::test_file("tests/scenarios/test-grouseticks.R")

### Description
# count GLMMs for the lme4 grouseticks data (ticks on red grouse chicks in 118
# broods over three years, with the altitude of the brood site): Poisson and
# negative-binomial (Poisson-gamma) models with random brood intercepts, a
# mean-difference year factor, and a standardized height (formula_scale_list).
# Compared with glmer (adaptive Gauss-Hermite quadrature) and glmer.nb, with
# brms under the same priors, and with reference marginal likelihoods that
# integrate each brood's intercept by adaptive Gauss-Hermite quadrature and
# the remaining parameters by importance sampling. Bayes factors for the year
# effect and for overdispersion (negative binomial against Poisson), and
# fixed-effect predictions on the log scale. Broods are nested in years, so
# the year coefficients and the brood intercepts are strongly correlated a
# posteriori and mix slowly in every brood parameterization; the centered
# parameterization with long chains is used (TODO: reparameterize the brood
# intercepts around their year means).

testthat::test_that("grouseticks Poisson and negative-binomial brood GLMMs", {

  data("grouseticks", package = "lme4")
  gt <- grouseticks
  gt$HEIGHT_s <- (gt$HEIGHT - mean(gt$HEIGHT)) / stats::sd(gt$HEIGHT)
  year_contrast <- contr.meandif(levels(gt$YEAR))
  gt$year_1 <- year_contrast[as.integer(gt$YEAR), 1L]
  gt$year_2 <- year_contrast[as.integer(gt$YEAR), 2L]

  poisson_likelihood <- "model{
    for(i in 1:N){
      TICKS[i] ~ dpois(exp(mu[i]))
    }
  }
  "
  nb_likelihood <- "model{
    for(i in 1:N){
      TICKS[i] ~ dnegbin(phi / (phi + exp(mu[i])), phi)
    }
  }
  "
  data_list <- list(TICKS = gt$TICKS, N = nrow(gt))
  poisson_log_posterior <- function(parameters, data) {
    sum(stats::dpois(data[["TICKS"]], exp(parameters[["mu"]]), log = TRUE))
  }
  nb_log_posterior <- function(parameters, data) {
    sum(stats::dnbinom(data[["TICKS"]], size = parameters[["phi"]], mu = exp(parameters[["mu"]]), log = TRUE))
  }
  # priors on the log scale: intercept ~ N(0, 2.5), height and the two year
  # mean-difference coefficients ~ N(0, 1), brood SD ~ half-normal(0, 1), and
  # the negative-binomial size phi ~ Gamma(2, 0.5)
  intercept_prior <- prior("normal", list(mean = 0, sd = 2.5))
  unit_prior      <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior   <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  year_prior      <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
  phi_prior       <- prior("gamma", list(shape = 2, rate = 0.5))

  ### glmer fits ----
  glmer_poisson <- lme4::glmer(TICKS ~ HEIGHT_s + year_1 + year_2 + (1 | BROOD), data = gt, family = poisson, nAGQ = 25)
  glmer_nb      <- lme4::glmer.nb(TICKS ~ HEIGHT_s + year_1 + year_2 + (1 | BROOD), data = gt)

  ### BayesTools fits ----
  fits <- list()
  fits[["poisson"]] <- scenario_fit("fit_poisson", {
    fit <- JAGS_fit(poisson_likelihood, data_list,
                    formula_list = list(mu = ~ HEIGHT + YEAR + (1 | BROOD)), formula_data_list = list(mu = gt),
                    formula_scale_list = list(mu = list(HEIGHT = TRUE)),
                    formula_prior_list = list(mu = list(intercept = intercept_prior, HEIGHT = unit_prior, YEAR = year_prior)),
                    formula_random_prior_list = list(mu = prior_random(BROOD = random_block(sd = unit_sd_prior, parameterization = "centered"))),
                    chains = 4, adapt = 1000, burnin = 5000, sample = 20000, seed = 1)
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = poisson_log_posterior, data = data_list, seed = 1))
  })
  fits[["poisson_no_year"]] <- scenario_fit("fit_poisson_no_year", {
    fit <- JAGS_fit(poisson_likelihood, data_list,
                    formula_list = list(mu = ~ HEIGHT + (1 | BROOD)), formula_data_list = list(mu = gt),
                    formula_scale_list = list(mu = list(HEIGHT = TRUE)),
                    formula_prior_list = list(mu = list(intercept = intercept_prior, HEIGHT = unit_prior)),
                    formula_random_prior_list = list(mu = prior_random(BROOD = random_block(sd = unit_sd_prior, parameterization = "centered"))),
                    chains = 4, adapt = 1000, burnin = 5000, sample = 20000, seed = 1)
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = poisson_log_posterior, data = data_list, seed = 1))
  })
  fits[["nb"]] <- scenario_fit("fit_nb", {
    fit <- JAGS_fit(nb_likelihood, data_list, prior_list = list(phi = phi_prior),
                    formula_list = list(mu = ~ HEIGHT + YEAR + (1 | BROOD)), formula_data_list = list(mu = gt),
                    formula_scale_list = list(mu = list(HEIGHT = TRUE)),
                    formula_prior_list = list(mu = list(intercept = intercept_prior, HEIGHT = unit_prior, YEAR = year_prior)),
                    formula_random_prior_list = list(mu = prior_random(BROOD = random_block(sd = unit_sd_prior, parameterization = "centered"))),
                    chains = 4, adapt = 1000, burnin = 5000, sample = 20000, seed = 1)
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = nb_log_posterior, data = data_list,
                                                  prior_list = list(phi = phi_prior), seed = 1))
  })
  fits[["nb_no_year"]] <- scenario_fit("fit_nb_no_year", {
    fit <- JAGS_fit(nb_likelihood, data_list, prior_list = list(phi = phi_prior),
                    formula_list = list(mu = ~ HEIGHT + (1 | BROOD)), formula_data_list = list(mu = gt),
                    formula_scale_list = list(mu = list(HEIGHT = TRUE)),
                    formula_prior_list = list(mu = list(intercept = intercept_prior, HEIGHT = unit_prior)),
                    formula_random_prior_list = list(mu = prior_random(BROOD = random_block(sd = unit_sd_prior, parameterization = "centered"))),
                    chains = 4, adapt = 1000, burnin = 5000, sample = 20000, seed = 1)
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = nb_log_posterior, data = data_list,
                                                  prior_list = list(phi = phi_prior), seed = 1))
  })

  ### brms fits ----
  brms_priors <- c(brms::prior(normal(0, 2.5), class = b, coef = Intercept), brms::prior(normal(0, 1), class = b),
                   brms::prior(normal(0, 1), class = sd))
  brms_marglik <- function(fit) {
    set.seed(1)
    brms::bridge_sampler(fit, silent = TRUE)
  }
  brms_fits <- list()
  brms_fits[["poisson"]] <- scenario_fit("brms_poisson", {
    fit <- brms::brm(TICKS ~ 0 + Intercept + HEIGHT_s + year_1 + year_2 + (1 | BROOD), data = gt, family = poisson(),
                     prior = brms_priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["poisson_no_year"]] <- scenario_fit("brms_poisson_no_year", {
    fit <- brms::brm(TICKS ~ 0 + Intercept + HEIGHT_s + (1 | BROOD), data = gt, family = poisson(),
                     prior = brms_priors, save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["nb"]] <- scenario_fit("brms_nb", {
    fit <- brms::brm(TICKS ~ 0 + Intercept + HEIGHT_s + year_1 + year_2 + (1 | BROOD), data = gt, family = brms::negbinomial(),
                     prior = c(brms_priors, brms::prior(gamma(2, 0.5), class = shape)),
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["nb_no_year"]] <- scenario_fit("brms_nb_no_year", {
    fit <- brms::brm(TICKS ~ 0 + Intercept + HEIGHT_s + (1 | BROOD), data = gt, family = brms::negbinomial(),
                     prior = c(brms_priors, brms::prior(gamma(2, 0.5), class = shape)),
                     save_pars = brms::save_pars(all = TRUE), backend = "rstan",
                     chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0)
    list(fit = fit, marglik = brms_marglik(fit))
  })

  ### reference marginal likelihoods ----
  # log p(y | beta, tau[, phi]) integrates each brood's intercept u by 20-node
  # Gauss-Hermite quadrature centred at its conditional mode, found by damped
  # Newton steps with the analytic derivatives of the concave log integrand,
  # and scaled by its curvature there; beta, log tau, and log phi are
  # importance-sampled from a multivariate t (5 df) fitted to the BayesTools
  # posterior draws, with the priors above (the half-normal and gamma
  # densities transformed to the log scale)
  reference <- scenario_fit("reference_logml", {
    reference_logml <- function(fit, with_year, nb, size = 10000L) {
      brood <- as.integer(factor(gt$BROOD))
      y <- gt$TICKS
      X <- cbind(1, gt$HEIGHT_s, if (with_year) cbind(gt$year_1, gt$year_2))
      nodes <- statmod::gauss.quad(20L, kind = "hermite")
      log_brood <- function(beta, tau, phi) {
        eta <- drop(X %*% beta)
        # the log integrand at one value of u per brood and node (columns), and
        # its first two derivatives at one value per brood
        log_f <- function(u) {
          lambda <- exp(eta + u[brood, , drop = FALSE])
          log_p <- if (nb) stats::dnbinom(y, size = phi, mu = lambda, log = TRUE) else stats::dpois(y, lambda, log = TRUE)
          rowsum(log_p, brood, reorder = TRUE) + stats::dnorm(u, 0, tau, log = TRUE)
        }
        derivatives <- function(u) {
          lambda <- exp(eta + u[brood])
          d1 <- if (nb) y - (y + phi) * lambda / (phi + lambda) else y - lambda
          d2 <- if (nb) -(y + phi) * phi * lambda / (phi + lambda)^2 else -lambda
          cbind(rowsum(d1, brood, reorder = TRUE) - u / tau^2, rowsum(d2, brood, reorder = TRUE) - 1 / tau^2)
        }
        u <- numeric(max(brood))
        for (iteration in 1:50) {
          d    <- derivatives(u)
          step <- pmax(pmin(-d[, 1L] / d[, 2L], 1), -1)
          u    <- u + step
          if (max(abs(step)) < 1e-10) break
        }
        s <- sqrt(-1 / derivatives(u)[, 2L])
        values <- log_f(u + sqrt(2) * s %*% t(nodes$nodes)) +
          matrix(log(nodes$weights) + nodes$nodes^2, max(brood), 20L, byrow = TRUE)
        top <- apply(values, 1L, max)
        sum(top + log(rowSums(exp(values - top))) + log(sqrt(2) * s))
      }
      draws <- JAGS_estimates_table(fit, return_samples = TRUE)
      U <- cbind(draws[, c("(mu) intercept", "(mu) HEIGHT", if (with_year) c("(mu) YEAR{1}", "(mu) YEAR{2}"))],
                 log(draws[, "(mu) sd(intercept)"]), if (nb) log(draws[, "phi"]))
      k <- ncol(X)
      set.seed(1)
      center <- colMeans(U)
      scale  <- 1.5 * stats::cov(U)
      proposals <- mvtnorm::rmvt(size, sigma = scale, df = 5, delta = center, type = "shifted")
      log_q     <- mvtnorm::dmvt(proposals, delta = center, sigma = scale, df = 5, log = TRUE, type = "shifted")
      log_w <- vapply(seq_len(size), function(i) {
        p   <- proposals[i, ]
        tau <- exp(p[[k + 1L]])
        phi <- if (nb) exp(p[[k + 2L]]) else NULL
        log_brood(p[seq_len(k)], tau, phi) +
          stats::dnorm(p[[1L]], 0, 2.5, log = TRUE) + sum(stats::dnorm(p[2:k], 0, 1, log = TRUE)) +
          log(2) + stats::dnorm(tau, 0, 1, log = TRUE) + p[[k + 1L]] +
          if (nb) stats::dgamma(phi, 2, 0.5, log = TRUE) + p[[k + 2L]] else 0
      }, numeric(1)) - log_q
      if (anyNA(log_w)) {
        stop("The reference importance weights are not all defined.", call. = FALSE)
      }
      weights <- exp(log_w - max(log_w))
      c(logml = max(log_w) + log(mean(weights)), relative_se = stats::sd(weights) / sqrt(size) / mean(weights))
    }
    list(
      poisson         = reference_logml(fits[["poisson"]][["fit"]], with_year = TRUE, nb = FALSE),
      poisson_no_year = reference_logml(fits[["poisson_no_year"]][["fit"]], with_year = FALSE, nb = FALSE),
      nb              = reference_logml(fits[["nb"]][["fit"]], with_year = TRUE, nb = TRUE),
      nb_no_year      = reference_logml(fits[["nb_no_year"]][["fit"]], with_year = FALSE, nb = TRUE)
    )
  })

  ### outputs ----
  scenario_text("estimates-poisson", JAGS_estimates_table(fits[["poisson"]][["fit"]]))
  scenario_text("estimates-nb",      JAGS_estimates_table(fits[["nb"]][["fit"]]))

  scenario_text("comparison-nb", {
    draws <- JAGS_estimates_table(fits[["nb"]][["fit"]], return_samples = TRUE)
    brms_draws <- as.matrix(brms_fits[["nb"]][["fit"]])
    parameters <- c("intercept", "height", "year 1", "year 2", "brood SD", "phi")
    glmer_values <- c(lme4::fixef(glmer_nb), attr(lme4::VarCorr(glmer_nb)[["BROOD"]], "stddev"), lme4::getME(glmer_nb, "glmer.nb.theta"))
    scenario_comparison_table(
      parameters,
      references = list(glmer.nb = stats::setNames(unname(glmer_values), parameters)),
      posteriors = list(
        BayesTools = scenario_draws(draws, stats::setNames(c("(mu) intercept", "(mu) HEIGHT", "(mu) YEAR{1}", "(mu) YEAR{2}",
                                                             "(mu) sd(intercept)", "phi"), parameters)),
        brms       = scenario_draws(brms_draws, stats::setNames(c("b_Intercept", "b_HEIGHT_s", "b_year_1", "b_year_2",
                                                                  "sd_BROOD__Intercept", "shape"), parameters))
      )
    )
  })
  scenario_text("comparison-poisson", {
    draws <- JAGS_estimates_table(fits[["poisson"]][["fit"]], return_samples = TRUE)
    brms_draws <- as.matrix(brms_fits[["poisson"]][["fit"]])
    parameters <- c("intercept", "height", "year 1", "year 2", "brood SD")
    glmer_values <- c(lme4::fixef(glmer_poisson), attr(lme4::VarCorr(glmer_poisson)[["BROOD"]], "stddev"))
    scenario_comparison_table(
      parameters,
      references = list("glmer (AGQ)" = stats::setNames(unname(glmer_values), parameters)),
      posteriors = list(
        BayesTools = scenario_draws(draws, stats::setNames(c("(mu) intercept", "(mu) HEIGHT", "(mu) YEAR{1}", "(mu) YEAR{2}",
                                                             "(mu) sd(intercept)"), parameters)),
        brms       = scenario_draws(brms_draws, stats::setNames(c("b_Intercept", "b_HEIGHT_s", "b_year_1", "b_year_2",
                                                                  "sd_BROOD__Intercept"), parameters))
      )
    )
  })

  scenario_text("bayes-factors", {
    log_bf <- function(fits, h1, h0) scenario_log_bf_bridge(fits[[h1]][["marglik"]], fits[[h0]][["marglik"]])
    reference_bf <- function(h1, h0) reference[[h1]][["logml"]] - reference[[h0]][["logml"]]
    reference_note <- function(h1, h0) {
      sprintf("reference relative SEs %.1f%%, %.1f%%", 100 * reference[[h1]][["relative_se"]], 100 * reference[[h0]][["relative_se"]])
    }
    rbind(
      scenario_bf_table("year (Poisson)", "no year", c(
        "importance sampling"    = reference_bf("poisson", "poisson_no_year"),
        "bridge sampling"        = log_bf(fits, "poisson", "poisson_no_year"),
        "bridge sampling (brms)" = log_bf(brms_fits, "poisson", "poisson_no_year")),
        note = reference_note("poisson", "poisson_no_year")),
      scenario_bf_table("year (negative binomial)", "no year", c(
        "importance sampling"    = reference_bf("nb", "nb_no_year"),
        "bridge sampling"        = log_bf(fits, "nb", "nb_no_year"),
        "bridge sampling (brms)" = log_bf(brms_fits, "nb", "nb_no_year")),
        note = reference_note("nb", "nb_no_year")),
      scenario_bf_table("negative binomial", "Poisson", c(
        "importance sampling"    = reference_bf("nb", "poisson"),
        "bridge sampling"        = log_bf(fits, "nb", "poisson"),
        "bridge sampling (brms)" = log_bf(brms_fits, "nb", "poisson")),
        note = reference_note("nb", "poisson"))
    )
  })

  scenario_text("predictions-nb", {
    grid <- expand.grid(YEAR = levels(gt$YEAR), HEIGHT = c(420, 480, 540))
    grid$YEAR     <- factor(grid$YEAR, levels = levels(gt$YEAR))
    grid$HEIGHT_s <- (grid$HEIGHT - mean(gt$HEIGHT)) / stats::sd(gt$HEIGHT)
    grid$year_1   <- year_contrast[as.integer(grid$YEAR), 1L]
    grid$year_2   <- year_contrast[as.integer(grid$YEAR), 2L]
    scenario_prediction_table(grid[, c("YEAR", "HEIGHT")], list(
      "glmer.nb"   = stats::predict(glmer_nb, newdata = grid, re.form = NA),
      "BayesTools" = JAGS_predict_formula(fits[["nb"]][["fit"]], "mu", data = grid[, c("YEAR", "HEIGHT")], formula_target = "fixed")$value,
      "brms"       = t(brms::posterior_linpred(brms_fits[["nb"]][["fit"]], newdata = grid, re_formula = NA))
    ))
  })

  scenario_plot("brood-intercepts", {
    # posterior mean brood intercepts of the negative-binomial model:
    # BayesTools against brms, and against the glmer.nb conditional modes
    samples <- JAGS_estimates_table(fits[["nb"]][["fit"]], return_samples = TRUE)
    brood_mean <- function(fit) {
      predictions <- JAGS_predict_formula(fit, "mu", data = gt[, c("YEAR", "HEIGHT", "BROOD")], formula_target = "conditional")$value
      fixed       <- JAGS_predict_formula(fit, "mu", data = gt[, c("YEAR", "HEIGHT", "BROOD")], formula_target = "fixed")$value
      tapply(rowMeans(predictions - fixed), gt$BROOD, mean)
    }
    bayestools <- brood_mean(fits[["nb"]][["fit"]])
    brms_ranef <- brms::ranef(brms_fits[["nb"]][["fit"]])[["BROOD"]][, "Estimate", "Intercept"]
    glmer_ranef <- lme4::ranef(glmer_nb)[["BROOD"]][, 1L]
    names(glmer_ranef) <- rownames(lme4::ranef(glmer_nb)[["BROOD"]])
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    scenario_agreement_plot(brms_ranef[names(bayestools)], bayestools, main = "brood intercepts", reference_label = "brms")
    scenario_agreement_plot(glmer_ranef[names(bayestools)], bayestools, main = "brood intercepts", reference_label = "glmer.nb")
  })
})
