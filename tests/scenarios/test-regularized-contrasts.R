if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized-contrasts")
source(file.path(.scenario_config()$root, "helpers", "regularized-display.R"), local = TRUE)
# testthat::test_file("tests/scenarios/test-regularized-contrasts.R")

### Monthly ozone in New York
# Compare shrinkage of May–September deviations and the wind slope under two
# factor codings. Ozone and wind are standardized, with residual SD fixed at 1.
testthat::test_that("Monthly ozone effects under centered shrinkage priors", {

  ozone <- na.omit(airquality[c("Ozone", "Month", "Wind")])
  ozone$Month <- factor(ozone$Month, labels = month.abb[5:9])
  ozone$Ozone <- as.numeric(scale(ozone$Ozone))
  ozone$Wind  <- as.numeric(scale(ozone$Wind))
  months <- levels(ozone$Month)
  data_list <- list(N = nrow(ozone), y = ozone$Ozone)

  # Construct Stan's two month designs independently of BayesTools.
  helmert <- contr.helmert(months)
  orthonormal <- sweep(helmert, 2, sqrt(colSums(helmert^2)), "/")
  meandif <- orthonormal / sqrt(1 - 1 / length(months))
  X_md <- cbind(meandif[as.integer(ozone$Month), ], Wind = ozone$Wind)
  X_on <- cbind(orthonormal[as.integer(ozone$Month), ], Wind = ozone$Wind)
  stan_md <- c(data_list, list(K = ncol(X_md), X = X_md))
  stan_on <- c(data_list, list(K = ncol(X_on), X = X_on))
  # Four month coefficients share one local horseshoe scale. Wind has its own.
  units <- list(J = 2L, unit = c(1L, 1L, 1L, 1L, 2L))
  horseshoe_md <- c(stan_md, units)
  horseshoe_on <- c(stan_on, units)

  likelihood <- "model {
    for (i in 1:N) { y[i] ~ dnorm(mu[i], 1) }
  }"
  log_likelihood <- function(parameters, data) {
    sum(dnorm(data$y, parameters$mu, sd = 1, log = TRUE))
  }
  intercept <- prior("normal", list(0, 1))
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))
  fit_version <- regularized_fit_version(ozone, data_list, likelihood, deparse(body(log_likelihood)), intercept, global)

  ### Mean-difference coding / Ridge ----
  # The intercept has its own Normal prior. All slopes share the global scale.
  fit_meandif_ridge <- scenario_fit("fit_meandif_ridge", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("ridge", scale = global, contrast = "meandif"))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 53214280)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_meandif_ridge <- fit_meandif_ridge$fit
  # First inspect the fitted coefficients, then express that same fit as monthly deviations.
  scenario_text("estimates-meandif_ridge-coded", JAGS_estimates_table(jags_meandif_ridge))
  scenario_text("estimates-meandif_ridge", JAGS_estimates_table(jags_meandif_ridge, transform_factors = TRUE))

  ref_meandif_ridge <- scenario_fit("ref_meandif_ridge", {
    fit <- rstan::sampling(scenario_stan_model("ridge.stan"), data = stan_md,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(53214281)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(stan_md, scenario_stan_version("ridge.stan")))
  stan_meandif_ridge <- ref_meandif_ridge$fit
  scenario_text("stan-estimates-meandif_ridge", round(rstan::summary(stan_meandif_ridge, pars = c("alpha", "beta", "tau"))$summary, 3))

  # BayesTools reports the level deviations directly. Stan's first four beta
  # draws are coded coefficients; multiplying by the contrast matrix returns
  # the five deviations from the equally weighted mean at average wind.
  effects <- c(months, "Wind")
  jags_names <- setNames(c(paste0("(mu) Month[dif: ", months, "]"), "(mu) Wind"), effects)
  stan_month_draws <- function(fit, contrasts) {
    beta <- rstan::extract(fit, pars = "beta")$beta
    deviations <- beta[, 1:4] %*% t(contrasts)
    colnames(deviations) <- months
    cbind(deviations, Wind = beta[, 5])
  }

  draws_meandif_ridge <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_meandif_ridge, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_meandif_ridge, meandif))
  scenario_comparison("comparison-meandif_ridge", effects, posteriors = draws_meandif_ridge,
                      main = "Ridge with mean-difference month contrasts", show_table = TRUE)

  # Inspect one coefficient in the native BayesTools plot and add the Stan density.
  scenario_plot("ridge-wind-posterior", {
    selection <- parameter_catalog_resolve(parameter_catalog(jags_meandif_ridge), "(mu) Wind")
    posterior <- parameter_mixed_posterior(jags_meandif_ridge, selection)
    plot_posterior(list(Wind = posterior), "Wind", prior = TRUE,
                   main = "Ozone and wind at a fixed month", xlim = c(-1.5, 1.5),
                   xlab = "Change in standardized ozone per SD of wind",
                   col = "#0072B2", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_meandif_ridge$RStan[, "Wind"]), col = "#D55E00", lty = 2, lwd = 2)
    legend("topright", c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
           col = c("grey70", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)
  })

  ### Mean-difference coding / Lasso ----
  # The month contrasts form one group-lasso term. Wind has a scalar Laplace prior.
  fit_meandif_lasso <- scenario_fit("fit_meandif_lasso", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("lasso", scale = global, contrast = "meandif"))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 55271485)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_meandif_lasso <- fit_meandif_lasso$fit
  scenario_text("estimates-meandif_lasso", JAGS_estimates_table(jags_meandif_lasso, transform_factors = TRUE))

  ref_meandif_lasso <- scenario_fit("ref_meandif_lasso", {
    fit <- rstan::sampling(scenario_stan_model("grouped-lasso.stan"), data = stan_md,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(55271486)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(stan_md, scenario_stan_version("grouped-lasso.stan")))
  stan_meandif_lasso <- ref_meandif_lasso$fit
  scenario_text("stan-estimates-meandif_lasso", round(rstan::summary(stan_meandif_lasso, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_meandif_lasso <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_meandif_lasso, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_meandif_lasso, meandif))
  scenario_comparison("comparison-meandif_lasso", effects, posteriors = draws_meandif_lasso,
                      main = "Group lasso with mean-difference month contrasts", show_table = TRUE)

  scenario_plot("meandif-lasso-july-posterior", {
    selection <- parameter_catalog_resolve(parameter_catalog(jags_meandif_lasso), "mu_Month[Jul]")
    posterior <- parameter_mixed_posterior(jags_meandif_lasso, selection)
    plot_posterior(list(Jul = posterior), "Jul", prior = TRUE,
                   main = "Group lasso for July ozone (mean-difference coding)", xlim = c(-1.5, 1.5),
                   xlab = "July deviation from the five-month mean (standardized ozone)",
                   col = "#0072B2", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_meandif_lasso$RStan[, "Jul"]), col = "#D55E00", lty = 2, lwd = 2)
    legend("topright", c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
           col = c("grey70", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)
  })

  ### Mean-difference coding / Horseshoe ----
  fit_meandif_horseshoe <- scenario_fit("fit_meandif_horseshoe", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("horseshoe", scale = global, contrast = "meandif"))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 148854667)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_meandif_horseshoe <- fit_meandif_horseshoe$fit
  scenario_text("estimates-meandif_horseshoe", JAGS_estimates_table(jags_meandif_horseshoe, transform_factors = TRUE))

  ref_meandif_horseshoe <- scenario_fit("ref_meandif_horseshoe", {
    fit <- rstan::sampling(scenario_stan_model("horseshoe.stan"), data = horseshoe_md,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(148854668)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(horseshoe_md, scenario_stan_version("horseshoe.stan")))
  stan_meandif_horseshoe <- ref_meandif_horseshoe$fit
  scenario_text("stan-estimates-meandif_horseshoe", round(rstan::summary(stan_meandif_horseshoe, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_meandif_horseshoe <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_meandif_horseshoe, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_meandif_horseshoe, meandif))
  scenario_comparison("comparison-meandif_horseshoe", effects, posteriors = draws_meandif_horseshoe,
                      main = "Horseshoe with mean-difference month contrasts", show_table = TRUE)

  # An additional view of the same fit overlays the five levels of Month.
  month_colours <- setNames(c("#0072B2", "#D55E00", "#009E73", "#A44A87", "#66552E"), months)
  scenario_plot("month-overlay-meandif-horseshoe", {
    posterior <- as_mixed_posteriors(jags_meandif_horseshoe, parameters = c("mu_Month", "mu_Wind"))
    old_par <- par(no.readonly = TRUE)
    layout(matrix(1:2, nrow = 1), widths = c(1.6, 1))
    par(mar = c(4, 3.8, 2.3, 0.8), oma = c(2, 0, 4, 0))

    plot_posterior(posterior, "mu_Month", prior = TRUE, main = "Month",
                   xlab = "Deviation from the five-month mean", ylab = "Density", xlim = c(-1.5, 1.5),
                   col = unname(month_colours), lwd = 2, dots_prior = list(col = "grey70", lwd = 2),
                   legend_labels = months, legend_position = "topright")
    for (month in months) {
      lines(density(draws_meandif_horseshoe$RStan[, month]), col = month_colours[[month]], lty = 2, lwd = 2)
    }

    plot_posterior(posterior, "mu_Wind", prior = TRUE, main = "Wind",
                   xlab = "Slope per SD of wind", ylab = "Density", xlim = c(-1.5, 1.5),
                   col = "#34495E", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_meandif_horseshoe$RStan[, "Wind"]), col = "#34495E", lty = 2, lwd = 2)
    mtext("Horseshoe with mean-difference month contrasts", side = 3, outer = TRUE, line = 1.8, font = 2)
    mtext("Grey: prior     Solid: BayesTools posterior     Dashed: RStan posterior",
          side = 3, outer = TRUE, line = 0.35, cex = 0.85)
    mtext("Effects are in standard deviations of ozone. The prior densities are unbounded at zero.",
          side = 1, outer = TRUE, line = 0.35, cex = 0.8)
    par(old_par)
  })

  ### Monthly ozone at average wind ----
  # These predictions include the intercept. Compare expected standardized
  # ozone in each month at Wind = 0, the mean wind speed.
  month_grid <- data.frame(Month = factor(months, levels = months), Wind = 0)
  monthly_ozone <- JAGS_predict_formula(jags_meandif_horseshoe, "mu", data = month_grid, formula_target = "fixed")
  stan_coef <- rstan::extract(stan_meandif_horseshoe, pars = c("alpha", "beta"))
  stan_ozone <- sweep(meandif %*% t(stan_coef$beta[, 1:4]), 2, stan_coef$alpha, "+")
  scenario_text("ozone-at-average-wind", scenario_prediction_table(month_grid["Month"],
    list(BayesTools = monthly_ozone$value, RStan = stan_ozone)))

  # Is expected ozone higher in July than in May at the same wind speed?
  # Symmetry gives equal prior odds, so the RStan posterior odds are an
  # independent estimate of the directional Bayes factor.
  july_vs_may <- hypothesis_BF(jags_meandif_horseshoe, "Month[Jul] > Month[May] vs Month[Jul] < Month[May]", seed = 1)
  scenario_text("july-versus-may", july_vs_may)
  july_higher <- stan_ozone[months == "Jul", ] > stan_ozone[months == "May", ]
  scenario_text("july-versus-may-comparison", data.frame(
    BayesTools_BF = round(as.data.frame(july_vs_may)$BF10, 2),
    RStan_BF = round(sum(july_higher) / sum(!july_higher), 2)))

  ### Mean-difference coding / Regularized horseshoe (finite slab) ----
  # The Student-t slab has scale 2 and four degrees of freedom.
  fit_meandif_regularized_horseshoe <- scenario_fit("fit_meandif_regularized_horseshoe", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("horseshoe", scale = global, contrast = "meandif", slab_scale = 2, slab_df = 4))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 1086261)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_meandif_regularized_horseshoe <- fit_meandif_regularized_horseshoe$fit
  scenario_text("estimates-meandif_regularized_horseshoe", JAGS_estimates_table(jags_meandif_regularized_horseshoe, transform_factors = TRUE))

  ref_meandif_regularized_horseshoe <- scenario_fit("ref_meandif_regularized_horseshoe", {
    fit <- rstan::sampling(scenario_stan_model("finite-slab-horseshoe.stan"), data = horseshoe_md,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(1086262)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(horseshoe_md, scenario_stan_version("finite-slab-horseshoe.stan")))
  stan_meandif_regularized_horseshoe <- ref_meandif_regularized_horseshoe$fit
  scenario_text("stan-estimates-meandif_regularized_horseshoe", round(rstan::summary(stan_meandif_regularized_horseshoe, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_meandif_regularized_horseshoe <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_meandif_regularized_horseshoe, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_meandif_regularized_horseshoe, meandif))
  scenario_comparison("comparison-meandif_regularized_horseshoe", effects, posteriors = draws_meandif_regularized_horseshoe,
                      main = "Finite-slab horseshoe with mean-difference month contrasts", show_table = TRUE)

  scenario_plot("meandif-regularized-horseshoe-july-posterior", {
    selection <- parameter_catalog_resolve(parameter_catalog(jags_meandif_regularized_horseshoe), "mu_Month[Jul]")
    posterior <- parameter_mixed_posterior(jags_meandif_regularized_horseshoe, selection)
    plot_posterior(list(Jul = posterior), "Jul", prior = TRUE,
                   main = "Regularized horseshoe for July ozone (mean-difference coding)", xlim = c(-1.5, 1.5),
                   xlab = "July deviation from the five-month mean (standardized ozone)",
                   col = "#0072B2", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_meandif_regularized_horseshoe$RStan[, "Jul"]), col = "#D55E00", lty = 2, lwd = 2)
    legend("topright", c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
           col = c("grey70", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)
  })

  ### Mean-difference coding / Evidence for each prior relative to ridge ----
  scenario_text("bridge-meandif_ridge", fit_meandif_ridge$marglik)
  scenario_text("stan-bridge-meandif_ridge", ref_meandif_ridge$marglik)
  scenario_text("model-comparison-meandif", regularized_evidence_table(list(
    "Ridge" = list(BayesTools = fit_meandif_ridge$marglik, RStan = ref_meandif_ridge$marglik),
    "Lasso" = list(BayesTools = fit_meandif_lasso$marglik, RStan = ref_meandif_lasso$marglik),
    "Horseshoe" = list(BayesTools = fit_meandif_horseshoe$marglik, RStan = ref_meandif_horseshoe$marglik),
    "Finite-slab horseshoe" = list(BayesTools = fit_meandif_regularized_horseshoe$marglik, RStan = ref_meandif_regularized_horseshoe$marglik))))

  ### Orthonormal coding / Ridge ----
  # The intercept has its own Normal prior. All slopes share the global scale.
  fit_orthonormal_ridge <- scenario_fit("fit_orthonormal_ridge", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("ridge", scale = global, contrast = "orthonormal"))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 109793227)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_orthonormal_ridge <- fit_orthonormal_ridge$fit
  scenario_text("estimates-orthonormal_ridge-coded", JAGS_estimates_table(jags_orthonormal_ridge))
  scenario_text("estimates-orthonormal_ridge", JAGS_estimates_table(jags_orthonormal_ridge, transform_factors = TRUE))

  ref_orthonormal_ridge <- scenario_fit("ref_orthonormal_ridge", {
    fit <- rstan::sampling(scenario_stan_model("ridge.stan"), data = stan_on,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(109793228)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(stan_on, scenario_stan_version("ridge.stan")))
  stan_orthonormal_ridge <- ref_orthonormal_ridge$fit
  scenario_text("stan-estimates-orthonormal_ridge", round(rstan::summary(stan_orthonormal_ridge, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_orthonormal_ridge <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_orthonormal_ridge, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_orthonormal_ridge, orthonormal))
  scenario_comparison("comparison-orthonormal_ridge", effects, posteriors = draws_orthonormal_ridge,
                      main = "Ridge with orthonormal month contrasts", show_table = TRUE)

  ### Orthonormal coding / Lasso ----
  # The month contrasts form one group-lasso term. Wind has a scalar Laplace prior.
  fit_orthonormal_lasso <- scenario_fit("fit_orthonormal_lasso", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("lasso", scale = global, contrast = "orthonormal"))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 230491521)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_orthonormal_lasso <- fit_orthonormal_lasso$fit
  scenario_text("estimates-orthonormal_lasso", JAGS_estimates_table(jags_orthonormal_lasso, transform_factors = TRUE))

  ref_orthonormal_lasso <- scenario_fit("ref_orthonormal_lasso", {
    fit <- rstan::sampling(scenario_stan_model("grouped-lasso.stan"), data = stan_on,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(230491522)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(stan_on, scenario_stan_version("grouped-lasso.stan")))
  stan_orthonormal_lasso <- ref_orthonormal_lasso$fit
  scenario_text("stan-estimates-orthonormal_lasso", round(rstan::summary(stan_orthonormal_lasso, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_orthonormal_lasso <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_orthonormal_lasso, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_orthonormal_lasso, orthonormal))
  scenario_comparison("comparison-orthonormal_lasso", effects, posteriors = draws_orthonormal_lasso,
                      main = "Group lasso with orthonormal month contrasts", show_table = TRUE)

  scenario_plot("orthonormal-lasso-july-posterior", {
    selection <- parameter_catalog_resolve(parameter_catalog(jags_orthonormal_lasso), "mu_Month[Jul]")
    posterior <- parameter_mixed_posterior(jags_orthonormal_lasso, selection)
    plot_posterior(list(Jul = posterior), "Jul", prior = TRUE,
                   main = "Group lasso for July ozone (orthonormal coding)", xlim = c(-1.5, 1.5),
                   xlab = "July deviation from the five-month mean (standardized ozone)",
                   col = "#0072B2", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_orthonormal_lasso$RStan[, "Jul"]), col = "#D55E00", lty = 2, lwd = 2)
    legend("topright", c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
           col = c("grey70", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)
  })

  ### Orthonormal coding / Horseshoe ----
  fit_orthonormal_horseshoe <- scenario_fit("fit_orthonormal_horseshoe", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("horseshoe", scale = global, contrast = "orthonormal"))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 99242205)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_orthonormal_horseshoe <- fit_orthonormal_horseshoe$fit
  scenario_text("estimates-orthonormal_horseshoe", JAGS_estimates_table(jags_orthonormal_horseshoe, transform_factors = TRUE))

  ref_orthonormal_horseshoe <- scenario_fit("ref_orthonormal_horseshoe", {
    fit <- rstan::sampling(scenario_stan_model("horseshoe.stan"), data = horseshoe_on,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(99242206)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(horseshoe_on, scenario_stan_version("horseshoe.stan")))
  stan_orthonormal_horseshoe <- ref_orthonormal_horseshoe$fit
  scenario_text("stan-estimates-orthonormal_horseshoe", round(rstan::summary(stan_orthonormal_horseshoe, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_orthonormal_horseshoe <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_orthonormal_horseshoe, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_orthonormal_horseshoe, orthonormal))
  scenario_comparison("comparison-orthonormal_horseshoe", effects, posteriors = draws_orthonormal_horseshoe,
                      main = "Horseshoe with orthonormal month contrasts", show_table = TRUE)

  ### Orthonormal coding / Regularized horseshoe (finite slab) ----
  # The Student-t slab has scale 2 and four degrees of freedom.
  fit_orthonormal_regularized_horseshoe <- scenario_fit("fit_orthonormal_regularized_horseshoe", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = ~ Month + Wind), formula_data_list = list(mu = ozone),
      formula_prior_list = list(mu = list(intercept = intercept,
        sh = prior_regularized("horseshoe", scale = global, contrast = "orthonormal", slab_scale = 2, slab_df = 4))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 73429)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 190355105)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_orthonormal_regularized_horseshoe <- fit_orthonormal_regularized_horseshoe$fit
  scenario_text("estimates-orthonormal_regularized_horseshoe", JAGS_estimates_table(jags_orthonormal_regularized_horseshoe, transform_factors = TRUE))

  ref_orthonormal_regularized_horseshoe <- scenario_fit("ref_orthonormal_regularized_horseshoe", {
    fit <- rstan::sampling(scenario_stan_model("finite-slab-horseshoe.stan"), data = horseshoe_on,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 6283,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(190355106)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(horseshoe_on, scenario_stan_version("finite-slab-horseshoe.stan")))
  stan_orthonormal_regularized_horseshoe <- ref_orthonormal_regularized_horseshoe$fit
  scenario_text("stan-estimates-orthonormal_regularized_horseshoe", round(rstan::summary(stan_orthonormal_regularized_horseshoe, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_orthonormal_regularized_horseshoe <- list(
    BayesTools = scenario_draws(JAGS_estimates_table(jags_orthonormal_regularized_horseshoe, transform_factors = TRUE, return_samples = TRUE), jags_names),
    RStan = stan_month_draws(stan_orthonormal_regularized_horseshoe, orthonormal))
  scenario_comparison("comparison-orthonormal_regularized_horseshoe", effects, posteriors = draws_orthonormal_regularized_horseshoe,
                      main = "Finite-slab horseshoe with orthonormal month contrasts", show_table = TRUE)

  scenario_plot("orthonormal-regularized-horseshoe-july-posterior", {
    selection <- parameter_catalog_resolve(parameter_catalog(jags_orthonormal_regularized_horseshoe), "mu_Month[Jul]")
    posterior <- parameter_mixed_posterior(jags_orthonormal_regularized_horseshoe, selection)
    plot_posterior(list(Jul = posterior), "Jul", prior = TRUE,
                   main = "Regularized horseshoe for July ozone (orthonormal coding)", xlim = c(-1.5, 1.5),
                   xlab = "July deviation from the five-month mean (standardized ozone)",
                   col = "#0072B2", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_orthonormal_regularized_horseshoe$RStan[, "Jul"]), col = "#D55E00", lty = 2, lwd = 2)
    legend("topright", c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
           col = c("grey70", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)
  })

  ### Orthonormal coding / Evidence for each prior relative to ridge ----
  scenario_text("bridge-orthonormal_ridge", fit_orthonormal_ridge$marglik)
  scenario_text("stan-bridge-orthonormal_ridge", ref_orthonormal_ridge$marglik)
  scenario_text("model-comparison-orthonormal", regularized_evidence_table(list(
    "Ridge" = list(BayesTools = fit_orthonormal_ridge$marglik, RStan = ref_orthonormal_ridge$marglik),
    "Lasso" = list(BayesTools = fit_orthonormal_lasso$marglik, RStan = ref_orthonormal_lasso$marglik),
    "Horseshoe" = list(BayesTools = fit_orthonormal_horseshoe$marglik, RStan = ref_orthonormal_horseshoe$marglik),
    "Finite-slab horseshoe" = list(BayesTools = fit_orthonormal_regularized_horseshoe$marglik, RStan = ref_orthonormal_regularized_horseshoe$marglik))))

  ### The same July–May difference under both codings ----
  # Subtract level effects within each draw to retain their correlation.
  july_minus_may <- function(draws) draws[, "Jul"] - draws[, "May"]
  scenario_comparison("july-may-by-coding", c("Mean difference", "Orthonormal"),
    main = "July minus May under the horseshoe prior", show_table = TRUE,
    posteriors = list(
      BayesTools = cbind("Mean difference" = july_minus_may(draws_meandif_horseshoe$BayesTools),
                         Orthonormal = july_minus_may(draws_orthonormal_horseshoe$BayesTools)),
      RStan = cbind("Mean difference" = july_minus_may(draws_meandif_horseshoe$RStan),
                    Orthonormal = july_minus_may(draws_orthonormal_horseshoe$RStan))))

  ### Prior and posterior distributions for all monthly deviations and the wind slope ----
  # Reuse the fitted objects above. Each figure uses the same coefficient axes.
  fits <- list(meandif_ridge = jags_meandif_ridge,
               meandif_lasso = jags_meandif_lasso,
               meandif_horseshoe = jags_meandif_horseshoe,
               meandif_regularized_horseshoe = jags_meandif_regularized_horseshoe,
               orthonormal_ridge = jags_orthonormal_ridge,
               orthonormal_lasso = jags_orthonormal_lasso,
               orthonormal_horseshoe = jags_orthonormal_horseshoe,
               orthonormal_regularized_horseshoe = jags_orthonormal_regularized_horseshoe)
  draws <- list(meandif_ridge = draws_meandif_ridge$RStan,
                meandif_lasso = draws_meandif_lasso$RStan,
                meandif_horseshoe = draws_meandif_horseshoe$RStan,
                meandif_regularized_horseshoe = draws_meandif_regularized_horseshoe$RStan,
                orthonormal_ridge = draws_orthonormal_ridge$RStan,
                orthonormal_lasso = draws_orthonormal_lasso$RStan,
                orthonormal_horseshoe = draws_orthonormal_horseshoe$RStan,
                orthonormal_regularized_horseshoe = draws_orthonormal_regularized_horseshoe$RStan)
  titles <- c(meandif_ridge = "Ridge with mean-difference month contrasts",
              meandif_lasso = "Group lasso with mean-difference month contrasts",
              meandif_horseshoe = "Horseshoe with mean-difference month contrasts",
              meandif_regularized_horseshoe = "Finite-slab horseshoe with mean-difference month contrasts",
              orthonormal_ridge = "Ridge with orthonormal month contrasts",
              orthonormal_lasso = "Group lasso with orthonormal month contrasts",
              orthonormal_horseshoe = "Horseshoe with orthonormal month contrasts",
              orthonormal_regularized_horseshoe = "Finite-slab horseshoe with orthonormal month contrasts")
  selectors <- setNames(c(paste0("mu_Month[", months, "]"), "mu_Wind"), effects)
  xlim <- c(-1.5, 1.5)
  for (model in names(fits)) {
    scenario_plot(paste0("distributions-", model), {
      catalog <- parameter_catalog(fits[[model]])
      posteriors <- lapply(selectors, function(term) {
        parameter_mixed_posterior(fits[[model]], parameter_catalog_resolve(catalog, term))
      })
      priors <- lapply(posteriors, posterior_metadata, field = "prior_density")
      stan_curves <- lapply(names(selectors), function(effect) density(draws[[model]][, effect], n = 301))
      # Mark an unbounded peak with an arrow from the nearest finite grid ordinate.
      unbounded <- vapply(priors, function(p) prior_density_ordinate(p, 0)$behavior == "infinite", logical(1))
      prior_height <- vapply(seq_along(priors), function(i) {
        x <- if (unbounded[i]) 0.01 else 0
        exp(prior_density_ordinate(priors[[i]], x)$log_density)
      }, numeric(1))
      ymax <- ceiling(1.12 * max(prior_height, vapply(stan_curves, function(d) max(d$y), numeric(1)),
                                 vapply(posteriors, function(q) max(density(q)$y), numeric(1))))
      old_par <- par(no.readonly = TRUE)
      regularized_plot_layout(titles[[model]])
      for (i in seq_along(selectors)) {
        effect <- names(selectors)[i]
        plot(priors[[i]], main = effect, xlab = "", ylab = "", xlim = xlim, ylim = c(0, ymax),
             n_points = 301, col = "grey70", lwd = 2)
        if (unbounded[i]) arrows(0, prior_height[i], 0, ymax, col = "grey70", lwd = 1.5, length = 0.065)
        plot_posterior(posteriors[i], effect, add = TRUE, xlim = xlim, n_points = 301, col = "#0072B2", lwd = 2)
        lines(stan_curves[[i]], col = "#D55E00", lty = 2, lwd = 2)
      }
      mtext("Density", side = 2, outer = TRUE, line = 0.3)
      mtext("Standardized effect", side = 1, outer = TRUE, line = 0)
      if (any(unbounded)) mtext("Arrows mark unbounded prior density at zero.", side = 1, outer = TRUE, line = 1.4, cex = 0.8)
      par(old_par)
    })
  }
})
