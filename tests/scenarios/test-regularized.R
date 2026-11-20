if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("regularized")
source(file.path(.scenario_config()$root, "helpers", "regularized-display.R"), local = TRUE)
# testthat::test_file("tests/scenarios/test-regularized.R")

### Fertility in 47 Swiss provinces
# How do the five slopes change under four shrinkage priors? Fit the same
# standardized Gaussian regression in both engines, with residual SD fixed at 1.
testthat::test_that("Shrinkage of Swiss fertility predictors", {

  dat <- as.data.frame(scale(datasets::swiss))
  slopes <- ~ Agriculture + Examination + Education + Catholic + Infant.Mortality
  X <- model.matrix(slopes, dat)[, -1]
  predictors <- colnames(X)
  data_list <- list(N = nrow(dat), y = dat$Fertility)
  stan_data <- c(data_list, list(K = ncol(X), X = X))
  # Each Swiss slope has its own local horseshoe scale.
  horseshoe_data <- c(stan_data, list(J = ncol(X), unit = seq_len(ncol(X))))

  likelihood <- "model {
    for (i in 1:N) { y[i] ~ dnorm(mu[i], 1) }
  }"
  log_likelihood <- function(parameters, data) {
    sum(dnorm(data$y, parameters$mu, sd = 1, log = TRUE))
  }
  intercept <- prior("normal", list(0, 1))
  global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))
  fit_version <- regularized_fit_version(dat, data_list, deparse(slopes), likelihood, deparse(body(log_likelihood)), intercept, global)

  ### Ridge ----
  # The intercept has its own Normal prior. All slopes share the global scale.
  fit_ridge <- scenario_fit("fit_ridge", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = intercept,
        shrinkage = prior_regularized("ridge", scale = global))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 238354864)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_ridge <- fit_ridge$fit
  scenario_text("estimates-ridge", JAGS_estimates_table(jags_ridge))

  ref_ridge <- scenario_fit("ref_ridge", {
    fit <- rstan::sampling(scenario_stan_model("ridge.stan"), data = stan_data,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(238354865)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(stan_data, scenario_stan_version("ridge.stan")))
  stan_ridge <- ref_ridge$fit
  scenario_text("stan-estimates-ridge", round(rstan::summary(stan_ridge, pars = c("alpha", "beta", "tau"))$summary, 3))

  # Compare like with like. Stan numbers the slopes in the column order of X.
  parameters <- c("Intercept", predictors, "Global scale")
  jags_names <- setNames(c("mu_intercept", paste0("mu_", predictors), "mu_shrinkage_global"), parameters)
  stan_names <- setNames(c("alpha", sprintf("beta[%d]", seq_along(predictors)), "tau"), parameters)

  draws_ridge <- list(
    BayesTools = scenario_draws(jags_ridge$mcmc, jags_names),
    RStan = scenario_draws(rstan::As.mcmc.list(stan_ridge), stan_names))
  scenario_comparison("comparison-ridge", parameters, posteriors = draws_ridge,
                      main = "Swiss fertility: ridge", show_table = TRUE)

  # Inspect one coefficient in the native BayesTools plot and add the Stan density.
  scenario_plot("ridge-education-posterior", {
    selection <- parameter_catalog_resolve(parameter_catalog(jags_ridge), "(mu) Education")
    posterior <- parameter_mixed_posterior(jags_ridge, selection)
    plot_posterior(list(Education = posterior), "Education", prior = TRUE,
                   main = "Swiss fertility and education", xlim = c(-1.5, 1.5),
                   xlab = "Change in standardized fertility per SD of education",
                   col = "#0072B2", lwd = 2, dots_prior = list(col = "grey70", lwd = 2))
    lines(density(draws_ridge$RStan[, "Education"]), col = "#D55E00", lty = 2, lwd = 2)
    legend("topright", c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
           col = c("grey70", "#0072B2", "#D55E00"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)
  })

  ### Lasso ----
  fit_lasso <- scenario_fit("fit_lasso", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = intercept,
        shrinkage = prior_regularized("lasso", scale = global))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 256602512)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_lasso <- fit_lasso$fit
  scenario_text("estimates-lasso", JAGS_estimates_table(jags_lasso))

  ref_lasso <- scenario_fit("ref_lasso", {
    fit <- rstan::sampling(scenario_stan_model("lasso.stan"), data = stan_data,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(256602513)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(stan_data, scenario_stan_version("lasso.stan")))
  stan_lasso <- ref_lasso$fit
  scenario_text("stan-estimates-lasso", round(rstan::summary(stan_lasso, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_lasso <- list(
    BayesTools = scenario_draws(jags_lasso$mcmc, jags_names),
    RStan = scenario_draws(rstan::As.mcmc.list(stan_lasso), stan_names))
  scenario_comparison("comparison-lasso", parameters, posteriors = draws_lasso,
                      main = "Swiss fertility: lasso", show_table = TRUE)

  ### Fertility at three levels of education ----
  # Hold the other four predictors at their means. Compare the expected
  # standardized fertility at education one SD below, at, and one SD above its mean.
  education_grid <- data.frame(Agriculture = 0, Examination = 0, Education = c(-1, 0, 1),
                               Catholic = 0, Infant.Mortality = 0)
  fertility <- JAGS_predict_formula(jags_lasso, "mu", data = education_grid, formula_target = "fixed")
  stan_coef <- rstan::extract(stan_lasso, pars = c("alpha", "beta"))
  stan_fertility <- rbind(stan_coef$alpha - stan_coef$beta[, 3], stan_coef$alpha,
                          stan_coef$alpha + stan_coef$beta[, 3])
  scenario_text("fertility-by-education", scenario_prediction_table(education_grid["Education"],
    list(BayesTools = fertility$value, RStan = stan_fertility)))

  # Do provinces with higher education have lower expected fertility at the
  # same values of the other predictors? The symmetric slope prior gives
  # equal prior odds for the two directions.
  education_direction <- hypothesis_BF(jags_lasso, "Education < 0 vs Education > 0", seed = 1)
  scenario_text("education-direction", education_direction)
  education_lower <- draws_lasso$RStan[, "Education"] < 0
  # Few draws have a positive slope, so small probability differences are magnified on the odds scale.
  scenario_text("education-direction-comparison", data.frame(
    engine = c("BayesTools", "RStan"),
    Pr_lower = round(c(mean(draws_lasso$BayesTools[, "Education"] < 0), mean(education_lower)), 3),
    BF = round(c(as.data.frame(education_direction)$BF10, sum(education_lower) / sum(!education_lower)), 2)))

  ### Horseshoe ----
  fit_horseshoe <- scenario_fit("fit_horseshoe", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = intercept,
        shrinkage = prior_regularized("horseshoe", scale = global))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 104941567)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_horseshoe <- fit_horseshoe$fit
  scenario_text("estimates-horseshoe", JAGS_estimates_table(jags_horseshoe))

  ref_horseshoe <- scenario_fit("ref_horseshoe", {
    fit <- rstan::sampling(scenario_stan_model("horseshoe.stan"), data = horseshoe_data,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(104941568)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(horseshoe_data, scenario_stan_version("horseshoe.stan")))
  stan_horseshoe <- ref_horseshoe$fit
  scenario_text("stan-estimates-horseshoe", round(rstan::summary(stan_horseshoe, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_horseshoe <- list(
    BayesTools = scenario_draws(jags_horseshoe$mcmc, jags_names),
    RStan = scenario_draws(rstan::As.mcmc.list(stan_horseshoe), stan_names))
  scenario_comparison("comparison-horseshoe", parameters, posteriors = draws_horseshoe,
                      main = "Swiss fertility: horseshoe", show_table = TRUE)

  ### Finite-slab horseshoe ----
  # The Student-t slab has scale 2 and four degrees of freedom.
  fit_regularized_horseshoe <- scenario_fit("fit_regularized_horseshoe", {
    fit <- JAGS_fit(likelihood, data_list, prior_list = list(),
      formula_list = list(mu = slopes), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = intercept,
        shrinkage = prior_regularized("horseshoe", scale = global, slab_scale = 2, slab_df = 4))),
      chains = 3, adapt = 1500, burnin = 3000, sample = 15000, seed = 89211)
    ml <- JAGS_bridgesampling(fit, log_posterior = log_likelihood, data = data_list, seed = 152666980)
    list(fit = fit, marglik = ml)
  }, cache_version = fit_version)
  jags_regularized_horseshoe <- fit_regularized_horseshoe$fit
  scenario_text("estimates-regularized_horseshoe", JAGS_estimates_table(jags_regularized_horseshoe))

  ref_regularized_horseshoe <- scenario_fit("ref_regularized_horseshoe", {
    fit <- rstan::sampling(scenario_stan_model("finite-slab-horseshoe.stan"), data = horseshoe_data,
      chains = 3, cores = 3, warmup = 1500, iter = 6500, seed = 26011,
      control = list(adapt_delta = 0.995, max_treedepth = 14), refresh = 0)
    set.seed(152666981)
    ml <- bridgesampling::bridge_sampler(fit, method = "normal", use_neff = TRUE, maxiter = 10000, silent = TRUE)
    list(fit = fit, marglik = ml)
  }, cache_version = regularized_fit_version(horseshoe_data, scenario_stan_version("finite-slab-horseshoe.stan")))
  stan_regularized_horseshoe <- ref_regularized_horseshoe$fit
  scenario_text("stan-estimates-regularized_horseshoe", round(rstan::summary(stan_regularized_horseshoe, pars = c("alpha", "beta", "tau"))$summary, 3))

  draws_regularized_horseshoe <- list(
    BayesTools = scenario_draws(jags_regularized_horseshoe$mcmc, jags_names),
    RStan = scenario_draws(rstan::As.mcmc.list(stan_regularized_horseshoe), stan_names))
  scenario_comparison("comparison-regularized_horseshoe", parameters, posteriors = draws_regularized_horseshoe,
                      main = "Swiss fertility: horseshoe with a finite slab", show_table = TRUE)

  ### Evidence for each prior relative to ridge ----
  # The Examination local scale can mix slowly in the finite-slab horseshoe.
  scenario_text("bridge-ridge", fit_ridge$marglik)
  scenario_text("stan-bridge-ridge", ref_ridge$marglik)
  scenario_text("model-comparison", regularized_evidence_table(list(
    "Ridge" = list(BayesTools = fit_ridge$marglik, RStan = ref_ridge$marglik),
    "Lasso" = list(BayesTools = fit_lasso$marglik, RStan = ref_lasso$marglik),
    "Horseshoe" = list(BayesTools = fit_horseshoe$marglik, RStan = ref_horseshoe$marglik),
    "Finite-slab horseshoe" = list(BayesTools = fit_regularized_horseshoe$marglik, RStan = ref_regularized_horseshoe$marglik))))

  # If each shrinkage prior is equally plausible before seeing the data,
  # how much posterior model probability does each receive?
  fitted_models <- list(Ridge = fit_ridge, Lasso = fit_lasso, Horseshoe = fit_horseshoe,
                        "Finite-slab horseshoe" = fit_regularized_horseshoe)
  models <- models_inference(lapply(fitted_models, function(x) {
    list(model = x$fit, marglik = x$marglik, prior_weights = 1)
  }))
  stan_logml <- c(ref_ridge$marglik$logml, ref_lasso$marglik$logml,
                  ref_horseshoe$marglik$logml, ref_regularized_horseshoe$marglik$logml)
  stan_weights <- exp(stan_logml - max(stan_logml))
  scenario_text("model-probabilities", data.frame(
    model = names(fitted_models), prior = 1 / length(models),
    BayesTools = round(vapply(models, function(x) x$inference$post_prob, numeric(1)), 3),
    RStan = round(stan_weights / sum(stan_weights), 3), row.names = NULL))

  ### Prior and posterior distributions for all five slopes ----
  # Reuse the fitted objects above. Each figure uses the same coefficient axes.
  fits <- list(ridge = jags_ridge,
               lasso = jags_lasso,
               horseshoe = jags_horseshoe,
               regularized_horseshoe = jags_regularized_horseshoe)
  draws <- list(ridge = draws_ridge$RStan,
                lasso = draws_lasso$RStan,
                horseshoe = draws_horseshoe$RStan,
                regularized_horseshoe = draws_regularized_horseshoe$RStan)
  titles <- c(ridge = "Swiss fertility: ridge",
              lasso = "Swiss fertility: lasso",
              horseshoe = "Swiss fertility: horseshoe",
              regularized_horseshoe = "Swiss fertility: horseshoe with a finite slab")
  selectors <- setNames(paste0("mu_", predictors), predictors)
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
