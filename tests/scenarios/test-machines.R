if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("machines")
# testthat::test_file("tests/scenarios/test-machines.R")

### Description
# factor random slopes and random-effect covariance structures for the nlme
# Machines data (6 workers each scoring 3 times on 3 machines): unstructured
# machine effects per worker, random intercepts with treatment and
# mean-difference machine slopes, diagonal, identity, compound-symmetry,
# heterogeneous compound-symmetry, and nested worker/machine effects, compared
# with lme4 structured covariances and with brms under matched priors where
# brms has the structure; bridge-sampling marginal likelihoods compare the
# structures, and predictions for a new worker compare the new-level policies.
# lme4 confirms two equivalences: the treatment-slope model reparameterizes the
# unstructured one, and compound symmetry with a positive correlation is the
# nested model (equal log-likelihoods)
#
# TODO: with six workers the fixed machine effects trade off with the means of
# the noncentered worker effects, so these multi-column blocks mix slowly
# (an ESS of about 40 from 16000 draws); the mean-centered parameterization is
# available only for a scalar random intercept. The fits use long thinned
# chains for now; consider other parameterizations for multi-column and
# structured blocks.

testthat::test_that("Machines factor random slopes and covariance structures", {

  data("Machines", package = "nlme")
  Machines <- as.data.frame(Machines)
  Machines$Worker  <- factor(as.character(Machines$Worker), levels = as.character(1:6))
  Machines$score_z <- as.numeric(scale(Machines$score))

  likelihood <- "model{
    for(i in 1:N){
      score_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  data_list <- list(score_z = Machines$score_z, N = nrow(Machines))
  log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["score_z"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  ### BayesTools fits ----
  # priors: intercept ~ N(0, 1), one N(0, 1) per treatment-contrast machine
  # coefficient, half-normal(0, 1) for sigma and the SDs, LKJ(1) for
  # unstructured correlations, and the default uniform(-1/2, 1) CS/HCS
  # correlation; long thinned chains (see the TODO above). Each block returns
  # the fit and its bridge-sampling marginal likelihood
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  fixed_priors  <- list(mu = list(
    intercept = prior("normal", list(mean = 0, sd = 1)),
    Machine   = prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
  ))
  marglik <- function(fit) {
    JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, prior_list = list(sigma = unit_sd_prior), seed = 1)
  }
  structures <- list()
  # one correlated coefficient per machine (no random intercept)
  structures[["us"]] <- scenario_fit("fit_us", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + us(0 + Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(
        sd = unit_sd_prior, cor = prior_lkj(eta = 1), contrasts = c(Machine = "independent")))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # a random intercept with correlated treatment-contrast machine slopes
  structures[["slope"]] <- scenario_fit("fit_slope", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + (1 + Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(
        sd = unit_sd_prior, cor = prior_lkj(eta = 1)))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # mean-difference machine slopes: one exchangeable unit with a shared SD
  structures[["meandif"]] <- scenario_fit("fit_meandif", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + (1 + Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(
        sd = unit_sd_prior, cor = prior_lkj(eta = 1), contrasts = c(Machine = "meandif")))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # independent machine coefficients with one SD each
  structures[["diag"]] <- scenario_fit("fit_diag", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + diag(0 + Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(
        sd = unit_sd_prior, contrasts = c(Machine = "independent")))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # independent machine coefficients with one shared SD
  structures[["id"]] <- scenario_fit("fit_id", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + id(0 + Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(
        sd = unit_sd_prior, contrasts = c(Machine = "independent")))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # compound symmetry over the machines: one SD and one correlation
  structures[["cs"]] <- scenario_fit("fit_cs", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + cs(Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(sd = unit_sd_prior))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # heterogeneous compound symmetry: one SD per machine and one correlation
  structures[["hcs"]] <- scenario_fit("fit_hcs", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + hcs(Machine | Worker)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(Worker = random_block(sd = unit_sd_prior))),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })
  # worker intercepts and worker-by-machine intercepts
  structures[["nested"]] <- scenario_fit("fit_nested", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = list(mu = ~ Machine + (1 | Worker/Machine)), formula_data_list = list(mu = Machines),
      formula_prior_list = fixed_priors,
      formula_random_prior_list = list(mu = prior_random(sd = unit_sd_prior)),
      chains = 4, adapt = 1000, burnin = 5000, sample = 10000, thin = 10, seed = 1
    )
    list(fit = fit, marglik = marglik(fit))
  })

  ### lme4 fits ----
  # a single structured term, so lme4 2.0 cannot misassign its structure
  lme4_fits <- list(
    us     = lme4::lmer(score_z ~ Machine + (0 + Machine | Worker),                  data = Machines, REML = FALSE),
    slope  = lme4::lmer(score_z ~ Machine + (1 + Machine | Worker),                  data = Machines, REML = FALSE),
    diag   = lme4::lmer(score_z ~ Machine + diag(0 + Machine | Worker),              data = Machines, REML = FALSE),
    id     = lme4::lmer(score_z ~ Machine + diag(0 + Machine | Worker, hom = TRUE),  data = Machines, REML = FALSE),
    cs     = lme4::lmer(score_z ~ Machine + cs(0 + Machine | Worker, hom = TRUE),    data = Machines, REML = FALSE),
    hcs    = lme4::lmer(score_z ~ Machine + cs(0 + Machine | Worker),                data = Machines, REML = FALSE),
    nested = lme4::lmer(score_z ~ Machine + (1 | Worker/Machine),                    data = Machines, REML = FALSE)
  )

  ### brms fits ----
  # the structures brms has, with the same priors (brms bounds the SD priors
  # at zero); identity machine effects within a worker are independent
  # worker-by-machine intercepts with one SD, (1 | Worker:Machine)
  brms_priors <- c(
    brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 1), class = sd),
    brms::prior(normal(0, 1), class = sigma)
  )
  brms_marglik <- function(fit) {
    set.seed(1)
    brms::bridge_sampler(fit, silent = TRUE)
  }
  brms_fits <- list()
  brms_fits[["us"]] <- scenario_fit("brms_us", {
    fit <- brms::brm(
      score_z ~ 0 + Intercept + Machine + (0 + Machine | Worker), data = Machines,
      prior = c(brms_priors, brms::prior(lkj(1), class = cor)), save_pars = brms::save_pars(all = TRUE),
      backend = "rstan", chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["slope"]] <- scenario_fit("brms_slope", {
    fit <- brms::brm(
      score_z ~ 0 + Intercept + Machine + (1 + Machine | Worker), data = Machines,
      prior = c(brms_priors, brms::prior(lkj(1), class = cor)), save_pars = brms::save_pars(all = TRUE),
      backend = "rstan", chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["diag"]] <- scenario_fit("brms_diag", {
    fit <- brms::brm(
      score_z ~ 0 + Intercept + Machine + (0 + Machine || Worker), data = Machines,
      prior = brms_priors, save_pars = brms::save_pars(all = TRUE),
      backend = "rstan", chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["id"]] <- scenario_fit("brms_id", {
    fit <- brms::brm(
      score_z ~ 0 + Intercept + Machine + (1 | Worker:Machine), data = Machines,
      prior = brms_priors, save_pars = brms::save_pars(all = TRUE),
      backend = "rstan", chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    list(fit = fit, marglik = brms_marglik(fit))
  })
  brms_fits[["nested"]] <- scenario_fit("brms_nested", {
    fit <- brms::brm(
      score_z ~ 0 + Intercept + Machine + (1 | Worker/Machine), data = Machines,
      prior = brms_priors, save_pars = brms::save_pars(all = TRUE),
      backend = "rstan", chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    list(fit = fit, marglik = brms_marglik(fit))
  })

  ### model summaries ----
  scenario_text("estimates-us",      JAGS_estimates_table(structures[["us"]][["fit"]]))
  scenario_text("estimates-slope",   JAGS_estimates_table(structures[["slope"]][["fit"]]))
  scenario_text("estimates-meandif", JAGS_estimates_table(structures[["meandif"]][["fit"]]))
  scenario_text("estimates-diag",    JAGS_estimates_table(structures[["diag"]][["fit"]]))
  scenario_text("estimates-id",      JAGS_estimates_table(structures[["id"]][["fit"]]))
  scenario_text("estimates-cs",      JAGS_estimates_table(structures[["cs"]][["fit"]]))
  scenario_text("estimates-hcs",     JAGS_estimates_table(structures[["hcs"]][["fit"]]))
  scenario_text("estimates-nested",  JAGS_estimates_table(structures[["nested"]][["fit"]]))

  ### exact marginal likelihoods ----
  # for the structures without an LKJ correlation matrix, the marginal
  # likelihood reduces to a low-dimensional integral over sigma, the SDs, and
  # the CS/HCS correlation (uniform(-1/2, 1) prior), which importance sampling
  # evaluates with an error far below that of bridge sampling
  exact_logml <- scenario_fit("exact_marglik", {
    same_worker  <- outer(Machines$Worker, Machines$Worker, "==") * 1
    same_cell    <- same_worker * outer(Machines$Machine, Machines$Machine, "==")
    machine_cell <- lapply(levels(Machines$Machine), function(m) same_cell * outer(Machines$Machine == m, Machines$Machine == m))
    machine_sd   <- function(sds) sds[as.integer(Machines$Machine)]
    residual     <- diag(nrow(Machines))
    # theta: sigma, the SDs, and the correlation, in the order of `columns`
    covariances <- list(
      id     = function(t) t[2]^2 * same_cell + t[1]^2 * residual,
      diag   = function(t) Reduce(`+`, Map(`*`, t[2:4]^2, machine_cell)) + t[1]^2 * residual,
      nested = function(t) t[2]^2 * same_worker + t[3]^2 * same_cell + t[1]^2 * residual,
      cs     = function(t) t[2]^2 * (t[3] * same_worker + (1 - t[3]) * same_cell) + t[1]^2 * residual,
      hcs    = function(t) outer(machine_sd(t[2:4]), machine_sd(t[2:4])) * (t[5] * same_worker + (1 - t[5]) * same_cell) +
        t[1]^2 * residual
    )
    machine_sds <- paste0("(mu) sd(Machine[", levels(Machines$Machine), "])")
    columns <- list(
      id = c("sigma", "(mu) sd"), diag = c("sigma", machine_sds),
      nested = c("sigma", "(mu) Worker: sd(intercept)", "(mu) Machine_Worker: sd(intercept)"),
      cs = c("sigma", "(mu) sd", "(mu) cor"), hcs = c("sigma", machine_sds, "(mu) cor")
    )
    set.seed(1)
    vapply(names(covariances), function(structure) {
      draws <- JAGS_estimates_table(structures[[structure]][["fit"]], return_samples = TRUE)[, columns[[structure]]]
      scenario_gaussian_logml(
        Machines$score_z, stats::model.matrix(~ Machine, Machines), covariances[[structure]], draws,
        correlation_bounds = list("(mu) cor" = c(-1 / 2, 1))
      )[["logml"]]
    }, numeric(1))
  })

  ### structure comparison ----
  # lme4 log-likelihoods (equal for us and slope, and for cs and nested), and
  # log marginal likelihoods (exact where available, bridge sampling in
  # BayesTools and brms) with log BF against the unstructured model. Bridge
  # sampling errors exceed their reported errors: the BayesTools HCS estimate
  # (seed 1) lies 0.26 below the exact value (other seeds give -11.84 to
  # -11.98) and the brms nested estimate 0.59 above it
  scenario_text("structure-comparison", {
    log_ml <- function(fits) vapply(fits, function(x) x[["marglik"]][["logml"]], numeric(1))
    bayes  <- log_ml(structures)
    brms   <- log_ml(brms_fits)
    out <- data.frame(
      structure      = names(structures),
      lme4_logLik    = vapply(names(structures), function(s) {
        if (s %in% names(lme4_fits)) as.numeric(stats::logLik(lme4_fits[[s]])) else NA_real_
      }, numeric(1)),
      exact_logml      = unname(exact_logml[names(structures)]),
      BayesTools_logml = bayes,
      brms_logml       = unname(brms[names(structures)]),
      BayesTools_log_BF_vs_us = bayes - bayes[["us"]],
      brms_log_BF_vs_us       = unname(brms[names(structures)]) - brms[["us"]],
      row.names = NULL
    )
    out[, -1] <- round(out[, -1], 3)
    out
  })

  ### lme4 and brms comparisons ----
  # each structure's parameters with BayesTools labels; lme4 values come from
  # the worker covariance matrix, brms columns are named per structure
  coefficient_labels <- c("(Intercept)" = "intercept", MachineA = "Machine[A]", MachineB = "Machine[B]", MachineC = "Machine[C]")
  lme4_values <- function(model) {
    fixed <- lme4::fixef(model)
    c(intercept = fixed[["(Intercept)"]], "Machine[B]" = fixed[["MachineB"]], "Machine[C]" = fixed[["MachineC"]],
      sigma = stats::sigma(model))
  }
  lme4_worker <- function(model) {
    covariance <- as.matrix(lme4::VarCorr(model)[["Worker"]])
    labels     <- coefficient_labels[rownames(covariance)]
    correlation <- stats::cov2cor(covariance)
    pairs <- which(lower.tri(correlation), arr.ind = TRUE)
    c(stats::setNames(sqrt(diag(covariance)), paste0("sd(", labels, ")")),
      stats::setNames(correlation[pairs], paste0("cor(", labels[pairs[, "col"]], ",", labels[pairs[, "row"]], ")")))
  }
  structure_references <- function(structure, lme4_random) {
    list(lme4 = c(lme4_values(lme4_fits[[structure]]), lme4_random))
  }
  structure_posteriors <- function(structure, parameters, brms_columns = NULL) {
    samples       <- JAGS_estimates_table(structures[[structure]][["fit"]], return_samples = TRUE)
    bayes_columns <- stats::setNames(ifelse(parameters == "sigma", "sigma", paste0("(mu) ", parameters)), parameters)
    fixed_brms    <- c(intercept = "b_Intercept", "Machine[B]" = "b_MachineB", "Machine[C]" = "b_MachineC", sigma = "sigma")
    c(
      list(BayesTools = scenario_draws(samples, bayes_columns)),
      if (!is.null(brms_columns)) list(brms = scenario_draws(as.matrix(brms_fits[[structure]][["fit"]]), c(fixed_brms, brms_columns)))
    )
  }
  fixed_parameters <- c("intercept", "Machine[B]", "Machine[C]", "sigma")
  machine_sd_brms  <- c("sd(Machine[A])" = "sd_Worker__MachineA", "sd(Machine[B])" = "sd_Worker__MachineB",
                        "sd(Machine[C])" = "sd_Worker__MachineC")

  us_random <- lme4_worker(lme4_fits[["us"]])
  scenario_comparison(
    "comparison-us", c(fixed_parameters, names(us_random)), main = "Unstructured machine effects",
    references = structure_references("us", us_random),
    posteriors = structure_posteriors("us", c(fixed_parameters, names(us_random)), c(
      machine_sd_brms,
      "cor(Machine[A],Machine[B])" = "cor_Worker__MachineA__MachineB", "cor(Machine[A],Machine[C])" = "cor_Worker__MachineA__MachineC",
      "cor(Machine[B],Machine[C])" = "cor_Worker__MachineB__MachineC"
    ))
  )
  slope_random <- lme4_worker(lme4_fits[["slope"]])
  scenario_comparison(
    "comparison-slope", c(fixed_parameters, names(slope_random)), main = "Random intercepts and treatment machine slopes",
    references = structure_references("slope", slope_random),
    posteriors = structure_posteriors("slope", c(fixed_parameters, names(slope_random)), c(
      "sd(intercept)" = "sd_Worker__Intercept", "sd(Machine[B])" = "sd_Worker__MachineB", "sd(Machine[C])" = "sd_Worker__MachineC",
      "cor(intercept,Machine[B])" = "cor_Worker__Intercept__MachineB", "cor(intercept,Machine[C])" = "cor_Worker__Intercept__MachineC",
      "cor(Machine[B],Machine[C])" = "cor_Worker__MachineB__MachineC"
    ))
  )
  diag_random <- lme4_worker(lme4_fits[["diag"]])[1:3]
  scenario_comparison(
    "comparison-diag", c(fixed_parameters, names(diag_random)), main = "Independent machine effects",
    references = structure_references("diag", diag_random),
    posteriors = structure_posteriors("diag", c(fixed_parameters, names(diag_random)), machine_sd_brms)
  )
  scenario_comparison(
    "comparison-id", c(fixed_parameters, "sd"), main = "Identity machine effects (one SD)",
    references = structure_references("id", c(sd = lme4_worker(lme4_fits[["id"]])[[1L]])),
    posteriors = structure_posteriors("id", c(fixed_parameters, "sd"), c(sd = "sd_Worker:Machine__Intercept"))
  )
  cs_worker <- lme4_worker(lme4_fits[["cs"]])
  scenario_comparison(
    "comparison-cs", c(fixed_parameters, "sd", "cor"), main = "Compound symmetry",
    references = structure_references("cs", c(sd = cs_worker[[1L]], cor = cs_worker[[4L]])),
    posteriors = structure_posteriors("cs", c(fixed_parameters, "sd", "cor"))
  )
  hcs_worker <- lme4_worker(lme4_fits[["hcs"]])
  scenario_comparison(
    "comparison-hcs", c(fixed_parameters, names(hcs_worker)[1:3], "cor"), main = "Heterogeneous compound symmetry",
    references = structure_references("hcs", c(hcs_worker[1:3], cor = hcs_worker[[4L]])),
    posteriors = structure_posteriors("hcs", c(fixed_parameters, names(hcs_worker)[1:3], "cor"))
  )
  nested_sds    <- vapply(lme4::VarCorr(lme4_fits[["nested"]]), attr, numeric(1), "stddev")
  nested_random <- c("Worker: sd(intercept)" = nested_sds[["Worker"]], "Machine_Worker: sd(intercept)" = nested_sds[["Machine:Worker"]])
  scenario_comparison(
    "comparison-nested", c(fixed_parameters, names(nested_random)), main = "Nested worker and worker-by-machine effects",
    references = structure_references("nested", nested_random),
    posteriors = structure_posteriors("nested", c(fixed_parameters, names(nested_random)), c(
      "Worker: sd(intercept)" = "sd_Worker__Intercept", "Machine_Worker: sd(intercept)" = "sd_Worker:Machine__Intercept"
    ))
  )

  ### predictions for a new worker ----
  # lme4 predicts a new worker at the population level, as the "zero" policy;
  # the "sample" policy draws the new worker's machine effects from each
  # structure's fitted covariance, as brms sample_new_levels = "gaussian"
  scenario_text("predictions-new-worker", {
    new_worker <- data.frame(Worker = factor("new"), Machine = factor(c("A", "B", "C")))
    predict_new <- function(structure, new_levels, seed = NULL) {
      JAGS_predict_formula(structures[[structure]][["fit"]], "mu", data = new_worker, formula_target = "conditional",
                           new_levels = new_levels, seed = seed)$value
    }
    brms_new <- function(structure) {
      t(brms::posterior_epred(brms_fits[[structure]][["fit"]], newdata = new_worker,
                              allow_new_levels = TRUE, sample_new_levels = "gaussian"))
    }
    scenario_prediction_table(new_worker[, "Machine", drop = FALSE], list(
      "lme4"          = stats::predict(lme4_fits[["us"]], newdata = new_worker, allow.new.levels = TRUE),
      "us zero"       = predict_new("us", "zero"),
      "us sample"     = predict_new("us", "sample", seed = 1),
      "us brms"       = brms_new("us"),
      "cs sample"     = predict_new("cs", "sample", seed = 1),
      "nested sample" = predict_new("nested", "sample", seed = 1),
      "nested brms"   = brms_new("nested")
    ))
  })

  ### fitted values ----
  fitted_mean <- function(structure) {
    rowMeans(JAGS_predict_formula(structures[[structure]][["fit"]], "mu", formula_target = "conditional")$value)
  }
  scenario_plot("fitted-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (structure in c("us", "cs", "hcs", "nested")) {
      scenario_agreement_plot(stats::fitted(lme4_fits[[structure]]), fitted_mean(structure), structure)
    }
  })
  scenario_plot("brms-fitted-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (structure in c("us", "slope", "diag", "nested")) {
      scenario_agreement_plot(colMeans(brms::posterior_epred(brms_fits[[structure]][["fit"]])), fitted_mean(structure),
                              structure, reference_label = "brms")
    }
  })
})
