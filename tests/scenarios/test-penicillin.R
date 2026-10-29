if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("penicillin")
# testthat::test_file("tests/scenarios/test-penicillin.R")

### Description
# crossed random intercepts for plates and samples of the lme4 Penicillin data
# (the crossed example of the RandomEffects vignette): independent half-normal
# SD priors and a variance allocation prior, compared with the lme4
# maximum-likelihood fit, including the plate and sample effects and the
# fitted values; the sample block is mean-centered

testthat::test_that("Penicillin crossed random intercepts", {

  data("Penicillin", package = "lme4")
  Penicillin$diameter_z <- as.numeric(scale(Penicillin$diameter))

  likelihood <- "model{
    for(i in 1:N){
      diameter_z[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  data_list     <- list(diameter_z = Penicillin$diameter_z, N = nrow(Penicillin))
  unit_prior    <- prior("normal", list(mean = 0, sd = 1))
  unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
  formula_list  <- list(mu = ~ 1 + (1 | plate) + (1 | sample))

  ### model fits ----
  fit_lme4 <- lme4::lmer(diameter_z ~ 1 + (1 | plate) + (1 | sample), data = Penicillin, REML = FALSE)

  fit_independent <- scenario_fit("fit_independent", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = formula_list, formula_data_list = list(mu = Penicillin),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      # with six samples, the intercept and the sample effects trade off; the
      # mean-centered sample block samples them around the intercept
      formula_random_prior_list = list(mu = prior_random(
        plate  = random_block(sd = unit_sd_prior),
        sample = random_block(sd = unit_sd_prior, parameterization = "mean_centered")
      )),
      chains = 4, adapt = 2000, burnin = 2000, sample = 4000, seed = 2
    )
  })
  fit_allocation <- scenario_fit("fit_allocation", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list, prior_list = list(sigma = unit_sd_prior),
      formula_list = formula_list, formula_data_list = list(mu = Penicillin),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(
        sample     = random_block(parameterization = "mean_centered"),
        allocation = random_variance_allocation(
          name    = "total_re",
          terms   = c(plate = "plate", sample = "sample"),
          sd      = unit_sd_prior,
          weights = prior("dirichlet", parameters = list(alpha = c(1, 1)))
        )
      )),
      chains = 4, adapt = 2000, burnin = 2000, sample = 4000, seed = 3
    )
  })

  ### model summaries ----
  fit_lme4
  scenario_text("estimates-independent",      JAGS_estimates_table(fit_independent))
  scenario_text("estimates-independent-full", JAGS_estimates_table(fit_independent, random_effects_summary = "full"))
  scenario_text("estimates-allocation",       JAGS_estimates_table(fit_allocation))
  scenario_text("estimates-allocation-full",  JAGS_estimates_table(fit_allocation, random_effects_summary = "full"))

  scenario_text("convergence", {
    convergence <- lapply(list(independent = fit_independent, allocation = fit_allocation), function(fit) {
      diagnostics <- attr(JAGS_check_convergence(fit), "diagnostics")
      data.frame(
        max_Rhat       = max(diagnostics[["Rhat"]], na.rm = TRUE),
        min_ESS        = min(diagnostics[["ESS"]], na.rm = TRUE),
        worst_Rhat     = diagnostics[["parameter"]][which.max(diagnostics[["Rhat"]])],
        worst_ESS      = diagnostics[["parameter"]][which.min(diagnostics[["ESS"]])]
      )
    })
    data.frame(fit = names(convergence), do.call(rbind, convergence), row.names = NULL)
  })

  ### lme4 comparison ----
  # the allocation quantities of lme4 and of the independent-prior fit are
  # derived from their two SDs
  samples_independent <- JAGS_estimates_table(fit_independent, random_effects_summary = "full", return_samples = TRUE)
  samples_allocation  <- JAGS_estimates_table(fit_allocation,  random_effects_summary = "full", return_samples = TRUE)
  allocation_summaries <- function(sd_plate, sd_sample) {
    variance_total <- sd_plate^2 + sd_sample^2
    cbind(sd_total = sqrt(variance_total), var_prop_plate = sd_plate^2 / variance_total)
  }

  scenario_text("lme4-comparison", {
    lme4_sd <- vapply(lme4::VarCorr(fit_lme4), attr, numeric(1), "stddev")
    lme4_values <- c(
      intercept    = lme4::fixef(fit_lme4)[["(Intercept)"]],
      sigma        = stats::sigma(fit_lme4),
      sd_plate     = lme4_sd[["plate"]],
      sd_sample    = lme4_sd[["sample"]],
      allocation_summaries(lme4_sd[["plate"]], lme4_sd[["sample"]])[1L, ]
    )
    posterior_summary <- function(draws) {
      c(estimate = mean(draws), median = stats::median(draws),
        lCI = unname(stats::quantile(draws, 0.025)), uCI = unname(stats::quantile(draws, 0.975)))
    }
    independent_draws <- cbind(
      samples_independent[, c("(mu) intercept", "sigma", "(mu) plate: sd(intercept)", "(mu) sample: sd(intercept)")],
      allocation_summaries(samples_independent[, "(mu) plate: sd(intercept)"], samples_independent[, "(mu) sample: sd(intercept)"])
    )
    allocation_draws <- samples_allocation[, c(
      "(mu) intercept", "sigma", "(mu) plate: sd(intercept)", "(mu) sample: sd(intercept)",
      "(mu) total_re: sd_total", "(mu) total_re: var_prop(plate)"
    )]
    comparison <- rbind(
      data.frame(fit = "lme4", parameter = names(lme4_values), estimate = unname(lme4_values),
                 median = NA_real_, lCI = NA_real_, uCI = NA_real_),
      data.frame(fit = "independent", parameter = names(lme4_values),
                 t(apply(independent_draws, 2, posterior_summary)), row.names = NULL),
      data.frame(fit = "allocation", parameter = names(lme4_values),
                 t(apply(allocation_draws, 2, posterior_summary)), row.names = NULL)
    )
    comparison[, 3:6] <- round(comparison[, 3:6], 3)
    comparison <- comparison[order(match(comparison$parameter, names(lme4_values))), ]
    rownames(comparison) <- NULL
    comparison
  })

  ### allocation identities ----
  # the allocation reports one total SD and its variance shares; they must
  # reproduce the block SDs draw by draw
  scenario_text("allocation-identities", {
    sd_plate  <- samples_allocation[, "(mu) plate: sd(intercept)"]
    sd_sample <- samples_allocation[, "(mu) sample: sd(intercept)"]
    sd_total  <- samples_allocation[, "(mu) total_re: sd_total"]
    prop      <- samples_allocation[, c("(mu) total_re: var_prop(plate)", "(mu) total_re: var_prop(sample)")]
    data.frame(
      identity      = c(
        "sd_total^2 = sd(plate)^2 + sd(sample)^2",
        "var_prop(plate) + var_prop(sample) = 1",
        "var_prop(plate) = sd(plate)^2 / sd_total^2"
      ),
      max_abs_error = signif(c(
        max(abs(sd_total^2 - sd_plate^2 - sd_sample^2)),
        max(abs(rowSums(prop) - 1)),
        max(abs(prop[, 1] - sd_plate^2 / sd_total^2))
      ), 3)
    )
  })

  ### random effects and fitted values ----
  # the effect of a level is its conditional prediction with only its block
  # minus the fixed part
  level_effects <- function(fit, block) {
    prediction <- JAGS_predict_formula(fit, "mu", formula_target = "conditional", blocks = block, components = TRUE)
    effects    <- rowMeans(prediction$value - prediction$components$fixed)
    levels     <- as.character(Penicillin[[block]])
    effects[match(levels(Penicillin[[block]]), levels)]
  }
  ranef_lme4 <- lme4::ranef(fit_lme4)

  scenario_text("ranef-comparison", {
    fits <- list(independent = fit_independent, allocation = fit_allocation)
    rows <- list()
    for (block in c("plate", "sample")) {
      reference <- ranef_lme4[[block]][levels(Penicillin[[block]]), "(Intercept)"]
      for (fit_name in names(fits)) {
        estimate <- level_effects(fits[[fit_name]], block)
        rows[[length(rows) + 1L]] <- data.frame(
          fit                = fit_name,
          block              = block,
          levels             = length(reference),
          correlation        = round(stats::cor(reference, estimate), 4),
          max_abs_difference = round(max(abs(estimate - reference)), 3),
          mean_difference    = round(mean(estimate - reference), 3)
        )
      }
    }
    do.call(rbind, rows)
  })

  scenario_plot("ranef-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    for (block in c("plate", "sample")) {
      reference <- ranef_lme4[[block]][levels(Penicillin[[block]]), "(Intercept)"]
      scenario_agreement_plot(reference, level_effects(fit_independent, block), paste0(block, " effects (independent)"))
      scenario_agreement_plot(reference, level_effects(fit_allocation,  block), paste0(block, " effects (allocation)"))
    }
  })

  scenario_plot("fitted-comparison", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    reference <- stats::fitted(fit_lme4)
    for (fit_name in c("independent", "allocation")) {
      fit <- if (fit_name == "independent") fit_independent else fit_allocation
      prediction <- JAGS_predict_formula(fit, "mu", formula_target = "conditional")
      scenario_agreement_plot(reference, rowMeans(prediction$value), paste0("Fitted values (", fit_name, ")"))
    }
  })

  ### prior and posterior distributions ----
  posterior_of <- function(fit, name) {
    selection <- parameter_catalog_resolve(parameter_catalog(fit), name)
    stats::setNames(list(parameter_mixed_posterior(fit, selection)), name)
  }
  plot_quantity <- function(fit, name) {
    plot_posterior(posterior_of(fit, name), name, prior = TRUE, main = name)
  }

  scenario_plot("posterior-independent", {
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    plot_quantity(fit_independent, "(mu) plate: sd(intercept)")
    plot_quantity(fit_independent, "(mu) sample: sd(intercept)")
    plot_quantity(fit_independent, "sigma")
  })

  scenario_plot("posterior-allocation", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    plot_quantity(fit_allocation, "(mu) total_re: sd_total")
    plot_quantity(fit_allocation, "(mu) total_re: var_prop(plate)")
    plot_quantity(fit_allocation, "(mu) plate: sd(intercept)")
    plot_quantity(fit_allocation, "(mu) sample: sd(intercept)")
  })
})
