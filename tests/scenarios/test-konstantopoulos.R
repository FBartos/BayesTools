if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("konstantopoulos")
# testthat::test_file("tests/scenarios/test-konstantopoulos.R")

### Description
# a three-level meta-analysis of the metadat dat.konstantopoulos2011 data (56
# standardized mean differences from schools within 11 districts) with known
# sampling variances: nested district and school random effects with
# independent half-normal SD priors and with a variance allocation prior,
# compared with metafor::rma.mv (ML) and with brms under matched priors; Bayes
# factors for the district level (bridge sampling in BayesTools and brms, and
# spike and slab), district and school effects against the metafor BLUPs and
# brms, and the pooled effect and a new school in a new district

testthat::test_that("Konstantopoulos three-level meta-analysis", {

  dat <- as.data.frame(metadat::dat.konstantopoulos2011)
  dat$district <- factor(dat$district)
  dat$school   <- factor(dat$school)
  dat$sei      <- sqrt(dat$vi)

  likelihood <- "model{
    for(i in 1:K){
      yi[i] ~ dnorm(mu[i], pow(sei[i], -2))
    }
  }
  "
  data_list <- list(yi = dat$yi, sei = dat$sei, K = nrow(dat))
  log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["yi"]], parameters[["mu"]], data[["sei"]], log = TRUE))
  }
  # priors on the SMD scale: pooled effect ~ N(0, 1), SDs ~ half-normal(0, 0.5)
  unit_prior <- prior("normal", list(mean = 0, sd = 1))
  sd_prior   <- prior("normal", list(mean = 0, sd = 0.5), truncation = list(0, Inf))

  ### metafor fit ----
  metafor_three <- metafor::rma.mv(yi, vi, random = ~ 1 | district/school, data = dat, method = "ML")

  ### BayesTools fits ----
  # with 11 districts the pooled effect and the district effects trade off;
  # the mean-centered district block samples them around the pooled effect
  three_level <- scenario_fit("fit_three_level", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ 1 + (1 | district/school)), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(
        district        = random_block(sd = sd_prior, parameterization = "mean_centered"),
        school_district = random_block(sd = sd_prior)
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  two_level <- scenario_fit("fit_two_level", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ 1 + (1 | district:school)), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(sd = sd_prior)),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  # the total heterogeneity SD with a Dirichlet(1, 1) split of its variance
  # between districts and schools
  fit_allocation <- scenario_fit("fit_allocation", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ 1 + (1 | district/school)), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(
        district   = random_block(parameterization = "mean_centered"),
        allocation = random_variance_allocation(
          name    = "heterogeneity",
          terms   = c(district = "district", school = "school_district"),
          sd      = sd_prior,
          weights = prior("dirichlet", parameters = list(alpha = c(1, 1)))
        )
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
    )
  })
  # a spike-and-slab prior on the district SD
  fit_district_spike <- scenario_fit("fit_district_spike", {
    JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ 1 + (1 | district/school)), formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(
        district        = random_block(sd = prior_spike_and_slab(prior_parameter = sd_prior)),
        school_district = random_block(sd = sd_prior)
      )),
      chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
    )
  })
  fit_three_level <- three_level[["fit"]]

  ### brms fits ----
  # `se(sei)` fixes the sampling SDs; the priors match, and brms bounds the SD
  # priors at zero
  brms_three_level <- scenario_fit("brms_three_level", {
    fit <- brms::brm(
      yi | se(sei) ~ 0 + Intercept + (1 | district/school), data = dat,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 0.5), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99)
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_two_level <- scenario_fit("brms_two_level", {
    fit <- brms::brm(
      yi | se(sei) ~ 0 + Intercept + (1 | district:school), data = dat,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 0.5), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99)
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })

  ### model summaries ----
  scenario_text("estimates-three-level",     JAGS_estimates_table(fit_three_level))
  scenario_text("estimates-allocation-full", JAGS_estimates_table(fit_allocation, random_effects_summary = "full"))

  ### metafor and brms comparison ----
  # the district share of the heterogeneity variance is derived from the two
  # SDs, except in the allocation fit, which samples it
  district_share <- function(sd_district, sd_school) {
    cbind(var_prop_district = sd_district^2 / (sd_district^2 + sd_school^2))
  }
  metafor_sd <- sqrt(metafor_three$sigma2)
  metafor_values <- c(
    intercept = as.numeric(stats::coef(metafor_three)), sd_district = metafor_sd[[1L]], sd_school = metafor_sd[[2L]],
    district_share(metafor_sd[[1L]], metafor_sd[[2L]])[1L, ]
  )
  with_share <- function(draws) cbind(draws, district_share(draws[, "sd_district"], draws[, "sd_school"]))
  bayes_columns <- c(intercept = "(mu) intercept", sd_district = "(mu) district: sd(intercept)",
                     sd_school = "(mu) school_district: sd(intercept)")
  scenario_comparison(
    "comparison-three-level", names(metafor_values), main = "Three-level meta-analysis",
    references = list(metafor = metafor_values),
    posteriors = list(
      BayesTools  = with_share(scenario_draws(JAGS_estimates_table(fit_three_level, return_samples = TRUE), bayes_columns)),
      brms        = with_share(scenario_draws(as.matrix(brms_three_level[["fit"]]), c(
        intercept = "b_Intercept", sd_district = "sd_district__Intercept", sd_school = "sd_district:school__Intercept"))),
      allocation  = scenario_draws(JAGS_estimates_table(fit_allocation, random_effects_summary = "full", return_samples = TRUE),
                                   c(bayes_columns, var_prop_district = "(mu) heterogeneity: var_prop(district)"))
    )
  )

  ### exact marginal likelihoods ----
  # with the pooled effect (N(0, 1)) and the random effects integrated out,
  # the marginal likelihoods are integrals over the half-normal(0, 0.5) SDs,
  # which importance sampling evaluates
  exact_logml <- scenario_fit("exact_marglik", {
    same_district <- outer(dat$district, dat$district, "==") * 1
    same_school   <- same_district * outer(dat$school, dat$school, "==")
    sampling      <- diag(dat$vi)
    pooled        <- matrix(1, nrow(dat))
    set.seed(1)
    c(
      three_level = scenario_gaussian_logml(
        dat$yi, pooled, function(t) t[1]^2 * same_district + t[2]^2 * same_school + sampling,
        JAGS_estimates_table(fit_three_level, return_samples = TRUE)[, c(
          "(mu) district: sd(intercept)", "(mu) school_district: sd(intercept)")],
        sd_prior_sd = 0.5
      )[["logml"]],
      two_level = scenario_gaussian_logml(
        dat$yi, pooled, function(t) t[1]^2 * same_school + sampling,
        JAGS_estimates_table(two_level[["fit"]], return_samples = TRUE)[, "(mu) sd(intercept)", drop = FALSE],
        sd_prior_sd = 0.5
      )[["logml"]]
    )
  })

  ### Bayes factors ----
  # the district level: three levels against schools only, as log BF10. No
  # draw of the spike-and-slab fit excludes the district SD, so its inclusion
  # BF is printed as the draw-count limit; it is not a valid lower bound here:
  # at the exact BF of about 1850 roughly 20 of the 40000 draws would exclude
  # the SD if the inclusion indicator mixed
  scenario_text("bayes-factor-district", {
    district_spike <- JAGS_inference_table(fit_district_spike)[["inclusion_BF"]]
    spike_bound    <- identical(attr(district_spike, "bound_operator"), ">")
    scenario_bf_table(
      H1 = "(1 | district/school)",
      H0 = "(1 | district:school)",
      log_BF10 = c(
        "exact"                  = exact_logml[["three_level"]] - exact_logml[["two_level"]],
        "bridge sampling"        = scenario_log_bf_bridge(three_level[["marglik"]], two_level[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_three_level[["marglik"]], brms_two_level[["marglik"]]),
        "spike and slab"         = log(as.numeric(district_spike))
      ),
      note = c("", "", "", if (spike_bound) "draw-count limit: no draw excludes the district SD" else "")
    )
  })
  scenario_text("inference-district-spike", JAGS_inference_table(fit_district_spike))

  ### pooled effect and a new school in a new district ----
  # metafor: the pooled estimate with its SE, and with the SD of a new
  # school's true effect (SE and both heterogeneity variances); BayesTools: the
  # fixed part and the "sample" policy for new district and school levels;
  # brms: re_formula = NA and sample_new_levels = "gaussian"
  scenario_text("predictions-new-school", {
    targets  <- data.frame(target = c("pooled effect", "new school in a new district"))
    new_data <- data.frame(district = factor("new"), school = factor("new"), sei = 0)
    pooled   <- JAGS_predict_formula(fit_three_level, "mu", data = new_data, formula_target = "fixed")$value
    new      <- JAGS_predict_formula(fit_three_level, "mu", data = new_data, formula_target = "conditional",
                                     new_levels = "sample", seed = 1)$value
    brms_pooled <- brms::posterior_epred(brms_three_level[["fit"]], newdata = new_data, re_formula = NA)
    brms_new    <- brms::posterior_epred(brms_three_level[["fit"]], newdata = new_data,
                                         allow_new_levels = TRUE, sample_new_levels = "gaussian")
    metafor_se <- metafor_three$se
    scenario_prediction_table(targets, list(
      "metafor"    = rep(as.numeric(stats::coef(metafor_three)), 2L),
      "metafor sd" = c(metafor_se, sqrt(metafor_se^2 + sum(metafor_three$sigma2))),
      "BayesTools" = rbind(pooled, new),
      "brms"       = rbind(t(brms_pooled), t(brms_new))
    ))
  })

  ### district and school effects ----
  # a level's effect is its conditional prediction with only its block minus
  # the fixed part; metafor BLUPs (ML) and brms posterior means as references
  level_effects <- function(block, groups) {
    prediction <- JAGS_predict_formula(fit_three_level, "mu", formula_target = "conditional", blocks = block, components = TRUE)
    tapply(rowMeans(prediction$value - prediction$components$fixed), groups, mean)
  }
  district_groups <- as.character(dat$district)
  school_groups   <- paste(dat$district, dat$school, sep = "/")
  scenario_plot("effects-comparison", {
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    blups  <- metafor::ranef(metafor_three)
    brms_effects <- brms::ranef(brms_three_level[["fit"]])
    district <- level_effects("district", district_groups)
    school   <- level_effects("school_district", school_groups)
    scenario_agreement_plot(blups[["district"]][names(district), "intrcpt"], unname(district),
                            "District effects", reference_label = "metafor")
    scenario_agreement_plot(blups[["district/school"]][names(school), "intrcpt"], unname(school),
                            "School effects", reference_label = "metafor")
    scenario_agreement_plot(brms_effects[["district"]][names(district), "Estimate", "Intercept"], unname(district),
                            "District effects", reference_label = "brms")
    scenario_agreement_plot(brms_effects[["district:school"]][sub("/", "_", names(school)), "Estimate", "Intercept"],
                            unname(school), "School effects", reference_label = "brms")
  })

  ### prior and posterior distributions ----
  scenario_plot("posterior-allocation", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    for (name in c("(mu) heterogeneity: sd_total", "(mu) heterogeneity: var_prop(district)")) {
      selection <- parameter_catalog_resolve(parameter_catalog(fit_allocation), name)
      plot_posterior(stats::setNames(list(parameter_mixed_posterior(fit_allocation, selection)), name), name,
                     prior = TRUE, main = name)
    }
  })
})
