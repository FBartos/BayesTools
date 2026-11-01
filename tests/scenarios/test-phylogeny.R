if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("phylogeny")
# testthat::test_file("tests/scenarios/test-phylogeny.R")

### Description
# a phylogenetic meta-analysis of the metadat dat.lim2014 offspring-size
# correlations (o_o_unadj: 170 Fisher-z correlations from 120 species in 125
# articles) with known sampling variances: random article and species
# intercepts and phylogenetic species effects whose correlation is known from
# the species tree (Grafen branch lengths) through a known group covariance,
# compared with metafor::rma.mv (ML) and with brms gr(cov = ) under matched
# priors; Bayes factors for the phylogenetic effects (bridge sampling in
# BayesTools and brms, and Savage-Dickey), and species effects against the
# metafor BLUPs and brms

testthat::test_that("Lim 2014 phylogenetic meta-analysis", {

  dat  <- as.data.frame(metafor::escalc(measure = "ZCOR", ri = ri, ni = ni, data = metadat::dat.lim2014[["o_o_unadj"]]))
  dat$article     <- factor(dat$article)
  dat$species_phy <- dat$species
  dat$sei         <- sqrt(dat$vi)
  phylogeny <- ape::vcv(ape::compute.brlen(metadat::dat.lim2014[["o_o_unadj_tree"]]), corr = TRUE)

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
  # priors on the Fisher-z scale: pooled effect ~ N(0, 1), SDs ~ half-normal(0, 0.5)
  unit_prior <- prior("normal", list(mean = 0, sd = 1))
  sd_prior   <- prior("normal", list(mean = 0, sd = 0.5), truncation = list(0, Inf))

  ### metafor fits ----
  metafor_phylo <- metafor::rma.mv(yi, vi, random = list(~ 1 | article, ~ 1 | species, ~ 1 | species_phy),
                                   R = list(species_phy = phylogeny), data = dat, method = "ML")
  metafor_no_phylo <- metafor::rma.mv(yi, vi, random = list(~ 1 | article, ~ 1 | species), data = dat, method = "ML")

  ### BayesTools fits ----
  # the phylogenetic effects share the tree's species correlation; their mean
  # trades off with the pooled effect, so the chains are long and thinned
  phylo <- scenario_fit("fit_phylo", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = random_effects_formula(
        ~ 1 + (1 | article) + (1 | species) + (1 | species_phy),
        group_covariance = list(species_phy = random_group_covariance(phylogeny))
      )),
      formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(sd = sd_prior)),
      chains = 4, adapt = 1000, burnin = 2000, sample = 10000, thin = 5, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  no_phylo <- scenario_fit("fit_no_phylo", {
    fit <- JAGS_fit(
      model_syntax = likelihood, data = data_list,
      formula_list = list(mu = ~ 1 + (1 | article) + (1 | species)),
      formula_data_list = list(mu = dat),
      formula_prior_list = list(mu = list(intercept = unit_prior)),
      formula_random_prior_list = list(mu = prior_random(sd = sd_prior)),
      chains = 4, adapt = 1000, burnin = 2000, sample = 10000, thin = 5, seed = 1
    )
    list(fit = fit, marglik = JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data_list, seed = 1))
  })
  fit_phylo <- phylo[["fit"]]

  ### brms fits ----
  brms_phylo <- scenario_fit("brms_phylo", {
    fit <- brms::brm(
      yi | se(sei) ~ 0 + Intercept + (1 | article) + (1 | species) + (1 | gr(species_phy, cov = phylogeny)),
      data = dat, data2 = list(phylogeny = phylogeny),
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 0.5), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12)
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })
  brms_no_phylo <- scenario_fit("brms_no_phylo", {
    fit <- brms::brm(
      yi | se(sei) ~ 0 + Intercept + (1 | article) + (1 | species), data = dat,
      prior = c(brms::prior(normal(0, 1), class = b), brms::prior(normal(0, 0.5), class = sd)),
      save_pars = brms::save_pars(all = TRUE), backend = "rstan",
      chains = 4, cores = 4, iter = 6000, warmup = 1000, seed = 1, refresh = 0,
      control = list(adapt_delta = 0.99)
    )
    set.seed(1)
    list(fit = fit, marglik = brms::bridge_sampler(fit, silent = TRUE))
  })

  ### model summaries ----
  scenario_text("estimates-phylo",    JAGS_estimates_table(fit_phylo))
  scenario_text("estimates-no-phylo", JAGS_estimates_table(no_phylo[["fit"]]))

  ### metafor and brms comparison ----
  # metafor's ML estimate of the species SD is 0 (a boundary)
  scenario_text("comparison-phylo", {
    metafor_sd <- sqrt(metafor_phylo$sigma2)
    metafor_values <- c(intercept = as.numeric(stats::coef(metafor_phylo)), sd_article = metafor_sd[[1L]],
                        sd_species = metafor_sd[[2L]], sd_phylogeny = metafor_sd[[3L]])
    scenario_comparison_table(
      names(metafor_values), references = list(metafor = metafor_values),
      posteriors = list(
        BayesTools = scenario_draws(JAGS_estimates_table(fit_phylo, return_samples = TRUE), c(
          intercept = "(mu) intercept", sd_article = "(mu) article: sd(intercept)",
          sd_species = "(mu) species: sd(intercept)", sd_phylogeny = "(mu) species_phy: sd(intercept)")),
        brms       = scenario_draws(as.matrix(brms_phylo[["fit"]]), c(
          intercept = "b_Intercept", sd_article = "sd_article__Intercept",
          sd_species = "sd_species__Intercept", sd_phylogeny = "sd_species_phy__Intercept"))
      )
    )
  })

  ### exact marginal likelihoods ----
  # with the pooled effect (N(0, 1)) and all random effects integrated out,
  # the marginal likelihoods are integrals over the half-normal(0, 0.5) SDs,
  # which importance sampling evaluates with an error far below that of bridge
  # sampling over the ~360 latent effects
  exact_logml <- scenario_fit("exact_marglik", {
    same_article <- outer(dat$article, dat$article, "==") * 1
    same_species <- outer(dat$species, dat$species, "==") * 1
    phylo_rows   <- phylogeny[dat$species_phy, dat$species_phy]
    sampling     <- diag(dat$vi)
    pooled       <- matrix(1, nrow(dat))
    sd_columns   <- c("(mu) article: sd(intercept)", "(mu) species: sd(intercept)")
    set.seed(1)
    c(
      phylo = scenario_gaussian_logml(
        dat$yi, pooled, function(t) t[1]^2 * same_article + t[2]^2 * same_species + t[3]^2 * phylo_rows + sampling,
        JAGS_estimates_table(fit_phylo, return_samples = TRUE)[, c(sd_columns, "(mu) species_phy: sd(intercept)")],
        sd_prior_sd = 0.5
      )[["logml"]],
      no_phylo = scenario_gaussian_logml(
        dat$yi, pooled, function(t) t[1]^2 * same_article + t[2]^2 * same_species + sampling,
        JAGS_estimates_table(no_phylo[["fit"]], return_samples = TRUE)[, sd_columns],
        sd_prior_sd = 0.5
      )[["logml"]]
    )
  })

  ### Bayes factors ----
  # the phylogenetic effects as log BF10: the exact reference, bridge sampling
  # against the model without them (whose brms estimate over the latent
  # effects varies by several tenths between refits), and Savage-Dickey for
  # their SD at 0 within the full fit
  scenario_text("bayes-factor-phylogeny", {
    scenario_bf_table(
      H1 = "phylogenetic species effects",
      H0 = "no phylogenetic effects",
      log_BF10 = c(
        "exact"                  = exact_logml[["phylo"]] - exact_logml[["no_phylo"]],
        "bridge sampling"        = scenario_log_bf_bridge(phylo[["marglik"]], no_phylo[["marglik"]]),
        "bridge sampling (brms)" = scenario_log_bf_bridge(brms_phylo[["marglik"]], brms_no_phylo[["marglik"]]),
        "Savage-Dickey"          = log(as.numeric(hypothesis_BF(fit_phylo, "species_phy: sd(intercept) = 0", seed = 1)[["BF"]]))
      )
    )
  })
  scenario_text("likelihood-ratio-metafor", {
    data.frame(model = c("phylogenetic", "no phylogeny"),
               logLik = round(c(as.numeric(stats::logLik(metafor_phylo)), as.numeric(stats::logLik(metafor_no_phylo))), 3))
  })

  ### species effects ----
  # a species' phylogenetic effect is its conditional prediction with only the
  # phylogenetic block minus the fixed part; metafor BLUPs and brms posterior
  # means as references. The effects trade off with the pooled effect (ESS of
  # about 300), so they share its Monte Carlo error: against brms they shift
  # together by the difference in the pooled effects
  scenario_plot("species-effects", {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    prediction <- JAGS_predict_formula(fit_phylo, "mu", formula_target = "conditional", blocks = "species_phy", components = TRUE)
    effects    <- tapply(rowMeans(prediction$value - prediction$components$fixed), dat$species_phy, mean)
    blups      <- metafor::ranef(metafor_phylo)[["species_phy"]]
    brms_effects <- brms::ranef(brms_phylo[["fit"]])[["species_phy"]]
    scenario_agreement_plot(blups[names(effects), "intrcpt"], unname(effects), "Phylogenetic effects",
                            reference_label = "metafor")
    scenario_agreement_plot(brms_effects[names(effects), "Estimate", "Intercept"], unname(effects), "Phylogenetic effects",
                            reference_label = "brms")
  })

  ### prior and posterior distributions ----
  scenario_plot("posterior-sd", {
    par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    for (name in c("(mu) article: sd(intercept)", "(mu) species: sd(intercept)", "(mu) species_phy: sd(intercept)")) {
      selection <- parameter_catalog_resolve(parameter_catalog(fit_phylo), name)
      plot_posterior(stats::setNames(list(parameter_mixed_posterior(fit_phylo, selection)), name), name,
                     prior = TRUE, main = name)
    }
  })
})
