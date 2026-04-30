# TEST FILE: Weightfunction prior redesign
# ============================================================================ #
#
# PURPOSE:
#   Contract tests for the unified prior_weightfunction() API, independent
#   weight priors, one/two-sided mapping, JAGS code generation, and bridge
#   sampling parameter extraction.
#
# TAGS: @priors, @jags, @weightfunctions
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("prior_weightfunction stores canonical geometry and weight priors", {

  wf <- prior_weightfunction(
    side = "one-sided",
    steps = c(.025, .05),
    weights = wf_cumulative(c(1, 2, 3)),
    prior_weights = 2
  )

  expect_true(is.prior.weightfunction(wf))
  expect_equal(wf$distribution, "weightfunction")
  expect_equal(wf$side, "one-sided")
  expect_equal(wf$steps, c(.025, .05))
  expect_equal(wf$bins$lower, c(0, .025, .05))
  expect_equal(wf$bins$upper, c(.025, .05, 1))
  expect_equal(wf$bins$reference, c(TRUE, FALSE, FALSE))
  expect_equal(wf$weights$type, "cumulative")
  expect_equal(wf$weights$alpha, c(1, 2, 3))
  expect_equal(wf$prior_weights, 2)

  samples <- rng(wf, 1000)
  expect_true(all(samples[,1] == 1))
})

test_that("weightfunction constructors validate independent scales", {

  expect_silent(wf_independent(prior("beta", list(1, 1))))
  expect_silent(wf_independent(prior("gamma", list(2, 1))))
  expect_error(
    wf_independent(prior("normal", list(0, 1))),
    "non-negative support"
  )

  expect_silent(wf_independent(
    prior("normal", list(0, 1)),
    scale = "log_omega"
  ))

  expect_error(
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(.9, .5))),
    "reference-bin"
  )
})

test_that("weightfunctions_mapping expands two-sided priors onto one-sided cuts", {

  one_sided <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_cumulative(c(1, 1, 1))
  )
  two_sided <- prior_weightfunction(
    "two-sided",
    c(.05),
    wf_fixed(c(1, .5))
  )

  expect_equal(
    weightfunctions_mapping(list(one_sided, two_sided), cuts_only = TRUE, one_sided = TRUE),
    c(0, .025, .05, .975, 1)
  )
  expect_equal(
    weightfunctions_mapping(list(one_sided, two_sided), one_sided = TRUE),
    list(c(1L, 2L, 3L, 3L), c(1L, 2L, 2L, 1L))
  )
})

test_that("JAGS generation uses component-local omega for bias mixtures", {

  fixed <- prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .5)))
  independent <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_independent(prior("beta", list(2, 3)))
  )
  bias <- prior_mixture(list(prior_none(), fixed, independent))

  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_match(syntax, "omega_component_1\\[1\\] <- 1")
  expect_match(syntax, "omega_local_component_2\\[2\\] <- 0.5")
  expect_match(syntax, "omega_local_component_3\\[2\\] ~ dbeta\\(2,3\\)")
  expect_match(syntax, "omega\\[1\\] <- omega_component_1\\[1\\] \\* equals\\(bias_indicator, 1\\)")
  expect_false(grepl("eta2omega", syntax, fixed = TRUE))
})

test_that("JAGS bridge helpers use natural latent weight parameters", {

  cumulative <- prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)))
  independent <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_independent(prior("beta", list(2, 3)))
  )
  independent_gamma <- prior_weightfunction(
    "one-sided",
    c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  log_independent <- prior_weightfunction(
    "one-sided",
    c(.05),
    wf_independent(
      prior("normal", list(0, 1)),
      "log_omega"
    )
  )

  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(cumulative)), paste0("eta[", 1:3, "]"))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(independent)), paste0("omega[", 2:3, "]"))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(independent_gamma)), "omega[2]")
  expect_equal(attr(.JAGS_bridgesampling_posterior_info.weightfunction(independent_gamma), "ub"), c("omega[2]" = Inf))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(log_independent)), "log_omega[2]")
  expect_equal(attr(.JAGS_bridgesampling_posterior_info.weightfunction(log_independent), "ub"), c("log_omega[2]" = Inf))

  samples <- c("eta[1]" = 1, "eta[2]" = 2, "eta[3]" = 3, "omega[2]" = 1.2, "omega[3]" = .1, "log_omega[2]" = .5)
  expect_equal(JAGS_marglik_parameters(samples, list(omega = cumulative))$omega, c(1, 5/6, 1/2))
  expect_equal(JAGS_marglik_parameters(samples, list(omega = independent))$omega, c(1, 1.2, .1))
  expect_equal(JAGS_marglik_parameters(samples, list(omega = independent_gamma))$omega, c(1, 1.2))
  expect_equal(
    JAGS_marglik_priors(samples, list(omega = independent_gamma)),
    mlpdf(independent_gamma$weights$prior, 1.2),
    tolerance = 1e-12
  )
  expect_equal(JAGS_marglik_parameters(samples, list(omega = log_independent))$omega, c(1, exp(.5)))
  expect_equal(
    JAGS_marglik_priors(samples, list(omega = log_independent)),
    mlpdf(log_independent$weights$prior, .5),
    tolerance = 1e-12
  )
})

test_that("heterogeneous bias mixtures map cumulative, omega, log-omega, fixed, and null weightfunctions", {

  bias <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))

  cuts <- weightfunctions_mapping(bias[sapply(bias, is.prior.weightfunction)], cuts_only = TRUE, one_sided = TRUE)
  mapping <- weightfunctions_mapping(bias[sapply(bias, is.prior.weightfunction)], one_sided = TRUE)

  expect_equal(cuts, c(0, .025, .05, .10, .975, 1))
  expect_equal(mapping[[1]], c(1L, 2L, 3L, 3L, 3L))
  expect_equal(mapping[[2]], c(1L, 1L, 2L, 3L, 3L))
  expect_equal(mapping[[3]], c(1L, 2L, 2L, 2L, 2L))
  expect_equal(mapping[[4]], c(1L, 2L, 2L, 2L, 1L))

  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_match(syntax, "eta_component_2\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax, "omega_local_component_3\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(syntax, "log_omega_component_4\\[2\\] ~ dnorm")
  expect_match(syntax, "omega_local_component_4\\[2\\] <- exp\\(log_omega_component_4\\[2\\]\\)")
  expect_match(syntax, "omega_local_component_5\\[2\\] <- 0.4")
  expect_match(syntax, "omega_component_5\\[5\\] <- omega_local_component_5\\[1\\]")
  expect_match(syntax, "omega\\[3\\] <- omega_component_1\\[3\\] \\* equals\\(bias_indicator, 1\\)")
  expect_false(grepl("eta2omega", syntax, fixed = TRUE))

  prior_samples <- .mix_priors.weightfunction(
    as.list(bias)[c(1, which(sapply(bias, is.prior.weightfunction)))],
    parameter = "omega",
    seed = 13,
    n_samples = 600
  )

  expect_equal(colnames(prior_samples), c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.1]", "omega[0.1,0.975]", "omega[0.975,1]"))
  expect_true(all(prior_samples[attr(prior_samples, "models_ind") == 1, ] == 1))
  expect_true(all(prior_samples[attr(prior_samples, "models_ind") == 5, "omega[0.025,0.05]"] == .4))
  expect_true(all(prior_samples[attr(prior_samples, "models_ind") == 5, "omega[0.975,1]"] == 1))
  expect_gt(mean(prior_samples[attr(prior_samples, "models_ind") == 3, "omega[0.05,0.1]"] > 1), .90)
  expect_gt(mean(prior_samples[attr(prior_samples, "models_ind") == 4, "omega[0.025,0.05]"] > 1), .95)
})

test_that("JAGS syntax and fitting allow independent omega weights above one", {

  skip_if_not_installed("rjags")
  skip_if_missing_fits(c("fit_wf_independent_gamma", "fit_wf_independent_log"))

  omega_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega")
  )

  omega_syntax <- JAGS_add_priors("model{}", list(omega = omega_prior))
  log_syntax   <- JAGS_add_priors("model{}", list(omega = log_prior))

  expect_match(omega_syntax, "omega\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(log_syntax, "log_omega\\[2\\] ~ dnorm")
  expect_match(log_syntax, "omega\\[2\\] <- exp\\(log_omega\\[2\\]\\)")

  omega_fit <- readRDS(file.path(temp_fits_dir, "fit_wf_independent_gamma.RDS"))
  log_fit   <- readRDS(file.path(temp_fits_dir, "fit_wf_independent_log.RDS"))

  omega_samples <- as.matrix(.fit_to_posterior(omega_fit))
  log_samples   <- as.matrix(.fit_to_posterior(log_fit))

  expect_true("omega[2]" %in% colnames(omega_samples))
  expect_true("omega[2]" %in% colnames(log_samples))
  expect_true("log_omega[2]" %in% colnames(log_samples))
  expect_gt(mean(omega_samples[, "omega[2]"] > 1), .90)
  expect_gt(mean(log_samples[, "omega[2]"] > 1), .90)
  expect_equal(
    unname(log_samples[, "omega[2]"]),
    unname(exp(log_samples[, "log_omega[2]"])),
    tolerance = 1e-8
  )
})

test_that("JAGS fits heterogeneous bias mixtures with omega and log-omega weights above one", {

  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_bias_heterogeneous_wf")

  fit <- readRDS(file.path(temp_fits_dir, "fit_bias_heterogeneous_wf.RDS"))

  posterior <- as.matrix(.fit_to_posterior(fit))
  expect_true(all(c("bias_indicator", paste0("omega[", 1:5, "]")) %in% colnames(posterior)))
  expect_true(all(1:5 %in% posterior[, "bias_indicator"]))

  indicator <- posterior[, "bias_indicator"]
  expect_true(all(posterior[indicator == 1, paste0("omega[", 1:5, "]")] == 1))
  expect_true(all(posterior[indicator == 2, paste0("omega[", 1:5, "]")] <= 1))
  expect_gt(mean(posterior[indicator == 3, "omega[3]"] > 1), .90)
  expect_gt(mean(posterior[indicator == 4, "omega[2]"] > 1), .95)
  expect_true(all(abs(posterior[indicator == 5, "omega[2]"] - .4) < 1e-8))
  expect_true(all(abs(posterior[indicator == 5, "omega[5]"] - 1) < 1e-8))

  mixed <- as_mixed_posteriors(fit, parameters = "bias", conditional = "omega")
  expect_equal(colnames(mixed$bias), c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.1]", "omega[0.1,0.975]", "omega[0.975,1]"))
  expect_true(all(attr(mixed$bias, "models_ind") %in% 2:5))
  expect_gt(mean(mixed$bias[attr(mixed$bias, "models_ind") == 3, "omega[0.05,0.1]"] > 1), .90)
  expect_gt(mean(mixed$bias[attr(mixed$bias, "models_ind") == 4, "omega[0.025,0.05]"] > 1), .95)
})

test_that("JAGS fits full bias mixtures with PET, PEESE, and heterogeneous weightfunctions", {

  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_bias_petpeese_heterogeneous_wf")

  bias <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, .4), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_PEESE("gamma", list(shape = 3, rate = 2), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))

  omega_names <- c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.1]", "omega[0.1,0.975]", "omega[0.975,1]")

  syntax <- JAGS_add_priors("model{}", list(bias = bias))
  expect_match(syntax, "PET <- PET_1 \\* equals\\(bias_indicator, 2\\)")
  expect_match(syntax, "eta_component_3\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax, "omega_local_component_4\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(syntax, "PEESE <- PEESE_1 \\* equals\\(bias_indicator, 5\\)")
  expect_match(syntax, "log_omega_component_6\\[2\\] ~ dnorm")
  expect_match(syntax, "omega_local_component_7\\[2\\] <- 0.4")
  expect_match(syntax, "omega_component_2\\[3\\] <- 1")
  expect_match(syntax, "omega_component_5\\[3\\] <- 1")
  expect_match(syntax, "omega\\[3\\] <- omega_component_1\\[3\\] \\* equals\\(bias_indicator, 1\\) \\+ omega_component_2\\[3\\] \\* equals\\(bias_indicator, 2\\)")

  fit <- readRDS(file.path(temp_fits_dir, "fit_bias_petpeese_heterogeneous_wf.RDS"))

  posterior <- as.matrix(.fit_to_posterior(fit))
  indicator <- posterior[, "bias_indicator"]
  expect_true(all(c("bias_indicator", paste0("omega[", 1:5, "]"), "PET", "PEESE") %in% colnames(posterior)))
  expect_true(all(1:7 %in% indicator))

  expect_true(all(posterior[indicator %in% c(1, 2, 5), paste0("omega[", 1:5, "]")] == 1))
  expect_true(all(posterior[indicator == 3, paste0("omega[", 1:5, "]")] <= 1))
  expect_gt(mean(posterior[indicator == 4, "omega[3]"] > 1), .90)
  expect_gt(mean(posterior[indicator == 6, "omega[2]"] > 1), .95)
  expect_true(all(abs(posterior[indicator == 7, "omega[2]"] - .4) < 1e-8))
  expect_true(all(abs(posterior[indicator == 7, "omega[5]"] - 1) < 1e-8))

  expect_true(all(abs(posterior[indicator != 2, "PET"]) < 1e-8))
  expect_true(any(posterior[indicator == 2, "PET"] > 0))
  expect_true(all(abs(posterior[indicator != 5, "PEESE"]) < 1e-8))
  expect_true(any(posterior[indicator == 5, "PEESE"] > 0))

  mixed_all <- as_mixed_posteriors(fit, parameters = "bias")
  expect_equal(colnames(mixed_all$bias), c(omega_names, "PET", "PEESE"))
  expect_true(all(1:7 %in% attr(mixed_all$bias, "models_ind")))

  mixed_omega <- as_mixed_posteriors(fit, parameters = "bias", conditional = "omega")
  expect_equal(colnames(mixed_omega$bias), omega_names)
  expect_true(all(attr(mixed_omega$bias, "models_ind") %in% c(3, 4, 6, 7)))
  expect_gt(mean(mixed_omega$bias[attr(mixed_omega$bias, "models_ind") == 4, "omega[0.05,0.1]"] > 1), .90)
  expect_gt(mean(mixed_omega$bias[attr(mixed_omega$bias, "models_ind") == 6, "omega[0.025,0.05]"] > 1), .95)

  mixed_pet <- as_mixed_posteriors(fit, parameters = "bias", conditional = "PET")
  expect_equal(colnames(mixed_pet$bias), "PET")
  expect_true(all(attr(mixed_pet$bias, "models_ind") == 2))
  expect_true(all(mixed_pet$bias[, "PET"] > 0))

  mixed_peese <- as_mixed_posteriors(fit, parameters = "bias", conditional = "PEESE")
  expect_equal(colnames(mixed_peese$bias), "PEESE")
  expect_true(all(attr(mixed_peese$bias, "models_ind") == 5))
  expect_true(all(mixed_peese$bias[, "PEESE"] > 0))

  mixed_petpeese <- as_mixed_posteriors(fit, parameters = "bias", conditional = "PETPEESE")
  expect_equal(colnames(mixed_petpeese$bias), c("PET", "PEESE"))
  expect_true(all(attr(mixed_petpeese$bias, "models_ind") %in% c(2, 5)))
  expect_true(all(mixed_petpeese$bias[attr(mixed_petpeese$bias, "models_ind") == 2, "PEESE"] == 0))
  expect_true(all(mixed_petpeese$bias[attr(mixed_petpeese$bias, "models_ind") == 5, "PET"] == 0))

  table_samples <- suppressWarnings(runjags_estimates_table(
    fit,
    conditional        = TRUE,
    return_samples     = TRUE,
    remove_diagnostics = TRUE
  ))
  expect_true(all(c(omega_names, "PET", "PEESE") %in% colnames(table_samples)))
  expect_true(any(is.na(table_samples[, "PET"])))
  expect_true(any(table_samples[, "PET"] > 0, na.rm = TRUE))
  expect_true(any(is.na(table_samples[, "PEESE"])))
  expect_true(any(table_samples[, "PEESE"] > 0, na.rm = TRUE))
  expect_true(any(is.na(table_samples[, "omega[0.05,0.1]"])))
  expect_true(any(table_samples[, "omega[0.05,0.1]"] > 1, na.rm = TRUE))
})

test_that("point(1) weightfunction null components are handled explicitly", {

  wf <- prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1)), prior_weights = 3)
  point_null <- prior("point", list(1), prior_weights = 1)

  mixed <- .mix_priors.weightfunction(
    list(point_null, wf),
    parameter = "omega",
    seed = 1,
    n_samples = 40
  )

  expect_equal(colnames(mixed), c("omega[0,0.05]", "omega[0.05,1]"))
  expect_equal(unname(mixed[attr(mixed, "models_ind") == 1, ]), matrix(1, nrow = 10, ncol = 2))
  expect_error(
    .mix_priors.weightfunction(
      list(prior("point", list(.5), prior_weights = 1), wf),
      parameter = "omega",
      n_samples = 40
    ),
    "point\\(1\\)/none null priors"
  )
})

test_that("omega diagnostics reject bias mixtures without weightfunctions", {

  fit <- structure(list(), class = c("runjags", "BayesTools_fit"))
  attr(fit, "prior_list") <- list(
    bias = prior_mixture(list(
      prior_PET("normal", list(0, 1), prior_weights = 1),
      prior_PEESE("normal", list(0, 1), prior_weights = 1)
    ))
  )

  expect_error(
    JAGS_diagnostics_density(fit, parameter = "omega"),
    "at least one weightfunction component"
  )
})
