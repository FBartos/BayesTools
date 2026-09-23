skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Hypothesis Bayes Factor Bridge Comparisons
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()


.hypothesis_treatment_level_marginal_for_test <- function(fit){

  model <- coda::as.mcmc(as.matrix(fit[["mcmc"]]))
  class(model) <- c("BayesTools_fit", class(model))
  attr(model, "prior_list") <- attr(fit, "prior_list", exact = TRUE)

  mixed <- as_mixed_posteriors(
    model,
    parameters      = "mu_x_fac2t",
    n_prior_samples = 20000
  )

  marginal_posterior(
    samples       = mixed,
    parameter     = "mu_x_fac2t",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 20000
  )
}


# Log posterior odds of `parameter > 0` in the fitted draws and their
# delta-method Monte Carlo standard error, 1 / sqrt(ESS p (1 - p)), from the
# effective sample size of the region indicator across chains. For the
# treatment contrast, level B minus level A is the coefficient itself.
.hypothesis_draws_log_odds_for_test <- function(fit, parameter){

  indicator <- coda::mcmc.list(lapply(fit[["mcmc"]], function(chain){
    coda::mcmc(as.numeric(chain[, parameter] > 0))
  }))
  p   <- mean(unlist(indicator))
  ess <- unname(coda::effectiveSize(indicator))

  list(
    log_odds = log(p / (1 - p)),
    mcse     = 1 / sqrt(ess * p * (1 - p))
  )
}


# Checks an inequality hypothesis BF against the ratio of the truncated-prior
# bridge marginal likelihoods (the encompassing-prior identity; the symmetric
# N(0, 1) treatment prior has prior odds 1). The fixture's log BF is small
# (0.22), so the direction is asserted explicitly: a reciprocal BF or BF = 1
# fails. The magnitude must agree within 3 Monte Carlo standard errors of the
# difference: the bridge MCSEs of both marginal likelihoods and the posterior
# draw error of the encompassing fit (together about 0.09 on the log scale).
# The hypothesis BF must also reproduce the log posterior odds of the same
# draws within 0.08, four standard errors of resampling 20000 posterior and
# 20000 prior draws (sqrt(1 / (20000 p (1 - p)) + 1 / (20000 / 4)) = 0.02).
.expect_inequality_BF_matches_bridge <- function(bf_hyp, fit, ml_positive, ml_negative){

  log_BF_hypothesis <- log(attr(bf_hyp, "raw_BF"))
  log_BF_bridge     <- ml_positive[["logml"]] - ml_negative[["logml"]]
  draws             <- .hypothesis_draws_log_odds_for_test(fit, "mu_x_fac2t")
  bridge_mcse       <- sqrt(sum(c(
    ml_positive[["repetitions"]][["mcse"]],
    ml_negative[["repetitions"]][["mcse"]]
  )^2))
  difference_mcse   <- sqrt(draws$mcse^2 + bridge_mcse^2)

  expect_gt(log_BF_bridge, 0)
  expect_gt(log_BF_hypothesis, 0)
  expect_lt(abs(log_BF_hypothesis - log_BF_bridge), 3 * difference_mcse)
  expect_lt(abs(log_BF_hypothesis - draws$log_odds), 0.08)
}


test_that("hypothesis_BF point-null agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  fit_null <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  post <- as.matrix(fit_alt[["mcmc"]])[, "m"]
  prior <- attr(fit_alt, "prior_list")[["m"]]

  bf_hyp <- hypothesis_BF(
    posterior  = post,
    prior      = prior,
    hypothesis = "m == 0"
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.25)
  expect_s3_class(fit_null, "BayesTools_fit")
})


test_that("hypothesis_BF explicit level point-null agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] = 0",
    seed       = 101
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.35)
})


test_that("hypothesis_BF treatment-level point contrast agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] - mu_x_fac2t[A] = 0",
    seed       = 102
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.45)
})


test_that("hypothesis_BF treatment-level inequality agrees with truncated-prior bridge BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))
  ml_negative <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_negative.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] > mu_x_fac2t[A] vs mu_x_fac2t[B] < mu_x_fac2t[A]",
    seed       = 103
  )

  .expect_inequality_BF_matches_bridge(bf_hyp, fit_alt, ml_positive, ml_negative)
})


test_that("hypothesis_BF transformed level inequality agrees with truncated-prior bridge BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))
  ml_negative <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_negative.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = paste(
      "exp(mu_x_fac2t[B]) > exp(mu_x_fac2t[A])",
      "vs",
      "exp(mu_x_fac2t[B]) < exp(mu_x_fac2t[A])"
    ),
    seed       = 104
  )

  .expect_inequality_BF_matches_bridge(bf_hyp, fit_alt, ml_positive, ml_negative)
})


test_that("hypothesis_BF point-vs-level-region agrees with truncated-prior bridge BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = paste(
      "mu_x_fac2t[B] - mu_x_fac2t[A] = 0",
      "vs",
      "mu_x_fac2t[B] - mu_x_fac2t[A] > 0"
    ),
    seed       = 105
  )
  bf_bridge <- exp(ml_null[["logml"]] - ml_positive[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.55)
})


test_that("hypothesis_BF nested coefficient agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  post <- as.matrix(fit_alt[["mcmc"]])[, "mu_x_fac2t"]
  prior <- attr(fit_alt, "prior_list")[["mu_x_fac2t"]]

  bf_hyp <- hypothesis_BF(
    posterior  = post,
    prior      = prior,
    hypothesis = "mu_x_fac2t == 0"
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.35)
})
