skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Hypothesis Bayes Factor Bridge Comparisons
# ============================================================================ #

source(testthat::test_path("common-functions.R"))


.hypothesis_bridge_cache_ready_for_test <- function(names){

  files <- c(
    file.path(temp_fits_dir, paste0(names, ".RDS")),
    file.path(temp_marglik_dir, paste0(names, ".RDS"))
  )
  if(all(file.exists(files))){
    return(TRUE)
  }

  testthat::fail(
    "Pre-fitted bridge fixtures not found. Run test-00-model-fits.R first."
  )
  return(FALSE)
}


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


test_that("hypothesis_BF point-null agrees with bridge-sampling model BF", {

  if(!.hypothesis_bridge_cache_ready_for_test(c("fit_simple_normal", "fit_simple_spike"))){
    return()
  }

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

  if(!.hypothesis_bridge_cache_ready_for_test(c("fit_formula_treatment", "fit_formula_simple"))){
    return()
  }

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

  if(!.hypothesis_bridge_cache_ready_for_test(c("fit_formula_treatment", "fit_formula_simple"))){
    return()
  }

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

  if(!.hypothesis_bridge_cache_ready_for_test(c(
    "fit_formula_treatment",
    "fit_formula_treatment_positive",
    "fit_formula_treatment_negative"
  ))){
    return()
  }

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))
  ml_negative <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_negative.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] > mu_x_fac2t[A] vs mu_x_fac2t[B] < mu_x_fac2t[A]",
    seed       = 103
  )
  bf_bridge <- exp(ml_positive[["logml"]] - ml_negative[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.45)
})


test_that("hypothesis_BF transformed level inequality agrees with truncated-prior bridge BF", {

  if(!.hypothesis_bridge_cache_ready_for_test(c(
    "fit_formula_treatment",
    "fit_formula_treatment_positive",
    "fit_formula_treatment_negative"
  ))){
    return()
  }

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
  bf_bridge <- exp(ml_positive[["logml"]] - ml_negative[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.45)
})


test_that("hypothesis_BF point-vs-level-region agrees with truncated-prior bridge BF", {

  if(!.hypothesis_bridge_cache_ready_for_test(c(
    "fit_formula_simple",
    "fit_formula_treatment",
    "fit_formula_treatment_positive"
  ))){
    return()
  }

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

  if(!.hypothesis_bridge_cache_ready_for_test(c("fit_formula_treatment", "fit_formula_simple"))){
    return()
  }

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
