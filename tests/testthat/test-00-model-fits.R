skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: Model Fits for Reuse Across Tests
# ============================================================================ #
#
# PURPOSE:
#   Centralized model fitting for all JAGS models used across the test suite.
#   Fitted models are saved to temp directory for reuse in other test files.
#   This reduces redundant MCMC sampling and speeds up the overall test suite.
#
# DEPENDENCIES:
#   - rjags, runjags, bridgesampling: For model fitting
#
# SKIP CONDITIONS:
#   - skip_on_cran(): Long-running model fitting
#   - skip_if_not_installed("rjags")
#
# MODELS/FIXTURES:
#   - Creates all pre-fitted models used by other test files
#   - Models saved to BAYESTOOLS_TEST_FITS_DIR environment variable
#   - Maintains model_registry.RDS with metadata
#
# TAGS: @slow, @JAGS, @model-fits
# ============================================================================ #

# This file contains all model fitting procedures used across the test suite.
# Fitted models are saved to a temporary directory for reuse in other tests.
# This reduces redundant MCMC sampling and speeds up the overall test suite.

skip_on_cran()
skip_if_not_installed("rjags")

# Load common test helpers
source(testthat::test_path("common-functions.R"))
fit_cache_catalog <- bayestools_required_fit_catalog()
skip_refit_if_cached(
  "model-fit",
  required_fits = fit_cache_catalog$model_name,
  required_margliks = fit_cache_catalog$model_name[fit_cache_catalog$has_marglik],
  registry_file = file.path(test_files_dir, "model_registry.RDS")
)
rm(fit_cache_catalog)

# Initialize model registry to track metadata about each fitted model
model_registry <- list()

# ============================================================================ #
# SECTION 1: SIMPLE PRIOR DISTRIBUTIONS
# ============================================================================ #
test_that("Simple prior models fit correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data <- list(
    x = rnorm(50, 0, .5),
    N = 50
  )

  # Model 1: Normal and truncated normal priors
  priors_simple_normal <- list(
    m = prior("normal", list(0, 1)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"

  fit_simple_normal <- JAGS_fit(model_syntax, data, priors_simple_normal,
                                chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)

  # Compute marginal likelihood for model averaging
  log_posterior_simple_normal <- function(parameters, data){
    sum(stats::dnorm(data$x, parameters[["m"]], parameters[["s"]], log = TRUE))
  }
  marglik_simple_normal <- JAGS_bridgesampling(fit_simple_normal,
                                               log_posterior = log_posterior_simple_normal,
                                               data = data, prior_list = priors_simple_normal)

  result <- save_fit(fit_simple_normal, "fit_simple_normal",
                     marglik = marglik_simple_normal,
                     simple_priors = TRUE,
                     note = "Normal and truncated normal priors with data")
  model_registry[["fit_simple_normal"]] <<- result$registry_entry
  fit_simple_normal <- result$fit

  # Model 2: Spike and normal priors (for model averaging)
  priors_simple_spike <- list(
    m = prior("spike", list(0)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )

  fit_simple_spike <- JAGS_fit(model_syntax, data, priors_simple_spike,
                               chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)

  # Compute marginal likelihood for model averaging
  marglik_simple_spike <- JAGS_bridgesampling(fit_simple_spike,
                                              log_posterior = log_posterior_simple_normal,
                                              data = data, prior_list = priors_simple_spike)

  result <- save_fit(fit_simple_spike, "fit_simple_spike",
                     marglik = marglik_simple_spike,
                     simple_priors = TRUE,
                     note = "Spike and truncated normal priors with data (for model averaging)")
  model_registry[["fit_simple_spike"]] <<- result$registry_entry
  fit_simple_spike <- result$fit

  # Model 3: Various prior distributions
  priors_various <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("t", list(0, .5, 5)),
    p4  = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    p5  = prior("gamma", list(2, 1)),
    p6  = prior("invgamma", list(3, 2), list(1, 3)),
    p7  = prior("exp", list(1.5)),
    p8  = prior("beta", list(3, 2)),
    p9  = prior("uniform", list(1, 5)),
    p10 = prior("point", list(1))
  )

  model_syntax_simple <- "model{}"

  fit_simple_various <- suppressWarnings(JAGS_fit(model_syntax_simple, data = NULL, prior_list = priors_various,
                                                  chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_simple_various, "fit_simple_various",
                     simple_priors = TRUE,
                     note = "Various univariate distributions: normal, lognormal, t, Cauchy, gamma, invgamma, exp, beta, uniform, point")
  model_registry[["fit_simple_various"]] <<- result$registry_entry
  fit_simple_various <- result$fit

  # Model 4: PET and PEESE priors
  priors_pub_bias <- list(
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
  )

  model_syntax_pb <- "model{}"

  fit_simple_pub_bias <- suppressWarnings(JAGS_fit(model_syntax_pb, data = NULL, prior_list = priors_pub_bias,
                                                   chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_simple_pub_bias, "fit_simple_pub_bias",
                     pub_bias_priors = TRUE,
                     note = "PET and PEESE priors for publication bias")
  model_registry[["fit_simple_pub_bias"]] <<- result$registry_entry
  fit_simple_pub_bias <- result$fit

  # Model 5: Test with thinning parameter
  priors_thin <- list(
    mu = prior("normal", list(0, 1))
  )
  model_syntax_thin <- "model{}"

  fit_simple_thin <- suppressWarnings(JAGS_fit(model_syntax_thin, data = NULL, prior_list = priors_thin,
                                               chains = 2, adapt = 100, burnin = 150, sample = 300, thin = 3, seed = 2))
  result <- save_fit(fit_simple_thin, "fit_simple_thin",
                     simple_priors = TRUE, thinning = TRUE,
                     note = "Simple normal prior with thinning parameter (thin=3)")
  model_registry[["fit_simple_thin"]] <<- result$registry_entry
  fit_simple_thin <- result$fit
})


# ============================================================================ #
# SECTION 1B: MODELS FOR SUMMARY TABLES TESTING
# ============================================================================ #
test_that("Summary tables models fit correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data_summary <- list(
    x = rnorm(20, 0, 1),
    N = 20
  )

  model_syntax_summary <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, 1)
      }
    }"

  # Log posterior for summary tables (constant, no data dependency)
  log_posterior_summary <- function(parameters, data){
    return(0)
  }

  # Model 1: Normal prior with prior_none weightfunction
  priors_summary0 <- list(
    m     = prior("normal", list(0, 1)),
    omega = prior_none()
  )

  fit_summary0 <- JAGS_fit(model_syntax_summary, data_summary, priors_summary0,
                           chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 0)

  marglik_summary0 <- JAGS_bridgesampling(fit_summary0,
                                          log_posterior = log_posterior_summary,
                                          data = data_summary, prior_list = priors_summary0)

  result <- save_fit(fit_summary0, "fit_summary0",
                     marglik = marglik_summary0,
                     simple_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Model for summary tables with no weightfunction")
  model_registry[["fit_summary0"]] <<- result$registry_entry
  fit_summary0 <- result$fit

  # Model 2: Normal prior with one-sided weightfunction (2 intervals)
  priors_summary1 <- list(
    m  = prior("normal", list(0, .5)),
    omega = prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))
  )

  fit_summary1 <- JAGS_fit(model_syntax_summary, data_summary, priors_summary1,
                           chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)

  marglik_summary1 <- JAGS_bridgesampling(fit_summary1,
                                          log_posterior = log_posterior_summary,
                                          data = data_summary, prior_list = priors_summary1)

  result <- save_fit(fit_summary1, "fit_summary1",
                     marglik = marglik_summary1,
                     simple_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Model for summary tables with one-sided weightfunction (cutpoint at .05)")
  model_registry[["fit_summary1"]] <<- result$registry_entry
  fit_summary1 <- result$fit

  # Model 3: Normal prior with one-sided weightfunction (3 intervals)
  priors_summary2 <- list(
    m  = prior("normal", list(0, .3)),
    omega = prior_weightfunction("one-sided", c(0.05, 0.50), wf_cumulative(c(1, 1, 1)))
  )

  fit_summary2 <- JAGS_fit(model_syntax_summary, data_summary, priors_summary2,
                           chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)

  marglik_summary2 <- JAGS_bridgesampling(fit_summary2,
                                          log_posterior = log_posterior_summary,
                                          data = data_summary, prior_list = priors_summary2)

  result <- save_fit(fit_summary2, "fit_summary2",
                     marglik = marglik_summary2,
                     simple_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Model for summary tables with one-sided weightfunction (cutpoints at .05, .50)")
  model_registry[["fit_summary2"]] <<- result$registry_entry
  fit_summary2 <- result$fit

  # Model 4: Normal prior with fixed weightfunction
  priors_summary3 <- list(
    m  = prior("normal", list(0, .3)),
    omega = prior_weightfunction("two-sided", 0.20, wf_fixed(c(1, .3)))
  )

  fit_summary3 <- JAGS_fit(model_syntax_summary, data_summary, priors_summary3,
                           chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)

  marglik_summary3 <- JAGS_bridgesampling(fit_summary3,
                                          log_posterior = log_posterior_summary,
                                          data = data_summary, prior_list = priors_summary3)

  result <- save_fit(fit_summary3, "fit_summary3",
                     marglik = marglik_summary3,
                     simple_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Model for summary tables with fixed weightfunction")
  model_registry[["fit_summary3"]] <<- result$registry_entry
  fit_summary3 <- result$fit
})


# ============================================================================ #
# SECTION 2: VECTOR PRIOR DISTRIBUTIONS
# ============================================================================ #
test_that("Vector prior models fit correctly", {

  skip_if_not_installed("rjags")

  # Multivariate normal
  priors_mnormal <- list(
    p1 = prior("mnormal", list(mean = 0, sd = 1, K = 3))
  )

  model_syntax_vec <- "model{}"

  fit_vector_mnormal <- suppressWarnings(JAGS_fit(model_syntax_vec, data = NULL, prior_list = priors_mnormal,
                                                  chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_vector_mnormal, "fit_vector_mnormal",
                     vector_priors = TRUE,
                     note = "Multivariate normal prior (K=3)")
  model_registry[["fit_vector_mnormal"]] <<- result$registry_entry
  fit_vector_mnormal <- result$fit

  # Multivariate cauchy
  priors_mcauchy <- list(
    p1 = prior("mcauchy", list(location = 0, scale = 1.5, K = 2))
  )

  model_syntax_mc <- "model{}"

  fit_vector_mcauchy <- suppressWarnings(JAGS_fit(model_syntax_mc, data = NULL, prior_list = priors_mcauchy,
                                                  chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2))
  result <- save_fit(fit_vector_mcauchy, "fit_vector_mcauchy",
                     vector_priors = TRUE,
                     note = "Multivariate Cauchy prior (K=2)")
  model_registry[["fit_vector_mcauchy"]] <<- result$registry_entry
  fit_vector_mcauchy <- result$fit

  # Multivariate t
  priors_mt <- list(
    p1 = prior("mt", list(location = 2, scale = 0.5, df = 5, K = 2))
  )

  model_syntax_mt <- "model{}"

  fit_vector_mt <- suppressWarnings(JAGS_fit(model_syntax_mt, data = NULL, prior_list = priors_mt,
                                             chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3))
  result <- save_fit(fit_vector_mt, "fit_vector_mt",
                     vector_priors = TRUE,
                     note = "Multivariate t prior with df=5 (K=2)")
  model_registry[["fit_vector_mt"]] <<- result$registry_entry
  fit_vector_mt <- result$fit
})


# ============================================================================ #

# SECTION 3: FACTOR PRIOR DISTRIBUTIONS
# ============================================================================ #
test_that("Factor prior models fit correctly", {

  skip_if_not_installed("rjags")

  # Orthonormal contrast
  priors_orthonormal <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_orthonormal[[1]], "levels") <- 3

  model_syntax_orth <- "model{}"

  fit_factor_orthonormal <- suppressWarnings(JAGS_fit(model_syntax_orth, data = NULL, prior_list = priors_orthonormal,
                                                      chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_factor_orthonormal, "fit_factor_orthonormal",
                     factor_priors = TRUE,
                     note = "Orthonormal contrast with 3 levels")
  model_registry[["fit_factor_orthonormal"]] <<- result$registry_entry
  fit_factor_orthonormal <- result$fit

  # Treatment contrast
  priors_treatment <- list(
    p1 = prior_factor("beta", list(alpha = 1, beta = 1), contrast = "treatment")
  )
  attr(priors_treatment[[1]], "levels") <- 2

  model_syntax_treat <- "model{}"

  fit_factor_treatment <- suppressWarnings(JAGS_fit(model_syntax_treat, data = NULL, prior_list = priors_treatment,
                                                    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2))
  result <- save_fit(fit_factor_treatment, "fit_factor_treatment",
                     factor_priors = TRUE,
                     note = "Treatment contrast with 2 levels and beta prior")
  model_registry[["fit_factor_treatment"]] <<- result$registry_entry
  fit_factor_treatment <- result$fit

  # Independent contrast
  priors_independent <- list(
    p1 = prior_factor("gamma", list(shape = 2, rate = 3), contrast = "independent")
  )
  attr(priors_independent[[1]], "levels") <- 3

  model_syntax_ind <- "model{}"

  fit_factor_independent <- suppressWarnings(JAGS_fit(model_syntax_ind, data = NULL, prior_list = priors_independent,
                                                      chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3))
  result <- save_fit(fit_factor_independent, "fit_factor_independent",
                     factor_priors = TRUE,
                     note = "Independent contrast with 3 levels and gamma prior")
  model_registry[["fit_factor_independent"]] <<- result$registry_entry
  fit_factor_independent <- result$fit

  # Meandif contrast
  priors_meandif <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 0.5), contrast = "meandif")
  )
  attr(priors_meandif[[1]], "levels") <- 3

  model_syntax_md <- "model{}"

  fit_factor_meandif <- suppressWarnings(JAGS_fit(model_syntax_md, data = NULL, prior_list = priors_meandif,
                                                  chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 4))
  result <- save_fit(fit_factor_meandif, "fit_factor_meandif",
                     factor_priors = TRUE,
                     note = "Meandif contrast with 3 levels")
  model_registry[["fit_factor_meandif"]] <<- result$registry_entry
  fit_factor_meandif <- result$fit
})


# ============================================================================ #
# SECTION 4: WEIGHTFUNCTION PRIORS
# ============================================================================ #
test_that("Weightfunction prior models fit correctly", {

  skip_if_not_installed("rjags")

  # One-sided weightfunction (2 intervals)
  priors_wf_onesided2 <- list(
    omega = prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1)))
  )

  model_syntax_wf1 <- "model{}"

  fit_weightfunction_onesided2 <- suppressWarnings(JAGS_fit(model_syntax_wf1, data = NULL, prior_list = priors_wf_onesided2,
                                                            chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_weightfunction_onesided2, "fit_weightfunction_onesided2",
                     weightfunction_priors = TRUE,
                     note = "One-sided weightfunction with 2 intervals (cutpoint at .05)")
  model_registry[["fit_weightfunction_onesided2"]] <<- result$registry_entry
  fit_weightfunction_onesided2 <- result$fit

  # One-sided weightfunction (3 intervals)
  priors_wf_onesided3 <- list(
    omega = prior_weightfunction("one-sided", c(.05, 0.10), wf_cumulative(c(1, 2, 3)))
  )

  model_syntax_wf2 <- "model{}"

  fit_weightfunction_onesided3 <- suppressWarnings(JAGS_fit(model_syntax_wf2, data = NULL, prior_list = priors_wf_onesided3,
                                                            chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2))
  result <- save_fit(fit_weightfunction_onesided3, "fit_weightfunction_onesided3",
                     weightfunction_priors = TRUE,
                     note = "One-sided weightfunction with 3 intervals (cutpoints at .05, .10)")
  model_registry[["fit_weightfunction_onesided3"]] <<- result$registry_entry
  fit_weightfunction_onesided3 <- result$fit

  # Two-sided weightfunction
  priors_wf_twosided <- list(
    omega = prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)))
  )

  model_syntax_wf3 <- "model{}"

  fit_weightfunction_twosided <- suppressWarnings(JAGS_fit(model_syntax_wf3, data = NULL, prior_list = priors_wf_twosided,
                                                           chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3))
  result <- save_fit(fit_weightfunction_twosided, "fit_weightfunction_twosided",
                     weightfunction_priors = TRUE,
                     note = "Two-sided weightfunction with cutpoint at .05")
  model_registry[["fit_weightfunction_twosided"]] <<- result$registry_entry
  fit_weightfunction_twosided <- result$fit

  # One-sided fixed weightfunction
  priors_wf_fixed <- list(
    omega = prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5)))
  )

  model_syntax_wf4 <- "model{}"

  fit_weightfunction_fixed <- suppressWarnings(JAGS_fit(model_syntax_wf4, data = NULL, prior_list = priors_wf_fixed,
                                                        chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 4))
  result <- save_fit(fit_weightfunction_fixed, "fit_weightfunction_fixed",
                     weightfunction_priors = TRUE,
                     note = "One-sided fixed weightfunction (weights: 1, .5)")
  model_registry[["fit_weightfunction_fixed"]] <<- result$registry_entry
  fit_weightfunction_fixed <- result$fit
})


# ============================================================================ #
# SECTION 4B: WEIGHTFUNCTION REDESIGN AND SELECTION KERNELS
# ============================================================================ #
test_that("Weightfunction redesign and selection-kernel models fit correctly", {

  skip_if_not_installed("rjags")

  # Independent gamma weightfunction
  omega_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  fit_wf_independent_gamma <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(omega = omega_prior),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 300,
    seed       = 11
  ))
  result <- save_fit(fit_wf_independent_gamma, "fit_wf_independent_gamma",
                     weightfunction_priors = TRUE,
                     note = "Independent gamma weightfunction for redesigned API tests")
  model_registry[["fit_wf_independent_gamma"]] <<- result$registry_entry
  fit_wf_independent_gamma <- result$fit

  # Independent log-omega weightfunction
  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega")
  )
  fit_wf_independent_log <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(omega = log_prior),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 300,
    seed       = 12
  ))
  result <- save_fit(fit_wf_independent_log, "fit_wf_independent_log",
                     weightfunction_priors = TRUE,
                     note = "Independent log-omega weightfunction for redesigned API tests")
  model_registry[["fit_wf_independent_log"]] <<- result$registry_entry
  fit_wf_independent_log <- result$fit

  # Heterogeneous bias mixture with cumulative, omega, log-omega, fixed, and null components
  bias_heterogeneous_wf <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))
  fit_bias_heterogeneous_wf <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(bias = bias_heterogeneous_wf),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 1000,
    seed       = 14
  ))
  result <- save_fit(fit_bias_heterogeneous_wf, "fit_bias_heterogeneous_wf",
                     pub_bias_priors = TRUE, weightfunction_priors = TRUE,
                     mixture_priors = TRUE,
                     note = "Heterogeneous bias mixture with cumulative, omega, log-omega, fixed, and null components")
  model_registry[["fit_bias_heterogeneous_wf"]] <<- result$registry_entry
  fit_bias_heterogeneous_wf <- result$fit

  # Full bias mixture with PET, PEESE, and heterogeneous weightfunctions
  bias_petpeese_heterogeneous_wf <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, .4), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_PEESE("gamma", list(shape = 3, rate = 2), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))
  fit_bias_petpeese_hetero_wf <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(bias = bias_petpeese_heterogeneous_wf),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 1200,
    seed       = 15
  ))
  result <- save_fit(fit_bias_petpeese_hetero_wf, "fit_bias_petpeese_hetero_wf",
                     pub_bias_priors = TRUE, weightfunction_priors = TRUE,
                     mixture_priors = TRUE,
                     note = "Full bias mixture with PET, PEESE, and heterogeneous weightfunctions")
  model_registry[["fit_bias_petpeese_hetero_wf"]] <<- result$registry_entry
  fit_bias_petpeese_hetero_wf <- result$fit

  # Ordinary mixture plus selection-kernel bias mixture for summary table tests
  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking  <- prior_phacking(form = "linear")
  bias_selection_kernel <- prior_mixture(list(prior_none(), phacking, prior_bias(selection, phacking)))
  mu_prior  <- prior_mixture(list(
    prior("point", list(0)),
    prior("normal", list(0, 1))
  ))
  fit_selection_kernel_summary <- suppressWarnings(JAGS_fit(
    "model{ y ~ dnorm(mu, 1) }",
    data       = list(y = 0),
    prior_list = list(mu = mu_prior, bias = bias_selection_kernel),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 100,
    seed       = 16,
    silent     = TRUE
  ))
  result <- save_fit(fit_selection_kernel_summary, "fit_selection_kernel_summary",
                     simple_priors = TRUE, pub_bias_priors = TRUE,
                     mixture_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Summary table fixture with ordinary mixture and selection-kernel bias mixture")
  model_registry[["fit_selection_kernel_summary"]] <<- result$registry_entry
  fit_selection_kernel_summary <- result$fit
})


# ============================================================================ #
# SECTION 5: SPIKE-AND-SLAB PRIORS
# ============================================================================ #
test_that("Spike-and-slab prior models fit correctly", {

  skip_if_not_installed("rjags")

  # Simple spike-and-slab
  priors_spike_slab_simple <- list(
    "mu" = prior_spike_and_slab(prior("normal", list(0, 1)),
                                prior_inclusion = prior("beta", list(1,1)))
  )

  model_syntax_ss1 <- "model{}"

  fit_spike_slab_simple <- suppressWarnings(JAGS_fit(model_syntax_ss1, data = NULL, prior_list = priors_spike_slab_simple,
                                                     chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_spike_slab_simple, "fit_spike_slab_simple",
                     spike_and_slab_priors = TRUE,
                     note = "Simple spike-and-slab with normal alternative and beta inclusion prior")
  model_registry[["fit_spike_slab_simple"]] <<- result$registry_entry
  fit_spike_slab_simple <- result$fit

  # Spike-and-slab with factor prior
  priors_spike_slab_factor <- list(
    "beta" = prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
                                  prior_inclusion = prior("beta", list(1,1)))
  )

  # Set levels attribute on the factor prior component within the spike_and_slab mixture
  # The spike_and_slab prior contains multiple components; we need to set levels on the factor component
  components <- attr(priors_spike_slab_factor$beta, "components")
  alternative_idx <- which(components == "alternative")
  # Set to 3 levels for a 3-level factor (A, B, C)
  attr(priors_spike_slab_factor$beta[[alternative_idx]], "levels") <- 3

  model_syntax_ss2 <- "model{}"

  fit_spike_slab_factor <- suppressWarnings(JAGS_fit(model_syntax_ss2, data = NULL, prior_list = priors_spike_slab_factor,
                                                     chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2))
  result <- save_fit(fit_spike_slab_factor, "fit_spike_slab_factor",
                     spike_and_slab_priors = TRUE, factor_priors = TRUE,
                     note = "Spike-and-slab with orthonormal factor prior (3 levels) as alternative")
  model_registry[["fit_spike_slab_factor"]] <<- result$registry_entry
  fit_spike_slab_factor <- result$fit
})


# ============================================================================ #
# SECTION 6: MIXTURE PRIORS
# ============================================================================ #
test_that("Mixture prior models fit correctly", {

  skip_if_not_installed("rjags")

  # Simple mixture
  priors_mixture_simple <- list(
    "mu" = prior_mixture(
      list(
        prior("normal", list(0,  1), prior_weights = 1),
        prior("normal", list(-3, 1), prior_weights = 5),
        prior("gamma",  list(5, 10), prior_weights = 1)
      ),
      is_null = c(T, F, T)
    )
  )

  model_syntax_mix1 <- "model{}"

  fit_mixture_simple <- suppressWarnings(JAGS_fit(model_syntax_mix1, data = NULL, prior_list = priors_mixture_simple,
                                                  chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_mixture_simple, "fit_mixture_simple",
                     mixture_priors = TRUE,
                     note = "Mixture of 3 components (2 normals, 1 gamma) with is_null flags")
  model_registry[["fit_mixture_simple"]] <<- result$registry_entry
  fit_mixture_simple <- result$fit

  # Mixture with components
  priors_mixture_components <- list(
    "beta" = prior_mixture(
      list(
        prior("normal", list(0,  1), prior_weights = 1),
        prior("normal", list(-3, 1), prior_weights = 5)
      ),
      components = c("b", "a")
    )
  )

  model_syntax_mix2 <- "model{}"

  fit_mixture_components <- suppressWarnings(JAGS_fit(model_syntax_mix2, data = NULL, prior_list = priors_mixture_components,
                                                      chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2))
  result <- save_fit(fit_mixture_components, "fit_mixture_components",
                     mixture_priors = TRUE,
                     note = "Mixture with named components (a, b)")
  model_registry[["fit_mixture_components"]] <<- result$registry_entry
  fit_mixture_components <- result$fit

  # Mixture with spike
  priors_mixture_spike <- list(
    "gamma" = prior_mixture(
      list(
        prior("spike", list(2)),
        prior("normal", list(-3, 1))
      )
    )
  )

  model_syntax_mix3 <- "model{}"

  fit_mixture_spike <- suppressWarnings(JAGS_fit(model_syntax_mix3, data = NULL, prior_list = priors_mixture_spike,
                                                 chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3))
  result <- save_fit(fit_mixture_spike, "fit_mixture_spike",
                     mixture_priors = TRUE,
                     note = "Mixture containing spike prior at value 2")
  model_registry[["fit_mixture_spike"]] <<- result$registry_entry
  fit_mixture_spike <- result$fit
})


# ============================================================================ #
# SECTION 7: FORMULA-BASED MODELS (SIMPLE REGRESSION)
# ============================================================================ #
test_that("Simple formula-based regression models fit correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data_formula <- data.frame(
    x_cont1 = rnorm(100),
    x_fac2t = factor(rep(c("A", "B"), 50), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(100, .4 * data_formula$x_cont1, 1),
    N = 100
  )

  # Simple linear regression
  formula_list_simple <- list(mu = ~ x_cont1)
  formula_data_list_simple <- list(mu = data_formula)
  formula_prior_list_simple <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1))
    )
  )
  prior_list_simple <- list(sigma = prior("lognormal", list(0, 1)))

  model_syntax_simple <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit_formula_simple <- JAGS_fit(
    model_syntax = model_syntax_simple, data = data, prior_list = prior_list_simple,
    formula_list = formula_list_simple, formula_data_list = formula_data_list_simple,
    formula_prior_list = formula_prior_list_simple,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)

  # Compute marginal likelihood for model averaging
  log_posterior_formula <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  marglik_formula_simple <- JAGS_bridgesampling(
    fit_formula_simple, log_posterior = log_posterior_formula, data = data,
    prior_list = prior_list_simple,
    formula_list = formula_list_simple, formula_data_list = formula_data_list_simple,
    formula_prior_list = formula_prior_list_simple)

  result <- save_fit(fit_formula_simple, "fit_formula_simple",
                     marglik = marglik_formula_simple,
                     formulas = TRUE, simple_priors = TRUE,
                     note = "Simple linear regression with continuous predictor")
  model_registry[["fit_formula_simple"]] <<- result$registry_entry
  fit_formula_simple <- result$fit

  # Regression with treatment factor
  formula_list_treatment <- list(mu = ~ x_cont1 + x_fac2t)
  formula_data_list_treatment <- list(mu = data_formula)
  formula_prior_list_treatment <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1)),
      "x_fac2t"   = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )

  fit_formula_treatment <- JAGS_fit(
    model_syntax = model_syntax_simple, data = data, prior_list = prior_list_simple,
    formula_list = formula_list_treatment, formula_data_list = formula_data_list_treatment,
    formula_prior_list = formula_prior_list_treatment,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)

  # Compute marginal likelihood for model averaging
  marglik_formula_treatment <- JAGS_bridgesampling(
    fit_formula_treatment, log_posterior = log_posterior_formula, data = data,
    prior_list = prior_list_simple,
    formula_list = formula_list_treatment, formula_data_list = formula_data_list_treatment,
    formula_prior_list = formula_prior_list_treatment)

  result <- save_fit(fit_formula_treatment, "fit_formula_treatment",
                     marglik = marglik_formula_treatment,
                     formulas = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Regression with continuous predictor and 2-level treatment factor")
  model_registry[["fit_formula_treatment"]] <<- result$registry_entry
  fit_formula_treatment <- result$fit

  # Regression with positive treatment factor
  formula_prior_list_treatment_positive <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1)),
      "x_fac2t"   = prior_factor(
        "normal",
        parameters = list(0, 1),
        truncation = list(0, Inf),
        contrast   = "treatment"
      )
    )
  )

  fit_formula_treatment_positive <- JAGS_fit(
    model_syntax = model_syntax_simple, data = data, prior_list = prior_list_simple,
    formula_list = formula_list_treatment, formula_data_list = formula_data_list_treatment,
    formula_prior_list = formula_prior_list_treatment_positive,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 4)

  marglik_formula_treatment_positive <- JAGS_bridgesampling(
    fit_formula_treatment_positive, log_posterior = log_posterior_formula, data = data,
    prior_list = prior_list_simple,
    formula_list = formula_list_treatment, formula_data_list = formula_data_list_treatment,
    formula_prior_list = formula_prior_list_treatment_positive)

  result <- save_fit(fit_formula_treatment_positive, "fit_formula_treatment_positive",
                     marglik = marglik_formula_treatment_positive,
                     formulas = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Regression with continuous predictor and positive-truncated 2-level treatment factor")
  model_registry[["fit_formula_treatment_positive"]] <<- result$registry_entry
  fit_formula_treatment_positive <- result$fit

  # Regression with negative treatment factor
  formula_prior_list_treatment_negative <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1)),
      "x_fac2t"   = prior_factor(
        "normal",
        parameters = list(0, 1),
        truncation = list(-Inf, 0),
        contrast   = "treatment"
      )
    )
  )

  fit_formula_treatment_negative <- JAGS_fit(
    model_syntax = model_syntax_simple, data = data, prior_list = prior_list_simple,
    formula_list = formula_list_treatment, formula_data_list = formula_data_list_treatment,
    formula_prior_list = formula_prior_list_treatment_negative,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 5)

  marglik_formula_treatment_negative <- JAGS_bridgesampling(
    fit_formula_treatment_negative, log_posterior = log_posterior_formula, data = data,
    prior_list = prior_list_simple,
    formula_list = formula_list_treatment, formula_data_list = formula_data_list_treatment,
    formula_prior_list = formula_prior_list_treatment_negative)

  result <- save_fit(fit_formula_treatment_negative, "fit_formula_treatment_negative",
                     marglik = marglik_formula_treatment_negative,
                     formulas = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Regression with continuous predictor and negative-truncated 2-level treatment factor")
  model_registry[["fit_formula_treatment_negative"]] <<- result$registry_entry
  fit_formula_treatment_negative <- result$fit

  # Regression with orthonormal factor
  formula_list_orthonormal <- list(mu = ~ x_cont1 + x_fac3o)
  formula_data_list_orthonormal <- list(mu = data_formula)
  formula_prior_list_orthonormal <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1)),
      "x_fac3o"   = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  fit_formula_orthonormal <- JAGS_fit(
    model_syntax = model_syntax_simple, data = data, prior_list = prior_list_simple,
    formula_list = formula_list_orthonormal, formula_data_list = formula_data_list_orthonormal,
    formula_prior_list = formula_prior_list_orthonormal,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3)

  # Compute marginal likelihood for model averaging
  marglik_formula_orthonormal <- JAGS_bridgesampling(
    fit_formula_orthonormal, log_posterior = log_posterior_formula, data = data,
    prior_list = prior_list_simple,
    formula_list = formula_list_orthonormal, formula_data_list = formula_data_list_orthonormal,
    formula_prior_list = formula_prior_list_orthonormal)

  result <- save_fit(fit_formula_orthonormal, "fit_formula_orthonormal",
                     marglik = marglik_formula_orthonormal,
                     formulas = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Regression with continuous predictor and 3-level orthonormal factor")
  model_registry[["fit_formula_orthonormal"]] <<- result$registry_entry
  fit_formula_orthonormal <- result$fit
})


# ============================================================================ #
# SECTION 8: FORMULA-BASED MODELS (INTERACTIONS)
# ============================================================================ #
test_that("Formula-based interaction models fit correctly", {

  skip_if_not_installed("rjags")

  set.seed(1)
  data_formula <- data.frame(
    x_cont1 = rnorm(100),
    x_cont2 = rnorm(100),
    x_fac2t = factor(rep(c("A", "B"), 50), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(100, .4 * data_formula$x_cont1 - 0.15 * data_formula$x_cont1 * data_formula$x_cont2, 1),
    N = 100
  )

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )
  prior_list <- list(sigma = prior("lognormal", list(0, 1)))

  # Continuous interaction
  formula_list_cont_int <- list(mu = ~ x_cont1 * x_cont2)
  formula_data_list_cont_int <- list(mu = data_formula)
  formula_prior_list_cont_int <- list(
    mu = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_cont2"         = prior("normal", list(0, 1)),
      "x_cont1:x_cont2" = prior("normal", list(0, 1))
    )
  )

  fit_formula_interaction_cont <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_cont_int, formula_data_list = formula_data_list_cont_int,
    formula_prior_list = formula_prior_list_cont_int,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)
  result <- save_fit(fit_formula_interaction_cont, "fit_formula_interaction_cont",
                     formulas = TRUE, interactions = TRUE, simple_priors = TRUE,
                     note = "Continuous-continuous interaction")
  model_registry[["fit_formula_interaction_cont"]] <<- result$registry_entry
  fit_formula_interaction_cont <- result$fit

  # Test standardization: manual vs automatic scaling
  # Create data with large scale differences (far from being scaled)
  set.seed(2)
  data_unscaled <- data.frame(
    x_cont1 = rnorm(100, mean = 1000, sd = 1000),  # Large scale
    x_cont2 = rnorm(100, mean = 0.5, sd = 0.01)    # Small scale
  )
  data_scale <- list(
    y = rnorm(100, 500 * data_unscaled$x_cont1 - 20 * data_unscaled$x_cont1 * data_unscaled$x_cont2, 1),
    N = 100
  )

  # Manual scaling: scale the data manually before fitting
  data_manual_scaled <- data_unscaled
  x_cont1_mean <- mean(data_unscaled$x_cont1)
  x_cont1_sd   <- sd(data_unscaled$x_cont1)
  x_cont2_mean <- mean(data_unscaled$x_cont2)
  x_cont2_sd   <- sd(data_unscaled$x_cont2)
  data_manual_scaled$x_cont1 <- (data_unscaled$x_cont1 - x_cont1_mean) / x_cont1_sd
  data_manual_scaled$x_cont2 <- (data_unscaled$x_cont2 - x_cont2_mean) / x_cont2_sd

  formula_list_scale <- list(mu = ~ x_cont1 * x_cont2)
  formula_prior_list_scale <- list(
    mu = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_cont2"         = prior("normal", list(0, 1)),
      "x_cont1:x_cont2" = prior("normal", list(0, 1))
    )
  )

  # Fit 1: Manual scaling
  formula_data_list_manual <- list(mu = data_manual_scaled)
  fit_formula_manual_scaled <- JAGS_fit(
    model_syntax = model_syntax, data = data_scale, prior_list = prior_list,
    formula_list = formula_list_scale, formula_data_list = formula_data_list_manual,
    formula_prior_list = formula_prior_list_scale,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  # Store scaling info as attribute for comparison
  attr(fit_formula_manual_scaled, "manual_scale") <- list(
    mu_x_cont1 = list(mean = x_cont1_mean, sd = x_cont1_sd),
    mu_x_cont2 = list(mean = x_cont2_mean, sd = x_cont2_sd)
  )
  result <- save_fit(fit_formula_manual_scaled, "fit_formula_manual_scaled",
                     formulas = TRUE, interactions = TRUE, simple_priors = TRUE,
                     note = "Manual scaling of continuous predictors")
  model_registry[["fit_formula_manual_scaled"]] <<- result$registry_entry
  fit_formula_manual_scaled <- result$fit

  # Fit 2: Automatic scaling
  formula_data_list_auto <- list(mu = data_unscaled)
  formula_scale_list_auto <- list(mu = list(x_cont1 = TRUE, x_cont2 = TRUE))
  fit_formula_auto_scaled <- JAGS_fit(
    model_syntax = model_syntax, data = data_scale, prior_list = prior_list,
    formula_list = formula_list_scale, formula_data_list = formula_data_list_auto,
    formula_prior_list = formula_prior_list_scale,
    formula_scale_list = formula_scale_list_auto,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  result <- save_fit(fit_formula_auto_scaled, "fit_formula_auto_scaled",
                     formulas = TRUE, interactions = TRUE, simple_priors = TRUE,
                     note = "Automatic scaling of continuous predictors")
  model_registry[["fit_formula_auto_scaled"]] <<- result$registry_entry
  fit_formula_auto_scaled <- result$fit

  # Continuous-factor interaction
  formula_list_mix_int <- list(mu = ~ x_cont1 * x_fac3o)
  formula_data_list_mix_int <- list(mu = data_formula)
  formula_prior_list_mix_int <- list(
    mu = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
      "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  fit_formula_interaction_mix <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_mix_int, formula_data_list = formula_data_list_mix_int,
    formula_prior_list = formula_prior_list_mix_int,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  result <- save_fit(fit_formula_interaction_mix, "fit_formula_interaction_mix",
                     formulas = TRUE, interactions = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Continuous-factor interaction with 3-level orthonormal factor")
  model_registry[["fit_formula_interaction_mix"]] <<- result$registry_entry
  fit_formula_interaction_mix <- result$fit

  # Continuous-factor interaction (Main effects only)
  formula_list_mix_main <- list(mu = ~ x_cont1 + x_fac3o)
  formula_data_list_mix_main <- list(mu = data_formula)
  formula_prior_list_mix_main <- list(
    mu = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  fit_formula_interaction_mix_main <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_mix_main, formula_data_list = formula_data_list_mix_main,
    formula_prior_list = formula_prior_list_mix_main,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  result <- save_fit(fit_formula_interaction_mix_main, "fit_formula_interaction_mix_main",
                     formulas = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Continuous-factor main effects only (for interaction test)")
  model_registry[["fit_formula_interaction_mix_main"]] <<- result$registry_entry
  fit_formula_interaction_mix_main <- result$fit

  # Factor-factor interaction
  formula_list_fac_int <- list(mu = ~ x_fac2t * x_fac3o)
  formula_data_list_fac_int <- list(mu = data_formula)
  formula_prior_list_fac_int <- list(
    mu = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac2t"         = prior_factor("normal", contrast = "treatment", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
      "x_fac2t:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  fit_formula_interaction_fac <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_fac_int, formula_data_list = formula_data_list_fac_int,
    formula_prior_list = formula_prior_list_fac_int,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3)
  result <- save_fit(fit_formula_interaction_fac, "fit_formula_interaction_fac",
                     formulas = TRUE, interactions = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Factor-factor interaction: 2-level treatment x 3-level orthonormal")
  model_registry[["fit_formula_interaction_fac"]] <<- result$registry_entry
  fit_formula_interaction_fac <- result$fit

  # Regression with prior_mixture for factor predictor
  # Testing mixture of spike and normal factor priors
  set.seed(1)
  data_formula_mix <- data.frame(
    x_cont  = rnorm(100),
    x_fac3t = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C"))
  )
  data_mix <- list(
    y = rnorm(100, 0.20 * data_formula_mix$x_cont, 1),
    N = 100
  )

  formula_list_factor_mix <- list(mu = ~ x_cont + x_fac3t)
  formula_data_list_factor_mix <- list(mu = data_formula_mix)
  formula_prior_list_factor_mix <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont"    = prior("normal", list(0, 1)),
      "x_fac3t"   = prior_mixture(list(
        prior("spike", list(0)),
        prior_factor("normal", list(0, 0.3), contrast = "treatment")
      ), is_null = c(TRUE, FALSE))
    )
  )

  fit_formula_factor_mixture <- JAGS_fit(
    model_syntax = model_syntax, data = data_mix, prior_list = prior_list,
    formula_list = formula_list_factor_mix, formula_data_list = formula_data_list_factor_mix,
    formula_prior_list = formula_prior_list_factor_mix,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 4)
  result <- save_fit(fit_formula_factor_mixture, "fit_formula_factor_mixture",
                     formulas = TRUE, mixture_priors = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Regression with mixture prior on 3-level treatment factor (spike vs normal)")
  model_registry[["fit_formula_factor_mixture"]] <<- result$registry_entry
  fit_formula_factor_mixture <- result$fit
})


# ============================================================================ #
# SECTION 9: FORMULA-BASED MODELS (MULTIPLE FORMULAS)
# ============================================================================ #
test_that("Multi-formula models fit correctly", {

  skip_if_not_installed("rjags")

  set.seed(1)
  data_formula <- data.frame(
    x_cont1 = rnorm(100),
    x_fac2t = factor(rep(c("A", "B"), 50), levels = c("A", "B"))
  )
  data_mu <- 0.20 * data_formula$x_cont1
  data_sigma <- 0.50 * exp(ifelse(data_formula$x_fac2t == "A", -0.5, 0.5))
  data <- list(
    y = rnorm(100, data_mu, data_sigma),
    N = 100
  )

  # Model with two formulas (mu and sigma)
  formula_list_multi <- list(
    mu        = ~ x_cont1,
    sigma_exp = ~ x_fac2t
  )
  formula_data_list_multi <- list(
    mu        = data_formula,
    sigma_exp = data_formula
  )
  formula_prior_list_multi <- list(
    mu         = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1))
    ),
    sigma_exp  = list(
      "intercept" = prior("spike", list(0)),
      "x_fac2t"   = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    )
  )
  prior_list_multi <- list(
    "sigma" = prior("normal", list(0, 5), list(0, Inf))
  )

  model_syntax_multi <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma * exp(sigma_exp[i]), 2))\n",
    "}\n",
    "}"
  )

  fit_formula_multi <- JAGS_fit(
    model_syntax = model_syntax_multi, data = data, prior_list = prior_list_multi,
    formula_list = formula_list_multi, formula_data_list = formula_data_list_multi,
    formula_prior_list = formula_prior_list_multi,
    chains = 2, adapt = 500, burnin = 500, sample = 500, seed = 1)
  result <- save_fit(fit_formula_multi, "fit_formula_multi",
                     formulas = TRUE, multi_formula = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Two formulas: mu (continuous) and sigma_exp (meandif factor)")
  model_registry[["fit_formula_multi"]] <<- result$registry_entry
  fit_formula_multi <- result$fit
})


# ============================================================================ #
# SECTION 10: RANDOM EFFECTS MODELS
# ============================================================================ #
test_that("Random effects models fit correctly", {

  skip_if_not_installed("rjags")

  set.seed(1)
  data_formula <- data.frame(
    x_cont1 = rnorm(100),
    x_fac3  = as.factor(sample(LETTERS[1:3], 100, replace = TRUE)),
    id      = factor(rep(LETTERS[1:10], 10))
  )
  id_values <- rnorm(10, 0, 0.5)
  names(id_values) <- LETTERS[1:10]

  data <- list(
    y = rnorm(100, 0.4 * data_formula$x_cont1 + id_values[data_formula$id]),
    N = 100
  )

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )
  prior_list <- list(sigma = prior("lognormal", list(0, 1)))

  # Random intercept only
  # Note: Using || for uncorrelated random effects (as opposed to | for correlated)
  formula_list_re_int <- list(mu = ~ 1 + (1 ||id))
  formula_data_list_re_int <- list(mu = data_formula)
  formula_prior_list_re_int <- list(
    mu = list(
      "intercept"    = prior("normal", list(0, 5))
    )
  )
  formula_random_prior_list_re_int <- list(
    mu = prior_random(id = random_block(sd = prior("normal", list(0, 1), list(0, 1))))
  )

  fit_random_intercept <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_int, formula_data_list = formula_data_list_re_int,
    formula_prior_list = formula_prior_list_re_int,
    formula_random_prior_list = formula_random_prior_list_re_int,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)
  result <- save_fit(fit_random_intercept, "fit_random_intercept",
                     formulas = TRUE, random_effects = TRUE, simple_priors = TRUE,
                     note = "Random intercept model (uncorrelated random effects)")
  model_registry[["fit_random_intercept"]] <<- result$registry_entry
  fit_random_intercept <- result$fit

  # Random slope (no intercept)
  formula_list_re_slope <- list(mu = ~ 1 + (0 + x_cont1 ||id))
  formula_data_list_re_slope <- list(mu = data_formula)
  formula_prior_list_re_slope <- list(
    mu = list(
      "intercept"  = prior("normal", list(0, 5))
    )
  )
  formula_random_prior_list_re_slope <- list(
    mu = prior_random(id = random_block(sd = prior("normal", list(0, 1), list(0, 1))))
  )

  fit_random_slope <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_slope, formula_data_list = formula_data_list_re_slope,
    formula_prior_list = formula_prior_list_re_slope,
    formula_random_prior_list = formula_random_prior_list_re_slope,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  result <- save_fit(fit_random_slope, "fit_random_slope",
                     formulas = TRUE, random_effects = TRUE, simple_priors = TRUE,
                     note = "Random slope for continuous predictor (no random intercept)")
  model_registry[["fit_random_slope"]] <<- result$registry_entry
  fit_random_slope <- result$fit

  # Random factor slope
  formula_list_re_fac <- list(mu = ~ 1 + x_cont1 + (x_fac3 ||id))
  formula_data_list_re_fac <- list(mu = data_formula)
  formula_prior_list_re_fac <- list(
    mu = list(
      "intercept"    = prior("normal", list(0, 5)),
      "x_cont1"      = prior("normal", list(0, 1))
    )
  )
  formula_random_prior_list_re_fac <- list(
    mu = prior_random(id = random_block(sd = prior("normal", list(0, 1), list(0, 1))))
  )

  fit_random_factor_slope <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_fac, formula_data_list = formula_data_list_re_fac,
    formula_prior_list = formula_prior_list_re_fac,
    formula_random_prior_list = formula_random_prior_list_re_fac,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3)
  result <- save_fit(fit_random_factor_slope, "fit_random_factor_slope",
                     formulas = TRUE, random_effects = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Random factor slopes with random intercept")
  model_registry[["fit_random_factor_slope"]] <<- result$registry_entry
  fit_random_factor_slope <- result$fit

  # Random factor slope with orthonormal contrast
  formula_list_re_fac <- list(mu = ~ 1 + x_fac3 + (x_fac3 ||id))
  formula_data_list_re_fac <- list(mu = data_formula)
  formula_prior_list_re_fac <- list(
    mu = list(
      "intercept"    = prior("normal", list(0, 5)),
      "x_fac3"       = prior_factor("mnormal", list(0, 1))
    )
  )
  formula_random_prior_list_re_fac <- list(
    mu = prior_random(id = random_block(sd = prior("normal", list(0, 1), list(0, 1))))
  )

  fit_random_factor_slope2 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_fac, formula_data_list = formula_data_list_re_fac,
    formula_prior_list = formula_prior_list_re_fac,
    formula_random_prior_list = formula_random_prior_list_re_fac,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3)
  result <- save_fit(fit_random_factor_slope2, "fit_random_factor_slope2",
                     formulas = TRUE, random_effects = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Random factor slopes with random intercept")
  model_registry[["fit_random_factor_slope2"]] <<- result$registry_entry
  fit_random_factor_slope2 <- result$fit


  # Random factor slope independent spike and slab contrast
  formula_list_re_fac <- list(mu = ~ -1 + x_fac3 + (x_fac3 - 1 ||id))
  formula_data_list_re_fac <- list(mu = data_formula)
  formula_prior_list_re_fac <- list(
    mu = list(
      "x_fac3"       = prior_factor("normal", list(0, 1), contrast = "independent")
    )
  )
  formula_random_prior_list_re_fac <- list(
    mu = prior_random(id = random_block(sd = prior_spike_and_slab(prior("normal", list(0, 1), list(0, 1)))))
  )

  fit_random_factor_slope3 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_fac, formula_data_list = formula_data_list_re_fac,
    formula_prior_list = formula_prior_list_re_fac,
    formula_random_prior_list = formula_random_prior_list_re_fac,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3)
  result <- save_fit(fit_random_factor_slope3, "fit_random_factor_slope3",
                     formulas = TRUE, random_effects = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Random factor slopes with random intercept")
  model_registry[["fit_random_factor_slope3"]] <<- result$registry_entry
  fit_random_factor_slope3 <- result$fit
})


# ============================================================================ #
# SECTION 11: SPIKE FACTOR PRIORS
# ============================================================================ #
test_that("Spike factor prior models fit correctly", {

  skip_if_not_installed("rjags")

  set.seed(1)
  data_formula <- data.frame(
    x_fac2i = factor(rep(c("A", "B"), 50), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C")),
    x_fac3t = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C")),
    x_fac3md = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C"))
  )
  data <- list(y = rnorm(100, 0, 1), N = 100)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )
  prior_list <- list(sigma = prior("lognormal", list(0, 1)))

  # Spike priors with different contrasts
  # Note: Using - 1 to remove the intercept since spike priors for independent factors
  # define all levels explicitly, and we're testing different contrast behaviors
  formula_list_spike <- list(mu = ~ x_fac2i + x_fac3o + x_fac3t + x_fac3md - 1)
  formula_data_list_spike <- list(mu = data_formula)
  formula_prior_list_spike <- list(
    mu = list(
      "x_fac2i"  = prior_factor("spike", contrast = "independent", list(1)),
      "x_fac3o"  = prior_factor("spike", contrast = "orthonormal", list(0)),
      "x_fac3t"  = prior_factor("spike", contrast = "treatment", list(2)),
      "x_fac3md" = prior_factor("spike", contrast = "meandif", list(0))
    )
  )

  fit_spike_factors <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_spike, formula_data_list = formula_data_list_spike,
    formula_prior_list = formula_prior_list_spike,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)
  result <- save_fit(fit_spike_factors, "fit_spike_factors",
                     formulas = TRUE, factor_priors = TRUE,
                     note = "Spike priors with all 4 contrast types: independent, orthonormal, treatment, meandif")
  model_registry[["fit_spike_factors"]] <<- result$registry_entry
  fit_spike_factors <- result$fit
})


# ============================================================================ #
# SECTION 12: JOINT MODELS (FORMULA + SPIKE-AND-SLAB + MIXTURE)
# ============================================================================ #
test_that("Joint complex models fit correctly", {

  skip_if_not_installed("rjags")

  set.seed(1)
  data_formula <- data.frame(
    x_cont1 = rnorm(100),
    x_fac2t = factor(rep(c("A", "B"), 50), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), length.out = 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(100, 0.20 * data_formula$x_cont1, 1),
    N = 100
  )

  # Model with mixture intercept, spike-and-slab continuous, spike-and-slab factor
  formula_list_joint <- list(mu = ~ x_cont1 + x_fac3t)
  formula_data_list_joint <- list(mu = data_formula)
  formula_prior_list_joint <- list(
    mu = list(
      "intercept" = prior_mixture(
        list(
          prior("spike",  list(0),       prior_weights = 2),
          prior("normal", list(-1, 0.5), prior_weights = 1),
          prior("normal", list( 1, 0.5), prior_weights = 1)
        ),
        is_null = c(T, F, F)
      ),
      "x_cont1" = prior_mixture(
        list(
          prior("spike",  list(0),    prior_weights = 1),
          prior("normal", list(0, 1), prior_weights = 1)
        ),
        is_null = c(T, F)
      ),
      "x_fac3t" = prior_spike_and_slab(
        prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
        prior_inclusion = prior("spike", list(0.5))
      )
    )
  )
  # Scale the continuous predictor by sigma (standard practice for hierarchical centering)
  attr(formula_prior_list_joint$mu$x_cont1, "multiply_by") <- "sigma"

  prior_list_joint <- list(
    "sigma" = prior_mixture(
      list(
        prior("normal",    list(0, 1), truncation = list(0, Inf)),
        prior("lognormal", list(0, 1))
      ),
      is_null = c(T, F)
    )
  )

  model_syntax_joint <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit_joint_complex <- JAGS_fit(
    model_syntax = model_syntax_joint, data = data, prior_list = prior_list_joint,
    formula_list = formula_list_joint, formula_data_list = formula_data_list_joint,
    formula_prior_list = formula_prior_list_joint,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)
  result <- save_fit(fit_joint_complex, "fit_joint_complex",
                     formulas = TRUE, mixture_priors = TRUE, spike_and_slab_priors = TRUE,
                     factor_priors = TRUE, simple_priors = TRUE,
                     note = "Complex model: mixture intercept, mixture sigma, spike-and-slab continuous, spike-and-slab factor")
  model_registry[["fit_joint_complex"]] <<- result$registry_entry
  fit_joint_complex <- result$fit
})


# ============================================================================ #
# SECTION 13: EXPRESSION PRIORS
# ============================================================================ #
test_that("Expression prior models fit correctly", {

  skip_if_not_installed("rjags")

  # Simple prior with expression
  priors_expr_simple <- list(
    x        = prior("normal",   list(0, expression(x_sigma))),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax_expr1 <- "model{}"

  fit_expression_simple <- suppressWarnings(JAGS_fit(model_syntax_expr1, data = NULL, prior_list = priors_expr_simple,
                                                     chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1))
  result <- save_fit(fit_expression_simple, "fit_expression_simple",
                     expression_priors = TRUE, simple_priors = TRUE,
                     note = "Normal prior with expression referencing another parameter (x_sigma)")
  model_registry[["fit_expression_simple"]] <<- result$registry_entry
  fit_expression_simple <- result$fit

  # Spike-and-slab with expression
  priors_expr_ss <- list(
    x        = prior_spike_and_slab(
      prior("normal", list(0, expression(x_sigma)))
    ),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax_expr2 <- "model{}"

  fit_expression_spike_slab <- suppressWarnings(JAGS_fit(model_syntax_expr2, data = NULL, prior_list = priors_expr_ss,
                                                         chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2))
  result <- save_fit(fit_expression_spike_slab, "fit_expression_spike_slab",
                     expression_priors = TRUE, spike_and_slab_priors = TRUE, simple_priors = TRUE,
                     note = "Spike-and-slab with expression in alternative prior")
  model_registry[["fit_expression_spike_slab"]] <<- result$registry_entry
  fit_expression_spike_slab <- result$fit

  # Mixture with expression
  priors_expr_mix <- list(
    x        = prior_mixture(list(
      prior("normal", list(0, expression(x_sigma))),
      prior("cauchy", list(0, 1))
    ), is_null = c(T, F)),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax_expr3 <- "model{}"

  fit_expression_mixture <- suppressWarnings(JAGS_fit(model_syntax_expr3, data = NULL, prior_list = priors_expr_mix,
                                                      chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3))
  result <- save_fit(fit_expression_mixture, "fit_expression_mixture",
                     expression_priors = TRUE, mixture_priors = TRUE, simple_priors = TRUE,
                     note = "Mixture prior with expression in one component")
  model_registry[["fit_expression_mixture"]] <<- result$registry_entry
  fit_expression_mixture <- result$fit
})


# ============================================================================ #
# SECTION 14: ADVANCED JAGS_FIT FEATURES
# ============================================================================ #
test_that("Advanced JAGS_fit features work correctly", {

  skip_if_not_installed("rjags")

  set.seed(1)
  data <- list(
    x = rnorm(20, 0, 1),
    N = 20
  )
  priors_list <- list(
    m = prior("normal", list(0, 1)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"

  # Test 1: add_parameters - monitoring additional parameters not in prior_list
  model_syntax_add_param <-
    "model
    {
      g ~ dnorm(0, 1)
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"

  log_posterior <- function(parameters, data){
    return(stats::dnorm(parameters[["g"]], log = TRUE))
    #return(sum(stats::dnorm(data$x, mean = parameters[["m"]], sd = parameters[["s"]], log = TRUE)))
  }
  add_l <- c("g" = -Inf)
  add_u <- c("g" =  Inf)

  fit_add_parameters <- JAGS_fit(model_syntax_add_param, data, priors_list,
                                 add_parameters = "g",
                                 chains = 2, adapt = 100, burnin = 100, sample = 300, seed = 1)
  marglik_fit_add_parameters <- JAGS_bridgesampling(
    fit                = fit_add_parameters,
    log_posterior      = log_posterior,
    data               = data,
    prior_list         = priors_list,
    add_parameters     = "g",
    add_bounds         = list("lb" = add_l, "ub" = add_u)
    )

  result <- save_fit(fit_add_parameters, "fit_add_parameters",
                     simple_priors = TRUE, add_parameters = TRUE,
                     note = "Model with additional monitored parameter 'g' not in prior_list")
  model_registry[["fit_add_parameters"]] <<- result$registry_entry
  fit_add_parameters <- result$fit

  # Verify that 'g' is in the output
  expect_true("g" %in% colnames(fit_add_parameters$mcmc[[1]]))
  expect_equal(ncol(fit_add_parameters$mcmc[[1]]), 3) # m, s, g

  # Test 2: autofit - automatic refitting until convergence
  # Using a model that requires more samples to converge
  priors_autofit <- list(
    m = prior("normal", list(0, 1))
  )
  data_autofit <- list(
    x = c(-500),
    N = 1
  )
  model_syntax_autofit <-
    "model
    {
      l = 1
      for(i in 1:N){
        x[i] ~ dt(m, pow(.3, -2), 1)
      }
    }"

  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  # First fit without autofit (should have poor convergence)
  fit_no_autofit <- JAGS_fit(model_syntax_autofit, data_autofit, priors_autofit,
                             autofit = FALSE,
                             chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 2)
  result <- save_fit(fit_no_autofit, "fit_no_autofit",
                     simple_priors = TRUE,
                     note = "Model without autofit (poor convergence expected)")
  model_registry[["fit_no_autofit"]] <<- result$registry_entry
  fit_no_autofit <- result$fit

  summary_no_autofit <- suppressWarnings(summary(fit_no_autofit))
  # Check that convergence is poor
  expect_true(summary_no_autofit[1,"MCerr"] > 0.069 || summary_no_autofit[1,"MC%ofSD"] > 8)

  # Now fit with autofit using max_error criterion
  fit_autofit_error <- JAGS_fit(model_syntax_autofit, data_autofit, priors_autofit,
                                autofit = TRUE,
                                autofit_control = list(max_error = 0.05, sample_extend = 100),
                                chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 2)
  result <- save_fit(fit_autofit_error, "fit_autofit_error",
                     simple_priors = TRUE, autofit = TRUE,
                     note = "Autofit with max_error criterion (< 0.05)")
  model_registry[["fit_autofit_error"]] <<- result$registry_entry
  fit_autofit_error <- result$fit

  summary_autofit_error <- summary(fit_autofit_error)
  # Should have better convergence
  expect_true(summary_autofit_error[1,"MCerr"] < 0.05)

  # Test autofit with min_ESS criterion
  fit_autofit_ess <- JAGS_fit(model_syntax_autofit, data_autofit, priors_autofit,
                              autofit = TRUE,
                              autofit_control = list(min_ESS = 200, sample_extend = 100),
                              chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 3)
  result <- save_fit(fit_autofit_ess, "fit_autofit_ess",
                     simple_priors = TRUE, autofit = TRUE,
                     note = "Autofit with min_ESS criterion (> 200)")
  model_registry[["fit_autofit_ess"]] <<- result$registry_entry
  fit_autofit_ess <- result$fit

  summary_autofit_ess <- summary(fit_autofit_ess)
  expect_true(summary_autofit_ess[1,"SSeff"] > 200)

  # Test 3: parallel - running chains in parallel
  # Note: parallel execution is tested but results should be the same as non-parallel
  fit_parallel <- JAGS_fit(model_syntax, data, priors_list,
                           parallel = TRUE, cores = 2,
                           chains = 2, adapt = 100, burnin = 100, sample = 300, seed = 4)
  result <- save_fit(fit_parallel, "fit_parallel",
                     simple_priors = TRUE, parallel = TRUE,
                     note = "Model fitted with parallel chains (cores=2)")
  model_registry[["fit_parallel"]] <<- result$registry_entry
  fit_parallel <- result$fit

  # Verify the fit worked and has the expected structure
  expect_equal(length(fit_parallel$mcmc), 2) # 2 chains
  expect_true(all(sapply(fit_parallel$mcmc, function(mcmc) ncol(mcmc) == 2))) # m and s
})


# ============================================================================ #
# SECTION 15: MODELS FOR MARGINAL DISTRIBUTION TESTING
# ============================================================================ #
# These models test marginal_posterior, ensemble_inference, and mix_posteriors
# with complex formulas including interactions and multiply_by scaling.
test_that("Marginal distribution models fit correctly", {

  skip_if_not_installed("rjags")

  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data_formula_marg <- data.frame(
    x_cont1  = rnorm(180),
    x_fac2t  = factor(rep(c("A", "B"), 90), levels = c("A", "B")),
    x_fac3md = factor(rep(c("A", "B", "C"), 60), levels = c("A", "B", "C"))
  )
  data_marg <- list(
    y = rnorm(180, 0.1, 0.5) + 0.5 + 0.20 * data_formula_marg$x_cont1 +
      ifelse(data_formula_marg$x_fac3md == "A", 0.15, ifelse(data_formula_marg$x_fac3md == "B", -0.15, 0)),
    N = 180
  )

  # Null model: spike priors on factor effects
  prior_list_marg_0 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("spike", contrast = "treatment", list(0)),
    "x_fac3md"         = prior_factor("spike", contrast = "meandif",   list(0)),
    "x_cont1:x_fac3md" = prior_factor("spike", contrast = "meandif",   list(0))
  )
  attr(prior_list_marg_0$x_cont1, "multiply_by") <- "sigma"

  # Alternative model: normal priors on factor effects
  prior_list_marg_1 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("normal",  contrast = "treatment", list(0, 1.00)),
    "x_fac3md"         = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25)),
    "x_cont1:x_fac3md" = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25))
  )
  attr(prior_list_marg_1$x_cont1, "multiply_by") <- "sigma"

  prior_list_marg <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 5))
  )
  model_syntax_marg <- paste0(
    "model{",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )
  log_posterior_marg <- function(parameters, data){
    return(sum(stats::dnorm(data$y, mean = parameters[["mu"]], sd = parameters[["sigma"]], log = TRUE)))
  }
  model_formula_marg <- list(mu = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md)

  # Fit null model
  fit_marginal_0 <- JAGS_fit(
    model_syntax = model_syntax_marg, data = data_marg,
    prior_list = prior_list_marg,
    formula_list       = model_formula_marg,
    formula_prior_list = list(mu = prior_list_marg_0),
    formula_data_list  = list(mu = data_formula_marg),
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)

  marglik_marginal_0 <- JAGS_bridgesampling(
    fit                = fit_marginal_0,
    log_posterior      = log_posterior_marg,
    data               = data_marg,
    prior_list         = prior_list_marg,
    formula_list       = model_formula_marg,
    formula_prior_list = list(mu = prior_list_marg_0),
    formula_data_list  = list(mu = data_formula_marg))

  result <- save_fit(fit_marginal_0, "fit_marginal_0",
                     marglik = marglik_marginal_0,
                     formulas = TRUE, factor_priors = TRUE, interactions = TRUE,
                     note = "Marginal dist null model: spike priors on factors with interaction and multiply_by")
  model_registry[["fit_marginal_0"]] <<- result$registry_entry
  fit_marginal_0 <- result$fit

  # Fit alternative model
  fit_marginal_1 <- JAGS_fit(
    model_syntax = model_syntax_marg, data = data_marg,
    prior_list = prior_list_marg,
    formula_list       = model_formula_marg,
    formula_prior_list = list(mu = prior_list_marg_1),
    formula_data_list  = list(mu = data_formula_marg),
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)

  marglik_marginal_1 <- JAGS_bridgesampling(
    fit                = fit_marginal_1,
    log_posterior      = log_posterior_marg,
    data               = data_marg,
    prior_list         = prior_list_marg,
    formula_list       = model_formula_marg,
    formula_prior_list = list(mu = prior_list_marg_1),
    formula_data_list  = list(mu = data_formula_marg))

  result <- save_fit(fit_marginal_1, "fit_marginal_1",
                     marglik = marglik_marginal_1,
                     formulas = TRUE, factor_priors = TRUE, interactions = TRUE,
                     note = "Marginal dist alt model: normal priors on factors with interaction and multiply_by")
  model_registry[["fit_marginal_1"]] <<- result$registry_entry
  fit_marginal_1 <- result$fit

  # Spike-and-slab/mixture model for marginal distributions
  prior_list_marg_ss <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior_mixture(list(
      prior("spike", list(0)),
      prior("normal", list(0, 1))
    ), is_null = c(T, F)),
    "x_fac2t"          = prior_spike_and_slab(prior_factor("normal",  contrast = "treatment", list(0, 1.00))),
    "x_fac3md"         = prior_spike_and_slab(prior_factor("mnormal", contrast = "meandif",   list(0, 0.25))),
    "x_cont1:x_fac3md" = prior_spike_and_slab(prior_factor("mnormal", contrast = "meandif",   list(0, 0.25)))
  )
  attr(prior_list_marg_ss$x_cont1, "multiply_by") <- "sigma"

  fit_marginal_ss <- JAGS_fit(
    model_syntax = model_syntax_marg, data = data_marg,
    prior_list = prior_list_marg,
    formula_list       = model_formula_marg,
    formula_prior_list = list(mu = prior_list_marg_ss),
    formula_data_list  = list(mu = data_formula_marg),
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 3)

  result <- save_fit(fit_marginal_ss, "fit_marginal_ss",
                     formulas = TRUE, factor_priors = TRUE, interactions = TRUE,
                     spike_and_slab_priors = TRUE, mixture_priors = TRUE,
                     note = "Marginal dist model: spike-and-slab and mixture priors with interaction and multiply_by")
  model_registry[["fit_marginal_ss"]] <<- result$registry_entry
  fit_marginal_ss <- result$fit
})


# ============================================================================ #
# SECTION: MODELS FOR ENSEMBLE PLOTS TESTING
# ============================================================================ #
test_that("PET-PEESE models fit correctly", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data <- NULL
  model_syntax <- "model{}"
  log_posterior <- function(parameters, data){ return(0) }

  # PET model
  priors_pet <- list(
    mu    = prior("spike", list(0)),
    PET   = prior_PET("normal", list(0, .2))
  )
  fit_pet <- suppressWarnings(JAGS_fit(model_syntax, data, priors_pet, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 0))
  marglik_pet <- JAGS_bridgesampling(fit_pet, log_posterior = log_posterior, data = data, prior_list = priors_pet)
  result <- save_fit(fit_pet, "fit_pet", marglik = marglik_pet, pub_bias_priors = TRUE, note = "PET prior only")
  model_registry[["fit_pet"]] <<- result$registry_entry

  # PEESE model
  priors_peese <- list(
    mu    = prior("spike", list(0)),
    PEESE = prior_PEESE("normal", list(0, .8))
  )
  fit_peese <- suppressWarnings(JAGS_fit(model_syntax, data, priors_peese, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  marglik_peese <- JAGS_bridgesampling(fit_peese, log_posterior = log_posterior, data = data, prior_list = priors_peese)
  result <- save_fit(fit_peese, "fit_peese", marglik = marglik_peese, pub_bias_priors = TRUE, note = "PEESE prior only")
  model_registry[["fit_peese"]] <<- result$registry_entry

  # Missing model (overwhelming)
  priors_missing <- list(
    mu = prior("normal", list(.2, .2), prior_weights = 4)
  )
  fit_missing <- suppressWarnings(JAGS_fit(model_syntax, data, priors_missing, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  marglik_missing <- JAGS_bridgesampling(fit_missing, log_posterior = log_posterior, data = data, prior_list = priors_missing)
  result <- save_fit(fit_missing, "fit_missing", marglik = marglik_missing, simple_priors = TRUE, note = "Overwhelming missing model")
  model_registry[["fit_missing"]] <<- result$registry_entry
})

test_that("Weightfunction models fit correctly", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data <- NULL
  model_syntax <- "model{}"
  log_posterior <- function(parameters, data){ return(0) }

  # One-sided
  priors_wf_onesided <- list(
    omega = prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 1)))
  )
  fit_wf_onesided <- suppressWarnings(JAGS_fit(model_syntax, data, priors_wf_onesided, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 0))
  marglik_wf_onesided <- JAGS_bridgesampling(fit_wf_onesided, log_posterior = log_posterior, data = data, prior_list = priors_wf_onesided)
  result <- save_fit(fit_wf_onesided, "fit_wf_onesided", marglik = marglik_wf_onesided, weightfunction_priors = TRUE, note = "One-sided weightfunction")
  model_registry[["fit_wf_onesided"]] <<- result$registry_entry

  # Two-sided
  priors_wf_twosided <- list(
    omega = prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)))
  )
  fit_wf_twosided <- suppressWarnings(JAGS_fit(model_syntax, data, priors_wf_twosided, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  marglik_wf_twosided <- JAGS_bridgesampling(fit_wf_twosided, log_posterior = log_posterior, data = data, prior_list = priors_wf_twosided)
  result <- save_fit(fit_wf_twosided, "fit_wf_twosided", marglik = marglik_wf_twosided, weightfunction_priors = TRUE, note = "Two-sided weightfunction")
  model_registry[["fit_wf_twosided"]] <<- result$registry_entry

  # Missing model for WF (overwhelming)
  priors_wf_missing <- list(
    mu = prior("normal", list(0, .8), prior_weights = 4)
  )
  fit_wf_missing <- suppressWarnings(JAGS_fit(model_syntax, data, priors_wf_missing, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  marglik_wf_missing <- JAGS_bridgesampling(fit_wf_missing, log_posterior = log_posterior, data = data, prior_list = priors_wf_missing)
  result <- save_fit(fit_wf_missing, "fit_wf_missing", marglik = marglik_wf_missing, simple_priors = TRUE, note = "Overwhelming missing model for WF")
  model_registry[["fit_wf_missing"]] <<- result$registry_entry
})

test_that("Orthonormal contrast models fit correctly", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)
  data_formula <- data.frame(
    x_fac3o = factor(rep(c("A", "B", "C"), 40), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(120, .4  + ifelse(data_formula$x_fac3o == "A", 0.0, ifelse(data_formula$x_fac3o == "B", -0.5, 0.5)), 1),
    N = 120
  )

  formula_list0 <- list(mu = ~ 1)
  formula_list1 <- list(mu = ~ x_fac3o)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 0.5))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit_orthonormal_0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  marglik_orthonormal_0 <- JAGS_bridgesampling(
    fit_orthonormal_0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  result <- save_fit(fit_orthonormal_0, "fit_orthonormal_0", marglik = marglik_orthonormal_0, formulas = TRUE, factor_priors = TRUE, note = "Orthonormal null model")
  model_registry[["fit_orthonormal_0"]] <<- result$registry_entry

  fit_orthonormal_1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)
  marglik_orthonormal_1 <- JAGS_bridgesampling(
    fit_orthonormal_1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)
  result <- save_fit(fit_orthonormal_1, "fit_orthonormal_1", marglik = marglik_orthonormal_1, formulas = TRUE, factor_priors = TRUE, note = "Orthonormal alternative model")
  model_registry[["fit_orthonormal_1"]] <<- result$registry_entry
})


# ============================================================================ #
# SECTION 2: COMPLEX MODELS FOR PLOTTING
# ============================================================================ #
test_that("Complex models for plotting fit correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_if_not_installed("RoBMA")
  requireNamespace("RoBMA", quietly = TRUE)

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(300),
    x_fac2t = factor(rep(c("A", "B"), 150), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, -0.15 + 0.20 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.2)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 300
  )

  # create model with mix of a formula and free parameters ---
  formula_list1 <- list(
    mu    = ~ x_cont1 + x_fac2t + x_fac3t
  )
  formula_data_list1 <- list(
    mu    = data_formula
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"  = prior_mixture(
        list(
          prior("spike",   list(0),       prior_weights = 2),
          prior("normal",  list(-1, 0.5), prior_weights = 1),
          prior("normal",  list( 1, 0.5), prior_weights = 1)
        ),
        is_null = c(TRUE, FALSE, FALSE)
      ),
      "x_cont1"    = prior_spike_and_slab(prior("normal",  list(0, 1), prior_weights = 1)),
      "x_fac2t"    = prior_mixture(list(
          prior("spike", list(0), prior_weights = 1),
          prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
        ),
        is_null = c(TRUE, FALSE)
      ),
      "x_fac3t"    = prior_mixture(list(
          prior("spike", list(0), prior_weights = 1),
          prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
        ),
        is_null = c(TRUE, FALSE)
      )
    )
  )

  attr(formula_prior_list1$mu$x_cont1, "multiply_by") <- "sigma"
  prior_list1 <- list(
    "sigma" = prior_mixture(
      list(
        prior("normal",    list(0, 1), truncation = list(0, Inf)),
        prior("lognormal", list(0, 1))
      ),
      components = c("normal", "lognormal")
    ),
    "bias"  = prior_mixture(list(
      prior_none(prior_weights = 1),
      prior_weightfunction("two-sided", c(0.05), wf_cumulative(c(1, 1)), prior_weights = 1/3),
      prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1)), prior_weights = 1/3),
      prior_PET("normal", list(0, 1), prior_weights = 1/3)
    ), is_null = c(TRUE, FALSE, FALSE, FALSE))
  )
  model_syntax1 <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit_complex_mixed <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1,
    chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)

  result <- save_fit(fit_complex_mixed, "fit_complex_mixed",
                     formulas = TRUE, mixture_priors = TRUE, spike_and_slab_priors = TRUE,
                     pub_bias_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Complex model with formula, mixtures, spike and slab, and publication bias")
  model_registry[["fit_complex_mixed"]] <<- result$registry_entry
  fit_complex_mixed <- result$fit

  expect_true(file.exists(file.path(temp_fits_dir, "fit_complex_mixed.RDS")))

  # Simple formula mixed model
  formula_list_simple_mixed <- list(
    mu    = ~ x_cont1 + x_fac2t + x_fac3t
  )
  formula_data_list_simple_mixed <- list(
    mu    = data_formula
  )
  formula_prior_list_simple_mixed <- list(
    mu    = list(
      "intercept"  = prior("normal",  list(-1, 0.5), prior_weights = 1),
      "x_cont1"    = prior("normal",  list(0, 1), prior_weights = 1),
      "x_fac2t"    = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x_fac3t"    = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    )
  )

  attr(formula_prior_list_simple_mixed$mu$x_cont1, "multiply_by") <- "sigma"
  prior_list_simple_mixed <- list(
    "sigma" =  prior("lognormal", list(0, 1))
  )
  model_syntax_simple_mixed <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit_simple_formula_mixed <- JAGS_fit(
    model_syntax = model_syntax_simple_mixed, data = data, prior_list = prior_list_simple_mixed,
    formula_list = formula_list_simple_mixed, formula_data_list = formula_data_list_simple_mixed, formula_prior_list = formula_prior_list_simple_mixed,
    chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)

  result <- save_fit(fit_simple_formula_mixed, "fit_simple_formula_mixed",
                     formulas = TRUE, factor_priors = TRUE, simple_priors = TRUE,
                     note = "Simple formula model with continuous, orthonormal factor, and meandif factor")
  model_registry[["fit_simple_formula_mixed"]] <<- result$registry_entry
  fit_simple_formula_mixed <- result$fit
})

# ============================================================================ #
# SECTION 3: COMPLEX BIAS ONLY MODEL FOR PLOTTING
# ============================================================================ #
test_that("Complex models for plotting fit correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_if_not_installed("RoBMA")
  requireNamespace("RoBMA", quietly = TRUE)

  set.seed(1)

  prior_list1 <- list(
    "mu"    = prior("gamma", list(3, 3)),
    "bias"  = prior_mixture(list(
      prior_none(prior_weights = 1),
      prior_weightfunction("two-sided", c(0.05), wf_cumulative(c(1, 1)), prior_weights = 1/3),
      prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1)), prior_weights = 1/3),
      prior_PET("normal", list(0, 1), prior_weights = 1/3),
      prior_PEESE("normal", list(0, 2), prior_weights = 1/3)
    ), is_null = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  )
  model_syntax1 <- "model{}"

  fit_complex_bias <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax1, data = NULL, prior_list = prior_list1,
    chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1))

  result <- save_fit(fit_complex_bias, "fit_complex_bias",
                     formulas = FALSE, mixture_priors = TRUE, spike_and_slab_priors = FALSE,
                     pub_bias_priors = TRUE, weightfunction_priors = TRUE,
                     note = "Model with complex publication bias mixture prior")
  model_registry[["fit_complex_bias"]] <<- result$registry_entry
  fit_complex_bias <- result$fit
})


# ============================================================================ #
# SECTION 4: DUAL PARAMETER REGRESSION WITH LOG(INTERCEPT) AND FORMULA_SCALE
# ============================================================================ #
test_that("Dual parameter regression with log(intercept) and formula_scale fits correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(1)

  # Generate data with heteroscedastic variance
  n <- 1000
  data_formula_dual <- data.frame(
    x_mu    = rnorm(n, mean = 5, sd = 2),
    x_sigma = rnorm(n, mean = 3, sd = 1.5)
  )

  # True parameters
  true_mu    <- 1 + 0.3 * data_formula_dual$x_mu
  true_sigma <- exp(log(0.5) - 0.2 * data_formula_dual$x_sigma)
  y <- rnorm(n, mean = true_mu, sd = true_sigma)

  data_dual <- list(y = y, N = n)

  # Formula for mu (standard intercept)
  formula_mu <- ~ x_mu

  # Formula for log_sigma with log(intercept) attribute
  formula_log_sigma <- ~ x_sigma
  attr(formula_log_sigma, "log(intercept)") <- TRUE

  formula_list_dual <- list(
    mu        = formula_mu,
    log_sigma = formula_log_sigma
  )

  formula_data_list_dual <- list(
    mu        = data_formula_dual,
    log_sigma = data_formula_dual
  )

  # Scale both continuous predictors
  formula_scale_list_dual <- list(
    mu        = list(x_mu = TRUE),
    log_sigma = list(x_sigma = TRUE)
  )

  formula_prior_list_dual <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 2)),
      "x_mu"      = prior("normal", list(0, 1))
    ),
    log_sigma = list(
      "intercept" = prior("lognormal", list(0, 0.5)),
      "x_sigma"   = prior("normal", list(0, 0.5))
    )
  )

  # Model syntax uses exp() on log_sigma to get positive sigma
  model_syntax_dual <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(exp(log_sigma[i]), 2))\n",
    "}\n",
    "}"
  )

  # Log posterior for marginal likelihood
  log_posterior_dual <- function(parameters, data){
    sigma <- exp(parameters[["log_sigma"]])
    sum(stats::dnorm(data$y, parameters[["mu"]], sigma, log = TRUE))
  }

  fit_dual_param_regression <- JAGS_fit(
    model_syntax       = model_syntax_dual,
    data               = data_dual,
    prior_list         = NULL,
    formula_list       = formula_list_dual,
    formula_data_list  = formula_data_list_dual,
    formula_prior_list = formula_prior_list_dual,
    formula_scale_list = formula_scale_list_dual,
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)

  marglik_dual_param_regression <- JAGS_bridgesampling(
    fit                = fit_dual_param_regression,
    log_posterior      = log_posterior_dual,
    data               = data_dual,
    prior_list         = NULL,
    formula_list       = formula_list_dual,
    formula_data_list  = formula_data_list_dual,
    formula_prior_list = formula_prior_list_dual,
    formula_scale_list = formula_scale_list_dual)

  result <- save_fit(fit_dual_param_regression, "fit_dual_param_regression",
                     marglik = marglik_dual_param_regression,
                     formulas = TRUE, simple_priors = TRUE,
                     note = "Dual parameter regression: mu and log_sigma with log(intercept) and formula_scale")
  model_registry[["fit_dual_param_regression"]] <<- result$registry_entry
  fit_dual_param_regression <- result$fit

  # Verify the model has the expected structure
  expect_true("mu_intercept" %in% colnames(fit_dual_param_regression$mcmc[[1]]))
  expect_true("mu_x_mu" %in% colnames(fit_dual_param_regression$mcmc[[1]]))
  expect_true("log_sigma_intercept" %in% colnames(fit_dual_param_regression$mcmc[[1]]))
  expect_true("log_sigma_x_sigma" %in% colnames(fit_dual_param_regression$mcmc[[1]]))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-fit-edge-cases.R
# ============================================================================ #

skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS Fit Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for JAGS fitting functions including input validation,
#   error handling, and boundary conditions.
#
# DEPENDENCIES:
#   - rjags: For JAGS model syntax generation and testing
#   - common-functions.R: REFERENCE_DIR, test_reference_text, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"): For all tests
#
# MODELS/FIXTURES:
#   - Some tests use pre-fitted models from test-00-model-fits.R
#
# TAGS: @edge-cases, @JAGS, @input-validation
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-fit-edge-cases")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: Input validation tests
# ============================================================================ #
test_that("JAGS_add_priors input validation works", {

  # Empty prior_list returns original syntax
  expect_equal(JAGS_add_priors("model{}", list()), "model{}")

  # prior_list must be a list of priors
  expect_error(JAGS_add_priors("model{}", list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_add_priors("model{}", prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that("JAGS_get_inits input validation works", {

  # Empty prior_list returns empty list
  expect_equal(JAGS_get_inits(list(), chains = 2, seed = 1), list())

  # Input validation
  expect_error(JAGS_get_inits(list(x = 1), chains = 2, seed = 1), "'prior_list' must be a list of priors.")
  expect_error(JAGS_get_inits(prior("normal", list(0, 1)), chains = 2, seed = 1), "'prior_list' must be a list of priors.")

})


test_that("JAGS_to_monitor input validation works", {

  # Empty prior_list returns empty string
  expect_equal(JAGS_to_monitor(list()), "")

  # Input validation
  expect_error(JAGS_to_monitor(list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_to_monitor(prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that(".check_JAGS_syntax validates syntax correctly", {

  # Test with valid syntax
  expect_silent(JAGS_add_priors("model{}", list(mu = prior("normal", list(0, 1)))))
  expect_equal(JAGS_add_priors(NULL, list()), "model{}")
  expect_match(JAGS_add_priors(NULL, list(mu = prior("normal", list(0, 1)))), "^model\\{")

  # Test with missing "model" keyword
  expect_error(
    JAGS_add_priors("invalid{}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing opening brace
  expect_error(
    JAGS_add_priors("model}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing closing brace
  expect_error(
    JAGS_add_priors("model{", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with non-character input
  expect_error(
    JAGS_add_priors(123, list(mu = prior("normal", list(0, 1)))),
    "must be a character"
  )

})


test_that("JAGS_extend error handling", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Test error when fit is not a BayesTools_fit
  expect_error(
    JAGS_extend(list(), autofit_control = list()),
    "'fit' must be a 'BayesTools_fit'"
  )

})


test_that("required packages are checked locally and on parallel workers", {

  expect_silent({
    local_loaded <- .JAGS_require_packages("stats")
  })
  expect_equal(unname(local_loaded), TRUE)
  expect_error(
    .JAGS_require_packages("BayesToolsMissingPackageForTest"),
    "Required packages are not available: 'BayesToolsMissingPackageForTest'.",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit("model{}", required_packages = NA_character_),
    "The 'required_packages' argument cannot contain NA/NaN values.",
    fixed = TRUE
  )

  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  expect_silent({
    worker_loaded <- .JAGS_require_packages("stats", cl)
  })
  expect_equal(unname(worker_loaded), TRUE)
  expect_error(
    .JAGS_require_packages(c("stats", "BayesToolsMissingPackageForTest"), cl),
    "Required packages are not available: 'BayesToolsMissingPackageForTest'.",
    fixed = TRUE
  )

})

test_that("JAGS_fit rejects unrecognized formula_scale_list entries", {

  skip_if_not_installed("runjags")

  data <- data.frame(x = c(-1, 0, 1))
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = ~ x),
      formula_data_list = list(mu = data),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )),
      formula_scale_list = list(muu = list(x = TRUE))
    ),
    "not recognized"
  )
})

test_that("JAGS_fit stores formula design metadata on fitted formula models", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7),
    x = c(-1, 0, 1, 2)
  )

  fit <- JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 1)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ x),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 123
  )

  design <- JAGS_formula_design(fit, "mu")

  expect_s3_class(fit, "BayesTools_fit")
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(colnames(design$model_matrix), c("(Intercept)", "x"))
  expect_equal(unname(design$model_matrix[, "x"]), df$x)
  expect_equal(design$jags_data_names$x, "mu_data_x")
  expect_identical(JAGS_formula_design(fit), list(mu = design))
})

test_that("JAGS_fit does not emit unused sampled data for marginalized random effects", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.3, -0.1),
    study = factor(c("s1", "s1", "s2", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"))
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  prior_random <- prior_random(
    study = random_block(sd = sd_prior),
    estimate = random_block(sd = sd_prior)
  )
  warnings <- character()

  fit <- withCallingHandlers(
    JAGS_fit(
      model_syntax = "model{
        for(i in 1:N_mu){
          y[i] ~ dnorm(mu[i], 1)
        }
      }",
      data = list(y = df$y),
      formula_list = list(
        mu = ~ 1 +
          random(1 | study, name = "study", covariance = "diag") +
          random(1 | estimate, name = "estimate", covariance = "diag")
      ),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 1))
      )),
      formula_random_prior_list = list(mu = prior_random),
      formula_random_effects_compile_list = list(
        mu = random_effects_compile(marginalized = "estimate")
      ),
      chains = 1,
      adapt = 100,
      burnin = 100,
      sample = 100,
      silent = FALSE,
      seed = 123
    ),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  design <- JAGS_formula_design(fit, "mu")
  estimate_term <- design$random_effects[[which(vapply(
    design$random_effects,
    function(random_term) identical(random_term$block_name, "estimate"),
    logical(1)
  ))]]

  expect_equal(warnings, character())
  expect_equal(estimate_term$compile_mode, "marginalized")
  expect_equal(estimate_term$jags_data_names, character())
  expect_false(grepl("mu__xREx__estimate_xRE_DATAx", attr(fit, "model"), fixed = TRUE))
  expect_false(grepl("mu__xREx__estimate_xRE_MAPx", attr(fit, "model"), fixed = TRUE))
})

test_that("JAGS_fit runs dummy structured random-effect formula models", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  sd_prior <- prior("normal", list(0, 0.5), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  model_syntax <- "model{
    for(i in 1:N_mu){
      y[i] ~ dnorm(mu[i], 4)
    }
  }"

  df_cs <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3, 0.5, 0.9),
    x = c(-1, 0, 1, 2, -2, 0.5, 1.5, 2.5),
    idx = factor(rep(c("t1", "t2"), 4), levels = c("t1", "t2")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 2))
  )

  fit_cs <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_cs$y),
    formula_list = list(mu = ~ 1 + x + cs(idx | id)),
    formula_data_list = list(mu = df_cs),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = rho_prior,
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 101
  ))

  expect_s3_class(fit_cs, "BayesTools_fit")
  expect_equal(attr(fit_cs, "jags_modules"), character())
  expect_equal(JAGS_formula_design(fit_cs, "mu")$random_effects[[1]]$structure, "cs")
  expect_true("mu__xREx__id_rho" %in% colnames(as.matrix(fit_cs$mcmc)))

  df_ar1 <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3, 0.5, 0.9, 0.2, 0.6, 0.8, 1.0),
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )

  fit_ar1 <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_ar1$y),
    formula_list = list(mu = ~ 1 + ar1(f | id)),
    formula_data_list = list(mu = df_ar1),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = rho_prior,
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 102
  ))

  expect_s3_class(fit_ar1, "BayesTools_fit")
  ar1_design <- JAGS_formula_design(fit_ar1, "mu")$random_effects[[1]]
  expect_equal(ar1_design$structure, "ar1")
  expect_equal(ar1_design$column_names, c("fa", "fb", "fc"))
  expect_true("mu__xREx__id_rho" %in% colnames(as.matrix(fit_ar1$mcmc)))

  df_us <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    x = c(-1, 0, 1, 2, -2, 0.5),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit_us <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_us$y),
    formula_list = list(mu = ~ 1 + x + (1 + x | id)),
    formula_data_list = list(mu = df_us),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 1, include_correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 103
  ))

  expect_s3_class(fit_us, "BayesTools_fit")
  expect_equal(attr(fit_us, "jags_modules"), "BayesTools")
  expect_equal(JAGS_formula_design(fit_us, "mu")$random_effects[[1]]$structure, "us")
  expect_true(any(grepl("mu__xREx__id_xRE_CORx_L", colnames(as.matrix(fit_us$mcmc)), fixed = TRUE)))
})

test_that("JAGS_fit monitors random coefficients for observed-level prediction", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + diag(1 | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(coefficients = TRUE, correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 103
  ))

  posterior <- as.matrix(fit$mcmc)
  expect_true(any(grepl("mu__xREx__id_xRE_COEFx", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list")
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS_fit predicts observed random effects from latent monitors", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + diag(1 | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(coefficients = FALSE, correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 104
  ))

  posterior <- as.matrix(fit$mcmc)
  expect_true(any(grepl("mu__xREx__id_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("mu__xREx__id_xRE_COEFx", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list")
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("centered and noncentered fits preserve substantive output schemas", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.4, -0.1, 0.2, 0.5, 0.8, 0.1, 0.3, 0.6, 0.9, 1.1,
          -0.3, 0, 0.1, 0.4, 0.7, 0.2, 0.5, 0.7, 1, 1.2),
    index = factor(rep(c("i1", "i2"), 10L)),
    id = factor(rep(paste0("g", seq_len(4L)), each = 5L))
  )
  fit_parameterization <- function(parameterization, seed){

    suppressWarnings(JAGS_fit(
      model_syntax = "model{
        for(i in 1:N_mu){
          y[i] ~ dnorm(mu[i], 4)
        }
      }",
      data = list(y = df$y),
      formula_list = list(mu = ~ 1 + cs(index | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          rho = prior("normal", list(0, 0.5)),
          parameterization = parameterization
        )
      )),
      chains = 1,
      adapt = 100,
      burnin = 200,
      sample = 400,
      silent = TRUE,
      seed = seed
    ))
  }

  noncentered <- fit_parameterization("noncentered", 207)
  centered    <- fit_parameterization("centered", 208)
  noncentered_draws <- as.matrix(noncentered$mcmc)
  centered_draws    <- as.matrix(centered$mcmc)
  centered_term <- JAGS_formula_design(centered, "mu")$random_effects[[1L]]
  sd_name       <- centered_term$sd_parameter_names[[1L]]
  rho_name      <- centered_term$correlation$rho_name

  expect_setequal(colnames(centered_draws), colnames(noncentered_draws))
  expect_equal(
    mean(centered_draws[, "mu_intercept"]),
    mean(noncentered_draws[, "mu_intercept"]),
    tolerance = 0.2
  )
  expect_equal(
    mean(centered_draws[, sd_name]),
    mean(noncentered_draws[, sd_name]),
    tolerance = 0.25
  )
  expect_equal(
    mean(centered_draws[, rho_name]),
    mean(noncentered_draws[, rho_name]),
    tolerance = 0.25
  )
  expect_identical(centered_term$parameterization_resolved, "centered")
})

test_that("group-local structured fits monitor only active latent cells", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  K <- 40L
  df <- data.frame(
    y = seq(-0.5, 0.5, length.out = K),
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(8L)), length.out = K))
  )
  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + cs(index | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        rho = prior("point", list(location = 0.2))
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 209
  ))
  term <- JAGS_formula_design(fit, "mu")$random_effects[[1L]]
  posterior <- as.matrix(fit$mcmc)
  z_names <- grep("_xRE_Zx[", colnames(posterior), fixed = TRUE, value = TRUE)

  expect_s3_class(
    term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_length(z_names, K)
  expect_lt(length(z_names), term$n_groups * term$n_columns)
  expect_false(any(grepl("_xRE_CORx_L[", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("_xRE_CORx_R[", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + cs(index | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list")
  )
  expect_equal(dim(prediction), c(K, nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS transformed scalar rho remains representably inside support", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, -0.1, 0.2, 0.5),
    index = factor(rep(c("i1", "i2", "i3"), 2L)),
    id = factor(rep(c("g1", "g2"), each = 3L))
  )
  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + ar1(index | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        rho = prior("point", list(location = 1e300))
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 210
  ))
  random_term <- JAGS_formula_design(fit, "mu")$random_effects[[1L]]
  interior <- .bt_random_effect_representable_rho_bounds(
    random_term$correlation$bounds,
    random_term$structure
  )
  rho <- as.matrix(fit$mcmc)[, random_term$correlation$rho_name]

  expect_true(all(is.finite(rho)))
  expect_equal(unname(rho), rep(interior[["upper"]], length(rho)), tolerance = 0)
  expect_lt(max(rho), random_term$correlation$bounds[["upper"]])
})

test_that("JAGS_fit predicts row-indexed external SD random effects from latent monitors", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    tau_data = c(0.5, 0.75, 1.0, 1.25)
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
        tau[i] <- tau_data[i]
      }
    }",
    data = list(y = df$y, tau_data = df$tau_data),
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      allocation = random_variance_allocation(
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )),
    add_parameters = c("tau", "mu"),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 105
  ))

  posterior <- as.matrix(fit$mcmc)
  tau_names <- paste0("tau[", seq_len(nrow(df)), "]")
  mu_names <- paste0("mu[", seq_len(nrow(df)), "]")

  expect_true(all(tau_names %in% colnames(posterior)))
  expect_true(all(mu_names %in% colnames(posterior)))
  expect_true(any(grepl("mu__xREx__study_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_true(any(grepl("mu__xREx__drug_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("_xRE_COEFx", colnames(posterior), fixed = TRUE)))
  expect_equal(
    unname(posterior[, tau_names, drop = FALSE]),
    unname(matrix(
      rep(df$tau_data, each = nrow(posterior)),
      nrow = nrow(posterior),
      ncol = nrow(df)
    )),
    tolerance = 1e-12
  )

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = formula,
    parameter = "mu",
    data = df,
    fitted_rows = seq_len(nrow(df)),
    prior_list = attr(fit, "prior_list")
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
  expect_equal(
    unname(t(prediction)),
    unname(posterior[, mu_names, drop = FALSE]),
    tolerance = 1e-8
  )
})

test_that("JAGS_fit stores formula design metadata on failed sampling objects", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(x = c(-1, 0, 1))

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      broken_node <- missing_node
    }",
    formula_list = list(mu = ~ x),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 123
  ))

  design <- JAGS_formula_design(fit, "mu")

  expect_s3_class(fit, "BayesTools_fit")
  expect_s3_class(fit, "error")
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(unname(design$model_matrix[, "x"]), df$x)
})


# ============================================================================ #
# SECTION 2: Convergence edge cases
# ============================================================================ #
test_that("autofit settings keep indicator checks off by default", {

  settings <- JAGS_check_and_list_autofit_settings(list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10
  ))
  expect_false(settings$check_indicators)

  settings <- JAGS_check_and_list_autofit_settings(list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10,
    check_indicators = TRUE
  ))
  expect_true(settings$check_indicators)

})


test_that("JAGS_check_convergence ignores indicator variables unless requested", {

  set.seed(1)
  chain_1 <- cbind(mu = rnorm(100), mu_indicator = rep(0, 100))
  chain_2 <- cbind(mu = rnorm(100), mu_indicator = rep(1, 100))
  fit <- list(
    mcmc         = coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2)),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"

  prior_list <- list(mu = prior("normal", list(0, 1)))

  expect_true(JAGS_check_convergence(
    fit,
    prior_list       = prior_list,
    max_Rhat        = 1.05,
    min_ESS         = NULL,
    max_error       = NULL,
    max_SD_error    = NULL,
    check_indicators = FALSE
  ))

  with_indicators <- JAGS_check_convergence(
    fit,
    prior_list       = prior_list,
    max_Rhat        = 1.05,
    min_ESS         = NULL,
    max_error       = NULL,
    max_SD_error    = NULL,
    check_indicators = TRUE
  )
  expect_false(with_indicators)
  expect_match(attr(with_indicators, "errors"), "R-hat")

})

test_that("JAGS_check_convergence ignores add_parameters without priors", {

  set.seed(2)
  chain_1 <- cbind(mu = rnorm(100), "aux[1]" = rep(0, 100))
  chain_2 <- cbind(mu = rnorm(100), "aux[1]" = rep(1, 100))
  fit <- list(
    mcmc         = coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2)),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"

  prior_list <- list(mu = prior("normal", list(0, 1)))

  without_aux <- JAGS_check_convergence(
    fit,
    prior_list    = prior_list,
    max_Rhat     = 1.05,
    min_ESS      = NULL,
    max_error    = NULL,
    max_SD_error = NULL
  )
  expect_false(without_aux)

  expect_true(JAGS_check_convergence(
    fit,
    prior_list     = prior_list,
    max_Rhat      = 1.05,
    min_ESS       = NULL,
    max_error     = NULL,
    max_SD_error  = NULL,
    add_parameters = "aux"
  ))
})


test_that("JAGS_check_convergence marks single-chain R-hat as not assessable", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  old_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.runjags = old_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.runjags = TRUE)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 1,  # Single chain - R-hat cannot be computed
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent.jags = TRUE
  ))

  convergence <- expect_silent(
    JAGS_check_convergence(
      fit,
      prior_list = prior_list,
      max_Rhat = 1.05,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL
    )
  )
  expect_false(convergence)
  expect_match(attr(convergence, "errors"), "R-hat.*not assessable")

})


test_that("JAGS_check_convergence handles ESS and error checks", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  old_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.runjags = old_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.runjags = TRUE)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 2,
    adapt = 50,
    burnin = 50,
    sample = 50,  # Small sample for testing convergence failures
    silent.jags = TRUE
  ))

  # Test with very strict ESS requirement (should fail)
  result_ess <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = 10000, max_error = NULL, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_ess)
  expect_true(!is.null(attr(result_ess, "errors")))

  # Test with very strict error requirement
  result_err <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = 0.00001, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_err)

  # Test with very strict SD error requirement
  result_sd <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = 0.00001, fail_fast = FALSE)
  expect_false(result_sd)

})


# ============================================================================ #
# SECTION 3: JAGS_fit with is_JASP mode
# ============================================================================ #
test_that("JAGS_fit works with is_JASP mode", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Simple model for testing is_JASP mode
  set.seed(1)
  data <- list(
    y = rnorm(20, 0.5, 1),
    N = 20
  )

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1), list(0, Inf))
  )

  model_syntax <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu, 1/pow(sigma, 2))
    }
  }"

  # Mock JASP progress bar functions (they should be skipped if not available)
  # The is_JASP mode should work but simply skip progress bars if functions don't exist
  fit_jasp <- capture.output(tryCatch({
    suppressWarnings(JAGS_fit(
      model_syntax = model_syntax,
      data = data,
      prior_list = prior_list,
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      seed = 1,
      silent = TRUE,
      is_JASP = TRUE,
      is_JASP_prefix = "Test"
    ))
  }, error = function(e) {
    # If JASP functions don't exist, this should still produce a fit
    # or fail gracefully
    if (grepl("JASP", e$message)) {
      skip("JASP progress bar functions not available")
    }
    stop(e)
  }))

  test_reference_text(paste0(fit_jasp, collapse = ","), "fit_jasp.txt")

})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-fit-lm-oracles.R
# ============================================================================ #

skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS Fit LM Oracles
# ============================================================================ #
#
# PURPOSE:
#   Fit-profile semantic oracle tests for Gaussian JAGS formula models. These
#   compare automatic formula scaling and manual standardization against lm()
#   predictions on the original data scale, and compare a known-sigma Gaussian
#   model against its closed-form posterior.
#
# TAGS: @fit, @JAGS, @formula, @standardization
# ============================================================================ #

skip_on_cran()
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

.fit_gaussian_formula_oracle <- function(data, formula_data, formula_scale_list = NULL, seed = 1L) {
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], pow(sigma, -2))\n",
    "}\n",
    "}"
  )

  JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = data$y, N = nrow(data)),
    prior_list = list(
      sigma = prior("lognormal", list(0, 1))
    ),
    formula_list = list(mu = ~ x1 * x2),
    formula_data_list = list(mu = formula_data),
    formula_prior_list = list(
      mu = list(
        "intercept" = prior("normal", list(0, 10)),
        "x1" = prior("normal", list(0, 5)),
        "x2" = prior("normal", list(0, 5)),
        "x1:x2" = prior("normal", list(0, 5))
      )
    ),
    formula_scale_list = formula_scale_list,
    chains = 2,
    adapt = 250,
    burnin = 500,
    sample = 1200,
    seed = seed,
    silent = TRUE
  )
}

.fit_known_sigma_gaussian_formula_oracle <- function(data, sigma_known, seed = 1L) {
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], pow(sigma_known, -2))\n",
    "}\n",
    "}"
  )

  JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = data$y, N = nrow(data), sigma_known = sigma_known),
    prior_list = NULL,
    formula_list = list(mu = ~ x1 * x2),
    formula_data_list = list(mu = data[c("x1", "x2")]),
    formula_prior_list = list(
      mu = list(
        "intercept" = prior("normal", list(0, 10)),
        "x1" = prior("normal", list(0, 5)),
        "x2" = prior("normal", list(0, 5)),
        "x1:x2" = prior("normal", list(0, 5))
      )
    ),
    formula_scale_list = list(mu = list(x1 = TRUE, x2 = TRUE)),
    chains = 2,
    adapt = 500,
    burnin = 1000,
    sample = 3000,
    seed = seed,
    silent = TRUE
  )
}

test_that("Gaussian JAGS formula fit agrees with lm oracle after scaling", {
  data <- bayestools_oracle_gaussian_regression_data()
  manual_formula_data <- bayestools_manual_scaled_data(data[c("x1", "x2")], c("x1", "x2"))
  manual_scale <- attr(manual_formula_data, "manual_scale")

  fit_manual <- .fit_gaussian_formula_oracle(
    data = data,
    formula_data = manual_formula_data,
    formula_scale_list = NULL,
    seed = 7701L
  )
  attr(fit_manual, "manual_scale") <- manual_scale

  fit_auto <- .fit_gaussian_formula_oracle(
    data = data,
    formula_data = data[c("x1", "x2")],
    formula_scale_list = list(mu = list(x1 = TRUE, x2 = TRUE)),
    seed = 7701L
  )

  expect_formula_scale_equal(attr(fit_auto, "formula_scale")$mu, manual_scale)

  scaled_parameters <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")
  posterior_manual <- as.matrix(suppressWarnings(coda::as.mcmc(fit_manual)))
  posterior_auto <- as.matrix(suppressWarnings(coda::as.mcmc(fit_auto)))

  expect_posterior_summary_close(
    posterior_auto[, scaled_parameters, drop = FALSE],
    colMeans(posterior_manual[, scaled_parameters, drop = FALSE]),
    tolerance = 0.05
  )

  posterior_original <- transform_scale_samples(
    posterior_auto[, scaled_parameters, drop = FALSE],
    attr(fit_auto, "formula_scale")
  )
  posterior_manual_original <- transform_scale_samples(
    posterior_manual[, scaled_parameters, drop = FALSE],
    list(mu = manual_scale)
  )

  lm_fit <- stats::lm(y ~ x1 * x2, data = data)
  lm_coef <- stats::coef(lm_fit)
  names(lm_coef) <- bayestools_lm_coef_to_jags_names(names(lm_coef))

  expect_posterior_summary_close(
    posterior_original[, names(lm_coef), drop = FALSE],
    lm_coef,
    tolerance = 0.30
  )

  newdata <- data.frame(
    x1 = stats::quantile(data$x1, probs = c(0.10, 0.50, 0.90), names = FALSE),
    x2 = stats::quantile(data$x2, probs = c(0.20, 0.60, 0.80), names = FALSE)
  )
  prediction_design <- stats::model.matrix(~ x1 * x2, data = newdata)
  prediction_names <- bayestools_lm_coef_to_jags_names(colnames(prediction_design))
  prediction_auto <- as.numeric(prediction_design %*% colMeans(
    posterior_original[, prediction_names, drop = FALSE]
  ))
  prediction_manual <- as.numeric(prediction_design %*% colMeans(
    posterior_manual_original[, prediction_names, drop = FALSE]
  ))
  expect_equal(prediction_auto, prediction_manual, tolerance = 0.15)

  expect_lm_predictions_equal(
    lm_fit = lm_fit,
    posterior = posterior_original,
    newdata = newdata,
    formula = ~ x1 * x2,
    tolerance = 0.30
  )

  sigma_samples <- posterior_auto[, "sigma"]
  expect_equal(mean(sigma_samples), summary(lm_fit)$sigma, tolerance = 0.25)

  posterior_intervals <- apply(
    posterior_original[, names(lm_coef), drop = FALSE],
    2,
    stats::quantile,
    probs = c(0.005, 0.995)
  )
  expect_true(all(lm_coef >= posterior_intervals[1, ] & lm_coef <= posterior_intervals[2, ]))
})

test_that("known-sigma Gaussian JAGS formula fit matches closed-form posterior", {
  data <- bayestools_oracle_gaussian_regression_data()
  sigma_known <- 0.6
  parameters <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")

  fit <- .fit_known_sigma_gaussian_formula_oracle(
    data = data,
    sigma_known = sigma_known,
    seed = 7711L
  )

  posterior <- as.matrix(suppressWarnings(coda::as.mcmc(fit)))
  posterior <- posterior[, parameters, drop = FALSE]

  scale_info <- attr(fit, "formula_scale")$mu
  expected_scale <- bayestools_manual_scaled_data(data[c("x1", "x2")], c("x1", "x2"))
  expect_formula_scale_equal(scale_info, attr(expected_scale, "manual_scale"))

  design <- stats::model.matrix(~ x1 * x2, data = expected_scale)
  jags_design <- JAGS_formula_design(fit, "mu")$model_matrix
  expect_equal(colnames(design), c("(Intercept)", "x1", "x2", "x1:x2"))
  expect_equal(colnames(jags_design), c("(Intercept)", "x1", "x2", "x1__xXx__x2"))
  expect_equal(unname(jags_design), unname(design), tolerance = 1e-12)
  expect_setequal(colnames(posterior), parameters)

  oracle <- bayestools_gaussian_posterior_oracle(
    X = design,
    y = data$y,
    sigma = sigma_known,
    prior_mean = c(0, 0, 0, 0),
    prior_sd = c(10, 5, 5, 5)
  )
  names(oracle$mean) <- parameters
  names(oracle$sd) <- parameters

  expect_equal(colMeans(posterior), oracle$mean, tolerance = 0.07)
  expect_equal(apply(posterior, 2, stats::sd), oracle$sd, tolerance = 0.04)
  expect_equal(unname(stats::cov(posterior)), unname(oracle$cov), tolerance = 0.02)

  posterior_intervals <- apply(posterior, 2, stats::quantile, probs = c(0.01, 0.99))
  expect_true(all(oracle$mean >= posterior_intervals[1, ] & oracle$mean <= posterior_intervals[2, ]))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-lkj-cholesky-fit.R
# ============================================================================ #

skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS LKJ-Cholesky Fit Oracles
# ============================================================================ #
#
# PURPOSE:
#   Fit-profile MCMC checks for the package-shipped JAGS LKJ-Cholesky module.
#
# TAGS: @fit, @JAGS, @Stan, @LKJ, @Cholesky
# ============================================================================ #

skip_on_cran()
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

.fit_jags_lkj_cholesky_prior <- function(K, eta, sample = 3000, seed = 1L) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = K,
    eta = eta,
    include_correlation = TRUE,
    include_primitives = TRUE
  )

  user_silent.jags <- runjags::runjags.getOption("silent.jags")
  user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  fit <- suppressWarnings(runjags::run.jags(
    model = paste0("model{\n", module$syntax, "\n}"),
    monitor = module$monitor,
    n.chains = 2,
    adapt = 500,
    burnin = 500,
    sample = sample,
    thin = 1,
    method = "rjags",
    summarise = FALSE,
    plots = FALSE,
    inits = lapply(seq_len(2), function(i){
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  ))

  list(
    fit = fit,
    module = module,
    samples = as.matrix(fit$mcmc)
  )
}

.fit_jags_lkj_cholesky_prior_cached <- local({
  cache <- new.env(parent = emptyenv())

  function(K, eta, sample = 3000, seed = 1L) {
    key <- paste(K, eta, sample, seed, sep = "|")
    if(!exists(key, envir = cache, inherits = FALSE)){
      assign(
        key,
        .fit_jags_lkj_cholesky_prior(
          K = K,
          eta = eta,
          sample = sample,
          seed = seed
        ),
        envir = cache
      )
    }

    get(key, envir = cache, inherits = FALSE)
  }
})

.eval_jags_lkj_cholesky_transform <- function(u, K) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  K2 <- K * K
  syntax <- c(
    paste0("  Omega_L_flat[1:", K2, "] <- bt_lkj_cholesky(u, ", K, ")"),
    paste0("  Omega_R_flat[1:", K2, "] <- bt_lkj_corr(u, ", K, ")")
  )
  for(i in seq_len(K)){
    for(j in seq_len(K)){
      flat_index <- BayesTools:::.bt_lkj_cholesky_flat_index(i, j, K)
      syntax <- c(
        syntax,
        paste0("  Omega_L[", i, ",", j, "] <- Omega_L_flat[", flat_index, "]"),
        paste0("  Omega_R[", i, ",", j, "] <- Omega_R_flat[", flat_index, "]")
      )
    }
  }

  con <- textConnection(paste0("model{\n", paste(syntax, collapse = "\n"), "\n}\n"))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file = con,
    data = list(u = u),
    n.chains = 1,
    n.adapt = 0,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = c("Omega_L", "Omega_R"),
    n.iter = 1,
    quiet = TRUE,
    progress.bar = "none"
  )

  as.matrix(samples)
}

.jags_lkj_vector_column <- function(samples, base_name, index) {
  indexed_name <- paste0(base_name, "[", index, "]")
  if(indexed_name %in% colnames(samples)){
    return(indexed_name)
  }
  if(index == 1L && base_name %in% colnames(samples)){
    return(base_name)
  }
  stop("Missing monitored JAGS column: ", indexed_name, call. = FALSE)
}

.jags_lkj_matrix_draw <- function(samples, row, prefix, K) {
  out <- matrix(NA_real_, K, K)
  for(i in seq_len(K)){
    for(j in seq_len(K)){
      out[i, j] <- samples[row, paste0(prefix, "[", i, ",", j, "]")]
    }
  }
  out
}

.jags_lkj_offdiag_columns <- function(samples, prefix, K) {
  out <- character(0)
  for(row in 2:K){
    for(column in seq_len(row - 1L)){
      out <- c(out, paste0(prefix, "[", column, ",", row, "]"))
    }
  }
  if(!all(out %in% colnames(samples))){
    stop("Missing monitored JAGS correlation columns.", call. = FALSE)
  }

  out
}

.eval_jags_lkj_cpc_deviance <- function(u, alpha) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  con <- textConnection(paste0(
    "model{\n",
    "  u[1:", length(u), "] ~ dbt_lkj_cpc(alpha)\n",
    "}\n"
  ))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file = con,
    data = list(u = u, alpha = alpha),
    n.chains = 2,
    n.adapt = 0,
    quiet = TRUE
  )
  dic <- rjags::dic.samples(
    model = model,
    n.iter = 1,
    type = "pD",
    quiet = TRUE,
    progress.bar = "none"
  )

  list(
    deviance = as.numeric(dic[["deviance"]]),
    penalty = as.numeric(dic[["penalty"]])
  )
}

.lkj_bivariate_data <- function(N = 28L, rho = 0.55, seed = 4711L) {
  set.seed(seed)
  z_1 <- stats::rnorm(N)
  z_2 <- rho * z_1 + sqrt(1 - rho^2) * stats::rnorm(N)
  cbind(z_1, z_2)
}

.lkj_bivariate_loglik <- function(rho, y) {
  denom <- 1 - rho^2
  quad <- (y[, 1]^2 - 2 * rho * y[, 1] * y[, 2] + y[, 2]^2) / denom
  sum(-0.5 * log(denom) - 0.5 * quad)
}

.lkj_bivariate_grid_oracle <- function(y, eta, grid_size = 20001L) {
  rho <- seq(-.999, .999, length.out = grid_size)
  log_density <- vapply(rho, .lkj_bivariate_loglik, numeric(1), y = y) +
    (eta - 1) * log1p(-rho^2)
  log_density <- log_density - max(log_density)
  weights <- exp(log_density)
  weights <- weights / sum(weights)
  cdf <- cumsum(weights)
  posterior_mean <- sum(weights * rho)

  list(
    mean = posterior_mean,
    sd = sqrt(sum(weights * (rho - posterior_mean)^2)),
    quantiles = stats::approx(cdf, rho, xout = c(.1, .5, .9), ties = "ordered")$y
  )
}

.fit_jags_lkj_bivariate_posterior <- function(y, eta, sample = 5000,
                                              seed = 812L) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  model <- paste0(
    "model{\n",
    "  alpha[1] <- eta\n",
    "  u[1:1] ~ dbt_lkj_cpc(alpha)\n",
    "  R_flat[1:4] <- bt_lkj_corr(u, 2)\n",
    "  rho <- R_flat[2]\n",
    "  for(i in 1:N){\n",
    "    loglik[i] <- -0.5 * log(1 - pow(rho, 2)) - ",
    "(pow(y[i,1], 2) - 2 * rho * y[i,1] * y[i,2] + pow(y[i,2], 2)) / ",
    "(2 * (1 - pow(rho, 2)))\n",
    "    zeros[i] ~ dpois(C - loglik[i])\n",
    "  }\n",
    "}\n"
  )

  user_silent.jags <- runjags::runjags.getOption("silent.jags")
  user_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.jags = user_silent.jags, silent.runjags = user_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  fit <- suppressWarnings(runjags::run.jags(
    model = model,
    data = list(y = y, eta = eta, N = nrow(y), zeros = rep(0, nrow(y)), C = 100),
    monitor = "rho",
    n.chains = 2,
    adapt = 1000,
    burnin = 1000,
    sample = sample,
    thin = 1,
    method = "rjags",
    summarise = FALSE,
    plots = FALSE,
    inits = lapply(seq_len(2), function(i){
      list(u = c(0.5), .RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  ))

  as.numeric(as.matrix(fit$mcmc)[, "rho"])
}

test_that("JAGS LKJ-Cholesky module matches primitive and correlation marginals", {
  settings <- data.frame(
    K = c(2L, 3L, 4L, 5L),
    eta = c(0.35, 1, 2.5, 0.75),
    sample = c(3000L, 3000L, 3000L, 4000L)
  )

  for(setting_i in seq_len(nrow(settings))){
    K <- settings$K[setting_i]
    eta <- settings$eta[setting_i]
    result <- .fit_jags_lkj_cholesky_prior_cached(
      K = K,
      eta = eta,
      sample = settings$sample[setting_i],
      seed = 100 + setting_i
    )
    samples <- result$samples
    pairs <- result$module$pairs

    u_columns <- character(nrow(pairs))
    for(p in seq_len(nrow(pairs))){
      u_name <- .jags_lkj_vector_column(samples, "Omega_lkj_u", p)
      u_columns[p] <- u_name
      alpha <- pairs$alpha[p]
      expect_equal(mean(samples[, u_name]), 0.5, tolerance = 0.05)
      expect_true(abs(stats::var(samples[, u_name]) - 1 / (4 * (2 * alpha + 1))) < 0.012)
      expect_equal(
        stats::quantile(samples[, u_name], c(.1, .5, .9), names = FALSE),
        stats::qbeta(c(.1, .5, .9), alpha, alpha),
        tolerance = 0.06
      )
    }

    if(length(u_columns) >= 2L){
      primitive_cor <- stats::cor(samples[, u_columns, drop = FALSE])
      expect_lt(max(abs(primitive_cor[upper.tri(primitive_cor)])), 0.10)
    }

    rho_columns <- .jags_lkj_offdiag_columns(samples, "Omega_R", K)
    alpha_rho <- eta - 1 + K / 2
    rho_quantiles <- 2 * stats::qbeta(c(.1, .5, .9), alpha_rho, alpha_rho) - 1
    rho_means <- vapply(rho_columns, function(name) mean(samples[, name]), numeric(1))
    rho_sds <- vapply(rho_columns, function(name) stats::sd(samples[, name]), numeric(1))
    for(rho_name in rho_columns){
      expect_equal(mean(samples[, rho_name]), 0, tolerance = 0.045)
      expect_true(abs(stats::var(samples[, rho_name]) - 1 / (2 * alpha_rho + 1)) < 0.025)
      expect_equal(
        stats::quantile(samples[, rho_name], c(.1, .5, .9), names = FALSE),
        rho_quantiles,
        tolerance = 0.08
      )
    }
    expect_lt(max(rho_means) - min(rho_means), 0.08)
    expect_lt(max(rho_sds) - min(rho_sds), 0.08)

    check_rows <- unique(round(seq(1, nrow(samples), length.out = min(25, nrow(samples)))))
    for(row in check_rows){
      L <- .jags_lkj_matrix_draw(samples, row, "Omega_L", K)
      R <- .jags_lkj_matrix_draw(samples, row, "Omega_R", K)

      expect_equal(L[upper.tri(L)], rep(0, K * (K - 1) / 2), tolerance = 1e-10)
      expect_true(all(diag(L) > 0))
      expect_equal(rowSums(L^2), rep(1, K), tolerance = 1e-8)
      expect_equal(R, L %*% t(L), tolerance = 1e-8)
      expect_equal(diag(R), rep(1, K), tolerance = 1e-8)
      expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
    }
  }
})

test_that("compiled LKJ-Cholesky functions match hand-coded transform oracles", {
  settings <- list(
    list(
      K = 2L,
      u = c(0.31),
      expected_L = matrix(c(
        1, 0,
        -0.38, sqrt(1 - .38^2)
      ), nrow = 2, byrow = TRUE)
    ),
    list(
      K = 3L,
      u = c(0.61, 0.22, 0.74),
      expected_L = matrix(c(
        1, 0, 0,
        .22, sqrt(1 - .22^2), 0,
        -.56, .48 * sqrt(1 - .56^2), sqrt(1 - .56^2) * sqrt(1 - .48^2)
      ), nrow = 3, byrow = TRUE)
    )
  )

  for(setting in settings){
    samples <- .eval_jags_lkj_cholesky_transform(setting$u, setting$K)
    expected_L <- setting$expected_L
    expected_R <- expected_L %*% t(expected_L)
    observed_L <- .jags_lkj_matrix_draw(samples, 1, "Omega_L", setting$K)
    observed_R <- .jags_lkj_matrix_draw(samples, 1, "Omega_R", setting$K)

    expect_equal(observed_L, expected_L, tolerance = 1e-12)
    expect_equal(observed_R, expected_R, tolerance = 1e-12)
    expect_equal(observed_L[upper.tri(observed_L)], rep(0, setting$K * (setting$K - 1L) / 2L), tolerance = 1e-12)
    expect_true(all(diag(observed_L) > 0))
    expect_equal(rowSums(observed_L^2), rep(1, setting$K), tolerance = 1e-12)
    expect_equal(diag(observed_R), rep(1, setting$K), tolerance = 1e-12)
  }
})

test_that("compiled LKJ transforms evaluate boundary proposals", {
  settings <- list(
    list(K = 2L, u = -1e-8),
    list(K = 2L, u = 0),
    list(K = 2L, u = 1),
    list(K = 2L, u = 1 + 1e-8),
    list(K = 3L, u = c(0, 1, 0.5))
  )

  for(setting in settings){
    samples <- .eval_jags_lkj_cholesky_transform(setting$u, setting$K)
    L <- .jags_lkj_matrix_draw(samples, 1, "Omega_L", setting$K)
    R <- .jags_lkj_matrix_draw(samples, 1, "Omega_R", setting$K)

    expect_true(all(is.finite(L)))
    expect_true(all(is.finite(R)))
    expect_equal(R, L %*% t(L), tolerance = 1e-12)
    expect_equal(diag(R), rep(1, setting$K), tolerance = 1e-12)
  }
})

test_that("compiled LKJ CPC distribution logDensity matches beta density for observed nodes", {
  K <- 4L
  eta <- 0.75
  u <- c(0.08, 0.91, 0.47, 0.12, 0.63, 0.88)
  alpha <- BayesTools:::.bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)$alpha
  deviance <- .eval_jags_lkj_cpc_deviance(u = u, alpha = alpha)
  expected_deviance <- -2 * sum(stats::dbeta(u, alpha, alpha, log = TRUE))

  expect_equal(deviance$deviance, expected_deviance, tolerance = 1e-12)
  expect_equal(deviance$penalty, 0, tolerance = 1e-12)
})

test_that("JAGS LKJ posterior for bivariate normal correlation matches grid oracle", {
  eta <- 0.8
  y <- .lkj_bivariate_data()
  oracle <- .lkj_bivariate_grid_oracle(y = y, eta = eta)
  rho <- .fit_jags_lkj_bivariate_posterior(y = y, eta = eta, sample = 5000)

  expect_equal(mean(rho), oracle$mean, tolerance = 0.06)
  expect_equal(stats::sd(rho), oracle$sd, tolerance = 0.06)
  expect_equal(
    stats::quantile(rho, c(.1, .5, .9), names = FALSE),
    oracle$quantiles,
    tolerance = 0.09
  )
})

test_that("JAGS_fit loads BayesTools module from generated LKJ metadata", {
  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = 2,
    eta = 1.1,
    include_correlation = TRUE
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = paste0("model{\n", module$syntax, "\n}"),
    prior_list = NULL,
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    add_parameters = module$monitor,
    required_packages = module$required_packages,
    jags_modules = module$jags_module,
    silent = TRUE,
    seed = 10
  ))

  expect_s3_class(fit, "BayesTools_fit")
  expect_equal(attr(fit, "jags_modules"), "BayesTools")
  samples <- as.matrix(fit$mcmc)
  expect_true(all(c("Omega_L[1,1]", "Omega_R[1,2]") %in% colnames(samples)))
})

test_that("one-dimensional generated LKJ module runs without native LKJ calls", {
  module <- JAGS_lkj_corr_cholesky(
    name = "One",
    K = 1,
    eta = 3,
    include_correlation = TRUE
  )

  expect_false(grepl("dbt_lkj_cpc", module$syntax, fixed = TRUE))
  expect_false(grepl("bt_lkj_cholesky", module$syntax, fixed = TRUE))
  expect_false(grepl("bt_lkj_corr", module$syntax, fixed = TRUE))

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = paste0("model{\n", module$syntax, "\n}"),
    prior_list = NULL,
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    add_parameters = module$monitor,
    required_packages = module$required_packages,
    jags_modules = module$jags_module,
    silent = TRUE,
    seed = 11
  ))

  expect_s3_class(fit, "BayesTools_fit")
  samples <- as.matrix(fit$mcmc)
  expect_equal(colnames(samples), c("One_L", "One_R"))
  expect_equal(unname(samples[, "One_L"]), rep(1, nrow(samples)))
  expect_equal(unname(samples[, "One_R"]), rep(1, nrow(samples)))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-marglik.R
# ============================================================================ #

skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS Marginal Likelihood Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for JAGS marginal likelihood computation functions.
#   Uses simple models where the log marginal likelihood is known to be 0
#   (for prior samples, the marginal likelihood for any proper prior is 1).
#
# DEPENDENCIES:
#   - rjags: For JAGS model fitting
#   - bridgesampling: For marginal likelihood computation
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags")
#   - Note: Creates fresh models, does not need pre-fitted models
#
# MODELS/FIXTURES:
#   - Creates models with known analytical marginal likelihoods for validation
#
# TAGS: @evaluation, @JAGS, @marginal-likelihood
# ============================================================================ #

# Load common test helpers
source(testthat::test_path("common-functions.R"))

make_bridge_random_fixture <- function(formula = ~ 1 + x + us(1 + x | id),
                                       data = NULL,
                                       prior_list = NULL,
                                       prior_random_list = NULL){

  if(is.null(data)){
    data <- data.frame(
      x = c(-1, 0, 1, 2),
      id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
    )
  }
  if(is.null(prior_list)){
    prior_list <- list(intercept = prior("normal", list(0, 1)))
    if("x" %in% all.vars(formula)){
      prior_list$x <- prior("normal", list(0, 1))
    }
  }
  if(is.null(prior_random_list)){
    prior_random_list <- prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1, include_correlation = FALSE)
      )
    )
  }

  list(
    data = data,
    formula = formula,
    prior_list = prior_list,
    prior_random_list = prior_random_list,
    result = JAGS_formula(
      formula = formula,
      parameter = "mu",
      data = data,
      prior_list = prior_list,
      prior_random = prior_random_list
    )
  )
}

expect_formula_random_prior_only_bridge <- function(formula, data, prior_list,
                                                    prior_random_list,
                                                    n_iter = 8000,
                                                    tolerance = 0.08,
                                                    seed = 1,
                                                    maxiter = 2000){

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(seed)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = seed),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = n_iter,
    quiet = TRUE,
    progress.bar = "none"
  )
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = maxiter
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = tolerance)

  invisible(list(
    marglik = marglik,
    formula_result = formula_result,
    samples = samples
  ))
}

test_that("direct composed bias priors support bridge-sampling helpers", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))
  phacking  <- prior_phacking(form = "linear", alpha = prior("beta", list(2, 3)))
  bias      <- prior_bias(selection = selection, phacking = phacking)

  posterior <- matrix(
    c(
      1.5, 2.5, .2,
      1.1, 2.1, .4
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(posterior) <- c("eta[1]", "eta[2]", "alpha")

  prepared <- JAGS_bridgesampling_posterior(posterior, list(bias = bias))
  expect_equal(colnames(prepared), c("eta[1]", "eta[2]", "alpha"))
  expect_equal(attr(prepared, "lb"), c("eta[1]" = 0, "eta[2]" = 0, "alpha" = 0))
  expect_equal(attr(prepared, "ub"), c("eta[1]" = Inf, "eta[2]" = Inf, "alpha" = 1))

  samples <- posterior[1, ]
  expected_prior_density <-
    sum(stats::dgamma(samples[c("eta[1]", "eta[2]")], shape = c(1, 2), rate = 1, log = TRUE)) +
    stats::dbeta(samples[["alpha"]], 2, 3, log = TRUE)
  expect_equal(
    JAGS_marglik_priors(samples, list(bias = bias)),
    expected_prior_density,
    tolerance = 1e-12
  )

  parameters <- JAGS_marglik_parameters(samples, list(bias = bias))
  constants <- phack_backend_constants(phacking$form, phacking$source, phacking$destination, target = phacking$target)
  expect_equal(parameters$omega, c(1, samples[["eta[2]"]] / sum(samples[c("eta[1]", "eta[2]")])))
  expect_equal(parameters$alpha, samples[["alpha"]])
  expect_equal(
    parameters$pi_null,
    samples[["alpha"]] * constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    parameters$beta_null,
    samples[["alpha"]] * constants$beta_null_per_alpha,
    tolerance = 1e-12
  )
})

test_that("p-hacking bridge helpers support point and inverse-gamma alpha priors", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))

  point_alpha <- prior_phacking(
    form  = "linear",
    alpha = prior("point", list(.25))
  )
  point_bias <- prior_bias(selection = selection, phacking = point_alpha)
  point_posterior <- matrix(
    c(
      1.5, 2.5,
      1.1, 2.1
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(point_posterior) <- c("eta[1]", "eta[2]")

  point_prepared <- JAGS_bridgesampling_posterior(point_posterior, list(bias = point_bias))
  expect_equal(colnames(point_prepared), c("eta[1]", "eta[2]"))
  expect_equal(attr(point_prepared, "lb"), c("eta[1]" = 0, "eta[2]" = 0))
  expect_equal(attr(point_prepared, "ub"), c("eta[1]" = Inf, "eta[2]" = Inf))

  point_samples <- point_posterior[1, ]
  expect_equal(
    JAGS_marglik_priors(point_samples, list(bias = point_bias)),
    sum(stats::dgamma(point_samples[c("eta[1]", "eta[2]")], shape = c(1, 2), rate = 1, log = TRUE)),
    tolerance = 1e-12
  )

  point_parameters <- JAGS_marglik_parameters(point_samples, list(bias = point_bias))
  point_constants <- phack_backend_constants(
    point_alpha$form, point_alpha$source, point_alpha$destination,
    target = point_alpha$target
  )
  expect_equal(point_parameters$alpha, .25)
  expect_equal(
    point_parameters$pi_null,
    .25 * point_constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    point_parameters$beta_null,
    .25 * point_constants$beta_null_per_alpha,
    tolerance = 1e-12
  )

  invgamma_alpha_prior <- prior("invgamma", list(shape = 3, scale = .4), list(0, 1))
  invgamma_alpha <- prior_phacking(
    form  = "linear",
    alpha = invgamma_alpha_prior
  )
  invgamma_bias <- prior_bias(selection = selection, phacking = invgamma_alpha)
  invgamma_posterior <- matrix(
    c(
      1.5, 2.5, 0.4,
      1.1, 2.1, 0.5
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(invgamma_posterior) <- c("eta[1]", "eta[2]", "alpha")

  invgamma_prepared <- JAGS_bridgesampling_posterior(invgamma_posterior, list(bias = invgamma_bias))
  expect_equal(colnames(invgamma_prepared), c("eta[1]", "eta[2]", "alpha"))
  expect_equal(attr(invgamma_prepared, "lb"), c("eta[1]" = 0, "eta[2]" = 0, "alpha" = 0))
  expect_equal(attr(invgamma_prepared, "ub"), c("eta[1]" = Inf, "eta[2]" = Inf, "alpha" = 1))

  invgamma_samples <- invgamma_posterior[1, ]
  expected_invgamma_prior_density <-
    sum(stats::dgamma(invgamma_samples[c("eta[1]", "eta[2]")], shape = c(1, 2), rate = 1, log = TRUE)) +
    lpdf(invgamma_alpha_prior, invgamma_samples[["alpha"]])
  expect_equal(
    JAGS_marglik_priors(invgamma_samples, list(bias = invgamma_bias)),
    expected_invgamma_prior_density,
    tolerance = 1e-12
  )

  invgamma_parameters <- JAGS_marglik_parameters(invgamma_samples, list(bias = invgamma_bias))
  invgamma_constants <- phack_backend_constants(
    invgamma_alpha$form, invgamma_alpha$source, invgamma_alpha$destination,
    target = invgamma_alpha$target
  )
  expect_equal(invgamma_parameters$alpha, invgamma_samples[["alpha"]])
  expect_equal(
    invgamma_parameters$pi_null,
    invgamma_samples[["alpha"]] * invgamma_constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    invgamma_parameters$beta_null,
    invgamma_samples[["alpha"]] * invgamma_constants$beta_null_per_alpha,
    tolerance = 1e-12
  )

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma test.
  legacy_invgamma_posterior <- invgamma_posterior
  colnames(legacy_invgamma_posterior)[3] <- "inv_alpha"
  legacy_invgamma_posterior[, "inv_alpha"] <- 1 / legacy_invgamma_posterior[, "inv_alpha"]
  legacy_prepared <- JAGS_bridgesampling_posterior(legacy_invgamma_posterior, list(bias = invgamma_bias))
  expect_equal(colnames(legacy_prepared), c("eta[1]", "eta[2]", "alpha"))
  expect_equal(legacy_prepared[, "alpha"], invgamma_posterior[, "alpha"])
})

test_that("bias mixtures fail explicitly in bridge-sampling helpers", {

  bias <- prior_mixture(list(
    prior_none(),
    prior_phacking(form = "linear")
  ))
  samples <- c("bias_indicator" = 1, "alpha" = .2)
  posterior <- matrix(samples, nrow = 1)

  expect_error(
    JAGS_bridgesampling_posterior(posterior, list(bias = bias)),
    "bias mixture priors"
  )
  expect_error(
    JAGS_marglik_priors(samples, list(bias = bias)),
    "bias mixture priors"
  )
  expect_error(
    JAGS_marglik_parameters(samples, list(bias = bias)),
    "bias mixture priors"
  )
})

# This file tests the JAGS marginal likelihood computation functions
# It uses simple models where the log marginal likelihood is known to be 0
# (for prior samples, the marginal likelihood for any proper prior is 1, log(1) = 0)
# More complex consistency tests (e.g., including formulas etc part of `test-00-model-fits.R`)

test_that("JAGS model functions work (simple)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("normal", list(0, 1), list(1, Inf)),
    p3  = prior("lognormal", list(0, .5)),
    p4  = prior("t", list(0, .5, 5)),
    p5  = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    p6  = prior("gamma", list(2, 1)),
    p7  = prior("invgamma", list(3, 2), list(1, 3)),
    p8  = prior("exp", list(1.5)),
    p9  = prior("beta", list(3, 2)),
    p10 = prior("uniform", list(1, 5)),
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
    #p12 = prior("bernoulli", list(0.75)) discrete priors are not supported with bridgesampling
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

# skip the rest as it takes too long
skip_on_cran()

test_that("JAGS model functions work (vector)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("mnormal", list(mean = 0, sd = 1, K = 3),),
    p2  = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    p3  = prior("mt",      list(location = 2, scale = 0.5, df = 5, K = 2))
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 5*1e-2) # the mCauchy is a bit more variable
  }

})

test_that("JAGS model functions work (factor)", {

  skip_if_not_installed("rjags")
  all_priors   <- list(
    p1  = prior_factor("mnorm", list(mean = 0, sd = 1),    contrast = "orthonormal"),
    p2  = prior_factor("beta",  list(alpha = 1, beta = 1), contrast = "treatment"),
    p3  = prior_factor("beta",  list(alpha = 2, beta = 2), contrast = "treatment"),
    p4  = prior_factor("gamma",   list(shape = 2, rate = 3), contrast = "independent"),
    p5  = prior_factor("uniform", list(a = -0.5, b = 1.5),   contrast = "independent"),
    p6  = prior_factor("mnorm", list(mean = 0, sd = 1),     contrast = "meandif")
  )

  # add levels
  attr(all_priors[[1]], "levels") <- 3
  attr(all_priors[[2]], "levels") <- 2
  attr(all_priors[[3]], "levels") <- 3
  attr(all_priors[[4]], "levels") <- 1
  attr(all_priors[[5]], "levels") <- 3
  attr(all_priors[[6]], "levels") <- 3
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS marginal-likelihood helpers reject spike-and-slab priors explicitly", {

  prior_list <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    )
  )
  posterior <- matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "theta"))
  samples <- c(theta = 0, theta_inclusion = .5)

  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list),
    "spike and slab priors is not implemented"
  )
  expect_error(
    JAGS_marglik_priors(samples, prior_list),
    "prior mixture priors is not implemented"
  )
  expect_error(
    JAGS_marglik_parameters(samples, prior_list),
    "prior mixture priors is not implemented"
  )
})

test_that("JAGS model functions work (weightfunctions)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1))),
    prior_weightfunction("one-sided", c(.05, 0.10), wf_cumulative(c(1, 2, 3))),
    prior_weightfunction("one-sided", c(.05, 0.60), wf_independent(prior("beta", list(1, 1)))),
    prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)))
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    names(prior_list) <- "omega"
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (spikes)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1    = prior("spike", list(1)),
    p2.2  = prior_factor("spike", list(location = 2), contrast = "treatment"),
    p3.2  = prior_factor("spike", list(location = 3), contrast = "independent"),
    p4.2  = prior_factor("spike", list(location = 0), contrast = "orthonormal"),
    p5.2  = prior_factor("spike", list(location = 0), contrast = "meandif"),
    p2.5  = prior_factor("spike", list(location = 2), contrast = "treatment"),
    p3.5  = prior_factor("spike", list(location = 3), contrast = "independent"),
    p4.5  = prior_factor("spike", list(location = 0), contrast = "orthonormal"),
    p5.5  = prior_factor("spike", list(location = 0), contrast = "meandif")
  )
  attr(all_priors$p2.2, "levels") <- 2
  attr(all_priors$p3.2, "levels") <- 2
  attr(all_priors$p4.2, "levels") <- 2
  attr(all_priors$p5.2, "levels") <- 2
  attr(all_priors$p2.5, "levels") <- 2
  attr(all_priors$p3.5, "levels") <- 2
  attr(all_priors$p4.5, "levels") <- 2
  attr(all_priors$p5.5, "levels") <- 2
  nuisance_prior <- list(sigma = prior("normal", list(0, 1)))
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- c(all_priors[i], nuisance_prior)
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("bridge sampling object function works",{

  marglik0 <- bridgesampling_object()
  marglik1 <- bridgesampling_object(1)

  expect_equal(marglik0$logml, -Inf)
  expect_equal(marglik1$logml, 1)
  expect_s3_class(marglik0, "BayesTools_marglik")

})

test_that("JAGS marglik with formula works", {

  # Test marginal likelihood computation with formula interface
  # Uses intercept-only formula with various priors
  # When sampling from prior and computing marglik, the result should be ~0 (log(1))

  skip_if_not_installed("rjags")

  # Simple data for the formula
  set.seed(1)
  df_test <- data.frame(x = rnorm(10))
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Create formula prior list with intercept only
  prior_list <- list(
    "intercept" = prior("gamma",  list(2, 2)),
    "x"         = prior("normal", list(0, 1))
  )

  # Process formula to get JAGS syntax
  formula_result <- JAGS_formula(~ 1 + x, parameter = "mu", data = df_test, prior_list = prior_list)

  # Build JAGS model with formula priors
  model_syntax <- JAGS_add_priors("model{}", formula_result$prior_list)
  monitor      <- JAGS_to_monitor(formula_result$prior_list)
  inits        <- JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1)

  # Sample from prior using JAGS
  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  # Compute marginal likelihood using formula interface
  marglik <- JAGS_bridgesampling(
    fit                = samples,
    log_posterior      = log_posterior,
    data               = list(),
    prior_list         = NULL,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = list(mu = df_test),
    formula_prior_list = list(mu = prior_list)
  )

  expect_equal(marglik$logml, 0, tolerance = 1e-3)
})

test_that("JAGS marglik with exp(intercept) formula works", {

  # Test marginal likelihood computation with formula interface
  # Uses intercept-only formula with various priors
  # When sampling from prior and computing marglik, the result should be ~0 (log(1))

  skip_if_not_installed("rjags")

  # Simple data for the formula
  set.seed(1)
  df_test <- data.frame(x = rnorm(10))
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Create formula prior list with intercept only
  prior_list <- list(
    "intercept" = prior("gamma",  list(2, 2)),
    "x"         = prior("normal", list(0, 1))
  )

  # Process formula to get JAGS syntax
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- TRUE
  formula_result <- JAGS_formula(formula, parameter = "mu", data = df_test, prior_list = prior_list)
  expect_equal(formula_result$formula_syntax, "for(i in 1:N_mu){\n  mu[i] = log(mu_intercept) + mu_x * mu_data_x[i]\n}\n")

  # Build JAGS model with formula priors
  model_syntax <- JAGS_add_priors("model{}", formula_result$prior_list)
  monitor      <- JAGS_to_monitor(formula_result$prior_list)
  inits        <- JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1)

  # Sample from prior using JAGS
  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  # Compute marginal likelihood using formula interface
  marglik <- JAGS_bridgesampling(
    fit                = samples,
    log_posterior      = log_posterior,
    data               = list(),
    prior_list         = NULL,
    formula_list       = list(mu = formula),
    formula_data_list  = list(mu = df_test),
    formula_prior_list = list(mu = prior_list)
  )

  expect_equal(marglik$logml, 0, tolerance = 1e-3)
})

test_that("JAGS bridgesampling infers formula scaling metadata from fits", {

  df_test <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x"         = prior("normal", list(0, 1))
  )

  scaled_formula <- JAGS_formula(
    ~ 1 + x,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  fit <- matrix(c(0, 1), ncol = 2)
  attr(fit, "formula_scale") <- list(mu = scaled_formula$formula_scale)

  inferred_scale <- BayesTools:::.JAGS_formula_scale_list_from_fit(fit, "mu")
  expect_equal(inferred_scale, list(mu = list(x = TRUE)))

  bridge_formula <- JAGS_formula(
    ~ 1 + x,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    formula_scale = inferred_scale$mu
  )

  expect_equal(bridge_formula$data$mu_data_x, scaled_formula$data$mu_data_x)
})

test_that("JAGS bridgesampling posterior supports add-only parameters", {

  posterior <- matrix(
    c(
      -1, 0.2,
       0, 0.5,
       1, 0.8
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(posterior) <- c("x", "prob")

  info <- BayesTools:::.JAGS_bridgesampling_posterior_info(NULL)
  expect_type(info, "character")
  expect_equal(length(info), 0L)
  expect_equal(attr(info, "lb"), numeric())
  expect_equal(attr(info, "ub"), numeric())

  result_null <- JAGS_bridgesampling_posterior(
    posterior = posterior,
    prior_list = NULL,
    add_parameters = c("x", "prob"),
    add_bounds = list(
      lb = c(x = -Inf, prob = 0),
      ub = c(x = Inf, prob = 1)
    )
  )
  expect_equal(colnames(result_null), c("x", "prob"))
  expect_equal(attr(result_null, "lb"), c(x = -Inf, prob = 0))
  expect_equal(attr(result_null, "ub"), c(x = Inf, prob = 1))

  result_empty <- JAGS_bridgesampling_posterior(
    posterior = posterior,
    prior_list = list(),
    add_parameters = "x",
    add_bounds = list(lb = c(x = -Inf), ub = c(x = Inf))
  )
  expect_equal(colnames(result_empty), "x")
  expect_equal(attr(result_empty, "lb"), c(x = -Inf))
  expect_equal(attr(result_empty, "ub"), c(x = Inf))
})

test_that("JAGS bridgesampling passes requested bridge context to callback", {

  skip_if_not_installed("bridgesampling")
  skip_if_not_installed("coda")

  set.seed(1)
  posterior <- matrix(
    rnorm(2000),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )
  posterior <- coda::as.mcmc(posterior)
  seen <- new.env(parent = emptyenv())
  log_posterior <- function(parameters, data, bridge_context){
    if(!exists("context", envir = seen, inherits = FALSE)){
      seen$context <- bridge_context
      seen$mu <- parameters$mu
    }
    expect_s3_class(bridge_context, "BayesTools_bridge_context")
    expect_true("mu" %in% names(bridge_context$state))
    expect_equal(bridge_context$nodes[["mu"]], parameters$mu)
    0
  }

  marglik <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = log_posterior,
    data = list(),
    prior_list = list(mu = prior("normal", list(0, 1))),
    bridge_context = TRUE,
    maxiter = 1000
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_s3_class(seen$context, "BayesTools_bridge_context")
  expect_equal(seen$context$nodes[["mu"]], seen$mu)
})

test_that("JAGS bridgesampling validates rebuilt formula random design metadata", {

  fixture <- make_bridge_random_fixture()
  fitted <- list(mu = fixture$result$formula_design)
  rebuilt <- list(mu = fixture$result$formula_design)

  expect_silent(BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, rebuilt))

  changed <- rebuilt
  changed$mu$random_effects[[1]]$structure <- "diag"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "covariance structure",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$structure <- NULL
  changed$mu$random_effects[[1]]$covariance <- "us"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "missing canonical 'random_term\\$structure'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$homogeneous_sd <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$correlation <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "missing canonical 'random_term\\$correlation'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$block_name <- "other"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "block names",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$group_levels <- rev(changed$mu$random_effects[[1]]$group_levels)
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "group levels",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$column_names[1] <- "changed"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "column names",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$n_columns <- changed$mu$random_effects[[1]]$n_columns + 1L
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "group or column counts",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$model_matrix <- changed$mu$random_effects[[1]]$model_matrix[-1, , drop = FALSE]
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "model matrix shape or columns",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$group_map <- rev(changed$mu$random_effects[[1]]$group_map)
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(fitted, changed),
    "group map",
    fixed = TRUE
  )

  structured_fixture <- make_bridge_random_fixture(
    formula = ~ 1 + ar1(idx | id),
    data = data.frame(
      idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
      id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  structured_fitted <- list(mu = structured_fixture$result$formula_design)
  structured_changed <- structured_fitted
  structured_changed$mu$random_effects[[1]]$structured_index$label <- "changed"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      structured_fitted,
      structured_changed
    ),
    "structured random-effect index metadata",
    fixed = TRUE
  )

  allocation_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )
  allocation_formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  allocation_prior_list <- list(intercept = prior("normal", list(0, 1)))
  allocation_fixture <- make_bridge_random_fixture(
    formula = allocation_formula,
    data = allocation_data,
    prior_list = allocation_prior_list,
    prior_random_list = prior_random(
      allocation = random_variance_allocation(
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  independent_fixture <- make_bridge_random_fixture(
    formula = allocation_formula,
    data = allocation_data,
    prior_list = allocation_prior_list,
    prior_random_list = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2))),
      drug = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      list(mu = allocation_fixture$result$formula_design),
      list(mu = independent_fixture$result$formula_design)
    ),
    "scale/allocation",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling errors on fitted/rebuilt random design mismatches", {

  fixture <- make_bridge_random_fixture()
  fit <- coda::mcmc(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "dummy")))
  attr(fit, "formula_design") <- list(mu = fixture$result$formula_design)

  mismatch_data <- fixture$data
  mismatch_data$id <- factor(as.character(mismatch_data$id), levels = c("b", "a"))

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = NULL,
      formula_list = list(mu = fixture$formula),
      formula_data_list = list(mu = mismatch_data),
      formula_prior_list = list(mu = fixture$prior_list),
      formula_random_prior_list = list(mu = fixture$prior_random_list),
      maxiter = 10
    ),
    "original formula source data differ",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling supports formula random effects through prior_random", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  if(!isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE))){
    skip("BayesTools JAGS module is not available")
  }

  df_test <- data.frame(
    x = c(-0.5, 0.5),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 + x + us(1 + x | id)
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior_lkj(eta = 1, include_correlation = FALSE)
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_true("mu__xREx__id_xRE_CORx_lkj_u[1]" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 5000,
    quiet = TRUE,
    progress.bar = "none"
  )
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling supports continuous-time CAR formula random effects", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 + car(0 + time | id)
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      rho = prior("normal", list(0, 0.5))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_equal(formula_result$formula_design$random_effects[[1]]$structure, "car")
  expect_equal(formula_result$formula_design$random_effects[[1]]$correlation$bounds, c(lower = 0, upper = 1))
  expect_true("mu__xREx__id_rho" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("centered continuous-time CAR syntax samples sequential conditionals", {

  skip_if_not_installed("rjags")

  df_test <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = df_test,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        rho = prior("normal", list(0, 0.5)),
        monitor = random_monitor(coefficients = TRUE),
        parameterization = "centered"
      )
    )
  )

  expect_false(grepl(
    "dmnorm.vcov",
    formula_result$formula_syntax,
    fixed = TRUE
  ))
  expect_match(
    formula_result$formula_syntax,
    "pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[2], 1)",
    fixed = TRUE
  )

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters,
    "mu__xREx__id_xRE_CAR_PHIX",
    "mu__xREx__id_xRE_CAR_INNOV_VARx"
  ))
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 1, seed = 31),
    n.chains = 1,
    n.adapt = 100,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 250,
    quiet = TRUE,
    progress.bar = "none"
  )
  draws <- as.matrix(samples[[1L]])
  coefficient_columns <- grep(
    "mu__xREx__id_xRE_COEFx[",
    colnames(draws),
    fixed = TRUE
  )
  transition_columns <- c(
    grep(
      "mu__xREx__id_xRE_CAR_PHIX[",
      colnames(draws),
      fixed = TRUE
    ),
    grep(
      "mu__xREx__id_xRE_CAR_INNOV_VARx[",
      colnames(draws),
      fixed = TRUE
    )
  )

  expect_length(coefficient_columns, 6L)
  expect_length(transition_columns, 4L)
  expect_true(all(is.finite(draws[, coefficient_columns, drop = FALSE])))
  expect_true(all(is.finite(draws[, transition_columns, drop = FALSE])))
  expect_true(all(apply(
    draws[, coefficient_columns, drop = FALSE],
    2L,
    stats::sd
  ) > 0))
  expect_true(all(apply(
    draws[, transition_columns, drop = FALSE],
    2L,
    stats::sd
  ) > 0))
})

test_that("JAGS bridgesampling supports Dirichlet variance-allocation random effects", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    allocation = random_variance_allocation(
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_true("mu__xRE_ALLOCx_allocation__weight" %in% names(formula_result$prior_list))
  expect_true("mu__xREx__study_xRE_Zx" %in% formula_result$add_parameters)
  expect_true("mu__xREx__drug_xRE_Zx" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling reconstructs row-indexed external SD sources from values", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    tau_factor = c(0.5, 0.75, 1.0, 1.25)
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  tau_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      data$tau_factor[seq_len(n_rows)]
    }
  )
  prior_random_list <- prior_random(
    allocation = random_variance_allocation(
      sd_source = random_sd_source(tau_source),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )
  no_values_prior_random_list <- prior_random(
    allocation = random_variance_allocation(
      sd_source = random_sd_source("tau", shape = "row"),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  no_values_formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = no_values_prior_random_list
  )
  model_syntax <- JAGS_add_priors(
    paste0(
      "model{\n",
      "for(i in 1:N_mu){\n",
      "  tau[i] = tau_factor[i]\n",
      "}\n",
      formula_result$formula_syntax,
      "\n}"
    ),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(11)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = c(formula_result$data, list(tau_factor = df_test$tau_factor)),
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 11),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  expect_false(any(grepl("^tau\\[", colnames(as.matrix(samples)))))

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(tau_factor = df_test$tau_factor),
    prior_list = NULL,
    maxiter = 1000
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)

  graft_samples <- samples
  attr(graft_samples, "formula_design") <- list(mu = no_values_formula_result$formula_design)
  expect_error(
    JAGS_bridgesampling(
      fit = graft_samples,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = NULL,
      formula_list = list(mu = formula),
      formula_data_list = list(mu = df_test),
      formula_prior_list = list(mu = prior_list),
      formula_random_prior_list = list(mu = prior_random_list),
      maxiter = 1000
    ),
    "scale/allocation metadata differ",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling gives unit marglik for prior-only random-effect settings", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  continuous_data <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"))
  )
  factor_data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 3), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 3), levels = c("g1", "g2", "g3"))
  )

  diag_heterogeneous <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + diag(1 + x | id, hom = FALSE),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        terms = list(x = prior("gamma", list(3, 1)))
      )
    ),
    seed = 10
  )
  expect_false(diag_heterogeneous$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    diag_heterogeneous$formula_result$formula_design$random_effects[[1]]$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_x")
  )

  id_homogeneous <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + id(1 + x | id),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        monitor = random_monitor(
          latent = TRUE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    ),
    seed = 11
  )
  expect_true(id_homogeneous$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    unique(id_homogeneous$formula_result$formula_design$random_effects[[1]]$sd_parameter_names),
    "mu__xREx__id_sd"
  )

  lkj_module <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + us(1 + x | id),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(
          eta = 2,
          include_correlation = FALSE
        )
      )
    ),
    seed = 12
  )
  expect_equal(lkj_module$formula_result$jags_modules, "BayesTools")
  expect_true(any(grepl("_xRE_CORx_lkj_u", lkj_module$formula_result$add_parameters, fixed = TRUE)))

  cs_fisher_z <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + cs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = prior("normal", list(0, 0.5))
      )
    ),
    seed = 13
  )
  expect_true(cs_fisher_z$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    cs_fisher_z$formula_result$formula_design$random_effects[[1]]$correlation$rho_scale,
    "fisher_z"
  )
  expect_true("mu__xREx__id_rho_z" %in% names(cs_fisher_z$formula_result$prior_list))

  hcs_logit <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + hcs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5)),
          rho_scale = "logit"
        )
      )
    ),
    seed = 14
  )
  expect_false(hcs_logit$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_true("mu__xREx__id_rho_logit" %in% names(hcs_logit$formula_result$prior_list))
  expect_equal(
    hcs_logit$formula_result$formula_design$random_effects[[1]]$correlation$rho_scale,
    "logit"
  )

  ar1_fisher_z <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + ar1(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = prior("normal", list(0, 0.5))
      )
    ),
    seed = 15
  )
  expect_equal(ar1_fisher_z$formula_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_true("mu__xREx__id_rho_z" %in% names(ar1_fisher_z$formula_result$prior_list))

  har_raw <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + har(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5), truncation = list(lower = -1, upper = 1)),
          rho_scale = "rho"
        )
      )
    ),
    seed = 16
  )
  expect_equal(har_raw$formula_result$formula_design$random_effects[[1]]$structure, "har")
  expect_false(har_raw$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_true("mu__xREx__id_rho" %in% names(har_raw$formula_result$prior_list))
  expect_false("mu__xREx__id_rho_z" %in% names(har_raw$formula_result$prior_list))

  nested_allocation <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    data = data.frame(
      study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
      paper = factor(c("p1", "p2", "p1", "p2"), levels = c("p1", "p2")),
      drug = factor(c("d1", "d1", "d2", "d2"), levels = c("d1", "d2"))
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        weights = prior("dirichlet", list(alpha = c(3, 2)))
      )
    ),
    n_iter = 10000,
    seed = 17
  )
  expect_true("mu__xRE_ALLOCx_total_re__weight" %in% names(nested_allocation$formula_result$prior_list))
  expect_true("mu__xRE_ALLOCx_nested_split__weight" %in% names(nested_allocation$formula_result$prior_list))
  expect_equal(
    length(nested_allocation$formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors),
    2L
  )

  sd_leaf_allocation <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + hcs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      allocation = random_variance_allocation(
        name = "leaf_alloc",
        terms = "id",
        target = "sd_component",
        scale = "mean_variance",
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = c(2, 3, 4)))
      ),
      id = random_block(rho = prior("normal", list(0, 0.5)))
    ),
    n_iter = 10000,
    seed = 18
  )
  expect_true("mu__xRE_ALLOCx_leaf_alloc__weight" %in% names(sd_leaf_allocation$formula_result$prior_list))
  expect_equal(
    sd_leaf_allocation$formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$target,
    "sd_component"
  )
  expect_equal(
    sd_leaf_allocation$formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$scale,
    "mean_variance"
  )
})

test_that("JAGS formula marglik reconstructs inverse-gamma terms on natural scale", {

  samples <- c(
    "mu_intercept" = 0.5,
    "mu_x"         = 0.25
  )
  formula_data_list <- list(
    mu = list(
      N_mu      = 2,
      mu_data_x = c(10, 20)
    )
  )
  formula_prior_list <- list(
    mu = list(
      mu_intercept = prior("invgamma", list(2, 1)),
      mu_x         = prior("invgamma", list(2, 1))
    )
  )

  parameters <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(parameters$mu, c(0.5 + 0.25 * 10, 0.5 + 0.25 * 20))

  formula_log_intercept <- ~ 1 + x
  attr(formula_log_intercept, "log(intercept)") <- TRUE
  parameters_log <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = formula_log_intercept),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(parameters_log$mu, c(log(0.5) + 0.25 * 10, log(0.5) + 0.25 * 20))

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma test.
  legacy_samples <- c(
    "inv_mu_intercept" = 2,
    "inv_mu_x"         = 4
  )
  legacy_parameters <- JAGS_marglik_parameters_formula(
    samples            = legacy_samples,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(legacy_parameters$mu, parameters$mu)
})


test_that("JAGS marglik reconstructs indexed factor inverse-gamma parameters", {
  theta_prior <- prior_factor("invgamma", list(2, 1), contrast = "independent")
  attr(theta_prior, "levels") <- 2

  parameters <- BayesTools:::.JAGS_marglik_parameters.factor(
    samples = c("theta[1]" = 0.5, "theta[2]" = 0.25),
    prior = theta_prior,
    parameter_name = "theta"
  )

  expect_equal(parameters$theta, c(0.5, 0.25))

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma test.
  legacy_parameters <- BayesTools:::.JAGS_marglik_parameters.factor(
    samples = c("inv_theta[1]" = 2, "inv_theta[2]" = 4),
    prior = theta_prior,
    parameter_name = "theta"
  )

  expect_equal(legacy_parameters$theta, parameters$theta)
})


test_that("JAGS formula marglik preserves predictor names containing _data", {
  samples <- c(
    "mu_intercept" = 1,
    "mu_x_data"   = 2
  )
  formula_data_list <- list(
    mu = list(
      N_mu           = 2,
      mu_data_x_data = c(10, 20)
    )
  )
  formula_prior_list <- list(
    mu = list(
      mu_intercept = prior("normal", list(0, 1)),
      mu_x_data    = prior("normal", list(0, 1))
    )
  )

  parameters <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = ~ 1 + x_data),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list()
  )

  expect_equal(parameters$mu, c(1 + 2 * 10, 1 + 2 * 20))
})


# Targeted tests for uncovered code paths in JAGS-marglik.R

test_that("JAGS_bridgesampling_posterior input validation works", {

  posterior <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(posterior) <- c("mu", "sigma", "x")

  # Input validation errors

  expect_error(JAGS_bridgesampling_posterior(data.frame(x = 1), prior_list = NULL), "'posterior' must be a matrix")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = "x"), "'prior_list' must be a list.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = 1), "'add_parameters' must be a character")
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_bounds = list(lb = -Inf, ub = Inf)
    ),
    "requires at least one 'add_parameters'",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = character(),
      add_bounds = list(lb = numeric(), ub = numeric())
    ),
    "requires at least one 'add_parameters'",
    fixed = TRUE
  )
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = "x"), "'add_bounds' must be a list")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(a = 1)), "'add_bounds' must contain lower and upper bounds")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = c("x", "y"), add_bounds = list(lb = 0, ub = 1)), "'lb' and 'ub' must have the same length")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(lb = "a", ub = "b")), "'lb' and 'ub' must be numeric")
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = -Inf, ub = Inf)
    ),
    "names must be unique and match",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = stats::setNames(-Inf, "wrong"), ub = stats::setNames(Inf, "x"))
    ),
    "names must match",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = c(x = 1), ub = c(x = 0))
    ),
    "smaller than upper",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = list(mu = prior("normal", list(0, 1))),
      add_parameters = "mu",
      add_bounds = list(lb = c(mu = -Inf), ub = c(mu = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = list(sigma = prior("invgamma", list(2, 1))),
      add_parameters = "sigma",
      add_bounds = list(lb = c(sigma = 0), ub = c(sigma = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_add_parameters_not_formula(
      add_parameters = "mu",
      formula_design_list = list(mu = list()),
      formula_prior_list = list()
    ),
    "BayesTools-owned formula parameter",
    fixed = TRUE
  )

  dirichlet_posterior <- matrix(
    c(1, 3, 0.25),
    nrow = 1,
    dimnames = list(NULL, c("prior_par_eta_w[1]", "prior_par_eta_w[2]", "w[1]"))
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      dirichlet_posterior,
      prior_list = list(w = prior("dirichlet", list(alpha = c(1, 1)))),
      add_parameters = "w[1]",
      add_bounds = list(lb = c("w[1]" = 0), ub = c("w[1]" = 1))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )

  # Unsupported prior types
  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list = list(p1 = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1))))),
    "spike and slab"
  )
  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list = list(p1 = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1))), is_null = c(TRUE, FALSE)))),
    "prior mixture"
  )

  # Missing parameters
  posterior_small <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(posterior_small) <- c("a", "b")
  expect_error(JAGS_bridgesampling_posterior(posterior_small, prior_list = list(x = prior("normal", list(0, 1)))), "'posterior' does not contain all")

  # Successful case with add_parameters
  result <- JAGS_bridgesampling_posterior(posterior, prior_list = list(mu = prior("normal", list(0, 1))), add_parameters = "x", add_bounds = list(lb = c(x = -Inf), ub = c(x = Inf)))
  expect_true(is.matrix(result))
  expect_true("x" %in% colnames(result))

})

test_that("JAGS_marglik_priors input validation and edge cases work", {

  # Empty prior_list contributes zero log prior density

  expect_equal(JAGS_marglik_priors(list(), prior_list = list()), 0)
  expect_equal(JAGS_marglik_priors_formula(list(), formula_prior_list = list(mu = list())), 0)

  # Input validation
  expect_error(JAGS_marglik_priors(list(), prior_list = "x"), "'prior_list' must be a list.")
  expect_error(JAGS_marglik_priors(list(), prior_list = prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")
  expect_error(JAGS_marglik_priors(list(), prior_list = list(x = 1)), "'prior_list' must be a list of priors.")

})

test_that("JAGS_marglik_parameters input validation and edge cases work", {

  # Test: empty prior_list returns empty list
  result <- JAGS_marglik_parameters(list(), prior_list = list())
  expect_equal(result, list())

  # Test: prior_list must be a list
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = "not_a_list"),
    "'prior_list' must be a list."
  )

  # Test: prior_list must be a list of priors (single prior passed)
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = prior("normal", list(0, 1))),
    "'prior_list' must be a list of priors."
  )

  # Test: prior_list must be a list of priors (non-prior elements)
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = list(x = 1)),
    "'prior_list' must be a list of priors."
  )

})

test_that(".fit_to_posterior handles different input types", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)

  # mcmc.list (rjags::coda.samples)
  samples_mcmc_list <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 100, quiet = TRUE, progress.bar = "none")
  marglik <- JAGS_bridgesampling(samples_mcmc_list, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik, "BayesTools_marglik")

  # mcmc (coda::as.mcmc)
  samples_mcmc <- coda::as.mcmc(samples_mcmc_list[[1]])
  marglik_mcmc <- JAGS_bridgesampling(samples_mcmc, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_mcmc, "BayesTools_marglik")

  # Error for unsupported input
  expect_error(JAGS_bridgesampling("bad_input", prior_list = prior_list, data = list(), log_posterior = log_posterior), "not implemented")

})

test_that(".fit_to_posterior handles jags.samples output", {

  skip_if_not_installed("rjags")

  # Scalar parameter
  prior_list <- list(mu = prior("normal", list(0, 1)), sigma = prior("gamma", list(1, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples_jags <- rjags::jags.samples(model = model, variable.names = monitor, n.iter = 100, progress.bar = "none")
  marglik_jags <- JAGS_bridgesampling(samples_jags, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_jags, "BayesTools_marglik")

})

test_that(".fit_to_posterior handles vector parameters in jags.samples", {

  skip_if_not_installed("rjags")

  # Vector parameter (K > 1)
  prior_list <- list(p = prior("mnormal", list(mean = 0, sd = 1, K = 3)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples_jags <- rjags::jags.samples(model = model, variable.names = monitor, n.iter = 100, progress.bar = "none")
  marglik_jags <- JAGS_bridgesampling(samples_jags, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_jags, "BayesTools_marglik")

})

test_that("JAGS_bridgesampling handles runjags output", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  log_posterior <- STANDARD_LOG_POSTERIOR
  old_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.runjags = old_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.runjags = TRUE)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 2,
    adapt = 100,
    burnin = 100,
    sample = 500,
    silent.jags = TRUE,
    modules = "glm"
  ))

  marglik <- JAGS_bridgesampling(fit, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.1)

})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-nonlocal-fit.R
# ============================================================================ #

skip_if_not_test_profile("fit")

expect_nonlocal_prior_only_samples <- function(prior, samples, tolerance = .08) {
  expect_true(all(is.finite(samples)))
  expect_true(all(samples >= prior$truncation[["lower"]]))
  expect_true(all(samples <= prior$truncation[["upper"]]))

  probs <- c(.1, .25, .5, .75, .9)
  quantiles <- quant(prior, probs)
  sample_cdf <- vapply(quantiles, function(x) mean(samples <= x), numeric(1))
  expect_equal(sample_cdf, probs, tolerance = tolerance)

  prior_mean <- mean(prior)
  prior_sd <- sd(prior)
  if (is.finite(prior_mean) && is.finite(prior_sd)) {
    expect_equal(mean(samples), prior_mean, tolerance = tolerance)
    expect_equal(stats::sd(samples), prior_sd, tolerance = tolerance)
  }
}

test_that("BayesTools JAGS module initializes truncated nonlocal priors", {
  skip_if_not_installed("rjags")
  skip_on_cran()

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    moment = prior(
      "moment",
      list(mode = .5),
      truncation = list(lower = -.1, upper = .1)
    ),
    invmoment = prior(
      "invmoment",
      list(mode = .5, df = 6),
      truncation = list(lower = -.1, upper = .1)
    )
  )

  for(prior_i in priors){
    syntax <- JAGS_add_priors("model{}", list(theta = prior_i))
    expect_silent(local({
      con <- textConnection(syntax)
      on.exit(close(con), add = TRUE)
      rjags::jags.model(
        file     = con,
        data     = list(),
        n.chains = 1,
        n.adapt  = 0,
        quiet    = TRUE
      )
    }))
  }
})

test_that("BayesTools JAGS module samples nonlocal priors", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  skip_if_not_installed("bridgesampling")
  skip_on_cran()

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    moment = prior("moment", list(mode = .5)),
    invmoment = prior("invmoment", list(mode = .5, df = 6)),
    moment_truncated = prior("moment", list(mode = .5), truncation = list(lower = -Inf, upper = 0)),
    invmoment_truncated = prior("invmoment", list(mode = .5, df = 6), truncation = list(lower = -Inf, upper = 0))
  )

  for(prior_name in names(priors)){
    prior_list <- list(theta = priors[[prior_name]])
    fit <- suppressWarnings(JAGS_fit(
      model_syntax = "model{}",
      data = NULL,
      prior_list = prior_list,
      chains = 2,
      adapt = 250,
      burnin = 250,
      sample = 4000,
      silent = TRUE,
      seed = 1
    ))

    expect_s3_class(fit, "BayesTools_fit")
    expect_true("BayesTools" %in% attr(fit, "jags_modules"))
    samples <- as.matrix(fit$mcmc)
    expect_true("theta" %in% colnames(samples))
    expect_nonlocal_prior_only_samples(priors[[prior_name]], samples[, "theta"])

    marglik <- JAGS_bridgesampling(
      fit = fit,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = prior_list,
      maxiter = 2000
    )
    expect_s3_class(marglik, "BayesTools_marglik")
    expect_equal(marglik$logml, 0, tolerance = .08)
  }
})

test_that("fully structural fits retain deterministic draw geometry", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  fit <- JAGS_fit(
    model_syntax = "model{}",
    prior_list = list(theta = prior("point", list(0))),
    chains = 2,
    adapt = 100,
    burnin = 100,
    sample = 100,
    silent = TRUE,
    seed = 1
  )

  registry <- JAGS_parameter_registry(fit)
  expect_identical(
    registry$role[registry$name == .bt_backend_anchor_name],
    "backend_anchor"
  )
  expect_true(registry$internal[registry$name == .bt_backend_anchor_name])
  expect_identical(registry$monitor_status[registry$name == "theta"], "structural")
  expect_identical(registry$fixed_value[registry$name == "theta"], 0)

  catalog <- parameter_catalog(fit)
  expect_false(.bt_backend_anchor_name %in% catalog$quantities$canonical_name)
  theta <- catalog$quantities[catalog$quantities$canonical_name == "theta", ]
  expect_identical(theta$status, "structural")
  expect_identical(theta$fixed_value, 0)
  expect_identical(JAGS_fit_contract(fit)$parameter_catalog_version, 1L)

  geometry <- JAGS_draw_geometry(fit)
  expect_identical(geometry$chains$iterations, c(100L, 100L))
  expect_identical(geometry$total_draws, 200L)

  draws <- JAGS_materialize_draws(fit)
  expect_identical(colnames(draws[[1L]]), "theta")
  expect_identical(as.numeric(draws[[1L]][, "theta"]), rep(0, 100))
  expect_false(.bt_backend_anchor_name %in% colnames(draws[[1L]]))

  empty_draws <- JAGS_materialize_draws(fit, character())
  expect_identical(dim(empty_draws[[1L]]), c(100L, 0L))
  expect_identical(attr(empty_draws[[1L]], "mcpar"), c(201, 300, 1))

  extended <- JAGS_extend(
    fit,
    autofit_control = list(
      max_Rhat = NULL,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL,
      max_time = list(time = 30, unit = "secs"),
      sample_extend = 100,
      restarts = 1,
      max_extend = 1
    ),
    silent = TRUE
  )
  extended_geometry <- JAGS_draw_geometry(extended)
  expect_identical(extended_geometry$chains$iterations, c(200L, 200L))
  expect_identical(extended_geometry$total_draws, 400L)
  expect_identical(extended_geometry$chains$end, c(400L, 400L))
  expect_identical(parameter_catalog(extended), catalog)
})

# ============================================================================ #
# SAVE MODEL REGISTRY
# ============================================================================ #
# Convert the model registry list to a data frame for easy inspection and querying
test_that("Model registry is created and saved", {

  skip_on_cran()

  # Combine all registry entries into a single data frame
  model_registry_df <- do.call(rbind, model_registry)
  rownames(model_registry_df) <- NULL

  # Save the registry alongside the fitted models
  registry_file <- file.path(test_files_dir, "model_registry.RDS")
  saveRDS(model_registry_df, registry_file)

  # Verify registry was created
  expect_true(file.exists(registry_file))
  expect_s3_class(model_registry_df, "data.frame")
  expect_true(nrow(model_registry_df) > 0)

  expected_fit_files <- file.path(temp_fits_dir, paste0(model_registry_df$model_name, ".RDS"))
  expect_true(all(file.exists(expected_fit_files)))

  marglik_names <- model_registry_df$model_name[model_registry_df$has_marglik]
  expected_marglik_files <- file.path(temp_marglik_dir, paste0(marglik_names, ".RDS"))
  expect_true(all(file.exists(expected_marglik_files)))

  mark_refit_cache_complete(
    "model-fit",
    required_fits = model_registry_df$model_name,
    required_margliks = marglik_names,
    registry_file = registry_file
  )
})
