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
skip_refit_if_cached(
  "model-fit",
  required_fits = c(
    "fit_wf_independent_gamma",
    "fit_wf_independent_log",
    "fit_bias_heterogeneous_wf",
    "fit_bias_petpeese_heterogeneous_wf",
    "fit_selection_kernel_summary"
  )
)

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
    prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1)))
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
    prior_weightfunction("one-sided", c(.05, 0.10), wf_cumulative(c(1, 2, 3)))
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
    prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)))
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
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5)))
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
  fit_bias_petpeese_heterogeneous_wf <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(bias = bias_petpeese_heterogeneous_wf),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 1200,
    seed       = 15
  ))
  result <- save_fit(fit_bias_petpeese_heterogeneous_wf, "fit_bias_petpeese_heterogeneous_wf",
                     pub_bias_priors = TRUE, weightfunction_priors = TRUE,
                     mixture_priors = TRUE,
                     note = "Full bias mixture with PET, PEESE, and heterogeneous weightfunctions")
  model_registry[["fit_bias_petpeese_heterogeneous_wf"]] <<- result$registry_entry
  fit_bias_petpeese_heterogeneous_wf <- result$fit

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
      "intercept"    = prior("normal", list(0, 5)),
      "intercept|id" = prior("normal", list(0, 1), list(0, 1))
    )
  )

  fit_random_intercept <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_int, formula_data_list = formula_data_list_re_int,
    formula_prior_list = formula_prior_list_re_int,
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
      "intercept"  = prior("normal", list(0, 5)),
      "x_cont1|id" = prior("normal", list(0, 1), list(0, 1))
    )
  )

  fit_random_slope <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_slope, formula_data_list = formula_data_list_re_slope,
    formula_prior_list = formula_prior_list_re_slope,
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
      "x_cont1"      = prior("normal", list(0, 1)),
      "intercept|id" = prior("normal", list(0, 1), list(0, 1)),
      "x_fac3|id"    = prior("normal", list(0, 1), list(0, 1))
    )
  )

  fit_random_factor_slope <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_fac, formula_data_list = formula_data_list_re_fac,
    formula_prior_list = formula_prior_list_re_fac,
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
      "x_fac3"       = prior_factor("mnormal", list(0, 1)),
      "intercept|id" = prior("normal", list(0, 1), list(0, 1)),
      "x_fac3|id"    = prior("normal", list(0, 1), list(0, 1))
    )
  )

  fit_random_factor_slope2 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_fac, formula_data_list = formula_data_list_re_fac,
    formula_prior_list = formula_prior_list_re_fac,
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
      "x_fac3"       = prior_factor("normal", list(0, 1), contrast = "independent"),
      "x_fac3|id"    = prior_spike_and_slab(prior("normal", list(0, 1), list(0, 1)))
    )
  )

  fit_random_factor_slope3 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list_re_fac, formula_data_list = formula_data_list_re_fac,
    formula_prior_list = formula_prior_list_re_fac,
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
  skip_on_os(c("mac", "linux", "solaris"))
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
})
