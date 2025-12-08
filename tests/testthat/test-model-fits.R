context("Model fits for reuse across tests")

# This file contains all model fitting procedures used across the test suite.
# Fitted models are saved to a temporary directory for reuse in other tests.
# This reduces redundant MCMC sampling and speeds up the overall test suite.

skip_on_cran()

# Setup directory for saving fitted models
temp_fits_dir <- file.path(tempdir(), "BayesTools_test_fits")
dir.create(temp_fits_dir, showWarnings = FALSE, recursive = TRUE)
# Set environment variable so other test files can locate pre-fitted models
Sys.setenv(BAYESTOOLS_TEST_FITS_DIR = temp_fits_dir)

# Helper function to save fitted models
save_fit <- function(fit, name) {
  saveRDS(fit, file = file.path(temp_fits_dir, paste0(name, ".RDS")))
}

# ==============================================================================
# SECTION 1: SIMPLE PRIOR DISTRIBUTIONS
# ==============================================================================
test_that("Simple prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  
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
  save_fit(fit_simple_normal, "fit_simple_normal")
  
  # Model 2: Various prior distributions
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
  model_syntax_simple <- JAGS_add_priors(model_syntax_simple, priors_various)
  monitor <- JAGS_to_monitor(priors_various)
  inits <- JAGS_get_inits(priors_various, chains = 2, seed = 1)
  
  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax_simple), 
                              inits = inits, n.chains = 2, quiet = TRUE)
  fit_simple_various <- rjags::coda.samples(model = model, variable.names = monitor, 
                                             n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_simple_various, "fit_simple_various")
  
  # Model 3: PET and PEESE priors
  priors_pub_bias <- list(
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
  )
  
  model_syntax_pb <- JAGS_add_priors("model{}", priors_pub_bias)
  monitor_pb <- JAGS_to_monitor(priors_pub_bias)
  inits_pb <- JAGS_get_inits(priors_pub_bias, chains = 2, seed = 1)
  
  set.seed(1)
  model_pb <- rjags::jags.model(file = textConnection(model_syntax_pb), 
                                 inits = inits_pb, n.chains = 2, quiet = TRUE)
  fit_simple_pub_bias <- rjags::coda.samples(model = model_pb, variable.names = monitor_pb, 
                                              n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_simple_pub_bias, "fit_simple_pub_bias")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_simple_normal.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_simple_various.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_simple_pub_bias.RDS")))
})


# ==============================================================================
# SECTION 2: VECTOR PRIOR DISTRIBUTIONS
# ==============================================================================
test_that("Vector prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  
  # Multivariate normal
  priors_mnormal <- list(
    p1 = prior("mnormal", list(mean = 0, sd = 1, K = 3))
  )
  
  model_syntax_vec <- JAGS_add_priors("model{}", priors_mnormal)
  monitor_vec <- JAGS_to_monitor(priors_mnormal)
  inits_vec <- JAGS_get_inits(priors_mnormal, chains = 2, seed = 1)
  
  set.seed(1)
  model_vec <- rjags::jags.model(file = textConnection(model_syntax_vec), 
                                  inits = inits_vec, n.chains = 2, quiet = TRUE)
  fit_vector_mnormal <- rjags::coda.samples(model = model_vec, variable.names = monitor_vec, 
                                             n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_vector_mnormal, "fit_vector_mnormal")
  
  # Multivariate cauchy
  priors_mcauchy <- list(
    p1 = prior("mcauchy", list(location = 0, scale = 1.5, K = 2))
  )
  
  model_syntax_mc <- JAGS_add_priors("model{}", priors_mcauchy)
  monitor_mc <- JAGS_to_monitor(priors_mcauchy)
  inits_mc <- JAGS_get_inits(priors_mcauchy, chains = 2, seed = 1)
  
  set.seed(1)
  model_mc <- rjags::jags.model(file = textConnection(model_syntax_mc), 
                                 inits = inits_mc, n.chains = 2, quiet = TRUE)
  fit_vector_mcauchy <- rjags::coda.samples(model = model_mc, variable.names = monitor_mc, 
                                             n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_vector_mcauchy, "fit_vector_mcauchy")
  
  # Multivariate t
  priors_mt <- list(
    p1 = prior("mt", list(location = 2, scale = 0.5, df = 5, K = 2))
  )
  
  model_syntax_mt <- JAGS_add_priors("model{}", priors_mt)
  monitor_mt <- JAGS_to_monitor(priors_mt)
  inits_mt <- JAGS_get_inits(priors_mt, chains = 2, seed = 1)
  
  set.seed(1)
  model_mt <- rjags::jags.model(file = textConnection(model_syntax_mt), 
                                 inits = inits_mt, n.chains = 2, quiet = TRUE)
  fit_vector_mt <- rjags::coda.samples(model = model_mt, variable.names = monitor_mt, 
                                        n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_vector_mt, "fit_vector_mt")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_vector_mnormal.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_vector_mcauchy.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_vector_mt.RDS")))
})


# ==============================================================================
# SECTION 3: FACTOR PRIOR DISTRIBUTIONS
# ==============================================================================
test_that("Factor prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))  # multivariate sampling differences
  
  # Orthonormal contrast
  priors_orthonormal <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_orthonormal[[1]], "levels") <- 3
  
  model_syntax_orth <- JAGS_add_priors("model{}", priors_orthonormal)
  monitor_orth <- JAGS_to_monitor(priors_orthonormal)
  inits_orth <- JAGS_get_inits(priors_orthonormal, chains = 2, seed = 1)
  
  set.seed(1)
  model_orth <- rjags::jags.model(file = textConnection(model_syntax_orth), 
                                   inits = inits_orth, n.chains = 2, quiet = TRUE)
  fit_factor_orthonormal <- rjags::coda.samples(model = model_orth, variable.names = monitor_orth, 
                                                 n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_factor_orthonormal, "fit_factor_orthonormal")
  
  # Treatment contrast
  priors_treatment <- list(
    p1 = prior_factor("beta", list(alpha = 1, beta = 1), contrast = "treatment")
  )
  attr(priors_treatment[[1]], "levels") <- 2
  
  model_syntax_treat <- JAGS_add_priors("model{}", priors_treatment)
  monitor_treat <- JAGS_to_monitor(priors_treatment)
  inits_treat <- JAGS_get_inits(priors_treatment, chains = 2, seed = 1)
  
  set.seed(1)
  model_treat <- rjags::jags.model(file = textConnection(model_syntax_treat), 
                                    inits = inits_treat, n.chains = 2, quiet = TRUE)
  fit_factor_treatment <- rjags::coda.samples(model = model_treat, variable.names = monitor_treat, 
                                               n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_factor_treatment, "fit_factor_treatment")
  
  # Independent contrast
  priors_independent <- list(
    p1 = prior_factor("gamma", list(shape = 2, rate = 3), contrast = "independent")
  )
  attr(priors_independent[[1]], "levels") <- 3
  
  model_syntax_ind <- JAGS_add_priors("model{}", priors_independent)
  monitor_ind <- JAGS_to_monitor(priors_independent)
  inits_ind <- JAGS_get_inits(priors_independent, chains = 2, seed = 1)
  
  set.seed(1)
  model_ind <- rjags::jags.model(file = textConnection(model_syntax_ind), 
                                  inits = inits_ind, n.chains = 2, quiet = TRUE)
  fit_factor_independent <- rjags::coda.samples(model = model_ind, variable.names = monitor_ind, 
                                                 n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_factor_independent, "fit_factor_independent")
  
  # Meandif contrast
  priors_meandif <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 0.5), contrast = "meandif")
  )
  attr(priors_meandif[[1]], "levels") <- 3
  
  model_syntax_md <- JAGS_add_priors("model{}", priors_meandif)
  monitor_md <- JAGS_to_monitor(priors_meandif)
  inits_md <- JAGS_get_inits(priors_meandif, chains = 2, seed = 1)
  
  set.seed(1)
  model_md <- rjags::jags.model(file = textConnection(model_syntax_md), 
                                 inits = inits_md, n.chains = 2, quiet = TRUE)
  fit_factor_meandif <- rjags::coda.samples(model = model_md, variable.names = monitor_md, 
                                             n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_factor_meandif, "fit_factor_meandif")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_factor_treatment.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_factor_independent.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_factor_meandif.RDS")))
})


# ==============================================================================
# SECTION 4: WEIGHTFUNCTION PRIORS
# ==============================================================================
test_that("Weightfunction prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  
  # One-sided weightfunction (2 intervals)
  priors_wf_onesided2 <- list(
    prior_weightfunction("one.sided", list(c(.05), c(1, 1)))
  )
  
  model_syntax_wf1 <- JAGS_add_priors("model{}", priors_wf_onesided2)
  monitor_wf1 <- JAGS_to_monitor(priors_wf_onesided2)
  inits_wf1 <- JAGS_get_inits(priors_wf_onesided2, chains = 2, seed = 1)
  
  set.seed(1)
  model_wf1 <- rjags::jags.model(file = textConnection(model_syntax_wf1), 
                                  inits = inits_wf1, n.chains = 2, quiet = TRUE)
  fit_weightfunction_onesided2 <- rjags::coda.samples(model = model_wf1, variable.names = monitor_wf1, 
                                                       n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_weightfunction_onesided2, "fit_weightfunction_onesided2")
  
  # One-sided weightfunction (3 intervals)
  priors_wf_onesided3 <- list(
    prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 2, 3)))
  )
  
  model_syntax_wf2 <- JAGS_add_priors("model{}", priors_wf_onesided3)
  monitor_wf2 <- JAGS_to_monitor(priors_wf_onesided3)
  inits_wf2 <- JAGS_get_inits(priors_wf_onesided3, chains = 2, seed = 1)
  
  set.seed(1)
  model_wf2 <- rjags::jags.model(file = textConnection(model_syntax_wf2), 
                                  inits = inits_wf2, n.chains = 2, quiet = TRUE)
  fit_weightfunction_onesided3 <- rjags::coda.samples(model = model_wf2, variable.names = monitor_wf2, 
                                                       n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_weightfunction_onesided3, "fit_weightfunction_onesided3")
  
  # Two-sided weightfunction
  priors_wf_twosided <- list(
    prior_weightfunction("two.sided", list(c(.05), c(1, 1)))
  )
  
  model_syntax_wf3 <- JAGS_add_priors("model{}", priors_wf_twosided)
  monitor_wf3 <- JAGS_to_monitor(priors_wf_twosided)
  inits_wf3 <- JAGS_get_inits(priors_wf_twosided, chains = 2, seed = 1)
  
  set.seed(1)
  model_wf3 <- rjags::jags.model(file = textConnection(model_syntax_wf3), 
                                  inits = inits_wf3, n.chains = 2, quiet = TRUE)
  fit_weightfunction_twosided <- rjags::coda.samples(model = model_wf3, variable.names = monitor_wf3, 
                                                      n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_weightfunction_twosided, "fit_weightfunction_twosided")
  
  # One-sided fixed weightfunction
  priors_wf_fixed <- list(
    prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)))
  )
  
  model_syntax_wf4 <- JAGS_add_priors("model{}", priors_wf_fixed)
  monitor_wf4 <- JAGS_to_monitor(priors_wf_fixed)
  inits_wf4 <- JAGS_get_inits(priors_wf_fixed, chains = 2, seed = 1)
  
  set.seed(1)
  model_wf4 <- rjags::jags.model(file = textConnection(model_syntax_wf4), 
                                  inits = inits_wf4, n.chains = 2, quiet = TRUE)
  fit_weightfunction_fixed <- rjags::coda.samples(model = model_wf4, variable.names = monitor_wf4, 
                                                   n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_weightfunction_fixed, "fit_weightfunction_fixed")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_weightfunction_onesided2.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_weightfunction_onesided3.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_weightfunction_twosided.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_weightfunction_fixed.RDS")))
})


# ==============================================================================
# SECTION 5: SPIKE-AND-SLAB PRIORS
# ==============================================================================
test_that("Spike-and-slab prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  
  # Simple spike-and-slab
  priors_spike_slab_simple <- list(
    "mu" = prior_spike_and_slab(prior("normal", list(0, 1)), 
                                 prior_inclusion = prior("beta", list(1,1)))
  )
  
  model_syntax_ss1 <- JAGS_add_priors("model{}", priors_spike_slab_simple)
  monitor_ss1 <- JAGS_to_monitor(priors_spike_slab_simple)
  inits_ss1 <- JAGS_get_inits(priors_spike_slab_simple, chains = 2, seed = 1)
  
  set.seed(1)
  model_ss1 <- rjags::jags.model(file = textConnection(model_syntax_ss1), 
                                  inits = inits_ss1, n.chains = 2, quiet = TRUE)
  fit_spike_slab_simple <- rjags::coda.samples(model = model_ss1, variable.names = monitor_ss1, 
                                                n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_spike_slab_simple, "fit_spike_slab_simple")
  
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
  
  model_syntax_ss2 <- JAGS_add_priors("model{}", priors_spike_slab_factor)
  monitor_ss2 <- JAGS_to_monitor(priors_spike_slab_factor)
  inits_ss2 <- JAGS_get_inits(priors_spike_slab_factor, chains = 2, seed = 1)
  
  set.seed(1)
  model_ss2 <- rjags::jags.model(file = textConnection(model_syntax_ss2), 
                                  inits = inits_ss2, n.chains = 2, quiet = TRUE)
  fit_spike_slab_factor <- rjags::coda.samples(model = model_ss2, variable.names = monitor_ss2, 
                                                n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_spike_slab_factor, "fit_spike_slab_factor")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_spike_slab_simple.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_spike_slab_factor.RDS")))
})


# ==============================================================================
# SECTION 6: MIXTURE PRIORS
# ==============================================================================
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
  
  model_syntax_mix1 <- JAGS_add_priors("model{}", priors_mixture_simple)
  monitor_mix1 <- JAGS_to_monitor(priors_mixture_simple)
  inits_mix1 <- JAGS_get_inits(priors_mixture_simple, chains = 2, seed = 1)
  
  set.seed(1)
  model_mix1 <- rjags::jags.model(file = textConnection(model_syntax_mix1), 
                                   inits = inits_mix1, n.chains = 2, quiet = TRUE)
  fit_mixture_simple <- rjags::coda.samples(model = model_mix1, variable.names = monitor_mix1, 
                                             n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_mixture_simple, "fit_mixture_simple")
  
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
  
  model_syntax_mix2 <- JAGS_add_priors("model{}", priors_mixture_components)
  monitor_mix2 <- JAGS_to_monitor(priors_mixture_components)
  inits_mix2 <- JAGS_get_inits(priors_mixture_components, chains = 2, seed = 1)
  
  set.seed(1)
  model_mix2 <- rjags::jags.model(file = textConnection(model_syntax_mix2), 
                                   inits = inits_mix2, n.chains = 2, quiet = TRUE)
  fit_mixture_components <- rjags::coda.samples(model = model_mix2, variable.names = monitor_mix2, 
                                                 n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_mixture_components, "fit_mixture_components")
  
  # Mixture with spike
  priors_mixture_spike <- list(
    "gamma" = prior_mixture(
      list(
        prior("spike", list(2)),
        prior("normal", list(-3, 1))
      )
    )
  )
  
  model_syntax_mix3 <- JAGS_add_priors("model{}", priors_mixture_spike)
  monitor_mix3 <- JAGS_to_monitor(priors_mixture_spike)
  inits_mix3 <- JAGS_get_inits(priors_mixture_spike, chains = 2, seed = 1)
  
  set.seed(1)
  model_mix3 <- rjags::jags.model(file = textConnection(model_syntax_mix3), 
                                   inits = inits_mix3, n.chains = 2, quiet = TRUE)
  fit_mixture_spike <- rjags::coda.samples(model = model_mix3, variable.names = monitor_mix3, 
                                            n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_mixture_spike, "fit_mixture_spike")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_mixture_simple.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_mixture_components.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_mixture_spike.RDS")))
})


# ==============================================================================
# SECTION 7: FORMULA-BASED MODELS (SIMPLE REGRESSION)
# ==============================================================================
test_that("Simple formula-based regression models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))
  
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
  save_fit(fit_formula_simple, "fit_formula_simple")
  
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
  save_fit(fit_formula_treatment, "fit_formula_treatment")
  
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
  save_fit(fit_formula_orthonormal, "fit_formula_orthonormal")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_simple.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_treatment.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_orthonormal.RDS")))
})


# ==============================================================================
# SECTION 8: FORMULA-BASED MODELS (INTERACTIONS)
# ==============================================================================
test_that("Formula-based interaction models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))
  
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
  save_fit(fit_formula_interaction_cont, "fit_formula_interaction_cont")
  
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
  save_fit(fit_formula_interaction_mix, "fit_formula_interaction_mix")
  
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
  save_fit(fit_formula_interaction_fac, "fit_formula_interaction_fac")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_interaction_mix.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_interaction_fac.RDS")))
})


# ==============================================================================
# SECTION 9: FORMULA-BASED MODELS (MULTIPLE FORMULAS)
# ==============================================================================
test_that("Multi-formula models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))
  
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
    chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 1)
  save_fit(fit_formula_multi, "fit_formula_multi")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_formula_multi.RDS")))
})


# ==============================================================================
# SECTION 10: RANDOM EFFECTS MODELS
# ==============================================================================
test_that("Random effects models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))
  
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
  save_fit(fit_random_intercept, "fit_random_intercept")
  
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
  save_fit(fit_random_slope, "fit_random_slope")
  
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
  save_fit(fit_random_factor_slope, "fit_random_factor_slope")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_random_intercept.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_random_slope.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_random_factor_slope.RDS")))
})


# ==============================================================================
# SECTION 11: SPIKE FACTOR PRIORS
# ==============================================================================
test_that("Spike factor prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))
  
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
  save_fit(fit_spike_factors, "fit_spike_factors")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_spike_factors.RDS")))
})


# ==============================================================================
# SECTION 12: JOINT MODELS (FORMULA + SPIKE-AND-SLAB + MIXTURE)
# ==============================================================================
test_that("Joint complex models fit correctly", {
  
  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris"))
  
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
  save_fit(fit_joint_complex, "fit_joint_complex")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_joint_complex.RDS")))
})


# ==============================================================================
# SECTION 13: EXPRESSION PRIORS
# ==============================================================================
test_that("Expression prior models fit correctly", {
  
  skip_if_not_installed("rjags")
  
  # Simple prior with expression
  priors_expr_simple <- list(
    x        = prior("normal",   list(0, expression(x_sigma))),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )
  
  model_syntax_expr1 <- JAGS_add_priors("model{}", priors_expr_simple)
  monitor_expr1 <- JAGS_to_monitor(priors_expr_simple)
  inits_expr1 <- JAGS_get_inits(priors_expr_simple, chains = 2, seed = 1)
  
  set.seed(1)
  model_expr1 <- rjags::jags.model(file = textConnection(model_syntax_expr1), 
                                    inits = inits_expr1, n.chains = 2, quiet = TRUE)
  fit_expression_simple <- rjags::coda.samples(model = model_expr1, variable.names = monitor_expr1, 
                                                n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_expression_simple, "fit_expression_simple")
  
  # Spike-and-slab with expression
  priors_expr_ss <- list(
    x        = prior_spike_and_slab(
      prior("normal", list(0, expression(x_sigma)))
    ),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )
  
  model_syntax_expr2 <- JAGS_add_priors("model{}", priors_expr_ss)
  monitor_expr2 <- JAGS_to_monitor(priors_expr_ss)
  inits_expr2 <- JAGS_get_inits(priors_expr_ss, chains = 2, seed = 1)
  
  set.seed(1)
  model_expr2 <- rjags::jags.model(file = textConnection(model_syntax_expr2), 
                                    inits = inits_expr2, n.chains = 2, quiet = TRUE)
  fit_expression_spike_slab <- rjags::coda.samples(model = model_expr2, variable.names = monitor_expr2, 
                                                    n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_expression_spike_slab, "fit_expression_spike_slab")
  
  # Mixture with expression
  priors_expr_mix <- list(
    x        = prior_mixture(list(
      prior("normal", list(0, expression(x_sigma))),
      prior("cauchy", list(0, 1))
    ), is_null = c(T, F)),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )
  
  model_syntax_expr3 <- JAGS_add_priors("model{}", priors_expr_mix)
  monitor_expr3 <- JAGS_to_monitor(priors_expr_mix)
  inits_expr3 <- JAGS_get_inits(priors_expr_mix, chains = 2, seed = 1)
  
  set.seed(1)
  model_expr3 <- rjags::jags.model(file = textConnection(model_syntax_expr3), 
                                    inits = inits_expr3, n.chains = 2, quiet = TRUE)
  fit_expression_mixture <- rjags::coda.samples(model = model_expr3, variable.names = monitor_expr3, 
                                                 n.iter = 1000, quiet = TRUE, progress.bar = "none")
  save_fit(fit_expression_mixture, "fit_expression_mixture")
  
  expect_true(file.exists(file.path(temp_fits_dir, "fit_expression_simple.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_expression_spike_slab.RDS")))
  expect_true(file.exists(file.path(temp_fits_dir, "fit_expression_mixture.RDS")))
})
