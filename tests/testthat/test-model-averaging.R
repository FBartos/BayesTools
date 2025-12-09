context("Model-averaging functions")

# ==============================================================================
# SECTION 1: BASIC MODEL-AVERAGING FUNCTIONS (NO JAGS FITS)
# ==============================================================================
test_that("Model-averaging functions work", {

  expect_equal(compute_inference(c(1,1), c(1, 1))$prior_probs, c(0.5, 0.5))
  expect_equal(compute_inference(c(1,1), c(1, 1))$post_probs,  c(0.5, 0.5))
  expect_equal(compute_inference(c(1,1), c(1, 1))$BF,          Inf)
  expect_equal(attr(compute_inference(c(1,1), c(1, 1)), "is_null"), c(FALSE, FALSE))

  expect_equal(compute_inference(c(1,4),   c(1, 1))$prior_probs,    c(0.2, 0.8))
  expect_equal(compute_inference(c(1,1,3), c(1, 1, 1))$prior_probs, c(0.2, 0.2, 0.6))
  expect_equal(compute_inference(c(1,1,4), c(1, 1, 1), c(F, T, F), conditional = TRUE)$prior_probs, c(0.2, 0, 0.8))

  expect_equal(compute_inference(c(1,4),   c(1, 1))$post_probs,    c(0.2, 0.8))
  expect_equal(compute_inference(c(1,1,3), c(1, 1, 1))$post_probs, c(0.2, 0.2, 0.6))
  expect_equal(compute_inference(c(1,1,4), c(1, 1, 1), c(F, T, F), conditional = TRUE)$post_probs, c(0.2, 0, 0.8))
  expect_equal(attr(compute_inference(c(1,1,4), c(1, 1, 1), c(2)), "is_null"), c(F, T, F))

  # automatically tests inclusion_bf as well
  expect_equal(compute_inference(c(1,1), c(1, 1), 1)$BF, 1)
  expect_equal(compute_inference(c(1,1), c(1, 2), c(F, T))$BF, exp(1-2))
  expect_equal(compute_inference(c(1,1,1), c(1, 1, 1), c(F, T, F))$BF, 1)
  expect_equal(compute_inference(c(1,1,1), c(1, 2, 1), c(F, T, F))$BF, exp(1-2))

  # and check BF formatting
  expect_equivalent(format_BF(c(0, 1, 2, Inf)), c(0, 1, 2, Inf))
  expect_equivalent(format_BF(c(0, 1, 2, Inf), BF01 = TRUE), 1/c(0, 1, 2, Inf))
  expect_equivalent(format_BF(c(0, 1, 2, Inf), logBF = TRUE), log(c(0, 1, 2, Inf)))
  expect_equivalent(format_BF(c(0, 1, 2, Inf), BF01 = TRUE, logBF = TRUE), log(1/c(0, 1, 2, Inf)))
  expect_equal(attr(format_BF(1), "name"), "BF")
  expect_equal(attr(format_BF(1, logBF = TRUE), "name"), "log(BF)")
  expect_equal(attr(format_BF(1, BF01 = TRUE, logBF = TRUE), "name"), "log(1/BF)")

  # additional BF checks
  expect_equal(inclusion_BF(prior_probs = c(.5, .5), post_probs = c(.5, .5), is_null = c(T, F)), 1)
  expect_equal(inclusion_BF(prior_probs = c(.5, .5), post_probs = c(.75, .25), is_null = c(T, F)), 1/3)
  expect_equal(inclusion_BF(prior_probs = c(.25, .25, .25, .25), post_probs = c(.75, 0, .25, 0), is_null = c(T, T, F, F)), 1/3)
  expect_equal(inclusion_BF(prior_probs = c(.25, .25, .25, .25), post_probs = c(.65, .10, .20, 0.05), is_null = c(T, T, F, F)), 1/3)
  expect_equal(inclusion_BF(prior_probs = c(1, 0), post_probs = c(1, 0), is_null = c(T, F)), 0)
  expect_equal(inclusion_BF(prior_probs = c(1, 0), post_probs = c(1, 0), is_null = c(F, T)), Inf)

  # test the marglik versions of BF
  temp_prior_probs <- 1:6/sum(1:6)
  temp_margliks    <- -2:3
  temp_post_probs  <- bridgesampling::post_prob(temp_margliks, prior_prob = temp_prior_probs)
  expect_equal(
    inclusion_BF(prior_probs = temp_prior_probs, post_probs = temp_post_probs, is_null = rep(c(T, F), 3)),
    inclusion_BF(prior_probs = temp_prior_probs, margliks = temp_margliks, is_null = rep(c(T, F), 3))
  )

  # check for over/underflow
  temp_prior_probs <- 1:6/sum(1:6)
  temp_margliks    <- c(-2:2, 100)
  temp_post_probs  <- bridgesampling::post_prob(temp_margliks, prior_prob = temp_prior_probs)
  expect_true(is.infinite(inclusion_BF(prior_probs = temp_prior_probs, post_probs = temp_post_probs, is_null = rep(c(T, F), 3))))
  expect_false(is.infinite(inclusion_BF(prior_probs = temp_prior_probs, margliks = temp_margliks, is_null = rep(c(T, F), 3))))

  # additional omega mapping checks
  expect_equal(weightfunctions_mapping(prior_list = list(
    prior_none(),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1),       steps = c(0.05)),        prior_weights = 1/2),
    prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1),    steps = c(0.05, 0.10)),  prior_weights = 1/2)
  )), list(NULL, c(2, 1, 1), c(3, 2, 1)))
})


# ==============================================================================
# SECTION 2: JAGS MODEL-AVERAGING WITH PREFITTED MODELS
# ==============================================================================
# Skip on CRAN as these tests use pre-fitted models
skip_on_cran()

# Get the directory where prefitted models are stored
# First check environment variable, then fall back to standard temp directory
temp_fits_dir <- Sys.getenv("BAYESTOOLS_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  temp_fits_dir <- file.path(tempdir(), "BayesTools_test_fits")
}
if (!dir.exists(temp_fits_dir)) {
  skip("Pre-fitted models not available. Run test-00-model-fits.R first.")
}

test_that("JAGS model-averaging with simple priors", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Load pre-fitted models
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  fit_simple_various <- readRDS(file.path(temp_fits_dir, "fit_simple_various.RDS"))
  
  # Test data from fit_simple_normal
  set.seed(1)
  data <- list(
    x = rnorm(50, 0, .5),
    N = 50
  )
  
  # Define log posterior for bridgesampling
  log_posterior_normal <- function(parameters, data){
    sum(stats::dnorm(data$x, parameters[["m"]], parameters[["s"]], log = TRUE))
  }
  
  # Get priors for the models
  priors_normal <- list(
    m = prior("normal", list(0, 1)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  
  # Create an alternative model with spike on m
  priors_spike <- list(
    m = prior("spike", list(0)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"
  
  fit_spike <- JAGS_fit(model_syntax, data, priors_spike, 
                        chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  
  # Get marginal likelihoods
  marglik_normal <- JAGS_bridgesampling(fit_simple_normal, 
                                         log_posterior = log_posterior_normal, 
                                         data = data, prior_list = priors_normal)
  marglik_spike <- JAGS_bridgesampling(fit_spike, 
                                        log_posterior = log_posterior_normal, 
                                        data = data, prior_list = priors_spike)
  
  # Create model list
  models <- list(
    list(fit = fit_spike, marglik = marglik_spike, prior_weights = 1),
    list(fit = fit_simple_normal, marglik = marglik_normal, prior_weights = 1)
  )
  
  # Test ensemble inference
  inference <- ensemble_inference(model_list = models, parameters = c("m", "s"), 
                                  is_null_list = list("m" = 1, "s" = 0), conditional = FALSE)
  inference_conditional <- ensemble_inference(model_list = models, parameters = c("m", "s"), 
                                               is_null_list = list("m" = 1, "s" = 0), conditional = TRUE)
  
  # Test mix posteriors
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "s"), 
                                      is_null_list = list("m" = 1, "s" = 0), seed = 1)
  mixed_posteriors_conditional <- mix_posteriors(model_list = models, parameters = c("m", "s"), 
                                                   is_null_list = list("m" = 1, "s" = 0), conditional = TRUE, seed = 1)
  
  # Checks
  expect_true(is.list(inference))
  expect_true(all(c("m", "s") %in% names(inference)))
  expect_true(is.numeric(inference$m$BF))
  expect_true(is.numeric(inference$s$BF))
  expect_equal(length(mixed_posteriors$m), length(mixed_posteriors$s))
  expect_true(mean(mixed_posteriors$m == 0) > 0) # Some spike samples
  
  # Visual check
  vdiffr::expect_doppelganger("model-averaging-simple-priors", function(){
    par(mfrow = c(2, 2))
    hist(mixed_posteriors$m, main = "model-averaged (m)")
    hist(mixed_posteriors_conditional$m, main = "conditional (m)")
    hist(mixed_posteriors$s, main = "model-averaged (s)")
    hist(mixed_posteriors_conditional$s, main = "conditional (s)")
  })
})

test_that("JAGS model-averaging with weightfunction priors", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Load pre-fitted weightfunction models
  fit_wf_none <- readRDS(file.path(temp_fits_dir, "fit_weightfunction_none.RDS"))
  fit_wf_onesided2 <- readRDS(file.path(temp_fits_dir, "fit_weightfunction_onesided2.RDS"))
  fit_wf_onesided3 <- readRDS(file.path(temp_fits_dir, "fit_weightfunction_onesided3.RDS"))
  fit_wf_twosided <- readRDS(file.path(temp_fits_dir, "fit_weightfunction_twosided.RDS"))
  
  # Create dummy data and log posterior (weightfunctions don't use data)
  data <- list()
  log_posterior <- function(parameters, data) { return(0) }
  
  # Define priors
  priors_none <- list(omega = prior_none())
  priors_onesided2 <- list(prior_weightfunction("one.sided", list(c(.05), c(1, 1))))
  priors_onesided3 <- list(prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 2, 3))))
  priors_twosided <- list(prior_weightfunction("two.sided", list(c(.05), c(1, 1))))
  
  # Get marginal likelihoods
  marglik_none <- JAGS_bridgesampling(fit_wf_none, log_posterior = log_posterior, 
                                       data = data, prior_list = priors_none)
  marglik_onesided2 <- JAGS_bridgesampling(fit_wf_onesided2, log_posterior = log_posterior, 
                                            data = data, prior_list = priors_onesided2)
  marglik_onesided3 <- JAGS_bridgesampling(fit_wf_onesided3, log_posterior = log_posterior, 
                                            data = data, prior_list = priors_onesided3)
  marglik_twosided <- JAGS_bridgesampling(fit_wf_twosided, log_posterior = log_posterior, 
                                           data = data, prior_list = priors_twosided)
  
  # Test coefficient mapping
  expect_equal(
    weightfunctions_mapping(list(priors_none$omega, priors_onesided2[[1]], priors_onesided3[[1]])), 
    list(NULL, c(2, 1, 1), c(3, 2, 1))
  )
  
  # Create model list
  models <- list(
    list(fit = fit_wf_none, marglik = marglik_none, prior_weights = 1),
    list(fit = fit_wf_onesided2, marglik = marglik_onesided2, prior_weights = 1),
    list(fit = fit_wf_onesided3, marglik = marglik_onesided3, prior_weights = 1),
    list(fit = fit_wf_twosided, marglik = marglik_twosided, prior_weights = 1)
  )
  
  # Test ensemble inference
  inference <- ensemble_inference(model_list = models, parameters = c("omega"), 
                                  is_null_list = list("omega" = 1), conditional = FALSE)
  
  # Test mix posteriors
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("omega"), 
                                      is_null_list = list("omega" = 1), seed = 1)
  mixed_posteriors_conditional <- mix_posteriors(model_list = models, parameters = c("omega"), 
                                                   is_null_list = list("omega" = 1), conditional = TRUE, seed = 1)
  
  # Checks
  expect_true(is.list(inference))
  expect_true("omega" %in% names(inference))
  expect_true(is.numeric(inference$omega$BF))
  expect_true(is.matrix(mixed_posteriors$omega))
  expect_true(all(mixed_posteriors$omega[1,] == 1)) # First row should be all 1s (no weightfunction)
  
  # Visual check
  vdiffr::expect_doppelganger("model-averaging-weightfunctions", function(){
    par(mfrow = c(2, 3))
    for(i in 1:min(3, ncol(mixed_posteriors$omega))) {
      hist(mixed_posteriors$omega[,i], main = paste("averaged omega", i), 
           xlab = colnames(mixed_posteriors$omega)[i])
    }
    for(i in 1:min(3, ncol(mixed_posteriors_conditional$omega))) {
      hist(mixed_posteriors_conditional$omega[,i], main = paste("conditional omega", i), 
           xlab = colnames(mixed_posteriors_conditional$omega)[i])
    }
  })
})

test_that("JAGS model-averaging with mixture priors", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Load pre-fitted mixture models
  fit_mixture_simple <- readRDS(file.path(temp_fits_dir, "fit_mixture_simple.RDS"))
  fit_mixture_components <- readRDS(file.path(temp_fits_dir, "fit_mixture_components.RDS"))
  fit_mixture_spike <- readRDS(file.path(temp_fits_dir, "fit_mixture_spike.RDS"))
  
  # Create dummy data and log posterior
  data <- list()
  log_posterior <- function(parameters, data) { return(0) }
  
  # Define priors matching the fitted models
  priors_simple <- list(
    "mu" = prior_mixture(
      list(
        prior("normal", list(0,  1), prior_weights = 1),
        prior("normal", list(-3, 1), prior_weights = 5),
        prior("gamma",  list(5, 10), prior_weights = 1)
      ),
      is_null = c(T, F, T)
    )
  )
  
  priors_components <- list(
    "beta" = prior_mixture(
      list(
        prior("normal", list(0,  1), prior_weights = 1),
        prior("normal", list(-3, 1), prior_weights = 5)
      ),
      components = c("b", "a")
    )
  )
  
  priors_spike <- list(
    "gamma" = prior_mixture(
      list(
        prior("spike", list(2)),
        prior("normal", list(-3, 1))
      )
    )
  )
  
  # Get marginal likelihoods
  marglik_simple <- JAGS_bridgesampling(fit_mixture_simple, log_posterior = log_posterior, 
                                         data = data, prior_list = priors_simple)
  marglik_components <- JAGS_bridgesampling(fit_mixture_components, log_posterior = log_posterior, 
                                             data = data, prior_list = priors_components)
  marglik_spike <- JAGS_bridgesampling(fit_mixture_spike, log_posterior = log_posterior, 
                                        data = data, prior_list = priors_spike)
  
  # Create model list
  models <- list(
    list(fit = fit_mixture_simple, marglik = marglik_simple, prior_weights = 1),
    list(fit = fit_mixture_components, marglik = marglik_components, prior_weights = 1),
    list(fit = fit_mixture_spike, marglik = marglik_spike, prior_weights = 1)
  )
  
  # Test models_inference
  models_with_inference <- models_inference(models)
  
  # Checks
  expect_true(length(models_with_inference) == 3)
  expect_true(all(sapply(models_with_inference, function(m) "inference" %in% names(m))))
  expect_true(all(sapply(models_with_inference, function(m) is.list(m$inference))))
})

test_that("JAGS model-averaging with spike-and-slab priors", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Load pre-fitted spike-and-slab models
  fit_spike_slab_simple <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_simple.RDS"))
  fit_spike_slab_factor <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_factor.RDS"))
  
  # Create dummy data and log posterior
  data <- list()
  log_posterior <- function(parameters, data) { return(0) }
  
  # Define priors
  priors_simple <- list(
    "mu" = prior_spike_and_slab(prior("normal", list(0, 1)),
                                 prior_inclusion = prior("beta", list(1,1)))
  )
  
  priors_factor <- list(
    "beta" = prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
                                   prior_inclusion = prior("beta", list(1,1)))
  )
  components <- attr(priors_factor$beta, "components")
  alternative_idx <- which(components == "alternative")
  attr(priors_factor$beta[[alternative_idx]], "levels") <- 3
  
  # Get marginal likelihoods
  marglik_simple <- JAGS_bridgesampling(fit_spike_slab_simple, log_posterior = log_posterior, 
                                         data = data, prior_list = priors_simple)
  marglik_factor <- JAGS_bridgesampling(fit_spike_slab_factor, log_posterior = log_posterior, 
                                         data = data, prior_list = priors_factor)
  
  # Create model list
  models <- list(
    list(fit = fit_spike_slab_simple, marglik = marglik_simple, prior_weights = 1),
    list(fit = fit_spike_slab_factor, marglik = marglik_factor, prior_weights = 1)
  )
  
  # Test models_inference
  models_with_inference <- models_inference(models)
  
  # Checks
  expect_true(length(models_with_inference) == 2)
  expect_true(all(sapply(models_with_inference, function(m) "inference" %in% names(m))))
})

test_that("JAGS model-averaging with formula models", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Load pre-fitted formula models
  fit_formula_simple <- readRDS(file.path(temp_fits_dir, "fit_formula_simple.RDS"))
  fit_formula_treatment <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  fit_formula_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_formula_orthonormal.RDS"))
  
  # Recreate data
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
  
  # Define log posterior
  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  
  # Define priors and formulas
  formula_list_simple <- list(mu = ~ x_cont1)
  formula_data_list_simple <- list(mu = data_formula)
  formula_prior_list_simple <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1))
    )
  )
  prior_list_simple <- list(sigma = prior("lognormal", list(0, 1)))
  
  formula_list_treatment <- list(mu = ~ x_cont1 + x_fac2t)
  formula_data_list_treatment <- list(mu = data_formula)
  formula_prior_list_treatment <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1)),
      "x_fac2t"   = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )
  
  formula_list_orthonormal <- list(mu = ~ x_cont1 + x_fac3o)
  formula_data_list_orthonormal <- list(mu = data_formula)
  formula_prior_list_orthonormal <- list(
    mu = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior("normal", list(0, 1)),
      "x_fac3o"   = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )
  
  # Get marginal likelihoods
  marglik_simple <- JAGS_bridgesampling(
    fit_formula_simple, log_posterior = log_posterior, data = data, 
    prior_list = prior_list_simple,
    formula_list = formula_list_simple, 
    formula_data_list = formula_data_list_simple, 
    formula_prior_list = formula_prior_list_simple)
  
  marglik_treatment <- JAGS_bridgesampling(
    fit_formula_treatment, log_posterior = log_posterior, data = data, 
    prior_list = prior_list_simple,
    formula_list = formula_list_treatment, 
    formula_data_list = formula_data_list_treatment, 
    formula_prior_list = formula_prior_list_treatment)
  
  marglik_orthonormal <- JAGS_bridgesampling(
    fit_formula_orthonormal, log_posterior = log_posterior, data = data, 
    prior_list = prior_list_simple,
    formula_list = formula_list_orthonormal, 
    formula_data_list = formula_data_list_orthonormal, 
    formula_prior_list = formula_prior_list_orthonormal)
  
  # Create model list
  models <- list(
    list(fit = fit_formula_simple, marglik = marglik_simple, prior_weights = 1),
    list(fit = fit_formula_treatment, marglik = marglik_treatment, prior_weights = 1),
    list(fit = fit_formula_orthonormal, marglik = marglik_orthonormal, prior_weights = 1)
  )
  
  # Test ensemble inference
  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3o"),
    is_null_list = list(
      "mu_x_cont1" = c(FALSE, FALSE, FALSE),
      "mu_x_fac2t" = c(TRUE,  FALSE, TRUE),
      "mu_x_fac3o" = c(TRUE,  TRUE,  FALSE)
    ),
    conditional = FALSE)
  
  # Test mix posteriors
  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3o"),
    is_null_list = list(
      "mu_x_cont1" = c(FALSE, FALSE, FALSE),
      "mu_x_fac2t" = c(TRUE,  FALSE, TRUE),
      "mu_x_fac3o" = c(TRUE,  TRUE,  FALSE)
    ),
    seed = 1, n_samples = 1000)
  
  # Checks
  expect_true(is.list(inference))
  expect_true(all(c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3o") %in% names(inference)))
  expect_true(is.numeric(inference$mu_x_cont1$BF))
  expect_true(is.numeric(inference$mu_x_fac2t$BF))
  expect_true(is.numeric(inference$mu_x_fac3o$BF))
  expect_equal(length(mixed_posteriors$mu_x_cont1), 1000)
  
  # Visual check
  vdiffr::expect_doppelganger("model-averaging-formulas", function(){
    par(mfrow = c(2, 2))
    hist(mixed_posteriors$mu_x_cont1, main = "mu_x_cont1")
    hist(mixed_posteriors$mu_x_fac2t, main = "mu_x_fac2t")
    if(is.matrix(mixed_posteriors$mu_x_fac3o)) {
      hist(mixed_posteriors$mu_x_fac3o[,1], main = "mu_x_fac3o[1]")
      hist(mixed_posteriors$mu_x_fac3o[,2], main = "mu_x_fac3o[2]")
    }
  })
})


# ==============================================================================
# SECTION 3: PRINT TESTS WITH TXT FILE COMPARISON
# ==============================================================================
test_that("Model-averaging print output matches expected format", {
  
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  
  # Create a list to store fits for printing
  fits <- list()
  
  # Load some pre-fitted models
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  fit_simple_various <- readRDS(file.path(temp_fits_dir, "fit_simple_various.RDS"))
  
  # Create test data and models for averaging
  set.seed(1)
  data <- list(
    x = rnorm(50, 0, .5),
    N = 50
  )
  
  log_posterior_normal <- function(parameters, data){
    sum(stats::dnorm(data$x, parameters[["m"]], parameters[["s"]], log = TRUE))
  }
  
  priors_normal <- list(
    m = prior("normal", list(0, 1)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  
  priors_spike <- list(
    m = prior("spike", list(0)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"
  
  fit_spike <- JAGS_fit(model_syntax, data, priors_spike, 
                        chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
  
  # Get marginal likelihoods
  marglik_normal <- JAGS_bridgesampling(fit_simple_normal, 
                                         log_posterior = log_posterior_normal, 
                                         data = data, prior_list = priors_normal)
  marglik_spike <- JAGS_bridgesampling(fit_spike, 
                                        log_posterior = log_posterior_normal, 
                                        data = data, prior_list = priors_spike)
  
  # Create model list
  models <- list(
    list(fit = fit_spike, marglik = marglik_spike, prior_weights = 1),
    list(fit = fit_simple_normal, marglik = marglik_normal, prior_weights = 1)
  )
  
  # Generate ensemble inference and mixed posteriors
  inference <- ensemble_inference(model_list = models, parameters = c("m", "s"), 
                                  is_null_list = list("m" = 1, "s" = 0), conditional = FALSE)
  
  # Store inference objects for print testing
  fits[[1]] <- inference
  
  # Compare printed output with saved files
  for(i in 1:length(fits)){
    saved_file <- file.path("../results/print", paste0(i, ".txt"))
    if(file.exists(saved_file)) {
      expect_equal(
        capture_output_lines(fits[[i]], print = TRUE, width = 150),
        read.table(file = saved_file, header = FALSE, blank.lines.skip = FALSE)[,1])
    } else {
      # If file doesn't exist yet, skip the test (will be generated by UPDATE_OUTPUT script)
      skip(paste("Print test file", i, "not yet generated. Run UPDATE_OUTPUT section."))
    }
  }
})


# ==============================================================================
# OUTPUT GENERATION SCRIPT (UPDATE_OUTPUT = FALSE)
# ==============================================================================
# This section generates the initial print test files
# Set UPDATE_OUTPUT = TRUE to regenerate the expected output files
UPDATE_OUTPUT <- FALSE

if(UPDATE_OUTPUT && !identical(Sys.getenv("CI"), "true")) {
  
  test_that("Generate print output files", {
    
    skip_if_not_installed("rjags")
    skip_if_not_installed("bridgesampling")
    skip_on_cran()
    skip_on_ci()
    
    # Ensure print directory exists
    print_dir <- file.path("tests/results/print")
    if(!dir.exists(print_dir)) {
      dir.create(print_dir, recursive = TRUE)
    }
    
    # Create a list to store fits for printing
    fits <- list()
    
    # Load some pre-fitted models
    fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
    fit_simple_various <- readRDS(file.path(temp_fits_dir, "fit_simple_various.RDS"))
    
    # Create test data and models for averaging
    set.seed(1)
    data <- list(
      x = rnorm(50, 0, .5),
      N = 50
    )
    
    log_posterior_normal <- function(parameters, data){
      sum(stats::dnorm(data$x, parameters[["m"]], parameters[["s"]], log = TRUE))
    }
    
    priors_normal <- list(
      m = prior("normal", list(0, 1)),
      s = prior("normal", list(0, 1), list(0, Inf))
    )
    
    priors_spike <- list(
      m = prior("spike", list(0)),
      s = prior("normal", list(0, 1), list(0, Inf))
    )
    
    model_syntax <-
      "model
      {
        for(i in 1:N){
          x[i] ~ dnorm(m, pow(s, -2))
        }
      }"
    
    fit_spike <- JAGS_fit(model_syntax, data, priors_spike, 
                          chains = 2, adapt = 100, burnin = 150, sample = 500, seed = 2)
    
    # Get marginal likelihoods
    marglik_normal <- JAGS_bridgesampling(fit_simple_normal, 
                                           log_posterior = log_posterior_normal, 
                                           data = data, prior_list = priors_normal)
    marglik_spike <- JAGS_bridgesampling(fit_spike, 
                                          log_posterior = log_posterior_normal, 
                                          data = data, prior_list = priors_spike)
    
    # Create model list
    models <- list(
      list(fit = fit_spike, marglik = marglik_spike, prior_weights = 1),
      list(fit = fit_simple_normal, marglik = marglik_normal, prior_weights = 1)
    )
    
    # Generate ensemble inference
    inference <- ensemble_inference(model_list = models, parameters = c("m", "s"), 
                                    is_null_list = list("m" = 1, "s" = 0), conditional = FALSE)
    
    # Store inference objects for print testing
    fits[[1]] <- inference
    
    # Generate print files
    for(i in seq_along(fits)){
      write.table(capture_output_lines(fits[[i]], print = TRUE, width = 150), 
                  file = file.path(print_dir, paste0(i, ".txt")), 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
    message("Print output files generated in ", print_dir)
  })
}
