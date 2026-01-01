# ============================================================================ #
# TEST FILE: Posterior Density Extraction Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for posterior extraction helper functions including
#   .extract_posterior_samples and .remove_auxiliary_parameters.
#
# DEPENDENCIES:
#   - rjags, runjags, coda: For JAGS model handling
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"), skip_if_not_installed("runjags")
#   - Note: Creates mock objects, does not need pre-fitted models
#
# MODELS/FIXTURES:
#   - Creates mock runjags objects for testing
#
# TAGS: @evaluation, @JAGS, @posterior-extraction
# ============================================================================ #

# Tests for posterior extraction helper functions
test_that(".extract_posterior_samples extracts samples correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  
  # Load runjags to ensure S3 methods are registered
  library(runjags)

  # Create a proper runjags object structure for testing
  # The runjags package has an as.mcmc method that handles mcmc.list objects
  set.seed(123)
  mcmc1 <- coda::mcmc(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu")),
                      start = 1, end = 100, thin = 1)
  mcmc2 <- coda::mcmc(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu")),
                      start = 1, end = 100, thin = 1)
  mcmc_list <- coda::mcmc.list(mcmc1, mcmc2)
  
  # Create a minimal runjags object
  fit <- structure(
    list(mcmc = mcmc_list),
    class = c("runjags", "list")
  )

  # Test matrix extraction (as_list = FALSE)
  # This calls coda::as.mcmc on the runjags object which returns an mcmc object
  samples_matrix <- BayesTools:::.extract_posterior_samples(fit, as_list = FALSE)
  # mcmc objects inherit from matrix
  expect_true(inherits(samples_matrix, "mcmc"))
  expect_equal(ncol(samples_matrix), 1)
  expect_true("mu" %in% colnames(samples_matrix))
  expect_equal(nrow(samples_matrix), 200)  # 100 samples x 2 chains merged

  # Test list extraction (as_list = TRUE)
  samples_list <- BayesTools:::.extract_posterior_samples(fit, as_list = TRUE)
  expect_true(inherits(samples_list, "mcmc.list"))
  expect_equal(length(samples_list), 2) # 2 chains
})


test_that(".remove_auxiliary_parameters removes invgamma support", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with invgamma support parameter
  model_samples <- matrix(rnorm(100), ncol = 2)
  colnames(model_samples) <- c("sigma", "inv_sigma")

  prior_list <- list(
    sigma = prior("invgamma", list(1, 1))
  )

  result <- BayesTools:::.remove_auxiliary_parameters(model_samples, prior_list, NULL)

  expect_false("inv_sigma" %in% colnames(result$model_samples))
  expect_true("sigma" %in% colnames(result$model_samples))
})


test_that(".process_spike_and_slab handles conditional samples", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with spike and slab
  model_samples <- matrix(c(
    rnorm(50, 0, 1),  # mu values
    rep(1, 50)        # indicator (all in slab)
  ), ncol = 2)
  colnames(model_samples) <- c("mu", "mu_indicator")

  prior_list <- list(
    mu = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("spike", list(0.5))
    )
  )

  result <- BayesTools:::.process_spike_and_slab(
    model_samples, prior_list, "mu",
    conditional = TRUE, remove_inclusion = FALSE, warnings = NULL
  )

  expect_true("mu (inclusion)" %in% colnames(result$model_samples))
  expect_false("mu_indicator" %in% colnames(result$model_samples))
  expect_true(is.prior.simple(result$prior_list$mu))
})


test_that(".apply_parameter_transformations applies transformations", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples
  model_samples <- matrix(rnorm(100, 0, 1), ncol = 1)
  colnames(model_samples) <- "mu"

  prior_list <- list(
    mu = prior("normal", list(0, 1))
  )

  # Apply exp transformation
  transformations <- list(
    mu = list(fun = exp, arg = list())
  )

  result <- BayesTools:::.apply_parameter_transformations(
    model_samples, transformations, prior_list
  )

  expect_true(all(result[, "mu"] > 0))  # exp makes all values positive
  expect_equal(ncol(result), 1)
})


test_that(".rename_factor_levels renames treatment factors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with factor
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior_obj, "levels") <- 4  # 4 levels total (treatment has K-1 parameters for K levels)
  attr(prior_obj, "level_names") <- c("A", "B", "C", "D")  # Should be a vector, not a list
  
  prior_list <- list(group = prior_obj)

  result <- BayesTools:::.rename_factor_levels(model_samples, prior_list)

  expect_true("group[B]" %in% colnames(result))
  expect_true("group[C]" %in% colnames(result))
  expect_true("group[D]" %in% colnames(result))
  expect_false("group[1]" %in% colnames(result))
})


test_that(".transform_factor_contrasts transforms orthonormal to differences", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with orthonormal contrasts
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  attr(prior_obj, "levels") <- 4  # 4 levels total (orthonormal has K-1 parameters for K levels)
  attr(prior_obj, "level_names") <- c("A", "B", "C", "D")  # Should be a vector, not a list
  
  prior_list <- list(group = prior_obj)

  expect_message(
    result <- BayesTools:::.transform_factor_contrasts(
      model_samples, prior_list, transform_factors = TRUE
    ),
    "transformation was applied"
  )

  # Should have 4 columns after transformation (one per level)
  expect_equal(ncol(result), 4)
  expect_true(any(grepl("dif:", colnames(result))))
})


test_that(".filter_parameters removes spike at 0 priors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("point", list(0)),  # spike at 0
    tau = prior("normal", list(1, 1))
  )

  # With remove_spike_0 = TRUE
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = TRUE)
  expect_true("delta" %in% result)
  expect_false("mu" %in% result)
  expect_false("tau" %in% result)

  # With remove_spike_0 = FALSE
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = FALSE)
  expect_equal(length(result), 0)
})


test_that(".filter_parameters removes character specified parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1))
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = c("mu", "tau"), remove_spike_0 = FALSE)
  expect_true("mu" %in% result)
  expect_true("tau" %in% result)
  expect_false("delta" %in% result)
})


test_that(".filter_parameters removes non-formula parameters when TRUE", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with formula attributes
  prior_formula <- prior("normal", list(0, 1))
  attr(prior_formula, "parameter") <- "y"

  prior_list <- list(
    intercept = prior_formula,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = TRUE, remove_spike_0 = FALSE)
  expect_true("sigma" %in% result)
  expect_false("intercept" %in% result)
})


test_that(".filter_parameters removes formula-specific parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_formulas = "y", remove_spike_0 = FALSE)
  expect_true("intercept_y" %in% result)
  expect_true("slope_y" %in% result)
  expect_false("intercept_x" %in% result)
  expect_false("sigma" %in% result)
})


test_that(".filter_parameters keeps only specified parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1))
  )

  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "mu", remove_spike_0 = FALSE)
  expect_false("mu" %in% result)
  expect_true("delta" %in% result)
  expect_true("tau" %in% result)
})


test_that(".filter_parameters keeps only specified formulas", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, keep_formulas = "y", remove_spike_0 = FALSE)
  expect_false("intercept_y" %in% result)
  expect_false("slope_y" %in% result)
  expect_true("intercept_x" %in% result)
  expect_true("sigma" %in% result)
})


test_that(".filter_parameters combines keep_parameters and keep_formulas", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  # Keep formula "y" and parameter "sigma"
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "sigma", keep_formulas = "y", remove_spike_0 = FALSE)
  expect_false("intercept_y" %in% result)
  expect_false("slope_y" %in% result)
  expect_false("sigma" %in% result)
  expect_true("intercept_x" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias is removed", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component (simulating bias)
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  # When bias is removed, PET should also be removed
  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias contains PEESE", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PEESE component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PEESE("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PEESE" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias contains weightfunction", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with weightfunction component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one.sided", list(c(0.05), c(1, 1)), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("omega" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters keeps bias-related parameters when bias is kept", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1)),
    bias = bias_prior
  )

  # When only bias is kept, mu and tau should be removed, but PET should be kept
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "bias", remove_spike_0 = FALSE)
  expect_false("bias" %in% result)
  expect_false("PET" %in% result)
  expect_true("mu" %in% result)
  expect_true("tau" %in% result)
})


test_that(".filter_parameters handles non-mixture bias priors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a single PET prior named bias (not a mixture)
  bias_prior <- prior_PET("normal", list(0, 1))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias is not in keep list", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1)),
    bias = bias_prior
  )

  # When only mu is kept, bias should be removed along with PET
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "mu", remove_spike_0 = FALSE)
  expect_false("mu" %in% result)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_true("tau" %in% result)
})


test_that("helper functions work with runjags estimates extraction", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Test the helper functions with mock data (not full integration)
  # This tests that our refactored code correctly uses the helpers
  
  # Create mock posterior samples
  set.seed(123)
  model_samples <- matrix(rnorm(200), ncol = 2, dimnames = list(NULL, c("mu", "inv_sigma")))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("invgamma", list(1, 1))
  )
  
  # Test that remove_auxiliary_parameters helper works
  cleaned <- BayesTools:::.remove_auxiliary_parameters(model_samples, prior_list, NULL)
  
  # Should remove inv_sigma
  expect_false("inv_sigma" %in% colnames(cleaned$model_samples))
  expect_true("mu" %in% colnames(cleaned$model_samples))
  expect_equal(ncol(cleaned$model_samples), 1)
})
