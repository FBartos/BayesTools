# Tests for posterior extraction helper functions
test_that(".extract_posterior_samples extracts samples correctly", {

  skip_on_cran()
  skip_on_ci()
  skip_on_os("linux")

  # Simple model
  data <- data.frame(y = rnorm(50))
  fit <- JAGS_fit(
    model_syntax = "model{\n for(i in 1:N){\n y[i] ~ dnorm(mu, 1)\n }\n mu ~ dnorm(0, 1)\n}",
    data = list(y = data$y, N = nrow(data)),
    prior_list = list(mu = prior("normal", list(0, 1))),
    chains = 2,
    sample = 100,
    burnin = 100,
    adapt = 100,
    silent = TRUE
  )

  # Test matrix extraction
  samples_matrix <- BayesTools:::.extract_posterior_samples(fit, as_list = FALSE)
  expect_true(is.matrix(samples_matrix))
  expect_equal(ncol(samples_matrix), 1)
  expect_true("mu" %in% colnames(samples_matrix))

  # Test list extraction
  samples_list <- BayesTools:::.extract_posterior_samples(fit, as_list = TRUE)
  expect_true(inherits(samples_list, "mcmc.list"))
  expect_equal(length(samples_list), 2) # 2 chains
})


test_that(".remove_auxiliary_parameters removes invgamma support", {
  skip_on_cran()
  skip_on_ci()
  skip_on_os("linux")

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
  skip_on_ci()
  skip_on_os("linux")

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
  skip_on_ci()
  skip_on_os("linux")

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
  skip_on_ci()
  skip_on_os("linux")

  # Create mock samples with factor
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior_obj, "levels") <- c("A", "B", "C", "D")
  
  prior_list <- list(group = prior_obj)

  result <- BayesTools:::.rename_factor_levels(model_samples, prior_list)

  expect_true("group[B]" %in% colnames(result))
  expect_true("group[C]" %in% colnames(result))
  expect_true("group[D]" %in% colnames(result))
  expect_false("group[1]" %in% colnames(result))
})


test_that(".transform_factor_contrasts transforms orthonormal to differences", {
  skip_on_cran()
  skip_on_ci()
  skip_on_os("linux")

  # Create mock samples with orthonormal contrasts
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("normal", list(0, 1), contrast = "orthonormal")
  attr(prior_obj, "levels") <- c("A", "B", "C", "D")
  
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


test_that("helper functions integrate correctly in runjags_estimates_table", {
  skip_on_cran()
  skip_on_ci()
  skip_on_os("linux")

  # Fit a simple model
  data <- data.frame(y = rnorm(50))
  fit <- JAGS_fit(
    model_syntax = "model{\n for(i in 1:N){\n y[i] ~ dnorm(mu, 1)\n }\n mu ~ dnorm(0, 1)\n}",
    data = list(y = data$y, N = nrow(data)),
    prior_list = list(mu = prior("normal", list(0, 1))),
    chains = 2,
    sample = 100,
    burnin = 100,
    adapt = 100,
    silent = TRUE
  )

  # Test that runjags_estimates_table still works
  summary <- runjags_estimates_table(fit)

  expect_true(inherits(summary, "BayesTools_table"))
  expect_true(inherits(summary, "BayesTools_runjags_summary"))
  expect_true("mu" %in% rownames(summary))
  expect_equal(ncol(summary), 9)  # Mean, SD, lCI, Median, uCI, MCMC_error, MCMC_SD_error, ESS, R_hat
})
