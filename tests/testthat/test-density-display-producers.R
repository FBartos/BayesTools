skip_if_not_test_profile("fixture")
source(testthat::test_path("common-functions.R"))

test_that("actual scalar producers match precomputed plot identity before transformations", {

  skip_if_missing_fits("fit_simple_normal")
  fit <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  samples <- as_mixed_posteriors(fit, parameters = "m")
  grid <- seq(-4, 4, length.out = 101)
  curve <- posterior_density_attribute(grid, stats::dnorm(grid), "reference", "precomputed", parameter = "m")
  posterior_metadata(samples$m, "posterior_density") <- curve
  data <- .plot_data_samples.simple(samples, "m", 101, "lin", list(a = 10, b = 2), FALSE, "precomputed")
  expect_equal(data$density$x, 10 + 2 * grid)
  expect_equal(data$density$y, stats::dnorm(grid) / 2)
  wrong <- curve
  wrong$parameter <- "another_target"
  posterior_metadata(samples$m, "posterior_density") <- wrong
  expect_error(plot_posterior(samples, "m", prior = FALSE, density_method = "precomputed"),
    "selected parameter and condition", fixed = TRUE)
  expect_no_error(plot_posterior(samples, "m", prior = FALSE, density_method = "KDE"))
})

test_that("actual transformed scalar producers attach their declared source prior law", {

  skip_if_missing_fits("fit_simple_normal")
  fit <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  samples <- as_mixed_posteriors(fit, parameters = "m")
  original_prior <- attr(samples$m, "prior_list")
  original_prior <- if(is.prior(original_prior)) original_prior else original_prior[[1L]]
  posterior_metadata(samples$m, "prior_density") <- NULL
  set.seed(18)
  rng_before <- .Random.seed
  mapped <- posterior_transform(samples, "lin", list(a = 10, b = 2))
  expect_identical(.Random.seed, rng_before)
  law <- posterior_metadata(mapped$m, "prior_density")
  expect_s3_class(law, "prior_linear_density")
  expected <- prior_density_ordinate(original_prior, 0)
  actual <- prior_density_ordinate(law, 10)
  expect_equal(actual$log_density, expected$log_density - log(2), tolerance = 1e-12)
  # The attached law is consumed independently of plotting.
  marginal <- marginal_posterior(samples, "m", prior_samples = TRUE)
  transformed_marginal <- posterior_transform(marginal, "lin", list(a = 10, b = 2))
  expect_s3_class(Savage_Dickey_BF(transformed_marginal, null_hypothesis = 10,
    normal_approximation = TRUE), "BayesTools_BF")
})
