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
