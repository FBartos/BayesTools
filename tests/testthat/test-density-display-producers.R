skip_if_not_test_profile("fixture")
source(testthat::test_path("common-functions.R"))

test_that("named marginal BF rows retain actual producer fields when inference order changes", {

  skip_if_missing_fits("fit_factor_orthonormal")
  fit <- readRDS(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS"))
  mixed <- as_mixed_posteriors(fit, parameters = "p1")
  marginal <- marginal_posterior(mixed, "p1", prior_samples = TRUE, use_formula = FALSE)
  producers <- Savage_Dickey_BF(marginal, normal_approximation = TRUE, silent = TRUE)
  names(producers) <- names(marginal)
  samples <- list(p1 = marginal)
  for(settings in list(list(), list(BF01 = TRUE), list(logBF = TRUE))){
    reference <- do.call(marginal_estimates_table,
      c(list(samples = samples, inference = list(p1 = producers), parameters = "p1"), settings))
    reordered <- do.call(marginal_estimates_table,
      c(list(samples = samples, inference = list(p1 = rev(producers)), parameters = "p1"), settings))
    expect_equal(as.data.frame(reordered), as.data.frame(reference))
    expect_identical(attr(reordered, "raw_log_BF"), attr(reference, "raw_log_BF"))
    expect_identical(attr(reordered, "numerical_diagnostics"), attr(reference, "numerical_diagnostics"))
    expect_identical(attr(reordered, "warnings"), attr(reference, "warnings"))
  }
  missing <- producers[-1L]
  expect_error(marginal_estimates_table(samples, list(p1 = missing), "p1"),
    "cannot be aligned", fixed = TRUE)
})

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

test_that("actual publication-weight producers own fixed and reference scalar marginals", {

  fixture_names <- c("fit_weightfunction_fixed", "fit_weightfunction_onesided2", "fit_wf_independent_gamma", "fit_bias_petpeese_hetero_wf")
  skip_if_missing_fits(fixture_names)
  fits <- lapply(fixture_names, function(name) readRDS(file.path(temp_fits_dir, paste0(name, ".RDS"))))
  for(i in seq_len(3L)){
    samples <- as_mixed_posteriors(fits[[i]], parameters = "omega")
    atoms <- posterior_metadata(samples$omega, "atoms")
    expect_identical(names(atoms$marginals), colnames(samples$omega))
    expect_equal(unname(atoms$marginals[[1L]]$locations[, 1L]), 1)
    expect_equal(atoms$marginals[[1L]]$mass, 1)
    if(i == 1L){
      expect_equal(unname(atoms$marginals[[2L]]$locations[, 1L]), .5)
      mapped <- posterior_transform(samples, "lin", list(a = 0, b = 2))
      expect_equal(unname(posterior_metadata(mapped$omega, "atoms")$marginals[[2L]]$locations[, 1L]), 1)
      expect_error(plot_posterior(mapped, "omega", individual = TRUE, prior = TRUE),
        class = "BayesTools_formula_prior_density_unavailable")
    }else expect_length(atoms$marginals[[2L]]$mass, 0L)
  }
  bias <- as_mixed_posteriors(fits[[4L]], parameters = "bias")
  atoms <- posterior_metadata(bias$bias, "atoms")
  omega <- grep("^omega\\[", colnames(bias$bias), value = TRUE)
  expect_true(all(!vapply(atoms$marginals[omega], is.null, logical(1))))
  expect_true(all(vapply(atoms$marginals[omega], function(marginal){
    !anyDuplicated(marginal$locations[, 1L])
  }, logical(1))))
  expect_equal(unname(atoms$marginals[[omega[[1L]]]]$locations[, 1L]), 1)
  expect_equal(atoms$marginals[[omega[[1L]]]]$mass, 1)
  pair <- .model_probability_pair(atoms$component_probabilities,
    log(atoms$component_probabilities), "component", "raw")
  for(marginal in atoms$marginals[omega]){
    expect_identical(marginal$component_probabilities, atoms$component_probabilities)
    expect_identical(marginal$component_log_probabilities, pair$logs)
    expect_identical(marginal$model_probability_declaration, pair$declaration)
  }
  simplified <- .simplify_as_mixed_posterior_bias(bias, "omega")
  expect_identical(posterior_metadata(simplified$omega, "atoms")$marginals, atoms$marginals[omega])
  expect_identical(vapply(attr(simplified$omega, "prior_list", exact = TRUE), .prior_model_weight, numeric(1)),
    vapply(.bias_samples_prior_list(bias), .prior_model_weight, numeric(1)))
})
