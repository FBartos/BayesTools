skip_if_not_test_profile("fixture")
source(testthat::test_path("common-functions.R"))
skip_if_missing_fits(c("fit_formula_multiplier_state", "fit_formula_multiplier_zero",
  "fit_formula_multiplier_missing"))

test_that("actual multiplier fits preserve predictors, gates and physical rows", {

  fit <- readRDS(file.path(temp_fits_dir, "fit_formula_multiplier_state.RDS"))
  source <- as.matrix(.fit_to_posterior(fit))
  ordinary <- as_mixed_posteriors(fit, "sigma", transform_scaled = TRUE)
  expect_identical(as.numeric(ordinary$sigma), unname(source[, "sigma"]))
  expect_null(.bt_formula_state_get(ordinary))
  original <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"), transform_scaled = TRUE)
  fitted <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x"))
  expect_equal(as.numeric(original$mu_intercept), 5 - 4 * source[, "sigma"], tolerance = 1e-12)
  expect_equal(as.numeric(original$mu_x), rep(.2, nrow(source)), tolerance = 1e-14)
  for(at in c(0, 1)){
    fitted_row <- marginal_posterior(fitted, "mu_intercept", ~ x, list(x = at), prior_samples = TRUE)$intercept
    original_row <- marginal_posterior(original, "mu_intercept", ~ x,
      list(x = 20 + 10 * at), prior_samples = TRUE)$intercept
    expect_equal(as.numeric(fitted_row), source[, paste0("mu[", at + 2, "]")], tolerance = 1e-12)
    expect_identical(as.numeric(fitted_row), as.numeric(original_row))
    expect_identical(.posterior_atoms_get(fitted_row), .posterior_atoms_get(original_row))
  }
  row <- marginal_posterior(fitted, "mu_intercept", ~ x, list(x = 1), prior_samples = TRUE)$intercept
  expect_identical(.posterior_atoms_get(row)$mass, mean(source[, "sigma_indicator"] == 1))
  expect_identical(.posterior_support_get(row)$bounds, c(5, Inf))
  expect_identical(.posterior_components_get(row)$index, as.integer(source[, "sigma_indicator"]))
  conditional <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "sigma"), conditional = "sigma")
  selected <- conditional[c("mu_intercept", "mu_x")]
  positive <- which(source[, "sigma_indicator"] == 2)
  expect_identical(.bt_formula_state_get(selected)$draw_index, positive)
  row <- marginal_posterior(selected, "mu_intercept", ~ x, list(x = 1))$intercept
  expect_equal(as.numeric(row), 5 + 2 * source[positive, "sigma"], tolerance = 1e-12)
  expect_length(.posterior_atoms_get(row)$mass, 0L)
  view <- JAGS_with_draws(fit, fit$mcmc)
  expect_identical(as.numeric(as_mixed_posteriors(view, "mu_intercept", transform_scaled = TRUE)$mu_intercept),
    as.numeric(original$mu_intercept))
  expect_error(JAGS_extend(view), class = "BayesTools_draws_view_unavailable")
})

test_that("actual fixed formula retention excludes unrelated random sources", {

  skip_if_missing_fits("fit_label_random_slope")
  fit <- readRDS(file.path(temp_fits_dir, "fit_label_random_slope.RDS"))
  for(original in c(FALSE, TRUE)){
    mixed <- as_mixed_posteriors(fit, names(attr(fit, "prior_list", exact = TRUE)), transform_scaled = original)
    state <- .bt_formula_state_get(mixed)
    expect_false(any(grepl("__xREx__|__xRE_ALLOCx|__xRE_SUMMARY__", state$models[[1L]]$required)))
    expect_false(any(grepl("__xREx__|__xRE_ALLOCx|__xRE_SUMMARY__", colnames(state$values))))
  }
})

test_that("actual unallocated and missing-coefficient models retain declared masses", {

  names <- c("fit_formula_multiplier_state", "fit_formula_multiplier_zero", "fit_formula_multiplier_missing")
  fits <- lapply(names, function(name) readRDS(file.path(temp_fits_dir, paste0(name, ".RDS"))))
  probabilities <- c(.95, .025, .025)
  models <- lapply(seq_along(fits), function(i){
    list(fit = fits[[i]], marglik = bridgesampling_object(0), prior_weights = probabilities[[i]])
  })
  mixed <- mix_posteriors(models, c("mu_intercept", "mu_x"),
    list(c(FALSE, FALSE, FALSE), c(FALSE, FALSE, TRUE)), seed = 74, n_samples = 2)
  state <- .bt_formula_state_get(mixed)
  expect_identical(state$model, c(1L, 1L))
  expect_equal(state$posterior_model_probabilities, probabilities, tolerance = 1e-15)
  row <- marginal_posterior(mixed, "mu_intercept", ~ x, list(x = 1), prior_samples = TRUE)$intercept
  source <- as.matrix(.fit_to_posterior(fits[[1L]]))
  expected_mass <- .95 * mean(source[, "sigma_indicator"] == 1) + .05
  expect_equal(.posterior_atoms_get(row)$mass, expected_mass, tolerance = 1e-15)
  law <- posterior_metadata(row, "prior_density")
  expect_equal(law$points$p[law$points$x == 5], .95 * .25 + .05, tolerance = 1e-15)
})
