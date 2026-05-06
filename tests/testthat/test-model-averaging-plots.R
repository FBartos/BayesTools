skip_if_not_test_profile(c("visual", "visual-fixture"))

# ============================================================================ #
# TEST FILE: Model Averaging Plots
# ============================================================================ #
#
# PURPOSE:
#   Tests for plot_prior_list, plot_posterior, plot_models, and related
#   visualization functions in model averaging.
#   Semantic plot-data companions live in
#   test-model-averaging-plots-edge-cases.R; this file is for presentation
#   regressions that still need rendered visual coverage.
#
# DEPENDENCIES:
#   - vdiffr: Visual regression testing
#   - rjags: For tests using pre-fitted models
#   - common-functions.R: temp_fits_dir, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("vdiffr"): For visual tests
#   - skip_if_not_visual_fixture_tests(): For visual tests using pre-fitted models
#   - skip_if_no_fits(): For tests using pre-fitted models
#
# MODELS/FIXTURES:
#   - fit_summary0, fit_summary1, fit_orthonormal_0, fit_orthonormal_1
#
# TAGS: @evaluation, @visual, @plots, @model-averaging
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

skip_if_not_installed("vdiffr")


# ============================================================================ #
# SECTION 1: plot_prior_list basic tests
# ============================================================================ #
test_that("plot_prior_list handles simple cases", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Test with a single normal prior
  prior_list_normal <- list(
    p1 = prior("normal", list(0, 1))
  )

  vdiffr::expect_doppelganger("plot-prior-list-single-normal", function() {
    plot_prior_list(prior_list_normal)
  })

  # Test with multiple priors
  prior_list_multi <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(1, 0.5)),
    p3 = prior("cauchy", list(0, 1))
  )

  vdiffr::expect_doppelganger("plot-prior-list-multi", function() {
    plot_prior_list(prior_list_multi)
  })

  # Test with gamma prior
  prior_list_gamma <- list(
    p1 = prior("gamma", list(2, 1))
  )

  vdiffr::expect_doppelganger("plot-prior-list-gamma", function() {
    plot_prior_list(prior_list_gamma)
  })

})


test_that("plot_prior_list handles orthonormal priors", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Create orthonormal factor prior
  prior_orth <- prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  attr(prior_orth, "levels") <- 3

  prior_list <- list(p1 = prior_orth)

  # Base plot
  vdiffr::expect_doppelganger("plot-prior-list-orthonormal-base", function() {
    plot_prior_list(prior_list)
  })

  # ggplot
  vdiffr::expect_doppelganger("plot-prior-list-orthonormal-ggplot", {
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  vdiffr::expect_doppelganger("plot-prior-list-orthonormal2-ggplot", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lty = 2, linetype = 2)
  })

  # Create orthonormal factor prior with spike
  prior_orth0 <- prior_factor("spike", list(0), contrast = "orthonormal")
  attr(prior_orth0, "levels") <- 3

  prior_list0 <- list(p1 = prior_orth0)

  # Base plot
  vdiffr::expect_doppelganger("plot-prior-list-orthonormal-spike", function() {
    plot_prior_list(prior_list0)
  })

  prior_list2 <- list(
    p1 = prior_orth,
    p2 = prior_orth0
  )

  # Base plot with transformation
  vdiffr::expect_doppelganger("plot-prior-list-orthonormal-spike-and-slab", function() {
    suppressMessages(plot_prior_list(prior_list2, transformation = "exp", transformation_settings = TRUE, xlim = c(0.01, 5)))
  })

})


test_that("plot_prior_list handles meandif priors", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Create meandif factor prior
  prior_md <- prior_factor("mnorm", list(mean = 0, sd = 0.5), contrast = "meandif")
  attr(prior_md, "levels") <- 3

  prior_list <- list(p1 = prior_md)

  # Base plot
  vdiffr::expect_doppelganger("plot-prior-list-meandif-base", function() {
    plot_prior_list(prior_list)
  })

  # ggplot
  vdiffr::expect_doppelganger("plot-prior-list-meandif-ggplot", {
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

})


test_that("plot_prior_list handles weightfunction priors", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Create one-sided weightfunction prior
  wf_prior <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))

  prior_list_wf <- list(p1 = wf_prior)

  vdiffr::expect_doppelganger("plot-prior-list-weightfunction", function() {
    plot_prior_list(prior_list_wf)
  })

  # Test ggplot version
  vdiffr::expect_doppelganger("plot-prior-list-weightfunction-ggplot", {
    plot_prior_list(prior_list_wf, plot_type = "ggplot")
  })

})


test_that("scale_y2 is handled correctly for mixed distributions", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Create a list with both continuous and point priors
  prior_list <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("spike", list(0))
  )

  # Base plot should handle dual y-axis
  vdiffr::expect_doppelganger("plot-prior-list-dual-axis", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })

  # ggplot should handle it differently
  vdiffr::expect_doppelganger("plot-prior-list-dual-axis-ggplot", {
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

})


# ============================================================================ #
# SECTION 2: lines_prior_list tests
# ============================================================================ #
test_that("lines_prior_list handles various configurations", {
  skip_if_not_visual_tests()
  set.seed(1)
  prior_list <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 2))
  )

  # Test adding lines to existing plot
  vdiffr::expect_doppelganger("lines-prior-list-add", function() {
    plot(NULL, xlim = c(-5, 5), ylim = c(0, 0.5), xlab = "", ylab = "")
    lines_prior_list(prior_list, col = "red", lwd = 2)
  })

  # Test with custom xlim
  vdiffr::expect_doppelganger("lines-prior-list-xlim", function() {
    plot(NULL, xlim = c(-3, 3), ylim = c(0, 0.5), xlab = "", ylab = "")
    lines_prior_list(prior_list, xlim = c(-3, 3), col = "blue")
  })

})


# ============================================================================ #
# SECTION 3: geom_prior_list tests
# ============================================================================ #
test_that("geom_prior_list handles various configurations", {
  skip_if_not_visual_tests()
  set.seed(1)
  prior_list <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("spike", list(0.5))
  )

  # Test adding to ggplot
  vdiffr::expect_doppelganger("geom-prior-list-add", {
    ggplot2::ggplot() +
      ggplot2::xlim(-4, 4) +
      ggplot2::ylim(0, 1) +
      geom_prior_list(prior_list, col = "red")
  })

})


# ============================================================================ #
# SECTION 4: plot_posterior tests
# ============================================================================ #
test_that("plot_posterior handles various sample types", {
  set.seed(1)
  skip_if_not_visual_fixture_tests()
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load fits
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE), "omega" = c(TRUE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  # Test simple posterior plot
  vdiffr::expect_doppelganger("plot-posterior-simple", function() {
    plot_posterior(mixed_posteriors, "m")
  })

  # Test with prior overlay
  vdiffr::expect_doppelganger("plot-posterior-with-prior", function() {
    plot_posterior(mixed_posteriors, "m", prior = TRUE)
  })

  # Test ggplot version
  vdiffr::expect_doppelganger("plot-posterior-ggplot", {
    plot_posterior(mixed_posteriors, "m", plot_type = "ggplot")
  })

  # Test with custom xlim
  vdiffr::expect_doppelganger("plot-posterior-xlim", function() {
    plot_posterior(mixed_posteriors, "m", xlim = c(-2, 2))
  })

})


test_that("plot_posterior handles weightfunction posteriors", {
  set.seed(1)
  skip_if_not_visual_fixture_tests()
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load fits
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1)
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE), "omega" = c(TRUE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  # Test weightfunction posterior plot
  vdiffr::expect_doppelganger("plot-posterior-omega", function() {
    plot_posterior(mixed_posteriors, "omega", n_points = 50, n_samples = 500)
  })

})


# ============================================================================ #
# SECTION 5: plot_models tests
# ============================================================================ #
test_that("plot_models handles various configurations", {
  set.seed(1)
  skip_if_not_visual_fixture_tests()
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load fits
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary0)),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary1))
  )
  models <- models_inference(models)

  inference <- ensemble_inference(
    model_list = models,
    parameters = c("m"),
    is_null_list = list("m" = c(FALSE, FALSE))
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = c("m"),
    is_null_list = list("m" = c(FALSE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  # Test basic plot_models
  vdiffr::expect_doppelganger("plot-models-basic", function() {
    plot_models(models, mixed_posteriors, inference, "m")
  })

  # Test ggplot version
  vdiffr::expect_doppelganger("plot-models-ggplot", {
    plot_models(models, mixed_posteriors, inference, "m", plot_type = "ggplot")
  })

})


test_that("plot_models handles order argument", {
  set.seed(1)
  skip_if_not_visual_fixture_tests()
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load fits
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary0)),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary1))
  )
  models <- models_inference(models)

  inference <- ensemble_inference(
    model_list = models,
    parameters = c("m"),
    is_null_list = list("m" = c(FALSE, FALSE))
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = c("m"),
    is_null_list = list("m" = c(FALSE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  # Test with order = decreasing by estimate
  vdiffr::expect_doppelganger("plot-models-order-decreasing-estimate", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", order = list("decreasing", "estimate"))
  })

  # Test with order = increasing by BF
  vdiffr::expect_doppelganger("plot-models-order-increasing-bf", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", order = list("decreasing", "BF"))
  })

  # Test with order = decreasing by probability
  vdiffr::expect_doppelganger("plot-models-order-decreasing-prob", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", order = list("decreasing", "probability"))
  })

  # Test with transformation
  vdiffr::expect_doppelganger("plot-models-order-trans", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", transformation = "exp")
  })

  # Test with transformation and prior
  vdiffr::expect_doppelganger("plot-models-order-trans-prior", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", prior = TRUE, transformation = "exp")
  })

  # Test with transformation ggplot
  vdiffr::expect_doppelganger("plot-models-order-trans-ggplot", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", transformation = "exp", plot_type = "ggplot")
  })

  # Test with transformation and prior ggplot
  vdiffr::expect_doppelganger("plot-models-order-trans-prior-ggplot", function() {
    BayesTools::plot_models(models, mixed_posteriors, inference, "m", prior = TRUE, transformation = "exp", plot_type = "ggplot")
  })

})


test_that("plot_models handles orthonormal priors", {
  set.seed(1)
  skip_if_not_visual_fixture_tests()
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load orthonormal models with marginal likelihoods
  fit_orthonormal_0 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0.RDS"))
  marglik_orthonormal_0 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_0.RDS"))

  fit_orthonormal_1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik_orthonormal_1 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_1.RDS"))

  models <- list(
    list(fit = fit_orthonormal_0, marglik = marglik_orthonormal_0, prior_weights = 1,
         fit_summary = suppressMessages(runjags_estimates_table(fit_orthonormal_0, transform_factors = TRUE))),
    list(fit = fit_orthonormal_1, marglik = marglik_orthonormal_1, prior_weights = 1,
         fit_summary = suppressMessages(runjags_estimates_table(fit_orthonormal_1, transform_factors = TRUE)))
  )
  models <- models_inference(models)

  # Get factor parameter names from the model
  prior_list <- attr(fit_orthonormal_1, "prior_list")
  factor_params <- names(prior_list)[sapply(prior_list, is.prior.factor)]

  inference <- ensemble_inference(
    model_list = models,
    parameters = factor_params,
    is_null_list = setNames(list(c(TRUE, FALSE)), factor_params)
  )

  mixed_posteriors <- mix_posteriors(
    model_list = models,
    parameters = factor_params,
    is_null_list = setNames(list(c(TRUE, FALSE)), factor_params),
    seed = 1,
    n_samples = 1000
  )

  # Test with orthonormal priors
  vdiffr::expect_doppelganger("plot-models-orthonormal", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4), mfrow = c(3, 1))
    BayesTools::plot_models(models, mixed_posteriors, inference, factor_params)
  })

  vdiffr::expect_doppelganger("plot-models-orthonormal-2", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4), mfrow = c(3, 1))
    BayesTools::plot_models(models, mixed_posteriors, inference, factor_params, transformation = "exp")
  })

  vdiffr::expect_doppelganger("plot-models-orthonormal-3", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4), mfrow = c(3, 1))
    BayesTools::plot_models(models, mixed_posteriors, inference, factor_params, transformation = "exp", prior = TRUE)
  })

})


# ============================================================================ #
# SECTION 6: .plot_prior_list.factor tests
# ============================================================================ #
test_that(".plot_prior_list.factor handles point priors within factor", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Test factor prior with spike
  prior_spike <- prior("spike", list(0))
  prior_factor_treat <- prior_factor("normal", list(mean = 0, sd = 1), contrast = "treatment")
  attr(prior_factor_treat, "levels") <- 3

  prior_list <- list(p1 = prior_spike, p2 = prior_factor_treat)

  vdiffr::expect_doppelganger("plot-factor-with-spike", function() {
    plot_prior_list(prior_list)
  })

  vdiffr::expect_doppelganger("plot-factor-with-spike-trans", function() {
    plot_prior_list(prior_list, transformation = "tanh")
  })

  vdiffr::expect_doppelganger("plot-factor-with-spike-trans-settings", function() {
    plot_prior_list(prior_list, transformation = "tanh", transformation_settings = T, xlim = c(-0.5, 0.5))
  })

})


test_that(".plot_prior_list.factor handles transformation", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Create treatment factor prior with normal distribution
  prior_treat <- prior_factor("normal", list(mean = 0, sd = 0.5), contrast = "treatment")
  attr(prior_treat, "levels") <- 3

  prior_list <- list(p1 = prior_treat)

  vdiffr::expect_doppelganger("plot-factor-transformation", function() {
    plot_prior_list(prior_list, transformation = "exp")
  })

})

# ============================================================================ #
# SECTION: get_scale_transformation with plotting
# ============================================================================ #

test_that("exp_lin transformation functions are defined correctly", {
  skip_if_not_visual_tests()
  # Test that exp_lin transformation is correctly defined
  # (used for log-intercept unscaling)

  # Get the transformation functions
  trans_funcs <- BayesTools:::.density.prior_transformation_functions("exp_lin")

  # Verify the functions exist
  expect_true(is.function(trans_funcs$fun))
  expect_true(is.function(trans_funcs$inv))
  expect_true(is.function(trans_funcs$jac))

  # Test the transformation: exp(a + b * log(x))
  x <- 2
  a <- 0.5
  b <- 1.3

  # fun: compare against the closed form
  expected <- exp(a + b * log(x))
  expect_equal(trans_funcs$fun(x, a, b), expected)

  # inv: should reverse the transformation
  y <- trans_funcs$fun(x, a, b)
  expect_equal(trans_funcs$inv(y, a, b), x, tolerance = 1e-10)

  # jac: derivative of the transformation at the original value
  expect_equal(trans_funcs$jac(x, a, b), b * y / x)
  expect_equal(trans_funcs$jac(x, 0, 1), 1)
})

test_that("linear transformation matches expected behavior", {
  skip_if_not_visual_tests()
  set.seed(1)

  # Create a normal prior
  prior_list <- list(p1 = prior("normal", list(0, 1)))

  # Apply linear transformation: a + b*x with a=0, b=0.5
  # This should compress the distribution by half
  vdiffr::expect_doppelganger("plot-normal-lin-compress", function() {
    plot_prior_list(prior_list,
                    transformation = "lin",
                    transformation_arguments = list(a = 0, b = 0.5))
  })

  # Apply linear transformation with offset: a + b*x with a=2, b=0.5
  # This should compress and shift
  vdiffr::expect_doppelganger("plot-normal-lin-shift-compress", function() {
    plot_prior_list(prior_list,
                    transformation = "lin",
                    transformation_arguments = list(a = 2, b = 0.5))
  })
})


# ============================================================================ #
# SECTION: transform_scaled visual tests
# ============================================================================ #
# These tests use pre-fitted regression models with formula_scale to visually
# verify that the transform_scaled feature correctly transforms priors and
# posteriors from standardized to original scale.

.make_fake_scaled_formula_fit <- function(prior_list, posterior_samples, formula_scale) {
  fit <- coda::mcmc(posterior_samples)
  class(fit) <- c("BayesTools_fit", class(fit))

  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- formula_scale

  return(fit)
}

.integrate_density_mass <- function(x, y) {
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}


test_that("transform_scaled is auto-detected from samples attribute", {
  skip_if_not_visual_fixture_tests()
  skip_if_no_fits()

  # Load a model with formula_scale
  fit_path <- file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS")
  skip_if_not(file.exists(fit_path), "Pre-fitted model not available")

  fit <- readRDS(fit_path)

  # Extract with transform_scaled = TRUE
  samples_scaled <- as_mixed_posteriors(fit, parameters = "mu_intercept", transform_scaled = TRUE)

  # Verify the attribute is set

  expect_true(isTRUE(attr(samples_scaled, "transform_scaled")))
  expect_false(is.null(attr(samples_scaled, "prior_densities")))
  expect_null(attr(samples_scaled, "prior_samples"))

  # Extract without transform_scaled
  samples_unscaled <- as_mixed_posteriors(fit, parameters = "mu_intercept", transform_scaled = FALSE)

  # Verify the attribute is NOT set
  expect_null(attr(samples_unscaled, "transform_scaled"))
})


test_that("linear prior density helper preserves transformed point masses", {
  skip_if_not_visual_tests()
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    BayesTools:::.prior_linear_density_point(0),
    n_points = 128
  )

  expect_true(any(vapply(plot_data, inherits, logical(1), what = "density.prior.point")))
  expect_false(any(vapply(plot_data, inherits, logical(1), what = "density.prior.simple")))
  expect_equal(plot_data[[1]]$x, 0)
  expect_equal(plot_data[[1]]$y, 1)
})

test_that("transform_scaled helper preserves total mass for mixed point-continuous priors", {
  skip_if_not_visual_tests()
  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(1, 0.25), prior_weights = 4)
      ),
      is_null = c(TRUE, FALSE)
    )),
    weights = c(x = 1),
    n_grid = 512
  )
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    prior_density,
    n_points = 512
  )

  density_data <- plot_data[[which(vapply(plot_data, inherits, logical(1), what = "density.prior.simple"))]]
  point_mass <- sum(vapply(
    plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.point")],
    function(item) item$y,
    numeric(1)
  ))

  expect_equal(.integrate_density_mass(density_data$x, density_data$y) + point_mass, 1, tolerance = 0.03)
})

test_that("transform_scaled helper preserves total mass under user transformations", {
  skip_if_not_visual_tests()
  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 3),
        prior("normal", list(1.5, 0.4), prior_weights = 17)
      ),
      is_null = c(TRUE, FALSE)
    )),
    weights = c(x = 1),
    n_grid = 512
  )
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    prior_density,
    n_points = 512,
    transformation = "lin",
    transformation_arguments = list(a = 2, b = 0.5)
  )

  density_data <- plot_data[[which(vapply(plot_data, inherits, logical(1), what = "density.prior.simple"))]]
  point_data <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.point")]
  point_mass <- sum(vapply(point_data, function(item) item$y, numeric(1)))

  expect_equal(point_data[[1]]$x, 2)
  expect_equal(.integrate_density_mass(density_data$x, density_data$y) + point_mass, 1, tolerance = 0.03)
})

test_that("plot_posterior errors when transformed prior densities are missing", {
  skip_if_not_visual_tests()
  sample_entry <- structure(
    rnorm(64),
    class = c("mixed_posteriors.formula", "mixed_posteriors.simple", "mixed_posteriors"),
    formula_parameter = "mu",
    prior_list = list(prior("normal", list(0, 1)))
  )
  samples <- structure(
    list(mu_x1 = sample_entry),
    class = c("as_mixed_posteriors", "mixed_posteriors"),
    transform_scaled = TRUE
  )

  expect_error(
    plot_posterior(samples, "mu_x1", prior = TRUE, plot_type = "ggplot"),
    "no prior densities found"
  )
})


test_that("transform_scaled visual: spike prior remains atomic", {
  skip_if_not_visual_tests()
  skip_on_cran()

  set.seed(123)

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x1 = prior("spike", list(0))
  )
  for(i in seq_along(prior_list)){
    attr(prior_list[[i]], "parameter") <- "mu"
  }

  formula_scale <- list(mu = list(mu_x1 = list(mean = 5, sd = 2)))
  attr(formula_scale$mu, "log_intercept") <- FALSE

  posterior_samples <- cbind(
    mu_intercept = rnorm(4000, 0, 1),
    mu_x1 = rnorm(4000, 0.25, 0.20)
  )

  fit <- .make_fake_scaled_formula_fit(prior_list, posterior_samples, formula_scale)
  samples_scaled <- as_mixed_posteriors(
    fit,
    parameters = "mu_x1",
    transform_scaled = TRUE,
    n_prior_samples = 20000
  )

  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(samples_scaled, "prior_densities")$mu_x1, 0),
    1
  )

  vdiffr::expect_doppelganger("transform-scaled-spike-only-prior", function() {
    plot_posterior(
      samples_scaled,
      "mu_x1",
      prior = TRUE,
      main = "Spike Prior (Original Scale)",
      dots_prior = list(col = "grey"),
      xlim = c(-1, 1)
    )
  })
})


test_that("transform_scaled visual: conditional mixture prior removes spike", {
  skip_if_not_visual_tests()
  skip_on_cran()

  set.seed(124)

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x1 = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(0.6, 0.4), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE)
    )
  )
  for(i in seq_along(prior_list)){
    attr(prior_list[[i]], "parameter") <- "mu"
  }

  formula_scale <- list(mu = list(mu_x1 = list(mean = 5, sd = 2)))
  attr(formula_scale$mu, "log_intercept") <- FALSE

  indicator <- sample(c(1L, 2L), size = 4000, replace = TRUE)
  posterior_samples <- cbind(
    mu_intercept = rnorm(4000, 0, 1),
    mu_x1 = ifelse(indicator == 1L, 0, rnorm(4000, 0.6, 0.4)),
    mu_x1_indicator = indicator
  )

  fit <- .make_fake_scaled_formula_fit(prior_list, posterior_samples, formula_scale)
  samples_unconditional <- as_mixed_posteriors(
    fit,
    parameters = "mu_x1",
    transform_scaled = TRUE,
    n_prior_samples = 20000
  )
  samples_conditional <- as_mixed_posteriors(
    fit,
    parameters = "mu_x1",
    conditional = "mu_x1",
    transform_scaled = TRUE,
    n_prior_samples = 20000
  )

  expect_gt(
    BayesTools:::.prior_linear_density_point_mass(attr(samples_unconditional, "prior_densities")$mu_x1, 0),
    0
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(samples_conditional, "prior_densities")$mu_x1, 0),
    0
  )

  vdiffr::expect_doppelganger("transform-scaled-conditional-mixture-prior", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))

    plot_posterior(
      samples_unconditional,
      "mu_x1",
      prior = TRUE,
      main = "Mixture Prior (Unconditional)",
      dots_prior = list(col = "grey"),
      xlim = c(-1, 2)
    )

    plot_posterior(
      samples_conditional,
      "mu_x1",
      prior = TRUE,
      main = "Mixture Prior (Conditional Alt)",
      dots_prior = list(col = "grey"),
      xlim = c(-1, 2)
    )
  })
})


test_that("transform_scaled visual: auto-scaled continuous predictors intercept", {
  skip_if_not_visual_fixture_tests()
  skip_on_cran()
  skip_if_no_fits()

  fit_path <- file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS")
  skip_if_not(file.exists(fit_path), "Pre-fitted model not available")

  fit <- readRDS(fit_path)

  # Extract posteriors with and without transform_scaled
  samples_scaled <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
                                        transform_scaled = TRUE)
  samples_unscaled <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
                                          transform_scaled = FALSE)

  # Visual test: intercept - scaled (left) vs original (right)
  vdiffr::expect_doppelganger("transform-scaled-intercept-comparison", function() {
    par(mfrow = c(1, 2))

    # Left: Standardized scale
    plot_posterior(samples_unscaled, "mu_intercept", prior = TRUE,
                   main = "Intercept (Standardized Scale)", dots_prior = list(col = "grey"))

    # Right: Original scale (auto-detected from samples)
    plot_posterior(samples_scaled, "mu_intercept", prior = TRUE,
                   main = "Intercept (Original Scale)", dots_prior = list(col = "grey"))
  })
})


test_that("transform_scaled visual: auto-scaled continuous predictor coefficient", {
  skip_if_not_visual_fixture_tests()
  skip_on_cran()
  skip_if_no_fits()

  fit_path <- file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS")
  skip_if_not(file.exists(fit_path), "Pre-fitted model not available")

  fit      <- readRDS(fit_path)

  samples_scaled <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
                                        transform_scaled = TRUE, n_prior_samples = 1e5)
  samples_unscaled <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
                                          transform_scaled = FALSE)

  # Visual test: coefficient x_cont1 - scaled (left) vs original (right)
  vdiffr::expect_doppelganger("transform-scaled-coef-x_cont1-comparison", function() {
    par(mfrow = c(1, 2))

    # Left: Standardized scale
    plot_posterior(samples_unscaled, "mu_x_cont1", prior = TRUE,
                   main = "x_cont1 (Standardized Scale)", dots_prior = list(col = "grey"))

    # Right: Original scale (auto-detected from samples)
    plot_posterior(samples_scaled, "mu_x_cont1", prior = TRUE,
                   main = "x_cont1 (Original Scale)", dots_prior = list(col = "grey"))
  })

  # Visual test: coefficient x_cont2
  vdiffr::expect_doppelganger("transform-scaled-coef-x_cont2-comparison", function() {
    par(mfrow = c(1, 2))

    plot_posterior(samples_unscaled, "mu_x_cont2", prior = TRUE,
                   main = "x_cont2 (Standardized Scale)", dots_prior = list(col = "grey"))

    plot_posterior(samples_scaled, "mu_x_cont2", prior = TRUE,
                   main = "x_cont2 (Original Scale)", dots_prior = list(col = "grey"))
  })
})


test_that("transform_scaled visual: all parameters side-by-side", {
  skip_if_not_visual_fixture_tests()
  skip_on_cran()
  skip_if_no_fits()

  fit_path <- file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS")
  skip_if_not(file.exists(fit_path), "Pre-fitted model not available")

  fit      <- readRDS(fit_path)

  samples_scaled <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
                                        transform_scaled = TRUE)
  samples_unscaled <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_cont2"),
                                          transform_scaled = FALSE)

  # Visual test: 3x2 grid showing all parameters
  vdiffr::expect_doppelganger("transform-scaled-all-params-grid", function() {
    par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

    # Row 1: Intercept
    plot_posterior(samples_unscaled, "mu_intercept", prior = TRUE, main = "Intercept (Scaled)", dots_prior = list(col = "grey"))
    plot_posterior(samples_scaled, "mu_intercept", prior = TRUE, main = "Intercept (Original)", dots_prior = list(col = "grey"))

    # Row 2: x_cont1
    plot_posterior(samples_unscaled, "mu_x_cont1", prior = TRUE, main = "x_cont1 (Scaled)", dots_prior = list(col = "grey"))
    plot_posterior(samples_scaled, "mu_x_cont1", prior = TRUE, main = "x_cont1 (Original)", dots_prior = list(col = "grey"))

    # Row 3: x_cont2
    plot_posterior(samples_unscaled, "mu_x_cont2", prior = TRUE, main = "x_cont2 (Scaled)", dots_prior = list(col = "grey"))
    plot_posterior(samples_scaled, "mu_x_cont2", prior = TRUE, main = "x_cont2 (Original)", dots_prior = list(col = "grey"))
  })
})


test_that("transform_scaled visual: dual parameter regression with log(intercept)", {
  skip_if_not_visual_fixture_tests()
  skip_on_cran()
  skip_if_no_fits()

  fit_path <- file.path(temp_fits_dir, "fit_dual_param_regression.RDS")
  fit <- readRDS(fit_path)

  # Get available mu parameters (those with formula_scale applied)
  params <- names(attr(fit, "prior_list"))
  sigma_params <- params[grepl("^log_sigma_", params)]

  samples_scaled <- as_mixed_posteriors(fit, parameters = sigma_params, transform_scaled = TRUE)
  samples_unscaled <- as_mixed_posteriors(fit, parameters = sigma_params, transform_scaled = FALSE)

  # Visual test: intercept for dual-parameter model
  vdiffr::expect_doppelganger("transform-scaled-dual-param-intercept", function() {
    par(mfrow = c(2, 2))
    plot_posterior(samples_unscaled, "log_sigma_intercept", prior = TRUE, main = "Dual: Intercept (Scaled)", dots_prior = list(col = "grey"), xlim = c(0, 1))
    plot_posterior(samples_scaled, "log_sigma_intercept", prior = TRUE, main = "Dual: Intercept (Original)", dots_prior = list(col = "grey"), xlim = c(0, 1))
    plot_posterior(samples_unscaled, "log_sigma_x_sigma", prior = TRUE, main = "Dual: Slope (Scaled)", dots_prior = list(col = "grey"), xlim = c(-1, 1))
    plot_posterior(samples_scaled, "log_sigma_x_sigma", prior = TRUE, main = "Dual: Slope (Original)", dots_prior = list(col = "grey"), xlim = c(-1, 1))
  })
})

