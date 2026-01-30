# ============================================================================ #
# TEST FILE: Model Averaging Plots
# ============================================================================ #
#
# PURPOSE:
#   Tests for plot_prior_list, plot_posterior, plot_models, and related
#   visualization functions in model averaging.
#
# DEPENDENCIES:
#   - vdiffr: Visual regression testing
#   - rjags: For tests using pre-fitted models
#   - common-functions.R: temp_fits_dir, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("vdiffr"): For visual tests
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
  set.seed(1)

  # Create one-sided weightfunction prior
  wf_prior <- prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))

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
  b <- 1

  # fun: exp(0.5 + 1 * log(2)) = exp(0.5 + 0.693) = exp(1.193) ≈ 3.30
  expected <- exp(a + b * log(x))
  expect_equal(trans_funcs$fun(x, a, b), expected)

  # inv: should reverse the transformation
  y <- trans_funcs$fun(x, a, b)
  expect_equal(trans_funcs$inv(y, a, b), x, tolerance = 1e-10)

  # jac: 1 / (b * x)
  expect_equal(trans_funcs$jac(x, a, b), 1 / (b * x))
})


test_that("get_scale_transformation prior/posterior plots with cached fit", {
  set.seed(1)
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load the auto-scaled fit
  fit <- readRDS(file.path(temp_fits_dir, "fit_formula_auto_scaled.RDS"))
  formula_scale <- attr(fit, "formula_scale")

  # Create priors matching those used in test-00-model-fits.R
  prior_x_cont   <- prior("normal", list(0, 1))
  prior_x_int    <- prior("normal", list(0, 1))
  prior_intercept <- prior("normal", list(0, 5))

  # Get transformations for the scaled parameters
  # get_scale_transformation now returns MARGINAL transformation by default,

  # using L2 norm of transformation coefficients for correct prior visualization
  trans_x1  <- get_scale_transformation(fit, "mu_x_cont1", formula_scale)
  trans_int <- get_scale_transformation(fit, "mu_x_cont1__xXx__x_cont2", formula_scale)
  trans_intercept <- get_scale_transformation(fit, "mu_intercept", formula_scale)

  # Extract posterior samples and transform back to original scale
  posterior <- as.matrix(coda::as.mcmc.list(fit))
  posterior_transformed <- transform_scale_samples(posterior, formula_scale)

  posterior_x1_unscaled <- posterior_transformed[, "mu_x_cont1"]
  posterior_int_unscaled <- posterior_transformed[, "mu_x_cont1__xXx__x_cont2"]
  posterior_intercept_unscaled <- posterior_transformed[, "mu_intercept"]

  # 1x2 layout for x_cont1: scaled | unscaled
  vdiffr::expect_doppelganger("scale-transform-x1-1x2", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    # Scaled panel: histogram of scaled posterior with posterior density and prior
    hist(posterior[, "mu_x_cont1"], breaks = 50, probability = TRUE,
         main = "x_cont1 (scaled)", xlab = "")
    lines(stats::density(posterior[, "mu_x_cont1"]), col = "blue", lwd = 2)
    lines(prior_x_cont, col = "red", lwd = 2)

    # Unscaled panel: histogram of unscaled posterior with density and transformed prior
    hist(posterior_x1_unscaled, breaks = 50, probability = TRUE,
         main = "x_cont1 (unscaled)", xlab = "")
    lines(stats::density(posterior_x1_unscaled), col = "blue", lwd = 2)  # Use transformed posterior
    lines(stats::density(
      trans_x1$transformation_arguments$a + trans_x1$transformation_arguments$b * posterior[, "mu_x_cont1"]),
      col = "blue", lwd = 2, lty = 2)
    lines(prior_x_cont,
          transformation = trans_x1$transformation,
          transformation_arguments = trans_x1$transformation_arguments,
          col = "red", lwd = 2)
  })

  # 1x2 layout for interaction: scaled | unscaled
  vdiffr::expect_doppelganger("scale-transform-interaction-1x2", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    # Scaled panel: histogram of scaled posterior with posterior density and prior
    hist(posterior[, "mu_x_cont1__xXx__x_cont2"], breaks = 30, probability = TRUE,
         main = "x1:x2 (scaled)", xlab = "")
    lines(stats::density(posterior[, "mu_x_cont1__xXx__x_cont2"]), col = "blue", lwd = 2)
    lines(prior_x_int, col = "red", lwd = 2)

    # Unscaled panel: histogram of unscaled posterior with density and transformed prior
    hist(posterior_int_unscaled, breaks = 30, probability = TRUE,
         main = "x1:x2 (unscaled)", xlab = "")
    lines(stats::density(posterior_int_unscaled), col = "blue", lwd = 2)
    lines(prior_x_int,
          transformation = trans_int$transformation,
          transformation_arguments = trans_int$transformation_arguments,
          col = "red", lwd = 2)
  })

  # 1x2 layout for intercept: scaled | unscaled
  # Intercept uses MARGINAL transformation (L2 norm of entire row)
  vdiffr::expect_doppelganger("scale-transform-intercept-1x2", function() {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    # Scaled panel: histogram of scaled posterior with posterior density and prior
    hist(posterior[, "mu_intercept"], breaks = 30, probability = TRUE,
         main = "intercept (scaled)", xlab = "")
    lines(stats::density(posterior[, "mu_intercept"]), col = "blue", lwd = 2)
    lines(prior_intercept, col = "red", lwd = 2)

    # Unscaled panel: histogram of unscaled posterior with density and transformed prior
    hist(posterior_intercept_unscaled, breaks = 30, probability = TRUE,
         main = "intercept (unscaled)", xlab = "")
    lines(stats::density(posterior_intercept_unscaled), col = "blue", lwd = 2)
    lines(prior_intercept,
          transformation = trans_intercept$transformation,
          transformation_arguments = trans_intercept$transformation_arguments,
          col = "red", lwd = 2)
  })
})


test_that("get_scale_transformation with dual parameter model (log intercept)", {
  set.seed(1)
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  # Load the dual parameter regression fit (has log(intercept) for log_sigma)
  fit <- readRDS(file.path(temp_fits_dir, "fit_dual_param_regression.RDS"))
  formula_scale <- attr(fit, "formula_scale")

  # Create priors matching those used in test-00-model-fits.R
  prior_mu_x      <- prior("normal", list(0, 1))
  prior_ls_x      <- prior("normal", list(0, 0.5))

  # Get transformations for mu (standard intercept) - fit must be first argument
  trans_mu_x <- get_scale_transformation(fit, "mu_x_mu", formula_scale)
  trans_mu_int <- get_scale_transformation(fit, "mu_intercept", formula_scale)

  # Get transformations for log_sigma (log intercept)
  trans_ls_x <- get_scale_transformation(fit, "log_sigma_x_sigma", formula_scale)
  trans_ls_int <- get_scale_transformation(fit, "log_sigma_intercept", formula_scale)

  # Verify log_sigma intercept uses exp_lin transformation
  expect_equal(trans_ls_int$transformation, "exp_lin")

  # Extract posterior samples
  posterior <- as.matrix(coda::as.mcmc.list(fit))

  # Plot mu coefficient: standardized vs unscaled
  vdiffr::expect_doppelganger("dual-mu-x-standardized", function() {
    hist(posterior[, "mu_x_mu"], breaks = 30, probability = TRUE,
         main = "mu_x_mu (standardized)", xlab = "")
    lines(prior_mu_x, col = "red", lwd = 2)
  })

  posterior_mu_x_unscaled <- trans_mu_x$transformation_arguments$a +
    trans_mu_x$transformation_arguments$b * posterior[, "mu_x_mu"]

  vdiffr::expect_doppelganger("dual-mu-x-unscaled", function() {
    hist(posterior_mu_x_unscaled, breaks = 30, probability = TRUE,
         main = "mu_x_mu (unscaled)", xlab = "")
    lines(prior_mu_x,
          transformation = trans_mu_x$transformation,
          transformation_arguments = trans_mu_x$transformation_arguments,
          col = "red", lwd = 2)
  })

  # Plot log_sigma coefficient: standardized vs unscaled
  vdiffr::expect_doppelganger("dual-log-sigma-x-standardized", function() {
    hist(posterior[, "log_sigma_x_sigma"], breaks = 30, probability = TRUE,
         main = "log_sigma_x_sigma (standardized)", xlab = "")
    lines(prior_ls_x, col = "red", lwd = 2)
  })

  posterior_ls_x_unscaled <- trans_ls_x$transformation_arguments$a +
    trans_ls_x$transformation_arguments$b * posterior[, "log_sigma_x_sigma"]

  vdiffr::expect_doppelganger("dual-log-sigma-x-unscaled", function() {
    hist(posterior_ls_x_unscaled, breaks = 30, probability = TRUE,
         main = "log_sigma_x_sigma (unscaled)", xlab = "")
    lines(prior_ls_x,
          transformation = trans_ls_x$transformation,
          transformation_arguments = trans_ls_x$transformation_arguments,
          col = "red", lwd = 2)
  })
})


test_that("linear transformation matches expected behavior", {
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
