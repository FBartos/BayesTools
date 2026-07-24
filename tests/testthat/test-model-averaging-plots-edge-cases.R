skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Model Averaging Plots Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for plot functions including input validation and
#   error handling for invalid prior configurations.
#
# DEPENDENCIES:
#   - None (pure R testing)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @edge-cases, @plots, @input-validation
# ============================================================================ #

.model_plot_density_mass <- function(density) {
  sum(diff(density$x) * (head(density$y, -1) + tail(density$y, -1)) / 2)
}


# ============================================================================ #
# SECTION 1: plot_prior_list input validation
# ============================================================================ #
test_that("plot_prior_list rejects non-list input", {

  expect_error(
    plot_prior_list(prior("normal", list(0, 1))),
    "must be a list of priors"
  )

})


test_that("plot_prior_list rejects PET-PEESE without prior_list_mu", {

  pet_list <- list(
    p1 = prior_PET("normal", list(0, 1))
  )
  expect_error(
    plot_prior_list(pet_list),
    "prior_list_mu"
  )

})


test_that("plot_prior_list rejects prior_list_mu when not needed", {

  simple_list <- list(
    p1 = prior("normal", list(0, 1))
  )
  expect_error(
    plot_prior_list(simple_list, prior_list_mu = list(prior("spike", list(0)))),
    "prior_list_mu"
  )

})

.test_plot_models_fixture <- function(is_null = c(FALSE, FALSE)) {
  make_summary <- function(parameter, mean, lo, hi) {
    out <- data.frame(
      Mean = mean,
      "0.025" = lo,
      "0.975" = hi,
      check.names = FALSE
    )
    class(out) <- c("BayesTools_table", "BayesTools_runjags_summary", "data.frame")
    attr(out, "parameters") <- parameter
    out
  }

  make_model <- function(mean, lo, hi, prior_prob, post_prob, BF, prior) {
    fit <- list()
    attr(fit, "prior_list") <- list(theta = prior)
    list(
      fit = fit,
      fit_summary = make_summary("theta", mean, lo, hi),
      inference = list(prior_prob = prior_prob, post_prob = post_prob, inclusion_BF = BF)
    )
  }

  models <- list(
    make_model(1, 0, 2, .25, .30, 1.2, prior("point", list(location = 0))),
    make_model(2, 1, 3, .75, .70, 0.8, prior("normal", list(mean = 0, sd = 1)))
  )
  theta <- structure(
    c(0, 1, 2, 3),
    class = "mixed_posteriors",
    prior_list = lapply(models, function(model) attr(model$fit, "prior_list")$theta)
  )
  samples <- list(theta = theta)
  inference <- list(theta = structure(list(), is_null = is_null))

  list(
    models = models,
    models_summary = lapply(models, `[[`, "fit_summary"),
    models_inference = lapply(models, `[[`, "inference"),
    prior_list = attr(theta, "prior_list"),
    samples = samples,
    theta = theta,
    inference = inference,
    total_inference = inference$theta
  )
}

test_that("plot_models rejects invalid ordering options before rendering", {

  fixture <- .test_plot_models_fixture()

  expect_error(
    plot_models(
      fixture$models,
      fixture$samples,
      fixture$inference,
      "theta",
      plot_type = "ggplot",
      order = list("sideways", "model")
    ),
    "order\\[\\[1\\]\\]"
  )
  expect_error(
    plot_models(
      fixture$models,
      fixture$samples,
      fixture$inference,
      "theta",
      plot_type = "ggplot",
      order = list("decreasing", "weight")
    ),
    "order\\[\\[2\\]\\]"
  )

})

test_that("plot_models helper rejects unknown non-point parameters", {

  fixture <- .test_plot_models_fixture()

  expect_error(
    BayesTools:::.plot_models_data_posterior(
      models_summary = fixture$models_summary,
      parameter = "missing",
      prior_list = fixture$prior_list,
      models_inference = fixture$models_inference
    ),
    "Posterior distribution summary for 'missing' is not available.",
    fixed = TRUE
  )

})

test_that("plot_models conditional view rejects empty model sets", {

  fixture <- .test_plot_models_fixture(is_null = c(TRUE, TRUE))

  expect_error(
    BayesTools:::.plot_models.simple(
      models_summary = fixture$models_summary,
      models_inference = fixture$models_inference,
      total_inference = fixture$total_inference,
      total_samples = fixture$theta,
      prior_list = fixture$prior_list,
      parameter = "theta",
      par_name = NULL,
      plot_type = "ggplot",
      prior = TRUE,
      conditional = TRUE,
      order = NULL,
      transformation = NULL,
      transformation_arguments = NULL,
      transformation_settings = FALSE
    )
  )

})

test_that("plot_models rejects bad transformations and invalid style switches", {

  fixture <- .test_plot_models_fixture()

  expect_error(
    BayesTools:::.plot_models.simple(
      models_summary = fixture$models_summary,
      models_inference = fixture$models_inference,
      total_inference = fixture$total_inference,
      total_samples = fixture$theta,
      prior_list = fixture$prior_list,
      parameter = "theta",
      par_name = NULL,
      plot_type = "ggplot",
      prior = TRUE,
      conditional = FALSE,
      order = NULL,
      transformation = list(fun = identity, inv = identity),
      transformation_arguments = NULL,
      transformation_settings = FALSE
    ),
    "Transformation must be either",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.plot_models.simple(
      models_summary = fixture$models_summary,
      models_inference = fixture$models_inference,
      total_inference = fixture$total_inference,
      total_samples = fixture$theta,
      prior_list = fixture$prior_list,
      parameter = "theta",
      par_name = NULL,
      plot_type = "ggplot",
      prior = TRUE,
      conditional = FALSE,
      order = NULL,
      transformation = NULL,
      transformation_arguments = NULL,
      transformation_settings = FALSE,
      show_updating = "yes"
    ),
    "invalid 'y' type",
    fixed = TRUE
  )

})

test_that("plot_models ggplot layer data follows requested estimate ordering", {

  skip_if_not_installed("ggplot2")

  make_summary <- function(mean, lo, hi) {
    out <- data.frame(
      Mean = mean,
      "0.025" = lo,
      "0.975" = hi,
      check.names = FALSE
    )
    class(out) <- c("BayesTools_table", "BayesTools_runjags_summary", "data.frame")
    attr(out, "parameters") <- "theta"
    out
  }

  make_model <- function(mean, lo, hi, prior_prob, post_prob, BF, prior_location) {
    fit <- list()
    attr(fit, "prior_list") <- list(theta = prior("point", list(location = prior_location)))
    list(
      fit = fit,
      fit_summary = make_summary(mean, lo, hi),
      inference = list(prior_prob = prior_prob, post_prob = post_prob, inclusion_BF = BF)
    )
  }

  models <- list(
    make_model(2, 1, 3, .2, .7, 4.0, 0),
    make_model(-1, -2, 0, .8, .3, 0.5, 1)
  )

  theta <- structure(
    c(-2, 0, 1, 3),
    class = "mixed_posteriors",
    prior_list = lapply(models, function(model) attr(model$fit, "prior_list")$theta)
  )
  samples <- list(theta = theta)
  inference <- list(theta = structure(list(), is_null = c(FALSE, FALSE)))

  plot <- plot_models(
    models,
    samples,
    inference,
    "theta",
    plot_type = "ggplot",
    prior = TRUE,
    order = list("decreasing", "estimate")
  )

  built <- ggplot2::ggplot_build(plot)
  posterior_points <- built$data[[2]]
  prior_points <- built$data[[4]]
  overall_polygon <- built$data[[5]]

  expect_equal(posterior_points$x, c(2, -1), tolerance = 1e-12)
  expect_equal(posterior_points$y, c(3, 4), tolerance = 1e-12)
  expect_gt(posterior_points$size[1], posterior_points$size[2])

  expect_equal(prior_points$x, c(0, 1), tolerance = 1e-12)
  expect_equal(prior_points$y, c(3.25, 4.25), tolerance = 1e-12)
  expect_gt(prior_points$size[2], prior_points$size[1])

  expect_equal(
    overall_polygon$x,
    c(
      unname(stats::quantile(theta, .025)),
      mean(theta),
      unname(stats::quantile(theta, .975)),
      mean(theta)
    ),
    tolerance = 1e-12
  )
})

test_that("plot_models conditional view filters nulls before probability ordering", {

  skip_if_not_installed("ggplot2")

  make_summary <- function(parameter, mean, lo, hi) {
    out <- data.frame(
      Mean = mean,
      "0.025" = lo,
      "0.975" = hi,
      check.names = FALSE
    )
    attr(out, "parameters") <- parameter
    out
  }

  models_summary <- list(
    make_summary("theta", 1, 0, 2),
    make_summary("theta", 20, 19, 21),
    make_summary("theta", 3, 2, 4),
    make_summary("theta", 4, 3, 5)
  )
  models_inference <- list(
    list(prior_prob = .25, post_prob = .10, inclusion_BF = 0.4),
    list(prior_prob = .25, post_prob = .99, inclusion_BF = 99),
    list(prior_prob = .25, post_prob = .20, inclusion_BF = 0.8),
    list(prior_prob = .25, post_prob = .80, inclusion_BF = 3.2)
  )
  prior_list <- list(
    prior("point", list(location = 1)),
    prior("point", list(location = 20)),
    prior("point", list(location = 3)),
    prior("point", list(location = 4))
  )
  total_inference <- structure(list(), is_null = c(FALSE, TRUE, FALSE, FALSE))

  plot <- BayesTools:::.plot_models.simple(
    models_summary = models_summary,
    models_inference = models_inference,
    total_inference = total_inference,
    total_samples = c(1, 3, 4, 4),
    prior_list = prior_list,
    parameter = "theta",
    par_name = NULL,
    plot_type = "ggplot",
    prior = TRUE,
    conditional = TRUE,
    order = list("decreasing", "probability"),
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  built <- ggplot2::ggplot_build(plot)
  posterior_points <- built$data[[2]]
  prior_points <- built$data[[4]]

  expect_equal(posterior_points$x, c(4, 3, 1), tolerance = 1e-12)
  expect_false(any(posterior_points$x == 20))
  expect_equal(posterior_points$y, c(3, 4, 5), tolerance = 1e-12)
  expect_gt(posterior_points$size[1], posterior_points$size[2])
  expect_gt(posterior_points$size[2], posterior_points$size[3])

  expect_equal(prior_points$x, c(4, 3, 1), tolerance = 1e-12)
  expect_false(any(prior_points$x == 20))
  expect_equal(prior_points$y, c(3.25, 4.25, 5.25), tolerance = 1e-12)
})

test_that("posterior plot data separates spike mass from continuous samples", {
  theta <- c(rep(0, 3), seq_len(7))
  attr(theta, "models_ind") <- c(rep(1, 3), rep(2, 7))
  attr(theta, "prior_list") <- list(
    prior("point", list(location = 0)),
    prior("normal", list(mean = 0, sd = 1))
  )

  plot_data <- BayesTools:::.plot_data_samples.simple(
    samples = list(theta = theta),
    parameter = "theta",
    n_points = 128,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(names(plot_data), c("density", "points1"))
  expect_equal(plot_data$points1$x, 0)
  expect_equal(plot_data$points1$y, .3)
  expect_equal(plot_data$density$samples, seq_len(7))

  dx <- plot_data$density$x[2] - plot_data$density$x[1]
  expect_equal(sum(plot_data$density$y) * dx, .7, tolerance = .08)
})

test_that("bounded posterior KDE reflects support and keeps spike mass separate", {
  theta <- c(rep(0, 20), seq(.005, .995, length.out = 80))
  attr(theta, "models_ind") <- c(rep(1, 20), rep(2, 80))
  attr(theta, "prior_list") <- list(
    prior("point", list(location = 0)),
    prior("uniform", list(0, 1))
  )

  plot_data <- BayesTools:::.plot_data_samples.simple(
    samples = list(theta = theta),
    parameter = "theta",
    n_points = 512,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(plot_data$points1$x, 0)
  expect_equal(plot_data$points1$y, .2)
  expect_false(any(plot_data$density$samples == 0))
  expect_equal(plot_data$density$samples, theta[21:100])
  expect_equal(range(plot_data$density$x), c(0, 1))
  expect_equal(plot_data$density$x[1:2], c(0, 0), tolerance = 1e-12)
  expect_equal(plot_data$density$y[1], 0)
  expect_gt(plot_data$density$y[2], 0)
  expect_equal(tail(plot_data$density$x, 2), c(1, 1), tolerance = 1e-12)
  expect_gt(tail(plot_data$density$y, 2)[1], 0)
  expect_equal(tail(plot_data$density$y, 1), 0)
  expect_equal(.model_plot_density_mass(plot_data$density), .8, tolerance = .03)
  expect_true(isTRUE(attr(plot_data$density, "boundary_reflection")))
})

test_that("conditional posterior prior overlay uses conditioned spike-and-slab slab", {
  prior_list <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(mean = 1, sd = 0.2)),
      prior_inclusion = prior("point", list(location = 0.5))
    )
  )

  posterior <- cbind(
    theta           = c(rep(0, 50), seq(0.5, 1.5, length.out = 50)),
    theta_indicator = c(rep(0, 50), rep(1, 50))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list

  samples <- as_mixed_posteriors(
    model       = fit,
    parameters  = "theta",
    conditional = "theta"
  )
  prior_list_plot <- BayesTools:::.simplify_prior_list(attr(samples$theta, "prior_list"))

  expect_true(BayesTools:::.plot_data_prior_should_use_context(
    samples = samples,
    parameter = "theta",
    transform_scaled = FALSE,
    prior_list = prior_list_plot
  ))

  plot_data_prior <- BayesTools:::.plot_data_prior_density_context(
    prior_density_context = attr(samples, "prior_density_context"),
    samples = samples,
    parameter = "theta",
    prior_list = prior_list_plot,
    n_points = 128,
    x_range = NULL,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(names(plot_data_prior), "density")
  expect_false(any(vapply(plot_data_prior, inherits, logical(1), what = "density.prior.point")))
  expect_gt(max(plot_data_prior$density$y), 0)

  expect_s3_class(
    plot_posterior(samples, "theta", plot_type = "ggplot", prior = TRUE),
    "ggplot"
  )
})

test_that("plot_posterior handles attached point priors outside xlim", {

  theta <- structure(
    rep(0, 64),
    class = c("mixed_posteriors.simple", "mixed_posteriors"),
    prior_list = list(prior("spike", list(0))),
    models_ind = rep(1L, 64),
    prior_density = BayesTools:::.prior_linear_density_point(0)
  )
  samples <- list(theta = theta)

  expect_null(BayesTools:::.plot_data_attached_prior_density(
    samples = samples,
    parameter = "theta",
    n_points = 32,
    x_range = c(10, 20)
  ))
  expect_s3_class(
    plot_posterior(
      samples,
      "theta",
      prior = TRUE,
      plot_type = "ggplot",
      xlim = c(10, 20),
      n_points = 32
    ),
    "ggplot"
  )
})

test_that("conditional PET/PEESE prior overlays do not reintroduce excluded bias branches", {
  prior_list <- list(
    bias = prior_mixture(
      list(
        prior_none(prior_weights = 1),
        prior_PET("normal", list(mean = 0, sd = 1), prior_weights = 1),
        prior_PEESE("normal", list(mean = 0, sd = 1), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE, FALSE)
    )
  )

  posterior <- cbind(
    bias_indicator = c(rep(1, 20), rep(2, 20), rep(3, 20)),
    PET            = c(rep(0, 20), seq(0.1, 1.0, length.out = 20), rep(0, 20)),
    PEESE          = c(rep(0, 20), rep(0, 20), seq(0.2, 1.2, length.out = 20))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list

  samples <- as_mixed_posteriors(
    model       = fit,
    parameters  = "bias",
    conditional = "PEESE",
    force_plots = TRUE
  )
  simplified <- BayesTools:::.simplify_as_mixed_posterior_bias(samples, "PEESE")
  prior_list_plot <- BayesTools:::.simplify_prior_list(attr(simplified$PEESE, "prior_list"))

  point_prior <- vapply(prior_list_plot, is.prior.point, logical(1))
  expect_true(all(vapply(prior_list_plot[point_prior], BayesTools:::.prior_model_weight, numeric(1)) == 0))

  plot_data <- BayesTools:::.plot_data_samples.simple(
    samples = simplified,
    parameter = "PEESE",
    n_points = 128,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )
  expect_equal(names(plot_data), "density")
  expect_equal(plot_data$density$samples, posterior[posterior[, "bias_indicator"] == 3, "PEESE"])

  plot_data_prior <- BayesTools:::.plot_data_prior_list.simple(
    prior_list = prior_list_plot,
    x_seq = NULL,
    x_range = NULL,
    x_range_quant = NULL,
    n_points = 128,
    n_samples = 10000,
    force_samples = FALSE,
    individual = FALSE,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(names(plot_data_prior), "density")
  expect_false(any(vapply(plot_data_prior, inherits, logical(1), what = "density.prior.point")))

  expect_s3_class(
    plot_posterior(samples, "PEESE", plot_type = "ggplot", prior = TRUE, individual = TRUE),
    "ggplot"
  )
})

test_that("omega posterior KDE excludes exact one samples and reflects bounded support", {
  continuous_samples <- seq(.005, .995, length.out = 75)
  omega_samples <- cbind(
    "omega[0,0.05]" = rep(1, 100),
    "omega[0.05,1]" = c(rep(1, 25), continuous_samples)
  )
  attr(omega_samples, "prior_list") <- list(
    prior_weightfunction(
      "one-sided",
      c(.05),
      wf_independent(prior("beta", list(1, 1)))
    )
  )
  attr(omega_samples, "models_ind") <- rep(1, nrow(omega_samples))

  plot_data <- BayesTools:::.plot_data_samples.weightparameter(
    samples = list(omega = omega_samples),
    parameter = "omega[0.05,1]",
    n_points = 512
  )

  expect_equal(plot_data$points1$x, 1)
  expect_equal(plot_data$points1$y, .25)
  expect_false(any(plot_data$density$samples == 1))
  expect_equal(plot_data$density$samples, continuous_samples)
  expect_equal(range(plot_data$density$x), c(0, 1))
  expect_equal(plot_data$density$x[1:2], c(0, 0), tolerance = 1e-12)
  expect_equal(plot_data$density$y[1], 0)
  expect_gt(plot_data$density$y[2], 0)
  expect_equal(tail(plot_data$density$x, 2), c(1, 1), tolerance = 1e-12)
  expect_gt(tail(plot_data$density$y, 2)[1], 0)
  expect_equal(tail(plot_data$density$y, 1), 0)
  expect_equal(.model_plot_density_mass(plot_data$density), .75, tolerance = .03)
  expect_true(isTRUE(attr(plot_data$density, "boundary_reflection")))
})

test_that("omega prior and posterior curves integrate their continuous masses", {
  prior_list <- list(
    bias = prior_mixture(
      list(
        prior_none(prior_weights = 1),
        prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)), prior_weights = 1 / 3),
        prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 1, 1)), prior_weights = 1 / 3),
        prior_PET("normal", list(0, 1), prior_weights = 1 / 3)
      ),
      is_null = c(TRUE, FALSE, FALSE, FALSE)
    )
  )
  parameter <- "omega[0.05,0.975]"

  prior_data <- BayesTools:::.plot_data_prior_list.weightparameter(
    prior_list$bias,
    parameter = parameter,
    n_points = 2048,
    n_samples = 1
  )
  prior_point_mass <- sum(vapply(prior_data, function(component) {
    if(inherits(component, "density.prior.point")) component$y else 0
  }, numeric(1)))

  n_samples <- 4000
  models_ind <- c(rep(1, 500), rep(2, 1500), rep(3, 1500), rep(4, 500))
  omega_samples <- matrix(
    1,
    nrow = n_samples,
    ncol = 4,
    dimnames = list(NULL, c(
      "omega[0,0.025]",
      "omega[0.025,0.05]",
      "omega[0.05,0.975]",
      "omega[0.975,1]"
    ))
  )
  omega_samples[models_ind == 2, parameter] <- seq(.005, .995, length.out = sum(models_ind == 2))
  omega_samples[models_ind == 3, parameter] <- seq(.995, .005, length.out = sum(models_ind == 3))
  attr(omega_samples, "prior_list") <- prior_list$bias
  attr(omega_samples, "models_ind") <- models_ind

  posterior_data <- BayesTools:::.plot_data_samples.weightparameter(
    samples = list(bias = omega_samples),
    parameter = parameter,
    n_points = 512
  )

  expect_equal(.model_plot_density_mass(prior_data$density), 1 / 3, tolerance = 1e-3)
  expect_equal(prior_point_mass, 2 / 3, tolerance = 1e-8)
  expect_equal(.model_plot_density_mass(posterior_data$density), .75, tolerance = .03)
  expect_equal(posterior_data$points1$y, .25, tolerance = 1e-8)
  expect_equal(.model_plot_density_mass(prior_data$density) + prior_point_mass, 1, tolerance = 1e-3)
  expect_equal(.model_plot_density_mass(posterior_data$density) + posterior_data$points1$y, 1, tolerance = .03)
  expect_gt(.model_plot_density_mass(posterior_data$density), .model_plot_density_mass(prior_data$density))
})

test_that("posterior plot data uses stored posterior density when available", {
  theta <- seq(-2, 2, length.out = 40)
  stored_x <- seq(-3, 3, length.out = 61)
  stored_y <- stats::dnorm(stored_x, mean = .25, sd = .9)
  attr(theta, "models_ind") <- rep(1, length(theta))
  attr(theta, "prior_list") <- list(prior("normal", list(mean = 0, sd = 1)))
  attr(theta, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

  plot_data <- BayesTools:::.plot_data_samples.simple(
    samples = list(theta = theta),
    parameter = "theta",
    n_points = 16,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE,
    density_method = "precomputed"
  )

  expect_equal(plot_data$density$x, stored_x)
  expect_equal(plot_data$density$y, stored_y)
  expect_equal(attr(plot_data$density, "posterior_density_method"), "iwmde")
})

test_that("posterior plot data does not add sample spikes to stored full density", {
  theta <- c(rep(0, 25), seq(-2, 2, length.out = 75))
  stored_x <- seq(-2, 2, length.out = 51)
  stored_y <- stats::dnorm(stored_x)
  attr(theta, "models_ind") <- c(rep(1, 25), rep(2, 75))
  attr(theta, "prior_list") <- list(
    prior("point", list(location = 0)),
    prior("normal", list(mean = 0, sd = 1))
  )
  attr(theta, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

  expect_warning(
    plot_data <- BayesTools:::.plot_data_samples.simple(
      samples                  = list(theta = theta),
      parameter                = "theta",
      n_points                 = 16,
      transformation           = NULL,
      transformation_arguments = NULL,
      transformation_settings  = FALSE,
      density_method           = "precomputed"
    ),
    "does not declare 'point_masses'",
    fixed = TRUE
  )

  expect_equal(plot_data$density$x, stored_x)
  expect_equal(
    length(plot_data[vapply(plot_data, inherits, logical(1), "density.prior.point")]),
    0L
  )
})

test_that("posterior plot data ignores stored density by default", {
  theta <- seq(-2, 2, length.out = 40)
  stored_x <- seq(-3, 3, length.out = 61)
  stored_y <- stats::dnorm(stored_x, mean = .25, sd = .9)
  attr(theta, "models_ind") <- rep(1, length(theta))
  attr(theta, "prior_list") <- list(prior("normal", list(mean = 0, sd = 1)))
  attr(theta, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

  plot_data <- BayesTools:::.plot_data_samples.simple(
    samples = list(theta = theta),
    parameter = "theta",
    n_points = 16,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(length(plot_data$density$x), 16)
  expect_false(identical(plot_data$density$x, stored_x))
  expect_null(attr(plot_data$density, "posterior_density_method"))
})

test_that("plot-model posterior data falls back to point priors for absent rows", {
  make_summary <- function(parameters, mean, lo, hi) {
    out <- data.frame(
      Mean = mean,
      "0.025" = lo,
      "0.975" = hi,
      check.names = FALSE
    )
    attr(out, "parameters") <- parameters
    out
  }

  models_summary <- list(
    make_summary("other", 99, 98, 100),
    make_summary("theta", 2, 1, 3)
  )
  prior_list <- list(
    prior("point", list(location = 0)),
    prior("normal", list(mean = 0, sd = 1))
  )
  models_inference <- list(
    list(prior_prob = .2, post_prob = .7, inclusion_BF = 4),
    list(prior_prob = .8, post_prob = .3, inclusion_BF = .5)
  )

  plot_data <- BayesTools:::.plot_models_data_posterior(
    models_summary = models_summary,
    parameter = "theta",
    prior_list = prior_list,
    models_inference = models_inference
  )

  expect_equal(plot_data$model, 1:2)
  expect_equal(plot_data$y, c(0, 2))
  expect_equal(plot_data$y_lCI, c(0, 1))
  expect_equal(plot_data$y_uCI, c(0, 3))
  expect_equal(plot_data$prior_prob, c(.2, .8))
  expect_equal(plot_data$post_prob, c(.7, .3))
  expect_equal(plot_data$BF, c(4, .5))
})

test_that("PET-PEESE prior plot data uses deterministic linear-combination summaries", {

  prior_list <- list(
    PET = prior_PET("normal", list(0, 1), truncation = list(-Inf, Inf))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )

  plot_data <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = c(0, 0.5, 1),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 3,
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )

  expect_null(plot_data$samples)
  expect_equal(plot_data$y, c(0, 0, 0), tolerance = 1e-8)
  expect_equal(plot_data$y_lCI, stats::qnorm(.025) * c(0, 0.5, 1), tolerance = 0.02)
  expect_equal(plot_data$y_uCI, stats::qnorm(.975) * c(0, 0.5, 1), tolerance = 0.02)
  expect_equal(
    attr(plot_data, "y_range"),
    range(plot_data$y, plot_data$y_lCI, plot_data$y_uCI)
  )
})

test_that("PET-PEESE public prior plots include uncertainty ribbons in y limits", {

  pet_prior <- prior_PET(
    "normal",
    list(0, 1),
    truncation = list(-Inf, Inf)
  )
  mu_prior <- prior("spike", list(0))

  list_plot <- plot_prior_list(
    prior_list = list(pet_prior),
    prior_list_mu = list(mu_prior),
    plot_type = "ggplot",
    xlim = c(0, 1),
    n_points = 3
  )
  list_ribbon <- list_plot$layers[[1]]$data$y
  list_limits <- list_plot$scales$get_scales("y")$limits

  expect_lte(min(list_limits), min(list_ribbon))
  expect_gte(max(list_limits), max(list_ribbon))

  single_plot <- plot(
    pet_prior,
    plot_type = "ggplot",
    xlim = c(0, 1),
    n_points = 3
  )
  single_ribbon <- single_plot$layers[[1]]$data$y
  single_limits <- single_plot$scales$get_scales("y")$limits

  expect_lte(min(single_limits), min(single_ribbon))
  expect_gte(max(single_limits), max(single_ribbon))
})

test_that("PET-PEESE public prior plots order ribbons after decreasing transformations", {

  x_seq <- c(0, 0.5, 1)
  bias_priors <- list(
    PET = prior_PET(
      "normal",
      list(0, 1),
      truncation = list(-Inf, Inf)
    ),
    PEESE = prior_PEESE(
      "normal",
      list(0, 1),
      truncation = list(-Inf, Inf)
    )
  )
  se_scale <- list(PET = x_seq, PEESE = x_seq^2)

  for(prior_name in names(bias_priors)){
    prior_plot <- plot(
      bias_priors[[prior_name]],
      plot_type = "ggplot",
      x_seq = x_seq,
      transformation = "lin",
      transformation_arguments = list(a = 2, b = -3)
    )

    ribbon <- prior_plot$layers[[1]]$data
    n_points <- length(x_seq)
    lower <- ribbon$y[seq_len(n_points)]
    upper <- rev(tail(ribbon$y, n_points))

    expect_true(all(lower <= upper))
    expect_equal(
      lower,
      (2 - 3 * stats::qnorm(.975)) * se_scale[[prior_name]]
    )
    expect_equal(
      upper,
      (2 - 3 * stats::qnorm(.025)) * se_scale[[prior_name]]
    )
  }
})

test_that("plot_prior_list keeps retained components paired with their weights", {

  point_priors <- list(
    prior("point", list(100), prior_weights = 1),
    prior("point", list(0), prior_weights = 9999)
  )
  point_plot <- plot_prior_list(
    point_priors,
    plot_type = "ggplot",
    n_samples = 10000
  )
  point_data <- point_plot$layers[[1]]$data

  expect_equal(point_data$x, 0)
  expect_equal(point_data$yend, .9999)

  continuous_priors <- list(
    prior("normal", list(100, 1), prior_weights = 1),
    prior("normal", list(0, 1), prior_weights = 9999)
  )
  density_plot <- plot_prior_list(
    continuous_priors,
    plot_type = "ggplot",
    xlim = c(-5, 105),
    n_points = 1101,
    n_samples = 10000
  )
  density_data <- density_plot$layers[[1]]$data

  expect_equal(density_data$x[which.max(density_data$y)], 0)
  expect_equal(max(density_data$y), stats::dnorm(0) * .9999, tolerance = 1e-10)
  expect_equal(.model_plot_density_mass(density_data), .9999, tolerance = 1e-5)
})

test_that("factor ggplot prior point layers use point plot data", {

  prior_list <- list(
    beta = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE)
    )
  )
  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = prior_list,
    weights    = c(beta = 1),
    n_grid     = 256
  )
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    density,
    n_points   = 64,
    factor     = TRUE,
    level      = 1,
    level_name = "A"
  )
  plot <- BayesTools:::.plot_prior_list.factor(
    plot_data = plot_data,
    plot_type = "ggplot",
    hardcode  = TRUE,
    legend    = FALSE
  )

  segment_layers <- vapply(plot[["layers"]], function(layer) {
    inherits(layer[["geom"]], "GeomSegment")
  }, logical(1))
  segment_data <- plot[["layers"]][[which(segment_layers)]][["data"]]

  expect_equal(sum(segment_layers), 1)
  expect_equal(NROW(segment_data), 1)
  expect_equal(segment_data[["x"]], 0, tolerance = 1e-8)
})

test_that("factor plot data normalization preserves unnamed points", {

  density_data <- BayesTools:::.prior_linear_density_to_plot_data(
    BayesTools:::.prior_linear_density_point(1),
    n_points   = 64,
    factor     = TRUE,
    level      = 1,
    level_name = "A"
  )[["points1"]]
  class(density_data) <- setdiff(class(density_data), "density.prior.point")

  point_data <- list(
    call    = call("density", "point"),
    bw      = NULL,
    n       = 64,
    x       = 0,
    y       = 0.5,
    samples = NULL
  )
  class(point_data) <- c("density", "density.prior", "density.prior.point")
  attr(point_data, "x_range") <- 0
  attr(point_data, "y_range") <- c(0, 0.5)

  normalized <- BayesTools:::.plot_prior_factor_normalize_data(
    list(density = density_data, point = point_data)
  )

  expect_equal(normalized[["level_names"]], "A")
  expect_equal(attr(normalized[["densities"]][[1]], "level_id"), 1)
  expect_equal(attr(normalized[["points"]][[1]], "level_id"), NA_integer_)
})

test_that("factor point layers use matching level colors when metadata is present", {

  plot_data <- c(
    BayesTools:::.prior_linear_density_to_plot_data(
      BayesTools:::.prior_linear_density_point(0),
      n_points   = 64,
      factor     = TRUE,
      level      = 1,
      level_name = "A"
    ),
    BayesTools:::.prior_linear_density_to_plot_data(
      BayesTools:::.prior_linear_density_point(1),
      n_points   = 64,
      factor     = TRUE,
      level      = 2,
      level_name = "B"
    )
  )
  plot <- BayesTools:::.plot_prior_list.factor(
    plot_data = plot_data,
    plot_type = "ggplot",
    col       = c("red", "blue"),
    hardcode  = TRUE
  )

  segment_layers <- vapply(plot[["layers"]], function(layer) {
    inherits(layer[["geom"]], "GeomSegment")
  }, logical(1))
  segments <- plot[["layers"]][segment_layers]

  expect_equal(length(segments), 2L)
  expect_equal(
    unname(vapply(segments, function(layer) layer[["aes_params"]][["colour"]], character(1))),
    c("red", "blue")
  )
  expect_equal(
    unname(vapply(segments, function(layer) layer[["data"]][["x"]], numeric(1))),
    c(0, 1),
    tolerance = 1e-8
  )
})

test_that("factor point layers reuse single style values for all levels", {

  plot_data <- c(
    BayesTools:::.prior_linear_density_to_plot_data(
      BayesTools:::.prior_linear_density_point(0),
      n_points   = 64,
      factor     = TRUE,
      level      = 1,
      level_name = "A"
    ),
    BayesTools:::.prior_linear_density_to_plot_data(
      BayesTools:::.prior_linear_density_point(1),
      n_points   = 64,
      factor     = TRUE,
      level      = 2,
      level_name = "B"
    )
  )
  plot <- BayesTools:::.plot_prior_list.factor(
    plot_data = plot_data,
    plot_type = "ggplot",
    col       = "red",
    lty       = 2,
    hardcode  = TRUE
  )

  segment_layers <- vapply(plot[["layers"]], function(layer) {
    inherits(layer[["geom"]], "GeomSegment")
  }, logical(1))
  segments <- plot[["layers"]][segment_layers]

  expect_equal(length(segments), 2L)
  expect_equal(
    unname(vapply(segments, function(layer) layer[["aes_params"]][["colour"]], character(1))),
    c("red", "red")
  )
  expect_equal(
    unname(vapply(segments, function(layer) layer[["aes_params"]][["linetype"]], numeric(1))),
    c(2, 2)
  )
})

test_that("factor plots map component-expanded prior styles back to levels", {

  prior_list <- list(
    beta = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE)
    )
  )
  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = prior_list,
    weights    = c(beta = 1),
    n_grid     = 256
  )
  plot_data <- c(
    BayesTools:::.prior_linear_density_to_plot_data(
      density,
      n_points   = 64,
      factor     = TRUE,
      level      = 1,
      level_name = "A"
    ),
    BayesTools:::.prior_linear_density_to_plot_data(
      density,
      n_points   = 64,
      factor     = TRUE,
      level      = 2,
      level_name = "B"
    )
  )
  plot <- BayesTools:::.plot_prior_list.factor(
    plot_data = plot_data,
    plot_type = "ggplot",
    col       = c("red", "red", "blue", "blue"),
    lty       = c(2, 2, 3, 3),
    hardcode  = TRUE
  )

  segment_layers <- vapply(plot[["layers"]], function(layer) {
    inherits(layer[["geom"]], "GeomSegment")
  }, logical(1))
  line_layers <- vapply(plot[["layers"]], function(layer) {
    inherits(layer[["geom"]], "GeomLine")
  }, logical(1))
  segments <- plot[["layers"]][segment_layers]
  lines    <- plot[["layers"]][line_layers]

  expect_equal(
    unname(vapply(segments, function(layer) layer[["aes_params"]][["colour"]], character(1))),
    c("red", "blue")
  )
  expect_equal(
    unname(vapply(segments, function(layer) layer[["aes_params"]][["linetype"]], numeric(1))),
    c(2, 3)
  )
  expect_equal(
    unname(vapply(lines, function(layer) layer[["aes_params"]][["colour"]], character(1))),
    c("red", "blue")
  )
  expect_equal(
    unname(vapply(lines, function(layer) layer[["aes_params"]][["linetype"]], numeric(1))),
    c(2, 3)
  )
})

.capture_graphics_legend_calls_for_test <- function(code){

  code <- substitute(code)
  env  <- parent.frame()
  calls <- list()

  testthat::local_mocked_bindings(
    legend = function(...){
      calls[[length(calls) + 1L]] <<- list(...)
      invisible(NULL)
    },
    .package = "graphics"
  )

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  eval(code, env)
  calls
}

.factor_plot_data_for_legend_test <- function(){

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(0, 1))),
    weights    = c(theta = 1),
    n_grid     = 256
  )

  c(
    BayesTools:::.prior_linear_density_to_plot_data(
      density,
      n_points   = 64,
      factor     = TRUE,
      level      = 1,
      level_name = "theta[dif: A]"
    ),
    BayesTools:::.prior_linear_density_to_plot_data(
      density,
      n_points   = 64,
      factor     = TRUE,
      level      = 2,
      level_name = "theta[dif: B]"
    )
  )
}

.factor_marginal_samples_for_legend_test <- function(){

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(0, 1))),
    weights    = c(theta = 1),
    n_grid     = 256
  )

  set.seed(917)
  level_1 <- stats::rnorm(80, -0.25, 0.4)
  level_2 <- stats::rnorm(80,  0.35, 0.4)
  attr(level_1, "prior_density") <- density
  attr(level_2, "prior_density") <- density

  theta <- list(A = level_1, B = level_2)
  class(theta) <- c("marginal_posterior.factor", "marginal_posterior")

  list(theta = theta)
}

test_that("factor legend suppression does not collapse level styles", {

  skip_if_not_installed("ggplot2")

  plot <- BayesTools:::.plot_prior_list.factor(
    plot_data = .factor_plot_data_for_legend_test(),
    plot_type = "ggplot",
    legend    = FALSE
  )
  built <- ggplot2::ggplot_build(plot)

  line_layers <- vapply(plot[["layers"]], function(layer) {
    inherits(layer[["geom"]], "GeomLine")
  }, logical(1))
  colours <- unique(unlist(lapply(built[["data"]], function(layer_data){
    if("colour" %in% names(layer_data)) layer_data[["colour"]]
  })))

  expect_true(all(vapply(plot[["layers"]][line_layers], function(layer) identical(layer[["show.legend"]], FALSE), logical(1))))
  expect_gte(length(colours), 2L)
})

test_that("standalone factor prior helper draws one owned base legend", {

  calls <- .capture_graphics_legend_calls_for_test(
    BayesTools:::.plot_prior_list.factor(
      plot_data = .factor_plot_data_for_legend_test(),
      plot_type = "base"
    )
  )

  expect_length(calls, 1L)
})

test_that("plot_marginal owns factor legends in base overlays", {

  samples <- .factor_marginal_samples_for_legend_test()

  expect_length(
    .capture_graphics_legend_calls_for_test(
      plot_marginal(samples, "theta", prior = TRUE, n_points = 64)
    ),
    1L
  )
  expect_length(
    .capture_graphics_legend_calls_for_test(
      plot_marginal(samples, "theta", prior = TRUE, n_points = 64, legend = FALSE, dots_prior = list(legend = TRUE))
    ),
    0L
  )
  expect_length(
    .capture_graphics_legend_calls_for_test(
      plot_marginal(samples, "theta", prior = FALSE, n_points = 64)
    ),
    1L
  )

  calls <- .capture_graphics_legend_calls_for_test(
    plot_marginal(
      samples,
      "theta",
      prior = TRUE,
      n_points = 64,
      legend_title = "tailor",
      legend_labels = c("0", "1"),
      legend_position = "bottomleft"
    )
  )

  expect_length(calls, 1L)
  expect_equal(calls[[1L]][[1L]], "bottomleft")
  expect_equal(calls[[1L]][["legend"]], c("0", "1"))
  expect_equal(calls[[1L]][["title"]], "tailor")
})

test_that("plot_marginal factor legend metadata is used in ggplot scales", {

  skip_if_not_installed("ggplot2")

  plot <- plot_marginal(
    .factor_marginal_samples_for_legend_test(),
    "theta",
    prior = TRUE,
    plot_type = "ggplot",
    n_points = 64,
    legend_title = "tailor",
    legend_labels = c("0", "1"),
    legend_position = "bottom"
  )
  built <- ggplot2::ggplot_build(plot)
  colour_scale <- built[["plot"]][["scales"]]$get_scales("colour")

  expect_equal(colour_scale[["name"]], "tailor")
  expect_equal(colour_scale$get_labels(), c("0", "1"))
  expect_equal(built[["plot"]][["theme"]][["legend.position"]], "bottom")
})

test_that("omega plot helpers use selection components from composed bias priors", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking  <- prior_phacking(form = "linear")
  bias      <- prior_bias(selection = selection, phacking = phacking)

  model <- matrix(
    c(
      1.0, .5, .2, .010, 1,
      1.0, .5, .4, .020, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(model) <- c("omega[1]", "omega[2]", "alpha", "pi_null", "phack_kind")
  class(model) <- c("matrix", "BayesTools_fit")
  attr(model, "prior_list") <- list(bias = bias)

  mixed <- as_mixed_posteriors(model, parameters = "bias")

  plot_data <- BayesTools:::.plot_data_samples.weightfunction(
    mixed,
    x_seq         = NULL,
    x_range       = NULL,
    x_range_quant = NULL,
    n_points      = 32
  )
  expect_equal(plot_data$x, c(0, .025, .025, 1))
  expect_equal(unname(plot_data$y), c(1, 1, .5, .5), tolerance = 1e-8)

  omega_samples <- BayesTools:::.simplify_as_mixed_posterior_bias(mixed, "omega")
  prior_data <- BayesTools:::.plot_data_prior_list.weightparameter(
    attr(omega_samples$omega, "prior_list"),
    parameter = "omega[0.025,1]",
    n_points  = 32,
    n_samples = 100
  )
  expect_equal(prior_data$points1$x, .5, tolerance = 1e-8)
})

test_that("omega plot helpers preserve composed selection branches in bias mixtures", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking  <- prior_phacking(form = "linear")
  bias      <- prior_mixture(
    list(prior_none(), phacking, prior_bias(selection = selection, phacking = phacking)),
    is_null = c(TRUE, FALSE, FALSE)
  )

  model <- matrix(
    c(
      1, 1.0, 1.0, 0.0, 0.000, 0,
      2, 1.0, 1.0, 0.2, 0.010, 1,
      3, 1.0, 0.5, 0.4, 0.020, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(model) <- c("bias_indicator", "omega[1]", "omega[2]", "alpha", "pi_null", "phack_kind")
  class(model) <- c("matrix", "BayesTools_fit")
  attr(model, "prior_list") <- list(bias = bias)

  mixed <- as_mixed_posteriors(model, parameters = "bias")
  omega_samples <- BayesTools:::.simplify_as_mixed_posterior_bias(mixed, "omega")

  plot_data <- BayesTools:::.plot_data_samples.weightfunction(
    omega_samples,
    x_seq         = NULL,
    x_range       = NULL,
    x_range_quant = NULL,
    n_points      = 32
  )
  expect_equal(plot_data$x, c(0, .025, .025, 1))

  prior_data <- BayesTools:::.plot_data_prior_list.weightparameter(
    attr(omega_samples$omega, "prior_list"),
    parameter = "omega[0.025,1]",
    n_points  = 32,
    n_samples = 100
  )
  expect_equal(prior_data$points1$x, .5, tolerance = 1e-8)
})

.test_weightfunction_posterior_samples <- function(){

  weight_prior <- prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5)))
  omega_samples <- matrix(c(rep(1, 20), rep(.5, 20)), ncol = 2)
  colnames(omega_samples) <- c("omega[0,0.05]", "omega[0.05,1]")
  attr(omega_samples, "prior_list") <- list(weight_prior)
  attr(omega_samples, "models_ind") <- rep(1, nrow(omega_samples))
  class(omega_samples) <- c("matrix", "mixed_posteriors")

  samples <- list(omega = omega_samples)
  class(samples) <- c("as_mixed_posteriors", "mixed_posteriors")
  samples
}

test_that("plot_posterior draws observed p-values on ggplot weightfunction plots", {

  samples <- .test_weightfunction_posterior_samples()
  p <- c(.025, .05, .8)

  plot <- plot_posterior(
    samples,
    "weightfunction",
    plot_type = "ggplot",
    data = data.frame(p = p),
    show_data = TRUE,
    rescale_x = TRUE,
    dots_data = list(color = "red", alpha = .5, linewidth = .4, rug_side = "top", rug_height = .02)
  )

  rug_layers <- vapply(plot$layers, function(layer) inherits(layer$geom, "GeomRug"), logical(1))
  rug_layer <- plot$layers[[which(rug_layers)]]
  expected_p <- stats::approx(c(0, .05, 1), c(0, .5, 1), xout = p, rule = 2)$y

  expect_s3_class(plot, "ggplot")
  expect_equal(sum(rug_layers), 1L)
  expect_equal(rug_layer$data$p, expected_p, tolerance = 1e-8)
  expect_equal(rug_layer$geom_params$sides, "t")
})

test_that("plot_posterior draws observed p-values on base weightfunction plots", {

  samples <- .test_weightfunction_posterior_samples()

  device_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(device_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_posterior(
    samples,
    "weightfunction",
    plot_type = "base",
    data = c(.01, .2, .8),
    show_data = TRUE,
    rescale_x = FALSE,
    dots_data = list(col = "blue", lwd = 2, side = "top", height = .02)
  ))
})

test_that("plot_posterior validates observed p-value rug inputs", {

  samples <- .test_weightfunction_posterior_samples()

  expect_error(
    plot_posterior(samples, "weightfunction", show_data = TRUE),
    "'data' must be supplied when 'show_data = TRUE'.",
    fixed = TRUE
  )
  expect_error(
    plot_posterior(samples, "mu", show_data = TRUE, data = .1),
    "'show_data' is currently supported only for weightfunction posterior plots.",
    fixed = TRUE
  )
  expect_error(
    plot_posterior(samples, "weightfunction", show_data = TRUE, data = .1, individual = TRUE),
    "'show_data' is supported only for the full weightfunction posterior plot.",
    fixed = TRUE
  )
  expect_error(
    plot_posterior(samples, "weightfunction", show_data = TRUE, data = data.frame(q = .1)),
    "'data' must be a numeric vector of p-values or a data frame with a 'p' column.",
    fixed = TRUE
  )
  expect_error(
    plot_posterior(samples, "weightfunction", show_data = TRUE, data = c(.1, 1.2)),
    "The 'data' must be equal or lower than 1.",
    fixed = TRUE
  )
})

test_that("PET-PEESE prior plot data uses CDF quantiles for half-Cauchy PET", {

  x_seq <- c(0, 0.25, 0.5, 1)
  prior_list <- list(
    PET = prior_PET("cauchy", list(0, 1))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )

  positive <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu,
    effect_direction         = "positive"
  )

  expected_positive <- function(p) x_seq * tan(pi * p / 2)
  expect_null(positive$samples)
  expect_equal(positive$y_lCI, expected_positive(.025), tolerance = 1e-4)
  expect_equal(positive$y,     expected_positive(.500), tolerance = 1e-4)
  expect_equal(positive$y_uCI, expected_positive(.975), tolerance = 1e-4)

  negative <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu,
    effect_direction         = "negative"
  )

  expect_equal(negative$y_lCI, -positive$y_uCI, tolerance = 1e-4)
  expect_equal(negative$y,     -positive$y,     tolerance = 1e-4)
  expect_equal(negative$y_uCI, -positive$y_lCI, tolerance = 1e-4)
})

test_that("PET-PEESE prior plot data uses CDF quantiles for PET and PEESE mixtures", {

  x_seq <- c(0, 0.25, 0.5, 1)
  prior_list <- list(
    PET   = prior_PET("cauchy", list(0, 1)),
    PEESE = prior_PEESE("cauchy", list(0, 5))
  )
  prior_list_mu <- list(
    mu1 = prior("spike", list(0)),
    mu2 = prior("spike", list(0))
  )

  half_cauchy_mixture_quantile <- function(se, p){
    if(se == 0){
      return(0)
    }

    a <- se
    b <- 5 * se^2
    if(isTRUE(all.equal(p, .5))){
      return(sqrt(a * b))
    }

    tangent <- tan(pi * p)
    (-(a + b) + sign(tangent) * sqrt((a + b)^2 + 4 * tangent^2 * a * b)) /
      (2 * tangent)
  }

  plot_data <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )

  expect_null(plot_data$samples)
  expect_equal(plot_data$y_lCI, vapply(x_seq, half_cauchy_mixture_quantile, numeric(1), p = .025), tolerance = 5e-4)
  expect_equal(plot_data$y,     vapply(x_seq, half_cauchy_mixture_quantile, numeric(1), p = .500), tolerance = 5e-4)
  expect_equal(plot_data$y_uCI, vapply(x_seq, half_cauchy_mixture_quantile, numeric(1), p = .975), tolerance = 5e-4)

  dense <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = seq(0, 1, length.out = 101),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 101,
    n_samples                = 1000,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )
  expect_lt(max(abs(diff(diff(dense$y_uCI)))), 0.02)
})

test_that("PET-PEESE prior plot data transforms ordered CDF quantiles", {

  prior_list <- list(
    PET = prior_PET("normal", list(0, 1), truncation = list(-Inf, Inf))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )
  x_seq <- c(0, 0.5, 1)

  exp_plot <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = "exp",
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )
  expect_equal(exp_plot$y,     exp(rep(0, length(x_seq))), tolerance = 1e-8)
  expect_equal(exp_plot$y_lCI, exp(stats::qnorm(.025) * x_seq), tolerance = 1e-4)
  expect_equal(exp_plot$y_uCI, exp(stats::qnorm(.975) * x_seq), tolerance = 1e-4)

  lin_plot <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = length(x_seq),
    n_samples                = 1000,
    transformation           = "lin",
    transformation_arguments = list(a = 2, b = -3),
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )
  expect_equal(lin_plot$y, rep(2, length(x_seq)), tolerance = 1e-8)
  expect_equal(lin_plot$y_lCI, 2 - 3 * stats::qnorm(.975) * x_seq, tolerance = 1e-4)
  expect_equal(lin_plot$y_uCI, 2 - 3 * stats::qnorm(.025) * x_seq, tolerance = 1e-4)
})

test_that("PET-PEESE prior plot data falls back to samples for custom transformations", {

  set.seed(1)
  prior_list <- list(
    PET = prior_PET("normal", list(0, 1), truncation = list(-Inf, Inf))
  )
  prior_list_mu <- list(
    mu = prior("spike", list(0))
  )
  custom_transformation <- list(
    fun = function(x) x^2,
    inv = sqrt,
    jac = function(x) 1 / (2 * sqrt(x))
  )

  plot_data <- BayesTools:::.plot_data_prior_list.PETPEESE(
    prior_list               = prior_list,
    x_seq                    = c(0, 0.5, 1),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 3,
    n_samples                = 100,
    transformation           = custom_transformation,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    prior_list_mu            = prior_list_mu
  )

  expect_false(is.null(plot_data$samples))
  expect_equal(dim(plot_data$samples), c(100, 3))
  expect_true(all(plot_data$y_lCI <= plot_data$y))
  expect_true(all(plot_data$y <= plot_data$y_uCI))
  expect_equal(
    attr(plot_data, "y_range"),
    range(plot_data$y, plot_data$y_lCI, plot_data$y_uCI)
  )
})

test_that("PET-PEESE posterior plot data does not recycle coefficient rows", {

  make_samples <- function(x, prior) {
    attr(x, "prior_list") <- list(prior)
    attr(x, "models_ind") <- rep(1, length(x))
    class(x) <- c("mixed_posteriors", "mixed_posteriors.simple")
    x
  }

  samples <- list(
    mu    = make_samples(c(0, 1), prior("normal", list(0, 1))),
    PET   = make_samples(c(2, 3), prior_PET("normal", list(0, 1))),
    PEESE = make_samples(c(0, 0), prior("point", list(0)))
  )
  class(samples) <- c("mixed_posteriors", "list")

  plot_data <- BayesTools:::.plot_data_samples.PETPEESE(
    samples                  = samples,
    x_seq                    = c(0, 1),
    x_range                  = c(0, 1),
    x_range_quant            = NULL,
    n_points                 = 2,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE
  )

  expect_equal(dim(plot_data$samples), c(2, 2))
  expect_equal(plot_data$y, c(0.5, 3))
})

test_that("PET-PEESE posterior plot data rejects unrelated bias columns", {

  bias_samples <- cbind(
    omega = c(0.5, 0.75),
    alpha = c(1, 2)
  )

  expect_error(
    BayesTools:::.plot_data_samples.PETPEESE(
      samples                  = list(mu = c(0, 1), bias = bias_samples),
      x_seq                    = c(0, 1),
      x_range                  = c(0, 1),
      x_range_quant            = NULL,
      n_points                 = 2,
      transformation           = NULL,
      transformation_arguments = NULL,
      transformation_settings  = FALSE
    ),
    "At least one 'PET' or 'PEESE' model needs to be specified.",
    fixed = TRUE
  )
})

test_that("PET-PEESE posterior plot data honors negative effect direction", {
  samples <- list(
    mu    = c(0, 1, 2, 3),
    PET   = c(1, 2, 3, 4),
    PEESE = c(0.5, 0, -0.5, 1)
  )
  x_seq <- c(0, .5, 1)

  plot_data <- BayesTools:::.plot_data_samples.PETPEESE(
    samples = samples,
    x_seq = x_seq,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = length(x_seq),
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE,
    effect_direction = "negative"
  )

  expected_samples <- sapply(x_seq, function(x) samples$mu - samples$PET * x - samples$PEESE * x^2)
  expected_quantiles <- apply(expected_samples, 2, stats::quantile, probs = c(.500, .025, .975), names = FALSE)

  expect_equal(plot_data$samples, expected_samples)
  expect_equal(plot_data$y, expected_quantiles[1,])
  expect_equal(plot_data$y_lCI, expected_quantiles[2,])
  expect_equal(plot_data$y_uCI, expected_quantiles[3,])
  expect_equal(
    attr(plot_data, "y_range"),
    range(plot_data$y, plot_data$y_lCI, plot_data$y_uCI)
  )

  intercept_fallback <- BayesTools:::.plot_data_samples.PETPEESE(
    samples = list(mu_intercept = samples$mu, PET = samples$PET, PEESE = samples$PEESE),
    x_seq = x_seq,
    x_range = c(0, 1),
    x_range_quant = NULL,
    n_points = length(x_seq),
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE,
    effect_direction = "negative"
  )
  expect_equal(intercept_fallback$samples, plot_data$samples)
  expect_equal(intercept_fallback$y, plot_data$y)
})

test_that("factor posterior density curves keep continuous mass scale", {

  set.seed(1)
  n_samples <- 5000
  samples <- cbind(
    `mu_alloc[random]`     = stats::rnorm(n_samples, -0.5, 0.8),
    `mu_alloc[systematic]` = stats::rnorm(n_samples,  0.5, 0.8)
  )

  prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior, "levels") <- 3
  attr(prior, "level_names") <- c("alternate", "random", "systematic")

  attr(samples, "prior_list") <- prior
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples = list(mu_alloc = samples),
    parameter = "mu_alloc",
    n_points = 512,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.simple")]
  areas <- vapply(density_entries, function(density){
    sum(density$y) * (density$x[2] - density$x[1])
  }, numeric(1))

  expect_equal(unname(areas), c(1, 1), tolerance = 0.08)
})

test_that("bounded factor posterior KDE reflects each level support", {
  n_samples <- 4000
  samples <- cbind(
    `theta[A]` = seq(.005, .995, length.out = n_samples),
    `theta[B]` = rev(seq(.005, .995, length.out = n_samples))
  )

  prior <- prior_factor("uniform", list(0, 1), contrast = "independent")
  attr(prior, "levels") <- 2
  attr(prior, "level_names") <- c("A", "B")

  attr(samples, "prior_list") <- prior
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector", "matrix")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples = list(theta = samples),
    parameter = "theta",
    n_points = 512,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.simple")]
  areas <- vapply(density_entries, .model_plot_density_mass, numeric(1))
  ranges <- t(vapply(density_entries, function(density) range(density$x), numeric(2)))

  expect_length(density_entries, 2)
  expect_equal(unname(vapply(density_entries, attr, numeric(1), which = "level")), 1:2)
  expect_equal(unname(vapply(density_entries, attr, character(1), which = "level_name")), colnames(samples))
  expect_equal(unname(areas), c(1, 1), tolerance = .03)
  expect_equal(unname(ranges), matrix(c(0, 1, 0, 1), ncol = 2, byrow = TRUE))
  expect_true(all(vapply(density_entries, function(density){
    isTRUE(attr(density, "boundary_reflection"))
  }, logical(1))))
})

test_that("transformed bounded factor posterior KDE avoids singular endpoints", {
  n_samples <- 4000
  samples <- cbind(
    `theta[A]` = seq(.005, .995, length.out = n_samples),
    `theta[B]` = rev(seq(.005, .995, length.out = n_samples))
  )

  prior <- prior_factor("uniform", list(0, 1), contrast = "independent")
  attr(prior, "levels") <- 2
  attr(prior, "level_names") <- c("A", "B")

  attr(samples, "prior_list") <- prior
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector", "matrix")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples = list(theta = samples),
    parameter = "theta",
    n_points = 512,
    transformation = "exp_lin",
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.simple")]

  expect_length(density_entries, 2)
  expect_true(all(vapply(density_entries, function(density){
    all(is.finite(density$x)) && all(is.finite(density$y))
  }, logical(1))))
  expect_true(all(vapply(density_entries, function(density){
    min(density$x) > 0 && max(density$x) < 1
  }, logical(1))))
  expect_true(all(vapply(density_entries, function(density){
    isTRUE(attr(density, "boundary_reflection"))
  }, logical(1))))
})

test_that("factor posterior plot data uses level-matched stored densities", {

  n_samples <- 100
  samples <- cbind(
    `mu_alloc[random]`     = seq(-1, 0, length.out = n_samples),
    `mu_alloc[systematic]` = seq(1, 2, length.out = n_samples)
  )

  prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior, "levels") <- 3
  attr(prior, "level_names") <- c("alternate", "random", "systematic")

  stored_random_x <- seq(-2, 1, length.out = 31)
  stored_systematic_x <- seq(0, 3, length.out = 31)
  attr(samples, "prior_list") <- prior
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "posterior_density") <- list(
    random = list(
      parameter = "random",
      x         = stored_random_x,
      y         = stats::dnorm(stored_random_x, mean = -.5, sd = .5),
      method    = "iwmde"
    ),
    systematic = list(
      parameter = "systematic",
      x         = stored_systematic_x,
      y         = stats::dnorm(stored_systematic_x, mean = 1.5, sd = .5),
      method    = "iwmde"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples                  = list(mu_alloc = samples),
    parameter                = "mu_alloc",
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  expect_equal(plot_data$density1$x, stored_random_x)
  expect_equal(plot_data$density2$x, stored_systematic_x)
  expect_equal(attr(plot_data$density1, "posterior_density_method"), "iwmde")
  expect_equal(attr(plot_data$density2, "posterior_density_method"), "iwmde")
})

test_that("factor posterior plot data uses stored point masses once", {

  n_samples <- 100
  samples <- cbind(
    `mu_alloc[random]`     = c(rep(0, 40), seq(-1, 0, length.out = 60)),
    `mu_alloc[systematic]` = c(rep(0, 40), seq(1, 2, length.out = 60))
  )

  prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior, "levels") <- 3
  attr(prior, "level_names") <- c("alternate", "random", "systematic")

  stored_x <- seq(-2, 2, length.out = 31)
  attr(samples, "prior_list") <- list(
    prior_factor("point", list(location = 0), contrast = "treatment"),
    prior
  )
  attr(samples, "models_ind") <- c(rep(1, 40), rep(2, 60))
  attr(samples, "posterior_density") <- list(
    random = list(
      parameter    = "random",
      x            = stored_x,
      y            = stats::dnorm(stored_x),
      point_masses = data.frame(x = 0, mass = .25),
      method       = "iwmde"
    ),
    systematic = list(
      parameter    = "systematic",
      x            = stored_x,
      y            = stats::dnorm(stored_x),
      point_masses = data.frame(x = 0, mass = .25),
      method       = "iwmde"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples                  = list(mu_alloc = samples),
    parameter                = "mu_alloc",
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  point_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.point")]

  expect_false(any(grepl("^points", names(plot_data))))
  expect_equal(length(point_entries), 2L)
  expect_equal(
    unname(vapply(point_entries, function(point) point[["y"]], numeric(1))),
    c(.25, .25)
  )
})

test_that("factor posterior plot data keeps fallback spikes per level", {

  n_samples <- 100
  samples <- cbind(
    `mu_alloc[random]`     = c(rep(0, 40), seq(-1, 0, length.out = 60)),
    `mu_alloc[systematic]` = c(rep(0, 40), seq(1, 2, length.out = 60))
  )

  prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  point_prior <- prior_factor("point", list(location = 0), contrast = "treatment")
  attr(prior, "levels") <- 3
  attr(prior, "level_names") <- c("alternate", "random", "systematic")
  attr(point_prior, "levels") <- 3
  attr(point_prior, "level_names") <- c("alternate", "random", "systematic")

  stored_x <- seq(-2, 2, length.out = 31)
  attr(samples, "prior_list") <- list(point_prior, prior)
  attr(samples, "models_ind") <- c(rep(1, 40), rep(2, 60))
  attr(samples, "posterior_density") <- list(
    random = list(
      parameter    = "random",
      x            = stored_x,
      y            = stats::dnorm(stored_x),
      point_masses = data.frame(x = 0, mass = .25),
      method       = "iwmde"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples                  = list(mu_alloc = samples),
    parameter                = "mu_alloc",
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  point_entries <- plot_data[vapply(plot_data, inherits, logical(1), "density.prior.point")]
  point_levels <- vapply(point_entries, function(point) attr(point, "level"), numeric(1))

  expect_setequal(unname(point_levels), c(1, 2))
  expect_equal(point_entries[[which(point_levels == 1)]]$y, .25)
  expect_equal(point_entries[[which(point_levels == 2)]]$y, .40)
})

test_that("factor posterior plot data matches interaction cell aliases", {

  n_samples <- 80
  samples <- cbind(
    `mu_a[A]__xXx__b[B]` = seq(-1, 0, length.out = n_samples),
    `mu_a[C]__xXx__b[D]` = seq(1, 2, length.out = n_samples)
  )

  prior <- prior_factor("normal", list(0, 1), contrast = "independent")
  attr(prior, "levels") <- 2
  attr(prior, "level_names") <- list(a = c("A", "C"), b = c("B", "D"))
  attr(prior, "factor_cell_names") <- c("A, B", "C, D")

  stored_ab_x <- seq(-2, 1, length.out = 31)
  stored_cd_x <- seq(0, 3, length.out = 31)
  attr(samples, "prior_list") <- prior
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "level_names") <- attr(prior, "level_names")
  attr(samples, "factor_cell_names") <- attr(prior, "factor_cell_names")
  attr(samples, "posterior_density") <- list(
    "A, B" = list(
      parameter = "A, B",
      x         = stored_ab_x,
      y         = stats::dnorm(stored_ab_x, mean = -.5, sd = .5),
      method    = "iwmde"
    ),
    "C, D" = list(
      parameter = "C, D",
      x         = stored_cd_x,
      y         = stats::dnorm(stored_cd_x, mean = 1.5, sd = .5),
      method    = "iwmde"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples                  = list(mu_a__xXx__b = samples),
    parameter                = "mu_a__xXx__b",
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  expect_equal(plot_data$density1$x, stored_ab_x)
  expect_equal(plot_data$density2$x, stored_cd_x)
  expect_equal(attr(plot_data$density1, "posterior_density_method"), "iwmde")
  expect_equal(attr(plot_data$density2, "posterior_density_method"), "iwmde")
})

test_that("density_method is named-only on exported plot APIs", {

  expect_gt(
    match("density_method", names(formals(plot_posterior))),
    match("...", names(formals(plot_posterior)))
  )
  expect_gt(
    match("density_method", names(formals(plot_marginal))),
    match("...", names(formals(plot_marginal)))
  )
})

test_that("factor posterior plot data aggregates duplicate point-mass models", {

  samples_matrix <- matrix(
    c(
      rep(0, 10),
      rep(2, 10)
    ),
    ncol = 2,
    dimnames = list(NULL, c("theta[A]", "theta[B]"))
  )
  attr(samples_matrix, "models_ind") <- c(rep(1, 2), rep(2, 3), rep(3, 5))
  attr(samples_matrix, "prior_list") <- list(
    prior_factor("point", list(location = 0), contrast = "treatment"),
    prior_factor("point", list(location = 0), contrast = "treatment"),
    prior_factor("point", list(location = 2), contrast = "treatment")
  )
  class(samples_matrix) <- c(
    "mixed_posteriors",
    "mixed_posteriors.factor",
    "mixed_posteriors.vector",
    "matrix"
  )

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples = list(theta = samples_matrix),
    parameter = "theta",
    n_points = 64,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(names(plot_data), c("points1", "points2"))
  expect_equal(plot_data$points1$x, 0)
  expect_equal(plot_data$points1$y, .5)
  expect_equal(plot_data$points2$x, 2)
  expect_equal(plot_data$points2$y, .5)
  expect_true(all(vapply(plot_data, inherits, logical(1), what = "density.prior.point")))
})

.scaled_multi_factor_interaction_samples <- function(n_posterior = 20, n_prior = 128){

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3")),
    x = seq_len(4)
  )
  formula_result <- JAGS_formula(
    formula = ~ x * a * b,
    parameter = "mu",
    data = df,
    formula_scale = list(x = TRUE),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      a = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "x:a" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x:b" = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x:a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  prior_list <- formula_result$prior_list
  posterior_columns <- unlist(lapply(names(prior_list), function(parameter) {
    if(is.prior.factor(prior_list[[parameter]])){
      BayesTools:::.JAGS_prior_factor_names(parameter, prior_list[[parameter]])
    }else{
      parameter
    }
  }), use.names = FALSE)

  set.seed(515)
  posterior <- matrix(rnorm(n_posterior * length(posterior_columns)), nrow = n_posterior)
  colnames(posterior) <- posterior_columns
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  samples <- as_mixed_posteriors(
    fit,
    parameters = "mu_x__xXx__a__xXx__b",
    transform_scaled = TRUE,
    n_prior_samples = n_prior
  )

  list(samples = samples, prior_list = prior_list)
}

.scaled_default_treatment_interaction_samples <- function(n_posterior = 20, n_prior = 128){

  df <- data.frame(
    alloc = factor(
      rep(c("alternate", "random", "systematic"), each = 4),
      levels = c("alternate", "random", "systematic")
    ),
    year = seq(1960, 1971, length.out = 12)
  )
  formula_result <- JAGS_formula(
    formula       = ~ alloc * year,
    parameter     = "mu",
    data          = df,
    formula_scale = list(year = TRUE),
    prior_list    = list(
      "__default_continuous" = prior("normal", list(0, 1)),
      "__default_factor"     = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  prior_list <- formula_result$prior_list
  posterior_columns <- unlist(lapply(names(prior_list), function(parameter) {
    if(is.prior.factor(prior_list[[parameter]])){
      BayesTools:::.JAGS_prior_factor_names(parameter, prior_list[[parameter]])
    }else{
      parameter
    }
  }), use.names = FALSE)

  set.seed(616)
  posterior <- matrix(rnorm(n_posterior * length(posterior_columns)), nrow = n_posterior)
  colnames(posterior) <- posterior_columns
  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  samples <- as_mixed_posteriors(
    fit,
    parameters = "mu_alloc",
    transform_scaled = TRUE,
    n_prior_samples = n_prior
  )

  list(samples = samples, prior_list = prior_list)
}

test_that("transformed factor prior plot data uses indexed prior density columns", {

  fixture <- .scaled_multi_factor_interaction_samples()
  samples <- fixture$samples
  prior_list <- fixture$prior_list

  plot_data <- BayesTools:::.plot_data_prior_factor_density_transformed(
    prior_density_context = attr(samples, "prior_density_context"),
    samples = samples,
    parameter = "mu_x__xXx__a__xXx__b",
    prior_list = attr(samples$mu_x__xXx__a__xXx__b, "prior_list"),
    n_points = 32
  )

  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.factor")]
  expect_equal(length(density_entries), 6L)
  expect_equal(
    unname(vapply(density_entries, attr, character(1), which = "level_name")),
    paste0(
      "mu_x__xXx__a[dif: ",
      rep(c("a1", "a2"), times = 3),
      "]__xXx__b[dif: ",
      rep(c("b1", "b2", "b3"), each = 2),
      "]"
    )
  )
})

test_that("transformed treatment factor prior plot data omits reference coefficient", {

  fixture <- .scaled_default_treatment_interaction_samples()
  samples <- fixture$samples

  factor_weights <- BayesTools:::.prior_factor_level_weight_matrix(
    sample_metadata = samples$mu_alloc,
    parameter       = "mu_alloc",
    samples         = samples
  )
  expect_equal(rownames(factor_weights), colnames(samples$mu_alloc))
  expect_equal(colnames(factor_weights), c("mu_alloc[1]", "mu_alloc[2]"))
  expect_false(any(rowSums(abs(factor_weights)) == 0))

  plot_data <- BayesTools:::.plot_data_prior_factor_density_transformed(
    prior_density_context = attr(samples, "prior_density_context"),
    samples               = samples,
    parameter             = "mu_alloc",
    prior_list            = attr(samples$mu_alloc, "prior_list"),
    n_points              = 32
  )

  expect_equal(length(plot_data), 2L)
  expect_false(any(vapply(plot_data, inherits, logical(1), what = "density.prior.point")))
  expect_equal(
    unname(vapply(plot_data, attr, character(1), which = "level_name")),
    colnames(samples$mu_alloc)
  )
})

test_that("plot_posterior handles transformed scaled factor-interaction priors", {

  skip_if_not_installed("ggplot2")

  fixture <- .scaled_multi_factor_interaction_samples(n_posterior = 24, n_prior = 128)

  expect_silent({
    plot <- plot_posterior(
      samples = fixture$samples,
      parameter = "mu_x__xXx__a__xXx__b",
      plot_type = "ggplot",
      prior = TRUE,
      n_points = 32
    )
  })

  expect_s3_class(plot, "ggplot")
  expect_gt(length(plot$layers), 1)
})

test_that("plot_posterior owns factor legends in transformed prior overlays", {

  fixture <- .scaled_default_treatment_interaction_samples(n_posterior = 20, n_prior = 128)

  expect_length(
    .capture_graphics_legend_calls_for_test(
      plot_posterior(fixture$samples, "mu_alloc", prior = TRUE, n_points = 32)
    ),
    1L
  )
  expect_length(
    .capture_graphics_legend_calls_for_test(
      plot_posterior(
        fixture$samples,
        "mu_alloc",
        prior = TRUE,
        n_points = 32,
        legend = FALSE,
        dots_prior = list(legend = TRUE)
      )
    ),
    0L
  )
})
