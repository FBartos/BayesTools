skip_if_not_test_profile(c("unit", "visual-fixture"))

# ============================================================================ #
# TEST FILE: JAGS Marginal Distributions
# ============================================================================ #
#
# PURPOSE:
#   Tests for marginal_posterior, ensemble_inference, mix_posteriors,
#   and related functions. Uses pre-fitted models from test-00-model-fits.R.
#
# DEPENDENCIES:
#   - rjags, bridgesampling: JAGS model fitting and marginal likelihood
#   - common-functions.R: temp_fits_dir, skip_if_no_fits, test_reference_table
#
# SKIP CONDITIONS:
#   - skip_if_no_fits(): Pre-fitted models required
#   - skip_if_not_installed("rjags"), skip_if_not_installed("bridgesampling")
#   - skip_on_os(): Multivariate sampling differs across OSes (meandif priors)
#
# MODELS/FIXTURES:
#   - fit_marginal_0, fit_marginal_1
#
# TAGS: @evaluation, @JAGS, @model-averaging, @marginal
# ============================================================================ #

# Reference directory for table outputs
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-marginal-distributions")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

.expect_marginal_table_current_inputs <- function(table, samples, inference,
                                                  parameters,
                                                  probs = c(0.025, 0.5, 0.975)){

  draw_groups <- list()
  inference_groups <- list()
  for(parameter in parameters){
    parameter_samples <- samples[[parameter]]
    parameter_inference <- inference[[parameter]]
    if(is.list(parameter_samples)){
      parameter_draws <- lapply(parameter_samples, as.numeric)
    }else{
      parameter_draws <- list(as.numeric(parameter_samples))
    }
    if(is.list(parameter_samples) && length(parameter_samples) > 1L){
      parameter_inferences <- lapply(
        seq_along(parameter_draws),
        function(i) parameter_inference[[i]]
      )
    }else{
      parameter_inferences <- list(parameter_inference[[1L]])
    }
    draw_groups <- c(draw_groups, parameter_draws)
    inference_groups <- c(inference_groups, parameter_inferences)
  }

  estimate_names <- c("Mean", "SD", as.character(probs))
  expected_estimates <- t(vapply(draw_groups, function(draws){
    c(
      Mean = mean(draws),
      SD = stats::sd(draws),
      vapply(
        probs,
        function(prob) unname(stats::quantile(draws, probs = prob)),
        numeric(1)
      )
    )
  }, numeric(length(estimate_names))))
  colnames(expected_estimates) <- estimate_names

  expect_equal(
    unname(as.matrix(table[, estimate_names, drop = FALSE])),
    unname(expected_estimates),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(table$inclusion_BF),
    vapply(inference_groups, as.numeric, numeric(1)),
    tolerance = 1e-12
  )

  expected_BF_errors <- vapply(inference_groups, function(x){
    error <- attr(x, "BF_error_percent")
    if(is.null(error)){
      return(NA_real_)
    }
    error <- as.numeric(error)[1L]
    if(!is.finite(error) || error < 0){
      return(NA_real_)
    }
    error
  }, numeric(1))
  if("BF_error_percent" %in% names(table)){
    expect_equal(
      as.numeric(table$BF_error_percent),
      expected_BF_errors,
      tolerance = 1e-12
    )
  }else{
    expect_false(any(is.finite(expected_BF_errors)))
  }
}

.plot_prior_density_for_test <- function(x, main = "", xlim = NULL, ylim = NULL, add = FALSE,
                                         lty = 1, col = graphics::par("fg"), ...){
  prior_density <- .bt_meta_get(x, "prior_density")
  if(is.null(prior_density)){
    stop("The object does not contain a deterministic prior density.", call. = FALSE)
  }
  plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    prior_density,
    n_points = 512,
    x_range  = xlim
  )

  if(length(plot_data) == 0){
    if(is.null(xlim)) xlim <- c(-1, 1)
    if(is.null(ylim)) ylim <- c(0, 1)
  }else{
    if(is.null(xlim)){
      xlim <- range(unlist(lapply(plot_data, function(d) attr(d, "x_range"))), finite = TRUE)
      if(!all(is.finite(xlim)) || diff(xlim) <= 0) xlim <- xlim + c(-1, 1)
    }
    if(is.null(ylim)){
      y_max <- max(unlist(lapply(plot_data, function(d) attr(d, "y_range")[2])), na.rm = TRUE)
      ylim <- c(0, if(is.finite(y_max) && y_max > 0) y_max else 1)
    }
  }

  if(!add){
    graphics::plot(NA, xlim = xlim, ylim = ylim, main = main, xlab = "", ylab = "Density", ...)
  }
  for(d in plot_data){
    if(inherits(d, "density.prior.point")){
      BayesTools:::.lines.prior.point(d, scale_y2 = 1, lty = lty, col = col)
    }else{
      graphics::lines(d$x, d$y, lty = lty, col = col)
    }
  }
  invisible(plot_data)
}

.posterior_density_for_test <- function(x, y, method = "iwmde",
                                        density_method = "precomputed", ...){
  posterior_density_attribute(x = x, y = y, method = method,
                              density_method = density_method, ...)
}

.posterior_ordinate_for_test <- function(value, ordinate, method = "qCMDE",
                                         density_method = "precomputed", ...){
  posterior_ordinate_attribute(value = value, ordinate = ordinate,
                               method = method,
                               density_method = density_method, ...)
}

.marginal_posterior_with_prior_density_for_test <- function(samples, prior_density) {
  class(samples) <- c("marginal_posterior.simple", "marginal_posterior", class(samples))
  samples <- .bt_meta_set(samples, "prior_density", prior_density)
  samples <- .bt_meta_set(samples, "atoms", posterior_atom_attribute())
  samples
}

.marginal_semantic_fixture_for_test <- function(){

  df <- expand.grid(
    x_cont1  = c(-1, 0, 1),
    x_fac2t  = factor(c("A", "B"), levels = c("A", "B")),
    x_fac3md = factor(c("A", "B", "C"), levels = c("A", "B", "C"))
  )

  formula_result <- JAGS_formula(
    formula = ~ x_cont1 + x_fac2t + x_cont1 * x_fac3md,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept        = prior("normal", list(0, 1)),
      x_cont1          = prior("normal", list(0, 1)),
      x_fac2t          = prior_factor("normal", list(0, 1), contrast = "treatment"),
      x_fac3md         = prior_factor("mnormal", list(0, 0.25), contrast = "meandif"),
      "x_cont1:x_fac3md" = prior_factor("mnormal", list(0, 0.25), contrast = "meandif")
    )
  )

  posterior <- cbind(
    mu_intercept = seq(-0.4, 0.5, length.out = 10),
    mu_x_cont1 = seq(-0.2, 0.7, length.out = 10),
    mu_x_fac2t = seq(0.1, 1.0, length.out = 10),
    `mu_x_fac3md[1]` = seq(-0.6, 0.3, length.out = 10),
    `mu_x_fac3md[2]` = seq(0.2, 1.1, length.out = 10),
    `mu_x_cont1__xXx__x_fac3md[1]` = seq(-1.0, -0.1, length.out = 10),
    `mu_x_cont1__xXx__x_fac3md[2]` = seq(0.3, 1.2, length.out = 10)
  )

  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(
    fit,
    parameters = c(
      "mu_intercept",
      "mu_x_cont1",
      "mu_x_fac2t",
      "mu_x_fac3md",
      "mu_x_cont1__xXx__x_fac3md"
    ),
    n_prior_samples = 128
  )

  list(
    data = df,
    posterior = posterior,
    prior_list = formula_result$prior_list,
    samples = samples
  )
}

.expect_density_area_for_test <- function(plot_data, expected = 1, tolerance = 0.08){
  dx <- diff(plot_data$x)
  area <- sum(dx * (plot_data$y[-1] + plot_data$y[-length(plot_data$y)]) / 2)
  expect_equal(area, expected, tolerance = tolerance)
}

.density_component_mass_for_test <- function(plot_data) {
  if(inherits(plot_data, "density.prior.point")){
    return(sum(plot_data$y))
  }

  dx <- diff(plot_data$x)
  sum(dx * (plot_data$y[-1] + plot_data$y[-length(plot_data$y)]) / 2)
}

.density_mass_by_level_for_test <- function(plot_data) {
  level_names <- vapply(plot_data, function(component) {
    level_name <- attr(component, "level_name")
    if(is.null(level_name)) "__all__" else level_name
  }, character(1))

  vapply(unique(level_names), function(level_name) {
    sum(vapply(
      plot_data[level_names == level_name],
      .density_component_mass_for_test,
      numeric(1)
    ))
  }, numeric(1))
}

.expect_density_mass_by_level_for_test <- function(plot_data, expected = 1, tolerance = 0.08) {
  masses <- .density_mass_by_level_for_test(plot_data)
  expect_true(all(is.finite(masses)))
  expect_true(all(masses >= 0))
  expect_equal(unname(masses), rep(expected, length(masses)), tolerance = tolerance)
  invisible(masses)
}

test_that("posterior density method helpers validate density sources", {

  expect_equal(posterior_density_method_match("KDE"), "KDE")
  expect_equal(posterior_density_method_match("precomputed"), "precomputed")
  expect_error(
    posterior_density_method_match("bad", name = "density_method"),
    "density_method"
  )

  expect_false(posterior_density_method_uses_precomputed("KDE"))
  expect_true(posterior_density_method_uses_precomputed("precomputed"))
  # estimator names of other packages are not BayesTools density methods:
  # those packages map them to "precomputed" themselves
  expect_false(posterior_density_method_uses_precomputed("qCMDE"))
  expect_false(posterior_density_method_uses_precomputed("IWMDE"))
  expect_error(posterior_density_method_uses_precomputed(c("precomputed", "KDE")),
    "The 'method' argument must have length '1'.", fixed = TRUE)
})

test_that("Savage_Dickey_BF uses prior density over normal posterior height", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior_z <- seq(-4, 4, length.out = 1001)
  posterior <- 0.4 + (posterior_z - mean(posterior_z)) / stats::sd(posterior_z) * 1.3
  posterior <- .marginal_posterior_with_prior_density_for_test(posterior, prior_density)

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) /
    stats::dnorm(0, mean = mean(posterior), sd = stats::sd(posterior))

  out <- Savage_Dickey_BF(posterior, null_hypothesis = 0, normal_approximation = TRUE, silent = TRUE)
  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "posterior_density_source"), "normal")
})

test_that("Savage_Dickey_BF uses stored posterior density when available", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  stored_x <- seq(-4, 4, length.out = 401)
  stored_y <- stats::dnorm(stored_x, mean = 0.4, sd = 1.2)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  ))

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) /
    stats::approx(stored_x, stored_y, xout = 0)[["y"]]

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )
  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "posterior_density_source"), "precomputed")
})

test_that("Savage_Dickey_BF rejects mismatched direct precomputed attributes", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  attr(posterior, "parameter") <- "theta"

  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    parameter = "phi",
    value     = 0,
    ordinate  = .5,
    method    = "qCMDE"
  ))
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      silent               = TRUE,
      density_method       = "precomputed"
    ),
    "requires valid posterior ordinate or posterior density metadata",
    fixed = TRUE
  )

  posterior <- .bt_meta_set(posterior, "posterior_ordinate", NULL)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    parameter = "phi",
    x         = seq(-1, 1, length.out = 101),
    y         = rep(.5, 101),
    method    = "qCMDE"
  ))
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      silent               = TRUE,
      density_method       = "precomputed"
    ),
    "requires valid posterior ordinate or posterior density metadata",
    fixed = TRUE
  )
})

test_that("posterior density and ordinate constructors create reusable attributes", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )

  ordinate <- posterior_ordinate_attribute(
    value          = 0,
    ordinate       = .5,
    method         = "iwmde",
    density_method = "IWMDE",
    diagnostics    = list(BF_error_percent = 2.5),
    parameter      = "theta"
  )
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", ordinate)

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis = 0,
    silent          = TRUE,
    density_method  = "precomputed"
  )
  expect_equal(attr(out, "posterior_density_source"), "precomputed")
  expect_equal(attr(out, "BF_error_percent"), 2.5)
  expect_true(posterior_ordinate_has_value(ordinate, 0))
  expect_true(posterior_ordinate_supports_bf(ordinate))
  expect_false(posterior_ordinate_supports_bf(
    ordinate,
    validator = function(x) FALSE
  ))

  second <- posterior_ordinate_attribute(
    value          = 1,
    ordinate       = .25,
    method         = "iwmde",
    density_method = "IWMDE"
  )
  appended <- posterior_ordinate_append(ordinate, second)
  expect_true(posterior_ordinate_has_value(appended, 0))
  expect_true(posterior_ordinate_has_value(appended, 1))
  expect_true(posterior_ordinate_supports_bf(
    appended,
    validator = function(x) !is.null(x[["diagnostics"]])
  ))
  expect_error(
    posterior_ordinate_append(NULL, list(value = 0)),
    "valid posterior ordinate"
  )
  expect_error(
    posterior_ordinate_attribute(
      value          = c(0, 0),
      ordinate       = c(.5, .6),
      method         = "iwmde",
      density_method = "IWMDE"
    ),
    "unique"
  )
  expect_error(
    posterior_ordinate_append(ordinate, ordinate),
    "duplicate"
  )
  expect_error(
    posterior_ordinate_attribute(
      value          = 0,
      ordinate       = .5,
      method         = "iwmde",
      density_method = "IWMDE",
      x              = 0
    ),
    "reserved fields"
  )
  expect_error(
    posterior_ordinate_attribute(
      value          = 0,
      ordinate       = .5,
      method         = "iwmde",
      density_method = "IWMDE",
      parameter      = "theta",
      parameter      = "mu"
    ),
    "unique, nonmissing names",
    fixed = TRUE
  )

  density <- posterior_density_attribute(
    x              = seq(-2, 2, length.out = 201),
    y              = .99 * stats::dnorm(seq(-2, 2, length.out = 201)),
    method         = "iwmde",
    density_method = "IWMDE",
    support        = posterior_support_attribute(c(-Inf, Inf), points = 0),
    parameter      = "theta"
  )
  parsed_density <- BayesTools:::.posterior_density_from_attribute(density)
  expect_equal(parsed_density[["method"]], "iwmde")
  expect_equal(parsed_density[["support"]][["bounds"]], c(-Inf, Inf))
  expect_equal(parsed_density[["support"]][["points"]], 0)
  # the point masses of the posterior are declared as atoms of the draws
  expect_false(any(c("point_masses", "point_masses_declared") %in% names(parsed_density)))
  expect_error(
    posterior_density_attribute(
      x              = seq(-2, 2, length.out = 201),
      y              = .99 * stats::dnorm(seq(-2, 2, length.out = 201)),
      method         = "iwmde",
      density_method = "IWMDE",
      point_masses   = data.frame(x = 0, mass = .01)
    ),
    "Posterior densities do not carry 'point_masses'",
    fixed = TRUE
  )
  expect_error(
    posterior_density_attribute(
      x              = c(0, 1, Inf),
      y              = c(1, 1, 1),
      method         = "iwmde",
      density_method = "IWMDE"
    ),
    "grid values must be finite",
    fixed = TRUE
  )
  expect_error(
    posterior_density_attribute(
      x              = 0:1,
      y              = c(1, 1),
      method         = "iwmde",
      density_method = "IWMDE",
      parameter      = "theta",
      parameter      = "mu"
    ),
    "unique, nonmissing names",
    fixed = TRUE
  )
  expect_error(
    posterior_density_attribute(
      x              = 0:1,
      y              = c(1, 1),
      method         = "iwmde",
      density_method = "IWMDE",
      density        = .posterior_density_for_test(x = 0:1, y = c(1, 1))
    ),
    "reserved fields"
  )
  for(raw_support in list(c(0, 1), list(bounds = c(0, 1)), list(lower = 0, upper = 1))){
    expect_error(
      posterior_density_attribute(
        x              = 0:1,
        y              = c(1, 1),
        method         = "iwmde",
        density_method = "IWMDE",
        support        = raw_support
      ),
      "'support' must be created with 'posterior_support_attribute()'.",
      fixed = TRUE
    )
  }
})

test_that("posterior_support_attribute validates and classifies support", {

  interval <- posterior_support_attribute(c(0, Inf), source = "test")
  expect_s3_class(interval, "BayesTools_posterior_support")
  expect_equal(interval$bounds, c(0, Inf))
  expect_equal(interval$type, "interval")
  expect_true(interval$exact)
  expect_equal(interval$source, "test")

  points <- posterior_support_attribute(c(0, 1), points = c(0, 1))
  expect_equal(points$type, "points")
  expect_false(BayesTools:::.posterior_support_contains_value(points, .5))

  mixed <- posterior_support_attribute(c(0, 1), points = 0, type = "mixed",
                                       exact = FALSE)
  expect_equal(mixed$type, "mixed")
  expect_false(mixed$exact)

  expect_error(posterior_support_attribute(c(1, 0)),
               "'bounds' must be ordered", fixed = TRUE)
  expect_error(posterior_support_attribute(c(0, NA)),
               "cannot contain NA/NaN", fixed = TRUE)
  expect_error(posterior_support_attribute(0),
               "must have length '2'", fixed = TRUE)
  expect_error(posterior_support_attribute(c(0, 1), points = Inf),
               "'points' must be finite.", fixed = TRUE)
  expect_error(posterior_support_attribute(c(0, 1), type = "points"),
               "requires 'points'", fixed = TRUE)
  expect_error(posterior_support_attribute(c(0, 1), type = "ray"),
               "not recognized", fixed = TRUE)

  # an unclassed list of supports keyed by column carries no single support
  expect_null(BayesTools:::.posterior_support_from_attribute(list(
    theta = interval
  )))
  modified <- interval
  modified$bounds <- c(2, 1)
  expect_error(
    BayesTools:::.posterior_support_from_attribute(modified),
    "Posterior support metadata is invalid",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF validates scalar options and rejects invalid precomputed source", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )

  expect_error(
    Savage_Dickey_BF(posterior, null_hypothesis = Inf, silent = TRUE),
    "must be finite",
    fixed = TRUE
  )
  expect_error(
    Savage_Dickey_BF(posterior, null_hypothesis = NA_real_, silent = TRUE),
    "cannot contain NA/NaN",
    fixed = TRUE
  )
  expect_error(
    Savage_Dickey_BF(posterior, normal_approximation = NA, silent = TRUE),
    "cannot contain NA/NaN",
    fixed = TRUE
  )
  expect_error(
    Savage_Dickey_BF(posterior, silent = NA),
    "cannot contain NA/NaN",
    fixed = TRUE
  )

  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(2, 3, length.out = 101),
    y      = rep(1, 101),
    method = "qCMDE"
  ))
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      silent               = TRUE,
      density_method       = "precomputed"
    ),
    "Stored posterior density does not span",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF reports stored density BF error only for matched nulls", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  stored_x <- seq(-4, 4, length.out = 401)
  stored_y <- stats::dnorm(stored_x, mean = 0.4, sd = 1.2)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x           = stored_x,
    y           = stored_y,
    method      = "iwmde",
    diagnostics = list(bf_relative_mcse = .123)
  ))

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_null(attr(out, "BF_error_percent"))

  stored_density <- .bt_meta_get(posterior, "posterior_density")
  stored_density$diagnostics <- list(
    bf_value         = 0,
    bf_relative_mcse = .123
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", stored_density)
  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(attr(out, "BF_error_percent"), 12.3)
})

test_that("Savage_Dickey_BF prefers matching stored ordinates", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(2, 3, length.out = 101),
    y      = rep(100, 101),
    method = "plot-only"
  ))
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value       = c(0, .5),
    ordinate    = c(.25, .5),
    method      = "qCMDE",
    diagnostics = list(relative_mcse = c(.1, .2))
  ))

  expected <- BayesTools:::.prior_linear_density_height(prior_density, .5) / .5
  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = .5,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "BF_error_percent"), 20)
})

test_that("stored posterior ordinate parser reads constructor attributes only", {

  ordinate <- .posterior_ordinate_for_test(
    value       = c(0, .5),
    ordinate    = c(.25, .5),
    diagnostics = list(relative_mcse = c(.1, .2))
  )
  parsed <- BayesTools:::.posterior_ordinate_from_attribute(
    ordinate,
    null_hypothesis = .5
  )
  expect_equal(parsed$x, .5)
  expect_equal(parsed$y, .5)
  expect_equal(parsed$diagnostics$relative_mcse, .2)
  expect_null(BayesTools:::.posterior_ordinate_from_attribute(
    ordinate,
    null_hypothesis = 1
  ))

  # raw lists and data frames are not posterior ordinate metadata
  raw_message <- paste0(
    "Posterior ordinate metadata must be created with ",
    "'posterior_ordinate_attribute()' or 'posterior_ordinate_append()'."
  )
  expect_error(
    BayesTools:::.posterior_ordinate_from_attribute(
      list(value = 0, ordinate = .25, method = "qCMDE"),
      null_hypothesis = 0
    ),
    raw_message,
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.posterior_ordinate_from_attribute(
      list(
        ordinates = data.frame(value = c(0, .5), ordinate = c(.25, .5)),
        method    = "qCMDE"
      ),
      null_hypothesis = .5
    ),
    raw_message,
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.posterior_ordinate_from_attribute(
      data.frame(x = 0, y = .25),
      null_hypothesis = 0
    ),
    raw_message,
    fixed = TRUE
  )

  # ambiguous or non-positive ordinates cannot be constructed
  expect_error(
    .posterior_ordinate_for_test(value = c(0, 0), ordinate = c(.25, .30)),
    "unique",
    fixed = TRUE
  )
  expect_error(
    .posterior_ordinate_for_test(value = 0, ordinate = 0),
    "finite and positive",
    fixed = TRUE
  )

  # a modified attribute that no longer satisfies the schema is rejected
  modified <- .posterior_ordinate_for_test(value = 0, ordinate = .25)
  modified$value <- c(0, 0)
  expect_error(
    BayesTools:::.posterior_ordinate_from_attribute(modified, 0),
    "Posterior ordinate metadata is invalid",
    fixed = TRUE
  )
})

test_that("stored posterior ordinate parser keeps diagnostics aligned", {

  parsed <- BayesTools:::.posterior_ordinate_from_attribute(
    .posterior_ordinate_for_test(
      value       = c(1, 0),
      ordinate    = c(.25, .50),
      diagnostics = list(
        relative_mcse    = c(.9, .1),
        BF_error_percent = c(90, 10),
        estimator        = "q_grid_cmde"
      )
    ),
    null_hypothesis = 0
  )

  expect_equal(parsed$x, 0)
  expect_equal(parsed$y, .50)
  expect_equal(parsed$diagnostics$relative_mcse, .1)
  expect_equal(parsed$diagnostics$BF_error_percent, 10)
  expect_equal(parsed$diagnostics$estimator, "q_grid_cmde")
})

test_that("stored posterior ordinate attachment preserves multiple null values", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- BayesTools:::.posterior_ordinate_attach(
    samples = posterior,
    sources = list(list(
      .posterior_ordinate_for_test(parameter = "theta", value = 0,  ordinate = .25),
      .posterior_ordinate_for_test(parameter = "theta", value = .5, ordinate = .50)
    )),
    parameter = "theta"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, .5) / .50
  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = .5,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
})

test_that("Savage_Dickey_BF ignores stored posterior density by default", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(-4, 4, length.out = 401),
    y      = rep(100, 401),
    method = "iwmde"
  ))

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) /
    BayesTools:::.Savage_Dickey_BF.kd(posterior, 0)

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE
  )
  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "posterior_density_source"), "KDE")
})

test_that("Savage_Dickey_BF does not infer exact support from bounded prior grids", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )

  # the null 0 lies below the draws: an unreflected kernel-tail ordinate
  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) /
    as.numeric(BayesTools:::.Savage_Dickey_BF.kd(posterior, 0, warn_extrapolation = FALSE))

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE
  )

  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "posterior_density_source"), "KDE")
  expect_null(attr(out, "posterior_density_boundary_reflection", exact = TRUE))
})

test_that("Savage_Dickey_BF uses exact posterior support for KDE fallback", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", posterior_support_attribute(c(0, 1)))

  posterior_height <- BayesTools:::.Savage_Dickey_BF.kd(posterior, 0)
  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) /
    as.numeric(posterior_height)

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE
  )

  expect_true(attr(posterior_height, "boundary_reflection"))
  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "posterior_density_source"), "KDE")
  expect_true(attr(out, "posterior_density_boundary_reflection"))
  expect_equal(attr(out, "posterior_density_support"), c(0, 1))

  # support metadata must come from the constructor
  for(raw_support in list(c(0, 1), list(lower = 0, upper = 1, exact = TRUE),
                          list(bounds = c(0, 1)))){
    expect_error(
      .bt_meta_set(posterior, "support", raw_support),
      "Posterior support metadata must be created with 'posterior_support_attribute()'.",
      fixed = TRUE
    )
  }
})

test_that("Savage_Dickey_BF rejects stored density support when grid misses the null", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x       = seq(.25, .75, length.out = 101),
    y       = rep(1, 101),
    method  = "iwmde",
    support = posterior_support_attribute(c(0, 1))
  ))

  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      density_method       = "precomputed"
    ),
    "Stored posterior density does not span",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF lets exact support override stale precomputed heights", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(0, 1, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", posterior_support_attribute(c(0, 1)))
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value    = -.5,
    ordinate = .5,
    method   = "qCMDE"
  ))

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = -.5,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )
  expect_equal(as.numeric(out), Inf)
  expect_equal(attr(out, "posterior_density_source"), "exact_support_exclusion")
  expect_true(attr(out, "posterior_density_fallback"))
  expect_match(
    attr(out, "posterior_density_fallback_warnings"),
    "posterior support excludes the null hypothesis"
  )

  posterior <- .bt_meta_set(posterior, "posterior_ordinate", NULL)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x       = seq(-1, 1, length.out = 101),
    y       = rep(.5, 101),
    method  = "qCMDE",
    support = posterior_support_attribute(c(0, 1))
  ))
  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = -.5,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )
  expect_equal(as.numeric(out), Inf)
  expect_equal(attr(out, "posterior_density_source"), "exact_support_exclusion")
  expect_true(attr(out, "posterior_density_fallback"))
  expect_equal(attr(out, "posterior_density_support"), c(0, 1))
  expect_match(
    attr(out, "posterior_density_fallback_warnings"),
    "stored posterior density support excludes the null hypothesis"
  )
})

test_that("Savage_Dickey_BF uses matched density support to reject stale ordinates", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(0, 1, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value    = -.5,
    ordinate = .5,
    method   = "stale-ordinate"
  ))
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x       = seq(0, 1, length.out = 101),
    y       = rep(.5, 101),
    method  = "support-guard",
    support = posterior_support_attribute(c(0, 1))
  ))

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = -.5,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out), Inf)
  expect_equal(attr(out, "posterior_density_source"), "exact_support_exclusion")
  expect_true(attr(out, "posterior_density_fallback"))
  expect_equal(attr(out, "posterior_density_support"), c(0, 1))
  expect_match(
    attr(out, "posterior_density_fallback_warnings"),
    "Ignoring the precomputed posterior ordinate",
    fixed = TRUE
  )
  expect_match(
    attr(out, "posterior_density_fallback_warnings"),
    "stored posterior density support excludes the null hypothesis"
  )
})

test_that("Savage_Dickey_BF ignores non-exact or incompatible support metadata", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )

  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(0, 1), exact = FALSE, source = "test"))
  expected <- BayesTools:::.prior_linear_density_height(prior_density, .5) /
    BayesTools:::.Savage_Dickey_BF.kd(posterior, .5)
  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = .5,
    normal_approximation = FALSE,
    silent               = TRUE
  )
  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_null(attr(out, "posterior_density_boundary_reflection", exact = TRUE))

  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(.2, .8), source = "test"))
  expect_warning(
    out <- Savage_Dickey_BF(
      posterior,
      null_hypothesis      = .5,
      normal_approximation = FALSE
    ),
    "Exact posterior support metadata is incompatible",
    fixed = TRUE
  )
  expect_null(attr(out, "posterior_density_boundary_reflection", exact = TRUE))
})

test_that("posterior support distinguishes point support from interval support", {

  point_support <- BayesTools:::.posterior_support_new(
    c(0, 1),
    points = c(0, 1),
    type   = "points"
  )
  expect_false(BayesTools:::.posterior_support_contains_value(point_support, .5))
  expect_true(BayesTools:::.posterior_support_contains_value(point_support, 1))

  point_support_info <- BayesTools:::.posterior_support_for_kde(
    c(0, 1, 0, 1),
    support = point_support
  )
  expect_null(point_support_info[["bounds"]])
  expect_match(
    point_support_info[["warning"]],
    "continuous interval",
    fixed = TRUE
  )

  mixed_support <- BayesTools:::.posterior_support_new(
    c(0, 1),
    points = 0,
    type   = "mixed"
  )
  expect_true(BayesTools:::.posterior_support_contains_value(mixed_support, .5))
})

test_that("Savage_Dickey_BF reports incompatible support on precomputed paths", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(.2, .8), source = "test"))
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", posterior_ordinate_attribute(
    value          = .5,
    ordinate       = 1,
    method         = "qCMDE",
    density_method = "precomputed"
  ))

  expect_warning(
    out <- Savage_Dickey_BF(
      posterior,
      null_hypothesis = .5,
      density_method  = "precomputed"
    ),
    "Exact posterior support metadata is incompatible",
    fixed = TRUE
  )
  expect_equal(attr(out, "posterior_density_source"), "precomputed")
  expect_match(
    attr(out, "warnings"),
    "Exact posterior support metadata is incompatible",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF ignores incompatible support before excluding nulls", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-1, 1, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(0, 1), source = "test"))
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", posterior_ordinate_attribute(
    value          = -.5,
    ordinate       = .5,
    method         = "qCMDE",
    density_method = "precomputed"
  ))

  expect_warning(
    out <- Savage_Dickey_BF(
      posterior,
      null_hypothesis = -.5,
      density_method  = "precomputed"
    ),
    "Exact posterior support metadata is incompatible",
    fixed = TRUE
  )
  expect_true(is.finite(out))
  expect_equal(attr(out, "posterior_density_source"), "precomputed")
})

test_that("marginal_posterior propagates exact scalar support from mixed samples", {

  theta_prior <- prior("beta", list(alpha = 1, beta = 1))
  theta <- seq(.001, .999, length.out = 101)
  class(theta) <- c("mixed_posteriors", "mixed_posteriors.simple", class(theta))
  theta <- .bt_meta_set(theta, "draw_index", seq_along(theta))
  theta <- .bt_draws_set_component(theta, source = "model", component = rep(1, length(theta)))
  attr(theta, "parameter") <- "theta"
  attr(theta, "prior_list") <- theta_prior
  theta <- BayesTools:::.posterior_support_set_from_prior_list(theta, theta_prior)

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")
  samples <- .bt_meta_set(samples, "prior_context", BayesTools:::.prior_density_build_context(
    prior_list   = list(theta = theta_prior),
    column_names = "theta",
    n_grid       = 1024
  ))

  marginal <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = TRUE
  )
  marginal_no_prior <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = FALSE
  )
  transformed <- marginal_posterior(
    samples,
    parameter                = "theta",
    prior_samples            = TRUE,
    transformation           = "lin",
    transformation_arguments = list(a = 1, b = 2)
  )

  expect_equal(BayesTools:::.posterior_support_bounds(marginal), c(0, 1))
  expect_equal(BayesTools:::.posterior_support_bounds(marginal_no_prior), c(0, 1))
  expect_equal(BayesTools:::.posterior_support_bounds(transformed), c(1, 3))

  theta <- .bt_meta_set(theta, "support", BayesTools:::.posterior_support_new(c(.2, .8), source = "posterior"))
  samples[["theta"]] <- theta
  marginal_existing_support <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = TRUE
  )
  expect_equal(
    BayesTools:::.posterior_support_bounds(marginal_existing_support),
    c(.2, .8)
  )

  support <- BayesTools:::.posterior_support_new(c(0, 1))
  expect_null(BayesTools:::.posterior_support_transform(
    support,
    "lin",
    list(a = 1, b = 0)
  ))
  expect_null(BayesTools:::.posterior_support_transform(
    support,
    "exp_lin",
    list(a = 1, b = 0)
  ))

  connected_support <- BayesTools:::.posterior_support_union(list(
    posterior_support_attribute(c(0, 1)),
    posterior_support_attribute(c(1, 2))
  ))
  disjoint_support <- BayesTools:::.posterior_support_union(list(
    posterior_support_attribute(c(0, 1)),
    posterior_support_attribute(c(2, 3))
  ))
  expect_true(connected_support$exact)
  expect_false(disjoint_support$exact)
})

test_that("marginal_posterior preserves an attached simple prior density", {

  theta <- seq(.1, .9, length.out = 101)
  class(theta) <- c("mixed_posteriors", "mixed_posteriors.simple", class(theta))
  theta <- .bt_meta_set(theta, "draw_index", seq_along(theta))
  theta <- .bt_draws_set_component(theta, source = "model", component = rep(1, length(theta)))
  attr(theta, "parameter")  <- "theta"
  attr(theta, "prior_list") <- prior_none()
  stored_prior <- prior("uniform", list(a = 0, b = 1))
  theta <- .bt_meta_set(theta, "prior_density", stored_prior)

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")

  marginal <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = TRUE
  )
  transformed <- marginal_posterior(
    samples,
    parameter                = "theta",
    prior_samples            = TRUE,
    transformation           = "lin",
    transformation_arguments = list(a = 1, b = 2)
  )

  expect_identical(
    .bt_meta_get(marginal, "prior_density"),
    stored_prior
  )
  expect_false(identical(
    .bt_meta_get(transformed, "prior_density"),
    stored_prior
  ))
})

test_that("marginal_posterior infers support from the current prior context", {

  raw_prior <- prior("beta", list(alpha = 1, beta = 1))
  transformed_prior <- prior("uniform", list(10, 20))
  theta <- seq(11, 19, length.out = 51)
  class(theta) <- c("mixed_posteriors", "mixed_posteriors.simple", class(theta))
  theta <- .bt_meta_set(theta, "draw_index", seq_along(theta))
  theta <- .bt_draws_set_component(theta, source = "model", component = rep(1, length(theta)))
  attr(theta, "parameter") <- "theta"
  attr(theta, "prior_list") <- raw_prior

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")
  samples <- .bt_meta_set(samples, "transform_scaled", TRUE)
  samples <- .bt_meta_set(samples, "prior_context", BayesTools:::.prior_density_context(
    prior_list   = list(theta = transformed_prior),
    column_names = "theta",
    n_grid       = 64
  ))

  marginal <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = FALSE
  )

  expect_equal(BayesTools:::.posterior_support_bounds(marginal), c(10, 20))
})

test_that("marginal_posterior rebuilds conditional context for support", {

  theta_prior <- prior_spike_and_slab(
    prior("uniform", list(10, 20)),
    prior_inclusion = prior("point", list(location = .5))
  )
  prior_list <- list(theta = theta_prior)
  condition_event <- BayesTools:::.condition_event(
    prior_list        = prior_list,
    conditional       = "theta",
    conditional_rule  = "AND"
  )
  theta <- seq(11, 19, length.out = 51)
  class(theta) <- c("mixed_posteriors", "mixed_posteriors.simple", class(theta))
  theta <- .bt_meta_set(theta, "draw_index", seq_along(theta))
  theta <- .bt_draws_set_component(theta, source = "model", component = rep(1, length(theta)))
  attr(theta, "parameter") <- "theta"
  attr(theta, "prior_list") <- theta_prior
  theta <- BayesTools:::.posterior_support_set_from_prior_list(theta, theta_prior)
  theta <- BayesTools:::.condition_event_set_attributes(theta, condition_event)

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")
  samples <- .bt_meta_set(samples, "prior_context", BayesTools:::.prior_density_context(
    prior_list   = prior_list,
    column_names = "theta",
    n_grid       = 64
  ))

  marginal <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = FALSE
  )

  expect_equal(BayesTools:::.posterior_support_bounds(marginal), c(10, 20))
  expect_equal(.bt_meta_condition(marginal, "conditional"), "theta")
  expect_equal(.bt_meta_condition(marginal, "condition_key"), condition_event[["condition_key"]])
})

test_that("support propagation ignores zero-weight and stale raw components", {

  null_prior <- BayesTools:::.set_prior_model_weight(
    prior("point", list(location = 0)),
    0
  )
  slab_prior <- BayesTools:::.set_prior_model_weight(
    prior("uniform", list(10, 20)),
    1
  )
  support <- BayesTools:::.posterior_support_from_prior_list(
    list(null_prior, slab_prior)
  )

  expect_equal(support$bounds, c(10, 20))

  mixture_prior <- prior_mixture(
    list(
      prior("point", list(location = 0)),
      prior("uniform", list(10, 20))
    ),
    is_null = c(TRUE, FALSE)
  )
  attr(mixture_prior, "prior_weights") <- c(0, 1)
  mixture_context <- BayesTools:::.prior_density_context(
    prior_list   = list(theta = mixture_prior),
    column_names = "theta"
  )
  mixture_support <- BayesTools:::.posterior_support_from_prior_context_weights(
    mixture_context,
    c(theta = 1)
  )

  expect_equal(mixture_support$bounds, c(10, 20))

  zero_weight_null <- BayesTools:::.set_prior_model_weight(prior_none(), 0)
  weightfunction_prior <- BayesTools:::.set_prior_model_weight(
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .4))),
    1
  )
  omega_context <- list(
    names   = c("omega[0,0.05]", "omega[0.05,1]"),
    mapping = list(c(NA_integer_, NA_integer_), c(1L, 2L))
  )
  omega_support <- BayesTools:::.posterior_support_weightfunction_columns(
    list(zero_weight_null, weightfunction_prior),
    omega_context
  )

  expect_equal(omega_support[["omega[0.05,1]"]]$bounds, c(.4, .4))

  samples <- matrix(1, nrow = 2, ncol = 2)
  colnames(samples) <- c("theta", "display_theta")
  samples <- .bt_meta_set(samples, "support", list(
    theta = BayesTools:::.posterior_support_new(
      c(0, 20),
      source = "prior_list"
    ),
    display_theta = BayesTools:::.posterior_support_new(
      c(0, 20),
      source = "prior_list"
    ),
    posterior_only = BayesTools:::.posterior_support_new(
      c(-1, 1),
      source = "posterior"
    )
  ))
  context <- BayesTools:::.prior_density_context(
    prior_list   = list(theta = slab_prior),
    column_names = "theta"
  )

  refreshed <- BayesTools:::.posterior_support_set_from_prior_context(
    samples,
    context
  )
  refreshed_support <- .bt_meta_get(refreshed, "support")

  expect_equal(BayesTools:::.posterior_support_bounds(refreshed, "theta"), c(10, 20))
  expect_null(refreshed_support[["display_theta"]])
  expect_equal(
    BayesTools:::.posterior_support_bounds(refreshed, "posterior_only"),
    c(-1, 1)
  )
})

test_that("formula marginal support is propagated without prior densities", {

  theta_prior <- prior(
    "beta",
    list(alpha = 1, beta = 1),
    truncation = list(lower = .2, upper = .8)
  )
  theta <- seq(.25, .75, length.out = 101)
  class(theta) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "mixed_posteriors.formula",
    class(theta)
  )
  theta <- .bt_meta_set(theta, "draw_index", seq_along(theta))
  theta <- .bt_draws_set_component(theta, source = "model", component = rep(1, length(theta)))
  attr(theta, "parameter") <- "mu_x"
  theta <- .bt_meta_set(theta, "formula_parameter", "mu")
  attr(theta, "prior_list") <- theta_prior

  samples <- list(mu_x = theta)
  class(samples) <- c("mixed_posteriors", "list")

  marginal_no_prior <- marginal_posterior(
    samples,
    parameter     = "mu_x",
    formula       = y ~ 0 + x,
    prior_samples = FALSE
  )
  marginal_with_prior <- marginal_posterior(
    samples,
    parameter     = "mu_x",
    formula       = y ~ 0 + x,
    prior_samples = TRUE
  )

  expected_bounds <- list(
    "-1SD" = c(-.8, -.2),
    "0SD"  = c(0, 0),
    "1SD"  = c(.2, .8)
  )
  for(level in names(expected_bounds)){
    expect_equal(
      BayesTools:::.posterior_support_bounds(marginal_no_prior[[level]]),
      expected_bounds[[level]]
    )
    expect_equal(
      BayesTools:::.posterior_support_bounds(marginal_with_prior[[level]]),
      expected_bounds[[level]]
    )
  }
})

test_that("formula marginal_posterior attaches matched top-level precomputed metadata", {

  mu_intercept <- rep(0, 51)
  class(mu_intercept) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "mixed_posteriors.formula",
    class(mu_intercept)
  )
  mu_intercept <- .bt_meta_set(mu_intercept, "draw_index", seq_along(mu_intercept))
  mu_intercept <- .bt_draws_set_component(mu_intercept, source = "model", component = rep(1, length(mu_intercept)))
  attr(mu_intercept, "parameter") <- "mu_intercept"
  mu_intercept <- .bt_meta_set(mu_intercept, "formula_parameter", "mu")
  attr(mu_intercept, "prior_list") <- prior("normal", list(0, 1))

  mu_x <- seq(-1, 1, length.out = 51)
  class(mu_x) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "mixed_posteriors.formula",
    class(mu_x)
  )
  mu_x <- .bt_meta_set(mu_x, "draw_index", seq_along(mu_x))
  mu_x <- .bt_draws_set_component(mu_x, source = "model", component = rep(1, length(mu_x)))
  attr(mu_x, "parameter") <- "mu_x"
  mu_x <- .bt_meta_set(mu_x, "formula_parameter", "mu")
  attr(mu_x, "prior_list") <- prior("normal", list(0, 1))

  samples <- list(mu_intercept = mu_intercept, mu_x = mu_x)
  class(samples) <- c("mixed_posteriors", "list")
  samples <- .bt_meta_set(samples, "posterior_density", list(
    one_sd = posterior_density_attribute(
      x         = seq(-1, 1, length.out = 101),
      y         = rep(.5, 101),
      method    = "formula-density",
      density_method = "precomputed",
      parameter = "mu_x[1SD]"
    )
  ))

  marginal <- marginal_posterior(
    samples,
    parameter     = "mu_x",
    formula       = y ~ x,
    prior_samples = FALSE
  )
  transformed <- marginal_posterior(
    samples,
    parameter                = "mu_x",
    formula                  = y ~ x,
    prior_samples            = FALSE,
    transformation           = "lin",
    transformation_arguments = list(a = 0, b = 2)
  )

  expect_equal(
    .bt_meta_get(marginal[["1SD"]], "posterior_density")[["method"]],
    "formula-density"
  )
  expect_null(.bt_meta_get(marginal[["0SD"]], "posterior_density"))
  # the stored density follows the transformation 2 x: locations doubled,
  # heights halved (the Jacobian)
  transformed_density <- .bt_meta_get(transformed[["1SD"]], "posterior_density")
  expect_equal(transformed_density$x, 2 * seq(-1, 1, length.out = 101))
  expect_equal(transformed_density$y, rep(.25, 101))
  expect_identical(transformed_density$method, "formula-density")
})

test_that("spike-and-slab posterior constructors attach support metadata", {

  slab <- prior(
    "beta",
    list(alpha = 2, beta = 2),
    truncation = list(lower = .2, upper = .8)
  )
  spike_slab <- prior_spike_and_slab(
    slab,
    prior_inclusion = prior("point", list(location = .5))
  )

  posterior_samples <- BayesTools:::.as_mixed_posteriors.spike_and_slab(
    cbind(theta = seq(.2, .8, length.out = 100),
          theta_indicator = rep(c(0, 1), 50)),
    spike_slab,
    parameter = "theta"
  )
  posterior_support <- BayesTools:::.posterior_support_get(posterior_samples)

  expect_equal(posterior_support$bounds, c(0, .8))
  expect_true(0 %in% posterior_support$points)
})

test_that("simplex posterior support uses component and convex-hull bounds", {

  simplex_prior <- prior("dirichlet", list(alpha = c(2, 3, 5)))

  component_support <- BayesTools:::.posterior_support_from_prior(simplex_prior)
  expect_equal(component_support$bounds, c(0, 1))
  expect_true(component_support$exact)

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(w = simplex_prior),
    column_names = paste0("w[", 1:3, "]")
  )

  all_weight_support <- BayesTools:::.posterior_support_from_prior_context_weights(
    context,
    c("w[1]" = 2, "w[2]" = 5, "w[3]" = 7)
  )
  expect_equal(all_weight_support$bounds, c(2, 7))

  partial_weight_support <- BayesTools:::.posterior_support_from_prior_context_weights(
    context,
    c("w[1]" = 2, "w[3]" = 5)
  )
  expect_equal(partial_weight_support$bounds, c(0, 5))
})

test_that("vector point support remains exact point support", {

  vector_point_prior <- prior("mpoint", list(location = 2, K = 3))

  component_support <- BayesTools:::.posterior_support_from_prior(vector_point_prior)
  expect_equal(component_support$bounds, c(2, 2))
  expect_equal(component_support$points, 2)
  expect_equal(component_support$type, "points")

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(theta = vector_point_prior),
    column_names = paste0("theta[", 1:3, "]")
  )
  linear_support <- BayesTools:::.posterior_support_from_prior_context_weights(
    context,
    c("theta[1]" = 2, "theta[2]" = -0.5)
  )

  expect_equal(linear_support$bounds, c(3, 3))
  expect_equal(linear_support$points, 3)
  expect_equal(linear_support$type, "points")
})

test_that("top-level marginal metadata does not replace child-specific metadata", {

  child <- 1:10
  class(child) <- c("marginal_posterior.simple", class(child))
  attr(child, "level_name") <- "A"
  child <- .bt_meta_set(child, "posterior_density", posterior_density_attribute(
    x              = seq(-1, 1, length.out = 11),
    y              = rep(.5, 11),
    method         = "child-density",
    density_method = "precomputed",
    parameter      = "theta[A]"
  ))

  marginal <- list(A = child)
  class(marginal) <- c("marginal_posterior.factor", "list")
  samples <- .bt_meta_update(
    structure(list(), class = c("mixed_posteriors", "list")),
    posterior_density = list(
      posterior_density_attribute(
        x              = seq(-1, 1, length.out = 11),
        y              = rep(.5, 11),
        method         = "top-density",
        density_method = "precomputed",
        parameter      = "theta[A]"
      )
    )
  )

  out <- BayesTools:::.marginal_posterior_attach_precomputed_metadata(
    marginal  = marginal,
    samples   = samples,
    parameter = "theta"
  )

  expect_equal(
    .bt_meta_get(out[["A"]], "posterior_density")[["method"]],
    "child-density"
  )
})

test_that("Savage_Dickey_BF rejects stored density missing the null", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-1, 1, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(.5, 1, length.out = 101),
    y      = rep(1, 101),
    method = "iwmde"
  ))

  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      density_method       = "precomputed"
    ),
    "Stored posterior density does not span",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF diagnoses invalid precomputed metadata", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-1, 1, length.out = 301),
    prior_density
  )

  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis = 0,
      density_method  = "precomputed",
      silent          = TRUE
    ),
    "requires valid posterior ordinate or posterior density metadata",
    fixed = TRUE
  )

  # attributes modified into an invalid state stop with the reason
  invalid_density <- .posterior_density_for_test(
    x      = c(-1, 1),
    y      = c(1, 1),
    method = "invalid-density"
  )
  invalid_density$x <- 0
  invalid_density$y <- 1
  posterior <- .bt_meta_set(posterior, "posterior_density", invalid_density)
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis = 0,
      density_method  = "precomputed"
    ),
    "Posterior density metadata is invalid: the density grid needs at least two points",
    fixed = TRUE
  )

  # raw lists are rejected when they are attached
  expect_error(
    .bt_meta_set(posterior, "posterior_density", list(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(1, 101),
      method = "raw-list"
    )),
    "Posterior density metadata must be created with 'posterior_density_attribute()'.",
    fixed = TRUE
  )

  posterior <- .bt_meta_set(posterior, "posterior_density", NULL)
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value    = 1,
    ordinate = .5,
    method   = "wrong-null"
  ))
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis = 0,
      density_method  = "precomputed"
    ),
    "Precomputed posterior ordinate metadata is present but invalid",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF rejects stored density with zero null height", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-1, 1, length.out = 301),
    prior_density
  )
  stored_x <- seq(-1, 1, length.out = 101)
  stored_y <- abs(stored_x)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  ))

  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      density_method       = "precomputed"
    ),
    "zero or non-finite height"
  )
})

test_that("stored posterior density parser rejects degenerate grids", {

  parsed <- BayesTools:::.posterior_density_from_attribute(.posterior_density_for_test(
    x      = c(0, 0, 1),
    y      = c(1, 3, 2),
    method = "iwmde"
  ))
  expect_equal(parsed$x, c(0, 1))
  expect_equal(unname(parsed$y), c(2, 2))

  grid_message <- "Posterior density attribute must contain a valid positive density grid."
  expect_error(
    .posterior_density_for_test(x = c(0, 1), y = c(0, 0)),
    grid_message,
    fixed = TRUE
  )
  expect_error(
    .posterior_density_for_test(x = c(1, 1), y = c(1, 2)),
    grid_message,
    fixed = TRUE
  )

  # a grid made degenerate after construction is rejected, not dropped
  modified <- .posterior_density_for_test(x = c(0, 1), y = c(1, 2))
  modified$y <- c(0, 0)
  expect_error(
    BayesTools:::.posterior_density_from_attribute(modified),
    "Posterior density metadata is invalid",
    fixed = TRUE
  )

  # raw lists are not posterior density metadata
  expect_error(
    BayesTools:::.posterior_density_from_attribute(list(
      x = c(0, 1), y = c(1, 2), method = "iwmde"
    )),
    "Posterior density metadata must be created with 'posterior_density_attribute()'.",
    fixed = TRUE
  )
  for(raw_density in list(
    list(density = list(x = c(0, 1), y = c(1, 2)), estimator = "iwmde"),
    data.frame(x = c(0, 1), y = c(1, 2))
  )){
    expect_error(
      BayesTools:::.posterior_density_from_attribute(raw_density),
      "must be created with 'posterior_density_attribute()'",
      fixed = TRUE
    )
  }
  # containers of attributes carry no single density
  expect_null(BayesTools:::.posterior_density_from_attribute(list(
    theta = .posterior_density_for_test(x = c(0, 1), y = c(1, 2))
  )))
})

test_that("stored posterior density selection validates names and conditionals", {

  density <- .posterior_density_for_test(
    parameter        = "theta",
    conditional      = c("b", "a"),
    conditional_rule = "OR",
    x                = seq(-1, 1, length.out = 11),
    y                = rep(1, 11),
    method           = "iwmde"
  )
  sources <- list(
    theta = density,
    phi = modifyList(density, list(parameter = "phi"))
  )

  expect_equal(
    BayesTools:::.posterior_density_from_sources(
      sources          = list(sources),
      aliases          = "theta",
      conditional      = c("a", "b"),
      conditional_rule = "OR"
    )$method,
    "iwmde"
  )
  expect_null(BayesTools:::.posterior_density_from_sources(
    sources          = list(sources),
    aliases          = "theta",
    conditional      = c("a", "b"),
    conditional_rule = "AND"
  ))
  density[["condition_key"]] <- BayesTools:::.condition_event_key(c("a", "b"), "OR")
  sources <- list(theta = density)
  expect_equal(
    BayesTools:::.posterior_density_from_sources(
      sources          = list(sources),
      aliases          = "theta",
      conditional      = c("a", "b"),
      conditional_rule = "OR",
      condition_key    = BayesTools:::.condition_event_key(c("a", "b"), "OR")
    )$method,
    "iwmde"
  )
  expect_null(BayesTools:::.posterior_density_from_sources(
    sources          = list(sources),
    aliases          = "theta",
    conditional      = c("a", "b"),
    conditional_rule = "AND",
    condition_key    = BayesTools:::.condition_event_key(c("a", "b"), "AND")
  ))
  alias_density <- density
  alias_density[["condition_key"]] <- NULL
  alias_density[["conditional"]] <- "phacking"
  alias_density[["conditional_rule"]] <- "AND"
  expect_equal(
    BayesTools:::.posterior_density_from_sources(
      sources          = list(list(theta = alias_density)),
      aliases          = "theta",
      conditional      = "alpha",
      conditional_rule = "AND"
    )$method,
    "iwmde"
  )
  expect_null(BayesTools:::.posterior_density_from_sources(
    sources = list(sources),
    aliases = "missing"
  ))
})

test_that("stored posterior ordinate selection continues past empty alias branches", {

  alias_branch <- .posterior_ordinate_for_test(
    parameter = "phi",
    value     = 0,
    ordinate  = 100,
    method    = "wrong-parameter"
  )
  valid_branch <- .posterior_ordinate_for_test(
    parameter = "theta",
    value     = 0,
    ordinate  = .5,
    method    = "valid-ordinate"
  )

  matched <- BayesTools:::.posterior_ordinate_from_sources(
    sources = list(list(theta = alias_branch, fallback = valid_branch)),
    aliases = "theta"
  )

  expect_equal(matched[["method"]], "valid-ordinate")
  expect_equal(
    BayesTools:::.posterior_ordinate_from_attribute(matched, 0)[["y"]],
    .5
  )
})

test_that("posterior atom declarations validate and aggregate point masses", {

  atoms <- posterior_atom_attribute(list(
    x    = c(0, 0, 1),
    mass = c(.2, .3, .4)
  ))

  expect_equal(as.numeric(atoms$locations[, 1L]), c(0, 1))
  expect_equal(atoms$mass, c(.5, .4))

  for(point_masses in list(
    list(x = c(0, 1, 2), mass = c(.6, .5, 2)),
    list(location = c(0, 1), p = c(.2, .3)),
    data.frame(location = 0, mass = .2)
  )){
    expect_error(
      posterior_atom_attribute(point_masses),
      "Posterior 'point_masses' metadata is invalid.",
      fixed = TRUE
    )
  }
})

test_that("Savage_Dickey_BF ignores stored density range for normal approximation", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.5, 1, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(-1, 1, length.out = 101),
    y      = rep(1, 101),
    method = "iwmde"
  ))

  expect_warning(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis      = 0,
      normal_approximation = TRUE,
      density_method       = "precomputed"
    ),
    "Posterior samples do not span"
  )
})

test_that("Savage_Dickey_BF uses top-level stored density for list posteriors", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  stored_x <- seq(-4, 4, length.out = 401)
  stored_y <- stats::dnorm(stored_x, mean = .4, sd = 1.2)
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(
    level = .posterior_density_for_test(
      x      = stored_x,
      y      = stored_y,
      method = "iwmde"
    )
  ))

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) /
    stats::approx(stored_x, stored_y, xout = 0)[["y"]]
  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)
  expect_equal(attr(out[["level"]], "posterior_density_source"), "precomputed")

  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", NULL)
  posterior_list <- .bt_meta_set(posterior_list, "posterior_densities", list(
    list(
      level = .posterior_density_for_test(
        x      = stored_x,
        y      = stored_y,
        method = "iwmde"
      )
    )
  ))
  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)
  expect_equal(attr(out[["level"]], "posterior_density_source"), "precomputed")
})

test_that("Savage_Dickey_BF respects child conditionals for top-level sources", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(
    posterior, "condition", list(conditional = "theta", conditional_rule = "AND")
  )
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(
    level = .posterior_density_for_test(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(100, 101),
      method = "stale"
    ),
    matching = .posterior_density_for_test(
      parameter   = "level",
      conditional = "theta",
      x           = seq(-1, 1, length.out = 101),
      y           = rep(.50, 101),
      method      = "qCMDE"
    )
  ))

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50
  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)

  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", NULL)
  posterior_list <- .bt_meta_set(posterior_list, "posterior_ordinates", list(list(
    level = .posterior_ordinate_for_test(
      value    = 0,
      ordinate = 100,
      method   = "stale"
    ),
    matching = .posterior_ordinate_for_test(
      parameter   = "level",
      conditional = "theta",
      value       = 0,
      ordinate    = .50,
      method      = "qCMDE"
    )
  )))

  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)
})

test_that("Savage_Dickey_BF revalidates positional top-level sources", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(
    posterior, "condition", list(conditional = "theta", conditional_rule = "AND")
  )
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")

  mismatched_density <- .posterior_density_for_test(
    parameter   = "level",
    conditional = "phi",
    x           = seq(-1, 1, length.out = 101),
    y           = rep(100, 101),
    method      = "stale-density"
  )
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(mismatched_density))

  expect_null(BayesTools:::.posterior_density_child_attributes(posterior_list)[[1]])

  matching_density <- mismatched_density
  matching_density[["conditional"]] <- "theta"
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(matching_density))
  expect_equal(
    BayesTools:::.posterior_density_child_attributes(posterior_list)[[1]][["method"]],
    "stale-density"
  )

  mismatched_ordinate <- .posterior_ordinate_for_test(
    parameter   = "level",
    conditional = "phi",
    value       = 0,
    ordinate    = 100,
    method      = "stale-ordinate"
  )
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", NULL)
  posterior_list <- .bt_meta_set(posterior_list, "posterior_ordinate", list(mismatched_ordinate))

  expect_null(BayesTools:::.posterior_ordinate_child_attributes(
    posterior_list,
    null_hypothesis = 0
  )[[1]])

  matching_ordinate <- mismatched_ordinate
  matching_ordinate[["conditional"]] <- "theta"
  posterior_list <- .bt_meta_set(posterior_list, "posterior_ordinate", list(matching_ordinate))
  expect_equal(
    BayesTools:::.posterior_ordinate_child_attributes(
      posterior_list,
      null_hypothesis = 0
    )[[1]][["method"]],
    "stale-ordinate"
  )
})

test_that("Savage_Dickey_BF accepts multi-ordinate posterior attributes", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", posterior_ordinate_append(
    .posterior_ordinate_for_test(
      value       = -.50,
      ordinate    = 100,
      parameter   = "theta",
      diagnostics = list(relative_mcse = .10)
    ),
    .posterior_ordinate_for_test(
      value       = 0,
      ordinate    = .50,
      parameter   = "theta",
      diagnostics = list(relative_mcse = .20)
    )
  ))

  expect_true(BayesTools:::.posterior_ordinate_has_data(
    .bt_meta_get(posterior, "posterior_ordinate")
  ))

  matched_source <- BayesTools:::.posterior_ordinate_from_sources(
    sources         = list(.bt_meta_get(posterior, "posterior_ordinate")),
    aliases         = "theta",
    null_hypothesis = 0
  )
  expect_s3_class(matched_source, "BayesTools_posterior_ordinates")
  expect_length(matched_source$ordinates, 2L)
  expect_equal(vapply(matched_source$ordinates, `[[`, numeric(1), "value"), c(-.5, 0))

  parsed <- BayesTools:::.posterior_ordinate_from_attribute(matched_source, 0)
  expect_equal(parsed[["method"]], "qCMDE")
  expect_equal(parsed[["diagnostics"]][["relative_mcse"]], .20)

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis = 0,
    silent          = TRUE,
    density_method  = "precomputed"
  )
  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50

  expect_equal(as.numeric(out), expected, tolerance = 1e-12)
  expect_equal(attr(out, "BF_error_percent"), 20)
})

test_that("Savage_Dickey_BF gives valid child sources precedence over top-level sources", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(-1, 1, length.out = 101),
    y      = rep(.50, 101),
    method = "child-density"
  ))
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value    = 0,
    ordinate = .50,
    method   = "child-ordinate"
  ))
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(
    level = .posterior_density_for_test(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(100, 101),
      method = "top-density"
    )
  ))
  posterior_list <- .bt_meta_set(posterior_list, "posterior_ordinate", list(
    level = .posterior_ordinate_for_test(
      value    = 0,
      ordinate = 100,
      method   = "top-ordinate"
    )
  ))

  child_density <- BayesTools:::.posterior_density_child_attributes(posterior_list)[[1]]
  child_ordinate <- BayesTools:::.posterior_ordinate_child_attributes(
    posterior_list,
    null_hypothesis = 0
  )[[1]]
  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50
  expect_equal(child_density[["method"]], "child-density")
  expect_equal(child_ordinate[["method"]], "child-ordinate")
  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)
})

test_that("Savage_Dickey_BF replaces null-unusable child sources and rejects invalid ones", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value    = 1,
    ordinate = 100,
    method   = "wrong-null"
  ))
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  posterior_list <- .bt_meta_set(posterior_list, "posterior_ordinate", list(
    level = .posterior_ordinate_for_test(
      value    = 0,
      ordinate = .50,
      method   = "top-ordinate"
    )
  ))

  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50
  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)

  # an invalid child attribute is an error, not a reason to use another source
  invalid_density <- .posterior_density_for_test(
    x      = c(-1, 1),
    y      = c(100, 100),
    method = "invalid-child-density"
  )
  invalid_density$x <- 0
  invalid_density$y <- 100
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", NULL)
  posterior <- .bt_meta_set(posterior, "posterior_density", invalid_density)
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(
    level = .posterior_density_for_test(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(.50, 101),
      method = "top-density"
    )
  ))

  expect_error(
    Savage_Dickey_BF(
      posterior_list,
      null_hypothesis      = 0,
      normal_approximation = FALSE,
      silent               = TRUE,
      density_method       = "precomputed"
    ),
    "Posterior density metadata is invalid",
    fixed = TRUE
  )
})

test_that("Savage_Dickey_BF replaces null-unusable child density with top-level density", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(1, 2, length.out = 101),
    y      = rep(100, 101),
    method = "child-misses-null"
  ))
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  posterior_list <- .bt_meta_set(posterior_list, "posterior_density", list(
    level = .posterior_density_for_test(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(.50, 101),
      method = "top-density"
    )
  ))

  child_density <- BayesTools:::.posterior_density_child_attributes(
    posterior_list,
    null_hypothesis = 0
  )[[1]]
  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50
  expect_equal(child_density[["method"]], "top-density")
  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)
  expect_equal(attr(out[["level"]], "posterior_density_source"), "precomputed")
})

test_that("Savage_Dickey_BF uses declarations rather than posterior-null clusters", {

  continuous_prior <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior_cluster <- .marginal_posterior_with_prior_density_for_test(
    c(rep(0, 8), seq(-2, 2, length.out = 92)),
    continuous_prior
  )

  expect_silent(
    Savage_Dickey_BF(
      posterior_cluster,
      null_hypothesis = 0,
      normal_approximation = TRUE,
      silent = TRUE
    )
  )

  posterior_without_declaration <- posterior_cluster
  posterior_without_declaration <- .bt_meta_set(posterior_without_declaration, "atoms", NULL)
  expect_error(
    Savage_Dickey_BF(
      posterior_without_declaration,
      null_hypothesis = 0,
      normal_approximation = TRUE
    ),
    "explicit atom/no-atom declaration",
    fixed = TRUE
  )

  posterior_cluster <- .bt_meta_set(posterior_cluster, "atoms", posterior_atom_attribute(
    list(x = 0, mass = .08)
  ))
  expect_error(
    Savage_Dickey_BF(
      posterior_cluster,
      null_hypothesis = 0,
      normal_approximation = TRUE
    ),
    "declared point mass"
  )

  point_prior <- BayesTools:::.prior_linear_density_point(0)
  posterior_continuous <- .marginal_posterior_with_prior_density_for_test(
    seq(-2, 2, length.out = 101),
    point_prior
  )

  expect_error(
    Savage_Dickey_BF(posterior_continuous, null_hypothesis = 0, normal_approximation = TRUE),
    "point mass in the prior"
  )

  posterior_with_stored_point <- .marginal_posterior_with_prior_density_for_test(
    seq(-2, 2, length.out = 101),
    continuous_prior
  )
  posterior_with_stored_point <- .bt_meta_set(posterior_with_stored_point, "posterior_density", .posterior_density_for_test(
    x            = seq(-2, 2, length.out = 101),
    y            = .8 * stats::dnorm(seq(-2, 2, length.out = 101))
  ))
  posterior_with_stored_point <- .bt_meta_set(posterior_with_stored_point, "atoms", posterior_atom_attribute(list(x = 0, mass = .2)))

  expect_error(
    Savage_Dickey_BF(
      posterior_with_stored_point,
      null_hypothesis = 0,
      density_method = "precomputed"
    ),
    "declared point mass"
  )
})

test_that("Savage_Dickey_BF excludes off-null atoms from the continuous ordinate", {

  continuous_prior <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  continuous_draws <- stats::rnorm(800, mean = 0, sd = 1)
  spike_draws <- rep(2, 200)
  posterior <- .marginal_posterior_with_prior_density_for_test(
    c(continuous_draws, spike_draws),
    continuous_prior
  )
  posterior <- .bt_meta_set(posterior, "atoms", posterior_atom_attribute(
    list(x = 2, mass = 0.2)
  ))

  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis = 2,
      normal_approximation = TRUE,
      silent = TRUE
    ),
    "declared point mass"
  )

  with_atoms <- Savage_Dickey_BF(
    posterior,
    null_hypothesis = 0,
    normal_approximation = TRUE,
    silent = TRUE
  )
  continuous_only <- .marginal_posterior_with_prior_density_for_test(
    continuous_draws,
    continuous_prior
  )
  expected <- Savage_Dickey_BF(
    continuous_only,
    null_hypothesis = 0,
    normal_approximation = TRUE,
    silent = TRUE
  )
  # Continuous-only BF uses density of continuous draws; mixed BF scales that
  # continuous ordinate by continuous mass 0.8, so BF is larger by 1/0.8.
  expect_equal(as.numeric(with_atoms), as.numeric(expected) / 0.8, tolerance = 1e-10)

  continuous_info <- BayesTools:::.Savage_Dickey_BF.continuous_posterior(
    posterior,
    .bt_meta_get(posterior, "atoms")
  )
  expect_equal(continuous_info$continuous_mass, 0.8)
  expect_false(any(continuous_info$samples == 2))
})

test_that("Savage_Dickey_BF diagnoses zero prior density at point null", {

  # A zero prior density at the null (Beta(2, 2) at 0) makes the density
  # ratio a 0/0 limit that the posterior kernel estimate cannot estimate: a
  # scalar call stops with the classed condition of hypothesis_BF(), and a
  # level of a list posterior gets an NA Bayes factor with the reason.
  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 2, beta = 2))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(0, 1), source = "test"))

  condition <- tryCatch(Savage_Dickey_BF(posterior, null_hypothesis = 0, silent = TRUE),
                        error = function(e) e)
  expect_s3_class(condition, "BayesTools_zero_ordinate")
  expect_s3_class(condition, "BayesTools_hypothesis_ordinate")
  expect_identical(
    conditionMessage(condition),
    "Prior density at point hypothesis 'parameter = 0' is zero, so the Savage-Dickey density ratio is undefined."
  )

  levels <- list(A = posterior, B = posterior)
  class(levels) <- c("marginal_posterior", "list")
  out <- Savage_Dickey_BF(levels, null_hypothesis = 0, silent = TRUE)
  expect_true(is.na(out[["A"]]))
  expect_identical(
    attr(out[["A"]], "warnings"),
    "The prior density at the null hypothesis value is zero. The Savage-Dickey Bayes factor is undefined."
  )
  expect_identical(attr(out[["A"]], "posterior_density_source"), "zero_prior_ordinate")
})

test_that("Savage_Dickey_BF applies the exactness rule of point hypotheses", {

  # The prior ordinate at the null follows the rule of hypothesis_BF() point
  # hypotheses with its condition classes (each also of class
  # BayesTools_hypothesis_ordinate).
  set.seed(21)
  draws <- stats::rnorm(4000, .3, .5)
  bf_condition <- function(prior_density, null_hypothesis, samples = draws){
    tryCatch(
      Savage_Dickey_BF(.marginal_posterior_with_prior_density_for_test(samples, prior_density),
                       null_hypothesis = null_hypothesis, silent = TRUE),
      error = function(e) e
    )
  }
  expect_ordinate_class <- function(condition, class){
    expect_s3_class(condition, class)
    expect_s3_class(condition, "BayesTools_hypothesis_ordinate")
  }

  # three t terms have no structural route (a numerical grid height before)
  t3 <- prior("t", list(0, 1, 3))
  three_t <- BayesTools:::.prior_linear_combination_density(
    list(a = t3, b = t3, c = t3), c(a = 1, b = 1, c = 1)
  )
  expect_ordinate_class(bf_condition(three_t, 0), "BayesTools_inexact_ordinate")

  # N(0, 1e-3) * Cauchy(0, 1): its quadrature at .3 is rejected by its
  # diagnostics (an unclassed error before)
  priors <- list(beta = prior("normal", list(0, 1e-3)), sigma = prior("cauchy", list(0, 1)))
  attr(priors$beta, "multiply_by") <- "sigma"
  rejected <- bf_condition(BayesTools:::.prior_linear_combination_density(priors, c(beta = 1)), .3)
  expect_ordinate_class(rejected, "BayesTools_inexact_ordinate")
  expect_match(conditionMessage(rejected), "rejected by its diagnostics", fixed = TRUE)

  # a density grid without provenance (an unclassed error before)
  grid <- structure(list(density = list(x = c(-3, 0, 3), y = c(0, 1 / 3, 0), mass = 1),
                         points = data.frame(x = numeric(), p = numeric())),
                    class = c("prior_linear_density", "prior_density"))
  expect_ordinate_class(bf_condition(grid, 0), "BayesTools_inexact_ordinate")

  # A positive power preserves the zero atom independently of the height at 1.
  nonnegative_spike <- prior_spike_and_slab(prior("normal", list(0, 1), list(0, Inf)),
                                            prior_inclusion = prior("spike", list(.5)))
  undefined <- BayesTools:::.prior_linear_combination_density(
    list(x = nonnegative_spike), c(x = 1), output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 0, b = 2)
  )
  transformed_BF <- bf_condition(undefined, 1, abs(draws))
  expect_true(is.numeric(transformed_BF) && is.finite(transformed_BF) && transformed_BF > 0)
  expect_equal(exp(prior_density_ordinate(undefined, 1)$log_density), stats::dnorm(1) / 2, tolerance = 1e-12)

  # a prior point mass at the null, as before
  spike <- prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("spike", list(.5)))
  expect_ordinate_class(
    bf_condition(BayesTools:::.prior_linear_combination_density(list(x = spike), c(x = 1)), 0),
    "BayesTools_point_mass_at_null"
  )

  # the infinite gamma(1/2) density at 0 (Inf with a warning before)
  expect_ordinate_class(
    bf_condition(
      BayesTools:::.prior_linear_combination_density(list(x = prior("gamma", list(.5, 1))), c(x = 1)),
      0, abs(draws)
    ),
    "BayesTools_infinite_ordinate"
  )
})

test_that("plot_marginal uses stored posterior density when available", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 512
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    stats::rnorm(100, 0, 1),
    prior_density
  )
  stored_x <- seq(-2, 2, length.out = 51)
  stored_y <- stats::dnorm(stored_x, sd = .8)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  ))

  plot_data <- BayesTools:::.plot_data_marginal_samples(
    samples                  = list(theta = posterior),
    parameter                = "theta",
    prior                    = FALSE,
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  expect_equal(plot_data[["density1"]][["x"]], stored_x)
  expect_equal(plot_data[["density1"]][["y"]], stored_y)
  expect_equal(attr(plot_data[["density1"]], "posterior_density_method"), "iwmde")
})

test_that("plot_marginal draws declared atoms with a stored continuous density", {

  posterior <- .marginal_posterior_with_prior_density_for_test(
    c(rep(0, 25), seq(-2, 2, length.out = 75)),
    BayesTools:::.prior_linear_density_point(0, p = .25)
  )
  posterior <- .bt_meta_set(posterior, "atoms", posterior_atom_attribute(
    data.frame(x = 0, mass = .25)
  ))
  stored_x <- seq(-2, 2, length.out = 51)
  stored_y <- .75 * stats::dnorm(stored_x)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  ))

  # the stored density is the continuous part and the declared atom its
  # point mass
  expect_no_warning(
    plot_data <- BayesTools:::.plot_data_marginal_samples(
      samples                  = list(theta = posterior),
      parameter                = "theta",
      prior                    = FALSE,
      n_points                 = 16,
      transformation           = NULL,
      transformation_arguments = NULL,
      transformation_settings  = FALSE,
      density_method           = "precomputed"
    )
  )

  expect_equal(plot_data$density1$x, stored_x)
  expect_equal(plot_data$density1$y, stored_y)
  point_entries <- plot_data[vapply(plot_data, inherits, logical(1), "density.prior.point")]
  expect_equal(length(point_entries), 1L)
  expect_equal(point_entries[[1L]]$x, 0)
  expect_equal(point_entries[[1L]]$y, .25)
})

test_that("plot_marginal accepts posterior density diagnostics", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 512
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    c(rep(0, 20), stats::rnorm(80, 0, 1)),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "atoms", posterior_atom_attribute(list(x = 0, mass = .2)))
  stored_x <- seq(-2, 2, length.out = 51)
  stored_y <- .8 * stats::dnorm(stored_x, sd = .8)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    parameter    = "theta",
    x            = stored_x,
    y            = stored_y,
    method       = "iwmde",
    diagnostics  = list(min_ess = 40)
  ))

  plot_data <- BayesTools:::.plot_data_marginal_samples(
    samples                  = list(theta = posterior),
    parameter                = "theta",
    prior                    = FALSE,
    n_points                 = 16,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE,
    density_method           = "precomputed"
  )

  expect_equal(plot_data[["density1"]][["x"]], stored_x)
  expect_equal(plot_data[["density1"]][["y"]], stored_y)
  expect_equal(attr(plot_data[["density1"]], "posterior_density_method"), "iwmde")
  expect_equal(attr(plot_data[["density1"]], "posterior_density_diagnostics")$min_ess, 40)
  point_data <- plot_data[vapply(plot_data, inherits, logical(1), "density.prior.point")][[1]]
  expect_equal(point_data[["x"]], 0)
  expect_equal(point_data[["y"]], .2)
})

test_that("marginal_posterior handles direct multi-factor transformed interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  interaction_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(interaction_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")
  class(interaction_samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")
  attr(interaction_samples, "levels")            <- BayesTools:::.get_prior_factor_levels(interaction_prior)
  attr(interaction_samples, "level_names")       <- attr(interaction_prior, "level_names")
  attr(interaction_samples, "interaction")       <- TRUE
  attr(interaction_samples, "interaction_terms") <- attr(interaction_prior, "interaction_terms")
  attr(interaction_samples, "term_components")   <- attr(interaction_prior, "term_components")
  attr(interaction_samples, "factor_terms")      <- attr(interaction_prior, "factor_terms")
  attr(interaction_samples, "factor_contrasts")  <- attr(interaction_prior, "factor_contrasts")
  attr(interaction_samples, "factor_design")     <- attr(interaction_prior, "factor_design")
  attr(interaction_samples, "factor_cell_names") <- attr(interaction_prior, "factor_cell_names")
  attr(interaction_samples, "orthonormal")       <- TRUE
  attr(interaction_samples, "meandif")           <- FALSE
  attr(interaction_samples, "treatment")         <- FALSE
  attr(interaction_samples, "independent")       <- FALSE
  attr(interaction_samples, "prior_list")        <- interaction_prior

  samples <- list(mu_a__xXx__b = interaction_samples)
  class(samples) <- c("as_mixed_posteriors", "mixed_posteriors")

  marginal <- marginal_posterior(samples, "mu_a__xXx__b", use_formula = FALSE)
  expected <- interaction_samples %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(names(marginal), attr(interaction_prior, "factor_cell_names"))
  expect_equal(attr(marginal, "level_names"), attr(interaction_prior, "factor_cell_names"))
  expect_equal(as.numeric(marginal[[1]]), as.numeric(expected[, 1]))
  expect_equal(as.numeric(marginal[[6]]), as.numeric(expected[, 6]))

  marginal_with_prior <- marginal_posterior(
    samples = samples,
    parameter = "mu_a__xXx__b",
    prior_samples = TRUE,
    use_formula = FALSE,
    n_samples = 32
  )

  expect_equal(names(marginal_with_prior), attr(interaction_prior, "factor_cell_names"))
  expect_true(all(vapply(marginal_with_prior, function(x) {
    inherits(.bt_meta_get(x, "prior_density"), "prior_linear_density")
  }, logical(1))))
})

test_that("marginal_posterior uses transformed treatment metadata for simple factor priors", {

  df <- data.frame(
    x_fac2t = factor(c("A", "B", "A", "B"), levels = c("A", "B"))
  )
  formula_result <- JAGS_formula(
    formula = ~ x_fac2t,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x_fac2t   = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )

  posterior <- matrix(rnorm(20), ncol = 1)
  colnames(posterior) <- "mu_x_fac2t"
  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(
    fit,
    parameters = "mu_x_fac2t"
  )
  marginal <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_x_fac2t",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 32
  )

  expect_equal(names(marginal), c("A", "B"))
  expect_true(inherits(.bt_meta_get(marginal[["A"]], "prior_density"), "prior_linear_density"))
  expect_true(inherits(.bt_meta_get(marginal[["B"]], "prior_density"), "prior_linear_density"))
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(.bt_meta_get(marginal[["A"]], "prior_density"), 0),
    1
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(.bt_meta_get(marginal[["B"]], "prior_density"), 0),
    0
  )
})

test_that("marginal_posterior handles as_mixed_posteriors multi-factor interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_a__xXx__b[", 1:2, "]")
  fit <- coda::mcmc(complete_test_posterior(posterior, formula_result$prior_list))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(fit, parameters = "mu_a__xXx__b")
  marginal <- marginal_posterior(
    samples = samples,
    parameter = "mu_a__xXx__b",
    prior_samples = TRUE,
    use_formula = FALSE,
    n_samples = 32
  )
  expected <- posterior %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(names(marginal), attr(interaction_prior, "factor_cell_names"))
  expect_equal(as.numeric(marginal[[1]]), as.numeric(expected[, 1]))
  expect_equal(as.numeric(marginal[[6]]), as.numeric(expected[, 6]))
  expect_true(all(vapply(marginal, function(x) {
    inherits(.bt_meta_get(x, "prior_density"), "prior_linear_density")
  }, logical(1))))
})

test_that("marginal_posterior handles one-coefficient as_mixed_posteriors interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2"), levels = c("b1", "b2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  posterior <- matrix(seq_len(10), nrow = 10, ncol = 1)
  colnames(posterior) <- "mu_a__xXx__b"
  fit <- coda::mcmc(complete_test_posterior(posterior, formula_result$prior_list))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(fit, parameters = "mu_a__xXx__b")
  marginal <- marginal_posterior(samples, "mu_a__xXx__b", use_formula = FALSE)
  expected <- posterior %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(names(marginal), attr(interaction_prior, "factor_cell_names"))
  expect_equal(as.numeric(marginal[[1]]), as.numeric(expected[, 1]))
  expect_equal(as.numeric(marginal[[4]]), as.numeric(expected[, 4]))
})

test_that("marginal_posterior aligns selected formula levels and at expansions", {

  fixture <- .marginal_semantic_fixture_for_test()

  marginal <- marginal_posterior(
    samples       = fixture$samples,
    parameter     = "mu_x_fac3md",
    at            = list(x_cont1 = 1, x_fac2t = c("A", "B")),
    formula       = ~ x_cont1 + x_fac2t + x_cont1 * x_fac3md,
    prior_samples = TRUE,
    n_samples     = 128
  )

  expect_s3_class(marginal, "marginal_posterior")
  expect_equal(attr(marginal, "parameter"), "mu_x_fac3md")
  expect_equal(names(marginal), c("A", "B", "C"))
  expect_equal(attr(marginal, "level_names"), c("A", "B", "C"))
  expected_level_at <- data.frame(x_fac3md = factor(c("A", "B", "C"), levels = c("A", "B", "C")))
  actual_level_at <- attr(marginal, "level_at")
  attr(actual_level_at, "out.attrs") <- NULL
  expect_equal(actual_level_at, expected_level_at)

  expected_meandif <- fixture$posterior[, c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  expected_interaction <- fixture$posterior[, c(
    "mu_x_cont1__xXx__x_fac3md[1]",
    "mu_x_cont1__xXx__x_fac3md[2]"
  )] %*% t(contr.meandif(1:3))

  for(i in seq_along(marginal)){
    expected_a <- fixture$posterior[, "mu_intercept"] +
      fixture$posterior[, "mu_x_cont1"] +
      expected_meandif[, i] +
      expected_interaction[, i]
    expected_b <- expected_a + fixture$posterior[, "mu_x_fac2t"]

    expect_equal(dim(marginal[[i]]), c(2, nrow(fixture$posterior)))
    expect_equal(as.numeric(marginal[[i]][1, ]), as.numeric(expected_a))
    expect_equal(as.numeric(marginal[[i]][2, ]), as.numeric(expected_b))
    expect_equal(attr(marginal[[i]], "data")$x_cont1, c(1, 1))
    expect_equal(as.character(attr(marginal[[i]], "data")$x_fac2t), c("A", "B"))
    expect_equal(as.character(attr(marginal[[i]], "data")$x_fac3md), rep(names(marginal)[i], 2))
    expect_s3_class(.bt_meta_get(marginal[[i]], "prior_density"), "prior_linear_density")
  }
})

test_that("marginal_posterior transformation preserves level names and transformed values", {

  fixture <- .marginal_semantic_fixture_for_test()

  marginal_raw <- marginal_posterior(
    samples       = fixture$samples,
    parameter     = "mu_x_cont1",
    formula       = ~ x_cont1 + x_fac2t + x_cont1 * x_fac3md,
    prior_samples = TRUE,
    n_samples     = 128
  )
  marginal_exp <- marginal_posterior(
    samples        = fixture$samples,
    parameter      = "mu_x_cont1",
    formula        = ~ x_cont1 + x_fac2t + x_cont1 * x_fac3md,
    transformation = "exp",
    prior_samples  = TRUE,
    n_samples      = 128
  )

  expect_equal(names(marginal_exp), names(marginal_raw))
  expect_equal(attr(marginal_exp, "level_at"), attr(marginal_raw, "level_at"))
  for(level in names(marginal_raw)){
    expect_equal(as.numeric(marginal_exp[[level]]), exp(as.numeric(marginal_raw[[level]])))

    exp_prior_data <- BayesTools:::.prior_linear_density_to_plot_data(
      .bt_meta_get(marginal_exp[[level]], "prior_density"),
      n_points = 128
    )
    expect_true(all(vapply(exp_prior_data, function(d) all(d$x > 0), logical(1))))
  }
})

test_that("marginal posterior plot data separates densities, point masses, labels, and selected parameters", {

  fixture <- .marginal_semantic_fixture_for_test()
  marginal <- marginal_posterior(
    samples       = fixture$samples,
    parameter     = "mu_x_fac2t",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 128
  )
  samples <- list(mu_x_fac2t = marginal)

  expect_error(suppressWarnings(plot_marginal(samples, parameter = "not_here")))

  plot_data <- BayesTools:::.plot_data_marginal_samples(
    samples,
    parameter = "mu_x_fac2t",
    prior = TRUE,
    n_points = 128,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )

  expect_equal(sort(unique(vapply(plot_data, attr, character(1), which = "level_name"))), c("A", "B"))

  point_data <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.point")]
  density_data <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.simple")]
  expect_length(point_data, 1)
  expect_equal(attr(point_data[[1]], "level_name"), "A")
  expect_equal(point_data[[1]]$x, 0)
  expect_equal(point_data[[1]]$y, 1)
  expect_length(density_data, 1)
  expect_equal(attr(density_data[[1]], "level_name"), "B")
  .expect_density_area_for_test(density_data[[1]])

  prior_data <- unlist(lapply(plot_data, attr, which = "prior"), recursive = FALSE)
  expect_true(any(vapply(prior_data, inherits, logical(1), what = "density.prior.point")))
  expect_true(any(vapply(prior_data, inherits, logical(1), what = "density.prior.simple")))

  expect_error(
    plot_marginal(samples, parameter = "mu_x_fac2t", prior = TRUE, transformation = "not-a-transformation"),
    "not recognized|must be"
  )
  expect_error(
    plot_marginal(list(mu_x_fac2t = marginal_posterior(fixture$samples, "mu_x_fac2t", use_formula = FALSE)),
                  parameter = "mu_x_fac2t", prior = TRUE),
    "'samples' did not contain prior densities"
  )
  expect_error(
    plot_marginal(samples, parameter = "mu_x_fac2t", plot_type = "lattice"),
    "'plot_type'"
  )
})

test_that("marginal density plot data preserves probability mass by level", {

  fixture <- .marginal_semantic_fixture_for_test()

  factor_marginal <- marginal_posterior(
    samples       = fixture$samples,
    parameter     = "mu_x_fac2t",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 512
  )
  factor_plot_data <- BayesTools:::.plot_data_marginal_samples(
    list(mu_x_fac2t = factor_marginal),
    parameter = "mu_x_fac2t",
    prior = TRUE,
    n_points = 512,
    transformation = NULL,
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )
  factor_prior_data <- unlist(lapply(factor_plot_data, attr, which = "prior"), recursive = FALSE)

  expect_equal(sort(names(.expect_density_mass_by_level_for_test(factor_plot_data))), c("A", "B"))
  expect_equal(sort(names(.expect_density_mass_by_level_for_test(factor_prior_data))), c("A", "B"))

  conditional_marginal <- marginal_posterior(
    samples       = fixture$samples,
    parameter     = "mu_x_fac3md",
    at            = list(x_cont1 = 1, x_fac2t = c("A", "B")),
    formula       = ~ x_cont1 + x_fac2t + x_cont1 * x_fac3md,
    prior_samples = TRUE,
    n_samples     = 512
  )
  conditional_prior_data <- unlist(lapply(names(conditional_marginal), function(level_name) {
    BayesTools:::.prior_linear_density_to_plot_data(
      .bt_meta_get(conditional_marginal[[level_name]], "prior_density"),
      n_points = 512,
      factor = TRUE,
      level_name = level_name
    )
  }), recursive = FALSE)

  expect_equal(sort(names(.expect_density_mass_by_level_for_test(conditional_prior_data))), c("A", "B", "C"))

  transformed_marginal <- .marginal_posterior_with_prior_density_for_test(
    stats::qnorm(seq(0.001, 0.999, length.out = 1000)),
    BayesTools:::.prior_linear_combination_density(
      prior_list = list(theta = prior("normal", list(0, 1))),
      weights    = c(theta = 1),
      n_grid     = 2048
    )
  )
  transformed_plot_data <- BayesTools:::.plot_data_marginal_samples(
    list(theta = transformed_marginal),
    parameter = "theta",
    prior = TRUE,
    n_points = 512,
    transformation = "exp",
    transformation_arguments = NULL,
    transformation_settings = FALSE
  )
  transformed_prior_data <- unlist(lapply(transformed_plot_data, attr, which = "prior"), recursive = FALSE)

  .expect_density_mass_by_level_for_test(transformed_plot_data, tolerance = 0.15)
  .expect_density_mass_by_level_for_test(transformed_prior_data, tolerance = 0.15)
})

test_that("ggplot marginal output carries selected levels in rendered data", {

  skip_if_not_installed("ggplot2")

  fixture <- .marginal_semantic_fixture_for_test()
  marginal <- marginal_posterior(
    samples       = fixture$samples,
    parameter     = "mu_x_fac2t",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 128
  )

  plot <- plot_marginal(
    list(mu_x_fac2t = marginal),
    parameter = "mu_x_fac2t",
    plot_type = "ggplot",
    prior = TRUE,
    par_name = "fac2t",
    n_points = 128
  )
  built <- ggplot2::ggplot_build(plot)
  layer_data <- do.call(rbind, lapply(built$data, function(x) {
    if(all(c("x", "y") %in% names(x))) x[, intersect(c("x", "y", "colour", "linetype", "group"), names(x)), drop = FALSE]
  }))

  expect_s3_class(plot, "ggplot")
  expect_true(length(built$data) >= 2)
  expect_true(any(abs(layer_data$x) < sqrt(.Machine$double.eps)))
  expect_true(length(unique(stats::na.omit(layer_data$colour))) >= 2)
})

test_that("marginal_estimates_table reports exact level summaries and Bayes factors", {

  samples <- list(
    theta = list(
      low  = c(1, 2, 3, 4),
      high = c(10, 20, 30, 40)
    )
  )
  inference <- list(theta = list(low = 2, high = 0.5))

  table <- marginal_estimates_table(
    samples    = samples,
    inference  = inference,
    parameters = "theta",
    probs      = c(0.25, 0.5, 0.75)
  )

  expected_samples <- do.call(rbind, lapply(samples$theta, function(x) {
    c(
      Mean = mean(x),
      SD = stats::sd(x),
      "0.25" = unname(stats::quantile(x, 0.25)),
      "0.5" = unname(stats::quantile(x, 0.5)),
      "0.75" = unname(stats::quantile(x, 0.75))
    )
  }))

  expect_s3_class(table, "BayesTools_table")
  expect_equal(rownames(table), c("theta[low]", "theta[high]"))
  expect_equal(colnames(table), c("Mean", "SD", "0.25", "0.5", "0.75", "inclusion_BF"))
  expect_equal(unname(as.matrix(table[, 1:5])), unname(expected_samples), tolerance = 1e-12)
  expect_equal(as.numeric(table$inclusion_BF), c(2, 0.5), tolerance = 1e-12)
  expect_equal(attr(table, "type"), c(rep("estimate", 5), "inclusion_BF"))
  expect_true(attr(table, "rownames"))

  log_table <- marginal_estimates_table(
    samples    = samples,
    inference  = inference,
    parameters = "theta",
    probs      = c(0.25),
    logBF      = TRUE,
    BF01       = TRUE
  )
  expect_equal(as.numeric(log_table$inclusion_BF), log(c(1 / 2, 1 / 0.5)), tolerance = 1e-12)
  expect_equal(attr(log_table$inclusion_BF, "name"), "log(Exclusion BF)")
})

test_that("marginal_estimates_table keeps Bayes factor warnings of scalar parameters", {

  # Savage_Dickey_BF() returns a bare numeric for a scalar marginal posterior;
  # its warnings must reach the table like those of formula levels
  bf_theta <- structure(4, warnings = "Scalar warning.", BF_error_percent = 2)
  bf_level <- structure(.5, warnings = "Level warning.")
  table <- marginal_estimates_table(
    samples    = list(theta = c(1, 2, 3, 4), mu = list(intercept = c(2, 3, 4, 5)),
                      gamma = list(A = c(0, 1, 2, 3))),
    inference  = list(theta = bf_theta, mu = list(intercept = bf_level),
                      gamma = list(A = bf_level)),
    parameters = c("theta", "mu", "gamma")
  )

  expect_equal(as.numeric(table$inclusion_BF), c(4, .5, .5))
  expect_equal(as.numeric(table$BF_error_percent), c(2, NA, NA))
  expect_equal(
    attr(table, "warnings"),
    c("theta: Scalar warning.", "mu: Level warning.", "gamma[A]: Level warning.")
  )
})

# ============================================================================ #
# SECTION: Marginal posterior regressions (review round 3)
# ============================================================================ #

.mock_mixing_fit_for_marginal <- function(samples, prior_list){
  samples <- coda::mcmc(as.matrix(samples))
  fit <- structure(
    list(
      mcmc = coda::mcmc.list(samples),
      sample = nrow(samples),
      summary.pars = list(mutate = NULL),
      monitor = colnames(samples)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- prior_list
  attach_test_parameter_map(fit)
}

.prior_height_for_test <- function(x, value){
  ordinate <- prior_density_ordinate(.bt_meta_get(x, "prior_density"), value)
  exp(ordinate$log_density)
}

# Grid-based prior heights are checked against the adaptive evaluation's own
# documented error bound; exact ordinates must match to rounding error.
.expect_prior_height_for_test <- function(x, value, expected){
  height <- BayesTools:::.prior_linear_density_height(.bt_meta_get(x, "prior_density"), value)
  error_bound <- attr(height, "adaptive_evaluation")$error_bound
  if(is.null(error_bound)){
    error_bound <- 1e-8 * max(1, abs(expected))
  }
  expect_lte(abs(as.numeric(height) - expected), error_bound)
}

.collect_warnings_for_test <- function(expr){
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w){
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = messages)
}

test_that("Savage-Dickey BFs with the null outside the posterior draws warn once per parameter and level", {

  extrapolation <- paste0(
    "Posterior samples do not span both sides of the null hypothesis. ",
    "The posterior density at the null hypothesis is an extrapolation from ",
    "Gaussian kernel tails; the Bayes factor is not reliable evidence."
  )

  # model-averaged simple parameters: 'mu' lies above the null, 'tau' spans it
  set.seed(11)
  n <- 1000
  draws <- cbind(mu = 1.5 + .25 * stats::rnorm(n), tau = stats::rnorm(n))
  prior_list <- list(mu = prior("normal", list(0, 1)), tau = prior("normal", list(0, 1)))
  models <- lapply(1:2, function(i) list(
    fit = .mock_mixing_fit_for_marginal(draws, prior_list),
    marglik = bridgesampling_object(0), prior_weights = 1
  ))
  inference <- .collect_warnings_for_test(marginal_inference(
    models, marginal_parameters = c("mu", "tau"), parameters = c("mu", "tau"),
    is_null_list = list(mu = c(FALSE, FALSE), tau = c(FALSE, FALSE)),
    formula = NULL, n_samples = n, seed = 1
  ))
  expect_gt(min(inference$value$conditional$mu), 0)
  expect_identical(inference$warnings, paste0("mu: ", extrapolation))
  bf_mu <- inference$value$inference$mu
  expect_true(is.finite(bf_mu) && bf_mu > 1)
  expect_identical(attr(bf_mu, "warnings"), extrapolation)
  expect_null(attr(inference$value$inference$tau, "warnings"))
  table <- marginal_estimates_table(
    inference$value$conditional, inference$value$inference, c("mu", "tau")
  )
  expect_identical(attr(table, "warnings"), paste0("mu: ", extrapolation))

  # a null within the draws is an ordinary KDE ordinate
  expect_no_warning(
    bf_inside <- Savage_Dickey_BF(inference$value$conditional$mu, null_hypothesis = 1.5)
  )
  expect_null(attr(bf_inside, "warnings"))

  # formula levels: one warning per level, labelled as in the summary table
  formula_result <- JAGS_formula(
    ~ x, "mu", data = data.frame(x = c(-1, 0, 1)),
    prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)))
  )
  posterior <- cbind(
    mu_intercept = 1.5 + .2 * stats::rnorm(n),
    mu_x         = .3 + .05 * stats::rnorm(n)
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)
  level_inference <- .collect_warnings_for_test(as_marginal_inference(
    fit, marginal_parameters = "mu_x", parameters = c("mu_intercept", "mu_x"),
    conditional_list = list(mu_x = NULL), conditional_rule = "AND",
    formula = ~ x, n_samples = n
  ))
  # warnings name the rendered level labels of the marginal means
  level_warnings <- paste0("(mu) x[", c("-1SD", "0SD", "1SD"), "]: ", extrapolation)
  expect_identical(level_inference$warnings, level_warnings)
  level_BFs <- unlist(level_inference$value$inference$mu_x)
  expect_true(all(is.finite(level_BFs) & level_BFs > 1))
  level_table <- marginal_estimates_table(
    level_inference$value$conditional, level_inference$value$inference, "mu_x"
  )
  expect_identical(attr(level_table, "warnings"), level_warnings)

  # the list method labels the same levels
  direct <- .collect_warnings_for_test(
    Savage_Dickey_BF(level_inference$value$conditional$mu_x)
  )
  expect_identical(direct$warnings, level_warnings)
  expect_identical(unlist(direct$value), level_BFs)
})

# Mean and variance of a Gaussian KDE ordinate at x (reflected at a finite
# lower support bound) for n draws resampled with replacement from n_source
# draws of 'density', by quadrature: the reference for the per-component
# Savage-Dickey ordinates.
.reflected_kde_moments_for_test <- function(x, h, density, lower = -Inf, n, n_source){

  kernel <- function(y){
    k <- stats::dnorm(x, mean = y, sd = h)
    if(is.finite(lower)){
      k <- k + stats::dnorm(x, mean = 2 * lower - y, sd = h)
    }
    k
  }
  centres <- if(is.finite(lower)) c(x, 2 * lower - x) else x
  from <- max(lower, min(centres) - 12 * h)
  to <- max(centres) + 12 * h
  first <- stats::integrate(function(y) kernel(y) * density(y), from, to,
                            rel.tol = 1e-10, subdivisions = 500L)$value
  second <- stats::integrate(function(y) kernel(y)^2 * density(y), from, to,
                             rel.tol = 1e-10, subdivisions = 500L)$value
  # resampling n of n_source iid draws: Var(mean) = sigma^2 (1/n + 1/n_source)
  c(mean = first, variance = (second - first^2) * (1 / n + 1 / n_source))
}

.component_mixture_for_test <- function(n_source, seed){

  set.seed(seed)
  bounded <- prior("normal", list(.5, 1), list(0, Inf))
  data <- data.frame(t = factor(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi")))
  prior_lists <- lapply(list(prior("normal", list(0, 1)), bounded), function(intercept){
    JAGS_formula(~ 1 + t, "mu", data = data, prior_list = list(
      intercept = intercept,
      t         = prior_factor("normal", list(0, 1), contrast = "treatment")
    ))$prior_list
  })
  # prior-only draws with equal marginal likelihoods: the posterior is the prior
  draws <- list(
    cbind(mu_intercept = stats::rnorm(n_source)),
    cbind(mu_intercept = rng(bounded, n_source))
  )
  draws <- lapply(draws, function(x){
    cbind(x, "mu_t[1]" = stats::rnorm(n_source), "mu_t[2]" = stats::rnorm(n_source))
  })
  models <- lapply(1:2, function(i) list(
    fit = .mock_mixing_fit_for_marginal(draws[[i]], prior_lists[[i]]),
    marglik = bridgesampling_object(0), prior_weights = 1
  ))
  list(models = models, bounded = bounded)
}

test_that("Savage-Dickey mixes per-model ordinates when model supports differ", {

  n_source <- 20000
  n_mixed  <- 20000
  fixture <- .component_mixture_for_test(n_source, seed = 5)
  mixed <- mix_posteriors(
    fixture$models, parameters = c("mu_intercept", "mu_t"),
    is_null_list = list(mu_intercept = c(FALSE, FALSE), mu_t = c(FALSE, FALSE)),
    seed = 1, n_samples = n_mixed
  )
  simple <- marginal_posterior(mixed, "mu_intercept", use_formula = FALSE, prior_samples = TRUE)
  levels <- marginal_posterior(mixed, "mu_t", formula = ~ 1 + t, prior_samples = TRUE)
  levels <- lapply(levels, function(level){
    class(level) <- c(class(level), "marginal_posterior")
    level
  })
  models_ind <- .bt_meta_get(mixed$mu_intercept, "component")
  intercept <- as.numeric(mixed$mu_intercept)

  # model 1: N(0, 1) on the real line; model 2: N(0.5, 1) truncated to [0, Inf)
  densities <- list(
    function(y) stats::dnorm(y),
    function(y) stats::dnorm(y, .5, 1) / stats::pnorm(0, .5, 1, lower.tail = FALSE)
  )
  lower <- c(-Inf, 0)
  expect_equal(
    lapply(.bt_meta_get(simple, "components")$supports, `[[`, "bounds"),
    list(c(-Inf, Inf), c(0, Inf))
  )

  for(null in c(0, .001, .05, -.5)){
    # the true posterior (and prior) ordinate, one-sided on model 2's bound
    truth <- .5 * densities[[1]](null) + .5 * if(null >= 0) densities[[2]](null) else 0

    # reference: sum_m (n_m / n) E[f_m(null)], its Monte Carlo variance
    moments <- vapply(1:2, function(m){
      draws_m <- intercept[models_ind == m]
      if(null < lower[m]) return(c(mean = 0, variance = 0))
      .reflected_kde_moments_for_test(
        null, h = stats::bw.nrd0(draws_m), density = densities[[m]],
        lower = lower[m], n = length(draws_m), n_source = n_source
      )
    }, numeric(2))
    shares <- tabulate(models_ind, 2) / length(models_ind)
    expected <- sum(shares * moments["mean", ])
    mc_sd <- sqrt(sum(shares^2 * moments["variance", ]))

    for(marginal in list(simple, levels[["lo"]])){
      bf <- Savage_Dickey_BF(marginal, null_hypothesis = null)
      prior_height <- BayesTools:::.prior_linear_density_height(.bt_meta_get(marginal, "prior_density"), null)
      expect_equal(as.numeric(prior_height), truth, tolerance = 1e-6)
      posterior_height <- as.numeric(prior_height) / as.numeric(bf)
      # the ordinate is the per-model KDE mixture (4 Monte Carlo SD)
      expect_lte(abs(posterior_height - expected), 4 * mc_sd)
      # BF = 1 up to the reflected-KDE smoothing bias and 4 Monte Carlo SD
      # (at null 0: bias 1.6%, SD 2.4%; the pooled KDE gave BF 1.37, 11.7 SD off)
      expect_lte(abs(log(as.numeric(bf))), abs(log(expected / truth)) + 4 * mc_sd / expected)
      components <- attr(bf, "posterior_density_components")
      expect_equal(components$n, tabulate(models_ind, 2))
      expect_true(all(components$ordinate[lower > null] == 0))
    }
  }

  # direct point hypotheses use the same ordinate
  bf <- Savage_Dickey_BF(simple, null_hypothesis = .001)
  hypothesis <- hypothesis_BF(simple, hypothesis = "mu_intercept = 0.001", columns = "all")
  expect_equal(
    hypothesis$posterior,
    as.numeric(BayesTools:::.prior_linear_density_height(.bt_meta_get(simple, "prior_density"), .001)) /
      as.numeric(bf),
    tolerance = 1e-12
  )

  # levels whose supports agree across models keep the pooled estimate
  for(level in c("mid", "hi")){
    pooled <- levels[[level]]
    pooled <- .bt_meta_set(pooled, "components", NULL)
    bf <- Savage_Dickey_BF(levels[[level]], null_hypothesis = .05)
    expect_identical(bf, Savage_Dickey_BF(pooled, null_hypothesis = .05))
    expect_null(attr(bf, "posterior_density_components"))
  }
})

test_that("Savage-Dickey keeps the pooled ordinate for shared or single-model supports", {

  set.seed(8)
  n <- 4000
  bounded_1 <- prior("normal", list(0, 1), list(0, Inf))
  bounded_2 <- prior("normal", list(.5, 1), list(0, Inf))
  pooled_BF <- function(models){
    mixed <- mix_posteriors(
      models, parameters = "mu", is_null_list = list(mu = rep(FALSE, length(models))),
      seed = 1, n_samples = n
    )
    marginal <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
    expect_s3_class(.bt_meta_get(marginal, "components"), "BayesTools_posterior_components")
    stripped <- marginal
    stripped <- .bt_meta_set(stripped, "components", NULL)
    for(null in c(0, .3)){
      bf <- Savage_Dickey_BF(marginal, null_hypothesis = null, silent = TRUE)
      expect_identical(bf, Savage_Dickey_BF(stripped, null_hypothesis = null, silent = TRUE))
      expect_null(attr(bf, "posterior_density_components"))
    }
  }

  # both models truncated to [0, Inf)
  pooled_BF(list(
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu = rng(bounded_1, n)), list(mu = bounded_1)),
         marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu = rng(bounded_2, n)), list(mu = bounded_2)),
         marglik = bridgesampling_object(0), prior_weights = 1)
  ))
  # a point-null model contributes atoms only: one continuous component
  pooled_BF(list(
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu = rep(.5, n)), list(mu = prior("spike", list(.5)))),
         marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu = rng(bounded_2, n)), list(mu = bounded_2)),
         marglik = bridgesampling_object(0), prior_weights = 1)
  ))
})

test_that("Savage-Dickey extrapolation warnings use the components supporting the null", {

  n <- 4000
  unbounded <- prior("normal", list(0, 1))
  shifted <- prior("normal", list(6, 1), list(5, Inf))
  models <- list(
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu = stats::qnorm(stats::ppoints(n))), list(mu = unbounded)),
         marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu = 5 + abs(stats::qnorm(stats::ppoints(n)))), list(mu = shifted)),
         marglik = bridgesampling_object(0), prior_weights = 1)
  )
  mixed <- mix_posteriors(
    models, parameters = "mu", is_null_list = list(mu = c(FALSE, FALSE)),
    seed = 1, n_samples = n
  )
  marginal <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  draws <- as.numeric(marginal)
  models_ind <- .bt_meta_get(mixed$mu, "component")
  expect_gt(min(draws[models_ind == 2]), 5)

  # 4.5 lies within the pooled draws, but only model 1 supports it and its
  # draws end below 4.5: the ordinate is model 1's kernel tail
  expect_lt(max(draws[models_ind == 1]), 4.5)
  warnings <- .collect_warnings_for_test(Savage_Dickey_BF(marginal, null_hypothesis = 4.5))
  expect_match(warnings$warnings, "^mu: Posterior samples do not span both sides", all = TRUE)
  expect_length(warnings$warnings, 1L)
  expect_true(is.finite(warnings$value) && warnings$value > 0)
  expect_equal(attr(warnings$value, "posterior_density_components")$ordinate[2], 0)

  # a null spanned by the supporting component's draws is an ordinary ordinate
  expect_no_warning(bf <- Savage_Dickey_BF(marginal, null_hypothesis = -.5))
  expect_null(attr(bf, "warnings"))
})

.single_fit_for_test <- function(posterior, prior_list){

  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- prior_list
  fit <- attach_test_parameter_map(fit)
  fit
}

# Checks a Savage-Dickey BF of prior-only draws (true BF 1) against the
# quadrature reference of the per-component KDE ordinate: the ordinate lies
# within 4 Monte Carlo SD of its expectation, and the BF within the KDE
# smoothing bias plus 4 Monte Carlo SD (plus the prior ordinate's numerical
# error) of 1. 'components' lists the density, lower support bound and draws of
# each component; the draws are used without resampling.
.expect_component_ordinate_for_test <- function(bf, marginal, null, truth, components){

  n_total <- sum(vapply(components, function(component) length(component$draws), numeric(1)))
  moments <- vapply(components, function(component){
    if(null < component$lower){
      return(c(mean = 0, variance = 0))
    }
    .reflected_kde_moments_for_test(
      null, h = stats::bw.nrd0(component$draws), density = component$density,
      lower = component$lower, n = length(component$draws), n_source = Inf
    )
  }, numeric(2))
  shares <- vapply(components, function(component) length(component$draws), numeric(1)) / n_total
  expected <- sum(shares * moments["mean", ])
  mc_sd <- sqrt(sum(shares^2 * moments["variance", ]))

  prior_height <- as.numeric(BayesTools:::.prior_linear_density_height(.bt_meta_get(marginal, "prior_density"), null))
  prior_error <- abs(prior_height / truth - 1)
  expect_lt(prior_error, 1e-4)
  posterior_height <- prior_height / as.numeric(bf)
  expect_lte(abs(posterior_height - expected), 4 * mc_sd)
  expect_lte(abs(log(as.numeric(bf))),
             abs(log(expected / truth)) + 4 * mc_sd / expected + prior_error)
}

test_that("Savage-Dickey mixes per-component ordinates of single-fit mixture priors", {

  set.seed(12)
  n <- 20000
  unbounded <- prior("normal", list(0, 1))
  bounded   <- prior("normal", list(.5, 1), list(0, Inf))
  density_bounded <- function(y) ifelse(y >= 0, stats::dnorm(y, .5, 1) / stats::pnorm(.5), 0)

  # simple parameter: mu ~ mixture(N(0, 1), N(0.5, 1)T[0, Inf)), prior-only draws
  indicator <- sample(1:2, n, TRUE)
  mu <- ifelse(indicator == 1L, rng(unbounded, n), rng(bounded, n))
  mixed <- as_mixed_posteriors(
    .single_fit_for_test(cbind(mu = mu, mu_indicator = indicator),
                         list(mu = prior_mixture(list(unbounded, bounded), is_null = c(FALSE, FALSE)))),
    parameters = "mu"
  )
  marginal <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  components <- .bt_meta_get(marginal, "components")
  expect_setequal(components$keys[, "mu"], c(1, 2))
  expect_equal(
    lapply(components$supports, `[[`, "bounds")[match(c(1, 2), components$keys[, "mu"])],
    list(c(-Inf, Inf), c(0, Inf))
  )
  for(null in c(0, .001, .05, -.5)){
    bf <- Savage_Dickey_BF(marginal, null_hypothesis = null)
    .expect_component_ordinate_for_test(
      bf, marginal, null,
      truth = .5 * stats::dnorm(null) + .5 * density_bounded(null),
      components = list(
        list(density = stats::dnorm, lower = -Inf, draws = mu[indicator == 1L]),
        list(density = density_bounded, lower = 0, draws = mu[indicator == 2L])
      )
    )
  }

  # formula level mu = intercept + x with intercept ~ mixture(N(0, 1),
  # N(0.5, 1)T[0, Inf)) and x ~ mixture(spike(0), N(0, 1)): four components,
  # one of them (bounded intercept, spike slope) on [0, Inf)
  formula_result <- JAGS_formula(
    ~ x, "mu", data = data.frame(x = c(-1, 0, 1)),
    prior_list = list(
      intercept = prior_mixture(list(unbounded, bounded), is_null = c(FALSE, FALSE)),
      x         = prior_mixture(list(prior("spike", list(0)), unbounded), is_null = c(TRUE, FALSE))
    )
  )
  intercept_indicator <- sample(1:2, n, TRUE)
  slope_indicator     <- sample(1:2, n, TRUE)
  intercept <- ifelse(intercept_indicator == 1L, rng(unbounded, n), rng(bounded, n))
  slope     <- ifelse(slope_indicator == 1L, 0, stats::rnorm(n))
  mixed <- as_mixed_posteriors(
    .single_fit_for_test(
      cbind(mu_intercept = intercept, mu_x = slope,
            mu_intercept_indicator = intercept_indicator, mu_x_indicator = slope_indicator),
      formula_result$prior_list
    ),
    parameters = c("mu_intercept", "mu_x")
  )
  levels <- marginal_posterior(mixed, "mu_x", formula = ~ x, prior_samples = TRUE)
  level <- levels[["1SD"]]
  class(level) <- c(class(level), "marginal_posterior")
  keys <- .bt_meta_get(level, "components")$keys
  expect_setequal(paste(keys[, "mu_intercept"], keys[, "mu_x"]), c("1 1", "1 2", "2 1", "2 2"))

  # bounded intercept + N(0, 1) slope: N(0.5, sqrt(2)) * P(T >= 0 | T + Z)
  density_sum <- function(q){
    stats::dnorm(q, .5, sqrt(2)) * stats::pnorm(sqrt(2) * (q + .5) / 2) / stats::pnorm(.5)
  }
  tuple <- paste(intercept_indicator, slope_indicator)
  level_draws <- intercept + slope
  # the prior ordinate at the density jump itself (null 0) is unavailable
  # (adaptive grid refinement cannot converge across the jump)
  for(null in c(.001, .05, -.5)){
    bf <- Savage_Dickey_BF(level, null_hypothesis = null)
    .expect_component_ordinate_for_test(
      bf, level, null,
      truth = .25 * (stats::dnorm(null) + stats::dnorm(null, 0, sqrt(2)) +
                       density_bounded(null) + density_sum(null)),
      components = list(
        list(density = stats::dnorm, lower = -Inf, draws = level_draws[tuple == "1 1"]),
        list(density = function(y) stats::dnorm(y, 0, sqrt(2)), lower = -Inf, draws = level_draws[tuple == "1 2"]),
        list(density = density_bounded, lower = 0, draws = level_draws[tuple == "2 1"]),
        list(density = density_sum, lower = -Inf, draws = level_draws[tuple == "2 2"])
      )
    )
  }
})

test_that("Savage-Dickey keeps the pooled ordinate of single-fit mixtures with shared supports", {

  set.seed(13)
  n <- 4000
  expect_pooled <- function(marginal, nulls = c(0, .3)){
    stripped <- marginal
    stripped <- .bt_meta_set(stripped, "components", NULL)
    for(null in nulls){
      bf <- Savage_Dickey_BF(marginal, null_hypothesis = null, silent = TRUE)
      expect_identical(bf, Savage_Dickey_BF(stripped, null_hypothesis = null, silent = TRUE))
      expect_null(attr(bf, "posterior_density_components"))
    }
  }

  # spike-and-slab: the spike draws are atoms, one continuous component
  included <- stats::rbinom(n, 1, .5)
  spike_and_slab <- as_mixed_posteriors(
    .single_fit_for_test(cbind(mu = included * stats::rnorm(n, .2), mu_indicator = included),
                         list(mu = prior_spike_and_slab(prior("normal", list(0, 1))))),
    parameters = "mu"
  )
  marginal <- marginal_posterior(spike_and_slab, "mu", prior_samples = TRUE)
  expect_s3_class(.bt_meta_get(marginal, "components"), "BayesTools_posterior_components")
  expect_pooled(marginal, nulls = c(.3, -.4))

  # RoBMA-like effect: mixture of a null spike and a normal alternative
  component <- sample(1:2, n, TRUE)
  robma_like <- as_mixed_posteriors(
    .single_fit_for_test(
      cbind(mu = ifelse(component == 1L, 0, stats::rnorm(n, .2)), mu_indicator = component),
      list(mu = prior_mixture(list(prior("spike", list(0)), prior("normal", list(0, 1))),
                              is_null = c(TRUE, FALSE)))
    ),
    parameters = "mu"
  )
  expect_pooled(marginal_posterior(robma_like, "mu", prior_samples = TRUE), nulls = c(.3, -.4))

  # formula levels whose components all live on the real line
  formula_result <- JAGS_formula(
    ~ x, "mu", data = data.frame(x = c(-1, 0, 1)),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior_spike_and_slab(prior("normal", list(0, 1)))
    )
  )
  levels <- marginal_posterior(
    as_mixed_posteriors(
      .single_fit_for_test(
        cbind(mu_intercept = stats::rnorm(n), mu_x = included * stats::rnorm(n), mu_x_indicator = included),
        formula_result$prior_list
      ),
      parameters = c("mu_intercept", "mu_x")
    ),
    "mu_x", formula = ~ x, prior_samples = TRUE
  )
  level <- levels[["1SD"]]
  class(level) <- c(class(level), "marginal_posterior")
  expect_s3_class(.bt_meta_get(level, "components"), "BayesTools_posterior_components")
  expect_pooled(level)
})

.half_normal_levels_for_test <- function(n, intercept_prior = NULL, seed = 4){

  set.seed(seed)
  half_normal <- prior("normal", list(0, 1), list(0, Inf))
  if(is.null(intercept_prior)){
    intercept_prior <- half_normal
  }
  formula_result <- JAGS_formula(
    ~ f, "mu", data = data.frame(f = factor(c("A", "B", "C"), levels = c("A", "B", "C"))),
    prior_list = list(
      intercept = intercept_prior,
      f         = prior_factor("normal", list(0, 1), list(0, Inf), contrast = "treatment")
    )
  )
  list(formula_result = formula_result, half_normal = half_normal)
}

test_that("linear level hypotheses use the exact support of the combination", {

  n <- 20000
  fixture <- .half_normal_levels_for_test(n)
  # prior-only draws: intercept and both treatment coefficients half-normal
  posterior <- cbind(
    mu_intercept = abs(stats::rnorm(n)),
    "mu_f[1]"    = abs(stats::rnorm(n)),
    "mu_f[2]"    = abs(stats::rnorm(n))
  )
  mixed <- as_mixed_posteriors(
    .single_fit_for_test(posterior, fixture$formula_result$prior_list),
    parameters = c("mu_intercept", "mu_f")
  )
  levels <- marginal_posterior(mixed, "mu_f", formula = ~ f, prior_samples = TRUE)
  half_normal <- function(y) ifelse(y >= 0, 2 * stats::dnorm(y), 0)

  # B - A is the half-normal coefficient: support [0, Inf) although each level
  # alone is supported on [0, Inf); the null 0 lies on the bound
  difference <- hypothesis_BF(levels, hypothesis = "mu_f[B] - mu_f[A] = 0", columns = "all", seed = 1)
  expect_identical(difference$method, "Savage-Dickey")
  expect_equal(as.numeric(difference$prior), 2 * stats::dnorm(0), tolerance = 1e-6)
  moments <- .reflected_kde_moments_for_test(
    0, h = stats::bw.nrd0(posterior[, "mu_f[1]"]), density = half_normal,
    lower = 0, n = n, n_source = Inf
  )
  # the reflected ordinate (the unreflected KDE gave about half of it)
  expect_lte(abs(difference$posterior - moments[["mean"]]), 4 * sqrt(moments[["variance"]]))

  # 2 A + B = 3 intercept + coefficient: [0, Inf); prior by quadrature
  combination <- hypothesis_BF(levels, hypothesis = "2*mu_f[A] + mu_f[B] = 1", columns = "all", seed = 1)
  expect_identical(combination$method, "Savage-Dickey")
  convolution <- stats::integrate(
    function(t) half_normal(t / 3) / 3 * half_normal(1 - t), 0, 1, rel.tol = 1e-10
  )$value
  expect_equal(as.numeric(combination$prior), convolution, tolerance = 1e-4)

  # an unbounded linear expression also uses the exact joint-prior density:
  # B - C = f1 - f2 is the difference of two half-normals, with density
  # int 4 phi(a)^2 da over a > 0 = 1 / sqrt(pi) at 0; the posterior ordinate
  # is the (unreflected) kernel sum of the expression draws
  draws <- data.frame(check.names = FALSE,
                      "mu_f[B]" = as.numeric(levels[["B"]]), "mu_f[C]" = as.numeric(levels[["C"]]),
                      "mu_f[A]" = as.numeric(levels[["A"]]))
  unbounded <- hypothesis_BF(levels, hypothesis = "mu_f[B] - mu_f[C] = 0", columns = "all", seed = 1)
  expect_identical(unbounded$method, "Savage-Dickey")
  expect_equal(as.numeric(unbounded$prior), 1 / sqrt(pi), tolerance = 1e-8)
  expect_equal(
    as.numeric(unbounded$posterior),
    as.numeric(BayesTools:::.hypothesis_sample_density_height(
      draws[["mu_f[B]"]] - draws[["mu_f[C]"]], 0, "posterior"
    )),
    tolerance = 1e-12
  )

  # nonlinear expressions of a deterministic prior have no exact prior
  # ordinate: a kernel estimate from sampled prior draws is refused
  expect_error(
    hypothesis_BF(levels, hypothesis = "exp(mu_f[B]) - exp(mu_f[C]) = 0", seed = 1),
    class = "BayesTools_inexact_ordinate"
  )
  # abs(B - A - 2) equals the linear 2 - (B - A) at every probe point, but its
  # value 0.5 is also reached at B - A = 2.5 (prior density 2 phi(1.5) +
  # 2 phi(2.5), not 2 phi(1.5) alone): it is not linear
  expect_error(
    hypothesis_BF(levels, hypothesis = "abs(mu_f[B] - mu_f[A] - 2) = 0.5", seed = 1),
    class = "BayesTools_inexact_ordinate"
  )
  expect_null(BayesTools:::.hypothesis_linear_coefficients(
    str2lang("abs(`mu_f[B]` - `mu_f[A]` - 2)"), c("mu_f[B]", "mu_f[A]"), draws
  ))
  # functions of constants keep a linear form
  expect_equal(
    BayesTools:::.hypothesis_linear_coefficients(
      str2lang("exp(0) * (`mu_f[B]` - `mu_f[A]`) / abs(-2) + 1"), c("mu_f[B]", "mu_f[A]"), draws
    ),
    list(constant = 1, coefficients = c("mu_f[B]" = .5, "mu_f[A]" = -.5))
  )
})

test_that("linear level hypotheses mix per-component ordinates of mixture terms", {

  n <- 8000
  bounded <- prior("normal", list(.5, 1), list(0, Inf))
  fixture <- .half_normal_levels_for_test(
    n, intercept_prior = prior_mixture(list(prior("normal", list(0, 1)), bounded),
                                       is_null = c(FALSE, FALSE)),
    seed = 6
  )
  indicator <- sample(1:2, n, TRUE)
  posterior <- cbind(
    mu_intercept = ifelse(indicator == 1L, stats::rnorm(n), rng(bounded, n)),
    "mu_f[1]"    = abs(stats::rnorm(n)),
    "mu_f[2]"    = abs(stats::rnorm(n)),
    mu_intercept_indicator = indicator
  )
  mixed <- as_mixed_posteriors(
    .single_fit_for_test(posterior, fixture$formula_result$prior_list),
    parameters = c("mu_intercept", "mu_f")
  )
  levels <- marginal_posterior(mixed, "mu_f", formula = ~ f, prior_samples = TRUE)

  # A + B = 2 intercept + coefficient: the real line with the unbounded
  # intercept component, [0, Inf) with the bounded one
  out <- hypothesis_BF(levels, hypothesis = "mu_f[A] + mu_f[B] = 0.2", columns = "all", seed = 1)
  expect_identical(out$method, "Savage-Dickey")
  values <- 2 * posterior[, "mu_intercept"] + posterior[, "mu_f[1]"]
  heights <- c(
    BayesTools:::.Savage_Dickey_BF.kd(values[indicator == 1L], .2, support = c(-Inf, Inf)),
    BayesTools:::.Savage_Dickey_BF.kd(values[indicator == 2L], .2, support = c(0, Inf))
  )
  expect_equal(as.numeric(out$posterior), sum(tabulate(indicator, 2) / n * heights), tolerance = 1e-12)
})

test_that("linear level hypotheses treat a structurally fixed level as a constant", {

  # Treatment levels on the coefficient scale: the reference level A is
  # structurally 0 (prior and posterior point mass), so B - A is B, with the
  # prior density phi(0) at 0 and the same Bayes factor as B = 0.
  set.seed(7)
  n <- 4000
  formula_result <- JAGS_formula(
    ~ f, "mu", data = data.frame(f = factor(c("A", "B", "C"), levels = c("A", "B", "C"))),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f         = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  posterior <- cbind(mu_intercept = stats::rnorm(n), "mu_f[1]" = stats::rnorm(n, .4, .3),
                     "mu_f[2]" = stats::rnorm(n, -.2, .3))
  mixed <- as_mixed_posteriors(
    .single_fit_for_test(posterior, formula_result$prior_list),
    parameters = c("mu_intercept", "mu_f")
  )
  levels <- marginal_posterior(mixed, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  contrast <- hypothesis_BF(levels, hypothesis = "mu_f[B] - mu_f[A] = 0", columns = "all", seed = 1)
  direct <- hypothesis_BF(levels, hypothesis = "mu_f[B] = 0", columns = "all", seed = 1)
  expect_identical(contrast$method, "Savage-Dickey")
  expect_equal(as.numeric(contrast$prior), stats::dnorm(0), tolerance = 1e-12)
  expect_equal(attr(contrast, "raw_BF"), attr(direct, "raw_BF"), tolerance = 1e-10)
})

.treatment_factor_prior_for_test <- function(sd){

  treatment <- prior_factor_levels(
    prior_factor("normal", list(0, sd), contrast = "treatment"),
    c("A", "B", "C")
  )
  # JAGS_fit() binds the factor metadata of fitted prior lists to their names
  BayesTools:::.complete_factor_metadata_prior_list(list(mu_f = treatment))$mu_f
}

test_that("marginal inference gives levels fixed at the null an NA Bayes factor with its reason", {

  fixed <- paste0(
    "The posterior is fixed at the null hypothesis value. The ",
    "Savage-Dickey Bayes factor is undefined."
  )
  point_mass <- paste0(
    "The posterior contains a declared point mass at the exact null ",
    "hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
  )
  set.seed(3)
  n <- 2000
  factor_draws <- function(){
    cbind("mu_f[1]" = stats::rnorm(n, .3, .4), "mu_f[2]" = stats::rnorm(n, .6, .4))
  }

  # treatment-coded coefficients: the reference level A is fixed at 0
  models <- lapply(c(1, .5), function(sd) list(
    fit = .mock_mixing_fit_for_marginal(factor_draws(), list(mu_f = .treatment_factor_prior_for_test(sd))),
    marglik = bridgesampling_object(0), prior_weights = 1
  ))
  inference <- .collect_warnings_for_test(marginal_inference(
    models, marginal_parameters = "mu_f", parameters = "mu_f",
    is_null_list = list(mu_f = c(FALSE, FALSE)), formula = NULL, n_samples = n, seed = 1
  ))
  expect_identical(inference$warnings, paste0("mu_f[A]: ", fixed))
  bf <- inference$value$inference$mu_f
  expect_identical(names(bf), c("A", "B", "C"))
  expect_true(is.na(bf[["A"]]))
  expect_identical(attr(bf[["A"]], "warnings"), fixed)
  levels <- lapply(inference$value$conditional$mu_f, function(level){
    class(level) <- c(class(level), "marginal_posterior")
    level
  })
  # the other levels match their direct Savage-Dickey Bayes factors
  for(level in c("B", "C")){
    expect_true(is.finite(bf[[level]]))
    expect_identical(bf[[level]], Savage_Dickey_BF(levels[[level]]))
  }
  # a scalar call keeps its error
  expect_error(Savage_Dickey_BF(levels[["A"]]), point_mass, fixed = TRUE)

  table <- marginal_estimates_table(
    inference$value$conditional, inference$value$inference, "mu_f"
  )
  expect_equal(as.numeric(table$inclusion_BF), c(NA, as.numeric(bf[["B"]]), as.numeric(bf[["C"]])))
  expect_identical(attr(table, "warnings"), paste0("mu_f[A]: ", fixed))

  # a partial point mass at the null: NA with the point-mass reason
  partial <- inference$value$conditional$mu_f
  partial[["B"]] <- .bt_meta_set(partial[["B"]], "atoms", posterior_atom_attribute(data.frame(x = 0, mass = .5)))
  partial_BF <- Savage_Dickey_BF(partial, silent = TRUE)
  expect_true(is.na(partial_BF[["B"]]))
  expect_identical(attr(partial_BF[["B"]], "warnings"), point_mass)
  expect_identical(partial_BF[["C"]], bf[["C"]])

  # as_marginal_inference on one fit
  fit <- coda::mcmc(factor_draws())
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(mu_f = .treatment_factor_prior_for_test(1))
  fit <- attach_test_parameter_map(fit)
  single <- .collect_warnings_for_test(as_marginal_inference(
    fit, marginal_parameters = "mu_f", parameters = "mu_f",
    conditional_list = list(mu_f = NULL), conditional_rule = "AND",
    formula = NULL, n_samples = n
  ))
  expect_identical(single$warnings, paste0("mu_f[A]: ", fixed))
  single_BF <- single$value$inference$mu_f
  expect_true(is.na(single_BF[["A"]]))
  expect_true(all(is.finite(unlist(single_BF[c("B", "C")]))))
})

test_that("marginal inference gives levels with a zero or infinite prior ordinate an NA Bayes factor", {

  # An ordered factor with the default flat Dirichlet allocation on the
  # coefficient scale (no intercept), averaged over two models whose totals
  # are N(0, 1) and N(0, .5): the first level is fixed at 0, the second is
  # total * S with S ~ Beta(1, 1), whose prior density at 0 is infinite in
  # both models (E[1 / S] diverges), and the third is the total, with the
  # regular prior density (dnorm(0) + dnorm(0, 0, .5)) / 2 at 0. The second
  # level gets an NA Bayes factor with the prior-ordinate reason, the first
  # the fixed-level reason, and the third its ordinary Savage-Dickey Bayes
  # factor.
  infinite <- paste0(
    "The prior density at the null hypothesis value is infinite. The ",
    "Savage-Dickey Bayes factor is undefined."
  )
  fixed <- paste0(
    "The posterior is fixed at the null hypothesis value. The ",
    "Savage-Dickey Bayes factor is undefined."
  )
  ordered_prior <- function(sd){
    total <- prior_factor_levels(
      prior_ordered(prior("normal", list(0, sd))),
      c("low", "mid", "high")
    )
    BayesTools:::.complete_factor_metadata_prior_list(list(mu_f = total))$mu_f
  }
  set.seed(5)
  n <- 2000
  models <- lapply(c(1, .5), function(sd){
    prior <- ordered_prior(sd)
    total <- stats::rnorm(n,.15,.1)
    gamma <- matrix(stats::rexp(2*n),n)
    spec <- .bt_ordered_spec("mu_f",prior)
    colnames(gamma) <- spec$allocations[[1L]]$gamma_coordinates
    shares <- gamma/rowSums(gamma)
    draws <- cbind("mu_f[1]"=total*shares[,1],"mu_f[2]"=total*shares[,2],
      mu_f_ordered_total=total,gamma)
    list(fit=.mock_mixing_fit_for_marginal(draws,list(mu_f=prior)),
      marglik=bridgesampling_object(0),prior_weights=1)
  })
  inference <- .collect_warnings_for_test(marginal_inference(
    models, marginal_parameters = "mu_f", parameters = "mu_f",
    is_null_list = list(mu_f = c(FALSE, FALSE)), formula = NULL, n_samples = n, seed = 1
  ))
  bf <- inference$value$inference$mu_f
  expect_identical(names(bf), c("low", "mid", "high"))
  expect_true(is.na(bf[["low"]]))
  expect_identical(attr(bf[["low"]], "warnings"), fixed)
  expect_true(is.na(bf[["mid"]]))
  expect_identical(attr(bf[["mid"]], "warnings"), infinite)
  expect_identical(attr(bf[["mid"]], "posterior_density_source"), "infinite_prior_ordinate")
  # the levels of the transformed ordered contrast are labelled by their cells
  expect_identical(inference$warnings, paste0("mu_f[", c("low", "mid"), "]: ", c(fixed, infinite)))

  high <- inference$value$conditional$mu_f[["high"]]
  class(high) <- c(class(high), "marginal_posterior")
  expect_equal(exp(prior_density_ordinate(.bt_meta_get(high, "prior_density"), 0)$log_density),
               (stats::dnorm(0) + stats::dnorm(0, 0, .5)) / 2, tolerance = 1e-14)
  expect_true(is.finite(bf[["high"]]))
  expect_identical(bf[["high"]], Savage_Dickey_BF(high))
  # a scalar call on the level with the infinite prior ordinate stops
  mid <- inference$value$conditional$mu_f[["mid"]]
  class(mid) <- c(class(mid), "marginal_posterior")
  condition <- tryCatch(Savage_Dickey_BF(mid, silent = TRUE), error = function(e) e)
  expect_s3_class(condition, "BayesTools_infinite_ordinate")

  table <- marginal_estimates_table(
    inference$value$conditional, inference$value$inference, "mu_f"
  )
  expect_true(all(is.na(table$inclusion_BF[1:2])))
  expect_equal(as.numeric(table$inclusion_BF[3]), as.numeric(bf[["high"]]))
  expect_identical(attr(table, "warnings"), paste0("mu_f[", c("low", "mid"), "]: ", c(fixed, infinite)))
})

test_that("use_formula = FALSE prior densities ignore the coefficient's own multiply_by", {

  # JAGS monitors the raw coefficient; 'multiply_by' only scales the linear
  # predictor, so the coefficient prior is its declared prior (exact normal
  # ordinates; the product sigma * beta would be singular at zero).
  set.seed(1)
  n <- 200
  sigma_prior <- prior("normal", list(0, 1), list(0, Inf))
  beta_prior_1 <- prior("normal", list(0, 1))
  beta_prior_2 <- prior("normal", list(0, 2))
  attr(beta_prior_1, "multiply_by") <- "sigma"
  attr(beta_prior_2, "multiply_by") <- "sigma"
  models <- lapply(list(beta_prior_1, beta_prior_2), function(beta_prior){
    list(
      fit = .mock_mixing_fit_for_marginal(
        cbind(mu_x = rnorm(n), sigma = abs(rnorm(n, 2, .1))),
        list(mu_x = beta_prior, sigma = sigma_prior)
      ),
      marglik = bridgesampling_object(0),
      prior_weights = 1
    )
  })
  mixed <- mix_posteriors(
    models,
    parameters   = c("mu_x", "sigma"),
    is_null_list = list(mu_x = c(FALSE, FALSE), sigma = c(FALSE, FALSE)),
    seed         = 1,
    n_samples    = n
  )
  marginal <- marginal_posterior(mixed, "mu_x", use_formula = FALSE, prior_samples = TRUE)

  for(value in c(0, 1)){
    expect_equal(
      .prior_height_for_test(marginal, value),
      .5 * stats::dnorm(value) + .5 * stats::dnorm(value, 0, 2),
      tolerance = 1e-8
    )
  }
  expect_true(is.finite(Savage_Dickey_BF(marginal, silent = TRUE)))

  # Conditional single-model context: included component only.
  spike_slab <- prior_spike_and_slab(prior("normal", list(0, 1)))
  attr(spike_slab, "multiply_by") <- "sigma"
  indicator <- rep(c(0, 1), length.out = n)
  fit <- coda::mcmc(cbind(
    mu_x = indicator * rnorm(n),
    mu_x_indicator = indicator,
    sigma = abs(rnorm(n, 2, .1))
  ))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(mu_x = spike_slab, sigma = sigma_prior)
  fit <- attach_test_parameter_map(fit)
  conditional <- as_mixed_posteriors(
    fit,
    parameters = c("mu_x", "sigma"),
    conditional = "mu_x",
    n_prior_samples = 1000
  )
  conditional_marginal <- marginal_posterior(
    conditional, "mu_x", use_formula = FALSE, prior_samples = TRUE
  )
  expect_equal(.prior_height_for_test(conditional_marginal, 0), stats::dnorm(0), tolerance = 1e-8)
  expect_equal(.prior_height_for_test(conditional_marginal, 1), stats::dnorm(1), tolerance = 1e-8)
})

test_that("transform_scaled slope stays raw and intercept uses compiled products", {

  set.seed(1)
  data <- data.frame(x = rnorm(50, 3, 2))
  x_prior <- prior_spike_and_slab(prior("normal", list(0, 1)))
  attr(x_prior, "multiply_by") <- "sigma"
  formula_result <- JAGS_formula(
    ~ x, parameter = "mu", data = data,
    prior_list = list(intercept = prior("normal", list(0, 1)), x = x_prior),
    formula_scale = list(x = TRUE)
  )
  scale <- formula_result$formula_scale[["mu_x"]]
  n <- 400
  indicator <- rep(c(0, 1), length.out = n)
  posterior <- cbind(
    mu_intercept = rnorm(n), mu_x = indicator * rnorm(n),
    mu_x_indicator = indicator, sigma = rlnorm(n)
  )
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = n
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- c(formula_result$prior_list, list(sigma = prior("lognormal", list(0, 1))))
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(
    fit, c("mu_intercept", "mu_x", "sigma"), transform_scaled = TRUE, n_prior_samples = 2000
  )
  prior_densities <- .bt_meta_get(samples, "prior_densities")
  ratio <- scale$mean / scale$sd
  # original-scale raw slope b / s: 1/2 point at 0 + 1/2 N(0, 1 / s)
  expect_equal(BayesTools:::.prior_linear_density_point_mass(prior_densities$mu_x, 0), .5, tolerance = 1e-12)
  for(value in c(.2, 1)){
    .expect_prior_height_for_test(
      .bt_meta_update(
        structure(0),
        prior_density = prior_densities$mu_x
      ),
      value, .5 * stats::dnorm(value, 0, 1 / scale$sd)
    )
  }
  # The original intercept is b0 - (m / s) b sigma. Integrate independently
  # over log(sigma) ~ N(0, 1), alongside the spike component's N(0, 1).
  for(value in c(.5, -1)){
    .expect_prior_height_for_test(
      .bt_meta_update(
        structure(0),
        prior_density = prior_densities$mu_intercept
      ),
      value, .5 * stats::dnorm(value) + .5 * stats::integrate(function(log_sigma){
        stats::dnorm(value, 0, sqrt(1 + ratio^2 * exp(2 * log_sigma))) * stats::dnorm(log_sigma)
      }, -Inf, Inf, rel.tol = 1e-11, subdivisions = 500L)$value
    )
  }
  # linear-predictor targets keep the formula prior's multiply_by
  expect_identical(
    attr(.bt_meta_get(samples, "prior_context")$prior_list$mu_x, "multiply_by"),
    "sigma"
  )
})

test_that("marginal_posterior uses log(intercept) for log-intercept formulas", {

  log_formula <- ~ x
  attr(log_formula, "log(intercept)") <- TRUE
  formula_result <- JAGS_formula(
    formula = log_formula,
    parameter = "ls",
    data = data.frame(x = c(-1, 0, 1)),
    prior_list = list(
      intercept = prior("lognormal", list(0, .5)),
      x         = prior("normal", list(0, .5))
    )
  )
  set.seed(2)
  n <- 100
  posterior <- cbind(ls_intercept = stats::rlnorm(n, -1, .2), ls_x = stats::rnorm(n, -.2, .1))
  make_fit <- function(design){
    fit <- coda::mcmc(posterior)
    class(fit) <- c("BayesTools_fit", class(fit))
    attr(fit, "prior_list") <- formula_result$prior_list
    fit <- attach_test_parameter_map(fit)
    # Only the persisted log(intercept) flag of the design is under test.
    attr(fit, "formula_design") <- design
    fit
  }
  expected <- lapply(c(-1, 0, 1), function(x){
    log(posterior[, "ls_intercept"]) + x * posterior[, "ls_x"]
  })

  # persisted fitted-design metadata marks the formula as log(intercept)
  samples <- as_mixed_posteriors(
    make_fit(list(ls = list(log_intercept = TRUE))),
    parameters = c("ls_intercept", "ls_x"),
    n_prior_samples = 1000
  )
  marginal <- marginal_posterior(samples, "ls_x", formula = ~ x, prior_samples = TRUE)
  expect_equal(unname(lapply(marginal, as.numeric)), expected, tolerance = 1e-12)
  # log(intercept) ~ N(0, .5) and x ~ N(0, .5): the level at x is N(0, .5 * sqrt(1 + x^2))
  .expect_prior_height_for_test(marginal[["1SD"]], -1, stats::dnorm(-1, 0, .5 * sqrt(2)))
  .expect_prior_height_for_test(marginal[["0SD"]], -1, stats::dnorm(-1, 0, .5))

  # an explicit formula attribute without persisted metadata
  attributed <- marginal_posterior(
    as_mixed_posteriors(make_fit(NULL), parameters = c("ls_intercept", "ls_x")),
    "ls_x",
    formula = log_formula
  )
  expect_equal(unname(lapply(attributed, as.numeric)), expected, tolerance = 1e-12)

  linear_formula <- ~ x
  attr(linear_formula, "log(intercept)") <- FALSE
  expect_error(
    marginal_posterior(samples, "ls_x", formula = linear_formula),
    "does not match the fitted formula",
    fixed = TRUE
  )

  # model-averaged draws carry the fitted-design flag from every model
  models <- lapply(1:2, function(i){
    fit <- .mock_mixing_fit_for_marginal(posterior, formula_result$prior_list)
    attr(fit, "formula_design") <- list(ls = list(log_intercept = TRUE))
    list(fit = fit, marglik = bridgesampling_object(0), prior_weights = 1)
  })
  mixed <- mix_posteriors(
    models,
    parameters   = c("ls_intercept", "ls_x"),
    is_null_list = list(ls_intercept = c(FALSE, FALSE), ls_x = c(FALSE, FALSE)),
    seed         = 1,
    n_samples    = 50
  )
  mixed_marginal <- marginal_posterior(mixed, "ls_x", formula = ~ x)
  rows <- .bt_meta_get(mixed$ls_x, "draw_index")
  expect_equal(
    as.numeric(mixed_marginal[["1SD"]]),
    log(posterior[rows, "ls_intercept"]) + posterior[rows, "ls_x"],
    tolerance = 1e-12
  )
})

test_that("marginal_posterior without prior samples declares scaled log-intercept metadata", {

  log_formula <- ~ x
  attr(log_formula, "log(intercept)") <- TRUE
  set.seed(5)
  formula_result <- JAGS_formula(
    log_formula, "ls", data = data.frame(x = rnorm(50, 3, 1.5)),
    prior_list = list(intercept = prior("lognormal", list(0, .5)), x = prior("normal", list(0, .5))),
    formula_scale = list(x = TRUE)
  )
  n <- 100
  posterior <- cbind(ls_intercept = stats::rlnorm(n, -1, .2), ls_x = stats::rnorm(n, -.2, .1))
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = n
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(ls = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(ls = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(
    fit, parameters = c("ls_intercept", "ls_x"),
    transform_scaled = TRUE, n_prior_samples = 500
  )
  original <- BayesTools:::.bt_transform_scale_posterior(posterior, list(ls = formula_result$formula_scale))

  levels <- marginal_posterior(samples, "ls_x", formula = ~ x, prior_samples = FALSE)
  expect_equal(
    unname(lapply(levels, as.numeric)),
    lapply(c(-1, 0, 1), function(x) log(original[, "ls_intercept"]) + x * original[, "ls_x"]),
    tolerance = 1e-12
  )
  # the scaled log-intercept levels are linear in (log(b0), b1) through the log
  # source: log(b0) + k b1 is normal (lognormal intercept, normal slope), so
  # the exact support is the real line and the declared atoms are none
  for(level in names(levels)){
    expect_true(posterior_atoms_free(levels[[level]]))
    support <- .bt_meta_get(levels[[level]], "support")
    expect_identical(support$bounds, c(-Inf, Inf))
    expect_true(support$exact)
  }

  intercept <- marginal_posterior(samples, "ls_intercept", formula = ~ x, prior_samples = FALSE)
  expect_equal(as.numeric(intercept[["intercept"]]), log(original[, "ls_intercept"]), tolerance = 1e-12)
})

test_that("intercept-only formulas give the intercept level", {

  set.seed(3)
  n <- 100L
  truncated <- prior("normal", list(.5, 1), list(0, Inf))
  data <- data.frame(y = 1:3)
  intercept_only <- function(intercept){
    JAGS_formula(~ 1, "mu", data = data, prior_list = list(intercept = intercept))$prior_list
  }
  # N(0, 1) and N(0.5, 1)T(0, Inf) mixed with equal weights; the truncated
  # component's ordinate at its bound is the one-sided limit.
  height_0 <- .5 * stats::dnorm(0) + .5 * stats::dnorm(0, .5) / stats::pnorm(.5)

  models <- list(
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu_intercept = stats::rnorm(n)),
                                             intercept_only(prior("normal", list(0, 1)))),
         marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(cbind(mu_intercept = rng(truncated, n)),
                                             intercept_only(truncated)),
         marglik = bridgesampling_object(0), prior_weights = 1)
  )
  mixed <- mix_posteriors(models, parameters = "mu_intercept",
                          is_null_list = list(mu_intercept = c(FALSE, FALSE)),
                          seed = 1, n_samples = n)
  levels <- marginal_posterior(mixed, "mu_intercept", formula = ~ 1, prior_samples = TRUE)
  expect_identical(names(levels), "intercept")
  expect_equal(as.numeric(levels[["intercept"]]), as.numeric(mixed$mu_intercept))
  expect_equal(.prior_height_for_test(levels[["intercept"]], 0), height_0, tolerance = 1e-10)

  indicator <- sample(1:2, n, TRUE)
  draws <- ifelse(indicator == 1L, stats::rnorm(n), rng(truncated, n))
  fit <- coda::mcmc(cbind(mu_intercept = draws, mu_intercept_indicator = indicator))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- intercept_only(prior_mixture(
    list(prior("normal", list(0, 1)), truncated), is_null = c(FALSE, FALSE)
  ))
  fit <- attach_test_parameter_map(fit)
  single <- marginal_posterior(as_mixed_posteriors(fit, parameters = "mu_intercept"),
                               "mu_intercept", formula = ~ 1, prior_samples = TRUE)
  expect_identical(names(single), "intercept")
  expect_equal(as.numeric(single[["intercept"]]), draws)
  expect_equal(.prior_height_for_test(single[["intercept"]], 0), height_0, tolerance = 1e-10)
})

.ordered_prior_for_test <- function(total, allocation = c(.4, .6)){
  data <- data.frame(f = ordered(c("low", "mid", "high"), levels = c("low", "mid", "high")))
  JAGS_formula(
    ~ f, "mu", data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f         = prior_ordered(total, allocation = allocation)
    )
  )$prior_list$mu_f
}

test_that("transform_scaled factor atoms are rescaled for level-labelled columns", {

  set.seed(6)
  data <- data.frame(
    x = rnorm(40, 5, 2),
    f = factor(rep(c("a", "b", "c"), length.out = 40), levels = c("a", "b", "c"))
  )
  formula_result <- JAGS_formula(
    ~ x * f, "mu", data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 5)),
      x         = prior("normal", list(0, 1)),
      f         = prior_spike_and_slab(prior_factor("normal", list(0, 1), contrast = "treatment")),
      "x:f"     = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    formula_scale = list(x = TRUE)
  )
  n <- 200
  indicator <- rep(c(0, 1), length.out = n)
  posterior <- cbind(
    mu_intercept = rnorm(n), mu_x = rnorm(n),
    "mu_f[1]" = indicator * rnorm(n), "mu_f[2]" = indicator * rnorm(n),
    mu_f_indicator = indicator,
    "mu_x__xXx__f[1]" = rnorm(n), "mu_x__xXx__f[2]" = rnorm(n)
  )
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = n
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(
    fit, c("mu_intercept", "mu_x", "mu_f", "mu_x__xXx__f"), transform_scaled = TRUE
  )
  # original-scale f = f* - (m / s) (x:f)* is continuous even in spike draws
  expect_equal(colnames(samples$mu_f), c("mu_f[b]", "mu_f[c]"))
  expect_true(all(unclass(samples$mu_f) != 0))
  atoms <- .bt_meta_get(samples$mu_f, "atoms")
  expect_true(atoms$declared)
  expect_length(atoms$mass, 0L)
  expect_equal(colnames(atoms$locations), c("mu_f[b]", "mu_f[c]"))

  marginal <- marginal_posterior(samples, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  for(level in c("b", "c")){
    level_posterior <- marginal[[level]]
    class(level_posterior) <- c(class(level_posterior), "marginal_posterior")
    expect_true(is.finite(Savage_Dickey_BF(level_posterior, silent = TRUE)))
  }
})

test_that("marginal posteriors of 0.3.0 objects ask for recomputation", {

  fixture <- .marginal_semantic_fixture_for_test()
  samples <- fixture$samples
  # BayesTools 0.3.0 mixed factor posteriors did not record the ordered flag
  attr(samples$mu_x_fac2t, "ordered") <- NULL

  expect_error(
    marginal_posterior(samples, "mu_x_fac2t", formula = ~ x_cont1 + x_fac2t + x_cont1 * x_fac3md),
    "lack the factor metadata recorded by the current BayesTools version (missing: 'ordered')",
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
  # an interaction-only factor predictor without its fitted levels and contrast
  condition <- tryCatch(
    .marginal_posterior_interaction_predictor_info(
      predictor   = "g",
      parameter   = "mu_x",
      priors_info = list(mu_x__xXx__g = list(
        term_components = c("x", "g"),
        factor_terms    = "g"
      ))
    ),
    error = identity
  )
  expect_identical(
    class(condition),
    c("BayesTools_refit_required", "error", "condition")
  )
  expect_identical(
    conditionMessage(condition),
    paste0(
      "The mixed posterior samples lack the fitted levels or contrast of the ",
      "factor predictor 'g'. Recreate them with mix_posteriors() or ",
      "as_mixed_posteriors() from models fitted with the current version."
    )
  )

  # 0.3.0 marginal posteriors carry no atom declaration
  legacy <- marginal_posterior(fixture$samples, "mu_x_cont1", use_formula = FALSE, prior_samples = TRUE)
  legacy <- .bt_meta_set(legacy, "atoms", NULL)
  expect_error(
    Savage_Dickey_BF(legacy, silent = TRUE),
    "Marginal posteriors created by BayesTools 0.3.0 do not record it",
    fixed = TRUE
  )
})

test_that("point-mass metadata merge atoms by exact location", {

  location <- 0.1 + 0.2
  expect_false(location == 0.3)
  atoms <- posterior_atom_attribute(data.frame(x = c(location, location, 0.3), mass = c(.2, .2, .1)))
  expect_identical(atoms$locations[, 1], c(0.3, location))
  expect_equal(atoms$mass, c(.1, .4), tolerance = 1e-15)

  # Savage-Dickey removes the draws of the merged atom from the continuous part
  set.seed(9)
  posterior <- c(rep(location, 400), rnorm(600, 1, .3))
  class(posterior) <- c("marginal_posterior.simple", "marginal_posterior", class(posterior))
  posterior <- .bt_meta_set(posterior, "prior_density", prior("normal", list(0, 1)))
  posterior <- .bt_meta_set(posterior, "atoms", posterior_atom_attribute(
    data.frame(x = c(location, location), mass = c(.2, .2))
  ))
  continuous <- BayesTools:::.Savage_Dickey_BF.continuous_posterior(
    posterior, BayesTools:::.posterior_atoms_get(posterior)
  )
  expect_length(continuous$samples, 600L)
  expect_equal(continuous$continuous_mass, .6, tolerance = 1e-12)

  density <- BayesTools:::.posterior_density_from_attribute(.posterior_density_for_test(
    x = c(0, location, 0.3, location, 1),
    y = c(1, 2, 3, 4, 5)
  ))
  expect_identical(density$x, c(0, 0.3, location, 1))
  expect_equal(density$y, c(1, 3, 3, 5))
})

test_that("factor terms omitted by a mixed model are zero on every coefficient column", {

  data <- data.frame(t = factor(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi")))
  formula_result <- JAGS_formula(
    ~ 1 + t, "mu", data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      t         = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  set.seed(7)
  n <- 200
  alternative <- cbind(mu_intercept = rnorm(n), "mu_t[1]" = rnorm(n, .4), "mu_t[2]" = rnorm(n, 1))
  null <- cbind(mu_intercept = rnorm(n))
  models <- list(
    list(fit = .mock_mixing_fit_for_marginal(alternative, formula_result$prior_list),
         marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(null, formula_result$prior_list["mu_intercept"]),
         marglik = bridgesampling_object(0), prior_weights = 1)
  )
  mixed <- mix_posteriors(
    models, parameters = c("mu_intercept", "mu_t"),
    is_null_list = list(mu_intercept = c(FALSE, FALSE), mu_t = c(FALSE, TRUE)),
    seed = 1, n_samples = n
  )

  levels <- marginal_posterior(mixed, "mu_t", formula = ~ 1 + t)
  expect_equal(as.numeric(levels[["mid"]]), as.numeric(mixed$mu_intercept) + unclass(mixed$mu_t)[, 1],
               tolerance = 1e-12)

  # the omitted term contributes a point at zero with prior model probability 1/2
  coefficients <- marginal_posterior(mixed, "mu_t", use_formula = FALSE, prior_samples = TRUE)
  for(level in c("mid", "hi")){
    prior_density <- .bt_meta_get(coefficients[[level]], "prior_density")
    expect_equal(prior_density_ordinate(prior_density, 0)$point_mass, .5, tolerance = 1e-12)
    .expect_prior_height_for_test(coefficients[[level]], .5, .5 * stats::dnorm(.5))
    expect_equal(.bt_meta_get(coefficients[[level]], "atoms")$mass, .5, tolerance = 1e-12)
  }
  formula_levels <- marginal_posterior(mixed, "mu_t", formula = ~ 1 + t, prior_samples = TRUE)
  # level mid = intercept + coefficient: N(0, 1) + {0 or N(0, 1)}
  .expect_prior_height_for_test(
    formula_levels[["mid"]], .5,
    .5 * stats::dnorm(.5) + .5 * stats::dnorm(.5, 0, sqrt(2))
  )
})

test_that("terms with unknown support leave level support unknown", {

  data <- data.frame(f = ordered(c("low", "mid", "high"), levels = c("low", "mid", "high")))
  formula_result <- JAGS_formula(
    ~ f, "mu", data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1), list(0, Inf)),
      f         = prior_ordered(prior("normal", list(0, 1)), allocation = c(.4, .6))
    )
  )
  total <- seq(-3,2,length.out=20)
  posterior <- cbind(mu_intercept = seq(.1, 2, length.out = 20),
                     "mu_f[1]" = .4*total,"mu_f[2]" = .6*total,
                     mu_f_ordered_total=total)
  fit <- coda::mcmc(posterior)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)
  samples <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_f"))

  levels <- marginal_posterior(samples, "mu_f", formula = ~ f)
  # the reference level is the truncated intercept alone; the other levels add
  # an ordered term whose support is not derived, so their support is unknown
  expect_equal(BayesTools:::.posterior_support_bounds(levels[["low"]]), c(0, Inf))
  expect_null(.bt_meta_get(levels[["mid"]], "support"))
  expect_null(.bt_meta_get(levels[["high"]], "support"))

  known <- BayesTools:::.posterior_support_new(c(0, Inf))
  expect_null(BayesTools:::.posterior_support_sum(list(known, NULL)))
  expect_null(BayesTools:::.posterior_support_union(list(known, NULL)))
  expect_equal(BayesTools:::.posterior_support_sum(list(known, known))$bounds, c(0, Inf))
})

test_that("mixed ordered spike-and-slab totals declare their within-model spike", {

  ordered_prior <- .ordered_prior_for_test(prior_spike_and_slab(prior("normal", list(0, 1))))
  set.seed(3)
  n <- 400
  indicators <- list(rbinom(n, 1, .3), rbinom(n, 1, .6))
  models <- lapply(indicators, function(indicator){
    variable <- rnorm(n,.3,.1)
    total <- indicator * variable
    posterior <- cbind(
      "mu_f[1]" = .4 * total,
      "mu_f[2]" = .6 * total,
      "mu_f_ordered_total" = total,
      "mu_f_ordered_total_variable" = variable,
      "mu_f_ordered_total_indicator" = indicator
    )
    list(
      fit = .mock_mixing_fit_for_marginal(posterior, list(mu_f = ordered_prior)),
      marglik = bridgesampling_object(0),
      prior_weights = 1
    )
  })
  mixed <- mix_posteriors(
    models, parameters = "mu_f", is_null_list = list(mu_f = c(FALSE, FALSE)),
    seed = 1, n_samples = n
  )

  # posterior model probabilities are 1/2; each spike carries the model's
  # posterior exclusion probability
  atoms <- .bt_meta_get(mixed$mu_f, "atoms")
  expect_equal(unname(atoms$locations), matrix(0, 2, 2))
  expect_equal(atoms$mass, .5 * vapply(indicators, function(x) mean(x == 0), numeric(1)), tolerance = 1e-12)

  marginal <- marginal_posterior(mixed, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  for(level in c("mid", "high")){
    level_atoms <- .bt_meta_get(marginal[[level]], "atoms")
    expect_equal(unname(level_atoms$locations[, 1]), 0)
    expect_equal(level_atoms$mass, sum(atoms$mass), tolerance = 1e-12)
    # the spike-and-slab total excludes the effect with prior probability 1/2
    ordinate <- prior_density_ordinate(.bt_meta_get(marginal[[level]], "prior_density"), 0)
    expect_equal(ordinate$point_mass, .5, tolerance = 1e-12)
  }

  unmonitored <- models[[1]]
  unmonitored$fit <- .mock_mixing_fit_for_marginal(
    as.matrix(unmonitored$fit$mcmc[[1]])[, c("mu_f[1]", "mu_f[2]")],
    list(mu_f = ordered_prior)
  )
  expect_error(
    mix_posteriors(
      list(unmonitored, models[[2]]), parameters = "mu_f",
      is_null_list = list(mu_f = c(FALSE, FALSE)), seed = 1, n_samples = n
    ),
    "required total-prior indicator",
    fixed = TRUE
  )
})

test_that("ordered mixture totals with a spike(0) component declare their point mass", {

  total <- prior_mixture(
    list(prior("spike", list(0)), prior("normal", list(0, 1))),
    is_null = c(TRUE, FALSE)
  )
  ordered_prior <- .ordered_prior_for_test(total)
  mixture_posterior <- function(indicator){
    total_draws <- ifelse(indicator == 1L, 0, rnorm(length(indicator)))
    cbind(
      "mu_f[1]" = .4 * total_draws,
      "mu_f[2]" = .6 * total_draws,
      "mu_f_ordered_total" = total_draws,
      "mu_f_ordered_total_indicator" = indicator
    )
  }

  # model averaging: model probability 1/2 times the posterior spike fraction
  set.seed(11)
  n <- 400
  indicators <- list(sample(1:2, n, TRUE, c(.3, .7)), sample(1:2, n, TRUE, c(.6, .4)))
  models <- lapply(indicators, function(indicator){
    list(
      fit = .mock_mixing_fit_for_marginal(mixture_posterior(indicator), list(mu_f = ordered_prior)),
      marglik = bridgesampling_object(0),
      prior_weights = 1
    )
  })
  mixed <- mix_posteriors(
    models, parameters = "mu_f", is_null_list = list(mu_f = c(FALSE, FALSE)),
    seed = 1, n_samples = n
  )
  atoms <- .bt_meta_get(mixed$mu_f, "atoms")
  expect_equal(unname(atoms$locations), matrix(0, 2, 2))
  expect_equal(atoms$mass, .5 * vapply(indicators, function(x) mean(x == 1L), numeric(1)), tolerance = 1e-12)

  # prior-only single model: the declared mass is the observed zero fraction
  # and the Savage-Dickey ratio of prior draws is 1 up to KDE error
  # (relative sd ~2%, bias ~1% for 10,000 continuous draws; |log BF| < 0.1)
  set.seed(12)
  n <- 20000
  prior_draws <- mixture_posterior(sample(1:2, n, TRUE))
  fit <- coda::mcmc(prior_draws)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(mu_f = ordered_prior)
  fit <- attach_test_parameter_map(fit)
  samples <- as_mixed_posteriors(fit, parameters = "mu_f", n_prior_samples = 1000)
  expect_equal(
    .bt_meta_get(samples$mu_f, "atoms")$mass,
    mean(prior_draws[, "mu_f[1]"] == 0),
    tolerance = 1e-12
  )

  marginal <- marginal_posterior(samples, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  for(level in c("mid", "high")){
    level_posterior <- marginal[[level]]
    expect_equal(
      .bt_meta_get(level_posterior, "atoms")$mass,
      mean(level_posterior == 0),
      tolerance = 1e-12
    )
    expect_equal(prior_density_ordinate(.bt_meta_get(level_posterior, "prior_density"), 0)$point_mass, .5,
                 tolerance = 1e-12)
    class(level_posterior) <- c(class(level_posterior), "marginal_posterior")
    for(null_hypothesis in c(.05, -.3)){
      BF <- Savage_Dickey_BF(level_posterior, null_hypothesis = null_hypothesis, silent = TRUE)
      expect_lt(abs(log(as.numeric(BF))), .1)
    }
  }
})

test_that("mixed formula levels declare within-model ordered-total spikes", {

  data <- data.frame(f = ordered(c("low", "mid", "high"), levels = c("low", "mid", "high")))
  formula_priors <- function(total){
    JAGS_formula(
      ~ f, "mu", data = data,
      prior_list = list(
        intercept = prior("spike", list(.3)),
        f         = prior_ordered(total, allocation = c(.4, .6))
      )
    )$prior_list
  }
  spike_and_slab_priors <- formula_priors(prior_spike_and_slab(prior("normal", list(0, 1))))
  mixture_priors <- formula_priors(prior_mixture(
    list(prior("spike", list(0)), prior("normal", list(0, 1))),
    is_null = c(TRUE, FALSE)
  ))
  ordered_draws <- function(total, indicator){
    cbind(
      mu_intercept = .3,
      "mu_f[1]" = .4 * total,
      "mu_f[2]" = .6 * total,
      "mu_f_ordered_total" = total,
      "mu_f_ordered_total_indicator" = indicator
    )
  }

  # prior-only draws of three models with the point intercept .3: a
  # spike-and-slab total (spike: indicator 0), a mixture total (spike:
  # component 1), each zero with probability 1/2, and a model without 'f'
  set.seed(21)
  n <- 30000
  ss_indicator  <- stats::rbinom(n, 1, .5)
  mix_indicator <- sample(1:2, n, TRUE)
  models <- list(
    list(fit = .mock_mixing_fit_for_marginal(
      ordered_draws(ss_indicator * stats::rnorm(n), ss_indicator), spike_and_slab_priors),
      marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(
      ordered_draws(ifelse(mix_indicator == 1L, 0, stats::rnorm(n)), mix_indicator), mixture_priors),
      marglik = bridgesampling_object(0), prior_weights = 1),
    list(fit = .mock_mixing_fit_for_marginal(
      cbind(mu_intercept = rep(.3, n)), spike_and_slab_priors["mu_intercept"]),
      marglik = bridgesampling_object(0), prior_weights = 1)
  )
  mixed <- mix_posteriors(
    models, parameters = c("mu_intercept", "mu_f"),
    is_null_list = list(mu_intercept = c(FALSE, FALSE, FALSE), mu_f = c(FALSE, FALSE, TRUE)),
    seed = 1, n_samples = n
  )

  # per-draw total components follow the mixture draws (NA without a spiked
  # total): the spike is component 2 of the spike-and-slab total and
  # component 1 of the mixture total
  models_ind <- .bt_meta_get(mixed$mu_f, "component")
  total_indicator <- .bt_meta_get(mixed$mu_f, "ordered_total_component")
  expect_length(total_indicator, n)
  expect_true(all(is.na(total_indicator[models_ind == 3L])))
  excluded <- ifelse(models_ind == 1L, total_indicator == 2L, total_indicator == 1L)
  # posterior model probability 1/3 times the within-model zero frequency
  expected_mass <- (mean(excluded[models_ind == 1L]) + mean(excluded[models_ind == 2L]) + 1) / 3

  levels <- marginal_posterior(mixed, "mu_f", formula = ~ f, prior_samples = TRUE)
  for(level in c("mid", "high")){
    level_posterior <- levels[[level]]
    atoms <- .bt_meta_get(level_posterior, "atoms")
    expect_equal(unname(atoms$locations[, 1]), .3)
    expect_equal(atoms$mass, expected_mass, tolerance = 1e-12)
    # the observed share of draws at .3 differs only through the multinomial
    # model counts (|difference| ~ 4e-4 here)
    expect_lt(abs(atoms$mass - mean(level_posterior == .3)), 2e-3)
    expect_equal(prior_density_ordinate(.bt_meta_get(level_posterior, "prior_density"), .3)$point_mass, 2 / 3,
                 tolerance = 1e-12)
    # prior-only draws: the Savage-Dickey ratio is 1 up to KDE error (10,000
    # continuous draws: relative sd ~2%, bias ~1-2%; |log BF| < 0.1)
    class(level_posterior) <- c(class(level_posterior), "marginal_posterior")
    for(null_hypothesis in c(.35, .8)){
      BF <- Savage_Dickey_BF(level_posterior, null_hypothesis = null_hypothesis, silent = TRUE)
      expect_lt(abs(log(as.numeric(BF))), .1)
    }
  }
})

test_that("ordered point(0) totals are structural zero coefficients", {

  alternative_prior <- .ordered_prior_for_test(prior("normal", list(0, 1)))
  null_prior <- .ordered_prior_for_test(prior("point", list(0)))
  set.seed(4)
  n <- 400
  total <- rnorm(n, .3, .1)
  alternative <- cbind("mu_f[1]" = .4 * total, "mu_f[2]" = .6 * total,mu_f_ordered_total=total)
  null <- cbind("mu_f[1]" = rep(0, n), "mu_f[2]" = rep(0, n))
  mixed <- mix_posteriors(
    list(
      list(fit = .mock_mixing_fit_for_marginal(alternative, list(mu_f = alternative_prior)),
           marglik = bridgesampling_object(0), prior_weights = 1),
      list(fit = .mock_mixing_fit_for_marginal(null, list(mu_f = null_prior)),
           marglik = bridgesampling_object(log(3)), prior_weights = 1)
    ),
    parameters = "mu_f", is_null_list = list(mu_f = c(FALSE, TRUE)),
    seed = 1, n_samples = n
  )
  atoms <- .bt_meta_get(mixed$mu_f, "atoms")
  expect_equal(unname(atoms$locations), matrix(0, 1, 2))
  expect_equal(atoms$mass, .75, tolerance = 1e-12)

  marginal <- marginal_posterior(mixed, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  # level mid = .4 * total and level high = total, each with prior null mass 1/2
  for(level in c("mid", "high")){
    scale <- if(level == "mid") .4 else 1
    prior_density <- .bt_meta_get(marginal[[level]], "prior_density")
    expect_equal(prior_density_ordinate(prior_density, 0)$point_mass, .5, tolerance = 1e-12)
    .expect_prior_height_for_test(marginal[[level]], .2, .5 * stats::dnorm(.2, 0, scale))
    expect_equal(.bt_meta_get(marginal[[level]], "atoms")$mass, .75, tolerance = 1e-12)
  }

  fit <- coda::mcmc(null)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(mu_f = null_prior)
  fit <- attach_test_parameter_map(fit)
  single <- as_mixed_posteriors(fit, parameters = "mu_f")
  expect_equal(.bt_meta_get(single$mu_f, "atoms")$mass, 1)
  single_marginal <- marginal_posterior(single, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  expect_equal(
    prior_density_ordinate(.bt_meta_get(single_marginal[["high"]], "prior_density"), 0)$point_mass,
    1
  )
})

test_that("Savage-Dickey posterior ordinates are exact reflected kernel sums", {

  # The posterior ordinate at the null is the Gaussian kernel sum of the
  # continuous draws with bandwidth bw.nrd0, reflected at exact support
  # bounds. The 512-point binned (unbounded) or gridded (bounded) estimate it
  # replaces spans the whole draw range; with long right tails its spacing
  # exceeds the bandwidth, and it was +43% on average (up to +272%) for the
  # lognormal posterior below at 1 with 1e4 draws, and +137% (up to +630%) for
  # the inverse-gamma tail. References: the kernel sum written out, and the
  # mean and variance of the reflected KDE under the true density by
  # quadrature (.reflected_kde_moments_for_test()); the sums lie within 4 of
  # their Monte Carlo SDs of that mean, and the mean itself is within 3% of
  # the true density.
  cases <- list(
    list(prior = prior("lognormal", list(0, 1.5)), draw = function(n) exp(stats::rnorm(n, 0, 1.5)),
         density = function(y) stats::dlnorm(y, 0, 1.5), null = 1, seed = 11),
    list(prior = prior("invgamma", list(1.5, 1)), draw = function(n) 1 / stats::rgamma(n, 1.5, 1),
         density = function(y) ifelse(y > 0, stats::dgamma(1 / y, 1.5, 1) / y^2, 0), null = .5, seed = 12)
  )
  for(case in cases){
    set.seed(case$seed)
    draws <- case$draw(6e4)
    prior_density <- BayesTools:::.prior_linear_combination_density(
      prior_list = list(theta = case$prior), weights = c(theta = 1)
    )
    posterior <- .marginal_posterior_with_prior_density_for_test(draws, prior_density)
    posterior <- .bt_meta_set(posterior, "support", .posterior_support_new(c(0, Inf)))
    bandwidth <- stats::bw.nrd0(draws)
    kernel_sum <- mean(stats::dnorm(case$null, draws, bandwidth)) +
      mean(stats::dnorm(case$null, -draws, bandwidth))
    bf <- Savage_Dickey_BF(posterior, null_hypothesis = case$null, silent = TRUE)
    prior_height <- as.numeric(BayesTools:::.prior_linear_density_height(prior_density, case$null))
    expect_equal(prior_height / as.numeric(bf), kernel_sum, tolerance = 1e-12)
    expect_true(attr(bf, "posterior_density_boundary_reflection"))
    moments <- .reflected_kde_moments_for_test(
      case$null, h = bandwidth, density = case$density, lower = 0, n = length(draws), n_source = Inf
    )
    expect_lte(abs(kernel_sum - moments[["mean"]]), 4 * sqrt(moments[["variance"]]))
    expect_lt(abs(moments[["mean"]] / case$density(case$null) - 1), .03)

    # the raw-draw hypothesis height is the unreflected kernel sum
    expect_equal(
      as.numeric(BayesTools:::.hypothesis_sample_density_height(draws, case$null, "posterior")),
      mean(stats::dnorm(case$null, draws, bandwidth)),
      tolerance = 1e-12
    )
  }

  # Draws about 100 bandwidths above the support bound 0: the reflected
  # kernel sum at the bound (log-sum-exp reference) is below the double range,
  # so the ordinate is 0 and the Bayes factor +Inf.
  set.seed(13)
  draws <- stats::rnorm(1e4, .5, .02)
  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(0, 1), list(0, Inf))), weights = c(theta = 1)
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(draws, prior_density)
  posterior <- .bt_meta_set(posterior, "support", .posterior_support_new(c(0, Inf)))
  bandwidth <- stats::bw.nrd0(draws)
  log_kernels <- stats::dnorm(draws / bandwidth, log = TRUE)
  log_ordinate <- log(2) + max(log_kernels) + log(sum(exp(log_kernels - max(log_kernels)))) -
    log(length(draws) * bandwidth)
  expect_lt(log_ordinate, log(.Machine$double.xmin))
  bf <- Savage_Dickey_BF(posterior, null_hypothesis = 0, silent = TRUE)
  expect_identical(as.numeric(bf), Inf)
})

test_that("Savage-Dickey prior support is the exact support, not the grid range", {

  # The support of the prior measure is exact, not the range of its
  # numerical grid: a normal prior's grid ends near +/-3.7 SD (tail
  # probability 1e-4), but a null at 6 lies inside its support.
  set.seed(14)
  normal <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(0, 1))), weights = c(theta = 1)
  )
  expect_lt(max(normal$density$x), 6)
  posterior <- .marginal_posterior_with_prior_density_for_test(stats::rnorm(2e4, 6, .5), normal)
  bf <- Savage_Dickey_BF(posterior, null_hypothesis = 6, silent = TRUE)
  expect_null(attr(bf, "warnings"))

  # a null outside the support of a half-normal prior has a zero prior
  # ordinate: the density ratio is undefined and the call stops (the former
  # "does not span both sides" warning came with a Bayes factor of 0)
  half_normal <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(0, 1), list(0, Inf))), weights = c(theta = 1)
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(stats::rnorm(2e4, -1, .5), half_normal)
  condition <- tryCatch(Savage_Dickey_BF(posterior, null_hypothesis = -1, silent = TRUE),
                        error = function(e) e)
  expect_s3_class(condition, "BayesTools_zero_ordinate")
})

test_that("N07 point slab atoms use the canonical component indices", {

  slab <- prior_spike_and_slab(prior("point", list(3)))
  atoms <- BayesTools:::.posterior_atoms_from_components(
    slab, c(1L, 1L, 1L, 2L), 1L, "x"
  )
  ordering <- order(atoms$locations[, 1L])
  expect_identical(as.numeric(atoms$locations[ordering, 1L]), c(0, 3))
  expect_equal(atoms$mass[ordering], c(.25, .75), tolerance = 0)
  expect_true(atoms$declared)
})

test_that("N07 all-included point slab retains its declared atom", {

  atoms <- BayesTools:::.posterior_atoms_from_components(
    prior_spike_and_slab(prior("point", list(3))), rep(1L, 4), 1L, "x"
  )
  expect_identical(as.numeric(atoms$locations), 3)
  expect_identical(atoms$mass, 1)
})

test_that("N07 continuous slab has only the excluded atom and zero slabs coalesce", {

  component <- c(1L, 1L, 1L, 2L)
  continuous <- BayesTools:::.posterior_atoms_from_components(
    prior_spike_and_slab(prior("normal", list(0, 1))), component, 1L, "x"
  )
  expect_identical(as.numeric(continuous$locations), 0)
  expect_identical(continuous$mass, .25)
  zero <- BayesTools:::.posterior_atoms_from_components(
    prior_spike_and_slab(prior("point", list(0))), component, 1L, "x"
  )
  expect_identical(as.numeric(zero$locations), 0)
  expect_identical(zero$mass, 1)
})

test_that("N07 missing component provenance is refused rather than reconstructed", {

  expect_error(
    BayesTools:::.posterior_atoms_from_components(
      prior_spike_and_slab(prior("point", list(3))), NULL, 1L, "x"
    ),
    "Posterior component indices must index the components of the mixture or spike-and-slab prior.",
    fixed = TRUE
  )
})


test_that("N08 nearby unequal simplex coefficients retain their exact interval", {

  priors <- list(w = prior("dirichlet", list(alpha = c(2, 3))))
  support <- BayesTools:::.posterior_support_from_prior_list_weights(
    priors, c("w[1]" = 1, "w[2]" = 1 + 1e-9)
  )
  expect_identical(support$type, "interval")
  expect_identical(support$bounds, c(1, 1 + 1e-9))
  expect_true(BayesTools:::.posterior_support_contains_value(support, 1 + 5e-10))
  expect_true(support$exact)
})

test_that("N08 equal simplex coefficients and omitted coefficients keep their meaning", {

  priors <- list(w = prior("dirichlet", list(alpha = c(2, 3))))
  equal <- BayesTools:::.posterior_support_from_prior_list_weights(
    priors, c("w[1]" = 1, "w[2]" = 1)
  )
  expect_identical(equal$type, "points")
  expect_identical(equal$points, 1)
  zero <- BayesTools:::.posterior_support_from_prior_list_weights(
    priors, c("w[1]" = 0, "w[2]" = 0)
  )
  expect_identical(zero$points, 0)
  omitted <- BayesTools:::.posterior_support_from_prior_list_weights(priors, c("w[2]" = 5))
  expect_identical(omitted$type, "interval")
  expect_identical(omitted$bounds, c(0, 5))
})

# File-level skips: All remaining tests in this file require pre-fitted models
skip_if_not_visual_fixture_tests()
skip_if_no_fits()
skip_if_not_installed("rjags")
skip_if_not_installed("bridgesampling")

test_that("Marginal distribution prior and posterior functions work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  set.seed(1)

  # Load pre-fitted marginal distribution models
  fit0     <- readRDS(file.path(temp_fits_dir, "fit_marginal_0.RDS"))
  fit1     <- readRDS(file.path(temp_fits_dir, "fit_marginal_1.RDS"))
  marglik0 <- readRDS(file.path(temp_marglik_dir, "fit_marginal_0.RDS"))
  marglik1 <- readRDS(file.path(temp_marglik_dir, "fit_marginal_1.RDS"))

  # Define prior lists (needed for manual mixing validation and prior densities)
  prior_list_0 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("spike", contrast = "treatment", list(0)),
    "x_fac3md"         = prior_factor("spike", contrast = "meandif",   list(0)),
    "x_cont1:x_fac3md" = prior_factor("spike", contrast = "meandif",   list(0))
  )
  prior_list_1 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("normal",  contrast = "treatment", list(0, 1.00)),
    "x_fac3md"         = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25)),
    "x_cont1:x_fac3md" = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25))
  )
  prior_list <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(prior_list_0$x_cont1, "multiply_by") <- "sigma"
  attr(prior_list_1$x_cont1, "multiply_by") <- "sigma"

  # make the mixing equal
  marglik1$logml <- marglik0$logml

  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    is_null_list = list(
      "sigma"                     = c(FALSE, FALSE),
      "mu_intercept"              = c(FALSE, FALSE),
      "mu_x_cont1"                = c(FALSE, FALSE),
      "mu_x_fac2t"                = c(TRUE, FALSE),
      "mu_x_fac3md"               = c(TRUE, FALSE),
      "mu_x_cont1__xXx__x_fac3md" = c(TRUE, FALSE)
    ),
    conditional  = FALSE)
  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    is_null_list = list(
      "sigma"                     = c(FALSE, FALSE),
      "mu_intercept"              = c(FALSE, FALSE),
      "mu_x_cont1"                = c(FALSE, FALSE),
      "mu_x_fac2t"                = c(TRUE, FALSE),
      "mu_x_fac3md"               = c(TRUE, FALSE),
      "mu_x_cont1__xXx__x_fac3md" = c(TRUE, FALSE)
    ),
    seed         = 1,
    conditional  = FALSE
  )

  # manual mixing
  posterior_manual0 <- suppressWarnings(coda::as.mcmc(fit0))
  posterior_manual1 <- suppressWarnings(coda::as.mcmc(fit1))
  add_missing_null_columns <- function(posterior, columns) {
    posterior <- as.matrix(posterior)
    missing <- setdiff(columns, colnames(posterior))
    if (length(missing) > 0L) {
      posterior <- cbind(
        posterior,
        matrix(
          0,
          nrow = nrow(posterior),
          ncol = length(missing),
          dimnames = list(NULL, missing)
        )
      )
    }
    posterior
  }
  null_formula_columns <- c(
    "mu_x_fac2t",
    "mu_x_fac3md[1]",
    "mu_x_fac3md[2]",
    "mu_x_cont1__xXx__x_fac3md[1]",
    "mu_x_cont1__xXx__x_fac3md[2]"
  )
  posterior_manual0 <- add_missing_null_columns(posterior_manual0, null_formula_columns)
  posterior_manual1 <- add_missing_null_columns(posterior_manual1, null_formula_columns)

  ### test error checks ----
  expect_error(marginal_posterior(
    samples           = list(posterior_manual0),
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md),
    "'samples' must be a be an object generated by 'mix_posteriors' function.")
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = "~ x_cont1 + x_fac2t + x_cont1*x_fac3md"),
    "'formula' must be a formula")
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = c(x_fac2t = NA),
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md),
    "'at' must be a list")
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + not_here),
    "The posterior samples for the 'mu_not_here' term is missing in the samples.")
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "not_here",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "The 'not_here' values are not recognized by the 'parameter' argument.")
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(mu_x_cont1 = 1),
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "The following values passed via the 'at' argument do not correspond to the specified model: 'mu_x_cont1'"
    )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(x_cont1 = 1),
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "Values of the parameter of interested cannot be specified via the 'at' argument."
  )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(x_fac2t = "D"),
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "Levels specified in the 'x_fac2t' factor variable do not match the levels used for model specification."
  )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(x_fac2t = NA),
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "Unspecified levels in the 'x_fac2t' factor",
  )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(x_cont1 = "A"),
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "Nonnumeric values in the 'x_cont1' continuous variable."
  )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(x_cont1 = NA),
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE),
    "Unspecified levels in the 'x_cont1' variable"
  )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    use_formula       = FALSE),
    "'formula' is supposed to be NULL when dealing with simple posteriors"
  )
  expect_error(marginal_posterior(
    samples           = mixed_posteriors,
    at                = list(x_cont1 = NA),
    parameter         = "mu_x_fac2t",
    use_formula       = FALSE),
    "'at' is supposed to be NULL when dealing with simple posteriors"
  )


  ### simple: continuous parameter ----
  marg_post_sigma <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "sigma",
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-simple-con", function(){
    hist(marg_post_sigma, freq = FALSE, main = "marginal posterior sigma")
    lines(density(c(posterior_manual0[,"sigma"], posterior_manual1[,"sigma"])))
  })

  vdiffr::expect_doppelganger("marginal-simple-con-p", function(){
    .plot_prior_density_for_test(marg_post_sigma, main = "marginal prior sigma")
    lines(density(prior_list$sigma))
  })


  ### simple: factor ----
  marg_post_simple_x_fac2t <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    prior_samples     = TRUE,
    use_formula       = FALSE)

  vdiffr::expect_doppelganger("marginal-simple-fac", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(marg_post_simple_x_fac2t[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")

    hist(marg_post_simple_x_fac2t[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_x_fac2t"], posterior_manual1[,"mu_x_fac2t"])))

  })

  vdiffr::expect_doppelganger("marginal-simple-fac-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    .plot_prior_density_for_test(marg_post_simple_x_fac2t[["A"]], main = "marg_post_x_fac2t = A")

    .plot_prior_density_for_test(marg_post_simple_x_fac2t[["B"]], main = "marg_post_x_fac2t = B")
    curve(dnorm(x, 0, 1)/2, add = T)

  })


  ### formula: intercept ----
  marg_post_int <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_intercept",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-form-int", function(){
    hist(marg_post_int[["intercept"]], freq = FALSE, main = "marginal posterior intercept")
    lines(density(c(posterior_manual0[,"mu_intercept"], posterior_manual1[,"mu_intercept"] )))
  })

  vdiffr::expect_doppelganger("marginal-form-int-p", function(){
    .plot_prior_density_for_test(marg_post_int[["intercept"]], main = "marginal prior intercept")
    lines(prior_list_0$intercept)
  })


  ### formula: continuous parameter (-+1SD) ----
  marg_post_x_cont1 <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-form-con", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_cont1[["-1SD"]], freq = FALSE, main = "marginal posterior x_cont1\n(-1)")
    lines(density(c(posterior_manual0[,"mu_intercept"] + -1 * posterior_manual0[,"mu_x_cont1"] * posterior_manual0[,"sigma"],
                    posterior_manual1[,"mu_intercept"] + -1 * posterior_manual1[,"mu_x_cont1"] * posterior_manual1[,"sigma"])))

    hist(marg_post_x_cont1[["0SD"]], freq = FALSE, main = "marginal posterior x_cont1\n(0)")
    lines(density(c(posterior_manual0[,"mu_intercept"] + 0 * posterior_manual0[,"mu_x_cont1"] * posterior_manual0[,"sigma"],
                    posterior_manual1[,"mu_intercept"] + 0 * posterior_manual1[,"mu_x_cont1"] * posterior_manual1[,"sigma"])))

    hist(marg_post_x_cont1[["1SD"]], freq = FALSE, main = "marginal posterior x_cont1\n(1)")
    lines(density(c(posterior_manual0[,"mu_intercept"] + 1 * posterior_manual0[,"mu_x_cont1"] * posterior_manual0[,"sigma"],
                    posterior_manual1[,"mu_intercept"] + 1 * posterior_manual1[,"mu_x_cont1"] * posterior_manual1[,"sigma"])))

  })

  vdiffr::expect_doppelganger("marginal-form-con-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_cont1[["-1SD"]], main = "marginal prior x_cont1\n(-1)", xlim = c(-10, 10))
    .plot_prior_density_for_test(marg_post_x_cont1[["0SD"]],  main = "marginal prior x_cont1\n(0)",  xlim = c(-10, 10))
    .plot_prior_density_for_test(marg_post_x_cont1[["1SD"]],  main = "marginal prior x_cont1\n(1)",  xlim = c(-10, 10))

  })


  ### formula: treatment factor ----
  marg_post_x_fac2t <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-form-fac.t", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(marg_post_x_fac2t[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")
    lines(density(c(posterior_manual0[,"mu_intercept"], posterior_manual1[,"mu_intercept"])))

    hist(marg_post_x_fac2t[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0[,"mu_x_fac2t"], posterior_manual1[,"mu_intercept"] + posterior_manual1[,"mu_x_fac2t"])))

  })

  vdiffr::expect_doppelganger("marginal-form-fac.t-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    .plot_prior_density_for_test(marg_post_x_fac2t[["A"]], main = "marginal prior x_fac2t = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x), add = TRUE)
    .plot_prior_density_for_test(marg_post_x_fac2t[["B"]], main = "marginal prior x_fac2t = B",  xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 1^2)) / 2), add = TRUE)
  })


  ### formula: meandif factor ----
  marg_post_x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  posterior_manual0.md <- posterior_manual0[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.md <- posterior_manual1[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))

  vdiffr::expect_doppelganger("marginal-form-fac.md", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_fac3md[["A"]], freq = FALSE, main = "marg_post_x_fac3md = A", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1])))

    hist(marg_post_x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,2], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,2])))

    hist(marg_post_x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,3], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,3])))
  })

  vdiffr::expect_doppelganger("marginal-form-fac.md-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_fac3md[["A"]], main = "marginal prior x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 0.25^2)) / 2), add = TRUE)
    .plot_prior_density_for_test(marg_post_x_fac3md[["B"]], main = "marginal prior x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 0.25^2)) / 2), add = TRUE)
    .plot_prior_density_for_test(marg_post_x_fac3md[["C"]], main = "marginal prior x_fac3md = C", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 0.25^2)) / 2), add = TRUE)
  })


  ### formula: meandif factor interaction ----
  marg_post_x_cont1__xXx__x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1__xXx__x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  posterior_manual0.md <- matrix(posterior_manual0[, "mu_intercept"], ncol = 9, nrow = nrow(posterior_manual0)) +
    posterior_manual0[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% do.call(cbind, lapply(1:3, function(i) t(contr.meandif(1:3)))) +
    matrix(posterior_manual0[, "sigma"], ncol = 9, nrow = nrow(posterior_manual0)) * posterior_manual0[, "mu_x_cont1"] %*% t(c(-1, -1, -1, 0, 0, 0, 1, 1, 1)) +
    posterior_manual0[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% (do.call(cbind, lapply(1:3, function(i) t(contr.meandif(1:3)))) * matrix(c(-1, -1, -1, 0, 0, 0, 1, 1, 1), ncol = 9, nrow = 2, byrow = TRUE))
  posterior_manual1.md <- matrix(posterior_manual1[, "mu_intercept"], ncol = 9, nrow = nrow(posterior_manual1)) +
    posterior_manual1[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% do.call(cbind, lapply(1:3, function(i) t(contr.meandif(1:3)))) +
    matrix(posterior_manual1[, "sigma"], ncol = 9, nrow = nrow(posterior_manual1)) * posterior_manual1[, "mu_x_cont1"] %*% t(c(-1, -1, -1, 0, 0, 0, 1, 1, 1)) +
    posterior_manual1[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% (do.call(cbind, lapply(1:3, function(i) t(contr.meandif(1:3)))) * matrix(c(-1, -1, -1, 0, 0, 0, 1, 1, 1), ncol = 9, nrow = 2, byrow = TRUE))
  posterior_manual.md <- rbind(posterior_manual0.md, posterior_manual1.md)

  vdiffr::expect_doppelganger("marginal-form-fac.mdi", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(3, 3))
    hist(marg_post_x_cont1__xXx__x_fac3md[["-1SD, A"]], freq = FALSE, main = "x_cont1 = -1\nmarg_post_x_fac3md = A", breaks = 20)
    lines(density(posterior_manual.md[,1]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["-1SD, B"]], freq = FALSE, main = "x_cont1 = -1\nmarg_post_x_fac3md = B", breaks = 20)
    lines(density(posterior_manual.md[,2]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["-1SD, C"]], freq = FALSE, main = "x_cont1 = -1\nmarg_post_x_fac3md = B", breaks = 20)
    lines(density(posterior_manual.md[,3]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["0SD, A"]], freq = FALSE, main = "x_cont1 = 0\nmarg_post_x_fac3md = A", breaks = 20)
    lines(density(posterior_manual.md[,4]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["0SD, B"]], freq = FALSE, main = "x_cont1 = 0\nmarg_post_x_fac3md = B", breaks = 20)
    lines(density(posterior_manual.md[,5]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["0SD, C"]], freq = FALSE, main = "x_cont1 = 0\nmarg_post_x_fac3md = B", breaks = 20)
    lines(density(posterior_manual.md[,6]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["1SD, A"]], freq = FALSE, main = "x_cont1 = 1\nmarg_post_x_fac3md = A", breaks = 20)
    lines(density(posterior_manual.md[,7]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["1SD, B"]], freq = FALSE, main = "x_cont1 = 1\nmarg_post_x_fac3md = B", breaks = 20)
    lines(density(posterior_manual.md[,8]))

    hist(marg_post_x_cont1__xXx__x_fac3md[["1SD, C"]], freq = FALSE, main = "x_cont1 = 1\nmarg_post_x_fac3md = B", breaks = 20)
    lines(density(posterior_manual.md[,9]))

  })

  vdiffr::expect_doppelganger("marginal-form-fac.mdi-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(3, 3))
    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["-1SD, A"]], main = "x_cont1 = -1\nmarg_post_x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["-1SD, B"]], main = "x_cont1 = -1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["-1SD, C"]], main = "x_cont1 = -1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["0SD, A"]], main = "x_cont1 = 0\nmarg_post_x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 0 + 0.25^2 + 0)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["0SD, B"]], main = "x_cont1 = 0\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 0 + 0.25^2 + 0)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["0SD, C"]], main = "x_cont1 = 0\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 0 + 0.25^2 + 0)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["1SD, A"]], main = "x_cont1 = 1\nmarg_post_x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["1SD, B"]], main = "x_cont1 = 1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["1SD, C"]], main = "x_cont1 = 1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)) , add = TRUE)

  })


  ### formula: meandif factor + at specification ----
  marg_post_x_fac3md_AT <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    at                = list(
      x_cont1 = 1,
      x_fac2t = c("A", "B")
    ),
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  posterior_manual0.md  <- posterior_manual0[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.md  <- posterior_manual1[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual0.mdi <- posterior_manual0[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.mdi <- posterior_manual1[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% t(contr.meandif(1:3))

  vdiffr::expect_doppelganger("marginal-form-fac.md-at", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(3, 2))
    hist(marg_post_x_fac3md_AT[["A"]][1,], freq = FALSE, main = "marg_post_x_fac3md = A | 1,A", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1] + posterior_manual0[,"sigma"] * posterior_manual0[,"mu_x_cont1"] + posterior_manual0.mdi[,1],
                    posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1] + posterior_manual1[,"sigma"] * posterior_manual1[,"mu_x_cont1"] + posterior_manual1.mdi[,1])))

    hist(marg_post_x_fac3md_AT[["A"]][2,], freq = FALSE, main = "marg_post_x_fac3md = A | 1,B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1] + posterior_manual0[,"sigma"] * posterior_manual0[,"mu_x_cont1"] + posterior_manual0[,"mu_x_fac2t"] + posterior_manual0.mdi[,1],
                    posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1] + posterior_manual1[,"sigma"] * posterior_manual1[,"mu_x_cont1"] + posterior_manual1[,"mu_x_fac2t"] + posterior_manual1.mdi[,1])))

    hist(marg_post_x_fac3md_AT[["B"]][1,], freq = FALSE, main = "marg_post_x_fac3md = B | 1,A", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,2] + posterior_manual0[,"sigma"] * posterior_manual0[,"mu_x_cont1"] + posterior_manual0.mdi[,2],
                    posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,2] + posterior_manual1[,"sigma"] * posterior_manual1[,"mu_x_cont1"] + posterior_manual1.mdi[,2])))

    hist(marg_post_x_fac3md_AT[["B"]][2,], freq = FALSE, main = "marg_post_x_fac3md = B | 1,B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,2] + posterior_manual0[,"sigma"] * posterior_manual0[,"mu_x_cont1"] + posterior_manual0[,"mu_x_fac2t"] + posterior_manual0.mdi[,2],
                    posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,2] + posterior_manual1[,"sigma"] * posterior_manual1[,"mu_x_cont1"] + posterior_manual1[,"mu_x_fac2t"] + posterior_manual1.mdi[,2])))

    hist(marg_post_x_fac3md_AT[["C"]][1,], freq = FALSE, main = "marg_post_x_fac3md = C | 1,A", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,3] + posterior_manual0[,"sigma"] * posterior_manual0[,"mu_x_cont1"] + posterior_manual0.mdi[,3],
                    posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,3] + posterior_manual1[,"sigma"] * posterior_manual1[,"mu_x_cont1"] + posterior_manual1.mdi[,3])))

    hist(marg_post_x_fac3md_AT[["C"]][2,], freq = FALSE, main = "marg_post_x_fac3md = C | 1,B", breaks = 20)
    lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,3] + posterior_manual0[,"sigma"] * posterior_manual0[,"mu_x_cont1"] + posterior_manual0[,"mu_x_fac2t"] + posterior_manual0.mdi[,3],
                    posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,3] + posterior_manual1[,"sigma"] * posterior_manual1[,"mu_x_cont1"] + posterior_manual1[,"mu_x_fac2t"] + posterior_manual1.mdi[,3])))

  })

  ### formula: transformation ----
  marg_post_x_cont1.exp <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    transformation    = "exp",
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-form-con-exp", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_cont1.exp[["-1SD"]], freq = FALSE, main = "exp marginal posterior x_cont1\n(-1)")
    lines(density(exp(marg_post_x_cont1[["-1SD"]])))

    hist(marg_post_x_cont1.exp[["0SD"]], freq = FALSE, main = "exp marginal posterior x_cont1\n(0)")
    lines(density(exp(marg_post_x_cont1[["0SD"]])))

    hist(marg_post_x_cont1.exp[["1SD"]], freq = FALSE, main = "exp marginal posterior x_cont1\n(1)")
    lines(density(exp(marg_post_x_cont1[["1SD"]])))

  })

  vdiffr::expect_doppelganger("marginal-form-con-p-exp", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_cont1.exp[["-1SD"]], main = "marginal prior x_cont1\n(-1)", xlim = c(0, 10))
    .plot_prior_density_for_test(marg_post_x_cont1.exp[["0SD"]],  main = "marginal prior x_cont1\n(0)",  xlim = c(0, 10))
    .plot_prior_density_for_test(marg_post_x_cont1.exp[["1SD"]],  main = "marginal prior x_cont1\n(1)",  xlim = c(0, 10))
  })

  ### Savage-Dickey BFs ----
  # Smoke tests on mixed_posteriors (model-averaged). Savage-Dickey for a
  # spike-and-slab / mixture parameter is defined on the continuous
  # (conditional-on-inclusion) posterior; use as_marginal_inference() for that.
  # test input
  expect_error(Savage_Dickey_BF(list(posterior_manual0)), "'Savage_Dickey_BF' requires an object of class 'marginal_posterior'.")
  expect_error(Savage_Dickey_BF(marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "sigma",
    prior_samples     = FALSE)), "there are no prior densities for the posterior distribution")

  # simple restricted prior: the posterior draws of sigma lie about 108
  # bandwidths above the null (and support bound) 0, so the exact reflected
  # kernel sum at 0, 2 sum_i phi(x_i / h) / (n h), is about exp(-5828) (the
  # log-sum-exp reference below) and underflows: the Bayes factor exceeds the
  # double range. (The earlier 512-point grid KDE reported 1.5e15.)
  BF.marg_post_sigma <- Savage_Dickey_BF(marg_post_sigma, silent = TRUE)
  sigma_draws <- as.numeric(marg_post_sigma)
  bandwidth <- stats::bw.nrd0(sigma_draws)
  log_kernels <- stats::dnorm(sigma_draws / bandwidth, log = TRUE)
  log_ordinate <- log(2) + max(log_kernels) + log(sum(exp(log_kernels - max(log_kernels)))) -
    log(length(sigma_draws) * bandwidth)
  expect_lt(log_ordinate, log(.Machine$double.xmin))
  expect_identical(as.numeric(BF.marg_post_sigma), Inf)
  expect_null(attr(BF.marg_post_sigma, "warnings"))
  expect_true(isTRUE(attr(BF.marg_post_sigma, "posterior_density_boundary_reflection")))
  expect_equal(attr(BF.marg_post_sigma, "posterior_density_support"), c(0, 5))

  # simple factor: the reference level A is fixed at 0 and level B has the
  # null model's point mass at 0, so both levels are NA with their reasons
  BF.marg_post_simple_x_fac2t <- Savage_Dickey_BF(marg_post_simple_x_fac2t, silent = TRUE)
  expect_true(all(is.na(unlist(BF.marg_post_simple_x_fac2t))))
  expect_identical(
    attr(BF.marg_post_simple_x_fac2t[["A"]], "warnings"),
    "The posterior is fixed at the null hypothesis value. The Savage-Dickey Bayes factor is undefined."
  )
  expect_identical(
    attr(BF.marg_post_simple_x_fac2t[["B"]], "warnings"),
    "The posterior contains a declared point mass at the exact null hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
  )
  marg_post_simple_x_fac2t_B <- marg_post_simple_x_fac2t[["B"]]
  class(marg_post_simple_x_fac2t_B) <- c(class(marg_post_simple_x_fac2t_B), "marginal_posterior")
  expect_error(
    Savage_Dickey_BF(marg_post_simple_x_fac2t_B),
    "exact null hypothesis value",
    fixed = TRUE
  )


  BF.marg_post_x_fac3md <- Savage_Dickey_BF(marg_post_x_fac3md, silent = TRUE)
  BF.marg_post_x_fac3md_values <- unlist(BF.marg_post_x_fac3md, use.names = FALSE)
  expect_true(all(is.finite(BF.marg_post_x_fac3md_values)))
  expect_gt(min(BF.marg_post_x_fac3md_values), 1e50)
  expect_equal(attr(BF.marg_post_x_fac3md[["A"]], "warnings"),
               "Posterior samples do not span both sides of the null hypothesis. The posterior density at the null hypothesis is an extrapolation from Gaussian kernel tails; the Bayes factor is not reliable evidence.")

  BF2.marg_post_x_fac3md <- suppressWarnings(Savage_Dickey_BF(marg_post_x_fac3md, null_hypothesis = 0.5))
  BF2.marg_post_x_fac3md_values <- unlist(BF2.marg_post_x_fac3md, use.names = FALSE)
  expect_true(all(is.finite(BF2.marg_post_x_fac3md_values)))
  expect_true(all(BF2.marg_post_x_fac3md_values > 0))
  expect_gt(BF2.marg_post_x_fac3md_values[1], 1)
  expect_true(all(BF2.marg_post_x_fac3md_values[-1] < 1))

  BF2_normal.marg_post_x_fac3md <- suppressWarnings(Savage_Dickey_BF(marg_post_x_fac3md, null_hypothesis = 0.5, normal_approximation = TRUE))
  BF2_normal.marg_post_x_fac3md_values <- unlist(
    BF2_normal.marg_post_x_fac3md,
    use.names = FALSE
  )
  expect_true(all(is.finite(BF2_normal.marg_post_x_fac3md_values)))
  expect_true(all(BF2_normal.marg_post_x_fac3md_values > 0))
  expect_true(all(BF2_normal.marg_post_x_fac3md_values < 1))
  expect_true(all(
    BF2.marg_post_x_fac3md_values >
      BF2_normal.marg_post_x_fac3md_values
  ))


  ### marginal_inference ----
  set.seed(1)
  out <- marginal_inference(
    model_list          = models,
    marginal_parameters = c("mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    parameters          = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    is_null_list        = list(
      "sigma"                     = c(FALSE, FALSE),
      "mu_intercept"              = c(FALSE, FALSE),
      "mu_x_cont1"                = c(FALSE, FALSE),
      "mu_x_fac2t"                = c(TRUE, FALSE),
      "mu_x_fac3md"               = c(TRUE, FALSE),
      "mu_x_cont1__xXx__x_fac3md" = c(TRUE, FALSE)
    ),
    formula      =  ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    silent       = TRUE
  )

  # test samples against previously generated ones
  vdiffr::expect_doppelganger("marginal_inference-cont",     function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_cont1[["-1SD"]], freq = FALSE, main = "mu_x_cont1 = -1SD", breaks = 20)
    lines(density(out$averaged$mu_x_cont1[["-1SD"]]))

    hist(marg_post_x_cont1[["0SD"]], freq = FALSE, main = "mu_x_cont1 = 0SD", breaks = 20)
    lines(density(out$averaged$mu_x_cont1[["0SD"]]))

    hist(marg_post_x_cont1[["1SD"]], freq = FALSE, main = "mu_x_cont1 = +1SD", breaks = 20)
    lines(density(out$averaged$mu_x_cont1[["1SD"]]))

  })
  vdiffr::expect_doppelganger("marginal_inference-cont-p",   function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_cont1[["-1SD"]], main = "mu_x_cont1 = -1SD", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_cont1[["-1SD"]], add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_cont1[["0SD"]], main = "mu_x_cont1 = 0SD", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_cont1[["0SD"]], add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_cont1[["1SD"]], main = "mu_x_cont1 = 1SD", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_cont1[["1SD"]], add = TRUE, lty = 2)
  })
  vdiffr::expect_doppelganger("marginal_inference-fac.md",   function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_fac3md[["A"]], freq = FALSE, main = "marg_post_x_fac3md = A", breaks = 20)
    lines(density(out$averaged$mu_x_fac3md$A))

    hist(marg_post_x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
    lines(density(out$averaged$mu_x_fac3md$B))

    hist(marg_post_x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
    lines(density(out$averaged$mu_x_fac3md$C))

  })
  vdiffr::expect_doppelganger("marginal_inference-fac.md-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_fac3md[["A"]], main = "marginal prior x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_fac3md$A, add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_fac3md[["B"]], main = "marginal prior x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_fac3md$B, add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_fac3md[["C"]], main = "marginal prior x_fac3md = C", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_fac3md$C, add = TRUE, lty = 2)
  })
  # the previous BFs were based on model-averaged posteriors so they won't match

  # test summary table
  marginal_parameters <- c(
    "mu_intercept",
    "mu_x_cont1",
    "mu_x_fac2t",
    "mu_x_fac3md",
    "mu_x_cont1__xXx__x_fac3md"
  )
  marginal_table <- marginal_estimates_table(
    out$conditional,
    out$inference,
    parameters = marginal_parameters
  )
  test_reference_table_stochastic(
    marginal_table,
    "marginal_estimates_table_model_avg.txt",
    info_msg = "marginal_estimates_table for model averaging"
  )
  .expect_marginal_table_current_inputs(
    table = marginal_table,
    samples = out$conditional,
    inference = out$inference,
    parameters = marginal_parameters
  )

  # plots
  vdiffr::expect_doppelganger("plot_marginal-mu_x_fac2t-1", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t")})
  vdiffr::expect_doppelganger("plot_marginal-mu_x_fac2t-2", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", par_name = "fac2t", lwd = 2)})
  vdiffr::expect_doppelganger("plot_marginal-mu_x_fac2t-3", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2))})
  vdiffr::expect_doppelganger("plot_marginal-mu_x_fac2t-4", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1))})
  vdiffr::expect_doppelganger("plot_marginal-mu_x_fac2t-5", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2), transformation = "exp", xlim = c(0, 5), transformation_settings = T)})

  vdiffr::expect_doppelganger("ggplot_marginal-mu_x_fac2t-1", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t"))
  vdiffr::expect_doppelganger("ggplot_marginal-mu_x_fac2t-2", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t", par_name = "fac2t", lwd = 2))
  vdiffr::expect_doppelganger("ggplot_marginal-mu_x_fac2t-3", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2)))
  vdiffr::expect_doppelganger("ggplot_marginal-mu_x_fac2t-4", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1)))

  vdiffr::expect_doppelganger("plot_marginal-mu_x_cont1", function(){plot_marginal(out$conditional, parameter = "mu_x_cont1", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1))})
  vdiffr::expect_doppelganger("ggplot_marginal-mu_x_cont1", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_cont1", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1)))

  vdiffr::expect_doppelganger("plot_marginal-mu_x_fac3md", function(){plot_marginal(out$averaged, parameter = "mu_x_fac3md", prior = TRUE, dots_prior = list(lty = 2), xlim = c(-1, 1))})
  vdiffr::expect_doppelganger("ggplot_marginal-mu_x_fac3md", plot_marginal(out$averaged, plot_type = "ggplot", parameter = "mu_x_fac3md", prior = TRUE, dots_prior = list(lty = 2), xlim = c(-1, 1)))

  vdiffr::expect_doppelganger("plot_marginal-int", plot_marginal(out$averaged, plot_type = "ggplot", parameter = "mu_intercept", prior = TRUE, dots_prior = list(lty = 2), xlim = c(-1, 1)))

})

test_that("Marginal distributions with spike and slab and mixture priors work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Load pre-fitted spike-and-slab model
  fit <- readRDS(file.path(temp_fits_dir, "fit_marginal_ss.RDS"))

  # Define prior lists (needed for prior density validation in marginal_posterior)
  prior_pars <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior_mixture(list(
      prior("spike", list(0)),
      prior("normal", list(0, 1))
    ), is_null = c(T, F)),
    "x_fac2t"          = prior_spike_and_slab(prior_factor("normal",  contrast = "treatment", list(0, 1.00))),
    "x_fac3md"         = prior_spike_and_slab(prior_factor("mnormal", contrast = "meandif",   list(0, 0.25))),
    "x_cont1:x_fac3md" = prior_spike_and_slab(prior_factor("mnormal", contrast = "meandif",   list(0, 0.25)))
  )
  prior_list <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(prior_pars$x_cont1, "multiply_by") <- "sigma"

  mixed_posteriors <- as_mixed_posteriors(
    model        = fit,
    parameters   = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md")
  )

  ### simple: continuous parameter ----
  marg_post_sigma <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "sigma",
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-simple-con", function(){
    hist(marg_post_sigma, freq = FALSE, main = "marginal posterior sigma")
  })

  vdiffr::expect_doppelganger("marginal-ss-simple-con-p", function(){
    .plot_prior_density_for_test(marg_post_sigma, main = "marginal prior sigma")
    lines(density(prior_list$sigma))
  })


  ### simple: factor ----
  marg_post_simple_x_fac2t <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    prior_samples     = TRUE,
    use_formula       = FALSE)

  vdiffr::expect_doppelganger("marginal-ss-simple-fac", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(marg_post_simple_x_fac2t[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")

    hist(marg_post_simple_x_fac2t[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  })

  vdiffr::expect_doppelganger("marginal-ss-simple-fac-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    .plot_prior_density_for_test(marg_post_simple_x_fac2t[["A"]], main = "marg_post_x_fac2t = A")

    .plot_prior_density_for_test(marg_post_simple_x_fac2t[["B"]], main = "marg_post_x_fac2t = B")
    curve(dnorm(x, 0, 1)/2, add = T)

  })


  ### formula: intercept ----
  marg_post_int <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_intercept",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-form-int", function(){
    hist(marg_post_int[["intercept"]], freq = FALSE, main = "marginal posterior intercept")
  })

  vdiffr::expect_doppelganger("marginal-ss-form-int-p", function(){
    .plot_prior_density_for_test(marg_post_int[["intercept"]], main = "marginal prior intercept")
    lines(prior_pars$intercept)
  })


  ### formula: continuous parameter (-+1SD) ----
  marg_post_x_cont1 <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-form-con", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_cont1[["-1SD"]], freq = FALSE, main = "marginal posterior x_cont1\n(-1)")
    hist(marg_post_x_cont1[["0SD"]], freq = FALSE, main = "marginal posterior x_cont1\n(0)")
    hist(marg_post_x_cont1[["1SD"]], freq = FALSE, main = "marginal posterior x_cont1\n(1)")

  })

  vdiffr::expect_doppelganger("marginal-ss-form-con-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_cont1[["-1SD"]], main = "marginal prior x_cont1\n(-1)", xlim = c(-10, 10))
    .plot_prior_density_for_test(marg_post_x_cont1[["0SD"]],  main = "marginal prior x_cont1\n(0)",  xlim = c(-10, 10))
    .plot_prior_density_for_test(marg_post_x_cont1[["1SD"]],  main = "marginal prior x_cont1\n(1)",  xlim = c(-10, 10))

  })


  ### formula: treatment factor ----
  marg_post_x_fac2t <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-form-fac.t", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(marg_post_x_fac2t[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")
    hist(marg_post_x_fac2t[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)

  })

  vdiffr::expect_doppelganger("marginal-ss-form-fac.t-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    .plot_prior_density_for_test(marg_post_x_fac2t[["A"]], main = "marginal prior x_fac2t = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x), add = TRUE)
    .plot_prior_density_for_test(marg_post_x_fac2t[["B"]], main = "marginal prior x_fac2t = B",  xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 1^2)) / 2), add = TRUE)
  })


  ### formula: meandif factor ----
  marg_post_x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-form-fac.md", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_fac3md[["A"]], freq = FALSE, main = "marg_post_x_fac3md = A", breaks = 20)
    hist(marg_post_x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
    hist(marg_post_x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  })

  vdiffr::expect_doppelganger("marginal-ss-form-fac.md-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_fac3md[["A"]], main = "marginal prior x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 0.25^2)) / 2), add = TRUE)
    .plot_prior_density_for_test(marg_post_x_fac3md[["B"]], main = "marginal prior x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 0.25^2)) / 2), add = TRUE)
    .plot_prior_density_for_test(marg_post_x_fac3md[["C"]], main = "marginal prior x_fac3md = C", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, (sqrt(1^2 + 0^2) + sqrt(1^2 + 0.25^2)) / 2), add = TRUE)
  })


  ### formula: meandif factor interaction ----
  marg_post_x_cont1__xXx__x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1__xXx__x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-form-fac.mdi", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(3, 3))
    hist(marg_post_x_cont1__xXx__x_fac3md[["-1SD, A"]], freq = FALSE, main = "x_cont1 = -1\nmarg_post_x_fac3md = A", breaks = 20)
    hist(marg_post_x_cont1__xXx__x_fac3md[["-1SD, B"]], freq = FALSE, main = "x_cont1 = -1\nmarg_post_x_fac3md = B", breaks = 20)
    hist(marg_post_x_cont1__xXx__x_fac3md[["-1SD, C"]], freq = FALSE, main = "x_cont1 = -1\nmarg_post_x_fac3md = B", breaks = 20)

    hist(marg_post_x_cont1__xXx__x_fac3md[["0SD, A"]], freq = FALSE, main = "x_cont1 = 0\nmarg_post_x_fac3md = A", breaks = 20)
    hist(marg_post_x_cont1__xXx__x_fac3md[["0SD, B"]], freq = FALSE, main = "x_cont1 = 0\nmarg_post_x_fac3md = B", breaks = 20)
    hist(marg_post_x_cont1__xXx__x_fac3md[["0SD, C"]], freq = FALSE, main = "x_cont1 = 0\nmarg_post_x_fac3md = B", breaks = 20)

    hist(marg_post_x_cont1__xXx__x_fac3md[["1SD, A"]], freq = FALSE, main = "x_cont1 = 1\nmarg_post_x_fac3md = A", breaks = 20)
    hist(marg_post_x_cont1__xXx__x_fac3md[["1SD, B"]], freq = FALSE, main = "x_cont1 = 1\nmarg_post_x_fac3md = B", breaks = 20)
    hist(marg_post_x_cont1__xXx__x_fac3md[["1SD, C"]], freq = FALSE, main = "x_cont1 = 1\nmarg_post_x_fac3md = B", breaks = 20)

  })

  vdiffr::expect_doppelganger("marginal-ss-form-fac.mdi-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(3, 3))
    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["-1SD, A"]], main = "x_cont1 = -1\nmarg_post_x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(0.5*sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)^2 + 0.5*sqrt(0.25^2 + 0.25^2))), add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["-1SD, B"]], main = "x_cont1 = -1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(0.5*sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)^2 + 0.5*sqrt(0.25^2 + 0.25^2))), add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["-1SD, C"]], main = "x_cont1 = -1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(0.5*sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)^2 + 0.5*sqrt(0.25^2 + 0.25^2))), add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["0SD, A"]], main = "x_cont1 = 0\nmarg_post_x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 0 + 0.25^2 + 0)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["0SD, B"]], main = "x_cont1 = 0\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 0 + 0.25^2 + 0)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["0SD, C"]], main = "x_cont1 = 0\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(1^2 + 0 + 0.25^2 + 0)) , add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["1SD, A"]], main = "x_cont1 = 1\nmarg_post_x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(0.5*sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)^2 + 0.5*sqrt(0.25^2 + 0.25^2))), add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["1SD, B"]], main = "x_cont1 = 1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(0.5*sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)^2 + 0.5*sqrt(0.25^2 + 0.25^2))), add = TRUE)

    .plot_prior_density_for_test(marg_post_x_cont1__xXx__x_fac3md[["1SD, C"]], main = "x_cont1 = 1\nmarg_post_x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    curve(dnorm(x, 0, sqrt(0.5*sqrt(1^2 + 1^2 + 0.25^2 + 0.25^2)^2 + 0.5*sqrt(0.25^2 + 0.25^2))), add = TRUE)

  })


  ### formula: meandif factor + at specification ----
  marg_post_x_fac3md_AT <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    at                = list(
      x_cont1 = 1,
      x_fac2t = c("A", "B")
    ),
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)


  vdiffr::expect_doppelganger("marginal-ss-form-fac.md-at", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(3, 2))
    hist(marg_post_x_fac3md_AT[["A"]][1,], freq = FALSE, main = "marg_post_x_fac3md = A | 1,A", breaks = 20)

    hist(marg_post_x_fac3md_AT[["A"]][2,], freq = FALSE, main = "marg_post_x_fac3md = A | 1,B", breaks = 20)

    hist(marg_post_x_fac3md_AT[["B"]][1,], freq = FALSE, main = "marg_post_x_fac3md = B | 1,A", breaks = 20)

    hist(marg_post_x_fac3md_AT[["B"]][2,], freq = FALSE, main = "marg_post_x_fac3md = B | 1,B", breaks = 20)

    hist(marg_post_x_fac3md_AT[["C"]][1,], freq = FALSE, main = "marg_post_x_fac3md = C | 1,A", breaks = 20)

    hist(marg_post_x_fac3md_AT[["C"]][2,], freq = FALSE, main = "marg_post_x_fac3md = C | 1,B", breaks = 20)

  })

  ### formula: transformation ----
  marg_post_x_cont1.exp <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    transformation    = "exp",
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-form-con-exp", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_cont1.exp[["-1SD"]], freq = FALSE, main = "exp marginal posterior x_cont1\n(-1)")
    lines(density(exp(marg_post_x_cont1[["-1SD"]])))

    hist(marg_post_x_cont1.exp[["0SD"]], freq = FALSE, main = "exp marginal posterior x_cont1\n(0)")
    lines(density(exp(marg_post_x_cont1[["0SD"]])))

    hist(marg_post_x_cont1.exp[["1SD"]], freq = FALSE, main = "exp marginal posterior x_cont1\n(1)")
    lines(density(exp(marg_post_x_cont1[["1SD"]])))

  })

  vdiffr::expect_doppelganger("marginal-ss-form-con-p-exp", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_cont1.exp[["-1SD"]], main = "marginal prior x_cont1\n(-1)", xlim = c(0, 10))
    .plot_prior_density_for_test(marg_post_x_cont1.exp[["0SD"]],  main = "marginal prior x_cont1\n(0)",  xlim = c(0, 10))
    .plot_prior_density_for_test(marg_post_x_cont1.exp[["1SD"]],  main = "marginal prior x_cont1\n(1)",  xlim = c(0, 10))
  })

  ### conditional marginal samples ----
  mixed_posteriors <- as_mixed_posteriors(
    model        = fit,
    parameters   = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    conditional  = c("mu_x_cont1", "mu_x_fac3md"),
    conditional_rule = "AND"
  )
  marg_post_sigma <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  vdiffr::expect_doppelganger("marginal-ss-cond-fac", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(2, 3))
    hist(marg_post_sigma[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")
    hist(marg_post_sigma[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B")
    hist(marg_post_sigma[["C"]], freq = FALSE, main = "marg_post_x_fac2t = C")

    .plot_prior_density_for_test(marg_post_sigma[["A"]], main = "marg_post_x_fac2t = A")
    curve(dnorm(x, 0, sqrt(1^2 + 0.25^2)), add = TRUE)
    .plot_prior_density_for_test(marg_post_sigma[["B"]], main = "marg_post_x_fac2t = B")
    curve(dnorm(x, 0, sqrt(1^2 + 0.25^2)), add = TRUE)
    .plot_prior_density_for_test(marg_post_sigma[["C"]], main = "marg_post_x_fac2t = C")
    curve(dnorm(x, 0, sqrt(1^2 + 0.25^2)), add = TRUE)
  })

  ### marginal_inference ----
  out <- as_marginal_inference(
    model               = fit,
    parameters          = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    marginal_parameters = c("mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    conditional_list    = list(
      "mu_intercept"               = c(),
      "mu_x_cont1"                 = c("mu_x_cont1"),
      "mu_x_fac2t"                 = c("mu_x_cont1", "mu_x_fac2t"),
      "mu_x_fac3md"                = c("mu_x_fac3md"),
      "mu_x_cont1__xXx__x_fac3md"  = c("mu_x_fac2t", "mu_x_fac3md","mu_x_cont1__xXx__x_fac3md")
    ),
    conditional_rule    = "OR",
    formula      =  ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    silent       = TRUE
  )

  # test samples against previously generated ones
  vdiffr::expect_doppelganger("marginal_inference-ss-cont",     function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_cont1[["-1SD"]], freq = FALSE, main = "mu_x_cont1 = -1SD", breaks = 20)
    lines(density(out$averaged$mu_x_cont1[["-1SD"]]))

    hist(marg_post_x_cont1[["0SD"]], freq = FALSE, main = "mu_x_cont1 = 0SD", breaks = 20)
    lines(density(out$averaged$mu_x_cont1[["0SD"]]))

    hist(marg_post_x_cont1[["1SD"]], freq = FALSE, main = "mu_x_cont1 = +1SD", breaks = 20)
    lines(density(out$averaged$mu_x_cont1[["1SD"]]))

  })
  vdiffr::expect_doppelganger("marginal_inference-ss-cont-p",   function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_cont1[["-1SD"]], main = "mu_x_cont1 = -1SD", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_cont1[["-1SD"]], add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_cont1[["0SD"]], main = "mu_x_cont1 = 0SD", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_cont1[["0SD"]], add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_cont1[["1SD"]], main = "mu_x_cont1 = 1SD", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_cont1[["1SD"]], add = TRUE, lty = 2)
  })
  vdiffr::expect_doppelganger("marginal_inference-ss-fac.md",   function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(marg_post_x_fac3md[["A"]], freq = FALSE, main = "marg_post_x_fac3md = A", breaks = 20)
    lines(density(out$averaged$mu_x_fac3md$A))

    hist(marg_post_x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
    lines(density(out$averaged$mu_x_fac3md$B))

    hist(marg_post_x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac3md = C", breaks = 20)
    lines(density(out$averaged$mu_x_fac3md$C))

  })
  vdiffr::expect_doppelganger("marginal_inference-ss-fac.md-p", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    .plot_prior_density_for_test(marg_post_x_fac3md[["A"]], main = "marginal prior x_fac3md = A", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_fac3md$A, add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_fac3md[["B"]], main = "marginal prior x_fac3md = B", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_fac3md$B, add = TRUE, lty = 2)
    .plot_prior_density_for_test(marg_post_x_fac3md[["C"]], main = "marginal prior x_fac3md = C", xlim = c(-5, 5), ylim = c(0, 0.4))
    .plot_prior_density_for_test(out$averaged$mu_x_fac3md$C, add = TRUE, lty = 2)
  })
  # the previous BFs were based on model-averaged posteriors so they won't match

  # test summary table (note that these differ from the first set of tests because of the different model settings)
  marginal_parameters <- c(
    "mu_intercept",
    "mu_x_cont1",
    "mu_x_fac2t",
    "mu_x_fac3md",
    "mu_x_cont1__xXx__x_fac3md"
  )
  marginal_table <- marginal_estimates_table(
    out$conditional,
    out$inference,
    parameters = marginal_parameters
  )
  test_reference_table_stochastic(
    marginal_table,
    "marginal_estimates_table_spike_slab.txt",
    info_msg = "marginal_estimates_table for spike-and-slab"
  )
  .expect_marginal_table_current_inputs(
    table = marginal_table,
    samples = out$conditional,
    inference = out$inference,
    parameters = marginal_parameters
  )

  # plots
  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_fac2t-1", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t")})
  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_fac2t-2", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", par_name = "fac2t", lwd = 2)})
  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_fac2t-3", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2))})
  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_fac2t-4", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1))})
  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_fac2t-5", function(){plot_marginal(out$conditional, parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2), transformation = "exp", xlim = c(0, 5), transformation_settings = T)})

  vdiffr::expect_doppelganger("ggplot_marginal-ss-mu_x_fac2t-1", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t"))
  vdiffr::expect_doppelganger("ggplot_marginal-ss-mu_x_fac2t-2", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t", par_name = "fac2t", lwd = 2))
  vdiffr::expect_doppelganger("ggplot_marginal-ss-mu_x_fac2t-3", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2)))
  vdiffr::expect_doppelganger("ggplot_marginal-ss-mu_x_fac2t-4", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_fac2t", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1)))

  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_cont1", function(){plot_marginal(out$conditional, parameter = "mu_x_cont1", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1))})
  vdiffr::expect_doppelganger("ggplot_marginal-ss-mu_x_cont1", plot_marginal(out$conditional, plot_type = "ggplot", parameter = "mu_x_cont1", prior = TRUE, dots_prior = list(lty = 2), xlim = c(0, 1)))

  vdiffr::expect_doppelganger("plot_marginal-ss-mu_x_fac3md", function(){plot_marginal(out$averaged, parameter = "mu_x_fac3md", prior = TRUE, dots_prior = list(lty = 2), xlim = c(-1, 1))})
  vdiffr::expect_doppelganger("ggplot_marginal-ss-mu_x_fac3md", plot_marginal(out$averaged, plot_type = "ggplot", parameter = "mu_x_fac3md", prior = TRUE, dots_prior = list(lty = 2), xlim = c(-1, 1)))

  vdiffr::expect_doppelganger("plot_marginal-ss-int", plot_marginal(out$averaged, plot_type = "ggplot", parameter = "mu_intercept", prior = TRUE, dots_prior = list(lty = 2), xlim = c(-1, 1)))

})


test_that("Marginal distributions with one-sided weightfunction model work", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # Load pre-fitted one-sided weightfunction model
  fit_wf <- readRDS(file.path(temp_fits_dir, "fit_wf_onesided.RDS"))

  mixed_posteriors <- as_mixed_posteriors(
    model        = fit_wf,
    parameters   = "omega"
  )

  # Visual tests for weightfunction posteriors
  vdiffr::expect_doppelganger("marginal-wf-onesided-hist", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(mixed_posteriors$omega[,1], freq = FALSE, main = "omega[0,0.025]", breaks = 50, xlim = c(0, 1))
    hist(mixed_posteriors$omega[,2], freq = FALSE, main = "omega[0.025,1]", breaks = 50, xlim = c(0, 1))
  })

})


test_that("Marginal distributions with independent factor model work", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Load pre-fitted independent factor model
  fit_ind <- readRDS(file.path(temp_fits_dir, "fit_factor_independent.RDS"))

  mixed_posteriors <- as_mixed_posteriors(
    model        = fit_ind,
    parameters   = "p1"
  )
  marginal_posteriors <- marginal_posterior(mixed_posteriors, parameter = "p1", prior_samples = TRUE)

  # Visual tests for independent factor posteriors
  vdiffr::expect_doppelganger("marginal-factor-independent-hist", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(mixed_posteriors$p1[,1], freq = FALSE, main = "p1[1] (level 1)", breaks = 50)
    lines(density(marginal_posteriors[[1]]))
    .plot_prior_density_for_test(marginal_posteriors[[1]], add = TRUE, lty = 2)
    hist(mixed_posteriors$p1[,2], freq = FALSE, main = "p1[2] (level 2)", breaks = 50)
    lines(density(marginal_posteriors[[2]]))
    .plot_prior_density_for_test(marginal_posteriors[[2]], add = TRUE, lty = 2)
    hist(mixed_posteriors$p1[,3], freq = FALSE, main = "p1[3] (level 3)", breaks = 50)
    lines(density(marginal_posteriors[[3]]))
    .plot_prior_density_for_test(marginal_posteriors[[3]], add = TRUE, lty = 2)
  })

})


test_that("transform_scaled raw coefficients of fit_complex_mixed use their own prior", {

  skip_if_not_installed("rjags")

  # mu_x_cont1 ~ N(0, 1) x spike(0.5) with multiply_by = "sigma"; the monitored
  # coefficient is the raw node, so its prior is 1/2 at 0 + 1/2 N(0, 1)
  # (the fit has no formula scaling, so transform_scaled leaves it unchanged)
  fit <- readRDS(file.path(temp_fits_dir, "fit_complex_mixed.RDS"))
  expect_identical(attr(attr(fit, "prior_list")$mu_x_cont1, "multiply_by"), "sigma")

  samples <- as_mixed_posteriors(
    fit, parameters = c("mu_intercept", "mu_x_cont1", "sigma"), transform_scaled = TRUE
  )
  marginal <- marginal_posterior(samples, "mu_x_cont1", use_formula = FALSE, prior_samples = TRUE)
  prior_density <- .bt_meta_get(marginal, "prior_density")
  expect_equal(prior_density_ordinate(prior_density, 0)$point_mass, .5, tolerance = 1e-12)
  for(value in c(0, .5)){
    .expect_prior_height_for_test(marginal, value, .5 * stats::dnorm(value))
  }
})
