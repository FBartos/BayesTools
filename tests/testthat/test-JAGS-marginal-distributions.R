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

.plot_prior_density_for_test <- function(x, main = "", xlim = NULL, ylim = NULL, add = FALSE,
                                         lty = 1, col = graphics::par("fg"), ...){
  prior_density <- attr(x, "prior_density")
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

.marginal_posterior_with_prior_density_for_test <- function(samples, prior_density) {
  class(samples) <- c("marginal_posterior.simple", "marginal_posterior", class(samples))
  attr(samples, "prior_density") <- prior_density
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
  expect_true(posterior_density_method_uses_precomputed("qCMDE"))
  expect_true(posterior_density_method_uses_precomputed("IWMDE"))
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
  attr(posterior, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

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

  attr(posterior, "posterior_ordinate") <- list(
    parameter = "phi",
    value     = 0,
    ordinate  = .5,
    method    = "qCMDE"
  )
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

  attr(posterior, "posterior_ordinate") <- NULL
  attr(posterior, "posterior_density") <- list(
    parameter = "phi",
    x         = seq(-1, 1, length.out = 101),
    y         = rep(.5, 101),
    method    = "qCMDE"
  )
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
  attr(posterior, "posterior_ordinate") <- ordinate

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
    y              = stats::dnorm(seq(-2, 2, length.out = 201)),
    method         = "iwmde",
    density_method = "IWMDE",
    point_masses   = data.frame(x = 0, mass = .01),
    support        = list(bounds = c(-Inf, Inf), points = 0),
    parameter      = "theta"
  )
  parsed_density <- BayesTools:::.posterior_density_from_attribute(density)
  expect_equal(parsed_density[["method"]], "iwmde")
  expect_equal(parsed_density[["point_masses"]][["mass"]], .01)
  expect_true(BayesTools:::.posterior_density_point_masses_declared(parsed_density))
  expect_equal(parsed_density[["support"]][["bounds"]], c(-Inf, Inf))
  expect_equal(parsed_density[["support"]][["points"]], 0)
  density_without_points <- posterior_density_attribute(
    x              = seq(-2, 2, length.out = 201),
    y              = stats::dnorm(seq(-2, 2, length.out = 201)),
    method         = "iwmde",
    density_method = "IWMDE"
  )
  expect_false(BayesTools:::.posterior_density_point_masses_declared(
    BayesTools:::.posterior_density_from_attribute(density_without_points)
  ))
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
      point_masses   = data.frame(x = 0, mass = 1.2)
    ),
    "point_masses"
  )
  expect_error(
    posterior_density_attribute(
      x              = 0:1,
      y              = c(1, 1),
      method         = "iwmde",
      density_method = "IWMDE",
      point_masses   = data.frame(x = c(0, 1), mass = c(.6, .5))
    ),
    "point_masses"
  )
  expect_error(
    posterior_density_attribute(
      x              = 0:1,
      y              = c(1, 1),
      method         = "iwmde",
      density_method = "IWMDE",
      density        = list(x = 0:1, y = c(1, 1))
    ),
    "reserved fields"
  )
  expect_error(
    posterior_density_attribute(
      x              = 0:1,
      y              = c(1, 1),
      method         = "iwmde",
      density_method = "IWMDE",
      support        = c(1, 0)
    ),
    "support metadata"
  )
  expect_warning(
    expect_error(
      posterior_density_attribute(
        x              = 0:1,
        y              = c(1, 1),
        method         = "iwmde",
        density_method = "IWMDE",
        support        = list(bounds = "bad")
      ),
      "support metadata"
    ),
    NA
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

  attr(posterior, "posterior_density") <- list(
    x      = seq(2, 3, length.out = 101),
    y      = rep(1, 101),
    method = "qCMDE"
  )
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
  attr(posterior, "posterior_density") <- list(
    x           = stored_x,
    y           = stored_y,
    method      = "iwmde",
    diagnostics = list(bf_relative_mcse = .123)
  )

  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_null(attr(out, "BF_error_percent"))

  attr(posterior, "posterior_density")$diagnostics <- list(
    bf_value         = 0,
    bf_relative_mcse = .123
  )
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
  attr(posterior, "posterior_density") <- list(
    x      = seq(2, 3, length.out = 101),
    y      = rep(100, 101),
    method = "plot-only"
  )
  attr(posterior, "posterior_ordinate") <- list(
    value       = c(0, .5),
    ordinate    = c(.25, .5),
    method      = "qCMDE",
    diagnostics = list(relative_mcse = c(.1, .2))
  )

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

test_that("stored posterior ordinate parser rejects ambiguous values", {

  parsed <- BayesTools:::.posterior_ordinate_from_attribute(
    list(
      ordinates = data.frame(
        value         = c(0, .5),
        ordinate      = c(.25, .5),
        relative_mcse = c(.1, .2)
      ),
      method = "qCMDE"
    ),
    null_hypothesis = .5
  )
  expect_equal(parsed$x, .5)
  expect_equal(parsed$y, .5)
  expect_equal(parsed$diagnostics$relative_mcse, .2)

  parsed <- BayesTools:::.posterior_ordinate_from_attribute(
    list(
      ordinates = data.frame(
        value            = c(0, .5),
        ordinate         = c(.25, .5),
        relative_mcse    = c(.1, .2),
        BF_error_percent = c(10, 20)
      ),
      diagnostics = list(relative_mcse = .9, estimator = "q_grid_cmde"),
      method      = "qCMDE"
    ),
    null_hypothesis = .5
  )
  expect_equal(parsed$diagnostics$relative_mcse, .2)
  expect_equal(parsed$diagnostics$BF_error_percent, 20)
  expect_equal(parsed$diagnostics$estimator, "q_grid_cmde")

  expect_null(BayesTools:::.posterior_ordinate_from_attribute(
    list(value = c(0, 0), ordinate = c(.25, .30)),
    null_hypothesis = 0
  ))
  expect_null(BayesTools:::.posterior_ordinate_from_attribute(
    list(value = 0, ordinate = 0),
    null_hypothesis = 0
  ))
  expect_null(BayesTools:::.posterior_ordinate_from_attribute(
    list(value = 1, ordinate = .25),
    null_hypothesis = 0
  ))
})

test_that("stored posterior ordinate parser keeps diagnostics aligned", {

  parsed <- BayesTools:::.posterior_ordinate_from_attribute(
    list(
      value       = c(NA, 0),
      ordinate    = c(.25, .50),
      method      = "qCMDE",
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
      list(parameter = "theta", value = 0,  ordinate = .25, method = "qCMDE"),
      list(parameter = "theta", value = .5, ordinate = .50, method = "qCMDE")
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
  attr(posterior, "posterior_density") <- list(
    x      = seq(-4, 4, length.out = 401),
    y      = rep(100, 401),
    method = "iwmde"
  )

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
  attr(posterior, "posterior_support") <- c(0, 1)

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

  attr(posterior, "posterior_support") <- list(lower = 0, upper = 1, exact = TRUE)
  out <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE
  )
  expect_true(attr(out, "posterior_density_boundary_reflection"))
  expect_equal(attr(out, "posterior_density_support"), c(0, 1))
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
  attr(posterior, "posterior_density") <- list(
    x       = seq(.25, .75, length.out = 101),
    y       = rep(1, 101),
    method  = "iwmde",
    support = list(bounds = c(0, 1), exact = TRUE)
  )

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
  attr(posterior, "posterior_support") <- c(0, 1)
  attr(posterior, "posterior_ordinate") <- list(
    value    = -.5,
    ordinate = .5,
    method   = "qCMDE"
  )

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

  attr(posterior, "posterior_ordinate") <- NULL
  attr(posterior, "posterior_density") <- list(
    x       = seq(-1, 1, length.out = 101),
    y       = rep(.5, 101),
    method  = "qCMDE",
    support = list(lower = 0, upper = 1, exact = TRUE)
  )
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
  attr(posterior, "posterior_ordinate") <- list(
    value    = -.5,
    ordinate = .5,
    method   = "stale-ordinate"
  )
  attr(posterior, "posterior_density") <- list(
    x       = seq(0, 1, length.out = 101),
    y       = rep(.5, 101),
    method  = "support-guard",
    support = list(lower = 0, upper = 1, exact = TRUE)
  )

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

  attr(posterior, "posterior_support") <-
    BayesTools:::.posterior_support_new(c(0, 1), exact = FALSE, source = "test")
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

  attr(posterior, "posterior_support") <-
    BayesTools:::.posterior_support_new(c(.2, .8), source = "test")
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
  attr(posterior, "posterior_support") <-
    BayesTools:::.posterior_support_new(c(.2, .8), source = "test")
  attr(posterior, "posterior_ordinate") <- posterior_ordinate_attribute(
    value          = .5,
    ordinate       = 1,
    method         = "qCMDE",
    density_method = "precomputed"
  )

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
  attr(posterior, "posterior_support") <-
    BayesTools:::.posterior_support_new(c(0, 1), source = "test")
  attr(posterior, "posterior_ordinate") <- posterior_ordinate_attribute(
    value          = -.5,
    ordinate       = .5,
    method         = "qCMDE",
    density_method = "precomputed"
  )

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
  attr(theta, "sample_ind") <- seq_along(theta)
  attr(theta, "models_ind") <- rep(1, length(theta))
  attr(theta, "parameter") <- "theta"
  attr(theta, "prior_list") <- theta_prior
  theta <- BayesTools:::.posterior_support_set_from_prior_list(theta, theta_prior)

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")
  attr(samples, "prior_density_context") <- BayesTools:::.prior_density_build_context(
    prior_list   = list(theta = theta_prior),
    column_names = "theta",
    n_grid       = 1024
  )

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

  attr(theta, "posterior_support") <-
    BayesTools:::.posterior_support_new(c(.2, .8), source = "posterior")
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

  connected_support <- BayesTools:::.posterior_support_union(list(c(0, 1), c(1, 2)))
  disjoint_support <- BayesTools:::.posterior_support_union(list(c(0, 1), c(2, 3)))
  expect_true(connected_support$exact)
  expect_false(disjoint_support$exact)
})

test_that("marginal_posterior infers support from the current prior context", {

  raw_prior <- prior("beta", list(alpha = 1, beta = 1))
  transformed_prior <- prior("uniform", list(10, 20))
  theta <- seq(11, 19, length.out = 51)
  class(theta) <- c("mixed_posteriors", "mixed_posteriors.simple", class(theta))
  attr(theta, "sample_ind") <- seq_along(theta)
  attr(theta, "models_ind") <- rep(1, length(theta))
  attr(theta, "parameter") <- "theta"
  attr(theta, "prior_list") <- raw_prior

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")
  attr(samples, "transform_scaled") <- TRUE
  attr(samples, "prior_density_context") <- BayesTools:::.prior_density_context(
    prior_list   = list(theta = transformed_prior),
    column_names = "theta",
    n_grid       = 64
  )

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
  attr(theta, "sample_ind") <- seq_along(theta)
  attr(theta, "models_ind") <- rep(1, length(theta))
  attr(theta, "parameter") <- "theta"
  attr(theta, "prior_list") <- theta_prior
  theta <- BayesTools:::.posterior_support_set_from_prior_list(theta, theta_prior)
  theta <- BayesTools:::.condition_event_set_attributes(theta, condition_event)

  samples <- list(theta = theta)
  class(samples) <- c("mixed_posteriors", "list")
  attr(samples, "prior_density_context") <- BayesTools:::.prior_density_context(
    prior_list   = prior_list,
    column_names = "theta",
    n_grid       = 64
  )

  marginal <- marginal_posterior(
    samples,
    parameter     = "theta",
    prior_samples = FALSE
  )

  expect_equal(BayesTools:::.posterior_support_bounds(marginal), c(10, 20))
  expect_equal(attr(marginal, "conditional", exact = TRUE), "theta")
  expect_equal(attr(marginal, "condition_key", exact = TRUE), condition_event[["condition_key"]])
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
  attr(samples, "posterior_support") <- list(
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
  )
  context <- BayesTools:::.prior_density_context(
    prior_list   = list(theta = slab_prior),
    column_names = "theta"
  )

  refreshed <- BayesTools:::.posterior_support_set_from_prior_context(
    samples,
    context
  )
  refreshed_support <- attr(refreshed, "posterior_support", exact = TRUE)

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
  attr(theta, "sample_ind") <- seq_along(theta)
  attr(theta, "models_ind") <- rep(1, length(theta))
  attr(theta, "parameter") <- "mu_x"
  attr(theta, "formula_parameter") <- "mu"
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
  attr(mu_intercept, "sample_ind") <- seq_along(mu_intercept)
  attr(mu_intercept, "models_ind") <- rep(1, length(mu_intercept))
  attr(mu_intercept, "parameter") <- "mu_intercept"
  attr(mu_intercept, "formula_parameter") <- "mu"
  attr(mu_intercept, "prior_list") <- prior("normal", list(0, 1))

  mu_x <- seq(-1, 1, length.out = 51)
  class(mu_x) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "mixed_posteriors.formula",
    class(mu_x)
  )
  attr(mu_x, "sample_ind") <- seq_along(mu_x)
  attr(mu_x, "models_ind") <- rep(1, length(mu_x))
  attr(mu_x, "parameter") <- "mu_x"
  attr(mu_x, "formula_parameter") <- "mu"
  attr(mu_x, "prior_list") <- prior("normal", list(0, 1))

  samples <- list(mu_intercept = mu_intercept, mu_x = mu_x)
  class(samples) <- c("mixed_posteriors", "list")
  attr(samples, "posterior_density") <- list(
    one_sd = posterior_density_attribute(
      x         = seq(-1, 1, length.out = 101),
      y         = rep(.5, 101),
      method    = "formula-density",
      density_method = "precomputed",
      parameter = "mu_x[1SD]"
    )
  )

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
    attr(marginal[["1SD"]], "posterior_density", exact = TRUE)[["method"]],
    "formula-density"
  )
  expect_null(attr(marginal[["0SD"]], "posterior_density", exact = TRUE))
  expect_null(attr(transformed[["1SD"]], "posterior_density", exact = TRUE))
})

test_that("spike-and-slab and vector constructors attach support metadata", {

  slab <- prior(
    "beta",
    list(alpha = 2, beta = 2),
    truncation = list(lower = .2, upper = .8)
  )
  spike_slab <- prior_spike_and_slab(
    slab,
    prior_inclusion = prior("point", list(location = .5))
  )

  prior_samples <- BayesTools:::.as_mixed_priors.spike_and_slab(
    spike_slab,
    parameter = "theta",
    n_samples = 200
  )
  posterior_samples <- BayesTools:::.as_mixed_posteriors.spike_and_slab(
    cbind(theta = seq(.2, .8, length.out = 100),
          theta_indicator = rep(c(0, 1), 50)),
    spike_slab,
    parameter = "theta"
  )
  prior_support <- BayesTools:::.posterior_support_get(prior_samples)
  posterior_support <- BayesTools:::.posterior_support_get(posterior_samples)

  expect_equal(prior_support$bounds, c(0, .8))
  expect_equal(posterior_support$bounds, c(0, .8))
  expect_true(0 %in% prior_support$points)
  expect_true(0 %in% posterior_support$points)
  expect_false(prior_support$exact)

  vector_prior <- prior("mnormal", list(mean = 0, sd = 1, K = 2))
  vector_samples <- BayesTools:::.as_mixed_priors.vector(
    vector_prior,
    parameter = "theta",
    n_samples = 25
  )
  vector_support <- attr(vector_samples, "posterior_support", exact = TRUE)

  expect_equal(names(vector_support), colnames(vector_samples))
  expect_true(all(vapply(
    vector_support,
    function(support) identical(support$bounds, c(-Inf, Inf)),
    logical(1)
  )))
})

test_that("simplex posterior support uses component and convex-hull bounds", {

  simplex_prior <- prior("dirichlet", list(alpha = c(2, 3, 5)))

  component_support <- BayesTools:::.posterior_support_from_prior(simplex_prior)
  expect_equal(component_support$bounds, c(0, 1))
  expect_true(component_support$exact)

  simplex_samples <- BayesTools:::.as_mixed_priors.vector(
    simplex_prior,
    parameter = "w",
    n_samples = 25
  )
  simplex_support <- attr(simplex_samples, "posterior_support", exact = TRUE)
  expect_true(all(vapply(
    simplex_support,
    function(support) identical(support$bounds, c(0, 1)),
    logical(1)
  )))

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
  attr(child, "posterior_density") <- posterior_density_attribute(
    x              = seq(-1, 1, length.out = 11),
    y              = rep(.5, 11),
    method         = "child-density",
    density_method = "precomputed",
    parameter      = "theta[A]"
  )

  marginal <- list(A = child)
  class(marginal) <- c("marginal_posterior.factor", "list")
  samples <- structure(
    list(),
    class = c("mixed_posteriors", "list"),
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
    attr(out[["A"]], "posterior_density", exact = TRUE)[["method"]],
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
  attr(posterior, "posterior_density") <- list(
    x      = seq(.5, 1, length.out = 101),
    y      = rep(1, 101),
    method = "iwmde"
  )

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

  attr(posterior, "posterior_density") <- list(
    x      = 0,
    y      = 1,
    method = "invalid-density"
  )
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis = 0,
      density_method  = "precomputed"
    ),
    "Precomputed posterior density metadata is present but invalid",
    fixed = TRUE
  )

  attr(posterior, "posterior_density") <- list(
    x            = seq(-1, 1, length.out = 101),
    y            = rep(1, 101),
    method       = "invalid-point-mass",
    point_masses = list(location = 0, p = 1.2)
  )
  expect_error(
    Savage_Dickey_BF(
      posterior,
      null_hypothesis = 0,
      density_method  = "precomputed"
    ),
    "Precomputed posterior density metadata is present but invalid",
    fixed = TRUE
  )

  attr(posterior, "posterior_density") <- NULL
  attr(posterior, "posterior_ordinate") <- list(
    value    = 1,
    ordinate = .5,
    method   = "wrong-null"
  )
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
  attr(posterior, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

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

  parsed <- BayesTools:::.posterior_density_from_attribute(list(
    x      = c(0, 0, 1),
    y      = c(1, 3, 2),
    method = "iwmde"
  ))
  expect_equal(parsed$x, c(0, 1))
  expect_equal(unname(parsed$y), c(2, 2))

  expect_null(BayesTools:::.posterior_density_from_attribute(list(
    x = c(0, 1),
    y = c(0, 0)
  )))
  expect_null(BayesTools:::.posterior_density_from_attribute(list(
    x = c(1, 1),
    y = c(1, 2)
  )))
})

test_that("stored posterior density selection validates names and conditionals", {

  density <- list(
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

  alias_branch <- list(
    parameter = "phi",
    value     = 0,
    ordinate  = 100,
    method    = "wrong-parameter"
  )
  valid_branch <- list(
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

test_that("stored posterior density parser validates and aggregates point masses", {

  parsed <- BayesTools:::.posterior_density_from_attribute(list(
    x = seq(-1, 1, length.out = 11),
    y = rep(1, 11),
    point_masses = list(
      location = c(0, 0, 1),
      p        = c(.2, .3, .4)
    )
  ))

  expect_equal(parsed$point_masses$x, c(0, 1))
  expect_equal(parsed$point_masses$mass, c(.5, .4))

  parsed_invalid <- BayesTools:::.posterior_density_from_attribute(list(
    x = seq(-1, 1, length.out = 11),
    y = rep(1, 11),
    point_masses = list(
      x    = c(0, 1, 2),
      mass = c(.6, .5, 2)
    )
  ))
  expect_null(parsed_invalid)
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
  attr(posterior, "posterior_density") <- list(
    x      = seq(-1, 1, length.out = 101),
    y      = rep(1, 101),
    method = "iwmde"
  )

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
  attr(posterior_list, "posterior_density") <- list(
    level = list(
      x      = stored_x,
      y      = stored_y,
      method = "iwmde"
    )
  )

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

  attr(posterior_list, "posterior_density") <- NULL
  attr(posterior_list, "posterior_densities") <- list(
    list(
      level = list(
        x      = stored_x,
        y      = stored_y,
        method = "iwmde"
      )
    )
  )
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
  attr(posterior, "conditional") <- "theta"
  attr(posterior, "conditional_rule") <- "AND"
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  attr(posterior_list, "posterior_density") <- list(
    level = list(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(100, 101),
      method = "stale"
    ),
    matching = list(
      parameter   = "level",
      conditional = "theta",
      x           = seq(-1, 1, length.out = 101),
      y           = rep(.50, 101),
      method      = "qCMDE"
    )
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50
  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)

  attr(posterior_list, "posterior_density") <- NULL
  attr(posterior_list, "posterior_ordinates") <- list(list(
    level = list(
      value    = 0,
      ordinate = 100,
      method   = "stale"
    ),
    matching = list(
      parameter   = "level",
      conditional = "theta",
      value       = 0,
      ordinate    = .50,
      method      = "qCMDE"
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
  attr(posterior, "conditional") <- "theta"
  attr(posterior, "conditional_rule") <- "AND"
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")

  mismatched_density <- list(
    parameter   = "level",
    conditional = "phi",
    x           = seq(-1, 1, length.out = 101),
    y           = rep(100, 101),
    method      = "stale-density"
  )
  attr(posterior_list, "posterior_density") <- list(mismatched_density)

  expect_null(BayesTools:::.posterior_density_child_attributes(posterior_list)[[1]])

  matching_density <- mismatched_density
  matching_density[["conditional"]] <- "theta"
  attr(posterior_list, "posterior_density") <- list(matching_density)
  expect_equal(
    BayesTools:::.posterior_density_child_attributes(posterior_list)[[1]][["method"]],
    "stale-density"
  )

  mismatched_ordinate <- list(
    parameter   = "level",
    conditional = "phi",
    value       = 0,
    ordinate    = 100,
    method      = "stale-ordinate"
  )
  attr(posterior_list, "posterior_density") <- NULL
  attr(posterior_list, "posterior_ordinate") <- list(mismatched_ordinate)

  expect_null(BayesTools:::.posterior_ordinate_child_attributes(
    posterior_list,
    null_hypothesis = 0
  )[[1]])

  matching_ordinate <- mismatched_ordinate
  matching_ordinate[["conditional"]] <- "theta"
  attr(posterior_list, "posterior_ordinate") <- list(matching_ordinate)
  expect_equal(
    BayesTools:::.posterior_ordinate_child_attributes(
      posterior_list,
      null_hypothesis = 0
    )[[1]][["method"]],
    "stale-ordinate"
  )
})

test_that("Savage_Dickey_BF accepts list-of-record posterior ordinates", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  attr(posterior, "posterior_ordinate") <- list(
    parameter   = "theta",
    method      = "qCMDE",
    diagnostics = list(relative_mcse = c(.10, .20)),
    ordinates   = list(
      list(value = -.50, ordinate = 100),
      list(value = 0, ordinate = .50)
    )
  )

  expect_true(BayesTools:::.posterior_ordinate_has_data(
    attr(posterior, "posterior_ordinate")
  ))

  matched_source <- BayesTools:::.posterior_ordinate_from_sources(
    sources         = list(attr(posterior, "posterior_ordinate")),
    aliases         = "theta",
    null_hypothesis = 0
  )
  expect_equal(matched_source[["method"]], "qCMDE")

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
  attr(posterior, "posterior_density") <- list(
    x      = seq(-1, 1, length.out = 101),
    y      = rep(.50, 101),
    method = "child-density"
  )
  attr(posterior, "posterior_ordinate") <- list(
    value    = 0,
    ordinate = .50,
    method   = "child-ordinate"
  )
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  attr(posterior_list, "posterior_density") <- list(
    level = list(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(100, 101),
      method = "top-density"
    )
  )
  attr(posterior_list, "posterior_ordinate") <- list(
    level = list(
      value    = 0,
      ordinate = 100,
      method   = "top-ordinate"
    )
  )

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

test_that("Savage_Dickey_BF replaces invalid child sources with valid top-level sources", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  attr(posterior, "posterior_ordinate") <- list(
    value    = 1,
    ordinate = 100,
    method   = "wrong-null"
  )
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  attr(posterior_list, "posterior_ordinate") <- list(
    level = list(
      value    = 0,
      ordinate = .50,
      method   = "top-ordinate"
    )
  )

  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / .50
  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)

  attr(posterior, "posterior_ordinate") <- NULL
  attr(posterior, "posterior_density") <- list(
    x      = 0,
    y      = 100,
    method = "invalid-child-density"
  )
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  attr(posterior_list, "posterior_density") <- list(
    level = list(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(.50, 101),
      method = "top-density"
    )
  )

  out <- Savage_Dickey_BF(
    posterior_list,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE,
    density_method       = "precomputed"
  )

  expect_equal(as.numeric(out[["level"]]), expected, tolerance = 1e-12)
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
  attr(posterior, "posterior_density") <- list(
    x      = seq(1, 2, length.out = 101),
    y      = rep(100, 101),
    method = "child-misses-null"
  )
  posterior_list <- list(level = posterior)
  class(posterior_list) <- c("marginal_posterior", "list")
  attr(posterior_list, "posterior_density") <- list(
    level = list(
      x      = seq(-1, 1, length.out = 101),
      y      = rep(.50, 101),
      method = "top-density"
    )
  )

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

test_that("Savage_Dickey_BF errors for point mass and posterior-null clusters", {

  continuous_prior <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior_cluster <- .marginal_posterior_with_prior_density_for_test(
    c(rep(0, 8), seq(-2, 2, length.out = 92)),
    continuous_prior
  )

  expect_error(
    Savage_Dickey_BF(posterior_cluster, null_hypothesis = 0, normal_approximation = TRUE),
    "posterior samples at the exact null hypothesis"
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
  attr(posterior_with_stored_point, "posterior_density") <- list(
    x            = seq(-2, 2, length.out = 101),
    y            = stats::dnorm(seq(-2, 2, length.out = 101)),
    point_masses = list(location = 0, p = .2)
  )

  expect_error(
    Savage_Dickey_BF(
      posterior_with_stored_point,
      null_hypothesis = 0,
      density_method  = "precomputed"
    ),
    "Stored posterior density contains"
  )
})

test_that("Savage_Dickey_BF diagnoses zero prior density at point null", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 2, beta = 2))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior <- .marginal_posterior_with_prior_density_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  attr(posterior, "posterior_support") <-
    BayesTools:::.posterior_support_new(c(0, 1), source = "test")

  out <- Savage_Dickey_BF(posterior, null_hypothesis = 0, silent = TRUE)
  expect_equal(as.numeric(out), 0)
  expect_match(
    attr(out, "warnings"),
    "Prior density at the null hypothesis value is zero or non-finite",
    fixed = TRUE
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
  attr(posterior, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

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

test_that("plot_marginal does not add sample spikes to stored full density", {

  posterior <- .marginal_posterior_with_prior_density_for_test(
    c(rep(0, 25), seq(-2, 2, length.out = 75)),
    BayesTools:::.prior_linear_density_point(0, p = .25)
  )
  stored_x <- seq(-2, 2, length.out = 51)
  stored_y <- stats::dnorm(stored_x)
  attr(posterior, "posterior_density") <- list(
    x      = stored_x,
    y      = stored_y,
    method = "iwmde"
  )

  expect_warning(
    plot_data <- BayesTools:::.plot_data_marginal_samples(
      samples                  = list(theta = posterior),
      parameter                = "theta",
      prior                    = FALSE,
      n_points                 = 16,
      transformation           = NULL,
      transformation_arguments = NULL,
      transformation_settings  = FALSE,
      density_method           = "precomputed"
    ),
    "does not declare 'point_masses'",
    fixed = TRUE
  )

  expect_equal(plot_data$density1$x, stored_x)
  expect_equal(
    length(plot_data[vapply(plot_data, inherits, logical(1), "density.prior.point")]),
    0L
  )
})

test_that("plot_marginal accepts posterior density diagnostics", {

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
  attr(posterior, "posterior_density") <- list(
    parameter    = "theta",
    status       = "ok",
    point_masses = list(location = 0, p = .2),
    x            = stored_x,
    y            = stored_y,
    method       = "iwmde",
    diagnostics = list(min_ess = 40)
  )

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
    inherits(attr(x, "prior_density"), "prior_linear_density")
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
  expect_true(inherits(attr(marginal[["A"]], "prior_density"), "prior_linear_density"))
  expect_true(inherits(attr(marginal[["B"]], "prior_density"), "prior_linear_density"))
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(marginal[["A"]], "prior_density"), 0),
    1
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(marginal[["B"]], "prior_density"), 0),
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
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list

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
    inherits(attr(x, "prior_density"), "prior_linear_density")
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
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list

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
    expect_s3_class(attr(marginal[[i]], "prior_density"), "prior_linear_density")
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
      attr(marginal_exp[[level]], "prior_density"),
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
      attr(conditional_marginal[[level_name]], "prior_density"),
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

# File-level skips: All remaining tests in this file require pre-fitted models
skip_if_not_visual_fixture_tests()
skip_if_no_fits()
skip_if_not_installed("rjags")
skip_if_not_installed("bridgesampling")

test_that("Marginal distribution prior and posterior functions work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes

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
  # (uses model-averaged posteriors rather than conditional ones -- which would be correct)
  # test input
  expect_error(Savage_Dickey_BF(list(posterior_manual0)), "'BF_savage_dickey' function requires an object of class 'marginal_posteriors'")
  expect_error(Savage_Dickey_BF(marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "sigma",
    prior_samples     = FALSE)), "there are no prior densities for the posterior distribution")

  # simple restricted prior
  BF.marg_post_sigma <- Savage_Dickey_BF(marg_post_sigma, silent = TRUE)
  expect_true(is.finite(BF.marg_post_sigma))
  expect_gt(BF.marg_post_sigma, 1e10)
  expect_null(attr(BF.marg_post_sigma, "warnings"))
  expect_true(isTRUE(attr(BF.marg_post_sigma, "posterior_density_boundary_reflection")))
  expect_equal(attr(BF.marg_post_sigma, "posterior_density_support"), c(0, 5))

  # simple factor
  expect_error(
    Savage_Dickey_BF(marg_post_simple_x_fac2t),
    "exact null hypothesis value|point mass in the prior"
  )
  marg_post_simple_x_fac2t_B <- marg_post_simple_x_fac2t[["B"]]
  class(marg_post_simple_x_fac2t_B) <- c(class(marg_post_simple_x_fac2t_B), "marginal_posterior")
  expect_error(
    Savage_Dickey_BF(marg_post_simple_x_fac2t_B),
    "exact null hypothesis value",
    fixed = TRUE
  )


  BF.marg_post_x_fac3md <- Savage_Dickey_BF(marg_post_x_fac3md, silent = TRUE)
  expect_equal(BF.marg_post_x_fac3md, list("A" = Inf, "B" = Inf, "C" = Inf), ignore_attr = TRUE)
  expect_equal(attr(BF.marg_post_x_fac3md[["A"]], "warnings"),
               "Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.")

  BF2.marg_post_x_fac3md <- suppressWarnings(Savage_Dickey_BF(marg_post_x_fac3md, null_hypothesis = 0.5))
  expect_equal(BF2.marg_post_x_fac3md, list("A" = 4.062, "B" = 0.1404, "C" = 0.153), tolerance = 5e-3, ignore_attr = TRUE)

  BF2.marg_post_x_fac3md <- suppressWarnings(Savage_Dickey_BF(marg_post_x_fac3md, null_hypothesis = 0.5, normal_approximation = TRUE))
  expect_equal(BF2.marg_post_x_fac3md, list("A" = 0.614, "B" = 0.0997, "C" = 0.1237), tolerance = 5e-3, ignore_attr = TRUE)


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
  test_reference_table(
    marginal_estimates_table(out$conditional, out$inference, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md")),
    "marginal_estimates_table_model_avg.txt",
    info_msg = "marginal_estimates_table for model averaging"
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

test_that("Marginal distribution prior functions work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  set.seed(1)

  ### independent prior distribution ----
  priors <- list(
      prior_factor("spike",  list(0), contrast = "independent"),
      prior_factor("normal", list(0, .3), contrast = "independent"),
      prior_factor("normal", list(2, .3), contrast = "independent")
  )
  attr(priors[[1]], "levels") <- 3
  attr(priors[[2]], "levels") <- 3
  attr(priors[[3]], "levels") <- 3
  temp_prior <- BayesTools:::.mix_priors.factor(priors, "mu", seed = NULL, n_samples = 10000)


  vdiffr::expect_doppelganger("marginal-prior-ind", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 3))
    hist(temp_prior[,1], freq = FALSE, main = "marginal prior independent (1)", breaks = 50)
    hist(temp_prior[,2], freq = FALSE, main = "marginal prior independent (2)", breaks = 50)
    hist(temp_prior[,3], freq = FALSE, main = "marginal prior independent (3)", breaks = 50)

  })

  ### 3 level treatment prior distribution ----
  priors <- list(
    prior_factor("spike",  list(0),     contrast = "treatment"),
    prior_factor("normal", list(2, .3), contrast = "treatment")
  )
  attr(priors[[1]], "levels") <- 3
  attr(priors[[2]], "levels") <- 3
  temp_prior <- BayesTools:::.mix_priors.factor(priors, "mu", seed = NULL, n_samples = 10000)


  vdiffr::expect_doppelganger("marginal-prior-trt", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(temp_prior[,1], freq = FALSE, main = "marginal prior treatment (1)", breaks = 50)
    hist(temp_prior[,2], freq = FALSE, main = "marginal prior treatment (2)", breaks = 50)

  })

  ### weightfunction prior distribution ----
  priors <- list(
    prior_weightfunction("one-sided", c(0.05, 0.50), wf_fixed(c(1, 1, 1))),
    prior_weightfunction("one-sided", c(0.10), wf_cumulative(c(1, 1)))
  )
  temp_prior <- BayesTools:::.mix_priors.weightfunction(priors, "mu", seed = NULL, n_samples = 10000)

  vdiffr::expect_doppelganger("marginal-prior-weightfunction", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 4))
    hist(temp_prior[,1], freq = FALSE, main = "marginal prior weightfunction (1)", breaks = 50)
    hist(temp_prior[,2], freq = FALSE, main = "marginal prior weightfunction (2)", breaks = 50)
    hist(temp_prior[,3], freq = FALSE, main = "marginal prior weightfunction (3)", breaks = 50)
    hist(temp_prior[,4], freq = FALSE, main = "marginal prior weightfunction (4)", breaks = 50)

  })

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
  test_reference_table_numeric(
    marginal_estimates_table(out$conditional, out$inference, parameters = c("mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md")),
    "marginal_estimates_table_spike_slab.txt",
    tolerance = 1e-2,
    info_msg = "marginal_estimates_table for spike-and-slab"
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

  # Not implemented for weightfunctions
  #  marginal_posterior(mixed_posteriors, parameter = "omega", prior_samples = TRUE)
  temp_samples <- .as_mixed_priors.weightfunction(attr(fit_wf, "prior_list")[[1]], parameter = "omega")

  # Visual tests for weightfunction posteriors
  vdiffr::expect_doppelganger("marginal-wf-onesided-hist", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))

    par(mfrow = c(1, 2))
    hist(mixed_posteriors$omega[,1], freq = FALSE, main = "omega[0,0.025]", breaks = 50, xlim = c(0, 1))
    # The first one-sided weight is the fixed reference bin; do not smooth a point mass.
    if(stats::sd(temp_samples[,1]) > sqrt(.Machine$double.eps)){
      lines(density(temp_samples[,1]))
    }
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

