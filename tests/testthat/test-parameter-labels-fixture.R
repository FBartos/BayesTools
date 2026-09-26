skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Parameter labels of fitted models
# ============================================================================ #
#
# PURPOSE:
#   Labels, catalog quantities, and original-scale transforms of the rows and
#   mixed-posterior columns of real fits: log intercepts, LKJ primitives,
#   ordered-factor shares, and random-effect SDs of standardized predictors.
#
# MODELS/FIXTURES:
#   - fit_label_* models from test-00-model-fits.R
#
# TAGS: @fixture, @JAGS, @labels
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

.label_cached_fit <- function(name){

  skip_if_not_installed("rjags")
  skip_if_missing_fits(name)
  readRDS(file.path(temp_fits_dir, paste0(name, ".RDS")))
}

test_that("estimates tables render the rows of a log-intercept formula from their parts", {

  # the rows of a log-intercept (scale) formula on the original scale are
  # rendered from the table's parts, also under another formula parameter
  fit <- .label_cached_fit("fit_label_log_intercept")
  table <- JAGS_estimates_table(fit, transform_scaled = TRUE)
  quantities <- attr(table, "quantities")
  expect_identical(rownames(table), c("(mu) exp(intercept)", "(mu) x"))
  expect_identical(parameter_labels(quantities, style = "table"), rownames(table))
  expect_identical(quantities$label_parts[[1L]]$transformation, "exp")
  # the exponentiated log intercept names the intercept (a catalog alias)
  catalog <- parameter_catalog(fit)
  expect_identical(
    quantities$quantity_id,
    c(parameter_catalog_resolve(catalog, "(mu) exp(intercept)")$quantity_id,
      parameter_catalog_resolve(catalog, "mu_x")$quantity_id)
  )
  expect_identical(
    parameter_labels(
      .bt_label_parts_update(quantities$label_parts, formula_parameter = "log(tau)"),
      style = "table"
    ),
    c("(log(tau)) exp(intercept)", "(log(tau)) x")
  )
  # the exponentiated log intercept keeps the output transformations of its
  # values (e.g. of posterior_transform())
  lin_part <- .bt_label_parts_update(quantities$label_parts[1L], transformation = c("none", "lin"))
  expect_identical(
    .bt_label_parts_log_intercept(lin_part, attr(fit, "formula_scale"))[[1L]]$transformation,
    c("exp", "lin")
  )
})

test_that("raw rows of LKJ primitives are rendered backend coordinates", {

  fit <- .label_cached_fit("fit_label_lkj")
  raw <- JAGS_estimates_table(fit, random_effects_summary = "raw",
                              remove_diagnostics = TRUE)
  primitives <- c(
    "(mu) lkj_u(intercept,x | g)",
    "(mu) lkj_u(intercept,z | g)",
    "(mu) lkj_u(x,z | g)"
  )
  # each primitive row summarizes its own backend coordinate
  expect_true(all(primitives %in% rownames(raw)))
  expect_false(any(grepl("lkj_u[", rownames(raw), fixed = TRUE)))
  draws <- do.call(rbind, lapply(fit$mcmc, as.matrix))
  expect_equal(
    raw[primitives, "Mean"],
    unname(colMeans(draws[, paste0("mu__xREx__g_xRE_CORx_lkj_u[", 1:3, "]")])),
    tolerance = 1e-12
  )
  expect_identical(
    rownames(JAGS_estimates_table(fit, random_effects_summary = "raw",
                                  formula_prefix = FALSE,
                                  remove_diagnostics = TRUE))[
      match(primitives, rownames(raw))
    ],
    sub("(mu) ", "", primitives, fixed = TRUE)
  )
  # raw rows are backend coordinates, not catalog quantities
  catalog <- parameter_catalog(fit)
  for(row in primitives){
    expect_error(
      parameter_catalog_resolve(catalog, row),
      class = "BayesTools_parameter_not_found"
    )
  }
})

test_that("semantic tables omit the allocation shares of ordered-factor priors", {

  fit <- .label_cached_fit("fit_label_ordered")
  coordinates <- parameter_coordinates(fit)
  shares <- coordinates$coordinate_name[
    coordinates$internal & coordinates$role == "parameter"
  ]
  expect_length(shares, 2L)
  expected <- c("(mu) intercept", "(mu) f[mid]", "(mu) f{2}", "(mu) f_ordered_total")
  for(mode in c("standard", "full")){
    expect_identical(
      rownames(JAGS_estimates_table(fit, random_effects_summary = mode,
                                    remove_diagnostics = TRUE)),
      expected,
      info = mode
    )
  }
  # raw tables show every backend coordinate
  expect_identical(
    rownames(JAGS_estimates_table(fit, random_effects_summary = "raw",
                                  remove_diagnostics = TRUE)),
    c(expected, shares)
  )
})

test_that("original-scale random-effect SDs require the fitted design", {

  fit <- .label_cached_fit("fit_label_random_slope")
  fitted_scale <- attr(fit, "formula_scale")
  hand_built <- list(mu = list(mu_x = list(
    mean = fitted_scale$mu$mu_x$mean,
    sd   = fitted_scale$mu$mu_x$sd
  )))
  parameters <- c("mu__xREx__g_intercept", "mu__xREx__g_x")
  mixed <- as_mixed_posteriors(fit, parameters = parameters)

  # the SD columns alone: standardization information without the fitted
  # design (and the random-effect structure it comes with) stops instead of
  # leaving the SDs standardized
  expect_error(
    ensemble_estimates_table(mixed, parameters = parameters,
                             transform_scaled = TRUE, formula_scale = hand_built),
    class = "BayesTools_formula_transform_unavailable"
  )
  # the fitted object's formula_scale transforms them as the model table does
  ensemble <- ensemble_estimates_table(mixed, parameters = parameters,
                                       transform_scaled = TRUE,
                                       formula_scale = fitted_scale)
  scaled_table <- JAGS_estimates_table(fit, transform_scaled = TRUE,
                                       remove_diagnostics = TRUE)
  expect_equal(ensemble[, "Mean"], scaled_table[rownames(ensemble), "Mean"],
               tolerance = 1e-10)
})

test_that("standardization passed for a fitted model takes its fitted structure", {

  # the same means and SDs as the fitted formula_scale, built by hand: the
  # fitted design, the random-effect SD structure, and the log intercept come
  # from the fitted model
  hand_built <- function(fit){
    fitted_scale <- attr(fit, "formula_scale")
    list(mu = list(mu_x = list(
      mean = fitted_scale$mu$mu_x$mean,
      sd   = fitted_scale$mu$mu_x$sd
    )))
  }
  expect_same_transforms <- function(fit){
    expect_identical(
      transform_scale_samples(fit, formula_scale = hand_built(fit)),
      transform_scale_samples(fit)
    )
    expect_identical(
      transform_prior_samples(fit, n_samples = 200, seed = 1,
                              formula_scale = hand_built(fit)),
      transform_prior_samples(fit, n_samples = 200, seed = 1)
    )
  }

  # a random slope of a standardized predictor
  expect_same_transforms(.label_cached_fit("fit_label_random_slope"))

  # a log intercept
  fit <- .label_cached_fit("fit_label_log_intercept")
  expect_true(isTRUE(attr(attr(fit, "formula_scale")$mu, "log_intercept")))
  expect_same_transforms(fit)
})

test_that("random-effect SD mixed columns carry the catalog's random-effect labels", {

  # a random treatment slope: every SD coordinate is its catalog SD quantity
  fit <- .label_cached_fit("fit_label_random_factor_slope")
  parameters <- c("mu__xREx__g_intercept", "mu__xREx__g_f")
  mixed <- as_mixed_posteriors(fit, parameters = parameters)
  catalog <- parameter_catalog(fit)
  quantities <- rbind(
    posterior_metadata(mixed$mu__xREx__g_intercept, "quantities"),
    posterior_metadata(mixed$mu__xREx__g_f, "quantities")
  )
  labels <- c("(mu) sd(intercept)", "(mu) sd(f[b])", "(mu) sd(f[c])")
  expect_identical(parameter_labels(quantities, "table"), labels)
  expect_identical(
    quantities$quantity_id,
    vapply(labels, function(label){
      parameter_catalog_resolve(catalog, label)$quantities$quantity_id
    }, character(1), USE.NAMES = FALSE)
  )
  # the mixed columns keep their names
  expect_identical(colnames(mixed$mu__xREx__g_f), c("mu__xREx__g_f[b]", "mu__xREx__g_f[c]"))
  ensemble <- ensemble_estimates_table(mixed, parameters = parameters)
  expect_identical(rownames(ensemble), labels)
  standard <- JAGS_estimates_table(fit, remove_diagnostics = TRUE)
  expect_equal(ensemble[labels, "Mean"], standard[labels, "Mean"], tolerance = 1e-10)

  # a random slope of a standardized predictor: the SD columns are transformed
  # by the random-effect covariance transform, not the fixed-effect design
  fit <- .label_cached_fit("fit_label_random_slope")
  parameters <- c("mu_intercept", "mu_x", "mu__xREx__g_intercept", "mu__xREx__g_x")
  mixed <- as_mixed_posteriors(fit, parameters = parameters)
  labels <- c("(mu) sd(intercept)", "(mu) sd(x)")
  ensemble <- ensemble_estimates_table(
    mixed, parameters = parameters,
    transform_scaled = TRUE, formula_scale = attr(fit, "formula_scale")
  )
  expect_identical(rownames(ensemble)[3:4], labels)
  scaled_table <- JAGS_estimates_table(fit, transform_scaled = TRUE, remove_diagnostics = TRUE)
  expect_equal(ensemble[labels, "Mean"], scaled_table[labels, "Mean"], tolerance = 1e-10)
  expect_equal(
    ensemble[c("(mu) intercept", "(mu) x"), "Mean"],
    scaled_table[c("(mu) intercept", "(mu) x"), "Mean"],
    tolerance = 1e-10
  )
})

test_that("random-effect SD columns describe the scale of the draws they hold", {

  # a standardized random slope: the catalog SD quantities are original-scale
  # functions of the block's fitted coordinates
  fit <- .label_cached_fit("fit_label_random_slope")
  formula_scale <- attr(fit, "formula_scale")
  catalog <- parameter_catalog(fit)
  parameters <- c("mu__xREx__g_intercept", "mu__xREx__g_x")
  sd_quantities <- function(samples){
    rbind(
      posterior_metadata(samples$mu__xREx__g_intercept, "quantities"),
      posterior_metadata(samples$mu__xREx__g_x, "quantities")
    )
  }
  fitted_labels   <- c("(mu) g: sd(intercept)", "(mu) g: sd(x)")
  original_labels <- c("(mu) sd(intercept)", "(mu) sd(x)")
  original_ids <- vapply(original_labels, function(label){
    parameter_catalog_resolve(catalog, label)$quantities$quantity_id
  }, character(1), USE.NAMES = FALSE)
  raw_table    <- JAGS_estimates_table(fit, random_effects_summary = "raw",
                                       remove_diagnostics = TRUE)
  scaled_table <- JAGS_estimates_table(fit, transform_scaled = TRUE,
                                       remove_diagnostics = TRUE)

  # on the fitted scale the columns are their fitted coordinates: labelled and
  # valued as the raw table shows the same draws, and no catalog quantity
  mixed <- as_mixed_posteriors(fit, parameters = parameters)
  expect_identical(parameter_labels(sd_quantities(mixed), "table"), fitted_labels)
  expect_identical(sd_quantities(mixed)$quantity_id, c("", ""))
  ensemble <- ensemble_estimates_table(mixed, parameters = parameters)
  expect_identical(rownames(ensemble), fitted_labels)
  expect_equal(ensemble[fitted_labels, "Mean"], raw_table[fitted_labels, "Mean"],
               tolerance = 1e-10)

  # transformed to the original scale they are the catalog's SD quantities,
  # labelled and valued as the model table shows them
  transformed <- .transform_scale_samples_list(mixed, formula_scale)
  expect_identical(parameter_labels(sd_quantities(transformed), "table"), original_labels)
  expect_identical(sd_quantities(transformed)$quantity_id, original_ids)
  ensemble <- ensemble_estimates_table(mixed, parameters = parameters,
                                       transform_scaled = TRUE,
                                       formula_scale = formula_scale)
  expect_identical(rownames(ensemble), original_labels)
  expect_equal(ensemble[original_labels, "Mean"], scaled_table[original_labels, "Mean"],
               tolerance = 1e-10)

  # draws created on the original scale carry the original-scale quantities
  mixed_original <- as_mixed_posteriors(fit, parameters = parameters,
                                        transform_scaled = TRUE)
  expect_identical(parameter_labels(sd_quantities(mixed_original), "table"),
                   original_labels)
  expect_identical(sd_quantities(mixed_original)$quantity_id, original_ids)
  expect_equal(
    vapply(parameters, function(parameter) mean(mixed_original[[parameter]]),
           numeric(1), USE.NAMES = FALSE),
    scaled_table[original_labels, "Mean"],
    tolerance = 1e-10
  )

  # a random slope without its intercept: the original-scale SD is a
  # one-to-one transform of the coordinate, which is not that quantity either
  fit <- .label_cached_fit("fit_label_random_slope_only")
  mixed <- as_mixed_posteriors(fit, parameters = "mu__xREx__g_x")
  expect_identical(
    rownames(ensemble_estimates_table(mixed, parameters = "mu__xREx__g_x")),
    "(mu) g: sd(x)"
  )
  expect_identical(
    rownames(ensemble_estimates_table(mixed, parameters = "mu__xREx__g_x",
                                      transform_scaled = TRUE,
                                      formula_scale = attr(fit, "formula_scale"))),
    "(mu) sd(x)"
  )
})

test_that("original-scale random-effect SDs require the random structure of the passed scale", {

  transform_reason <- function(expr){
    condition <- tryCatch(expr, BayesTools_formula_transform_unavailable = function(e) e)
    expect_s3_class(condition, "BayesTools_formula_transform_unavailable")
    condition$reason
  }
  fit_slope     <- .label_cached_fit("fit_label_random_slope")
  fit_intercept <- .label_cached_fit("fit_label_random_intercept")
  scale_slope     <- attr(fit_slope, "formula_scale")
  scale_intercept <- attr(fit_intercept, "formula_scale")
  parameters <- c("mu_intercept", "mu_x", "mu__xREx__g_intercept", "mu__xREx__g_x")
  mixed <- as_mixed_posteriors(fit_slope, parameters = parameters)

  # the formula_scale of the model without the random slope does not contain
  # its SD, which would stay standardized (0.615 for 0.792 and 0.204 for 0.121)
  expect_identical(transform_reason(ensemble_estimates_table(
    mixed, parameters = parameters, transform_scaled = TRUE,
    formula_scale = scale_intercept
  )), "random_effects_outside_structure")
  # the intercept SD alone is an SD of that structure, but of a block it
  # leaves unchanged, while the samples' original-scale SD differs
  intercept_parameters <- c("mu_intercept", "mu_x", "mu__xREx__g_intercept")
  mixed_intercept_sd <- as_mixed_posteriors(fit_slope, parameters = intercept_parameters)
  expect_identical(transform_reason(ensemble_estimates_table(
    mixed_intercept_sd, parameters = intercept_parameters, transform_scaled = TRUE,
    formula_scale = scale_intercept
  )), "random_effect_structure_differs")
  # and the structure of the slope model needs the slope SD with it
  expect_identical(transform_reason(ensemble_estimates_table(
    mixed_intercept_sd, parameters = intercept_parameters, transform_scaled = TRUE,
    formula_scale = scale_slope
  )), "random_effect_block_incomplete")

  # matching structures transform as the model tables do
  ensemble <- ensemble_estimates_table(mixed, parameters = parameters,
                                       transform_scaled = TRUE,
                                       formula_scale = scale_slope)
  scaled_table <- JAGS_estimates_table(fit_slope, transform_scaled = TRUE,
                                       remove_diagnostics = TRUE)
  expect_equal(ensemble[, "Mean"], scaled_table[rownames(ensemble), "Mean"],
               tolerance = 1e-10)
  mixed_intercept <- as_mixed_posteriors(fit_intercept, parameters = intercept_parameters)
  ensemble <- ensemble_estimates_table(mixed_intercept, parameters = intercept_parameters,
                                       transform_scaled = TRUE,
                                       formula_scale = scale_intercept)
  scaled_table <- JAGS_estimates_table(fit_intercept, transform_scaled = TRUE,
                                       remove_diagnostics = TRUE)
  expect_identical(rownames(ensemble)[3L], "(mu) sd(intercept)")
  expect_equal(ensemble[, "Mean"], scaled_table[rownames(ensemble), "Mean"],
               tolerance = 1e-10)

  # a variance-allocation model: its intercept SDs are allocation-derived
  # nodes, which its structure contains, but not the random slope
  fit_allocation <- .label_cached_fit("fit_label_allocation")
  scale_allocation <- attr(fit_allocation, "formula_scale")
  expect_identical(transform_reason(ensemble_estimates_table(
    mixed, parameters = parameters, transform_scaled = TRUE,
    formula_scale = scale_allocation
  )), "random_effects_outside_structure")
  # its own samples (coefficients and allocation coordinates) transform with
  # its own structure as before: the allocation coordinates are not
  # standardized, the coefficients as in the model table
  allocation_parameters <- names(attr(fit_allocation, "prior_list"))
  mixed_allocation <- as_mixed_posteriors(fit_allocation, parameters = allocation_parameters)
  untransformed <- ensemble_estimates_table(mixed_allocation, parameters = allocation_parameters)
  ensemble <- ensemble_estimates_table(mixed_allocation, parameters = allocation_parameters,
                                       transform_scaled = TRUE,
                                       formula_scale = scale_allocation)
  scaled_table <- JAGS_estimates_table(fit_allocation, transform_scaled = TRUE,
                                       remove_diagnostics = TRUE)
  coefficients <- c("(mu) intercept", "(mu) x")
  expect_equal(ensemble[coefficients, "Mean"], scaled_table[coefficients, "Mean"],
               tolerance = 1e-10)
  allocation_rows <- setdiff(rownames(ensemble), coefficients)
  expect_true(length(allocation_rows) > 0L)
  expect_equal(ensemble[allocation_rows, "Mean"], untransformed[allocation_rows, "Mean"],
               tolerance = 1e-12)
})
