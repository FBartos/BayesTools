skip_if_not_test_profile("unit")

.formula_coefficient_density_fit <- function(formula_result,
                                             sampled_columns){

  design <- formula_result$formula_design
  parameter <- design$parameter
  formula_design <- stats::setNames(list(design), parameter)
  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  formula_scale <- NULL
  if("formula_scale" %in% names(formula_result)){
    formula_scale <- stats::setNames(
      list(formula_result$formula_scale),
      parameter
    )
    attr(fit, "formula_scale") <- formula_scale
  }
  attr(fit, "formula_design") <- formula_design
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = sampled_columns,
    prior_list = formula_result$prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  .bt_attach_fit_contract(fit)
}

.formula_coefficient_source_names <- function(formula_result){

  unique(unlist(Map(
    .prior_linear_prior_columns,
    names(formula_result$prior_list),
    formula_result$prior_list
  ), use.names = FALSE))
}

test_that("formula coefficient transforms expose the sample transformation", {

  formula_result <- JAGS_formula(
    formula = ~ 1 + x * z,
    parameter = "mu",
    data = data.frame(x = c(2, 4, 6), z = c(-1, 1, 3)),
    prior_list = list(
      intercept = prior("point", list(3)),
      x = prior("normal", list(0, 1)),
      z = prior("normal", list(0, 1)),
      `x:z` = prior("normal", list(0, 1))
    ),
    formula_scale = TRUE
  )
  source_names <- .formula_coefficient_source_names(formula_result)
  sampled_names <- setdiff(source_names, "mu_intercept")
  fit <- .formula_coefficient_density_fit(formula_result, sampled_names)
  transform <- JAGS_formula_coefficient_transform(fit, "mu")

  expect_s3_class(transform, "BayesTools_formula_coefficient_transform")
  expect_identical(transform$schema_version, 1L)
  expect_identical(transform$formula_design_version, 4L)
  expect_identical(transform$parameter_map_version, 4L)
  expect_identical(transform$source_names, source_names)
  expect_identical(transform$target_names, source_names)
  expect_identical(
    transform$matrix["mu_intercept", ],
    c(mu_intercept = 1, mu_x = -2, mu_z = -0.5,
      mu_x__xXx__z = 1)
  )
  expect_identical(
    transform$matrix["mu_x", ],
    c(mu_intercept = 0, mu_x = 0.5, mu_z = 0,
      mu_x__xXx__z = -0.25)
  )
  expect_identical(
    transform$dependencies$source[
      transform$dependencies$target == "mu_x"
    ],
    c("mu_x", "mu_x__xXx__z")
  )
  expect_identical(
    transform$sources$monitor_status,
    c("structural", "sampled", "sampled", "sampled")
  )
  expect_identical(
    transform$targets$structural_status,
    rep("dependent", 4L)
  )

  samples <- matrix(
    c(0.2, -0.5, 0.75, -1, 2, 0.4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, sampled_names)
  )
  source_samples <- cbind(mu_intercept = 3, samples)
  source_samples <- source_samples[, transform$source_names, drop = FALSE]
  expected <- source_samples %*% t(transform$matrix)
  transformed <- transform_scale_samples(
    samples,
    formula_scale = list(mu = formula_result$formula_scale)
  )
  expect_equal(
    transformed[, transform$target_names, drop = FALSE],
    expected,
    tolerance = 1e-14
  )
})

.formula_coefficient_sample_fit <- function(formula_result, samples){

  design <- formula_result$formula_design
  parameter <- design$parameter
  formula_design <- stats::setNames(list(design), parameter)
  formula_scale <- stats::setNames(
    list(formula_result$formula_scale),
    parameter
  )
  fit <- coda::mcmc(samples)
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_scale") <- formula_scale
  attr(fit, "formula_design") <- formula_design
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = colnames(samples),
    prior_list = formula_result$prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  .bt_attach_fit_contract(fit)
}

.formula_coefficient_design_matrix <- function(formula_result){

  data_names <- names(formula_result$data)[
    startsWith(names(formula_result$data), "mu_data_")
  ]
  cbind(1, do.call(cbind, lapply(formula_result$data[data_names], as.matrix)))
}

test_that("formula coefficient transforms follow the fitted design for nested slopes", {

  # ~ f/x fits one slope per level (full indicator coding of f:x) next to a
  # treatment-coded f, so f:x column k is not coded like f column k.
  data <- data.frame(
    x = c(1, 3, 7, 2, 6, 11, 4, 5, 9),
    f = factor(rep(c("A", "B", "C"), each = 3L), levels = c("A", "B", "C"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    f = prior_factor("normal", list(0, 1), contrast = "treatment"),
    "f:x" = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  scaled <- JAGS_formula(~ f/x, "mu", data, prior_list,
                         formula_scale = list(x = TRUE))
  original <- JAGS_formula(~ f/x, "mu", data, prior_list)
  source_names <- .formula_coefficient_source_names(scaled)
  expect_identical(
    source_names,
    c("mu_intercept", "mu_f[1]", "mu_f[2]",
      "mu_f__xXx__x[1]", "mu_f__xXx__x[2]", "mu_f__xXx__x[3]")
  )
  fit <- .formula_coefficient_density_fit(scaled, source_names)
  transform <- JAGS_formula_coefficient_transform(fit, "mu")

  # Analytic map: a_A = b0 - g_A m / s, (f=B) = b_B - (g_B - g_A) m / s,
  # (f=C) = b_C - (g_C - g_A) m / s, slopes g_k / s.
  m <- mean(data$x)
  s <- stats::sd(data$x)
  expected <- matrix(0, 6L, 6L, dimnames = list(source_names, source_names))
  expected["mu_intercept", c("mu_intercept", "mu_f__xXx__x[1]")] <- c(1, -m / s)
  expected["mu_f[1]", c("mu_f[1]", "mu_f__xXx__x[1]", "mu_f__xXx__x[2]")] <-
    c(1, m / s, -m / s)
  expected["mu_f[2]", c("mu_f[2]", "mu_f__xXx__x[1]", "mu_f__xXx__x[3]")] <-
    c(1, m / s, -m / s)
  for(k in 1:3){
    slope <- paste0("mu_f__xXx__x[", k, "]")
    expected[slope, slope] <- 1 / s
  }
  expect_equal(transform$matrix, expected, tolerance = 1e-12)
  expect_identical(transform$matrix != 0, expected != 0)

  # The fitted linear predictor is reproduced on the original scale; the
  # name-paired map of the review scenario missed it by more than 2.
  coefficients <- rbind(
    c(0.3, -0.7, 1.1, 0.25, -0.4, 0.9),
    c(-1.2, 0.5, 0.2, -0.6, 0.35, 0.15)
  )
  colnames(coefficients) <- source_names
  expect_equal(
    .formula_coefficient_design_matrix(original) %*%
      (transform$matrix %*% t(coefficients)),
    .formula_coefficient_design_matrix(scaled) %*% t(coefficients),
    tolerance = 1e-12
  )

  # Posterior transformation of a fitted object uses the same design map.
  sample_fit <- .formula_coefficient_sample_fit(scaled, coefficients)
  expect_equal(
    transform_scale_samples(sample_fit)[, source_names],
    coefficients %*% t(transform$matrix),
    tolerance = 1e-12
  )
})

test_that("JAGS_formula formula-scale metadata carry the fitted design", {

  data <- data.frame(
    x = c(1, 3, 7, 2, 6, 11, 4, 5, 9),
    f = factor(rep(c("A", "B", "C"), each = 3L), levels = c("A", "B", "C"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    f = prior_factor("normal", list(0, 1), contrast = "treatment"),
    "f:x" = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  scaled <- JAGS_formula(~ f/x, "mu", data, prior_list,
                         formula_scale = list(x = TRUE))
  original <- JAGS_formula(~ f/x, "mu", data, prior_list)
  expect_false(is.null(attr(scaled$formula_scale, "unscale_design")))
  # the design keeps its fitted formula-scale copy (compared on bridge rebuilds)
  expect_null(attr(scaled$formula_design$formula_scale, "unscale_design"))
  expect_null(original$formula_scale)

  # formula-scale-only consumers: a posterior matrix and the stored metadata
  source_names <- .formula_coefficient_source_names(scaled)
  coefficients <- rbind(
    c(0.3, -0.7, 1.1, 0.25, -0.4, 0.9),
    c(-1.2, 0.5, 0.2, -0.6, 0.35, 0.15)
  )
  colnames(coefficients) <- source_names
  transformed <- transform_scale_samples(
    coefficients,
    formula_scale = list(mu = scaled$formula_scale)
  )
  expect_equal(
    .formula_coefficient_design_matrix(original) %*% t(transformed[, source_names]),
    .formula_coefficient_design_matrix(scaled) %*% t(coefficients),
    tolerance = 1e-12
  )

  # Formula-scale metadata stored before the design was attached keep the
  # name-paired map without error.
  legacy_scale <- scaled$formula_scale
  attr(legacy_scale, "unscale_design") <- NULL
  expect_identical(
    transform_scale_samples(coefficients, formula_scale = list(mu = legacy_scale)),
    coefficients %*% t(.build_unscale_matrix_by_names(
      source_names, legacy_scale, "mu", require_closure = FALSE
    ))
  )
})

test_that("design-derived unscaling keeps the coding of logical predictors", {

  # Unstandardized logical predictors are indicator-coded ('lgTRUE'); the
  # synthetic verification design must reproduce that coding. Standardized
  # logical predictors are numeric in the fitted design.
  data <- data.frame(
    x = c(1, 3, 7, 2, 6, 11, 4, 5, 9, 8),
    lg = rep(c(TRUE, FALSE), 5L)
  )
  normal <- prior("normal", list(0, 1))
  prior_list <- list(intercept = normal, x = normal, lg = normal, "x:lg" = normal)
  coefficients <- rbind(
    c(0.3, -0.7, 1.1, 0.25),
    c(-1.2, 0.5, 0.2, -0.6)
  )

  for(formula_scale in list(list(x = TRUE), TRUE)){
    scaled <- JAGS_formula(~ x * lg, "mu", data, prior_list,
                           formula_scale = formula_scale)
    original <- JAGS_formula(~ x * lg, "mu", data, prior_list)
    source_names <- .formula_coefficient_source_names(scaled)
    colnames(coefficients) <- source_names
    transformed <- transform_scale_samples(
      coefficients,
      formula_scale = list(mu = scaled$formula_scale)
    )
    expect_equal(
      .formula_coefficient_design_matrix(original) %*% t(transformed[, source_names]),
      .formula_coefficient_design_matrix(scaled) %*% t(coefficients),
      tolerance = 1e-12
    )
    transform <- JAGS_formula_coefficient_transform(
      .formula_coefficient_density_fit(scaled, source_names),
      "mu"
    )
    expect_equal(
      transform$matrix,
      .build_unscale_matrix_by_names(
        source_names, scaled$formula_scale, "mu", require_closure = FALSE
      ),
      tolerance = 1e-14
    )
  }
})

test_that("formula coefficient transforms refuse terms whose centering is not representable", {

  data <- data.frame(
    x = c(1, 3, 7, 2, 6, 11, 4, 5, 9),
    f = factor(rep(c("A", "B", "C"), each = 3L), levels = c("A", "B", "C")),
    d = c(0, 1, 0, 1, 1, 0, 0, 1, 1)
  )
  cases <- list(
    list(
      formula = ~ x + x:f,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        "x:f" = prior_factor("normal", list(0, 1), contrast = "treatment")
      ),
      term = "x:f"
    ),
    list(
      formula = ~ x + x:d,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        "x:d" = prior("normal", list(0, 1))
      ),
      term = "x:d"
    )
  )

  for(case in cases){
    formula_result <- JAGS_formula(
      case$formula, "mu", data, case$prior_list,
      formula_scale = list(x = TRUE)
    )
    source_names <- .formula_coefficient_source_names(formula_result)
    fit <- .formula_coefficient_density_fit(formula_result, source_names)
    error <- tryCatch(JAGS_formula_coefficient_transform(fit, "mu"),
                      error = identity)
    expect_s3_class(error, "BayesTools_formula_transform_unavailable")
    expect_identical(error$reason, "original_scale_not_representable")
    expect_identical(error$terms, case$term)
    expect_match(conditionMessage(error), case$term, fixed = TRUE)
    expect_match(conditionMessage(error), "does not contain", fixed = TRUE)

    samples <- matrix(
      c(0.5, 0.2, -0.3, 0.4, 1, -0.2, 0.1, 0.3)[seq_len(2L * length(source_names))],
      nrow = 2L,
      dimnames = list(NULL, source_names)
    )
    sample_fit <- .formula_coefficient_sample_fit(formula_result, samples)
    expect_error(
      transform_scale_samples(sample_fit),
      "does not contain"
    )
  }
})

test_that("design-derived unscaling accepts exact maps under extreme centering", {

  # Centering products up to (m / s)^3 ~ 2e8: the exact name-paired map
  # reproduces the standardized design only to ~1e-7 in floating point, which
  # an absolute 1e-8 bound would reject as "not representable".
  offsets <- c(-0.9, -0.5, -0.2, 0, 0.1, 0.3, 0.6, 1.0, -0.7, 0.4, -0.1, 0.8)
  data <- data.frame(
    year = 2010 + 0.5 * offsets,
    x = 100 + rev(offsets),
    z = 100 + offsets[c(4:12, 1:3)]
  )
  normal <- prior("normal", list(0, 1))
  scale <- list(year = TRUE, x = TRUE, z = TRUE)
  full_priors <- list(
    intercept = normal, year = normal, x = normal, z = normal,
    "year:x" = normal, "year:z" = normal, "x:z" = normal, "year:x:z" = normal
  )
  formula_result <- JAGS_formula(~ year * x * z, "mu", data, full_priors,
                                 formula_scale = scale)
  source_names <- .formula_coefficient_source_names(formula_result)
  fit <- .formula_coefficient_density_fit(formula_result, source_names)
  transform <- JAGS_formula_coefficient_transform(fit, "mu")

  expect_identical(
    transform$matrix,
    .build_unscale_matrix_by_names(
      source_names, formula_result$formula_scale, "mu",
      require_closure = FALSE
    )
  )
  original <- JAGS_formula(~ year * x * z, "mu", data, full_priors)
  coefficients <- c(0.3, -0.7, 1.1, 0.25, -0.4, 0.9, 0.2, -0.15)
  # absolute error of X_o (A b): about eps * |X_o| |A| |b| < 1e-6 here
  expect_equal(
    .formula_coefficient_design_matrix(original) %*%
      (transform$matrix %*% coefficients),
    .formula_coefficient_design_matrix(formula_result) %*% coefficients,
    tolerance = 1e-6
  )

  # The least-squares path (nested slopes, m / s ~ 7000) stays exact and keeps
  # its structural zeros.
  nested_data <- data.frame(
    x = 2010 + 0.5 * offsets[1:9],
    f = factor(rep(c("A", "B", "C"), each = 3L), levels = c("A", "B", "C"))
  )
  nested_priors <- list(
    intercept = normal,
    f = prior_factor("normal", list(0, 1), contrast = "treatment"),
    "f:x" = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  nested <- JAGS_formula(~ f/x, "mu", nested_data, nested_priors,
                         formula_scale = list(x = TRUE))
  nested_names <- .formula_coefficient_source_names(nested)
  nested_transform <- JAGS_formula_coefficient_transform(
    .formula_coefficient_density_fit(nested, nested_names),
    "mu"
  )
  m <- mean(nested_data$x)
  s <- stats::sd(nested_data$x)
  expected <- matrix(0, 6L, 6L, dimnames = list(nested_names, nested_names))
  expected["mu_intercept", c("mu_intercept", "mu_f__xXx__x[1]")] <- c(1, -m / s)
  expected["mu_f[1]", c("mu_f[1]", "mu_f__xXx__x[1]", "mu_f__xXx__x[2]")] <-
    c(1, m / s, -m / s)
  expected["mu_f[2]", c("mu_f[2]", "mu_f__xXx__x[1]", "mu_f__xXx__x[3]")] <-
    c(1, m / s, -m / s)
  for(k in 1:3){
    slope <- paste0("mu_f__xXx__x[", k, "]")
    expected[slope, slope] <- 1 / s
  }
  expect_equal(nested_transform$matrix, expected, tolerance = 1e-10)
  expect_identical(nested_transform$matrix != 0, expected != 0)

  # Omitting x:z, which centering year in year:x:z induces, stays rejected
  # at the same scaling.
  reduced_priors <- full_priors[names(full_priors) != "x:z"]
  reduced <- JAGS_formula(~ year * x * z - x:z, "mu", data, reduced_priors,
                          formula_scale = scale)
  reduced_fit <- .formula_coefficient_density_fit(
    reduced,
    .formula_coefficient_source_names(reduced)
  )
  error <- tryCatch(JAGS_formula_coefficient_transform(reduced_fit, "mu"),
                    error = identity)
  expect_s3_class(error, "BayesTools_formula_transform_unavailable")
  expect_identical(error$reason, "original_scale_not_representable")
  expect_identical(error$terms, "year:x:z")
})

test_that("design-derived unscaling applies only to fitted coefficient coordinates", {

  data <- data.frame(
    x = c(1, 3, 7, 2, 6, 11, 4, 5, 9),
    f = factor(rep(c("a", "b", "c"), 3L), levels = c("a", "b", "c"))
  )
  formula_result <- JAGS_formula(
    ~ x * f, "mu", data,
    list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "x:f" = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    formula_scale = list(x = TRUE)
  )
  with_design <- list(mu = formula_result$formula_scale)
  expect_false(is.null(attr(with_design$mu, "unscale_design")))
  without_design <- with_design
  attr(without_design$mu, "unscale_design") <- NULL

  # Level-wise summaries (one column per level, or level-labelled columns)
  # are not the fitted coefficient vector and keep the name-paired map.
  level_columns <- list(
    matrix(c(0, 0.4, -0.2, 0, 0.1, 0.3), nrow = 2, byrow = TRUE,
           dimnames = list(NULL, paste0("mu_x__xXx__f[", 1:3, "]"))),
    cbind(`mu_f[b]` = c(0.5, 1), `mu_f[c]` = c(-1, 0.2),
          `mu_x__xXx__f[b]` = c(0.3, 0.1), `mu_x__xXx__f[c]` = c(-0.4, 0.2))
  )
  for(samples in level_columns){
    expect_identical(
      transform_scale_samples(samples, formula_scale = with_design),
      transform_scale_samples(samples, formula_scale = without_design)
    )
  }

  # The fitted coefficient coordinates use the verified design map, which
  # equals the name-paired map for this crossed formula.
  coefficients <- matrix(
    c(1, 0.5, 0.2, -0.3, 0.4, -0.1),
    nrow = 1,
    dimnames = list(NULL, .formula_coefficient_source_names(formula_result))
  )
  expect_equal(
    transform_scale_samples(coefficients, formula_scale = with_design),
    transform_scale_samples(coefficients, formula_scale = without_design),
    tolerance = 1e-14
  )
})

test_that("formula prior densities of raw coefficients ignore their multiply_by", {

  # As in fixture fit_complex_mixed: the monitored mu_x node is the raw
  # coefficient with a spike-and-slab N(0, 1) x Spike(0.5) prior; 'sigma'
  # multiplies only its linear-predictor contribution.
  x_prior <- prior_spike_and_slab(prior("normal", list(0, 1), prior_weights = 1))
  attr(x_prior, "multiply_by") <- "sigma"
  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = data.frame(x = c(-1, 0.5, 2)),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = x_prior
    )
  )
  expect_identical(attr(formula_result$prior_list$mu_x, "multiply_by"), "sigma")
  source_names <- .formula_coefficient_source_names(formula_result)
  fit <- .formula_coefficient_density_fit(formula_result, source_names)
  attr(fit, "prior_list") <- c(
    attr(fit, "prior_list"),
    list(sigma = prior("lognormal", list(0, 1)))
  )

  density <- JAGS_formula_prior_density(fit, parameter = "mu", target = "mu_x")
  continuous <- density$density
  expect_equal(continuous$mass, 0.5, tolerance = 1e-12)
  expect_equal(continuous$mass * max(continuous$y), 0.5 * stats::dnorm(0),
               tolerance = 1e-3)
  expect_lt(max(abs(continuous$x)), 10)

  ordinate <- prior_density_ordinate(density, 1)
  expect_identical(ordinate$behavior, "regular")
  expect_true(ordinate$exact)
  expect_equal(ordinate$log_density, log(0.5) + stats::dnorm(1, log = TRUE),
               tolerance = 1e-12)
  null_ordinate <- prior_density_ordinate(density, 0)
  expect_identical(null_ordinate$behavior, "point_mass")
  expect_equal(null_ordinate$point_mass, 0.5, tolerance = 1e-12)

  # A supplied (e.g., conditional or model-mixture) context is treated alike.
  context <- .prior_density_build_context(
    prior_list = attr(fit, "prior_list"),
    column_names = c(source_names, "sigma")
  )
  supplied <- JAGS_formula_prior_density(
    fit, parameter = "mu", target = "mu_x", context = context
  )
  expect_equal(
    prior_density_ordinate(supplied, 1)$log_density,
    log(0.5) + stats::dnorm(1, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("formula coefficient transforms require current linked schemas", {

  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = data.frame(x = c(-1, 0, 1)),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )
  )
  source_names <- .formula_coefficient_source_names(formula_result)
  fit <- .formula_coefficient_density_fit(formula_result, source_names)
  transform <- JAGS_formula_coefficient_transform(fit, "mu")

  stale_design <- transform
  stale_design$formula_design_version <-
    stale_design$formula_design_version - 1L
  expect_error(
    BayesTools:::.bt_validate_formula_coefficient_transform(stale_design),
    "missing or unsupported"
  )

  unknown_map <- transform
  unknown_map$parameter_map_version <- unknown_map$parameter_map_version + 1L
  expect_error(
    BayesTools:::.bt_validate_formula_coefficient_transform(unknown_map),
    "missing or unsupported"
  )
})

test_that("formula coefficient transforms exclude random-effect priors", {

  data <- data.frame(
    x = c(-1, 0, 1, -1, 0, 1),
    id = factor(rep(c("a", "b"), each = 3L))
  )
  compile_policies <- list(
    sampled = NULL,
    marginalized = random_effects_compile(marginalized = "block")
  )

  for(policy in compile_policies){
    formula_result <- JAGS_formula(
      formula = ~ 1 + x +
        random(1 | id, name = "block", covariance = "diag"),
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      formula_scale = TRUE,
      prior_random = prior_random(
        block = random_block(sd = prior("gamma", list(2, 2)))
      ),
      random_effects_compile = policy
    )
    all_prior_coordinates <- .formula_coefficient_source_names(formula_result)
    expect_true(any(grepl("__xREx__", all_prior_coordinates, fixed = TRUE)))
    fit <- .formula_coefficient_density_fit(
      formula_result,
      all_prior_coordinates
    )

    transform <- JAGS_formula_coefficient_transform(fit, "mu")
    expect_identical(
      transform$source_names,
      c("mu_intercept", "mu_x")
    )
    expect_false(any(grepl("__xREx__", transform$source_names, fixed = TRUE)))

    density <- JAGS_formula_prior_density(
      fit,
      parameter = "mu",
      target = "mu_x"
    )
    ordinate <- prior_density_ordinate(density, 0)
    expect_identical(ordinate$behavior, "regular")
    expect_equal(
      ordinate$log_density,
      stats::dnorm(
        0,
        sd = abs(transform$matrix["mu_x", "mu_x"]),
        log = TRUE
      ),
      tolerance = 1e-12
    )
  }
})

test_that("formula prior density names a single fixed source", {

  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 | id, name = "block", covariance = "diag"),
    parameter = "mu",
    data = data.frame(id = factor(c("a", "a", "b", "b"))),
    prior_list = list(
      intercept = prior("normal", list(0, 2))
    ),
    prior_random = prior_random(
      block = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  all_prior_coordinates <- .formula_coefficient_source_names(formula_result)
  fit <- .formula_coefficient_density_fit(
    formula_result,
    all_prior_coordinates
  )

  density <- JAGS_formula_prior_density(
    fit,
    parameter = "mu",
    target = "mu_intercept"
  )
  ordinate <- prior_density_ordinate(density, 0)
  expect_identical(ordinate$behavior, "regular")
  expect_equal(
    ordinate$log_density,
    stats::dnorm(0, sd = 2, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("formula prior densities distinguish structural and dependent targets", {

  continuous_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = prior("point", list(3)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = TRUE
  )
  continuous_fit <- .formula_coefficient_density_fit(
    continuous_result,
    "mu_x"
  )
  continuous_transform <- JAGS_formula_coefficient_transform(
    continuous_fit,
    "mu"
  )
  expect_identical(
    continuous_transform$targets$structural_status,
    c("dependent", "dependent")
  )

  continuous_density <- JAGS_formula_prior_density(
    continuous_fit,
    parameter = "mu",
    target = "mu_intercept"
  )
  continuous_ordinate <- prior_density_ordinate(continuous_density, 3)
  expect_identical(continuous_ordinate$behavior, "regular")
  expect_identical(continuous_ordinate$method, "scalar_affine")
  expect_true(continuous_ordinate$exact)
  expect_equal(
    continuous_ordinate$log_density,
    stats::dnorm(3, mean = 3, sd = 2, log = TRUE),
    tolerance = 1e-12
  )

  fixed_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = prior("point", list(3)),
      x = prior("point", list(0.5))
    ),
    formula_scale = TRUE
  )
  fixed_fit <- .formula_coefficient_density_fit(
    fixed_result,
    character()
  )
  fixed_transform <- JAGS_formula_coefficient_transform(fixed_fit, "mu")
  expect_identical(
    fixed_transform$targets$structural_status,
    c("structural", "structural")
  )
  expect_identical(fixed_transform$targets$fixed_value, c(2, 0.25))

  fixed_density <- JAGS_formula_prior_density(
    fixed_fit,
    parameter = "mu",
    target = "mu_intercept"
  )
  fixed_ordinate <- prior_density_ordinate(fixed_density, 2)
  expect_identical(fixed_ordinate$behavior, "point_mass")
  expect_identical(fixed_ordinate$point_mass, 1)
  expect_true(fixed_ordinate$exact)
})

test_that("log-intercept formula densities apply source and output Jacobians", {

  log_formula <- ~ 1 + x
  attr(log_formula, "log(intercept)") <- TRUE
  formula_result <- JAGS_formula(
    formula = log_formula,
    parameter = "mu",
    data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = prior("point", list(2)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = TRUE
  )
  fit <- .formula_coefficient_density_fit(formula_result, "mu_x")
  transform <- JAGS_formula_coefficient_transform(fit, "mu")

  expect_identical(
    transform$source_transforms,
    c(mu_intercept = "log", mu_x = "identity")
  )
  expect_identical(
    transform$output_transforms,
    c(mu_intercept = "exp", mu_x = "identity")
  )

  density <- JAGS_formula_prior_density(
    fit,
    parameter = "mu",
    target = "mu_intercept"
  )
  ordinate <- prior_density_ordinate(density, 1)
  expect_identical(ordinate$behavior, "regular")
  expect_identical(ordinate$method, "named_transform")
  expect_true(ordinate$exact)
  expect_equal(
    ordinate$log_density,
    stats::dlnorm(1, meanlog = log(2), sdlog = 2, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("unscaled log-intercepts retain their positive-scale transform", {

  log_formula <- ~ 1
  attr(log_formula, "log(intercept)") <- TRUE
  formula_result <- JAGS_formula(
    formula = log_formula,
    parameter = "log_tau",
    data = data.frame(row = seq_len(3L)),
    prior_list = list(
      intercept = prior("gamma", list(2, 2))
    ),
    formula_scale = TRUE
  )
  expect_null(formula_result$formula_scale)
  expect_true(formula_result$formula_design$log_intercept)

  fit <- .formula_coefficient_density_fit(
    formula_result,
    "log_tau_intercept"
  )
  transform <- JAGS_formula_coefficient_transform(fit, "log_tau")

  expect_identical(
    transform$source_transforms,
    c(log_tau_intercept = "log")
  )
  expect_identical(
    transform$output_transforms,
    c(log_tau_intercept = "exp")
  )
  expect_identical(
    transform$matrix,
    matrix(
      1,
      nrow = 1L,
      dimnames = list("log_tau_intercept", "log_tau_intercept")
    )
  )

  positive_draws <- matrix(
    c(0.25, 1, 4),
    ncol = 1L,
    dimnames = list(NULL, "log_tau_intercept")
  )
  expect_equal(
    .bt_apply_formula_coefficient_transform(positive_draws, transform),
    positive_draws,
    tolerance = 1e-14
  )

  density <- JAGS_formula_prior_density(
    fit,
    parameter = "log_tau",
    target = "log_tau_intercept"
  )
  ordinate <- prior_density_ordinate(density, 1.5)
  expect_identical(ordinate$behavior, "regular")
  expect_identical(ordinate$method, "named_transform")
  expect_true(ordinate$exact)
  expect_equal(
    ordinate$log_density,
    stats::dgamma(1.5, shape = 2, rate = 2, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("formula prior densities preserve model-mixture atoms and fail closed", {

  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = prior("point", list(0)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = TRUE
  )
  fit <- .formula_coefficient_density_fit(formula_result, "mu_x")
  mixture_context <- .prior_density_model_mixture_context(
    prior_list = list(
      mu_intercept = list(
        prior("point", list(0), prior_weights = 1),
        prior("point", list(0), prior_weights = 3)
      ),
      mu_x = list(
        prior("point", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 3)
      )
    ),
    column_names = c("mu_intercept", "mu_x"),
    n_grid = 128
  )
  density <- JAGS_formula_prior_density(
    fit,
    parameter = "mu",
    target = "mu_intercept",
    context = mixture_context
  )
  ordinate <- prior_density_ordinate(density, 0)
  expect_identical(ordinate$behavior, "point_mass")
  expect_equal(ordinate$point_mass, 0.25)
  expect_equal(
    ordinate$log_density,
    log(0.75) + stats::dnorm(0, sd = 2, log = TRUE),
    tolerance = 1e-12
  )

  incomplete_context <- .prior_density_context(
    prior_list = list(mu_intercept = prior("point", list(0))),
    column_names = "mu_intercept"
  )
  missing_error <- tryCatch(
    JAGS_formula_prior_density(
      fit,
      parameter = "mu",
      target = "mu_intercept",
      context = incomplete_context
    ),
    error = identity
  )
  expect_s3_class(
    missing_error,
    "BayesTools_formula_prior_density_unavailable"
  )
  expect_identical(missing_error$reason, "missing_source_coordinates")
  expect_identical(missing_error$missing, "mu_x")

  unknown_error <- tryCatch(
    JAGS_formula_prior_density(fit, "mu", "unknown"),
    error = identity
  )
  expect_s3_class(
    unknown_error,
    "BayesTools_formula_prior_density_unavailable"
  )
  expect_identical(unknown_error$reason, "unknown_target")

  stale <- JAGS_formula_coefficient_transform(fit, "mu")
  stale$schema_version <- 0L
  expect_error(
    .bt_validate_formula_coefficient_transform(stale),
    "missing or unsupported"
  )
})
