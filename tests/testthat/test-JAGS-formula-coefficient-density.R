skip_if_not_test_profile("unit")

.formula_coefficient_density_fit <- function(formula_result,
                                             sampled_columns){

  design <- formula_result$formula_design
  parameter <- design$parameter
  formula_design <- stats::setNames(list(design), parameter)
  registry <- .bt_build_parameter_registry(
    columns = sampled_columns,
    prior_list = formula_result$prior_list,
    formula_design = formula_design
  )

  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  if("formula_scale" %in% names(formula_result)){
    attr(fit, "formula_scale") <- stats::setNames(
      list(formula_result$formula_scale),
      parameter
    )
  }
  attr(fit, "formula_design") <- formula_design
  attr(fit, "parameter_registry") <- registry
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
  expect_identical(transform$formula_design_version, 3L)
  expect_identical(transform$parameter_registry_version, 3L)
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
