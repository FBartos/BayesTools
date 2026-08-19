skip_if_not_test_profile("unit")

test_that("parameter map and coordinate schemas are explicit and versioned", {

  map_schema <- parameter_map_schema()
  schema <- parameter_coordinates_schema()
  expect_identical(map_schema$schema_version, 1L)
  expect_identical(
    schema$field,
    c(
      "coordinate_name", "monitor_name", "formula_parameter", "role",
      "random_block", "random_name", "term", "column", "index", "dimensions",
      "fitted_scale", "monitor_status", "fixed_value", "display_label",
      "random_grouping", "random_structure", "internal"
    )
  )
  expect_identical(schema$type[nrow(schema)], "logical")
  expect_identical(anyDuplicated(schema$field), 0L)
})

test_that("fitted coordinates classify concrete random coordinates exactly", {

  data <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = TRUE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    )
  )
  columns <- c(
    "mu_intercept",
    "mu__xREx__id_x",
    "mu__xREx__id_xRE_Zx[1,1]",
    "mu__xREx__id_xRE_COEFx[1,1]"
  )
  coordinates <- build_test_parameter_coordinates(
    columns = columns,
    monitor_names = c(
      "mu_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_Zx",
      "mu__xREx__id_xRE_COEFx"
    ),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design),
    formula_scale = list(mu = formula_result$formula_scale)
  )

  expect_s3_class(coordinates, "BayesTools_parameter_coordinates")
  expect_identical(coordinates$coordinate_name, columns)
  expect_identical(anyDuplicated(coordinates$coordinate_name), 0L)
  expect_identical(
    coordinates$role,
    c(
      "fixed_coefficient",
      "random_sd",
      "random_latent",
      "random_group_coefficient"
    )
  )
  expect_identical(
    coordinates$fitted_scale,
    c(
      "fitted_standardized",
      "fitted_covariance",
      "unit_latent",
      "fitted_standardized"
    )
  )
  expect_identical(coordinates$internal, c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(coordinates$random_block, c("", "id", "id", "id"))
  expect_identical(coordinates$column, c("mu_intercept", "x", "x", "x"))
  expect_identical(
    coordinates$display_label,
    c(
      "(mu) intercept",
      "(mu) id: sd(x)",
      "(mu) z(id[a], x)",
      "(mu) coef(id[a], x)"
    )
  )
})

test_that("coordinate map keeps LKJ primitives internal to their random block", {

  data <- data.frame(
    group = factor(
      c("sensitivity", "specificity", "sensitivity", "specificity"),
      levels = c("sensitivity", "specificity")
    ),
    study = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(0 + group | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(
        sd = prior("gamma", list(2, 2)),
        covariance = random_covariance(cor = prior_lkj(eta = 1)),
        monitor = random_monitor(lkj_primitives = TRUE),
        contrasts = c(group = "independent")
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  coordinate_names <- c(
    random_term$correlation$primitive_names,
    random_term$correlation$cpc_names
  )
  coordinates <- build_test_parameter_coordinates(
    columns = coordinate_names,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )

  expect_identical(
    random_term$sd_leaves$leaf_terms_by_column,
    c("group[sensitivity]", "group[specificity]")
  )
  expect_identical(random_term$contrast_owner, "random_block")
  expect_identical(coordinates$coordinate_name, coordinate_names)
  expect_identical(
    coordinates$role,
    rep("random_correlation_coordinate", length(coordinate_names))
  )
  expect_identical(coordinates$formula_parameter, rep("mu", length(coordinate_names)))
  expect_identical(coordinates$random_block, rep("study", length(coordinate_names)))
  expect_identical(coordinates$random_name, rep("study", length(coordinate_names)))
  expect_identical(coordinates$random_grouping, rep("study", length(coordinate_names)))
  expect_identical(coordinates$random_structure, rep("us", length(coordinate_names)))
  expect_identical(coordinates$fitted_scale, rep("unitless", length(coordinate_names)))
  expect_true(all(coordinates$internal))
})

test_that("Dirichlet auxiliary coordinates remain coordinate-only", {

  prior_list <- list(
    weights = prior("dirichlet", list(alpha = c(2, 3)))
  )
  columns <- c(
    "weights[1]",
    "weights[2]",
    "prior_par_eta_weights[1]",
    "prior_par_eta_weights[2]"
  )
  coordinates <- build_test_parameter_coordinates(
    columns = columns,
    prior_list = prior_list
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates = coordinates,
    prior_list = prior_list
  )

  expect_identical(
    coordinates$internal,
    c(FALSE, FALSE, TRUE, TRUE)
  )
  expect_identical(
    catalog$quantities$canonical_name,
    c("weights[1]", "weights[2]")
  )
})

test_that("formula metadata exposes exact LKJ primitive coordinate priors", {

  data <- data.frame(
    group = factor(rep(c("a", "b", "c"), 2L)),
    study = factor(rep(c("s1", "s2"), each = 3L))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(0 + group | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(
        sd = prior("gamma", list(2, 2)),
        covariance = random_covariance(cor = prior_lkj(eta = 2)),
        contrasts = c(group = "independent")
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  columns <- c(
    "mu_intercept",
    random_term$sd_parameter_names,
    random_term$correlation$primitive_names
  )
  posterior <- matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = columns,
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit, monitor_names = columns)

  coordinate_priors <- JAGS_formula_internal_coordinate_priors(fit)
  alpha <- vapply(
    coordinate_priors,
    function(x) x$parameters$alpha,
    numeric(1)
  )

  expect_identical(
    names(coordinate_priors),
    random_term$correlation$primitive_names
  )
  expect_equal(unname(alpha), c(2.5, 2.5, 2))
  expect_true(all(vapply(coordinate_priors, is.prior, logical(1))))
  expect_true(all(vapply(coordinate_priors, function(x) {
    identical(x$distribution, "beta") &&
      identical(x$parameters$alpha, x$parameters$beta)
  }, logical(1))))
})

test_that("coordinate map owns random SD spike-and-slab auxiliaries", {

  data <- data.frame(
    x_fac3 = factor(
      c("A", "B", "C", "A", "B", "C"),
      levels = c("A", "B", "C")
    ),
    id = factor(c("one", "one", "one", "two", "two", "two"))
  )
  formula_result <- JAGS_formula(
    formula = ~ -1 + x_fac3 + (x_fac3 - 1 || id),
    parameter = "mu",
    data = data,
    prior_list = list(
      x_fac3 = prior_factor(
        "normal",
        list(0, 1),
        contrast = "independent"
      )
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior_spike_and_slab(
          prior("normal", list(0, 1), list(0, 1))
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  sd_names <- random_term$sd_parameter_names
  sd_columns <- random_term$sd_leaves$column_names
  n_sd <- length(sd_names)
  expect_identical(
    random_term$sd_leaves$leaf_terms_by_column,
    paste0("x_fac3[", levels(data$x_fac3), "]")
  )
  sd_base <- unique(BayesTools:::.bt_parameter_coordinates_base(sd_names))
  expect_length(sd_base, 1L)

  columns <- c(
    paste0(sd_base, "_indicator"),
    paste0(sd_base, "_inclusion"),
    sd_names,
    paste0(sd_base, "_variable[", seq_along(sd_names), "]")
  )
  coordinates <- build_test_parameter_coordinates(
    columns = columns,
    monitor_names = c(
      paste0(sd_base, "_indicator"),
      paste0(sd_base, "_inclusion"),
      sd_base,
      paste0(sd_base, "_variable")
    ),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  coordinates <- coordinates[
    match(columns, coordinates$coordinate_name),
    ,
    drop = FALSE
  ]

  expect_identical(
    coordinates$role,
    c(
      "random_inclusion_indicator",
      "random_inclusion_probability",
      rep("random_sd", n_sd),
      rep("random_sd_variable", n_sd)
    )
  )
  expect_identical(coordinates$random_block, rep("id", length(columns)))
  expect_identical(coordinates$formula_parameter, rep("mu", length(columns)))
  expect_identical(
    coordinates$fitted_scale,
    c(
      "unitless",
      "unitless",
      rep("fitted_covariance", 2L * n_sd)
    )
  )
  expect_identical(
    coordinates$column,
    c("", "", sd_columns, sd_columns)
  )
  expect_identical(
    coordinates$internal,
    c(TRUE, TRUE, rep(FALSE, n_sd), rep(TRUE, n_sd))
  )
})

test_that("coordinate labels do not overwrite fixed scale formatting", {

  interaction_prior <- prior("normal", list(0, 1))
  attr(interaction_prior, "parameter") <- "mu"
  log_intercept_prior <- prior("normal", list(0, 1))
  attr(log_intercept_prior, "parameter") <- "log_sigma"
  prior_list <- list(
    mu_x__xXx__z = interaction_prior,
    log_sigma_intercept = log_intercept_prior
  )
  coordinate_names <- names(prior_list)
  coordinates <- build_test_parameter_coordinates(
    columns = coordinate_names,
    prior_list = prior_list
  )

  expect_identical(
    coordinates$display_label,
    c("(mu) x:z", "(log_sigma) intercept")
  )

  formatted_names <- c("(mu) x:z", "(log_sigma) exp(intercept)")
  expect_identical(
    BayesTools:::.bt_random_effect_summary_display_names(
      names = formatted_names,
      raw_names = coordinate_names,
      prior_list = prior_list,
      coordinates = coordinates
    ),
    formatted_names
  )
})

test_that("parameter-map accessors reject missing and malformed fitted objects", {

  samples <- coda::mcmc(matrix(
    1:4,
    ncol = 1L,
    dimnames = list(NULL, "theta")
  ))
  fit <- samples
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(theta = prior("normal", list(0, 1)))

  expect_error(
    parameter_map(fit),
    "Refit the model"
  )

  fit <- attach_test_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)
  map <- parameter_map(fit)
  coordinates <- parameter_coordinates(fit)
  catalog <- parameter_catalog(fit)
  expect_null(attr(fit, "parameter_registry", exact = TRUE))
  expect_null(attr(fit, "parameter_catalog", exact = TRUE))
  expect_identical(names(map), c(
    "schema_version", "coordinates", "quantities", "aliases"
  ))
  expect_identical(map$coordinates, coordinates)
  expect_identical(map$quantities, catalog$quantities)
  expect_identical(map$aliases, catalog$aliases)
  expect_identical(JAGS_fit_contract(fit)$parameter_map_version,
                   map$schema_version)
  expect_identical(coordinates$coordinate_name, "theta")
  expect_identical(coordinates$role, "parameter")

  malformed <- coordinates
  malformed$coordinate_name <- ""
  map$coordinates <- malformed
  attr(fit, "parameter_map") <- map
  expect_error(
    parameter_coordinates(fit),
    "unique, non-missing coordinate names"
  )
})

test_that("parameter map validates semantic dependencies atomically", {

  map <- .bt_build_parameter_map(columns = "theta")
  broken <- map
  broken$quantities$extraction_key[[1L]] <- list(
    type = "factor_level",
    dependencies = "missing_coordinate",
    weights = 1
  )

  expect_error(
    .bt_validate_parameter_map(broken),
    "unknown coordinate dependencies: 'missing_coordinate'"
  )
})

test_that("structural point parameters are registered when JAGS omits them", {

  coordinates <- build_test_parameter_coordinates(
    columns = "theta",
    prior_list = list(
      theta = prior("normal", list(0, 1)),
      fixed = prior("point", list(0))
    )
  )

  fixed <- coordinates[coordinates$coordinate_name == "fixed", , drop = FALSE]
  expect_equal(nrow(fixed), 1L)
  expect_identical(fixed$monitor_status, "structural")
  expect_identical(fixed$fixed_value, 0)
  expect_identical(fixed$role, "parameter")
  expect_false(fixed$internal)

  monitored <- build_test_parameter_coordinates(
    columns = c("theta", "fixed"),
    monitor_names = c("theta", "fixed"),
    prior_list = list(
      theta = prior("normal", list(0, 1)),
      fixed = prior("point", list(0))
    )
  )
  fixed_monitored <- monitored[
    monitored$coordinate_name == "fixed",
    ,
    drop = FALSE
  ]
  expect_identical(fixed_monitored$monitor_status, "structural")
  expect_identical(fixed_monitored$fixed_value, 0)
})

test_that("structural coordinates retain exact scalar and vector values", {

  factor_prior <- prior_factor(
    "point",
    list(location = -2),
    contrast = "treatment"
  )
  attr(factor_prior, "levels") <- 3L
  coordinates <- build_test_parameter_coordinates(
    columns = "theta",
    prior_list = list(
      theta = prior("normal", list(0, 1)),
      scalar = prior("point", list(3.5)),
      vector = prior("mpoint", list(location = 2, K = 3)),
      factor = factor_prior
    )
  )

  structural <- coordinates[coordinates$monitor_status == "structural", ]
  expect_identical(
    structural$coordinate_name,
    c("scalar", "vector[1]", "vector[2]", "vector[3]", "factor[1]", "factor[2]")
  )
  expect_identical(structural$fixed_value, c(3.5, 2, 2, 2, -2, -2))
  expect_true(all(is.na(coordinates$fixed_value[coordinates$monitor_status == "sampled"])))
})

test_that("coordinate map prevents random-block ownership collisions", {

  short_term <- list(
    parameter_stem = "mu__xREx__a",
    parameter = "mu",
    block_name = "a",
    group_label = "group_a",
    has_explicit_name = TRUE,
    structure = "diag",
    column_names = "intercept",
    sd_parameter_names = "mu__xREx__a_intercept",
    sd_leaves = structure(
      list(leaf_terms = c(mu__xREx__a_intercept = "intercept")),
      class = c("BayesTools_random_effect_sd_leaves", "list")
    ),
    group_levels = "one"
  )
  long_term <- list(
    parameter_stem = "mu__xREx__a_b",
    parameter = "mu",
    block_name = "a_b",
    group_label = "group_a_b",
    has_explicit_name = TRUE,
    structure = "diag",
    column_names = "intercept",
    sd_parameter_names = "mu__xREx__a_b_intercept",
    sd_leaves = structure(
      list(leaf_terms = c(mu__xREx__a_b_intercept = "intercept")),
      class = c("BayesTools_random_effect_sd_leaves", "list")
    ),
    group_levels = "one"
  )
  formula_design <- list(
    mu = structure(
      list(
        parameter = "mu",
        random_effects = list(short_term, long_term)
      ),
      class = c("BayesTools_formula_design", "list")
    )
  )
  columns <- c(
    "mu_intercept",
    "mu__xREx__a_intercept",
    "mu__xREx__a_b_intercept"
  )
  coordinates <- build_test_parameter_coordinates(
    columns = columns,
    formula_design = formula_design
  )

  expect_identical(coordinates$random_block, c("", "a", "a_b"))

  samples <- matrix(
    seq_len(6L),
    nrow = 2L,
    dimnames = list(NULL, columns)
  )
  removed <- BayesTools:::.bt_random_effect_summary_filter_raw_columns(
    model_samples = samples,
    coordinates = coordinates,
    remove_random_effects = "a"
  )
  expect_identical(
    colnames(removed),
    c("mu_intercept", "mu__xREx__a_b_intercept")
  )
})
