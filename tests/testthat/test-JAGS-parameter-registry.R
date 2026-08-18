skip_if_not_test_profile("unit")

test_that("parameter registry schema is explicit and versioned", {

  schema <- JAGS_parameter_registry_schema()
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

test_that("fitted registry classifies concrete random coordinates exactly", {

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
  registry <- build_test_parameter_registry(
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

  expect_s3_class(registry, "BayesTools_parameter_registry")
  expect_identical(attr(registry, "schema_version"), 4L)
  expect_identical(registry$coordinate_name, columns)
  expect_identical(anyDuplicated(registry$coordinate_name), 0L)
  expect_identical(
    registry$role,
    c(
      "fixed_coefficient",
      "random_sd",
      "random_latent",
      "random_group_coefficient"
    )
  )
  expect_identical(
    registry$fitted_scale,
    c(
      "fitted_standardized",
      "fitted_covariance",
      "unit_latent",
      "fitted_standardized"
    )
  )
  expect_identical(registry$internal, c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(registry$random_block, c("", "id", "id", "id"))
  expect_identical(registry$column, c("mu_intercept", "x", "x", "x"))
  expect_identical(
    registry$display_label,
    c(
      "(mu) intercept",
      "(mu) id: sd(x)",
      "(mu) z(id[a], x)",
      "(mu) coef(id[a], x)"
    )
  )
})

test_that("registry keeps LKJ primitive coordinates internal to their random block", {

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
  coordinates <- c(
    random_term$correlation$primitive_names,
    random_term$correlation$cpc_names
  )
  registry <- build_test_parameter_registry(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )

  expect_identical(
    random_term$sd_leaves$leaf_terms_by_column,
    c("group[sensitivity]", "group[specificity]")
  )
  expect_identical(random_term$contrast_owner, "random_block")
  expect_identical(registry$coordinate_name, coordinates)
  expect_identical(
    registry$role,
    rep("random_correlation_coordinate", length(coordinates))
  )
  expect_identical(registry$formula_parameter, rep("mu", length(coordinates)))
  expect_identical(registry$random_block, rep("study", length(coordinates)))
  expect_identical(registry$random_name, rep("study", length(coordinates)))
  expect_identical(registry$random_grouping, rep("study", length(coordinates)))
  expect_identical(registry$random_structure, rep("us", length(coordinates)))
  expect_identical(registry$fitted_scale, rep("unitless", length(coordinates)))
  expect_true(all(registry$internal))
})

test_that("registry owns random SD spike-and-slab auxiliaries", {

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
  sd_base <- unique(BayesTools:::.bt_parameter_registry_base(sd_names))
  expect_length(sd_base, 1L)

  columns <- c(
    paste0(sd_base, "_indicator"),
    paste0(sd_base, "_inclusion"),
    sd_names,
    paste0(sd_base, "_variable[", seq_along(sd_names), "]")
  )
  registry <- build_test_parameter_registry(
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
  registry <- registry[
    match(columns, registry$coordinate_name),
    ,
    drop = FALSE
  ]

  expect_identical(
    registry$role,
    c(
      "random_inclusion_indicator",
      "random_inclusion_probability",
      "random_sd",
      "random_sd",
      "random_sd_variable",
      "random_sd_variable"
    )
  )
  expect_identical(registry$random_block, rep("id", length(columns)))
  expect_identical(registry$formula_parameter, rep("mu", length(columns)))
  expect_identical(
    registry$fitted_scale,
    c(
      "unitless",
      "unitless",
      "fitted_covariance",
      "fitted_covariance",
      "fitted_covariance",
      "fitted_covariance"
    )
  )
  expect_identical(
    registry$column,
    c("", "", "x_fac3B", "x_fac3C", "x_fac3B", "x_fac3C")
  )
  expect_identical(
    registry$internal,
    c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
  )
})

test_that("registry display labels do not overwrite fixed scale formatting", {

  interaction_prior <- prior("normal", list(0, 1))
  attr(interaction_prior, "parameter") <- "mu"
  log_intercept_prior <- prior("normal", list(0, 1))
  attr(log_intercept_prior, "parameter") <- "log_sigma"
  prior_list <- list(
    mu_x__xXx__z = interaction_prior,
    log_sigma_intercept = log_intercept_prior
  )
  coordinate_names <- names(prior_list)
  registry <- build_test_parameter_registry(
    columns = coordinate_names,
    prior_list = prior_list
  )

  expect_identical(
    registry$display_label,
    c("(mu) x:z", "(log_sigma) intercept")
  )

  formatted_names <- c("(mu) x:z", "(log_sigma) exp(intercept)")
  expect_identical(
    BayesTools:::.bt_random_effect_summary_display_names(
      names = formatted_names,
      raw_names = coordinate_names,
      prior_list = prior_list,
      parameter_registry = registry
    ),
    formatted_names
  )
})

test_that("registry accessor rejects unversioned and malformed fitted objects", {

  samples <- coda::mcmc(matrix(
    1:4,
    ncol = 1L,
    dimnames = list(NULL, "theta")
  ))
  fit <- samples
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(theta = prior("normal", list(0, 1)))

  expect_error(
    JAGS_parameter_registry(fit),
    "Refit the model"
  )

  fit <- attach_test_parameter_registry(fit)
  registry <- JAGS_parameter_registry(fit)
  expect_identical(registry$coordinate_name, "theta")
  expect_identical(registry$role, "parameter")

  malformed <- registry
  malformed$coordinate_name <- ""
  attr(fit, "parameter_registry") <- malformed
  expect_error(
    JAGS_parameter_registry(fit),
    "unique, non-missing coordinate names"
  )
})

test_that("structural point parameters are registered when JAGS omits them", {

  registry <- build_test_parameter_registry(
    columns = "theta",
    prior_list = list(
      theta = prior("normal", list(0, 1)),
      fixed = prior("point", list(0))
    )
  )

  fixed <- registry[registry$coordinate_name == "fixed", , drop = FALSE]
  expect_equal(nrow(fixed), 1L)
  expect_identical(fixed$monitor_status, "structural")
  expect_identical(fixed$fixed_value, 0)
  expect_identical(fixed$role, "parameter")
  expect_false(fixed$internal)

  monitored <- build_test_parameter_registry(
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

test_that("structural registry coordinates retain exact scalar and vector values", {

  factor_prior <- prior_factor(
    "point",
    list(location = -2),
    contrast = "treatment"
  )
  attr(factor_prior, "levels") <- 3L
  registry <- build_test_parameter_registry(
    columns = "theta",
    prior_list = list(
      theta = prior("normal", list(0, 1)),
      scalar = prior("point", list(3.5)),
      vector = prior("mpoint", list(location = 2, K = 3)),
      factor = factor_prior
    )
  )

  structural <- registry[registry$monitor_status == "structural", ]
  expect_identical(
    structural$coordinate_name,
    c("scalar", "vector[1]", "vector[2]", "vector[3]", "factor[1]", "factor[2]")
  )
  expect_identical(structural$fixed_value, c(3.5, 2, 2, 2, -2, -2))
  expect_true(all(is.na(registry$fixed_value[registry$monitor_status == "sampled"])))
})

test_that("registry prevents prefix-related random-block ownership collisions", {

  short_term <- list(
    parameter_stem = "mu__xREx__a",
    parameter = "mu",
    block_name = "a",
    group_label = "group_a",
    has_explicit_name = TRUE,
    structure = "diag",
    column_names = "intercept",
    sd_parameter_names = "mu__xREx__a_intercept",
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
  registry <- build_test_parameter_registry(
    columns = columns,
    formula_design = formula_design
  )

  expect_identical(registry$random_block, c("", "a", "a_b"))

  samples <- matrix(
    seq_len(6L),
    nrow = 2L,
    dimnames = list(NULL, columns)
  )
  removed <- BayesTools:::.bt_random_effect_summary_filter_raw_columns(
    model_samples = samples,
    parameter_registry = registry,
    remove_random_effects = "a"
  )
  expect_identical(
    colnames(removed),
    c("mu_intercept", "mu__xREx__a_b_intercept")
  )
})
