skip_if_not_test_profile("unit")

test_that("parameter registry schema is explicit and versioned", {

  schema <- JAGS_parameter_registry_schema()
  expect_identical(
    schema$field,
    c(
      "canonical_name", "monitor_name", "formula_parameter", "role",
      "random_block", "random_name", "term", "column", "index", "dimensions",
      "fitted_scale", "monitor_status", "display_label",
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
  expect_identical(attr(registry, "schema_version"), 1L)
  expect_identical(registry$canonical_name, columns)
  expect_identical(anyDuplicated(registry$canonical_name), 0L)
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
      "(mu) sd(x | id)",
      "(mu) z(id[a], x)",
      "(mu) coef(id[a], x)"
    )
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
  canonical_names <- names(prior_list)
  registry <- build_test_parameter_registry(
    columns = canonical_names,
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
      raw_names = canonical_names,
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
  expect_identical(registry$canonical_name, "theta")
  expect_identical(registry$role, "parameter")

  malformed <- registry
  malformed$canonical_name <- ""
  attr(fit, "parameter_registry") <- malformed
  expect_error(
    JAGS_parameter_registry(fit),
    "unique, non-missing canonical names"
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

  fixed <- registry[registry$canonical_name == "fixed", , drop = FALSE]
  expect_equal(nrow(fixed), 1L)
  expect_identical(fixed$monitor_status, "structural")
  expect_identical(fixed$role, "parameter")
  expect_false(fixed$internal)
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
