skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Random-effect design, grouping and structured terms
# ============================================================================ #
#
# PURPOSE:
#   Random-effect formula parsing, grouping maps, covariance specials, known group
#   covariance, and structured (level-indexed) random-effect terms.
# ============================================================================ #

test_that("named covariance triangles follow supplied rows after column alignment", {

  lower <- matrix(c(2, .2, .3, 0, 3, .4, 0, 0, 4), 3L,
    dimnames = list(letters[1:3], letters[1:3]))
  complete <- lower + t(lower) - diag(diag(lower))
  for(partial in list(lower, t(lower))){
    for(columns in list(1:3, c(3L, 1L, 2L))){
      actual <- random_group_covariance(partial[, columns], scale = "none")
      expect_equal(actual$covariance, complete, ignore_attr = TRUE)
    }
  }
  expect_error(random_group_covariance(lower[c(3L, 1L, 2L), ]),
    "must be exactly symmetric", fixed = TRUE)
  for(scale in c("cor", "none", "cor0", "cov0")){
    original <- random_group_covariance(complete, scale = scale)
    permuted <- random_group_covariance(complete[c(3L, 1L, 2L), c(2L, 3L, 1L)], scale = scale)
    expected <- BayesTools:::.bt_prepare_group_covariance_kernel(original, letters[1:3], "g")
    actual <- BayesTools:::.bt_prepare_group_covariance_kernel(permuted, letters[1:3], "g")
    expect_equal(actual$kernel, expected$kernel, tolerance = 0)
  }
  sparse <- complete
  sparse[1L, 3L] <- sparse[3L, 1L] <- 0
  expect_equal(random_group_covariance(sparse[c(3L, 1L, 2L), c(2L, 3L, 1L)])$covariance,
    sparse[c(3L, 1L, 2L), c(3L, 1L, 2L)], ignore_attr = TRUE)
  asymmetric <- complete
  asymmetric[1L, 2L] <- asymmetric[1L, 2L] + .Machine$double.eps
  expect_error(random_group_covariance(asymmetric[, c(3L, 1L, 2L)]),
    "must be exactly symmetric", fixed = TRUE)
})

test_that("group ownership retains exact numeric identity and safe display replay", {

  compile <- function(data, formula = ~ 1 + diag(1 | g), group_covariance = NULL){
    if(!is.null(group_covariance)){
      formula <- random_effects_formula(~ diag(1 | g), group_covariance = group_covariance)
    }
    JAGS_formula(formula, "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        sd = prior("normal", list(0, 1), list(0, Inf))
      ))
  }
  adjacent <- c(.3, .3 + .Machine$double.eps)
  data <- data.frame(g = adjacent)
  result <- compile(data)
  term <- result$formula_design$random_effects[[1L]]
  expect_identical(term$n_groups, 2L)
  expect_identical(term$group_map, 1:2)
  expect_identical(term$group_component_display_values$g, c("0.3", "0.3"))
  expect_false(anyDuplicated(term$group_levels) > 0L)
  replay <- function(term, data, allow_new = FALSE){
    BayesTools:::.bt_random_effect_prediction_data(term, data,
      group_data = data, allow_new_groups = allow_new)
  }
  expect_identical(replay(term, data[2L, , drop = FALSE])$group_map, 2L)
  expect_identical(replay(term, data[2:1, , drop = FALSE])$group_map, 2:1)
  expect_error(replay(term, data.frame(g = "0.3"), TRUE), "is ambiguous", fixed = TRUE)
  expect_identical(replay(term, data.frame(g = term$group_component_levels$g[2L]))$group_map, 2L)
  for(allow_new in c(FALSE, TRUE)){
    for(group_label in list(term$group_levels[2L], factor(term$group_levels[2L]))){
      public <- replay(term, data.frame(g = group_label), allow_new)
      expect_identical(public$group_map, 2L)
      expect_identical(public$group_levels, term$group_levels)
      expect_identical(public$group_tuple_keys, term$group_tuple_keys)
    }
  }
  expect_identical(replay(term, data.frame(g = .3 + 2 * .Machine$double.eps), TRUE)$group_map, 3L)

  ordinary <- compile(data.frame(g = c(2, 1, -0, 0)))$formula_design$random_effects[[1L]]
  expect_identical(ordinary$group_levels, c("0", "1", "2"))
  expect_identical(replay(ordinary, data.frame(g = c(0L, 1L, 2L)))$group_map, 1:3)
  expect_identical(replay(ordinary, data.frame(g = c("0", "1", "2")))$group_map, 1:3)
  character_owned <- compile(data.frame(g = c("1", "2")))$formula_design$random_effects[[1L]]
  expect_identical(replay(character_owned, data.frame(g = c(2, 1)))$group_map, 2:1)
  factor_owned <- compile(data.frame(g = factor(c("é", "a"), levels = c("é", "a", "unused"))))$formula_design$random_effects[[1L]]
  expect_identical(factor_owned$group_levels, c("é", "a", "unused"))
  expect_identical(factor_owned$group_observed_levels, c("é", "a"))
  expect_identical(replay(factor_owned, data.frame(g = "é"))$group_map, 1L)
  infinite <- compile(data.frame(g = c(Inf, -Inf)))$formula_design$random_effects[[1L]]
  expect_identical(replay(infinite, data.frame(g = c(-Inf, Inf)))$group_map, 1:2)

  nested_data <- expand.grid(g = adjacent, h = adjacent)
  nested <- compile(nested_data, ~ 1 + diag(1 | g / h))$formula_design$random_effects[[1L]]
  expect_identical(nested$n_groups, 4L)
  expect_identical(replay(nested, nested_data[4:1, ])$group_map, nested$group_map[4:1])
  old <- options(OutDec = ",")
  on.exit(options(old), add = TRUE)
  expect_identical(replay(term, data)$group_map, 1:2)
  expect_identical(BayesTools:::.bt_random_effect_numeric_identity(adjacent),
    c("0.29999999999999999", "0.30000000000000021"))
  expect_identical(BayesTools:::.bt_random_effect_structured_index_component(adjacent)$level_keys,
    substring(term$group_component_identity_keys$g, 3L))
  options(old)

  kernel <- diag(2)
  dimnames(kernel) <- list(term$group_levels, term$group_levels)
  known <- compile(data, group_covariance = random_group_covariance(kernel))
  expect_equal(known$formula_design$random_effects[[1L]]$group_covariance$kernel, kernel)
  known_term <- known$formula_design$random_effects[[1L]]
  for(allow_new in c(FALSE, TRUE)){
    for(group_label in list(known_term$group_levels[2L], factor(known_term$group_levels[2L]))){
      public <- replay(known_term, data.frame(g = group_label), allow_new)
      expect_identical(public$group_map, 2L)
      expect_identical(public$group_levels, known_term$group_levels)
      expect_identical(public$group_tuple_keys, known_term$group_tuple_keys)
    }
  }
  rounded <- matrix(1, 1L, 1L, dimnames = list("0.3", "0.3"))
  expect_error(compile(data, group_covariance = random_group_covariance(rounded)),
    "is missing fitted level(s)", fixed = TRUE)

  posterior <- matrix(1, 1L, 1L, dimnames = list(NULL, unique(term$sd_parameter_names)))
  legacy <- result$formula_design
  legacy$schema_version <- 4L
  expect_error(random_effects_marginal_vcov(legacy, posterior_samples = posterior,
    prior_list = result$prior_list, new_levels = "zero"), class = "BayesTools_refit_required")
  expect_error(random_effects_marginal_vcov(list(mu = legacy), posterior_samples = posterior,
    prior_list = result$prior_list), class = "BayesTools_refit_required")
  missing <- term
  missing$group_component_identity_keys <- NULL
  expect_error(replay(missing, data.frame(g = "unknown"), TRUE), class = "BayesTools_refit_required")
  malformed <- term
  malformed$group_component_identity_keys$g <- c("n:invalid", "n:invalid2")
  expect_error(replay(malformed, data, TRUE), class = "BayesTools_refit_required")
  malformed$group_component_identity_keys$g <- 1:2
  expect_error(replay(malformed, data, TRUE), class = "BayesTools_refit_required")
  malformed$group_component_identity_keys$g <- list("n:0.3", "n:0.4")
  expect_error(replay(malformed, data, TRUE), class = "BayesTools_refit_required")
})

test_that("AR1 and HAR warn without changing misleading numeric label order", {
  dat <- data.frame(g = factor(rep(c("a", "b"), each = 3)))
  compile <- function(wave, structure){
    dat$wave <- wave
    JAGS_formula(
      stats::as.formula(paste0("~ ", structure, "(wave | g)")), "mu", dat,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("point", list(location = 1)))
    )$formula_design$random_effects[[1]]
  }
  for(structure in c("ar1", "har")){
    for(wave in list(rep(c("1", "2", "10"), 2), factor(rep(c("1", "2", "10"), 2)))){
      expect_warning(block <- compile(wave, structure), "order: 1, 10, 2", fixed = TRUE)
      expect_identical(block$column_names, c("wave1", "wave10", "wave2"))
    }
    expect_no_warning(block <- compile(rep(c(1, 2, 10), 2), structure))
    expect_identical(block$column_names, c("wave1", "wave2", "wave10"))
    expect_no_warning(block <- compile(ordered(rep(c("1", "2", "10"), 2),
      levels = c("1", "10", "2")), structure))
    expect_identical(block$column_names, c("wave1", "wave10", "wave2"))
  }
})

test_that("random name maps distinguish prefix-related block names", {
  dat <- data.frame(g = factor(c("a", "b", "a", "b")))
  result <- JAGS_formula(
    ~ random(1 | g, name = "study", covariance = "diag") +
      random(1 | g, name = "study_long", covariance = "diag") +
      random(1 | g, name = "study2", covariance = "diag"),
    "mu", dat, list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("point", list(location = 1)))
  )
  map <- result$formula_design$name_map
  for(block in result$formula_design$random_effects){
    row <- match(block$sd_parameter_names, map$jags_name)
    expect_false(anyNA(row))
    expect_true(all(map$kind[row] == "random"))
    expect_true(all(map$term[row] == block$block_name))
  }
})

test_that("CAR single-coordinate covariance and observed factor basis replay", {
  dat <- data.frame(g = factor(rep(c("a", "b"), each = 3)), time = 0)
  compile <- function(data){
    JAGS_formula(~ car(time | g), "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("point", list(location = 1))))
  }
  result <- compile(dat)
  block <- result$formula_design$random_effects[[1]]
  expect_null(block$correlation)
  expect_equal(block$car$time_values, 0)
  covariance <- JAGS_formula_random_marginal_covariance(
    result$formula_design, list(seq_len(nrow(dat))))
  expect_match(covariance$syntax, "_cor[1,1] = 1", fixed = TRUE)
  expect_equal(as.numeric(random_effects_correlation_draws(block, matrix(numeric(), 2, 0))), c(1, 1))

  dat$time <- ordered(rep(c("0", "1", "5"), 2), levels = c("0", "1", "5", "20"))
  result <- compile(dat)
  block <- result$formula_design$random_effects[[1]]
  expect_identical(block$xlevels$time, c("0", "1", "5"))
  expect_equal(unname(block$contrast_matrices$time), diag(3))
  replay <- .bt_random_effect_prediction_data(block, dat)
  expect_equal(replay$model_matrix, block$model_matrix)
})

test_that("structured index levels are identified by exact tuple and numeric keys", {
  compile <- function(formula, data){
    JAGS_formula(formula, "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("point", list(location = 1)))
    )$formula_design$random_effects[[1]]
  }

  # Character cells whose '.'-pasted labels coincide stay distinct.
  character_data <- data.frame(
    g = factor(rep(c("s1", "s2", "s3"), each = 2)),
    a = rep(c("x.y", "x"), 3),
    b = rep(c("z", "y.z"), 3)
  )
  block <- compile(~ 1 + cs(a + b | g), character_data)
  expect_identical(block$n_columns, 2L)
  expect_identical(block$column_names, c("a_bx_y.z", "a_bx.y_z"))
  expect_equal(
    block$model_matrix,
    cbind(rep(c(0, 1), 3), rep(c(1, 0), 3)),
    ignore_attr = TRUE
  )

  # Numeric cells (1.5, 2) and (1, 5.2) share the label "1.5.2".
  numeric_data <- data.frame(
    g = factor(rep(c("s1", "s2"), each = 3)),
    a = rep(c(1.5, 1, 2), 2),
    b = rep(c(2, 5.2, 2), 2)
  )
  block <- compile(~ 1 + cs(a + b | g), numeric_data)
  expect_identical(block$n_columns, 3L)
  expect_identical(block$column_names, c("a_b1_5.2", "a_b1.5_2", "a_b2_2"))
  expected <- cbind(
    numeric_data$a == 1,
    numeric_data$a == 1.5,
    numeric_data$a == 2
  ) * 1
  expect_equal(block$model_matrix, expected, ignore_attr = TRUE)
  replay <- .bt_random_effect_prediction_data(block, numeric_data)
  expect_equal(replay$model_matrix, block$model_matrix, ignore_attr = TRUE)
  new_cell <- .bt_random_effect_prediction_data(
    block,
    data.frame(g = factor("s1"), a = 1, b = 5.2)
  )
  expect_equal(
    new_cell$model_matrix,
    matrix(c(1, 0, 0), nrow = 1L),
    ignore_attr = TRUE
  )
  expect_error(
    .bt_random_effect_prediction_data(
      block,
      data.frame(g = factor("s1"), a = 5, b = 5)
    ),
    "do not match the levels used for model specification"
  )

  # Distinct doubles with equal 15-digit labels never merge or crash.
  float_data <- data.frame(
    g = factor(rep(c("s1", "s2"), each = 2)),
    t = c(0.3, 0.1 + 0.2, 0.3, 0.1 + 0.2)
  )
  for(structure in c("cs", "ar1")){
    block <- compile(
      stats::as.formula(paste0("~ 1 + ", structure, "(t | g)")),
      float_data
    )
    expect_identical(
      block$column_names,
      c("t0.29999999999999999", "t0.30000000000000004"),
      info = structure
    )
    expect_equal(
      block$model_matrix,
      cbind(c(1, 0, 1, 0), c(0, 1, 0, 1)),
      ignore_attr = TRUE,
      info = structure
    )
  }
  expect_error(
    .bt_random_effect_prediction_data(
      block,
      data.frame(g = factor("s1"), t = 0.3 + 1e-12)
    ),
    "do not match the levels used for model specification"
  )

  # Unambiguous labels are unchanged.
  block <- compile(
    ~ 1 + cs(a + b | g),
    data.frame(g = factor(rep(c("s1", "s2"), each = 2)), a = c(1, 2, 1, 2),
               b = c(0.5, 0.5, 3, 3))
  )
  expect_identical(
    block$column_names,
    c("a_b1.0.5", "a_b1.3", "a_b2.0.5", "a_b2.3")
  )
  # Index metadata without stored level keys ask for a refit.
  legacy <- block
  legacy$structured_index$level_keys <- NULL
  legacy$structured_index$levels <- NULL
  expect_error(
    .bt_random_effect_prediction_data(
      legacy,
      data.frame(g = factor("s2"), a = 2, b = 3)
    ),
    paste0(
      "Structured random-effect index metadata for '", block$structured_index$name,
      "' are missing the fitted level keys. Refit the model with this version ",
      "of BayesTools."
    ),
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
})

test_that("numeric index prediction falls back to unambiguous fitted labels", {
  compile <- function(formula, data){
    JAGS_formula(formula, "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("point", list(location = 1)))
    )$formula_design$random_effects[[1]]
  }
  indicator <- function(block, data){
    .bt_random_effect_prediction_data(block, data)$model_matrix
  }
  unit_row <- function(column, n_columns){
    matrix(replace(numeric(n_columns), column, 1), nrow = 1L)
  }

  # seq() yields 0.30000000000000004, displayed as "0.3".
  times <- seq(0, 0.5, by = 0.1)
  data <- data.frame(g = factor(rep(c("s1", "s2"), each = 6)), t = rep(times, 2))
  for(structure in c("cs", "ar1")){
    block <- compile(
      stats::as.formula(paste0("~ 1 + ", structure, "(t | g)")),
      data
    )
    expect_identical(block$column_names[4L], "t0.3", info = structure)
    # Exact key.
    expect_equal(
      indicator(block, data.frame(g = factor("s1"), t = times[4L])),
      unit_row(4L, 6L),
      ignore_attr = TRUE,
      info = structure
    )
    # A different double with the unique fitted label "0.3".
    expect_false(identical(0.3, times[4L]))
    expect_equal(
      indicator(block, data.frame(g = factor("s1"), t = c(0.3, 0.5))),
      rbind(unit_row(4L, 6L), unit_row(6L, 6L)),
      ignore_attr = TRUE,
      info = structure
    )
    expect_equal(
      indicator(block, data.frame(g = factor("s1"), t = "0.3")),
      unit_row(4L, 6L),
      ignore_attr = TRUE,
      info = structure
    )
    # A value whose label matches no fitted level is a new index level.
    expect_error(
      indicator(block, data.frame(g = factor("s1"), t = 0.35)),
      "do not match the levels used for model specification",
      info = structure
    )
  }

  # When the label identifies two fitted levels, it identifies neither.
  ambiguous <- compile(
    ~ 1 + cs(t | g),
    data.frame(g = factor(rep(c("s1", "s2"), each = 2)),
               t = c(0.3, 0.1 + 0.2, 0.3, 0.1 + 0.2))
  )
  near <- 0.3 + 2e-16
  expect_identical(as.character(near), "0.3")
  expect_false(near %in% c(0.3, 0.1 + 0.2))
  expect_error(
    indicator(ambiguous, data.frame(g = factor("s1"), t = near)),
    "do not match the levels used for model specification"
  )
  expect_equal(
    indicator(ambiguous, data.frame(g = factor("s1"), t = 0.1 + 0.2)),
    unit_row(2L, 2L),
    ignore_attr = TRUE
  )

  # Multi-variable cells resolve each numeric component the same way.
  cells <- compile(
    ~ 1 + cs(a + b | g),
    data.frame(g = factor(rep(c("s1", "s2"), each = 2)),
               a = rep(c(0.1 * 3, 1), 2), b = rep(c(2, 2), 2))
  )
  expect_identical(cells$column_names, c("a_b0.3.2", "a_b1.2"))
  expect_equal(
    indicator(cells, data.frame(g = factor("s1"), a = 0.3, b = 2)),
    unit_row(1L, 2L),
    ignore_attr = TRUE
  )
})

test_that("expression-only subtraction preserves structural zero intercept", {
  expect_equal(.remove_expressions(~ expression(x) - 1), ~ 1 - 1,
    ignore_formula_env = TRUE)
  result <- JAGS_formula(~ expression(x[i]) - 1, "mu", data.frame(x = 1:3),
                        list(intercept = prior("point", list(location = 0))))
  expect_true(is.prior.point(result$prior_list$mu_intercept))
  expect_equal(result$prior_list$mu_intercept$parameters$location, 0)
  expect_identical(result$formula_design$transformed_terms, list("x[i]"))
})

test_that("cumulative-level factor metadata retains all and only its levels", {
  p <- prior_ordered(prior("normal", list(0, 1)), contrast = "cumulative_levels")
  attr(p, "levels") <- 4L
  attr(p, "factor_terms") <- "f"
  attr(p, "factor_contrasts") <- c(f = "contr.ordered_cumulative_levels")
  expect_identical(.factor_level_list(p)$f, as.character(1:4))
  expect_equal(dim(.factor_term_design_from_metadata(p)$design), c(4L, 4L))
  samples <- matrix(1:8, 2, 4)
  class(samples) <- c("mixed_posteriors.factor", class(samples))
  attr(samples, "ordered") <- TRUE
  attr(samples, "levels") <- 4L
  attr(samples, "factor_terms") <- "f"
  attr(samples, "factor_contrasts") <- c(f = "contr.ordered_cumulative_levels")
  transformed <- transform_factor_samples(list(f = samples))$f
  expect_equal(dim(transformed), c(2L, 4L))
  expect_identical(attr(transformed, "level_names"), as.character(1:4))
})

test_that("one-component random-effect gates do not create artificial weights", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"))
  )
  result <- JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "component",
        terms = c(study = "study"),
        sd = prior("gamma", list(2, 2)),
        inclusion = list(study = prior("spike", list(location = 0.5)))
      )
    )
  )

  allocation <- result$formula_design$random_effects[[1]]$
    sd_binding$allocations[[1]]
  expect_true(allocation$gate_only)
  expect_null(allocation$weight_name)
  expect_false(any(grepl("__weight", names(result$prior_list), fixed = TRUE)))
  expect_match(
    result$formula_syntax,
    paste0(
      "mu__xREx__study_intercept = ",
      "mu__xRE_ALLOCx_component__allocation_sd * ",
      "mu__xRE_ALLOCx_component__include_study_indicator"
    ),
    fixed = TRUE
  )

  samples <- matrix(
    c(1, 0, 2, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        allocation$source_node,
        allocation$inclusion$study$indicator_name
      )
    )
  )
  sd_draws <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = samples,
    prior_list = result$prior_list
  )
  expect_equal(sd_draws[, 1], c(0, 2))

  fit <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(samples)),
      sample = nrow(samples)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  fit <- attach_test_parameter_map(fit)
  catalog <- parameter_catalog(fit)
  inclusion <- parameter_catalog_resolve(
    catalog,
    "(mu) component: inclusion(study)"
  )
  sd <- parameter_catalog_resolve(catalog, "(mu) sd(intercept)")
  expect_equal(as.numeric(as.matrix(parameter_draws(fit, inclusion))), c(0, 1))
  expect_equal(as.numeric(as.matrix(parameter_draws(fit, sd))), c(0, 2))
  expect_false(any(catalog$quantities$quantity == "var_prop"))

  inference <- JAGS_inference_table(fit)
  expect_equal(rownames(inference), "(mu) component: inclusion(study)")
  expect_equal(as.numeric(inference$prior_prob), 0.5)
  expect_equal(as.numeric(inference$post_prob), 0.5)
  expect_equal(as.numeric(inference$inclusion_BF), 1)
  expect_identical(attr(inference, "parameter_roles"), "random_inclusion")

  conditional_samples <- JAGS_estimates_table(
    fit,
    keep_parameters = "random",
    random_effects_summary = "standard",
    conditional = TRUE,
    remove_inclusion = TRUE,
    return_samples = TRUE
  )
  expect_equal(
    conditional_samples[, "(mu) sd(intercept)"],
    c(NA_real_, 2)
  )
})

test_that("random formulas reject transformations and persist tuple ordering", {
  for(random_formula in list(
    ~ 1 + diag(1 + I(x^2) | g1),
    ~ 1 + diag(1 + scale(x) | g1),
    ~ 1 + diag(1 + custom_transform(x) | g1)
  )){
    expect_error(
      BayesTools:::.bt_parse_random_effects(random_formula),
      "Create the transformed predictor as an explicit data column",
      fixed = TRUE
    )
  }

  data <- data.frame(
    x = rep(c(-1, 1), 4),
    g1 = factor(rep(c("b", "a", "b", "a"), 2), levels = c("a", "b")),
    g2 = factor(rep(c("y", "y", "x", "x"), 2), levels = c("x", "y"))
  )
  sd_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  result <- JAGS_formula(
    ~ 1 + diag(1 + x | g1:g2),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g1_g2 = random_block(sd = sd_prior)
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]
  expected_group <- interaction(
    data$g1,
    data$g2,
    drop = TRUE,
    sep = ":",
    lex.order = TRUE
  )

  expect_identical(random_term$group_levels, levels(expected_group))
  expect_identical(random_term$group_map, as.integer(expected_group))
  expect_identical(random_term$group_components, c("g1", "g2"))
  expect_identical(
    random_term$group_tuples,
    matrix(
      c("a", "x", "a", "y", "b", "x", "b", "y"),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("g1", "g2"))
    )
  )
  expect_identical(
    unname(random_term$group_tuple_index[random_term$group_tuple_keys]),
    seq_along(random_term$group_tuple_keys)
  )

  shuffled <- data[c(4, 1, 3, 2), , drop = FALSE]
  prediction <- BayesTools:::.bt_random_effect_prediction_data(
    random_term,
    shuffled,
    group_data = shuffled
  )
  expected_shuffled <- interaction(
    shuffled$g1,
    shuffled$g2,
    drop = TRUE,
    sep = ":",
    lex.order = TRUE
  )
  expect_identical(
    prediction$group_map,
    match(as.character(expected_shuffled), random_term$group_levels)
  )

  nested_result <- JAGS_formula(
    ~ 1 + diag(1 | g1 / g2),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = sd_prior)
  )
  nested_term <- nested_result$formula_design$random_effects[[1L]]
  expect_identical(nested_term$group_components, c("g1", "g2"))
  expect_identical(nested_term$group_levels, levels(expected_group))
  expect_identical(nested_term$group_map, as.integer(expected_group))
  expect_identical(
    nested_term$group_tuples,
    matrix(
      c("a", "x", "a", "y", "b", "x", "b", "y"),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("g1", "g2"))
    )
  )

  collision_data <- data.frame(
    g1 = factor(c("a", "a:b"), levels = c("a", "a:b")),
    g2 = factor(c("b:c", "c"), levels = c("b:c", "c"))
  )
  collision_result <- JAGS_formula(
    ~ 1 + diag(1 | g1:g2),
    "mu",
    collision_data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g1_g2 = random_block(sd = sd_prior)
    )
  )
  collision_term <- collision_result$formula_design$random_effects[[1L]]
  expect_identical(collision_term$n_groups, 2L)
  expect_true(anyDuplicated(collision_term$group_labels) > 0L)
  expect_false(anyDuplicated(collision_term$group_tuple_keys) > 0L)
  expect_false(anyDuplicated(collision_term$group_levels) > 0L)
  expect_identical(
    BayesTools:::.bt_random_group_unique_labels(
      "a:b:c",
      "3:a:b|1:c",
      existing = "a:b:c"
    ),
    "a:b:c [3:a:b|1:c]"
  )
})

test_that("nested grouping expands consistently across covariance specials and wrappers", {

  block_labels <- function(formula){
    terms <- BayesTools:::.bt_parse_random_effects(formula)$terms
    list(
      labels     = vapply(terms, function(t) t$group_label, character(1)),
      structures = vapply(terms, function(t) t$structure, character(1)),
      blocks     = vapply(terms, function(t) t$block_name, character(1)),
      nesting    = lapply(terms, function(t) t$group_nesting_components)
    )
  }

  # The plain-bar reference expansion (already supported) defines the canonical
  # order/interaction convention every wrapped structure must match.
  reference <- block_labels(~ 1 + (1 | g1 / g2))
  expect_equal(reference$labels, c("g2:g1", "g1"))
  expect_equal(reference$blocks, c("g2_g1", "g1"))
  expect_identical(reference$nesting, list(c("g1", "g2"), "g1"))
  parenthesized_reference <- block_labels(~ 1 + random(1 | (g1 / g2)))
  expect_equal(parenthesized_reference$labels, c("g2:g1", "g1"))
  expect_equal(parenthesized_reference$blocks, c("g2_g1", "g1"))
  expect_identical(parenthesized_reference$nesting, reference$nesting)
  parenthesized_special <- block_labels(~ 1 + diag(1 | (g1 / g2)))
  expect_equal(parenthesized_special$labels, c("g2:g1", "g1"))
  expect_equal(parenthesized_special$blocks, c("g2_g1", "g1"))
  expect_identical(parenthesized_special$nesting, reference$nesting)

  for(structure in c("diag", "id", "us", "un", "cs", "hcs", "ar1", "ar", "har", "car")){
    parsed <- block_labels(stats::as.formula(
      paste0("~ 1 + x + ", structure, "(1 + x | g1 / g2)")
    ))
    expect_equal(parsed$labels, c("g2:g1", "g1"),
                 info = paste("labels for", structure))
    expect_equal(parsed$blocks, c("g2_g1", "g1"),
                 info = paste("blocks for", structure))
    expect_equal(length(unique(parsed$structures)), 1L,
                 info = paste("single structure for", structure))
  }

  for(wrapper in c("random", "re")){
    parsed <- block_labels(stats::as.formula(
      paste0("~ 1 + ", wrapper, "(1 | g1 / g2)")
    ))
    expect_equal(parsed$labels, c("g2:g1", "g1"))
    expect_equal(parsed$blocks, c("g2_g1", "g1"))
    expect_equal(parsed$structures, c("us", "us"))
  }

  # Three-level nesting expands to one block per cumulative level.
  three <- block_labels(~ 1 + diag(1 | g1 / g2 / g3))
  expect_equal(three$labels, c("g3:g2:g1", "g2:g1", "g1"))
  expect_equal(three$blocks, c("g3_g2_g1", "g2_g1", "g1"))
  expect_identical(
    three$nesting,
    list(c("g1", "g2", "g3"), c("g1", "g2"), "g1")
  )

  # Formula order is preserved when wrapped and plain terms are mixed.
  mixed <- block_labels(~ 1 + x + diag(1 | g1 / g2) + random(1 | h))
  expect_equal(mixed$blocks, c("g2_g1", "g1", "h"))
  expect_equal(mixed$structures, c("diag", "diag", "us"))

  # Named covariance / hom arguments propagate to every expanded sub-term.
  propagated <- block_labels(~ 1 + x + random(1 + x | g1 / g2, covariance = "cs"))
  expect_equal(propagated$structures, c("cs", "cs"))
  hom_propagated <- BayesTools:::.bt_parse_random_effects(~ 1 + cs(1 + x | g1 / g2, hom = TRUE))$terms
  expect_true(all(vapply(hom_propagated, function(t) isTRUE(t$hom), logical(1))))
  wrapper_hcs <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + random(1 + x | g1, covariance = "cs", hom = FALSE)
  )$terms[[1]]
  expect_equal(wrapper_hcs$structure, "hcs")
  expect_null(wrapper_hcs$hom)
  wrapper_har <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + re(1 + x | g1, covariance = "ar1", hom = FALSE)
  )$terms[[1]]
  expect_equal(wrapper_har$structure, "har")
  expect_null(wrapper_har$hom)
  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + random(cs(1 + x | g1))),
    "Do not wrap covariance-special random-effect calls",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + re(ar1(1 + x | g1), name = "time")),
    "Do not wrap covariance-special random-effect calls",
    fixed = TRUE
  )
  expect_error(
    random_effects_formula(~ foo() + random(1 | g1)),
    "only random-effect terms",
    fixed = TRUE
  )

  # '||' under a wrapper expands and still resolves to diagonal blocks, matching
  # the plain-bar '||' nesting expansion.
  double_bar <- block_labels(~ 1 + x + random(1 + x || g1 / g2))
  expect_equal(double_bar$structures, c("diag", "diag"))
  expect_equal(double_bar$blocks, c("g2_g1", "g1"))

  # An explicit block name becomes a per-level prefix on each expanded block,
  # preserving the user's naming intent while keeping the blocks distinct.
  named_two <- block_labels(~ 1 + random(1 | g1 / g2, name = "foo"))
  expect_equal(named_two$blocks, c("foo_g2_g1", "foo_g1"))
  named_special <- block_labels(~ 1 + diag(1 | g1 / g2, name = "foo"))
  expect_equal(named_special$blocks, c("foo_g2_g1", "foo_g1"))
  named_three <- block_labels(~ 1 + random(1 | g1 / g2 / g3, name = "spatial"))
  expect_equal(named_three$blocks, c("spatial_g3_g2_g1", "spatial_g2_g1", "spatial_g1"))

  # End-to-end: numeric grouping IDs must build two nested blocks, not a single
  # block keyed on the literal 'g1 / g2' quotient.
  df_numeric <- data.frame(
    g1 = rep(1:4, each = 10),               # 4 outer groups
    g2 = rep(rep(1:2, each = 5), 4),        # 2 inner-of-g1 -> 8 (g2:g1) groups
    g3 = rep(1:5, 8)                        # 5 inner-of-(g2:g1) -> 40 (g3:g2:g1) groups
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  map_level_counts <- function(formula){
    built <- JAGS_formula(
      formula,
      parameter    = "mu",
      data         = df_numeric,
      prior_list   = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = sd_prior)
    )
    map_columns <- grep("xRE_MAPx", names(built$data), value = TRUE)
    sort(vapply(map_columns, function(m) length(unique(built$data[[m]])), integer(1)))
  }

  # Two-level nesting: inner (g2:g1) has 8 groups, outer (g1) has 4.
  for(formula in list(
    ~ 1 + (1 | g1 / g2),
    ~ 1 + diag(1 | g1 / g2),
    ~ 1 + random(1 | g1 / g2),
    ~ 1 + diag(1 | (g1 / g2)),
    ~ 1 + random(1 | (g1 / g2))
  )){
    expect_equal(unname(map_level_counts(formula)), c(4L, 8L))
  }

  # Three-level nesting: 40 / 8 / 4 across the cumulative levels.
  for(formula in list(~ 1 + (1 | g1 / g2 / g3), ~ 1 + diag(1 | g1 / g2 / g3), ~ 1 + random(1 | g1 / g2 / g3))){
    expect_equal(unname(map_level_counts(formula)), c(4L, 8L, 40L))
  }
})

test_that("a plain double-bar term with nested grouping expands like random()", {

  block_labels <- function(formula){
    terms <- BayesTools:::.bt_parse_random_effects(formula)$terms
    list(
      labels     = vapply(terms, function(t) t$group_label, character(1)),
      structures = vapply(terms, function(t) t$structure, character(1)),
      blocks     = vapply(terms, function(t) t$block_name, character(1)),
      nesting    = lapply(terms, function(t) t$group_nesting_components)
    )
  }

  for(term in c("1 + x", "1", "0 + x")){
    plain <- block_labels(stats::as.formula(paste0("~ 1 + x + (", term, " || g1 / g2)")))
    wrapped <- block_labels(stats::as.formula(paste0("~ 1 + x + random(", term, " || g1 / g2)")))
    expect_identical(plain$blocks, c("g2_g1", "g1"), info = term)
    expect_identical(plain$labels, c("g2:g1", "g1"), info = term)
    expect_identical(plain$structures, c("diag", "diag"), info = term)
    expect_identical(plain$nesting, list(c("g1", "g2"), "g1"), info = term)
    expect_identical(plain, wrapped, info = term)
  }

  expected_formula <- "diag(1 + x | g2:g1, name = \"g2_g1\") + diag(1 + x | g1, name = \"g1\")"
  expect_identical(
    deparse(random_effects_formula(~ (1 + x || g1 / g2))$formula[[2L]]),
    expected_formula
  )
  expect_identical(
    deparse(random_effects_formula(~ random(1 + x || g1 / g2))$formula[[2L]]),
    expected_formula
  )

  set.seed(1)
  df <- data.frame(
    x = rnorm(16),
    g1 = factor(rep(1:2, each = 8)),
    g2 = factor(rep(1:4, each = 4))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  compile <- function(formula){
    JAGS_formula(formula, "mu", df, list(intercept = prior("normal", list(0, 1))),
                 prior_random = prior_random(sd = sd_prior))
  }
  plain <- compile(~ 1 + (1 + x || g1 / g2))
  wrapped <- compile(~ 1 + random(1 + x || g1 / g2))
  blocks <- plain$formula_design$random_effects
  expect_identical(
    vapply(blocks, function(block) block$block_name, character(1)),
    c("g2_g1", "g1")
  )
  expect_identical(
    vapply(blocks, function(block) block$structure, character(1)),
    c("diag", "diag")
  )
  expect_identical(plain$formula_syntax, wrapped$formula_syntax)
  expect_identical(names(plain$prior_list), names(wrapped$prior_list))
  expect_identical(
    vapply(compile(~ 1 + (1 || g1 / g2))$formula_design$random_effects,
           function(block) block$block_name, character(1)),
    c("g2_g1", "g1")
  )
})

test_that("covariance specials preserve parenthesized bar terms", {

  diag_parsed <- BayesTools:::.bt_parse_random_effects(~ diag((1 | g)))
  expect_length(diag_parsed$terms, 1L)
  expect_equal(diag_parsed$terms[[1L]]$structure, "diag")
  expect_equal(diag_parsed$terms[[1L]]$term_formula, ~ 1, ignore_formula_env = TRUE)
  expect_equal(diag_parsed$fixed_formula, ~ 1, ignore_formula_env = TRUE)

  cs_parsed <- BayesTools:::.bt_parse_random_effects(~ cs((x | g)))
  expect_length(cs_parsed$terms, 1L)
  expect_equal(cs_parsed$terms[[1L]]$structure, "cs")
  expect_equal(cs_parsed$terms[[1L]]$term_formula, ~ x, ignore_formula_env = TRUE)
  expect_equal(cs_parsed$fixed_formula, ~ 1, ignore_formula_env = TRUE)

  named_diag <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + diag((1 | g), name = "diag_block")
  )
  expect_equal(named_diag$terms[[1L]]$block_name, "diag_block")
  expect_equal(named_diag$terms[[1L]]$structure, "diag")
  expect_equal(named_diag$fixed_formula, ~ 1, ignore_formula_env = TRUE)

  named_cs <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + cs((x | g), name = "cs_block", hom = TRUE)
  )
  expect_equal(named_cs$terms[[1L]]$block_name, "cs_block")
  expect_equal(named_cs$terms[[1L]]$structure, "cs")
  expect_true(named_cs$terms[[1L]]$hom)
  expect_equal(named_cs$fixed_formula, ~ 1, ignore_formula_env = TRUE)
})

test_that("random-effect formula lists retain component hierarchy", {

  df <- data.frame(
    district = factor(c("d1", "d1", "d2", "d2")),
    school   = factor(c("s1", "s2", "s1", "s2")),
    id       = factor(c("i1", "i1", "i2", "i2")),
    study    = factor(c("a", "a", "b", "b"))
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )

  random_effects <- random_effects_formula(
    list(
      nested = ~ 1 | district / school,
      study  = ~ 1 | study
    )
  )
  expect_equal(
    vapply(random_effects$terms, `[[`, character(1), "block_name"),
    c("nested_school_district", "nested_district", "study")
  )
  expect_equal(
    random_effects$components,
    list(nested = c("nested_school_district", "nested_district"), study = "study")
  )

  expect_equal(
    attr(random_effects$formula, "random_components"),
    random_effects$components
  )
  expect_equal(
    vapply(attr(random_effects$formula, "random_terms"), `[[`, character(1), "block_name"),
    c("nested_school_district", "nested_district", "study")
  )

  explicit_named_component <- random_effects_formula(
    list(study = ~ random(1 | id, name = "id"))
  )
  expect_equal(
    vapply(explicit_named_component$terms, `[[`, character(1), "block_name"),
    "id"
  )
  expect_equal(explicit_named_component$components, list(study = "id"))
  expect_equal(explicit_named_component$terms[[1]]$component_label, "study")
  expect_true(explicit_named_component$terms[[1]]$component_visible)

  explicit_named_result <- JAGS_formula(
    formula = explicit_named_component$formula,
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(explicit_named_result$formula_design$random_effects[[1]]$block_name, "id")
  expect_error(
    JAGS_formula(
      formula = explicit_named_component$formula,
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(study = random_block(sd = sd_prior))
    ),
    "block override names were not found",
    fixed = TRUE
  )

  random_prior <- prior_random(
    random_variance_allocation(
      name  = "random_total",
      terms = c(nested = "nested", study = "study"),
      sd    = sd_prior
    ),
    random_variance_allocation(
      name   = "nested_split",
      terms  = c(school_district = "nested_school_district", district = "nested_district"),
      parent = allocation_ref("random_total", "nested")
    )
  )
  allocations <- BayesTools:::.bt_random_allocation_list(random_prior$allocation)
  expect_length(allocations, 2L)
  expect_equal(allocations[[1]]$terms, c(nested = "nested", study = "study"))
  expect_equal(allocations[[2]]$parent$allocation, "random_total")
  expect_equal(allocations[[2]]$parent$component, "nested")
  expect_equal(
    allocations[[2]]$terms,
    c(school_district = "nested_school_district", district = "nested_district")
  )

  result <- JAGS_formula(
    formula      = random_effects$formula,
    parameter    = "mu",
    data         = df,
    prior_list   = list(intercept = prior("normal", list(0, 1))),
    prior_random = random_prior
  )
  expect_match(result$formula_syntax, "random_total__component_nested_sd", fixed = TRUE)
  expect_equal(length(result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors), 2L)
  expect_equal(length(result$formula_design$random_effects[[3]]$sd_binding$allocations[[1L]]$factors), 1L)

  expect_error(
    random_effects_formula(
      list(
        `study-component` = ~ 1 | study,
        study_component   = ~ 1 | district
      )
    ),
    "unique after sanitization"
  )
  expect_error(
    random_effects_formula(
      ~ us(1 | study, name = "shared") + us(1 | district, name = "shared")
    ),
    "block names must be unique"
  )
  expect_error(
    random_effects_formula(~ x + (1 | study)),
    "only random-effect terms"
  )
})

test_that("random group covariance constructor validates and scales kernels", {

  K <- matrix(
    c(4, 1, 1,
      1, 9, 2,
      1, 2, 16),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("b", "a", "extra"), c("b", "a", "extra"))
  )
  none <- random_group_covariance(K, scale = "none")
  expect_s3_class(none, "random_group_covariance")
  expect_equal(none$covariance, K, ignore_attr = TRUE)
  expect_equal(none$scale, "none")
  expect_identical(none$triangle_source, "complete")

  lower <- matrix(
    c(2, .4, 0, 3),
    nrow = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  lower_covariance <- random_group_covariance(lower, scale = "none")
  expect_identical(lower_covariance$triangle_source, "lower")
  expect_equal(
    lower_covariance$covariance,
    matrix(
      c(2, .4, .4, 3),
      nrow = 2,
      dimnames = dimnames(lower)
    ),
    ignore_attr = TRUE
  )

  upper <- t(lower)
  upper_covariance <- random_group_covariance(upper, scale = "none")
  expect_identical(upper_covariance$triangle_source, "upper")
  expect_equal(
    upper_covariance$covariance,
    lower_covariance$covariance,
    ignore_attr = TRUE
  )

  asymmetric <- matrix(
    c(2, .4, .4 + .Machine$double.eps, 3),
    nrow = 2,
    dimnames = dimnames(lower)
  )
  expect_error(
    random_group_covariance(asymmetric, scale = "none"),
    "must be exactly symmetric",
    fixed = TRUE
  )

  prepared_none <- BayesTools:::.bt_prepare_group_covariance_kernel(
    none,
    group_levels = c("a", "b"),
    block_name = "id"
  )
  expect_equal(prepared_none$levels, c("a", "b"))
  expect_equal(prepared_none$dropped_levels, "extra")
  expect_equal(prepared_none$kernel, K[c("a", "b"), c("a", "b")])
  expect_true(all(is.finite(prepared_none$precision)))
  expect_true(is.finite(prepared_none$log_det))

  K2 <- matrix(
    c(4, 1, 1, 9),
    nrow = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  prepared_cor <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K2, scale = "cor"),
    group_levels = c("a", "b"),
    block_name = "id"
  )
  expect_equal(prepared_cor$kernel, stats::cov2cor(K2), tolerance = 1e-12)

  prepared_cor0 <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K2, scale = "cor0"),
    group_levels = c("a", "b"),
    block_name = "id"
  )
  R <- stats::cov2cor(K2)
  expect_equal(prepared_cor0$kernel, (R - min(R)) / (1 - min(R)), tolerance = 1e-12)

  prepared_cov0 <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K2, scale = "cov0"),
    group_levels = c("a", "b"),
    block_name = "id"
  )
  expect_equal(prepared_cov0$kernel, K2 - min(K2), tolerance = 1e-12)

  expect_error(
    BayesTools:::.bt_prepare_group_covariance_kernel(
      random_group_covariance(K2, scale = "none"),
      group_levels = c("a", "missing"),
      block_name = "id"
    ),
    "missing fitted level",
    fixed = TRUE
  )
  singular <- matrix(
    c(1, 1, 1, 1),
    nrow = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  expect_error(
    BayesTools:::.bt_prepare_group_covariance_kernel(
      random_group_covariance(singular, scale = "none"),
      group_levels = c("a", "b"),
      block_name = "id"
    ),
    "positive definite",
    fixed = TRUE
  )

  nonfinite_precision <- diag(c(1e-310, 1))
  dimnames(nonfinite_precision) <- list(c("a", "b"), c("a", "b"))
  expect_error(
    BayesTools:::.bt_prepare_group_covariance_kernel(
      random_group_covariance(nonfinite_precision, scale = "none"),
      group_levels = c("a", "b"),
      block_name = "id"
    ),
    "finite precision and log-determinant metadata",
    fixed = TRUE
  )

  expect_error(
    random_group_covariance(matrix(1, 1, 1), scale = "none"),
    "row and column names",
    fixed = TRUE
  )
  bad <- K2
  dimnames(bad) <- list(c("a", "a"), c("a", "b"))
  expect_error(
    random_group_covariance(bad, scale = "none"),
    "row names must be unique",
    fixed = TRUE
  )
})

test_that("known heteroscedastic covariance scaling preserves exact symmetry", {

  # X_i = s_i * (sqrt(.4) * U + sqrt(.6) * E_i), independent standard
  # normals U/E_i, with s = (1.7, 2.1, 3.7). These literal covariances
  # therefore have independently known common correlation .4.
  K <- matrix(c(2.89, 1.428, 2.516,
                1.428, 4.41, 3.108,
                2.516, 3.108, 13.69), 3L, byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  expected_cor <- matrix(.4, 3L, 3L, dimnames = dimnames(K))
  diag(expected_cor) <- 1
  prepared <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K, scale = "cor"),
    group_levels = rownames(K), block_name = "study")
  expect_identical(prepared$kernel, t(prepared$kernel))
  expect_equal(prepared$kernel, expected_cor, tolerance = 1e-14)

  # This positive-definite Gram-plus-diagonal fixture made cov2cor() differ
  # across triangles by 1.11e-16 and rejected the public formula constructor.
  K <- tcrossprod(matrix(c(.1, .2, .3, .7, 1.1, .2, .4, 1.3),
    4L, 2L, byrow = TRUE)) + diag(.1, 4L)
  labels <- letters[1:4]
  dimnames(K) <- list(labels, labels)
  original <- K
  formula <- random_effects_formula(~ 1 | study,
    group_covariance = random_group_covariance(K, scale = "cor"))
  result <- JAGS_formula(formula = formula, parameter = "mu",
    data = data.frame(study = factor(rep(labels, each = 2L), levels = labels)),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("normal", list(0, 1),
      truncation = list(lower = 0, upper = Inf))))
  actual <- result$formula_design$random_effects[[1L]]$group_covariance$kernel
  expect_identical(actual, t(actual))
  expect_identical(dimnames(actual), dimnames(original))
  expect_equal(actual, stats::cov2cor(original), tolerance = 1e-12)
  expect_identical(K, original)
})

test_that("known group kernels require resolvable positive definiteness", {

  # Covariance of (U + W, U, W), with independent unit-variance U and W.
  # The first column is exactly the sum of the other two, even when raw
  # floating-point Cholesky returns three positive pivots.
  K <- matrix(c(2, 1, 1, 1, 1, 0, 1, 0, 1), 3L,
    dimnames = list(letters[1:3], letters[1:3]))
  message <- paste0("Known group covariance for random-effect block 'study' is unavailable ",
    "because positive definiteness after scaling cannot be resolved at working precision.")
  for(scale in c("none", "cor")){
    specification <- random_group_covariance(K, scale = scale)
    expect_error(BayesTools:::.bt_prepare_group_covariance_kernel(specification,
      group_levels = letters[1:3], block_name = "study"), message, fixed = TRUE)
  }
  # Subsetting remains prior to the support check: the fitted principal
  # submatrix can be positive definite even when unused levels are dependent.
  subset <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K, scale = "none"), group_levels = c("a", "b"),
    block_name = "study")
  expect_identical(subset$kernel, K[c("a", "b"), c("a", "b")])

  # Rescaling units must not turn small independent variances into null axes.
  independent <- diag(c(1e-24, 1, 1e12))
  dimnames(independent) <- dimnames(K)
  prepared <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(independent, scale = "none"),
    group_levels = letters[1:3], block_name = "study")
  expect_identical(prepared$kernel, independent)
  near_singular <- matrix(c(1, 1 - .Machine$double.eps,
    1 - .Machine$double.eps, 1), 2L,
    dimnames = list(letters[1:2], letters[1:2]))
  expect_error(BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(near_singular, scale = "none"),
    group_levels = letters[1:2], block_name = "study"), message, fixed = TRUE)
})

test_that("random_effects_formula attaches known group covariance by block names", {

  study_kernel <- diag(2)
  dimnames(study_kernel) <- list(c("a", "b"), c("a", "b"))
  district_kernel <- diag(2)
  dimnames(district_kernel) <- list(c("d1", "d2"), c("d1", "d2"))

  by_block <- random_effects_formula(
    ~ random(1 | study, name = "study_block") +
      random(1 | district, name = "district_block"),
    group_covariance = list(
      study_block = random_group_covariance(study_kernel, scale = "none"),
      district_block = district_kernel
    )
  )
  expect_equal(
    vapply(by_block$terms, function(term) term$group_covariance$scale, character(1)),
    c("none", "cor")
  )
  expect_equal(names(by_block$group_covariance), c("study_block", "district_block"))

  by_group <- random_effects_formula(
    ~ random(1 | study, name = "study_block"),
    group_covariance = list(study = random_group_covariance(study_kernel))
  )
  expect_s3_class(by_group$terms[[1]]$group_covariance, "random_group_covariance")

  expect_error(
    random_effects_formula(
      ~ random(1 | study, name = "a") + random(1 | district, name = "b"),
      group_covariance = random_group_covariance(study_kernel)
    ),
    "single random-effect block",
    fixed = TRUE
  )
  expect_error(
    random_effects_formula(
      ~ random(1 | study, name = "a"),
      group_covariance = list(missing = random_group_covariance(study_kernel))
    ),
    "does not match",
    fixed = TRUE
  )
})

test_that("random-effect design exposes grouping maps and correlated syntax", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    idx = factor(rep(c("t1", "t2"), 3), levels = c("t1", "t2")),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )
  random_term <- BayesTools:::.bt_parse_random_effects(~ 1 + x + (1 + x || id))$terms[[1]]

  expect_s3_class(random_term, "BayesTools_random_effect_term")
  expect_equal(random_term$group_label, "id")
  expect_equal(random_term$structure, "diag")
  expect_false("covariance" %in% names(random_term))
  expect_equal(attr(random_term, "structure"), "diag")
  expect_null(attr(random_term, "covariance"))
  expect_true(isTRUE(random_term$independent))
  expect_equal(random_term$term_formula, ~ 1 + x, ignore_formula_env = TRUE)
  expect_equal(attr(random_term, "grouping_factor"), "id")
  expect_true(isTRUE(attr(random_term, "independent")))

  reordered_wrapper <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x + random(name = "study", covariance = "diag", 1 + x | id)
  )$terms[[1]]
  expect_equal(reordered_wrapper$block_name, "study")
  expect_equal(reordered_wrapper$structure, "diag")
  expect_equal(reordered_wrapper$term_formula, ~ 1 + x, ignore_formula_env = TRUE)

  reordered_special <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x + diag(name = "study", 1 + x | id)
  )$terms[[1]]
  expect_equal(reordered_special$block_name, "study")
  expect_equal(reordered_special$structure, "diag")

  explicit_diag_double_bar <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x + random(1 + x || id, name = "study", covariance = "diag")
  )$terms[[1]]
  expect_equal(explicit_diag_double_bar$block_name, "study")
  expect_equal(explicit_diag_double_bar$structure, "diag")

  expect_error(
    BayesTools:::.bt_parse_random_effects(
      ~ 1 + x + random(1 + x || id, name = "study", covariance = "us")
    ),
    "cannot be combined with '||' syntax",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parse_random_effects(
      ~ 1 + x +
        random(1 | drug, name = "drug", covariance = "diag") +
        us(1 + x || id)
    ),
    "cannot be combined with '||' syntax",
    fixed = TRUE
  )

  bare_and_wrapper <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + (1 | id) + random(1 | id, name = "study", covariance = "diag")
  )
  expect_equal(
    vapply(bare_and_wrapper$terms, function(term) term$block_name, character(1)),
    c("id", "study")
  )
  expect_equal(
    vapply(bare_and_wrapper$terms, function(term) term$structure, character(1)),
    c("us", "diag")
  )

  bare_and_special <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + (1 | id) + diag(name = "study", 1 | id)
  )
  expect_equal(
    vapply(bare_and_special$terms, function(term) term$block_name, character(1)),
    c("id", "study")
  )
  expect_equal(
    vapply(bare_and_special$terms, function(term) term$structure, character(1)),
    c("us", "diag")
  )

  ordered_wrapper_terms <- BayesTools:::.bt_parse_random_effects(
    ~ 1 +
      random(1 | id, name = "study", covariance = "diag") +
      diag(1 | drug, name = "drug")
  )
  expect_equal(
    vapply(ordered_wrapper_terms$terms, function(term) term$block_name, character(1)),
    c("study", "drug")
  )

  reversed_wrapper_terms <- BayesTools:::.bt_parse_random_effects(
    ~ 1 +
      diag(1 | drug, name = "drug") +
      random(1 | id, name = "study", covariance = "diag")
  )
  expect_equal(
    vapply(reversed_wrapper_terms$terms, function(term) term$block_name, character(1)),
    c("drug", "study")
  )

  formula_env <- new.env(parent = baseenv())
  formula_env$make_group_for_random_effect_test <- function(x) paste0("g_", x)
  grouped_formula <- ~ 1 + random(1 | make_group_for_random_effect_test(id), name = "derived", covariance = "diag")
  environment(grouped_formula) <- formula_env
  expect_error(
    BayesTools:::.bt_parse_random_effects(grouped_formula),
    "grouping expressions must be variables",
    fixed = TRUE
  )

  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + random(1 | id, TRUE, name = "id")),
    "exactly one unnamed term",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + diag(1 | id, TRUE)),
    "exactly one unnamed term",
    fixed = TRUE
  )

  correlated <- BayesTools:::.JAGS_random_effect_formula(
    formula = BayesTools:::.bt_parse_random_effects(~ 1 + x + (1 + x | id))$terms[[1]],
    parameter = "mu",
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1),
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )
  )

  expect_equal(correlated$random_effect$structure, "us")
  expect_match(paste(correlated$random_syntax, collapse = "\n"), "dbt_lkj_cpc", fixed = TRUE)
  expect_equal(correlated$jags_modules, "BayesTools")
  expect_equal(
    correlated$add_parameters,
    c(
      "mu__xREx__id_xRE_CORx_L",
      "mu__xREx__id_xRE_CORx_R",
      "mu__xREx__id_xRE_CORx_lkj_u[1]"
    )
  )
})

test_that("independent random-effect transforms keep group-local JAGS dependencies", {

  df <- data.frame(study = rep(c("a", "b", "c"), each = 2L), x = rep(c(-1, 1), 3L))
  formulas <- list(~ 1 + (1 | study), ~ 1 + id(1 + x | study),
                   ~ 1 + diag(1 + x | study))
  for(formula in formulas){
    result <- JAGS_formula(
      formula = formula,
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("gamma", list(2, 2)))
    )
    expect_match(
      result$formula_syntax,
      "for(g in 1:3){\n     mu__xREx__study_xRE_COEFx[g,i] = mu__xREx__study_xRE_Zx[g,i] * mu__xREx__study_xRE_STDx[i]",
      fixed = TRUE
    )
    expect_false(grepl("_xRE_COEFx[1:3,", result$formula_syntax, fixed = TRUE))
  }
})

test_that("random-effect terms use structure as canonical internal metadata", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )

  parsed <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x +
      (1 + x | id) +
      random(1 + x || id, name = "independent", covariance = "diag") +
      ar1(x | id, name = "ordered")
  )

  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )

  terms <- c(parsed$terms, formula_result$formula_design$random_effects)
  expect_true(length(terms) > 0L)

  for (random_term in terms) {
    expect_s3_class(random_term, "BayesTools_random_effect_term")
    expect_true("structure" %in% names(random_term))
    expect_false("covariance" %in% names(random_term))
    expect_true(nzchar(random_term$structure))
    expect_equal(attr(random_term, "structure"), random_term$structure)
    expect_null(attr(random_term, "covariance"))
  }

  malformed_term <- terms[[1]]
  malformed_term$covariance <- malformed_term$structure
  malformed_term$structure <- NULL
  attr(malformed_term, "covariance") <- malformed_term$covariance
  attr(malformed_term, "structure") <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_structure(malformed_term),
    "missing canonical 'random_term\\$structure'"
  )
})

test_that("JAGS_formula random-effect design exposes grouping maps and public variable names", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  random_prior_diag <- prior_random(
    id = random_block(sd = prior("gamma", list(2, 2)))
  )
  random_prior_us <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior_lkj(eta = 1)
    )
  )

  result <- JAGS_formula(
    formula = ~ 1 + x + (1 + x || id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_diag
  )

  expected_map <- match(as.character(df$id), levels(df$id))

  expect_equal(result$data$mu__xREx__id_xRE_MAPx, expected_map)
  expect_equal(dim(result$data$mu__xREx__id_xRE_DATAx), c(6L, 2L))
  expect_equal(colnames(result$data$mu__xREx__id_xRE_DATAx), c("(Intercept)", "x"))
  expect_equal(unname(result$data$mu__xREx__id_xRE_DATAx[, "x"]), df$x)

  expect_equal(
    result$formula_design$jags_data_names[["__xREx__id"]],
    c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx")
  )
  expect_true(isTRUE(attr(result$formula_design$random_effects[[1]], "independent")))
  expect_equal(attr(result$formula_design$random_effects[[1]], "grouping_factor"), "id")
  expect_equal(result$formula_design$random_effects[[1]]$group_map, expected_map)
  expect_equal(result$formula_design$random_effects[[1]]$group_levels, levels(df$id))
  expect_equal(result$formula_design$random_effects[[1]]$column_names, c("(Intercept)", "x"))
  expect_equal(result$formula_design$random_effects[[1]]$jags_data_names, c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx"))
  expect_s3_class(
    result$formula_design$random_effects[[1]]$sd_leaves,
    "BayesTools_random_effect_sd_leaves"
  )
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_leaves$leaf_names_by_column,
    result$formula_design$random_effects[[1]]$sd_parameter_names
  )
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_leaves$leaf_terms_by_column,
    c("intercept", "x")
  )
  expect_equal(
    names(result$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__id_x")
  )
  expect_true(isTRUE(attr(result$prior_list$mu__xREx__id_x, "random_sd")))
  expect_equal(attr(result$prior_list$mu__xREx__id_x, "random_factor"), "id")
  expect_match(result$formula_syntax, "for\\(i in 1:3\\)")
  expect_match(result$formula_syntax, "xRE_Zx\\[i,j\\] ~ dnorm\\(0, 1\\)")

  diag_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_diag
  )

  expect_equal(
    diag_result$data$mu__xREx__id_xRE_MAPx,
    result$data$mu__xREx__id_xRE_MAPx
  )
  expect_equal(
    unname(diag_result$data$mu__xREx__id_xRE_DATAx),
    unname(result$data$mu__xREx__id_xRE_DATAx)
  )
  expect_equal(diag_result$formula_design$random_effects[[1]]$structure, "diag")

  us_result <- JAGS_formula(
    formula = ~ 1 + x + (1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_us
  )
  expect_equal(us_result$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(us_result$jags_modules, "BayesTools")
  expect_true(any(grepl("xRE_CORx_R", us_result$add_parameters, fixed = TRUE)))

  explicit_us_result <- JAGS_formula(
    formula = ~ 1 + x + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_us
  )
  expect_equal(explicit_us_result$formula_design$random_effects[[1]]$structure, "us")

  single_us_result <- JAGS_formula(
    formula = ~ 1 + (1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
  )
  expect_equal(single_us_result$jags_modules, character())
  expect_null(single_us_result$formula_design$random_effects[[1]]$correlation)
  expect_false(grepl("xRE_CORx", single_us_result$formula_syntax, fixed = TRUE))
  expect_error(
    JAGS_formula(
      formula = ~ 1 + (1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      ))
    ),
    "Single-column random-effect structure 'us' has no correlation parameter",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id, unknown = TRUE),
      parameter = "mu",
      data = df,
      prior_list = prior_list,
      prior_random = random_prior_diag
    ),
    "The 'diag' random-effect covariance structure does not support extra arguments yet.",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + random(1 | id, name = "xRE_Zx_block"),
      parameter = "mu",
      data = df,
      prior_list = prior_list,
      prior_random = prior_random(
        xRE_Zx_block = random_block(sd = prior("gamma", list(2, 2)))
      )
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + cs(1 + x | id, hom = FALSE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "Use 'hcs' for heteroscedastic standard deviations",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + ar1(1 + x | id, hom = FALSE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "Use 'har' for heteroscedastic standard deviations",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id, hom = NA),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "'hom' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + hcs(1 + x | id, hom = TRUE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "already heteroscedastic",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + us(1 + x | id, hom = TRUE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "Homogeneous random-effect standard deviations are not supported yet",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula reconstructs factor random-slope designs", {

  df <- data.frame(
    f = factor(c("a", "b", "c", "a"), levels = c("a", "b", "c")),
    id = factor(c("g1", "g1", "g2", "g2"), levels = c("g1", "g2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior_factor("gamma", list(2, 2), contrast = "treatment"),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  posterior <- matrix(
    c(1, 10, 20, 30, 40),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[1,2]",
      "mu__xREx__id_xRE_COEFx[2,1]",
      "mu__xREx__id_xRE_COEFx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = data.frame(
      f = c("c", "b"),
      id = factor(c("g1", "g2"), levels = c("g2", "g1"))
    ),
    prior_list = formula_result$prior_list,
    formula_target = "conditional"
  )

  expect_equal(unname(drop(prediction)), c(21, 31))
})

test_that("factor random slopes compile identically from data frames and tibbles", {

  skip_if_not_installed("tibble")

  df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2"), each = 3), levels = c("g1", "g2"))
  )
  compile_formula <- function(data){
    JAGS_formula(
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      formula_scale = TRUE,
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      )
    )
  }

  data_frame_result <- compile_formula(df)
  tibble_result <- compile_formula(tibble::as_tibble(df))
  data_frame_term <- data_frame_result$formula_design$random_effects[[1L]]
  tibble_term <- tibble_result$formula_design$random_effects[[1L]]

  expect_identical(tibble_term$model_terms_type, data_frame_term$model_terms_type)
  expect_identical(tibble_term$model_terms_type[["f"]], "factor")
  expect_identical(tibble_term$sd_parameter_names, data_frame_term$sd_parameter_names)
  expect_equal(
    tibble_result$data$mu__xREx__id_xRE_DATAx,
    data_frame_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_identical(names(tibble_result$prior_list), names(data_frame_result$prior_list))
  expect_length(tibble_result$formula_scale, 0L)
})

test_that("random-effect factor level names read tibble columns as vectors", {

  skip_if_not_installed("tibble")

  df <- data.frame(
    study = factor(rep(paste0("s", 1:6), each = 12)),
    t     = factor(rep(c("t1", "t2", "t3", "t4"), 18)),
    f     = factor(rep(c("x", "y", "z"), 24)),
    g     = factor(rep(c("p", "q"), each = 36))
  )
  tbl <- tibble::as_tibble(df)
  leaf_terms <- function(formula, data){
    compiled <- JAGS_formula(
      formula, "mu", data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("gamma", list(2, 2)))
    )
    unname(compiled$formula_design$random_effects[[1L]]$sd_leaves$leaf_terms)
  }

  expect_identical(
    leaf_terms(~ 1 + hcs(t | study), df),
    c("t[t1]", "t[t2]", "t[t3]", "t[t4]")
  )
  expect_identical(leaf_terms(~ 1 + us(1 + f | study), df), c("intercept", "f[y]", "f[z]"))
  for(formula in list(
    ~ 1 + hcs(t | study), ~ 1 + har(t | study), ~ 1 + us(1 + f | study),
    ~ 1 + us(1 + f * g | study)
  )){
    expect_identical(leaf_terms(formula, tbl), leaf_terms(formula, df), info = deparse(formula))
  }
})

test_that("prior_random maps to explicitly named random-effect blocks", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b")),
    drug = factor(c("d1", "d1", "d2", "d2", "d3", "d3"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))

  result <- JAGS_formula(
    formula = ~ 1 + x +
      random(1 + x | id, name = "study") +
      diag(0 + x | drug, name = "drug_slope"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      study = random_block(sd = sd_prior, cor = prior_lkj(eta = 2)),
      drug_slope = random_block(sd = sd_prior)
    )
  )

  expect_equal(
    names(result$prior_list),
    c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__study_intercept",
      "mu__xREx__study_x",
      "mu__xREx__drug_slope_x"
    )
  )
  expect_equal(
    names(result$formula_design$jags_data_names),
    c("x", "__xREx__study", "__xREx__drug_slope")
  )
  expect_equal(result$formula_design$random_effects[[1]]$block_name, "study")
  expect_equal(result$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(result$formula_design$random_effects[[2]]$block_name, "drug_slope")
  expect_equal(result$formula_design$random_effects[[2]]$structure, "diag")
  expect_equal(result$jags_modules, "BayesTools")
  expect_match(result$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)
  expect_false("backend" %in% names(result$formula_design$random_effects[[1]]$correlation))
  expect_false(any(grepl("dbeta", result$formula_syntax, fixed = TRUE)))

  inherited_covariance <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      cor = prior_lkj(eta = 3),
      study = random_block(sd = sd_prior)
    )
  )
  expect_equal(inherited_covariance$jags_modules, "BayesTools")
  expect_match(inherited_covariance$formula_syntax, "_lkj_alpha[1] <- 3", fixed = TRUE)
  expect_match(inherited_covariance$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)

  eta_override <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      cor = prior_lkj(eta = 3),
      study = random_block(
        sd = sd_prior,
        covariance = random_covariance(eta = 4)
      )
    )
  )
  expect_equal(eta_override$jags_modules, "BayesTools")
  expect_match(eta_override$formula_syntax, "_lkj_alpha[1] <- 4", fixed = TRUE)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + random(x | id, name = "study", covariance = "cs"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        covariance = random_covariance(cor = prior("normal", list(0, 0.5))),
        study = random_block(sd = sd_prior, covariance = random_covariance(eta = 2))
      )
    ),
    "uses a scalar correlation prior",
    fixed = TRUE
  )

  mixed_covariance <- JAGS_formula(
    formula = ~ 1 + x +
      random(1 + x | id, name = "study") +
      random(0 + x | drug, name = "drug_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      sd = sd_prior,
      cor = prior_lkj(eta = 3),
      drug_slope = random_block(covariance = random_covariance(structure = "diag"))
    )
  )
  expect_equal(mixed_covariance$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(mixed_covariance$formula_design$random_effects[[2]]$structure, "diag")
  expect_equal(mixed_covariance$jags_modules, "BayesTools")
  expect_match(mixed_covariance$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)
  expect_false(any(grepl("drug_slope_xRE_CORx", mixed_covariance$formula_syntax, fixed = TRUE)))

  wrapper_covariance <- JAGS_formula(
    formula = ~ 1 + x + random(0 + x | id, name = "id_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id_slope = random_block(sd = sd_prior)
    )
  )
  expect_equal(wrapper_covariance$formula_design$random_effects[[1]]$block_name, "id_slope")
  expect_equal(wrapper_covariance$formula_design$random_effects[[1]]$structure, "diag")
  expect_equal(names(wrapper_covariance$prior_list), c("mu_intercept", "mu_x", "mu__xREx__id_slope_x"))

  unnamed_blocks <- JAGS_formula(
    formula = ~ 1 + x + (1 | id) + (0 + x | drug),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(sd = sd_prior)
  )
  expect_equal(
    vapply(unnamed_blocks$formula_design$random_effects, `[[`, character(1), "block_name"),
    c("id", "drug")
  )
  expect_equal(
    names(unnamed_blocks$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__drug_x")
  )
  expect_no_error(
    JAGS_formula(
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(sd = sd_prior)
    )
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(study = random_block(sd = sd_prior))
    ),
    "block override names were not found",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        allocation = random_variance_allocation(name = "allocation",
          terms = c("study", "drug"),
          sd = sd_prior
        )
      )
    ),
    "Variance allocation priors require formula random-effect terms",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x +
        random(1 | id, name = "study") +
        random(0 + x | drug, name = "study"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(study = random_block(sd = sd_prior))
    ),
    "Random-effect block names must be unique",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x +
        random(1 | id, name = "study-id") +
        random(0 + x | drug, name = "study_id"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(study_id = random_block(sd = sd_prior))
    ),
    "Random-effect block names must be unique",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + random(1 | id, name = "!!!"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(sd = sd_prior)
    ),
    "must contain at least one letter",
    fixed = TRUE
  )
})

test_that("nested grouping formulas expand to evaluable random-effect blocks", {

  df <- data.frame(
    x = 1:6,
    study = factor(c("s1", "s1", "s1", "s2", "s2", "s2")),
    paper = factor(c("p1", "p1", "p2", "p1", "p2", "p2"))
  )

  result <- JAGS_formula(
    formula = ~ x +
      random(1 | paper:study, name = "paper_study", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      paper_study = random_block(sd = prior("gamma", list(2, 2))),
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  expect_equal(
    names(result$formula_design$jags_data_names),
    c("x", "__xREx__paper_study", "__xREx__study")
  )
  expect_equal(
    names(result$prior_list),
    c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__paper_study_intercept",
      "mu__xREx__study_intercept"
    )
  )
  expect_equal(result$formula_design$random_effects[[1]]$group_label, "paper:study")
  expect_equal(result$formula_design$random_effects[[1]]$block_name, "paper_study")
  expect_equal(result$formula_design$random_effects[[2]]$group_label, "study")
  expect_equal(result$formula_design$random_effects[[1]]$n_groups, 4L)
  expect_equal(result$formula_design$random_effects[[2]]$n_groups, 2L)

  char_df <- transform(
    df,
    study = as.character(study),
    paper = as.character(paper)
  )
  char_result <- JAGS_formula(
    formula = ~ x +
      random(1 | paper:study, name = "paper_study", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = char_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      paper_study = random_block(sd = prior("gamma", list(2, 2))),
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expect_equal(char_result$formula_design$random_effects[[1]]$group_label, "paper:study")
  expect_equal(char_result$formula_design$random_effects[[1]]$n_groups, 4L)
  expect_equal(char_result$formula_design$random_effects[[2]]$n_groups, 2L)
})

test_that("formula designs keep no environment of the function that built them", {

  # JAGS_fit() stores the design, standardization, and priors of every
  # formula. Built in a function that holds a large local object, none of
  # their formulas may keep that function's frame as its environment: a saved
  # fit would serialize it. (The returned fixed-effect 'formula' is not stored
  # and keeps the caller's environment.)
  formula_output <- function(local_object_size){
    local_object <- stats::runif(local_object_size)
    data <- data.frame(
      x = c(-1.2, 0.4, 0.9, -0.3, 1.5, -0.8, 0.2, 1.1),
      t = factor(rep(c("t1", "t2"), 4)),
      g = factor(rep(c("g1", "g2"), each = 4)),
      s = factor(rep(c("s1", "s2", "s3", "s4"), each = 2))
    )
    output <- JAGS_formula(
      formula = ~ 1 + x + us(1 + x | g) + ar1(t | s),
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      formula_scale = list(x = TRUE),
      prior_random = prior_random(
        g = random_block(sd = prior("gamma", list(2, 2)), cor = prior_lkj(eta = 1)),
        s = random_block(sd = prior("gamma", list(2, 2)), cor = prior("normal", list(0, 0.5)))
      )
    )
    output[c("formula_design", "formula_scale", "prior_list")]
  }
  # The size of 'small' is taken before 'large' exists: a captured frame
  # reaches this test's environment, and through it both objects.
  small <- formula_output(0)
  small_size <- length(serialize(small, NULL))
  large <- formula_output(1e6)

  expect_identical(stored_environment_paths(large, "stored"), character())
  expect_identical(
    vapply(large$formula_design$random_effects, function(random_term){
      identical(environment(random_term$term_formula), baseenv())
    }, logical(1)),
    c(TRUE, TRUE)
  )
  expect_identical(large, small)
  expect_identical(length(serialize(large, NULL)), small_size)
})

test_that("diag random-effect syntax handles slope-only and factor-slope designs", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )

  slope_result <- JAGS_formula(
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  expect_equal(.remove_random_effects(~ z + diag(0 + x | id)), ~ z, ignore_formula_env = TRUE)
  expect_equal(slope_result$formula_design$random_effects[[1]]$term_formula, ~ 0 + x, ignore_formula_env = TRUE)
  expect_equal(slope_result$formula_design$random_effects[[1]]$structure, "diag")
  expect_true(isTRUE(slope_result$formula_design$random_effects[[1]]$independent))
  expect_equal(slope_result$formula_design$random_effects[[1]]$n_columns, 1L)
  expect_equal(colnames(slope_result$data$mu__xREx__id_xRE_DATAx), "x")
  expect_equal(unname(slope_result$data$mu__xREx__id_xRE_DATAx[, "x"]), df$x)
  expect_equal(
    names(slope_result$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_x")
  )

  factor_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 2))
  )

  factor_result <- JAGS_formula(
    formula = ~ 1 + f + diag(0 + f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  expected_factor_df <- factor_df
  stats::contrasts(expected_factor_df$f) <- "contr.treatment"
  expected_factor_design <- stats::model.matrix(~ f, data = expected_factor_df)[, -1, drop = FALSE]

  expect_equal(colnames(factor_result$data$mu__xREx__id_xRE_DATAx), c("fb", "fc"))
  expect_equal(
    as.vector(factor_result$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_factor_design)
  )
  expect_equal(factor_result$formula_design$random_effects[[1]]$model_terms_type[["f"]], "factor")
  expect_equal(factor_result$formula_design$random_effects[[1]]$n_columns, 2L)
  expect_equal(
    names(factor_result$prior_list),
    c("mu_intercept", "mu_f", "mu__xREx__id_f")
  )
  expect_true(is.prior.treatment(factor_result$prior_list$mu__xREx__id_f))
  expect_equal(attr(factor_result$prior_list$mu__xREx__id_f, "random_factor"), "id")

  ordered_df <- factor_df
  ordered_df$f <- ordered(ordered_df$f, levels = c("a", "b", "c"))
  old_contrasts <- options("contrasts")
  on.exit(options(old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  ordered_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = ordered_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expected_ordered_df <- ordered_df
  stats::contrasts(expected_ordered_df$f) <- "contr.treatment"
  expected_ordered <- stats::model.matrix(~ f, data = expected_ordered_df)[, -1, drop = FALSE]
  expect_equal(colnames(ordered_result$data$mu__xREx__id_xRE_DATAx), c("fb", "fc"))
  expect_equal(
    as.vector(ordered_result$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_ordered)
  )
  expect_true(is.prior.treatment(ordered_result$prior_list$mu__xREx__id_f))

  random_only_orthonormal <- NULL
  expect_no_warning(
    random_only_orthonormal <- JAGS_formula(
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)), contrasts = c(f = "orthonormal"))
      )
    )
  )
  orthonormal_df <- factor_df
  stats::contrasts(orthonormal_df$f) <- "contr.orthonormal"
  expected_orthonormal <- stats::model.matrix(~ f, data = orthonormal_df)[, -1, drop = FALSE]

  expect_equal(
    as.vector(random_only_orthonormal$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_orthonormal),
    tolerance = 1e-12
  )
  expect_equal(unname(random_only_orthonormal$formula_design$random_effects[[1]]$contrast_matrices$f),
               unname(contr.orthonormal(levels(factor_df$f))), tolerance = 0)

  random_only_meandif <- NULL
  expect_no_warning(
    random_only_meandif <- JAGS_formula(
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)), contrasts = c(f = "meandif"))
      )
    )
  )
  meandif_df <- factor_df
  stats::contrasts(meandif_df$f) <- "contr.meandif"
  expected_meandif <- stats::model.matrix(~ f, data = meandif_df)[, -1, drop = FALSE]

  expect_equal(
    as.vector(random_only_meandif$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_meandif),
    tolerance = 1e-12
  )
  expect_equal(unname(random_only_meandif$formula_design$random_effects[[1]]$contrast_matrices$f),
               unname(contr.meandif(levels(factor_df$f))), tolerance = 0)

  interaction_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    g = factor(rep(c("u", "v"), each = 6), levels = c("u", "v")),
    id = factor(rep(c("s1", "s2", "s3", "s4"), each = 3))
  )
  interaction_result <- JAGS_formula(
    formula = ~ 1 + f * g + diag(0 + f:g | id),
    parameter = "mu",
    data = interaction_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment"),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expected_interaction_df <- interaction_df
  stats::contrasts(expected_interaction_df$f) <- "contr.treatment"
  stats::contrasts(expected_interaction_df$g) <- "contr.treatment"
  expected_interaction <- stats::model.matrix(~ f:g, data = expected_interaction_df)[, -1, drop = FALSE]
  colnames(expected_interaction) <- gsub(":", "__xXx__", colnames(expected_interaction), fixed = TRUE)

  expect_equal(
    colnames(interaction_result$data$mu__xREx__id_xRE_DATAx),
    colnames(expected_interaction)
  )
  expect_equal(
    as.vector(interaction_result$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_interaction),
    tolerance = 1e-12
  )
  expect_equal(
    attr(interaction_result$prior_list$mu__xREx__id_f__xXx__g, "level_names"),
    list(f = c("a", "b", "c"), g = c("u", "v"))
  )

  mixed_contrast_df <- interaction_df
  stats::contrasts(mixed_contrast_df$f) <- "contr.treatment"
  stats::contrasts(mixed_contrast_df$g) <- "contr.orthonormal"
  mixed_contrast_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f:g | id),
    parameter = "mu",
    data = mixed_contrast_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expected_mixed <- stats::model.matrix(
    ~ f:g,
    data = mixed_contrast_df
  )[, -1, drop = FALSE]
  mixed_term <- mixed_contrast_result$formula_design$random_effects[[1L]]
  expect_equal(
    as.vector(mixed_term$model_matrix),
    as.vector(expected_mixed),
    tolerance = 1e-12
  )
  expect_equal(
    mixed_term$contrast_matrices$f,
    stats::contrasts(mixed_contrast_df$f)
  )
  expect_equal(
    mixed_term$contrast_matrices$g,
    stats::contrasts(mixed_contrast_df$g)
  )
})

test_that("fixed and random blocks own independent concrete factor bases", {

  data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4)),
    id_treatment = factor(rep(c("s1", "s2", "s3", "s4"), each = 3)),
    id_meandif = factor(rep(c("q1", "q2", "q3"), 4))
  )
  inherited <- JAGS_formula(
    formula = ~ 0 + f + diag(0 + f | id_treatment),
    parameter = "mu",
    data = data,
    prior_list = list(
      f = prior_factor(
        "normal",
        list(0, 1),
        contrast = "independent"
      )
    ),
    prior_random = prior_random(
      id_treatment = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  inherited_term <- inherited$formula_design$random_effects[[1L]]

  expect_equal(
    inherited$formula_design$contrast_matrices$f,
    contr.independent(levels(data$f))
  )
  expect_equal(
    inherited_term$contrast_matrices$f,
    contr.independent(levels(data$f))
  )
  expect_equal(inherited_term$n_columns, nlevels(data$f))

  result <- JAGS_formula(
    formula = ~ f +
      diag(0 + f | id_treatment) +
      diag(0 + f | id_meandif),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor(
        "mnormal",
        list(0, 1),
        contrast = "orthonormal"
      )
    ),
    prior_random = prior_random(
      id_treatment = random_block(
        sd = prior("gamma", list(2, 2)),
        contrasts = c(f = "treatment")
      ),
      id_meandif = random_block(
        sd = prior("gamma", list(2, 2)),
        contrasts = c(f = "meandif")
      )
    )
  )

  fixed_matrix <- contr.orthonormal(levels(data$f))
  treatment_matrix <- stats::contr.treatment(levels(data$f))
  meandif_matrix <- contr.meandif(levels(data$f))
  treatment_term <- result$formula_design$random_effects[[1L]]
  meandif_term <- result$formula_design$random_effects[[2L]]

  expect_identical(result$formula_design$schema_version, 7L)
  expect_equal(result$formula_design$contrast_matrices$f, fixed_matrix)
  expect_equal(treatment_term$contrast_matrices$f, treatment_matrix)
  expect_equal(meandif_term$contrast_matrices$f, meandif_matrix)
  expect_identical(treatment_term$contrast_owner, "random_block")
  expect_identical(meandif_term$contrast_owner, "random_block")
  expect_false(isTRUE(all.equal(fixed_matrix, treatment_matrix)))
  expect_false(isTRUE(all.equal(fixed_matrix, meandif_matrix)))
  # The mean-difference columns are one exchangeable unit with one SD.
  expect_identical(
    meandif_term$sd_parameter_names,
    rep("mu__xREx__id_meandif_f", 2L)
  )
  expect_identical(
    meandif_term$sd_leaves$leaf_terms_by_column,
    rep("f", 2L)
  )

  treatment_prediction <- .bt_random_effect_prediction_data(
    treatment_term,
    data
  )
  meandif_prediction <- .bt_random_effect_prediction_data(
    meandif_term,
    data
  )
  expect_equal(
    treatment_prediction$model_matrix,
    treatment_term$model_matrix
  )
  expect_equal(
    meandif_prediction$model_matrix,
    meandif_term$model_matrix
  )

  missing_basis <- treatment_term
  missing_basis$contrast_matrices <- NULL
  expect_error(
    .bt_random_effect_prediction_data(missing_basis, data),
    "missing its owner-scoped concrete factor basis",
    fixed = TRUE
  )
})

test_that("random-block contrast specifications validate their design scope", {

  expect_error(
    random_block(contrasts = "treatment"),
    "must be a named character vector",
    fixed = TRUE
  )
  expect_error(
    random_block(contrasts = c(f = "sum")),
    "Unknown random-block contrast value",
    fixed = TRUE
  )

  data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2)),
    id = factor(rep(c("s1", "s2"), each = 3))
  )
  expect_error(
    JAGS_formula(
      ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          contrasts = c(unknown = "treatment")
        )
      )
    ),
    "contrast overrides reference predictors outside this design",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ 1 + cs(f | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          cor = prior("normal", list(0, 0.5)),
          contrasts = c(f = "orthonormal")
        )
      )
    ),
    "level basis is defined by the covariance structure",
    fixed = TRUE
  )

  # A random-only factor carrying a contrasts matrix or an unsupported contrast
  # name is refused naming the factor, instead of failing inside vapply().
  compile_with_attribute <- function(contrasts_attribute, formula, block){
    d <- data
    if(!is.null(contrasts_attribute)){
      stats::contrasts(d$f) <- contrasts_attribute
    }
    JAGS_formula(
      formula,
      parameter = "mu",
      data = d,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = block)
    )
  }
  gamma_block <- random_block(sd = prior("gamma", list(2, 2)))
  for(formula in list(~ 1 + diag(0 + f | id), ~ 1 + us(1 + f | id), ~ 1 + (1 + f | id))){
    for(matrix_contrasts in list(stats::contr.treatment(3), stats::contr.sum(3))){
      expect_error(
        compile_with_attribute(matrix_contrasts, formula, gamma_block),
        "Random-effect factor 'f' uses a contrasts matrix, which random effects do not support.",
        fixed = TRUE
      )
    }
    expect_error(
      compile_with_attribute("contr.sum", formula, gamma_block),
      "Random-effect factor 'f' uses the contrasts 'contr.sum', which random effects do not support.",
      fixed = TRUE
    )
    expect_error(
      compile_with_attribute(stats::contr.sum(3), formula, gamma_block),
      "random_block(contrasts = c(f = \"treatment\"))",
      fixed = TRUE
    )
    # F318 (2): a classed refusal naming the factor
    for(contrasts_attribute in list(stats::contr.sum(3), "contr.sum")){
      condition <- tryCatch(
        compile_with_attribute(contrasts_attribute, formula, gamma_block),
        error = identity
      )
      expect_s3_class(condition, "BayesTools_formula_transform_unavailable")
      expect_identical(condition$reason, "unsupported_random_factor_contrast")
      expect_identical(condition$predictor, "f")
    }
    # the message lists the accepted aliases, which compile
    expect_error(
      compile_with_attribute(stats::contr.sum(3), formula, gamma_block),
      paste0("\"treatment\" (also \"dummy\")", ", \"independent\", ",
             "\"orthonormal\", \"meandif\", \"cumulative\" and \"cumulative_levels\", ",
             "each also with the 'contr.' prefix; a factor already carrying the ",
             "'contr.mixed' coding is kept as it is."),
      fixed = TRUE
    )
    for(alias in c("dummy", "contr.treatment", "contr.orthonormal")){
      expect_no_error(compile_with_attribute(
        stats::contr.sum(3), formula,
        random_block(sd = prior("gamma", list(2, 2)), contrasts = c(f = alias))
      ))
    }
    # Without the attribute, with a supported name or with a block override the
    # same design compiles.
    expect_no_error(compile_with_attribute(NULL, formula, gamma_block))
    expect_no_error(compile_with_attribute("contr.treatment", formula, gamma_block))
    expect_no_error(compile_with_attribute(
      stats::contr.sum(3), formula,
      random_block(sd = prior("gamma", list(2, 2)), contrasts = c(f = "treatment"))
    ))
  }
  # Structures that own the level basis never read the factor's contrasts.
  cs_block <- random_block(sd = prior("gamma", list(2, 2)), cor = prior("normal", list(0, 0.5)))
  for(contrasts_attribute in list(stats::contr.sum(3), "contr.sum")){
    expect_no_error(compile_with_attribute(contrasts_attribute, ~ 1 + cs(f | id), cs_block))
  }
})

test_that("structured random-effect terms use level-indexed factor columns and scalar rho priors", {

  factor_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    g = factor(rep(c("u", "v"), 6), levels = c("u", "v")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )
  continuous_df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 0.5),
    id = factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"))
  )
  car_df <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))

  large_factor_df <- data.frame(
    f = factor(sprintf("level_%03d", seq_len(113L))),
    id = factor(rep(sprintf("group_%02d", seq_len(17L)), length.out = 113L))
  )
  large_result <- JAGS_formula(
      formula = ~ 1 + cs(f | id),
      parameter = "mu",
      data = large_factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    )
  large_term <- large_result$formula_design$random_effects[[1]]
  expect_s3_class(
    large_term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_equal(large_term$latent_layout$n_local, 113L)
  expect_equal(
    names(large_result$data),
    c("N_mu", "mu__xREx__id_xRE_MAPx", "mu__xREx__id_xRE_COLx")
  )
  expect_false(grepl("_xRE_CORx_L", large_result$formula_syntax, fixed = TRUE))
  expect_false(grepl("_xRE_CORx_R", large_result$formula_syntax, fixed = TRUE))

  coefficient_result <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = large_factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior("normal", list(0, 0.5)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE)
      )
    )
  )
  coefficient_term <- coefficient_result$formula_design$random_effects[[1]]
  expect_null(coefficient_term$latent_layout)
  expect_true("mu__xREx__id_xRE_COEFx" %in% coefficient_result$add_parameters)
  expect_match(
    coefficient_result$formula_syntax,
    "mu__xREx__id_xRE_Zx[i,j] ~ dnorm(0, 1)",
    fixed = TRUE
  )

  old_multiplier <- getOption("BayesTools.random_effects_complexity_multiplier")
  on.exit(options(BayesTools.random_effects_complexity_multiplier = old_multiplier), add = TRUE)
  options(BayesTools.random_effects_complexity_multiplier = 0)
  expect_error(
    BayesTools:::.bt_random_effect_check_dense_complexity(
      random_term = list(block_name = "id"),
      structure = "cs",
      n_groups = 17L,
      n_columns = 113L,
      n_rows = 113L,
      monitor_policy = random_monitor()
    ),
    "must be a positive numeric scalar",
    fixed = TRUE
  )
  options(BayesTools.random_effects_complexity_multiplier = old_multiplier)

  options(BayesTools.random_effects_complexity_multiplier = Inf)
  expect_null(
    BayesTools:::.bt_random_effect_check_dense_complexity(
      random_term = list(block_name = "id"),
      structure = "cs",
      n_groups = 17L,
      n_columns = 113L,
      n_rows = 113L,
      monitor_policy = random_monitor()
    )
  )
  options(BayesTools.random_effects_complexity_multiplier = old_multiplier)

  default_cs <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(names(default_cs$prior_list), c(
    "mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho"
  ))
  expect_equal(default_cs$prior_list$mu__xREx__id_rho$distribution, "uniform")
  expect_equal(
    default_cs$prior_list$mu__xREx__id_rho$parameters,
    list(a = -0.5, b = 1)
  )
  expect_equal(
    default_cs$formula_design$random_effects[[1]]$correlation$rho_scale,
    "rho"
  )

  fisher_cs <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(names(fisher_cs$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_equal(fisher_cs$prior_list$mu__xREx__id_rho_z$distribution, "normal")
  expect_equal(fisher_cs$prior_list$mu__xREx__id_rho_z$parameters, list(mean = 0, sd = 0.5))
  expect_true(fisher_cs$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(fisher_cs$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(fisher_cs$data$mu__xREx__id_xRE_DATAx)
  expect_equal(
    fisher_cs$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 3)
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + cs(1 | id),
      parameter = "mu",
      data = continuous_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "does not support explicit '1', '0', or '-1' terms",
    fixed = TRUE
  )

  id_result <- JAGS_formula(
    formula = ~ 1 + x + id(1 + x | id),
    parameter = "mu",
    data = continuous_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(names(id_result$prior_list), c("mu_intercept", "mu_x", "mu__xREx__id_sd"))
  expect_true(id_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    id_result$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 2)
  )

  diag_heterogeneous <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id, hom = FALSE),
    parameter = "mu",
    data = continuous_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(names(diag_heterogeneous$prior_list), c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__id_x"))
  expect_false(diag_heterogeneous$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    diag_heterogeneous$formula_design$random_effects[[1]]$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_x")
  )

  hom_env <- new.env(parent = globalenv())
  hom_env$hom_flag <- FALSE
  hom_formula <- stats::as.formula("~ 1 + x + diag(1 + x | id, hom = hom_flag)", env = hom_env)
  hom_from_env <- JAGS_formula(
    formula = hom_formula,
    parameter = "mu",
    data = continuous_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_false(hom_from_env$formula_design$random_effects[[1]]$homogeneous_sd)

  hcs_index <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    names(hcs_index$prior_list),
    c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho_z")
  )
  expect_false(hcs_index$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    hcs_index$formula_design$random_effects[[1]]$sd_parameter_names,
    paste0("mu__xREx__id_f[", 1:3, "]")
  )

  wrapper_hcs_result <- JAGS_formula(
    formula = ~ 1 + random(f | id, covariance = "cs", hom = FALSE),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(wrapper_hcs_result$formula_design$random_effects[[1]]$structure, "hcs")
  expect_false(wrapper_hcs_result$formula_design$random_effects[[1]]$homogeneous_sd)

  hcs_composite <- JAGS_formula(
    formula = ~ 1 + hcs(f + g | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$structured_index$variables, c("f", "g"))
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$structured_index$label, "f:g")
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$n_columns, 6L)
  expect_null(hcs_composite$data$mu__xREx__id_xRE_DATAx)
  expect_equal(
    hcs_composite$formula_design$random_effects[[1]]$sd_leaves$leaf_terms_by_column,
    paste0("f_g[", c("a.u", "a.v", "b.u", "b.v", "c.u", "c.v"), "]")
  )
  expect_equal(
    BayesTools:::.bt_random_effect_summary_sd_components(
      hcs_composite$formula_design$random_effects[[1]],
      paste0("mu__xREx__id_f_g[", 1:6, "]")
    ),
    paste0("f:g[", c("a.u", "a.v", "b.u", "b.v", "c.u", "c.v"), "]")
  )
  hcs_composite_raw <- .parameter_catalog_random_summary_samples(
    model_samples = matrix(
      seq_len(7),
      nrow = 1,
      dimnames = list(
        NULL,
        c(paste0("mu__xREx__id_f_g[", 1:6, "]"), "mu__xREx__id_rho")
      )
    ),
    prior_list = hcs_composite$prior_list,
    formula_design = list(mu = hcs_composite$formula_design),
    mode = "raw"
  )
  hcs_composite_raw_coordinates <- build_test_parameter_coordinates(
    columns = colnames(hcs_composite_raw$model_samples),
    prior_list = hcs_composite_raw$prior_list,
    formula_design = list(mu = hcs_composite$formula_design)
  )
  hcs_composite_raw_display <- hcs_composite_raw_coordinates$display_label
  expect_true("(mu) id: sd(f:g[a.u])" %in% hcs_composite_raw_display)
  expect_true("(mu) id: cor" %in% hcs_composite_raw_display)
  expect_false(any(grepl("f_g", hcs_composite_raw_display, fixed = TRUE)))
  expect_false(any(grepl("sd((mu)", hcs_composite_raw_display, fixed = TRUE)))

  hcs_composite_allocation <- JAGS_formula(
    formula = ~ 1 + hcs(f + g | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      random_variance_allocation(
        name = "leaf_alloc",
        terms = "id",
        target = "sd_component",
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = rep(1, 6)))
      ),
      id = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$parent_factors <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing 'allocation\\$parent_factors'"
  )

  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$n_targets <- 6.5
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$n_targets'"
  )

  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column <-
    malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column[-1L]
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "leaf_index_by_column"
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$n_targets <- 1L
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column <-
    rep(1L, malformed_random_term$n_columns)
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_names <- "only"
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_terms <- "only"
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$n_targets'"
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column <-
    rep(1L, malformed_random_term$n_columns)
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "must cover every target",
    fixed = TRUE
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_names[[2L]] <-
    malformed_random_term$sd_binding$allocations[[1L]]$leaf_names[[1L]]
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$leaf_names'"
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_terms <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$leaf_terms'"
  )

  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(ar1_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(names(ar1_result$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_true(ar1_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(ar1_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_equal(
    ar1_result$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 3)
  )
  expect_match(
    ar1_result$formula_syntax,
    "mu__xREx__id_xRE_UNIT_COEFx[g,i - 1]",
    fixed = TRUE
  )

  ar1_implicit_levels <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    ar1_implicit_levels$data$mu__xREx__id_xRE_DATAx,
    ar1_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_equal(
    ar1_implicit_levels$formula_design$random_effects[[1]]$term_formula,
    ~ f - 1,
    ignore_formula_env = TRUE
  )

  wrapper_har_result <- JAGS_formula(
    formula = ~ 1 + re(f | id, covariance = "ar1", hom = FALSE),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(wrapper_har_result$formula_design$random_effects[[1]]$structure, "har")
  expect_false(wrapper_har_result$formula_design$random_effects[[1]]$homogeneous_sd)

  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(hcs_result$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(hcs_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(names(hcs_result$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho_z"))
  expect_false(hcs_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    hcs_result$formula_design$random_effects[[1]]$sd_parameter_names,
    paste0("mu__xREx__id_f[", 1:3, "]")
  )
  expect_equal(hcs_result$prior_list$mu__xREx__id_rho_z$truncation$lower, atanh(-1 / 2), tolerance = 1e-12)
  expect_s3_class(hcs_result$prior_list$mu__xREx__id_f, "prior.independent")
  expect_equal(attr(hcs_result$prior_list$mu__xREx__id_f, "levels"), 3L)
  expect_equal(attr(hcs_result$prior_list$mu__xREx__id_f, "level_names"), c("a", "b", "c"))
  expect_equal(
    hcs_result$formula_design$random_effects[[1]]$sd_leaves$leaf_terms_by_column,
    c("f[a]", "f[b]", "f[c]")
  )
  hcs_implicit_levels <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    hcs_implicit_levels$data$mu__xREx__id_xRE_DATAx,
    hcs_result$data$mu__xREx__id_xRE_DATAx
  )

  har_implicit_levels <- JAGS_formula(
    formula = ~ 1 + har(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(har_implicit_levels$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(har_implicit_levels$data$mu__xREx__id_xRE_DATAx)
  expect_false(har_implicit_levels$formula_design$random_effects[[1]]$homogeneous_sd)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + hcs(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
          cor = prior("normal", list(0, 0.5))
        )
      )
    ),
    "Unsupported multivariate random-effect SD prior",
    fixed = TRUE
  )
  factor_sd <- prior_factor("gamma", list(2, 2), contrast = "treatment")
  expect_true(is.prior.factor(factor_sd))
  expect_false(is.prior.vector(factor_sd))
  expect_no_error(random_block(sd = factor_sd))
  expect_error(
    JAGS_formula(
      formula = ~ 1 + hcs(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = factor_sd, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "Random-effect SD prior for 'f' must be an ordinary scalar prior.",
    fixed = TRUE
  )
  hcs_mixture_sd <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior_mixture(list(prior("gamma", list(2, 2)), prior("point", list(0)))),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  expect_s3_class(hcs_mixture_sd$prior_list$mu__xREx__id_f, "prior.factor_mixture")

  raw_result <- NULL
  expect_warning(
    raw_result <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(
            cor = prior("normal", list(0, 0.5)),
            cor_scale = "cor"
          )
        )
      )
    ),
    "valid correlation range",
    fixed = TRUE
  )

  expect_equal(names(raw_result$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_equal(raw_result$prior_list$mu__xREx__id_rho$truncation, list(lower = -1, upper = 1))
  expect_equal(
    raw_result$formula_design$random_effects[[1]]$sd_parameter_names,
    paste0("mu__xREx__id_f[", 1:3, "]")
  )
  expect_false(grepl("rho_z", raw_result$formula_syntax, fixed = TRUE))
  expect_match(
    raw_result$formula_syntax,
    "mu__xREx__id_rho * mu__xREx__id_xRE_UNIT_COEFx[g,i - 1]",
    fixed = TRUE
  )

  raw_inherited <- NULL
  expect_no_warning(
    raw_inherited <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "cor"
        ),
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "har")
        )
      )
    )
  )
  expect_equal(names(raw_inherited$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_identical(raw_inherited$prior_list$mu__xREx__id_rho$distribution, "uniform")
  expect_equal(raw_inherited$prior_list$mu__xREx__id_rho$parameters, list(a = -1, b = 1))
  expect_false(grepl("rho_z", raw_inherited$formula_syntax, fixed = TRUE))

  raw_block_override <- NULL
  expect_no_warning(
    raw_block_override <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(cor = prior("normal", list(0, 0.5))),
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "har", cor_scale = "cor")
        )
      )
    )
  )
  expect_equal(names(raw_block_override$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_identical(raw_block_override$prior_list$mu__xREx__id_rho$distribution, "uniform")
  expect_equal(raw_block_override$prior_list$mu__xREx__id_rho$parameters, list(a = -1, b = 1))
  expect_false(grepl("rho_z", raw_block_override$formula_syntax, fixed = TRUE))

  raw_top_scale_block_rho <- NULL
  expect_no_warning(
    raw_top_scale_block_rho <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(cor_scale = "cor"),
        id = random_block(
          sd = sd_prior,
          cor = prior("normal", list(0, 0.5))
        )
      )
    )
  )
  expect_equal(names(raw_top_scale_block_rho$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho_z"))
  expect_match(raw_top_scale_block_rho$formula_syntax, "tanh(mu__xREx__id_rho_z)", fixed = TRUE)

  car_result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  car_term <- car_result$formula_design$random_effects[[1]]
  car_distance <- abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))

  expect_equal(car_term$column_names, c("time_0", "time_0p5", "time_2"))
  expect_null(car_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(names(car_result$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_true(car_term$homogeneous_sd)
  expect_equal(car_term$structure, "car")
  expect_equal(car_term$car$time_variable, "time")
  expect_equal(car_term$car$time_values, c(0, 0.5, 2))
  expect_null(car_term$car$distance_matrix)
  expect_equal(car_term$correlation$bounds, c(lower = 0, upper = 1))
  expect_null(car_term$correlation$distance_matrix)
  expect_equal(car_term$correlation$time_values, c(0, 0.5, 2))
  expect_equal(car_result$prior_list$mu__xREx__id_rho_z$truncation$lower, 0)
  expect_equal(car_result$prior_list$mu__xREx__id_rho_z$truncation$upper, Inf)
  expect_match(
    car_result$formula_syntax,
    "mu__xREx__id_rho <- tanh(mu__xREx__id_rho_z)",
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    "tanh(mu__xREx__id_rho_z)",
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_LOG_PHIX[2] <- ",
      "0.5 * log(mu__xREx__id_rho)"
    ),
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_LOG_PHIX[3] <- ",
      "1.5 * log(mu__xREx__id_rho)"
    ),
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_PHIX[2] <- ",
      "exp(mu__xREx__id_xRE_CAR_LOG_PHIX[2])"
    ),
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_INNOV_VARx[2] <- ",
      "pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[2], 1)"
    ),
    fixed = TRUE
  )
  expect_false(grepl(
    "pow(mu__xREx__id_rho",
    car_result$formula_syntax,
    fixed = TRUE
  ))
  expect_false(grepl(
    "1 - pow(mu__xREx__id_xRE_CAR_PHIX",
    car_result$formula_syntax,
    fixed = TRUE
  ))

  car_independent <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("point", list(location = 0)))
    )
  )
  car_independent_term <- car_independent$formula_design$random_effects[[1]]
  empty_posterior <- matrix(numeric(0), nrow = 2, ncol = 0)
  expect_equal(
    BayesTools:::.bt_random_effect_rho_draws(car_independent_term, empty_posterior),
    rep(
      BayesTools:::.bt_random_effect_representable_rho_bounds(
        car_independent_term$correlation$bounds,
        car_independent_term$structure
      )[["lower"]],
      2L
    )
  )
  expect_equal(
    bayestools_reference_random_effect_scalar_rho_support(
      empty_posterior,
      car_independent_term
    ),
    0
  )

  car_independent_raw <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("point", list(location = 0)),
          cor_scale = "cor"
        )
      )
    )
  )
  expect_equal(
    car_independent_raw$formula_design$random_effects[[1]]$correlation$sample_fixed,
    0
  )

  car_explicit_no_intercept <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    car_explicit_no_intercept$data$mu__xREx__id_xRE_DATAx,
    car_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_error(
    BayesTools:::.bt_random_effect_prediction_data(
      car_term,
      data.frame(time = 1, id = factor("a", levels = c("a", "b")))
    ),
    "New CAR time coordinate",
    fixed = TRUE
  )

  car_scale_request <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = TRUE,
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_length(car_scale_request$formula_scale, 0L)
  expect_identical(attr(car_scale_request$formula_scale, "unscale_design")$owner_scope, "compiler")
  expect_equal(car_scale_request$formula_design$random_effects[[1]]$car$time_values, c(0, 0.5, 2))

  fixed_car_scale <- JAGS_formula(
    formula = ~ 1 + time + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      time = prior("normal", list(0, 1))
    ),
    formula_scale = TRUE,
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(fixed_car_scale$formula_design$random_effects[[1]]$car$time_values, c(0, 0.5, 2))
  expect_null(fixed_car_scale$formula_design$random_effects[[1]]$car$distance_matrix)
  expect_equal(
    fixed_car_scale$data$mu__xREx__id_xRE_DATAx,
    car_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_equal(
    unname(fixed_car_scale$data$mu_data_time),
    unname(as.numeric(scale(car_df$time)))
  )
  expect_equal(names(fixed_car_scale$formula_scale), "mu_time")

  car_prediction_term <- fixed_car_scale$formula_design$random_effects[[1]]
  car_coefficient_names <- as.vector(BayesTools:::.bt_random_effect_coefficient_names(
    random_term = car_prediction_term,
    n_groups = length(car_prediction_term$group_levels),
    n_columns = car_prediction_term$n_columns
  ))
  car_prediction_posterior <- matrix(
    0,
    nrow = 1,
    ncol = 2L + length(car_coefficient_names),
    dimnames = list(NULL, c("mu_intercept", "mu_time", car_coefficient_names))
  )
  car_prediction_posterior[, "mu_intercept"] <- 1
  car_prediction_posterior[, "mu_time"] <- 0.25
  car_prediction_fit <- coda::mcmc(car_prediction_posterior)
  attr(car_prediction_fit, "formula_design") <- list(mu = fixed_car_scale$formula_design)
  attr(car_prediction_fit, "formula_scale") <- list(mu = fixed_car_scale$formula_scale)
  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      fit = car_prediction_fit,
      formula = ~ 1 + time + car(0 + time | id),
      parameter = "mu",
      data = car_df,
      prior_list = fixed_car_scale$prior_list,
      formula_target = "conditional"
    ))),
    1 + 0.25 * unname(as.numeric(scale(car_df$time))),
    tolerance = 1e-12
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = transform(car_df, time = c(0, NA, 2, 0, 0.5, 2)),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "must contain only finite values",
    fixed = TRUE
  )
  bad_ordered_car_df <- car_df
  bad_ordered_car_df$time <- ordered(
    c("early", "mid", "late", "early", "mid", "late"),
    levels = c("early", "mid", "late")
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = bad_ordered_car_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "numeric level labels",
    fixed = TRUE
  )

  ordered_car_df <- car_df
  ordered_car_df$time <- ordered(
    as.character(ordered_car_df$time),
    levels = c("0", "0.5", "2")
  )
  ordered_car_result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = ordered_car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  ordered_car_term <- ordered_car_result$formula_design$random_effects[[1]]
  expect_equal(ordered_car_term$car$time_values, c(0, 0.5, 2))
  expect_equal(
    ordered_car_result$data$mu__xREx__id_xRE_DATAx,
    car_result$data$mu__xREx__id_xRE_DATAx
  )

  car_logit <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "logit"
        )
      )
    )
  )
  expect_equal(names(car_logit$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_logit"))
  expect_match(car_logit$formula_syntax, "mu__xREx__id_rho <- 0 +", fixed = TRUE)
  expect_match(car_logit$formula_syntax, "* ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)

  ar_logit <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "logit"
        )
      )
    )
  )
  expect_equal(names(ar_logit$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_logit"))
  expect_match(ar_logit$formula_syntax, "mu__xREx__id_rho <- -1 + 2", fixed = TRUE)
  expect_match(ar_logit$formula_syntax, "* ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + x + time | id),
      parameter = "mu",
      data = transform(car_df, x = seq_len(nrow(car_df))),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "exactly one untransformed time variable",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = transform(car_df, time = factor(time)),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "must be numeric or an ordered factor",
    fixed = TRUE
  )
})

test_that("structured correlation Cholesky syntax exposes intended covariance patterns", {

  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  block_prior <- random_block(sd = sd_prior, cor = rho_prior)

  cs_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "cs",
    block_prior = block_prior,
    include_correlation = TRUE
  )
  expect_match(cs_module$syntax, "mu__xREx__id_xRE_CORx_R[1,2] <- mu__xREx__id_rho", fixed = TRUE)
  expect_match(cs_module$syntax, "mu__xREx__id_xRE_CORx_L[1,2] <- 0", fixed = TRUE)
  expect_match(cs_module$syntax, "mu__xREx__id_xRE_CORx_L[2,2] <- sqrt", fixed = TRUE)

  ar_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 4,
    structure = "ar1",
    block_prior = block_prior,
    include_correlation = TRUE
  )
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_R[1,3] <- pow(mu__xREx__id_rho, 2)", fixed = TRUE)
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_R[1,4] <- pow(mu__xREx__id_rho, 3)", fixed = TRUE)
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_L[1,4] <- 0", fixed = TRUE)
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_L[4,4] <- sqrt", fixed = TRUE)

  car_distance <- abs(outer(c(0, 1.5, 3), c(0, 1.5, 3), "-"))
  expect_error(
    BayesTools:::.bt_JAGS_structured_corr_cholesky(
      node_prefix = "mu__xREx__id",
      prior_prefix = "_xREx__id",
      K = 3,
      structure = "car",
      block_prior = block_prior,
      include_correlation = TRUE
    ),
    "Dense CAR Cholesky compilation is unsupported",
    fixed = TRUE
  )
  car_module <- BayesTools:::.bt_JAGS_structured_corr_direct(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "car",
    block_prior = block_prior,
    include_correlation = TRUE,
    distance_matrix = car_distance
  )
  expect_match(
    car_module$syntax,
    paste0(
      "mu__xREx__id_xRE_CORx_R[1,2] <- ",
      "exp(1.5 * log(mu__xREx__id_rho))"
    ),
    fixed = TRUE
  )
  expect_match(
    car_module$syntax,
    paste0(
      "mu__xREx__id_xRE_CORx_R[1,3] <- ",
      "exp(3 * log(mu__xREx__id_rho))"
    ),
    fixed = TRUE
  )
  expect_false(grepl(
    "pow(mu__xREx__id_rho",
    car_module$syntax,
    fixed = TRUE
  ))
  expect_false(grepl("_xRE_CORx_L", car_module$syntax, fixed = TRUE))
  expect_null(car_module$cholesky_name)
  expect_equal(car_module$prior_list$`_xREx__id_rho_z`$truncation$lower, 0)
  expect_equal(car_module$bridge$bounds, c(lower = 0, upper = 1))
  expect_equal(car_module$bridge$distance_matrix, car_distance)

  logit_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "cs",
    block_prior = random_block(
      sd = sd_prior,
      covariance = random_covariance(
        cor = rho_prior,
        cor_scale = "logit"
      )
    ),
    include_correlation = TRUE
  )
  expect_match(logit_module$syntax, "mu__xREx__id_rho <- -0.5 + 1.5", fixed = TRUE)
  expect_match(logit_module$syntax, "* ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)
  expect_true("_xREx__id_rho_logit" %in% names(logit_module$prior_list))
})

.jags_lme4_random_oracle_data <- function() {
  data.frame(
    y = c(1.5, 0.2, 2.1, 3.4, -0.7, 4.2),
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )
}

.lme4_lFormula_or_fallback <- function(formula, fallback_formula, data) {
  result <- try(lme4::lFormula(formula, data = data), silent = TRUE)
  if (inherits(result, "try-error")) {
    result <- lme4::lFormula(fallback_formula, data = data)
  }
  result
}

.expect_jags_random_design_matches_lme4 <- function(jags_result, lme4_result,
                                                    expected_re_columns) {
  bayes_random_data <- jags_result$data$mu__xREx__id_xRE_DATAx
  bayes_random_map <- jags_result$data$mu__xREx__id_xRE_MAPx
  lme4_re_terms <- lme4_result$reTrms
  lme4_group <- lme4_re_terms$flist$id
  lme4_Z <- t(as.matrix(lme4_re_terms$Zt))

  expect_equal(dim(jags_result$formula_design$model_matrix), dim(lme4_result$X))
  expect_equal(
    as.vector(jags_result$formula_design$model_matrix),
    as.vector(lme4_result$X)
  )
  expect_equal(colnames(jags_result$formula_design$model_matrix), colnames(lme4_result$X))

  expect_equal(bayes_random_map, unname(as.integer(lme4_group)))
  expect_equal(levels(lme4_group), levels(.jags_lme4_random_oracle_data()$id))
  expect_equal(attr(jags_result$formula_design$random_effects[[1]], "grouping_factor"), "id")
  expect_true(isTRUE(attr(jags_result$formula_design$random_effects[[1]], "independent")))
  expect_s3_class(jags_result$formula_design$random_effects[[1]], "BayesTools_random_effect_term")
  expect_equal(jags_result$formula_design$random_effects[[1]]$structure, "diag")
  expect_equal(
    jags_result$formula_design$jags_data_names[["__xREx__id"]],
    c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx")
  )

  expect_equal(dim(bayes_random_data), c(nrow(lme4_result$fr), length(expected_re_columns)))
  expect_equal(colnames(bayes_random_data), expected_re_columns)
  expect_equal(lme4_re_terms$cnms$id, expected_re_columns)
  expect_equal(dim(lme4_Z), c(nrow(lme4_result$fr), length(expected_re_columns) * nlevels(lme4_group)))

  for (column_i in seq_along(expected_re_columns)) {
    selected <- bayes_random_map + (column_i - 1L) * nlevels(lme4_group)
    expect_equal(unname(bayes_random_data[, column_i]), unname(lme4_Z[cbind(seq_len(nrow(lme4_Z)), selected)]))
  }
}

test_that("JAGS_formula independent random-effect design matches lme4 lFormula oracles", {

  skip_if_not_installed("lme4")

  df <- .jags_lme4_random_oracle_data()

  intercept_result <- JAGS_formula(
    formula = ~ x + (1 || id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  intercept_lme4 <- .lme4_lFormula_or_fallback(
    y ~ x + (1 || id),
    fallback_formula = y ~ x + (1 | id),
    data = df
  )

  .expect_jags_random_design_matches_lme4(
    jags_result = intercept_result,
    lme4_result = intercept_lme4,
    expected_re_columns = "(Intercept)"
  )

  slope_result <- JAGS_formula(
    formula = ~ x + (0 + x || id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  slope_lme4 <- lme4::lFormula(y ~ x + (0 + x || id), data = df)

  .expect_jags_random_design_matches_lme4(
    jags_result = slope_result,
    lme4_result = slope_lme4,
    expected_re_columns = "x"
  )

  df_us <- data.frame(
    y = seq_len(12) / 10,
    x = c(-1, 0, 1, 2, -2, 3, -0.5, 0.5, 1.5, -1.5, 2.5, -2.5),
    id = factor(rep(c("c", "a", "b"), 4), levels = c("c", "a", "b"))
  )
  us_result <- JAGS_formula(
    formula = ~ x + (1 + x | id),
    parameter = "mu",
    data = df_us,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )
  us_lme4 <- lme4::lFormula(y ~ x + (1 + x | id), data = df_us)
  us_map <- unname(as.integer(us_lme4$reTrms$flist$id))

  expect_equal(us_result$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(
    as.vector(us_result$formula_design$model_matrix),
    as.vector(us_lme4$X)
  )
  expect_equal(us_result$formula_design$random_effects[[1]]$column_names, us_lme4$reTrms$cnms$id)
  expect_equal(us_result$data$mu__xREx__id_xRE_MAPx, us_map)
  expect_equal(colnames(us_result$data$mu__xREx__id_xRE_DATAx), c("(Intercept)", "x"))
  expect_equal(unname(us_result$data$mu__xREx__id_xRE_DATAx[, "(Intercept)"]), rep(1, nrow(df_us)))
  expect_equal(unname(us_result$data$mu__xREx__id_xRE_DATAx[, "x"]), df_us$x)
  expect_match(us_result$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)

  df_ar1 <- data.frame(
    y = seq_len(12) / 10,
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )
  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = df_ar1,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  ar1_lme4 <- lme4::lFormula(y ~ 1 + ar1(0 + f | id), data = df_ar1)

  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names, ar1_lme4$reTrms$cnms$id)
  expect_equal(ar1_result$data$mu__xREx__id_xRE_MAPx, unname(as.integer(ar1_lme4$reTrms$flist$id)))
  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(ar1_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(ar1_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_true(ar1_result$formula_design$random_effects[[1]]$homogeneous_sd)
})

test_that("random slopes of a two-level factor name the scalar SD node of their prior", {

  testthat::skip_if_not_installed("runjags")

  # one contrast column is a scalar node of the SD prior, not an indexed element
  df <- data.frame(
    x = c(-1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2),
    f2 = factor(rep(c("a", "b"), 6), levels = c("a", "b")),
    f3 = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    id = factor(rep(paste0("s", 1:4), each = 3))
  )
  orthonormal_df <- df
  stats::contrasts(orthonormal_df$f2) <- "contr.orthonormal"
  priors <- list(intercept = prior("normal", list(0, 1)))
  sd_prior <- prior("gamma", list(2, 2))
  monitor <- random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)

  sd_prior_nodes <- function(result){
    term <- result$formula_design$random_effects[[1L]]
    owners <- unique(sub("\\[[0-9]+\\]$", "", term$sd_parameter_names))
    unlist(lapply(owners, function(owner){
      owner_prior <- result$prior_list[[owner]]
      if(is.prior.factor(owner_prior)){
        BayesTools:::.JAGS_prior_factor_names(owner, owner_prior)
      }else{
        owner
      }
    }), use.names = FALSE)
  }
  sd_prior_random <- function(correlated){
    prior_random(id = random_block(
      sd = sd_prior,
      cor = if(correlated) prior_lkj(eta = 1) else NULL,
      monitor = monitor
    ))
  }
  # the carrier monitors the nodes JAGS emits for the SD priors
  carrier_fit <- function(result, correlated){
    term <- result$formula_design$random_effects[[1L]]
    sd_names <- sd_prior_nodes(result)
    values <- c(mu_intercept = 0, stats::setNames(seq_along(sd_names) / 2, sd_names))
    if(correlated){
      n <- length(sd_names)
      cor_matrix <- diag(n)
      cor_matrix[1, 2] <- cor_matrix[2, 1] <- 0.4
      R_names <- outer(seq_len(n), seq_len(n), function(row, column){
        paste0(term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
      })
      L_names <- BayesTools:::.bt_random_effect_cholesky_names(term, n)
      values <- c(
        values,
        stats::setNames(as.vector(cor_matrix), as.vector(R_names)),
        stats::setNames(as.vector(t(chol(cor_matrix))), as.vector(L_names))
      )
    }
    posterior <- matrix(rep(values, 2L), nrow = 2L, byrow = TRUE,
                        dimnames = list(NULL, names(values)))
    fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)),
                summary.pars = list(mutate = NULL),
                monitor = colnames(posterior), sample = nrow(posterior))
    class(fit) <- c("runjags", "BayesTools_fit")
    attr(fit, "prior_list") <- result$prior_list
    attr(fit, "formula_design") <- list(mu = result$formula_design)
    attr(fit, "formula_scale") <- list(mu = result$formula_scale)
    attach_test_parameter_map(fit)
  }

  cases <- list(
    list(formula = ~ 1 + diag(1 + f2 | id), data = df, correlated = FALSE,
         names = c("mu__xREx__id_intercept", "mu__xREx__id_f2"),
         rows = c("(mu) sd(intercept)", "(mu) sd(f2[b])")),
    list(formula = ~ 1 + diag(0 + f2 | id), data = df, correlated = FALSE,
         names = "mu__xREx__id_f2", rows = "(mu) sd(f2[b])"),
    list(formula = ~ 1 + us(1 + f2 | id), data = df, correlated = TRUE,
         names = c("mu__xREx__id_intercept", "mu__xREx__id_f2"),
         rows = c("(mu) sd(intercept)", "(mu) sd(f2[b])", "(mu) cor(intercept,f2[b])")),
    list(formula = ~ 1 + diag(1 + f2 * x | id, hom = FALSE), data = df, correlated = FALSE,
         names = c("mu__xREx__id_intercept", "mu__xREx__id_f2", "mu__xREx__id_x",
                   "mu__xREx__id_f2__xXx__x"),
         rows = c("(mu) sd(intercept)", "(mu) sd(f2[b])", "(mu) sd(x)", "(mu) sd(f2[b]:x)")),
    list(formula = ~ 1 + diag(1 + f2 | id), data = orthonormal_df, correlated = FALSE,
         names = c("mu__xREx__id_intercept", "mu__xREx__id_f2"),
         # an orthonormal factor term is one unit, labelled by the term
         rows = c("(mu) sd(intercept)", "(mu) sd(f2)"))
  )
  for(case in cases){
    label <- paste(deparse(case$formula), "(", case$names[length(case$names)], ")")
    result <- JAGS_formula(case$formula, "mu", case$data, priors,
                           prior_random = sd_prior_random(case$correlated))
    term <- result$formula_design$random_effects[[1L]]
    expect_identical(unique(term$sd_parameter_names), case$names, label = label)
    expect_identical(unname(term$sd_leaves$leaf_names_by_column), unname(term$sd_parameter_names),
                     label = label)
    expect_identical(unique(term$sd_parameter_names), sd_prior_nodes(result), label = label)
    syntax_lines <- trimws(unlist(strsplit(result$formula_syntax, "\n", fixed = TRUE)))
    for(column in seq_along(case$names)){
      expect_true(
        paste0(term$parameter_stem, "_xRE_STDx[", column, "] = ", case$names[column]) %in% syntax_lines,
        label = paste(label, "column", column)
      )
    }
    expect_false(grepl("f2[1]", result$formula_syntax, fixed = TRUE), label = label)

    fit <- carrier_fit(result, case$correlated)
    table <- JAGS_estimates_table(fit, remove_diagnostics = TRUE)
    expect_true(all(case$rows %in% rownames(table)), label = label)
  }

  # factors with three levels keep their indexed nodes
  three <- JAGS_formula(~ 1 + diag(1 + f3 | id), "mu", df, priors,
                        prior_random = sd_prior_random(FALSE))
  expect_identical(
    three$formula_design$random_effects[[1L]]$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_f3[1]", "mu__xREx__id_f3[2]")
  )
  expect_identical(
    unique(three$formula_design$random_effects[[1L]]$sd_parameter_names),
    c("mu__xREx__id_intercept", sd_prior_nodes(three)[-1L])
  )
})

test_that("structured random-effect indexes with a single level compile", {

  testthat::skip_if_not_installed("runjags")

  # one index level is one indicator column and has no correlation parameter
  one_level <- data.frame(g = factor(c("a", "a", "b", "b")), time = 0)
  priors <- list(intercept = prior("normal", list(0, 1)))
  sd_prior <- prior("normal", list(0, 1), truncation = list(0, Inf))
  compile <- function(formula, data){
    JAGS_formula(formula, "mu", data, priors, prior_random = prior_random(
      sd = sd_prior, cor = prior("normal", list(0, 0.5))
    ))
  }

  for(structure in c("cs", "hcs", "ar1", "har")){
    result <- compile(stats::as.formula(paste0("~ 1 + ", structure, "(time | g)")), one_level)
    term <- result$formula_design$random_effects[[1L]]
    expect_identical(term$structure, structure, label = structure)
    expect_identical(term$n_columns, 1L, label = structure)
    expect_identical(term$column_names, "time0", label = structure)
    sd_name <- if(structure %in% c("cs", "ar1")) "mu__xREx__g_sd" else "mu__xREx__g_time"
    expect_identical(unique(term$sd_parameter_names), sd_name, label = structure)
    expect_identical(names(result$prior_list), c("mu_intercept", sd_name), label = structure)
    expect_false(grepl("rho|_CORx_|dbt_lkj", result$formula_syntax), label = structure)

    posterior <- matrix(c(0, 0.5), nrow = 2L, ncol = 2L, byrow = TRUE,
                        dimnames = list(NULL, c("mu_intercept", sd_name)))
    fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)),
                summary.pars = list(mutate = NULL),
                monitor = colnames(posterior), sample = nrow(posterior))
    class(fit) <- c("runjags", "BayesTools_fit")
    attr(fit, "prior_list") <- result$prior_list
    attr(fit, "formula_design") <- list(mu = result$formula_design)
    attr(fit, "formula_scale") <- list(mu = result$formula_scale)
    table <- JAGS_estimates_table(attach_test_parameter_map(fit), remove_diagnostics = TRUE)
    expect_identical(
      rownames(table),
      c("(mu) intercept", if(structure %in% c("cs", "ar1")) "(mu) sd" else "(mu) sd(time[0])"),
      label = structure
    )
  }

  # replayed data with a single index level builds the same indicator column
  replay <- BayesTools:::.bt_random_effect_design_matrix(
    ~ time - 1,
    data.frame(time = factor(rep("0", 3L))),
    preserve_no_intercept_contrasts = FALSE,
    structure = "hcs"
  )$model_matrix
  expect_identical(dim(replay), c(3L, 1L))
  expect_true(all(replay == 1))
  expect_identical(colnames(replay), "time0")

  # the indicator matrix equals the model.matrix() coding of a multi-level index
  index_data <- data.frame(
    time = factor(c("b", "a", "c", "a", "b"), levels = c("a", "b", "c"))
  )
  stats::contrasts(index_data$time) <- "contr.treatment"
  for(structure in c("cs", "hcs", "ar1", "har")){
    indicator <- BayesTools:::.bt_random_effect_design_matrix(
      ~ time - 1, index_data, preserve_no_intercept_contrasts = FALSE, structure = structure
    )$model_matrix
    reference <- BayesTools:::.bt_model_matrix(
      stats::model.frame(~ time - 1, index_data), formula = ~ time - 1, data = index_data
    )
    expect_identical(indicator, reference, label = structure)
  }
  multi_level <- compile(~ 1 + hcs(time | g), data.frame(
    g = factor(rep(c("a", "b"), each = 2L)), time = rep(1:2, 2L)
  ))
  expect_identical(multi_level$formula_design$random_effects[[1L]]$n_columns, 2L)
  expect_identical(
    multi_level$formula_design$random_effects[[1L]]$sd_parameter_names,
    paste0("mu__xREx__g_time[", 1:2, "]")
  )
})

test_that("generated grouping labels are unique against the final label set", {

  unique_labels <- BayesTools:::.bt_random_group_unique_labels

  # shared labels take their tuple key; unique labels are unchanged
  expect_identical(
    unique_labels(c("a", "b", "a"), c("k1", "k2", "k3")),
    c("a [k1]", "b", "a [k3]")
  )
  # a generated label that equals an ordinary label takes the smallest unused
  # integer suffix and leaves the ordinary label alone
  expect_identical(
    unique_labels(c("a", "a", "a [k1]"), c("k1", "k2", "k3")),
    c("a [k1] #2", "a [k2]", "a [k1]")
  )
  expect_identical(
    unique_labels(c("a", "a", "a [k1]", "a [k1] #2"), c("k1", "k2", "k3", "k4")),
    c("a [k1] #3", "a [k2]", "a [k1]", "a [k1] #2")
  )
  # two generated labels that coincide keep the first
  expect_identical(
    unique_labels(c("a [x]", "a [x]", "a", "a"), c("k", "k", "x", "y")),
    c("a [x] [k]", "a [x] [k] #2", "a [x]", "a [y]")
  )
  # and the same against the labels that exist already (new prediction groups)
  expect_identical(
    unique_labels("a", "k1", existing = c("a", "a [k1]")),
    "a [k1] #2"
  )
  expect_false(anyDuplicated(unique_labels(
    c("a", "a", "a [k1]", "a [k1] #2", "b"), c("k1", "k2", "k3", "k4", "k5"),
    existing = c("a [k2] #2")
  )) > 0L)

  # interaction groups whose labels contain the separator and a generated suffix
  collision_data <- data.frame(
    a = c("a", "a:b", "a"),
    b = c("b:c", "c", "b:c [3:c:a|5:c:b:c]")
  )
  built <- JAGS_formula(
    ~ 1 + diag(1 | a:b),
    "mu",
    collision_data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      a_b = random_block(sd = prior("normal", list(0, 1), truncation = list(0, Inf)))
    )
  )
  term <- built$formula_design$random_effects[[1L]]
  expect_identical(term$n_groups, 3L)
  expect_identical(
    term$group_levels,
    c("a:b:c [3:c:a|5:c:b:c] #2", "a:b:c [3:c:a|5:c:b:c]", "a:b:c [5:c:a:b|3:c:c]")
  )
  expect_identical(term$group_map, c(1L, 3L, 2L))

  # prediction replay of the fitted data matches every row to its fitted group
  replay <- BayesTools:::.bt_random_effect_prediction_data(
    term, collision_data, group_data = collision_data
  )
  expect_identical(replay$group_map, term$group_map)
  expect_identical(replay$group_levels, term$group_levels)
})

test_that("car() refuses an explicit intercept and accepts the redundant 0 and - 1", {

  car_data <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  compile <- function(formula){
    JAGS_formula(
      formula, "mu", car_data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("normal", list(0, 0.5))
      ))
    )
  }

  # 'car(1 + time | id)' used to drop the requested intercept without a message
  for(term in c("1 + time", "time + 1", "(1) + time", "time - 0", "0 + 1 + time")){
    for(formula in list(
      stats::as.formula(paste0("~ 1 + car(", term, " | id)")),
      stats::as.formula(paste0("~ 1 + random(", term, " | id, covariance = \"car\")"))
    )){
      expect_error(
        compile(formula),
        "CAR random-effect terms have no intercept, so an explicit '1' is not accepted",
        fixed = TRUE,
        info = deparse(formula)
      )
    }
  }
  expect_error(
    compile(~ 1 + car(1 + time | id)),
    "'0' and '- 1' are accepted and redundant.",
    fixed = TRUE
  )

  # the time variable is checked first
  expect_error(
    compile(~ 1 + car(1 | id)),
    "CAR random-effect terms require exactly one untransformed time variable",
    fixed = TRUE
  )
  expect_error(
    compile(~ 1 + car(1 + time + id | id)),
    "CAR random-effect terms require exactly one untransformed time variable",
    fixed = TRUE
  )

  # 0 and - 1 request what CAR already does
  reference <- compile(~ 1 + car(time | id))
  for(term in c("0 + time", "time - 1", "-1 + time", "time + 0", "0 + time - 1")){
    accepted <- compile(stats::as.formula(paste0("~ 1 + car(", term, " | id)")))
    expect_identical(accepted$formula_syntax, reference$formula_syntax, info = term)
    expect_identical(accepted$data, reference$data, info = term)
    expect_identical(
      accepted$formula_design$random_effects[[1L]]$column_names,
      reference$formula_design$random_effects[[1L]]$column_names,
      info = term
    )
  }
})

# The estimates table of a block with independent SDs, from a synthetic carrier
# fit that holds the SD nodes of the block's priors.
.oracle_independent_sd_table <- function(result){

  term <- result$formula_design$random_effects[[1L]]
  owners <- unique(sub("\\[[0-9]+\\]$", "", term$sd_parameter_names))
  nodes <- unlist(lapply(owners, function(owner){
    owner_prior <- result$prior_list[[owner]]
    if(is.prior.factor(owner_prior)){
      BayesTools:::.JAGS_prior_factor_names(owner, owner_prior)
    }else{
      owner
    }
  }), use.names = FALSE)
  values <- c(mu_intercept = 0, stats::setNames(seq_along(nodes) / 2, nodes))
  posterior <- matrix(rep(values, 2L), nrow = 2L, byrow = TRUE,
                      dimnames = list(NULL, names(values)))
  fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)),
              summary.pars = list(mutate = NULL),
              monitor = colnames(posterior), sample = nrow(posterior))
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "formula_scale") <- list(mu = result$formula_scale)

  JAGS_estimates_table(attach_test_parameter_map(fit), remove_diagnostics = TRUE)
}

test_that("random factor interactions match factors to terms by their exact names", {

  testthat::skip_if_not_installed("runjags")

  # factors named x_g and g (and g and age_g) coexist in one term
  df <- data.frame(
    x_g = factor(rep(c("a", "b"), 6)),
    g = factor(rep(c("p", "q", "r"), 4)),
    age_g = factor(rep(c("u", "v"), each = 6)),
    h = factor(rep(c("x", "y", "z"), 4)),
    id = factor(rep(paste0("s", 1:4), each = 3))
  )
  priors <- list(intercept = prior("normal", list(0, 1)))
  compile <- function(formula){
    JAGS_formula(formula, "mu", df, priors,
                 prior_random = prior_random(sd = prior("gamma", list(2, 2))))
  }
  leaf_terms <- function(result){
    result$formula_design$random_effects[[1L]]$sd_leaves$leaf_terms
  }

  suffix_first <- compile(~ 1 + diag(1 + x_g:g | id))
  expect_identical(
    unname(leaf_terms(suffix_first)),
    c("intercept", paste0("x_g[", c("a", "b"), "]__xXx__g[", rep(c("p", "q", "r"), each = 2L), "]"))
  )
  suffix_last <- compile(~ 1 + diag(1 + g:x_g | id))
  expect_identical(
    unname(leaf_terms(suffix_last)),
    c("intercept", paste0("g[", c("p", "q", "r"), "]__xXx__x_g[", rep(c("a", "b"), each = 3L), "]"))
  )
  expect_identical(
    unname(leaf_terms(compile(~ 1 + diag(1 + g * age_g | id, hom = FALSE)))),
    c("intercept", "g[q]", "g[r]", "age_g[v]", "g[q]__xXx__age_g[v]", "g[r]__xXx__age_g[v]")
  )
  expect_identical(
    unname(leaf_terms(compile(~ 1 + diag(1 + age_g * g | id, hom = FALSE)))),
    c("intercept", "age_g[v]", "g[q]", "g[r]", "age_g[v]__xXx__g[q]", "age_g[v]__xXx__g[r]")
  )
  expect_identical(
    unname(leaf_terms(compile(~ 1 + diag(1 + h:g | id))))[2:4],
    c("h[x]__xXx__g[p]", "h[y]__xXx__g[p]", "h[z]__xXx__g[p]")
  )

  # the coordinate labels of the same priors are formatted by the same rule
  expect_identical(
    rownames(.oracle_independent_sd_table(suffix_first)),
    c("(mu) intercept", "(mu) sd(intercept)",
      paste0("(mu) sd(x_g[", c("a", "b"), "]:g[", rep(c("p", "q", "r"), each = 2L), "])"))
  )
  expect_identical(
    rownames(.oracle_independent_sd_table(suffix_last)),
    c("(mu) intercept", "(mu) sd(intercept)",
      paste0("(mu) sd(g[", c("p", "q", "r"), "]:x_g[", rep(c("a", "b"), each = 3L), "])"))
  )

  # the formatter itself: exact names, the random-term prefix, no suffix matching
  format_names <- BayesTools:::.format_factor_level_parameter_names
  levels <- list(x_g = c("a", "b"), g = c("p", "q"))
  expect_identical(
    format_names("x_g__xXx__g", levels, 4L),
    c("x_g[a]__xXx__g[p]", "x_g[b]__xXx__g[p]", "x_g[a]__xXx__g[q]", "x_g[b]__xXx__g[q]")
  )
  expect_identical(
    format_names("mu__xREx__id_x_g__xXx__g", levels, 4L, prefix = "mu__xREx__id_"),
    c("mu__xREx__id_x_g[a]__xXx__g[p]", "mu__xREx__id_x_g[b]__xXx__g[p]",
      "mu__xREx__id_x_g[a]__xXx__g[q]", "mu__xREx__id_x_g[b]__xXx__g[q]")
  )
  expect_error(
    format_names("mu__xREx__id_x_g__xXx__g", levels, 4L),
    "factor metadata do not match",
    fixed = TRUE
  )
  expect_error(
    format_names("pre_x__xXx__g", list(x = c("a", "b"), g = c("p", "q")), 4L),
    "factor metadata do not match",
    fixed = TRUE
  )
})

test_that("SD components of mixed-contrast random interactions are level cells only when the column is that cell", {

  testthat::skip_if_not_installed("runjags")

  df <- data.frame(
    study = factor(rep(paste0("s", 1:6), each = 12)),
    f = factor(rep(c("x", "y", "z"), 24), levels = c("x", "y", "z")),
    g = factor(rep(c("p", "q"), each = 36), levels = c("p", "q"))
  )
  compile <- function(formula, f_contrast, g_contrast){
    data <- df
    stats::contrasts(data$f) <- f_contrast
    stats::contrasts(data$g) <- g_contrast
    JAGS_formula(
      formula, "mu", data, list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = prior("gamma", list(2, 2)))
    )
  }
  leaf_terms <- function(result){
    unname(result$formula_design$random_effects[[1L]]$sd_leaves$leaf_terms)
  }
  both <- ~ 1 + diag(1 + f * g | study, hom = FALSE)

  # an orthonormal or mean-difference factor is an exchangeable unit with one
  # SD; crossed with a treatment factor it gives one unit per partner column,
  # labelled by the partner's level
  treatment_orthonormal <- compile(both, "contr.treatment", "contr.orthonormal")
  expect_identical(
    leaf_terms(treatment_orthonormal),
    c("intercept", "f[y]", "f[z]", "g", "f[y]__xXx__g", "f[z]__xXx__g")
  )
  expect_identical(
    leaf_terms(compile(both, "contr.orthonormal", "contr.treatment")),
    c("intercept", "f", "g[q]", "f__xXx__g[q]")
  )
  expect_identical(
    leaf_terms(compile(both, "contr.meandif", "contr.treatment")),
    c("intercept", "f", "g[q]", "f__xXx__g[q]")
  )
  # the public labels follow the leaf terms
  expect_identical(
    rownames(.oracle_independent_sd_table(treatment_orthonormal)),
    c("(mu) intercept", "(mu) sd(intercept)", "(mu) sd(f[y])", "(mu) sd(f[z])",
      "(mu) sd(g)", "(mu) sd(f[y]:g)", "(mu) sd(f[z]:g)")
  )

  # treatment and treatment: the interaction columns are the (y, q) and (z, q) cells
  treatment_treatment <- compile(both, "contr.treatment", "contr.treatment")
  expect_identical(
    leaf_terms(treatment_treatment),
    c("intercept", "f[y]", "f[z]", "g[q]", "f[y]__xXx__g[q]", "f[z]__xXx__g[q]")
  )
  # independent f and treatment g: the interaction columns are direct cells too
  independent_treatment <- compile(both, "contr.independent", "contr.treatment")
  expect_identical(
    leaf_terms(independent_treatment),
    c("intercept", "f[x]", "f[y]", "f[z]", "g[q]",
      "f[x]__xXx__g[q]", "f[y]__xXx__g[q]", "f[z]__xXx__g[q]")
  )
  expect_identical(
    rownames(.oracle_independent_sd_table(independent_treatment)),
    c("(mu) intercept", "(mu) sd(intercept)", "(mu) sd(f[x])", "(mu) sd(f[y])",
      "(mu) sd(f[z])", "(mu) sd(g[q])", "(mu) sd(f[x]:g[q])", "(mu) sd(f[y]:g[q])",
      "(mu) sd(f[z]:g[q])")
  )
  # the f:g-only term is coded by level indicators
  only_interaction <- compile(~ 1 + diag(1 + f:g | study, hom = FALSE),
                              "contr.treatment", "contr.treatment")
  expect_identical(
    leaf_terms(only_interaction),
    c("intercept", paste0("f[", c("x", "y", "z"), "]__xXx__g[", rep(c("p", "q"), each = 3L), "]"))
  )

  # R codes the f:g-only term by level indicators whatever the declared
  # contrasts, so its columns are cells also for a mixed-contrast term; the
  # tag gate stays for terms whose margins are present
  cells <- c("intercept", paste0("f[", c("x", "y", "z"), "]__xXx__g[", rep(c("p", "q"), each = 3L), "]"))
  interaction_only <- ~ 1 + diag(1 + f:g | study, hom = FALSE)
  for(contrasts in list(c("contr.treatment", "contr.orthonormal"),
                        c("contr.orthonormal", "contr.treatment"),
                        c("contr.meandif", "contr.treatment"))){
    mixed_only <- compile(interaction_only, contrasts[1L], contrasts[2L])
    expect_identical(leaf_terms(mixed_only), cells, info = paste(contrasts, collapse = " x "))
    expect_identical(
      rownames(.oracle_independent_sd_table(mixed_only)),
      c("(mu) intercept", "(mu) sd(intercept)",
        paste0("(mu) sd(f[", c("x", "y", "z"), "]:g[", rep(c("p", "q"), each = 3L), "])")),
      info = paste(contrasts, collapse = " x ")
    )
  }
  # the factor of the interaction that is coded by its contrast is a unit per
  # level of the indicator-coded factor: orthonormal g with its margin f in
  # the block, orthonormal f with g
  expect_identical(
    leaf_terms(compile(~ 1 + diag(1 + f + f:g | study, hom = FALSE), "contr.treatment", "contr.orthonormal")),
    c("intercept", "f[y]", "f[z]", paste0("f[", c("x", "y", "z"), "]__xXx__g"))
  )
  expect_identical(
    leaf_terms(compile(~ 1 + diag(1 + g + f:g | study, hom = FALSE), "contr.orthonormal", "contr.treatment")),
    c("intercept", "g[q]", paste0("g[", c("p", "q"), "]__xXx__f"))
  )
  # while f, coded by indicators next to a treatment-coded g, gives cells
  expect_identical(
    leaf_terms(compile(~ 1 + diag(1 + f + f:g | study, hom = FALSE), "contr.orthonormal", "contr.treatment")),
    c("intercept", "f", paste0("f[", c("x", "y", "z"), "]__xXx__g[q]"))
  )
  # an orthonormal main-effect design is one unit
  expect_identical(
    leaf_terms(compile(~ 1 + diag(1 + g | study, hom = FALSE), "contr.treatment", "contr.orthonormal")),
    c("intercept", "g")
  )
})
