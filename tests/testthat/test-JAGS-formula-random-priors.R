skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Random-effect priors, SD sources and variance allocation
# ============================================================================ #
#
# PURPOSE:
#   prior_random() settings, external SD sources, and variance allocation priors.
# ============================================================================ #

test_that("prior_random rejects unsupported and ignored production settings", {

  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  lkj_prior <- prior_lkj(eta = 2)
  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    idx = factor(rep(c("t1", "t2"), 3), levels = c("t1", "t2")),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  group_as_data_frame <- function(x){
    data.frame(group = x)
  }

  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = transform(df, id = factor(c("b", NA, "b", "c", "a", "c"))),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "must not contain missing values",
    fixed = TRUE
  )
  explicit_na_id <- factor(
    c("b", NA, "b", "c", "a", "c"),
    levels = c("a", "b", "c", NA),
    exclude = NULL
  )
  expect_false(anyNA(explicit_na_id))
  expect_true(anyNA(levels(explicit_na_id)))
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = transform(df, id = explicit_na_id),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "must not contain missing values",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | group_as_data_frame(id), name = "bad_group", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(bad_group = random_block(sd = sd_prior))
    ),
    "grouping expressions must be variables",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | id[1], name = "bad_group", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(bad_group = random_block(sd = sd_prior))
    ),
    "grouping expressions must be variables",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | unknown_group_function(id), name = "bad_group", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(bad_group = random_block(sd = sd_prior))
    ),
    "grouping expressions must be variables",
    fixed = TRUE
  )
  expect_error(
    prior_random(random_block(sd = sd_prior)),
    "must be named",
    fixed = TRUE
  )
  expect_error(
    prior_random(
      study = random_block(sd = sd_prior),
      study = random_block(sd = sd_prior)
    ),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    prior_random(study = sd_prior),
    "must be created with random_block",
    fixed = TRUE
  )
  expect_s3_class(
    random_block(terms = list(intercept = sd_prior, x = random_block(sd = sd_prior))),
    "random_block"
  )
  expect_error(
    random_block(terms = list(sd_prior)),
    "must be named",
    fixed = TRUE
  )
  expect_error(
    random_block(terms = list(intercept = sd_prior, intercept = sd_prior)),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    random_block(terms = list(intercept = 1)),
    "must be a prior or random_block",
    fixed = TRUE
  )
  malformed_terms_block <- random_block(sd = sd_prior)
  malformed_terms_block$terms <- list(intercept = 1)
  expect_error(
    prior_random(id = malformed_terms_block),
    "must be a prior or random_block",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = prior("point", list(-1))))
    ),
    # random_block() rejects the negative SD support before the formula is built
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + missing_x | id),
      parameter = "mu",
      data = df["id"],
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "predictor variable is missing",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + x | id),
      parameter = "mu",
      data = transform(df, x = c(-1, NA, 1, 2, -2, 3)),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "missing predictor values",
    fixed = TRUE
  )

  allocation_prior <- random_variance_allocation(name = "allocation",
    terms = c("study", "drug"),
    sd = sd_prior,
    weights = prior("dirichlet", list(alpha = c(2, 3)))
  )
  expect_s3_class(allocation_prior, "random_variance_allocation")
  expect_s3_class(prior_random(allocation = allocation_prior), "prior_random")
  named_dot_allocation <- prior_random(total = allocation_prior)
  expect_equal(
    names(BayesTools:::.bt_random_allocation_list(named_dot_allocation$allocation)),
    "total"
  )
  expect_error(
    prior_random(`bad-name` = allocation_prior),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    prior_random(allocation = stats::setNames(list(allocation_prior), "bad-name")),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    prior_random(`xRE_ALLOCx` = allocation_prior),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    prior_random(allocation = stats::setNames(list(allocation_prior), "__xRE_SUMMARY__")),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "study"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = "study",
      sd = sd_prior
    ),
    "at least two resolved random-effect blocks",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug"),
      sd = sd_prior,
      weights = prior("beta", list(1, 1))
    ),
    "Dirichlet simplex priors",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = prior("spike", list(location = 0.5))
    ),
    "named list of scalar probability priors",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = list(study = 0.5)
    ),
    "must be a prior object",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = list(`bad-name` = prior("spike", list(location = 0.5)))
    ),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = list(study = prior("spike", list(location = 1.5)))
    ),
    "bounded within 0 and 1",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1, 1)))
    ),
    "dimension must match",
    fixed = TRUE
  )
  hidden_component_allocation <- random_variance_allocation(
    name            = "hidden_component",
    terms           = c(component_1 = "study"),
    sd              = sd_prior,
    inclusion       = list(
      component_1 = prior("spike", list(location = 0.5))
    ),
    component_names = ""
  )
  expect_identical(hidden_component_allocation$component_names, "")
  expect_identical(
    BayesTools:::.bt_random_effect_allocation_component_name(
      hidden_component_allocation,
      "component_1"
    ),
    ""
  )
  expect_error(
    random_variance_allocation(
      name            = "bad_hidden_component",
      terms           = c(study = "study", drug = "drug"),
      sd              = sd_prior,
      component_names = c("", "drug")
    ),
    "sole empty name",
    fixed = TRUE
  )
  positional_allocation <- random_variance_allocation(
    "allocation",
    c("study", "drug"),
    sd_prior,
    prior("dirichlet", list(alpha = c(1, 1)))
  )

  unnamed_single <- random_effects_formula(list(~ 1 | study))
  expect_equal(
    vapply(unnamed_single$terms, `[[`, character(1), "block_name"),
    "study"
  )
  expect_equal(unnamed_single$components, list(component_1 = "study"))
  expect_equal(unnamed_single$terms[[1L]]$component, "component 1")
  expect_false(unnamed_single$terms[[1L]]$component_visible)

  bare_single <- random_effects_formula(~ 1 | study)
  expect_equal(bare_single$components, list(component_1 = "study"))
  expect_equal(bare_single$terms[[1L]]$component, "component 1")
  expect_false(bare_single$terms[[1L]]$component_visible)

  unnamed_multiple <- random_effects_formula(
    list(~ 1 | study, ~ 1 | id)
  )
  expect_equal(
    vapply(unnamed_multiple$terms, `[[`, character(1), "block_name"),
    c("component_1", "component_2")
  )
  expect_equal(
    vapply(unnamed_multiple$terms, `[[`, character(1), "component"),
    c("component 1", "component 2")
  )
  expect_true(all(vapply(
    unnamed_multiple$terms,
    `[[`,
    logical(1),
    "component_visible"
  )))
  expect_s3_class(positional_allocation$weights, "prior.simplex")
  positional_named_allocation <- random_variance_allocation(
    "total_re",
    c("study", "drug"),
    sd_prior,
    prior("dirichlet", list(alpha = c(1, 1)))
  )
  expect_equal(positional_named_allocation$name, "total_re")
  positional_child_allocation <- random_variance_allocation(
    "nested_split",
    c("paper", "estimate"),
    NULL,
    prior("dirichlet", list(alpha = c(1, 1))),
    allocation_ref("total_re", "study")
  )
  expect_equal(positional_child_allocation$name, "nested_split")
  expect_s3_class(positional_child_allocation$parent, "random_allocation_ref")
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      name = "bad-name"
    ),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      name = "xRE_ALLOCx"
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      name = "__xRE_SUMMARY__"
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(xRE_ALLOCx = "study", drug = "drug"),
      sd = sd_prior
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    allocation_ref("xRE_ALLOCx", "study"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    allocation_ref("total_re", "__xRE_SUMMARY__"),
    "internally used",
    fixed = TRUE
  )
  corrupted_ref <- allocation_ref("total_re", "study")
  corrupted_ref$component <- "__xRE_SUMMARY__"
  expect_error(
    BayesTools:::.bt_check_random_allocation_ref(corrupted_ref),
    "internally used",
    fixed = TRUE
  )
  corrupted_ref$component <- "bad-name"
  expect_error(
    BayesTools:::.bt_check_random_allocation_ref(corrupted_ref),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(total = sd_prior)),
    "Variance allocation 'total' is invalid",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(allocation_prior, sd_prior)),
    "Variance allocation #2 is invalid",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(total = sd_prior)),
    "entries must be created with random_variance_allocation",
    fixed = TRUE
  )
  expect_false("allocation" %in% names(formals(random_block)))
  expect_false("allocation" %in% names(random_block()))
  expect_false("allocation" %in% names(
    BayesTools:::.bt_random_prior_for_block(prior_random(), "study")
  ))
  expect_s3_class(
    prior_random(sd = sd_prior, new_levels = random_new_levels(method = "zero")),
    "prior_random"
  )
  expect_error(
    random_covariance(cor = prior_factor("normal", list(0, 1), contrast = "treatment")),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(cor = prior("mnormal", list(mean = 0, sd = 1, K = 2))),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(cor = prior_none(), cor_scale = "logit"),
    "cannot use prior_none",
    fixed = TRUE
  )
  expect_error(
    random_covariance(
      cor = prior_mixture(list(
        prior_factor("normal", list(0, 1), contrast = "treatment"),
        prior("point", list(0))
      ))
    ),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(structure = "diag", eta = 2),
    "structure 'diag' has no correlation parameter",
    fixed = TRUE
  )
  expect_error(
    random_covariance(eta = 2, cor = lkj_prior),
    "explicit 'cor'",
    fixed = TRUE
  )
  expect_error(
    random_covariance(eta = 2, cor = rho_prior),
    "explicit 'cor'",
    fixed = TRUE
  )
  expect_s3_class(random_covariance(structure = "car"), "random_covariance")
  expect_error(
    random_covariance(structure = "car", cor = prior_lkj(eta = 1)),
    "uses a scalar correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + car(1 | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "CAR random-effect terms require exactly one untransformed time variable"
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + (1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "diag")
        )
      )
    ),
    "The formula owns the covariance structure",
    fixed = TRUE
  )
  default_us <- JAGS_formula(
    formula = ~ 1 + x + (1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(default_us$formula_design$random_effects[[1]]$correlation$eta, 1)
  expect_true("mu__xREx__id_xRE_CORx_lkj_u[1]" %in% default_us$add_parameters)
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(sd = sd_prior)
        )
      )
    ),
    "SD prior was supplied both",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + x | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = prior_factor("gamma", list(2, 2), contrast = "treatment"))
      )
    ),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = prior_mixture(list(prior("point", list(-1)), prior("gamma", list(2, 2)))))
      )
    ),
    # random_block() checks each mixture component's support at construction
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  zero_sd <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior("point", list(0))))
  )
  expect_s3_class(zero_sd$prior_list$mu__xREx__id_intercept, "prior.point")

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior, cor = lkj_prior))
    ),
    "structure 'diag' has no correlation parameter",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + (1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior, cor = rho_prior))
    ),
    "structure 'us' uses an LKJ correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + ar1(idx | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior, cor = lkj_prior))
    ),
    "structure 'ar1' uses a scalar correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + id(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        covariance = random_covariance(structure = "id"),
        id = random_block(sd = sd_prior, cor = rho_prior)
      )
    ),
    "random-effect block 'id' supplies a correlation prior, but structure 'id' has no correlation parameter.",
    fixed = TRUE
  )
})

test_that("parameter and random SD sources validate simple external references", {

  sd_prior <- prior("gamma", list(2, 2))

  source <- parameter_source("tau", shape = "row")
  expect_s3_class(source, "parameter_source")
  expect_equal(source$shape, "row")
  expect_false(any(c("index", "expression") %in% names(source)))
  expect_equal(
    BayesTools:::.bt_parameter_source_jags_expression(source, row_index = "i"),
    "tau[i]"
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_label(parameter_source("tau", shape = "row")),
    "tau[row]"
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_row_names(source, 2),
    c("tau[1]", "tau[2]")
  )
  expect_error(
    BayesTools:::.bt_parameter_source_row_names(source, 1.5),
    "integer vector",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_row_names(parameter_source("tau"), 2),
    "not row-shaped",
    fixed = TRUE
  )
  source_values <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$tau * data$scale[seq_len(n_rows)]
    },
    inputs = "tau"
  )
  expect_true(is.function(source_values$values))
  ellipsis_values <- parameter_source(
    "tau",
    shape = "row",
    values = function(...){
      args <- list(...)
      rep(args$parameters$tau, args$n_rows)
    },
    inputs = "tau"
  )
  expect_true(is.function(ellipsis_values$values))
  expect_error(
    parameter_source("tau", shape = "row", values = function(x)x),
    "must accept named arguments",
    fixed = TRUE
  )
  broken_values <- parameter_source("tau", shape = "row")
  broken_values$values <- function(x)x
  expect_error(
    BayesTools:::.bt_check_parameter_source(broken_values),
    "must accept named arguments",
    fixed = TRUE
  )

  sd_source <- random_sd_source(source)
  expect_s3_class(sd_source, "random_sd_source")
  expect_equal(sd_source$shape, "row")
  expect_false(any(c("index", "expression", "row_indexed", "values") %in% names(sd_source)))
  expect_true(BayesTools:::.bt_random_sd_source_is_row_indexed(sd_source))
  expect_equal(
    BayesTools:::.bt_random_sd_source_expression(sd_source, row_index = "i"),
    "tau[i]"
  )
  expect_error(
    BayesTools:::.bt_random_sd_source_expression(sd_source),
    "requires a row index",
    fixed = TRUE
  )
  expect_equal(
    BayesTools:::.bt_random_sd_source_expression(random_sd_source("tau")),
    "tau"
  )
  expect_true(is.function(BayesTools:::.bt_parameter_source_values_function(
    random_sd_source(source_values)
  )))
  source_posterior <- matrix(
    c(2, 4),
    ncol = 1,
    dimnames = list(NULL, "tau")
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_value_draws(
      source = source_values,
      n_rows = 2,
      posterior = source_posterior,
      data = list(scale = c(1, 3))
    ),
    matrix(
      c(2, 6, 4, 12),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("tau[1]", "tau[2]"))
    ),
    tolerance = 1e-12
  )
  dependent_parameters <- BayesTools:::.bt_parameter_source_forbid_formula_parameters(
    list(theta = 2, deterministic_formula = 3, log_sigma = 4),
    "log_sigma"
  )
  allowed_source <- parameter_source(
    "allowed",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters$theta + parameters[["deterministic_formula"]], n_rows)
    },
    inputs = c("theta", "deterministic_formula")
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_value_draws(
      source = allowed_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    matrix(
      rep(5, 4),
      nrow = 2,
      dimnames = list(NULL, c("allowed[1]", "allowed[2]"))
    )
  )
  forbidden_dollar_source <- parameter_source(
    "forbidden_dollar",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters$log_sigma, n_rows)
    },
    inputs = "log_sigma"
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = forbidden_dollar_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    "source 'forbidden_dollar\\[row\\]'.*formula parameter 'log_sigma'",
    perl = TRUE
  )
  forbidden_bracket_source <- parameter_source(
    "forbidden_bracket",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters[["log_sigma"]], n_rows)
    },
    inputs = "log_sigma"
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = forbidden_bracket_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    "source 'forbidden_bracket\\[row\\]'.*formula parameter 'log_sigma'",
    perl = TRUE
  )
  forbidden_subset_source <- parameter_source(
    "forbidden_subset",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters["log_sigma"][[1]], n_rows)
    },
    inputs = "log_sigma"
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = forbidden_subset_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    "source 'forbidden_subset\\[row\\]'.*formula parameter 'log_sigma'",
    perl = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) "bad", inputs = character()),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "must return a numeric vector",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) 1, inputs = character()),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "must return a numeric vector of length 2",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) c(1, NA_real_), inputs = character()),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "returned missing values",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) stop("boom"), inputs = character()),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "failed: boom",
    fixed = TRUE
  )
  source_data <- BayesTools:::.bt_JAGS_marglik_parameter_source_data(
    model_data = list(raw_x = c(1, 2)),
    formula_data = list(raw_x = c(1, 2), formula_only = c(3, 4)),
    design = structure(
      list(
        model_frame = data.frame(raw_x = c(-1, 1), design_only = c(5, 6)),
        source_data = data.frame(raw_x = c(1, 2), source_only = c(5, 6))
      ),
      class = "BayesTools_formula_design"
    )
  )
  expect_equal(source_data$raw_x, c(1, 2))
  expect_equal(source_data$source_only, c(5, 6))
  expect_false("formula_only" %in% names(source_data))
  expect_false("design_only" %in% names(source_data))
  expect_error(
    BayesTools:::.bt_JAGS_marglik_parameter_source_data(
      model_data = list(raw_x = c(9, 9)),
      formula_data = NULL,
      design = structure(
        list(source_data = data.frame(raw_x = c(1, 2))),
        class = "BayesTools_formula_design"
      )
    ),
    "conflict with the fitted source snapshot",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_parameter_source_data(
      model_data = list(raw_x = c(9, 9)),
      formula_data = list(raw_x = c(1, 2)),
      design = NULL
    ),
    "conflicting data for variable 'raw_x'",
    fixed = TRUE
  )
  scaled_source_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = data.frame(
      id = factor(c("a", "a", "b", "b")),
      x = c(10, 20, 30, 40)
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd_source = random_sd_source(parameter_source(
          "tau",
          shape = "row",
          values = function(parameters, data, n_rows){
            data$x[seq_len(n_rows)]
          },
          inputs = character()
        ))
      )
    )
  )
  expect_equal(scaled_source_result$formula_design$source_data$x, c(10, 20, 30, 40))
  expect_false(isTRUE(all.equal(
    scaled_source_result$formula_design$model_frame$x,
    scaled_source_result$formula_design$source_data$x
  )))

  expect_error(
    parameter_source("1tau"),
    "must start with a letter",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", shape = "rows"),
    "should be one of"
  )
  expect_error(
    parameter_source("tau", index = "i"),
    "unused argument",
    fixed = TRUE
  )
  expect_error(
    random_sd_source("tau", shape = "rows"),
    "should be one of"
  )
  expect_error(
    random_sd_source("tau", index = "i"),
    "unused argument",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", values = 1),
    "'values' must be NULL or a function",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", values = function(parameters, data, n_rows) 1),
    "only for row-shaped parameter sources",
    fixed = TRUE
  )
  expect_error(
    parameter_source("mu__xREx__study_intercept"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_sd_source("mu__xRE_ALLOCx_allocation__allocation_sd"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_sd_source("mu__xRE_SUMMARY__sd__study__intercept"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_sd_source(parameter_source("tau"), shape = "row"),
    "must not be supplied",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_sd_binding(
      source = NULL,
      application = "block"
    ),
    "require a shared 'source'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_sd_binding(
      source = NULL,
      application = "column"
    ),
    "require shared or per-column sources",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_sd_binding(
      sources_by_column = list(random_sd_source("tau", shape = "row")),
      application = "column"
    ),
    "Row-indexed per-column SD sources are not supported",
    fixed = TRUE
  )
  broken_sd_source <- BayesTools:::.bt_prior_owned_sd_source("tau")
  broken_sd_source$owned <- FALSE
  expect_error(
    BayesTools:::.bt_check_random_sd_binding_source(broken_sd_source),
    "inconsistent",
    fixed = TRUE
  )
  one_target_factor <- BayesTools:::.bt_random_variance_allocation_factor(
    weight_name = "w",
    index = 1L,
    scale = "total_variance",
    n_targets = 1L
  )
  expect_error(
    BayesTools:::.bt_check_random_variance_allocation_factor(one_target_factor),
    "n_targets",
    fixed = TRUE
  )
  broken_source <- parameter_source("tau")
  broken_source$shape <- "rows"
  expect_error(
    BayesTools:::.bt_check_parameter_source(broken_source),
    "metadata are inconsistent",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug"),
      sd = sd_prior,
      sd_source = random_sd_source("tau")
    ),
    "exactly one of 'sd' or 'sd_source'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug")
    ),
    "exactly one of 'sd' or 'sd_source'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      parent = allocation_ref("total", "study"),
      terms = c("paper", "estimate"),
      sd_source = random_sd_source("tau")
    ),
    "must not specify 'sd_source'",
    fixed = TRUE
  )
})

test_that("source data conflicts name the variable and the remedy, not a refit (F017)", {

  design <- structure(
    list(source_data = data.frame(raw_x = c(1, 2))),
    class = "BayesTools_formula_design"
  )
  error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_parameter_source_data(
      model_data = list(raw_x = c(9, 9)),
      formula_data = NULL,
      design = design
    ),
    error = function(e) e
  )
  expect_s3_class(error, "BayesTools_refit_required")
  expect_match(conditionMessage(error), "variable 'raw_x'", fixed = TRUE)
  expect_match(conditionMessage(error), "must equal the fitted source data", fixed = TRUE)
  expect_match(conditionMessage(error), "Rename the likelihood variable or pass the fitted values.", fixed = TRUE)
  expect_false(grepl("Refit the model", conditionMessage(error), fixed = TRUE))
})

test_that("variance allocation priors generate shared total SD and Dirichlet allocation syntax", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    paper = factor(c("p1", "p2", "p1", "p2", "p3", "p3")),
    drug = factor(c("a", "b", "a", "b", "a", "b")),
    x = c(-1, 0, 1, 2, -2, 3)
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  random_prior <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = random_prior
  )

  expect_equal(
    names(result$prior_list),
    c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation__allocation_sd",
      "mu__xRE_ALLOCx_allocation__weight"
    )
  )
  expect_s3_class(result$prior_list$mu__xRE_ALLOCx_allocation__weight, "prior.simplex")
  expect_match(
    result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )
  expect_match(
    result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  expect_equal(
    result$add_parameters,
    c(
      "mu__xREx__study_intercept",
      "mu__xREx__study_xRE_Zx",
      "mu__xREx__drug_intercept",
      "mu__xREx__drug_xRE_Zx"
    )
  )
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$weight_name,
    "mu__xRE_ALLOCx_allocation__weight"
  )
  expect_equal(result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$index, 1L)
  expect_equal(result$formula_design$random_effects[[2]]$sd_binding$allocations[[1L]]$index, 2L)
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_parameter_names,
    "mu__xREx__study_intercept"
  )

  allocation_inclusion_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "total_re",
        terms = c(study = "study", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3))),
        inclusion = list(study = prior("spike", list(location = 0.5)))
      )
    )
  )
  expect_true("mu__xRE_ALLOCx_total_re__include_study_prob" %in%
                names(allocation_inclusion_result$prior_list))
  expect_match(
    allocation_inclusion_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__include_study_indicator ~ dbern(mu__xRE_ALLOCx_total_re__include_study_prob)",
    fixed = TRUE
  )
  expect_match(
    allocation_inclusion_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * mu__xRE_ALLOCx_total_re__include_study_indicator * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    allocation_inclusion_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[2])",
    fixed = TRUE
  )
  expect_false(grepl(
    "include_drug",
    allocation_inclusion_result$formula_syntax,
    fixed = TRUE
  ))
  expect_true("mu__xRE_ALLOCx_total_re__include_study_indicator" %in%
                allocation_inclusion_result$add_parameters)
  allocation_inclusion_syntax <- JAGS_add_priors(
    "model{}",
    allocation_inclusion_result$prior_list
  )
  expect_match(
    allocation_inclusion_syntax,
    "mu__xRE_ALLOCx_total_re__include_study_prob = 0.5",
    fixed = TRUE
  )

  allocation_inclusion_posterior <- matrix(
    c(
      2, 1, 3, 1,
      2, 1, 3, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]",
        "mu__xRE_ALLOCx_total_re__include_study_indicator"
      )
    )
  )
  inclusion_study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = allocation_inclusion_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = allocation_inclusion_posterior,
    prior_list = allocation_inclusion_result$prior_list
  )
  inclusion_drug_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = allocation_inclusion_result$formula_design$random_effects[[2]],
    n_columns = 1,
    posterior = allocation_inclusion_posterior,
    prior_list = allocation_inclusion_result$prior_list
  )
  expect_equal(inclusion_study_sd[, 1], c(2 * sqrt(1 / 4), 0))
  expect_equal(inclusion_drug_sd[, 1], c(2 * sqrt(3 / 4), 2 * sqrt(3 / 4)))
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = allocation_inclusion_result$formula_design$random_effects[[1]],
      n_columns = 1,
      posterior = allocation_inclusion_posterior[
        ,
        colnames(allocation_inclusion_posterior) !=
          "mu__xRE_ALLOCx_total_re__include_study_indicator",
        drop = FALSE
      ],
      prior_list = allocation_inclusion_result$prior_list
    ),
    "missing Bernoulli indicator",
    fixed = TRUE
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_formula_allocation_inclusion_names(
      list(mu = allocation_inclusion_result$formula_design)
    ),
    "mu__xRE_ALLOCx_total_re__include_study_indicator"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_allocation_factor_metadata(
      allocation_inclusion_result$formula_design$random_effects[[1]]$
        sd_binding$allocations[[1L]]$factors[[1L]]
    )$inclusion_name,
    "mu__xRE_ALLOCx_total_re__include_study_indicator"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_allocation_metadata(
      allocation_inclusion_result$formula_design$random_effects[[1]]$
        sd_binding$allocations[[1L]]
    )$inclusion$study$indicator_name,
    "mu__xRE_ALLOCx_total_re__include_study_indicator"
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_formula_random_prior_evaluator(
      list(mu = allocation_inclusion_result$formula_design)
    ),
    "Bridge sampling for variance allocation inclusion gates",
    fixed = TRUE
  )
  expect_error(
    bayestools_reference_formula_random_log_prior(
      samples = allocation_inclusion_posterior[1, ],
      formula_design_list = list(mu = allocation_inclusion_result$formula_design)
    ),
    "Bridge sampling for variance allocation inclusion gates",
    fixed = TRUE
  )

  total_component_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "a",
        terms = c(total = "nested", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      random_variance_allocation(
        name = "split",
        terms = c("study", "paper"),
        parent = allocation_ref("a", "total"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  expect_match(
    total_component_result$formula_syntax,
    "mu__xRE_ALLOCx_a__component_total_sd = mu__xRE_ALLOCx_a__allocation_sd * sqrt(mu__xRE_ALLOCx_a__weight[1])",
    fixed = TRUE
  )
  expect_false(grepl(
    "mu__xRE_ALLOCx_a__allocation_sd = mu__xRE_ALLOCx_a__allocation_sd",
    total_component_result$formula_syntax,
    fixed = TRUE
  ))

  malformed_binding <- result$formula_design$random_effects[[1]]$sd_binding
  malformed_binding$factors <- list()
  expect_error(
    BayesTools:::.bt_check_random_sd_binding(malformed_binding),
    "binding$factors",
    fixed = TRUE
  )

  reversed_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = c("drug", "study"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_equal(reversed_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$index, 2L)
  expect_equal(reversed_result$formula_design$random_effects[[2]]$sd_binding$allocations[[1L]]$index, 1L)
  expect_match(
    reversed_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  expect_match(
    reversed_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )

  full_syntax <- JAGS_add_priors("model{}", result$prior_list)
  expect_match(
    full_syntax,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1] ~ dgamma(2, 1)",
    fixed = TRUE
  )
  expect_match(
    full_syntax,
    "mu__xRE_ALLOCx_allocation__weight[2] <- prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2] / sum(prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1:2])",
    fixed = TRUE
  )

  posterior <- matrix(
    c(
      2, 1, 3,
      4, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_allocation__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = posterior,
    prior_list = result$prior_list
  )
  drug_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = result$formula_design$random_effects[[2]],
    n_columns = 1,
    posterior = posterior,
    prior_list = result$prior_list
  )
  expect_equal(study_sd[, 1], c(2 * sqrt(1 / 4), 4 * sqrt(2 / 4)))
  expect_equal(drug_sd[, 1], c(2 * sqrt(3 / 4), 4 * sqrt(2 / 4)))
  reversed_study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = reversed_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = posterior,
    prior_list = reversed_result$prior_list
  )
  reversed_drug_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = reversed_result$formula_design$random_effects[[2]],
    n_columns = 1,
    posterior = posterior,
    prior_list = reversed_result$prior_list
  )
  expect_equal(reversed_study_sd[, 1], c(2 * sqrt(3 / 4), 4 * sqrt(2 / 4)))
  expect_equal(reversed_drug_sd[, 1], c(2 * sqrt(1 / 4), 4 * sqrt(2 / 4)))

  bridge_row <- posterior[1, ]
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = bridge_row,
      random_term = result$formula_design$random_effects[[1]],
      prior_list = result$prior_list
    ),
    2 * sqrt(1 / 4)
  )

  single_factor_df <- data.frame(
    f = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("d1", "d2", "d1", "d2"))
  )
  single_factor_result <- JAGS_formula(
    formula = ~ 1 +
      random(0 + f | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = single_factor_df,
    prior_list = fixed_priors,
    prior_random = random_prior
  )
  study_random <- single_factor_result$formula_design$random_effects[[1]]
  expect_equal(study_random$n_columns, 1L)
  expect_equal(study_random$sd_parameter_names, "mu__xREx__study_f[1]")
  expect_equal(study_random$sd_leaves$leaf_names_by_column, study_random$sd_parameter_names)
  expect_match(
    single_factor_result$formula_syntax,
    "mu__xREx__study_f[1] = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )

  random_bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = result$formula_design)
  )
  expect_true(all(grepl("_xRE_Zx", random_bridge$parameters, fixed = TRUE)))
  expect_false(any(grepl("_intercept", random_bridge$parameters, fixed = TRUE)))

  bridge_samples <- c(
    "mu_intercept" = 10,
    "mu__xRE_ALLOCx_allocation__allocation_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__study_xRE_Zx[3,1]" = 0.3,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  fixed_formula_priors <- BayesTools:::.bt_JAGS_marglik_formula_fixed_priors(
    result$prior_list,
    "mu"
  )
  expect_equal(names(fixed_formula_priors), "mu_intercept")

  reconstructed <- JAGS_marglik_parameters_formula(
    samples = bridge_samples,
    formula_list = list(mu = result$formula),
    formula_data_list = list(mu = result$data),
    formula_prior_list = list(mu = result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = result$formula_design)
  )
  expected_mu <- 10 +
    c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
    sqrt(3) * c(1, 2, 1, 2, 1, 2)
  expect_equal(reconstructed$mu, expected_mu, tolerance = 1e-12)

  expect_equal(
    JAGS_marglik_priors_formula(
      samples = bridge_samples,
      formula_prior_list = list(mu = result$prior_list)
    ),
    stats::dnorm(10, 0, 1, log = TRUE) +
      stats::dgamma(2, shape = 2, rate = 2, log = TRUE) +
      stats::dgamma(1, shape = 2, rate = 1, log = TRUE) +
    stats::dgamma(3, shape = 3, rate = 1, log = TRUE),
    tolerance = 1e-12
  )
  edge_bridge_samples <- bridge_samples
  edge_bridge_samples[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]"]] <- 0
  expect_equal(
    JAGS_marglik_priors_formula(
      samples = edge_bridge_samples,
      formula_prior_list = list(mu = result$prior_list)
    ),
    -Inf
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples                 = edge_bridge_samples,
      formula_list            = list(mu = result$formula),
      formula_data_list       = list(mu = result$data),
      formula_prior_list      = list(mu = result$prior_list),
      prior_list_parameters   = list(),
      formula_design_list     = list(mu = result$formula_design)
    ),
    "Bridge samples contain out-of-support positive auxiliary coordinate 'prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]'.",
    fixed = TRUE,
    class = "BayesTools_marglik_out_of_support"
  )
  expect_equal(
    bayestools_reference_formula_random_log_prior(
      bridge_samples,
      list(mu = result$formula_design)
    ),
    sum(stats::dnorm(bridge_samples[grepl("_xRE_Zx", names(bridge_samples), fixed = TRUE)], log = TRUE)),
    tolerance = 1e-12
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = bridge_samples[!grepl("prior_par_eta_", names(bridge_samples), fixed = TRUE)],
      random_term = result$formula_design$random_effects[[1]],
      prior_list = result$prior_list
    ),
    "missing Dirichlet allocation coordinates",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = bridge_samples[names(bridge_samples) != "mu__xRE_ALLOCx_allocation__allocation_sd"],
      random_term = result$formula_design$random_effects[[1]],
      prior_list = result$prior_list
    ),
    "does not contain all monitored formula prior parameters",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        allocation = random_variance_allocation(name = "allocation",
          terms = c("study", "missing"),
          sd = prior("gamma", list(2, 2))
        )
      )
    ),
    "unknown random-effect block",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        study = random_block(sd = prior("gamma", list(4, 1))),
        allocation = random_variance_allocation(name = "allocation",
          terms = c("study", "drug"),
          sd = prior("gamma", list(2, 2))
        )
      )
    ),
    "cannot supply block-specific SD",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 + x | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = random_prior
    ),
    "one SD component",
    fixed = TRUE
  )
})

test_that("bridge positive support is enforced before reconstruction", {

  dirichlet_prior <- prior("dirichlet", list(alpha = c(1, 2)))
  dirichlet_samples <- c(
    "prior_par_eta_w[1]" = 0,
    "prior_par_eta_w[2]" = 2
  )
  expect_equal(JAGS_marglik_priors(dirichlet_samples, list(w = dirichlet_prior)), -Inf)
  expect_error(
    JAGS_marglik_parameters(dirichlet_samples, list(w = dirichlet_prior)),
    "out-of-support positive auxiliary coordinate"
  )

  invgamma_prior <- prior("invgamma", list(1, 2))
  invgamma_samples <- c("sigma" = 0)
  expect_equal(JAGS_marglik_priors(invgamma_samples, list(sigma = invgamma_prior)), -Inf)
  expect_error(
    JAGS_marglik_parameters(invgamma_samples, list(sigma = invgamma_prior)),
    "out-of-support inverse-gamma coordinate"
  )

  formula_prior_list <- list(
    mu = list(mu_intercept = invgamma_prior)
  )
  expect_equal(
    JAGS_marglik_priors_formula(
      samples = c("mu_intercept" = 0),
      formula_prior_list = formula_prior_list
    ),
    -Inf
  )
  intercept_design <- JAGS_formula(
    ~ 1, "mu", data.frame(x = 1),
    list(intercept = invgamma_prior)
  )$formula_design
  expect_error(
    JAGS_marglik_parameters_formula(
      samples = c("mu_intercept" = 0),
      formula_list = list(mu = ~ 1),
      formula_data_list = list(mu = list(N_mu = 1)),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = list(mu = intercept_design)
    ),
    "out-of-support inverse-gamma coordinate"
  )

  # BayesTools 0.3.0 'inv_' precision coordinates are not read
  expect_error(
    JAGS_marglik_parameters(
      c("inv_sigma" = 0),
      list(sigma = invgamma_prior)
    ),
    "'samples' does not contain all monitored inverse-gamma prior parameters.",
    fixed = TRUE
  )
})

test_that("external variance allocation sources generate scalar and row-indexed syntax", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))

  direct_scalar_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      sd = prior("gamma", list(2, 2)),
      study = random_block(sd_source = random_sd_source("tau"))
    )
  )
  expect_equal(names(direct_scalar_result$prior_list), "mu_intercept")
  expect_match(
    direct_scalar_result$formula_syntax,
    "mu__xREx__study_xRE_STDx[1] = tau",
    fixed = TRUE
  )
  expect_equal(
    direct_scalar_result$formula_design$random_effects[[1]]$sd_binding$source$name,
    "tau"
  )
  expect_false(direct_scalar_result$formula_design$random_effects[[1]]$sd_binding$true_allocation)
  expect_true("tau" %in% direct_scalar_result$add_parameters)
  direct_scalar_coordinates <- BayesTools:::.bt_build_parameter_coordinates(
    columns = "tau",
    monitor_names = direct_scalar_result$add_parameters,
    prior_list = direct_scalar_result$prior_list,
    formula_design = list(mu = direct_scalar_result$formula_design)
  )
  expect_identical(direct_scalar_coordinates$role, "parameter")
  expect_identical(direct_scalar_coordinates$display_label, "tau")
  direct_scalar_summary <- .parameter_catalog_random_summary_samples(
    model_samples = matrix(
      2,
      nrow = 1,
      dimnames = list(NULL, "tau")
    ),
    prior_list = direct_scalar_result$prior_list,
    formula_design = list(mu = direct_scalar_result$formula_design)
  )
  expect_false(any(grepl("__xRE_SUMMARY__sd", colnames(direct_scalar_summary$model_samples), fixed = TRUE)))

  direct_row_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      study = random_block(sd_source = random_sd_source("tau", shape = "row"))
    )
  )
  expect_true(all(is.na(direct_row_result$formula_design$random_effects[[1]]$sd_parameter_names)))
  expect_match(
    direct_row_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * 1 * inprod(mu__xREx__study_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_equal(
    direct_row_result$formula_design$random_effects[[1]]$sd_binding$source$shape,
    "row"
  )

  expect_error(
    random_block(
      sd = prior("gamma", list(2, 2)),
      sd_source = random_sd_source("tau")
    ),
    "'sd_source' cannot be supplied together with block-local 'sd'",
    fixed = TRUE
  )

  scalar_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  expect_equal(
    names(scalar_result$prior_list),
    c("mu_intercept", "mu__xRE_ALLOCx_allocation__weight")
  )
  expect_match(
    scalar_result$formula_syntax,
    "mu__xREx__study_intercept = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )
  expect_match(
    scalar_result$formula_syntax,
    "mu__xREx__drug_intercept = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  expect_equal(
    scalar_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$source$name,
    "tau"
  )
  expect_equal(
    scalar_result$formula_design$random_effects[[1]]$sd_binding$source$shape,
    "scalar"
  )

  scalar_posterior <- matrix(
    c(
      2, 1, 3,
      4, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "tau",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  scalar_study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = scalar_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = scalar_posterior,
    prior_list = scalar_result$prior_list
  )
  expect_equal(scalar_study_sd[, 1], c(2 * sqrt(1 / 4), 4 * sqrt(2 / 4)))
  scalar_monitored_leaf_posterior <- matrix(
    c(
      0.5, 1.5, 1, 3,
      0.75, 1.25, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xREx__study_intercept",
        "mu__xREx__drug_intercept",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  scalar_study_leaf_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = scalar_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = scalar_monitored_leaf_posterior,
    prior_list = scalar_result$prior_list
  )
  expect_equal(scalar_study_leaf_sd[, 1], c(0.5, 0.75))
  scalar_bridge_samples <- c(
    "tau" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = scalar_bridge_samples,
      random_term = scalar_result$formula_design$random_effects[[1]],
      prior_list = scalar_result$prior_list
    ),
    c(0.1, 0.1, 0.2, 0.2),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_scale_metadata(scalar_result$formula_design$random_effects[[1]])$sd_binding$source$shape,
    "scalar"
  )

  scalar_sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_equal(
    names(scalar_sd_leaf_result$prior_list),
    c("mu_intercept", "mu__xRE_ALLOCx_allocation__weight")
  )
  expect_match(
    scalar_sd_leaf_result$formula_syntax,
    "mu__xREx__study_intercept = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )
  expect_match(
    scalar_sd_leaf_result$formula_syntax,
    "mu__xREx__study_x = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  scalar_sd_leaf_posterior <- matrix(
    c(2, 1, 3),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "tau",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = scalar_sd_leaf_result$formula_design$random_effects[[1]],
      n_columns = 2,
      posterior = scalar_sd_leaf_posterior,
      prior_list = scalar_sd_leaf_result$prior_list
    )[1, ],
    c(2 * sqrt(1 / 4), 2 * sqrt(3 / 4)),
    tolerance = 1e-12
  )

  row_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  expect_equal(
    names(row_result$prior_list),
    c("mu_intercept", "mu__xRE_ALLOCx_allocation__weight")
  )
  expect_true(all(is.na(row_result$formula_design$random_effects[[1]]$sd_parameter_names)))
  expect_equal(
    row_result$formula_design$random_effects[[1]]$sd_binding$source$shape,
    "row"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_scale_metadata(row_result$formula_design$random_effects[[1]])$sd_binding$source$shape,
    "row"
  )
  row_bridge_metadata <- BayesTools:::.bt_JAGS_bridge_scale_metadata(
    row_result$formula_design$random_effects[[1]]
  )
  expect_equal(row_bridge_metadata$sd_binding$source$shape, "row")
  expect_equal(row_bridge_metadata$sd_binding$source$source$shape, "row")
  expect_false(any(c("index", "expression", "row_indexed", "values") %in% names(row_bridge_metadata$sd_binding$source)))
  expect_equal(row_bridge_metadata$sd_binding$allocations[[1L]]$source$shape, "row")
  expect_false(any(c("source_index", "source_row_indexed") %in% names(row_bridge_metadata$sd_binding$allocations[[1L]])))
  expect_false(grepl("mu__xREx__study_intercept = tau[i]", row_result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__study_xRE_STDx", row_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * inprod(mu__xREx__study_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_equal(
    row_result$add_parameters,
    c("mu__xREx__study_xRE_Zx", "mu__xREx__drug_xRE_Zx")
  )
  row_id_result <- JAGS_formula(
    formula = ~ 1 +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_false(grepl("mu__xREx__study_xRE_STDx", row_id_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_id_result$formula_syntax,
    "for(g in 1:2){\n     mu__xREx__study_xRE_UNIT_COEFx[g,i] = mu__xREx__study_xRE_Zx[g,i]",
    fixed = TRUE
  )
  expect_match(
    row_id_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * inprod(mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i], 1:2]",
    fixed = TRUE
  )
  row_cs_result <- JAGS_formula(
    formula = ~ 1 +
      cs(f | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, f = factor(c("a", "b", "a", "b"))),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      study = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_false(grepl("mu__xREx__study_xRE_STDx", row_cs_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_cs_result$formula_syntax,
    "mu__xREx__study_xRE_UNIT_COEFx[g,i] <- mu__xREx__study_xRE_CS_PREFIXx[g,i - 1] + mu__xREx__study_xRE_CS_DIAGx[i] * mu__xREx__study_xRE_Zx[g,i]",
    fixed = TRUE
  )
  expect_match(
    row_cs_result$formula_syntax,
    paste0(
      "mu__xREx__study[i] = tau[i] * ",
      "sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * ",
      "mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i],",
      "mu__xREx__study_xRE_COLx[i]]"
    ),
    fixed = TRUE
  )
  expect_null(row_cs_result$data$mu__xREx__study_xRE_DATAx)

  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = scalar_posterior[1, -1],
      random_term = row_result$formula_design$random_effects[[1]],
      prior_list = row_result$prior_list
    ),
    "scalar random-effect SD values are not defined for row-indexed sources",
    fixed = TRUE
  )
  row_bridge_samples <- c(
    "mu_intercept" = 10,
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 1,
    "mu__xREx__study_xRE_Zx[2,1]" = 2,
    "mu__xREx__drug_xRE_Zx[1,1]" = 3,
    "mu__xREx__drug_xRE_Zx[2,1]" = 4
  )
  row_bridge_posterior <- matrix(
    row_bridge_samples,
    nrow = 1,
    dimnames = list(NULL, names(row_bridge_samples))
  )
  attr(row_bridge_posterior, "lb") <- stats::setNames(
    rep(-Inf, ncol(row_bridge_posterior)),
    colnames(row_bridge_posterior)
  )
  attr(row_bridge_posterior, "ub") <- stats::setNames(
    rep(Inf, ncol(row_bridge_posterior)),
    colnames(row_bridge_posterior)
  )
  attr(row_bridge_posterior, "lb")[paste0("tau[", 1:4, "]")] <- 0
  expect_silent(BayesTools:::.bt_JAGS_bridge_check_row_indexed_external_sd_sources(
    formula_design_list = list(mu = row_result$formula_design),
    bridgesampling_posterior = row_bridge_posterior
  ))
  bad_row_bridge_posterior <- row_bridge_posterior
  attr(bad_row_bridge_posterior, "lb")["tau[2]"] <- -Inf
  expect_error(
    BayesTools:::.bt_JAGS_bridge_check_row_indexed_external_sd_sources(
      formula_design_list = list(mu = row_result$formula_design),
      bridgesampling_posterior = bad_row_bridge_posterior
    ),
    "non-negative lower bounds",
    fixed = TRUE
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_bridge_samples,
      random_term = row_result$formula_design$random_effects[[1]],
      prior_list = row_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * 1,
      4 * sqrt(1 / 4) * 1,
      6 * sqrt(1 / 4) * 2,
      8 * sqrt(1 / 4) * 2
    ),
    tolerance = 1e-12
  )
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = row_bridge_samples,
      formula_list = list(mu = row_result$formula),
      formula_data_list = list(mu = row_result$data),
      formula_prior_list = list(mu = row_result$prior_list),
      prior_list_parameters = list(),
      formula_design_list = list(mu = row_result$formula_design)
    )$mu,
    10 + c(
      2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
      4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
      6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
      8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
    ),
    tolerance = 1e-12
  )
  value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$total_tau * data$tau_factor[seq_len(n_rows)]
    },
    inputs = "total_tau"
  )
  value_prior_random <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd_source = random_sd_source(value_source),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )
  row_values_data <- transform(
    df,
    tau_factor = c(1, 2, 3, 4)
  )
  row_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = row_values_data,
    prior_list = fixed_priors,
    prior_random = value_prior_random
  )
  row_value_samples <- row_bridge_samples[
    !grepl("^tau\\[", names(row_bridge_samples))
  ]
  row_expected_value <- 10 + c(
    2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
    4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
    6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
    8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
  )
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = row_value_samples,
      formula_list = list(mu = row_values_result$formula),
      formula_data_list = list(mu = row_values_result$data),
      formula_prior_list = list(mu = row_values_result$prior_list),
      prior_list_parameters = list(total_tau = 2),
      formula_design_list = list(mu = row_values_result$formula_design),
      model_data = list(tau_factor = c(1, 2, 3, 4))
    )$mu,
    row_expected_value,
    tolerance = 1e-12
  )
  total_tau_result <- JAGS_formula(
    formula = ~ 1,
    parameter = "total_tau",
    data = df,
    prior_list = fixed_priors
  )
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = c(row_value_samples, "total_tau_intercept" = 2),
      formula_list = list(
        mu = row_values_result$formula,
        total_tau = total_tau_result$formula
      ),
      formula_data_list = list(
        mu = row_values_result$data,
        total_tau = total_tau_result$data
      ),
      formula_prior_list = list(
        mu = row_values_result$prior_list,
        total_tau = total_tau_result$prior_list
      ),
      prior_list_parameters = list(),
      formula_design_list = list(
        mu = row_values_result$formula_design,
        total_tau = total_tau_result$formula_design
      ),
      model_data = list(tau_factor = c(1, 2, 3, 4))
    )$mu,
    row_expected_value,
    tolerance = 1e-12
  )
  bridge_preflight_fit <- matrix(
    row_value_samples,
    nrow = 1,
    dimnames = list(NULL, names(row_value_samples))
  )
  bridge_preflight_fit <- bridge_preflight_fit[, names(row_value_samples) != "mu_intercept", drop = FALSE]
  bridge_preflight_fit <- cbind(
    "mu_intercept" = 0,
    bridge_preflight_fit
  )
  bridge_preflight_fit <- coda::mcmc(bridge_preflight_fit)
  row_values_result$formula_design$formula_scale <- .bt_formula_scale_finalize(
    row_values_result$formula_scale, row_values_result$formula_design, row_values_result$prior_list)
  attr(bridge_preflight_fit, "formula_design") <- list(
    mu = row_values_result$formula_design
  )
  bridge_value_context <- BayesTools:::.bt_JAGS_bridge_formula_context(
    fit = bridge_preflight_fit,
    formula_list = list(mu = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag")),
    formula_data_list = list(mu = row_values_data),
    formula_prior_list = list(mu = fixed_priors),
    formula_scale_list = NULL,
    formula_random_prior_list = list(mu = value_prior_random)
  )
  expect_true(BayesTools:::.bt_parameter_source_has_values(
    BayesTools:::.bt_random_effect_row_indexed_source(
      bridge_value_context$formula_design_list$mu$random_effects[[1]]
    )
  ))
  expect_null(bridge_value_context$formula_data_list)
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = row_value_samples,
      formula_list = bridge_value_context$formula_list,
      formula_data_list = NULL,
      formula_prior_list = bridge_value_context$formula_prior_list,
      prior_list_parameters = list(total_tau = 2),
      formula_design_list = bridge_value_context$formula_design_list
    )$mu,
    row_expected_value,
    tolerance = 1e-12
  )
  ambiguous_row_value_samples <- c(
    row_value_samples,
    "tau[1]" = 999,
    "tau[2]" = 999,
    "tau[3]" = 999,
    "tau[4]" = 999
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = ambiguous_row_value_samples,
      random_term = row_values_result$formula_design$random_effects[[1]],
      prior_list = row_values_result$prior_list
    ),
    "Use one row-source reconstruction path only",
    fixed = TRUE
  )
  ambiguous_bridge_posterior <- matrix(
    ambiguous_row_value_samples,
    nrow = 1,
    dimnames = list(NULL, names(ambiguous_row_value_samples))
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_check_row_indexed_external_sd_sources(
      formula_design_list = list(mu = row_values_result$formula_design),
      bridgesampling_posterior = ambiguous_bridge_posterior
    ),
    "also appears as bridge parameter column",
    fixed = TRUE
  )
  bad_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      c(1, -1, 1, 1)[seq_len(n_rows)]
    },
    inputs = character()
  )
  bad_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(bad_value_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  row_value_posterior <- matrix(
    row_value_samples,
    nrow = 1,
    dimnames = list(NULL, names(row_value_samples))
  )
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_contribution_from_latent(
      random_term = bad_values_result$formula_design$random_effects[[1]],
      model_matrix = bad_values_result$formula_design$random_effects[[1]]$model_matrix,
      group_map = bad_values_result$formula_design$random_effects[[1]]$group_map,
      posterior = row_value_posterior,
      prior_list = bad_values_result$prior_list,
      context = "Prediction"
    ),
    "non-finite or negative source values",
    fixed = TRUE
  )
  negative_marglik_error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_value_samples,
      random_term = bad_values_result$formula_design$random_effects[[1]],
      prior_list = bad_values_result$prior_list
    ),
    error = function(e)e
  )
  expect_s3_class(negative_marglik_error, "BayesTools_marglik_out_of_support")
  nan_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      c(1, NaN, 1, 1)[seq_len(n_rows)]
    },
    inputs = character()
  )
  nan_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(nan_value_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  nan_marglik_error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_value_samples,
      random_term = nan_values_result$formula_design$random_effects[[1]],
      prior_list = nan_values_result$prior_list
    ),
    error = function(e)e
  )
  expect_s3_class(nan_marglik_error, "BayesTools_marglik_out_of_support")
  short_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      1
    },
    inputs = character()
  )
  short_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(short_value_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_value_samples,
      random_term = short_values_result$formula_design$random_effects[[1]],
      prior_list = short_values_result$prior_list
    ),
    "must return a numeric vector of length 4",
    fixed = TRUE
  )
  missing_source_fit <- bridge_preflight_fit
  row_result$formula_design$formula_scale <- .bt_formula_scale_finalize(
    row_result$formula_scale, row_result$formula_design, row_result$prior_list)
  attr(missing_source_fit, "formula_design") <- list(
    mu = row_result$formula_design
  )
  expect_error(
    JAGS_bridgesampling(
      fit = missing_source_fit,
      log_posterior = function(parameters, data) 0,
      data = list()
    ),
    "cannot reconstruct row-indexed external SD source 'tau[row]'",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        allocation = random_variance_allocation(name = "allocation",
          sd_source = random_sd_source("tau", shape = "row"),
          weights = prior("dirichlet", list(alpha = c(2, 3)))
        ),
        monitor = random_monitor(coefficients = TRUE)
      )
    ),
    "Group-level coefficient monitoring is not supported",
    fixed = TRUE
  )

  row_nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, paper = factor(c("p1", "p2", "p1", "p2"))),
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 3)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        weights = prior("dirichlet", list(alpha = c(3, 2)))
      )
    )
  )
  expect_false(grepl(
    "mu__xRE_ALLOCx_total_re__component_nested_sd = tau[i]",
    row_nested_result$formula_syntax,
    fixed = TRUE
  ))
  expect_match(
    row_nested_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1]) * inprod(mu__xREx__study_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_match(
    row_nested_result$formula_syntax,
    "mu__xREx__drug[i] = tau[i] * sqrt(mu__xRE_ALLOCx_total_re__weight[2]) * inprod(mu__xREx__drug_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_equal(length(row_nested_result$formula_design$random_effects[[1]]$sd_binding$factors), 2L)
  expect_equal(length(row_nested_result$formula_design$random_effects[[3]]$sd_binding$factors), 1L)
  row_nested_samples <- c(
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[1]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[2]" = 2,
    "mu__xREx__study_xRE_Zx[1,1]" = 1,
    "mu__xREx__study_xRE_Zx[2,1]" = 2,
    "mu__xREx__paper_xRE_Zx[1,1]" = 3,
    "mu__xREx__paper_xRE_Zx[2,1]" = 4,
    "mu__xREx__drug_xRE_Zx[1,1]" = 5,
    "mu__xREx__drug_xRE_Zx[2,1]" = 6
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_nested_samples,
      random_term = row_nested_result$formula_design$random_effects[[1]],
      prior_list = row_nested_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * sqrt(3 / 5) * 1,
      4 * sqrt(1 / 4) * sqrt(3 / 5) * 1,
      6 * sqrt(1 / 4) * sqrt(3 / 5) * 2,
      8 * sqrt(1 / 4) * sqrt(3 / 5) * 2
    ),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_nested_samples,
      random_term = row_nested_result$formula_design$random_effects[[2]],
      prior_list = row_nested_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * sqrt(2 / 5) * 3,
      4 * sqrt(1 / 4) * sqrt(2 / 5) * 4,
      6 * sqrt(1 / 4) * sqrt(2 / 5) * 3,
      8 * sqrt(1 / 4) * sqrt(2 / 5) * 4
    ),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_nested_samples,
      random_term = row_nested_result$formula_design$random_effects[[3]],
      prior_list = row_nested_result$prior_list
    ),
    c(
      2 * sqrt(3 / 4) * 5,
      4 * sqrt(3 / 4) * 6,
      6 * sqrt(3 / 4) * 5,
      8 * sqrt(3 / 4) * 6
    ),
    tolerance = 1e-12
  )

  row_sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  expect_true(all(is.na(row_sd_leaf_result$formula_design$random_effects[[1]]$sd_parameter_names)))
  expect_equal(row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$application, "column")
  expect_equal(length(row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$factors_by_column), 2L)
  expect_false(grepl("mu__xREx__study_intercept = tau[i]", row_sd_leaf_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_sd_leaf_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * (sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i],1] * mu__xREx__study_xRE_DATAx[i,1] + sqrt(mu__xRE_ALLOCx_allocation__weight[2]) * mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i],2] * mu__xREx__study_xRE_DATAx[i,2])",
    fixed = TRUE
  )
  row_sd_leaf_forced_latent_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      study = random_block(
        monitor = random_monitor(latent = FALSE, coefficients = FALSE)
      )
    )
  )
  expect_true("mu__xREx__study_xRE_Zx" %in%
                row_sd_leaf_forced_latent_result$add_parameters)
  expect_false(any(grepl(
    "mu__xREx__study_xRE_COEFx",
    row_sd_leaf_forced_latent_result$add_parameters,
    fixed = TRUE
  )))
  row_sd_leaf_samples <- c(
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 1,
    "mu__xREx__study_xRE_Zx[2,1]" = 2,
    "mu__xREx__study_xRE_Zx[1,2]" = 3,
    "mu__xREx__study_xRE_Zx[2,2]" = 4
  )
  row_sd_leaf_expected <- c(
    2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * -1),
    4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * 0),
    6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 1),
    8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 2)
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_samples,
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_result$prior_list
    ),
    row_sd_leaf_expected,
    tolerance = 1e-12
  )
  row_sd_leaf_ambiguous_samples <- c(
    row_sd_leaf_samples,
    "mu__xRE_ALLOCx_allocation__weight[1]" = 0.25,
    "mu__xRE_ALLOCx_allocation__weight[2]" = 0.75
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_ambiguous_samples,
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_result$prior_list
    ),
    "both normalized Dirichlet allocation coordinates",
    fixed = TRUE
  )
  expect_equal(
    unname(BayesTools:::.bt_random_effect_row_indexed_contribution_from_latent(
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      model_matrix = row_sd_leaf_result$formula_design$random_effects[[1]]$model_matrix,
      group_map = row_sd_leaf_result$formula_design$random_effects[[1]]$group_map,
      posterior = matrix(row_sd_leaf_samples, nrow = 1, dimnames = list(NULL, names(row_sd_leaf_samples))),
      prior_list = row_sd_leaf_result$prior_list
    )[, 1]),
    row_sd_leaf_expected,
    tolerance = 1e-12
  )
  row_sd_leaf_bad_eta <- row_sd_leaf_samples
  row_sd_leaf_bad_eta[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]"]] <- 0
  row_sd_leaf_bad_eta[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"]] <- 0
  row_sd_leaf_bad_eta_error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_bad_eta,
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_result$prior_list
    ),
    error = function(e)e
  )
  expect_s3_class(row_sd_leaf_bad_eta_error, "BayesTools_marglik_out_of_support")
  valid_weight_draws <- matrix(
    c(0.25, 0.75),
    nrow = 1,
    dimnames = list(NULL, c("w[1]", "w[2]"))
  )
  expect_equal(
    BayesTools:::.bt_random_effect_dirichlet_draws(
      parameter_name = "w",
      posterior = valid_weight_draws,
      prior_list = list(w = prior("dirichlet", list(alpha = c(1, 1))))
    ),
    valid_weight_draws
  )
  invalid_weight_error <- tryCatch(
    BayesTools:::.bt_random_effect_dirichlet_draws(
      parameter_name = "w",
      posterior = matrix(
        c(0.25, 0.90),
        nrow = 1,
        dimnames = list(NULL, c("w[1]", "w[2]"))
      ),
      prior_list = list(w = prior("dirichlet", list(alpha = c(1, 1))))
    ),
    error = function(e)e
  )
  expect_s3_class(
    invalid_weight_error,
    "BayesTools_random_effect_allocation_out_of_support"
  )
  expect_identical(class(invalid_weight_error), c("BayesTools_random_effect_allocation_out_of_support",
    "BayesTools_random_effects_error", "error", "condition"))
  expect_identical(invalid_weight_error$call, NULL)
  expect_type(invalid_weight_error$message, "character")
  expect_true(BayesTools:::.bt_JAGS_marglik_random_effect_support_error(invalid_weight_error))
  bad_row_sd_leaf_target <- row_sd_leaf_result$formula_design$random_effects[[1]]
  bad_row_sd_leaf_target$sd_binding$allocations[[1L]]$target <- "block"
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = bad_row_sd_leaf_target,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "missing 'allocation\\$factors'"
  )
  missing_allocation_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  missing_allocation_term$sd_binding$allocations <- list()
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = missing_allocation_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "exactly one allocation record",
    fixed = TRUE
  )
  extra_allocation_binding <- row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding
  extra_allocation_binding$allocations <- c(
    extra_allocation_binding$allocations,
    extra_allocation_binding$allocations
  )
  expect_error(
    BayesTools:::.bt_check_random_sd_binding(extra_allocation_binding),
    "exactly one allocation record",
    fixed = TRUE
  )
  stale_allocation_binding <- row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding
  stale_allocation_binding$true_allocation <- FALSE
  expect_error(
    BayesTools:::.bt_check_random_sd_binding(stale_allocation_binding),
    "must not contain allocation records",
    fixed = TRUE
  )
  bad_row_sd_leaf_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  bad_row_sd_leaf_term$sd_binding$application <- "block"
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = bad_row_sd_leaf_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "Block SD bindings must not use 'factors_by_column'",
    fixed = TRUE
  )
  empty_column_chain_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  empty_column_chain_term$sd_binding$factors_by_column <- list()
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = empty_column_chain_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "missing canonical 'binding\\$factors_by_column'"
  )
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_contribution_from_latent(
      random_term = empty_column_chain_term,
      model_matrix = empty_column_chain_term$model_matrix,
      group_map = empty_column_chain_term$group_map,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list
    ),
    "missing canonical 'binding\\$factors_by_column'"
  )
  mismatched_column_chain_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  mismatched_column_chain_term$sd_binding$factors_by_column[[1L]][[1L]]$index <- 2L
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = mismatched_column_chain_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "column factor chains do not match",
    fixed = TRUE
  )
  row_sd_leaf_mean_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "mean_variance",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_samples,
      random_term = row_sd_leaf_mean_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_mean_result$prior_list
    ),
    c(
      2 * (sqrt(2 * 1 / 4) * 1 + sqrt(2 * 3 / 4) * 3 * -1),
      4 * (sqrt(2 * 1 / 4) * 1 + sqrt(2 * 3 / 4) * 3 * 0),
      6 * (sqrt(2 * 1 / 4) * 2 + sqrt(2 * 3 / 4) * 4 * 1),
      8 * (sqrt(2 * 1 / 4) * 2 + sqrt(2 * 3 / 4) * 4 * 2)
    ),
    tolerance = 1e-12
  )

  row_child_sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "het", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(het = "het", simple = "drug"),
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 3)))
      ),
      random_variance_allocation(
        name = "het_sd",
        parent = allocation_ref("total_re", "het"),
        terms = "het",
        target = "sd_component",
        weights = prior("dirichlet", list(alpha = c(1, 3)))
      )
    )
  )
  expect_false(grepl("mu__xRE_ALLOCx_total_re__component_het_sd = tau[i]", row_child_sd_leaf_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_child_sd_leaf_result$formula_syntax,
    "mu__xREx__het[i] = tau[i] * (sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_het_sd__weight[1]) * mu__xREx__het_xRE_UNIT_COEFx[mu__xREx__het_xRE_MAPx[i],1]",
    fixed = TRUE
  )
  expect_equal(length(row_child_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$factors_by_column[[1L]]), 2L)
  row_child_sd_leaf_samples <- c(
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[2]" = 3,
    "mu__xREx__het_xRE_Zx[1,1]" = 1,
    "mu__xREx__het_xRE_Zx[2,1]" = 2,
    "mu__xREx__het_xRE_Zx[1,2]" = 3,
    "mu__xREx__het_xRE_Zx[2,2]" = 4,
    "mu__xREx__drug_xRE_Zx[1,1]" = 5,
    "mu__xREx__drug_xRE_Zx[2,1]" = 6
  )
  row_child_sd_leaf_expected <- c(
    2 * sqrt(1 / 4) * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * -1),
    4 * sqrt(1 / 4) * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * 0),
    6 * sqrt(1 / 4) * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 1),
    8 * sqrt(1 / 4) * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 2)
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_child_sd_leaf_samples,
      random_term = row_child_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_child_sd_leaf_result$prior_list
    ),
    row_child_sd_leaf_expected,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_child_sd_leaf_samples,
      random_term = row_child_sd_leaf_result$formula_design$random_effects[[2]],
      prior_list = row_child_sd_leaf_result$prior_list
    ),
    c(
      2 * sqrt(3 / 4) * 5,
      4 * sqrt(3 / 4) * 6,
      6 * sqrt(3 / 4) * 5,
      8 * sqrt(3 / 4) * 6
    ),
    tolerance = 1e-12
  )
  row_child_bridge_design <- list(mu = row_child_sd_leaf_result$formula_design)
  row_child_bridge_changed <- row_child_bridge_design
  row_child_bridge_changed_term <- row_child_bridge_changed$mu$random_effects[[1]]
  row_child_bridge_changed_allocation <- row_child_bridge_changed_term$sd_binding$allocations[[1L]]
  row_child_bridge_changed_allocation$parent_factors[[1L]]$index <- 2L
  row_child_bridge_changed_term$sd_binding$allocations[[1L]] <- row_child_bridge_changed_allocation
  row_child_bridge_changed$mu$random_effects[[1]] <- row_child_bridge_changed_term
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = row_child_bridge_design[["mu"]],
      rebuilt = row_child_bridge_changed[["mu"]]
    ),
    "binding$factors",
    fixed = TRUE
  )
})

test_that("variance allocation graph supports child block and SD-leaf allocations", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))

  df_nested <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    paper = factor(c("p1", "p2", "p1", "p2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  nested_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(nested = "nested", drug = "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      weights = prior("dirichlet", list(alpha = c(3, 2)))
    )
  )
  nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = nested_prior
  )

  expect_equal(
    names(nested_result$prior_list),
    c(
      "mu_intercept",
      "mu__xRE_ALLOCx_total_re__allocation_sd",
      "mu__xRE_ALLOCx_total_re__weight",
      "mu__xRE_ALLOCx_nested_split__weight"
    )
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__component_nested_sd = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[2])",
    fixed = TRUE
  )
  expect_equal(length(nested_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors), 2L)
  expect_equal(length(nested_result$formula_design$random_effects[[3]]$sd_binding$allocations[[1L]]$factors), 1L)

  nested_gate_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(nested = "nested", drug = "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1))),
      inclusion = list(nested = prior("spike", list(location = 0.5)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      weights = prior("dirichlet", list(alpha = c(3, 2)))
    )
  )
  nested_gate_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = nested_gate_prior
  )
  expect_match(
    nested_gate_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__component_nested_sd = mu__xRE_ALLOCx_total_re__allocation_sd * mu__xRE_ALLOCx_total_re__include_nested_indicator * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_gate_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * mu__xRE_ALLOCx_total_re__include_nested_indicator * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1])",
    fixed = TRUE
  )
  expect_true("mu__xRE_ALLOCx_total_re__include_nested_indicator" %in%
                nested_gate_result$add_parameters)
  nested_gate_posterior <- matrix(
    c(
      4, 1, 3, 3, 2, 1,
      4, 1, 3, 3, 2, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[2]",
        "mu__xRE_ALLOCx_total_re__include_nested_indicator"
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_gate_result$formula_design$random_effects[[1]],
      n_columns = 1,
      posterior = nested_gate_posterior,
      prior_list = nested_gate_result$prior_list
    )[, 1],
    c(4 * sqrt(1 / 4) * sqrt(3 / 5), 0),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_gate_result$formula_design$random_effects[[3]],
      n_columns = 1,
      posterior = nested_gate_posterior,
      prior_list = nested_gate_result$prior_list
    )[, 1],
    c(4 * sqrt(3 / 4), 4 * sqrt(3 / 4)),
    tolerance = 1e-12
  )

  nested_external_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(nested = "nested", drug = "drug"),
      sd_source = random_sd_source("tau"),
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      weights = prior("dirichlet", list(alpha = c(3, 2)))
    )
  )
  nested_external_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = nested_external_prior
  )
  expect_equal(
    names(nested_external_result$prior_list),
    c(
      "mu_intercept",
      "mu__xRE_ALLOCx_total_re__weight",
      "mu__xRE_ALLOCx_nested_split__weight"
    )
  )
  expect_match(
    nested_external_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__component_nested_sd = tau * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_external_result$formula_syntax,
    "mu__xREx__study_intercept = tau * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1])",
    fixed = TRUE
  )
  expect_equal(
    nested_external_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$source$kind,
    "external"
  )
  expect_equal(length(nested_external_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors), 2L)

  nested_posterior <- matrix(
    c(4, 1, 3, 3, 2),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[2]"
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_result$formula_design$random_effects[[1]],
      n_columns = 1,
      posterior = nested_posterior,
      prior_list = nested_result$prior_list
    )[, 1],
    4 * sqrt(1 / 4) * sqrt(3 / 5),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_result$formula_design$random_effects[[3]],
      n_columns = 1,
      posterior = nested_posterior,
      prior_list = nested_result$prior_list
    )[, 1],
    4 * sqrt(3 / 4),
    tolerance = 1e-12
  )

  df_sd <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    x = c(-1, 0, 1, 2)
  )
  sd_leaf_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(het = "het", simple = "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "het_sd",
      parent = allocation_ref("total_re", "het"),
      terms = "het",
      target = "sd_component",
      scale = "mean_variance",
      weights = prior("dirichlet", list(alpha = c(2, 2)))
    )
  )
  sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "het", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_sd,
    prior_list = fixed_priors,
    prior_random = sd_leaf_prior
  )

  expect_match(
    sd_leaf_result$formula_syntax,
    "mu__xREx__het_intercept = mu__xRE_ALLOCx_total_re__component_het_sd * sqrt(2 * mu__xRE_ALLOCx_het_sd__weight[1])",
    fixed = TRUE
  )
  expect_match(
    sd_leaf_result$formula_syntax,
    "mu__xREx__het_x = mu__xRE_ALLOCx_total_re__component_het_sd * sqrt(2 * mu__xRE_ALLOCx_het_sd__weight[2])",
    fixed = TRUE
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$target,
    "sd_component"
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$leaf_terms,
    c(mu__xREx__het_intercept = "intercept", mu__xREx__het_x = "x")
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_leaves$leaf_names,
    names(sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$leaf_terms)
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_leaves$leaf_index_by_column,
    c(1L, 2L)
  )

  sd_leaf_samples <- c(
    "mu__xRE_ALLOCx_total_re__allocation_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[2]" = 3
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = sd_leaf_samples,
      random_term = sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = sd_leaf_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * sqrt(2 * 1 / 4),
      2 * sqrt(1 / 4) * sqrt(2 * 3 / 4)
    ),
    tolerance = 1e-12
  )
})

test_that("variance allocation graph rejects ambiguous or conflicting specifications", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    x = c(-1, 0, 1, 2)
  )

  expect_error(
    random_variance_allocation(
      name = "bad_child",
      parent = allocation_ref("total_re", "study"),
      terms = "study",
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "must not specify 'sd'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      name = "bad_mean",
      terms = c("study", "drug"),
      sd = sd_prior,
      scale = "mean_variance"
    ),
    "target = \"sd_component\"",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      name = "bad_sd_terms",
      terms = c("study", "drug"),
      sd = sd_prior,
      target = "sd_component"
    ),
    "exactly one",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1))),
          inclusion = list(missing = prior("spike", list(location = 0.5)))
        )
      )
    ),
    "unknown component label",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "study_sd",
          terms = "study",
          sd = sd_prior,
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1))),
          inclusion = list(study = prior("spike", list(location = 0.5)))
        )
      )
    ),
    "supports target = \"block\"",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "either all named or all unnamed",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(name = "allocation",
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "second",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "without explicit 'terms'",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "duplicate",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "duplicate",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "labels must be unique",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(symbolic = "symbolic", drug = "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "Unknown targets are allowed only when consumed by a child allocation",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "child_first",
          parent = allocation_ref("total_re", "study"),
          terms = "study",
          target = "sd_component",
          scale = "mean_variance",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "total_re",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "must be defined before child allocation",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "bad_ref",
          parent = allocation_ref("total_re", "missing"),
          terms = "study",
          target = "sd_component",
          scale = "mean_variance",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "does not contain component 'missing'",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split_1",
          parent = allocation_ref("total_re", "study"),
          terms = "study",
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split_2",
          parent = allocation_ref("total_re", "study"),
          terms = "drug",
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "can be consumed by only one child",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split",
          parent = allocation_ref("total_re", "study"),
          terms = c("study", "drug"),
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "also names a random-effect block",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_variance_allocation_component_labels(c("a-b", "a_b")),
    "unique after sanitization",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "study_sd",
          terms = "study",
          sd = sd_prior,
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1, 1)))
        )
      )
    ),
    "dimension must match",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "study_sd",
          terms = "study",
          sd = sd_prior,
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "at least two resolved SD components",
    fixed = TRUE
  )
})
