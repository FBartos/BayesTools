skip_if_not_test_profile("unit")

.re_compile_sd_prior <- function(){
  prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
}

.re_compile_data <- function(){
  data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"),
                      levels = c("e1", "e2", "e3", "e4")),
    x = c(0, 1, 0, 1)
  )
}

.re_compile_formula <- function(){
  ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | estimate, name = "estimate", covariance = "diag")
}

.re_compile_prior_random <- function(monitor = random_monitor()){
  sd_prior <- .re_compile_sd_prior()
  prior_random(
    study = random_block(sd = sd_prior, monitor = monitor),
    estimate = random_block(sd = sd_prior, monitor = monitor)
  )
}

.re_compile_result <- function(compile = NULL, prior_random = .re_compile_prior_random()){
  JAGS_formula(
    formula = .re_compile_formula(),
    parameter = "mu",
    data = .re_compile_data(),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random,
    random_effects_compile = compile
  )
}

.re_compile_block_names <- function(random_effects){
  vapply(random_effects, `[[`, character(1), "block_name")
}

.re_compile_known_group_kernel <- function(){
  levels <- c("s2", "extra", "s1", "s3")
  K <- diag(c(1.0, 2.5, 1.5, 2.0))
  K[lower.tri(K)] <- c(.20, .10, .05, .15, .08, .12)
  K <- K + t(K) - diag(diag(K))
  dimnames(K) <- list(levels, levels)
  K
}

.re_compile_terms_by_mode <- function(design, mode){
  design$random_effects[
    vapply(design$random_effects, function(random_term){
      identical(random_term$compile_mode, mode)
    }, logical(1))
  ]
}

test_that("random_effects_compile validates block requests", {

  compile <- random_effects_compile(marginalized = "estimate")
  expect_true(inherits(compile, "random_effects_compile"))
  expect_equal(compile$marginalized, "estimate")
  expect_named(compile, "marginalized")
  expect_false("sampled" %in% names(formals(random_effects_compile)))

  expect_error(
    random_effects_compile(marginalized = 1),
    "must be a character vector",
    fixed = TRUE
  )
  expect_error(
    random_effects_compile(marginalized = c("study", "study")),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    random_effects_compile(marginalized = NA_character_),
    "cannot contain NA",
    fixed = TRUE
  )
  malformed_resolved <- BayesTools:::.bt_random_effects_compile_resolved(
    sampled = "study",
    marginalized = "study",
    mode = c(study = "sampled")
  )
  expect_error(
    BayesTools:::.bt_check_random_effects_compile(malformed_resolved),
    "cannot be both sampled and marginalized",
    fixed = TRUE
  )
  malformed_mode <- BayesTools:::.bt_random_effects_compile_resolved(
    sampled = "study",
    marginalized = character(),
    mode = c(study = "invalid")
  )
  expect_error(
    BayesTools:::.bt_check_random_effects_compile(malformed_mode),
    "must be 'sampled' or 'marginalized'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_effects_compile_mode(malformed_mode, "study"),
    "must be 'sampled' or 'marginalized'",
    fixed = TRUE
  )
  inconsistent_mode <- BayesTools:::.bt_random_effects_compile_resolved(
    sampled = "study",
    marginalized = "estimate",
    mode = c(study = "marginalized", estimate = "sampled")
  )
  expect_error(
    BayesTools:::.bt_check_random_effects_compile(inconsistent_mode),
    "inconsistent with the sampled and marginalized blocks",
    fixed = TRUE
  )
  expect_error(
    .re_compile_result(random_effects_compile(marginalized = "missing")),
    "unknown random-effect block",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1,
      parameter = "mu",
      data = .re_compile_data(),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      random_effects_compile = random_effects_compile(marginalized = "missing")
    ),
    "unknown random-effect block",
    fixed = TRUE
  )
  expect_false(exists(
    "random_effects_marginal_variance_expression",
    envir = asNamespace("BayesTools"),
    inherits = FALSE
  ))
})

test_that("default random-effect compilation remains all sampled", {

  result <- .re_compile_result()
  design <- result$formula_design

  expect_equal(.re_compile_block_names(design$random_effects), c("study", "estimate"))
  expect_false("random_effects_all" %in% names(design))
  expect_false("marginalized_random_effects" %in% names(design))
  expect_equal(design$random_effects_compile$sampled, c("study", "estimate"))
  expect_null(design$random_effects_compile$marginalized)
  expect_length(.re_compile_terms_by_mode(design, "marginalized"), 0L)
  expect_equal(
    design$random_effects_compile$mode,
    c(study = "sampled", estimate = "sampled")
  )
  expect_match(result$formula_syntax, "mu__xREx__study_xRE_Zx", fixed = TRUE)
  expect_match(result$formula_syntax, "mu__xREx__estimate_xRE_Zx", fixed = TRUE)
  expect_match(result$formula_syntax, "mu__xREx__estimate\\[i\\] =")
})

test_that("known group covariance compiles as sampled random-intercept kernel", {

  data <- data.frame(
    study = factor(c("s2", "s1", "s3", "s2"), levels = c("s2", "s1", "s3"))
  )
  K <- .re_compile_known_group_kernel()
  random_effects <- random_effects_formula(
    ~ 1 | study,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  result <- JAGS_formula(
    formula = random_effects,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(
        sd = .re_compile_sd_prior(),
        monitor = random_monitor(latent = TRUE, coefficients = TRUE)
      )
    )
  )
  random_term <- result$formula_design$random_effects[[1]]
  expected_levels <- levels(data$study)
  expected_kernel <- K[expected_levels, expected_levels]

  expect_s3_class(random_term$group_covariance, "random_group_covariance_kernel")
  expect_equal(random_term$group_covariance$levels, expected_levels)
  expect_equal(random_term$group_covariance$dropped_levels, "extra")
  expect_equal(random_term$group_covariance$kernel, expected_kernel)
  expect_equal(
    result$data[["mu__xREx__study_xRE_GROUP_PRECx"]],
    chol2inv(chol(expected_kernel)),
    tolerance = 1e-12
  )
  expect_equal(
    result$data[["mu__xREx__study_xRE_GROUP_MUx"]],
    rep(0, length(expected_levels))
  )
  expect_match(
    result$formula_syntax,
    "mu__xREx__study_xRE_GROUP_Zx[1:3,j] ~ dmnorm",
    fixed = TRUE
  )
  expect_match(
    result$formula_syntax,
    "mu__xREx__study_xRE_Zx[i,j] = mu__xREx__study_xRE_GROUP_Zx[i,j]",
    fixed = TRUE
  )
  expect_false(grepl("mu__xREx__study_xRE_PRECx", result$formula_syntax, fixed = TRUE))
  expect_true("mu__xREx__study_xRE_Zx" %in% result$add_parameters)
  expect_true("mu__xREx__study_xRE_COEFx" %in% result$add_parameters)
})

test_that("known group covariance compiles as marginalized random-intercept metadata", {

  data <- data.frame(
    study = factor(c("s2", "s1", "s3", "s2"), levels = c("s2", "s1", "s3"))
  )
  K <- .re_compile_known_group_kernel()
  random_effects <- random_effects_formula(
    ~ 1 | study,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  result <- JAGS_formula(
    formula = random_effects,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(
        sd = .re_compile_sd_prior(),
        monitor = random_monitor(latent = TRUE, coefficients = TRUE)
      )
    ),
    random_effects_compile = random_effects_compile(marginalized = "study")
  )
  random_term <- result$formula_design$random_effects[[1]]
  expected_levels <- levels(data$study)
  expected_kernel <- K[expected_levels, expected_levels]

  expect_true(all(c(
    "block_name", "group_label", "group_levels", "group_map", "model_matrix",
    "n_groups", "n_columns", "structure", "homogeneous_sd",
    "sd_parameter_names", "sd_binding", "compile_mode", "group_covariance"
  ) %in% names(random_term)))
  expect_equal(random_term$compile_mode, "marginalized")
  expect_s3_class(random_term$group_covariance, "random_group_covariance_kernel")
  expect_equal(random_term$group_covariance$levels, expected_levels)
  expect_equal(random_term$group_covariance$dropped_levels, "extra")
  expect_equal(random_term$group_covariance$kernel, expected_kernel)
  expect_equal(random_term$group_map, c(1, 2, 3, 1))
  expect_equal(dim(random_term$model_matrix), c(4L, 1L))
  expect_equal(random_term$n_groups, 3L)
  expect_equal(random_term$n_columns, 1L)
  expect_equal(random_term$jags_data_names, character())

  expect_false("mu__xREx__study_xRE_GROUP_PRECx" %in% names(result$data))
  expect_false("mu__xREx__study_xRE_GROUP_MUx" %in% names(result$data))
  expect_false("mu__xREx__study_xRE_DATAx" %in% names(result$data))
  expect_false("mu__xREx__study_xRE_MAPx" %in% names(result$data))
  expect_false(grepl("mu__xREx__study_xRE_GROUP_Zx", result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__study_xRE_Zx", result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__study_xRE_COEFx", result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__study[i] =", result$formula_syntax, fixed = TRUE))
  expect_false(any(grepl("mu__xREx__study_xRE_Zx", result$add_parameters, fixed = TRUE)))
  expect_false(any(grepl("mu__xREx__study_xRE_COEFx", result$add_parameters, fixed = TRUE)))
})

test_that("known group covariance supports exact noncentered random slopes", {

  data <- data.frame(
    study = factor(c("s2", "s1", "s3", "s2"), levels = c("s2", "s1", "s3")),
    x = c(0, 1, 0, 1)
  )
  K <- .re_compile_known_group_kernel()
  random_intercept <- random_effects_formula(
    ~ 1 | study,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  random_slope <- random_effects_formula(
    ~ x | study,
    group_covariance = random_group_covariance(K, scale = "none")
  )

  slope_prior <- prior_random(
    study = random_block(
      sd = .re_compile_sd_prior(),
      cor = prior_lkj(eta = 2, include_primitives = TRUE)
    )
  )
  sampled_slope <- JAGS_formula(
    formula = random_slope,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = slope_prior
  )
  marginalized_slope <- JAGS_formula(
    formula = random_slope,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = slope_prior,
    random_effects_compile = random_effects_compile(marginalized = "study")
  )

  sampled_term <- sampled_slope$formula_design$random_effects[[1L]]
  marginalized_term <- marginalized_slope$formula_design$random_effects[[1L]]
  expect_equal(sampled_term$n_columns, 2L)
  expect_identical(sampled_term$parameterization_resolved, "noncentered")
  expect_s3_class(sampled_term$group_covariance, "random_group_covariance_kernel")
  expect_match(
    sampled_slope$formula_syntax,
    "mu__xREx__study_xRE_GROUP_Zx[1:3,j] ~ dmnorm",
    fixed = TRUE
  )
  expect_match(sampled_slope$formula_syntax, "for(j in 1:2)", fixed = TRUE)
  expect_equal(marginalized_term$n_columns, 2L)
  expect_identical(marginalized_term$compile_mode, "marginalized")
  expect_s3_class(
    marginalized_term$group_covariance,
    "random_group_covariance_kernel"
  )
  expect_false(grepl(
    "mu__xREx__study_xRE_GROUP_Zx",
    marginalized_slope$formula_syntax,
    fixed = TRUE
  ))

  expect_error(
    JAGS_formula(
      formula = random_slope,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        study = random_block(
          sd = .re_compile_sd_prior(),
          cor = prior_lkj(eta = 2, include_primitives = TRUE),
          parameterization = "centered"
        )
      )
    ),
    "requires the exact noncentered parameterization",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = random_intercept,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        study = random_block(sd_source = random_sd_source("tau", shape = "row"))
      )
    ),
    "row-indexed external SD sources",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = random_intercept,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        study = random_block(sd_source = random_sd_source("tau", shape = "row"))
      ),
      random_effects_compile = random_effects_compile(marginalized = "study")
    ),
    "row-indexed external SD sources",
    fixed = TRUE
  )
})

test_that("marginalized blocks keep metadata but omit latent mean nodes", {

  result <- .re_compile_result(
    random_effects_compile(marginalized = "estimate")
  )
  design <- result$formula_design

  sampled_terms <- .re_compile_terms_by_mode(design, "sampled")
  marginalized_terms <- .re_compile_terms_by_mode(design, "marginalized")
  expect_equal(.re_compile_block_names(design$random_effects), c("study", "estimate"))
  expect_equal(.re_compile_block_names(sampled_terms), "study")
  expect_equal(.re_compile_block_names(marginalized_terms), "estimate")
  expect_equal(design$random_effects_compile$sampled, "study")
  expect_equal(design$random_effects_compile$marginalized, "estimate")
  expect_equal(
    design$random_effects_compile$mode,
    c(study = "sampled", estimate = "marginalized")
  )
  expect_equal(marginalized_terms[[1]]$compile_mode, "marginalized")
  expect_equal(marginalized_terms[[1]]$sd_parameter_names,
               "mu__xREx__estimate_intercept")
  expect_equal(marginalized_terms[[1]]$jags_data_names, character())
  expect_equal(design$jags_data_names[["__xREx__estimate"]], character())
  expect_false("mu__xREx__estimate_xRE_DATAx" %in% names(result$data))
  expect_false("mu__xREx__estimate_xRE_MAPx" %in% names(result$data))
  expect_equal(dim(marginalized_terms[[1]]$model_matrix), c(4L, 1L))
  expect_equal(marginalized_terms[[1]]$group_map, seq_len(4L))
  expect_true(all(.re_compile_terms_by_mode(design, "sampled")[[1]]$jags_data_names %in% names(result$data)))

  expect_false(grepl("mu__xREx__estimate_xRE_Zx", result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__estimate_xRE_COEFx", result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__estimate[i] =", result$formula_syntax, fixed = TRUE))
  expect_match(result$formula_syntax, "mu__xREx__estimate_xRE_STDx", fixed = TRUE)
  expect_match(result$formula_syntax, "mu__xREx__study_xRE_Zx", fixed = TRUE)
})

test_that("marginalized blocks retain allocation and correlation metadata", {

  allocation_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(study = "study", estimate = "estimate"),
      sd = .re_compile_sd_prior(),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )
  allocation_result <- .re_compile_result(
    random_effects_compile(marginalized = "estimate"),
    prior_random = allocation_prior
  )
  marginalized <- .re_compile_terms_by_mode(
    allocation_result$formula_design,
    "marginalized"
  )[[1]]
  expect_true(isTRUE(marginalized$sd_binding$true_allocation))
  expect_equal(marginalized$sd_binding$allocations[[1L]]$label, "total_re")
  expect_equal(marginalized$sd_binding$allocations[[1L]]$terms,
               c(study = "study", estimate = "estimate"))
  expect_true("mu__xRE_ALLOCx_total_re__weight" %in%
                names(allocation_result$prior_list))

  hcs_data <- data.frame(
    time = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(time | id),
    parameter = "mu",
    data = hcs_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = .re_compile_sd_prior(),
        rho = prior("normal", list(mean = 0, sd = 0.5))
      )
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  hcs_term <- .re_compile_terms_by_mode(
    hcs_result$formula_design,
    "marginalized"
  )[[1]]
  expect_equal(hcs_term$structure, "hcs")
  expect_equal(hcs_term$correlation$type, "rho")
  expect_equal(length(hcs_term$sd_parameter_names), 2L)
  expect_false(grepl("mu__xREx__id_xRE_Zx", hcs_result$formula_syntax, fixed = TRUE))
  expect_match(hcs_result$formula_syntax, hcs_term$correlation$rho_name, fixed = TRUE)
  expect_null(hcs_term$correlation$cholesky_name)
  expect_null(hcs_term$correlation$correlation_name)
})

test_that("nested formula expansion resolves marginalized block names", {

  df <- data.frame(
    district = factor(c("d1", "d1", "d2", "d2")),
    school = factor(c("s1", "s2", "s3", "s4"))
  )
  random_effects <- random_effects_formula(
    list(nested = ~ 1 | district / school)
  )
  result <- JAGS_formula(
    formula = random_effects$formula,
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      nested_school_district = random_block(sd = .re_compile_sd_prior()),
      nested_district = random_block(sd = .re_compile_sd_prior())
    ),
    random_effects_compile = random_effects_compile(
      marginalized = "nested_school_district"
    )
  )

  expect_equal(
    .re_compile_block_names(result$formula_design$random_effects),
    c("nested_school_district", "nested_district")
  )
  expect_equal(
    .re_compile_block_names(.re_compile_terms_by_mode(
      result$formula_design,
      "marginalized"
    )),
    "nested_school_district"
  )
  expect_equal(
    .re_compile_block_names(.re_compile_terms_by_mode(
      result$formula_design,
      "sampled"
    )),
    "nested_district"
  )
})

test_that("formula evaluation adds sampled random effects only", {

  result <- .re_compile_result(
    random_effects_compile(marginalized = "estimate"),
    prior_random = .re_compile_prior_random(
      monitor = random_monitor(latent = FALSE, coefficients = TRUE)
    )
  )
  posterior <- matrix(
    c(
      10, 1, 3,
      20, 2, 4
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__study_xRE_COEFx[1,1]",
      "mu__xREx__study_xRE_COEFx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = .re_compile_formula(),
    parameter = "mu",
    data = .re_compile_data(),
    prior_list = result$prior_list
  )

  expect_equal(
    unname(prediction),
    cbind(c(11, 11, 13, 13), c(22, 22, 24, 24))
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + random(1 | missing, name = "missing", covariance = "diag"),
      parameter = "mu",
      data = transform(.re_compile_data(), missing = study),
      prior_list = result$prior_list
    ),
    "were not found in the fitted formula",
    fixed = TRUE
  )
})

test_that("bridge rebuild validation preserves compile policy", {

  result <- .re_compile_result(
    random_effects_compile(marginalized = "estimate")
  )
  all_sampled <- .re_compile_result()
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      list(mu = result$formula_design),
      list(mu = all_sampled$formula_design)
    ),
    "compile metadata",
    fixed = TRUE
  )

  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list(NULL, "dummy")
  ))
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  expect_error(
    BayesTools:::.bt_JAGS_bridge_formula_context(
      fit = fit,
      formula_list = list(mu = .re_compile_formula()),
      formula_data_list = list(mu = .re_compile_data()),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
      formula_scale_list = NULL,
      formula_random_prior_list = list(mu = .re_compile_prior_random()),
      formula_random_effects_compile_list = NULL
    ),
    "compile metadata",
    fixed = TRUE
  )

  context <- BayesTools:::.bt_JAGS_bridge_formula_context(
    fit = fit,
    formula_list = list(mu = .re_compile_formula()),
    formula_data_list = list(mu = .re_compile_data()),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
    formula_scale_list = NULL,
    formula_random_prior_list = list(mu = .re_compile_prior_random()),
    formula_random_effects_compile_list = list(
      mu = random_effects_compile(marginalized = "estimate")
    )
  )
  expect_equal(
    .re_compile_block_names(.re_compile_terms_by_mode(
      context$formula_design_list$mu,
      "marginalized"
    )),
    "estimate"
  )

  expect_error(
    BayesTools:::.bt_JAGS_bridge_formula_context(
      fit = fit,
      formula_list = list(mu = .re_compile_formula()),
      formula_data_list = list(mu = .re_compile_data()),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
      formula_scale_list = NULL,
      formula_random_prior_list = list(mu = .re_compile_prior_random()),
      formula_random_effects_compile_list = list(mu = list(marginalized = "estimate"))
    ),
    "random_effects_compile",
    fixed = TRUE
  )
})

test_that("bridge prior helpers use all structural metadata and sampled latent effects only", {

  hcs_data <- data.frame(
    time = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(time | id),
    parameter = "mu",
    data = hcs_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = .re_compile_sd_prior(),
        rho = prior("normal", list(mean = 0, sd = 0.5))
      )
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  rho_bridge <- BayesTools:::.bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    list(mu = hcs_result$formula_design)
  )
  expect_equal(rho_bridge$parameters, "mu__xREx__id_rho_z")

  samples <- c("mu__xREx__id_rho_z" = 0)
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      samples,
      list(mu = hcs_result$formula_design)
    ),
    0
  )
  samples_bad <- c("mu__xREx__id_rho_z" = NA_real_)
  expect_identical(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      samples_bad,
      list(mu = hcs_result$formula_design)
    ),
    -Inf
  )
})

test_that("bridge formula reconstruction ignores all-marginalized blocks in the mean", {

  result <- .re_compile_result(
    random_effects_compile(marginalized = c("study", "estimate"))
  )
  samples <- c(mu_intercept = 5)

  parameters <- JAGS_marglik_parameters_formula(
    samples = samples,
    formula_list = list(mu = result$formula),
    formula_data_list = list(mu = .re_compile_data()),
    formula_prior_list = list(mu = result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = result$formula_design)
  )

  expect_equal(parameters$mu, rep(5, nrow(.re_compile_data())))
})

test_that("JAGS_extend preserves marginalized random-effect metadata", {

  skip_if_not_installed("runjags")

  result <- .re_compile_result(
    random_effects_compile(marginalized = "estimate")
  )
  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "prior_list") <- list()
  attr(fit, "model_syntax") <- "model{}"
  attr(fit, "required_packages") <- character()
  attr(fit, "jags_modules") <- character()
  attr(fit, "add_parameters") <- character()
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_registry") <- build_test_parameter_registry(character())

  extended <- JAGS_extend(
    fit,
    autofit_control = list(
      max_Rhat = 1.05,
      min_ESS = 1,
      max_error = 1,
      max_SD_error = 1,
      max_time = list(time = 0, unit = "secs"),
      sample_extend = 1,
      restarts = 1,
      max_extend = 1,
      check_indicators = FALSE
    ),
    silent = TRUE
  )

  expect_equal(
    attr(extended, "formula_design")$mu$random_effects_compile$mode,
    c(study = "sampled", estimate = "marginalized")
  )
  expect_equal(
    .re_compile_block_names(.re_compile_terms_by_mode(
      attr(extended, "formula_design")$mu,
      "marginalized"
    )),
    "estimate"
  )
})
