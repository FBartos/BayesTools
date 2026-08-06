skip_if_not_test_profile("unit")

.parameter_catalog_test_fit <- function(chains, prior_list,
                                        formula_design = NULL){

  fit <- chains
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- prior_list
  if(!is.null(formula_design)){
    attr(fit, "formula_design") <- formula_design
  }
  attr(fit, "parameter_registry") <- .bt_build_parameter_registry(
    columns = colnames(chains[[1L]]),
    prior_list = prior_list,
    formula_design = formula_design
  )
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_parameter_catalog(fit)
  .bt_attach_fit_contract(fit)
}

test_that("parameter catalog construction is metadata-only and versioned", {

  prior_list <- list(
    theta = prior("normal", list(0, 1)),
    fixed = prior("point", list(-3))
  )
  registry <- .bt_build_parameter_registry(
    columns = "theta",
    prior_list = prior_list
  )
  testthat::local_mocked_bindings(
    .extract_posterior_samples = function(...){
      stop("posterior extraction is forbidden")
    },
    .package = "BayesTools"
  )

  expect_silent(catalog <- .bt_build_parameter_catalog(
    registry = registry,
    prior_list = prior_list
  ))
  expect_s3_class(catalog, "BayesTools_parameter_catalog")
  expect_identical(catalog$schema_version, 1L)
  expect_identical(
    names(catalog$quantities),
    .bt_parameter_catalog_quantity_columns
  )
  fixed <- catalog$quantities[catalog$quantities$canonical_name == "fixed", ]
  expect_identical(fixed$status, "structural")
  expect_identical(fixed$fixed_value, -3)
  expect_identical(fixed$extraction_key[[1L]]$dependencies, "fixed")
})

test_that("factor catalog components preserve fitted level identities", {

  data <- data.frame(f = factor(c("a", "b", "c")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + f,
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  registry <- .bt_build_parameter_registry(
    columns = c("mu_intercept", "mu_f[1]", "mu_f[2]"),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    registry = registry,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )

  factor_rows <- catalog$quantities[
    catalog$quantities$canonical_name %in% c("mu_f[1]", "mu_f[2]"),
    ,
    drop = FALSE
  ]
  expect_identical(factor_rows$component, c("b", "c"))
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      alias = "f",
      namespace = "mu",
      component = "b"
    )$quantities$canonical_name,
    "mu_f[1]"
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      alias = "f[c]",
      namespace = "mu"
    )$quantities$canonical_name,
    "mu_f[2]"
  )

  resolved <- hypothesis_resolve(
    hypothesis_parse("f[b] > f[c]"),
    catalog,
    namespace = "mu"
  )
  occurrence_map <- unique(resolved$occurrences[
    c("symbol", "canonical_name", "component")
  ])
  expect_identical(
    occurrence_map$canonical_name,
    c("mu_f[1]", "mu_f[2]")
  )
  expect_identical(occurrence_map$component, c("b", "c"))
  expect_error(
    hypothesis_resolve(
      hypothesis_parse("f[b] > 0"),
      catalog,
      namespace = "mu",
      component = "c"
    ),
    "does not match the requested catalog component",
    fixed = TRUE
  )
})

test_that("catalog extensions preserve ambiguity until filtered", {

  registry <- .bt_build_parameter_registry(
    columns = "theta",
    prior_list = list(theta = prior("normal", list(0, 1)))
  )
  catalog <- .bt_build_parameter_catalog(registry)
  location <- .bt_parameter_catalog_quantity(
    canonical_name = "effect_location",
    namespace = "location",
    role = "coefficient",
    extraction_key = list(type = "robma", dependencies = character())
  )
  scale <- .bt_parameter_catalog_quantity(
    canonical_name = "effect_scale",
    namespace = "scale",
    role = "coefficient",
    extraction_key = list(type = "robma", dependencies = character())
  )
  quantities <- rbind(location, scale)
  quantities$provider <- "RoBMA"
  quantities$quantity_id <- c("RoBMA::location", "RoBMA::scale")
  aliases <- data.frame(
    alias = c("effect", "effect"),
    quantity_id = quantities$quantity_id,
    namespace = quantities$namespace,
    component = c("", ""),
    stringsAsFactors = FALSE
  )

  extended <- parameter_catalog_extend(
    catalog,
    quantities = quantities,
    aliases = aliases,
    provider = "RoBMA"
  )
  expect_error(
    parameter_catalog_resolve(extended, "effect"),
    class = "BayesTools_parameter_ambiguous"
  )
  resolved <- parameter_catalog_resolve(
    extended,
    "effect",
    namespace = "scale"
  )
  expect_identical(resolved$quantity_id, "RoBMA::scale")
  expect_error(
    parameter_catalog_resolve(extended, "missing"),
    class = "BayesTools_parameter_not_found"
  )
  theta_id <- catalog$quantities$quantity_id[
    catalog$quantities$canonical_name == "theta"
  ]
  alias_only <- parameter_catalog_extend(
    extended,
    quantities = .bt_parameter_catalog_empty_quantities(),
    aliases = data.frame(
      alias = "shared_theta",
      quantity_id = theta_id,
      namespace = "model",
      component = "",
      stringsAsFactors = FALSE
    ),
    provider = "RoBMA"
  )
  expect_identical(
    parameter_catalog_resolve(alias_only, "shared_theta")$quantity_id,
    theta_id
  )
  expect_identical(unserialize(serialize(extended, NULL)), extended)
  expect_error(
    parameter_catalog_extend(
      extended,
      quantities = quantities[1L, , drop = FALSE],
      aliases = aliases[1L, , drop = FALSE],
      provider = "RoBMA"
    ),
    "already exist"
  )
})

test_that("ordinary catalog extraction requests only the selected coordinate", {

  chains <- coda::mcmc.list(coda::mcmc(
    matrix(
      c(1, 10, 2, 20, 3, 30),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("theta", "phi"))
    ),
    start = 7,
    thin = 2
  ))
  fit <- .parameter_catalog_test_fit(
    chains,
    prior_list = list(
      theta = prior("normal", list(0, 1)),
      phi = prior("normal", list(0, 1))
    )
  )
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "theta")
  draws <- parameter_draws(fit, selection)

  expect_identical(colnames(draws[[1L]]), "theta")
  expect_identical(as.numeric(draws[[1L]][, "theta"]), c(1, 2, 3))
  expect_identical(attr(draws[[1L]], "mcpar"), c(7, 11, 2))
})

test_that("random summaries are cataloged and extracted from declared dependencies", {

  data <- data.frame(
    x = 1:4,
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj()
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  cholesky <- random_term$correlation$cholesky_name
  columns <- c(
    "mu_intercept",
    random_term$sd_parameter_names,
    paste0(cholesky, "[1,1]"),
    paste0(cholesky, "[2,1]"),
    paste0(cholesky, "[2,2]")
  )
  values <- cbind(
    0,
    1:3,
    4:6,
    1,
    c(0, .2, .4),
    sqrt(1 - c(0, .2, .4)^2)
  )
  colnames(values) <- columns
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values, start = 11)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )

  catalog <- parameter_catalog(fit)
  derived <- catalog$quantities[catalog$quantities$status == "derived", ]
  expect_setequal(
    derived$role,
    c("random_sd", "random_correlation")
  )
  sd_quantity <- derived[derived$role == "random_sd", , drop = FALSE][1L, ]
  expect_false("mu_intercept" %in%
                 sd_quantity$extraction_key[[1L]]$dependencies)

  observed <- NULL
  original <- .bt_parameter_draw_dependencies
  testthat::local_mocked_bindings(
    .bt_parameter_draw_dependencies = function(fit, dependencies){
      observed <<- dependencies
      original(fit, dependencies)
    },
    .package = "BayesTools"
  )
  selection <- parameter_catalog_resolve(
    catalog,
    sd_quantity$canonical_name
  )
  draws <- parameter_draws(fit, selection)

  expect_identical(observed, sd_quantity$extraction_key[[1L]]$dependencies)
  expect_identical(colnames(draws[[1L]]), sd_quantity$canonical_name)
  expect_identical(as.numeric(draws[[1L]][, 1L]), c(1, 2, 3))
})

test_that("declared variance allocations have metadata-only catalog rows", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  random_terms <- formula_result$formula_design$random_effects
  columns <- c(
    "mu_intercept",
    "mu__xRE_ALLOCx_allocation__total_sd",
    "mu__xRE_ALLOCx_allocation__weight[1]",
    "mu__xRE_ALLOCx_allocation__weight[2]",
    unlist(lapply(random_terms, `[[`, "sd_parameter_names"))
  )
  registry <- .bt_build_parameter_registry(
    columns = columns,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    registry = registry,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  derived <- catalog$quantities[catalog$quantities$status == "derived", ]

  expect_setequal(
    derived$role,
    c("random_sd_total", "random_sd", "random_var_frac")
  )
  study_sd_label <- "(mu) sd(intercept | study)"
  study_sd <- parameter_catalog_resolve(
    catalog,
    alias = study_sd_label,
    namespace = "mu"
  )
  expect_identical(study_sd$quantities$status, "derived")
  sampled_study_sd <- catalog$quantities[
    catalog$quantities$status == "sampled" &
      catalog$quantities$display_label == study_sd_label,
    ,
    drop = FALSE
  ]
  expect_identical(nrow(sampled_study_sd), 1L)
  sampled_selection <- parameter_catalog_resolve(
    catalog,
    alias = sampled_study_sd$canonical_name,
    namespace = "mu"
  )
  expect_identical(sampled_selection$quantities$status, "sampled")
  fractions <- derived[derived$role == "random_var_frac", ]
  expect_identical(fractions$component, c("study", "drug"))
  expect_true(all(vapply(fractions$extraction_key, function(key){
    all(c(
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]"
    ) %in% key$dependencies)
  }, logical(1))))
  study_fraction <- parameter_catalog_resolve(
    catalog,
    alias = "var_frac(allocation: study)",
    namespace = "mu"
  )
  expect_identical(
    study_fraction$quantities$canonical_name,
    fractions$canonical_name[fractions$component == "study"]
  )
})

test_that("malformed catalogs and stale selections fail closed", {

  registry <- .bt_build_parameter_registry(columns = "theta")
  catalog <- .bt_build_parameter_catalog(registry)
  selection <- parameter_catalog_resolve(catalog, "theta")

  broken <- catalog
  broken$schema_version <- 2L
  expect_error(
    .bt_validate_parameter_catalog(broken),
    "Refit or rebuild"
  )

  stale <- selection
  stale$quantity_id <- "BayesTools::different"
  expect_error(
    .bt_validate_parameter_selection(stale, catalog),
    "disagree"
  )
})
