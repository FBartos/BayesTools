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
  expect_identical(catalog$schema_version, 2L)
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
  reference <- parameter_catalog_resolve(
    catalog,
    alias = "f[a]",
    namespace = "mu"
  )$quantities
  expect_identical(reference$status, "structural")
  expect_identical(reference$fixed_value, 0)
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

test_that("factor catalog quantities reconstruct fitted term-level cells", {

  data <- data.frame(f = factor(c("a", "b", "c")))
  build_fit <- function(factor_prior_input, values){
    formula_result <- JAGS_formula(
      ~ 1 + f,
      "mu",
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        f = factor_prior_input
      )
    )
    factor_prior <- formula_result$prior_list$mu_f
    coordinates <- .JAGS_prior_factor_names("mu_f", factor_prior)
    chains <- coda::mcmc.list(coda::mcmc(
      cbind(mu_intercept = 0, values),
      start = 5L,
      thin = 2L
    ))
    colnames(chains[[1L]]) <- c("mu_intercept", coordinates)
    list(
      fit = .parameter_catalog_test_fit(
        chains,
        prior_list = formula_result$prior_list,
        formula_design = list(mu = formula_result$formula_design)
      ),
      prior = factor_prior
    )
  }

  treatment <- build_fit(
    prior_factor("normal", list(0, 1), contrast = "treatment"),
    matrix(c(1, 10, 2, 20), ncol = 2L, byrow = TRUE)
  )
  treatment_catalog <- parameter_catalog(treatment$fit)
  reference <- parameter_catalog_resolve(
    treatment_catalog,
    "f[a]",
    namespace = "mu"
  )
  reference_draws <- parameter_draws(treatment$fit, reference)
  expect_identical(as.numeric(reference_draws[[1L]][, 1L]), c(0, 0))
  expect_identical(attr(reference_draws[[1L]], "mcpar"), c(5, 7, 2))
  direct <- parameter_catalog_resolve(
    treatment_catalog,
    "f[b]",
    namespace = "mu"
  )
  expect_identical(direct$quantities$canonical_name, "mu_f[1]")
  expect_identical(
    as.numeric(parameter_draws(treatment$fit, direct)[[1L]][, 1L]),
    c(1, 2)
  )

  independent <- build_fit(
    prior_factor("normal", list(0, 1), contrast = "independent"),
    matrix(c(1, 10, 100, 2, 20, 200), ncol = 3L, byrow = TRUE)
  )
  independent_catalog <- parameter_catalog(independent$fit)
  for(level_i in seq_along(levels(data$f))){
    selection <- parameter_catalog_resolve(
      independent_catalog,
      paste0("f[", levels(data$f)[level_i], "]"),
      namespace = "mu"
    )
    expect_identical(
      selection$quantities$canonical_name,
      paste0("mu_f[", level_i, "]")
    )
  }

  for(contrast in c("orthonormal", "meandif")){
    transformed <- build_fit(
      prior_factor("mnormal", list(0, 1), contrast = contrast),
      matrix(c(1, 10, 2, 20), ncol = 2L, byrow = TRUE)
    )
    design <- .factor_term_design_from_metadata(transformed$prior)$design
    catalog <- parameter_catalog(transformed$fit)
    for(level_i in seq_along(levels(data$f))){
      selection <- parameter_catalog_resolve(
        catalog,
        paste0("f[", levels(data$f)[level_i], "]"),
        namespace = "mu"
      )
      expect_identical(selection$quantities$status, "derived")
      observed <- as.numeric(
        parameter_draws(transformed$fit, selection)[[1L]][, 1L]
      )
      expected <- as.vector(
        matrix(c(1, 10, 2, 20), ncol = 2L, byrow = TRUE) %*%
          design[level_i, ]
      )
      expect_equal(observed, expected, info = contrast)
    }
  }

  wrapped_priors <- list(
    spike_and_slab = prior_spike_and_slab(
      prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    ),
    mixture = prior_mixture(list(
      prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      prior_factor("mnormal", list(0, 2), contrast = "orthonormal")
    ))
  )
  for(wrapper in names(wrapped_priors)){
    transformed <- build_fit(
      wrapped_priors[[wrapper]],
      matrix(c(1, 10, 2, 20), ncol = 2L, byrow = TRUE)
    )
    catalog <- parameter_catalog(transformed$fit)
    selections <- lapply(levels(data$f), function(level){
      parameter_catalog_resolve(
        catalog,
        paste0("f[", level, "]"),
        namespace = "mu"
      )
    })
    expect_true(all(vapply(selections, function(selection){
      identical(selection$quantities$status, "derived")
    }, logical(1))), info = wrapper)
  }

  ordered <- build_fit(
    prior_ordered(prior("normal", list(0, 1))),
    matrix(c(1, 10, 2, 20), ncol = 2L, byrow = TRUE)
  )
  ordered_catalog <- parameter_catalog(ordered$fit)
  expect_identical(
    parameter_catalog_resolve(
      ordered_catalog,
      "f[a]",
      namespace = "mu"
    )$quantities$status,
    "structural"
  )
  expect_identical(
    parameter_catalog_resolve(
      ordered_catalog,
      "f[b]",
      namespace = "mu"
    )$quantities$canonical_name,
    "mu_f[1]"
  )
  ordered_c <- parameter_catalog_resolve(
    ordered_catalog,
    "f[c]",
    namespace = "mu"
  )
  expect_identical(ordered_c$quantities$status, "derived")
  expect_identical(
    as.numeric(parameter_draws(ordered$fit, ordered_c)[[1L]][, 1L]),
    c(11, 22)
  )

  ordered_levels <- build_fit(
    prior_ordered(
      prior("normal", list(0, 1)),
      contrast = "cumulative_levels"
    ),
    matrix(c(1, 10, 100, 2, 20, 200), ncol = 3L, byrow = TRUE)
  )
  ordered_levels_catalog <- parameter_catalog(ordered_levels$fit)
  ordered_levels_b <- parameter_catalog_resolve(
    ordered_levels_catalog,
    "f[b]",
    namespace = "mu"
  )
  ordered_levels_c <- parameter_catalog_resolve(
    ordered_levels_catalog,
    "f[c]",
    namespace = "mu"
  )
  expect_identical(
    as.numeric(
      parameter_draws(ordered_levels$fit, ordered_levels_b)[[1L]][, 1L]
    ),
    c(11, 22)
  )
  expect_identical(
    as.numeric(
      parameter_draws(ordered_levels$fit, ordered_levels_c)[[1L]][, 1L]
    ),
    c(111, 222)
  )
})

test_that("factor interaction cells use only their persisted term design", {

  data <- data.frame(
    f = factor(c("a", "b", "c", "a", "b", "c")),
    g = factor(c("u", "v", "u", "v", "u", "v"))
  )
  formula_result <- JAGS_formula(
    ~ f * g,
    "mu",
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  coordinates <- unlist(lapply(names(formula_result$prior_list), function(name){
    prior <- formula_result$prior_list[[name]]
    if(is.prior.factor(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  registry <- .bt_build_parameter_registry(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    registry,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  )

  structural <- parameter_catalog_resolve(
    catalog,
    "f:g",
    namespace = "mu",
    component = "f=a, g=u"
  )$quantities
  expect_identical(structural$status, "structural")
  expect_identical(structural$fixed_value, 0)

  derived <- parameter_catalog_resolve(
    catalog,
    "f:g",
    namespace = "mu",
    component = "f=b, g=v"
  )$quantities
  expect_identical(derived$status, "derived")
  expect_setequal(
    derived$extraction_key[[1L]]$dependencies,
    c("mu_f__xXx__g[1]", "mu_f__xXx__g[2]")
  )
  expect_false(any(c("mu_intercept", "mu_f[1]", "mu_g") %in%
                     derived$extraction_key[[1L]]$dependencies))

  resolved <- hypothesis_resolve(
    hypothesis_parse("`f:g[f=b, g=v]` > 0"),
    catalog,
    namespace = "mu"
  )
  expect_identical(unique(resolved$occurrences$component), "f=b, g=v")
})

test_that("factor interaction components quote ambiguous level delimiters", {

  data <- expand.grid(
    f = c("a", "a, g=v"),
    g = c("u", "v, g=u"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = TRUE
  )
  formula_result <- JAGS_formula(
    ~ f * g,
    "mu",
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment"),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_factor(
        "mnormal",
        list(0, 1),
        contrast = "orthonormal"
      )
    )
  )
  coordinates <- unlist(lapply(names(formula_result$prior_list), function(name){
    prior <- formula_result$prior_list[[name]]
    if(is.prior.factor(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  registry <- .bt_build_parameter_registry(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    registry,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  )

  first <- parameter_catalog_resolve(
    catalog,
    "f:g",
    namespace = "mu",
    component = "f=\"a, g=v\", g=u"
  )
  second <- parameter_catalog_resolve(
    catalog,
    "f:g",
    namespace = "mu",
    component = "f=a, g=\"v, g=u\""
  )
  expect_false(identical(first$quantity_id, second$quantity_id))
  expect_identical(
    first$quantities$component,
    "f=\"a, g=v\", g=u"
  )
  expect_identical(
    second$quantities$component,
    "f=a, g=\"v, g=u\""
  )
  resolved <- hypothesis_resolve(
    hypothesis_parse("`f:g[f=\"a, g=v\", g=u]` > 0"),
    catalog,
    namespace = "mu"
  )
  expect_identical(
    unique(resolved$occurrences$quantity_id),
    first$quantity_id
  )
})

test_that("factor components escape hypothesis syntax characters", {

  data <- data.frame(f = factor(c("", "a]b")))
  formula_result <- JAGS_formula(
    ~ 1 + f,
    "mu",
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  prior <- formula_result$prior_list$mu_f
  coordinates <- c(
    "mu_intercept",
    .JAGS_prior_factor_names("mu_f", prior)
  )
  registry <- .bt_build_parameter_registry(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    registry,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  )

  empty <- parameter_catalog_resolve(
    catalog,
    "f[\"\"]",
    namespace = "mu"
  )
  bracket <- parameter_catalog_resolve(
    catalog,
    "f[a%5Db]",
    namespace = "mu"
  )
  expect_identical(empty$quantities$component, "\"\"")
  expect_identical(bracket$quantities$component, "a%5Db")
  resolved <- hypothesis_resolve(
    hypothesis_parse("`f[a%5Db]` > `f[\"\"]`"),
    catalog,
    namespace = "mu"
  )
  expect_setequal(
    unique(resolved$occurrences$quantity_id),
    c(empty$quantity_id, bracket$quantity_id)
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
  broken$schema_version <- 3L
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
