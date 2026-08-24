skip_if_not_test_profile("unit")

test_that("parameter catalog construction is metadata-only and versioned", {

  prior_list <- list(
    theta = prior("normal", list(0, 1)),
    fixed = prior("point", list(-3))
  )
  coordinates <- .bt_build_parameter_coordinates(
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
    coordinates = coordinates,
    prior_list = prior_list
  ))
  expect_s3_class(catalog, "BayesTools_parameter_catalog")
  expect_identical(catalog$schema_version, 3L)
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
  coordinates <- .bt_build_parameter_coordinates(
    columns = c("mu_intercept", "mu_f[1]", "mu_f[2]"),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )

  factor_rows <- catalog$quantities[
    catalog$quantities$canonical_name %in% c("mu_f[1]", "mu_f[2]"),
    ,
    drop = FALSE
  ]
  expect_identical(factor_rows$component, c("b", "c"))
  expect_identical(factor_rows$display_label, c("(mu) f[b]", "(mu) f[c]"))
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
  expect_identical(reference$display_label, "(mu) f[a]")
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

  incomplete <- coordinates[coordinates$coordinate_name != "mu_f[2]", , drop = FALSE]
  expect_error(
    .bt_build_parameter_catalog(
      coordinates = incomplete,
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design)
    ),
    "factor coordinates are missing or malformed"
  )
})

test_that("factor cell mapping is limited to fixed formula terms", {

  ordinary <- prior_factor(
    "mnormal",
    list(0, 1),
    contrast = "orthonormal"
  )
  attr(ordinary, "levels") <- 3L
  ordinary_coordinates <- .bt_build_parameter_coordinates(
    .JAGS_prior_factor_names("p1", ordinary),
    prior_list = list(p1 = ordinary)
  )
  ordinary_catalog <- .bt_build_parameter_catalog(
    ordinary_coordinates,
    prior_list = list(p1 = ordinary)
  )
  expect_false(any(vapply(
    ordinary_catalog$quantities$extraction_key,
    function(key) identical(key$type, "factor_level"),
    logical(1)
  )))

  data <- data.frame(
    f = factor(c("a", "b", "c")),
    id = factor(c("g1", "g2", "g1"))
  )
  formula_result <- JAGS_formula(
    ~ 1 + (f || id),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, 1)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  random_coordinates <- .bt_build_parameter_coordinates(
    c("mu_intercept", random_term$sd_parameter_names),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  expect_silent(.bt_build_parameter_catalog(
    random_coordinates,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  ))
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
  coordinates <- .bt_build_parameter_coordinates(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
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
    hypothesis_parse("f:g[f=b, g=v] > 0"),
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
  coordinates <- .bt_build_parameter_coordinates(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
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
  coordinates <- .bt_build_parameter_coordinates(
    columns = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
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

test_that("factor components preserve boundary whitespace", {

  data <- data.frame(f = factor(c(" a", "a ")))
  formula_result <- JAGS_formula(
    ~ 1 + f,
    "mu",
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "independent")
    )
  )
  prior <- formula_result$prior_list$mu_f
  coordinates <- .bt_build_parameter_coordinates(
    columns = c("mu_intercept", .JAGS_prior_factor_names("mu_f", prior)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  )

  leading <- parameter_catalog_resolve(catalog, "f[\" a\"]", "mu")
  trailing <- parameter_catalog_resolve(catalog, "f[\"a \"]", "mu")
  expect_identical(leading$quantities$component, "\" a\"")
  expect_identical(trailing$quantities$component, "\"a \"")
  resolved <- hypothesis_resolve(
    hypothesis_parse("\`f[\" a\"]\` > \`f[\"a \"]\`"),
    catalog,
    namespace = "mu"
  )
  expect_setequal(
    unique(resolved$occurrences$quantity_id),
    c(leading$quantity_id, trailing$quantity_id)
  )
})

test_that("catalog extensions preserve ambiguity until filtered", {

  coordinates <- .bt_build_parameter_coordinates(
    columns = "theta",
    prior_list = list(theta = prior("normal", list(0, 1)))
  )
  catalog <- .bt_build_parameter_catalog(coordinates)
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
    simplified = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  expect_error(
    parameter_catalog_extend(
      catalog,
      quantities = location,
      aliases = .bt_parameter_catalog_empty_aliases(),
      provider = "BayesTools"
    ),
    "reserved"
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
      simplified = FALSE,
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
  derived <- catalog$quantities[
    catalog$quantities$role == "random_correlation",
  ]
  expect_identical(unique(derived$status), "sampled")
  expect_false(any(catalog$quantities$internal))
  correlation <- derived[1L, , drop = FALSE]
  expect_identical(correlation$source_type, "one_to_one_transform")
  expect_identical(
    correlation$extraction_key[[1L]]$source_parameter,
    random_term$correlation$primitive_names
  )
  expect_identical(
    correlation$extraction_key[[1L]]$source_transform,
    "lkj2"
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "(mu) cor(intercept,x)",
      namespace = "mu"
    )$quantity_id,
    correlation$quantity_id
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "cor(intercept,x)",
      namespace = "mu"
    )$quantity_id,
    correlation$quantity_id
  )
  expect_setequal(
    correlation$extraction_key[[1L]]$dependencies,
    c(
      paste0(cholesky, "[1,1]"),
      paste0(cholesky, "[2,1]"),
      paste0(cholesky, "[2,2]")
    )
  )
  expect_false(any(c("mu_intercept", random_term$sd_parameter_names) %in%
                     correlation$extraction_key[[1L]]$dependencies))

  sd_label <- "(mu) sd(intercept)"
  expect_identical(
    sum(catalog$quantities$display_label == sd_label),
    1L
  )
  sd_quantity <- parameter_catalog_resolve(
    catalog,
    sd_label,
    namespace = "mu"
  )$quantities
  expect_identical(sd_quantity$status, "sampled")
  sd_draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(catalog, sd_label, namespace = "mu")
  )
  expect_identical(as.numeric(sd_draws[[1L]][, 1L]), c(1, 2, 3))
  var_label <- "(mu) var(intercept)"
  var_quantity <- parameter_catalog_resolve(
    catalog,
    var_label,
    namespace = "mu"
  )$quantities
  expect_identical(var_quantity$role, "random_var")
  var_draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(catalog, var_label, namespace = "mu")
  )
  expect_identical(as.numeric(var_draws[[1L]][, 1L]), c(1, 4, 9))

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
    correlation$canonical_name
  )
  transform <- parameter_transform(fit, selection)
  expect_identical(
    transform,
    list(type = "affine", offset = -1, scale = 2)
  )
  expect_equal(
    parameter_transform_forward(c(.25, .75), transform),
    c(-.5, .5)
  )
  draws <- parameter_draws(fit, selection)

  expect_identical(observed, correlation$extraction_key[[1L]]$dependencies)
  expect_identical(colnames(draws[[1L]]), correlation$canonical_name)
  expect_equal(as.numeric(draws[[1L]][, 1L]), c(0, .2, .4))
})

test_that("semantic parameter transforms own scalar transform algebra", {

  transforms <- list(
    identity = list(type = "identity"),
    affine = list(type = "affine", offset = -1, scale = 2),
    tanh = list(type = "tanh"),
    bounded_logit = list(type = "bounded_logit", lower = -.5, upper = 1),
    sqrt_scale = list(type = "sqrt_scale", scale = 4),
    square = list(type = "square")
  )
  source <- list(
    identity = c(-1, 1),
    affine = c(.25, .75),
    tanh = c(-.5, .5),
    bounded_logit = c(-1, 1),
    sqrt_scale = c(.25, 1),
    square = c(1, 2)
  )
  for(name in names(transforms)){
    transformed <- parameter_transform_forward(source[[name]], transforms[[name]])
    expect_equal(
      parameter_transform_inverse(transformed, transforms[[name]]),
      source[[name]],
      info = name
    )
    expect_true(
      all(parameter_transform_jacobian(source[[name]], transforms[[name]]) > 0),
      info = name
    )
  }

  expect_error(
    parameter_transform_forward(1, list(type = "affine", offset = 0, scale = 0)),
    "Unsupported semantic parameter transform",
    fixed = TRUE
  )

  sqrt_transform <- list(type = "sqrt_scale", scale = 4)
  expect_warning(
    expect_equal(
      parameter_transform_forward(c(-1, 0, 1), sqrt_transform),
      c(NaN, 0, 2)
    ),
    NA
  )
  expect_warning(
    expect_equal(
      parameter_transform_jacobian(c(-1, 0, 1), sqrt_transform),
      c(NaN, Inf, 1)
    ),
    NA
  )
})

test_that("prior sampling includes stored LKJ primitive coordinates", {

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
        cor = prior_lkj(eta = 2)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  primitive_name <- random_term$correlation$primitive_names
  columns <- c(
    "mu_intercept",
    random_term$sd_parameter_names,
    primitive_name
  )
  fit <- coda::mcmc.list(coda::mcmc(matrix(
    0,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  samples <- transform_prior_samples(fit, n_samples = 5000L, seed = 914L)

  expect_true(primitive_name %in% colnames(samples))
  expect_true(all(samples[, primitive_name] > 0))
  expect_true(all(samples[, primitive_name] < 1))
  expect_equal(mean(samples[, primitive_name]), 0.5, tolerance = 0.02)
  expect_lt(abs(stats::var(samples[, primitive_name]) - 0.05), 0.005)
})

test_that("structured correlation aliases expose shared and pairwise semantics", {

  data <- data.frame(
    outcome = factor(c("sensitivity", "specificity", "sensitivity")),
    study = factor(c("a", "a", "b"))
  )
  formula_result <- JAGS_formula(
    ~ 1 + hcs(outcome | study),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(
        sd = prior("gamma", list(2, 2)),
        covariance = random_covariance(
          cor = prior("uniform", list(-1, 1)),
          cor_scale = "cor"
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  rho_name <- random_term$correlation$rho_name
  values <- cbind(
    mu_intercept = 0,
    matrix(1, nrow = 2L, ncol = 2L),
    c(.1, .2)
  )
  colnames(values)[2:3] <- random_term$sd_parameter_names
  colnames(values)[4L] <- rho_name
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- parameter_catalog(fit)
  rho <- parameter_catalog_resolve(catalog, "cor", "mu")

  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "(mu) cor(outcome[sensitivity],outcome[specificity])",
      "mu"
    )$quantity_id,
    rho$quantity_id
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "cor(outcome[sensitivity],outcome[specificity])",
      "mu"
    )$quantity_id,
    rho$quantity_id
  )
  expect_error(
    parameter_catalog_resolve(catalog, "study: cor", "mu"),
    "No public parameter quantity matches"
  )
})

test_that("explicitly named one-entry random lists retain their public owner", {

  data <- data.frame(id = factor(c("a", "a", "b", "b")))
  formula_result <- JAGS_formula(
    random_effects_formula(list(study = ~ 1 | id)),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  coordinates <- .bt_build_parameter_coordinates(
    columns = random_term$sd_parameter_names,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )

  expect_setequal(
    catalog$quantities$canonical_name,
    c("(mu) study: sd(intercept)", "(mu) study: var(intercept)")
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "study: sd(intercept)",
      "mu"
    )$quantities$canonical_name,
    "(mu) study: sd(intercept)"
  )
  expect_error(
    parameter_catalog_resolve(catalog, "sd(intercept)", "mu"),
    "No public parameter quantity matches"
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "study: sd",
      "mu",
      simplify_names = TRUE
    )$quantities$canonical_name,
    "(mu) study: sd(intercept)"
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "sd",
      "mu",
      simplify_names = TRUE
    )$quantities$canonical_name,
    "(mu) study: sd(intercept)"
  )
  expect_identical(
    catalog$quantities$display_label[
      catalog$quantities$quantity == "sd"
    ],
    "(mu) study: sd"
  )
  simplified_ast <- hypothesis_parse(
    "study: sd > 0",
    catalog = catalog,
    namespace = "mu",
    simplify_names = TRUE
  )
  expect_identical(
    unique(hypothesis_resolve(
      simplified_ast,
      catalog,
      namespace = "mu",
      simplify_names = TRUE
    )$occurrences$canonical_name),
    "(mu) study: sd(intercept)"
  )
})

test_that("random covariance families share one semantic naming grammar", {

  data <- data.frame(
    group = factor(c("a", "a", "b", "b")),
    level = factor(c("x", "y", "x", "y")),
    time = 1:4
  )
  prior_list <- list(intercept = prior("normal", list(0, 1)))
  random_prior <- prior_random(sd = prior("gamma", list(2, 2)))
  formulas <- list(
    id = ~ 1 + id(1 | group),
    diag = ~ 1 + diag(1 + time | group),
    us = ~ 1 + us(1 + time | group),
    cs = ~ 1 + cs(level | group),
    hcs = ~ 1 + hcs(level | group),
    ar1 = ~ 1 + ar1(level | group),
    har = ~ 1 + har(level | group),
    car = ~ 1 + car(time | group)
  )
  expected <- list(
    id = "(mu) sd",
    diag = c(
      "(mu) sd(intercept)",
      "(mu) sd(time)"
    ),
    us = c(
      "(mu) sd(intercept)",
      "(mu) sd(time)",
      "(mu) cor(intercept,time)"
    ),
    cs = c("(mu) sd", "(mu) cor"),
    hcs = c(
      "(mu) sd(level[x])",
      "(mu) sd(level[y])",
      "(mu) cor"
    ),
    ar1 = c("(mu) sd", "(mu) cor"),
    har = c(
      "(mu) sd(level[x])",
      "(mu) sd(level[y])",
      "(mu) cor"
    ),
    car = c("(mu) sd", "(mu) cor")
  )

  catalogs <- lapply(formulas, function(formula){
    formula_result <- JAGS_formula(
      formula = formula,
      parameter = "mu",
      data = data,
      prior_list = prior_list,
      prior_random = random_prior
    )
    random_term <- formula_result$formula_design$random_effects[[1L]]
    correlation <- random_term$correlation
    correlation_coordinates <- character()
    if(!is.null(correlation) && identical(correlation$type, "rho")){
      correlation_coordinates <- unique(c(
        correlation$rho_name,
        correlation$sample_name
      ))
    }
    if(!is.null(correlation) && identical(correlation$type, "lkj")){
      indices <- which(
        lower.tri(matrix(0, random_term$n_columns, random_term$n_columns),
                  diag = TRUE),
        arr.ind = TRUE
      )
      correlation_coordinates <- paste0(
        correlation$cholesky_name,
        "[", indices[, "row"], ",", indices[, "col"], "]"
      )
    }
    columns <- unique(c(
      "mu_intercept",
      random_term$sd_parameter_names,
      correlation_coordinates
    ))
    coordinates <- .bt_build_parameter_coordinates(
      columns = columns,
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design)
    )
    .bt_build_parameter_catalog(
      coordinates = coordinates,
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design)
    )
  })

  for(name in names(catalogs)){
    catalog <- catalogs[[name]]
    random <- catalog$quantities[
      catalog$quantities$role %in% c("random_sd", "random_correlation"),
      ,
      drop = FALSE
    ]
    expect_setequal(random$canonical_name, expected[[name]])
    expect_true(all(random$status == "sampled"))
    for(canonical_name in random$canonical_name){
      prefixless <- sub("^\\(mu\\) ", "", canonical_name)
      expect_identical(
        parameter_catalog_resolve(catalog, prefixless, "mu")$quantity_id,
        random$quantity_id[random$canonical_name == canonical_name]
      )
    }
  }

  expect_identical(
    parameter_catalog_resolve(
      catalogs$cs,
      "cor(level[x],level[y])",
      "mu"
    )$quantities$canonical_name,
    "(mu) cor"
  )
  expect_error(
    parameter_catalog_resolve(
      catalogs$ar1,
      "cor(level[x],level[y])",
      "mu"
    ),
    "No public parameter quantity matches"
  )
  expect_error(parameter_catalog_resolve(catalogs$us, "cor", "mu"),
               "No public parameter quantity matches")
})

test_that("transformed random summaries hide fitted-scale implementation rows", {

  data <- data.frame(
    x = 1:4,
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    ~ 1 + diag(0 + x | id),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  sd_name <- random_term$sd_parameter_names
  values <- cbind(mu_intercept = 0, c(2, 4))
  colnames(values)[2L] <- sd_name
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design),
    formula_scale = list(mu = formula_result$formula_scale)
  )
  catalog <- parameter_catalog(fit)

  label <- "(mu) sd(x)"
  quantity <- parameter_catalog_resolve(catalog, label, "mu")$quantities
  expect_identical(quantity$status, "sampled")
  expect_identical(quantity$source_type, "one_to_one_transform")
  expect_equal(
    quantity$extraction_key[[1L]]$source_scale,
    1 / stats::sd(data$x)
  )
  expect_identical(quantity$extraction_key[[1L]]$dependencies, sd_name)
  expect_false(sd_name %in% catalog$quantities$canonical_name)
  draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(catalog, label, "mu")
  )
  expect_equal(
    as.numeric(draws[[1L]][, 1L]),
    c(2, 4) / stats::sd(data$x)
  )
})

test_that("transformed correlations declare SD and Cholesky inputs", {

  data <- data.frame(
    x = 1:4,
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    ~ 1 + us(1 + x | id),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj()
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  cholesky <- random_term$correlation$cholesky_name
  sd_values <- cbind(c(1, 2), c(2, 3))
  rho <- c(.2, .4)
  values <- cbind(
    mu_intercept = 0,
    sd_values,
    1,
    rho,
    sqrt(1 - rho^2)
  )
  colnames(values) <- c(
    "mu_intercept",
    random_term$sd_parameter_names,
    paste0(cholesky, "[1,1]"),
    paste0(cholesky, "[2,1]"),
    paste0(cholesky, "[2,2]")
  )
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design),
    formula_scale = list(mu = formula_result$formula_scale)
  )
  catalog <- parameter_catalog(fit)
  correlation <- catalog$quantities[
    catalog$quantities$role == "random_correlation" &
      catalog$quantities$status == "sampled",
    ,
    drop = FALSE
  ]
  expect_identical(nrow(correlation), 1L)
  expect_setequal(
    correlation$extraction_key[[1L]]$dependencies,
    colnames(values)[-1L]
  )

  mean_x <- mean(data$x)
  sd_x <- stats::sd(data$x)
  fitted_covariance <- rho * sd_values[, 1L] * sd_values[, 2L]
  intercept_variance <- sd_values[, 1L]^2 +
    (mean_x / sd_x)^2 * sd_values[, 2L]^2 -
    2 * (mean_x / sd_x) * fitted_covariance
  slope_variance <- sd_values[, 2L]^2 / sd_x^2
  covariance <- fitted_covariance / sd_x -
    mean_x * sd_values[, 2L]^2 / sd_x^2
  expected <- covariance / sqrt(intercept_variance * slope_variance)
  draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(catalog, correlation$canonical_name)
  )
  expect_equal(as.numeric(draws[[1L]][, 1L]), expected)
})

test_that("identity random summaries preserve structural provenance", {

  data <- data.frame(id = factor(c("a", "b")))
  formula_result <- JAGS_formula(
    ~ 1 + (1 | id),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("point", list(location = 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  coordinates <- .bt_build_parameter_coordinates(
    columns = c("mu_intercept", random_term$sd_parameter_names),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
    formula_result$prior_list,
    list(mu = formula_result$formula_design)
  )

  label <- "(mu) sd(intercept)"
  quantity <- parameter_catalog_resolve(
    catalog,
    label,
    namespace = "mu"
  )$quantities
  expect_identical(quantity$status, "structural")
  expect_identical(quantity$fixed_value, 2)
  expect_identical(
    catalog$quantities$display_label[
      catalog$quantities$canonical_name == label
    ],
    "(mu) sd"
  )
  expect_false(any(
    catalog$quantities$role == "random_sd" &
      catalog$quantities$status == "sampled"
  ))
  expect_identical(
    parameter_catalog_resolve(catalog, "intercept", "mu")$quantities$role,
    "fixed_coefficient"
  )
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
      allocation = random_variance_allocation(name = "allocation",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  random_terms <- formula_result$formula_design$random_effects
  allocation_columns <- c(
    "mu__xRE_ALLOCx_allocation__allocation_sd",
    "mu__xRE_ALLOCx_allocation__weight[1]",
    "mu__xRE_ALLOCx_allocation__weight[2]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
  )
  columns <- c(
    "mu_intercept",
    allocation_columns,
    unlist(lapply(random_terms, `[[`, "sd_parameter_names"))
  )
  coordinates <- .bt_build_parameter_coordinates(
    columns = columns,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  derived <- catalog$quantities[
    startsWith(catalog$quantities$role, "random_"),
  ]

  expect_true(all(coordinates$internal[
    coordinates$coordinate_name %in% allocation_columns
  ]))
  expect_false(any(
    allocation_columns %in% catalog$quantities$canonical_name
  ))

  expect_true("random_var_prop" %in% derived$role)
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      "allocation: sd_total",
      "mu"
    )$quantities$role,
    "random_sd_total"
  )
  allocation_var <- parameter_catalog_resolve(
    catalog,
    "allocation: var_total",
    "mu"
  )$quantities
  expect_identical(allocation_var$role, "random_var_total")
  expect_identical(allocation_var$scale_role, "total")
  expect_identical(allocation_var$source_type, "one_to_one_transform")
  expect_identical(
    allocation_var$extraction_key[[1L]]$source_transform,
    "square"
  )
  allocation_sd <- parameter_catalog_resolve(
    catalog,
    "allocation: sd_total",
    "mu"
  )$quantities
  expect_identical(allocation_sd$scale_role, "total")
  expect_identical(allocation_var$parent_quantity_id, allocation_sd$quantity_id)
  study_sd_label <- "(mu) study: sd(intercept)"
  study_sd <- parameter_catalog_resolve(
    catalog,
    alias = study_sd_label,
    namespace = "mu"
  )
  expect_identical(study_sd$quantities$status, "sampled")
  expect_identical(
    catalog$quantities$display_label[
      catalog$quantities$canonical_name == study_sd_label
    ],
    "(mu) study: sd"
  )
  expect_error(
    parameter_catalog_resolve(
      catalog,
      alias = "sd",
      namespace = "mu",
      simplify_names = TRUE
    ),
    class = "BayesTools_parameter_ambiguous"
  )
  fractions <- derived[derived$role == "random_var_prop", ]
  expect_identical(fractions$component, c("study", "drug"))
  expect_true(all(vapply(fractions$extraction_key, function(key){
    identical(
      key$dependencies,
      c(
        "mu__xRE_ALLOCx_allocation__weight[1]",
        "mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  }, logical(1))))
  study_fraction <- parameter_catalog_resolve(
    catalog,
    alias = "allocation: var_prop(study)",
    namespace = "mu"
  )
  expect_identical(
    study_fraction$quantities$canonical_name,
    fractions$canonical_name[fractions$component == "study"]
  )

  unprefixed_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name            = "internal_allocation",
        display_name    = "",
        terms           = c(component_1 = "study", component_2 = "drug"),
        component_names = c("component 1", "component 2"),
        sd              = prior("gamma", list(2, 2))
      )
    )
  )
  unprefixed_coordinates <- .bt_build_parameter_coordinates(
    columns = c(
      "mu_intercept",
      "mu__xRE_ALLOCx_internal_allocation__allocation_sd",
      "mu__xRE_ALLOCx_internal_allocation__weight[1]",
      "mu__xRE_ALLOCx_internal_allocation__weight[2]"
    ),
    prior_list = unprefixed_result$prior_list,
    formula_design = list(mu = unprefixed_result$formula_design)
  )
  unprefixed_catalog <- .bt_build_parameter_catalog(
    coordinates = unprefixed_coordinates,
    prior_list = unprefixed_result$prior_list,
    formula_design = list(mu = unprefixed_result$formula_design)
  )
  expect_identical(
    parameter_catalog_resolve(
      unprefixed_catalog,
      "sd_total",
      "mu"
    )$quantities$canonical_name,
    "(mu) sd_total"
  )
  expect_identical(
    parameter_catalog_resolve(
      unprefixed_catalog,
      "var_prop(component 1)",
      "mu"
    )$quantities$component,
    "component 1"
  )
  unprefixed_hypothesis <- hypothesis_parse(
    "var_prop(component 1) = 0",
    catalog = unprefixed_catalog
  )
  expect_identical(
    unique(hypothesis_resolve(
      unprefixed_hypothesis,
      unprefixed_catalog,
      namespace = "mu"
    )$occurrences$quantity_id),
    parameter_catalog_resolve(
      unprefixed_catalog,
      "var_prop(component 1)",
      "mu"
    )$quantity_id
  )
})

test_that("allocation quantities expose deterministic induced prior densities", {

  data <- data.frame(
    group = factor(c("a", "b", "a", "b")),
    study = factor(c("s1", "s1", "s2", "s2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 0 + group + us(0 + group | study),
    parameter = "mu",
    data = data,
    prior_list = list(
      group = prior_factor(
        "normal",
        list(mean = 0, sd = 1),
        contrast = "independent"
      )
    ),
    prior_random = prior_random(
      study = random_block(contrasts = c(group = "independent")),
      random_variance_allocation(
        name = "heterogeneity",
        display_name = "",
        terms = "study",
        target = "sd_component",
        scale = "mean_variance",
        sd = prior(
          "normal",
          list(mean = 0, sd = 1),
          truncation = list(0, Inf)
        ),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  allocation  <- random_term$sd_binding$allocations[[1L]]
  cholesky    <- random_term$correlation$cholesky_name
  columns <- c(
    "mu_group[1]", "mu_group[2]",
    allocation$source_node,
    paste0(allocation$weight_name, "[", 1:2, "]"),
    random_term$correlation$primitive_names,
    paste0(cholesky, c("[1,1]", "[2,1]", "[1,2]", "[2,2]"))
  )
  values <- rbind(
    c(0, 0, 1, .4, .6, .25, 1, -.5, 0, sqrt(.75)),
    c(0, 0, 1, .4, .6, .25, 1, -.5, 0, sqrt(.75))
  )
  colnames(values) <- columns
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- parameter_catalog(fit)
  component_selection <- parameter_catalog_resolve(
    catalog,
    "sd(group[a])",
    "mu"
  )
  multiplier_selection <- parameter_catalog_resolve(
    catalog,
    "var_mult(group[a])",
    "mu"
  )
  sd_mult_selection <- parameter_catalog_resolve(
    catalog,
    "sd_mult(group[a])",
    "mu"
  )
  expect_error(
    parameter_catalog_resolve(
      catalog,
      paste0("var_", "ratio(group[a])"),
      "mu"
    ),
    "No public parameter quantity matches",
    fixed = TRUE
  )
  expect_error(
    parameter_catalog_resolve(
      catalog,
      paste0("sd_", "ratio(group[a])"),
      "mu"
    ),
    "No public parameter quantity matches",
    fixed = TRUE
  )
  common_variance_selection <- parameter_catalog_resolve(
    catalog,
    "var_common",
    "mu"
  )
  component_density <- parameter_prior_density(
    fit,
    component_selection,
    n_grid = 512L
  )
  repeated_density <- parameter_prior_density(
    fit,
    component_selection,
    n_grid = 512L
  )
  multiplier_density <- parameter_prior_density(
    fit,
    multiplier_selection,
    n_grid = 512L
  )
  sd_mult_density <- parameter_prior_density(
    fit,
    sd_mult_selection,
    n_grid = 512L
  )
  common_variance_density <- parameter_prior_density(
    fit,
    common_variance_selection,
    n_grid = 512L
  )
  x <- component_density$density$x
  y <- component_density$density$y
  second_moment <- sum(diff(x) *
    (head(x^2 * y, -1L) + tail(x^2 * y, -1L)) / 2)

  expect_s3_class(component_density, "prior_linear_density")
  expect_identical(component_density$density, repeated_density$density)
  expect_equal(second_moment, 2 * (2 / 5), tolerance = .03)
  expect_equal(
    .prior_linear_density_height(multiplier_density, .5),
    stats::dbeta(.25, 2, 3) / 2,
    tolerance = .02
  )
  expect_equal(
    .prior_linear_density_height(sd_mult_density, 1),
    stats::dbeta(.5, 2, 3),
    tolerance = .02
  )
  common_variance_interior <- prior_density_ordinate(
    common_variance_density,
    .5
  )
  expect_identical(common_variance_interior$behavior, "regular")
  expect_true(common_variance_interior$exact)
  expect_equal(
    common_variance_interior$log_density,
    stats::dchisq(.5, df = 1, log = TRUE)
  )
  common_variance_boundary <- prior_density_ordinate(
    common_variance_density,
    0
  )
  expect_identical(common_variance_boundary$behavior, "infinite")
  expect_true(common_variance_boundary$exact)
  standard <- .bt_parameter_catalog_random_summary_samples(
    fit = fit,
    model_samples = values,
    prior_list = formula_result$prior_list,
    coordinates = parameter_coordinates(fit),
    mode = "standard"
  )$model_samples
  semantic <- intersect(
    colnames(standard),
    catalog$quantities$canonical_name[startsWith(
      catalog$quantities$role,
      "random_"
    )]
  )
  expect_identical(
    sub("^\\(mu\\) ", "", semantic),
    c(
      "sd_common",
      "var_mult(group[a])",
      "var_mult(group[b])",
      "cor(group[a],group[b])"
    )
  )
})


test_that("block allocations expose deterministic component-SD priors", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    observation = factor(seq_len(4L))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | observation, name = "observation", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      random_variance_allocation(
        name            = "heterogeneity",
        display_name    = "",
        terms           = c(study = "study", observation = "observation"),
        component_names = c("study", "observation"),
        sd = prior(
          "normal",
          list(mean = 0, sd = 1),
          truncation = list(0, Inf)
        ),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  allocation <- formula_result$formula_design$random_allocations[[1L]]
  columns <- c(
    "mu_intercept",
    allocation$source_node,
    paste0(allocation$weight_name, "[", 1:2, "]")
  )
  values <- matrix(
    rep(c(0, 1, .4, .6), 2L),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, columns)
  )
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  selection <- parameter_catalog_resolve(
    parameter_catalog(fit),
    "study: sd(intercept)",
    "mu"
  )
  density <- parameter_prior_density(fit, selection, n_grid = 512L)
  x <- density$density$x
  y <- density$density$y
  second_moment <- sum(diff(x) *
    (head(x^2 * y, -1L) + tail(x^2 * y, -1L)) / 2)

  expect_s3_class(density, "prior_linear_density")
  expect_equal(second_moment, 2 / 5, tolerance = .03)
})


test_that("unnamed local allocations retain owners with multiple blocks", {

  data <- data.frame(
    study = factor(c("a", "a", "b", "b")),
    site  = factor(c("x", "y", "x", "y")),
    z     = c(-1, 0, 1, 2)
  )
  scale_prior  <- prior("gamma", list(2, 2))
  weight_prior <- prior("dirichlet", list(alpha = c(1, 1)))
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + z | study, name = "study", covariance = "diag") +
      random(1 + z | site, name = "site", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      random_variance_allocation(
        name = "study_allocation", display_name = "", terms = "study",
        target = "sd_component", scale = "mean_variance",
        sd = scale_prior, weights = weight_prior
      ),
      random_variance_allocation(
        name = "site_allocation", display_name = "", terms = "site",
        target = "sd_component", scale = "mean_variance",
        sd = scale_prior, weights = weight_prior
      )
    )
  )
  allocations <- unlist(lapply(
    formula_result$formula_design$random_effects,
    function(term) term$sd_binding$allocations
  ), recursive = FALSE)
  columns <- unique(c(
    "mu_intercept",
    unlist(lapply(
      formula_result$formula_design$random_effects,
      `[[`,
      "sd_parameter_names"
    )),
    unlist(lapply(allocations, function(allocation){
      c(
        allocation$source_node,
        paste0(allocation$weight_name, "[", seq_len(allocation$n_targets), "]")
      )
    }))
  ))
  coordinates <- .bt_build_parameter_coordinates(
    columns = columns,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates = coordinates,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  common_sd <- catalog$quantities[
    catalog$quantities$quantity == "sd_common",
    ,
    drop = FALSE
  ]

  expect_setequal(
    common_sd$canonical_name,
    c("(mu) study: sd_common", "(mu) site: sd_common")
  )
  expect_error(
    parameter_catalog_resolve(catalog, "sd_common", "mu"),
    "No public parameter quantity matches"
  )
  expect_identical(
    parameter_catalog_resolve(catalog, "study: sd_common", "mu")$quantity_id,
    common_sd$quantity_id[common_sd$owner_name == "study"]
  )
})

test_that("malformed catalogs and stale selections fail closed", {

  coordinates <- .bt_build_parameter_coordinates(columns = "theta")
  catalog <- .bt_build_parameter_catalog(coordinates)
  selection <- parameter_catalog_resolve(catalog, "theta")

  broken <- catalog
  broken$schema_version <- 7L
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

  spoofed <- catalog
  spoofed$quantities$extraction_key[[1L]] <- list(
    type = "bogus",
    dependencies = "theta"
  )
  expect_error(
    .bt_validate_parameter_catalog(spoofed),
    "extraction keys are malformed"
  )
})
