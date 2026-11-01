skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Exported Schemas, Aliases, and Small Helpers
# ============================================================================ #
#
# PURPOSE:
#   Relationship tests for exports that no other test file references: schema
#   descriptions against the objects they describe, aliases against the
#   functions they alias, and small helpers against their documented behavior.
#
# DEPENDENCIES:
#   - helper-parameter-catalog.R: fit-free fitted object for the catalog schema
#
# SKIP CONDITIONS:
#   - None (no fits, no external data)
#
# TAGS: @evaluation, @schemas, @aliases
# ============================================================================ #


# Whether 'value' has the type that a schema row declares.
.schema_value_has_type <- function(value, type){

  switch(
    type,
    integer    = is.integer(value),
    numeric    = is.numeric(value),
    character  = is.character(value),
    logical    = is.logical(value),
    list       = is.list(value),
    data.frame = is.data.frame(value),
    stop("Unknown schema type '", type, "'.", call. = FALSE)
  )
}

# The fields of 'schema' whose element in 'object' does not have the declared type.
.schema_type_problems <- function(schema, object){

  held <- vapply(seq_len(nrow(schema)), function(i){
    .schema_value_has_type(object[[schema$field[i]]], schema$type[i])
  }, logical(1))

  schema$field[!held]
}

# Whether reading the fit contract of 'fit' stops with a refit request.
.contract_refused <- function(fit){

  inherits(tryCatch(JAGS_fit_contract(fit), error = identity), "BayesTools_refit_required")
}


# ---------------------------------------------------------------------------- #
# JAGS_draw_geometry_schema()
# ---------------------------------------------------------------------------- #

test_that("JAGS_draw_geometry_schema describes the fields of JAGS_draw_geometry()", {

  chains <- coda::mcmc.list(lapply(0:2, function(chain){
    coda::mcmc(
      matrix(chain * 8L + 1:8, nrow = 4, dimnames = list(NULL, c("theta", "private"))),
      start = 21, thin = 5
    )
  }))
  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "draw_geometry") <- .bt_draw_geometry_from_chains(chains)
  fit <- .bt_attach_fit_contract(fit)

  geometry <- JAGS_draw_geometry(fit)
  schema   <- JAGS_draw_geometry_schema()

  expect_s3_class(schema, "data.frame")
  expect_identical(names(schema), c("field", "type", "description"))
  expect_identical(anyDuplicated(schema$field), 0L)
  expect_true(all(nzchar(schema$description)))

  # the schema lists the fields of the geometry and of its per-chain table, each
  # group in the order of the object
  chain_fields <- names(geometry$chains)
  expect_setequal(schema$field, c(names(geometry), chain_fields))
  expect_identical(schema$field[schema$field %in% names(geometry)], names(geometry))
  expect_identical(schema$field[schema$field %in% chain_fields], chain_fields)

  # and declares the type of each one
  top_level <- schema$field %in% names(geometry)
  expect_identical(.schema_type_problems(schema[top_level, ], geometry), character())
  expect_identical(.schema_type_problems(schema[!top_level, ], geometry$chains), character())

  # the described values: one row per chain, timing as recorded by coda, and the
  # rows of each chain after chain-major concatenation of the draws
  mcpar <- t(vapply(chains, function(chain) attr(chain, "mcpar"), numeric(3L)))
  expect_identical(nrow(geometry$chains), length(chains))
  expect_equal(geometry$chains$start, mcpar[, 1L])
  expect_equal(geometry$chains$end,   mcpar[, 2L])
  expect_equal(geometry$chains$thin,  mcpar[, 3L])
  expect_equal(geometry$chains$iterations, vapply(chains, nrow, numeric(1L)))

  pooled <- do.call(rbind, lapply(chains, as.matrix))
  expect_equal(geometry$total_draws, nrow(pooled))
  expect_identical(geometry$chain_order, geometry$chains$chain)
  for(i in seq_along(chains)){
    rows <- geometry$chains$draw_start[i]:geometry$chains$draw_end[i]
    expect_equal(unname(pooled[rows, , drop = FALSE]), unname(as.matrix(chains[[i]])))
  }
})


# ---------------------------------------------------------------------------- #
# JAGS_fit_contract_schema()
# ---------------------------------------------------------------------------- #

test_that("JAGS_fit_contract_schema describes the contract that the validator enforces", {

  fit      <- .bt_attach_fit_contract(structure(list(), class = "BayesTools_fit"))
  contract <- JAGS_fit_contract(fit)
  schema   <- JAGS_fit_contract_schema()

  expect_s3_class(schema, "data.frame")
  expect_identical(names(schema), c("field", "type", "description"))
  expect_identical(anyDuplicated(schema$field), 0L)
  expect_true(all(nzchar(schema$description)))

  # the schema lists the fields of the contract, in its order, with their types
  expect_identical(schema$field, names(contract))
  expect_identical(schema$field[1L], "schema_version")
  expect_identical(schema$type, rep("integer", nrow(schema)))
  expect_identical(.schema_type_problems(schema, contract), character())

  # every other field is the version of a component that the validator accepts
  versions   <- setdiff(schema$field, "schema_version")
  components <- sub("_version$", "", versions)
  expect_identical(JAGS_validate_fit_contract(fit, requires = components), contract)
  expect_error(
    JAGS_validate_fit_contract(fit, requires = "not_a_component"),
    "Unknown fitted-contract component: 'not_a_component'.",
    fixed = TRUE
  )

  # a stale version of a component is refused for that component only
  for(i in seq_along(versions)){
    stale <- fit
    attr(stale, "fit_contract")[[versions[i]]] <- 0L
    refused <- tryCatch(
      JAGS_validate_fit_contract(stale, requires = components[i]),
      error = identity
    )
    expect_s3_class(refused, "BayesTools_refit_required")
    expect_identical(
      JAGS_validate_fit_contract(stale, requires = components[-i]),
      JAGS_fit_contract(stale)
    )
  }

  # and the contract object itself needs every described field as an integer
  incomplete <- vapply(schema$field, function(field){
    broken <- fit
    attr(broken, "fit_contract")[[field]] <- NULL
    .contract_refused(broken)
  }, logical(1))
  mistyped <- vapply(schema$field, function(field){
    broken <- fit
    attr(broken, "fit_contract")[[field]] <- as.numeric(contract[[field]])
    .contract_refused(broken)
  }, logical(1))
  expect_true(all(incomplete))
  expect_true(all(mistyped))
  expect_false(.contract_refused(fit))
})


# ---------------------------------------------------------------------------- #
# parameter_catalog_schema()
# ---------------------------------------------------------------------------- #

test_that("parameter_catalog_schema describes the tables of parameter_catalog()", {

  chains <- coda::mcmc.list(coda::mcmc(matrix(
    0.5, nrow = 2L, ncol = 2L, dimnames = list(NULL, c("theta", "tau"))
  )))
  sampled <- list(
    theta = prior("normal", list(0, 1)),
    tau   = prior("gamma", list(2, 1))
  )
  prior_lists <- list(
    sampled    = sampled,
    structural = c(sampled, list(fixed = prior("point", list(-3))))
  )

  schema <- parameter_catalog_schema()
  expect_identical(names(schema), c("schema_version", "quantities", "aliases"))
  for(table in c("quantities", "aliases")){
    expect_identical(names(schema[[table]]), c("field", "type", "description"), info = table)
    expect_identical(anyDuplicated(schema[[table]]$field), 0L, info = table)
    expect_true(all(nzchar(schema[[table]]$description)), info = table)
  }

  catalogs <- list()
  for(name in names(prior_lists)){
    fit     <- .parameter_catalog_test_fit(chains, prior_lists[[name]])
    catalog <- parameter_catalog(fit)
    catalogs[[name]] <- catalog

    # the versions and the columns of both tables
    expect_identical(schema$schema_version, catalog$schema_version, info = name)
    expect_identical(schema$schema_version, parameter_map(fit)$schema_version, info = name)
    expect_identical(schema$quantities$field, names(catalog$quantities), info = name)
    expect_identical(schema$aliases$field, names(catalog$aliases), info = name)

    # and their declared types, also for the columns that hold only NA or lists
    expect_identical(.schema_type_problems(schema$quantities, catalog$quantities), character(), info = name)
    expect_identical(.schema_type_problems(schema$aliases, catalog$aliases), character(), info = name)

    # the described values: a status from the declared set, and a fixed value
    # that is exact for a structural quantity and NA otherwise
    quantities <- catalog$quantities
    expect_true(all(quantities$status %in% c("sampled", "structural", "derived", "unavailable")), info = name)
    expect_identical(!is.na(quantities$fixed_value), quantities$status == "structural", info = name)
    expect_identical(anyDuplicated(quantities$quantity_id), 0L, info = name)
    expect_true(all(catalog$aliases$quantity_id %in% quantities$quantity_id), info = name)
  }

  structural <- catalogs$structural$quantities
  expect_identical(structural$fixed_value[structural$canonical_name == "fixed"], -3)
  expect_identical(structural$status[structural$canonical_name == "fixed"], "structural")
  expect_true(all(is.na(catalogs$sampled$quantities$fixed_value)))

  # the parameter map documents the same tables
  map_schema <- parameter_map_schema()
  expect_identical(map_schema$schema_version, schema$schema_version)
  expect_identical(map_schema$quantities, schema$quantities)
  expect_identical(map_schema$aliases, schema$aliases)
})


# ---------------------------------------------------------------------------- #
# JAGS_summary_table(), JAGS_estimates_empty_table(), JAGS_inference_empty_table()
# ---------------------------------------------------------------------------- #

# Two fit-free models that differ in the prior of 'm' (a normal prior and a
# spike at zero), with their model-averaging inference.
.alias_models_for_test <- function(){

  prior_lists <- list(
    list(
      m = prior("normal", list(0, 1)),
      s = prior("normal", list(0, 1), truncation = list(lower = 0))
    ),
    list(
      m = prior("point", list(0)),
      s = prior("normal", list(0, 1), truncation = list(lower = 0))
    )
  )
  models <- lapply(seq_along(prior_lists), function(i){
    fit <- list()
    attr(fit, "prior_list") <- prior_lists[[i]]
    list(
      fit           = fit,
      marglik       = .bt_marglik_manual_result(c(-1, -2)[i]),
      prior_weights = 1
    )
  })

  models_inference(models)
}

test_that("JAGS_summary_table gives the model summary table", {

  models <- .alias_models_for_test()
  argument_sets <- list(
    list(),
    list(remove_spike_0 = FALSE),
    list(short_name = TRUE, formula_prefix = FALSE),
    list(
      model_description = list(Type = "spike"),
      title = "Title", footnotes = "Note", warnings = "Warning"
    ),
    list(remove_parameters = "s")
  )

  expect_identical(formals(JAGS_summary_table), formals(model_summary_table))
  for(model in models){
    for(arguments in argument_sets){
      expect_identical(
        do.call(JAGS_summary_table, c(list(model), arguments)),
        do.call(model_summary_table, c(list(model), arguments))
      )
    }
  }

  # the summary reports the model and the priors that are not point zero
  table <- JAGS_summary_table(models[[2L]], title = "Title")
  expect_s3_class(table, "BayesTools_table")
  expect_identical(attr(table, "type"), c("string_left", "string", "string", "prior"))
  expect_identical(attr(table, "title"), "Title")
  expect_identical(
    trimws(unlist(table[, 1L])),
    c("Model", "Prior prob.", "log(marglik)", "Post. prob.", "Inclusion BF")
  )
  expect_identical(trimws(unlist(table[, 2L]))[1L], "2")
  expect_identical(
    trimws(unlist(table[, 4L]))[1:2],
    c("Parameter prior distributions", "s ~ Normal(0, 1)[0, Inf]")
  )
  expect_true(any(grepl("m ~ Spike(0)", unlist(JAGS_summary_table(models[[2L]], remove_spike_0 = FALSE)[, 4L]), fixed = TRUE)))
  expect_false(any(grepl("s ~", unlist(JAGS_summary_table(models[[2L]], remove_parameters = "s")[, 4L]), fixed = TRUE)))
})

test_that("JAGS_estimates_empty_table and JAGS_inference_empty_table give the empty runjags tables", {

  # the defaults of the diagnostic columns come from options
  withr::local_options(list(
    BayesTools.JAGS_estimates_diagnostic_columns = NULL,
    BayesTools.JAGS_BF_diagnostic_columns        = NULL
  ))

  expect_identical(formals(JAGS_estimates_empty_table), formals(runjags_estimates_empty_table))
  expect_identical(formals(JAGS_inference_empty_table), formals(runjags_inference_empty_table))

  estimates_arguments <- list(
    list(),
    list(probs = c(0.1, 0.9)),
    list(title = "Title", footnotes = "Note", warnings = "Warning"),
    list(remove_diagnostics = TRUE),
    list(diagnostic_columns = c("ESS", "R_hat")),
    list(probs = 0.5, diagnostic_columns = TRUE)
  )
  for(arguments in estimates_arguments){
    expect_identical(
      do.call(JAGS_estimates_empty_table, arguments),
      do.call(runjags_estimates_empty_table, arguments)
    )
  }
  inference_arguments <- list(
    list(),
    list(title = "Title", footnotes = "Note", warnings = "Warning"),
    list(logBF = TRUE),
    list(BF01 = TRUE),
    list(BF_diagnostics = TRUE, logBF = TRUE, BF01 = TRUE),
    list(BF_diagnostic_columns = "ESS")
  )
  for(arguments in inference_arguments){
    expect_identical(
      do.call(JAGS_inference_empty_table, arguments),
      do.call(runjags_inference_empty_table, arguments)
    )
  }

  # the tables have no rows and the columns that the arguments request
  estimates <- JAGS_estimates_empty_table(
    probs = c(0.1, 0.9), diagnostic_columns = c("ESS", "R_hat"), title = "Title"
  )
  expect_s3_class(estimates, "BayesTools_table")
  expect_s3_class(estimates, "BayesTools_runjags_summary")
  expect_identical(nrow(estimates), 0L)
  expect_identical(colnames(estimates), c("Mean", "SD", "0.1", "0.9", "ESS", "R_hat"))
  expect_identical(attr(estimates, "type"), c(rep("estimate", 4L), "ESS", "R_hat"))
  expect_identical(attr(estimates, "title"), "Title")

  inference <- JAGS_inference_empty_table(logBF = TRUE, BF01 = TRUE, BF_diagnostics = TRUE)
  expect_s3_class(inference, "BayesTools_table")
  expect_s3_class(inference, "BayesTools_runjags_inference")
  expect_identical(nrow(inference), 0L)
  expect_identical(
    colnames(inference),
    c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error_percent")
  )
  expect_identical(attr(inference[["inclusion_BF"]], "name"), "log(Exclusion BF)")
  expect_identical(attr(inference[["BF_error_percent"]], "name"), "error%(Exclusion BF)")
})


# ---------------------------------------------------------------------------- #
# hypothesis_alias_label(), hypothesis_unique_named_warnings()
# ---------------------------------------------------------------------------- #

test_that("hypothesis_alias_label names a parameter by its first alias other than itself", {

  # the first alias of the parameter in the order of the list
  expect_identical(hypothesis_alias_label(list(effect = "mu", location = "mu"), "mu"), "effect")
  # the parameter's own name is no alias, wherever it stands
  expect_identical(hypothesis_alias_label(list(mu = "mu", effect = "mu"), "mu"), "effect")
  expect_identical(hypothesis_alias_label(list(effect = "mu", mu = "mu"), "mu"), "effect")
  # aliases of other parameters are ignored
  expect_identical(hypothesis_alias_label(list(spread = "tau", effect = "mu", sd = "tau"), "mu"), "effect")
  expect_identical(hypothesis_alias_label(list(spread = "tau", effect = "mu", sd = "tau"), "tau"), "spread")

  # duplicated aliases name the parameter once
  expect_identical(hypothesis_alias_label(list(effect = "mu", effect = "mu"), "mu"), "effect")
  expect_identical(hypothesis_alias_label(list(mu = "mu", mu = "mu", effect = "mu"), "mu"), "effect")

  # without an alias, the parameter itself
  expect_identical(hypothesis_alias_label(list(), "mu"), "mu")
  expect_identical(hypothesis_alias_label(list(mu = "mu"), "mu"), "mu")
  expect_identical(hypothesis_alias_label(list(spread = "tau"), "mu"), "mu")
  # unnamed and empty names are no alias
  expect_identical(hypothesis_alias_label(list("mu"), "mu"), "mu")
  expect_identical(hypothesis_alias_label(stats::setNames(list("mu", "mu"), c("", "effect")), "mu"), "effect")
  expect_identical(hypothesis_alias_label(stats::setNames(list("mu"), ""), "mu"), "mu")

  # the aliases and the parameter are validated
  expect_error(hypothesis_alias_label(c(effect = "mu"), "mu"), "'aliases' must be a list.", fixed = TRUE)
  expect_error(hypothesis_alias_label(list(effect = "mu"), c("mu", "tau")), "length '1'")
  expect_error(hypothesis_alias_label(list(effect = "mu"), NA_character_), "NA")
  expect_error(hypothesis_alias_label(list(effect = "mu"), 1), "character")
})

test_that("hypothesis_unique_named_warnings drops repeated warnings of the same row", {

  # the same text of the same row once, in the order of first occurrence and with
  # the names of the kept warnings
  warnings <- c(
    row1 = "first", row1 = "first", row2 = "first", row1 = "second", row2 = "first"
  )
  expect_identical(
    hypothesis_unique_named_warnings(warnings),
    c(row1 = "first", row2 = "first", row1 = "second")
  )
  # warnings that are already unique are returned as they are
  unique_warnings <- c(row1 = "first", row2 = "second", row3 = "third")
  expect_identical(hypothesis_unique_named_warnings(unique_warnings), unique_warnings)
  expect_identical(
    hypothesis_unique_named_warnings(hypothesis_unique_named_warnings(warnings)),
    hypothesis_unique_named_warnings(warnings)
  )

  # unnamed warnings are compared by their text; a name makes a warning another one
  expect_identical(hypothesis_unique_named_warnings(c("a", "a", "b", "a")), c("a", "b"))
  expect_identical(hypothesis_unique_named_warnings(c("a", row = "a")), c("a", row = "a"))
  expect_identical(
    hypothesis_unique_named_warnings(c(row = "a", "a", row = "a")),
    c(row = "a", "a")
  )

  # missing and empty warnings are dropped, and nothing left is NULL
  expect_identical(
    hypothesis_unique_named_warnings(c(row1 = "", row2 = NA, row3 = "kept", row4 = "kept")),
    c(row3 = "kept", row4 = "kept")
  )
  expect_null(hypothesis_unique_named_warnings(c(row1 = "", row2 = NA)))
  expect_null(hypothesis_unique_named_warnings(character()))
  expect_null(hypothesis_unique_named_warnings(NULL))
})


# ---------------------------------------------------------------------------- #
# prior_informed_medicine_names
# ---------------------------------------------------------------------------- #

test_that("prior_informed_medicine_names lists the subfields that prior_informed() accepts", {

  subfields <- prior_informed_medicine_names
  types     <- c("smd", "logOR", "logRR", "RD", "logHR")

  expect_type(subfields, "character")
  expect_false(anyNA(subfields))
  expect_true(all(nzchar(subfields)))
  expect_identical(anyDuplicated(subfields), 0L)
  # prior_informed() matches names after dropping case and punctuation, so the
  # listed names have to stay distinguishable after it
  expect_identical(anyDuplicated(.prior_clean_input_name(subfields)), 0L)
  expect_true("Cochrane" %in% subfields)
  # the psychology priors are other names that need no parameter
  expect_false(any(c("van Erp", "Oosterwijk") %in% subfields))

  # a listed subfield has an effect prior for a type, or is refused as
  # unavailable for it, and no subfield is unavailable for every type
  outcome <- function(subfield, type){
    prior <- tryCatch(
      prior_informed(subfield, parameter = "effect", type = type),
      error = identity
    )
    if(is.prior(prior)){
      return("prior")
    }
    if(grepl("not available for type", conditionMessage(prior), fixed = TRUE)){
      return("unavailable")
    }
    conditionMessage(prior)
  }
  outcomes <- vapply(types, function(type){
    vapply(subfields, outcome, character(1L), type = type)
  }, character(length(subfields)))

  expect_identical(setdiff(as.vector(outcomes), c("prior", "unavailable")), character())
  expect_identical(subfields[rowSums(outcomes == "prior") == 0L], character())
  # the whole database has every type, the subfields of hazard ratios only that
  expect_true(all(outcomes["Cochrane", ] == "prior"))
  expect_identical(subfields[outcomes[, "logHR"] == "prior"], "Cochrane")
  heterogeneity <- lapply(types, function(type){
    prior_informed("Cochrane", parameter = "heterogeneity", type = type)
  })
  expect_true(all(vapply(heterogeneity, is.prior, logical(1L))))

  # punctuation of a listed name does not have to be typed
  punctuated <- grep("[,;/-]", subfields, value = TRUE)
  expect_gt(length(punctuated), 0L)
  expect_identical(
    lapply(punctuated, prior_informed, parameter = "effect", type = "logOR"),
    lapply(tolower(gsub("[,;/-]", "", punctuated)), prior_informed, parameter = "effect", type = "logOR")
  )

  # an unlisted name is unknown
  expect_error(prior_informed("Not a subfield", parameter = "effect"), "unknown prior 'name'")
})


# ---------------------------------------------------------------------------- #
# rtwo.sided_fixed()
# ---------------------------------------------------------------------------- #

test_that("rtwo.sided_fixed returns the fixed weights as every draw", {

  omega <- c(1, 0.5, 1.5)

  # every draw is the weights, and no random number is used
  set.seed(1)
  seed  <- .Random.seed
  draws <- rtwo.sided_fixed(4, omega = omega)
  expect_identical(.Random.seed, seed)
  expect_identical(draws, matrix(rep(omega, each = 4L), nrow = 4L))
  expect_identical(rtwo.sided_fixed(1, omega = omega), matrix(omega, nrow = 1L))

  # a matrix of weights has one row per draw
  omega_rows <- matrix(c(1, 0.3, 0.8, 1, 0.6, 0.2), nrow = 2L, byrow = TRUE)
  expect_identical(rtwo.sided_fixed(2, omega = omega_rows), omega_rows)
  expect_identical(
    rtwo.sided_fixed(3, omega = omega_rows[1L, , drop = FALSE]),
    matrix(rep(omega_rows[1L, ], each = 3L), nrow = 3L)
  )

  # the draws of the one-sided function, of the quantile function of the same
  # fixed weights at any probability, and of the fixed two-sided prior
  expect_identical(draws, rone.sided_fixed(4, omega = omega))
  expect_equal(draws[1L, , drop = TRUE], mqtwo.sided_fixed(0.37, omega = omega)[1L, ])
  prior <- prior_weightfunction("two-sided", c(0.01, 0.05), wf_fixed(omega))
  expect_equal(unname(rng(prior, 4)), draws)

  # the weights are validated as for the one-sided function
  expect_error(
    rtwo.sided_fixed(3, omega = omega_rows),
    "Incompatible dimensions of requested number of samples and 'omega'.",
    fixed = TRUE
  )
  expect_error(rtwo.sided_fixed(2, omega = c(0.9, 0.5)), "reference-bin")
  expect_error(rtwo.sided_fixed(2, omega = c(1, -0.5)), "'omega' must be non-negative.", fixed = TRUE)
  expect_error(rtwo.sided_fixed(2, omega = 1), "'omega' must be a vector of length at least 2.", fixed = TRUE)
  expect_error(rtwo.sided_fixed(2, omega = "1"), "'omega' must be a numeric vector or a matrix.", fixed = TRUE)
  expect_error(rtwo.sided_fixed(0, omega = omega), "'n' must be equal or higher than 1.", fixed = TRUE)
})


# ---------------------------------------------------------------------------- #
# transform_orthonormal_samples()
# ---------------------------------------------------------------------------- #

# Coefficient draws of a factor term with an orthonormal or mean-difference
# contrast, with the metadata that the mixed posterior of such a factor carries.
.contrast_samples_for_test <- function(levels, n, seed, contrast = "orthonormal", term = "f"){

  n_coefficients <- length(levels) - 1L
  coefficients <- withr::with_seed(seed, matrix(
    stats::rnorm(n * n_coefficients),
    nrow = n,
    dimnames = list(NULL, paste0(term, "[", seq_len(n_coefficients), "]"))
  ))
  samples <- coefficients
  class(samples) <- c("mixed_posteriors.factor", class(samples))
  attr(samples, contrast)            <- TRUE
  attr(samples, "level_names")       <- levels
  attr(samples, "factor_terms")      <- term
  attr(samples, "factor_contrasts")  <- stats::setNames(paste0("contr.", contrast), term)
  attr(samples, "levels")            <- n_coefficients

  list(coefficients = coefficients, samples = samples)
}

test_that("transform_orthonormal_samples maps the coefficients to the differences of levels from the mean", {

  for(levels in list(c("a", "b"), c("a", "b", "c", "d"))){
    n     <- 40L
    draws <- .contrast_samples_for_test(levels, n, seed = length(levels))
    cells <- transform_orthonormal_samples(list(f = draws$samples))$f

    # the cells are the coefficients combined by the rows of the orthonormal contrast
    contrast <- contr.orthonormal(levels)
    expect_identical(dim(cells), c(n, length(levels)))
    expect_equal_each(
      as.numeric(cells),
      as.numeric(draws$coefficients %*% t(contrast)),
      tolerance = 1e-12
    )

    # differences from the mean sum to zero, and the orthonormal columns of the
    # contrast recover the coefficients from them
    expect_equal_each(rowSums(cells), rep(0, n), tolerance = 1e-12)
    expect_equal_each(
      as.numeric(matrix(as.numeric(cells), nrow = n) %*% contrast),
      as.numeric(draws$coefficients),
      tolerance = 1e-10
    )

    # the cells are named and labelled by the levels, and keep the factor class
    expect_identical(colnames(cells), paste0("f[dif: ", levels, "]"))
    expect_identical(attr(cells, "level_names"), levels)
    expect_s3_class(cells, "mixed_posteriors.factor")
    expect_s3_class(cells, "mixed_posteriors.orthonormal_transformed")
  }
})

test_that("transform_orthonormal_samples transforms only orthonormal factor samples, once", {

  orthonormal <- .contrast_samples_for_test(c("a", "b", "c"), 20L, seed = 2L)
  meandif     <- .contrast_samples_for_test(c("x", "y", "z"), 20L, seed = 3L, contrast = "meandif", term = "g")
  samples     <- list(f = orthonormal$samples, g = meandif$samples, tau = c(0.1, 0.2, 0.3))

  transformed <- transform_orthonormal_samples(samples)
  expect_identical(names(transformed), names(samples))
  expect_identical(transformed$g, samples$g)
  expect_identical(transformed$tau, samples$tau)
  expect_identical(dim(transformed$f), c(20L, 3L))
  expect_false(identical(transformed$f, samples$f))
  expect_identical(transform_orthonormal_samples(transformed), transformed)

  # no samples, nothing to transform; anything but a list is refused
  expect_null(transform_orthonormal_samples(NULL))
  expect_identical(transform_orthonormal_samples(list()), list())
  expect_error(transform_orthonormal_samples(1:3), "'samples' argument must be a list", fixed = TRUE)
})

test_that("orthonormal and mean-difference cells differ by the scale of their contrasts", {

  levels      <- c("a", "b", "c", "d", "e")
  orthonormal <- .contrast_samples_for_test(levels, 30L, seed = 4L)
  meandif     <- .contrast_samples_for_test(levels, 30L, seed = 4L, contrast = "meandif", term = "g")
  expect_identical(unname(orthonormal$coefficients), unname(meandif$coefficients))

  orthonormal_cells <- transform_orthonormal_samples(list(f = orthonormal$samples))$f
  meandif_cells     <- transform_meandif_samples(list(g = meandif$samples))$g

  # contr.meandif(K) is contr.orthonormal(K) / sqrt(1 - 1 / K)
  expect_equal_each(
    as.numeric(meandif_cells) * sqrt(1 - 1 / length(levels)),
    as.numeric(orthonormal_cells),
    tolerance = 1e-12
  )

  # the factor transformation applies each one to the samples of its contrast
  both <- transform_factor_samples(list(f = orthonormal$samples, g = meandif$samples))
  expect_equal_each(as.numeric(both$f), as.numeric(orthonormal_cells), tolerance = 1e-12)
  expect_equal_each(as.numeric(both$g), as.numeric(meandif_cells), tolerance = 1e-12)
  expect_s3_class(both$f, "mixed_posteriors.orthonormal_transformed")
  expect_s3_class(both$g, "mixed_posteriors.meandif_transformed")
})
