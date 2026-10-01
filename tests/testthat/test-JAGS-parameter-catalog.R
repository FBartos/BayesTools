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
  expect_identical(catalog$schema_version, 9L)
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
    catalog$quantities$canonical_name %in% c("mu_f[b]", "mu_f[c]"),
    ,
    drop = FALSE
  ]
  expect_identical(factor_rows$component, c("b", "c"))
  expect_identical(factor_rows$display_label, c("(mu) f[b]", "(mu) f[c]"))
  expect_identical(
    lapply(factor_rows$extraction_key, `[[`, "dependencies"),
    list("mu_f[1]", "mu_f[2]")
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      alias = "f",
      namespace = "mu",
      component = "b"
    )$quantities$canonical_name,
    "mu_f[b]"
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog,
      alias = "f[c]",
      namespace = "mu"
    )$quantities$canonical_name,
    "mu_f[c]"
  )
  # Backend coordinate names are not selectors of factor levels.
  for(coordinate in c("mu_f[1]", "mu_f[2]")){
    expect_error(
      parameter_catalog_resolve(catalog, alias = coordinate),
      class = "BayesTools_parameter_not_found"
    )
  }

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
    c("mu_f[b]", "mu_f[c]")
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
  condition <- tryCatch(
    hypothesis_resolve(
      hypothesis_parse("f[b] > 0"),
      catalog,
      namespace = "mu",
      component = "c"
    ),
    error = identity
  )
  expect_identical(
    class(condition),
    c("BayesTools_hypothesis_component_mismatch",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  expect_identical(
    conditionMessage(condition),
    paste0(
      "The level in hypothesis symbol 'f[b]' does not match the requested ",
      "catalog component 'c'."
    )
  )
  expect_identical(condition$symbol, "f[b]")
  expect_identical(condition$component, "c")

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

test_that("factor cell mapping covers fixed factor priors, not random-effect priors", {

  # Ordinary factor priors label their levels 1..K by construction, so
  # `p1[k]` is level k and never the k-th backend coordinate.
  ordinary <- prior_factor(
    "mnormal",
    list(0, 1),
    contrast = "orthonormal"
  )
  ordinary <- prior_factor_levels(ordinary, 3L)
  ordinary_coordinates <- .bt_build_parameter_coordinates(
    .JAGS_prior_factor_names("p1", ordinary),
    prior_list = list(p1 = ordinary)
  )
  ordinary_catalog <- .bt_build_parameter_catalog(
    ordinary_coordinates,
    prior_list = list(p1 = ordinary)
  )
  expect_true(all(vapply(
    ordinary_catalog$quantities$extraction_key,
    function(key) identical(key$type, "factor_level"),
    logical(1)
  )))
  expect_setequal(
    ordinary_catalog$quantities$canonical_name,
    c("p1[1]", "p1[2]", "p1[3]", "p1{1}", "p1{2}")
  )
  ordinary_design <- .factor_term_design_from_metadata(
    .complete_factor_metadata(ordinary, "p1")
  )$design
  for(level in 1:3){
    level_key <- parameter_catalog_resolve(
      ordinary_catalog,
      paste0("p1[", level, "]")
    )$quantities$extraction_key[[1L]]
    expected_weights <- ordinary_design[level, ]
    expect_identical(level_key$dependencies,
                     c("p1[1]", "p1[2]")[expected_weights != 0])
    expect_identical(level_key$weights,
                     expected_weights[expected_weights != 0])
  }
  treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")
  treatment <- prior_factor_levels(treatment, 3L)
  treatment_catalog <- .bt_build_parameter_catalog(
    .bt_build_parameter_coordinates(
      .JAGS_prior_factor_names("p1", treatment),
      prior_list = list(p1 = treatment)
    ),
    prior_list = list(p1 = treatment)
  )
  expect_identical(
    parameter_catalog_resolve(treatment_catalog, "p1[1]")$quantities$source_type,
    "structural_zero"
  )
  expect_identical(
    parameter_catalog_resolve(
      treatment_catalog,
      "p1[2]"
    )$quantities$extraction_key[[1L]]$dependencies,
    "p1[1]"
  )
  expect_identical(
    parameter_catalog_resolve(
      treatment_catalog,
      "p1[3]"
    )$quantities$extraction_key[[1L]]$dependencies,
    "p1[2]"
  )

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
  expect_identical(reference$quantities$status, "structural")
  expect_identical(reference$quantities$fixed_value, 0)
  expect_identical(reference$quantities$source_type, "structural_zero")
  direct <- parameter_catalog_resolve(
    treatment_catalog,
    "f[b]",
    namespace = "mu"
  )
  expect_identical(direct$quantities$canonical_name, "mu_f[b]")
  expect_identical(direct$quantities$source_type, "identity")
  expect_identical(direct$quantities$status, "sampled")
  expect_identical(direct$quantities$extraction_key[[1L]]$dependencies, "mu_f[1]")
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
      paste0("mu_f[", levels(data$f)[level_i], "]")
    )
    expect_identical(
      selection$quantities$extraction_key[[1L]]$dependencies,
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
    # Contrast coefficients are `f{j}`, never a bracketed level selector.
    coefficient_values <- matrix(c(1, 10, 2, 20), ncol = 2L, byrow = TRUE)
    for(coefficient in 1:2){
      selection <- parameter_catalog_resolve(
        catalog,
        paste0("f{", coefficient, "}"),
        namespace = "mu"
      )
      expect_identical(selection$quantities$canonical_name,
                       paste0("mu_f{", coefficient, "}"), info = contrast)
      expect_identical(selection$quantities$status, "sampled", info = contrast)
      expect_identical(
        as.numeric(parameter_draws(transformed$fit, selection)[[1L]][, 1L]),
        coefficient_values[, coefficient],
        info = contrast
      )
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
    "mu_f[b]"
  )
  # Later ordered coordinates are increments, not levels: coefficient 2.
  expect_identical(
    parameter_catalog_resolve(
      ordered_catalog,
      "f{2}",
      namespace = "mu"
    )$quantities$extraction_key[[1L]]$dependencies,
    "mu_f[2]"
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
  # The first cumulative-levels coordinate is structurally level "a"; the
  # later coordinates are increments, coefficients 2 and 3.
  ordered_levels_a <- parameter_catalog_resolve(
    ordered_levels_catalog,
    "f[a]",
    namespace = "mu"
  )$quantities
  expect_identical(ordered_levels_a$status, "sampled")
  expect_identical(ordered_levels_a$extraction_key[[1L]]$dependencies, "mu_f[1]")
  expect_setequal(
    ordered_levels_catalog$quantities$canonical_name[
      ordered_levels_catalog$quantities$term == "f"
    ],
    c("mu_f[a]", "mu_f[b]", "mu_f[c]", "mu_f{2}", "mu_f{3}")
  )
})

# A formula factor fit with deterministic synthetic draws for the label
# contract: level selectors, coefficient selectors, and table rows.
.label_contract_fit <- function(levels, contrast, seed = 1L){

  factor_prior <- switch(
    contrast,
    treatment   = prior_factor("normal", list(0, 1), contrast = "treatment"),
    independent = prior_factor("normal", list(0, 1), contrast = "independent"),
    meandif     = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
    orthonormal = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
    ordered     = prior_ordered(prior("normal", list(0, 1)))
  )
  g <- factor(rep(levels, 3L), levels = levels)
  if(identical(contrast, "ordered")){
    g <- ordered(g, levels = levels)
  }
  data <- data.frame(g = g)
  formula <- if(identical(contrast, "independent")) ~ 0 + g else ~ 1 + g
  priors <- list(g = factor_prior)
  if(!identical(contrast, "independent")){
    priors <- c(list(intercept = prior("normal", list(0, 1))), priors)
  }
  result <- JAGS_formula(formula, "mu", data, priors)
  columns <- unlist(lapply(names(result$prior_list), function(name){
    prior <- result$prior_list[[name]]
    if(.bt_prior_is_factor_family(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  set.seed(seed)
  samples <- matrix(
    stats::rnorm(20L * length(columns)),
    nrow = 20L,
    dimnames = list(NULL, columns)
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = 20L),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attach_test_parameter_map(fit)
}

.label_contract_level_sets <- function(){
  list(
    numeric = c(5, 10, 20),
    index = 1:4,
    character = c("a", "b", "c")
  )
}

test_that("factor selectors name level labels, never coordinate positions", {

  contrasts <- c("treatment", "independent", "meandif", "orthonormal",
                 "ordered")
  for(contrast in contrasts){
    for(level_set in names(.label_contract_level_sets())){
      levels <- as.character(.label_contract_level_sets()[[level_set]])
      info <- paste(contrast, level_set)
      fit <- .label_contract_fit(levels, contrast)
      catalog <- parameter_catalog(fit)
      samples <- as.matrix(fit$mcmc)

      # Independent reference: marginal_posterior() evaluates each level
      # through the fitted model matrix; without the intercept it is the
      # level-L term value.
      mixed <- as_mixed_posteriors(
        fit,
        parameters = names(attr(fit, "prior_list"))
      )
      marginal <- marginal_posterior(
        mixed,
        parameter = "mu_g",
        formula = ~ g
      )
      # every level, position and alias of the unit is checked; the failures
      # are collected and asserted once per unit
      problems <- expectation_problems({
        expect_identical(names(marginal), levels, info = info)
        intercept <- if("mu_intercept" %in% colnames(samples)){
          samples[, "mu_intercept"]
        }else{
          0
        }

        level_checks <- lapply(levels, function(level){
          selectors <- c(
            paste0("g[", level, "]"),
            paste0("mu_g[", level, "]"),
            paste0("(mu) g[", level, "]"),
            paste0("mu_g[dif: ", level, "]"),
            paste0("g[dif: ", level, "]"),
            paste0("(mu) g[dif: ", level, "]")
          )
          selections <- lapply(selectors, function(selector){
            parameter_catalog_resolve(catalog, selector)
          })
          ids <- vapply(selections, `[[`, character(1), "quantity_id")
          quantity <- selections[[1L]]$quantities
          hypothesis <- hypothesis_resolve(
            hypothesis_parse(paste0("mu_g[", level, "] > g[", level, "]")),
            catalog
          )
          list(
            one_quantity        = identical(unique(ids), ids[[1L]]),
            canonical_name      = quantity$canonical_name,
            component           = quantity$component,
            draws               = as.numeric(as.matrix(parameter_draws(fit, selections[[1L]]))),
            hypothesis_quantity = identical(unique(hypothesis$occurrences$quantity_id), ids[[1L]])
          )
        })
        # every selector of a level is the one quantity of the level
        expect_true(
          all(vapply(level_checks, `[[`, logical(1), "one_quantity")),
          info = info
        )
        expect_identical(
          vapply(level_checks, `[[`, character(1), "canonical_name"),
          paste0("mu_g[", levels, "]"),
          info = info
        )
        expect_identical(
          vapply(level_checks, `[[`, character(1), "component"),
          levels,
          info = info
        )
        # Exact linear algebra on the same draws: equality up to rounding.
        for(k in seq_along(levels)){
          expect_equal(
            level_checks[[k]]$draws,
            as.numeric(marginal[[levels[[k]]]]) - intercept,
            tolerance = 1e-12,
            info = paste(info, levels[[k]])
          )
        }
        expect_true(
          all(vapply(level_checks, `[[`, logical(1), "hypothesis_quantity")),
          info = info
        )

        # Positions that are not level labels select nothing.
        positions <- setdiff(as.character(1:5), levels)
        position_selectors <- as.vector(rbind(
          paste0("mu_g[", positions, "]"),
          paste0("g[", positions, "]")
        ))
        not_found <- vapply(position_selectors, function(selector){
          inherits(
            tryCatch(parameter_catalog_resolve(catalog, selector), error = identity),
            "BayesTools_parameter_not_found"
          )
        }, logical(1))
        expect_true(
          all(not_found),
          info = paste(info, "positions found:", paste(position_selectors[!not_found], collapse = ", "))
        )

        # No alias other than the shared term name is ambiguous.
        aliases <- setdiff(unique(catalog$aliases$alias), "g")
        # every error counts, also one without a message
        refusals <- as.character(unlist(lapply(aliases, function(alias){
          error <- tryCatch({
            parameter_catalog_resolve(catalog, alias)
            NULL
          }, error = identity)
          if(is.null(error)) NULL else paste0(alias, ": ", conditionMessage(error))
        })))
        expect_identical(refusals, character(), info = info)
      })
      expect_identical(problems, character(), info = info)
    }
  }
})

# The selector, parsing, and resolution refusals of a contrast-coefficient
# selector `form` whose coordinate is the level `level_form`.
.expect_level_coefficient_refused <- function(catalog, form, level_form,
                                              coefficient, term, info){

  message <- paste0(
    "Selector '", form, "' is unavailable: coefficient ", coefficient,
    " of factor term '", term, "' is a level, not a contrast coefficient. ",
    "Select it by its level label, '", level_form, "'."
  )
  # the condition of code that must fail, NULL when it does not
  refusal <- function(code){
    tryCatch({
      code
      NULL
    }, error = identity)
  }
  refused <- function(condition, exact){
    inherits(condition, "BayesTools_selector_unavailable") &&
      if(exact) identical(conditionMessage(condition), message)
      else grepl(message, conditionMessage(condition), fixed = TRUE)
  }
  error <- refusal(parameter_catalog_resolve(catalog, form))
  # parsing against the catalog refuses it instead of failing to parse, and a
  # quoted selector parses and is refused by catalog resolution
  parse_error <- refusal(hypothesis_parse(paste(form, "> 0"), catalog = catalog))
  quoted_error <- refusal(
    hypothesis_resolve(hypothesis_parse(paste0("`", form, "` > 0")), catalog)
  )
  expect_identical(
    c(resolve = refused(error, TRUE), parse = refused(parse_error, FALSE),
      quoted = refused(quoted_error, FALSE)),
    c(resolve = TRUE, parse = TRUE, quoted = TRUE),
    info = info
  )
  expect_true(
    inherits(error, "BayesTools_hypothesis_target") &&
      inherits(error, "BayesTools_parameter_resolution_error"),
    info = info
  )
  expect_identical(
    list(selector = error$selector, level = error$level, quantity_id = error$quantity_id),
    list(selector = form, level = level_form,
         quantity_id = parameter_catalog_resolve(catalog, level_form)$quantity_id),
    info = info
  )
}

test_that("contrast-coefficient selectors of level coordinates name the level", {

  # Treatment and independent coordinates, and the first ordered coordinate,
  # are level cells: `g{j}` names the level of coordinate j.
  levels <- c("5", "10", "20")
  coordinate_levels <- list(
    treatment   = c("10", "20"),
    independent = c("5", "10", "20"),
    ordered     = "10"
  )
  for(contrast in names(coordinate_levels)){
    catalog <- parameter_catalog(.label_contract_fit(levels, contrast))
    for(j in seq_along(coordinate_levels[[contrast]])){
      level <- coordinate_levels[[contrast]][[j]]
      forms <- c(
        paste0("g{", j, "}"),
        paste0("mu_g{", j, "}"),
        paste0("(mu) g{", j, "}")
      )
      level_forms <- c(
        paste0("g[", level, "]"),
        paste0("mu_g[", level, "]"),
        paste0("(mu) g[", level, "]")
      )
      for(k in seq_along(forms)){
        .expect_level_coefficient_refused(
          catalog, forms[[k]], level_forms[[k]], coefficient = j, term = "g",
          info = paste(contrast, forms[[k]])
        )
      }
    }
    # a coefficient beyond the coordinates selects nothing
    expect_error(
      parameter_catalog_resolve(catalog, paste0("g{", length(levels) + 1L, "}")),
      class = "BayesTools_parameter_not_found"
    )
  }

  # Later ordered increments and mean-difference and orthonormal coordinates
  # are contrast coefficients and stay selectable.
  ordered <- parameter_catalog(.label_contract_fit(levels, "ordered"))
  expect_identical(
    parameter_catalog_resolve(ordered, "g{2}")$quantities$canonical_name,
    "mu_g{2}"
  )
  for(contrast in c("meandif", "orthonormal")){
    catalog <- parameter_catalog(.label_contract_fit(levels, contrast))
    for(form in c("g{1}", "mu_g{2}", "(mu) g{1}")){
      expect_no_error(parameter_catalog_resolve(catalog, form))
      expect_no_error(hypothesis_parse(paste(form, "> 0"), catalog = catalog))
    }
  }

  # Interaction cells and ordinary factor priors (levels 1..K).
  data <- data.frame(
    g1 = factor(rep(levels, 4L), levels = levels),
    g2 = factor(rep(c("a", "b"), each = 6L))
  )
  treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")
  formula_result <- JAGS_formula(~ 1 + g1 * g2, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    g1        = treatment,
    g2        = treatment,
    "g1:g2"   = treatment
  ))
  prior_list <- formula_result$prior_list
  prior_list$p <- prior_factor_levels(treatment, 3L)
  columns <- unlist(lapply(names(prior_list), function(name){
    prior <- prior_list[[name]]
    if(.bt_prior_is_factor_family(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(matrix(
      0.1, nrow = 2L, ncol = length(columns),
      dimnames = list(NULL, columns)
    ))),
    prior_list = prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- parameter_catalog(fit)
  .expect_level_coefficient_refused(
    catalog, "g1:g2{2}", "g1[20]:g2[b]", coefficient = 2L, term = "g1:g2",
    info = "interaction"
  )
  .expect_level_coefficient_refused(
    catalog, "(mu) g1:g2{1}", "(mu) g1[10]:g2[b]", coefficient = 1L,
    term = "g1:g2", info = "prefixed interaction"
  )
  .expect_level_coefficient_refused(
    catalog, "p{1}", "p[2]", coefficient = 1L, term = "p",
    info = "ordinary factor prior"
  )

  # The same factor term in two namespaces: a level coordinate in 'mu'
  # (treatment) and a contrast coefficient in 'sigma' (meandif). Within 'mu',
  # `g{1}` is refused although it names the 'sigma' coefficient elsewhere;
  # without a namespace filter it keeps resolving to that coefficient.
  data <- data.frame(g = factor(rep(levels, 4L), levels = levels))
  formula_results <- list(
    mu = JAGS_formula(~ 1 + g, "mu", data, list(
      intercept = prior("normal", list(0, 1)),
      g         = treatment
    )),
    sigma = JAGS_formula(~ 1 + g, "sigma", data, list(
      intercept = prior("normal", list(0, 1)),
      g         = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    ))
  )
  prior_list <- do.call(c, unname(lapply(formula_results, `[[`, "prior_list")))
  columns <- unlist(lapply(names(prior_list), function(name){
    prior <- prior_list[[name]]
    if(.bt_prior_is_factor_family(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(matrix(
      0.1, nrow = 2L, ncol = length(columns),
      dimnames = list(NULL, columns)
    ))),
    prior_list = prior_list,
    formula_design = lapply(formula_results, `[[`, "formula_design")
  )
  catalog <- parameter_catalog(fit)
  message <- paste0(
    "Selector 'g{1}' is unavailable: coefficient 1 of factor term 'g' is a ",
    "level, not a contrast coefficient. Select it by its level label, ",
    "'g[10]'."
  )
  error <- expect_error(
    parameter_catalog_resolve(catalog, "g{1}", namespace = "mu"),
    message,
    fixed = TRUE,
    class = "BayesTools_selector_unavailable"
  )
  expect_identical(
    error$quantity_id,
    parameter_catalog_resolve(catalog, "g[10]", namespace = "mu")$quantity_id
  )
  expect_error(
    hypothesis_parse("g{1} > 0", catalog = catalog, namespace = "mu"),
    message,
    fixed = TRUE,
    class = "BayesTools_selector_unavailable"
  )
  .expect_level_coefficient_refused(
    catalog, "(mu) g{1}", "(mu) g[10]", coefficient = 1L, term = "g",
    info = "prefixed selector beside another namespace"
  )
  for(namespace in list(NULL, "sigma")){
    expect_identical(
      parameter_catalog_resolve(catalog, "g{1}", namespace = namespace)$
        quantities$canonical_name,
      "sigma_g{1}"
    )
    expect_no_error(
      hypothesis_parse("g{1} > 0", catalog = catalog, namespace = namespace)
    )
  }

  # Without a catalog the selector is not hypothesis syntax.
  error <- expect_error(hypothesis_parse("g{1} > 0"))
  expect_false(inherits(error, "BayesTools_selector_unavailable"))

  # Every refused selector contains `{`: parsing and resolving text without
  # one never builds the refused selectors (their cost grows with the level
  # cells of the catalog).
  local_mocked_bindings(
    .bt_parameter_catalog_level_coefficient_selectors = function(catalog){
      stop("the refused selectors were built", call. = FALSE)
    }
  )
  expect_no_error(hypothesis_parse("(mu) g[10] > 0", catalog = catalog))
  expect_error(
    parameter_catalog_resolve(catalog, "g[99]"),
    class = "BayesTools_parameter_not_found"
  )
})

# Every displayed row of the factor term resolves to the catalog quantity
# whose draws produced it: the row mean is the mean of the selected draws.
# 'draw_means' holds the mean of the draws of each quantity already extracted
# (the rows of the tables of one fit select the same few quantities).
.expect_factor_rows_resolve <- function(table, fit, catalog, info,
                                        draw_means = new.env(parent = emptyenv())){

  rows <- setdiff(rownames(table), c("(mu) intercept", "intercept"))
  expect_true(length(rows) > 0L, info = info)
  selections <- lapply(rows, function(row) parameter_catalog_resolve(catalog, row))
  expect_identical(
    vapply(selections, function(selection) selection$quantities$term, character(1)),
    rep("g", length(rows)),
    info = info
  )
  # The table summarizes the same draws: agreement up to rounding.
  expect_equal_each(
    vapply(selections, function(selection){
      id <- selection$quantity_id
      if(is.null(draw_means[[id]])){
        draw_means[[id]] <- mean(as.matrix(parameter_draws(fit, selection)))
      }
      draw_means[[id]]
    }, numeric(1)),
    table[rows, "Mean"],
    tolerance = 1e-10,
    info = info
  )
}

test_that("displayed factor rows resolve to the quantities that produced them", {

  contrasts <- c("treatment", "independent", "meandif", "orthonormal",
                 "ordered")
  for(contrast in contrasts){
    for(level_set in names(.label_contract_level_sets())){
      levels <- as.character(.label_contract_level_sets()[[level_set]])
      info <- paste(contrast, level_set)
      fit <- .label_contract_fit(levels, contrast)
      catalog <- parameter_catalog(fit)
      mixed <- as_mixed_posteriors(
        fit,
        parameters = names(attr(fit, "prior_list"))
      )
      # every row and label of the unit is checked; the failures are
      # collected and asserted once per unit
      draw_means <- new.env(parent = emptyenv())
      problems <- expectation_problems({
        for(transform in c(FALSE, TRUE)){
          for(prefix in c(TRUE, FALSE)){
            model_table <- JAGS_estimates_table(
              fit,
              transform_factors = transform,
              formula_prefix = prefix
            )
            ensemble_table <- ensemble_estimates_table(
              mixed,
              parameters = "mu_g",
              transform_factors = transform,
              formula_prefix = prefix
            )
            # A bracketed row always names a level label.
            rows <- c(rownames(model_table), rownames(ensemble_table))
            bracket_content <- sub("^[^[]*\\[(dif: )?(.*)\\]$", "\\2",
                                   grep("\\[", rows, value = TRUE))
            expect_true(all(bracket_content %in% levels),
                        info = paste(info, transform, prefix))
            .expect_factor_rows_resolve(
              model_table, fit, catalog,
              paste(info, "model", transform, prefix),
              draw_means
            )
            .expect_factor_rows_resolve(
              ensemble_table, fit, catalog,
              paste(info, "ensemble", transform, prefix),
              draw_means
            )
          }
        }

        # Coordinate display labels and mixed columns name the quantity that
        # is exactly that coordinate: its level cell or contrast coefficient.
        coordinates <- parameter_coordinates(fit)
        coordinates <- coordinates[coordinates$term == "g", , drop = FALSE]
        keys <- lapply(coordinates$display_label, function(label){
          parameter_catalog_resolve(catalog, label)$quantities$extraction_key[[1L]]
        })
        expect_identical(
          vapply(keys, function(key) key$dependencies, character(1)),
          coordinates$coordinate_name,
          info = info
        )
        expect_true(
          all(vapply(keys, function(key) identical(key$weights, 1), logical(1))),
          info = info
        )
        mixed_columns <- colnames(mixed[["mu_g"]])
        expect_identical(
          format_parameter_names(mixed_columns, formula_parameters = "mu"),
          coordinates$display_label,
          info = info
        )
      })
      expect_identical(problems, character(), info = info)
    }
  }
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
  # Transformed summaries name the same cell per factor.
  for(label in c("mu_f[dif: b]__xXx__g[dif: v]", "(mu) f[dif: b]:g[dif: v]",
                 "f[dif: b]:g[dif: v]")){
    expect_identical(
      parameter_catalog_resolve(catalog, label)$quantity_id,
      derived$quantity_id,
      info = label
    )
  }
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

test_that("estimates-table row labels resolve to their catalog quantities", {

  data <- data.frame(
    x = seq(1, 9, length.out = 24),
    h = factor(rep(c("lo", "hi"), 12), levels = c("lo", "hi")),
    g = factor(rep(c("u", "v", "w", "q"), 6), levels = c("u", "v", "w", "q")),
    f = factor(rep(c("a", "b", "c"), each = 8))
  )
  table_fit <- function(formula, prior_list, formula_scale = NULL){
    result <- JAGS_formula(
      formula, "mu", data, prior_list,
      formula_scale = formula_scale
    )
    columns <- unlist(lapply(names(result$prior_list), function(name){
      prior <- result$prior_list[[name]]
      if(is.prior.factor(prior)){
        .JAGS_prior_factor_names(name, prior)
      }else{
        name
      }
    }), use.names = FALSE)
    samples <- matrix(
      seq_len(2L * length(columns)) / 10,
      nrow = 2L,
      dimnames = list(NULL, columns)
    )
    fit <- structure(
      list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = 2L),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(fit, "prior_list") <- result$prior_list
    attr(fit, "formula_design") <- list(mu = result$formula_design)
    if(!is.null(result$formula_scale)){
      attr(fit, "formula_scale") <- list(mu = result$formula_scale)
    }
    attach_test_parameter_map(fit)
  }
  expect_rows_resolve <- function(fit, expected, transform_scaled = FALSE){
    table <- JAGS_estimates_table(fit, transform_scaled = transform_scaled)
    catalog <- parameter_catalog(fit)
    for(row in names(expected)){
      expect_true(row %in% rownames(table), info = row)
      selection <- parameter_catalog_resolve(catalog, row)
      expect_identical(
        selection$quantities$canonical_name,
        expected[[row]],
        info = row
      )
      if(!transform_scaled){
        expect_equal(
          mean(as.matrix(parameter_draws(fit, selection))),
          table[row, "Mean"],
          info = row
        )
      }
    }
  }
  normal <- prior("normal", list(0, 1))

  # Two-level treatment interaction: one unindexed coefficient, labelled by
  # the level cell it is.
  expect_rows_resolve(
    table_fit(~ x * h, list(
      intercept = normal, x = normal,
      h = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "x:h" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )),
    c("(mu) x:h[hi]" = "mu_x__xXx__h[hi]", "(mu) h[hi]" = "mu_h[hi]")
  )
  # Mean-difference coefficients whose contrast row is a unit vector are
  # contrast coefficients: their rows are `{j}`, and positional `[j]` labels
  # never resolve.
  meandif_fit <- table_fit(~ x * g, list(
    intercept = normal, x = normal,
    g = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
    "x:g" = prior_factor("mnormal", list(0, 1), contrast = "meandif")
  ))
  expect_rows_resolve(
    meandif_fit,
    stats::setNames(
      c(paste0("mu_g{", 1:3, "}"), paste0("mu_x__xXx__g{", 1:3, "}")),
      c(paste0("(mu) g{", 1:3, "}"), paste0("(mu) x:g{", 1:3, "}"))
    )
  )
  meandif_catalog <- parameter_catalog(meandif_fit)
  for(j in 1:3){
    for(term in c("g", "x:g")){
      expect_error(
        parameter_catalog_resolve(meandif_catalog, paste0("(mu) ", term, "[", j, "]")),
        class = "BayesTools_parameter_not_found"
      )
      coefficient <- parameter_catalog_resolve(
        meandif_catalog,
        paste0("(mu) ", term, "{", j, "}")
      )$quantities
      expect_identical(
        coefficient$extraction_key[[1L]]$dependencies,
        paste0("mu_", sub(":", "__xXx__", term, fixed = TRUE), "[", j, "]")
      )
    }
  }
  # Treatment-by-treatment interaction cells.
  expect_rows_resolve(
    table_fit(~ f * g, list(
      intercept = normal,
      f = prior_factor("normal", list(0, 1), contrast = "treatment"),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )),
    c(
      "(mu) f[b]:g[v]" = "mu_f__xXx__g[f=b, g=v]",
      "(mu) f[c]:g[v]" = "mu_f__xXx__g[f=c, g=v]",
      "(mu) f[c]:g[q]" = "mu_f__xXx__g[f=c, g=q]"
    )
  )
  # Log-intercept label of transformed tables.
  log_formula <- ~ x
  attr(log_formula, "log(intercept)") <- TRUE
  expect_rows_resolve(
    table_fit(
      log_formula,
      list(intercept = prior("lognormal", list(0, 0.5)), x = normal),
      formula_scale = list(x = TRUE)
    ),
    c("(mu) exp(intercept)" = "mu_intercept", "(mu) x" = "mu_x"),
    transform_scaled = TRUE
  )

  # Table labels never take over a selector of another quantity: with
  # index-like level names, '(mu) g[2]' keeps naming its level cell.
  data$g <- factor(rep(1:4, 6))
  fit <- table_fit(~ g, list(
    intercept = normal,
    g = prior_factor("mnormal", list(0, 1), contrast = "meandif")
  ))
  catalog <- parameter_catalog(fit)
  base_aliases <- .bt_parameter_catalog_aliases(
    catalog$quantities,
    attr(fit, "formula_design")
  )
  owners <- function(aliases, keys){
    lapply(keys, function(key){
      sort(unique(aliases$quantity_id[
        paste(aliases$alias, aliases$namespace) == key
      ]))
    })
  }
  base_keys <- unique(paste(base_aliases$alias, base_aliases$namespace))
  expect_identical(
    owners(catalog$aliases, base_keys),
    owners(base_aliases, base_keys)
  )
  expect_identical(
    parameter_catalog_resolve(catalog, "(mu) g[2]")$quantities$component,
    "2"
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

# A synthetic fit of the block us(1 + x1 + ... + x[K - 1] | id) with an
# LKJ(eta) correlation prior, monitoring the fitted priors and LKJ primitives.
.lkj_block_fit <- function(K, eta,
                           sd = prior("normal", list(0, 1), list(0, Inf)),
                           formula_scale = NULL){

  set.seed(K)
  predictors <- paste0("x", seq_len(K - 1L))
  data <- as.data.frame(stats::setNames(
    lapply(predictors, function(predictor) stats::rnorm(8L, 1)),
    predictors
  ))
  data$id <- factor(rep(c("a", "b", "c", "d"), each = 2L))
  formula_result <- JAGS_formula(
    formula = stats::as.formula(paste0(
      "~ 1 + us(1 + ", paste(predictors, collapse = " + "), " | id)"
    )),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = formula_scale,
    prior_random = prior_random(
      id = random_block(sd = sd, cor = prior_lkj(eta = eta))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  list(
    formula_result = formula_result,
    random_term    = random_term,
    fit            = .prior_monitor_test_fit(formula_result, unique(c(
      JAGS_to_monitor(formula_result$prior_list),
      random_term$correlation$primitive_names
    )))
  )
}

# The prior density is exactly the LKJ marginal: exact ordinates equal to
# dbeta((r + 1) / 2, shape, shape) / 2 (tolerance: rounding of the Beta
# density, applied to each ordinate).
.expect_lkj_marginal_density <- function(density, shape, info){

  expect_s3_class(density, "prior_linear_density")
  expect_identical(
    attr(density, "parameter_prior_density")$source,
    "fitted_parameter_map",
    info = info
  )
  values <- c(-0.95, -0.4, 0, 0.3, 0.9)
  ordinates <- lapply(values, function(value) prior_density_ordinate(density, value))
  expect_true(
    all(vapply(ordinates, function(ordinate) isTRUE(ordinate$exact), logical(1))),
    info = info
  )
  reference <- stats::dbeta((values + 1) / 2, shape, shape) / 2
  density_values <- exp(vapply(ordinates, function(ordinate) ordinate$log_density, numeric(1)))
  expect_equal_each(density_values, reference, tolerance = 1e-12, info = info)
}

test_that("LKJ correlations have the exact LKJ marginal prior density", {

  skip_on_cran()
  # Reference (Lewandowski, Kurowicka and Joe, 2009): every off-diagonal
  # correlation r of an LKJ(eta) K x K correlation matrix has
  # (r + 1) / 2 ~ Beta(eta - 1 + K / 2, eta - 1 + K / 2). The prior draws of
  # the primitives (Beta draws with the shapes of the emitted JAGS syntax, as
  # the JAGS module samples them) pass through the catalog's correlation
  # evaluator, and each of 20 equal-width bins on (-1, 1) holds its exact Beta
  # probability p within 4 binomial standard errors sqrt(p (1 - p) / n).
  # The standardized deviation does not depend on n: with n = 1e5 the smallest
  # expected bin count is 288 (the edge bins of the shape 2.5 cases), the exact
  # Binomial(n, p) probability of a 4-SE excursion is at most 6.7e-5 per bin,
  # and the 340 bins expect 0.022 false alarms (0.022 with n = 1e6). A 4-SE
  # deviation is 4.2-7.0% of the probability of a central bin (1.3-2.2% with
  # n = 1e6).
  # The cases cover each dimension K = 2, 3, 5 and eta below, at, and above 1,
  # and the three boundary classes of the Beta margin: shape 0.5 (U-shaped,
  # singular at +-1), shape 1 (flat), and shape 2 and above (vanishing at +-1);
  # K = 5 at the default eta = 1 checks all ten entries, K = 3 at eta = 0.5 and
  # 2 the entries of every row on both sides of eta = 1. The exact margin
  # depends on K and eta only through its shape; K sets the chain of canonical
  # partial correlations, so every dimension is checked with all its entries.
  n <- 1e5
  breaks <- seq(-1, 1, length.out = 21L)
  cases <- data.frame(K = c(2L, 3L, 3L, 5L), eta = c(0.5, 0.5, 2, 1))
  for(case in seq_len(nrow(cases))){
    K <- cases$K[[case]]
    eta <- cases$eta[[case]]
    info <- paste0("K = ", K, ", eta = ", eta)
    shape <- eta - 1 + K / 2
    block <- .lkj_block_fit(K, eta)
    fit <- block$fit

    # every check of the unit is made; the failures are collected and
    # asserted once per unit
    problems <- expectation_problems({
      # The emitted canonical partial correlations of the first row, which are
      # the correlations r[1, j], have this Beta shape.
      syntax <- block$formula_result$formula_syntax
      emitted <- regmatches(
        syntax,
        gregexpr("lkj_alpha\\[[0-9]+\\] <- [0-9.]+", syntax)
      )[[1L]]
      alpha <- as.numeric(sub("^.* <- ", "", emitted))
      expect_length(alpha, K * (K - 1L) / 2L)
      first_row <- vapply(2:K, function(j) (j - 1L) * (j - 2L) / 2L + 1L,
                          numeric(1))
      expect_equal(alpha[first_row], rep(shape, K - 1L), info = info)

      catalog <- parameter_catalog(fit)
      correlations <- catalog$quantities[
        catalog$quantities$role == "random_correlation", , drop = FALSE
      ]
      expect_identical(nrow(correlations), as.integer(K * (K - 1L) / 2L),
                       info = info)
      raw <- transform_prior_samples(fit, n_samples = n,
                                     seed = as.integer(10L * K + 2 * eta))
      all_pairs <- .bt_random_effect_summary_correlation_samples(
        block$random_term,
        raw
      )$values
      expected <- diff(stats::pbeta((breaks + 1) / 2, shape, shape))
      for(i in seq_len(nrow(correlations))){
        name <- correlations$canonical_name[[i]]
        selection <- parameter_catalog_resolve(catalog, name)
        index <- correlations$extraction_key[[i]]$index
        # the evaluated pair is the catalog quantity's draws
        expect_true(
          identical(
            as.numeric(as.matrix(parameter_draws(
              fit,
              selection,
              model_samples = raw[1:100, , drop = FALSE]
            ))),
            all_pairs[1:100, index]
          ),
          info = paste(info, name, "the draws are the evaluated pair")
        )
        .expect_lkj_marginal_density(
          parameter_prior_density(fit, selection),
          shape,
          info = paste(info, name)
        )
        observed <- tabulate(
          findInterval(all_pairs[, index], breaks, rightmost.closed = TRUE,
                       all.inside = TRUE),
          20L
        ) / n
        expect_lt(
          max(abs(observed - expected) / sqrt(expected * (1 - expected) / n)),
          4,
          label = paste(info, name, "largest bin deviation in binomial SE")
        )
      }
    })
    expect_identical(problems, character(), info = info)
  }
})

test_that("LKJ correlation priors are exact where the correlation is an LKJ entry", {

  eta <- 1.5
  shape <- eta - 1 + 3 / 2
  names <- c("(mu) cor(intercept,x1)", "(mu) cor(intercept,x2)",
             "(mu) cor(x1,x2)")
  density_of <- function(fit, name){
    parameter_prior_density(
      fit,
      parameter_catalog_resolve(parameter_catalog(fit), name)
    )
  }

  # Scaled predictors: the original-scale slopes are the fitted slopes divided
  # by the predictor SDs, so cor(x1,x2) is the fitted correlation; the
  # intercept of centred predictors mixes the fitted intercept and slopes, so
  # its correlations combine the SDs and have no exact density.
  scaled <- .lkj_block_fit(3L, eta, formula_scale = list(x1 = TRUE, x2 = TRUE))
  .expect_lkj_marginal_density(density_of(scaled$fit, names[[3L]]), shape,
                               info = "scaled")
  expect_null(density_of(scaled$fit, names[[1L]]))
  expect_null(density_of(scaled$fit, names[[2L]]))
  # Standardized summaries drop the formula scale: every correlation is then
  # the fitted-scale LKJ entry.
  standardized <- scaled$fit
  attr(standardized, "formula_scale") <- list()
  for(name in names){
    .expect_lkj_marginal_density(density_of(standardized, name), shape,
                                 info = paste("standardized", name))
  }

  # Gated SDs (spike-and-slab priors with a point at zero): the primitives are
  # a priori independent of the SDs and their gates.
  gated_sd <- prior_spike_and_slab(
    prior("gamma", list(2, 2)),
    prior_inclusion = prior("spike", list(0.5))
  )
  gated <- .lkj_block_fit(3L, eta, sd = gated_sd)
  for(name in names){
    .expect_lkj_marginal_density(density_of(gated$fit, name), shape,
                                 info = paste("gated", name))
  }
  gated_scaled <- .lkj_block_fit(3L, eta, sd = gated_sd,
                                 formula_scale = list(x1 = TRUE, x2 = TRUE))
  .expect_lkj_marginal_density(density_of(gated_scaled$fit, names[[3L]]),
                               shape, info = "gated scaled")
  expect_null(density_of(gated_scaled$fit, names[[1L]]))

  # The original-scale cor(x1,x2) of the gated scaled block is the fitted
  # correlation wherever both SDs are positive and undefined elsewhere, so its
  # prior given definedness is the LKJ marginal: 10 bins within 4 binomial SE
  # of the Beta probabilities.
  n <- 20000L
  raw <- transform_prior_samples(gated_scaled$fit, n_samples = n, seed = 71L,
                                 formula_scale = list())
  selection <- parameter_catalog_resolve(
    parameter_catalog(gated_scaled$fit),
    names[[3L]]
  )
  original <- as.numeric(as.matrix(parameter_draws(
    gated_scaled$fit,
    selection,
    model_samples = raw
  )))
  fitted <- .bt_random_effect_summary_correlation_samples(
    gated_scaled$random_term,
    raw
  )$values[, 3L]
  sd_names <- gated_scaled$random_term$sd_parameter_names
  defined <- raw[, sd_names[[2L]]] > 0 & raw[, sd_names[[3L]]] > 0
  expect_true(any(!defined))
  expect_identical(is.na(original), !defined)
  # rounding of the covariance transform only
  expect_equal(original[defined], fitted[defined], tolerance = 1e-12)
  breaks <- seq(-1, 1, length.out = 11L)
  expected <- diff(stats::pbeta((breaks + 1) / 2, shape, shape))
  observed <- tabulate(
    findInterval(original[defined], breaks, rightmost.closed = TRUE,
                 all.inside = TRUE),
    10L
  ) / sum(defined)
  expect_lt(
    max(abs(observed - expected) /
          sqrt(expected * (1 - expected) / sum(defined))),
    4
  )
})

test_that("posterior plots draw the exact LKJ prior of us() correlations", {

  block <- .lkj_block_fit(3L, 2)
  # prior draws stand in for posterior draws of the block
  raw <- transform_prior_samples(block$fit, n_samples = 500L, seed = 5L)
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(raw)),
    prior_list     = block$formula_result$prior_list,
    formula_design = list(mu = block$formula_result$formula_design)
  )
  name <- "(mu) cor(intercept,x1)"
  posterior <- parameter_mixed_posterior(
    fit,
    parameter_catalog_resolve(parameter_catalog(fit), name)
  )
  expect_s3_class(posterior_metadata(posterior, "prior_density"),
                  "prior_linear_density")
  samples <- stats::setNames(list(posterior), name)
  plot <- plot_posterior(samples, name, prior = TRUE, plot_type = "ggplot")
  prior_layer <- ggplot2::ggplot_build(plot)$data[[1L]]
  interior <- prior_layer$x > -1 & prior_layer$x < 1
  expect_gt(sum(interior), 100L)
  # K = 3, eta = 2: (r + 1) / 2 ~ Beta(2.5, 2.5)
  expect_equal(
    prior_layer$y[interior],
    stats::dbeta((prior_layer$x[interior] + 1) / 2, 2.5, 2.5) / 2,
    tolerance = 1e-10
  )
})

test_that("posterior plots draw the posterior alone when the draws have no prior density", {

  skip_if_not(capabilities("png"))
  # The original-scale cor(intercept,x1) of centred scaled predictors has no
  # prior density: its draws declare prior_none(). A correlation has no point
  # mass, so its atoms are declared none.
  block <- .lkj_block_fit(3L, 2, formula_scale = list(x1 = TRUE, x2 = TRUE))
  raw <- transform_prior_samples(block$fit, n_samples = 500L, seed = 5L,
                                 formula_scale = list())
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(raw)),
    prior_list     = block$formula_result$prior_list,
    formula_design = list(mu = block$formula_result$formula_design),
    formula_scale  = list(mu = block$formula_result$formula_scale)
  )
  name <- "(mu) cor(intercept,x1)"
  posterior <- parameter_mixed_posterior(
    fit,
    parameter_catalog_resolve(parameter_catalog(fit), name)
  )
  expect_null(posterior_metadata(posterior, "prior_density"))

  # Draws without an atom declaration stop before any plot is drawn, without
  # announcing an omitted prior curve.
  undeclared <- posterior
  posterior_metadata(undeclared, "atoms") <- NULL
  for(plot_fun in list(plot_posterior, plot_marginal)){
    expect_no_warning(
      expect_error(
        plot_fun(stats::setNames(list(undeclared), name), name, prior = TRUE,
                 plot_type = "ggplot"),
        "Posterior atom status is unknown",
        fixed = TRUE
      ),
      class = "BayesTools_prior_curve_unavailable"
    )
  }

  posterior_metadata(posterior, "atoms") <- posterior_atom_attribute()
  samples <- stats::setNames(list(posterior), name)
  message <- paste0(
    "The prior density curve of '", name, "' is unavailable: its posterior ",
    "draws carry no prior density, as the quantity has no deterministic ",
    "prior-density route. The prior curve is omitted from the plot."
  )

  for(plot_function in c("plot_posterior", "plot_marginal")){
    plot_fun <- get(plot_function)
    # ggplot: the posterior layers of the plot without the prior
    plot <- NULL
    expect_warning(
      plot <- plot_fun(samples, name, prior = TRUE, plot_type = "ggplot"),
      message,
      fixed = TRUE,
      class = "BayesTools_prior_curve_unavailable"
    )
    posterior_only <- plot_fun(samples, name, prior = FALSE, plot_type = "ggplot")
    expect_identical(
      ggplot2::ggplot_build(plot)$data,
      ggplot2::ggplot_build(posterior_only)$data,
      info = plot_function
    )

    # base graphics: the same drawing as without the prior
    draw <- function(prior){
      file <- tempfile(fileext = ".png")
      grDevices::png(file, width = 480, height = 360)
      if(prior){
        expect_warning(
          plot_fun(samples, name, prior = TRUE),
          message,
          fixed = TRUE,
          class = "BayesTools_prior_curve_unavailable"
        )
      }else{
        plot_fun(samples, name, prior = FALSE)
      }
      grDevices::dev.off()
      on.exit(unlink(file))
      unname(tools::md5sum(file))
    }
    expect_identical(draw(TRUE), draw(FALSE), info = plot_function)
  }
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

test_that("prior draws carry the allocation-derived SDs of correlated random slopes", {

  data <- data.frame(
    x = c(1, 4, 6, 2, 8, 5, 3, 7),
    g = factor(rep(c("a", "b", "c", "d"), each = 2L))
  )
  eta <- 2
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + us(1 + x | g),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "het",
        terms = "g",
        target = "sd_component",
        scale = "mean_variance",
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      g = random_block(cor = prior_lkj(eta = eta))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  allocation <- random_term$sd_binding$allocations[[1L]]
  source_name <- allocation$source_node
  weight_names <- paste0(allocation$weight_name, "[", 1:2, "]")
  # The monitored SDs are deterministic nodes sd * sqrt(2 * w[k]).
  sd_names <- random_term$sd_parameter_names
  R_names <- .prior_monitor_matrix_names(random_term$correlation$correlation_name, 2L)
  columns <- c(
    "mu_intercept", "mu_x", source_name, weight_names,
    paste0(.JAGS_prior_dirichlet_eta_name(allocation$weight_name), "[", 1:2, "]"),
    sd_names,
    .prior_monitor_matrix_names(random_term$correlation$cholesky_name, 2L),
    R_names,
    random_term$correlation$primitive_names,
    as.vector(.bt_random_effect_latent_names(random_term, random_term$n_groups, 2L))
  )
  fit <- .prior_monitor_test_fit(formula_result, columns)

  n <- 10000L
  raw <- transform_prior_samples(fit, n_samples = n, seed = 58L,
                                 formula_scale = list())
  expect_true(all(sd_names %in% colnames(raw)))
  for(k in 1:2){
    expect_identical(
      unname(raw[, sd_names[[k]]]),
      unname(raw[, source_name] * sqrt(2 * raw[, weight_names[[k]]]))
    )
  }

  # Every public quantity, including the original-scale correlation that
  # depends on the SD monitors, is available from the prior draws.
  draws <- .prior_monitor_catalog_draws(fit, raw)
  expect_true("(mu) cor(intercept,x)" %in% names(draws))
  for(name in names(draws)){
    values <- as.numeric(as.matrix(draws[[name]]))
    expect_length(values, n)
    expect_true(all(is.finite(values)), info = name)
  }

  # Original-scale correlation of u0 - u1 m / s and u1 / s, from the fitted
  # correlation and SDs.
  scale_info <- formula_result$formula_scale$mu_x
  ratio <- scale_info$mean / scale_info$sd
  rho <- raw[, R_names[[2L]]]
  t0 <- raw[, sd_names[[1L]]]
  t1 <- raw[, sd_names[[2L]]]
  expected_cor <- (rho * t0 * t1 - ratio * t1^2) /
    sqrt((t0^2 - 2 * ratio * rho * t0 * t1 + ratio^2 * t1^2) * t1^2)
  expect_equal(
    unname(as.numeric(as.matrix(draws[["(mu) cor(intercept,x)"]]))),
    unname(expected_cor),
    tolerance = 1e-10
  )

  # The fitted-scale K = 2 correlation follows the LKJ(eta) marginal:
  # (r + 1) / 2 ~ Beta(eta, eta). The prior is exact, so with a fixed seed the
  # one-sample KS p-value is a single uniform draw; alpha = 0.01.
  expect_gt(
    stats::ks.test((rho + 1) / 2, "pbeta", eta, eta)$p.value,
    0.01
  )

  # The original-scale prior draws unscale the SD monitors and the
  # correlation together, matching the catalog quantities.
  transformed <- transform_prior_samples(fit, n_samples = n, seed = 58L)
  for(k in 1:2){
    expect_equal(
      unname(transformed[, sd_names[[k]]]),
      as.numeric(as.matrix(draws[[paste0("(mu) sd(", c("intercept", "x")[[k]], ")")]])),
      tolerance = 1e-12
    )
  }
  expect_equal(
    unname(transformed[, R_names[[2L]]]),
    unname(expected_cor),
    tolerance = 1e-10
  )
})

test_that("prior draws carry SDs of child component allocations", {

  data <- data.frame(
    x = c(1, 4, 6, 2, 8, 5, 3, 7),
    s = factor(rep(c("a", "b", "c", "d"), each = 2L)),
    d = factor(rep(c("u", "v"), 4L))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + us(1 + x | s) +
      random(1 | d, name = "d", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(s = "s", d = "d"),
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      random_variance_allocation(
        name = "s_components",
        parent = allocation_ref("total_re", "s"),
        terms = "s",
        target = "sd_component",
        scale = "mean_variance",
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  terms <- formula_result$formula_design$random_effects
  s_term <- terms[[which(vapply(terms, `[[`, character(1), "block_name") == "s")]]
  d_term <- terms[[which(vapply(terms, `[[`, character(1), "block_name") == "d")]]
  source_name <- "mu__xRE_ALLOCx_total_re__allocation_sd"
  total_weights <- paste0("mu__xRE_ALLOCx_total_re__weight[", 1:2, "]")
  child_weights <- paste0("mu__xRE_ALLOCx_s_components__weight[", 1:2, "]")
  columns <- c(
    "mu_intercept", "mu_x", source_name, total_weights, child_weights,
    s_term$sd_parameter_names,
    .prior_monitor_matrix_names(s_term$correlation$cholesky_name, 2L),
    .prior_monitor_matrix_names(s_term$correlation$correlation_name, 2L),
    s_term$correlation$primitive_names,
    d_term$sd_parameter_names
  )
  fit <- .prior_monitor_test_fit(formula_result, columns)

  raw <- transform_prior_samples(fit, n_samples = 2000L, seed = 59L,
                                 formula_scale = list())
  # s: sd * sqrt(w_total[1]) * sqrt(2 * w_child[k]); d: sd * sqrt(w_total[2]).
  s_source <- raw[, source_name] * sqrt(raw[, total_weights[[1L]]])
  for(k in 1:2){
    expect_identical(
      unname(raw[, s_term$sd_parameter_names[[k]]]),
      unname(s_source * sqrt(2 * raw[, child_weights[[k]]]))
    )
  }
  expect_identical(
    unname(raw[, d_term$sd_parameter_names]),
    unname(raw[, source_name] * sqrt(raw[, total_weights[[2L]]]))
  )

  draws <- .prior_monitor_catalog_draws(fit, raw)
  expect_true("(mu) s: cor(intercept,x)" %in% names(draws))
  for(name in names(draws)){
    expect_true(all(is.finite(as.numeric(as.matrix(draws[[name]])))),
                info = name)
  }
})

test_that("prior draws carry transformed scalar correlations and LKJ partial correlations", {

  data <- data.frame(
    x = c(1, 4, 6, 2, 8, 5, 3, 7, 2, 6),
    t = factor(rep(c("t1", "t2", "t3", "t4", "t5"), 2L)),
    g = factor(rep(c("a", "b"), each = 5L))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(0, Inf))

  # HCS with SD-component allocation (indexed SD monitors) and a Fisher-z
  # correlation: rho <- tanh(rho_z).
  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(t | g),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "comp",
        terms = "g",
        target = "sd_component",
        scale = "mean_variance",
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = rep(1, 5L)))
      ),
      g = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  hcs_term <- hcs_result$formula_design$random_effects[[1L]]
  hcs_rho <- hcs_term$correlation
  expect_identical(hcs_rho$rho_scale, "fisher_z")
  weight_names <- paste0("mu__xRE_ALLOCx_comp__weight[", 1:5, "]")
  hcs_fit <- .prior_monitor_test_fit(hcs_result, c(
    "mu_intercept", "mu__xRE_ALLOCx_comp__allocation_sd", weight_names,
    hcs_rho$sample_name, hcs_term$sd_parameter_names, hcs_rho$rho_name
  ))
  hcs_raw <- transform_prior_samples(hcs_fit, n_samples = 2000L, seed = 60L)
  for(k in 1:5){
    expect_identical(
      unname(hcs_raw[, hcs_term$sd_parameter_names[[k]]]),
      unname(hcs_raw[, "mu__xRE_ALLOCx_comp__allocation_sd"] *
               sqrt(5 * hcs_raw[, weight_names[[k]]]))
    )
  }
  expect_identical(
    unname(hcs_raw[, hcs_rho$rho_name]),
    unname(tanh(hcs_raw[, hcs_rho$sample_name]))
  )
  hcs_draws <- .prior_monitor_catalog_draws(hcs_fit, hcs_raw)
  expect_identical(
    as.numeric(as.matrix(hcs_draws[["(mu) cor"]])),
    unname(hcs_raw[, hcs_rho$rho_name])
  )
  for(name in names(hcs_draws)){
    expect_true(all(is.finite(as.numeric(as.matrix(hcs_draws[[name]])))),
                info = name)
  }

  # AR1 with a logit-scale correlation: rho <- -1 + 2 * ilogit(rho_logit).
  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(t | g),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("normal", list(0, 1)),
          cor_scale = "logit"
        )
      )
    )
  )
  ar1_term <- ar1_result$formula_design$random_effects[[1L]]
  ar1_rho <- ar1_term$correlation
  ar1_fit <- .prior_monitor_test_fit(ar1_result, c(
    "mu_intercept", unique(ar1_term$sd_parameter_names),
    ar1_rho$sample_name, ar1_rho$rho_name
  ))
  ar1_raw <- transform_prior_samples(ar1_fit, n_samples = 2000L, seed = 61L)
  expect_identical(
    unname(ar1_raw[, ar1_rho$rho_name]),
    unname(-1 + 2 * stats::plogis(ar1_raw[, ar1_rho$sample_name]))
  )
  ar1_draws <- .prior_monitor_catalog_draws(ar1_fit, ar1_raw)
  expect_identical(
    as.numeric(as.matrix(ar1_draws[["(mu) cor"]])),
    unname(ar1_raw[, ar1_rho$rho_name])
  )

  # Monitored LKJ partial correlations: cpc[p] <- 2 * u[p] - 1.
  lkj_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | g),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 2, include_primitives = TRUE)
      )
    )
  )
  lkj_correlation <- lkj_result$formula_design$random_effects[[1L]]$correlation
  expect_length(lkj_correlation$cpc_names, 1L)
  lkj_fit <- .prior_monitor_test_fit(lkj_result, c(
    "mu_intercept",
    lkj_result$formula_design$random_effects[[1L]]$sd_parameter_names,
    .prior_monitor_matrix_names(lkj_correlation$cholesky_name, 2L),
    .prior_monitor_matrix_names(lkj_correlation$correlation_name, 2L),
    lkj_correlation$primitive_names,
    lkj_correlation$cpc_names
  ))
  lkj_raw <- transform_prior_samples(lkj_fit, n_samples = 2000L, seed = 62L)
  expect_identical(
    unname(lkj_raw[, lkj_correlation$cpc_names]),
    unname(2 * lkj_raw[, lkj_correlation$primitive_names] - 1)
  )
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

  # every alias is the table label of its label parts
  expect_identical(parameter_labels(catalog$aliases, style = "table"), catalog$aliases$alias)
  # the pairwise alias renders under a caller vocabulary (cor -> rho), and its
  # rendering under the catalog's names selects the shared correlation
  pair <- catalog$aliases[
    catalog$aliases$alias == "(mu) cor(outcome[sensitivity],outcome[specificity])", ,
    drop = FALSE
  ]
  expect_identical(nrow(pair), 1L)
  expect_identical(pair$quantity_id, rho$quantity_id)
  expect_identical(pair$label_parts[[1L]]$random$arguments,
                   c("outcome[sensitivity]", "outcome[specificity]"))
  vocabulary <- c(sd = "tau", cor = "rho")
  expect_identical(
    parameter_labels(pair, style = "table", vocabulary = vocabulary),
    "(mu) rho(outcome[sensitivity],outcome[specificity])"
  )
  expect_identical(
    parameter_labels(pair$label_parts, style = "table", formula_prefix = FALSE,
                     vocabulary = vocabulary),
    "rho(outcome[sensitivity],outcome[specificity])"
  )
  expect_identical(
    parameter_catalog_resolve(
      catalog, parameter_labels(pair, style = "table"), "mu"
    )$quantity_id,
    rho$quantity_id
  )
  # the vocabulary renames random-effect quantities only
  sd_aliases <- catalog$aliases[grepl("(^|\\) )sd\\(", catalog$aliases$alias), , drop = FALSE]
  expect_true(nrow(sd_aliases) > 0L)
  expect_identical(
    parameter_labels(sd_aliases, style = "table", vocabulary = vocabulary),
    sub("(^|\\) )sd\\(", "\\1tau(", sd_aliases$alias)
  )
  intercept <- catalog$aliases[catalog$aliases$alias == "mu_intercept", , drop = FALSE]
  expect_identical(parameter_labels(intercept, style = "table", vocabulary = vocabulary),
                   "mu_intercept")
  expect_error(parameter_labels(pair, vocabulary = c("rho")), "'vocabulary' must be a named")
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

  testthat::local_mocked_bindings(
    .bt_random_effect_summary_sd_samples = function(...){
      stop("Invalid random-SD scaling metadata.", call. = FALSE)
    },
    .package = "BayesTools"
  )
  expect_error(
    .bt_build_parameter_catalog(
      coordinates = parameter_coordinates(fit),
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design),
      formula_scale = list(mu = formula_result$formula_scale)
    ),
    "Invalid random-SD scaling metadata.",
    fixed = TRUE
  )
})

test_that("variances of one-to-one scaled SDs are one-to-one transforms", {

  data <- data.frame(
    x = c(1, 2, 4, 7),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    ~ 1 + (0 + x | id),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    )
  )
  sd_name <- formula_result$formula_design$random_effects[[1L]]$sd_parameter_names
  source_values <- c(0.5, 1.5)
  values <- cbind(mu_intercept = 0, source_values)
  colnames(values)[2L] <- sd_name
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design),
    formula_scale = list(mu = formula_result$formula_scale)
  )
  catalog <- parameter_catalog(fit)
  scale <- 1 / stats::sd(data$x)
  sd_selection <- parameter_catalog_resolve(catalog, "(mu) sd(x)", "mu")
  var_selection <- parameter_catalog_resolve(catalog, "(mu) var(x)", "mu")
  sd_key <- sd_selection$quantities$extraction_key[[1L]]
  var_key <- var_selection$quantities$extraction_key[[1L]]

  expect_identical(var_selection$quantities$source_type, "one_to_one_transform")
  expect_identical(var_key$source_type, "one_to_one_transform")
  expect_identical(var_key$source_parameter, sd_name)
  expect_identical(var_key$source_parameter, sd_key$source_parameter)
  expect_identical(var_key$source_transform, "random_var")
  expect_equal(var_key$source_scale, scale)
  expect_identical(var_key$source_scale, sd_key$source_scale)

  var_draws <- as.numeric(parameter_draws(fit, var_selection)[[1L]][, 1L])
  expect_equal(var_draws, (scale * source_values)^2, tolerance = 1e-12)
  transform <- parameter_transform(fit, var_selection)
  expect_identical(transform$type, "square")
  expect_equal(transform$scale, scale)
  expect_equal(
    parameter_transform_forward(source_values, transform),
    var_draws,
    tolerance = 1e-12
  )
  expect_equal(
    parameter_transform_inverse(var_draws, transform),
    source_values,
    tolerance = 1e-12
  )
  step <- 1e-6
  expect_equal(
    parameter_transform_jacobian(source_values, transform),
    (parameter_transform_forward(source_values + step, transform) -
       parameter_transform_forward(source_values - step, transform)) /
      (2 * step),
    tolerance = 1e-8
  )
  # An unscaled square keeps its meaning.
  expect_identical(
    parameter_transform_forward(3, list(type = "square")),
    9
  )
  expect_error(
    parameter_transform_forward(3, list(type = "square", scale = -1)),
    "Unsupported semantic parameter transform."
  )

  # Prior of (s * sd)^2 for sd ~ Normal+(0, 1): dnorm(sqrt(v), 0, s) / sqrt(v).
  # Linear interpolation on the default 4096-knot grid is accurate to about
  # 1e-3 (relative) at these ordinates, as for the unscaled variance.
  density <- parameter_prior_density(fit, var_selection)
  expect_s3_class(density, "prior_linear_density")
  ordinates <- c(0.05, 0.2, 0.5)
  expect_equal(
    stats::approx(density$density$x, density$density$y, ordinates)$y,
    stats::dnorm(sqrt(ordinates), 0, scale) / sqrt(ordinates),
    tolerance = 5e-3
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

test_that("one-to-one random summaries take status from their own source", {

  set.seed(1)
  data <- data.frame(id = factor(rep(c("a", "b", "c"), each = 4)),
                     x = stats::rnorm(12))
  formula_result <- JAGS_formula(
    ~ 1 + x + random(1 + x | id, name = "id", covariance = "diag"),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1)),
         x = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(
      sd = prior("normal", list(0, 1), list(0, Inf)),
      terms = list(intercept = prior("point", list(location = 0.5)))
    ))
  )
  formula_design <- list(mu = formula_result$formula_design)
  random_term <- formula_result$formula_design$random_effects[[1L]]
  sd_names <- unique(random_term$sd_parameter_names)
  coordinates <- .bt_build_parameter_coordinates(
    columns = c("mu_intercept", "mu_x", "mu__xREx__id_x"),
    prior_list = formula_result$prior_list,
    formula_design = formula_design
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates,
    formula_result$prior_list,
    formula_design
  )
  quantity <- function(label){
    parameter_catalog_resolve(catalog, label, namespace = "mu")$quantities
  }

  fixed_sd <- quantity("(mu) sd(intercept)")
  expect_identical(fixed_sd$status, "structural")
  expect_identical(fixed_sd$fixed_value, 0.5)
  expect_identical(fixed_sd$source_type, "identity")
  fixed_var <- quantity("(mu) var(intercept)")
  expect_identical(fixed_var$status, "structural")
  expect_identical(fixed_var$fixed_value, 0.25)
  expect_identical(fixed_var$source_type, "one_to_one_transform")
  for(label in c("(mu) sd(x)", "(mu) var(x)")){
    expect_identical(quantity(label)$status, "sampled", info = label)
    expect_identical(quantity(label)$fixed_value, NA_real_, info = label)
  }
  # The block evaluator still declares every SD coordinate of the block.
  for(label in c("(mu) sd(intercept)", "(mu) var(intercept)",
                 "(mu) sd(x)", "(mu) var(x)")){
    expect_setequal(
      quantity(label)$extraction_key[[1L]]$dependencies,
      sd_names
    )
  }

  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(cbind(
      mu_intercept = c(0.1, 0.2),
      mu_x = c(0.3, 0.4),
      mu__xREx__id_x = c(1.5, 2)
    ))),
    prior_list = formula_result$prior_list,
    formula_design = formula_design
  )
  fit_catalog <- parameter_catalog(fit)
  draws <- function(label){
    as.numeric(parameter_draws(
      fit,
      parameter_catalog_resolve(fit_catalog, label, namespace = "mu")
    )[[1L]][, 1L])
  }
  expect_identical(draws("(mu) sd(intercept)"), c(0.5, 0.5))
  expect_identical(draws("(mu) var(intercept)"), c(0.25, 0.25))
  expect_identical(draws("(mu) sd(x)"), c(1.5, 2))
  expect_identical(draws("(mu) var(x)"), c(2.25, 4))
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

test_that("allocation parent links are scoped by formula parameter", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  build <- function(parameter, sd = NULL, sd_source = NULL){
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = parameter,
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        allocation = random_variance_allocation(
          name = "allocation",
          sd = sd,
          sd_source = sd_source,
          weights = prior("dirichlet", list(alpha = c(2, 3)))
        )
      )
    )
  }
  columns <- function(parameter, scale = TRUE){
    c(
      paste0(parameter, "_intercept"),
      if(scale) paste0(parameter, "__xRE_ALLOCx_allocation__allocation_sd"),
      paste0(parameter, "__xRE_ALLOCx_allocation__weight[", 1:2, "]")
    )
  }
  catalog_for <- function(mu, sigma, sigma_scale){
    prior_list <- c(mu$prior_list, sigma$prior_list)
    formula_design <- list(
      mu = mu$formula_design,
      sigma = sigma$formula_design
    )
    coordinates <- .bt_build_parameter_coordinates(
      columns = c(columns("mu"), columns("sigma", sigma_scale)),
      prior_list = prior_list,
      formula_design = formula_design
    )
    .bt_build_parameter_catalog(coordinates, prior_list, formula_design)
  }
  parent_of <- function(catalog, label, namespace){
    parameter_catalog_resolve(catalog, label, namespace)$quantities$parent_quantity_id
  }
  id_of <- function(catalog, label, namespace){
    parameter_catalog_resolve(catalog, label, namespace)$quantity_id
  }
  children <- c(
    "allocation: var_total",
    "allocation: var_prop(study)",
    "allocation: var_prop(drug)"
  )

  # Same allocation label in both formulas: each links to its own total.
  catalog <- catalog_for(
    build("mu", sd = prior("gamma", list(2, 2))),
    build("sigma", sd = prior("gamma", list(3, 3))),
    sigma_scale = TRUE
  )
  for(namespace in c("mu", "sigma")){
    for(label in children){
      expect_identical(
        parent_of(catalog, label, namespace),
        id_of(catalog, "allocation: sd_total", namespace),
        info = paste(namespace, label)
      )
    }
  }

  # A formula without an allocation total never links to another formula's.
  catalog <- catalog_for(
    build("mu", sd = prior("gamma", list(2, 2))),
    build("sigma", sd_source = random_sd_source("tau", shape = "row")),
    sigma_scale = FALSE
  )
  for(label in children[-1L]){
    expect_identical(parent_of(catalog, label, "sigma"), "", info = label)
    expect_identical(
      parent_of(catalog, label, "mu"),
      id_of(catalog, "allocation: sd_total", "mu"),
      info = label
    )
  }
})

test_that("coordinates and factor levels expose their prior densities", {

  # Treatment levels are their N(0, 1) coefficients (the reference level is
  # structurally 0); mean-difference levels are sum_j w_j beta_j of iid
  # N(0, 0.5) coordinates, i.e. N(0, 0.5 * sqrt(sum(w^2))). A coefficient's
  # 'multiply_by' scale enters the linear predictor, not its prior density.
  data <- data.frame(t = factor(c("a", "b", "c", "a")), m = factor(c("x", "y", "z", "z")),
                     x = c(-1, 0, 1, .5))
  formula_result <- JAGS_formula(~ t + m + x, "mu", data = data, prior_list = list(
    intercept = prior("normal", list(0, 2)),
    t = prior_factor("normal", list(0, 1), contrast = "treatment"),
    m = prior_factor("mnormal", list(0, .5), contrast = "meandif"),
    x = prior("normal", list(0, .3))
  ))
  prior_list <- formula_result$prior_list
  prior_list$sigma <- prior("normal", list(0, 1), list(0, Inf))
  attr(prior_list$mu_x, "multiply_by") <- "sigma"
  columns <- c("mu_intercept", "mu_t[1]", "mu_t[2]", "mu_m[1]", "mu_m[2]", "mu_x", "sigma")
  values <- matrix(.5, nrow = 2L, ncol = length(columns), dimnames = list(NULL, columns))
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- parameter_catalog(fit)
  density_at <- function(name, value){
    density <- parameter_prior_density(fit, parameter_catalog_resolve(catalog, name, "mu"))
    expect_s3_class(density, "prior_linear_density")
    ordinate <- prior_density_ordinate(density, value)
    expect_true(ordinate$exact)
    ordinate
  }
  for(value in c(-.4, .3)){
    expect_equal(exp(density_at("mu_intercept", value)$log_density), stats::dnorm(value, 0, 2),
                 tolerance = 1e-14)
    expect_equal(exp(density_at("mu_x", value)$log_density), stats::dnorm(value, 0, .3),
                 tolerance = 1e-14)
    expect_equal(exp(density_at("mu_t[b]", value)$log_density), stats::dnorm(value),
                 tolerance = 1e-14)
    for(level in c("x", "y", "z")){
      key <- parameter_catalog_resolve(catalog, paste0("mu_m[", level, "]"), "mu")$quantities$extraction_key[[1L]]
      expect_equal(exp(density_at(paste0("mu_m[", level, "]"), value)$log_density),
                   stats::dnorm(value, 0, .5 * sqrt(sum(key$weights^2))), tolerance = 1e-12)
    }
  }
  reference <- density_at("mu_t[a]", 0)
  expect_identical(reference$behavior, "point_mass")
  expect_equal(reference$point_mass, 1)
})

test_that("coordinate prior densities are NULL only without an owning prior", {

  columns <- c("mu", "tau")
  values <- matrix(.5, nrow = 2L, ncol = length(columns), dimnames = list(NULL, columns))
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(values)),
    prior_list = list(
      mu  = prior("normal", list(0, 1)),
      tau = prior("normal", list(0, 1), list(0, Inf))
    )
  )
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "mu")
  with_mu_prior <- function(value){
    prior_list <- attr(fit, "prior_list", exact = TRUE)
    prior_list["mu"] <- list(value)
    attr(fit, "prior_list") <- prior_list
    fit
  }

  # a coordinate that no fitted prior owns has no prior density by rule
  unowned <- fit
  attr(unowned, "prior_list") <- attr(fit, "prior_list", exact = TRUE)["tau"]
  expect_null(parameter_prior_density(unowned, selection))

  # an owning entry that is not a BayesTools prior stops (the failed context
  # build was previously reported as an unavailable density)
  for(unsupported in list(
    structure(list(distribution = "normal"), class = "not_a_prior"),
    "normal",
    list(1, 2)
  )){
    expect_error(
      parameter_prior_density(with_mu_prior(unsupported), selection),
      "The prior density of 'mu' is unavailable: the prior distribution of 'mu' is not a BayesTools prior.",
      fixed = TRUE
    )
  }

  # model lists of BayesTools priors remain supported owners
  density <- parameter_prior_density(
    with_mu_prior(list(prior("normal", list(0, 1)), prior("normal", list(0, 2)))),
    selection
  )
  expect_s3_class(density, "prior_linear_density")
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
  # The transformed densities record their deterministic provenance, so
  # heights are exact: var_mult = 2 w (affine) and sd_mult = sqrt(2 w)
  # (exp_lin) with w ~ Beta(2, 3).
  for(density in list(multiplier_density, sd_mult_density, common_variance_density)){
    expect_false(is.null(attr(density, "adaptive_evaluation", exact = TRUE)))
  }
  expect_equal(
    .prior_linear_density_height(multiplier_density, .5),
    stats::dbeta(.25, 2, 3) / 2,
    tolerance = 1e-12
  )
  expect_equal(
    .prior_linear_density_height(sd_mult_density, 1),
    stats::dbeta(.5, 2, 3),
    tolerance = 1e-12
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
  proportion <- parameter_catalog_resolve(
    parameter_catalog(fit), "var_prop(study)", "mu"
  )
  proportion_density <- parameter_prior_density(fit, proportion, n_grid = 512L)
  expect_s3_class(proportion_density, "prior_linear_density")
  expect_equal(prior_density_ordinate(proportion_density, .4)$log_density,
               stats::dbeta(.4, 2, 3, log = TRUE))
})

test_that("sparse allocation margins keep boundary-singular prior densities", {

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
        weights = prior("dirichlet", list(alpha = c(.5, 1.5)))
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
  catalog <- parameter_catalog(fit)

  # var_prop(study) is the Beta(0.5, 1.5) margin, infinite at zero.
  proportion <- parameter_catalog_resolve(catalog, "var_prop(study)", "mu")
  proportion_density <- parameter_prior_density(fit, proportion)
  x  <- proportion_density$density$x
  y  <- proportion_density$density$y
  dx <- x[2] - x[1]
  expect_true(all(is.finite(y)))
  # Exact CDF cell mass at the singular bound; the tolerance covers the grid
  # renormalisation of the O(sqrt(dx)) neighbouring-cell midpoint error.
  expect_equal(y[1] * dx, stats::pbeta(dx / 2, .5, 1.5), tolerance = 5e-3)
  expect_equal(
    .prior_linear_density_height(proportion_density, .4),
    stats::dbeta(.4, .5, 1.5)
  )
  expect_identical(prior_density_ordinate(proportion_density, 0)$behavior, "infinite")

  # The component SD is total SD * sqrt(var_prop): E[sd^2] = 1 * .5 / 2. Its
  # exact density (the scale prior times the square root of the Beta(0.5, 1.5)
  # share) integrates to it: the second moment is P(sd^2 > t) integrated over
  # t, from the exact region probabilities (the half-normal scale prior leaves
  # less than 1e-11 beyond t = 50). The integral runs over u = t^(1/4): P(sd^2 > t)
  # leaves 1 at t = 0 with an infinite slope (the density of sd is infinite at
  # zero), which the factor dt = 4 u^3 du flattens, so the quadrature needs 63
  # of the 399 evaluations of the integral over t.
  component <- parameter_catalog_resolve(catalog, "study: sd(intercept)", "mu")
  component_density <- parameter_prior_density(fit, component)
  expect_identical(attr(component_density, "adaptive_evaluation")$kind, "allocation_product")
  expect_identical(prior_density_ordinate(component_density, 0)$behavior, "infinite")
  second_moment <- stats::integrate(Vectorize(function(u){
    region <- list(intervals = matrix(c(u^2, Inf), 1L), indicator = function(x) x > u^2)
    4 * u^3 * as.numeric(.prior_linear_density_region_probability(component_density, region))
  }), 0, 50^(1 / 4), rel.tol = 1e-8)$value
  expect_equal(second_moment, .25, tolerance = 1e-6)
})


test_that("shared-gate proportions use their declared conditional Dirichlet prior", {

  make_fit <- function(independent = FALSE, probability = .5,
                       scale_prior = prior("gamma", list(2, 2))){

    split <- random_variance_allocation(
      name = "split", terms = c(study = "study", esid = "esid"),
      parent = if(!independent) allocation_ref("root", "component") else NULL,
      sd = if(independent) scale_prior else NULL,
      weights = prior("dirichlet", list(alpha = c(2, 3))),
      inclusion = if(independent) list(
        study = prior("spike", list(location = probability)),
        esid = prior("spike", list(location = probability))
      ) else NULL
    )
    allocations <- if(independent) list(split) else list(
      random_variance_allocation(
        name = "root", terms = c(component = "both"),
        sd = scale_prior,
        inclusion = list(component = prior("spike", list(location = probability)))
      ),
      split
    )
    result <- JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | esid, name = "esid", covariance = "diag"),
      parameter = "mu",
      data = data.frame(study = factor(c("a", "a", "b", "b")),
                        esid = factor(seq_len(4L))),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(allocation = allocations)
    )
    allocation <- result$formula_design$random_allocations$split
    gates <- .bt_random_effect_summary_allocation_gate_names(allocation)
    columns <- c("mu_intercept", allocation$source$name,
                  paste0(allocation$weight_name, "[", 1:2, "]"), gates)
    samples <- matrix(.5, nrow = 2L, ncol = length(columns),
                       dimnames = list(NULL, columns))
    if(is.prior.point(scale_prior)){
      samples[, allocation$source$name] <- scale_prior$parameters$location
    }
    samples[, gates] <- if(probability %in% c(0, 1)) probability else c(0, 1)
    .parameter_catalog_test_fit(
      coda::mcmc.list(coda::mcmc(samples)),
      prior_list = result$prior_list,
      formula_design = list(mu = result$formula_design)
    )
  }

  fit <- make_fit()
  points <- c(0, .1, .5, .9, 1)
  set.seed(1)
  before <- .Random.seed
  for(i in 1:2){
    selection <- parameter_catalog_resolve(
      parameter_catalog(fit),
      paste0("split: var_prop(", c("study", "esid")[[i]], ")"), "mu"
    )
    expect_null(parameter_transform(fit, selection))
    density <- parameter_prior_density(fit, selection, n_grid = 512L)
    expect_s3_class(density, "prior_linear_density")
    ordinates <- vapply(points, function(value){
      prior_density_ordinate(density, value)$log_density
    }, numeric(1))
    expect_equal(ordinates, stats::dbeta(points, c(2, 3)[[i]],
                                       c(3, 2)[[i]], log = TRUE))
  }
  expect_identical(.Random.seed, before)

  fixed_positive <- make_fit(probability = 1,
                             scale_prior = prior("point", list(location = .7)))
  selection <- parameter_catalog_resolve(
    parameter_catalog(fixed_positive), "split: var_prop(study)", "mu"
  )
  density <- parameter_prior_density(fixed_positive, selection)
  expect_equal(prior_density_ordinate(density, .4)$log_density,
               stats::dbeta(.4, 2, 3, log = TRUE))

  independent <- make_fit(independent = TRUE)
  independent_points <- c(.1, .4, .9)
  for(i in 1:2){
    selection <- parameter_catalog_resolve(
      parameter_catalog(independent),
      paste0("split: var_prop(", c("study", "esid")[[i]], ")"), "mu"
    )
    density <- parameter_prior_density(independent, selection, n_grid = 512L)
    expect_s3_class(density, "prior_linear_density")
    expect_equal(
      .prior_linear_density_point_mass(density, 0),
      1 / 3,
      tolerance = 1e-8
    )
    expect_equal(
      .prior_linear_density_point_mass(density, 1),
      1 / 3,
      tolerance = 1e-8
    )
    alpha_i <- c(2, 3)[[i]]
    beta_i <- c(3, 2)[[i]]
    ordinates <- vapply(independent_points, function(value){
      prior_density_ordinate(density, value)$log_density
    }, numeric(1))
    expect_equal(
      ordinates,
      log((1 / 3) * stats::dbeta(independent_points, alpha_i, beta_i)),
      tolerance = 1e-3
    )
    # the mixture of atoms and Beta components has an exact ordinate
    expect_true(all(vapply(independent_points, function(value){
      prior_density_ordinate(density, value)$exact
    }, logical(1))))
    expect_equal(
      vapply(independent_points, function(value){
        .prior_linear_density_height(density, value)
      }, numeric(1)),
      (1 / 3) * stats::dbeta(independent_points, alpha_i, beta_i),
      tolerance = 1e-12
    )
  }

  # Independently gated totals and component SDs: the realized total SD is
  # sd * sqrt(g1 w + g2 (1 - w)) and the study SD sd * sqrt(w) * g1, with
  # sd ~ gamma(2, 2), w ~ Beta(2, 3) and gates g ~ Bernoulli(1/2). The atoms
  # at 0 (all gates off: 1/4; study gate off: 1/2) are exact, and so are the
  # continuous parts (the allocation product measure): their distribution
  # functions (exact region probabilities) are checked against the Monte
  # Carlo distribution function within 4 MC SD.
  set.seed(2)
  n <- 2e5
  scale_draws <- stats::rgamma(n, 2, 2)
  share_draws <- stats::rbeta(n, 2, 3)
  gate_1 <- stats::rbinom(n, 1, .5)
  gate_2 <- stats::rbinom(n, 1, .5)
  total_sd <- scale_draws * sqrt(gate_1 * share_draws + gate_2 * (1 - share_draws))
  study_sd <- scale_draws * sqrt(share_draws) * gate_1
  references <- list(
    "split: sd_total" = list(draws = total_sd, zero = .25),
    "split: var_total" = list(draws = total_sd^2, zero = .25),
    "study: sd(intercept)" = list(draws = study_sd, zero = .5),
    "study: var(intercept)" = list(draws = study_sd^2, zero = .5)
  )
  for(name in names(references)){
    density <- parameter_prior_density(
      independent, parameter_catalog_resolve(parameter_catalog(independent), name, "mu")
    )
    expect_s3_class(density, "prior_linear_density")
    expect_equal(.prior_linear_density_point_mass(density, 0), references[[name]]$zero,
                 tolerance = 1e-12, info = name)
    expect_identical(attr(density, "adaptive_evaluation", exact = TRUE)$kind,
                     "allocation_product")
    for(value in c(.1, .3, .6, 1, 2)){
      reference <- mean(references[[name]]$draws <= value)
      region <- list(intervals = matrix(c(-Inf, value), 1L), indicator = function(x) x <= value)
      expect_lte(abs(as.numeric(.prior_linear_density_region_probability(density, region)) - reference),
                 4 * sqrt(reference * (1 - reference) / n))
      expect_true(prior_density_ordinate(density, value)$exact)
    }
  }

  # a model-averaged (mixture) scale prior keeps the gated proportion density
  mixture_scale <- make_fit(independent = TRUE, scale_prior = prior_mixture(list(
    prior("gamma", list(2, 2), prior_weights = 3),
    prior("gamma", list(3, 1), prior_weights = 2)
  ), is_null = c(FALSE, FALSE)))
  selection <- parameter_catalog_resolve(
    parameter_catalog(mixture_scale), "split: var_prop(study)", "mu"
  )
  density <- parameter_prior_density(mixture_scale, selection)
  expect_equal(.prior_linear_density_point_mass(density, 0), 1 / 3, tolerance = 1e-12)
  expect_equal(.prior_linear_density_point_mass(density, 1), 1 / 3, tolerance = 1e-12)
  expect_equal(.prior_linear_density_height(density, .4), stats::dbeta(.4, 2, 3) / 3,
               tolerance = 1e-12)

  independent_fixed <- make_fit(independent = TRUE, probability = 1)
  selection <- parameter_catalog_resolve(
    parameter_catalog(independent_fixed), "split: var_prop(study)", "mu"
  )
  density <- parameter_prior_density(independent_fixed, selection)
  expect_equal(prior_density_ordinate(density, .4)$log_density,
               stats::dbeta(.4, 2, 3, log = TRUE))
  expect_equal(.prior_linear_density_point_mass(density, 0), 0)
  expect_equal(.prior_linear_density_point_mass(density, 1), 0)

  for(unavailable in list(make_fit(probability = 0),
                          make_fit(independent = TRUE, probability = 0),
                          make_fit(scale_prior = prior("point", list(location = 0))))){
    selection <- parameter_catalog_resolve(
      parameter_catalog(unavailable), "split: var_prop(study)", "mu"
    )
    expect_null(parameter_prior_density(unavailable, selection))
  }
})

test_that("independently gated proportion mixtures fold fixed inclusion gates", {

  # One free gate, one always-on component, one always-off component:
  # the Beta second parameter is the always-on alpha only.
  density <- BayesTools:::.bt_parameter_prior_density_gated_var_prop_mixture(
    alpha = c(2, 3, 4),
    index = 1L,
    probability = c(0.5, 1, 0),
    n_grid = 256L,
    tail_prob = 1e-4
  )
  expect_s3_class(density, "prior_linear_density")
  expect_equal(.prior_linear_density_point_mass(density, 0), 0.5, tolerance = 1e-8)
  expect_equal(.prior_linear_density_point_mass(density, 1), 0)
  expect_equal(
    prior_density_ordinate(density, 0.4)$log_density,
    log(0.5 * stats::dbeta(0.4, 2, 3)),
    tolerance = 1e-3
  )
  expect_equal(
    .prior_linear_density_height(density, 0.4),
    0.5 * stats::dbeta(0.4, 2, 3),
    tolerance = 1e-12
  )

  too_many <- rep(0.5, 22L)
  expect_error(
    BayesTools:::.bt_parameter_prior_density_gated_var_prop_mixture(
      alpha = rep(1, 22L),
      index = 1L,
      probability = too_many,
      n_grid = 32L,
      tail_prob = 1e-4
    ),
    "Independently gated variance-proportion prior density is unavailable for more than 20 free inclusion gates.",
    fixed = TRUE
  )
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

test_that("coordinate quantities declare the support of their priors", {

  prior_list <- list(
    theta = prior("normal", list(0, 1), list(0, Inf)),
    phi   = prior_spike_and_slab(prior("normal", list(0, 1)))
  )
  coordinates <- .bt_build_parameter_coordinates(
    columns = c("theta", "phi", "extra"),
    prior_list = prior_list
  )
  catalog <- .bt_build_parameter_catalog(coordinates, prior_list = prior_list)
  quantities <- catalog$quantities
  support <- stats::setNames(quantities$support, quantities$canonical_name)
  expect_equal(support$theta$bounds, c(0, Inf))
  expect_true(support$theta$exact)
  # a spike-and-slab coefficient: the slab's interval and the spike's point
  expect_equal(support$phi$bounds, c(-Inf, Inf))
  expect_identical(support$phi$points, 0)
  # a coordinate without an owning prior has no declared support
  expect_null(support$extra)
  expect_true(all(quantities$definedness == "always"))

  # provider rows must declare both columns
  provider_rows <- quantities[1L, , drop = FALSE]
  provider_rows$definedness <- NA_character_
  expect_error(
    BayesTools:::.bt_validate_parameter_catalog_tables(provider_rows, catalog$aliases[0, ]),
    "malformed field types or missing metadata"
  )
})

test_that("malformed catalogs and stale selections fail closed", {

  coordinates <- .bt_build_parameter_coordinates(columns = "theta")
  catalog <- .bt_build_parameter_catalog(coordinates)
  selection <- parameter_catalog_resolve(catalog, "theta")

  broken <- catalog
  broken$schema_version <- BayesTools:::.bt_parameter_map_version + 1L
  expect_error(
    .bt_validate_parameter_catalog(broken),
    "Refit or rebuild",
    class = "BayesTools_refit_required"
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

  duplicated_name <- catalog
  duplicated_name$quantities <- rbind(
    duplicated_name$quantities,
    duplicated_name$quantities
  )
  duplicated_name$quantities$quantity_id[2L] <- "BayesTools::duplicate"
  # identical in canonical name, namespace, and component: the resolver could
  # not tell these apart, so construction still refuses them
  expect_error(
    .bt_validate_parameter_catalog(duplicated_name),
    "cannot be resolved"
  )
})

test_that("gated totals include the all-off zero and leave var_prop undefined", {

  data <- data.frame(
    study = factor(c("a", "a", "b", "b")),
    esid = factor(seq_len(4L))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "split",
        terms = c(study = "study", esid = "esid"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3))),
        inclusion = list(
          study = prior("spike", list(location = 0.5)),
          esid = prior("spike", list(location = 0.5))
        )
      )
    )
  )
  allocation <- formula_result$formula_design$random_allocations[[1L]]
  gates <- .bt_random_effect_summary_allocation_gate_names(allocation)
  columns <- c(
    "mu_intercept",
    allocation$source$name,
    paste0(allocation$weight_name, "[", 1:2, "]"),
    gates
  )
  samples <- matrix(
    c(
      0, 2, 0.4, 0.6, 0, 0,
      0, 2, 0.4, 0.6, 1, 1
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, columns)
  )
  fit <- .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(samples)),
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design)
  )
  catalog <- parameter_catalog(fit)
  coordinates <- parameter_coordinates(fit)

  expect_true(all(coordinates$internal[coordinates$coordinate_name %in% gates]))
  expect_false(any(gates %in% catalog$quantities$canonical_name))
  expect_false(any(gates %in% catalog$aliases$alias))
  expect_false(any(gates %in% colnames(JAGS_materialize_draws(fit)[[1L]])))

  sd_total <- parameter_catalog_resolve(catalog, "split: sd_total", "mu")
  var_total <- parameter_catalog_resolve(catalog, "split: var_total", "mu")
  var_prop <- parameter_catalog_resolve(catalog, "split: var_prop(study)", "mu")
  expect_identical(sd_total$quantities$source_type, "composite")
  expect_identical(var_prop$quantities$source_type, "composite")
  expect_identical(
    as.numeric(parameter_draws(fit, sd_total)[[1L]][, 1L]),
    c(0, 2)
  )
  expect_identical(
    as.numeric(parameter_draws(fit, var_total)[[1L]][, 1L]),
    c(0, 4)
  )
  expect_identical(
    as.numeric(parameter_draws(fit, var_prop)[[1L]][, 1L]),
    c(NA_real_, 0.4)
  )

  mixed <- sd_total
  intercept <- parameter_catalog_resolve(catalog, "mu_intercept")
  mixed$quantity_id <- c(intercept$quantity_id, sd_total$quantity_id)
  mixed$quantities <- rbind(intercept$quantities, sd_total$quantities)
  expect_error(
    parameter_draws(fit, mixed),
    "Mixed or multiple derived selections are not supported in one extraction call.",
    fixed = TRUE
  )
})

test_that("gated totals without a source coordinate are unavailable", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  build_fit <- function(sd_source, monitor_source){
    formula_result <- JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        allocation = random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd_source = sd_source,
          weights = prior("dirichlet", list(alpha = c(2, 2))),
          inclusion = list(study = prior("beta", list(2, 2)))
        )
      )
    )
    allocation <- formula_result$formula_design$random_allocations[[1L]]
    gates <- .bt_random_effect_summary_allocation_gate_names(allocation)
    expect_length(gates, 1L)
    samples <- cbind(
      mu_intercept = c(0, 0),
      tau = c(2, 2),
      matrix(c(0.4, 0.6, 0.4, 0.6), nrow = 2L, byrow = TRUE,
             dimnames = list(NULL, paste0(allocation$weight_name, "[", 1:2, "]"))),
      matrix(c(0, 1), ncol = 1L, dimnames = list(NULL, gates))
    )
    if(!monitor_source){
      samples <- samples[, colnames(samples) != "tau", drop = FALSE]
    }
    fit <- structure(
      list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(fit, "prior_list") <- formula_result$prior_list
    attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
    attach_test_parameter_map(fit)
  }
  totals <- c("(mu) total_re: sd_total", "(mu) total_re: var_total")
  proportions <- c(
    "(mu) total_re: var_prop(study)",
    "(mu) total_re: var_prop(drug)"
  )
  draws <- function(fit, label){
    as.numeric(parameter_draws(
      fit,
      parameter_catalog_resolve(parameter_catalog(fit), label, "mu")
    )[[1L]][, 1L])
  }

  unavailable <- list(
    row = build_fit(random_sd_source("tau", shape = "row"), FALSE),
    scalar_unmonitored = build_fit(random_sd_source("tau"), FALSE)
  )
  for(case in names(unavailable)){
    fit <- unavailable[[case]]
    quantities <- parameter_catalog(fit)$quantities
    expect_false(any(totals %in% quantities$canonical_name), info = case)
    expect_true(all(proportions %in% quantities$canonical_name), info = case)
    estimates <- JAGS_estimates_table(fit, return_samples = TRUE)
    expect_true(all(proportions %in% colnames(estimates)), info = case)
    expect_false(any(totals %in% colnames(estimates)), info = case)
    # var_prop is the realized share: the gated study component is off in
    # the first draw.
    expect_equal(draws(fit, proportions[1L]), c(0, 0.4), info = case)
    expect_equal(draws(fit, proportions[2L]), c(1, 0.6), info = case)
    expect_error(
      random_effects_summary_posterior(fit, "sd_total"),
      "No random-effect sd_total summaries are available",
      fixed = TRUE,
      info = case
    )
  }

  # A monitored scalar source keeps the realized total with its gate and
  # weight dependencies: tau * sqrt(I_study * w_study + w_drug).
  fit <- build_fit(random_sd_source("tau"), TRUE)
  sd_total <- parameter_catalog_resolve(parameter_catalog(fit), totals[1L], "mu")
  expect_true("tau" %in% sd_total$quantities$extraction_key[[1L]]$dependencies)
  expect_equal(draws(fit, totals[1L]), 2 * sqrt(c(0.6, 1)), tolerance = 1e-12)
  expect_equal(draws(fit, totals[2L]), 4 * c(0.6, 1), tolerance = 1e-12)
})

test_that("block SDs scaled by an unmonitored allocation source are unavailable", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  build_fit <- function(monitor_source, inclusion){
    formula_result <- JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        allocation = random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd_source = random_sd_source("tau"),
          weights = prior("dirichlet", list(alpha = c(2, 2))),
          inclusion = inclusion
        )
      )
    )
    allocation <- formula_result$formula_design$random_allocations[[1L]]
    gates <- .bt_random_effect_summary_allocation_gate_names(allocation)
    samples <- cbind(
      mu_intercept = c(0, 0),
      tau = c(2, 2),
      matrix(c(0.4, 0.6, 0.4, 0.6), nrow = 2L, byrow = TRUE,
             dimnames = list(NULL, paste0(allocation$weight_name, "[", 1:2, "]")))
    )
    if(length(gates) > 0L){
      samples <- cbind(
        samples,
        matrix(c(0, 1), ncol = 1L, dimnames = list(NULL, gates))
      )
    }
    if(!monitor_source){
      samples <- samples[, colnames(samples) != "tau", drop = FALSE]
    }
    fit <- structure(
      list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(fit, "prior_list") <- formula_result$prior_list
    attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
    attach_test_parameter_map(fit)
  }
  block_sds <- c(
    "(mu) study: sd(intercept)", "(mu) study: var(intercept)",
    "(mu) drug: sd(intercept)", "(mu) drug: var(intercept)"
  )
  gated <- list(study = prior("beta", list(2, 2)))

  for(inclusion in list(NULL, gated)){
    info <- if(is.null(inclusion)) "ungated" else "gated"
    fit <- build_fit(monitor_source = FALSE, inclusion = inclusion)
    catalog <- parameter_catalog(fit)
    quantities <- catalog$quantities[
      catalog$quantities$canonical_name %in% block_sds, ,
      drop = FALSE
    ]
    expect_setequal(quantities$canonical_name, block_sds)
    expect_true(all(quantities$status == "unavailable"), info = info)
    expect_true(all(is.na(quantities$fixed_value)), info = info)
    estimates <- JAGS_estimates_table(
      fit,
      random_effects_summary = "full",
      return_samples = TRUE
    )
    expect_false(any(block_sds %in% colnames(estimates)), info = info)
    expect_true(
      "(mu) total_re: var_prop(study)" %in% colnames(estimates),
      info = info
    )
    expect_error(
      parameter_draws(
        fit,
        parameter_catalog_resolve(catalog, block_sds[1L], "mu")
      ),
      paste0(
        "The parameter quantity '(mu) study: sd(intercept)' is unavailable ",
        "in this fit: its source coordinates are not part of the posterior ",
        "draws. Refit the model with those coordinates monitored."
      ),
      fixed = TRUE,
      class = "BayesTools_refit_monitoring",
      info = info
    )
  }

  # A monitored source keeps the realized component SDs tau * sqrt(I * w).
  fit <- build_fit(monitor_source = TRUE, inclusion = gated)
  catalog <- parameter_catalog(fit)
  draws <- function(label){
    selection <- parameter_catalog_resolve(catalog, label, "mu")
    expect_identical(selection$quantities$status, "sampled")
    as.numeric(parameter_draws(fit, selection)[[1L]][, 1L])
  }
  expect_equal(draws(block_sds[1L]), 2 * sqrt(0.4) * c(0, 1), tolerance = 1e-12)
  expect_equal(draws(block_sds[3L]), 2 * sqrt(c(0.6, 0.6)), tolerance = 1e-12)
})

test_that("known group-covariance scale remains sd/var rather than sd_mult", {

  data <- data.frame(
    id = factor(c("a", "b", "a", "c"), levels = c("a", "b", "c"))
  )
  kernel <- matrix(
    c(2, .4, .2,
      .4, 3, .5,
      .2, .5, 4),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  formula_result <- JAGS_formula(
    formula = random_effects_formula(
      ~ 1 | id,
      group_covariance = random_group_covariance(kernel, scale = "none")
    ),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 1)))
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
  random <- catalog$quantities[
    startsWith(catalog$quantities$role, "random_"),
    ,
    drop = FALSE
  ]

  expect_setequal(random$quantity, c("sd", "var"))
  expect_setequal(
    random$canonical_name,
    c("(mu) sd(intercept)", "(mu) var(intercept)")
  )
  expect_false(any(random$quantity %in% c("sd_mult", "var_mult")))
  expect_error(
    parameter_catalog_resolve(catalog, "sd_mult", "mu"),
    "No public parameter quantity matches"
  )
})


test_that("extended catalogs may reuse a canonical name across components", {

  # `canonical_name` is a selector, not a key. A provider extending the catalog
  # routinely describes the same underlying term as BayesTools does - a factor
  # moderator gets both the backend coefficient and the extending provider's
  # grouped view - so the two legitimately share a public name and differ by
  # component. Rejecting that at construction made the whole catalog
  # unbuildable, with an error telling the user to refit, which cannot help.
  coordinates <- .bt_build_parameter_coordinates(
    columns = "mu_group",
    prior_list = list(mu_group = prior("normal", list(0, 1)))
  )
  catalog <- .bt_build_parameter_catalog(coordinates)
  canonical <- catalog$quantities$canonical_name[[1L]]
  namespace <- catalog$quantities$namespace[[1L]]

  shared_name <- .bt_parameter_catalog_quantity(
    canonical_name = canonical,
    namespace = namespace,
    role = "formula_coefficient_group",
    extraction_key = list(type = "robma", dependencies = character())
  )
  shared_name$component <- "mods"
  shared_name$provider <- "RoBMA"
  shared_name$quantity_id <- "RoBMA::grouped"

  extended <- parameter_catalog_extend(
    catalog,
    quantities = shared_name,
    aliases = .bt_parameter_catalog_empty_aliases(),
    provider = "RoBMA"
  )
  expect_true(canonical %in% extended$quantities$canonical_name)
  expect_identical(sum(extended$quantities$canonical_name == canonical), 2L)

  # the collision surfaces where it is actionable - at resolution, as a typed
  # ambiguity the caller can narrow - rather than at construction
  expect_error(
    parameter_catalog_resolve(extended, canonical, namespace),
    class = "BayesTools_parameter_ambiguous"
  )
  expect_identical(
    parameter_catalog_resolve(
      extended, canonical, namespace, component = "mods"
    )$quantity_id,
    "RoBMA::grouped"
  )

  # genuinely indistinguishable rows stay refused
  indistinguishable <- shared_name
  indistinguishable$component <- catalog$quantities$component[[1L]]
  indistinguishable$quantity_id <- "RoBMA::indistinguishable"
  expect_error(
    parameter_catalog_extend(
      catalog,
      quantities = indistinguishable,
      aliases = .bt_parameter_catalog_empty_aliases(),
      provider = "RoBMA"
    ),
    "cannot be resolved"
  )
})

# Reference: the validator without its memo, which is the pre-memo behaviour.
test_that("catalog tables are validated once per content and modifications are rechecked", {

  memo_reset <- function(){
    .BayesTools_private$content_memo <- NULL
  }
  memo_reset()
  withr::defer(memo_reset())
  original <- .bt_validate_parameter_catalog_tables_uncached
  calls <- 0L
  testthat::local_mocked_bindings(
    .bt_validate_parameter_catalog_tables_uncached = function(...){
      calls <<- calls + 1L
      original(...)
    },
    .package = "BayesTools"
  )

  prior_list <- list(
    memo_a = prior("normal", list(0, 1)),
    memo_b = prior("gamma", list(2, 2))
  )
  coordinates <- .bt_build_parameter_coordinates(
    columns = c("memo_a", "memo_b"),
    prior_list = prior_list
  )
  catalog <- .bt_build_parameter_catalog(coordinates, prior_list = prior_list)

  # the first validation runs the checks, later ones on the same object and on
  # an equal copy of it (as after saving and loading) do not
  expect_true(.bt_validate_parameter_catalog(catalog))
  expect_identical(calls, 1L)
  for(i in 1:3){
    expect_true(.bt_validate_parameter_catalog(catalog))
  }
  rebuilt <- unserialize(serialize(catalog, NULL))
  expect_true(.bt_validate_parameter_catalog(rebuilt))
  wrapper <- .bt_parameter_map_catalog(.bt_parameter_map_new(
    coordinates, catalog$quantities, catalog$aliases
  ))
  expect_true(.bt_validate_parameter_catalog(wrapper))
  expect_identical(calls, 1L)

  # a modification is never covered by the validated original: an invalid one
  # fails on every call (a failure is not recorded), a valid one is checked once
  spoofed <- catalog
  spoofed$quantities$extraction_key[[1L]] <- list(type = "bogus", dependencies = "memo_a")
  for(i in 1:2){
    expect_error(.bt_validate_parameter_catalog(spoofed), "extraction keys are malformed")
  }
  expect_identical(calls, 3L)
  relabelled <- catalog
  relabelled$quantities$display_label[[1L]] <- "a different label"
  expect_true(.bt_validate_parameter_catalog(relabelled))
  expect_true(.bt_validate_parameter_catalog(relabelled))
  expect_identical(calls, 4L)
  expect_true(.bt_validate_parameter_catalog(catalog))
  expect_identical(calls, 4L)

  # tables edited in place are as unrecognised as a rebuilt catalog
  edited <- catalog
  edited$aliases$quantity_id[[1L]] <- "BayesTools::unknown"
  expect_error(.bt_validate_parameter_catalog(edited), "unknown quantity IDs")
  expect_error(.bt_validate_parameter_catalog(edited), "unknown quantity IDs")

  # the memo and the unmemoized validator agree on the message of every refusal
  duplicated <- catalog
  duplicated$quantities <- rbind(duplicated$quantities, duplicated$quantities)
  duplicated$quantities$quantity_id[3:4] <- paste0("BayesTools::duplicate", 1:2)
  missing_metadata <- catalog
  missing_metadata$quantities$definedness[[1L]] <- NA_character_
  refusals <- list(
    spoofed = spoofed,
    edited = edited,
    duplicated = duplicated,
    missing_metadata = missing_metadata
  )
  for(name in names(refusals)){
    refused <- refusals[[name]]
    memoized <- tryCatch(.bt_validate_parameter_catalog(refused), error = identity)
    unmemoized <- tryCatch(
      original(refused$quantities, refused$aliases),
      error = identity
    )
    expect_s3_class(memoized, "error")
    expect_identical(conditionMessage(memoized), conditionMessage(unmemoized))
    expect_identical(class(memoized), class(unmemoized))
  }
  expect_true(original(catalog$quantities, catalog$aliases))
})

test_that("the validation memo is bounded and keeps selections separate from catalogs", {

  memo_reset <- function(){
    .BayesTools_private$content_memo <- NULL
  }
  memo_reset()
  withr::defer(memo_reset())
  original <- .bt_validate_parameter_catalog_tables_uncached
  calls <- 0L
  testthat::local_mocked_bindings(
    .bt_validate_parameter_catalog_tables_uncached = function(...){
      calls <<- calls + 1L
      original(...)
    },
    .package = "BayesTools"
  )

  catalogs <- lapply(seq_len(20L), function(i){
    name <- paste0("bounded_", i)
    prior_list <- stats::setNames(list(prior("normal", list(0, 1))), name)
    .bt_build_parameter_catalog(
      .bt_build_parameter_coordinates(columns = name, prior_list = prior_list),
      prior_list = prior_list
    )
  })
  # construction validates each catalog; start counting from an empty memo
  memo_reset()
  calls <- 0L
  for(catalog in catalogs){
    .bt_validate_parameter_catalog(catalog)
  }
  expect_identical(calls, 20L)
  limit <- .bt_content_memo_limit()
  expect_identical(length(.BayesTools_private$content_memo$parameter_catalog_tables), limit)
  # the newest entries are kept, the oldest ones are checked again
  .bt_validate_parameter_catalog(catalogs[[20L]])
  expect_identical(calls, 20L)
  .bt_validate_parameter_catalog(catalogs[[1L]])
  expect_identical(calls, 21L)

  # a selection is validated once per content, in its own memo
  selection <- parameter_catalog_resolve(catalogs[[5L]], "bounded_5")
  calls_before <- calls
  for(i in 1:3){
    .bt_validate_parameter_selection(selection, catalog = catalogs[[5L]])
  }
  expect_identical(calls - calls_before, 0L)
  expect_identical(
    length(.BayesTools_private$content_memo$parameter_selection_tables),
    1L
  )
  stale <- selection
  stale$quantities$display_label[[1L]] <- "tampered"
  expect_error(
    .bt_validate_parameter_selection(stale, catalog = catalogs[[5L]]),
    "stale or does not belong"
  )
})

# The catalog builders assemble their tables without the data frame operations
# they used to apply row by row (one-row tables, `[<-`, `rbind()`, row
# subsetting). Every replaced step is compared with the step it replaces.
catalog_builder_test_inputs <- function(){

  data <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    f = factor(c("a", "b", "c", "a", "b", "c")),
    id = factor(c("a", "a", "b", "b", "c", "c"))
  )
  formula_result <- JAGS_formula(
    ~ 1 + x + f + us(1 + x | id),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)), cor = prior_lkj())
    ),
    formula_scale = list(x = TRUE)
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  cholesky <- random_term$correlation$cholesky_name
  columns <- c(
    "mu_intercept", "mu_x",
    .JAGS_prior_factor_names("mu_f", formula_result$prior_list$mu_f),
    random_term$sd_parameter_names,
    paste0(cholesky, c("[1,1]", "[2,1]", "[2,2]"))
  )
  prior_list <- formula_result$prior_list
  formula_design <- list(mu = formula_result$formula_design)
  formula_scale <- list(mu = formula_result$formula_scale)
  coordinates <- .bt_build_parameter_coordinates(
    columns = columns,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  list(
    coordinates = coordinates,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale,
    catalog = .bt_build_parameter_catalog(
      coordinates, prior_list, formula_design, formula_scale
    )
  )
}

test_that("catalog quantity rows and their table equal the data frame operations they replace", {

  inputs <- catalog_builder_test_inputs()
  quantities <- inputs$catalog$quantities
  expect_gt(nrow(quantities), 8L)
  expect_true(all(c("factor_level", "coordinate", "random_summary") %in%
                    vapply(quantities$extraction_key, `[[`, character(1), "type")))

  # a row built directly is the row the assignment of its values builds
  reference_row <- function(row){
    out <- .bt_parameter_catalog_empty_quantities_build()
    scalar_columns <- setdiff(
      names(out), c("arguments", "support", "label_parts", "extraction_key")
    )
    out[1L, scalar_columns] <- unname(lapply(scalar_columns, function(name) row[[name]]))
    out$arguments <- I(list(row$arguments[[1L]]))
    out$support <- I(list(NULL))
    out$label_parts <- I(list(row$label_parts[[1L]]))
    out$extraction_key <- I(list(row$extraction_key[[1L]]))
    out
  }
  parts <- .bt_label_parts("theta", selector = "theta")
  rows <- list(
    .bt_parameter_catalog_quantity(
      "theta", "model", "parameter", label_parts = parts, status = "sampled",
      source_type = "identity",
      extraction_key = list(type = "coordinate", dependencies = "theta")
    ),
    .bt_parameter_catalog_quantity(
      "fixed", "model", "parameter", status = "structural", fixed_value = 0L,
      internal = TRUE, arguments = c("a", "b"), source_type = "structural_zero",
      extraction_key = list(type = "factor_level", dependencies = character(),
                            weights = numeric())
    ),
    .bt_parameter_catalog_quantity(
      "extension", "mu", "parameter", formula_parameter = "mu", term = "x",
      component = "{1}", fixed_value = NA, source_type = "composite",
      extraction_key = list(type = "coordinate", dependencies = "x")
    )
  )
  for(row in rows){
    expect_identical(row, reference_row(row))
    expect_identical(names(row), .bt_parameter_catalog_quantity_columns)
  }
  expect_identical(rows[[2L]]$fixed_value, 0)

  # the table of rows is their rbind()
  reference_bind <- function(rows){
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  expect_identical(.bt_parameter_catalog_bind_quantities(rows), reference_bind(rows))
  expect_identical(.bt_parameter_catalog_bind_quantities(rows[1L]), reference_bind(rows[1L]))
  subset_rows <- lapply(seq_len(nrow(quantities)), function(i) quantities[i, , drop = FALSE])
  expect_identical(.bt_parameter_catalog_bind_quantities(subset_rows), reference_bind(subset_rows))
  expect_identical(.bt_parameter_catalog_bind_quantities(subset_rows), quantities)

  # rows read as the fields of the one-row table: the (row, field) pairs whose
  # value differs from the one-row table's are named
  differing_fields <- function(table){
    lists <- .bt_parameter_catalog_rows(table)
    expect_length(lists, nrow(table))
    as.character(unlist(lapply(seq_along(lists), function(i){
      reference <- table[i, , drop = FALSE]
      if(!identical(names(lists[[i]]), names(reference))){
        return(paste(i, "names"))
      }
      differs <- !vapply(names(reference), function(field){
        identical(lists[[i]][[field]], reference[[field]])
      }, logical(1))
      if(any(differs)) paste(i, names(reference)[differs])
    })))
  }
  expect_identical(differing_fields(quantities), character())
  expect_identical(differing_fields(inputs$coordinates), character())
})

test_that("catalog quantity rows refuse fields that do not hold one value", {

  # The assignment into the empty table refused a field without a value; the
  # directly built row must refuse it too instead of carrying a column of
  # another length (which the table validator does not see).
  build <- function(...){
    .bt_parameter_catalog_quantity(
      "theta", "model", "parameter", status = "sampled",
      source_type = "identity",
      extraction_key = list(type = "coordinate", dependencies = "theta"),
      ...
    )
  }
  expect_identical(nrow(build()), 1L)
  for(field in list(list(component = NULL), list(component = character()),
                    list(term = character()), list(fixed_value = numeric()),
                    list(owner_name = c("a", "b")), list(internal = c(TRUE, FALSE)))){
    expect_error(
      do.call(build, field),
      paste0("Parameter catalog quantity fields must each hold one value: '",
             names(field), "'."),
      fixed = TRUE
    )
  }
})

test_that("catalog support and definedness share the prior-list work without changing it", {

  inputs <- catalog_builder_test_inputs()
  quantities <- inputs$catalog$quantities
  object <- structure(
    list(),
    prior_list = inputs$prior_list,
    formula_design = inputs$formula_design,
    formula_scale = inputs$formula_scale
  )
  # reference: one data frame row and one full prior-list pass per quantity
  reference <- quantities
  reference$support <- I(lapply(seq_len(nrow(reference)), function(i){
    .bt_parameter_catalog_quantity_support(object, quantities[i, , drop = FALSE])
  }))
  reference$definedness <- vapply(seq_len(nrow(reference)), function(i){
    .bt_parameter_catalog_quantity_definedness(object, quantities[i, , drop = FALSE])
  }, character(1))
  expect_identical(
    .bt_parameter_catalog_add_support(
      quantities, inputs$prior_list, inputs$formula_design, inputs$formula_scale
    ),
    reference
  )
  expect_identical(quantities$support, reference$support)
  expect_true(any(!vapply(reference$support, is.null, logical(1))))
  expect_true(any(quantities$definedness == "correlation"))

  # the shared context is what a single call computes
  context <- .bt_parameter_catalog_linear_support_context(inputs$prior_list)
  weights <- c(mu_intercept = 1, mu_x = -2)
  expect_identical(
    .bt_parameter_catalog_linear_support(inputs$prior_list, weights, context = context),
    .bt_parameter_catalog_linear_support(inputs$prior_list, weights)
  )
  expect_null(.bt_parameter_catalog_linear_support(inputs$prior_list, c(unowned = 1), context = context))
  expect_null(.bt_parameter_catalog_linear_support_context(list()))
})

test_that("catalog aliases equal the per-quantity data frame assembly they replace", {

  inputs <- catalog_builder_test_inputs()
  quantities <- inputs$catalog$quantities

  reference_alias_parts <- function(alias, candidates = list()){
    for(parts in candidates){
      if(identical(.bt_label(parts, style = "table"), alias)){
        return(parts)
      }
    }
    .bt_label_parts(alias, selector = alias)
  }
  reference_label_aliases <- function(quantity, formula_scale = NULL){
    parts <- quantity$label_parts[[1L]]
    if(is.null(parts)){
      return(list(values = character(), parts = list()))
    }
    original_scale <- .bt_label_parts_log_intercept(parts, formula_scale)[[1L]]
    renderings <- list(
      .bt_parameter_catalog_alias_rendering(parts, formula_prefix = TRUE),
      .bt_parameter_catalog_alias_rendering(parts, formula_prefix = FALSE),
      .bt_parameter_catalog_alias_rendering(original_scale, formula_prefix = TRUE),
      .bt_parameter_catalog_alias_rendering(original_scale, formula_prefix = FALSE)
    )
    if(.bt_parameter_catalog_is_factor_quantity(quantity) &&
       length(parts$levels) > 0L){
      dif <- .bt_label_parts_update(parts, transformation = "dif")[[1L]]
      selector <- .bt_label(dif, style = "selector")
      renderings <- c(
        renderings,
        list(
          reference_alias_parts(selector),
          .bt_parameter_catalog_alias_rendering(dif, formula_prefix = TRUE),
          .bt_parameter_catalog_alias_rendering(dif, formula_prefix = FALSE)
        )
      )
    }
    .bt_parameter_catalog_rendered_aliases(renderings)
  }
  reference_aliases <- function(quantities, formula_design, formula_scale){
    public <- quantities[!quantities$internal, , drop = FALSE]
    rows <- list()
    secondary <- list()
    add_aliases <- function(quantity, aliases, simplified){
      keep <- !is.na(aliases$values) & nzchar(aliases$values) &
        !duplicated(aliases$values)
      values <- aliases$values[keep]
      if(length(values) == 0L){
        return(invisible(NULL))
      }
      row <- data.frame(
        alias = values,
        quantity_id = rep(quantity$quantity_id, length(values)),
        namespace = rep(quantity$namespace, length(values)),
        component = rep(quantity$component, length(values)),
        simplified = rep(simplified, length(values)),
        stringsAsFactors = FALSE
      )
      row$label_parts <- I(unname(aliases$parts[keep]))
      rows[[length(rows) + 1L]] <<- row
      invisible(NULL)
    }
    for(i in seq_len(nrow(public))){
      quantity <- public[i, , drop = FALSE]
      parts <- quantity$label_parts[[1L]]
      if(startsWith(quantity$role, "random_")){
        aliases <- .bt_parameter_catalog_rendered_aliases(c(
          if(!is.null(parts)){
            list(.bt_parameter_catalog_alias_rendering(parts, formula_prefix = FALSE))
          },
          .bt_parameter_catalog_random_correlation_aliases(quantity, formula_design)
        ))
      }else{
        structured <- if(!is.null(parts)){
          list(parts, .bt_parameter_catalog_alias_rendering(parts, formula_prefix = FALSE))
        }else{
          list()
        }
        values <- unique(c(
          quantity$canonical_name,
          quantity$display_label,
          quantity$term,
          if(nzchar(quantity$term) && nzchar(quantity$component) &&
             (identical(quantity$role, "fixed_coefficient") ||
                .bt_parameter_catalog_is_factor_quantity(quantity))){
            .bt_parameter_catalog_level_alias(quantity$term, quantity$component)
          }else{
            character()
          }
        ))
        values <- values[!is.na(values) & nzchar(values)]
        aliases <- list(
          values = values,
          parts  = lapply(values, reference_alias_parts, candidates = structured)
        )
        label_aliases <- reference_label_aliases(quantity, formula_scale)
        new_labels <- !label_aliases$values %in% values
        if(any(new_labels)){
          secondary_rows <- data.frame(
            alias = label_aliases$values[new_labels],
            quantity_id = quantity$quantity_id,
            stringsAsFactors = FALSE
          )
          secondary_rows$label_parts <- I(unname(label_aliases$parts[new_labels]))
          secondary[[length(secondary) + 1L]] <- secondary_rows
        }
      }
      add_aliases(quantity, aliases, simplified = FALSE)
      if(startsWith(quantity$role, "random_")){
        add_aliases(
          quantity,
          .bt_parameter_catalog_random_simplified_aliases(quantity),
          simplified = TRUE
        )
      }
    }
    out <- do.call(rbind, rows)
    out <- out[!duplicated(out[setdiff(names(out), "label_parts")]), , drop = FALSE]
    rownames(out) <- NULL
    .bt_parameter_catalog_secondary_aliases(
      out,
      public,
      secondary = if(length(secondary) > 0L) do.call(rbind, secondary)
    )
  }

  # the fitted scaling, and the same catalog with a log intercept, whose
  # original-scale label differs from the fitted one
  log_scale <- inputs$formula_scale
  attr(log_scale$mu, "log_intercept") <- TRUE
  for(formula_scale in list(inputs$formula_scale, log_scale)){
    aliases <- .bt_parameter_catalog_aliases(quantities, inputs$formula_design, formula_scale)
    expect_gt(nrow(aliases), nrow(quantities))
    expect_identical(
      aliases,
      reference_aliases(quantities, inputs$formula_design, formula_scale)
    )
  }
  expect_identical(inputs$catalog$aliases, reference_aliases(
    quantities, inputs$formula_design, inputs$formula_scale
  ))
  expect_false(identical(
    .bt_parameter_catalog_aliases(quantities, inputs$formula_design, log_scale),
    inputs$catalog$aliases
  ))

  # the precomputed candidate labels select the candidate the loop selects
  parts <- quantities$label_parts[[which(!vapply(quantities$label_parts, is.null, logical(1)))[[1L]]]]
  candidates <- list(parts, .bt_parameter_catalog_alias_rendering(parts, formula_prefix = FALSE))
  labels <- vapply(candidates, .bt_label, character(1), style = "table")
  for(alias in c(labels, "no such label")){
    expect_identical(
      .bt_parameter_catalog_alias_parts(alias, candidates, labels),
      reference_alias_parts(alias, candidates)
    )
    expect_identical(
      .bt_parameter_catalog_alias_parts(alias, candidates),
      reference_alias_parts(alias, candidates)
    )
  }
  expect_identical(
    .bt_parameter_catalog_alias_parts("no candidates"),
    reference_alias_parts("no candidates")
  )
})
