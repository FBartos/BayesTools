skip_if_not_test_profile("unit")

.random_update_test_fit <- function(formula, data){

  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("gamma", list(2, 2)))
  )
  term <- result$formula_design$random_effects[[1L]]
  correlation <- term$correlation
  correlation_columns <- if(is.null(correlation)){
    character()
  }else if(identical(correlation$type, "lkj")){
    correlation$primitive_names
  }else{
    correlation$rho_name
  }
  columns <- unique(c(
    "mu_intercept",
    term$sd_parameter_names,
    correlation_columns
  ))
  draws <- matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}


# The factor form of an affine allocation-derived plan: the update of a block
# that the allocation scale drives only in part, which the grid evaluates. The
# plans of the fits below are affine because the scale drives every column.
.random_update_test_factor_form <- function(fit, plan){

  term <- Filter(
    function(term) identical(term$block_name, plan$blocks[[1L]]),
    attr(fit, "formula_design")$mu$random_effects
  )[[1L]]
  variance <- identical(plan$coefficient_transform, list(type = "identity"))
  plan$family <- "factor"
  plan$update <- "column_scale"
  plan$coefficient_scale_transform <- if(variance) list(type = "sqrt") else list(type = "identity")
  plan$coefficient_transform <- NULL
  plan$driven_columns <- seq_len(term$n_columns)

  plan
}


.random_update_test_plan <- function(fit, role, component = NULL){

  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities[
    catalog$quantities$role == role,
    ,
    drop = FALSE
  ]
  if(!is.null(component)){
    quantities <- quantities[quantities$component == component, , drop = FALSE]
  }
  stopifnot(nrow(quantities) == 1L)
  selection <- parameter_catalog_resolve(
    catalog,
    quantities$canonical_name[[1L]],
    "mu"
  )
  random_effects_marginal_update_plan(fit, selection)
}


test_that("marginal update plans propagate invalid transformation metadata", {

  fit <- .random_update_test_fit(
    ~ 1 + (1 | study),
    data.frame(study = factor(c("a", "a", "b", "b")))
  )
  testthat::local_mocked_bindings(
    .bt_parameter_transform_from_quantity = function(...){
      stop("Invalid source transformation metadata.", call. = FALSE)
    },
    .package = "BayesTools"
  )
  expect_error(.random_update_test_plan(fit, "random_sd"),
               "Invalid source transformation metadata.", fixed = TRUE)
})


.random_update_test_allocation_fit <- function(
    formula = ~ 1 + random(1 + x | id, name = "study", covariance = "diag")){

  data <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    x = c(-1, 0, 1, 2),
    index = c(1, 2, 1, 2)
  )
  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "total_variance",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  weight <- "mu__xRE_ALLOCx_allocation__weight"
  term <- result$formula_design$random_effects[[1L]]
  columns <- c(
    term$correlation$rho_name,
    "mu_intercept", "tau",
    result$formula_design$random_effects[[1L]]$sd_parameter_names,
    paste0(weight, "[", 1:2, "]"),
    paste0("prior_par_eta_", weight, "[", 1:2, "]")
  )
  draws <- matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}


.random_update_test_generic_fit <- function(formula, data, prior_random,
                                            prior_list, formula_scale = NULL,
                                            extra_columns = character(),
                                            values = NULL){

  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = prior_list,
    formula_scale = formula_scale,
    prior_random = prior_random
  )
  columns <- character()
  for(name in names(result$prior_list)){
    current <- result$prior_list[[name]]
    if(is.prior.point(current)){
      next
    }
    if(is.prior.simplex(current)){
      K <- current$parameters$K
      columns <- c(
        columns,
        paste0(name, "[", seq_len(K), "]"),
        paste0("prior_par_eta_", name, "[", seq_len(K), "]")
      )
    }else{
      columns <- c(columns, name)
    }
  }
  for(term in result$formula_design$random_effects){
    correlation <- term$correlation
    if(!is.null(correlation)){
      columns <- c(columns, if(identical(correlation$type, "lkj")){
        correlation$primitive_names
      }else{
        correlation$sample_name
      })
    }
  }
  columns <- unique(c(columns, extra_columns))
  draws <- matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  if(length(values) > 0L){
    stopifnot(all(names(values) %in% columns))
    draws[, names(values)] <- rep(values, each = nrow(draws))
  }
  scale <- if(is.null(result$formula_scale)){
    NULL
  }else{
    list(mu = result$formula_scale)
  }
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "formula_scale") <- scale
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design),
    formula_scale = scale
  )
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}


test_that("composite random SD sources have no single-coordinate update", {

  data <- data.frame(
    id = factor(rep(letters[1:4], each = 3L)),
    x = c(1.2, 3.5, 4.1, 2.2, 5.3, 0.4, 3.3, 2.9, 6.1, 1.7, 4.4, 2.6)
  )
  fixed <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  sd_prior <- prior_random(sd = prior("gamma", list(2, 2)))
  formulas <- list(
    diag = ~ 1 + x + (1 + x || id),
    us = ~ 1 + x + (1 + x | id)
  )

  # With a scaled slope, the public intercept SD mixes the fitted SDs (and
  # the correlation), so no single source coordinate carries the update. The
  # public slope SD is the fitted slope SD over sd(x): a one-to-one rescaling
  # of its own coordinate, so its update is in that coordinate and the plan
  # carries the transform from the coordinate to the public quantity.
  scale <- 1 / stats::sd(data$x)
  for(structure in names(formulas)){
    fit <- .random_update_test_generic_fit(
      formulas[[structure]], data, sd_prior, fixed,
      formula_scale = list(x = TRUE)
    )
    catalog <- parameter_catalog(fit)
    quantities <- catalog$quantities[
      catalog$quantities$role %in% c("random_sd", "random_var"), ,
      drop = FALSE
    ]
    expect_identical(nrow(quantities), 4L, info = structure)
    expect_identical(
      quantities$source_type[quantities$component == "intercept"],
      c("composite", "composite"),
      info = structure
    )
    expect_identical(
      quantities$source_type[quantities$component == "x"],
      c("one_to_one_transform", "one_to_one_transform"),
      info = structure
    )
    for(name in quantities$canonical_name[quantities$component == "intercept"]){
      plan <- random_effects_marginal_update_plan(
        fit,
        parameter_catalog_resolve(catalog, name, "mu")
      )
      expect_identical(plan$family, "unsupported", info = name)
      expect_identical(plan$reason, "composite_sd_source", info = name)
    }
    for(role in c("random_sd", "random_var")){
      plan <- .random_update_test_plan(fit, role, "x")
      expect_identical(plan$family, if(structure == "diag") "affine" else "factor",
                       info = paste(structure, role))
      expect_identical(plan$update, "column_scale", info = paste(structure, role))
      expect_identical(plan$source_parameter, "mu__xREx__id_x",
                       info = paste(structure, role))
      expect_identical(
        plan$source_transform_spec,
        if(role == "random_sd") {
          list(type = "affine", offset = 0, scale = scale)
        }else{
          list(type = "square", scale = scale)
        },
        info = paste(structure, role)
      )
    }
  }
  # the factor update of the correlated block replaces the slope SD
  # coordinate: the candidate scales are the fitted values
  us_fit <- .random_update_test_generic_fit(
    formulas$us, data, sd_prior, fixed, formula_scale = list(x = TRUE)
  )
  us_plan <- .random_update_test_plan(us_fit, "random_sd", "x")
  us_grid <- random_effects_marginal_update_grid(
    us_fit, us_plan, values = c(0.3, 0.6)
  )
  expect_equal(us_grid$candidate_scale[, 1L], c(0.3, 0.6))

  unscaled <- .random_update_test_generic_fit(
    formulas$us, data, sd_prior, fixed
  )
  plan <- .random_update_test_plan(unscaled, "random_sd", "x")
  expect_identical(plan$family, "factor")
  expect_identical(plan$source_parameter, "mu__xREx__id_x")
  grid <- random_effects_marginal_update_grid(
    unscaled, plan, values = c(0.3, 0.6)
  )
  expect_equal(grid$candidate_scale[, 1L], c(0.3, 0.6))
})


test_that("scaled allocation-derived SDs have no quantity-scale factor grid", {

  data <- data.frame(
    id = factor(rep(c("a", "b", "c"), each = 4L)),
    x = c(1.2, 3.5, 4.1, 2.2, 5.3, 0.4, 3.3, 2.9, 6.1, 1.7, 4.4, 2.6)
  )
  fixed <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  allocation <- prior_random(
    allocation = random_variance_allocation(
      name = "allocation",
      terms = "study",
      target = "sd_component",
      scale = "total_variance",
      sd_source = random_sd_source("tau"),
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    )
  )
  formula <- ~ 1 + x + random(1 + x | id, name = "study")
  unscaled <- .random_update_test_generic_fit(
    formula, data, allocation, fixed, extra_columns = "tau"
  )
  scaled <- .random_update_test_generic_fit(
    formula, data, allocation, fixed,
    formula_scale = list(x = TRUE), extra_columns = "tau"
  )

  for(component in c("intercept", "x")){
    for(role in c("random_sd", "random_var")){
      plan <- .random_update_test_plan(unscaled, role, component)
      expect_identical(plan$family, "affine")
      expect_identical(plan$coefficient_input, "quantity")
      plan <- .random_update_test_plan(scaled, role, component)
      expect_identical(plan$family, "unsupported")
      expect_identical(plan$reason, "scaled_component_sd")
    }
  }

  # Without formula scaling the public SD is the fitted coefficient scale
  # that a quantity-scale grid replaces.
  plan <- .random_update_test_factor_form(
    unscaled, .random_update_test_plan(unscaled, "random_sd", "x")
  )
  grid <- random_effects_marginal_update_grid(unscaled, plan, values = 1)
  public <- parameter_draws(
    unscaled,
    parameter_catalog_resolve(
      parameter_catalog(unscaled),
      plan$canonical_name,
      "mu"
    )
  )
  expect_equal(
    grid$coefficient_scale[, plan$component_index],
    unname(as.numeric(as.matrix(public)[, 1L])),
    tolerance = 1e-12
  )
})


test_that("scaled allocation-derived SDs have no quantity-scale affine plan", {

  data <- data.frame(
    id = factor(rep(c("a", "b", "c"), each = 4L)),
    x = c(1.2, 3.5, 4.1, 2.2, 5.3, 0.4, 3.3, 2.9, 6.1, 1.7, 4.4, 2.6),
    z = c(0.2, -0.5, 0.1, 0.4, -0.3, 0.6, -0.1, 0.3, -0.2, 0.5, 0.0, -0.4)
  )
  fixed <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1)),
    z = prior("normal", list(0, 1))
  )
  allocation <- prior_random(
    allocation = random_variance_allocation(
      name = "allocation",
      terms = "study",
      target = "sd_component",
      scale = "total_variance",
      sd_source = random_sd_source("tau"),
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    )
  )
  fit <- function(formula, formula_scale = NULL){
    .random_update_test_generic_fit(
      formula, data, allocation, fixed,
      formula_scale = formula_scale, extra_columns = "tau"
    )
  }
  scaled_slope <- ~ 1 + x + z + random(1 + x | id, name = "study", covariance = "diag")
  unscaled_slope <- ~ 1 + x + z + random(1 + z | id, name = "study", covariance = "diag")
  cases <- list(
    unscaled = list(fit = fit(scaled_slope), component = "x"),
    scale_free = list(
      fit = fit(unscaled_slope, formula_scale = list(x = TRUE)),
      component = "z"
    ),
    scaled = list(
      fit = fit(scaled_slope, formula_scale = list(x = TRUE)),
      component = "x"
    )
  )

  for(case in names(cases)){
    for(component in c("intercept", cases[[case]]$component)){
      for(role in c("random_sd", "random_var")){
        plan <- .random_update_test_plan(cases[[case]]$fit, role, component)
        info <- paste(case, role, component)
        if(identical(case, "scaled")){
          # The public SDs are on the original predictor scale: the intercept
          # SD mixes both fitted SDs and the slope SD is divided by sd(x).
          expect_identical(plan$family, "unsupported", info = info)
          expect_identical(plan$reason, "scaled_component_sd", info = info)
        }else{
          expect_identical(plan$family, "affine", info = info)
          expect_identical(plan$update, "scale", info = info)
          expect_identical(plan$coefficient_input, "quantity", info = info)
        }
      }
    }
  }
})


test_that("allocation plans report Dirichlet weight coordinates", {

  fit <- .random_update_test_allocation_fit()
  weight <- "mu__xRE_ALLOCx_allocation__weight"
  for(index in 1:2){
    component <- c("intercept", "x")[[index]]
    for(role in c("random_var_prop", "random_sd_mult")){
      plan <- .random_update_test_plan(fit, role, component)
      info <- paste(role, component)
      expect_identical(plan$family, "affine", info = info)
      expect_identical(plan$update, "allocation", info = info)
      expect_identical(plan$source_parameter, weight, info = info)
      expect_identical(plan$allocation$weight_name, weight, info = info)
      expect_equal(plan$allocation$index, index, info = info)
      expect_identical(plan$allocation$n_targets, 2L, info = info)
    }
  }
})


test_that("source-coordinate factor and Markov plans declare their input", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 2L)),
    time = rep(c(0, 1, 3), 2L)
  )
  hcs_sd <- .random_update_test_plan(
    .random_update_test_fit(~ 1 + hcs(f | id), data),
    "random_sd",
    "f[x]"
  )
  us_sd <- .random_update_test_plan(
    .random_update_test_fit(~ 1 + us(1 + f | id), data),
    "random_sd",
    "intercept"
  )
  ar1_cor <- .random_update_test_plan(
    .random_update_test_fit(~ 1 + ar1(f | id), data),
    "random_correlation"
  )
  car_cor <- .random_update_test_plan(
    .random_update_test_fit(~ 1 + car(time | id), data),
    "random_correlation"
  )

  for(plan in list(hcs_sd, us_sd)){
    expect_identical(plan$family, "factor")
    expect_identical(plan$coefficient_input, "source")
  }
  for(plan in list(ar1_cor, car_cor)){
    expect_identical(plan$family, "markov")
    expect_identical(plan$coefficient_input, "source")
  }
})


test_that("a source-coordinate plan lists every block that uses the source", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    other = factor(c("a", "b", "a", "b")),
    third = factor(c("u", "v", "v", "u")),
    x = c(0, 1, 0, 1)
  )
  tau <- random_sd_source("tau")
  allocation <- function(name, term){
    random_variance_allocation(
      name = name,
      terms = term,
      target = "sd_component",
      sd_source = tau,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    )
  }
  plan_for <- function(fit, quantity){
    catalog <- parameter_catalog(fit)
    random_effects_marginal_update_plan(
      fit,
      parameter_catalog_resolve(catalog, quantity)
    )
  }
  change <- function(fit, blocks, tau_value){
    posterior <- as.matrix(fit)
    reference <- random_effects_marginal_vcov(
      fit, parameter = "mu", blocks = blocks, posterior_samples = posterior
    )$samples[1L, , ]
    posterior[, "tau"] <- tau_value
    updated <- random_effects_marginal_vcov(
      fit, parameter = "mu", blocks = blocks, posterior_samples = posterior
    )$samples[1L, , ]
    unname(updated - reference)
  }

  # two diagonal random-slope blocks whose SD-component allocations share
  # random_sd_source("tau")
  fit <- .random_update_test_generic_fit(
    ~ 1 + random(1 + x | id, name = "first", covariance = "diag") +
      random(1 + x | other, name = "second", covariance = "diag"),
    df,
    prior_random(allocation = list(
      allocation("alloc1", "first"),
      allocation("alloc2", "second")
    )),
    list(intercept = prior("normal", list(0, 1))),
    extra_columns = "tau"
  )
  total <- plan_for(fit, "(mu) alloc1: sd_total")
  expect_identical(total$family, "affine")
  expect_identical(total$coefficient_input, "source")
  expect_identical(total$source_parameter, "tau")
  expect_identical(total$blocks, c("first", "second"))
  expect_identical(plan_for(fit, "(mu) alloc1: var_total")$blocks, c("first", "second"))
  expect_identical(plan_for(fit, "(mu) alloc2: sd_total")$blocks, c("second", "first"))

  # the plan's blocks reproduce the full change of the covariance
  expected <- 0.375 * (1 + outer(df$x, df$x)) * (
    outer(as.integer(df$id), as.integer(df$id), "==") +
      outer(as.integer(df$other), as.integer(df$other), "==")
  )
  full <- change(fit, blocks = NULL, tau_value = 1)
  expect_equal(full, unname(expected), tolerance = 1e-12)
  expect_equal(change(fit, total$blocks, tau_value = 1), full, tolerance = 1e-12)

  # a weight coordinate belongs to its own allocation, and a quantity-scale
  # component SD moves the shared scale, which reaches every block that uses it
  proportion <- plan_for(fit, "(mu) alloc1: var_prop(x)")
  expect_identical(proportion$source_parameter, "mu__xRE_ALLOCx_alloc1__weight")
  expect_identical(proportion$blocks, "first")
  component <- plan_for(fit, "(mu) first: sd(x)")
  expect_identical(component$coefficient_input, "quantity")
  expect_identical(component$blocks, c("first", "second"))

  # a third block bound directly to the same source is listed as well
  direct <- .random_update_test_generic_fit(
    ~ 1 + random(1 + x | id, name = "first", covariance = "diag") +
      random(1 + x | other, name = "second", covariance = "diag") +
      random(1 | third, name = "third", covariance = "diag"),
    df,
    prior_random(
      third = random_block(sd_source = tau),
      allocation = list(allocation("alloc1", "first"), allocation("alloc2", "second"))
    ),
    list(intercept = prior("normal", list(0, 1))),
    extra_columns = "tau"
  )
  direct_total <- plan_for(direct, "(mu) alloc1: sd_total")
  expect_identical(direct_total$blocks, c("first", "second", "third"))
  expect_equal(
    change(direct, direct_total$blocks, tau_value = 1),
    change(direct, NULL, tau_value = 1),
    tolerance = 1e-12
  )
  expect_gt(
    max(abs(change(direct, "first", tau_value = 1) - change(direct, NULL, tau_value = 1))),
    0.1
  )
  expect_identical(
    plan_for(direct, "(mu) first: sd(x)")$blocks,
    c("first", "second", "third")
  )
})

test_that("interaction SD labels keep backslashes and levels that differ by a character", {

  # The levels "a\b" and "ab" give two distinct cell labels with their
  # characters intact (they were once spliced into a sub() replacement).
  backslash <- intToUtf8(92L)
  lev_f <- c(paste0("a", backslash, "b"), "ab")
  data <- expand.grid(
    f = factor(lev_f, levels = lev_f), h = factor(c("u", "v")),
    id = factor(c("g1", "g2")), KEEP.OUT.ATTRS = FALSE
  )
  result <- JAGS_formula(
    ~ 1 + diag(0 + f:h | id), "mu", data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      sd = prior("gamma", list(2, 2)),
      id = random_block(contrasts = c(f = "independent", h = "independent"))
    )
  )
  term <- result$formula_design$random_effects[[1L]]
  columns <- unique(c("mu_intercept", term$sd_parameter_names))
  draws <- matrix(0.5, 2L, length(columns), dimnames = list(NULL, columns))
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns, prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_fit_contract(.bt_attach_draw_geometry(fit))
  names <- parameter_catalog(fit)$quantities$canonical_name
  sd_names <- names[startsWith(names, "(mu) sd(")]
  expect_length(sd_names, 4L)
  expect_false(anyDuplicated(sd_names) > 0L)
  expect_true(paste0("(mu) sd(f[a", backslash, "b]:h[u])") %in% sd_names)
  expect_true("(mu) sd(f[ab]:h[u])" %in% sd_names)
})

test_that("affine update plans are affine in their declared coefficient", {

  # A plan that claims `A + h(value) B` is checked against the covariance, not
  # only against its labels: the divided differences of C(value) in h(value)
  # agree, and the declared blocks change.
  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 2L))
  )
  h <- function(transform, value){
    switch(
      transform$type,
      identity = value,
      square = value^2,
      stop("unexpected coefficient transform")
    )
  }
  check_affine <- function(fit, plan){
    posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))[1L, , drop = FALSE]
    expect_identical(plan$family, "affine")
    covariance <- function(value){
      draw <- posterior
      draw[, plan$source_parameter] <- value
      random_effects_marginal_vcov(fit, posterior_samples = draw)$samples[1L, , ]
    }
    values <- c(0.15, 0.4, 0.8)
    C <- lapply(values, covariance)
    coefficient <- h(plan$coefficient_transform, values)
    first <- (C[[2L]] - C[[1L]]) / (coefficient[[2L]] - coefficient[[1L]])
    second <- (C[[3L]] - C[[1L]]) / (coefficient[[3L]] - coefficient[[1L]])
    expect_equal(first, second, tolerance = 1e-12)
    expect_gt(max(abs(first)), 0)
  }
  diagonal <- .random_update_test_fit(~ 1 + diag(1 + f | id), data)
  check_affine(diagonal, .random_update_test_plan(diagonal, "random_sd", "intercept"))
  hcs_fit <- .random_update_test_fit(~ 1 + hcs(f | id), data)
  check_affine(hcs_fit, .random_update_test_plan(hcs_fit, "random_correlation"))
})


test_that("random covariance updates are classified from formula metadata", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 2L))
  )
  diagonal <- .random_update_test_fit(~ 1 + diag(1 + f | id), data)
  hcs_fit  <- .random_update_test_fit(~ 1 + hcs(f | id), data)
  ar1_fit  <- .random_update_test_fit(~ 1 + ar1(f | id), data)
  ar1_two_fit <- .random_update_test_fit(
    ~ 1 + ar1(f | id),
    droplevels(data[data$f != "z", , drop = FALSE])
  )

  diagonal_sd <- .random_update_test_plan(
    diagonal,
    "random_sd",
    "intercept"
  )
  hcs_sd <- .random_update_test_plan(hcs_fit, "random_sd", "f[x]")
  hcs_cor <- .random_update_test_plan(hcs_fit, "random_correlation")
  ar1_cor <- .random_update_test_plan(ar1_fit, "random_correlation")
  ar1_two_cor <- .random_update_test_plan(
    ar1_two_fit,
    "random_correlation"
  )

  expect_identical(diagonal_sd$family, "affine")
  expect_identical(diagonal_sd$update, "column_scale")
  expect_identical(diagonal_sd$coefficient_input, "source")
  expect_identical(hcs_sd$family, "factor")
  expect_identical(hcs_cor$family, "affine")
  expect_identical(hcs_cor$update, "correlation")
  expect_identical(hcs_cor$coefficient_input, "source")
  expect_identical(ar1_cor$family, "markov")
  expect_identical(ar1_two_cor$family, "affine")
})


test_that("non-affine update grids reproduce compiled factor states", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 2L))
  )
  cases <- list(
    factor = list(
      fit = .random_update_test_fit(~ 1 + hcs(f | id), data),
      role = "random_sd",
      component = "f[x]",
      values = c(0.25, 0.8)
    ),
    markov = list(
      fit = .random_update_test_fit(~ 1 + ar1(f | id), data),
      role = "random_correlation",
      component = NULL,
      values = c(-0.25, 0.4)
    )
  )

  for(family in names(cases)){
    current <- cases[[family]]
    fit     <- current$fit
    plan <- .random_update_test_plan(
      fit,
      current$role,
      current$component
    )
    posterior <- as.matrix(fit)
    grid <- random_effects_marginal_update_grid(
      fit = fit,
      update = plan,
      values = current$values,
      posterior_samples = posterior
    )
    expect_s3_class(
      grid,
      "BayesTools_random_effects_marginal_update_grid"
    )
    expect_identical(grid$family, family)

    for(value_i in seq_along(current$values)){
      candidate <- posterior
      candidate[, plan$source_parameter] <- current$values[[value_i]]
      direct <- random_effects_marginal_factor_states(
        fit = fit,
        parameter = "mu",
        posterior_samples = candidate,
        prior_list = attr(fit, "prior_list"),
        blocks = plan$blocks,
        row_blocks = list(seq_len(nrow(data)))
      )
      for(draw in seq_len(nrow(posterior))){
        if(identical(family, "factor")){
          scale <- grid$coefficient_scale[draw, ]
          scale[[grid$component_index]] <-
            grid$candidate_scale[value_i, draw]
          cholesky <- matrix(
            grid$coefficient_cholesky[draw, , ],
            nrow = length(scale),
            ncol = length(scale)
          )
        }else{
          scale <- grid$coefficient_scale[draw, ]
          cholesky <- matrix(
            grid$candidate_cholesky[value_i, , ],
            nrow = length(scale),
            ncol = length(scale)
          )
          expect_equal(
            direct$factor_states[[draw]][[1L]]$markov_transition,
            grid$candidate_transition[value_i, ],
            tolerance = 1e-12
          )
          expect_equal(
            direct$factor_states[[draw]][[1L]]$markov_innovation_variance,
            grid$candidate_innovation_variance[value_i, ],
            tolerance = 1e-12
          )
        }
        expect_equal(
          direct$factor_states[[draw]][[1L]]$coefficient_factor,
          cholesky * scale,
          tolerance = 1e-12
        )
      }
    }
  }
})


test_that("every random covariance structure declares an exact update route", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 2L)),
    time = rep(1:3, 2L)
  )
  formulas <- list(
    id = ~ 1 + id(1 + f | id),
    diag = ~ 1 + diag(1 + f | id),
    us = ~ 1 + us(1 + f | id),
    cs = ~ 1 + cs(f | id),
    hcs = ~ 1 + hcs(f | id),
    ar1 = ~ 1 + ar1(f | id),
    har = ~ 1 + har(f | id),
    car = ~ 1 + car(time | id)
  )

  for(structure in names(formulas)){
    fit <- .random_update_test_fit(formulas[[structure]], data)
    catalog <- parameter_catalog(fit)
    quantities <- catalog$quantities[
      catalog$quantities$role %in% c("random_sd", "random_correlation"),
      ,
      drop = FALSE
    ]
    plans <- lapply(seq_len(nrow(quantities)), function(index){
      selection <- parameter_catalog_resolve(
        catalog,
        quantities$canonical_name[[index]],
        "mu"
      )
      random_effects_marginal_update_plan(fit, selection)
    })
    sd_plans <- plans[quantities$role == "random_sd"]
    expect_true(
      length(sd_plans) > 0L &&
        all(vapply(
          sd_plans,
          function(plan) plan$family %in% c("affine", "factor"),
          logical(1)
        )),
      info = structure
    )
    correlation_plans <- plans[quantities$role == "random_correlation"]
    if(identical(structure, "us")){
      expect_true(all(vapply(
        correlation_plans,
        function(plan) identical(plan$family, "unsupported"),
        logical(1)
      )))
    }else{
      expect_true(all(vapply(
        correlation_plans,
        function(plan) plan$family %in% c("affine", "markov"),
        logical(1)
      )), info = structure)
    }

    non_affine <- Filter(
      function(plan) plan$family %in% c("factor", "markov"),
      plans
    )
    for(update in non_affine){
      values <- if(identical(update$family, "factor")){
        c(0.25, 0.75)
      }else{
        c(0.2, 0.4)
      }
      grid <- random_effects_marginal_update_grid(
        fit = fit,
        update = update,
        values = values,
        posterior_samples = as.matrix(fit)
      )
      expect_identical(grid$family, update$family, info = structure)
      expect_identical(nrow(grid$coefficient_scale), 2L, info = structure)
    }
  }
})


test_that("single direct random scale exposes its exact row covariance basis", {

  data <- data.frame(id = factor(c("a", "a", "b", "c")))
  fit <- .random_update_test_fit(~ 1 + (1 | id), data)
  plan <- .random_update_test_plan(fit, "random_sd")
  var_plan <- .random_update_test_plan(fit, "random_var")

  expect_identical(plan$family, "affine")
  expect_identical(plan$coefficient_transform, list(type = "square"))
  expect_identical(var_plan$family, "affine")
  expect_identical(var_plan$coefficient_input, "source")
  expect_identical(var_plan$coefficient_transform, list(type = "square"))
  expect_equal(
    plan$invariant_covariance$update_covariance,
    outer(as.integer(data$id), as.integer(data$id), "==") * 1
  )
  expect_equal(plan$invariant_covariance$base_covariance, matrix(0, 4L, 4L))
})


test_that("variance allocations declare exact scalar covariance inputs", {

  fit <- .random_update_test_allocation_fit()
  total_sd <- .random_update_test_plan(fit, "random_sd_total")
  total_var <- .random_update_test_plan(fit, "random_var_total")
  proportion <- .random_update_test_plan(
    fit,
    "random_var_prop",
    "intercept"
  )
  component_sd <- .random_update_test_plan(
    fit,
    "random_sd",
    "intercept"
  )
  component_var <- .random_update_test_plan(
    fit,
    "random_var",
    "intercept"
  )

  expect_identical(total_sd$family, "affine")
  expect_identical(total_sd$coefficient_input, "source")
  expect_identical(total_sd$coefficient_transform, list(type = "square"))
  expect_identical(total_var$coefficient_transform, list(type = "square"))
  expect_identical(proportion$update, "allocation")
  expect_identical(proportion$coefficient_input, "source")
  expect_identical(proportion$coefficient_transform, list(type = "identity"))
  expect_identical(component_sd$update, "scale")
  expect_identical(component_sd$coefficient_input, "quantity")
  expect_identical(component_sd$coefficient_transform, list(type = "square"))
  expect_identical(component_var$update, "scale")
  expect_identical(component_var$coefficient_input, "quantity")
  expect_identical(component_var$coefficient_transform, list(type = "identity"))
})



test_that("correlated component allocations retain nonlinear covariances", {

  fit <- .random_update_test_allocation_fit(
    ~ 1 + har(index | id, name = "study")
  )
  plan <- .random_update_test_plan(fit, "random_var_prop", "index[1]")
  expect_identical(plan$family, "unsupported")
  expect_identical(plan$reason, "correlated_component_allocation")
  expect_identical(.random_update_test_plan(fit, "random_sd_total")$family,
                   "affine")

  posterior <- as.matrix(fit)[rep(1L, 3L), , drop = FALSE]
  weight <- "mu__xRE_ALLOCx_allocation__weight"
  values <- c(.1, .5, .9)
  posterior[, paste0(weight, "[1]")] <- values
  posterior[, paste0(weight, "[2]")] <- 1 - values
  covariance <- random_effects_marginal_vcov(
    fit,
    parameter = "mu",
    posterior_samples = posterior
  )$samples
  for(i in seq_along(values)){
    w <- values[[i]]
    block <- .5^2 * matrix(
      c(w, .5 * sqrt(w * (1 - w)), .5 * sqrt(w * (1 - w)), 1 - w),
      2L, 2L
    )
    expected <- matrix(0, 4L, 4L)
    expected[1:2, 1:2] <- block
    expected[3:4, 3:4] <- block
    expect_equal(unname(covariance[i, , ]), expected, tolerance = 1e-14)
  }
  expect_gt(abs(covariance[2L, 1L, 2L] -
                  mean(covariance[c(1L, 3L), 1L, 2L])), .01)

  # the allocation scale drives both columns of the block: the update is the
  # affine scale update on this correlated structure; the grid evaluates its
  # factor form
  component_sd <- .random_update_test_plan(fit, "random_sd", "index[1]")
  expect_identical(component_sd$family, "affine")
  expect_identical(component_sd$update, "scale")
  expect_identical(component_sd$coefficient_input, "quantity")
  expect_identical(component_sd$blocks, "study")
  component_sd <- .random_update_test_factor_form(fit, component_sd)
  grid <- random_effects_marginal_update_grid(
    fit = fit,
    update = component_sd,
    values = c(0.25, 0.75),
    posterior_samples = as.matrix(fit)
  )
  expect_identical(grid$family, "factor")
  expect_equal(grid$candidate_scale[, 1L], c(0.25, 0.75))
  expect_equal(grid$candidate_scale[, 2L], c(0.25, 0.75))

  component_var <- .random_update_test_plan(fit, "random_var", "index[1]")
  expect_identical(component_var$family, "affine")
  expect_identical(component_var$coefficient_transform, list(type = "identity"))
  component_var <- .random_update_test_factor_form(fit, component_var)
  variance_grid <- random_effects_marginal_update_grid(
    fit = fit,
    update = component_var,
    values = c(0.25, 0.75)^2
  )
  expect_identical(component_var$coefficient_input, "quantity")
  expect_identical(component_var$coefficient_scale_transform, list(type = "sqrt"))
  expect_equal(variance_grid$candidate_scale, grid$candidate_scale)

  # setting the SD of the first component to v moves the allocation scale
  # (tau = v / sqrt(w_1)) with the weights fixed, so both column SDs are
  # multiplied by v / SD_1 (here 0.5 sqrt(1/2) each)
  expect_identical(component_sd$held_fixed, "allocation_weights")
  expect_identical(component_var$held_fixed, "allocation_weights")
  current_sd <- 0.5 / sqrt(2)
  expect_equal(unname(grid$coefficient_scale[1L, ]), rep(current_sd, 2L))
  expect_equal(variance_grid$column_factor, grid$column_factor)
  for(i in seq_len(nrow(grid$candidate_scale))){
    factor <- grid$candidate_scale[i, 1L] / current_sd
    expect_equal(grid$column_factor[i, ], rep(factor, 2L), tolerance = 1e-14)
    scales <- current_sd * factor * c(1, 1)
    expected <- kronecker(diag(2), tcrossprod(scales) *
      matrix(c(1, 0.5, 0.5, 1), nrow = 2L))
    coefficient_factor <- diag(current_sd * grid$column_factor[i, ]) %*%
      variance_grid$coefficient_cholesky[1L, , ]
    expect_equal(kronecker(diag(2), tcrossprod(coefficient_factor)),
                 expected, tolerance = 1e-14)

    candidate <- as.matrix(fit)
    candidate[, "tau"] <- grid$candidate_scale[i, 1L] / sqrt(0.5)
    direct <- random_effects_marginal_vcov(
      fit, parameter = "mu", posterior_samples = candidate
    )$samples
    expect_equal(unname(direct[1L, , ]), expected, tolerance = 1e-14)
  }
  for(update in list(component_sd, component_var)){
    expect_error(
      random_effects_marginal_update_grid(fit, update, values = -0.2),
      "must be non-negative"
    )
    zero <- random_effects_marginal_update_grid(fit, update, values = 0)
    expect_equal(zero$candidate_scale, matrix(0, nrow = 1L, ncol = 2L))
  }
})


test_that("nested aggregate allocations use their public covariance scale", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    paper = factor(c("p1", "p2", "p1", "p2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        weights = prior("dirichlet", list(alpha = c(3, 2)))
      )
    )
  )
  allocation_columns <- function(name){
    weight <- paste0("mu__xRE_ALLOCx_", name, "__weight")
    c(
      paste0(weight, "[", 1:2, "]"),
      paste0("prior_par_eta_", weight, "[", 1:2, "]")
    )
  }
  columns <- c(
    "mu_intercept",
    "mu__xRE_ALLOCx_total_re__allocation_sd",
    allocation_columns("total_re"),
    allocation_columns("nested_split")
  )
  draws <- matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)

  nested_sd <- .random_update_test_plan(
    fit,
    "random_sd_total",
    "nested_split"
  )
  nested_var <- .random_update_test_plan(
    fit,
    "random_var_total",
    "nested_split"
  )
  outer_sd <- .random_update_test_plan(
    fit,
    "random_sd_total",
    "total_re"
  )

  expect_identical(outer_sd$blocks, c("study", "paper", "drug"))
  # the nested aggregate moves the root scale with the weights fixed, which
  # reaches the other component of the root as well
  expect_identical(nested_sd$blocks, c("study", "paper", "drug"))
  expect_identical(nested_sd$coefficient_input, "quantity")
  expect_identical(nested_sd$held_fixed, "allocation_weights")
  expect_identical(nested_sd$coefficient_transform, list(type = "square"))
  expect_identical(nested_var$coefficient_input, "quantity")
  expect_identical(nested_var$coefficient_transform, list(type = "identity"))
})


test_that("shared independent-coefficient scale exposes an exact dense basis", {

  data <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    x = c(-1, 0, 0.5, 2)
  )
  fit <- .random_update_test_fit(~ 1 + id(1 + x | id), data)
  plan <- .random_update_test_plan(fit, "random_sd")
  design <- attr(fit, "formula_design")$mu$random_effects[[1L]]
  same_group <- outer(as.integer(data$id), as.integer(data$id), "==") * 1
  expected <- same_group * tcrossprod(design$model_matrix)

  expect_identical(plan$family, "affine")
  expect_equal(plan$invariant_covariance$update_covariance, unname(expected))
})


test_that("the multi-column update plan checks its dense covariance against the memory ceiling", {

  # F286: the multi-column branch forms rows x rows copies (group kernel,
  # design cross-product, product) that nothing bounded
  option <- "BayesTools.random_effects_memory_limit_bytes"
  old_limit <- getOption(option)
  on.exit(do.call(options, stats::setNames(list(old_limit), option)), add = TRUE)
  options(stats::setNames(list(Inf), option))

  data <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    x = c(-1, 0, 0.5, 2)
  )
  fit <- .random_update_test_fit(~ 1 + id(1 + x | id), data)
  reference <- .random_update_test_plan(fit, "random_sd")
  expect_identical(reference$family, "affine")

  # 4 rows, one block: the kept block matrix (8 * 16) and the working copies (36 * 16)
  peak <- 8 * 16 + 36 * 16
  options(stats::setNames(list(peak - 1), option))
  condition <- tryCatch(.random_update_test_plan(fit, "random_sd"), error = identity)
  expect_s3_class(condition, "BayesTools_random_effects_memory_limit")
  expect_s3_class(condition, "BayesTools_random_effects_error")
  expect_identical(condition$operation, "random-effect marginal update plan")
  expect_identical(condition$peak_bytes, peak)
  expect_identical(condition$limit_bytes, peak - 1)
  expect_match(
    conditionMessage(condition),
    paste0("Estimated peak memory for random-effect marginal update plan for block 'id' ",
           "exceeds option 'BayesTools.random_effects_memory_limit_bytes'. ",
           "Dimensions: 4 rows x 4 rows; 1 block;"),
    fixed = TRUE
  )
  expect_match(conditionMessage(condition), "fit the model to fewer rows", fixed = TRUE)

  # at the estimate the plan is unchanged
  options(stats::setNames(list(peak), option))
  expect_identical(.random_update_test_plan(fit, "random_sd"), reference)

  estimate <- BayesTools:::.bt_random_effect_update_plan_memory_estimate(4L, 1L)
  expect_identical(
    estimate$components,
    c(block_covariances = 8 * 16, assembly_working = 36 * 16)
  )
  expect_identical(estimate$peak_bytes, peak)
  # every block keeps its own matrix until the blocks are summed
  three <- BayesTools:::.bt_random_effect_update_plan_memory_estimate(4L, 3L)
  expect_identical(three$components[["block_covariances"]], 3 * 8 * 16)
  expect_match(three$dimensions, "3 blocks", fixed = TRUE)
  expect_identical(
    BayesTools:::.bt_random_effect_update_plan_memory_estimate(1e200, 2L)$peak_bytes,
    Inf
  )

  # the guard is the shared one: its other sites carry the class as well
  options(stats::setNames(list(100), option))
  expect_error(
    BayesTools:::.bt_random_effect_check_memory(estimate, alternative = "Use fewer rows."),
    class = "BayesTools_random_effects_memory_limit"
  )
})


test_that("known group covariance remains an exact affine scale basis", {

  data <- data.frame(id = factor(c("b", "a", "c", "b")))
  kernel <- matrix(
    c(4, 1, 0.5, 1, 9, 2, 0.5, 2, 16),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  formula <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(kernel, scale = "none")
  )
  fit <- .random_update_test_fit(formula, data)
  plan <- .random_update_test_plan(fit, "random_sd")
  expected <- kernel[as.character(data$id), as.character(data$id)]

  expect_identical(plan$family, "affine")
  expect_equal(
    plan$invariant_covariance$update_covariance,
    unname(expected)
  )
})


test_that("LKJ correlation updates distinguish scalar and composite paths", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    x = c(-1, 0, 1, -0.5, 0.5, 2),
    z = c(0, 1, 0, 1, 0, 1)
  )
  two <- .random_update_test_fit(~ 1 + us(1 + x | id), data)
  three <- .random_update_test_fit(~ 1 + us(1 + x + z | id), data)
  two_plan <- .random_update_test_plan(two, "random_correlation")

  catalog <- parameter_catalog(three)
  correlations <- catalog$quantities[
    catalog$quantities$role == "random_correlation",
    ,
    drop = FALSE
  ]
  selection <- parameter_catalog_resolve(
    catalog,
    correlations$canonical_name[[1L]],
    "mu"
  )
  three_plan <- random_effects_marginal_update_plan(three, selection)

  expect_identical(two_plan$family, "affine")
  expect_identical(three_plan$family, "unsupported")
  expect_identical(three_plan$reason, "non_scalar_covariance_path")
})


.random_update_test_block_allocation_fit <- function(gates = "shared"){

  data <- data.frame(
    study = factor(c("a", "a", "b", "b")),
    esid = factor(seq_len(4L))
  )
  if(identical(gates, "shared")){
    allocations <- list(
      random_variance_allocation(
        name = "root",
        terms = c(component = "component_gated"),
        sd = prior("gamma", list(2, 2)),
        inclusion = list(component = prior("spike", list(location = 0.5)))
      ),
      random_variance_allocation(
        name = "split",
        terms = c(study = "study", esid = "esid"),
        parent = allocation_ref("root", "component"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  }else{
    allocations <- random_variance_allocation(
      name = "split",
      terms = c(study = "study", esid = "esid"),
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(1, 1))),
      inclusion = if(identical(gates, "independent")) list(
        study = prior("spike", list(location = 0.5)),
        esid = prior("spike", list(location = 0.5))
      ) else NULL
    )
  }
  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(allocation = allocations)
  )
  allocation <- result$formula_design$random_allocations$split
  source <- allocation$source$name
  weight <- allocation$weight_name
  indicators <- .bt_random_effect_summary_allocation_gate_names(allocation)
  columns <- c("mu_intercept", source, paste0(weight, "[", 1:2, "]"),
               indicators)
  samples <- matrix(
    0.5, nrow = 2L, ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  samples[, source] <- c(0.7, 0.9)
  samples[, paste0(weight, "[1]")] <- c(0.25, 0.4)
  samples[, paste0(weight, "[2]")] <- c(0.75, 0.6)
  if(length(indicators) > 0L){
    samples[, indicators] <- 1
    samples[1L, indicators] <- 0
  }
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)
  list(fit = fit, samples = samples, data = data, source = source,
       weight = weight, indicators = indicators)
}


test_that("shared parent gates declare conditional block-allocation updates", {

  fixture <- .random_update_test_block_allocation_fit()
  selectors <- c("(mu) split: sd_total", "(mu) esid: sd(intercept)",
                 "(mu) split: var_prop(esid)")
  coefficient_inputs <- c("source", "quantity", "source")
  transforms <- c("square", "square", "identity")
  values <- c(0.2, 0.5, 0.8)
  study_covariance <- outer(fixture$data$study, fixture$data$study, "==") * 1
  estimate_covariance <- diag(nrow(fixture$data))

  for(i in seq_along(selectors)){
    selection <- parameter_catalog_resolve(parameter_catalog(fixture$fit),
                                           selectors[[i]])
    update <- random_effects_marginal_update_plan(fixture$fit, selection)
    expect_identical(update$family, "unsupported")
    expect_identical(update$conditional$required_active, fixture$indicators)
    plan <- update$conditional$plan
    expect_s3_class(plan, "BayesTools_random_effects_marginal_update_plan")
    expect_identical(plan$family, "affine")
    expect_identical(plan$coefficient_input, coefficient_inputs[[i]])
    expect_identical(plan$coefficient_transform, list(type = transforms[[i]]))
    expect_setequal(plan$blocks, c("study", "esid"))
    expect_identical(plan$source_parameter,
                     if(i == 3L) fixture$weight else fixture$source)
    expect_identical(plan$allocation$weight_name, fixture$weight)
    expect_identical(plan$allocation$n_targets, 2L)

    samples <- fixture$samples[rep(2L, length(values)), , drop = FALSE]
    source <- samples[1L, fixture$source]
    proportion <- samples[1L, paste0(fixture$weight, "[2]")]
    if(i == 1L){
      samples[, fixture$source] <- values
      expected <- lapply(values, function(value){
        value^2 * ((1 - proportion) * study_covariance +
                     proportion * estimate_covariance)
      })
    }else if(i == 2L){
      samples[, fixture$source] <- values / sqrt(proportion)
      expected <- lapply(values, function(value){
        value^2 * ((1 - proportion) / proportion * study_covariance +
                     estimate_covariance)
      })
    }else{
      samples[, paste0(fixture$weight, "[1]")] <- 1 - values
      samples[, paste0(fixture$weight, "[2]")] <- values
      expected <- lapply(values, function(value){
        source^2 * ((1 - value) * study_covariance +
                      value * estimate_covariance)
      })
    }
    covariance <- random_effects_marginal_vcov(
      fixture$fit, parameter = "mu", posterior_samples = samples
    )$samples
    for(row in seq_along(values)){
      expect_equal(unname(covariance[row, , ]), unname(expected[[row]]),
                   tolerance = 1e-12)
    }
    samples[, fixture$indicators] <- 0
    inactive <- random_effects_marginal_vcov(
      fixture$fit, parameter = "mu", posterior_samples = samples
    )$samples
    expect_equal(unname(inactive), array(0, dim = dim(inactive)))
  }
})


test_that("block component scales resolve uniquely without enabling independent gates", {

  fixture <- .random_update_test_block_allocation_fit("none")
  selection <- parameter_catalog_resolve(
    parameter_catalog(fixture$fit), "(mu) esid: sd(intercept)"
  )
  plan <- random_effects_marginal_update_plan(fixture$fit, selection)
  expect_identical(plan$family, "affine")
  expect_null(plan$conditional)
  expect_identical(plan$coefficient_input, "quantity")
  expect_identical(plan$coefficient_transform, list(type = "square"))
  expect_identical(plan$source_parameter, fixture$source)
  expect_identical(plan$allocation$index, 2L)

  independent <- .random_update_test_block_allocation_fit("independent")
  for(selector in c("(mu) split: sd_total", "(mu) esid: sd(intercept)",
                    "(mu) split: var_prop(esid)")){
    selection <- parameter_catalog_resolve(parameter_catalog(independent$fit),
                                           selector)
    plan <- random_effects_marginal_update_plan(independent$fit, selection)
    expect_identical(plan$family, "unsupported")
    expect_null(plan$conditional)
  }
})

# Allocation fit of one random block with distinct weights and a non-zero
# correlation, for the covariance structures that resolve several SD components.
.random_update_test_structure_fit <- function(structure){

  data <- data.frame(
    id = factor(rep(c("a", "b", "c"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 3L)),
    x = rep(c(-1, 0.5, 2), 3L)
  )
  formula <- switch(
    structure,
    diag = ~ 1 + random(1 + x | id, name = "study", covariance = "diag"),
    us = ~ 1 + random(1 + x | id, name = "study", covariance = "us"),
    hcs = ~ 1 + hcs(f | id, name = "study"),
    har = ~ 1 + har(f | id, name = "study")
  )
  weights <- if(structure %in% c("diag", "us")) c(0.36, 0.64) else c(0.2, 0.3, 0.5)
  K <- length(weights)
  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "total_variance",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = rep(1, K)))
      )
    )
  )
  term <- result$formula_design$random_effects[[1L]]
  weight <- "mu__xRE_ALLOCx_allocation__weight"
  correlation <- term$correlation
  correlation_values <- if(is.null(correlation)){
    numeric()
  }else if(identical(correlation$type, "lkj")){
    stats::setNames(rep(0.8, length(correlation$primitive_names)), correlation$primitive_names)
  }else{
    stats::setNames(atanh(0.6), correlation$sample_name)
  }
  values <- c(
    mu_intercept = 0,
    tau = 0.5,
    stats::setNames(0.5 * sqrt(weights), unique(term$sd_parameter_names)),
    stats::setNames(weights, paste0(weight, "[", seq_len(K), "]")),
    stats::setNames(rep(1, K), paste0("prior_par_eta_", weight, "[", seq_len(K), "]")),
    correlation_values
  )
  draws <- matrix(rep(values, 2L), nrow = 2L, byrow = TRUE,
                  dimnames = list(NULL, names(values)))
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = colnames(draws),
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_draw_geometry(fit)

  list(fit = .bt_attach_fit_contract(fit), weights = weights, tau = 0.5)
}

test_that("an allocation-derived SD update holds the weights fixed and moves the allocation scale", {

  # F092: the same candidate SD implies the same covariance on every
  # structure: tau' = v / sqrt(w_k) with the weights and the correlation fixed
  v <- 0.25
  for(structure in c("diag", "us", "hcs", "har")){
    case <- .random_update_test_structure_fit(structure)
    fit <- case$fit
    weights <- case$weights
    current_sd <- case$tau * sqrt(weights[1L])
    posterior <- as.matrix(fit)
    catalog <- parameter_catalog(fit)
    sd_rows <- catalog$quantities[catalog$quantities$role == "random_sd", , drop = FALSE]
    plan <- random_effects_marginal_update_plan(
      fit, parameter_catalog_resolve(catalog, sd_rows$canonical_name[[1L]], "mu")
    )
    variance_plan <- random_effects_marginal_update_plan(
      fit, parameter_catalog_resolve(
        catalog, catalog$quantities$canonical_name[
          catalog$quantities$role == "random_var"][[1L]], "mu")
    )
    expect_identical(plan$held_fixed, "allocation_weights", info = structure)
    expect_identical(variance_plan$held_fixed, "allocation_weights", info = structure)
    expect_identical(plan$coefficient_input, "quantity", info = structure)
    # the allocation scale drives every column of the single block: the update
    # is affine in the selected SD on every structure
    expect_identical(plan$family, "affine", info = structure)
    expect_identical(plan$blocks, "study", info = structure)

    # the covariance recomputed from (tau' = v / sqrt(w_1), w)
    candidate <- posterior
    candidate[, "tau"] <- v / sqrt(weights[1L])
    expected <- random_effects_marginal_vcov(
      fit, parameter = "mu", posterior_samples = candidate
    )$samples[1L, , ]
    current <- random_effects_marginal_vcov(
      fit, parameter = "mu", posterior_samples = posterior
    )$samples[1L, , ]
    # the sibling SDs move with the selected one
    expect_equal(unname(expected), unname(current) * (v / current_sd)^2,
                 tolerance = 1e-12, info = structure)

    # A + h(v) B with h the square: the block contribution scales by c^2
    expect_identical(plan$coefficient_transform, list(type = "square"), info = structure)
    if(identical(structure, "diag")){
      next
    }

    # the factor form of the same update (a block that the scale drives only in
    # part), evaluated through the grid with every column driven
    plan <- .random_update_test_factor_form(fit, plan)
    grid <- random_effects_marginal_update_grid(
      fit, plan, values = c(v, 0), posterior_samples = posterior
    )
    expect_identical(grid$held_fixed, "allocation_weights", info = structure)
    expect_identical(grid$driven_columns, seq_len(length(weights)), info = structure)
    expect_equal(grid$candidate_scale[, 1L], c(v, 0), info = structure)
    expect_equal(unname(grid$coefficient_scale[1L, ]), case$tau * sqrt(weights),
                 tolerance = 1e-12, info = structure)
    expect_equal(grid$column_factor[, 1L], c(v, 0) / current_sd,
                 tolerance = 1e-12, info = structure)
    term <- attr(fit, "formula_design")$mu$random_effects[[1L]]
    for(value_i in 1:2){
      scale <- grid$coefficient_scale[1L, ] * grid$column_factor[value_i, 1L]
      factor <- matrix(grid$coefficient_cholesky[1L, , ], length(scale)) * scale
      covariance <- outer(term$group_map, term$group_map, "==") *
        (term$model_matrix %*% tcrossprod(factor) %*% t(term$model_matrix))
      reference <- if(value_i == 1L) expected else {
        candidate[, "tau"] <- 0
        random_effects_marginal_vcov(
          fit, parameter = "mu", posterior_samples = candidate
        )$samples[1L, , ]
      }
      expect_equal(unname(covariance), unname(reference), tolerance = 1e-12,
                   info = paste(structure, "value", value_i))
    }

    # a draw whose selected SD is zero has no finite factor
    zero <- posterior
    zero[, paste0("mu__xRE_ALLOCx_allocation__weight[1]")] <- 0
    zero[, "mu__xRE_ALLOCx_allocation__weight[2]"] <-
      zero[, "mu__xRE_ALLOCx_allocation__weight[2]"] + weights[1L]
    zero[, unique(term$sd_parameter_names)[1L]] <- 0
    zero_grid <- random_effects_marginal_update_grid(
      fit, plan, values = c(v, 0), posterior_samples = zero
    )
    expect_true(all(is.na(zero_grid$column_factor)), info = structure)
    expect_equal(zero_grid$candidate_scale[, 1L], c(v, 0), info = structure)
  }

  # a direct per-column SD holds its siblings fixed and has no held_fixed field
  direct <- .random_update_test_fit(
    ~ 1 + hcs(f | id), data.frame(id = factor(rep(c("a", "b"), each = 3L)),
                                  f = factor(rep(c("x", "y", "z"), 2L)))
  )
  direct_plan <- .random_update_test_plan(direct, "random_sd", "f[x]")
  expect_null(direct_plan$held_fixed)
  expect_null(random_effects_marginal_update_grid(
    direct, direct_plan, values = 0.3, posterior_samples = as.matrix(direct)
  )$column_factor)
})


# Split-root layout of RoBMA: a root allocation over two blocks and one
# `<block>_components` SD-component allocation per block, plus a block `site`
# with its own SD that the root scale does not reach. The values give every
# component a distinct SD, SD_k = tau * sqrt(w_block) * sqrt(K u_k).
.random_update_test_split_root_fit <- function(study_structure){

  data <- data.frame(
    id = factor(rep(c("a", "b", "c"), each = 3L)),
    other = factor(rep(c("u", "v", "w"), 3L)),
    site = factor(rep(c("k", "l", "m"), 3L)),
    x = rep(c(-1, 0.5, 2), 3L)
  )
  formula <- if(identical(study_structure, "us")){
    ~ 1 + random(1 + x | id, name = "study", covariance = "us") +
      random(1 + x | other, name = "paper", covariance = "diag") +
      random(1 | site, name = "site", covariance = "diag")
  }else{
    ~ 1 + random(1 + x | id, name = "study", covariance = "diag") +
      random(1 + x | other, name = "paper", covariance = "diag") +
      random(1 | site, name = "site", covariance = "diag")
  }
  dirichlet <- prior("dirichlet", list(alpha = c(1, 1)))
  prior_random <- prior_random(site = random_block(sd = prior("gamma", list(2, 2))), allocation = list(
    random_variance_allocation(
      name = "root", terms = c(s = "study", p = "paper"),
      sd = prior("gamma", list(2, 2)), weights = dirichlet
    ),
    random_variance_allocation(
      name = "study_components", terms = "study",
      parent = allocation_ref("root", "s"), weights = dirichlet,
      target = "sd_component", scale = "mean_variance"
    ),
    random_variance_allocation(
      name = "paper_components", terms = "paper",
      parent = allocation_ref("root", "p"), weights = dirichlet,
      target = "sd_component", scale = "mean_variance"
    )
  ))
  tau <- 0.8
  root <- c(0.35, 0.65)
  study <- c(0.3, 0.7)
  paper <- c(0.6, 0.4)
  weights <- function(name, w){
    stats::setNames(
      c(w, rep(1, length(w))),
      paste0(
        rep(c("", "prior_par_eta_"), each = length(w)),
        "mu__xRE_ALLOCx_", name, "__weight[", seq_along(w), "]"
      )
    )
  }
  values <- c(
    stats::setNames(tau, "mu__xRE_ALLOCx_root__allocation_sd"),
    weights("root", root),
    weights("study_components", study),
    weights("paper_components", paper)
  )
  fit <- .random_update_test_generic_fit(
    formula, data, prior_random, list(intercept = prior("normal", list(0, 1))),
    values = if(identical(study_structure, "us")){
      c(values, "mu__xREx__study_xRE_CORx_lkj_u[1]" = 0.8)
    }else{
      values
    }
  )
  list(
    fit = fit,
    source = "mu__xRE_ALLOCx_root__allocation_sd",
    tau = tau,
    sd = list(
      study = tau * sqrt(root[[1L]]) * sqrt(2 * study),
      paper = tau * sqrt(root[[2L]]) * sqrt(2 * paper)
    )
  )
}


test_that("an allocation-derived SD update lists every block that the moved scale reaches", {

  # F092: the moved scale is the root allocation scale, which every block of the
  # split-root layout derives from; with the weights fixed the covariance of all
  # reached blocks scales by c^2 = (v / SD_k)^2, on any covariance structure
  for(study_structure in c("diag", "us")){
    case <- .random_update_test_split_root_fit(study_structure)
    fit <- case$fit
    posterior <- as.matrix(fit)
    catalog <- parameter_catalog(fit)
    vcov_of <- function(samples, blocks = NULL){
      random_effects_marginal_vcov(
        fit, parameter = "mu", blocks = blocks, posterior_samples = samples
      )$samples[1L, , ]
    }
    current <- vcov_of(posterior)
    cases <- list(
      list("(mu) study: sd(intercept)", "study", 1L, FALSE),
      list("(mu) study: sd(x)", "study", 2L, FALSE),
      list("(mu) paper: sd(intercept)", "paper", 1L, FALSE),
      list("(mu) paper: sd(x)", "paper", 2L, FALSE),
      list("(mu) study: var(x)", "study", 2L, TRUE),
      list("(mu) paper: var(intercept)", "paper", 1L, TRUE)
    )
    for(case_i in cases){
      info <- paste(study_structure, case_i[[1L]])
      block <- case_i[[2L]]
      current_sd <- case$sd[[block]][[case_i[[3L]]]]
      plan <- random_effects_marginal_update_plan(
        fit, parameter_catalog_resolve(catalog, case_i[[1L]])
      )
      expect_identical(plan$family, "affine", info = info)
      expect_identical(plan$update, "scale", info = info)
      expect_identical(plan$held_fixed, "allocation_weights", info = info)
      expect_identical(plan$coefficient_input, "quantity", info = info)
      expect_identical(
        plan$coefficient_transform,
        if(case_i[[4L]]) list(type = "identity") else list(type = "square"),
        info = info
      )
      # the selected block first, then every block that the scale reaches
      expect_identical(plan$blocks, c(block, setdiff(c("study", "paper"), block)),
                       info = info)

      # candidate SD value v (variance v^2): the moved scale is tau' = tau v / SD_k
      v <- 0.5 * current_sd
      candidate <- posterior
      candidate[, case$source] <- case$tau * v / current_sd
      moved <- vcov_of(candidate)
      c2 <- (v / current_sd)^2
      expect_equal(moved, current + (c2 - 1) * vcov_of(posterior, plan$blocks),
                   tolerance = 1e-12, info = info)
      for(reached in plan$blocks){
        expect_equal(
          vcov_of(candidate, reached), vcov_of(posterior, reached) * c2,
          tolerance = 1e-12, info = paste(info, reached)
        )
      }
      # the block that does not use the scale stays
      expect_equal(vcov_of(candidate, "site"), vcov_of(posterior, "site"),
                   tolerance = 1e-12, info = info)
      # the selected block alone does not reproduce the update
      only_block <- current + (vcov_of(posterior, block)) * (c2 - 1)
      expect_gt(max(abs(only_block - moved)), 1e-3)
    }
  }
})


test_that("a scale shared by several blocks reaches all of them in an allocation-derived SD update", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    other = factor(c("a", "b", "a", "b")),
    third = factor(c("u", "v", "v", "u")),
    x = c(0, 1, 0, 1)
  )
  tau <- random_sd_source("tau")
  allocation <- function(name, term){
    random_variance_allocation(
      name = name, terms = term, target = "sd_component", sd_source = tau,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    )
  }
  plan_for <- function(fit, quantity){
    random_effects_marginal_update_plan(
      fit, parameter_catalog_resolve(parameter_catalog(fit), quantity)
    )
  }
  values <- c(
    "tau" = 0.9,
    "mu__xRE_ALLOCx_alloc1__weight[1]" = 0.25,
    "mu__xRE_ALLOCx_alloc1__weight[2]" = 0.75,
    "mu__xRE_ALLOCx_alloc2__weight[1]" = 0.6,
    "mu__xRE_ALLOCx_alloc2__weight[2]" = 0.4
  )
  for(structure in c("diag", "us")){
    # two blocks bound to random_sd_source("tau") through SD-component allocations
    # and a third block bound to it directly
    fit <- .random_update_test_generic_fit(
      stats::as.formula(paste0(
        "~ 1 + random(1 + x | id, name = 'first', covariance = '", structure, "') +",
        "random(1 + x | other, name = 'second', covariance = 'diag') +",
        "random(1 | third, name = 'third', covariance = 'diag')"
      )),
      df,
      prior_random(
        third = random_block(sd_source = tau),
        allocation = list(allocation("alloc1", "first"), allocation("alloc2", "second"))
      ),
      list(intercept = prior("normal", list(0, 1))),
      extra_columns = "tau",
      values = values
    )
    posterior <- as.matrix(fit)
    catalog <- parameter_catalog(fit)
    vcov_of <- function(samples, blocks = NULL){
      random_effects_marginal_vcov(
        fit, parameter = "mu", blocks = blocks, posterior_samples = samples
      )$samples[1L, , ]
    }
    current <- vcov_of(posterior)
    for(quantity in c("(mu) first: sd(x)", "(mu) second: sd(intercept)", "(mu) first: var(intercept)")){
      info <- paste(structure, quantity)
      plan <- plan_for(fit, quantity)
      expect_identical(plan$family, "affine", info = info)
      expect_identical(plan$held_fixed, "allocation_weights", info = info)
      expect_identical(plan$coefficient_input, "quantity", info = info)
      expect_setequal(plan$blocks, c("first", "second", "third"))
      expect_identical(plan$blocks[[1L]], sub("^\\(mu\\) (.*): .*$", "\\1", quantity), info = info)

      # SD_k = tau sqrt(w_k) (total variance) of the selected component
      block <- plan$blocks[[1L]]
      column <- if(grepl("intercept", quantity, fixed = TRUE)) 1L else 2L
      current_sd <- 0.9 * sqrt(values[[paste0(
        "mu__xRE_ALLOCx_", if(identical(block, "first")) "alloc1" else "alloc2",
        "__weight[", column, "]"
      )]])
      v <- 1.3 * current_sd
      candidate <- posterior
      candidate[, "tau"] <- 0.9 * v / current_sd
      expect_equal(vcov_of(candidate), current * (v / current_sd)^2,
                   tolerance = 1e-12, info = info)
      expect_equal(
        vcov_of(candidate, plan$blocks), vcov_of(posterior, plan$blocks) * (v / current_sd)^2,
        tolerance = 1e-12, info = info
      )
    }
  }
})


test_that("a scale that reaches several blocks only in part has no allocation-derived SD update", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    other = factor(c("a", "b", "a", "b")),
    x = c(0, 1, 0, 1)
  )
  tau <- random_sd_source("tau")
  allocation <- function(name, term){
    random_variance_allocation(
      name = name, terms = term, target = "sd_component", sd_source = tau,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    )
  }
  fit <- .random_update_test_generic_fit(
    ~ 1 + random(1 + x | id, name = "first", covariance = "diag") +
      random(1 + x | other, name = "second", covariance = "diag"),
    df,
    prior_random(allocation = list(allocation("alloc1", "first"), allocation("alloc2", "second"))),
    list(intercept = prior("normal", list(0, 1))),
    extra_columns = "tau"
  )
  plan_for <- function(fit, quantity){
    random_effects_marginal_update_plan(
      fit, parameter_catalog_resolve(parameter_catalog(fit), quantity)
    )
  }
  with_binding <- function(fit, block, binding){
    design <- attr(fit, "formula_design")
    terms <- design$mu$random_effects
    index <- which(vapply(terms, function(term) identical(term$block_name, block), logical(1)))
    design$mu$random_effects[[index]]$sd_binding <- binding
    attr(fit, "formula_design") <- design
    fit
  }
  owned <- function(name) BayesTools:::.bt_prior_owned_sd_source(name)
  column_binding <- function(sources){
    BayesTools:::.bt_random_sd_binding(
      sources_by_column = sources, application = "column",
      factors_by_column = rep(list(list()), length(sources))
    )
  }

  # the second block uses the scale for its intercept only: the scale reaches
  # two blocks and only part of the second
  partial <- with_binding(
    fit, "second", column_binding(list(tau, owned("second_other")))
  )
  refused <- plan_for(partial, "(mu) first: sd(x)")
  expect_identical(refused$family, "unsupported")
  expect_identical(refused$reason, "partial_scale_reach")
  expect_identical(refused$blocks, character())
  expect_match(refused$message, "only part", fixed = TRUE)

  # the same layout without the partial block is a full reach of both blocks
  expect_identical(plan_for(fit, "(mu) first: sd(x)")$family, "affine")

  # the scale reaches only the selected block and drives only some of its
  # columns: the factor family with the driven columns
  alone <- with_binding(
    with_binding(fit, "second", column_binding(list(owned("second_a"), owned("second_b")))),
    "first", {
      binding <- attr(fit, "formula_design")$mu$random_effects[[1L]]$sd_binding
      binding$sources_by_column <- list(tau, owned("first_other"))
      binding
    }
  )
  factor_plan <- plan_for(alone, "(mu) first: sd(intercept)")
  expect_identical(factor_plan$family, "factor")
  expect_identical(factor_plan$update, "column_scale")
  expect_identical(factor_plan$blocks, "first")
  expect_identical(factor_plan$held_fixed, "allocation_weights")
  expect_identical(factor_plan$driven_columns, 1L)
  expect_identical(factor_plan$component_index, 1L)
  expect_identical(factor_plan$coefficient_input, "quantity")

  # a selected column that the scale does not drive has no factor update: the
  # scale would not move it
  undriven <- plan_for(alone, "(mu) first: sd(x)")
  expect_identical(undriven$family, "unsupported")
  expect_identical(undriven$reason, "selected_column_not_driven")
  expect_identical(undriven$blocks, character())
  expect_match(undriven$message, "does not drive", fixed = TRUE)
  expect_identical(plan_for(alone, "(mu) first: var(x)")$reason, "selected_column_not_driven")

  # an aggregate has no selected column: a partly driven block refuses it
  aggregate <- plan_for(partial, "(mu) alloc1: sd_total")
  expect_identical(aggregate$family, "unsupported")
  expect_identical(aggregate$reason, "partial_scale_reach")
})


# The plan of a quantity that is the allocation scale times weights, against the
# covariance recomputed from the moved scale with every weight fixed: the blocks
# of the plan scale by the squared factor of the quantity, the others stay.
.random_update_test_scale_reach <- function(fit, name, source, factor, unreached,
                                            info = name){

  posterior <- as.matrix(fit)
  selection <- parameter_catalog_resolve(parameter_catalog(fit), name)
  plan <- random_effects_marginal_update_plan(fit, selection)
  vcov_of <- function(samples, blocks = NULL){
    random_effects_marginal_vcov(
      fit, parameter = "mu", blocks = blocks, posterior_samples = samples
    )$samples[1L, , ]
  }
  all_blocks <- vapply(
    attr(fit, "formula_design")$mu$random_effects,
    function(term) term$block_name, character(1)
  )
  candidate <- posterior
  candidate[, source] <- posterior[, source] * factor

  # the selected quantity is proportional to the scale (an SD) or to its square
  variance <- grepl(": var", name, fixed = TRUE)
  current <- as.numeric(as.matrix(parameter_draws(fit, selection)))[1L]
  moved_value <- as.numeric(as.matrix(
    parameter_draws(fit, selection, model_samples = candidate)
  ))[1L]
  expect_equal(moved_value, current * if(variance) factor^2 else factor,
               tolerance = 1e-12, info = info)

  expect_identical(plan$family, "affine", info = info)
  expect_identical(plan$update, "scale", info = info)
  expect_identical(sort(setdiff(all_blocks, plan$blocks)), sort(unreached), info = info)
  expect_equal(
    vcov_of(candidate),
    vcov_of(posterior) + (factor^2 - 1) * vcov_of(posterior, plan$blocks),
    tolerance = 1e-12, info = info
  )
  for(reached in plan$blocks){
    expect_equal(
      vcov_of(candidate, reached), vcov_of(posterior, reached) * factor^2,
      tolerance = 1e-12, info = paste(info, reached)
    )
  }
  if(length(unreached) > 0L){
    expect_equal(vcov_of(candidate, unreached), vcov_of(posterior, unreached),
                 tolerance = 1e-12, info = info)
  }

  plan
}


test_that("the aggregates of nested allocations move the root scale and list every block that uses it", {

  # F092: sd_common/var_common of RoBMA's `<block>_components` allocations and
  # the totals of a nested block split are the root scale times weights; moving
  # them with the weights fixed moves the root scale, which reaches the sibling
  # block, and not the block `site` that has its own SD
  for(study_structure in c("diag", "us")){
    case <- .random_update_test_split_root_fit(study_structure)
    for(allocation in c("study_components", "paper_components")){
      for(quantity in c("sd_common", "var_common")){
        name <- paste0("(mu) ", allocation, ": ", quantity)
        info <- paste(study_structure, name)
        plan <- .random_update_test_scale_reach(
          case$fit, name, case$source, factor = 0.5, unreached = "site", info = info
        )
        block <- sub("_components$", "", allocation)
        expect_identical(plan$blocks, c(block, setdiff(c("study", "paper"), block)), info = info)
        expect_identical(plan$coefficient_input, "quantity", info = info)
        expect_identical(plan$held_fixed, "allocation_weights", info = info)
        expect_identical(
          plan$coefficient_transform,
          if(identical(quantity, "var_common")) list(type = "identity") else list(type = "square"),
          info = info
        )
      }
    }
  }

  # a block split nested in a root allocation, and the root itself
  data <- data.frame(
    id = factor(rep(c("a", "b", "c"), each = 4L)),
    other = factor(rep(c("u", "v", "w", "z"), 3L)),
    drug = factor(rep(c("p", "q"), 6L)),
    site = factor(rep(c("k", "l", "m"), 4L))
  )
  dirichlet <- function(alpha) prior("dirichlet", list(alpha = alpha))
  fit <- .random_update_test_generic_fit(
    ~ 1 + random(1 | id, name = "study", covariance = "diag") +
      random(1 | other, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag") +
      random(1 | site, name = "site", covariance = "diag"),
    data,
    prior_random(
      site = random_block(sd = prior("gamma", list(2, 2))),
      allocation = list(
        random_variance_allocation(
          name = "total_re", terms = c(nested = "nested", drug = "drug"),
          sd = prior("gamma", list(2, 2)), weights = dirichlet(c(1, 1))
        ),
        random_variance_allocation(
          name = "nested_split", parent = allocation_ref("total_re", "nested"),
          terms = c(study = "study", paper = "paper"), weights = dirichlet(c(3, 2))
        )
      )
    ),
    list(intercept = prior("normal", list(0, 1))),
    values = c(
      "mu__xRE_ALLOCx_total_re__allocation_sd" = 0.8,
      "mu__xRE_ALLOCx_total_re__weight[1]" = 0.35,
      "mu__xRE_ALLOCx_total_re__weight[2]" = 0.65,
      "mu__xRE_ALLOCx_nested_split__weight[1]" = 0.3,
      "mu__xRE_ALLOCx_nested_split__weight[2]" = 0.7
    )
  )
  source <- "mu__xRE_ALLOCx_total_re__allocation_sd"
  for(quantity in c("sd_total", "var_total")){
    nested <- .random_update_test_scale_reach(
      fit, paste0("(mu) nested_split: ", quantity), source, factor = 1.6,
      unreached = "site"
    )
    expect_identical(nested$blocks, c("study", "paper", "drug"))
    expect_identical(nested$coefficient_input, "quantity")
    expect_identical(nested$held_fixed, "allocation_weights")
    root <- .random_update_test_scale_reach(
      fit, paste0("(mu) total_re: ", quantity), source, factor = 1.6,
      unreached = "site"
    )
    expect_identical(root$blocks, c("study", "paper", "drug"))
    expect_identical(root$coefficient_input, "source")
  }
})
