test_that("variance-allocation gate priors sample independent indicators", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug  = factor(c("a", "b", "a", "b"))
  )
  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name      = "total",
        terms     = c(study = "study", drug = "drug"),
        sd        = prior("gamma", list(2, 2)),
        weights   = prior("dirichlet", list(alpha = c(1, 1))),
        inclusion = list(
          study = prior("spike", list(location = 0.3)),
          drug  = prior("spike", list(location = 0.7))
        )
      )
    )
  )

  allocation <- result$formula_design$random_effects[[1L]]$
    sd_binding$allocations[[1L]]
  indicators <- vapply(
    allocation$inclusion,
    `[[`,
    character(1),
    "indicator_name"
  )
  posterior_columns <- c(
    "mu_intercept",
    allocation$source_node,
    paste0(allocation$weight_name, "[", 1:2, "]"),
    vapply(allocation$inclusion, `[[`, character(1), "prob_name"),
    indicators
  )
  fit <- coda::mcmc.list(coda::mcmc(matrix(
    0,
    nrow = 2L,
    ncol = length(posterior_columns),
    dimnames = list(NULL, posterior_columns)
  )))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list")     <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  samples_a <- transform_prior_samples(fit, n_samples = 10000L, seed = 913L)
  samples_b <- transform_prior_samples(fit, n_samples = 10000L, seed = 913L)
  expect_identical(samples_a, samples_b)
  expect_true(all(indicators %in% colnames(samples_a)))
  expect_true(all(samples_a[, indicators] %in% c(0, 1)))
  expect_equal(unname(colMeans(samples_a[, indicators])), c(0.3, 0.7),
               tolerance = 0.015)
  expect_lt(abs(stats::cor(samples_a[, indicators])[1L, 2L]), 0.03)

  fit <- attach_test_parameter_map(fit)
  selection <- parameter_catalog_resolve(
    parameter_catalog(fit),
    "(mu) total: inclusion(study)"
  )
  draws <- parameter_draws(
    fit,
    selection,
    model_samples = samples_a
  )
  expect_equal(as.numeric(as.matrix(draws)), samples_a[, indicators[[1L]]])

  weight_columns <- paste0(allocation$weight_name, "[", 1:2, "]")
  gated_weights <- samples_a[, weight_columns, drop = FALSE] *
    samples_a[, indicators, drop = FALSE]
  realized_fraction <- rowSums(gated_weights)
  source <- samples_a[, allocation$source_node]

  total_selection <- parameter_catalog_resolve(
    parameter_catalog(fit),
    "(mu) total: sd_total"
  )
  total_draws <- as.numeric(as.matrix(parameter_draws(
    fit,
    total_selection,
    model_samples = samples_a
  )))
  expect_equal(
    total_draws,
    source * sqrt(realized_fraction),
    tolerance = 1e-12
  )
  expect_true(all(total_draws[realized_fraction == 0] == 0))
  expect_true(all(c(weight_columns, indicators) %in%
                    total_selection$quantities$extraction_key[[1L]]$dependencies))
  expect_identical(total_selection$quantities$source_type, "composite")
  total_plan <- random_effects_marginal_update_plan(fit, total_selection)
  expect_identical(total_plan$family, "unsupported")
  expect_identical(total_plan$reason, "gated_realized_allocation")

  proportion_draws <- lapply(c("study", "drug"), function(component){
    proportion_selection <- parameter_catalog_resolve(
      parameter_catalog(fit),
      paste0("(mu) total: var_prop(", component, ")")
    )
    expect_true(all(c(weight_columns, indicators) %in%
                      proportion_selection$quantities$extraction_key[[1L]]$dependencies))
    expect_identical(proportion_selection$quantities$source_type, "composite")
    proportion_plan <- random_effects_marginal_update_plan(
      fit,
      proportion_selection
    )
    expect_identical(proportion_plan$family, "unsupported")
    as.numeric(as.matrix(parameter_draws(
      fit,
      proportion_selection,
      model_samples = samples_a
    )))
  })
  active <- realized_fraction > 0
  expect_true(all(vapply(proportion_draws, function(x){
    all(is.na(x[!active]))
  }, logical(1))))
  expect_equal(
    unname(do.call(cbind, proportion_draws)[active, , drop = FALSE]),
    unname(gated_weights[active, , drop = FALSE] /
             realized_fraction[active]),
    tolerance = 1e-12
  )
  expect_equal(
    rowSums(do.call(cbind, proportion_draws)[active, , drop = FALSE]),
    rep(1, sum(active)),
    tolerance = 1e-12
  )
})


test_that("a consumed unary root gate retains canonical allocation size", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug  = factor(c("a", "b", "a", "b"))
  )
  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = list(
        random_variance_allocation(
          name      = "root",
          terms     = c(component = "component_gated"),
          sd        = prior("gamma", list(2, 2)),
          inclusion = list(
            component = prior("spike", list(location = 0.5))
          )
        ),
        random_variance_allocation(
          name    = "split",
          terms   = c(study = "study", drug = "drug"),
          parent  = allocation_ref("root", "component"),
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    )
  )

  root <- result$formula_design$random_allocations$root
  split <- result$formula_design$random_allocations$split
  expect_true(root$gate_only)
  expect_identical(root$n_targets, 1L)
  expect_identical(root$components$component$factors[[1L]]$n_targets, 1L)

  source_name <- root$source_node
  gate_name   <- root$inclusion$component$indicator_name
  weight_name <- split$weight_name
  samples <- matrix(
    c(
      0, 1, 0, 0.2, 0.8,
      0, 2, 1, 0.7, 0.3
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      source_name,
      gate_name,
      paste0(weight_name, "[1]"),
      paste0(weight_name, "[2]")
    ))
  )
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

  sd_total <- random_effects_summary_posterior(
    fit,
    summary = "sd_total"
  )
  expect_equal(as.numeric(sd_total[[1L]]), c(0, 2))
})


test_that("a sole hidden allocation omits its public component argument", {

  data <- data.frame(study = factor(c("s1", "s1", "s2", "s2")))
  result <- JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name            = "component",
        display_name    = "",
        terms           = c(component_1 = "study"),
        component_names = "",
        sd              = prior("gamma", list(2, 2)),
        inclusion       = list(
          component_1 = prior("spike", list(location = 0.5))
        )
      )
    )
  )
  allocation  <- result$formula_design$random_allocations$component
  source_name <- allocation$source_node
  gate_name   <- allocation$inclusion$component_1$indicator_name
  samples <- matrix(
    c(0, 1, 0, 0, 2, 1),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", source_name, gate_name))
  )
  fit <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(samples)),
      sample = nrow(samples)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  fit <- attach_test_parameter_map(fit)

  selection <- parameter_catalog_resolve(
    parameter_catalog(fit),
    "(mu) inclusion"
  )
  expect_identical(selection$quantities$arguments[[1L]], "")
  expect_false(any(grepl(
    "component 1",
    parameter_catalog(fit)$quantities$canonical_name,
    fixed = TRUE
  )))
})


test_that("allocation slab auxiliaries remain private coordinates", {

  data <- data.frame(study = factor(c("s1", "s1", "s2", "s2")))
  slab <- prior_mixture(
    prior_list = list(
      prior("gamma", list(2, 2)),
      prior("gamma", list(3, 2))
    ),
    components = c("alternative", "alternative")
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
        sd = slab,
        inclusion = list(study = prior("spike", list(location = 0.5)))
      )
    )
  )
  allocation <- result$formula_design$random_allocations$component
  source_name <- allocation$source_node
  slab_indicator <- paste0(source_name, "_indicator")
  gate_indicator <- allocation$inclusion$study$indicator_name
  samples <- matrix(
    c(
      0, 1, 1, 0,
      0, 2, 2, 1
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", source_name, slab_indicator, gate_indicator)
    )
  )
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

  coordinates <- parameter_coordinates(fit)
  product_indicators <- c(slab_indicator, gate_indicator)
  indicator_rows <- match(product_indicators, coordinates$coordinate_name)
  source_row <- match(source_name, coordinates$coordinate_name)
  expect_false(anyNA(c(source_row, indicator_rows)))
  expect_true(coordinates$internal[[source_row]])
  expect_false(any(coordinates$internal[indicator_rows]))
  expect_true(all(coordinates$role[indicator_rows] == "allocation"))

  catalog <- parameter_catalog(fit)
  expect_false(any(c(source_name, product_indicators) %in%
                     catalog$quantities$canonical_name))
  expect_no_error(parameter_catalog_resolve(
    catalog,
    "(mu) component: inclusion(study)"
  ))
  expect_no_error(parameter_catalog_resolve(catalog, "(mu) sd(intercept)"))
  expect_true(all(product_indicators %in% colnames(as.matrix(
    JAGS_materialize_draws(fit)
  ))))

  inference <- JAGS_inference_table(fit)
  expect_identical(
    attr(inference, "parameter_roles"),
    c("random_slab", "random_inclusion")
  )

  summary_samples <- JAGS_estimates_table(
    fit,
    keep_parameters = "random",
    random_effects_summary = "standard",
    remove_inclusion = TRUE,
    return_samples = TRUE
  )
  expect_identical(colnames(summary_samples), "(mu) sd(intercept)")
  expect_equal(summary_samples[, 1L], c(0, 2))
})
