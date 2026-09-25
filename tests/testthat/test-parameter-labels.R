skip_if_not_test_profile("unit")

# Synthetic fitted object of a formula: seeded placeholder draws for every
# fitted coordinate, with the parameter map of a real fit.
.label_test_fit <- function(formula, data, prior_list, parameter = "mu",
                            formula_scale = NULL, seed = 1L, n = 20L){

  formula_result <- JAGS_formula(
    formula       = formula,
    parameter     = parameter,
    data          = data,
    prior_list    = prior_list,
    formula_scale = formula_scale
  )
  columns <- unlist(lapply(names(formula_result$prior_list), function(name){
    prior <- formula_result$prior_list[[name]]
    if(.bt_prior_is_factor_family(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  set.seed(seed)
  draws <- matrix(
    stats::rnorm(n * length(columns)),
    nrow = n,
    dimnames = list(NULL, columns)
  )
  scale <- formula_result$formula_scale
  fit <- structure(
    list(
      mcmc         = coda::mcmc.list(coda::mcmc(draws)),
      sample       = n,
      summary.pars = list(mutate = NULL),
      monitor      = columns
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- stats::setNames(
    list(formula_result$formula_design),
    parameter
  )
  if(length(scale) > 0L){
    attr(fit, "formula_scale") <- stats::setNames(list(scale), parameter)
  }
  fit <- .bt_attach_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}

.label_test_factor_prior <- function(contrast){

  switch(
    contrast,
    treatment   = prior_factor("normal", list(0, 1), contrast = "treatment"),
    independent = prior_factor("normal", list(0, 1), contrast = "independent"),
    meandif     = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
    orthonormal = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
    ordered     = prior_ordered(prior("normal", list(0, 1)))
  )
}

.label_test_data <- function(levels, contrast, h_levels = NULL){

  g <- factor(rep(levels, 4L), levels = levels)
  if(identical(contrast, "ordered")){
    g <- ordered(g, levels = levels)
  }
  data <- data.frame(g = g, x = seq(-1, 1, length.out = length(g)))
  if(!is.null(h_levels)){
    data$h <- factor(rep(h_levels, length.out = nrow(data)), levels = h_levels)
  }
  data
}

test_that("the label renderer renders every style from structured parts", {

  cell <- .bt_label_parts(
    components = c("g", "x"),
    formula_parameter = "mu",
    levels = c(g = "b")
  )
  expect_identical(.bt_label(cell, "selector"), "mu_g__xXx__x[b]")
  expect_identical(.bt_label(cell, "table"), "(mu) g[b]:x")
  expect_identical(.bt_label(cell, "table", formula_prefix = FALSE), "g[b]:x")
  expect_identical(.bt_label(cell, "plot"), "b")
  expect_identical(.bt_label(cell, "warning"), "(mu) g[b]:x")

  interaction <- .bt_label_parts(
    components = c("g", "h"),
    formula_parameter = "mu",
    levels = c(g = "a,b", h = "[u]")
  )
  expect_identical(
    .bt_label(interaction, "selector"),
    "mu_g__xXx__h[g=\"a,b\", h=%5Bu%5D]"
  )
  expect_identical(.bt_label(interaction, "table"), "(mu) g[a,b]:h[[u]]")
  expect_identical(.bt_label(interaction, "plot"), "a,b, [u]")

  dif <- .bt_label_parts_update(interaction, transformation = "dif")[[1L]]
  expect_identical(
    .bt_label(dif, "selector"),
    "mu_g[dif: a,b]__xXx__h[dif: [u]]"
  )
  expect_identical(.bt_label(dif, "table"), "(mu) g[dif: a,b]:h[dif: [u]]")
  # legends show the level text only, without the transformation marker
  expect_identical(.bt_label(dif, "plot"), "a,b, [u]")

  coefficient <- .bt_label_parts(
    components = c("g", "x"),
    formula_parameter = "mu",
    coefficient = 2L
  )
  expect_identical(.bt_label(coefficient, "selector"), "mu_g__xXx__x{2}")
  expect_identical(.bt_label(coefficient, "table"), "(mu) g:x{2}")
  expect_identical(.bt_label(coefficient, "plot"), "{2}")

  intercept <- .bt_label_parts("intercept", formula_parameter = "mu")
  expect_identical(.bt_label(intercept, "selector"), "mu_intercept")
  expect_identical(
    .bt_label(.bt_label_parts_update(intercept, transformation = "exp"), "table"),
    "(mu) exp(intercept)"
  )
  inclusion <- .bt_label_parts_update(intercept, inclusion = "")[[1L]]
  expect_identical(.bt_label(inclusion, "table"), "(mu) intercept (inclusion)")
  expect_identical(
    .bt_label(.bt_label_parts_update(intercept, inclusion = "alt"), "table"),
    "(mu) intercept (inclusion: alt)"
  )

  marginal <- .bt_label_parts(
    components = c("x", "g"),
    formula_parameter = "mu",
    levels = c(x = "-1SD", g = "b"),
    marginal = TRUE
  )
  expect_identical(.bt_label(marginal, "table"), "(mu) x:g[-1SD, b]")
  expect_identical(.bt_label(marginal, "plot"), "-1SD, b")
  expect_identical(.bt_label(marginal, "selector"), "mu_x__xXx__g[-1SD, b]")

  random <- .bt_label_parts(
    components = "id",
    formula_parameter = "mu",
    random = list(
      owner = "", quantity = "sd", arguments = "intercept",
      display_arguments = character()
    )
  )
  expect_identical(.bt_label(random, "selector"), "(mu) sd(intercept)")
  expect_identical(.bt_label(random, "table", formula_prefix = FALSE),
                   "sd(intercept)")
  expect_identical(.bt_label(random, "table", simplify = TRUE), "(mu) sd")

  plain <- .bt_label_parts("omega[0,0.05]")
  expect_identical(
    .bt_label(list(plain, plain), "table"),
    rep("omega[0,0.05]", 2L)
  )

  expect_error(
    .bt_label_parts(c("g", "x"), levels = c(z = "a")),
    "label parts are malformed"
  )
  expect_error(
    .bt_label_parts("g", levels = c(g = "a"), coefficient = 1L),
    "label parts are malformed"
  )
})

test_that("catalog level tokens round trip through the one codec", {

  levels <- c(
    "a", "5", "-1", "1.5", "a b", "a,b", "a=b", " a", "b ", "",
    "x]y", "[z]", "w{1}", "}{", "100%", "%25", "back\\slash", "q\"uote",
    "tick`", "tab\there", "line\nbreak", "\u00e9t\u00e9"
  )
  tokens <- .bt_label_token(levels)
  expect_identical(.bt_label_token_decode(tokens), levels)
  expect_false(any(grepl("[][{}`\\\\]", tokens)))
  # tokens with interaction separators are quoted, so cells split uniquely
  expect_true(all(startsWith(tokens[levels %in% c("a,b", "a=b", " a", "b ", "")], "\"")))
  expect_identical(anyDuplicated(tokens), 0L)
})

test_that("every catalog quantity renders its canonical name and resolves by its labels", {

  level_sets <- list(
    c("5", "10", "20"),
    c("1", "2", "3", "4"),
    c("lo", "hi"),
    c("a b", "a,b", "a=b"),
    c("x]y", "[z]", "w{1}")
  )
  checked <- 0L
  for(contrast in c("treatment", "meandif", "orthonormal", "ordered")){
    for(levels in level_sets){
      data <- .label_test_data(levels, contrast, h_levels = c("u", "v"))
      fit <- .label_test_fit(
        ~ g * x + g * h,
        data,
        list(
          intercept = prior("normal", list(0, 1)),
          g         = .label_test_factor_prior(contrast),
          x         = prior("normal", list(0, 1)),
          h         = .label_test_factor_prior("treatment"),
          "g:x"     = .label_test_factor_prior(contrast),
          "g:h"     = .label_test_factor_prior(contrast)
        )
      )
      catalog <- parameter_catalog(fit)
      quantities <- catalog$quantities[!catalog$quantities$internal, , drop = FALSE]
      info <- paste(contrast, paste(levels, collapse = "|"))

      expect_identical(
        parameter_labels(quantities, "selector"),
        quantities$canonical_name,
        info = info
      )
      expect_identical(
        parameter_labels(quantities, "table", simplify = TRUE),
        quantities$display_label,
        info = info
      )
      for(prefix in c(TRUE, FALSE)){
        labels <- parameter_labels(quantities, "table", formula_prefix = prefix)
        for(i in seq_along(labels)){
          resolved <- parameter_catalog_resolve(catalog, labels[[i]])
          expect_identical(
            resolved$quantity_id,
            quantities$quantity_id[[i]],
            info = paste(info, labels[[i]])
          )
          checked <- checked + 1L
        }
      }
      # the level names and rows of transformed contrasts select the levels
      cells <- vapply(quantities$label_parts, function(parts){
        length(parts$levels) > 0L
      }, logical(1))
      dif <- .bt_label_parts_update(
        unclass(quantities$label_parts)[cells],
        transformation = "dif"
      )
      for(style in c("selector", "table")){
        dif_labels <- .bt_label(dif, style)
        for(i in seq_along(dif_labels)){
          expect_identical(
            parameter_catalog_resolve(catalog, dif_labels[[i]])$quantity_id,
            quantities$quantity_id[cells][[i]],
            info = paste(info, dif_labels[[i]])
          )
        }
      }
    }
  }
  expect_gt(checked, 500L)
})

test_that("two-level interaction cells are labelled by their level", {

  data <- .label_test_data(c("lo", "hi"), "treatment")
  fit <- .label_test_fit(
    ~ g * x,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("treatment"),
      x         = prior("normal", list(0, 1)),
      "g:x"     = .label_test_factor_prior("treatment")
    )
  )
  catalog <- parameter_catalog(fit)
  cell <- parameter_catalog_resolve(catalog, "(mu) g[hi]:x")
  expect_identical(cell$quantities$canonical_name, "mu_g__xXx__x[hi]")
  expect_identical(cell$quantities$display_label, "(mu) g[hi]:x")
  expect_identical(
    parameter_coordinates(fit)$display_label[
      parameter_coordinates(fit)$coordinate_name == "mu_g__xXx__x"
    ],
    "(mu) g[hi]:x"
  )
})

test_that("mixed posteriors carry the fitted coordinate and catalog quantity of every column", {

  data <- .label_test_data(c("a", "b", "c"), "meandif", h_levels = c("u", "v"))
  fit <- .label_test_fit(
    ~ g * x + h,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("meandif"),
      x         = prior("normal", list(0, 1)),
      h         = .label_test_factor_prior("treatment"),
      "g:x"     = .label_test_factor_prior("meandif")
    )
  )
  catalog <- parameter_catalog(fit)
  draws <- as.matrix(fit$mcmc)
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  for(parameter in names(mixed)){
    samples <- mixed[[parameter]]
    quantities <- posterior_metadata(samples, "quantities")
    expect_identical(
      quantities$column,
      if(is.null(dim(samples))) parameter else colnames(samples),
      info = parameter
    )
    # every column is exactly its fitted coordinate ...
    expect_true(all(lengths(quantities$dependencies) == 1L), info = parameter)
    expect_equal(
      matrix(as.numeric(samples), nrow = NROW(samples)),
      unname(draws[, unlist(quantities$dependencies), drop = FALSE]),
      info = parameter
    )
    # ... and the catalog quantity that is that coordinate
    rows <- match(quantities$quantity_id, catalog$quantities$quantity_id)
    expect_false(anyNA(rows), info = parameter)
    expect_identical(
      lapply(catalog$quantities$extraction_key[rows], `[[`, "dependencies"),
      unclass(quantities$dependencies),
      info = parameter
    )
  }

  # transformed levels are the design combinations of the fitted coordinates
  transformed <- transform_factor_samples(mixed)
  for(parameter in c("mu_g", "mu_g__xXx__x")){
    quantities <- posterior_metadata(transformed[[parameter]], "quantities")
    expect_identical(quantities$column, colnames(transformed[[parameter]]))
    for(i in seq_len(nrow(quantities))){
      expect_equal(
        as.numeric(transformed[[parameter]][, i]),
        as.numeric(
          draws[, quantities$dependencies[[i]], drop = FALSE] %*%
            quantities$weights[[i]]
        ),
        tolerance = 1e-12,
        info = paste(parameter, i)
      )
    }
  }
})

# Original-scale coefficients of the formula g + g:x (treatment g, level-wise
# slopes, x standardized with mean m and sd s): at level k the standardized
# linear predictor b0 + g_k + c_k (x - m) / s equals
# (b0 + g_k - c_k m / s) + (c_k / s) x, so the intercept is b0 - c_1 m / s,
# level k adds g_k - (c_k - c_1) m / s, and the slopes are c_k / s.
.label_test_level_slope_original <- function(intercept, levels, slopes, m, s){

  list(
    intercept = intercept - slopes[, 1L] * m / s,
    levels    = levels - (slopes[, -1L, drop = FALSE] - slopes[, 1L]) * m / s,
    slopes    = slopes / s
  )
}

test_that("original-scale ensemble tables transform mixed columns through the fitted design", {

  set.seed(11)
  data <- data.frame(
    g = factor(rep(c("1", "2", "3"), 8), levels = c("1", "2", "3")),
    x = stats::rnorm(24, 3, 2)
  )
  priors <- list(
    intercept = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("treatment"),
    "g:x"     = .label_test_factor_prior("independent")
  )
  fit <- .label_test_fit(~ g + g:x, data, priors, formula_scale = list(x = TRUE))
  formula_scale <- attr(fit, "formula_scale")
  m <- formula_scale$mu$mu_x$mean
  s <- formula_scale$mu$mu_x$sd
  draws <- as.matrix(fit$mcmc)
  expected <- .label_test_level_slope_original(
    intercept = draws[, "mu_intercept"],
    levels    = draws[, c("mu_g[1]", "mu_g[2]"), drop = FALSE],
    slopes    = draws[, paste0("mu_g__xXx__x[", 1:3, "]"), drop = FALSE],
    m         = m,
    s         = s
  )

  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  table <- ensemble_estimates_table(
    mixed,
    parameters       = names(mixed),
    transform_scaled = TRUE,
    formula_scale    = formula_scale
  )
  expect_equal(table["(mu) intercept", "Mean"], mean(expected$intercept), tolerance = 1e-10)
  expect_equal(
    table[c("(mu) g[2]", "(mu) g[3]"), "Mean"],
    unname(colMeans(expected$levels)),
    tolerance = 1e-10
  )
  expect_equal(
    table[c("(mu) g[1]:x", "(mu) g[2]:x", "(mu) g[3]:x"), "Mean"],
    unname(colMeans(expected$slopes)),
    tolerance = 1e-10
  )
  # the same numbers as the design-verified model table
  reference <- JAGS_estimates_table(fit, transform_scaled = TRUE, remove_diagnostics = TRUE)
  expect_equal(
    table[rownames(reference), "Mean"],
    reference[, "Mean"],
    tolerance = 1e-10
  )

  # mixtures of models are transformed the same way, draw by draw
  second <- .label_test_fit(~ g + g:x, data, priors, formula_scale = list(x = TRUE),
                            seed = 2L)
  mixed_models <- mix_posteriors(
    model_list = list(
      list(fit = fit,    marglik = bridgesampling_object(-10), prior_weights = 1),
      list(fit = second, marglik = bridgesampling_object(-10.5), prior_weights = 1)
    ),
    parameters   = names(mixed),
    is_null_list = stats::setNames(rep(list(c(FALSE, FALSE)), 3L), names(mixed)),
    seed         = 1,
    n_samples    = 200
  )
  mixed_draws <- function(parameter){
    samples <- mixed_models[[parameter]]
    quantities <- posterior_metadata(samples, "quantities")
    out <- matrix(as.numeric(samples), nrow = NROW(samples))
    colnames(out) <- unlist(quantities$dependencies)
    out
  }
  expected_models <- .label_test_level_slope_original(
    intercept = mixed_draws("mu_intercept")[, 1L],
    levels    = mixed_draws("mu_g")[, c("mu_g[1]", "mu_g[2]"), drop = FALSE],
    slopes    = mixed_draws("mu_g__xXx__x"),
    m         = m,
    s         = s
  )
  table_models <- ensemble_estimates_table(
    mixed_models,
    parameters       = names(mixed_models),
    transform_scaled = TRUE,
    formula_scale    = formula_scale
  )
  expect_equal(
    table_models[c("(mu) intercept", "(mu) g[2]", "(mu) g[3]"), "Mean"],
    c(mean(expected_models$intercept), unname(colMeans(expected_models$levels))),
    tolerance = 1e-10
  )
  expect_equal(
    table_models[c("(mu) g[1]:x", "(mu) g[2]:x", "(mu) g[3]:x"), "Mean"],
    unname(colMeans(expected_models$slopes)),
    tolerance = 1e-10
  )

  # columns that do not identify their fitted coordinates are refused
  unlabelled <- mixed
  unlabelled$mu_g <- .bt_meta_set(unlabelled$mu_g, "quantities", NULL)
  expect_error(
    ensemble_estimates_table(
      unlabelled,
      parameters       = names(unlabelled),
      transform_scaled = TRUE,
      formula_scale    = formula_scale
    ),
    "do not identify their fitted coordinates"
  )
  expect_error(
    .build_unscale_matrix(
      c("mu_intercept", "mu_g[2]", "mu_g[3]"),
      formula_scale$mu,
      prefix = "mu"
    ),
    "not the fitted coefficient coordinates"
  )
})

test_that("original-scale marginal tables keep the marginal predictions", {

  set.seed(3)
  data <- data.frame(
    x = stats::rnorm(30, 5, 2),
    g = factor(rep(c("a", "b"), 15))
  )
  fit <- .label_test_fit(
    ~ x + g,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("treatment")
    ),
    formula_scale = list(x = TRUE),
    n = 50L
  )
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  marginal <- list(mu_x = marginal_posterior(mixed, "mu_x", formula = ~ x + g))
  inference <- list(mu_x = stats::setNames(
    as.list(rep(1, length(marginal$mu_x))),
    names(marginal$mu_x)
  ))
  # the predictions at x = m + k s do not depend on the standardization
  standardized <- marginal_estimates_table(marginal, inference, parameters = "mu_x")
  original <- marginal_estimates_table(
    marginal, inference,
    parameters       = "mu_x",
    transform_scaled = TRUE,
    formula_scale    = attr(fit, "formula_scale")
  )
  expect_equal(original[, "Mean"], standardized[, "Mean"], tolerance = 1e-12)
  expect_equal(
    standardized[, "Mean"],
    unname(vapply(marginal$mu_x, mean, numeric(1))),
    tolerance = 1e-12
  )
})

test_that("mixed-posterior columns and ensemble rows are catalog labels of their quantities", {

  for(contrast in c("treatment", "meandif", "ordered")){
    data <- .label_test_data(c("a", "b", "c"), contrast, h_levels = c("u", "v"))
    fit <- .label_test_fit(
      ~ g * x + g * h,
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        g         = .label_test_factor_prior(contrast),
        x         = prior("normal", list(0, 1)),
        h         = .label_test_factor_prior("treatment"),
        "g:x"     = .label_test_factor_prior(contrast),
        "g:h"     = .label_test_factor_prior(contrast)
      )
    )
    catalog <- parameter_catalog(fit)
    mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
    for(parameter in names(mixed)){
      quantities <- posterior_metadata(mixed[[parameter]], "quantities")
      # every mixed column is the canonical name of its catalog quantity
      for(i in seq_len(nrow(quantities))){
        expect_identical(
          parameter_catalog_resolve(catalog, quantities$column[[i]])$quantity_id,
          quantities$quantity_id[[i]],
          info = paste(contrast, quantities$column[[i]])
        )
      }
    }
    # every ensemble row, untransformed and transformed, selects its quantity
    for(transform_factors in c(FALSE, TRUE)){
      table <- ensemble_estimates_table(
        mixed,
        parameters        = names(mixed),
        transform_factors = transform_factors
      )
      rows <- rownames(table)
      rows <- rows[rows != "(mu) intercept"]
      for(row in rows){
        selection <- parameter_catalog_resolve(catalog, row)
        expect_equal(
          mean(as.matrix(parameter_draws(fit, selection))),
          table[row, "Mean"],
          tolerance = 1e-10,
          info = paste(contrast, transform_factors, row)
        )
      }
    }
  }
})

test_that("formula prefixes come from the formula parameter, not from name prefixes", {

  data <- data.frame(
    x        = seq(-1, 1, length.out = 12),
    z        = stats::rnorm(12),
    mu_income = seq(0, 1, length.out = 12)
  )
  mu <- JAGS_formula(
    ~ x + mu_income, "mu", data,
    list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1)),
      mu_income = prior("normal", list(0, 1))
    )
  )
  mu_tau <- JAGS_formula(
    ~ z, "mu_tau", data,
    list(
      intercept = prior("normal", list(0, 1)),
      z         = prior("normal", list(0, 1))
    )
  )
  prior_list <- c(mu$prior_list, mu_tau$prior_list)
  set.seed(4)
  draws <- matrix(
    stats::rnorm(20 * length(prior_list)),
    nrow = 20,
    dimnames = list(NULL, names(prior_list))
  )
  fit <- structure(
    list(
      mcmc         = coda::mcmc.list(coda::mcmc(draws)),
      sample       = 20L,
      summary.pars = list(mutate = NULL),
      monitor      = names(prior_list)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_design") <- list(
    mu     = mu$formula_design,
    mu_tau = mu_tau$formula_design
  )
  fit <- .bt_attach_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)

  expected <- c(
    "(mu) intercept", "(mu) x", "(mu) mu_income",
    "(mu_tau) intercept", "(mu_tau) z"
  )
  mixed <- as_mixed_posteriors(fit, parameters = names(prior_list))
  expect_identical(
    rownames(ensemble_estimates_table(mixed, parameters = names(prior_list))),
    expected
  )
  expect_identical(
    rownames(ensemble_estimates_table(
      mixed,
      parameters     = c("mu_x", "mu_mu_income", "mu_tau_z"),
      formula_prefix = FALSE
    )),
    c("x", "mu_income", "z")
  )
  catalog <- parameter_catalog(fit)
  expect_identical(
    parameter_catalog_resolve(catalog, "(mu_tau) z")$quantities$canonical_name,
    "mu_tau_z"
  )

  models <- list(list(
    fit       = fit,
    inference = list(m_number = 1, marglik = 0, prior_prob = 1,
                     post_prob = 1, inclusion_BF = 1)
  ))
  expect_identical(
    colnames(ensemble_summary_table(models, names(prior_list)))[2:6],
    expected
  )
})

test_that("model tables, diagnostics, and ensembles label each quantity alike", {

  for(contrast in c("treatment", "meandif")){
    data <- .label_test_data(c("lo", "hi"), contrast, h_levels = c("u", "v"))
    fit <- .label_test_fit(
      ~ g * x + g * h,
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        g         = .label_test_factor_prior(contrast),
        x         = prior("normal", list(0, 1)),
        h         = .label_test_factor_prior("treatment"),
        "g:x"     = .label_test_factor_prior(contrast),
        "g:h"     = .label_test_factor_prior(contrast)
      )
    )
    catalog <- parameter_catalog(fit)
    mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
    for(transform_factors in c(FALSE, TRUE)){
      for(formula_prefix in c(TRUE, FALSE)){
        info <- paste(contrast, transform_factors, formula_prefix)
        table <- JAGS_estimates_table(
          fit,
          transform_factors = transform_factors,
          formula_prefix    = formula_prefix
        )
        rows <- setdiff(rownames(table), c("(mu) intercept", "intercept"))
        # every row selects the quantity whose draws it summarizes
        for(row in rows){
          selection <- parameter_catalog_resolve(catalog, row)
          expect_equal(
            mean(as.matrix(parameter_draws(fit, selection))),
            table[row, "Mean"],
            tolerance = 1e-10,
            info = paste(info, row)
          )
        }
        # one label per quantity on the model and ensemble routes
        ensemble <- ensemble_estimates_table(
          mixed,
          parameters        = names(mixed),
          transform_factors = transform_factors,
          formula_prefix    = formula_prefix
        )
        expect_identical(rownames(ensemble), rownames(table), info = info)
        expect_equal(ensemble[, "Mean"], table[, "Mean"], tolerance = 1e-10,
                     info = info)
      }
    }
    # the two-level interaction cells name their level
    table <- JAGS_estimates_table(fit)
    if(contrast == "treatment"){
      expect_true(all(c("(mu) g[hi]:x", "(mu) g[hi]:h[v]") %in% rownames(table)))
    }
    # diagnostic titles are the same labels
    plot_data <- .diagnostics_plot_data(
      fit               = fit,
      parameter         = "mu_g__xXx__x",
      prior_list        = attr(fit, "prior_list"),
      transformations   = NULL,
      transform_factors = FALSE
    )
    expect_identical(
      .bt_label(attr(plot_data, "label_parts"), style = "table"),
      grep("^\\(mu\\) g.*:x", rownames(table), value = TRUE),
      info = contrast
    )
  }
})

test_that("model tables and the catalog take formula prefixes from the formula parameter", {

  data <- data.frame(
    x = seq(-1, 1, length.out = 12),
    z = stats::rnorm(12)
  )
  mu <- JAGS_formula(
    ~ x, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)))
  )
  mu_tau <- JAGS_formula(
    ~ z, "mu_tau", data,
    list(intercept = prior("normal", list(0, 1)), z = prior("normal", list(0, 1)))
  )
  prior_list <- c(mu$prior_list, mu_tau$prior_list)
  set.seed(5)
  draws <- matrix(
    stats::rnorm(20 * length(prior_list)),
    nrow = 20,
    dimnames = list(NULL, names(prior_list))
  )
  fit <- structure(
    list(
      mcmc         = coda::mcmc.list(coda::mcmc(draws)),
      sample       = 20L,
      summary.pars = list(mutate = NULL),
      monitor      = names(prior_list)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_design") <- list(
    mu     = mu$formula_design,
    mu_tau = mu_tau$formula_design
  )
  fit <- .bt_attach_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)

  expect_identical(
    rownames(JAGS_estimates_table(fit)),
    c("(mu) intercept", "(mu) x", "(mu_tau) intercept", "(mu_tau) z")
  )
  catalog <- parameter_catalog(fit)
  expect_identical(
    parameter_catalog_resolve(catalog, "(mu_tau) z")$quantities$canonical_name,
    "mu_tau_z"
  )
  expect_error(
    parameter_catalog_resolve(catalog, "(mu) tau_z"),
    class = "BayesTools_parameter_not_found"
  )
})

test_that("marginal names, rows, and warnings are rendered labels of the marginal means", {

  data <- .label_test_data(c("a", "b"), "treatment")
  fit <- .label_test_fit(
    ~ x * g,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("treatment"),
      "x:g"     = .label_test_factor_prior("treatment")
    )
  )
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  marginal <- marginal_posterior(mixed, "mu_x__xXx__g", formula = ~ x * g)
  expected_names <- paste0(rep(c("-1SD", "0SD", "1SD"), 2L), ", ",
                           rep(c("a", "b"), each = 3L))
  expect_identical(names(marginal), expected_names)
  for(level in names(marginal)){
    parts <- posterior_metadata(marginal[[level]], "quantities")$label_parts[[1L]]
    expect_true(parts$marginal)
    expect_identical(.bt_label(parts, "plot"), level)
  }

  inference <- list(mu_x__xXx__g = stats::setNames(
    lapply(names(marginal), function(level){
      structure(1, warnings = "check this level")
    }),
    names(marginal)
  ))
  table <- marginal_estimates_table(
    list(mu_x__xXx__g = marginal),
    inference,
    parameters = "mu_x__xXx__g"
  )
  expect_identical(
    rownames(table),
    paste0("(mu) x:g[", expected_names, "]")
  )
  # the warnings name the rows they belong to
  expect_identical(
    attr(table, "warnings"),
    paste0(rownames(table), ": check this level")
  )
  table <- marginal_estimates_table(
    list(mu_x__xXx__g = marginal),
    inference,
    parameters     = "mu_x__xXx__g",
    formula_prefix = FALSE
  )
  expect_identical(
    attr(table, "warnings"),
    paste0("x:g[", expected_names, "]: check this level")
  )
})

# The labels of the discrete (legend) scales of a ggplot.
.label_test_legends <- function(plot){

  built  <- ggplot2::ggplot_build(plot)
  scales <- Filter(function(scale) inherits(scale, "ScaleDiscrete"), built$plot$scales$scales)
  unique(lapply(scales, function(scale) scale$get_labels()))
}

test_that("factor plot legends are the level text of each level cell", {

  legends <- function(samples, parameter){
    plot_data <- .plot_data_samples.factor(
      samples,
      parameter                = parameter,
      n_points                 = 64,
      transformation           = NULL,
      transformation_arguments = NULL,
      transformation_settings  = FALSE
    )
    .plot_prior_factor_normalize_data(plot_data)$level_names
  }

  # treatment x treatment: every cell names the levels of both factors
  data <- .label_test_data(c("a", "b", "c"), "treatment")
  data$g <- factor(rep(c("a", "b", "c"), length.out = nrow(data)))
  data$h <- factor(rep(c("u", "v", "w"), each = 3L, length.out = nrow(data)))
  fit <- .label_test_fit(
    ~ g * h,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("treatment"),
      h         = .label_test_factor_prior("treatment"),
      "g:h"     = .label_test_factor_prior("treatment")
    ),
    n = 200L
  )
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  cells <- c("b, v", "c, v", "b, w", "c, w")
  expect_identical(legends(mixed, "mu_g__xXx__h"), cells)
  expect_identical(
    .label_test_legends(plot_posterior(mixed, "mu_g__xXx__h", plot_type = "ggplot")),
    list(cells)
  )
  expect_identical(
    .label_test_legends(plot_posterior(mixed, "mu_g__xXx__h", plot_type = "ggplot", prior = TRUE)),
    list(cells)
  )

  # mean-difference levels: the level text without the transformation
  # marker or padding, with brackets and braces kept whole
  data <- .label_test_data(c("x]y", "[z]", "w{1}"), "meandif", h_levels = c("u", "v"))
  fit <- .label_test_fit(
    ~ g * h,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("meandif"),
      h         = .label_test_factor_prior("meandif"),
      "g:h"     = .label_test_factor_prior("meandif")
    ),
    n = 200L
  )
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  transformed <- transform_factor_samples(mixed)
  levels <- c("x]y", "[z]", "w{1}")
  cells  <- c("x]y, u", "[z], u", "w{1}, u", "x]y, v", "[z], v", "w{1}, v")
  expect_identical(legends(mixed, "mu_g"), levels)
  expect_identical(legends(transformed, "mu_g"), levels)
  expect_identical(legends(transformed, "mu_g__xXx__h"), cells)
  expect_identical(
    .label_test_legends(plot_posterior(
      mixed, "mu_g__xXx__h", plot_type = "ggplot", transform_factors = TRUE
    )),
    list(cells)
  )
  expect_identical(
    .label_test_legends(plot_posterior(
      mixed, "mu_g", plot_type = "ggplot", transform_factors = TRUE, prior = TRUE
    )),
    list(levels)
  )
})

test_that("ordered factors plot their level priors with the posterior", {

  data <- .label_test_data(c("lo", "mid", "hi"), "ordered")
  fit <- .label_test_fit(
    ~ 1 + g,
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("ordered")
    ),
    n = 200L
  )
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  # the first level is fixed at zero by the ordered contrast; the priors are
  # drawn for the plotted levels
  expect_identical(
    .label_test_legends(plot_posterior(mixed, "mu_g", plot_type = "ggplot", prior = TRUE)),
    list(c("mid", "hi"))
  )
  expect_identical(
    .label_test_legends(plot_posterior(
      mixed, "mu_g", plot_type = "ggplot", prior = TRUE, transform_factors = TRUE
    )),
    list(c("mid", "hi"))
  )
})
