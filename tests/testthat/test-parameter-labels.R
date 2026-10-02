skip_if_not_test_profile("unit")

# Synthetic fitted object of a formula: seeded placeholder draws for every
# fitted coordinate, with the parameter map of a real fit.
.label_test_fit <- function(formula, data, prior_list, parameter = "mu",
                            formula_scale = NULL, seed = 1L, n = 20L,
                            prior_random = NULL){

  formula_result <- JAGS_formula(
    formula       = formula,
    parameter     = parameter,
    data          = data,
    prior_list    = prior_list,
    formula_scale = formula_scale,
    prior_random  = prior_random
  )
  columns <- unlist(lapply(names(formula_result$prior_list), function(name){
    prior <- formula_result$prior_list[[name]]
    if(is.prior.spike_and_slab(prior)){
      # the spike-and-slab nodes: indicator, inclusion, value, and variable
      variable <- .get_spike_and_slab_variable(prior)
      names_of <- function(node){
        if(.bt_prior_is_factor_family(variable)) .JAGS_prior_factor_names(node, variable) else node
      }
      return(c(
        paste0(name, c("_indicator", "_inclusion")),
        names_of(name),
        names_of(paste0(name, "_variable"))
      ))
    }
    if(.bt_prior_is_factor_family(prior)){
      .JAGS_prior_factor_names(name, prior)
    }else if(is.prior.mixture(prior)){
      c(name, paste0(name, "_indicator"))
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
  # indicators: 0/1 for spike-and-slab priors, component indices of mixtures
  for(column in columns[endsWith(columns, "_indicator")]){
    owner <- formula_result$prior_list[[sub("_indicator$", "", column)]]
    draws[, column] <- if(is.prior.spike_and_slab(owner)){
      stats::rbinom(n, 1L, 0.5)
    }else{
      sample.int(length(owner), n, replace = TRUE)
    }
  }
  draws[, endsWith(columns, "_inclusion")] <- 0.5
  # Ordered label fixtures declare independent total draws and equal primitive
  # Gamma coordinates. Their increments are the total divided equally over
  # each ordered axis, so source data and coefficient data describe one model.
  for(name in names(formula_result$prior_list)){
    prior <- formula_result$prior_list[[name]]
    if(!is.prior.ordered(prior)) next
    spec <- .bt_ordered_spec(name, prior)
    totals <- matrix(stats::rnorm(n * length(spec$total_names)), n,
      dimnames = list(NULL, spec$total_names))
    draws <- cbind(draws, totals)
    for(record in spec$allocations){
      if(!identical(record$spec$type, "dirichlet")) next
      draws <- cbind(draws, matrix(1, n, record$dim,
        dimnames = list(NULL, record$gamma_coordinates)))
    }
    dimensions <- vapply(spec$metadata$ordered_terms, function(factor){
      .prior_ordered_allocation_for_coefficient(spec$metadata, factor, 1L)$dim
    }, integer(1))
    draws[, spec$coefficient_names] <- totals[, spec$slice_index, drop = FALSE] / prod(dimensions)
  }
  columns <- colnames(draws)
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
      # every alias is the table label of its label parts
      expect_identical(
        parameter_labels(catalog$aliases, "table"),
        catalog$aliases$alias,
        info = info
      )
      # the quantity every label resolves to, compared for all labels at once
      # (the names of the compared vectors are the labels)
      resolved_ids <- function(labels){
        stats::setNames(vapply(labels, function(label){
          parameter_catalog_resolve(catalog, label)$quantity_id
        }, character(1), USE.NAMES = FALSE), labels)
      }
      for(prefix in c(TRUE, FALSE)){
        labels <- parameter_labels(quantities, "table", formula_prefix = prefix)
        expect_identical(
          resolved_ids(labels),
          stats::setNames(quantities$quantity_id, labels),
          info = info
        )
        checked <- checked + length(labels)
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
        expect_identical(
          resolved_ids(dif_labels),
          stats::setNames(quantities$quantity_id[cells], dif_labels),
          info = paste(info, style)
        )
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

test_that("original-scale mixtures are transformed only through a design containing their terms", {

  # models whose formulas differ: g * x and g + x (no g:x); the g:x columns of
  # the mixture are coefficients of the larger design only
  set.seed(12)
  data <- data.frame(
    g = factor(rep(c("a", "b", "c"), 8), levels = c("a", "b", "c")),
    x = stats::rnorm(24, 3, 2)
  )
  full <- .label_test_fit(~ g * x, data, list(
    intercept = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("treatment"),
    x         = prior("normal", list(0, 1)),
    "g:x"     = .label_test_factor_prior("treatment")
  ), formula_scale = list(x = TRUE), n = 50L)
  smaller <- .label_test_fit(~ g + x, data, list(
    intercept = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("treatment"),
    x         = prior("normal", list(0, 1))
  ), formula_scale = list(x = TRUE), seed = 2L, n = 50L)
  parameters <- c("mu_intercept", "mu_g", "mu_x", "mu_g__xXx__x")
  mixed <- mix_posteriors(
    model_list = list(
      list(fit = full,    marglik = bridgesampling_object(-10),   prior_weights = 1),
      list(fit = smaller, marglik = bridgesampling_object(-10.2), prior_weights = 1)
    ),
    parameters   = parameters,
    is_null_list = list(
      mu_intercept = c(FALSE, FALSE), mu_g = c(FALSE, FALSE),
      mu_x = c(FALSE, FALSE), mu_g__xXx__x = c(FALSE, TRUE)
    ),
    seed      = 1,
    n_samples = 200
  )

  # through the larger design: at level k the standardized predictor
  # b0 + g_k + (bx + c_k) (x - m) / s gives the intercept b0 - bx m / s, the
  # levels g_k - c_k m / s, and the slopes bx / s and c_k / s
  formula_scale <- attr(full, "formula_scale")
  m <- formula_scale$mu$mu_x$mean
  s <- formula_scale$mu$mu_x$sd
  b0 <- as.numeric(mixed$mu_intercept)
  bx <- as.numeric(mixed$mu_x)
  g  <- matrix(as.numeric(mixed$mu_g), ncol = 2L)
  gx <- matrix(as.numeric(mixed$mu_g__xXx__x), ncol = 2L)
  table <- ensemble_estimates_table(
    mixed, parameters = parameters, transform_scaled = TRUE,
    formula_scale = formula_scale
  )
  expect_equal(
    table[, "Mean"],
    c(mean(b0 - bx * m / s), colMeans(g - gx * m / s), mean(bx / s), colMeans(gx / s)),
    tolerance = 1e-10
  )

  # the smaller model's design does not contain g:x: no original-scale values
  expect_error(
    ensemble_estimates_table(
      mixed, parameters = parameters, transform_scaled = TRUE,
      formula_scale = attr(smaller, "formula_scale")
    ),
    paste0(
      "Cannot transform 'mu_g__xXx__x[1]', 'mu_g__xXx__x[2]' to the original ",
      "predictor scale: the fitted design in 'formula_scale' of formula ",
      "parameter 'mu' does not contain these coefficients."
    ),
    fixed = TRUE,
    class = "BayesTools_formula_transform_unavailable"
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
    # Effect and total rows select catalog quantities. Sampled share labels
    # agree with the source rows of the corresponding single-model table.
    for(transform_factors in c(FALSE, TRUE)){
      table <- ensemble_estimates_table(
        mixed,
        parameters        = names(mixed),
        transform_factors = transform_factors
      )
      rows <- rownames(table)
      rows <- rows[rows != "(mu) intercept"]
      model_table <- runjags_estimates_table(fit, transform_factors = transform_factors)
      for(row in rows){
        if(contrast == "ordered" && !transform_factors && grepl("allocation", row, fixed = TRUE)){
          expect_true(row %in% rownames(model_table), info = row)
          expect_equal(table[row, "Mean"], model_table[row, "Mean"], tolerance = 1e-10, info = row)
          next
        }
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

test_that("estimates tables carry the quantity and label parts of every row", {

  # a mean-difference factor: its transformed levels are the [dif: ] rows
  data <- .label_test_data(c("lo", "mid", "hi"), "meandif")
  fit <- .label_test_fit(~ g, data, list(
    intercept = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("meandif")
  ))
  catalog <- parameter_catalog(fit)
  table <- JAGS_estimates_table(fit, transform_factors = TRUE)
  quantities <- attr(table, "quantities")
  expect_identical(names(quantities), c("row", "quantity_id", "label_parts"))
  expect_identical(quantities$row, rownames(table))
  expect_identical(parameter_labels(quantities, style = "table"), rownames(table))
  expect_identical(
    rownames(table)[startsWith(rownames(table), "(mu) g[")],
    c("(mu) g[dif: lo]", "(mu) g[dif: mid]", "(mu) g[dif: hi]")
  )
  row <- match("(mu) g[dif: mid]", rownames(table))
  expect_identical(quantities$label_parts[[row]]$transformation, "dif")
  expect_identical(quantities$label_parts[[row]]$levels, c(g = "mid"))
  # the row is the level quantity, and renders without the formula prefix
  expect_identical(
    quantities$quantity_id[[row]],
    parameter_catalog_resolve(catalog, "mu_g[mid]")$quantity_id
  )
  expect_identical(
    parameter_labels(quantities[row, ], style = "table", formula_prefix = FALSE),
    "g[dif: mid]"
  )
  expect_identical(
    quantities$quantity_id[[match("(mu) intercept", rownames(table))]],
    parameter_catalog_resolve(catalog, "mu_intercept")$quantity_id
  )
  # subsetting keeps the quantities of the kept rows
  subset <- table[c(row, 1L), ]
  expect_identical(attr(subset, "quantities")$row, rownames(subset))
  expect_identical(attr(update(table, remove_parameters = "(mu) g[dif: lo]"), "quantities")$row,
                   setdiff(rownames(table), "(mu) g[dif: lo]"))
  # the ensemble route carries the parts of its rows, and the transformed
  # levels name the catalog level quantities that hold their values, as the
  # model table does
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
  ensemble <- ensemble_estimates_table(mixed, parameters = names(mixed),
                                       transform_factors = TRUE)
  expect_identical(
    parameter_labels(attr(ensemble, "quantities"), style = "table"),
    rownames(ensemble)
  )
  expect_identical(attr(ensemble, "quantities")$row, rownames(table))
  expect_identical(attr(ensemble, "quantities")$quantity_id, quantities$quantity_id)
  expect_true(all(nzchar(attr(ensemble, "quantities")$quantity_id)))
  mid <- parameter_catalog_resolve(catalog, "mu_g[mid]")
  expect_equal(ensemble["(mu) g[dif: mid]", "Mean"],
               mean(as.matrix(parameter_draws(fit, mid))), tolerance = 1e-12)

  # rows whose values 'transformations' changed hold no catalog quantity
  # (as the draws of posterior_transform()): their label parts record the
  # transformation, their labels are unchanged, and their quantity id is ""
  expect_message(
    transformed <- JAGS_estimates_table(fit, transform_factors = TRUE, transformations = list(
      mu_intercept = list(fun = exp),
      mu_g         = list(fun = function(x) 2 * x)
    )),
    "The transformation was applied to the differences from the mean"
  )
  transformed_quantities <- attr(transformed, "quantities")
  expect_identical(rownames(transformed), rownames(table))
  expect_identical(parameter_labels(transformed_quantities, style = "table"), rownames(transformed))
  intercept_row <- match("(mu) intercept", rownames(transformed))
  expect_equal(
    transformed[intercept_row, "Mean"],
    mean(exp(as.matrix(parameter_draws(fit, parameter_catalog_resolve(catalog, "mu_intercept"))))),
    tolerance = 1e-12
  )
  expect_equal(transformed[row, "Mean"], 2 * table[row, "Mean"], tolerance = 1e-12)
  expect_identical(transformed_quantities$quantity_id[c(intercept_row, row)], c("", ""))
  expect_identical(transformed_quantities$label_parts[[intercept_row]]$transformation, c("none", "custom"))
  expect_identical(transformed_quantities$label_parts[[row]]$transformation, c("dif", "custom"))
  # coefficients transformed before the contrast transformation (a treatment
  # factor: its coordinates are level quantities)
  treatment_data <- .label_test_data(c("lo", "mid", "hi"), "treatment")
  treatment_fit <- .label_test_fit(~ g, treatment_data, list(
    intercept = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("treatment")
  ))
  treatment_table <- JAGS_estimates_table(treatment_fit)
  treatment_transformed <- JAGS_estimates_table(
    treatment_fit, transformations = list(mu_g = list(fun = exp))
  )
  level_rows <- startsWith(rownames(treatment_table), "(mu) g[")
  expect_true(any(level_rows))
  expect_true(all(nzchar(attr(treatment_table, "quantities")$quantity_id[level_rows])))
  expect_identical(
    attr(treatment_transformed, "quantities")$quantity_id,
    ifelse(level_rows, "", attr(treatment_table, "quantities")$quantity_id)
  )
  expect_identical(rownames(treatment_transformed), rownames(treatment_table))
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

  # the intercept's marginal mean (its draws are a one-row matrix of
  # predictions) is one quantity, labelled by its term
  intercept <- marginal_posterior(mixed, "mu_intercept", formula = ~ x * g)
  expect_identical(
    parameter_labels(posterior_metadata(intercept$intercept, "quantities"), "table"),
    "(mu) intercept"
  )
  table <- marginal_estimates_table(
    list(mu_intercept = intercept),
    list(mu_intercept = list(intercept = structure(1, warnings = "check the intercept"))),
    parameters = "mu_intercept"
  )
  expect_identical(rownames(table), "(mu) intercept")
  expect_equal(table[, "Mean"], mean(intercept$intercept), tolerance = 1e-12)
  expect_identical(attr(table, "warnings"), "(mu) intercept: check the intercept")
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

test_that("hypotheses reference levels in the catalog's escaped level form", {

  # every escaped character of the level-token codec decodes back
  levels <- c("w{1}", "(0,1]", "a%b", "[z]", "q\"r", " edge ")
  expect_identical(.bt_label_token_decode(.bt_label_token(levels)), levels)
  expect_identical(
    .hypothesis_match_level_names(.bt_label_token(levels), levels, "mu"),
    levels
  )

  context <- .prior_density_context(
    prior_list   = list(
      a = prior("normal", list(mean = 0, sd = 1)),
      b = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("a", "b"),
    n_grid       = 128
  )
  make_posterior <- function(level_names){
    set.seed(2)
    posterior <- list(
      .bt_meta_update(
        structure(stats::rnorm(4000, 0.5, 0.2), class = c("marginal_posterior.simple", "numeric")),
        linear_weights = c(a = 1, b = 0),
        atoms = posterior_atom_attribute()
      ),
      .bt_meta_update(
        structure(stats::rnorm(4000, 0, 0.2), class = c("marginal_posterior.simple", "numeric")),
        linear_weights = c(a = 0, b = 1),
        atoms = posterior_atom_attribute()
      )
    )
    names(posterior) <- level_names
    class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
    attr(posterior, "parameter") <- "mu"
    .bt_meta_set(posterior, "prior_context", context)
  }

  reference <- hypothesis_BF(
    make_posterior(c("A", "B")),
    hypothesis = "mu[A] > mu[B]",
    seed       = 1,
    columns    = "all"
  )
  posterior <- make_posterior(c("w{1}", "v}"))
  for(hypothesis in c(
    "`mu[w{1}]` > `mu[v}]`",
    "mu[w%7B1%7D] > mu[v%7D]",
    "mu[\"w%7B1%7D\"] > `mu[v}]`"
  )){
    out <- hypothesis_BF(posterior, hypothesis = hypothesis, seed = 1,
                         columns = "all")
    expect_equal(attr(out, "raw_BF"), attr(reference, "raw_BF"),
                 tolerance = 1e-12, info = hypothesis)
  }
})

# A synthetic two-chain fit of the columns of 'chains' under 'prior_list'.
.label_test_mock_fit <- function(chains, prior_list){

  fit <- list(
    mcmc         = do.call(coda::mcmc.list, lapply(chains, coda::mcmc)),
    sample       = nrow(chains[[1L]]),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  fit <- .bt_attach_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}

test_that("estimates tables blank the diagnostics of declared structural weight bins", {

  set.seed(22)
  n <- 80
  weight_chain <- function(omega_2){
    cbind(mu = stats::rnorm(n), "omega[1]" = 1, "omega[2]" = omega_2)
  }
  diagnostics <- c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")

  # estimated weights: the reference bin is the only structural constant
  estimated <- .label_test_mock_fit(
    list(weight_chain(stats::rbeta(n, 4, 2)), weight_chain(stats::rbeta(n, 4, 2))),
    list(
      mu    = prior("normal", list(0, 1)),
      omega = prior_weightfunction("one-sided", .05, wf_cumulative(c(2, 4)))
    )
  )
  table <- JAGS_estimates_table(estimated)
  expect_identical(rownames(table), c("mu", "omega[0,0.05]", "omega[0.05,1]"))
  expect_true(all(is.na(unlist(table["omega[0,0.05]", diagnostics]))))
  expect_false(anyNA(unlist(table[c("mu", "omega[0.05,1]"), c("MCMC_error", "ESS", "R_hat")])))

  # fixed weights: every bin is a declared constant
  fixed <- .label_test_mock_fit(
    list(weight_chain(.5), weight_chain(.5)),
    list(
      mu    = prior("normal", list(0, 1)),
      omega = prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
    )
  )
  table <- JAGS_estimates_table(fixed)
  expect_equal(table[["Mean"]], c(mean(c(fixed$mcmc[[1]][, "mu"], fixed$mcmc[[2]][, "mu"])), 1, .5))
  expect_true(all(is.na(unlist(table[c("omega[0,0.05]", "omega[0.05,1]"), diagnostics]))))
  expect_false(anyNA(unlist(table["mu", diagnostics])))

  # a publication-bias mixture: the bins every branch fixes (the reference
  # bin and its two-sided mirror) are the structural constants
  mixture_chain <- function(){
    cbind(
      mu             = stats::rnorm(n),
      bias_indicator = rep(1:3, length.out = n),
      "omega[1]"     = 1,
      "omega[2]"     = stats::rbeta(n, 4, 2),
      "omega[3]"     = stats::rbeta(n, 4, 2),
      "omega[4]"     = stats::rbeta(n, 4, 2),
      "omega[5]"     = 1
    )
  }
  mixture <- .label_test_mock_fit(
    list(mixture_chain(), mixture_chain()),
    list(
      mu   = prior("normal", list(0, 1)),
      bias = prior_mixture(list(
        prior_none(prior_weights = 1),
        prior_weightfunction("two-sided", .05, wf_cumulative(c(1, 1)), prior_weights = 1),
        prior_weightfunction("two-sided", c(.05, .10), wf_cumulative(c(1, 1, 1)), prior_weights = 1)
      ), is_null = c(TRUE, FALSE, FALSE))
    )
  )
  table <- JAGS_estimates_table(mixture)
  bins <- grep("^omega", rownames(table), value = TRUE)
  expect_identical(
    bins,
    c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.95]",
      "omega[0.95,0.975]", "omega[0.975,1]")
  )
  expect_true(all(is.na(unlist(table[bins[c(1, 5)], diagnostics]))))
  expect_false(anyNA(unlist(table[bins[2:4], c("MCMC_error", "ESS", "R_hat")])))
})

test_that("inclusion rows of names containing 'inclusion' are formatted from the marker only", {

  expect_identical(
    format_parameter_names(
      c("mu__xREx__id_xinclusion", "mu__xREx__id_x(inclusion)"),
      formula_random = "id"
    ),
    c("sd(mu_xinclusion|id)", "mu_x|id (inclusion)")
  )
})

test_that("the formula-coordinate encoding is internal to the name map", {

  exports <- getNamespaceExports("BayesTools")
  expect_false(any(c(
    "JAGS_parameter_encode", "JAGS_parameter_decode",
    "JAGS_parameter_encoding_schema"
  ) %in% exports))
  expect_true("JAGS_formula_name_map" %in% exports)
})

test_that("formulas generating the same JAGS node name are reported by their terms", {

  syntax <- "model{\nfor(i in 1:N){\n  y[i] ~ dnorm(a[i] + a_b[i], 4)\n}\n}"
  data <- data.frame(b_x = seq(-1, 1, length.out = 10), x = seq(1, -1, length.out = 10))
  fit <- function(prior_list = NULL){
    JAGS_fit(
      model_syntax       = syntax,
      data               = list(y = rep(0, 10), N = 10),
      prior_list         = prior_list,
      formula_list       = list(a = ~ 0 + b_x, a_b = ~ 0 + x),
      formula_data_list  = list(a = data, a_b = data),
      formula_prior_list = list(
        a   = list(b_x = prior("normal", list(0, 1))),
        a_b = list(x = prior("normal", list(0, 1)))
      ),
      chains = 1, adapt = 100, burnin = 100, sample = 100, silent = TRUE
    )
  }
  expect_error(
    fit(),
    paste0(
      "formula parameter 'a' (term 'b_x') and formula parameter 'a_b' ",
      "(term 'x') define the same JAGS node 'a_b_x'. Rename a formula ",
      "parameter or predictor so that the node names differ."
    ),
    fixed = TRUE
  )

  # a formula node that is also a model parameter
  expect_error(
    JAGS_fit(
      model_syntax       = "model{\nfor(i in 1:N){\n  y[i] ~ dnorm(a[i], a_x)\n}\n}",
      data               = list(y = rep(0, 10), N = 10),
      prior_list         = list(a_x = prior("gamma", list(1, 1))),
      formula_list       = list(a = ~ 0 + x),
      formula_data_list  = list(a = data),
      formula_prior_list = list(a = list(x = prior("normal", list(0, 1)))),
      chains = 1, adapt = 100, burnin = 100, sample = 100, silent = TRUE
    ),
    "formula parameter 'a' (term 'x') and 'prior_list' define the same JAGS node 'a_x'.",
    fixed = TRUE
  )
})

test_that("mixed columns of spike-and-slab and random-effect factor priors are named by their terms", {

  # the variable part of a spike-and-slab formula term keeps the formula
  # parameter of the term
  data <- .label_test_data(c("a", "b", "c"), "meandif")
  fit <- .label_test_fit(~ x * g, data, list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("meandif"),
    "x:g"     = prior_spike_and_slab(
      .label_test_factor_prior("meandif"),
      prior_inclusion = prior("beta", list(1, 1))
    )
  ))
  mixed <- as_mixed_posteriors(fit, parameters = "mu_x__xXx__g")
  expect_identical(colnames(mixed$mu_x__xXx__g), c("mu_x__xXx__g{1}", "mu_x__xXx__g{2}"))
  expect_identical(
    parameter_labels(posterior_metadata(mixed$mu_x__xXx__g, "quantities"), "table"),
    c("(mu) x:g{1}", "(mu) x:g{2}")
  )

  # the SD coordinates of a random treatment slope are its level cells
  data <- .label_test_data(c("a", "b", "c"), "treatment")
  data$id <- factor(rep(c("s1", "s2", "s3", "s4"), each = 3L))
  fit <- .label_test_fit(
    ~ 1 + g + (1 + g || id),
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      g         = .label_test_factor_prior("treatment")
    ),
    prior_random = prior_random(id = random_block(
      sd = prior("normal", list(0, 1), list(0, Inf))
    ))
  )
  mixed <- as_mixed_posteriors(fit, parameters = "mu__xREx__id_g")
  expect_identical(colnames(mixed$mu__xREx__id_g), c("mu__xREx__id_g[b]", "mu__xREx__id_g[c]"))
  expect_identical(
    .label_test_legends(plot_posterior(mixed, "mu__xREx__id_g", plot_type = "ggplot")),
    list(c("b", "c"))
  )

  # the random-effect SD priors are labelled as the SDs of their block's terms
  # wherever a prior-list entry is labelled as a whole
  sd_priors <- c("mu__xREx__id_intercept", "mu__xREx__id_g")
  expect_true(all(sd_priors %in% names(attr(fit, "prior_list"))))
  models <- list(list(
    fit       = fit,
    inference = list(m_number = 1, marglik = 0, prior_prob = 1,
                     post_prob = 1, inclusion_BF = 1)
  ))
  expect_identical(
    colnames(ensemble_summary_table(models, sd_priors))[2:3],
    c("(mu) id: sd(intercept)", "(mu) id: sd(g)")
  )
  expect_identical(
    .bt_label(
      .bt_label_parts_update(
        .bt_label_parts_term("mu__xREx__id_g", attr(fit, "prior_list")$mu__xREx__id_g),
        inclusion = ""
      ),
      style = "table", formula_prefix = FALSE
    ),
    "id: sd(g) (inclusion)"
  )

  # the variable coordinates of a spike-and-slab factor term keep their
  # coordinate index, one label per coordinate
  data <- .label_test_data(c("a", "b", "c"), "treatment")
  fit <- .label_test_fit(~ x * g, data, list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1)),
    g         = prior_spike_and_slab(
      .label_test_factor_prior("treatment"),
      prior_inclusion = prior("beta", list(1, 1))
    ),
    "x:g"     = .label_test_factor_prior("treatment")
  ))
  coordinates <- parameter_coordinates(fit)
  variable <- coordinates$display_label[
    startsWith(coordinates$coordinate_name, "mu_g_variable")
  ]
  expect_identical(variable, c("(mu) g_variable[1]", "(mu) g_variable[2]"))
})
test_that("original-scale transforms require the fitted design of the formula", {

  set.seed(11)
  data <- data.frame(
    g = factor(rep(c("1", "2", "3"), 8), levels = c("1", "2", "3")),
    x = stats::rnorm(24, 3, 2)
  )
  fit <- .label_test_fit(~ g + g:x, data, list(
    intercept = prior("normal", list(0, 1)),
    g         = .label_test_factor_prior("treatment"),
    "g:x"     = .label_test_factor_prior("independent")
  ), formula_scale = list(x = TRUE))
  mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))

  # the same standardization information, without the fitted design
  formula_scale <- attr(fit, "formula_scale")
  hand_built <- list(mu = list(mu_x = list(
    mean = formula_scale$mu$mu_x$mean,
    sd   = formula_scale$mu$mu_x$sd
  )))
  message <- paste0(
    "Cannot transform the coefficients of formula parameter 'mu' to the ",
    "original predictor scale: 'formula_scale' does not carry the fitted ",
    "design of the formula. Use the 'formula_scale' attribute of the fitted ",
    "model (attr(fit, \"formula_scale\")), which carries it."
  )
  expect_error(
    ensemble_estimates_table(mixed, parameters = names(mixed),
                             transform_scaled = TRUE, formula_scale = hand_built),
    class = "BayesTools_formula_transform_unavailable"
  )
  expect_error(
    ensemble_estimates_table(mixed, parameters = names(mixed),
                             transform_scaled = TRUE, formula_scale = hand_built),
    message,
    fixed = TRUE
  )
  draws <- as.matrix(fit$mcmc)
  expect_error(
    .bt_transform_scale_posterior(draws, hand_built),
    class = "BayesTools_formula_transform_unavailable"
  )
  # the fitted object's formula_scale carries the design
  expect_s3_class(
    ensemble_estimates_table(mixed, parameters = names(mixed),
                             transform_scaled = TRUE, formula_scale = formula_scale),
    "BayesTools_table"
  )
})

test_that("inference rows and Bayes factor warnings are rendered labels", {

  data <- data.frame(x = seq(-1, 1, length.out = 12), z = sin(seq_len(12)))
  fit <- .label_test_fit(~ x + z, data, list(
    intercept = prior("normal", list(0, 1)),
    x = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(.5))
    ),
    z = prior_mixture(
      list(prior("normal", list(0, 1)), prior("normal", list(0, 2))),
      components = c("narrow", "wide")
    )
  ), n = 20L)

  for(formula_prefix in c(TRUE, FALSE)){
    table <- runjags_inference_table(fit, BF_diagnostics = TRUE,
                                     formula_prefix = formula_prefix)
    rows <- c("(mu) x", "(mu) z[narrow]", "(mu) z[wide]")
    if(!formula_prefix){
      rows <- sub("(mu) ", "", rows, fixed = TRUE)
    }
    # a mixture component is named after the term, as in estimates tables
    expect_identical(rownames(table), rows)
    # the Bayes factor MC-error warnings name the rows they belong to
    warnings <- attr(table, "warnings")
    expect_identical(names(warnings), rows)
    expect_true(all(startsWith(
      warnings,
      paste0("Bayes factor MC error for ", rows, " is based on")
    )))
  }
})

test_that("publication-weight bin rows are catalog selectors of their bins", {

  set.seed(31)
  n <- 80
  weights <- function() pmin(pmax(stats::rbeta(n, 4, 2), 1e-3), 1 - 1e-3)
  round_trip <- function(fit, expected_rows){
    catalog <- parameter_catalog(fit)
    table <- JAGS_estimates_table(fit, remove_diagnostics = TRUE)
    draws <- do.call(rbind, lapply(fit$mcmc, as.matrix))
    rows <- rownames(table)[startsWith(rownames(table), "omega[")]
    expect_identical(rows, expected_rows)
    for(row in rows){
      selection <- parameter_catalog_resolve(catalog, row)
      # the interval form is the display label; the index form stays the
      # canonical name of the same quantity
      expect_identical(selection$quantities$display_label, row)
      index_form <- selection$quantities$canonical_name
      expect_true(grepl("^omega[[][0-9]+[]]$", index_form), info = row)
      expect_identical(
        parameter_catalog_resolve(catalog, index_form)$quantities$quantity_id,
        selection$quantities$quantity_id
      )
      expect_equal(mean(draws[, index_form]), table[row, "Mean"],
                   tolerance = 1e-10, info = row)
    }
    # the mixed-posterior columns name the catalog bins holding their values
    mixed <- as_mixed_posteriors(fit, parameters = names(attr(fit, "prior_list")))
    ensemble <- ensemble_estimates_table(mixed, parameters = names(mixed))
    model_ids <- attr(table, "quantities")$quantity_id[match(rows, rownames(table))]
    ensemble_ids <- attr(ensemble, "quantities")$quantity_id[match(rows, rownames(ensemble))]
    expect_identical(ensemble_ids, model_ids)
    expect_true(all(nzchar(ensemble_ids)))
  }

  # a one-sided weight function
  chain <- function() cbind(mu = stats::rnorm(n), "omega[1]" = 1, "omega[2]" = weights())
  round_trip(
    .label_test_mock_fit(list(chain(), chain()), list(
      mu    = prior("normal", list(0, 1)),
      omega = prior_weightfunction("one-sided", .05, wf_cumulative(c(2, 4)))
    )),
    c("omega[0,0.05]", "omega[0.05,1]")
  )

  # a composed two-sided bias prior on the one-sided cut grid
  chain <- function() cbind(
    mu = stats::rnorm(n), "omega[1]" = 1, "omega[2]" = weights(), "omega[3]" = 1
  )
  round_trip(
    .label_test_mock_fit(list(chain(), chain()), list(
      mu   = prior("normal", list(0, 1)),
      bias = prior_bias(selection = prior_weightfunction(
        "two-sided", .05, wf_cumulative(c(1, 1))
      ))
    )),
    c("omega[0,0.025]", "omega[0.025,0.975]", "omega[0.975,1]")
  )

  # a publication-bias mixture of two-sided weight functions
  chain <- function() cbind(
    mu = stats::rnorm(n), bias_indicator = rep(1:3, length.out = n),
    "omega[1]" = 1, "omega[2]" = weights(), "omega[3]" = weights(),
    "omega[4]" = weights(), "omega[5]" = 1
  )
  round_trip(
    .label_test_mock_fit(list(chain(), chain()), list(
      mu   = prior("normal", list(0, 1)),
      bias = prior_mixture(list(
        prior_none(prior_weights = 1),
        prior_weightfunction("two-sided", .05, wf_cumulative(c(1, 1)), prior_weights = 1),
        prior_weightfunction("two-sided", c(.05, .10), wf_cumulative(c(1, 1, 1)), prior_weights = 1)
      ), is_null = c(TRUE, FALSE, FALSE))
    )),
    c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.95]",
      "omega[0.95,0.975]", "omega[0.975,1]")
  )
})

test_that("selection-prior coordinates declare the constants of their branches as atoms", {

  set.seed(3)
  n <- 60
  atoms_of <- function(fit, name){
    selection <- parameter_catalog_resolve(parameter_catalog(fit), name)
    atoms <- posterior_metadata(parameter_mixed_posterior(fit, selection), "atoms")
    list(x = as.numeric(atoms$locations[, 1L]), mass = atoms$mass)
  }

  # a single one-sided weight function: the reference bin is the constant 1
  # (atom of mass 1), the estimated bin is continuous
  chain <- function() cbind(mu = stats::rnorm(n), "omega[1]" = 1, "omega[2]" = stats::rbeta(n, 4, 2))
  single <- .label_test_mock_fit(list(chain(), chain()), list(
    mu    = prior("normal", list(0, 1)),
    omega = prior_weightfunction("one-sided", .05, wf_cumulative(c(2, 4)))
  ))
  expect_identical(atoms_of(single, "omega[1]"), list(x = 1, mass = 1))
  expect_identical(atoms_of(single, "omega[2]"), list(x = numeric(), mass = numeric()))
  states <- parameter_gate_states(single, parameter_catalog_resolve(parameter_catalog(single), "omega[1]"))
  expect_identical(states$prior_atoms, data.frame(x = 1, mass = 1))

  # a publication-bias mixture of no bias (weight 2), a two-sided weight
  # function, PET, and p-hacking (weight 1 each): the weights are 1 in the
  # branches without a selection and on the weight function's reference
  # bins, PET is 0 outside its branch, and alpha is 0 without p-hacking; the
  # atom masses are the shares of the bias_indicator draws in those branches
  # and the prior masses the normalized branch weights
  chain <- function(){
    indicator <- rep(1:4, length.out = n)
    cbind(
      mu = stats::rnorm(n), bias_indicator = indicator,
      "omega[1]" = 1, "omega[2]" = ifelse(indicator == 2, stats::rbeta(n, 4, 2), 1),
      "omega[3]" = 1,
      PET = ifelse(indicator == 3, abs(stats::rnorm(n)), 0),
      alpha = ifelse(indicator == 4, stats::rbeta(n, 2, 2), 0)
    )
  }
  chains <- list(chain(), chain())
  mixture <- .label_test_mock_fit(chains, list(
    mu   = prior("normal", list(0, 1)),
    bias = prior_mixture(list(
      prior_none(prior_weights = 2),
      prior_weightfunction("two-sided", .05, wf_cumulative(c(1, 1)), prior_weights = 1),
      prior_PET("normal", list(0, 1), list(0, Inf), prior_weights = 1),
      prior_phacking(form = "linear", prior_weights = 1)
    ), is_null = c(TRUE, FALSE, FALSE, FALSE))
  ))
  indicator <- unlist(lapply(chains, function(x) x[, "bias_indicator"]))
  share <- function(branches) mean(indicator %in% branches)
  expect_identical(atoms_of(mixture, "omega[1]"), list(x = 1, mass = 1))
  expect_identical(atoms_of(mixture, "omega[2]"), list(x = 1, mass = share(c(1, 3, 4))))
  expect_identical(atoms_of(mixture, "omega[3]"), list(x = 1, mass = 1))
  expect_identical(atoms_of(mixture, "PET"), list(x = 0, mass = share(c(1, 2, 4))))
  expect_identical(atoms_of(mixture, "alpha"), list(x = 0, mass = share(1:3)))
  # the declared masses are the shares of the draws on the atoms
  draws <- do.call(rbind, chains)
  expect_equal(mean(draws[, "omega[2]"] == 1), share(c(1, 3, 4)))
  expect_equal(mean(draws[, "PET"] == 0), share(c(1, 2, 4)))
  catalog <- parameter_catalog(mixture)
  prior_atoms <- function(name){
    parameter_gate_states(mixture, parameter_catalog_resolve(catalog, name))$prior_atoms
  }
  expect_equal(prior_atoms("omega[2]"), data.frame(x = 1, mass = 4 / 5))
  expect_equal(prior_atoms("PET"), data.frame(x = 0, mass = 4 / 5))
  expect_equal(prior_atoms("alpha"), data.frame(x = 0, mass = 4 / 5))
  expect_identical(
    parameter_gate_states(mixture, parameter_catalog_resolve(catalog, "PET"))$atom,
    ifelse(indicator == 3, NA_real_, 0)
  )
})

# Reference: the validator without its memo, which is the pre-memo behaviour.
test_that("label parts are validated once per content and modifications are rechecked", {

  .BayesTools_private$content_memo <- NULL
  withr::defer(.BayesTools_private$content_memo <- NULL)
  original <- .bt_validate_label_parts_uncached
  calls <- 0L
  testthat::local_mocked_bindings(
    .bt_validate_label_parts_uncached = function(parts){
      calls <<- calls + 1L
      original(parts)
    },
    .package = "BayesTools"
  )

  # building the parts validates them; rendering them again does not
  parts <- .bt_label_parts(c("mu", "x"), formula_parameter = "mu", levels = c(mu = "a"))
  expect_identical(calls, 1L)
  labels <- vapply(c("selector", "table", "plot", "warning"), function(style){
    .bt_label(parts, style = style)
  }, character(1))
  expect_identical(calls, 1L)
  expect_identical(.bt_label(unserialize(serialize(parts, NULL)), style = "table"), labels[["table"]])
  expect_identical(calls, 1L)

  # every other content is checked (once), and a refusal is never remembered
  updated <- .bt_label_parts_update(parts, transformation = "dif")[[1L]]
  expect_identical(calls, 2L)
  .bt_label(updated)
  expect_identical(calls, 2L)
  malformed <- list(
    parts_field = {
      out <- parts
      out$components <- character()
      out
    },
    transformation = {
      out <- parts
      out$transformation <- "unknown"
      out
    },
    random = {
      out <- parts
      out$random <- list(owner = "g")
      out
    },
    class = unclass(parts)
  )
  for(name in names(malformed)){
    for(i in 1:2){
      memoized <- tryCatch(.bt_validate_label_parts(malformed[[name]]), error = identity)
      unmemoized <- tryCatch(original(malformed[[name]]), error = identity)
      expect_s3_class(memoized, "error")
      expect_identical(conditionMessage(memoized), conditionMessage(unmemoized))
    }
  }
  expect_identical(calls, 2L + 2L * length(malformed))
  # the parts that were valid before the refusals still are, without a check
  before <- calls
  .bt_label(parts)
  expect_identical(calls, before)
})
