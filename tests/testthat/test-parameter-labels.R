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
  .parameter_catalog_test_fit(
    coda::mcmc.list(coda::mcmc(draws)),
    prior_list     = formula_result$prior_list,
    formula_design = stats::setNames(list(formula_result$formula_design), parameter),
    formula_scale  = if(length(scale) > 0L) stats::setNames(list(scale), parameter)
  )
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
    "tick`", "tab\there", "line\nbreak", "été"
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
