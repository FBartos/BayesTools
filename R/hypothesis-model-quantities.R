# ============================================================================ #
# hypothesis-model-quantities.R
# ============================================================================ #
#
# The eligibility table of hypothesis_quantities(): which tests the engine
# runs on each quantity of a model, rendered from the plans of the statements
# '<q> = <null>', '<q> > <null>', and '<level> = <other level>'.
#
# ============================================================================ #


# The table of a model without quantities that the engine tests: the columns
# of hypothesis_quantities() without rows.
.hypothesis_model_quantities_empty <- function() {

  data.frame(
    alias                 = character(),
    parameter             = character(),
    component             = character(),
    term                  = character(),
    bracket               = character(),
    point_test            = logical(),
    direction_test        = logical(),
    contrast_test         = logical(),
    point_test_methods    = character(),
    contrast_test_methods = character(),
    reason                = character(),
    stringsAsFactors      = FALSE
  )
}


# One row of eligibility columns.
.hypothesis_model_quantities_row <- function(bracket = FALSE, point_test = FALSE,
                                       direction_test = FALSE,
                                       contrast_test = NA,
                                       point_test_methods = "",
                                       contrast_test_methods = NA_character_,
                                       reason = "") {

  data.frame(
    bracket               = bracket,
    point_test            = point_test,
    direction_test        = direction_test,
    contrast_test         = contrast_test,
    point_test_methods    = point_test_methods,
    contrast_test_methods = contrast_test_methods,
    reason                = reason,
    stringsAsFactors      = FALSE
  )
}


# The rendered statements of a catalog entry and their plans. Statements
# name the quantity as .hypothesis_model_quantities_reference_root() does.
.hypothesis_model_quantities_plans <- function(provider, entry, metadata,
                                               cache, transform_scaled) {

  # Statements are planned as hypothesis_BF_model() plans them.
  plan_of <- function(statement) {
    .hypothesis_model_plans(
      provider         = provider,
      hypothesis       = statement,
      component        = entry[["component"]],
      transform_scaled = transform_scaled,
      metadata         = metadata,
      cache            = cache
    )[[1L]]
  }
  root <- .hypothesis_model_quantities_reference_root(metadata, entry)
  reference <- function(level = NULL) {
    if (is.null(level)) {
      return(paste0("`", root, "`"))
    }
    paste0("`", root, "[", level, "]`")
  }

  if (!identical(entry[["role"]], "formula_coefficient_group")) {
    point <- .hypothesis_model_quantities_point_plan(plan_of, reference(), provider)
    return(list(
      point    = list(point),
      region   = list(plan_of(paste0(
        reference(), " > ", .hypothesis_model_quantities_region_value(point)
      ))),
      contrast = NULL
    ))
  }

  levels <- .hypothesis_model_plan_term_levels(metadata, entry)
  tested <- levels[["level"]][!levels[["fixed"]]]
  point  <- lapply(tested, function(level) {
    .hypothesis_model_quantities_point_plan(plan_of, reference(level), provider)
  })
  region <- lapply(seq_along(tested), function(i) {
    plan_of(paste0(
      reference(tested[[i]]), " > ",
      .hypothesis_model_quantities_region_value(point[[i]])
    ))
  })
  pairs <- if (nrow(levels) > 1L) utils::combn(levels[["level"]], 2L) else NULL
  contrast <- lapply(seq_len(NCOL(pairs)), function(i) {
    plan_of(paste0(reference(pairs[1L, i]), " = ", reference(pairs[2L, i])))
  })

  list(point = point, region = region, contrast = contrast)
}


# The name of a catalog entry's quantity in its rendered statements.
# Random-effect quantities and factor terms are named by their parameter (the
# level references of factor terms are resolved across components). Other
# quantities are named by their displayed alias within the entry's component,
# unless that alias also names another quantity of the component (e.g. the
# 'intercept' of several scale formulas): such quantities are named by their
# catalog selector.
.hypothesis_model_quantities_reference_root <- function(metadata, entry) {

  parameter <- entry[["parameter"]]
  if (identical(entry[["component"]], "random") ||
      identical(entry[["role"]], "formula_coefficient_group")) {
    return(parameter)
  }

  entries <- metadata[["entries"]]
  aliases <- entries[["aliases"]][[
    match(entry[["quantity_id"]], entries[["quantity_id"]])
  ]]
  aliases <- as.list(stats::setNames(rep(parameter, length(aliases)), aliases))
  label   <- .hypothesis_model_alias_label(aliases, parameter)
  if (.hypothesis_model_quantities_alias_resolves(
    catalog     = metadata[["catalog"]],
    alias       = label,
    component   = entry[["component"]],
    quantity_id = entry[["quantity_id"]]
  )) {
    return(label)
  }

  quantities <- metadata[["catalog"]][["quantities"]]
  parameter_labels(
    quantities[
      match(entry[["quantity_id"]], quantities[["quantity_id"]]), ,
      drop = FALSE
    ],
    style = "selector"
  )
}


# Whether 'alias' names the quantity 'quantity_id' within 'component': the
# resolution the engine applies to the statements of that component. An
# alias that also names another quantity of the component is ambiguous there
# (e.g. the 'intercept' of several scale formulas).
.hypothesis_model_quantities_alias_resolves <- function(catalog, alias, component,
                                                  quantity_id) {

  tryCatch(
    identical(
      parameter_catalog_resolve(
        catalog        = catalog,
        alias          = alias,
        component      = component,
        simplify_names = TRUE
      )[["quantity_id"]],
      quantity_id
    ),
    BayesTools_parameter_ambiguous = function(error) FALSE
  )
}


# The plan of '<q> = 0', or of '<q> = <interior value>' when the value 0
# itself is refused (a prior point mass, a support boundary without a
# regular prior ordinate, or a value outside the support): the statement
# rendered for a quantity's point tests.
.hypothesis_model_quantities_point_plan <- function(plan_of, reference, provider) {

  plan    <- plan_of(paste0(reference, " = 0"))
  refusal <- .hypothesis_model_quantities_status(plan, provider, "KDE")
  if (!is.null(refusal) && isTRUE(refusal[["value_specific"]])) {
    value <- .hypothesis_model_quantities_interior_value(plan[["support"]])
    plan  <- plan_of(paste0(reference, " = ", format(value, digits = 15L)))
    plan[["null_value"]] <- value
    return(plan)
  }
  plan[["null_value"]] <- 0

  plan
}


# The value of the rendered region statement '<q> > <value>': the null value
# of the point statement, or an interior value when that null value is a
# support bound (a region beyond a bound has no complementary prior mass).
.hypothesis_model_quantities_region_value <- function(point) {

  value   <- point[["null_value"]]
  support <- point[["support"]]
  if (!is.null(support) && (value <= support[[1L]] || value >= support[[2L]])) {
    value <- .hypothesis_model_quantities_interior_value(support)
  }

  format(value, digits = 15L)
}


# A value inside the support bounds 'support' (c(lower, upper), or NULL for
# an unbounded support): the midpoint of a bounded support, one unit inside a
# half-bounded support, and 1 otherwise.
.hypothesis_model_quantities_interior_value <- function(support) {

  bounds <- if (is.null(support)) c(-Inf, Inf) else as.numeric(support)
  if (all(is.finite(bounds))) {
    return(mean(bounds))
  }
  if (is.finite(bounds[[1L]])) {
    return(bounds[[1L]] + 1)
  }
  if (is.finite(bounds[[2L]])) {
    return(bounds[[2L]] - 1)
  }

  1
}


.hypothesis_model_quantities_render <- function(provider, entry, metadata,
                                                cache, transform_scaled) {

  plans <- .hypothesis_model_quantities_plans(
    provider, entry, metadata, cache, transform_scaled
  )
  .hypothesis_model_quantities_render_plans(
    plans    = plans,
    bracket  = identical(entry[["role"]], "formula_coefficient_group"),
    provider = provider
  )
}


# The status of the named density method for a plan, with the default
# settings of the method (computed on every call).
.hypothesis_model_quantities_status <- function(plan, provider, method) {

  backend <- .hypothesis_model_backend(method)

  .hypothesis_model_plan_status(
    plan, provider, backend,
    control = .hypothesis_model_backend_control(backend, NULL)
  )
}


# The methods that hypothesis_quantities() lists: the registered backends that
# are advertised and offered for the model (the normal approximation is a
# rough check of near-normal interior ordinates and is not advertised).
.hypothesis_model_quantities_methods <- function(provider) {

  backends <- .hypothesis_model_backends()
  allowed  <- .hypothesis_model_density_methods(provider)
  names(backends)[vapply(backends, function(backend) {
    isTRUE(backend[["advertised"]]) && backend[["name"]] %in% allowed
  }, logical(1))]
}


# The eligibility columns of the plans of one quantity: the methods that
# evaluate every point (contrast) plan, whether every region plan runs, and
# the distinct refusals of the plans.
.hypothesis_model_quantities_render_plans <- function(plans, bracket, provider) {

  advertised <- .hypothesis_model_quantities_methods(provider)
  methods_of <- function(plans) {
    if (length(plans) == 0L) {
      return(character())
    }
    advertised[vapply(advertised, function(method) {
      all(vapply(plans, function(plan) {
        is.null(.hypothesis_model_quantities_status(plan, provider, method))
      }, logical(1)))
    }, logical(1))]
  }
  point_methods    <- methods_of(plans[["point"]])
  contrast_methods <- if (is.null(plans[["contrast"]])) {
    NULL
  } else {
    methods_of(plans[["contrast"]])
  }
  direction_test <- length(plans[["region"]]) > 0L && all(vapply(
    plans[["region"]],
    function(plan) is.null(.hypothesis_model_quantities_status(plan, provider, "KDE")),
    logical(1)
  ))
  reasons <- unlist(lapply(
    c(plans[["point"]], plans[["region"]], plans[["contrast"]]),
    function(plan) {
      vapply(advertised, function(method) {
        refusal <- .hypothesis_model_quantities_status(plan, provider, method)
        if (is.null(refusal)) NA_character_ else refusal[["reason"]]
      }, character(1))
    }
  ), use.names = FALSE)
  reasons <- unique(reasons[!is.na(reasons)])

  .hypothesis_model_quantities_row(
    bracket               = bracket,
    point_test            = length(point_methods) > 0L,
    direction_test        = direction_test,
    contrast_test         = if (is.null(contrast_methods)) {
      NA
    } else {
      length(contrast_methods) > 0L
    },
    point_test_methods    = paste(point_methods, collapse = ", "),
    contrast_test_methods = if (is.null(contrast_methods)) {
      NA_character_
    } else {
      paste(contrast_methods, collapse = ", ")
    },
    reason                = paste(reasons, collapse = " ")
  )
}


# The eligibility table of the marginal means stored by a provider of kind
# "marginal_means": one row per level of each term.
.hypothesis_model_quantities_marginal <- function(provider) {

  term_map <- provider[["metadata"]]()[["term_map"]]
  averaged <- provider[["inference"]][["averaged"]]
  if (is.null(averaged)) {
    averaged <- provider[["inference"]][["conditional"]]
  }
  rows  <- list()
  cache <- .hypothesis_model_plan_cache(provider)
  for (i in seq_len(nrow(term_map))) {

    parameter <- term_map[["parameter"]][[i]]
    levels <- names(averaged[[parameter]])
    if (is.null(levels) ||
        length(levels) != length(averaged[[parameter]]) ||
        any(!nzchar(levels))) {
      stop(
        "Grouped marginal means must have non-empty level names.",
        call. = FALSE
      )
    }
    for (level in levels) {
      rows[[length(rows) + 1L]] <- cbind(
        .hypothesis_model_quantities_marginal_row(term_map[i, , drop = FALSE], level),
        .hypothesis_model_quantities_marginal_render(
          provider    = provider,
          parameter = parameter,
          level     = level,
          levels    = levels,
          cache     = cache
        )
      )
    }
  }
  if (length(rows) == 0L) {
    return(.hypothesis_model_quantities_empty())
  }
  out <- do.call(rbind, rows)
  out <- out[, c(
    "alias", "parameter", "component", "term", "bracket", "point_test",
    "direction_test", "contrast_test", "point_test_methods",
    "contrast_test_methods", "reason"
  ), drop = FALSE]
  rownames(out) <- NULL

  return(out)
}


.hypothesis_model_quantities_marginal_row <- function(term_row, level) {

  parameter <- term_row[["parameter"]][[1L]]
  return(data.frame(
    alias      = term_row[["term"]][[1L]],
    parameter  = parameter,
    component  = "marginal_means",
    term       = term_row[["term"]][[1L]],
    bracket    = paste0(parameter, "[", level, "]"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
}


# The plans of one marginal-mean level: its point and region statements and
# its contrasts with the other levels.
.hypothesis_model_quantities_marginal_render <- function(provider, parameter, level,
                                                   levels, cache) {

  plan_of <- function(statement) {
    .hypothesis_model_plan_marginal_means(
      provider    = provider,
      statement = statement,
      parameter = parameter,
      cache     = cache
    )
  }
  reference <- function(level) {
    paste0("`", parameter, "[", level, "]`")
  }
  point  <- .hypothesis_model_quantities_point_plan(plan_of, reference(level), provider)
  plans  <- list(
    point    = list(point),
    region   = list(plan_of(paste0(
      reference(level), " > ", .hypothesis_model_quantities_region_value(point)
    ))),
    contrast = if (length(levels) > 1L) {
      lapply(setdiff(levels, level), function(other) {
        plan_of(paste0(reference(level), " = ", reference(other)))
      })
    }
  )
  out <- .hypothesis_model_quantities_render_plans(plans, bracket = FALSE, provider)
  out[["bracket"]] <- NULL

  out
}


#' @title Available Hypothesis Quantities
#'
#' @description Lists the parameter names and aliases that
#' [hypothesis_BF()] accepts for a fitted model, with the tests that it runs
#' on each of them.
#'
#' @param object a fitted model with a [hypothesis_provider()] method, for
#'   example a \code{BayesTools_fit}.
#' @param transform_scaled whether the rendered statements test coefficients of
#'   scaled predictors on the original scale (\code{TRUE}) or the fitted scale
#'   (\code{FALSE}). The default \code{NULL} is the default of the provider
#'   (\code{FALSE} for a \code{BayesTools_fit}).
#' @param ... unused; unused arguments warn.
#'
#' @details The eligibility columns render the plans that
#' [hypothesis_BF_model()] executes (see its details) for the statements
#' \code{<q> = <null>} (point tests), \code{<q> > <null>} (directional and
#' interval tests), and, for factor terms, \code{<level> = <other level>}
#' (level contrasts) of every quantity \code{q}. The null value is 0; when 0
#' is a prior point mass, a support boundary without a regular prior
#' ordinate, or outside the support (for example, a gated random-effect
#' standard deviation or a variance proportion), the point test is rendered
#' at an interior value of the support instead (the midpoint of a bounded
#' support, one unit inside a half-bounded support, and 1 otherwise): point
#' hypotheses at such values follow the per-value rules of
#' [hypothesis_BF_model()]. A factor term lists the point-test methods
#' available for all of its levels that the contrast does not fix, and the
#' contrast-test methods available for all pairs of its levels. Region tests
#' do not use a density method. The methods are the registered density
#' backends that are advertised ([hypothesis_density_backends()]; the normal
#' approximation is a rough check and is not listed).
#'
#' @return A data frame with the columns \code{alias} (a name that
#' [hypothesis_BF()] resolves to the row's quantity with the row's
#' \code{component}; a name shared by several quantities of a component is not
#' listed), \code{parameter}, \code{component}, \code{term}, \code{bracket} (a
#' \code{parameter[level]} template for fitted factor terms, or a concrete
#' level selector for each marginal-mean row), \code{point_test},
#' \code{direction_test}, \code{contrast_test} (\code{NA} for quantities
#' without levels), \code{point_test_methods}, \code{contrast_test_methods},
#' and \code{reason} (the refusals of the rendered statements).
#'
#' @seealso [hypothesis_BF()], [hypothesis_BF_model()]
#' @export
hypothesis_quantities <- function(object, ...){

  UseMethod("hypothesis_quantities")
}


#' @rdname hypothesis_quantities
#' @exportS3Method hypothesis_quantities default
hypothesis_quantities.default <- function(object, transform_scaled = NULL, ...){

  .hypothesis_model_warn_unused_dots(
    dots    = list(...),
    allowed = character(),
    caller  = "hypothesis_quantities()"
  )
  provider <- hypothesis_provider(object)
  .hypothesis_model_validate_provider(provider)
  if(is.null(transform_scaled)){
    transform_scaled <- isTRUE(provider[["transform_scaled_default"]])
  }
  check_bool(transform_scaled, "transform_scaled", allow_NA = FALSE)
  if(identical(provider[["kind"]], "marginal_means")){
    return(.hypothesis_model_quantities_marginal(provider))
  }

  .hypothesis_model_quantities_model(provider, transform_scaled)
}


#' @rdname hypothesis_quantities
#' @exportS3Method hypothesis_quantities BayesTools_fit
hypothesis_quantities.BayesTools_fit <- function(object, transform_scaled = NULL, ...){

  hypothesis_quantities.default(object, transform_scaled = transform_scaled, ...)
}


# The entries of the metadata of a provider, one row per alias.
.hypothesis_model_quantities_catalog <- function(metadata){

  entries <- metadata[["entries"]]
  rows <- lapply(seq_len(nrow(entries)), function(i){
    aliases <- entries[["aliases"]][[i]]
    data.frame(
      alias             = aliases,
      parameter         = entries[["parameter"]][[i]],
      component         = entries[["component"]][[i]],
      term              = entries[["term"]][[i]],
      source            = entries[["source"]][[i]],
      formula_parameter = entries[["formula_parameter"]][[i]],
      quantity_id       = entries[["quantity_id"]][[i]],
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(NULL)
  }
  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL

  out
}


.hypothesis_model_quantities_model <- function(provider, transform_scaled){

  metadata <- provider[["metadata"]]()
  catalog  <- .hypothesis_model_quantities_catalog(metadata)
  if (is.null(catalog)) {
    return(.hypothesis_model_quantities_empty())
  }
  catalog  <- catalog[
    !catalog[["component"]] %in% names(provider[["components"]][["refused"]]),
    , drop = FALSE
  ]
  if (nrow(catalog) == 0L) {
    return(.hypothesis_model_quantities_empty())
  }
  out  <- catalog[, c("alias", "parameter", "component", "term"), drop = FALSE]
  keys <- catalog[["quantity_id"]]
  unique_rows <- match(unique(keys), keys)

  entries <- metadata[["entries"]]
  cache   <- .hypothesis_model_plan_cache(provider)
  rows <- lapply(unique_rows, function(i){
    entry <- as.list(entries[
      match(catalog[["quantity_id"]][[i]], entries[["quantity_id"]]),
      setdiff(names(entries), "aliases"),
      drop = FALSE
    ])
    .hypothesis_model_quantities_render(
      provider         = provider,
      entry            = entry,
      metadata         = metadata,
      cache            = cache,
      transform_scaled = transform_scaled
    )
  })
  rows <- do.call(rbind, rows)[match(keys, unique(keys)), , drop = FALSE]
  out  <- cbind(out, rows)
  # Only aliases that name their own quantity within its component, as the
  # engine resolves them; ambiguous aliases are not listed.
  resolves <- vapply(seq_len(nrow(catalog)), function(i){
    .hypothesis_model_quantities_alias_resolves(
      catalog     = metadata[["catalog"]],
      alias       = catalog[["alias"]][[i]],
      component   = catalog[["component"]][[i]],
      quantity_id = catalog[["quantity_id"]][[i]]
    )
  }, logical(1))
  out <- out[resolves, , drop = FALSE]
  out[["bracket"]] <- as.character(ifelse(
    out[["bracket"]],
    paste0(out[["parameter"]], "[level]"),
    NA_character_
  ))
  out <- out[, c(
    "alias", "parameter", "component", "term", "bracket", "point_test",
    "direction_test", "contrast_test", "point_test_methods",
    "contrast_test_methods", "reason"
  ), drop = FALSE]
  rownames(out) <- NULL

  out
}
