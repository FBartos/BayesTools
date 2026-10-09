# ============================================================================ #
# posterior-transform.R
# ============================================================================ #
#
# A monotone transformation of posterior draws together with the metadata
# that describe them (R/draws-metadata.R): values, exact supports, atoms,
# mixture-component supports, prior densities, precomputed posterior
# densities and ordinates, linear-combination weights, and the label parts of
# the quantities. marginal_posterior(transformation = ) builds the
# untransformed marginal posterior and transforms it here.
#
# ============================================================================ #

#' @title Transform posterior draws with their metadata
#'
#' @description Applies a strictly monotone transformation to BayesTools
#' posterior draws (an element returned by [marginal_posterior()],
#' [parameter_mixed_posterior()], [as_mixed_posteriors()], or
#' [mix_posteriors()], or a list of them) and transforms the draw metadata
#' ([posterior_metadata()]) with the values, so that inference, plots, and
#' summaries of the transformed draws remain valid. Arithmetic on draws
#' returns plain numeric draws without metadata; this function is the route
#' for transformed posterior distributions.
#'
#' @param x BayesTools posterior draws with metadata, or a list of them.
#' @param transformation a transformation in the form accepted by
#' [marginal_posterior()]: \code{"lin"} (\code{a + b * x}), \code{"exp_lin"}
#' (\code{exp(a + b * log(x))}), \code{"exp"}, \code{"tanh"}, or a list of
#' functions \code{fun}, \code{inv}, and \code{jac} (the transformation, its
#' inverse, and its derivative), optionally with \code{output_support}.
#' @param transformation_arguments optional named list of arguments of the
#' transformation (\code{a} and \code{b} of \code{"lin"} and
#' \code{"exp_lin"}).
#'
#' @details The transformation must be strictly monotone and invertible.
#' \code{"lin"} and \code{"exp_lin"} with \code{b = 0} are constant, and a
#' transformation given as functions must have a finite, nonzero derivative
#' of one sign at every draw and atom; otherwise the function stops with an
#' error of class \code{BayesTools_nonmonotone_transformation}. Draws that
#' the transformation maps to non-finite values (for example, negative draws
#' of \code{"exp_lin"}) stop with class
#' \code{BayesTools_transformation_domain}. Both classes have the parent class
#' \code{BayesTools_transformation}.
#' Finite builtin images that round to mathematical limits (exponential zero,
#' tanh endpoints, or positive-power zero from a strictly positive source) stop
#' with \code{BayesTools_transformation_image_unavailable}, also inheriting
#' \code{BayesTools_transformation}, \code{error} and \code{condition}. Its
#' \code{call} is \code{NULL}; \code{transformation} names the builtin or
#' \code{"custom"}, and \code{source_values} and \code{images} contain matched
#' failing finite entries. Declared point and support mapping uses the same
#' condition for nonrepresentable images. True source zero under a positive
#' power and infinite mathematical support limits remain valid. No clipping or
#' probability-mass repair is performed.
#'
#' The metadata are transformed as follows:
#' A missing canonical prior law is first built from the complete declared
#' source context and recipe, preserving its conditioning and coefficient or
#' formula-contribution space. It is then transformed once. An attached law
#' on the current output scale takes precedence over the original source
#' context in plots. Known unsupported or numerical prior-law routes retain
#' usable posterior draws and declare the prior density unavailable; requesting
#' a prior overlay raises \code{BayesTools_formula_prior_density_unavailable}.
#' Transformed publication-bias scalar prior overlays require a supported
#' current-scale law and otherwise have the same explicit limitation.
#' \describe{
#'   \item{\code{support}}{the bounds and points are mapped, and the bounds
#'   of decreasing transformations swapped. The support of a transformation
#'   given as functions is unknown and removed.}
#'   \item{\code{atoms}}{the locations are mapped; the masses are kept.}
#'   \item{\code{prior_density}, \code{prior_densities}}{change of variables
#'   through the transformation machinery of the prior densities: a density
#'   built from a prior-density context or a linear combination of priors is
#'   rebuilt with the transformation as its output transformation (the
#'   density [marginal_posterior()] returns with \code{transformation}), a
#'   density that already has an output transformation or is an allocation
#'   product is transformed on top of its recorded provenance, and a density
#'   grid without provenance stays a grid for plotting only.}
#'   \item{\code{posterior_density}, \code{posterior_ordinate}}{the locations
#'   are mapped and the heights divided by the absolute derivative of the
#'   transformation; a density or ordinate whose transformed heights are not
#'   finite and positive (e.g. at a point where the derivative vanishes) is
#'   removed. Diagnostics are kept; their null-value entries
#'   (\code{bf_value}, \code{value}, \code{null_hypothesis}) are mapped.}
#'   \item{\code{components}}{the per-component supports are transformed.}
#'   \item{\code{linear_weights}}{\code{"lin"} scales the weights by
#'   \code{b} and adds \code{a} to the offset (\code{linear_offset}); any
#'   other transformation removes them and records its name in the joint
#'   prior transformation.}
#'   \item{\code{quantities}}{the name of the transformation (\code{"custom"}
#'   for functions) is appended to the \code{transformation} of the label
#'   parts; the columns are no longer the catalog quantities (the quantity id
#'   is \code{""}) and declare fitted coordinates only for \code{"lin"} with
#'   \code{a = 0}, whose weights are scaled by \code{b}. The rendered labels
#'   are unchanged.}
#' }
#' The conditioning (\code{condition}), undefined draws, and the other
#' metadata are kept, the applied transformation is appended to the
#' \code{output_transformations} metadata (whether or not the draws have a
#' \code{quantities} table), and the metadata record the fingerprint of the
#' transformed values. The \code{prior_list} attribute of mixed posteriors
#' describes their untransformed values: [marginal_posterior()] refuses
#' transformed mixed posteriors whose prior it would build from it; use its
#' \code{transformation} argument instead.
#'
#' @return \code{x} with transformed values and metadata.
#'
#' @seealso [marginal_posterior()], [posterior_metadata()]
#'
#' @examples
#' draws <- structure(stats::rnorm(100), class = c("marginal_posterior.simple", "marginal_posterior"))
#' posterior_metadata(draws, "support") <- posterior_support_attribute(c(-Inf, Inf))
#' posterior_metadata(draws, "atoms") <- posterior_atom_attribute()
#' transformed <- posterior_transform(draws, "exp")
#' posterior_metadata(transformed, "support")
#'
#' @export
posterior_transform <- function(x, transformation, transformation_arguments = NULL){

  if(missing(transformation) || is.null(transformation)){
    stop("'transformation' must be specified.", call. = FALSE)
  }
  .check_transformation_input(transformation, transformation_arguments, FALSE)
  map <- .bt_posterior_transformation(transformation, transformation_arguments)

  .bt_posterior_transform(x, map)
}

# The transformation 'map' of posterior_transform(): its functions, its name
# ("custom" for a list of functions), and its direction (1 increasing, -1
# decreasing, NA when only the draws can tell).
.bt_posterior_transformation <- function(transformation,
                                         transformation_arguments = NULL){

  functions <- .density.prior_transformation_functions(transformation)
  if(is.null(functions)){
    stop(
      "Unknown transformation '", transformation, "'. Use 'lin', 'exp_lin', ",
      "'tanh', 'exp', or a list of functions 'fun', 'inv', and 'jac'.",
      call. = FALSE
    )
  }
  name <- if(is.character(transformation)) transformation else "custom"
  direction <- NA_real_
  if(name %in% c("lin", "exp_lin")){
    a <- .posterior_support_transform_argument(transformation_arguments, "a", 0)
    b <- .posterior_support_transform_argument(transformation_arguments, "b", 1)
    if(!is.finite(a) || !is.finite(b)){
      stop("The arguments 'a' and 'b' of the '", name, "' transformation must be ",
           "finite numbers.", call. = FALSE)
    }
    if(b == 0){
      .bt_posterior_transformation_stop(
        "nonmonotone",
        paste0("the '", name, "' transformation with 'b = 0' is constant")
      )
    }
    direction <- sign(b)
  }else if(name %in% c("exp", "tanh")){
    direction <- 1
  }

  list(
    transformation = transformation,
    arguments      = transformation_arguments,
    name           = name,
    direction      = direction,
    fun            = function(x){
      .density.prior_transformation_x(x, transformation, transformation_arguments)
    },
    jac            = function(x){
      arguments <- c(list(x = x), transformation_arguments)
      out <- do.call(functions$jac, arguments)
      rep_len(as.numeric(out), length(x))
    }
  )
}

.bt_posterior_transformation_stop <- function(type, detail){

  if(identical(type, "nonmonotone")){
    message <- paste0(
      "The transformation of the posterior draws is unavailable: it must be ",
      "strictly monotone and invertible, but ", detail, "."
    )
    class <- "BayesTools_nonmonotone_transformation"
  }else{
    message <- paste0(
      "The transformation of the posterior draws is unavailable: ", detail, "."
    )
    class <- "BayesTools_transformation_domain"
  }
  stop(errorCondition(
    message,
    class = c(class, "BayesTools_transformation"),
    call  = NULL
  ))
}

# Checks the transformation on the finite values it is applied to (draws and
# atom locations): the images must be finite, and a transformation given as
# functions must have a finite nonzero derivative of one sign there.
.bt_posterior_transformation_check <- function(map, values){

  values <- values[is.finite(values)]
  if(length(values) == 0L){
    return(invisible(map))
  }
  images <- suppressWarnings(map$fun(values))
  if(any(!is.finite(images))){
    .bt_posterior_transformation_stop(
      "domain",
      paste0(
        "the '", map$name, "' transformation maps the value ",
        format(values[!is.finite(images)][[1L]]), " to a non-finite value"
      )
    )
  }
  arguments <- .density.prior_transformation_named_arguments(map$transformation, map$arguments)
  bad <- .density.prior_transformation_image_bad(values, images, map$name, arguments)
  if(any(bad)){
    .density.prior_transformation_image_stop(map$transformation, values[bad], images[bad])
  }
  if(is.na(map$direction)){
    jacobian <- suppressWarnings(map$jac(values))
    if(any(!is.finite(jacobian)) || any(jacobian == 0) ||
       length(unique(sign(jacobian))) > 1L){
      .bt_posterior_transformation_stop(
        "nonmonotone",
        "its derivative 'jac' is not finite and nonzero with one sign on the draws"
      )
    }
  }

  invisible(map)
}

# Transforms draws or a list of draws.
.bt_posterior_transform <- function(x, map){

  if(is.list(x) && !is.numeric(x)){
    x_attributes <- attributes(x)
    for(i in seq_along(x)){
      name <- if(is.null(names(x))) NULL else names(x)[[i]]
      x[[i]] <- .bt_posterior_transform_draws(x[[i]], map, parent = x, parameter = name)
    }
    attributes(x) <- x_attributes
    # list-level metadata (e.g. sources of precomputed posterior densities)
    fields <- .bt_meta_get_fields(x, .bt_meta_fields())
    updates <- .bt_posterior_transform_fields(fields, map)
    if(length(updates) > 0L){
      x <- .bt_meta_assign(x, updates)
    }
    return(x)
  }

  .bt_posterior_transform_draws(x, map)
}

.bt_posterior_transform_source_law <- function(x, parent = NULL, parameter = NULL){

  fields <- .bt_meta_get_fields(x, .bt_meta_fields())
  if(!is.null(fields$prior_density) || !is.null(fields$prior_densities) ||
     (!is.null(fields$measure_unavailable) && any(fields$measure_unavailable$measure == "prior_density"))) return(x)
  if(.plot_data_samples_without_prior(x)) return(x)
  parent_fields <- if(is.null(parent)) list() else .bt_meta_get_fields(parent, .bt_meta_fields())
  promise_fields <- c("prior_context", "linear_weights", "formula_state", "formula_parameter", "formula_scale",
    "ordered_source", "model_probabilities")
  bare <- is.null(attr(x, "prior_list", exact = TRUE)) &&
    (is.null(parent) || is.null(attr(parent, "prior_list", exact = TRUE))) &&
    all(vapply(c(fields[promise_fields], parent_fields[promise_fields]), is.null, logical(1)))
  if(bare) return(x)
  if(inherits(x, c("mixed_posteriors.weightfunction", "mixed_posteriors.bias"))){
    return(.bt_formula_measure_mark(x, colnames(x), "prior_density",
      "The transformed publication-bias scalar prior law is unavailable.",
      cause = "unsupported_contribution_measure"))
  }
  if(length(fields$output_transformations)){
    .bt_formula_density_stop("The transformed source has no current-scale prior law.",
      target = parameter, reason = "structural_target_law_unavailable")
  }
  if(is.null(parameter)) parameter <- attr(x, "parameter", exact = TRUE)
  if(is.null(parameter) && !is.null(fields$quantities) && nrow(fields$quantities) == 1L){
    parameter <- fields$quantities$column[[1L]]
  }
  if(is.null(parent)){
    parent <- stats::setNames(list(x), if(is.null(parameter)) "value" else parameter)
    if(!is.null(fields$prior_context)) parent <- .bt_meta_set(parent, "prior_context", fields$prior_context)
    if(!is.null(fields$condition)) parent <- .bt_meta_set(parent, "condition", fields$condition)
    priors <- attr(x, "prior_list", exact = TRUE)
    if(is.prior(priors) || (is.list(priors) && all(vapply(priors, is.prior, logical(1))))){
      if(is.null(parameter) && is.null(fields$prior_context)){
        .bt_formula_density_stop("The detached source prior has no declared target name.", reason = "unknown_target")
      }
      if(!is.null(parameter)) priors <- stats::setNames(list(priors), parameter)
      attr(parent, "prior_list") <- priors
    }
  }
  context <- if(!is.null(fields$prior_context)) fields$prior_context else parent_fields$prior_context
  priors <- attr(parent, "prior_list", exact = TRUE)
  if(is.null(priors) && !is.null(context)) priors <- context$prior_list
  if(is.null(priors) && is.null(context)){
    child_priors <- attr(x, "prior_list", exact = TRUE)
    if(!is.null(child_priors)){
      if(is.null(parameter) || length(parameter) != 1L){
        .bt_formula_density_stop("The declared child prior has no unique target name.", reason = "unknown_target")
      }
      priors <- stats::setNames(list(child_priors), parameter)
      attr(parent, "prior_list") <- priors
    }
  }
  if(is.null(priors) && is.null(context)){
    .bt_formula_density_stop("The source has no complete declared prior context.",
      target = parameter, reason = "missing_source_context")
  }
  if(!is.null(priors) && all(vapply(priors, is.prior.none, logical(1)))) return(x)
  build <- function(){
    source <- if(is.null(fields$prior_context)) parent else .bt_meta_set(parent, "prior_context", fields$prior_context)
    raw_coefficients <- if(is.null(fields$linear_weights)){
      !isTRUE(fields$transform_scaled)
    }else identical(.bt_linear_weight_space(x), "coefficient")
    context <- .marginal_posterior_prior_density_context(source, priors,
      column_names = context$column_names,
      n_samples = if(is.null(context$n_grid)) .prior_linear_density_default_grid() else context$n_grid,
      condition_source = x, raw_coefficients = raw_coefficients)
    weights <- fields$linear_weights
    if(is.null(weights) && is.matrix(x) && inherits(x, "mixed_posteriors.factor")){
      target_parameter <- attr(x, "parameter", exact = TRUE)
      if(is.null(target_parameter)) target_parameter <- parameter
      factor_weights <- .prior_factor_level_weight_matrix(x, target_parameter, parent)
      if(!all(colnames(factor_weights) %in% context$column_names)){
        .bt_formula_density_stop("The factor prior recipe is not aligned with its source context.",
          target = target_parameter, reason = "unknown_target")
      }
      weights <- matrix(0, nrow(factor_weights), length(context$column_names),
        dimnames = list(colnames(x), context$column_names))
      weights[, colnames(factor_weights)] <- factor_weights
    }
    if(is.null(weights)){
      columns <- if(is.matrix(x)) colnames(x) else parameter
      if(is.null(columns) || !all(columns %in% context$column_names)){
        .bt_formula_density_stop("The source prior context has no declared target recipe.",
          target = parameter, reason = "unknown_target")
      }
      weights <- diag(length(context$column_names))[match(columns, context$column_names), , drop = FALSE]
      colnames(weights) <- context$column_names
      rownames(weights) <- columns
    }
    if(is.null(dim(weights))) weights <- matrix(weights, nrow = 1L,
      dimnames = list(parameter, names(weights)))
    if(!is.matrix(x) && nrow(weights) > 1L){
      density <- .prior_density_from_context_rows(context, weights)
      if(!is.null(fields$linear_offset) && fields$linear_offset != 0){
        density <- .prior_density_output_transform(density,
          .bt_posterior_transformation("lin", list(a = fields$linear_offset, b = 1)))
      }
      return(.bt_meta_set(x, "prior_density", density))
    }
    if(is.matrix(x) && nrow(weights) != ncol(x)){
      stop("The declared source prior recipe does not align with its draw columns.", call. = FALSE)
    }
    densities <- lapply(seq_len(nrow(weights)), function(i){
      target <- if(nrow(weights) == 1L && !is.null(parameter)){
        .marginal_posterior_simple_target(context, parameter)
      }else list(source_transforms = NULL, output_transformation = NULL)
      density <- .prior_density_from_context(context, stats::setNames(weights[i, ], colnames(weights)),
        source_transforms = target$source_transforms, output_transformation = target$output_transformation)
      if(!is.null(fields$linear_offset) && fields$linear_offset != 0){
        density <- .prior_density_output_transform(density,
          .bt_posterior_transformation("lin", list(a = fields$linear_offset, b = 1)))
      }
      density
    })
    if(nrow(weights) == 1L && !is.matrix(x)){
      .bt_meta_set(x, "prior_density", densities[[1L]])
    }else{
      names(densities) <- if(is.null(rownames(weights))) colnames(x) else rownames(weights)
      .bt_meta_set(x, "prior_densities", densities)
    }
  }
  tryCatch(build(),
    BayesTools_formula_measure_unavailable = function(condition){
      if(!condition$reason %in% setdiff(.bt_formula_measure_causes, "missing_multiplier_law")) stop(condition)
      .bt_formula_measure_mark(x, if(is.null(parameter)) colnames(x) else parameter,
        "prior_density", conditionMessage(condition), cause = condition$reason, diagnostics = condition$diagnostics)
    },
    BayesTools_numerical_condition = function(condition){
      .bt_formula_measure_mark(x, if(is.null(parameter)) colnames(x) else parameter,
        "prior_density", conditionMessage(condition), cause = "numerical_scale_unavailable")
    })
}

.bt_posterior_transform_draws <- function(x, map, parent = NULL, parameter = NULL){

  if(!.bt_meta_is_draws(x) ||
     (is.null(.bt_meta_container(x)) &&
      !inherits(x, c("mixed_posteriors", "marginal_posterior")))){
    .bt_draws_stop_plain(
      "'posterior_transform' requires BayesTools posterior draws, not plain numeric draws"
    )
  }

  fields <- .bt_meta_get_fields(x, .bt_meta_fields())
  atoms <- .posterior_atoms_from_attribute(fields$atoms)
  .bt_posterior_transformation_check(
    map,
    c(as.numeric(x), if(!is.null(atoms)) as.numeric(atoms$locations))
  )

  if(is.null(parameter)) parameter <- attr(x, "parameter", exact = TRUE)
  if(is.null(parameter) && is.matrix(x)) parameter <- colnames(x)
  x <- .bt_posterior_transform_source_law(x, parent, parameter)
  fields <- .bt_meta_get_fields(x, .bt_meta_fields())

  out <- .bt_draws_transform_values(x, map$fun)
  updates <- .bt_posterior_transform_fields(fields, map,
    parameter = parameter)
  # the draws record the transformations applied to their values
  updates["output_transformations"] <- list(c(fields$output_transformations, map$name))
  out <- .bt_meta_assign(out, updates)
  out
}

# The transformed values of the metadata fields present in 'fields' (a named
# list of fields; NULL values are absent fields). NULL values in the result
# remove a field.
.bt_posterior_transform_fields <- function(fields, map, parameter = NULL){

  updates <- list()
  set <- function(field, value){
    updates[field] <<- list(value)
  }
  transform_prior <- function(density, column = parameter){
    mark <- function(condition, cause){
      if(is.null(column) && !is.null(fields$quantities)) column <- fields$quantities$column
      if(is.null(column)) stop(condition)
      holder <- .bt_meta_set(numeric(), "measure_unavailable",
        if(!is.null(updates$measure_unavailable)) updates$measure_unavailable else fields$measure_unavailable)
      for(target in column) holder <- .bt_formula_measure_mark(holder, target, "prior_density",
        conditionMessage(condition), cause = cause, diagnostics = condition$diagnostics)
      set("measure_unavailable", .bt_meta_get(holder, "measure_unavailable"))
      NULL
    }
    tryCatch(.prior_density_output_transform(density, map),
      BayesTools_formula_measure_unavailable = function(condition){
        if(!condition$reason %in% setdiff(.bt_formula_measure_causes, "missing_multiplier_law")) stop(condition)
        mark(condition, condition$reason)
      },
      BayesTools_numerical_condition = function(condition){
        mark(condition, "numerical_scale_unavailable")
      })
  }
  transform_support <- function(support){
    .posterior_support_transform(support, map$transformation, map$arguments)
  }
  if(!is.null(fields$hypothesis_evaluation)) set("hypothesis_evaluation", NULL)

  if(!is.null(fields$support)){
    set("support", if(.posterior_metadata_is_container(fields$support)){
      lapply(fields$support, transform_support)
    }else{
      transform_support(fields$support)
    })
  }
  if(!is.null(fields$atoms)){
    set("atoms", .posterior_atoms_transform(
      fields$atoms, map$transformation, map$arguments
    ))
  }
  if(!is.null(fields$ordered_source)){
    source <- fields$ordered_source
    source$view_transformations <- c(source$view_transformations,list(list(
      transformation=if(is.character(map$transformation)) map$transformation else "unavailable",
      arguments=map$arguments)))
    set("ordered_source",source)
  }
  if(!is.null(fields$components)){
    components <- fields$components
    set("components", .posterior_components_new(
      index    = components$index,
      supports = lapply(components$supports, transform_support),
      keys     = components$keys
    ))
  }
  if(!is.null(fields$prior_density)){
    set("prior_density", transform_prior(fields$prior_density))
  }
  if(!is.null(fields$prior_densities)){
    densities <- fields$prior_densities
    for(i in seq_along(densities)){
      if(!is.null(densities[[i]])){
        column <- if(is.null(names(densities))) parameter else names(densities)[[i]]
        densities[i] <- list(transform_prior(densities[[i]], column))
      }
    }
    set("prior_densities", densities)
  }
  for(field in c("posterior_density", "posterior_densities")){
    if(!is.null(fields[[field]])){
      set(field, .posterior_density_output_transform(fields[[field]], map))
    }
  }
  for(field in c("posterior_ordinate", "posterior_ordinates")){
    if(!is.null(fields[[field]])){
      set(field, .posterior_ordinate_output_transform(fields[[field]], map))
    }
  }
  if(!is.null(fields$quantities)){
    set("quantities", .bt_quantities_output_transform(fields$quantities, map))
  }
  if(!is.null(fields$original_scale_quantities)){
    # fitted-scale replacements do not apply to transformed values
    set("original_scale_quantities", NULL)
  }
  if(!is.null(fields$linear_weights)){
    if(identical(map$name, "lin")){
      a <- .posterior_support_transform_argument(map$arguments, "a", 0)
      b <- .posterior_support_transform_argument(map$arguments, "b", 1)
      offset <- if(is.null(fields$linear_offset)) 0 else fields$linear_offset
      set("linear_weights", b * fields$linear_weights)
      set("linear_offset", a + b * offset)
    }else{
      set("linear_weights", NULL)
      set("linear_offset", NULL)
      set("joint_prior_transformation", map$name)
    }
  }

  updates
}

# The column table of transformed draws: the transformation is appended to
# the label parts, the columns are no catalog quantities, and only 'lin' with
# 'a = 0' keeps them linear in the fitted coordinates.
.bt_quantities_output_transform <- function(quantities, map){

  parts <- lapply(unclass(quantities$label_parts), function(part){
    part$transformation <- c(part$transformation, map$name)
    .bt_validate_label_parts(part)
    part
  })
  n <- nrow(quantities)
  linear <- identical(map$name, "lin") &&
    .posterior_support_transform_argument(map$arguments, "a", 0) == 0
  b <- .posterior_support_transform_argument(map$arguments, "b", 1)

  .bt_draws_quantity_table(
    columns      = quantities$column,
    quantity_ids = rep("", n),
    dependencies = if(linear) unclass(quantities$dependencies) else rep(list(character()), n),
    weights      = if(linear){
      lapply(unclass(quantities$weights), function(weights) b * weights)
    }else{
      rep(list(numeric()), n)
    },
    label_parts  = parts
  )
}

# The prior density of transformed values (change of variables). Densities
# whose builder records no output transformation are rebuilt with it (the
# density a transformation passed to the builder gives); densities with an
# output transformation or another recorded provenance (allocation products)
# get the transformation on top of it ('output_transformation' provenance);
# a grid without provenance is transformed for plotting only.
.prior_density_output_transform <- function(density, map){

  if(is.prior(density)){
    return(.prior_linear_combination_density(
      prior_list                      = list(value = density),
      weights                         = c(value = 1),
      output_transformation           = map$transformation,
      output_transformation_arguments = map$arguments
    ))
  }

  adaptive <- attr(density, "adaptive_evaluation", exact = TRUE)
  builders <- list(
    linear_combination   = .prior_linear_combination_density,
    density_context      = .prior_density_from_context,
    density_context_rows = .prior_density_from_context_rows
  )
  if(is.list(adaptive) && is.character(adaptive$kind) &&
     length(adaptive$kind) == 1L && adaptive$kind %in% names(builders) &&
     is.null(adaptive$arguments$output_transformation)){
    arguments <- adaptive$arguments
    arguments$output_transformation           <- map$transformation
    arguments$output_transformation_arguments <- map$arguments
    return(do.call(builders[[adaptive$kind]], arguments))
  }

  out <- .prior_linear_density_transform(
    density,
    map$transformation,
    map$arguments,
    n_grid = density$n_grid
  )
  for(attribute in c("product_grid_resolution", "provenance_unavailable")){
    attr(out, attribute) <- attr(density, attribute, exact = TRUE)
  }
  if(is.list(adaptive)){
    attr(out, "adaptive_evaluation") <- list(
      kind      = "output_transformation",
      arguments = list(
        source                          = adaptive,
        output_transformation           = map$transformation,
        output_transformation_arguments = map$arguments
      )
    )
  }
  out
}

# Precomputed posterior densities of transformed values: locations mapped and
# heights divided by the absolute derivative; NULL when the transformed grid
# is not finite (a vanishing or infinite derivative).
.posterior_density_output_transform <- function(density, map){

  kind <- .posterior_density_kind(density)
  if(identical(kind, "null")){
    return(NULL)
  }
  if(identical(kind, "container")){
    out <- lapply(density, .posterior_density_output_transform, map = map)
    out <- out[!vapply(out, is.null, logical(1))]
    return(if(length(out) > 0L) out)
  }

  density <- .posterior_density_from_attribute(density)
  x <- suppressWarnings(map$fun(density$x))
  y <- density$y / abs(suppressWarnings(map$jac(density$x)))
  if(any(!is.finite(x)) || any(!is.finite(y))){
    return(NULL)
  }
  order_x <- order(x)
  density$x <- x[order_x]
  density$y <- y[order_x]
  density["support"] <- list(.posterior_support_transform(
    density$support, map$transformation, map$arguments
  ))
  density["diagnostics"] <- list(.bt_posterior_transform_diagnostics(
    density$diagnostics, map
  ))
  if(is.character(.posterior_density_normalize(density))){
    return(NULL)
  }

  density
}

# Precomputed posterior ordinates of transformed values: the values mapped
# and the ordinates divided by the absolute derivative there; an ordinate
# whose transformed values or heights are not finite and positive is removed.
.posterior_ordinate_output_transform <- function(ordinate, map){

  kind <- .posterior_ordinate_kind(ordinate)
  if(identical(kind, "null")){
    return(NULL)
  }
  if(identical(kind, "container")){
    out <- lapply(ordinate, .posterior_ordinate_output_transform, map = map)
    out <- out[!vapply(out, is.null, logical(1))]
    return(if(length(out) > 0L) out)
  }
  if(identical(kind, "ordinates")){
    entries <- lapply(
      .posterior_ordinate_entries(ordinate),
      .posterior_ordinate_output_transform,
      map = map
    )
    entries <- entries[!vapply(entries, is.null, logical(1))]
    out <- NULL
    for(entry in entries){
      out <- posterior_ordinate_append(out, entry)
    }
    return(out)
  }

  values <- .posterior_ordinate_values(ordinate)
  value <- suppressWarnings(map$fun(values$value))
  jacobian <- abs(suppressWarnings(map$jac(values$value)))
  log_height <- if(is.null(values$log_ordinate)) NULL else values$log_ordinate - log(jacobian)
  height <- if(is.null(log_height)) values$ordinate / jacobian else exp(log_height)
  if(any(!is.finite(value)) ||
     (is.null(log_height) && (any(!is.finite(height)) || any(height <= 0))) ||
     (!is.null(log_height) && any(!is.finite(log_height))) ||
     anyDuplicated(value)){
    return(NULL)
  }
  ordinate$value <- value
  ordinate$ordinate <- height
  if(!is.null(log_height)) ordinate$log_ordinate <- log_height
  ordinate["diagnostics"] <- list(.bt_posterior_transform_diagnostics(
    ordinate$diagnostics, map
  ))

  ordinate
}

# Estimator diagnostics of transformed densities and ordinates: the entries
# that identify the null value of a stored Bayes-factor evaluation are mapped
# with the values; the other diagnostics are kept as estimated.
.bt_posterior_transform_diagnostics <- function(diagnostics, map){

  if(is.null(diagnostics) || !is.list(diagnostics)){
    return(diagnostics)
  }
  for(name in intersect(c("bf_value", "value", "null_hypothesis"), names(diagnostics))){
    value <- diagnostics[[name]]
    if(is.numeric(value)){
      diagnostics[[name]] <- suppressWarnings(map$fun(value))
    }
  }

  diagnostics
}

# The transformations applied to the values of draws by
# posterior_transform() (their 'output_transformations' metadata; empty for
# untransformed draws).
.bt_draws_output_transformations <- function(x){

  transformations <- .bt_meta_get(x, "output_transformations")
  if(is.null(transformations)) character() else transformations
}
