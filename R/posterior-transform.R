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
#'
#' The metadata are transformed as follows:
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
      x[[i]] <- .bt_posterior_transform_draws(x[[i]], map)
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

.bt_posterior_transform_draws <- function(x, map){

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

  out <- .bt_draws_transform_values(x, map$fun)
  updates <- .bt_posterior_transform_fields(fields, map)
  # the draws record the transformations applied to their values
  updates["output_transformations"] <- list(c(fields$output_transformations, map$name))
  out <- .bt_meta_assign(out, updates)
  out
}

# The transformed values of the metadata fields present in 'fields' (a named
# list of fields; NULL values are absent fields). NULL values in the result
# remove a field.
.bt_posterior_transform_fields <- function(fields, map){

  updates <- list()
  set <- function(field, value){
    updates[field] <<- list(value)
  }
  transform_support <- function(support){
    .posterior_support_transform(support, map$transformation, map$arguments)
  }

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
  if(!is.null(fields$components)){
    components <- fields$components
    set("components", .posterior_components_new(
      index    = components$index,
      supports = lapply(components$supports, transform_support),
      keys     = components$keys
    ))
  }
  if(!is.null(fields$prior_density)){
    set("prior_density", .prior_density_output_transform(fields$prior_density, map))
  }
  if(!is.null(fields$prior_densities)){
    densities <- fields$prior_densities
    for(i in seq_along(densities)){
      if(!is.null(densities[[i]])){
        densities[i] <- list(.prior_density_output_transform(densities[[i]], map))
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
# or point masses are not finite (a vanishing or infinite derivative).
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
  point_masses <- density$point_masses
  point_masses$x <- suppressWarnings(map$fun(point_masses$x))
  if(any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(point_masses$x))){
    return(NULL)
  }
  order_x <- order(x)
  density$x <- x[order_x]
  density$y <- y[order_x]
  density$point_masses <- point_masses[order(point_masses$x), , drop = FALSE]
  rownames(density$point_masses) <- NULL
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
  height <- values$ordinate / abs(suppressWarnings(map$jac(values$value)))
  if(any(!is.finite(value)) || any(!is.finite(height)) || any(height <= 0) ||
     anyDuplicated(value)){
    return(NULL)
  }
  ordinate$value <- value
  ordinate$ordinate <- height
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
