# ============================================================================ #
# posterior-support.R
# ============================================================================ #
#
# Internal helpers for exact scalar support metadata on posterior samples.
#
# ============================================================================ #

.posterior_support_new <- function(bounds, points = numeric(), exact = TRUE,
                                   source = NULL, type = NULL){

  bounds <- .posterior_support_parse_bounds(bounds)
  if(is.null(bounds)){
    stop("'bounds' must be an ordered numeric vector of length 2.", call. = FALSE)
  }
  points <- .posterior_support_parse_points(points)
  if(is.null(points)){
    stop("'points' must be a numeric vector.", call. = FALSE)
  }
  type <- .posterior_support_parse_type(type, bounds, points)
  if(is.null(type)){
    stop("'type' must be one of 'interval', 'points', or 'mixed'.",
         call. = FALSE)
  }

  out <- list(
    bounds = bounds,
    points = points,
    exact  = isTRUE(exact),
    source = source,
    type   = type
  )
  class(out) <- c("BayesTools_posterior_support", "list")

  out
}

.posterior_support_parse_bounds <- function(bounds){

  if(!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds)){
    return(NULL)
  }

  bounds <- as.numeric(bounds)
  if(bounds[1] > bounds[2]){
    return(NULL)
  }

  bounds
}

.posterior_support_parse_points <- function(points){

  if(is.null(points)){
    return(numeric())
  }
  if(!is.numeric(points)){
    return(NULL)
  }

  points <- as.numeric(points)
  points <- points[is.finite(points)]
  unique(points)
}

.posterior_support_parse_type <- function(type, bounds, points){

  if(is.null(type)){
    if(length(points) > 0L &&
       all(is.finite(bounds)) &&
       isTRUE(all.equal(range(points), bounds))){
      return("points")
    }
    return("interval")
  }

  if(!is.character(type) || length(type) != 1L || is.na(type)){
    return(NULL)
  }
  type <- tryCatch(
    match.arg(type, c("interval", "points", "mixed")),
    error = function(e) NULL
  )
  if(is.null(type)){
    return(NULL)
  }

  if(identical(type, "points") && length(points) == 0L){
    return(NULL)
  }

  type
}

.posterior_support_has_interval <- function(support){

  support <- .posterior_support_from_attribute(support)
  if(is.null(support)){
    return(FALSE)
  }

  is.null(support$type) ||
    support$type %in% c("interval", "mixed")
}

.posterior_support_from_attribute <- function(support, exact = NULL){

  if(is.null(support)){
    return(NULL)
  }

  if(inherits(support, "BayesTools_posterior_support")){
    if(!is.null(exact)){
      support$exact <- isTRUE(exact)
    }
    support$type <- .posterior_support_parse_type(
      support$type,
      support$bounds,
      support$points
    )
    if(is.null(support$type)){
      return(NULL)
    }
    return(support)
  }

  support_exact <- if(is.null(exact)) TRUE else isTRUE(exact)
  support_source <- NULL
  support_points <- numeric()
  support_type <- NULL

  if(is.numeric(support) && length(support) == 2L){
    bounds <- support
  }else if(is.list(support)){
    if(!is.null(support[["exact"]])){
      support_exact <- isTRUE(support[["exact"]])
    }
    if(!is.null(support[["source"]])){
      support_source <- support[["source"]]
    }
    if(!is.null(support[["points"]])){
      support_points <- support[["points"]]
    }
    if(!is.null(support[["type"]])){
      support_type <- support[["type"]]
    }else if(!is.null(support[["support_type"]])){
      support_type <- support[["support_type"]]
    }

    if(!is.null(support[["bounds"]])){
      bounds <- support[["bounds"]]
    }else if(!is.null(support[["lower"]]) && !is.null(support[["upper"]])){
      bounds <- c(support[["lower"]], support[["upper"]])
    }else{
      return(NULL)
    }
  }else{
    return(NULL)
  }

  bounds <- .posterior_support_parse_bounds(bounds)
  support_points <- .posterior_support_parse_points(support_points)
  if(is.null(bounds) || is.null(support_points)){
    return(NULL)
  }

  .posterior_support_new(
    bounds = bounds,
    points = support_points,
    exact  = support_exact,
    source = support_source,
    type   = support_type
  )
}

.posterior_support_source_is_raw_prior <- function(source){

  if(is.null(source)){
    return(FALSE)
  }
  source <- as.character(source)
  source <- source[!is.na(source) & nzchar(source)]
  if(length(source) == 0L){
    return(FALSE)
  }

  all(source %in% c("prior", "prior_list", "weightfunction"))
}

.posterior_support_is_raw_prior <- function(support){

  direct_support <- .posterior_support_from_attribute(support)
  if(!is.null(direct_support)){
    return(.posterior_support_source_is_raw_prior(direct_support$source))
  }
  if(!is.list(support)){
    return(FALSE)
  }

  support_list <- lapply(support, .posterior_support_from_attribute)
  support_list <- support_list[!vapply(support_list, is.null, logical(1))]
  if(length(support_list) == 0L){
    return(FALSE)
  }

  all(vapply(
    support_list,
    function(x) .posterior_support_source_is_raw_prior(x$source),
    logical(1)
  ))
}

.posterior_support_set <- function(x, support){

  if(is.null(support)){
    attr(x, "posterior_support") <- NULL
    return(x)
  }

  attr(x, "posterior_support") <- support
  x
}

.posterior_support_drop <- function(x, recursive = FALSE){

  attr(x, "posterior_support") <- NULL

  if(isTRUE(recursive) && is.list(x)){
    x_attributes <- attributes(x)
    for(i in seq_along(x)){
      x[[i]] <- .posterior_support_drop(x[[i]], recursive = TRUE)
    }
    attributes(x) <- x_attributes
    attr(x, "posterior_support") <- NULL
  }

  x
}

.posterior_support_get <- function(x, name = NULL){

  support <- attr(x, "posterior_support", exact = TRUE)
  if(is.null(support)){
    return(NULL)
  }

  if(!is.null(name) && !inherits(support, "BayesTools_posterior_support")){
    if(!is.null(support[[name]])){
      return(.posterior_support_from_attribute(support[[name]]))
    }
    return(NULL)
  }

  .posterior_support_from_attribute(support)
}

.posterior_support_bounds <- function(x, name = NULL, exact_only = TRUE,
                                      interval_only = FALSE){

  support <- .posterior_support_get(x, name = name)
  if(is.null(support)){
    return(NULL)
  }
  if(isTRUE(exact_only) && !isTRUE(support$exact)){
    return(NULL)
  }
  if(isTRUE(interval_only) && !.posterior_support_has_interval(support)){
    return(NULL)
  }

  support$bounds
}

.posterior_support_bounds_contains_value <- function(bounds, value,
                                                     tolerance = sqrt(.Machine$double.eps)){

  bounds <- .posterior_support_parse_bounds(bounds)
  if(is.null(bounds) || length(value) != 1L || !is.finite(value)){
    return(NA)
  }

  scale <- max(1, abs(value), abs(bounds[is.finite(bounds)]))
  tol <- tolerance * scale

  lower_contains <- !is.finite(bounds[1]) || value >= bounds[1] - tol
  upper_contains <- !is.finite(bounds[2]) || value <= bounds[2] + tol

  lower_contains && upper_contains
}

.posterior_support_contains_value <- function(support, value,
                                              exact_only = TRUE,
                                              tolerance = sqrt(.Machine$double.eps)){

  support <- .posterior_support_from_attribute(support)
  if(is.null(support)){
    return(NA)
  }
  if(isTRUE(exact_only) && !isTRUE(support$exact)){
    return(NA)
  }

  if(!.posterior_support_has_interval(support)){
    points <- support$points
    if(length(points) == 0L || length(value) != 1L || !is.finite(value)){
      return(NA)
    }
    scale <- max(1, abs(value), abs(points))
    tol <- tolerance * scale
    return(any(abs(points - value) <= tol))
  }

  .posterior_support_bounds_contains_value(
    support$bounds,
    value     = value,
    tolerance = tolerance
  )
}

.posterior_support_excludes_value <- function(support, value,
                                              exact_only = TRUE,
                                              tolerance = sqrt(.Machine$double.eps)){

  identical(
    .posterior_support_contains_value(
      support,
      value      = value,
      exact_only = exact_only,
      tolerance  = tolerance
    ),
    FALSE
  )
}

.posterior_support_bounds_connected <- function(bounds){

  if(is.null(bounds) || nrow(bounds) <= 1L){
    return(TRUE)
  }

  bounds <- bounds[order(bounds[, 1], bounds[, 2]), , drop = FALSE]
  current_upper <- bounds[1, 2]

  for(i in 2:nrow(bounds)){
    if(is.finite(current_upper) && is.finite(bounds[i, 1])){
      scale <- max(1, abs(current_upper), abs(bounds[i, 1]))
      tol <- sqrt(.Machine$double.eps) * scale
      if(bounds[i, 1] > current_upper + tol){
        return(FALSE)
      }
    }else if(is.infinite(current_upper) && current_upper < 0 &&
             !(is.infinite(bounds[i, 1]) && bounds[i, 1] < 0)){
      return(FALSE)
    }
    current_upper <- max(current_upper, bounds[i, 2])
  }

  TRUE
}

.posterior_support_union <- function(supports, source = NULL){

  supports <- lapply(supports, .posterior_support_from_attribute)
  supports <- supports[!vapply(supports, is.null, logical(1))]
  if(length(supports) == 0L){
    return(NULL)
  }

  bounds <- do.call(rbind, lapply(supports, function(support) support$bounds))
  points <- unique(unlist(lapply(supports, function(support) support$points),
                          use.names = FALSE))
  has_interval <- any(vapply(supports, .posterior_support_has_interval, logical(1)))
  has_points <- length(points) > 0L
  type <- if(has_interval && has_points){
    "mixed"
  }else if(has_interval){
    "interval"
  }else{
    "points"
  }
  exact <- all(vapply(supports, function(support) isTRUE(support$exact), logical(1)))
  if(has_interval){
    exact <- exact && .posterior_support_bounds_connected(bounds)
  }

  .posterior_support_new(
    bounds = c(min(bounds[, 1], na.rm = TRUE), max(bounds[, 2], na.rm = TRUE)),
    points = points,
    exact  = exact,
    source = source,
    type   = type
  )
}

.posterior_support_sum <- function(supports, source = NULL){

  supports <- lapply(supports, .posterior_support_from_attribute)
  supports <- supports[!vapply(supports, is.null, logical(1))]
  if(length(supports) == 0L){
    return(NULL)
  }

  bounds <- do.call(rbind, lapply(supports, function(support) support$bounds))
  points <- numeric()
  point_lists <- lapply(supports, function(support) support$points)
  if(all(vapply(point_lists, length, integer(1)) > 0L)){
    point_grid <- expand.grid(point_lists)
    points <- rowSums(point_grid)
  }
  has_interval <- any(vapply(supports, .posterior_support_has_interval, logical(1)))
  has_points <- length(points) > 0L
  type <- if(has_interval && has_points){
    "mixed"
  }else if(has_interval){
    "interval"
  }else{
    "points"
  }
  exact <- all(vapply(supports, function(support) isTRUE(support$exact), logical(1)))
  if(has_interval && any(vapply(
    supports,
    function(support) !.posterior_support_has_interval(support) &&
      length(support$points) > 1L,
    logical(1)
  ))){
    exact <- FALSE
  }

  .posterior_support_new(
    bounds = c(sum(bounds[, 1]), sum(bounds[, 2])),
    points = points,
    exact  = exact,
    source = source,
    type   = type
  )
}

.posterior_support_scale <- function(support, weight, source = NULL){

  support <- .posterior_support_from_attribute(support)
  if(is.null(support)){
    return(NULL)
  }

  if(abs(weight) <= .prior_linear_density_zero_tol()){
    return(.posterior_support_new(c(0, 0), points = 0, source = source,
                                  type = "points"))
  }

  bounds <- weight * support$bounds
  points <- weight * support$points

  .posterior_support_new(
    bounds = range(bounds),
    points = points,
    exact  = support$exact,
    source = source,
    type   = support$type
  )
}

.posterior_support_transform <- function(support, transformation,
                                         transformation_arguments = NULL){

  support <- .posterior_support_from_attribute(support)
  if(is.null(support) || is.null(transformation)){
    return(support)
  }

  if(!is.character(transformation) || length(transformation) != 1L){
    return(NULL)
  }
  if(!transformation %in% c("lin", "exp_lin", "tanh", "exp")){
    return(NULL)
  }
  if(!.posterior_support_transform_is_valid(transformation, transformation_arguments)){
    return(NULL)
  }

  transformed_bounds <- suppressWarnings(.density.prior_transformation_x(
    support$bounds,
    transformation,
    transformation_arguments
  ))
  if(length(transformed_bounds) != 2L || anyNA(transformed_bounds)){
    return(NULL)
  }

  transformed_points <- numeric()
  if(length(support$points) > 0L){
    transformed_points <- suppressWarnings(.density.prior_transformation_x(
      support$points,
      transformation,
      transformation_arguments
    ))
    transformed_points <- transformed_points[is.finite(transformed_points)]
  }

  .posterior_support_new(
    bounds = range(transformed_bounds),
    points = transformed_points,
    exact  = support$exact,
    source = support$source,
    type   = support$type
  )
}

.posterior_support_transform_is_valid <- function(transformation,
                                                  transformation_arguments = NULL){

  if(transformation %in% c("lin", "exp_lin")){
    b <- .posterior_support_transform_argument(
      transformation_arguments,
      name    = "b",
      default = 1
    )
    if(!is.finite(b) || abs(b) <= .prior_linear_density_zero_tol()){
      return(FALSE)
    }
  }

  TRUE
}

.posterior_support_transform_argument <- function(transformation_arguments,
                                                  name,
                                                  default){

  if(is.null(transformation_arguments) ||
     is.null(transformation_arguments[[name]])){
    return(default)
  }

  value <- transformation_arguments[[name]]
  if(!is.numeric(value) || length(value) != 1L || anyNA(value)){
    return(NA_real_)
  }

  as.numeric(value)
}

.posterior_support_point <- function(value, source = NULL){

  .posterior_support_new(
    c(value, value),
    points = value,
    source = source,
    type   = "points"
  )
}

.posterior_support_positive_prior_indices <- function(priors){

  if(is.null(priors) || !is.list(priors)){
    return(NULL)
  }

  weights <- NULL
  if(is.prior.mixture(priors)){
    weights <- attr(priors, "prior_weights", exact = TRUE)
  }else if(!is.prior(priors)){
    weights <- vapply(priors, function(prior){
      if(!is.prior(prior)){
        return(NA_real_)
      }
      weight <- .prior_model_weight(prior)
      if(!is.numeric(weight) || length(weight) != 1L){
        return(NA_real_)
      }
      as.numeric(weight)
    }, numeric(1))
  }

  if(is.null(weights) || length(weights) != length(priors)){
    return(seq_along(priors))
  }

  which(is.finite(weights) & weights > 0)
}

.posterior_support_from_prior <- function(prior, source = "prior"){

  if(is.null(prior)){
    return(NULL)
  }
  if(is.prior.none(prior)){
    return(.posterior_support_new(c(0, 0), points = 0, source = source,
                                  type = "points"))
  }
  if(is.prior.spike_and_slab(prior)){
    variable_support <- .posterior_support_from_prior(
      .get_spike_and_slab_variable(prior),
      source = source
    )
    return(.posterior_support_union(list(
      .posterior_support_new(c(0, 0), points = 0, source = source,
                             type = "points"),
      variable_support
    ), source = source))
  }
  if(is.prior.mixture(prior)){
    positive_indices <- .posterior_support_positive_prior_indices(prior)
    return(.posterior_support_union(
      lapply(prior[positive_indices], .posterior_support_from_prior, source = source),
      source = source
    ))
  }
  if(is.prior.point(prior) && !is.prior.vector(prior)){
    location <- prior$parameters[["location"]]
    return(.posterior_support_new(
      c(location, location),
      points = location,
      source = source,
      type   = "points"
    ))
  }
  if(is.prior.discrete(prior)){
    support <- switch(
      prior[["distribution"]],
      "bernoulli" = c(0, 1),
      numeric()
    )
    support <- support[
      support >= prior$truncation[["lower"]] &
        support <= prior$truncation[["upper"]]
    ]
    if(length(support) == 0L){
      return(NULL)
    }
    return(.posterior_support_new(
      range(support),
      points = support,
      source = source,
      type   = "points"
    ))
  }
  if(is.prior.weightfunction(prior)){
    return(.posterior_support_union(
      lapply(
        .weightfunction_marginal_components(prior),
        function(component) .posterior_support_from_weightfunction_component(component, source)
      ),
      source = source
    ))
  }
  if(is_prior_phacking(prior) || is_prior_bias(prior)){
    return(NULL)
  }
  if(is.prior.simple(prior)){
    bounds <- c(prior$truncation[["lower"]], prior$truncation[["upper"]])
    return(.posterior_support_new(bounds, source = source, type = "interval"))
  }
  if(is.prior.simplex(prior)){
    return(.posterior_support_new(c(0, 1), source = source,
                                  type = "interval"))
  }
  if(is.prior.vector(prior)){
    if(is.prior.point(prior)){
      location <- prior$parameters[["location"]]
      if(!is.numeric(location) || length(location) != 1L){
        return(NULL)
      }
      return(.posterior_support_point(location, source = source))
    }
    return(.posterior_support_new(c(-Inf, Inf), source = source,
                                  type = "interval"))
  }

  NULL
}

.posterior_support_from_prior_list <- function(priors, source = "prior_list"){

  if(is.null(priors)){
    return(NULL)
  }
  if(is.prior(priors)){
    return(.posterior_support_from_prior(priors, source = source))
  }
  if(!is.list(priors)){
    return(NULL)
  }

  positive_indices <- .posterior_support_positive_prior_indices(priors)
  .posterior_support_union(
    lapply(priors[positive_indices], .posterior_support_from_prior, source = source),
    source = source
  )
}

.posterior_support_set_from_prior_list <- function(samples, priors,
                                                   source = "prior_list"){

  support <- .posterior_support_from_prior_list(priors, source = source)
  .posterior_support_set(samples, support)
}

.posterior_support_set_columns_from_prior_list <- function(samples, priors,
                                                           column_names = colnames(samples),
                                                           source = "prior_list"){

  support <- .posterior_support_from_prior_list(priors, source = source)
  if(is.null(support)){
    return(.posterior_support_set(samples, NULL))
  }
  if(is.null(column_names)){
    return(.posterior_support_set(samples, support))
  }

  attr(samples, "posterior_support") <- stats::setNames(
    rep(list(support), length(column_names)),
    column_names
  )
  samples
}

.posterior_support_from_weightfunction_component <- function(component,
                                                             source = "weightfunction"){

  bounds <- .density.prior_weightfunction_component_bounds(component)
  points <- if(identical(component$type, "point")) component$location else numeric()

  .posterior_support_new(
    bounds = bounds,
    points = points,
    source = source,
    type   = if(identical(component$type, "point")) "points" else "interval"
  )
}

.posterior_support_weightfunction_columns <- function(priors, omega_context,
                                                      source = "weightfunction"){

  if(is.null(omega_context)){
    return(NULL)
  }
  if(is.prior(priors)){
    priors <- list(priors)
  }
  if(!is.list(priors)){
    return(NULL)
  }

  omega_names <- omega_context$names
  out <- vector("list", length(omega_names))
  names(out) <- omega_names

  for(bin_i in seq_along(omega_names)){
    bin_supports <- list()
    for(prior_i in .posterior_support_positive_prior_indices(priors)){
      prior <- priors[[prior_i]]
      if(.is_prior_weightfunction_null(prior)){
        bin_supports[[length(bin_supports) + 1L]] <-
          .posterior_support_new(c(1, 1), points = 1, source = source,
                                 type = "points")
      }else if(is.prior.weightfunction(prior)){
        components <- .weightfunction_marginal_components(prior)
        component_i <- omega_context$mapping[[prior_i]][bin_i]
        if(length(component_i) == 1L && !is.na(component_i)){
          bin_supports[[length(bin_supports) + 1L]] <-
            .posterior_support_from_weightfunction_component(
              components[[component_i]],
              source = source
            )
        }
      }
    }
    out[[bin_i]] <- .posterior_support_union(bin_supports, source = source)
  }

  out
}

.posterior_support_set_weightfunction_columns <- function(samples, priors,
                                                          omega_context,
                                                          source = "weightfunction"){

  support <- .posterior_support_weightfunction_columns(
    priors        = priors,
    omega_context = omega_context,
    source        = source
  )
  if(is.null(support)){
    return(samples)
  }

  .posterior_support_set(samples, support)
}

.posterior_support_group_full_weights <- function(group){

  prior <- group$prior
  K <- prior$parameters[["K"]]
  if(!is.numeric(K) || length(K) != 1L || is.na(K) || K < 1L){
    return(NULL)
  }
  K <- as.integer(K)

  indices <- group$indices
  weights <- group$weights
  if(length(indices) != length(weights) || anyNA(indices) ||
     any(indices < 1L | indices > K)){
    return(NULL)
  }

  full_weights <- numeric(K)
  full_weights[indices] <- as.numeric(weights)
  full_weights
}

.posterior_support_group_simplex <- function(group, source = "linear_prior"){

  full_weights <- .posterior_support_group_full_weights(group)
  if(is.null(full_weights) || any(!is.finite(full_weights))){
    return(NULL)
  }

  bounds <- range(full_weights)
  if(isTRUE(all.equal(bounds[1], bounds[2]))){
    return(.posterior_support_point(bounds[1], source = source))
  }

  .posterior_support_new(bounds, source = source, type = "interval")
}

.posterior_support_group_vector_point <- function(group, source = "linear_prior"){

  location <- group$prior$parameters[["location"]]
  if(!is.numeric(location) || anyNA(location)){
    return(NULL)
  }

  weights <- as.numeric(group$weights)
  indices <- group$indices
  if(length(weights) != length(indices) || anyNA(indices)){
    return(NULL)
  }

  if(length(location) == 1L){
    value <- sum(weights) * location
  }else{
    if(any(indices < 1L | indices > length(location))){
      return(NULL)
    }
    value <- sum(weights * location[indices])
  }

  .posterior_support_point(value, source = source)
}

.posterior_support_group_linear <- function(group, source = "linear_prior"){

  prior <- group$prior
  weights <- group$weights

  if(is.prior.none(prior)){
    return(.posterior_support_new(c(0, 0), points = 0, source = source,
                                  type = "points"))
  }
  if(is.prior.spike_and_slab(prior)){
    variable_group <- group
    variable_group$prior <- .get_spike_and_slab_variable(prior)
    return(.posterior_support_union(list(
      .posterior_support_new(c(0, 0), points = 0, source = source,
                             type = "points"),
      .posterior_support_group_linear(variable_group, source = source)
    ), source = source))
  }
  if(is.prior.mixture(prior)){
    positive_indices <- .posterior_support_positive_prior_indices(prior)
    return(.posterior_support_union(lapply(prior[positive_indices], function(component){
      component_group <- group
      component_group$prior <- component
      .posterior_support_group_linear(component_group, source = source)
    }), source = source))
  }
  if(is.prior.vector(prior) && !is.prior.treatment(prior) && !is.prior.independent(prior)){
    if(is.prior.simplex(prior)){
      return(.posterior_support_group_simplex(group, source = source))
    }
    if(is.prior.point(prior)){
      return(.posterior_support_group_vector_point(group, source = source))
    }
    return(.posterior_support_new(c(-Inf, Inf), source = source,
                                  type = "interval"))
  }

  .posterior_support_sum(lapply(weights, function(weight){
    .posterior_support_scale(
      .posterior_support_from_prior(prior, source = source),
      weight,
      source = source
    )
  }), source = source)
}

.posterior_support_from_prior_list_weights <- function(prior_list, weights,
                                                       source = "linear_prior"){

  if(is.null(names(weights))){
    return(NULL)
  }
  weights <- weights[is.finite(weights)]
  weights <- weights[abs(weights) > .prior_linear_density_zero_tol()]
  if(length(weights) == 0L){
    return(.posterior_support_new(c(0, 0), points = 0, source = source,
                                  type = "points"))
  }

  groups <- .prior_linear_weight_groups(prior_list, weights)
  .posterior_support_sum(
    lapply(groups, .posterior_support_group_linear, source = source),
    source = source
  )
}

.posterior_support_from_prior_context_weights <- function(context, weights,
                                                          output_transformation = NULL,
                                                          output_transformation_arguments = NULL){

  if(is.null(context)){
    return(NULL)
  }

  if(!is.null(dim(weights))){
    supports <- lapply(seq_len(nrow(weights)), function(row_i){
      .posterior_support_from_prior_context_weights(
        context                         = context,
        weights                         = weights[row_i, ],
        output_transformation           = output_transformation,
        output_transformation_arguments = output_transformation_arguments
      )
    })
    return(.posterior_support_union(supports, source = "linear_prior_rows"))
  }

  if(inherits(context, "prior_density_context")){
    support <- .posterior_support_from_prior_list_weights(
      context$prior_list,
      .prior_density_context_standardized_weights(context, weights)
    )
  }else if(inherits(context, "prior_density_model_mixture_context")){
    model_indices <- which(is.finite(context$model_weights) & context$model_weights > 0)
    supports <- lapply(model_indices, function(model_i){
      model_prior_list <- lapply(context$prior_list, function(parameter_priors){
        if(is.prior(parameter_priors)){
          return(parameter_priors)
        }
        parameter_priors[[model_i]]
      })
      names(model_prior_list) <- names(context$prior_list)
      for(parameter in names(model_prior_list)){
        if(is.null(model_prior_list[[parameter]])){
          model_prior_list[[parameter]] <- prior("point", list(location = 0))
        }
      }
      .posterior_support_from_prior_list_weights(model_prior_list, weights)
    })
    support <- .posterior_support_union(supports, source = "linear_prior_models")
  }else if(inherits(context, "prior_density_conditional_context")){
    model_indices <- which(is.finite(context$model_weights) & context$model_weights > 0)
    supports <- lapply(model_indices, function(model_i){
      prior_list <- context$prior_lists[[model_i]]
      if(!is.null(context$formula_scale) && length(context$formula_scale) > 0L){
        component_context <- .prior_density_context(
          prior_list    = prior_list,
          column_names  = context$column_names,
          formula_scale = context$formula_scale,
          n_grid        = context$n_grid,
          tail_prob     = context$tail_prob
        )
        return(.posterior_support_from_prior_context_weights(component_context, weights))
      }
      .posterior_support_from_prior_list_weights(prior_list, weights)
    })
    support <- .posterior_support_union(supports, source = "linear_prior_conditioned")
  }else{
    return(NULL)
  }

  .posterior_support_transform(
    support,
    transformation          = output_transformation,
    transformation_arguments = output_transformation_arguments
  )
}

.posterior_support_set_from_prior_context <- function(samples, context,
                                                      parameter = attr(samples, "parameter", exact = TRUE)){

  if(is.null(context)){
    return(samples)
  }
  column_names <- context$column_names
  if(is.null(column_names)){
    return(samples)
  }

  if(!is.null(dim(samples))){
    sample_columns <- colnames(samples)
    if(is.null(sample_columns)){
      return(samples)
    }
    matched_columns <- sample_columns[sample_columns %in% column_names]
    if(length(matched_columns) == 0L){
      return(samples)
    }

    support <- vector("list", length(matched_columns))
    names(support) <- matched_columns
    for(column in matched_columns){
      weights <- rep(0, length(column_names))
      names(weights) <- column_names
      weights[[column]] <- 1
      support[[column]] <- tryCatch(
        .posterior_support_from_prior_context_weights(context, weights),
        error = function(e) NULL
      )
    }
    support <- support[!vapply(support, is.null, logical(1))]
    if(length(support) == 0L){
      return(samples)
    }
    existing_support <- attr(samples, "posterior_support", exact = TRUE)
    if(is.list(existing_support) &&
       !inherits(existing_support, "BayesTools_posterior_support") &&
       !is.null(names(existing_support))){
      existing_support <- existing_support[
        !names(existing_support) %in% names(support)
      ]
      existing_support <- existing_support[
        !vapply(existing_support, .posterior_support_is_raw_prior, logical(1))
      ]
      support <- c(existing_support, support)
    }
    attr(samples, "posterior_support") <- support
    return(samples)
  }

  if(is.null(parameter) || !parameter %in% column_names){
    return(samples)
  }

  weights <- rep(0, length(column_names))
  names(weights) <- column_names
  weights[[parameter]] <- 1
  support <- tryCatch(
    .posterior_support_from_prior_context_weights(context, weights),
    error = function(e) NULL
  )
  if(is.null(support)){
    return(samples)
  }
  .posterior_support_set(
    samples,
    support
  )
}

.posterior_support_for_kde <- function(samples, support = NULL,
                                       tolerance = sqrt(.Machine$double.eps)){

  if(is.null(support)){
    support <- .posterior_support_get(samples)
  }else{
    support <- .posterior_support_from_attribute(support)
  }
  if(is.null(support) || !isTRUE(support$exact)){
    return(list(bounds = NULL, warning = NULL))
  }
  if(!.posterior_support_has_interval(support)){
    return(list(
      bounds = NULL,
      warning = "Exact posterior support metadata does not describe a continuous interval. Falling back to the standard kernel density estimate."
    ))
  }

  bounds <- support$bounds
  sample_values <- as.numeric(samples)
  sample_values <- sample_values[is.finite(sample_values)]
  if(length(sample_values) > 0L){
    scale <- max(1, abs(bounds[is.finite(bounds)]), abs(sample_values))
    tol <- tolerance * scale
    outside <- rep(FALSE, length(sample_values))
    if(is.finite(bounds[1])){
      outside <- outside | sample_values < bounds[1] - tol
    }
    if(is.finite(bounds[2])){
      outside <- outside | sample_values > bounds[2] + tol
    }
    if(any(outside)){
      return(list(
        bounds = NULL,
        warning = "Exact posterior support metadata is incompatible with the posterior samples. Falling back to the standard kernel density estimate."
      ))
    }
  }

  list(bounds = bounds, warning = NULL)
}
