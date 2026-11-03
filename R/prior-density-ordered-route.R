# Routes of ordered terms that are not one allocation share of one slice
# (R/prior-density-ordered-product.R holds the leaf they build).

# The positive-probability components of an ordered total (recursively through
# mixtures and spike-and-slab priors): their probabilities and priors; NULL
# when the weights are not exact numeric values.
.prior_ordered_total_components <- function(total){

  if(!is.prior.mixture(total) && !is.prior.spike_and_slab(total)){
    return(list(probabilities = 1, priors = list(total)))
  }
  probabilities <- .prior_density_ordinate_mixture_weights(total)
  if(is.null(probabilities)){
    return(NULL)
  }
  out <- list(probabilities = numeric(), priors = list())
  for(i in which(probabilities > 0)){
    component <- .prior_ordered_total_components(total[[i]])
    if(is.null(component)){
      return(NULL)
    }
    out$probabilities <- c(out$probabilities, probabilities[[i]] * component$probabilities)
    out$priors <- c(out$priors, component$priors)
  }
  out
}

.prior_ordered_product_normal <- function(prior){

  is.prior.simple(prior) && identical(prior$distribution, "normal") &&
    .prior_density_ordinate_parameters_numeric(prior) &&
    identical(prior$truncation$lower, -Inf) && identical(prior$truncation$upper, Inf) &&
    length(prior$parameters$mean) == 1L && length(prior$parameters$sd) == 1L &&
    is.finite(prior$parameters$mean) && is.finite(prior$parameters$sd) && prior$parameters$sd > 0
}

# Route of a linear combination that involves ordered terms whose levels are
# not one allocation share of one slice (several slices, several allocations,
# shared allocations, other coefficients): the ordered-product leaf, or a
# refusal ('unknown' route with its reason). Mixture and spike-and-slab totals
# are expanded into one route per component combination (a shared indicator is
# one component per term, never one per slice); the terms that are not ordered
# must be points or normal terms.
.prior_density_route_ordered_product <- function(prior_list, weights, source_transforms, n_grid){

  refuse <- function(reason, cause = NULL){
    provenance <- list(kind = "unsupported_provenance", ordered = TRUE)
    if(!is.null(cause)){
      provenance$cause <- cause
    }
    .prior_density_route_unknown(reason = reason, provenance = provenance)
  }
  # priors with expression parameters have no R-side law
  expression_sources <- .prior_linear_expression_sources(prior_list, weights)
  if(length(expression_sources) > 0L){
    return(.prior_density_route_unknown(
      reason     = .prior_linear_expression_reason(expression_sources),
      provenance = list(kind = "expression_prior_law_unavailable", sources = expression_sources)
    ))
  }
  groups <- .prior_linear_weight_groups(prior_list, weights)
  is_ordered <- vapply(groups, function(group) is.prior.ordered(group$prior), logical(1))
  ordered <- groups[is_ordered]
  others <- groups[!is_ordered]
  if(any(vapply(groups, function(group) !is.null(attr(group$prior, "multiply_by", exact = TRUE)), logical(1)))){
    return(refuse("Scaled ordered-prior terms are not structurally classified."))
  }
  ordered_columns <- unlist(lapply(ordered, function(group) names(group$weights)), use.names = FALSE)
  if(any(!is.na(source_transforms[ordered_columns]))){
    return(refuse("The source transformation of an ordered term is not supported structurally."))
  }

  # the terms added to the ordered terms: mixtures first, then points and normal terms
  mixtures <- names(others)[vapply(others, function(group){
    is.prior.mixture(group$prior) || is.prior.spike_and_slab(group$prior)
  }, logical(1))]
  if(length(mixtures) > 0L){
    expansion <- .prior_density_route_mixture_expansion(prior_list, mixtures, n_grid)
    if(is.null(expansion)){
      return(refuse("The mixture expansion of the terms added to the ordered terms exceeds the leaf cap.",
                    "leaf_cap"))
    }
    return(.prior_density_route_expand(expansion, weights, source_transforms, n_grid))
  }
  offset <- 0
  non_point <- character()
  for(group in others){
    point <- .prior_density_ordinate_point_group_location(group, source_transforms)
    if(length(point) == 1L){
      if(is.na(point)){
        return(refuse("A point term added to the ordered terms is not a numeric location."))
      }
      offset <- offset + point
    }else{
      non_point <- c(non_point, group$parameter)
    }
  }
  additive_mean <- 0
  additive_sd <- 0
  if(length(non_point) > 0L){
    columns <- unlist(lapply(others[non_point], function(group) names(group$weights)), use.names = FALSE)
    normal <- .prior_density_ordinate_linear_normal(
      prior_list[non_point], weights[columns], source_transforms[columns], 0
    )
    if(is.null(normal) || !identical(normal$method, "linear_normal")){
      return(refuse(paste0(
        "The terms added to ordered terms that have no one-share representation must be normal ",
        "terms or points."
      )))
    }
    additive_mean <- normal$provenance$mean
    additive_sd <- normal$provenance$sd
  }

  slices <- tryCatch(
    lapply(ordered, function(group) .prior_ordered_linear_slices(group$prior, group$weights, group$indices)),
    error = function(e) e
  )
  if(inherits(slices, "error")){
    return(refuse(conditionMessage(slices)))
  }
  expansions <- lapply(ordered, function(group) .prior_ordered_total_components(group$prior$total))
  if(any(vapply(expansions, is.null, logical(1)))){
    return(refuse("The weights of an ordered mixture total are not exact numeric values."))
  }
  for(group in ordered){
    total <- group$prior$total
    metadata <- .prior_ordered_metadata(.prior_ordered_default_bound(group$prior))
    if(metadata$theta_dim > 1L && is.prior.mixture(total) && !is.prior.spike_and_slab(total)){
      return(refuse("A mixture total of several slices has no structural density route."))
    }
  }
  combinations <- expand.grid(lapply(expansions, function(expansion) seq_along(expansion$probabilities)),
                              KEEP.OUT.ATTRS = FALSE)
  if(nrow(combinations) > .prior_density_route_leaf_cap(n_grid)){
    return(refuse("The mixture expansion of the ordered totals exceeds the leaf cap.", "leaf_cap"))
  }
  routes <- lapply(seq_len(nrow(combinations)), function(row){
    components <- lapply(seq_along(expansions), function(g){
      expansions[[g]]$priors[[combinations[row, g]]]
    })
    .prior_ordered_product_leaf(
      ordered = ordered, slices = slices, components = components,
      others = others, prior_list = prior_list, weights = weights,
      source_transforms = source_transforms, offset = offset,
      additive_mean = additive_mean, additive_sd = additive_sd, n_grid = n_grid
    )
  })
  probabilities <- vapply(seq_len(nrow(combinations)), function(row){
    prod(vapply(seq_along(expansions), function(g){
      expansions[[g]]$probabilities[[combinations[row, g]]]
    }, numeric(1)))
  }, numeric(1))
  if(length(routes) == 1L){
    return(routes[[1L]])
  }
  .prior_density_route_mixture(components = routes, weights = probabilities)
}

# The route of one combination of total components.
.prior_ordered_product_leaf <- function(ordered, slices, components, others, prior_list, weights,
                                        source_transforms, offset, additive_mean, additive_sd, n_grid){

  refuse <- function(reason){
    .prior_density_route_unknown(
      reason = reason, provenance = list(kind = "unsupported_provenance", ordered = TRUE)
    )
  }
  terms <- list()
  for(g in seq_along(ordered)){
    total <- components[[g]]
    kind <- if(is.prior.none(total)){
      "point"
    }else if(is.prior.point(total)){
      if(.is_prior_expression(total) || !is.numeric(total$parameters$location) ||
         !is.finite(total$parameters$location)){
        return(refuse("The point total of an ordered term is not a numeric location."))
      }
      "point"
    }else if(.prior_density_simple_continuous(total)){
      "random"
    }else{
      return(refuse("The total of an ordered term in this combination must be a simple continuous prior or a point."))
    }
    location <- if(is.prior.none(total)) 0 else if(identical(kind, "point")) total$parameters$location else NA_real_
    normal <- identical(kind, "random") && .prior_ordered_product_normal(total)
    for(slice in slices[[g]]){
      if(all(slice$tensor$coefficients == 0) || (identical(kind, "point") && location == 0)){
        next
      }
      terms[[length(terms) + 1L]] <- list(
        kind = kind, total = total, location = location, tensor = slice$tensor,
        label = paste0(ordered[[g]]$parameter, "|", slice$slice),
        mean = if(normal) total$parameters$mean,
        sd = if(normal) total$parameters$sd
      )
    }
  }
  other_list <- prior_list[names(others)]
  other_weights <- weights[unlist(lapply(others, function(group) names(group$weights)), use.names = FALSE)]
  other_transforms <- source_transforms[names(other_weights)]
  if(length(terms) == 0L){
    return(.prior_density_route_linear(other_list, other_weights, other_transforms, n_grid))
  }
  structure <- .prior_ordered_product_structure(terms)
  p <- length(structure$coordinates)

  # no allocation stays random: independent totals with fixed coefficients
  if(p == 0L){
    for(j in seq_along(terms)){
      name <- paste0(".ordered_total[", terms[[j]]$label, "]#", j)
      other_list[[name]] <- terms[[j]]$total
      other_weights[[name]] <- terms[[j]]$tensor$coefficients[[1L]]
      other_transforms[[name]] <- NA_character_
    }
    return(.prior_density_route_linear(other_list, other_weights, other_transforms[names(other_weights)], n_grid))
  }

  for(j in seq_along(terms)){
    terms[[j]]$class_tensor <- .prior_ordered_product_class_tensor(terms[[j]], structure)
  }
  random <- which(vapply(terms, function(term) identical(term$kind, "random"), logical(1)))
  points <- which(vapply(terms, function(term) identical(term$kind, "point"), logical(1)))
  all_normal <- all(vapply(terms[random], function(term) !is.null(term$mean), logical(1)))
  mode <- if(length(random) == 0L){
    if(additive_sd > 0) "normal" else "deterministic"
  }else if(length(random) == 1L && additive_sd == 0){
    "scale"
  }else if(all_normal){
    "normal"
  }else{
    return(refuse(paste0(
      "Several totals in one ordered combination (slices or terms) have a closed-form law only when ",
      "they are normal priors, as do totals added to normal terms."
    )))
  }
  spec <- list(
    mode = mode, terms = terms, structure = structure, offset = offset,
    random_terms = random, point_terms = points,
    additive_mean = additive_mean, additive_sd = additive_sd
  )
  if(identical(mode, "deterministic")){
    if(p != 1L){
      return(refuse(paste0(
        "The levels of ordered terms with point totals are a nonlinear function of ", p,
        " allocation shares, whose density has no structural route."
      )))
    }
    state <- .prior_ordered_product_state(spec, matrix(c(0, 1), ncol = 1L))
    coordinate <- structure$coordinates[[1L]]
    scale <- state$shift[[2L]] - state$shift[[1L]]
    if(scale == 0){
      return(.prior_density_route_atom(
        locations = state$shift[[1L]], probability = 1,
        provenance = list(kind = "scalar_affine", offset = state$shift[[1L]], scale = 0)
      ))
    }
    return(list(
      type = "scalar", prior = prior("beta", list(alpha = coordinate$a, beta = coordinate$b)),
      offset = state$shift[[1L]], scale = scale, source_transform = NULL
    ))
  }
  if(length(.prior_ordered_product_sizes(p)) < 2L){
    return(refuse(paste0(
      "The combination involves ", p, " independent allocation coordinates, more than the tensor ",
      "quadrature of the product admits."
    )))
  }
  if(identical(mode, "scale")){
    spec$log_density <- .prior_simple_lpdf_evaluator(terms[[random]]$total, purpose = "natural_integral")
  }
  if(.prior_ordered_product_bounded_total(spec)){
    exponents <- .prior_ordered_product_bound_exponents(terms[[random]]$total)
    if(anyNA(exponents) || any(exponents < .prior_ordered_product_exponent_limit())){
      return(refuse(paste0(
        "The density of the total of an ordered term at a bound of its support is infinite with ",
        if(anyNA(exponents)) "an unknown exponent" else paste0("the exponent ", format(min(exponents))),
        ", which the pieces next to the bound do not integrate to the accuracy of an ordinate (the ",
        "exponent must be at least ", format(.prior_ordered_product_exponent_limit()), ")."
      )))
    }
  }
  if(.prior_ordered_product_bounded_total(spec) && p > .prior_ordered_product_piece_coordinates()){
    return(refuse(paste0(
      "A total with a bounded support in a product of ", p, " allocation coordinates has no verified ",
      "quadrature: its density has kinks where the argument crosses a bound, which the nested pieces ",
      "split for at most ", .prior_ordered_product_piece_coordinates(), " coordinates."
    )))
  }
  spec$hull <- .prior_ordered_product_hull(spec)
  list(type = "ordered_product", spec = spec, n_grid = n_grid)
}

# The refusal reason of a route that contains a leaf without a structural
# representation, NULL otherwise.
.prior_ordered_route_refusal <- function(route){

  if(is.null(route)){
    return("The ordered combination has no structural route.")
  }
  if(identical(route$type, "unknown")){
    return(if(is.null(route$reason)) "The ordered combination has no structural route." else route$reason)
  }
  if(identical(route$type, "mixture")){
    for(component in route$components[route$weights > 0]){
      reason <- .prior_ordered_route_refusal(component)
      if(!is.null(reason)){
        return(reason)
      }
    }
  }
  NULL
}

# The atoms of a route of an ordered combination (a data frame of locations and
# probabilities): the atom leaves of the mixture, weighted.
.prior_ordered_route_points <- function(route){

  points <- .prior_linear_density_empty_points()
  if(identical(route$type, "atom")){
    return(data.frame(x = route$locations, p = route$probability))
  }
  if(identical(route$type, "mixture")){
    for(i in which(route$weights > 0)){
      component <- .prior_ordered_route_points(route$components[[i]])
      if(nrow(component) > 0L){
        points <- rbind(points, data.frame(x = component$x, p = component$p * route$weights[[i]]))
      }
    }
    if(nrow(points) > 1L){
      points <- data.frame(x = sort(unique(points$x)),
                           p = as.vector(tapply(points$p, points$x, sum)))
    }
  }
  points
}

# The route of the combination of one ordered term alone (for ranges and
# grids), refused with the formula-measure condition when it has none.
.prior_ordered_group_route <- function(ordered_prior, weights, indices, n_grid){

  ordered_prior <- .prior_ordered_default_bound(ordered_prior)
  parameter <- .prior_ordered_metadata(ordered_prior)$parameter_name
  columns <- .prior_linear_prior_columns(parameter, ordered_prior)
  if(is.null(names(weights))){
    names(weights) <- columns[indices]
  }
  if(!all(names(weights) %in% columns)){
    stop("Internal ordered combination weights do not match the columns of the term.", call. = FALSE)
  }
  route <- .prior_density_route_linear(stats::setNames(list(ordered_prior), parameter), weights, NULL, n_grid)
  reason <- .prior_ordered_route_refusal(route)
  if(!is.null(reason)){
    .prior_ordered_density_unavailable(reason)
  }
  route
}

# Numerical range of a combination of an ordered term that is not one
# allocation share: the sum over its slices of the products of the range of the
# total and the range of the slice's coefficient. The coefficient is multilinear
# in the allocations, each on its simplex, so it attains its extremes at the
# vertices, where it equals an entry of the slice's reduced tensor: its range is
# the range of the tensor's coefficients (no enumeration of the shares).
.prior_ordered_product_range <- function(ordered_prior, weights, indices, tail_prob){

  slices <- .prior_ordered_linear_slices(ordered_prior, weights, indices)
  if(length(slices) == 0L){
    return(c(0, 0))
  }
  .prior_ordered_total_linear_measure_check(ordered_prior$total)
  total_range <- .prior_linear_group_range(
    list(prior = ordered_prior$total, weights = c(.ordered_total = 1), indices = 1L),
    tail_prob = tail_prob, source_transforms = c(.ordered_total = NA_character_)
  )
  ranges <- vapply(slices, function(slice){
    range(.prior_ordered_product_bound_products(range(slice$tensor$coefficients), total_range))
  }, numeric(2))
  c(sum(ranges[1L, ]), sum(ranges[2L, ]))
}

# The grid of the combination of one ordered term that is not one allocation
# share: its route's continuous density over the numerical range, with the
# atoms of the route (as the level of a Beta share, .prior_ordered_linear_distribution()).
.prior_ordered_product_distribution <- function(ordered_prior, weights, indices, n_grid, tail_prob){

  route <- .prior_ordered_group_route(ordered_prior, weights, indices, n_grid)
  # the density at every grid value is a tensor quadrature over the shares: the
  # grid of a product of three or more coordinates is coarser (display only)
  if(.prior_ordered_route_dimension(route) >= 3L){
    n_grid <- min(n_grid, 257L)
  }
  out <- .prior_linear_density_route_product(
    route  = route,
    range  = .prior_ordered_product_range(ordered_prior, weights, indices, tail_prob),
    points = .prior_ordered_route_points(route),
    n_grid = n_grid
  )
  if(is.null(out)){
    strong <- .prior_ordered_route_strong_singularity(route)
    .prior_ordered_density_unavailable(if(is.null(strong)){
      "its density has no resolved values over the numerical range."
    }else{
      paste0(strong, ".")
    })
  }
  attr(out, "ordered_measure") <- list(
    method = "analytic_components",
    total = ordered_prior$total,
    allocation = .prior_ordered_measure_allocation(ordered_prior, weights, indices),
    coefficient_weights = stats::setNames(as.numeric(weights), names(weights))
  )
  out
}

# The largest number of allocation coordinates of an ordered-product leaf of a route.
.prior_ordered_route_dimension <- function(route){

  if(identical(route$type, "ordered_product")){
    return(length(route$spec$structure$coordinates))
  }
  if(identical(route$type, "mixture")){
    return(max(0L, vapply(route$components[route$weights > 0], .prior_ordered_route_dimension, integer(1))))
  }
  0L
}

# The reason that an ordered-product leaf of a route has a strong singularity
# (NULL without one).
.prior_ordered_route_strong_singularity <- function(route){

  if(identical(route$type, "ordered_product")){
    return(.prior_ordered_product_strong_singularity(route$spec))
  }
  if(identical(route$type, "mixture")){
    for(component in route$components[route$weights > 0]){
      reason <- .prior_ordered_route_strong_singularity(component)
      if(!is.null(reason)){
        return(reason)
      }
    }
  }
  NULL
}
