.prior_linear_is_factor_prior <- function(prior){
  is.prior.factor(prior) ||
    inherits(prior, "prior.factor_mixture") ||
    inherits(prior, "prior.factor_spike_and_slab")
}

.prior_linear_prior_dimension <- function(prior){

  if(.prior_linear_is_factor_prior(prior)){
    if(!is.null(attr(prior, "K")) && !is.na(attr(prior, "K"))){
      return(attr(prior, "K"))
    }
    if(is.prior.mixture(prior)){
      factor_components <- prior[vapply(prior, function(p) is.prior.factor(p) || inherits(p, "prior.factor_spike_and_slab"), logical(1))]
      if(length(factor_components) > 0){
        return(.get_prior_factor_levels(factor_components[[1]]))
      }
    }
    return(.get_prior_factor_levels(prior))
  }

  if(is.prior.vector(prior)){
    return(prior$parameters[["K"]])
  }

  1L
}

.prior_linear_prior_columns <- function(parameter, prior){

  K <- .prior_linear_prior_dimension(prior)
  if(length(K) != 1 || is.na(K) || K <= 1){
    return(parameter)
  }

  paste0(parameter, "[", seq_len(K), "]")
}

.prior_linear_representative_prior <- function(prior){

  if(is.prior(prior)){
    return(prior)
  }

  if(is.list(prior)){
    prior <- prior[vapply(prior, is.prior, logical(1))]
    if(length(prior) > 0){
      return(prior[[1]])
    }
  }

  NULL
}

.prior_linear_active_parameters <- function(prior_list, weights){

  active <- character()
  if(is.null(names(weights))){
    return(active)
  }

  weights <- weights[is.finite(weights)]
  weights <- weights[weights != 0]
  if(length(weights) == 0){
    return(active)
  }

  for(parameter in names(prior_list)){
    prior <- .prior_linear_representative_prior(prior_list[[parameter]])
    if(is.null(prior)){
      next
    }
    columns <- .prior_linear_prior_columns(parameter, prior)
    present <- intersect(columns, names(weights))
    if(length(present) > 0){
      active <- c(active, parameter)
    }
  }

  active
}

.prior_linear_active_conditionals <- function(prior_list, weights, conditional){

  .condition_event_active_labels(
    prior_list  = prior_list,
    weights     = weights,
    conditional = conditional
  )
}

.prior_linear_weight_groups <- function(prior_list, weights){

  if(is.null(names(weights))){
    stop("'weights' must be a named numeric vector.", call. = FALSE)
  }

  groups <- list()
  matched <- rep(FALSE, length(weights))
  names(matched) <- names(weights)

  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    if(is.null(prior)){
      next
    }
    if(!is.prior(prior)){
      stop("All entries of 'prior_list' must be prior objects.", call. = FALSE)
    }

    columns <- .prior_linear_prior_columns(parameter, prior)
    present <- intersect(columns, names(weights))
    present <- present[weights[present] != 0]
    if(length(present) == 0){
      next
    }

    groups[[parameter]] <- list(
      parameter = parameter,
      prior     = prior,
      columns   = columns,
      weights   = weights[present],
      indices   = match(present, columns)
    )
    matched[present] <- TRUE
  }

  unmatched <- names(weights)[weights != 0 & !matched]
  if(length(unmatched) > 0){
    stop(
      "No prior distribution was found for coefficient column(s): '",
      paste0(unmatched, collapse = "', '"), "'.",
      call. = FALSE
    )
  }

  groups
}

.prior_linear_scalar_range <- function(prior, weight, tail_prob, source_transform = NULL){

  source_transform <- .prior_linear_source_transform(source_transform)

  if(weight == 0 || is.prior.none(prior)){
    return(c(0, 0))
  }

  if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    if(identical(source_transform, "log")){
      if(location <= 0){
        stop("A log-transformed prior source must have positive support.", call. = FALSE)
      }
      location <- log(location)
    }
    return(rep(weight * location, 2))
  }

  if(is.prior.discrete(prior)){
    support <- switch(
      prior[["distribution"]],
      "bernoulli" = c(0, 1),
      stop("Unsupported discrete prior distribution for linear-combination densities.", call. = FALSE)
    )
    if(identical(source_transform, "log")){
      if(any(support <= 0)){
        stop("A log-transformed prior source must have positive support.", call. = FALSE)
      }
      support <- log(support)
    }
    return(range(weight * support))
  }

  if(identical(source_transform, "log")){
    prior_range <- quant(prior, c(tail_prob, 1 - tail_prob))
    if(any(!is.finite(prior_range)) || any(prior_range <= 0)){
      prior_range <- range(prior, quantiles = tail_prob)
      prior_range[prior_range <= 0] <- min(prior_range[prior_range > 0], na.rm = TRUE)
    }
    if(any(!is.finite(prior_range)) || any(prior_range <= 0)){
      stop("A log-transformed prior source must have positive support.", call. = FALSE)
    }
    prior_range <- log(prior_range)
  }else{
    prior_range <- range(prior, quantiles = tail_prob)
  }

  range(weight * prior_range)
}

.prior_linear_vector_scalar_prior <- function(prior, weights){

  weights <- as.numeric(weights)
  norm_weight <- sqrt(sum(weights^2))

  if(norm_weight == 0){
    return(prior("point", list(location = 0)))
  }

  if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    return(prior("point", list(location = sum(weights) * location)))
  }

  switch(
    prior[["distribution"]],
    "mnormal" = prior(
      "normal",
      list(
        mean = sum(weights) * prior$parameters[["mean"]],
        sd   = norm_weight * prior$parameters[["sd"]]
      )
    ),
    "mt" = prior(
      "t",
      list(
        location = sum(weights) * prior$parameters[["location"]],
        scale    = norm_weight * prior$parameters[["scale"]],
        df       = prior$parameters[["df"]]
      )
    ),
    "mpoint" = prior(
      "point",
      list(location = sum(weights) * prior$parameters[["location"]])
    ),
    stop("Unsupported vector prior distribution for linear-combination densities.", call. = FALSE)
  )
}

.prior_linear_group_range <- function(group, tail_prob, source_transforms = NULL){

  prior <- group$prior
  weights <- group$weights
  source_transforms <- source_transforms[names(weights)]

  if(is.prior.none(prior)){
    return(c(0, 0))
  }

  if(is.prior.ordered(prior)){
    return(.prior_ordered_linear_range(
      ordered_prior = prior,
      weights = weights,
      indices = group$indices,
      tail_prob = tail_prob
    ))
  }

  if(is.prior.spike_and_slab(prior)){
    variable_prior <- .get_spike_and_slab_variable(prior)
    variable_group <- group
    variable_group$prior <- variable_prior
    variable_range <- .prior_linear_group_range(variable_group, tail_prob, source_transforms)
    return(range(c(0, variable_range)))
  }

  if(is.prior.mixture(prior)){
    ranges <- do.call(rbind, lapply(prior, function(component){
      component_group <- group
      component_group$prior <- component
      .prior_linear_group_range(component_group, tail_prob, source_transforms)
    }))
    return(range(ranges))
  }

  if(is.prior.vector(prior) && !is.prior.treatment(prior) && !is.prior.independent(prior)){
    if(any(!is.na(source_transforms))){
      stop("Source transformations are only supported for scalar prior terms.", call. = FALSE)
    }
    scalar_prior <- .prior_linear_vector_scalar_prior(prior, weights)
    return(.prior_linear_scalar_range(scalar_prior, 1, tail_prob))
  }

  ranges <- do.call(rbind, Map(function(weight, source_transform){
    .prior_linear_scalar_range(prior, weight, tail_prob, source_transform)
  }, weights, source_transforms))

  c(sum(ranges[, 1]), sum(ranges[, 2]))
}

.prior_linear_scalar_distribution <- function(prior, weight, dx, tail_prob, source_transform = NULL, n_grid = NULL){

  source_transform <- .prior_linear_source_transform(source_transform)

  if(weight == 0 || is.prior.none(prior)){
    return(.prior_linear_density_point(0))
  }

  if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    if(identical(source_transform, "log")){
      if(location <= 0){
        stop("A log-transformed prior source must have positive support.", call. = FALSE)
      }
      location <- log(location)
    }
    return(.prior_linear_density_point(weight * location))
  }

  if(is.prior.discrete(prior)){
    support <- switch(
      prior[["distribution"]],
      "bernoulli" = c(0, 1),
      stop("Unsupported discrete prior distribution for linear-combination densities.", call. = FALSE)
    )
    probs <- mpdf(prior, support)
    if(identical(source_transform, "log")){
      if(any(support <= 0)){
        stop("A log-transformed prior source must have positive support.", call. = FALSE)
      }
      support <- log(support)
    }
    points <- data.frame(x = weight * support, p = probs / sum(probs))
    return(.prior_linear_density_coalesce(points = points, dx = dx, n_grid = n_grid))
  }

  x_range <- .prior_linear_scalar_range(prior, weight, tail_prob, source_transform)
  if(x_range[1] == x_range[2]){
    return(.prior_linear_density_point(x_range[1]))
  }

  x <- seq(x_range[1], x_range[2], by = dx)
  if(length(x) < 3){
    x <- seq(x_range[1], x_range[2], length.out = 3)
  }

  source_x <- x / weight
  if(identical(source_transform, "log")){
    original_x <- exp(source_x)
    y <- mpdf(prior, original_x) * original_x / abs(weight)
  }else{
    y <- mpdf(prior, source_x) / abs(weight)
  }

  if(any(!is.finite(y))){
    stop(
      "A continuous prior component produced a non-finite density on its ",
      "numerical grid; the value cannot be replaced by zero faithfully.",
      call. = FALSE
    )
  }
  area <- sum(y) * (x[2] - x[1])
  if(!is.finite(area) || area <= 0){
    stop("A continuous prior component evaluated to zero mass on its grid.", call. = FALSE)
  }
  y <- y / area

  .prior_linear_density_coalesce(
    densities = list(list(x = x, y = y, mass = 1)),
    dx = dx,
    n_grid = n_grid
  )
}

.prior_linear_group_distribution <- function(group, dx, tail_prob, source_transforms = NULL, n_grid = NULL){

  prior <- group$prior
  weights <- group$weights
  source_transforms <- source_transforms[names(weights)]

  if(is.prior.none(prior)){
    return(.prior_linear_density_point(0))
  }

  if(is.prior.ordered(prior)){
    return(.prior_ordered_linear_distribution(
      ordered_prior = prior,
      weights = weights,
      indices = group$indices,
      dx = dx,
      n_grid = n_grid,
      tail_prob = tail_prob
    ))
  }

  if(is.prior.spike_and_slab(prior)){
    variable_prior  <- .get_spike_and_slab_variable(prior)
    inclusion_prior <- .get_spike_and_slab_inclusion(prior)
    inclusion <- mean(inclusion_prior)
    inclusion <- min(max(inclusion, 0), 1)

    variable_group <- group
    variable_group$prior <- variable_prior

    return(.prior_linear_density_mix(
      dists = list(
        .prior_linear_group_distribution(variable_group, dx, tail_prob, source_transforms, n_grid),
        .prior_linear_density_point(0)
      ),
      weights = c(inclusion, 1 - inclusion),
      dx = dx,
      n_grid = n_grid
    ))
  }

  if(is.prior.mixture(prior)){
    prior_weights <- attr(prior, "prior_weights")
    prior_weights <- prior_weights / sum(prior_weights)
    dists <- lapply(prior, function(component){
      component_group <- group
      component_group$prior <- component
      .prior_linear_group_distribution(component_group, dx, tail_prob, source_transforms, n_grid)
    })
    return(.prior_linear_density_mix(dists, prior_weights, dx = dx, n_grid = n_grid))
  }

  if(is.prior.vector(prior) && !is.prior.treatment(prior) && !is.prior.independent(prior)){
    if(any(!is.na(source_transforms))){
      stop("Source transformations are only supported for scalar prior terms.", call. = FALSE)
    }
    scalar_prior <- .prior_linear_vector_scalar_prior(prior, weights)
    return(.prior_linear_scalar_distribution(scalar_prior, 1, dx, tail_prob, n_grid = n_grid))
  }

  dist <- .prior_linear_density_point(0)
  for(i in seq_along(weights)){
    source_dist <- .prior_linear_scalar_distribution(
      prior             = prior,
      weight            = weights[i],
      dx                = dx,
      tail_prob         = tail_prob,
      source_transform  = source_transforms[i],
      n_grid            = n_grid
    )
    dist <- .prior_linear_density_convolve(dist, source_dist, dx)
  }

  dist
}

.prior_ordered_total_linear_distribution <- function(total, dx, n_grid,
                                                      tail_prob){

  total_group <- list(
    prior = total,
    weights = c(.ordered_total = 1),
    indices = 1L
  )
  .prior_linear_group_distribution(
    group = total_group,
    dx = dx,
    tail_prob = tail_prob,
    source_transforms = c(.ordered_total = NA_character_),
    n_grid = n_grid
  )
}

.prior_ordered_linear_range <- function(ordered_prior, weights, indices,
                                        tail_prob){

  ordered_prior <- .prior_ordered_default_bound(ordered_prior)
  total_group <- list(
    prior = ordered_prior$total,
    weights = c(.ordered_total = 1),
    indices = 1L
  )
  total_range <- .prior_linear_group_range(
    total_group,
    tail_prob = tail_prob,
    source_transforms = c(.ordered_total = NA_character_)
  )
  scale <- max(c(1, abs(weights), diff(total_range)), na.rm = TRUE)
  n_grid <- .prior_linear_density_default_grid()
  multiplier <- .prior_ordered_linear_multiplier(
    ordered_prior = ordered_prior,
    weights = weights,
    indices = indices,
    dx = scale / max(1, n_grid - 1L),
    n_grid = n_grid,
    tail_prob = tail_prob
  )
  multiplier_range <- .prior_linear_density_range(multiplier)
  products <- as.vector(outer(total_range, multiplier_range, `*`))
  products <- products[is.finite(products)]
  if(length(products) == 0L){
    return(c(0, 0))
  }
  range(products)
}

.prior_ordered_linear_multiplier <- function(ordered_prior, weights, indices,
                                             dx, n_grid, tail_prob){

  ordered_prior <- .prior_ordered_default_bound(ordered_prior)
  metadata <- .prior_ordered_metadata(ordered_prior)
  if(length(metadata$ordered_terms) != 1L ||
     metadata$theta_dim != 1L ||
     length(metadata$allocations) != 1L){
    stop(
      "Mixed-measure ordered densities currently require one ordered term ",
      "and one scalar total. Split the interaction into explicitly named ",
      "terms before requesting its marginal density.",
      call. = FALSE
    )
  }

  coefficient_weights <- numeric(metadata$coefficient_dim)
  coefficient_weights[indices] <- weights
  record <- metadata$allocations[[1L]]
  increments <- metadata$coefficient_grid[[metadata$ordered_terms[[1L]]]]
  allocation_weights <- numeric(length(record$spec$weights))
  if(identical(record$spec$type, "dirichlet")){
    allocation_weights <- numeric(length(record$spec$alpha))
  }
  for(i in seq_along(coefficient_weights)){
    allocation_weights[increments[[i]]] <-
      allocation_weights[increments[[i]]] + coefficient_weights[[i]]
  }

  if(identical(record$spec$type, "fixed")){
    scale <- sum(allocation_weights * record$spec$weights)
    return(.prior_linear_density_point(scale))
  }
  if(!identical(record$spec$type, "dirichlet")){
    stop("Unsupported ordered allocation specification.", call. = FALSE)
  }

  unique_weights <- unique(allocation_weights)
  if(length(unique_weights) == 1L){
    return(.prior_linear_density_point(unique_weights[[1L]]))
  }

  nonzero <- allocation_weights != 0
  nonzero_values <- unique(allocation_weights[nonzero])
  if(length(nonzero_values) != 1L){
    stop(
      "This ordered-prior linear combination is not a level or a single ",
      "allocation subset, so its Dirichlet multiplier has no beta reduction. ",
      "Request factor-level marginals instead.",
      call. = FALSE
    )
  }

  scale <- nonzero_values[[1L]]
  alpha_selected <- sum(record$spec$alpha[nonzero])
  alpha_remaining <- sum(record$spec$alpha[!nonzero])
  if(alpha_selected == 0){
    return(.prior_linear_density_point(0))
  }
  if(alpha_remaining == 0){
    return(.prior_linear_density_point(scale))
  }

  beta_prior <- prior(
    "beta",
    list(alpha = alpha_selected, beta = alpha_remaining)
  )
  beta_group <- list(
    prior = beta_prior,
    weights = c(.ordered_allocation = scale),
    indices = 1L
  )
  .prior_linear_group_distribution(
    group = beta_group,
    dx = dx,
    tail_prob = tail_prob,
    source_transforms = c(.ordered_allocation = NA_character_),
    n_grid = n_grid
  )
}

.prior_ordered_linear_distribution <- function(ordered_prior, weights, indices,
                                               dx = NA_real_, n_grid = NULL,
                                               tail_prob = .prior_linear_density_tail_prob()){

  if(is.null(n_grid)){
    n_grid <- .prior_linear_density_default_grid()
  }
  multiplier <- .prior_ordered_linear_multiplier(
    ordered_prior = ordered_prior,
    weights = weights,
    indices = indices,
    dx = dx,
    n_grid = n_grid,
    tail_prob = tail_prob
  )
  total <- .prior_ordered_total_linear_distribution(
    total = ordered_prior$total,
    dx = dx,
    n_grid = n_grid,
    tail_prob = tail_prob
  )

  out <- .prior_linear_density_product(
    total,
    multiplier,
    n_grid = n_grid
  )
  attr(out, "ordered_measure") <- list(
    method = "analytic_components",
    total = ordered_prior$total,
    allocation = .prior_ordered_metadata(
      .prior_ordered_default_bound(ordered_prior)
    )$allocations[[1L]]$spec,
    coefficient_weights = stats::setNames(
      as.numeric(weights),
      names(weights)
    )
  )
  out
}

.prior_linear_density_transform <- function(dist, transformation, transformation_arguments = NULL, n_grid = NULL){

  if(is.null(transformation)){
    return(dist)
  }

  densities <- list()
  if(!is.null(dist$density) && dist$density$mass > 0){
    x_old <- dist$density$x
    y_old <- dist$density$y
    x_new <- .density.prior_transformation_x(x_old, transformation, transformation_arguments)
    y_new <- .density.prior_transformation_y(x_new, y_old, transformation, transformation_arguments)

    if(any(!is.finite(x_new)) || any(!is.finite(y_new))){
      stop(
        "The requested prior-density transformation produced non-finite ",
        "grid values.",
        call. = FALSE
      )
    }

    if(length(x_new) >= 2){
      ord <- order(x_new)
      x_new <- x_new[ord]
      y_new <- y_new[ord]

      keep <- c(TRUE, diff(x_new) > .prior_linear_density_grid_tol() * pmax(1, abs(x_new[-length(x_new)])))
      x_new <- x_new[keep]
      y_new <- y_new[keep]

      area <- if(length(x_new) > 1){
        sum(diff(x_new) * (head(y_new, -1) + tail(y_new, -1)) / 2)
      }else{
        0
      }
      if(is.finite(area) && area > 0){
        y_new <- y_new / area
      }
      densities[[1]] <- list(x = x_new, y = y_new, mass = dist$density$mass)
    }
  }

  points <- dist$points
  if(!is.null(points) && nrow(points) > 0){
    points$x <- .density.prior_transformation_x(points$x, transformation, transformation_arguments)
    if(any(!is.finite(points$x))){
      stop(
        "The requested prior-density transformation produced a non-finite ",
        "point-mass location.",
        call. = FALSE
      )
    }
  }

  if(length(densities) == 0){
    return(.prior_linear_density_coalesce(
      points = points,
      dx     = NA_real_,
      n_grid = if(is.null(n_grid)) dist$n_grid else n_grid
    ))
  }

  out <- list(
    density = densities[[1]],
    points  = .prior_linear_density_aggregate_points(points, NA_real_),
    n_grid  = if(is.null(n_grid)) dist$n_grid else n_grid
  )
  class(out) <- c("prior_linear_density", "prior_density")
  .prior_linear_density_normalize(out, warn = TRUE)
}

.prior_linear_split_multiply_groups <- function(prior_list, weights){

  additive_weights <- weights
  product_groups <- list()

  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    if(is.null(prior)){
      next
    }
    if(!is.prior(prior)){
      stop("All entries of 'prior_list' must be prior objects.", call. = FALSE)
    }

    columns <- .prior_linear_prior_columns(parameter, prior)
    present <- intersect(columns, names(additive_weights))
    present <- present[additive_weights[present] != 0]
    if(length(present) == 0){
      next
    }

    multiply_by <- attr(prior, "multiply_by")
    if(is.null(multiply_by)){
      next
    }

    if(is.numeric(multiply_by)){
      if(length(multiply_by) != 1L){
        stop("Numeric 'multiply_by' values must be scalar for deterministic prior densities.", call. = FALSE)
      }
      additive_weights[present] <- additive_weights[present] * multiply_by
      next
    }

    if(!is.character(multiply_by) || length(multiply_by) != 1L){
      stop("'multiply_by' must be either a scalar numeric value or a scalar parameter name.", call. = FALSE)
    }

    if(is.null(product_groups[[multiply_by]])){
      product_groups[[multiply_by]] <- list(
        multiplier = multiply_by,
        prior_list = list(),
        weights    = numeric()
      )
    }

    product_groups[[multiply_by]]$prior_list[[parameter]] <- prior
    product_groups[[multiply_by]]$weights <- c(
      product_groups[[multiply_by]]$weights,
      additive_weights[present]
    )
    additive_weights[present] <- 0
  }

  list(
    additive_weights = additive_weights,
    product_groups   = product_groups
  )
}

.prior_linear_additive_combination_density <- function(prior_list, weights,
                                                       n_grid = .prior_linear_density_default_grid(),
                                                       tail_prob = .prior_linear_density_tail_prob(),
                                                       source_transforms = NULL){

  weights <- weights[is.finite(weights)]
  weights <- weights[weights != 0]

  if(length(weights) == 0){
    return(.prior_linear_density_point(0))
  }

  groups <- .prior_linear_weight_groups(prior_list, weights)
  if(length(groups) == 0){
    return(.prior_linear_density_point(0))
  }

  if(is.null(source_transforms)){
    source_transforms <- rep(NA_character_, length(weights))
    names(source_transforms) <- names(weights)
  }else{
    source_transforms <- source_transforms[names(weights)]
    source_transforms[is.na(source_transforms)] <- NA_character_
  }

  group_ranges <- do.call(rbind, lapply(groups, .prior_linear_group_range,
                                        tail_prob = tail_prob,
                                        source_transforms = source_transforms))
  target_range <- c(sum(group_ranges[, 1]), sum(group_ranges[, 2]))
  target_width <- diff(target_range)
  dx <- target_width / max(1, n_grid - 1)
  if(!is.finite(dx) || dx <= 0){
    dx <- 1
  }

  dist <- .prior_linear_density_point(0)
  for(group in groups){
    group_dist <- .prior_linear_group_distribution(
      group             = group,
      dx                = dx,
      tail_prob         = tail_prob,
      source_transforms = source_transforms,
      n_grid            = n_grid
    )
    dist <- .prior_linear_density_convolve(dist, group_dist, dx)
  }

  attr(dist, "weights") <- weights
  return(dist)
}

.prior_linear_combination_density <- function(prior_list, weights,
                                              n_grid = .prior_linear_density_default_grid(),
                                              tail_prob = .prior_linear_density_tail_prob(),
                                              source_transforms = NULL,
                                              output_transformation = NULL,
                                              output_transformation_arguments = NULL,
                                              .record_evaluation = TRUE){

  check_list(prior_list, "prior_list")
  if(is.null(weights) || length(weights) == 0){
    weights <- numeric()
  }else{
    check_real(weights, "weights", check_length = 0, allow_NA = FALSE)
    if(any(!is.finite(weights))){
      stop("The 'weights' argument must contain only finite values.", call. = FALSE)
    }
  }
  check_int(n_grid, "n_grid", lower = 16)
  check_real(tail_prob, "tail_prob", lower = 0, upper = 0.5, allow_bound = FALSE)

  weights <- weights[weights != 0]

  if(length(weights) == 0){
    dist <- .prior_linear_density_point(0)
    return(.prior_linear_density_transform(dist, output_transformation, output_transformation_arguments, n_grid))
  }

  if(is.null(source_transforms)){
    source_transforms <- rep(NA_character_, length(weights))
    names(source_transforms) <- names(weights)
  }else{
    source_transforms <- source_transforms[names(weights)]
    source_transforms[is.na(source_transforms)] <- NA_character_
  }

  split <- .prior_linear_split_multiply_groups(prior_list, weights)
  components <- list(
    .prior_linear_additive_combination_density(
      prior_list         = prior_list,
      weights            = split$additive_weights,
      n_grid             = n_grid,
      tail_prob          = tail_prob,
      source_transforms  = source_transforms
    )
  )
  singular_density_points <- numeric()

  if(length(split$product_groups) > 0){
    for(product_group in split$product_groups){
      multiplier <- product_group$multiplier
      if(!multiplier %in% names(prior_list)){
        stop("No prior distribution was found for 'multiply_by' parameter '", multiplier, "'.", call. = FALSE)
      }
      multiplier_prior <- prior_list[[multiplier]]
      if(is.null(multiplier_prior) || !is.prior(multiplier_prior)){
        stop("The 'multiply_by' parameter '", multiplier, "' does not have a supported prior distribution.", call. = FALSE)
      }
      if(!is.null(attr(multiplier_prior, "multiply_by"))){
        stop("Nested 'multiply_by' prior densities are not supported.", call. = FALSE)
      }
      if(.prior_linear_prior_dimension(multiplier_prior) != 1L){
        stop("The 'multiply_by' parameter '", multiplier, "' must have a scalar prior distribution.", call. = FALSE)
      }

      linear_dist <- .prior_linear_additive_combination_density(
        prior_list         = product_group$prior_list,
        weights            = product_group$weights,
        n_grid             = n_grid,
        tail_prob          = tail_prob,
        source_transforms  = source_transforms
      )

      multiplier_weights <- 1
      names(multiplier_weights) <- multiplier
      multiplier_dist <- .prior_linear_additive_combination_density(
        prior_list         = prior_list[multiplier],
        weights            = multiplier_weights,
        n_grid             = n_grid,
        tail_prob          = tail_prob,
        source_transforms  = source_transforms
      )
      if(.prior_linear_density_grid_height(linear_dist, 0) > 0 &&
         .prior_linear_density_grid_height(multiplier_dist, 0) > 0){
        singular_density_points <- c(singular_density_points, 0)
      }

      components[[length(components) + 1L]] <- .prior_linear_density_product(
        linear_dist,
        multiplier_dist,
        n_grid = n_grid
      )
    }
  }

  dist <- .prior_linear_density_sum_independent(components, n_grid = n_grid)
  dist <- .prior_linear_density_transform(dist, output_transformation,
                                          output_transformation_arguments, n_grid)
  attr(dist, "weights") <- weights
  if(length(singular_density_points) > 0L){
    attr(dist, "singular_density_points") <-
      sort(unique(singular_density_points))
  }
  if(isTRUE(.record_evaluation)){
    attr(dist, "adaptive_evaluation") <- list(
      kind = "linear_combination",
      arguments = list(
        prior_list = prior_list,
        weights = weights,
        n_grid = n_grid,
        tail_prob = tail_prob,
        source_transforms = source_transforms,
        output_transformation = output_transformation,
        output_transformation_arguments = output_transformation_arguments
      )
    )
    attr(dist, "numerical_diagnostics") <- list(
      n_grid = n_grid,
      tail_probability_per_source = tail_prob,
      intended_captured_probability_per_continuous_source =
        max(0, 1 - 2 * tail_prob),
      numerical_range = .prior_linear_density_range(dist),
      grid_normalization =
        attr(dist, "grid_normalization", exact = TRUE),
      fft_clipping =
        attr(dist, "fft_clipping", exact = TRUE),
      adaptive_evaluation = TRUE
    )
  }

  if(length(weights) == 1L && is.null(output_transformation)){
    parameter <- names(weights)
    scalar_prior <- if(length(parameter) == 1L) prior_list[[parameter]] else NULL
    source_transform <- source_transforms[parameter]
    if(is.prior.simple(scalar_prior) &&
       !is.prior.mixture(scalar_prior) &&
       !is.prior.spike_and_slab(scalar_prior) &&
       !is.prior.point(scalar_prior) &&
       !is.prior.discrete(scalar_prior) &&
       is.null(attr(scalar_prior, "multiply_by", exact = TRUE)) &&
       (length(source_transform) == 0L || is.na(source_transform))){
      weight <- unname(weights[[1L]])
      attr(dist, "density_evaluator") <- local({
        prior_value <- scalar_prior
        weight_value <- weight
        function(value){
          mpdf(prior_value, value / weight_value) / abs(weight_value)
        }
      })
    }
  }
  return(.prior_linear_density_normalize(dist, warn = TRUE))
}

.prior_linear_density_refinement_tolerance <- function(){

  list(relative = 1e-4, absolute = 1e-12)
}

.prior_linear_density_grid_height <- function(x, value){

  height <- 0
  if(!is.null(x$density) &&
     value >= min(x$density$x) && value <= max(x$density$x)){
    height <- stats::approx(
      x$density$x,
      x$density$y * x$density$mass,
      xout = value,
      yleft = 0,
      yright = 0
    )$y
  }
  height
}

.prior_linear_density_refinements <- function(x, max_refinements = 4L){

  context <- attr(x, "adaptive_evaluation", exact = TRUE)
  if(is.null(context) ||
     !context$kind %in%
       c("linear_combination", "density_context", "density_context_rows")){
    return(list())
  }

  arguments <- context$arguments
  refinements <- vector("list", max_refinements)
  for(i in seq_len(max_refinements)){
    if(identical(context$kind, "linear_combination")){
      arguments$n_grid <- min(
        max(as.integer(arguments$n_grid) * 2L, 4096L),
        32768L
      )
      arguments$tail_prob <- max(arguments$tail_prob / 1000, 1e-12)
    }else{
      arguments$context$n_grid <- min(
        max(as.integer(arguments$context$n_grid) * 2L, 4096L),
        32768L
      )
      arguments$context$tail_prob <- max(
        arguments$context$tail_prob / 1000,
        1e-12
      )
    }
    refined_arguments <- arguments
    refined_arguments$.record_evaluation <- FALSE
    refinements[[i]] <- if(identical(context$kind, "linear_combination")){
      do.call(.prior_linear_combination_density, refined_arguments)
    }else if(identical(context$kind, "density_context_rows")){
      do.call(.prior_density_from_context_rows, refined_arguments)
    }else{
      do.call(.prior_density_from_context, refined_arguments)
    }
    attr(refinements[[i]], "refinement_settings") <- list(
      n_grid = if(identical(context$kind, "linear_combination")){
        arguments$n_grid
      }else{
        arguments$context$n_grid
      },
      tail_prob = if(identical(context$kind, "linear_combination")){
        arguments$tail_prob
      }else{
        arguments$context$tail_prob
      }
    )
  }

  refinements
}

.prior_linear_density_to_plot_data <- function(x, n_points = 1000, x_range = NULL,
                                               transformation = NULL,
                                               transformation_arguments = NULL,
                                               transformation_settings = FALSE,
                                               factor = FALSE,
                                               level = NULL,
                                               level_name = NULL){

  if(!inherits(x, "prior_linear_density")){
    stop("'x' must be a prior linear density object.", call. = FALSE)
  }

  check_int(n_points, "n_points", lower = 2)
  check_real(x_range, "x_range", check_length = 2, allow_NULL = TRUE)
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)

  dist <- x
  if(!is.null(transformation) && transformation_settings && !is.null(x_range)){
    x_range <- .density.prior_transformation_inv_x(x_range, transformation, transformation_arguments)
  }

  out <- list()

  if(!is.null(dist$density) && dist$density$mass > 0){
    if(is.null(x_range)){
      x_den <- seq(min(dist$density$x), max(dist$density$x), length.out = n_points)
    }else{
      x_den <- seq(x_range[1], x_range[2], length.out = n_points)
    }
    y_den <- stats::approx(dist$density$x, dist$density$y, xout = x_den, yleft = 0, yright = 0)$y
    y_den <- y_den * dist$density$mass

    if(!is.null(transformation)){
      x_den <- .density.prior_transformation_x(x_den, transformation, transformation_arguments)
      y_den <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
    }

    out_den <- list(
      call    = call("density", "linear-combination prior"),
      bw      = NULL,
      n       = n_points,
      x       = x_den,
      y       = y_den,
      samples = NULL
    )
    class(out_den) <- c("density", "density.prior", "density.prior.simple",
                        if(factor) "density.prior.factor")
    attr(out_den, "x_range") <- range(x_den)
    attr(out_den, "y_range") <- c(0, max(y_den, 0, na.rm = TRUE))
    if(!is.null(level)) attr(out_den, "level") <- level
    if(!is.null(level_name)) attr(out_den, "level_name") <- level_name
    out[["density"]] <- out_den
  }

  points <- dist$points
  if(!is.null(points) && nrow(points) > 0){
    points <- points[points$p > 0, , drop = FALSE]
    if(!is.null(x_range)){
      points <- points[points$x >= min(x_range) & points$x <= max(x_range), , drop = FALSE]
    }
    if(nrow(points) > 0 && !is.null(transformation)){
      points$x <- .density.prior_transformation_x(points$x, transformation, transformation_arguments)
    }

    for(i in seq_len(nrow(points))){
      out_point <- list(
        call    = call("density", paste0("point", i)),
        bw      = NULL,
        n       = n_points,
        x       = points$x[i],
        y       = points$p[i],
        samples = NULL
      )
      class(out_point) <- c("density", "density.prior", "density.prior.point",
                            if(factor) "density.prior.factor")
      attr(out_point, "x_range") <- range(points$x)
      attr(out_point, "y_range") <- c(0, max(points$p))
      if(!is.null(level)) attr(out_point, "level") <- level
      if(!is.null(level_name)) attr(out_point, "level_name") <- level_name
      out[[paste0("points", i)]] <- out_point
    }
  }

  return(out)
}

.prior_linear_density_height <- function(x, value){

  if(!inherits(x, "prior_linear_density")){
    stop("'x' must be a prior linear density object.", call. = FALSE)
  }

  singular_points <- attr(x, "singular_density_points", exact = TRUE)
  if(length(value) == 1L && value %in% singular_points){
    return(Inf)
  }

  evaluator <- attr(x, "density_evaluator", exact = TRUE)
  if(is.function(evaluator)){
    height <- evaluator(value)
    if(!is.numeric(height) || length(height) != length(value) || anyNA(height)){
      stop("The analytic prior density evaluator returned invalid values.", call. = FALSE)
    }
    return(height)
  }

  if(length(value) != 1L || !is.finite(value)){
    stop("Adaptive prior-density evaluation requires one finite ordinate.",
         call. = FALSE)
  }

  height <- .prior_linear_density_grid_height(x, value)
  refinements <- .prior_linear_density_refinements(x)
  if(length(refinements) == 0L){
    if(!is.null(x$density) &&
       (value < min(x$density$x) || value > max(x$density$x))){
      stop(
        "The requested ordinate is outside the numerical approximation range, ",
        "and the density has no provenance for adaptive extension.",
        call. = FALSE
      )
    }
    return(height)
  }

  tolerance <- .prior_linear_density_refinement_tolerance()
  previous <- height
  for(i in seq_along(refinements)){
    refined <- refinements[[i]]
    current <- .prior_linear_density_grid_height(refined, value)
    change <- abs(current - previous)
    bound <- tolerance$absolute +
      tolerance$relative * max(abs(current), abs(previous))
    inside <- !is.null(refined$density) &&
      value >= min(refined$density$x) &&
      value <= max(refined$density$x)
    if(isTRUE(inside) && is.finite(current) && change <= bound){
      settings <- attr(refined, "refinement_settings", exact = TRUE)
      attr(current, "adaptive_evaluation") <- c(
        settings,
        list(
          refinements = i,
          absolute_change = change,
          error_bound = bound,
          numerical_range = .prior_linear_density_range(refined),
          converged = TRUE
        )
      )
      return(current)
    }
    previous <- current
  }

  final_range <- .prior_linear_density_range(
    refinements[[length(refinements)]]
  )
  if(value < final_range[1L] || value > final_range[2L]){
    stop(
      "The requested ordinate remains outside the numerical approximation ",
      "range after adaptive extension.",
      call. = FALSE
    )
  }
  stop(
    "Adaptive prior-density evaluation did not converge within the documented ",
    "grid-refinement error criterion.",
    call. = FALSE
  )
}

.prior_linear_density_point_mass <- function(x, value){

  if(!inherits(x, "prior_linear_density") || is.null(x$points) || nrow(x$points) == 0){
    return(0)
  }

  # Exact IEEE equality by design: point locations are not snapped or
  # coalesced to nearby grid / cut values.
  sum(x$points$p[x$points$x == value])
}
