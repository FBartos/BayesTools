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
  source_transforms <- if(is.null(source_transforms)){
    rep(NA_character_, length(weights))
  }else{
    source_transforms[names(weights)]
  }

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
    stop(
      "Continuous prior density is unavailable because its numerical range ",
      "collapses to one representable value. Center or rescale the modeled ",
      "quantity and its prior parameters before evaluating this density.",
      call. = FALSE
    )
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

  # An integrable singularity at a support bound (e.g., gamma or beta shape
  # below one) has no finite ordinate; its knot carries the exact mass of the
  # grid cell instead.
  infinite <- is.infinite(y) & y > 0
  if(any(infinite)){
    y[infinite] <- .prior_linear_scalar_cell_density(
      prior            = prior,
      x                = x[infinite],
      dx               = x[2] - x[1],
      weight           = weight,
      source_transform = source_transform
    )
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

.prior_linear_scalar_cell_density <- function(prior, x, dx, weight, source_transform = NULL){

  # Average density of weight * source over [x - dx / 2, x + dx / 2]. The
  # prior CDF restricts the cell to the prior support.
  source_bounds <- cbind(x - dx / 2, x + dx / 2) / weight
  if(identical(source_transform, "log")){
    source_bounds <- exp(source_bounds)
  }
  mass <- abs(
    mcdf(prior, source_bounds[, 2L]) - mcdf(prior, source_bounds[, 1L])
  )
  as.numeric(mass) / dx
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

.prior_ordered_linear_share <- function(ordered_prior, weights, indices){

  # Allocation share multiplying the ordered total in a level combination:
  # a fixed point or 'scale' times a Beta(alpha[1], alpha[2]) variable.
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
    return(list(type = "point", scale = scale))
  }
  if(!identical(record$spec$type, "dirichlet")){
    stop("Unsupported ordered allocation specification.", call. = FALSE)
  }

  unique_weights <- unique(allocation_weights)
  if(length(unique_weights) == 1L){
    return(list(type = "point", scale = unique_weights[[1L]]))
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
    return(list(type = "point", scale = 0))
  }
  if(alpha_remaining == 0){
    return(list(type = "point", scale = scale))
  }

  list(type = "beta", scale = scale, alpha = c(alpha_selected, alpha_remaining))
}

.prior_ordered_linear_multiplier <- function(ordered_prior, weights, indices,
                                             dx, n_grid, tail_prob){

  share <- .prior_ordered_linear_share(ordered_prior, weights, indices)
  if(identical(share$type, "point")){
    return(.prior_linear_density_point(share$scale))
  }

  beta_prior <- prior(
    "beta",
    list(alpha = share$alpha[[1L]], beta = share$alpha[[2L]])
  )
  beta_group <- list(
    prior = beta_prior,
    weights = c(.ordered_allocation = share$scale),
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

  if(is.character(transformation) && length(transformation) == 1L &&
     transformation %in% c("lin", "exp_lin") &&
     (is.null(transformation_arguments) || is.list(transformation_arguments))){
    a <- if(is.null(transformation_arguments[["a"]])){
      0
    }else{
      transformation_arguments[["a"]]
    }
    b <- if(is.null(transformation_arguments[["b"]])){
      1
    }else{
      transformation_arguments[["b"]]
    }
    if(is.numeric(a) && length(a) == 1L && is.finite(a) &&
       is.numeric(b) && length(b) == 1L && is.finite(b) && b == 0){
      location <- if(identical(transformation, "lin")) a else exp(a)
      if(!is.finite(location) || location == 0 &&
         identical(transformation, "exp_lin")){
        stop(
          "The constant prior-density transformation is not representable.",
          call. = FALSE
        )
      }
      return(.prior_linear_density_point(location))
    }
  }

  densities <- list()
  if(!is.null(dist$density) && dist$density$mass > 0){
    transformed <- .density.prior_transformation_grid(
      dist$density$x,
      dist$density$y,
      transformation,
      transformation_arguments
    )
    x_new <- transformed$x[!transformed$drop]
    y_new <- transformed$y[!transformed$drop]
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

      # Keep every distinct image knot. A tolerance on the spacing would drop
      # the dense knots near a saturating boundary (exp near zero, tanh near
      # +/-1), where the transformed density is largest.
      keep <- c(TRUE, diff(x_new) > 0)
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
    dist <- .prior_linear_density_convolve(
      dist,
      .prior_linear_density_on_spacing(group_dist, dx, n_grid),
      dx
    )
  }

  attr(dist, "weights") <- weights
  return(dist)
}

.prior_linear_density_on_spacing <- function(dist, dx, n_grid = NULL){

  # The FFT convolution assumes both continuous parts share the spacing 'dx'.
  # Groups built by the product quadrature (ordered totals times allocation
  # shares) carry their own grid and are resampled, preserving their mass.
  group_dx <- .prior_linear_density_dx(dist)
  if(!is.finite(group_dx) || !is.finite(dx) ||
     abs(group_dx - dx) <= .prior_linear_density_grid_tol() * dx){
    return(dist)
  }
  .prior_linear_density_regrid(dist, dx, n_grid)
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

.prior_conditional_normal_groups <- function(prior_list, split){

  if(length(split$product_groups) != 1L){
    return(NULL)
  }
  product <- split$product_groups[[1L]]
  multiplier <- prior_list[[product$multiplier]]
  if(!is.prior.simple(multiplier) || is.prior.point(multiplier) ||
     is.prior.discrete(multiplier) || is.prior.mixture(multiplier) ||
     is.prior.spike_and_slab(multiplier) ||
     !is.null(attr(multiplier, "multiply_by", exact = TRUE)) ||
     .prior_linear_prior_dimension(multiplier) != 1L){
    return(NULL)
  }
  bounds <- unlist(multiplier$truncation[c("lower", "upper")], use.names = FALSE)
  if(length(bounds) != 2L || anyNA(bounds) || bounds[1L] >= bounds[2L]){
    return(NULL)
  }
  additive_groups <- .prior_linear_weight_groups(prior_list, split$additive_weights)
  product_groups <- .prior_linear_weight_groups(product$prior_list, product$weights)
  if(length(intersect(names(additive_groups), names(product_groups))) > 0L ||
     product$multiplier %in% c(names(additive_groups), names(product_groups))){
    return(NULL)
  }
  list(multiplier = multiplier, bounds = bounds,
       additive = additive_groups, multiplied = product_groups)
}

.prior_conditional_normal_expansion <- function(prior_list, split, source_transforms, n_grid){

  groups <- .prior_conditional_normal_groups(prior_list, split)
  if(is.null(groups)) return(NULL)
  max_leaves <- floor(n_grid / .prior_conditional_normal_initial_evaluations(groups$bounds))
  if(max_leaves < 1) return(NULL)
  inspect_group <- function(group, limit){
    prior <- group$prior
    if(is.prior.mixture(prior) || is.prior.spike_and_slab(prior)){
      probabilities <- .prior_density_ordinate_mixture_weights(prior)
      if(is.null(probabilities)) return(NULL)
      indices <- which(probabilities > 0)
      count <- 0
      always_normal <- TRUE
      zero_points <- TRUE
      branch_counts <- numeric(length(indices))
      for(j in seq_along(indices)){
        child <- group
        child$prior <- .prior_density_copy_parent_attributes(prior[[indices[j]]], prior)
        if(!identical(attr(child$prior, "multiply_by", exact = TRUE),
                      attr(prior, "multiply_by", exact = TRUE))) return(NULL)
        info <- inspect_group(child, limit - count)
        if(is.null(info)) return(NULL)
        count <- count + info$count
        branch_counts[j] <- info$count
        always_normal <- always_normal && info$always_normal
        zero_points <- zero_points && info$zero_points
      }
      return(list(count = count, always_normal = always_normal,
                  zero_points = zero_points, indices = indices,
                  branch_counts = branch_counts))
    }
    if(limit < 1) return(NULL)
    normal <- .prior_density_ordinate_linear_normal(
      stats::setNames(list(prior), group$parameter), group$weights, source_transforms, 0
    )
    if(!is.null(normal) && isTRUE(normal$exact) && identical(normal$method, "linear_normal")){
      return(list(count = 1, always_normal = TRUE, zero_points = TRUE))
    }
    point <- .prior_density_ordinate_point_group_location(group, source_transforms)
    if(length(point) == 1L && is.finite(point)){
      return(list(count = 1, always_normal = FALSE, zero_points = point == 0))
    }
    NULL
  }
  all_groups <- c(groups$additive, groups$multiplied)
  counts <- list()
  total <- 1
  for(parameter in names(all_groups)){
    info <- inspect_group(all_groups[[parameter]], floor(max_leaves / total))
    if(is.null(info)) return(NULL)
    total <- total * info$count
    counts[[parameter]] <- info
  }
  additive_normal <- vapply(counts[names(groups$additive)], `[[`, logical(1), "always_normal")
  product_normal <- vapply(counts[names(groups$multiplied)], `[[`, logical(1), "always_normal")
  product_zero <- vapply(counts[names(groups$multiplied)], `[[`, logical(1), "zero_points")
  if(!any(additive_normal) || (!any(product_normal) && !all(product_zero))) return(NULL)
  mixture_names <- names(counts)[vapply(counts, function(x) !is.null(x$indices), logical(1))]
  if(length(mixture_names) == 0L) return(NULL)
  parameter <- mixture_names[[1L]]
  first <- counts[[parameter]]
  list(parameter = parameter, indices = first$indices,
       budgets = floor(n_grid / total) * (total / first$count) * first$branch_counts)
}

.prior_conditional_normal_spec <- function(prior_list, split, source_transforms){

  groups <- .prior_conditional_normal_groups(prior_list, split)
  if(is.null(groups)) return(NULL)
  product <- split$product_groups[[1L]]
  additive <- .prior_density_ordinate_linear_normal(
    prior_list, split$additive_weights, source_transforms, 0
  )
  multiplied <- .prior_density_ordinate_linear_normal(
    product$prior_list, product$weights, source_transforms, 0
  )
  if(is.null(additive) || is.null(multiplied) ||
     !identical(additive$method, "linear_normal") ||
     !identical(multiplied$method, "linear_normal") ||
     !isTRUE(additive$exact) || !isTRUE(multiplied$exact) ||
     !is.finite(additive$provenance$sd) || additive$provenance$sd <= 0 ||
     !is.finite(multiplied$provenance$sd) || multiplied$provenance$sd <= 0){
    return(NULL)
  }
  list(
    additive_mean = additive$provenance$mean,
    additive_sd = additive$provenance$sd,
    product_mean = multiplied$provenance$mean,
    product_sd = multiplied$provenance$sd,
    multiplier = groups$multiplier,
    bounds = groups$bounds,
    sources = list(additive = names(groups$additive),
                   multiplied = names(groups$multiplied), multiplier = product$multiplier)
  )
}

.prior_conditional_normal_ordinate <- function(spec, value, n_grid){

  # Finite bounds use 21-point rules; infinite bounds use 15-point rules,
  # with paired evaluations for two-sided infinite intervals.
  initial_evaluations <- .prior_conditional_normal_initial_evaluations(spec$bounds)
  subdivisions <- floor((n_grid + initial_evaluations) / (2 * initial_evaluations))
  tolerance <- .prior_linear_density_refinement_tolerance()
  evaluations <- 0L
  integrand <- function(multiplier){
    if(evaluations + length(multiplier) > n_grid){
      stop("the integration evaluation budget was exhausted", call. = FALSE)
    }
    evaluations <<- evaluations + length(multiplier)
    product_sd <- abs(multiplier) * spec$product_sd
    scale <- pmax(spec$additive_sd, product_sd)
    conditional_sd <- scale * sqrt((spec$additive_sd / scale)^2 + (product_sd / scale)^2)
    conditional_mean <- spec$additive_mean + multiplier * spec$product_mean
    exp(stats::dnorm(value, conditional_mean, conditional_sd, log = TRUE) +
          lpdf(spec$multiplier, multiplier))
  }
  integral <- if(!is.finite(subdivisions) || subdivisions < 1){
    list(value = NA_real_, abs.error = NA_real_, subdivisions = 0L,
         message = paste0("fewer than ", initial_evaluations,
                          " integration evaluations are available"))
  }else tryCatch(
    stats::integrate(integrand, spec$bounds[1L], spec$bounds[2L],
                     subdivisions = subdivisions, rel.tol = tolerance$relative,
                     abs.tol = tolerance$absolute, stop.on.error = FALSE),
    error = function(e){
      list(value = NA_real_, abs.error = NA_real_, subdivisions = 0L,
           message = conditionMessage(e))
    }
  )
  if(identical(integral$message, "OK") && isTRUE(integral$value == 0)){
    integral$message <- "zero ordinate for a structurally positive density"
  }
  bound <- tolerance$absolute + tolerance$relative * abs(integral$value)
  accepted <- identical(integral$message, "OK") && is.finite(integral$value) &&
    integral$value > 0 && is.finite(integral$abs.error) && integral$abs.error <= bound
  if(!isTRUE(accepted)){
    integral$value <- NA_real_
  }
  .prior_density_ordinate_result(
    value = value, behavior = "regular",
    log_density = log(integral$value), exact = TRUE,
    method = "conditional_normal_mixture",
    provenance = list(
      kind = "conditional_normal_mixture",
      additive = c(mean = spec$additive_mean, sd = spec$additive_sd),
      multiplied = c(mean = spec$product_mean, sd = spec$product_sd),
      multiplier = .prior_density_ordinate_prior_provenance(spec$multiplier),
      independent_sources = spec$sources,
      structural_regularity = "positive_variance_gaussian_convolution",
      integration = list(
        kind = "conditional_normal_mixture", exact = FALSE,
        absolute_error = integral$abs.error, error_bound = bound,
        evaluations = evaluations,
        budget = n_grid, converged = isTRUE(accepted), message = integral$message
      )
    )
  )
}

.prior_conditional_normal_initial_evaluations <- function(bounds){

  if(all(is.finite(bounds))) return(21L)
  if(all(is.infinite(bounds))) return(30L)
  15L
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

.prior_linear_density_refinement <- function(x){

  context <- attr(x, "adaptive_evaluation", exact = TRUE)
  if(is.null(context) ||
     !context$kind %in%
       c("linear_combination", "density_context", "density_context_rows")){
    return(NULL)
  }

  arguments <- context$arguments
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
  if(isTRUE(all.equal(arguments, context$arguments, tolerance = 0))){
    return(NULL)
  }
  refined_arguments <- arguments
  refined_arguments$.record_evaluation <- FALSE
  refined <- if(identical(context$kind, "linear_combination")){
    do.call(.prior_linear_combination_density, refined_arguments)
  }else if(identical(context$kind, "density_context_rows")){
    do.call(.prior_density_from_context_rows, refined_arguments)
  }else{
    do.call(.prior_density_from_context, refined_arguments)
  }
  context$arguments <- arguments
  attr(refined, "adaptive_evaluation") <- context
  attr(refined, "refinement_settings") <- list(
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

  refined
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

  dist                <- x
  transformed_x_range <- NULL
  if(!is.null(transformation) && transformation_settings && !is.null(x_range)){
    transformed_x_range <- x_range
    x_range <- NULL
  }

  out <- list()

  if(!is.null(dist$density) && dist$density$mass > 0){
    if(!is.null(transformed_x_range)){
      x_den <- seq(transformed_x_range[1], transformed_x_range[2], length.out = n_points)
      x_raw <- suppressWarnings(.density.prior_transformation_inv_grid(
        x_den,
        transformation,
        transformation_arguments
      ))
    }else if(is.null(x_range)){
      x_raw <- seq(min(dist$density$x), max(dist$density$x), length.out = n_points)
    }else{
      x_raw <- seq(x_range[1], x_range[2], length.out = n_points)
    }

    finite_raw <- is.finite(x_raw)
    y_den      <- rep(NA_real_, length(x_raw))
    evaluator  <- attr(dist, "density_evaluator", exact = TRUE)
    if(any(finite_raw) && is.function(evaluator)){
      evaluated <- evaluator(x_raw[finite_raw])
      if(!is.numeric(evaluated) || length(evaluated) != sum(finite_raw) ||
         anyNA(evaluated)){
        stop("The analytic prior density evaluator returned invalid values.",
             call. = FALSE)
      }
      y_den[finite_raw] <- evaluated
    }else if(any(finite_raw)){
      y_den[finite_raw] <- stats::approx(
        dist$density$x,
        dist$density$y,
        xout   = x_raw[finite_raw],
        yleft  = 0,
        yright = 0
      )$y
    }
    y_den <- y_den * dist$density$mass

    if(!is.null(transformation)){
      if(is.null(transformed_x_range)){
        x_den <- .density.prior_transformation_x(
          x_raw,
          transformation,
          transformation_arguments
        )
      }
      y_transformed <- rep(NA_real_, length(y_den))
      y_transformed[finite_raw] <- .density.prior_transformation_y(
        x_den[finite_raw],
        y_den[finite_raw],
        transformation,
        transformation_arguments
      )
      y_den <- y_transformed
    }else{
      x_den <- x_raw
    }

    finite <- is.finite(x_den) & is.finite(y_den)
    x_den  <- x_den[finite]
    y_den  <- y_den[finite]

    if(length(x_den) > 0L){
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
      attr(out_den, "x_range") <- if(is.null(transformed_x_range)) range(x_den) else transformed_x_range
      attr(out_den, "y_range") <- c(0, max(y_den, 0, na.rm = TRUE))
      if(!is.null(level)) attr(out_den, "level") <- level
      if(!is.null(level_name)) attr(out_den, "level_name") <- level_name
      out[["density"]] <- out_den
    }
  }

  points <- dist$points
  if(!is.null(points) && nrow(points) > 0){
    points <- points[points$p > 0, , drop = FALSE]
    if(!is.null(x_range) && is.null(transformed_x_range)){
      points <- points[points$x >= min(x_range) & points$x <= max(x_range), , drop = FALSE]
    }
    if(nrow(points) > 0 && !is.null(transformation)){
      points$x <- .density.prior_transformation_x(points$x, transformation, transformation_arguments)
    }
    if(nrow(points) > 0 && !is.null(transformed_x_range)){
      points <- points[
        points$x >= min(transformed_x_range) &
        points$x <= max(transformed_x_range),
        ,
        drop = FALSE
      ]
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
      attr(out_point, "x_range") <- if(is.null(transformed_x_range)) range(points$x) else transformed_x_range
      attr(out_point, "y_range") <- c(0, max(points$p))
      if(!is.null(level)) attr(out_point, "level") <- level
      if(!is.null(level_name)) attr(out_point, "level_name") <- level_name
      out[[paste0("points", i)]] <- out_point
    }
  }

  return(out)
}

.prior_linear_group_support_hull <- function(group, source_transforms = NULL){

  # Closed interval containing the support of one weighted prior group, or
  # NULL when that support is not known exactly from the prior definition.
  prior   <- group$prior
  weights <- group$weights
  transforms <- if(is.null(source_transforms)){
    rep(NA_character_, length(weights))
  }else{
    unname(source_transforms[names(weights)])
  }

  if(is.prior.none(prior)){
    return(c(0, 0))
  }
  if(is.prior.ordered(prior) || is.prior.weightfunction(prior) ||
     is_prior_phacking(prior) || is_prior_bias(prior)){
    return(NULL)
  }
  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    components <- if(is.prior.spike_and_slab(prior)){
      list(.get_spike_and_slab_variable(prior), prior("point", list(location = 0)))
    }else{
      probabilities <- .prior_density_ordinate_mixture_weights(prior)
      if(is.null(probabilities)) prior else prior[probabilities > 0]
    }
    hulls <- lapply(components, function(component){
      component_group <- group
      component_group$prior <- component
      .prior_linear_group_support_hull(component_group, source_transforms)
    })
    if(length(hulls) == 0L || any(vapply(hulls, is.null, logical(1)))){
      return(NULL)
    }
    return(range(unlist(hulls)))
  }
  if(is.prior.vector(prior) && !is.prior.treatment(prior) && !is.prior.independent(prior)){
    if(any(!is.na(transforms))){
      return(NULL)
    }
    if(identical(prior$distribution, "mpoint")){
      location <- prior$parameters[["location"]]
      if(!is.numeric(location) || length(location) != 1L || !is.finite(location)){
        return(NULL)
      }
      return(rep(sum(weights) * location, 2L))
    }
    if(prior$distribution %in% c("mnormal", "mt")){
      return(c(-Inf, Inf))
    }
    return(NULL)
  }
  if(!is.prior.simple(prior)){
    return(NULL)
  }

  bounds <- if(is.prior.point(prior)){
    rep(prior$parameters[["location"]], 2L)
  }else if(is.prior.discrete(prior)){
    range(.prior_simple_truncated_discrete(prior)$support)
  }else{
    c(prior$truncation[["lower"]], prior$truncation[["upper"]])
  }
  if(!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds)){
    return(NULL)
  }

  hull <- c(0, 0)
  for(i in seq_along(weights)){
    source_bounds <- bounds
    if(!is.na(transforms[i])){
      if(!identical(transforms[i], "log") || source_bounds[1L] < 0){
        return(NULL)
      }
      source_bounds <- c(
        if(source_bounds[1L] == 0) -Inf else log(source_bounds[1L]),
        log(source_bounds[2L])
      )
    }
    hull <- hull + range(weights[[i]] * source_bounds)
  }
  hull
}

.prior_linear_combination_support_hull <- function(prior_list, weights,
                                                   source_transforms = NULL){

  weights <- weights[is.finite(weights) & weights != 0]
  if(length(weights) == 0L){
    return(c(0, 0))
  }
  split <- tryCatch(
    .prior_linear_split_multiply_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(split) || length(split$product_groups) > 0L){
    return(NULL)
  }
  additive <- split$additive_weights[split$additive_weights != 0]
  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, additive),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }
  hull <- c(0, 0)
  for(group in groups){
    group_hull <- .prior_linear_group_support_hull(group, source_transforms)
    if(is.null(group_hull)){
      return(NULL)
    }
    hull <- hull + group_hull
  }
  hull
}

.prior_linear_context_support_hull <- function(context, weights,
                                               source_transforms = NULL){

  union_hull <- function(hulls){
    if(length(hulls) == 0L || any(vapply(hulls, is.null, logical(1)))){
      return(NULL)
    }
    range(unlist(hulls))
  }

  if(inherits(context, "prior_density_context")){
    standardized <- tryCatch(
      .prior_density_context_standardized_weights(context, weights),
      error = function(e) NULL
    )
    if(is.null(standardized)){
      return(NULL)
    }
    return(.prior_linear_combination_support_hull(
      context$prior_list,
      standardized,
      if(is.null(source_transforms)) NULL else source_transforms[names(standardized)]
    ))
  }
  if(inherits(context, "prior_density_model_mixture_context")){
    models <- which(is.finite(context$model_weights) & context$model_weights > 0)
    return(union_hull(lapply(models, function(model_i){
      model_prior_list <- lapply(context$prior_list, function(parameter_priors){
        if(is.prior(parameter_priors)) parameter_priors else parameter_priors[[model_i]]
      })
      names(model_prior_list) <- names(context$prior_list)
      for(parameter in names(model_prior_list)){
        if(is.null(model_prior_list[[parameter]])){
          model_prior_list[[parameter]] <- prior("point", list(location = 0))
        }
      }
      .prior_linear_combination_support_hull(model_prior_list, weights, source_transforms)
    })))
  }
  if(inherits(context, "prior_density_conditional_context")){
    models <- which(is.finite(context$model_weights) & context$model_weights > 0)
    return(union_hull(lapply(models, function(model_i){
      prior_list <- context$prior_lists[[model_i]]
      if(!is.null(context$formula_scale) && length(context$formula_scale) > 0L){
        component_context <- .prior_density_context(
          prior_list    = prior_list,
          column_names  = context$column_names,
          formula_scale = context$formula_scale,
          n_grid        = context$n_grid,
          tail_prob     = context$tail_prob
        )
        return(.prior_linear_context_support_hull(component_context, weights, source_transforms))
      }
      .prior_linear_combination_support_hull(prior_list, weights, source_transforms)
    })))
  }
  NULL
}

.prior_linear_density_support_hull <- function(adaptive){

  # Closed interval containing the support of a recorded prior density, from
  # the prior definitions and weights only; NULL when it is not known exactly.
  if(!is.list(adaptive) || !is.character(adaptive$kind) ||
     length(adaptive$kind) != 1L || !is.list(adaptive$arguments)){
    return(NULL)
  }
  arguments <- adaptive$arguments
  hull <- switch(
    adaptive$kind,
    "linear_combination" = .prior_linear_combination_support_hull(
      arguments$prior_list,
      arguments$weights,
      arguments$source_transforms
    ),
    "density_context" = .prior_linear_context_support_hull(
      arguments$context,
      arguments$weights,
      arguments$source_transforms
    ),
    "density_context_rows" = {
      weights <- arguments$weights
      if(is.null(dim(weights))){
        .prior_linear_context_support_hull(
          arguments$context, weights, arguments$source_transforms
        )
      }else{
        hulls <- lapply(seq_len(nrow(weights)), function(row_i){
          .prior_linear_context_support_hull(
            arguments$context, weights[row_i, ], arguments$source_transforms
          )
        })
        if(length(hulls) == 0L || any(vapply(hulls, is.null, logical(1)))){
          NULL
        }else{
          range(unlist(hulls))
        }
      }
    },
    NULL
  )
  if(is.null(hull) || anyNA(hull)){
    return(NULL)
  }

  transformation <- arguments$output_transformation
  if(is.null(transformation)){
    return(hull)
  }
  if(!is.character(transformation) || length(transformation) != 1L ||
     !transformation %in% c("lin", "exp", "exp_lin", "tanh")){
    return(NULL)
  }
  transformation_arguments <- arguments$output_transformation_arguments
  if(transformation %in% c("lin", "exp_lin")){
    b <- transformation_arguments[["b"]]
    if(!is.null(b) && (!is.numeric(b) || length(b) != 1L || !is.finite(b) || b == 0)){
      return(NULL)
    }
  }
  if(identical(transformation, "exp_lin") && hull[1L] < 0){
    return(NULL)
  }
  mapped <- suppressWarnings(.density.prior_transformation_x(
    hull, transformation, transformation_arguments
  ))
  if(anyNA(mapped)){
    return(NULL)
  }
  range(mapped)
}

.prior_linear_density_height <- function(x, value){

  if(!inherits(x, "prior_linear_density")){
    stop("'x' must be a prior linear density object.", call. = FALSE)
  }

  ordinate <- NULL
  if(length(value) == 1L && is.finite(value)){
    ordinate <- .prior_density_ordinate_from_adaptive(
      attr(x, "adaptive_evaluation", exact = TRUE), value
    )
    continuous <- if(is.null(ordinate)){
      NULL
    }else{
      .prior_density_ordinate_continuous_behavior(ordinate)
    }
    if(identical(continuous, "regular") &&
       .prior_density_ordinate_has_quadrature(ordinate$provenance)){
      integration <- .prior_density_ordinate_integration(ordinate$provenance)
      if(is.na(ordinate$log_density) || !isTRUE(integration$converged)){
        stop("Conditional-normal prior density was rejected by diagnostics: integration reported '",
             integration$message, "' with absolute error ", format(integration$absolute_error),
             ". Inspect the prior specification and increase 'n_samples' for marginal inference.",
             call. = FALSE)
      }
      height <- exp(ordinate$log_density)
      attr(height, "numerical_diagnostics") <- integration
      return(height)
    }
    # Exact structural ordinates (including exact finite mixtures) are used
    # directly; grid refinement cannot converge across a density jump.
    if(!is.null(ordinate) && isTRUE(ordinate$exact)){
      if(identical(continuous, "infinite")){
        return(Inf)
      }
      if(identical(continuous, "zero")){
        return(0)
      }
      if(identical(continuous, "regular") && !is.na(ordinate$log_density)){
        return(exp(ordinate$log_density))
      }
    }
  }

  singular_points <- attr(x, "singular_density_points", exact = TRUE)
  if(length(value) == 1L && value %in% singular_points){
    if(!is.null(ordinate) && isTRUE(ordinate$exact) &&
       .prior_density_ordinate_continuous_behavior(ordinate) %in% c("regular", "zero") &&
       !is.na(ordinate$log_density)){
      return(exp(ordinate$log_density))
    }
    stop("The prior density at the flagged product ordinate is unavailable from supported deterministic provenance. Inspect the prior specification or use a supported prior-density evaluator.",
         call. = FALSE)
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

  support <- .prior_linear_density_support_hull(
    attr(x, "adaptive_evaluation", exact = TRUE)
  )
  if(!is.null(support) && (value < support[1L] || value > support[2L])){
    return(0)
  }

  height <- .prior_linear_density_grid_height(x, value)
  refined <- .prior_linear_density_refinement(x)
  if(is.null(refined)){
    if(!is.null(attr(x, "adaptive_evaluation", exact = TRUE))){
      stop(
        "Adaptive prior-density evaluation did not converge within the documented ",
        "grid-refinement error criterion.",
        call. = FALSE
      )
    }
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
  for(i in seq_len(4L)){
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
    if(i < 4L){
      next_refined <- .prior_linear_density_refinement(refined)
      if(is.null(next_refined)){
        break
      }
      refined <- next_refined
    }
  }

  final_range <- .prior_linear_density_range(refined)
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
