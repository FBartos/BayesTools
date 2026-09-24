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

  source_density <- function(x){
    source_x <- x / weight
    if(identical(source_transform, "log")){
      original_x <- exp(source_x)
      mpdf(prior, original_x) * original_x / abs(weight)
    }else{
      mpdf(prior, source_x) / abs(weight)
    }
  }

  # 'seq(by = dx)' can stop one rounding error short of a range end that is a
  # multiple of dx. When that end is a support bound with an infinite density,
  # the knot belongs on the bound, where it carries its cell mass below;
  # left inside, it would carry a finite but arbitrarily large ordinate.
  last <- length(x)
  if(x[last] < x_range[2] &&
     x_range[2] - x[last] <= 64 * .Machine$double.eps * max(abs(x_range))){
    bound_density <- as.numeric(source_density(x_range[2]))
    if(isTRUE(is.infinite(bound_density) && bound_density > 0)){
      x[last] <- x_range[2]
    }
  }

  y <- source_density(x)

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
    out <- .prior_linear_density_coalesce(
      points = points,
      dx     = NA_real_,
      n_grid = if(is.null(n_grid)) dist$n_grid else n_grid
    )
    attr(out, "grid_resolution") <- attr(dist, "grid_resolution", exact = TRUE)
    return(out)
  }

  out <- list(
    density = densities[[1]],
    points  = .prior_linear_density_aggregate_points(points, NA_real_),
    n_grid  = if(is.null(n_grid)) dist$n_grid else n_grid
  )
  class(out) <- c("prior_linear_density", "prior_density")
  # Adaptive refinement halves the spacing of the source (linear-predictor)
  # grid, which the transformed knots no longer show.
  attr(out, "grid_resolution") <- attr(dist, "grid_resolution", exact = TRUE)
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

.prior_linear_group_robust_scale <- function(group, source_transforms = NULL){

  # Smallest robust scale (IQR / 1.349 in output units) of the continuous
  # sources of one weighted prior group; Inf when the group has none. Ordered
  # groups are resolved by their own product grid and are not included.
  prior   <- group$prior
  weights <- group$weights
  transforms <- if(is.null(source_transforms)){
    rep(NA_character_, length(weights))
  }else{
    unname(source_transforms[names(weights)])
  }

  if(is.prior.none(prior) || is.prior.ordered(prior) ||
     is.prior.point(prior) || is.prior.discrete(prior)){
    return(Inf)
  }
  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    components <- if(is.prior.spike_and_slab(prior)){
      list(.get_spike_and_slab_variable(prior))
    }else{
      probabilities <- .prior_density_ordinate_mixture_weights(prior)
      if(is.null(probabilities)) prior else prior[probabilities > 0]
    }
    scales <- vapply(components, function(component){
      component_group <- group
      component_group$prior <- component
      .prior_linear_group_robust_scale(component_group, source_transforms)
    }, numeric(1))
    return(min(c(Inf, scales)))
  }
  if(is.prior.vector(prior) && !is.prior.treatment(prior) && !is.prior.independent(prior)){
    scalar_prior <- tryCatch(
      .prior_linear_vector_scalar_prior(prior, weights),
      error = function(e) NULL
    )
    if(is.null(scalar_prior)){
      return(Inf)
    }
    return(.prior_linear_group_robust_scale(
      list(prior = scalar_prior, weights = c(.vector = 1))
    ))
  }
  if(!is.prior.simple(prior)){
    return(Inf)
  }

  quartiles <- mquant(prior, c(.25, .75))
  scales <- vapply(seq_along(weights), function(i){
    source_quartiles <- if(identical(transforms[i], "log")) log(quartiles) else quartiles
    abs(weights[[i]]) * (source_quartiles[2] - source_quartiles[1]) / 1.349
  }, numeric(1))
  scales <- scales[is.finite(scales) & scales > 0]
  if(length(scales) == 0L) Inf else min(scales)
}

.prior_linear_additive_combination_density <- function(prior_list, weights,
                                                       n_grid = .prior_linear_density_default_grid(),
                                                       tail_prob = .prior_linear_density_tail_prob(),
                                                       source_transforms = NULL,
                                                       grid_spacing = NULL){

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
  # Resolve the narrowest continuous source even when a heavy-tailed source
  # sets the range.
  robust_scale <- min(c(Inf, vapply(groups, .prior_linear_group_robust_scale,
                                    numeric(1), source_transforms = source_transforms)))
  resolution <- .prior_linear_density_resolution(
    width        = target_width,
    n_grid       = n_grid,
    scale        = robust_scale,
    grid_spacing = grid_spacing
  )
  n_grid <- resolution$n_grid
  dx <- resolution$dx
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
  attr(dist, "grid_resolution") <- c(spacing = dx, n_grid = n_grid)
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
                                              grid_spacing = NULL,
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
  check_real(grid_spacing, "grid_spacing", lower = 0, allow_bound = FALSE, allow_NULL = TRUE)

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
      source_transforms  = source_transforms,
      grid_spacing       = grid_spacing
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
      # The product and the multiplier's own term share one random variable,
      # so they cannot be convolved as independent components.
      multiplier_columns <- .prior_linear_prior_columns(multiplier, multiplier_prior)
      if(any(split$additive_weights[intersect(multiplier_columns, names(split$additive_weights))] != 0)){
        stop(
          "The prior density of this linear combination is unavailable because '",
          multiplier, "' enters it both as the 'multiply_by' scale of other ",
          "coefficients and with its own weight, which makes the terms dependent. ",
          "Evaluate the terms separately.",
          call. = FALSE
        )
      }

      linear_dist <- .prior_linear_additive_combination_density(
        prior_list         = product_group$prior_list,
        weights            = product_group$weights,
        n_grid             = n_grid,
        tail_prob          = tail_prob,
        source_transforms  = source_transforms,
        grid_spacing       = grid_spacing
      )

      multiplier_weights <- 1
      names(multiplier_weights) <- multiplier
      multiplier_dist <- .prior_linear_additive_combination_density(
        prior_list         = prior_list[multiplier],
        weights            = multiplier_weights,
        n_grid             = n_grid,
        tail_prob          = tail_prob,
        source_transforms  = source_transforms,
        grid_spacing       = grid_spacing
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

  dist <- .prior_linear_density_sum_independent(
    components,
    n_grid       = max(c(n_grid, vapply(components, function(component){
      as.integer(component$n_grid)
    }, integer(1)))),
    grid_spacing = grid_spacing
  )
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
        output_transformation_arguments = output_transformation_arguments,
        grid_spacing = grid_spacing
      )
    )
    attr(dist, "numerical_diagnostics") <- list(
      n_grid = n_grid,
      grid_resolution = attr(dist, "grid_resolution", exact = TRUE),
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

# A simple continuous scalar prior with a numeric interval support: the
# multiplier of a scale mixture or a term of a two-term convolution.
.prior_density_simple_continuous <- function(prior){

  if(!is.prior(prior) || !is.prior.simple(prior) || is.prior.point(prior) ||
     is.prior.discrete(prior) || is.prior.mixture(prior) ||
     is.prior.spike_and_slab(prior) || is.prior.vector(prior) ||
     .is_prior_expression(prior) ||
     !is.null(attr(prior, "multiply_by", exact = TRUE)) ||
     !.prior_density_ordinate_parameters_numeric(prior)){
    return(FALSE)
  }
  bounds <- unlist(prior$truncation[c("lower", "upper")], use.names = FALSE)
  is.numeric(bounds) && length(bounds) == 2L && !anyNA(bounds) && bounds[1L] < bounds[2L]
}

# Conditional-normal specification of a product term X = a + b * s: an
# additive part a (normal terms and points; a_s = 0 when it is deterministic),
# a multiplied normal part b ~ N(b_m, b_s) with b_s > 0 and a simple continuous
# multiplier s. Given s, X is normal with mean a_m + b_m s and SD
# sqrt(a_s^2 + b_s^2 s^2); with a_s = 0 it is a pure scale mixture. NULL for
# other products.
.prior_conditional_normal_spec <- function(prior_list, split, source_transforms){

  if(length(split$product_groups) != 1L){
    return(NULL)
  }
  product <- split$product_groups[[1L]]
  multiplier <- prior_list[[product$multiplier]]
  if(!.prior_density_simple_continuous(multiplier) ||
     .prior_linear_prior_dimension(multiplier) != 1L){
    return(NULL)
  }
  bounds <- unlist(multiplier$truncation[c("lower", "upper")], use.names = FALSE)
  additive_groups <- .prior_linear_weight_groups(prior_list, split$additive_weights)
  product_groups <- .prior_linear_weight_groups(product$prior_list, product$weights)
  if(length(intersect(names(additive_groups), names(product_groups))) > 0L ||
     product$multiplier %in% c(names(additive_groups), names(product_groups))){
    return(NULL)
  }
  multiplied <- .prior_density_ordinate_linear_normal(
    product$prior_list, product$weights, source_transforms, 0
  )
  if(is.null(multiplied) || !identical(multiplied$method, "linear_normal") ||
     !isTRUE(multiplied$exact) ||
     !is.finite(multiplied$provenance$sd) || multiplied$provenance$sd <= 0){
    return(NULL)
  }
  additive <- .prior_density_ordinate_linear_normal(
    prior_list, split$additive_weights, source_transforms, 0
  )
  if(is.null(additive)){
    # a deterministic additive part (points only): a pure scale mixture
    additive_mean <- .prior_density_ordinate_deterministic_offset(
      prior_list, split$additive_weights, source_transforms
    )
    if(is.null(additive_mean) || !is.finite(additive_mean)){
      return(NULL)
    }
    additive_sd <- 0
  }else{
    if(!identical(additive$method, "linear_normal") || !isTRUE(additive$exact) ||
       !is.finite(additive$provenance$sd) || additive$provenance$sd <= 0){
      return(NULL)
    }
    additive_mean <- additive$provenance$mean
    additive_sd <- additive$provenance$sd
  }
  list(
    additive_mean = additive_mean,
    additive_sd = additive_sd,
    product_mean = multiplied$provenance$mean,
    product_sd = multiplied$provenance$sd,
    multiplier = multiplier,
    bounds = bounds,
    sources = list(additive = names(additive_groups),
                   multiplied = names(product_groups), multiplier = product$multiplier)
  )
}

.prior_conditional_normal_ordinate <- function(spec, value, n_grid){

  # A pure scale mixture (a_s = 0) at its offset a_m is classified from the
  # multiplier's declared behavior at zero.
  if(spec$additive_sd == 0 && value == spec$additive_mean){
    return(.prior_conditional_normal_offset_ordinate(spec, value, n_grid))
  }
  # The integral runs over the other term's (the multiplier's) declared
  # support, split at breakpoints so that no piece is dominated by mass that
  # its initial quadrature rule cannot see (a narrow Gaussian peak, or a
  # concentrated other term far from zero). Each piece is an independent
  # integral with the full evaluation budget and its own diagnostics; the
  # ordinate is their sum, and the acceptance criterion applies to the total.
  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  integrand <- function(multiplier){
    conditional <- .prior_conditional_normal_moments(spec, multiplier)
    out <- exp(stats::dnorm(value, conditional$mean, conditional$sd, log = TRUE) +
                 multiplier_lpdf(multiplier))
    # away from the offset, a pure scale mixture has no Gaussian mass at
    # the value where the multiplier is zero (the integrand's limit)
    if(spec$additive_sd == 0){
      out[conditional$sd == 0] <- 0
    }
    out
  }
  integral <- .prior_conditional_normal_quadrature(
    integrand, .prior_conditional_normal_breakpoints(spec, value), n_grid,
    zero_message = "zero ordinate for a structurally positive density"
  )
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
      structural_regularity = if(spec$additive_sd > 0){
        "positive_variance_gaussian_convolution"
      }else{
        "scale_mixture_away_from_offset"
      },
      integration = integral$integration
    )
  )
}

# Ordinate of a pure scale mixture X = a_m + b * s (b ~ N(b_m, b_s)) at its
# offset a_m. Given s, the normal density at a_m is phi(b_m / b_s) / (b_s |s|),
# so f(a_m) = phi(b_m / b_s) / b_s * E[1 / |s|]. The expectation is finite
# exactly when the multiplier's density vanishes at zero (a density that
# vanishes like |s|^p, p > 0, at zero, or zero outside its support); a positive
# or infinite multiplier density at zero makes it infinite.
.prior_conditional_normal_offset_ordinate <- function(spec, value, n_grid){

  multiplier_zero <- .prior_density_ordinate_primitive(spec$multiplier, 0)
  multiplier_behavior <- .prior_density_ordinate_continuous_behavior(multiplier_zero)
  provenance <- list(
    kind = "conditional_normal_mixture",
    additive = c(mean = spec$additive_mean, sd = spec$additive_sd),
    multiplied = c(mean = spec$product_mean, sd = spec$product_sd),
    multiplier = .prior_density_ordinate_prior_provenance(spec$multiplier),
    independent_sources = spec$sources,
    structural_regularity = "scale_mixture_offset",
    multiplier_at_zero = multiplier_behavior
  )
  if(multiplier_behavior %in% c("regular", "infinite")){
    provenance$kind <- "product_singularity"
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "infinite",
      log_density = Inf,
      exact       = TRUE,
      method      = "conditional_normal_mixture",
      reason      = paste0(
        "A normal term multiplied by a scale whose density at zero is positive ",
        "or infinite has an infinite density at the requested value."
      ),
      provenance  = provenance
    ))
  }
  if(!identical(multiplier_behavior, "zero")){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      provenance  = provenance
    ))
  }
  moment <- .prior_density_inverse_moment(spec$multiplier, n_grid)
  provenance$inverse_moment <- moment[c("value", "method")]
  if(!is.null(moment$integration)){
    provenance$integration <- moment$integration
  }
  .prior_density_ordinate_result(
    value       = value,
    behavior    = "regular",
    log_density = stats::dnorm(spec$product_mean / spec$product_sd, log = TRUE) -
      log(spec$product_sd) + log(moment$value),
    exact       = TRUE,
    method      = "conditional_normal_mixture",
    provenance  = provenance
  )
}

# E[1 / |s|] of a simple continuous prior whose density vanishes at zero:
# closed forms for the untruncated gamma (shape > 1), inverse-gamma, lognormal
# and beta (alpha > 1) families, otherwise the quadrature of f(s) / |s| over
# the support, split at zero and the prior's quantiles, with the acceptance
# criterion of the conditional-normal quadrature (value NA when rejected).
.prior_density_inverse_moment <- function(prior, n_grid){

  family <- prior$distribution
  parameters <- prior$parameters
  lower <- prior$truncation$lower
  upper <- prior$truncation$upper
  closed <- if(identical(family, "gamma") && lower == 0 && upper == Inf &&
               parameters$shape > 1){
    parameters$rate / (parameters$shape - 1)
  }else if(identical(family, "invgamma") && lower == 0 && upper == Inf){
    parameters$shape / parameters$scale
  }else if(identical(family, "lognormal") && lower == 0 && upper == Inf){
    exp(-parameters$meanlog + parameters$sdlog^2 / 2)
  }else if(identical(family, "beta") && lower == 0 && upper == 1 &&
           parameters$alpha > 1){
    (parameters$alpha + parameters$beta - 1) / (parameters$alpha - 1)
  }else{
    NULL
  }
  if(!is.null(closed)){
    return(list(value = closed, method = "closed_form", integration = NULL))
  }

  spec <- list(additive_mean = 0, additive_sd = 1, product_mean = 0,
               product_sd = 0, multiplier = prior, bounds = c(lower, upper))
  points <- .prior_conditional_normal_breakpoints(spec, 0, extra = 0)
  prior_lpdf <- .prior_simple_lpdf_evaluator(prior)
  integral <- .prior_conditional_normal_quadrature(
    function(s){
      out <- exp(prior_lpdf(s) - log(abs(s)))
      out[s == 0] <- 0
      out
    },
    points, n_grid,
    zero_message = "zero inverse moment of a structurally positive density"
  )
  list(value = integral$value, method = "quadrature", integration = integral$integration)
}

# Region probability P(X in region) of the same conditional-normal mixture:
# the 1-D integral over the multiplier (the other term of a Gaussian
# convolution) s of its density times the Gaussian probability of the region's
# disjoint intervals, int f(s) sum_j [Phi((u_j - mu(s)) / sd(s)) -
# Phi((l_j - mu(s)) / sd(s))] ds (infinite bounds allowed). It uses the
# ordinate's quadrature: the same breakpoints, with a Gaussian-peak window
# around the location s* of every finite interval endpoint under the
# ordinate's peak guard (so always for Gaussian convolutions; the Gaussian
# interval probability changes over about one local SD there), the full budget
# per piece, and the acceptance criterion on the total; a region with a
# positive Gaussian variance never has zero probability, so an exactly zero
# total is rejected as for the ordinate (also when the probability underflows,
# or when it comes only from the pieces evaluated as c = 0 below).
# Unlike the ordinate, whose Gaussian factor vanishes away from its peak, the
# integrand on a piece that the region covers entirely is the other term's
# density itself, which QUADPACK integrates poorly in heavy tails. For a
# Gaussian convolution X = G + w T, G ~ N(m, s), a piece [a, b] that lies
# entirely outside every window t*_e +- 10 s / |w| of the finite region
# endpoints e is therefore evaluated exactly as c * (F(b) - F(a)) with T's
# declared distribution function F: the mean m + w t of G is then at least
# 10 SDs from every region bound, so the Gaussian region probability is
# c = 1 (mean inside the region) or c = 0 (outside) to within
# P(|G - m| > 10 s) = 2 * Phi(-10) ~= 1.5e-23, and the piece's absolute error
# is at most 1.5e-23 * (F(b) - F(a)). The window bounds are breakpoints, so
# only pieces inside the windows use QUADPACK. Known limitation: conditional-
# normal scale mixtures (b_s > 0) integrate every piece; for heavy-tailed
# multipliers (e.g. half-Cauchy or inverse-gamma(1)), QUADPACK flags the
# infinite end piece of a region with an infinite bound ("roundoff error",
# "probably divergent") and the probability stops.
.prior_conditional_normal_region <- function(spec, intervals, n_grid){

  lower <- intervals[, 1L]
  upper <- intervals[, 2L]
  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  integrand <- function(multiplier){
    conditional <- .prior_conditional_normal_moments(spec, multiplier)
    probability <- 0
    for(i in seq_along(lower)){
      probability <- probability + .prior_normal_interval_probability(
        lower[i], upper[i], conditional$mean, conditional$sd
      )
    }
    log_density <- multiplier_lpdf(multiplier)
    out <- numeric(length(multiplier))
    positive <- probability > 0
    out[positive] <- exp(log(probability[positive]) + log_density[positive])
    out
  }
  endpoints <- c(lower, upper)
  endpoints <- unique(endpoints[is.finite(endpoints)])
  exact_piece <- NULL
  if(isTRUE(spec$product_sd == 0) && isTRUE(spec$product_mean != 0) &&
     length(endpoints) > 0L){
    # the same arithmetic as the peak windows of the breakpoints
    centre <- (endpoints - spec$additive_mean) / spec$product_mean
    width <- sqrt(spec$additive_sd^2 + (spec$product_sd * centre)^2) / abs(spec$product_mean)
    window_lower <- centre + (-10) * width
    window_upper <- centre + 10 * width
    exact_piece <- function(a, b){
      if(!all(b <= window_lower | a >= window_upper)){
        return(NULL)
      }
      interior <- if(is.finite(a) && is.finite(b)){
        a / 2 + b / 2
      }else if(is.finite(b)){
        b - 1 - abs(b)
      }else if(is.finite(a)){
        a + 1 + abs(a)
      }else{
        0
      }
      mean <- spec$additive_mean + interior * spec$product_mean
      mass <- .prior_region_prior_mass(spec$multiplier, matrix(c(a, b), 1L))
      list(
        value       = if(any(mean > lower & mean < upper)) mass else 0,
        abs.error   = 2 * stats::pnorm(-10) * mass,
        message     = "OK",
        evaluations = 0L
      )
    }
  }
  .prior_conditional_normal_quadrature(
    integrand, .prior_conditional_normal_breakpoints(spec, endpoints), n_grid,
    zero_message = "zero probability for a structurally positive region",
    exact_piece = exact_piece
  )
}

# Scale products X = c + w * L * s of a simple continuous term L (not a
# full-support normal, which is the conditional-normal route) and a simple
# continuous multiplier s: the ordered-level route with a non-normal total
# (L = total, s = its Beta allocation share) and 'multiply_by' products of a
# non-normal coefficient prior. Away from the offset c the density is the 1-D
# integral f(x) = int f_s(s) f_L((x - c) / (w s)) / |w s| ds over the
# multiplier's support, evaluated by the conditional-normal quadrature (split
# at the multiplier's bounds and quantiles, at its zero, and at the images
# s = (x - c) / (w q) of L's quantiles and finite bounds q); the offset is
# classified from the declared behaviors of L and s at zero.
.prior_scale_product_spec <- function(offset, scale, factor, multiplier, sources){

  list(
    offset     = offset,
    scale      = scale,
    factor     = factor,
    multiplier = multiplier,
    bounds     = unlist(multiplier$truncation[c("lower", "upper")], use.names = FALSE),
    sources    = sources
  )
}

# Quadrature settings of a scale product: the multiplier's breakpoints, no
# Gaussian peak window.
.prior_scale_product_breakpoints <- function(spec, distances){

  factor_bounds <- unlist(spec$factor$truncation[c("lower", "upper")], use.names = FALSE)
  quantiles <- tryCatch(
    suppressWarnings(as.numeric(quant(
      spec$factor, c(1e-6, 1e-3, .02, .25, .5, .75, .98, 1 - 1e-3, 1 - 1e-6)
    ))),
    error = function(e) numeric()
  )
  targets <- c(quantiles, factor_bounds)
  targets <- targets[is.finite(targets) & targets != 0]
  images <- as.vector(outer(distances, targets, `/`))
  .prior_conditional_normal_breakpoints(
    list(additive_mean = 0, additive_sd = 1, product_mean = 0, product_sd = 0,
         multiplier = spec$multiplier, bounds = spec$bounds),
    value = 0,
    extra = c(0, images[is.finite(images)])
  )
}

# Closed interval containing the support of c + w * L * s.
.prior_scale_product_hull <- function(spec){

  factor_bounds <- unlist(spec$factor$truncation[c("lower", "upper")], use.names = FALSE)
  products <- as.vector(outer(factor_bounds, spec$bounds, function(a, b){
    ifelse(a == 0 | b == 0, 0, a * b)
  }))
  spec$offset + sort(spec$scale * range(products))
}

.prior_scale_product_provenance <- function(spec){

  list(
    kind                = "scale_mixture",
    offset              = spec$offset,
    scale               = spec$scale,
    factor              = .prior_density_ordinate_prior_provenance(spec$factor),
    multiplier          = .prior_density_ordinate_prior_provenance(spec$multiplier),
    independent_sources = spec$sources
  )
}

.prior_scale_product_ordinate <- function(spec, value, n_grid){

  provenance <- .prior_scale_product_provenance(spec)
  hull <- .prior_scale_product_hull(spec)
  if(value < hull[1L] || value > hull[2L]){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "zero",
      log_density = -Inf,
      exact       = TRUE,
      method      = "scale_mixture",
      reason      = "The requested value is outside the prior support.",
      provenance  = provenance
    ))
  }
  if(value == spec$offset){
    return(.prior_scale_product_offset_ordinate(spec, value, n_grid, provenance))
  }
  if(value == hull[1L] || value == hull[2L]){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "The limit at a finite bound of a product support is not structurally classified.",
      provenance  = provenance
    ))
  }

  distance <- (value - spec$offset) / spec$scale
  factor_lpdf <- .prior_simple_lpdf_evaluator(spec$factor)
  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  integrand <- function(multiplier){
    out <- exp(factor_lpdf(distance / multiplier) +
                 multiplier_lpdf(multiplier) - log(abs(spec$scale * multiplier)))
    out[multiplier == 0] <- 0
    out
  }
  integral <- .prior_conditional_normal_quadrature(
    integrand, .prior_scale_product_breakpoints(spec, distance), n_grid,
    zero_message = "zero ordinate for a structurally positive density",
    kind = "scale_mixture"
  )
  provenance$structural_regularity <- "scale_mixture_inside_support"
  provenance$integration <- integral$integration
  .prior_density_ordinate_result(
    value = value, behavior = "regular",
    log_density = log(integral$value), exact = TRUE,
    method = "scale_mixture", provenance = provenance
  )
}

# The offset c of c + w * L * s: with f_L and f_s the declared densities at
# zero (one-sided limits at a support bound), the density at c is infinite
# when either is infinite or both are positive, f_L(0) E[1 / |s|] / |w| when
# only f_s vanishes there, f_s(0) E[1 / |L|] / |w| when only f_L vanishes
# there, and zero when both vanish. A term bounded at zero combined with a
# two-sided other term makes the offset a density jump, which is not
# classified.
.prior_scale_product_offset_ordinate <- function(spec, value, n_grid, provenance){

  factor_zero <- .prior_density_ordinate_primitive(spec$factor, 0)
  multiplier_zero <- .prior_density_ordinate_primitive(spec$multiplier, 0)
  behaviors <- c(
    factor     = .prior_density_ordinate_continuous_behavior(factor_zero),
    multiplier = .prior_density_ordinate_continuous_behavior(multiplier_zero)
  )
  provenance$structural_regularity <- "scale_mixture_offset"
  provenance$behaviors_at_zero <- behaviors
  unknown <- function(reason){
    .prior_density_ordinate_result(
      value = value, behavior = "unknown", log_density = NA_real_,
      exact = FALSE, method = "unsupported_provenance", reason = reason,
      provenance = provenance
    )
  }
  if(any(!behaviors %in% c("regular", "zero", "infinite"))){
    return(unknown(NULL))
  }
  one_sided <- function(prior){
    bounds <- unlist(prior$truncation[c("lower", "upper")], use.names = FALSE)
    any(bounds == 0)
  }
  two_sided <- function(prior){
    bounds <- unlist(prior$truncation[c("lower", "upper")], use.names = FALSE)
    bounds[1L] < 0 && bounds[2L] > 0
  }
  if((one_sided(spec$factor) && two_sided(spec$multiplier) && behaviors[["factor"]] != "zero") ||
     (one_sided(spec$multiplier) && two_sided(spec$factor) && behaviors[["multiplier"]] != "zero")){
    return(unknown("The density of the product jumps at the requested value."))
  }
  if(any(behaviors == "infinite") || all(behaviors == "regular")){
    return(.prior_density_ordinate_result(
      value = value, behavior = "infinite", log_density = Inf, exact = TRUE,
      method = "scale_mixture",
      reason = paste0(
        "A product of independent continuous terms whose densities at zero ",
        "are both positive, or one of them infinite, has an infinite density ",
        "at the requested value."
      ),
      provenance = provenance
    ))
  }
  if(all(behaviors == "zero")){
    return(.prior_density_ordinate_result(
      value = value, behavior = "zero", log_density = -Inf, exact = TRUE,
      method = "scale_mixture", provenance = provenance
    ))
  }
  if(behaviors[["factor"]] == "regular"){
    density_zero <- factor_zero$log_density
    moment <- .prior_density_inverse_moment(spec$multiplier, n_grid)
  }else{
    density_zero <- multiplier_zero$log_density
    moment <- .prior_density_inverse_moment(spec$factor, n_grid)
  }
  provenance$inverse_moment <- moment[c("value", "method")]
  if(!is.null(moment$integration)){
    provenance$integration <- moment$integration
  }
  .prior_density_ordinate_result(
    value = value, behavior = "regular",
    log_density = density_zero + log(moment$value) - log(abs(spec$scale)),
    exact = TRUE, method = "scale_mixture", provenance = provenance
  )
}

# P(lower < L < upper) for a simple continuous prior, vectorized over the
# bounds; upper-tail probabilities are used above the median.
.prior_scalar_interval_probability <- function(prior, lower, upper){

  lower_cdf <- cdf(prior, lower)
  out <- cdf(prior, upper) - lower_cdf
  upper_tail <- !is.na(lower_cdf) & lower_cdf > .5
  if(any(upper_tail)){
    out[upper_tail] <- ccdf(prior, lower[upper_tail]) - ccdf(prior, upper[upper_tail])
  }
  pmax(out, 0)
}

# Region probability of c + w * L * s: the 1-D integral over s of its density
# times P(c + w L s in region), with the ordinate's breakpoints at the images
# of every finite region endpoint.
.prior_scale_product_region <- function(spec, intervals, n_grid){

  lower <- intervals[, 1L]
  upper <- intervals[, 2L]
  offset_inside <- as.numeric(any(spec$offset > lower & spec$offset < upper))
  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  integrand <- function(multiplier){
    probability <- numeric(length(multiplier))
    zero <- multiplier == 0
    s <- multiplier[!zero]
    for(i in seq_along(lower)){
      a <- (lower[i] - spec$offset) / (spec$scale * s)
      b <- (upper[i] - spec$offset) / (spec$scale * s)
      flip <- spec$scale * s < 0
      probability[!zero] <- probability[!zero] + .prior_scalar_interval_probability(
        spec$factor, ifelse(flip, b, a), ifelse(flip, a, b)
      )
    }
    probability[zero] <- offset_inside
    out <- numeric(length(multiplier))
    positive <- probability > 0
    out[positive] <- exp(log(probability[positive]) + multiplier_lpdf(multiplier[positive]))
    out
  }
  endpoints <- c(lower, upper)
  endpoints <- unique(endpoints[is.finite(endpoints)])
  .prior_conditional_normal_quadrature(
    integrand,
    .prior_scale_product_breakpoints(spec, (endpoints - spec$offset) / spec$scale),
    n_grid,
    zero_message = "zero probability for a structurally positive region",
    kind = "scale_mixture"
  )
}

# Mean and SD of the conditional normal N(a_m + b_m s, sqrt(a_s^2 + b_s^2 s^2))
# at multiplier values s, with the SD computed without overflow.
.prior_conditional_normal_moments <- function(spec, multiplier){

  product_sd <- abs(multiplier) * spec$product_sd
  if(spec$additive_sd == 0){
    # a pure scale mixture: the conditional SD is that of the multiplied term
    return(list(mean = spec$additive_mean + multiplier * spec$product_mean,
                sd   = product_sd))
  }
  scale <- pmax(spec$additive_sd, product_sd)
  list(
    mean = spec$additive_mean + multiplier * spec$product_mean,
    sd   = scale * sqrt((spec$additive_sd / scale)^2 + (product_sd / scale)^2)
  )
}

# P(lower < Z < upper) for Z ~ N(mean, sd), vectorized over 'mean' and 'sd';
# upper-tail probabilities are used when the interval lies above the mean, so
# small probabilities in either tail keep their relative precision. A zero SD
# (the multiplier's zero of a pure scale mixture) is the point mass at 'mean'.
.prior_normal_interval_probability <- function(lower, upper, mean, sd){

  z_lower <- (lower - mean) / sd
  z_upper <- (upper - mean) / sd
  out <- stats::pnorm(z_upper) - stats::pnorm(z_lower)
  upper_tail <- !is.na(z_lower) & z_lower > 0
  out[upper_tail] <-
    stats::pnorm(z_lower[upper_tail], lower.tail = FALSE) -
    stats::pnorm(z_upper[upper_tail], lower.tail = FALSE)
  degenerate <- !is.na(sd) & sd == 0
  if(any(degenerate)){
    mean <- rep_len(mean, length(out))
    out[degenerate] <- as.numeric(mean[degenerate] > lower & mean[degenerate] < upper)
  }
  pmax(out, 0)
}

# Sum of the budgeted QUADPACK pieces between consecutive 'points'. Each piece
# is an independent integral with the full budget 'n_grid' and its own
# diagnostics; the value and absolute error are summed, and the acceptance
# criterion (all pieces converged, abs. error <= 1e-12 + 1e-4 * value, value
# positive) applies to the total. A rejected total has value NA. An optional
# 'exact_piece(lower, upper)' returns a piece evaluated without quadrature
# (value, abs.error, message, evaluations) or NULL for a QUADPACK piece.
.prior_conditional_normal_quadrature <- function(integrand, points, n_grid,
                                                 zero_message,
                                                 exact_piece = NULL,
                                                 kind = "conditional_normal_mixture"){

  tolerance <- .prior_linear_density_refinement_tolerance()
  n_pieces <- length(points) - 1L
  exact <- rep(FALSE, n_pieces)
  pieces <- lapply(seq_len(n_pieces), function(i){
    if(!is.null(exact_piece)){
      piece <- exact_piece(points[i], points[i + 1L])
      if(!is.null(piece)){
        exact[i] <<- TRUE
        return(piece)
      }
    }
    .prior_conditional_normal_piece(
      integrand, points[i], points[i + 1L], n_grid,
      relative = tolerance$relative, absolute = tolerance$absolute / n_pieces
    )
  })
  piece_values <- vapply(pieces, `[[`, numeric(1), "value")
  piece_errors <- vapply(pieces, `[[`, numeric(1), "abs.error")
  piece_messages <- vapply(pieces, `[[`, character(1), "message")
  integral <- list(
    value     = sum(piece_values),
    abs.error = sum(piece_errors),
    message   = if(all(piece_messages == "OK")){
      "OK"
    }else{
      paste(unique(piece_messages[piece_messages != "OK"]), collapse = "; ")
    }
  )
  if(identical(integral$message, "OK") && isTRUE(integral$value == 0)){
    integral$message <- zero_message
  }
  bound <- tolerance$absolute + tolerance$relative * abs(integral$value)
  accepted <- identical(integral$message, "OK") && is.finite(integral$value) &&
    integral$value > 0 && is.finite(integral$abs.error) && integral$abs.error <= bound
  if(!isTRUE(accepted)){
    integral$value <- NA_real_
  }
  integration <- list(
    kind = kind, exact = FALSE,
    absolute_error = integral$abs.error, error_bound = bound,
    evaluations = sum(vapply(pieces, `[[`, integer(1), "evaluations")),
    budget = n_grid, converged = isTRUE(accepted), message = integral$message,
    breakpoints = points,
    piece_evaluations = vapply(pieces, `[[`, integer(1), "evaluations"),
    piece_absolute_errors = piece_errors
  )
  if(!is.null(exact_piece)){
    integration$exact_pieces <- exact
  }
  list(value = integral$value, integration = integration)
}

# Breakpoints of the conditional-normal integral over the multiplier (or the
# other term of a Gaussian convolution) s:
# * the support bounds of s and quantiles of its declared prior;
# * the location peak of the conditional normal N(value; a_m + b_m s,
#   sqrt(a_s^2 + b_s^2 s^2)) at s* = (value - a_m) / b_m and s* +- k w
#   (k = 1, 3, 10), w = sqrt(a_s^2 + b_s^2 s*^2) / |b_m|, when the multiplied
#   SD is at most half of the absolute value of its mean (b_s <= |b_m| / 2).
#   The standardized distance |value - a_m - b_m s| / sqrt(a_s^2 + b_s^2 s^2)
#   then still grows away from s* (to |b_m| / b_s >= 2 far from it), so the
#   integrand is concentrated around s*; with a tighter guard (1/10), pieces
#   missed such peaks just beyond it. For b_s much larger than |b_m| the
#   window is not a peak (its points lie far outside the mass and made the
#   pieces miss mass elsewhere), and no peak points are used.
#   A Gaussian convolution (b_s = 0, b_m = w) always has its Gaussian peak
#   u* = (value - m) / w with width s / |w|.
# Only points strictly inside the open support with a finite density are
# kept. Next to a finite bound where the density is infinite, QUADPACK
# integrates a piece ending very close to the bound as if the singularity
# were at that end, counting the mass between the bound and the end twice:
# there the extreme (1e-6) quantile is not used, and between the bound and the
# quartile on that side a point is kept only if its distance to the bound is
# at least 1e-3 of the next kept point's distance. Gaussian-peak points at
# least one local SD w from the bound are exempt: they resolve a real peak,
# which the piece from the bound would otherwise miss; peak points closer
# than w to the bound (rounding residues of s* - k w) are still dropped.
# Points closer than min(1e-9 * max(1, |endpoints|), w / 2), but at least
# 16 * eps * max(1, |endpoints|), are merged (the bounds are kept): points of
# different sources that nearly coincide (e.g. a peak point a rounding error
# from a bound) would leave a piece only a few ulps wide, which QUADPACK
# cannot resolve, while the peak points, at least w apart, are never merged.
# Known limitations: pure scale mixtures (b_m = 0) with heavy-tailed
# multipliers can stop as non-convergent (mostly with a small multiplied SD);
# a flagged piece stops the ordinate even when its value is far below the
# absolute tolerance (e.g. a far tail piece reported as "probably divergent");
# an exactly zero integral is rejected (a missed peak looks the same), which
# also stops a mixture with such a component; very narrow Gaussian peaks can
# stop as non-convergent where QUADPACK reaches the floating-point resolution
# (w = 1e-10 next to a singular bound at 1, or w of a few hundred ulps); and
# the scale peak s ~ |value - a_m| / b_s of a value far in the multiplier's
# heavy tail (beyond its extreme quantiles) is not a breakpoint, also with
# b_m != 0 when b_s is not small against |b_m|, so its mass can be missed
# without a convergence failure. Such ordinates are small at ordinary scales,
# but the missed fraction does not depend on the units of the value.
# Region probabilities pass several values (the finite region endpoints): each
# gets its own peak window, each peak point keeps its own local SD for the
# exemption next to a singular bound, and the merge width is capped by the
# smallest local SD, so no peak point of any window is merged.
# A pure scale mixture (a_s = 0) also splits at the multiplier's zero and at
# |value - a_m| / b_s * (1/10, 1, 10) on both sides of it: the integrand rises
# from zero at s = 0 to the scale peak near |s| = |value - a_m| / b_s (exactly
# there for b_m = 0) and then decays like f(s) / |s|. 'extra' points (e.g. the
# images of another term's quantiles) are added like quantiles.
.prior_conditional_normal_breakpoints <- function(spec, value, extra = numeric()){

  lower <- spec$bounds[1L]
  upper <- spec$bounds[2L]
  inner <- numeric()
  is_peak <- logical()
  widths <- numeric()
  if(isTRUE(spec$product_mean != 0) &&
     isTRUE(spec$product_sd <= abs(spec$product_mean) / 2)){
    for(peak_value in value){
      centre <- (peak_value - spec$additive_mean) / spec$product_mean
      width <- sqrt(spec$additive_sd^2 + (spec$product_sd * centre)^2) / abs(spec$product_mean)
      inner <- c(inner, centre, centre + as.vector(outer(c(-1, 1), c(1, 3, 10))) * width)
      is_peak <- c(is_peak, rep(TRUE, 7L))
      widths <- c(widths, rep(if(is.finite(width)) width else Inf, 7L))
    }
  }
  if(isTRUE(spec$additive_sd == 0) && isTRUE(spec$product_sd > 0)){
    scale_points <- unlist(lapply(value, function(peak_value){
      distance <- abs(peak_value - spec$additive_mean) / spec$product_sd
      if(!is.finite(distance) || distance <= 0){
        return(numeric())
      }
      as.vector(outer(c(-1, 1), distance * c(.1, 1, 10)))
    }))
    extra <- c(extra, 0, scale_points)
  }
  if(length(extra) > 0L){
    inner <- c(inner, extra)
    is_peak <- c(is_peak, rep(FALSE, length(extra)))
    widths <- c(widths, rep(Inf, length(extra)))
  }
  peak_width <- min(c(Inf, widths))
  singular <- vapply(c(lower, upper), function(bound){
    is.finite(bound) && isTRUE(is.infinite(suppressWarnings(exp(lpdf(spec$multiplier, bound)))))
  }, logical(1))
  probabilities <- c(if(!singular[1L]) 1e-6, 1e-3, .02, .25, .5, .75, .98, 1 - 1e-3,
                     if(!singular[2L]) 1 - 1e-6)
  quantiles <- tryCatch(
    suppressWarnings(quant(spec$multiplier, probabilities)),
    error = function(e) numeric()
  )
  inner <- c(inner, as.numeric(quantiles))
  is_peak <- c(is_peak, rep(FALSE, length(quantiles)))
  widths <- c(widths, rep(Inf, length(quantiles)))
  inside <- is.finite(inner) & inner > lower & inner < upper
  inner <- inner[inside]
  is_peak <- is_peak[inside]
  widths <- widths[inside]
  if(length(inner) > 0L){
    density <- suppressWarnings(exp(lpdf(spec$multiplier, inner)))
    inner <- inner[is.finite(density)]
    is_peak <- is_peak[is.finite(density)]
    widths <- widths[is.finite(density)]
  }
  # sorted and unique; a peak point equal to a quantile stays a peak point
  sorted <- order(inner)
  inner <- inner[sorted]
  is_peak <- is_peak[sorted]
  widths <- widths[sorted]
  first <- !duplicated(inner)
  inner <- inner[first]
  is_peak <- is_peak[first]
  widths <- widths[first]

  if(any(singular)){
    quartiles <- tryCatch(
      suppressWarnings(as.numeric(quant(spec$multiplier, c(.25, .75)))),
      error = function(e) c(NA_real_, NA_real_)
    )
    for(side in which(singular)){
      bound <- c(lower, upper)[side]
      start <- abs(quartiles[side] - bound)
      if(!isTRUE(is.finite(start))){
        next
      }
      distance <- abs(inner - bound)
      keep <- rep(TRUE, length(inner))
      last <- start
      for(i in order(distance, decreasing = TRUE)){
        if(distance[i] >= start){
          next
        }
        exempt <- is_peak[i] && distance[i] >= widths[i]
        if(exempt || distance[i] >= 1e-3 * last){
          last <- distance[i]
        }else{
          keep[i] <- FALSE
        }
      }
      inner <- inner[keep]
      is_peak <- is_peak[keep]
      widths <- widths[keep]
    }
  }

  minimum_width <- function(a, b){
    scale <- max(1, abs(c(a, b))[is.finite(c(a, b))])
    max(16 * .Machine$double.eps * scale, min(1e-9 * scale, peak_width / 2))
  }
  points <- lower
  for(point in inner){
    if(point - points[length(points)] >= minimum_width(points[length(points)], point)){
      points <- c(points, point)
    }
  }
  if(length(points) > 1L && upper - points[length(points)] < minimum_width(points[length(points)], upper)){
    points <- points[-length(points)]
  }
  c(points, upper)
}

# One budgeted QUADPACK piece. Finite pieces use 21-point rules; infinite ends
# use 15-point rules, paired for two-sided infinite pieces. With k evaluations
# per interval, m intervals use k * (2 * m - 1) evaluations; 'max_intervals'
# is the largest m within the budget. QUADPACK reports "maximum number of
# subdivisions reached" whenever the interval count equals 'subdivisions', even
# for a converged result, so it gets one more interval and the evaluation cap
# enforces the budget.
.prior_conditional_normal_piece <- function(integrand, lower, upper, n_grid,
                                            relative, absolute){

  initial_evaluations <- .prior_conditional_normal_initial_evaluations(c(lower, upper))
  max_intervals <- floor((n_grid + initial_evaluations) / (2 * initial_evaluations))
  evaluations <- 0L
  budgeted <- function(x){
    if(evaluations + length(x) > n_grid){
      stop("the integration evaluation budget was exhausted", call. = FALSE)
    }
    evaluations <<- evaluations + length(x)
    integrand(x)
  }
  integral <- if(!is.finite(max_intervals) || max_intervals < 1){
    list(value = NA_real_, abs.error = NA_real_,
         message = paste0("fewer than ", initial_evaluations,
                          " integration evaluations are available"))
  }else tryCatch(
    stats::integrate(budgeted, lower, upper,
                     subdivisions = max_intervals + 1L,
                     rel.tol = relative, abs.tol = absolute, stop.on.error = FALSE),
    error = function(e){
      list(value = NA_real_, abs.error = NA_real_, message = conditionMessage(e))
    }
  )
  list(
    value       = as.numeric(integral$value),
    abs.error   = as.numeric(integral$abs.error),
    message     = as.character(integral$message),
    evaluations = evaluations
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

.prior_linear_density_recipe_width <- function(context, tail_prob){

  # Width of the source-grid range that a recorded prior density would span
  # with 'tail_prob' omitted per source, from the prior ranges alone; NULL when
  # it cannot be determined.
  combination_range <- function(prior_list, weights, source_transforms){
    weights <- weights[is.finite(weights) & weights != 0]
    if(length(weights) == 0L){
      return(c(0, 0))
    }
    groups <- .prior_linear_weight_groups(prior_list, weights)
    if(is.null(source_transforms)){
      source_transforms <- stats::setNames(rep(NA_character_, length(weights)), names(weights))
    }
    ranges <- lapply(groups, .prior_linear_group_range, tail_prob = tail_prob,
                     source_transforms = source_transforms)
    # A 'multiply_by' scale widens the range of the coefficients it scales.
    multipliers <- unique(unlist(lapply(groups, function(group){
      multiply_by <- attr(group$prior, "multiply_by", exact = TRUE)
      if(is.character(multiply_by) && length(multiply_by) == 1L) multiply_by else NULL
    })))
    for(multiplier in setdiff(multipliers, names(groups))){
      ranges[[length(ranges) + 1L]] <- .prior_linear_group_range(
        list(prior = prior_list[[multiplier]], weights = stats::setNames(1, multiplier), indices = 1L),
        tail_prob = tail_prob
      )
    }
    c(sum(vapply(ranges, `[`, numeric(1), 1L)), sum(vapply(ranges, `[`, numeric(1), 2L)))
  }
  context_range <- function(density_context, weights, source_transforms){
    if(inherits(density_context, "prior_density_context")){
      standardized <- .prior_density_context_standardized_weights(density_context, weights)
      return(combination_range(
        density_context$prior_list,
        standardized,
        if(is.null(source_transforms)) NULL else source_transforms[names(standardized)]
      ))
    }
    if(inherits(density_context, "prior_density_model_mixture_context")){
      models <- which(density_context$model_weights > 0)
      return(range(unlist(lapply(models, function(model_i){
        combination_range(
          .prior_density_model_prior_list(density_context$prior_list, model_i),
          weights, source_transforms
        )
      }))))
    }
    if(inherits(density_context, "prior_density_conditional_context")){
      models <- which(density_context$model_weights > 0)
      return(range(unlist(lapply(density_context$prior_lists[models], function(prior_list){
        if(!is.null(density_context$formula_scale) && length(density_context$formula_scale) > 0L){
          return(context_range(
            .prior_density_context(prior_list, density_context$column_names,
                                   density_context$formula_scale),
            weights, source_transforms
          ))
        }
        combination_range(prior_list, weights, source_transforms)
      }))))
    }
    NULL
  }

  arguments <- context$arguments
  out <- tryCatch({
    bounds <- if(identical(context$kind, "linear_combination")){
      combination_range(arguments$prior_list, arguments$weights, arguments$source_transforms)
    }else if(identical(context$kind, "density_context_rows") && !is.null(dim(arguments$weights))){
      range(unlist(lapply(seq_len(nrow(arguments$weights)), function(row_i){
        context_range(arguments$context, arguments$weights[row_i, ], arguments$source_transforms)
      })))
    }else{
      context_range(arguments$context, arguments$weights, arguments$source_transforms)
    }
    diff(bounds)
  }, error = function(e) NULL)
  if(length(out) != 1L || !is.finite(out) || out <= 0){
    return(NULL)
  }
  out
}

.prior_linear_density_refined_tail <- function(context, tail_prob){

  # The omitted tail mass always shrinks, so its truncation bias stays visible
  # to the convergence check; the largest reduction whose range at most
  # doubles is preferred, otherwise the smallest one.
  fallback <- max(tail_prob / 10, 1e-12)
  width <- .prior_linear_density_recipe_width(context, tail_prob)
  if(is.null(width)){
    return(fallback)
  }
  for(factor in c(1000, 100)){
    candidate <- max(tail_prob / factor, 1e-12)
    if(candidate >= fallback){
      break
    }
    candidate_width <- .prior_linear_density_recipe_width(context, candidate)
    if(!is.null(candidate_width) && candidate_width <= 2 * width){
      return(candidate)
    }
  }
  fallback
}

.prior_linear_density_refinement <- function(x){

  context <- attr(x, "adaptive_evaluation", exact = TRUE)
  if(is.null(context) ||
     !context$kind %in%
       c("linear_combination", "density_context", "density_context_rows")){
    return(NULL)
  }

  # Each refinement strictly halves the spacing of the source grid and omits
  # less tail probability by the largest factor (1000, 100, 10) whose range at
  # most doubles, and by 10 when none does. When the halved spacing needs more
  # knots than the grid limit, no refinement is available and callers report
  # non-convergence.
  spacing <- attr(x, "grid_resolution", exact = TRUE)[["spacing"]]
  if(is.null(spacing) || !is.finite(spacing) || spacing <= 0){
    spacing <- .prior_linear_density_dx(x)
  }
  if(!is.finite(spacing) || spacing <= 0){
    return(NULL)
  }
  arguments <- context$arguments
  tail_prob <- if(identical(context$kind, "linear_combination")){
    arguments$tail_prob
  }else{
    arguments$context$tail_prob
  }
  tail_prob <- .prior_linear_density_refined_tail(context, tail_prob)
  if(identical(context$kind, "linear_combination")){
    arguments$tail_prob    <- tail_prob
    arguments$grid_spacing <- spacing / 2
  }else{
    arguments$context$tail_prob    <- tail_prob
    arguments$context$grid_spacing <- spacing / 2
  }
  refined_arguments <- arguments
  refined_arguments$.record_evaluation <- FALSE
  refined <- tryCatch(
    if(identical(context$kind, "linear_combination")){
      do.call(.prior_linear_combination_density, refined_arguments)
    }else if(identical(context$kind, "density_context_rows")){
      do.call(.prior_density_from_context_rows, refined_arguments)
    }else{
      do.call(.prior_density_from_context, refined_arguments)
    },
    BayesTools_prior_grid_limit = function(e) NULL
  )
  if(is.null(refined)){
    return(NULL)
  }
  context$arguments <- arguments
  attr(refined, "adaptive_evaluation") <- context
  attr(refined, "refinement_settings") <- list(
    n_grid = as.integer(attr(refined, "grid_resolution", exact = TRUE)[["n_grid"]]),
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

    # The continuous density is evaluated on its structural route (closed
    # forms, and the ordinate's quadrature at each plotted value); only a
    # combination without a structural route, or a density without recorded
    # provenance, interpolates its numerical grid.
    finite_raw <- is.finite(x_raw)
    y_den      <- rep(NA_real_, length(x_raw))
    route      <- .prior_density_route_from_adaptive(
      attr(dist, "adaptive_evaluation", exact = TRUE)
    )
    evaluator  <- attr(dist, "density_evaluator", exact = TRUE)
    if(any(finite_raw) && !is.null(route) && !identical(route$type, "unknown")){
      y_den[finite_raw] <- .prior_density_route_density(route, x_raw[finite_raw])
    }else if(any(finite_raw) && is.null(route) && is.function(evaluator)){
      evaluated <- evaluator(x_raw[finite_raw])
      if(!is.numeric(evaluated) || length(evaluated) != sum(finite_raw) ||
         anyNA(evaluated)){
        stop("The analytic prior density evaluator returned invalid values.",
             call. = FALSE)
      }
      y_den[finite_raw] <- evaluated * dist$density$mass
    }else if(any(finite_raw)){
      y_den[finite_raw] <- stats::approx(
        dist$density$x,
        dist$density$y,
        xout   = x_raw[finite_raw],
        yleft  = 0,
        yright = 0
      )$y * dist$density$mass
    }

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
      .prior_linear_combination_support_hull(
        .prior_density_model_prior_list(context$prior_list, model_i),
        weights, source_transforms
      )
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

# Height from an exact or regular structural ordinate, or NULL. Quadrature
# ordinates must have converged; exact ordinates (including exact finite
# mixtures) are used directly, since grid refinement cannot converge across a
# density jump.
.prior_linear_density_exact_height <- function(ordinate){

  if(is.null(ordinate)){
    return(NULL)
  }
  continuous <- .prior_density_ordinate_continuous_behavior(ordinate)
  if(identical(continuous, "regular") &&
     .prior_density_ordinate_has_quadrature(ordinate$provenance)){
    integration <- .prior_density_ordinate_integration(ordinate$provenance)
    if(is.na(ordinate$log_density) || !isTRUE(integration$converged)){
      subject <- if(identical(integration$kind, "scale_mixture")){
        "Scale-mixture prior density"
      }else{
        "Conditional-normal prior density"
      }
      stop(subject, " was rejected by diagnostics: integration reported '",
           integration$message, "' with absolute error ", format(integration$absolute_error),
           ". Inspect the prior specification and increase 'n_samples' for marginal inference.",
           call. = FALSE)
    }
    height <- exp(ordinate$log_density)
    attr(height, "numerical_diagnostics") <- integration
    return(height)
  }
  if(isTRUE(ordinate$exact)){
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
  NULL
}

# Height of a finite mixture whose ordinate is not exact (a model or
# conditional mixture, or mixture and spike-and-slab terms of one linear
# combination). Components with an exact or regular ordinate contribute it
# directly; every other component gets its own adaptive grid, so no grid spans
# a density jump between components. The component grids are refined in
# lockstep and the documented criterion is applied to the mixture height with
# the weighted absolute changes of the components, sum_k w_k |dH_k| <=
# 1e-12 + 1e-4 * H, so changes cannot cancel between components and a
# component whose own density is about zero at the value does not block
# convergence. As for single grids, the criterion is a refinement-change
# criterion: grid-based heights of components with singular source densities
# remain approximate within it. Component grids start at the mixture grid's source
# spacing ('grid_spacing'), so they are never coarser than the grid they
# replace. NULL when the density is not such a mixture, is a row mixture, or
# carries an output transformation.
.prior_linear_density_component_height <- function(adaptive, value,
                                                   grid_spacing = NULL){

  if(!is.list(adaptive) || !is.list(adaptive$arguments) ||
     !is.null(adaptive$arguments$output_transformation)){
    return(NULL)
  }
  if(!is.numeric(grid_spacing) || length(grid_spacing) != 1L ||
     !is.finite(grid_spacing) || grid_spacing <= 0){
    grid_spacing <- NULL
  }
  arguments <- adaptive$arguments
  if(identical(adaptive$kind, "linear_combination")){
    if(!is.list(arguments$prior_list) || !is.numeric(arguments$weights) ||
       is.null(names(arguments$weights))){
      return(NULL)
    }
    terms <- .prior_linear_mixture_terms(
      prior_list        = arguments$prior_list,
      weights           = arguments$weights,
      source_transforms = arguments$source_transforms,
      value             = value,
      n_grid            = if(is.null(arguments$n_grid)) .prior_linear_density_default_grid() else arguments$n_grid,
      tail_prob         = if(is.null(arguments$tail_prob)) .prior_linear_density_tail_prob() else arguments$tail_prob,
      grid_spacing      = grid_spacing,
      top_level         = TRUE
    )
    return(.prior_linear_mixture_lockstep_height(terms, value))
  }
  if(!identical(adaptive$kind, "density_context")){
    return(NULL)
  }
  context <- arguments$context
  weights <- arguments$weights
  source_transforms <- arguments$source_transforms

  context_terms <- function(component_context, top_level){
    standardized <- tryCatch(
      .prior_density_context_standardized_weights(component_context, weights),
      error = function(e) NULL
    )
    if(is.null(standardized)){
      return(NULL)
    }
    .prior_linear_mixture_terms(
      prior_list        = component_context$prior_list,
      weights           = standardized,
      source_transforms = if(is.null(source_transforms)) NULL else source_transforms[names(standardized)],
      value             = value,
      n_grid            = component_context$n_grid,
      tail_prob         = component_context$tail_prob,
      grid_spacing      = grid_spacing,
      top_level         = top_level
    )
  }
  if(inherits(context, "prior_density_context")){
    return(.prior_linear_mixture_lockstep_height(context_terms(context, TRUE), value))
  }

  # with a single positive-weight component, that component is the density
  indices <- which(context$model_weights > 0)
  single <- length(indices) < 2L
  if(inherits(context, "prior_density_model_mixture_context")){
    components <- lapply(indices, function(model_i){
      .prior_linear_mixture_terms(
        .prior_density_model_prior_list(context$prior_list, model_i),
        weights, source_transforms, value,
        context$n_grid, context$tail_prob, grid_spacing, single
      )
    })
  }else if(inherits(context, "prior_density_conditional_context")){
    components <- lapply(context$prior_lists[indices], function(prior_list){
      if(!is.null(context$formula_scale) && length(context$formula_scale) > 0L){
        component_context <- .prior_density_context(
          prior_list, context$column_names, context$formula_scale,
          context$n_grid, context$tail_prob
        )
        component_context$grid_spacing <- context$grid_spacing
        return(context_terms(component_context, single))
      }
      .prior_linear_mixture_terms(
        prior_list, weights, source_transforms, value,
        context$n_grid, context$tail_prob, grid_spacing, single
      )
    })
  }else{
    return(NULL)
  }
  .prior_linear_mixture_lockstep_height(
    .prior_linear_weight_mixture_terms(components, context$model_weights[indices]),
    value
  )
}

# Mixture terms of one linear combination: a list of components, each with its
# mixture weight and either its exact/regular height ('exact', 'quadrature')
# or its own density for grid evaluation ('grid'). NULL for a top-level
# combination that is not a mixture (its own grid applies).
.prior_linear_mixture_terms <- function(prior_list, weights, source_transforms,
                                        value, n_grid, tail_prob, grid_spacing,
                                        top_level){

  weights <- weights[weights != 0]
  ordinate <- .prior_density_ordinate_linear_base(
    prior_list, weights, source_transforms, value, n_grid
  )
  exact <- .prior_linear_density_exact_height(ordinate)
  if(!is.null(exact)){
    quadrature <- .prior_density_ordinate_has_quadrature(ordinate$provenance)
    return(list(list(
      weight = 1,
      method = if(quadrature) "quadrature" else "exact",
      height = as.numeric(exact)
    )))
  }

  active <- .prior_linear_active_parameters(prior_list, weights)
  multipliers <- unlist(lapply(prior_list[active], function(prior){
    multiply_by <- attr(prior, "multiply_by", exact = TRUE)
    if(is.character(multiply_by) && length(multiply_by) == 1L) multiply_by else NULL
  }))
  plan <- .prior_density_route_mixture_expansion(prior_list, c(active, multipliers), n_grid)
  if(!is.null(plan)){
    components <- lapply(plan$prior_lists, function(component_priors){
      .prior_linear_mixture_terms(
        component_priors, weights, source_transforms, value,
        n_grid, tail_prob, grid_spacing, FALSE
      )
    })
    return(.prior_linear_weight_mixture_terms(components, plan$probabilities))
  }
  if(isTRUE(top_level)){
    return(NULL)
  }
  .prior_linear_density_check_grid(.prior_density_route_linear(
    prior_list, weights, source_transforms, n_grid
  ))

  # the mixture spacing may need more knots than this component's own grid
  # admits; the component then starts at its own spacing
  build <- function(spacing){
    .prior_linear_combination_density(
      prior_list        = prior_list,
      weights           = weights,
      n_grid            = n_grid,
      tail_prob         = tail_prob,
      source_transforms = source_transforms,
      grid_spacing      = spacing
    )
  }
  density <- tryCatch(build(grid_spacing), BayesTools_prior_grid_limit = function(e) NULL)
  if(is.null(density)){
    density <- build(NULL)
  }
  list(list(weight = 1, method = "grid", density = density))
}

.prior_linear_weight_mixture_terms <- function(components, weights){

  if(length(components) == 0L || any(vapply(components, is.null, logical(1)))){
    return(NULL)
  }
  weights <- weights / sum(weights)
  unlist(lapply(seq_along(components), function(i){
    lapply(components[[i]], function(term){
      term$weight <- term$weight * weights[[i]]
      term
    })
  }), recursive = FALSE)
}

.prior_linear_mixture_lockstep_height <- function(terms, value){

  if(is.null(terms) || length(terms) == 0L){
    return(NULL)
  }
  # grid components outside their support hull contribute exactly zero, and
  # scalar components carry an analytic density evaluator
  for(i in seq_along(terms)){
    if(!identical(terms[[i]]$method, "grid")){
      next
    }
    density <- terms[[i]]$density
    support <- .prior_linear_density_support_hull(
      attr(density, "adaptive_evaluation", exact = TRUE)
    )
    if(!is.null(support) && (value < support[1L] || value > support[2L])){
      terms[[i]] <- list(weight = terms[[i]]$weight, method = "exact", height = 0)
      next
    }
    if(value %in% attr(density, "singular_density_points", exact = TRUE)){
      stop("The prior density at the flagged product ordinate is unavailable from supported deterministic provenance. Inspect the prior specification or use a supported prior-density evaluator.",
           call. = FALSE)
    }
    evaluator <- attr(density, "density_evaluator", exact = TRUE)
    if(is.function(evaluator)){
      terms[[i]] <- list(weight = terms[[i]]$weight, method = "exact",
                         height = as.numeric(evaluator(value)))
    }
  }

  weights <- vapply(terms, `[[`, numeric(1), "weight")
  grid <- vapply(terms, function(term) identical(term$method, "grid"), logical(1))
  fixed_heights <- vapply(terms[!grid], `[[`, numeric(1), "height")
  if(any(is.infinite(fixed_heights) & weights[!grid] > 0)){
    return(Inf)
  }
  fixed <- sum(weights[!grid] * fixed_heights)
  diagnostics <- function(grid_heights){
    heights <- numeric(length(terms))
    heights[!grid] <- fixed_heights
    heights[grid] <- grid_heights
    lapply(seq_along(terms), function(i){
      list(weight = unname(weights[i]), method = terms[[i]]$method,
           height = unname(heights[i]))
    })
  }
  if(!any(grid)){
    height <- fixed
    attr(height, "component_heights") <- diagnostics(numeric())
    return(height)
  }

  densities <- lapply(terms[grid], `[[`, "density")
  grid_heights <- function(densities){
    vapply(densities, .prior_linear_density_grid_height, numeric(1), value = value)
  }
  mixture_height <- function(heights) fixed + sum(weights[grid] * heights)

  tolerance <- .prior_linear_density_refinement_tolerance()
  previous_heights <- grid_heights(densities)
  previous <- mixture_height(previous_heights)
  refined <- lapply(densities, .prior_linear_density_refinement)
  if(any(vapply(refined, is.null, logical(1)))){
    stop(
      "Adaptive prior-density evaluation did not converge within the documented ",
      "grid-refinement error criterion.",
      call. = FALSE
    )
  }
  for(i in seq_len(4L)){
    current_heights <- grid_heights(refined)
    current <- mixture_height(current_heights)
    change <- sum(weights[grid] * abs(current_heights - previous_heights))
    bound <- tolerance$absolute +
      tolerance$relative * max(abs(current), abs(previous))
    inside <- all(vapply(refined, function(density){
      !is.null(density$density) &&
        value >= min(density$density$x) && value <= max(density$density$x)
    }, logical(1)))
    if(isTRUE(inside) && is.finite(current) && change <= bound){
      attr(current, "adaptive_evaluation") <- list(
        components      = sum(grid),
        refinements     = i,
        absolute_change = change,
        error_bound     = bound,
        converged       = TRUE
      )
      attr(current, "component_heights") <- diagnostics(current_heights)
      return(current)
    }
    previous <- current
    previous_heights <- current_heights
    if(i < 4L){
      next_refined <- lapply(refined, .prior_linear_density_refinement)
      if(any(vapply(next_refined, is.null, logical(1)))){
        break
      }
      refined <- next_refined
    }
  }

  outside <- vapply(refined, function(density){
    final_range <- .prior_linear_density_range(density)
    value < final_range[1L] || value > final_range[2L]
  }, logical(1))
  if(any(outside)){
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

.prior_linear_density_height <- function(x, value){

  if(!inherits(x, "prior_linear_density")){
    stop("'x' must be a prior linear density object.", call. = FALSE)
  }

  ordinate <- NULL
  if(length(value) == 1L && is.finite(value)){
    ordinate <- .prior_density_ordinate_from_adaptive(
      attr(x, "adaptive_evaluation", exact = TRUE), value
    )
    exact <- .prior_linear_density_exact_height(ordinate)
    if(!is.null(exact)){
      return(exact)
    }
    components <- .prior_linear_density_component_height(
      attr(x, "adaptive_evaluation", exact = TRUE), value,
      grid_spacing = attr(x, "grid_resolution", exact = TRUE)[["spacing"]]
    )
    if(!is.null(components)){
      return(components)
    }
    .prior_linear_density_check_grid(.prior_density_route_from_adaptive(
      attr(x, "adaptive_evaluation", exact = TRUE)
    ))
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
