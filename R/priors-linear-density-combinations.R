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
  share <- .prior_ordered_linear_share(ordered_prior, weights, indices)
  if(identical(share$type,"point") && share$scale==0) return(c(0,0))
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
  share_range <- .prior_ordered_linear_share_range(
    share,
    tail_prob
  )
  products <- as.vector(outer(total_range, share_range, `*`))
  products <- products[is.finite(products)]
  if(length(products) == 0L){
    return(c(0, 0))
  }
  range(products)
}

# Numerical range of an allocation share (.prior_ordered_linear_share()): its
# scale for a fixed share, and the range of the scaled Beta term otherwise
# (from the prior alone, without a grid of the share).
.prior_ordered_linear_share_range <- function(share, tail_prob){

  if(identical(share$type, "point")){
    return(rep(share$scale, 2L))
  }
  .prior_linear_scalar_range(
    prior("beta", list(alpha = share$alpha[[1L]], beta = share$alpha[[2L]])),
    share$scale,
    tail_prob
  )
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

# Prior terms of an ordered level (or allocation subset): the ordered total
# with the share's scale as weight and, for a Beta(alpha_1, alpha_2) share,
# the share as its 'multiply_by' scale (named 'share_name', listed before the
# total 'total_name').
.prior_ordered_share_terms <- function(total, share, total_name, share_name){

  prior_list <- list()
  if(identical(share$type, "beta")){
    attr(total, "multiply_by") <- share_name
    prior_list[[share_name]] <- prior(
      "beta",
      list(alpha = share$alpha[[1L]], beta = share$alpha[[2L]])
    )
    attr(prior_list[[share_name]], "ordered_allocation") <- TRUE
  }
  prior_list[[total_name]] <- total
  list(prior_list = prior_list, total_name = total_name, weight = share$scale)
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

# Density of an ordered level (or allocation subset): the total times its
# allocation share, constructed from the level's structural route when it
# has one (.prior_linear_density_route_product()), and from the capped
# product grid of the total and the share only otherwise.
.prior_ordered_linear_distribution <- function(ordered_prior, weights, indices,
                                               dx = NA_real_, n_grid = NULL,
                                               tail_prob = .prior_linear_density_tail_prob()){

  if(is.null(n_grid)){
    n_grid <- .prior_linear_density_default_grid()
  }
  share <- .prior_ordered_linear_share(ordered_prior, weights, indices)
  if(identical(share$type,"point") && share$scale==0){
    out <- .prior_linear_density_point(0)
    attr(out,"ordered_measure") <- list(method="analytic_components",total=ordered_prior$total,
      allocation=.prior_ordered_metadata(.prior_ordered_default_bound(ordered_prior))$allocations[[1L]]$spec,
      coefficient_weights=stats::setNames(as.numeric(weights),names(weights)))
    return(out)
  }
  total <- .prior_ordered_total_linear_distribution(
    total = ordered_prior$total,
    dx = dx,
    n_grid = n_grid,
    tail_prob = tail_prob
  )

  # a Beta share: the product from its structural route, over the product of
  # the total's range and the share's range, with the total's zero atom
  out <- NULL
  if(identical(share$type, "beta")){
    terms <- .prior_ordered_share_terms(
      total      = ordered_prior$total,
      share      = share,
      total_name = ".ordered_total",
      share_name = ".ordered_share"
    )
    route <- .prior_density_route_linear(
      terms$prior_list, stats::setNames(terms$weight, terms$total_name), NULL, n_grid
    )
    out <- .prior_linear_density_route_product(
      route  = route,
      range  = range(outer(
        .prior_linear_density_range(total),
        .prior_ordered_linear_share_range(share, tail_prob)
      )),
      points = .prior_linear_density_product_atoms(
        total$points, .prior_linear_density_continuous_mass(total),
        .prior_linear_density_empty_points(), 1
      ),
      n_grid = n_grid
    )
  }
  if(is.null(out)){
    multiplier <- .prior_ordered_linear_multiplier(
      ordered_prior = ordered_prior,
      weights = weights,
      indices = indices,
      dx = dx,
      n_grid = n_grid,
      tail_prob = tail_prob
    )
    out <- .prior_linear_density_product(
      total,
      multiplier,
      n_grid = n_grid
    )
  }
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
  # the grid and its attributes (a deferred row-mixture grid is built)
  dist <- .prior_linear_density_materialize(dist)

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
  group_dists <- list()
  for(group in groups){
    group_dist <- .prior_linear_group_distribution(
      group             = group,
      dx                = dx,
      tail_prob         = tail_prob,
      source_transforms = source_transforms,
      n_grid            = n_grid
    )
    group_dists[[length(group_dists) + 1L]] <- group_dist
    dist <- .prior_linear_density_convolve(
      dist,
      .prior_linear_density_on_spacing(group_dist, dx, n_grid),
      dx
    )
  }

  attr(dist, "weights") <- weights
  attr(dist, "grid_resolution") <- c(spacing = dx, n_grid = n_grid)
  attr(dist, "product_grid_resolution") <- .prior_linear_density_merge_resolution(group_dists)
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

      # the product from its structural route; the capped product grid only
      # where the product has none
      product <- .prior_linear_density_route_product(
        route  = .prior_density_route_linear(
          prior_list, product_group$weights,
          source_transforms[names(product_group$weights)], n_grid
        ),
        range  = .prior_linear_density_product_range(linear_dist, multiplier_dist),
        points = .prior_linear_density_product_atoms(
          linear_dist$points, .prior_linear_density_continuous_mass(linear_dist),
          multiplier_dist$points, .prior_linear_density_continuous_mass(multiplier_dist)
        ),
        n_grid = n_grid
      )
      if(is.null(product)){
        product <- .prior_linear_density_product(
          linear_dist,
          multiplier_dist,
          n_grid = n_grid
        )
      }
      components[[length(components) + 1L]] <- product
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
  # whether the grids of product components resolve their densities
  # (.prior_linear_density_route_product()); plots of combinations without
  # a structural route omit a curve built from an unresolved one
  attr(dist, "product_grid_resolution") <- .prior_linear_density_merge_resolution(components)
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
  return(.prior_linear_density_normalize(dist, warn = TRUE))
}

# Acceptance of numerical prior heights and region probabilities is purely
# relative (at most 1e-4 of the value), so a small density (a far tail) is not
# accepted on an absolute floor; 'quadrature_floor' is only the absolute
# target of the first QUADPACK pass (see .prior_conditional_normal_quadrature).
.prior_linear_density_refinement_tolerance <- function(){

  list(relative = 1e-4, quadrature_floor = 1e-12)
}

# The one grid-refinement loop of prior-density heights and region
# probabilities. 'densities' are numerical grids with recorded provenance,
# 'evaluate' the functional of one grid (a height at a value or a region
# probability) and 'covers' whether a grid covers the evaluation. Each round
# refines every grid (halving its source spacing and omitting less tail
# probability) and applies the documented criterion to the combined value
# fixed + sum_k w_k g_k with the weighted absolute changes of the grids,
# sum_k w_k |g_k - g_k'| <= 1e-4 * max(|current|, |previous|) (purely
# relative), so changes cannot cancel between grids; at most four refinements
# follow the initial grids. Returns the converged values, their combination,
# the refined grids and the criterion; 'converged = FALSE' with the last
# grids; or NULL when a grid cannot be refined at all.
.prior_linear_density_refine_grids <- function(densities, weights, evaluate,
                                               covers = function(density) TRUE,
                                               fixed = 0){

  combine <- function(values) fixed + sum(weights * values)
  tolerance <- .prior_linear_density_refinement_tolerance()
  previous <- vapply(densities, evaluate, numeric(1))
  refined <- lapply(densities, .prior_linear_density_refinement)
  if(any(vapply(refined, is.null, logical(1)))){
    return(NULL)
  }
  for(i in seq_len(4L)){
    current <- vapply(refined, evaluate, numeric(1))
    total <- combine(current)
    change <- sum(weights * abs(current - previous))
    bound <- tolerance$relative * max(abs(total), abs(combine(previous)))
    if(all(vapply(refined, covers, logical(1))) && is.finite(total) &&
       change <= bound){
      return(list(
        converged       = TRUE,
        values          = current,
        total           = total,
        densities       = refined,
        refinements     = i,
        absolute_change = change,
        error_bound     = bound
      ))
    }
    previous <- current
    if(i < 4L){
      next_refined <- lapply(refined, .prior_linear_density_refinement)
      if(any(vapply(next_refined, is.null, logical(1)))){
        break
      }
      refined <- next_refined
    }
  }
  list(converged = FALSE, densities = refined)
}

# Whether a density grid covers 'value'.
.prior_linear_density_covers <- function(density, value){

  !is.null(density$density) &&
    value >= min(density$density$x) && value <= max(density$density$x)
}

# Stops after an unconverged refinement: an ordinate outside the final grids'
# range, or non-convergence of the 'quantity' ("density" or "probability").
.prior_linear_density_stop_refinement <- function(refinement, value = NULL,
                                                  quantity = "density"){

  if(!is.null(refinement) && !is.null(value)){
    outside <- vapply(refinement$densities, function(density){
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
  }
  stop(
    "Adaptive prior-", quantity, " evaluation did not converge within the ",
    "documented grid-refinement error criterion.",
    call. = FALSE
  )
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

# Integrand of the conditional-normal ordinate over the multiplier s at the
# values 'value' (vectorized over pairs): the conditional normal density of the
# value given s times the multiplier's density. Away from the offset, a pure
# scale mixture has no Gaussian mass at the value where the multiplier is zero
# (the integrand's limit).
.prior_conditional_normal_integrand <- function(spec){

  parts <- .prior_conditional_normal_integrand_parts(spec)
  function(multiplier, value){
    parts$value(parts$shared(multiplier), value)
  }
}

# The integrand in two parts: 'shared(multiplier)' holds the terms that do not
# depend on the value (the conditional moments and the multiplier's log
# density), and 'value(shared, value)' combines them with the value in the
# same arithmetic, so a batched quadrature evaluates the shared terms once per
# distinct node.
.prior_conditional_normal_integrand_parts <- function(spec){

  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  list(
    shared = function(multiplier){
      conditional <- .prior_conditional_normal_moments(spec, multiplier)
      list(mean = conditional$mean, sd = conditional$sd,
           log_density = multiplier_lpdf(multiplier))
    },
    value = function(shared, value){
      out <- exp(stats::dnorm(value, shared$mean, shared$sd, log = TRUE) +
                   shared$log_density)
      if(spec$additive_sd == 0){
        out[shared$sd == 0] <- 0
      }
      out
    }
  )
}

# Values of 'x' where the conditional-normal density is not the regular
# integral: the offset of a pure scale mixture.
.prior_conditional_normal_special <- function(spec, x){

  spec$additive_sd == 0 & x == spec$additive_mean
}

# Whether the conditional-normal quadrature at the values 'x' keeps full
# precision. A pure scale mixture a_m + b s (b ~ N(b_m, b_s)) is the scale
# product with a normal factor: its distance x - a_m from the offset and the
# distances standardized by b_s and b_m (the scale peak near
# s = |x - a_m| / b_s, the location peak near (x - a_m) / b_m) must be
# representable at full precision (.prior_density_affine_full_precision()),
# as for the scale-product leaf. With an additive normal term (a_s > 0) the
# value enters only through Gaussian kernels of x - a_m - b_m s with SD at
# least a_s, whose relative change under an absolute rounding of at most the
# smallest subnormal is negligible, so every value keeps full precision.
.prior_conditional_normal_full_precision <- function(spec, x){

  if(spec$additive_sd > 0){
    return(rep(TRUE, length(x)))
  }
  .prior_density_affine_full_precision(
    x, spec$additive_mean, c(spec$product_sd, spec$product_mean)
  )
}

# Batched plan of the conditional-normal density at the values 'x'
# (.prior_density_route_quadrature_density()): the values classified by the
# ordinate itself ('special'), and for the others the ordinate's integrand and
# breakpoints; 'batch' is FALSE when the integrand has an integrable
# singularity, i.e. the multiplier's density is infinite at a finite bound
# (except at zero for a pure scale mixture, where the Gaussian factor
# vanishes faster than any power). The integrand and breakpoints of such a
# plan are included only with 'singular' (see
# .prior_density_route_quadrature_density()), and never when the multiplier
# has a strong singularity at any bound (.prior_density_strong_singularity()).
.prior_conditional_normal_density_plan <- function(spec, x, singular = FALSE){

  # values without full precision take the ordinate, which has no value there
  special <- .prior_conditional_normal_special(spec, x) |
    !.prior_conditional_normal_full_precision(spec, x)
  setup <- .prior_conditional_normal_breakpoint_setup(spec$multiplier, spec$bounds)
  bounds <- c(setup$lower, setup$upper)
  batch <- !any(setup$singular & !(spec$additive_sd == 0 & bounds == 0))
  if(!batch && (!isTRUE(singular) ||
                .prior_density_strong_singularity(spec$multiplier, bounds))){
    return(list(batch = FALSE, special = special))
  }
  regular <- x[!special]
  parts <- .prior_conditional_normal_integrand_parts(spec)
  list(
    batch       = batch,
    special     = special,
    integrand   = list(
      shared = parts$shared,
      value  = function(shared, index) parts$value(shared, regular[index])
    ),
    breakpoints = .prior_conditional_normal_breakpoints_values(spec, regular, setup)
  )
}

.prior_conditional_normal_ordinate <- function(spec, value, n_grid){

  # A pure scale mixture (a_s = 0) at its offset a_m is classified from the
  # multiplier's declared behavior at zero.
  if(.prior_conditional_normal_special(spec, value)){
    return(.prior_conditional_normal_offset_ordinate(spec, value, n_grid))
  }
  if(!.prior_conditional_normal_full_precision(spec, value)){
    return(.prior_density_ordinate_imprecise(
      value, "The distance of the requested value from the scale mixture's offset",
      "conditional_normal_mixture",
      list(kind = "conditional_normal_mixture",
           additive = c(mean = spec$additive_mean, sd = spec$additive_sd),
           multiplied = c(mean = spec$product_mean, sd = spec$product_sd),
           multiplier = .prior_density_ordinate_prior_provenance(spec$multiplier),
           independent_sources = spec$sources,
           structural_regularity = "scale_mixture_away_from_offset")
    ))
  }
  # The integral runs over the other term's (the multiplier's) declared
  # support, split at breakpoints so that no piece is dominated by mass that
  # its initial quadrature rule cannot see (a narrow Gaussian peak, or a
  # concentrated other term far from zero). Each piece is an independent
  # integral with the full evaluation budget and its own diagnostics; the
  # ordinate is their sum, and the acceptance criterion applies to the total.
  integrand <- .prior_conditional_normal_integrand(spec)
  integral <- .prior_conditional_normal_quadrature(
    function(multiplier) integrand(multiplier, value),
    .prior_conditional_normal_breakpoints(spec, value), n_grid,
    zero_message = "zero ordinate for a structurally positive density"
  )
  .prior_density_quadrature_ordinate(
    value, integral, "conditional_normal_mixture",
    list(
      kind = "conditional_normal_mixture",
      additive = c(mean = spec$additive_mean, sd = spec$additive_sd),
      multiplied = c(mean = spec$product_mean, sd = spec$product_sd),
      multiplier = .prior_density_ordinate_prior_provenance(spec$multiplier),
      independent_sources = spec$sources,
      structural_regularity = if(spec$additive_sd > 0){
        "positive_variance_gaussian_convolution"
      }else{
        "scale_mixture_away_from_offset"
      }
    )
  )
}

# The regular ordinate of a quadrature leaf (conditional-normal mixture, scale
# product, two-term convolution) from its integral, with the integration
# record in its provenance. An accepted integral value below
# .Machine$double.xmin is subnormal and rounded to a multiple of the smallest
# subnormal, so its log has no value there (the full-precision rule,
# .prior_density_full_precision(): the log image of a half-Cauchy product was
# off by 1.2e-2 at z = 370, where the product's density is 4e-322); plotted
# densities keep it as a display estimate.
.prior_density_quadrature_ordinate <- function(value, integral, method, provenance){

  provenance$integration <- integral$integration
  if(isTRUE(integral$value > 0) && !.prior_density_full_precision(integral$value)){
    provenance$integration$estimate <- integral$value
    return(.prior_density_ordinate_imprecise(
      value, "The density at the requested value", method, provenance
    ))
  }
  .prior_density_ordinate_result(
    value = value, behavior = "regular",
    log_density = log(integral$value), exact = TRUE,
    method = method, provenance = provenance
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

# Scale products X = c + w * L * m(s) of a simple continuous term L (not a
# full-support normal, which is the conditional-normal route) and a simple
# continuous multiplier s mapped by m: the ordered-level route with a
# non-normal total (L = total, s = its Beta allocation share, m the identity),
# 'multiply_by' products of a non-normal coefficient prior (m the identity),
# and allocation-derived random-effect SDs (L = the scale prior, s = a Beta
# allocation share, m(s) = sqrt(k s), the 'sqrt' map with scale k). Away from
# the offset c the density is the 1-D integral
# f(x) = int f_s(s) f_L((x - c) / (w m(s))) / |w m(s)| ds over the
# multiplier's support, evaluated by the conditional-normal quadrature (split
# at the multiplier's bounds and quantiles, at its zero, and at the images
# s = m^-1((x - c) / (w q)) of L's quantiles and finite bounds q); the offset
# is classified from the declared behaviors of L and m(s) at zero.
.prior_scale_product_spec <- function(offset, scale, factor, multiplier, sources,
                                      map = NULL){

  if(!is.null(map)){
    beta_share <- identical(multiplier$distribution, "beta") &&
      isTRUE(multiplier$truncation$lower == 0) &&
      isTRUE(multiplier$truncation$upper == 1)
    if(!identical(map$type, "sqrt") || !is.numeric(map$scale) ||
       length(map$scale) != 1L || !is.finite(map$scale) || map$scale <= 0 ||
       !beta_share){
      stop("Scale-product multiplier maps are square roots of scaled Beta shares.",
           call. = FALSE)
    }
  }
  list(
    offset     = offset,
    scale      = scale,
    factor     = factor,
    multiplier = multiplier,
    bounds     = unlist(multiplier$truncation[c("lower", "upper")], use.names = FALSE),
    sources    = sources,
    map        = map
  )
}

# The mapped multiplier m(s) of a scale product (the identity without a map)
# and its inverse on the mapped values (NA where no share maps to them).
.prior_scale_product_map <- function(spec, s){

  if(is.null(spec$map)){
    return(s)
  }
  sqrt(spec$map$scale * s)
}

.prior_scale_product_map_inverse <- function(spec, v){

  if(is.null(spec$map)){
    return(v)
  }
  out <- v^2 / spec$map$scale
  out[!is.finite(v) | v < 0] <- NA_real_
  out
}

# Quadrature settings of a scale product: the multiplier's breakpoints, no
# Gaussian peak window. 'setup' holds the value-independent parts.
.prior_scale_product_breakpoint_setup <- function(spec){

  factor_bounds <- unlist(spec$factor$truncation[c("lower", "upper")], use.names = FALSE)
  quantiles <- tryCatch(
    suppressWarnings(as.numeric(quant(
      spec$factor, c(1e-6, 1e-3, .02, .25, .5, .75, .98, 1 - 1e-3, 1 - 1e-6)
    ))),
    error = function(e) numeric()
  )
  targets <- c(quantiles, factor_bounds)
  list(
    targets    = targets[is.finite(targets) & targets != 0],
    spec       = list(additive_mean = 0, additive_sd = 1, product_mean = 0, product_sd = 0,
                      multiplier = spec$multiplier, bounds = spec$bounds),
    multiplier = .prior_conditional_normal_breakpoint_setup(spec$multiplier, spec$bounds)
  )
}

.prior_scale_product_breakpoints <- function(spec, distances,
                                             setup = .prior_scale_product_breakpoint_setup(spec)){

  images <- .prior_scale_product_map_inverse(
    spec, as.vector(outer(distances, setup$targets, `/`))
  )
  .prior_conditional_normal_breakpoints(
    setup$spec,
    value = 0,
    extra = c(0, images[is.finite(images)]),
    setup = setup$multiplier
  )
}

# .prior_scale_product_breakpoints() of each of the single 'distances',
# computed together (.prior_conditional_normal_breakpoints_values()).
.prior_scale_product_breakpoints_values <- function(spec, distances, setup){

  images <- .prior_scale_product_map_inverse(
    spec, as.vector(outer(distances, setup$targets, `/`))
  )
  images <- matrix(images, nrow = length(distances))
  images[!is.finite(images)] <- NA_real_
  .prior_conditional_normal_breakpoints_values(
    setup$spec,
    values = rep(0, length(distances)),
    setup  = setup$multiplier,
    extra  = cbind(0, images)
  )
}

# Integrand of the scale-product ordinate over the multiplier s at the
# standardized distances (x - c) / w (vectorized over pairs).
.prior_scale_product_integrand <- function(spec){

  parts <- .prior_scale_product_integrand_parts(spec)
  function(multiplier, distance){
    parts$value(parts$shared(multiplier), distance)
  }
}

# The integrand in a part shared by all distances (the mapped multiplier, its
# log density and the log Jacobian) and a part combining it with the distance
# (as .prior_conditional_normal_integrand_parts()).
.prior_scale_product_integrand_parts <- function(spec){

  factor_lpdf <- .prior_simple_lpdf_evaluator(spec$factor)
  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  list(
    shared = function(multiplier){
      mapped <- .prior_scale_product_map(spec, multiplier)
      list(mapped = mapped, log_density = multiplier_lpdf(multiplier),
           log_jacobian = log(abs(spec$scale * mapped)))
    },
    value = function(shared, distance){
      out <- exp(factor_lpdf(distance / shared$mapped) +
                   shared$log_density - shared$log_jacobian)
      out[shared$mapped == 0] <- 0
      out
    }
  )
}

# Batched plan of the scale-product density at the values 'x' (as
# .prior_conditional_normal_density_plan()): values outside the support hull
# have a zero density ('zero'), and its bounds and the offset are classified
# by the ordinate, as are values whose distance from the offset is not
# representable at full precision (the ordinate has no value there). The
# integrand is singular ('batch' FALSE; its integrand and
# breakpoints only with 'singular', and never when the multiplier or the
# factor has a strong singularity at any bound,
# .prior_density_strong_singularity()) where the multiplier's density is
# infinite at a nonzero finite bound, and at the image of a nonzero finite
# bound where the factor's density is infinite.
.prior_scale_product_density_plan <- function(spec, x, singular = FALSE){

  hull <- .prior_scale_product_hull(spec)
  zero <- x < hull[1L] | x > hull[2L]
  special <- !zero & (x == hull[1L] | x == hull[2L] | x == spec$offset |
                        !.prior_density_affine_full_precision(x, spec$offset, spec$scale))
  setup <- .prior_scale_product_breakpoint_setup(spec)
  factor_bounds <- unlist(spec$factor$truncation[c("lower", "upper")], use.names = FALSE)
  # with a square-root map, a share density that is infinite at zero is not
  # cancelled by the factor's tail (a heavy-tailed factor leaves s^(a - 1/2))
  batch <- !any(setup$multiplier$singular & (spec$bounds != 0 | !is.null(spec$map))) &&
    !any(.prior_density_singular_bounds(spec$factor) & factor_bounds != 0)
  if(!batch && (!isTRUE(singular) ||
                .prior_density_strong_singularity(spec$multiplier, spec$bounds) ||
                .prior_density_strong_singularity(spec$factor, factor_bounds))){
    return(list(batch = FALSE, zero = zero, special = special))
  }
  distance <- (x[!zero & !special] - spec$offset) / spec$scale
  parts <- .prior_scale_product_integrand_parts(spec)
  list(
    batch       = batch,
    zero        = zero,
    special     = special,
    integrand   = list(
      shared = parts$shared,
      value  = function(shared, index) parts$value(shared, distance[index])
    ),
    breakpoints = .prior_scale_product_breakpoints_values(spec, distance, setup)
  )
}

# Closed interval containing the support of c + w * L * m(s).
.prior_scale_product_hull <- function(spec){

  factor_bounds <- unlist(spec$factor$truncation[c("lower", "upper")], use.names = FALSE)
  products <- as.vector(outer(factor_bounds, .prior_scale_product_map(spec, spec$bounds), function(a, b){
    ifelse(a == 0 | b == 0, 0, a * b)
  }))
  spec$offset + sort(spec$scale * range(products))
}

.prior_scale_product_provenance <- function(spec){

  out <- list(
    kind                = "scale_mixture",
    offset              = spec$offset,
    scale               = spec$scale,
    factor              = .prior_density_ordinate_prior_provenance(spec$factor),
    multiplier          = .prior_density_ordinate_prior_provenance(spec$multiplier),
    independent_sources = spec$sources
  )
  if(!is.null(spec$map)){
    out$multiplier_map <- spec$map
  }
  out
}

# Structural provenance of a scale-product route (named transformations of
# the route read its support and its behavior at the offset).
.prior_scale_product_route_provenance <- function(spec){

  provenance <- .prior_scale_product_provenance(spec)
  provenance$support <- .prior_scale_product_hull(spec)
  provenance$offset_behavior <- .prior_scale_product_offset_behavior(spec)$behavior
  provenance
}

# Behavior of the mapped multiplier m(s) at zero (the one-sided limit at a
# bound): the share's own classification for the identity map, and for
# m(s) = sqrt(k s) of a Beta(a, b) share, whose density near zero is
# 2 v^(2a - 1) / (k^a B(a, b)), zero for a > 1/2, 2 / (sqrt(k) B(1/2, b)) for
# a = 1/2 and infinite for a < 1/2.
.prior_scale_product_multiplier_zero <- function(spec){

  if(is.null(spec$map)){
    return(.prior_density_ordinate_primitive(spec$multiplier, 0))
  }
  alpha <- spec$multiplier$parameters$alpha
  beta <- spec$multiplier$parameters$beta
  behavior <- if(alpha > 1 / 2) "zero" else if(alpha < 1 / 2) "infinite" else "regular"
  list(
    behavior    = behavior,
    log_density = switch(
      behavior,
      "zero"     = -Inf,
      "infinite" = Inf,
      "regular"  = log(2) - log(spec$map$scale) / 2 - lbeta(1 / 2, beta)
    )
  )
}

# E[1 / |m(s)|] of the mapped multiplier when its density vanishes at zero:
# E[(k s)^(-1/2)] = B(a - 1/2, b) / (sqrt(k) B(a, b)) for a Beta(a, b) share
# (a > 1/2) and the multiplier's own inverse moment for the identity map.
.prior_scale_product_inverse_moment <- function(spec, n_grid){

  if(is.null(spec$map)){
    return(.prior_density_inverse_moment(spec$multiplier, n_grid))
  }
  alpha <- spec$multiplier$parameters$alpha
  beta <- spec$multiplier$parameters$beta
  list(
    value  = exp(lbeta(alpha - 1 / 2, beta) - lbeta(alpha, beta) -
                   log(spec$map$scale) / 2),
    method = "closed_form",
    integration = NULL
  )
}

# The scale-product quadrature at a value keeps full double precision only when
# the distance x - c from the offset and the standardized distance (x - c) / w
# are representable at full precision (.prior_density_affine_full_precision()):
# at a subnormal distance the factor's argument (x - c) / (w m(s)) of the
# integrand is rounded to a multiple of the smallest subnormal, which moves the
# factor's density where it varies near zero (for a gamma(2, 4) factor and a
# lognormal multiplier, by 1.4e-4 relative at 1e-320 and 7e-2 at 4.9e-324).
# Such a value has no ordinate value, whichever route reaches the leaf
# (products, ordered levels, allocations, transformations, log images).
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

  provenance$structural_regularity <- "scale_mixture_inside_support"
  if(!.prior_density_affine_full_precision(value, spec$offset, spec$scale)){
    return(.prior_density_ordinate_imprecise(
      value, "The distance of the requested value from the product's offset",
      "scale_mixture", provenance
    ))
  }
  distance <- (value - spec$offset) / spec$scale
  integrand <- .prior_scale_product_integrand(spec)
  integral <- .prior_conditional_normal_quadrature(
    function(multiplier) integrand(multiplier, distance),
    .prior_scale_product_breakpoints(spec, distance), n_grid,
    zero_message = "zero ordinate for a structurally positive density",
    kind = "scale_mixture"
  )
  .prior_density_quadrature_ordinate(value, integral, "scale_mixture", provenance)
}

# The offset c of c + w * L * s: with f_L and f_s the declared densities at
# zero (one-sided limits at a support bound), the density at c is infinite
# when either is infinite or both are positive, f_L(0) E[1 / |s|] / |w| when
# only f_s vanishes there, f_s(0) E[1 / |L|] / |w| when only f_L vanishes
# there, and zero when both vanish. A term bounded at zero with a positive
# finite limit there, combined with a two-sided other term whose density
# vanishes at zero, makes the offset a density jump (the one-sided limits
# weight E[1 / |s|] over the two signs of the other term separately), which
# is not classified. With a two-sided other term whose density at zero is
# positive or infinite both one-sided limits are infinite (e.g. an ordered
# level of a t total with a Beta(1, b) share), so no jump occurs.
# Classification of the offset c of a scale product from the declared
# behaviors of L and m(s) at zero (see .prior_scale_product_offset_ordinate()):
# 'behavior' is "infinite", "zero", "regular", or "unknown" (with 'reason').
.prior_scale_product_offset_behavior <- function(spec){

  factor_zero <- .prior_density_ordinate_primitive(spec$factor, 0)
  multiplier_zero <- .prior_scale_product_multiplier_zero(spec)
  behaviors <- c(
    factor     = .prior_density_ordinate_continuous_behavior(factor_zero),
    multiplier = .prior_density_ordinate_continuous_behavior(multiplier_zero)
  )
  out <- list(behaviors = behaviors, factor_zero = factor_zero,
              multiplier_zero = multiplier_zero, reason = NULL)
  if(any(!behaviors %in% c("regular", "zero", "infinite"))){
    out$behavior <- "unknown"
    return(out)
  }
  bounds_of <- function(prior, mapped){
    bounds <- unlist(prior$truncation[c("lower", "upper")], use.names = FALSE)
    if(mapped) .prior_scale_product_map(spec, bounds) else bounds
  }
  one_sided <- function(bounds) any(bounds == 0)
  two_sided <- function(bounds) bounds[1L] < 0 && bounds[2L] > 0
  factor_bounds <- bounds_of(spec$factor, FALSE)
  multiplier_bounds <- bounds_of(spec$multiplier, TRUE)
  if((one_sided(factor_bounds) && two_sided(multiplier_bounds) &&
      behaviors[["factor"]] == "regular" && behaviors[["multiplier"]] == "zero") ||
     (one_sided(multiplier_bounds) && two_sided(factor_bounds) &&
      behaviors[["multiplier"]] == "regular" && behaviors[["factor"]] == "zero")){
    out$behavior <- "unknown"
    out$reason <- "The density of the product jumps at the requested value."
    return(out)
  }
  out$behavior <- if(any(behaviors == "infinite") || all(behaviors == "regular")){
    "infinite"
  }else if(all(behaviors == "zero")){
    "zero"
  }else{
    "regular"
  }
  out
}

.prior_scale_product_offset_ordinate <- function(spec, value, n_grid, provenance){

  offset <- .prior_scale_product_offset_behavior(spec)
  behaviors <- offset$behaviors
  provenance$structural_regularity <- "scale_mixture_offset"
  provenance$behaviors_at_zero <- behaviors
  unknown <- function(reason){
    .prior_density_ordinate_result(
      value = value, behavior = "unknown", log_density = NA_real_,
      exact = FALSE, method = "unsupported_provenance", reason = reason,
      provenance = provenance
    )
  }
  if(identical(offset$behavior, "unknown")){
    return(unknown(offset$reason))
  }
  if(identical(offset$behavior, "infinite")){
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
  if(identical(offset$behavior, "zero")){
    return(.prior_density_ordinate_result(
      value = value, behavior = "zero", log_density = -Inf, exact = TRUE,
      method = "scale_mixture", provenance = provenance
    ))
  }
  if(behaviors[["factor"]] == "regular"){
    density_zero <- offset$factor_zero$log_density
    moment <- .prior_scale_product_inverse_moment(spec, n_grid)
  }else{
    density_zero <- offset$multiplier_zero$log_density
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

# Region probability of c + w * L * m(s): the 1-D integral over s of its
# density times P(c + w L m(s) in region), with the ordinate's breakpoints at
# the images of every finite region endpoint.
.prior_scale_product_region <- function(spec, intervals, n_grid){

  lower <- intervals[, 1L]
  upper <- intervals[, 2L]
  offset_inside <- as.numeric(any(spec$offset > lower & spec$offset < upper))
  multiplier_lpdf <- .prior_simple_lpdf_evaluator(spec$multiplier)
  integrand <- function(multiplier){
    mapped <- .prior_scale_product_map(spec, multiplier)
    probability <- numeric(length(multiplier))
    zero <- mapped == 0
    s <- mapped[!zero]
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

# Two-term convolutions X = c + w_A A + w_B B of simple continuous scalar
# terms without a Gaussian part (a Gaussian term plus one other term is the
# conditional-normal Gaussian convolution): the density is the 1-D integral
# f(x) = int f_A(a) f_B((x - c - w_A a) / w_B) / |w_B| da over A's support,
# evaluated by the conditional-normal quadrature split at A's bounds and
# quantiles and at the images a = (x - c - w_B q) / w_A of B's quantiles and
# finite bounds q. Values where a finite bound of A meets a finite bound of B
# are classified from the densities' exponents at those bounds (1 for a
# positive finite density, > 1 for a vanishing one, and the shape parameter
# of a gamma or beta density that is infinite there): with e = p_A + p_B - 1,
# an end of the support has a zero density for e > 0 and an infinite one for
# e < 0, and an inner meeting point has an infinite density for e <= 0.
.prior_convolution_spec <- function(offset, terms, sources){

  list(
    offset  = offset,
    first   = terms[[1L]]$prior,
    weight  = terms[[1L]]$weight,
    second  = terms[[2L]]$prior,
    other   = terms[[2L]]$weight,
    bounds  = unlist(terms[[1L]]$prior$truncation[c("lower", "upper")], use.names = FALSE),
    sources = sources
  )
}

.prior_convolution_provenance <- function(spec){

  list(
    kind                = "convolution",
    offset              = spec$offset,
    weights             = c(spec$weight, spec$other),
    terms               = list(
      .prior_density_ordinate_prior_provenance(spec$first),
      .prior_density_ordinate_prior_provenance(spec$second)
    ),
    independent_sources = spec$sources
  )
}

# Closed interval containing the support of c + w_A A + w_B B.
.prior_convolution_hull <- function(spec){

  second_bounds <- unlist(spec$second$truncation[c("lower", "upper")], use.names = FALSE)
  spec$offset + range(spec$weight * spec$bounds) + range(spec$other * second_bounds)
}

# Exponent of a simple continuous prior's density at a finite support bound:
# 1 when positive and finite, 2 (any value above 1 classifies alike) when it
# vanishes, the shape of an infinite gamma or beta density, NA otherwise.
.prior_density_bound_exponent <- function(prior, bound){

  behavior <- .prior_density_ordinate_continuous_behavior(
    .prior_density_ordinate_primitive(prior, bound)
  )
  if(identical(behavior, "regular")){
    return(1)
  }
  if(identical(behavior, "zero")){
    return(2)
  }
  if(!identical(behavior, "infinite")){
    return(NA_real_)
  }
  parameters <- prior$parameters
  if(identical(prior$distribution, "gamma") && bound == 0){
    return(parameters$shape)
  }
  if(identical(prior$distribution, "beta") && bound == 0){
    return(parameters$alpha)
  }
  if(identical(prior$distribution, "beta") && bound == 1){
    return(parameters$beta)
  }
  NA_real_
}

# Quadrature settings of a two-term convolution over the first term: its
# breakpoints and the images of the second term's quantiles and bounds.
# 'setup' holds the value-independent parts.
.prior_convolution_breakpoint_setup <- function(spec){

  second_bounds <- unlist(spec$second$truncation[c("lower", "upper")], use.names = FALSE)
  quantiles <- tryCatch(
    suppressWarnings(as.numeric(quant(
      spec$second, c(1e-6, 1e-3, .02, .25, .5, .75, .98, 1 - 1e-3, 1 - 1e-6)
    ))),
    error = function(e) numeric()
  )
  targets <- c(quantiles, second_bounds)
  list(
    targets = targets[is.finite(targets)],
    spec    = list(additive_mean = 0, additive_sd = 1, product_mean = 0, product_sd = 0,
                   multiplier = spec$first, bounds = spec$bounds),
    first   = .prior_conditional_normal_breakpoint_setup(spec$first, spec$bounds)
  )
}

.prior_convolution_breakpoints <- function(spec, distances,
                                           setup = .prior_convolution_breakpoint_setup(spec)){

  images <- as.vector(outer(distances, setup$targets, function(d, q) (d - spec$other * q) / spec$weight))
  .prior_conditional_normal_breakpoints(
    setup$spec,
    value = 0,
    extra = images[is.finite(images)],
    setup = setup$first
  )
}

# Integrand of the convolution ordinate over the first term at the distances
# x - c (vectorized over pairs).
.prior_convolution_integrand <- function(spec){

  first_lpdf <- .prior_simple_lpdf_evaluator(spec$first)
  second_lpdf <- .prior_simple_lpdf_evaluator(spec$second)
  function(first, distance){
    exp(first_lpdf(first) +
          second_lpdf((distance - spec$weight * first) / spec$other) -
          log(abs(spec$other)))
  }
}

# Whether the convolution quadrature at the values 'x' keeps full precision:
# the distance d = x - c from the offset enters the other term's argument
# (d - w_A a) / w_B and the breakpoints (d - w_B q) / w_A, so d, d / w_A and
# d / w_B must be representable at full precision
# (.prior_density_affine_full_precision()). Next to a meeting point of bounds
# at the offset, a subnormal d leaves the integrand a subnormal support.
.prior_convolution_full_precision <- function(spec, x){

  .prior_density_affine_full_precision(x, spec$offset, c(spec$weight, spec$other))
}

# Meeting points c + w_A a + w_B b of finite support bounds of both terms.
.prior_convolution_meeting_points <- function(spec){

  second_bounds <- unlist(spec$second$truncation[c("lower", "upper")], use.names = FALSE)
  first_bounds <- spec$bounds[is.finite(spec$bounds)]
  second_bounds <- second_bounds[is.finite(second_bounds)]
  as.vector(outer(first_bounds, second_bounds, function(a, b){
    spec$offset + spec$weight * a + spec$other * b
  }))
}

# Batched plan of the convolution density at the values 'x' (as
# .prior_conditional_normal_density_plan()): values outside the support hull
# have a zero density ('zero'), and meeting points of support bounds are
# classified by the ordinate. The integrand is singular where either term's
# density is infinite at a finite bound.
.prior_convolution_density_plan <- function(spec, x){

  hull <- .prior_convolution_hull(spec)
  zero <- x < hull[1L] | x > hull[2L]
  special <- rep(FALSE, length(x))
  for(meeting in .prior_convolution_meeting_points(spec)){
    special <- special | vapply(x, function(value){
      isTRUE(.prior_density_ordinate_endpoint_matches(meeting, value))
    }, logical(1))
  }
  # values without full precision take the ordinate, which has no value there
  special <- (special | !.prior_convolution_full_precision(spec, x)) & !zero
  setup <- .prior_convolution_breakpoint_setup(spec)
  batch <- !any(setup$first$singular) && !any(.prior_density_singular_bounds(spec$second))
  if(!batch){
    return(list(batch = FALSE, zero = zero, special = special))
  }
  distance <- x[!zero & !special] - spec$offset
  integrand <- .prior_convolution_integrand(spec)
  list(
    batch       = TRUE,
    zero        = zero,
    special     = special,
    integrand   = function(first, index) integrand(first, distance[index]),
    breakpoints = lapply(distance, .prior_convolution_breakpoints, spec = spec, setup = setup)
  )
}

.prior_convolution_ordinate <- function(spec, value, n_grid){

  provenance <- .prior_convolution_provenance(spec)
  hull <- .prior_convolution_hull(spec)
  result <- function(behavior, log_density, method = "convolution", reason = NULL,
                     exact = TRUE){
    .prior_density_ordinate_result(
      value = value, behavior = behavior, log_density = log_density,
      exact = exact, method = method, reason = reason, provenance = provenance
    )
  }
  if(value < hull[1L] || value > hull[2L]){
    return(result("zero", -Inf, reason = "The requested value is outside the prior support."))
  }

  # meeting points of finite bounds of both terms
  second_bounds <- unlist(spec$second$truncation[c("lower", "upper")], use.names = FALSE)
  first_bounds <- spec$bounds[is.finite(spec$bounds)]
  second_bounds <- second_bounds[is.finite(second_bounds)]
  exponent <- Inf
  for(a in first_bounds){
    for(b in second_bounds){
      meeting <- spec$offset + spec$weight * a + spec$other * b
      if(isTRUE(.prior_density_ordinate_endpoint_matches(meeting, value))){
        exponent <- min(exponent, .prior_density_bound_exponent(spec$first, a) +
                          .prior_density_bound_exponent(spec$second, b) - 1, na.rm = FALSE)
      }
    }
  }
  if(is.na(exponent)){
    return(result("unknown", NA_real_, method = "unsupported_provenance", exact = FALSE,
                  reason = "The density where two support bounds meet is not structurally classified."))
  }
  if(is.finite(exponent)){
    at_end <- value == hull[1L] || value == hull[2L]
    if(at_end && exponent > 0){
      return(result("zero", -Inf))
    }
    if(exponent < 0 || (!at_end && exponent == 0)){
      return(result("infinite", Inf, reason = paste0(
        "The densities of the two terms are singular where their support ",
        "bounds meet, which makes the density of the sum infinite at the ",
        "requested value."
      )))
    }
    if(at_end){
      return(result("unknown", NA_real_, method = "unsupported_provenance", exact = FALSE,
                    reason = "The positive limit at an end of the support is not structurally classified."))
    }
  }

  if(!.prior_convolution_full_precision(spec, value)){
    return(.prior_density_ordinate_imprecise(
      value, "The distance of the requested value from the convolution's offset",
      "convolution", provenance
    ))
  }
  distance <- value - spec$offset
  integrand <- .prior_convolution_integrand(spec)
  integral <- .prior_conditional_normal_quadrature(
    function(first) integrand(first, distance),
    .prior_convolution_breakpoints(spec, distance), n_grid,
    zero_message = "zero ordinate for a structurally positive density",
    kind = "convolution"
  )
  .prior_density_quadrature_ordinate(value, integral, "convolution", provenance)
}

# Region probability of c + w_A A + w_B B: the 1-D integral over A of its
# density times P(c + w_A a + w_B B in region), with the ordinate's
# breakpoints at the images of every finite region endpoint.
.prior_convolution_region <- function(spec, intervals, n_grid){

  lower <- intervals[, 1L]
  upper <- intervals[, 2L]
  first_lpdf <- .prior_simple_lpdf_evaluator(spec$first)
  integrand <- function(first){
    probability <- numeric(length(first))
    for(i in seq_along(lower)){
      a <- (lower[i] - spec$offset - spec$weight * first) / spec$other
      b <- (upper[i] - spec$offset - spec$weight * first) / spec$other
      if(spec$other < 0){
        probability <- probability + .prior_scalar_interval_probability(spec$second, b, a)
      }else{
        probability <- probability + .prior_scalar_interval_probability(spec$second, a, b)
      }
    }
    out <- numeric(length(first))
    positive <- probability > 0
    out[positive] <- exp(log(probability[positive]) + first_lpdf(first[positive]))
    out
  }
  endpoints <- c(lower, upper)
  endpoints <- unique(endpoints[is.finite(endpoints)])
  .prior_conditional_normal_quadrature(
    integrand,
    .prior_convolution_breakpoints(spec, endpoints - spec$offset),
    n_grid,
    zero_message = "zero probability for a structurally positive region",
    kind = "convolution"
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
# criterion (all pieces converged, value positive, and the purely relative
# abs. error <= 1e-4 * value) applies to the total. The pieces first run with
# QUADPACK's relative target 1e-4 and an absolute floor of 1e-12 in total; a
# total that misses the relative criterion (a small density, e.g. a far tail,
# where the floor stopped QUADPACK early) is refined once against its own
# value, each piece with its full budget again. A rejected total has value NA
# (a total that only misses the relative criterion keeps a display
# 'estimate'). An optional 'exact_piece(lower, upper)' returns a piece
# evaluated without quadrature (value, abs.error, message, evaluations) or
# NULL for a QUADPACK piece; in the refinement an exact piece evaluated as 0
# whose error bound is too large for the total is integrated instead.
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
      relative = tolerance$relative, absolute = tolerance$quadrature_floor / n_pieces
    )
  })
  total <- function(pieces){
    piece_messages <- vapply(pieces, `[[`, character(1), "message")
    integral <- list(
      value     = sum(vapply(pieces, `[[`, numeric(1), "value")),
      abs.error = sum(vapply(pieces, `[[`, numeric(1), "abs.error")),
      message   = if(all(piece_messages == "OK")){
        "OK"
      }else{
        paste(unique(piece_messages[piece_messages != "OK"]), collapse = "; ")
      }
    )
    if(identical(integral$message, "OK") && isTRUE(integral$value == 0)){
      integral$message <- zero_message
    }
    integral
  }
  # acceptance is purely relative: the reported error of the total at most
  # 1e-4 of its value (the integrands are nonnegative, so the pieces add up)
  accept <- function(integral){
    identical(integral$message, "OK") && is.finite(integral$value) &&
      integral$value > 0 && is.finite(integral$abs.error) &&
      integral$abs.error <= tolerance$relative * integral$value
  }
  integral <- total(pieces)
  evaluations <- sum(vapply(pieces, `[[`, integer(1), "evaluations"))
  refined <- FALSE
  if(!accept(integral) && identical(integral$message, "OK") &&
     is.finite(integral$value) && integral$value > 0 &&
     is.finite(integral$abs.error)){
    # the absolute floor of the first pass stopped QUADPACK before the
    # relative criterion (a small total): every quadrature piece is refined
    # against the total (per piece max(rel / 2 * |piece|, rel * total /
    # (2 n)) sums to at most rel * total), with its full budget again. An
    # exact piece evaluated as 0 whose error bound (a fixed fraction of its
    # mass) exceeds that share is integrated as well; exact pieces evaluated
    # as their mass have a bound relative to their own value and are kept.
    refined <- TRUE
    target <- tolerance$relative * integral$value / (2 * n_pieces)
    quadrature <- !exact | vapply(pieces, function(piece){
      piece$value == 0 && piece$abs.error > target
    }, logical(1))
    pieces[quadrature] <- lapply(which(quadrature), function(i){
      .prior_conditional_normal_piece(
        integrand, points[i], points[i + 1L], n_grid,
        relative = tolerance$relative / 2, absolute = target
      )
    })
    exact <- exact & !quadrature
    integral <- total(pieces)
    evaluations <- evaluations +
      sum(vapply(pieces[quadrature], `[[`, integer(1), "evaluations"))
  }
  bound <- tolerance$relative * abs(integral$value)
  relative_error <- integral$abs.error / abs(integral$value)
  accepted <- accept(integral)
  # a total that only misses the relative criterion is kept as a display
  # estimate (plotted curves); it is never an ordinate or a probability
  estimate <- if(!accepted && identical(integral$message, "OK") &&
                 is.finite(integral$value) && integral$value > 0){
    integral$value
  }
  if(!isTRUE(accepted)){
    integral$value <- NA_real_
  }
  integration <- list(
    kind = kind, exact = FALSE,
    absolute_error = integral$abs.error, error_bound = bound,
    evaluations = evaluations,
    budget = n_grid, converged = isTRUE(accepted), message = integral$message,
    refined = refined,
    breakpoints = points,
    piece_evaluations = vapply(pieces, `[[`, integer(1), "evaluations"),
    piece_absolute_errors = vapply(pieces, `[[`, numeric(1), "abs.error")
  )
  if(!accepted && identical(integral$message, "OK")){
    # the observed metric, not the criterion (public-API message rules); the
    # criterion stays in 'error_bound'
    integration$message <- paste0(
      "relative error estimate ", format(signif(relative_error, 3))
    )
    integration$estimate <- estimate
  }
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
# images of another term's quantiles) are added like quantiles. 'setup' holds
# the value-independent parts (support bounds, singular bounds and quantiles of
# the multiplier), so that batched densities compute them once.
.prior_conditional_normal_breakpoint_setup <- function(multiplier, bounds){

  lower <- bounds[1L]
  upper <- bounds[2L]
  singular <- .prior_density_singular_bounds(multiplier, c(lower, upper))
  probabilities <- c(if(!singular[1L]) 1e-6, 1e-3, .02, .25, .5, .75, .98, 1 - 1e-3,
                     if(!singular[2L]) 1 - 1e-6)
  quantiles <- tryCatch(
    suppressWarnings(quant(multiplier, probabilities)),
    error = function(e) numeric()
  )
  quartiles <- if(any(singular)){
    tryCatch(
      suppressWarnings(as.numeric(quant(multiplier, c(.25, .75)))),
      error = function(e) c(NA_real_, NA_real_)
    )
  }
  list(lower = lower, upper = upper, singular = singular,
       quantiles = as.numeric(quantiles), quartiles = quartiles)
}

.prior_conditional_normal_breakpoints <- function(spec, value, extra = numeric(),
                                                  setup = NULL){

  if(is.null(setup)){
    setup <- .prior_conditional_normal_breakpoint_setup(spec$multiplier, spec$bounds)
  }
  lower <- setup$lower
  upper <- setup$upper
  singular <- setup$singular
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
  quantiles <- setup$quantiles
  inner <- c(inner, quantiles)
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
    quartiles <- setup$quartiles
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

# The breakpoints of .prior_conditional_normal_breakpoints() for each of the
# single values 'values' with the same 'setup', computed for all values
# together (batched densities evaluate hundreds of values): element j of the
# returned list is identical to
# .prior_conditional_normal_breakpoints(spec, values[j], extra[j, ], setup)
# with the NA entries of row j of the matrix 'extra' left out. The points of
# all values are filtered, sorted and deduplicated in one pass, and the two
# sequential rules (the thinning next to a singular bound and the merging of
# close points) run over the values' points in the same order, one position
# at a time for all values.
.prior_conditional_normal_breakpoints_values <- function(spec, values, setup, extra = NULL){

  n <- length(values)
  if(n == 0L){
    return(list())
  }
  lower <- setup$lower
  upper <- setup$upper
  columns <- list()
  peak <- logical()
  peak_width <- rep(Inf, n)
  if(isTRUE(spec$product_mean != 0) &&
     isTRUE(spec$product_sd <= abs(spec$product_mean) / 2)){
    centre <- (values - spec$additive_mean) / spec$product_mean
    width <- sqrt(spec$additive_sd^2 + (spec$product_sd * centre)^2) / abs(spec$product_mean)
    offsets <- as.vector(outer(c(-1, 1), c(1, 3, 10)))
    columns <- c(list(centre), lapply(offsets, function(offset) centre + offset * width))
    peak <- rep(TRUE, 7L)
    peak_width <- ifelse(is.finite(width), width, Inf)
  }
  if(!is.null(extra)){
    columns <- c(columns, lapply(seq_len(ncol(extra)), function(j) extra[, j]))
    peak <- c(peak, rep(FALSE, ncol(extra)))
  }
  if(isTRUE(spec$additive_sd == 0) && isTRUE(spec$product_sd > 0)){
    distance <- abs(values - spec$additive_mean) / spec$product_sd
    valid <- is.finite(distance) & distance > 0
    scale_points <- lapply(c(.1, 1, 10), function(multiple){
      point <- distance * multiple
      point[!valid] <- NA_real_
      list(-1 * point, 1 * point)
    })
    columns <- c(columns, list(rep(0, n)), unlist(scale_points, recursive = FALSE))
    peak <- c(peak, rep(FALSE, 7L))
  }
  quantiles <- setup$quantiles
  columns <- c(columns, lapply(quantiles, rep, times = n))
  peak <- c(peak, rep(FALSE, length(quantiles)))

  # the candidate points of each value in the order of the single-value
  # function (value by value), without the absent (NA) extra points
  k <- length(columns)
  point <- if(k > 0L) as.vector(t(do.call(cbind, columns))) else numeric()
  index <- rep(seq_len(n), each = k)
  is_peak <- rep(peak, times = n)
  widths <- ifelse(is_peak, peak_width[index], Inf)
  inside <- is.finite(point) & point > lower & point < upper
  point <- point[inside]
  index <- index[inside]
  is_peak <- is_peak[inside]
  widths <- widths[inside]
  if(length(point) > 0L){
    density <- suppressWarnings(exp(lpdf(spec$multiplier, point)))
    finite <- is.finite(density)
    point <- point[finite]
    index <- index[finite]
    is_peak <- is_peak[finite]
    widths <- widths[finite]
  }
  # sorted and unique per value (both orders are stable, so a peak point
  # equal to a quantile stays a peak point)
  sorted <- order(index, point)
  point <- point[sorted]
  index <- index[sorted]
  is_peak <- is_peak[sorted]
  widths <- widths[sorted]
  m <- length(point)
  first <- if(m > 0L) c(TRUE, index[-1L] != index[-m] | point[-1L] != point[-m]) else logical()
  point <- point[first]
  index <- index[first]
  is_peak <- is_peak[first]
  widths <- widths[first]

  # the rank of each point in 'order' within its value, as a matrix of point
  # positions with one row per value
  positions <- function(order){
    counts <- tabulate(index, nbins = n)
    out <- matrix(NA_integer_, n, max(c(0L, counts)))
    if(length(order) > 0L){
      out[cbind(index[order], sequence(counts[unique(index[order])]))] <- order
    }
    out
  }

  if(any(setup$singular)){
    for(side in which(setup$singular)){
      bound <- c(lower, upper)[side]
      start <- abs(setup$quartiles[side] - bound)
      if(!isTRUE(is.finite(start)) || length(point) == 0L){
        next
      }
      distance <- abs(point - bound)
      keep <- rep(TRUE, length(point))
      last <- rep(start, n)
      visit <- positions(order(index, distance, decreasing = c(FALSE, TRUE), method = "radix"))
      for(rank in seq_len(ncol(visit))){
        i <- visit[, rank]
        i <- i[!is.na(i)]
        i <- i[distance[i] < start]
        exempt <- is_peak[i] & distance[i] >= widths[i]
        take <- exempt | distance[i] >= 1e-3 * last[index[i]]
        last[index[i][take]] <- distance[i][take]
        keep[i[!take]] <- FALSE
      }
      point <- point[keep]
      index <- index[keep]
      is_peak <- is_peak[keep]
      widths <- widths[keep]
    }
  }

  minimum_width <- function(a, b, peak_width){
    scale <- pmax(1, ifelse(is.finite(a), abs(a), 1), ifelse(is.finite(b), abs(b), 1))
    pmax(16 * .Machine$double.eps * scale, pmin(1e-9 * scale, peak_width / 2))
  }
  keep <- rep(FALSE, length(point))
  last <- rep(lower, n)
  last_kept <- rep(NA_integer_, n)
  visit <- positions(seq_along(point))
  for(rank in seq_len(ncol(visit))){
    i <- visit[, rank]
    i <- i[!is.na(i)]
    value_index <- index[i]
    take <- point[i] - last[value_index] >=
      minimum_width(last[value_index], point[i], peak_width[value_index])
    keep[i[take]] <- TRUE
    last[value_index[take]] <- point[i][take]
    last_kept[value_index[take]] <- i[take]
  }
  closing <- !is.na(last_kept)
  closing[closing] <- upper - point[last_kept[closing]] <
    minimum_width(point[last_kept[closing]], upper, peak_width[closing])
  keep[last_kept[closing]] <- FALSE

  point <- point[keep]
  index <- index[keep]
  unname(split(
    c(rep(lower, n), point, rep(upper, n)),
    factor(c(seq_len(n), index, seq_len(n)), levels = seq_len(n))
  ))
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
      standardized <- .prior_density_context_standardized_weights(density_context, weights,
                                                                  source_transforms)
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

  x <- .prior_linear_density_materialize(x)
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

  # A curve without a structural route that relies on a numerical grid with an
  # unresolved product component (a heavy-tailed factor; see
  # .prior_linear_density_route_product()) is omitted with a classed warning;
  # the atoms are still drawn.
  draw_curve <- !is.null(dist$density) && dist$density$mass > 0
  route <- NULL
  if(draw_curve){
    route <- .prior_density_route_from_adaptive(
      attr(dist, "adaptive_evaluation", exact = TRUE)
    )
    if(!identical(route$type, "unknown") && .prior_density_route_has_leaf(route, "unknown")){
      route <- .prior_density_route_with_grids(route)
    }
    unresolved <- if(identical(route$type, "unknown")){
      isFALSE(attr(dist, "product_grid_resolution", exact = TRUE)$resolved)
    }else{
      .prior_density_route_unresolved_products(route)
    }
    if(draw_curve && unresolved){
      .prior_linear_density_warn_curve_unavailable()
      draw_curve <- FALSE
    }
  }

  if(draw_curve){
    # The continuous density is evaluated on its structural route (closed
    # forms, and quadrature leaves by one batched quadrature over the plotted
    # values); only a combination without a structural route, or a density
    # without recorded provenance, interpolates its numerical grid. A route
    # with quadrature leaves is plotted on at most
    # .prior_linear_density_display_size() equally spaced values plus the
    # values it must include (support bounds and jumps, with a point just
    # outside each, offsets and other peaks, and atoms), and the vertex of a
    # parabola through each local maximum and its neighbours. Every curve
    # drops to zero at a finite support bound of the continuous density that
    # lies in the plotted range and has a positive density
    # (.prior_linear_density_edge_bounds()).
    structural <- !is.null(route) && !identical(route$type, "unknown")
    display <- structural && .prior_density_route_has_quadrature(route)
    size <- if(display) min(n_points, .prior_linear_density_display_size()) else n_points

    grid <- function(size){
      if(!is.null(transformed_x_range)){
        den <- seq(transformed_x_range[1], transformed_x_range[2], length.out = size)
        return(list(raw = suppressWarnings(.density.prior_transformation_inv_grid(
          den, transformation, transformation_arguments
        )), den = den))
      }
      limits <- if(is.null(x_range)) range(dist$density$x) else x_range
      list(raw = seq(limits[1], limits[2], length.out = size), den = NULL)
    }

    unresolved_values <- numeric()
    evaluate <- function(raw){
      finite_raw <- is.finite(raw)
      y <- rep(NA_real_, length(raw))
      if(any(finite_raw) && structural){
        y[finite_raw] <- .prior_density_route_density(route, raw[finite_raw])
        unresolved_values <<- unique(c(unresolved_values, raw[finite_raw & is.na(y)]))
      }else if(any(finite_raw)){
        y[finite_raw] <- stats::approx(
          dist$density$x,
          dist$density$y,
          xout   = raw[finite_raw],
          yleft  = 0,
          yright = 0
        )$y * dist$density$mass
      }
      y
    }
    plotted <- function(raw, den, y){
      if(is.null(transformation)){
        return(list(x = raw, y = y))
      }
      if(is.null(den)){
        den <- .density.prior_transformation_x(raw, transformation, transformation_arguments)
      }
      finite_raw <- is.finite(raw)
      y_transformed <- rep(NA_real_, length(y))
      y_transformed[finite_raw] <- .density.prior_transformation_y(
        den[finite_raw],
        y[finite_raw],
        transformation,
        transformation_arguments
      )
      # the image of the source value 0 under exp_lin takes the analytic
      # limit of density.prior() (.density.prior_transformation_grid()), so a
      # support edge at 0 is mapped as there
      at_zero <- finite_raw & raw == 0
      if(is.character(transformation) && identical(transformation, "exp_lin") && any(at_zero)){
        limit <- .density.prior_transformation_grid(
          raw[at_zero], y[at_zero], transformation, transformation_arguments
        )
        y_transformed[at_zero] <- ifelse(limit$drop, NA_real_, limit$y)
      }
      list(x = den, y = y_transformed)
    }

    # raw values in plotting order, with their plotted coordinates when the
    # plotted range is on the transformed scale
    points <- grid(size)
    edges  <- .prior_linear_density_edge_bounds(route, dist, structural, points$raw, evaluate)
    extra  <- unique(c(
      if(display) .prior_linear_density_display_values(route, dist, points$raw),
      edges$interior
    ))
    if(length(extra) > 0L){
      points <- grid(max(ceiling(size / 2), size - length(extra)))
      points <- .prior_linear_density_add_display_values(
        points, extra, transformation, transformation_arguments, transformed_x_range
      )
    }

    # a support bound with a positive density inside it repeats its value with
    # density 0 on its outer side (the vertical edge of density.prior()); a
    # mapped bound whose route density is 0 by rounding takes the ordinate's
    # limit inside the support
    y <- .prior_linear_density_edge_limits(points$raw, evaluate(points$raw), edges)
    repeats <- .prior_linear_density_edge_repeats(points$raw, edges)
    if(!is.null(repeats)){
      points$raw <- points$raw[repeats$index]
      if(!is.null(points$den)){
        points$den <- points$den[repeats$index]
      }
      y <- y[repeats$index]
      y[repeats$zero] <- 0
    }
    curve <- plotted(points$raw, points$den, y)
    if(display){
      curve <- .prior_linear_density_refine_peaks(
        curve, evaluate, plotted, transformation, transformation_arguments
      )
    }
    x_den <- curve$x
    y_den <- curve$y

    finite <- is.finite(x_den) & is.finite(y_den)
    x_den  <- x_den[finite]
    y_den  <- y_den[finite]
    if(length(unresolved_values) > 0L){
      .prior_linear_density_warn_curve_unavailable(
        unresolved_values = unresolved_values, partial = length(x_den) > 0L
      )
    }

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

# Warning of a partly or wholly unavailable prior curve.
.prior_linear_density_warn_curve_unavailable <- function(unresolved_values = NULL, partial = FALSE){

  warning(structure(
    class = c("BayesTools_prior_curve_unavailable", "BayesTools_plot_condition",
              "warning", "condition"),
    list(
      message = if(!is.null(unresolved_values)) paste0(
        "The prior density curve is ", if(partial) "partially unavailable" else "unavailable",
        ": numerical evaluations were unresolved at ", length(unresolved_values),
        if(length(unresolved_values) == 1L) " evaluation coordinate on the source scale. " else
          " evaluation coordinates on the source scale. ",
        if(partial) "Available curve points and declared atoms are retained; " else
          "No continuous curve points are available; declared atoms are retained; ",
        "use 'prior = FALSE' to draw the posterior alone."
      ) else paste0(
        "The prior density curve is unavailable: this prior-density ",
        "combination has no exact route, and its numerical grid cannot ",
        "resolve the scale of a heavy-tailed product term (e.g. a Cauchy ",
        "ordered total or 'multiply_by' factor). The prior curve is omitted ",
        "from the plot; plot the terms separately."
      ),
      call = NULL,
      unresolved_values = unresolved_values
    )
  ))
}

# Number of equally spaced values of a plotted density with quadrature leaves.
.prior_linear_density_display_size <- function(){

  200L
}

# Values within the range of the raw plotting values 'raw' that a plotted
# density of 'route' must include: its display points, its support bounds and
# the atoms of 'dist'. The density at a bound is the one-sided limit inside
# the support, so the value 1e-6 of the plotted range outside the bound (below
# a lower bound, above an upper bound) draws the jump; the value as far
# inside it would repeat the value at the bound and is added only where the
# density at the bound is numerically unavailable. An infinite density stays
# off the plotted grid without an extra nearby point whose arbitrary height
# would dominate the display scale.
.prior_linear_density_display_values <- function(route, dist, raw){

  raw <- raw[is.finite(raw)]
  if(length(raw) == 0L){
    return(numeric())
  }
  limits <- range(raw)
  special <- .prior_density_route_display_points(route)
  atoms <- if(!is.null(dist$points) && nrow(dist$points) > 0L){
    dist$points$x[dist$points$p > 0]
  }
  delta <- 1e-6 * diff(limits)
  bounds <- c(special$lower, special$upper)
  inside <- c(special$lower + delta, special$upper - delta)
  undrawn <- inside >= limits[1L] & inside <= limits[2L]
  if(any(undrawn)){
    undrawn[undrawn] <- is.na(.prior_density_route_density(route, bounds[undrawn]))
  }
  values <- c(special$points, atoms, bounds,
              special$lower - delta, special$upper + delta, inside[undrawn])
  unique(values[is.finite(values) & values >= limits[1L] & values <= limits[2L]])
}

# Finite support bounds of the continuous density at which a plotted curve
# drops to zero (the edge of density.prior()): bounds of the structural route
# (or, for a curve interpolated from a numerical grid, the exact support hull of
# its provenance) within the range of the raw plotting values 'raw' (up to a
# relative 1e-9 of it, where the range of a numerical grid ends at the bound
# up to rounding), with a positive finite density as 'evaluate'
# returns it and exactly zero density 1e-6 of the range outside the bound. A
# bound at which the density is zero or unavailable, that lies outside the
# range, or inside another component's support (an interior jump to a positive
# density, drawn by the display values) is not an edge. 'interior' are the
# edge bounds strictly inside the range, which the plotted values must include;
# 'tolerance' matches a bound to the nearest plotted value. A mapped bound
# whose inverse image rounds just outside the source support has the route
# density 0; the density there is the one-sided limit inside the support that
# the ordinate takes (.prior_density_ordinate_endpoint_matches()), so such a
# bound is an edge when that limit is positive, and 'snap' lists these bounds
# and their limits, the plotted values at them.
.prior_linear_density_edge_bounds <- function(route, dist, structural, raw, evaluate){

  none <- list(lower = numeric(), upper = numeric(), interior = numeric(), tolerance = 0,
               snap = list(bound = numeric(), y = numeric()))
  raw <- raw[is.finite(raw)]
  if(length(raw) < 2L){
    return(none)
  }
  limits <- range(raw)
  width <- diff(limits)
  if(!is.finite(width) || width <= 0){
    return(none)
  }

  candidates <- if(structural){
    bounds <- .prior_density_route_display_points(route)
    list(lower = bounds$lower, upper = bounds$upper)
  }else{
    hull <- .prior_linear_density_support_hull(attr(dist, "adaptive_evaluation", exact = TRUE))
    list(lower = hull[1L], upper = hull[2L])
  }
  tolerance <- 1e-9 * width
  density <- function(values){
    out <- tryCatch(suppressWarnings(evaluate(values)), error = function(e) NULL)
    if(length(out) == length(values)) out else rep(NA_real_, length(values))
  }
  # the one-sided limit inside the support at bounds where the route density
  # is 0: the continuous density of the route's ordinate
  limit <- function(bounds){
    vapply(bounds, function(bound){
      if(!structural){
        return(NA_real_)
      }
      tryCatch(
        suppressWarnings(.prior_density_ordinate_height_value(
          .prior_density_route_ordinate(route, bound)
        )),
        error = function(e) NA_real_
      )
    }, numeric(1))
  }
  edge <- function(bounds, side){
    bounds <- unique(bounds[is.finite(bounds) &
                              bounds >= limits[1L] - tolerance & bounds <= limits[2L] + tolerance])
    if(length(bounds) == 0L){
      return(list(bounds = numeric(), y = numeric(), snapped = logical()))
    }
    inside  <- density(bounds)
    outside <- density(bounds + side * 1e-6 * width)
    snapped <- !is.na(inside) & inside == 0
    if(any(snapped)){
      inside[snapped] <- limit(bounds[snapped])
    }
    keep <- is.finite(inside) & inside > 0 & !is.na(outside) & outside == 0
    list(bounds = bounds[keep], y = inside[keep], snapped = snapped[keep])
  }

  lower <- edge(candidates$lower, -1)
  upper <- edge(candidates$upper, 1)
  bounds  <- c(lower$bounds, upper$bounds)
  snapped <- c(lower$snapped, upper$snapped)
  list(lower = lower$bounds, upper = upper$bounds,
       interior = bounds[bounds > limits[1L] & bounds < limits[2L]],
       tolerance = tolerance,
       snap = list(bound = bounds[snapped], y = c(lower$y, upper$y)[snapped]))
}

# The plotted values 'y' at the raw values 'raw' with the one-sided limit
# inside the support (edges$snap, .prior_linear_density_edge_bounds()) at each
# bound whose own route density is 0.
.prior_linear_density_edge_limits <- function(raw, y, edges){

  snap <- edges$snap
  if(length(snap$bound) == 0L){
    return(y)
  }
  positions <- .prior_linear_density_edge_positions(raw, snap$bound, edges$tolerance)
  use <- which(!is.na(positions))
  use <- use[vapply(use, function(i) isTRUE(y[positions[i]] == 0), logical(1))]
  y[positions[use]] <- snap$y[use]
  y
}

# The position in 'raw' of the finite value nearest to each of 'values' within
# 'tolerance' (NA where there is none).
.prior_linear_density_edge_positions <- function(raw, values, tolerance){

  finite <- which(is.finite(raw))
  vapply(values, function(value){
    distance <- abs(raw[finite] - value)
    nearest  <- which.min(distance)
    if(length(nearest) == 1L && distance[nearest] <= tolerance) finite[nearest] else NA_integer_
  }, integer(1))
}

# The plotted raw values 'raw' with the plotted value of each edge bound of
# 'edges' (.prior_linear_density_edge_bounds()) repeated: the indices of the
# repeated vector, and the positions of the repeats drawn at density 0, each on
# the outer side of its bound (before the bound's own value where the raw
# values increase, after it where they decrease). NULL without an edge.
.prior_linear_density_edge_repeats <- function(raw, edges){

  finite <- which(is.finite(raw))
  if(length(finite) < 2L || length(c(edges$lower, edges$upper)) == 0L){
    return(NULL)
  }
  ascending <- raw[finite[length(finite)]] > raw[finite[1L]]
  locate <- function(values){
    positions <- .prior_linear_density_edge_positions(raw, values, edges$tolerance)
    unique(positions[!is.na(positions)])
  }
  lower  <- locate(edges$lower)
  upper  <- locate(edges$upper)
  before <- if(ascending) lower else upper
  after  <- setdiff(if(ascending) upper else lower, before)
  if(length(before) + length(after) == 0L){
    return(NULL)
  }
  copies <- rep(1L, length(raw))
  copies[c(before, after)] <- 2L
  first <- cumsum(c(1L, copies))[seq_along(raw)]
  list(index = rep(seq_along(raw), copies),
       zero  = c(first[before], first[after] + 1L))
}

# Plotting values 'points' (raw values and, for a range on the transformed
# scale, their plotted coordinates) with the raw values 'extra' added, in
# plotting order.
.prior_linear_density_add_display_values <- function(points, extra, transformation,
                                                     transformation_arguments,
                                                     transformed_x_range){

  if(is.null(points$den)){
    return(list(raw = sort(unique(c(points$raw, extra))), den = NULL))
  }
  den_extra <- suppressWarnings(.density.prior_transformation_x(
    extra, transformation, transformation_arguments
  ))
  keep <- is.finite(den_extra) &
    den_extra >= min(transformed_x_range) & den_extra <= max(transformed_x_range)
  den <- c(points$den, den_extra[keep])
  raw <- c(points$raw, extra[keep])
  ordered <- order(den)
  unique_den <- !duplicated(den[ordered])
  list(raw = raw[ordered][unique_den], den = den[ordered][unique_den])
}

# The plotted curve with the vertex of the parabola through each interior local
# maximum and its two neighbours added (in plotted coordinates), evaluated by
# 'evaluate' (raw values) and mapped by 'plotted'.
.prior_linear_density_refine_peaks <- function(curve, evaluate, plotted,
                                               transformation, transformation_arguments){

  x <- curve$x
  y <- curve$y
  n <- length(x)
  if(n < 3L){
    return(curve)
  }
  i <- seq(2L, n - 1L)
  # a value repeated at a support edge has no parabola
  peak <- i[is.finite(x[i - 1L]) & is.finite(x[i + 1L]) &
              is.finite(y[i - 1L]) & is.finite(y[i]) & is.finite(y[i + 1L]) &
              x[i] != x[i - 1L] & x[i] != x[i + 1L] &
              y[i] > y[i - 1L] & y[i] >= y[i + 1L]]
  if(length(peak) == 0L){
    return(curve)
  }
  x0 <- x[peak - 1L]
  x1 <- x[peak]
  x2 <- x[peak + 1L]
  y0 <- y[peak - 1L]
  y1 <- y[peak]
  y2 <- y[peak + 1L]
  numerator <- (x1 - x0)^2 * (y1 - y2) - (x1 - x2)^2 * (y1 - y0)
  denominator <- (x1 - x0) * (y1 - y2) - (x1 - x2) * (y1 - y0)
  vertex <- x1 - numerator / (2 * denominator)
  vertex <- unique(vertex[is.finite(vertex) & vertex > x0 & vertex < x2 & vertex != x1])
  if(length(vertex) == 0L){
    return(curve)
  }
  raw <- if(is.null(transformation)) vertex else suppressWarnings(
    .density.prior_transformation_inv_grid(vertex, transformation, transformation_arguments)
  )
  added <- plotted(raw, if(is.null(transformation)) NULL else vertex, evaluate(raw))
  x <- c(x, added$x)
  y <- c(y, added$y)
  ordered <- order(x)
  list(x = x[ordered], y = y[ordered])
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
      .prior_density_context_standardized_weights(context, weights, source_transforms),
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
    "allocation_product" = .prior_allocation_product_hull(arguments),
    # the support of the transformed source, mapped below
    "output_transformation" = .prior_linear_density_support_hull(arguments$source),
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
      }else if(identical(integration$kind, "convolution")){
        "Convolution prior density"
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
# 1e-4 * H (purely relative), so changes cannot cancel between components
# and a component whose own density is about zero at the value does not block
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
      .prior_density_context_standardized_weights(component_context, weights, source_transforms),
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
  # grid components outside their support hull contribute exactly zero
  for(i in seq_along(terms)){
    if(!identical(terms[[i]]$method, "grid")){
      next
    }
    support <- .prior_linear_density_support_hull(
      attr(terms[[i]]$density, "adaptive_evaluation", exact = TRUE)
    )
    if(!is.null(support) && (value < support[1L] || value > support[2L])){
      terms[[i]] <- list(weight = terms[[i]]$weight, method = "exact", height = 0)
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
  refinement <- .prior_linear_density_refine_grids(
    densities, weights[grid],
    evaluate = function(density) .prior_linear_density_grid_height(density, value),
    covers   = function(density) .prior_linear_density_covers(density, value),
    fixed    = fixed
  )
  if(!isTRUE(refinement$converged)){
    .prior_linear_density_stop_refinement(refinement, value)
  }
  height <- refinement$total
  attr(height, "adaptive_evaluation") <- list(
    components      = sum(grid),
    refinements     = refinement$refinements,
    absolute_change = refinement$absolute_change,
    error_bound     = refinement$error_bound,
    converged       = TRUE
  )
  attr(height, "component_heights") <- diagnostics(refinement$values)
  height
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
    # the component and grid heights read the grid and its resolution
    x <- .prior_linear_density_materialize(x)
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


  if(length(value) != 1L || !is.finite(value)){
    stop("Adaptive prior-density evaluation requires one finite ordinate.",
         call. = FALSE)
  }

  # A numerical grid without deterministic provenance (the prior measure it
  # approximates) has no error control and cannot be refined: it is used for
  # plots only. Without a continuous part the continuous height is zero.
  if(is.null(attr(x, "adaptive_evaluation", exact = TRUE))){
    if(is.null(x$density) || !isTRUE(x$density$mass > 0)){
      return(0)
    }
    .prior_linear_density_stop_no_provenance("height")
  }

  support <- .prior_linear_density_support_hull(
    attr(x, "adaptive_evaluation", exact = TRUE)
  )
  if(!is.null(support) && (value < support[1L] || value > support[2L])){
    return(0)
  }

  refinement <- .prior_linear_density_refine_grids(
    list(x), 1,
    evaluate = function(density) .prior_linear_density_grid_height(density, value),
    covers   = function(density) .prior_linear_density_covers(density, value)
  )
  if(!isTRUE(refinement$converged)){
    .prior_linear_density_stop_refinement(refinement, value)
  }
  height <- refinement$total
  refined <- refinement$densities[[1L]]
  attr(height, "adaptive_evaluation") <- c(
    attr(refined, "refinement_settings", exact = TRUE),
    list(
      refinements = refinement$refinements,
      absolute_change = refinement$absolute_change,
      error_bound = refinement$error_bound,
      numerical_range = .prior_linear_density_range(refined),
      converged = TRUE
    )
  )
  height
}

.prior_linear_density_stop_no_provenance <- function(quantity){

  stop(
    "The prior density has no deterministic provenance, so its ", quantity,
    " is unavailable: a numerical density grid without provenance is used ",
    "only for plotting. Construct the density with the BayesTools ",
    "prior-density builders.",
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
