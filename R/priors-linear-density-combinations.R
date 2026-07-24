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
  weights <- weights[abs(weights) > .prior_linear_density_zero_tol()]
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
    present <- present[abs(weights[present]) > .prior_linear_density_zero_tol()]
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

  unmatched <- names(weights)[abs(weights) > .prior_linear_density_zero_tol() & !matched]
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

  if(abs(weight) <= .prior_linear_density_zero_tol() || is.prior.none(prior)){
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

  if(norm_weight <= .prior_linear_density_zero_tol()){
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

  if(abs(weight) <= .prior_linear_density_zero_tol() || is.prior.none(prior)){
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
  if(isTRUE(all.equal(x_range[1], x_range[2]))){
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

  y[!is.finite(y)] <- 0
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

    keep <- is.finite(x_new) & is.finite(y_new)
    x_new <- x_new[keep]
    y_new <- y_new[keep]

    if(length(x_new) >= 2){
      ord <- order(x_new)
      x_new <- x_new[ord]
      y_new <- y_new[ord]

      keep <- c(TRUE, diff(x_new) > .prior_linear_density_zero_tol() * pmax(1, abs(x_new[-length(x_new)])))
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
    points <- points[is.finite(points$x), , drop = FALSE]
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
  .prior_linear_density_normalize(out)
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
    present <- present[abs(additive_weights[present]) > .prior_linear_density_zero_tol()]
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
  weights <- weights[abs(weights) > .prior_linear_density_zero_tol()]

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
                                              output_transformation_arguments = NULL){

  check_list(prior_list, "prior_list")
  if(is.null(weights) || length(weights) == 0){
    weights <- numeric()
  }else{
    check_real(weights, "weights", check_length = 0)
  }
  check_int(n_grid, "n_grid", lower = 16)
  check_real(tail_prob, "tail_prob", lower = 0, upper = 0.5, allow_bound = FALSE)

  weights <- weights[is.finite(weights)]
  weights <- weights[abs(weights) > .prior_linear_density_zero_tol()]

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
  return(dist)
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

  if(!is.null(dist$density) && dist$density$mass > .prior_linear_density_zero_tol()){
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
    points <- points[points$p > .prior_linear_density_zero_tol(), , drop = FALSE]
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

  height <- 0
  if(!is.null(x$density) && value >= min(x$density$x) && value <= max(x$density$x)){
    height <- height + stats::approx(
      x$density$x,
      x$density$y * x$density$mass,
      xout = value,
      yleft = 0,
      yright = 0
    )$y
  }

  height
}

.prior_linear_density_point_mass <- function(x, value){

  if(!inherits(x, "prior_linear_density") || is.null(x$points) || nrow(x$points) == 0){
    return(0)
  }

  sum(x$points$p[abs(x$points$x - value) <= .prior_linear_density_zero_tol()])
}
