.prior_linear_density_default_grid <- function(){
  4096L
}

.prior_linear_density_tail_prob <- function(){
  1e-4
}

.prior_linear_density_product_grid <- function(){
  1024L
}

.prior_linear_density_zero_tol <- function(){
  sqrt(.Machine$double.eps)
}

.prior_linear_density_empty_points <- function(){
  data.frame(x = numeric(), p = numeric())
}

.prior_linear_density_point <- function(x, p = 1){
  out <- list(
    density = NULL,
    points  = data.frame(x = x, p = p),
    n_grid  = 1L
  )
  class(out) <- c("prior_linear_density", "prior_density")
  return(.prior_linear_density_normalize(out))
}

.prior_linear_density_normalize <- function(dist){

  point_mass <- if(!is.null(dist$points) && nrow(dist$points) > 0) sum(dist$points$p) else 0
  density_mass <- if(!is.null(dist$density)) dist$density$mass else 0
  total_mass <- point_mass + density_mass

  if(!is.finite(total_mass) || total_mass <= 0){
    stop("The computed prior density has zero total mass.", call. = FALSE)
  }

  if(!is.null(dist$density)){
    dist$density$mass <- dist$density$mass / total_mass
  }
  if(!is.null(dist$points) && nrow(dist$points) > 0){
    dist$points$p <- dist$points$p / total_mass
  }

  class(dist) <- unique(c(class(dist), "prior_linear_density", "prior_density"))
  return(dist)
}

.prior_linear_density_fft_convolve <- function(a, b){

  n <- length(a) + length(b) - 1L
  n_fft <- 2^ceiling(log2(n))
  out <- fft(
    fft(c(a, rep(0, n_fft - length(a)))) *
      fft(c(b, rep(0, n_fft - length(b)))),
    inverse = TRUE
  ) / n_fft

  pmax(Re(out[seq_len(n)]), 0)
}

.prior_linear_density_aggregate_points <- function(points, dx){

  if(is.null(points) || nrow(points) == 0){
    return(.prior_linear_density_empty_points())
  }

  points <- points[is.finite(points$x) & is.finite(points$p) & points$p > 0, , drop = FALSE]
  if(nrow(points) == 0){
    return(.prior_linear_density_empty_points())
  }

  if(is.finite(dx) && dx > 0){
    key <- as.character(round(points$x / dx))
  }else{
    key <- format(signif(points$x, 14), scientific = FALSE)
  }

  split_points <- split(points, key)
  out <- do.call(rbind, lapply(split_points, function(p){
    data.frame(
      x = stats::weighted.mean(p$x, p$p),
      p = sum(p$p)
    )
  }))
  rownames(out) <- NULL
  out <- out[order(out$x), , drop = FALSE]

  return(out)
}

.prior_linear_density_coalesce <- function(densities = list(), points = NULL, dx, n_grid = NULL){

  points <- .prior_linear_density_aggregate_points(points, dx)

  densities <- densities[!vapply(densities, is.null, logical(1))]
  densities <- densities[vapply(densities, function(d){
    !is.null(d$x) && !is.null(d$y) && is.finite(d$mass) && d$mass > 0
  }, logical(1))]

  density <- NULL
  if(length(densities) > 0){
    x_min <- min(vapply(densities, function(d) min(d$x), numeric(1)))
    x_max <- max(vapply(densities, function(d) max(d$x), numeric(1)))

    if(!is.finite(dx) || dx <= 0){
      if(is.null(n_grid)){
        n_grid <- .prior_linear_density_default_grid()
      }
      dx <- (x_max - x_min) / max(1, n_grid - 1)
    }

    if(!is.finite(dx) || dx <= 0 || isTRUE(all.equal(x_min, x_max))){
      x <- x_min
      y_mass <- sum(vapply(densities, function(d) d$mass, numeric(1)))
    }else{
      x <- seq(x_min, x_max, by = dx)
      if(length(x) < 2){
        x <- seq(x_min, x_max, length.out = 2)
      }
      y_mass <- numeric(length(x))
      for(d in densities){
        y_mass <- y_mass + d$mass * stats::approx(d$x, d$y, xout = x, yleft = 0, yright = 0)$y
      }
    }

    density_mass <- sum(vapply(densities, function(d) d$mass, numeric(1)))
    area <- if(length(x) > 1) sum(y_mass) * (x[2] - x[1]) else density_mass
    if(is.finite(area) && area > 0 && density_mass > 0){
      y <- y_mass / area
    }else{
      y <- y_mass
    }

    density <- list(
      x    = x,
      y    = y,
      mass = density_mass
    )
  }

  out <- list(
    density = density,
    points  = points,
    n_grid  = if(is.null(n_grid)) length(if(!is.null(density)) density$x else points$x) else n_grid
  )
  class(out) <- c("prior_linear_density", "prior_density")
  return(.prior_linear_density_normalize(out))
}

.prior_linear_density_convolve <- function(lhs, rhs, dx){

  densities <- list()
  points <- .prior_linear_density_empty_points()

  if(!is.null(lhs$points) && nrow(lhs$points) > 0 && !is.null(rhs$points) && nrow(rhs$points) > 0){
    point_grid <- merge(lhs$points, rhs$points, by = NULL)
    points <- rbind(points, data.frame(
      x = point_grid$x.x + point_grid$x.y,
      p = point_grid$p.x * point_grid$p.y
    ))
  }

  if(!is.null(lhs$density) && !is.null(rhs$density)){
    y <- .prior_linear_density_fft_convolve(lhs$density$y, rhs$density$y) * dx
    area <- sum(y) * dx
    if(is.finite(area) && area > 0){
      y <- y / area
    }
    x <- lhs$density$x[1] + rhs$density$x[1] + dx * (seq_along(y) - 1)
    densities[[length(densities) + 1L]] <- list(
      x    = x,
      y    = y,
      mass = lhs$density$mass * rhs$density$mass
    )
  }

  if(!is.null(lhs$density) && !is.null(rhs$points) && nrow(rhs$points) > 0){
    for(i in seq_len(nrow(rhs$points))){
      densities[[length(densities) + 1L]] <- list(
        x    = lhs$density$x + rhs$points$x[i],
        y    = lhs$density$y,
        mass = lhs$density$mass * rhs$points$p[i]
      )
    }
  }

  if(!is.null(rhs$density) && !is.null(lhs$points) && nrow(lhs$points) > 0){
    for(i in seq_len(nrow(lhs$points))){
      densities[[length(densities) + 1L]] <- list(
        x    = rhs$density$x + lhs$points$x[i],
        y    = rhs$density$y,
        mass = rhs$density$mass * lhs$points$p[i]
      )
    }
  }

  .prior_linear_density_coalesce(densities = densities, points = points, dx = dx,
                                 n_grid = max(lhs$n_grid, rhs$n_grid))
}

.prior_linear_density_mix <- function(dists, weights, dx, n_grid = NULL){

  weights <- weights / sum(weights)
  densities <- list()
  points <- .prior_linear_density_empty_points()

  for(i in seq_along(dists)){
    dist <- dists[[i]]
    w <- weights[i]

    if(!is.null(dist$density)){
      densities[[length(densities) + 1L]] <- list(
        x    = dist$density$x,
        y    = dist$density$y,
        mass = dist$density$mass * w
      )
    }
    if(!is.null(dist$points) && nrow(dist$points) > 0){
      temp_points <- dist$points
      temp_points$p <- temp_points$p * w
      points <- rbind(points, temp_points)
    }
  }

  .prior_linear_density_coalesce(densities = densities, points = points, dx = dx, n_grid = n_grid)
}

.prior_linear_density_range <- function(dist){

  values <- numeric()
  if(!is.null(dist$density) && length(dist$density$x) > 0){
    values <- c(values, range(dist$density$x, finite = TRUE))
  }
  if(!is.null(dist$points) && nrow(dist$points) > 0){
    values <- c(values, dist$points$x[dist$points$p > .prior_linear_density_zero_tol()])
  }
  values <- values[is.finite(values)]

  if(length(values) == 0){
    return(c(0, 0))
  }
  range(values)
}

.prior_linear_density_dx <- function(dist){

  if(!is.null(dist$density) && length(dist$density$x) > 1){
    dx <- median(diff(dist$density$x))
    if(is.finite(dx) && dx > 0){
      return(dx)
    }
  }
  NA_real_
}

.prior_linear_density_regrid <- function(dist, dx, n_grid = NULL){

  densities <- list()
  if(!is.null(dist$density)){
    densities[[1]] <- dist$density
  }

  .prior_linear_density_coalesce(
    densities = densities,
    points    = dist$points,
    dx        = dx,
    n_grid    = if(is.null(n_grid)) dist$n_grid else n_grid
  )
}

.prior_linear_density_sum_independent <- function(dists, n_grid = NULL){

  dists <- dists[!vapply(dists, is.null, logical(1))]
  if(length(dists) == 0){
    return(.prior_linear_density_point(0))
  }
  if(length(dists) == 1){
    return(dists[[1]])
  }

  ranges <- do.call(rbind, lapply(dists, .prior_linear_density_range))
  target_width <- sum(pmax(0, ranges[, 2] - ranges[, 1]))
  if(is.null(n_grid)){
    n_grid <- max(vapply(dists, function(dist) dist$n_grid, integer(1)))
  }
  dx <- target_width / max(1, n_grid - 1)
  if(!is.finite(dx) || dx <= 0){
    dx_values <- vapply(dists, .prior_linear_density_dx, numeric(1))
    dx_values <- dx_values[is.finite(dx_values) & dx_values > 0]
    dx <- if(length(dx_values) > 0) min(dx_values) else 1
  }
  if(!is.finite(dx) || dx <= 0){
    dx <- 1
  }

  dist <- .prior_linear_density_point(0)
  for(component in dists){
    dist <- .prior_linear_density_convolve(
      .prior_linear_density_regrid(dist, dx, n_grid),
      .prior_linear_density_regrid(component, dx, n_grid),
      dx
    )
  }

  dist
}

.prior_linear_density_scaled <- function(dist, scale, mass = 1, dx = NA_real_, n_grid = NULL){

  if(abs(scale) <= .prior_linear_density_zero_tol()){
    out <- .prior_linear_density_point(0)
    out$points$p <- mass
    return(out)
  }

  densities <- list()
  if(!is.null(dist$density) && dist$density$mass > 0){
    x <- dist$density$x * scale
    y <- dist$density$y / abs(scale)
    ord <- order(x)
    densities[[1]] <- list(
      x    = x[ord],
      y    = y[ord],
      mass = dist$density$mass
    )
  }

  points <- .prior_linear_density_empty_points()
  if(!is.null(dist$points) && nrow(dist$points) > 0){
    points <- data.frame(
      x = dist$points$x * scale,
      p = dist$points$p
    )
  }

  out <- .prior_linear_density_coalesce(
    densities = densities,
    points    = points,
    dx        = dx,
    n_grid    = if(is.null(n_grid)) dist$n_grid else n_grid
  )
  if(!is.null(out$density)){
    out$density$mass <- out$density$mass * mass
  }
  if(!is.null(out$points) && nrow(out$points) > 0){
    out$points$p <- out$points$p * mass
  }
  out
}

.prior_linear_density_product_range <- function(lhs, rhs){

  lhs_range <- .prior_linear_density_range(lhs)
  rhs_range <- .prior_linear_density_range(rhs)
  products <- as.vector(outer(lhs_range, rhs_range, `*`))

  lhs_points <- if(!is.null(lhs$points) && nrow(lhs$points) > 0) lhs$points$x else numeric()
  rhs_points <- if(!is.null(rhs$points) && nrow(rhs$points) > 0) rhs$points$x else numeric()
  if(length(lhs_points) > 0){
    products <- c(products, as.vector(outer(lhs_points, rhs_range, `*`)))
  }
  if(length(rhs_points) > 0){
    products <- c(products, as.vector(outer(lhs_range, rhs_points, `*`)))
  }
  if(length(lhs_points) > 0 && length(rhs_points) > 0){
    products <- c(products, as.vector(outer(lhs_points, rhs_points, `*`)))
  }

  products <- products[is.finite(products)]
  if(length(products) == 0){
    return(c(0, 0))
  }
  range(products)
}

.prior_linear_density_product <- function(lhs, rhs, n_grid = NULL){

  if(is.null(n_grid)){
    n_grid <- max(lhs$n_grid, rhs$n_grid)
  }
  n_grid <- min(max(16L, n_grid), .prior_linear_density_product_grid())

  densities <- list()
  points <- .prior_linear_density_empty_points()
  dx <- NA_real_

  if(!is.null(lhs$points) && nrow(lhs$points) > 0 && !is.null(rhs$points) && nrow(rhs$points) > 0){
    point_grid <- merge(lhs$points, rhs$points, by = NULL)
    points <- rbind(points, data.frame(
      x = point_grid$x.x * point_grid$x.y,
      p = point_grid$p.x * point_grid$p.y
    ))
  }

  if(!is.null(lhs$density) && !is.null(rhs$points) && nrow(rhs$points) > 0){
    for(i in seq_len(nrow(rhs$points))){
      scaled <- .prior_linear_density_scaled(lhs, rhs$points$x[i],
                                             mass = rhs$points$p[i],
                                             n_grid = n_grid)
      if(!is.null(scaled$density)){
        densities[[length(densities) + 1L]] <- scaled$density
      }
      if(!is.null(scaled$points) && nrow(scaled$points) > 0){
        points <- rbind(points, scaled$points)
      }
    }
  }

  if(!is.null(rhs$density) && !is.null(lhs$points) && nrow(lhs$points) > 0){
    for(i in seq_len(nrow(lhs$points))){
      scaled <- .prior_linear_density_scaled(rhs, lhs$points$x[i],
                                             mass = lhs$points$p[i],
                                             n_grid = n_grid)
      if(!is.null(scaled$density)){
        densities[[length(densities) + 1L]] <- scaled$density
      }
      if(!is.null(scaled$points) && nrow(scaled$points) > 0){
        points <- rbind(points, scaled$points)
      }
    }
  }

  if(!is.null(lhs$density) && !is.null(rhs$density)){
    z_range <- .prior_linear_density_product_range(lhs, rhs)
    if(isTRUE(all.equal(z_range[1], z_range[2]))){
      points <- rbind(points, data.frame(
        x = z_range[1],
        p = lhs$density$mass * rhs$density$mass
      ))
    }else{
      z <- seq(z_range[1], z_range[2], length.out = n_grid)
      dx <- z[2] - z[1]

      x <- lhs$density$x
      fx <- lhs$density$y
      dx_x <- .prior_linear_density_dx(lhs)
      zero_tol <- .prior_linear_density_zero_tol() * max(1, max(abs(x), na.rm = TRUE))
      keep <- is.finite(x) & is.finite(fx) & fx > 0 & abs(x) > zero_tol
      x <- x[keep]
      fx <- fx[keep]

      y <- numeric(length(z))
      if(length(x) > 0 && is.finite(dx_x) && dx_x > 0){
        integrand_weight <- fx / abs(x)
        chunk_size <- max(1L, floor(5e6 / length(x)))
        for(start in seq(1L, length(z), by = chunk_size)){
          end <- min(length(z), start + chunk_size - 1L)
          source <- outer(1 / x, z[start:end], `*`)
          fy <- stats::approx(
            rhs$density$x,
            rhs$density$y,
            xout   = as.vector(source),
            yleft  = 0,
            yright = 0
          )$y
          fy <- matrix(fy, nrow = length(x))
          y[start:end] <- colSums(fy * integrand_weight) * dx_x
        }
      }

      y[!is.finite(y)] <- 0
      area <- sum(y) * dx
      if(is.finite(area) && area > 0){
        y <- y / area
        densities[[length(densities) + 1L]] <- list(
          x    = z,
          y    = y,
          mass = lhs$density$mass * rhs$density$mass
        )
      }
    }
  }

  if(is.na(dx)){
    dx_values <- vapply(densities, function(d){
      if(length(d$x) > 1) median(diff(d$x)) else NA_real_
    }, numeric(1))
    dx_values <- dx_values[is.finite(dx_values) & dx_values > 0]
    dx <- if(length(dx_values) > 0) min(dx_values) else NA_real_
  }
  if(!is.finite(dx) || dx <= 0){
    dx <- NA_real_
  }

  .prior_linear_density_coalesce(
    densities = densities,
    points    = points,
    dx        = dx,
    n_grid    = n_grid
  )
}

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

  if(length(conditional) == 0){
    return(character())
  }

  active <- .prior_linear_active_parameters(prior_list, weights)
  conditional[conditional %in% active]
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

.prior_density_context <- function(prior_list, column_names, formula_scale = NULL,
                                   n_grid = .prior_linear_density_default_grid(),
                                   tail_prob = .prior_linear_density_tail_prob()){

  check_list(prior_list, "prior_list")
  check_char(column_names, "column_names", check_length = FALSE)
  check_list(formula_scale, "formula_scale", allow_NULL = TRUE)
  check_int(n_grid, "n_grid", lower = 16)
  check_real(tail_prob, "tail_prob", lower = 0, upper = 0.5, allow_bound = FALSE)

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    .check_formula_scale_info(formula_scale)
  }

  transforms <- list()
  if(!is.null(formula_scale) && length(formula_scale) > 0){
    for(param_name in names(formula_scale)){
      affected_cols <- grep(paste0("^", param_name, "_"), column_names, value = TRUE)
      if(length(affected_cols) == 0){
        next
      }
      transforms[[param_name]] <- list(
        columns       = affected_cols,
        matrix        = .build_unscale_matrix(affected_cols, formula_scale[[param_name]], param_name),
        log_intercept = isTRUE(attr(formula_scale[[param_name]], "log_intercept")),
        intercept     = paste0(param_name, "_intercept")
      )
    }
  }

  out <- list(
    prior_list    = prior_list,
    column_names  = column_names,
    formula_scale = formula_scale,
    transforms    = transforms,
    n_grid        = n_grid,
    tail_prob     = tail_prob
  )
  class(out) <- "prior_density_context"
  return(out)
}

.prior_density_context_standardized_weights <- function(context, weights){

  if(!inherits(context, "prior_density_context")){
    stop("'context' must be a prior density context.", call. = FALSE)
  }
  if(is.null(names(weights))){
    stop("'weights' must be a named numeric vector.", call. = FALSE)
  }

  weights <- weights[intersect(names(weights), context$column_names)]
  out <- rep(0, length(context$column_names))
  names(out) <- context$column_names

  used <- rep(FALSE, length(weights))
  names(used) <- names(weights)

  for(transform in context$transforms){
    cols <- intersect(transform$columns, names(weights))
    if(length(cols) == 0){
      next
    }

    if(transform$log_intercept &&
       transform$intercept %in% cols &&
       abs(weights[[transform$intercept]]) > .prior_linear_density_zero_tol()){
      stop(
        "Linear-combination prior densities with log-intercept scaling are only available ",
        "for the transformed intercept coefficient itself.",
        call. = FALSE
      )
    }

    temp_weights <- rep(0, length(transform$columns))
    names(temp_weights) <- transform$columns
    temp_weights[cols] <- weights[cols]
    out[transform$columns] <- out[transform$columns] +
      as.numeric(temp_weights %*% transform$matrix)
    used[cols] <- TRUE
  }

  remaining <- names(weights)[!used]
  if(length(remaining) > 0){
    out[remaining] <- out[remaining] + weights[remaining]
  }

  out[abs(out) > .prior_linear_density_zero_tol()]
}

.prior_density_context_density <- function(context, weights,
                                           source_transforms = NULL,
                                           output_transformation = NULL,
                                           output_transformation_arguments = NULL){

  standardized_weights <- .prior_density_context_standardized_weights(context, weights)

  if(!is.null(source_transforms)){
    source_transforms <- source_transforms[names(standardized_weights)]
  }

  .prior_linear_combination_density(
    prior_list                       = context$prior_list,
    weights                          = standardized_weights,
    n_grid                           = context$n_grid,
    tail_prob                        = context$tail_prob,
    source_transforms                = source_transforms,
    output_transformation            = output_transformation,
    output_transformation_arguments  = output_transformation_arguments
  )
}

.prior_density_model_mixture_context <- function(prior_list, column_names,
                                                  n_grid = .prior_linear_density_default_grid(),
                                                  tail_prob = .prior_linear_density_tail_prob()){

  check_list(prior_list, "prior_list")
  check_char(column_names, "column_names", check_length = FALSE)

  prior_weights <- do.call(cbind, lapply(prior_list, function(parameter_priors){
    if(is.prior(parameter_priors)){
      return(parameter_priors[["prior_weights"]])
    }
    sapply(parameter_priors, function(prior) prior[["prior_weights"]])
  }))

  if(!all(prior_weights[, 1] == prior_weights)){
    stop("The model prior distributions are not aligned across parameters.", call. = FALSE)
  }

  model_weights <- prior_weights[, 1]
  model_weights <- model_weights / sum(model_weights)

  out <- list(
    prior_list   = prior_list,
    column_names = column_names,
    model_weights = model_weights,
    n_grid       = n_grid,
    tail_prob    = tail_prob
  )
  class(out) <- "prior_density_model_mixture_context"
  return(out)
}

.prior_density_condition_component <- function(prior){

  if(is.prior.spike_and_slab(prior)){
    components <- attr(prior, "components")
    if(!all(components %in% c("null", "alternative"))){
      stop("conditional mixture posterior distributions are available only for 'null' and 'alternative' components", call. = FALSE)
    }

    inclusion <- mean(.get_spike_and_slab_inclusion(prior))
    inclusion <- min(max(inclusion, 0), 1)
    probabilities <- ifelse(components == "alternative", inclusion, 1 - inclusion)

    return(lapply(seq_along(prior), function(i){
      list(
        prior       = prior[[i]],
        probability = probabilities[i],
        alternative = components[i] == "alternative"
      )
    }))
  }

  if(is.prior.mixture(prior)){
    components <- attr(prior, "components")
    if(!all(components %in% c("null", "alternative"))){
      stop("conditional mixture posterior distributions are available only for 'null' and 'alternative' components", call. = FALSE)
    }

    prior_weights <- attr(prior, "prior_weights")
    prior_weights <- prior_weights / sum(prior_weights)

    return(lapply(seq_along(prior), function(i){
      list(
        prior       = prior[[i]],
        probability = prior_weights[i],
        alternative = components[i] == "alternative"
      )
    }))
  }

  list(list(
    prior       = prior,
    probability = 1,
    alternative = TRUE
  ))
}

.prior_density_copy_parent_attributes <- function(component, parent){

  parent_attributes <- attributes(parent)
  skip <- c("class", "names", "components", "prior_weights", "inclusion_prior")

  for(attribute in setdiff(names(parent_attributes), skip)){
    if(is.null(attr(component, attribute, exact = TRUE))){
      attr(component, attribute) <- parent_attributes[[attribute]]
    }
  }

  component
}

.prior_density_condition_models <- function(prior_list, conditional, conditional_rule){

  conditional <- unique(conditional[conditional %in% names(prior_list)])
  if(length(conditional) == 0){
    return(NULL)
  }

  options <- lapply(conditional, function(parameter){
    .prior_density_condition_component(prior_list[[parameter]])
  })
  names(options) <- conditional

  option_grid <- expand.grid(lapply(options, seq_along))
  keep <- logical(nrow(option_grid))
  model_weights <- numeric(nrow(option_grid))

  for(i in seq_len(nrow(option_grid))){
    alternatives <- logical(length(conditional))
    probabilities <- numeric(length(conditional))

    for(j in seq_along(conditional)){
      option <- options[[j]][[option_grid[i, j]]]
      alternatives[j] <- option$alternative
      probabilities[j] <- option$probability
    }

    keep[i] <- if(conditional_rule == "AND") all(alternatives) else any(alternatives)
    model_weights[i] <- prod(probabilities)
  }

  option_grid <- option_grid[keep & model_weights > 0, , drop = FALSE]
  model_weights <- model_weights[keep & model_weights > 0]

  if(nrow(option_grid) == 0){
    return(list(prior_lists = list(), weights = numeric()))
  }

  prior_lists <- lapply(seq_len(nrow(option_grid)), function(i){
    model_prior_list <- prior_list
    for(j in seq_along(conditional)){
      parameter <- conditional[j]
      component <- options[[j]][[option_grid[i, j]]]$prior
      model_prior_list[[parameter]] <- .prior_density_copy_parent_attributes(
        component = component,
        parent    = prior_list[[parameter]]
      )
    }
    model_prior_list
  })

  list(
    prior_lists = prior_lists,
    weights     = model_weights / sum(model_weights)
  )
}

.prior_density_conditional_context <- function(prior_list, column_names, conditional,
                                               conditional_rule = "AND", formula_scale = NULL,
                                               n_grid = .prior_linear_density_default_grid(),
                                               tail_prob = .prior_linear_density_tail_prob()){

  condition_models <- .prior_density_condition_models(prior_list, conditional, conditional_rule)
  if(is.null(condition_models)){
    return(.prior_density_context(prior_list, column_names, formula_scale, n_grid, tail_prob))
  }

  out <- list(
    prior_list      = prior_list,
    column_names    = column_names,
    formula_scale   = formula_scale,
    prior_lists     = condition_models$prior_lists,
    model_weights   = condition_models$weights,
    n_grid          = n_grid,
    tail_prob       = tail_prob
  )
  class(out) <- "prior_density_conditional_context"
  return(out)
}

.prior_density_model_mixture_density <- function(context, weights,
                                                  output_transformation = NULL,
                                                  output_transformation_arguments = NULL){

  if(!inherits(context, "prior_density_model_mixture_context")){
    stop("'context' must be a prior density model-mixture context.", call. = FALSE)
  }

  dists <- vector("list", length(context$model_weights))
  for(model_i in seq_along(context$model_weights)){
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

    dists[[model_i]] <- .prior_linear_combination_density(
      prior_list = model_prior_list,
      weights    = weights,
      n_grid     = context$n_grid,
      tail_prob  = context$tail_prob
    )
  }

  dx <- min(vapply(dists, function(dist){
    if(!is.null(dist$density) && length(dist$density$x) > 1){
      return(dist$density$x[2] - dist$density$x[1])
    }
    Inf
  }, numeric(1)))
  if(!is.finite(dx)){
    dx <- NA_real_
  }

  dist <- .prior_linear_density_mix(
    dists   = dists,
    weights = context$model_weights,
    dx      = dx,
    n_grid  = context$n_grid
  )

  .prior_linear_density_transform(dist, output_transformation,
                                  output_transformation_arguments,
                                  n_grid = context$n_grid)
}

.prior_density_conditional_context_density <- function(context, weights,
                                                       source_transforms = NULL,
                                                       output_transformation = NULL,
                                                       output_transformation_arguments = NULL){

  if(!inherits(context, "prior_density_conditional_context")){
    stop("'context' must be a conditional prior density context.", call. = FALSE)
  }

  if(length(context$prior_lists) == 0){
    return(.prior_linear_density_point(0))
  }

  dists <- lapply(context$prior_lists, function(prior_list){
    if(!is.null(context$formula_scale) && length(context$formula_scale) > 0){
      component_context <- .prior_density_context(
        prior_list    = prior_list,
        column_names  = context$column_names,
        formula_scale = context$formula_scale,
        n_grid        = context$n_grid,
        tail_prob     = context$tail_prob
      )
      return(.prior_density_context_density(
        context           = component_context,
        weights           = weights,
        source_transforms = source_transforms
      ))
    }

    .prior_linear_combination_density(
      prior_list        = prior_list,
      weights           = weights,
      n_grid            = context$n_grid,
      tail_prob         = context$tail_prob,
      source_transforms = source_transforms
    )
  })

  dx <- min(vapply(dists, function(dist){
    if(!is.null(dist$density) && length(dist$density$x) > 1){
      return(dist$density$x[2] - dist$density$x[1])
    }
    Inf
  }, numeric(1)))
  if(!is.finite(dx)){
    dx <- NA_real_
  }

  dist <- .prior_linear_density_mix(
    dists   = dists,
    weights = context$model_weights,
    dx      = dx,
    n_grid  = context$n_grid
  )

  .prior_linear_density_transform(dist, output_transformation,
                                  output_transformation_arguments,
                                  n_grid = context$n_grid)
}

.prior_density_build_context <- function(prior_list, column_names, formula_scale = NULL,
                                         n_grid = .prior_linear_density_default_grid(),
                                         tail_prob = .prior_linear_density_tail_prob(),
                                         conditional = NULL,
                                         conditional_rule = "AND"){

  if(length(conditional) > 0){
    return(.prior_density_conditional_context(
      prior_list       = prior_list,
      column_names     = column_names,
      conditional      = conditional,
      conditional_rule = conditional_rule,
      formula_scale    = formula_scale,
      n_grid           = n_grid,
      tail_prob        = tail_prob
    ))
  }

  if(all(vapply(prior_list, is.prior, logical(1)))){
    return(.prior_density_context(prior_list, column_names, formula_scale, n_grid, tail_prob))
  }

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    stop("Formula-scale prior densities for model-list mixtures are not implemented.", call. = FALSE)
  }

  .prior_density_model_mixture_context(prior_list, column_names, n_grid, tail_prob)
}

.prior_density_from_context <- function(context, weights,
                                        source_transforms = NULL,
                                        output_transformation = NULL,
                                        output_transformation_arguments = NULL){

  if(inherits(context, "prior_density_context")){
    return(.prior_density_context_density(
      context                         = context,
      weights                         = weights,
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    if(!is.null(source_transforms)){
      stop("Source transformations are not supported for model-list prior mixtures.", call. = FALSE)
    }
    return(.prior_density_model_mixture_density(
      context                         = context,
      weights                         = weights,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  if(inherits(context, "prior_density_conditional_context")){
    return(.prior_density_conditional_context_density(
      context                         = context,
      weights                         = weights,
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  stop("Unknown prior density context.", call. = FALSE)
}

.prior_density_coefficient_weights <- function(column_names, parameter){

  weights <- rep(0, length(column_names))
  names(weights) <- column_names
  if(parameter %in% names(weights)){
    weights[[parameter]] <- 1
  }
  weights
}

.generate_transformed_prior_densities <- function(prior_list, column_names, formula_scale = NULL,
                                                  conditional = NULL, conditional_rule = "AND",
                                                  n_grid = .prior_linear_density_default_grid(),
                                                  tail_prob = .prior_linear_density_tail_prob()){

  context <- .prior_density_build_context(
    prior_list       = prior_list,
    column_names     = column_names,
    formula_scale    = formula_scale,
    n_grid           = n_grid,
    tail_prob        = tail_prob,
    conditional      = conditional,
    conditional_rule = conditional_rule
  )

  prior_columns <- unlist(lapply(names(prior_list), function(parameter){
    .prior_linear_prior_columns(parameter, prior_list[[parameter]])
  }), use.names = FALSE)
  prior_columns <- intersect(prior_columns, column_names)

  out <- list()
  for(parameter in prior_columns){
    weights <- .prior_density_coefficient_weights(column_names, parameter)
    source_transforms <- NULL
    output_transformation <- NULL

    for(transform in context$transforms){
      if(transform$log_intercept && identical(parameter, transform$intercept)){
        weights <- rep(0, length(column_names))
        names(weights) <- column_names
        weights[transform$columns] <- transform$matrix[parameter, ]

        source_transforms <- rep(NA_character_, length(column_names))
        names(source_transforms) <- column_names
        source_transforms[[transform$intercept]] <- "log"
        output_transformation <- "exp"
        break
      }
    }

    if(!is.null(output_transformation)){
      out[[parameter]] <- .prior_linear_combination_density(
        prior_list             = prior_list,
        weights                = weights,
        n_grid                 = n_grid,
        tail_prob              = tail_prob,
        source_transforms      = source_transforms,
        output_transformation  = output_transformation
      )
    }else{
      out[[parameter]] <- .prior_density_from_context(
        context               = context,
        weights               = weights,
        source_transforms     = source_transforms,
        output_transformation = output_transformation
      )
    }
  }

  attr(out, "context") <- context
  class(out) <- c("prior_density_list", "list")
  return(out)
}

.prior_factor_level_weight_matrix <- function(sample_metadata, parameter, samples = NULL){

  metadata <- sample_metadata
  if(!is.null(samples)){
    metadata <- .add_factor_metadata_from_named_objects(metadata, parameter, samples)
  }

  if(isTRUE(attr(metadata, "orthonormal")) ||
     isTRUE(attr(metadata, "meandif")) ||
     inherits(metadata, "mixed_posteriors.treatment_transformed")){
    design_info <- .factor_term_design_from_metadata(metadata)
    weights <- design_info$design
    rownames(weights) <- .factor_contrast_parameter_names(
      parameter   = parameter,
      level_names = design_info$level_names,
      cell_names  = design_info$cell_names
    )
    prior <- attr(metadata, "prior_list")
    if(!is.null(prior) && !is.prior(prior)){
      prior_candidates <- prior[vapply(prior, is.prior, logical(1))]
      prior <- if(length(prior_candidates) > 0) prior_candidates[[1]] else NULL
    }
    prior_columns <- if(!is.null(prior) && is.prior(prior)){
      .prior_linear_prior_columns(parameter, prior)
    }else{
      NULL
    }
    colnames(weights) <- if(length(prior_columns) == ncol(weights)){
      prior_columns
    }else{
      colnames(as.matrix(sample_metadata))
    }
    return(weights)
  }

  raw_columns <- colnames(as.matrix(sample_metadata))
  weights <- diag(length(raw_columns))
  rownames(weights) <- raw_columns
  prior <- attr(metadata, "prior_list")
  if(!is.null(prior) && !is.prior(prior)){
    prior_candidates <- prior[vapply(prior, is.prior, logical(1))]
    prior <- if(length(prior_candidates) > 0) prior_candidates[[1]] else NULL
  }
  prior_columns <- if(!is.null(prior) && is.prior(prior)){
    .prior_linear_prior_columns(parameter, prior)
  }else{
    NULL
  }
  colnames(weights) <- if(length(prior_columns) == ncol(weights)){
    prior_columns
  }else{
    raw_columns
  }
  return(weights)
}
