.prior_linear_density_default_grid <- function(){
  4096L
}

.prior_linear_density_tail_prob <- function(){
  1e-4
}

.prior_linear_density_product_grid <- function(){
  1024L
}

.prior_linear_density_grid_tol <- function(){
  sqrt(.Machine$double.eps)
}
.prior_linear_source_transform <- function(source_transform){

  if(is.null(source_transform) || length(source_transform) == 0 || is.na(source_transform)){
    return(NULL)
  }

  unname(source_transform)
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

.prior_linear_density_normalize <- function(dist, warn = FALSE){

  point_mass <- if(!is.null(dist$points) && nrow(dist$points) > 0) sum(dist$points$p) else 0
  density_mass <- if(!is.null(dist$density)) dist$density$mass else 0
  total_mass <- point_mass + density_mass

  if(!is.finite(total_mass) || total_mass <= 0){
    stop("The computed prior density has zero total mass.", call. = FALSE)
  }

  # Intermediate mixture / product pieces may intentionally carry partial mass.
  # Callers of finished densities should pass warn = TRUE.
  if(isTRUE(warn) && abs(total_mass - 1) > 1e-6){
    warning(
      "Computed prior density mass was ",
      format(total_mass, digits = 8),
      "; renormalizing to 1.",
      call. = FALSE
    )
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

  out <- Re(out[seq_len(n)])
  error_bound <- 64 * .Machine$double.eps * max(1, log2(n_fft)) *
    max(1, max(abs(out)))
  minimum <- min(out)
  if(minimum < -error_bound){
    stop(
      "FFT convolution produced a negative density value beyond its ",
      "floating-point error bound.",
      call. = FALSE
    )
  }
  negative <- out < 0
  diagnostics <- list(
    error_bound = error_bound,
    minimum_unclipped_value = minimum,
    clipped_negative_sum = sum(-out[negative]),
    clipped_value_count = sum(negative)
  )
  out[negative] <- 0
  attr(out, "fft_clipping") <- diagnostics
  out
}

.prior_linear_density_aggregate_points <- function(points, dx){

  if(is.null(points) || nrow(points) == 0){
    return(.prior_linear_density_empty_points())
  }

  points <- points[is.finite(points$x) & is.finite(points$p) & points$p > 0, , drop = FALSE]
  if(nrow(points) == 0){
    return(.prior_linear_density_empty_points())
  }

  key <- sprintf("%a", points$x)

  split_points <- split(points, key)
  out <- do.call(rbind, lapply(split_points, function(p){
    data.frame(
      x = p$x[1L],
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
  grid_normalization <- NULL
  if(length(densities) > 0){
    x_min <- min(vapply(densities, function(d) min(d$x), numeric(1)))
    x_max <- max(vapply(densities, function(d) max(d$x), numeric(1)))

    if(!is.finite(dx) || dx <= 0){
      if(is.null(n_grid)){
        n_grid <- .prior_linear_density_default_grid()
      }
      dx <- (x_max - x_min) / max(1, n_grid - 1)
    }

    if(!is.finite(dx) || dx <= 0 || x_min == x_max){
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
      grid_normalization <- list(
        captured_numerical_mass = area,
        target_continuous_mass = density_mass,
        normalization_factor = density_mass / area
      )
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
  if(!is.null(grid_normalization)){
    attr(out, "grid_normalization") <- grid_normalization
  }
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
    y <- .prior_linear_density_fft_convolve(lhs$density$y, rhs$density$y)
    fft_clipping <- attr(y, "fft_clipping", exact = TRUE)
    y <- as.numeric(y) * dx
    if(!is.null(fft_clipping)){
      fft_clipping$clipped_negative_mass <- fft_clipping$clipped_negative_sum * dx
    }
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

  out <- .prior_linear_density_coalesce(
    densities = densities,
    points = points,
    dx = dx,
    n_grid = max(lhs$n_grid, rhs$n_grid)
  )
  inherited_clipping <- c(
    attr(lhs, "fft_clipping", exact = TRUE),
    attr(rhs, "fft_clipping", exact = TRUE)
  )
  if(exists("fft_clipping", inherits = FALSE) && !is.null(fft_clipping)){
    inherited_clipping <- c(inherited_clipping, list(fft_clipping))
  }
  if(length(inherited_clipping) > 0L){
    attr(out, "fft_clipping") <- inherited_clipping
  }
  out
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
    values <- c(values, dist$points$x[dist$points$p > 0])
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

  if(scale == 0){
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
    lhs_continuous <- lhs
    lhs_continuous$points <- .prior_linear_density_empty_points()
    for(i in seq_len(nrow(rhs$points))){
      scaled <- .prior_linear_density_scaled(lhs_continuous, rhs$points$x[i],
                                             mass = lhs$density$mass * rhs$points$p[i],
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
    rhs_continuous <- rhs
    rhs_continuous$points <- .prior_linear_density_empty_points()
    for(i in seq_len(nrow(lhs$points))){
      scaled <- .prior_linear_density_scaled(rhs_continuous, lhs$points$x[i],
                                             mass = rhs$density$mass * lhs$points$p[i],
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
    if(z_range[1] == z_range[2]){
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
      zero_tol <- .prior_linear_density_grid_tol() * max(1, max(abs(x), na.rm = TRUE))
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

      if(any(!is.finite(y))){
        stop(
          "Product-density quadrature produced non-finite values.",
          call. = FALSE
        )
      }
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
