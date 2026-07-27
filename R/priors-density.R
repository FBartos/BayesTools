#' @title Prior density
#'
#' @description Computes density of a prior
#' distribution across a range of values.
#'
#' @param x a prior
#' @param x_seq sequence of x coordinates
#' @param x_range vector of length two with
#' lower and upper range for the support
#' (used if \code{x_seq} is unspecified)
#' @param x_range_quant quantile used for
#' automatically obtaining \code{x_range}
#' if both \code{x_range} and \code{x_seq}
#' are unspecified. Defaults to \code{0.005}
#' for all but Cauchy, Student-t, Gamma, inverse-gamma,
#' moment, and inverse-moment distributions that use
#' \code{0.010}.
#' @param n_points number of equally spaced points
#' in the \code{x_range} if \code{x_seq} is unspecified
#' @param n_samples number of samples from the prior
#' distribution if the density cannot be obtained
#' analytically (or if samples are forced with
#' \code{force_samples = TRUE})
#' @param force_samples should prior be sampled instead
#' of obtaining analytic solution whenever possible
#' @param individual should individual densities be returned
#' (e.g., in case of weightfunction)
#' @param transformation transformation to be applied
#' to the prior distribution. Either a character
#' specifying one of the prepared transformations:
#' \describe{
#'   \item{lin}{linear transformation in form of \code{a + b*x}}
#'   \item{tanh}{hyperbolic tangent transformation}
#'   \item{exp}{exponential transformation}
#' }, or a list containing the transformation function \code{fun},
#' inverse transformation function \code{inv}, and derivative of the
#' transformation \code{jac}, evaluated on the original support. See examples
#' for details.
#' @param transformation_arguments a list with named arguments for
#' the \code{transformation}
#' @param transformation_settings boolean indicating whether the
#' settings the \code{x_seq} or \code{x_range} was specified on
#' the transformed support
#' @param truncate_end whether the density should be set to zero
#' for endpoints of truncated distributions
#' @param ... additional arguments
#'
#' @return \code{density.prior} returns an object of class 'density'. For
#' Dirichlet simplex priors it returns a named list of beta-marginal density
#' objects, one for each simplex coordinate. Ordered priors return one component
#' per factor level. When a level has both discrete and continuous probability,
#' its \code{atoms} table stores exact \code{location} and \code{mass} values,
#' while its \code{continuous} table stores a density already weighted to
#' integrate to the remaining continuous mass.
#'
#' @details Sample-based density estimates for continuous priors with finite
#' support use boundary-reflected kernel density estimates. The plotting range
#' controls the evaluation grid, while reflection is based on the prior's true
#' truncation bounds.
#'
#' @importFrom stats density
#' @seealso [prior()]
#' @rdname density.prior
#' @export
density.prior <- function(x,
                          x_seq = NULL, x_range = NULL, x_range_quant = NULL, n_points = 1000,
                          n_samples = 10000, force_samples = FALSE, individual = FALSE,
                          transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE, truncate_end = TRUE, ...){

  # input check
  .check_prior(x, "x")
  check_real(x_seq, "x_seq", check_length = 0, allow_NULL = TRUE)
  check_real(x_range, "x_range", check_length = 2, allow_NULL = TRUE)
  if(!is.null(x_range) && x_range[1] > x_range[2])
    stop("The lower range limit must be lower than the upper range limit.")
  check_real(x_range_quant, "x_range_quant", lower = 0, upper = 1, allow_NULL = TRUE)
  check_int(n_points, "n_points",  lower = 2)
  check_int(n_samples, "n_samples", lower = 1)
  check_bool(force_samples, "force_samples")
  check_bool(individual, "individual")
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)
  check_bool(truncate_end, "truncate_end")


  ### setting the range
  # get plotting range if not specified
  if(is.null(x_range)){
    if(!is.null(x_seq)){
      x_range <- range(x_seq)
    }else{
      if(is_prior_phacking(x) || is_prior_bias(x)){
        .selection_prior_stop_unsupported_generic("density", x)
      }else if(!individual & (is.prior.PET(x) | is.prior.PEESE(x))){
        x_range <- c(0, 1)
      }else if(!individual & is.prior.weightfunction(x)){
        x_range <- c(0, 1)
      }else if(is.prior.spike_and_slab(x)){
        x_range <- range(c(range(.get_spike_and_slab_variable(x), if(is.null(x_range_quant)) .range.prior_quantile_default(.get_spike_and_slab_variable(x)) else x_range_quant), 0))
      }else if(is.prior.ordered(x)){
        x_range <- .prior_ordered_range(x, quantiles = x_range_quant)
      }else if(is.prior.discrete(x)){
        x_range <- c(x[["truncation"]][["lower"]], x[["truncation"]][["upper"]])
      }else{
        x_range <- range(x, if(is.null(x_range_quant)) .range.prior_quantile_default(x) else x_range_quant)
      }
    }
  }

  # get the x_seq for plotting
  if(is.null(x_seq)){
    if(is.prior.discrete(x)){
      x_seq <- seq(x_range[1], x_range[2], by = 1)
    }else{
      x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
    }
  }

  # specify it on the transformed range if requested
  if(transformation_settings & !is.null(transformation)){
    x_seq   <- .density.prior_transformation_inv_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_inv_x(x_range, transformation, transformation_arguments)
  }


  # use the corresponding density subfunction
  if(is.prior.weightfunction(x)){
    out <- .density.prior.weightfunction(x, x_seq, x_range, n_points, n_samples, force_samples, individual)
  }else if(is_prior_phacking(x) || is_prior_bias(x)){
    .selection_prior_stop_unsupported_generic("density", x)
  }else if(is.prior.PET(x) | is.prior.PEESE(x)){
    out <- .density.prior.PETPEESE(x, x_seq, x_range, n_points, n_samples, force_samples, individual, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.spike_and_slab(x)){
    out <- .density.prior.spike_and_slab(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.point(x)){
    out <- .density.prior.point(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments)
  }else if(is.prior.ordered(x)){
    out <- .density.prior.ordered(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.orthonormal(x) | is.prior.meandif(x)){
    out <- .density.prior.orthonormal_or_meandif(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.simplex(x)){
    out <- .density.prior.simplex(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.simple(x)){
    out <- .density.prior.simple(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }

  if(!is.null(transformation)){
    attr(out, "transformation") <- transformation
  }

  return(out)
}

.density.prior.simple                 <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  boundary_reflection <- FALSE

  # get the samples to estimate density / obtain the density directly
  if(force_samples | .density.prior_need_samples(x)){
    x_sam <- rng(x, n_samples)
    if(is.prior.discrete(x)){
      x_seq <- unique(round(x_seq))
      x_den <- vapply(x_seq, function(x_i) mean(x_sam == x_i), numeric(1))
    }else{
      x_density <- .density_kde_boundary(
        x      = x_sam,
        n      = n_points,
        from   = x_range[1],
        to     = x_range[2],
        bounds = c(x$truncation[["lower"]], x$truncation[["upper"]])
      )
      x_seq               <- x_density$x
      x_den               <- x_density$y
      boundary_reflection <- isTRUE(attr(x_density, "boundary_reflection"))
    }
  }else{

    if(is.prior.discrete(x)){
      x_seq <- unique(round(x_seq))
    }

    x_den <- mpdf(x, x_seq)
    x_sam <- NULL
  }


  # set the endpoints to zero if they correspond to truncation
  if(truncate_end){
    if(isTRUE(all.equal(x$truncation[["lower"]], x_seq[1])) | x$truncation[["lower"]] >= x_seq[1]){
      x_den <- c(0, x_den)
      x_seq <- c(x_seq[1], x_seq)
    }
    if(isTRUE(all.equal(x$truncation[["upper"]], x_seq[length(x_seq)])) | x$truncation[["upper"]] <= x_seq[length(x_seq)]){
      x_den <- c(x_den, 0)
      x_seq <- c(x_seq, x_seq[length(x_seq)])
    }
  }



  # transform the output, if requested
  if(!is.null(transformation)){
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
    if(!is.null(x_sam)){
      x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
    }
    x_den   <- .density.prior_transformation_y(x_seq, x_den, transformation, transformation_arguments)
  }


  # create the output object
  out <- list(
    call    = call("density", print(x, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_den,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", "density.prior.simple")
  attr(out, "x_range") <- x_range
  attr(out, "y_range") <- c(0, max(x_den))
  if(boundary_reflection){
    attr(out, "boundary_reflection") <- TRUE
  }

  return(out)
}

.density.prior.ordered                <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  mixed <- .density.prior.ordered_mixed(
    x = x,
    x_seq = x_seq,
    x_range = x_range,
    n_points = n_points,
    n_samples = n_samples,
    force_samples = force_samples,
    transformation = transformation,
    transformation_arguments = transformation_arguments
  )
  if(!is.null(mixed)){
    return(mixed)
  }

  if(!force_samples && is.null(transformation)){
    direct <- .density.prior.ordered_direct(x, x_seq, n_points, transformation, transformation_arguments, truncate_end)
    if(!is.null(direct)){
      return(direct)
    }
  }

  samples <- rng(x, n_samples, transform_factor_samples = TRUE)
  if(!is.matrix(samples)){
    samples <- matrix(samples, ncol = 1L)
  }

  out <- vector("list", ncol(samples))
  names(out) <- colnames(samples)

  for(i in seq_len(ncol(samples))){
    component_samples <- samples[, i]

    if(all(component_samples == component_samples[1])){
      out[[i]] <- .density.prior.point(
        prior("point", list(location = component_samples[1])),
        x_seq,
        x_range,
        n_points,
        n_samples,
        force_samples = TRUE,
        transformation,
        transformation_arguments
      )
    }else{
      x_density <- .density_kde_boundary(
        x    = component_samples,
        n    = n_points,
        from = x_range[1],
        to   = x_range[2]
      )
      x_values <- x_density$x
      y_values <- x_density$y

      if(!is.null(transformation)){
        x_values         <- .density.prior_transformation_x(x_values, transformation, transformation_arguments)
        component_samples <- .density.prior_transformation_x(component_samples, transformation, transformation_arguments)
        y_values         <- .density.prior_transformation_y(x_values, y_values, transformation, transformation_arguments)
      }

      out[[i]] <- list(
        call    = call("density", print(x, silent = TRUE)),
        bw      = x_density$bw,
        n       = n_points,
        x       = x_values,
        y       = y_values,
        samples = component_samples
      )
      class(out[[i]]) <- c("density", "density.prior", "density.prior.simple")
      attr(out[[i]], "x_range") <- range(x_values)
      attr(out[[i]], "y_range") <- range(y_values)
    }

    if(inherits(out[[i]], "density.prior.point")){
      attr(out[[i]], "x_range") <- range(c(attr(out[[i]], "x_range"), out[[i]]$x), na.rm = TRUE)
    }
    class(out[[i]]) <- c("density.prior.ordered_component", class(out[[i]]))
    attr(out[[i]], "component") <- i
    attr(out[[i]], "component_name") <- names(out)[i]
  }

  attr(out, "x_range")        <- range(unlist(lapply(out, attr, which = "x_range")), na.rm = TRUE)
  attr(out, "y_range")        <- range(unlist(lapply(out, attr, which = "y_range")), na.rm = TRUE)
  attr(out, "parameter_name") <- names(out)
  class(out) <- c("density.prior.ordered", "list")

  out
}

.density.prior.ordered_mixed <- function(x, x_seq, x_range, n_points,
                                         n_samples, force_samples,
                                         transformation,
                                         transformation_arguments){

  x <- .prior_ordered_default_bound(x)
  target_dx <- diff(range(x_seq)) / max(1, length(x_seq) - 1L)
  if(!is.finite(target_dx) || target_dx <= 0){
    target_dx <- NA_real_
  }
  total <- .prior_ordered_total_linear_distribution(
    total = x$total,
    dx = target_dx,
    n_grid = n_points,
    tail_prob = .prior_linear_density_tail_prob()
  )
  if(is.null(total) || is.null(total$density) ||
     is.null(total$points) || nrow(total$points) == 0L){
    return(NULL)
  }

  metadata <- .prior_ordered_metadata(x)
  if(length(metadata$ordered_terms) != 1L ||
     metadata$theta_dim != 1L ||
     length(metadata$allocations) != 1L){
    stop(
      "Mixed-measure ordered densities currently require one ordered term ",
      "and one scalar total. Split the interaction into explicitly named ",
      "terms before requesting its density.",
      call. = FALSE
    )
  }

  design_info <- .factor_term_design_from_metadata(x)
  level_names <- design_info[["cell_names"]]
  component_names <- .factor_contrast_parameter_names(
    parameter = metadata$parameter_name,
    level_names = .factor_level_list(x),
    cell_names = level_names
  )
  weights <- design_info$design
  colnames(weights) <- .JAGS_prior_factor_names(metadata$parameter_name, x)

  samples <- NULL
  if(force_samples){
    samples <- rng(x, n_samples, transform_factor_samples = TRUE)
    if(!is.matrix(samples)){
      samples <- matrix(samples, ncol = 1L)
    }
  }

  densities <- vector("list", nrow(weights))
  names(densities) <- component_names
  for(i in seq_len(nrow(weights))){
    dist <- .prior_ordered_linear_distribution(
      ordered_prior = x,
      weights = weights[i, ],
      indices = seq_len(ncol(weights)),
      dx = target_dx,
      n_grid = n_points,
      tail_prob = .prior_linear_density_tail_prob()
    )
    dist <- .density.prior.ordered_regrid_mixed(dist, x_seq)
    if(!is.null(transformation)){
      dist <- .prior_linear_density_transform(
        dist,
        transformation,
        transformation_arguments,
        n_grid = n_points
      )
    }
    component_samples <- if(is.null(samples)){
      NULL
    }else{
      values <- samples[, i]
      if(!is.null(transformation)){
        values <- .density.prior_transformation_x(
          values,
          transformation,
          transformation_arguments
        )
      }
      values
    }
    densities[[i]] <- .density.prior.ordered_mixed_component(
      dist = dist,
      prior = x,
      n_points = n_points,
      samples = component_samples,
      component = i,
      component_name = component_names[[i]],
      x_range = x_range,
      transformation = transformation,
      transformation_arguments = transformation_arguments
    )
  }

  attr(densities, "x_range") <- range(unlist(lapply(
    densities,
    attr,
    which = "x_range"
  )), na.rm = TRUE)
  attr(densities, "y_range") <- range(unlist(lapply(
    densities,
    attr,
    which = "y_range"
  )), na.rm = TRUE)
  attr(densities, "parameter_name") <- names(densities)
  attr(densities, "method") <- "analytic_mixed_measure"
  attr(densities, "measure_schema_version") <- 1L
  class(densities) <- c("density.prior.ordered", "list")
  densities
}

.density.prior.ordered_regrid_mixed <- function(dist, x_seq){

  if(is.null(dist$density) || dist$density$mass <= 0){
    return(dist)
  }
  y <- stats::approx(
    dist$density$x,
    dist$density$y,
    xout = x_seq,
    yleft = 0,
    yright = 0
  )$y
  area <- .density.prior.ordered_curve_integral(x_seq, y)
  if(!is.finite(area) || area <= 0){
    stop(
      "The continuous part of an ordered mixed-measure prior has no ",
      "numerical mass on the requested density grid.",
      call. = FALSE
    )
  }
  dist$density$x <- x_seq
  dist$density$y <- y / area
  dist$n_grid <- length(x_seq)
  attr(dist, "ordered_grid_diagnostics") <- list(
    captured_continuous_shape_integral = area
  )
  dist
}

.density.prior.ordered_curve_integral <- function(x, y){

  if(length(x) < 2L || length(y) != length(x)){
    return(0)
  }
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

.density.prior.ordered_mixed_component <- function(
    dist, prior, n_points, samples, component, component_name,
    x_range, transformation, transformation_arguments){

  atoms <- data.frame(location = numeric(), mass = numeric())
  if(!is.null(dist$points) && nrow(dist$points) > 0L){
    atoms <- data.frame(
      location = dist$points$x,
      mass = dist$points$p
    )
  }
  continuous <- NULL
  x_values <- numeric()
  y_values <- numeric()
  if(!is.null(dist$density) && dist$density$mass > 0){
    y_values <- dist$density$y * dist$density$mass
    continuous <- data.frame(
      x = dist$density$x,
      density = y_values
    )
    attr(continuous, "mass") <- dist$density$mass
    x_values <- continuous$x
  }
  if(length(x_values) == 0L){
    x_values <- atoms$location
    y_values <- atoms$mass
  }

  atom_mass <- sum(atoms$mass)
  continuous_mass <- if(is.null(continuous)) 0 else attr(continuous, "mass")
  continuous_integral <- if(is.null(continuous)){
    0
  }else{
    .density.prior.ordered_curve_integral(
      continuous$x,
      continuous$density
    )
  }
  mass_bound <- 128 * .Machine$double.eps
  if(abs(atom_mass + continuous_mass - 1) > mass_bound ||
     abs(continuous_integral - continuous_mass) > 1e-10){
    stop(
      "The ordered mixed-measure density did not preserve unit probability ",
      "mass.",
      call. = FALSE
    )
  }

  transformed_range <- x_range
  if(!is.null(transformation)){
    transformed_range <- .density.prior_transformation_x(
      transformed_range,
      transformation,
      transformation_arguments
    )
  }
  component_range <- range(c(
    transformed_range,
    atoms$location,
    if(is.null(continuous)) numeric() else continuous$x
  ), na.rm = TRUE)
  y_max <- max(c(0, atoms$mass, y_values), na.rm = TRUE)
  out <- list(
    call = call("density", print(prior, silent = TRUE)),
    bw = NULL,
    n = n_points,
    x = x_values,
    y = y_values,
    samples = samples,
    atoms = atoms,
    continuous = continuous,
    diagnostics = list(
      method = "analytic_components",
      atom_mass = atom_mass,
      continuous_mass = continuous_mass,
      continuous_integral = continuous_integral,
      grid = attr(dist, "ordered_grid_diagnostics", exact = TRUE),
      numerical = attr(dist, "numerical_diagnostics", exact = TRUE),
      ordered_measure = attr(dist, "ordered_measure", exact = TRUE)
    ),
    transformation = list(
      name = transformation,
      arguments = transformation_arguments
    )
  )
  component_class <- if(is.null(continuous)){
    "density.prior.point"
  }else{
    "density.prior.simple"
  }
  class(out) <- c(
    "density.prior.ordered_component",
    "density.prior.mixed_measure",
    "density",
    "density.prior",
    component_class
  )
  attr(out, "x_range") <- component_range
  attr(out, "y_range") <- c(0, y_max)
  attr(out, "component") <- component
  attr(out, "component_name") <- component_name
  attr(out, "measure_schema_version") <- 1L
  out
}

.density.prior.ordered_direct         <- function(x, x_seq, n_points, transformation, transformation_arguments, truncate_end){

  x <- .prior_ordered_default_bound(x)
  metadata <- .prior_ordered_metadata(x)

  if(length(metadata$ordered_terms) != 1L || metadata$theta_dim != 1L ||
     length(metadata$allocations) != 1L || !is.prior.simple(x$total) ||
     is.prior.discrete(x$total) || is.prior.point(x$total)){
    return(NULL)
  }

  record <- metadata$allocations[[1]]
  if(!record$spec$type %in% c("fixed", "dirichlet")){
    return(NULL)
  }

  level_names <- .factor_term_design_from_metadata(x)[["cell_names"]]
  component_names <- .factor_contrast_parameter_names(
    parameter = metadata$parameter_name,
    level_names = .factor_level_list(x),
    cell_names = level_names
  )

  densities <- vector("list", length(component_names))
  names(densities) <- component_names

  if(identical(record$spec$type, "fixed")){
    cumulative <- if(identical(x$contrast, "cumulative")){
      c(0, cumsum(record$spec$weights))
    }else{
      cumsum(record$spec$weights)
    }
    for(i in seq_along(cumulative)){
      densities[[i]] <- .density.prior.ordered_scaled_total(
        total = x$total,
        scale = cumulative[[i]],
        x_seq = x_seq,
        n_points = n_points
      )
    }
  }else{
    alpha <- record$spec$alpha
    D <- length(alpha)
    for(i in seq_along(densities)){
      m <- if(identical(x$contrast, "cumulative")) i - 1L else i
      if(m == 0L){
        densities[[i]] <- .density.prior.point(prior("point", list(location = 0)), x_seq, range(x_seq), n_points, n_samples = 1L, force_samples = FALSE, transformation = NULL, transformation_arguments = NULL)
      }else if(m == D){
        densities[[i]] <- .density.prior.simple(x$total, x_seq, range(x_seq), n_points, n_samples = 1L, force_samples = FALSE, transformation = NULL, transformation_arguments = NULL, truncate_end = truncate_end)
      }else{
        densities[[i]] <- .density.prior.ordered_dirichlet_product(
          total = x$total,
          alpha1 = sum(alpha[seq_len(m)]),
          alpha2 = sum(alpha[(m + 1L):D]),
          x_seq = x_seq,
          n_points = n_points
        )
        if(is.null(densities[[i]])){
          return(NULL)
        }
      }
    }
  }

  for(i in seq_along(densities)){
    if(inherits(densities[[i]], "density.prior.point")){
      attr(densities[[i]], "x_range") <- range(c(attr(densities[[i]], "x_range"), densities[[i]]$x), na.rm = TRUE)
    }
    attr(densities[[i]], "component") <- i
    attr(densities[[i]], "component_name") <- names(densities)[i]
    class(densities[[i]]) <- c("density.prior.ordered_component", class(densities[[i]]))
  }

  attr(densities, "x_range")        <- range(unlist(lapply(densities, attr, which = "x_range")), na.rm = TRUE)
  attr(densities, "y_range")        <- range(unlist(lapply(densities, attr, which = "y_range")), na.rm = TRUE)
  attr(densities, "parameter_name") <- names(densities)
  attr(densities, "method")         <- "direct"
  class(densities) <- c("density.prior.ordered", "list")

  densities
}

.density.prior.ordered_scaled_total   <- function(total, scale, x_seq, n_points){

  if(scale == 0){
    return(.density.prior.point(
      prior("point", list(location = 0)),
      x_seq,
      range(x_seq),
      n_points,
      n_samples = 1L,
      force_samples = FALSE,
      transformation = NULL,
      transformation_arguments = NULL
    ))
  }

  y <- pdf(total, x_seq / scale) / abs(scale)
  out <- list(
    call    = call("density", print(total, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = y,
    samples = NULL
  )
  class(out) <- c("density", "density.prior", "density.prior.simple")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- c(0, max(y, na.rm = TRUE))
  out
}

.density.prior.ordered_dirichlet_product <- function(total, alpha1, alpha2, x_seq, n_points){

  rel_tol <- 1e-7
  abs_tol <- 1e-10
  quadrature <- vector("list", length(x_seq))
  y <- vapply(seq_along(x_seq), function(i){
    x_value <- x_seq[i]
    density_at_zero <- if(x_value == 0) pdf(total, 0) else NA_real_

    if(x_value == 0 && is.finite(density_at_zero) && density_at_zero > 0){
      if(alpha1 <= 1){
        quadrature[[i]] <<- list(
          value = Inf,
          abs.error = 0,
          message = "analytic singular boundary"
        )
        return(Inf)
      }
      value <- density_at_zero * (alpha1 + alpha2 - 1) / (alpha1 - 1)
      quadrature[[i]] <<- list(
        value = value,
        abs.error = 0,
        message = "analytic boundary"
      )
      return(value)
    }

    integration <- tryCatch(
      stats::integrate(
        function(logit_c){
          c_value <- stats::plogis(logit_c)
          out <- numeric(length(c_value))
          interior <- c_value > 0 & c_value < 1
          if(any(interior)){
            log_integrand <- lpdf(total, x_value / c_value[interior]) +
              stats::dbeta(
                c_value[interior],
                shape1 = alpha1,
                shape2 = alpha2,
                log = TRUE
              ) +
              log1p(-c_value[interior])
            out[interior] <- exp(log_integrand)
          }
          out
        },
        lower = -Inf,
        upper = Inf,
        subdivisions = 500L,
        rel.tol = rel_tol,
        abs.tol = abs_tol,
        stop.on.error = FALSE
      ),
      error = function(e) e
    )
    if(inherits(integration, "error") ||
       !identical(integration$message, "OK") ||
       !is.finite(integration$value) ||
       !is.finite(integration$abs.error) ||
       integration$abs.error > max(abs_tol, rel_tol * abs(integration$value))){
      detail <- if(inherits(integration, "error")){
        conditionMessage(integration)
      }else{
        paste0(
          integration$message,
          "; absolute error ",
          format(integration$abs.error, digits = 6)
        )
      }
      stop(
        "Ordered-prior product density quadrature failed at x = ",
        format(x_value, digits = 17), ": ", detail, ".",
        call. = FALSE
      )
    }
    quadrature[[i]] <<- list(
      value = integration$value,
      abs.error = integration$abs.error,
      message = integration$message
    )
    integration$value
  }, numeric(1))

  if(anyNA(y)){
    stop("Ordered-prior product density quadrature returned missing values.",
         call. = FALSE)
  }

  out <- list(
    call    = call("density", print(total, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = y,
    samples = NULL
  )
  class(out) <- c("density", "density.prior", "density.prior.simple")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- c(0, max(y, na.rm = TRUE))
  attr(out, "quadrature") <- quadrature
  attr(out, "quadrature_tolerance") <- c(
    relative = rel_tol,
    absolute = abs_tol
  )
  out
}

.density.prior.simplex                <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  if(!identical(x[["distribution"]], "dirichlet")){
    stop("Only Dirichlet simplex prior densities are supported.", call. = FALSE)
  }

  alpha  <- x$parameters[["alpha"]]
  alpha0 <- sum(alpha)

  out <- vector("list", length(alpha))
  names(out) <- paste0("V", seq_along(alpha))

  for(i in seq_along(alpha)){
    component_prior <- prior(
      "beta",
      list(alpha = alpha[i], beta = alpha0 - alpha[i])
    )
    out[[i]] <- .density.prior.simple(
      component_prior,
      x_seq,
      x_range,
      n_points,
      n_samples,
      force_samples,
      transformation,
      transformation_arguments,
      truncate_end
    )
    attr(out[[i]], "component")      <- i
    attr(out[[i]], "component_name") <- names(out)[i]
    class(out[[i]]) <- c("density.prior.simplex_component", class(out[[i]]))
  }

  attr(out, "x_range")        <- range(unlist(lapply(out, attr, which = "x_range")), na.rm = TRUE)
  attr(out, "y_range")        <- range(unlist(lapply(out, attr, which = "y_range")), na.rm = TRUE)
  attr(out, "parameter_name") <- names(out)
  class(out) <- c("density.prior.simplex", "list")

  return(out)
}
.density.prior.point                  <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments){

  # return the samples if requested
  if(force_samples){
    x_sam <- rng(x, n_samples)
  }else{
    x_sam <- NULL
  }

  x_seq <- x$parameters[["location"]]
  x_den <- 1


  # transform the output, if requested
  if(!is.null(transformation)){
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
    if(!is.null(x_sam)){
      x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
    }
  }


  # create the output object
  out <- list(
    call    = call("density", print(x, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_den,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", "density.prior.point", if(is.prior.orthonormal(x)) "density.prior.orthonormal" else if(is.prior.meandif(x)) "density.prior.meandif")
  attr(out, "x_range") <- x_range
  attr(out, "y_range") <- c(0, max(x_den))

  return(out)
}
.density.prior.weightfunction         <- function(x, x_seq, x_range, n_points, n_samples, force_samples, individual){

  # create either distribution for the individual weights or the whole weightfunction
  if(individual){

    out <- list()
    out_types <- .density.prior_type(x)
    components <- .weightfunction_marginal_components(x)
    component_bounds <- lapply(components, .density.prior_weightfunction_component_bounds)
    density_boundary_reflection <- rep(FALSE, length(components))

    if(force_samples | .density.prior_need_samples(x)){
      x_sam <- rng(x, n_samples)
      density_ind <- which(out_types != "point")
      densities <- lapply(density_ind, function(i){
        .density_kde_boundary(
          x      = x_sam[,i],
          n      = n_points,
          from   = x_range[1],
          to     = x_range[2],
          bounds = component_bounds[[i]]
        )
      })
      if(length(densities) > 0L){
        x_seq <- densities[[1]]$x
        x_den <- matrix(0, nrow = length(x_seq), ncol = ncol(x_sam))
        for(j in seq_along(density_ind)){
          x_den[,density_ind[j]] <- densities[[j]]$y
          density_boundary_reflection[density_ind[j]] <- isTRUE(attr(densities[[j]], "boundary_reflection"))
        }
      }else{
        x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
        x_den <- matrix(0, nrow = length(x_seq), ncol = ncol(x_sam))
      }
    }else{
      x_den <- mpdf(x, x_seq)
      x_sam <- NULL
    }

    for(i in 1:ncol(x_den)){

      temp_samples <- if(is.null(x_sam)) NULL else x_sam[,i]
      temp_x_seq <- x_seq
      temp_y_den <- x_den[,i]
      temp_bounds <- component_bounds[[i]]

      if(out_types[i] != "point"){
        at_lower_bound <- is.finite(temp_bounds[1]) &&
          (isTRUE(all.equal(temp_bounds[1], temp_x_seq[1])) | temp_bounds[1] >= temp_x_seq[1])
        at_upper_bound <- is.finite(temp_bounds[2]) &&
          (isTRUE(all.equal(temp_bounds[2], temp_x_seq[length(temp_x_seq)])) | temp_bounds[2] <= temp_x_seq[length(temp_x_seq)])

        if(at_lower_bound){
          temp_y_den <- c(0, temp_y_den)
          temp_x_seq <- c(temp_x_seq[1], temp_x_seq)
        }
        if(at_upper_bound){
          temp_y_den <- c(temp_y_den, 0)
          temp_x_seq <- c(temp_x_seq, temp_x_seq[length(temp_x_seq)])
        }
      }

      # create the output object
      if(out_types[i] == "point"){
        temp_out <- list(
          call    = call("density", print(x, silent = TRUE)),
          bw      = NULL,
          n       = n_points,
          x       = components[[i]]$location,
          y       = 1,
          samples = temp_samples
        )
      }else{
        temp_out <- list(
          call    = call("density", print(x, silent = TRUE)),
          bw      = NULL,
          n       = n_points,
          x       = temp_x_seq,
          y       = temp_y_den,
          samples = temp_samples
        )
      }


      class(temp_out) <- c("density", "density.prior", paste0("density.prior.",out_types[i]))
      attr(temp_out, "x_range") <- x_range
      attr(temp_out, "y_range") <- if(out_types[i] == "point") c(0, 1) else c(0, max(temp_y_den))
      attr(temp_out, "steps")   <- c(x$bins$lower[i], x$bins$upper[i])
      if(density_boundary_reflection[i]){
        attr(temp_out, "boundary_reflection") <- TRUE
      }

      out[[i]] <- temp_out
    }

  }else{

    # weightfunction specific stuff
    x_seq     <- .weightfunction_local_cuts(x)
    x_seq_rep <- c(1, sort(rep(2:(length(x_seq)-1), 2)) ,length(x_seq))
    x_val_rep <- sort(rep(1:(length(x_seq)-1), 2))
    if(force_samples | .density.prior_need_samples(x)){
      x_sam  <- rng(x, n_samples)
      x_lCI  <- apply(x_sam, 2, stats::quantile, probs = .025)
      x_uCI  <- apply(x_sam, 2, stats::quantile, probs = .975)
      x_mean <- apply(x_sam, 2, mean)
    }else{
      x_sam  <- NULL
      x_lCI  <- mquant(x, .025)
      x_uCI  <- mquant(x, .975)
      x_mean <- mean(x)
    }

    out <- list(
      call    = call("density", print(x, silent = TRUE)),
      bw      = NULL,
      n       = n_points,
      x       = x_seq[x_seq_rep],
      y       = x_mean[x_val_rep],
      y_lCI   = x_lCI[x_val_rep],
      y_uCI   = x_uCI[x_val_rep],
      samples = x_sam
    )


    class(out) <- c("density", "density.prior", "density.prior.weightfunction")
    attr(out, "x_range") <- c(0, 1)
    attr(out, "y_range") <- c(0, max(1, x_mean, x_lCI, x_uCI, na.rm = TRUE))
  }

  return(out)
}
.density.prior.PETPEESE               <- function(x, x_seq, x_range, n_points, n_samples, force_samples, individual, transformation, transformation_arguments, truncate_end){

  # create either distribution for the parameter or the PET/PEESE function
  if(individual){

    if(is.prior.point(x)){
      out <- .density.prior.point(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments)
    }else if(is.prior.simple(x)){
      out <- .density.prior.simple(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
    }

  }else{

    # get the samples to estimate density / obtain the density directly
    if(force_samples | .density.prior_need_samples(x)){
      x_sam  <- rng(x, n_samples)
      x_med  <- stats::quantile(x_sam, .500)
      x_lCI  <- stats::quantile(x_sam, .025)
      x_uCI  <- stats::quantile(x_sam, .975)
    }else{
      x_med  <- quant(x, .500)
      x_lCI  <- quant(x, .025)
      x_uCI  <- quant(x, .975)
      x_sam  <- NULL
    }


    # transform the output, if requested
    if(!is.null(transformation)){
      x_med   <- .density.prior_transformation_x(x_med,   transformation, transformation_arguments)
      x_lCI   <- .density.prior_transformation_x(x_lCI,   transformation, transformation_arguments)
      x_uCI   <- .density.prior_transformation_x(x_uCI,   transformation, transformation_arguments)
      x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
      if(!is.null(x_sam)){
        x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
      }
    }


    # compute the PET/PEESE
    if(is.prior.PET(x)){
      y_med   = x_med  * x_seq
      y_lCI   = x_lCI  * x_seq
      y_uCI   = x_uCI  * x_seq
    }else if(is.prior.PEESE(x)){
      y_med   = x_med  * x_seq^2
      y_lCI   = x_lCI  * x_seq^2
      y_uCI   = x_uCI  * x_seq^2
    }

    y_lower <- pmin(y_lCI, y_uCI)
    y_upper <- pmax(y_lCI, y_uCI)
    y_lCI   <- y_lower
    y_uCI   <- y_upper


    out <- list(
      call    = call("density", print(x, silent = TRUE)),
      bw      = NULL,
      n       = n_points,
      x       = x_seq,
      y       = y_med,
      y_lCI   = y_lCI,
      y_uCI   = y_uCI,
      samples = x_sam
    )


    class(out) <- c("density", "density.prior", if(is.prior.PET(x)) "density.prior.PET" else if(is.prior.PEESE(x)) "density.prior.PEESE")
    attr(out, "x_range") <- range(x_seq)
    attr(out, "y_range") <- range(out$y, out$y_lCI, out$y_uCI)
  }

  return(out)
}
.density.prior.orthonormal_or_meandif <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  boundary_reflection <- FALSE

  # get the samples to estimate density / obtain the density directly
  if(force_samples | .density.prior_need_samples(x)){

    if(is.na(x$parameters[["K"]]) && !is.null(attr(x, "levels"))){
      x$parameters[["K"]] <- .get_prior_factor_levels(x)
    }else if(is.na(x$parameters[["K"]])){
      x$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    x_sam <- rng(x, n_samples)
    x_sam <- as.vector(x_sam)
    x_density <- .density_kde_boundary(
      x      = x_sam,
      n      = n_points,
      from   = x_range[1],
      to     = x_range[2],
      bounds = c(x$truncation[["lower"]], x$truncation[["upper"]])
    )
    x_seq               <- x_density$x
    x_den               <- x_density$y
    boundary_reflection <- isTRUE(attr(x_density, "boundary_reflection"))

  }else{
    x_den <- mpdf(x, x_seq)
    x_sam <- NULL
  }


  # set the endpoints to zero if they correspond to truncation
  if(truncate_end){
    if(isTRUE(all.equal(x$truncation[["lower"]], x_seq[1])) | x$truncation[["lower"]] >= x_seq[1]){
      x_den <- c(0, x_den)
      x_seq <- c(x_seq[1], x_seq)
    }
    if(isTRUE(all.equal(x$truncation[["upper"]], x_seq[length(x_seq)])) | x$truncation[["upper"]] <= x_seq[length(x_seq)]){
      x_den <- c(x_den, 0)
      x_seq <- c(x_seq, x_seq[length(x_seq)])
    }
  }


  # transform the output, if requested
  if(!is.null(transformation)){
    message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the orthonormal/meandif contrasts to the differences from the mean.")
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
    if(!is.null(x_sam)){
      x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
    }
    x_den   <- .density.prior_transformation_y(x_seq, x_den, transformation, transformation_arguments)
  }


  # create the output object
  out <- list(
    call    = call("density", print(x, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_den,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", if(is.prior.orthonormal(x)) "density.prior.orthonormal" else if(is.prior.meandif(x)) "density.prior.meandif")
  attr(out, "x_range") <- x_range
  attr(out, "y_range") <- c(0, max(x_den))
  if(boundary_reflection){
    attr(out, "boundary_reflection") <- TRUE
  }

  return(out)
}
.density.prior.spike_and_slab         <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  density_variable  <- .density.prior.simple(.get_spike_and_slab_variable(x), x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  density_inclusion <- .density.prior.point(prior(distribution = "spike", parameters = list(location = 0)), x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments)

  inclusion_prob <- mean(.get_spike_and_slab_inclusion(x))
  density_variable$y  <- density_variable[["y"]]  * inclusion_prob
  density_inclusion$y <- density_inclusion[["y"]] * (1 - inclusion_prob)

  attr(density_variable,  "y_range") <- attr(density_variable, "y_range")  * inclusion_prob
  attr(density_inclusion, "y_range") <- attr(density_inclusion, "y_range") * (1 - inclusion_prob)

  # create the output object
  out <- list(
    call      = call("density", print(x, silent = TRUE)),
    variable  = density_variable,
    inclusion = density_inclusion
  )


  class(out) <- c("density", "density.prior.spike_and_slab")
  attr(out, "x_range") <- range(c(attr(density_variable, "x_range"), attr(density_inclusion, "x_range")))
  attr(out, "y_range_variable")  <- attr(density_variable,  "y_range")
  attr(out, "y_range_inclusion") <- attr(density_inclusion, "y_range")

  return(out)
}


#' @title Prior range
#'
#' @description Computes range of a prior
#' distribution (if the prior distribution is
#' unbounded range from \code{quantiles} to
#' \code{1 -quantiles}) is returned.
#'
#' @param x a prior
#' @param quantiles quantile to be returned in
#' case of unbounded distribution.
#' @param ... additional arguments
#' @param na.rm unused
#'
#' @return \code{range.prior} returns a numeric vector of
#' length with a plotting range of a prior distribution.
#'
#' @seealso [prior()]
#' @rdname range.prior
#' @export
range.prior  <- function(x, quantiles = NULL, ..., na.rm = FALSE){

  .check_prior(x)
  if(!is.null(quantiles)){
    check_real(quantiles, "quantiles", upper = 0.5, allow_bound = FALSE)
  }else{
    quantiles <- .range.prior_quantile_default(x)
  }


  x_range <- c(NA, NA)

  if(is.prior.weightfunction(x)){
    return(.weightfunction_range(x, quantiles))
  }
  if(is_prior_phacking(x) || is_prior_bias(x)){
    .selection_prior_stop_unsupported_generic("range", x)
  }
  if(is.prior.ordered(x)){
    return(.prior_ordered_range(x, quantiles = quantiles))
  }

  if(is.infinite(x[["truncation"]][["lower"]])){
    x_range[1] <- min(mquant(x, quantiles), na.rm = TRUE)
  }else{
    x_range[1] <- x[["truncation"]][["lower"]]
  }

  if(is.infinite(x[["truncation"]][["upper"]])){
    x_range[2] <- max(mquant(x, 1 - quantiles), na.rm = TRUE)
  }else{
    x_range[2] <- x[["truncation"]][["upper"]]
  }

  return(x_range)
}



# helper functions
.density_kde_gaussian_height <- function(x, value, bw,
                                         bounds = c(-Inf, Inf)){

  if(length(x) == 0L || length(value) != 1L || !is.finite(value) ||
     length(bw) != 1L || !is.finite(bw) || bw <= 0){
    return(NA_real_)
  }
  if(value < bounds[1L] || value > bounds[2L]){
    return(0)
  }

  height <- mean(stats::dnorm(value, mean = x, sd = bw))
  if(is.finite(bounds[1L])){
    height <- height +
      mean(stats::dnorm(value, mean = 2 * bounds[1L] - x, sd = bw))
  }
  if(is.finite(bounds[2L])){
    height <- height +
      mean(stats::dnorm(value, mean = 2 * bounds[2L] - x, sd = bw))
  }

  sample_range <- range(x)
  if(value < sample_range[1L] || value > sample_range[2L]){
    attr(height, "kde_extrapolation") <- list(
      value = value,
      sample_range = sample_range,
      bandwidth = bw,
      distance_bandwidths = if(value < sample_range[1L]){
        (sample_range[1L] - value) / bw
      }else{
        (value - sample_range[2L]) / bw
      },
      boundary_reflection = any(is.finite(bounds))
    )
  }

  height
}


.density_kde_boundary       <- function(x, n, from = NULL, to = NULL, bounds = c(-Inf, Inf), na.rm = FALSE, ...){

  if(!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds)){
    stop("'bounds' must be a numeric vector of length 2.", call. = FALSE)
  }
  if(bounds[1] > bounds[2]){
    stop("The lower boundary must be lower than or equal to the upper boundary.", call. = FALSE)
  }

  # Estimate the bandwidth on the original sample. The reflected sample is only
  # used to correct leakage at the true support boundaries.
  dots <- list(...)
  density_args <- c(list(x = x, n = n, na.rm = na.rm), dots)
  if(!is.null(from)){
    density_args$from <- from
  }
  if(!is.null(to)){
    density_args$to <- to
  }

  density_base <- do.call(stats::density, density_args)
  attr(density_base, "boundary_reflection") <- FALSE

  if(!any(is.finite(bounds))){
    return(density_base)
  }

  # 'from' and 'to' define only the evaluation grid. 'bounds' defines the true
  # support and therefore where reflection is applied.
  x_ref <- x
  weights_ref <- dots[["weights"]]
  if(na.rm){
    keep <- !is.na(x_ref)
    if(!is.null(weights_ref)){
      keep <- keep & !is.na(weights_ref)
    }
    x_ref <- x_ref[keep]
    if(!is.null(weights_ref)){
      weights_ref <- weights_ref[keep]
    }
  }

  if(is.null(weights_ref)){
    weights_ref <- rep(1 / length(x_ref), length(x_ref))
  }

  x_reflected <- x_ref
  weights_reflected <- weights_ref
  if(is.finite(bounds[1])){
    x_reflected <- c(x_reflected, 2 * bounds[1] - x_ref)
    weights_reflected <- c(weights_reflected, weights_ref)
  }
  if(is.finite(bounds[2])){
    x_reflected <- c(x_reflected, 2 * bounds[2] - x_ref)
    weights_reflected <- c(weights_reflected, weights_ref)
  }

  density_reflected_args <- dots
  density_reflected_args[c("adjust", "bw", "subdensity", "warnWbw", "weights", "width")] <- NULL
  density_reflected_args <- c(
    list(
      x          = x_reflected,
      bw         = density_base$bw,
      n          = n,
      from       = if(is.null(from)) density_base$x[1] else from,
      to         = if(is.null(to)) density_base$x[length(density_base$x)] else to,
      weights    = weights_reflected,
      subdensity = TRUE,
      warnWbw    = FALSE,
      na.rm      = FALSE
    ),
    density_reflected_args
  )

  density_reflected <- do.call(stats::density, density_reflected_args)
  if(is.finite(bounds[1])){
    density_reflected$y[density_reflected$x < bounds[1]] <- 0
  }
  if(is.finite(bounds[2])){
    density_reflected$y[density_reflected$x > bounds[2]] <- 0
  }
  attr(density_reflected, "boundary_reflection") <- TRUE

  return(density_reflected)
}
.density.prior_weightfunction_component_bounds <- function(component){

  switch(
    component$type,
    "point" = c(component$location, component$location),
    "beta"  = c(0, 1),
    "prior" = {
      lower <- component$prior$truncation[["lower"]]
      upper <- component$prior$truncation[["upper"]]
      if(component$scale == "omega"){
        c(lower, upper)
      }else{
        c(exp(lower), exp(upper))
      }
    },
    "one_minus_product_beta" = c(0, 1)
  )
}
.density.prior_weightfunction_components_bounds <- function(components){

  bounds <- do.call(rbind, lapply(components, .density.prior_weightfunction_component_bounds))
  c(min(bounds[,1], na.rm = TRUE), max(bounds[,2], na.rm = TRUE))
}
.density.prior_need_samples   <- function(prior){

  return(FALSE)
}
.density.prior_type           <- function(prior){
  if(is.prior.point(prior)){
    return("point")
  }else if(is.prior.simple(prior)){
    return("simple")
  }else if(is.prior.weightfunction(prior)){
    components <- .weightfunction_marginal_components(prior)
    return(vapply(components, function(component){
      if(component$type == "point") "point" else "simple"
    }, character(1)))
  }else if(is_prior_phacking(prior) || is_prior_bias(prior)){
    .selection_prior_stop_unsupported_generic("density", prior)
  }else if(is.prior.orthonormal(prior)){
    return("orthonormal")
  }else if(is.prior.meandif(prior)){
    return("meandif")
  }
}
.range.prior_quantile_default <- function(prior){

  switch(
    prior[["distribution"]],
    "normal"    = .005,
    "lognormal" = .005,
    "t"         = .010,
    "gamma"     = .010,
    "invgamma"  = .010,
    "moment"    = .010,
    "invmoment" = .010,
    "beta"      = .005,
    "exp"       = .005,
    "uniform"   = .005,
    "point"     = .005,
    "weightfunction" = .005,
    "mnormal"    = .005,
    "mt"         = .010,
    "dirichlet"  = .005
  )

}

# transformation functions
.density.prior_transformation_x         <- function(x, transformation, transformation_arguments = NULL){

  arg <- list(x = x)
  for(i in seq_along(transformation_arguments)){
    arg[[names(transformation_arguments)[i]]] <- transformation_arguments[[i]]
  }

  do.call(.density.prior_transformation_functions(transformation)$fun, arg)
}
.density.prior_transformation_inv_x     <- function(x, transformation, transformation_arguments = NULL){

  arg <- list(x = x)
  for(i in seq_along(transformation_arguments)){
    arg[[names(transformation_arguments)[i]]] <- transformation_arguments[[i]]
  }

  do.call(.density.prior_transformation_functions(transformation)$inv, arg)
}
.density.prior_transformation_y         <- function(x, y, transformation, transformation_arguments = NULL){

  x_inv <- .density.prior_transformation_inv_x(x, transformation, transformation_arguments)
  arg <- list(x = x_inv)
  for(i in seq_along(transformation_arguments)){
    arg[[names(transformation_arguments)[i]]] <- transformation_arguments[[i]]
  }

  y / abs(do.call(.density.prior_transformation_functions(transformation)$jac, arg))
}
.density.prior_transformation_functions <- function(transformation){

  if(is.character(transformation) & length(transformation) == 1){

    return(switch(
      transformation,
      "lin" = list(
        fun = function(x, a = 0, b = 1)a + b * x,
        inv = function(x, a = 0, b = 1)(x - a) / b,
        jac = function(x, a = 0, b = 1)b
      ),
      "exp_lin" = list(
        # Exponential-linear transformation: exp(a + b * log(x))
        # Used for log-intercept unscaling where: intercept_orig = exp(log(intercept_z) * b + a)
        # When a = 0 and b = 1, this is identity: exp(log(x)) = x
        fun = function(x, a = 0, b = 1) exp(a + b * log(x)),
        inv = function(x, a = 0, b = 1) exp((log(x) - a) / b),
        jac = function(x, a = 0, b = 1) b * exp(a + b * log(x)) / x
      ),
      "tanh" = list(
        fun = tanh,
        inv = atanh,
        jac = function(x)1 - tanh(x)^2
      ),
      "exp"  = list(
        fun = exp,
        inv = log,
        jac = exp
      )
    ))

  }else if(is.list(transformation) & length(transformation) == 3 & all(names(transformation) %in% c("fun", "inv", "jac"))){

    return(transformation)

  }else{

    stop("Transformation must be either a character vector of length 1 corresponding to one of known transformations ('lin' = linear, 'exp_lin' = exponential-linear for log-intercept, 'tanh' = hyperbolic tangent, 'exp' = exponential) or a list of three functions (fun = transformation function, inv = inverse transformation, jac = derivative of the transformation).")

  }

}
