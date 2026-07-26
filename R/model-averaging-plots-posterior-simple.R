.plot_data_prior_factor_density_transformed <- function(prior_density_context, samples, parameter, prior_list, n_points, x_range = NULL,
                                                        transformation = NULL, transformation_arguments = NULL,
                                                        transformation_settings = FALSE){

  if(is.null(samples[[parameter]]) || !inherits(samples[[parameter]], "mixed_posteriors.factor")){
    return(NULL)
  }

  factor_weights <- .prior_factor_level_weight_matrix(
    sample_metadata = samples[[parameter]],
    parameter       = parameter,
    samples         = samples
  )

  plot_data <- list()
  for(level_i in seq_len(nrow(factor_weights))){
    weights <- rep(0, length(prior_density_context$column_names))
    names(weights) <- prior_density_context$column_names
    weights[colnames(factor_weights)] <- factor_weights[level_i, ]

    level_density <- .prior_density_from_context(prior_density_context, weights)
    level_plot_data <- .prior_linear_density_to_plot_data(
      level_density,
      n_points                  = n_points,
      x_range                   = x_range,
      transformation            = transformation,
      transformation_arguments  = transformation_arguments,
      transformation_settings   = transformation_settings,
      factor                    = TRUE,
      level                     = level_i,
      level_name                = rownames(factor_weights)[level_i]
    )

    for(data_i in seq_along(level_plot_data)){
      plot_data[[paste0("level", level_i, "_", names(level_plot_data)[data_i])]] <- level_plot_data[[data_i]]
    }
  }

  return(plot_data)
}

.plot_data_prior_has_condition <- function(samples, parameter){

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = samples[[parameter]]
  )

  length(.posterior_density_normalize_condition(condition_metadata[["conditional"]])) > 0L
}

.plot_data_prior_should_use_context <- function(samples, parameter, transform_scaled, prior_list){

  if(.plot_data_prior_has_condition(samples, parameter)){
    return(TRUE)
  }

  isTRUE(transform_scaled) && any(vapply(prior_list, is.prior.factor, logical(1)))
}

.plot_data_prior_context_for_plot <- function(prior_density_context, samples, parameter, n_points){

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = samples[[parameter]]
  )

  if(!is.null(prior_density_context) &&
     .marginal_posterior_context_matches_condition(prior_density_context, condition_metadata)){
    return(prior_density_context)
  }

  prior_list <- attr(samples, "prior_list", exact = TRUE)
  if(is.null(prior_list)){
    return(NULL)
  }

  column_names <- NULL
  if(!is.null(prior_density_context)){
    column_names <- prior_density_context[["column_names"]]
  }

  .marginal_posterior_prior_density_context(
    samples          = samples,
    prior_list       = prior_list,
    column_names     = column_names,
    n_samples        = max(16L, n_points),
    allow_failure    = TRUE,
    condition_source = samples[[parameter]]
  )
}

.plot_data_prior_density_context <- function(prior_density_context, samples, parameter, prior_list, n_points, x_range = NULL,
                                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE){

  prior_density_context <- .plot_data_prior_context_for_plot(
    prior_density_context = prior_density_context,
    samples               = samples,
    parameter             = parameter,
    n_points              = n_points
  )
  if(is.null(prior_density_context)){
    return(NULL)
  }

  if(any(vapply(prior_list, is.prior.factor, logical(1)))){
    return(.plot_data_prior_factor_density_transformed(
      prior_density_context      = prior_density_context,
      samples                   = samples,
      parameter                 = parameter,
      prior_list                = prior_list,
      n_points                  = n_points,
      x_range                   = x_range,
      transformation            = transformation,
      transformation_arguments  = transformation_arguments,
      transformation_settings   = transformation_settings
    ))
  }

  column_names <- prior_density_context[["column_names"]]
  if(is.null(column_names) || !parameter %in% column_names){
    return(NULL)
  }

  weights <- rep(0, length(column_names))
  names(weights) <- column_names
  weights[parameter] <- 1

  prior_density <- tryCatch(
    .prior_density_from_context(prior_density_context, weights),
    error = function(e) NULL
  )
  if(is.null(prior_density)){
    return(NULL)
  }

  .prior_linear_density_to_plot_data(
    prior_density,
    n_points                  = n_points,
    x_range                   = x_range,
    transformation            = transformation,
    transformation_arguments  = transformation_arguments,
    transformation_settings   = transformation_settings
  )
}

.plot_data_attached_prior_density <- function(samples, parameter, n_points,
                                              x_range = NULL,
                                              transformation = NULL,
                                              transformation_arguments = NULL,
                                              transformation_settings = FALSE){

  if(is.null(samples[[parameter]])){
    return(NULL)
  }
  prior_density <- attr(samples[[parameter]], "prior_density", exact = TRUE)
  if(!inherits(prior_density, "prior_linear_density")){
    return(NULL)
  }

  plot_data <- .prior_linear_density_to_plot_data(
    prior_density,
    n_points                  = n_points,
    x_range                   = x_range,
    transformation            = transformation,
    transformation_arguments  = transformation_arguments,
    transformation_settings   = transformation_settings
  )

  if(length(plot_data) == 0L){
    return(NULL)
  }

  return(plot_data)
}

.plot_data_samples_prior_bounds <- function(prior_list, factor_contrasts = FALSE){

  prior_list_simple <- prior_list[!vapply(prior_list, is.prior.point, logical(1))]
  if(length(prior_list_simple) == 0L){
    return(c(-Inf, Inf))
  }

  if(factor_contrasts && any(vapply(prior_list_simple, function(p) is.prior.orthonormal(p) || is.prior.meandif(p), logical(1)))){
    return(c(-Inf, Inf))
  }

  lower <- vapply(prior_list_simple, function(p) p$truncation[["lower"]], numeric(1))
  upper <- vapply(prior_list_simple, function(p) p$truncation[["upper"]], numeric(1))

  c(min(lower), max(upper))
}
.plot_data_samples_density_range <- function(bounds, transformation = NULL){

  from <- if(!is.infinite(bounds[1])){
    if(is.null(transformation)){
      bounds[1]
    }else{
      .representable_interior_value(bounds[1], 1)
    }
  }else{
    NULL
  }
  to <- if(!is.infinite(bounds[2])){
    if(is.null(transformation)){
      bounds[2]
    }else{
      .representable_interior_value(bounds[2], -1)
    }
  }else{
    NULL
  }

  if(!is.null(from) && !is.null(to) && from >= to){
    from <- bounds[1]
    to   <- bounds[2]
  }

  list(from = from, to = to)
}
.plot_data_density_add_boundary_zeros <- function(x_den, y_den, bounds){

  if(is.null(bounds) || !is.numeric(bounds) || length(bounds) != 2L ||
     anyNA(bounds) || length(x_den) == 0L || length(y_den) == 0L){
    return(list(x = x_den, y = y_den))
  }

  if(is.finite(bounds[1]) &&
     (isTRUE(all.equal(bounds[1], x_den[1])) || bounds[1] >= x_den[1])){
    x_den <- c(x_den[1], x_den)
    y_den <- c(0, y_den)
  }

  if(is.finite(bounds[2]) &&
     (isTRUE(all.equal(bounds[2], x_den[length(x_den)])) || bounds[2] <= x_den[length(x_den)])){
    x_den <- c(x_den, x_den[length(x_den)])
    y_den <- c(y_den, 0)
  }

  list(x = x_den, y = y_den)
}

.plot_data_warn_missing_stored_point_masses <- function(){

  warning(
    "Stored posterior density does not declare 'point_masses'; sample-derived point masses are not added. Provide explicit 'point_masses' when atomic posterior mass should be shown.",
    call. = FALSE
  )
}

.plot_data_stored_point_masses <- function(posterior_density, transformation = NULL,
                                           transformation_arguments = NULL){

  point_masses <- posterior_density[["point_masses"]]
  x_points <- point_masses[["x"]]
  y_points <- point_masses[["mass"]]
  if(length(y_points) == 0L){
    return(list(x = NULL, y = NULL))
  }
  if(!is.null(transformation)){
    x_points <- .density.prior_transformation_x(
      x_points,
      transformation,
      transformation_arguments
    )
  }

  list(x = x_points, y = y_points)
}

.plot_data_factor_sample_points_for_levels <- function(sample_point_data, levels,
                                                       level_names){

  out <- list()
  if(length(sample_point_data) == 0L || length(levels) == 0L){
    return(out)
  }

  for(level in levels){
    for(point_i in seq_along(sample_point_data)){
      point_data <- sample_point_data[[point_i]]
      attr(point_data, "level") <- level
      if(length(level_names) >= level){
        attr(point_data, "level_name") <- level_names[[level]]
      }
      out[[paste0("points", level, "_", point_i)]] <- point_data
    }
  }

  out
}
.plot_data_samples.simple         <- function(samples, parameter, n_points, transformation, transformation_arguments, transformation_settings,
                                             density_method = c("KDE", "precomputed")){

  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)
  density_method <- .posterior_density_method(density_method)

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL
  boundary_reflection <- FALSE

  # extract the relevant data
  samples    <- samples[[parameter]]
  prior_list <- attr(samples, "prior_list")
  posterior_density <- .posterior_density_for_method(attr(samples, "posterior_density"), density_method)
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  # deal with spikes
  if(any(sapply(prior_list, is.prior.point))){

    # aggregate samples across spikes
    spikes_simplified <- .simplify_spike_samples(samples, prior_list)

    if(nrow(spikes_simplified) > 0){
      x_points <- spikes_simplified[,"location"]
      y_points <- spikes_simplified[,"probability"]
    }else{
      x_points <- NULL
      y_points <- NULL
    }

    # apply transformations
    if(!is.null(transformation)){
      x_points <- .density.prior_transformation_x(x_points, transformation, transformation_arguments)
    }
  }

  # deal with the densities
  if(any(!sapply(prior_list, is.prior.point))){

    samples_density   <- samples[attr(samples, "models_ind") %in% which(!sapply(prior_list, is.prior.point))]

    if(!is.null(posterior_density)){

      sample_points_available <- !is.null(y_points)
      x_points <- NULL
      y_points <- NULL
      if(.posterior_density_point_masses_declared(posterior_density)){
        stored_points <- .plot_data_stored_point_masses(
          posterior_density,
          transformation,
          transformation_arguments
        )
        x_points <- stored_points[["x"]]
        y_points <- stored_points[["y"]]
      }else if(sample_points_available){
        .plot_data_warn_missing_stored_point_masses()
      }

      x_den <- posterior_density[["x"]]
      y_den <- posterior_density[["y"]]

      if(!is.null(transformation)){
        x_den   <- .density.prior_transformation_x(x_den, transformation, transformation_arguments)
        y_den   <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
        samples_density <- .density.prior_transformation_x(samples_density, transformation, transformation_arguments)
      }

    }else if(length(samples_density) > 0){

      # Keep evaluation range separate from true support so bounded KDEs can
      # reflect at the support while avoiding transformation singularities.
      density_bounds <- .posterior_support_bounds(samples, interval_only = TRUE)
      if(is.null(density_bounds)){
        density_bounds <- .plot_data_samples_prior_bounds(prior_list)
      }
      density_range  <- .plot_data_samples_density_range(density_bounds, transformation)

      # get the density estimate
      density_continuous <- .density_kde_boundary(
        x      = samples_density,
        n      = n_points,
        from   = density_range[["from"]],
        to     = density_range[["to"]],
        bounds = density_bounds
      )
      x_den <- density_continuous$x
      y_den <- density_continuous$y * (length(samples_density) / length(samples))
      boundary_reflection <- isTRUE(attr(density_continuous, "boundary_reflection"))

      # apply transformations
      if(!is.null(transformation)){
        x_den   <- .density.prior_transformation_x(x_den,   transformation, transformation_arguments)
        y_den   <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
        samples_density <- .density.prior_transformation_x(samples_density,   transformation, transformation_arguments)
      }

      if(boundary_reflection){
        plot_bounds <- density_bounds
        if(!is.null(transformation)){
          plot_bounds <- .density.prior_transformation_x(plot_bounds, transformation, transformation_arguments)
        }
        density_plot <- .plot_data_density_add_boundary_zeros(x_den, y_den, plot_bounds)
        x_den <- density_plot$x
        y_den <- density_plot$y
      }

    }
  }


  # create the output object
  out <- list()

  # add continuous densities
  if(!is.null(y_den)){
    out_den    <- list(
      call    = call("density", "mixed samples"),
      bw      = NULL,
      n       = n_points,
      x       = x_den,
      y       = y_den,
      samples = samples_density
    )

    class(out_den) <- c("density", "density.prior", "density.prior.simple")
    attr(out_den, "x_range") <- range(x_den)
    attr(out_den, "y_range") <- c(0, max(y_den))
    if(!is.null(posterior_density)){
      attr(out_den, "posterior_density_method") <- posterior_density[["method"]]
      attr(out_den, "posterior_density_diagnostics") <- posterior_density[["diagnostics"]]
    }
    if(boundary_reflection){
      attr(out_den, "boundary_reflection") <- TRUE
    }

    out[["density"]] <- out_den
  }

  # add spikes
  if(!is.null(y_points)){
    for(i in seq_along(y_points)){
      temp_points <- list(
        call    = call("density", paste0("point", i)),
        bw      = NULL,
        n       = n_points,
        x       = x_points[i],
        y       = y_points[i],
        samples = NULL
      )

      class(temp_points) <- c("density", "density.prior", "density.prior.point")
      attr(temp_points, "x_range") <- range(x_points[i])
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))

      out[[paste0("points",i)]] <- temp_points
    }
  }

  return(out)
}
