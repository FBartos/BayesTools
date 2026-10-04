.plot_data_prior_factor_density_transformed <- function(prior_density_context, samples, parameter, prior_list, n_points, x_range = NULL,
                                                        transformation = NULL, transformation_arguments = NULL,
                                                        transformation_settings = FALSE){

  if(is.null(samples[[parameter]]) || !inherits(samples[[parameter]], "mixed_posteriors.factor")){
    return(NULL)
  }

  # The posterior of an ordered factor is plotted on its levels, without the
  # level that the contrast design fixes at zero; its priors follow.
  sample_metadata <- samples[[parameter]]
  structural_levels <- NULL
  if(any(vapply(prior_list, is.prior.ordered, logical(1)))){
    if(!inherits(sample_metadata, "mixed_posteriors.ordered_transformed")){
      sample_metadata <- transform_factor_samples(samples)[[parameter]]
    }
    design <- as.matrix(.factor_term_design_from_metadata(sample_metadata)$design)
    structural_levels <- rowSums(design != 0) == 0
  }
  factor_weights <- .prior_factor_level_weight_matrix(
    sample_metadata = sample_metadata,
    parameter       = parameter,
    samples         = samples
  )
  level_legends <- .bt_label(
    .bt_draws_label_parts(sample_metadata, parameter),
    style = "plot"
  )
  if(length(level_legends) != nrow(factor_weights)){
    level_legends <- rownames(factor_weights)
  }
  if(!is.null(structural_levels) && length(structural_levels) == nrow(factor_weights)){
    factor_weights <- factor_weights[!structural_levels, , drop = FALSE]
    level_legends <- level_legends[!structural_levels]
  }

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
      attr(level_plot_data[[data_i]], "level_legend") <- level_legends[level_i]
      plot_data[[paste0("level", level_i, "_", names(level_plot_data)[data_i])]] <- level_plot_data[[data_i]]
    }
  }

  .plot_factor_level_universe(plot_data, rownames(factor_weights), level_legends,
    ordered = any(vapply(prior_list, is.prior.ordered, logical(1))))
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
  # the levels of an ordered factor have different priors
  if(any(vapply(prior_list, is.prior.ordered, logical(1)))){
    return(TRUE)
  }

  isTRUE(transform_scaled) && any(vapply(prior_list, is.prior.factor, logical(1)))
}

.plot_data_prior_context_for_plot <- function(prior_density_context, samples, parameter, n_points){

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = samples[[parameter]]
  )

  # Plotted posterior samples are monitored coefficient columns (the JAGS
  # nodes); 'multiply_by' scales only their linear-predictor contribution.
  prior_list <- attr(samples, "prior_list", exact = TRUE)
  raw_priors <- .plot_data_prior_list_without_multiply_by(prior_list)

  if(!raw_priors$changed && !is.null(prior_density_context) &&
     .marginal_posterior_context_matches_condition(prior_density_context, condition_metadata)){
    return(prior_density_context)
  }

  if(is.null(prior_list)){
    return(NULL)
  }

  column_names <- NULL
  if(!is.null(prior_density_context)){
    column_names <- prior_density_context[["column_names"]]
  }

  if(!raw_priors$changed){
    return(.marginal_posterior_prior_density_context(
      samples          = samples,
      prior_list       = prior_list,
      column_names     = column_names,
      n_samples        = max(16L, n_points),
      condition_source = samples[[parameter]]
    ))
  }

  if(is.null(column_names)){
    column_names <- unique(unlist(lapply(names(prior_list), function(parameter_name){
      parameter_prior <- prior_list[[parameter_name]]
      if(!is.prior(parameter_prior)){
        parameter_prior <- parameter_prior[[1]]
      }
      .prior_linear_prior_columns(parameter_name, parameter_prior)
    }), use.names = FALSE))
  }
  n_grid <- if(!is.null(prior_density_context[["n_grid"]])){
    prior_density_context[["n_grid"]]
  }else{
    max(16L, n_points)
  }
  formula_scale <- if(isTRUE(.bt_meta_get(samples, "transform_scaled"))){
    .bt_meta_get(samples, "formula_scale")
  }else{
    NULL
  }

  .prior_density_build_context(
    prior_list       = raw_priors$prior_list,
    column_names     = column_names,
    formula_scale    = formula_scale,
    n_grid           = n_grid,
    conditional      = condition_metadata[["conditional"]],
    conditional_rule = condition_metadata[["conditional_rule"]],
    condition_event  = condition_metadata[["condition_event"]]
  )
}

.plot_data_prior_list_without_multiply_by <- function(prior_list){

  changed <- FALSE
  drop_multiply_by <- function(prior){
    if(!is.null(attr(prior, "multiply_by", exact = TRUE))){
      attr(prior, "multiply_by") <- NULL
      changed <<- TRUE
    }
    prior
  }

  if(is.list(prior_list) && !is.prior(prior_list)){
    for(i in seq_along(prior_list)){
      if(!is.prior(prior_list[[i]])){
        next
      }
      prior <- drop_multiply_by(prior_list[[i]])
      if(is.prior.mixture(prior)){
        for(k in seq_along(prior)){
          if(!is.null(attr(prior[[k]], "multiply_by", exact = TRUE))){
            prior[[k]] <- drop_multiply_by(prior[[k]])
          }
        }
      }
      prior_list[[i]] <- prior
    }
  }

  list(prior_list = prior_list, changed = changed)
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

  prior_density <- .prior_density_from_context(prior_density_context, weights)

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
  prior_density <- .bt_meta_get(samples[[parameter]], "prior_density")
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

# Whether posterior draws declare that they have no prior: a prior list of
# prior_none(), as for parameter_mixed_posterior() draws of a quantity
# without a prior density. Their prior curve is unavailable, not missing.
.plot_data_samples_without_prior <- function(x){

  prior_list <- attr(x, "prior_list", exact = TRUE)
  if(is.prior(prior_list)){
    return(is.prior.none(prior_list))
  }
  is.list(prior_list) && length(prior_list) > 0L &&
    all(vapply(prior_list, is.prior.none, logical(1)))
}

# Warning of an omitted prior curve of draws without a prior density (the
# family of .prior_linear_density_warn_curve_unavailable()).
.plot_data_warn_prior_curve_unavailable <- function(parameter){

  warning(structure(
    class = c("BayesTools_prior_curve_unavailable", "BayesTools_plot_condition",
              "warning", "condition"),
    list(
      message = paste0(
        "The prior density curve of '", parameter, "' is unavailable: its ",
        "posterior draws carry no prior density, as the quantity has no ",
        "deterministic prior-density route. The prior curve is omitted from ",
        "the plot."
      ),
      call = NULL
    )
  ))
}

.plot_data_samples_prior_bounds <- function(prior_list, factor_contrasts = FALSE){

  prior_list_simple <- prior_list[!vapply(prior_list, function(prior){
    is.prior.point(prior) || is.prior.none(prior)
  }, logical(1))]
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
  .bt_ordered_source_require_measure(samples)
  prior_list <- attr(samples, "prior_list")
  posterior_density <- .posterior_density_for_method(.bt_meta_get(samples, "posterior_density"), density_method)
  posterior_atoms <- .posterior_atoms_get(samples)
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  if(is.null(posterior_atoms)){
    .plot_data_stop_unknown_atoms()
  }
  posterior_atoms <- .posterior_atoms_for_column(
    posterior_atoms,
    if(ncol(posterior_atoms$locations) == 1L) 1L else parameter
  )
  if(is.null(posterior_atoms)){
    stop("Simple posterior plotting is unavailable because atom metadata do not identify the requested parameter.", call. = FALSE)
  }
  continuous <- .Savage_Dickey_BF.continuous_posterior(samples, posterior_atoms)
  samples_density <- as.numeric(continuous$samples)
  continuous_mass <- continuous$continuous_mass
  if(nrow(posterior_atoms$locations) > 0L){
    x_points <- as.numeric(posterior_atoms$locations[, 1L])
    y_points <- posterior_atoms$mass
  }
  if(!is.null(x_points) && !is.null(transformation)){
    x_points <- .density.prior_transformation_x(x_points, transformation, transformation_arguments)
  }

  # deal with the densities
  if(!is.null(posterior_density) || continuous_mass > 0){

    if(!is.null(posterior_density)){

      # the stored density is the continuous part; the point masses are the
      # declared posterior atoms
      x_den <- posterior_density[["x"]]
      y_den <- posterior_density[["y"]]

      if(!is.null(transformation)){
        x_den   <- .density.prior_transformation_x(x_den, transformation, transformation_arguments)
        y_den   <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
        samples_density <- .density.prior_transformation_x(samples_density, transformation, transformation_arguments)
      }

    }else if(length(samples_density) < 2L || diff(range(samples_density)) == 0){
      stop(
        "Posterior density is unavailable for declared continuous samples with fewer than two distinct values. ",
        "Provide valid 'posterior_density' metadata (posterior_metadata()) and set 'density_method' to 'precomputed'.",
        call. = FALSE
      )
    }else{

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
      y_den <- density_continuous$y * continuous_mass
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
