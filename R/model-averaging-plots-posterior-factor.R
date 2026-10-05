.plot_data_samples.weightparameter<- function(samples, parameter, n_points){

  check_list(samples, "samples", check_names = "omega", allow_other = TRUE)
  if(!is.null(samples[["omega"]])){
    samples <- samples[["omega"]]
  }else if(!is.null(samples[["bias"]])){
    samples <- samples[["bias"]]
  }else{
    stop("No 'omega' or 'bias' samples found.")
  }

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL
  boundary_reflection <- FALSE

  # extract the relevant data
  prior_list <- attr(samples, "prior_list")
  # the component of each draw: its model in an ensemble, or its component of
  # the weightfunction prior list
  draw_component <- .bt_draws_component(samples)
  posterior_atom_metadata <- .posterior_atoms_get(samples)
  if(is.null(posterior_atom_metadata)){
    .plot_data_stop_unknown_atoms()
  }
  samples    <- samples[,parameter]
  n_samples_total <- length(samples)
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  # One component per original prior entry: the model indicators and the
  # recorded component probabilities index the unmerged prior list.
  context <- .weightfunction_prior_list_context(prior_list, merge = FALSE)
  parameter_ind <- match(parameter, context$omega_names)
  if(is.na(parameter_ind)){
    stop(
      "Weightfunction posterior plotting is unavailable for '", parameter,
      "' because it is not a weight of the weightfunction prior list.",
      call. = FALSE
    )
  }
  components <- .weightfunction_prior_entry_components(
    context,
    parameter_ind
  )

  component_probabilities <- posterior_atom_metadata$component_probabilities
  if(!is.null(component_probabilities)){
    if(length(component_probabilities) != length(components)){
      stop("Recorded weightfunction component probabilities do not match the weightfunction prior list.", call. = FALSE)
    }
    component_probabilities <- component_probabilities /
      sum(component_probabilities)
  }else{
    if(length(draw_component) != n_samples_total ||
       !all(draw_component %in% seq_along(components))){
      stop("Weightfunction model indicators do not match the weightfunction prior list.", call. = FALSE)
    }
    component_probabilities <- tabulate(
      draw_component,
      nbins = length(components)
    ) / n_samples_total
  }
  for(i in seq_along(components)){
    components[[i]]$weight <- component_probabilities[i]
  }
  component_index <- which(component_probabilities > 0)
  components      <- components[component_index]

  point_components <- vapply(
    components,
    function(component) identical(component$type, "point"),
    logical(1)
  )
  if(any(point_components)){
    point_locations <- vapply(
      components[point_components],
      `[[`,
      numeric(1),
      "location"
    )
    point_masses <- vapply(
      components[point_components],
      `[[`,
      numeric(1),
      "weight"
    )
    point_keys <- sprintf("%a", point_locations)
    x_points <- unname(vapply(
      split(point_locations, point_keys),
      function(x) x[1L],
      numeric(1)
    ))
    y_points <- unname(vapply(
      split(point_masses, point_keys),
      sum,
      numeric(1)
    ))
  }
  continuous_component_mass <- sum(vapply(
    components[!point_components],
    `[[`,
    numeric(1),
    "weight"
  ))

  # deal with the densities
  if(any(!point_components)){

    continuous_components <- component_index[!point_components]
    samples_density <- samples[draw_component %in% continuous_components]

    if(length(samples_density) > 0){

      # Use component support for reflection and a display range for evaluation.
      density_components <- components[!point_components]
      density_bounds <- .density.prior_weightfunction_components_bounds(
        density_components
      )
      density_range <- .weightfunction_components_range(
        density_components,
        samples = samples_density
      )

      # get the density estimate
      density_continuous <- .density_kde_boundary(
        x      = samples_density,
        n      = n_points,
        from   = density_range[1],
        to     = density_range[2],
        bounds = density_bounds
      )
      x_den <- density_continuous$x
      y_den <- density_continuous$y * continuous_component_mass
      boundary_reflection <- isTRUE(attr(density_continuous, "boundary_reflection"))

      if(boundary_reflection){
        density_plot <- .plot_data_density_add_boundary_zeros(x_den, y_den, density_bounds)
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
    attr(out_den, "parameter") <- parameter
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
      attr(temp_points, "x_range") <- if(!is.null(x_den)) range(x_den) else range(c(0, 1, x_points[i]))
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))
      attr(temp_points, "parameter") <- parameter

      out[[paste0("points",i)]] <- temp_points
    }
  }

  return(out)
}
.plot_data_samples.factor         <- function(samples, parameter, n_points, transformation, transformation_arguments, transformation_settings,
                                             density_method = c("KDE", "precomputed")){

  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)
  density_method <- .posterior_density_method(density_method)

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL
  sample_point_data <- list()

  # transform & extract the relevant data
  prior_list <- attr(samples[[parameter]], "prior_list")
  posterior_density_sources <- .posterior_density_sources(samples, samples[[parameter]])
  posterior_density_conditional <- .bt_meta_condition(samples[[parameter]], "conditional")
  posterior_density_conditional_rule <- .bt_meta_condition(samples[[parameter]], "conditional_rule")
  posterior_density_condition_key <- .bt_meta_condition(samples[[parameter]], "condition_key")
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  if(any(sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x)))){
    samples  <- transform_factor_samples(samples)
    if(!is.null(transformation)){
      message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the meandif/orthonormal contrasts to the differences from the mean.")
    }
  }else if(any(sapply(prior_list, is.prior.ordered))){
    # Ordered increments are not levels: plot the level effects instead, as
    # for mean-difference contrasts, omitting the structurally zero level.
    samples <- transform_factor_samples(samples)
    samples[[parameter]] <- .transformed_factor_drop_structural_levels(
      samples[[parameter]]
    )
  }

  samples    <- samples[[parameter]]
  # the legend label of every level: the level text of its level cell
  level_parts   <- .bt_draws_label_parts(samples, parameter)
  level_legends <- .bt_label(level_parts, style = "plot")

  # Declared posterior atoms are authoritative (as for simple parameters):
  # each level uses the point masses and continuous mass of its own column.
  column_atoms     <- .plot_data_factor_column_atoms(samples)
  level_point_data <- NULL

  # create the output object
  out <- list()

  # deal with spikes
  column_points <- lapply(column_atoms, function(atoms){
    list(x = as.numeric(atoms$locations[, 1L]), y = atoms$mass)
  })
  if(all(vapply(column_points, identical, logical(1), y = column_points[[1L]]))){
    x_points <- column_points[[1L]]$x
    y_points <- column_points[[1L]]$y
  }else{
    level_point_data <- lapply(column_points, function(points){
      .plot_data_factor_points(points$x, points$y, n_points, transformation, transformation_arguments)
    })
  }
  if(length(y_points) > 0L){
    sample_point_data <- .plot_data_factor_points(x_points, y_points, n_points, transformation, transformation_arguments)
  }

  # deal with the densities; a stored density of a level is its continuous
  # part, the level's point masses are its declared posterior atoms
  samples_density <- samples

  if(nrow(samples_density) > 0){
    for(i in 1:ncol(samples_density)){

      continuous    <- .Savage_Dickey_BF.continuous_posterior(samples[,i], column_atoms[[i]])
      level_samples <- as.numeric(continuous$samples)
      level_mass    <- continuous$continuous_mass
      if(level_mass <= .Machine$double.eps * max(8, length(column_atoms[[i]]$mass))){
        next
      }

      boundary_reflection <- FALSE
      density_aliases <- .plot_data_factor_density_aliases(
        parameter   = parameter,
        samples     = samples,
        sample_name = colnames(samples_density)[i],
        level_i     = i,
        level_parts = level_parts[[i]]
      )
      posterior_density <- .posterior_density_for_method(
        .posterior_density_from_sources(
          sources          = posterior_density_sources,
          aliases          = density_aliases,
          conditional      = posterior_density_conditional,
          conditional_rule = posterior_density_conditional_rule,
          condition_key    = posterior_density_condition_key
        ),
        density_method
      )

      if(!is.null(posterior_density)){

        x_den <- posterior_density[["x"]]
        y_den <- posterior_density[["y"]]
        if(!is.null(transformation)){
          x_den <- .density.prior_transformation_x(x_den, transformation, transformation_arguments)
          y_den <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
          level_samples <- .density.prior_transformation_x(level_samples, transformation, transformation_arguments)
        }

      }else if(length(level_samples) < 2L || diff(range(level_samples)) == 0){
        stop(
          "Posterior density is unavailable for declared continuous samples with fewer than two distinct values. ",
          "Provide valid 'posterior_density' metadata (posterior_metadata()) and set 'density_method' to 'precomputed'.",
          call. = FALSE
        )
      }else{

        # Factor contrasts may be transformed before plotting; in that case
        # the marginal differences no longer have the original bounded support.
        density_bounds <- .posterior_support_bounds(
          samples,
          name = colnames(samples_density)[i],
          interval_only = TRUE
        )
        if(is.null(density_bounds)){
          density_bounds <- .plot_data_samples_prior_bounds(prior_list, factor_contrasts = TRUE)
        }
        density_range  <- .plot_data_samples_density_range(density_bounds, transformation)

        # get the density estimate
        density_continuous <- .density_kde_boundary(
          x      = level_samples,
          n      = n_points,
          from   = density_range[["from"]],
          to     = density_range[["to"]],
          bounds = density_bounds
        )
        x_den <- density_continuous$x
        y_den <- density_continuous$y * level_mass
        boundary_reflection <- isTRUE(attr(density_continuous, "boundary_reflection"))

        # apply transformations
        if(!is.null(transformation)){
          x_den   <- .density.prior_transformation_x(x_den,   transformation, transformation_arguments)
          y_den   <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
          level_samples <- .density.prior_transformation_x(level_samples, transformation, transformation_arguments)
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

      out_den    <- list(
        call    = call("density", "mixed samples"),
        bw      = NULL,
        n       = n_points,
        x       = x_den,
        y       = y_den,
        samples = level_samples
      )

      class(out_den) <- c("density", "density.prior", "density.prior.factor", "density.prior.simple")
      attr(out_den, "x_range")    <- range(x_den)
      attr(out_den, "y_range")    <- c(0, max(y_den))
      attr(out_den, "level")      <- i
      attr(out_den, "level_name") <- colnames(samples_density)[i]
      attr(out_den, "level_legend") <- level_legends[i]
      if(!is.null(posterior_density)){
        attr(out_den, "posterior_density_method") <- posterior_density[["method"]]
        attr(out_den, "posterior_density_diagnostics") <- posterior_density[["diagnostics"]]
      }
      if(boundary_reflection){
        attr(out_den, "boundary_reflection") <- TRUE
      }

      out[[paste0("density", i)]] <- out_den

    }
  }

  if(length(sample_point_data) > 0L){
    out <- c(sample_point_data, out)
  }
  if(!is.null(level_point_data)){
    level_points <- list()
    for(level in seq_along(level_point_data)){
      for(point_i in seq_along(level_point_data[[level]])){
        point_data <- level_point_data[[level]][[point_i]]
        attr(point_data, "level")      <- level
        attr(point_data, "level_name") <- colnames(samples)[level]
        attr(point_data, "level_legend") <- level_legends[level]
        level_points[[paste0("points", level, "_", point_i)]] <- point_data
      }
    }
    out <- c(level_points, out)
  }

  .plot_factor_level_universe(out, colnames(samples), level_legends,
    ordered = any(vapply(prior_list, is.prior.ordered, logical(1))))
}
.plot_data_factor_column_atoms <- function(samples){

  .bt_ordered_source_require_measure(samples)
  posterior_atoms <- .posterior_atoms_get(samples)
  if(is.null(posterior_atoms)){
    .plot_data_stop_unknown_atoms()
  }

  n_atom_columns <- ncol(posterior_atoms$locations)
  if(!n_atom_columns %in% c(1L, ncol(samples))){
    stop("Factor posterior plotting is unavailable because atom metadata do not match the factor columns.", call. = FALSE)
  }

  lapply(seq_len(ncol(samples)), function(i){
    column_atoms <- .posterior_atoms_for_column(
      posterior_atoms,
      if(n_atom_columns == 1L) 1L else i
    )
    if(is.null(column_atoms)){
      stop("Factor posterior plotting is unavailable because atom metadata do not identify the factor columns.", call. = FALSE)
    }
    column_atoms
  })
}
.plot_data_factor_points <- function(x_points, y_points, n_points, transformation = NULL,
                                     transformation_arguments = NULL){

  if(!is.null(transformation) && length(x_points) > 0L){
    x_points <- .density.prior_transformation_checked_x(x_points, transformation, transformation_arguments)
  }

  out <- list()
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

  out
}

# Names under which a precomputed posterior density of one factor level may be
# stored: the column name, and the level text of its level cell (per factor
# and joined over factors, as is and after 'dif: ' for transformed contrast
# levels), bare and after the parameter, from the level's label parts.
.plot_data_factor_density_aliases <- function(parameter, samples, sample_name, level_i,
                                              level_parts = NULL){

  aliases <- .posterior_density_aliases(sample_name)

  if(!is.null(level_parts) && length(level_parts$levels) > 0L){
    level_text <- unname(level_parts$levels)
    if(identical(.bt_label_relation(level_parts), "dif")){
      level_text <- c(level_text, paste0("dif: ", level_text))
    }
    cell_text <- paste0(unname(level_parts$levels), collapse = ", ")
    aliases <- .posterior_density_aliases(
      aliases,
      level_text,
      paste0(parameter, "[", level_text, "]"),
      cell_text,
      paste0(parameter, "[", cell_text, "]")
    )
  }

  level_names <- attr(samples, "level_names", exact = TRUE)
  if(is.list(level_names)){
    level_names <- .factor_cell_labels(level_names)
  }
  if(length(level_names) == ncol(samples)){
    aliases <- .posterior_density_aliases(
      aliases,
      level_names[[level_i]],
      paste0(parameter, "[", level_names[[level_i]], "]")
    )
  }

  factor_cell_names <- attr(samples, "factor_cell_names", exact = TRUE)
  if(length(factor_cell_names) == ncol(samples)){
    aliases <- .posterior_density_aliases(
      aliases,
      factor_cell_names[[level_i]],
      paste0(parameter, "[", factor_cell_names[[level_i]], "]")
    )
  }

  return(aliases)
}
