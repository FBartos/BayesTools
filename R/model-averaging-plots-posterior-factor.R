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
  models_ind <- attr(samples, "models_ind")
  posterior_atoms <- .posterior_atoms_for_column(
    .posterior_atoms_get(samples),
    parameter
  )
  posterior_atom_metadata <- .posterior_atoms_get(samples)
  samples    <- samples[,parameter]
  n_samples_total <- length(samples)
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  context <- .weightfunction_prior_list_context(prior_list)
  parameter_ind <- match(parameter, context$omega_names)
  components <- NULL
  continuous_component_mass <- NULL

  if(!is.na(parameter_ind)){
    components <- .weightfunction_prior_marginal_components(
      context,
      parameter_ind
    )

    component_probabilities <- if(!is.null(posterior_atom_metadata)){
      posterior_atom_metadata$component_probabilities
    }else{
      NULL
    }
    if(!is.null(component_probabilities) &&
       length(component_probabilities) == length(components)){
      component_probabilities <- component_probabilities /
        sum(component_probabilities)
      for(i in seq_along(components)){
        components[[i]]$weight <- component_probabilities[i]
      }
    }else if(length(models_ind) == n_samples_total &&
             length(components) > 0L &&
             all(models_ind %in% seq_along(components))){
      component_probabilities <- tabulate(
        models_ind,
        nbins = length(components)
      ) / n_samples_total
      for(i in seq_along(components)){
        components[[i]]$weight <- component_probabilities[i]
      }
    }

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
  }else if(!is.null(posterior_atoms) &&
           nrow(posterior_atoms$locations) > 0L){
    # Retain explicit atom metadata for nonstandard legacy weight names.
    x_points <- posterior_atoms$locations[, 1L]
    y_points <- posterior_atoms$mass
  }

  # deal with the densities
  has_continuous_component <- if(is.null(components)){
    !all(sapply(prior_list, \(x) is.prior.point(x) || is.prior.none(x)))
  }else{
    any(vapply(
      components,
      function(component) !identical(component$type, "point"),
      logical(1)
    ))
  }
  if(has_continuous_component){

    if(is.null(components)){
      continuous_components <- which(!sapply(prior_list, is.prior.point))
    }else{
      continuous_components <- which(vapply(
        components,
        function(component) !identical(component$type, "point"),
        logical(1)
      ))
    }
    samples_density <- samples[models_ind %in% continuous_components]

    if(length(samples_density) > 0){

      # Use component support for reflection and a display range for evaluation.
      if(is.na(parameter_ind)){
        density_bounds <- c(-Inf, Inf)
        density_range <- range(c(0, 1, samples_density), finite = TRUE)
      }else{
        density_components <- components[continuous_components]
        density_bounds <- .density.prior_weightfunction_components_bounds(
          density_components
        )
        density_range <- .weightfunction_components_range(
          density_components,
          samples = samples_density
        )
      }

      # get the density estimate
      density_continuous <- .density_kde_boundary(
        x      = samples_density,
        n      = n_points,
        from   = density_range[1],
        to     = density_range[2],
        bounds = density_bounds
      )
      x_den <- density_continuous$x
      density_mass <- if(is.null(continuous_component_mass)){
        length(samples_density) / n_samples_total
      }else{
        continuous_component_mass
      }
      y_den <- density_continuous$y * density_mass
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
  sample_points_suppressed_by_level <- logical()
  missing_stored_point_masses_warning <- FALSE

  # transform & extract the relevant data
  prior_list <- attr(samples[[parameter]], "prior_list")
  posterior_density_sources <- .posterior_density_sources(samples, samples[[parameter]])
  posterior_density_conditional <- attr(samples[[parameter]], "conditional", exact = TRUE)
  posterior_density_conditional_rule <- attr(samples[[parameter]], "conditional_rule", exact = TRUE)
  posterior_density_condition_key <- attr(samples[[parameter]], "condition_key", exact = TRUE)
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  if(any(sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x)))){
    samples  <- transform_factor_samples(samples)
    if(!is.null(transformation)){
      message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the meandif/orthonormal contrasts to the differences from the mean.")
    }
  }

  samples    <- samples[[parameter]]

  # create the output object
  out <- list()

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

      sample_point_data[[paste0("points",i)]] <- temp_points
    }
  }

  # deal with the densities
  if(any(!sapply(prior_list, is.prior.point))){

    samples_density <- samples[attr(samples, "models_ind") %in% which(!sapply(prior_list, is.prior.point)),,drop=FALSE]
    sample_points_suppressed_by_level <- rep(FALSE, ncol(samples_density))

    if(nrow(samples_density) > 0){
      for(i in 1:ncol(samples_density)){

        boundary_reflection <- FALSE
        density_aliases <- .plot_data_factor_density_aliases(
          parameter   = parameter,
          samples     = samples,
          sample_name = colnames(samples_density)[i],
          level_i     = i
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

          if(.posterior_density_point_masses_declared(posterior_density)){
            sample_points_suppressed_by_level[i] <- TRUE
            stored_points <- .plot_data_stored_point_masses(
              posterior_density,
              transformation,
              transformation_arguments
            )
            x_points_i <- stored_points[["x"]]
            y_points_i <- stored_points[["y"]]
            for(point_i in seq_along(y_points_i)){
              temp_points <- list(
                call    = call("density", paste0("point", point_i)),
                bw      = NULL,
                n       = n_points,
                x       = x_points_i[point_i],
                y       = y_points_i[point_i],
                samples = NULL
              )

              class(temp_points) <- c("density", "density.prior", "density.prior.point")
              attr(temp_points, "x_range")    <- range(x_points_i[point_i])
              attr(temp_points, "y_range")    <- c(0, max(y_points_i[point_i]))
              attr(temp_points, "level")      <- i
              attr(temp_points, "level_name") <- colnames(samples_density)[i]

              out[[paste0("density", i, "_points", point_i)]] <- temp_points
            }
          }else if(length(sample_point_data) > 0L){
            sample_points_suppressed_by_level[i] <- TRUE
            if(!missing_stored_point_masses_warning){
              .plot_data_warn_missing_stored_point_masses()
              missing_stored_point_masses_warning <- TRUE
            }
          }

          x_den <- posterior_density[["x"]]
          y_den <- posterior_density[["y"]]
          if(!is.null(transformation)){
            x_den <- .density.prior_transformation_x(x_den, transformation, transformation_arguments)
            y_den <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
            samples_density[,i] <- .density.prior_transformation_x(samples_density[,i], transformation, transformation_arguments)
          }

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
            x      = samples_density[,i],
            n      = n_points,
            from   = density_range[["from"]],
            to     = density_range[["to"]],
            bounds = density_bounds
          )
          x_den <- density_continuous$x
          y_den <- density_continuous$y * (nrow(samples_density) / nrow(samples))
          boundary_reflection <- isTRUE(attr(density_continuous, "boundary_reflection"))

          # apply transformations
          if(!is.null(transformation)){
            x_den   <- .density.prior_transformation_x(x_den,   transformation, transformation_arguments)
            y_den   <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
            samples_density[,i] <- .density.prior_transformation_x(samples_density[,i],   transformation, transformation_arguments)
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
          samples = samples_density
        )

        class(out_den) <- c("density", "density.prior", "density.prior.factor", "density.prior.simple")
        attr(out_den, "x_range")    <- range(x_den)
        attr(out_den, "y_range")    <- c(0, max(y_den))
        attr(out_den, "level")      <- i
        attr(out_den, "level_name") <- colnames(samples_density)[i]
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
  }

  if(length(sample_point_data) > 0L){
    if(length(sample_points_suppressed_by_level) == 0L ||
       !any(sample_points_suppressed_by_level)){
      out <- c(sample_point_data, out)
    }else{
      out <- c(
        .plot_data_factor_sample_points_for_levels(
          sample_point_data,
          which(!sample_points_suppressed_by_level),
          colnames(samples_density)
        ),
        out
      )
    }
  }



  return(out)
}

.plot_data_factor_density_aliases <- function(parameter, samples, sample_name, level_i){

  bracket_matches <- regmatches(sample_name, gregexpr("\\[[^]]+\\]", sample_name))[[1]]
  bracket_aliases <- gsub("^\\[|\\]$", "", bracket_matches)
  aliases <- .posterior_density_aliases(sample_name)

  if(length(bracket_aliases) > 0L){
    aliases <- .posterior_density_aliases(
      aliases,
      bracket_aliases,
      paste0(parameter, "[", bracket_aliases, "]")
    )
  }
  if(length(bracket_aliases) > 1L){
    cell_alias <- paste0(bracket_aliases, collapse = ", ")
    aliases <- .posterior_density_aliases(
      aliases,
      cell_alias,
      paste0(parameter, "[", cell_alias, "]")
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
