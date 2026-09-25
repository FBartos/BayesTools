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
  sample_points_suppressed_by_level <- logical()
  missing_stored_point_masses_warning <- FALSE

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
    samples[[parameter]] <- .plot_data_factor_drop_structural_levels(
      samples[[parameter]]
    )
  }

  samples    <- samples[[parameter]]

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

  # deal with the densities
  samples_density <- samples
  sample_points_suppressed_by_level <- rep(FALSE, ncol(samples_density))

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
        }else if(length(sample_point_data) > 0L || length(level_point_data[[i]]) > 0L){
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
  if(!is.null(level_point_data)){
    levels <- if(length(sample_points_suppressed_by_level) == 0L){
      seq_along(level_point_data)
    }else{
      which(!sample_points_suppressed_by_level)
    }
    level_points <- list()
    for(level in levels){
      for(point_i in seq_along(level_point_data[[level]])){
        point_data <- level_point_data[[level]][[point_i]]
        attr(point_data, "level")      <- level
        attr(point_data, "level_name") <- colnames(samples)[level]
        level_points[[paste0("points", level, "_", point_i)]] <- point_data
      }
    }
    out <- c(level_points, out)
  }

  return(out)
}
# Drops the transformed level columns that the persisted contrast design fixes
# at zero (the reference level of a cumulative ordered contrast); they carry
# no posterior density. Identified from the design, never from the draws.
.plot_data_factor_drop_structural_levels <- function(samples){

  design <- .factor_term_design_from_metadata(samples)$design
  if(nrow(as.matrix(design)) != ncol(samples)){
    stop("The factor design metadata do not match the transformed factor levels.",
         call. = FALSE)
  }
  keep <- rowSums(as.matrix(design) != 0) > 0
  if(all(keep)){
    return(samples)
  }

  out <- samples[, keep, drop = FALSE]
  attributes_kept <- attributes(samples)
  attributes_kept <- attributes_kept[!names(attributes_kept) %in% c(
    "dim", "dimnames", "names", "level_names", "factor_cell_names"
  )]
  attributes(out) <- c(attributes(out), attributes_kept)
  out <- .bt_meta_refresh(out)
  out <- .bt_meta_set(out, "atoms", NULL)
  for(name in c("level_names", "factor_cell_names")){
    value <- attr(samples, name, exact = TRUE)
    if(!is.null(value) && !is.list(value) && length(value) == length(keep)){
      attr(out, name) <- value[keep]
    }
  }
  posterior_atoms <- .posterior_atoms_get(samples)
  if(!is.null(posterior_atoms)){
    out <- .posterior_atoms_set(
      out,
      .posterior_atoms_linear_transform(
        posterior_atoms,
        diag(length(keep))[keep, , drop = FALSE],
        column_names = colnames(out)
      )
    )
  }
  out
}
.plot_data_factor_column_atoms <- function(samples){

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
    x_points <- .density.prior_transformation_x(x_points, transformation, transformation_arguments)
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
