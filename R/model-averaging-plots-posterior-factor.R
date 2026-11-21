.plot_data_samples.weightparameter<- function(samples, parameter, n_points,
                                              transformation = NULL, transformation_arguments = NULL,
                                              density_method = c("KDE", "precomputed")){

  check_list(samples, "samples", check_names = "omega", allow_other = TRUE)
  density_method <- .posterior_density_method(density_method)
  parent <- samples
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

  # Scalar declarations, including an explicit empty marginal, determine
  # the measure. Repeated posterior values never introduce an atom.
  posterior_atom_metadata <- .posterior_atoms_get(samples, allow_partial = TRUE)
  if(is.null(posterior_atom_metadata)){
    .plot_data_stop_unknown_atoms()
  }
  atoms <- .posterior_atoms_for_column(posterior_atom_metadata, parameter)
  if(is.null(atoms)){
    column <- if(is.character(parameter)) parameter else colnames(samples)[parameter]
    .bt_formula_measure_check(samples, "atoms", column)
    .plot_data_stop_unknown_atoms()
  }
  density_bounds <- .posterior_support_bounds(samples, name = parameter, interval_only = TRUE)
  support_hull   <- .posterior_support_bounds(samples, name = parameter, exact_only = FALSE, interval_only = TRUE)
  parameter_label <- if(is.character(parameter)) parameter else colnames(samples)[parameter]
  # an explicit precomputed request resolves the stored curve of this omega
  # while the parent samples and their condition are available
  posterior_density <- if(identical(density_method, "precomputed")){
    .plot_data_selected_posterior_density(parent, samples,
      if(is.character(parameter)) parameter else colnames(samples)[parameter])
  }else NULL
  samples <- samples[, parameter]
  continuous <- .Savage_Dickey_BF.continuous_posterior(samples, atoms)
  samples_density <- as.numeric(continuous$samples)
  continuous_component_mass <- continuous$continuous_mass
  .plot_data_require_precomputed_density(posterior_density, density_method, continuous_component_mass)
  if(nrow(atoms$locations) > 0L){
    x_points <- as.numeric(atoms$locations[, 1L])
    y_points <- atoms$mass
  }
  if(!is.null(x_points) && !is.null(transformation)){
    x_points <- .density.prior_transformation_checked_x(x_points, transformation, transformation_arguments)
  }

  # deal with the densities
  if(!is.null(posterior_density) || continuous_component_mass > 0){

    curve <- .plot_data_posterior_curve(
      posterior_density        = posterior_density,
      draws                    = samples_density,
      mass                     = continuous_component_mass,
      n_points                 = n_points,
      bounds                   = function(){
        # Use component support for reflection and a display range for evaluation.
        if(is.null(density_bounds)) c(-Inf, Inf) else density_bounds
      },
      transformation           = transformation,
      transformation_arguments = transformation_arguments,
      range_bounds             = support_hull
    )
    if(is.null(curve)){
      # no kernel density: the declared point masses are drawn, the curve is omitted
      if(is.null(y_points)){
        stop("Weightfunction posterior density is unavailable for declared continuous samples with fewer than two distinct values.", call. = FALSE)
      }
      .plot_data_warn_posterior_curve_unavailable(parameter_label, continuous_component_mass, samples_density)
    }else{
      x_den               <- curve$x
      y_den               <- curve$y
      samples_density     <- curve$samples
      boundary_reflection <- curve$boundary_reflection
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
    # without a continuous part, the display range is the unit interval of omega
    point_range <- c(0, 1)
    if(is.null(x_den) && !is.null(transformation)){
      point_range <- range(.density.prior_transformation_x(point_range, transformation, transformation_arguments))
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
      attr(temp_points, "x_range") <- if(!is.null(x_den)) range(x_den) else range(c(point_range[is.finite(point_range)], x_points[i]))
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))
      attr(temp_points, "parameter") <- parameter

      out[[paste0("points",i)]] <- temp_points
    }
  }

  return(out)
}
.plot_data_samples.factor         <- function(samples, parameter, n_points, transformation, transformation_arguments, transformation_settings,
                                             density_method = c("KDE", "precomputed"), level = NULL){

  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)
  density_method <- .posterior_density_method(density_method)

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL
  sample_point_data <- list()

  # transform & extract the relevant data
  prior_list <- .plot_factor_priors(samples[[parameter]])
  posterior_density_sources <- .posterior_density_sources(samples, samples[[parameter]])
  density_condition <- .marginal_posterior_condition_metadata(samples, samples[[parameter]])
  posterior_density_conditional <- density_condition$conditional
  posterior_density_conditional_rule <- density_condition$conditional_rule
  posterior_density_condition_key <- density_condition$condition_key

  samples <- .plot_factor_samples(samples, parameter, prior_list, transformation)
  # the legend label of every level: the level text of its level cell
  level_parts   <- .bt_draws_label_parts(samples, parameter)
  level_legends <- .bt_label(level_parts, style = "plot")
  selected <- .plot_factor_level_indices(level, level_legends)

  # Declared posterior atoms are authoritative (as for simple parameters):
  # each level uses the point masses and continuous mass of its own column.
  column_atoms     <- .plot_data_factor_column_atoms(samples, selected)
  level_point_data <- NULL

  # create the output object
  out <- list()

  # deal with spikes
  column_points <- lapply(column_atoms, function(atoms){
    if(!is.null(atoms)) list(x = as.numeric(atoms$locations[, 1L]), y = atoms$mass)
  })
  if(all(vapply(column_points, identical, logical(1), y = column_points[[1L]]))){
    x_points <- column_points[[1L]]$x
    y_points <- column_points[[1L]]$y
  }else{
    level_point_data <- lapply(column_points[selected], function(points){
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
    for(i in selected){

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
      .plot_data_require_precomputed_density(posterior_density, density_method, level_mass)

      curve <- .plot_data_posterior_curve(
        posterior_density        = posterior_density,
        draws                    = level_samples,
        mass                     = level_mass,
        n_points                 = n_points,
        bounds                   = function(){
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
          density_bounds
        },
        transformation           = transformation,
        transformation_arguments = transformation_arguments,
        range_bounds             = .posterior_support_bounds(
          samples,
          name = colnames(samples_density)[i],
          exact_only = FALSE,
          interval_only = TRUE
        )
      )
      if(is.null(curve)){
        # no kernel density: the declared point masses of the level are drawn,
        # the curve of the level is omitted
        if(length(column_atoms[[i]]$mass) == 0L) .plot_data_stop_posterior_curve_unavailable()
        .plot_data_warn_posterior_curve_unavailable(
          colnames(samples_density)[i], level_mass, level_samples)
        next
      }
      x_den               <- curve$x
      y_den               <- curve$y
      level_samples       <- curve$samples
      boundary_reflection <- curve$boundary_reflection

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
    for(k in seq_along(selected)){
      level <- selected[[k]]
      for(point_i in seq_along(level_point_data[[k]])){
        point_data <- level_point_data[[k]][[point_i]]
        attr(point_data, "level")      <- level
        attr(point_data, "level_name") <- colnames(samples)[level]
        attr(point_data, "level_legend") <- level_legends[level]
        level_points[[paste0("points", level, "_", point_i)]] <- point_data
      }
    }
    out <- c(level_points, out)
  }

  out <- .plot_factor_level_universe(out, colnames(samples), level_legends,
    ordered = any(vapply(prior_list, is.prior.ordered, logical(1))))
  if(!is.null(level)) attr(out, "visible_levels") <- selected
  out
}
#' Resolve plotted columns of a factor posterior
#'
#' @inheritParams plot_posterior
#' @param parameter name of the factor-valued element in \code{samples}.
#' @return A character vector of posterior column names on the fitted plotting
#' scale, in display order. Centered contrasts use their level deviations.
#' Ordered contrasts omit the structural reference.
#' @seealso [plot_posterior()]
#' @export
posterior_factor_columns <- function(samples, parameter, level = NULL){

  check_char(parameter, "parameter")
  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)
  if(!inherits(samples[[parameter]], "mixed_posteriors.factor")){
    stop("'parameter' must identify a factor-valued posterior.", call. = FALSE)
  }
  prior_list <- .plot_factor_priors(samples[[parameter]])
  draws <- .plot_factor_samples(samples, parameter, prior_list)
  labels <- .bt_label(.bt_draws_label_parts(draws, parameter), style = "plot")
  colnames(draws)[.plot_factor_level_indices(level, labels)]
}

.plot_factor_priors <- function(draws){

  priors <- attr(draws, "prior_list", exact = TRUE)
  if(is.prior(priors) && !is.prior.mixture(priors) && !is.prior.spike_and_slab(priors)){
    return(list(priors))
  }
  priors
}

.plot_factor_samples <- function(samples, parameter, prior_list, transformation = NULL){

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

  samples[[parameter]]
}

.plot_factor_visible_data <- function(data, visible){

  # Unindexed components are shared by every level (e.g., a common spike).
  if(is.null(data)) return(NULL)
  keep <- vapply(data, function(component){
    index <- attr(component, "level", exact = TRUE)
    is.null(index) || index %in% visible
  }, logical(1))
  result <- data[keep]
  attr(result, "factor_level_universe") <- attr(data, "factor_level_universe", exact = TRUE)
  attr(result, "visible_levels") <- visible
  result
}

.plot_posterior_factor_levels <- function(data, prior_data, level, individual,
                                          plot_type, par_name, dots, dots_prior,
                                          legend, legend_title, legend_labels, legend_position){

  labels <- attr(data, "factor_level_universe", exact = TRUE)$label
  selected <- .plot_factor_level_indices(level, labels)
  dots$plot_type <- plot_type
  dots$par_name <- par_name
  dots$legend <- legend
  dots$legend_title <- legend_title
  dots$legend_labels <- legend_labels
  dots$legend_position <- legend_position

  if(!is.null(prior_data)){
    dots_prior <- do.call(.transfer_dots, c(list(dots_prior), dots))
    dots_prior <- .plot_ordered_overlay_styles(dots_prior, data, dots)
    dots_prior$plot_type <- plot_type
    dots_prior$par_name <- par_name
    dots_prior$legend <- FALSE
    dots_prior$legend_labels <- legend_labels
    dots_prior$hardcode <- plot_type == "ggplot"
  }

  draw <- function(visible, overlay_scale = NULL){
    posterior <- .plot_factor_visible_data(data, visible)
    if(length(visible) == 1L || length(visible) < length(labels)){
      # A shared spike is drawn once, in the first visible level's style.
      universe <- attr(data, "factor_level_universe", exact = TRUE)
      for(i in seq_along(posterior)){
        if(inherits(posterior[[i]], "density.prior.point") &&
           is.null(attr(posterior[[i]], "level", exact = TRUE))){
          attr(posterior[[i]], "level") <- visible[[1L]]
          attr(posterior[[i]], "level_name") <- universe$name[visible[[1L]]]
        }
      }
    }
    dots$plot_data <- posterior
    prior <- .plot_factor_visible_data(prior_data, visible)
    axes <- c(prior, posterior)
    if(isTRUE(dots$add) && plot_type == "ggplot"){
      if(is.null(overlay_scale)){
        # The receiving ggplot supplies its probability-axis mapping at addition.
        return(.bt_ggplot_prior_overlay(function(scale) draw(visible, scale), axes,
          scale_y2 = .plot_scale_y2_resolve(axes, dots)))
      }
      dots$hardcode <- TRUE
      dots$legend <- FALSE
      dots$scale_y2 <- overlay_scale
      dots_prior$scale_y2 <- overlay_scale
    }
    if(!length(prior)){
      plot <- do.call(.plot_prior_list.factor, dots)
      if(!is.null(overlay_scale)) return(plot$geoms)
      return(plot)
    }
    dots_prior$axis_data <- axes
    dots_prior$plot_data <- prior
    dots_prior$add <- isTRUE(dots$add)
    # A fresh figure scales the combined measures. A base overlay takes
    # its scaling from the active figure through the renderer.
    if(!isTRUE(dots$add) ||
       (plot_type == "base" && is.null(.plot_scale_y2_state_current()))){
      scale <- .plot_scale_y2_resolve(axes, dots)
      dots_prior$scale_y2 <- scale
      dots$scale_y2 <- scale
    }
    render_prior <- if(any(vapply(prior, inherits, logical(1), "density.prior.factor"))){
      .plot_prior_list.factor
    }else .plot_prior_list.both
    plot <- do.call(render_prior, dots_prior)
    dots$add <- TRUE
    posterior_plot <- do.call(.plot_prior_list.factor, dots)
    if(!is.null(overlay_scale)) return(c(plot$geoms, posterior_plot$geoms))
    if(plot_type == "ggplot") plot + posterior_plot else invisible(NULL)
  }

  if(!individual) return(draw(selected))
  plots <- lapply(selected, draw)
  stats::setNames(plots, labels[selected])
}

.plot_factor_level_indices <- function(level, labels){

  if(is.null(level)) return(seq_along(labels))
  if(!is.character(level) || !length(level) || anyNA(level) || anyDuplicated(level)){
    stop("'level' must contain one or more distinct factor level names.", call. = FALSE)
  }
  indices <- match(level, labels)
  if(anyNA(indices)){
    stop("Unknown factor level ", paste(shQuote(level[is.na(indices)]), collapse = ", "),
         ". Available levels are ", paste(shQuote(labels), collapse = ", "), ".", call. = FALSE)
  }
  if(any(labels[duplicated(labels)] %in% level)){
    stop("The requested factor level labels are ambiguous.", call. = FALSE)
  }
  indices
}

.plot_data_factor_column_atoms <- function(samples, selected = seq_len(ncol(samples))){

  .bt_ordered_source_require_measure(samples)
  posterior_atoms <- .posterior_atoms_get(samples, allow_partial = TRUE)
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
    if(is.null(column_atoms) && i %in% selected){
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
