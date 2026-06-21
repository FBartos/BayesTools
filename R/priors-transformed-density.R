.generate_transformed_prior_densities <- function(prior_list, column_names, formula_scale = NULL,
                                                  conditional = NULL, conditional_rule = "AND",
                                                  condition_event = NULL,
                                                  n_grid = .prior_linear_density_default_grid(),
                                                  tail_prob = .prior_linear_density_tail_prob()){

  context <- .prior_density_build_context(
    prior_list       = prior_list,
    column_names     = column_names,
    formula_scale    = formula_scale,
    n_grid           = n_grid,
    tail_prob        = tail_prob,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_event  = condition_event
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

#' @title Plot Transformed Prior Densities
#'
#' @description Plots a prior density after applying the same deterministic
#' coefficient transformation used for formula-scale back-transforms. The helper
#' supports continuous densities and point-mass priors in the same plot.
#'
#' @param prior_list named list of prior distributions.
#' @param column_names character vector of coefficient column names defining the
#' fitted parameter space.
#' @param formula_scale optional nested formula-scale metadata, in the same shape
#' as \code{attr(fit, "formula_scale")}.
#' @param parameter coefficient name to plot.
#' @param n_points number of plotting points.
#' @param x_range optional plotting range on the untransformed scale.
#' @param transformation optional output transformation passed to the prior
#' plotting machinery.
#' @param transformation_arguments optional list of transformation arguments.
#' @param transformation_settings whether \code{x_range} is supplied on the
#' transformed scale.
#' @param plot_type either \code{"base"} or \code{"ggplot"}.
#' @param par_name optional plotting label.
#' @param ... additional plotting arguments.
#'
#' @return For \code{plot_type = "base"}, invisibly returns base-plot metadata.
#' For \code{plot_type = "ggplot"}, returns a \code{ggplot} object. Returns
#' \code{NULL} when the requested parameter is an identity transform and no
#' output transformation was requested.
#'
#' @seealso [plot_prior_list()] [transform_scale_samples()]
#' @export
plot_transformed_prior <- function(prior_list, column_names, formula_scale = NULL,
                                   parameter, n_points = 1000, x_range = NULL,
                                   transformation = NULL,
                                   transformation_arguments = NULL,
                                   transformation_settings = FALSE,
                                   plot_type = c("base", "ggplot"),
                                   par_name = NULL, ...){

  check_list(prior_list, "prior_list")
  check_char(column_names, "column_names", check_length = FALSE)
  check_list(formula_scale, "formula_scale", allow_NULL = TRUE)
  check_char(parameter, "parameter")
  check_int(n_points, "n_points", lower = 2)
  check_real(x_range, "x_range", check_length = 2, allow_NULL = TRUE)
  plot_type <- match.arg(plot_type)
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)

  if(!parameter %in% column_names){
    stop("Parameter '", parameter, "' was not found in 'column_names'.", call. = FALSE)
  }

  if(is.null(transformation) &&
     .plot_transformed_prior_is_identity(prior_list, column_names, formula_scale, parameter)){
    return(NULL)
  }

  prior_densities <- .generate_transformed_prior_densities(
    prior_list    = prior_list,
    column_names  = column_names,
    formula_scale = formula_scale
  )

  if(!parameter %in% names(prior_densities)){
    return(NULL)
  }

  plot_data <- .prior_linear_density_to_plot_data(
    prior_densities[[parameter]],
    n_points                 = n_points,
    x_range                  = x_range,
    transformation           = transformation,
    transformation_arguments = transformation_arguments,
    transformation_settings  = transformation_settings
  )

  .plot_prior_list.both(
    plot_data = plot_data,
    plot_type = plot_type,
    par_name  = if(is.null(par_name)) parameter else par_name,
    ...
  )
}

.plot_transformed_prior_is_identity <- function(prior_list, column_names, formula_scale, parameter){

  if(is.null(formula_scale) || length(formula_scale) == 0L){
    return(TRUE)
  }

  context <- .prior_density_build_context(
    prior_list    = prior_list,
    column_names  = column_names,
    formula_scale = formula_scale
  )

  if(!inherits(context, "prior_density_context")){
    return(FALSE)
  }

  weights <- .prior_density_coefficient_weights(column_names, parameter)

  for(transform in context$transforms){
    if(transform$log_intercept && identical(parameter, transform$intercept)){
      return(FALSE)
    }
  }

  transformed_weights <- .prior_density_context_standardized_weights(context, weights)
  weights <- weights[abs(weights) > .prior_linear_density_zero_tol()]

  all_names <- union(names(weights), names(transformed_weights))
  raw <- rep(0, length(all_names))
  names(raw) <- all_names
  transformed <- raw
  raw[names(weights)] <- weights
  transformed[names(transformed_weights)] <- transformed_weights

  isTRUE(all.equal(raw, transformed, tolerance = .prior_linear_density_zero_tol()))
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
