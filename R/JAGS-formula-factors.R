.get_parameter_scaling_factor_matrix <- function(term, prior_list, posterior, nrow, ncol){

  if(!is.null(attr(prior_list[[term]], "multiply_by"))){
    if(is.numeric(attr(prior_list[[term]], "multiply_by"))){
      temp_multiply_by <- matrix(attr(prior_list[[term]], "multiply_by"), nrow = nrow, ncol = ncol)
    }else{
      temp_multiply_by <- matrix(posterior[,JAGS_parameter_names(attr(prior_list[[term]], "multiply_by"))], nrow = nrow, ncol = ncol, byrow = TRUE)
    }
  }else{
    temp_multiply_by <- matrix(1, nrow = nrow, ncol = ncol)
  }

  return(temp_multiply_by)
}

.factor_level_list <- function(x){

  level_names <- attr(x, "level_names")
  if(is.null(level_names)){
    factor_terms <- attr(x, "factor_terms")
    if(!is.null(factor_terms) && length(factor_terms) > 1){
      return(NULL)
    }

    n_levels <- attr(x, "levels")
    if(is.null(n_levels)){
      return(NULL)
    }

    if(is.prior.factor(x)){
      level_names <- .get_prior_factor_level_names(x)
    }else if(isTRUE(attr(x, "independent")) ||
             identical(unname(attr(x, "factor_contrasts")),
                       "contr.ordered_cumulative_levels")){
      level_names <- seq_len(n_levels)
    }else{
      level_names <- seq_len(n_levels + 1)
    }
  }

  if(is.list(level_names)){
    factor_terms <- attr(x, "factor_terms")
    if(is.null(factor_terms)){
      factor_terms <- names(level_names)
    }
    if(is.null(factor_terms) || any(!nzchar(factor_terms))){
      factor_terms <- paste0("factor", seq_along(level_names))
    }
    level_names <- level_names[factor_terms]
    names(level_names) <- factor_terms
  }else{
    factor_terms <- attr(x, "factor_terms")
    if(is.null(factor_terms) || length(factor_terms) != 1){
      factor_terms <- ".factor"
    }
    level_names <- setNames(list(level_names), factor_terms)
  }

  level_names <- lapply(level_names, as.character)
  return(level_names)
}

.factor_cell_grid <- function(level_names){

  if(is.null(level_names) || length(level_names) == 0){
    return(data.frame())
  }

  expand.grid(level_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

.factor_cell_labels <- function(level_names){

  level_grid <- .factor_cell_grid(level_names)
  if(ncol(level_grid) == 0){
    return(character(0))
  }

  if(ncol(level_grid) == 1){
    return(as.character(level_grid[[1]]))
  }

  apply(level_grid, 1, function(level_row) {
    paste0(names(level_row), "=", unname(level_row), collapse = ", ")
  })
}

.factor_object_contrast_name <- function(x){

  if(is.prior.independent(x) || isTRUE(attr(x, "independent"))){
    return("contr.independent")
  }else if(is.prior.ordered(x)){
    return(.prior_ordered_contrast_name(x$contrast))
  }else if(is.prior.treatment(x) || isTRUE(attr(x, "treatment"))){
    return("contr.treatment")
  }else if(is.prior.orthonormal(x) || isTRUE(attr(x, "orthonormal"))){
    return("contr.orthonormal")
  }else if(is.prior.meandif(x) || isTRUE(attr(x, "meandif"))){
    return("contr.meandif")
  }

  return(NULL)
}

.factor_contrast_matrix <- function(level_names, contrast){

  switch(
    contrast,
    "contr.treatment"   = stats::contr.treatment(level_names),
    "contr.independent" = contr.independent(level_names),
    "contr.orthonormal" = contr.orthonormal(level_names),
    "contr.meandif"     = contr.meandif(level_names),
    "contr.ordered_cumulative" = contr.ordered_cumulative(level_names),
    "contr.ordered_cumulative_levels" = contr.ordered_cumulative_levels(level_names),
    stop("Unsupported factor contrast '", contrast, "'.", call. = FALSE)
  )
}

.bt_concrete_factor_contrasts <- function(data, factor_names,
                                          context = "Factor design"){

  if(length(factor_names) == 0L){
    return(list())
  }
  out <- lapply(factor_names, function(factor_name){
    if(!factor_name %in% names(data) || !is.factor(data[[factor_name]])){
      stop(
        context, " is missing fitted factor metadata for '",
        factor_name, "'.",
        call. = FALSE
      )
    }
    contrast_matrix <- tryCatch(
      stats::contrasts(data[[factor_name]], contrasts = TRUE),
      error = function(e){
        stop(
          context, " could not resolve the concrete contrast matrix for '",
          factor_name, "': ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    if(is.null(contrast_matrix)){
      stop(
        context, " has no concrete contrast matrix for '",
        factor_name, "'.",
        call. = FALSE
      )
    }
    contrast_matrix <- as.matrix(contrast_matrix)
    if(nrow(contrast_matrix) != nlevels(data[[factor_name]]) ||
       any(!is.finite(contrast_matrix))){
      stop(
        context, " has an invalid concrete contrast matrix for '",
        factor_name, "'.",
        call. = FALSE
      )
    }
    contrast_matrix
  })
  names(out) <- factor_names

  out
}

.bt_model_matrix <- function(model_frame, formula, data = model_frame){

  factor_names <- names(model_frame)[vapply(model_frame, is.factor, logical(1))]
  supported_contrasts <- c(
    "contr.treatment",
    "contr.independent",
    "contr.orthonormal",
    "contr.meandif",
    "contr.ordered_cumulative",
    "contr.ordered_cumulative_levels"
  )

  symbolic_contrasts <- lapply(factor_names, function(factor_name){
    attr(model_frame[[factor_name]], "contrasts", exact = TRUE)
  })
  names(symbolic_contrasts) <- factor_names

  resolved_names <- factor_names[vapply(symbolic_contrasts, function(contrast){
    is.character(contrast) && length(contrast) > 0L &&
      contrast[[1L]] %in% supported_contrasts
  }, logical(1))]
  concrete_names <- factor_names[vapply(
    symbolic_contrasts,
    is.matrix,
    logical(1)
  )]
  contrast_names <- factor_names[
    factor_names %in% c(resolved_names, concrete_names)
  ]
  contrasts_arg <- lapply(contrast_names, function(factor_name){
    contrast <- symbolic_contrasts[[factor_name]]
    if(is.matrix(contrast)){
      return(contrast)
    }
    .factor_contrast_matrix(
      levels(model_frame[[factor_name]]),
      contrast[[1L]]
    )
  })
  names(contrasts_arg) <- contrast_names

  model_matrix <- stats::model.matrix(
    model_frame,
    formula       = formula,
    data          = data,
    contrasts.arg = if(length(contrasts_arg) == 0L) NULL else contrasts_arg
  )

  # Keep the public design metadata symbolic while resolving package-defined
  # contrast functions without relying on the user's search path.
  if(length(resolved_names) > 0L){
    model_contrasts <- attr(model_matrix, "contrasts", exact = TRUE)
    for(factor_name in resolved_names){
      model_contrasts[[factor_name]] <- symbolic_contrasts[[factor_name]][[1L]]
    }
    attr(model_matrix, "contrasts") <- model_contrasts
  }

  return(model_matrix)
}
.bt_validate_model_matrix_finite <- function(model_matrix, context){

  if(any(!is.finite(model_matrix))){
    stop(context, " design matrix contains non-finite values.", call. = FALSE)
  }

  invisible(TRUE)
}

.complete_factor_metadata_prior_list <- function(prior_list){

  if(is.null(prior_list) || length(prior_list) == 0){
    return(prior_list)
  }

  for(parameter in names(prior_list)){
    prior_list[[parameter]] <- .complete_factor_metadata(prior_list[[parameter]], parameter)
  }

  return(prior_list)
}

.copy_missing_factor_metadata <- function(target, source){

  for(metadata_name in .bt_factor_metadata_names){
    if(is.null(attr(target, metadata_name, exact = TRUE)) &&
       !is.null(attr(source, metadata_name, exact = TRUE))){
      attr(target, metadata_name) <- attr(source, metadata_name, exact = TRUE)
    }
  }

  return(target)
}

# Binds a factor prior to the parameter it is fitted as. Every factor prior
# must carry its complete factor metadata (formula factor terms, random-effect
# SD factor priors, and prior_factor_levels() set it); the levels of a factor
# are never inferred from its coefficient count. A factor set outside a
# formula is named after the parameter, and ordered priors are bound to the
# parameter's nodes.
.complete_factor_metadata <- function(x, parameter = NULL){

  if((is.prior.mixture(x) || is.prior.spike_and_slab(x)) && length(x) > 0){
    prior_components <- vapply(x, is.prior, logical(1))
    for(component_i in which(prior_components)){
      x[[component_i]] <- .complete_factor_metadata(x[[component_i]], parameter)
    }

    factor_components <- which(vapply(x, is.prior.factor, logical(1)))
    if(length(factor_components) > 0){
      x <- .copy_missing_factor_metadata(x, x[[factor_components[[1]]]])
      x <- .bt_factor_prior_bind_term(x, parameter)
    }
  }

  if(!is.prior.factor(x)){
    return(x)
  }

  if(!.bt_factor_metadata_complete(x)){
    .bt_stop_incomplete_factor_metadata(parameter)
  }
  x <- .bt_factor_prior_bind_term(x, parameter)

  if(is.prior.ordered(x)){
    x <- .bt_bind_ordered_prior_metadata(x, parameter)
  }

  return(x)
}

.add_factor_metadata_from_named_objects <- function(x, parameter, objects){

  level_names <- .factor_level_list(x)
  if(is.null(level_names)){
    return(x)
  }

  factor_terms <- names(level_names)
  if(is.null(attr(x, "factor_terms"))){
    attr(x, "factor_terms") <- factor_terms
  }

  factor_contrasts <- attr(x, "factor_contrasts")
  if(is.null(factor_contrasts)){
    return(x)
  }else{
    factor_contrasts <- as.character(factor_contrasts)
    if(is.null(names(factor_contrasts))){
      names(factor_contrasts) <- factor_terms[seq_along(factor_contrasts)]
    }
    factor_contrasts <- factor_contrasts[factor_terms]
  }

  attr(x, "factor_contrasts") <- factor_contrasts
  return(x)
}

.factor_term_design_from_formula <- function(formula, data, predictors, predictors_type, term_index, term_components, factor_terms, has_intercept){

  level_names <- lapply(factor_terms, function(factor_term) levels(data[[factor_term]]))
  names(level_names) <- factor_terms
  cell_grid <- .factor_cell_grid(level_names)

  grid_data <- data[rep(1, nrow(cell_grid)), predictors, drop = FALSE]
  rownames(grid_data) <- NULL

  for(predictor in predictors){
    if(predictors_type[[predictor]] == "factor"){
      predictor_values <- if(predictor %in% factor_terms){
        cell_grid[[predictor]]
      }else{
        rep(levels(data[[predictor]])[1], nrow(cell_grid))
      }
      grid_data[[predictor]] <- factor(predictor_values, levels = levels(data[[predictor]]))
      attr(grid_data[[predictor]], "contrasts") <-
        attr(data[[predictor]], "contrasts")
    }else{
      grid_data[[predictor]] <- if(predictor %in% term_components) 1 else 0
    }
  }

  grid_model_frame <- stats::model.frame(formula, data = grid_data)
  grid_model_matrix <- .bt_model_matrix(grid_model_frame, formula = formula, data = grid_data)
  grid_terms_indexes <- attr(grid_model_matrix, "assign")
  if(has_intercept){
    grid_terms_indexes <- grid_terms_indexes + 1
    grid_terms_indexes[1] <- 0
  }

  design <- grid_model_matrix[, grid_terms_indexes == term_index, drop = FALSE]

  return(list(
    design     = unname(design),
    cell_grid  = cell_grid,
    cell_names = .factor_cell_labels(level_names),
    level_names = level_names
  ))
}

.factor_term_design_from_metadata <- function(x){

  factor_design <- attr(x, "factor_design")
  if(!is.null(factor_design)){
    factor_design <- as.matrix(factor_design)
    cell_names <- attr(x, "factor_cell_names")
    level_names <- .factor_level_list(x)
    if(is.null(cell_names)){
      if(is.null(level_names)){
        stop("Factor level names are missing and the factor contrast cannot be transformed.", call. = FALSE)
      }
      cell_names <- .factor_cell_labels(level_names)
    }
    return(list(
      design     = factor_design,
      cell_names = cell_names,
      level_names = level_names
    ))
  }

  level_names <- .factor_level_list(x)
  if(is.null(level_names)){
    stop("Factor level names are missing and the factor contrast cannot be transformed.", call. = FALSE)
  }

  factor_terms <- names(level_names)
  factor_contrasts <- attr(x, "factor_contrasts")
  if(is.null(factor_contrasts)){
    stop("Factor contrast metadata is missing and cannot be inferred.", call. = FALSE)
  }else{
    factor_contrasts <- as.character(factor_contrasts)
    if(is.null(names(factor_contrasts))){
      names(factor_contrasts) <- factor_terms[seq_along(factor_contrasts)]
    }
    factor_contrasts <- factor_contrasts[factor_terms]
  }

  if(any(is.na(factor_contrasts))){
    stop("Factor contrast metadata is incomplete and cannot be inferred.", call. = FALSE)
  }

  contrast_matrices <- lapply(factor_terms, function(factor_term) {
    .factor_contrast_matrix(level_names[[factor_term]], factor_contrasts[[factor_term]])
  })
  names(contrast_matrices) <- factor_terms

  level_grid <- expand.grid(lapply(level_names, seq_along), KEEP.OUT.ATTRS = FALSE)
  coef_grid <- expand.grid(lapply(contrast_matrices, function(contrast_matrix) seq_len(ncol(contrast_matrix))), KEEP.OUT.ATTRS = FALSE)

  design <- matrix(NA_real_, nrow = nrow(level_grid), ncol = nrow(coef_grid))
  for(row_i in seq_len(nrow(level_grid))){
    for(col_i in seq_len(nrow(coef_grid))){
      design[row_i, col_i] <- prod(vapply(factor_terms, function(factor_term) {
        contrast_matrices[[factor_term]][level_grid[[factor_term]][row_i], coef_grid[[factor_term]][col_i]]
      }, numeric(1)))
    }
  }

  return(list(
    design     = design,
    cell_names = .factor_cell_labels(level_names),
    level_names = level_names
  ))
}

.transform_factor_contrast_samples <- function(coefficient_samples, metadata, parameter, transformed_class){

  if(!is.matrix(coefficient_samples)){
    coefficient_samples <- matrix(coefficient_samples, ncol = 1)
  }

  design_info <- .factor_term_design_from_metadata(metadata)
  design <- design_info[["design"]]

  if(ncol(coefficient_samples) != ncol(design)){
    stop(
      "The factor contrast design for '", parameter, "' has ", ncol(design),
      " coefficient columns, but the samples contain ", ncol(coefficient_samples), ".",
      call. = FALSE
    )
  }

  transformed_samples <- coefficient_samples %*% t(design)
  # transformed contrast levels are named by their level cells (treatment
  # levels are the level effects themselves)
  level_parts <- .bt_label_factor_level_parts(
    parameter         = parameter,
    x                 = metadata,
    transformation    = if(identical(transformed_class, "mixed_posteriors.treatment_transformed")){
      "none"
    }else{
      "dif"
    },
    formula_parameter = .transformed_factor_formula_parameter(
      coefficient_samples,
      metadata
    )
  )
  colnames(transformed_samples) <- .bt_label(level_parts, style = "selector")

  transformed_quantities <- .transformed_factor_quantities(
    coefficient_samples = coefficient_samples,
    level_parts         = level_parts,
    design              = design,
    columns             = colnames(transformed_samples)
  )
  posterior_atoms <- .posterior_atoms_get(coefficient_samples)
  old_attributes <- attributes(coefficient_samples)
  old_class <- class(coefficient_samples)
  old_attributes <- old_attributes[
    !names(old_attributes) %in% c(
      "dim", "dimnames", "names", "class", "level_names"
    )
  ]
  attributes(transformed_samples) <- c(attributes(transformed_samples), old_attributes)
  transformed_samples <- .bt_meta_refresh(transformed_samples)
  # the coefficient supports and atoms do not describe the transformed levels
  transformed_samples <- .bt_meta_set(transformed_samples, "support", NULL)
  transformed_samples <- .bt_meta_set(transformed_samples, "atoms", NULL)
  transformed_samples <- .bt_meta_set(transformed_samples, "quantities", transformed_quantities)
  transformed_samples <- .bt_meta_set(transformed_samples, "level_quantities", NULL)
  # the level names of every factor of the term (one factor: its levels, the
  # cell names of the transformed columns)
  attr(transformed_samples, "level_names")       <- if(length(design_info[["level_names"]]) == 1L){
    design_info[["cell_names"]]
  }else{
    design_info[["level_names"]]
  }
  attr(transformed_samples, "factor_cell_names") <- design_info[["cell_names"]]
  if(!is.null(posterior_atoms)){
    transformed_samples <- .posterior_atoms_set(
      transformed_samples,
      .posterior_atoms_linear_transform(
        posterior_atoms,
        design,
        column_names = colnames(transformed_samples)
      )
    )
  }
  class(transformed_samples) <- unique(c(old_class, class(transformed_samples), transformed_class))

  return(transformed_samples)
}

# The formula parameter of transformed factor levels: that of the
# coefficient columns' label parts, of the factor prior, or of the draws.
.transformed_factor_formula_parameter <- function(coefficient_samples, metadata){

  coefficient_quantities <- .bt_meta_get(coefficient_samples, "quantities")
  if(!is.null(coefficient_quantities) && nrow(coefficient_quantities) > 0L){
    return(coefficient_quantities$label_parts[[1L]]$formula_parameter)
  }
  formula_parameter <- .bt_label_formula_parameter(metadata)
  if(nzchar(formula_parameter)){
    return(formula_parameter)
  }

  .bt_label_formula_parameter(coefficient_samples)
}

# The column table of transformed factor levels: each level cell is the
# linear combination of the coefficient columns' fitted coordinates given by
# its design row, and the catalog quantity of that level cell (declared by the
# producer of the draws in 'level_quantities') holds its values. NULL when the
# coefficient columns do not identify their fitted coordinates.
.transformed_factor_quantities <- function(coefficient_samples, level_parts,
                                           design, columns){

  # the coefficient columns are the fitted coordinates of the term in
  # coordinate order (their own names without column metadata)
  coefficient_quantities <- .bt_meta_get(coefficient_samples, "quantities")
  coordinates <- if(is.null(coefficient_quantities)){
    colnames(coefficient_samples)
  }else if(nrow(coefficient_quantities) == ncol(design)){
    vapply(seq_len(nrow(coefficient_quantities)), function(i){
      dependencies <- coefficient_quantities$dependencies[[i]]
      weights <- coefficient_quantities$weights[[i]]
      if(length(dependencies) == 1L && isTRUE(weights == 1)) dependencies else NA_character_
    }, character(1))
  }
  if(length(coordinates) != ncol(design) || anyNA(coordinates) ||
     length(level_parts) != nrow(design)){
    return(NULL)
  }
  level_quantities <- .bt_meta_get(coefficient_samples, "level_quantities")
  quantity_ids <- rep("", nrow(design))
  if(!is.null(level_quantities)){
    cells <- .bt_label(
      .bt_label_parts_update(level_parts, transformation = "none"),
      style = "selector"
    )
    rows <- match(cells, level_quantities$column)
    quantity_ids[!is.na(rows)] <- level_quantities$quantity_id[rows[!is.na(rows)]]
  }
  .bt_draws_quantity_table(
    columns      = columns,
    quantity_ids = quantity_ids,
    dependencies = lapply(seq_len(nrow(design)), function(cell){
      coordinates[design[cell, ] != 0]
    }),
    weights      = lapply(seq_len(nrow(design)), function(cell){
      unname(as.numeric(design[cell, design[cell, ] != 0]))
    }),
    label_parts  = level_parts
  )
}

#' @title Transform factor posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on meandif/orthonormal prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#' @return \code{transform_meandif_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [transform_meandif_samples] [transform_meandif_samples] [transform_orthonormal_samples]
#'
#' @export
transform_factor_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  samples <- transform_meandif_samples(samples)
  samples <- transform_orthonormal_samples(samples)
  samples <- transform_ordered_samples(samples)

  return(samples)
}

#' @title Transform meandif posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on meandif prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#' @return \code{transform_meandif_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [contr.meandif]
#'
#' @export
transform_meandif_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.meandif_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && isTRUE(attr(samples[[i]], "meandif"))){

      meandif_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = meandif_samples,
        metadata            = meandif_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.meandif_transformed"
      )
    }
  }

  return(samples)
}

#' @title Transform orthonomal posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on orthonormal prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#' @return \code{transform_orthonormal_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [contr.orthonormal]
#'
#' @export
transform_orthonormal_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.orthonormal_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && isTRUE(attr(samples[[i]], "orthonormal"))){

      orthonormal_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = orthonormal_samples,
        metadata            = orthonormal_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.orthonormal_transformed"
      )
    }
  }

  return(samples)
}

# not part of transform factor samples (as it's usefull only for marginal effects)
transform_treatment_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.treatment_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && isTRUE(attr(samples[[i]], "treatment"))){

      treatment_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = treatment_samples,
        metadata            = treatment_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.treatment_transformed"
      )
    }
  }

  return(samples)
}

transform_ordered_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.ordered_transformed") &&
       inherits(samples[[i]], "mixed_posteriors.factor") &&
       isTRUE(attr(samples[[i]], "ordered", exact = TRUE))){

      ordered_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = ordered_samples,
        metadata            = ordered_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.ordered_transformed"
      )
    }
  }

  return(samples)
}
