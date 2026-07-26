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
    }else if(isTRUE(attr(x, "independent"))){
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

.factor_contrast_parameter_names <- function(parameter, level_names, cell_names){

  if(is.null(level_names)){
    return(paste0(parameter, "[dif: ", cell_names, "]"))
  }

  level_grid <- .factor_cell_grid(level_names)
  if(nrow(level_grid) != length(cell_names)){
    return(paste0(parameter, "[dif: ", cell_names, "]"))
  }

  if(length(level_names) == 1){
    return(paste0(parameter, "[dif: ", level_grid[[1]], "]"))
  }

  parameter_terms <- strsplit(parameter, "__xXx__", fixed = TRUE)[[1]]
  factor_terms <- names(level_names)
  factor_positions <- vapply(factor_terms, function(factor_term) {
    factor_position <- which(
      parameter_terms == factor_term |
        endsWith(parameter_terms, paste0("_", factor_term))
    )

    if(length(factor_position) == 1){
      return(factor_position)
    }

    return(NA_integer_)
  }, integer(1))

  if(!all(!is.na(factor_positions))){
    return(paste0(parameter, "[dif: ", cell_names, "]"))
  }

  vapply(seq_len(nrow(level_grid)), function(level_i) {
    formatted_terms <- parameter_terms

    for(factor_term in factor_terms){
      factor_position <- factor_positions[[factor_term]]
      prefix_length <- nchar(formatted_terms[[factor_position]]) - nchar(factor_term)
      formatted_terms[[factor_position]] <- paste0(
        substr(formatted_terms[[factor_position]], 1, prefix_length),
        factor_term,
        "[dif: ",
        level_grid[[factor_term]][level_i],
        "]"
      )
    }

    paste0(formatted_terms, collapse = "__xXx__")
  }, character(1))
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

  metadata_names <- c(
    "levels",
    "coefficient_dim",
    "level_names",
    "interaction",
    "interaction_terms",
    "term_components",
    "factor_terms",
    "factor_contrasts",
    "factor_design",
    "factor_cell_names",
    "ordered_metadata"
  )

  for(metadata_name in metadata_names){
    if(is.null(attr(target, metadata_name, exact = TRUE)) &&
       !is.null(attr(source, metadata_name, exact = TRUE))){
      attr(target, metadata_name) <- attr(source, metadata_name, exact = TRUE)
    }
  }

  return(target)
}

.complete_factor_metadata <- function(x, parameter = NULL){

  if((is.prior.mixture(x) || is.prior.spike_and_slab(x)) && length(x) > 0){
    prior_components <- vapply(x, is.prior, logical(1))
    for(component_i in which(prior_components)){
      x[[component_i]] <- .complete_factor_metadata(x[[component_i]], parameter)
    }

    factor_components <- which(vapply(x, is.prior.factor, logical(1)))
    if(length(factor_components) > 0){
      x <- .copy_missing_factor_metadata(x, x[[factor_components[[1]]]])
    }
  }

  if(!is.prior.factor(x)){
    return(x)
  }

  level_names <- .factor_level_list(x)
  if(is.null(level_names) || length(level_names) != 1L){
    return(x)
  }

  factor_terms <- attr(x, "factor_terms", exact = TRUE)
  if(is.null(factor_terms) || length(factor_terms) != 1L ||
     anyNA(factor_terms) || !nzchar(factor_terms)){
    factor_terms <- if(!is.null(parameter) && nzchar(parameter)){
      parameter
    }else{
      ".factor"
    }
    attr(x, "factor_terms") <- factor_terms
  }

  contrast <- .factor_object_contrast_name(x)
  if(is.null(contrast)){
    return(x)
  }

  factor_contrasts <- attr(x, "factor_contrasts", exact = TRUE)
  if(is.null(factor_contrasts)){
    factor_contrasts <- stats::setNames(contrast, factor_terms)
  }else{
    factor_contrasts <- as.character(factor_contrasts)
    if(is.null(names(factor_contrasts))){
      names(factor_contrasts) <- factor_terms[seq_along(factor_contrasts)]
    }
    factor_contrasts <- factor_contrasts[factor_terms]
    if(any(is.na(factor_contrasts))){
      factor_contrasts[is.na(factor_contrasts)] <- contrast
    }
  }
  attr(x, "factor_contrasts") <- factor_contrasts

  if(is.null(attr(x, "factor_design", exact = TRUE)) ||
     is.null(attr(x, "factor_cell_names", exact = TRUE))){
    design_info <- .factor_term_design_from_metadata(x)
    attr(x, "factor_design")     <- design_info[["design"]]
    attr(x, "factor_cell_names") <- design_info[["cell_names"]]
  }

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
  colnames(transformed_samples) <- .factor_contrast_parameter_names(
    parameter = parameter,
    level_names = design_info[["level_names"]],
    cell_names = design_info[["cell_names"]]
  )

  posterior_atoms <- .posterior_atoms_get(coefficient_samples)
  old_attributes <- attributes(coefficient_samples)
  old_class <- class(coefficient_samples)
  old_attributes <- old_attributes[
    !names(old_attributes) %in% c(
      "dim", "dimnames", "names", "class", "level_names",
      "posterior_support", "posterior_atoms"
    )
  ]
  attributes(transformed_samples) <- c(attributes(transformed_samples), old_attributes)
  attr(transformed_samples, "level_names")       <- design_info[["cell_names"]]
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
       isTRUE(attr(samples[[i]], "ordered"))){

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
