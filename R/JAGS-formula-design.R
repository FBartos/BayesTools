.formula_scale_point_terms <- function(prior_list, parameter, model_terms){

  fixed_names <- paste0(parameter, "_", model_terms)
  fixed_names <- fixed_names[fixed_names %in% names(prior_list)]
  if(length(fixed_names) == 0L){
    return(stats::setNames(numeric(), character()))
  }

  point_terms <- numeric()
  for(fixed_name in fixed_names){
    this_prior <- prior_list[[fixed_name]]
    if(!is.prior.point(this_prior)){
      next
    }
    location <- this_prior$parameters[["location"]]
    if(!is.numeric(location) || length(location) != 1L ||
       is.na(location) || !is.finite(location)){
      next
    }
    point_terms[[fixed_name]] <- location
  }

  point_terms
}

.JAGS_formula_design_object <- function(parameter, formula, log_intercept,
                                        model_frame, source_data, model_matrix,
                                        raw_column_names, column_names,
                                        predictors, predictors_type,
                                        model_terms, model_terms_type,
                                        prior_list, formula_scale,
                                        expressions, random_effects,
                                        random_effects_compile = NULL,
                                        jags_data_names,
                                        random_allocations = list(),
                                        random_effects_interface = NULL){

  formula_terms <- stats::terms(formula)
  attr(formula_terms, ".Environment") <- emptyenv()
  formula_output <- formula
  environment(formula_output) <- emptyenv()
  model_frame_output <- model_frame
  attr(model_frame_output, "terms") <- formula_terms

  qr_info <- qr(model_matrix)
  aliased <- rep(FALSE, ncol(model_matrix))
  if(qr_info$rank < ncol(model_matrix)){
    aliased[qr_info$pivot[(qr_info$rank + 1L):ncol(model_matrix)]] <- TRUE
  }
  names(aliased) <- colnames(model_matrix)

  factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  xlevels <- lapply(factor_predictors, function(predictor){
    if(predictor %in% names(model_frame) && is.factor(model_frame[[predictor]])){
      levels(model_frame[[predictor]])
    }else{
      NULL
    }
  })
  names(xlevels) <- factor_predictors
  xlevels <- xlevels[!vapply(xlevels, is.null, logical(1))]

  if(is.null(random_effects_compile)){
    random_effects_compile <- .bt_random_effects_compile_resolved(
      sampled = vapply(random_effects, function(term) term$block_name, character(1)),
      marginalized = character(),
      mode = .bt_random_effects_compile_modes_from_terms(random_effects)
    )
  }

  out <- list(
    parameter          = parameter,
    formula            = formula_output,
    log_intercept      = isTRUE(log_intercept),
    model_frame        = model_frame_output,
    source_data        = source_data,
    model_matrix       = model_matrix,
    column_names       = column_names,
    raw_column_names   = raw_column_names,
    assign             = attr(model_matrix, "assign"),
    terms              = formula_terms,
    contrasts          = attr(model_matrix, "contrasts"),
    xlevels            = xlevels,
    predictors         = predictors,
    predictor_types    = predictors_type,
    model_terms        = model_terms,
    model_terms_type   = model_terms_type,
    prior_list         = prior_list,
    formula_scale      = formula_scale,
    rank               = qr_info$rank,
    qr_pivot           = qr_info$pivot,
    aliased            = aliased,
    transformed_terms  = expressions,
    random_effects     = random_effects,
    random_effects_compile = random_effects_compile,
    jags_data_names    = jags_data_names,
    random_allocations = random_allocations,
    random_effects_interface = random_effects_interface
  )
  class(out) <- c("BayesTools_formula_design", "list")

  return(out)
}

.bt_random_effects_predictor_types <- function(random_effects, data){

  if(length(random_effects) == 0L){
    return(stats::setNames(character(), character()))
  }

  scale_terms <- random_effects[
    !vapply(random_effects, function(random_term){
      .bt_random_effect_structure(random_term) %in% c("cs", "hcs", "ar1", "car", "har")
    }, logical(1))
  ]
  if(length(scale_terms) == 0L){
    return(stats::setNames(character(), character()))
  }

  predictors <- unique(unlist(lapply(scale_terms, function(random_term){
    formula_terms <- stats::terms(random_term$term_formula)
    as.character(attr(formula_terms, "variables"))[-1L]
  }), use.names = FALSE))
  if(length(predictors) == 0L){
    return(stats::setNames(character(), character()))
  }

  missing_predictors <- predictors[!predictors %in% colnames(data)]
  if(length(missing_predictors) > 0L){
    stop(
      paste0(
        "The ",
        paste0("'", missing_predictors, "'", collapse = ", "),
        " predictor variable is missing in the data set."
      ),
      call. = FALSE
    )
  }

  vapply(predictors, function(predictor){
    if(is.factor(data[, predictor]) || is.character(data[, predictor])){
      "factor"
    }else{
      "continuous"
    }
  }, character(1))
}

.bt_merge_predictor_types <- function(...){

  predictor_types <- list(...)
  predictor_types <- predictor_types[vapply(predictor_types, length, integer(1)) > 0L]
  if(length(predictor_types) == 0L){
    return(stats::setNames(character(), character()))
  }

  out <- predictor_types[[1L]]
  if(length(predictor_types) == 1L){
    return(out)
  }

  for(i in seq(2L, length(predictor_types))){
    current <- predictor_types[[i]]
    conflicts <- intersect(names(out), names(current))
    conflicts <- conflicts[out[conflicts] != current[conflicts]]
    if(length(conflicts) > 0L){
      stop(
        "Predictor type conflicts between fixed and random-effect formulas for: ",
        paste(conflicts, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    out[setdiff(names(current), names(out))] <- current[setdiff(names(current), names(out))]
  }

  out
}

.bt_should_scale_predictor <- function(formula_scale, predictor){

  if(is.logical(formula_scale) && length(formula_scale) == 1L){
    return(isTRUE(formula_scale))
  }
  if(is.list(formula_scale) && !is.null(formula_scale[[predictor]])){
    return(isTRUE(formula_scale[[predictor]]))
  }

  FALSE
}

.bt_validate_formula_scale <- function(formula_scale, predictor_types){

  if(is.null(formula_scale)){
    return(invisible(NULL))
  }

  if(is.logical(formula_scale)){
    check_bool(formula_scale, "formula_scale", allow_NA = FALSE)
    return(invisible(NULL))
  }

  if(!is.list(formula_scale)){
    stop("'formula_scale' must be NULL, a single TRUE or FALSE, or a named list.", call. = FALSE)
  }

  if(is.null(names(formula_scale)) || anyNA(names(formula_scale)) || any(names(formula_scale) == "")){
    stop("'formula_scale' must be a named list.", call. = FALSE)
  }

  if(any(duplicated(names(formula_scale)))){
    stop("'formula_scale' names must be unique.", call. = FALSE)
  }

  unknown_predictors <- setdiff(names(formula_scale), names(predictor_types))
  if(length(unknown_predictors) > 0L){
    stop(
      "The '",
      paste0(unknown_predictors, collapse = "', '"),
      "' entries in 'formula_scale' are not predictor variables in the formula.",
      call. = FALSE
    )
  }

  for(predictor in names(formula_scale)){
    check_bool(
      formula_scale[[predictor]],
      paste0("formula_scale[['", predictor, "']]"),
      allow_NA = FALSE
    )
  }

  scaled_predictors <- names(formula_scale)[vapply(formula_scale, isTRUE, logical(1))]
  noncontinuous_predictors <- scaled_predictors[predictor_types[scaled_predictors] != "continuous"]
  if(length(noncontinuous_predictors) > 0L){
    stop(
      "Only continuous predictors can be standardized; '",
      paste0(noncontinuous_predictors, collapse = "', '"),
      "' in 'formula_scale' is not continuous.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' @title Extract Fitted JAGS Formula Design Metadata
#'
#' @description Returns the fitted formula design metadata stored by
#' [JAGS_fit()]. The design contains the processed formula, fitted model frame,
#' exact model matrix used for JAGS data construction, JAGS-safe coefficient
#' names, contrast and factor-level metadata, rank diagnostics, prior metadata,
#' and formula-scale information.
#'
#' @param fit a fitted object returned by [JAGS_fit()].
#' @param parameter optional formula parameter name. If \code{NULL}, all stored
#' formula designs are returned.
#'
#' @return A named list of formula designs, or one formula design when
#' \code{parameter} is supplied. Returns \code{NULL} when no formula design
#' metadata is stored on \code{fit}.
#'
#' @seealso [JAGS_fit()] [JAGS_formula()]
#' @export
JAGS_formula_design <- function(fit, parameter = NULL){

  check_char(parameter, "parameter", allow_NULL = TRUE)

  formula_design <- attr(fit, "formula_design")
  if(is.null(formula_design)){
    return(NULL)
  }

  if(is.null(parameter)){
    return(formula_design)
  }

  if(!parameter %in% names(formula_design)){
    stop("Formula design for parameter '", parameter, "' was not found.", call. = FALSE)
  }

  return(formula_design[[parameter]])
}

