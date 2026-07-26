#' @title Evaluate JAGS formula using posterior samples
#'
#' @description Evaluates a JAGS formula on a posterior distribution obtained
#' from a fitted model. Formula random effects can be evaluated for existing
#' grouping levels when either standardized latent random effects and covariance
#' hyperparameters were monitored via \code{random_monitor(latent = TRUE)}, or
#' the group-level coefficients were monitored via
#' \code{random_monitor(coefficients = TRUE)}.
#' Row-indexed external random-effect SD sources, such as
#' \code{random_sd_source("tau", shape = "row")}, are evaluated from latent
#' random effects only. Model generation automatically monitors the required
#' latent effects for these blocks. The posterior samples must either contain
#' source columns named \code{tau[1]}, ..., \code{tau[N]}. When prediction
#' \code{data} are supplied, \code{fitted_rows} must explicitly map each
#' prediction row to its fitted observation index. Alternatively, the source
#' must provide a
#' \code{parameter_source()} \code{values} function for reconstructing row-wise
#' source values from the posterior samples and supplied prediction data.
#' Literal \code{expression()} terms cannot be reconstructed automatically.
#' Replaying a fitted formula that contains them produces an error rather than
#' silently omitting their contribution. An explicit expression-free formula
#' can still be supplied to evaluate a selected subset of the fitted formula.
#' Inline transformations, offsets, dot expansion, and arbitrary calls are not
#' supported. Create transformed predictors as explicit columns in \code{data}.
#'
#' @param fit model fitted with either \link[runjags]{runjags} posterior
#' samples obtained with \link[rjags]{rjags-package}
#' @param formula formula specifying the right hand side of the assignment (the
#' left hand side is ignored). If `NULL`, the fitted formula stored in
#' `formula_design` metadata is used. If the formula has a
#' \code{"log(intercept)"} attribute set to \code{TRUE}, the intercept values
#' will be log-transformed before computing the linear predictor.
#' @param parameter name of the parameter created with the formula
#' @param data data.frame containing predictors included in the formula. If
#' `NULL`, versioned original-scale fitted source data from `formula_design`
#' metadata are used. Fits without that metadata must be refitted.
#' @param fitted_rows optional integer vector mapping supplied prediction rows
#' to fitted observation indices. It is required whenever `data` is supplied
#' and a selected random-effect block uses a posterior-indexed row source.
#' Reordering and duplicate indices are supported. Callback-computed row
#' sources do not use this mapping.
#' @param prior_list named list of prior distribution of parameters specified
#' within the \code{formula}. If `NULL`, fitted priors from `formula_design`
#' metadata are used.
#' @param formula_target optional formula prediction target. `NULL` preserves
#' the historical safety behavior. `"fixed"` evaluates only the fixed formula
#' contribution. `"conditional"` evaluates fixed effects plus fitted or
#' explicitly generated random-effect contributions.
#' @param blocks optional random-effect block names used with
#' `formula_target = "conditional"`.
#' @param new_levels optional new-level policy used only with
#' `formula_target = "conditional"`. Use a `random_new_levels()` object or one
#' of `"error"`, `"zero"`, or `"sample"`.
#'
#'
#' @return \code{JAGS_evaluate_formula} returns a matrix of the evaluated posterior samples on
#' the supplied data.
#'
#' @seealso [JAGS_fit()] [JAGS_formula()]
#' @export
JAGS_evaluate_formula <- function(fit, formula = NULL, parameter,
                                  data = NULL, prior_list = NULL,
                                  formula_target = NULL, blocks = NULL,
                                  new_levels = NULL, fitted_rows = NULL){

  check_char(parameter, "parameter", allow_NA = FALSE)
  .bt_check_jags_node_name(parameter, "parameter")
  data_supplied <- !is.null(data)
  if(!is.null(fitted_rows) && !data_supplied){
    stop("'fitted_rows' can be supplied only with 'data'.", call. = FALSE)
  }
  formula_target <- .bt_formula_prediction_target(
    formula_target,
    allow_marginal = FALSE,
    context = "JAGS_evaluate_formula()"
  )
  if(!is.null(blocks)){
    check_char(blocks, "blocks", check_length = 0, allow_NA = FALSE)
    if(anyDuplicated(blocks)){
      stop("'blocks' must be unique.", call. = FALSE)
    }
  }
  if(!is.null(blocks) && !identical(formula_target, "conditional")){
    stop("'blocks' can be used only with formula_target = \"conditional\".", call. = FALSE)
  }
  if(!is.null(new_levels) && !identical(formula_target, "conditional")){
    stop("'new_levels' can be used only with formula_target = \"conditional\".", call. = FALSE)
  }
  if(!is.null(new_levels)){
    new_levels <- .bt_random_new_levels_resolve(new_levels)
  }
  replay_fitted_formula <- is.null(formula)
  fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
  resolved_inputs <- .bt_JAGS_evaluate_formula_resolve_inputs(
    fit = fit,
    formula = formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    fitted_design = fitted_design
  )
  formula <- resolved_inputs$formula
  data <- resolved_inputs$data
  prior_list <- resolved_inputs$prior_list
  fitted_design <- resolved_inputs$fitted_design

  if(!inherits(formula, "formula"))
    stop("'formula' must be a formula", call. = FALSE)
  formula_expressions <- .extract_expressions(formula)
  fitted_expressions <- fitted_design$transformed_terms
  if(length(formula_expressions) > 0L ||
     (replay_fitted_formula && length(fitted_expressions) > 0L)){
    stop(
      "JAGS_evaluate_formula() cannot evaluate literal expression() terms. ",
      "Supply an explicit expression-free formula to evaluate a selected subset.",
      call. = FALSE
    )
  }
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  if(!is.null(fitted_rows)){
    check_int(
      fitted_rows,
      "fitted_rows",
      lower = 1L,
      check_length = nrow(data),
      allow_NA = FALSE
    )
  }
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")

  # extract the posterior distribution
  posterior <- as.matrix(.fit_to_posterior(fit))
  .bt_JAGS_evaluate_formula_validate_posterior_names(posterior)

  # remove the specified response (would crash the model.frame if not included)
  formula <- .remove_response(formula)
  .bt_validate_formula_replay_grammar(formula)
  formula_has_random <- .has_random_effects(formula)
  fitted_has_random <- !is.null(fitted_design) &&
    .bt_formula_design_has_any_random_effects(fitted_design)
  if(!is.null(blocks) && !fitted_has_random){
    stop(
      "The fitted formula for parameter '", parameter,
      "' does not include random-effect blocks.",
      call. = FALSE
    )
  }
  if(identical(formula_target, "fixed")){
    formula <- .remove_random_effects(formula)
  }else if(formula_has_random ||
           (identical(formula_target, "conditional") && fitted_has_random)){
    return(.bt_JAGS_evaluate_formula_with_random_effects(
      fit = fit,
      formula = formula,
      parameter = parameter,
      data = data,
      prior_list = prior_list,
      posterior = posterior,
      formula_target = formula_target,
      blocks = blocks,
      new_levels = new_levels,
      fitted_rows = fitted_rows,
      data_supplied = data_supplied,
      replay_fitted_formula = replay_fitted_formula
    ))
  }
  if(is.null(formula_target) &&
     !inherits(fitted_design, "try-error") &&
     .bt_formula_design_has_sampled_random_effects(fitted_design)){
    stop(
      "The fitted formula for parameter '", parameter,
      "' includes random effects. JAGS_evaluate_formula() cannot currently evaluate ",
      "random-effect fits without silently dropping group-level contributions.",
      call. = FALSE
    )
  }
  log_intercept <- isTRUE(attr(formula, "log(intercept)"))
  if(attr(stats::terms(formula), "intercept") == 0){
    formula <- formula_add_intercept(formula)
    if(log_intercept){
      attr(formula, "log(intercept)") <- TRUE
    }
  }

  # select priors corresponding to the prior distribution
  prior_parameter <- sapply(prior_list, function(p) if(is.null(attr(p, "parameter"))) "__none" else attr(p, "parameter"))
  if(!any(parameter %in% unique(prior_parameter)))
    stop("The specified parameter '", parameter, "' was not used in any of the prior distributions.")
  prior_list_formula <- prior_list[prior_parameter == parameter]
  names(prior_list_formula) <- format_parameter_names(names(prior_list_formula), formula_parameters = parameter, formula_prefix = FALSE)

  # extract the terms information from the formula
  formula_terms    <- stats::terms(formula)
  has_intercept    <- attr(formula_terms, "intercept") == 1
  predictors       <- as.character(attr(formula_terms, "variables"))[-1]
  model_terms      <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))

  # check that all predictors have data and prior distribution
  if(!all(predictors %in% colnames(data)))
    stop(paste0("The ", paste0("'", predictors[!predictors %in% colnames(data)], "'", collapse = ", ")," predictor variable is missing in the data."))
  missing_terms <- model_terms[!model_terms %in% names(prior_list_formula)]
  if(length(missing_terms) > 0L)
    stop(paste0("The prior distribution for the ", paste0("'", missing_terms, "'", collapse = ", ")," term is missing in the prior_list."))
  if(log_intercept){
    .bt_validate_formula_log_intercept_prior(prior_list_formula)
  }

  # obtain predictors characteristics -- based on prior distributions used to fit the original model
  # (i.e., do not truest the supplied data -- probably passed by the user)
  model_terms_type <- sapply(model_terms, function(model_term){
    if(model_term == "intercept"){
      return("continuous")
    }else if(is.prior.factor(prior_list_formula[[model_term]]) || inherits(prior_list_formula[[model_term]], "prior.factor_mixture") || inherits(prior_list_formula[[model_term]], "prior.factor_spike_and_slab")){
      return("factor")
    }else if(is.prior.simple(prior_list_formula[[model_term]]) || inherits(prior_list_formula[[model_term]], "prior.simple_mixture") || inherits(prior_list_formula[[model_term]], "prior.simple_spike_and_slab")){
      return("continuous")
    } else {
      stop(paste0("Unrecognized prior distribution for the '", model_term, "' term."))
    }
  })
  predictors_type <- .bt_JAGS_evaluate_predictor_types(
    predictors = predictors,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    prior_list = prior_list_formula,
    fitted_design = fitted_design
  )

  # check that passed data correspond to the specified priors (factor levels etc...) and set the proper contrasts
  if(any(predictors_type == "factor")){

    # check the proper data input for each factor prior
    for(factor in names(predictors_type[predictors_type == "factor"])){

      factor_metadata <- .bt_JAGS_evaluate_factor_metadata(
        predictor = factor,
        prior_list = prior_list_formula,
        fitted_design = fitted_design
      )
      .bt_validate_categorical_level_names(
        factor_metadata$levels,
        factor,
        context = "Fitted factor metadata"
      )
      .bt_validate_categorical_values(
        data[[factor]],
        factor,
        context = "Factor predictor"
      )
      observed_levels <- unique(as.character(data[[factor]]))
      observed_levels <- observed_levels[!is.na(observed_levels)]
      if(any(!observed_levels %in% factor_metadata$levels)){
        stop(paste0("Levels specified in the '", factor, "' factor variable do not match the levels used for model specification."))
      }

      data[[factor]] <- if(factor_metadata$ordered){
        ordered(data[[factor]], levels = factor_metadata$levels)
      }else{
        factor(data[[factor]], levels = factor_metadata$levels)
      }
      stats::contrasts(data[[factor]]) <- factor_metadata$contrast
    }
  }
  if(any(predictors_type == "continuous")){

    # check the proper data input for each continuous prior
    for(continuous in names(predictors_type[predictors_type == "continuous"])){

      # select the corresponding prior in the variable
      this_prior <- prior_list_formula[[continuous]]

      if(is.prior.factor(this_prior)|| is.prior.discrete(this_prior) || is.prior.PET(this_prior) || is.prior.PEESE(this_prior) || is.prior.weightfunction(this_prior)){
        stop(paste0("Unsupported prior distribution defined for '", continuous, "' continuous variable. See '?prior' for details."))
      }
    }

    data <- .bt_apply_formula_scale_to_data(
      fit = fit,
      parameter = parameter,
      data = data,
      predictors_type = predictors_type
    )
  }

  # get the design matrix
  model_frame  <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.pass),
    error = function(e){
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  if(anyNA(model_frame)){
    stop("Formula predictors contain missing values.", call. = FALSE)
  }
  model_matrix <- .bt_model_matrix(model_frame, formula = formula, data = data)
  .bt_validate_model_matrix_finite(model_matrix, "Formula")

  ### evaluate the design matrix on the samples -> output[data, posterior]
  if(has_intercept){

    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    # check for scaling factors
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term = "intercept", prior_list = prior_list_formula, posterior = posterior, nrow = nrow(data), ncol = nrow(posterior))

    # get intercept values and apply log() transformation if log(intercept) attribute is set
    if(is.prior.point(prior_list_formula[["intercept"]])){
      intercept_values <- rep(
        prior_list_formula[["intercept"]]$parameters[["location"]],
        nrow(posterior)
      )
    }else{
      intercept_values <- posterior[, JAGS_parameter_names("intercept", formula_parameter = parameter)]
    }
    if(log_intercept){
      intercept_values <- log(intercept_values)
    }
    output           <- temp_multiply_by * matrix(intercept_values, nrow = nrow(data), ncol = nrow(posterior), byrow = TRUE)

  }else{

    terms_indexes    <- attr(model_matrix, "assign")
    output           <- matrix(0, nrow = nrow(data), ncol = nrow(posterior))

  }

  # add remaining terms (omitting the intercept indexed as NA)
  for(i in unique(terms_indexes[terms_indexes > 0])){

    # subset the model matrix
    temp_data <- model_matrix[,terms_indexes == i,drop = FALSE]

    # get the posterior (unless point prior was used)
    if(is.prior.point(prior_list_formula[[model_terms[i]]])){
      temp_posterior <- matrix(
        prior_list_formula[[model_terms[i]]]$parameters[["location"]],
        nrow = nrow(posterior),
        ncol = if(model_terms_type[i] == "factor") .get_prior_factor_levels(prior_list_formula[[model_terms[i]]]) else 1
      )
    }else{
      temp_posterior <- posterior[,paste0(
        JAGS_parameter_names(model_terms[i], formula_parameter = parameter),
        if(model_terms_type[i] == "factor" && .get_prior_factor_levels(prior_list_formula[[model_terms[i]]]) > 1) paste0("[", 1:.get_prior_factor_levels(prior_list_formula[[model_terms[i]]]), "]"))
        ,drop = FALSE]
    }

    # check for scaling factors
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term = model_terms[i], prior_list = prior_list_formula, posterior = posterior, nrow = nrow(data), ncol = nrow(posterior))

    output <- output + temp_multiply_by * (temp_data %*% t(temp_posterior))

  }

  return(output)
}

.bt_JAGS_evaluate_predictor_types <- function(predictors, model_terms,
                                              model_terms_type, prior_list,
                                              fitted_design){

  if(length(predictors) == 0L){
    return(stats::setNames(character(), character()))
  }

  predictors_type <- stats::setNames(rep(NA_character_, length(predictors)), predictors)
  design_types <- fitted_design$predictor_types
  if(is.character(design_types) && !is.null(names(design_types))){
    matched <- intersect(predictors, names(design_types))
    predictors_type[matched] <- design_types[matched]
  }

  unresolved <- names(predictors_type)[is.na(predictors_type)]
  for(predictor in unresolved){
    main_term <- which(model_terms == predictor)
    if(length(main_term) == 1L){
      predictors_type[[predictor]] <- model_terms_type[[main_term]]
      next
    }

    containing_terms <- vapply(
      model_terms,
      function(model_term){
        predictor %in% strsplit(model_term, ":", fixed = TRUE)[[1L]]
      },
      logical(1)
    )
    candidate_terms <- model_terms[containing_terms]
    candidate_terms <- candidate_terms[candidate_terms != "intercept"]
    if(length(candidate_terms) == 0L){
      next
    }

    candidate_priors <- prior_list[candidate_terms]
    factor_terms <- unique(unlist(lapply(
      candidate_priors,
      function(this_prior) attr(this_prior, "factor_terms", exact = TRUE)
    ), use.names = FALSE))
    candidate_types <- model_terms_type[match(candidate_terms, model_terms)]
    if(any(candidate_types == "factor") && length(factor_terms) == 0L){
      next
    }
    if(predictor %in% factor_terms){
      predictors_type[[predictor]] <- "factor"
    }else{
      predictors_type[[predictor]] <- "continuous"
    }
  }

  if(anyNA(predictors_type) ||
     any(!predictors_type %in% c("continuous", "factor"))){
    invalid <- names(predictors_type)[
      is.na(predictors_type) |
        !predictors_type %in% c("continuous", "factor")
    ]
    stop(
      "Could not determine the fitted type of predictor(s): ",
      paste0("'", invalid, "'", collapse = ", "),
      ". Supply fit metadata created by JAGS_formula().",
      call. = FALSE
    )
  }

  predictors_type
}

.bt_JAGS_evaluate_factor_metadata <- function(predictor, prior_list,
                                              fitted_design){

  candidate_names <- names(prior_list)[vapply(
    seq_along(prior_list),
    function(i){
      identical(names(prior_list)[[i]], predictor) ||
        predictor %in% attr(prior_list[[i]], "factor_terms", exact = TRUE)
    },
    logical(1)
  )]
  candidate_priors <- prior_list[candidate_names]
  this_prior <- if(length(candidate_priors) > 0L) candidate_priors[[1L]] else NULL

  fitted_levels <- fitted_design$xlevels[[predictor]]
  if(is.null(fitted_levels)){
    fitted_levels <- NULL
    for(candidate_name in candidate_names){
      level_names <- attr(prior_list[[candidate_name]], "level_names", exact = TRUE)
      if(is.list(level_names)){
        level_names <- level_names[[predictor]]
      }else if(!identical(candidate_name, predictor)){
        level_names <- NULL
      }
      if(!is.null(level_names)){
        fitted_levels <- level_names
        break
      }
    }
  }
  if(is.null(fitted_levels) || length(fitted_levels) == 0L){
    stop(
      "Could not recover fitted levels for factor predictor '", predictor,
      "'. Supply fit metadata created by JAGS_formula().",
      call. = FALSE
    )
  }

  fitted_factor <- fitted_design$model_frame[[predictor]]
  ordered_factor <- if(is.factor(fitted_factor)){
    is.ordered(fitted_factor)
  }else{
    !is.null(this_prior) && is.prior.ordered(this_prior)
  }

  fitted_contrast <- fitted_design$contrast_matrices[[predictor]]
  if(is.null(fitted_contrast) && is.null(fitted_design) && !is.null(this_prior)){
    factor_contrasts <- attr(this_prior, "factor_contrasts", exact = TRUE)
    contrast_name <- if(
      !is.null(factor_contrasts) &&
      predictor %in% names(factor_contrasts)
    ){
      factor_contrasts[[predictor]]
    }else{
      .factor_object_contrast_name(this_prior)
    }
    if(!is.null(contrast_name)){
      fitted_contrast <- .factor_contrast_matrix(
        fitted_levels,
        contrast_name
      )
    }
  }
  if(is.null(fitted_contrast)){
    stop(
      "Could not recover the concrete fitted contrast matrix for factor predictor '",
      predictor, "'. Supply fit metadata created by JAGS_formula().",
      call. = FALSE
    )
  }

  list(
    levels = as.character(fitted_levels),
    ordered = ordered_factor,
    contrast = fitted_contrast
  )
}

.bt_formula_prediction_target <- function(formula_target,
                                          allow_marginal = TRUE,
                                          context = "Formula prediction"){

  if(is.null(formula_target)){
    return(NULL)
  }
  check_char(formula_target, "formula_target", check_length = 1,
             allow_NULL = FALSE, allow_NA = FALSE)
  allowed <- c("fixed", "conditional", if(isTRUE(allow_marginal)) "marginal")
  if(!formula_target %in% allowed){
    stop(
      context, " supports formula_target = ",
      paste0("'", allowed, "'", collapse = ", "),
      if(!isTRUE(allow_marginal)) ". Use JAGS_predict_formula() for formula_target = 'marginal'." else ".",
      call. = FALSE
    )
  }

  formula_target
}

.bt_JAGS_evaluate_formula_validate_posterior_names <- function(posterior){

  posterior_names <- colnames(posterior)
  if(is.null(posterior_names) || length(posterior_names) != ncol(posterior) ||
     anyNA(posterior_names) || any(!nzchar(posterior_names))){
    stop(
      "Posterior samples used by JAGS_evaluate_formula() must have non-empty column names.",
      call. = FALSE
    )
  }
  if(anyDuplicated(posterior_names)){
    duplicated_names <- unique(posterior_names[duplicated(posterior_names)])
    stop(
      "Posterior samples used by JAGS_evaluate_formula() must have unique ",
      "column names. Duplicated column(s): ",
      paste0(
        "'", duplicated_names[seq_len(min(4L, length(duplicated_names)))], "'",
        collapse = ", "
      ),
      if(length(duplicated_names) > 4L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  invisible(NULL)
}

.bt_JAGS_evaluate_formula_design <- function(fit, parameter){

  fitted_design <- try(JAGS_formula_design(fit, parameter), silent = TRUE)
  if(inherits(fitted_design, "try-error")){
    return(NULL)
  }

  fitted_design
}

.bt_JAGS_evaluate_formula_resolve_inputs <- function(fit, formula,
                                                     parameter, data,
                                                     prior_list,
                                                     fitted_design = NULL){

  if(is.null(fitted_design)){
    fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
  }
  if(!is.null(fitted_design)){
    .bt_validate_formula_design_replay_schema(
      fitted_design,
      context = "JAGS_evaluate_formula()"
    )
  }

  if(is.null(formula)){
    if(is.null(fitted_design)){
      stop("'formula' must be a formula.", call. = FALSE)
    }
    formula <- fitted_design$formula
    # Stored design formulas intentionally drop their original environment.
    # Use the user workspace as a pragmatic evaluation environment for replay.
    environment(formula) <- globalenv()
    if(isTRUE(fitted_design$log_intercept)){
      attr(formula, "log(intercept)") <- TRUE
    }
  }
  if(is.null(data)){
    if(is.null(fitted_design)){
      stop("'data' must be a data.frame.", call. = FALSE)
    }
    data <- fitted_design$source_data
  }
  if(is.null(prior_list)){
    fit_prior_list <- attr(fit, "prior_list", exact = TRUE)
    if(is.list(fit_prior_list) && length(fit_prior_list) > 0L){
      prior_list <- fit_prior_list
    }else if(!is.null(fitted_design) && is.list(fitted_design$prior_list)){
      prior_list <- fitted_design$prior_list
    }else{
      stop("'prior_list' must be supplied.", call. = FALSE)
    }
  }

  list(
    formula = formula,
    data = data,
    prior_list = prior_list,
    fitted_design = fitted_design
  )
}
