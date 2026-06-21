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
#' source columns named \code{tau[1]}, ..., \code{tau[nrow(data)]} aligned to
#' the supplied prediction rows, or the source must provide a
#' \code{parameter_source()} \code{values} function for reconstructing row-wise
#' source values from the posterior samples and supplied prediction data.
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
#' `NULL`, fitted source data from `formula_design` metadata are used.
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
                                  new_levels = NULL){

  check_char(parameter, "parameter")
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

  if(!is.language(formula))
    stop("'formula' must be a formula")
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")

  # extract the posterior distribution
  posterior <- as.matrix(.fit_to_posterior(fit))

  # remove the specified response (would crash the model.frame if not included)
  formula <- .remove_response(formula)
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
      new_levels = new_levels
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
  if(!all(model_terms %in% names(prior_list_formula)))
    stop(paste0("The prior distribution for the ", paste0("'", predictors[!model_terms %in% format_parameter_names(names(prior_list_formula), formula_parameters = parameter, formula_prefix = FALSE)], "'", collapse = ", ")," term is missing in the prior_list."))

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
  predictors_type <- model_terms_type[predictors]

  # check that passed data correspond to the specified priors (factor levels etc...) and set the proper contrasts
  if(any(predictors_type == "factor")){

    # check the proper data input for each factor prior
    for(factor in names(predictors_type[predictors_type == "factor"])){

      # select the corresponding prior in the variable
      this_prior <- prior_list_formula[[factor]]

      if(is.factor(data[,factor])){
        if(all(levels(data[,factor]) %in% .get_prior_factor_level_names(this_prior))){
          # either the formatting is correct, or the supplied levels are a subset of the original levels
          # reformat to check ordering and etc...
          data[,factor] <- factor(data[,factor], levels = .get_prior_factor_level_names(this_prior))
        }else{
          # there are some additional levels
          stop(paste0("Levels specified in the '", factor, "' factor variable do not match the levels used for model specification."))
        }
      }else if(all(unique(data[,factor]) %in% .get_prior_factor_level_names(this_prior))){
        # the variable was not passed as a factor but the values matches the factor levels
        data[,factor] <- factor(data[,factor], levels = .get_prior_factor_level_names(this_prior))
      }else{
        # there are some additional mismatching values
        stop(paste0("Levels specified in the '", factor, "' factor variable do not match the levels used for model specification."))
      }

      # set the contrast
      if(is.prior.orthonormal(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.orthonormal"
      }else if(is.prior.meandif(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.meandif"
      }else if(is.prior.independent(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.independent"
      }else if(is.prior.treatment(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.treatment"
      }
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
  model_matrix <- stats::model.matrix(model_frame, formula = formula, data = data)

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
    if(is.null(data)){
      data <- as.data.frame(fitted_design$model_frame)
    }
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

