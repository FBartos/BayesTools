#' @title Model-average marginal posterior distributions and
#' marginal Bayes factors
#'
#' @description Creates marginal model-averaged and conditional
#' posterior distributions based on a list of models, vector of parameters,
#' formula, and a list of indicators of the null or alternative hypothesis models
#' for each parameter. Computes inclusion Bayes factors for each
#' marginal estimate via a Savage-Dickey density approximation.
#'
#' @param marginal_parameters parameters for which the the marginal summary
#' should be created
#' @param parameters all parameters included in the model_list that are
#' relevant for the formula (all of which need to have specification of
#' \code{is_null_list})
#' @param seed seed for random number generation
#' @param density_method posterior density method used for Savage-Dickey Bayes
#' factors. Currently only \code{"KDE"} is supported by
#' \code{marginal_inference()} and \code{as_marginal_inference()} because
#' generic model-averaged marginal posteriors do not have a well-defined
#' precomputed density source.
#' @inheritParams ensemble_inference
#' @inheritParams marginal_posterior
#' @inheritParams Savage_Dickey_BF
#'
#' @return \code{marginal_inference} returns an object of class 'marginal_inference'.
#'
#' @seealso [ensemble_inference] [mix_posteriors] [BayesTools_ensemble_tables]
#'
#' @export
marginal_inference <- function(model_list, marginal_parameters, parameters, is_null_list, formula,
                               null_hypothesis = 0, normal_approximation = FALSE,
                               n_samples = 10000, seed = NULL, silent = FALSE,
                               density_method = "KDE"){

  # check input (majority of the checks performed within mix_posteriors)
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_char(marginal_parameters, "marginal_parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  if(!all(unlist(sapply(model_list, function(m) sapply(attr(m[["fit"]], "prior_list"), function(p) is.prior(p))))))
    stop("model_list:priors must contain 'BayesTools' priors")
  density_method <- .marginal_inference_density_method(density_method)


  # create one full model-averaged ensemble
  averaged_posterior <- mix_posteriors(
    model_list   = model_list,
    parameters   = parameters,
    is_null_list = is_null_list,
    seed         = seed,
    n_samples    = n_samples,
    conditional  = FALSE
  )

  # prepare output object
  out <- list(
    conditional = list(),
    averaged    = list(),
    inference   = list()
  )

  for(i in seq_along(marginal_parameters)){

    if(all(is_null_list[[marginal_parameters[i]]])){
      warning(paste0("parameter '", marginal_parameters[i], "' does not contain any alternative hypothesis models."), immediate. = TRUE, call. = FALSE)
      next
    }

    # obtain model-averaged posterior conditional on including the parameter of interest
    # (different from individual conditionals)
    temp_conditional_posterior <- mix_posteriors(
      model_list   = model_list[!is_null_list[[marginal_parameters[i]]]],
      parameters   = parameters,
      is_null_list = lapply(is_null_list, function(l) l[!is_null_list[[marginal_parameters[i]]]]),
      seed         = seed,
      n_samples    = n_samples,
      conditional  = FALSE
    )

    # compute the marginals
    out[["averaged"]][[marginal_parameters[i]]] <- marginal_posterior(
      samples           = averaged_posterior,
      parameter         = marginal_parameters[i],
      formula           = formula,
      prior_samples     = TRUE,
      n_samples         = n_samples
    )
    out[["conditional"]][[marginal_parameters[i]]] <- marginal_posterior(
      samples           = temp_conditional_posterior,
      parameter         = marginal_parameters[i],
      formula           = formula,
      prior_samples     = TRUE,
      n_samples         = n_samples
    )

    # and inclusion Bayes factor
    out[["inference"]][[marginal_parameters[i]]] <- Savage_Dickey_BF(
      posterior            = out[["conditional"]][[marginal_parameters[i]]],
      null_hypothesis      = null_hypothesis,
      normal_approximation = normal_approximation,
      silent               = silent,
      density_method       = density_method
    )
  }

  attr(out, "null_hypothesis")      <- null_hypothesis
  attr(out, "normal_approximation") <- normal_approximation
  attr(out, "density_method")        <- density_method
  class(out) <- c(class(out), "marginal_inference")
  return(out)
}


#' @title Model-average marginal posterior distributions and
#' marginal Bayes factors based on BayesTools JAGS model via \code{marginal_inference}
#'
#' @description Creates marginal model-averaged and conditional
#' posterior distributions based on a BayesTools JAGS model, vector of parameters,
#' formula, and a list of conditional specifications for each parameter.
#' Computes inclusion Bayes factors for each marginal estimate via a Savage-Dickey
#' density approximation.
#'
#' @param marginal_parameters parameters for which the the marginal summary
#' should be created
#' @param conditional_list list of conditional parameters for each marginal parameter
#' @param parameters all parameters included in the model_list that are
#' relevant for the formula (all of which need to have specification of
#' \code{is_null_list})
#' @param compute_BF whether to compute inclusion Bayes factors. When
#' \code{FALSE}, the averaged and conditional marginal posteriors are returned
#' and the \code{inference} list is empty.
#' @inheritParams as_mixed_posteriors
#' @inheritParams marginal_inference
#' @inheritParams Savage_Dickey_BF
#'
#' @details For \code{as_marginal_inference()}, \code{conditional_list} is
#' applied separately to each output marginal or level. A requested conditional
#' parameter is active only for levels whose linear combination has a nonzero
#' weight for that parameter. Requested conditionals with zero weight are ignored
#' for that level; if no requested conditionals are active, the level uses the
#' fully averaged posterior and prior context. Level comparisons require the
#' compared levels to use the same effective conditional subset and rule.
#'
#' @return \code{as_marginal_inference} returns an object of class 'marginal_inference'.
#'
#' @seealso [marginal_inference] [as_mixed_posteriors]
#'
#' @export
as_marginal_inference <- function(model, marginal_parameters, parameters, conditional_list, conditional_rule, formula,
                                  null_hypothesis = 0, normal_approximation = FALSE,
                                  n_samples = 10000, silent = FALSE, force_plots = FALSE,
                                  density_method = "KDE", compute_BF = TRUE){

  # check input (majority of the checks performed within mix_posteriors)
  # check input
  if(!inherits(model, "BayesTools_fit"))
    stop("'model' must be a 'BayesTools_fit'")
  check_char(parameters, "parameters", check_length = FALSE)
  check_char(marginal_parameters, "marginal_parameters", check_length = FALSE)
  check_list(conditional_list, "conditional_list", check_length = length(marginal_parameters))
  check_char(conditional_rule, "conditional_rule")
  check_bool(compute_BF, "compute_BF")
  density_method <- .marginal_inference_density_method(density_method)

  priors <- attr(model, "prior_list")


  # create one full model-averaged ensemble
  averaged_posterior <- as_mixed_posteriors(
    model        = model,
    parameters   = parameters
  )

  # prepare output object
  out <- list(
    conditional = list(),
    averaged    = list(),
    inference   = list()
  )

  for(i in seq_along(marginal_parameters)){

    check_char(
      conditional_list[[marginal_parameters[i]]],
      sprintf("conditional_list[[%1$s]]", marginal_parameters[i]),
      check_length = FALSE,
      allow_values = c(parameters, "PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null"),
      allow_NULL = TRUE
    )

    # compute the marginals
    out[["averaged"]][[marginal_parameters[i]]] <- marginal_posterior(
      samples           = averaged_posterior,
      parameter         = marginal_parameters[i],
      formula           = formula,
      prior_samples     = TRUE,
      n_samples         = n_samples
    )

    out[["conditional"]][[marginal_parameters[i]]] <- .marginal_inference_conditional_posterior(
      model              = model,
      parameters         = parameters,
      marginal_parameter = marginal_parameters[i],
      formula            = formula,
      averaged_marginal  = out[["averaged"]][[marginal_parameters[i]]],
      prior_list         = attr(averaged_posterior, "prior_list"),
      conditional        = conditional_list[[marginal_parameters[i]]],
      conditional_rule   = conditional_rule,
      n_samples          = n_samples,
      force_plots        = force_plots
    )

    if(length(out[["conditional"]][[marginal_parameters[i]]]) == 0){
      out[["averaged"]][[marginal_parameters[i]]] <- NULL
      next
    }

    if(compute_BF){
      # and inclusion Bayes factor
      out[["inference"]][[marginal_parameters[i]]] <- Savage_Dickey_BF(
        posterior            = out[["conditional"]][[marginal_parameters[i]]],
        null_hypothesis      = null_hypothesis,
        normal_approximation = normal_approximation,
        silent               = silent,
        density_method       = density_method
      )
    }
  }

  attr(out, "null_hypothesis")      <- null_hypothesis
  attr(out, "normal_approximation") <- normal_approximation
  attr(out, "density_method")        <- density_method
  class(out) <- c(class(out), "marginal_inference")
  return(out)
}

.marginal_inference_density_method <- function(density_method){

  density_method <- .posterior_density_method(density_method)
  if(identical(density_method, "precomputed")){
    stop(
      "'density_method = \"precomputed\"' is not supported by ",
      "'marginal_inference()' or 'as_marginal_inference()'. ",
      "Precomputed densities must be computed for the marginal posterior ",
      "itself and passed directly to 'Savage_Dickey_BF()'.",
      call. = FALSE
    )
  }

  return(density_method)
}

.marginal_inference_condition_key <- function(conditional, conditional_rule = "AND"){

  .condition_event_key(conditional, conditional_rule)
}

.marginal_inference_level_conditionals <- function(marginal, prior_list, conditional,
                                                   conditional_rule = "AND"){

  scalar_marginal <- !is.list(marginal)
  levels <- if(scalar_marginal) ".scalar" else names(marginal)
  conditionals <- lapply(levels, function(level){
    level_marginal <- if(scalar_marginal) marginal else marginal[[level]]
    weights <- attr(level_marginal, "linear_weights")
    if(is.null(weights)){
      return(conditional)
    }
    if(!is.null(dim(weights))){
      row_conditionals <- lapply(seq_len(nrow(weights)), function(row_i){
        .prior_linear_active_conditionals(
          prior_list  = prior_list,
          weights     = weights[row_i, ],
          conditional = conditional
        )
      })
      row_keys <- vapply(row_conditionals, .condition_labels_key, character(1))
      if(length(unique(row_keys)) > 1L){
        stop(
          "Row-varying active conditional sets are not supported for marginal inference.",
          call. = FALSE
        )
      }
      return(row_conditionals[[1]])
    }
    .prior_linear_active_conditionals(
      prior_list   = prior_list,
      weights      = weights,
      conditional  = conditional
    )
  })
  names(conditionals) <- levels

  conditionals
}

.marginal_inference_conditional_posterior <- function(model, parameters, marginal_parameter,
                                                      formula, averaged_marginal, prior_list,
                                                      conditional, conditional_rule, n_samples,
                                                      force_plots){

  scalar_marginal <- !is.list(averaged_marginal)
  if(scalar_marginal && length(conditional) == 0L){
    return(averaged_marginal)
  }

  level_conditionals <- .marginal_inference_level_conditionals(
    marginal    = averaged_marginal,
    prior_list  = prior_list,
    conditional = conditional,
    conditional_rule = conditional_rule
  )

  conditional_marginal <- averaged_marginal
  marginal_cache <- list()

  for(level in names(level_conditionals)){
    level_conditional <- level_conditionals[[level]]
    level_event <- .condition_event(
      prior_list        = prior_list,
      conditional       = level_conditional,
      conditional_rule  = conditional_rule
    )
    key <- level_event[["condition_key"]]

    if(scalar_marginal && length(level_conditional) == 0L){
      conditional_marginal <- .condition_event_set_attributes(
        averaged_marginal,
        level_event,
        effective = TRUE
      )
      next
    }

    if(is.null(marginal_cache[[key]])){
      conditional_posterior <- as_mixed_posteriors(
        model            = model,
        parameters       = parameters,
        conditional      = if(length(level_conditional) == 0) NULL else level_conditional,
        conditional_rule = conditional_rule,
        force_plots      = force_plots
      )

      if(length(conditional_posterior) == 0){
        return(list())
      }else{
        marginal_cache[[key]] <- marginal_posterior(
          samples       = conditional_posterior,
          parameter     = marginal_parameter,
          formula       = formula,
          prior_samples = TRUE,
          n_samples     = n_samples
        )
      }
    }

    if(length(marginal_cache[[key]]) == 0){
      return(list())
    }

    if(scalar_marginal){
      conditional_marginal <- .condition_event_set_attributes(
        marginal_cache[[key]],
        level_event,
        effective = TRUE
      )
    }else{
      conditional_marginal[[level]] <- marginal_cache[[key]][[level]]
      conditional_marginal[[level]] <- .condition_event_set_attributes(
        conditional_marginal[[level]],
        level_event,
        effective = TRUE
      )
    }
  }

  conditional_marginal
}
