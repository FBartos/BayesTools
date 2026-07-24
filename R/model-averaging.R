#' @title Compute posterior probabilities and inclusion Bayes factors
#'
#' @description Computes prior probabilities, posterior probabilities,
#' and inclusion Bayes factors based either on (1) a list of models,
#' vector of parameters, and a list of indicators the models represent the
#' null or alternative hypothesis for each parameter, (2) on prior
#' model odds, marginal likelihoods, and indicator whether
#' the models represent the null or alternative hypothesis, or (3) list of
#' models for each model.
#'
#' @param model_list list of models, each of which contains marginal
#' likelihood estimated with bridge sampling \code{marglik} and prior model
#' odds \code{prior_weights}
#' @param parameters vector of parameters names for which inference should
#' be drawn
#' @param is_null_list list with entries for each parameter carrying either
#' logical vector of indicators specifying whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis; use \code{0} or \code{integer(0)}
#' when no models are null)
#' @param prior_weights vector of prior model odds
#' @param margliks vector of marginal likelihoods
#' @param is_null logical vector of indicators specifying whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis; use \code{0} or \code{integer(0)}
#' when no models are null)
#' @param conditional whether prior and posterior model probabilities should
#' be returned only for the conditional model. Defaults to \code{FALSE}
#'
#'
#' @return \code{compute_inference} returns a named list of prior probabilities,
#' posterior probabilities, and Bayes factors, \code{ppoint} gives the
#' distribution function, \code{ensemble_inference} gives a list of named lists of
#' inferences for each parameter, and \code{models_inference} returns a list of
#' models, each expanded by the inference list.
#'
#' @seealso [mix_posteriors] [BayesTools_ensemble_tables]
#'
#' @name ensemble_inference
#' @export compute_inference
#' @export ensemble_inference
#' @export models_inference
NULL

#' @rdname ensemble_inference
compute_inference <- function(prior_weights, margliks, is_null = NULL, conditional = FALSE){

  check_real(prior_weights, "prior_weights", lower = 0, check_length = 0)
  check_real(margliks,   "margliks", check_length = length(prior_weights))
  check_bool(conditional, "conditional", allow_NA = FALSE)
  is_null <- .model_averaging_is_null(is_null, length(prior_weights))

  prior_probs <- .model_averaging_prior_probs(prior_weights)
  margliks    <- .model_averaging_margliks(margliks, prior_probs)
  post_probs  <- .model_averaging_post_probs(margliks, prior_probs)
  BF          <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)

  if(conditional){
    if(all(is_null))
      stop("Conditional inference requires at least one non-null model.", call. = FALSE)
    prior_probs <- .model_averaging_prior_probs(ifelse(is_null, 0, prior_weights))
    margliks    <- .model_averaging_margliks(margliks, prior_probs)
    post_probs  <- .model_averaging_post_probs(margliks, prior_probs)
  }

  output <- list(
    prior_probs = prior_probs,
    post_probs  = post_probs,
    BF          = BF
  )

  attr(output, "is_null")     <- is_null
  attr(output, "conditional") <- conditional
  class(output) <- c(class(output), "inference")

  return(output)
}

.model_averaging_prior_probs <- function(prior_weights){

  if(any(!is.finite(prior_weights))){
    stop("'prior_weights' must be finite.", call. = FALSE)
  }
  if(!any(prior_weights > 0)){
    stop("At least one prior model weight must be positive.", call. = FALSE)
  }

  scaled_weights <- prior_weights / max(prior_weights)
  scaled_weights / sum(scaled_weights)
}

.model_averaging_margliks <- function(margliks, prior_probs){

  if(any(is.infinite(margliks) & margliks > 0 & prior_probs > 0, na.rm = TRUE)){
    stop("Infinite positive marginal likelihoods are not supported.", call. = FALSE)
  }

  margliks[is.na(margliks)] <- -Inf

  if(!any(is.finite(margliks) & prior_probs > 0)){
    stop("No finite marginal likelihoods are available for models with positive prior probability.", call. = FALSE)
  }

  margliks
}

.model_averaging_post_probs <- function(margliks, prior_probs){

  unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
}

.posterior_mixture_sample_counts <- function(post_probs, n_samples){

  .mixture_sample_counts(post_probs, n_samples)
}

.mixture_sample_counts <- function(probs, n_samples, preserve_positive = TRUE){

  probs <- as.numeric(probs)
  if(length(probs) == 0L){
    return(integer())
  }

  if(length(n_samples) != 1L || !is.finite(n_samples) || n_samples < 0 || !.is.wholenumber(n_samples)){
    stop("'n_samples' must be a non-negative integer.", call. = FALSE)
  }
  n_samples <- as.integer(n_samples)

  if(any(!is.finite(probs)) || any(probs < 0)){
    stop("'probs' must contain non-negative finite values.", call. = FALSE)
  }
  if(sum(probs) <= 0){
    stop("At least one mixture probability must be positive.", call. = FALSE)
  }

  probs <- probs / sum(probs)
  raw_counts <- probs * n_samples
  counts <- floor(raw_counts)

  remaining <- n_samples - sum(counts)
  if(remaining > 0L){
    add_order <- order(-(raw_counts - counts), -probs, seq_along(probs))
    counts[add_order[seq_len(remaining)]] <- counts[add_order[seq_len(remaining)]] + 1L
  }

  if(preserve_positive && n_samples > 0L){
    positive <- probs > 0
    missing <- which(positive & counts == 0L)
    if(length(missing) > 0L && sum(positive) <= n_samples){
      for(m in missing){
        donors <- which(counts > 1L)
        if(length(donors) == 0L){
          break
        }

        donor_cost <- abs((counts[donors] - 1L) - raw_counts[donors]) - abs(counts[donors] - raw_counts[donors])
        donor <- donors[order(donor_cost, -counts[donors], -probs[donors], donors)[1L]]
        counts[donor] <- counts[donor] - 1L
        counts[m] <- 1L
      }
    }
  }

  as.integer(counts)
}

#' @rdname ensemble_inference
ensemble_inference <- function(model_list, parameters, is_null_list, conditional = FALSE){

  # check input
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("marglik", "prior_weights"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m)inherits(m[["marglik"]], what = "bridge"))))
    stop("model_list:marglik must contain 'bridgesampling' marginal likelihoods")
  sapply(model_list, function(m)check_real(m[["prior_weights"]], "model_list:prior_weights", lower = 0))


  # extract the object
  margliks      <- sapply(model_list, function(m) m[["marglik"]][["logml"]])
  prior_weights <- sapply(model_list, function(m) m[["prior_weights"]])

  out <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    out[[parameters[p]]] <- compute_inference(prior_weights = prior_weights, margliks = margliks, is_null = is_null_list[[p]], conditional = conditional)

    # add parameter names
    parameter_name    <- parameters[p]
    formula_parameter <- unique(unlist(lapply(model_list, function(m) attr(attr(m[["fit"]], "prior_list")[[parameters[p]]], "parameter"))))

    if(!is.null(unlist(formula_parameter))){
      parameter_name <- format_parameter_names(parameter_name, formula_parameters = formula_parameter, formula_prefix = TRUE)
      class(out[[parameters[p]]]) <- c(class(out[[parameters[p]]]), "inference.formula")
      attr(out[[parameters[p]]], "formula_parameter")  <- formula_parameter
    }
    attr(out[[parameters[p]]], "parameter_name")  <- parameter_name

  }

  attr(out, "conditional") <- conditional
  return(out)
}

#' @rdname ensemble_inference
models_inference   <- function(model_list){

  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("marglik", "prior_weights"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m)inherits(m[["marglik"]], what = "bridge"))))
    stop("model_list:marglik must contain 'bridgesampling' marginal likelihoods")
  sapply(model_list, function(m)check_real(m[["prior_weights"]], "model_list:prior_weights", lower = 0))

  margliks    <- sapply(model_list, function(model)model[["marglik"]][["logml"]])
  prior_weights  <- sapply(model_list, function(model)model[["prior_weights"]])
  prior_probs <- .model_averaging_prior_probs(prior_weights)
  margliks    <- .model_averaging_margliks(margliks, prior_probs)
  post_probs  <- .model_averaging_post_probs(margliks, prior_probs)
  incl_BF     <- sapply(seq_along(model_list), function(i){
    is_null <- rep(TRUE, length(model_list))
    is_null[i] <- FALSE
    return(inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null))
  })

  for(i in seq_along(model_list)){
    model_list[[i]][["inference"]] <- list(
      "m_number"     = i,
      "marglik"      = margliks[i],
      "prior_prob"   = prior_probs[i],
      "post_prob"    = post_probs[i],
      "inclusion_BF" = incl_BF[i]
    )
  }

  return(model_list)
}
