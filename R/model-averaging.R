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
#' corresponding to the null hypothesis)
#' @param prior_weights vector of prior model odds
#' @param margliks vector of marginal likelihoods
#' @param is_null logical vector of indicators specifying whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis)
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
  check_bool(conditional, "conditional")
  if(!is.null(is_null)){
    if(is.numeric(is_null)){
      check_int(is_null, "is_null", lower = 0, upper = length(prior_weights), check_length = 0)
      is_null <- c(1:length(prior_weights)) %in% is_null
    }else if(is.logical(is_null)){
      check_bool(is_null, "is_null", check_length = length(prior_weights))
    }else{
      stop("'is_null' argument must be either logical vector, integer vector, or NULL.")
    }
  }else{
    is_null <- rep(FALSE, length(prior_weights))
  }

  prior_probs <- prior_weights / sum(prior_weights)
  post_probs  <- unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
  BF          <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)

  if(conditional){
    prior_probs <- ifelse(is_null, 0, prior_weights) / sum(prior_weights[!is_null])
    post_probs  <- unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
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
  margliks    <- ifelse(is.na(margliks), -Inf, margliks)
  prior_weights  <- sapply(model_list, function(model)model[["prior_weights"]])
  prior_probs <- prior_weights / sum(prior_weights)
  post_probs  <- unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
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

#' @title Model-average posterior distributions
#'
#' @description Model-averages posterior distributions based on
#' a list of models, vector of parameters, and a list of
#' indicators of the null or alternative hypothesis models
#' for each parameter.
#'
#' @param seed integer specifying seed for sampling posteriors for
#' model averaging. Defaults to \code{NULL}.
#' @param n_samples number of samples to be drawn for the model-averaged
#' posterior distribution
#' @inheritParams ensemble_inference
#'
#' @return \code{mix_posteriors} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [ensemble_inference] [BayesTools_ensemble_tables] [as_mixed_posteriors]
#'
#' @name mix_posteriors
#' @export
mix_posteriors <- function(model_list, parameters, is_null_list, conditional = FALSE, seed = NULL, n_samples = 10000){

  # check input
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("fit", "marglik", "prior_weights"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m) inherits(m[["fit"]], what = "runjags")) | sapply(model_list, function(m)inherits(m[["fit"]], what = "stanfit")) | sapply(model_list, function(m)inherits(m[["fit"]], what = "null_model"))))
    stop("model_list:fit must contain 'runjags' or 'rstan' models")
  if(!all(sapply(model_list, function(m) inherits(m[["marglik"]], what = "bridge"))))
    stop("model_list:marglik must contain 'bridgesampling' marginal likelihoods")
  if(!all(unlist(sapply(model_list, function(m) sapply(attr(m[["fit"]], "prior_list"), function(p) is.prior(p))))))
    stop("model_list:priors must contain 'BayesTools' priors")
  sapply(model_list, function(m) check_real(m[["prior_weights"]], "model_list:prior_weights", lower = 0))


  # extract the object
  fits           <- lapply(model_list, function(m) m[["fit"]])
  margliks       <- sapply(model_list, function(m) m[["marglik"]][["logml"]])
  priors         <- lapply(model_list, function(m) attr(m[["fit"]], "prior_list"))
  formula_priors <- lapply(model_list, function(m) m[["formula_priors"]])
  prior_weights  <- sapply(model_list, function(m) m[["prior_weights"]])

  inference  <- ensemble_inference(model_list, parameters, is_null_list, conditional)
  out <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter    <- parameters[p]
    temp_inference    <- inference[[temp_parameter]]
    temp_priors       <- lapply(priors, function(p) p[[temp_parameter]])

    if(any(sapply(temp_priors, is.prior.weightfunction)) && all(sapply(temp_priors, is.prior.weightfunction) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.prior.none) | sapply(temp_priors, is.null))){
      # weightfunctions:

      # replace missing priors with default prior: none
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior_none()
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        temp_priors[[i]][["prior_weights"]] <- temp_inference$prior_probs[i]
      }

      out[[temp_parameter]] <- .mix_posteriors.weightfunction(fits, temp_priors, temp_parameter, temp_inference$post_probs, seed, n_samples)

    }else if(any(sapply(temp_priors, is.prior.factor)) && all(sapply(temp_priors, is.prior.factor) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # factor priors

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        temp_priors[[i]][["prior_weights"]] <- temp_inference$prior_probs[i]
      }

      out[[temp_parameter]] <- .mix_posteriors.factor(fits, temp_priors, temp_parameter, temp_inference$post_probs, seed, n_samples)

    }else if(any(sapply(temp_priors, is.prior.vector)) && all(sapply(temp_priors, is.prior.vector) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # vector priors:

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        temp_priors[[i]][["prior_weights"]] <- temp_inference$prior_probs[i]
      }

      out[[temp_parameter]] <- .mix_posteriors.vector(fits, temp_priors, temp_parameter, temp_inference$post_probs, seed, n_samples)

    }else if(all(sapply(temp_priors, is.prior.simple) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # simple priors:

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        temp_priors[[i]][["prior_weights"]] <- temp_inference$prior_probs[i]
      }

      out[[temp_parameter]] <- .mix_posteriors.simple(fits, temp_priors, temp_parameter, temp_inference$post_probs, seed, n_samples)

    }else{
      stop("The posterior samples cannot be mixed: unsupported mixture of prior distributions.")
    }

    # add formula relevant information
    if(!is.null(unique(unlist(lapply(temp_priors, attr, which = "parameter"))))){
      class(out[[temp_parameter]]) <- c(class(out[[temp_parameter]]), "mixed_posteriors.formula")
      attr(out[[temp_parameter]], "formula_parameter")  <- unique(unlist(lapply(temp_priors, attr, which = "parameter")))
    }

  }

  class(out) <- c(class(out), "mixed_posteriors")
  return(out)
}

.mix_posteriors.simple         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.simple) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of simple priors")

  # gather and check compatibility of prior distributions
  priors_info <- lapply(priors, function(p){
    if(is.prior.point(p) | is.prior.none(p)){
      return(FALSE)
    }else{
      list(
        "interaction"       = .is_prior_interaction(p),
        "interaction_terms" = attr(p, "interaction_terms")
      )
    }
  })
  priors_info <- priors_info[!sapply(priors_info, isFALSE)]
  if(length(priors_info) >= 2 && any(!unlist(lapply(priors_info, function(i) all.equal(i, priors_info[[1]]))))){
    stop("non-matching prior factor type specifications")
  }else if(length(priors_info) != 0){
    priors_info <- priors_info[[1]]
  }



  # set seed at the beginning makes sure that the samples of different parameters from the same models retain their correlation
  if(!is.null(seed)){
    set.seed(seed)
  }else{
    set.seed(1)
  }

  # prepare output objects
  samples <- NULL
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  for(i in seq_along(fits)[round(post_probs * n_samples) > 1]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- suppressWarnings(coda::as.mcmc(fits[[i]]))
      if(!is.matrix(model_samples)){
        # deal with automatic coercion into a vector in case of a single predictor
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- fits[[i]]$monitor
      }
    }else if(inherits(fits[[i]], "stanfit")){
      .check_rstan()
      model_samples <- .extract_stan(fits[[i]])
    }

    # sample indexes
    temp_ind <- sample(nrow(model_samples), round(n_samples * post_probs[i]), replace = TRUE)

    if(is.prior.point(priors[[i]])){
      # not sampling the priors as the samples would be already transformed
      samples <- c(samples, rep(priors[[i]]$parameters[["location"]], length(temp_ind)))
    }else{
      samples <- c(samples, model_samples[temp_ind, parameter])
    }

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  samples <- unname(samples)
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  attr(samples, "interaction")       <- if(length(priors_info) == 0) FALSE else priors_info[["interaction"]]
  attr(samples, "interaction_terms") <- priors_info[["interaction_terms"]]
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.simple")

  return(samples)
}
.mix_posteriors.vector         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.vector) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of vector priors")

  # set seed at the beginning makes sure that the samples of different parameters from the same models retain their correlation
  if(!is.null(seed)){
    set.seed(seed)
  }else{
    set.seed(1)
  }

  # prepare output objects
  K <- unique(sapply(priors[sapply(priors, is.prior.vector)], function(p) p$parameters[["K"]]))
  if(length(K) != 1)
    stop("all vector priors must be of the same length")

  samples    <- matrix(nrow = 0, ncol = K)
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  for(i in seq_along(fits)[round(post_probs * n_samples) > 1]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- suppressWarnings(coda::as.mcmc(fits[[i]]))
      if(!is.matrix(model_samples)){
        # deal with automatic coercion into a vector in case of a single predictor
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- fits[[i]]$monitor
      }
    }else if(inherits(fits[[i]], "stanfit")){
      .check_rstan()
      model_samples <- .extract_stan(fits[[i]])
    }

    # sample indexes
    temp_ind <- sample(nrow(model_samples), round(n_samples * post_probs[i]), replace = TRUE)

    if(is.prior.point(priors[[i]])){
      # not sampling the priors as the samples would be already transformed
      samples <- rbind(samples, matrix(priors[[i]]$parameters[["location"]], nrow = length(temp_ind), ncol = K))
    }else if(K == 1){
      samples <- rbind(samples, matrix(model_samples[temp_ind, parameter], nrow = length(temp_ind), ncol = K))
    }else{
      samples <- rbind(samples, model_samples[temp_ind, paste0(parameter,"[",1:K,"]")])
    }

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  rownames(samples) <- NULL
  colnames(samples) <- paste0(parameter,"[",1:K,"]")
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.vector")

  return(samples)
}
.mix_posteriors.factor         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.factor) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of factor priors")

  # check the prior levels
  levels <- unique(sapply(priors[sapply(priors, is.prior.factor)], .get_prior_factor_levels))
  if(length(levels) != 1)
    stop("all factor priors must be of the same number of levels")

  # gather and check compatibility of prior distributions
  priors_info <- lapply(priors, function(p){
    if(is.prior.point(p) | is.prior.none(p)){
      return(FALSE)
    }else if(is.prior.factor(p)){
      return(list(
        "levels"            = .get_prior_factor_levels(p),
        "level_names"       = .get_prior_factor_level_names(p),
        "interaction"       = .is_prior_interaction(p),
        "interaction_terms" = attr(p, "interaction_terms"),
        "treatment"         = is.prior.treatment(p),
        "independent"       = is.prior.independent(p),
        "orthonormal"       = is.prior.orthonormal(p),
        "meandif"           = is.prior.meandif(p)
      ))
    }else{
      stop("unsupported prior type")
    }
  })
  priors_info <- priors_info[!sapply(priors_info, isFALSE)]
  if(length(priors_info) >= 2 && any(!unlist(lapply(priors_info, function(i) all.equal(i, priors_info[[1]]))))){
    stop("non-matching prior factor type specifications")
  }
  if(length(priors_info) != 0){
    priors_info <- priors_info[[1]]
  }


  if(priors_info[["treatment"]]){

    if(levels == 1){

      samples <- .mix_posteriors.simple(fits, priors, parameter, post_probs, seed, n_samples)

      sample_ind <- attr(samples, "sample_ind")
      models_ind <- attr(samples, "models_ind")

      samples <- matrix(samples, ncol = 1)

    }else{

      # keep the same seed across levels
      if(is.null(seed)){
        seed <- sample(666666, 1)
      }

      samples <- lapply(1:levels, function(i) .mix_posteriors.simple(fits, priors, paste0(parameter, "[", i, "]"), post_probs, seed, n_samples))

      sample_ind <- attr(samples[[1]], "sample_ind")
      models_ind <- attr(samples[[1]], "models_ind")

      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- paste0(parameter,"[",priors_info$level_names[-1],"]")
    attr(samples, "sample_ind") <- sample_ind
    attr(samples, "models_ind") <- models_ind
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- priors
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(priors_info[["independent"]]){

    if(levels == 1){

      samples <- .mix_posteriors.simple(fits, priors, parameter, post_probs, seed, n_samples)

      sample_ind <- attr(samples, "sample_ind")
      models_ind <- attr(samples, "models_ind")

      samples <- matrix(samples, ncol = 1)

    }else{

      # keep the same seed across levels
      if(is.null(seed)){
        seed <- sample(666666, 1)
      }

      samples <- lapply(1:levels, function(i) .mix_posteriors.simple(fits, priors, paste0(parameter, "[", i, "]"), post_probs, seed, n_samples))

      sample_ind <- attr(samples[[1]], "sample_ind")
      models_ind <- attr(samples[[1]], "models_ind")

      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- paste0(parameter,"[",priors_info$level_names,"]")
    attr(samples, "sample_ind") <- sample_ind
    attr(samples, "models_ind") <- models_ind
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- priors
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(priors_info[["orthonormal"]] | priors_info[["meandif"]]){

    for(i in seq_along(priors)){
      if(is.prior.factor(priors[[i]])){
        priors[[i]]$parameters[["K"]] <- levels
      }
    }

    samples <- .mix_posteriors.vector(fits, priors, parameter, post_probs, seed, n_samples)
    class(samples) <- c(class(samples), "mixed_posteriors.factor")

  }

  attr(samples, "levels")            <- priors_info[["levels"]]
  attr(samples, "level_names")       <- priors_info[["level_names"]]
  attr(samples, "interaction")       <- if(length(priors_info) == 0) FALSE else priors_info[["interaction"]]
  attr(samples, "interaction_terms") <- priors_info[["interaction_terms"]]
  attr(samples, "treatment")         <- priors_info[["treatment"]]
  attr(samples, "independent")       <- priors_info[["independent"]]
  attr(samples, "orthonormal")       <- priors_info[["orthonormal"]]
  attr(samples, "meandif")           <- priors_info[["meandif"]]

  return(samples)
}
.mix_posteriors.weightfunction <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.weightfunction) | sapply(priors, is.prior.point) | sapply(priors, is.prior.none)))
    stop("'priors' must be a list of weightfunction priors distributions")


  # set seed at the beginning makes sure that the samples of different parameters from the same models retain their correlation
  if(!is.null(seed)){
    set.seed(seed)
  }else{
    set.seed(1)
  }

  # obtain mapping for the weight coefficients
  omega_mapping <- weightfunctions_mapping(priors)
  omega_cuts    <- weightfunctions_mapping(priors, cuts_only = TRUE)
  omega_names   <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))

  # prepare output objects
  samples    <- matrix(nrow = 0, ncol = length(omega_cuts) - 1)
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  for(i in seq_along(fits)[round(post_probs * n_samples) > 1]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- suppressWarnings(coda::as.mcmc(fits[[i]]))
      if(!is.matrix(model_samples)){
        # deal with automatic coercion into a vector in case of a single predictor
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- fits[[i]]$monitor
      }
    }else if(inherits(fits[[i]], "stanfit")){
      .check_rstan()
      model_samples <- .extract_stan(fits[[i]])
    }

    # sample indexes
    temp_ind <- sample(nrow(model_samples), round(n_samples * post_probs[i]), replace = TRUE)

    if(is.prior.none(priors[[i]])){
      samples <- rbind(samples, matrix(1, ncol = length(omega_cuts) - 1, nrow = length(temp_ind)))
    }else{
      samples <- rbind(samples, model_samples[temp_ind, paste0("omega[",omega_mapping[[i]],"]")])
    }

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}

#' @title Export BayesTools JAGS model posterior distribution as model-average posterior distributions via \code{mix_posteriors}
#'
#' @description Creates a model-averages posterior distributions on a single
#' model that allows mimicking the [mix_posteriors] functionality. This function
#' is useful when the model-averaged ensemble is based on [prior_spike_and_slab]
#' or [prior_mixture] priors - the model-averaging is done within the model.
#'
#' @param model model fit via the [JAGS_fit] function
#' @param conditional a character vector of parameters to be conditioned on
#' @param conditional_rule a character string specifying the rule for conditioning.
#' Either "AND" or "OR". Defaults to "AND".
#' @param force_plots temporal argument allowing to generate conditional posterior samples
#' suitable for prior and posterior plots. Only available when conditioning on a
#' single parameter.
#' @inheritParams ensemble_inference
#' @inheritParams mix_posteriors
#'
#' @return \code{as_mix_posteriors} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors]
#'
#' @name as_mixed_posteriors
#' @export
as_mixed_posteriors <- function(model, parameters, conditional = NULL, conditional_rule = "AND", force_plots = FALSE){

  # check input
  if(!inherits(model, "BayesTools_fit"))
    stop("'model' must be a 'BayesTools_fit'")
  check_char(parameters, "parameters", check_length = FALSE)
  check_char(conditional, "conditional", check_length = FALSE, allow_values = c(parameters, "PET", "PEESE", "PETPEESE", "omega"), allow_NULL = TRUE)
  check_char(conditional_rule, "conditional_rule", allow_values = c("AND", "OR"))

  # extract the list of priors
  priors <- attr(model, "prior_list")

  # extract the samples
  model_samples <- suppressWarnings(coda::as.mcmc(model))
  if(!is.matrix(model_samples)){
    # deal with automatic coercion into a vector in case of a single predictor
    model_samples <- matrix(model_samples, ncol = 1)
    colnames(model_samples) <- model$monitor
  }

  # apply conditioning
  if(length(conditional) > 0){

    # subset the posterior distribution
    conditioning_samples <- do.call(cbind, lapply(conditional, function(parameter){

      # special cases for PET / PEESE / PET-PEESE / weightfunctions
      if(parameter == "PET" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_PET <- sapply(priors[["bias"]], is.prior.PET)
        return(model_samples[, "bias_indicator"] %in% which(is_PET))
      }
      if(parameter == "PEESE" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
        return(model_samples[, "bias_indicator"] %in% which(is_PEESE))
      }
      if(parameter == "PETPEESE" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_PET   <- sapply(priors[["bias"]], is.prior.PET)
        is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
        return(model_samples[, "bias_indicator"] %in% which(is_PET | is_PEESE))
      }
      if(parameter == "omega" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_weightfunction <- sapply(priors[["bias"]], is.prior.weightfunction)
        return(model_samples[, "bias_indicator"] %in% which(is_weightfunction))
      }

      # normal cases
      temp_prior <- priors[[parameter]]

      if(is.prior.spike_and_slab(temp_prior)){

        return(model_samples[,paste0(parameter, "_indicator")] == 1)

      }else if(is.prior.mixture(temp_prior)){

        components <- attr(temp_prior, "components")
        if(!all(components %in% c("null", "alternative")))
          stop("conditional mixture posterior distributions are available only for 'null' and 'alternative' components")

        return(model_samples[,paste0(parameter, "_indicator")] %in% which(components == "alternative"))

      }else{

        warning(sprintf("The parameter '%s' is not a conditional parameter.", parameter), call. = FALSE, immediate. = TRUE)
        return(rep(TRUE, nrow(model_samples)))

      }
    }))
    conditioning_samples <- apply(conditioning_samples, 1, ifelse(conditional_rule == "AND", all, any))

    if(sum(conditioning_samples) == 0){
      warning("No samples left after conditioning.", call. = FALSE, immediate. = TRUE)
      return(list())
    }


    model_samples <- model_samples[conditioning_samples,,drop=FALSE]

    # set prior weights to 0 for null distributions
    # TODO: this needs to be implemented for enabling of the conditional mixture posterior distributions when more than one components is present
    # (e.g., conditional marginal and posterior plots)
    # the current workaround is suitable only for a single parameters (to produce averaged prior and posterior plots)
    if(length(conditional) == 1 && length(parameters) == 1 && conditional == parameters && force_plots){

      # special cases for PET / PEESE / PET-PEESE / weightfunctions
      if(conditional == "PET" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_PET <- sapply(priors[["bias"]], is.prior.PET)
        for(i in seq(along = is_PET)){
          if(!is_PET[i]){
            priors[[parameters]][[i]][["prior_weights"]] <- 0
          }
        }
      }else if(conditional == "PEESE" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        for(i in seq(along = is_PEESE)){
          if(!is_PEESE[i]){
            priors[[parameters]][[i]][["prior_weights"]] <- 0
          }
        }
      }else if(conditional == "PETPEESE" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_PET   <- sapply(priors[["bias"]], is.prior.PET)
        is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
        for(i in seq(along = is_PET)){
          if(!(is_PET[i] || is_PEESE[i])){
            priors[[parameters]][[i]][["prior_weights"]] <- 0
          }
        }
      }else if(conditional == "omega" && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])){
        is_weightfunction <- sapply(priors[["bias"]], is.prior.weightfunction)
        for(i in seq(along = is_weightfunction)){
          if(!is_weightfunction[i]){
            priors[[parameters]][[i]][["prior_weights"]] <- 0
          }
        }
      }else if(is.prior.spike_and_slab(priors[[parameters]])){
        priors[[parameters]][["inclusion"]] <- prior("spike", list(1))
      }else if(is.prior.mixture(priors[[parameters]])){

        components <- attr(priors[[parameters]], "components")

        attr(priors[[parameters]], "prior_weights")[which(components == "null")] <- 0
        for(i in seq(along = components)){
          if(components[i] == "null"){
            priors[[parameters]][[i]][["prior_weights"]] <- 0
          }
        }
      }
    }
  }


  out    <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter <- parameters[p]
    temp_prior     <- priors[[temp_parameter]]

    if(is.prior.spike_and_slab(temp_prior)){
      # spike and slab priors
      out[[temp_parameter]] <- .as_mixed_posteriors.spike_and_slab(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.mixture(temp_prior)){
      # mixture priors
      out[[temp_parameter]] <- .as_mixed_posteriors.mixture(model_samples, temp_prior, temp_parameter, conditional)

    }else if(is.prior.weightfunction(temp_prior)){
      # weight functions
      out[[temp_parameter]] <- .as_mixed_posteriors.weightfunction(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.factor(temp_prior)){
      # factor priors
      out[[temp_parameter]] <- .as_mixed_posteriors.factor(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.vector(temp_prior)){
      # vector priors
      out[[temp_parameter]] <- .as_mixed_posteriors.vector(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.simple(temp_prior)){
      # simple priors
      out[[temp_parameter]] <- .as_mixed_posteriors.simple(model_samples, temp_prior, temp_parameter)

    }else{
      stop("The posterior samples cannot be mixed: unsupported prior distributions.")
    }

    # add formula relevant information
    if(!is.null(attr(temp_prior, which = "parameter"))){
      class(out[[temp_parameter]]) <- c(class(out[[temp_parameter]]), "mixed_posteriors.formula")
      attr(out[[temp_parameter]], "formula_parameter")  <- attr(temp_prior, which = "parameter")
    }

    # add conditioning information
    attr(out[[temp_parameter]], "conditional")      <- conditional
    attr(out[[temp_parameter]], "conditional_rule") <- conditional_rule

  }

  attr(out, "prior_list")       <- priors
  attr(out, "conditional")      <- conditional
  attr(out, "conditional_rule") <- conditional_rule
  class(out) <- c(class(out), "as_mixed_posteriors", "mixed_posteriors")
  return(out)
}

.as_mixed_posteriors.simple         <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  # gather information about the prior distribution
  prior_info <- list(
    "interaction"       = .is_prior_interaction(prior),
    "interaction_terms" = attr(prior, "interaction_terms")
  )

  # prepare output objects
  samples <- model_samples[, parameter]

  # format the output
  samples <- unname(samples)
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  attr(samples, "interaction")       <- if(length(prior_info) == 0) FALSE else prior_info[["interaction"]]
  attr(samples, "interaction_terms") <- prior_info[["interaction_terms"]]
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.simple")

  return(samples)
}
.as_mixed_posteriors.vector         <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  # gather information about the prior distribution
  K <- prior$parameter[["K"]]
  if(length(K) != 1)
    stop("all vector prior must be of the same length")

  # prepare output objects
  if(K == 1){
    samples <- model_samples[, parameter, drop = FALSE]
  }else{
    samples <- model_samples[, paste0(parameter,"[",1:K,"]"), drop = FALSE]
  }

  rownames(samples) <- NULL
  colnames(samples) <- paste0(parameter,"[",1:K,"]")
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.vector")

  return(samples)
}
.as_mixed_posteriors.factor         <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  # gather information about the prior distribution
  prior_info <- list(
    "levels"            = .get_prior_factor_levels(prior),
    "level_names"       = .get_prior_factor_level_names(prior),
    "interaction"       = .is_prior_interaction(prior),
    "interaction_terms" = attr(prior, "interaction_terms"),
    "treatment"         = is.prior.treatment(prior),
    "independent"       = is.prior.independent(prior),
    "orthonormal"       = is.prior.orthonormal(prior),
    "meandif"           = is.prior.meandif(prior)
  )


  if(prior_info[["treatment"]]){

    if(prior_info[["levels"]] == 1){

      samples <- .as_mixed_posteriors.simple(model_samples, prior, parameter)
      samples <- matrix(samples, ncol = 1)

    }else{

      samples <- lapply(1:prior_info[["levels"]], function(i) .as_mixed_posteriors.simple(model_samples, prior, paste0(parameter, "[", i, "]")))
      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- paste0(parameter,"[",prior_info[["level_names"]][-1],"]")
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- FALSE
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["independent"]]){

    if(prior_info[["levels"]] == 1){

      samples <- .as_mixed_posteriors.simple(model_samples, prior, parameter)
      samples <- matrix(samples, ncol = 1)

    }else{

      samples <- lapply(1:prior_info[["levels"]], function(i) .as_mixed_posteriors.simple(model_samples, prior, paste0(parameter, "[", i, "]")))
      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- paste0(parameter,"[",prior_info[["level_names"]],"]")
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- FALSE
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["orthonormal"]] | prior_info[["meandif"]]){

    prior$parameter[["K"]] <- prior_info[["levels"]]
    samples <- .as_mixed_posteriors.vector(model_samples, prior, parameter)
    class(samples) <- c(class(samples), "mixed_posteriors.factor")

  }

  attr(samples, "levels")            <- prior_info[["levels"]]
  attr(samples, "level_names")       <- prior_info[["level_names"]]
  attr(samples, "interaction")       <- if(length(prior_info) == 0) FALSE else prior_info[["interaction"]]
  attr(samples, "interaction_terms") <- prior_info[["interaction_terms"]]
  attr(samples, "treatment")         <- prior_info[["treatment"]]
  attr(samples, "independent")       <- prior_info[["independent"]]
  attr(samples, "orthonormal")       <- prior_info[["orthonormal"]]
  attr(samples, "meandif")           <- prior_info[["meandif"]]

  return(samples)
}
.as_mixed_posteriors.weightfunction <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  # obtain mapping for the weight coefficients
  omega_mapping <- weightfunctions_mapping(list(prior))
  omega_cuts    <- weightfunctions_mapping(list(prior), cuts_only = TRUE)
  omega_names   <- sapply(1:(length(omega_cuts)-1), function(i) paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))

  # prepare output objects
  samples <- model_samples[, sapply(1:(length(omega_cuts)-1), function(i) paste0("omega[",i,"]"))]

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}
.as_mixed_posteriors.spike_and_slab <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  prior_variable <- prior[["variable"]]

  # prepare output objects
  if(is.prior.factor(prior_variable)){

    samples <- .as_mixed_posteriors.factor(model_samples, prior_variable, parameter)
    attr(samples, "models_ind") <- as.vector(model_samples[,paste0(parameter, "_indicator")])

  }else if(is.prior.simple(prior_variable)){

    samples <- .as_mixed_posteriors.simple(model_samples, prior_variable, parameter)
    attr(samples, "models_ind") <- as.vector(model_samples[,paste0(parameter, "_indicator")])

  }

  class(samples) <- c("mixed_posteriors.spike_and_slab", class(samples))
  attr(samples, "prior_list") <- prior

  return(samples)
}
.as_mixed_posteriors.mixture        <- function(model_samples, prior, parameter, conditional){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)


  # prepare output objects
  if(inherits(prior, "prior.bias_mixture")){

    is_PET            <- sapply(prior, is.prior.PET)
    is_PEESE          <- sapply(prior, is.prior.PEESE)
    is_weightfunction <- sapply(prior, is.prior.weightfunction)

    # prepare weightfunction parameter names
    if(any(is_weightfunction)){
      omega_mapping <- weightfunctions_mapping(prior[is_weightfunction], one_sided = TRUE)
      omega_cuts    <- weightfunctions_mapping(prior[is_weightfunction], cuts_only = TRUE, one_sided = TRUE)
      omega_names   <- sapply(1:(length(omega_cuts)-1), function(i) paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
      omega_par     <- rev(sapply(1:(length(omega_cuts)-1), function(i) paste0("omega[",i,"]")))
    }

    # deal with conditional parameters
    if(length(conditional) > 0 && any(c("PET", "PEESE", "PETPEESE", "omega") %in% conditional)){

      if("omega" %in% conditional){
        out_names <- omega_names
        par_names <- omega_par
      }else if("PETPEESE" %in% conditional){
        out_names <- par_names <- c("PET", "PEESE")
      }else if("PET" %in% conditional){
        out_names <- par_names <- "PET"
      }else if("PEESE" %in% conditional){
        out_names <- par_names <- "PEESE"
      }

    }else{

      out_names <- NULL
      par_names <- NULL

      if(any(is_weightfunction)){
        out_names <- c(out_names, omega_names)
        par_names <- c(par_names, omega_par)
      }
      if(any(is_PET)){
        out_names <- c(out_names, "PET")
        par_names <- c(par_names, "PET")
      }
      if(any(is_PEESE)){
        out_names <- c(out_names, "PEESE")
        par_names <- c(par_names, "PEESE")
      }
    }

    # select samples
    samples   <- model_samples[, par_names,drop=FALSE]
    indicator <- model_samples[,paste0(parameter, "_indicator")]

    rownames(samples) <- NULL
    colnames(samples) <- out_names
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- as.vector(indicator)
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.bias")

  }else{

    if(inherits(prior, "prior.simple_mixture")){
      samples <- .as_mixed_posteriors.simple(model_samples, prior, parameter)
    }else if(inherits(prior, "prior.factor_mixture")){
      samples <- .as_mixed_posteriors.factor(model_samples, prior, parameter)
    }
    attr(samples, "models_ind") <- as.vector(model_samples[,paste0(parameter, "_indicator")])

  }

  class(samples) <- c("mixed_posteriors.mixture", class(samples))
  attr(samples, "prior_list") <- prior

  return(samples)
}

#' @title Compute inclusion Bayes factors
#'
#' @description Computes inclusion Bayes factors based on prior model probabilities,
#' posterior model probabilities (or marginal likelihoods), and indicator whether
#' the models represent the null or alternative hypothesis.
#'
#' @param prior_probs vector of prior model probabilities
#' @param post_probs vector of posterior model probabilities
#' @param margliks vector of marginal likelihoods.
#' @param is_null logical vector of indicators whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis)
#'
#' @details Supplying \code{margliks} as the input is preferred since it is better at dealing with
#' under/overflow (posterior probabilities are very close to either 0 or 1). In case that both the
#' \code{post_probs} and \code{margliks} are supplied, the results are based on \code{margliks}.
#'
#' @return \code{inclusion_BF} returns a Bayes factor.
#'
#' @export
inclusion_BF         <- function(prior_probs, post_probs, margliks, is_null){


  if(is.numeric(is_null)){
    check_int(is_null, "is_null", lower = 1, upper = length(prior_probs), check_length = 0)
    is_null <- c(1:length(prior_probs)) %in% is_null
  }else if(is.logical(is_null)){
    check_bool(is_null, "is_null", check_length = length(prior_probs))
  }else{
    stop("'is_null' argument must be either logical vector, integer vector, or NULL.")
  }

  # deal with all null or alternative scenarios
  if(all(!is_null)){
    return(Inf)
  }else if(all(is_null)){
    return(0)
  }

  if(!missing(prior_probs) && !missing(margliks)){
    return(.inclusion_BF.margliks(prior_probs = prior_probs, margliks = margliks, is_null = is_null))
  }else if(!missing(prior_probs) && !missing(post_probs)){
    return(.inclusion_BF.probs(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null))
  }else{
    stop("'prior_probs' and either 'post_probs' or 'marglik' must be specified.")
  }
}

.inclusion_BF.probs    <- function(prior_probs, post_probs, is_null){

  check_real(prior_probs, "prior_probs", lower = 0, upper = 1, check_length = 0)
  check_real(post_probs,  "post_probs", lower = 0, upper = 1, check_length = length(prior_probs))

  if(isTRUE(all.equal(sum(post_probs[!is_null]), 1)) | isTRUE(all.equal(sum(prior_probs[!is_null]), 1))){
    return(Inf)
  }else if(isTRUE(all.equal(sum(post_probs[is_null]), 1)) | isTRUE(all.equal(sum(prior_probs[is_null]), 1))){
    return(0)
  }else{
    return(
      (sum(post_probs[!is_null]) / sum(post_probs[is_null]))  /
        (sum(prior_probs[!is_null]) / sum(prior_probs[is_null]))
    )
  }
}
.inclusion_BF.margliks <- function(prior_probs, margliks, is_null){

  check_real(prior_probs, "prior_probs", lower = 0, upper = 1, check_length = 0)
  check_real(margliks,  "margliks", check_length = length(prior_probs))

  # substract the max marglikto remove problems with overflow
  margliks <- margliks - max(margliks)

  return(
    (sum(exp(margliks[!is_null]) * prior_probs[!is_null]) / sum(exp(margliks[is_null]) * prior_probs[is_null])) /
      (sum(prior_probs[!is_null]) / sum(prior_probs[is_null]))
  )
}


#' @title Create coefficient mapping between multiple weightfunctions
#'
#' @description Creates coefficients mapping between multiple weightfunctions.
#'
#' @param prior_list list of prior distributions
#' @param cuts_only whether only p-value cuts should be returned
#' @param one_sided force one-sided output
#'
#' @return \code{weightfunctions_mapping} returns a list of indices
#' mapping the publication weights omega from the individual weightfunctions
#' into a joint weightfunction.
#'
#' @export
weightfunctions_mapping <- function(prior_list, cuts_only = FALSE, one_sided = FALSE){

  # check input
  if(!all(sapply(prior_list, is.prior.weightfunction) | sapply(prior_list, is.prior.point) | sapply(prior_list, is.prior.none)))
    stop("'priors' must be a list of weightfunction priors distributions")
  check_bool(cuts_only, "cuts_only")
  check_bool(one_sided, "one_sided")

  # extract cuts and types
  priors_cuts <- lapply(prior_list, function(prior)rev(prior[["parameters"]][["steps"]]))
  priors_type <- sapply(prior_list, function(prior)prior[["distribution"]])


  # get new cutpoint appropriate cut-points
  priors_cuts_new <- priors_cuts
  if(one_sided || any(grepl("one.sided", priors_type))){

    # translate two.sided into one.sided
    for(p in seq_along(priors_type)){
      if(grepl("two.sided", priors_type[p])){
        priors_cuts_new[[p]] <- c(priors_cuts[[p]]/2, 1 - rev(priors_cuts[[p]])/2)
      }
    }
  }

  # combine the steps
  all_cuts <- c(0, sort(unique(unlist(priors_cuts_new))), 1)

  # return the naming for summary function if only asked for labels
  if(cuts_only){
    return(all_cuts)
  }


  # get lower and upper bounds + indicies
  omega_ind <- list()
  p_bound   <- list()
  for(p in seq_along(priors_type)){
    if(!is.null(priors_cuts_new[[p]])){

      p_bound[[p]] <- list(
        l = c(0, priors_cuts_new[[p]]),
        u = c(priors_cuts_new[[p]], 1))

      if(one_sided || any(grepl("one.sided", priors_type))){
        if(grepl("two.sided", priors_type[p])){
          omega_ind[[p]] <- rev(c( (length(priors_cuts[[p]]) + 1):2, 1:(length(priors_cuts[[p]]) + 1) ))
        }else if(grepl("one.sided", priors_type[p])){
          omega_ind[[p]] <- rev(1:(length(priors_cuts[[p]]) + 1))
        }
      }else{
        if(grepl("two.sided", priors_type[p])){
          omega_ind[[p]] <- rev(1:(length(priors_cuts[[p]]) + 1))
        }
      }

    }
  }

  # create mapping to weights
  omega_mapping <- list()
  for(p in seq_along(priors_type)){
    if(is.prior.weightfunction(prior_list[[p]])){
      omega_mapping[[p]] <- sapply(1:(length(all_cuts)-1), function(i)
        omega_ind[[p]][all_cuts[i] >= p_bound[[p]]$l & all_cuts[i+1] <= p_bound[[p]]$u]
      )
    }
  }


  return(omega_mapping)
}
