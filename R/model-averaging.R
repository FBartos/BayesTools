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
      parameter_name <- gsub(paste0(formula_parameter, "_"), paste0("(", formula_parameter, ") "), parameter_name)
      parameter_name <- gsub("__xXx__", ":", parameter_name)
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
#' indicators the models represent the null or alternative hypothesis
#' for each parameter.
#'
#' @param seed integer specifying seed for sampling posteriors for
#' model averaging. Defaults to \code{1}.
#' @param n_samples number of samples to be drawn for the model-averaged
#' posterior distribution
#' @inheritParams ensemble_inference
#'
#' @return \code{mix_posteriors} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [ensemble_inference] [BayesTools_ensemble_tables]
#'
#' @name mix_posteriors
#' @export
mix_posteriors <- function(model_list, parameters, is_null_list, conditional = FALSE, seed = NULL, n_samples = 10000){

  # check input
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  check_real(seed, "seed", allow_NULL = TRUE)
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

  return(out)
}

.mix_posteriors.simple         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.simple) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of simple priors")


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
      samples <- c(samples, rng(priors[[i]], length(temp_ind)))
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
      samples <- rbind(samples, matrix(rng(priors[[i]], 1), nrow = length(temp_ind), ncol = K))
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
        "levels"      = .get_prior_factor_levels(p),
        "level_names" = .get_prior_factor_level_names(p),
        "interaction" = .is_prior_interaction(p),
        "treatment"   = is.prior.dummy(p),
        "independent" = is.prior.independent(p),
        "orthonormal" = is.prior.orthonormal(p),
        "meandif"     = is.prior.meandif(p)
      ))
    }else{
      stop("unsupported prior type")
    }
  })
  priors_info <- priors_info[!sapply(priors_info, isFALSE)]
  if(length(priors_info) >= 2 && any(!unlist(lapply(priors_info, function(i) all.equal(i, priors_info[[1]]))))){
    stop("non-matching prior factor type specifications")
  }
  priors_info <- priors_info[[1]]

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

  attr(samples, "levels")      <- priors_info[["levels"]]
  attr(samples, "level_names") <- priors_info[["level_names"]]
  attr(samples, "interaction") <- priors_info[["interaction"]]
  attr(samples, "treatment")   <- priors_info[["treatment"]]
  attr(samples, "independent") <- priors_info[["independent"]]
  attr(samples, "orthonormal") <- priors_info[["orthonormal"]]
  attr(samples, "meandif")     <- priors_info[["meandif"]]

  return(samples)
}
.mix_posteriors.weightfunction <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
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
#'
#' @return \code{weightfunctions_mapping} returns a list of indices
#' mapping the publication weights omega from the individual weightfunctions
#' into a joint weightfunction.
#'
#' @export
weightfunctions_mapping <- function(prior_list, cuts_only = FALSE){

  # check input
  if(!all(sapply(prior_list, is.prior.weightfunction) | sapply(prior_list, is.prior.point) | sapply(prior_list, is.prior.none)))
    stop("'priors' must be a list of weightfunction priors distributions")
  check_bool(cuts_only, "cuts_only")


  # extract cuts and types
  priors_cuts <- lapply(prior_list, function(prior)rev(prior[["parameters"]][["steps"]]))
  priors_type <- sapply(prior_list, function(prior)prior[["distribution"]])


  # get new cutpoint appropriate cut-points
  priors_cuts_new <- priors_cuts
  if(any(grepl("one.sided", priors_type))){

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

      if(any(grepl("one.sided", priors_type))){
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
