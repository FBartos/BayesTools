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
#' odds \code{prior_odds}
#' @param parameters vector of parameters names for which inference should
#' be drawn
#' @param is_null_list list with entries for each parameter carrying either
#' logical vector of indicators specifying whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis)
#' @param prior_odds vector of prior model odds
#' @param margliks vector of marginal likelihoods
#' @param is_null logical vector of indicators specifying whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis)
#' @param conditional whether prior and posterior model probabilities should
#' be returned only for the conditional model. Defaults to \code{FALSE}
#'
#' @name ensemble_inference
#' @export compute_inference
#' @export ensemble_inference
#' @export models_inference
#'
#' @seealso [mix_posteriors] [BayesTools_ensemble_tables]
NULL

#' @rdname ensemble_inference
compute_inference <- function(prior_odds, margliks, is_null = NULL, conditional = FALSE){

  check_real(prior_odds, "prior_odds", lower = 0, check_length = 0)
  check_real(margliks,   "margliks", check_length = length(prior_odds))
  check_bool(conditional, "conditional")
  if(!is.null(is_null)){
    if(is.numeric(is_null)){
      check_int(is_null, "is_null", lower = 0, upper = length(prior_odds), check_length = 0)
      is_null <- c(1:length(prior_odds)) %in% is_null
    }else if(is.logical(is_null)){
      check_bool(is_null, "is_null", check_length = length(prior_odds))
    }else{
      stop("'is_null' argument must be either logical vector, integer vector, or NULL.")
    }
  }else{
    is_null <- rep(FALSE, length(prior_odds))
  }

  prior_probs <- prior_odds / sum(prior_odds)
  post_probs  <- unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
  BF          <- inclusion_BF(prior_probs, post_probs, is_null)

  if(conditional){
    prior_probs <- ifelse(is_null, 0, prior_odds) / sum(prior_odds[!is_null])
    post_probs  <- unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
  }

  output <- list(
    prior_probs = prior_probs,
    post_probs  = post_probs,
    BF          = BF
  )
  attr(output, "is_null")    <- is_null
  attr(output, "condtional") <- conditional

  return(output)
}

#' @rdname ensemble_inference
ensemble_inference <- function(model_list, parameters, is_null_list, conditional = FALSE){

  # check input
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("marglik", "prior_odds"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m)inherits(m[["marglik"]], what = "bridge"))))
    stop("model_list:marglik must contain 'bridgesampling' marginal likelihoods")
  sapply(model_list, function(m)check_real(m[["prior_odds"]], "model_list:prior_odds", lower = 0))


  # extract the object
  margliks   <- sapply(model_list, function(m)m[["marglik"]][["logml"]])
  prior_odds <- sapply(model_list, function(m)m[["prior_odds"]])

  out <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter <- parameters[p]
    temp_is_null   <- is_null_list[[p]]

    out[[temp_parameter]] <- compute_inference(prior_odds = prior_odds, margliks = margliks, is_null = temp_is_null, conditional = conditional)
  }

  attr(out, "conditional") <- conditional
  return(out)
}

#' @rdname ensemble_inference
models_inference   <- function(model_list){

  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("marglik", "prior_odds"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m)inherits(m[["marglik"]], what = "bridge"))))
    stop("model_list:marglik must contain 'bridgesampling' marginal likelihoods")
  sapply(model_list, function(m)check_real(m[["prior_odds"]], "model_list:prior_odds", lower = 0))

  margliks    <- sapply(model_list, function(model)model[["marglik"]][["logml"]])
  margliks    <- ifelse(is.na(margliks), -Inf, margliks)
  prior_odds  <- sapply(model_list, function(model)model[["prior_odds"]])
  prior_probs <- prior_odds / sum(prior_odds)
  post_probs  <- unname(bridgesampling::post_prob(margliks, prior_prob = prior_probs))
  incl_BF     <- sapply(seq_along(model_list), function(i) (post_probs[i] / (1 - post_probs[i])) / (prior_probs[i] / (1 - prior_probs[i])))

  for(i in seq_along(model_list)){
    model_list[[i]][["inference"]] <- list(
      "m_number"     = i,
      "marglik"      = margliks[i],
      "prior_prob"   = prior_probs[i],
      "post_probs"   = post_probs[i],
      "inclusion_BF" = incl_BF[i]
    )
  }

  return(model_list)
}

#' @title Model-average posterior distributions
#'
#' @description Model-averages posterior distributions based either
#'
#' on (1) a list of models, vector of parameters, and a list of
#' indicators the models represent the null or alternative hypothesis
#' for each parameter, or (2) on model fits, list of priors for a given
#' parameter, parameter name, and posterior probabilities
#'
#' @param fits list with either runjags or rstan fitted models
#' @param priors list of prior distributions for the model-averaged
#' parameter for each passed fit
#' @param parameter name of the parameter to be model-averaged
#' @param post_probs vector of posterior probabilities for each fit
#' @param seed integer specifying seed for sampling posteriors for
#' model averaging. Defaults to \code{1}.
#' @param n_samples number of samples to be drawn for the model-averaged
#' posterior distribution
#'
#' @inheritParams ensemble_inference
#' @name mix_posteriors
#' @export mix_posteriors
#' @export mix_posteriors.simple
#' @export mix_posteriors.weightfunction
#'
#' @seealso [ensemble_inference] [BayesTools_ensemble_tables]
NULL

#' @rdname mix_posteriors
mix_posteriors <- function(model_list, parameters, is_null_list, conditional = FALSE, seed = NULL, n_samples = 10000){

  # check input
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  check_real(seed, "seed", allow_NULL = TRUE)
  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("fit", "marglik", "priors", "prior_odds"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m)inherits(m[["fit"]], what = "runjags")) | sapply(model_list, function(m)inherits(m[["fit"]], what = "rstan")) | sapply(model_list, function(m)inherits(m[["fit"]], what = "null_model"))))
    stop("model_list:fit must contain 'runjags' or 'rstan' models")
  if(!all(sapply(model_list, function(m)inherits(m[["marglik"]], what = "bridge"))))
    stop("model_list:marglik must contain 'bridgesampling' marginal likelihoods")
  if(!all(unlist(sapply(model_list, function(m)sapply(m[["priors"]], function(p)is.prior(p))))))
    stop("model_list:priors must contain 'BayesTools' priors")
  sapply(model_list, function(m)check_real(m[["prior_odds"]], "model_list:prior_odds", lower = 0))


  # extract the object
  fits       <- lapply(model_list, function(m)m[["fit"]])
  margliks   <- sapply(model_list, function(m)m[["marglik"]][["logml"]])
  priors     <- lapply(model_list, function(m)m[["priors"]])
  prior_odds <- sapply(model_list, function(m)m[["prior_odds"]])

  inference  <- ensemble_inference(model_list, parameters, is_null_list, conditional)
  out <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter    <- parameters[p]
    temp_inference    <- inference[[temp_parameter]]
    temp_priors       <- lapply(priors, function(p)p[[temp_parameter]])

    if(all(sapply(temp_priors, is.prior.simple) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }
      out[[temp_parameter]] <- mix_posteriors.simple(fits, temp_priors, temp_parameter, temp_inference$post_probs, seed, n_samples)

    }else if(all(sapply(temp_priors, is.prior.weightfunction) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.prior.none) | sapply(temp_priors, is.null))){

      # replace missing priors with default prior: none
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior_none()
        }
      }
      out[[temp_parameter]] <- mix_posteriors.weightfunction(fits, temp_priors, temp_parameter, temp_inference$post_probs, seed, n_samples)

    }

  }

  return(out)
}

#' @rdname mix_posteriors
mix_posteriors.simple         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "rstan") | sapply(fits, inherits, what = "null_model")))
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
    }else if(inherits(fits[[i]], "rstan")){
      model_samples <- rstan::extract(fits[[i]])
      par_names     <- names(model_samples)
      par_dims      <- sapply(model_samples, function(s)if(is.matrix(s)) ncol(s) else 1)
      par_names     <- unlist(sapply(seq_along(par_names), function(p){
        if(par_dims[p] == 1){
          return(par_names[p])
        }else{
          return(paste0(par_names[p], "[", 1:par_dims[p],"]"))
        }
      }))
      model_samples <- data.frame(do.call(cbind, model_samples))
      colnames(model_samples) <- par_names
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
  return(samples)
}

#' @rdname mix_posteriors
mix_posteriors.weightfunction <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000){

  # check input
  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  check_real(seed, "seed", allow_NULL = TRUE)
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "rstan") | sapply(fits, inherits, what = "null_model")))
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
    }else if(inherits(fits[[i]], "rstan")){
      model_samples <- rstan::extract(fits[[i]])
      par_names     <- names(model_samples)
      par_dims      <- sapply(model_samples, function(s)if(is.matrix(s)) ncol(s) else 1)
      par_names     <- unlist(sapply(seq_along(par_names), function(p){
        if(par_dims[p] == 1){
          return(par_names[p])
        }else{
          return(paste0(par_names[p], "[", 1:par_dims[p],"]"))
        }
      }))
      model_samples <- data.frame(do.call(cbind, model_samples))
      colnames(model_samples) <- par_names
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
  return(samples)
}



#' @title Compute inclusion Bayes factors
#'
#' @description Computes inclusion Bayes factors based on prior and posterior
#' model probabilities and indicator whether the models represent the null or
#' alternative hypothesis
#'
#' @param prior_probs vector of prior model probabilities
#' @param post_probs vector of posterior model probabilities
#' @param is_null logical vector of indicators whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis)
#'
#' @export
inclusion_BF     <- function(prior_probs, post_probs, is_null){

  check_real(prior_probs, "prior_probs", lower = 0, upper = 1, check_length = 0)
  check_real(post_probs,  "post_probs", lower = 0, upper = 1, check_length = length(prior_probs))
  if(is.numeric(is_null)){
    check_int(is_null, "is_null", lower = 1, upper = length(prior_probs), check_length = 0)
    is_null <- c(1:length(prior_probs)) %in% is_null
  }else if(is.logical(is_null)){
    check_bool(is_null, "is_null", check_length = length(prior_probs))
  }else{
    stop("'is_null' argument must be either logical vector, integer vector, or NULL.")
  }

  if(all(!is_null)){
    return(Inf)
  }else if(all(is_null)){
    return(0)
  }else{
    return(
      (sum(post_probs[!is_null]) / sum(post_probs[is_null]))  /
        (sum(prior_probs[!is_null]) / sum(prior_probs[is_null]))
    )
  }
}


#' @title Create coefficient mapping between multiple weightfunctions
#'
#' @description Creates coefficients mapping between multiple weightfunctions.
#'
#' @param prior_list list of prior distributions
#' @param cuts_only whether only p-value cuts should be returned
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

      if(grepl("two.sided", priors_type[p])){
        omega_ind[[p]] <- rev(c( (length(priors_cuts[[p]]) + 1):2, 1:(length(priors_cuts[[p]]) + 1) ))
      }else if(grepl("one.sided", priors_type[p])){
        omega_ind[[p]] <- rev(1:(length(priors_cuts[[p]]) + 1))
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


#' @title Format Bayes factor
#'
#' @description Formats Bayes factor
#'
#' @param BF Bayes factor(s)
#' @param logBF log(BF)
#' @param BF01 1/BF
#'
#' @export
format_BF <- function(BF, logBF = FALSE, BF01 = FALSE){

  check_real(BF, "BF", lower = 0, check_length = FALSE)
  check_bool(logBF, "logBF")
  check_bool(BF01,  "BF01")

  name <- "BF10"

  if(BF01){
    BF   <- 1/BF
    name <- "BF01"
  }
  if(logBF){
    BF   <- log(BF)
    name <- paste0("log(", name, ")")
  }

  attr(BF, "name")  <- name
  attr(BF, "logBF") <- logBF
  attr(BF, "BF01")  <- BF01

  return(BF)
}
