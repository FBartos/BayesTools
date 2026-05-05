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
  if(sum(prior_weights) <= 0){
    stop("At least one prior model weight must be positive.", call. = FALSE)
  }

  prior_weights / sum(prior_weights)
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

  # Use one shared sampling seed so formula terms are drawn from aligned
  # posterior rows across parameters.
  if(!is.null(seed)){
    set.seed(seed)
    common_sample_seed <- seed
  }else{
    common_sample_seed <- sample(.Machine$integer.max, 1)
  }

  out <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter    <- parameters[p]
    temp_inference    <- inference[[temp_parameter]]
    temp_priors       <- lapply(priors, function(p) p[[temp_parameter]])

    if(any(sapply(temp_priors, is.prior.weightfunction)) && all(sapply(temp_priors, is.prior.weightfunction) | sapply(temp_priors, .is_prior_weightfunction_null) | sapply(temp_priors, is.null))){
      # weightfunctions:

      # replace missing priors with default prior: none
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior_none()
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        temp_priors[[i]] <- .set_prior_model_weight(temp_priors[[i]], temp_inference$prior_probs[i])
      }

      out[[temp_parameter]] <- .mix_posteriors.weightfunction(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples)

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
        temp_priors[[i]] <- .set_prior_model_weight(temp_priors[[i]], temp_inference$prior_probs[i])
      }

      out[[temp_parameter]] <- .mix_posteriors.factor(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples)

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
        temp_priors[[i]] <- .set_prior_model_weight(temp_priors[[i]], temp_inference$prior_probs[i])
      }

      out[[temp_parameter]] <- .mix_posteriors.vector(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples)

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
        temp_priors[[i]] <- .set_prior_model_weight(temp_priors[[i]], temp_inference$prior_probs[i])
      }

      out[[temp_parameter]] <- .mix_posteriors.simple(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples)

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



  # When seed is NULL, continue using the current RNG state so mix_posteriors()
  # can coordinate draws across parameters.
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  samples <- NULL
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
  for(i in seq_along(fits)[sample_counts > 0]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
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
    temp_ind <- sample(nrow(model_samples), sample_counts[i], replace = TRUE)

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

  # When seed is NULL, continue using the current RNG state so mix_posteriors()
  # can coordinate draws across parameters.
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  K <- unique(sapply(priors[sapply(priors, is.prior.vector)], function(p) p$parameters[["K"]]))
  if(length(K) != 1)
    stop("all vector priors must be of the same length")

  samples    <- matrix(nrow = 0, ncol = K)
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
  for(i in seq_along(fits)[sample_counts > 0]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
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
    temp_ind <- sample(nrow(model_samples), sample_counts[i], replace = TRUE)

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
        "term_components"   = attr(p, "term_components"),
        "factor_terms"      = attr(p, "factor_terms"),
        "factor_contrasts"  = attr(p, "factor_contrasts"),
        "factor_design"     = attr(p, "factor_design"),
        "factor_cell_names" = attr(p, "factor_cell_names"),
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
  attr(samples, "term_components")   <- priors_info[["term_components"]]
  attr(samples, "factor_terms")      <- priors_info[["factor_terms"]]
  attr(samples, "factor_contrasts")  <- priors_info[["factor_contrasts"]]
  attr(samples, "factor_design")     <- priors_info[["factor_design"]]
  attr(samples, "factor_cell_names") <- priors_info[["factor_cell_names"]]
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
  if(!all(sapply(priors, is.prior.weightfunction) | sapply(priors, .is_prior_weightfunction_null)))
    stop("'priors' must be a list of weightfunction priors or point(1)/none null priors")


  # When seed is NULL, continue using the current RNG state so mix_posteriors()
  # can coordinate draws across parameters.
  if(!is.null(seed)){
    set.seed(seed)
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
  sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
  for(i in seq_along(fits)[sample_counts > 0]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
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
    temp_ind <- sample(nrow(model_samples), sample_counts[i], replace = TRUE)

    if(.is_prior_weightfunction_null(priors[[i]])){
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
#' @param transform_scaled whether to transform samples from standardized (scaled) to
#' original (unscaled) scale. When \code{TRUE}, posterior samples are
#' transformed, and the result can be directly passed to [plot_posterior] which will
#' automatically detect the transformation and use transformed deterministic prior densities.
#' Requires a model fitted with \code{formula_scale_list}. Defaults to \code{FALSE}.
#' @param n_prior_samples controls the numerical grid used for transformed
#' prior densities when \code{transform_scaled = TRUE}. Defaults to 10000.
#' @inheritParams ensemble_inference
#'
#' @return \code{as_mix_posteriors} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors]
#'
#' @name as_mixed_posteriors
#' @export
as_mixed_posteriors <- function(model, parameters, conditional = NULL, conditional_rule = "AND", force_plots = FALSE,
                                 transform_scaled = FALSE, n_prior_samples = 10000){

  # check input
  if(!inherits(model, "BayesTools_fit"))
    stop("'model' must be a 'BayesTools_fit'")
  check_char(parameters, "parameters", check_length = FALSE)
  check_char(conditional, "conditional", check_length = FALSE, allow_values = c(parameters, "PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null"), allow_NULL = TRUE)
  check_char(conditional_rule, "conditional_rule", allow_values = c("AND", "OR"))
  check_bool(transform_scaled, "transform_scaled")
  check_int(n_prior_samples, "n_prior_samples", lower = 1)

  # extract the list of priors
  priors <- attr(model, "prior_list")
  prior_density_priors <- priors

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

      # special cases for PET / PEESE / PET-PEESE / selection / p-hacking within the bias parameter
      if (!is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])) {
        bias_branch_info <- .selection_prior_branch_info(priors[["bias"]])
        if(parameter == "PET"){
          is_PET <- sapply(priors[["bias"]], is.prior.PET)
          return(model_samples[, "bias_indicator"] %in% which(is_PET))
        }
        if(parameter == "PEESE"){
          is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
          return(model_samples[, "bias_indicator"] %in% which(is_PEESE))
        }
        if(parameter == "PETPEESE"){
          is_PET   <- sapply(priors[["bias"]], is.prior.PET)
          is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
          return(model_samples[, "bias_indicator"] %in% which(is_PET | is_PEESE))
        }
        if(parameter == "omega"){
          has_selection <- vapply(bias_branch_info, function(x) !is.null(x$selection), logical(1))
          return(model_samples[, "bias_indicator"] %in% which(has_selection))
        }
        if(parameter %in% c("phacking", "alpha", "pi_null")){
          has_phacking <- vapply(bias_branch_info, function(x) !is.null(x$phacking), logical(1))
          return(model_samples[, "bias_indicator"] %in% which(has_phacking))
        }
      }

      if(parameter == "omega" && any(vapply(priors, function(x){
        is.prior.weightfunction(x) || (is_prior_bias(x) && !is.null(x$selection))
      }, logical(1)))){
        return(rep(TRUE, nrow(model_samples)))
      }

      if(parameter %in% c("phacking", "alpha", "pi_null") && any(vapply(priors, function(x){
        is_prior_phacking(x) || (is_prior_bias(x) && !is.null(x$phacking))
      }, logical(1)))){
        return(rep(TRUE, nrow(model_samples)))
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

        warning(sprintf("The parameter '%s' is not a conditional parameter. All samples are assumed to compe from the conditional posterior distribution.", parameter), call. = FALSE, immediate. = TRUE)
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
    if(length(conditional) == 1){

      if (conditional %in% c("bias", "PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null") && !is.null(priors[["bias"]]) && is.prior.mixture(priors[["bias"]])) {

        # special cases for PET / PEESE / PET-PEESE / selection / p-hacking
        bias_branch_info <- .selection_prior_branch_info(priors[["bias"]])
        if(conditional == "bias"){
          components <- attr(priors[["bias"]], "components")
          if(is.null(components) || length(components) != length(priors[["bias"]])){
            is_null <- sapply(priors[["bias"]], is.prior.none)
          }else{
            is_null <- components == "null"
          }
          for(i in seq(along = is_null)){
            if(is_null[i]){
              priors[["bias"]][[i]][["prior_weights"]] <- 0
            }
          }
        }else if(conditional == "PET"){
          is_PET <- sapply(priors[["bias"]], is.prior.PET)
          for(i in seq(along = is_PET)){
            if(!is_PET[i]){
              priors[["bias"]][[i]][["prior_weights"]] <- 0
            }
          }
        }else if(conditional == "PEESE"){
          is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
          for(i in seq(along = is_PEESE)){
            if(!is_PEESE[i]){
              priors[["bias"]][[i]][["prior_weights"]] <- 0
            }
          }
        }else if(conditional == "PETPEESE"){
          is_PET   <- sapply(priors[["bias"]], is.prior.PET)
          is_PEESE <- sapply(priors[["bias"]], is.prior.PEESE)
          for(i in seq(along = is_PET)){
            if(!(is_PET[i] || is_PEESE[i])){
              priors[["bias"]][[i]][["prior_weights"]] <- 0
            }
          }
        }else if(conditional == "omega"){
          has_selection <- vapply(bias_branch_info, function(x) !is.null(x$selection), logical(1))
          for(i in seq(along = has_selection)){
            if(!has_selection[i]){
              priors[["bias"]][[i]][["prior_weights"]] <- 0
            }
          }
        }else if(conditional %in% c("phacking", "alpha", "pi_null")){
          has_phacking <- vapply(bias_branch_info, function(x) !is.null(x$phacking), logical(1))
          for(i in seq(along = has_phacking)){
            if(!has_phacking[i]){
              priors[["bias"]][[i]][["prior_weights"]] <- 0
            }
          }
        }

        # propagate the prior weights to the mixture prior itself
        attr(priors[["bias"]], "prior_weights") <- sapply(priors[["bias"]], \(x) x[["prior_weights"]])

      }else if(is.prior.mixture(priors[[conditional]])){

        components <- attr(priors[[conditional]], "components")
        for(i in seq(along = components)){
          if(components[i] == "null"){
            priors[[conditional]][[i]][["prior_weights"]] <- 0
          }
        }

        # propagate the prior weights to the mixture prior itself
        attr(priors[[conditional]], "prior_weights") <- sapply(priors[[conditional]], \(x) x[["prior_weights"]])
      }

    }
  }

  # extract formula_scale early for transform_scaled support
  formula_scale <- attr(model, "formula_scale")

  # apply scale transformation to posterior samples if requested
  if(transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0){
    model_samples <- transform_scale_samples(model_samples, formula_scale)
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

    }else if(is_prior_phacking(temp_prior)){
      # p-hacking priors
      out[[temp_parameter]] <- .as_mixed_posteriors.phacking(model_samples, temp_prior, temp_parameter)

    }else if(is_prior_bias(temp_prior)){
      # composed publication-bias priors
      out[[temp_parameter]] <- .as_mixed_posteriors.bias(model_samples, temp_prior, temp_parameter, conditional)

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

  # propagate formula_scale attribute for transform_scaled support
  if(!is.null(formula_scale)){
    attr(out, "formula_scale") <- formula_scale
  }

  prior_context_priors      <- prior_density_priors
  prior_context_conditional <- conditional
  special_conditionals      <- c("PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null")
  if(length(conditional) == 1 &&
     conditional %in% special_conditionals &&
     !conditional %in% names(prior_density_priors)){
    prior_context_priors      <- priors
    prior_context_conditional <- NULL
  }

  # generate and store transformed prior densities if requested
  if(transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0){
    prior_densities <- .generate_transformed_prior_densities(
      prior_list       = prior_context_priors,
      column_names     = colnames(model_samples),
      n_grid           = n_prior_samples,
      formula_scale    = formula_scale,
      conditional      = prior_context_conditional,
      conditional_rule = conditional_rule
    )
    attr(out, "prior_densities")       <- prior_densities
    attr(out, "prior_density_context") <- attr(prior_densities, "context")
    attr(out, "transform_scaled")      <- TRUE
  }else{
    attr(out, "prior_density_context") <- .prior_density_build_context(
      prior_list       = prior_context_priors,
      column_names     = colnames(model_samples),
      n_grid           = n_prior_samples,
      conditional      = prior_context_conditional,
      conditional_rule = conditional_rule
    )
  }

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
  attr(samples, "models_ind") <- rep(1, length(samples))
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
  attr(samples, "models_ind") <- rep(1, nrow(samples))
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
    "term_components"   = attr(prior, "term_components"),
    "factor_terms"      = attr(prior, "factor_terms"),
    "factor_contrasts"  = attr(prior, "factor_contrasts"),
    "factor_design"     = attr(prior, "factor_design"),
    "factor_cell_names" = attr(prior, "factor_cell_names"),
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

    level_names <- prior_info[["level_names"]]
    if(is.list(level_names)){
      level_names <- lapply(level_names, function(x) x[-1])
    }else{
      level_names <- level_names[-1]
    }

    rownames(samples) <- NULL
    colnames(samples) <- .format_factor_level_parameter_names(parameter, level_names, ncol(samples))
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- rep(1, nrow(samples))
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
    colnames(samples) <- .format_factor_level_parameter_names(parameter, prior_info[["level_names"]], ncol(samples))
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- rep(1, nrow(samples))
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
  attr(samples, "term_components")   <- prior_info[["term_components"]]
  attr(samples, "factor_terms")      <- prior_info[["factor_terms"]]
  attr(samples, "factor_contrasts")  <- prior_info[["factor_contrasts"]]
  attr(samples, "factor_design")     <- prior_info[["factor_design"]]
  attr(samples, "factor_cell_names") <- prior_info[["factor_cell_names"]]
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
  omega_info    <- .weightfunction_mapping_info(list(prior))
  omega_mapping <- omega_info$mapping
  omega_names   <- omega_info$names
  omega_par     <- omega_info$pars

  # prepare output objects
  samples <- model_samples[, omega_par, drop = FALSE]

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  samples <- .weightfunction_set_omega_context(samples, omega_info)
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}
.as_mixed_posteriors.phacking      <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  par_names <- intersect(.phacking_report_parameter(prior), colnames(model_samples))
  samples   <- model_samples[, par_names, drop = FALSE]

  rownames(samples) <- NULL
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.phacking")

  return(samples)
}
.as_mixed_posteriors.bias          <- function(model_samples, prior, parameter, conditional){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  branch_info   <- .selection_prior_branch_info(prior)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

  out_names <- NULL
  par_names <- NULL

  if(any(has_selection)){
    selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
    omega_info       <- .weightfunction_mapping_info(selection_priors, one_sided = TRUE)
    omega_names      <- omega_info$names
    omega_par        <- omega_info$pars
  }
  if(any(has_phacking)){
    phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
    phacking_par    <- .selection_phacking_report_parameters(phacking_priors)
    phacking_names  <- phacking_par
  }

  if(length(conditional) > 0 && any(c("omega", "phacking", "alpha", "pi_null") %in% conditional)){
    if("omega" %in% conditional && any(has_selection)){
      out_names <- c(out_names, omega_names)
      par_names <- c(par_names, omega_par)
    }
    if("phacking" %in% conditional && any(has_phacking)){
      out_names <- c(out_names, phacking_names)
      par_names <- c(par_names, phacking_par)
    }
    if("alpha" %in% conditional && any(has_phacking)){
      out_names <- c(out_names, "alpha")
      par_names <- c(par_names, "alpha")
    }
    if("pi_null" %in% conditional && any(has_phacking)){
      out_names <- c(out_names, "pi_null")
      par_names <- c(par_names, "pi_null")
    }
  }else{
    if(any(has_selection)){
      out_names <- c(out_names, omega_names)
      par_names <- c(par_names, omega_par)
    }
    if(any(has_phacking)){
      out_names <- c(out_names, phacking_names)
      par_names <- c(par_names, phacking_par)
    }
  }

  if(is.null(par_names)){
    par_names <- character()
    out_names <- character()
  }
  keep_unique <- !duplicated(par_names)
  par_names <- par_names[keep_unique]
  out_names <- out_names[keep_unique]
  keep_par <- par_names %in% colnames(model_samples)
  par_names <- par_names[keep_par]
  out_names <- out_names[keep_par]
  samples   <- model_samples[, par_names, drop = FALSE]

  rownames(samples) <- NULL
  colnames(samples) <- out_names
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  if(any(has_selection)){
    samples <- .weightfunction_set_omega_context(samples, omega_info)
  }
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.bias")

  return(samples)
}
.as_mixed_posteriors.spike_and_slab <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  prior_variable <- .get_spike_and_slab_variable(prior)

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
    branch_info       <- .selection_prior_branch_info(prior)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))

    # prepare weightfunction parameter names
    if(any(has_selection)){
      selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
      omega_info    <- .weightfunction_mapping_info(selection_priors, one_sided = TRUE)
      omega_mapping <- omega_info$mapping
      omega_names   <- omega_info$names
      omega_par     <- omega_info$pars
    }
    if(any(has_phacking)){
      phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
      phacking_par    <- .selection_phacking_report_parameters(phacking_priors)
      phacking_names  <- phacking_par
    }

    # deal with conditional parameters
    if(length(conditional) > 0 && any(c("PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null") %in% conditional)){

      out_names <- NULL
      par_names <- NULL

      if("omega" %in% conditional && any(has_selection)){
        out_names <- c(out_names, omega_names)
        par_names <- c(par_names, omega_par)
      }
      if("phacking" %in% conditional && any(has_phacking)){
        out_names <- c(out_names, phacking_names)
        par_names <- c(par_names, phacking_par)
      }
      if("alpha" %in% conditional && any(has_phacking)){
        out_names <- c(out_names, "alpha")
        par_names <- c(par_names, "alpha")
      }
      if("pi_null" %in% conditional && any(has_phacking)){
        out_names <- c(out_names, "pi_null")
        par_names <- c(par_names, "pi_null")
      }
      if("PETPEESE" %in% conditional){
        # subset in case only PET/PEESE is supplied
        out_names <- c(out_names, colnames(model_samples)[colnames(model_samples) %in% c("PET", "PEESE")])
        par_names <- c(par_names, colnames(model_samples)[colnames(model_samples) %in% c("PET", "PEESE")])
      }
      if("PET" %in% conditional && any(is_PET)){
        out_names <- c(out_names, "PET")
        par_names <- c(par_names, "PET")
      }
      if("PEESE" %in% conditional && any(is_PEESE)){
        out_names <- c(out_names, "PEESE")
        par_names <- c(par_names, "PEESE")
      }

    }else{

      out_names <- NULL
      par_names <- NULL

      if(any(has_selection)){
        out_names <- c(out_names, omega_names)
        par_names <- c(par_names, omega_par)
      }
      if(any(has_phacking)){
        out_names <- c(out_names, phacking_names)
        par_names <- c(par_names, phacking_par)
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
    if(is.null(par_names)){
      par_names <- character()
      out_names <- character()
    }
    keep_unique <- !duplicated(par_names)
    par_names <- par_names[keep_unique]
    out_names <- out_names[keep_unique]
    keep_par <- par_names %in% colnames(model_samples)
    par_names <- par_names[keep_par]
    out_names <- out_names[keep_par]
    samples   <- model_samples[, par_names,drop=FALSE]
    indicator <- model_samples[,paste0(parameter, "_indicator")]

    rownames(samples) <- NULL
    colnames(samples) <- out_names
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- as.vector(indicator)
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    if(any(has_selection)){
      samples <- .weightfunction_set_omega_context(samples, omega_info)
    }
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
#' If the prior probability of either the null or alternative hypothesis is
#' zero, the Bayes factor is undefined and \code{NA} is returned.
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

  prior_alt  <- sum(prior_probs[!is_null])
  prior_null <- sum(prior_probs[is_null])
  post_alt   <- sum(post_probs[!is_null])
  post_null  <- sum(post_probs[is_null])

  if(isTRUE(all.equal(prior_alt, 0)) || isTRUE(all.equal(prior_null, 0))){
    return(NA_real_)
  }

  if(isTRUE(all.equal(post_alt, 1))){
    return(Inf)
  }else if(isTRUE(all.equal(post_null, 1))){
    return(0)
  }else{
    return(
      (post_alt / post_null) / (prior_alt / prior_null)
    )
  }
}
.inclusion_BF.margliks <- function(prior_probs, margliks, is_null){

  check_real(prior_probs, "prior_probs", lower = 0, upper = 1, check_length = 0)
  check_real(margliks,  "margliks", check_length = length(prior_probs))

  prior_alt  <- sum(prior_probs[!is_null])
  prior_null <- sum(prior_probs[is_null])

  if(isTRUE(all.equal(prior_alt, 0)) || isTRUE(all.equal(prior_null, 0))){
    return(NA_real_)
  }

  margliks <- .model_averaging_margliks(margliks, prior_probs)

  active <- prior_probs > 0 & is.finite(margliks)
  if(!any(active & !is_null)){
    return(0)
  }
  if(!any(active & is_null)){
    return(Inf)
  }

  # subtract the max among positive-prior finite models to avoid overflow.
  margliks <- margliks - max(margliks[active])

  alt_ind  <- active & !is_null
  null_ind <- active & is_null

  alt_marginal  <- sum(exp(margliks[alt_ind])  * prior_probs[alt_ind])
  null_marginal <- sum(exp(margliks[null_ind]) * prior_probs[null_ind])

  if(alt_marginal == 0 && null_marginal == 0){
    return(NaN)
  }
  if(alt_marginal == 0){
    return(0)
  }
  if(null_marginal == 0){
    return(Inf)
  }

  return(
    (alt_marginal / null_marginal) / (prior_alt / prior_null)
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
  if(!all(sapply(prior_list, is.prior.weightfunction) | sapply(prior_list, .is_prior_weightfunction_null)))
    stop("'priors' must be a list of weightfunction priors or point(1)/none null priors")
  check_bool(cuts_only, "cuts_only")
  check_bool(one_sided, "one_sided")

  force_one_sided <- one_sided || any(sapply(prior_list, function(prior){
    is.prior.weightfunction(prior) && prior$side == "one-sided"
  }))

  prior_expansions <- lapply(prior_list, function(prior){
    if(!is.prior.weightfunction(prior)){
      return(NULL)
    }
    .weightfunction_mapping_expansion(prior, force_one_sided)
  })

  all_cuts <- .weightfunction_unique_cuts(unlist(lapply(prior_expansions, function(expansion){
    if(is.null(expansion)) NULL else expansion$cuts
  })))
  if(length(all_cuts) == 0L){
    all_cuts <- c(0, 1)
  }

  # return the naming for summary function if only asked for labels
  if(cuts_only){
    return(all_cuts)
  }

  # create mapping to weights
  omega_mapping <- list()
  for(p in seq_along(prior_list)){
    if(is.prior.weightfunction(prior_list[[p]])){
      expansion <- prior_expansions[[p]]
      omega_mapping[[p]] <- expansion$index[.weightfunction_global_bin_indices(all_cuts, expansion)]
    }
  }


  return(omega_mapping)
}

.weightfunction_mapping_info <- function(prior_list, one_sided = FALSE){

  cuts <- weightfunctions_mapping(prior_list, cuts_only = TRUE, one_sided = one_sided)
  list(
    mapping = weightfunctions_mapping(prior_list, one_sided = one_sided),
    cuts    = cuts,
    names   = .weightfunction_omega_names(cuts),
    pars    = paste0("omega[", seq_len(length(cuts) - 1L), "]"),
    one_sided = one_sided
  )
}

.weightfunction_set_omega_context <- function(samples, omega_context){

  if(is.null(omega_context)){
    return(samples)
  }

  attr(samples, "omega_context") <- omega_context
  prior_list <- attr(samples, "prior_list")
  if(!is.null(prior_list)){
    attr(prior_list, "omega_context") <- omega_context
    attr(samples, "prior_list") <- prior_list
  }

  samples
}

.weightfunction_omega_names <- function(cuts){
  sapply(seq_len(length(cuts) - 1L), function(i){
    paste0("omega[", cuts[i], ",", cuts[i + 1L], "]")
  })
}

.weightfunction_unique_cuts <- function(cuts, tolerance = sqrt(.Machine$double.eps)){

  cuts <- sort(cuts)
  if(length(cuts) <= 1L){
    return(cuts)
  }

  out <- cuts[1L]
  for(cut in cuts[-1L]){
    if(abs(cut - out[length(out)]) <= tolerance){
      if(abs(cut) <= tolerance || abs(out[length(out)]) <= tolerance){
        out[length(out)] <- 0
      }else if(abs(cut - 1) <= tolerance || abs(out[length(out)] - 1) <= tolerance){
        out[length(out)] <- 1
      }else{
        out[length(out)] <- min(cut, out[length(out)])
      }
    }else{
      out <- c(out, cut)
    }
  }

  out
}

.weightfunction_global_bin_indices <- function(global_cuts, expansion, tolerance = sqrt(.Machine$double.eps)){

  vapply(seq_len(length(global_cuts) - 1L), function(i){
    ind <- which(
      global_cuts[i] >= expansion$lower - tolerance &
        global_cuts[i + 1L] <= expansion$upper + tolerance
    )
    if(length(ind) != 1L){
      stop("Could not map global weightfunction bin to a local bin.", call. = FALSE)
    }
    ind
  }, integer(1))
}

.weightfunction_mapping_expansion <- function(prior, force_one_sided = FALSE){

  if(prior$side == "two-sided" && force_one_sided){
    J <- .weightfunction_n_bins(prior)
    cuts <- c(0, prior$steps / 2, 1 - rev(prior$steps) / 2, 1)
    index <- c(seq_len(J), seq.int(J - 1L, 1L))
  }else{
    cuts <- .weightfunction_local_cuts(prior)
    index <- seq_len(.weightfunction_n_bins(prior))
  }

  list(
    cuts  = cuts,
    lower = cuts[-length(cuts)],
    upper = cuts[-1],
    index = index
  )
}
