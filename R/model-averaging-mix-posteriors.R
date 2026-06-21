#' @title Model-average posterior distributions
#'
#' @description Model-averages posterior distributions based on
#' a list of models, vector of parameters, and a list of
#' indicators of the null or alternative hypothesis models
#' for each parameter.
#'
#' @details For simplex parameters such as Dirichlet priors, absent model
#' parameters are not filled with an implicit spike at zero because the all-zero
#' vector is off the simplex. Every model with positive prior or posterior
#' probability must provide an explicit compatible simplex prior or point prior
#' on the simplex.
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

    if(all(sapply(temp_priors, is.null))){
      stop(
        "The parameter '", temp_parameter,
        "' is not available in any model prior list.",
        call. = FALSE
      )
    }

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

      .mix_posteriors_validate_simplex_priors(
        temp_priors,
        temp_parameter,
        temp_inference$prior_probs,
        temp_inference$post_probs
      )

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

.mix_posteriors_is_dirichlet_simplex <- function(prior){
  is.prior.simplex(prior) && identical(prior[["distribution"]], "dirichlet")
}

.mix_posteriors_is_simplex_point <- function(prior, K){

  if(!is.prior.point(prior)){
    return(FALSE)
  }

  location <- prior$parameters[["location"]]
  if(is.numeric(location) && length(location) == 1L){
    location <- rep(location, K)
  }
  is.numeric(location) &&
    length(location) == K &&
    all(is.finite(location)) &&
    all(location >= 0) &&
    isTRUE(all.equal(sum(location), 1, tolerance = 1e-8))
}

.mix_posteriors_validate_simplex_priors <- function(priors, parameter, prior_probs, post_probs){

  simplex <- vapply(priors, .mix_posteriors_is_dirichlet_simplex, logical(1))
  if(!any(simplex)){
    return(invisible(NULL))
  }

  K <- unique(vapply(priors[simplex], function(prior){
    prior$parameters[["K"]]
  }, numeric(1)))
  if(length(K) != 1L){
    stop(
      "The simplex parameter '", parameter,
      "' cannot be mixed across Dirichlet priors with different dimensions.",
      call. = FALSE
    )
  }

  contributing <- prior_probs > 0 | post_probs > 0
  missing <- vapply(priors, is.null, logical(1))
  if(any(missing & contributing)){
    stop(
      "The simplex parameter '", parameter,
      "' has no implicit spike-at-zero null. Every model with positive prior ",
      "or posterior probability must provide an explicit compatible Dirichlet ",
      "simplex prior or point prior on the simplex.",
      call. = FALSE
    )
  }

  compatible <- vapply(priors, function(prior){
    if(is.null(prior)){
      return(TRUE)
    }
    .mix_posteriors_is_dirichlet_simplex(prior) ||
      .mix_posteriors_is_simplex_point(prior, K)
  }, logical(1))

  if(any(!compatible & contributing)){
    stop(
      "The simplex parameter '", parameter,
      "' can only mix Dirichlet simplex priors with the same dimension or ",
      "explicit point priors on the simplex.",
      call. = FALSE
    )
  }

  invisible(NULL)
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
  samples <- .posterior_support_set_from_prior_list(samples, priors)
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
      samples <- rbind(
        samples,
        matrix(
          rep(priors[[i]]$parameters[["location"]], times = length(temp_ind)),
          nrow = length(temp_ind),
          ncol = K,
          byrow = TRUE
        )
      )
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
  samples <- .posterior_support_set_columns_from_prior_list(samples, priors)
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
  priors <- .complete_factor_metadata_prior_list(priors)

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

  if(isTRUE(priors_info[["treatment"]]) || isTRUE(priors_info[["independent"]])){
    factor_support <- .posterior_support_from_prior_list(priors)
    if(!is.null(factor_support) && !is.null(colnames(samples))){
      attr(samples, "posterior_support") <- stats::setNames(
        rep(list(factor_support), ncol(samples)),
        colnames(samples)
      )
    }
  }

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
  samples <- .posterior_support_set_weightfunction_columns(
    samples,
    priors,
    .weightfunction_mapping_info(priors)
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}

