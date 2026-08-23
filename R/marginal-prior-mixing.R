.mix_priors                <- function(prior_list, seed = NULL, n_samples = 10000){

  check_list(prior_list, "prior_list")
  for(i in seq_along(prior_list)){
    if(any(!sapply(prior_list[[i]], is.prior)))
      stop("'prior_list' must be a list of prior distributions")
  }
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  ### get model indices
  prior_weights <- do.call(cbind, lapply(seq_along(prior_list), function(i){
    prior_weights <- sapply(prior_list[[i]], .prior_model_weight)
    prior_weights <- prior_weights / sum(prior_weights)
    return(prior_weights)
  }))
  if(!all(prior_weights[,1] == prior_weights))
    stop("the prior samples are not alligned across models/draws")
  prior_weights <- prior_weights[,1]

  # set seed only once at the beginning -- not in the individual draws as the priors will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }
  sample_counts <- .prior_mixture_sample_counts(prior_weights, n_samples)

  ### adapted from 'mix_posteriors'
  parameters <- names(prior_list)
  out        <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter    <- parameters[p]
    temp_priors       <- prior_list[[temp_parameter]]

    if(any(sapply(temp_priors, is.prior.weightfunction)) && all(sapply(temp_priors, is.prior.weightfunction) | sapply(temp_priors, .is_prior_weightfunction_null) | sapply(temp_priors, is.null))){
      # weightfunctions:

      # replace missing priors with default prior: none
      for(i in 1:length(temp_priors)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior_none(prior_weights = prior_weights[i])
        }
      }

      out[[temp_parameter]] <- .mix_priors.weightfunction(
        temp_priors, temp_parameter, NULL, n_samples, sample_counts
      )

    }else if(any(sapply(temp_priors, is.prior.factor)) && all(sapply(temp_priors, is.prior.factor) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # factor priors

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(temp_priors)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0), prior_weights = prior_weights[i])
        }
      }

      out[[temp_parameter]] <- .mix_priors.factor(
        temp_priors, temp_parameter, NULL, n_samples, sample_counts
      )

    }else if(any(sapply(temp_priors, is.prior.vector)) && all(sapply(temp_priors, is.prior.vector) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # vector priors:

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(temp_priors)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0), prior_weights = prior_weights[i])
        }
      }

      out[[temp_parameter]] <- .mix_priors.vector(
        temp_priors, temp_parameter, NULL, n_samples, sample_counts
      )

    }else if(all(sapply(temp_priors, is.prior.simple) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # simple priors:

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(temp_priors)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0), prior_weights = prior_weights[i])
        }
      }

      out[[temp_parameter]] <- .mix_priors.simple(
        temp_priors, temp_parameter, NULL, n_samples, sample_counts
      )

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
.mix_priors.simple <- function(priors, parameter, seed = NULL, n_samples = 10000,
                               sample_counts = NULL){

  # check input
  check_list(priors, "priors")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(priors, is.prior.simple) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of simple priors")

  # get prior model probabilities
  prior_probs <- sapply(priors, .prior_model_weight)
  prior_probs <- prior_probs / sum(prior_probs)

  # do not set seed when sampling multiple priors for the same model -- they will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  samples <- NULL
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  if(is.null(sample_counts)){
    sample_counts <- .prior_mixture_sample_counts(prior_probs, n_samples)
  }
  for(i in seq_along(priors)[sample_counts > 0]){

    # sample indexes
    temp_ind <- seq_len(sample_counts[i])

    # sample prior
    samples <- c(samples, rng(priors[[i]], length(temp_ind), transform_factor_samples = FALSE))

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  # assure the correct number of samples
  samples    <- samples[1:n_samples]
  sample_ind <- sample_ind[1:n_samples]
  models_ind <- models_ind[1:n_samples]

  samples <- unname(samples)
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  samples <- .posterior_support_set_from_prior_list(samples, priors)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors, prior_probs,
      n_columns = 1L,
      column_names = parameter,
      source = "prior_model_probabilities"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.simple")

  return(samples)
}
.mix_priors.vector <- function(priors, parameter, seed = NULL, n_samples = 10000,
                               sample_counts = NULL){

  # check input
  check_list(priors, "priors")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(priors, is.prior.vector) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of vector priors")

  # get prior model probabilities
  prior_probs <- sapply(priors, .prior_model_weight)
  prior_probs <- prior_probs / sum(prior_probs)

  # do not set seed when sampling multiple priors for the same model -- they will end up completely correlated
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
  if(is.null(sample_counts)){
    sample_counts <- .prior_mixture_sample_counts(prior_probs, n_samples)
  }
  for(i in seq_along(priors)[sample_counts > 0]){

    # sample indexes
    temp_ind <- seq_len(sample_counts[i])

    if(is.prior.point(priors[[i]]) & is.prior.simple(priors[[i]])){
      # not sampling the priors in case they were imputed (missing dimensions)
      samples <- rbind(samples, matrix(priors[[i]]$parameters[["location"]], nrow = length(temp_ind), ncol = K))
    }else if(K == 1){
      samples <- rbind(samples, matrix(rng(priors[[i]], length(temp_ind), transform_factor_samples = FALSE), nrow = length(temp_ind), ncol = K))
    }else{
      samples <- rbind(samples, rng(priors[[i]], length(temp_ind), transform_factor_samples = FALSE))
    }

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  # assure the correct number of samples
  samples    <- samples[1:n_samples,,drop=FALSE]
  sample_ind <- sample_ind[1:n_samples]
  models_ind <- models_ind[1:n_samples]

  rownames(samples) <- NULL
  colnames(samples) <- paste0(parameter,"[",1:K,"]")
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  samples <- .posterior_support_set_columns_from_prior_list(samples, priors)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors, prior_probs,
      n_columns = K,
      column_names = colnames(samples),
      source = "prior_model_probabilities"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.vector")

  return(samples)
}
.mix_priors.factor <- function(priors, parameter, seed = NULL, n_samples = 10000,
                               sample_counts = NULL){

  # check input
  check_list(priors, "priors")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(priors, is.prior.factor) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of factor priors")

  for(i in seq_along(priors)){
    if(is.prior.ordered(priors[[i]])){
      priors[[i]] <- .prior_ordered_default_bound(priors[[i]], parameter)
    }
  }

  # get prior model probabilities
  prior_probs <- sapply(priors, .prior_model_weight)
  prior_probs <- prior_probs / sum(prior_probs)

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
        "interaction_terms" = attr(p, "interaction_terms"),
        "term_components"   = attr(p, "term_components"),
        "factor_terms"      = attr(p, "factor_terms"),
        "factor_contrasts"  = attr(p, "factor_contrasts"),
        "factor_design"     = attr(p, "factor_design"),
        "factor_cell_names" = attr(p, "factor_cell_names"),
        "treatment"   = is.prior.treatment(p),
        "independent" = is.prior.independent(p),
        "orthonormal" = is.prior.orthonormal(p),
        "meandif"     = is.prior.meandif(p),
        "ordered"     = is.prior.ordered(p)
      ))
    }else{
      stop("unsupported prior type")
    }
  })
  priors_info <- priors_info[!sapply(priors_info, isFALSE)]
  if(length(priors_info) >= 2 && any(!vapply(
    priors_info,
    function(i) isTRUE(all.equal(i, priors_info[[1]])),
    logical(1)
  ))){
    stop("non-matching prior factor type specifications")
  }
  priors_info <- priors_info[[1]]

  # Draw the model allocation once. All coefficients of a factor must use the
  # same model assignment so their rows remain joint draws.
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(sample_counts)){
    sample_counts <- .prior_mixture_sample_counts(prior_probs, n_samples)
  }

  if(priors_info[["ordered"]]){

    ordered_prior <- priors[[which(vapply(priors, is.prior.ordered, logical(1)))[1]]]
    coefficient_names <- .JAGS_prior_factor_names(parameter, ordered_prior)
    samples    <- matrix(nrow = 0, ncol = levels)
    sample_ind <- NULL
    models_ind <- NULL

    for(i in seq_along(priors)[sample_counts > 0]){

      temp_ind <- seq_len(sample_counts[i])

      if(is.prior.point(priors[[i]])){
        temp_samples <- matrix(
          priors[[i]]$parameters[["location"]],
          nrow = length(temp_ind),
          ncol = levels
        )
      }else{
        temp_samples <- rng(
          priors[[i]],
          length(temp_ind),
          transform_factor_samples = FALSE,
          quantity = "coefficient"
        )
        if(!is.matrix(temp_samples)){
          temp_samples <- matrix(temp_samples, nrow = length(temp_ind))
        }
        if(nrow(temp_samples) != length(temp_ind) || ncol(temp_samples) != levels){
          stop(
            "Ordered prior coefficient samples for '", parameter,
            "' do not match the bound factor design.",
            call. = FALSE
          )
        }
      }

      samples <- rbind(samples, temp_samples)
      sample_ind <- c(sample_ind, temp_ind)
      models_ind <- c(models_ind, rep(i, length(temp_ind)))
    }

    rownames(samples) <- NULL
    colnames(samples) <- coefficient_names
    attr(samples, "sample_ind") <- sample_ind
    attr(samples, "models_ind") <- models_ind
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- priors
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(priors_info[["treatment"]]){

    if(levels == 1){

      samples <- .mix_priors.simple(
        priors, parameter, NULL, n_samples, sample_counts
      )

      sample_ind <- attr(samples, "sample_ind")
      models_ind <- attr(samples, "models_ind")

      samples <- matrix(samples, ncol = 1)

    }else{

      samples <- lapply(1:levels, function(i) .mix_priors.simple(
        priors, paste0(parameter, "[", i, "]"), NULL, n_samples, sample_counts
      ))

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

      samples <- .mix_priors.simple(
        priors, parameter, NULL, n_samples, sample_counts
      )

      sample_ind <- attr(samples, "sample_ind")
      models_ind <- attr(samples, "models_ind")

      samples <- matrix(samples, ncol = 1)

    }else{

      samples <- lapply(1:levels, function(i) .mix_priors.simple(
        priors, paste0(parameter, "[", i, "]"), NULL, n_samples, sample_counts
      ))

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

    samples <- .mix_priors.vector(
      priors, parameter, NULL, n_samples, sample_counts
    )
    class(samples) <- c(class(samples), "mixed_posteriors.factor")

  }

  attr(samples, "levels")      <- priors_info[["levels"]]
  attr(samples, "level_names") <- priors_info[["level_names"]]
  attr(samples, "interaction") <- priors_info[["interaction"]]
  attr(samples, "interaction_terms") <- priors_info[["interaction_terms"]]
  attr(samples, "term_components")   <- priors_info[["term_components"]]
  attr(samples, "factor_terms")      <- priors_info[["factor_terms"]]
  attr(samples, "factor_contrasts")  <- priors_info[["factor_contrasts"]]
  attr(samples, "factor_design")     <- priors_info[["factor_design"]]
  attr(samples, "factor_cell_names") <- priors_info[["factor_cell_names"]]
  attr(samples, "treatment")   <- priors_info[["treatment"]]
  attr(samples, "independent") <- priors_info[["independent"]]
  attr(samples, "orthonormal") <- priors_info[["orthonormal"]]
  attr(samples, "meandif")     <- priors_info[["meandif"]]
  attr(samples, "ordered")     <- priors_info[["ordered"]]
  if(isTRUE(priors_info[["ordered"]])){
    ordered_prior <- priors[[which(vapply(priors, is.prior.ordered, logical(1)))[1]]]
    attr(samples, "ordered_metadata") <- attr(ordered_prior, "ordered_metadata")
  }

  if(isTRUE(priors_info[["treatment"]]) || isTRUE(priors_info[["independent"]])){
    factor_support <- .posterior_support_from_prior_list(priors)
    if(!is.null(factor_support) && !is.null(colnames(samples))){
      attr(samples, "posterior_support") <- stats::setNames(
        rep(list(factor_support), ncol(samples)),
        colnames(samples)
      )
    }
  }

  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors, prior_probs,
      n_columns = ncol(samples),
      column_names = colnames(samples),
      source = "prior_model_probabilities"
    )
  )

  return(samples)
}

.prior_mixture_sample_counts <- function(prior_probs, n_samples){

  .mixture_sample_counts(prior_probs, n_samples)
}

.mix_priors.weightfunction <- function(
    priors, parameter, seed = NULL, n_samples = 10000,
    sample_counts = NULL){

  # check input
  check_list(priors, "priors")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(priors, is.prior.weightfunction) | sapply(priors, .is_prior_weightfunction_null)))
    stop("'priors' must be a list of weightfunction priors or point(1)/none null priors")

  # get prior model probabilities
  prior_probs <- sapply(priors, .prior_model_weight)
  prior_probs <- prior_probs / sum(prior_probs)

  # do not set seed when sampling multiple priors for the same model -- they will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }

  # obtain mapping for the weight coefficients
  omega_info    <- .weightfunction_mapping_info(priors)
  omega_mapping <- omega_info$mapping
  omega_cuts    <- omega_info$cuts
  omega_names   <- omega_info$names

  # prepare output objects
  samples    <- matrix(nrow = 0, ncol = length(omega_cuts) - 1)
  sample_ind <- NULL
  models_ind <- NULL

  # mix samples
  if(is.null(sample_counts)){
    sample_counts <- .prior_mixture_sample_counts(prior_probs, n_samples)
  }
  for(i in seq_along(priors)[sample_counts > 0]){

    # sample indexes
    temp_ind <- seq_len(sample_counts[i])

    if(.is_prior_weightfunction_null(priors[[i]])){
      samples <- rbind(samples, matrix(1, ncol = length(omega_cuts) - 1, nrow = length(temp_ind)))
    }else{
      # create temp samples so names can be matched by mapping
      temp_samples <- rng(priors[[i]], length(temp_ind))
      colnames(temp_samples) <- paste0("omega[",1:ncol(temp_samples),"]")
      samples <- rbind(samples, temp_samples[, paste0("omega[",omega_mapping[[i]],"]")])
    }

    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  # assure the correct number of samples
  samples    <- samples[seq_len(n_samples),,drop=FALSE]
  sample_ind <- sample_ind[seq_len(n_samples)]
  models_ind <- models_ind[seq_len(n_samples)]

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  samples <- .weightfunction_set_omega_context(samples, omega_info)
  samples <- .posterior_support_set_weightfunction_columns(samples, priors, omega_info)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors, prior_probs,
      n_columns = ncol(samples),
      column_names = colnames(samples),
      source = "prior_model_probabilities",
      null_location = 1
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}

.as_mixed_priors            <- function(prior_list, seed = NULL, n_samples = 10000, conditional = NULL, conditional_rule = NULL){

  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of prior distributions")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(is.null(conditional_rule)){
    conditional_rule <- "AND"
  }
  condition_event <- .condition_event(
    prior_list        = prior_list,
    conditional       = conditional,
    conditional_rule  = conditional_rule
  )

  # set seed only once at the beginning -- not in the individual draws as the priors will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }

  # adapted from 'as_mixed_posteriors'
  parameters <- names(prior_list)
  out        <- list()

  # estimate the number of necessary samples for conditioning
  if(length(condition_event[["conditional"]]) > 0){

    condition_models <- .condition_event_model_options(prior_list, condition_event)
    conditioning_probability <- if(is.null(condition_models)){
      1
    }else{
      condition_models[["event_probability"]]
    }
    if(!is.finite(conditioning_probability) || conditioning_probability <= 0){
      stop("No prior models remain after applying the conditional event.", call. = FALSE)
    }
    all_alternative <- conditioning_probability == 1

    # multiply by 1.25 to ensure that the requested number of samples is reached
    requested_samples <- n_samples
    n_samples <- round(n_samples / conditioning_probability * 1.25)
  }

  # create the samples
  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter <- parameters[p]
    temp_prior     <- prior_list[[temp_parameter]]

    if(is.prior.spike_and_slab(temp_prior)){
      # spike and slab priors
      out[[temp_parameter]] <- .as_mixed_priors.spike_and_slab(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is.prior.mixture(temp_prior)){
      # mixture priors
      out[[temp_parameter]] <- .as_mixed_priors.mixture(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is.prior.weightfunction(temp_prior)){
      # weightfunctions:
      out[[temp_parameter]] <- .as_mixed_priors.weightfunction(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is_prior_phacking(temp_prior)){
      # p-hacking priors:
      out[[temp_parameter]] <- .as_mixed_priors.phacking(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is_prior_bias(temp_prior)){
      # composed publication-bias priors:
      out[[temp_parameter]] <- .as_mixed_priors.bias(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is.prior.factor(temp_prior)){
      # factor priors
      out[[temp_parameter]] <- .as_mixed_priors.factor(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is.prior.vector(temp_prior)){
      # vector priors:
      out[[temp_parameter]] <- .as_mixed_priors.vector(temp_prior, temp_parameter, NULL, n_samples)

    }else if(is.prior.simple(temp_prior)){
      # simple priors:
      out[[temp_parameter]] <- .as_mixed_priors.simple(temp_prior, temp_parameter, NULL, n_samples)

    }else{
      stop("The posterior samples cannot be mixed: unsupported mixture of prior distributions.")
    }

    # add formula relevant information
    if(!is.null(attr(temp_prior, which = "parameter"))){
      class(out[[temp_parameter]]) <- c(class(out[[temp_parameter]]), "mixed_posteriors.formula")
      attr(out[[temp_parameter]], "formula_parameter")  <- attr(temp_prior, which = "parameter")
    }
  }


  # perform conditioning (and copy back with attributes)
  if(length(condition_event[["conditional"]]) > 0){

    # obtain the indicator samples
    indicator_samples <- matrix(nrow = n_samples, ncol = 0L)
    for(parameter in parameters){
      models_ind <- attr(out[[parameter]], "models_ind")
      if(length(models_ind) == n_samples && !identical(models_ind, rep(FALSE, n_samples))){
        indicator_samples <- cbind(indicator_samples, models_ind)
        colnames(indicator_samples)[ncol(indicator_samples)] <- paste0(parameter, "_indicator")
        if(parameter == "bias"){
          indicator_samples <- cbind(indicator_samples, models_ind)
          colnames(indicator_samples)[ncol(indicator_samples)] <- "bias_indicator"
        }
      }
    }
    conditioning_samples <- .condition_event_posterior_mask(
      event         = condition_event,
      prior_list    = prior_list,
      model_samples = indicator_samples
    )

    # check enough samples were drawn (if too many remove the extra ones)
    if(sum(conditioning_samples) < requested_samples){
      warning(sprintf("Only %d samples were drawn from the prior distributions due to conditioning.", sum(conditioning_samples)))
    }else{
      conditioning_samples[which(conditioning_samples)[-(1:requested_samples)]] <- FALSE
    }

    # select the conditional samples (and copy attributes)
    for(p in seq_along(parameters)){
      temp <- attributes(out[[parameters[p]]])
      if(is.null(dim(out[[parameters[p]]]))){
        out[[parameters[p]]] <- out[[parameters[p]]][conditioning_samples]
        attributes(out[[parameters[p]]]) <- c(attributes(out[[parameters[p]]]), temp)
        attr(out[[parameters[p]]], "models_ind") <- attr(out[[parameters[p]]], "models_ind")[conditioning_samples]
      }else{
        out[[parameters[p]]] <- out[[parameters[p]]][conditioning_samples,,drop=FALSE]
        attributes(out[[parameters[p]]]) <- c(attributes(out[[parameters[p]]])[!names(attributes(out[[parameters[p]]])) %in% c("dimnames")], temp[!names(temp) %in% c("dim")])
        attr(out[[parameters[p]]], "models_ind") <- attr(out[[parameters[p]]], "models_ind")[conditioning_samples]
      }
      out[[parameters[p]]] <- .posterior_atoms_refresh_from_prior(
        out[[parameters[p]]],
        prior_list[[parameters[p]]],
        parameters[p]
      )
    }

    # put a check whether all samples were conditional
    attr(out, "all_alternative") <- all_alternative
  }

  for(parameter in names(out)){
    out[[parameter]] <- .condition_event_set_attributes(
      out[[parameter]],
      condition_event
    )
  }
  out <- .condition_event_set_attributes(out, condition_event)

  return(out)
}
.as_mixed_priors.simple         <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  # do not set seed when sampling multiple prior for the same model -- they will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  samples <- rng(prior, n_samples, transform_factor_samples = FALSE)

  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  samples <- .posterior_support_set_from_prior_list(samples, prior)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      prior, 1,
      n_columns = 1L,
      column_names = parameter,
      source = "single_prior_structure"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.simple")

  return(samples)
}
.as_mixed_priors.vector         <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  # do not set seed when sampling multiple prior for the same model -- they will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  K <- prior$parameters[["K"]]

  if(is.prior.point(prior) & is.prior.simple(prior)){
    # not sampling the prior in case they were imputed (missing dimensions)
    samples <- matrix(prior$parameters[["location"]], nrow = n_samples, ncol = K)
  }else if(K == 1){
    samples <- matrix(rng(prior, n_samples, transform_factor_samples = FALSE), nrow = n_samples, ncol = K)
  }else{
    samples <- rng(prior, n_samples, transform_factor_samples = FALSE)
  }

  rownames(samples) <- NULL
  colnames(samples) <- paste0(parameter,"[",1:K,"]")
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  samples <- .posterior_support_set_columns_from_prior_list(samples, prior)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      prior, 1,
      n_columns = K,
      column_names = colnames(samples),
      source = "single_prior_structure"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.vector")

  return(samples)
}
.as_mixed_priors.factor         <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  if(is.prior.ordered(prior)){
    prior <- .prior_ordered_default_bound(prior, parameter)
  }

  # check the prior levels
  levels <- .get_prior_factor_levels(prior)

  # gather and check compatibility of prior distributions
  prior_info <- list(
    "levels"      = .get_prior_factor_levels(prior),
    "level_names" = .get_prior_factor_level_names(prior),
    "interaction" = .is_prior_interaction(prior),
    "interaction_terms" = attr(prior, "interaction_terms"),
    "term_components"   = attr(prior, "term_components"),
    "factor_terms"      = attr(prior, "factor_terms"),
    "factor_contrasts"  = attr(prior, "factor_contrasts"),
    "factor_design"     = attr(prior, "factor_design"),
    "factor_cell_names" = attr(prior, "factor_cell_names"),
    "treatment"   = is.prior.treatment(prior),
    "independent" = is.prior.independent(prior),
    "orthonormal" = is.prior.orthonormal(prior),
    "meandif"     = is.prior.meandif(prior),
    "ordered"     = is.prior.ordered(prior)
  )

  if(prior_info[["ordered"]]){

    if(!is.null(seed)){
      set.seed(seed)
    }

    samples <- rng(
      prior,
      n_samples,
      transform_factor_samples = FALSE,
      quantity = "coefficient"
    )
    ordered_total_indicator <- attr(
      samples,
      "ordered_total_indicator",
      exact = TRUE
    )
    if(!is.matrix(samples)){
      samples <- matrix(samples, nrow = n_samples)
    }
    if(nrow(samples) != n_samples || ncol(samples) != levels){
      stop(
        "Ordered prior coefficient samples for '", parameter,
        "' do not match the bound factor design.",
        call. = FALSE
      )
    }

    rownames(samples) <- NULL
    colnames(samples) <- .JAGS_prior_factor_names(parameter, prior)
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- FALSE
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    if(!is.null(ordered_total_indicator)){
      attr(samples, "ordered_total_indicator") <- ordered_total_indicator
    }
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["treatment"]]){

    if(levels == 1){

      samples <- .as_mixed_priors.simple(prior, parameter, seed, n_samples)
      samples <- matrix(samples, ncol = 1)

    }else{

      if(!is.null(seed)){
        set.seed(seed)
      }

      samples <- lapply(1:levels, function(i) .as_mixed_priors.simple(prior, paste0(parameter, "[", i, "]"), NULL, n_samples))
      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- paste0(parameter,"[",prior_info$level_names[-1],"]")
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- FALSE
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["independent"]]){

    if(levels == 1){

      samples <- .as_mixed_priors.simple(prior, parameter, seed, n_samples)
      samples <- matrix(samples, ncol = 1)

    }else{

      if(!is.null(seed)){
        set.seed(seed)
      }

      samples <- lapply(1:levels, function(i) .as_mixed_priors.simple(prior, paste0(parameter, "[", i, "]"), NULL, n_samples))
      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- paste0(parameter,"[",prior_info$level_names,"]")
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- FALSE
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["orthonormal"]] | prior_info[["meandif"]]){

    prior$parameters[["K"]] <- levels
    samples <- .as_mixed_priors.vector(prior, parameter, seed, n_samples)
    class(samples) <- c(class(samples), "mixed_posteriors.factor")

  }

  attr(samples, "levels")      <- prior_info[["levels"]]
  attr(samples, "level_names") <- prior_info[["level_names"]]
  attr(samples, "interaction") <- prior_info[["interaction"]]
  attr(samples, "interaction_terms") <- prior_info[["interaction_terms"]]
  attr(samples, "term_components")   <- prior_info[["term_components"]]
  attr(samples, "factor_terms")      <- prior_info[["factor_terms"]]
  attr(samples, "factor_contrasts")  <- prior_info[["factor_contrasts"]]
  attr(samples, "factor_design")     <- prior_info[["factor_design"]]
  attr(samples, "factor_cell_names") <- prior_info[["factor_cell_names"]]
  attr(samples, "treatment")   <- prior_info[["treatment"]]
  attr(samples, "independent") <- prior_info[["independent"]]
  attr(samples, "orthonormal") <- prior_info[["orthonormal"]]
  attr(samples, "meandif")     <- prior_info[["meandif"]]
  attr(samples, "ordered")     <- prior_info[["ordered"]]
  attr(samples, "ordered_metadata") <- attr(prior, "ordered_metadata")

  if(isTRUE(prior_info[["treatment"]]) || isTRUE(prior_info[["independent"]])){
    factor_support <- .posterior_support_from_prior_list(prior)
    if(!is.null(factor_support) && !is.null(colnames(samples))){
      attr(samples, "posterior_support") <- stats::setNames(
        rep(list(factor_support), ncol(samples)),
        colnames(samples)
      )
    }
  }

  ordered_atoms <- .posterior_atoms_from_ordered_total(
    prior,
    n_columns = ncol(samples),
    column_names = colnames(samples),
    source = "ordered_total_prior"
  )
  samples <- .posterior_atoms_set(
    samples,
    if(is.null(ordered_atoms)){
      .posterior_atoms_from_priors(
        prior, 1,
        n_columns = ncol(samples),
        column_names = colnames(samples),
        source = "single_prior_structure"
      )
    }else{
      ordered_atoms
    }
  )

  return(samples)
}
.as_mixed_priors.weightfunction <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  # do not set seed when sampling multiple prior for the same model -- they will end up completely correlated
  if(!is.null(seed)){
    set.seed(seed)
  }

  # obtain mapping for the weight coefficients
  omega_cuts    <- weightfunctions_mapping(list(prior), cuts_only = TRUE)
  omega_names   <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))

  # prepare output objects
  samples <- rng(prior, n_samples)

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  omega_info <- .weightfunction_mapping_info(list(prior))
  samples <- .weightfunction_set_omega_context(samples, omega_info)
  samples <- .posterior_support_set_weightfunction_columns(samples, prior, omega_info)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      prior, 1,
      n_columns = ncol(samples),
      column_names = colnames(samples),
      source = "single_prior_structure",
      null_location = 1
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}
.as_mixed_priors.phacking <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  if(!is.null(seed)){
    set.seed(seed)
  }

  samples <- rng(prior, n_samples)
  par_names <- intersect(.phacking_report_parameter(prior), colnames(samples))
  samples <- samples[, par_names, drop = FALSE]

  rownames(samples) <- NULL
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.phacking")

  return(samples)
}
.as_mixed_priors.bias <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  if(!is.null(seed)){
    set.seed(seed)
  }

  spec          <- selection_backend_spec(prior)
  branch_info   <- .selection_prior_branch_info(prior)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

  samples <- rng(prior, n_samples)

  out_names <- character()
  par_names <- character()

  if(any(has_selection)){
    omega_cuts  <- spec$step$breaks
    omega_names <- sapply(seq_len(length(omega_cuts) - 1L), function(i) paste0("omega[", omega_cuts[i], ",", omega_cuts[i + 1L], "]"))
    omega_par   <- paste0("omega[", seq_len(length(omega_cuts) - 1L), "]")
    out_names   <- c(out_names, omega_names)
    par_names   <- c(par_names, omega_par)
  }
  if(any(has_phacking)){
    phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
    phacking_names <- .selection_phacking_report_parameters(phacking_priors)
    out_names <- c(out_names, phacking_names)
    par_names <- c(par_names, phacking_names)
  }

  keep <- par_names %in% colnames(samples)
  samples <- samples[, par_names[keep], drop = FALSE]
  colnames(samples) <- out_names[keep]

  rownames(samples) <- NULL
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- FALSE
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.bias")

  return(samples)
}
.as_mixed_priors.spike_and_slab <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  # set the optional seed once and let nested samplers continue the same stream
  if(!is.null(seed)){
    set.seed(seed)
  }

  prior_variable   <- .get_spike_and_slab_variable(prior)
  prior_inclusion  <- .get_spike_and_slab_inclusion(prior)

  inclusion <- stats::rbinom(n_samples, size = 1, prob = rng(prior_inclusion, n_samples))

  if(is.prior.factor(prior_variable)){

    samples <- .as_mixed_priors.factor(prior_variable, parameter, NULL, n_samples)

  }else if(is.prior.simple(prior_variable)){

    samples <- .as_mixed_priors.simple(prior_variable, parameter, NULL, n_samples)

  }

  # merge with names and attributes
  samples       <- samples * inclusion

  class(samples) <- c(class(samples), "mixed_posteriors.spike_and_slab")
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- inclusion
  attr(samples, "prior_list") <- prior
  if(!is.null(dim(samples))){
    samples <- .posterior_support_set_columns_from_prior_list(samples, prior)
  }else{
    samples <- .posterior_support_set_from_prior_list(samples, prior)
  }
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_indicator(
      prior = prior,
      indicator = inclusion,
      n_columns = if(is.null(dim(samples))) 1L else ncol(samples),
      column_names = if(is.null(dim(samples))) parameter else colnames(samples),
      spike_and_slab = TRUE
    )
  )

  return(samples)
}
.as_mixed_priors.mixture        <- function(prior, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(prior, "prior")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  is_PET            <- sapply(prior, is.prior.PET)
  is_PEESE          <- sapply(prior, is.prior.PEESE)
  is_weightfunction <- sapply(prior, is.prior.weightfunction)
  is_phacking       <- sapply(prior, is_prior_phacking)
  is_bias           <- sapply(prior, is_prior_bias)

  if(any(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias)){

    temp_samples <- .mix_priors.bias(prior, parameter = parameter, seed = seed, n_samples = n_samples)

  }else{

    is_factor <- sapply(prior, is.prior.factor)

    if(any(is_factor)){
      temp_samples <- .mix_priors.factor(prior, parameter = parameter, seed = seed, n_samples = n_samples)
    }else{
      temp_samples <- .mix_priors.simple(prior, parameter = parameter, seed = seed, n_samples = n_samples)
    }

  }

  # the samples parameters need to be randomly shuffled
  # (the  .mix_priors.XXX functions generate the samples model by model to keep bridge-sampling model-averaging consistent structure,
  #  this however does not apply to the spike an slab priors)
  random_ind <- sample(n_samples)
  if(is.null(dim(temp_samples))){
    samples <- temp_samples[random_ind]
  }else{
    samples <- temp_samples[random_ind,,drop=FALSE]
  }
  attributes(samples) <- attributes(temp_samples)
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- attr(samples, "models_ind")[random_ind]

  # append classes and priors
  class(samples) <- c(class(samples), "mixed_posteriors.mixture")
  attr(samples, "prior_list") <- prior
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_indicator(
      prior = prior,
      indicator = attr(samples, "models_ind"),
      n_columns = if(is.null(dim(samples))) 1L else ncol(samples),
      column_names = if(is.null(dim(samples))) parameter else colnames(samples)
    )
  )

  return(samples)
}
.mix_priors.bias <- function(priors, parameter, seed = NULL, n_samples = 10000){

  # check input
  check_list(priors, "priors")
  check_char(parameter, "parameter")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")

  allowed <- sapply(priors, function(prior){
    is.prior.none(prior) || is.prior.PET(prior) || is.prior.PEESE(prior) ||
      is.prior.weightfunction(prior) || is_prior_phacking(prior) || is_prior_bias(prior)
  })
  if(!all(allowed)){
    stop("'priors' must be a list of publication-bias priors.", call. = FALSE)
  }

  spec        <- selection_backend_spec(priors)
  branch_info <- lapply(priors, .selection_branch_info)

  is_PET       <- sapply(priors, is.prior.PET)
  is_PEESE     <- sapply(priors, is.prior.PEESE)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

  prior_probs <- sapply(priors, .prior_model_weight)
  prior_probs <- prior_probs / sum(prior_probs)

  if(!is.null(seed)){
    set.seed(seed)
  }

  out_names <- character()
  if(any(has_selection)){
    omega_cuts  <- spec$step$breaks
    omega_names <- .weightfunction_omega_names(omega_cuts)
    out_names   <- c(out_names, omega_names)

    selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
    omega_mapping    <- .weightfunction_mapping_info(selection_priors, one_sided = TRUE)$mapping
    selection_index  <- which(has_selection)
  }
  if(any(has_phacking)){
    phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
    phacking_names  <- .selection_phacking_report_parameters(phacking_priors)
    out_names <- c(out_names, phacking_names)
  }
  if(any(is_PET)){
    out_names <- c(out_names, "PET")
  }
  if(any(is_PEESE)){
    out_names <- c(out_names, "PEESE")
  }

  samples    <- matrix(nrow = 0, ncol = length(out_names))
  colnames(samples) <- out_names
  sample_ind <- NULL
  models_ind <- NULL

  sample_counts <- .prior_mixture_sample_counts(prior_probs, n_samples)
  for(i in seq_along(priors)[sample_counts > 0]){

    temp_ind <- seq_len(sample_counts[i])
    temp_samples <- matrix(0, nrow = length(temp_ind), ncol = length(out_names))
    colnames(temp_samples) <- out_names

    if(any(has_selection)){
      temp_samples[, omega_names] <- 1
    }

    if(has_selection[i]){
      selection_i <- match(i, selection_index)
      selection_samples <- rng(branch_info[[i]]$selection, length(temp_ind))
      temp_samples[, omega_names] <- selection_samples[, paste0("omega[", omega_mapping[[selection_i]], "]"), drop = FALSE]
    }

    if(has_phacking[i]){
      phacking_samples <- rng(branch_info[[i]]$phacking, length(temp_ind))
      phacking_names <- .phacking_report_parameter(branch_info[[i]]$phacking)
      temp_samples[, phacking_names] <- phacking_samples[, phacking_names, drop = FALSE]
    }

    if(is_PET[i]){
      temp_samples[, "PET"] <- rng(priors[[i]], length(temp_ind), transform_factor_samples = FALSE)
    }
    if(is_PEESE[i]){
      temp_samples[, "PEESE"] <- rng(priors[[i]], length(temp_ind), transform_factor_samples = FALSE)
    }

    samples    <- rbind(samples, temp_samples)
    sample_ind <- c(sample_ind, temp_ind)
    models_ind <- c(models_ind, rep(i, length(temp_ind)))
  }

  samples    <- samples[seq_len(n_samples), , drop = FALSE]
  sample_ind <- sample_ind[seq_len(n_samples)]
  models_ind <- models_ind[seq_len(n_samples)]

  rownames(samples) <- NULL
  attr(samples, "sample_ind") <- sample_ind
  attr(samples, "models_ind") <- models_ind
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  if(any(has_selection)){
    samples <- .weightfunction_set_omega_context(samples, list(
      mapping   = NULL,
      cuts      = omega_cuts,
      names     = omega_names,
      pars      = paste0("omega[", seq_len(length(omega_cuts) - 1L), "]"),
      one_sided = TRUE
    ))
  }
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.bias")

  return(samples)
}
