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
#' corresponding to the null hypothesis; use \code{0} or \code{integer(0)}
#' when no models are null)
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


  is_null <- .model_averaging_is_null(is_null, length(prior_probs))

  if(!missing(prior_probs) && !missing(margliks)){
    return(.inclusion_BF.margliks(prior_probs = prior_probs, margliks = margliks, is_null = is_null))
  }else if(!missing(prior_probs) && !missing(post_probs)){
    return(.inclusion_BF.probs(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null))
  }else{
    stop("'prior_probs' and either 'post_probs' or 'marglik' must be specified.")
  }
}

.model_averaging_is_null <- function(is_null, n_models){

  if(is.null(is_null)){
    return(rep(FALSE, n_models))
  }

  if(is.numeric(is_null)){
    check_int(is_null, "is_null", lower = 0, upper = n_models, check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
    if(length(is_null) == 0L || (length(is_null) == 1L && is_null == 0)){
      return(rep(FALSE, n_models))
    }
    if(any(is_null == 0)){
      stop("'is_null' can contain 0 only when no null models are specified.", call. = FALSE)
    }
    return(seq_len(n_models) %in% is_null)
  }

  if(is.logical(is_null)){
    check_bool(is_null, "is_null", check_length = n_models, allow_NA = FALSE)
    return(is_null)
  }

  stop("'is_null' argument must be either logical vector, integer vector, or NULL.", call. = FALSE)
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
