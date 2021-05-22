#' @title Compute posterior probabilities and inclusion Bayes factors
#'
#' @description Computes prior probabilities, posterior probabilities,
#' and inclusion Bayes factors based on prior model odds, marginal
#' likelihoods, and indicator whether the models represent the null or
#' alternative hypothesis
#'
#' @param prior_odds vector of prior model odds
#' @param margliks vector of marginal likelihoods
#' @param is_null logical vector of indicators whether the model corresponds
#' to the null or alternative hypothesis (or an integer vector indexing models
#' corresponding to the null hypothesis)
#' @param conditional whether prior and posterior model probabilities should
#' be returned only for the conditional model. Defaults to \code{FALSE}
#'
#' @export
models_inference <- function(prior_odds, margliks, is_null = NULL, conditional = FALSE){

  check_real(prior_odds, "prior_odds", lower = 0, check_length = 0)
  check_real(margliks,   "margliks", check_length = length(prior_odds))
  check_bool(conditional, "conditional")
  if(!is.null(is_null)){
    if(is.numeric(is_null)){
      check_int(is_null, "is_null", lower = 1, upper = length(prior_odds), check_length = 0)
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
