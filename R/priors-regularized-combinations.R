# The stochastic scale coordinates used by a coefficient, including the
# alternatives of a gated coefficient that have positive prior probability.
.bt_regularized_shared_scales <- function(prior){

  if(inherits(prior, "prior.regularized_term")){
    spec <- prior$regularization
    return(c(spec$global_name, spec$slab_name))
  }
  if(is.prior.mixture(prior)){
    probabilities <- .prior_density_ordinate_mixture_weights(prior)
    if(!is.null(probabilities)) prior <- prior[probabilities > 0]
    return(unique(unlist(lapply(prior, .bt_regularized_shared_scales), use.names = FALSE)))
  }
  character()
}

# A common scale makes even uncorrelated signed coefficients dependent.
# For ridge and lasso, condition on it by factoring it out of the sum:
# sum(w_j beta_j) = tau * sum(w_j z_j) with independent normal/Laplace z_j.
# The conditional Normal route supplies ridge ordinates. Lasso uses the
# existing product grid for display; a general product has no certified
# structural ordinate.
.bt_regularized_linear_plan <- function(prior_list, weights, source_transforms = NULL){

  if(!length(weights) || !any(vapply(prior_list, .bt_has_regularized_term, logical(1)))) return(NULL)
  entries <- .prior_linear_active_entries(prior_list, weights)
  coordinates <- list()
  multipliers <- character()
  for(name in names(entries)){
    entry <- entries[[name]]
    multiplier <- attr(entry$prior, "multiply_by", exact = TRUE)
    if(is.numeric(multiplier) && isTRUE(multiplier == 0)) next
    if(is.character(multiplier)){
      law <- prior_list[[multiplier]]
      if(is.prior.none(law) || (is.prior.point(law) && isTRUE(law$parameters$location == 0))) next
      multipliers <- c(multipliers, multiplier)
    }
    coordinates[[name]] <- entry$present
  }
  for(name in intersect(unique(multipliers), names(prior_list))){
    coordinates[[name]] <- union(coordinates[[name]], name)
  }
  owners <- names(coordinates)
  scales <- lapply(prior_list[owners], .bt_regularized_shared_scales)
  for(name in intersect(owners, unlist(scales, use.names = FALSE))){
    scales[[name]] <- unique(c(scales[[name]], name))
  }
  keys <- unlist(rep(scales, lengths(coordinates)), use.names = FALSE)
  shared <- unique(keys[duplicated(keys)])
  if(!length(shared)) return(NULL)
  affected <- owners[vapply(scales, function(keys) any(keys %in% shared), logical(1))]
  factorable <- vapply(affected, function(name){
    prior <- prior_list[[name]]
    inherits(prior, "prior.regularized_term") &&
      prior$regularization$family %in% c("ridge", "lasso") &&
      is.null(attr(prior, "multiply_by", exact = TRUE)) && !name %in% multipliers &&
      !any(!is.na(source_transforms[coordinates[[name]]]))
  }, logical(1))
  if(!all(factorable)){
    message <- paste0("The prior density of this combination is unavailable because its coefficients ",
      "share shrinkage scales (", paste(shared, collapse = ", "), "). Use joint prior draws.")
    # A nonzero coefficient of an ungated term has full real support and
    # a continuous density conditional on all scales and mixture gates.
    # A JAGS expression or a coefficient-valued multiplier could use that
    # same coefficient a second time and cancel this contribution.
    atom_free <- !length(multipliers) && !any(!is.na(source_transforms)) &&
      any(vapply(prior_list[affected], inherits, logical(1), "prior.regularized_term")) &&
      !any(vapply(prior_list[owners], .bt_prior_has_expression, logical(1)))
    return(list(reason = message, sources = affected, atom_free = atom_free))
  }
  for(name in affected){
    prior <- prior_list[[name]]
    spec <- prior$regularization
    base <- if(identical(spec$family, "ridge")) prior("normal", list(0, spec$multiplier)) else
      prior("laplace", list(0, spec$multiplier))
    prior$distribution <- base$distribution
    prior$parameters <- base$parameters
    prior$regularization <- NULL
    class(prior) <- setdiff(class(prior), "prior.regularized_term")
    attr(prior, "multiply_by") <- spec$global_name
    prior_list[[name]] <- prior
    prior_list[[spec$global_name]] <- spec$scale
  }
  list(prior_list = prior_list)
}
