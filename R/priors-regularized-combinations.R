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
# Ridge uses the conditional Normal route. Two Laplace terms have an
# analytic conditional convolution. Wider shared-lasso grids have no
# verified display resolution.
.bt_regularized_linear_plan <- function(prior_list, weights, source_transforms = NULL){

  if(!length(weights) || !any(vapply(prior_list, .bt_has_regularized_term, logical(1)))) return(NULL)
  entries <- .prior_linear_active_entries(prior_list, weights)
  centered <- .bt_regularized_centered_plan(prior_list, entries, weights, source_transforms)
  if(!is.null(centered)) return(centered)
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
      is.null(prior$regularization$group_lasso_dimension) &&
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
  # In fitted coefficient units, the sum has the distribution of the
  # larger absolute weight times (L1 + r*L2), r the smaller-to-larger
  # ratio. Both independent Laplace directions use the same global scale.
  if(length(weights) == 2L && length(shared) == 1L &&
     all(vapply(prior_list[affected], function(p) p$regularization$family == "lasso", logical(1)))){
    name <- affected[[1L]]
    projection <- prior_list[[name]]
    width <- unlist(lapply(affected, function(owner){
      .prior_linear_fold_scale(abs(weights[coordinates[[owner]]]), prior_list[[owner]]$regularization$multiplier)
    }), use.names = FALSE)
    projection$regularization$contrast_ratio <- min(width) / max(width)
    if(projection$regularization$contrast_ratio == 0) .regularized_numerical_unavailable("contrast scale ratio")
    projection$regularization$dimension <- 1L
    projection$regularization$multiplier <- 1
    projection$regularization$unit_index <- NULL
    projection$regularization$unit_centered <- NULL
    projection$regularization$unit_dimension <- NULL
    projection$regularization$variance_names <- NULL
    class(projection) <- c("prior", "prior.regularized_term", "prior.simple")
    attributes(projection) <- attributes(projection)[c("names", "class")]
    return(list(prior_list = stats::setNames(list(projection), name),
                weights = stats::setNames(max(width), name)))
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
    if(spec$family == "lasso") attr(prior, "regularized_lasso_direction") <- TRUE
    prior_list[[name]] <- prior
    prior_list[[spec$global_name]] <- spec$scale
  }
  list(prior_list = prior_list)
}
