# Internal random-effect parameterization policy helpers.

.bt_random_effect_parameterization_requested <- function(block_prior){

  parameterization <- block_prior$parameterization
  if(is.null(parameterization)){
    parameterization <- "noncentered"
  }
  check_char(
    parameterization,
    "parameterization",
    allow_values = c("noncentered", "centered", "auto"),
    allow_NA = FALSE
  )

  parameterization
}

.bt_random_effect_prior_has_zero_atom <- function(prior){

  if(is.null(prior) || is.prior.none(prior)){
    return(FALSE)
  }
  if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    return(length(location) == 1L && !is.na(location) && location == 0)
  }
  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    return(any(vapply(prior, .bt_random_effect_prior_has_zero_atom, logical(1))))
  }
  if(is.prior.ordered(prior)){
    if(.bt_random_effect_prior_has_zero_atom(prior$total)){
      return(TRUE)
    }
    allocation <- prior$allocation
    return(
      is.list(allocation) && identical(allocation$type, "fixed") &&
        any(allocation$weights == 0)
    )
  }

  FALSE
}

.bt_random_effect_centered_eligibility <- function(block_prior, prior_list,
                                                   sd_binding,
                                                   row_indexed_external_sd){

  if(isTRUE(row_indexed_external_sd)){
    return(list(ok = FALSE, reason = "row-indexed external SD source"))
  }
  if(!is.null(block_prior$sd_source) ||
     .bt_random_sd_binding_has_external_source(sd_binding)){
    return(list(ok = FALSE, reason = "external SD source without a positive-support contract"))
  }
  binding_sources <- if(is.null(sd_binding)){
    list()
  }else{
    c(list(sd_binding$source), sd_binding$sources_by_column)
  }
  binding_sources <- binding_sources[!vapply(binding_sources, is.null, logical(1))]
  binding_priors <- lapply(binding_sources, function(source) source$prior)
  binding_priors <- binding_priors[!vapply(binding_priors, is.null, logical(1))]
  scale_priors <- c(prior_list, binding_priors)
  if(any(vapply(scale_priors, .bt_random_effect_prior_has_zero_atom, logical(1)))){
    return(list(ok = FALSE, reason = "SD prior with an atom at zero"))
  }
  if(!is.null(sd_binding)){
    factors <- c(
      sd_binding$factors,
      unlist(sd_binding$factors_by_column, recursive = FALSE)
    )
    factor_inclusion <- any(vapply(factors, function(factor){
      is.list(factor) && !is.null(factor$inclusion_name)
    }, logical(1)))
    allocations <- sd_binding$allocations
    allocation_inclusion <- any(vapply(allocations, function(allocation){
      is.list(allocation) && length(allocation$inclusion) > 0L
    }, logical(1)))
    if(factor_inclusion || allocation_inclusion){
      return(list(ok = FALSE, reason = "variance-allocation inclusion gate"))
    }
  }

  list(ok = TRUE, reason = "strictly nondegenerate random-effect scales")
}

.bt_random_effect_auto_centered_design <- function(model_matrix, group_map,
                                                   n_groups = max(group_map),
                                                   max_columns = 8L,
                                                   min_effective_n = 5,
                                                   min_rcond = 1e-4){

  K <- ncol(model_matrix)
  if(K > max_columns){
    return(list(ok = FALSE, reason = paste0("more than ", max_columns, " columns")))
  }

  observed_groups <- sort(unique(as.integer(group_map)))
  if(!identical(observed_groups, seq_len(as.integer(n_groups)))){
    return(list(ok = FALSE, reason = "one or more grouping levels are unobserved"))
  }

  for(group in observed_groups){
    X <- model_matrix[group_map == group, , drop = FALSE]
    column_scale <- apply(abs(X), 2L, max)
    X_scaled <- sweep(X, 2L, ifelse(column_scale > 0, column_scale, 1), "/")
    squared <- colSums(X_scaled^2)
    fourth  <- colSums(X_scaled^4)
    effective_n <- ifelse(fourth > 0, squared^2 / fourth, 0)
    if(any(!is.finite(effective_n)) || any(effective_n < min_effective_n)){
      return(list(ok = FALSE, reason = "insufficient within-group information"))
    }
    if(K > 1L){
      X_normalized <- sweep(X_scaled, 2L, sqrt(squared), "/")
      condition <- rcond(crossprod(X_normalized))
      if(!is.finite(condition) || condition < min_rcond){
        return(list(ok = FALSE, reason = "rank-deficient or ill-conditioned within-group design"))
      }
    }
  }

  list(
    ok = TRUE,
    reason = "within-group replication and conditioning favor centered parameterization"
  )
}

.bt_random_effect_resolve_parameterization <- function(block_prior, prior_list,
                                                       sd_binding,
                                                       row_indexed_external_sd,
                                                       model_matrix, group_map,
                                                       n_groups = max(group_map),
                                                       compile_mode,
                                                       block_name = NULL){

  requested <- .bt_random_effect_parameterization_requested(block_prior)
  if(!identical(compile_mode, "sampled")){
    return(list(
      requested = requested,
      resolved = "marginalized",
      reason = "random effect is analytically marginalized"
    ))
  }

  eligibility <- .bt_random_effect_centered_eligibility(
    block_prior = block_prior,
    prior_list = prior_list,
    sd_binding = sd_binding,
    row_indexed_external_sd = row_indexed_external_sd
  )
  if(identical(requested, "centered") && !isTRUE(eligibility$ok)){
    block_label <- if(is.null(block_name)) "unknown" else block_name
    stop(
      "Centered parameterization is not available for random-effect block '",
      block_label, "': ", eligibility$reason, ".",
      call. = FALSE
    )
  }
  if(identical(requested, "noncentered")){
    return(list(
      requested = requested,
      resolved = "noncentered",
      reason = "explicit noncentered parameterization"
    ))
  }
  if(identical(requested, "centered")){
    return(list(
      requested = requested,
      resolved = "centered",
      reason = "explicit centered parameterization"
    ))
  }
  if(!isTRUE(eligibility$ok)){
    return(list(
      requested = requested,
      resolved = "noncentered",
      reason = eligibility$reason
    ))
  }

  design <- .bt_random_effect_auto_centered_design(
    model_matrix = model_matrix,
    group_map = group_map,
    n_groups = n_groups
  )
  list(
    requested = requested,
    resolved = if(isTRUE(design$ok)) "centered" else "noncentered",
    reason = design$reason
  )
}
