.bt_random_effect_latent_names <- function(random_term, n_groups, n_columns){

  check_int(n_groups, "n_groups", lower = 0, allow_NA = FALSE)
  check_int(n_columns, "n_columns", lower = 0, allow_NA = FALSE)
  layout <- random_term$latent_layout
  if(inherits(layout, "BayesTools_random_effect_structured_local_layout")){
    node_names <- layout$node_names
    if(!is.character(node_names) ||
       length(node_names) != layout$n_local ||
       anyNA(node_names) ||
       any(!nzchar(node_names)) ||
       anyDuplicated(node_names)){
      stop(
        "Random-effect local latent metadata",
        .bt_random_effect_metadata_block_detail(random_term),
        " must contain one unique, non-missing, non-empty character ",
        "node name per local latent cell.",
        call. = FALSE
      )
    }
    local_group <- layout$local_group
    local_column <- layout$local_column
    valid_local_indices <- is.numeric(local_group) &&
      is.numeric(local_column) &&
      length(local_group) == layout$n_local &&
      length(local_column) == layout$n_local &&
      !anyNA(local_group) &&
      !anyNA(local_column) &&
      all(is.finite(local_group)) &&
      all(is.finite(local_column)) &&
      all(local_group == as.integer(local_group)) &&
      all(local_column == as.integer(local_column)) &&
      all(local_group >= 1L) &&
      all(local_column >= 1L) &&
      all(local_group <= layout$n_groups) &&
      all(local_column <= layout$global_n_columns)
    if(!isTRUE(valid_local_indices)){
      stop(
        "Random-effect local latent metadata",
        .bt_random_effect_metadata_block_detail(random_term),
        " must contain one positive group and column index per local latent ",
        "cell within the stored dimensions.",
        call. = FALSE
      )
    }
    keep <- local_group <= n_groups & local_column <= n_columns
    return(node_names[keep])
  }

  outer(
    seq_len(n_groups),
    seq_len(n_columns),
    Vectorize(function(group, column){
      paste0(random_term$parameter_stem, "_xRE_Zx[", group, ",", column, "]")
    })
  )
}

.bt_random_effect_sd_draws <- function(random_term, n_columns, posterior,
                                       prior_list){

  posterior <- .bt_random_effect_marginal_covariance_validate_posterior(posterior)
  if(!is.null(random_term$sd_binding) &&
     isTRUE(random_term$sd_binding$true_allocation)){
    .bt_check_random_sd_binding(random_term$sd_binding)
    allocation <- random_term$sd_binding$allocations[[1L]]
    target <- .bt_random_effect_allocation_target_metadata(allocation)
    .bt_random_effect_allocation_scale_metadata(allocation)
    if(identical(target, "sd_component")){
      .bt_check_random_sd_component_binding(
        binding = random_term$sd_binding,
        n_columns = n_columns,
        context = "Random-effect allocation metadata"
      )
    }else{
      .bt_random_effect_allocation_factors_metadata(allocation)
    }
    allocated <- .bt_random_effect_allocated_sd_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior,
      prior_list = prior_list
    )
    if(!is.null(allocated)){
      return(allocated)
    }
  }

  sd_names <- random_term$sd_parameter_names
  if(is.null(sd_names) || length(sd_names) != n_columns || any(is.na(sd_names))){
    return(NULL)
  }

  out <- matrix(NA_real_, nrow = nrow(posterior), ncol = n_columns)
  for(column in seq_len(n_columns)){
    values <- .bt_random_effect_parameter_draws(
      parameter_name = sd_names[column],
      posterior = posterior,
      prior_list = prior_list
    )
    if(is.null(values)){
      return(NULL)
    }
    out[, column] <- values
  }

  out
}

.bt_random_effect_allocation_target_metadata <- function(
    allocation,
    context = "Random-effect allocation metadata"){

  target <- allocation$target
  if(is.character(target) && length(target) == 1L &&
     !is.na(target) && target %in% c("block", "sd_component")){
    return(target)
  }

  stop(
    context,
    " are missing canonical 'allocation$target'.",
    call. = FALSE
  )
}

.bt_random_effect_allocation_scale_metadata <- function(
    allocation,
    context = "Random-effect allocation metadata"){

  scale <- allocation$scale
  if(is.character(scale) && length(scale) == 1L && !is.na(scale) &&
     scale %in% c("total_variance", "mean_variance")){
    return(scale)
  }

  stop(
    context,
    " are missing canonical 'allocation$scale'.",
    call. = FALSE
  )
}

.bt_random_effect_allocation_factors_metadata <- function(
    allocation,
    context = "Random-effect allocation metadata"){

  factors <- allocation$factors
  if(is.list(factors)){
    return(factors)
  }

  stop(
    context,
    " are missing canonical 'allocation$factors'.",
    call. = FALSE
  )
}

.bt_random_effect_allocation_parent_factors_metadata <- function(
    allocation,
    context = "Random-effect allocation metadata"){

  parent_factors <- allocation$parent_factors
  if(is.list(parent_factors)){
    return(parent_factors)
  }

  stop(
    context,
    " are missing canonical 'allocation$parent_factors'.",
    call. = FALSE
  )
}

.bt_random_effect_parameter_draws <- function(parameter_name, posterior,
                                             prior_list){

  if(parameter_name %in% colnames(posterior)){
    return(posterior[, parameter_name])
  }

  prior_name <- sub("\\[[0-9]+\\]$", "", parameter_name)
  if(!prior_name %in% names(prior_list)){
    return(NULL)
  }
  prior <- prior_list[[prior_name]]
  if(!is.prior.point(prior)){
    return(NULL)
  }

  location <- prior$parameters[["location"]]
  if(length(location) != 1L || is.na(location)){
    return(NULL)
  }

  rep(location, nrow(posterior))
}

.bt_random_effect_allocated_sd_draws <- function(random_term, n_columns,
                                                 posterior, prior_list){

  binding <- random_term$sd_binding
  if(is.null(binding) || !isTRUE(binding$true_allocation) ||
     length(binding$allocations) == 0L){
    return(NULL)
  }
  allocation <- binding$allocations[[1L]]
  target <- .bt_random_effect_allocation_target_metadata(allocation)
  scale <- .bt_random_effect_allocation_scale_metadata(allocation)
  factors <- if(identical(target, "sd_component")){
    .bt_random_effect_allocation_parent_factors_metadata(allocation)
  }else{
    .bt_random_effect_allocation_factors_metadata(allocation)
  }

  base <- .bt_random_effect_parameter_draws(
    parameter_name = .bt_random_sd_binding_source_name(allocation$source),
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(base)){
    return(NULL)
  }

  base <- .bt_random_effect_apply_allocation_factors(
    base = base,
    factors = factors,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(base)){
    return(NULL)
  }

  if(!identical(target, "sd_component")){
    return(matrix(base, nrow = nrow(posterior), ncol = n_columns))
  }

  sd_component_metadata <- .bt_check_random_sd_component_allocation(
    allocation = allocation,
    n_columns = n_columns,
    context = "Random-effect allocation metadata"
  )
  weights <- .bt_random_effect_dirichlet_draws(
    parameter_name = allocation$weight_name,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(weights)){
    return(NULL)
  }
  leaf_index <- sd_component_metadata$leaf_index_by_column
  K <- sd_component_metadata$n_targets
  if(ncol(weights) != K){
    stop(
      "Random-effect allocation metadata for '",
      allocation$weight_name,
      "' expected ", K,
      " Dirichlet coordinate(s), but found ", ncol(weights), ".",
      call. = FALSE
    )
  }
  out <- matrix(NA_real_, nrow = nrow(posterior), ncol = n_columns)
  for(column in seq_len(n_columns)){
    out[, column] <- base * .bt_random_effect_allocation_multiplier(
      weights = weights[, leaf_index[column]],
      scale = scale,
      n_targets = K
    )
  }

  out
}

.bt_random_effect_apply_allocation_factors <- function(base, factors,
                                                       posterior,
                                                       prior_list){

  factor_plan <- .bt_random_effect_allocation_factor_plan(factors)
  .bt_random_effect_apply_allocation_factor_plan(
    base = base,
    factor_plan = factor_plan,
    posterior = posterior,
    prior_list = prior_list
  )
}

.bt_random_effect_allocation_factor_plan <- function(factors){

  factor_plan <- attr(
    factors,
    "BayesTools_random_effect_allocation_factor_plan",
    exact = TRUE
  )
  if(!is.null(factor_plan)){
    return(factor_plan)
  }

  .bt_random_effect_compile_allocation_factor_plan(factors)
}

.bt_random_effect_allocation_factors_with_plan <- function(factors){

  attr(
    factors,
    "BayesTools_random_effect_allocation_factor_plan"
  ) <- .bt_random_effect_compile_allocation_factor_plan(factors)

  factors
}

.bt_random_effect_compile_allocation_factor_plan <- function(factors){

  if(!is.list(factors)){
    stop(
      "Random-effect allocation metadata are missing canonical 'allocation$factors'.",
      call. = FALSE
    )
  }

  factor_plan <- vector("list", length(factors))
  for(factor_i in seq_along(factors)){
    factor <- factors[[factor_i]]
    if(!is.list(factor)){
      stop(
        "Random-effect allocation factor metadata are missing canonical fields.",
        call. = FALSE
      )
    }
    if(!is.character(factor$weight_name) || length(factor$weight_name) != 1L ||
       is.na(factor$weight_name) || !nzchar(factor$weight_name)){
      stop(
        "Random-effect allocation factor metadata are missing canonical 'weight_name'.",
        call. = FALSE
      )
    }
    if(!is.numeric(factor$index) || length(factor$index) != 1L ||
       is.na(factor$index) || factor$index != as.integer(factor$index) ||
       factor$index < 1L){
      stop(
        "Random-effect allocation factor metadata are missing canonical 'index'.",
        call. = FALSE
      )
    }
    if(!is.numeric(factor$n_targets) || length(factor$n_targets) != 1L ||
       is.na(factor$n_targets) || factor$n_targets != as.integer(factor$n_targets) ||
       factor$n_targets < 2L){
      stop(
        "Random-effect allocation factor metadata are missing canonical 'n_targets'.",
        call. = FALSE
      )
    }
    factor_plan[[factor_i]] <- list(
      weight_name = factor$weight_name,
      index = factor$index,
      scale = .bt_random_effect_allocation_scale_metadata(
        factor,
        context = "Random-effect allocation factor metadata"
      ),
      n_targets = factor$n_targets,
      inclusion_name = factor$inclusion_name
    )
    if(!is.null(factor_plan[[factor_i]]$inclusion_name) &&
       (!is.character(factor_plan[[factor_i]]$inclusion_name) ||
        length(factor_plan[[factor_i]]$inclusion_name) != 1L ||
        is.na(factor_plan[[factor_i]]$inclusion_name) ||
        !nzchar(factor_plan[[factor_i]]$inclusion_name))){
      stop(
        "Random-effect allocation factor metadata are missing canonical 'inclusion_name'.",
        call. = FALSE
      )
    }
  }

  factor_plan
}

.bt_random_effect_apply_allocation_factor_plan <- function(base, factor_plan,
                                                           posterior,
                                                           prior_list){

  if(length(factor_plan) == 0L){
    return(base)
  }

  out <- base
  for(factor in factor_plan){
    scale <- .bt_random_effect_allocation_scale_metadata(
      factor,
      context = "Random-effect allocation factor metadata"
    )
    weights <- .bt_random_effect_dirichlet_draws(
      parameter_name = factor$weight_name,
      posterior = posterior,
      prior_list = prior_list
    )
    if(is.null(weights)){
      return(NULL)
    }
    if(ncol(weights) != factor$n_targets){
      stop(
        "Random-effect allocation factor metadata for '",
        factor$weight_name,
        "' expected ", factor$n_targets,
        " Dirichlet coordinate(s), but found ", ncol(weights), ".",
        call. = FALSE
      )
    }
    if(factor$index > ncol(weights)){
      stop(
        "Random-effect allocation factor metadata for '",
        factor$weight_name,
        "' reference coordinate ", factor$index,
        ", but only ", ncol(weights), " coordinate(s) are available.",
        call. = FALSE
      )
    }
    gate <- .bt_random_effect_allocation_gate_draws(
      parameter_name = factor$inclusion_name,
      posterior = posterior
    )
    if(is.null(gate)){
      stop(
        "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
        factor$inclusion_name,
        "'.",
        call. = FALSE
      )
    }
    out <- out * .bt_random_effect_allocation_multiplier(
      weights = weights[, factor$index],
      scale = scale,
      n_targets = factor$n_targets
    ) * gate
  }

  out
}

.bt_random_effect_allocation_gate_draws <- function(parameter_name, posterior){

  if(is.null(parameter_name)){
    return(rep(1, nrow(posterior)))
  }
  if(!parameter_name %in% colnames(posterior)){
    return(NULL)
  }

  values <- as.numeric(posterior[, parameter_name])
  invalid <- !is.finite(values) | !(values %in% c(0, 1))
  if(any(invalid)){
    .bt_random_effect_allocation_out_of_support(
      "Random-effect allocation inclusion samples for '",
      parameter_name,
      "' must be Bernoulli indicators."
    )
  }

  values
}

.bt_random_effect_allocation_factor_chain_label <- function(factors){

  if(!is.list(factors) || length(factors) == 0L){
    return("")
  }

  labels <- vapply(factors, function(factor){
    if(!is.list(factor) ||
       !is.character(factor$weight_name) ||
       length(factor$weight_name) != 1L ||
       is.na(factor$weight_name) ||
       !is.numeric(factor$index) ||
       length(factor$index) != 1L ||
       is.na(factor$index)){
      return("<malformed>")
    }
    label <- paste0(factor$weight_name, "[", factor$index, "]")
    if(is.character(factor$inclusion_name) &&
       length(factor$inclusion_name) == 1L &&
       !is.na(factor$inclusion_name) &&
       nzchar(factor$inclusion_name)){
      label <- paste0(label, " * ", factor$inclusion_name)
    }
    label
  }, character(1))

  paste0(" (", paste(labels, collapse = " -> "), ")")
}

.bt_random_effect_allocation_multiplier <- function(weights, scale, n_targets){

  if(identical(scale, "mean_variance")){
    if(!is.numeric(n_targets) || length(n_targets) != 1L ||
       is.na(n_targets) || n_targets < 1L){
      stop(
        "Random-effect allocation metadata are missing canonical 'allocation$n_targets'.",
        call. = FALSE
      )
    }
    return(sqrt(n_targets * weights))
  }
  if(identical(scale, "total_variance")){
    return(sqrt(weights))
  }

  stop(
    "Random-effect allocation metadata are missing canonical 'allocation$scale'.",
    call. = FALSE
  )
}

.bt_random_effect_dirichlet_draws <- function(parameter_name, posterior,
                                             prior_list){

  if(!parameter_name %in% names(prior_list)){
    return(NULL)
  }
  prior <- prior_list[[parameter_name]]
  if(!is.prior.simplex(prior) || !identical(prior$distribution, "dirichlet")){
    return(NULL)
  }

  K <- prior$parameters[["K"]]
  cache <- .bt_random_effect_dirichlet_draw_cache(posterior)
  weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
  if(all(weight_names %in% colnames(posterior))){
    cache_key <- .bt_random_effect_dirichlet_cache_key(parameter_name, K, "weights")
    if(!is.null(cache) && exists(cache_key, envir = cache, inherits = FALSE)){
      return(get(cache_key, envir = cache, inherits = FALSE))
    }
    weights <- .bt_random_effect_validate_dirichlet_weights(
      weights = posterior[, weight_names, drop = FALSE],
      parameter_name = parameter_name
    )
    .bt_random_effect_dirichlet_cache_assign(cache, cache_key, weights)
    return(weights)
  }

  eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(K), "]")
  if(!all(eta_names %in% colnames(posterior))){
    return(NULL)
  }

  cache_key <- .bt_random_effect_dirichlet_cache_key(parameter_name, K, "eta")
  if(!is.null(cache) && exists(cache_key, envir = cache, inherits = FALSE)){
    return(get(cache_key, envir = cache, inherits = FALSE))
  }

  eta <- posterior[, eta_names, drop = FALSE]
  invalid <- !is.finite(eta) | eta <= 0
  if(any(invalid)){
    invalid_column <- col(eta)[which(invalid)[1L]]
    .bt_random_effect_allocation_out_of_support(
      "Random-effect Dirichlet allocation auxiliary samples must be positive for '",
      eta_names[invalid_column],
      "'."
    )
  }
  weights <- .bt_random_effect_validate_dirichlet_weights(
    weights = eta / rowSums(eta),
    parameter_name = parameter_name
  )
  .bt_random_effect_dirichlet_cache_return(
    weights = weights,
    cache = cache,
    cache_key = cache_key
  )
}

.bt_random_effect_dirichlet_draw_cache <- function(posterior){

  cache <- attr(
    posterior,
    "BayesTools_random_effect_dirichlet_draw_cache",
    exact = TRUE
  )
  if(is.environment(cache)){
    return(cache)
  }

  NULL
}

.bt_random_effect_dirichlet_cache_key <- function(parameter_name, K, mode){

  paste0(parameter_name, "\r", K, "\r", mode)
}

.bt_random_effect_dirichlet_cache_assign <- function(cache, cache_key, weights){

  if(!is.null(cache)){
    assign(cache_key, weights, envir = cache)
  }

  invisible(weights)
}

.bt_random_effect_dirichlet_cache_return <- function(weights, cache, cache_key){

  .bt_random_effect_dirichlet_cache_assign(
    cache = cache,
    cache_key = cache_key,
    weights = weights
  )

  weights
}

.bt_random_effect_validate_dirichlet_weights <- function(weights,
                                                         parameter_name){

  if(!is.matrix(weights) || ncol(weights) < 2L){
    .bt_random_effect_allocation_out_of_support(
      "Random-effect Dirichlet allocation samples for '",
      parameter_name,
      "' must be a matrix with at least two coordinates."
    )
  }
  invalid <- !is.finite(weights) | weights < 0
  if(any(invalid)){
    invalid_column <- col(weights)[which(invalid)[1L]]
    invalid_name <- colnames(weights)[invalid_column]
    if(is.null(invalid_name) || is.na(invalid_name) || !nzchar(invalid_name)){
      invalid_name <- paste0(parameter_name, "[", invalid_column, "]")
    }
    .bt_random_effect_allocation_out_of_support(
      "Random-effect Dirichlet allocation samples must be finite and non-negative for '",
      invalid_name,
      "'."
    )
  }
  canonical <- tryCatch(
    .canonicalize_simplex(
      weights,
      name = paste0("Dirichlet allocation samples for '", parameter_name, "'"),
      diagnostics = TRUE
    ),
    error = function(e) e
  )
  if(inherits(canonical, "error")){
    .bt_random_effect_allocation_out_of_support(
      "Random-effect Dirichlet allocation samples for '",
      parameter_name, "' are not on the simplex: ",
      conditionMessage(canonical)
    )
  }

  weights <- canonical$values
  if(canonical$diagnostics$max_correction > 0){
    attr(weights, "simplex_canonicalization") <- canonical$diagnostics
  }
  weights
}


.bt_random_effect_allocation_out_of_support <- function(...){

  stop(structure(
    list(
      message = paste0(...),
      call = NULL
    ),
    class = c(
      "BayesTools_random_effect_allocation_out_of_support",
      "error",
      "condition"
    )
  ))
}
