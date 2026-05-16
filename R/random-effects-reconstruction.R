# Internal random-effect posterior reconstruction helpers.

.bt_try_random_effect_contribution_from_latent <- function(random_term,
                                                          model_matrix,
                                                          group_map,
                                                          posterior,
                                                          prior_list){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  n_groups <- length(random_term$group_levels)

  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(!all(as.vector(z_names) %in% colnames(posterior))){
    return(NULL)
  }

  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(sd_draws)){
    return(NULL)
  }

  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )
  if(is.null(cholesky)){
    return(NULL)
  }

  z_draws <- lapply(seq_len(n_columns), function(latent_column){
    posterior[, z_names[, latent_column], drop = FALSE]
  })

  .bt_random_effect_contribution_from_latent_draws(
    model_matrix = model_matrix,
    group_map = group_map,
    z_draws = z_draws,
    sd_draws = sd_draws,
    cholesky = cholesky
  )
}

.bt_random_effect_contribution_from_latent_draws <- function(model_matrix,
                                                             group_map,
                                                             z_draws,
                                                             sd_draws,
                                                             cholesky){

  n_draws <- nrow(sd_draws)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  n_groups <- ncol(z_draws[[1L]])
  output <- matrix(0, nrow = n_rows, ncol = n_draws)

  for(column in seq_len(n_columns)){
    coefficient_matrix <- matrix(0, nrow = n_draws, ncol = n_groups)
    for(latent_column in seq_len(n_columns)){
      L_ij <- cholesky[, column, latent_column]
      if(all(L_ij == 0)){
        next
      }
      z_matrix <- z_draws[[latent_column]]
      coefficient_matrix <- coefficient_matrix +
        z_matrix * matrix(L_ij, nrow = n_draws, ncol = n_groups)
    }
    coefficient_matrix <- coefficient_matrix *
      matrix(sd_draws[, column], nrow = n_draws, ncol = n_groups)
    output <- output +
      t(coefficient_matrix[, group_map, drop = FALSE]) *
      matrix(model_matrix[, column], nrow = n_rows, ncol = n_draws)
  }

  output
}

.bt_random_effect_latent_names <- function(random_term, n_groups, n_columns){

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

  sd_names <- random_term$sd_parameter_names
  if(is.null(sd_names) || length(sd_names) != n_columns || any(is.na(sd_names))){
    return(NULL)
  }

  if(!is.null(random_term$allocation)){
    .bt_random_effect_allocation_components_metadata(random_term$allocation)
    .bt_random_effect_allocation_scale_metadata(random_term$allocation)
    .bt_random_effect_allocation_factors_metadata(random_term$allocation)
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

.bt_random_effect_allocation_components_metadata <- function(
    allocation,
    context = "Random-effect allocation metadata"){

  components <- allocation$components
  if(is.character(components) && length(components) == 1L &&
     !is.na(components) && components %in% c("block", "sd")){
    return(components)
  }

  stop(
    context,
    " are missing canonical 'allocation$components'.",
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

  allocation <- random_term$allocation
  if(is.null(allocation)){
    return(NULL)
  }
  components <- .bt_random_effect_allocation_components_metadata(allocation)
  scale <- .bt_random_effect_allocation_scale_metadata(allocation)
  factors <- .bt_random_effect_allocation_factors_metadata(allocation)

  base <- .bt_random_effect_parameter_draws(
    parameter_name = allocation$source_name,
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

  if(!identical(components, "sd")){
    return(matrix(base, nrow = nrow(posterior), ncol = n_columns))
  }

  weights <- .bt_random_effect_dirichlet_draws(
    parameter_name = allocation$weight_name,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(weights)){
    return(NULL)
  }
  leaf_index <- allocation$leaf_index_by_column
  if(!is.numeric(leaf_index) || length(leaf_index) != n_columns ||
     any(is.na(leaf_index))){
    stop(
      "Random-effect allocation metadata are missing canonical 'allocation$leaf_index_by_column'.",
      call. = FALSE
    )
  }
  K <- allocation$n_targets
  if(!is.numeric(K) || length(K) != 1L || is.na(K) || K < 1L){
    stop(
      "Random-effect allocation metadata are missing canonical 'allocation$n_targets'.",
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

  if(!is.list(factors)){
    stop(
      "Random-effect allocation metadata are missing canonical 'allocation$factors'.",
      call. = FALSE
    )
  }
  if(length(factors) == 0L){
    return(base)
  }
  out <- base
  for(factor in factors){
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
       is.na(factor$index)){
      stop(
        "Random-effect allocation factor metadata are missing canonical 'index'.",
        call. = FALSE
      )
    }
    scale <- .bt_random_effect_allocation_scale_metadata(
      factor,
      context = "Random-effect allocation factor metadata"
    )
    if(!is.numeric(factor$n_targets) || length(factor$n_targets) != 1L ||
       is.na(factor$n_targets) || factor$n_targets < 1L){
      stop(
        "Random-effect allocation factor metadata are missing canonical 'n_targets'.",
        call. = FALSE
      )
    }
    weights <- .bt_random_effect_dirichlet_draws(
      parameter_name = factor$weight_name,
      posterior = posterior,
      prior_list = prior_list
    )
    if(is.null(weights)){
      return(NULL)
    }
    out <- out * .bt_random_effect_allocation_multiplier(
      weights = weights[, factor$index],
      scale = scale,
      n_targets = factor$n_targets
    )
  }

  out
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
  weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
  if(all(weight_names %in% colnames(posterior))){
    return(posterior[, weight_names, drop = FALSE])
  }

  eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(K), "]")
  if(!all(eta_names %in% colnames(posterior))){
    return(NULL)
  }

  eta <- posterior[, eta_names, drop = FALSE]
  invalid <- !is.finite(eta) | eta <= 0
  if(any(invalid)){
    invalid_column <- col(eta)[which(invalid)[1L]]
    stop(
      "Random-effect Dirichlet allocation auxiliary samples must be positive for '",
      eta_names[invalid_column],
      "'.",
      call. = FALSE
    )
  }
  eta / rowSums(eta)
}

.bt_random_effect_cholesky_draws <- function(random_term, n_columns,
                                            posterior){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect posterior reconstruction metadata"
  )
  if(structure %in% c("diag", "id") || n_columns == 1L){
    out <- array(0, dim = c(nrow(posterior), n_columns, n_columns))
    for(column in seq_len(n_columns)){
      out[, column, column] <- 1
    }
    return(out)
  }
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Random-effect posterior reconstruction metadata"
  )

  L_names <- .bt_random_effect_cholesky_names(
    random_term = random_term,
    n_columns = n_columns
  )
  if(all(as.vector(L_names) %in% colnames(posterior))){
    out <- array(NA_real_, dim = c(nrow(posterior), n_columns, n_columns))
    for(row in seq_len(n_columns)){
      for(column in seq_len(n_columns)){
        out[, row, column] <- posterior[, L_names[row, column]]
      }
    }
    return(out)
  }

  if(identical(structure, "us")){
    u_names <- .bt_random_effect_lkj_primitive_names(random_term, n_columns)
    if(all(u_names %in% colnames(posterior))){
      return(.bt_lkj_cholesky_cpc_u_to_L(
        posterior[, u_names, drop = FALSE],
        K = n_columns
      ))
    }
  }

  if(structure %in% c("cs", "hcs", "ar1", "car", "har")){
    rho <- .bt_random_effect_rho_draws(random_term, posterior)
    if(!is.null(rho)){
      out <- array(NA_real_, dim = c(nrow(posterior), n_columns, n_columns))
      for(draw in seq_len(nrow(posterior))){
        R <- .bt_random_effect_structured_correlation_matrix(
          structure = structure,
          K = n_columns,
          rho = rho[draw],
          distance_matrix = if(identical(structure, "car")) correlation$distance_matrix else NULL
        )
        out[draw, , ] <- t(chol(R))
      }
      return(out)
    }
  }

  NULL
}

.bt_random_effect_rho_draws <- function(random_term, posterior){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect posterior reconstruction metadata"
  )
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Random-effect posterior reconstruction metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    return(NULL)
  }
  if(!is.character(correlation$rho_name) || length(correlation$rho_name) != 1L ||
     is.na(correlation$rho_name) || !nzchar(correlation$rho_name)){
    stop(
      "Random-effect posterior reconstruction metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical 'random_term$correlation$rho_name'.",
      call. = FALSE
    )
  }
  if(!is.character(correlation$sample_name) || length(correlation$sample_name) != 1L ||
     is.na(correlation$sample_name) || !nzchar(correlation$sample_name)){
    stop(
      "Random-effect posterior reconstruction metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical 'random_term$correlation$sample_name'.",
      call. = FALSE
    )
  }

  rho_scale <- .bt_random_effect_rho_scale_metadata(correlation, random_term)
  if(!identical(rho_scale, "rho") &&
     correlation$sample_name %in% colnames(posterior)){
    rho <- .bt_random_effect_transform_rho(
      posterior[, correlation$sample_name],
      correlation = correlation,
      random_term = random_term
    )
  }else if(correlation$rho_name %in% colnames(posterior)){
    rho <- posterior[, correlation$rho_name]
  }else if(correlation$sample_name %in% colnames(posterior)){
    rho <- .bt_random_effect_transform_rho(
      posterior[, correlation$sample_name],
      correlation = correlation,
      random_term = random_term
    )
  }else{
    sample_fixed <- .bt_random_effect_rho_fixed_sample_metadata(
      correlation,
      random_term
    )
    if(is.null(sample_fixed)){
      return(NULL)
    }
    rho <- .bt_random_effect_transform_rho(
      rep(sample_fixed, nrow(posterior)),
      correlation = correlation,
      random_term = random_term
    )
  }

  bounds <- .bt_random_effect_rho_bounds_metadata(correlation, random_term)
  if(any(rho <= bounds[["lower"]] | rho >= bounds[["upper"]])){
    return(NULL)
  }

  rho
}

.bt_random_effect_transform_rho <- function(value, correlation,
                                            random_term = NULL){

  rho_scale <- .bt_random_effect_rho_scale_metadata(correlation, random_term)
  if(identical(rho_scale, "fisher_z")){
    return(tanh(value))
  }
  if(identical(rho_scale, "logit")){
    bounds <- .bt_random_effect_rho_bounds_metadata(correlation, random_term)
    return(bounds[["lower"]] + (bounds[["upper"]] - bounds[["lower"]]) * stats::plogis(value))
  }

  value
}

.bt_random_effect_rho_scale_metadata <- function(correlation,
                                                 random_term = NULL){

  rho_scale <- correlation$rho_scale
  if(!is.character(rho_scale) || length(rho_scale) != 1L ||
     is.na(rho_scale) || !rho_scale %in% c("fisher_z", "logit", "rho")){
    stop(
      "Random-effect posterior reconstruction metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical 'random_term$correlation$rho_scale'.",
      call. = FALSE
    )
  }

  rho_scale
}

.bt_random_effect_rho_fixed_sample_metadata <- function(correlation,
                                                        random_term = NULL){

  sample_fixed <- correlation$sample_fixed
  if(is.null(sample_fixed)){
    return(NULL)
  }
  if(!is.numeric(sample_fixed) || length(sample_fixed) != 1L ||
     is.na(sample_fixed)){
    stop(
      "Random-effect posterior reconstruction metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation$sample_fixed'.",
      call. = FALSE
    )
  }

  sample_fixed
}

.bt_random_effect_rho_bounds_metadata <- function(correlation,
                                                  random_term = NULL){

  bounds <- correlation$bounds
  if((is.list(bounds) || is.numeric(bounds)) &&
     is.numeric(bounds[["lower"]]) && length(bounds[["lower"]]) == 1L &&
     !is.na(bounds[["lower"]]) &&
     is.numeric(bounds[["upper"]]) && length(bounds[["upper"]]) == 1L &&
     !is.na(bounds[["upper"]]) &&
     bounds[["lower"]] < bounds[["upper"]]){
    return(bounds)
  }

  stop(
    "Random-effect posterior reconstruction metadata",
    .bt_random_effect_metadata_block_detail(random_term),
    " are missing canonical 'random_term$correlation$bounds'.",
    call. = FALSE
  )
}

.bt_random_effect_structured_correlation_matrix <- function(structure, K, rho,
                                                           distance_matrix = NULL){

  R <- matrix(NA_real_, nrow = K, ncol = K)
  if(identical(structure, "car")){
    distance_matrix <- .bt_random_effect_validate_car_distance_matrix(distance_matrix, K)
  }

  for(row in seq_len(K)){
    for(column in seq_len(K)){
      R[row, column] <- if(row == column){
        1
      }else if(structure %in% c("cs", "hcs")){
        rho
      }else if(identical(structure, "car")){
        rho^distance_matrix[row, column]
      }else{
        rho^abs(row - column)
      }
    }
  }

  R
}

.bt_random_effect_cholesky_names <- function(random_term, n_columns){

  outer(
    seq_len(n_columns),
    seq_len(n_columns),
    Vectorize(function(row, column){
      paste0(random_term$parameter_stem, "_xRE_CORx_L[", row, ",", column, "]")
    })
  )
}

.bt_random_effect_lkj_primitive_names <- function(
    random_term,
    n_columns,
    context = "Random-effect posterior reconstruction metadata"){

  n_pairs <- n_columns * (n_columns - 1L) / 2L
  if(n_pairs < 1L){
    return(character(0))
  }

  correlation <- .bt_random_effect_correlation_metadata(
    random_term = random_term,
    structure = "us",
    context = context
  )
  if(is.null(correlation) || !identical(correlation$type, "lkj")){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical LKJ 'random_term$correlation'.",
      call. = FALSE
    )
  }

  primitive_names <- correlation$primitive_names
  if(!is.character(primitive_names) || length(primitive_names) != n_pairs ||
     any(is.na(primitive_names)) || any(!nzchar(primitive_names))){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical 'random_term$correlation$primitive_names'.",
      call. = FALSE
    )
  }

  primitive_bounds <- correlation$primitive_bounds
  if(!is.list(primitive_bounds) ||
     !all(c("lb", "ub") %in% names(primitive_bounds)) ||
     !is.numeric(primitive_bounds$lb) ||
     !is.numeric(primitive_bounds$ub) ||
     length(primitive_bounds$lb) != n_pairs ||
     length(primitive_bounds$ub) != n_pairs ||
     !identical(names(primitive_bounds$lb), primitive_names) ||
     !identical(names(primitive_bounds$ub), primitive_names) ||
     any(is.na(primitive_bounds$lb)) ||
     any(is.na(primitive_bounds$ub)) ||
     any(primitive_bounds$lb != 0) ||
     any(primitive_bounds$ub != 1)){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical 'random_term$correlation$primitive_bounds'.",
      call. = FALSE
    )
  }

  primitive_names
}

.bt_random_effect_coefficient_names <- function(random_term, n_groups, n_columns){

  outer(
    seq_len(n_groups),
    seq_len(n_columns),
    Vectorize(function(group, column){
      paste0(random_term$parameter_stem, "_xRE_COEFx[", group, ",", column, "]")
    })
  )
}
