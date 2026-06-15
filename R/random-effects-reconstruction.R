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

.bt_random_effect_row_indexed_contribution_from_latent <- function(
    random_term, model_matrix, group_map, posterior, prior_list,
    data = NULL, parameters = NULL, context = "Prediction"){

  source_draws <- .bt_random_effect_row_indexed_source_draws(
    random_term = random_term,
    n_rows = nrow(model_matrix),
    posterior = posterior,
    data = data,
    parameters = parameters,
    context = context
  )
  if(any(!is.finite(source_draws) | source_draws < 0)){
    stop(
      context, " with row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' returned non-finite or negative source values.",
      call. = FALSE
    )
  }
  column_allocation_draws <- .bt_random_effect_row_indexed_column_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    n_columns = ncol(model_matrix)
  )
  if(!is.null(column_allocation_draws)){
    if(any(!is.finite(column_allocation_draws) | column_allocation_draws < 0)){
      stop(
        context, " with row-indexed external SD source '",
        .bt_random_effect_external_sd_source_label(random_term),
        "' returned non-finite or negative allocation values.",
        call. = FALSE
      )
    }
    unit_columns <- .bt_try_random_effect_unit_column_contributions_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior
    )
    if(is.null(unit_columns)){
      stop(
        "Random-effect block '", random_term$block_name,
        "' with row-indexed external SD source '",
        .bt_random_effect_external_sd_source_label(random_term),
        "' cannot be reconstructed from the posterior samples. Refit with ",
        "random_monitor(latent = TRUE) and monitor the required correlation ",
        "coordinates before using JAGS_evaluate_formula() with random effects.",
        call. = FALSE
      )
    }
    return(.bt_random_effect_apply_row_indexed_source_to_unit_columns(
      unit_columns = unit_columns,
      source_draws = source_draws,
      allocation_draws = column_allocation_draws
    ))
  }
  allocation_draws <- .bt_random_effect_row_indexed_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list
  )
  if(any(!is.finite(allocation_draws) | allocation_draws < 0)){
    stop(
      context, " with row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' returned non-finite or negative allocation values.",
      call. = FALSE
    )
  }
  unit_contribution <- .bt_try_random_effect_unit_contribution_from_latent(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    posterior = posterior
  )
  if(is.null(unit_contribution)){
    stop(
      "Random-effect block '", random_term$block_name,
      "' with row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' cannot be reconstructed from the posterior samples. Refit with ",
      "random_monitor(latent = TRUE) and monitor the required correlation ",
      "coordinates before using JAGS_evaluate_formula() with random effects.",
      call. = FALSE
    )
  }

  unit_contribution *
    t(source_draws) *
    matrix(allocation_draws, nrow = nrow(model_matrix), ncol = nrow(posterior),
           byrow = TRUE)
}

.bt_try_random_effect_unit_contribution_from_latent <- function(random_term,
                                                               model_matrix,
                                                               group_map,
                                                               posterior){

  unit_columns <- .bt_try_random_effect_unit_column_contributions_from_latent(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    posterior = posterior
  )
  if(is.null(unit_columns)){
    return(NULL)
  }

  rowSums(unit_columns, dims = 2L)
}

.bt_try_random_effect_unit_column_contributions_from_latent <- function(random_term,
                                                                       model_matrix,
                                                                       group_map,
                                                                       posterior){

  n_draws <- nrow(posterior)
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
  sd_draws <- matrix(1, nrow = n_draws, ncol = n_columns)

  .bt_random_effect_column_contributions_from_latent_draws(
    model_matrix = model_matrix,
    group_map = group_map,
    z_draws = z_draws,
    sd_draws = sd_draws,
    cholesky = cholesky
  )
}

.bt_random_effect_row_indexed_source_draws <- function(random_term, n_rows,
                                                       posterior,
                                                       data = NULL,
                                                       parameters = NULL,
                                                       context = "Prediction"){

  source <- .bt_random_effect_row_indexed_source(random_term)
  source_names <- .bt_parameter_source_row_names(source$source, n_rows)
  source_values <- .bt_parameter_source_value_draws(
    source = source$source,
    n_rows = n_rows,
    posterior = posterior,
    data = data,
    parameters = parameters,
    context = context
  )
  if(!is.null(source_values)){
    return(source_values)
  }

  if(all(source_names %in% colnames(posterior))){
    return(posterior[, source_names, drop = FALSE])
  }

  missing <- !source_names %in% colnames(posterior)

  stop(
    context, " with row-indexed external SD source '",
    .bt_random_effect_external_sd_source_label(random_term),
    "' is missing values for prediction row(s): ",
    paste(which(missing), collapse = ", "),
    ". Expected posterior column(s): ",
    paste0("'", source_names[missing][seq_len(min(3L, sum(missing)))], "'",
           collapse = ", "),
    if(sum(missing) > 3L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_random_effect_row_indexed_source <- function(random_term){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    source <- NULL
  }else{
    .bt_check_random_sd_binding(binding)
    source <- binding$source
  }
  if(is.null(source) || !inherits(source, "random_sd_source")){
    stop(
      "Random-effect row-indexed external SD metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical source information.",
      call. = FALSE
    )
  }
  .bt_check_random_sd_source(source)
  if(!.bt_random_sd_source_is_row_indexed(source)){
    stop(
      "Random-effect row-indexed external SD metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical source information.",
      call. = FALSE
    )
  }

  source
}

.bt_random_effect_row_indexed_allocation_draws <- function(random_term,
                                                          posterior,
                                                          prior_list){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    stop(
      "Random-effect row-indexed external SD metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical SD binding information.",
      call. = FALSE
    )
  }
  .bt_check_random_sd_binding(binding)
  if(identical(binding$application, "column")){
    stop(
      "Random-effect row-indexed external SD metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " use column-specific SD allocation metadata and cannot be reconstructed with block-level allocation.",
      call. = FALSE
    )
  }
  factors <- .bt_random_effect_allocation_factors_metadata(binding)
  out <- .bt_random_effect_apply_allocation_factors(
    base = rep(1, nrow(posterior)),
    factors = factors,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(out)){
    allocation_label <- if(length(binding$allocations) > 0L){
      binding$allocations[[1L]]$label
    }else{
      "<direct>"
    }
    stop(
      "Prediction with row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' requires Dirichlet allocation coordinates for allocation '",
      allocation_label,
      "'.",
      call. = FALSE
    )
  }

  out
}

.bt_random_effect_row_indexed_column_allocation_draws <- function(random_term,
                                                                 posterior,
                                                                 prior_list,
                                                                 n_columns){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    return(NULL)
  }
  .bt_check_random_sd_binding(binding)
  if(isTRUE(binding$true_allocation)){
    allocation_target <- .bt_random_effect_allocation_target_metadata(
      binding$allocations[[1L]],
      context = paste0(
        "Random-effect row-indexed external SD metadata",
        .bt_random_effect_metadata_block_detail(random_term)
      )
    )
    if(identical(allocation_target, "sd_component")){
      .bt_check_random_sd_component_binding(
        binding = binding,
        n_columns = n_columns,
        context = paste0(
          "Random-effect row-indexed external SD metadata",
          .bt_random_effect_metadata_block_detail(random_term)
        )
      )
    }else if(length(binding$factors_by_column) > 0L){
      stop(
        "Random-effect row-indexed external SD metadata",
        .bt_random_effect_metadata_block_detail(random_term),
        " with column factor chains require 'allocation$target' to be 'sd_component'.",
        call. = FALSE
      )
    }
  }
  if(length(binding$factors_by_column) == 0L){
    return(NULL)
  }
  if(!identical(binding$application, "column")){
    stop(
      "Random-effect row-indexed external SD metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " contain column factor chains for a non-column SD binding.",
      call. = FALSE
    )
  }
  if(length(binding$factors_by_column) != n_columns){
    stop(
      "Random-effect row-indexed external SD metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " do not match the number of random-effect columns.",
      call. = FALSE
    )
  }

  out <- matrix(NA_real_, nrow = nrow(posterior), ncol = n_columns)
  for(column in seq_len(n_columns)){
    values <- .bt_random_effect_apply_allocation_factors(
      base = rep(1, nrow(posterior)),
      factors = binding$factors_by_column[[column]],
      posterior = posterior,
      prior_list = prior_list
    )
    if(is.null(values)){
      stop(
        "Prediction with row-indexed external SD source '",
        .bt_random_effect_external_sd_source_label(random_term),
        "' requires Dirichlet allocation coordinates for column ",
        column,
        " factor chain",
        .bt_random_effect_allocation_factor_chain_label(
          binding$factors_by_column[[column]]
        ),
        "'.",
        call. = FALSE
      )
    }
    out[, column] <- values
  }

  out
}

.bt_random_effect_contribution_from_latent_draws <- function(model_matrix,
                                                             group_map,
                                                             z_draws,
                                                             sd_draws,
                                                             cholesky){

  column_contributions <- .bt_random_effect_column_contributions_from_latent_draws(
    model_matrix = model_matrix,
    group_map = group_map,
    z_draws = z_draws,
    sd_draws = sd_draws,
    cholesky = cholesky
  )

  rowSums(column_contributions, dims = 2L)
}

.bt_random_effect_column_contributions_from_latent_draws <- function(model_matrix,
                                                                     group_map,
                                                                     z_draws,
                                                                     sd_draws,
                                                                     cholesky){

  n_draws <- nrow(sd_draws)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  n_groups <- ncol(z_draws[[1L]])
  output <- array(0, dim = c(n_rows, n_draws, n_columns))

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
    output[, , column] <-
      t(coefficient_matrix[, group_map, drop = FALSE]) *
      matrix(model_matrix[, column], nrow = n_rows, ncol = n_draws)
  }

  output
}

.bt_random_effect_apply_row_indexed_source_to_unit_columns <- function(
    unit_columns, source_draws, allocation_draws){

  if(length(dim(unit_columns)) != 3L){
    stop("'unit_columns' must be a three-dimensional array.", call. = FALSE)
  }
  n_rows <- dim(unit_columns)[1L]
  n_draws <- dim(unit_columns)[2L]
  n_columns <- dim(unit_columns)[3L]
  if(!is.matrix(source_draws) || nrow(source_draws) != n_draws ||
     ncol(source_draws) != n_rows){
    stop("Row-indexed source draws do not match random-effect contribution dimensions.", call. = FALSE)
  }
  if(is.null(dim(allocation_draws))){
    allocation_draws <- matrix(allocation_draws, ncol = 1L)
  }
  if(!is.matrix(allocation_draws) || nrow(allocation_draws) != n_draws ||
     !ncol(allocation_draws) %in% c(1L, n_columns)){
    stop("Row-indexed allocation draws do not match random-effect contribution dimensions.", call. = FALSE)
  }

  source_matrix <- t(source_draws)
  out <- matrix(0, nrow = n_rows, ncol = n_draws)
  for(column in seq_len(n_columns)){
    allocation_column <- if(ncol(allocation_draws) == 1L) 1L else column
    out <- out +
      unit_columns[, , column] *
      source_matrix *
      matrix(allocation_draws[, allocation_column], nrow = n_rows,
             ncol = n_draws, byrow = TRUE)
  }

  out
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
      n_targets = factor$n_targets
    )
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
    out <- out * .bt_random_effect_allocation_multiplier(
      weights = weights[, factor$index],
      scale = scale,
      n_targets = factor$n_targets
    )
  }

  out
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
    paste0(factor$weight_name, "[", factor$index, "]")
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
  row_sums <- rowSums(weights)
  if(any(!is.finite(row_sums) | abs(row_sums - 1) > 1e-8)){
    .bt_random_effect_allocation_out_of_support(
      "Random-effect Dirichlet allocation samples for '",
      parameter_name,
      "' must sum to one."
    )
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

.bt_random_effect_cholesky_draws <- function(random_term, n_columns,
                                            posterior,
                                            sample_space = FALSE){

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
    rho <- .bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior,
      sample_space = sample_space
    )
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

.bt_random_effect_rho_draws <- function(random_term, posterior,
                                        missing = c("null", "error"),
                                        out_of_support = c("null", "error"),
                                        sample_space = FALSE,
                                        context = "Random-effect posterior reconstruction metadata"){

  missing <- match.arg(missing)
  out_of_support <- match.arg(out_of_support)

  structure <- .bt_random_effect_structure(
    random_term,
    context = context
  )
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = context
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    return(NULL)
  }
  if(!is.character(correlation$rho_name) || length(correlation$rho_name) != 1L ||
     is.na(correlation$rho_name) || !nzchar(correlation$rho_name)){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical 'random_term$correlation$rho_name'.",
      call. = FALSE
    )
  }
  if(!is.character(correlation$sample_name) || length(correlation$sample_name) != 1L ||
     is.na(correlation$sample_name) || !nzchar(correlation$sample_name)){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical 'random_term$correlation$sample_name'.",
      call. = FALSE
    )
  }

  rho_scale <- .bt_random_effect_rho_scale_metadata(
    correlation,
    random_term,
    context = context
  )
  if(!identical(rho_scale, "rho") &&
     correlation$sample_name %in% colnames(posterior)){
    sample_value <- posterior[, correlation$sample_name]
    rho_source <- "sample"
    rho <- .bt_random_effect_transform_rho(
      sample_value,
      correlation = correlation,
      random_term = random_term,
      context = context
    )
  }else if(correlation$rho_name %in% colnames(posterior)){
    sample_value <- NULL
    rho_source <- "rho"
    rho <- posterior[, correlation$rho_name]
  }else if(correlation$sample_name %in% colnames(posterior)){
    sample_value <- posterior[, correlation$sample_name]
    rho_source <- "sample"
    rho <- .bt_random_effect_transform_rho(
      sample_value,
      correlation = correlation,
      random_term = random_term,
      context = context
    )
  }else{
    sample_fixed <- .bt_random_effect_rho_fixed_sample_metadata(
      correlation,
      random_term,
      context = context
    )
    if(is.null(sample_fixed)){
      if(identical(missing, "error")){
        .bt_random_effect_missing_rho_draws_stop(
          random_term = random_term,
          correlation = correlation,
          context = context
        )
      }
      return(NULL)
    }
    sample_value <- rep(sample_fixed, nrow(posterior))
    rho_source <- "fixed_sample"
    rho <- .bt_random_effect_transform_rho(
      sample_value,
      correlation = correlation,
      random_term = random_term,
      context = context
    )
  }

  bounds <- .bt_random_effect_rho_bounds_metadata(
    correlation,
    random_term,
    context = context
  )
  if(isTRUE(sample_space) &&
     identical(rho_source, "sample") &&
     !identical(rho_scale, "rho")){
    sample_bounds <- .bt_random_effect_rho_sample_bounds(
      correlation = correlation,
      random_term = random_term,
      context = context
    )
    invalid_sample <- .bt_random_effect_rho_outside_support(
      sample_value,
      sample_bounds,
      structure
    )
    if(any(invalid_sample)){
      if(identical(out_of_support, "error")){
        .bt_random_effect_rho_out_of_support_draw_stop(
          random_term = random_term,
          rho = rho,
          invalid = invalid_sample,
          bounds = bounds,
          structure = structure,
          context = context
        )
      }
      return(NULL)
    }
    rho <- .bt_random_effect_clamp_sample_space_rho(
      rho = rho,
      bounds = bounds,
      structure = structure
    )
  }
  invalid <- .bt_random_effect_rho_outside_support(rho, bounds, structure)
  if(any(invalid)){
    if(identical(out_of_support, "error")){
      .bt_random_effect_rho_out_of_support_draw_stop(
        random_term = random_term,
        rho = rho,
        invalid = invalid,
        bounds = bounds,
        structure = structure,
        context = context
      )
    }
    return(NULL)
  }

  rho
}

.bt_random_effect_rho_sample_bounds <- function(correlation,
                                                random_term = NULL,
                                                context = "Random-effect posterior reconstruction metadata"){

  rho_scale <- .bt_random_effect_rho_scale_metadata(
    correlation,
    random_term,
    context = context
  )
  bounds <- .bt_random_effect_rho_bounds_metadata(
    correlation,
    random_term,
    context = context
  )

  if(identical(rho_scale, "fisher_z")){
    lower <- if(bounds[["lower"]] <= -1) -Inf else atanh(bounds[["lower"]])
    upper <- if(bounds[["upper"]] >= 1)  Inf else atanh(bounds[["upper"]])
  }else if(identical(rho_scale, "logit")){
    lower <- -Inf
    upper <-  Inf
  }else{
    lower <- bounds[["lower"]]
    upper <- bounds[["upper"]]
  }

  c(lower = lower, upper = upper)
}

.bt_random_effect_clamp_sample_space_rho <- function(rho, bounds, structure){

  interval_width <- bounds[["upper"]] - bounds[["lower"]]
  if(!is.finite(interval_width) || interval_width <= 0){
    return(rho)
  }
  eps <- .Machine$double.eps * max(1, abs(bounds[["lower"]]), abs(bounds[["upper"]]))
  lower <- bounds[["lower"]] + eps
  upper <- bounds[["upper"]] - eps

  if(!.bt_random_effect_rho_lower_inclusive(structure)){
    rho <- ifelse(is.finite(rho) & rho <= bounds[["lower"]], lower, rho)
  }
  rho <- ifelse(is.finite(rho) & rho >= bounds[["upper"]], upper, rho)

  rho
}

.bt_random_effect_missing_rho_draws_stop <- function(random_term, correlation,
                                                     context){

  sample_names <- unique(c(correlation$rho_name, correlation$sample_name))
  sample_names <- sample_names[!is.na(sample_names) & nzchar(sample_names)]

  stop(
    context,
    " samples are missing canonical scalar correlation coordinates for block '",
    random_term$block_name,
    "'. Expected posterior column(s): ",
    paste0("'", sample_names, "'", collapse = ", "),
    ", or fixed rho metadata.",
    call. = FALSE
  )
}

.bt_random_effect_rho_out_of_support_draw_stop <- function(random_term, rho,
                                                           invalid, bounds,
                                                           structure,
                                                           context){

  draw <- which(invalid)[1L]
  interval <- paste0(
    if(.bt_random_effect_rho_lower_inclusive(structure)) "[" else "(",
    bounds[["lower"]],
    ", ",
    bounds[["upper"]],
    ")"
  )

  stop(
    context,
    " scalar correlation samples for block '",
    random_term$block_name,
    "' contain an out-of-support draw at row ",
    draw,
    " (rho = ",
    format(rho[draw], digits = 6),
    "). Expected rho in ",
    interval,
    ".",
    call. = FALSE
  )
}

.bt_random_effect_transform_rho <- function(value, correlation,
                                            random_term = NULL,
                                            context = "Random-effect posterior reconstruction metadata"){

  rho_scale <- .bt_random_effect_rho_scale_metadata(correlation, random_term, context)
  if(identical(rho_scale, "fisher_z")){
    return(tanh(value))
  }
  if(identical(rho_scale, "logit")){
    bounds <- .bt_random_effect_rho_bounds_metadata(correlation, random_term, context)
    return(bounds[["lower"]] + (bounds[["upper"]] - bounds[["lower"]]) * stats::plogis(value))
  }

  value
}

.bt_random_effect_rho_scale_metadata <- function(correlation,
                                                 random_term = NULL,
                                                 context = "Random-effect posterior reconstruction metadata"){

  rho_scale <- correlation$rho_scale
  if(!is.character(rho_scale) || length(rho_scale) != 1L ||
     is.na(rho_scale) || !rho_scale %in% c("fisher_z", "logit", "rho")){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical 'random_term$correlation$rho_scale'.",
      call. = FALSE
    )
  }

  rho_scale
}

.bt_random_effect_rho_fixed_sample_metadata <- function(correlation,
                                                        random_term = NULL,
                                                        context = "Random-effect posterior reconstruction metadata"){

  sample_fixed <- correlation$sample_fixed
  if(is.null(sample_fixed)){
    return(NULL)
  }
  if(!is.numeric(sample_fixed) || length(sample_fixed) != 1L ||
     is.na(sample_fixed)){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation$sample_fixed'.",
      call. = FALSE
    )
  }

  sample_fixed
}

.bt_random_effect_rho_bounds_metadata <- function(correlation,
                                                  random_term = NULL,
                                                  context = "Random-effect posterior reconstruction metadata"){

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
    context,
    .bt_random_effect_metadata_block_detail(random_term),
    " are missing canonical 'random_term$correlation$bounds'.",
    call. = FALSE
  )
}

.bt_random_effect_rho_outside_support <- function(rho, bounds, structure){

  lower_outside <- if(.bt_random_effect_rho_lower_inclusive(structure)){
    rho < bounds[["lower"]]
  }else{
    rho <= bounds[["lower"]]
  }

  !is.finite(rho) | lower_outside | rho >= bounds[["upper"]]
}

.bt_random_effect_rho_lower_inclusive <- function(structure){
  identical(structure, "car")
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
