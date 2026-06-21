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

