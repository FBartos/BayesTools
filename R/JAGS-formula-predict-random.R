.bt_apply_formula_scale_to_data <- function(fit, parameter, data,
                                            predictors_type){

  continuous_predictors <- names(predictors_type[predictors_type == "continuous"])

  formula_scale <- attr(fit, "formula_scale", exact = TRUE)
  param_scale <- if(is.list(formula_scale)){
    formula_scale[[parameter]]
  }else{
    NULL
  }
  if(is.null(param_scale)){
    fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
    if(!is.null(fitted_design) && is.list(fitted_design$formula_scale)){
      param_scale <- fitted_design$formula_scale
    }
  }
  if(is.null(param_scale)){
    return(data)
  }

  scaled_predictors <- .formula_scale_strip_prefix(names(param_scale), parameter)
  continuous_predictors <- unique(c(continuous_predictors, scaled_predictors))
  if(length(continuous_predictors) == 0L){
    return(data)
  }

  for(continuous in continuous_predictors){
    if(!continuous %in% colnames(data)){
      next
    }
    scaled_name <- paste0(parameter, "_", continuous)
    if(scaled_name %in% names(param_scale)){
      scale_info <- param_scale[[scaled_name]]
      data[, continuous] <- (data[, continuous] - scale_info$mean) / scale_info$sd
    }
  }

  data
}

.bt_JAGS_evaluate_formula_with_random_effects <- function(fit, formula,
                                                          parameter, data,
                                                          prior_list,
                                                          posterior,
                                                          formula_target = NULL,
                                                          blocks = NULL,
                                                          new_levels = NULL,
                                                          fitted_rows = NULL,
                                                          data_supplied = FALSE,
                                                          replay_fitted_formula = FALSE,
                                                          expressions_to_eval = list()){

  fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
  if(is.null(fitted_design)){
    stop(
      "JAGS_evaluate_formula() needs fitted formula design metadata to evaluate random effects. ",
      "Use a fit produced by JAGS_fit() with formula_list.",
      call. = FALSE
    )
  }
  if(!.bt_formula_design_has_any_random_effects(fitted_design)){
    stop(
      "The supplied formula contains random effects, but the fitted formula for parameter '",
      parameter, "' does not.",
      call. = FALSE
    )
  }

  random_terms <- if(.has_random_effects(formula)){
    .bt_parse_random_effects(formula)$terms
  }else{
    list()
  }
  if(length(random_terms) > 0L){
    .bt_validate_random_effect_prediction_terms(
      requested = random_terms,
      fitted = .bt_formula_design_random_effects(fitted_design),
      parameter = parameter
    )
  }
  selected_random_effects <- .bt_JAGS_evaluate_formula_random_effect_terms(
    fitted_design = fitted_design,
    requested_terms = random_terms,
    formula_target = formula_target,
    blocks = blocks,
    replay_fitted_formula = replay_fitted_formula
  )

  fixed_formula <- .remove_expressions(.remove_random_effects(formula))
  output <- JAGS_evaluate_formula(
    fit = fit,
    formula = fixed_formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    formula_target = "fixed"
  )
  random_data <- .bt_apply_formula_scale_to_data(
    fit = fit,
    parameter = parameter,
    data = data,
    predictors_type = fitted_design$predictor_types
  )

  for(random_term in selected_random_effects){
    random_structure <- .bt_random_effect_structure(
      random_term,
      context = "Random-effect prediction metadata"
    )
    output <- output + .bt_JAGS_evaluate_random_effect_term(
      random_term = random_term,
      data = if(random_structure %in% c("cs", "hcs", "ar1", "car", "har")) data else random_data,
      group_data = data,
      posterior = posterior,
      prior_list = prior_list,
      new_levels = new_levels,
      fitted_rows = fitted_rows,
      data_supplied = data_supplied
    )
  }
  if(length(expressions_to_eval) > 0L){
    expression_data <- .bt_formula_expression_merge_data(
      data,
      if(!data_supplied) fitted_design$expression_data else NULL,
      context = paste0(
        "JAGS_evaluate_formula() for parameter '", parameter, "'"
      )
    )
    output <- output + .bt_formula_expression_contribution_matrix(
      expressions = expressions_to_eval,
      data = expression_data,
      n_rows = nrow(data),
      n_draws = nrow(posterior),
      context = paste0(
        "JAGS_evaluate_formula() for parameter '", parameter, "'"
      ),
      samples = posterior
    )
  }

  output
}

.bt_JAGS_evaluate_formula_random_effect_terms <- function(fitted_design,
                                                          requested_terms,
                                                          formula_target,
                                                          blocks,
                                                          replay_fitted_formula = FALSE){

  fitted_terms <- .bt_formula_design_random_effects(fitted_design)
  fitted_names <- vapply(fitted_terms, `[[`, character(1), "block_name")
  if(!is.null(blocks)){
    unknown <- setdiff(blocks, fitted_names)
    if(length(unknown) > 0L){
      stop(
        "Unknown random-effect block(s): ",
        paste(unknown, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if(is.null(formula_target) && isTRUE(replay_fitted_formula)){
    return(.bt_formula_design_sampled_random_effects(fitted_design))
  }

  requested_names <- vapply(requested_terms, `[[`, character(1), "block_name")
  selected_names <- if(!is.null(blocks)){
    blocks
  }else if(length(requested_names) > 0L){
    requested_names
  }else{
    fitted_names
  }
  selected_terms <- fitted_terms[match(selected_names, fitted_names)]
  modes <- vapply(selected_terms, .bt_random_effect_term_compile_mode, character(1))
  if(is.null(formula_target)){
    return(selected_terms[modes == "sampled"])
  }

  marginalized <- selected_names[modes == "marginalized"]
  if(length(marginalized) > 0L){
    stop(
      "JAGS_evaluate_formula() cannot condition on random-effect block(s) ",
      "compiled as marginalized: ",
      paste(marginalized, collapse = ", "),
      ". Use formula_target = \"marginal\" with JAGS_predict_formula() or refit ",
      "with the block(s) sampled.",
      call. = FALSE
    )
  }

  selected_terms
}

.bt_validate_random_effect_prediction_terms <- function(requested, fitted,
                                                        parameter){

  requested_names <- vapply(requested, function(term) term$block_name, character(1))
  if(anyDuplicated(requested_names)){
    stop(
      "Random-effect block names in the supplied formula must be unique.",
      call. = FALSE
    )
  }
  fitted_names <- vapply(fitted, function(term) term$block_name, character(1))
  unknown_names <- setdiff(requested_names, fitted_names)
  if(length(unknown_names) > 0L){
    stop(
      "Random-effect block(s) in the supplied formula do not match the fitted formula for parameter '",
      parameter, "'; block(s) were not found in the fitted formula: ",
      paste(unknown_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for(block in requested_names){
    requested_term <- requested[[match(block, requested_names)]]
    fitted_term <- fitted[[match(block, fitted_names)]]
    requested_structure <- .bt_random_effect_structure(
      requested_term,
      context = "Random-effect prediction metadata"
    )
    fitted_structure <- .bt_random_effect_structure(
      fitted_term,
      context = "Random-effect prediction metadata"
    )
    .bt_validate_random_effect_term_supported(requested_term)
    requested_homogeneous <- .bt_random_effect_homogeneous_sd(
      requested_term,
      requested_structure
    )
    fitted_homogeneous <- .bt_random_effect_homogeneous_sd_metadata(
      fitted_term,
      context = "Random-effect prediction metadata"
    )
    if(!identical(requested_structure, fitted_structure) ||
       !identical(requested_homogeneous, fitted_homogeneous) ||
       !identical(.bt_deparse_expr(requested_term$expr), .bt_deparse_expr(fitted_term$expr)) ||
       !identical(requested_term$group_label, fitted_term$group_label)){
      stop(
        "Random-effect block '", block,
        "' in the supplied formula does not match the fitted formula.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_JAGS_evaluate_random_effect_term <- function(random_term, data, posterior,
                                                prior_list,
                                                group_data = data,
                                                new_levels = NULL,
                                                fitted_rows = NULL,
                                                data_supplied = FALSE){

  new_levels <- .bt_random_effect_new_levels_policy(random_term, new_levels)
  prediction <- .bt_random_effect_prediction_data(
    random_term,
    data,
    group_data = group_data,
    allow_new_groups = isTRUE(new_levels$allow),
    context = "JAGS_evaluate_formula()"
  )
  model_matrix <- prediction$model_matrix
  group_map <- prediction$group_map
  n_rows <- nrow(model_matrix)
  n_draws <- nrow(posterior)
  .bt_random_effect_check_memory(
    estimate = .bt_random_effect_output_memory_estimate(
      operation = "conditional random-effect prediction",
      n_rows = n_rows,
      n_draws = n_draws
    ),
    block_name = random_term$block_name,
    alternative = paste0(
      "Reduce the number of prediction rows or posterior draws, or select ",
      "fewer random-effect blocks. Raise the option (or set it to Inf) only ",
      "after verifying the operation's memory budget."
    )
  )
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  fitted_n_groups <- length(random_term$group_levels)
  new_row <- group_map > fitted_n_groups
  prediction_rows <- if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    .bt_random_effect_prediction_fitted_rows(
      random_term = random_term,
      n_rows = n_rows,
      data_supplied = data_supplied,
      fitted_rows = fitted_rows,
      new_row = new_row,
      context = "JAGS_evaluate_formula()"
    )
  }else{
    NULL
  }

  if(any(!new_row)){
    existing_rows <- which(!new_row)
    output[existing_rows, ] <- .bt_JAGS_evaluate_random_effect_term_existing(
      random_term = random_term,
      model_matrix = model_matrix[existing_rows, , drop = FALSE],
      group_map = group_map[existing_rows],
      posterior = posterior,
      prior_list = prior_list,
      group_data = group_data[existing_rows, , drop = FALSE],
      prediction_rows = if(is.null(prediction_rows)){
        existing_rows
      }else{
        prediction_rows[existing_rows]
      }
    )
  }
  if(any(new_row) && identical(new_levels$method, "sample")){
    output <- output + .bt_random_effect_new_level_contribution_sample(
      random_term = random_term,
      prediction = prediction,
      new_row = new_row,
      posterior = posterior,
      prior_list = prior_list,
      source_data = group_data
    )
  }

  output
}

.bt_JAGS_evaluate_random_effect_term_existing <- function(random_term,
                                                          model_matrix,
                                                          group_map,
                                                          posterior,
                                                          prior_list,
                                                          group_data,
                                                          prediction_rows = NULL){

  n_columns <- ncol(model_matrix)
  n_groups <- length(random_term$group_levels)

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(
      .bt_random_effect_row_indexed_contribution_from_latent(
        random_term = random_term,
        model_matrix = model_matrix,
        group_map = group_map,
        posterior = posterior,
        prior_list = prior_list,
        data = group_data,
        prediction_rows = prediction_rows,
        context = "Prediction"
      )
    )
  }

  coefficient_names <- .bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(all(as.vector(coefficient_names) %in% colnames(posterior))){
    return(.bt_random_effect_contribution_from_coefficients(
      model_matrix = model_matrix,
      group_map = group_map,
      coefficient_names = coefficient_names,
      posterior = posterior
    ))
  }

  latent_contribution <- .bt_try_random_effect_contribution_from_latent(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    posterior = posterior,
    prior_list = prior_list
  )
  if(!is.null(latent_contribution)){
    return(latent_contribution)
  }

  stop(
    "Random-effect coefficients for block '", random_term$block_name,
    "' cannot be reconstructed from the posterior samples. Refit with ",
    "random_monitor(latent = TRUE) or random_monitor(coefficients = TRUE) ",
    "for that random-effect block before using JAGS_evaluate_formula() ",
    "with random effects.",
    call. = FALSE
  )
}

.bt_random_effect_contribution_from_coefficients <- function(model_matrix,
                                                            group_map,
                                                            coefficient_names,
                                                            posterior){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  for(column in seq_len(n_columns)){
    coefficient_matrix <- posterior[, coefficient_names[, column], drop = FALSE]
    output <- output +
      t(coefficient_matrix[, group_map, drop = FALSE]) *
      matrix(model_matrix[, column], nrow = n_rows, ncol = n_draws)
  }

  output
}

.bt_random_effect_new_level_contribution_sample <- function(random_term,
                                                            prediction,
                                                            new_row,
                                                            posterior,
                                                            prior_list,
                                                            source_data){

  .bt_random_effect_group_contribution_sample(
    random_term = random_term,
    model_matrix = prediction$model_matrix,
    group_map = prediction$group_map,
    rows = which(new_row),
    posterior = posterior,
    prior_list = prior_list,
    source_data = source_data
  )
}

.bt_random_effect_group_contribution_sample <- function(random_term,
                                                        model_matrix,
                                                        group_map,
                                                        rows = seq_len(nrow(model_matrix)),
                                                        posterior,
                                                        prior_list,
                                                        source_data,
                                                        prediction_rows = NULL){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  if(length(rows) == 0L){
    return(output)
  }

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(.bt_random_effect_group_contribution_sample_row_indexed(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      prior_list = prior_list,
      source_data = source_data,
      prediction_rows = prediction_rows
    ))
  }
  if(.bt_random_effect_has_known_group_covariance(random_term)){
    return(.bt_random_effect_group_contribution_sample_known_group_covariance(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      prior_list = prior_list
    ))
  }

  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(sd_draws)){
    .bt_random_effect_marginal_covariance_missing_sd_stop(
      random_term = random_term,
      n_columns = n_columns
    )
  }
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = sd_draws,
    n_draws = n_draws,
    n_columns = n_columns,
    label = "SD",
    random_term = random_term,
    nonnegative = TRUE,
    context = "Random-effect prediction"
  )
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect prediction metadata"
  )
  if(n_columns == 1L || structure %in% c("diag", "id")){
    return(.bt_random_effect_group_contribution_sample_independent(
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      column_scale_draws = sd_draws
    ))
  }
  rho_draws <- if(n_columns > 1L &&
                    structure %in% c("cs", "hcs", "ar1", "car", "har")){
    .bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior,
      missing = "null",
      context = "Random-effect prediction metadata"
    )
  }else{
    NULL
  }
  if(n_columns > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har") &&
     !is.null(rho_draws)){
    return(.bt_random_effect_group_contribution_sample_structured(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      sd_draws = sd_draws
    ))
  }
  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )
  if(is.null(cholesky)){
    .bt_random_effect_marginal_covariance_missing_correlation_stop(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
  }
  .bt_random_effect_marginal_covariance_validate_correlation_cholesky(
    cholesky = cholesky,
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    factor <- sweep(
      cholesky[draw, , ],
      MARGIN = 1L,
      STATS = sd_draws[draw, ],
      FUN = "*"
    )
    effects <- .bt_random_effect_mvn_group_draws_from_factor(
      factor = factor,
      n_groups = length(groups)
    )
    output[rows, draw] <- rowSums(
      model_matrix[rows, , drop = FALSE] *
        effects[group_index, , drop = FALSE]
    )
  }

  output
}

.bt_random_effect_group_contribution_sample_known_group_covariance <- function(
    random_term,
    model_matrix,
    group_map,
    rows,
    posterior,
    prior_list){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  if(ncol(model_matrix) != 1L){
    stop(
      "Random-effect prediction for block '", random_term$block_name,
      "' with known group covariance supports one random-effect column only.",
      call. = FALSE
    )
  }

  group_covariance <- .bt_random_effect_known_group_covariance(
    random_term,
    context = "Random-effect prediction"
  )
  if(any(group_map[rows] > length(group_covariance$levels))){
    stop(
      "Random-effect prediction for block '", random_term$block_name,
      "' cannot sample new levels with known group covariance.",
      call. = FALSE
    )
  }
  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = 1L,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(sd_draws)){
    .bt_random_effect_marginal_covariance_missing_sd_stop(
      random_term = random_term,
      n_columns = 1L
    )
  }
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = sd_draws,
    n_draws = n_draws,
    n_columns = 1L,
    label = "SD",
    random_term = random_term,
    nonnegative = TRUE,
    context = "Random-effect prediction"
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  kernel <- group_covariance$kernel[groups, groups, drop = FALSE]
  factor <- t(chol(kernel))
  group_effects <- .bt_random_effect_mvn_group_draws_from_factor(
    factor = factor,
    n_groups = n_draws
  )
  group_effects <- group_effects * sd_draws[, 1L]
  output[rows, ] <- t(group_effects[, group_index, drop = FALSE]) *
    matrix(
      model_matrix[rows, 1L],
      nrow = length(rows),
      ncol = n_draws
    )

  output
}

# Sample scalar-structured new groups from only their requested columns.
.bt_random_effect_group_contribution_sample_structured <- function(
    random_term,
    model_matrix,
    group_map,
    rows,
    posterior,
    sd_draws){

  n_draws   <- nrow(posterior)
  n_rows    <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output    <- matrix(0, nrow = n_rows, ncol = n_draws)
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect prediction metadata"
  )
  rho_draws <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "error",
    out_of_support = "error",
    context = "Random-effect prediction metadata"
  )
  column_coordinates <- if(identical(structure, "car")){
    correlation <- .bt_random_effect_correlation_metadata(
      random_term = random_term,
      structure = structure,
      context = "Random-effect prediction metadata"
    )
    .bt_random_effect_car_time_values(
      random_term = random_term,
      correlation = correlation,
      n_columns = n_columns,
      context = "Random-effect prediction metadata"
    )
  }else{
    seq_len(n_columns)
  }
  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = n_columns,
    column_coordinates = column_coordinates
  )
  for(rho in unique(rho_draws)){
    .bt_random_effect_structured_local_check_rho(
      structure = structure,
      rho = rho,
      global_n_columns = n_columns
    )
  }

  groups <- sort(unique(group_map[rows]))
  group_rows <- lapply(groups, function(group){
    rows[group_map[rows] == group]
  })
  group_columns <- lapply(group_rows, function(group_rows_i){
    columns <- unname(which(colSums(
      model_matrix[group_rows_i, , drop = FALSE] != 0
    ) > 0L))
    columns[order(coordinates[columns], columns)]
  })
  nonempty <- lengths(group_columns) > 0L
  if(!any(nonempty)){
    return(output)
  }

  group_rows    <- group_rows[nonempty]
  group_columns <- group_columns[nonempty]
  active_groups <- groups[nonempty]
  subset_keys   <- vapply(group_columns, paste, collapse = ",", character(1))
  unique_keys   <- unique(subset_keys)
  subset_groups <- lapply(unique_keys, function(key) which(subset_keys == key))
  subset_columns <- lapply(subset_groups, function(group_index){
    group_columns[[group_index[1L]]]
  })

  for(draw in seq_len(n_draws)){
    for(subset in seq_along(subset_columns)){
      columns     <- subset_columns[[subset]]
      group_index <- subset_groups[[subset]]
      z <- matrix(
        stats::rnorm(length(group_index) * length(columns)),
        nrow = length(group_index),
        ncol = length(columns)
      )
      effects <- matrix(NA_real_, nrow = length(group_index),
                        ncol = length(columns))
      for(index in seq_along(group_index)){
        group_position <- group_index[index]
        effects[index, ] <- .bt_random_effect_prediction_structured_subset_transform(
          structure = structure,
          columns = columns,
          latent = z[index, ],
          rho = rho_draws[draw],
          coordinates = coordinates,
          context = paste0(
            "Random-effect new-group prediction",
            .bt_random_effect_metadata_block_detail(random_term),
            ", posterior draw ", draw,
            ", group ", active_groups[group_position]
          )
        )
      }
      effects <- sweep(
        effects,
        MARGIN = 2L,
        STATS = sd_draws[draw, columns],
        FUN = "*"
      )

      for(index in seq_along(group_index)){
        group_position <- group_index[index]
        rows_i <- group_rows[[group_position]]
        output[rows_i, draw] <- drop(
          model_matrix[rows_i, columns, drop = FALSE] %*%
            effects[index, ]
        )
      }
    }
  }

  output
}

# Sample independent new groups without an identity covariance decomposition.
.bt_random_effect_group_contribution_sample_independent <- function(
    model_matrix,
    group_map,
    rows,
    posterior,
    column_scale_draws = NULL,
    row_scale_draws = NULL,
    draw_scale = NULL){

  n_draws   <- nrow(posterior)
  n_rows    <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output    <- matrix(0, nrow = n_rows, ncol = n_draws)
  if(length(rows) == 0L){
    return(output)
  }

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    effects <- matrix(
      stats::rnorm(length(groups) * n_columns),
      nrow = length(groups),
      ncol = n_columns
    )
    if(!is.null(column_scale_draws)){
      effects <- sweep(
        effects,
        MARGIN = 2L,
        STATS = column_scale_draws[draw, ],
        FUN = "*"
      )
    }
    contribution <- rowSums(
      model_matrix[rows, , drop = FALSE] *
        effects[group_index, , drop = FALSE]
    )
    if(!is.null(row_scale_draws)){
      contribution <- contribution * row_scale_draws[draw, rows]
    }
    if(!is.null(draw_scale)){
      contribution <- contribution * draw_scale[draw]
    }
    output[rows, draw] <- contribution
  }

  output
}

# Apply the exact scalar-structure recurrence to one validated principal subset.
.bt_random_effect_prediction_structured_subset_transform <- function(
    structure,
    columns,
    latent,
    rho,
    coordinates,
    context = NULL){

  .bt_random_effect_structured_subset_transform(
    structure = structure,
    columns = columns,
    latent = latent,
    rho = rho,
    global_n_columns = length(coordinates),
    column_coordinates = coordinates,
    context = context
  )
}

.bt_random_effect_group_contribution_sample_row_indexed <- function(
    random_term,
    model_matrix,
    group_map,
    rows,
    posterior,
    prior_list,
    source_data,
    prediction_rows = NULL){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)

  source_draws <- .bt_random_effect_row_indexed_source_draws(
    random_term = random_term,
    n_rows = n_rows,
    posterior = posterior,
    data = source_data,
    prediction_rows = prediction_rows,
    context = "Prediction"
  )
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = source_draws,
    n_draws = n_draws,
    n_columns = n_rows,
    label = "row-indexed SD source",
    random_term = random_term,
    nonnegative = TRUE,
    context = "Random-effect prediction"
  )
  column_allocation <- .bt_random_effect_row_indexed_column_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    n_columns = n_columns
  )
  if(is.null(column_allocation)){
    allocation <- .bt_random_effect_row_indexed_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list
    )
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = matrix(allocation, ncol = 1L),
      n_draws = n_draws,
      n_columns = 1L,
      label = "row-indexed SD allocation",
      random_term = random_term,
      nonnegative = TRUE,
      context = "Random-effect prediction"
    )
  }else{
    allocation <- NULL
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = column_allocation,
      n_draws = n_draws,
      n_columns = n_columns,
      label = "row-indexed column SD allocation",
      random_term = random_term,
      nonnegative = TRUE,
      context = "Random-effect prediction"
    )
  }
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect prediction metadata"
  )
  if(n_columns == 1L || structure %in% c("diag", "id")){
    if(is.null(column_allocation)){
      return(.bt_random_effect_group_contribution_sample_independent(
        model_matrix = model_matrix,
        group_map = group_map,
        rows = rows,
        posterior = posterior,
        row_scale_draws = source_draws,
        draw_scale = allocation
      ))
    }
    return(.bt_random_effect_group_contribution_sample_independent(
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      column_scale_draws = column_allocation,
      row_scale_draws = source_draws
    ))
  }
  rho_draws <- if(n_columns > 1L &&
                    structure %in% c("cs", "hcs", "ar1", "car", "har")){
    .bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior,
      missing = "null",
      context = "Random-effect prediction metadata"
    )
  }else{
    NULL
  }
  if(n_columns > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har") &&
     !is.null(rho_draws)){
    unit_contribution <- .bt_random_effect_group_contribution_sample_structured(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      sd_draws = matrix(1, nrow = n_draws, ncol = n_columns)
    )
    if(is.null(column_allocation)){
      allocation_matrix <- matrix(
        allocation,
        nrow = n_rows,
        ncol = n_draws,
        byrow = TRUE
      )
    }else{
      row_column <- .bt_random_effect_structured_indicator_columns(
        model_matrix,
        context = "Prediction for a row-indexed scalar-structured random effect"
      )
      allocation_matrix <- t(column_allocation[, row_column, drop = FALSE])
    }
    return(unit_contribution * t(source_draws) * allocation_matrix)
  }
  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )
  if(is.null(cholesky)){
    .bt_random_effect_marginal_covariance_missing_correlation_stop(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
  }
  .bt_random_effect_marginal_covariance_validate_correlation_cholesky(
    cholesky = cholesky,
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    effects <- .bt_random_effect_mvn_group_draws_from_factor(
      factor = cholesky[draw, , ],
      n_groups = length(groups)
    )
    Z <- model_matrix[rows, , drop = FALSE] *
      matrix(source_draws[draw, rows], nrow = length(rows), ncol = n_columns)
    if(is.null(column_allocation)){
      Z <- Z * allocation[draw]
    }else{
      Z <- Z * matrix(
        column_allocation[draw, ],
        nrow = length(rows),
        ncol = n_columns,
        byrow = TRUE
      )
    }
    output[rows, draw] <- rowSums(Z * effects[group_index, , drop = FALSE])
  }

  output
}

.bt_random_effect_mvn_group_draws_from_factor <- function(factor, n_groups){

  if(n_groups == 0L){
    return(matrix(numeric(), nrow = 0L, ncol = ncol(factor)))
  }
  if(!is.matrix(factor) || !is.numeric(factor) ||
     nrow(factor) != ncol(factor) || any(!is.finite(factor))){
    stop("Random-effect prediction factor must be a finite numeric square matrix.",
         call. = FALSE)
  }
  z <- matrix(
    stats::rnorm(n_groups * ncol(factor)),
    nrow = n_groups,
    ncol = ncol(factor)
  )

  unname(z %*% t(factor))
}

.bt_random_effect_prediction_data <- function(random_term, data,
                                              group_data = data,
                                              allow_new_groups = FALSE,
                                              context = "JAGS_evaluate_formula()"){

  prediction_data <- data
  random_structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect prediction metadata"
  )
  prediction_data <- .bt_random_effect_prediction_structured_index_data(
    random_term,
    prediction_data
  )
  factor_levels <- random_term$xlevels
  if(!is.null(factor_levels) && length(factor_levels) > 0L){
    for(factor_name in names(factor_levels)){
      .bt_validate_categorical_level_names(
        factor_levels[[factor_name]],
        factor_name,
        context = "Fitted random-effect factor metadata"
      )
      if(!factor_name %in% names(prediction_data)){
        stop(
          "The '", factor_name,
          "' predictor needed for random-effect prediction is missing in the data.",
          call. = FALSE
        )
      }
      .bt_validate_categorical_values(
        prediction_data[[factor_name]],
        factor_name,
        context = "Random-effect factor predictor"
      )
      if(is.factor(prediction_data[[factor_name]])){
        observed_levels <- unique(as.character(prediction_data[[factor_name]]))
        if(!all(observed_levels %in% factor_levels[[factor_name]])){
          stop(
            "Levels specified in the '", factor_name,
            "' factor variable do not match the levels used for model specification.",
            call. = FALSE
          )
        }
        prediction_data[[factor_name]] <- factor(
          as.character(prediction_data[[factor_name]]),
          levels = factor_levels[[factor_name]]
        )
      }else if(all(unique(prediction_data[[factor_name]]) %in% factor_levels[[factor_name]])){
        prediction_data[[factor_name]] <- factor(
          prediction_data[[factor_name]],
          levels = factor_levels[[factor_name]]
        )
      }else{
        stop(
          "Levels specified in the '", factor_name,
          "' factor variable do not match the levels used for model specification.",
          call. = FALSE
        )
      }
    }
  }

  contrast_owner <- random_term$contrast_owner
  contrast_matrices <- random_term$contrast_matrices
  if(is.null(contrast_owner) ||
     !contrast_owner %in% c("random_block", "structure") ||
     is.null(contrast_matrices)){
    stop(
      "Random-effect prediction metadata for block '",
      random_term$block_name,
      "' is missing its owner-scoped concrete factor basis. Refit the model ",
      "with this version of BayesTools.",
      call. = FALSE
    )
  }
  if(identical(contrast_owner, "random_block")){
    for(factor_name in names(contrast_matrices)){
      if(factor_name %in% names(prediction_data) &&
         is.factor(prediction_data[[factor_name]])){
        attr(prediction_data[[factor_name]], "contrasts") <-
          contrast_matrices[[factor_name]]
      }
    }
  }

  .bt_random_effect_check_memory(
    estimate = .bt_random_effect_design_memory_estimate(
      n_rows = nrow(prediction_data),
      n_columns = random_term$n_columns,
      n_groups = random_term$n_groups +
        if(isTRUE(allow_new_groups)) nrow(prediction_data) else 0L,
      monitor_policy = random_term$monitor,
      structure = random_structure,
      compile_mode = random_term$compile_mode
    ),
    block_name = random_term$block_name,
    alternative = paste0(
      "Reduce the number of prediction rows or random-effect columns, or ",
      "request fewer random-effect blocks. Raise the option (or set it to ",
      "Inf) only after verifying the operation's memory budget."
    )
  )

  random_design <- .bt_random_effect_design_matrix(
    random_term$term_formula,
    prediction_data,
    preserve_no_intercept_contrasts = !random_structure %in% c("cs", "hcs", "ar1", "car", "har"),
    structure = random_structure,
    car_time_values = if(identical(random_structure, "car")) random_term$car$time_values else NULL,
    block_name = random_term$block_name
  )
  model_matrix <- random_design$model_matrix
  attr(model_matrix, "contrasts") <- random_term$contrasts
  colnames(model_matrix) <- gsub(":", "__xXx__", colnames(model_matrix))

  if(!identical(colnames(model_matrix), random_term$column_names)){
    stop(
      "Random-effect design columns for block '", random_term$block_name,
      "' do not match the fitted formula.",
      call. = FALSE
    )
  }

  grouping_observations <- .bt_random_group_observations(
    random_term,
    group_data
  )
  if(nrow(grouping_observations$tuple_values) != nrow(model_matrix)){
    stop(
      "Random-effect grouping data for block '", random_term$block_name,
      "' must have one value per prediction row.",
      call. = FALSE
    )
  }
  if(is.null(random_term$group_components) ||
     is.null(random_term$group_tuple_keys) ||
     is.null(random_term$group_tuple_index)){
    stop(
      "Random-effect prediction metadata for block '", random_term$block_name,
      "' is missing the fitted grouping tuple map. Refit the model with this ",
      "version of BayesTools.",
      call. = FALSE
    )
  }
  if(!identical(
    unname(grouping_observations$component_names),
    unname(random_term$group_components)
  )){
    stop(
      "Random-effect grouping components for block '", random_term$block_name,
      "' do not match the fitted formula.",
      call. = FALSE
    )
  }
  if(!is.null(random_term$group_component_levels)){
    for(component_name in names(random_term$group_component_levels)){
      .bt_validate_categorical_level_names(
        random_term$group_component_levels[[component_name]],
        component_name,
        context = "Fitted random-effect grouping metadata"
      )
    }
  }

  group_levels <- random_term$group_levels
  group_tuple_keys <- random_term$group_tuple_keys
  group_tuple_index <- random_term$group_tuple_index
  group_map <- unname(group_tuple_index[grouping_observations$tuple_keys])
  if(any(is.na(group_map))){
    new_keys <- unique(grouping_observations$tuple_keys[is.na(group_map)])
    new_rows <- match(new_keys, grouping_observations$tuple_keys)
    new_groups <- grouping_observations$display_labels[new_rows]
    if(.bt_random_effect_has_known_group_covariance(random_term)){
      stop(
        "New random-effect level(s) for block '", random_term$block_name,
        "' are not supported with known group covariance: ",
        paste(new_groups, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if(isTRUE(allow_new_groups)){
      new_groups <- .bt_random_group_unique_labels(
        new_groups,
        new_keys,
        existing = group_levels
      )
      group_levels <- c(group_levels, new_groups)
      new_indices <- seq.int(
        length(group_tuple_keys) + 1L,
        length(group_tuple_keys) + length(new_keys)
      )
      group_tuple_keys <- c(group_tuple_keys, new_keys)
      group_tuple_index <- c(
        group_tuple_index,
        stats::setNames(new_indices, new_keys)
      )
      group_map <- unname(
        group_tuple_index[grouping_observations$tuple_keys]
      )
      return(list(
        model_matrix = model_matrix,
        group_map = group_map,
        group_levels = group_levels,
        group_tuple_keys = group_tuple_keys,
        group_tuple_index = group_tuple_index
      ))
    }
    stop(
      "New random-effect level(s) for block '", random_term$block_name,
      "' are not supported by ", context, ": ",
      paste(new_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    model_matrix = model_matrix,
    group_map = group_map,
    group_levels = group_levels,
    group_tuple_keys = group_tuple_keys,
    group_tuple_index = group_tuple_index
  )
}

.bt_random_effect_prediction_structured_index_data <- function(random_term,
                                                               data){

  index <- random_term$structured_index
  if(is.null(index)){
    return(data)
  }

  missing_variables <- index$variables[!index$variables %in% names(data)]
  if(length(missing_variables) > 0L){
    stop(
      "The ",
      paste0("'", missing_variables, "'", collapse = ", "),
      " structured random-effect index variable is missing in the data.",
      call. = FALSE
    )
  }

  data[[index$name]] <- .bt_random_effect_structured_index_values(
    data,
    index$variables
  )
  data
}
