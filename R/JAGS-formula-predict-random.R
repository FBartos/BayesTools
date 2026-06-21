.bt_apply_formula_scale_to_data <- function(fit, parameter, data,
                                            predictors_type){

  continuous_predictors <- names(predictors_type[predictors_type == "continuous"])

  formula_scale <- attr(fit, "formula_scale")
  if(is.null(formula_scale)){
    return(data)
  }

  param_scale <- formula_scale[[parameter]]
  if(is.null(param_scale)){
    return(data)
  }

  scaled_predictors <- sub(paste0("^", parameter, "_"), "", names(param_scale))
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
                                                          new_levels = NULL){

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
    parameter = parameter
  )

  fixed_formula <- .remove_random_effects(formula)
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
      new_levels = new_levels
    )
  }

  output
}

.bt_JAGS_evaluate_formula_random_effect_terms <- function(fitted_design,
                                                          requested_terms,
                                                          formula_target,
                                                          blocks,
                                                          parameter){

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

  if(!identical(formula_target, "conditional")){
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
  marginalized <- selected_names[modes == "marginalized"]
  if(length(marginalized) > 0L){
    stop(
      "JAGS_evaluate_formula() cannot use formula_target = \"conditional\" for ",
      "random-effect block(s) compiled as marginalized: ",
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
                                                new_levels = NULL){

  new_levels <- .bt_random_effect_new_levels_policy(random_term, new_levels)
  prediction <- .bt_random_effect_prediction_data(
    random_term,
    data,
    group_data = group_data,
    allow_new_groups = isTRUE(new_levels$allow)
  )
  model_matrix <- prediction$model_matrix
  group_map <- prediction$group_map
  n_rows <- nrow(model_matrix)
  n_draws <- nrow(posterior)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  fitted_n_groups <- length(random_term$group_levels)
  new_row <- group_map > fitted_n_groups

  if(any(!new_row)){
    existing_rows <- which(!new_row)
    output[existing_rows, ] <- .bt_JAGS_evaluate_random_effect_term_existing(
      random_term = random_term,
      model_matrix = model_matrix[existing_rows, , drop = FALSE],
      group_map = group_map[existing_rows],
      posterior = posterior,
      prior_list = prior_list,
      group_data = group_data[existing_rows, , drop = FALSE]
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
                                                          group_data){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
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
                                                        source_data){

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
      source_data = source_data
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
  correlation <- .bt_random_effect_marginal_covariance_correlation_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    G <- matrix(
      correlation[draw, , ],
      nrow = n_columns,
      ncol = n_columns
    ) * tcrossprod(sd_draws[draw, ])
    effects <- .bt_random_effect_mvn_group_draws(
      covariance = G,
      n_groups = length(groups)
    )
    output[rows, draw] <- rowSums(
      model_matrix[rows, , drop = FALSE] *
        effects[group_index, , drop = FALSE]
    )
  }

  output
}

.bt_random_effect_group_contribution_sample_row_indexed <- function(
    random_term,
    model_matrix,
    group_map,
    rows,
    posterior,
    prior_list,
    source_data){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)

  source_draws <- .bt_random_effect_row_indexed_source_draws(
    random_term = random_term,
    n_rows = n_rows,
    posterior = posterior,
    data = source_data,
    context = "Prediction"
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
  }else{
    allocation <- NULL
  }
  correlation <- .bt_random_effect_marginal_covariance_correlation_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    effects <- .bt_random_effect_mvn_group_draws(
      covariance = matrix(correlation[draw, , ], nrow = n_columns, ncol = n_columns),
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

.bt_random_effect_mvn_group_draws <- function(covariance, n_groups){

  if(n_groups == 0L){
    return(matrix(numeric(), nrow = 0L, ncol = ncol(covariance)))
  }
  decomposition <- eigen(covariance, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  transform <- decomposition$vectors %*%
    (sqrt(values) * t(decomposition$vectors))
  z <- matrix(
    stats::rnorm(n_groups * ncol(covariance)),
    nrow = n_groups,
    ncol = ncol(covariance)
  )

  z %*% transform
}

.bt_random_effect_prediction_data <- function(random_term, data,
                                              group_data = data,
                                              allow_new_groups = FALSE){

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
      if(!factor_name %in% names(prediction_data)){
        stop(
          "The '", factor_name,
          "' predictor needed for random-effect prediction is missing in the data.",
          call. = FALSE
        )
      }
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

  contrasts <- random_term$contrasts
  if(!is.null(contrasts) && length(contrasts) > 0L){
    for(factor_name in names(contrasts)){
      if(factor_name %in% names(prediction_data) && is.factor(prediction_data[[factor_name]])){
        stats::contrasts(prediction_data[[factor_name]]) <- contrasts[[factor_name]]
      }
    }
  }

  random_design <- .bt_random_effect_design_matrix(
    random_term$term_formula,
    prediction_data,
    preserve_no_intercept_contrasts = !random_structure %in% c("cs", "hcs", "ar1", "car", "har"),
    structure = random_structure,
    car_time_values = if(identical(random_structure, "car")) random_term$car$time_values else NULL,
    block_name = random_term$block_name
  )
  model_matrix <- random_design$model_matrix
  colnames(model_matrix) <- gsub(":", "__xXx__", colnames(model_matrix))

  if(!identical(colnames(model_matrix), random_term$column_names)){
    stop(
      "Random-effect design columns for block '", random_term$block_name,
      "' do not match the fitted formula.",
      call. = FALSE
    )
  }

  grouping_values <- .bt_random_group_values(random_term, group_data)
  if(length(grouping_values) != nrow(model_matrix)){
    stop(
      "Random-effect grouping data for block '", random_term$block_name,
      "' must have one value per prediction row.",
      call. = FALSE
    )
  }
  group_levels <- random_term$group_levels
  group_map <- match(as.character(grouping_values), group_levels)
  if(any(is.na(group_map))){
    new_groups <- unique(as.character(grouping_values)[is.na(group_map)])
    if(isTRUE(allow_new_groups)){
      group_levels <- c(group_levels, new_groups)
      group_map <- match(as.character(grouping_values), group_levels)
      return(list(
        model_matrix = model_matrix,
        group_map = group_map,
        group_levels = group_levels
      ))
    }
    stop(
      "New random-effect level(s) for block '", random_term$block_name,
      "' are not supported by JAGS_evaluate_formula(): ",
      paste(new_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    model_matrix = model_matrix,
    group_map = group_map,
    group_levels = group_levels
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
