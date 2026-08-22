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

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect posterior reconstruction metadata"
  )
  if(n_columns > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har")){
    sd_draws <- .bt_random_effect_sd_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list,
      n_columns = n_columns
    )
    if(is.null(sd_draws)){
      return(NULL)
    }
    structured_contribution <- .bt_random_effect_structured_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      scale_draws = sd_draws
    )
    if(!is.null(structured_contribution)){
      return(structured_contribution)
    }
  }

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

  if(structure %in% c("diag", "id")){
    return(.bt_random_effect_independent_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      column_scale_draws = sd_draws
    ))
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

.bt_random_effect_structured_contribution_from_latent <- function(
    random_term, model_matrix, group_map, posterior, scale_draws){

  if(inherits(
    random_term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )){
    return(.bt_random_effect_structured_local_contribution(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      scale_draws = scale_draws
    ))
  }

  .bt_random_effect_structured_dense_contribution(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    posterior = posterior,
    scale_draws = scale_draws
  )
}

.bt_random_effect_structured_dense_contribution <- function(
    random_term, model_matrix, group_map, posterior, scale_draws){

  n_columns <- ncol(model_matrix)
  n_groups  <- length(random_term$group_levels)
  z_names   <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(!all(as.vector(z_names) %in% colnames(posterior))){
    return(NULL)
  }
  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior
  )
  if(is.null(rho)){
    return(NULL)
  }

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect posterior reconstruction metadata"
  )
  column_coordinates <- if(identical(structure, "car")){
    correlation <- .bt_random_effect_correlation_metadata(
      random_term = random_term,
      structure = structure,
      context = "Random-effect posterior reconstruction metadata"
    )
    .bt_random_effect_car_time_values(
      random_term = random_term,
      correlation = correlation,
      n_columns = n_columns,
      context = "Random-effect posterior reconstruction metadata"
    )
  }else{
    seq_len(n_columns)
  }
  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = n_columns,
    column_coordinates = column_coordinates
  )
  row_column <- .bt_random_effect_structured_indicator_columns(
    model_matrix,
    context = "Prediction for a scalar-structured random effect"
  )
  n_draws <- nrow(posterior)
  latent <- do.call(cbind, lapply(seq_len(n_columns), function(column){
    as.vector(unname(posterior[, z_names[, column], drop = FALSE]))
  }))
  transition_context <- function(index){
    draw  <- (index - 1L) %% n_draws + 1L
    group <- (index - 1L) %/% n_draws + 1L
    paste0(
      "Random-effect posterior reconstruction",
      .bt_random_effect_metadata_block_detail(random_term),
      ", posterior draw ", draw, ", group ", group
    )
  }
  unit <- .bt_random_effect_structured_subset_transform_draws(
    structure = structure,
    columns = seq_len(n_columns),
    latent = latent,
    rho = rep(rho, times = n_groups),
    global_n_columns = n_columns,
    column_coordinates = coordinates,
    context = transition_context
  )
  out <- matrix(NA_real_, nrow = nrow(model_matrix), ncol = n_draws)
  for(column in seq_len(n_columns)){
    rows <- which(row_column == column)
    if(length(rows) == 0L){
      next
    }
    coefficient <- matrix(
      unit[, column] * rep(scale_draws[, column], times = n_groups),
      nrow = n_draws,
      ncol = n_groups
    )
    out[rows, ] <- t(coefficient[, group_map[rows], drop = FALSE])
  }

  out
}

.bt_random_effect_structured_indicator_columns <- function(model_matrix,
                                                           context){

  nonzero <- model_matrix != 0
  if(any(rowSums(nonzero) != 1L)){
    stop(
      context, " requires one index level per row.",
      call. = FALSE
    )
  }
  row_column <- max.col(nonzero, ties.method = "first")
  selected <- model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)]
  if(any(selected != 1)){
    stop(
      context, " requires unit index indicators.",
      call. = FALSE
    )
  }

  row_column
}

.bt_random_effect_structured_local_contribution <- function(random_term,
                                                           model_matrix,
                                                           group_map,
                                                           posterior,
                                                           scale_draws){

  layout  <- random_term$latent_layout
  z_names <- layout$node_names
  if(!all(z_names %in% colnames(posterior))){
    return(NULL)
  }
  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior
  )
  if(is.null(rho)){
    return(NULL)
  }

  out <- matrix(
    NA_real_,
    nrow = nrow(model_matrix),
    ncol = nrow(posterior)
  )
  row_column <- .bt_random_effect_structured_indicator_columns(
    model_matrix,
    context = "Prediction for a group-local structured random effect"
  )
  requested_columns <- lapply(seq_len(layout$n_groups), function(group){
    unique(row_column[group_map == group])
  })
  has_missing_columns <- any(vapply(seq_len(layout$n_groups), function(group){
    length(setdiff(requested_columns[[group]], layout$group_columns[[group]])) > 0L
  }, logical(1)))
  if(!has_missing_columns){
    for(group in seq_len(layout$n_groups)){
      rows <- which(group_map == group)
      if(length(rows) == 0L){
        next
      }
      observed_columns <- layout$group_columns[[group]]
      names <- .bt_random_effect_structured_local_node_names(
        parameter_stem = random_term$parameter_stem,
        group = rep(group, length(observed_columns)),
        column = observed_columns
      )
      transition_context <- function(draw){
        paste0(
          "Random-effect posterior reconstruction",
          .bt_random_effect_metadata_block_detail(random_term),
          ", posterior draw ", draw, ", group ", group
        )
      }
      unit <- .bt_random_effect_structured_subset_transform_draws(
        structure = layout$structure,
        columns = observed_columns,
        latent = unname(posterior[, names, drop = FALSE]),
        rho = rho,
        global_n_columns = layout$global_n_columns,
        column_coordinates = layout$column_coordinates,
        context = transition_context
      )
      observed_coefficients <- unit *
        scale_draws[, observed_columns, drop = FALSE]
      requested_index <- match(row_column[rows], observed_columns)
      out[rows, ] <- t(observed_coefficients[, requested_index, drop = FALSE])
    }
    return(out)
  }

  for(draw in seq_len(nrow(posterior))){
    group_coefficients <- vector("list", layout$n_groups)
    for(group in seq_len(layout$n_groups)){
      observed_columns <- layout$group_columns[[group]]
      transition_context <- paste0(
        "Random-effect posterior reconstruction",
        .bt_random_effect_metadata_block_detail(random_term),
        ", posterior draw ", draw, ", group ", group
      )
      names <- .bt_random_effect_structured_local_node_names(
        parameter_stem = random_term$parameter_stem,
        group = rep(group, length(observed_columns)),
        column = observed_columns
      )
      unit <- .bt_random_effect_structured_subset_transform(
        structure = layout$structure,
        columns = observed_columns,
        latent = as.numeric(posterior[draw, names]),
        rho = rho[draw],
        global_n_columns = layout$global_n_columns,
        column_coordinates = layout$column_coordinates,
        context = transition_context
      )
      observed_coefficients <- unit * scale_draws[draw, observed_columns]
      requested_columns <- unique(row_column[group_map == group])
      missing_columns <- setdiff(requested_columns, observed_columns)
      values <- stats::setNames(observed_coefficients, observed_columns)
      if(length(missing_columns) > 0L){
        values <- c(values, .bt_random_effect_structured_local_conditional_draw(
          structure = layout$structure,
          observed_columns = observed_columns,
          missing_columns = missing_columns,
          observed_unit = unit,
          sd = scale_draws[draw, ],
          rho = rho[draw],
          global_n_columns = layout$global_n_columns,
          column_coordinates = layout$column_coordinates,
          context = transition_context
        ))
      }
      group_coefficients[[group]] <- values
    }
    out[, draw] <- vapply(seq_len(nrow(model_matrix)), function(row){
      group_coefficients[[group_map[row]]][as.character(row_column[row])]
    }, numeric(1))
  }

  out
}

.bt_random_effect_structured_local_conditional_draw <- function(
    structure, observed_columns, missing_columns, observed_unit,
    sd, rho, global_n_columns, column_coordinates, context = NULL){

  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  check_int(global_n_columns, "global_n_columns", lower = 1, allow_NA = FALSE)
  columns <- c(observed_columns, missing_columns)
  if(!is.numeric(observed_columns) || length(observed_columns) < 1L ||
     !is.numeric(missing_columns) || length(missing_columns) < 1L ||
     any(is.na(columns)) || any(columns != as.integer(columns)) ||
     any(columns < 1L) || any(columns > global_n_columns) ||
     anyDuplicated(columns)){
    stop(
      "Observed and missing structured columns must be disjoint positive integer indices.",
      call. = FALSE
    )
  }
  if(!is.numeric(observed_unit) ||
     length(observed_unit) != length(observed_columns) ||
     any(!is.finite(observed_unit))){
    stop("'observed_unit' must contain one finite value per observed column.",
         call. = FALSE)
  }
  if(!is.numeric(sd) || length(sd) != global_n_columns ||
     any(!is.finite(sd)) || any(sd < 0)){
    stop("'sd' must contain one finite non-negative value per global column.",
         call. = FALSE)
  }
  .bt_random_effect_structured_local_check_rho(
    structure = structure,
    rho = rho,
    global_n_columns = global_n_columns
  )

  innovations <- stats::rnorm(length(missing_columns))
  if(structure %in% c("cs", "hcs")){
    draw <- .bt_random_effect_cs_conditional_unit_draw(
      observed_unit = observed_unit,
      rho = rho,
      n_missing = length(missing_columns),
      innovations = innovations
    )
  }else{
    coordinates <- .bt_random_effect_structured_local_coordinates(
      structure = structure,
      n_columns = global_n_columns,
      column_coordinates = column_coordinates
    )
    draw <- .bt_random_effect_markov_conditional_unit_draw(
      observed_columns = observed_columns,
      missing_columns = missing_columns,
      observed_unit = observed_unit,
      rho = rho,
      column_coordinates = coordinates,
      innovations = innovations,
      context = context
    )
  }

  stats::setNames(draw * sd[missing_columns], missing_columns)
}

# Sample a CS/HCS conditional block using its diagonal-plus-rank-one form.
.bt_random_effect_cs_conditional_unit_draw <- function(
    observed_unit, rho, n_missing, innovations){

  check_int(n_missing, "n_missing", lower = 1, allow_NA = FALSE)
  if(!is.numeric(innovations) || length(innovations) != n_missing ||
     any(!is.finite(innovations))){
    stop("'innovations' must contain one finite value per missing column.",
         call. = FALSE)
  }

  n_observed <- length(observed_unit)
  residual_variance <- 1 - rho
  observed_denominator <- 1 + (n_observed - 1L) * rho
  conditional_common <- rho * residual_variance / observed_denominator
  conditional_diagonal <- residual_variance + conditional_common
  conditional_mean <- rho * sum(observed_unit) / observed_denominator
  if(!is.finite(conditional_diagonal) || conditional_diagonal < 0){
    stop("Conditional CS/HCS random-effect variance is not non-negative.",
         call. = FALSE)
  }
  if(n_missing == 1L){
    return(conditional_mean + sqrt(conditional_diagonal) * innovations)
  }

  conditional_rho <- conditional_common / conditional_diagonal
  standardized <- .bt_random_effect_structured_subset_transform(
    structure = "cs",
    columns = seq_len(n_missing),
    latent = innovations,
    rho = conditional_rho,
    global_n_columns = n_missing
  )
  conditional_mean + sqrt(conditional_diagonal) * standardized
}

# Sample AR1/HAR/CAR missing coordinates through exact Gaussian Markov bridges.
.bt_random_effect_markov_conditional_unit_draw <- function(
    observed_columns, missing_columns, observed_unit, rho,
    column_coordinates, innovations, context = NULL){

  columns <- c(observed_columns, missing_columns)
  coordinates <- column_coordinates[columns]
  order_index <- order(coordinates, columns)
  ordered_columns <- columns[order_index]
  ordered_coordinates <- coordinates[order_index]
  if(any(diff(ordered_coordinates) <= 0)){
    stop("Structured conditional columns must have unique ordered coordinates.",
         call. = FALSE)
  }

  is_observed <- ordered_columns %in% observed_columns
  values <- rep(NA_real_, length(ordered_columns))
  values[is_observed] <- observed_unit[
    match(ordered_columns[is_observed], observed_columns)
  ]
  observed_positions <- which(is_observed)
  innovations_by_column <- stats::setNames(innovations, missing_columns)

  first_observed <- observed_positions[1L]
  if(first_observed > 1L){
    for(position in seq.int(first_observed - 1L, 1L)){
      transition <- .bt_random_effect_markov_transition(
        rho = rho,
        left_coordinate = ordered_coordinates[position],
        right_coordinate = ordered_coordinates[position + 1L],
        context = context
      )
      values[position] <- transition$phi * values[position + 1L] +
        sqrt(transition$innovation_variance) *
        innovations_by_column[[as.character(ordered_columns[position])]]
    }
  }

  if(length(observed_positions) > 1L){
    for(interval in seq_len(length(observed_positions) - 1L)){
      left <- observed_positions[interval]
      right <- observed_positions[interval + 1L]
      if(right - left <= 1L){
        next
      }
      current <- left
      for(position in seq.int(left + 1L, right - 1L)){
        left_transition <- .bt_random_effect_markov_transition(
          rho = rho,
          left_coordinate = ordered_coordinates[current],
          right_coordinate = ordered_coordinates[position],
          context = context
        )
        right_transition <- .bt_random_effect_markov_transition(
          rho = rho,
          left_coordinate = ordered_coordinates[position],
          right_coordinate = ordered_coordinates[right],
          context = context
        )
        span_transition <- .bt_random_effect_markov_transition(
          rho = rho,
          left_coordinate = ordered_coordinates[current],
          right_coordinate = ordered_coordinates[right],
          context = context
        )
        right_ratio <- right_transition$innovation_variance /
          span_transition$innovation_variance
        left_ratio <- left_transition$innovation_variance /
          span_transition$innovation_variance
        conditional_mean <-
          left_transition$phi * right_ratio * values[current] +
          right_transition$phi * left_ratio * values[right]
        conditional_variance <-
          left_transition$innovation_variance * right_ratio
        if(!is.finite(conditional_variance) || conditional_variance <= 0){
          bridge_values <- format(
            c(
              ordered_coordinates[current],
              ordered_coordinates[position],
              ordered_coordinates[right],
              ordered_coordinates[position] - ordered_coordinates[current],
              ordered_coordinates[right] - ordered_coordinates[position],
              ordered_coordinates[right] - ordered_coordinates[current],
              rho
            ),
            digits = 17L,
            scientific = TRUE,
            trim = TRUE
          )
          error_label <- .bt_random_effect_markov_error_label(
            context,
            operation = "bridge"
          )
          stop(
            error_label, " at coordinate ", bridge_values[2L],
            " between coordinates ", bridge_values[1L], " and ",
            bridge_values[3L],
            " has a non-positive or non-finite conditional variance ",
            "(rho = ", bridge_values[7L],
            ", left gap = ", bridge_values[4L],
            ", right gap = ", bridge_values[5L],
            ", span gap = ", bridge_values[6L],
            "). The requested coordinate/time resolution is not ",
            "representable for this rho.",
            call. = FALSE
          )
        }
        values[position] <- conditional_mean +
          sqrt(conditional_variance) *
          innovations_by_column[[as.character(ordered_columns[position])]]
        current <- position
      }
    }
  }

  last_observed <- observed_positions[length(observed_positions)]
  if(last_observed < length(ordered_columns)){
    for(position in seq.int(last_observed + 1L, length(ordered_columns))){
      transition <- .bt_random_effect_markov_transition(
        rho = rho,
        left_coordinate = ordered_coordinates[position - 1L],
        right_coordinate = ordered_coordinates[position],
        context = context
      )
      values[position] <- transition$phi * values[position - 1L] +
        sqrt(transition$innovation_variance) *
        innovations_by_column[[as.character(ordered_columns[position])]]
    }
  }

  values[match(missing_columns, ordered_columns)]
}

.bt_random_effect_row_indexed_contribution_from_latent <- function(
    random_term, model_matrix, group_map, posterior, prior_list,
    data = NULL, parameters = NULL, prediction_rows = NULL,
    context = "Prediction"){

  source_draws <- .bt_random_effect_row_indexed_source_draws(
    random_term = random_term,
    n_rows = nrow(model_matrix),
    posterior = posterior,
    data = data,
    parameters = parameters,
    prediction_rows = prediction_rows,
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
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect posterior reconstruction metadata"
  )
  if(structure %in% c("diag", "id")){
    if(!is.null(column_allocation_draws)){
      if(any(!is.finite(column_allocation_draws) |
             column_allocation_draws < 0)){
        stop(
          context, " with row-indexed external SD source '",
          .bt_random_effect_external_sd_source_label(random_term),
          "' returned non-finite or negative allocation values.",
          call. = FALSE
        )
      }
      contribution <- .bt_random_effect_independent_contribution_from_latent(
        random_term = random_term,
        model_matrix = model_matrix,
        group_map = group_map,
        posterior = posterior,
        column_scale_draws = column_allocation_draws,
        row_scale_draws = source_draws
      )
      if(is.null(contribution)){
        .bt_random_effect_row_indexed_reconstruction_missing_stop(random_term)
      }
      return(contribution)
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
    contribution <- .bt_random_effect_independent_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      row_scale_draws = source_draws,
      draw_scale = allocation_draws
    )
    if(is.null(contribution)){
      .bt_random_effect_row_indexed_reconstruction_missing_stop(random_term)
    }
    return(contribution)
  }
  if(ncol(model_matrix) > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har")){
    if(!is.null(column_allocation_draws) &&
       any(!is.finite(column_allocation_draws) |
           column_allocation_draws < 0)){
      stop(
        context, " with row-indexed external SD source '",
        .bt_random_effect_external_sd_source_label(random_term),
        "' returned non-finite or negative allocation values.",
        call. = FALSE
      )
    }
    allocation_draws <- if(is.null(column_allocation_draws)){
      .bt_random_effect_row_indexed_allocation_draws(
        random_term = random_term,
        posterior = posterior,
        prior_list = prior_list
      )
    }else{
      NULL
    }
    if(!is.null(allocation_draws) &&
       any(!is.finite(allocation_draws) | allocation_draws < 0)){
      stop(
        context, " with row-indexed external SD source '",
        .bt_random_effect_external_sd_source_label(random_term),
        "' returned non-finite or negative allocation values.",
        call. = FALSE
      )
    }
    unit_contribution <- .bt_random_effect_structured_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      scale_draws = matrix(
        1,
        nrow = nrow(posterior),
        ncol = ncol(model_matrix)
      )
    )
    if(!is.null(unit_contribution)){
      if(is.null(column_allocation_draws)){
        allocation_matrix <- matrix(
          allocation_draws,
          nrow = nrow(model_matrix),
          ncol = nrow(posterior),
          byrow = TRUE
        )
      }else{
        row_column <- .bt_random_effect_structured_indicator_columns(
          model_matrix,
          context = "Prediction for a row-indexed scalar-structured random effect"
        )
        allocation_matrix <- t(
          column_allocation_draws[, row_column, drop = FALSE]
        )
      }
      return(unit_contribution * t(source_draws) * allocation_matrix)
    }
  }
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

.bt_random_effect_row_indexed_reconstruction_missing_stop <- function(random_term){

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

# Reconstruct independent random effects without materializing identity matrices.
.bt_random_effect_independent_contribution_from_latent <- function(
    random_term,
    model_matrix,
    group_map,
    posterior,
    column_scale_draws = NULL,
    row_scale_draws = NULL,
    draw_scale = NULL){

  n_draws   <- nrow(posterior)
  n_rows    <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  n_groups  <- length(random_term$group_levels)
  z_names   <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(!all(as.vector(z_names) %in% colnames(posterior))){
    return(NULL)
  }
  if(!is.null(column_scale_draws) &&
     (!is.matrix(column_scale_draws) ||
      !identical(dim(column_scale_draws), c(n_draws, n_columns)))){
    stop("Independent random-effect column scales do not match posterior dimensions.",
         call. = FALSE)
  }
  if(!is.null(row_scale_draws) &&
     (!is.matrix(row_scale_draws) ||
      !identical(dim(row_scale_draws), c(n_draws, n_rows)))){
    stop("Independent random-effect row scales do not match posterior dimensions.",
         call. = FALSE)
  }
  if(!is.null(draw_scale) &&
     (!is.numeric(draw_scale) || length(draw_scale) != n_draws)){
    stop("Independent random-effect draw scales do not match posterior dimensions.",
         call. = FALSE)
  }

  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  for(column in seq_len(n_columns)){
    coefficient <- unname(posterior[, z_names[, column], drop = FALSE])
    if(!is.null(column_scale_draws)){
      coefficient <- coefficient * matrix(
        column_scale_draws[, column],
        nrow = n_draws,
        ncol = n_groups
      )
    }
    contribution <- t(coefficient[, group_map, drop = FALSE]) *
      unname(model_matrix[, column])
    output <- output + contribution
  }
  if(!is.null(row_scale_draws)){
    output <- output * t(unname(row_scale_draws))
  }
  if(!is.null(draw_scale)){
    output <- output * matrix(
      draw_scale,
      nrow = n_rows,
      ncol = n_draws,
      byrow = TRUE
    )
  }
  dimnames(output) <- NULL

  output
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
                                                       prediction_rows = NULL,
                                                       context = "Prediction"){

  source <- .bt_random_effect_row_indexed_source(random_term)
  if(is.null(prediction_rows)){
    prediction_rows <- seq_len(n_rows)
  }
  if(!is.numeric(prediction_rows) || length(prediction_rows) != n_rows ||
     anyNA(prediction_rows) ||
     any(prediction_rows != as.integer(prediction_rows)) ||
     any(prediction_rows < 1L)){
    stop(
      context, " with row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' received invalid prediction-row indices.",
      call. = FALSE
    )
  }
  prediction_rows <- as.integer(prediction_rows)
  source_names <- .bt_parameter_source_row_names(
    source$source,
    max(prediction_rows)
  )[prediction_rows]
  source_values <- .bt_parameter_source_value_draws(
    source = source$source,
    n_rows = n_rows,
    posterior = posterior,
    data = data,
    parameters = parameters,
    context = context
  )
  if(!is.null(source_values)){
    colnames(source_values) <- source_names
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
    paste(prediction_rows[missing], collapse = ", "),
    ". Expected posterior column(s): ",
    paste0("'", source_names[missing][seq_len(min(3L, sum(missing)))], "'",
           collapse = ", "),
    if(sum(missing) > 3L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_random_effect_prediction_fitted_rows <- function(
    random_term,
    n_rows,
    data_supplied,
    fitted_rows = NULL,
    new_row = NULL,
    context = "Prediction"){

  source <- .bt_random_effect_row_indexed_source(random_term)
  if(.bt_parameter_source_has_values(source$source)){
    return(NULL)
  }

  if(isTRUE(data_supplied) && is.null(fitted_rows)){
    stop(
      context, " with posterior-indexed row source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' for block '", random_term$block_name,
      "' requires an explicit 'fitted_rows' mapping whenever 'data' is supplied.",
      call. = FALSE
    )
  }
  if(is.null(fitted_rows)){
    fitted_rows <- seq_len(n_rows)
  }

  fitted_n_rows <- nrow(random_term$model_matrix)
  check_int(
    fitted_rows,
    "fitted_rows",
    lower = 1L,
    upper = fitted_n_rows,
    check_length = n_rows,
    allow_NA = FALSE,
    call = paste0(context, ": ")
  )
  fitted_rows <- as.integer(fitted_rows)

  if(!is.null(new_row) && any(new_row)){
    stop(
      context, " with posterior-indexed row source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' for block '", random_term$block_name,
      "' cannot evaluate new observation rows. Supply a ",
      "parameter_source(..., values = ...) callback to compute the source ",
      "for arbitrary prediction rows.",
      call. = FALSE
    )
  }

  fitted_rows
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
