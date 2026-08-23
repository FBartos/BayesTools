# Group-local layouts for one-hot structured random effects.
.bt_random_effect_structured_local_layout <- function(model_matrix, group_map,
                                                      structure,
                                                      parameter_stem = NULL,
                                                      column_coordinates = NULL,
                                                      n_groups = NULL,
                                                      exact_indicator = FALSE){

  check_bool(exact_indicator, "exact_indicator", allow_NA = FALSE)
  if(!isTRUE(exact_indicator)){
    return(NULL)
  }

  if(!is.matrix(model_matrix) || !is.numeric(model_matrix) ||
     nrow(model_matrix) < 1L || ncol(model_matrix) < 1L){
    stop("Structured random-effect model matrix must be a non-empty numeric matrix.",
         call. = FALSE)
  }
  if(any(!is.finite(model_matrix))){
    stop("Structured random-effect model matrix must contain only finite values.",
         call. = FALSE)
  }
  if(!is.numeric(group_map) || length(group_map) != nrow(model_matrix) ||
     any(is.na(group_map)) || any(!is.finite(group_map)) ||
     any(group_map != as.integer(group_map)) || any(group_map < 1L)){
    stop(
      "Structured random-effect group map must contain one positive integer per model-matrix row.",
      call. = FALSE
    )
  }
  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  group_map <- as.integer(group_map)
  if(is.null(n_groups)){
    n_groups <- max(group_map)
  }
  check_int(n_groups, "n_groups", lower = max(group_map), allow_NA = FALSE)
  groups <- seq_len(n_groups)
  observed_groups <- sort(unique(group_map))
  if(any(!observed_groups %in% groups)){
    stop("Structured random-effect group map exceeds 'n_groups'.", call. = FALSE)
  }

  if(any(!model_matrix %in% c(0, 1))){
    stop(
      "A design declared as an exact indicator matrix contains a non-binary value.",
      call. = FALSE
    )
  }
  nonzero <- model_matrix == 1
  if(any(rowSums(nonzero) != 1L)){
    stop(
      "Group-local structured compilation requires exactly one non-zero design value per row.",
      call. = FALSE
    )
  }
  row_column <- max.col(nonzero, ties.method = "first")
  selected   <- model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)]
  if(any(selected != 1)){
    stop(
      "Group-local structured compilation requires the non-zero design value in every row to equal one.",
      call. = FALSE
    )
  }

  n_columns <- ncol(model_matrix)
  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = n_columns,
    column_coordinates = column_coordinates
  )
  group_columns <- lapply(groups, function(group){
    columns <- unique(row_column[group_map == group])
    columns[order(coordinates[columns], columns)]
  })
  names(group_columns) <- as.character(groups)

  row_local <- integer(length(group_map))
  for(group in observed_groups){
    rows <- which(group_map == group)
    row_local[rows] <- match(row_column[rows], group_columns[[group]])
  }

  local_group  <- rep(groups, lengths(group_columns))
  local_column <- unlist(group_columns, use.names = FALSE)
  local_index  <- unlist(lapply(lengths(group_columns), seq_len), use.names = FALSE)
  node_names   <- if(is.null(parameter_stem)){
    character()
  }else{
    .bt_random_effect_structured_local_node_names(
      parameter_stem = parameter_stem,
      group = local_group,
      column = local_column
    )
  }

  out <- list(
    type = "group_local",
    structure = structure,
    global_n_columns = n_columns,
    n_groups = length(groups),
    all_groups_observed = length(observed_groups) == length(groups),
    n_local = length(local_column),
    group_columns = group_columns,
    group_coordinates = lapply(group_columns, function(columns) coordinates[columns]),
    row_column = row_column,
    row_local = row_local,
    local_group = local_group,
    local_column = local_column,
    local_index = local_index,
    column_coordinates = coordinates,
    node_names = node_names
  )
  class(out) <- c("BayesTools_random_effect_structured_local_layout", "list")

  out
}

# Normalize the scalar-correlation structures supported by local compilation.
.bt_random_effect_structured_local_normalize_structure <- function(structure){

  check_char(
    structure,
    "structure",
    allow_values = c("cs", "hcs", "ar", "ar1", "har", "car",
                     "CS", "HCS", "AR", "AR1", "HAR", "CAR"),
    allow_NA = FALSE
  )
  structure <- tolower(structure)
  if(identical(structure, "ar")){
    structure <- "ar1"
  }

  structure
}

# Resolve the global index coordinate attached to every structured column.
.bt_random_effect_structured_local_coordinates <- function(
    structure, n_columns, column_coordinates = NULL){

  if(is.null(column_coordinates)){
    if(identical(structure, "car")){
      stop(
        "Group-local CAR compilation requires one numeric coordinate per global column.",
        call. = FALSE
      )
    }
    return(seq_len(n_columns))
  }
  if(!is.numeric(column_coordinates) || length(column_coordinates) != n_columns ||
     any(is.na(column_coordinates)) || any(!is.finite(column_coordinates)) ||
     anyDuplicated(column_coordinates)){
    stop(
      "Structured random-effect column coordinates must be finite, unique, and match the global column count.",
      call. = FALSE
    )
  }
  if(structure %in% c("ar1", "har") &&
     any(column_coordinates != as.integer(column_coordinates))){
    stop(
      "AR1/HAR column coordinates must be integer ordered-index positions.",
      call. = FALSE
    )
  }

  as.numeric(column_coordinates)
}

# Generate sparse latent-node names without changing the historical convention.
.bt_random_effect_structured_local_node_names <- function(parameter_stem,
                                                          group, column){

  check_char(parameter_stem, "parameter_stem", allow_NA = FALSE)
  if(length(parameter_stem) != 1L || !nzchar(parameter_stem)){
    stop("'parameter_stem' must be one non-empty string.", call. = FALSE)
  }
  if(!is.numeric(group) || !is.numeric(column) || length(group) != length(column) ||
     any(is.na(group)) || any(is.na(column)) ||
     any(group != as.integer(group)) || any(column != as.integer(column)) ||
     any(group < 1L) || any(column < 1L)){
    stop("Local random-effect group and column indices must be positive integers.",
         call. = FALSE)
  }

  paste0(
    parameter_stem,
    "_xRE_Zx[", as.integer(group), ",", as.integer(column), "]"
  )
}

# Validate rho against the full global covariance dimension.
.bt_random_effect_structured_local_check_rho <- function(structure, rho,
                                                         global_n_columns){

  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  check_int(global_n_columns, "global_n_columns", lower = 1, allow_NA = FALSE)
  check_real(rho, "rho", allow_NA = FALSE)
  if(length(rho) != 1L || !is.finite(rho)){
    stop("'rho' must be one finite numeric value.", call. = FALSE)
  }
  .bt_random_effect_structured_check_rho_draws(
    structure        = structure,
    rho              = rho,
    global_n_columns = global_n_columns
  )
}

# Validate a vector of scalar-correlation draws against the global structure.
.bt_random_effect_structured_check_rho_draws <- function(
    structure, rho, global_n_columns){

  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  check_int(global_n_columns, "global_n_columns", lower = 1, allow_NA = FALSE)
  if(!is.numeric(rho) || length(rho) < 1L || any(!is.finite(rho))){
    stop("'rho' must contain one or more finite numeric values.", call. = FALSE)
  }
  if(global_n_columns == 1L){
    return(invisible(TRUE))
  }

  bounds <- .bt_random_effect_structured_rho_bounds(
    K = global_n_columns,
    structure = structure
  )
  lower_valid <- if(identical(structure, "car")){
    rho >= bounds[["lower"]]
  }else{
    rho > bounds[["lower"]]
  }
  valid <- lower_valid & rho < bounds[["upper"]]
  if(any(!valid)){
    interval <- paste0(
      if(identical(structure, "car")) "[" else "(",
      bounds[["lower"]], ", ", bounds[["upper"]], ")"
    )
    stop(
      "'rho' is outside the global ", toupper(structure),
      " support for K = ", global_n_columns, ": expected ", interval, ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Format the caller identity attached to structured Markov numerical failures.
.bt_random_effect_markov_error_label <- function(context, operation){

  if(is.null(context)){
    return(paste0("Structured Markov ", operation))
  }
  if(!is.character(context) || length(context) != 1L ||
     is.na(context) || !nzchar(context)){
    stop("'context' must be NULL or one non-empty string.", call. = FALSE)
  }

  paste0(context, ": structured Markov ", operation)
}

.bt_random_effect_batch_context <- function(context, row){

  if(is.function(context)){
    context <- context(row)
  }
  if(is.null(context)){
    return(NULL)
  }
  if(!is.character(context) || length(context) != 1L ||
     is.na(context) || !nzchar(context)){
    stop(
      "Batch 'context' must be NULL, one non-empty string, or a function returning one non-empty string.",
      call. = FALSE
    )
  }

  context
}

# Return stable Markov transitions for a vector of rho draws.
.bt_random_effect_markov_transition_draws <- function(
    rho, left_coordinate, right_coordinate, context = NULL){

  if(!is.numeric(rho) || length(rho) < 1L || any(!is.finite(rho))){
    stop("'rho' must contain one or more finite numeric values.", call. = FALSE)
  }
  if(!is.numeric(left_coordinate) || length(left_coordinate) != 1L ||
     !is.finite(left_coordinate) ||
     !is.numeric(right_coordinate) || length(right_coordinate) != 1L ||
     !is.finite(right_coordinate)){
    stop(
      "Structured Markov transition coordinates must be finite numeric values.",
      call. = FALSE
    )
  }

  gap <- right_coordinate - left_coordinate
  if(!is.finite(gap) || gap <= 0){
    transition_values <- format(
      c(left_coordinate, right_coordinate, gap, rho[[1L]]),
      digits = 17L,
      scientific = TRUE,
      trim = TRUE
    )
    error_label <- .bt_random_effect_markov_error_label(
      .bt_random_effect_batch_context(context, 1L),
      operation = "transition"
    )
    stop(
      error_label, " from coordinate ",
      transition_values[1L], " to ", transition_values[2L],
      " has a non-positive or non-finite gap (rho = ",
      transition_values[4L], ", gap = ",
      transition_values[3L], ").",
      call. = FALSE
    )
  }

  negative <- rho < 0
  if(any(negative) && gap != floor(gap)){
    row <- which(negative)[[1L]]
    transition_values <- format(
      c(left_coordinate, right_coordinate, gap, rho[[row]]),
      digits = 17L,
      scientific = TRUE,
      trim = TRUE
    )
    error_label <- .bt_random_effect_markov_error_label(
      .bt_random_effect_batch_context(context, row),
      operation = "transition"
    )
    stop(
      error_label, " with a negative correlation requires an integer ",
      "coordinate gap; transition from coordinate ",
      transition_values[1L], " to ", transition_values[2L],
      " has rho = ", transition_values[4L],
      " and gap = ", transition_values[3L], ".",
      call. = FALSE
    )
  }

  log_phi <- rep(-Inf, length(rho))
  nonzero <- rho != 0
  log_phi[nonzero] <- gap * base::log(abs(rho[nonzero]))
  phi <- base::exp(log_phi)
  if(gap %% 2 != 0){
    phi[negative] <- -phi[negative]
  }
  innovation_variance <- -base::expm1(2 * log_phi)
  invalid <- !is.finite(innovation_variance) | innovation_variance <= 0
  if(any(invalid)){
    row <- which(invalid)[[1L]]
    transition_values <- format(
      c(left_coordinate, right_coordinate, gap, rho[[row]]),
      digits = 17L,
      scientific = TRUE,
      trim = TRUE
    )
    error_label <- .bt_random_effect_markov_error_label(
      .bt_random_effect_batch_context(context, row),
      operation = "transition"
    )
    stop(
      error_label, " from coordinate ",
      transition_values[1L], " to ", transition_values[2L],
      " has a non-positive or non-finite innovation variance (rho = ",
      transition_values[4L], ", gap = ",
      transition_values[3L], "). The requested coordinate/time resolution ",
      "is not representable for this rho.",
      call. = FALSE
    )
  }

  list(
    log_phi = log_phi,
    phi = phi,
    innovation_variance = innovation_variance
  )
}

# Return a stable Markov transition over one increasing coordinate interval.
.bt_random_effect_markov_transition <- function(
    rho, left_coordinate, right_coordinate, context = NULL){

  if(!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)){
    stop("'rho' must be one finite numeric value.", call. = FALSE)
  }

  transition <- .bt_random_effect_markov_transition_draws(
    rho = rho,
    left_coordinate = left_coordinate,
    right_coordinate = right_coordinate,
    context = context
  )
  lapply(transition, `[[`, 1L)
}

# Exact Cholesky factor for a principal structured-correlation block.
.bt_random_effect_structured_subset_cholesky <- function(
    structure, columns, rho, global_n_columns,
    column_coordinates = NULL, context = NULL){

  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  check_int(global_n_columns, "global_n_columns", lower = 1, allow_NA = FALSE)
  if(!is.numeric(columns) || length(columns) < 1L || any(is.na(columns)) ||
     any(columns != as.integer(columns)) || any(columns < 1L) ||
     any(columns > global_n_columns) || anyDuplicated(columns)){
    stop("'columns' must contain unique global column indices.", call. = FALSE)
  }
  columns <- as.integer(columns)
  .bt_random_effect_structured_local_check_rho(
    structure = structure,
    rho = rho,
    global_n_columns = global_n_columns
  )
  if(length(columns) == 1L){
    if(identical(structure, "car")){
      .bt_random_effect_structured_local_coordinates(
        structure = structure,
        n_columns = global_n_columns,
        column_coordinates = column_coordinates
      )
    }
    return(matrix(1, nrow = 1L, ncol = 1L))
  }

  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = global_n_columns,
    column_coordinates = column_coordinates
  )[columns]
  if(any(diff(coordinates) <= 0)){
    stop(
      "Structured subset columns must be ordered by increasing index coordinate.",
      call. = FALSE
    )
  }

  if(structure %in% c("cs", "hcs")){
    return(.bt_random_effect_cs_subset_cholesky(length(columns), rho))
  }

  .bt_random_effect_markov_subset_cholesky(
    coordinates = coordinates,
    rho = rho,
    context = context
  )
}

# Exact principal structured-correlation block without global materialization.
.bt_random_effect_structured_subset_correlation <- function(
    structure, columns, rho, global_n_columns,
    column_coordinates = NULL, context = NULL){

  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  check_int(global_n_columns, "global_n_columns", lower = 1, allow_NA = FALSE)
  if(!is.numeric(columns) || length(columns) < 1L || any(is.na(columns)) ||
     any(columns != as.integer(columns)) || any(columns < 1L) ||
     any(columns > global_n_columns) || anyDuplicated(columns)){
    stop("'columns' must contain unique global column indices.", call. = FALSE)
  }
  columns <- as.integer(columns)
  .bt_random_effect_structured_local_check_rho(
    structure = structure,
    rho = rho,
    global_n_columns = global_n_columns
  )
  if(structure %in% c("cs", "hcs")){
    out <- matrix(rho, nrow = length(columns), ncol = length(columns))
    diag(out) <- 1
    return(out)
  }

  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = global_n_columns,
    column_coordinates = column_coordinates
  )[columns]
  coordinate_order <- order(coordinates)
  sorted_correlation <- tcrossprod(
    .bt_random_effect_markov_subset_cholesky(
      coordinates = coordinates[coordinate_order],
      rho = rho,
      context = context
    )
  )
  original_order <- order(coordinate_order)
  correlation <- sorted_correlation[
    original_order,
    original_order,
    drop = FALSE
  ]
  diag(correlation) <- 1

  correlation
}

# Apply a principal structured Cholesky recurrence to multiple draws.
.bt_random_effect_structured_subset_transform_draws <- function(
    structure, columns, latent, rho, global_n_columns,
    column_coordinates = NULL, context = NULL){

  structure <- .bt_random_effect_structured_local_normalize_structure(structure)
  check_int(global_n_columns, "global_n_columns", lower = 1, allow_NA = FALSE)
  if(!is.numeric(columns) || length(columns) < 1L || any(is.na(columns)) ||
     any(columns != as.integer(columns)) || any(columns < 1L) ||
     any(columns > global_n_columns) || anyDuplicated(columns)){
    stop("'columns' must contain unique global column indices.", call. = FALSE)
  }
  columns <- as.integer(columns)
  if(!is.matrix(latent) || !is.numeric(latent) || nrow(latent) < 1L ||
     ncol(latent) != length(columns) || any(!is.finite(latent))){
    stop(
      "'latent' must be a finite numeric matrix with one column per structured column.",
      call. = FALSE
    )
  }
  if(!is.numeric(rho) || length(rho) != nrow(latent) || any(!is.finite(rho))){
    stop("'rho' must contain one finite value per latent row.", call. = FALSE)
  }
  .bt_random_effect_structured_check_rho_draws(
    structure = structure,
    rho = rho,
    global_n_columns = global_n_columns
  )
  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = global_n_columns,
    column_coordinates = column_coordinates
  )[columns]
  if(any(diff(coordinates) <= 0)){
    stop(
      "Structured subset columns must be ordered by increasing index coordinate.",
      call. = FALSE
    )
  }

  out <- matrix(NA_real_, nrow = nrow(latent), ncol = length(columns))
  out[, 1L] <- latent[, 1L]
  if(length(columns) == 1L){
    return(out)
  }

  if(structure %in% c("cs", "hcs")){
    prefix <- rho * latent[, 1L]
    for(index in 2:length(columns)){
      diagonal <- sqrt(
        (1 - rho) * (1 + (index - 1L) * rho) /
          (1 + (index - 2L) * rho)
      )
      update <- rho * sqrt(
        (1 - rho) /
          ((1 + (index - 2L) * rho) * (1 + (index - 1L) * rho))
      )
      out[, index] <- prefix + diagonal * latent[, index]
      prefix <- prefix + update * latent[, index]
    }
    return(out)
  }

  for(index in 2:length(columns)){
    transition <- .bt_random_effect_markov_transition_draws(
      rho = rho,
      left_coordinate = coordinates[index - 1L],
      right_coordinate = coordinates[index],
      context = context
    )
    out[, index] <- transition$phi * out[, index - 1L] +
      sqrt(transition$innovation_variance) * latent[, index]
  }

  out
}

# Scalar wrapper retained for conditional simulation and independent callers.
.bt_random_effect_structured_subset_transform <- function(
    structure, columns, latent, rho, global_n_columns,
    column_coordinates = NULL, context = NULL){

  if(!is.numeric(latent) || length(latent) != length(columns) ||
     any(!is.finite(latent))){
    stop("'latent' must contain one finite value per structured column.",
         call. = FALSE)
  }
  .bt_random_effect_structured_local_check_rho(
    structure = structure,
    rho = rho,
    global_n_columns = global_n_columns
  )
  out <- .bt_random_effect_structured_subset_transform_draws(
    structure = structure,
    columns = columns,
    latent = matrix(latent, nrow = 1L),
    rho = rho,
    global_n_columns = global_n_columns,
    column_coordinates = column_coordinates,
    context = context
  )

  as.numeric(out[1L, ])
}

# Closed-form equicorrelation Cholesky factor, including admissible negative rho.
.bt_random_effect_cs_subset_cholesky <- function(n_columns, rho){

  L <- matrix(0, nrow = n_columns, ncol = n_columns)
  L[1L, 1L] <- 1
  if(n_columns == 1L){
    return(L)
  }

  for(row in 2:n_columns){
    L[row, 1L] <- rho
    if(row > 2L){
      for(column in 2:(row - 1L)){
        L[row, column] <- rho * sqrt(
          (1 - rho) /
            ((1 + (column - 2L) * rho) * (1 + (column - 1L) * rho))
        )
      }
    }
    L[row, row] <- sqrt(
      (1 - rho) * (1 + (row - 1L) * rho) /
        (1 + (row - 2L) * rho)
    )
  }

  L
}

# Markov Cholesky recurrence for irregular AR1/CAR principal blocks.
.bt_random_effect_markov_subset_cholesky <- function(
    coordinates, rho, context = NULL){

  n_columns <- length(coordinates)
  L <- matrix(0, nrow = n_columns, ncol = n_columns)
  L[1L, 1L] <- 1
  if(n_columns == 1L){
    return(L)
  }

  for(row in 2:n_columns){
    transition <- .bt_random_effect_markov_transition(
      rho = rho,
      left_coordinate = coordinates[row - 1L],
      right_coordinate = coordinates[row],
      context = context
    )
    L[row, ] <- transition$phi * L[row - 1L, ]
    L[row, row] <- sqrt(transition$innovation_variance)
  }

  L
}

.bt_random_effect_native_structured_cholesky <- function(
    structure, rho, coordinates){

  .Call(
    "BayesTools_structured_cholesky",
    as.numeric(rho),
    as.numeric(coordinates),
    structure %in% c("cs", "hcs"),
    PACKAGE = "BayesTools"
  )
}

