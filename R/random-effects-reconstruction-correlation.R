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
  if(structure %in% c("cs", "hcs", "ar1", "car", "har")){
    correlation <- .bt_random_effect_correlation_metadata(
      random_term,
      structure = structure,
      context = "Random-effect posterior reconstruction metadata"
    )
    rho <- .bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior,
      missing = "error",
      out_of_support = "error"
    )
    column_coordinates <- if(identical(structure, "car")){
      distance_matrix <- .bt_random_effect_correlation_draws_distance(
        random_term = random_term,
        correlation = correlation,
        structure = structure,
        n_columns = n_columns,
        context = "Random-effect posterior reconstruction metadata"
      )
      distance_matrix[1L, ]
    }else{
      NULL
    }
    columns <- seq_len(n_columns)
    out <- array(NA_real_, dim = c(nrow(posterior), n_columns, n_columns))
    for(draw in seq_len(nrow(posterior))){
      out[draw, , ] <- .bt_random_effect_structured_subset_cholesky(
        structure = structure,
        columns = columns,
        rho = rho[draw],
        global_n_columns = n_columns,
        column_coordinates = column_coordinates
      )
    }
    return(out)
  }

  if(identical(structure, "us")){
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

    u_names <- .bt_random_effect_lkj_primitive_names(random_term, n_columns)
    if(all(u_names %in% colnames(posterior))){
      return(.bt_lkj_cholesky_cpc_u_to_L(
        posterior[, u_names, drop = FALSE],
        K = n_columns
      ))
    }
  }

  NULL
}

.bt_random_effect_rho_draws <- function(random_term, posterior,
                                        missing = c("null", "error"),
                                        out_of_support = c("null", "error"),
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
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " does not define canonical scalar correlation metadata.",
      call. = FALSE
    )
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
  if(rho_source %in% c("sample", "fixed_sample") &&
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

  if(!is.numeric(value) || any(!is.finite(value))){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " scalar correlation coordinates must be finite.",
      call. = FALSE
    )
  }

  rho_scale <- .bt_random_effect_rho_scale_metadata(correlation, random_term, context)
  bounds <- .bt_random_effect_rho_bounds_metadata(correlation, random_term, context)
  structure <- .bt_random_effect_structure(random_term, context = context)
  interior <- .bt_random_effect_representable_rho_bounds(bounds, structure)
  if(identical(rho_scale, "fisher_z")){
    return(pmax(interior[["lower"]], pmin(interior[["upper"]], tanh(value))))
  }
  if(identical(rho_scale, "logit")){
    return(
      interior[["lower"]] +
        (interior[["upper"]] - interior[["lower"]]) * stats::plogis(value)
    )
  }

  value
}

.bt_random_effect_representable_rho_bounds <- function(bounds, structure){

  interior <- c(
    lower = if(.bt_random_effect_rho_lower_inclusive(structure)){
      bounds[["lower"]]
    }else{
      .bt_random_effect_representable_rho_neighbor(
        bounds[["lower"]],
        direction = 1
      )
    },
    upper = .bt_random_effect_representable_rho_neighbor(
      bounds[["upper"]],
      direction = -1
    )
  )
  if(!all(is.finite(interior)) ||
     (!.bt_random_effect_rho_lower_inclusive(structure) &&
      interior[["lower"]] <= bounds[["lower"]]) ||
     (.bt_random_effect_rho_lower_inclusive(structure) &&
      interior[["lower"]] != bounds[["lower"]]) ||
     interior[["upper"]] >= bounds[["upper"]] ||
     interior[["lower"]] >= interior[["upper"]]){
    stop("Scalar correlation bounds have no representable interior.", call. = FALSE)
  }

  interior
}

.bt_random_effect_representable_rho_neighbor <- function(bound, direction){

  smallest <- .Machine$double.xmin * .Machine$double.eps
  margin <- .Machine$double.eps / 2 * abs(bound)
  if(!is.finite(margin) || margin == 0){
    margin <- smallest
  }
  for(i in seq_len(64L)){
    candidate <- bound + direction * margin
    if(is.finite(candidate) && candidate != bound){
      return(candidate)
    }
    margin <- margin * 2
  }

  stop("Scalar correlation bound has no representable interior neighbor.",
       call. = FALSE)
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
  if(anyDuplicated(primitive_names)){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " must define unique 'random_term$correlation$primitive_names'.",
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
