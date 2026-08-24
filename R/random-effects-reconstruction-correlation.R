.bt_random_effect_cholesky_draws <- function(random_term, n_columns,
                                            posterior){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect posterior reconstruction metadata"
  )
  cholesky_evaluator <- .bt_random_effect_compile_cholesky_evaluator(
    random_term,
    n_columns = n_columns,
    structure = structure,
    posterior_names = colnames(posterior)
  )
  if(is.null(cholesky_evaluator)){
    return(NULL)
  }
  cholesky_evaluator(posterior)
}

.bt_random_effect_compile_cholesky_evaluator <- function(
    random_term, n_columns, structure = NULL, posterior_names = NULL){

  if(is.null(structure)){
    structure <- .bt_random_effect_structure(
      random_term,
      context = "Random-effect posterior reconstruction metadata"
    )
  }
  if(structure %in% c("diag", "id") || n_columns <= 1L){
    force(n_columns)
    return(function(posterior){
      out <- array(0, dim = c(nrow(posterior), n_columns, n_columns))
      for(column in seq_len(n_columns)){
        out[, column, column] <- 1
      }
      out
    })
  }

  if(identical(structure, "us")){
    primitive_names <- .bt_random_effect_lkj_primitive_names(
      random_term,
      n_columns
    )
    cholesky_names <- as.vector(.bt_random_effect_cholesky_names(
      random_term = random_term,
      n_columns = n_columns
    ))
    fixed_posterior_names <- !is.null(posterior_names)
    cached_posterior_names <- posterior_names
    primitive_indices <- if(fixed_posterior_names){
      match(primitive_names, posterior_names)
    }else{
      rep(NA_integer_, length(primitive_names))
    }
    cholesky_indices <- if(fixed_posterior_names){
      match(cholesky_names, posterior_names)
    }else{
      rep(NA_integer_, length(cholesky_names))
    }
    force(primitive_names)
    force(cholesky_names)
    force(n_columns)

    return(function(posterior){
      if(!fixed_posterior_names){
        current_names <- colnames(posterior)
        if(!identical(current_names, cached_posterior_names)){
          primitive_indices <<- match(primitive_names, current_names)
          cholesky_indices  <<- match(cholesky_names, current_names)
          cached_posterior_names <<- current_names
        }
      }
      if(length(primitive_names) > 0L && !anyNA(primitive_indices)){
        return(.bt_lkj_cholesky_cpc_u_to_L(
          posterior[, primitive_indices, drop = FALSE],
          K = n_columns
        ))
      }
      if(length(cholesky_names) > 0L && !anyNA(cholesky_indices)){
        values <- posterior[, cholesky_indices, drop = FALSE]
        return(array(
          unname(values),
          dim = c(nrow(posterior), n_columns, n_columns)
        ))
      }
      NULL
    })
  }

  if(!structure %in% c("cs", "hcs", "ar1", "car", "har")){
    return(NULL)
  }

  context <- "Random-effect posterior reconstruction metadata"
  rho_plan <- .bt_random_effect_compile_rho_draw_plan(
    random_term = random_term,
    context = context
  )
  correlation <- rho_plan$correlation
  coordinates <- if(identical(structure, "car")){
    .bt_random_effect_car_time_values(
      random_term = random_term,
      correlation = correlation,
      n_columns = n_columns,
      context = context
    )
  }else{
    seq_len(n_columns)
  }
  coordinates <- .bt_random_effect_structured_local_coordinates(
    structure = structure,
    n_columns = n_columns,
    column_coordinates = coordinates
  )
  structure_bounds <- .bt_random_effect_structured_rho_bounds(
    K = n_columns,
    structure = structure
  )
  rho_evaluator <- .bt_random_effect_compile_rho_draw_evaluator(
    random_term = random_term,
    missing = "error",
    out_of_support = "error",
    context = context,
    plan = rho_plan,
    posterior_names = posterior_names
  )
  force(random_term)
  force(n_columns)
  force(structure)
  force(rho_plan)
  force(coordinates)
  force(structure_bounds)
  force(rho_evaluator)

  function(posterior){
    rho <- rho_evaluator(posterior)
    invalid <- .bt_random_effect_rho_outside_support(
      rho,
      bounds = structure_bounds,
      structure = structure
    )
    if(any(invalid)){
      .bt_random_effect_structured_local_check_rho(
        structure = structure,
        rho = rho[which(invalid)[1L]],
        global_n_columns = n_columns
      )
    }

    tryCatch(
      .bt_random_effect_native_structured_cholesky(
        structure = structure,
        rho = rho,
        coordinates = coordinates
      ),
      error = function(error){
        for(draw in seq_len(nrow(posterior))){
          transition_context <- paste0(
            "Random-effect Cholesky reconstruction",
            .bt_random_effect_metadata_block_detail(random_term),
            ", posterior draw ", draw
          )
          .bt_random_effect_structured_subset_cholesky(
            structure = structure,
            columns = seq_len(n_columns),
            rho = rho[draw],
            global_n_columns = n_columns,
            column_coordinates = coordinates,
            context = transition_context
          )
        }
        stop(conditionMessage(error), call. = FALSE)
      }
    )
  }
}

.bt_random_effect_compile_rho_draw_evaluator <- function(
    random_term,
    missing = c("null", "error"),
    out_of_support = c("null", "error"),
    context = "Random-effect posterior reconstruction metadata",
    plan = NULL,
    posterior_names = NULL){

  missing        <- match.arg(missing)
  out_of_support <- match.arg(out_of_support)
  if(is.null(plan)){
    plan <- .bt_random_effect_compile_rho_draw_plan(
      random_term = random_term,
      context = context
    )
  }
  structure   <- plan$structure
  correlation <- plan$correlation
  rho_scale   <- plan$rho_scale
  bounds      <- plan$bounds
  sample_fixed <- plan$sample_fixed
  sample_bounds <- plan$sample_bounds
  sample_name <- correlation$sample_name
  rho_name    <- correlation$rho_name
  fixed_posterior_names <- !is.null(posterior_names)
  cached_posterior_names <- posterior_names
  sample_index <- if(fixed_posterior_names){
    match(sample_name, posterior_names)
  }else{
    NA_integer_
  }
  rho_index <- if(fixed_posterior_names){
    match(rho_name, posterior_names)
  }else{
    NA_integer_
  }
  force(random_term)
  force(plan)
  force(missing)
  force(out_of_support)
  force(context)

  function(posterior){
    if(!fixed_posterior_names){
      current_names <- colnames(posterior)
      if(!identical(current_names, cached_posterior_names)){
        sample_index <<- match(sample_name, current_names)
        rho_index    <<- match(rho_name, current_names)
        cached_posterior_names <<- current_names
      }
    }

    sample_value <- NULL
    if(!identical(rho_scale, "rho") && !is.na(sample_index)){
      sample_value <- posterior[, sample_index]
      rho_source <- "sample"
      rho <- .bt_random_effect_transform_rho(
        sample_value,
        correlation = correlation,
        random_term = random_term,
        context = context,
        plan = plan
      )
    }else if(!is.na(rho_index)){
      rho_source <- "rho"
      rho <- posterior[, rho_index]
    }else if(!is.na(sample_index)){
      sample_value <- posterior[, sample_index]
      rho_source <- "sample"
      rho <- .bt_random_effect_transform_rho(
        sample_value,
        correlation = correlation,
        random_term = random_term,
        context = context,
        plan = plan
      )
    }else if(!is.null(sample_fixed)){
      sample_value <- rep(sample_fixed, nrow(posterior))
      rho_source <- "fixed_sample"
      rho <- .bt_random_effect_transform_rho(
        sample_value,
        correlation = correlation,
        random_term = random_term,
        context = context,
        plan = plan
      )
    }else{
      if(identical(missing, "error")){
        .bt_random_effect_missing_rho_draws_stop(
          random_term = random_term,
          correlation = correlation,
          context = context
        )
      }
      return(NULL)
    }

    if(rho_source %in% c("sample", "fixed_sample") &&
       !identical(rho_scale, "rho")){
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
}

.bt_random_effect_rho_draws <- function(random_term, posterior,
                                        missing = c("null", "error"),
                                        out_of_support = c("null", "error"),
                                        context = "Random-effect posterior reconstruction metadata",
                                        plan = NULL){

  missing <- match.arg(missing)
  out_of_support <- match.arg(out_of_support)
  if(is.null(plan)){
    plan <- .bt_random_effect_compile_rho_draw_plan(
      random_term = random_term,
      context = context
    )
  }
  structure   <- plan$structure
  correlation <- plan$correlation
  rho_scale   <- plan$rho_scale
  if(!identical(rho_scale, "rho") &&
     correlation$sample_name %in% colnames(posterior)){
    sample_value <- posterior[, correlation$sample_name]
    rho_source <- "sample"
    rho <- .bt_random_effect_transform_rho(
      sample_value,
      correlation = correlation,
      random_term = random_term,
      context = context,
      plan = plan
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
      context = context,
      plan = plan
    )
  }else{
    sample_fixed <- plan$sample_fixed
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
      context = context,
      plan = plan
    )
  }

  bounds <- plan$bounds
  if(rho_source %in% c("sample", "fixed_sample") &&
     !identical(rho_scale, "rho")){
    sample_bounds <- plan$sample_bounds
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

.bt_random_effect_compile_rho_draw_plan <- function(
    random_term,
    context = "Random-effect posterior reconstruction metadata"){

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
  bounds <- .bt_random_effect_rho_bounds_metadata(
    correlation,
    random_term,
    context = context
  )
  sample_fixed <- .bt_random_effect_rho_fixed_sample_metadata(
    correlation,
    random_term,
    context = context
  )
  sample_bounds <- if(identical(rho_scale, "rho")){
    NULL
  }else{
    .bt_random_effect_rho_sample_bounds(
      correlation = correlation,
      random_term = random_term,
      context = context
    )
  }
  list(
    structure = structure,
    correlation = correlation,
    rho_scale = rho_scale,
    bounds = bounds,
    sample_fixed = sample_fixed,
    sample_bounds = sample_bounds
  )
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
                                            context = "Random-effect posterior reconstruction metadata",
                                            plan = NULL){

  if(!is.numeric(value) || any(!is.finite(value))){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " scalar correlation coordinates must be finite.",
      call. = FALSE
    )
  }

  if(is.null(plan)){
    rho_scale <- .bt_random_effect_rho_scale_metadata(
      correlation,
      random_term,
      context
    )
    bounds <- .bt_random_effect_rho_bounds_metadata(
      correlation,
      random_term,
      context
    )
  }else{
    rho_scale <- plan$rho_scale
    bounds    <- plan$bounds
  }
  if(identical(rho_scale, "fisher_z")){
    return(tanh(value))
  }
  if(identical(rho_scale, "logit")){
    return(
      bounds[["lower"]] +
        (bounds[["upper"]] - bounds[["lower"]]) * stats::plogis(value)
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
