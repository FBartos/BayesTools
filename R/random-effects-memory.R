.bt_random_effect_memory_limit_bytes <- function(){

  limit <- getOption(
    "BayesTools.random_effects_memory_limit_bytes",
    16 * 1024^3
  )
  if(!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
     limit <= 0){
    stop(
      "Option 'BayesTools.random_effects_memory_limit_bytes' must be a ",
      "positive numeric scalar or Inf.",
      call. = FALSE
    )
  }

  as.numeric(limit)
}

.bt_random_effect_memory_product <- function(...){

  values <- as.numeric(c(...))
  if(anyNA(values) || any(values < 0) || any(!is.finite(values))){
    return(Inf)
  }
  if(any(values == 0)){
    return(0)
  }
  log_product <- sum(log(values))
  if(log_product > log(.Machine$double.xmax)){
    return(Inf)
  }

  prod(values)
}

.bt_random_effect_design_memory_estimate <- function(
    n_rows,
    n_columns,
    n_groups,
    monitor_policy,
    structure,
    compile_mode){

  design_payload <- .bt_random_effect_memory_product(
    8,
    n_rows,
    n_columns
  )
  design_working <- 4 * design_payload
  sampled <- identical(compile_mode, "sampled")
  group_state_copies <- if(sampled){
    2 +
      as.integer(isTRUE(monitor_policy$latent)) +
      as.integer(isTRUE(monitor_policy$coefficients))
  }else{
    0
  }
  group_state <- .bt_random_effect_memory_product(
    8,
    n_groups,
    n_columns,
    group_state_copies
  )
  covariance_copies <- if(structure %in% c(
    "us", "cs", "hcs", "ar1", "car", "har"
  )){
    3
  }else{
    0
  }
  covariance_state <- .bt_random_effect_memory_product(
    8,
    n_columns,
    n_columns,
    covariance_copies
  )
  row_maps <- .bt_random_effect_memory_product(4, n_rows, 4)
  peak <- sum(
    design_working,
    group_state,
    covariance_state,
    row_maps
  )
  if(!is.finite(peak)){
    peak <- Inf
  }

  list(
    operation = "random-effect design construction",
    dimensions = paste0(
      format(n_rows, scientific = FALSE, trim = TRUE),
      " rows x ",
      format(n_columns, scientific = FALSE, trim = TRUE),
      " columns; ",
      format(n_groups, scientific = FALSE, trim = TRUE),
      " groups"
    ),
    payload_bytes = design_payload,
    peak_bytes = peak,
    components = c(
      design_working = design_working,
      group_state = group_state,
      covariance_state = covariance_state,
      row_maps = row_maps
    )
  )
}

.bt_random_effect_output_memory_estimate <- function(
    operation,
    n_rows,
    n_draws,
    covariance = FALSE,
    diagonal_only = FALSE){

  entries <- if(isTRUE(covariance) && !isTRUE(diagonal_only)){
    .bt_random_effect_memory_product(n_draws, n_rows, n_rows)
  }else{
    .bt_random_effect_memory_product(n_draws, n_rows)
  }
  payload <- .bt_random_effect_memory_product(8, entries)
  peak <- 3 * payload
  if(!is.finite(peak)){
    peak <- Inf
  }

  list(
    operation = operation,
    dimensions = if(isTRUE(covariance) && !isTRUE(diagonal_only)){
      paste0(
        format(n_draws, scientific = FALSE, trim = TRUE),
        " draws x ",
        format(n_rows, scientific = FALSE, trim = TRUE),
        " rows x ",
        format(n_rows, scientific = FALSE, trim = TRUE),
        " rows"
      )
    }else{
      paste0(
        format(n_rows, scientific = FALSE, trim = TRUE),
        " rows x ",
        format(n_draws, scientific = FALSE, trim = TRUE),
        " draws"
      )
    },
    payload_bytes = payload,
    peak_bytes = peak,
    components = c(output_working = peak)
  )
}

.bt_random_effect_check_memory <- function(estimate, block_name = NULL,
                                           alternative){

  limit <- .bt_random_effect_memory_limit_bytes()
  if(is.infinite(limit) || estimate$peak_bytes <= limit){
    return(invisible(estimate))
  }

  block_detail <- if(is.null(block_name)){
    ""
  }else{
    paste0(" for block '", block_name, "'")
  }
  stop(
    "Estimated peak memory for ", estimate$operation, block_detail,
    " exceeds option 'BayesTools.random_effects_memory_limit_bytes'. ",
    "Dimensions: ", estimate$dimensions,
    "; primary payload: ",
    .bt_random_effect_format_bytes(estimate$payload_bytes),
    "; conservative peak estimate: ",
    .bt_random_effect_format_bytes(estimate$peak_bytes),
    "; configured ceiling: ",
    .bt_random_effect_format_bytes(limit),
    ". This guard cannot promise that smaller allocations will succeed. ",
    alternative,
    call. = FALSE
  )
}

.bt_random_effect_format_bytes <- function(bytes){

  if(is.infinite(bytes)){
    return("Inf bytes")
  }
  units <- c("bytes", "KiB", "MiB", "GiB", "TiB")
  unit <- min(
    length(units),
    max(1L, floor(log(max(bytes, 1), base = 1024)) + 1L)
  )
  scaled <- bytes / 1024^(unit - 1L)

  paste0(
    format(bytes, scientific = FALSE, trim = TRUE, big.mark = ","),
    " bytes (",
    format(round(scaled, 2L), nsmall = if(scaled < 10) 2L else 1L,
           trim = TRUE),
    " ",
    units[unit],
    ")"
  )
}

.bt_random_effect_design_preflight <- function(
    formula,
    data,
    preserve_no_intercept_contrasts,
    structure,
    block_name){

  if(identical(structure, "car")){
    formula_terms <- stats::terms(formula)
    predictors <- as.character(attr(formula_terms, "variables"))[-1L]
    if(length(predictors) != 1L || !predictors %in% names(data)){
      return(list(n_rows = nrow(data), n_columns = NULL))
    }
    time_values <- .bt_random_effect_car_observed_time(
      data[[predictors]],
      predictors
    )
    return(list(
      n_rows = nrow(data),
      n_columns = length(unique(time_values))
    ))
  }

  prototype <- data[1L, , drop = FALSE]
  prototype_design <- .bt_random_effect_design_matrix(
    formula = formula,
    data = prototype,
    preserve_no_intercept_contrasts = preserve_no_intercept_contrasts,
    structure = structure,
    block_name = block_name
  )

  list(
    n_rows = nrow(data),
    n_columns = ncol(prototype_design$model_matrix)
  )
}
