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

# The block design evaluated at every level cell of a random factor term (the
# factor design of the term's SD factor prior): the dense cells x block
# columns design, with the working copies of design construction, and the
# data frame of the cells.
.bt_random_effect_factor_design_memory_estimate <- function(n_cells,
                                                            n_columns,
                                                            n_predictors){

  design_payload <- .bt_random_effect_memory_product(8, n_cells, n_columns)
  design_working <- 4 * design_payload
  cell_data <- .bt_random_effect_memory_product(8, n_cells, max(n_predictors, 1))
  peak <- sum(design_working, cell_data)
  if(!is.finite(peak)){
    peak <- Inf
  }

  list(
    operation = "random-effect factor-level design construction",
    dimensions = paste0(
      format(n_cells, scientific = FALSE, trim = TRUE),
      " level cells x ",
      format(n_columns, scientific = FALSE, trim = TRUE),
      " columns"
    ),
    payload_bytes = design_payload,
    peak_bytes = peak,
    components = c(design_working = design_working, cell_data = cell_data)
  )
}

# The marginal covariance output and, when the block is described, the working
# arrays of its coefficient columns (`n_columns` is the block's number of
# random-effect columns):
# - sd_draws_working: the SD draw matrix (draws x columns) and its validation
#   copy, read for every structure;
# - coefficient_factor_working / coefficient_contraction_working:
#   `factor_structure` blocks (a correlation structure other than ID, DIAG and
#   CAR) hold the Cholesky and correlation draws (2 x draws x columns^2) and a
#   rows x columns design contraction (4 copies), as `car_columns` does for CAR.
#   The Cholesky reader reshapes its column slice in place and the validation
#   reads one column at a time, so no further full copy of either exists;
# - known_group_working: a known group kernel forms, for the full output, the
#   rows x rows kernel expansion (K[group_map, group_map]), its lower-triangle
#   mask and three rows x rows covariance copies per draw (44 bytes per
#   entry); its diagonal-only output keeps two rows x columns contractions.
.bt_random_effect_output_memory_estimate <- function(
    operation,
    n_rows,
    n_draws,
    covariance = FALSE,
    diagonal_only = FALSE,
    car_columns = NULL,
    n_columns = NULL,
    factor_structure = FALSE,
    known_group = FALSE,
    column_contributions = FALSE,
    group_blocks = NULL){

  entries <- if(isTRUE(covariance) && !isTRUE(diagonal_only)){
    .bt_random_effect_memory_product(n_draws, n_rows, n_rows)
  }else{
    .bt_random_effect_memory_product(n_draws, n_rows)
  }
  payload <- .bt_random_effect_memory_product(8, entries)
  components <- c(output_working = 3 * payload)
  if(!is.null(group_blocks) && isTRUE(covariance) && !isTRUE(diagonal_only)){
    # one rows x rows group matrix per block is formed before the draws loop
    components <- c(components,
      group_covariance_working = .bt_random_effect_memory_product(8, group_blocks, n_rows, n_rows))
  }
  if(!is.null(car_columns)){
    components <- c(components,
      car_factor_working = .bt_random_effect_memory_product(2, 8, n_draws, car_columns, car_columns),
      car_contraction_working = .bt_random_effect_memory_product(4, 8, n_rows, car_columns))
  }
  if(!is.null(n_columns)){
    components <- c(components,
      sd_draws_working = .bt_random_effect_memory_product(2, 8, n_draws, n_columns))
    if(isTRUE(factor_structure)){
      components <- c(components,
        coefficient_factor_working = .bt_random_effect_memory_product(2, 8, n_draws, n_columns, n_columns),
        coefficient_contraction_working = .bt_random_effect_memory_product(4, 8, n_rows, n_columns))
    }
    if(isTRUE(column_contributions)){
      # The rows x draws x columns array of the per-column contributions.
      components <- c(components,
        column_contribution_working = .bt_random_effect_memory_product(8, n_rows, n_draws, n_columns))
    }
    if(isTRUE(known_group)){
      components <- c(components,
        known_group_working = if(isTRUE(covariance) && !isTRUE(diagonal_only)){
          .bt_random_effect_memory_product(5 * 8 + 4, n_rows, n_rows)
        }else{
          .bt_random_effect_memory_product(2, 8, n_rows, n_columns)
        })
    }
  }
  peak <- sum(components)
  if(!is.finite(peak)){
    peak <- Inf
  }

  list(
    operation = operation,
    dimensions = paste0(if(isTRUE(covariance) && !isTRUE(diagonal_only)){
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
    }, if(!is.null(car_columns)){
      paste0(
        "; ", format(car_columns, scientific = FALSE, trim = TRUE),
        " coefficient/time columns"
      )
    }else if(!is.null(n_columns)){
      paste0(
        "; ", format(n_columns, scientific = FALSE, trim = TRUE),
        " coefficient columns"
      )
    }),
    payload_bytes = payload,
    peak_bytes = peak,
    components = components
  )
}

# The rows x rows covariances of the multi-column marginal update plan
# (.bt_random_effect_marginal_update_invariant_covariance()). One dense matrix
# per block is kept until the blocks are summed (8 bytes per entry). Assembling
# a block forms the group indicator (4 bytes), its double copy, the design
# cross-product and their product (8 bytes each) before the temporaries are
# collected, and the summed basis comes with a zero base matrix of the same
# size: 36 bytes per entry of working copies. (The expanded known-group kernel
# of a block with a known group covariance replaces the indicator and its
# copy, which is smaller.) The estimate is conservative against the garbage
# collector's high-water mark measured for n = 600 and 1200 rows.
.bt_random_effect_update_plan_memory_estimate <- function(n_rows, n_blocks){

  entries <- .bt_random_effect_memory_product(n_rows, n_rows)
  payload <- .bt_random_effect_memory_product(8, entries)
  components <- c(
    block_covariances = .bt_random_effect_memory_product(8, n_blocks, entries),
    assembly_working = .bt_random_effect_memory_product(36, entries)
  )
  peak <- sum(components)
  if(!is.finite(peak)){
    peak <- Inf
  }

  list(
    operation = "random-effect marginal update plan",
    dimensions = paste0(
      format(n_rows, scientific = FALSE, trim = TRUE), " rows x ",
      format(n_rows, scientific = FALSE, trim = TRUE), " rows; ",
      format(n_blocks, scientific = FALSE, trim = TRUE),
      if(n_blocks == 1L) " block" else " blocks"
    ),
    payload_bytes = payload,
    peak_bytes = peak,
    components = components
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
  stop(errorCondition(
    message = paste0(
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
      alternative
    ),
    class = c("BayesTools_random_effects_memory_limit",
              "BayesTools_random_effects_error"),
    call = NULL,
    operation = estimate$operation,
    peak_bytes = estimate$peak_bytes,
    limit_bytes = limit
  ))
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
