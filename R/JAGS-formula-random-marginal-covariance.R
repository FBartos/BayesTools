#' Compile formula random effects as marginal JAGS covariance nodes
#'
#' @description
#' Compiles marginalized formula random-effect metadata into deterministic JAGS
#' nodes containing observation-level covariance contributions. The compiler is
#' intentionally likelihood agnostic: callers provide dependency-preserving row
#' blocks and add the returned covariance to their sampling covariance.
#'
#' @param formula_design a `BayesTools_formula_design` object created by
#'   [JAGS_formula()].
#' @param row_blocks list of integer vectors partitioning the fitted rows.
#' @param prefix valid JAGS node prefix used for generated nodes and data.
#'
#' @return A list with `syntax`, `data`, `row_blocks`, and one lower-triangle
#'   covariance node name per row block in `lower_names`.
#'
#' @export
JAGS_formula_random_marginal_covariance <- function(
    formula_design, row_blocks, prefix = "random_marginal_covariance"){

  if(!inherits(formula_design, "BayesTools_formula_design")){
    stop("'formula_design' must be a BayesTools formula design.", call. = FALSE)
  }
  check_char(prefix, "prefix", check_length = 1, allow_NA = FALSE)
  if(!grepl("^[A-Za-z][A-Za-z0-9_]*$", prefix)){
    stop("'prefix' must be a valid JAGS node prefix.", call. = FALSE)
  }

  random_effects <- formula_design$random_effects
  if(is.null(random_effects)){
    random_effects <- list()
  }
  n_rows <- nrow(formula_design$model_matrix)
  if(length(random_effects) > 0L){
    n_rows <- nrow(random_effects[[1L]]$model_matrix)
  }
  row_blocks <- .bt_JAGS_formula_random_covariance_row_blocks(
    row_blocks = row_blocks,
    n_rows = n_rows
  )

  if(length(random_effects) == 0L){
    return(structure(list(
      syntax = "",
      data = list(),
      row_blocks = row_blocks,
      lower_names = rep(NA_character_, length(row_blocks)),
      term_names = character()
    ), class = c("BayesTools_JAGS_random_marginal_covariance", "list")))
  }

  compile_modes <- vapply(
    random_effects,
    .bt_random_effect_term_compile_mode,
    character(1)
  )
  if(any(compile_modes != "marginalized")){
    stop(
      "Every formula random-effect block must be compiled as 'marginalized' ",
      "before constructing marginal JAGS covariance nodes.",
      call. = FALSE
    )
  }
  .bt_JAGS_bridge_validate_marginal_random_row_blocks(
    random_effects = random_effects,
    row_blocks = row_blocks,
    parameter = formula_design$parameter
  )

  syntax <- character()
  data <- list()
  term_names <- character(length(random_effects))

  for(term_index in seq_along(random_effects)){
    random_term <- random_effects[[term_index]]
    term_prefix <- paste0(prefix, "_term_", term_index)
    term_names[[term_index]] <- random_term$block_name
    term <- .bt_JAGS_formula_random_covariance_term(
      random_term = random_term,
      prefix = term_prefix,
      n_rows = n_rows
    )
    syntax <- c(syntax, term$syntax)
    data <- c(data, term$data)
  }

  lower_names <- character(length(row_blocks))
  for(block_index in seq_along(row_blocks)){
    rows <- row_blocks[[block_index]]
    pairs <- .bt_lower_triangle_pairs(rows)
    block_prefix <- paste0(prefix, "_block_", block_index)
    row_1_name <- paste0(block_prefix, "_row_1")
    row_2_name <- paste0(block_prefix, "_row_2")
    data[[row_1_name]] <- pairs$row_1
    data[[row_2_name]] <- pairs$row_2

    contribution_names <- character(length(random_effects))
    for(term_index in seq_along(random_effects)){
      random_term <- random_effects[[term_index]]
      term_prefix <- paste0(prefix, "_term_", term_index)
      contribution_name <- paste0(block_prefix, "_term_", term_index)
      contribution_names[[term_index]] <- contribution_name
      group_name <- paste0(contribution_name, "_group")
      data[[group_name]] <- .bt_JAGS_formula_random_covariance_group_values(
        random_term = random_term,
        row_1 = pairs$row_1,
        row_2 = pairs$row_2
      )
      n_columns <- random_term$n_columns
      syntax <- c(syntax, paste0(
        "for(l in 1:", nrow(pairs), "){\n",
        "  for(a in 1:", n_columns, "){\n",
        "    ", contribution_name, "_tmp[l,a] = ",
        term_prefix, "_basis[", row_1_name, "[l],a] * inprod(",
        term_prefix, "_cor[a,1:", n_columns, "], ",
        term_prefix, "_basis[", row_2_name, "[l],1:", n_columns, "])\n",
        "  }\n",
        "  ", contribution_name, "[l] = ", group_name,
        "[l] * sum(", contribution_name, "_tmp[l,1:", n_columns, "])\n",
        "}\n"
      ))
    }

    lower_name <- paste0(block_prefix, "_lower")
    lower_names[[block_index]] <- lower_name
    syntax <- c(syntax, paste0(
      "for(l in 1:", nrow(pairs), "){\n",
      "  ", lower_name, "[l] = ",
      if(length(contribution_names) == 1L){
        paste0(contribution_names[[1L]], "[l]")
      }else{
        paste0("sum(", block_prefix, "_terms[l,1:",
               length(contribution_names), "])")
      },
      "\n}\n"
    ))
    if(length(contribution_names) > 1L){
      assignments <- vapply(seq_along(contribution_names), function(term_index){
        paste0(
          "for(l in 1:", nrow(pairs), "){ ", block_prefix,
          "_terms[l,", term_index, "] = ", contribution_names[[term_index]],
          "[l] }\n"
        )
      }, character(1))
      syntax <- append(
        syntax,
        assignments,
        after = length(syntax) - 1L
      )
    }
  }

  structure(list(
    syntax = paste0(syntax, collapse = ""),
    data = data,
    row_blocks = row_blocks,
    lower_names = lower_names,
    term_names = term_names
  ), class = c("BayesTools_JAGS_random_marginal_covariance", "list"))
}


.bt_JAGS_formula_random_covariance_row_blocks <- function(row_blocks, n_rows){

  if(!is.list(row_blocks) || length(row_blocks) == 0L){
    stop("'row_blocks' must be a non-empty list of row indices.", call. = FALSE)
  }
  normalized <- lapply(row_blocks, function(rows){
    if(!is.numeric(rows) || length(rows) == 0L || anyNA(rows) ||
       any(!is.finite(rows)) || any(rows != as.integer(rows)) ||
       any(rows < 1L) || any(rows > n_rows) || anyDuplicated(rows)){
      stop("Each element of 'row_blocks' must contain unique fitted-row indices.",
           call. = FALSE)
    }
    as.integer(rows)
  })
  all_rows <- unlist(normalized, use.names = FALSE)
  if(length(all_rows) != n_rows || anyDuplicated(all_rows) ||
     !identical(sort(all_rows), seq_len(n_rows))){
    stop("'row_blocks' must partition every fitted row exactly once.",
         call. = FALSE)
  }
  normalized
}


.bt_JAGS_formula_random_covariance_term <- function(random_term, prefix,
                                                     n_rows){

  model_matrix <- random_term$model_matrix
  if(!is.numeric(model_matrix) || !is.matrix(model_matrix) ||
     nrow(model_matrix) != n_rows || any(!is.finite(model_matrix))){
    stop(
      "Random-effect block '", random_term$block_name,
      "' has invalid fitted model-matrix metadata.",
      call. = FALSE
    )
  }
  n_columns <- ncol(model_matrix)
  if(n_columns < 1L || !identical(as.integer(random_term$n_columns),
                                  as.integer(n_columns))){
    stop(
      "Random-effect block '", random_term$block_name,
      "' has inconsistent column metadata.",
      call. = FALSE
    )
  }

  data_name <- paste0(prefix, "_data")
  data <- stats::setNames(list(unname(model_matrix)), data_name)
  syntax <- character()

  for(column in seq_len(n_columns)){
    sd_expression <- .bt_JAGS_formula_random_covariance_sd_expression(
      random_term = random_term,
      column = column,
      row_index = "i"
    )
    syntax <- c(syntax, paste0(
      "for(i in 1:", n_rows, "){ ", prefix, "_basis[i,", column,
      "] = ", data_name, "[i,", column, "] * ", sd_expression, " }\n"
    ))
  }

  correlation <- .bt_JAGS_formula_random_covariance_correlation(
    random_term = random_term,
    prefix = prefix
  )
  syntax <- c(syntax, correlation$syntax)
  data <- c(data, correlation$data)

  list(syntax = syntax, data = data)
}


.bt_JAGS_formula_random_covariance_sd_expression <- function(
    random_term, column, row_index){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    stop(
      "Random-effect block '", random_term$block_name,
      "' is missing SD-binding metadata.",
      call. = FALSE
    )
  }
  if(!.bt_random_sd_binding_has_row_external_source(binding)){
    return(paste0(random_term$parameter_stem, "_xRE_STDx[", column, "]"))
  }

  source <- binding$source
  if(length(binding$sources_by_column) > 0L){
    source <- binding$sources_by_column[[column]]
  }
  source_expression <- .bt_random_sd_binding_source_jags_expression(
    source,
    row_index = row_index
  )
  factors <- if(identical(binding$application, "column")){
    binding$factors_by_column[[column]]
  }else{
    binding$factors
  }
  factor_expression <- .bt_random_sd_binding_factors_expression(factors)
  paste(
    c(source_expression,
      if(!identical(factor_expression, "1")) factor_expression),
    collapse = " * "
  )
}


.bt_JAGS_formula_random_covariance_correlation <- function(random_term,
                                                            prefix){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Marginal JAGS covariance compilation"
  )
  n_columns <- random_term$n_columns
  syntax <- character()
  data <- list()

  if(identical(structure, "us")){
    correlation <- random_term$correlation
    if(n_columns == 1L){
      syntax <- paste0(prefix, "_cor[1,1] = 1\n")
    }else if(is.list(correlation) &&
             is.character(correlation$correlation_name) &&
             length(correlation$correlation_name) == 1L &&
             nzchar(correlation$correlation_name)){
      syntax <- paste0(
        "for(a in 1:", n_columns, "){\n",
        "  for(c in 1:", n_columns, "){ ", prefix, "_cor[a,c] = ",
        correlation$correlation_name, "[a,c] }\n",
        "}\n"
      )
    }else{
      stop(
        "Unstructured random-effect block '", random_term$block_name,
        "' is missing correlation-matrix node metadata.",
        call. = FALSE
      )
    }
    return(list(syntax = syntax, data = data))
  }

  if(structure %in% c("id", "diag")){
    identity_name <- paste0(prefix, "_cor_data")
    data[[identity_name]] <- diag(n_columns)
    syntax <- paste0(
      "for(a in 1:", n_columns, "){\n",
      "  for(c in 1:", n_columns, "){ ", prefix, "_cor[a,c] = ",
      identity_name, "[a,c] }\n",
      "}\n"
    )
    return(list(syntax = syntax, data = data))
  }

  if(!structure %in% c("cs", "hcs", "ar1", "har", "car")){
    stop(
      "Random-effect structure '", structure, "' in block '",
      random_term$block_name,
      "' is unavailable for marginal JAGS covariance compilation.",
      call. = FALSE
    )
  }
  correlation <- random_term$correlation
  if(!is.list(correlation) || !is.character(correlation$rho_name) ||
     length(correlation$rho_name) != 1L || !nzchar(correlation$rho_name)){
    stop(
      "Structured random-effect block '", random_term$block_name,
      "' is missing scalar-correlation node metadata.",
      call. = FALSE
    )
  }
  distance <- if(structure %in% c("cs", "hcs")){
    outer(seq_len(n_columns), seq_len(n_columns), "!=") * 1
  }else if(identical(structure, "car")){
    time_values <- correlation$time_values
    if(!is.numeric(time_values) || length(time_values) != n_columns ||
       any(!is.finite(time_values))){
      stop(
        "Continuous autoregressive block '", random_term$block_name,
        "' is missing resolved time values.",
        call. = FALSE
      )
    }
    abs(outer(time_values, time_values, "-"))
  }else{
    abs(outer(seq_len(n_columns), seq_len(n_columns), "-"))
  }
  distance_name <- paste0(prefix, "_distance")
  data[[distance_name]] <- unname(distance)
  syntax <- paste0(
    "for(a in 1:", n_columns, "){\n",
    "  for(c in 1:", n_columns, "){ ", prefix, "_cor[a,c] = pow(",
    correlation$rho_name, ", ", distance_name, "[a,c]) }\n",
    "}\n"
  )
  list(syntax = syntax, data = data)
}


.bt_JAGS_formula_random_covariance_group_values <- function(
    random_term, row_1, row_2){

  group_map <- random_term$group_map
  n_rows <- nrow(random_term$model_matrix)
  if(!is.numeric(group_map) || length(group_map) != n_rows || anyNA(group_map) ||
     any(!is.finite(group_map)) || any(group_map != as.integer(group_map)) ||
     any(group_map < 1L)){
    stop(
      "Random-effect block '", random_term$block_name,
      "' has invalid grouping metadata.",
      call. = FALSE
    )
  }
  group_covariance <- random_term$group_covariance
  if(is.list(group_covariance) && identical(group_covariance$type, "known")){
    kernel <- group_covariance$kernel
    if(!is.numeric(kernel) || !is.matrix(kernel) ||
       nrow(kernel) != ncol(kernel) || any(!is.finite(kernel)) ||
       any(group_map > nrow(kernel))){
      stop(
        "Random-effect block '", random_term$block_name,
        "' has invalid known group-covariance metadata.",
        call. = FALSE
      )
    }
    return(unname(kernel[cbind(group_map[row_1], group_map[row_2])]))
  }
  as.numeric(group_map[row_1] == group_map[row_2])
}
