# formula helper functions
.remove_response        <- function(formula){
  # removes response from the expression
  # (prevents crash on formula evaluations)
  if(attr(stats::terms(formula), "response")  == 1){
    formula[2] <- NULL
  }
  return(formula)
}
.extract_expressions    <- function(formula){
  # extract all expressions from the formula

  return(.bt_extract_expression_bodies(.bt_formula_rhs(formula)))
}
.clean_from_expression  <- function(x){
  # expression to character

  return(sub("expression\\((.*)\\)", "\\1", x))
}
.remove_expressions     <- function(formula){
  # remove all expressions from the formula

  rhs_index <- .bt_formula_rhs_index(formula)
  rhs <- .bt_remove_expression_terms(formula[[rhs_index]])
  if(is.null(rhs)){
    rhs <- 1
  }
  formula[[rhs_index]] <- rhs

  return(formula)
}
.bt_formula_rhs_index <- function(formula){
  if(length(formula) == 3L) 3L else 2L
}
.bt_formula_rhs <- function(formula){
  formula[[.bt_formula_rhs_index(formula)]]
}
.bt_validate_fixed_formula_grammar <- function(formula){

  validate_expression <- function(expression){
    if(is.symbol(expression)){
      if(identical(as.character(expression), ".")){
        stop(
          "Unsupported fixed-formula term '.'. Dot expansion is not supported; ",
          "list each data-frame column explicitly.",
          call. = FALSE
        )
      }
      return(invisible(TRUE))
    }
    if(is.numeric(expression) && length(expression) == 1L &&
       is.finite(expression) && expression %in% c(0, 1)){
      return(invisible(TRUE))
    }
    if(is.call(expression)){
      call_name <- if(is.symbol(expression[[1L]])){
        as.character(expression[[1L]])
      }else{
        ""
      }
      if(identical(call_name, "^") && length(expression) == 3L){
        degree <- expression[[3L]]
        if(is.numeric(degree) && length(degree) == 1L && is.finite(degree) &&
           degree >= 2 && degree <= .Machine$integer.max && degree == floor(degree)){
          validate_expression(expression[[2L]])
          return(invisible(TRUE))
        }
      }
      if(call_name %in% c("+", "-", "*", ":", "/", "^", "(")){
        for(argument in as.list(expression)[-1L]){
          validate_expression(argument)
        }
        return(invisible(TRUE))
      }

      expression_label <- .bt_deparse_expr(expression)
      if(identical(call_name, "offset")){
        stop(
          "Unsupported fixed-formula call '", expression_label,
          "'. offset() is not supported; use expression(...) for an explicit ",
          "JAGS-scale offset.",
          call. = FALSE
        )
      }
      stop(
        "Unsupported fixed-formula call '", expression_label,
        "'. Create the transformed value as a data-frame column and reference ",
        "that column by name.",
        call. = FALSE
      )
    }

    stop(
      "Unsupported fixed-formula expression '",
      .bt_deparse_expr(expression),
      "'. Use literal data-frame column names.",
      call. = FALSE
    )
  }

  validate_expression(.bt_formula_rhs(formula))
  invisible(TRUE)
}
.bt_validate_formula_replay_grammar <- function(formula){

  if(!inherits(formula, "formula")){
    stop("'formula' must be a formula.", call. = FALSE)
  }

  validation_formula <- .remove_response(formula)
  validation_formula <- .remove_expressions(validation_formula)
  validation_formula <- .remove_random_effects(validation_formula)
  .bt_validate_fixed_formula_grammar(validation_formula)
}
.bt_is_expression_call <- function(x){
  is.call(x) && identical(as.character(x[[1L]]), "expression")
}
.bt_contains_expression_call <- function(x){
  if(.bt_is_expression_call(x)){
    return(TRUE)
  }
  if(is.call(x) || is.pairlist(x)){
    return(any(vapply(as.list(x)[-1], .bt_contains_expression_call, logical(1))))
  }

  FALSE
}
.bt_extract_expression_bodies <- function(x){
  if(.bt_is_expression_call(x)){
    return(lapply(as.list(x)[-1], function(expression_body){
      paste0(deparse(expression_body), collapse = " ")
    }))
  }
  if(is.call(x) || is.pairlist(x)){
    return(unlist(lapply(as.list(x)[-1], .bt_extract_expression_bodies), recursive = FALSE))
  }

  list()
}
.bt_formula_expression_label <- function(expression_body){

  if(is.character(expression_body)){
    return(paste(expression_body, collapse = " "))
  }
  paste(deparse(expression_body), collapse = " ")
}
.bt_formula_expression_stop <- function(expression_label, detail = NULL){

  stop(
    "expression() term '", expression_label, "' is not replayable",
    if(is.null(detail)) "." else paste0(": ", detail, "."),
    " Use numeric constants, replayable data, 'i' row indexing, sampled scalar ",
    "or one-dimensional indexed parameters, arithmetic operators, or abs(), ",
    "exp(), log(), and sqrt().",
    call. = FALSE
  )
}
.bt_validate_formula_expression_node <- function(node, expression_label){

  if(is.numeric(node) && length(node) == 1L && is.finite(node)){
    return(invisible(TRUE))
  }
  if(is.symbol(node)){
    symbol <- as.character(node)
    if(!identical(symbol, "i")){
      valid_name <- tryCatch({
        .bt_check_jags_node_name(symbol, "expression symbol")
        TRUE
      }, error = function(e) FALSE)
      if(!isTRUE(valid_name)){
        .bt_formula_expression_stop(
          expression_label,
          paste0("invalid JAGS symbol '", symbol, "'")
        )
      }
    }
    return(invisible(TRUE))
  }
  if(!is.call(node) || !is.symbol(node[[1L]])){
    .bt_formula_expression_stop(expression_label, "unsupported syntax")
  }

  call_name <- as.character(node[[1L]])
  arguments <- as.list(node)[-1L]
  valid_arity <- switch(
    call_name,
    "(" = length(arguments) == 1L,
    "+" = length(arguments) %in% c(1L, 2L),
    "-" = length(arguments) %in% c(1L, 2L),
    "*" = length(arguments) == 2L,
    "/" = length(arguments) == 2L,
    "^" = length(arguments) == 2L,
    "[" = length(arguments) >= 2L,
    "abs" = length(arguments) == 1L,
    "exp" = length(arguments) == 1L,
    "log" = length(arguments) == 1L,
    "sqrt" = length(arguments) == 1L,
    FALSE
  )
  if(!isTRUE(valid_arity)){
    .bt_formula_expression_stop(
      expression_label,
      paste0("unsupported call '", call_name, "'")
    )
  }
  for(argument in arguments){
    .bt_validate_formula_expression_node(
      argument,
      expression_label
    )
  }
  invisible(TRUE)
}
.bt_parse_formula_expression <- function(expression_body){

  expression_label <- .bt_formula_expression_label(expression_body)
  parsed <- tryCatch(
    parse(text = expression_label),
    error = function(e) e
  )
  if(inherits(parsed, "error") || length(parsed) != 1L){
    .bt_formula_expression_stop(expression_label, "invalid expression syntax")
  }
  .bt_validate_formula_expression_node(
    parsed[[1L]],
    expression_label
  )
  parsed[[1L]]
}
.bt_formula_expression_symbols <- function(node){

  if(is.symbol(node)){
    return(as.character(node))
  }
  if(!is.call(node)){
    return(character())
  }
  unique(unlist(lapply(
    as.list(node)[-1L],
    .bt_formula_expression_symbols
  ), use.names = FALSE))
}
.bt_validate_formula_expression_parameter_index <- function(node,
                                                            parameter_names,
                                                            expression_label){

  if(!is.call(node)){
    return(invisible(TRUE))
  }
  call_name <- if(is.symbol(node[[1L]])) as.character(node[[1L]]) else ""
  arguments <- as.list(node)[-1L]
  if(identical(call_name, "[") && is.symbol(arguments[[1L]]) &&
     as.character(arguments[[1L]]) %in% parameter_names &&
     length(arguments) != 2L){
    .bt_formula_expression_stop(
      expression_label,
      paste0(
        "sampled parameter '", as.character(arguments[[1L]]),
        "' must use exactly one index"
      )
    )
  }
  for(argument in arguments){
    .bt_validate_formula_expression_parameter_index(
      argument,
      parameter_names,
      expression_label
    )
  }
  invisible(TRUE)
}
.bt_formula_expression_specs <- function(expressions, data_names = character(),
                                         parameter_names = character(),
                                         allow_unresolved = FALSE){

  lapply(expressions, function(expression_body){
    label <- .bt_formula_expression_label(expression_body)
    parsed <- .bt_parse_formula_expression(label)
    dependencies <- setdiff(.bt_formula_expression_symbols(parsed), "i")
    overlap <- intersect(
      dependencies,
      intersect(data_names, parameter_names)
    )
    if(length(overlap) > 0L){
      .bt_formula_expression_stop(
        label,
        paste0(
          "dependency ", paste0("'", overlap, "'", collapse = ", "),
          " is both data and a parameter"
        )
      )
    }
    data_dependencies <- dependencies[dependencies %in% data_names]
    parameter_dependencies <- dependencies[dependencies %in% parameter_names]
    .bt_validate_formula_expression_parameter_index(
      parsed,
      parameter_dependencies,
      label
    )
    unresolved_dependencies <- setdiff(
      dependencies,
      c(data_dependencies, parameter_dependencies)
    )
    if(length(unresolved_dependencies) > 0L && !isTRUE(allow_unresolved)){
      .bt_formula_expression_stop(
        label,
        paste0(
          "unknown replay dependency ",
          paste0("'", unresolved_dependencies, "'", collapse = ", ")
        )
      )
    }
    list(
      label = label,
      parsed = parsed,
      dependencies = dependencies,
      data_dependencies = data_dependencies,
      parameter_dependencies = parameter_dependencies,
      unresolved_dependencies = unresolved_dependencies
    )
  })
}
.bt_formula_expression_specs_valid <- function(specs){

  is.list(specs) && all(vapply(specs, function(spec){
    is.list(spec) && is.character(spec$label) && length(spec$label) == 1L &&
      (is.language(spec$parsed) ||
       (is.numeric(spec$parsed) && length(spec$parsed) == 1L &&
        is.finite(spec$parsed))) &&
      is.character(spec$dependencies) &&
      is.character(spec$data_dependencies) &&
      is.character(spec$parameter_dependencies) &&
      is.character(spec$unresolved_dependencies)
  }, logical(1)))
}
.bt_validate_formula_expressions <- function(expressions, data,
                                             parameter_names = character(),
                                             allow_unresolved = FALSE){

  if(length(expressions) == 0L){
    return(list())
  }
  data_names <- names(data)
  if(is.null(data_names) || any(!nzchar(data_names)) || anyDuplicated(data_names)){
    stop("Formula expression source data must have unique, nonempty names.",
         call. = FALSE)
  }
  if("i" %in% data_names){
    stop(
      "Formula expression source data cannot contain a column named 'i' ",
      "because it is reserved for JAGS-style row indexing.",
      call. = FALSE
    )
  }
  .bt_formula_expression_specs(
    expressions = expressions,
    data_names = data_names,
    parameter_names = parameter_names,
    allow_unresolved = allow_unresolved
  )
}
.bt_formula_expression_data_list <- function(data, context){

  if(is.null(data)){
    return(list())
  }
  if(is.data.frame(data)){
    return(as.list(data))
  }
  if(is.list(data)){
    return(data)
  }
  stop(context, " must be a data.frame or named list.", call. = FALSE)
}
.bt_formula_expression_validate_data_value <- function(value, name, n_rows,
                                                       context){

  if(!is.numeric(value) && !is.integer(value) && !is.logical(value)){
    stop(context, " data dependency '", name, "' must be numeric or logical.",
         call. = FALSE)
  }
  if(anyNA(value) || any(!is.finite(value))){
    stop(context, " data dependency '", name, "' must be finite.", call. = FALSE)
  }
  dimensions <- dim(value)
  row_aligned <- if(is.null(dimensions)){
    length(value) %in% c(1L, n_rows)
  }else{
    length(dimensions) > 0L && dimensions[1L] == n_rows
  }
  if(!isTRUE(row_aligned)){
    stop(
      context, " data dependency '", name,
      "' must be scalar or have ", n_rows, " rows.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
.bt_formula_expression_data <- function(specs, formula_data, model_data = NULL,
                                        n_rows, context){

  dependencies <- unique(unlist(lapply(
    specs,
    `[[`,
    "data_dependencies"
  ), use.names = FALSE))
  if(length(dependencies) == 0L){
    return(list())
  }
  formula_data <- .bt_formula_expression_data_list(formula_data, context)
  model_data <- .bt_formula_expression_data_list(model_data, context)
  overlap <- intersect(names(formula_data), names(model_data))
  for(name in intersect(overlap, dependencies)){
    if(!identical(formula_data[[name]], model_data[[name]])){
      stop(
        context, " received conflicting formula and model data for '", name,
        "'.",
        call. = FALSE
      )
    }
  }
  combined <- formula_data
  combined[setdiff(names(model_data), names(combined))] <-
    model_data[setdiff(names(model_data), names(combined))]
  missing <- setdiff(dependencies, names(combined))
  if(length(missing) > 0L){
    stop(
      context, " is missing data dependency ",
      paste0("'", missing, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  out <- combined[dependencies]
  for(name in names(out)){
    .bt_formula_expression_validate_data_value(
      out[[name]],
      name,
      n_rows,
      context
    )
  }
  out
}
.bt_formula_expression_merge_data <- function(data, stored_data = NULL,
                                              context = "Formula expression"){

  out <- .bt_formula_expression_data_list(data, context)
  stored <- .bt_formula_expression_data_list(stored_data, context)
  overlap <- intersect(names(out), names(stored))
  for(name in overlap){
    if(!identical(out[[name]], stored[[name]])){
      stop(context, " received conflicting values for data dependency '", name,
           "'.", call. = FALSE)
    }
  }
  out[setdiff(names(stored), names(out))] <- stored[setdiff(names(stored), names(out))]
  out
}
.bt_formula_expression_merge_jags_data <- function(generated, supplied,
                                                   context){

  if(is.null(generated)){
    generated <- list()
  }
  if(is.null(supplied)){
    supplied <- list()
  }
  overlap <- intersect(names(generated), names(supplied))
  for(name in overlap){
    if(!identical(generated[[name]], supplied[[name]])){
      stop(context, " contains conflicting JAGS data named '", name, "'.",
           call. = FALSE)
    }
  }
  generated[setdiff(names(supplied), names(generated))] <-
    supplied[setdiff(names(supplied), names(generated))]
  generated
}
.bt_formula_expression_finalize_design <- function(design, formula_data,
                                                   model_data,
                                                   parameter_names,
                                                   forbidden_parameters,
                                                   context){

  expressions <- design$transformed_terms
  if(length(expressions) == 0L){
    design$expression_specs <- list()
    design$expression_data <- list()
    return(design)
  }
  formula_data_list <- .bt_formula_expression_data_list(formula_data, context)
  model_data_list <- .bt_formula_expression_data_list(model_data, context)
  data_names <- unique(c(names(formula_data_list), names(model_data_list)))
  parameter_names <- unique(.bt_parameter_coordinates_base(parameter_names))
  specs <- .bt_formula_expression_specs(
    expressions = expressions,
    data_names = data_names,
    parameter_names = parameter_names,
    allow_unresolved = FALSE
  )
  expression_parameters <- unique(unlist(lapply(
    specs,
    `[[`,
    "parameter_dependencies"
  ), use.names = FALSE))
  forbidden <- intersect(expression_parameters, forbidden_parameters)
  if(length(forbidden) > 0L){
    stop(
      context, " cannot replay formula-output dependency ",
      paste0("'", forbidden, "'", collapse = ", "),
      "; cross-formula and self-referential expression dependencies are not ",
      "supported.",
      call. = FALSE
    )
  }
  design$expression_specs <- specs
  design$expression_data <- .bt_formula_expression_data(
    specs = specs,
    formula_data = formula_data_list,
    model_data = model_data_list,
    n_rows = nrow(design$source_data),
    context = context
  )
  design
}
.bt_formula_expression_sample_names <- function(samples){

  if(is.matrix(samples) || is.data.frame(samples)){
    return(colnames(samples))
  }
  names(samples)
}
.bt_formula_expression_parameter_roots <- function(samples,
                                                   parameters = NULL){

  sample_names <- .bt_formula_expression_sample_names(samples)
  sample_roots <- if(is.null(sample_names)) character() else
    .bt_parameter_coordinates_base(sample_names)
  unique(c(sample_roots, names(parameters)))
}
.bt_formula_expression_indexed_draws <- function(values, indices){

  out <- list(
    values = unname(as.matrix(values)),
    indices = indices,
    size = max(indices)
  )
  class(out) <- c("BayesTools_formula_expression_indexed_draws", "list")
  out
}
.bt_formula_expression_parameter_value <- function(draws, draw){

  if(inherits(draws, "BayesTools_formula_expression_indexed_draws")){
    value <- rep.int(NA_real_, draws$size)
    value[draws$indices] <- draws$values[draw, ]
    return(value)
  }
  unname(draws[draw, ])
}
.bt_formula_expression_resolve_specs <- function(expressions, data, samples,
                                                 parameters = NULL){

  if(.bt_formula_expression_specs_valid(expressions)){
    if(any(vapply(expressions, function(spec){
      length(spec$unresolved_dependencies) > 0L
    }, logical(1)))){
      expressions <- vapply(expressions, `[[`, character(1), "label")
    }else{
      return(expressions)
    }
  }
  data_names <- names(.bt_formula_expression_data_list(
    data,
    "Formula expression source data"
  ))
  .bt_formula_expression_specs(
    expressions = expressions,
    data_names = data_names,
    parameter_names = .bt_formula_expression_parameter_roots(
      samples,
      parameters
    ),
    allow_unresolved = FALSE
  )
}
.bt_formula_expression_validate_replay_data <- function(specs, data, n_rows,
                                                        context){

  dependencies <- unique(unlist(lapply(
    specs,
    `[[`,
    "data_dependencies"
  ), use.names = FALSE))
  missing <- setdiff(dependencies, names(data))
  if(length(missing) > 0L){
    stop(
      context, " is missing expression data dependency ",
      paste0("'", missing, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  for(name in dependencies){
    .bt_formula_expression_validate_data_value(
      data[[name]],
      name,
      n_rows,
      context
    )
  }
  invisible(TRUE)
}
.bt_formula_expression_parameter_draws <- function(samples, parameters,
                                                   parameter, n_draws,
                                                   context){

  if(!is.null(parameters[[parameter]])){
    if(n_draws != 1L){
      stop(context, " internal parameter reconstruction requires one draw.",
           call. = FALSE)
    }
    value <- parameters[[parameter]]
    if(!is.numeric(value) || is.matrix(value) || length(dim(value)) > 1L){
      stop(context, " parameter '", parameter,
           "' must be scalar or one-dimensional.", call. = FALSE)
    }
    return(matrix(unname(value), nrow = 1L))
  }

  if(is.matrix(samples) || is.data.frame(samples)){
    sample_names <- colnames(samples)
    direct <- !is.null(sample_names) && parameter %in% sample_names
    indexed <- JAGS_indexed_parameter_matrix(samples, parameter)
    if(isTRUE(direct) && !is.null(indexed)){
      stop(context, " parameter '", parameter,
           "' has both scalar and indexed coordinates.", call. = FALSE)
    }
    if(isTRUE(direct)){
      return(matrix(samples[, parameter], ncol = 1L))
    }
    if(!is.null(indexed)){
      indices <- .JAGS_indexed_parameter_indices(colnames(indexed), parameter)
      return(.bt_formula_expression_indexed_draws(indexed, indices))
    }
  }else{
    sample_names <- names(samples)
    direct <- !is.null(sample_names) && parameter %in% sample_names
    indexed <- JAGS_indexed_parameter_vector(samples, parameter)
    if(isTRUE(direct) && length(indexed) > 0L){
      stop(context, " parameter '", parameter,
           "' has both scalar and indexed coordinates.", call. = FALSE)
    }
    if(isTRUE(direct)){
      return(matrix(samples[[parameter]], nrow = 1L))
    }
    if(length(indexed) > 0L){
      indices <- .JAGS_indexed_parameter_indices(names(indexed), parameter)
      return(.bt_formula_expression_indexed_draws(
        matrix(unname(indexed), nrow = 1L),
        indices
      ))
    }
  }
  stop(
    context, " cannot reconstruct expression parameter '", parameter,
    "' from posterior or bridge coordinates.",
    call. = FALSE
  )
}
.bt_formula_expression_index <- function(x, ...){

  indices <- list(...)
  if(length(indices) <= 1L){
    return(base::`[`(x, ...))
  }
  dimensions <- dim(x)
  if(length(dimensions) != length(indices)){
    stop("Indexed expression must supply one index for every array dimension.",
         call. = FALSE)
  }
  index_lengths <- lengths(indices)
  common_length <- max(index_lengths)
  if(common_length == 0L ||
     any(!index_lengths %in% c(1L, common_length))){
    stop("Indexed expression indices must have length 1 or a common nonzero length.",
         call. = FALSE)
  }
  for(index_i in seq_along(indices)){
    index <- indices[[index_i]]
    if((!is.numeric(index) && !is.logical(index)) || is.complex(index) ||
       !is.null(dim(index)) || any(!is.finite(index)) ||
       any(index <= 0 | index != floor(index) | index > dimensions[[index_i]])){
      stop("Indexed expression coordinates must be finite positive integers within their array dimensions.",
           call. = FALSE)
    }
    index <- as.numeric(index)
    indices[[index_i]] <- if(length(index) == 1L){
      rep.int(index, common_length)
    }else{
      index
    }
  }
  x[do.call(cbind, indices)]
}
.bt_formula_expression_eval <- function(spec, data, n_rows,
                                        parameter_values = list(), context,
                                        mode = "row"){

  env_data <- data
  if(identical(mode, "row")) env_data[["i"]] <- seq_len(n_rows)
  env_data[names(parameter_values)] <- parameter_values
  env_data[["["]] <- if(identical(mode, "point")) .bt_formula_point_index else .bt_formula_expression_index
  value <- tryCatch(
    eval(spec$parsed, envir = list2env(env_data, parent = baseenv())),
    error = function(e) e
  )
  if(inherits(value, "error")){
    if(identical(mode, "point") && inherits(value, "BayesTools_formula_point_unavailable")) stop(value)
    stop(
      context, " '", spec$label, "' could not be evaluated: ",
      conditionMessage(value),
      call. = FALSE
    )
  }
  if(!is.numeric(value) && !is.integer(value) && !is.logical(value)){
    stop(context, " '", spec$label, "' must evaluate to numeric values.",
         call. = FALSE)
  }
  if(identical(mode, "point") && !is.numeric(value)){
    .bt_formula_point_stop(context, "nonnumeric_point_expression", "the result must be numeric")
  }
  value <- as.numeric(value)
  if(identical(mode, "point") && length(value) != 1L){
    .bt_formula_point_stop(context, "nonscalar_point_expression", "the result must be one finite scalar")
  }
  if(length(value) == 1L){
    value <- rep.int(value, n_rows)
  }
  if(length(value) != n_rows){
    stop(
      context, " '", spec$label,
      "' must evaluate to length 1 or ", n_rows, ".",
      call. = FALSE
    )
  }
  if(any(!is.finite(value))){
    if(identical(mode, "point")) .bt_formula_point_stop(context,
      "nonfinite_point_expression", "the expression produced a nonfinite value")
    stop(context, " '", spec$label, "' produced non-finite values.",
         call. = FALSE)
  }
  value
}
.bt_formula_expression_row_values <- function(expressions, data, n_rows,
                                              context = "Formula expression",
                                              samples = NULL,
                                              parameters = NULL){

  if(length(expressions) == 0L){
    return(rep.int(0, n_rows))
  }
  check_int(n_rows, "n_rows", lower = 0, allow_NA = FALSE)
  if(n_rows == 0L){
    return(numeric())
  }
  env_data <- .bt_formula_expression_data_list(data, context)
  if("i" %in% names(env_data)){
    stop(
      context,
      " source data cannot contain a column named 'i' because it is reserved ",
      "for JAGS-style row indexing inside expression() terms.",
      call. = FALSE
    )
  }
  specs <- .bt_formula_expression_resolve_specs(
    expressions,
    env_data,
    samples,
    parameters
  )
  .bt_formula_expression_validate_replay_data(
    specs,
    env_data,
    n_rows,
    context
  )

  total <- rep.int(0, n_rows)
  for(spec in specs){
    parameter_values <- lapply(spec$parameter_dependencies, function(parameter){
      values <- .bt_formula_expression_parameter_draws(
        samples,
        parameters,
        parameter,
        n_draws = 1L,
        context = context
      )
      .bt_formula_expression_parameter_value(values, 1L)
    })
    names(parameter_values) <- spec$parameter_dependencies
    total <- total + .bt_formula_expression_eval(
      spec,
      env_data,
      n_rows,
      parameter_values,
      context
    )
  }

  total
}

.bt_formula_expression_contribution_matrix <- function(expressions, data,
                                                       n_rows, n_draws,
                                                       context = "Formula expression",
                                                       samples = NULL,
                                                       parameters = NULL){

  env_data <- .bt_formula_expression_data_list(data, context)
  specs <- .bt_formula_expression_resolve_specs(
    expressions,
    env_data,
    samples,
    parameters
  )
  .bt_formula_expression_validate_replay_data(
    specs,
    env_data,
    n_rows,
    context
  )
  parameter_names <- unique(unlist(lapply(
    specs,
    `[[`,
    "parameter_dependencies"
  ), use.names = FALSE))
  parameter_draws <- lapply(parameter_names, function(parameter){
    .bt_formula_expression_parameter_draws(
      samples,
      parameters,
      parameter,
      n_draws,
      context
    )
  })
  names(parameter_draws) <- parameter_names

  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  for(spec in specs){
    if(length(spec$parameter_dependencies) == 0L){
      value <- .bt_formula_expression_eval(
        spec,
        env_data,
        n_rows,
        context = context
      )
      output <- output + matrix(value, nrow = n_rows, ncol = n_draws)
      next
    }
    for(draw in seq_len(n_draws)){
      parameter_values <- lapply(
        spec$parameter_dependencies,
        function(parameter){
          .bt_formula_expression_parameter_value(
            parameter_draws[[parameter]],
            draw
          )
        }
      )
      names(parameter_values) <- spec$parameter_dependencies
      output[, draw] <- output[, draw] + .bt_formula_expression_eval(
        spec,
        env_data,
        n_rows,
        parameter_values,
        context
      )
    }
  }
  output
}

.bt_resolve_formula_expression_terms <- function(formula, fitted_design,
                                                 replay_fitted_formula){

  if(isTRUE(replay_fitted_formula)){
    terms <- fitted_design$expression_specs
    if(is.null(terms)){
      terms <- fitted_design$transformed_terms
    }
    if(is.null(terms)){
      return(list())
    }
    return(terms)
  }
  if(is.null(formula)){
    return(list())
  }
  .extract_expressions(formula)
}

.bt_remove_expression_terms <- function(x){
  if(.bt_is_expression_call(x)){
    return(NULL)
  }

  if(is.call(x)){
    call_name <- if(is.symbol(x[[1L]])) as.character(x[[1L]]) else ""
    if(call_name == "+" && length(x) == 3L){
      lhs <- .bt_remove_expression_terms(x[[2L]])
      rhs <- .bt_remove_expression_terms(x[[3L]])
      if(is.null(lhs)){
        return(rhs)
      }
      if(is.null(rhs)){
        return(lhs)
      }
      return(call("+", lhs, rhs))
    }
    if(call_name == "-" && length(x) == 3L){
      if(.bt_contains_expression_call(x[[3L]])){
        stop("expression() terms must be additive formula terms.", call. = FALSE)
      }
      lhs <- .bt_remove_expression_terms(x[[2L]])
      rhs <- x[[3L]]
      if(is.null(lhs)){
        lhs <- 1
      }
      return(call("-", lhs, rhs))
    }
    if(.bt_contains_expression_call(x)){
      stop("expression() terms must be additive formula terms.", call. = FALSE)
    }
  }

  x
}
.bt_formula_random_terms <- function(formula){

  if(inherits(formula, "BayesTools_random_effects")){
    return(formula$terms)
  }
  random_terms <- attr(formula, "random_terms", exact = TRUE)
  if(is.list(random_terms)){
    return(random_terms)
  }

  .bt_parse_random_effects(formula)$terms
}
.bt_formula_random_formula <- function(formula){

  if(inherits(formula, "BayesTools_random_effects")){
    return(formula$formula)
  }

  formula
}
.bt_formula_preserve_random_terms <- function(formula, random_terms){

  if(is.list(random_terms)){
    attr(formula, "random_terms") <- random_terms
  }
  formula
}
.has_random_effects     <- function(formula){
  return(length(.bt_formula_random_terms(formula)) > 0L)
}
.remove_random_effects  <- function(formula){
  return(.bt_fixed_formula(formula))
}
.get_grouping_factor    <- function(x){
  has_grouping            <- grepl("\\|", x)
  grouping                <- rep("", length(x))
  grouping[has_grouping]  <- trimws(sub(".*\\|\\s*", "", x[has_grouping]))
  return(grouping)
}
.JAGS_formula_default_prior_names <- function(){
  c("__default_continuous", "__default_factor")
}
.JAGS_formula_is_lazy_default_prior <- function(x){
  is.function(x) && typeof(x) == "closure" && length(formals(x)) == 0L
}
.JAGS_formula_check_prior_list <- function(prior_list){

  default_prior_names <- .JAGS_formula_default_prior_names()
  prior_names         <- names(prior_list)

  for(i in seq_along(prior_list)){
    prior_name  <- if(is.null(prior_names)) "" else prior_names[[i]]
    prior_value <- prior_list[[i]]

    if(prior_name %in% default_prior_names){
      if(is.prior(prior_value) || .JAGS_formula_is_lazy_default_prior(prior_value)){
        next
      }
      if(is.function(prior_value)){
        stop(
          paste0(
            "The 'prior_list[[\"", prior_name, "\"]]' entry must be a prior object ",
            "or a zero-argument function returning a prior object."
          ),
          call. = FALSE
        )
      }
    }

    if(!is.prior(prior_value)){
      stop("'prior_list' must be a list of priors.", call. = FALSE)
    }
  }

  return()
}
.JAGS_formula_resolve_default_prior <- function(default_prior, default_name){

  if(.JAGS_formula_is_lazy_default_prior(default_prior)){
    default_prior <- default_prior()
    if(!is.prior(default_prior)){
      stop(
        paste0(
          "The 'prior_list[[\"", default_name, "\"]]' lazy default must return ",
          "a BayesTools prior object."
        ),
        call. = FALSE
      )
    }
  }

  return(default_prior)
}
.JAGS_formula_canonicalize_none_prior <- function(prior_object){

  if(!is.prior.none(prior_object)){
    return(prior_object)
  }

  output <- prior(
    "point",
    list(location = 0),
    prior_weights = .prior_model_weight(prior_object)
  )
  prior_attributes <- attributes(prior_object)
  metadata_names <- setdiff(names(prior_attributes), c("names", "class"))
  for(metadata_name in metadata_names){
    attr(output, metadata_name) <- prior_attributes[[metadata_name]]
  }

  output
}
.bt_validate_formula_term_priors <- function(prior_list, model_terms,
                                             model_terms_type){

  # The linear predictor enters the intercept unscaled; a 'multiply_by'
  # multiplier applies only to the coefficients of formula terms.
  if(!is.null(attr(prior_list[["intercept"]], "multiply_by", exact = TRUE))){
    stop(
      "The intercept prior 'prior_list[[\"intercept\"]]' has a 'multiply_by' ",
      "attribute, but intercept priors cannot be scaled. Remove 'multiply_by' ",
      "from the intercept prior; it is supported only for the priors of ",
      "formula terms.",
      call. = FALSE
    )
  }

  for(model_term in model_terms){
    this_prior <- prior_list[[model_term]]
    term_type <- model_terms_type[[model_term]]
    factor_prior <- .bt_prior_is_factor_family(this_prior)

    if(factor_prior){
      .validate_centered_factor_prior(
        this_prior,
        paste0("prior_list[[\"", model_term, "\"]]")
      )
    }
    if(identical(term_type, "factor") && !factor_prior){
      stop(
        "Unsupported prior distribution defined for '", model_term,
        "' factor variable. See '?prior_factor' for details.",
        call. = FALSE
      )
    }
    if(identical(term_type, "continuous") &&
       (factor_prior || is.prior.discrete(this_prior) ||
        is.prior.PET(this_prior) || is.prior.PEESE(this_prior) ||
        is.prior.weightfunction(this_prior))){
      stop(
        "Unsupported prior distribution defined for '", model_term,
        "' continuous variable. See '?prior' for details.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
# A formula term codes a factor by level indicators instead of its contrast
# when the term without that factor is not in the formula (the value 2 in the
# terms "factors" attribute that stats::model.matrix() follows), e.g., `g` in
# the `g:x` term of `~ g + g:x`. The term then has one coefficient per level
# of that factor. Mean-difference and orthonormal priors are defined on
# contrast coefficients, so on such a term they would be applied to the level
# coefficients as independent priors; they are rejected instead. A point mass
# at zero (such as the null-hypothesis spike of a model-averaged term) fixes
# every coefficient at zero in any basis and is accepted.
.bt_validate_indicator_coded_factor_priors <- function(prior_list, model_terms,
                                                       model_terms_type,
                                                       predictors_type,
                                                       term_factors){

  factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  for(model_term in model_terms[model_terms_type[model_terms] == "factor"]){
    this_prior <- prior_list[[model_term]]
    contrast <- if(is.prior.meandif(this_prior)){
      "meandif"
    }else if(is.prior.orthonormal(this_prior)){
      "orthonormal"
    }
    zero_prior <- if(is.prior.mixture(this_prior)){
      all(vapply(this_prior, .posterior_atoms_is_zero_point, logical(1)))
    }else{
      .posterior_atoms_is_zero_point(this_prior)
    }
    if(is.null(contrast) || zero_prior){
      next
    }

    components <- rownames(term_factors)[term_factors[, model_term] > 0]
    indicator_factors <- components[
      components %in% factor_predictors &
        term_factors[components, model_term] == 2
    ]
    if(length(indicator_factors) == 0L){
      next
    }

    missing_terms <- vapply(indicator_factors, function(indicator_factor){
      paste0(setdiff(components, indicator_factor), collapse = ":")
    }, character(1))
    stop(
      "The '", contrast, "' prior of the factor term '", model_term,
      "' is unavailable: the formula has no term ",
      paste0("'", missing_terms, "'", collapse = " or "),
      ", so '", model_term, "' codes ",
      paste0("'", indicator_factors, "'", collapse = " and "),
      " by level indicators and has one coefficient per level instead of '",
      contrast, "' contrast coefficients. Add ",
      paste0("'", missing_terms, "'", collapse = " and "),
      " to the formula to keep the '", contrast, "' contrast, or use ",
      "prior_factor(contrast = \"independent\") for one independent ",
      "coefficient per level.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
.bt_validate_formula_reconstruction_prior <- function(prior_object,
                                                       prior_name){

  if(is.prior.point(prior_object) ||
     is.prior.factor(prior_object) ||
     is.prior.simple(prior_object)){
    return(invisible(TRUE))
  }

  stop(
    "Unsupported formula reconstruction prior for '", prior_name,
    "'. Formula metadata must contain a canonical simple or factor prior.",
    call. = FALSE
  )
}
.bt_validate_formula_log_intercept_prior <- function(prior_list,
                                                      parameter = NULL){

  intercept_name <- if(is.null(parameter)){
    "intercept"
  }else{
    paste0(parameter, "_intercept")
  }
  if(!intercept_name %in% names(prior_list)){
    stop(
      "A formula using log(intercept) must define a prior for '",
      intercept_name, "'.",
      call. = FALSE
    )
  }

  .validate_strictly_positive_prior(
    prior_list[[intercept_name]],
    if(is.null(parameter)){
      "prior_list[[\"intercept\"]]"
    }else{
      paste0("formula_prior_list[[\"", intercept_name, "\"]]")
    }
  )
}
#' @title Add an Intercept to a Formula
#'
#' @description Converts a no-intercept formula to the corresponding formula
#' with an intercept while preserving the formula environment. Additive,
#' parenthesized, and unary-plus no-intercept encodings such as \code{- 1},
#' \code{+ 0}, and \code{0 +} are removed without editing transformed calls
#' such as \code{I(x - 1)} or \code{offset(x - 1)}.
#'
#' @param formula a formula object.
#'
#' @return A formula object with an intercept.
#'
#' @export
formula_add_intercept <- function(formula){

  if(!inherits(formula, "formula")){
    stop("'formula' must be a formula.", call. = FALSE)
  }

  if(attr(stats::terms(formula), "intercept") == 1L){
    return(formula)
  }

  formula_env   <- environment(formula)
  formula_attrs <- attributes(formula)
  rhs_index     <- if(length(formula) == 3L) 3L else 2L
  rhs           <- .formula_strip_no_intercept(formula[[rhs_index]])

  if(is.null(rhs)){
    rhs <- 1
  }

  out <- formula
  out[[rhs_index]] <- rhs
  environment(out) <- formula_env

  for(attribute in setdiff(names(formula_attrs), c("class", ".Environment", "names"))){
    attr(out, attribute) <- formula_attrs[[attribute]]
  }

  if(attr(stats::terms(out), "intercept") == 0L){
    out[[rhs_index]] <- call("+", 1, out[[rhs_index]])
    environment(out) <- formula_env
  }

  return(out)
}
.formula_strip_no_intercept <- function(expr){

  if(is.call(expr) && length(expr) == 2L &&
     (identical(expr[[1L]], as.name("(")) ||
      identical(expr[[1L]], as.name("+")))){
    return(.formula_strip_no_intercept(expr[[2L]]))
  }

  if(.formula_is_no_intercept_additive_term(expr)){
    return(NULL)
  }

  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    lhs <- .formula_strip_no_intercept(expr[[2L]])
    rhs <- .formula_strip_no_intercept(expr[[3L]])

    if(is.null(lhs)){
      return(rhs)
    }
    if(is.null(rhs)){
      return(lhs)
    }
    return(call("+", lhs, rhs))
  }

  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 3L &&
     .formula_is_numeric_constant(expr[[3L]], 1)){
    return(.formula_strip_no_intercept(expr[[2L]]))
  }

  return(expr)
}

.formula_is_no_intercept_additive_term <- function(expr){

  .formula_is_numeric_constant(expr, 0) || .formula_is_negative_one(expr)
}

.formula_is_negative_one <- function(expr){

  (is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), -1)) ||
    (is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L &&
       .formula_is_numeric_constant(expr[[2L]], 1))
}

.formula_is_numeric_constant <- function(expr, value){

  if(is.call(expr) && length(expr) == 2L &&
     (identical(expr[[1L]], as.name("(")) ||
      identical(expr[[1L]], as.name("+")))){
    return(.formula_is_numeric_constant(expr[[2L]], value))
  }

  is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), as.numeric(value))
}
