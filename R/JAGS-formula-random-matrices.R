.bt_random_effect_design_matrix <- function(formula, data,
                                            preserve_no_intercept_contrasts = TRUE,
                                            structure = NULL,
                                            car_time_values = NULL,
                                            block_name = NULL){

  if(identical(structure, "car")){
    return(.bt_random_effect_car_design_matrix(
      formula = formula,
      data = data,
      car_time_values = car_time_values
    ))
  }

  formula_terms <- stats::terms(formula)
  has_intercept <- attr(formula_terms, "intercept") == 1L
  matrix_formula <- formula
  if(!has_intercept && isTRUE(preserve_no_intercept_contrasts)){
    matrix_formula <- formula_add_intercept(formula)
  }

  label <- "Random-effect term"
  if(!is.null(block_name) && nzchar(block_name)){
    label <- paste0("Random-effect block '", block_name, "'")
  }
  model_frame <- stats::model.frame(matrix_formula, data = data)
  if(nrow(model_frame) != nrow(data) || anyNA(model_frame)){
    stop(
      label,
      " contains missing predictor values; random-effect design matrices must have one row per data row.",
      call. = FALSE
    )
  }
  model_matrix <- .bt_model_matrix(model_frame, formula = matrix_formula, data = data)
  .bt_validate_model_matrix_finite(model_matrix, label)

  if(!has_intercept && isTRUE(preserve_no_intercept_contrasts)){
    intercept_column <- which(colnames(model_matrix) == "(Intercept)")
    if(length(intercept_column) == 1L){
      assign    <- attr(model_matrix, "assign")
      contrasts <- attr(model_matrix, "contrasts", exact = TRUE)
      model_matrix <- model_matrix[, -intercept_column, drop = FALSE]
      attr(model_matrix, "assign") <- assign[-intercept_column]
      attr(model_matrix, "contrasts") <- contrasts
    }
  }

  list(
    model_frame = model_frame,
    model_matrix = model_matrix,
    exact_indicator = structure %in% c("cs", "hcs", "ar1", "har")
  )
}

.bt_random_effect_normalize_structured_formula <- function(formula, data,
                                                           structure){

  out <- list(formula = formula, data = data, index = NULL)
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har")){
    return(out)
  }

  if(identical(structure, "car")){
    out$formula <- .bt_random_effect_formula_without_intercept(formula)
    return(out)
  }

  index_variables <- .bt_random_effect_structured_index_variables(formula, structure)
  if(structure %in% c("ar1", "har") && length(index_variables) > 1L){
    stop(
      "The '", structure,
      "' random-effect covariance structure requires a single ordered index variable. ",
      "Create an explicit ordered index column before calling '", structure, "()'.",
      call. = FALSE
    )
  }

  missing_variables <- index_variables[!index_variables %in% names(data)]
  if(length(missing_variables) > 0L){
    stop(
      paste0(
        "The ",
        paste0("'", missing_variables, "'", collapse = ", "),
        " structured random-effect index variable is missing in the data set."
      ),
      call. = FALSE
    )
  }

  index_name <- .bt_random_effect_structured_index_name(index_variables)
  index_levels <- .bt_random_effect_structured_index_resolve(data, index_variables)
  data[[index_name]] <- .bt_random_effect_structured_index_factor(index_levels)
  if(structure %in% c("ar1", "har")){
    .bt_random_effect_warn_index_order(
      x = out$data[[index_variables]],
      resolved_levels = levels(data[[index_name]]),
      variable = index_variables,
      structure = structure
    )
  }
  out$data <- data
  out$formula <- stats::as.formula(
    call("~", call("-", as.name(index_name), 1)),
    env = environment(formula)
  )
  out$index <- list(
    variables = index_variables,
    name = index_name,
    label = paste(index_variables, collapse = ":"),
    structure = structure,
    levels = index_levels$labels,
    level_keys = index_levels$level_keys,
    component_levels = index_levels$components
  )

  out
}

.bt_random_effect_formula_without_intercept <- function(formula){

  rhs_index <- if(length(formula) == 3L) 3L else 2L
  stats::as.formula(
    call("~", call("-", formula[[rhs_index]], 1)),
    env = environment(formula)
  )
}

.bt_random_effect_structured_index_variables <- function(formula, structure){

  rhs_index <- if(length(formula) == 3L) 3L else 2L
  rhs <- formula[[rhs_index]]
  usage <- if(structure %in% c("cs", "hcs")){
    paste0("Use '", structure, "(index | group)' or '", structure, "(index1 + index2 | group)'.")
  }else{
    paste0("Use '", structure, "(index | group)'.")
  }
  if(.bt_random_effect_structured_index_has_intercept_control(rhs)){
    stop(
      "The '", structure,
      "' random-effect covariance structure uses index variables and does not support explicit ",
      "'1', '0', or '-1' terms. ", usage,
      call. = FALSE
    )
  }

  variables <- .bt_random_effect_structured_index_plus_terms(rhs, structure)
  if(length(variables) == 0L || any(!nzchar(variables))){
    stop(
      "The '", structure,
      "' random-effect covariance structure requires at least one index variable.",
      call. = FALSE
    )
  }
  if(anyDuplicated(variables)){
    stop(
      "Structured random-effect index variables must be unique.",
      call. = FALSE
    )
  }

  variables
}

.bt_random_effect_structured_index_plus_terms <- function(expr, structure){

  if(is.symbol(expr)){
    return(as.character(expr))
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    return(c(
      .bt_random_effect_structured_index_plus_terms(expr[[2L]], structure),
      .bt_random_effect_structured_index_plus_terms(expr[[3L]], structure)
    ))
  }

  stop(
    "The '", structure,
    "' random-effect covariance structure supports index variables",
    if(structure %in% c("cs", "hcs")) " separated by '+'." else ".",
    " ",
    if(structure %in% c("cs", "hcs")){
      paste0("Use '", structure, "(index | group)' or '", structure, "(index1 + index2 | group)'.")
    }else{
      paste0("Use '", structure, "(index | group)'.")
    },
    call. = FALSE
  )
}

.bt_random_effect_structured_index_has_intercept_control <- function(expr){

  if(is.numeric(expr) && length(expr) == 1L && expr %in% c(0, 1)){
    return(TRUE)
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L){
    return(.bt_random_effect_structured_index_has_intercept_control(expr[[2L]]))
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 3L){
    return(
      .bt_random_effect_structured_index_has_intercept_control(expr[[2L]]) ||
        .bt_random_effect_structured_index_has_intercept_control(expr[[3L]])
    )
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    return(
      .bt_random_effect_structured_index_has_intercept_control(expr[[2L]]) ||
        .bt_random_effect_structured_index_has_intercept_control(expr[[3L]])
    )
  }

  FALSE
}

.bt_random_effect_structured_index_name <- function(variables){

  if(length(variables) == 1L){
    return(variables)
  }

  .bt_random_effect_sanitize_name(paste(variables, collapse = "_"))
}

# Structured index levels are identified by exact keys: factor, character, and
# logical values by their labels, numeric values by their exact double value,
# and multi-variable cells by the length-prefixed tuple of component keys.
# Display labels follow the factor()/interaction(lex.order = TRUE) labels and
# are changed only where distinct keys would otherwise share a label.
.bt_random_effect_structured_index_values <- function(data, variables,
                                                      index = NULL){

  if(is.null(index) || is.null(index$level_keys) || is.null(index$levels) ||
     is.null(index$component_levels)){
    # Fits without stored index keys replay the index by its display labels.
    return(.bt_random_effect_structured_index_factor(
      .bt_random_effect_structured_index_resolve(data, variables)
    ))
  }

  level_keys <- index$level_keys
  levels <- index$levels
  component_levels <- index$component_levels
  if(!is.character(level_keys) || !is.character(levels) ||
     length(level_keys) != length(levels) || anyNA(level_keys) ||
     anyNA(levels) || anyDuplicated(level_keys) || anyDuplicated(levels) ||
     !is.list(component_levels) ||
     length(component_levels) != length(variables)){
    stop(
      "Structured random-effect index metadata for '", index$name,
      "' are malformed. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  row_keys <- .bt_random_effect_structured_index_prediction_keys(
    data = data,
    variables = variables,
    component_levels = component_levels
  )
  level_index <- match(row_keys, level_keys)
  if(anyNA(level_index[!is.na(row_keys)])){
    stop(
      "Levels specified in the '", index$name,
      "' factor variable do not match the levels used for model specification.",
      call. = FALSE
    )
  }
  factor(levels[level_index], levels = levels)
}

# Prediction keys of index values: the exact key when it is a fitted level;
# otherwise the fitted level whose display label (or, for a numeric fitted
# level, its 15-significant-digit label) equals the value's label, when that
# label identifies exactly one fitted level. Anything else keeps its own key
# and is a new index level.
.bt_random_effect_structured_index_prediction_keys <- function(
    data, variables, component_levels){

  component_keys <- lapply(seq_along(variables), function(i){
    x <- data[[variables[[i]]]]
    fitted <- component_levels[[i]]
    keys <- .bt_random_effect_structured_index_component(
      x,
      variables[[i]]
    )$row_keys
    unmatched <- which(!is.na(keys) & !keys %in% fitted$level_keys)
    if(length(unmatched) == 0L){
      return(keys)
    }
    fitted_short_labels <- if(isTRUE(fitted$numeric)){
      as.character(as.numeric(fitted$level_keys))
    }else{
      fitted$labels
    }
    value_labels <- as.character(x)[unmatched]
    for(j in seq_along(unmatched)){
      candidates <- which(
        fitted$labels == value_labels[[j]] |
          fitted_short_labels == value_labels[[j]]
      )
      if(length(candidates) == 1L){
        keys[[unmatched[[j]]]] <- fitted$level_keys[[candidates]]
      }
    }
    keys
  })
  if(length(component_keys) == 1L){
    return(component_keys[[1L]])
  }

  key_matrix <- do.call(cbind, component_keys)
  missing_rows <- apply(is.na(key_matrix), 1L, any)
  row_keys <- rep(NA_character_, nrow(key_matrix))
  row_keys[!missing_rows] <- apply(
    key_matrix[!missing_rows, , drop = FALSE],
    1L,
    .bt_random_group_tuple_key
  )
  row_keys
}

.bt_random_effect_structured_index_factor <- function(resolved){

  factor(
    resolved$labels[match(resolved$row_keys, resolved$level_keys)],
    levels = resolved$labels
  )
}

.bt_random_effect_structured_index_resolve <- function(data, variables){

  components <- lapply(variables, function(variable){
    .bt_random_effect_structured_index_component(data[[variable]], variable)
  })
  component_levels <- lapply(components, `[`,
                             c("level_keys", "labels", "numeric"))
  if(length(components) == 1L){
    return(c(
      components[[1L]][c("row_keys", "level_keys", "labels")],
      list(components = component_levels)
    ))
  }

  row_key_matrix <- do.call(cbind, lapply(components, `[[`, "row_keys"))
  codes <- do.call(cbind, lapply(components, function(component){
    match(component$row_keys, component$level_keys)
  }))
  missing_rows <- apply(is.na(row_key_matrix), 1L, any)
  row_keys <- rep(NA_character_, nrow(row_key_matrix))
  row_keys[!missing_rows] <- apply(
    row_key_matrix[!missing_rows, , drop = FALSE],
    1L,
    .bt_random_group_tuple_key
  )

  # Observed cells only, ordered with the first index variable varying
  # slowest, as interaction(drop = TRUE, lex.order = TRUE).
  first_rows <- which(!duplicated(row_keys) & !missing_rows)
  first_rows <- first_rows[do.call(
    order,
    unname(as.data.frame(codes[first_rows, , drop = FALSE]))
  )]
  label_parts <- do.call(cbind, lapply(seq_along(components), function(i){
    components[[i]]$labels[codes[first_rows, i]]
  }))
  labels <- .bt_random_effect_structured_index_tuple_labels(
    label_parts = label_parts,
    component_labels = lapply(components, `[[`, "labels"),
    variables = variables
  )

  list(
    row_keys = row_keys,
    level_keys = row_keys[first_rows],
    labels = labels,
    components = component_levels
  )
}

.bt_random_effect_structured_index_tuple_labels <- function(label_parts,
                                                            component_labels,
                                                            variables){

  if(nrow(label_parts) == 0L){
    return(character())
  }
  labels <- apply(label_parts, 1L, paste, collapse = ".")
  if(!anyDuplicated(labels)){
    return(labels)
  }

  # A separator absent from every component label keeps the cell labels
  # injective. It excludes characters reserved by coefficient names and
  # semantic parameter labels.
  used_labels <- unlist(component_labels, use.names = FALSE)
  for(separator in c("_", "-", "~", "/", "&", "#", "@")){
    if(!any(grepl(separator, used_labels, fixed = TRUE))){
      return(apply(label_parts, 1L, paste, collapse = separator))
    }
  }

  stop(
    "The structured random-effect index levels of ",
    paste0("'", variables, "'", collapse = ", "),
    " cannot be represented by unique labels. Recode the index levels.",
    call. = FALSE
  )
}

.bt_random_effect_structured_index_component <- function(x, variable = ""){

  if(is.factor(x)){
    labels <- levels(x)
    return(list(
      row_keys = as.character(x),
      level_keys = labels,
      labels = labels,
      numeric = FALSE
    ))
  }
  if(is.numeric(x)){
    x <- as.numeric(x)
    values <- sort(unique(x))
    # The exact double identifies the level; signed zeros are one value.
    values[values == 0] <- 0
    level_keys <- sprintf("%.17g", values)
    labels <- as.character(values)
    ambiguous <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    labels[ambiguous] <- level_keys[ambiguous]
    if(anyDuplicated(labels) || anyDuplicated(level_keys)){
      stop(
        "The structured random-effect index variable '", variable,
        "' cannot be represented by unique level labels.",
        call. = FALSE
      )
    }
    return(list(
      row_keys = level_keys[match(x, values)],
      level_keys = level_keys,
      labels = labels,
      numeric = TRUE
    ))
  }

  values <- sort(unique(x))
  labels <- as.character(values)
  row_keys <- as.character(x)
  row_keys[is.na(x)] <- NA_character_
  list(
    row_keys = row_keys,
    level_keys = labels,
    labels = labels,
    numeric = FALSE
  )
}

.bt_random_effect_warn_index_order <- function(x, resolved_levels, variable,
                                                structure){

  if(is.ordered(x) || !(is.factor(x) || is.character(x))){
    return(invisible(NULL))
  }
  numeric_levels <- suppressWarnings(as.numeric(resolved_levels))
  if(any(!is.finite(numeric_levels)) || !is.unsorted(numeric_levels)){
    return(invisible(NULL))
  }

  warning(
    "The '", structure, "' index variable '", variable, "' uses ",
    if(is.factor(x)) "declared factor level" else "sorted character",
    " order: ", paste(resolved_levels, collapse = ", "),
    ". Its numeric labels are not in increasing numeric order. ",
    "Use numeric values for numeric ordering, or an ordered factor with ",
    "explicit levels to confirm the intended order. The existing order is preserved.",
    call. = FALSE
  )
  invisible(NULL)
}

.bt_random_effect_car_design_matrix <- function(formula, data,
                                                car_time_values = NULL){

  formula_terms <- stats::terms(formula)
  has_intercept <- attr(formula_terms, "intercept") == 1L
  if(has_intercept){
    stop(
      "CAR random-effect terms require exactly one untransformed time variable. Use 'car(time | group)'.",
      call. = FALSE
    )
  }

  term_labels <- attr(formula_terms, "term.labels")
  predictors <- as.character(attr(formula_terms, "variables"))[-1L]
  if(length(term_labels) != 1L || length(predictors) != 1L ||
     !identical(term_labels, predictors)){
    stop(
      "CAR random-effect terms require exactly one untransformed time variable. Use 'car(time | group)'.",
      call. = FALSE
    )
  }

  time_name <- predictors
  if(!time_name %in% names(data)){
    stop("The '", time_name, "' CAR time variable is missing in the data set.", call. = FALSE)
  }

  observed_time <- .bt_random_effect_car_observed_time(data[[time_name]], time_name)
  if(is.null(car_time_values)){
    car_time_values <- sort(unique(observed_time))
  }else{
    car_time_values <- .bt_random_effect_car_reference_time_values(car_time_values)
    new_time <- observed_time[!observed_time %in% car_time_values]
    if(length(new_time) > 0L){
      stop(
        "New CAR time coordinate(s) for random-effect prediction are not supported: ",
        paste(sort(unique(new_time)), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  time_index <- match(observed_time, car_time_values)
  model_matrix <- matrix(0, nrow = length(observed_time), ncol = length(car_time_values))
  if(length(observed_time) > 0L){
    model_matrix[cbind(seq_along(observed_time), time_index)] <- 1
  }
  colnames(model_matrix) <- paste0(time_name, .bt_random_effect_car_time_suffix(car_time_values))
  attr(model_matrix, "assign") <- rep(1L, ncol(model_matrix))
  attr(model_matrix, "contrasts") <- NULL

  model_frame <- data[time_name]
  attr(model_frame, "terms") <- formula_terms

  list(
    model_frame = model_frame,
    model_matrix = model_matrix,
    exact_indicator = TRUE,
    car = list(
      time_variable = time_name,
      time_values = car_time_values
    )
  )
}

.bt_random_effect_car_observed_time <- function(x, time_name){

  if(is.numeric(x) || is.integer(x)){
    out <- as.numeric(x)
  }else if(is.ordered(x)){
    level_values <- suppressWarnings(as.numeric(as.character(levels(x))))
    if(any(is.na(level_values))){
      stop(
        "Ordered factor CAR time variable '", time_name,
        "' must have numeric level labels.",
        call. = FALSE
      )
    }
    out <- level_values[match(as.character(x), levels(x))]
  }else{
    stop(
      "CAR time variable '", time_name,
      "' must be numeric or an ordered factor with numeric level labels.",
      call. = FALSE
    )
  }

  if(any(is.na(out)) || any(!is.finite(out))){
    stop("CAR time variable '", time_name, "' must contain only finite values.", call. = FALSE)
  }

  out
}

.bt_random_effect_car_reference_time_values <- function(x){

  if(!is.numeric(x) && !is.integer(x)){
    stop("CAR reference time values must be numeric.", call. = FALSE)
  }
  x <- as.numeric(x)
  if(length(x) == 0L || any(is.na(x)) || any(!is.finite(x))){
    stop("CAR reference time values must be finite and non-empty.", call. = FALSE)
  }
  if(anyDuplicated(x)){
    stop("CAR reference time values must be unique.", call. = FALSE)
  }
  sort(x)
}

.bt_random_effect_car_time_suffix <- function(x){

  labels <- vapply(x, function(value){
    format(value, scientific = FALSE, trim = TRUE, digits = 17)
  }, character(1))
  labels <- gsub("-", "m", labels, fixed = TRUE)
  labels <- gsub(".", "p", labels, fixed = TRUE)
  labels <- gsub("[^A-Za-z0-9_]", "_", labels)
  labels <- gsub("_+", "_", labels)
  labels <- gsub("^_|_$", "", labels)
  if(anyDuplicated(labels)){
    stop(
      "CAR time values cannot be represented by unique random-effect labels.",
      call. = FALSE
    )
  }
  paste0("_", labels)
}
