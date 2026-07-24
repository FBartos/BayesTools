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
    model_matrix = model_matrix
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
  data[[index_name]] <- .bt_random_effect_structured_index_values(data, index_variables)
  out$data <- data
  out$formula <- stats::as.formula(
    call("~", call("-", as.name(index_name), 1)),
    env = environment(formula)
  )
  out$index <- list(
    variables = index_variables,
    name = index_name,
    label = paste(index_variables, collapse = ":"),
    structure = structure
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

.bt_random_effect_structured_index_values <- function(data, variables){

  values <- lapply(variables, function(variable){
    .bt_random_effect_structured_index_component(data[[variable]])
  })

  if(length(values) == 1L){
    return(values[[1L]])
  }

  do.call(interaction, c(values, list(drop = TRUE, lex.order = TRUE)))
}

.bt_random_effect_structured_index_component <- function(x){

  if(is.factor(x)){
    return(x)
  }
  if(is.character(x)){
    return(factor(x, levels = sort(unique(x))))
  }
  if(is.numeric(x) || is.integer(x) || is.logical(x)){
    return(factor(x, levels = sort(unique(x))))
  }

  factor(x, levels = sort(unique(x)))
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
