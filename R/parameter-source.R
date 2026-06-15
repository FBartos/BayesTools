#' External JAGS parameter source
#'
#' @description
#' `parameter_source()` describes an existing JAGS node that another generated
#' model component may reference without creating or owning a prior for it.
#'
#' @param name JAGS node name.
#' @param shape source shape. `"scalar"` references one existing JAGS node;
#'   `"row"` references a vector node aligned to the formula rows. Generated
#'   row-shaped random-effect syntax indexes the node inside the observation
#'   loop.
#' @param values optional function for R-side reconstruction of row-shaped
#'   sources. When supplied with `shape = "row"`, it must accept named
#'   arguments `parameters`, `data`, and `n_rows`, or include `...`, and return
#'   a numeric vector of length `n_rows`. This is used by prediction and bridge-sampling
#'   reconstruction when the row-shaped source is not present as explicit
#'   posterior columns. For formula bridge reconstruction, `data` contains raw
#'   row-aligned formula/model data rather than standardized design matrices.
#'   Scalar sources do not accept `values`.
#'
#' @return A list-like `parameter_source` object.
#' @export
parameter_source <- function(name, shape = c("scalar", "row"),
                             values = NULL){

  check_char(name, "name", allow_NA = FALSE)
  .bt_check_external_parameter_source_name(name)
  .bt_check_jags_node_name(name, "name")
  shape <- match.arg(shape)
  if(!is.null(values) && !is.function(values)){
    stop("'values' must be NULL or a function.", call. = FALSE)
  }
  .bt_check_parameter_source_values_function(values, "values")
  if(identical(shape, "scalar") && !is.null(values)){
    stop("'values' is supported only for row-shaped parameter sources.", call. = FALSE)
  }

  out <- list(
    name = name,
    shape = shape,
    values = values
  )
  class(out) <- c("parameter_source", "list")

  out
}

.bt_check_parameter_source <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!inherits(x, "parameter_source")){
    stop("'source' must be created with parameter_source().", call. = FALSE)
  }
  check_char(x$name, "source$name", allow_NA = FALSE)
  .bt_check_external_parameter_source_name(x$name)
  .bt_check_jags_node_name(x$name, "source$name")
  if(!is.character(x$shape) || length(x$shape) != 1L ||
     is.na(x$shape) || !x$shape %in% c("scalar", "row")){
    stop("'source' metadata are inconsistent.", call. = FALSE)
  }
  if(!is.null(x$values) && !is.function(x$values)){
    stop("'source$values' must be NULL or a function.", call. = FALSE)
  }
  .bt_check_parameter_source_values_function(x$values, "source$values")
  if(identical(x$shape, "scalar") && !is.null(x$values)){
    stop("'source$values' is supported only for row-shaped parameter sources.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_parameter_source_values_function <- function(values, name){

  if(is.null(values)){
    return(invisible(TRUE))
  }
  if(!is.function(values)){
    stop("'", name, "' must be NULL or a function.", call. = FALSE)
  }
  value_formals <- formals(values)
  if(is.null(value_formals)){
    stop(
      "'", name,
      "' must be an R function that accepts named arguments 'parameters', 'data', and 'n_rows', or includes '...'.",
      call. = FALSE
    )
  }
  formal_names <- names(value_formals)
  if("..." %in% formal_names){
    return(invisible(TRUE))
  }
  required_names <- c("parameters", "data", "n_rows")
  missing_names <- setdiff(required_names, formal_names)
  if(length(missing_names) > 0L){
    stop(
      "'", name,
      "' must accept named arguments 'parameters', 'data', and 'n_rows', or include '...'.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_parameter_source_shape <- function(source){

  if(!is.character(source$shape) || length(source$shape) != 1L ||
     is.na(source$shape) || !source$shape %in% c("scalar", "row")){
    stop("'source' metadata are inconsistent.", call. = FALSE)
  }

  source$shape
}

.bt_parameter_source_is_row <- function(source){

  identical(.bt_parameter_source_shape(source), "row")
}

.bt_parameter_source_jags_expression <- function(source, row_index = NULL){

  .bt_check_parameter_source(source)
  if(!.bt_parameter_source_is_row(source)){
    return(source$name)
  }
  if(is.null(row_index)){
    stop(
      "Row-shaped parameter source '", source$name,
      "' requires a row index in generated JAGS syntax.",
      call. = FALSE
    )
  }
  check_char(row_index, "row_index", allow_NA = FALSE)
  .bt_check_jags_simple_index(row_index, "row_index")

  paste0(source$name, "[", row_index, "]")
}

.bt_parameter_source_label <- function(source){

  if(is.null(source)){
    return(NA_character_)
  }
  if(.bt_parameter_source_is_row(source)){
    return(paste0(source$name, "[row]"))
  }

  source$name
}

.bt_parameter_source_row_names <- function(source, n_rows){

  .bt_check_parameter_source(source)
  if(!.bt_parameter_source_is_row(source)){
    stop(
      "Parameter source '", source$name, "' is not row-shaped.",
      call. = FALSE
    )
  }
  check_int(n_rows, "n_rows", lower = 1, allow_NA = FALSE)

  paste0(source$name, "[", seq_len(n_rows), "]")
}

.bt_check_jags_node_name <- function(x, name){

  if(!grepl("^[A-Za-z][A-Za-z0-9_]*$", x)){
    stop(
      "'", name, "' must start with a letter and contain only letters, numbers, and underscores.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_check_external_parameter_source_name <- function(x){

  .bt_validate_random_effect_reserved_name(
    x,
    context = "external parameter sources"
  )

  invisible(TRUE)
}

.bt_check_jags_simple_index <- function(x, name){

  if(!grepl("^([A-Za-z][A-Za-z0-9_]*|[1-9][0-9]*)$", x)){
    stop(
      "'", name, "' must be a simple JAGS index name or positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @rdname parameter_source
#' @description
#' `random_sd_source()` marks a parameter source as an external random-effect
#' SD source. It is intended for variance allocations that should use an
#' already-defined JAGS node, such as scalar `tau` or row-shaped `tau`, instead
#' of creating a prior-owned total SD parameter.
#'
#' @param source JAGS node name or `parameter_source()` object. When `source`
#'   is a `parameter_source()` object, `shape` must not be supplied separately;
#'   the shape is inherited from `source`.
#'
#' @return A list-like `random_sd_source` object.
#' @export
random_sd_source <- function(source, shape = c("scalar", "row")){

  if(inherits(source, "parameter_source")){
    if(!missing(shape)){
      stop("'shape' must not be supplied when 'source' is already a parameter_source().", call. = FALSE)
    }
    parameter <- source
  }else{
    shape <- match.arg(shape)
    parameter <- parameter_source(source, shape = shape)
  }
  .bt_check_parameter_source(parameter)

  out <- list(
    name = parameter$name,
    shape = parameter$shape,
    source = parameter,
    kind = "external",
    owned = FALSE
  )
  class(out) <- c("random_sd_source", "list")

  out
}

.bt_check_random_sd_source <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!inherits(x, "random_sd_source")){
    stop("'sd_source' must be created with random_sd_source().", call. = FALSE)
  }
  .bt_check_parameter_source(x$source)
  if(!identical(x$kind, "external") || !identical(x$owned, FALSE)){
    stop("'sd_source' metadata are inconsistent.", call. = FALSE)
  }
  if(!identical(x$name, x$source$name) ||
     !identical(x$shape, x$source$shape)){
    stop("'sd_source' metadata are inconsistent.", call. = FALSE)
  }
  if(!is.null(x$source$values) && !is.function(x$source$values)){
    stop("'sd_source$source$values' must be NULL or a function.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_random_sd_source_expression <- function(x, row_index = NULL){

  .bt_check_random_sd_source(x)
  .bt_parameter_source_jags_expression(x$source, row_index = row_index)
}

.bt_random_sd_source_label <- function(x){

  .bt_check_random_sd_source(x)
  .bt_parameter_source_label(x$source)
}

.bt_random_sd_source_is_row_indexed <- function(x){

  .bt_check_random_sd_source(x)
  .bt_parameter_source_is_row(x$source)
}

.bt_parameter_source_values_function <- function(source){

  if(inherits(source, "random_sd_source")){
    .bt_check_random_sd_source(source)
    return(source$source$values)
  }
  .bt_check_parameter_source(source)

  source$values
}

.bt_parameter_source_has_values <- function(source){

  !is.null(.bt_parameter_source_values_function(source))
}

.bt_parameter_source_draw_parameters <- function(posterior, draw,
                                                parameters = NULL){

  out <- as.list(stats::setNames(
    as.numeric(posterior[draw, , drop = TRUE]),
    colnames(posterior)
  ))
  if(!is.null(parameters)){
    if(!is.list(parameters)){
      stop("'parameters' must be a list.", call. = FALSE)
    }
    out[names(parameters)] <- parameters
  }

  out
}

.bt_parameter_source_value_draws <- function(source, n_rows, posterior,
                                            data = NULL,
                                            parameters = NULL,
                                            context = "Parameter source"){

  .bt_check_parameter_source(source)
  values_function <- .bt_parameter_source_values_function(source)
  if(is.null(values_function)){
    return(NULL)
  }
  if(!is.matrix(posterior)){
    stop("'posterior' must be a matrix.", call. = FALSE)
  }
  check_int(n_rows, "n_rows", lower = 1, allow_NA = FALSE)

  out <- matrix(NA_real_, nrow = nrow(posterior), ncol = n_rows)
  for(draw in seq_len(nrow(posterior))){
    draw_parameters <- .bt_parameter_source_draw_parameters(
      posterior = posterior,
      draw = draw,
      parameters = parameters
    )
    values <- tryCatch(
      values_function(
        parameters = draw_parameters,
        data = data,
        n_rows = n_rows
      ),
      error = function(e)e
    )
    if(inherits(values, "error")){
      stop(
        context, " for source '", .bt_parameter_source_label(source),
        "' failed: ", conditionMessage(values),
        call. = FALSE
      )
    }
    if(!is.numeric(values)){
      stop(
        context, " for source '", .bt_parameter_source_label(source),
        "' must return a numeric vector.",
        call. = FALSE
      )
    }
    values <- as.numeric(values)
    if(length(values) != n_rows){
      stop(
        context, " for source '", .bt_parameter_source_label(source),
        "' must return a numeric vector of length ", n_rows, ".",
        call. = FALSE
      )
    }
    if(anyNA(values)){
      stop(
        context, " for source '", .bt_parameter_source_label(source),
        "' returned missing values.",
        call. = FALSE
      )
    }
    out[draw, ] <- values
  }

  colnames(out) <- .bt_parameter_source_row_names(source, n_rows)
  out
}
