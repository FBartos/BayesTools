#' @title JAGS parameter-name helpers
#'
#' @description Format formula-derived parameter names and extract exact indexed
#' JAGS parameter columns.
#'
#' @param parameters a vector of parameter names
#' @param formula_parameter a formula parameter prefix name
#' @param formula_parameters a vector of formula parameter prefix names
#' @param formula_random a vector of random effects grouping factors
#' @param formula_prefix whether the \code{formula_parameters} names should be
#' kept. Defaults to \code{TRUE}.
#' @param formula_scale optional nested list containing scaling info. When provided,
#' intercepts from parameters with \code{log_intercept = TRUE} attribute will be
#' renamed to \code{exp(intercept)}.
#' @param x character vector to escape for regular-expression matching.
#' @param columns character vector of column names.
#' @param parameter exact unindexed JAGS parameter name.
#' @param samples matrix, data frame, \code{mcmc}, or \code{mcmc.list} object.
#' @param drop_missing whether missing indexed columns should return
#' \code{NULL}. If \code{FALSE}, a zero-column matrix with the input row count
#' is returned.
#' @param row named row vector, list, or one-row data frame.
#'
#' @examples
#' format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"),
#'                        formula_parameters = "mu")
#'
#' @return \code{format_parameter_names()}, \code{JAGS_parameter_names()}, and
#' \code{JAGS_regex_escape()} return character vectors.
#' \code{JAGS_indexed_parameter_columns()} returns a logical vector.
#' \code{JAGS_indexed_parameter_matrix()} returns a matrix of exact indexed
#' parameter columns sorted by numeric index, or \code{NULL} when none are
#' present and \code{drop_missing = TRUE}.
#' \code{JAGS_indexed_parameter_vector()} returns indexed row values sorted by
#' numeric index.
#'
#' @export format_parameter_names
#' @export JAGS_parameter_names
#' @export JAGS_regex_escape
#' @export JAGS_indexed_parameter_columns
#' @export JAGS_indexed_parameter_matrix
#' @export JAGS_indexed_parameter_vector
#' @name parameter_names
NULL

#' @rdname parameter_names
format_parameter_names <- function(parameters, formula_parameters = NULL, formula_random = NULL, formula_prefix = TRUE, formula_scale = NULL){

  check_char(parameters, "parameters", check_length = FALSE)
  check_char(formula_random, "formula_random", check_length = FALSE, allow_NULL = TRUE)
  check_char(formula_parameters, "formula_parameters", check_length = FALSE, allow_NULL = TRUE)
  check_bool(formula_prefix, "formula_prefix")
  check_list(formula_scale, "formula_scale", allow_NULL = TRUE)

  # rename intercept to exp(intercept) for parameters with log_intercept attribute
  if(!is.null(formula_scale)){
    for(param_name in names(formula_scale)){
      if(isTRUE(attr(formula_scale[[param_name]], "log_intercept"))){
        intercept_name <- paste0(param_name, "_intercept")
        if(intercept_name %in% parameters){
          parameters[parameters == intercept_name] <- paste0(param_name, "_exp(intercept)")
        }
      }
    }
  }

  for(i in seq_along(formula_parameters)){
    formula_prefix_pattern <- paste0(formula_parameters[i], "_")
    matching_parameters <- grepl(
      formula_prefix_pattern,
      parameters,
      fixed = TRUE
    )
    parameters[matching_parameters] <- gsub(
      formula_prefix_pattern,
      if(formula_prefix) paste0("(", formula_parameters[i], ") ") else "",
      parameters[matching_parameters],
      fixed = TRUE
    )
  }

  for(i in seq_along(formula_random)){
    temp_which <- grepl(
      paste0("_xREx__", formula_random[i], "_"),
      parameters,
      fixed = TRUE
    )
    temp_incl  <- grepl("(inclusion)", parameters)
    parameters[temp_which] <- gsub(
      paste0("_xREx__", formula_random[i], "_"),
      "",
      parameters[temp_which],
      fixed = TRUE
    )
    if(any(temp_which &  temp_incl)){
      parameters[temp_which &  temp_incl] <- paste0(gsub("(inclusion)", "", parameters[temp_which & temp_incl], fixed = TRUE), "|", formula_random[i], " (inclusion)")
    }
    if(any(temp_which & !temp_incl)){
      parameters[temp_which & !temp_incl] <- paste0("sd(", parameters[temp_which & !temp_incl], "|", formula_random[i], ")")
    }
  }

  parameters[grep("__xXx__", parameters)] <- gsub("__xXx__", ":", parameters[grep("__xXx__", parameters)])

  return(parameters)
}
#' @rdname parameter_names
JAGS_parameter_names   <- function(parameters, formula_parameter = NULL){

  check_char(parameters, "parameters", check_length = FALSE)
  check_char(formula_parameter, "formula_parameter", check_length = TRUE, allow_NULL = TRUE)

  if(!is.null(formula_parameter)){
    parameters <- paste0(formula_parameter, "_", parameters)
  }
  parameters <- gsub(":", "__xXx__", parameters)

  return(parameters)
}

#' @rdname parameter_names
JAGS_regex_escape <- function(x){

  if(!is.character(x) || !is.vector(x)){
    stop("The 'x' argument must be a character vector.", call. = FALSE)
  }
  if(anyNA(x)){
    stop("The 'x' argument cannot contain NA/NaN values.", call. = FALSE)
  }
  if(length(x) == 0L){
    return(character())
  }

  special <- c("\\", ".", "|", "(", ")", "[", "]", "{", "}", "^", "$",
               "*", "+", "?")
  out <- vapply(x, function(x_i){
    chars <- strsplit(x_i, "", fixed = TRUE)[[1L]]
    chars <- vapply(chars, function(ch){
      if(ch %in% special){
        return(paste0("\\", ch))
      }
      ch
    }, character(1))
    paste0(chars, collapse = "")
  }, character(1))

  return(unname(out))
}

#' @rdname parameter_names
JAGS_indexed_parameter_columns <- function(columns, parameter){

  .JAGS_check_columns(columns)
  .JAGS_check_indexed_parameter(parameter)

  return(!is.na(.JAGS_indexed_parameter_indices(columns, parameter)))
}

#' @rdname parameter_names
JAGS_indexed_parameter_matrix <- function(samples, parameter,
                                         drop_missing = TRUE){

  check_bool(drop_missing, "drop_missing", allow_NA = FALSE)
  .JAGS_check_indexed_parameter(parameter)

  if(inherits(samples, "mcmc.list") || inherits(samples, "mcmc")){
    samples <- as.matrix(samples)
  }
  if(!is.matrix(samples) && !is.data.frame(samples)){
    stop("'samples' must be a matrix, data frame, mcmc, or mcmc.list object.",
         call. = FALSE)
  }
  samples_is_data_frame <- is.data.frame(samples)

  columns <- colnames(samples)
  if(is.null(columns)){
    if(drop_missing){
      return(NULL)
    }
    out <- samples[, integer(0), drop = FALSE]
    return(if(samples_is_data_frame) as.matrix(out) else out)
  }

  selected <- .JAGS_indexed_parameter_sorted_columns(columns, parameter)
  if(length(selected) == 0L){
    if(drop_missing){
      return(NULL)
    }
    out <- samples[, integer(0), drop = FALSE]
    return(if(samples_is_data_frame) as.matrix(out) else out)
  }

  out <- samples[, selected, drop = FALSE]
  return(if(samples_is_data_frame) as.matrix(out) else out)
}

#' @rdname parameter_names
JAGS_indexed_parameter_vector <- function(row, parameter){

  .JAGS_check_indexed_parameter(parameter)

  if(is.data.frame(row)){
    if(nrow(row) != 1L){
      stop("'row' data frames must contain exactly one row.", call. = FALSE)
    }
    row_names <- names(row)
  }else{
    row_names <- names(row)
  }
  if(is.null(row_names)){
    stop("'row' must be a named vector or one-row data frame.", call. = FALSE)
  }

  selected <- .JAGS_indexed_parameter_sorted_columns(row_names, parameter)
  if(is.data.frame(row)){
    return(unlist(row[1L, selected, drop = FALSE], use.names = TRUE))
  }
  if(is.list(row) && !is.atomic(row)){
    return(unlist(row[selected], use.names = TRUE))
  }
  if(length(selected) == 0L){
    return(row[integer(0)])
  }

  return(row[selected])
}

.JAGS_check_columns <- function(columns){

  if(!is.character(columns) || !is.vector(columns)){
    stop("The 'columns' argument must be a character vector.",
         call. = FALSE)
  }
  if(anyNA(columns)){
    stop("The 'columns' argument cannot contain NA/NaN values.",
         call. = FALSE)
  }

  return(invisible(TRUE))
}

.JAGS_check_indexed_parameter <- function(parameter){

  check_char(parameter, "parameter", check_length = 1, allow_NA = FALSE)
  if(!nzchar(parameter)){
    stop("The 'parameter' argument must not be empty.", call. = FALSE)
  }

  return(invisible(TRUE))
}

.JAGS_indexed_parameter_indices <- function(columns, parameter){

  pattern <- paste0("^", JAGS_regex_escape(parameter), "\\[([1-9][0-9]*)\\]$")
  matches <- regexec(pattern, columns, perl = TRUE)
  parts <- regmatches(columns, matches)
  out <- rep(NA_integer_, length(columns))
  has_match <- vapply(parts, length, integer(1)) == 2L
  index_text <- vapply(
    parts[has_match],
    `[[`,
    character(1),
    2L
  )
  index_values <- suppressWarnings(as.numeric(index_text))
  if(any(!is.finite(index_values)) ||
     any(index_values > .Machine$integer.max)){
    stop(
      "Indexed JAGS parameter columns contain an index outside the supported integer range.",
      call. = FALSE
    )
  }
  out[has_match] <- as.integer(index_values)

  return(out)
}

.JAGS_indexed_parameter_sorted_columns <- function(columns, parameter){

  indices <- .JAGS_indexed_parameter_indices(columns, parameter)
  keep <- which(!is.na(indices))
  if(length(keep) == 0L){
    return(integer())
  }
  duplicate_indices <- unique(indices[keep][duplicated(indices[keep])])
  if(length(duplicate_indices) > 0L){
    stop(
      "Indexed JAGS parameter '", parameter,
      "' contains duplicate index",
      if(length(duplicate_indices) > 1L) " values: " else ": ",
      paste(duplicate_indices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  return(keep[order(indices[keep])])
}

.JAGS_prior_factor_names <- function(parameter, prior){

  levels <- .get_prior_factor_levels(prior)
  if((is.null(levels) || length(levels) == 0L || is.na(levels)) && "K" %in% names(prior[["parameters"]])){
    levels <- prior[["parameters"]][["K"]]
  }
  if(is.null(levels) || length(levels) == 0L || is.na(levels)){
    stop("Factor-prior dimensions must be available before constructing JAGS parameter names.", call. = FALSE)
  }

  if(levels == 1){
    par_names <- parameter
  }else{
    par_names <- paste0(parameter, "[", 1:levels, "]")
  }

  return(par_names)
}
