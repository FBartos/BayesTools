#' @title Print a BayesTools table
#'
#' @param x a BayesTools_values_tables
#' @param ... additional arguments.
#'
#' @return \code{print.BayesTools_table} returns \code{NULL}.
#'
#' @exportS3Method
print.BayesTools_table <- function(x, ...){

  # print formatting
  for(i in seq_along(attr(x, "type"))){
    values         <- x[[i]]
    colnames(x)[i] <- .format_column_names(colnames(x)[i], attr(x, "type")[i], values)
    x[[i]]         <- .format_column(values, attr(x, "type")[i], attr(x, "n_models")[i])
  }

  # print title
  if(!is.null(attr(x, "title"))){
    cat(paste0(attr(x, "title"), "\n"))
  }

  # print the table
  print_rownames <- attr(x, "rownames")
  if(is.null(print_rownames)){
    print_rownames <- TRUE
  }
  print(
    as.data.frame(x),
    quote = FALSE,
    right = TRUE,
    row.names = print_rownames
  )

  # print footnotes
  for(i in seq_along(attr(x, "footnotes"))){
    cat(paste0(attr(x, "footnotes")[i], "\n"))
  }

  # print warnings in red
  for(i in seq_along(attr(x, "warnings"))){
    cat(paste0("\033[0;31m", attr(x, "warnings")[i], "\033[0m\n"))
  }

  return(invisible())
}

#' @title Format Bayes factor
#'
#' @description Formats Bayes factor
#'
#' @param BF Bayes factor(s)
#' @param logBF log(BF)
#' @param BF01 1/BF
#' @param inclusion whether the Bayes factor is an inclusion BF (for naming purposes)
#'
#' @return \code{format_BF} returns a formatted Bayes factor.
#'
#' @export
format_BF <- function(BF, logBF = FALSE, BF01 = FALSE, inclusion = FALSE){

  if(!is.numeric(BF)){
    check_real(BF, "BF", lower = 0, check_length = FALSE, allow_NA = TRUE)
  }
  BF_names <- names(BF)
  bound_operator <- .standardize_BF_bound_operator(attr(BF, "bound_operator"), length(BF))
  BF <- as.numeric(BF)
  names(BF) <- BF_names
  check_real(BF, "BF", lower = 0, check_length = FALSE, allow_NA = TRUE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01,  "BF01",  allow_NA = FALSE)

  if(BF01){
    BF   <- 1/BF
    bound_operator <- .invert_BF_bound_operator(bound_operator)
    name <- ifelse(inclusion, "Exclusion BF", "1/BF")
  }else{
    name <- ifelse(inclusion, "Inclusion BF", "BF")
  }

  if(logBF){
    BF   <- log(BF)
    name <- paste0("log(", name, ")")
  }

  attr(BF, "name")  <- name
  attr(BF, "logBF") <- logBF
  attr(BF, "BF01")  <- BF01
  attr(BF, "bound_operator") <- bound_operator
  if(any(!is.na(bound_operator))){
    class(BF) <- unique(c("BayesTools_BF", class(BF)))
  }

  return(BF)
}

#' @export
`[.BayesTools_BF` <- function(x, i, ...){

  bound_operator <- attr(x, "bound_operator")
  out <- NextMethod("[")

  attr(out, "name")  <- attr(x, "name")
  attr(out, "logBF") <- attr(x, "logBF")
  attr(out, "BF01")  <- attr(x, "BF01")
  if(!is.null(bound_operator)){
    attr(out, "bound_operator") <- if(missing(i)) bound_operator else bound_operator[i]
  }
  if(any(!is.na(attr(out, "bound_operator")))){
    class(out) <- unique(c("BayesTools_BF", class(out)))
  }

  return(out)
}

#' @export
`[.BayesTools_table` <- function(x, i, j, ..., drop = TRUE){

  original_names <- names(x)
  original_col_attributes <- lapply(x, attributes)

  out <- NextMethod("[")
  if(!is.data.frame(out)){
    return(out)
  }

  out_names <- names(out)
  source_cols <- .match_table_columns(out_names, original_names)
  valid_cols  <- !is.na(source_cols)

  if(length(source_cols) > 0 && any(valid_cols)){
    for(k in which(valid_cols)){
      out[[k]] <- .restore_table_column_attributes(out[[k]], original_col_attributes[[source_cols[k]]])
    }
  }

  type <- attr(x, "type")
  if(length(type) == length(original_names)){
    attr(out, "type") <- type[source_cols]
  }

  n_models <- attr(x, "n_models")
  if(length(n_models) == length(original_names)){
    attr(out, "n_models") <- n_models[source_cols]
  }

  for(attribute in c("title", "footnotes", "rownames")){
    attr(out, attribute) <- attr(x, attribute)
  }

  selected_parameters <- .subset_table_parameters(x, out)
  if(!is.null(selected_parameters)){
    attr(out, "parameters") <- selected_parameters
  }
  attr(out, "warnings") <- .subset_table_warnings(attr(x, "warnings"), selected_parameters, rownames(out))
  out <- .subset_table_hypothesis_attributes(x, out)

  out
}

.BF_column_name <- function(logBF = FALSE, BF01 = FALSE, inclusion = FALSE){

  if(BF01){
    name <- ifelse(inclusion, "Exclusion BF", "1/BF")
  }else{
    name <- ifelse(inclusion, "Inclusion BF", "BF")
  }

  if(logBF){
    name <- paste0("log(", name, ")")
  }

  return(name)
}

.runjags_indicator_list <- function(model_samples_list, indicator_name, indicator_values){

  lapply(model_samples_list, function(samples){
    samples <- as.matrix(samples)
    if(!indicator_name %in% colnames(samples)){
      stop(paste0("The '", indicator_name, "' indicator was not found in the posterior samples."), call. = FALSE)
    }
    as.numeric(samples[, indicator_name] %in% indicator_values)
  })
}

.indicator_MCMC_summary <- function(indicator_list){

  indicator_list <- lapply(indicator_list, as.numeric)
  indicator      <- unlist(indicator_list, use.names = FALSE)

  indicator_mcmc <- coda::as.mcmc.list(lapply(indicator_list, function(x){
    coda::as.mcmc(matrix(x, ncol = 1, dimnames = list(NULL, "indicator")))
  }))

  mcmc_summary <- summary(indicator_mcmc, quantiles = NULL)$statistics
  if(is.null(dim(mcmc_summary))){
    mcmc_summary <- t(mcmc_summary)
  }

  ESS <- coda::effectiveSize(indicator_mcmc)[1]
  if(is.nan(ESS)){
    ESS <- 0
  }

  return(list(
    post_prob       = mean(indicator),
    MCMC_error      = mcmc_summary[1,"Time-series SE"],
    MCMC_SD_error   = mcmc_summary[1,"Time-series SE"] / mcmc_summary[1,"SD"],
    ESS             = ESS,
    visits          = sum(indicator),
    n_samples       = length(indicator)
  ))
}

.indicator_BF_diagnostics <- function(indicator_list, prior_prob, BF){

  diagnostics <- .indicator_MCMC_summary(indicator_list)

  if(diagnostics$post_prob <= 0 || diagnostics$post_prob >= 1 ||
     prior_prob <= 0 || prior_prob >= 1 ||
     !is.finite(BF) || !is.finite(diagnostics$MCMC_error)){
    BF_error_percent <- NA_real_
  }else{
    log_BF_error <- diagnostics$MCMC_error / (diagnostics$post_prob * (1 - diagnostics$post_prob))
    BF_error_percent <- 100 * log_BF_error
  }

  diagnostics$BF_error_percent <- BF_error_percent

  return(diagnostics)
}

.indicator_BF_diagnostic_row <- function(diagnostics){

  data.frame(
    ESS              = diagnostics$ESS,
    MCMC_error       = diagnostics$MCMC_error,
    BF_error_percent = diagnostics$BF_error_percent
  )
}

.indicator_BF_reporting_value <- function(BF, post_prob, prior_prob, n_samples){

  out <- list(value = BF, operator = NA_character_)
  if(length(post_prob) != 1L || length(prior_prob) != 1L || length(n_samples) != 1L){
    return(out)
  }
  if(!is.finite(post_prob) || !is.finite(prior_prob) || !is.finite(n_samples) ||
     prior_prob <= 0 || prior_prob >= 1 || n_samples <= 1){
    return(out)
  }

  prior_odds <- prior_prob / (1 - prior_prob)
  if(post_prob >= 1){
    out$value    <- (n_samples - 1) / prior_odds
    out$operator <- ">"
  }else if(post_prob <= 0){
    out$value    <- 1 / ((n_samples - 1) * prior_odds)
    out$operator <- "<"
  }

  return(out)
}

.indicator_BF_warnings <- function(parameter, diagnostics, severe_limit = 20, warning_limit = 100){

  minority_visits <- min(diagnostics$visits, diagnostics$n_samples - diagnostics$visits)

  if(minority_visits == 0){
    return()
  }else if(minority_visits < severe_limit){
    return(stats::setNames(
      sprintf("Bayes factor MC error for %s is based on only %i posterior samples from the less frequent model.", parameter, minority_visits),
      parameter
    ))
  }else if(minority_visits < warning_limit){
    return(stats::setNames(
      sprintf("Bayes factor MC error for %s is based on %i posterior samples from the less frequent model.", parameter, minority_visits),
      parameter
    ))
  }else{
    return()
  }
}

.BF_error_column_name <- function(BF01 = FALSE){

  "error%(Inclusion BF)"
}

.standardize_BF_bound_operator <- function(operator, n){

  if(is.null(operator)){
    return(rep(NA_character_, n))
  }
  operator <- as.character(operator)
  if(length(operator) != n){
    operator <- rep_len(operator, n)
  }
  invalid <- !is.na(operator) & !operator %in% c("<", ">")
  if(any(invalid)){
    stop("BF bound operators must be '<' or '>'.", call. = FALSE)
  }

  return(operator)
}

.invert_BF_bound_operator <- function(operator){

  operator <- as.character(operator)
  operator[operator == "<"] <- "TEMP_GT"
  operator[operator == ">"] <- "<"
  operator[operator == "TEMP_GT"] <- ">"

  return(operator)
}

.match_table_columns <- function(output_names, input_names){

  if(length(output_names) == 0){
    return(integer())
  }

  used <- rep(FALSE, length(input_names))
  vapply(output_names, function(output_name){
    matches <- which(input_names == output_name & !used)
    if(length(matches) == 0){
      return(NA_integer_)
    }
    used[matches[1]] <<- TRUE
    matches[1]
  }, integer(1))
}

.restore_table_column_attributes <- function(column, source_attributes){

  source_attributes <- source_attributes[!names(source_attributes) %in% c("names", "dim", "dimnames")]
  for(attribute in names(source_attributes)){
    attr(column, attribute) <- source_attributes[[attribute]]
  }

  column
}

.subset_table_parameters <- function(table, output){

  parameters <- attr(table, "parameters")
  if(is.null(parameters) || length(parameters) != nrow(table)){
    return(NULL)
  }

  row_indices <- match(rownames(output), rownames(table))
  if(any(is.na(row_indices))){
    return(NULL)
  }

  parameters[row_indices]
}

.subset_table_warnings <- function(warnings, selected_parameters = NULL, selected_rows = NULL){

  if(is.null(warnings) || length(warnings) == 0){
    return(warnings)
  }

  warning_names <- names(warnings)
  if(is.null(warning_names)){
    return(warnings)
  }
  warning_names[is.na(warning_names)] <- ""
  if(!any(nzchar(warning_names))){
    return(warnings)
  }

  selected <- unique(c(selected_parameters, selected_rows))
  if(length(selected) == 0){
    return(warnings[!nzchar(warning_names)])
  }

  keep <- !nzchar(warning_names) | warning_names %in% selected
  warnings[keep]
}

.subset_table_hypothesis_attributes <- function(table, output){

  if(!inherits(table, "BayesTools_hypothesis_BF")){
    return(output)
  }

  raw_BF <- attr(table, "raw_BF")
  attr(output, "raw_BF") <- NULL
  if(!is.null(raw_BF) && length(raw_BF) == nrow(table)){
    row_indices <- match(rownames(output), rownames(table))
    if(length(row_indices) == nrow(output) &&
       !any(is.na(row_indices)) &&
       !anyDuplicated(row_indices)){
      attr(output, "raw_BF") <- raw_BF[row_indices]
    }
  }

  output
}

.format_BF_column <- function(x){

  out <- format(round(x, digits = 3), nsmall = 3)
  bound_operator <- .standardize_BF_bound_operator(attr(x, "bound_operator"), length(x))
  has_bound <- !is.na(bound_operator) & !is.na(x)
  if(any(has_bound)){
    out[has_bound] <- paste0(bound_operator[has_bound], trimws(out[has_bound]))
  }

  return(out)
}

.copy_BayesTools_table_attributes <- function(new_table, table, type){

  class(new_table)        <- class(table)
  attr(new_table, "type") <- type

  copied_attributes <- setdiff(
    names(attributes(table)),
    c("class", "type", "names", "dim", "dimnames")
  )
  for(a in copied_attributes){
    attr(new_table, a) <- attr(table, a)
  }

  new_table
}

#' @title Adds column to BayesTools table
#'
#' @description Adds column to a BayesTools table while not
#' breaking formatting, attributes, etc...
#'
#' @param table BayesTools table
#' @param column_title title of the new column
#' @param column_values values of the new column
#' @param column_position position of the new column (defaults to \code{NULL} which
#' appends the column to the end)
#' @param column_type type of values of the new column table (important for formatting,
#' defaults to \code{NULL} = the function tries to guess numeric / character based on the
#' \code{column_values} but many more specific types are available)
#'
#' @return returns an object of 'BayesTools_table' class.
#'
#' @export
add_column <- function(table, column_title, column_values, column_position = NULL, column_type = NULL){

  if(!inherits(table, "BayesTools_table"))
    stop("The 'table' must be of class 'BayesTools_table'.")
  check_char(column_title, "column_title")
  if(!(is.vector(column_values) | is.numeric(column_values)) || length(column_values) != nrow(table))
    stop("The 'column_values' must be a vector of the same length as has the table rows.")
  check_int(column_position, "column_position", allow_NULL = TRUE, lower = 0, upper = ncol(table) + 1)
  .check_table_types(column_type, "column_type", allow_NULL = TRUE)

  # fill defaults
  if(is.null(column_position)){
    column_position <- ncol(table) + 1
  }
  if(is.null(column_type)){
    if(is.numeric(column_values) && !all(is.na(column_values)) && all(.is.wholenumber(column_values, na.rm = TRUE))){
      column_type <- "integer"
    }else if(is.numeric(column_values)){
      column_type <- "estimate"
    }else if(is.character(column_values)){
      column_type <- "string"
    }else{
      stop("The 'column_type' could not be guessed. Please, supply it manually.")
    }
  }

  before <- if(column_position > 1) seq_len(column_position - 1) else integer(0)
  after  <- if(column_position <= ncol(table)) seq.int(max(column_position, 1L), ncol(table)) else integer(0)

  new_table <- table
  new_table[[ncol(table) + 1L]] <- column_values
  new_table <- new_table[, c(before, ncol(table) + 1L, after), drop = FALSE]

  names(new_table) <- c(
    if(length(before) > 0) colnames(table)[before],
    column_title,
    if(length(after) > 0) colnames(table)[after]
  )

  .copy_BayesTools_table_attributes(new_table, table, c(
    if(length(before) > 0) attr(table, "type")[before],
    column_type,
    if(length(after) > 0) attr(table, "type")[after]
  ))
}

#' @title Removes column to BayesTools table
#'
#' @description Removes column to a BayesTools table while not
#' breaking formatting, attributes, etc...
#'
#' @param table BayesTools table
#' @param column_position position of the to be removed column (defaults to \code{NULL} which
#' removes the last column)
#'
#' @return returns an object of 'BayesTools_table' class.
#'
#' @export
remove_column <- function(table, column_position = NULL){

  if(!inherits(table, "BayesTools_table"))
    stop("The 'table' must be of class 'BayesTools_table'.")
  check_int(column_position, "column_position", allow_NULL = TRUE, lower = 1, upper = ncol(table))

  # fill defaults
  if(is.null(column_position)){
    column_position <- ncol(table)
  }

  keep <- setdiff(seq_len(ncol(table)), column_position)
  new_table <- table[, keep, drop = FALSE]

  # transfer column names
  colnames(new_table) <- colnames(table)[keep]

  .copy_BayesTools_table_attributes(new_table, table, attr(table, "type")[keep])
}

#' @title Updates BayesTools table
#'
#' @description Updates BayesTools table while not breaking formatting, attributes, etc...
#'
#' @param object a BayesTools table
#' @param title title of the table
#' @param footnotes add footnotes to the table
#' @param warnings add warnings of the table
#' @param remove_parameters remove parameters from the table
#' @param logBF whether to format Bayes factors as log(BF)
#' @param BF01 whether to format Bayes factors as 1/BF
#' @param ... additional arguments.
#'
#' @return returns an object of 'BayesTools_table' class.
#' @export
update.BayesTools_table <- function(object, title = NULL, footnotes = NULL, warnings = NULL, remove_parameters = NULL, logBF = FALSE, BF01 = FALSE, ...){

  check_char(title, "title", allow_NULL = TRUE)
  check_char(footnotes, "footnotes", check_length = 0, allow_NULL = TRUE)
  check_char(warnings, "warnings", check_length = 0, allow_NULL = TRUE)
  check_char(remove_parameters, "remove_parameters", check_length = 0, allow_NULL = TRUE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01,  "BF01",  allow_NA = FALSE)

  if(!is.null(footnotes)){
    attr(object, "footnotes") <- c(attr(object, "footnotes"), footnotes)
  }

  if(!is.null(warnings)){
    attr(object, "warnings")  <- c(attr(object, "warnings"), warnings)
  }

  if(!is.null(remove_parameters)){
    object <- object[!rownames(object) %in% remove_parameters,,drop=FALSE]
  }

  if(!is.null(title)){
    attr(object, "title") <- title
  }

  BF_types <- attr(object, "type") %in% c("BF", "inclusion_BF")
  if(any(BF_types)){
    for(BF_col in which(BF_types)){
      BF_values <- object[[BF_col]]
      raw_BF    <- as.numeric(BF_values)
      bound_operator <- .standardize_BF_bound_operator(attr(BF_values, "bound_operator"), length(BF_values))
      if(isTRUE(attr(BF_values, "logBF"))){
        raw_BF <- exp(raw_BF)
      }
      if(isTRUE(attr(BF_values, "BF01"))){
        raw_BF <- 1 / raw_BF
        bound_operator <- .invert_BF_bound_operator(bound_operator)
      }
      attr(raw_BF, "bound_operator") <- bound_operator
      object[[BF_col]] <- format_BF(
        raw_BF,
        logBF     = logBF,
        BF01      = BF01,
        inclusion = identical(attr(object, "type")[BF_col], "inclusion_BF")
      )
    }
  }
  if(any(attr(object, "type") == "BF_error")){
    BF_error_cols <- which(attr(object, "type") == "BF_error")
    if(length(BF_error_cols) == 1){
      BF_cols <- which(attr(object, "type") %in% c("BF", "inclusion_BF"))
      if(length(BF_cols) == 1 && identical(attr(object, "type")[BF_cols], "BF")){
        attr(object[[BF_error_cols]], "name") <- "error%(BF)"
      }else{
        attr(object[[BF_error_cols]], "name") <- .BF_error_column_name(BF01)
      }
    }
  }
  attr(object, "logBF") <- logBF
  attr(object, "BF01")  <- BF01

  return(object)
}
