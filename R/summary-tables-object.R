#' @title Print a BayesTools table
#'
#' @param x a BayesTools_values_tables
#' @param ... additional arguments.
#' @details Finite positive Bayes-factor bounds that would round to zero at
#' three decimals are printed as reciprocals (for example \code{<1/10000}),
#' or as \code{exp()} of their logarithm when the reciprocal overflows.
#' Log-scale bounds that would round to zero keep three significant digits,
#' whichever sign they have (for example \code{>-1e-05}). This
#' changes only their presentation; numeric values, BF direction and bound
#' operators are retained.
#' A validated finite-log endpoint whose natural value in the displayed
#' orientation is not representable (\code{0} or \code{Inf}) is printed as
#' \code{exp(log_endpoint)}.
#'
#' @return \code{print.BayesTools_table} returns \code{NULL}.
#'
#' @exportS3Method
print.BayesTools_table <- function(x, ...){

  # print formatting
  # n_models is row-aligned (one denominator per table row), so the Models
  # column must see the full vector rather than n_models[[column]].
  for(i in seq_along(attr(x, "type"))){
    values         <- x[[i]]
    colnames(x)[i] <- .format_column_names(colnames(x)[i], attr(x, "type")[i], values)
    x[[i]]         <- .format_column(values, attr(x, "type")[i], attr(x, "n_models"))
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
  class(x) <- "data.frame"
  print(
    x,
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
#' @details Valid live Bayes-factor carriers retain their canonical BF10
#' logarithms when reformatted into either direction and scale. This preserves
#' finite logs and representable reciprocals when the natural display is zero
#' or infinite. Replaced or arithmetically modified vectors follow ordinary
#' numeric formatting.
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

  canonical <- .BF_carrier_log(BF)
  if(!is.null(canonical)){
    check_bool(logBF, "logBF", allow_NA = FALSE)
    check_bool(BF01, "BF01", allow_NA = FALSE)
    check_bool(inclusion, "inclusion", allow_NA = FALSE)
    names(canonical) <- names(BF)
    bound <- .standardize_BF_bound_operator(attr(BF, "bound_operator"), length(BF))
    if(isTRUE(attr(BF, "BF01", exact = TRUE))) bound <- .invert_BF_bound_operator(bound)
    return(.format_BF_from_log(canonical, logBF, BF01, inclusion,
      bound_operator = bound, diagnostics = attr(BF, "numerical_diagnostics", exact = TRUE)))
  }

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

# Format Bayes factors supplied as natural-log BF values in the BF10
# orientation. Log-scale outputs are computed in log space, so Bayes factors
# beyond the double range keep finite log values instead of passing through
# exp() to 0 or Inf. Natural-scale outputs use 'BF' when it is supplied.
.format_BF_from_log <- function(log_BF, logBF = FALSE, BF01 = FALSE, inclusion = FALSE,
                                bound_operator = NULL, BF = NULL, diagnostics = NULL){

  BF_names <- names(log_BF)
  log_BF <- as.numeric(log_BF)
  BF <- if(is.null(BF)) exp(log_BF) else as.numeric(BF)
  names(BF) <- BF_names
  attr(BF, "bound_operator") <- bound_operator
  out <- format_BF(BF, logBF = logBF, BF01 = BF01, inclusion = inclusion)
  if(logBF){
    out[] <- if(BF01) -log_BF else log_BF
  }else if(BF01){
    hidden <- is.finite(log_BF) & (BF == 0 | is.infinite(BF))
    out[hidden] <- exp(-log_BF[hidden])
  }
  attr(out, "canonical_log_BF") <- log_BF
  attr(out, "producer_values") <- as.numeric(out)
  attr(out, "carrier_logBF") <- logBF
  attr(out, "carrier_BF01") <- BF01
  attr(out, "numerical_diagnostics") <- diagnostics
  class(out) <- unique(c("BayesTools_BF", class(out)))
  out
}

#' @title Format Bayes factors from natural-log values
#'
#' @description Builds a Bayes factor from natural-logarithm BF10 values without
#' passing them through \code{exp()}. Output on the logarithmic scale
#' (\code{logBF = TRUE}) and in the inverse orientation (\code{BF01 = TRUE}) is
#' computed in log space, so Bayes factors beyond the double range (an
#' \code{exp(log_BF)} equal to \code{0} or \code{Inf}) keep finite logarithms
#' and representable reciprocals. [format_BF()] cannot recover them from the
#' natural-scale numbers.
#'
#' @details The result is the Bayes-factor carrier that [format_BF()] returns
#' when it receives such a carrier: reformatting it with [format_BF()] keeps
#' the canonical logarithms. A plain numeric vector given to [format_BF()]
#' returns an unclassed vector (it has the class \code{BayesTools_BF} only
#' while a non-missing bound operator is present). The carrier stores the
#' logarithms in private attributes; [BF_carrier_drop()] removes them.
#' Subsetting keeps the logarithms of the selected values, arithmetic on the
#' carrier returns ordinary numbers, and replacing values removes the carrier
#' attributes.
#'
#' @param log_BF numeric vector of natural-logarithm BF10 values (the Bayes
#' factor in favor of the alternative hypothesis). Names are kept. Missing
#' values (\code{NA}, \code{NaN}) are propagated.
#' @param logBF whether the logarithm of the Bayes factor is returned.
#' @param BF01 whether the inverse Bayes factor (the BF01 orientation) is
#' returned.
#' @param inclusion whether the Bayes factor is an inclusion Bayes factor (for
#' naming purposes).
#' @param bound_operator \code{NULL} (all values exact), or a character vector
#' with one element per value of \code{log_BF}: \code{NA} for an exact value,
#' \code{"<"} or \code{">"} when the value is an upper or lower bound of the
#' Bayes factor. The operators refer to the BF10 values in \code{log_BF} and
#' are stored in the orientation of the output (inverted when
#' \code{BF01 = TRUE}, see [invert_BF_bound_operator()]).
#' @param diagnostics \code{NULL}, or a list with one element per value of
#' \code{log_BF} (elements can be \code{NULL}) holding numerical diagnostics of
#' the individual Bayes factors, stored with the carrier.
#'
#' @return \code{format_BF_from_log()} returns a numeric vector of class
#' \code{c("BayesTools_BF", "numeric")} with the names of \code{log_BF}. Its
#' public attributes are \code{name} (for example \code{"BF"}, \code{"1/BF"},
#' \code{"Inclusion BF"}, \code{"Exclusion BF"}, wrapped in \code{log()} when
#' \code{logBF = TRUE}), \code{logBF}, \code{BF01}, and \code{bound_operator}.
#'
#' @seealso [format_BF()] [BF_carrier_drop()] [invert_BF_bound_operator()]
#'
#' @examples
#' # a Bayes factor beyond the double range keeps its finite logarithm
#' BF <- format_BF_from_log(c(a = 2, b = 800))
#' format_BF(BF, logBF = TRUE)
#' format_BF(BF, BF01 = TRUE, logBF = TRUE)
#'
#' @export
format_BF_from_log <- function(log_BF, logBF = FALSE, BF01 = FALSE, inclusion = FALSE,
                               bound_operator = NULL, diagnostics = NULL){

  check_real(log_BF, "log_BF", check_length = 0, allow_NA = TRUE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01, "BF01", allow_NA = FALSE)
  check_bool(inclusion, "inclusion", allow_NA = FALSE)
  if(!is.null(bound_operator)){
    .BF_bound_operator_check(bound_operator, "bound_operator")
    if(length(bound_operator) != length(log_BF)){
      stop("The 'bound_operator' argument must have one element per value of 'log_BF'.", call. = FALSE)
    }
  }
  if(!is.null(diagnostics) && (!is.list(diagnostics) || length(diagnostics) != length(log_BF))){
    stop("The 'diagnostics' argument must be NULL or a list with one element per value of 'log_BF'.", call. = FALSE)
  }

  .format_BF_from_log(log_BF, logBF, BF01, inclusion, bound_operator,
    BF = NULL, diagnostics = diagnostics)
}

.BF_carrier_attributes <- function(){

  c("canonical_log_BF", "producer_values", "carrier_logBF", "carrier_BF01", "numerical_diagnostics")
}

.BF_carrier_log <- function(x){

  logs <- attr(x, "canonical_log_BF", exact = TRUE)
  if(!inherits(x, "BayesTools_BF") || length(logs) != length(x) ||
     !identical(attr(x, "producer_values", exact = TRUE), as.numeric(x)) ||
     !identical(attr(x, "carrier_logBF", exact = TRUE), attr(x, "logBF", exact = TRUE)) ||
     !identical(attr(x, "carrier_BF01", exact = TRUE), attr(x, "BF01", exact = TRUE))) return(NULL)
  logs
}

.BF_carrier_drop <- function(x){

  for(attribute in .BF_carrier_attributes()) attr(x, attribute) <- NULL
  for(attribute in c("log_BF", "log_prior_height", "log_posterior_height", "orientation")) attr(x, attribute) <- NULL
  if(!any(!is.na(attr(x, "bound_operator")))) class(x) <- setdiff(class(x), "BayesTools_BF")
  x
}

#' @title Remove the carrier attributes of a Bayes factor
#'
#' @description Removes the private attributes with which a Bayes-factor
#' carrier (see [format_BF_from_log()]) retains its canonical logarithms, so
#' that the values can be stored in a plain data frame or compared as ordinary
#' numbers.
#'
#' @details The names and the public attributes \code{name}, \code{logBF},
#' \code{BF01}, and \code{bound_operator} are kept. The class
#' \code{BayesTools_BF} is kept only while a non-missing bound operator remains,
#' as for the result of [format_BF()] on a plain vector. A numeric vector
#' without carrier attributes is returned with the same rule applied to its
#' class.
#'
#' @param x a numeric vector, typically a Bayes factor created by
#' [format_BF_from_log()] or [format_BF()].
#'
#' @return \code{BF_carrier_drop()} returns \code{x} without its carrier
#' attributes.
#'
#' @seealso [format_BF_from_log()] [format_BF()]
#'
#' @examples
#' BF <- format_BF_from_log(c(a = 2, b = 800))
#' names(attributes(BF))
#' names(attributes(BF_carrier_drop(BF)))
#'
#' @export
BF_carrier_drop <- function(x){

  if(!is.numeric(x)){
    stop("The 'x' argument must be a numeric vector.", call. = FALSE)
  }

  .BF_carrier_drop(x)
}

#' @export
`[.BayesTools_BF` <- function(x, i, ...){

  logs <- .BF_carrier_log(x)
  positions <- stats::setNames(seq_along(x), names(x))
  index <- if(missing(i)) positions else positions[i]
  index <- unname(index)
  out <- NextMethod("[")
  for(attribute in c("name", "logBF", "BF01")) attr(out, attribute) <- attr(x, attribute)
  bound <- attr(x, "bound_operator")
  if(!is.null(bound)) attr(out, "bound_operator") <- bound[index]
  out <- .BF_carrier_drop(out)
  if(!is.null(logs)){
    attr(out, "canonical_log_BF") <- logs[index]
    attr(out, "producer_values") <- as.numeric(out)
    attr(out, "carrier_logBF") <- attr(x, "logBF")
    attr(out, "carrier_BF01") <- attr(x, "BF01")
    diagnostics <- attr(x, "numerical_diagnostics")
    if(length(diagnostics) == length(x)) attr(out, "numerical_diagnostics") <- diagnostics[index]
    for(attribute in c("log_BF", "log_prior_height", "log_posterior_height")){
      values <- attr(x, attribute, exact = TRUE)
      if(length(values) == length(x)) attr(out, attribute) <- values[index]
    }
    attr(out, "orientation") <- attr(x, "orientation", exact = TRUE)
    class(out) <- unique(c("BayesTools_BF", class(out)))
  }
  if(any(!is.na(attr(out, "bound_operator", exact = TRUE)))){
    class(out) <- unique(c("BayesTools_BF", class(out)))
  }
  out
}

#' @export
`[<-.BayesTools_BF` <- function(x, i, ..., value){

  .BF_carrier_drop(NextMethod("[<-"))
}

#' @export
`[[<-.BayesTools_BF` <- function(x, i, ..., value){

  .BF_carrier_drop(NextMethod("[[<-"))
}

#' @export
Ops.BayesTools_BF <- function(e1, e2){

  if(inherits(e1, "BayesTools_BF")){
    values <- as.numeric(e1)
    names(values) <- names(e1)
    e1 <- values
  }
  if(!missing(e2) && inherits(e2, "BayesTools_BF")){
    values <- as.numeric(e2)
    names(values) <- names(e2)
    e2 <- values
  }
  if(missing(e2)) do.call(.Generic, list(e1)) else do.call(.Generic, list(e1, e2))
}

#' @export
Math.BayesTools_BF <- function(x, ...){

  values <- as.numeric(x)
  names(values) <- names(x)
  do.call(.Generic, c(list(values), list(...)))
}

#' @export
`[.BayesTools_table` <- function(x, i, j, ..., drop = TRUE){

  original_names <- names(x)
  original_col_attributes <- lapply(x, attributes)
  row_indices <- if(missing(i) || (missing(j) && (nargs() - !missing(drop)) < 3L)) seq_len(nrow(x)) else {
    if(is.character(i)) pmatch(i, rownames(x), duplicates.ok = TRUE) else seq_len(nrow(x))[i]
  }

  out <- NextMethod("[")
  if(!is.data.frame(out)){
    return(out)
  }

  # the source column of every selected column, resolved from the selector
  # (before the data frame makes repeated names unique): a column selected
  # twice keeps its own carrier and metadata
  column_positions <- stats::setNames(seq_along(original_names), original_names)
  select_columns   <- if((nargs() - !missing(drop)) < 3L) !missing(i) else !missing(j)
  source_cols <- if(!select_columns){
    unname(column_positions)
  }else if((nargs() - !missing(drop)) < 3L){
    unname(column_positions[i])
  }else{
    unname(column_positions[j])
  }
  valid_cols  <- !is.na(source_cols)

  if(length(source_cols) > 0 && any(valid_cols)){
    for(k in which(valid_cols)){
      source <- x[[source_cols[k]]]
      if(inherits(source, "BayesTools_BF") && length(row_indices) == nrow(out)){
        out[[k]] <- source[row_indices]
      }else{
        out[[k]] <- .restore_table_column_attributes(out[[k]], original_col_attributes[[source_cols[k]]])
      }
    }
  }

  type <- attr(x, "type")
  if(length(type) == length(original_names)){
    attr(out, "type") <- type[source_cols]
  }

  n_models <- attr(x, "n_models")
  if(!is.null(n_models) && length(n_models) == nrow(x)){
    if(length(row_indices) == nrow(out)){
      attr(out, "n_models") <- n_models[row_indices]
    }
  }

  for(attribute in c("title", "footnotes", "rownames")){
    attr(out, attribute) <- attr(x, attribute)
  }
  attr(out, "footnotes") <- .subset_table_row_footnotes(
    attr(x, "footnotes"),
    table_rows = rownames(x),
    selected_rows = rownames(out)
  )

  selected_parameters <- .subset_table_parameters(x, out, row_indices)
  if(!is.null(selected_parameters)){
    attr(out, "parameters") <- selected_parameters
  }
  attr(out, "quantities") <- .subset_table_quantities(x, out, row_indices)
  attr(out, "warnings") <- .subset_table_warnings(attr(x, "warnings"), selected_parameters, rownames(out))
  out <- .subset_table_hypothesis_attributes(x, out, row_indices)

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

  if(BF01) "error%(Exclusion BF)" else "error%(Inclusion BF)"
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

#' @title Invert the orientation of Bayes factor bound operators
#'
#' @description Swaps \code{"<"} and \code{">"}: the bound operators of a Bayes
#' factor in the BF10 orientation become those of its reciprocal in the BF01
#' orientation, and the other way round. Missing operators (exact values) are
#' kept.
#'
#' @param operator a character vector of bound operators (\code{"<"}, \code{">"}
#' or \code{NA}), or a vector of only \code{NA} values. Other values are
#' rejected.
#'
#' @return \code{invert_BF_bound_operator()} returns a character vector of the
#' length of \code{operator} (without names) with the operators swapped.
#'
#' @seealso [format_BF_from_log()] [format_BF()]
#'
#' @examples
#' invert_BF_bound_operator(c("<", ">", NA))
#'
#' @export
invert_BF_bound_operator <- function(operator){

  .BF_bound_operator_check(operator, "operator")

  .invert_BF_bound_operator(operator)
}

# Bound operators of the public Bayes factor helpers: a character vector (or a
# vector of only NA) whose non-missing values are '<' or '>'.
.BF_bound_operator_check <- function(operator, name){

  if(!(is.character(operator) || (is.logical(operator) && all(is.na(operator))))){
    stop("The '", name, "' argument must be a character vector of '<' and '>' bound operators (NA for exact values).",
         call. = FALSE)
  }
  invalid <- !is.na(operator) & !operator %in% c("<", ">")
  if(any(invalid)){
    stop("The '", paste0(unique(operator[invalid]), collapse = "', '"), "' values are not recognized by the '",
         name, "' argument. Bound operators must be '<', '>', or NA.", call. = FALSE)
  }

  invisible(TRUE)
}

.restore_table_column_attributes <- function(column, source_attributes){

  source_attributes <- source_attributes[!names(source_attributes) %in%
    c("names", "dim", "dimnames", "bound_operator", .BF_carrier_attributes())]
  for(attribute in names(source_attributes)){
    attr(column, attribute) <- source_attributes[[attribute]]
  }

  column
}

.subset_table_parameters <- function(table, output, row_indices = NULL){

  parameters <- attr(table, "parameters")
  if(is.null(parameters) || length(parameters) != nrow(table)){
    return(NULL)
  }

  if(is.null(row_indices)) row_indices <- match(rownames(output), rownames(table))
  if(length(row_indices) != nrow(output)){
    return(NULL)
  }

  parameters[row_indices]
}

# The per-row quantity table of the rows kept by subsetting (NULL when the
# table has none or the rows cannot be matched).
.subset_table_quantities <- function(table, output, row_indices = NULL){

  quantities <- attr(table, "quantities", exact = TRUE)
  if(is.null(quantities) || nrow(quantities) != nrow(table)){
    return(NULL)
  }
  if(is.null(row_indices)) row_indices <- match(rownames(output), rownames(table))
  if(length(row_indices) != nrow(output)){
    return(NULL)
  }
  out <- quantities[row_indices, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The per-row quantity table of an estimates table (its 'quantities'
# attribute): the row labels, the parameter_catalog() quantity id of each row
# ("" when the row is not a catalog quantity), and the label parts the row
# label is rendered from.
.bt_table_quantities <- function(rows, label_parts, quantity_ids = NULL){

  if(is.null(quantity_ids)){
    quantity_ids <- rep("", length(rows))
  }
  if(length(label_parts) != length(rows) || length(quantity_ids) != length(rows)){
    stop("The label parts of the table rows do not match its rows.", call. = FALSE)
  }
  out <- data.frame(
    row         = as.character(rows),
    quantity_id = as.character(quantity_ids),
    stringsAsFactors = FALSE
  )
  out$label_parts <- I(unname(label_parts))
  out
}

# The catalog quantity ids of label parts: the quantity whose canonical name
# is the rendered selector, otherwise the quantity of a unique exact alias
# (e.g. a transformed factor level `<parameter>[dif: level]`), in the
# namespace of the parts' formula parameter; "" when there is none, and for
# parts that record an output transformation of the values (the values are
# no catalog quantity).
.bt_table_quantity_ids <- function(label_parts, catalog){

  out <- rep("", length(label_parts))
  if(is.null(catalog) || length(label_parts) == 0L){
    return(out)
  }
  selectors <- .bt_label(label_parts, style = "selector")
  namespaces <- vapply(label_parts, function(part){
    if(nzchar(part$formula_parameter)) part$formula_parameter else "model"
  }, character(1))
  out <- .bt_draws_quantity_ids(label_parts, catalog)
  missing <- !nzchar(out)
  aliases <- catalog$aliases[!catalog$aliases$simplified, , drop = FALSE]
  if(any(missing) && nrow(aliases) > 0L){
    pairs <- unique(data.frame(
      key = paste(aliases$alias, aliases$namespace, sep = "\r"),
      id  = aliases$quantity_id,
      stringsAsFactors = FALSE
    ))
    pairs <- pairs[!pairs$key %in% pairs$key[duplicated(pairs$key)], , drop = FALSE]
    matched <- pairs$id[match(paste(selectors, namespaces, sep = "\r")[missing], pairs$key)]
    out[missing] <- ifelse(is.na(matched), "", matched)
  }
  transformed <- vapply(label_parts, function(part){
    length(.bt_label_output_transformation(part)) > 0L
  }, logical(1))
  out[transformed] <- ""

  out
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

# Footnotes named after a table row (such as the share of draws in which an
# original-scale correlation is defined) describe that row and are dropped
# with it. Unnamed footnotes and footnotes named otherwise are kept.
.subset_table_row_footnotes <- function(footnotes, table_rows, selected_rows){

  footnote_names <- names(footnotes)
  if(length(footnotes) == 0L || is.null(footnote_names)){
    return(footnotes)
  }

  dropped <- footnote_names %in% setdiff(table_rows, selected_rows)
  if(!any(dropped)){
    return(footnotes)
  }
  footnotes <- footnotes[!dropped]
  if(length(footnotes) == 0L){
    return(NULL)
  }
  if(!any(nzchar(names(footnotes)))){
    footnotes <- unname(footnotes)
  }

  footnotes
}

.subset_table_hypothesis_attributes <- function(table, output, row_indices = NULL){

  explicit_rows <- !is.null(row_indices)
  if(is.null(row_indices)) row_indices <- match(rownames(output), rownames(table))
  for(attribute in c("raw_BF", "raw_log_BF", "numerical_diagnostics", "prior_numerical_diagnostics")){
    values <- attr(table, attribute, exact = TRUE)
    attr(output, attribute) <- NULL
    if(attribute %in% c("raw_BF", "raw_log_BF") &&
       !any(attr(output, "type") %in% c("BF", "inclusion_BF"))) next
    if(length(values) == nrow(table) && length(row_indices) == nrow(output) &&
       (explicit_rows || !anyNA(row_indices))){
      attr(output, attribute) <- values[row_indices]
    }
  }
  output
}


# A positive Bayes factor that rounds to zero at the printed digits is shown as
# its reciprocal (1/10000), mirroring large values; a reciprocal that overflows
# is shown as exp(<log>) with the natural logarithm of the value.
.format_BF_small_value <- function(BF, log_BF){

  reciprocal <- 1 / BF
  if(is.finite(reciprocal)){
    paste0("1/", format(round(reciprocal), trim = TRUE))
  }else{
    paste0("exp(", format(log_BF, trim = TRUE), ")")
  }
}

.format_BF_column <- function(x){

  out <- format(round(x, digits = 3), nsmall = 3)
  bound_operator <- .standardize_BF_bound_operator(attr(x, "bound_operator"), length(x))
  log_scale <- isTRUE(attr(x, "logBF"))
  # A bound that would print as 0.000: a log-BF bound of either sign, and a
  # positive natural-scale bound.
  small_bound <- !is.na(bound_operator) & is.finite(x) & round(x, 3) == 0 & (x > 0 | (log_scale & x < 0))
  logs <- .BF_carrier_log(x)
  if(any(small_bound) && log_scale){
    # A small log-BF bound keeps significant digits.
    out[small_bound] <- format(x[small_bound], digits = 3, scientific = TRUE, trim = TRUE)
  }else if(any(small_bound)){
    # Shown as reciprocals (1/10000), mirroring large bounds (>10000.000);
    # an overflowing reciprocal shows the (canonical) logarithm instead.
    small_logs <- if(!is.null(logs) && !isTRUE(attr(x, "logBF"))){
      (if(isTRUE(attr(x, "BF01"))) -logs else logs)[small_bound]
    }else{
      log(x[small_bound])
    }
    out[small_bound] <- vapply(seq_len(sum(small_bound)), function(i){
      .format_BF_small_value(x[small_bound][i], small_logs[i])
    }, character(1))
  }
  if(!is.null(logs) && !log_scale){
    # A bound whose oriented natural value is 0 or Inf, while its logarithm is
    # finite, is shown as exp(<log>).
    oriented_logs <- if(isTRUE(attr(x, "BF01"))) -logs else logs
    hidden_bound <- !is.na(bound_operator) & is.finite(oriented_logs) & !is.na(x) & (x == 0 | is.infinite(x))
    if(any(hidden_bound, na.rm = TRUE)){
      out[hidden_bound] <- paste0("exp(", format(oriented_logs[hidden_bound], trim = TRUE), ")")
    }
  }
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

  .subset_table_hypothesis_attributes(table, new_table, seq_len(nrow(table)))
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
      canonical <- .BF_carrier_log(BF_values)
      bound_operator <- .standardize_BF_bound_operator(attr(BF_values, "bound_operator"), length(BF_values))
      inclusion <- identical(attr(object, "type")[BF_col], "inclusion_BF")
      if(!is.null(canonical)){
        if(isTRUE(attr(BF_values, "BF01"))) bound_operator <- .invert_BF_bound_operator(bound_operator)
        object[[BF_col]] <- .format_BF_from_log(canonical, logBF, BF01, inclusion,
          bound_operator = bound_operator, diagnostics = attr(BF_values, "numerical_diagnostics"))
        next
      }
      if(isTRUE(attr(BF_values, "logBF"))){
        # stay in log space so that log-scale values survive re-formatting
        if(isTRUE(attr(BF_values, "BF01"))){
          raw_BF <- -raw_BF
          bound_operator <- .invert_BF_bound_operator(bound_operator)
        }
        object[[BF_col]] <- .format_BF_from_log(
          raw_BF,
          logBF          = logBF,
          BF01           = BF01,
          inclusion      = inclusion,
          bound_operator = bound_operator
        )
        next
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
        inclusion = inclusion
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

# Plain coercions keep the existing layout and strip hidden BF carriers.
# Hypothesis tables retain their established long-frame export.
#' @export
as.data.frame.BayesTools_table <- function(x, row.names = NULL, optional = FALSE, ...){

  if(inherits(x, "BayesTools_hypothesis_BF")){
    return(as.data.frame.BayesTools_hypothesis_BF(x, row.names, optional, ...))
  }
  .BayesTools_table_plain_frame(x, row.names = row.names, optional = optional, ...)
}

# The displayed columns of a table as a plain data frame at full precision.
.BayesTools_table_plain_frame <- function(x, row.names = NULL, optional = FALSE, ...){

  out <- as.data.frame.data.frame(x, row.names = row.names, optional = optional, ...)
  for(column in names(out)){
    if(inherits(out[[column]], "BayesTools_BF")) out[[column]] <- .BF_carrier_drop(out[[column]])
  }
  for(attribute in c("raw_log_BF", "numerical_diagnostics", "prior_numerical_diagnostics")){
    attr(out, attribute) <- NULL
  }
  out
}
