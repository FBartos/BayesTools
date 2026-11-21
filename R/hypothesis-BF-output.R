

.hypothesis_BF_evidence_columns <- c(
  "prior_density", "posterior_density",
  "prior_prob_alt", "posterior_prob_alt", "prior_prob_null", "posterior_prob_null"
)

.hypothesis_BF_row <- function(quantity, result) {

  row <- data.frame(
    Alternative  = result[["alternative"]],
    Null         = result[["null"]],
    BF           = as.numeric(result[["BF"]]),
    BF_error     = as.numeric(result[["BF_error"]]),
    method       = result[["method"]],
    warning      = .hypothesis_collapse_warning(result[["warning"]]),
    row.names    = quantity[["label"]],
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
  for(column in .hypothesis_BF_evidence_columns){
    value <- result[[column]]
    row[[column]] <- if(is.null(value)) NA_real_ else as.numeric(value)
  }
  row
}


.hypothesis_BF_row_names <- function(labels, statements) {

  # Rows of several statements on one quantity carry the statement number
  # ("mu (1)", "mu (2)"); rbind()'s "mu1" reads like another parameter, and
  # "mu.1" would collide with the names `[.data.frame` gives duplicated rows,
  # which row-name based subsetting of the table relies on.
  repeated <- labels %in% labels[duplicated(labels)]
  labels[repeated] <- paste0(labels[repeated], " (", statements[repeated], ")")

  make.unique(labels, sep = " ")
}


# Resolve the schema once per request, before splitting a fit into groups.
# Even an unavailable ordinate keeps its columns when a point was requested.
.hypothesis_BF_output_columns <- function(hypothesis, include_densities,
                                          include_probabilities, include_method) {

  check_bool(include_densities, "include_densities", allow_NA = FALSE)
  check_bool(include_probabilities, "include_probabilities", allow_NA = FALSE)
  check_bool(include_method, "include_method", allow_NA = FALSE)

  columns <- c("Alternative", "Null", "BF", "BF_error")
  types <- unlist(lapply(hypothesis$statements, function(statement) {
    c(statement$left$type, statement$right$type)
  }), use.names = FALSE)
  if(include_densities && "point" %in% types){
    columns <- c(columns, "prior_density", "posterior_density")
  }
  if(include_probabilities && any(types %in% c("region", "unconstrained"))){
    columns <- c(columns, "prior_prob_alt", "posterior_prob_alt", "prior_prob_null", "posterior_prob_null")
  }
  if(include_method) columns <- c(columns, "method")
  columns
}


.hypothesis_BF_table_types <- function(columns) {

  type <- c(
    Alternative = "hypothesis_label",
    Null        = "hypothesis_label",
    BF          = "BF",
    BF_error    = "BF_error",
    method      = "string"
  )
  type[.hypothesis_BF_evidence_columns] <- "estimate"

  unname(type[columns])
}


.hypothesis_BF_table_warnings <- function(out) {

  warnings <- out[["warning"]]
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if(length(warnings) == 0L){
    return(NULL)
  }

  names(warnings) <- rownames(out)[!is.na(out[["warning"]]) &
    nzchar(out[["warning"]])]

  warnings
}


.hypothesis_trapz <- function(x, y) {

  if(length(x) < 2L){
    return(0)
  }

  sum(diff(x) * (y[-1L] + y[-length(y)]) / 2)
}


.hypothesis_number_label <- function(x) {

  # Keep familiar short labels when they preserve the exact parsed value; a
  # value whose fixed notation is long (1e-300, 1e22) uses the shortest
  # scientific notation that preserves it.
  x <- as.numeric(x)
  fixed <- NULL
  for(digits in 7:17){
    label <- format(x, digits = digits, trim = TRUE, scientific = FALSE,
                    decimal.mark = ".")
    if(identical(.hypothesis_read_number(label), x)){
      fixed <- label
      break
    }
  }
  if(!is.null(fixed) && nchar(fixed) <= 20L){
    return(fixed)
  }
  for(digits in 1:17){
    scientific <- format(x, digits = digits, trim = TRUE, scientific = TRUE,
                         decimal.mark = ".")
    if(identical(.hypothesis_read_number(scientific), x)){
      return(scientific)
    }
  }
  if(!is.null(fixed)){
    return(fixed)
  }
  .hypothesis_exact_number_label(x)
}

# Hypothesis text is read by R's parser, which reads a number as as.numeric()
# does, so a label must read back with as.numeric() on this platform.
.hypothesis_read_number <- function(text) {

  as.numeric(text)
}

# R reads decimal text by accumulating its digits and the power of ten in long
# double. Where long double is double (Apple silicon), some doubles have no
# decimal text of up to 17 significant digits that reads back as them; their
# hexadecimal literal reads back exactly on every platform. R cannot read the
# hexadecimal literal of a subnormal value (its power of two overflows), so
# such a value, a representability limit rather than a realistic hypothesis,
# is refused.
.hypothesis_exact_number_label <- function(x) {

  if(x != 0 && abs(x) < .Machine$double.xmin){
    stop(errorCondition(
      paste0("The hypothesis value ", sprintf("%.17g", x), " is subnormal, and ",
             "no text of it reads back as the same double on this platform. ",
             "Restate the hypothesis with a value that is 0 or at least ",
             "'.Machine$double.xmin' in magnitude."),
      class = "BayesTools_hypothesis_statement", call = NULL
    ))
  }
  sprintf("%a", x)
}


.hypothesis_collapse_warning <- function(warning) {

  warning <- unlist(warning, use.names = FALSE)
  warning <- warning[!is.na(warning) & nzchar(warning)]
  if(length(warning) == 0L){
    return(NA_character_)
  }

  paste(unique(warning), collapse = " ")
}
