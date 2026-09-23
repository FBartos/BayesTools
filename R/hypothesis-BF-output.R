

.hypothesis_BF_row <- function(quantity, result) {

  data.frame(
    Alternative  = result[["alternative"]],
    Null         = result[["null"]],
    BF           = result[["BF"]],
    BF_error     = result[["BF_error"]],
    prior        = result[["prior"]],
    posterior    = result[["posterior"]],
    method       = result[["method"]],
    warning      = .hypothesis_collapse_warning(result[["warning"]]),
    row.names    = quantity[["label"]],
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
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


.hypothesis_BF_output_columns <- function(columns) {

  default <- c("Alternative", "Null", "BF", "BF_error")
  extra   <- c("prior", "posterior", "method")

  if(any(columns == "all")){
    if(length(columns) > 1L){
      stop("'columns = \"all\"' cannot be combined with other column names.",
           call. = FALSE)
    }
    return(c(default, extra))
  }
  if(any(columns == "default")){
    columns <- setdiff(columns, "default")
    if(length(columns) == 0L){
      return(default)
    }
  }

  aliases <- c(
    Alternative             = "Alternative",
    alternative             = "Alternative",
    Null                    = "Null",
    null                    = "Null",
    BF                      = "BF",
    BF_error                = "BF_error",
    `error%(BF)`            = "BF_error",
    prior                   = "prior",
    Prior                   = "prior",
    posterior               = "posterior",
    Posterior               = "posterior",
    method                  = "method",
    Method                  = "method",
    computation_method      = "method",
    `computation method`    = "method"
  )
  mapped <- aliases[columns]
  if(any(is.na(mapped))){
    stop("Unknown 'columns' value: ",
         paste0("'", columns[is.na(mapped)], "'", collapse = ", "),
         ".", call. = FALSE)
  }

  unique(c(default, unname(mapped)))
}


.hypothesis_BF_table_types <- function(columns) {

  type <- c(
    Alternative = "hypothesis_label",
    Null        = "hypothesis_label",
    BF          = "BF",
    BF_error    = "BF_error",
    prior       = "estimate",
    posterior   = "estimate",
    method      = "string"
  )

  unname(type[columns])
}


.hypothesis_BF_table_footnotes <- function(columns) {

  if(any(c("prior", "posterior") %in% columns)){
    return(paste0(
      "Note: 'prior' and 'posterior' are diagnostic values, not always ",
      "probabilities. Point tests report density heights, region tests report ",
      "odds, and transitive point-vs-region tests report NA."
    ))
  }

  return(NULL)
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

  # Keep familiar short labels when they preserve the exact parsed value.
  for(digits in 7:17){
    label <- format(x, digits = digits, trim = TRUE, scientific = FALSE,
                    decimal.mark = ".")
    if(identical(as.numeric(label), as.numeric(x))){
      return(label)
    }
  }
  label
}


.hypothesis_collapse_warning <- function(warning) {

  warning <- unlist(warning, use.names = FALSE)
  warning <- warning[!is.na(warning) & nzchar(warning)]
  if(length(warning) == 0L){
    return(NA_character_)
  }

  paste(unique(warning), collapse = " ")
}
