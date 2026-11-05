

.hypothesis_BF_row <- function(quantity, result) {

  data.frame(
    Alternative  = result[["alternative"]],
    Null         = result[["null"]],
    BF           = as.numeric(result[["BF"]]),
    BF_error     = as.numeric(result[["BF_error"]]),
    prior        = as.numeric(result[["prior"]]),
    posterior    = as.numeric(result[["posterior"]]),
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

  unknown <- setdiff(columns, c(default, extra))
  if(length(unknown) > 0L){
    stop("Unknown 'columns' value: ",
         paste0("'", unknown, "'", collapse = ", "),
         ". Columns are selected by their names: ",
         paste0("'", c(default, extra), "'", collapse = ", "),
         ", or 'default' or 'all'.", call. = FALSE)
  }

  unique(c(default, columns))
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
# hexadecimal literal reads back exactly on every platform. R scales a
# hexadecimal literal by a power of two that overflows below 2^-1022, so a
# subnormal value is written as its multiple of 2^-1074 on the 2^-1022 scale.
.hypothesis_exact_number_label <- function(x) {

  if(x != 0 && abs(x) < .Machine$double.xmin){
    units  <- abs(x) * 2^1022 * 2^52
    digits <- (units %/% 16^(12:0)) %% 16
    return(paste0(if(x < 0) "-", "0x0.",
                  paste(c(0:9, letters[1:6])[digits + 1], collapse = ""), "p-1022"))
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
