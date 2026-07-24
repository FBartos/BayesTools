#' @title Hypothesis parser helpers
#'
#' @description Utilities for parsing point-null hypothesis references and
#' level references in the same syntax accepted by \code{\link{hypothesis_BF}}.
#'
#' @param hypothesis character vector with hypothesis statements.
#' @param allow_compound whether compound left-hand side expressions such as
#' \code{"theta + 0 = 0"} should be returned with \code{direct = FALSE}
#' instead of rejected.
#' @param text character vector with possible \code{parameter[level]}
#' references.
#'
#' @return \code{hypothesis_parse_point_reference()} returns a data frame with
#' columns \code{hypothesis}, \code{side}, \code{symbol}, \code{parameter},
#' \code{level}, \code{value}, \code{operator}, and \code{direct}.
#' \code{hypothesis_parse_level_reference()} returns a data frame with columns
#' \code{input}, \code{symbol}, \code{parameter}, \code{level}, and
#' \code{direct}. \code{hypothesis_normalize_level_references()} returns a
#' character vector.
#'
#' @export
hypothesis_parse_point_reference <- function(hypothesis,
                                             allow_compound = TRUE){

  check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  check_bool(allow_compound, "allow_compound", allow_NA = FALSE)

  rows <- list()
  for(i in seq_along(hypothesis)){
    parsed <- .parse_hypothesis_BF(hypothesis[[i]])
    side_names <- if(isTRUE(parsed[["explicit"]])) c("left", "right") else "left"
    for(side_name in side_names){
      side <- parsed[[side_name]]
      if(!side[["type"]] %in% c("point", "not_point")){
        next
      }
      row <- .hypothesis_point_reference_row(
        hypothesis = hypothesis[[i]],
        side_name  = side_name,
        side       = side
      )
      if(!isTRUE(allow_compound) && !isTRUE(row[["direct"]])){
        stop(
          "Point hypothesis side '", side[["label"]],
          "' is not a direct parameter or level reference.",
          call. = FALSE
        )
      }
      rows[[length(rows) + 1L]] <- row
    }
  }

  if(length(rows) == 0L){
    return(.hypothesis_empty_point_reference_table())
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  return(out)
}


#' @rdname hypothesis_parse_point_reference
#' @export
hypothesis_parse_level_reference <- function(text){

  check_char(text, "text", check_length = 0, allow_NA = FALSE)

  rows <- lapply(text, function(text_i){
    normalized <- .hypothesis_normalize_level_references(text_i)
    symbol <- tryCatch(
      .hypothesis_direct_symbol(normalized),
      error = function(e) NULL
    )
    ref <- .hypothesis_parse_level_symbol(symbol)

    data.frame(
      input     = text_i,
      symbol    = if(is.null(symbol)) NA_character_ else symbol,
      parameter = ref[["parameter"]],
      level     = ref[["level"]],
      direct    = !is.null(symbol) && isTRUE(ref[["direct"]]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  return(out)
}


#' @rdname hypothesis_parse_point_reference
#' @export
hypothesis_normalize_level_references <- function(text){

  check_char(text, "text", check_length = 0, allow_NA = FALSE)

  return(unname(vapply(
    text,
    .hypothesis_normalize_level_references,
    character(1)
  )))
}


.hypothesis_empty_point_reference_table <- function(){

  data.frame(
    hypothesis = character(),
    side       = character(),
    symbol     = character(),
    parameter  = character(),
    level      = character(),
    value      = numeric(),
    operator   = character(),
    direct     = logical(),
    stringsAsFactors = FALSE
  )
}


.hypothesis_point_reference_row <- function(hypothesis, side_name, side){

  symbol <- .hypothesis_direct_symbol(side[["expr"]])
  ref <- .hypothesis_parse_level_symbol(symbol)
  direct <- !is.null(symbol)
  parameter <- if(isTRUE(ref[["direct"]])){
    ref[["parameter"]]
  }else if(direct){
    symbol
  }else{
    NA_character_
  }

  data.frame(
    hypothesis = hypothesis,
    side       = side_name,
    symbol     = if(is.null(symbol)) NA_character_ else symbol,
    parameter  = parameter,
    level      = ref[["level"]],
    value      = side[["value"]],
    operator   = if(identical(side[["type"]], "point")) "=" else "!=",
    direct     = direct,
    stringsAsFactors = FALSE
  )
}


.hypothesis_parse_level_symbol <- function(symbol){

  empty <- list(
    parameter = NA_character_,
    level     = NA_character_,
    direct    = FALSE
  )
  if(is.null(symbol) || !nzchar(symbol)){
    return(empty)
  }

  ref <- regexec("^([^\\[]+)\\[([^\\]]+)\\]$", symbol, perl = TRUE)
  match <- regmatches(symbol, ref)[[1L]]
  if(length(match) != 3L){
    return(empty)
  }

  return(list(
    parameter = match[[2L]],
    level     = trimws(match[[3L]]),
    direct    = TRUE
  ))
}


.parse_hypothesis_BF <- function(hypothesis) {

  hypothesis <- .hypothesis_normalize_level_references(hypothesis)
  .hypothesis_reject_assignment_arrow(hypothesis)
  parts <- .hypothesis_split_vs(hypothesis)
  left  <- .parse_hypothesis_side(parts[[1L]])
  right <- if(length(parts) == 2L) {
    .parse_hypothesis_side(parts[[2L]])
  }else{
    .hypothesis_complement_side(left)
  }

  out <- list(
    input    = hypothesis,
    left     = left,
    right    = right,
    explicit = length(parts) == 2L
  )
  class(out) <- "BayesTools_hypothesis_BF_parsed"

  return(out)
}


.hypothesis_reject_assignment_arrow <- function(text) {

  chars <- strsplit(text, "", fixed = TRUE)[[1L]]
  in_backtick <- FALSE

  i <- 1L
  while(i <= length(chars)){
    ch <- chars[[i]]
    if(ch == "`"){
      in_backtick <- !in_backtick
    }else if(!in_backtick && ch == "<" && i < length(chars) &&
             chars[[i + 1L]] == "-"){
      stop("Hypothesis uses '<-'. Use '=' or '==' for equality hypotheses.",
           call. = FALSE)
    }
    i <- i + 1L
  }

  invisible(TRUE)
}


.hypothesis_normalize_level_references <- function(text) {

  ends_in_backtick <- endsWith(text, "`")
  pieces <- strsplit(text, "`", fixed = TRUE)[[1L]]
  if(ends_in_backtick){
    pieces <- c(pieces, "")
  }
  if(length(pieces) == 0L){
    return(text)
  }

  for(i in seq_along(pieces)){
    if(i %% 2L == 1L){
      pieces[[i]] <- gsub(
        "\\b([A-Za-z.][A-Za-z0-9._]*)\\s*\\[\\s*([^\\]\\[]+)\\s*\\]",
        "`\\1[\\2]`",
        pieces[[i]],
        perl = TRUE
      )
      pieces[[i]] <- gsub(
        "`([^`\\[]+)\\[\\s*([^\\]\\[]*\\S)\\s*\\]`",
        "`\\1[\\2]`",
        pieces[[i]],
        perl = TRUE
      )
    }
  }

  paste(pieces, collapse = "`")
}


.hypothesis_split_vs <- function(hypothesis) {

  chars <- strsplit(hypothesis, "", fixed = TRUE)[[1L]]
  in_backtick <- FALSE
  depth <- 0L
  found <- integer()

  i <- 1L
  while(i <= length(chars)){
    ch <- chars[[i]]
    if(ch == "`"){
      in_backtick <- !in_backtick
      i <- i + 1L
      next
    }
    if(!in_backtick){
      if(ch == "("){
        depth <- depth + 1L
      }else if(ch == ")"){
        depth <- max(0L, depth - 1L)
      }else if(depth == 0L && i > 1L && i < length(chars) &&
               tolower(ch) == "v" && tolower(chars[[i + 1L]]) == "s" &&
               grepl("\\s", chars[[i - 1L]]) &&
               i + 2L <= length(chars) && grepl("\\s", chars[[i + 2L]])){
        found <- c(found, i)
        i <- i + 2L
        next
      }
    }
    i <- i + 1L
  }

  if(length(found) == 0L){
    parts <- hypothesis
  }else if(length(found) == 1L){
    parts <- c(
      substr(hypothesis, 1L, found[[1L]] - 1L),
      substr(hypothesis, found[[1L]] + 2L, nchar(hypothesis))
    )
  }else{
    parts <- character()
  }
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]

  if(length(parts) == 0L || length(parts) > 2L){
    stop("Hypothesis must contain one statement or one explicit 'vs' comparison.",
         call. = FALSE)
  }

  return(parts)
}


.parse_hypothesis_side <- function(side) {

  side <- trimws(side)
  if(!nzchar(side)){
    stop("Empty hypothesis side.", call. = FALSE)
  }

  if(.hypothesis_has_boolean(side)){
    expr <- .hypothesis_parse_expression(side)
    if(.hypothesis_condition_has_equality(expr)){
      stop("Equality constraints cannot be combined with '&' or '|'.",
           call. = FALSE)
    }
    .hypothesis_validate_expression(expr, condition = TRUE)
    return(list(
      type      = "region",
      label     = .hypothesis_display_text(side),
      condition = side,
      simple    = FALSE
    ))
  }

  relation <- .hypothesis_find_relation(side)
  if(is.null(relation)){
    stop("Hypothesis side must contain one relation operator.", call. = FALSE)
  }

  lhs <- trimws(substr(side, 1L, relation[["start"]] - 1L))
  rhs <- trimws(substr(side, relation[["end"]] + 1L, nchar(side)))
  op  <- relation[["operator"]]

  if(!nzchar(lhs) || !nzchar(rhs)){
    stop("Hypothesis relation must have both left and right sides.",
         call. = FALSE)
  }

  if(op %in% c("=", "==")){
    .hypothesis_validate_expression(.hypothesis_parse_expression(lhs),
                                    condition = FALSE)
    value <- .hypothesis_parse_number(rhs)
    return(list(
      type  = "point",
      label = paste(.hypothesis_display_text(lhs), "=", .hypothesis_display_text(rhs)),
      expr  = lhs,
      value = value
    ))
  }
  if(op == "!="){
    .hypothesis_validate_expression(.hypothesis_parse_expression(lhs),
                                    condition = FALSE)
    value <- .hypothesis_parse_number(rhs)
    return(list(
      type  = "not_point",
      label = paste(.hypothesis_display_text(lhs), "!=", .hypothesis_display_text(rhs)),
      expr  = lhs,
      value = value
    ))
  }

  .hypothesis_validate_expression(.hypothesis_parse_expression(lhs),
                                  condition = FALSE)
  .hypothesis_validate_expression(.hypothesis_parse_expression(rhs),
                                  condition = FALSE)
  condition <- paste(lhs, op, rhs)
  .hypothesis_validate_expression(.hypothesis_parse_expression(condition),
                                  condition = TRUE)
  rhs_symbols <- .hypothesis_expression_symbols(.hypothesis_parse_expression(rhs))
  return(list(
    type      = "region",
    label     = .hypothesis_display_text(condition),
    condition = condition,
    expr      = lhs,
    value     = if(length(rhs_symbols) == 0L) .hypothesis_parse_number(rhs) else NULL,
    operator  = op,
    rhs       = rhs,
    simple    = TRUE
  ))
}


.hypothesis_complement_side <- function(side) {

  if(identical(side[["type"]], "point")){
    out <- side
    out[["type"]]  <- "not_point"
    out[["label"]] <- paste(
      .hypothesis_display_text(side[["expr"]]),
      "!=",
      .hypothesis_number_label(side[["value"]])
    )
    return(out)
  }

  if(identical(side[["type"]], "not_point")){
    out <- side
    out[["type"]]  <- "point"
    out[["label"]] <- paste(
      .hypothesis_display_text(side[["expr"]]),
      "=",
      .hypothesis_number_label(side[["value"]])
    )
    return(out)
  }

  if(!is.null(side[["simple"]]) && isTRUE(side[["simple"]])){
    complement <- switch(
      side[["operator"]],
      ">"  = "<=",
      ">=" = "<",
      "<"  = ">=",
      "<=" = ">"
    )
    condition <- paste(side[["expr"]], complement, side[["rhs"]])
    return(list(
      type      = "region",
      label     = .hypothesis_display_text(condition),
      condition = condition,
      expr      = side[["expr"]],
      value     = side[["value"]],
      operator  = complement,
      rhs       = side[["rhs"]],
      simple    = TRUE
    ))
  }

  return(list(
    type       = "region",
    label      = paste0("not (", side[["label"]], ")"),
    condition  = side[["condition"]],
    complement = TRUE,
    simple     = FALSE
  ))
}


.hypothesis_display_text <- function(text) {

  gsub("`([^`]+)`", "\\1", text, perl = TRUE)
}


.hypothesis_find_relation <- function(text) {

  chars       <- strsplit(text, "", fixed = TRUE)[[1L]]
  in_backtick <- FALSE
  depth       <- 0L
  found       <- list()
  two_char    <- c("<=", ">=", "==", "!=")
  one_char    <- c("<", ">", "=")

  i <- 1L
  while(i <= length(chars)){
    ch <- chars[[i]]
    if(ch == "`"){
      in_backtick <- !in_backtick
      i <- i + 1L
      next
    }
    if(!in_backtick){
      if(ch == "("){
        depth <- depth + 1L
      }else if(ch == ")"){
        depth <- max(0L, depth - 1L)
      }else if(depth == 0L){
        op <- NULL
        if(i < length(chars)){
          candidate <- paste0(chars[[i]], chars[[i + 1L]])
          if(candidate %in% two_char){
            op <- candidate
          }
        }
        if(is.null(op) && ch %in% one_char){
          op <- ch
        }
        if(!is.null(op)){
          found[[length(found) + 1L]] <- list(
            operator = op,
            start    = i,
            end      = i + nchar(op) - 1L
          )
          i <- i + nchar(op)
          next
        }
      }
    }
    i <- i + 1L
  }

  if(length(found) == 0L){
    return(NULL)
  }
  if(length(found) > 1L){
    stop("Hypothesis side must contain exactly one top-level relation operator.",
         call. = FALSE)
  }

  return(found[[1L]])
}


.hypothesis_has_boolean <- function(text) {

  chars       <- strsplit(text, "", fixed = TRUE)[[1L]]
  in_backtick <- FALSE
  for(ch in chars){
    if(ch == "`"){
      in_backtick <- !in_backtick
    }else if(!in_backtick && ch %in% c("&", "|")){
      return(TRUE)
    }
  }

  return(FALSE)
}


.hypothesis_parse_number <- function(text) {

  expr <- .hypothesis_parse_expression(text)
  names <- .hypothesis_expression_symbols(expr)
  if(length(names) > 0L){
    stop("Right side of a point/one-sided hypothesis must be a numeric value.",
         call. = FALSE)
  }

  value <- eval(expr, envir = .hypothesis_eval_parent())
  check_real(value, "hypothesis value", check_length = 1, allow_NA = FALSE)
  if(!is.finite(value)){
    stop("Hypothesis value must be finite.", call. = FALSE)
  }

  return(value)
}


.hypothesis_parse_expression <- function(text) {

  parsed <- tryCatch(parse(text = text, keep.source = FALSE),
                     error = function(e)e)
  if(inherits(parsed, "error") || length(parsed) != 1L){
    stop("Could not parse hypothesis expression '", text, "'.", call. = FALSE)
  }

  return(parsed[[1L]])
}


.hypothesis_eval_parent <- local({
  env <- new.env(parent = baseenv())
  env[["plogis"]] <- stats::plogis
  env[["qlogis"]] <- stats::qlogis
  function(){
    env
  }
})


.hypothesis_validate_expression <- function(expr, condition) {

  if(is.numeric(expr) || is.integer(expr)){
    return(invisible(TRUE))
  }
  if(is.logical(expr)){
    return(invisible(TRUE))
  }
  if(is.name(expr)){
    name <- as.character(expr)
    if(name %in% c("Inf", "NaN", "NA", "TRUE", "FALSE")){
      return(invisible(TRUE))
    }
    return(invisible(TRUE))
  }
  if(!is.call(expr)){
    stop("Unsupported hypothesis expression.", call. = FALSE)
  }

  fun <- as.character(expr[[1L]])
  allowed_arithmetic <- c("(", "+", "-", "*", "/", "^")
  allowed_functions  <- c("abs", "exp", "log", "sqrt", "plogis", "qlogis")
  allowed_condition  <- c("<", "<=", ">", ">=", "&", "|", "!")
  allowed <- c(allowed_arithmetic, allowed_functions)
  if(condition){
    allowed <- c(allowed, allowed_condition)
  }

  if(!fun %in% allowed){
    stop("Unsupported hypothesis expression operator or function '", fun, "'.",
         call. = FALSE)
  }
  for(i in seq.int(2L, length(expr))){
    .hypothesis_validate_expression(expr[[i]], condition = condition)
  }

  return(invisible(TRUE))
}


.hypothesis_condition_has_equality <- function(expr) {

  if(!is.call(expr)){
    return(FALSE)
  }
  fun <- as.character(expr[[1L]])
  if(fun %in% c("=", "==", "!=")){
    return(TRUE)
  }

  return(any(vapply(as.list(expr[-1L]), .hypothesis_condition_has_equality,
                    logical(1))))
}


.hypothesis_expression_symbols <- function(expr) {

  names <- all.names(expr, functions = TRUE, unique = TRUE)
  blocked <- c(
    "(", "+", "-", "*", "/", "^", "<", "<=", ">", ">=", "=", "==", "!=",
    "&", "|", "!", "abs", "exp", "log", "sqrt", "plogis", "qlogis",
    "Inf", "NaN", "NA", "TRUE", "FALSE"
  )

  return(setdiff(names, blocked))
}


.hypothesis_parsed_symbols <- function(parsed) {

  symbols <- character()
  for(side_name in c("left", "right")){
    side <- parsed[[side_name]]
    if(!is.null(side[["expr"]])){
      symbols <- c(symbols, .hypothesis_expression_symbols(
        .hypothesis_parse_expression(side[["expr"]])
      ))
    }
    if(!is.null(side[["condition"]])){
      symbols <- c(symbols, .hypothesis_expression_symbols(
        .hypothesis_parse_expression(side[["condition"]])
      ))
    }
  }

  return(unique(symbols))
}


.hypothesis_all_symbols <- function(parsed) {

  unique(unlist(lapply(parsed, .hypothesis_parsed_symbols), use.names = FALSE))
}


.hypothesis_single_symbol <- function(parsed) {

  symbols <- .hypothesis_all_symbols(parsed)
  if(length(symbols) == 1L){
    return(symbols)
  }

  return(NULL)
}
