#' @title Hypothesis parser helpers
#'
#' @description Utilities for parsing point-null hypothesis references and
#' level references in the same syntax accepted by \code{\link{hypothesis_BF}}.
#'
#' @param hypothesis character vector with hypothesis statements or a validated
#'   `BayesTools_hypothesis_ast`.
#' @param allow_compound whether compound left-hand side expressions such as
#' \code{"theta + 0 = 0"} should be returned with \code{direct = FALSE}
#' instead of rejected.
#' @param text character vector with possible \code{parameter[level]}
#' references. Formula interaction names such as \code{factor:moderator[level]}
#' may be supplied without backticks.
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

  check_bool(allow_compound, "allow_compound", allow_NA = FALSE)
  if(inherits(hypothesis, "BayesTools_hypothesis_ast")){
    ast <- hypothesis
    hypothesis_text <- hypothesis_render(ast)
  }else{
    check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
    hypothesis_text <- hypothesis
    ast <- hypothesis_parse(hypothesis)
  }
  rows <- list()
  for(i in seq_along(ast$statements)){
    statement <- ast$statements[[i]]
    side_names <- if(isTRUE(statement$explicit)) c("left", "right") else "left"
    for(side_name in side_names){
      side <- statement[[side_name]]
      if(!side[["type"]] %in% c("point", "not_point")){
        next
      }
      row <- .bt_hypothesis_point_reference_row(
        hypothesis = hypothesis_text[[i]],
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
    node <- tryCatch(
      .bt_hypothesis_ast_node(.hypothesis_parse_expression(normalized)),
      error = function(e) NULL
    )
    symbol <- if(is.null(node)){
      NULL
    }else if(identical(node$type, "symbol")){
      node$name
    }else if(identical(node$type, "level_reference")){
      paste0(node$parameter, "[", node$level, "]")
    }else{
      NULL
    }
    direct <- !is.null(node) && identical(node$type, "level_reference")

    data.frame(
      input     = text_i,
      symbol    = if(is.null(symbol)) NA_character_ else symbol,
      parameter = if(direct) node$parameter else NA_character_,
      level     = if(direct) node$level else NA_character_,
      direct    = direct,
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


.hypothesis_parse_statement_spec <- function(hypothesis) {

  hypothesis <- .hypothesis_normalize_level_references(hypothesis)
  .hypothesis_reject_assignment_arrow(hypothesis)
  parts <- .hypothesis_split_vs(hypothesis)
  left  <- .parse_hypothesis_side(parts[[1L]])
  right <- if(length(parts) == 2L) {
    .parse_hypothesis_side(parts[[2L]])
  }else{
    .hypothesis_complement_side(left)
  }

  list(
    input    = hypothesis,
    left     = left,
    right    = right,
    explicit = length(parts) == 2L
  )
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
        "\\b([A-Za-z.][A-Za-z0-9._]*(?::[A-Za-z.][A-Za-z0-9._]*)*)\\s*\\[\\s*([^\\]\\[]+)\\s*\\]",
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

  relation <- if(.hypothesis_has_boolean(side)){
    NULL
  }else{
    .hypothesis_find_relation(side)
  }
  if(!is.null(relation) &&
     relation[["operator"]] %in% c("=", "==", "!=")){
    return(.hypothesis_parse_point_side(side, relation))
  }

  expr <- .hypothesis_parse_expression(side)
  if(.hypothesis_condition_has_equality(expr)){
    if(.hypothesis_region_has_negated_equality(expr)){
      stop(
        "Point equalities cannot be negated. Use an explicit point ",
        "hypothesis or its top-level '!=' complement.",
        call. = FALSE
      )
    }
    stop(
      "Equality constraints cannot be combined with '&' or '|', and point ",
      "equalities cannot be parenthesized inside a region expression.",
      call. = FALSE
    )
  }
  .hypothesis_validate_expression(expr, condition = TRUE)

  simple_expr <- .hypothesis_unwrap_parentheses(expr)
  simple_fun <- .hypothesis_call_name(simple_expr)
  if(!is.null(simple_fun) &&
     simple_fun %in% c("<", "<=", ">", ">=")){
    lhs_expr <- simple_expr[[2L]]
    rhs_expr <- simple_expr[[3L]]
    lhs_symbols <- .hypothesis_expression_symbols(lhs_expr)
    rhs_symbols <- .hypothesis_expression_symbols(rhs_expr)
    if(length(lhs_symbols) == 0L && length(rhs_symbols) > 0L){
      simple_fun <- switch(
        simple_fun,
        ">"  = "<",
        ">=" = "<=",
        "<"  = ">",
        "<=" = ">="
      )
      swap_expr <- lhs_expr
      lhs_expr  <- rhs_expr
      rhs_expr  <- swap_expr
      expr <- as.call(list(as.name(simple_fun), lhs_expr, rhs_expr))
      rhs_symbols <- character()
    }
    lhs <- .hypothesis_expression_text(lhs_expr)
    rhs <- .hypothesis_expression_text(rhs_expr)
    return(list(
      type      = "region",
      label     = .hypothesis_display_text(side),
      condition = .hypothesis_expression_text(expr),
      condition_expression = expr,
      expr      = lhs,
      expression = lhs_expr,
      value     = if(length(rhs_symbols) == 0L){
        .hypothesis_parse_number(rhs)
      }else{
        NULL
      },
      operator  = simple_fun,
      rhs       = rhs,
      rhs_expression = rhs_expr,
      simple    = TRUE
    ))
  }

  list(
    type      = "region",
    label     = .hypothesis_display_text(side),
    condition = side,
    condition_expression = expr,
    simple    = FALSE
  )
}


.bt_hypothesis_point_reference_row <- function(hypothesis, side_name, side){

  node <- side$expression
  direct <- node$type %in% c("symbol", "level_reference")
  symbol <- if(!direct){
    NA_character_
  }else if(identical(node$type, "level_reference")){
    paste0(node$parameter, "[", node$level, "]")
  }else{
    node$name
  }
  parameter <- if(!direct){
    NA_character_
  }else if(identical(node$type, "level_reference")){
    node$parameter
  }else{
    node$name
  }
  level <- if(identical(node$type, "level_reference")){
    node$level
  }else{
    NA_character_
  }

  data.frame(
    hypothesis = hypothesis,
    side = side_name,
    symbol = symbol,
    parameter = parameter,
    level = level,
    value = side$value,
    operator = if(identical(side$type, "point")) "=" else "!=",
    direct = direct,
    stringsAsFactors = FALSE
  )
}

.hypothesis_parse_point_side <- function(side, relation){

  lhs <- trimws(substr(side, 1L, relation[["start"]] - 1L))
  rhs <- trimws(substr(side, relation[["end"]] + 1L, nchar(side)))
  op <- relation[["operator"]]
  if(!nzchar(lhs) || !nzchar(rhs)){
    stop(
      "Point hypothesis relation must have both left and right sides.",
      call. = FALSE
    )
  }

  lhs_expression <- .hypothesis_parse_expression(lhs)
  .hypothesis_validate_expression(lhs_expression, condition = FALSE)
  rhs_expression <- .hypothesis_parse_expression(rhs)
  lhs_symbols <- .hypothesis_expression_symbols(lhs_expression)
  rhs_symbols <- .hypothesis_expression_symbols(rhs_expression)
  if(length(rhs_symbols) == 0L){
    expression <- lhs_expression
    value <- .hypothesis_parse_point_value(rhs)
  }else if(length(lhs_symbols) == 0L){
    .hypothesis_validate_expression(rhs_expression, condition = FALSE)
    expression <- rhs_expression
    value <- .hypothesis_parse_point_value(lhs)
  }else{
    .hypothesis_validate_expression(rhs_expression, condition = FALSE)
    expression <- as.call(list(
      as.name("-"), lhs_expression, rhs_expression
    ))
    value <- 0
  }
  list(
    type  = if(op == "!=") "not_point" else "point",
    label = paste(
      .hypothesis_display_text(lhs),
      if(op == "!=") "!=" else "=",
      .hypothesis_display_text(rhs)
    ),
    expr  = .hypothesis_expression_text(expression),
    expression = expression,
    value = value
  )
}

.hypothesis_unwrap_parentheses <- function(expr){

  while(identical(.hypothesis_call_name(expr), "(")){
    expr <- expr[[2L]]
  }
  expr
}

.hypothesis_call_name <- function(expr){

  if(!is.call(expr) || !is.name(expr[[1L]])){
    return(NULL)
  }
  as.character(expr[[1L]])
}

.hypothesis_expression_text <- function(expr){

  if(is.name(expr)){
    name <- .hypothesis_decode_escaped_constant(as.character(expr))
    if(name %in% .hypothesis_escaped_constant_names() ||
       !identical(make.names(name), name)){
      name <- gsub("`", "\\`", name, fixed = TRUE)
      return(paste0("`", name, "`"))
    }
    return(name)
  }
  expr <- .hypothesis_restore_escaped_constants(expr)
  paste(deparse(expr, width.cutoff = 500L), collapse = "")
}

.hypothesis_restore_escaped_constants <- function(expr){

  if(is.name(expr)){
    return(as.name(
      .hypothesis_decode_escaped_constant(as.character(expr))
    ))
  }
  if(!is.call(expr)){
    return(expr)
  }
  as.call(lapply(as.list(expr), .hypothesis_restore_escaped_constants))
}

.hypothesis_region_has_negated_equality <- function(expr, negated = FALSE){

  if(!is.call(expr)){
    return(FALSE)
  }
  fun <- .hypothesis_call_name(expr)
  if(identical(fun, "!")){
    return(.hypothesis_region_has_negated_equality(expr[[2L]], TRUE))
  }
  if(!is.null(fun) && fun %in% c("=", "==", "!=")){
    return(negated)
  }
  any(vapply(
    as.list(expr[-1L]),
    .hypothesis_region_has_negated_equality,
    logical(1),
    negated = negated
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
    condition_expression <- as.call(list(
      as.name(complement),
      side[["expression"]],
      side[["rhs_expression"]]
    ))
    return(list(
      type      = "region",
      label     = .hypothesis_display_text(condition),
      condition = condition,
      condition_expression = condition_expression,
      expr      = side[["expr"]],
      expression = side[["expression"]],
      value     = side[["value"]],
      operator  = complement,
      rhs       = side[["rhs"]],
      rhs_expression = side[["rhs_expression"]],
      simple    = TRUE
    ))
  }

  condition_expression <- as.call(list(
    as.name("!"),
    as.call(list(as.name("("), side[["condition_expression"]]))
  ))
  return(list(
    type       = "region",
    label      = paste0("not (", side[["label"]], ")"),
    condition  = side[["condition"]],
    condition_expression = condition_expression,
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
  .hypothesis_validate_expression(expr, condition = FALSE)

  value <- eval(expr, envir = .hypothesis_eval_parent())
  check_real(value, "hypothesis value", check_length = 1, allow_NA = FALSE)
  if(!is.finite(value)){
    stop("Hypothesis value must be finite.", call. = FALSE)
  }

  return(value)
}

.hypothesis_parse_point_value <- function(text){

  expr <- .hypothesis_parse_expression(text)
  numeric_literal <- function(x){
    (is.numeric(x) || is.integer(x)) && length(x) == 1L && is.finite(x)
  }
  valid_literal <- function(x){
    if(numeric_literal(x)){
      return(TRUE)
    }
    fun <- .hypothesis_call_name(x)
    if(!is.null(fun) && fun %in% c("+", "-") &&
       length(x) == 2L){
      return(numeric_literal(x[[2L]]))
    }
    FALSE
  }
  if(!valid_literal(expr)){
    stop(
      "The right side of a point hypothesis must be a numeric value written ",
      "as one finite literal with an optional unary sign.",
      call. = FALSE
    )
  }

  value <- eval(expr, envir = .hypothesis_eval_parent())
  if(!is.numeric(value) || length(value) != 1L || !is.finite(value)){
    stop("Point hypothesis value must be finite.", call. = FALSE)
  }
  as.numeric(value)
}

.hypothesis_escaped_constant_names <- function(){
  c("Inf", "NaN", "NA", "TRUE", "FALSE")
}

.hypothesis_escaped_constant_symbol <- function(name){
  paste0(".BayesTools_escaped_constant_", name)
}

.hypothesis_protect_escaped_constants <- function(text){

  for(name in .hypothesis_escaped_constant_names()){
    text <- gsub(
      paste0("`", name, "`"),
      .hypothesis_escaped_constant_symbol(name),
      text,
      fixed = TRUE
    )
  }
  text
}

.hypothesis_decode_escaped_constant <- function(name){

  symbols <- vapply(
    .hypothesis_escaped_constant_names(),
    .hypothesis_escaped_constant_symbol,
    character(1)
  )
  index <- match(name, symbols)
  if(is.na(index)){
    return(name)
  }
  .hypothesis_escaped_constant_names()[[index]]
}


.hypothesis_parse_expression <- function(text) {

  if(is.name(text) || is.call(text) || is.numeric(text) || is.integer(text) ||
     is.logical(text)){
    return(text)
  }
  if(!is.character(text) || length(text) != 1L || is.na(text)){
    stop("Hypothesis expression must be scalar text or a parsed language object.",
         call. = FALSE)
  }

  protected_text <- .hypothesis_protect_escaped_constants(text)
  parsed <- tryCatch(parse(text = protected_text, keep.source = FALSE),
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

  if(condition){
    .hypothesis_validate_region(expr)
  }else{
    .hypothesis_validate_arithmetic(expr)
  }

  return(invisible(TRUE))
}

.hypothesis_validate_arithmetic <- function(expr){

  if(is.numeric(expr) || is.integer(expr)){
    if(length(expr) != 1L || !is.finite(expr)){
      stop(
        "Hypothesis arithmetic requires finite numeric literals.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }
  if(is.logical(expr)){
    literal <- paste(deparse(expr, width.cutoff = 500L), collapse = "")
    stop(
      "Unescaped reserved literal '", literal,
      "' is not a hypothesis parameter or finite numeric value.",
      call. = FALSE
    )
  }
  if(is.name(expr)){
    return(invisible(TRUE))
  }
  if(!is.call(expr)){
    stop("Unsupported hypothesis arithmetic expression.", call. = FALSE)
  }

  call_head <- expr[[1L]]
  if(!is.name(call_head)){
    call_label <- paste(deparse(call_head, width.cutoff = 500L), collapse = "")
    stop("Unsupported hypothesis expression call '", call_label, "'.",
         call. = FALSE)
  }
  fun <- as.character(call_head)
  allowed_functions <- c("abs", "exp", "log", "sqrt", "plogis", "qlogis")
  if(!fun %in% c("(", "+", "-", "*", "/", "^", allowed_functions)){
    stop("Unsupported hypothesis expression operator or function '", fun, "'.",
         call. = FALSE)
  }

  n_arguments <- length(expr) - 1L
  valid_arity <- switch(
    fun,
    "(" = n_arguments == 1L,
    "+" = n_arguments %in% c(1L, 2L),
    "-" = n_arguments %in% c(1L, 2L),
    "*" = n_arguments == 2L,
    "/" = n_arguments == 2L,
    "^" = n_arguments == 2L,
    n_arguments == 1L
  )
  if(!valid_arity){
    if(n_arguments == 0L){
      stop(
        "Hypothesis expression call '", fun,
        "' requires at least one argument.",
        call. = FALSE
      )
    }
    stop(
      "Hypothesis expression call '", fun,
      "' has an unsupported number of arguments.",
      call. = FALSE
    )
  }
  for(i in seq.int(2L, length(expr))){
    .hypothesis_validate_arithmetic(expr[[i]])
  }

  if(length(.hypothesis_expression_symbols(expr)) == 0L){
    value <- tryCatch(
      suppressWarnings(eval(expr, envir = .hypothesis_eval_parent())),
      error = function(e) NULL
    )
    if(is.null(value) || !is.numeric(value) || length(value) != 1L ||
       !is.finite(value)){
      stop(
        "Constant hypothesis arithmetic must evaluate to one finite numeric ",
        "value.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.hypothesis_validate_region <- function(expr){

  if(!is.call(expr)){
    stop(
      "A region hypothesis must be a comparison, optionally parenthesized, ",
      "negated, or combined with '&' or '|'.",
      call. = FALSE
    )
  }
  call_head <- expr[[1L]]
  if(!is.name(call_head)){
    stop("Unsupported region hypothesis expression.", call. = FALSE)
  }
  fun <- as.character(call_head)
  n_arguments <- length(expr) - 1L

  if(fun == "("){
    if(n_arguments != 1L){
      stop("Parenthesized regions require exactly one expression.",
           call. = FALSE)
    }
    .hypothesis_validate_region(expr[[2L]])
    return(invisible(TRUE))
  }
  if(fun == "!"){
    if(n_arguments != 1L){
      stop("Region negation requires exactly one region.", call. = FALSE)
    }
    .hypothesis_validate_region(expr[[2L]])
    return(invisible(TRUE))
  }
  if(fun %in% c("&", "|")){
    if(n_arguments != 2L){
      stop(
        "Region boolean operators require exactly two region expressions.",
        call. = FALSE
      )
    }
    .hypothesis_validate_region(expr[[2L]])
    .hypothesis_validate_region(expr[[3L]])
    return(invisible(TRUE))
  }
  if(fun %in% c("<", "<=", ">", ">=")){
    if(n_arguments != 2L){
      stop("Region relations require exactly two arithmetic expressions.",
           call. = FALSE)
    }
    .hypothesis_validate_arithmetic(expr[[2L]])
    .hypothesis_validate_arithmetic(expr[[3L]])
    return(invisible(TRUE))
  }

  stop(
    "Unsupported region operator '", fun,
    "'. Use '<', '<=', '>', or '>=' and combine regions with '&', '|', or '!'.",
    call. = FALSE
  )
}


.hypothesis_condition_has_equality <- function(expr) {

  if(!is.call(expr)){
    return(FALSE)
  }
  fun <- .hypothesis_call_name(expr)
  if(!is.null(fun) && fun %in% c("=", "==", "!=")){
    return(TRUE)
  }

  return(any(vapply(as.list(expr[-1L]), .hypothesis_condition_has_equality,
                    logical(1))))
}


.hypothesis_expression_symbols <- function(expr) {

  collect <- function(node){
    if(is.name(node)){
      return(.hypothesis_decode_escaped_constant(as.character(node)))
    }
    if(!is.call(node) || length(node) == 1L){
      return(character())
    }
    unlist(lapply(as.list(node[-1L]), collect), use.names = FALSE)
  }

  unique(collect(expr))
}


.hypothesis_statement_symbols <- function(statement) {

  symbols <- character()
  for(side_name in c("left", "right")){
    nodes <- .bt_hypothesis_symbol_nodes(statement[[side_name]]$expression)
    symbols <- c(symbols, vapply(nodes, function(node) {
      if(identical(node$type, "level_reference")) {
        paste0(node$parameter, "[", node$level, "]")
      } else {
        node$name
      }
    }, character(1)))
  }

  return(unique(symbols))
}


.hypothesis_all_symbols <- function(statements) {

  unique(unlist(
    lapply(statements, .hypothesis_statement_symbols),
    use.names = FALSE
  ))
}


.hypothesis_single_symbol <- function(statements) {

  symbols <- .hypothesis_all_symbols(statements)
  if(length(symbols) == 1L){
    return(symbols)
  }

  return(NULL)
}
