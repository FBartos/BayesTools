# Public versioned hypothesis syntax tree.

.bt_hypothesis_ast_version <- 1L
.bt_hypothesis_resolution_version <- 1L

#' Versioned hypothesis syntax trees
#'
#' @description
#' `hypothesis_parse()` parses the restricted language used by
#' [hypothesis_BF()] into a validated, serializable syntax tree.
#' `hypothesis_render()` renders that tree with token-safe quoting.
#' `hypothesis_symbols()` returns exact parameter roots or an occurrence table.
#'
#' `hypothesis_rewrite()` replaces only exact symbol roots. Function names,
#' literals, level labels, and longer identifiers are never edited.
#' `hypothesis_resolve()` resolves every occurrence through
#' [parameter_catalog_resolve()] without accessing posterior draws.
#'
#' @param hypothesis character vector of hypothesis statements.
#' @param ast a `BayesTools_hypothesis_ast` object.
#' @param occurrences whether to return one row per symbol occurrence.
#' @param mapping named character vector from exact old roots to new roots.
#' @param catalog optional `BayesTools_parameter_catalog`. When supplied to
#'   `hypothesis_parse()`, exact non-syntactic public aliases are recognized
#'   without requiring the caller to add backticks.
#' @param namespace optional exact catalog namespace filter.
#' @param component optional exact catalog component filter for unqualified
#'   symbols. A level-qualified symbol such as `term[level]` supplies its own
#'   per-occurrence component and must agree with this value when both are used.
#' @param simplify_names whether to recognize centrally generated simplified
#'   random-effect aliases. Defaults to `FALSE`.
#'
#' @return `hypothesis_parse()` and `hypothesis_rewrite()` return a
#' `BayesTools_hypothesis_ast`. `hypothesis_render()` returns character text.
#' `hypothesis_symbols()` returns unique roots or an occurrence data frame.
#' `hypothesis_resolve()` returns a versioned
#' `BayesTools_hypothesis_resolution` containing one row per occurrence.
#'
#' @export hypothesis_parse
#' @export hypothesis_ast_schema
#' @export hypothesis_render
#' @export hypothesis_symbols
#' @export hypothesis_rewrite
#' @export hypothesis_resolve
#' @name hypothesis_ast
NULL

#' @rdname hypothesis_ast
hypothesis_parse <- function(hypothesis, catalog = NULL, namespace = NULL,
                             component = NULL, simplify_names = FALSE){

  check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  if(length(hypothesis) == 0L){
    stop("'hypothesis' must contain at least one statement.", call. = FALSE)
  }
  if(!is.null(catalog)){
    hypothesis <- .bt_hypothesis_quote_catalog_aliases(
      hypothesis = hypothesis,
      catalog    = catalog,
      namespace  = namespace,
      component  = component,
      simplify_names = simplify_names
    )
  }else{
    check_char(namespace, "namespace", check_length = 1L, allow_NULL = TRUE,
               allow_NA = FALSE)
    check_char(component, "component", check_length = 1L, allow_NULL = TRUE,
               allow_NA = FALSE)
    check_bool(simplify_names, "simplify_names", allow_NA = FALSE)
    if(!is.null(namespace) || !is.null(component)){
      stop("'namespace' and 'component' require 'catalog'.", call. = FALSE)
    }
  }
  statements <- lapply(hypothesis, function(statement){
    .bt_hypothesis_ast_statement_from_spec(
      .hypothesis_parse_statement_spec(statement)
    )
  })
  out <- list(
    schema_version = .bt_hypothesis_ast_version,
    statements = statements
  )
  class(out) <- c("BayesTools_hypothesis_ast", "list")
  .bt_validate_hypothesis_ast(out)
  out
}

.bt_hypothesis_quote_catalog_aliases <- function(
    hypothesis, catalog, namespace, component, simplify_names){

  .bt_validate_parameter_catalog(catalog)
  check_char(namespace, "namespace", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(component, "component", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)

  quantities <- catalog$quantities
  public <- !quantities$internal
  canonical <- public
  if(!is.null(namespace)){
    canonical <- canonical & quantities$namespace == namespace
  }
  if(!is.null(component)){
    canonical <- canonical & quantities$component == component
  }

  aliases <- catalog$aliases
  alias_rows <- aliases$quantity_id %in% quantities$quantity_id[public] &
    (!aliases$simplified | simplify_names)
  if(!is.null(namespace)){
    alias_rows <- alias_rows & aliases$namespace == namespace
  }
  if(!is.null(component)){
    alias_rows <- alias_rows & aliases$component == component
  }
  names <- unique(c(
    quantities$canonical_name[canonical],
    aliases$alias[alias_rows]
  ))
  names <- names[
    nzchar(names) &
      make.names(names) != names
  ]
  if(length(names) == 0L){
    return(hypothesis)
  }
  names <- names[order(nchar(names), decreasing = TRUE)]

  vapply(
    hypothesis,
    .bt_hypothesis_quote_statement_aliases,
    character(1),
    aliases = names,
    USE.NAMES = FALSE
  )
}

.bt_hypothesis_quote_statement_aliases <- function(statement, aliases){

  n <- nchar(statement, type = "chars")
  if(n == 0L){
    return(statement)
  }
  out <- character()
  i <- 1L
  quote <- ""
  escaped <- FALSE
  while(i <= n){
    current <- substr(statement, i, i)
    if(nzchar(quote)){
      out <- c(out, current)
      if(escaped){
        escaped <- FALSE
      }else if(identical(current, "\\")){
        escaped <- TRUE
      }else if(identical(current, quote)){
        quote <- ""
      }
      i <- i + 1L
      next
    }
    if(current %in% c("`", "'", "\"")){
      quote <- current
      out <- c(out, current)
      i <- i + 1L
      next
    }

    matched <- NULL
    for(alias in aliases){
      end <- i + nchar(alias, type = "chars") - 1L
      if(end <= n && identical(substr(statement, i, end), alias) &&
         .bt_hypothesis_alias_boundary(statement, i, end)){
        matched <- alias
        break
      }
    }
    if(!is.null(matched)){
      out <- c(out, encodeString(matched, quote = "`"))
      i <- i + nchar(matched, type = "chars")
      next
    }
    out <- c(out, current)
    i <- i + 1L
  }

  paste0(out, collapse = "")
}

.bt_hypothesis_alias_boundary <- function(statement, start, end){

  identifier <- "[[:alnum:]_.]"
  first <- substr(statement, start, start)
  last  <- substr(statement, end, end)
  before <- if(start > 1L) substr(statement, start - 1L, start - 1L) else ""
  after  <- if(end < nchar(statement, type = "chars")) {
    substr(statement, end + 1L, end + 1L)
  }else{
    ""
  }
  left_ok <- !grepl(identifier, first) || !grepl(identifier, before)
  right_ok <- !grepl(identifier, last) || !grepl(identifier, after)

  left_ok && right_ok
}

#' @rdname hypothesis_ast
hypothesis_ast_schema <- function(){

  data.frame(
    object = c(
      rep("ast", 2L), rep("statement", 4L), rep("side", 5L),
      rep("node", 4L), rep("resolution", 4L)
    ),
    field = c(
      "schema_version", "statements",
      "source", "explicit", "left", "right",
      "type", "source", "label", "expression", "value",
      "type", "source", "operator/name", "children/value",
      "schema_version", "ast_schema_version", "parameter_map_version",
      "occurrences"
    ),
    type = c(
      "integer", "list", "character", "logical", "side", "side",
      "character", "character", "character", "node", "numeric/NULL",
      "character", "character", "character/NULL", "node/list/scalar",
      "integer", "integer", "integer", "data.frame"
    ),
    description = c(
      "Hypothesis AST schema version.",
      "One entry per input hypothesis statement.",
      "Normalized source label.",
      "Whether the statement contained an explicit top-level 'vs'.",
      "Alternative-side node.",
      "Null-side node, explicit or generated as a complement.",
      "Point, not_point, or region.",
      "Stable side source label.",
      "Stable display label.",
      "Arithmetic expression or region condition node.",
      "Finite point value; NULL for regions.",
      "Literal, symbol, level_reference, parentheses, arithmetic, function, comparison, boolean, or negation.",
      "Stable node source label.",
      "Operator or whitelisted function name when applicable.",
      "Typed node payload.",
      "Resolution schema version.",
      "AST schema version used for resolution.",
      "Parameter-map schema version used for resolution.",
      "Occurrence-to-quantity mapping table."
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname hypothesis_ast
hypothesis_render <- function(ast){

  .bt_validate_hypothesis_ast(ast)
  vapply(ast$statements, .bt_hypothesis_render_statement, character(1))
}

#' @rdname hypothesis_ast
hypothesis_symbols <- function(ast, occurrences = FALSE){

  .bt_validate_hypothesis_ast(ast)
  check_bool(occurrences, "occurrences", allow_NA = FALSE)
  rows <- list()
  occurrence_i <- 0L
  for(statement_i in seq_along(ast$statements)){
    statement <- ast$statements[[statement_i]]
    for(side_name in c("left", "right")){
      nodes <- .bt_hypothesis_symbol_nodes(statement[[side_name]]$expression)
      for(node in nodes){
        occurrence_i <- occurrence_i + 1L
        level <- if(identical(node$type, "level_reference")){
          node$level
        }else{
          NA_character_
        }
        parameter <- if(identical(node$type, "level_reference")){
          node$parameter
        }else{
          node$name
        }
        symbol <- if(identical(node$type, "level_reference")){
          paste0(node$parameter, "[", node$level, "]")
        }else{
          node$name
        }
        rows[[occurrence_i]] <- data.frame(
          statement = as.integer(statement_i),
          side = side_name,
          occurrence = as.integer(occurrence_i),
          symbol = symbol,
          parameter = parameter,
          level = level,
          source = node$source,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if(length(rows) == 0L){
    out <- data.frame(
      statement = integer(), side = character(), occurrence = integer(),
      symbol = character(), parameter = character(), level = character(),
      source = character(), stringsAsFactors = FALSE
    )
  }else{
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
  }
  if(isTRUE(occurrences)){
    return(out)
  }
  unique(out$parameter)
}

#' @rdname hypothesis_ast
hypothesis_rewrite <- function(ast, mapping){

  .bt_validate_hypothesis_ast(ast)
  if(!is.character(mapping) || anyNA(mapping) || any(!nzchar(mapping)) ||
     is.null(names(mapping)) || anyNA(names(mapping)) ||
     any(!nzchar(names(mapping))) || anyDuplicated(names(mapping))){
    stop("'mapping' must be a named non-missing character vector with unique non-empty source names.",
         call. = FALSE)
  }
  if(length(mapping) == 0L){
    return(ast)
  }
  round_trips <- vapply(
    unname(mapping),
    .bt_hypothesis_root_round_trips,
    logical(1)
  )
  if(any(!round_trips)){
    stop(
      "Rewrite targets must be parameter roots that preserve their exact identity through rendering and parsing: ",
      paste0("'", unname(mapping[!round_trips]), "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  roots <- hypothesis_symbols(ast)
  unknown <- setdiff(names(mapping), roots)
  if(length(unknown) > 0L){
    stop(
      "Rewrite mapping references unknown hypothesis symbol",
      if(length(unknown) > 1L) "s: " else ": ",
      paste0("'", unknown, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  rewritten_roots <- roots
  matches <- match(roots, names(mapping))
  replace <- !is.na(matches)
  rewritten_roots[replace] <- unname(mapping[matches[replace]])
  if(anyDuplicated(rewritten_roots)){
    stop("Rewrite mapping creates duplicate or colliding hypothesis symbols.",
         call. = FALSE)
  }

  out <- ast
  out$statements <- lapply(out$statements, function(statement){
    statement$left$expression <- .bt_hypothesis_rewrite_node(
      statement$left$expression,
      mapping
    )
    statement$right$expression <- .bt_hypothesis_rewrite_node(
      statement$right$expression,
      mapping
    )
    statement$left <- .bt_hypothesis_refresh_side(statement$left)
    statement$right <- .bt_hypothesis_refresh_side(statement$right)
    statement$source <- .bt_hypothesis_render_statement(statement)
    statement
  })
  .bt_validate_hypothesis_ast(out)
  out
}

#' @rdname hypothesis_ast
hypothesis_resolve <- function(ast, catalog, namespace = NULL,
                               component = NULL, simplify_names = FALSE){

  .bt_validate_hypothesis_ast(ast)
  .bt_validate_parameter_catalog(catalog)
  check_char(namespace, "namespace", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(component, "component", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)
  occurrences <- hypothesis_symbols(ast, occurrences = TRUE)
  if(nrow(occurrences) == 0L){
    stop("The hypothesis contains no parameter symbols to resolve.",
         call. = FALSE)
  }
  resolved <- vector("list", nrow(occurrences))
  for(i in seq_len(nrow(occurrences))){
    occurrence_component <- component
    if(!is.na(occurrences$level[i])){
      if(!is.null(component) && !identical(component, occurrences$level[i])){
        stop(
          "The level in hypothesis symbol '", occurrences$symbol[i],
          "' does not match the requested catalog component '", component,
          "'.",
          call. = FALSE
        )
      }
      occurrence_component <- occurrences$level[i]
    }
    selection <- parameter_catalog_resolve(
      catalog,
      alias = occurrences$parameter[i],
      namespace = namespace,
      component = occurrence_component,
      simplify_names = simplify_names
    )
    quantity <- selection$quantities
    resolved[[i]] <- data.frame(
      statement = occurrences$statement[i],
      side = occurrences$side[i],
      occurrence = occurrences$occurrence[i],
      symbol = occurrences$symbol[i],
      parameter = occurrences$parameter[i],
      level = occurrences$level[i],
      quantity_id = quantity$quantity_id,
      canonical_name = quantity$canonical_name,
      namespace = quantity$namespace,
      component = quantity$component,
      stringsAsFactors = FALSE
    )
  }
  out <- list(
    schema_version = .bt_hypothesis_resolution_version,
    ast_schema_version = ast$schema_version,
    parameter_map_version = catalog$schema_version,
    occurrences = do.call(rbind, resolved)
  )
  rownames(out$occurrences) <- NULL
  class(out) <- c("BayesTools_hypothesis_resolution", "list")
  out
}

.bt_hypothesis_root_round_trips <- function(root){

  node <- list(type = "symbol", source = root, name = root)
  reparsed <- tryCatch({
    rendered <- .bt_hypothesis_render_node(node)
    normalized <- .hypothesis_normalize_level_references(rendered)
    .bt_hypothesis_ast_node(.hypothesis_parse_expression(normalized))
  }, error = function(e){
    NULL
  })
  !is.null(reparsed) && identical(reparsed$type, "symbol") &&
    identical(reparsed$name, root)
}

.bt_hypothesis_ast_statement_from_spec <- function(spec){

  left <- .bt_hypothesis_ast_side(spec$left)
  right <- .bt_hypothesis_ast_side(spec$right)
  source <- if(isTRUE(spec$explicit)){
    paste(.bt_hypothesis_render_side(left), "vs",
          .bt_hypothesis_render_side(right))
  }else{
    .bt_hypothesis_render_side(left)
  }
  list(
    source = source,
    explicit = isTRUE(spec$explicit),
    left = left,
    right = right
  )
}

.bt_hypothesis_ast_side <- function(side){

  if(side$type %in% c("point", "not_point")){
    expression <- side$expression
    if(is.null(expression)){
      expression <- .hypothesis_parse_expression(side$expr)
    }
    out <- list(
      type = side$type,
      source = side$label,
      label = side$label,
      expression = .bt_hypothesis_ast_node(expression),
      value = as.numeric(side$value)
    )
  }else{
    expression <- side$condition_expression
    if(is.null(expression)){
      expression <- .hypothesis_parse_expression(side$condition)
      if(isTRUE(side$complement)){
        expression <- as.call(list(
          as.name("!"),
          as.call(list(as.name("("), expression))
        ))
      }
    }
    out <- list(
      type = "region",
      source = side$label,
      label = side$label,
      expression = .bt_hypothesis_ast_node(expression),
      value = NULL
    )
  }
  out$source <- .bt_hypothesis_render_side(out)
  out
}

.bt_hypothesis_ast_node <- function(expr){

  source <- .hypothesis_expression_text(expr)
  if(is.numeric(expr) || is.integer(expr)){
    return(list(type = "literal", source = source, value = as.numeric(expr)))
  }
  if(is.name(expr)){
    name <- .hypothesis_decode_escaped_constant(as.character(expr))
    reference <- .hypothesis_parse_level_symbol(name)
    if(isTRUE(reference$direct)){
      return(list(
        type = "level_reference",
        source = source,
        parameter = reference$parameter,
        level = reference$level
      ))
    }
    return(list(type = "symbol", source = source, name = name))
  }
  fun <- .hypothesis_call_name(expr)
  if(identical(fun, "(")){
    return(list(
      type = "parentheses",
      source = source,
      expression = .bt_hypothesis_ast_node(expr[[2L]])
    ))
  }
  if(identical(fun, "!")){
    return(list(
      type = "negation",
      source = source,
      expression = .bt_hypothesis_ast_node(expr[[2L]])
    ))
  }
  if(fun %in% c("&", "|")){
    return(list(
      type = "boolean",
      source = source,
      operator = fun,
      left = .bt_hypothesis_ast_node(expr[[2L]]),
      right = .bt_hypothesis_ast_node(expr[[3L]])
    ))
  }
  if(fun %in% c("<", "<=", ">", ">=")){
    return(list(
      type = "comparison",
      source = source,
      operator = fun,
      left = .bt_hypothesis_ast_node(expr[[2L]]),
      right = .bt_hypothesis_ast_node(expr[[3L]])
    ))
  }
  if(fun %in% c("+", "-", "*", "/", "^")){
    return(list(
      type = "arithmetic",
      source = source,
      operator = fun,
      arguments = lapply(as.list(expr[-1L]), .bt_hypothesis_ast_node)
    ))
  }
  list(
    type = "function",
    source = source,
    name = fun,
    arguments = lapply(as.list(expr[-1L]), .bt_hypothesis_ast_node)
  )
}

.bt_hypothesis_node_language <- function(node){

  switch(
    node$type,
    literal = node$value,
    symbol = as.name(if(node$name %in% .hypothesis_escaped_constant_names()){
      .hypothesis_escaped_constant_symbol(node$name)
    }else{
      node$name
    }),
    level_reference = as.name(paste0(
      node$parameter, "[", node$level, "]"
    )),
    parentheses = as.call(list(
      as.name("("),
      .bt_hypothesis_node_language(node$expression)
    )),
    negation = as.call(list(
      as.name("!"),
      .bt_hypothesis_node_language(node$expression)
    )),
    boolean = as.call(list(
      as.name(node$operator),
      .bt_hypothesis_node_language(node$left),
      .bt_hypothesis_node_language(node$right)
    )),
    comparison = as.call(list(
      as.name(node$operator),
      .bt_hypothesis_node_language(node$left),
      .bt_hypothesis_node_language(node$right)
    )),
    arithmetic = as.call(c(
      list(as.name(node$operator)),
      lapply(node$arguments, .bt_hypothesis_node_language)
    )),
    "function" = as.call(c(
      list(as.name(node$name)),
      lapply(node$arguments, .bt_hypothesis_node_language)
    )),
    stop("Unsupported hypothesis AST node type.", call. = FALSE)
  )
}

.bt_hypothesis_render_node <- function(node){

  .hypothesis_expression_text(.bt_hypothesis_node_language(node))
}

.bt_hypothesis_render_side <- function(side){

  expression <- .bt_hypothesis_render_node(side$expression)
  if(identical(side$type, "point")){
    return(paste(expression, "=", .hypothesis_number_label(side$value)))
  }
  if(identical(side$type, "not_point")){
    return(paste(expression, "!=", .hypothesis_number_label(side$value)))
  }
  expression
}

.bt_hypothesis_render_statement <- function(statement){

  left <- .bt_hypothesis_render_side(statement$left)
  if(isTRUE(statement$explicit)){
    return(paste(left, "vs", .bt_hypothesis_render_side(statement$right)))
  }
  left
}

.bt_hypothesis_refresh_side <- function(side){

  rendered <- .bt_hypothesis_render_side(side)
  side$source <- rendered
  side$label <- .hypothesis_display_text(rendered)
  side
}

.bt_hypothesis_symbol_nodes <- function(node){

  if(node$type %in% c("symbol", "level_reference")){
    return(list(node))
  }
  if(node$type %in% c("parentheses", "negation")){
    return(.bt_hypothesis_symbol_nodes(node$expression))
  }
  if(node$type %in% c("boolean", "comparison")){
    return(c(
      .bt_hypothesis_symbol_nodes(node$left),
      .bt_hypothesis_symbol_nodes(node$right)
    ))
  }
  if(node$type %in% c("arithmetic", "function")){
    return(unlist(lapply(node$arguments, .bt_hypothesis_symbol_nodes),
                  recursive = FALSE))
  }
  list()
}

.bt_hypothesis_rewrite_node <- function(node, mapping){

  if(identical(node$type, "symbol")){
    if(node$name %in% names(mapping)){
      node$name <- unname(mapping[[node$name]])
    }
  }else if(identical(node$type, "level_reference")){
    if(node$parameter %in% names(mapping)){
      node$parameter <- unname(mapping[[node$parameter]])
    }
  }else if(node$type %in% c("parentheses", "negation")){
    node$expression <- .bt_hypothesis_rewrite_node(node$expression, mapping)
  }else if(node$type %in% c("boolean", "comparison")){
    node$left <- .bt_hypothesis_rewrite_node(node$left, mapping)
    node$right <- .bt_hypothesis_rewrite_node(node$right, mapping)
  }else if(node$type %in% c("arithmetic", "function")){
    node$arguments <- lapply(
      node$arguments,
      .bt_hypothesis_rewrite_node,
      mapping = mapping
    )
  }
  node$source <- .bt_hypothesis_render_node(node)
  node
}

.bt_validate_hypothesis_ast <- function(ast){

  valid <- inherits(ast, "BayesTools_hypothesis_ast") && is.list(ast) &&
    identical(names(ast), c("schema_version", "statements")) &&
    identical(ast$schema_version, .bt_hypothesis_ast_version) &&
    is.list(ast$statements) && length(ast$statements) > 0L
  if(!valid){
    stop("Hypothesis AST metadata are missing or unsupported. Parse the hypothesis again with this version of BayesTools.",
         call. = FALSE)
  }
  for(statement in ast$statements){
    valid_statement <- is.list(statement) && identical(
      names(statement),
      c("source", "explicit", "left", "right")
    ) && is.character(statement$source) && length(statement$source) == 1L &&
      !is.na(statement$source) && nzchar(statement$source) &&
      is.logical(statement$explicit) && length(statement$explicit) == 1L &&
      !is.na(statement$explicit)
    if(!valid_statement){
      stop("Hypothesis AST contains a malformed statement. Parse the hypothesis again.",
           call. = FALSE)
    }
    .bt_validate_hypothesis_side(statement$left)
    .bt_validate_hypothesis_side(statement$right)
    if(!identical(statement$source,
                  .bt_hypothesis_render_statement(statement))){
      stop("Hypothesis AST statement source is inconsistent with its nodes.",
           call. = FALSE)
    }
    if(!isTRUE(statement$explicit)){
      expected <- .bt_hypothesis_complement_side(statement$left)
      if(!identical(statement$right$type, expected$type) ||
         !identical(.bt_hypothesis_render_side(statement$right),
                    .bt_hypothesis_render_side(expected))){
        stop("Hypothesis AST implicit comparison has an invalid complement.",
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.bt_hypothesis_complement_side <- function(side) {

  if(side$type %in% c("point", "not_point")) {
    out <- side
    out$type <- if(identical(side$type, "point")) "not_point" else "point"
    return(.bt_hypothesis_refresh_side(out))
  }

  expression <- .bt_hypothesis_node_language(side$expression)
  simple <- .hypothesis_unwrap_parentheses(expression)
  operator <- .hypothesis_call_name(simple)
  if(!is.null(operator) && operator %in% c("<", "<=", ">", ">=")) {
    complement <- switch(
      operator,
      ">"  = "<=",
      ">=" = "<",
      "<"  = ">=",
      "<=" = ">"
    )
    expression <- as.call(list(
      as.name(complement),
      simple[[2L]],
      simple[[3L]]
    ))
  } else {
    expression <- as.call(list(
      as.name("!"),
      as.call(list(as.name("("), expression))
    ))
  }
  out <- side
  out$expression <- .bt_hypothesis_ast_node(expression)
  .bt_hypothesis_refresh_side(out)
}

.bt_validate_hypothesis_side <- function(side){

  valid <- is.list(side) && identical(
    names(side),
    c("type", "source", "label", "expression", "value")
  ) && is.character(side$type) && length(side$type) == 1L &&
    side$type %in% c("point", "not_point", "region") &&
    is.character(side$source) && length(side$source) == 1L &&
    !is.na(side$source) && nzchar(side$source) &&
    is.character(side$label) && length(side$label) == 1L &&
    !is.na(side$label) && nzchar(side$label)
  if(!valid){
    stop("Hypothesis AST contains a malformed side. Parse the hypothesis again.",
         call. = FALSE)
  }
  .bt_validate_hypothesis_node(side$expression)
  expression <- .bt_hypothesis_node_language(side$expression)
  if(side$type %in% c("point", "not_point")){
    if(!is.numeric(side$value) || length(side$value) != 1L ||
       !is.finite(side$value)){
      stop("Hypothesis AST point node has no finite scalar value.",
           call. = FALSE)
    }
    .hypothesis_validate_expression(expression, condition = FALSE)
  }else{
    if(!is.null(side$value)){
      stop("Hypothesis AST region node must not contain a point value.",
           call. = FALSE)
    }
    .hypothesis_validate_expression(expression, condition = TRUE)
  }
  if(!identical(side$source, .bt_hypothesis_render_side(side))){
    stop("Hypothesis AST side source is inconsistent with its expression.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.bt_validate_hypothesis_node <- function(node){

  if(!is.list(node) || !is.character(node$type) || length(node$type) != 1L ||
     !is.character(node$source) || length(node$source) != 1L ||
     is.na(node$source) || !nzchar(node$source)){
    stop("Hypothesis AST contains a malformed expression node. Parse the hypothesis again.",
         call. = FALSE)
  }
  expected <- switch(
    node$type,
    literal = c("type", "source", "value"),
    symbol = c("type", "source", "name"),
    level_reference = c("type", "source", "parameter", "level"),
    parentheses = c("type", "source", "expression"),
    negation = c("type", "source", "expression"),
    boolean = c("type", "source", "operator", "left", "right"),
    comparison = c("type", "source", "operator", "left", "right"),
    arithmetic = c("type", "source", "operator", "arguments"),
    "function" = c("type", "source", "name", "arguments"),
    NULL
  )
  if(is.null(expected) || !identical(names(node), expected)){
    stop("Hypothesis AST contains an unsupported expression node. Parse the hypothesis again.",
         call. = FALSE)
  }
  if(identical(node$type, "literal")){
    if(!is.numeric(node$value) || length(node$value) != 1L ||
       !is.finite(node$value)){
      stop("Hypothesis AST literals must be finite numeric scalars.",
           call. = FALSE)
    }
  }else if(identical(node$type, "symbol")){
    if(!is.character(node$name) || length(node$name) != 1L ||
       is.na(node$name) || !nzchar(node$name)){
      stop("Hypothesis AST symbol name is malformed.", call. = FALSE)
    }
  }else if(identical(node$type, "level_reference")){
    values <- c(node$parameter, node$level)
    if(!is.character(values) || length(values) != 2L || anyNA(values) ||
       any(!nzchar(values))){
      stop("Hypothesis AST level reference is malformed.", call. = FALSE)
    }
  }else if(node$type %in% c("parentheses", "negation")){
    .bt_validate_hypothesis_node(node$expression)
  }else if(node$type %in% c("boolean", "comparison")){
    allowed <- if(identical(node$type, "boolean")){
      c("&", "|")
    }else{
      c("<", "<=", ">", ">=")
    }
    if(!is.character(node$operator) || length(node$operator) != 1L ||
       !node$operator %in% allowed){
      stop("Hypothesis AST operator is malformed.", call. = FALSE)
    }
    .bt_validate_hypothesis_node(node$left)
    .bt_validate_hypothesis_node(node$right)
  }else{
    operator <- if(identical(node$type, "arithmetic")){
      node$operator
    }else{
      node$name
    }
    allowed <- if(identical(node$type, "arithmetic")){
      c("+", "-", "*", "/", "^")
    }else{
      c("abs", "exp", "log", "sqrt", "plogis", "qlogis")
    }
    if(!is.character(operator) || length(operator) != 1L ||
       !operator %in% allowed || !is.list(node$arguments) ||
       length(node$arguments) == 0L){
      stop("Hypothesis AST call node is malformed.", call. = FALSE)
    }
    for(argument in node$arguments){
      .bt_validate_hypothesis_node(argument)
    }
  }
  if(!identical(node$source, .bt_hypothesis_render_node(node))){
    stop("Hypothesis AST node source is inconsistent with its expression.",
         call. = FALSE)
  }
  invisible(TRUE)
}
