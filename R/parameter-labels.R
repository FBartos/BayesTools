# One renderer for parameter labels.
#
# Every parameter label that is shown or typed - catalog selectors and
# aliases, table rows, mixed-posterior column names, plot legends, and warning
# text - is rendered by .bt_label() from structured label parts. The parts of
# fitted quantities are stored on the parameter catalog ('label_parts'), the
# parts of mixed and marginal posterior columns in their draw metadata
# ('quantities'). Consumers read the parts; they never parse label text.

.bt_label_styles <- c("selector", "table", "plot", "warning")
.bt_label_transformations <- c("none", "dif", "exp")

# Structured label parts of one quantity.
#
# - formula_parameter: owning formula output parameter, or "".
# - components: the term components (formula term `g:x` has c("g", "x")); a
#   quantity without formula term structure has one component, its name.
# - levels: named character, the level label of each factor component of a
#   level cell (or an estimated marginal mean); empty otherwise.
# - coefficient: contrast coefficient index `j` of a coordinate that is not a
#   level cell (`term{j}`), or NA.
# - transformation: "none", "dif" (a transformed contrast level, a difference
#   from the mean), or "exp" (the exponentiated log intercept).
# - component: a mixture component shown after the term (`term[component]`),
#   or "".
# - inclusion: NA, or the inclusion row of the quantity: "" for a
#   spike-and-slab inclusion, otherwise the mixture component.
# - random: NULL, or the random-effect semantic name fields `owner`,
#   `quantity`, `arguments`, and `display_arguments`.
# - marginal: whether the parts describe an estimated marginal mean at
#   `levels` (a different estimand than the coefficient level cell).
# - selector: the exact selector when it is not derived from the other parts
#   (fitted coordinate names), or "".
.bt_label_parts <- function(components, formula_parameter = "",
                            levels = character(), coefficient = NA_integer_,
                            transformation = "none", component = "",
                            inclusion = NA_character_, random = NULL,
                            marginal = FALSE, selector = ""){

  parts <- list(
    formula_parameter = as.character(formula_parameter),
    components        = as.character(components),
    levels            = if(length(levels) == 0L){
      stats::setNames(character(), character())
    }else{
      stats::setNames(as.character(levels), names(levels))
    },
    coefficient       = as.integer(coefficient),
    transformation    = transformation,
    component         = as.character(component),
    inclusion         = as.character(inclusion),
    random            = random,
    marginal          = marginal,
    selector          = as.character(selector)
  )
  class(parts) <- c("BayesTools_label_parts", "list")
  .bt_validate_label_parts(parts)
  parts
}

.bt_label_parts_fields <- c(
  "formula_parameter", "components", "levels", "coefficient",
  "transformation", "component", "inclusion", "random", "marginal", "selector"
)

.bt_validate_label_parts <- function(parts){

  scalar_character <- function(x){
    is.character(x) && length(x) == 1L && !is.na(x)
  }
  valid <- inherits(parts, "BayesTools_label_parts") &&
    identical(names(parts), .bt_label_parts_fields) &&
    scalar_character(parts$formula_parameter) &&
    is.character(parts$components) && length(parts$components) > 0L &&
    !anyNA(parts$components) && all(nzchar(parts$components)) &&
    is.character(parts$levels) && !anyNA(parts$levels) &&
    (length(parts$levels) == 0L || (
      !is.null(names(parts$levels)) && !anyNA(names(parts$levels)) &&
        !anyDuplicated(names(parts$levels)) &&
        all(names(parts$levels) %in% parts$components))) &&
    is.integer(parts$coefficient) && length(parts$coefficient) == 1L &&
    (is.na(parts$coefficient) || parts$coefficient >= 1L) &&
    (is.na(parts$coefficient) || length(parts$levels) == 0L) &&
    scalar_character(parts$transformation) &&
    parts$transformation %in% .bt_label_transformations &&
    scalar_character(parts$component) &&
    is.character(parts$inclusion) && length(parts$inclusion) == 1L &&
    is.logical(parts$marginal) && length(parts$marginal) == 1L &&
    !is.na(parts$marginal) &&
    scalar_character(parts$selector) &&
    (is.null(parts$random) || (
      is.list(parts$random) &&
        identical(names(parts$random),
                  c("owner", "quantity", "arguments", "display_arguments")) &&
        scalar_character(parts$random$owner) &&
        scalar_character(parts$random$quantity) &&
        is.character(parts$random$arguments) &&
        !anyNA(parts$random$arguments) &&
        is.character(parts$random$display_arguments) &&
        !anyNA(parts$random$display_arguments)))
  if(!isTRUE(valid)){
    stop(
      "Parameter label parts are malformed. Refit or rebuild the catalog with this version of BayesTools.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_label_parts_list <- function(parts){

  if(inherits(parts, "BayesTools_label_parts")){
    return(list(parts))
  }
  if(!is.list(parts)){
    stop("Parameter label parts must be a list of label parts.", call. = FALSE)
  }

  parts
}

# Modify fields of label parts (vectorized over a list of parts).
.bt_label_parts_update <- function(parts, ...){

  fields <- list(...)
  lapply(.bt_label_parts_list(parts), function(part){
    for(field in names(fields)){
      part[field] <- list(fields[[field]])
    }
    if("coefficient" %in% names(fields)){
      part$coefficient <- as.integer(part$coefficient)
    }
    .bt_validate_label_parts(part)
    part
  })
}

#' @title Parameter labels
#'
#' @description Renders the labels of parameter quantities from their
#' structured label parts: catalog quantities (the `label_parts` column of
#' [parameter_catalog()]) and the columns of mixed and marginal posterior
#' draws (the `quantities` draw metadata, see [posterior_metadata()]). The same
#' renderer produces the catalog aliases, the rows of the summary tables, the
#' column names of mixed posteriors, plot legends, and warning text, so a label
#' shown by BayesTools is always a label of the same quantity everywhere.
#'
#' @param x a parameter catalog (or its `quantities` table), a parameter
#' selection returned by [parameter_catalog_resolve()], a data frame with a
#' `label_parts` column, or a list of label parts.
#' @param style one of `"selector"` (the exact selector, i.e., the canonical
#' name of a catalog quantity), `"table"` (the summary-table row label),
#' `"plot"` (the plot-legend label: the level text of a level cell), or
#' `"warning"` (the row label used in warnings and footnotes).
#' @param formula_prefix whether table and warning labels show the formula
#' parameter, as in `(mu) x`. Defaults to `TRUE`.
#' @param simplify whether random-effect labels use their simplified display
#' arguments (for example, `sd` for a sole random intercept). Defaults to
#' `FALSE`.
#'
#' @details Square brackets after a factor term hold a level label (`g[b]`);
#' interaction cells name the level of every factor (`g[b]:h[v]`, `g[b]:x`);
#' contrast coefficients that are not level cells are written with curly
#' braces (`g{1}`). Selectors percent-escape syntax-sensitive level characters
#' and quote interaction tokens as the catalog does; table, plot, and warning
#' labels show the level labels as they are.
#'
#' @return a character vector of labels.
#'
#' @seealso [parameter_catalog()]
#'
#' @export
parameter_labels <- function(x, style = c("selector", "table", "plot", "warning"),
                             formula_prefix = TRUE, simplify = FALSE){

  style <- match.arg(style)
  check_bool(formula_prefix, "formula_prefix", allow_NA = FALSE)
  check_bool(simplify, "simplify", allow_NA = FALSE)

  if(inherits(x, "BayesTools_parameter_catalog")){
    .bt_validate_parameter_catalog(x)
    x <- x$quantities
  }else if(inherits(x, "BayesTools_parameter_selection")){
    .bt_validate_parameter_selection(x)
    x <- x$quantities
  }
  if(is.data.frame(x)){
    if(!"label_parts" %in% names(x)){
      stop("The 'x' argument does not contain a 'label_parts' column.",
           call. = FALSE)
    }
    missing <- vapply(x$label_parts, is.null, logical(1))
    out <- character(nrow(x))
    if(any(missing) && "display_label" %in% names(x)){
      # extension providers may describe their quantities by a label only
      out[missing] <- x$display_label[missing]
    }else if(any(missing)){
      stop("The 'x' argument contains quantities without label parts.",
           call. = FALSE)
    }
    if(any(!missing)){
      out[!missing] <- .bt_label(
        unclass(x$label_parts)[!missing],
        style          = style,
        formula_prefix = formula_prefix,
        simplify       = simplify
      )
    }
    return(out)
  }

  .bt_label(
    x,
    style          = style,
    formula_prefix = formula_prefix,
    simplify       = simplify
  )
}

# Render labels of label parts (one label per parts object).
.bt_label <- function(parts, style = "table", formula_prefix = TRUE,
                      simplify = FALSE){

  style <- match.arg(style, .bt_label_styles)
  parts <- .bt_label_parts_list(parts)
  out <- vapply(parts, function(part){
    .bt_validate_label_parts(part)
    switch(
      style,
      selector = .bt_label_selector(part),
      table    = .bt_label_table(part, formula_prefix, simplify),
      plot     = .bt_label_plot(part, simplify),
      warning  = .bt_label_table(part, formula_prefix, simplify)
    )
  }, character(1))

  unname(out)
}

.bt_label_prefix <- function(formula_parameter, formula_prefix){

  if(isTRUE(formula_prefix) && nzchar(formula_parameter)){
    paste0("(", formula_parameter, ") ")
  }else{
    ""
  }
}

.bt_label_jags_base <- function(part){

  base <- paste(part$components, collapse = "__xXx__")
  if(nzchar(part$formula_parameter)){
    base <- paste0(part$formula_parameter, "_", base)
  }
  base
}

.bt_label_suffix <- function(part){

  suffix <- ""
  if(nzchar(part$component)){
    suffix <- paste0(suffix, "[", part$component, "]")
  }
  if(!is.na(part$inclusion)){
    suffix <- paste0(
      suffix,
      if(nzchar(part$inclusion)){
        paste0(" (inclusion: ", part$inclusion, ")")
      }else{
        " (inclusion)"
      }
    )
  }
  suffix
}

# The exact selector: the canonical name of a catalog quantity, and the
# column name of mixed posterior draws.
.bt_label_selector <- function(part){

  if(nzchar(part$selector)){
    return(paste0(part$selector, .bt_label_suffix(part)))
  }
  if(!is.null(part$random)){
    return(.bt_random_effect_semantic_name(
      parameter      = part$formula_parameter,
      owner          = part$random$owner,
      quantity       = part$random$quantity,
      arguments      = part$random$arguments,
      formula_prefix = TRUE
    ))
  }

  base <- .bt_label_jags_base(part)
  label <- if(identical(part$transformation, "exp")){
    paste0(
      if(nzchar(part$formula_parameter)) paste0(part$formula_parameter, "_"),
      "exp(", paste(part$components, collapse = "__xXx__"), ")"
    )
  }else if(!is.na(part$coefficient)){
    paste0(base, "{", part$coefficient, "}")
  }else if(length(part$levels) == 0L){
    base
  }else if(part$marginal){
    paste0(base, "[", paste(part$levels, collapse = ", "), "]")
  }else if(identical(part$transformation, "dif")){
    # the level name of a transformed contrast marks every factor component
    components <- vapply(part$components, function(component){
      if(component %in% names(part$levels)){
        paste0(component, "[dif: ", part$levels[[component]], "]")
      }else{
        component
      }
    }, character(1))
    paste0(
      if(nzchar(part$formula_parameter)) paste0(part$formula_parameter, "_"),
      paste(components, collapse = "__xXx__")
    )
  }else{
    paste0(base, "[", .bt_label_cell_token(part$levels), "]")
  }

  paste0(label, .bt_label_suffix(part))
}

# The catalog token of a level cell: the escaped level of a single factor, or
# the escaped `factor=level` pairs of an interaction.
.bt_label_cell_token <- function(levels){

  if(length(levels) == 1L){
    return(.bt_label_token(levels[[1L]]))
  }
  paste0(
    .bt_label_token(names(levels)), "=", .bt_label_token(unname(levels)),
    collapse = ", "
  )
}

.bt_label_term <- function(part, dif = FALSE){

  components <- vapply(part$components, function(component){
    if(component %in% names(part$levels) && !part$marginal){
      paste0(
        component, "[",
        if(dif) "dif: " else "",
        part$levels[[component]], "]"
      )
    }else{
      component
    }
  }, character(1))
  paste(components, collapse = ":")
}

.bt_label_table <- function(part, formula_prefix, simplify){

  if(!is.null(part$random)){
    return(.bt_random_effect_semantic_name(
      parameter      = part$formula_parameter,
      owner          = part$random$owner,
      quantity       = part$random$quantity,
      arguments      = if(simplify){
        part$random$display_arguments
      }else{
        part$random$arguments
      },
      formula_prefix = formula_prefix
    ))
  }

  prefix <- .bt_label_prefix(part$formula_parameter, formula_prefix)
  label <- if(identical(part$transformation, "exp")){
    paste0("exp(", paste(part$components, collapse = ":"), ")")
  }else if(!is.na(part$coefficient)){
    paste0(paste(part$components, collapse = ":"), "{", part$coefficient, "}")
  }else if(part$marginal && length(part$levels) > 0L){
    paste0(
      paste(part$components, collapse = ":"),
      "[", paste(part$levels, collapse = ", "), "]"
    )
  }else{
    .bt_label_term(part, dif = identical(part$transformation, "dif"))
  }

  paste0(prefix, label, .bt_label_suffix(part))
}

# Plot legends distinguish the levels (cells) of one plotted term: the level
# text of every factor component, `b` or `b, v`.
.bt_label_plot <- function(part, simplify){

  if(length(part$levels) > 0L){
    return(paste(part$levels, collapse = ", "))
  }
  if(!is.na(part$coefficient)){
    return(paste0("{", part$coefficient, "}"))
  }

  .bt_label_table(part, formula_prefix = FALSE, simplify = simplify)
}

# Codec of catalog level tokens --------------------------------------------
#
# Level labels keep their text; the syntax-sensitive characters of selectors
# (square and curly brackets, backticks, quotes, backslashes, the escape
# character itself, and control characters) are percent-escaped. Tokens with
# ',' or '=' or edge whitespace are additionally quoted so that interaction
# tokens split unambiguously. .bt_label_token_decode() inverts the encoding.

.bt_label_token_escapes <- c(
  "%"  = "%25",
  "["  = "%5B",
  "]"  = "%5D",
  "`"  = "%60",
  "\\" = "%5C",
  "\"" = "%22",
  "{"  = "%7B",
  "}"  = "%7D",
  "\r" = "%0D",
  "\n" = "%0A",
  "\t" = "%09",
  "\f" = "%0C",
  "\b" = "%08",
  "\a" = "%07",
  "\v" = "%0B"
)

.bt_label_token <- function(x){

  vapply(as.character(x), function(x_i){
    for(character_i in names(.bt_label_token_escapes)){
      x_i <- gsub(character_i, .bt_label_token_escapes[[character_i]], x_i,
                  fixed = TRUE)
    }
    needs_quotes <- !nzchar(x_i) ||
      grepl("^[[:space:]]|[[:space:]]$", x_i) ||
      grepl(",", x_i, fixed = TRUE) ||
      grepl("=", x_i, fixed = TRUE)
    if(needs_quotes){
      encodeString(x_i, quote = "\"")
    }else{
      x_i
    }
  }, character(1), USE.NAMES = FALSE)
}

.bt_label_token_decode <- function(x){

  vapply(as.character(x), function(token){
    if(nchar(token) >= 2L && startsWith(token, "\"") && endsWith(token, "\"")){
      parsed <- tryCatch(
        parse(text = token, keep.source = FALSE),
        error = function(e) NULL
      )
      if(length(parsed) == 1L && is.character(parsed[[1L]]) &&
         length(parsed[[1L]]) == 1L){
        token <- parsed[[1L]]
      }
    }
    escapes <- stats::setNames(
      names(.bt_label_token_escapes),
      unname(.bt_label_token_escapes)
    )
    matches <- gregexpr(paste(names(escapes), collapse = "|"), token)
    regmatches(token, matches) <- lapply(
      regmatches(token, matches),
      function(escaped) unname(escapes[escaped])
    )
    token
  }, character(1), USE.NAMES = FALSE)
}

# Label parts of fixed factor terms ------------------------------------------

# The term components of a factor prior: the formula term components, or the
# prior name for an ordinary factor prior.
.bt_label_factor_components <- function(parameter, prior, term = ""){

  components <- attr(prior, "term_components", exact = TRUE)
  if(is.null(components)){
    components <- attr(prior, "interaction_terms", exact = TRUE)
  }
  if(is.null(components) || length(components) == 0L){
    components <- if(nzchar(term)) term else parameter
  }

  as.character(components)
}

# Label parts of every level cell of a fixed factor term, in design-row order,
# and of its fitted coordinates (a structural level cell, or contrast
# coefficient `j`), from the persisted factor metadata.
.bt_label_parts_factor <- function(parameter, prior, formula_parameter = "",
                                   term = ""){

  prior <- .complete_factor_metadata(prior, parameter)
  design_info <- .factor_term_design_from_metadata(prior)
  level_names <- design_info$level_names
  if(is.null(level_names) || length(level_names) == 0L){
    stop(
      "Factor levels of '", parameter, "' are missing. Refit the model with ",
      "this version of BayesTools.",
      call. = FALSE
    )
  }
  components <- .bt_label_factor_components(parameter, prior, term)
  if(!nzchar(formula_parameter) && length(components) == 1L &&
     !identical(components, parameter)){
    components <- parameter
  }
  if(!all(names(level_names) %in% components)){
    if(length(level_names) == 1L && length(components) == 1L){
      names(level_names) <- components
    }else{
      stop(
        "Factor metadata of '", parameter, "' do not identify the factor ",
        "components of its term. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
  }
  grid <- .factor_cell_grid(level_names)
  design <- as.matrix(design_info$design)
  if(nrow(grid) != nrow(design)){
    stop(
      "Factor metadata of '", parameter, "' do not identify every level cell. ",
      "Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  cells <- lapply(seq_len(nrow(grid)), function(cell){
    .bt_label_parts(
      components        = components,
      formula_parameter = formula_parameter,
      levels            = stats::setNames(
        vapply(grid, function(column) as.character(column[[cell]]), character(1)),
        names(grid)
      )
    )
  })
  direct <- .bt_factor_direct_cells(prior, design)
  coordinates <- lapply(seq_len(ncol(design)), function(coordinate){
    if(!is.na(direct[coordinate])){
      return(cells[[direct[coordinate]]])
    }
    .bt_label_parts(
      components        = components,
      formula_parameter = formula_parameter,
      coefficient       = coordinate
    )
  })

  list(
    cells       = cells,
    coordinates = coordinates,
    direct      = direct
  )
}
