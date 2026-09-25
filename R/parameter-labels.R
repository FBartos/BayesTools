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

# The term components of a factor prior: the formula term components (from
# the prior's term metadata, or from its name when it has none), or the prior
# name for an ordinary factor prior.
.bt_label_factor_components <- function(parameter, prior, term = "",
                                        formula_parameter = ""){

  components <- attr(prior, "term_components", exact = TRUE)
  if(is.null(components)){
    components <- attr(prior, "interaction_terms", exact = TRUE)
  }
  if(is.null(components) || length(components) == 0L){
    components <- if(nzchar(formula_parameter) &&
                     startsWith(parameter, paste0(formula_parameter, "_"))){
      .bt_label_parts_coefficient(parameter, formula_parameter)$components
    }else if(nzchar(term)){
      term
    }else{
      parameter
    }
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
  components <- .bt_label_factor_components(
    parameter,
    prior,
    term,
    formula_parameter = formula_parameter
  )
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

# Label parts of prior-list entries ------------------------------------------

# The fitted coordinate names of one prior-list entry and their label parts,
# in coordinate order: the level cells and contrast coefficients of a fixed
# factor prior, the term of a formula coefficient, and the parameter name of
# any other parameter. The formula parameter is the prior's own 'parameter'
# attribute; names are never matched against other formula parameters.
.bt_label_parts_prior <- function(parameter, prior){

  formula_parameter <- .bt_label_formula_parameter(prior)
  factor_prior <- .bt_parameter_catalog_factor_prior(parameter, prior)
  if(!is.null(factor_prior) &&
     !is.null(attr(factor_prior, "factor_design", exact = TRUE))){
    return(list(
      coordinates = .JAGS_prior_factor_names(parameter, factor_prior),
      parts       = .bt_label_parts_factor(
        parameter         = parameter,
        prior             = factor_prior,
        formula_parameter = formula_parameter
      )$coordinates
    ))
  }
  if(is.prior.vector(prior) || .bt_prior_is_factor_family(prior)){
    coordinates <- .JAGS_prior_factor_names(parameter, prior)
    return(list(
      coordinates = coordinates,
      parts       = lapply(coordinates, function(coordinate){
        .bt_label_parts(coordinate, selector = coordinate)
      })
    ))
  }

  list(
    coordinates = parameter,
    parts       = list(.bt_label_parts_coefficient(
      parameter         = parameter,
      formula_parameter = formula_parameter,
      interaction_terms = if(.is_prior_interaction(prior)){
        attr(prior, "interaction_terms", exact = TRUE)
      }
    ))
  )
}

# The formula parameter owning a prior (its 'parameter' attribute) or owning
# mixed or marginal draws (their 'formula_parameter' draw metadata); "" for
# any other object.
.bt_label_formula_parameter <- function(x){

  formula_parameter <- if(is.prior(x)){
    attr(x, "parameter", exact = TRUE)
  }else{
    .bt_meta_get(x, "formula_parameter")
  }
  if(is.character(formula_parameter) && length(formula_parameter) == 1L &&
     !is.na(formula_parameter)){
    formula_parameter
  }else{
    ""
  }
}

# Label parts of one fitted coefficient named by JAGS_parameter_names(): the
# formula parameter, then the term with interactions encoded. The formula
# parameter is the owner's own; names are never matched against other formula
# parameters.
.bt_label_parts_coefficient <- function(parameter, formula_parameter = "",
                                        interaction_terms = NULL){

  components <- parameter
  if(nzchar(formula_parameter)){
    stem <- paste0(formula_parameter, "_")
    if(!startsWith(parameter, stem) || nchar(parameter) <= nchar(stem)){
      stop(
        "The formula coefficient '", parameter, "' is not named after its ",
        "formula parameter '", formula_parameter, "'.",
        call. = FALSE
      )
    }
    components <- if(length(interaction_terms) > 0L){
      as.character(interaction_terms)
    }else{
      strsplit(
        substring(parameter, nchar(stem) + 1L),
        "__xXx__",
        fixed = TRUE
      )[[1L]]
    }
  }

  .bt_label_parts(
    components        = components,
    formula_parameter = formula_parameter,
    selector          = parameter
  )
}

# Label parts of a prior-list entry as a whole term (no level or contrast
# coefficient), such as its inclusion row or its column in prior summaries.
.bt_label_parts_term <- function(parameter, prior){

  formula_parameter <- .bt_label_formula_parameter(prior)
  factor_prior <- .bt_parameter_catalog_factor_prior(parameter, prior)
  if(!is.null(factor_prior)){
    components <- .bt_label_factor_components(
      parameter,
      factor_prior,
      formula_parameter = formula_parameter
    )
    if(!nzchar(formula_parameter)){
      components <- parameter
    }
    return(.bt_label_parts(
      components        = components,
      formula_parameter = formula_parameter,
      selector          = parameter
    ))
  }
  if(!nzchar(formula_parameter) || is.prior.vector(prior)){
    return(.bt_label_parts(parameter, selector = parameter))
  }

  .bt_label_parts_coefficient(
    parameter         = parameter,
    formula_parameter = formula_parameter,
    interaction_terms = if(.is_prior_interaction(prior)){
      attr(prior, "interaction_terms", exact = TRUE)
    }
  )
}

# Column names of the fitted coordinates of a prior-list entry: their
# selectors, i.e., the canonical catalog names of the level cells and
# contrast coefficients of factor priors.
.bt_label_prior_column_names <- function(parameter, prior){

  .bt_label(.bt_label_parts_prior(parameter, prior)$parts, style = "selector")
}

# Label parts of every level cell of a factor term (a factor prior, or mixed
# draws carrying the factor metadata), in design-row order, as transformed
# contrast levels ("dif") or level effects ("none").
.bt_label_factor_level_parts <- function(parameter, x, transformation = "dif",
                                         formula_parameter = NULL){

  if(is.null(formula_parameter)){
    formula_parameter <- .bt_label_formula_parameter(x)
  }
  .bt_label_parts_update(
    .bt_label_parts_factor(
      parameter         = parameter,
      prior             = x,
      formula_parameter = formula_parameter
    )$cells,
    transformation = transformation
  )
}

# Label parts of the columns of mixed draws 'x' of the element 'name': their
# 'quantities' draw metadata, otherwise the element's own coefficient (vector
# draws) or its column names (other draws).
.bt_draws_label_parts <- function(x, name){

  quantities <- .bt_draws_quantities(x)
  if(!is.null(quantities)){
    return(unclass(quantities$label_parts))
  }
  if(is.null(dim(x))){
    formula_parameter <- .bt_label_formula_parameter(x)
    if(nzchar(formula_parameter) &&
       startsWith(name, paste0(formula_parameter, "_"))){
      return(list(.bt_label_parts_coefficient(
        parameter         = name,
        formula_parameter = formula_parameter,
        interaction_terms = attr(x, "interaction_terms", exact = TRUE)
      )))
    }
    return(list(.bt_label_parts(name, selector = name)))
  }
  columns <- colnames(x)
  if(is.null(columns)){
    columns <- if(ncol(x) == 1L) name else paste0(name, "[", seq_len(ncol(x)), "]")
  }
  formula_parameter <- .bt_label_formula_parameter(x)

  # factor draws: the level cells and coefficients their factor metadata
  # names (columns are matched to the rendered selectors, never parsed)
  factor_parts <- .bt_draws_factor_column_parts(x, name, columns, formula_parameter)
  if(!is.null(factor_parts)){
    return(factor_parts)
  }

  # columns named after the formula parameter keep their term text
  lapply(columns, function(column){
    if(nzchar(formula_parameter) &&
       startsWith(column, paste0(formula_parameter, "_"))){
      .bt_label_parts_coefficient(column, formula_parameter)
    }else{
      .bt_label_parts(column, selector = column)
    }
  })
}

# Label parts of the columns of factor draws without 'quantities' metadata
# (draws assembled outside the BayesTools producers): the level cells
# (as level effects or transformed contrast levels) and contrast coefficients
# of the factor metadata the draws carry, or of a factor prior in their
# 'prior_list', whose selectors are exactly the column names. NULL when no
# factor metadata names every column.
.bt_draws_factor_column_parts <- function(x, name, columns, formula_parameter){

  prior_list <- attr(x, "prior_list", exact = TRUE)
  if(is.prior(prior_list)){
    prior_list <- list(prior_list)
  }
  candidates <- c(list(x), if(is.list(prior_list)) prior_list)
  candidates <- Filter(function(candidate){
    !is.null(attr(candidate, "level_names", exact = TRUE)) &&
      (!is.prior(candidate) || is.prior.factor(candidate))
  }, candidates)

  for(candidate in candidates){
    # metadata that does not describe the factor's level cells names nothing
    parts <- tryCatch(
      .bt_label_parts_factor(name, candidate, formula_parameter),
      error = function(e) NULL
    )
    if(is.null(parts)){
      next
    }
    lookup <- c(
      parts$coordinates,
      parts$cells,
      .bt_label_parts_update(parts$cells, transformation = "dif")
    )
    rows <- match(columns, .bt_label(lookup, style = "selector"))
    if(!anyNA(rows)){
      return(lookup[rows])
    }
  }
  NULL
}

# Label parts of the log intercept of a formula shown as the exponentiated
# original-scale intercept: tables with transform_scaled = TRUE label the
# intercept of a log(intercept) formula exp(intercept).
.bt_label_parts_log_intercept <- function(parts, formula_scale){

  lapply(.bt_label_parts_list(parts), function(part){
    if(is.null(formula_scale) || !nzchar(part$formula_parameter) ||
       !identical(part$components, "intercept") ||
       length(part$levels) > 0L || !is.na(part$coefficient) ||
       !isTRUE(attr(formula_scale[[part$formula_parameter]], "log_intercept"))){
      return(part)
    }
    .bt_label_parts_update(part, transformation = "exp")[[1L]]
  })
}

# Draw-metadata column tables ('quantities') -------------------------------

.bt_draws_quantity_table <- function(columns, quantity_ids, dependencies,
                                     weights, label_parts){

  out <- data.frame(
    column      = as.character(columns),
    quantity_id = as.character(quantity_ids),
    stringsAsFactors = FALSE
  )
  out$dependencies <- I(lapply(dependencies, as.character))
  out$weights      <- I(lapply(weights, as.numeric))
  out$label_parts  <- I(label_parts)
  out
}

# The catalog quantity ids of label parts rendered to canonical names (""
# when the catalog has no such quantity).
.bt_draws_quantity_ids <- function(parts, catalog){

  out <- rep("", length(parts))
  if(is.null(catalog) || length(parts) == 0L){
    return(out)
  }
  quantities <- catalog$quantities
  quantities <- quantities[quantities$provider == "BayesTools", , drop = FALSE]
  selectors <- .bt_label(parts, style = "selector")
  namespaces <- vapply(parts, function(part){
    if(nzchar(part$formula_parameter)) part$formula_parameter else "model"
  }, character(1))
  rows <- match(
    paste(selectors, namespaces, sep = "\r"),
    paste(quantities$canonical_name, quantities$namespace, sep = "\r")
  )
  out[!is.na(rows)] <- quantities$quantity_id[rows[!is.na(rows)]]
  out
}

# The column table of mixed-posterior draws of one prior-list entry whose
# columns are the prior's fitted coordinates in coordinate order.
.bt_mixed_quantities <- function(parameter, prior, columns, catalog = NULL){

  prior_parts <- .bt_label_parts_prior(parameter, prior)
  if(length(prior_parts$coordinates) != length(columns)){
    stop(
      "The mixed posterior columns of '", parameter, "' do not match its ",
      "fitted coordinates.",
      call. = FALSE
    )
  }
  .bt_draws_quantity_table(
    columns      = columns,
    quantity_ids = .bt_draws_quantity_ids(prior_parts$parts, catalog),
    dependencies = as.list(prior_parts$coordinates),
    weights      = rep(list(1), length(columns)),
    label_parts  = prior_parts$parts
  )
}

# The column table of draws whose columns have no label structure beyond
# their names (weight-function and publication-bias columns): each column is
# the given fitted coordinate, or no coordinate when the column is a mixture
# of different coordinates across models.
.bt_verbatim_quantities <- function(columns, coordinates = NULL){

  .bt_draws_quantity_table(
    columns      = columns,
    quantity_ids = rep("", length(columns)),
    dependencies = if(is.null(coordinates)){
      rep(list(character()), length(columns))
    }else{
      as.list(coordinates)
    },
    weights      = if(is.null(coordinates)){
      rep(list(numeric()), length(columns))
    }else{
      rep(list(1), length(columns))
    },
    label_parts  = lapply(columns, function(column){
      .bt_label_parts(column, selector = column)
    })
  )
}

# The column table of draws 'x', aligned with its columns (one row for vector
# draws), or NULL when the draws carry none.
.bt_draws_quantities <- function(x){

  quantities <- .bt_meta_get(x, "quantities")
  if(is.null(quantities)){
    return(NULL)
  }
  columns <- if(is.null(dim(x))) NULL else colnames(x)
  n_columns <- if(is.null(dim(x))) 1L else ncol(x)
  if(nrow(quantities) != n_columns ||
     (!is.null(columns) && !identical(quantities$column, columns))){
    stop(
      "The draw metadata 'quantities' do not describe the columns of the draws.",
      call. = FALSE
    )
  }
  quantities
}

# The fitted coordinate of every column of draws 'x' (NA for columns that are
# not a single fitted coordinate), from its 'quantities' metadata. Vector
# draws without that metadata are the parameter 'parameter' itself.
.bt_draws_coordinate_columns <- function(x, parameter){

  quantities <- .bt_draws_quantities(x)
  if(is.null(quantities)){
    if(is.null(dim(x))){
      return(parameter)
    }
    stop(
      "The posterior samples of '", parameter, "' do not identify their fitted ",
      "coordinates (draw metadata 'quantities'). Create them with ",
      "as_mixed_posteriors() or mix_posteriors().",
      call. = FALSE
    )
  }
  vapply(seq_len(nrow(quantities)), function(i){
    dependencies <- quantities$dependencies[[i]]
    weights <- quantities$weights[[i]]
    if(length(dependencies) == 1L && isTRUE(weights == 1)){
      dependencies
    }else{
      NA_character_
    }
  }, character(1))
}

# Label parts of the levels of a marginal posterior of 'parameter' (a list of
# levels, or the draws themselves for a simple parameter): their 'quantities'
# draw metadata, otherwise the parameter's term followed by the level name.
.bt_marginal_level_parts <- function(x, parameter){

  levels <- if(is.list(x) && !is.numeric(x)) x else list(x)
  level_names <- names(levels)
  if(is.null(level_names)){
    level_names <- rep("", length(levels))
  }
  formula_parameter <- .bt_label_formula_parameter(x)
  lapply(seq_along(levels), function(i){
    quantities <- .bt_draws_quantities(levels[[i]])
    if(!is.null(quantities)){
      return(quantities$label_parts[[1L]])
    }
    term <- if(nzchar(formula_parameter) &&
               startsWith(parameter, paste0(formula_parameter, "_"))){
      paste(
        .bt_label_parts_coefficient(parameter, formula_parameter)$components,
        collapse = ":"
      )
    }else{
      parameter
    }
    level <- if(nzchar(level_names[[i]]) && !identical(level_names[[i]], "intercept")){
      paste0("[", level_names[[i]], "]")
    }else{
      ""
    }
    .bt_label_parts(
      components        = paste0(term, level),
      formula_parameter = formula_parameter,
      selector          = paste0(parameter, level)
    )
  })
}

# Label parts of summary columns ------------------------------------------------

# Label parts of the columns of a summary built from a fit: the columns the
# summary created ('created', keyed by column: inclusion rows, mixture
# components, transformed factor levels), fitted coordinates (the coordinate
# label parts of the parameter map; of the prior list for objects without
# one), catalog quantities (random-effect summaries), and otherwise the column
# name itself.
.bt_estimates_column_parts <- function(columns, fit, created = list(),
                                       coordinates = NULL){

  out <- vector("list", length(columns))
  from_created <- columns %in% names(created)
  out[from_created] <- created[columns[from_created]]
  has_map <- !is.null(attr(fit, "parameter_map", exact = TRUE))
  prior_list <- attr(fit, "prior_list", exact = TRUE)

  if(has_map){
    if(is.null(coordinates)){
      coordinates <- parameter_coordinates(fit)
    }
    coordinate_rows <- match(columns, coordinates$coordinate_name)
    from_coordinates <- !from_created & !is.na(coordinate_rows)
    if(any(from_coordinates)){
      out[from_coordinates] <- .bt_parameter_coordinates_row_label_parts(
        coordinates    = coordinates[coordinate_rows[from_coordinates], , drop = FALSE],
        prior_list     = prior_list,
        formula_design = attr(fit, "formula_design", exact = TRUE)
      )
    }
  }else if(length(prior_list) > 0L){
    for(parameter in names(prior_list)){
      prior_parts <- .bt_label_parts_prior(parameter, prior_list[[parameter]])
      rows <- match(prior_parts$coordinates, columns)
      free <- !is.na(rows)
      free[free] <- vapply(out[rows[free]], is.null, logical(1))
      out[rows[free]] <- prior_parts$parts[free]
    }
  }

  remaining <- vapply(out, is.null, logical(1))
  if(any(remaining) && has_map){
    quantities <- parameter_catalog(fit)$quantities
    quantities <- quantities[quantities$provider == "BayesTools", , drop = FALSE]
    catalog_rows <- match(columns, quantities$canonical_name)
    from_catalog <- remaining & !is.na(catalog_rows)
    out[from_catalog] <- unclass(quantities$label_parts)[catalog_rows[from_catalog]]
    remaining <- remaining & !from_catalog
  }
  out[remaining] <- lapply(columns[remaining], function(column){
    .bt_label_parts(column, selector = column)
  })

  out
}
