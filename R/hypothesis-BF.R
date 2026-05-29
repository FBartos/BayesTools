# ============================================================================ #
# hypothesis-BF.R
# ============================================================================ #
#
# Generic Bayes factors for scalar prior/posterior quantities.
#
# ============================================================================ #


#' @title Hypothesis Bayes Factors
#'
#' @description Computes Bayes factors for scalar hypotheses written as
#' expressions, for example \code{"theta = 0"}, \code{"theta > 0"}, or
#' \code{"theta = 0 vs theta > 0"}. Factor/posterior-list levels can be
#' referenced as \code{"mu_alloc[alternate] > mu_alloc[random]"} when the
#' marginal posterior carries joint prior information. Point-null hypotheses use
#' \code{\link{Savage_Dickey_BF}} when a \code{marginal_posterior} object is
#' supplied, including precomputed qCMDE/IWMDE posterior ordinates when
#' \code{density_method = "precomputed"}. Level-specific
#' \code{marginal_posterior.*} subclass vectors with attached density
#' attributes are accepted as marginal posterior inputs.
#'
#' @param posterior posterior draws, a \code{marginal_posterior}, a
#' \code{marginal_inference} object, or a data frame/matrix of posterior draws.
#' @param prior prior draws for numeric/data-frame inputs. Ignored when
#' \code{posterior} already contains deterministic prior density information.
#' @param hypothesis character vector with scalar hypothesis statements.
#' @param parameter optional scalar quantity name for numeric vectors or
#' \code{marginal_posterior} objects.
#' @param logBF whether to display the Bayes factor on the log scale.
#' @param BF01 whether to display the inverse Bayes factor.
#' @param seed optional seed used only by downstream helpers that sample.
#' @param density_method posterior density source for point-null tests.
#' \code{"KDE"} uses kernel density estimates. \code{"precomputed"} uses valid
#' \code{posterior_ordinate} or \code{posterior_density} attributes before
#' falling back to KDE.
#' @param columns output columns. \code{"default"} returns \code{Alternative},
#' \code{Null}, \code{BF}, and \code{BF_error}. \code{"all"} also returns
#' \code{prior}, \code{posterior}, and \code{method} columns. The
#' \code{prior} and \code{posterior} columns are diagnostics, not always
#' probabilities: point tests report density heights, region tests report odds,
#' and transitive point-vs-region tests report \code{NA}. The default
#' \code{BF_error} column is printed as \code{error\%(BF)}.
#' @param ... unused.
#'
#' @return A BayesTools table of class \code{BayesTools_hypothesis_BF}. The
#' \code{BF_error} column reports approximate relative Monte Carlo error
#' percentage when available. Region odds errors are computed on
#' \code{log(BF)} from prior/posterior region indicators; point-vs-region
#' errors combine the available point-density and region-mass errors on the
#' \code{log(BF)} scale.
#'
#' @export
hypothesis_BF <- function(posterior, prior = NULL, hypothesis, parameter = NULL,
                          logBF = FALSE, BF01 = FALSE, seed = NULL,
                          density_method = c("KDE", "precomputed"),
                          columns = "default", ...) {

  check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  check_char(parameter, "parameter", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01, "BF01", allow_NA = FALSE)
  check_real(seed, "seed", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  density_method <- .posterior_density_method(density_method)
  columns        <- .hypothesis_BF_output_columns(columns)

  if(!is.null(seed)){
    set.seed(seed)
  }

  parsed <- lapply(hypothesis, .parse_hypothesis_BF)
  quantities <- .as_hypothesis_quantities(
    posterior  = posterior,
    prior      = prior,
    parsed     = parsed,
    parameter  = parameter
  )

  rows <- list()
  row_i <- 1L
  for(hyp_i in seq_along(parsed)){
    for(quantity_i in seq_along(quantities)){
      result <- .hypothesis_BF_compute(
        quantity       = quantities[[quantity_i]],
        parsed         = parsed[[hyp_i]],
        density_method = density_method
      )
      rows[[row_i]] <- .hypothesis_BF_row(
        quantity = quantities[[quantity_i]],
        parsed   = parsed[[hyp_i]],
        result   = result
      )
      row_i <- row_i + 1L
    }
  }

  out <- do.call(rbind, rows)
  raw_BF <- out[["BF"]]
  out[["BF"]] <- format_BF(raw_BF, logBF = logBF, BF01 = BF01)
  attr(out[["BF_error"]], "name") <- "error%(BF)"

  warnings <- .hypothesis_BF_table_warnings(out)
  out      <- out[, columns, drop = FALSE]

  attr(out, "raw_BF")   <- raw_BF
  attr(out, "parsed")   <- parsed
  attr(out, "logBF")    <- logBF
  attr(out, "BF01")     <- BF01
  attr(out, "type")      <- .hypothesis_BF_table_types(colnames(out))
  attr(out, "footnotes") <- .hypothesis_BF_table_footnotes(colnames(out))
  attr(out, "warnings")  <- warnings
  attr(out, "rownames")  <- TRUE
  class(out) <- c("BayesTools_table", "BayesTools_hypothesis_BF", "data.frame")

  return(out)
}


#' @export
print.BayesTools_hypothesis_BF <- function(x, ...) {

  if(!inherits(x, "BayesTools_table")){
    class(x) <- c("BayesTools_table", class(x))
  }
  print(x, ...)

  return(invisible(x))
}


.parse_hypothesis_BF <- function(hypothesis) {

  hypothesis <- .hypothesis_normalize_level_references(hypothesis)
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


.hypothesis_normalize_level_references <- function(text) {

  pieces <- strsplit(text, "`", fixed = TRUE)[[1L]]
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
    }
  }

  paste(pieces, collapse = "`")
}


.hypothesis_split_vs <- function(hypothesis) {

  parts <- strsplit(hypothesis, "\\s+[Vv][Ss]\\s+", perl = TRUE)[[1L]]
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

  value <- eval(expr, envir = baseenv())
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


.as_hypothesis_quantities <- function(posterior, prior, parsed, parameter) {

  if(inherits(posterior, "marginal_inference")){
    return(.as_hypothesis_quantities_marginal_inference(
      posterior = posterior,
      parsed    = parsed,
      parameter = parameter
    ))
  }

  if(.hypothesis_inherits_marginal_posterior(posterior)){
    return(.as_hypothesis_quantities_marginal_posterior(
      posterior = posterior,
      parsed    = parsed,
      parameter = parameter
    ))
  }

  if(is.data.frame(posterior) || is.matrix(posterior)){
    if(is.null(prior)){
      stop("Prior draws are required for data-frame posterior inputs.",
           call. = FALSE)
    }
    return(list(.hypothesis_quantity_from_draws(
      posterior = as.data.frame(posterior, check.names = FALSE),
      prior     = as.data.frame(prior, check.names = FALSE),
      label     = "draws"
    )))
  }

  if(is.numeric(posterior)){
    if(is.null(prior)){
      stop("Prior draws are required for numeric posterior inputs.",
           call. = FALSE)
    }
    check_real(posterior, "posterior", check_length = 0, allow_NA = FALSE)
    if(is.null(parameter)){
      parameter <- .hypothesis_single_symbol(parsed)
    }
    if(is.null(parameter)){
      stop("The 'parameter' argument is required for numeric draws unless ",
           "the hypothesis contains exactly one named quantity.", call. = FALSE)
    }

    prior_info <- .hypothesis_prior_input_to_draws(
      prior     = prior,
      parameter = parameter,
      n         = max(length(posterior), 10000L)
    )
    posterior_df <- data.frame(posterior, check.names = FALSE)
    prior_df     <- data.frame(prior_info[["draws"]], check.names = FALSE)
    names(posterior_df) <- parameter
    names(prior_df)     <- parameter

    return(list(.hypothesis_quantity_from_draws(
      posterior     = posterior_df,
      prior         = prior_df,
      label         = parameter,
      parameter     = parameter,
      prior_density = prior_info[["density"]]
    )))
  }

  stop("Unsupported posterior input for hypothesis_BF().", call. = FALSE)
}


.hypothesis_inherits_marginal_posterior <- function(x) {

  inherits(x, "marginal_posterior") ||
    any(grepl("^marginal_posterior\\.", class(x)))
}


.as_hypothesis_quantities_marginal_inference <- function(posterior, parsed,
                                                         parameter) {

  if(is.null(posterior[["conditional"]])){
    stop("'marginal_inference' object does not contain conditional draws.",
         call. = FALSE)
  }

  available <- names(posterior[["conditional"]])
  if(is.null(parameter)){
    level_refs <- .hypothesis_level_references(parsed)
    symbols <- .hypothesis_all_symbols(parsed)
    matches <- intersect(symbols, available)
    if(nrow(level_refs) > 0L &&
       length(unique(level_refs[["parameter"]])) == 1L &&
       unique(level_refs[["parameter"]]) %in% available){
      parameter <- unique(level_refs[["parameter"]])
    }else if(length(matches) == 1L){
      parameter <- matches
    }else if(length(matches) == 0L && length(available) == 1L){
      parameter <- available
    }else{
      stop("Specify 'parameter' for this marginal_inference object.",
           call. = FALSE)
    }
  }
  if(!parameter %in% available){
    stop("Parameter '", parameter, "' is not available in 'posterior'.",
         call. = FALSE)
  }

  return(.as_hypothesis_quantities_marginal_posterior(
    posterior = posterior[["conditional"]][[parameter]],
    parsed    = parsed,
    parameter = parameter
  ))
}


.as_hypothesis_quantities_marginal_posterior <- function(posterior, parsed,
                                                         parameter) {

  level_refs <- .hypothesis_level_references(parsed, parameter)
  symbol_parameter <- .hypothesis_single_symbol(parsed)
  if(is.null(parameter) && nrow(level_refs) > 0L &&
     length(unique(level_refs[["parameter"]])) == 1L){
    parameter <- unique(level_refs[["parameter"]])
  }
  if(is.null(parameter)){
    parameter <- symbol_parameter
  }
  if(is.null(parameter)){
    parameter <- attr(posterior, "parameter", exact = TRUE)
  }
  if(is.null(parameter)){
    stop("The 'parameter' argument is required for this marginal posterior.",
         call. = FALSE)
  }

  if(is.list(posterior)){
    level_refs <- .hypothesis_level_references(parsed, parameter)
    if(nrow(level_refs) > 0L){
      return(list(.hypothesis_quantity_from_marginal_posterior_levels(
        posterior  = posterior,
        parameter  = parameter,
        level_refs = level_refs,
        parsed     = parsed
      )))
    }

    out <- lapply(seq_along(posterior), function(i){
      level <- names(posterior)[i]
      if(is.null(level) || !nzchar(level)){
        level <- attr(posterior[[i]], "level", exact = TRUE)
      }
      if(is.null(level) || !nzchar(level)){
        level <- attr(posterior[[i]], "level_name", exact = TRUE)
      }
      label <- if(is.null(level) || !nzchar(level)){
        parameter
      }else{
        paste0(parameter, "[", level, "]")
      }
      .hypothesis_quantity_from_marginal_posterior(
        posterior = posterior[[i]],
        parameter = parameter,
        label     = label
      )
    })
    return(out)
  }

  return(list(.hypothesis_quantity_from_marginal_posterior(
    posterior = posterior,
    parameter = parameter,
    label     = parameter
  )))
}


.hypothesis_level_references <- function(parsed, parameter = NULL) {

  symbols <- .hypothesis_all_symbols(parsed)
  refs <- regexec("^([^\\[]+)\\[([^\\]]+)\\]$", symbols, perl = TRUE)
  matches <- regmatches(symbols, refs)
  has_match <- vapply(matches, length, integer(1)) == 3L
  if(!any(has_match)){
    return(data.frame(
      symbol    = character(),
      parameter = character(),
      level     = character(),
      stringsAsFactors = FALSE
    ))
  }

  matches <- matches[has_match]
  out <- data.frame(
    symbol    = vapply(matches, `[[`, character(1), 1L),
    parameter = vapply(matches, `[[`, character(1), 2L),
    level     = trimws(vapply(matches, `[[`, character(1), 3L)),
    stringsAsFactors = FALSE
  )

  if(!is.null(parameter)){
    out <- out[out[["parameter"]] == parameter, , drop = FALSE]
  }

  return(unique(out))
}


.hypothesis_prior_input_to_draws <- function(prior, parameter, n) {

  prior_density <- NULL
  if(inherits(prior, "prior")){
    prior_density <- tryCatch(
      .prior_linear_combination_density(
        prior_list = setNames(list(prior), parameter),
        weights    = setNames(1, parameter),
        n_grid     = .prior_linear_density_default_grid()
      ),
      error = function(e) NULL
    )
    prior <- rng(prior, n)
    if(is.matrix(prior) || is.data.frame(prior)){
      if(ncol(prior) != 1L){
        stop("Prior object must generate scalar draws for numeric posterior inputs.",
             call. = FALSE)
      }
      prior <- prior[, 1L]
    }
  }

  check_real(prior, "prior", check_length = 0, allow_NA = FALSE)

  return(list(
    draws   = prior,
    density = prior_density
  ))
}


.hypothesis_quantity_from_draws <- function(posterior, prior, label,
                                            parameter = NULL,
                                            prior_density = NULL) {

  out <- list(
    label              = label,
    parameter          = parameter,
    posterior_draws    = posterior,
    prior_draws        = prior,
    posterior_marginal = NULL,
    posterior_marginals = NULL,
    prior_densities    = NULL,
    prior_density      = prior_density
  )
  class(out) <- "BayesTools_hypothesis_quantity"

  return(out)
}


.hypothesis_quantity_from_marginal_posterior_levels <- function(posterior,
                                                                parameter,
                                                                level_refs,
                                                                parsed) {

  levels <- unique(level_refs[["level"]])
  available <- names(posterior)
  missing <- setdiff(levels, available)
  if(length(missing) > 0L){
    stop("Hypothesis references unknown level '",
         paste(missing, collapse = "', '"), "' for parameter '", parameter, "'.",
         call. = FALSE)
  }

  posterior_draws <- lapply(levels, function(level) as.numeric(posterior[[level]]))
  n_draws <- vapply(posterior_draws, length, integer(1))
  if(length(unique(n_draws)) != 1L){
    stop("Level comparisons require equal-length posterior draws.",
         call. = FALSE)
  }

  .hypothesis_validate_level_conditionals(posterior, parameter, levels)

  posterior_df <- as.data.frame(posterior_draws, check.names = FALSE)
  names(posterior_df) <- paste0(parameter, "[", levels, "]")

  prior_df <- NULL
  if(.hypothesis_level_hypotheses_need_prior_draws(parsed, level_refs)){
    prior_df <- .hypothesis_prior_draws_from_marginal_levels(
      posterior = posterior,
      parameter = parameter,
      levels    = levels,
      n         = max(nrow(posterior_df), 10000L)
    )
  }

  posterior_marginals <- lapply(levels, function(level) {
    .hypothesis_marginal_child(posterior[[level]])
  })
  names(posterior_marginals) <- names(posterior_df)
  prior_densities <- lapply(levels, function(level) {
    attr(posterior[[level]], "prior_density", exact = TRUE)
  })
  names(prior_densities) <- names(posterior_df)

  out <- list(
    label              = parameter,
    parameter          = parameter,
    posterior_draws    = posterior_df,
    prior_draws        = prior_df,
    posterior_marginal = NULL,
    posterior_marginals = posterior_marginals,
    prior_densities    = prior_densities,
    prior_density      = NULL
  )
  class(out) <- "BayesTools_hypothesis_quantity"

  return(out)
}


.hypothesis_marginal_child <- function(posterior) {

  if(.hypothesis_inherits_marginal_posterior(posterior) &&
     !inherits(posterior, "marginal_posterior")){
    class(posterior) <- unique(c(class(posterior), "marginal_posterior"))
  }

  return(posterior)
}


.hypothesis_validate_level_conditionals <- function(posterior, parameter,
                                                    levels) {

  conditionals <- lapply(levels, function(level) {
    attr(posterior[[level]], "effective_conditional", exact = TRUE)
  })
  has_conditional <- vapply(conditionals, function(x) !is.null(x), logical(1))
  if(!any(has_conditional)){
    return(invisible(TRUE))
  }

  keys <- vapply(conditionals, .hypothesis_conditional_key, character(1))
  if(length(unique(keys)) > 1L){
    stop(
      "Level comparison for parameter '", parameter,
      "' uses different conditional posterior subsets. Use averaged marginals ",
      "or compare levels with identical conditionals.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


.hypothesis_conditional_key <- function(conditional) {

  if(is.null(conditional)){
    return("<averaged>")
  }

  paste(sort(unique(as.character(conditional))), collapse = "\r")
}


.hypothesis_level_hypotheses_need_prior_draws <- function(parsed, level_refs) {

  valid_symbols <- unique(level_refs[["symbol"]])
  for(hypothesis in parsed){
    if(!.hypothesis_sides_point_complement(hypothesis[["left"]],
                                           hypothesis[["right"]]) &&
       !.hypothesis_sides_point_complement(hypothesis[["right"]],
                                           hypothesis[["left"]])){
      return(TRUE)
    }

    point_side <- if(identical(hypothesis[["left"]][["type"]], "point")){
      hypothesis[["left"]]
    }else{
      hypothesis[["right"]]
    }
    symbol <- .hypothesis_direct_symbol(point_side[["expr"]])
    if(is.null(symbol) || !symbol %in% valid_symbols){
      return(TRUE)
    }
  }

  return(FALSE)
}


.hypothesis_prior_draws_from_marginal_levels <- function(posterior, parameter,
                                                         levels, n) {

  attr_prior_draws <- attr(posterior, "prior_draws", exact = TRUE)
  columns <- paste0(parameter, "[", levels, "]")
  if(!is.null(attr_prior_draws)){
    attr_prior_draws <- as.data.frame(attr_prior_draws, check.names = FALSE)
    missing <- setdiff(columns, names(attr_prior_draws))
    if(length(missing) == 0L){
      return(attr_prior_draws[, columns, drop = FALSE])
    }
  }

  child_prior_draws <- lapply(levels, function(level) {
    attr(posterior[[level]], "prior_draws", exact = TRUE)
  })
  if(all(!vapply(child_prior_draws, is.null, logical(1)))){
    n_draws <- vapply(child_prior_draws, length, integer(1))
    if(length(unique(n_draws)) != 1L){
      stop("Level comparisons require equal-length prior draws.",
           call. = FALSE)
    }
    out <- as.data.frame(child_prior_draws, check.names = FALSE)
    names(out) <- columns
    return(out)
  }

  context <- attr(posterior, "prior_density_context", exact = TRUE)
  if(is.null(context)){
    stop("Joint prior information is required for level-comparison hypotheses.",
         call. = FALSE)
  }
  if(!.hypothesis_is_prior_density_context(context)){
    stop("Invalid joint prior information for level-comparison hypotheses.",
         call. = FALSE)
  }

  prior_matrix <- .hypothesis_prior_samples_from_context(context, n)
  level_weights <- lapply(levels, function(level) {
    weights <- attr(posterior[[level]], "linear_weights", exact = TRUE)
    if(is.null(weights)){
      stop("Linear prior weights are missing for level '", level, "'.",
           call. = FALSE)
    }
    .hypothesis_prepare_level_weights(weights)
  })
  names(level_weights) <- levels
  row_i <- .hypothesis_level_weight_rows(level_weights, nrow(prior_matrix))

  prior_values <- lapply(levels, function(level){
    .hypothesis_apply_level_weights(prior_matrix, level_weights[[level]], row_i)
  })

  out <- as.data.frame(prior_values, check.names = FALSE)
  names(out) <- columns

  return(out)
}


.hypothesis_is_prior_density_context <- function(context) {

  inherits(context, "prior_density_context") ||
    inherits(context, "prior_density_model_mixture_context") ||
    inherits(context, "prior_density_conditional_context")
}


.hypothesis_prior_samples_from_context <- function(context, n) {

  if(inherits(context, "prior_density_context")){
    samples <- .generate_transformed_prior_samples(
      prior_list    = context[["prior_list"]],
      column_names  = context[["column_names"]],
      n_samples     = n,
      formula_scale = context[["formula_scale"]]
    )
    return(.hypothesis_complete_prior_sample_matrix(
      samples      = samples,
      column_names = context[["column_names"]],
      n            = n
    ))
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    return(.hypothesis_prior_samples_from_model_mixture_context(context, n))
  }

  if(inherits(context, "prior_density_conditional_context")){
    return(.hypothesis_prior_samples_from_conditional_context(context, n))
  }

  stop("Unknown prior density context.", call. = FALSE)
}


.hypothesis_prior_samples_from_model_mixture_context <- function(context, n) {

  model_i <- sample(seq_along(context[["model_weights"]]), size = n,
                    replace = TRUE, prob = context[["model_weights"]])
  out <- .hypothesis_empty_prior_sample_matrix(context[["column_names"]], n)

  for(i in unique(model_i)){
    rows <- which(model_i == i)
    model_prior_list <- .hypothesis_model_mixture_prior_list(context, i)
    samples <- .generate_transformed_prior_samples(
      prior_list   = model_prior_list,
      column_names = context[["column_names"]],
      n_samples    = length(rows)
    )
    out[rows, ] <- .hypothesis_complete_prior_sample_matrix(
      samples      = samples,
      column_names = context[["column_names"]],
      n            = length(rows)
    )
  }

  return(out)
}


.hypothesis_model_mixture_prior_list <- function(context, model_i) {

  model_prior_list <- lapply(context[["prior_list"]], function(parameter_priors) {
    if(is.prior(parameter_priors)){
      return(parameter_priors)
    }
    parameter_priors[[model_i]]
  })
  names(model_prior_list) <- names(context[["prior_list"]])

  for(parameter in names(model_prior_list)){
    if(is.null(model_prior_list[[parameter]])){
      model_prior_list[[parameter]] <- prior("point", list(location = 0))
    }
  }

  return(model_prior_list)
}


.hypothesis_prior_samples_from_conditional_context <- function(context, n) {

  out <- .hypothesis_empty_prior_sample_matrix(context[["column_names"]], n)
  if(length(context[["prior_lists"]]) == 0L){
    return(out)
  }

  model_i <- sample(seq_along(context[["model_weights"]]), size = n,
                    replace = TRUE, prob = context[["model_weights"]])
  for(i in unique(model_i)){
    rows <- which(model_i == i)
    samples <- .generate_transformed_prior_samples(
      prior_list    = context[["prior_lists"]][[i]],
      column_names  = context[["column_names"]],
      n_samples     = length(rows),
      formula_scale = context[["formula_scale"]]
    )
    out[rows, ] <- .hypothesis_complete_prior_sample_matrix(
      samples      = samples,
      column_names = context[["column_names"]],
      n            = length(rows)
    )
  }

  return(out)
}


.hypothesis_empty_prior_sample_matrix <- function(column_names, n) {

  matrix(
    0,
    nrow = n,
    ncol = length(column_names),
    dimnames = list(NULL, column_names)
  )
}


.hypothesis_complete_prior_sample_matrix <- function(samples, column_names, n) {

  samples <- as.matrix(samples)
  out <- .hypothesis_empty_prior_sample_matrix(column_names, n)
  columns <- intersect(column_names, colnames(samples))
  if(length(columns) > 0L){
    out[, columns] <- samples[, columns, drop = FALSE]
  }

  return(out)
}


.hypothesis_prepare_level_weights <- function(weights) {

  if(is.null(dim(weights))){
    weight_names <- names(weights)
    weights <- matrix(weights, nrow = 1L)
    colnames(weights) <- weight_names
  }else{
    weights <- as.matrix(weights)
  }
  if(is.null(colnames(weights))){
    stop("Linear prior weights must be named.", call. = FALSE)
  }

  weights
}


.hypothesis_level_weight_rows <- function(level_weights, n) {

  row_counts <- vapply(level_weights, nrow, integer(1))
  row_counts <- row_counts[row_counts > 1L]
  if(length(row_counts) == 0L){
    return(NULL)
  }
  if(length(unique(row_counts)) != 1L){
    stop("Level comparisons with row-varying prior weights require matching ",
         "weight rows across referenced levels.", call. = FALSE)
  }

  sample.int(row_counts[[1L]], size = n, replace = TRUE)
}


.hypothesis_apply_level_weights <- function(samples, weights, row_i = NULL) {

  samples <- as.matrix(samples)
  weights <- .hypothesis_prepare_level_weights(weights)
  nonzero_columns <- colnames(weights)[colSums(abs(weights), na.rm = TRUE) >
                                          .prior_linear_density_zero_tol()]
  missing <- setdiff(nonzero_columns, colnames(samples))
  if(length(missing) > 0L){
    stop("Linear prior weights reference columns not available in the joint ",
         "prior context: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }

  columns <- intersect(colnames(weights), colnames(samples))
  if(length(columns) == 0L){
    return(rep(0, nrow(samples)))
  }
  if(nrow(weights) == 1L){
    return(as.numeric(samples[, columns, drop = FALSE] %*%
                        as.numeric(weights[1L, columns])))
  }

  if(is.null(row_i)){
    row_i <- sample.int(nrow(weights), size = nrow(samples), replace = TRUE)
  }
  if(any(row_i > nrow(weights))){
    stop("Shared level-weight rows exceed available weight rows.",
         call. = FALSE)
  }

  rowSums(samples[, columns, drop = FALSE] *
            weights[row_i, columns, drop = FALSE])
}


.hypothesis_quantity_from_marginal_posterior <- function(posterior, parameter,
                                                         label) {

  if(.hypothesis_inherits_marginal_posterior(posterior) &&
     !inherits(posterior, "marginal_posterior")){
    class(posterior) <- unique(c(class(posterior), "marginal_posterior"))
  }

  posterior_df <- data.frame(as.numeric(posterior), check.names = FALSE)
  names(posterior_df) <- parameter

  out <- list(
    label              = label,
    parameter          = parameter,
    posterior_draws    = posterior_df,
    prior_draws        = NULL,
    posterior_marginal = posterior,
    posterior_marginals = NULL,
    prior_densities    = NULL,
    prior_density      = attr(posterior, "prior_density", exact = TRUE)
  )
  class(out) <- "BayesTools_hypothesis_quantity"

  return(out)
}


.hypothesis_BF_compute <- function(quantity, parsed, density_method) {

  left  <- parsed[["left"]]
  right <- parsed[["right"]]
  explicit <- isTRUE(parsed[["explicit"]])

  if(.hypothesis_sides_point_complement(left, right)){
    if(explicit){
      return(.hypothesis_BF_result_labels(.hypothesis_point_BF(
        quantity       = quantity,
        side           = left,
        density_method = density_method,
        inverse        = FALSE
      ), alternative = left, null = right))
    }
    return(.hypothesis_BF_result_labels(.hypothesis_point_BF(
      quantity       = quantity,
      side           = left,
      density_method = density_method,
      inverse        = TRUE
    ), alternative = right, null = left))
  }
  if(.hypothesis_sides_point_complement(right, left)){
    return(.hypothesis_BF_result_labels(.hypothesis_point_BF(
      quantity       = quantity,
      side           = right,
      density_method = density_method,
      inverse        = TRUE
    ), alternative = left, null = right))
  }

  if(identical(left[["type"]], "region") &&
     identical(right[["type"]], "region")){
    return(.hypothesis_BF_result_labels(
      .hypothesis_region_odds_BF(quantity, left, right),
      alternative = left,
      null        = right
    ))
  }

  if(identical(left[["type"]], "point") &&
     identical(right[["type"]], "region")){
    return(.hypothesis_BF_result_labels(.hypothesis_transitive_BF(
      quantity       = quantity,
      point_side     = left,
      region_side    = right,
      density_method = density_method,
      inverse        = FALSE
    ), alternative = left, null = right))
  }
  if(identical(left[["type"]], "region") &&
     identical(right[["type"]], "point")){
    return(.hypothesis_BF_result_labels(.hypothesis_transitive_BF(
      quantity       = quantity,
      point_side     = right,
      region_side    = left,
      density_method = density_method,
      inverse        = TRUE
    ), alternative = left, null = right))
  }

  stop("Unsupported hypothesis comparison.", call. = FALSE)
}


.hypothesis_BF_result_labels <- function(result, alternative, null) {

  result[["alternative"]] <- alternative[["label"]]
  result[["null"]]        <- null[["label"]]

  return(result)
}


.hypothesis_sides_point_complement <- function(point_side, other_side) {

  if(!identical(point_side[["type"]], "point") ||
     !identical(other_side[["type"]], "not_point")){
    return(FALSE)
  }

  identical(point_side[["expr"]], other_side[["expr"]]) &&
    isTRUE(all.equal(point_side[["value"]], other_side[["value"]]))
}


.hypothesis_point_BF <- function(quantity, side, density_method,
                                 inverse = FALSE) {

  marginal <- .hypothesis_point_marginal(quantity, side)
  if(!is.null(marginal)){
    inclusion_BF <- Savage_Dickey_BF(
      posterior            = marginal[["posterior"]],
      null_hypothesis      = side[["value"]],
      normal_approximation = FALSE,
      silent               = TRUE,
      density_method       = density_method
    )
    prior_value <- .hypothesis_prior_density_height(
      marginal[["prior_density"]],
      side[["value"]]
    )
    posterior_value <- prior_value / as.numeric(inclusion_BF)
    BF <- posterior_value / prior_value
    warning <- attr(inclusion_BF, "warnings", exact = TRUE)
    BF_error <- attr(inclusion_BF, "BF_error_percent", exact = TRUE)
    method <- if(identical(density_method, "precomputed")){
      "Savage-Dickey (precomputed)"
    }else{
      "Savage-Dickey"
    }
  }else{
    posterior <- .hypothesis_eval_expression(
      side[["expr"]],
      quantity[["posterior_draws"]]
    )
    if(!is.null(quantity[["prior_density"]]) &&
       .hypothesis_expression_is_parameter(side[["expr"]], quantity[["parameter"]])){
      prior_value <- .hypothesis_prior_density_height(
        quantity[["prior_density"]],
        side[["value"]]
      )
    }else{
      prior <- .hypothesis_eval_expression(
        side[["expr"]],
        .hypothesis_prior_draws(quantity)
      )
      prior_value <- .hypothesis_sample_density_height(prior, side[["value"]],
                                                       "prior")
    }
    posterior_value <- .hypothesis_sample_density_height(
      posterior, side[["value"]], "posterior"
    )
    BF <- posterior_value / prior_value
    warning <- NULL
    BF_error <- NA_real_
    method <- "kernel Savage-Dickey"
  }

  if(inverse){
    BF <- 1 / as.numeric(BF)
  }

  return(list(
    BF        = as.numeric(BF),
    prior     = prior_value,
    posterior = posterior_value,
    method    = method,
    BF_error  = if(is.null(BF_error)) NA_real_ else as.numeric(BF_error),
    warning   = .hypothesis_collapse_warning(warning)
  ))
}


.hypothesis_point_marginal <- function(quantity, side) {

  if(!is.null(quantity[["posterior_marginal"]]) &&
     .hypothesis_expression_is_parameter(side[["expr"]],
                                         quantity[["parameter"]])){
    return(list(
      posterior     = quantity[["posterior_marginal"]],
      prior_density = quantity[["prior_density"]]
    ))
  }

  symbol <- .hypothesis_direct_symbol(side[["expr"]])
  if(is.null(symbol) || is.null(quantity[["posterior_marginals"]]) ||
     !symbol %in% names(quantity[["posterior_marginals"]])){
    return(NULL)
  }

  return(list(
    posterior     = quantity[["posterior_marginals"]][[symbol]],
    prior_density = quantity[["prior_densities"]][[symbol]]
  ))
}


.hypothesis_region_odds_BF <- function(quantity, left, right) {

  prior_left      <- .hypothesis_region_mass(quantity, left, prior = TRUE)
  prior_right     <- .hypothesis_region_mass(quantity, right, prior = TRUE)
  posterior_left  <- .hypothesis_region_mass(quantity, left, prior = FALSE)
  posterior_right <- .hypothesis_region_mass(quantity, right, prior = FALSE)

  .hypothesis_check_prior_mass(prior_left, left[["label"]])
  .hypothesis_check_prior_mass(prior_right, right[["label"]])

  BF       <- (posterior_left / posterior_right) / (prior_left / prior_right)
  BF_error <- .hypothesis_region_odds_BF_error_percent(quantity, left, right)
  warning  <- NULL
  if(posterior_left == 0 || posterior_right == 0){
    warning <- "Posterior region mass is zero; reported BF is boundary-valued."
  }

  return(list(
    BF        = BF,
    prior     = prior_left / prior_right,
    posterior = posterior_left / posterior_right,
    method    = "prior-posterior odds",
    BF_error  = BF_error,
    warning   = warning
  ))
}


.hypothesis_transitive_BF <- function(quantity, point_side, region_side,
                                      density_method, inverse) {

  if(!.hypothesis_point_region_compatible(point_side, region_side)){
    stop("Point-vs-region hypotheses must use the same scalar expression.",
         call. = FALSE)
  }

  point_BF <- .hypothesis_point_BF(
    quantity       = quantity,
    side           = point_side,
    density_method = density_method,
    inverse        = FALSE
  )
  region_prior     <- .hypothesis_region_mass(quantity, region_side, prior = TRUE)
  region_posterior <- .hypothesis_region_mass(quantity, region_side, prior = FALSE)
  .hypothesis_check_prior_mass(region_prior, region_side[["label"]])

  region_BF       <- region_posterior / region_prior
  region_BF_error <- .hypothesis_region_BF_error_percent(quantity, region_side)
  BF              <- point_BF[["BF"]] / region_BF
  if(inverse){
    BF <- 1 / BF
  }

  warning <- point_BF[["warning"]]
  if(region_posterior == 0){
    warning <- .hypothesis_collapse_warning(c(
      warning,
      "Posterior region mass is zero; reported BF is boundary-valued."
    ))
  }

  return(list(
    BF        = BF,
    prior     = NA_real_,
    posterior = NA_real_,
    method    = "transitive Savage-Dickey",
    BF_error  = .hypothesis_combine_BF_error_percent(
      point_BF[["BF_error"]],
      region_BF_error
    ),
    warning   = warning
  ))
}


.hypothesis_region_odds_BF_error_percent <- function(quantity, left, right) {

  posterior_var <- .hypothesis_region_log_odds_mc_var(
    quantity = quantity,
    left     = left,
    right    = right,
    prior    = FALSE
  )
  prior_var     <- .hypothesis_region_log_odds_mc_var(
    quantity = quantity,
    left     = left,
    right    = right,
    prior    = TRUE
  )

  .hypothesis_log_BF_error_percent(c(posterior_var, prior_var))
}


.hypothesis_region_BF_error_percent <- function(quantity, side) {

  posterior_var <- .hypothesis_region_log_mass_mc_var(
    quantity = quantity,
    side     = side,
    prior    = FALSE
  )
  prior_var     <- .hypothesis_region_log_mass_mc_var(
    quantity = quantity,
    side     = side,
    prior    = TRUE
  )

  .hypothesis_log_BF_error_percent(c(posterior_var, prior_var))
}


.hypothesis_region_log_odds_mc_var <- function(quantity, left, right, prior) {

  if(prior && is.null(quantity[["prior_draws"]])){
    return(0)
  }

  draws <- if(prior){
    quantity[["prior_draws"]]
  }else{
    quantity[["posterior_draws"]]
  }
  left_values  <- .hypothesis_draw_region_indicator(left, draws)
  right_values <- .hypothesis_draw_region_indicator(right, draws)

  .hypothesis_log_odds_indicator_mc_var(left_values, right_values)
}


.hypothesis_region_log_mass_mc_var <- function(quantity, side, prior) {

  if(prior && is.null(quantity[["prior_draws"]])){
    return(0)
  }

  draws <- if(prior){
    quantity[["prior_draws"]]
  }else{
    quantity[["posterior_draws"]]
  }
  values <- .hypothesis_draw_region_indicator(side, draws)

  .hypothesis_log_prob_indicator_mc_var(values)
}


.hypothesis_log_odds_indicator_mc_var <- function(left, right) {

  n <- length(left)
  if(n <= 1L || length(right) != n){
    return(NA_real_)
  }

  left    <- as.numeric(left)
  right   <- as.numeric(right)
  p_left  <- mean(left)
  p_right <- mean(right)
  if(!is.finite(p_left) || !is.finite(p_right) ||
     p_left <= 0 || p_right <= 0){
    return(NA_real_)
  }

  log_var <- stats::var(left) / (n * p_left^2) +
    stats::var(right) / (n * p_right^2) -
    2 * stats::cov(left, right) / (n * p_left * p_right)

  .hypothesis_normalize_log_BF_var(log_var)
}


.hypothesis_log_prob_indicator_mc_var <- function(values) {

  n <- length(values)
  if(n <= 1L){
    return(NA_real_)
  }

  values <- as.numeric(values)
  p      <- mean(values)
  if(!is.finite(p) || p <= 0){
    return(NA_real_)
  }

  log_var <- stats::var(values) / (n * p^2)

  .hypothesis_normalize_log_BF_var(log_var)
}


.hypothesis_normalize_log_BF_var <- function(log_var) {

  if(!is.finite(log_var)){
    return(NA_real_)
  }
  if(log_var < 0 && log_var > -sqrt(.Machine$double.eps)){
    log_var <- 0
  }
  if(log_var < 0){
    return(NA_real_)
  }

  return(log_var)
}


.hypothesis_log_BF_error_percent <- function(log_var) {

  log_var <- as.numeric(log_var)
  if(any(!is.finite(log_var) | log_var < 0)){
    return(NA_real_)
  }

  total <- sum(log_var)
  if(!is.finite(total) || total < 0){
    return(NA_real_)
  }

  return(100 * sqrt(total))
}


.hypothesis_combine_BF_error_percent <- function(...) {

  BF_error <- as.numeric(unlist(list(...), use.names = FALSE))
  if(length(BF_error) == 0L ||
     any(!is.finite(BF_error) | BF_error < 0)){
    return(NA_real_)
  }

  return(100 * sqrt(sum((BF_error / 100)^2)))
}


.hypothesis_point_region_compatible <- function(point_side, region_side) {

  point_key <- .hypothesis_expression_key(point_side[["expr"]])
  if(!is.null(region_side[["expr"]])){
    return(identical(point_key, .hypothesis_expression_key(region_side[["expr"]])))
  }

  region_exprs <- .hypothesis_region_scalar_expressions(
    .hypothesis_parse_expression(region_side[["condition"]])
  )
  if(length(region_exprs) == 0L || anyNA(region_exprs)){
    return(FALSE)
  }

  identical(unique(region_exprs), point_key)
}


.hypothesis_expression_key <- function(text) {

  paste(deparse(.hypothesis_parse_expression(text), width.cutoff = 500L),
        collapse = "")
}


.hypothesis_region_scalar_expressions <- function(expr) {

  if(!is.call(expr)){
    return(NA_character_)
  }

  fun <- as.character(expr[[1L]])
  if(fun %in% c("&", "|")){
    return(unlist(lapply(as.list(expr[-1L]),
                         .hypothesis_region_scalar_expressions),
                  use.names = FALSE))
  }
  if(fun == "!"){
    return(.hypothesis_region_scalar_expressions(expr[[2L]]))
  }
  if(fun %in% c("<", "<=", ">", ">=")){
    lhs <- expr[[2L]]
    rhs <- expr[[3L]]
    lhs_symbols <- .hypothesis_expression_symbols(lhs)
    rhs_symbols <- .hypothesis_expression_symbols(rhs)

    if(length(lhs_symbols) > 0L && length(rhs_symbols) == 0L){
      return(paste(deparse(lhs, width.cutoff = 500L), collapse = ""))
    }
    if(length(lhs_symbols) == 0L && length(rhs_symbols) > 0L){
      return(paste(deparse(rhs, width.cutoff = 500L), collapse = ""))
    }
  }

  return(NA_character_)
}


.hypothesis_direct_symbol <- function(expr_text) {

  expr <- .hypothesis_parse_expression(expr_text)
  if(is.name(expr)){
    return(as.character(expr))
  }

  return(NULL)
}


.hypothesis_region_mass <- function(quantity, side, prior) {

  if(prior && is.null(quantity[["prior_draws"]])){
    if(is.null(quantity[["prior_density"]])){
      stop("Prior information is required for region hypotheses.",
           call. = FALSE)
    }
    mass <- .hypothesis_prior_density_prob(
      prior_density = quantity[["prior_density"]],
      side          = side,
      parameter     = quantity[["parameter"]]
    )
    if(!is.null(side[["complement"]]) && isTRUE(side[["complement"]])){
      mass <- 1 - mass
    }
  }else{
    draws <- if(prior){
      quantity[["prior_draws"]]
    }else{
      quantity[["posterior_draws"]]
    }
    mass <- .hypothesis_draw_region_mass(side, draws)
  }

  return(mass)
}


.hypothesis_draw_region_mass <- function(side, draws) {

  values <- .hypothesis_draw_region_indicator(side, draws)

  return(mean(values))
}


.hypothesis_draw_region_indicator <- function(side, draws) {

  values <- .hypothesis_eval_condition(side[["condition"]], draws)
  if(!is.null(side[["complement"]]) && isTRUE(side[["complement"]])){
    values <- !values
  }

  return(values)
}


.hypothesis_prior_density_prob <- function(prior_density, side, parameter) {

  if(is.null(parameter)){
    stop("A quantity name is required for prior-density region tests.",
         call. = FALSE)
  }
  if(!inherits(prior_density, "prior_linear_density")){
    stop("Prior density is not a deterministic scalar prior density.",
         call. = FALSE)
  }

  prob <- 0
  if(!is.null(prior_density[["density"]])){
    density <- prior_density[["density"]]
    x <- density[["x"]]
    y <- density[["y"]] * density[["mass"]]
    draws <- data.frame(x, check.names = FALSE)
    names(draws) <- parameter
    inside <- .hypothesis_eval_condition(side[["condition"]], draws)
    if(length(x) > 1L){
      prob <- prob + .hypothesis_trapz(x, y * as.numeric(inside))
    }
  }

  points <- prior_density[["points"]]
  if(!is.null(points) && nrow(points) > 0L){
    draws <- data.frame(points[["x"]], check.names = FALSE)
    names(draws) <- parameter
    inside <- .hypothesis_eval_condition(side[["condition"]], draws)
    prob <- prob + sum(points[["p"]][inside])
  }

  prob <- max(0, min(1, prob))

  return(prob)
}


.hypothesis_prior_draws <- function(quantity) {

  if(is.null(quantity[["prior_draws"]])){
    stop("Prior draws are required for this hypothesis expression.",
         call. = FALSE)
  }

  return(quantity[["prior_draws"]])
}


.hypothesis_eval_expression <- function(text, draws) {

  expr <- .hypothesis_parse_expression(text)
  .hypothesis_validate_expression(expr, condition = FALSE)

  draws <- as.data.frame(draws, check.names = FALSE)
  missing <- setdiff(.hypothesis_expression_symbols(expr), names(draws))
  if(length(missing) > 0L){
    stop("Hypothesis expression references unknown quantity '",
         paste(missing, collapse = "', '"), "'.", call. = FALSE)
  }

  env <- list2env(as.list(draws), parent = baseenv())
  values <- eval(expr, envir = env)
  check_real(values, "hypothesis expression", check_length = 0,
             allow_NA = FALSE)
  if(any(!is.finite(values))){
    stop("Hypothesis expression produced non-finite values.", call. = FALSE)
  }

  return(as.numeric(values))
}


.hypothesis_eval_condition <- function(text, draws) {

  expr <- .hypothesis_parse_expression(text)
  .hypothesis_validate_expression(expr, condition = TRUE)

  draws <- as.data.frame(draws, check.names = FALSE)
  missing <- setdiff(.hypothesis_expression_symbols(expr), names(draws))
  if(length(missing) > 0L){
    stop("Hypothesis expression references unknown quantity '",
         paste(missing, collapse = "', '"), "'.", call. = FALSE)
  }

  env <- list2env(as.list(draws), parent = baseenv())
  values <- eval(expr, envir = env)
  if(!is.logical(values)){
    stop("Region hypothesis must evaluate to logical values.", call. = FALSE)
  }
  if(anyNA(values)){
    stop("Region hypothesis produced missing values.", call. = FALSE)
  }

  return(values)
}


.hypothesis_expression_is_parameter <- function(expr_text, parameter) {

  identical(.hypothesis_direct_symbol(expr_text), parameter)
}


.hypothesis_sample_density_height <- function(samples, value, label) {

  sample_sd <- stats::sd(samples)
  if(length(samples) < 2L || !is.finite(sample_sd) || sample_sd <= 0){
    stop("Cannot estimate ", label, " density from degenerate samples.",
         call. = FALSE)
  }

  if(value < min(samples) || value > max(samples)){
    if(label == "prior"){
      stop("Prior samples do not span the point hypothesis.", call. = FALSE)
    }
    return(0)
  }

  density <- stats::density(samples)
  height <- stats::approx(
    density[["x"]],
    density[["y"]],
    xout  = value,
    yleft = 0,
    yright = 0
  )[["y"]]

  if(!is.finite(height) || height < 0){
    stop("Could not estimate ", label, " density at the point hypothesis.",
         call. = FALSE)
  }

  return(height)
}


.hypothesis_prior_density_height <- function(prior_density, value) {

  if(is.null(prior_density)){
    stop("Prior density is required for point hypotheses.", call. = FALSE)
  }

  .prior_linear_density_height(prior_density, value)
}


.hypothesis_check_prior_mass <- function(mass, label) {

  if(!is.finite(mass) || mass <= 0 || mass >= 1){
    stop("Prior region mass for hypothesis '", label,
         "' is zero, one, or non-finite.", call. = FALSE)
  }

  return(invisible(TRUE))
}


.hypothesis_BF_row <- function(quantity, parsed, result) {

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

  format(x, trim = TRUE, scientific = FALSE)
}


.hypothesis_collapse_warning <- function(warning) {

  warning <- unlist(warning, use.names = FALSE)
  warning <- warning[!is.na(warning) & nzchar(warning)]
  if(length(warning) == 0L){
    return(NA_character_)
  }

  paste(unique(warning), collapse = " ")
}
