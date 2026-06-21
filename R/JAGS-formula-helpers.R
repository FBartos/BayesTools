# formula helper functions
.remove_response        <- function(formula){
  # removes response from the expression
  # (prevents crash on formula evaluations)
  if(attr(stats::terms(formula), "response")  == 1){
    formula[2] <- NULL
  }
  return(formula)
}
.has_expression         <- function(formula){
  # check if there is any expression in the formula
  return(.bt_contains_expression_call(.bt_formula_rhs(formula)))
}
.extract_expressions    <- function(formula){
  # extract all expressions from the formula

  return(.bt_extract_expression_bodies(.bt_formula_rhs(formula)))
}
.clean_from_expression  <- function(x){
  # expression to character

  return(sub("expression\\((.*)\\)", "\\1", x))
}
.remove_expressions     <- function(formula){
  # remove all expressions from the formula

  rhs_index <- .bt_formula_rhs_index(formula)
  rhs <- .bt_remove_expression_terms(formula[[rhs_index]])
  if(is.null(rhs)){
    rhs <- 1
  }
  formula[[rhs_index]] <- rhs

  return(formula)
}
.bt_formula_rhs_index <- function(formula){
  if(length(formula) == 3L) 3L else 2L
}
.bt_formula_rhs <- function(formula){
  formula[[.bt_formula_rhs_index(formula)]]
}
.bt_is_expression_call <- function(x){
  is.call(x) && identical(as.character(x[[1L]]), "expression")
}
.bt_contains_expression_call <- function(x){
  if(.bt_is_expression_call(x)){
    return(TRUE)
  }
  if(is.call(x) || is.pairlist(x)){
    return(any(vapply(as.list(x)[-1], .bt_contains_expression_call, logical(1))))
  }

  FALSE
}
.bt_extract_expression_bodies <- function(x){
  if(.bt_is_expression_call(x)){
    return(lapply(as.list(x)[-1], function(expression_body){
      paste0(deparse(expression_body), collapse = " ")
    }))
  }
  if(is.call(x) || is.pairlist(x)){
    return(unlist(lapply(as.list(x)[-1], .bt_extract_expression_bodies), recursive = FALSE))
  }

  list()
}
.bt_remove_expression_terms <- function(x){
  if(.bt_is_expression_call(x)){
    return(NULL)
  }

  if(is.call(x)){
    call_name <- as.character(x[[1L]])
    if(call_name == "+" && length(x) == 3L){
      lhs <- .bt_remove_expression_terms(x[[2L]])
      rhs <- .bt_remove_expression_terms(x[[3L]])
      if(is.null(lhs)){
        return(rhs)
      }
      if(is.null(rhs)){
        return(lhs)
      }
      return(call("+", lhs, rhs))
    }
    if(call_name == "-" && length(x) == 3L){
      lhs <- .bt_remove_expression_terms(x[[2L]])
      rhs <- .bt_remove_expression_terms(x[[3L]])
      if(is.null(lhs) && is.null(rhs)){
        return(NULL)
      }
      if(is.null(rhs)){
        return(lhs)
      }
      if(is.null(lhs)){
        stop("expression() terms must be additive formula terms.", call. = FALSE)
      }
      return(call("-", lhs, rhs))
    }
    if(.bt_contains_expression_call(x)){
      stop("expression() terms must be additive formula terms.", call. = FALSE)
    }
  }

  x
}
.has_random_effects     <- function(formula){
  return(length(.bt_parse_random_effects(formula)$terms) > 0L)
}
.remove_random_effects  <- function(formula){
  return(.bt_fixed_formula(formula))
}
.get_grouping_factor    <- function(x){
  has_grouping            <- grepl("\\|", x)
  grouping                <- rep("", length(x))
  grouping[has_grouping]  <- trimws(sub(".*\\|\\s*", "", x[has_grouping]))
  return(grouping)
}
.remove_grouping_factor <- function(formula){
  return(trimws(sub("\\|.*$", "", formula)))
}
#' @title Add an Intercept to a Formula
#'
#' @description Converts a no-intercept formula to the corresponding formula
#' with an intercept while preserving the formula environment. Top-level
#' no-intercept encodings such as \code{- 1}, \code{+ 0}, and \code{0 +} are
#' removed without editing transformed calls such as \code{I(x - 1)} or
#' \code{offset(x - 1)}.
#'
#' @param formula a formula object.
#'
#' @return A formula object with an intercept.
#'
#' @export
formula_add_intercept <- function(formula){

  if(!inherits(formula, "formula")){
    stop("'formula' must be a formula.", call. = FALSE)
  }

  if(attr(stats::terms(formula), "intercept") == 1L){
    return(formula)
  }

  formula_env   <- environment(formula)
  formula_attrs <- attributes(formula)
  rhs_index     <- if(length(formula) == 3L) 3L else 2L
  rhs           <- .formula_strip_no_intercept(formula[[rhs_index]])

  if(is.null(rhs)){
    rhs <- 1
  }

  out <- formula
  out[[rhs_index]] <- rhs
  environment(out) <- formula_env

  for(attribute in setdiff(names(formula_attrs), c("class", ".Environment", "names"))){
    attr(out, attribute) <- formula_attrs[[attribute]]
  }

  if(attr(stats::terms(out), "intercept") == 0L){
    out[[rhs_index]] <- call("+", 1, out[[rhs_index]])
    environment(out) <- formula_env
  }

  return(out)
}
.add_intercept_to_formula <- formula_add_intercept

.formula_strip_no_intercept <- function(expr){

  if(.formula_is_no_intercept_additive_term(expr)){
    return(NULL)
  }

  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    lhs <- .formula_strip_no_intercept(expr[[2L]])
    rhs <- .formula_strip_no_intercept(expr[[3L]])

    if(is.null(lhs)){
      return(rhs)
    }
    if(is.null(rhs)){
      return(lhs)
    }
    return(call("+", lhs, rhs))
  }

  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 3L &&
     .formula_is_numeric_constant(expr[[3L]], 1)){
    return(.formula_strip_no_intercept(expr[[2L]]))
  }

  return(expr)
}

.formula_is_no_intercept_additive_term <- function(expr){

  .formula_is_numeric_constant(expr, 0) || .formula_is_negative_one(expr)
}

.formula_is_negative_one <- function(expr){

  (is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), -1)) ||
    (is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L &&
       .formula_is_numeric_constant(expr[[2L]], 1))
}

.formula_is_numeric_constant <- function(expr, value){

  if(is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L){
    return(.formula_is_numeric_constant(expr[[2L]], value))
  }

  is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), as.numeric(value))
}

