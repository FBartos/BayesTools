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
.bt_validate_fixed_formula_grammar <- function(formula){

  validate_expression <- function(expression){
    if(is.symbol(expression)){
      if(identical(as.character(expression), ".")){
        stop(
          "Unsupported fixed-formula term '.'. Dot expansion is not supported; ",
          "list each data-frame column explicitly.",
          call. = FALSE
        )
      }
      return(invisible(TRUE))
    }
    if(is.numeric(expression) && length(expression) == 1L &&
       is.finite(expression) && expression %in% c(0, 1)){
      return(invisible(TRUE))
    }
    if(is.call(expression)){
      call_name <- if(is.symbol(expression[[1L]])){
        as.character(expression[[1L]])
      }else{
        ""
      }
      if(call_name %in% c("+", "-", "*", ":", "/", "^", "(")){
        for(argument in as.list(expression)[-1L]){
          validate_expression(argument)
        }
        return(invisible(TRUE))
      }

      expression_label <- .bt_deparse_expr(expression)
      if(identical(call_name, "offset")){
        stop(
          "Unsupported fixed-formula call '", expression_label,
          "'. offset() is not supported; use expression(...) for an explicit ",
          "JAGS-scale offset.",
          call. = FALSE
        )
      }
      stop(
        "Unsupported fixed-formula call '", expression_label,
        "'. Create the transformed value as a data-frame column and reference ",
        "that column by name.",
        call. = FALSE
      )
    }

    stop(
      "Unsupported fixed-formula expression '",
      .bt_deparse_expr(expression),
      "'. Use literal data-frame column names.",
      call. = FALSE
    )
  }

  validate_expression(.bt_formula_rhs(formula))
  invisible(TRUE)
}
.bt_validate_formula_replay_grammar <- function(formula){

  if(!inherits(formula, "formula")){
    stop("'formula' must be a formula.", call. = FALSE)
  }

  validation_formula <- .remove_response(formula)
  validation_formula <- .remove_expressions(validation_formula)
  validation_formula <- .remove_random_effects(validation_formula)
  .bt_validate_fixed_formula_grammar(validation_formula)
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
    call_name <- if(is.symbol(x[[1L]])) as.character(x[[1L]]) else ""
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
      if(.bt_contains_expression_call(x[[3L]])){
        stop("expression() terms must be additive formula terms.", call. = FALSE)
      }
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
.bt_formula_random_terms <- function(formula){

  if(inherits(formula, "BayesTools_random_effects")){
    return(formula$terms)
  }
  random_terms <- attr(formula, "random_terms", exact = TRUE)
  if(is.list(random_terms)){
    return(random_terms)
  }

  .bt_parse_random_effects(formula)$terms
}
.bt_formula_random_formula <- function(formula){

  if(inherits(formula, "BayesTools_random_effects")){
    return(formula$formula)
  }

  formula
}
.bt_formula_preserve_random_terms <- function(formula, random_terms){

  if(is.list(random_terms)){
    attr(formula, "random_terms") <- random_terms
  }
  formula
}
.has_random_effects     <- function(formula){
  return(length(.bt_formula_random_terms(formula)) > 0L)
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
.JAGS_formula_default_prior_names <- function(){
  c("__default_continuous", "__default_factor")
}
.JAGS_formula_is_lazy_default_prior <- function(x){
  is.function(x) && typeof(x) == "closure" && length(formals(x)) == 0L
}
.JAGS_formula_check_prior_list <- function(prior_list){

  default_prior_names <- .JAGS_formula_default_prior_names()
  prior_names         <- names(prior_list)

  for(i in seq_along(prior_list)){
    prior_name  <- if(is.null(prior_names)) "" else prior_names[[i]]
    prior_value <- prior_list[[i]]

    if(prior_name %in% default_prior_names){
      if(is.prior(prior_value) || .JAGS_formula_is_lazy_default_prior(prior_value)){
        next
      }
      if(is.function(prior_value)){
        stop(
          paste0(
            "The 'prior_list[[\"", prior_name, "\"]]' entry must be a prior object ",
            "or a zero-argument function returning a prior object."
          ),
          call. = FALSE
        )
      }
    }

    if(!is.prior(prior_value)){
      stop("'prior_list' must be a list of priors.", call. = FALSE)
    }
  }

  return()
}
.JAGS_formula_resolve_default_prior <- function(default_prior, default_name){

  if(.JAGS_formula_is_lazy_default_prior(default_prior)){
    default_prior <- default_prior()
    if(!is.prior(default_prior)){
      stop(
        paste0(
          "The 'prior_list[[\"", default_name, "\"]]' lazy default must return ",
          "a BayesTools prior object."
        ),
        call. = FALSE
      )
    }
  }

  return(default_prior)
}
.JAGS_formula_canonicalize_none_prior <- function(prior_object){

  if(!is.prior.none(prior_object)){
    return(prior_object)
  }

  output <- prior(
    "point",
    list(location = 0),
    prior_weights = .prior_model_weight(prior_object)
  )
  prior_attributes <- attributes(prior_object)
  metadata_names <- setdiff(names(prior_attributes), c("names", "class"))
  for(metadata_name in metadata_names){
    attr(output, metadata_name) <- prior_attributes[[metadata_name]]
  }

  output
}
.bt_formula_prior_is_factor <- function(x){

  is.prior.factor(x) ||
    inherits(x, "prior.factor_mixture") ||
    inherits(x, "prior.factor_spike_and_slab")
}
.bt_validate_formula_term_priors <- function(prior_list, model_terms,
                                             model_terms_type){

  for(model_term in model_terms){
    this_prior <- prior_list[[model_term]]
    term_type <- model_terms_type[[model_term]]
    factor_prior <- .bt_formula_prior_is_factor(this_prior)

    if(factor_prior){
      .validate_centered_factor_prior(
        this_prior,
        paste0("prior_list[[\"", model_term, "\"]]")
      )
    }
    if(identical(term_type, "factor") && !factor_prior){
      stop(
        "Unsupported prior distribution defined for '", model_term,
        "' factor variable. See '?prior_factor' for details.",
        call. = FALSE
      )
    }
    if(identical(term_type, "continuous") &&
       (factor_prior || is.prior.discrete(this_prior) ||
        is.prior.PET(this_prior) || is.prior.PEESE(this_prior) ||
        is.prior.weightfunction(this_prior))){
      stop(
        "Unsupported prior distribution defined for '", model_term,
        "' continuous variable. See '?prior' for details.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
.bt_validate_formula_reconstruction_prior <- function(prior_object,
                                                       prior_name){

  if(is.prior.point(prior_object) ||
     is.prior.factor(prior_object) ||
     is.prior.simple(prior_object)){
    return(invisible(TRUE))
  }

  stop(
    "Unsupported formula reconstruction prior for '", prior_name,
    "'. Formula metadata must contain a canonical simple or factor prior.",
    call. = FALSE
  )
}
.bt_validate_formula_log_intercept_prior <- function(prior_list,
                                                      parameter = NULL){

  intercept_name <- if(is.null(parameter)){
    "intercept"
  }else{
    paste0(parameter, "_intercept")
  }
  if(!intercept_name %in% names(prior_list)){
    stop(
      "A formula using log(intercept) must define a prior for '",
      intercept_name, "'.",
      call. = FALSE
    )
  }

  .validate_strictly_positive_prior(
    prior_list[[intercept_name]],
    if(is.null(parameter)){
      "prior_list[[\"intercept\"]]"
    }else{
      paste0("formula_prior_list[[\"", intercept_name, "\"]]")
    }
  )
}
.remove_grouping_factor <- function(formula){
  return(trimws(sub("\\|.*$", "", formula)))
}
#' @title Add an Intercept to a Formula
#'
#' @description Converts a no-intercept formula to the corresponding formula
#' with an intercept while preserving the formula environment. Additive,
#' parenthesized, and unary-plus no-intercept encodings such as \code{- 1},
#' \code{+ 0}, and \code{0 +} are removed without editing transformed calls
#' such as \code{I(x - 1)} or \code{offset(x - 1)}.
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

  if(is.call(expr) && length(expr) == 2L &&
     (identical(expr[[1L]], as.name("(")) ||
      identical(expr[[1L]], as.name("+")))){
    return(.formula_strip_no_intercept(expr[[2L]]))
  }

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

  if(is.call(expr) && length(expr) == 2L &&
     (identical(expr[[1L]], as.name("(")) ||
      identical(expr[[1L]], as.name("+")))){
    return(.formula_is_numeric_constant(expr[[2L]], value))
  }

  is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), as.numeric(value))
}
