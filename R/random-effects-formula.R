# Helpers for random-effect formula/list normalization.

#' Normalize random-effect formula input
#'
#' @description
#' Normalizes a formula or a list of formulas containing random-effect terms
#' into a single formula and parsed BayesTools random-effect terms.
#'
#' @details
#' List names define top-level random-effect components used by allocation and
#' summary metadata. Missing names are replaced by `"Component 1"`,
#' `"Component 2"`, and so on, then sanitized for JAGS node names. Unnamed
#' random-effect terms inside a list entry are auto-named from the component
#' label. Explicit `name =` arguments define the final random-effect block name
#' and are the names targeted by `prior_random()` block overrides.
#'
#' @param random a formula or non-empty list of formulas.
#' @param envir environment used for the returned formula.
#'
#' @return A `BayesTools_random_effects` object with `formula`, `terms`, and
#'   `components` fields.
#'
#' @export
random_effects_formula <- function(random, envir = parent.frame()){

  formulas <- .bt_random_effects_input_formulas(random)
  is_list_input <- is.list(random) && !inherits(random, "formula")
  component_names <- if(is_list_input){
    .bt_random_effects_input_names(random, length(formulas))
  }else{
    rep(NA_character_, length(formulas))
  }
  component_labels <- rep(NA_character_, length(formulas))
  if(is_list_input){
    component_labels <- vapply(
      component_names,
      .bt_random_effect_sanitize_name,
      character(1)
    )
    if(anyDuplicated(component_labels)){
      stop(
        "Random-effect list names must be unique after sanitization.",
        call. = FALSE
      )
    }
  }

  normalized_terms <- list()
  normalized_calls <- list()
  components <- list()

  for(i in seq_along(formulas)){
    formula <- formulas[[i]]
    if(length(formula) == 3L){
      warning(
        "The left-hand side of random-effect formulas is ignored.",
        call. = FALSE
      )
      formula <- formula[-2L]
    }

    parsed <- .bt_parse_random_effects(formula)
    terms <- parsed$terms
    if(length(terms) == 0L){
      stop("Random-effect formulas must contain at least one random-effect term.", call. = FALSE)
    }
    .bt_random_effects_validate_random_only(parsed)

    component <- component_names[[i]]
    component_label <- component_labels[[i]]
    if(!is.na(component_label)){
      terms <- .bt_random_effects_apply_component_name(
        terms = terms,
        component = component,
        component_label = component_label
      )
      components[[component_label]] <- vapply(terms, `[[`, character(1), "block_name")
    }

    for(term in terms){
      normalized_terms[[length(normalized_terms) + 1L]] <- term
      normalized_calls[[length(normalized_calls) + 1L]] <- .bt_random_effect_named_special_call(term)
    }
  }
  .bt_validate_random_effect_block_names(normalized_terms)

  formula <- stats::as.formula(
    call("~", .bt_random_effect_plus_calls(normalized_calls)),
    env = envir
  )
  attr(formula, "random_components") <- components
  attr(formula, "random_terms") <- normalized_terms

  out <- list(
    formula = formula,
    terms = normalized_terms,
    components = components
  )
  class(out) <- c("BayesTools_random_effects", "list")

  out
}

.bt_random_effects_input_formulas <- function(random){

  if(inherits(random, "formula")){
    return(list(random))
  }
  if(is.list(random) && length(random) > 0L &&
     all(vapply(random, inherits, logical(1), what = "formula"))){
    return(random)
  }

  stop(
    "'random' must be a formula or a non-empty list of formulas.",
    call. = FALSE
  )
}

.bt_random_effects_input_names <- function(random, n){

  random_names <- names(random)
  if(is.null(random_names) || length(random_names) != n){
    random_names <- rep("", n)
  }
  random_names[is.na(random_names)] <- ""
  missing_names <- !nzchar(random_names)
  random_names[missing_names] <- paste0("Component ", which(missing_names))

  random_names
}

.bt_random_effects_validate_random_only <- function(parsed){

  fixed_terms <- attr(stats::terms(parsed$fixed_formula), "term.labels")
  if(length(fixed_terms) > 0L){
    stop(
      "Random-effect formulas must contain only random-effect terms.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_effects_apply_component_name <- function(terms, component,
                                                   component_label){

  explicit <- vapply(terms, function(term){
    isTRUE(term$has_explicit_name)
  }, logical(1))

  original_blocks <- vapply(terms, `[[`, character(1), "block_name")
  for(i in seq_along(terms)){
    child_label <- original_blocks[[i]]
    if(!explicit[[i]]){
      terms[[i]]$block_name <- if(length(terms) == 1L){
        component_label
      }else{
        paste0(component_label, "_", child_label)
      }
    }
    terms[[i]]$component <- component
    terms[[i]]$component_label <- component_label
    terms[[i]]$component_child_label <- if(length(terms) == 1L){
      component_label
    }else{
      child_label
    }
    attr(terms[[i]], "random_block") <- terms[[i]]$block_name
  }

  terms
}

.bt_random_effect_named_special_call <- function(term){

  call_name <- switch(
    term$structure,
    id = "id",
    diag = "diag",
    cs = "cs",
    hcs = "hcs",
    ar1 = "ar1",
    car = "car",
    har = "har",
    us = "us",
    term$structure
  )

  call_args <- list(as.name(call_name), term$bar_call)
  call_args$name <- term$block_name
  if(!is.null(term$hom)){
    call_args$hom <- term$hom
  }

  as.call(call_args)
}

.bt_random_effect_plus_calls <- function(calls){

  if(length(calls) == 0L){
    stop("At least one random-effect call is required.", call. = FALSE)
  }
  if(length(calls) == 1L){
    return(calls[[1L]])
  }

  out <- calls[[1L]]
  for(i in seq.int(2L, length(calls))){
    out <- call("+", out, calls[[i]])
  }

  out
}
