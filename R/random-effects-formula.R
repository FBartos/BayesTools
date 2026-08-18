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
#' Random-effect terms have two distinct left-side grammars:
#'
#' * `id()`, `diag()`, and `us()` / `un()` use a coefficient formula. The
#'   `1`, `0`, and `-1` terms control the random intercept; continuous slopes,
#'   factor slopes, and interactions are supported. Plain `(expr | group)`
#'   defaults to `us()`, and `expr || group` to `diag()`.
#' * `cs()` / `hcs()`, `ar1()` / `ar()` / `har()`, and `car()` use a
#'   structure-owned index specification. They reject explicit intercept
#'   controls and do not accept `random_block(contrasts = ...)`.
#'
#' Random-coefficient factor coding comes from the block's resolved contrast
#' metadata. It reuses an already established fixed-factor contrast by default;
#' `random_block(contrasts = ...)` is the explicit per-block override. Removing
#' the random intercept does not by itself request one coefficient per factor
#' level.
#'
#' The discrete index structures `cs()`, `hcs()`, `ar1()`, and `har()` accept
#' factor, character, numeric/integer, or logical columns. Existing factor
#' levels are preserved; otherwise sorted unique values define levels and
#' AR(1) order. CS/HCS can combine several columns with `+`; AR1/HAR require
#' exactly one. CAR requires one finite numeric/integer column, or an ordered
#' factor with numeric labels, and uses actual coordinate distances.
#'
#' Random-effect predictors must be literal data-column names combined with
#' standard formula operators. Inline transformations and arbitrary calls are
#' rejected; create transformed predictors as explicit columns first. Grouping
#' terms support variables, `:` interactions, and `/` nesting. Grouping
#' interactions follow base-R/lme4 lexicographic tuple order, and the fitted
#' tuple-to-index map is stored for replay.
#' Categorical predictor and grouping levels must not contain BayesTools'
#' reserved internal tokens, such as `__xXx__`.
#'
#' @param random a formula or non-empty list of formulas.
#' @param envir environment used for the returned formula.
#' @param group_covariance optional known group covariance. Use one
#'   [random_group_covariance()] object or numeric matrix for a single random
#'   block, or a named list of objects/matrices for multiple random blocks.
#'   Names should match resolved random-effect block names, or grouping labels
#'   when the grouping label is unambiguous.
#'
#' @return A `BayesTools_random_effects` object with `formula`, `terms`, and
#'   `components` fields.
#'
#' @export
random_effects_formula <- function(random, envir = parent.frame(),
                                   group_covariance = NULL){

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
  normalized_terms <- .bt_random_effects_apply_group_covariance(
    terms = normalized_terms,
    group_covariance = group_covariance
  )

  formula <- stats::as.formula(
    call("~", .bt_random_effect_plus_calls(normalized_calls)),
    env = envir
  )
  attr(formula, "random_components") <- components
  attr(formula, "random_terms") <- normalized_terms

  out <- list(
    formula = formula,
    terms = normalized_terms,
    components = components,
    group_covariance = .bt_random_effects_group_covariance_metadata(
      normalized_terms
    )
  )
  class(out) <- c("BayesTools_random_effects", "list")

  out
}

.bt_random_effects_apply_group_covariance <- function(terms,
                                                      group_covariance = NULL){

  if(is.null(group_covariance)){
    return(terms)
  }
  entries <- .bt_random_effects_group_covariance_entries(group_covariance)
  entry_names <- names(entries)
  if(is.null(entry_names)){
    entry_names <- rep("", length(entries))
  }
  entry_names[is.na(entry_names)] <- ""

  if(length(entries) == 1L && !nzchar(entry_names[[1L]])){
    if(length(terms) != 1L){
      stop(
        "Unnamed 'group_covariance' is only supported for a single random-effect block.",
        call. = FALSE
      )
    }
    terms[[1L]]$group_covariance <- entries[[1L]]
    return(terms)
  }
  if(any(!nzchar(entry_names))){
    stop(
      "'group_covariance' must be named when multiple entries are supplied.",
      call. = FALSE
    )
  }
  if(anyDuplicated(entry_names)){
    stop("'group_covariance' names must be unique.", call. = FALSE)
  }

  matched <- rep(FALSE, length(terms))
  for(i in seq_along(entries)){
    index <- .bt_random_effects_group_covariance_match(
      terms = terms,
      name = entry_names[[i]]
    )
    if(matched[[index]]){
      stop(
        "Random-effect block '",
        terms[[index]]$block_name,
        "' has more than one known group covariance entry.",
        call. = FALSE
      )
    }
    terms[[index]]$group_covariance <- entries[[i]]
    matched[[index]] <- TRUE
  }

  terms
}

.bt_random_effects_group_covariance_entries <- function(group_covariance){

  if(.bt_is_random_group_covariance(group_covariance) ||
     is.matrix(group_covariance) ||
     is.data.frame(group_covariance)){
    return(list(.bt_as_random_group_covariance(group_covariance)))
  }
  if(!is.list(group_covariance) || length(group_covariance) == 0L){
    stop(
      "'group_covariance' must be a random_group_covariance() object, a numeric matrix, or a named list.",
      call. = FALSE
    )
  }

  lapply(group_covariance, .bt_as_random_group_covariance)
}

.bt_random_effects_group_covariance_match <- function(terms, name){

  block_names <- vapply(terms, `[[`, character(1), "block_name")
  block_matches <- which(block_names == name)
  if(length(block_matches) == 1L){
    return(block_matches)
  }
  if(length(block_matches) > 1L){
    stop(
      "'group_covariance' name '",
      name,
      "' matches multiple random-effect blocks.",
      call. = FALSE
    )
  }

  group_labels <- vapply(terms, `[[`, character(1), "group_label")
  group_matches <- which(group_labels == name)
  if(length(group_matches) == 0L){
    sanitized <- vapply(
      group_labels,
      .bt_random_effect_sanitize_name,
      character(1)
    )
    group_matches <- which(sanitized == name)
  }
  if(length(group_matches) == 1L){
    return(group_matches)
  }
  if(length(group_matches) > 1L){
    stop(
      "'group_covariance' name '",
      name,
      "' matches an ambiguous random-effect grouping label.",
      call. = FALSE
    )
  }

  stop(
    "'group_covariance' name '",
    name,
    "' does not match any random-effect block or unambiguous grouping label.",
    call. = FALSE
  )
}

.bt_random_effects_group_covariance_metadata <- function(terms){

  out <- lapply(terms, function(term){
    if(is.null(term$group_covariance)){
      return(NULL)
    }
    list(
      block_name = term$block_name,
      group_label = term$group_label,
      scale = term$group_covariance$scale,
      levels = rownames(term$group_covariance$covariance)
    )
  })
  names(out) <- vapply(terms, `[[`, character(1), "block_name")
  out[!vapply(out, is.null, logical(1))]
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
