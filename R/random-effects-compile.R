#' Random-effect compilation policy
#'
#' @description
#' `random_effects_compile()` specifies which resolved formula random-effect
#' blocks are compiled as structural marginalized blocks. Blocks not listed in
#' `marginalized` are compiled as ordinary sampled random effects. Marginalized
#' blocks keep their SD, allocation, correlation, covariance-structure, and
#' metadata, but do not generate latent group-level coefficients and do not
#' contribute to the formula linear predictor. Callers that marginalize a block
#' are responsible for using the stored metadata and generated hyperparameter
#' nodes in their likelihood.
#'
#' @param marginalized optional character vector of resolved random-effect
#'   block names that should be compiled structurally without latent
#'   group-level coefficients.
#'
#' @return A list-like S3 object containing the requested `marginalized` block
#'   names, used by [JAGS_formula()], [JAGS_fit()], and [JAGS_bridgesampling()].
#'   Resolved formula-design metadata additionally records the complete
#'   `sampled`, `marginalized`, and named `mode` outputs.
#'
#' @examples
#' random_effects_compile()
#' random_effects_compile(marginalized = "study_estimate")
#'
#' @export
random_effects_compile <- function(marginalized = NULL){

  marginalized <- .bt_validate_random_effects_compile_names(marginalized, "marginalized")

  out <- list(
    marginalized = marginalized
  )
  class(out) <- c("random_effects_compile", "list")

  out
}

.bt_validate_random_effects_compile_names <- function(x, name){

  if(is.null(x)){
    return(NULL)
  }
  check_char(x, name, check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  if(length(x) == 0L){
    return(character())
  }
  if(any(!nzchar(x))){
    stop("'", name, "' cannot contain empty random-effect block names.", call. = FALSE)
  }
  if(anyDuplicated(x)){
    stop("'", name, "' random-effect block names must be unique.", call. = FALSE)
  }

  x
}

is.random_effects_compile <- function(x){
  inherits(x, "random_effects_compile")
}

.bt_check_random_effects_compile <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!is.random_effects_compile(x)){
    stop("'random_effects_compile' must be created with random_effects_compile().", call. = FALSE)
  }

  .bt_validate_random_effects_compile_names(x$marginalized, "marginalized")
  if(inherits(x, "random_effects_compile_resolved")){
    .bt_validate_random_effects_compile_names(x$sampled, "sampled")
    overlap <- intersect(x$sampled, x$marginalized)
    if(length(overlap) > 0L){
      stop(
        "Resolved random-effect block(s) cannot be both sampled and marginalized: ",
        paste(overlap, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_resolve_random_effects_compile <- function(random_effects,
                                               random_effects_compile = NULL){

  .bt_check_random_effects_compile(random_effects_compile, allow_NULL = TRUE)

  block_names <- vapply(random_effects, function(term) term$block_name, character(1))
  requested_marginalized <- if(is.null(random_effects_compile)) NULL else random_effects_compile$marginalized
  unknown <- setdiff(requested_marginalized, block_names)
  if(length(unknown) > 0L){
    stop(
      "random_effects_compile() contains unknown random-effect block(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if(length(block_names) == 0L){
    return(.bt_random_effects_compile_resolved(
      sampled = character(),
      marginalized = character(),
      mode = stats::setNames(character(), character())
    ))
  }
  if(anyDuplicated(block_names)){
    stop(
      "Random-effect block names must be unique after formula parsing.",
      call. = FALSE
    )
  }

  marginalized <- intersect(block_names, requested_marginalized)
  sampled <- setdiff(block_names, marginalized)
  mode <- stats::setNames(rep("sampled", length(block_names)), block_names)
  mode[marginalized] <- "marginalized"

  .bt_random_effects_compile_resolved(
    sampled = sampled,
    marginalized = marginalized,
    mode = mode
  )
}

.bt_random_effects_compile_resolved <- function(sampled, marginalized, mode){

  out <- list(
    sampled = sampled,
    marginalized = marginalized,
    mode = mode
  )
  class(out) <- c("random_effects_compile_resolved", "random_effects_compile", "list")

  out
}

.bt_random_effects_compile_mode <- function(random_effects_compile, block_name){

  if(is.null(random_effects_compile) || is.null(random_effects_compile$mode)){
    return("sampled")
  }
  mode <- random_effects_compile$mode[[block_name]]
  if(is.null(mode) || length(mode) != 1L || is.na(mode) || !nzchar(mode)){
    return("sampled")
  }
  mode
}

.bt_random_effect_term_compile_mode <- function(random_term){

  mode <- random_term$compile_mode
  if(is.null(mode)){
    mode <- attr(random_term, "compile_mode", exact = TRUE)
  }
  if(is.character(mode) && length(mode) == 1L && !is.na(mode) &&
     mode %in% c("sampled", "marginalized")){
    return(mode)
  }

  "sampled"
}

.bt_random_effects_compile_modes_from_terms <- function(random_effects){

  if(length(random_effects) == 0L){
    return(stats::setNames(character(), character()))
  }

  modes <- vapply(random_effects, .bt_random_effect_term_compile_mode, character(1))
  names(modes) <- vapply(random_effects, function(term) term$block_name, character(1))
  modes
}

.bt_formula_design_random_effects <- function(design){

  if(!is.null(design) && inherits(design, "BayesTools_formula_design") &&
     is.list(design$random_effects)){
    return(design$random_effects)
  }

  list()
}

.bt_formula_design_random_effects_by_mode <- function(design, mode){

  random_effects <- .bt_formula_design_random_effects(design)
  if(length(random_effects) == 0L){
    return(list())
  }
  modes <- .bt_random_effects_compile_modes_from_terms(random_effects)
  random_effects[modes == mode]
}

.bt_formula_design_sampled_random_effects <- function(design){

  .bt_formula_design_random_effects_by_mode(design, "sampled")
}

.bt_formula_design_has_any_random_effects <- function(design){

  !is.null(design) &&
    inherits(design, "BayesTools_formula_design") &&
    length(.bt_formula_design_random_effects(design)) > 0L
}

.bt_formula_design_has_sampled_random_effects <- function(design){

  !is.null(design) &&
    inherits(design, "BayesTools_formula_design") &&
    length(.bt_formula_design_sampled_random_effects(design)) > 0L
}

.bt_formula_design_set_random_effects <- function(design, random_effects){

  design$random_effects <- random_effects
  modes <- .bt_random_effects_compile_modes_from_terms(random_effects)
  blocks <- names(modes)
  design$random_effects_compile <- .bt_random_effects_compile_resolved(
    sampled = blocks[modes == "sampled"],
    marginalized = blocks[modes == "marginalized"],
    mode = modes
  )
  design
}
