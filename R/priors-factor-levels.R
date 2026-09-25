#' @title Levels of a factor prior distribution
#'
#' @description Sets the factor levels of a factor prior distribution used
#' outside a formula, e.g., in the \code{prior_list} of [JAGS_fit()]. The
#' prior receives the complete factor metadata that the priors of formula
#' factor terms carry: the number and names of the levels, the contrast, and
#' the design that maps its coefficients to the levels. [JAGS_fit()] requires
#' this metadata of every factor prior outside a formula and stops without it.
#'
#' @param prior a factor prior distribution ([prior_factor()] or
#' [prior_ordered()]), or a mixture or spike-and-slab prior of factor prior
#' distributions.
#' @param levels the number of factor levels, or a character vector with the
#' level names.
#'
#' @details The contrast of the prior determines the number of coefficients:
#' one per level for \code{contrast = "independent"} and ordered priors with
#' \code{contrast = "cumulative_levels"}, and one per level after the first
#' otherwise. Levels given by their number are named \code{"1"}, \code{"2"},
#' and so on. The factor is named after the parameter under which the prior is
#' fitted.
#'
#' @return The prior distribution with its factor levels.
#'
#' @examples
#' p <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
#' p <- prior_factor_levels(p, c("low", "medium", "high"))
#'
#' @seealso [prior_factor()] [prior_ordered()] [JAGS_fit()]
#' @export
prior_factor_levels <- function(prior, levels){

  if(!is.prior(prior) || !.bt_prior_is_factor_family(prior)){
    stop("'prior' must be a factor prior distribution.", call. = FALSE)
  }
  if(is.character(levels)){
    check_char(levels, "levels", check_length = 0, allow_NA = FALSE)
    if(any(!nzchar(levels)) || anyDuplicated(levels) > 0L){
      stop("The 'levels' argument must contain unique, nonempty level names.", call. = FALSE)
    }
    level_names <- levels
  }else{
    check_int(levels, "levels", lower = 1, allow_NA = FALSE)
    level_names <- as.character(seq_len(levels))
  }

  .bt_factor_prior_set_levels(prior, level_names)
}

# Attributes holding the factor metadata of a factor prior.
.bt_factor_metadata_names <- c(
  "levels",
  "coefficient_dim",
  "level_names",
  "interaction",
  "interaction_terms",
  "term_components",
  "factor_terms",
  "factor_contrasts",
  "factor_design",
  "factor_cell_names",
  "ordered_metadata"
)

# The factor metadata that every factor prior carries once its levels are
# known: formula factor terms set it from the data, prior_factor_levels() from
# the given levels, and random-effect SD factor priors from their design.
.bt_factor_metadata_required <- c(
  "levels",
  "level_names",
  "factor_terms",
  "factor_contrasts",
  "factor_design",
  "factor_cell_names"
)

# Factor term of a factor prior whose levels were set outside a formula; the
# factor is named after the parameter when the prior is bound to it.
.bt_factor_placeholder_term <- ".factor"

.bt_factor_prior_set_levels <- function(prior, level_names){

  if(is.prior.mixture(prior) || is.prior.spike_and_slab(prior)){
    components <- which(vapply(prior, is.prior.factor, logical(1)))
    if(length(components) == 0L){
      stop("'prior' must be a factor prior distribution.", call. = FALSE)
    }
    for(component in components){
      prior[[component]] <- .bt_factor_prior_set_levels(prior[[component]], level_names)
    }
    return(.bt_factor_prior_copy_metadata(prior, prior[[components[[1L]]]]))
  }

  contrast <- .factor_object_contrast_name(prior)
  if(is.null(contrast)){
    stop("'prior' must be a factor prior distribution with a contrast.", call. = FALSE)
  }
  for(name in .bt_factor_metadata_names){
    attr(prior, name) <- NULL
  }
  term <- .bt_factor_placeholder_term
  attr(prior, "levels")           <- length(level_names)
  attr(prior, "level_names")      <- level_names
  attr(prior, "factor_terms")     <- term
  attr(prior, "factor_contrasts") <- stats::setNames(contrast, term)
  if(!isTRUE(.get_prior_factor_levels(prior) >= 1)){
    stop(
      "A factor prior with '", contrast, "' contrasts requires at least two levels.",
      call. = FALSE
    )
  }
  design_info <- .factor_term_design_from_metadata(prior)
  attr(prior, "factor_design")     <- design_info[["design"]]
  attr(prior, "factor_cell_names") <- design_info[["cell_names"]]
  if(is.prior.ordered(prior)){
    prior <- .bt_bind_ordered_prior_metadata(prior, ".ordered")
  }

  prior
}

# Sets the factor metadata of a factor prior (and of the factor components of
# a mixture or spike-and-slab factor prior) from a known coefficient design:
# 'level_names' is a character vector for one factor or a list of them named
# by factor, 'factor_contrasts' names the contrast of every factor, and the
# rows of 'design' are the level cells in expand.grid() order.
.bt_factor_prior_set_design <- function(prior, level_names, factor_contrasts,
                                        design){

  level_list <- if(is.list(level_names)) level_names else
    stats::setNames(list(level_names), names(factor_contrasts))
  if(!identical(names(level_list), names(factor_contrasts))){
    stop("The factor contrasts do not match the factor levels.", call. = FALSE)
  }
  cell_names <- .factor_cell_labels(lapply(level_list, as.character))
  if(!is.matrix(design) || nrow(design) != length(cell_names)){
    stop("The factor design does not match the factor levels.", call. = FALSE)
  }
  set_design <- function(x){

    attr(x, "level_names")       <- level_names
    attr(x, "factor_terms")      <- names(factor_contrasts)
    attr(x, "factor_contrasts")  <- factor_contrasts
    attr(x, "factor_design")     <- unname(design)
    attr(x, "factor_cell_names") <- cell_names
    x
  }

  if(is.prior.mixture(prior) || is.prior.spike_and_slab(prior)){
    for(component in which(vapply(prior, is.prior.factor, logical(1)))){
      prior[[component]] <- set_design(prior[[component]])
    }
  }
  set_design(prior)
}

# 'container' with the factor metadata of its factor 'component'.
.bt_factor_prior_copy_metadata <- function(container, component){

  for(name in .bt_factor_metadata_names){
    attr(container, name) <- attr(component, name, exact = TRUE)
  }
  container
}

.bt_factor_metadata_complete <- function(x){

  present <- vapply(.bt_factor_metadata_required, function(name){
    !is.null(attr(x, name, exact = TRUE))
  }, logical(1))
  if(!all(present)){
    return(FALSE)
  }
  design <- as.matrix(attr(x, "factor_design", exact = TRUE))
  nrow(design) == length(attr(x, "factor_cell_names", exact = TRUE)) &&
    isTRUE(ncol(design) == .get_prior_factor_levels(x))
}

.bt_stop_incomplete_factor_metadata <- function(parameter){

  stop(
    "The factor prior ",
    if(!is.null(parameter) && nzchar(parameter)) paste0("of '", parameter, "' "),
    "has no complete factor-level metadata. Set its levels with ",
    "'prior_factor_levels()', or specify the factor in a formula.",
    call. = FALSE
  )
}

# Names the placeholder factor of a prior whose levels were set outside a
# formula after the parameter it is fitted as.
.bt_factor_prior_bind_term <- function(x, parameter){

  if(is.null(parameter) || !nzchar(parameter) ||
     !identical(unname(attr(x, "factor_terms", exact = TRUE)), .bt_factor_placeholder_term)){
    return(x)
  }
  attr(x, "factor_terms") <- parameter
  factor_contrasts <- attr(x, "factor_contrasts", exact = TRUE)
  names(factor_contrasts) <- parameter
  attr(x, "factor_contrasts") <- factor_contrasts
  x
}
