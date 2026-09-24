# ============================================================================ #
# draws-metadata.R
# ============================================================================ #
#
# One metadata container for posterior and prior draws: mixed and marginal
# posteriors, their elements, and extracted parameter draws carry their
# supports, atoms, mixture components, undefined draws, prior densities and
# their context, precomputed posterior densities and ordinates, formula
# flags, linear weights, and conditioning in the single validated attribute
# 'bayestools_meta'. Every read and write goes through .bt_meta_get() and
# .bt_meta_set() (exact field access, validated values); posterior_metadata()
# exposes the fields downstream packages attach or read. Free attributes of
# earlier development versions are not read.
#
# ============================================================================ #

.bt_meta_attribute <- "bayestools_meta"

# Free attributes that held draw metadata before the container; they must not
# be used on draws (see the unit lint test).
.bt_meta_legacy_names <- c(
  "posterior_support", "posterior_atoms", "posterior_components",
  "models_ind", "sample_ind", "ordered_total_indicator", "undefined_draws",
  "prior_density", "prior_density_context", "prior_densities",
  "posterior_density", "posterior_densities", "posterior_ordinate",
  "posterior_ordinates", "formula_parameter", "formula_log_intercept",
  "formula_scale", "transform_scaled", "linear_weights", "linear_offset",
  "joint_prior_transformation", "conditional", "conditional_rule",
  "condition_key", "condition_event", "resolved_condition_event",
  "effective_conditional", "effective_conditional_rule"
)

# Elements of the 'condition' field (the conditioning of the draws).
.bt_meta_condition_names <- c(
  "conditional", "conditional_rule", "condition_key", "condition_event",
  "resolved_condition_event", "effective_conditional",
  "effective_conditional_rule"
)

# Validators of the draw-metadata fields: each returns NULL for a valid value
# and the reason otherwise.
.bt_meta_validators <- function(){

  list(
    support = function(value){
      # a single support or a list of them keyed by column; malformed support
      # stops with the constructor message
      supports <- if(.posterior_metadata_is_container(value)) value else list(value)
      for(support in supports){
        if(!is.null(support) && is.null(.posterior_support_from_attribute(support))){
          return("it must be created with 'posterior_support_attribute()'")
        }
      }
      NULL
    },
    atoms = function(value){
      if(inherits(value, "BayesTools_posterior_atoms") &&
         !is.null(.posterior_atoms_from_attribute(value))){
        return(NULL)
      }
      "it must be created with 'posterior_atom_attribute()'"
    },
    components = function(value){
      if(inherits(value, "BayesTools_posterior_components")) NULL else
        "it must be a posterior component table"
    },
    models_ind = function(value){
      if(is.numeric(value) && !anyNA(value)) NULL else
        "it must be a numeric vector without missing values"
    },
    sample_ind = function(value){
      if(is.numeric(value) || identical(value, FALSE)) NULL else
        "it must be a numeric vector"
    },
    ordered_total_indicator = function(value){
      if(is.numeric(value) || (is.logical(value) && all(is.na(value)))) NULL else
        "it must be a numeric vector"
    },
    undefined_draws = function(value){
      if(is.character(value) && !anyNA(value)) NULL else
        "it must be a character vector naming the undefined quantities"
    },
    prior_density = function(value){
      if(inherits(value, "prior_density") || is.prior(value)) NULL else
        "it must be a prior density or a prior distribution"
    },
    prior_context = function(value){
      if(.bt_meta_is_prior_context(value)) NULL else
        "it must be a prior-density context"
    },
    prior_densities = function(value){
      if(is.list(value) && (!is.object(value) || inherits(value, "prior_density_list")) &&
         all(vapply(value, function(density){
           is.null(density) || inherits(density, "prior_density") || is.prior(density)
         }, logical(1)))){
        return(NULL)
      }
      "it must be a list of prior densities"
    },
    posterior_density = function(value){
      .posterior_density_kind(value)
      NULL
    },
    posterior_densities = function(value){
      if(!is.list(value) || is.object(value)){
        return("it must be a list of posterior density metadata")
      }
      lapply(value, .posterior_density_kind)
      NULL
    },
    posterior_ordinate = function(value){
      .posterior_ordinate_kind(value)
      NULL
    },
    posterior_ordinates = function(value){
      if(!is.list(value) || is.object(value)){
        return("it must be a list of posterior ordinate metadata")
      }
      lapply(value, .posterior_ordinate_kind)
      NULL
    },
    formula_parameter = function(value){
      if(is.character(value) && length(value) > 0L && !anyNA(value)) NULL else
        "it must name the formula parameter"
    },
    log_intercept = function(value){
      if(is.logical(value) && length(value) == 1L) NULL else
        "it must be TRUE, FALSE, or NA"
    },
    formula_scale = function(value){
      if(is.list(value)) NULL else "it must be a list of formula scaling information"
    },
    transform_scaled = function(value){
      if(isTRUE(value) || isFALSE(value)) NULL else "it must be TRUE or FALSE"
    },
    condition = function(value){
      if(!is.list(value) || is.null(names(value)) ||
         !all(names(value) %in% .bt_meta_condition_names)){
        return("it must be a named list of condition fields")
      }
      for(name in c("conditional", "effective_conditional")){
        if(!is.null(value[[name]]) && !is.character(value[[name]])){
          return(paste0("'", name, "' must be a character vector"))
        }
      }
      for(name in c("conditional_rule", "effective_conditional_rule")){
        if(!is.null(value[[name]]) &&
           !(is.character(value[[name]]) && length(value[[name]]) == 1L &&
             value[[name]] %in% c("AND", "OR"))){
          return(paste0("'", name, "' must be 'AND' or 'OR'"))
        }
      }
      if(!is.null(value[["condition_key"]]) &&
         !(is.character(value[["condition_key"]]) &&
           length(value[["condition_key"]]) == 1L)){
        return("'condition_key' must be a character string")
      }
      for(name in c("condition_event", "resolved_condition_event")){
        if(!is.null(value[[name]]) && !is.list(value[[name]])){
          return(paste0("'", name, "' must be a condition event"))
        }
      }
      NULL
    },
    linear_weights = function(value){
      if(is.numeric(value)) NULL else "it must be a numeric vector or matrix"
    },
    linear_offset = function(value){
      if(is.numeric(value) && length(value) == 1L) NULL else
        "it must be a number"
    },
    joint_prior_transformation = function(value){
      if(is.character(value) && length(value) == 1L) NULL else
        "it must be a character string"
    }
  )
}

.bt_meta_is_prior_context <- function(x){

  inherits(x, "prior_density_context") ||
    inherits(x, "prior_density_model_mixture_context") ||
    inherits(x, "prior_density_conditional_context")
}

.bt_meta_field_names <- c(
  "support", "atoms", "components", "models_ind", "sample_ind",
  "ordered_total_indicator", "undefined_draws", "prior_density",
  "prior_context", "prior_densities", "posterior_density",
  "posterior_densities", "posterior_ordinate", "posterior_ordinates",
  "formula_parameter", "log_intercept", "formula_scale", "transform_scaled",
  "condition", "linear_weights", "linear_offset", "joint_prior_transformation"
)

.bt_meta_fields <- function(){
  .bt_meta_field_names
}

.bt_meta_check_field <- function(field){

  if(!is.character(field) || length(field) != 1L || is.na(field) ||
     !field %in% .bt_meta_fields()){
    stop(
      "'", paste(field, collapse = ", "), "' is not a draw-metadata field.",
      call. = FALSE
    )
  }
  invisible(field)
}

.bt_meta_validate <- function(field, value){

  reason <- .bt_meta_validators()[[field]](value)
  if(!is.null(reason)){
    stop("Draw metadata '", field, "' is invalid: ", reason, ".", call. = FALSE)
  }
  invisible(value)
}

# The draw-metadata container of 'x' (NULL when 'x' carries none).
.bt_meta_container <- function(x){

  meta <- attr(x, .bt_meta_attribute, exact = TRUE)
  if(is.null(meta)){
    return(NULL)
  }
  if(!inherits(meta, "BayesTools_draw_metadata") || !is.list(meta) ||
     is.null(names(meta)) || !all(names(meta) %in% .bt_meta_fields())){
    stop("The draw-metadata container is invalid.", call. = FALSE)
  }
  meta
}

# The value of a draw-metadata field, or NULL when it is absent.
.bt_meta_get <- function(x, field){

  .bt_meta_check_field(field)
  meta <- .bt_meta_container(x)
  if(is.null(meta)){
    return(NULL)
  }
  meta[[field]]
}

# 'x' with the draw-metadata field set to 'value'; NULL removes the field.
.bt_meta_set <- function(x, field, value){

  .bt_meta_check_field(field)
  meta <- .bt_meta_container(x)
  if(is.null(value) && is.null(meta)){
    return(x)
  }
  if(!is.null(value)){
    .bt_meta_validate(field, value)
  }
  if(is.null(meta)){
    meta <- structure(list(), names = character(), class = c("BayesTools_draw_metadata", "list"))
  }
  meta[[field]] <- value
  attr(x, .bt_meta_attribute) <- if(length(meta) > 0L) meta
  x
}

# 'x' with several draw-metadata fields set, given as named arguments.
.bt_meta_update <- function(x, ...){

  fields <- list(...)
  if(length(fields) > 0L &&
     (is.null(names(fields)) || any(!nzchar(names(fields))))){
    stop("Draw-metadata fields must be named.", call. = FALSE)
  }
  for(field in names(fields)){
    x <- .bt_meta_set(x, field, fields[[field]])
  }
  x
}

# One element of the 'condition' field.
.bt_meta_condition <- function(x, name){

  if(!name %in% .bt_meta_condition_names){
    stop("'", name, "' is not a draw-conditioning field.", call. = FALSE)
  }
  .bt_meta_get(x, "condition")[[name]]
}

# Fields of the public accessor posterior_metadata().
.bt_meta_public_fields <- c(
  "support", "atoms", "undefined_draws", "prior_density", "prior_context",
  "posterior_density", "posterior_densities", "posterior_ordinate",
  "posterior_ordinates", "condition"
)

#' @title Metadata of BayesTools posterior draws
#'
#' @description Reads and sets the metadata that BayesTools attaches to
#' posterior draws: the elements and lists returned by [mix_posteriors()],
#' [as_mixed_posteriors()], [marginal_posterior()], and
#' [random_effects_summary_posterior()], and the draws of [parameter_draws()].
#' All metadata are stored in one attribute and validated when they are set;
#' use this accessor instead of setting attributes.
#'
#' @param x posterior draws (an element or a list of them).
#' @param field metadata field:
#' \describe{
#'   \item{\code{"support"}}{exact support, created with
#'   [posterior_support_attribute()] (or a list of such objects keyed by
#'   column).}
#'   \item{\code{"atoms"}}{declared posterior point masses, created with
#'   [posterior_atom_attribute()]; point masses are read only from this
#'   field.}
#'   \item{\code{"undefined_draws"}}{a character vector naming quantities
#'   whose draws may be \code{NA} because the quantity is undefined in those
#'   draws (e.g. \code{"correlation"}), as returned by [parameter_draws()].}
#'   \item{\code{"prior_density"}}{the prior density of the quantity.}
#'   \item{\code{"prior_context"}}{the joint prior-density context of the
#'   draws.}
#'   \item{\code{"posterior_density"}, \code{"posterior_densities"}}{
#'   precomputed posterior densities created with
#'   [posterior_density_attribute()].}
#'   \item{\code{"posterior_ordinate"}, \code{"posterior_ordinates"}}{
#'   precomputed posterior ordinates created with
#'   [posterior_ordinate_attribute()] or [posterior_ordinate_append()].}
#'   \item{\code{"condition"}}{the conditioning of the draws: a list with
#'   \code{conditional}, \code{conditional_rule}, \code{condition_key},
#'   \code{condition_event}, \code{resolved_condition_event}, and, for
#'   levels whose conditioning was resolved per level,
#'   \code{effective_conditional} and \code{effective_conditional_rule}.}
#' }
#' @param value the new value of the field; \code{NULL} removes it.
#'
#' @return \code{posterior_metadata()} returns the value of the field or
#' \code{NULL}; the replacement form returns \code{x} with the field set.
#'
#' @examples
#' draws <- stats::rnorm(100)
#' posterior_metadata(draws, "atoms") <- posterior_atom_attribute()
#' posterior_metadata(draws, "atoms")
#'
#' @export
posterior_metadata <- function(x, field){

  .bt_meta_check_public_field(field)
  .bt_meta_get(x, field)
}

#' @rdname posterior_metadata
#' @export
`posterior_metadata<-` <- function(x, field, value){

  .bt_meta_check_public_field(field)
  .bt_meta_set(x, field, value)
}

.bt_meta_check_public_field <- function(field){

  check_char(field, "field", allow_values = .bt_meta_public_fields)
  invisible(field)
}
