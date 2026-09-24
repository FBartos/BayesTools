# ============================================================================ #
# draws-metadata.R
# ============================================================================ #
#
# Accessors for the metadata of posterior and prior draws (mixed and marginal
# posteriors, their elements, and extracted parameter draws): supports, atoms,
# mixture components, undefined draws, prior densities and their context,
# precomputed posterior densities and ordinates, formula flags, linear
# weights, and conditioning. Every read and write of these fields goes
# through .bt_meta_get()/.bt_meta_set(), which access them exactly (never by
# partial attribute-name matching) and validate the field name.
#
# ============================================================================ #

# Draw-metadata fields and the attribute that stores each of them.
.bt_meta_attribute_names <- c(
  support                    = "posterior_support",
  atoms                      = "posterior_atoms",
  components                 = "posterior_components",
  models_ind                 = "models_ind",
  sample_ind                 = "sample_ind",
  ordered_total_indicator    = "ordered_total_indicator",
  undefined_draws            = "undefined_draws",
  prior_density              = "prior_density",
  prior_context              = "prior_density_context",
  prior_densities            = "prior_densities",
  posterior_density          = "posterior_density",
  posterior_densities        = "posterior_densities",
  posterior_ordinate         = "posterior_ordinate",
  posterior_ordinates        = "posterior_ordinates",
  formula_parameter          = "formula_parameter",
  log_intercept              = "formula_log_intercept",
  formula_scale              = "formula_scale",
  transform_scaled           = "transform_scaled",
  linear_weights             = "linear_weights",
  linear_offset              = "linear_offset",
  joint_prior_transformation = "joint_prior_transformation"
)

# Elements of the 'condition' field (the conditioning of the draws).
.bt_meta_condition_names <- c(
  "conditional", "conditional_rule", "condition_key", "condition_event",
  "resolved_condition_event", "effective_conditional",
  "effective_conditional_rule"
)

.bt_meta_fields <- function(){
  c(names(.bt_meta_attribute_names), "condition")
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

# The value of a draw-metadata field, or NULL when it is absent.
.bt_meta_get <- function(x, field){

  .bt_meta_check_field(field)
  if(identical(field, "condition")){
    values <- lapply(.bt_meta_condition_names, function(name){
      attr(x, name, exact = TRUE)
    })
    names(values) <- .bt_meta_condition_names
    values <- values[!vapply(values, is.null, logical(1))]
    if(length(values) == 0L){
      return(NULL)
    }
    return(values)
  }

  attr(x, .bt_meta_attribute_names[[field]], exact = TRUE)
}

# 'x' with the draw-metadata field set to 'value'; NULL removes the field.
.bt_meta_set <- function(x, field, value){

  .bt_meta_check_field(field)
  if(identical(field, "condition")){
    if(!is.null(value) &&
       (!is.list(value) || is.null(names(value)) ||
        !all(names(value) %in% .bt_meta_condition_names))){
      stop("Draw conditioning metadata must be a named list of condition fields.",
           call. = FALSE)
    }
    for(name in .bt_meta_condition_names){
      attr(x, name) <- value[[name]]
    }
    return(x)
  }

  attr(x, .bt_meta_attribute_names[[field]]) <- value
  x
}

# One element of the 'condition' field.
.bt_meta_condition <- function(x, name){

  if(!name %in% .bt_meta_condition_names){
    stop("'", name, "' is not a draw-conditioning field.", call. = FALSE)
  }
  .bt_meta_get(x, "condition")[[name]]
}

# 'to' with every draw-metadata field of 'from'.
.bt_meta_copy <- function(to, from){

  for(field in .bt_meta_fields()){
    to <- .bt_meta_set(to, field, .bt_meta_get(from, field))
  }
  to
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
