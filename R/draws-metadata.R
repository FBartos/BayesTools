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
  "models_ind", "sample_ind", "ordered_total_indicator",
  "ordered_total_component", "undefined_draws",
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
  "effective_conditional_rule", "averaged"
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
      .posterior_atoms_from_attribute(value)
      NULL
    },
    components = function(value){
      if(inherits(value, "BayesTools_posterior_components")) NULL else
        "it must be a posterior component table"
    },
    component = function(value){
      if(.bt_meta_is_index(value)) NULL else
        "it must be a vector of positive integer component indices"
    },
    component_source = function(value){
      if(is.character(value) && length(value) == 1L &&
         value %in% .bt_meta_component_sources){
        return(NULL)
      }
      "it must be 'model', 'mixture', or 'spike_and_slab'"
    },
    draw_index = function(value){
      if(.bt_meta_is_index(value)) NULL else
        "it must be a vector of positive integer draw indices"
    },
    ordered_total_component = function(value){
      if(.bt_meta_is_index(value, allow_NA = TRUE)) NULL else
        "it must hold positive integer component indices (NA for models without a total spike)"
    },
    ordered_source = .bt_ordered_source_validate,
    formula_state = .bt_formula_state_validate,
    measure_unavailable = function(value){
      columns <- c("column", "measure", "reason")
      valid_cause <- !"cause" %in% names(value) || (is.character(value$cause) &&
        all(is.na(value$cause) | value$cause %in% .bt_formula_measure_causes))
      valid_diagnostics <- !"diagnostics" %in% names(value) || (is.list(value$diagnostics) &&
        all(vapply(value$diagnostics, .bt_formula_measure_diagnostics_valid, logical(1))))
      if(is.data.frame(value) && identical(class(value), "data.frame") &&
         identical(names(value)[seq_along(columns)], columns) &&
         all(names(value) %in% c(columns, "cause", "diagnostics")) &&
         !anyDuplicated(names(value)) && valid_cause && valid_diagnostics &&
         all(vapply(value[columns], function(x) is.character(x) && !anyNA(x) && all(nzchar(x)), logical(1))) &&
         all(value$measure %in% c("prior_density", "atoms", "support")) &&
         !anyDuplicated(value[c("column", "measure")])) return(NULL)
      "it must be a plain column/measure/reason table with unique target measures"
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
      .model_probability_context_validate(value)
      if(.bt_meta_is_prior_context(value)) NULL else
        "it must be a prior-density context"
    },
    model_probabilities = function(value){
      if(!is.list(value) || !identical(names(value), c("prior", "posterior"))) return("it must contain paired prior and posterior model probabilities")
      for(pair in value){
        if(!is.list(pair) || !identical(names(pair), c("probabilities", "logs", "declaration"))) return("its model probability pairs are malformed")
        .model_probability_validate(pair$probabilities, pair$logs, pair$declaration, normalized = TRUE)
      }
      NULL
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
      if(!is.null(value[["averaged"]]) &&
         !(isTRUE(value[["averaged"]]) || isFALSE(value[["averaged"]]))){
        return("'averaged' must be TRUE or FALSE")
      }
      NULL
    },
    linear_weights = function(value){
      if(is.numeric(value) && all(is.finite(value))) NULL else "it must be a finite numeric vector or matrix"
    },
    linear_weight_space = function(value){
      if(is.character(value) && length(value) == 1L && !is.na(value) &&
         value %in% c("coefficient", "formula_contribution")) NULL else
        "it must be 'coefficient' or 'formula_contribution'"
    },
    linear_offset = function(value){
      if(is.numeric(value) && length(value) == 1L) NULL else
        "it must be a number"
    },
    hypothesis_evaluation = function(value){
      if(!is.list(value) || !is.numeric(value$numerator) || any(!is.finite(value$numerator)) ||
         !is.numeric(value$divisor) || length(value$divisor) != 1L ||
         !is.finite(value$divisor) || value$divisor <= 0 ||
         !is.numeric(value$weights) || is.null(names(value$weights)) ||
         any(!is.finite(value$weights)) || !is.numeric(value$offset) ||
         length(value$offset) != 1L || !is.finite(value$offset)) return("invalid affine evaluation coordinates")
      NULL
    },
    joint_prior_transformation = function(value){
      if(is.character(value) && length(value) == 1L) NULL else
        "it must be a character string"
    },
    quantities = function(value){
      .bt_meta_quantities_reason(value)
    },
    original_scale_quantities = function(value){
      .bt_meta_quantities_reason(value)
    },
    level_quantities = function(value){
      .bt_meta_quantities_reason(value)
    },
    output_transformations = function(value){
      if(is.character(value) && length(value) > 0L && !anyNA(value) &&
         all(value %in% .bt_label_output_transformations)){
        return(NULL)
      }
      paste0(
        "it must name the transformations applied to the draws (",
        paste0("'", .bt_label_output_transformations, "'", collapse = ", "), ")"
      )
    }
  )
}

# The internal 'level_quantities' field of factor draws holds the catalog
# quantities of the level cells of the term (a column table whose columns are
# the cells' canonical names), declared by the producers of mixed posteriors
# so that the transformed levels (transform_factor_samples()) name the catalog
# quantity holding their values.
#
# The 'quantities' field: one row per draw column (the element itself for
# vector draws) with the column name, the catalog quantity id ("" for columns
# that are not catalog quantities), the fitted coordinates the column is a
# linear function of ('dependencies', with 'weights'), and the column's label
# parts. Consumers map columns to fitted coordinates and render labels from
# it, never from the column names. The quantity id and label parts describe
# the values the column holds; the internal 'original_scale_quantities' field
# (the same table, for the columns of fitted-scale mixed draws that are other
# quantities once transformed to the original predictor scale) replaces them
# when the draws are transformed.
.bt_meta_quantities_columns <- c(
  "column", "quantity_id", "dependencies", "weights", "label_parts"
)

.bt_meta_quantities_reason <- function(value){

  if(!is.data.frame(value) ||
     !identical(names(value), .bt_meta_quantities_columns)){
    return(paste0(
      "it must be a column table with the columns ",
      paste0("'", .bt_meta_quantities_columns, "'", collapse = ", ")
    ))
  }
  if(!is.character(value$column) || anyNA(value$column) ||
     any(!nzchar(value$column)) || anyDuplicated(value$column) ||
     !is.character(value$quantity_id) || anyNA(value$quantity_id)){
    return("its 'column' and 'quantity_id' must be character (unique, non-missing column names)")
  }
  valid_rows <- is.list(value$dependencies) && is.list(value$weights) &&
    is.list(value$label_parts) &&
    all(vapply(seq_len(nrow(value)), function(i){
      dependencies <- value$dependencies[[i]]
      weights <- value$weights[[i]]
      is.character(dependencies) && !anyNA(dependencies) &&
        is.numeric(weights) && !anyNA(weights) &&
        length(weights) == length(dependencies) &&
        inherits(value$label_parts[[i]], "BayesTools_label_parts")
    }, logical(1)))
  if(!valid_rows){
    return("its 'dependencies', 'weights', and 'label_parts' must describe every column")
  }

  NULL
}

# Sources of the per-draw component index: the model of a model-averaged
# ensemble (mix_posteriors()), or the component of a mixture or spike-and-slab
# prior of a single fit (as_mixed_posteriors()).
.bt_meta_component_sources <- c("model", "mixture", "spike_and_slab")

.bt_meta_is_index <- function(x, allow_NA = FALSE){

  if(is.logical(x) && allow_NA && all(is.na(x))){
    return(TRUE)
  }
  if(!is.numeric(x) || (!allow_NA && anyNA(x))){
    return(FALSE)
  }
  x <- x[!is.na(x)]
  all(x >= 1 & x == round(x))
}

.bt_meta_is_prior_context <- function(x){

  inherits(x, "prior_density_context") ||
    inherits(x, "prior_density_model_mixture_context") ||
    inherits(x, "prior_density_conditional_context")
}

.bt_meta_field_names <- c(
  "support", "atoms", "components", "component", "component_source",
  "draw_index", "ordered_total_component", "ordered_source", "undefined_draws", "prior_density",
  "formula_state", "measure_unavailable", "model_probabilities",
  "prior_context", "prior_densities", "posterior_density",
  "posterior_densities", "posterior_ordinate", "posterior_ordinates",
  "formula_parameter", "log_intercept", "formula_scale", "transform_scaled",
  "condition", "linear_weights", "linear_weight_space", "linear_offset", "joint_prior_transformation", "hypothesis_evaluation",
  "quantities", "original_scale_quantities", "level_quantities",
  "output_transformations"
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
     is.null(names(meta)) ||
     !all(names(meta) %in% c(.bt_meta_fields(), .bt_meta_fingerprint_field))){
    stop("The draw-metadata container is invalid.", call. = FALSE)
  }
  meta
}

# The value of a draw-metadata field, or NULL when it is absent. The metadata
# of draws must describe their current values (.bt_meta_check_current()).
.bt_meta_get <- function(x, field){

  .bt_meta_check_field(field)
  value <- .bt_meta_current_container(x)[[field]]
  if(!is.null(value) && field %in% c("formula_state", "ordered_source", "prior_context", "measure_unavailable", "model_probabilities")){
    .bt_validate_once(paste0("draw_metadata_", field), value, function() .bt_meta_validate(field, value))
  }
  value
}

# The values of several draw-metadata fields (a list named by 'fields', NULL
# for absent fields), with one check that the metadata describe the draws.
.bt_meta_get_fields <- function(x, fields){

  for(field in fields){
    .bt_meta_check_field(field)
  }
  meta <- .bt_meta_current_container(x)
  stats::setNames(lapply(fields, function(field){
    value <- meta[[field]]
    if(!is.null(value) && field %in% c("formula_state", "ordered_source", "prior_context", "measure_unavailable", "model_probabilities")){
      .bt_validate_once(paste0("draw_metadata_", field), value, function() .bt_meta_validate(field, value))
    }
    value
  }), fields)
}

# The container of 'x' (NULL when 'x' carries none), checked to describe the
# current values of draws.
.bt_meta_current_container <- function(x){

  meta <- .bt_meta_container(x)
  if(!is.null(meta) && identical(meta$component_source, "model") && is.null(meta$model_probabilities)){
    .bt_stop_refit_required("Model probability ownership is missing. Recreate these mixed posteriors from the original fits with the current BayesTools version.")
  }
  if(!is.null(meta) && .bt_meta_is_draws(x)){
    .bt_meta_check_current(meta, .bt_meta_fingerprint(x))
  }
  meta
}

# 'x' with the draw-metadata field set to 'value'; NULL removes the field.
# Setting a field of draws whose values changed after their metadata were
# attached stops: producers that change the values of draws rebuild the
# container first (.bt_meta_refresh()).
.bt_meta_set <- function(x, field, value){

  .bt_meta_assign(x, stats::setNames(list(value), field))
}

# 'x' with the fields of the named list 'fields' set (NULL removes a field),
# with one check and one fingerprint of the draws.
.bt_meta_assign <- function(x, fields){

  if("linear_weights" %in% names(fields) && !"linear_weight_space" %in% names(fields)){
    context <- fields$prior_context
    if(is.null(context)) context <- .bt_meta_get(x, "prior_context")
    fields["linear_weight_space"] <- list(if(is.null(fields$linear_weights)) NULL else
      if(identical(context$linear_weight_space, "formula_contribution")) "formula_contribution" else "coefficient")
  }
  if(!is.null(fields$measure_unavailable)){
    .bt_meta_validate("measure_unavailable", fields$measure_unavailable)
    if(nrow(fields$measure_unavailable) == 0L) fields["measure_unavailable"] <- list(NULL)
  }
  for(field in names(fields)){
    .bt_meta_check_field(field)
    if(!is.null(fields[[field]])){
      .bt_meta_validate(field, fields[[field]])
      if(.bt_meta_is_draws(x) && field %in% c("component", "draw_index", "ordered_total_component", "ordered_source", "formula_state")){
        rows <- if(identical(field, "ordered_source")) nrow(fields[[field]]$primitives) else
          if(identical(field, "formula_state")) nrow(fields[[field]]$values) else NROW(fields[[field]])
        if(rows != NROW(x)) stop("Draw metadata '", field, "' must have one row per draw.", call. = FALSE)
      }
      if(.bt_meta_is_draws(x) && identical(field, "measure_unavailable")){
        quantities <- .bt_meta_get(x, "quantities")
        columns <- if(is.matrix(x)) colnames(x) else
          if(!is.null(quantities) && nrow(quantities) == 1L) quantities$column else attr(x, "parameter", exact = TRUE)
        if(is.null(columns)) columns <- .bt_meta_get(x, "quantities")$column
        if(!is.null(columns) && any(!fields[[field]]$column %in% columns)) stop(
          "Unavailable measure keys must identify current draw columns.", call. = FALSE)
      }
      if(.bt_meta_is_draws(x) && identical(field, "hypothesis_evaluation") &&
         !identical(as.numeric(x), fields[[field]]$numerator / fields[[field]]$divisor)){
        stop("Affine hypothesis evaluation coordinates do not match the reported draws.", call. = FALSE)
      }
    }
  }
  meta <- .bt_meta_container(x)
  if(is.null(meta) && all(vapply(fields, is.null, logical(1)))){
    return(x)
  }
  fingerprint <- if(.bt_meta_is_draws(x)) .bt_meta_fingerprint(x)
  if(is.null(meta)){
    meta <- structure(list(), names = character(), class = c("BayesTools_draw_metadata", "list"))
  }else if(!is.null(fingerprint)){
    .bt_meta_check_current(meta, fingerprint)
  }
  for(field in names(fields)){
    meta[[field]] <- fields[[field]]
  }
  .bt_meta_write(x, meta, fingerprint)
}

# One row-subsetting operation for every row-aligned metadata field. Producers
# that change values separately retain fitted ordered primitives through this
# helper, rather than copying an old row count onto new draws.
.bt_draws_subset_rows <- function(x, rows){

  meta <- .bt_meta_current_container(x)
  out <- if(is.null(dim(x))) .bt_draws_plain(x)[rows] else .bt_draws_plain(x)[rows, , drop = FALSE]
  attributes_to_keep <- attributes(x)[setdiff(names(attributes(x)), c("names", "dim", "dimnames", .bt_meta_attribute))]
  attributes(out) <- c(attributes(out), attributes_to_keep)
  if(is.null(meta)) return(out)
  for(field in c("component", "draw_index", "ordered_total_component")){
    if(!is.null(meta[[field]])) meta[[field]] <- if(is.null(dim(meta[[field]]))) meta[[field]][rows] else meta[[field]][rows, , drop = FALSE]
  }
  if(!is.null(meta$ordered_source)) meta$ordered_source <- .bt_ordered_source_subset(meta$ordered_source, rows)
  if(!is.null(meta$formula_state)) meta$formula_state <- .bt_formula_state_subset(meta$formula_state, rows)
  if(!is.null(meta$hypothesis_evaluation)){
    meta$hypothesis_evaluation$numerator <- meta$hypothesis_evaluation$numerator[rows]
  }
  if(!is.null(meta$components) && length(meta$components$index)==NROW(x)){
    meta$components <- .posterior_components_new(meta$components$index[rows],meta$components$supports,meta$components$keys)
  }
  if(is.matrix(meta$linear_weights) && nrow(meta$linear_weights)==NROW(x)){
    meta$linear_weights <- meta$linear_weights[rows,,drop=FALSE]
  }
  out <- .bt_meta_write(out, meta, if(.bt_meta_is_draws(out)) .bt_meta_fingerprint(out))
  source <- meta$ordered_source
  if(!is.null(source$projection_design)){
    out <- .bt_meta_set(out,"atoms",NULL)
    columns <- if(is.null(dim(out))) rownames(source$projection_design) else colnames(out)
    out <- .bt_ordered_source_semantics(out,diag(length(columns)),columns)
  }
  out
}

# Draw metadata go stale when the values of the draws change after the
# metadata were attached (e.g., by 'x[] <- ', 'x[i] <- ', or pmin(), which
# keep the attributes of 'x'). The container of draws therefore stores a
# fingerprint of the values it describes: their number, the number of missing
# values, and the sums of the observed values and of the observed values
# weighted by their positions. The fingerprint is computed in one native pass
# (src/r-draw-fingerprint.c) or, when the package's native routines are not
# loaded (the package loads without its DLL when JAGS cannot be located), by
# an R evaluator of the same definition (.bt_meta_fingerprint_r()).
.bt_meta_fingerprint_field <- "fingerprint"

.bt_meta_is_draws <- function(x){

  is.atomic(x) && typeof(x) %in% c("double", "integer", "logical")
}

# The fingerprint of the values of draws 'x' ('value', stored in the
# container) with the scales of its rounding-error bound ('scale').
.bt_meta_fingerprint <- function(x){

  out <- if(.bt_meta_fingerprint_native()){
    .Call("BayesTools_draw_fingerprint", x, PACKAGE = "BayesTools")
  }else{
    .bt_meta_fingerprint_r(x)
  }
  list(
    value = c(length = out[[1L]], missing = out[[2L]], sum = out[[3L]], weighted_sum = out[[4L]]),
    scale = c(sum = out[[5L]], weighted_sum = out[[6L]])
  )
}

# Whether draw fingerprints are computed by the native pass: the native
# routines are loaded and the R evaluator is not forced (the internal switch
# .BayesTools_private$draw_fingerprint_r, used by the tests).
.bt_meta_fingerprint_native <- function(){

  !isTRUE(.BayesTools_private$draw_fingerprint_r) &&
    isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))
}

# The R evaluator of the native fingerprint pass: the same six values, with
# the same fixed summation order. Value i goes to lane i mod 4 (row of a
# 4-row matrix, padded with zeros, which add nothing to a lane sum),
# rowSums() adds the observed values of each lane in increasing i, and the
# lanes are combined as (p0 + p1) + (p2 + p3). The lane sums accumulate in
# long double where R has it, so they can differ from the native double lanes
# in the last bits, within the rounding-error bound of
# .bt_meta_fingerprint_matches(); they are identical where long double is
# double, and for integer and logical draws (exact sums). A lane sum that
# overflows double precision only transiently (draws near the largest
# double) can be infinite in one evaluator and finite in the other.
.bt_meta_fingerprint_r <- function(x){

  if(!.bt_meta_is_draws(x)){
    stop("Draw fingerprints require numeric or logical values.", call. = FALSE)
  }
  n <- length(x)
  values <- as.double(x)
  any_missing <- anyNA(values)
  padded <- 4 * ceiling(n / 4)
  lanes <- c(values, rep(0, padded - n))
  dim(lanes) <- c(4L, padded / 4)
  weighted <- lanes * seq_len(padded)
  combine <- function(lane_values){
    lane_sums <- rowSums(lane_values, na.rm = any_missing)
    (lane_sums[[1L]] + lane_sums[[2L]]) + (lane_sums[[3L]] + lane_sums[[4L]])
  }
  c(
    n,
    if(any_missing) sum(is.na(values)) else 0,
    combine(lanes),
    combine(weighted),
    combine(abs(lanes)),
    combine(abs(weighted))
  )
}

# Whether a stored fingerprint describes the values of the current one. The
# number of values and of missing values must be equal. The native pass sums
# in double precision in a fixed order, four interleaved partial sums of at
# most m = ceiling(n / 4) summands combined pairwise, so a computed sum is
# within gamma_k * sum(|t|) of the exact one, with k = m + 2 (the weighted
# sum's products are rounded too), gamma_k = k u / (1 - k u), u = 2^-53, and
# t the summands. Two computations of the same values - on the same build
# bit-identical, on other platforms possibly with fused multiply-adds - thus
# agree within 2 gamma_k (1 + gamma_k) times the computed sum(|t|), the last
# factor bounding the rounding of that scale. Differences within this bound
# (for example changes of single draws far below their magnitude) are not
# detected.
.bt_meta_fingerprint_matches <- function(stored, current){

  if(!is.numeric(stored) || !identical(names(stored), names(current$value)) ||
     !identical(stored[c("length", "missing")], current$value[c("length", "missing")])){
    return(FALSE)
  }
  if(identical(stored, current$value)){
    return(TRUE)
  }
  k <- ceiling(current$value[["length"]] / 4) + 2
  gamma <- k * .Machine$double.eps / 2
  gamma <- gamma / (1 - gamma)
  all(vapply(c("sum", "weighted_sum"), function(name){
    stored_sum  <- stored[[name]]
    current_sum <- current$value[[name]]
    if(!is.finite(stored_sum) || !is.finite(current_sum)){
      return(identical(stored_sum, current_sum))
    }
    abs(stored_sum - current_sum) <= 2 * gamma * (1 + gamma) * current$scale[[name]]
  }, logical(1)))
}

.bt_meta_check_current <- function(meta, current){

  if(!.bt_meta_fingerprint_matches(meta[[.bt_meta_fingerprint_field]], current)){
    stop(errorCondition(
      paste0(
        "The metadata of these posterior draws are unavailable: their values ",
        "changed after the metadata were attached (for example by 'x[] <- ', ",
        "'x[i] <- ', or 'pmin()'), so the supports, atoms, and prior densities ",
        "no longer describe them. Use 'posterior_transform()' (or ",
        "'marginal_posterior(transformation = )') for transformed posterior ",
        "distributions."
      ),
      class = c("BayesTools_stale_metadata", "BayesTools_metadata"),
      call = NULL
    ))
  }
  invisible(TRUE)
}

# 'x' with the container 'meta' (and, for draws, the fingerprint of their
# values); a container without fields is removed.
.bt_meta_write <- function(x, meta, fingerprint = NULL){

  meta[[.bt_meta_fingerprint_field]] <- NULL
  if(length(meta) == 0L){
    attr(x, .bt_meta_attribute) <- NULL
    return(x)
  }
  if(!is.null(fingerprint)){
    meta[[.bt_meta_fingerprint_field]] <- fingerprint$value
  }
  attr(x, .bt_meta_attribute) <- meta
  x
}

# Rebuilds the container of draws whose values a producer changed together
# with the metadata that describe them (the producer transforms, replaces, or
# removes the affected fields).
.bt_meta_refresh <- function(x){

  meta <- .bt_meta_container(x)
  if(is.null(meta)){
    return(x)
  }
  .bt_meta_write(x, meta, if(.bt_meta_is_draws(x)) .bt_meta_fingerprint(x))
}

# 'x' with several draw-metadata fields set, given as named arguments (one
# check and one fingerprint of the draws).
.bt_meta_update <- function(x, ...){

  fields <- list(...)
  if(length(fields) == 0L){
    return(x)
  }
  if(is.null(names(fields)) || any(!nzchar(names(fields)))){
    stop("Draw-metadata fields must be named.", call. = FALSE)
  }
  .bt_meta_assign(x, fields)
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
  "support", "atoms", "undefined_draws", "prior_density", "prior_densities",
  "prior_context", "posterior_density", "posterior_densities",
  "posterior_ordinate", "posterior_ordinates", "condition", "linear_weights",
  "quantities", "output_transformations", "measure_unavailable", "linear_weight_space",
  "ordered_source"
)

.bt_formula_state_validate <- function(value){

  if(!is.list(value) || !identical(names(value), c("schema_version", "models", "model", "draw_index", "values", "posterior_model_probabilities", "posterior_log_model_probabilities")) ||
     !identical(value$schema_version, 2L) || !is.list(value$models) || !length(value$models) ||
     !.bt_meta_is_index(value$model) || !.bt_meta_is_index(value$draw_index) ||
     length(value$model) != length(value$draw_index) ||
     !is.matrix(value$values) || !is.numeric(value$values) ||
     nrow(value$values) != length(value$model) ||
     (ncol(value$values) > 0L && is.null(colnames(value$values))) || anyDuplicated(colnames(value$values)) ||
     any(value$model > length(value$models)) || !is.numeric(value$posterior_model_probabilities) ||
     length(value$posterior_model_probabilities) != length(value$models) ||
     any(!is.finite(value$posterior_model_probabilities)) || any(value$posterior_model_probabilities < 0) ||
     abs(sum(value$posterior_model_probabilities) - 1) > 1e-12){
    return("it must be a versioned formula state with aligned rows and posterior model probabilities")
  }
  declaration <- attr(value, "model_probability_declaration", exact = TRUE)
  .model_probability_validate(value$posterior_model_probabilities,
    value$posterior_log_model_probabilities, declaration$posterior, normalized = TRUE)
  prior_probs <- vapply(value$models, `[[`, numeric(1), "prior_probability")
  prior_logs <- vapply(value$models, `[[`, numeric(1), "prior_log_probability")
  .model_probability_validate(prior_probs, prior_logs, declaration$prior, normalized = TRUE)
  for(model in seq_along(value$models)){
    record <- value$models[[model]]
    if(!is.list(record) || !is.list(record$prior_list) || !is.list(record$formula_scale) ||
       !is.character(record$required) || anyNA(record$required) || anyDuplicated(record$required)){
      return("formula model declarations are incomplete")
    }
    rows <- value$model == model
    if(any(rows) && (!all(record$required %in% colnames(value$values)) ||
                    any(!is.finite(value$values[rows, record$required, drop = FALSE])))){
      return("required own-model formula states must be finite and present")
    }
  }
  NULL
}

.bt_formula_state_subset <- function(value, rows){

  value$model <- value$model[rows]
  value$draw_index <- value$draw_index[rows]
  value$values <- value$values[rows, , drop = FALSE]
  value
}

.bt_formula_measure_diagnostics_valid <- function(value){

  if(is.null(value)) return(TRUE)
  if(!is.list(value) || is.object(value) || !identical(names(value),
    c("model_indices", "log_prior_probabilities", "log_posterior_probabilities", "stage"))) return(FALSE)
  indices <- value$model_indices
  valid_logs <- function(logs){
    is.null(logs) || (is.numeric(logs) && length(logs) == length(indices) &&
      !anyNA(logs) && all(is.finite(logs) | logs == -Inf))
  }
  is.integer(indices) && length(indices) > 0L && !anyNA(indices) &&
    all(indices > 0L) && !anyDuplicated(indices) &&
    (!is.null(value$log_prior_probabilities) || !is.null(value$log_posterior_probabilities)) &&
    valid_logs(value$log_prior_probabilities) && valid_logs(value$log_posterior_probabilities) &&
    is.character(value$stage) && length(value$stage) == 1L &&
    !is.na(value$stage) && value$stage %in% c("prior", "posterior", "conditional_prior", "conditional_posterior", "event", "component")
}

.bt_formula_measure_check <- function(x, measure, column = NULL){

  unavailable <- .bt_meta_get(x, "measure_unavailable")
  if(is.null(unavailable)) return(invisible(TRUE))
  if(is.null(column)){
    quantities <- .bt_meta_get(x, "quantities")
    column <- if(is.matrix(x)) colnames(x) else
      if(!is.null(quantities) && nrow(quantities) == 1L) quantities$column else attr(x, "parameter", exact = TRUE)
    if(is.null(column)){
      quantities <- .bt_meta_get(x, "quantities")
      column <- quantities$column
    }
  }
  selected <- unavailable$measure == measure & unavailable$column %in% column
  if(!any(selected)) return(invisible(TRUE))
  entry <- unavailable[which(selected)[1L], , drop = FALSE]
  cause <- if("cause" %in% names(entry) && !is.na(entry$cause)) entry$cause else entry$reason
  diagnostics <- if("diagnostics" %in% names(entry)) entry$diagnostics[[1L]] else NULL
  if(identical(measure, "prior_density")) .bt_formula_density_stop(
    paste0("Prior density for '", entry$column, "' is unavailable: ", entry$reason, "."),
    target = entry$column, reason = cause, detail = entry$reason, diagnostics = diagnostics)
  stop(errorCondition(paste0("Formula ", measure, " for '", entry$column,
    "' are unavailable: ", entry$reason, "."), call = NULL,
    class = c(paste0("BayesTools_formula_", measure, "_unavailable"),
      if(cause %in% .bt_formula_measure_causes) "BayesTools_formula_measure_unavailable"),
    target = entry$column, reason = cause, detail = entry$reason, diagnostics = diagnostics))
}

.bt_formula_measure_mark <- function(x, column, measure, reason, cause = NULL, diagnostics = NULL){

  unavailable <- .bt_meta_get(x, "measure_unavailable")
  entry <- data.frame(column = column, measure = measure, reason = reason, stringsAsFactors = FALSE)
  if(!is.null(cause)) entry$cause <- cause
  if(!is.null(diagnostics)) entry$diagnostics <- list(diagnostics)
  if(!is.null(unavailable) && "cause" %in% names(unavailable) && !"cause" %in% names(entry)) entry$cause <- NA_character_
  if(!is.null(unavailable) && "cause" %in% names(entry) && !"cause" %in% names(unavailable)) unavailable$cause <- NA_character_
  if(!is.null(unavailable) && "diagnostics" %in% names(unavailable) && !"diagnostics" %in% names(entry)) entry$diagnostics <- list(NULL)
  if(!is.null(unavailable) && "diagnostics" %in% names(entry) && !"diagnostics" %in% names(unavailable)) unavailable$diagnostics <- rep(list(NULL), nrow(unavailable))
  if(!is.null(unavailable)) unavailable <- unavailable[!(unavailable$column == column & unavailable$measure == measure), , drop = FALSE]
  .bt_meta_set(x, "measure_unavailable", rbind(unavailable, entry))
}

.bt_formula_measure_linear <- function(source, target, design){

  unavailable <- .bt_meta_get(source, "measure_unavailable")
  if(is.null(unavailable)) return(target)
  source_columns <- colnames(source)
  target_columns <- colnames(target)
  if(ncol(design) != length(source_columns) || nrow(design) != length(target_columns)){
    stop("Unavailable measure columns do not align with their declared linear design.", call. = FALSE)
  }
  target <- .bt_meta_set(target, "measure_unavailable", NULL)
  for(row in seq_len(nrow(design))){
    active <- source_columns[design[row, ] != 0]
    entries <- unavailable[unavailable$column %in% active, , drop = FALSE]
    if(nrow(entries)) for(measure in unique(entries$measure)){
      target <- .bt_formula_measure_mark(target, target_columns[[row]], measure,
        entries$reason[entries$measure == measure][[1L]],
        cause = if("cause" %in% names(entries)) entries$cause[entries$measure == measure][[1L]],
        diagnostics = if("diagnostics" %in% names(entries)) entries$diagnostics[entries$measure == measure][[1L]])
    }
  }
  target
}

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
#'   \item{\code{"prior_densities"}}{a list of prior densities keyed by
#'   parameter, attached to the list returned by [as_mixed_posteriors()]
#'   with \code{transform_scaled = TRUE} for the requested parameters.}
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
#'   \code{condition_event}, \code{resolved_condition_event},
#'   \code{averaged} (\code{TRUE} when the draws are not conditioned on any
#'   event, i.e., the unconditional model-averaged posterior; \code{FALSE}
#'   otherwise), and, for levels whose conditioning was resolved per level,
#'   \code{effective_conditional} and \code{effective_conditional_rule}.
#'   Read \code{averaged} instead of comparing \code{condition_key} with a
#'   literal key.}
#'   \item{\code{"measure_unavailable"}}{a plain data frame with exact columns
#'   \code{column}, \code{measure}, and \code{reason}, with an optional declared
#'   \code{cause} enum; absent causes use an already enumerated reason.
#'   An optional \code{diagnostics} list column retains compact model probability
#'   diagnostics. Each non-NULL cell is a plain list with \code{model_indices},
#'   \code{log_prior_probabilities}, \code{log_posterior_probabilities}, and
#'   \code{stage}; supplied log vectors align with positive unique integer
#'   model indices. A missing log vector is NULL. Stage identifies prior,
#'   posterior, conditional prior/posterior, event, or component probabilities.
#'   The cause \code{numerical_model_probability_unavailable} identifies a
#'   declared law whose active model probabilities cannot be represented at
#'   full precision. Invalid metadata requires recomputation/refitting.
#'   Human reason text alone does not classify a measure refusal. Tables use
#'   unique column/measure pairs. The measures are \code{prior_density}, \code{atoms}, or \code{support}.
#'   Entries prevent stale prior/density/support fallback while numeric draws
#'   remain usable. Atom/support refusals use
#'   \code{BayesTools_formula_atoms_unavailable} and
#'   \code{BayesTools_formula_support_unavailable}, inheriting
#'   \code{BayesTools_formula_measure_unavailable}.}
#'   \item{\code{"linear_weight_space"}}{\code{"coefficient"} for coefficient
#'   recipes, or \code{"formula_contribution"} for fitted design rows and compiled
#'   multipliers. Contribution weights never pass through coefficient C/A again.
#'   Old ambiguous stored contexts must be recreated.}
#'   \item{\code{"linear_weights"}}{the weights of the fitted coordinates
#'   that form a level of a [marginal_posterior()] (a named numeric vector,
#'   or a matrix with one row per draw).}
#'   \item{\code{"ordered_source"}}{internal retained fitted totals, normalized
#'   allocations, model provenance and source-row indices. Producers validate
#'   its row count and preserve the primitive values independently of effect
#'   transformations. Recreate missing sources from the fitted models; shares
#'   are never reconstructed from increments divided by totals.}
#'   \item{\code{"quantities"}}{the quantity of every column of mixed and
#'   marginal posterior draws: a data frame with one row per column (one row
#'   for vector draws) and the columns \code{column} (the column name),
#'   \code{quantity_id} (the [parameter_catalog()] quantity id, or
#'   \code{""} for columns that are not catalog quantities),
#'   \code{dependencies} and \code{weights} (the fitted coordinates the
#'   column is a linear function of), and \code{label_parts} (the label
#'   parts rendered by [parameter_labels()]). Set by [as_mixed_posteriors()],
#'   [mix_posteriors()], [transform_factor_samples()], and
#'   [marginal_posterior()] (whose estimated marginal means are predictions
#'   and declare no fitted coordinates); summaries map columns to fitted
#'   coordinates and render their labels from it.}
#'   \item{\code{"output_transformations"}}{the transformations applied to
#'   the values of the draws by [posterior_transform()] (and
#'   [marginal_posterior()] with \code{transformation}), in the order they
#'   were applied: \code{"lin"}, \code{"exp_lin"}, \code{"tanh"},
#'   \code{"exp"}, or \code{"custom"} for a transformation given as
#'   functions; absent for untransformed draws.}
#' }
#' @param value the new value of the field; \code{NULL} removes it.
#'
#' @details Arithmetic and mathematical functions of posterior draws
#' (\code{Ops} and \code{Math} group generics), \code{c()},
#' \code{as.numeric()}, and subsetting return plain numeric draws without
#' metadata, because supports, atoms, and prior densities do not follow the
#' transformed values. Functions that need the metadata (e.g.
#' [Savage_Dickey_BF()], [plot_posterior()], [marginal_posterior()]) stop on
#' such draws; [posterior_transform()] (and [marginal_posterior()] with
#' \code{transformation}) transforms the draws together with their metadata.
#'
#' Operations that replace values but keep the attributes of the draws
#' (\code{x[] <- }, \code{x[i] <- }, \code{pmin()}, \code{pmax()}) leave
#' metadata that no longer describe the draws. The metadata record a
#' fingerprint of the values they describe (their number, the number of
#' missing values, and two sums of the observed values), and reading or
#' setting the metadata of such draws, here or in any function that uses
#' them, stops with an error of class \code{BayesTools_stale_metadata}
#' (parent class \code{BayesTools_metadata}). Subsetting a list of mixed
#' posteriors with \code{[} keeps the list's metadata, with the prior
#' densities of the omitted parameters removed.
#'
#' The fingerprint cannot detect every change: the two sums are compared
#' within their rounding-error bound, so changes that keep both sums (e.g.,
#' adding d, -2d, and d to three consecutive draws), changes of single draws
#' smaller than about \eqn{5.6 \cdot 10^{-17} n^2} times the draws' mean
#' absolute value for \eqn{n} draws (about \eqn{6 \cdot 10^{-5}} for
#' \eqn{10^6} draws), and reorderings that move draws by only a few positions
#' can keep the metadata. The bound grows with the number of draws: such
#' changes are detected at some positions or for fewer draws (e.g., swapping
#' two adjacent standard-normal draws is typically detected for \eqn{10^5}
#' draws but not for \eqn{10^6}).
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

# Plain numeric draws: 'x' without its class and metadata (dimensions and
# names are kept).
.bt_draws_plain <- function(x){

  keep <- intersect(names(attributes(x)), c("dim", "dimnames", "names"))
  attributes(x) <- attributes(x)[keep]
  x
}

# Arithmetic and mathematical functions of posterior draws return plain
# numeric draws: supports, atoms, prior densities and the other draw metadata
# do not follow the transformed values. posterior_transform() transforms the
# draws together with their metadata.
.bt_draws_ops <- function(e1, e2){

  if(missing(e2)){
    return(get(.Generic)(.bt_draws_plain(e1)))
  }
  get(.Generic)(.bt_draws_plain(e1), .bt_draws_plain(e2))
}

.bt_draws_math <- function(x, ...){

  get(.Generic)(.bt_draws_plain(x), ...)
}

# One method object per group generic for every draw class, so that
# operations combining draws of different classes dispatch to it (R uses the
# internal method, which keeps attributes, when the two methods differ).
#' @exportS3Method Ops mixed_posteriors
Ops.mixed_posteriors <- .bt_draws_ops
#' @exportS3Method Ops marginal_posterior
Ops.marginal_posterior <- .bt_draws_ops
#' @exportS3Method Ops marginal_posterior.simple
Ops.marginal_posterior.simple <- .bt_draws_ops
#' @exportS3Method Ops marginal_posterior.factor
Ops.marginal_posterior.factor <- .bt_draws_ops
#' @exportS3Method Math mixed_posteriors
Math.mixed_posteriors <- .bt_draws_math
#' @exportS3Method Math marginal_posterior
Math.marginal_posterior <- .bt_draws_math
#' @exportS3Method Math marginal_posterior.simple
Math.marginal_posterior.simple <- .bt_draws_math
#' @exportS3Method Math marginal_posterior.factor
Math.marginal_posterior.factor <- .bt_draws_math

# Subsetting a list of mixed posteriors keeps the list and its metadata: the
# per-element prior densities are subset to the kept elements, and the joint
# fields (prior context, scaling, conditioning, posterior density sources)
# are kept. Subsetting draws returns plain numeric draws.
#' @exportS3Method "[" mixed_posteriors
`[.mixed_posteriors` <- function(x, ...){

  out <- NextMethod()
  if(!is.list(x)){
    return(out)
  }
  kept <- names(out)
  x_attributes <- attributes(x)
  x_attributes <- x_attributes[setdiff(names(x_attributes), c("names", .bt_meta_attribute))]
  x_attributes[["names"]] <- kept
  attributes(out) <- x_attributes

  meta <- .bt_meta_container(x)
  if(is.null(meta)){
    return(out)
  }
  prior_densities <- meta[["prior_densities"]]
  if(!is.null(prior_densities)){
    keep <- vapply(names(prior_densities), function(key){
      any(key == kept | startsWith(key, paste0(kept, "[")))
    }, logical(1))
    subset <- unclass(prior_densities)[keep]
    density_attributes <- attributes(prior_densities)
    density_attributes[["names"]] <- names(subset)
    attributes(subset) <- density_attributes
    meta[["prior_densities"]] <- if(any(keep)) subset
  }
  if(!is.null(meta$measure_unavailable)){
    keep <- vapply(meta$measure_unavailable$column, function(column){
      any(column == kept | startsWith(column, paste0(kept, "[")))
    }, logical(1))
    meta$measure_unavailable <- if(any(keep)) meta$measure_unavailable[keep, , drop = FALSE]
  }
  .bt_meta_write(out, meta)
}

.bt_linear_weight_space <- function(x){

  if(is.null(.bt_meta_get(x, "linear_weights"))) return(NULL)
  space <- .bt_meta_get(x, "linear_weight_space")
  if(is.null(space)) .bt_stop_refit_required(
    "Stored linear prior weights have no supported coordinate space. Recreate the posterior with this version of BayesTools.")
  context <- .bt_meta_get(x, "prior_context")
  if(length(context$transforms)) .bt_formula_context_check(context)
  if(identical(space, "formula_contribution") && (!is.null(context$formula_scale) && length(context$formula_scale))){
    .bt_stop_refit_required("Formula contribution weights cannot carry a coefficient unscaling map. Recreate the posterior with this version of BayesTools.")
  }
  space
}

# 'fun' applied to the values of draws 'x', keeping every attribute of 'x'
# (for producers that transform the metadata explicitly).
.bt_draws_transform_values <- function(x, fun){

  meta <- .bt_meta_container(x)
  if(!is.null(meta) && .bt_meta_is_draws(x)){
    .bt_meta_check_current(meta, .bt_meta_fingerprint(x))
  }
  out <- fun(.bt_draws_plain(x))
  attributes(out) <- attributes(x)
  .bt_meta_refresh(out)
}

# Error for plain numeric draws passed where BayesTools posterior draws with
# their metadata are required.
.bt_draws_stop_plain <- function(what){

  stop(
    what, ": arithmetic, mathematical functions, and subsetting of posterior ",
    "draws return plain numeric draws without their supports, atoms, and ",
    "prior densities. Use posterior_transform() (or ",
    "marginal_posterior(transformation = )) for transformed posterior ",
    "distributions.",
    call. = FALSE
  )
}

# Component indices (into the declared component list of a mixture or
# spike-and-slab prior) of fitted indicator draws: the 'dcat' index of a
# mixture, and for a spike-and-slab prior the 'dbern' inclusion indicator
# (1 = slab) mapped to the positions of its slab and spike components.
.bt_component_from_indicator <- function(prior, indicator){

  indicator <- as.numeric(indicator)
  if(anyNA(indicator)){
    stop("Mixture component indicator draws must not be missing.", call. = FALSE)
  }
  if(is.prior.spike_and_slab(prior)){
    if(any(!indicator %in% c(0, 1))){
      stop("Spike-and-slab indicator draws must be 0 or 1.", call. = FALSE)
    }
    components <- attr(prior, "components", exact = TRUE)
    slab  <- which(components == "alternative")
    spike <- which(components == "null")
    return(ifelse(indicator == 1, slab, spike))
  }
  if(is.prior.mixture(prior)){
    if(any(!indicator %in% seq_along(prior))){
      stop("Mixture indicator draws must index the mixture components.", call. = FALSE)
    }
    return(as.integer(indicator))
  }

  stop("Component indicators require a mixture or spike-and-slab prior.", call. = FALSE)
}

# Whether 'component' indexes the spike (the 'null' component) of a
# spike-and-slab prior.
.bt_component_is_spike <- function(prior, component){

  identical(attr(prior, "components", exact = TRUE)[component], "null")
}

# 'x' with its per-draw component index and the source of that index.
.bt_draws_set_component <- function(x, component, source){

  .bt_meta_update(x, component = as.integer(component), component_source = source)
}

# The per-draw component index of draws. By explicit rule, draws without a
# mixture (parameters of a single fit without mixture priors) form one
# component.
.bt_draws_component <- function(x, n = NROW(x)){

  component <- .bt_meta_get(x, "component")
  if(is.null(component)){
    return(rep(1L, n))
  }
  component
}

# The model of each draw of a model-averaged ensemble (mix_posteriors()), or
# NULL for draws of a single fit.
.bt_draws_model_component <- function(x){

  if(!identical(.bt_meta_get(x, "component_source"), "model")){
    return(NULL)
  }
  .bt_meta_get(x, "component")
}
