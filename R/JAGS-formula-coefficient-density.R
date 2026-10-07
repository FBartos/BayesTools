# Versioned fixed-coefficient transforms and induced prior densities.

.bt_formula_coefficient_transform_version <- 3L

#' Formula coefficient transformations and induced prior densities
#'
#' @description
#' `JAGS_formula_coefficient_transform()` exposes the exact fixed-coefficient
#' transformation used by [transform_scale_samples()]. It records fitted
#' sources, original-scale targets, named source/output transforms, exact
#' nonzero dependencies, and parameter-map structural metadata.
#'
#' `JAGS_formula_prior_density()` constructs the induced marginal prior measure
#' for one target, or for a weighted combination of targets, through
#' BayesTools' deterministic prior-density context and linear-density algebra.
#' The result can be passed directly to [prior_density_ordinate()]. Fitted
#' coefficient priors are raw; `multiply_by` belongs to their compiled formula
#' contribution. Original raw coefficients use
#' \eqn{z'_j = \sum_k A_{jk} z_k m_k / m_j} on one immutable fitted state.
#' Identical canonical multipliers cancel before division, including zero;
#' declared zero contributions require no irrelevant state. A different needed
#' zero denominator refuses the requested numeric vector, without row repair.
#'
#' Finite `matrix` rows are static raw maps. Dynamic rows are entirely `NA` and
#' retain finite `basis_matrix`, canonical multipliers, true state dependencies,
#' and a `raw_affine`, `contribution_affine`, or `unavailable` prior recipe.
#' Compatible weighted requests share one context; declared point siblings may
#' add a constant offset. Incompatible recipes refuse with
#' `"incompatible_prior_recipes"`; unsupported stochastic denominators refuse
#' with `"state_dependent_map"`. Unsupported dependence/products remain explicit.
#'
#' `JAGS_formula_internal_coordinate_priors()` returns exact scalar priors for
#' stochastic formula coordinates that are intentionally absent from the
#' ordinary fitted `prior_list`. Currently these are the independent beta
#' primitives used by compiled LKJ random-effect blocks.
#'
#' @param fit fitted object created by [JAGS_fit()].
#' @param parameter scalar formula parameter name.
#' @param target_scale requested coefficient scale. The current schema supports
#'   only `"original"`.
#' @param target exact target coordinate from the transformation object.
#'   Supply exactly one of `target` and `weights`.
#' @param context optional BayesTools prior-density context. This is useful for
#'   product-space or conditional prior mixtures; when `NULL`, the fitted
#'   `prior_list` is used.
#' @param weights named numeric vector of weights over `target_names` of the
#'   transformation object; the density is that of the weighted combination
#'   of those targets. Supply exactly one of `target` and `weights`.
#'
#' @return `JAGS_formula_coefficient_transform()` returns a
#' `BayesTools_formula_coefficient_transform` list (schema version 3). Its
#' `targets` data frame has one row per original-scale target with the
#' `structural_status`, `fixed_value` and `reason` of the target, its
#' `map_type` (`"identity"`: the target is its own fitted source, including a
#' log-transformed source with an exp output; `"affine"`: a linear
#' combination of fitted sources; `"exp_affine"`: exp of a linear combination
#' of identity and log-transformed sources; `"state_dependent"`: a map requiring
#' fitted multiplier states; `"unsupported"`: any other map),
#' and the `support` of the map (a list column of `c(lower, upper)`:
#' `c(0, Inf)` for maps with an exp output and `c(-Inf, Inf)` otherwise; the
#' prior support of a target can be narrower).
#' `JAGS_formula_coefficient_transform_schema()` returns field descriptions.
#' `JAGS_formula_prior_density()` returns a `prior_linear_density` accepted by
#' [prior_density_ordinate()]. An unavailable density stops with an error of
#' class `BayesTools_formula_prior_density_unavailable` (also
#' `BayesTools_formula_transform_unavailable`) whose field `reason` names the
#' cause, e.g. `"unknown_target"`, `"missing_source_coordinates"`, or
#' `"nonlinear_map"` for a weighted combination with a target whose map is not
#' linear in the fitted coefficients. Declared measure limitations additionally
#' inherit `BayesTools_formula_measure_unavailable`; generic context/program
#' failures and numeric transformation failures retain their strict boundary.
#' Numeric refusal reasons include `"zero_multiplier_denominator"`,
#' `"missing_multiplier_state"`, and `"nonfinite_transform"`. Formula-design 6,
#' unscale-design 2, and transform 3 require refitting older stored formats;
#' readers never reconstruct missing fitted declaration ownership.
#' A requested parameter without a persisted formula design stops with
#' `BayesTools_formula_transform_unavailable` and reason
#' `"missing_formula_design"`, also when requested through either prior-density
#' route.
#' `JAGS_formula_internal_coordinate_priors()` returns a uniquely named list of
#' scalar [prior()] objects keyed by concrete fitted coordinate.
#'
#' @export JAGS_formula_coefficient_transform
#' @export JAGS_formula_coefficient_transform_schema
#' @export JAGS_formula_prior_density
#' @export JAGS_formula_internal_coordinate_priors
#' @name JAGS_formula_coefficient_transform
NULL

#' @rdname JAGS_formula_coefficient_transform
JAGS_formula_internal_coordinate_priors <- function(fit){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  JAGS_validate_fit_contract(
    fit,
    requires = c("formula_design", "parameter_map")
  )

  designs <- JAGS_formula_design(fit)
  coordinates <- parameter_coordinates(fit)
  out <- list()
  for(parameter in names(designs)){
    design <- designs[[parameter]]
    for(random_term in .bt_formula_design_random_effects(design)){
      correlation <- random_term$correlation
      if(!is.list(correlation) || !identical(correlation$type, "lkj")){
        next
      }

      K <- random_term$n_columns
      primitive_names <- .bt_random_effect_lkj_primitive_names(
        random_term,
        K,
        context = "Internal formula-coordinate prior metadata"
      )
      alpha <- .bt_lkj_cholesky_alpha(K = K, eta = correlation$eta)
      if(length(primitive_names) != length(alpha)){
        stop(
          "Stored LKJ primitive metadata do not match the compiled random-effect dimension.",
          call. = FALSE
        )
      }
      for(i in seq_along(primitive_names)){
        coordinate_name <- primitive_names[[i]]
        coordinate <- coordinates[
          coordinates$coordinate_name == coordinate_name,
          ,
          drop = FALSE
        ]
        valid <- nrow(coordinate) == 1L &&
          identical(coordinate$role, "random_correlation_coordinate") &&
          isTRUE(coordinate$internal) &&
          identical(coordinate$monitor_status, "sampled")
        if(!valid){
          .bt_stop_refit_required(
            "LKJ primitive coordinate '", coordinate_name,
            "' is missing from the fitted parameter map. Refit the model with the current BayesTools version."
          )
        }
        if(coordinate_name %in% names(out)){
          stop(
            "Formula metadata contain duplicate internal coordinate '",
            coordinate_name, "'.",
            call. = FALSE
          )
        }
        out[[coordinate_name]] <- prior(
          "beta",
          parameters = list(alpha = alpha[[i]], beta = alpha[[i]])
        )
      }
    }
  }

  out
}

#' @rdname JAGS_formula_coefficient_transform
JAGS_formula_coefficient_transform <- function(
    fit, parameter, target_scale = "original"){

  check_char(parameter, "parameter", check_length = 1L, allow_NA = FALSE)
  check_char(target_scale, "target_scale", check_length = 1L,
             allow_NA = FALSE)
  if(!identical(target_scale, "original")){
    .bt_formula_transform_stop(
      "Only target_scale = 'original' is supported by the current coefficient-transform schema.",
      parameter = parameter,
      reason = "unsupported_target_scale"
    )
  }
  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  JAGS_validate_fit_contract(
    fit,
    requires = c("formula_design", "parameter_map")
  )
  designs <- JAGS_formula_design(fit)
  design <- designs[[parameter, exact = TRUE]]
  if(is.null(design)){
    .bt_formula_transform_stop(
      paste0("Formula design for parameter '", parameter, "' is unavailable."),
      parameter = parameter,
      reason = "missing_formula_design"
    )
  }
  .bt_validate_formula_design_replay_schema(
    design,
    context = paste0("Coefficient transform for parameter '", parameter, "'")
  )
  .bt_formula_scale_finalized_check(design$formula_scale, parameter, require_owner = TRUE)
  if(!identical(attr(design$formula_scale, "unscale_design", exact = TRUE),
                attr(attr(fit, "formula_scale", exact = TRUE)[[parameter]], "unscale_design", exact = TRUE))){
    .bt_stop_refit_required("Formula declaration carriers disagree. Refit the model with this version of BayesTools.")
  }
  coordinates <- parameter_coordinates(fit)
  sources <- .bt_formula_coefficient_sources(design, coordinates, parameter)

  .bt_formula_coefficient_transform(
    source_names = sources$source,
    formula_scale = .bt_formula_scale_with_unscale_design(
      design$formula_scale,
      design
    ),
    log_intercept = design$log_intercept,
    parameter = parameter,
    target_scale = target_scale,
    source_metadata = sources[, c("source", "monitor_status", "fixed_value"),
                              drop = FALSE],
    formula_design_version = design$schema_version,
    parameter_map_version = parameter_map(fit)$schema_version
  )
}

#' @rdname JAGS_formula_coefficient_transform
JAGS_formula_coefficient_transform_schema <- function(){

  data.frame(
    field = c(
      "schema_version", "formula_design_version",
      "parameter_map_version", "parameter", "target_scale",
      "source_names", "target_names", "matrix", "basis_matrix", "multipliers",
      "state_constants", "state_dependencies", "prior_recipes", "source_transforms",
      "output_transforms", "dependencies", "sources", "targets"
    ),
    type = c(
      rep("integer", 3L), rep("character", 4L), "numeric matrix",
      "numeric matrix", "named list", "named numeric", "data.frame", "named list",
      rep("named character", 2L), rep("data.frame", 3L)
    ),
    description = c(
      "Coefficient-transform schema version.",
      "Formula-design schema version used to construct the transform.",
      "Parameter-map schema version used for source status.",
      "Formula output parameter.",
      "Requested target coefficient scale.",
      "Ordered fitted-coordinate source names.",
      "Ordered original-coordinate target names.",
      "Static target-by-source raw coefficient matrix; dynamic rows are entirely NA.",
      "Finite fitted-to-original design geometry on the same axes as matrix.",
      "Compiler-owned constant or exact backend-state multiplier declarations.",
      "Declared literal coefficient points and scalar model-data/ordinary point constants.",
      "True nonconstant coefficient/multiplier state dependencies of dynamic targets.",
      "Per-target raw_affine, contribution_affine or unavailable prior-law recipes.",
      "Named identity/log transform for every fitted source.",
      "Named identity/exp transform for every target.",
      "One row per exact nonzero target/source coefficient.",
      "Registry-linked source monitor status and fixed value.",
      paste0(
        "Target structural/dependent status, exact fixed value, map type ",
        "(identity, affine, exp_affine, or unsupported), and the support of ",
        "the map as a list column of c(lower, upper)."
      )
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname JAGS_formula_coefficient_transform
JAGS_formula_prior_density <- function(
    fit, parameter, target = NULL, target_scale = "original", context = NULL,
    weights = NULL){

  if(is.null(target) == is.null(weights)){
    stop("Supply exactly one of 'target' and 'weights'.", call. = FALSE)
  }
  if(!is.null(weights)){
    return(.bt_formula_prior_density_weights(
      fit          = fit,
      parameter    = parameter,
      weights      = weights,
      target_scale = target_scale,
      context      = context
    ))
  }
  check_char(target, "target", check_length = 1L, allow_NA = FALSE)
  transform <- JAGS_formula_coefficient_transform(
    fit,
    parameter = parameter,
    target_scale = target_scale
  )
  target_i <- match(target, transform$target_names)
  if(is.na(target_i)){
    .bt_formula_density_stop(
      paste0("Target coefficient '", target, "' is not available for formula parameter '",
             parameter, "'."),
      parameter = parameter,
      target = target,
      reason = "unknown_target",
      available = transform$target_names
    )
  }
  target_metadata <- transform$targets[target_i, , drop = FALSE]
  if(identical(target_metadata$structural_status, "unavailable")){
    .bt_formula_density_stop(
      paste0("Prior density for target coefficient '", target,
             "' is structurally unavailable."),
      parameter = parameter,
      target = target,
      reason = "structural_target_law_unavailable",
      detail = target_metadata$reason
    )
  }

  recipe <- .bt_formula_prior_recipe_weights(transform, stats::setNames(1, target))
  weights <- recipe$weights
  density_context <- .bt_formula_prior_density_context(
    fit, context = context, source_names = transform$source_names,
    parameter = parameter, target = target, transform = transform, recipe = recipe$type
  )
  if(identical(recipe$type, "contribution_affine")) .bt_formula_require_multiplier_laws(density_context, weights, parameter)
  missing <- setdiff(names(weights), density_context$column_names)
  if(length(missing) > 0L){
    .bt_formula_density_stop(
      paste0(
        "Prior-density context is missing fitted source coordinate",
        if(length(missing) > 1L) "s " else " ",
        paste0("'", missing, "'", collapse = ", "), "."
      ),
      parameter = parameter,
      target = target,
      reason = "missing_source_coordinates",
      missing = missing
    )
  }

  source_transforms <- transform$source_transforms[names(weights)]
  source_transforms[source_transforms == "identity"] <- NA_character_
  if(all(is.na(source_transforms))){
    source_transforms <- NULL
  }
  output_transform <- transform$output_transforms[[target]]
  if(identical(output_transform, "identity")){
    output_transform <- NULL
  }

  density <- tryCatch(
    .prior_density_from_context(
      context = density_context,
      weights = weights,
      source_transforms = source_transforms,
      output_transformation = output_transform
    ),
    error = function(e){
      if(inherits(e, "BayesTools_formula_measure_unavailable")) stop(e)
      .bt_formula_density_stop(
        paste0("Prior density for target coefficient '", target,
               "' is unavailable: ", conditionMessage(e)),
        parameter = parameter,
        target = target,
        reason = "density_context_error",
        parent = e
      )
    }
  )
  attr(density, "formula_coefficient_transform") <- transform
  attr(density, "formula_coefficient_target") <- target
  density
}

# The original-scale prior density of the weighted combination sum_t v_t T_t
# of targets whose maps are identity or affine in the fitted sources: the
# joint prior-density context at the combined source weights v %*% M (a
# log-transformed source with an exp output is its own target, so its weight
# stays on the positive source). Other maps are not linear in the sources and
# their combinations are unavailable.
.bt_formula_prior_density_weights <- function(fit, parameter, weights,
                                              target_scale, context){

  if(!is.numeric(weights) || length(weights) == 0L || is.null(names(weights)) ||
     anyNA(names(weights)) || any(!nzchar(names(weights))) ||
     anyDuplicated(names(weights)) || any(!is.finite(weights))){
    stop("'weights' must be a named numeric vector of finite target weights.",
         call. = FALSE)
  }
  transform <- JAGS_formula_coefficient_transform(
    fit,
    parameter = parameter,
    target_scale = target_scale
  )
  label <- paste0(
    paste0(format(unname(weights)), " * ", names(weights)),
    collapse = " + "
  )
  unknown <- setdiff(names(weights), transform$target_names)
  if(length(unknown) > 0L){
    .bt_formula_density_stop(
      paste0(
        "Target coefficient", if(length(unknown) > 1L) "s " else " ",
        paste0("'", unknown, "'", collapse = ", "),
        if(length(unknown) > 1L) " are" else " is",
        " not available for formula parameter '", parameter, "'."
      ),
      parameter = parameter,
      target = label,
      reason = "unknown_target",
      available = transform$target_names
    )
  }
  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    stop("'weights' must contain at least one nonzero target weight.",
         call. = FALSE)
  }
  targets <- transform$targets[
    match(names(weights), transform$target_names), , drop = FALSE
  ]
  unavailable <- targets$structural_status == "unavailable"
  if(any(unavailable)){
    .bt_formula_density_stop(
      paste0("Prior density for target coefficient '",
             targets$target[unavailable][[1L]],
             "' is structurally unavailable."),
      parameter = parameter,
      target = label,
      reason = "structural_target_law_unavailable",
      detail = targets$reason[unavailable][[1L]]
    )
  }
  nonlinear <- targets$map_type %in% c("exp_affine", "unsupported") |
    (targets$map_type == "state_dependent" &
       transform$output_transforms[targets$target] != "identity")
  if(any(nonlinear)) .bt_formula_density_stop(
    paste0("Prior density of a weighted combination of target coefficients is ",
      "unavailable: the map of '", targets$target[nonlinear][[1L]],
      "' from the fitted coefficients is ", targets$map_type[nonlinear][[1L]], ", not linear."),
    parameter = parameter, target = label, reason = "nonlinear_map")
  recipe <- .bt_formula_prior_recipe_weights(transform, weights)
  source_weights <- recipe$weights
  density_context <- .bt_formula_prior_density_context(
    fit, context = context, source_names = transform$source_names,
    parameter = parameter, target = label, transform = transform, recipe = recipe$type
  )
  if(identical(recipe$type, "contribution_affine")) .bt_formula_require_multiplier_laws(density_context, source_weights, parameter)
  missing <- setdiff(names(source_weights), density_context$column_names)
  if(length(missing) > 0L){
    .bt_formula_density_stop(
      paste0(
        "Prior-density context is missing fitted source coordinate",
        if(length(missing) > 1L) "s " else " ",
        paste0("'", missing, "'", collapse = ", "), "."
      ),
      parameter = parameter,
      target = label,
      reason = "missing_source_coordinates",
      missing = missing
    )
  }

  density <- tryCatch(
    .prior_density_from_context(
      context = density_context,
      weights = source_weights,
      output_transformation = if(isTRUE(recipe$offset != 0)) "lin" else NULL,
      output_transformation_arguments = if(isTRUE(recipe$offset != 0)) list(a = recipe$offset, b = 1) else NULL
    ),
    error = function(e){
      if(inherits(e, "BayesTools_formula_measure_unavailable")) stop(e)
      .bt_formula_density_stop(
        paste0("Prior density of the weighted target combination is ",
               "unavailable: ", conditionMessage(e)),
        parameter = parameter,
        target = label,
        reason = "density_context_error",
        parent = e
      )
    }
  )
  attr(density, "formula_coefficient_transform") <- transform
  attr(density, "formula_coefficient_weights") <- weights
  density
}

.bt_formula_coefficient_sources <- function(design, coordinates, parameter){

  prior_list <- design$prior_list
  if(!is.list(prior_list) || is.null(names(prior_list)) ||
     anyNA(names(prior_list)) || any(!nzchar(names(prior_list))) ||
     anyDuplicated(names(prior_list)) ||
     any(!vapply(prior_list, is.prior, logical(1)))){
    .bt_formula_transform_stop(
      paste0("Formula prior metadata for parameter '", parameter,
             "' are incomplete or unsupported."),
      parameter = parameter,
      reason = "unsupported_formula_priors"
    )
  }
  prior_coordinates <- unique(unlist(Map(
    .prior_linear_prior_columns,
    names(prior_list),
    prior_list
  ), use.names = FALSE))
  rows <- coordinates$formula_parameter == parameter &
    coordinates$role == "fixed_coefficient" & !coordinates$internal
  source_coordinates <- coordinates[rows, , drop = FALSE]
  actual <- source_coordinates$coordinate_name
  expected <- prior_coordinates[prior_coordinates %in% actual]
  if(length(expected) == 0L || !setequal(expected, actual)){
    .bt_formula_transform_stop(
      paste0("Formula coefficient sources for parameter '", parameter,
             "' disagree with the parameter-coordinate table."),
      parameter = parameter,
      reason = "coordinate_source_mismatch",
      expected = expected,
      registered = actual,
      prior_coordinates = prior_coordinates
    )
  }
  source_coordinates <- source_coordinates[match(expected, actual), , drop = FALSE]
  data.frame(
    source = source_coordinates$coordinate_name,
    monitor_status = source_coordinates$monitor_status,
    fixed_value = source_coordinates$fixed_value,
    stringsAsFactors = FALSE
  )
}

# The transform is a function of its arguments only, and the prior-density
# contexts, the scale transformations of posteriors and the atoms of one
# original-scale request rebuild the same transform again and again (one
# context per parameter); a request for arguments identical() to earlier ones
# returns the kept transform (see .bt_content_memo()).
.bt_formula_coefficient_transform <- function(
    source_names, formula_scale, parameter, target_scale = "original",
    log_intercept = FALSE,
    source_metadata = NULL,
    formula_design_version = .bt_formula_design_schema_version(),
    parameter_map_version = .bt_parameter_map_version){

  .bt_content_memo(
    "formula_coefficient_transform",
    list(source_names, formula_scale, parameter, target_scale, log_intercept,
         source_metadata, formula_design_version, parameter_map_version),
    function() .bt_formula_coefficient_transform_uncached(
      source_names, formula_scale, parameter, target_scale, log_intercept,
      source_metadata, formula_design_version, parameter_map_version
    )
  )
}

.bt_formula_coefficient_transform_uncached <- function(
    source_names, formula_scale, parameter, target_scale = "original",
    log_intercept = FALSE,
    source_metadata = NULL,
    formula_design_version = .bt_formula_design_schema_version(),
    parameter_map_version = .bt_parameter_map_version){

  check_char(source_names, "source_names", check_length = FALSE,
             allow_NA = FALSE)
  check_char(parameter, "parameter", check_length = 1L, allow_NA = FALSE)
  if(length(source_names) == 0L || anyDuplicated(source_names)){
    stop("'source_names' must contain unique fitted coordinates.",
         call. = FALSE)
  }
  if(!is.null(formula_scale) && length(formula_scale) > 0L){
    scale_list <- stats::setNames(list(formula_scale), parameter)
    .check_formula_scale_info(scale_list)
    matrix <- .build_unscale_matrix(
      source_names,
      formula_scale,
      prefix = parameter
    )
  }else{
    matrix <- diag(length(source_names))
    rownames(matrix) <- colnames(matrix) <- source_names
  }

  source_transforms <- stats::setNames(
    rep("identity", length(source_names)),
    source_names
  )
  output_transforms <- stats::setNames(
    rep("identity", length(source_names)),
    source_names
  )
  log_intercept <- isTRUE(log_intercept) ||
    (!is.null(formula_scale) &&
       length(formula_scale) > 0L &&
       isTRUE(attr(formula_scale, "log_intercept")))
  intercept <- paste0(parameter, "_intercept")
  if(log_intercept && intercept %in% source_names){
    source_transforms[[intercept]] <- "log"
    output_transforms[[intercept]] <- "exp"
  }

  if(is.null(source_metadata)){
    source_metadata <- data.frame(
      source = source_names,
      monitor_status = rep("unknown", length(source_names)),
      fixed_value = rep(NA_real_, length(source_names)),
      stringsAsFactors = FALSE
    )
  }
  source_metadata <- source_metadata[match(source_names, source_metadata$source),
                                     , drop = FALSE]
  spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  if(!is.null(spec)) .bt_formula_unscale_design_spec_check(spec, parameter)
  sources <- data.frame(
    source = source_names,
    monitor_status = source_metadata$monitor_status,
    fixed_value = as.numeric(source_metadata$fixed_value),
    source_transform = unname(source_transforms),
    stringsAsFactors = FALSE
  )
  basis_matrix <- matrix
  spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  multipliers <- stats::setNames(lapply(source_names, function(source){
    declaration <- spec$multipliers[[source]]
    if(is.null(declaration)) list(type = "constant", value = 1) else declaration
  }), source_names)
  state_constants <- .bt_formula_state_constants(formula_scale)
  declared <- match(sources$source, names(state_constants))
  fixed <- !is.na(declared)
  sources$monitor_status[fixed] <- "structural"
  sources$fixed_value[fixed] <- state_constants[declared[fixed]]
  analysis <- .bt_formula_coefficient_analysis(basis_matrix, multipliers,
    state_constants, source_transforms, sources, output_transforms)
  matrix <- analysis$matrix
  dependencies <- .bt_formula_coefficient_dependencies(matrix)
  targets <- analysis$targets

  out <- list(
    schema_version = .bt_formula_coefficient_transform_version,
    formula_design_version = as.integer(formula_design_version),
    parameter_map_version = as.integer(parameter_map_version),
    parameter = parameter,
    target_scale = target_scale,
    source_names = source_names,
    target_names = rownames(matrix),
    matrix = matrix,
    basis_matrix = basis_matrix,
    multipliers = multipliers,
    state_constants = state_constants,
    state_dependencies = analysis$state_dependencies,
    prior_recipes = analysis$prior_recipes,
    source_transforms = source_transforms,
    output_transforms = output_transforms,
    dependencies = dependencies,
    sources = sources,
    targets = targets
  )
  class(out) <- c("BayesTools_formula_coefficient_transform", "list")
  .bt_validate_formula_coefficient_transform(out)
  out
}

.bt_formula_coefficient_dependencies <- function(matrix){

  rows <- list()
  row_i <- 0L
  for(target_i in seq_len(nrow(matrix))){
    source_i <- which(matrix[target_i, ] != 0)
    for(i in source_i){
      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        target = rownames(matrix)[target_i],
        source = colnames(matrix)[i],
        coefficient = unname(matrix[target_i, i]),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.bt_formula_multiplier_constant <- function(declaration, constants){

  if(identical(declaration$type, "constant")) return(declaration$value)
  if(declaration$name %in% names(constants)) return(unname(constants[[declaration$name]]))
  NULL
}

# Retain the original grouping when its intermediate products are normal.
# A ratio-first grouping can recover a representable final value when that
# grouping loses range. Only final results may be subnormal; neither order
# certifies zero from nonzero operands.
.bt_formula_coefficient_ratio <- function(basis, numerator, denominator,
                                           source, target, value = NULL){

  normal <- function(x) is.finite(x) & abs(x) >= .Machine$double.xmin
  final <- function(x) is.finite(x) & x != 0
  if(is.null(value)){
    product <- basis * numerator
    direct <- product / denominator
    valid <- normal(product) & final(direct)
    ratio <- numerator / denominator
    alternate <- ratio * basis
    alternate_valid <- normal(ratio) & final(alternate)
    zero <- numerator == 0
  }else{
    first <- basis * value
    product <- first * numerator
    direct <- product / denominator
    valid <- normal(first) & normal(product) & final(direct)
    ratio <- numerator / denominator
    second <- ratio * value
    alternate <- second * basis
    alternate_valid <- normal(ratio) & normal(second) & final(alternate)
    zero <- value == 0 | numerator == 0
  }
  resolved <- valid | alternate_valid | zero
  if(any(!resolved)) .bt_formula_transform_stop(
    "Original coefficient ratio arithmetic is numerically unavailable.",
    source = source, target = target, reason = "nonfinite_transform",
    observed = list(indices = which(!resolved), basis = basis, value = value,
      numerator = numerator, denominator = denominator, product = product,
      direct = direct, ratio = ratio, alternate = alternate))
  out <- direct
  out[!valid] <- rep_len(alternate, length(out))[!valid]
  out[zero] <- 0
  out
}

.bt_formula_coefficient_term <- function(source, target, basis, multipliers,
                                         constants, source_transforms){

  numerator <- multipliers[[source]]
  denominator <- multipliers[[target]]
  cancel <- identical(numerator, denominator)
  constant_numerator <- .bt_formula_multiplier_constant(numerator, constants)
  constant_denominator <- .bt_formula_multiplier_constant(denominator, constants)
  zero_source <- source %in% names(constants) && identical(source_transforms[[source]], "identity") &&
    isTRUE(constants[[source]] == 0)
  zero <- zero_source || (!cancel && !is.null(constant_numerator) && constant_numerator == 0)
  static <- zero || cancel || (!is.null(constant_numerator) &&
                                !is.null(constant_denominator) && constant_denominator != 0)
  weight <- if(zero) 0 else if(cancel) basis else if(static)
    .bt_formula_coefficient_ratio(basis, constant_numerator, constant_denominator, source, target) else NA_real_
  list(zero = zero, cancel = cancel, static = static, weight = weight, basis = unname(basis))
}

.bt_formula_coefficient_analysis <- function(basis, multipliers, constants,
                                             source_transforms, sources,
                                             output_transforms){

  matrix <- basis * 0
  dependencies <- list()
  recipes <- list()
  dynamic <- stats::setNames(rep(FALSE, nrow(basis)), rownames(basis))
  for(target in rownames(basis)){
    for(source in colnames(basis)[basis[target, ] != 0]){
      term <- .bt_formula_coefficient_term(source, target, basis[target, source],
        multipliers, constants, source_transforms)
      if(term$zero) next
      if(!source %in% names(constants)) dependencies[[length(dependencies) + 1L]] <-
        data.frame(target = target, source = source, role = "coefficient", stringsAsFactors = FALSE)
      if(!term$cancel){
        for(declaration in multipliers[c(source, target)]){
          if(is.null(.bt_formula_multiplier_constant(declaration, constants))){
            dependencies[[length(dependencies) + 1L]] <- data.frame(
              target = target, source = declaration$name, role = "multiplier", stringsAsFactors = FALSE)
          }
        }
      }
      matrix[target, source] <- term$weight
      dynamic[[target]] <- dynamic[[target]] || !term$static
    }
    if(dynamic[[target]]) matrix[target, ] <- NA_real_
    denominator <- .bt_formula_multiplier_constant(multipliers[[target]], constants)
    recipes[[target]] <- if(!dynamic[[target]]){
      list(type = "raw_affine", weights = stats::setNames(as.numeric(matrix[target, , drop = FALSE]), colnames(matrix)), reason = NA_character_)
    }else if(!is.null(denominator) && denominator != 0){
      weights <- stats::setNames(numeric(ncol(basis)), colnames(basis))
      numerical_scale <- FALSE
      for(source in colnames(basis)[basis[target, ] != 0]){
        if(.bt_formula_coefficient_term(source, target, basis[target, source], multipliers,
                                        constants, source_transforms)$zero) next
        weight <- basis[target, source] / denominator
        if(!is.finite(weight) || weight == 0) numerical_scale <- TRUE else weights[[source]] <- weight
      }
      if(numerical_scale) list(type = "unavailable", weights = stats::setNames(numeric(), character()),
        reason = "numerical_scale_unavailable") else
          list(type = "contribution_affine", weights = weights, reason = NA_character_)
    }else{
      list(type = "unavailable", weights = stats::setNames(numeric(), character()),
        reason = "state_dependent_map")
    }
  }
  finite_matrix <- matrix
  finite_matrix[is.na(finite_matrix)] <- 0
  targets <- .bt_formula_coefficient_targets(finite_matrix, sources, output_transforms)
  targets$map_type[dynamic] <- "state_dependent"
  targets$structural_status[dynamic] <- "dependent"
  targets$fixed_value[dynamic] <- NA_real_
  targets$reason[dynamic] <- NA_character_
  state_dependencies <- if(length(dependencies)) unique(do.call(rbind, dependencies)) else
    data.frame(target = character(), source = character(), role = character(), stringsAsFactors = FALSE)
  state_dependencies <- state_dependencies[state_dependencies$target %in% names(dynamic)[dynamic], , drop = FALSE]
  rownames(state_dependencies) <- NULL
  list(matrix = matrix, state_dependencies = state_dependencies, prior_recipes = recipes, targets = targets)
}

# Static consumers must select requested nonzero rows before accessing C.
.bt_formula_static_rows <- function(transform, targets){

  .bt_validate_formula_coefficient_transform(transform)
  if(any(!targets %in% transform$target_names)) stop("Unknown coefficient target.", call. = FALSE)
  rows <- transform$matrix[targets, , drop = FALSE]
  if(any(!is.finite(rows))) .bt_formula_density_stop(
    "A static coefficient projection is unavailable for a state-dependent map.",
    parameter = transform$parameter, target = targets, reason = "state_dependent_map")
  rows
}

.bt_formula_coefficient_targets <- function(matrix, sources,
                                             output_transforms){

  rows <- vector("list", nrow(matrix))
  for(target_i in seq_len(nrow(matrix))){
    dependencies <- which(matrix[target_i, ] != 0)
    statuses <- sources$monitor_status[dependencies]
    status <- if(any(statuses == "unavailable")){
      "unavailable"
    }else if(all(statuses == "structural")){
      "structural"
    }else{
      "dependent"
    }
    fixed_value <- NA_real_
    reason <- NA_character_
    if(identical(status, "structural")){
      fixed_result <- tryCatch(
        list(
          value = .bt_formula_coefficient_fixed_value(
            matrix[target_i, dependencies],
            sources[dependencies, , drop = FALSE],
            output_transforms[[target_i]]
          ),
          reason = NA_character_
        ),
        error = function(e) list(
          value = NA_real_,
          reason = conditionMessage(e)
        )
      )
      fixed_value <- fixed_result$value
      reason <- fixed_result$reason
      if(!is.finite(fixed_value)){
        status <- "unavailable"
      }
    }else if(identical(status, "unavailable")){
      reason <- "At least one nonzero fitted source is unavailable."
    }
    rows[[target_i]] <- data.frame(
      target = rownames(matrix)[target_i],
      structural_status = status,
      fixed_value = fixed_value,
      reason = reason,
      map_type = .bt_formula_coefficient_map_type(
        weights = stats::setNames(
          matrix[target_i, dependencies],
          colnames(matrix)[dependencies]
        ),
        target = rownames(matrix)[target_i],
        source_transforms = sources$source_transform[dependencies],
        output_transform = output_transforms[[target_i]]
      ),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  # the image of each target's map: exp outputs are positive
  out$support <- unname(lapply(output_transforms, function(output_transform){
    if(identical(output_transform, "exp")) c(0, Inf) else c(-Inf, Inf)
  }))
  out
}

# The type of the map from the fitted sources to one original-scale target:
# "identity" (the target is its own source: an identity source and output,
# or a log source with an exp output), "affine" (a linear combination of
# identity sources), "exp_affine" (exp of a linear combination of identity
# and log sources), or "unsupported" (a log source without an exp output).
.bt_formula_coefficient_map_type <- function(weights, target,
                                             source_transforms,
                                             output_transform){

  sources <- names(weights)
  unit <- length(weights) == 1L && identical(sources, target) &&
    isTRUE(unname(weights) == 1)
  if(unit && ((identical(unname(source_transforms), "identity") &&
               identical(output_transform, "identity")) ||
              (identical(unname(source_transforms), "log") &&
               identical(output_transform, "exp")))){
    return("identity")
  }
  if(identical(output_transform, "identity") &&
     all(source_transforms == "identity")){
    return("affine")
  }
  if(identical(output_transform, "exp") &&
     all(source_transforms %in% c("identity", "log"))){
    return("exp_affine")
  }

  "unsupported"
}

.bt_formula_coefficient_fixed_value <- function(weights, sources,
                                                output_transform){

  values <- sources$fixed_value
  log_sources <- sources$source_transform == "log"
  if(any(log_sources)){
    if(any(values[log_sources] <= 0)){
      stop("A structurally fixed log-transformed source is not positive.",
           call. = FALSE)
    }
    values[log_sources] <- log(values[log_sources])
  }
  value <- sum(unname(weights) * values)
  if(identical(output_transform, "exp")){
    value <- exp(value)
  }
  if(!is.finite(value)){
    stop("The structurally fixed target value is not finite.", call. = FALSE)
  }
  value
}

.bt_apply_formula_coefficient_transform <- function(samples, transform, targets = transform$target_names){

  .bt_validate_formula_coefficient_transform(transform)
  samples <- as.matrix(samples)
  if(!is.numeric(samples) || is.null(colnames(samples)) || anyDuplicated(colnames(samples))){
    stop("Fitted coefficient state must be a numeric matrix with unique columns.", call. = FALSE)
  }
  if(any(!targets %in% transform$target_names)) stop("Unknown requested coefficient target.", call. = FALSE)
  needed <- character()
  for(target_name in targets){
    for(source_name in colnames(transform$basis_matrix)[transform$basis_matrix[target_name, ] != 0]){
      term <- .bt_formula_coefficient_term(source_name, target_name,
        transform$basis_matrix[target_name, source_name], transform$multipliers,
        transform$state_constants, transform$source_transforms)
      if(term$zero) next
      needed <- union(needed, source_name)
      if(!term$cancel) for(declaration in transform$multipliers[c(source_name, target_name)]){
        if(identical(declaration$type, "state")) needed <- union(needed, declaration$name)
      }
    }
  }
  state <- .bt_formula_materialize_state(samples, transform$state_constants, validate = needed)
  target <- matrix(0, nrow(state), length(targets), dimnames = list(rownames(state), targets))
  read_state <- function(name, role){
    if(!name %in% colnames(state)) .bt_formula_transform_stop(
      paste0("Required fitted state '", name, "' is unavailable."),
      parameter = transform$parameter, reason = "missing_multiplier_state", state = name, role = role)
    value <- state[, name]
    if(any(!is.finite(value))) .bt_formula_transform_stop(
      paste0("Required fitted state '", name, "' is nonfinite."),
      parameter = transform$parameter, reason = "nonfinite_transform", state = name, role = role)
    value
  }
  multiplier <- function(declaration){
    value <- .bt_formula_multiplier_constant(declaration, transform$state_constants)
    if(!is.null(value)) return(value)
    read_state(declaration$name, "multiplier")
  }
  for(name in targets){
    static_weights <- stats::setNames(as.numeric(transform$matrix[name, , drop = FALSE]), transform$source_names)
    if(all(is.finite(static_weights))){
      active <- names(static_weights)[static_weights != 0]
      source <- matrix(0, nrow(state), length(active), dimnames = list(NULL, active))
      for(coordinate in active){
        value <- read_state(coordinate, "coefficient")
        if(identical(transform$source_transforms[[coordinate]], "log")){
          if(any(value <= 0)) .bt_formula_transform_stop(
            "Log-transformed fitted coefficient samples must be positive.",
            parameter = transform$parameter, reason = "nonfinite_transform", state = coordinate)
          value <- log(value)
        }
        source[, coordinate] <- value
      }
      if(length(active)) target[, name] <- source %*% static_weights[active]
    }else{
    for(source in colnames(transform$basis_matrix)[transform$basis_matrix[name, ] != 0]){
      term <- .bt_formula_coefficient_term(source, name, transform$basis_matrix[name, source],
        transform$multipliers, transform$state_constants, transform$source_transforms)
      if(isTRUE(term$zero)) next
      value <- read_state(source, "coefficient")
      if(identical(transform$source_transforms[[source]], "log")){
        if(any(value <= 0)) .bt_formula_transform_stop(
          "Log-transformed fitted coefficient samples must be positive.",
          parameter = transform$parameter, reason = "nonfinite_transform", state = source)
        value <- log(value)
      }
      if(isTRUE(term$cancel)){
        contribution <- term$basis * value
        if(any(value != 0 & contribution == 0)) .bt_formula_transform_stop(
          "Original coefficient arithmetic is numerically unavailable.",
          parameter = transform$parameter, source = source, target = name,
          reason = "nonfinite_transform", observed = list(basis = term$basis, value = value,
            contribution = contribution))
      }else{
        numerator <- multiplier(transform$multipliers[[source]])
        denominator <- multiplier(transform$multipliers[[name]])
        if(any(denominator == 0)) .bt_formula_transform_stop(
          paste0("Original coefficient '", name, "' has a zero multiplier denominator."),
          parameter = transform$parameter, target = name, reason = "zero_multiplier_denominator")
        contribution <- .bt_formula_coefficient_ratio(term$basis, numerator, denominator,
          source, name, value = value)
      }
      if(any(!is.finite(contribution))) .bt_formula_transform_stop(
        "Original coefficient arithmetic is nonfinite.", parameter = transform$parameter,
        target = name, reason = "nonfinite_transform")
      target[, name] <- target[, name] + contribution
    }
    }
    if(identical(transform$output_transforms[[name]], "exp")){
      target[, name] <- exp(target[, name])
      if(any(target[, name] <= 0)) .bt_formula_transform_stop(
        "Original logged-intercept output is unrepresentable.", parameter = transform$parameter,
        target = name, reason = "nonfinite_transform")
    }
    if(any(!is.finite(target[, name]))) .bt_formula_transform_stop(
      "Original coefficient arithmetic is nonfinite.", parameter = transform$parameter,
      target = name, reason = "nonfinite_transform")
  }
  absent <- setdiff(targets, colnames(samples))
  if(length(absent)) samples <- cbind(samples, target[, absent, drop = FALSE])
  samples[, targets] <- target
  samples
}

.bt_validate_formula_coefficient_transform <- function(transform){

  expected_names <- c(
    "schema_version", "formula_design_version",
    "parameter_map_version", "parameter", "target_scale",
    "source_names", "target_names", "matrix", "basis_matrix", "multipliers",
    "state_constants", "state_dependencies", "prior_recipes", "source_transforms",
    "output_transforms", "dependencies", "sources", "targets"
  )
  valid <- inherits(transform, "BayesTools_formula_coefficient_transform") &&
    is.list(transform) && identical(names(transform), expected_names) &&
    identical(transform$schema_version,
              .bt_formula_coefficient_transform_version) &&
    is.integer(transform$formula_design_version) &&
    length(transform$formula_design_version) == 1L &&
    !is.na(transform$formula_design_version) &&
    identical(
      transform$formula_design_version,
      .bt_formula_design_schema_version()
    ) &&
    is.integer(transform$parameter_map_version) &&
    length(transform$parameter_map_version) == 1L &&
    !is.na(transform$parameter_map_version) &&
    identical(transform$parameter_map_version, .bt_parameter_map_version) &&
    is.character(transform$parameter) && length(transform$parameter) == 1L &&
    !is.na(transform$parameter) && nzchar(transform$parameter) &&
    identical(transform$target_scale, "original") &&
    is.character(transform$source_names) &&
    length(transform$source_names) > 0L &&
    !anyNA(transform$source_names) && !anyDuplicated(transform$source_names) &&
    identical(transform$target_names, transform$source_names) &&
    is.matrix(transform$matrix) && is.numeric(transform$matrix) &&
    all(apply(transform$matrix, 1L, function(row) all(is.finite(row)) || all(is.na(row)))) &&
    is.matrix(transform$basis_matrix) && is.numeric(transform$basis_matrix) &&
    all(is.finite(transform$basis_matrix)) &&
    identical(dimnames(transform$basis_matrix), dimnames(transform$matrix)) &&
    is.list(transform$multipliers) && identical(names(transform$multipliers), transform$source_names) &&
    all(vapply(transform$multipliers, .bt_formula_multiplier_valid, logical(1))) &&
    is.numeric(transform$state_constants) && all(is.finite(transform$state_constants)) &&
    (length(transform$state_constants) == 0L || (!is.null(names(transform$state_constants)) &&
      !anyNA(names(transform$state_constants)) && all(nzchar(names(transform$state_constants))))) &&
    !anyDuplicated(names(transform$state_constants)) &&
    identical(dimnames(transform$matrix),
              list(transform$target_names, transform$source_names)) &&
    is.character(transform$source_transforms) &&
    length(transform$source_transforms) == length(transform$source_names) &&
    !anyNA(transform$source_transforms) &&
    identical(names(transform$source_transforms), transform$source_names) &&
    all(transform$source_transforms %in% c("identity", "log")) &&
    is.character(transform$output_transforms) &&
    length(transform$output_transforms) == length(transform$target_names) &&
    !anyNA(transform$output_transforms) &&
    identical(names(transform$output_transforms), transform$target_names) &&
    all(transform$output_transforms %in% c("identity", "exp"))
  if(!valid){
    .bt_stop_refit_required("Formula coefficient transform metadata are missing or unsupported. Refit the model with this version of BayesTools.")
  }

  valid_sources <- is.data.frame(transform$sources) && identical(
    names(transform$sources),
    c("source", "monitor_status", "fixed_value", "source_transform")
  ) && identical(transform$sources$source, transform$source_names) &&
    all(transform$sources$monitor_status %in%
          c("sampled", "structural", "unavailable", "unknown")) &&
    is.numeric(transform$sources$fixed_value) &&
    all(is.na(transform$sources$fixed_value[
      transform$sources$monitor_status != "structural"
    ])) && all(is.finite(transform$sources$fixed_value[
      transform$sources$monitor_status == "structural"
    ])) && identical(transform$sources$source_transform,
                     unname(transform$source_transforms))
  expected_dependencies <- .bt_formula_coefficient_dependencies(
    transform$matrix
  )
  analysis <- .bt_formula_coefficient_analysis(transform$basis_matrix,
    transform$multipliers, transform$state_constants, transform$source_transforms,
    transform$sources, transform$output_transforms)
  expected_targets <- analysis$targets
  if(!valid_sources || !identical(transform$matrix, analysis$matrix) ||
     !identical(transform$state_dependencies, analysis$state_dependencies) ||
     !identical(transform$prior_recipes, analysis$prior_recipes) ||
     !identical(transform$dependencies,
                                  expected_dependencies) ||
     !identical(transform$targets, expected_targets)){
    .bt_stop_refit_required("Formula coefficient transform structural metadata are malformed. Refit the model with this version of BayesTools.")
  }
  invisible(TRUE)
}

.bt_formula_prior_density_context <- function(
    fit, context, source_names, parameter, target,
    transform = NULL, recipe = "raw_affine"){

  if(is.null(context)){
    prior_list <- attr(fit, "prior_list", exact = TRUE)
    if(is.null(prior_list)){
      .bt_formula_density_stop(
        "The fitted object has no prior metadata for coefficient-density construction.",
        parameter = parameter,
        target = target,
        reason = "missing_prior_list"
      )
    }
    context <- tryCatch(
      .prior_density_build_context(
        prior_list = prior_list,
        column_names = source_names
      ),
      error = function(e){
        .bt_formula_density_stop(
          paste0("Could not construct the fitted prior-density context: ",
                 conditionMessage(e)),
          parameter = parameter,
          target = target,
          reason = "unsupported_fit_prior_context",
          parent = e
        )
      }
    )
  }
  if(!.bt_formula_prior_density_context_valid(context)){
    .bt_formula_density_stop(
      "'context' must be a valid BayesTools prior-density context.",
      parameter = parameter,
      target = target,
      reason = "invalid_prior_density_context"
    )
  }

  # The public transform already maps original targets to fitted sources.
  # Remove any older context-owned unscaling map to avoid applying it twice.
  out <- context
  if(inherits(out, "prior_density_context")){
    out$formula_scale <- NULL
    out$transforms <- list()
  }else if(inherits(out, "prior_density_conditional_context")){
    out$formula_scale <- NULL
  }

  # The sources are the raw monitored coefficients (the JAGS nodes). A
  # coefficient's 'multiply_by' scales only its linear-predictor contribution,
  # so it must not enter the density of the coefficient itself.
  out <- .bt_formula_prior_density_context_raw_sources(out, source_names)
  if(identical(recipe, "contribution_affine")){
    out <- .bt_formula_prior_density_context_contributions(out, transform)
  }
  out
}

.bt_formula_prior_recipe_weights <- function(transform, weights){

  weights <- weights[weights != 0]
  recipes <- transform$prior_recipes[names(weights)]
  types <- vapply(recipes, `[[`, character(1), "type")
  if(any(types == "unavailable")) .bt_formula_density_stop(
    paste0("The requested original coefficient prior law is unavailable: ",
      recipes[[which(types == "unavailable")[[1L]]]]$reason, "."),
    parameter = transform$parameter, target = names(weights),
    reason = recipes[[which(types == "unavailable")[[1L]]]]$reason)
  # A static identity sibling can be expressed in the contribution space if
  # its own multiplier is an authoritative nonzero constant.
  type <- if(all(types == "raw_affine")) "raw_affine" else "contribution_affine"
  offset <- 0
  if(length(unique(types)) > 1L){
    for(target in names(recipes)[types == "raw_affine"]){
      denominator <- .bt_formula_multiplier_constant(transform$multipliers[[target]], transform$state_constants)
      if(is.null(denominator) || denominator == 0){
        target_metadata <- transform$targets[match(target, transform$targets$target), , drop = FALSE]
        if(identical(target_metadata$structural_status, "structural") && is.finite(target_metadata$fixed_value)){
          offset <- offset + weights[[target]] * target_metadata$fixed_value
          recipes[[target]]$weights <- stats::setNames(rep(0, length(transform$source_names)), transform$source_names)
          next
        }
        .bt_formula_density_stop("The requested coefficient laws have incompatible prior recipes.",
          parameter = transform$parameter, target = names(weights), reason = "incompatible_prior_recipes")
      }
      converted <- stats::setNames(numeric(length(transform$source_names)), transform$source_names)
      for(source in transform$source_names[transform$basis_matrix[target, ] != 0]){
        if(.bt_formula_coefficient_term(source, target, transform$basis_matrix[target, source],
          transform$multipliers, transform$state_constants, transform$source_transforms)$zero) next
        converted[[source]] <- transform$basis_matrix[target, source] / denominator
        if(!is.finite(converted[[source]]) || converted[[source]] == 0) .bt_formula_density_stop(
          "The requested coefficient prior law is unavailable because its contribution scale lost representable range.",
          parameter = transform$parameter, target = names(weights), source = source,
          reason = "numerical_scale_unavailable")
      }
      recipes[[target]]$weights <- converted
    }
  }
  out <- stats::setNames(rep(0, length(transform$source_names)), transform$source_names)
  for(target in names(weights)) out <- out + weights[[target]] * recipes[[target]]$weights
  list(type = type, weights = out[out != 0], offset = offset)
}

.bt_formula_prior_density_context_contributions <- function(context, transform){

  canonical <- function(priors){
    for(owner in names(priors)){
      prior <- priors[[owner]]
      if(!is.prior(prior)) next
      prior <- .bt_prior_without_multiply_by(prior)
      columns <- intersect(.prior_linear_prior_columns(owner, prior), transform$source_names)
      if(length(columns)){
        declaration <- transform$multipliers[[columns[[1L]]]]
        if(!all(vapply(transform$multipliers[columns], identical, logical(1), declaration))){
          .bt_formula_density_stop("One prior owner has incompatible compiled multipliers.",
            parameter = transform$parameter, reason = "incompatible_prior_recipes")
        }
        value <- .bt_formula_multiplier_constant(declaration, transform$state_constants)
        attr(prior, "multiply_by") <- if(is.prior.point(prior) &&
          is.numeric(prior$parameters$location) && isTRUE(prior$parameters$location == 0)) NULL else
          if(is.null(value)) declaration$name else if(value == 1) NULL else value
      }
      priors[[owner]] <- prior
    }
    priors
  }
  if(inherits(context, "prior_density_model_mixture_context")){
    for(model in seq_along(context$model_weights)){
      leaf <- canonical(.prior_density_model_prior_list(context$prior_list, model))
      for(owner in names(leaf)) context$prior_list[[owner]][[model]] <- leaf[[owner]]
    }
  }else context$prior_list <- canonical(context$prior_list)
  if(is.list(context$prior_lists)) context$prior_lists <- lapply(context$prior_lists, canonical)
  context$linear_weight_space <- "formula_contribution"
  context
}

.bt_formula_prior_density_context_raw_sources <- function(context,
                                                           source_names){

  source_priors <- unique(sub("\\[[^]]+\\]$", "", source_names))
  strip_list <- function(prior_list){
    if(!is.list(prior_list) || is.null(names(prior_list))){
      return(prior_list)
    }
    for(name in intersect(names(prior_list), source_priors)){
      prior <- prior_list[[name]]
      if(is.prior(prior)){
        prior_list[[name]] <- .bt_prior_without_multiply_by(prior)
      }else if(is.list(prior)){
        # model-mixture contexts hold one prior per model
        for(model_i in seq_along(prior)){
          if(is.prior(prior[[model_i]])){
            prior[[model_i]] <- .bt_prior_without_multiply_by(prior[[model_i]])
          }
        }
        prior_list[[name]] <- prior
      }
    }
    prior_list
  }

  context$prior_list <- strip_list(context$prior_list)
  if(inherits(context, "prior_density_conditional_context") &&
     is.list(context$prior_lists)){
    context$prior_lists <- lapply(context$prior_lists, strip_list)
  }

  context
}

.bt_prior_without_multiply_by <- function(prior){

  attr(prior, "multiply_by") <- NULL
  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    for(component_i in seq_along(prior)){
      if(is.prior(prior[[component_i]])){
        prior[[component_i]] <- .bt_prior_without_multiply_by(prior[[component_i]])
      }
    }
  }

  prior
}

.bt_formula_prior_density_context_valid <- function(context){

  known_class <- inherits(context, "prior_density_context") ||
    inherits(context, "prior_density_model_mixture_context") ||
    inherits(context, "prior_density_conditional_context")
  known_class && is.list(context) &&
    is.character(context$column_names) && length(context$column_names) > 0L &&
    !anyNA(context$column_names) && !anyDuplicated(context$column_names) &&
    is.numeric(context$n_grid) && length(context$n_grid) == 1L &&
    is.finite(context$n_grid) && context$n_grid >= 16 &&
    is.numeric(context$tail_prob) && length(context$tail_prob) == 1L &&
    is.finite(context$tail_prob) && context$tail_prob > 0 &&
    context$tail_prob < 0.5
}

.bt_formula_transform_stop <- function(message, ...){

  condition <- structure(
    c(list(message = message, call = NULL), list(...)),
    class = c("BayesTools_formula_transform_unavailable", "error",
              "condition")
  )
  stop(condition)
}

.bt_formula_raw_children <- function(values, priors){

  out <- list()
  for(owner in names(priors)){
    prior <- priors[[owner]]
    out[[owner]] <- if(is.prior.spike_and_slab(prior)){
      .as_mixed_posteriors.spike_and_slab(values, prior, owner)
    }else if(is.prior.mixture(prior)){
      .as_mixed_posteriors.mixture(values, prior, owner, character())
    }else if(is.prior.factor(prior)){
      .as_mixed_posteriors.factor(values, prior, owner)
    }else if(is.prior.vector(prior)){
      .as_mixed_posteriors.vector(values, prior, owner)
    }else if(is.prior.simple(prior)){
      .as_mixed_posteriors.simple(values, prior, owner)
    }else NULL
  }
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors")
  out
}

.bt_formula_state_new <- function(fit, values, owners, draw_index = seq_len(nrow(values)),
                                   condition_event = NULL){

  designs <- attr(fit, "formula_design", exact = TRUE)
  priors <- attr(fit, "prior_list", exact = TRUE)
  selected <- Filter(function(design){
    any(owners %in% intersect(names(design$prior_list), paste0(design$parameter, "_", design$model_terms)))
  }, designs)
  if(!length(selected)) return(NULL)
  scales <- lapply(selected, `[[`, "formula_scale")
  required <- retained_owners <- character()
  for(prefix in names(selected)){
    scale <- scales[[prefix]]
    .bt_formula_scale_finalized_check(scale, prefix, require_owner = TRUE)
    spec <- attr(scale, "unscale_design", exact = TRUE)
    columns <- names(spec$multipliers)
    selected_columns <- unique(unlist(lapply(intersect(owners, names(selected[[prefix]]$prior_list)), function(owner){
      .prior_linear_prior_columns(owner, priors[[owner]])
    }), use.names = FALSE))
    selected_columns <- intersect(selected_columns, columns)
    transform <- .bt_formula_coefficient_transform(columns, scale, prefix,
      log_intercept = selected[[prefix]]$log_intercept)
    numeric_dependencies <- transform$dependencies$source[
      transform$dependencies$target %in% selected_columns]
    numeric_dependencies <- union(numeric_dependencies, transform$state_dependencies$source[
      transform$state_dependencies$target %in% selected_columns])
    required <- union(required, numeric_dependencies)
    active <- columns[colSums(abs(transform$basis_matrix[intersect(selected_columns, columns), , drop = FALSE])) != 0]
    active <- union(active, selected_columns)
    constants <- .bt_formula_state_constants(scale)
    values <- .bt_formula_materialize_state(values, constants)
    for(column in active){
      zero <- column %in% names(constants) && constants[[column]] == 0 &&
        transform$source_transforms[[column]] == "identity"
      numerator <- .bt_formula_multiplier_constant(spec$multipliers[[column]], constants)
      if(zero || (!is.null(numerator) && numerator == 0)) next
      required <- union(required, column)
      declaration <- spec$multipliers[[column]]
      if(identical(declaration$type, "state")) required <- union(required, declaration$name)
    }
    retained_owners <- union(retained_owners, names(priors)[vapply(names(priors), function(owner){
      any(active %in% .prior_linear_prior_columns(owner, priors[[owner]])) || owner %in% required
    }, logical(1))])
    for(owner in names(priors)[vapply(names(priors), function(owner){
      is.prior.ordered(priors[[owner]]) && any(required %in% .prior_linear_prior_columns(owner, priors[[owner]]))
    }, logical(1))]){
      ordered_spec <- .bt_ordered_spec(owner, priors[[owner]])
      replay <- .bt_deterministic_node_evaluate(.bt_dnode_ordered_coefficients(ordered_spec),
        .bt_deterministic_lookup(values, priors))
      if(is.null(replay)) .bt_ordered_stop(paste0("Ordered coefficient sources for '", owner, "' are unavailable."))
      values[, ordered_spec$coefficient_names] <- replay
    }
    values <- .bt_formula_materialize_state(values, constants, validate = required)
  }
  retained_owners <- union(retained_owners, intersect(condition_event$conditional, names(priors)))
  # Gate rows come from the same original eligible population as the output.
  gate_owners <- retained_owners[vapply(priors[retained_owners], function(prior){
    is.prior.mixture(prior) || is.prior.spike_and_slab(prior)
  }, logical(1))]
  gate_columns <- if(length(gate_owners)) paste0(gate_owners, "_indicator") else character()
  required <- union(required, intersect(gate_columns, colnames(values)))
  missing <- setdiff(required, colnames(values))
  if(length(missing)) .bt_formula_transform_stop(
    "Selected formula sources are unavailable in the supplied fitted state.",
    reason = "missing_multiplier_state", missing = missing)
  if(any(!is.finite(values[, required, drop = FALSE]))) .bt_formula_transform_stop(
    "Selected formula sources are nonfinite in the supplied fitted state.", reason = "nonfinite_transform")
  record_priors <- priors[retained_owners]
  gate_plan <- if(all(gate_columns %in% colnames(values))){
    raw <- .bt_formula_raw_children(values, record_priors)
    .posterior_atoms_formula_plan(raw, record_priors)
  }else NULL
  if(!is.null(gate_plan)) gate_plan$draw_index <- as.integer(draw_index)
  record <- list(prior_list = record_priors, formula_scale = scales,
    required = required, gate_plan = gate_plan, eligible_n = nrow(values), prior_probability = 1,
    condition_event = condition_event)
  state <- list(schema_version = 1L, models = list(record),
    model = rep(1L, nrow(values)), draw_index = as.integer(draw_index),
    values = values[, required, drop = FALSE], posterior_model_probabilities = 1)
  reason <- .bt_formula_state_validate(state)
  if(!is.null(reason)) stop(reason, call. = FALSE)
  state
}

.bt_formula_state_attach <- function(samples, state){

  if(is.null(state)) return(samples)
  samples <- .bt_meta_set(samples, "formula_state", state)
  for(owner in names(samples)){
    if(!is.null(.bt_meta_get(samples[[owner]], "formula_parameter"))){
      samples[[owner]] <- .bt_meta_set(samples[[owner]], "formula_state", state)
    }
  }
  samples
}

.bt_formula_state_get <- function(samples, parameter = NULL){

  if(!is.null(parameter)){
    child <- .bt_meta_get(samples[[parameter]], "formula_state")
    if(!is.null(child)) return(child)
  }
  state <- .bt_meta_get(samples, "formula_state")
  children <- Filter(Negate(is.null), lapply(samples, function(x) .bt_meta_get(x, "formula_state")))
  if(length(children) && all(vapply(children, function(child){
    identical(child$model, children[[1L]]$model) && identical(child$draw_index, children[[1L]]$draw_index)
  }, logical(1)))) state <- children[[1L]]
  state
}

.bt_formula_dynamic_quantities <- function(x, transform, coordinates){

  dynamic <- transform$targets$map_type[match(coordinates, transform$targets$target)] == "state_dependent"
  for(field in c("quantities", "original_scale_quantities")){
    quantities <- .bt_meta_get(x, field)
    if(is.null(quantities) || nrow(quantities) != length(coordinates)) next
    for(row in which(dynamic)){
      quantities$quantity_id[[row]] <- ""
      quantities$dependencies[[row]] <- character()
      quantities$weights[[row]] <- numeric()
    }
    x <- .bt_meta_set(x, field, quantities)
  }
  x
}

# Reuse marginal_posterior's row construction and ordered term mapping. Only
# the caller's effective numeric scale changes the at data; declarations and
# contrast ownership come from the immutable fitted owner.
.bt_formula_fitted_rows <- function(formula, data, record, prefix, original = FALSE){

  if(isTRUE(record$absent)) return(matrix(0, nrow(data), length(record$zero_coordinates),
    dimnames = list(NULL, record$zero_coordinates)))
  scale <- record$formula_scale[[prefix]]
  spec <- attr(scale, "unscale_design", exact = TRUE)
  .bt_formula_unscale_design_spec_check(spec, prefix)
  fitted_data <- data
  for(variable in if(original) intersect(spec$continuous, names(fitted_data)) else character()){
    numeric_scale <- scale[[paste0(prefix, "_", variable)]]
    if(!is.null(numeric_scale)) fitted_data[[variable]] <-
      (fitted_data[[variable]] - numeric_scale$mean) / numeric_scale$sd
  }
  for(variable in intersect(names(spec$factor_levels), names(fitted_data))){
    fitted_data[[variable]] <- factor(fitted_data[[variable]],
      levels = spec$factor_levels[[variable]], ordered = spec$factor_ordered[[variable]])
    fitted_data[[variable]] <- stats::`contrasts<-`(fitted_data[[variable]],
      how.many = ncol(spec$contrast_matrices[[variable]]), value = spec$contrast_matrices[[variable]])
  }
  frame <- stats::model.frame(formula, data = fitted_data, na.action = NULL)
  matrix <- .bt_model_matrix(frame, data = frame, formula = formula)
  matrix[is.na(matrix)] <- 0
  terms <- stats::terms(formula)
  labels <- attr(terms, "term.labels")
  assign <- attr(matrix, "assign")
  columns <- names(spec$multipliers)
  out <- matrix(0, nrow(data), length(columns), dimnames = list(NULL, columns))
  intercept <- paste0(prefix, "_intercept")
  if(attr(terms, "intercept") == 1L && intercept %in% columns) out[, intercept] <- 1
  for(i in seq_along(labels)){
    owner <- JAGS_parameter_names(labels[[i]], formula_parameter = prefix)
    if(owner %in% record$zero_owners) next
    prior <- record$prior_list[[owner]]
    if(is.null(prior)) .bt_formula_transform_stop(
      "The formula row needs an omitted coefficient owner.", reason = "missing_multiplier_state", state = owner)
    info <- lapply(c("ordered", "factor_terms", "factor_design", "level_names", "levels"), function(name) attr(prior, name, exact = TRUE))
    names(info) <- c("ordered", "factor_terms", "factor_design", "level_names", "levels")
    info$ordered <- is.prior.ordered(prior)
    if(is.prior.factor(prior)) info$levels <- .get_prior_factor_levels(prior)
    term_data <- .marginal_posterior_term_data(matrix, assign, i, fitted_data, info, owner)
    term_columns <- .prior_linear_prior_columns(owner, prior)
    if(ncol(term_data) != length(term_columns)) stop("Persisted formula row columns disagree with their coefficient owner.", call. = FALSE)
    out[, term_columns] <- term_data
  }
  out
}

.bt_formula_predictor_state <- function(state, formula, data, prefix, log_intercept, original = FALSE){

  out <- matrix(0, nrow(data), nrow(state$values))
  for(model in unique(state$model)){
    rows <- which(state$model == model)
    record <- state$models[[model]]
    if(isTRUE(record$absent)){
      if(log_intercept) .bt_formula_transform_stop(
        "An absent formula model has no declared positive log-intercept source.", reason = "missing_multiplier_state")
      next
    }
    weights <- .bt_formula_fitted_rows(formula, data, record, prefix, original)
    scale <- record$formula_scale[[prefix]]
    spec <- attr(scale, "unscale_design", exact = TRUE)
    values <- .bt_formula_materialize_state(state$values[rows, , drop = FALSE], .bt_formula_state_constants(scale))
    for(column in colnames(weights)[colSums(weights != 0) > 0]){
      constant <- .bt_formula_multiplier_constant(spec$multipliers[[column]], .bt_formula_state_constants(scale))
      point <- .bt_formula_state_constants(scale)[column]
      if((!is.null(constant) && constant == 0) ||
         (length(point) && !is.na(point) && point == 0 && !identical(column, paste0(prefix, "_intercept")))) next
      if(!column %in% colnames(values)) .bt_formula_transform_stop(
        "A selected formula coefficient state is unavailable.", reason = "missing_multiplier_state", state = column)
      value <- values[, column]
      if(log_intercept && identical(column, paste0(prefix, "_intercept"))) value <- log(value)
      if(is.null(constant)){
        name <- spec$multipliers[[column]]$name
        if(!name %in% colnames(values)) .bt_formula_transform_stop(
          "A selected formula multiplier state is unavailable.", reason = "missing_multiplier_state", state = name)
        multiplier <- values[, name]
      }else multiplier <- constant
      contribution <- weights[, column, drop = FALSE] %*% matrix(value * multiplier, nrow = 1L)
      if(any(!is.finite(contribution))) .bt_formula_transform_stop(
        "Formula contribution arithmetic is nonfinite.", reason = "nonfinite_transform")
      out[, rows] <- out[, rows, drop = FALSE] + contribution
    }
  }
  if(any(!is.finite(out))) .bt_formula_transform_stop("Formula contribution arithmetic is nonfinite.", reason = "nonfinite_transform")
  out
}

.bt_formula_contribution_context <- function(record, prefix, n_grid = .prior_linear_density_default_grid(),
                                              priors = record$prior_list, active_columns = NULL,
                                              condition_event = record$condition_event){

  if(isTRUE(record$absent)){
    context <- .prior_density_context(record$prior_list, record$zero_coordinates, n_grid = n_grid)
    context$linear_weight_space <- "formula_contribution"
    return(context)
  }
  scale <- record$formula_scale[[prefix]]
  spec <- attr(scale, "unscale_design", exact = TRUE)
  if(is.null(active_columns)) active_columns <- names(spec$multipliers)
  for(column in intersect(active_columns, names(spec$multipliers))){
    declaration <- spec$multipliers[[column]]
    owner <- names(priors)[vapply(names(priors), function(owner){
      column %in% .prior_linear_prior_columns(owner, priors[[owner]])
    }, logical(1))]
    if(length(owner) == 1L && is.prior.point(priors[[owner]]) &&
       is.numeric(priors[[owner]]$parameters$location) && isTRUE(priors[[owner]]$parameters$location == 0)) next
    if(identical(declaration$type, "state") &&
       is.null(.bt_formula_multiplier_constant(declaration, .bt_formula_state_constants(scale))) &&
       is.null(priors[[declaration$name]])) .bt_formula_density_stop(
      "A compiled formula multiplier has no supplied prior law.",
      parameter = prefix, reason = "missing_multiplier_law", state = declaration$name)
  }
  context <- .prior_density_build_context(priors, names(spec$multipliers), n_grid = n_grid,
    conditional = condition_event$conditional, conditional_rule = if(is.null(condition_event$conditional_rule)) "AND" else condition_event$conditional_rule,
    condition_event = condition_event)
  context <- .bt_formula_prior_density_context_contributions(context,
    list(parameter = prefix, source_names = names(spec$multipliers), multipliers = spec$multipliers,
         state_constants = .bt_formula_state_constants(scale)))
  context$formula_scale <- NULL
  context$transforms <- list()
  context
}

.bt_formula_route_atom_certificate <- function(route){

  if(route$type %in% c("log_scale_product", "truncated_normal_convolution")){
    return(list(type = "atom_free", location = NULL))
  }
  if(identical(route$type, "conditional_normal")){
    spec <- route$spec
    if(spec$additive_sd > 0 || (spec$product_sd > 0 && .prior_density_simple_continuous(spec$multiplier))){
      return(list(type = "atom_free", location = NULL))
    }
  }
  provenance <- .prior_density_route_provenance(route)
  point <- .prior_density_ordinate_provenance_constant(provenance)
  if(!is.null(point)) return(list(type = "point", location = point))
  atoms <- .prior_density_ordinate_provenance_atoms(provenance)
  if(!is.null(atoms) && length(atoms) == 0L) return(list(type = "atom_free", location = NULL))
  list(type = "unavailable", location = NULL)
}

.bt_formula_route_numerical_scale_reason <- function(route){

  if(identical(route$type, "unknown") &&
     identical(route$provenance$kind, "numerical_scale_unavailable")) return(route$reason)
  if(identical(route$type, "transform")) return(.bt_formula_route_numerical_scale_reason(route$source))
  if(identical(route$type, "mixture")){
    reasons <- lapply(route$components[route$weights > 0], .bt_formula_route_numerical_scale_reason)
    reasons <- Filter(Negate(is.null), reasons)
    if(length(reasons)) return(reasons[[1L]])
  }
  NULL
}

.bt_formula_route_support <- function(route){

  if(identical(route$type, "log_scale_product")){
    return(.posterior_support_log(.bt_formula_route_support(route$product)))
  }
  if(identical(route$type, "truncated_normal_convolution")){
    return(.posterior_support_new(c(-Inf, Inf), source = "formula_contribution"))
  }
  if(identical(route$type, "mixture")){
    supports <- lapply(route$components[route$weights > 0], .bt_formula_route_support)
    if(any(vapply(supports, is.null, logical(1)))) return(NULL)
    return(.posterior_support_union(supports, source = "formula_contribution"))
  }
  if(identical(route$type, "transform")){
    return(.posterior_support_transform(.bt_formula_route_support(route$source),
      route$transformation, route$arguments))
  }
  if(identical(route$type, "conditional_normal")){
    certificate <- .bt_formula_route_atom_certificate(route)
    if(!identical(certificate$type, "unavailable")) return(.posterior_support_new(c(-Inf, Inf), source = "formula_contribution"))
  }
  provenance <- .prior_density_route_provenance(route)
  support <- .prior_density_ordinate_provenance_support(provenance)
  if(is.null(support)) return(NULL)
  points <- .prior_density_ordinate_provenance_atoms(provenance)
  .posterior_support_new(support, points = if(is.null(points)) numeric() else points,
    source = "formula_contribution")
}

.bt_formula_contribution_metadata <- function(marginal, state, formula, data, prefix,
                                               column, source_transforms = NULL,
                                               prior_samples = FALSE, n_grid = .prior_linear_density_default_grid(),
                                               coefficient = NULL, original = FALSE, prior_density = NULL){

  laws <- model_supports <- model_atoms <- row_weights <- vector("list", length(state$models))
  prior_available <- atoms_available <- support_available <- TRUE
  prior_reason <- atom_reason <- support_reason <- "unsupported_contribution_measure"
  prior_cause <- atom_cause <- support_cause <- NULL
  component_supports <- list()
  component_keys <- list()
  n_rows <- if(is.null(coefficient)) nrow(data) else 1L
  component_index <- rep(NA_integer_, nrow(state$values) * n_rows)
  for(model in seq_along(state$models)){
    record <- state$models[[model]]
    output_transform <- NULL
    if(is.null(coefficient)){
      weights <- .bt_formula_fitted_rows(formula, data, record, prefix, original)
      recipe_type <- "contribution_affine"
      transform <- NULL
    }else if(coefficient %in% record$zero_coordinates){
      columns <- if(isTRUE(record$absent)) record$zero_coordinates else
        names(attr(record$formula_scale[[prefix]], "unscale_design", exact = TRUE)$multipliers)
      weights <- matrix(0, 1L, length(columns), dimnames = list(NULL, columns))
      recipe_type <- "raw_affine"
      transform <- NULL
    }else{
      scale <- record$formula_scale[[prefix]]
      spec <- attr(scale, "unscale_design", exact = TRUE)
      transform <- .bt_formula_coefficient_transform(names(spec$multipliers), scale, prefix,
        log_intercept = isTRUE(attr(scale, "log_intercept", exact = TRUE)))
      recipe <- transform$prior_recipes[[coefficient]]
      if(is.null(recipe)) stop("Unknown coefficient measure target.", call. = FALSE)
      if(identical(recipe$type, "unavailable")){
        if(record$prior_probability > 0){ prior_available <- FALSE; prior_reason <- recipe$reason }
        if(state$posterior_model_probabilities[[model]] > 0){
          atoms_available <- support_available <- FALSE
          atom_reason <- support_reason <- recipe$reason
        }
        next
      }
      weights <- matrix(recipe$weights, nrow = 1L, dimnames = list(NULL, names(recipe$weights)))
      recipe_type <- recipe$type
      source_transforms <- transform$source_transforms[transform$source_transforms == "log"]
      if(length(source_transforms) == 0L) source_transforms <- NULL
      if(identical(transform$output_transforms[[coefficient]], "exp")) output_transform <- "exp"
    }
    row_weights[[model]] <- weights
    context_builder <- function(priors, conditional = TRUE){
      if(identical(recipe_type, "contribution_affine")) return(.bt_formula_contribution_context(record, prefix, n_grid, priors,
        colnames(weights)[colSums(weights != 0) > 0], condition_event = if(conditional) record$condition_event))
      event <- if(conditional) record$condition_event
      .bt_formula_prior_density_context_raw_sources(
        .prior_density_build_context(priors, colnames(weights), n_grid = n_grid,
          conditional = event$conditional, conditional_rule = if(is.null(event$conditional_rule)) "AND" else event$conditional_rule,
          condition_event = event), colnames(weights))
    }
    build_route <- function(context, row){
      .prior_density_route_context(context, weights[row, ], source_transforms,
        output_transform, NULL)
    }
    context <- tryCatch(context_builder(record$prior_list),
      BayesTools_formula_measure_unavailable = function(e) e)
    if(inherits(context, "BayesTools_formula_measure_unavailable")){
      detail <- if(is.null(context$detail)) conditionMessage(context) else context$detail
      if(record$prior_probability > 0){
        prior_available <- FALSE
        prior_reason <- detail
        prior_cause <- context$reason
      }
      if(state$posterior_model_probabilities[[model]] > 0){
        atoms_available <- support_available <- FALSE
        atom_reason <- support_reason <- detail
        atom_cause <- support_cause <- context$reason
      }
      next
    }
    routes <- lapply(seq_len(nrow(weights)), function(row){
      build_route(context, row)
    })
    numerical_reasons <- Filter(Negate(is.null), lapply(routes, .bt_formula_route_numerical_scale_reason))
    if(record$prior_probability > 0 && length(numerical_reasons)){
      prior_available <- FALSE
      prior_reason <- numerical_reasons[[1L]]
      prior_cause <- "numerical_scale_unavailable"
    }else if(record$prior_probability > 0 && any(vapply(routes, function(route) identical(route$type, "unknown"), logical(1)))){
      prior_available <- FALSE
    }else if(prior_samples && record$prior_probability > 0){
      laws[[model]] <- if(length(state$models) == 1L && !is.null(prior_density)) prior_density else
        .prior_density_from_context_rows(context, weights, source_transforms = source_transforms,
          output_transformation = output_transform)
    }
    if(state$posterior_model_probabilities[[model]] <= 0) next
    plan <- record$gate_plan
    if(is.null(plan)){
      atoms_available <- support_available <- FALSE
      atom_reason <- support_reason <- "missing_joint_component_states"
      next
    }
    locations <- masses <- numeric()
    supports <- list()
    for(component in which(plan$probabilities > 0)){
      this_component_supports <- list()
      component_priors <- record$prior_list
      for(owner in colnames(plan$components)) component_priors[[owner]] <- .posterior_atoms_component_prior(
        component_priors[[owner]], plan$components[component, owner], FALSE,
        .posterior_atoms_plan_total_indicator(plan, component, owner))
      component_context <- context_builder(component_priors, conditional = FALSE)
      for(row in seq_len(nrow(weights))){
        route <- build_route(component_context, row)
        numerical_reason <- .bt_formula_route_numerical_scale_reason(route)
        if(!is.null(numerical_reason)){
          atom_reason <- support_reason <- numerical_reason
          atom_cause <- support_cause <- "numerical_scale_unavailable"
        }
        certificate <- .bt_formula_route_atom_certificate(route)
        if(identical(certificate$type, "unavailable")) atoms_available <- FALSE
        if(identical(certificate$type, "point")){
          locations <- c(locations, certificate$location)
          masses <- c(masses, state$posterior_model_probabilities[[model]] * plan$probabilities[[component]] / nrow(weights))
        }
        support <- .bt_formula_route_support(route)
        if(is.null(support)) support_available <- FALSE else{
          supports[[length(supports) + 1L]] <- support
          this_component_supports[[length(this_component_supports) + 1L]] <- support
        }
      }
      component_supports[length(component_supports) + 1L] <- list(if(length(this_component_supports) == nrow(weights))
        .posterior_support_union(this_component_supports, source = "formula_joint_component") else NULL)
      global_component <- length(component_supports)
      key <- stats::setNames(as.integer(plan$components[component, ]), colnames(plan$components))
      if(length(state$models) > 1L) key <- c(.component = component)
      component_keys[[global_component]] <- if(length(state$models) > 1L) c(.model = model, key) else key
      model_rows <- which(state$model == model)
      original_rows <- if(isTRUE(record$absent)) rep(1L, length(model_rows)) else
        match(state$draw_index[model_rows], plan$draw_index)
      if(anyNA(original_rows)) stop("Formula component rows do not belong to their original eligible population.", call. = FALSE)
      selected <- model_rows[plan$index[original_rows] == component]
      if(length(selected)){
        output_rows <- unlist(lapply(selected, function(row) (row - 1L) * n_rows + seq_len(n_rows)), use.names = FALSE)
        component_index[output_rows] <- global_component
      }
    }
    model_atoms[[model]] <- data.frame(x = locations, mass = masses)
    if(length(supports)) model_supports[[model]] <- .posterior_support_union(supports, source = "formula_contribution")
  }
  marginal <- .bt_meta_assign(marginal, list(support = NULL, atoms = NULL, components = NULL,
    prior_density = NULL, prior_context = NULL, linear_weights = NULL, linear_weight_space = NULL))
  if(prior_available && prior_samples){
    positive <- which(vapply(state$models, function(record) record$prior_probability > 0, logical(1)))
    probabilities <- vapply(state$models[positive], `[[`, numeric(1), "prior_probability")
    probabilities <- probabilities / sum(probabilities)
    dx <- min(vapply(laws[positive], function(law){
      if(!is.null(law$density) && length(law$density$x) > 1L) diff(law$density$x)[[1L]] else Inf
    }, numeric(1)))
    if(!is.finite(dx)) dx <- NA_real_
    law <- if(length(positive) == 1L) laws[[positive]] else
      .prior_linear_density_mix(laws[positive], probabilities, dx = dx, n_grid = n_grid)
    if(length(positive) > 1L) attr(law, "adaptive_evaluation") <- list(kind = "density_mixture",
      arguments = list(dists = laws[positive], weights = probabilities, n_grid = n_grid))
    marginal <- .bt_meta_set(marginal, "prior_density", law)
  }else if(!prior_available) marginal <- .bt_formula_measure_mark(marginal, column, "prior_density", prior_reason, cause = prior_cause)
  if(atoms_available){
    points <- .posterior_atoms_point_mass_table(do.call(rbind, model_atoms))
    atoms <- .posterior_atoms_new(locations = matrix(points$x, ncol = 1L), mass = points$mass,
      column_names = column, source = "formula_joint_component_frequencies")
    marginal <- .posterior_atoms_set(marginal, atoms)
  }else marginal <- .bt_formula_measure_mark(marginal, column, "atoms", atom_reason, cause = atom_cause)
  if(support_available){
    marginal <- .posterior_support_set(marginal, .posterior_support_union(
      model_supports[state$posterior_model_probabilities > 0], source = "formula_contribution"))
  }else marginal <- .bt_formula_measure_mark(marginal, column, "support", support_reason, cause = support_cause)
  if(!anyNA(component_index) && length(component_supports) &&
     (length(state$models) > 1L || any(vapply(state$models, function(record){
       any(vapply(record$prior_list, .posterior_components_is_mixture, logical(1)))
     }, logical(1))))){
    marginal <- .posterior_components_set(marginal,
      .posterior_components_new(component_index, component_supports, do.call(rbind, component_keys)))
  }
  common <- row_weights[[1L]]
  if(is.null(coefficient) && prior_available && all(vapply(row_weights, identical, logical(1), common)) && length(state$models) == 1L){
    context <- .bt_formula_contribution_context(state$models[[1L]], prefix, n_grid,
      active_columns = colnames(common)[colSums(common != 0) > 0])
    marginal <- .bt_meta_assign(marginal, list(prior_context = context, linear_weights = common,
      linear_weight_space = "formula_contribution"))
  }
  marginal
}

.bt_formula_density_stop <- function(message, ...){

  fields <- list(...)
  measure <- fields$reason %in% .bt_formula_measure_causes
  condition <- structure(
    c(list(message = message, call = NULL), fields),
    class = c(
      "BayesTools_formula_prior_density_unavailable",
      if(isTRUE(measure)) "BayesTools_formula_measure_unavailable",
      "BayesTools_formula_transform_unavailable",
      "error",
      "condition"
    )
  )
  stop(condition)
}

.bt_formula_measure_causes <- c("state_dependent_map", "nonlinear_map",
  "incompatible_prior_recipes", "missing_multiplier_law", "numerical_scale_unavailable",
  "structural_target_law_unavailable", "unsupported_contribution_measure")
