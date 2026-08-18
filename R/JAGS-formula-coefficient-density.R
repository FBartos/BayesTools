# Versioned fixed-coefficient transforms and induced prior densities.

.bt_formula_coefficient_transform_version <- 1L

#' Formula coefficient transformations and induced prior densities
#'
#' @description
#' `JAGS_formula_coefficient_transform()` exposes the exact fixed-coefficient
#' transformation used by [transform_scale_samples()]. It records fitted
#' sources, original-scale targets, named source/output transforms, exact
#' nonzero dependencies, and parameter-map structural metadata.
#'
#' `JAGS_formula_prior_density()` constructs the induced marginal prior measure
#' for one target through BayesTools' deterministic prior-density context and
#' linear-density algebra. The result can be passed directly to
#' [prior_density_ordinate()].
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
#' @param context optional BayesTools prior-density context. This is useful for
#'   product-space or conditional prior mixtures; when `NULL`, the fitted
#'   `prior_list` is used.
#'
#' @return `JAGS_formula_coefficient_transform()` returns a
#' `BayesTools_formula_coefficient_transform` list.
#' `JAGS_formula_coefficient_transform_schema()` returns field descriptions.
#' `JAGS_formula_prior_density()` returns a `prior_linear_density` accepted by
#' [prior_density_ordinate()].
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
          stop(
            "LKJ primitive coordinate '", coordinate_name,
            "' is missing from the fitted parameter map. Refit the model with the current BayesTools version.",
            call. = FALSE
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
  design <- JAGS_formula_design(fit, parameter = parameter)
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
  coordinates <- parameter_coordinates(fit)
  sources <- .bt_formula_coefficient_sources(design, coordinates, parameter)

  .bt_formula_coefficient_transform(
    source_names = sources$source,
    formula_scale = design$formula_scale,
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
      "source_names", "target_names", "matrix", "source_transforms",
      "output_transforms", "dependencies", "sources", "targets"
    ),
    type = c(
      rep("integer", 3L), rep("character", 4L), "numeric matrix",
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
      "Exact target-by-source additive coefficient matrix.",
      "Named identity/log transform for every fitted source.",
      "Named identity/exp transform for every target.",
      "One row per exact nonzero target/source coefficient.",
      "Registry-linked source monitor status and fixed value.",
      "Target structural/dependent status and exact fixed value."
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname JAGS_formula_coefficient_transform
JAGS_formula_prior_density <- function(
    fit, parameter, target, target_scale = "original", context = NULL){

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
      reason = target_metadata$reason
    )
  }

  weights <- stats::setNames(
    as.numeric(transform$matrix[target_i, , drop = FALSE]),
    colnames(transform$matrix)
  )
  weights <- weights[weights != 0]
  density_context <- .bt_formula_prior_density_context(
    fit,
    context = context,
    source_names = transform$source_names,
    parameter = parameter,
    target = target
  )
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

.bt_formula_coefficient_transform <- function(
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
  sources <- data.frame(
    source = source_names,
    monitor_status = source_metadata$monitor_status,
    fixed_value = as.numeric(source_metadata$fixed_value),
    source_transform = unname(source_transforms),
    stringsAsFactors = FALSE
  )
  dependencies <- .bt_formula_coefficient_dependencies(matrix)
  targets <- .bt_formula_coefficient_targets(
    matrix,
    sources,
    output_transforms
  )

  out <- list(
    schema_version = .bt_formula_coefficient_transform_version,
    formula_design_version = as.integer(formula_design_version),
    parameter_map_version = as.integer(parameter_map_version),
    parameter = parameter,
    target_scale = target_scale,
    source_names = source_names,
    target_names = rownames(matrix),
    matrix = matrix,
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
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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

.bt_apply_formula_coefficient_transform <- function(samples, transform){

  .bt_validate_formula_coefficient_transform(transform)
  samples <- as.matrix(samples)
  missing <- setdiff(transform$source_names, colnames(samples))
  if(length(missing) > 0L){
    stop(
      "Coefficient samples are missing fitted source coordinate",
      if(length(missing) > 1L) "s " else " ",
      paste0("'", missing, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  source <- samples[, transform$source_names, drop = FALSE]
  log_sources <- transform$source_transforms == "log"
  if(any(log_sources)){
    if(any(source[, log_sources, drop = FALSE] <= 0)){
      stop("Log-transformed fitted coefficient samples must be positive.",
           call. = FALSE)
    }
    source[, log_sources] <- log(source[, log_sources, drop = FALSE])
  }
  target <- source %*% t(transform$matrix)
  exp_targets <- transform$output_transforms == "exp"
  if(any(exp_targets)){
    target[, exp_targets] <- exp(target[, exp_targets, drop = FALSE])
  }
  samples[, transform$target_names] <- target
  samples
}

.bt_validate_formula_coefficient_transform <- function(transform){

  expected_names <- c(
    "schema_version", "formula_design_version",
    "parameter_map_version", "parameter", "target_scale",
    "source_names", "target_names", "matrix", "source_transforms",
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
    all(is.finite(transform$matrix)) &&
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
    stop("Formula coefficient transform metadata are missing or unsupported. Rebuild the transform with this version of BayesTools.",
         call. = FALSE)
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
  expected_targets <- .bt_formula_coefficient_targets(
    transform$matrix,
    transform$sources,
    transform$output_transforms
  )
  if(!valid_sources || !identical(transform$dependencies,
                                  expected_dependencies) ||
     !identical(transform$targets, expected_targets)){
    stop("Formula coefficient transform structural metadata are malformed. Rebuild the transform with this version of BayesTools.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.bt_formula_prior_density_context <- function(
    fit, context, source_names, parameter, target){

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
  out
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

.bt_formula_density_stop <- function(message, ...){

  condition <- structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(
      "BayesTools_formula_prior_density_unavailable",
      "BayesTools_formula_transform_unavailable",
      "error",
      "condition"
    )
  )
  stop(condition)
}
