# Exact fixed-formula predictor dependencies and update bases.

#' Fitted-formula coordinate dependencies
#'
#' @description
#' `JAGS_formula_coordinate_dependencies()` reports every fitted formula
#' predictor that structurally depends on the supplied posterior coordinates.
#' It uses the persisted formula designs and parameter map, including direct
#' coefficients, formula expressions, coefficient multipliers, and external
#' random-effect SD sources. A row-source callback is reported as opaque for
#' every queried coordinate because its R code may inspect arbitrary posterior
#' parameters.
#'
#' @param fit fitted object created by [JAGS_fit()].
#' @param coordinates unique fitted coordinate names.
#'
#' @return A data frame with `coordinate_name`, `formula_parameter`, and
#'   `dependency_type` columns. An empty data frame means that none of the
#'   fitted formula predictors depends on the supplied coordinates.
#'
#' @seealso [JAGS_formula_design()] [parameter_map()]
#' @export
JAGS_formula_coordinate_dependencies <- function(fit, coordinates){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  JAGS_validate_fit_contract(
    fit,
    requires = c("formula_design", "parameter_map")
  )
  if(!is.character(coordinates) || length(coordinates) < 1L ||
     anyNA(coordinates) || any(!nzchar(coordinates)) ||
     anyDuplicated(coordinates)){
    stop("'coordinates' must contain unique non-empty fitted coordinate names.",
         call. = FALSE)
  }

  coordinate_map <- parameter_coordinates(fit)
  unknown <- setdiff(coordinates, coordinate_map$coordinate_name)
  if(length(unknown) > 0L){
    stop(
      "Formula dependency lookup references unknown fitted coordinate",
      if(length(unknown) > 1L) "s " else " ",
      paste0("'", unknown, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  designs <- JAGS_formula_design(fit)
  rows <- list()
  row_i <- 0L
  for(parameter in names(designs)){
    design <- designs[[parameter]]
    .bt_validate_formula_design_replay_schema(
      design,
      context = paste0("Formula dependencies for parameter '", parameter, "'")
    )
    dependencies <- .bt_formula_predictor_design_dependencies(design)
    opaque <- .bt_formula_predictor_has_opaque_random_dependency(design)
    for(coordinate_name in coordinates){
      coordinate_row <- coordinate_map[
        coordinate_map$coordinate_name == coordinate_name,
        ,
        drop = FALSE
      ]
      types <- character()
      if(identical(coordinate_row$role, "fixed_coefficient") &&
         identical(coordinate_row$formula_parameter, parameter)){
        types <- c(types, "coefficient")
      }
      coordinate_base <- .bt_parameter_coordinates_base(coordinate_name)
      for(type in names(dependencies)){
        dependency_bases <- .bt_parameter_coordinates_base(dependencies[[type]])
        if(coordinate_base %in% dependency_bases){
          types <- c(types, type)
        }
      }
      if(opaque){
        types <- c(types, "opaque_random_sd_callback")
      }
      for(type in unique(types)){
        row_i <- row_i + 1L
        rows[[row_i]] <- data.frame(
          coordinate_name = coordinate_name,
          formula_parameter = parameter,
          dependency_type = type,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if(length(rows) == 0L){
    return(data.frame(
      coordinate_name = character(),
      formula_parameter = character(),
      dependency_type = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Exact fitted-formula predictor update basis
#'
#' @description
#' `JAGS_formula_predictor_basis()` uses the fitted parameter map and persisted
#' formula design to construct the exact observation-level change in a fixed
#' linear predictor for a supplied direction through fitted coefficient
#' coordinates. The result is conditional on all coordinates outside that
#' direction, including fitted random effects.
#'
#' The accessor reports `"affine"` only when the stored metadata proves that
#' the selected update is additive on the fitted predictor scale. Logged
#' intercepts and coordinates reused by formula expressions, coefficient
#' multipliers, or random-effect scale sources are reported as non-affine or
#' unsupported so callers can retain a generic evaluator.
#'
#' @param fit fitted object created by [JAGS_fit()].
#' @param directions a fully named numeric vector, or a numeric matrix with
#'   unique coefficient-coordinate column names. Matrix rows define separate
#'   posterior-row update directions.
#' @param posterior_samples optional posterior sample matrix aligned with the
#'   rows of `directions`. It is required when a selected formula term has a
#'   parameter-valued `multiply_by` attribute; when omitted, draws are obtained
#'   from `fit`.
#'
#' @return A `BayesTools_formula_predictor_basis` list with `status`, `reason`,
#' formula `parameter`, selected `coordinates`, and an `S` by `N` numeric
#' `basis`. `basis` is present only when `status` is `"affine"`.
#'
#' @seealso [JAGS_formula_design()] [parameter_map()]
#' @export
JAGS_formula_predictor_basis <- function(fit, directions,
                                         posterior_samples = NULL){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  JAGS_validate_fit_contract(
    fit,
    requires = c("formula_design", "parameter_map")
  )
  directions <- .bt_formula_predictor_directions(directions)
  coordinate_names <- colnames(directions)
  coordinates <- parameter_coordinates(fit)
  coordinate_rows <- match(coordinate_names, coordinates$coordinate_name)
  if(anyNA(coordinate_rows)){
    missing <- coordinate_names[is.na(coordinate_rows)]
    stop(
      "Formula predictor directions reference unknown fitted coordinate",
      if(length(missing) > 1L) "s " else " ",
      paste0("'", missing, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  selected <- coordinates[coordinate_rows, , drop = FALSE]
  if(any(selected$role != "fixed_coefficient") || any(selected$internal)){
    return(.bt_formula_predictor_basis_result(
      status = "unsupported",
      reason = paste(
        "Selected directions include coordinates that are not public",
        "fixed-formula coefficients."
      ),
      coordinates = coordinate_names
    ))
  }
  parameters <- unique(selected$formula_parameter)
  if(length(parameters) != 1L || is.na(parameters) || !nzchar(parameters)){
    return(.bt_formula_predictor_basis_result(
      status = "unsupported",
      reason = "Selected coefficient coordinates do not have one formula owner.",
      coordinates = coordinate_names
    ))
  }
  parameter <- parameters[[1L]]
  design <- JAGS_formula_design(fit, parameter = parameter)
  .bt_validate_formula_design_replay_schema(
    design,
    context = paste0("Predictor basis for parameter '", parameter, "'")
  )
  coordinate_map <- .bt_formula_predictor_coordinate_map(
    design = design,
    coordinates = selected
  )
  moving <- vapply(seq_len(ncol(directions)), function(i){
    any(directions[, i] != 0)
  }, logical(1))
  moving_coordinates <- coordinate_names[moving]

  intercept <- paste0(parameter, "_intercept")
  if(isTRUE(design$log_intercept) && intercept %in% moving_coordinates){
    return(.bt_formula_predictor_basis_result(
      status = "non_affine",
      reason = "The selected fitted intercept enters the predictor through log().",
      parameter = parameter,
      coordinates = coordinate_names
    ))
  }
  if(length(moving_coordinates) > 0L){
    dependencies <- JAGS_formula_coordinate_dependencies(
      fit,
      coordinates = moving_coordinates
    )
    indirect <- dependencies$dependency_type != "coefficient"
    if(any(indirect)){
      dependency_types <- unique(dependencies$dependency_type[indirect])
      reason <- .bt_formula_predictor_dependency_reason(dependency_types)
      return(.bt_formula_predictor_basis_result(
        status = "unsupported",
        reason = reason,
        parameter = parameter,
        coordinates = coordinate_names
      ))
    }
  }

  multiplier_names <- unique(coordinate_map$multiplier[
    !is.na(coordinate_map$multiplier) & nzchar(coordinate_map$multiplier)
  ])
  if(length(multiplier_names) > 0L && is.null(posterior_samples)){
    posterior_samples <- as.matrix(.fit_to_posterior(fit))
  }else if(!is.null(posterior_samples)){
    posterior_samples <- as.matrix(posterior_samples)
  }
  if(!is.null(posterior_samples)){
    if(nrow(directions) == 1L && nrow(posterior_samples) > 1L){
      directions <- directions[rep(1L, nrow(posterior_samples)), , drop = FALSE]
    }else if(nrow(directions) != nrow(posterior_samples)){
      stop(
        "'directions' and 'posterior_samples' must have the same number of rows.",
        call. = FALSE
      )
    }
  }
  if(length(multiplier_names) > 0L){
    posterior_names <- colnames(posterior_samples)
    missing <- setdiff(multiplier_names, posterior_names)
    if(length(missing) > 0L){
      return(.bt_formula_predictor_basis_result(
        status = "unsupported",
        reason = paste0(
          "Posterior samples are missing formula multiplier coordinate",
          if(length(missing) > 1L) "s " else " ",
          paste0("'", missing, "'", collapse = ", "),
          "."
        ),
        parameter = parameter,
        coordinates = coordinate_names
      ))
    }
  }

  basis <- matrix(
    0,
    nrow = nrow(directions),
    ncol = nrow(design$model_matrix),
    dimnames = list(rownames(directions), rownames(design$model_matrix))
  )
  for(i in seq_len(nrow(coordinate_map))){
    design_column <- design$model_matrix[, coordinate_map$model_column[i]]
    multiplier <- coordinate_map$multiplier_value[i]
    if(!is.na(coordinate_map$multiplier[i])){
      multiplier <- posterior_samples[, coordinate_map$multiplier[i]]
    }
    if(any(!is.finite(multiplier))){
      return(.bt_formula_predictor_basis_result(
        status = "unsupported",
        reason = "Formula multiplier values are not finite.",
        parameter = parameter,
        coordinates = coordinate_names
      ))
    }
    basis <- basis + outer(
      directions[, i] * multiplier,
      as.numeric(design_column)
    )
  }

  .bt_formula_predictor_basis_result(
    status = "affine",
    reason = "",
    parameter = parameter,
    coordinates = coordinate_names,
    basis = basis
  )
}

.bt_formula_predictor_directions <- function(directions){

  if(is.numeric(directions) && is.null(dim(directions))){
    direction_names <- names(directions)
    directions <- matrix(
      as.numeric(directions),
      nrow = 1L,
      dimnames = list(NULL, direction_names)
    )
  }
  if(!is.matrix(directions) || !is.numeric(directions) ||
     nrow(directions) < 1L || ncol(directions) < 1L ||
     any(!is.finite(directions))){
    stop(
      "'directions' must be a finite named numeric vector or numeric matrix.",
      call. = FALSE
    )
  }
  coordinate_names <- colnames(directions)
  if(is.null(coordinate_names) || anyNA(coordinate_names) ||
     any(!nzchar(coordinate_names)) || anyDuplicated(coordinate_names)){
    stop("'directions' must have unique non-empty coordinate names.",
         call. = FALSE)
  }
  directions
}

.bt_formula_predictor_coordinate_map <- function(design, coordinates){

  name_map <- design$name_map
  .bt_validate_formula_name_map(name_map)
  fixed_map <- name_map[
    name_map$kind == "fixed" &
      name_map$formula_parameter == design$parameter,
    ,
    drop = FALSE
  ]
  semantic_terms <- gsub(
    "__xXx__",
    ":",
    design$model_terms,
    fixed = TRUE
  )
  rows <- vector("list", nrow(coordinates))
  for(i in seq_len(nrow(coordinates))){
    coordinate_name <- coordinates$coordinate_name[i]
    base_name <- .bt_parameter_coordinates_base(coordinate_name)
    map_row <- fixed_map[fixed_map$jags_name == base_name, , drop = FALSE]
    if(nrow(map_row) != 1L){
      stop(
        "Formula predictor metadata do not uniquely map coordinate '",
        coordinate_name, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    term_i <- match(map_row$term, semantic_terms)
    if(is.na(term_i)){
      stop(
        "Formula predictor metadata do not map term '", map_row$term,
        "' to the fitted design. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    term_columns <- which(design$assign == (term_i - 1L))
    coefficient_i <- .bt_formula_predictor_coordinate_index(
      coordinate = coordinates$index[i],
      n_columns = length(term_columns),
      coordinate_name = coordinate_name
    )
    model_column <- term_columns[coefficient_i]
    prior <- design$prior_list[[map_row$jags_name]]
    if(is.null(prior) || !is.prior(prior)){
      stop(
        "Formula predictor metadata are missing the prior for coordinate '",
        coordinate_name, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    multiply_by <- attr(prior, "multiply_by", exact = TRUE)
    multiplier <- NA_character_
    multiplier_value <- 1
    if(is.numeric(multiply_by)){
      if(length(multiply_by) != 1L || !is.finite(multiply_by)){
        stop("Numeric formula 'multiply_by' metadata must be one finite value.",
             call. = FALSE)
      }
      multiplier_value <- as.numeric(multiply_by)
    }else if(!is.null(multiply_by)){
      if(!is.character(multiply_by) || length(multiply_by) != 1L ||
         is.na(multiply_by) || !nzchar(multiply_by)){
        stop("Formula 'multiply_by' metadata are malformed.", call. = FALSE)
      }
      multiplier <- JAGS_parameter_names(multiply_by)
      if(length(multiplier) != 1L){
        stop("Formula 'multiply_by' metadata must resolve to one coordinate.",
             call. = FALSE)
      }
    }
    rows[[i]] <- data.frame(
      coordinate = coordinate_name,
      model_column = model_column,
      multiplier = multiplier,
      multiplier_value = multiplier_value,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.bt_formula_predictor_coordinate_index <- function(coordinate, n_columns,
                                                    coordinate_name){

  if(n_columns < 1L){
    stop(
      "Formula predictor metadata have no design column for coordinate '",
      coordinate_name, "'. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  if(!nzchar(coordinate)){
    if(n_columns == 1L){
      return(1L)
    }
  }else{
    parsed <- suppressWarnings(as.integer(coordinate))
    if(length(parsed) == 1L && !is.na(parsed) &&
       parsed >= 1L && parsed <= n_columns){
      return(parsed)
    }
  }
  stop(
    "Formula predictor coordinate '", coordinate_name,
    "' does not identify one fitted design column. Refit the model with this ",
    "version of BayesTools.",
    call. = FALSE
  )
}

.bt_formula_predictor_multiplier_dependencies <- function(design){

  dependencies <- lapply(design$prior_list, function(prior){
    multiply_by <- attr(prior, "multiply_by", exact = TRUE)
    if(is.character(multiply_by) && length(multiply_by) == 1L &&
       !is.na(multiply_by) && nzchar(multiply_by)){
      return(JAGS_parameter_names(multiply_by))
    }
    character()
  })
  unique(unlist(dependencies, use.names = FALSE))
}

.bt_formula_predictor_design_dependencies <- function(design){

  list(
    expression = unique(unlist(lapply(
      design$expression_specs,
      `[[`,
      "parameter_dependencies"
    ), use.names = FALSE)),
    multiplier = .bt_formula_predictor_multiplier_dependencies(design),
    random_sd_source = .bt_formula_predictor_random_dependencies(design)
  )
}

.bt_formula_predictor_has_opaque_random_dependency <- function(design){

  for(random_term in design$random_effects){
    binding <- random_term$sd_binding
    if(is.null(binding)){
      next
    }
    .bt_check_random_sd_binding(binding)
    sources <- c(list(binding$source), binding$sources_by_column)
    for(source in sources){
      if(.bt_random_sd_binding_source_is_external(source) &&
         .bt_parameter_source_has_values(source)){
        return(TRUE)
      }
    }
  }
  FALSE
}

.bt_formula_predictor_dependency_reason <- function(types){

  if("opaque_random_sd_callback" %in% types){
    return(paste(
      "A persisted random-effect SD callback has opaque posterior-coordinate",
      "dependencies."
    ))
  }
  if("expression" %in% types){
    return(paste(
      "Selected coefficient coordinates are reused by a persisted formula",
      "expression."
    ))
  }
  if("multiplier" %in% types){
    return(paste(
      "Selected coefficient coordinates are reused as formula-term",
      "multipliers."
    ))
  }
  if("random_sd_source" %in% types){
    return(paste(
      "Selected coefficient coordinates are reused as random-effect scale",
      "sources."
    ))
  }
  "Selected coefficient coordinates have unsupported formula dependencies."
}

.bt_formula_predictor_random_dependencies <- function(design){

  dependencies <- list()
  dependency_i <- 0L
  for(random_term in design$random_effects){
    binding <- random_term$sd_binding
    if(is.null(binding)){
      next
    }
    .bt_check_random_sd_binding(binding)
    sources <- c(list(binding$source), binding$sources_by_column)
    for(source in sources){
      if(!.bt_random_sd_binding_source_is_external(source)){
        next
      }
      dependency_i <- dependency_i + 1L
      dependencies[[dependency_i]] <-
        .bt_random_sd_binding_source_name(source)
    }
  }
  unique(unlist(dependencies, use.names = FALSE))
}

.bt_formula_predictor_basis_result <- function(
    status, reason, parameter = "", coordinates = character(), basis = NULL){

  out <- list(
    status = status,
    reason = reason,
    parameter = parameter,
    coordinates = coordinates,
    basis = basis
  )
  class(out) <- c("BayesTools_formula_predictor_basis", "list")
  out
}
