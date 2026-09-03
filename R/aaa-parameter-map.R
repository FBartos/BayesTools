# Authoritative fitted parameter map.

.bt_parameter_map_version <- 4L

#' Fitted parameter map and coordinate view
#'
#' @description
#' `parameter_map()` returns the single versioned parameter map stored by
#' [JAGS_fit()]. It contains three linked tables: concrete fitted
#' `coordinates`, public semantic `quantities`, and exact `aliases`.
#'
#' `parameter_coordinates()` returns the coordinate table from that map. Each
#' row describes one sampled, structural, or unavailable fitted coordinate.
#' Coordinate names are backend identities, not additional public aliases.
#'
#' `parameter_map_schema()` documents all three tables, while
#' `parameter_coordinates_schema()` documents the coordinate view alone.
#'
#' @param object a fitted object.
#' @param ... arguments for methods.
#'
#' @return `parameter_map()` returns a `BayesTools_parameter_map`.
#' `parameter_coordinates()` returns a `BayesTools_parameter_coordinates`
#' data frame. The schema helpers return data frames or a named list of them.
#'
#' @export parameter_map
#' @export parameter_map_schema
#' @export parameter_coordinates
#' @export parameter_coordinates_schema
#' @name parameter_map
NULL

#' @rdname parameter_map
parameter_map <- function(object, ...){

  UseMethod("parameter_map")
}

#' @rdname parameter_map
#' @exportS3Method parameter_map BayesTools_fit
parameter_map.BayesTools_fit <- function(object, ...){

  if(!inherits(object, "BayesTools_fit")){
    stop("'object' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  map <- attr(object, "parameter_map", exact = TRUE)
  if(is.null(map)){
    stop(
      "The fitted object does not contain parameter-map metadata. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  .bt_validate_parameter_map(map)
  map
}

#' @rdname parameter_map
parameter_coordinates <- function(object, ...){

  UseMethod("parameter_coordinates")
}

#' @rdname parameter_map
#' @exportS3Method parameter_coordinates BayesTools_fit
parameter_coordinates.BayesTools_fit <- function(object, ...){

  parameter_map(object)$coordinates
}

#' @rdname parameter_map
parameter_coordinates_schema <- function(){

  data.frame(
    field = .bt_parameter_coordinates_columns,
    type = c(
      rep("character", 12L),
      "numeric",
      rep("character", 3L),
      "logical"
    ),
    description = c(
      "Unique concrete posterior or structural parameter name.",
      "JAGS monitor node requested during fitting.",
      "Formula output parameter owning the coefficient or random block.",
      "Coordinate role such as fixed coefficient, random SD, latent variable, allocation, or ordinary parameter.",
      "Canonical random-effect block identifier.",
      "Public random-effect owner name, distinct from the grouping label.",
      "Formula term represented by the coordinate.",
      "Concrete fixed/random design-matrix column.",
      "Comma-separated concrete JAGS array indices.",
      "Fitted array dimensions joined by 'x'.",
      "Scale used by the fitted coordinate.",
      "Whether the coordinate is sampled, structural, or unavailable.",
      "Exact prior-specified value for a structural coordinate; otherwise NA.",
      "Default coordinate display label.",
      "Random grouping-variable label.",
      "Random covariance structure.",
      "Whether the coordinate is an implementation-level dependency."
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname parameter_map
parameter_map_schema <- function(){

  catalog <- parameter_catalog_schema()
  list(
    schema_version = .bt_parameter_map_version,
    coordinates = parameter_coordinates_schema(),
    quantities = catalog$quantities,
    aliases = catalog$aliases
  )
}

.bt_parameter_map_new <- function(coordinates, quantities, aliases){

  out <- list(
    schema_version = .bt_parameter_map_version,
    coordinates = coordinates,
    quantities = quantities,
    aliases = aliases
  )
  class(out) <- c("BayesTools_parameter_map", "list")
  .bt_validate_parameter_map(out)
  out
}

.bt_parameter_map_catalog <- function(map){

  out <- list(
    schema_version = map$schema_version,
    quantities = map$quantities,
    aliases = map$aliases
  )
  class(out) <- c("BayesTools_parameter_catalog", "list")
  out
}

.bt_validate_parameter_map <- function(map){

  valid <- inherits(map, "BayesTools_parameter_map") &&
    is.list(map) &&
    identical(
      names(map),
      c("schema_version", "coordinates", "quantities", "aliases")
    ) &&
    identical(map$schema_version, .bt_parameter_map_version)
  if(!valid){
    stop(
      "Parameter-map metadata are missing, malformed, or unsupported. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }

  .bt_validate_parameter_coordinates(map$coordinates)
  .bt_validate_parameter_catalog_tables(map$quantities, map$aliases)

  native <- map$quantities$provider == "BayesTools"
  dependencies <- unique(unlist(
    lapply(map$quantities$extraction_key[native], `[[`, "dependencies"),
    use.names = FALSE
  ))
  missing <- setdiff(dependencies, map$coordinates$coordinate_name)
  if(length(missing) > 0L){
    stop(
      "Parameter-map quantities reference unknown coordinate dependencies: ",
      paste0("'", missing, "'", collapse = ", "),
      ". Refit or rebuild the map with this version of BayesTools.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_build_parameter_map <- function(columns, monitor_names = columns,
                                    prior_list = NULL,
                                    formula_design = NULL,
                                    formula_scale = NULL,
                                    backend_anchor = NULL){

  coordinates <- .bt_build_parameter_coordinates(
    columns = columns,
    monitor_names = monitor_names,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale,
    backend_anchor = backend_anchor
  )
  catalog <- .bt_build_parameter_catalog(
    coordinates = coordinates,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  .bt_parameter_map_new(
    coordinates = coordinates,
    quantities = catalog$quantities,
    aliases = catalog$aliases
  )
}

.bt_attach_parameter_map <- function(fit, monitor_names = NULL){

  if(inherits(fit, "error")){
    return(fit)
  }
  samples <- .extract_posterior_samples(fit, as_list = FALSE)
  columns <- colnames(samples)
  if(is.null(monitor_names)){
    monitor_names <- if(is.list(fit) && !is.null(fit$monitor)){
      fit$monitor
    }else{
      columns
    }
  }
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    monitor_names = monitor_names,
    prior_list = attr(fit, "prior_list", exact = TRUE),
    formula_design = attr(fit, "formula_design", exact = TRUE),
    formula_scale = attr(fit, "formula_scale", exact = TRUE),
    backend_anchor = attr(fit, "backend_anchor", exact = TRUE)
  )
  fit
}
