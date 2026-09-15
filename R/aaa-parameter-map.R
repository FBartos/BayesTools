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
  cache <- .bt_parameter_map_cache(map)
  if(!.bt_parameter_map_cache_matches(map, cache)){
    .bt_validate_parameter_map(map)
    # The map this cache slot described has been replaced, so every derived
    # entry in it is stale, including entries this package did not create.
    .bt_parameter_map_cache_store(map, cache)
  }
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
  attr(out, "runtime_cache_id") <- .bt_parameter_map_cache_new_id()
  .bt_parameter_map_cache_store(out, .bt_parameter_map_cache(out))
  out
}


# Parameter-map runtime cache
#
# The cache lives in a session-local registry keyed by an id carried on the
# map, never in the map itself. A cache attached to the object would be
# serialized with every saved fit: it would inflate the file, and - worse -
# a fit reloaded months later would replay entries derived from whatever the
# packages looked like when it was saved. Keying by id means a reloaded fit
# simply misses and recomputes.
#
# The registry is bounded, because entries hold references to the map tables
# (so the validity check stays an O(1) pointer comparison) and to whatever
# consumers store.
.bt_parameter_map_cache_limit <- function(){
  64L
}

# Ids must never draw from the \R RNG: maps are built inside seeded fitting
# code, and consuming the stream there would move seeded results. A per-session
# stamp plus a counter is unique within a session and across sessions, so a
# reloaded fit cannot land on a live slot belonging to a different map.
.bt_parameter_map_cache_new_id <- function(){

  stamp <- .BayesTools_private$parameter_map_cache_session
  if(is.null(stamp)){
    stamp <- paste0(
      Sys.getpid(), "-",
      format(as.numeric(Sys.time()), digits = 15, scientific = FALSE)
    )
    .BayesTools_private$parameter_map_cache_session <- stamp
  }
  sequence <- .BayesTools_private$parameter_map_cache_sequence
  sequence <- if(is.null(sequence)) 1L else sequence + 1L
  .BayesTools_private$parameter_map_cache_sequence <- sequence

  paste0("map-", stamp, "-", sequence)
}

.bt_parameter_map_cache_registry <- function(){

  registry <- .BayesTools_private$parameter_map_cache
  if(!is.environment(registry)){
    registry <- new.env(parent = emptyenv())
    .BayesTools_private$parameter_map_cache <- registry
    .BayesTools_private$parameter_map_cache_order <- character()
  }
  registry
}

.bt_parameter_map_cache <- function(map){

  id <- attr(map, "runtime_cache_id", exact = TRUE)
  if(!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)){
    return(NULL)
  }
  registry <- .bt_parameter_map_cache_registry()

  # A hit does no bookkeeping at all. Reordering for recency would allocate a
  # character vector on every accessor call, and this is the hot path the
  # cache exists to keep cheap; with a bound of 64 slots and a handful of live
  # maps, insertion order evicts just as well.
  existing <- registry[[id]]
  if(!is.null(existing)){
    return(existing)
  }

  order <- .BayesTools_private$parameter_map_cache_order
  cache <- new.env(parent = emptyenv())
  cache$providers <- new.env(parent = emptyenv())
  registry[[id]] <- cache
  order <- c(order, id)
  if(length(order) > .bt_parameter_map_cache_limit()){
    evicted <- order[seq_len(length(order) - .bt_parameter_map_cache_limit())]
    rm(list = evicted, envir = registry)
    order <- setdiff(order, evicted)
  }
  .BayesTools_private$parameter_map_cache_order <- order

  cache
}

.bt_parameter_map_cache_matches <- function(map, cache){

  is.environment(cache) &&
    isTRUE(cache$validated) &&
    identical(cache$coordinates, map$coordinates) &&
    identical(cache$quantities, map$quantities) &&
    identical(cache$aliases, map$aliases)
}

.bt_parameter_map_cache_store <- function(map, cache){

  if(!is.environment(cache)){
    return(invisible(NULL))
  }
  cache$validated <- TRUE
  cache$coordinates <- map$coordinates
  cache$quantities <- map$quantities
  cache$aliases <- map$aliases
  # Consumer entries were derived from the map that was here before, so they
  # cannot survive it. BayesTools owns this environment and clears it whole.
  cache$providers <- new.env(parent = emptyenv())
  invisible(NULL)
}


#' @title Cache a value derived from a fitted parameter map
#'
#' @description Stores one value per provider against a fitted parameter map,
#' for the lifetime of the \R session. It exists so that packages building on
#' BayesTools can avoid recomputing map-derived metadata without inventing
#' their own storage inside fitted objects.
#'
#' The cache is keyed by the map \emph{and} by `key`, a value naming everything
#' else the cached result was derived from. Whenever `key` stops being
#' [identical()] to the stored one the value is recomputed, so a result that
#' also depends on data or priors stays correct when those change. Passing a
#' `key` that does not cover every input is the one way to use this
#' incorrectly.
#'
#' Entries never travel with a saved fit, and BayesTools discards every
#' provider's entries whenever the map's own tables are replaced.
#'
#' @param map a parameter map, as returned by [parameter_map()].
#' @param provider name of the calling package.
#' @param key a value identifying every input other than `map` that `compute`
#'   depends on. Compared with [identical()].
#' @param compute a function of no arguments returning the value to cache.
#'
#' @return The cached or freshly computed value of `compute()`.
#'
#' @seealso [parameter_map()]
#' @export
parameter_map_cache <- function(map, provider, key, compute){

  if(!inherits(map, "BayesTools_parameter_map")){
    stop("'map' must be a 'BayesTools_parameter_map' object.", call. = FALSE)
  }
  check_char(provider, "provider", check_length = 1L, allow_NA = FALSE)
  if(!nzchar(provider)){
    stop("'provider' must be a non-empty package name.", call. = FALSE)
  }
  if(!is.function(compute)){
    stop("'compute' must be a function of no arguments.", call. = FALSE)
  }

  cache <- .bt_parameter_map_cache(map)
  if(!is.environment(cache) || !is.environment(cache$providers)){
    return(compute())
  }

  entry <- cache$providers[[provider]]
  if(is.list(entry) && identical(entry$key, key)){
    return(entry$value)
  }

  value <- compute()
  cache$providers[[provider]] <- list(key = key, value = value)
  value
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
