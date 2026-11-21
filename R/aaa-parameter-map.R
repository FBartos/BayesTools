# Authoritative fitted parameter map.

.bt_parameter_map_version <- 23L

#' Fitted parameter map and coordinate view
#'
#' @description
#' `parameter_map()` returns the single versioned parameter map stored by
#' [JAGS_fit()]. It contains concrete fitted `coordinates`, public semantic `quantities`,
#' exact `aliases` and directed [parameter_relationships()].
#' A complete pre-R2 version 22 map is read as version 23. Its saved coordinates,
#' quantities and aliases are preserved, and its relationships are derived from
#' the fitted model declarations.
#'
#' `parameter_coordinates()` returns the coordinate table from that map. Each
#' row describes one sampled, structural, or unavailable fitted coordinate.
#' Coordinate names are backend identities, not additional public aliases.
#'
#' `parameter_map_schema()` documents all four tables, while
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
    .bt_stop_refit_required(
      "The fitted object does not contain parameter-map metadata. ",
      "Refit the model with the current BayesTools version."
    )
  }
  contract <- attr(object, "fit_contract", exact = TRUE)
  if(.bt_parameter_map_v22_pair(map, contract)){
    return(.bt_parameter_map_from_v22(object, map, contract))
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
      "logical",
      "character"
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
      "Whether the coordinate is an implementation-level dependency.",
      paste(
        "Convergence role set at fit time: sampled, indicator, structural,",
        "derived, or auxiliary (see JAGS_check_convergence())."
      )
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
    aliases = catalog$aliases,
    relationships = catalog$relationships
  )
}

# Schema 22 stored the same three parameter tables. Its pre-R2 bindings have
# all the information needed to add relationships. Keep the saved tables and
# validate that the declarations still describe their native quantities.
.bt_parameter_map_v22_pair <- function(map, contract){

  (is.list(map) && identical(map[["schema_version"]], 22L)) ||
    (is.list(contract) && identical(contract[["parameter_map_version"]], 22L))
}

.bt_parameter_map_from_v22 <- function(fit, map, contract){

  priors <- attr(fit, "prior_list", exact = TRUE)
  designs <- attr(fit, "formula_design", exact = TRUE)
  scales <- attr(fit, "formula_scale", exact = TRUE)

  .bt_content_memo("parameter_map_v22", list(map, contract, priors, designs, scales), function(){
    refuse <- function(){
      .bt_stop_refit_required("The saved parameter map or its fitted declarations are incomplete or incompatible. ",
                              "Refit the model with the current BayesTools version.")
    }
    .bt_validate_fit_contract_object(contract)
    # Only a matching pair of version-22 declarations can be upgraded.
    if(!is.list(map) || !identical(map[["schema_version"]], 22L) ||
       !identical(contract[["parameter_map_version"]], 22L)){
      .bt_unsupported_parameter_map()
    }
    supported <- .bt_supported_fit_contract()
    supported$parameter_map <- 22L
    if(!identical(unname(contract[-1L]), unname(supported)) ||
       !inherits(map, "BayesTools_parameter_map") ||
       !identical(names(map), c("schema_version", "coordinates", "quantities", "aliases"))) refuse()

    tryCatch({
      .bt_validate_parameter_coordinates(map$coordinates)
      .bt_validate_parameter_catalog_tables(map$quantities, map$aliases)
      native <- map$quantities[map$quantities$provider == "BayesTools", , drop = FALSE]
      has_R2 <- function(p){
        is.R2_budget(p) || is.prior.R2(p) || is.prior.R2_budget(p) ||
          !is.null(attr(p, "R2_budget", exact = TRUE)) ||
          (is.prior.mixture(p) && any(vapply(p, has_R2, logical(1))))
      }
      if(any(native$owner_type == "r2_budget") ||
         any(vapply(native$extraction_key, function(key) key$type %in% c("R2_share", "R2_open"), logical(1))) ||
         any(vapply(priors, has_R2, logical(1)))) refuse()
      if(length(priors) && (!is.list(priors) || is.null(names(priors)) || anyNA(names(priors)) ||
         any(!nzchar(names(priors))) || anyDuplicated(names(priors)) || !all(vapply(priors, is.prior, logical(1))))) refuse()

      formulas <- unique(c(native$formula_parameter, map$coordinates$formula_parameter))
      if(!is.null(designs) && !is.list(designs)) refuse()
      if(!all(formulas[nzchar(formulas)] %in% names(designs)) || anyDuplicated(names(designs)) ||
         (length(designs) && (is.null(names(designs)) || anyNA(names(designs)) || any(!nzchar(names(designs)))))) refuse()
      for(parameter in names(designs)){
        design <- designs[[parameter]]
        if(!identical(design$parameter, parameter) ||
           !is.list(design$prior_list) || !all(names(design$prior_list) %in% names(priors)) ||
           any(vapply(design$prior_list, has_R2, logical(1))) ||
           !is.list(design$random_effects) || !is.list(design$random_allocations) ||
           !identical(design$formula_scale, scales[[parameter]])) refuse()
        .bt_validate_formula_design_replay_schema(design)
        for(block in design$random_effects) .bt_check_random_sd_binding(block$sd_binding)
      }

      catalog <- .bt_build_parameter_catalog(map$coordinates, priors, designs, scales)
      expected <- catalog$quantities[match(native$quantity_id, catalog$quantities$quantity_id), , drop = FALSE]
      rownames(native) <- rownames(expected) <- NULL
      allocation_ids <- native$quantity_id[native$owner_type == "variance_allocation"]
      reported <- catalog$relationships$consumer_id[catalog$relationships$consumer_type == "quantity"]
      if(nrow(native) != nrow(catalog$quantities) || !identical(native, expected) ||
         !all(allocation_ids %in% reported)) refuse()
      .bt_parameter_map_new(map$coordinates, map$quantities, map$aliases, catalog$relationships)
    }, BayesTools_refit_required = function(e) stop(e), error = function(e) refuse())
  })
}

.bt_fit_with_parameter_map <- function(fit, map){

  contract <- attr(fit, "fit_contract", exact = TRUE)
  if(identical(contract[["parameter_map_version"]], 22L)){
    contract$parameter_map_version <- map$schema_version
    attr(fit, "fit_contract") <- contract
  }
  attr(fit, "parameter_map") <- map
  fit
}

.bt_parameter_map_new <- function(coordinates, quantities, aliases,
                                  relationships = .bt_empty_relationships()){

  out <- list(
    schema_version = .bt_parameter_map_version,
    coordinates = coordinates,
    quantities = quantities,
    aliases = aliases,
    relationships = relationships
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

  .bt_parameter_map_header_valid(map) &&
    is.environment(cache) &&
    isTRUE(cache$validated) &&
    identical(cache$coordinates, map$coordinates) &&
    identical(cache$quantities, map$quantities) &&
    identical(cache$aliases, map$aliases) &&
    identical(cache$relationships, map$relationships)
}

.bt_parameter_map_cache_store <- function(map, cache){

  if(!is.environment(cache)){
    return(invisible(NULL))
  }
  cache$validated <- TRUE
  cache$coordinates <- map$coordinates
  cache$quantities <- map$quantities
  cache$aliases <- map$aliases
  cache$relationships <- map$relationships
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
  if(!.bt_parameter_map_cache_matches(map, cache)){
    .bt_validate_parameter_map(map)
    .bt_parameter_map_cache_store(map, cache)
  }
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
    aliases = map$aliases,
    relationships = map$relationships
  )
  class(out) <- c("BayesTools_parameter_catalog", "list")
  out
}

.bt_parameter_map_header_valid <- function(map){

  inherits(map, "BayesTools_parameter_map") &&
    is.list(map) &&
    identical(
      names(map),
      c("schema_version", "coordinates", "quantities", "aliases", "relationships")
    ) &&
    identical(map$schema_version, .bt_parameter_map_version)
}

.bt_unsupported_parameter_map <- function(){

  .bt_stop_refit_required(
    "Parameter-map metadata are missing, malformed, or unsupported. ",
    "Refit the model with the current BayesTools version."
  )
}

.bt_validate_parameter_map <- function(map){

  if(!.bt_parameter_map_header_valid(map)) .bt_unsupported_parameter_map()

  .bt_validate_parameter_coordinates(map$coordinates)
  .bt_validate_parameter_catalog_tables(map$quantities, map$aliases)
  .bt_validate_relationships(map$relationships, map$quantities)

  native <- map$quantities$provider == "BayesTools"
  dependencies <- unique(unlist(
    lapply(map$quantities$extraction_key[native], `[[`, "dependencies"),
    use.names = FALSE
  ))
  missing <- setdiff(dependencies, map$coordinates$coordinate_name)
  if(length(missing) > 0L){
    .bt_stop_refit_required(
      "Parameter-map quantities reference unknown coordinate dependencies: ",
      paste0("'", missing, "'", collapse = ", "),
      ". Refit or rebuild the map with this version of BayesTools."
    )
  }

  invisible(TRUE)
}

.bt_build_parameter_map <- function(columns, monitor_names = columns,
                                    prior_list = NULL,
                                    formula_design = NULL,
                                    formula_scale = NULL,
                                    backend_anchor = NULL,
                                    add_parameters = NULL,
                                    model_syntax = NULL,
                                    data_names = NULL){

  coordinates <- .bt_build_parameter_coordinates(
    columns = columns,
    monitor_names = monitor_names,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale,
    backend_anchor = backend_anchor,
    add_parameters = add_parameters,
    model_syntax = model_syntax,
    data_names = data_names
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
    aliases = catalog$aliases,
    relationships = catalog$relationships
  )
}

# 'data_names' are the names of fully observed model data; they let the
# convergence roles treat monitored functions of data as constants.
.bt_attach_parameter_map <- function(fit, monitor_names = NULL,
                                     data_names = NULL,
                                     model_syntax = attr(fit, "model_syntax", exact = TRUE)){

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
    backend_anchor = attr(fit, "backend_anchor", exact = TRUE),
    add_parameters = attr(fit, "add_parameters", exact = TRUE),
    model_syntax = model_syntax,
    data_names = data_names
  )
  fit
}
