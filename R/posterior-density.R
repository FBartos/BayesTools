#' @title Posterior density method helpers
#'
#' @description Helpers for normalizing public posterior-density method
#' arguments and identifying methods that are expected to use precomputed
#' posterior density or ordinate attributes.
#'
#' @param method density method.
#' @param allowed character vector of allowed density methods.
#' @param name argument name used in error messages.
#'
#' @return \code{posterior_density_method_match()} returns a scalar character
#' value. \code{posterior_density_method_uses_precomputed()} returns a scalar
#' logical value.
#'
#' @export
posterior_density_method_match <- function(method,
                                           allowed = c("KDE", "precomputed"),
                                           name = "method"){

  check_char(method, "method", check_length = 0, allow_NA = FALSE)
  check_char(allowed, "allowed", check_length = 0, allow_NA = FALSE)
  check_char(name, "name", check_length = 1, allow_NA = FALSE)

  method <- tryCatch(match.arg(method, allowed), error = function(e) NULL)
  if(is.null(method)){
    stop(
      "The '", name, "' argument must be one of '",
      paste(allowed, collapse = "', '"),
      "'.",
      call. = FALSE
    )
  }

  return(method)
}

#' @rdname posterior_density_method_match
#' @export
posterior_density_method_uses_precomputed <- function(method){

  check_char(method, "method", check_length = 1, allow_NA = FALSE)

  return(tolower(method) %in% c("precomputed", "qcmde", "iwmde"))
}

.posterior_density_method <- function(density_method){

  return(posterior_density_method_match(
    density_method,
    allowed = c("KDE", "precomputed")
  ))
}


#' @title Posterior density and ordinate attributes
#'
#' @description Construct and inspect posterior-density and posterior-ordinate
#' attributes in the schema consumed by \code{\link{Savage_Dickey_BF}} and
#' \code{\link{hypothesis_BF}}.
#'
#' @details Additional metadata must be named and cannot replace structural
#' schema fields. Posterior ordinate values are unique so a point-null lookup
#' maps to at most one stored ordinate. When a validator is supplied to
#' \code{posterior_ordinate_supports_bf()}, it is applied to matching ordinate
#' entries rather than to the multi-ordinate container as a whole.
#'
#' @param x numeric density grid locations.
#' @param y numeric density grid heights.
#' @param method estimator label stored in the attribute.
#' @param density_method public density method label stored in the attribute.
#' @param diagnostics optional estimator diagnostics.
#' @param point_masses optional point-mass table/list with \code{x} and
#' \code{mass} entries.
#' @param ... additional named metadata fields, for example \code{parameter},
#' \code{conditional}, or \code{conditional_rule}.
#'
#' @return \code{posterior_density_attribute()} returns a list suitable for a
#' \code{posterior_density} attribute.
#'
#' @export
posterior_density_attribute <- function(x, y, method, density_method,
                                        diagnostics = NULL,
                                        point_masses = NULL, ...){

  check_real(x, "x", check_length = 0, allow_NA = FALSE)
  check_real(y, "y", check_length = 0, allow_NA = FALSE)
  check_char(method, "method", check_length = 1, allow_NA = FALSE)
  check_char(density_method, "density_method", check_length = 1,
             allow_NA = FALSE)
  if(length(x) != length(y)){
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }

  metadata <- list(...)
  metadata_names <- names(metadata)
  if(length(metadata) > 0L &&
     (is.null(metadata_names) || any(!nzchar(metadata_names)))){
    stop("Additional posterior density metadata must be named.",
         call. = FALSE)
  }
  reserved <- c(
    "status", "x", "y", "method", "density_method", "diagnostics",
    "point_masses", "density", "posterior_density", "posterior_densities",
    "densities", "estimator"
  )
  if(any(metadata_names %in% reserved)){
    stop("Additional posterior density metadata cannot replace reserved fields.",
         call. = FALSE)
  }

  out <- c(list(
    status         = "ok",
    x              = x,
    y              = y,
    method         = method,
    density_method = density_method,
    diagnostics    = diagnostics,
    point_masses   = point_masses
  ), metadata)
  class(out) <- c("BayesTools_posterior_density", "list")

  if(is.null(.posterior_density_from_attribute(out))){
    stop("Posterior density attribute must contain a valid positive density grid.",
         call. = FALSE)
  }

  return(out)
}


#' @rdname posterior_density_attribute
#' @param value numeric null-hypothesis value or values.
#' @param ordinate numeric posterior ordinate height or heights.
#'
#' @return \code{posterior_ordinate_attribute()} returns a list suitable for a
#' \code{posterior_ordinate} attribute.
#'
#' @export
posterior_ordinate_attribute <- function(value, ordinate, method,
                                         density_method,
                                         diagnostics = NULL, ...){

  check_real(value, "value", check_length = 0, allow_NA = FALSE)
  check_real(ordinate, "ordinate", check_length = 0, allow_NA = FALSE)
  check_char(method, "method", check_length = 1, allow_NA = FALSE)
  check_char(density_method, "density_method", check_length = 1,
             allow_NA = FALSE)
  if(length(value) != length(ordinate)){
    stop("'value' and 'ordinate' must have the same length.", call. = FALSE)
  }
  if(any(!is.finite(value)) || any(!is.finite(ordinate)) ||
     any(ordinate <= 0)){
    stop("Posterior ordinates must be finite and positive.", call. = FALSE)
  }
  if(length(.posterior_ordinate_duplicated_values(value)) > 0L){
    stop("Posterior ordinate values must be unique.", call. = FALSE)
  }

  metadata <- list(...)
  metadata_names <- names(metadata)
  if(length(metadata) > 0L &&
     (is.null(metadata_names) || any(!nzchar(metadata_names)))){
    stop("Additional posterior ordinate metadata must be named.",
         call. = FALSE)
  }
  reserved <- c(
    "status", "value", "ordinate", "method", "density_method",
    "diagnostics", "x", "y", "null_hypothesis", "height",
    "posterior_height", "ordinates", "posterior_ordinate",
    "posterior_ordinates", "estimator"
  )
  if(any(metadata_names %in% reserved)){
    stop("Additional posterior ordinate metadata cannot replace reserved fields.",
         call. = FALSE)
  }

  out <- c(list(
    status         = "ok",
    value          = value,
    ordinate       = ordinate,
    method         = method,
    density_method = density_method,
    diagnostics    = diagnostics
  ), metadata)
  class(out) <- c("BayesTools_posterior_ordinate", "list")

  if(!.posterior_ordinate_has_data(out)){
    stop("Posterior ordinate attribute must contain valid value and ordinate fields.",
         call. = FALSE)
  }

  return(out)
}


#' @rdname posterior_density_attribute
#' @param existing existing posterior-ordinate attribute or \code{NULL}.
#' @param ordinate posterior-ordinate attribute to append.
#'
#' @return \code{posterior_ordinate_append()} returns a posterior-ordinate
#' attribute container that may hold multiple ordinates.
#'
#' @export
posterior_ordinate_append <- function(existing, ordinate){

  if(is.null(existing)){
    if(is.null(ordinate)){
      return(NULL)
    }
    if(!.posterior_ordinate_has_data(ordinate)){
      stop("'ordinate' is not a valid posterior ordinate attribute.",
           call. = FALSE)
    }
    if(length(.posterior_ordinate_duplicated_values(
      .posterior_ordinate_value_candidates(ordinate)
    )) > 0L){
      stop("Posterior ordinate attributes cannot contain duplicate values.",
           call. = FALSE)
    }
    return(ordinate)
  }
  if(is.null(ordinate)){
    if(!.posterior_ordinate_has_data(existing)){
      stop("'existing' is not a valid posterior ordinate attribute.",
           call. = FALSE)
    }
    if(length(.posterior_ordinate_duplicated_values(
      .posterior_ordinate_value_candidates(existing)
    )) > 0L){
      stop("Posterior ordinate attributes cannot contain duplicate values.",
           call. = FALSE)
    }
    return(existing)
  }
  if(!.posterior_ordinate_has_data(existing)){
    stop("'existing' is not a valid posterior ordinate attribute.",
         call. = FALSE)
  }
  if(!.posterior_ordinate_has_data(ordinate)){
    stop("'ordinate' is not a valid posterior ordinate attribute.",
         call. = FALSE)
  }
  if(length(.posterior_ordinate_duplicated_values(c(
    .posterior_ordinate_value_candidates(existing),
    .posterior_ordinate_value_candidates(ordinate)
  ))) > 0L){
    stop("Posterior ordinate attributes cannot contain duplicate values.",
         call. = FALSE)
  }

  out <- list(
    status    = "ok",
    ordinates = c(
      .posterior_ordinate_entries(existing),
      .posterior_ordinate_entries(ordinate)
    )
  )
  class(out) <- c("BayesTools_posterior_ordinates", "list")

  return(out)
}


#' @rdname posterior_density_attribute
#' @param validator optional function for stricter method-specific BF support
#' checks. It is called only after the generic BayesTools schema is valid.
#'
#' @return \code{posterior_ordinate_supports_bf()} and
#' \code{posterior_ordinate_has_value()} return scalar logical values.
#'
#' @export
posterior_ordinate_supports_bf <- function(ordinate, validator = NULL){

  if(is.null(ordinate) || !.posterior_ordinate_has_data(ordinate)){
    return(FALSE)
  }
  if(!is.null(validator)){
    if(!is.function(validator)){
      stop("'validator' must be a function.", call. = FALSE)
    }
  }

  for(entry in .posterior_ordinate_entries(ordinate)){
    values <- .posterior_ordinate_value_candidates(entry)
    if(length(values) == 0L){
      next
    }
    supported <- vapply(
      values,
      function(value) {
        parsed <- .posterior_ordinate_from_attribute(entry, value)
        if(is.null(parsed)){
          return(FALSE)
        }
        if(is.null(validator)){
          return(TRUE)
        }
        isTRUE(validator(entry)) || isTRUE(validator(parsed))
      },
      logical(1)
    )
    if(any(supported)){
      return(TRUE)
    }
  }

  return(FALSE)
}


#' @rdname posterior_density_attribute
#' @export
posterior_ordinate_has_value <- function(ordinate, value){

  check_real(value, "value", check_length = 1, allow_NA = FALSE)
  if(!is.finite(value)){
    stop("'value' must be finite.", call. = FALSE)
  }

  return(!is.null(.posterior_ordinate_from_attribute(ordinate, value)))
}


.posterior_ordinate_entries <- function(ordinate){

  if(is.list(ordinate) &&
     !is.data.frame(ordinate[["ordinates"]]) &&
     is.list(ordinate[["ordinates"]]) &&
     is.null(ordinate[["ordinates"]][["x"]]) &&
     is.null(ordinate[["ordinates"]][["value"]]) &&
     is.null(ordinate[["ordinates"]][["null_hypothesis"]])){
    return(ordinate[["ordinates"]])
  }

  return(list(ordinate))
}


.posterior_ordinate_value_candidates <- function(ordinate){

  values <- numeric()
  if(is.null(ordinate)){
    return(values)
  }
  if(is.data.frame(ordinate)){
    value_name <- intersect(c("x", "value", "null_hypothesis"),
                            colnames(ordinate))[1L]
    if(!is.na(value_name)){
      values <- c(values, ordinate[[value_name]])
    }
    values <- suppressWarnings(as.numeric(values))
    values <- values[is.finite(values)]
    return(unique(values))
  }
  if(!is.list(ordinate)){
    return(values)
  }
  for(value_name in c("x", "value", "null_hypothesis")){
    if(!is.null(ordinate[[value_name]])){
      values <- c(values, ordinate[[value_name]])
    }
  }
  for(container_name in c("ordinate", "ordinates")){
    if(!is.null(ordinate[[container_name]]) &&
       (is.list(ordinate[[container_name]]) ||
        is.data.frame(ordinate[[container_name]]))){
      if(is.list(ordinate[[container_name]]) &&
         !is.data.frame(ordinate[[container_name]]) &&
         is.null(ordinate[[container_name]][["x"]]) &&
         is.null(ordinate[[container_name]][["value"]]) &&
         is.null(ordinate[[container_name]][["null_hypothesis"]])){
        for(entry in ordinate[[container_name]]){
          values <- c(values, .posterior_ordinate_value_candidates(entry))
        }
      }else{
        values <- c(
          values,
          .posterior_ordinate_value_candidates(ordinate[[container_name]])
        )
      }
    }
  }

  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]

  return(unique(values))
}

.posterior_ordinate_duplicated_values <- function(values){

  values <- suppressWarnings(as.numeric(values))
  values <- sort(values[is.finite(values)])
  if(length(values) <= 1L){
    return(numeric())
  }

  tolerance <- sqrt(.Machine$double.eps) *
    pmax(1, abs(values[-1L]), abs(values[-length(values)]))
  duplicated <- abs(diff(values)) <= tolerance

  return(unique(values[-1L][duplicated]))
}

.posterior_density_for_method <- function(posterior_density, density_method){

  if(identical(density_method, "precomputed")){
    return(.posterior_density_from_attribute(posterior_density))
  }

  return(NULL)
}

.posterior_density_aliases <- function(...){

  aliases <- unlist(list(...), use.names = FALSE)
  aliases <- as.character(aliases)
  aliases <- aliases[!is.na(aliases) & nzchar(aliases)]

  return(unique(aliases))
}

.posterior_density_metadata_values <- function(posterior_density, fields){

  values <- NULL
  if(is.list(posterior_density)){
    for(field in fields){
      if(!is.null(posterior_density[[field]])){
        values <- c(values, unlist(posterior_density[[field]], use.names = FALSE))
      }
    }
    for(container_name in c("density", "ordinate", "ordinates")){
      if(!is.null(posterior_density[[container_name]]) &&
         is.list(posterior_density[[container_name]])){
        for(field in fields){
          if(!is.null(posterior_density[[container_name]][[field]])){
            values <- c(values, unlist(posterior_density[[container_name]][[field]], use.names = FALSE))
          }
        }
      }
    }
  }
  if(is.data.frame(posterior_density)){
    for(field in fields){
      if(field %in% colnames(posterior_density)){
        values <- c(values, unique(posterior_density[[field]]))
      }
    }
  }

  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]

  return(unique(values))
}

.posterior_density_parameter_metadata <- function(posterior_density){

  return(.posterior_density_metadata_values(
    posterior_density,
    c("parameter", "parameters", "parameter_name", "name")
  ))
}

.posterior_density_conditional_metadata <- function(posterior_density){

  return(.posterior_density_metadata_values(
    posterior_density,
    c("conditional", "condition", "conditioned_on", "conditioning")
  ))
}

.posterior_density_conditional_rule_metadata <- function(posterior_density){

  rule <- .posterior_density_metadata_values(
    posterior_density,
    c("conditional_rule", "condition_rule", "conditioning_rule")
  )
  if(length(rule) == 0L){
    return(NULL)
  }

  return(rule[[1]])
}

.posterior_density_condition_key_metadata <- function(posterior_density){

  key <- .posterior_density_metadata_values(
    posterior_density,
    c("condition_key", "conditional_key")
  )
  if(length(key) == 0L){
    return(NULL)
  }

  key[[1]]
}

.posterior_density_normalize_condition <- function(conditional){

  conditional <- unlist(conditional, use.names = FALSE)
  conditional <- as.character(conditional)
  conditional <- conditional[!is.na(conditional) & nzchar(conditional)]

  return(unique(conditional))
}

.posterior_density_condition_matches <- function(posterior_density, conditional, conditional_rule,
                                                 condition_key = NULL){

  requested_conditional <- .posterior_density_normalize_condition(conditional)
  density_conditional   <- .posterior_density_normalize_condition(
    .posterior_density_conditional_metadata(posterior_density)
  )
  requested_key <- if(is.null(condition_key)){
    .condition_event_key(requested_conditional, conditional_rule)
  }else{
    as.character(condition_key)
  }
  density_key <- .posterior_density_condition_key_metadata(posterior_density)
  if(!is.null(density_key)){
    return(identical(as.character(density_key), requested_key))
  }

  if(length(requested_conditional) == 0L && length(density_conditional) == 0L){
    return(TRUE)
  }
  if(length(requested_conditional) == 0L || length(density_conditional) == 0L){
    return(FALSE)
  }

  density_rule <- .posterior_density_conditional_rule_metadata(posterior_density)
  if(length(requested_conditional) > 1L && is.null(density_rule)){
    return(FALSE)
  }
  if(!is.null(density_rule) && !identical(as.character(density_rule), as.character(conditional_rule))){
    return(FALSE)
  }
  if(is.null(density_rule)){
    density_rule <- conditional_rule
  }

  density_key <- .condition_event_key(density_conditional, density_rule)

  identical(requested_key, density_key)
}

.posterior_density_candidate_matches <- function(posterior_density, aliases, conditional, conditional_rule,
                                                 condition_key = NULL, allow_unlabeled = FALSE){

  if(is.null(.posterior_density_from_attribute(posterior_density))){
    return(FALSE)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_density)
  if(length(parameter_names) > 0L){
    if(!any(parameter_names %in% aliases)){
      return(FALSE)
    }
  }else if(!allow_unlabeled){
    return(FALSE)
  }

  return(.posterior_density_condition_matches(
    posterior_density,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key
  ))
}

.posterior_ordinate_candidate_matches <- function(posterior_ordinate, aliases, conditional, conditional_rule,
                                                  condition_key = NULL, allow_unlabeled = FALSE,
                                                  null_hypothesis = NULL){

  if(!.posterior_ordinate_has_data(posterior_ordinate)){
    return(FALSE)
  }
  if(!is.null(null_hypothesis) &&
     is.null(.posterior_ordinate_from_attribute(posterior_ordinate, null_hypothesis))){
    return(FALSE)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_ordinate)
  if(length(parameter_names) > 0L){
    if(!any(parameter_names %in% aliases)){
      return(FALSE)
    }
  }else if(!allow_unlabeled){
    return(FALSE)
  }

  return(.posterior_density_condition_matches(
    posterior_ordinate,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key
  ))
}

.posterior_ordinate_has_data <- function(posterior_ordinate){

  if(is.null(posterior_ordinate)){
    return(FALSE)
  }
  if(is.data.frame(posterior_ordinate)){
    return(
      any(c("x", "value", "null_hypothesis") %in% colnames(posterior_ordinate)) &&
        any(c("y", "ordinate", "height", "posterior_height") %in% colnames(posterior_ordinate))
    )
  }
  if(!is.list(posterior_ordinate)){
    return(FALSE)
  }
  if(!is.null(posterior_ordinate[["status"]]) &&
     !identical(posterior_ordinate[["status"]], "ok")){
    return(FALSE)
  }
  if((!is.null(posterior_ordinate[["x"]]) ||
      !is.null(posterior_ordinate[["value"]]) ||
      !is.null(posterior_ordinate[["null_hypothesis"]])) &&
     (!is.null(posterior_ordinate[["y"]]) ||
      !is.null(posterior_ordinate[["ordinate"]]) ||
      !is.null(posterior_ordinate[["height"]]) ||
      !is.null(posterior_ordinate[["posterior_height"]]))){
    return(TRUE)
  }
  if(!is.null(posterior_ordinate[["ordinate"]]) &&
     (is.list(posterior_ordinate[["ordinate"]]) ||
      is.data.frame(posterior_ordinate[["ordinate"]]))){
    return(.posterior_ordinate_has_data(posterior_ordinate[["ordinate"]]))
  }
  if(!is.null(posterior_ordinate[["ordinates"]]) &&
     (is.list(posterior_ordinate[["ordinates"]]) ||
      is.data.frame(posterior_ordinate[["ordinates"]]))){
    if(.posterior_ordinate_has_data(posterior_ordinate[["ordinates"]])){
      return(TRUE)
    }
    if(is.list(posterior_ordinate[["ordinates"]])){
      return(any(vapply(
        posterior_ordinate[["ordinates"]],
        .posterior_ordinate_has_data,
        logical(1)
      )))
    }
  }

  return(FALSE)
}

.posterior_density_from_sources <- function(sources, aliases, conditional = NULL, conditional_rule = "AND",
                                            condition_key = NULL, allow_unlabeled = FALSE){

  aliases <- .posterior_density_aliases(aliases)
  if(length(aliases) == 0L){
    return(NULL)
  }

  for(source in sources){
    out <- .posterior_density_from_source(
      source             = source,
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = FALSE,
      depth              = 0L
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_density_from_source <- function(source, aliases, conditional, conditional_rule,
                                           condition_key, allow_unlabeled, selected_by_name, depth){

  if(is.null(source) || depth > 4L){
    return(NULL)
  }

  if(.posterior_density_candidate_matches(
    source,
    aliases          = aliases,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key,
    allow_unlabeled  = allow_unlabeled || selected_by_name
  )){
    return(source)
  }

  if(!is.list(source)){
    return(NULL)
  }

  container_names <- c("posterior_density", "posterior_densities", "densities")
  for(container_name in container_names[container_names %in% names(source)]){
    out <- .posterior_density_from_source(
      source             = source[[container_name]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = selected_by_name,
        depth              = depth + 1L
    )
    if(!is.null(out)){
      return(out)
    }
  }

  source_names <- names(source)
  if(!is.null(source_names)){
    for(alias in aliases[aliases %in% source_names]){
      out <- .posterior_density_from_source(
        source             = source[[alias]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = TRUE,
        selected_by_name   = TRUE,
        depth              = depth + 1L
      )
      if(!is.null(out)){
        return(out)
      }
    }
  }

  for(i in seq_along(source)){
    out <- .posterior_density_from_source(
      source             = source[[i]],
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = selected_by_name,
      depth              = depth + 1L
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_density_sources <- function(...){

  objects <- list(...)
  sources <- list()
  for(object in objects){
    if(is.null(object)){
      next
    }
    object_sources <- list(
      attr(object, "posterior_density", exact = TRUE),
      attr(object, "posterior_densities", exact = TRUE)
    )
    object_sources <- object_sources[!vapply(object_sources, is.null, logical(1))]
    sources <- c(sources, object_sources)
  }

  return(sources)
}

.posterior_density_attach <- function(samples, sources, parameter, conditional = NULL, conditional_rule = "AND",
                                      condition_key = NULL, allow_unlabeled = FALSE){

  if(length(sources) == 0L){
    return(samples)
  }

  if(is.matrix(samples) || is.data.frame(samples)){
    density_list <- attr(samples, "posterior_density", exact = TRUE)
    if(is.null(density_list) || !is.list(density_list)){
      density_list <- list()
    }

    density_aliases <- list()
    sample_names <- colnames(samples)
    if(!is.null(sample_names)){
      for(sample_name in sample_names){
        density_aliases[[sample_name]] <- .posterior_density_aliases(sample_name)
      }
    }

    level_names <- attr(samples, "level_names", exact = TRUE)
    if(is.list(level_names)){
      level_names <- .factor_cell_labels(level_names)
    }
    factor_cell_names <- attr(samples, "factor_cell_names", exact = TRUE)
    level_names <- .posterior_density_aliases(level_names, factor_cell_names)
    for(level_name in level_names){
      density_aliases[[level_name]] <- .posterior_density_aliases(
        level_name,
        paste0(parameter, "[", level_name, "]")
      )
    }

    if(length(density_aliases) == 0L && ncol(samples) == 1L){
      density_aliases[[parameter]] <- .posterior_density_aliases(parameter)
    }

    for(density_name in names(density_aliases)){
      density <- .posterior_density_from_sources(
        sources            = sources,
        aliases            = density_aliases[[density_name]],
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled && length(density_aliases) == 1L
      )
      if(!is.null(density)){
        density_list[[density_name]] <- density
      }
    }

    if(length(density_list) > 0L){
      attr(samples, "posterior_density") <- density_list
    }

    return(samples)
  }

  density <- .posterior_density_from_sources(
    sources            = sources,
    aliases            = .posterior_density_aliases(parameter, attr(samples, "parameter", exact = TRUE)),
    conditional        = conditional,
    conditional_rule   = conditional_rule,
    condition_key      = condition_key,
    allow_unlabeled    = allow_unlabeled
  )
  if(!is.null(density)){
    attr(samples, "posterior_density") <- density
  }

  return(samples)
}

.posterior_ordinate_sources <- function(...){

  objects <- list(...)
  sources <- list()
  for(object in objects){
    if(is.null(object)){
      next
    }
    object_sources <- list(
      attr(object, "posterior_ordinate", exact = TRUE),
      attr(object, "posterior_ordinates", exact = TRUE)
    )
    object_sources <- object_sources[!vapply(object_sources, is.null, logical(1))]
    sources <- c(sources, object_sources)
  }

  return(sources)
}

.posterior_ordinate_from_sources <- function(sources, aliases, conditional = NULL, conditional_rule = "AND",
                                             condition_key = NULL, allow_unlabeled = FALSE,
                                             null_hypothesis = NULL){

  aliases <- .posterior_density_aliases(aliases)
  if(length(aliases) == 0L){
    return(NULL)
  }

  if(is.null(null_hypothesis)){
    matches <- list()
    for(source in sources){
      matches <- c(matches, .posterior_ordinate_collect_from_source(
        source             = source,
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = FALSE,
        depth              = 0L
      ))
    }
    matches <- .posterior_ordinate_unique_sources(matches)
    if(length(matches) == 0L){
      return(NULL)
    }
    if(length(matches) == 1L){
      return(matches[[1]])
    }

    return(list(ordinates = matches))
  }

  for(source in sources){
    out <- .posterior_ordinate_from_source(
      source             = source,
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = FALSE,
      depth              = 0L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_ordinate_collect_from_source <- function(source, aliases, conditional, conditional_rule,
                                                    condition_key, allow_unlabeled,
                                                    selected_by_name, depth){

  if(is.null(source) || depth > 4L){
    return(list())
  }

  out <- list()
  if(.posterior_ordinate_candidate_matches(
    source,
    aliases          = aliases,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key,
    allow_unlabeled  = allow_unlabeled || selected_by_name
  )){
    out <- c(out, list(source))
  }

  if(!is.list(source)){
    return(out)
  }

  container_names <- c("posterior_ordinate", "posterior_ordinates", "ordinates")
  for(container_name in container_names[container_names %in% names(source)]){
    out <- c(out, .posterior_ordinate_collect_from_source(
      source             = source[[container_name]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = selected_by_name,
        depth              = depth + 1L
    ))
  }

  source_names <- names(source)
  if(!is.null(source_names)){
    matched_names <- aliases[aliases %in% source_names]
    if(length(matched_names) > 0L){
      for(alias in matched_names){
        out <- c(out, .posterior_ordinate_collect_from_source(
          source             = source[[alias]],
          aliases            = aliases,
          conditional        = conditional,
          conditional_rule   = conditional_rule,
          condition_key      = condition_key,
          allow_unlabeled    = TRUE,
          selected_by_name   = TRUE,
          depth              = depth + 1L
        ))
      }

      return(out)
    }
  }

  for(i in seq_along(source)){
    out <- c(out, .posterior_ordinate_collect_from_source(
      source             = source[[i]],
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = selected_by_name,
      depth              = depth + 1L
    ))
  }

  return(out)
}

.posterior_ordinate_unique_sources <- function(sources){

  if(length(sources) <= 1L){
    return(sources)
  }

  keys <- vapply(sources, function(source) {
    paste(serialize(source, NULL, ascii = TRUE), collapse = "")
  }, character(1))

  return(sources[!duplicated(keys)])
}

.posterior_ordinate_from_source <- function(source, aliases, conditional, conditional_rule,
                                            condition_key, allow_unlabeled, selected_by_name, depth,
                                            null_hypothesis){

  if(is.null(source) || depth > 4L){
    return(NULL)
  }

  has_ordinate <- if(is.null(null_hypothesis)){
    !is.null(source)
  }else{
    !is.null(.posterior_ordinate_from_attribute(source, null_hypothesis))
  }
  if(has_ordinate &&
     .posterior_ordinate_candidate_matches(
       source,
     aliases          = aliases,
     conditional      = conditional,
     conditional_rule = conditional_rule,
     condition_key    = condition_key,
     allow_unlabeled  = allow_unlabeled || selected_by_name,
     null_hypothesis  = null_hypothesis
     )){
    return(source)
  }

  if(!is.list(source)){
    return(NULL)
  }

  container_names <- c("posterior_ordinate", "posterior_ordinates", "ordinates")
  for(container_name in container_names[container_names %in% names(source)]){
    out <- .posterior_ordinate_from_source(
      source             = source[[container_name]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = selected_by_name,
        depth              = depth + 1L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  source_names <- names(source)
  if(!is.null(source_names)){
    for(alias in aliases[aliases %in% source_names]){
      out <- .posterior_ordinate_from_source(
        source             = source[[alias]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = TRUE,
        selected_by_name   = TRUE,
        depth              = depth + 1L,
        null_hypothesis    = null_hypothesis
      )
      if(!is.null(out)){
        return(out)
      }
    }
  }

  for(i in seq_along(source)){
    out <- .posterior_ordinate_from_source(
      source             = source[[i]],
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = selected_by_name,
      depth              = depth + 1L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_ordinate_attach <- function(samples, sources, parameter, conditional = NULL, conditional_rule = "AND",
                                       condition_key = NULL, allow_unlabeled = FALSE){

  if(length(sources) == 0L){
    return(samples)
  }

  if(is.matrix(samples) || is.data.frame(samples)){
    ordinate_list <- attr(samples, "posterior_ordinate", exact = TRUE)
    if(is.null(ordinate_list) || !is.list(ordinate_list)){
      ordinate_list <- list()
    }

    ordinate_aliases <- list()
    sample_names <- colnames(samples)
    if(!is.null(sample_names)){
      for(sample_name in sample_names){
        ordinate_aliases[[sample_name]] <- .posterior_density_aliases(sample_name)
      }
    }

    level_names <- attr(samples, "level_names", exact = TRUE)
    if(is.list(level_names)){
      level_names <- .factor_cell_labels(level_names)
    }
    factor_cell_names <- attr(samples, "factor_cell_names", exact = TRUE)
    level_names <- .posterior_density_aliases(level_names, factor_cell_names)
    for(level_name in level_names){
      ordinate_aliases[[level_name]] <- .posterior_density_aliases(
        level_name,
        paste0(parameter, "[", level_name, "]")
      )
    }

    if(length(ordinate_aliases) == 0L && ncol(samples) == 1L){
      ordinate_aliases[[parameter]] <- .posterior_density_aliases(parameter)
    }

    for(ordinate_name in names(ordinate_aliases)){
      ordinate <- .posterior_ordinate_from_sources(
        sources            = sources,
        aliases            = ordinate_aliases[[ordinate_name]],
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled && length(ordinate_aliases) == 1L
      )
      if(!is.null(ordinate)){
        ordinate_list[[ordinate_name]] <- ordinate
      }
    }

    if(length(ordinate_list) > 0L){
      attr(samples, "posterior_ordinate") <- ordinate_list
    }

    return(samples)
  }

  ordinate <- .posterior_ordinate_from_sources(
    sources            = sources,
    aliases            = .posterior_density_aliases(parameter, attr(samples, "parameter", exact = TRUE)),
    conditional        = conditional,
    conditional_rule   = conditional_rule,
    condition_key      = condition_key,
    allow_unlabeled    = allow_unlabeled
  )
  if(!is.null(ordinate)){
    attr(samples, "posterior_ordinate") <- ordinate
  }

  return(samples)
}

.posterior_density_child_attributes <- function(samples){

  out <- lapply(samples, function(x) {
    density <- attr(x, "posterior_density", exact = TRUE)
    if(is.null(.posterior_density_from_attribute(density))){
      return(NULL)
    }
    if(!.posterior_density_condition_matches(
      density,
      conditional      = attr(x, "conditional", exact = TRUE),
      conditional_rule = attr(x, "conditional_rule", exact = TRUE),
      condition_key    = attr(x, "condition_key", exact = TRUE)
    )){
      return(NULL)
    }

    return(density)
  })
  top_sources <- list(
    attr(samples, "posterior_density", exact = TRUE),
    attr(samples, "posterior_densities", exact = TRUE)
  )
  top_sources <- top_sources[!vapply(top_sources, is.null, logical(1))]
  if(length(top_sources) > 0L){
    sample_names <- names(samples)
    for(i in seq_along(samples)){
      aliases <- .posterior_density_aliases(
        sample_names[[i]],
        attr(samples[[i]], "parameter", exact = TRUE),
        attr(samples[[i]], "level_name", exact = TRUE)
      )
      conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
      conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
      condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
      density <- .posterior_density_from_sources(
        sources          = top_sources,
        aliases          = aliases,
        conditional      = conditional,
        conditional_rule = conditional_rule,
        condition_key    = condition_key,
        allow_unlabeled  = FALSE
      )
      if(is.null(out[[i]]) && !is.null(density)){
        out[[i]] <- density
      }
    }
  }
  for(top_level in top_sources){
    if(is.list(top_level) && is.null(names(top_level)) && length(top_level) == length(samples)){
      missing <- vapply(out, is.null, logical(1))
      if(any(missing)){
        sample_names <- names(samples)
        for(i in which(missing)){
          aliases <- .posterior_density_aliases(
            if(!is.null(sample_names)) sample_names[[i]] else NULL,
            attr(samples[[i]], "parameter", exact = TRUE),
            attr(samples[[i]], "level_name", exact = TRUE)
          )
          conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
          conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
          condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
          if(.posterior_density_candidate_matches(
            top_level[[i]],
            aliases          = aliases,
            conditional      = conditional,
            conditional_rule = conditional_rule,
            condition_key    = condition_key,
            allow_unlabeled  = TRUE
          )){
            out[[i]] <- top_level[[i]]
          }
        }
        names(out) <- sample_names
      }
    }
  }

  return(out)
}

.posterior_ordinate_for_method <- function(posterior_ordinate, null_hypothesis,
                                           density_method){

  if(identical(density_method, "precomputed")){
    return(.posterior_ordinate_from_attribute(posterior_ordinate, null_hypothesis))
  }

  return(NULL)
}

.posterior_ordinate_child_attributes <- function(samples, null_hypothesis = NULL){

  out <- lapply(samples, function(x) {
    ordinate <- attr(x, "posterior_ordinate", exact = TRUE)
    if(is.null(null_hypothesis)){
      if(!.posterior_ordinate_has_data(ordinate)){
        return(NULL)
      }
    }else if(is.null(.posterior_ordinate_from_attribute(ordinate, null_hypothesis))){
      return(NULL)
    }
    if(!.posterior_density_condition_matches(
      ordinate,
      conditional      = attr(x, "conditional", exact = TRUE),
      conditional_rule = attr(x, "conditional_rule", exact = TRUE),
      condition_key    = attr(x, "condition_key", exact = TRUE)
    )){
      return(NULL)
    }

    return(ordinate)
  })
  top_sources <- list(
    attr(samples, "posterior_ordinate", exact = TRUE),
    attr(samples, "posterior_ordinates", exact = TRUE)
  )
  top_sources <- top_sources[!vapply(top_sources, is.null, logical(1))]
  if(length(top_sources) > 0L){
    sample_names <- names(samples)
    for(i in seq_along(samples)){
      aliases <- .posterior_density_aliases(
        sample_names[[i]],
        attr(samples[[i]], "parameter", exact = TRUE),
        attr(samples[[i]], "level_name", exact = TRUE)
      )
      conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
      conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
      condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
      ordinate <- .posterior_ordinate_from_sources(
        sources          = top_sources,
        aliases          = aliases,
        conditional      = conditional,
        conditional_rule = conditional_rule,
        condition_key    = condition_key,
        allow_unlabeled  = FALSE,
        null_hypothesis  = null_hypothesis
      )
      if(is.null(out[[i]]) && !is.null(ordinate)){
        out[[i]] <- ordinate
      }
    }
  }
  for(top_level in top_sources){
    if(is.list(top_level) && is.null(names(top_level)) && length(top_level) == length(samples)){
      missing <- vapply(out, is.null, logical(1))
      if(any(missing)){
        sample_names <- names(samples)
        for(i in which(missing)){
          aliases <- .posterior_density_aliases(
            if(!is.null(sample_names)) sample_names[[i]] else NULL,
            attr(samples[[i]], "parameter", exact = TRUE),
            attr(samples[[i]], "level_name", exact = TRUE)
          )
          conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
          conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
          condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
          if(.posterior_ordinate_candidate_matches(
            top_level[[i]],
            aliases          = aliases,
            conditional      = conditional,
            conditional_rule = conditional_rule,
            condition_key    = condition_key,
            allow_unlabeled  = TRUE,
            null_hypothesis  = null_hypothesis
          )){
            out[[i]] <- top_level[[i]]
          }
        }
        names(out) <- sample_names
      }
    }
  }

  return(out)
}

.posterior_ordinate_from_attribute <- function(posterior_ordinate,
                                               null_hypothesis){

  if(is.null(posterior_ordinate)){
    return(NULL)
  }

  if(is.list(posterior_ordinate) &&
     !is.null(posterior_ordinate[["status"]]) &&
     !identical(posterior_ordinate[["status"]], "ok")){
    return(NULL)
  }
  if(is.list(posterior_ordinate) &&
     !is.data.frame(posterior_ordinate[["ordinates"]]) &&
     is.list(posterior_ordinate[["ordinates"]]) &&
     is.null(posterior_ordinate[["ordinates"]][["x"]]) &&
     is.null(posterior_ordinate[["ordinates"]][["value"]]) &&
     is.null(posterior_ordinate[["ordinates"]][["null_hypothesis"]])){
    candidates <- lapply(
      posterior_ordinate[["ordinates"]],
      .posterior_ordinate_from_attribute,
      null_hypothesis = null_hypothesis
    )
    matched <- !vapply(candidates, is.null, logical(1))
    if(sum(matched) != 1L){
      return(NULL)
    }

    out <- candidates[[which(matched)]]
    if(is.null(out[["method"]])){
      if(!is.null(posterior_ordinate[["method"]])){
        out[["method"]] <- posterior_ordinate[["method"]]
      }else if(!is.null(posterior_ordinate[["estimator"]])){
        out[["method"]] <- posterior_ordinate[["estimator"]]
      }
    }
    if(!is.null(posterior_ordinate[["diagnostics"]])){
      diagnostics <- .posterior_ordinate_subset_diagnostics(
        posterior_ordinate[["diagnostics"]],
        which(matched)
      )
      out[["diagnostics"]] <- .posterior_ordinate_merge_diagnostics(
        diagnostics,
        out[["diagnostics"]]
      )
    }

    return(out)
  }

  source      <- posterior_ordinate
  method      <- NULL
  diagnostics <- NULL

  if(is.list(posterior_ordinate)){
    if(!is.null(posterior_ordinate[["method"]])){
      method <- posterior_ordinate[["method"]]
    }
    if(!is.null(posterior_ordinate[["estimator"]])){
      method <- posterior_ordinate[["estimator"]]
    }
    if(!is.null(posterior_ordinate[["diagnostics"]])){
      diagnostics <- posterior_ordinate[["diagnostics"]]
    }
    if(!is.null(posterior_ordinate[["ordinate"]]) &&
       (is.list(posterior_ordinate[["ordinate"]]) ||
        is.data.frame(posterior_ordinate[["ordinate"]]))){
      source <- posterior_ordinate[["ordinate"]]
    }
    if(!is.null(posterior_ordinate[["ordinates"]]) &&
       (is.list(posterior_ordinate[["ordinates"]]) ||
        is.data.frame(posterior_ordinate[["ordinates"]]))){
      source <- posterior_ordinate[["ordinates"]]
    }
  }

  source_diagnostics <- NULL
  if(is.data.frame(source)){
    x_name <- intersect(c("x", "value", "null_hypothesis"), colnames(source))[1]
    y_name <- intersect(c("y", "ordinate", "height", "posterior_height"), colnames(source))[1]
    if(is.na(x_name) || is.na(y_name)){
      return(NULL)
    }
    x <- source[[x_name]]
    y <- source[[y_name]]
    source_diagnostics <- source[, setdiff(colnames(source), c(x_name, y_name)), drop = FALSE]
  }else if(is.list(source)){
    if(is.null(method) && !is.null(source[["method"]])){
      method <- source[["method"]]
    }
    if(is.null(method) && !is.null(source[["estimator"]])){
      method <- source[["estimator"]]
    }
    if(!is.null(source[["x"]])){
      x <- source[["x"]]
    }else if(!is.null(source[["value"]])){
      x <- source[["value"]]
    }else if(!is.null(source[["null_hypothesis"]])){
      x <- source[["null_hypothesis"]]
    }else{
      return(NULL)
    }

    if(!is.null(source[["y"]])){
      y <- source[["y"]]
    }else if(!is.null(source[["ordinate"]])){
      y <- source[["ordinate"]]
    }else if(!is.null(source[["height"]])){
      y <- source[["height"]]
    }else if(!is.null(source[["posterior_height"]])){
      y <- source[["posterior_height"]]
    }else{
      return(NULL)
    }
    source_diagnostics <- source[setdiff(
      names(source),
      c("x", "value", "null_hypothesis", "y", "ordinate", "height",
        "posterior_height", "method", "estimator", "parameter", "parameters",
        "parameter_name", "name", "conditional", "condition", "conditioned_on",
        "conditioning", "conditional_rule", "condition_rule",
        "conditioning_rule", "condition_key", "conditional_key")
    )]
  }else{
    return(NULL)
  }

  x <- as.numeric(x)
  y <- as.numeric(y)
  if(length(x) != length(y)){
    return(NULL)
  }

  keep <- is.finite(x) & is.finite(y) & y > 0
  x <- x[keep]
  y <- y[keep]
  diagnostics <- .posterior_ordinate_subset_keep(diagnostics, keep)
  source_diagnostics <- .posterior_ordinate_subset_keep(source_diagnostics, keep)
  if(length(x) == 0L){
    return(NULL)
  }

  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(null_hypothesis))
  index <- which(abs(x - null_hypothesis) <= tolerance)
  if(length(index) != 1L){
    return(NULL)
  }
  if(is.null(diagnostics)){
    diagnostics <- source_diagnostics
  }
  diagnostics        <- .posterior_ordinate_subset_diagnostics(diagnostics, index)
  source_diagnostics <- .posterior_ordinate_subset_diagnostics(source_diagnostics, index)
  diagnostics        <- .posterior_ordinate_merge_diagnostics(diagnostics, source_diagnostics)

  return(list(
    x           = x[index],
    y           = y[index],
    method      = method,
    diagnostics = diagnostics
  ))
}

.posterior_ordinate_subset_keep <- function(diagnostics, keep){

  if(is.null(diagnostics)){
    return(NULL)
  }
  if(is.data.frame(diagnostics)){
    return(diagnostics[keep, , drop = FALSE])
  }
  if(!is.list(diagnostics)){
    return(diagnostics)
  }

  for(name in names(diagnostics)){
    value <- diagnostics[[name]]
    if(is.atomic(value) && length(value) == length(keep)){
      diagnostics[[name]] <- value[keep]
    }
  }

  return(diagnostics)
}

.posterior_ordinate_subset_diagnostics <- function(diagnostics, index){

  if(is.null(diagnostics) || !is.list(diagnostics)){
    return(diagnostics)
  }

  if(is.data.frame(diagnostics)){
    diagnostics <- as.list(diagnostics[index, , drop = FALSE])
  }

  for(name in names(diagnostics)){
    value <- diagnostics[[name]]
    if(is.atomic(value) && length(value) > 1L && length(value) >= index){
      diagnostics[[name]] <- value[index]
    }
  }

  return(diagnostics)
}

.posterior_ordinate_merge_diagnostics <- function(diagnostics, source_diagnostics){

  if(is.null(diagnostics)){
    return(source_diagnostics)
  }
  if(is.null(source_diagnostics)){
    return(diagnostics)
  }
  if(!is.list(diagnostics)){
    return(source_diagnostics)
  }
  if(!is.list(source_diagnostics)){
    return(diagnostics)
  }

  for(name in names(source_diagnostics)){
    diagnostics[[name]] <- source_diagnostics[[name]]
  }

  return(diagnostics)
}

.posterior_density_from_attribute <- function(posterior_density){

  if(is.null(posterior_density)){
    return(NULL)
  }

  if(is.list(posterior_density) &&
     !is.null(posterior_density[["status"]]) &&
     !identical(posterior_density[["status"]], "ok")){
    return(NULL)
  }

  source       <- posterior_density
  method       <- NULL
  diagnostics  <- NULL
  point_masses <- NULL

  if(is.list(posterior_density)){
    if(!is.null(posterior_density[["method"]])){
      method <- posterior_density[["method"]]
    }
    if(!is.null(posterior_density[["estimator"]])){
      method <- posterior_density[["estimator"]]
    }
    if(!is.null(posterior_density[["diagnostics"]])){
      diagnostics <- posterior_density[["diagnostics"]]
    }
    if(!is.null(posterior_density[["point_masses"]])){
      point_masses <- posterior_density[["point_masses"]]
    }
    if(!is.null(posterior_density[["density"]]) &&
       (is.list(posterior_density[["density"]]) ||
        is.data.frame(posterior_density[["density"]]))){
      source <- posterior_density[["density"]]
    }
  }

  if(is.data.frame(source)){
    if(!all(c("x", "y") %in% colnames(source))){
      return(NULL)
    }
    x <- source[["x"]]
    y <- source[["y"]]
  }else if(is.list(source) && !is.null(source[["x"]]) && !is.null(source[["y"]])){
    x <- source[["x"]]
    y <- source[["y"]]
    if(is.null(method) && !is.null(source[["method"]])){
      method <- source[["method"]]
    }
    if(is.null(method) && !is.null(source[["estimator"]])){
      method <- source[["estimator"]]
    }
  }else{
    return(NULL)
  }

  x <- as.numeric(x)
  y <- as.numeric(y)
  if(length(x) != length(y)){
    return(NULL)
  }

  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if(length(x) < 2L || any(y < 0) || !any(y > 0)){
    return(NULL)
  }

  order_x <- order(x)
  x <- x[order_x]
  y <- y[order_x]
  if(anyDuplicated(x)){
    y_by_x <- split(y, x)
    x      <- as.numeric(names(y_by_x))
    y      <- vapply(y_by_x, mean, numeric(1))

    order_x <- order(x)
    x <- x[order_x]
    y <- y[order_x]
  }
  if(length(x) < 2L || diff(range(x)) <= 0){
    return(NULL)
  }

  return(list(
    x            = x,
    y            = y,
    method       = method,
    diagnostics  = diagnostics,
    point_masses = .posterior_density_point_masses(point_masses)
  ))
}

.posterior_density_point_masses <- function(point_masses){

  empty <- data.frame(x = numeric(), mass = numeric())
  if(is.null(point_masses)){
    return(empty)
  }

  if(is.data.frame(point_masses)){
    if(!"x" %in% colnames(point_masses)){
      if("location" %in% colnames(point_masses)){
        point_masses[["x"]] <- point_masses[["location"]]
      }else{
        return(empty)
      }
    }
    if(!"mass" %in% colnames(point_masses)){
      if("p" %in% colnames(point_masses)){
        point_masses[["mass"]] <- point_masses[["p"]]
      }else{
        return(empty)
      }
    }
    out <- point_masses[, c("x", "mass"), drop = FALSE]
  }else if(is.list(point_masses) &&
           (!is.null(point_masses[["x"]]) || !is.null(point_masses[["location"]])) &&
           (!is.null(point_masses[["mass"]]) || !is.null(point_masses[["p"]]))){
    out <- data.frame(
      x    = if(!is.null(point_masses[["x"]])) point_masses[["x"]] else point_masses[["location"]],
      mass = if(!is.null(point_masses[["mass"]])) point_masses[["mass"]] else point_masses[["p"]]
    )
  }else{
    return(empty)
  }

  out[["x"]]    <- as.numeric(out[["x"]])
  out[["mass"]] <- as.numeric(out[["mass"]])
  keep <- is.finite(out[["x"]]) & is.finite(out[["mass"]]) & out[["mass"]] > 0
  out <- out[keep, , drop = FALSE]
  out <- out[out[["mass"]] <= 1, , drop = FALSE]
  if(nrow(out) > 0L && anyDuplicated(out[["x"]])){
    mass_by_x <- tapply(out[["mass"]], out[["x"]], sum)
    out <- data.frame(
      x    = as.numeric(names(mass_by_x)),
      mass = as.numeric(mass_by_x)
    )
    out <- out[order(out[["x"]]), , drop = FALSE]
  }
  if(nrow(out) > 0L && sum(out[["mass"]]) > 1 + sqrt(.Machine$double.eps)){
    return(empty)
  }
  rownames(out) <- NULL

  return(out)
}

.posterior_density_height <- function(posterior_density, null_hypothesis){

  posterior_density <- .posterior_density_from_attribute(posterior_density)
  if(is.null(posterior_density)){
    return(NA_real_)
  }

  x <- posterior_density[["x"]]
  y <- posterior_density[["y"]]
  if(null_hypothesis < min(x) || null_hypothesis > max(x)){
    return(0)
  }

  height <- stats::approx(
    x    = x,
    y    = y,
    xout = null_hypothesis,
    rule = 1,
    ties = mean
  )[["y"]]
  if(!is.finite(height)){
    return(0)
  }

  return(max(0, height))
}

.posterior_density_bf_error_percent <- function(posterior_density, null_hypothesis = NULL){

  posterior_density <- .posterior_density_from_attribute(posterior_density)
  if(is.null(posterior_density) || is.null(posterior_density[["diagnostics"]])){
    return(NA_real_)
  }

  diagnostics <- posterior_density[["diagnostics"]]
  if(!is.list(diagnostics)){
    return(NA_real_)
  }

  index <- .posterior_density_diagnostics_null_index(diagnostics, null_hypothesis)
  if(is.na(index)){
    return(NA_real_)
  }
  diagnostics <- .posterior_ordinate_subset_diagnostics(diagnostics, index)

  if(!is.null(diagnostics[["BF_error_percent"]])){
    BF_error_percent <- as.numeric(diagnostics[["BF_error_percent"]])[1]
    if(is.finite(BF_error_percent) && BF_error_percent >= 0){
      return(BF_error_percent)
    }
  }
  if(!is.null(diagnostics[["bf_error_percent"]])){
    BF_error_percent <- as.numeric(diagnostics[["bf_error_percent"]])[1]
    if(is.finite(BF_error_percent) && BF_error_percent >= 0){
      return(BF_error_percent)
    }
  }

  if(is.null(diagnostics[["bf_relative_mcse"]])){
    return(NA_real_)
  }

  relative_mcse <- as.numeric(diagnostics[["bf_relative_mcse"]])[1]
  if(!is.finite(relative_mcse) || relative_mcse < 0){
    return(NA_real_)
  }

  return(100 * relative_mcse)
}

.posterior_density_diagnostics_null_index <- function(diagnostics, null_hypothesis){

  if(is.null(null_hypothesis)){
    return(NA_integer_)
  }
  if(!is.list(diagnostics)){
    return(NA_integer_)
  }

  diagnostic_list <- diagnostics
  if(is.data.frame(diagnostics)){
    diagnostic_list <- as.list(diagnostics)
  }

  for(name in c("bf_value", "value", "null_hypothesis")){
    if(is.null(diagnostic_list[[name]])){
      next
    }
    values <- suppressWarnings(as.numeric(diagnostic_list[[name]]))
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(null_hypothesis))
    index <- which(is.finite(values) & abs(values - null_hypothesis) <= tolerance)
    if(length(index) != 1L){
      return(NA_integer_)
    }

    return(index)
  }

  return(NA_integer_)
}

.posterior_ordinate_bf_error_percent <- function(posterior_ordinate){

  if(is.null(posterior_ordinate) || is.null(posterior_ordinate[["diagnostics"]])){
    return(NA_real_)
  }

  diagnostics <- posterior_ordinate[["diagnostics"]]
  if(!is.list(diagnostics)){
    return(NA_real_)
  }

  for(name in c("BF_error_percent", "bf_error_percent")){
    if(!is.null(diagnostics[[name]])){
      BF_error_percent <- as.numeric(diagnostics[[name]])[1]
      if(is.finite(BF_error_percent) && BF_error_percent >= 0){
        return(BF_error_percent)
      }
    }
  }

  for(name in c("relative_mcse", "bf_relative_mcse")){
    if(!is.null(diagnostics[[name]])){
      relative_mcse <- as.numeric(diagnostics[[name]])[1]
      if(is.finite(relative_mcse) && relative_mcse >= 0){
        return(100 * relative_mcse)
      }
    }
  }

  return(NA_real_)
}
