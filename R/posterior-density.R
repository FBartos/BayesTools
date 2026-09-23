#' @title Posterior density method helpers
#'
#' @description Helpers for normalizing public posterior-density method
#' arguments and identifying methods that are expected to use precomputed
#' posterior density or ordinate attributes. Estimator labels such as
#' \code{"qCMDE"} and \code{"IWMDE"} are treated as precomputed metadata
#' sources by \code{posterior_density_method_uses_precomputed()}, but they are
#' not public \code{density_method} argument values.
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
#' entries rather than to the multi-ordinate container as a whole. The
#' \code{support} field is used to check whether exact support excludes a
#' requested point-null value and, when a kernel-density fallback must be
#' computed from posterior samples, to supply boundary reflection. It must
#' describe the true posterior support on the same scale as \code{x}, not the
#' finite range of an estimator grid. The stable public metadata schema is the
#' object returned by these constructors. Raw list and data-frame attributes
#' with equivalent fields are parsed for backwards compatibility, but extension
#' code should prefer the constructors.
#'
#' @param x numeric density grid locations.
#' @param y numeric density grid heights.
#' @param method estimator label stored in the attribute.
#' @param density_method public density method label stored in the attribute.
#' @param diagnostics optional estimator diagnostics.
#' @param point_masses optional point-mass table/list with \code{x} and
#' \code{mass} entries. Aliases \code{location} for \code{x} and \code{p} for
#' \code{mass} are accepted. Locations and masses must be finite, masses must
#' be positive, and the aggregated point mass cannot exceed one.
#' @param support optional exact support metadata for the density scale. Supply
#' a trusted numeric \code{c(lower, upper)} vector, or a list with
#' \code{bounds}, optional \code{points}, optional \code{type}, and optional
#' \code{exact} entries. The \code{type} entry may be \code{"interval"},
#' \code{"points"}, or \code{"mixed"}; KDE boundary reflection uses only
#' interval-capable support. Set \code{exact = FALSE} when the values are
#' plotting or integration limits rather than true support boundaries.
#' @param ... additional named metadata fields, for example \code{parameter},
#' \code{conditional}, or \code{conditional_rule}.
#'
#' @return \code{posterior_density_attribute()} returns a list suitable for a
#' \code{posterior_density} attribute.
#'
#' @export
posterior_density_attribute <- function(x, y, method, density_method,
                                        diagnostics = NULL,
                                        point_masses = NULL,
                                        support = NULL, ...){

  check_real(x, "x", check_length = 0, allow_NA = FALSE)
  check_real(y, "y", check_length = 0, allow_NA = FALSE)
  check_char(method, "method", check_length = 1, allow_NA = FALSE)
  check_char(density_method, "density_method", check_length = 1,
             allow_NA = FALSE)
  if(length(x) != length(y)){
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  if(any(!is.finite(x)) || any(!is.finite(y))){
    stop("Posterior density grid values must be finite.", call. = FALSE)
  }

  metadata <- list(...)
  metadata_names <- names(metadata)
  if(length(metadata) > 0L &&
     (is.null(metadata_names) || anyNA(metadata_names) ||
      any(!nzchar(metadata_names)) || anyDuplicated(metadata_names))){
    stop("Additional posterior density metadata must have unique, nonmissing names.",
         call. = FALSE)
  }
  reserved <- c(
    "status", "x", "y", "method", "density_method", "diagnostics",
    "point_masses", "point_masses_declared", "support",
    "posterior_support", "density",
    "posterior_density", "posterior_densities", "densities", "estimator"
  )
  if(any(metadata_names %in% reserved)){
    stop("Additional posterior density metadata cannot replace reserved fields.",
         call. = FALSE)
  }

  if(!is.null(support)){
    support <- .posterior_support_from_attribute(support)
    if(is.null(support)){
      stop("Posterior density support metadata is invalid.", call. = FALSE)
    }
  }

  point_masses_declared <- !is.null(point_masses)
  point_masses <- .posterior_density_point_masses(point_masses)
  if(is.null(point_masses)){
    stop("Posterior density 'point_masses' metadata is invalid.",
         call. = FALSE)
  }

  out <- c(list(
    status         = "ok",
    x              = x,
    y              = y,
    method         = method,
    density_method = density_method,
    diagnostics    = diagnostics,
    point_masses   = point_masses,
    point_masses_declared = point_masses_declared,
    support        = support
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
#' @param ordinate for \code{posterior_ordinate_attribute()}, numeric
#' posterior ordinate height or heights; for
#' \code{posterior_ordinate_append()}, a posterior-ordinate attribute to
#' append.
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
     (is.null(metadata_names) || anyNA(metadata_names) ||
      any(!nzchar(metadata_names)) || anyDuplicated(metadata_names))){
    stop("Additional posterior ordinate metadata must have unique, nonmissing names.",
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

  duplicated <- diff(values) == 0

  return(unique(values[-1L][duplicated]))
}
