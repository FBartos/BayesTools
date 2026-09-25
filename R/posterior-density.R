#' @title Posterior density method helpers
#'
#' @description Helpers for normalizing public posterior-density method
#' arguments and identifying the method that uses precomputed posterior
#' density or ordinate metadata. Only \code{"precomputed"} does; packages
#' with their own estimators (e.g., qCMDE or IWMDE) map their method names to
#' \code{"precomputed"} themselves.
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

  return(identical(tolower(method), "precomputed"))
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
#' finite range of an estimator grid.
#'
#' The objects returned by these constructors are the only accepted posterior
#' density and ordinate metadata: other objects set as \code{posterior_density}
#' or \code{posterior_ordinate} metadata of posterior draws
#' ([posterior_metadata()]) are rejected with an error. An unclassed list of
#' such objects, named by parameter or level, can hold the metadata of several
#' parameters or levels. The fields used to match the metadata to posterior
#' samples are \code{parameter}, \code{conditional}, \code{conditional_rule},
#' and \code{condition_key}.
#'
#' A stored density describes only the continuous part of the posterior: for
#' a posterior with point masses, its heights integrate to the continuous
#' mass. The point masses are the posterior atoms of the draws, which are
#' declared only in their \code{atoms} metadata with
#' [posterior_atom_attribute()]; densities do not carry point masses.
#'
#' @param x numeric density grid locations.
#' @param y numeric density grid heights.
#' @param method estimator label stored in the attribute.
#' @param density_method public density method label stored in the attribute.
#' @param diagnostics optional estimator diagnostics.
#' @param support optional exact support metadata for the density scale,
#' created with \code{\link{posterior_support_attribute}()}. KDE boundary
#' reflection uses only interval-capable support; support with
#' \code{exact = FALSE} (plotting or integration limits rather than true
#' support boundaries) is not used to exclude a null hypothesis.
#' @param ... additional named metadata fields, for example \code{parameter},
#' \code{conditional}, or \code{conditional_rule}.
#'
#' @return \code{posterior_density_attribute()} returns a list suitable for the
#' \code{posterior_density} metadata of posterior draws
#' ([posterior_metadata()]).
#'
#' @export
posterior_density_attribute <- function(x, y, method, density_method,
                                        diagnostics = NULL,
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
  if(any(metadata_names %in% c("point_masses", "point_masses_declared"))){
    stop(
      "Posterior densities do not carry 'point_masses': declare the point ",
      "masses of the posterior as the 'atoms' metadata of its draws with ",
      "'posterior_atom_attribute()'.",
      call. = FALSE
    )
  }
  reserved <- c(
    "status", "x", "y", "method", "density_method", "diagnostics",
    "support", "posterior_support", "density",
    "posterior_density", "posterior_densities", "densities", "estimator"
  )
  if(any(metadata_names %in% reserved)){
    stop("Additional posterior density metadata cannot replace reserved fields.",
         call. = FALSE)
  }

  if(!is.null(support)){
    if(!inherits(support, "BayesTools_posterior_support")){
      stop("'support' must be created with 'posterior_support_attribute()'.",
           call. = FALSE)
    }
    support <- .posterior_support_from_attribute(support)
  }

  out <- c(list(
    status         = "ok",
    x              = x,
    y              = y,
    method         = method,
    density_method = density_method,
    diagnostics    = diagnostics,
    support        = support
  ), metadata)
  class(out) <- c("BayesTools_posterior_density", "list")

  if(is.character(.posterior_density_normalize(out))){
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
#' @return \code{posterior_ordinate_attribute()} returns a list suitable for the
#' \code{posterior_ordinate} metadata of posterior draws
#' ([posterior_metadata()]).
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
    if(!.posterior_ordinate_is_attribute(ordinate)){
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
    if(!.posterior_ordinate_is_attribute(existing)){
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
  if(!.posterior_ordinate_is_attribute(existing)){
    stop("'existing' is not a valid posterior ordinate attribute.",
         call. = FALSE)
  }
  if(!.posterior_ordinate_is_attribute(ordinate)){
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


# The single-ordinate attributes of a posterior-ordinate attribute.
.posterior_ordinate_entries <- function(ordinate){

  kind <- .posterior_ordinate_kind(ordinate)
  if(identical(kind, "ordinate")){
    return(list(ordinate))
  }
  if(!identical(kind, "ordinates")){
    return(list())
  }

  entries <- ordinate[["ordinates"]]
  if(!identical(ordinate[["status"]], "ok") || !is.list(entries) ||
     is.object(entries) ||
     !all(vapply(entries, inherits, logical(1), what = "BayesTools_posterior_ordinate"))){
    stop(
      "Posterior ordinate metadata is invalid: a multi-ordinate attribute needs ",
      "status 'ok' and 'ordinates' created with 'posterior_ordinate_attribute()'.",
      call. = FALSE
    )
  }

  return(entries)
}


.posterior_ordinate_is_attribute <- function(ordinate){

  inherits(ordinate, c(
    "BayesTools_posterior_ordinate",
    "BayesTools_posterior_ordinates"
  )) && .posterior_ordinate_has_data(ordinate)
}


.posterior_ordinate_value_candidates <- function(ordinate){

  values <- unlist(lapply(
    .posterior_ordinate_entries(ordinate),
    function(entry) .posterior_ordinate_values(entry)[["value"]]
  ), use.names = FALSE)

  return(unique(as.numeric(values)))
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
