# Posterior density, ordinate, and support metadata are accepted only as the
# classed objects created by posterior_density_attribute(),
# posterior_ordinate_attribute() / posterior_ordinate_append(), and
# posterior_support_attribute(). An unclassed list whose elements are lists
# (such objects or nested lists of them) is a container keyed by parameter,
# level, or position; any other object is rejected.
.posterior_metadata_is_container <- function(x){

  is.list(x) && !is.object(x) &&
    all(vapply(x, function(element) is.null(element) || is.list(element), logical(1)))
}

.posterior_density_kind <- function(posterior_density){

  if(is.null(posterior_density)){
    return("null")
  }
  if(inherits(posterior_density, "BayesTools_posterior_density")){
    return("density")
  }
  if(.posterior_metadata_is_container(posterior_density)){
    return("container")
  }

  stop(
    "Posterior density metadata must be created with 'posterior_density_attribute()'.",
    call. = FALSE
  )
}

.posterior_ordinate_kind <- function(posterior_ordinate){

  if(is.null(posterior_ordinate)){
    return("null")
  }
  if(inherits(posterior_ordinate, "BayesTools_posterior_ordinate")){
    return("ordinate")
  }
  if(inherits(posterior_ordinate, "BayesTools_posterior_ordinates")){
    return("ordinates")
  }
  if(.posterior_metadata_is_container(posterior_ordinate)){
    return("container")
  }

  stop(
    "Posterior ordinate metadata must be created with ",
    "'posterior_ordinate_attribute()' or 'posterior_ordinate_append()'.",
    call. = FALSE
  )
}

.posterior_density_from_attribute <- function(posterior_density){

  if(!identical(.posterior_density_kind(posterior_density), "density")){
    return(NULL)
  }

  out <- .posterior_density_normalize(posterior_density)
  if(is.character(out)){
    stop("Posterior density metadata is invalid: ", out, ".", call. = FALSE)
  }

  return(out)
}

# Returns the density attribute with a sorted grid (heights at repeated
# locations averaged) and validated support, or a character reason when the
# attribute is invalid.
.posterior_density_normalize <- function(posterior_density){

  if(!identical(posterior_density[["status"]], "ok")){
    return("its status is not 'ok'")
  }

  x <- posterior_density[["x"]]
  y <- posterior_density[["y"]]
  if(!is.numeric(x) || !is.numeric(y) || length(x) != length(y)){
    return("'x' and 'y' must be numeric vectors of the same length")
  }
  if(any(!is.finite(x)) || any(!is.finite(y))){
    return("density grid values must be finite")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  if(length(x) < 2L || any(y < 0) || !any(y > 0)){
    return("the density grid needs at least two points, non-negative heights, and a positive height")
  }

  order_x <- order(x)
  x <- x[order_x]
  y <- y[order_x]
  if(anyDuplicated(x)){
    # group by exact value (character keys would round to 15 digits)
    unique_x <- unique(x)
    index    <- match(x, unique_x)
    y        <- as.numeric(rowsum(y, index, reorder = TRUE)) / tabulate(index)
    x        <- unique_x

    order_x <- order(x)
    x <- x[order_x]
    y <- y[order_x]
  }
  if(length(x) < 2L || diff(range(x)) <= 0){
    return("the density grid needs at least two distinct locations")
  }

  posterior_density[["x"]] <- x
  posterior_density[["y"]] <- y
  posterior_density["support"] <- list(
    .posterior_support_from_attribute(posterior_density[["support"]])
  )

  return(posterior_density)
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
    index <- which(is.finite(values) & values == null_hypothesis)
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
