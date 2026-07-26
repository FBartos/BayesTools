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
  point_masses_declared <- NULL
  support      <- NULL

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
    if(!is.null(posterior_density[["point_masses_declared"]])){
      point_masses_declared <- isTRUE(posterior_density[["point_masses_declared"]])
    }
    if(!is.null(posterior_density[["support"]])){
      support <- .posterior_support_from_attribute(
        posterior_density[["support"]],
        exact = posterior_density[["support_exact"]]
      )
    }else if(!is.null(posterior_density[["posterior_support"]])){
      support <- .posterior_support_from_attribute(
        posterior_density[["posterior_support"]],
        exact = posterior_density[["support_exact"]]
      )
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
    if(is.null(point_masses) && !is.null(source[["point_masses"]])){
      point_masses <- source[["point_masses"]]
    }
    if(is.null(support)){
      if(!is.null(source[["support"]])){
        support <- .posterior_support_from_attribute(
          source[["support"]],
          exact = source[["support_exact"]]
        )
      }else if(!is.null(source[["posterior_support"]])){
        support <- .posterior_support_from_attribute(
          source[["posterior_support"]],
          exact = source[["support_exact"]]
        )
      }
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

  if(is.null(point_masses_declared)){
    point_masses_declared <- !is.null(point_masses)
  }
  point_masses <- .posterior_density_point_masses(point_masses)
  if(is.null(point_masses)){
    return(NULL)
  }

  return(list(
    x                     = x,
    y                     = y,
    method                = method,
    diagnostics           = diagnostics,
    support               = support,
    point_masses          = point_masses,
    point_masses_declared = point_masses_declared
  ))
}

.posterior_density_point_masses <- function(point_masses){

  empty <- data.frame(x = numeric(), mass = numeric())
  if(is.null(point_masses)){
    return(empty)
  }

  if(is.data.frame(point_masses)){
    x_name <- if("x" %in% colnames(point_masses)){
      "x"
    }else if("location" %in% colnames(point_masses)){
      "location"
    }else{
      NULL
    }
    mass_name <- if("mass" %in% colnames(point_masses)){
      "mass"
    }else if("p" %in% colnames(point_masses)){
      "p"
    }else{
      NULL
    }
    if(is.null(x_name) || is.null(mass_name)){
      return(NULL)
    }
    x <- point_masses[[x_name]]
    mass <- point_masses[[mass_name]]
  }else if(is.list(point_masses) &&
            (!is.null(point_masses[["x"]]) || !is.null(point_masses[["location"]])) &&
            (!is.null(point_masses[["mass"]]) || !is.null(point_masses[["p"]]))){
    x <- if(!is.null(point_masses[["x"]])){
      point_masses[["x"]]
    }else{
      point_masses[["location"]]
    }
    mass <- if(!is.null(point_masses[["mass"]])){
      point_masses[["mass"]]
    }else{
      point_masses[["p"]]
    }
  }else{
    return(NULL)
  }

  if(length(x) != length(mass)){
    return(NULL)
  }
  if(length(x) == 0L){
    return(empty)
  }

  out <- data.frame(
    x    = suppressWarnings(as.numeric(x)),
    mass = suppressWarnings(as.numeric(mass))
  )
  if(any(!is.finite(out[["x"]])) ||
     any(!is.finite(out[["mass"]])) ||
     any(out[["mass"]] <= 0)){
    return(NULL)
  }

  if(nrow(out) > 0L && anyDuplicated(out[["x"]])){
    mass_by_x <- tapply(out[["mass"]], out[["x"]], sum)
    out <- data.frame(
      x    = as.numeric(names(mass_by_x)),
      mass = as.numeric(mass_by_x)
    )
    out <- out[order(out[["x"]]), , drop = FALSE]
  }
  mass_bound <- .Machine$double.eps * max(8, nrow(out))
  if(nrow(out) > 0L && sum(out[["mass"]]) > 1 + mass_bound){
    return(NULL)
  }
  rownames(out) <- NULL

  return(out)
}

.posterior_density_point_masses_declared <- function(posterior_density){

  isTRUE(posterior_density[["point_masses_declared"]])
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
