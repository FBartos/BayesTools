#' @title Compute Savage-Dickey inclusion Bayes factors
#'
#' @description Computes Savage-Dickey (density ratio) inclusion Bayes factors
#' based the change of height from prior to posterior distribution at the test value.
#'
#' @param posterior marginal posterior distribution generated via the
#' \code{marginal_posterior} function
#' @param null_hypothesis point null hypothesis to test. Defaults to \code{0}
#' @param normal_approximation whether the height of prior and posterior density should be
#' approximated via a normal distribution (rather than kernel density). Defaults to \code{FALSE}.
#' @param silent whether warnings should be returned silently. Defaults to \code{FALSE}
#' @param density_method density source for the posterior ordinate. \code{"KDE"}
#' computes a kernel density estimate, using boundary reflection when exact
#' posterior-support metadata is available. \code{"precomputed"} requires a
#' valid \code{posterior_ordinate} attribute when present, or otherwise a valid
#' \code{posterior_density} attribute.
#'
#' @details Marginal posterior vectors may carry a \code{posterior_ordinate}
#' attribute with exact \code{value} and \code{ordinate} entries. When
#' \code{density_method = "precomputed"} and
#' \code{normal_approximation = FALSE}, a matching ordinate is used for the
#' Savage-Dickey ratio. If no matching ordinate is available, a valid
#' \code{posterior_density} grid with \code{x} and \code{y} coordinates is
#' used. If neither valid precomputed source is available, an error is thrown.
#' Exact support metadata is also checked before accepting precomputed posterior
#' ordinates or densities; support is used to exclude the null only after it
#' is validated against the posterior samples. A stale precomputed value is
#' ignored when compatible exact support excludes the null; returned Bayes
#' factors label this path as exact support exclusion rather than as a KDE
#' estimate. Support is not inferred from prior-density grids,
#' posterior-density grids, or plotting ranges, which may be finite numerical
#' integration ranges rather than true support boundaries. Exact support
#' metadata is either a numeric
#' \code{c(lower, upper)} vector or a list with
#' \code{bounds = c(lower, upper)} and \code{exact = TRUE}; optional
#' \code{type = "points"} support is not treated as a continuous interval for
#' KDE boundary reflection. The same schema may be supplied as \code{support}
#' in a \code{posterior_density} attribute. When the prior density at the null
#' is zero or non-finite, the returned Bayes factor carries a warning because
#' the point-null Savage-Dickey ratio is not a regular density ratio at that
#' point.
#'
#' @return \code{Savage_Dickey_BF} returns a Bayes factor.
#'
#' @export
Savage_Dickey_BF <- function(posterior, null_hypothesis = 0, normal_approximation = FALSE, silent = FALSE,
                             density_method = c("KDE", "precomputed")){

  if(!inherits(posterior, "marginal_posterior"))
    stop("'BF_savage_dickey' function requires an object of class 'marginal_posteriors'")
  check_real(null_hypothesis, "null_hypothesis", allow_NA = FALSE)
  if(!is.finite(null_hypothesis)){
    stop("The 'null_hypothesis' argument must be finite.", call. = FALSE)
  }
  check_bool(normal_approximation, "normal_approximation", allow_NA = FALSE)
  check_bool(silent, "silent", allow_NA = FALSE)
  density_method <- .posterior_density_method(density_method)

  if(is.list(posterior)){
    bf <- list()
    for(i in seq_along(posterior)){
      posterior_i <- .posterior_precomputed_child(
        parent          = posterior,
        child           = posterior[[i]],
        index           = i,
        null_hypothesis = null_hypothesis,
        density_method  = density_method
      )
      bf[[i]] <- .Savage_Dickey_BF.fun(posterior_i, null_hypothesis, normal_approximation, silent, density_method)
    }
    names(bf) <- names(posterior)
  }else{
    bf <- .Savage_Dickey_BF.fun(posterior, null_hypothesis, normal_approximation, silent, density_method)
  }

  return(bf)
}

.Savage_Dickey_BF.fun    <- function(posterior, null_hypothesis, normal_approximation, silent, density_method){

  if(is.null(attr(posterior, "prior_density")))
    stop("there are no prior densities for the posterior distribution", call. = FALSE)

  prior <- attr(posterior, "prior_density")

  warnings <- NULL
  stored_posterior_density <- NULL
  stored_posterior_ordinate <- NULL
  posterior_density_source <- if(isTRUE(normal_approximation)) "normal" else "KDE"
  posterior_density_boundary_reflection <- FALSE
  posterior_density_support_bounds <- NULL
  posterior_density_fallback_warnings <- NULL
  BF_error_percent <- NA_real_
  posterior_ordinate_status <- NULL
  posterior_density_status <- NULL
  if(!normal_approximation){
    posterior_ordinate_status <- .posterior_ordinate_direct_status(
      posterior,
      null_hypothesis = null_hypothesis
    )
    posterior_density_status <- .posterior_density_direct_status(posterior)
    if(identical(density_method, "precomputed") &&
       isTRUE(posterior_ordinate_status[["valid"]])){
      stored_posterior_ordinate <- posterior_ordinate_status[["value"]]
    }
    if(identical(density_method, "precomputed") &&
       isTRUE(posterior_density_status[["valid"]])){
      stored_posterior_density <- posterior_density_status[["value"]]
    }
  }
  stored_posterior_density_support <- if(!is.null(stored_posterior_density)){
    stored_posterior_density[["support"]]
  }else{
    NULL
  }
  if(!isTRUE(normal_approximation) &&
     identical(density_method, "precomputed") &&
     is.null(stored_posterior_ordinate) &&
     is.null(stored_posterior_density)){
    if(!is.null(posterior_ordinate_status) &&
       isTRUE(posterior_ordinate_status[["present"]]) &&
       isTRUE(posterior_ordinate_status[["relevant"]]) &&
       !isTRUE(posterior_ordinate_status[["valid"]])){
      stop(
        "Precomputed posterior ordinate metadata is present but invalid ",
        "for the requested null hypothesis.",
        call. = FALSE
      )
    }
    if(!is.null(posterior_density_status) &&
       isTRUE(posterior_density_status[["present"]]) &&
       isTRUE(posterior_density_status[["relevant"]]) &&
       !isTRUE(posterior_density_status[["valid"]])){
      stop(
        "Precomputed posterior density metadata is present but invalid.",
        call. = FALSE
      )
    }
    stop(
      "'density_method = \"precomputed\"' requires valid posterior ordinate ",
      "or posterior density metadata for the requested null hypothesis.",
      call. = FALSE
    )
  }

  if(mean(posterior == null_hypothesis) > 0.05){
    stop(
      "There is a considerable cluster of posterior samples at the exact null hypothesis value. The Savage-Dickey density ratio is invalid.",
      call. = FALSE
    )
  }
  if(!is.null(stored_posterior_density) && nrow(stored_posterior_density[["point_masses"]]) > 0L){
    point_masses <- stored_posterior_density[["point_masses"]]
    point_tol <- sqrt(.Machine$double.eps) * max(1, abs(null_hypothesis))
    null_point_mass <- sum(point_masses[["mass"]][abs(point_masses[["x"]] - null_hypothesis) <= point_tol])
    if(null_point_mass > 0){
      stop(
        "Stored posterior density contains a point mass at the exact null hypothesis value. The Savage-Dickey density ratio is invalid.",
        call. = FALSE
      )
    }
  }
  if(.prior_linear_density_point_mass(prior, null_hypothesis) > 0){
    stop(
      "There is a point mass in the prior at the exact null hypothesis value. The Savage-Dickey density ratio is invalid.",
      call. = FALSE
    )
  }

  prior_range <- range(c(
    if(!is.null(prior$density)) prior$density$x else NULL,
    if(!is.null(prior$points) && nrow(prior$points) > 0) prior$points$x else NULL
  ))
  if(null_hypothesis < prior_range[1] || null_hypothesis > prior_range[2]){
    warnings <- c(warnings, "Prior density does not span both sides of the null hypothesis. Check whether the prior distribution contains the null hypothesis in the first place. The Savage-Dickey density ratio is likely to be invalid.")
  }
  posterior_range <- range(posterior)
  posterior_support_bounds <- .posterior_support_bounds(
    posterior,
    interval_only = TRUE
  )
  if(is.null(posterior_support_bounds)){
    stored_support <- .posterior_support_from_attribute(stored_posterior_density_support)
    if(!is.null(stored_support) && isTRUE(stored_support$exact) &&
       .posterior_support_has_interval(stored_support)){
      posterior_support_bounds <- stored_support$bounds
    }
  }
  null_at_support_boundary <- FALSE
  if(!is.null(posterior_support_bounds)){
    boundary_tol <- sqrt(.Machine$double.eps) * max(1, abs(null_hypothesis), abs(posterior_support_bounds[is.finite(posterior_support_bounds)]))
    null_at_support_boundary <- any(is.finite(posterior_support_bounds) & abs(null_hypothesis - posterior_support_bounds) <= boundary_tol)
  }
  if(!is.null(stored_posterior_density) && is.null(stored_posterior_ordinate)){
    posterior_range <- range(stored_posterior_density[["x"]], finite = TRUE)
    if(null_hypothesis < posterior_range[1] || null_hypothesis > posterior_range[2]){
      stop(
        "Stored posterior density does not span both sides of the null hypothesis.",
        call. = FALSE
      )
    }
  }
  if(is.null(stored_posterior_ordinate) &&
     (null_hypothesis < posterior_range[1] || null_hypothesis > posterior_range[2]) &&
     !isTRUE(null_at_support_boundary)){
    warnings <- c(warnings, "Posterior samples do not span both sides of the null hypothesis. The Savage-Dickey density ratio is likely to be overestimated.")
  }

  kde_height <- function(support = NULL){
    height <- .Savage_Dickey_BF.kd(
      samples         = posterior,
      null_hypothesis = null_hypothesis,
      support         = support
    )
    support_warning <- attr(height, "posterior_support_warning", exact = TRUE)
    if(!is.null(support_warning)){
      warnings <<- c(warnings, support_warning)
    }
    if(isTRUE(attr(height, "boundary_reflection", exact = TRUE))){
      posterior_density_boundary_reflection <<- TRUE
    }
    if(isTRUE(attr(height, "posterior_support_exclusion", exact = TRUE))){
      posterior_density_source <<- "exact_support_exclusion"
    }
    support_bounds <- attr(height, "posterior_support_bounds", exact = TRUE)
    if(!is.null(support_bounds)){
      posterior_density_support_bounds <<- support_bounds
    }
    height
  }

  if(normal_approximation){
    posterior_height <- .Savage_Dickey_BF.normal(posterior, null_hypothesis)
  }else if(!is.null(stored_posterior_ordinate)){
    support_exclusion <- .Savage_Dickey_BF.support_exclusion(
      posterior,
      null_hypothesis = null_hypothesis,
      source_support  = stored_posterior_density_support
    )
    warnings <- c(warnings, support_exclusion[["warnings"]])
    if(isTRUE(support_exclusion[["excluded"]])){
      fallback_warning <- paste0(
        "Exact ", support_exclusion[["source_label"]],
        " excludes the null hypothesis. Ignoring the precomputed posterior ordinate."
      )
      warnings <- c(warnings, fallback_warning)
      posterior_density_fallback_warnings <- c(posterior_density_fallback_warnings, fallback_warning)
      posterior_height <- 0
      posterior_density_source <- "exact_support_exclusion"
      posterior_density_support_bounds <- support_exclusion[["bounds"]]
    }else{
      posterior_height <- stored_posterior_ordinate[["y"]]
      posterior_density_source <- "precomputed"
      BF_error_percent <- .posterior_ordinate_bf_error_percent(stored_posterior_ordinate)
    }
  }else if(!is.null(stored_posterior_density)){
    support_exclusion <- .Savage_Dickey_BF.support_exclusion(
      posterior,
      null_hypothesis = null_hypothesis,
      source_support  = stored_posterior_density[["support"]]
    )
    warnings <- c(warnings, support_exclusion[["warnings"]])
    if(isTRUE(support_exclusion[["excluded"]])){
      fallback_warning <- paste0(
        "Exact ", support_exclusion[["source_label"]],
        " excludes the null hypothesis. Ignoring the precomputed posterior density."
      )
      warnings <- c(warnings, fallback_warning)
      posterior_density_fallback_warnings <- c(posterior_density_fallback_warnings, fallback_warning)
      posterior_height <- 0
      posterior_density_source <- "exact_support_exclusion"
      posterior_density_support_bounds <- support_exclusion[["bounds"]]
    }else{
      posterior_height <- .posterior_density_height(stored_posterior_density, null_hypothesis)
    }
    if(!isTRUE(support_exclusion[["excluded"]]) &&
       (!is.finite(posterior_height) || posterior_height <= 0)){
      stop(
        "Stored posterior density has zero or non-finite height at the null hypothesis.",
        call. = FALSE
      )
    }else if(!isTRUE(support_exclusion[["excluded"]])){
      posterior_density_source <- "precomputed"
      BF_error_percent <- .posterior_density_bf_error_percent(stored_posterior_density, null_hypothesis)
    }
  }else{
    posterior_height <- kde_height(stored_posterior_density_support)
  }
  prior_height <- .prior_linear_density_height(prior, null_hypothesis)
  if(!is.finite(prior_height) || prior_height <= 0){
    warnings <- c(
      warnings,
      "Prior density at the null hypothesis value is zero or non-finite. The Savage-Dickey density ratio is invalid."
    )
  }

  if(!silent && !is.null(warnings)){
    sapply(warnings, warning, call. = FALSE)
  }

  BF <- exp(log(prior_height) - log(posterior_height))

  if(!is.null(warnings)){
    attr(BF, "warnings") <- warnings
  }
  if(is.finite(BF_error_percent)){
    attr(BF, "BF_error_percent") <- BF_error_percent
  }
  attr(BF, "posterior_density_source") <- posterior_density_source
  if(isTRUE(posterior_density_boundary_reflection)){
    attr(BF, "posterior_density_boundary_reflection") <- TRUE
  }
  if(!is.null(posterior_density_support_bounds)){
    attr(BF, "posterior_density_support") <- posterior_density_support_bounds
  }
  if(length(posterior_density_fallback_warnings) > 0L){
    attr(BF, "posterior_density_fallback") <- TRUE
    attr(BF, "posterior_density_fallback_warnings") <- posterior_density_fallback_warnings
  }

  return(BF)
}
.Savage_Dickey_BF.normal <- function(samples, null_hypothesis){

  height <- stats::dnorm(null_hypothesis, mean = mean(samples), sd = stats::sd(samples))

  return(height)
}

.Savage_Dickey_BF.support_exclusion <- function(samples, null_hypothesis,
                                                 source_support = NULL){

  supports <- list(source_support, .posterior_support_get(samples))
  support_labels <- c("stored posterior density support", "posterior support")
  valid_bounds <- NULL
  support_warnings <- NULL

  for(i in seq_along(supports)){
    support <- supports[[i]]
    support <- .posterior_support_from_attribute(support)
    if(is.null(support) || !isTRUE(support$exact)){
      next
    }

    support_info <- .posterior_support_for_kde(samples, support = support)
    if(!is.null(support_info[["warning"]])){
      support_warnings <- c(support_warnings, support_info[["warning"]])
    }
    if(is.null(support_info[["bounds"]])){
      next
    }
    if(is.null(valid_bounds)){
      valid_bounds <- support_info[["bounds"]]
    }

    if(.posterior_support_excludes_value(
      .posterior_support_new(support_info[["bounds"]]),
      null_hypothesis
    )){
      return(list(
        excluded    = TRUE,
        bounds      = support_info[["bounds"]],
        source_label = support_labels[[i]],
        warnings    = support_warnings
      ))
    }
  }

  list(excluded = FALSE, bounds = valid_bounds, source_label = NULL,
       warnings = support_warnings)
}

.Savage_Dickey_BF.kd     <- function(samples, null_hypothesis, support = NULL){

  sample_values <- as.numeric(samples)
  sample_values <- sample_values[is.finite(sample_values)]
  support_info <- .posterior_support_for_kde(samples, support = support)
  support_bounds <- support_info[["bounds"]]

  if(!is.null(support_bounds)){
    if(null_hypothesis < support_bounds[1] || null_hypothesis > support_bounds[2]){
      height <- 0
      attr(height, "posterior_support_bounds") <- support_bounds
      attr(height, "posterior_support_exclusion") <- TRUE
      return(height)
    }

    if(any(is.finite(support_bounds))){
      density_args <- list(
        x      = sample_values,
        n      = 512L,
        bounds = support_bounds,
        na.rm  = TRUE
      )
      if(is.finite(support_bounds[1])){
        density_args[["from"]] <- support_bounds[1]
      }
      if(is.finite(support_bounds[2])){
        density_args[["to"]] <- support_bounds[2]
      }
      density_posterior <- do.call(.density_kde_boundary, density_args)
      height <- stats::approx(
        density_posterior$x,
        density_posterior$y,
        xout  = null_hypothesis,
        yleft = 0,
        yright = 0
      )[["y"]]
      attr(height, "boundary_reflection") <- isTRUE(attr(density_posterior, "boundary_reflection"))
      attr(height, "posterior_support_bounds") <- support_bounds
      return(height)
    }
  }

  if(null_hypothesis < min(sample_values) || null_hypothesis > max(sample_values)){
    height <- 0
  }else{
    density_posterior <- stats::density(sample_values)
    height <- stats::approx(
      density_posterior$x,
      density_posterior$y,
      xout   = null_hypothesis,
      yleft  = 0,
      yright = 0
    )[["y"]]
  }

  if(!is.null(support_info[["warning"]])){
    attr(height, "posterior_support_warning") <- support_info[["warning"]]
  }

  return(height)
}


