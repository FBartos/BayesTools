#' Shrink regression coefficients towards zero
#'
#' @description Ridge uses a zero-mean Normal prior and lasso uses a
#' zero-location Laplace prior. The horseshoe gives each fitted coefficient
#' its own standard half-Cauchy local scale. All members of a named group
#' share the global scale and, when specified, the slab variance.
#'
#' @param family One of `"ridge"`, `"lasso"`, and `"horseshoe"`.
#' @param scale A positive number fixing the global scale, or a proper scalar
#' continuous prior with positive support from which to learn it. The scale
#' is in fitted coefficient units. For lasso it is the Laplace scale (the
#' conditional standard deviation is \eqn{\sqrt{2}} times this scale).
#' @param terms Optional character vector of fixed-effect term labels.
#' Place a group under its own name in the formula prior list. Omitting
#' `terms` selects the non-intercept terms without a prior of their own.
#' Separate named groups have independent scales, even if they use the
#' same scale prior. A multi-column term uses its fitted coefficient basis.
#' @param slab_scale Positive scale of the Student-t slab for the regularized
#' horseshoe. The default `NULL` selects the classical horseshoe.
#' @param slab_df Positive degrees of freedom of the Student-t slab.
#' @param x An R object.
#'
#' @details Conditional on global scale \eqn{\tau} and local scale
#' \eqn{\lambda_j}, a horseshoe coefficient is Normal with mean zero and
#' standard deviation \eqn{s_j = \tau\lambda_j}. With a finite slab its
#' conditional variance is \eqn{c^2 s_j^2/(c^2+s_j^2)}, where the shared
#' \eqn{c^2} follows an inverse-gamma distribution with shape
#' \eqn{\nu/2} and scale \eqn{\nu s^2/2} for \eqn{\nu=} `slab_df`
#' and \eqn{s=} `slab_scale`.
#'
#' @return With numeric scale and omitted terms, ridge and lasso return
#' ordinary scalar priors accepted by [prior_spike_and_slab()] and
#' [prior_mixture()]. Otherwise the result declares a coefficient group for
#' [JAGS_formula()]. `is.prior.regularized()` identifies group declarations.
#'
#' @references \insertRef{piironen2017ejs_sparsity}{BayesTools}
#' @examples
#' prior_regularized("ridge", scale = 0.5)
#' prior_regularized("lasso", scale = 0.5)
#' global <- prior("cauchy", list(0, 0.2), truncation = list(0, Inf))
#' priors <- list(
#'   intercept = prior("normal", list(0, 2)),
#'   shrinkage = prior_regularized("horseshoe", scale = global,
#'                                 terms = c("wt", "hp"), slab_scale = 2)
#' )
#' JAGS_formula(~ wt + hp, "mu", mtcars, prior_list = priors)
#' @export
prior_regularized <- function(family, scale, terms = NULL, slab_scale = NULL, slab_df = 4){

  check_char(family, "family", allow_values = c("ridge", "lasso", "horseshoe"), allow_NA = FALSE)
  check_char(terms, "terms", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  if(!is.null(terms) && (!length(terms) || any(!nzchar(terms)) || anyDuplicated(terms) ||
                         any(terms %in% c("intercept", "(Intercept)")))){
    stop("Supply distinct non-intercept term labels in 'terms'.", call. = FALSE)
  }
  if(is.prior(scale)){
    if(!is.prior.simple(scale) || is.prior.factor(scale) || is.prior.discrete(scale) ||
       is.prior.PET(scale) || is.prior.PEESE(scale) || .is_prior_expression(scale) ||
       inherits(scale, "prior.regularized_term") || !is.null(attr(scale, "multiply_by", exact = TRUE))){
      stop("Supply a positive number or a continuous scalar prior on positive values in 'scale'.", call. = FALSE)
    }
    if(is.prior.point(scale)){
      scale <- scale$parameters$location
    }else if(!isTRUE(scale$truncation$lower >= 0 && scale$truncation$upper > scale$truncation$lower) ||
             any(!is.finite(unlist(scale$parameters)))){
      stop("The prior in 'scale' must have finite parameters and positive support.", call. = FALSE)
    }
  }
  if(!is.prior(scale)){
    .regularized_check_positive(scale, "scale")
  }
  if(!is.null(slab_scale)){
    if(family != "horseshoe"){
      stop("'slab_scale' applies only to the horseshoe.", call. = FALSE)
    }
    .regularized_check_positive(slab_scale, "slab_scale")
    .regularized_check_positive(slab_df, "slab_df")
  }else if(!missing(slab_df)){
    stop("Specify 'slab_scale' when supplying 'slab_df'.", call. = FALSE)
  }
  spec <- structure(list(family = family, terms = terms, scale = scale,
                         slab_scale = slab_scale, slab_df = if(!is.null(slab_scale)) slab_df),
                    class = c("prior", "prior.regularized"))
  if(is.numeric(scale) && is.null(terms) && family != "horseshoe"){
    scalar <- prior(if(family == "ridge") "normal" else "laplace", list(0, scale))
    attr(scalar, "regularized_spec") <- spec
    return(scalar)
  }
  spec
}

.regularized_check_positive <- function(x, name){

  check_real(x, name, lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if(!is.finite(x)){
    stop(paste0("The '", name, "' argument must be finite."), call. = FALSE)
  }
}

#' @rdname prior_regularized
#' @export
is.prior.regularized <- function(x){

  inherits(x, "prior.regularized")
}

.print_prior_regularized <- function(x, digits, short_name){

  bound <- inherits(x, "prior.regularized_term")
  spec <- if(bound) x$regularization else x
  label <- switch(spec$family, ridge = "Ridge", lasso = "Lasso",
    horseshoe = if(is.null(spec$slab_scale)) "Horseshoe" else "Regularized horseshoe")
  if(short_name) return(label)
  scale <- if(is.prior(spec$scale)) print(spec$scale, digits_estimates = digits, silent = TRUE) else
    format(signif(spec$scale, max(1, digits)))
  entries <- paste0("scale = ", scale)
  if(!is.null(spec$slab_scale)){
    entries <- c(entries, paste0("slab scale = ", spec$slab_scale), paste0("slab df = ", spec$slab_df))
  }
  members <- if(bound) paste0("group = ", spec$group) else
    paste0("terms = ", if(is.null(spec$terms)) "all unassigned slopes" else paste(spec$terms, collapse = ", "))
  paste0(label, "(", paste(c(entries, members), collapse = ", "), ")")
}
