.BayesTools_require_native_invgamma <- function(){

  if(isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))){
    return(invisible(TRUE))
  }

  .BayesTools_load_native_routines(pkgname = "BayesTools", warn = TRUE)
  if(!isTRUE(.BayesTools_native_routines_loaded(pkgname = "BayesTools"))){
    stop("BayesTools native inverse-gamma routines are not loaded.", call. = FALSE)
  }

  invisible(TRUE)
}

.dinvgamma_prior <- function(x, shape, scale, log = FALSE){

  .check_log(log)
  .BayesTools_require_native_invgamma()
  out <- .Call("BayesTools_invgamma_d", x, shape, scale, log, PACKAGE = "BayesTools")
  .prior_numerical_result(out, x, is.finite(shape) && shape > 0 && is.finite(scale) && scale > 0,
    "density", "invgamma", if(log) "log" else "natural",
    interior = is.finite(x) & x > 0, rng = FALSE)
}

.pinvgamma_prior <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_invgamma()
  out <- .Call("BayesTools_invgamma_p", q, shape, scale, lower.tail, log.p, PACKAGE = "BayesTools")
  .prior_numerical_result(out, q, is.finite(shape) && shape > 0 && is.finite(scale) && scale > 0,
    "distribution", "invgamma", if(log.p) "log" else "natural",
    interior = is.finite(q) & q > 0, rng = FALSE)
}

.qinvgamma_prior <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_invgamma()
  out <- .Call("BayesTools_invgamma_q", p, shape, scale, lower.tail, log.p, PACKAGE = "BayesTools")
  .prior_numerical_result(out, p, is.finite(shape) && shape > 0 && is.finite(scale) && scale > 0,
    "quantile", "invgamma", "natural",
    interior = if(log.p) is.finite(p) & p < 0 else p > 0 & p < 1, rng = FALSE)
}

.rinvgamma_prior <- function(n, shape, scale){

  .BayesTools_require_native_invgamma()
  out <- .Call("BayesTools_invgamma_r", as.integer(n), shape, scale, PACKAGE = "BayesTools")
  .prior_numerical_result(out, rep(0, length(out)), is.finite(shape) && shape > 0 && is.finite(scale) && scale > 0,
    "sampling", "invgamma", "natural",
    interior = rep(TRUE, length(out)), rng = TRUE)
}
