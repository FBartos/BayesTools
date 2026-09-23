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
  .Call("BayesTools_invgamma_d", x, shape, scale, log, PACKAGE = "BayesTools")
}

.pinvgamma_prior <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_invgamma()
  .Call("BayesTools_invgamma_p", q, shape, scale, lower.tail, log.p, PACKAGE = "BayesTools")
}

.qinvgamma_prior <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_invgamma()
  .Call("BayesTools_invgamma_q", p, shape, scale, lower.tail, log.p, PACKAGE = "BayesTools")
}

.rinvgamma_prior <- function(n, shape, scale){

  .BayesTools_require_native_invgamma()
  .Call("BayesTools_invgamma_r", as.integer(n), shape, scale, PACKAGE = "BayesTools")
}
