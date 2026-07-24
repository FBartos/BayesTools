### create generic method
#' @title Creates generics for common statistical functions
#'
#' @description Density (pdf / lpdf), distribution
#' function (cdf / ccdf), quantile function (quant),
#' random generation (rng), mean, standard deviation (sd),
#' and marginal variants of the functions (mpdf, mlpf, mcdf,
#' mccdf, mquant).
#'
#' @param x main argument
#' @param ... unused arguments
#'
#' @return \code{pdf} (\code{mpdf}) and \code{lpdf} (\code{mlpdf}) give
#' the (marginal) density and the log of (marginal) density,
#' \code{cdf} (\code{mcdf}) and \code{ccdf} (\code{mccdf}) give the
#' (marginal) distribution and the complement of (marginal) distribution function,
#' \code{quant} (\code{mquant}) give the (marginal) quantile function,
#' and \code{rng} generates random deviates for an object of class 'prior'.
#'
#' The \code{pdf} function proceeds to PDF graphics device if \code{x} is a character.
#'
#' @export rng
#' @export cdf
#' @export ccdf
#' @export quant
#' @export lpdf
#' @export pdf
#' @export mcdf
#' @export mccdf
#' @export mquant
#' @export mlpdf
#' @export mpdf
#' @importFrom grDevices pdf
#' @name prior_functions_methods
NULL

#' @rdname prior_functions_methods
rng    <- function(x, ...){
  UseMethod("rng")
}
#' @rdname prior_functions_methods
cdf    <- function(x, ...){
  UseMethod("cdf")
}
#' @rdname prior_functions_methods
ccdf   <- function(x, ...){
  UseMethod("ccdf")
}
#' @rdname prior_functions_methods
quant  <- function(x, ...){
  UseMethod("quant")
}
#' @rdname prior_functions_methods
lpdf   <- function(x, ...){
  UseMethod("lpdf")
}
#' @rdname prior_functions_methods
pdf    <- function(x, ...){
  UseMethod("pdf")
}
#' @rdname prior_functions_methods
mcdf   <- function(x, ...){
  UseMethod("mcdf")
}
#' @rdname prior_functions_methods
mccdf  <- function(x, ...){
  UseMethod("mccdf")
}
#' @rdname prior_functions_methods
mquant <- function(x, ...){
  UseMethod("mquant")
}
#' @rdname prior_functions_methods
mlpdf  <- function(x, ...){
  UseMethod("mlpdf")
}
#' @rdname prior_functions_methods
mpdf   <- function(x, ...){
  UseMethod("mpdf")
}

#' @export
pdf.default  <- function(x, ...){
  grDevices::pdf(x, ...)
}
