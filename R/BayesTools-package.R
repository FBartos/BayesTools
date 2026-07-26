##' BayesTools
##'
##' BayesTools: Provides tools for conducting Bayesian analyses. The package contains
##' functions for creating a wide range of prior distribution objects, mixing posterior
##' samples from JAGS and Stan models, plotting posterior distributions, and etc...
##' The tools for working with prior distribution span from visualization, generating JAGS
##' and bridgesampling syntax to basic functions such as rng, quantile, and distribution functions.
##'
##' Random-effect operations use the
##' \code{BayesTools.random_effects_memory_limit_bytes} option as a
##' deterministic pre-allocation guard. The default is \code{16 * 1024^3}
##' bytes (16 GiB). The estimate includes the primary numeric payload,
##' simultaneous working copies, row maps, group state, and covariance state
##' known before allocation. Set a different positive byte ceiling with
##' \code{options(BayesTools.random_effects_memory_limit_bytes = value)}.
##' Experts may explicitly use \code{Inf} to disable the guard after verifying
##' the operation's memory budget. Passing the guard does not promise that an
##' allocation will succeed because R, JAGS, and the operating system have
##' additional memory requirements.
##'
##'
##' @name BayesTools
##' @aliases BayesTools-package BayesTools
##' @docType package
##' @author František Bartoš \email{f.bartos96@@gmail.com}
##' @keywords package
##' @importFrom rlang .data
##' @importFrom Rdpack reprompt
##' @importFrom stats fft median setNames
##' @importFrom utils head tail
"_PACKAGE"
