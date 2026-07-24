#' @title Assess convergence of a runjags model
#'
#' @description Checks whether the supplied \link[runjags]{runjags-package} model
#' satisfied convergence criteria.
#' @param fit a runjags model
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#' @param max_Rhat maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}.
#' @param min_ESS minimum effective sample size. Defaults to \code{500}.
#' @param max_error maximum MCMC error. Defaults to \code{0.01}.
#' @param max_SD_error maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.
#' @param add_parameters vector of additional parameter names that should be used
#' (only allows removing last, fixed, omega element if omega is tracked manually).
#' @param fail_fast whether the function should stop after the first failed convergence check.
#' @param check_indicators whether model indicator variables should be included
#' in convergence checks. Defaults to \code{FALSE}.
#'
#' @examples \dontrun{
#' # simulate data
#' set.seed(1)
#' data <- list(
#'   x = rnorm(10),
#'   N = 10
#' )
#' data$x
#'
#' # define priors
#' priors_list <- list(mu = prior("normal", list(0, 1)))
#'
#' # define likelihood for the data
#' model_syntax <-
#'   "model{
#'     for(i in 1:N){
#'       x[i] ~ dnorm(mu, 1)
#'     }
#'   }"
#'
#' # fit the models
#' fit <- JAGS_fit(model_syntax, data, priors_list)
#' JAGS_check_convergence(fit, priors_list)
#' }
#' @return \code{JAGS_check_convergence} returns a boolean
#' indicating whether the model converged or not, with an
#' attribute 'errors' carrying the failed convergence checks (if any).
#'
#' @seealso [JAGS_fit()]
#' @export
JAGS_check_convergence <- function(fit, prior_list, max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, add_parameters = NULL, fail_fast = FALSE, check_indicators = FALSE){

  # check input
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  check_list(prior_list, "prior_list", allow_NULL = TRUE)
  if(!is.null(prior_list) && any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_real(max_Rhat,     "max_Rhat",     lower = 1, allow_NULL = TRUE)
  check_real(min_ESS,      "min_ESS",      lower = 0, allow_NULL = TRUE)
  check_real(max_error,    "max_error",    lower = 0, allow_NULL = TRUE)
  check_real(max_SD_error, "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE)
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE)
  check_bool(check_indicators, "check_indicators")

  # extract samples and parameter information
  mcmc_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  mcmc_samples      <- do.call(rbind, mcmc_samples_list)

  # build remove_parameters list: point priors, spike priors, indicators, inclusions
  remove_params <- c(
    # point priors
    names(prior_list)[sapply(prior_list, is.prior.point)],
    # mixture with single point prior
    names(prior_list)[sapply(prior_list, function(p) {
      is.prior.mixture(p) && length(p) == 1 && is.prior.point(p[[1]])
    })],
    # add_parameters that should be excluded
    add_parameters
  )

  # use helper to remove auxiliary parameters
  cleaned <- .remove_auxiliary_parameters(mcmc_samples, prior_list, remove_params)
  mcmc_samples <- cleaned$model_samples

  # remove auxiliary inclusion probabilities and, by default, model indicators
  indicator_cols <- grepl("_indicator(\\[[^]]+\\])?$", colnames(mcmc_samples))
  inclusion_cols <- grepl("_inclusion(\\[[^]]+\\])?$", colnames(mcmc_samples))
  mcmc_samples <- mcmc_samples[, !(inclusion_cols | (!check_indicators & indicator_cols)), drop = FALSE]

  if(ncol(mcmc_samples) == 0){
    return(TRUE)
  }

  # convert back to mcmc.list for convergence checks
  n_chains <- length(mcmc_samples_list)
  samples_per_chain <- nrow(mcmc_samples) / n_chains
  mcmc_samples_list_cleaned <- lapply(1:n_chains, function(i) {
    start_idx <- (i - 1) * samples_per_chain + 1
    end_idx <- i * samples_per_chain
    coda::as.mcmc(mcmc_samples[start_idx:end_idx, , drop = FALSE])
  })
  mcmc_samples <- coda::as.mcmc.list(mcmc_samples_list_cleaned)

  ### check the convergence
  fails <- NULL

  # assess R-hat
  if(!is.null(max_Rhat)){
    if(length(fit$mcmc) == 1){
      warning("Only one chain was run. R-hat cannot be computed.", immediate. = TRUE)
    }else{
      temp_Rhat <- coda::gelman.diag(mcmc_samples, multivariate = FALSE, autoburnin = FALSE)$psrf
      temp_Rhat[is.na(temp_Rhat)] <- 1
      temp_Rhat <- max(temp_Rhat)
      if(temp_Rhat > max_Rhat){
        fails <- c(fails, paste0("R-hat ", round(temp_Rhat, 3), " is larger than the set target (", max_Rhat, ")."))
        if(fail_fast){
          return(FALSE)
        }
      }
    }
  }

  if(!is.null(min_ESS)){
    temp_ESS <- coda::effectiveSize(mcmc_samples)
    temp_ESS[is.nan(temp_ESS) | temp_ESS == 0] <- Inf
    temp_ESS <- min(temp_ESS)
    if(temp_ESS < min_ESS){
      fails <- c(fails, paste0("ESS ", round(temp_ESS), " is lower than the set target (", min_ESS, ")."))
      if(fail_fast){
        return(FALSE)
      }
    }
  }

  # compute the MCMC error and & SD error
  if(!(is.null(max_error) && is.null(max_SD_error))){
    temp_summary <- summary(mcmc_samples, quantiles = NULL)$statistics
    if(is.null(dim(temp_summary))){
      temp_summary <- t(temp_summary)
    }
  }


  if(!is.null(max_error)){
    temp_error    <- temp_summary[,"Time-series SE"]
    temp_error[is.na(temp_error)] <- 0
    temp_error    <- max(temp_error)
    if(temp_error > max_error){
      fails <- c(fails, paste0("MCMC error ", round(temp_error, 5), " is larger than the set target (", max_error, ")."))
      if(fail_fast){
        return(FALSE)
      }
    }
  }

  if(!is.null(max_SD_error)){
    temp_error_SD <- temp_summary[,"Time-series SE"] / temp_summary[,"SD"]
    temp_error_SD[is.na(temp_error_SD)] <- 0
    temp_error_SD <- max(temp_error_SD)
    if(temp_error_SD > max_SD_error){
      fails <- c(fails, paste0("MCMC SD error ", round(temp_error_SD, 3), " is larger than the set target (", max_SD_error, ")."))
      if(fail_fast){
        return(FALSE)
      }
    }
  }

  converged <- length(fails) == 0
  attr(converged, "errors") <- fails
  return(converged)
}
