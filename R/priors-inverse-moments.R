.prior_inverse_moment_result <- function(log_value, method = "closed_form",
                                         error_estimate = 0, reason = NULL,
                                         integration = NULL, value = NULL){

  threshold <- log1p(.prior_linear_density_refinement_tolerance()$relative)
  available <- is.finite(log_value) && is.finite(error_estimate) &&
    error_estimate <= threshold && is.null(reason)
  if(is.null(reason) && !available){
    reason <- if(!is.finite(log_value)) "The log inverse moment is outside supported arithmetic range." else
      "The log inverse moment is insufficiently resolved for the ordinate accuracy criterion."
  }
  if(is.null(value)) value <- if(is.na(log_value)) NA_real_ else exp(log_value)
  list(value = value, log_value = log_value, method = method,
       available = available, reason = reason, error_estimate = error_estimate,
       error_estimate_kind = "operand floating-resolution diagnostic, not a special-function error bound",
       integration = integration)
}

.prior_inverse_moment_offset <- function(moment, addends){

  if(!isTRUE(moment$available)) return(moment)
  value <- sum(c(addends, moment$log_value))
  estimate <- moment$error_estimate +
    8 * .Machine$double.eps * sum(abs(c(addends, moment$log_value, value)))
  .prior_inverse_moment_result(value, moment$method, estimate,
                              integration = moment$integration)
}

.prior_inverse_moment_provenance <- function(moment){

  moment[c("value", "log_value", "method", "available", "reason", "error_estimate",
           "error_estimate_kind")]
}
