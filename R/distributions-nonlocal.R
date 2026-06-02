.nonlocal_validate_finite <- function(x, name){

  if(!is.finite(x)){
    stop(paste0("The '", name, "' must be finite."), call. = FALSE)
  }

  invisible(NULL)
}

.nonlocal_validate_mode_tau <- function(parameters, distribution){

  has_mode <- "mode" %in% names(parameters)
  has_tau  <- "tau"  %in% names(parameters)

  if(has_mode == has_tau){
    stop(paste0("The ", distribution, " prior distribution requires exactly one of 'mode' or 'tau'."), call. = FALSE)
  }

  invisible(NULL)
}

.nonlocal_validate_order <- function(order){

  check_int(order, "order", lower = 1, allow_NA = FALSE)
  .nonlocal_validate_finite(order, "order")

  as.integer(order)
}

.nonlocal_parameters_moment <- function(parameters){

  if(is.null(names(parameters))){
    if(length(parameters) != 1L){
      stop("Moment prior distribution positional input requires a single 'mode' parameter. Use named parameters for 'tau', 'order', or 'location'.", call. = FALSE)
    }
    names(parameters) <- "mode"
  }else{
    if(anyDuplicated(names(parameters)[nzchar(names(parameters))])){
      stop("Moment prior distribution parameters must not contain duplicate names.", call. = FALSE)
    }
    unsupported <- !names(parameters) %in% c("mode", "tau", "order", "location")
    if(any(unsupported)){
      stop(paste0("Parameters ", paste(paste0("'", names(parameters)[unsupported], "'"), collapse = ", "), " are not supported for a moment distribution."), call. = FALSE)
    }
    .nonlocal_validate_mode_tau(parameters, "moment")
    if(!"location" %in% names(parameters)){
      parameters[["location"]] <- 0
    }
  }

  if(!"order" %in% names(parameters)){
    parameters[["order"]] <- 1
  }
  if(!"location" %in% names(parameters)){
    parameters[["location"]] <- 0
  }
  .nonlocal_validate_mode_tau(parameters, "moment")

  order    <- .nonlocal_validate_order(parameters[["order"]])
  location <- parameters[["location"]]
  check_real(location, "location", allow_NA = FALSE)
  .nonlocal_validate_finite(location, "location")

  if("mode" %in% names(parameters)){
    mode <- parameters[["mode"]]
    check_real(mode, "mode", allow_NA = FALSE)
    .nonlocal_validate_finite(mode, "mode")
    mode <- abs(mode)
    if(mode == 0){
      stop("The 'mode' must be finite and nonzero after taking its absolute value.", call. = FALSE)
    }
    tau <- mode^2 / (2 * order)
  }else{
    tau <- parameters[["tau"]]
    check_real(tau, "tau", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
    .nonlocal_validate_finite(tau, "tau")
    mode <- sqrt(2 * order * tau)
  }

  list(
    mode     = mode,
    tau      = tau,
    order    = order,
    location = location
  )
}

.nonlocal_parameters_invmoment <- function(parameters){

  if(is.null(names(parameters))){
    if(length(parameters) == 2L){
      names(parameters) <- c("mode", "df")
    }else{
      stop("Inverse-moment prior distribution positional input requires 'mode' and 'df'. Use named parameters for 'tau', 'order', 'df'/'nu', or 'location'.", call. = FALSE)
    }
  }else{
    if(anyDuplicated(names(parameters)[nzchar(names(parameters))])){
      stop("Inverse-moment prior distribution parameters must not contain duplicate names.", call. = FALSE)
    }
    unsupported <- !names(parameters) %in% c("mode", "tau", "order", "df", "nu", "location")
    if(any(unsupported)){
      stop(paste0("Parameters ", paste(paste0("'", names(parameters)[unsupported], "'"), collapse = ", "), " are not supported for an inverse-moment distribution."), call. = FALSE)
    }
    .nonlocal_validate_mode_tau(parameters, "inverse-moment")
    if("df" %in% names(parameters) && "nu" %in% names(parameters)){
      stop("Use only one of 'df' or 'nu' for an inverse-moment prior distribution.", call. = FALSE)
    }
    if(!"df" %in% names(parameters) && "nu" %in% names(parameters)){
      parameters[["df"]] <- parameters[["nu"]]
      parameters[["nu"]] <- NULL
    }
    if(!"df" %in% names(parameters)){
      stop("Inverse-moment prior distribution requires a 'df' parameter.", call. = FALSE)
    }
    if(!"location" %in% names(parameters)){
      parameters[["location"]] <- 0
    }
  }

  if(!"order" %in% names(parameters)){
    parameters[["order"]] <- 1
  }
  if(!"location" %in% names(parameters)){
    parameters[["location"]] <- 0
  }
  .nonlocal_validate_mode_tau(parameters, "inverse-moment")

  order    <- .nonlocal_validate_order(parameters[["order"]])
  df       <- parameters[["df"]]
  location <- parameters[["location"]]

  check_real(df, "df", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  .nonlocal_validate_finite(df, "df")
  check_real(location, "location", allow_NA = FALSE)
  .nonlocal_validate_finite(location, "location")

  if("mode" %in% names(parameters)){
    mode <- parameters[["mode"]]
    check_real(mode, "mode", allow_NA = FALSE)
    .nonlocal_validate_finite(mode, "mode")
    mode <- abs(mode)
    if(mode == 0){
      stop("The 'mode' must be finite and nonzero after taking its absolute value.", call. = FALSE)
    }
    tau <- mode^2 * ((df + 1) / (2 * order))^(1 / order)
  }else{
    tau <- parameters[["tau"]]
    check_real(tau, "tau", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
    .nonlocal_validate_finite(tau, "tau")
    mode <- sqrt(tau) * (2 * order / (df + 1))^(1 / (2 * order))
  }

  list(
    mode     = mode,
    tau      = tau,
    order    = order,
    df       = df,
    location = location
  )
}

.dmoment_prior <- function(x, location, tau, order, log = FALSE){

  .check_log(log)
  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_moment_d", x, location, tau, as.numeric(order), log, PACKAGE = "BayesTools")
}

.pmoment_prior <- function(q, location, tau, order, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_moment_p", q, location, tau, as.numeric(order), lower.tail, log.p, PACKAGE = "BayesTools")
}

.qmoment_prior <- function(p, location, tau, order, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_moment_q", p, location, tau, as.numeric(order), lower.tail, log.p, PACKAGE = "BayesTools")
}

.rmoment_prior <- function(n, location, tau, order){

  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_moment_r", as.integer(n), location, tau, as.numeric(order), PACKAGE = "BayesTools")
}

.dinvmoment_prior <- function(x, location, tau, order, df, log = FALSE){

  .check_log(log)
  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_invmoment_d", x, location, tau, as.numeric(order), df, log, PACKAGE = "BayesTools")
}

.pinvmoment_prior <- function(q, location, tau, order, df, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_invmoment_p", q, location, tau, as.numeric(order), df, lower.tail, log.p, PACKAGE = "BayesTools")
}

.qinvmoment_prior <- function(p, location, tau, order, df, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_invmoment_q", p, location, tau, as.numeric(order), df, lower.tail, log.p, PACKAGE = "BayesTools")
}

.rinvmoment_prior <- function(n, location, tau, order, df){

  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_invmoment_r", as.integer(n), location, tau, as.numeric(order), df, PACKAGE = "BayesTools")
}
