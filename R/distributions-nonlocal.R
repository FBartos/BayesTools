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

  check_int(order, "order", lower = 1, upper = .Machine$integer.max,
            allow_NA = FALSE)
  .nonlocal_validate_finite(order, "order")

  order <- as.integer(round(order))
  if(is.na(order)){
    stop("The 'order' argument cannot be represented as an integer.",
         call. = FALSE)
  }

  order
}

.nonlocal_validate_tau <- function(tau, name){

  if(identical(name, "mode") && (!is.finite(tau) || tau <= 0)){
    stop("The supplied 'mode' yields a derived 'tau' outside representable positive range.", call. = FALSE)
  }
  check_real(tau, name, lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  .nonlocal_validate_finite(tau, name)

  invisible(NULL)
}

.nonlocal_validate_mode <- function(mode, name){

  check_real(mode, name, allow_NA = FALSE)
  .nonlocal_validate_finite(mode, name)
  mode <- abs(mode)
  if(mode == 0){
    stop(paste0(
      "The '", name, "' must be finite and nonzero after taking its absolute value."
    ), call. = FALSE)
  }

  mode
}

.nonlocal_validate_derived_mode <- function(mode, name){

  .nonlocal_validate_finite(mode, name)
  if(mode == 0){
    stop(paste0("The '", name, "' implies a zero mode."), call. = FALSE)
  }

  invisible(NULL)
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
    mode <- .nonlocal_validate_mode(parameters[["mode"]], "mode")
    tau <- mode^2 / (2 * order)
    if(!is.finite(tau) || mode^2 < .Machine$double.xmin){
      tau <- exp(2 * log(mode) - log(2) - log(order))
    }
    .nonlocal_validate_tau(tau, "mode")
  }else{
    tau <- parameters[["tau"]]
    .nonlocal_validate_tau(tau, "tau")
    scaled_tau <- 2 * order * tau
    mode <- if(is.finite(scaled_tau) && scaled_tau >= .Machine$double.xmin){
      sqrt(scaled_tau)
    }else{
      sqrt(tau) * sqrt(2 * order)
    }
    .nonlocal_validate_derived_mode(mode, "tau")
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
    mode <- .nonlocal_validate_mode(parameters[["mode"]], "mode")
    tau <- mode^2 * ((df + 1) / (2 * order))^(1 / order)
    if(!is.finite(tau) || mode^2 < .Machine$double.xmin){
      tau <- exp(2 * log(mode) + (log(df + 1) - log(2) - log(order)) / order)
    }
    .nonlocal_validate_tau(tau, "mode")
  }else{
    tau <- parameters[["tau"]]
    .nonlocal_validate_tau(tau, "tau")
    mode <- sqrt(tau) * (2 * order / (df + 1))^(1 / (2 * order))
    .nonlocal_validate_derived_mode(mode, "tau")
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
  out <- .Call("BayesTools_moment_d", x, location, tau, as.numeric(order), log, PACKAGE = "BayesTools")
  .prior_numerical_result(out, x, .prior_numerical_nonlocal_domain(location, tau, order),
    "density", "moment", if(log) "log" else "natural",
    interior = is.finite(x) & x != location, rng = FALSE)
}

.pmoment_prior <- function(q, location, tau, order, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_moment_p", q, location, tau, as.numeric(order), lower.tail, log.p, PACKAGE = "BayesTools")
  .prior_numerical_result(out, q, .prior_numerical_nonlocal_domain(location, tau, order),
    "distribution", "moment", if(log.p) "log" else "natural",
    interior = is.finite(q) & q != location, rng = FALSE)
}

.qmoment_prior <- function(p, location, tau, order, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_moment_q", p, location, tau, as.numeric(order), lower.tail, log.p, PACKAGE = "BayesTools")
  .prior_numerical_result(out, p, .prior_numerical_nonlocal_domain(location, tau, order),
    "quantile", "moment", "natural",
    interior = if(log.p) is.finite(p) & p < 0 else p > 0 & p < 1, rng = FALSE)
}

.rmoment_prior <- function(n, location, tau, order){

  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_moment_r", as.integer(n), location, tau, as.numeric(order), PACKAGE = "BayesTools")
  .prior_numerical_result(out, rep(0, length(out)), .prior_numerical_nonlocal_domain(location, tau, order),
    "sampling", "moment", "natural",
    interior = rep(TRUE, length(out)), rng = TRUE)
}

.dinvmoment_prior <- function(x, location, tau, order, df, log = FALSE){

  .check_log(log)
  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_invmoment_d", x, location, tau, as.numeric(order), df, log, PACKAGE = "BayesTools")
  .prior_numerical_result(out, x, .prior_numerical_nonlocal_domain(location, tau, order, df),
    "density", "invmoment", if(log) "log" else "natural",
    interior = is.finite(x) & x != location, rng = FALSE)
}

.pinvmoment_prior <- function(q, location, tau, order, df, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_invmoment_p", q, location, tau, as.numeric(order), df, lower.tail, log.p, PACKAGE = "BayesTools")
  .prior_numerical_result(out, q, .prior_numerical_nonlocal_domain(location, tau, order, df),
    "distribution", "invmoment", if(log.p) "log" else "natural",
    interior = is.finite(q) & q != location, rng = FALSE)
}

.qinvmoment_prior <- function(p, location, tau, order, df, lower.tail = TRUE, log.p = FALSE){

  .check_lower.tail(lower.tail)
  .check_log.p(log.p)
  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_invmoment_q", p, location, tau, as.numeric(order), df, lower.tail, log.p, PACKAGE = "BayesTools")
  .prior_numerical_result(out, p, .prior_numerical_nonlocal_domain(location, tau, order, df),
    "quantile", "invmoment", "natural",
    interior = if(log.p) is.finite(p) & p < 0 else p > 0 & p < 1, rng = FALSE)
}

.rinvmoment_prior <- function(n, location, tau, order, df){

  .BayesTools_require_native_nonlocal()
  out <- .Call("BayesTools_invmoment_r", as.integer(n), location, tau, as.numeric(order), df, PACKAGE = "BayesTools")
  .prior_numerical_result(out, rep(0, length(out)), .prior_numerical_nonlocal_domain(location, tau, order, df),
    "sampling", "invmoment", "natural",
    interior = rep(TRUE, length(out)), rng = TRUE)
}
