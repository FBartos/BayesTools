#' @title Create list of monitored parameters for 'JAGS' model
#'
#' @description Creates a vector of parameter names to be
#' monitored in a 'JAGS' model.
#'
#' @inheritParams JAGS_add_priors
#'
#' @return \code{JAGS_to_monitor} returns a character vector of
#' parameter names.
#'
#' @export
JAGS_to_monitor             <- function(prior_list){

  # return empty string in case that no prior was specified
  if(length(prior_list) == 0){
    return("")
  }

  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)
  .bt_validate_ordered_shared_allocations(prior_list)


  # add the monitored parameters
  monitor <- character()
  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.weightfunction(prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.phacking(prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.bias(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.mixture(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.factor(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      monitor <- c(monitor, .JAGS_monitor.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  monitor <- unique(monitor)

  if(length(monitor) == 0L){
    return("")
  }

  return(monitor)
}


.JAGS_monitor.simple         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!(is.prior.simple(prior) | is.prior.vector(prior) | is.prior.factor(prior)))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(prior[["distribution"]] %in% c("point", "mpoint")){
    monitor <- character()
  }else{
    monitor <- parameter_name
  }

  return(monitor)
}
.JAGS_monitor.vector         <- function(prior, parameter_name){

  monitor <- .JAGS_monitor.simple(prior, parameter_name)
  if(prior[["distribution"]] == "dirichlet"){
    monitor <- c(monitor, .JAGS_prior_dirichlet_eta_name(parameter_name))
  }

  return(monitor)
}
.JAGS_monitor.factor         <- function(prior, parameter_name){

  if(is.prior.ordered(prior)){
    return(.JAGS_monitor.ordered(prior, parameter_name))
  }

  monitor <- .JAGS_monitor.simple(prior, parameter_name)

  return(monitor)
}
.JAGS_monitor.ordered        <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.ordered(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  metadata <- .prior_ordered_metadata(prior)
  total_name <- .prior_ordered_total_name(parameter_name)

  monitor <- parameter_name

  if(is.prior.spike_and_slab(prior$total) && metadata$theta_dim > 1L){
    monitor <- c(
      monitor,
      paste0(total_name, "_indicator"),
      JAGS_to_monitor(setNames(list(.get_spike_and_slab_inclusion(prior$total)), paste0(total_name, "_inclusion"))),
      if(!is.prior.point(.get_spike_and_slab_variable(prior$total))) paste0(total_name, "_variable"),
      total_name
    )
  }else if(metadata$theta_dim == 1L){
    monitor <- c(monitor, JAGS_to_monitor(setNames(list(prior$total), total_name)))
  }else if(!is.prior.point(prior$total)){
    monitor <- c(monitor, total_name)
  }

  for(record in .prior_ordered_dirichlet_records(prior)){
    monitor <- c(monitor, .JAGS_prior_dirichlet_eta_name(record$node))
  }

  unique(monitor[nzchar(monitor)])
}
.JAGS_monitor.PP             <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    monitor <- .JAGS_monitor.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    monitor <- .JAGS_monitor.simple(prior, "PEESE")
  }

  return(monitor)
}
.JAGS_monitor.weightfunction <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  return(selection_backend_spec(prior)$monitor)
}
.JAGS_monitor_private.weightfunction <- function(prior){

  if(prior$weights$type == "cumulative"){
    return("eta")
  }
  if(prior$weights$type == "independent" && prior$weights$scale == "log_omega"){
    return("log_omega")
  }

  character()
}
.JAGS_monitor.phacking      <- function(prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)$monitor
}
.JAGS_monitor.bias          <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  selection_backend_spec(prior)$monitor
}
.JAGS_monitor.spike_and_slab <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  prior_variable  <- list(.get_spike_and_slab_variable(prior))
  prior_inclusion <- list(.get_spike_and_slab_inclusion(prior))
  names(prior_variable)  <- paste0(parameter_name, "_variable")
  names(prior_inclusion) <- paste0(parameter_name, "_inclusion")

  monitor <- c(
    paste0(parameter_name, "_indicator"),
    JAGS_to_monitor(prior_inclusion),
    parameter_name,
    JAGS_to_monitor(prior_variable)
  )

  return(monitor)
}
.JAGS_monitor.mixture        <- function(prior_list, parameter_name){

  .check_prior_list(prior_list, allow_expressions = TRUE)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_phacking       <- sapply(prior_list, is_prior_phacking)
    is_bias           <- sapply(prior_list, is_prior_bias)
    is_none           <- sapply(prior_list, is.prior.none)
    branch_info       <- lapply(prior_list, .selection_branch_info)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    monitor <- if(any(has_selection) || any(has_phacking)){
      selection_backend_spec(prior_list)$monitor
    }else{
      "bias_indicator"
    }

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      monitor <- c(monitor, "PET")
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      monitor <- c(monitor, "PEESE")
    }
  }else{

    monitor <- c(paste0(parameter_name, "_indicator"), parameter_name)
  }

  return(unique(monitor))
}
