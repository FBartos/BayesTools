#' @title Create initial values for 'JAGS' model
#'
#' @description Creates initial values for priors in
#' a 'JAGS' model.
#'
#' @param chains number of chains
#' @param seed seed for random number generation
#'
#' @inheritParams JAGS_add_priors
#'
#' @return \code{JAGS_add_priors} returns a list of JAGS
#' initial values.
#'
#' @export
JAGS_get_inits            <- function(prior_list, chains, seed){

  # return empty list in case that no prior was specified
  if(length(prior_list) == 0){
    return(list())
  }

  check_int(chains, "chains", lower = 1)
  check_real(seed, "seed", allow_NULL = TRUE)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)
  .bt_validate_ordered_shared_allocations(prior_list)


  # select seed at random if none was specified
  if(is.null(seed)){
    seed <- sample(666666, 1)
  }
  set.seed(seed)


  # create the starting values
  inits <- vector("list", chains)
  for(j in 1:chains){

    temp_inits <- .JAGS_get_inits.fun(prior_list)

    temp_inits[[".RNG.seed"]] <- seed + j
    temp_inits[[".RNG.name"]] <- if(chains > 4) "lecuyer::RngStream" else "base::Super-Duper"

    inits[[j]] <- temp_inits
  }

  return(inits)
}

.JAGS_get_inits.fun        <- function(prior_list){

  temp_inits <- list()
  ordered_allocation_keys <- character()

  for(i in seq_along(prior_list)){

    if(is.prior.point(prior_list[[i]])){

      next

    }else if(is.prior.weightfunction(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.weightfunction(prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.phacking(prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.bias(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.mixture(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.ordered(prior_list[[i]])){

      ordered_inits <- .JAGS_init.ordered(
        prior_list[[i]],
        names(prior_list)[i],
        emitted_allocations = ordered_allocation_keys
      )
      temp_inits <- c(temp_inits, ordered_inits[["inits"]])
      ordered_allocation_keys <- unique(c(
        ordered_allocation_keys,
        ordered_inits[["allocation_keys"]]
      ))

    }else if(is.prior.factor(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      temp_inits <- c(temp_inits, .JAGS_init.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(temp_inits)
}
.JAGS_init.simple          <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # no initialization for expression priors: require higher-level input
  if(.is_prior_expression(prior)){
    return()
  }

  if(prior[["distribution"]] == "point"){

    return()

  }else{
    init <- list()

    init[[parameter_name]] <- rng(prior, 1)
  }

  return(init)
}
.JAGS_positive_gamma_initialization <- function(shape, label){

  shape <- as.numeric(shape)
  if(length(shape) == 0L || any(!is.finite(shape)) || any(shape <= 0)){
    stop("Gamma initialization shapes must be finite and positive.",
         call. = FALSE)
  }

  values <- stats::rgamma(length(shape), shape = shape, rate = 1)
  invalid <- !is.finite(values) | values <= 0
  if(!any(invalid)){
    return(values)
  }

  warning(
    "RNG initialization failed for ", label,
    "; attempting a deterministic, order-one rescaling of the distribution medians.",
    call. = FALSE,
    immediate. = TRUE
  )
  medians <- stats::qgamma(0.5, shape = shape, rate = 1)
  if(all(medians == 0) && length(unique(shape)) == 1L){
    # Equal shapes have equal medians. Their common magnitude is immaterial to
    # the normalized simplex, so an order-one vector represents the exact
    # median proportions even when the common median underflows.
    medians[] <- 1
  }
  if(any(!is.finite(medians)) || any(medians <= 0)){
    stop(
      "RNG initialization failed for ", label,
      ", and the distribution medians are not representable as strictly ",
      "positive finite values.",
      call. = FALSE
    )
  }
  medians <- medians / max(medians)
  if(any(!is.finite(medians)) || any(medians <= 0)){
    stop(
      "The median fallback for ", label,
      " could not be rescaled to a strictly positive interior state.",
      call. = FALSE
    )
  }

  medians
}
.JAGS_binary_cumulative_initialization <- function(alpha, label){

  eta <- .JAGS_positive_gamma_initialization(
    shape = alpha,
    label = label
  )
  eta <- eta / max(eta)
  omega_ratio <- eta[2L] / sum(eta)
  if(!is.finite(omega_ratio) || omega_ratio <= 0 || omega_ratio >= 1){
    stop(
      "The initialization for ", label,
      " is not representable as a finite interior cumulative weight.",
      call. = FALSE
    )
  }

  omega_ratio
}

.JAGS_init.vector          <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # no initialization for expression priors: require higher-level input
  if(.is_prior_expression(prior)){
    return()
  }

  if(prior[["distribution"]] == "point"){

    return()

  }else{

    init <- list()


    if(prior[["distribution"]] == "dirichlet"){
      eta_init <- .JAGS_positive_gamma_initialization(
        shape = rep(
          prior$parameters[["alpha"]],
          length.out = prior$parameters[["K"]]
        ),
        label = paste0("Dirichlet prior '", parameter_name, "'")
      )
      init[[.JAGS_prior_dirichlet_eta_name(parameter_name)]] <- eta_init
    }else if(prior[["distribution"]] == "mt"){
      init[[paste0("prior_par_s_", parameter_name)]] <- rng(prior("gamma", list(shape = prior$parameters[["df"]]/2, rate = prior$parameters[["df"]]/2)), 1)
      init[[paste0("prior_par_z_", parameter_name)]] <- rng(prior("mnormal", list(mean = 0, sd = prior$parameters[["scale"]], K = prior$parameters[["K"]])), 1)[1,]
    }else{
      init[[parameter_name]] <- rng(prior, 1)[1,]
    }

  }

  return(init)
}
.JAGS_init.factor          <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(.get_prior_factor_levels(prior), "levels", lower = 1)

  # no initialization for expression priors: require higher-level input
  if(.is_prior_expression(prior)){
    return()
  }

  if(is.prior.ordered(prior)){

    init <- .JAGS_init.ordered(prior, parameter_name)[["inits"]]

  }else if(is.prior.treatment(prior) | is.prior.independent(prior)){

    init <- list()
    init[[parameter_name]] <- rng(prior, .get_prior_factor_levels(prior))

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)

    # remove the orthonormal/meandif class, otherwise samples from the transformed distributions are generated
    class(prior) <- class(prior)[!class(prior) %in% c("prior.orthonormal", "prior.meandif")]

    init <- .JAGS_init.vector(prior, parameter_name)

  }

  return(init)
}
.JAGS_init.ordered         <- function(prior, parameter_name, emitted_allocations = character()){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.ordered(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  metadata <- .prior_ordered_metadata(prior)
  total_name <- .prior_ordered_total_name(parameter_name)
  init <- .JAGS_init.ordered_total(prior$total, total_name, metadata$theta_dim)

  emitted_now <- character()
  for(record in .prior_ordered_dirichlet_records(prior)){
    if(record$key %in% emitted_allocations || record$key %in% emitted_now){
      next
    }
    eta_init <- .JAGS_positive_gamma_initialization(
      shape = rep(record$spec$alpha, length.out = record$dim),
      label = paste0("ordered allocation '", record$node, "'")
    )
    init[[.JAGS_prior_dirichlet_eta_name(record$node)]] <- eta_init
    emitted_now <- c(emitted_now, record$key)
  }

  list(inits = init, allocation_keys = emitted_now)
}

.JAGS_init.ordered_total   <- function(total, total_name, theta_dim){

  if(.is_prior_expression(total)){
    return(list())
  }

  if(is.prior.point(total)){
    return(list())
  }

  if(is.prior.spike_and_slab(total) && theta_dim > 1L){
    variable <- .get_spike_and_slab_variable(total)
    inclusion <- .get_spike_and_slab_inclusion(total)
    init <- list()
    if(!is.prior.point(variable)){
      init[[paste0(total_name, "_variable")]] <- rng(variable, theta_dim)
    }
    if(!is.prior.point(inclusion)){
      init[[paste0(total_name, "_inclusion")]] <- rng(inclusion, 1)
    }
    return(init)
  }

  if(theta_dim > 1L){
    init <- list()
    init[[total_name]] <- rng(total, theta_dim)
    return(init)
  }

  total_list <- list(total)
  names(total_list) <- total_name
  .JAGS_get_inits.fun(total_list)
}
.JAGS_init.PP              <- function(prior){

  .check_prior(prior)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    init <- .JAGS_init.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    init <- .JAGS_init.simple(prior, "PEESE")
  }

  return(init)
}
.JAGS_init.weightfunction  <- function(prior, component_id = NULL){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  if(is.null(component_id)){
    return(selection_backend_spec(prior)$init)
  }

  .selection_JAGS_init_weightfunction_component(prior, component_id)
}
.JAGS_init.phacking       <- function(prior, component_id = NULL){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  if(is.null(component_id)){
    return(selection_backend_spec(prior)$init)
  }

  alpha_name <- paste0("alpha_component_", component_id)
  init <- .JAGS_init.simple(prior$alpha, alpha_name)

  return(init)
}
.JAGS_init.bias           <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  return(selection_backend_spec(prior)$init)
}
.JAGS_init.spike_and_slab  <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")

  prior_variable        <- list(.get_spike_and_slab_variable(prior))
  names(prior_variable) <- paste0(parameter_name, "_variable")
  init <- .JAGS_get_inits.fun(prior_variable)

  if(!is.prior.point(.get_spike_and_slab_inclusion(prior))){
    init[[paste0(parameter_name, "_inclusion")]] <- rng(.get_spike_and_slab_inclusion(prior), 1)
  }


  return(init)
}
.JAGS_init.mixture         <- function(prior_list, parameter_name){

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

    init <- list()

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      named_prior_PET <- prior_list[[which(is_PET)]]
      class(named_prior_PET) <- class(named_prior_PET)[!class(named_prior_PET) %in% "prior.PET"]
      named_prior_PET <- list("PET_1" = named_prior_PET)

      init <- c(init, .JAGS_get_inits.fun(named_prior_PET))
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      named_prior_PEESE <- prior_list[[which(is_PEESE)]]
      class(named_prior_PEESE) <- class(named_prior_PEESE)[!class(named_prior_PEESE) %in% "prior.PEESE"]
      named_prior_PEESE <- list("PEESE_1" = named_prior_PEESE)

      init <- c(init, .JAGS_get_inits.fun(named_prior_PEESE))
    }
    if(any(has_selection) || any(has_phacking)){
      init <- c(init, selection_backend_spec(prior_list)$init)
    }else{
      init[["bias_indicator"]] <- rng(prior_list, 1, sample_components = TRUE)
    }

  }else{

    prior_components <- as.list(prior_list)
    class(prior_components) <- "list"
    names(prior_components) <- paste0(parameter_name, "_component_", seq_along(prior_components))

    init <- .JAGS_get_inits.fun(prior_components)
    init[[paste0(parameter_name, "_indicator")]] <- rng(prior_list, 1, sample_components = TRUE)

  }

  return(init)
}
