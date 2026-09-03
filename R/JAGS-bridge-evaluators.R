.bt_JAGS_bridge_compile_prior_list_evaluator <- function(prior_list,
                                                         emitted_allocations = character()){

  if(length(prior_list) == 0L){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list(),
      allocation_keys = unique(emitted_allocations)
    ))
  }

  if(!is.list(prior_list))
    stop("'prior_list' must be a list.")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)

  evaluators <- vector("list", length(prior_list))
  ordered_allocation_keys <- unique(emitted_allocations)
  for(i in seq_along(prior_list)){
    if(is.prior.ordered(prior_list[[i]])){
      evaluator <- .bt_JAGS_bridge_compile_ordered_evaluator(
        prior_list[[i]],
        names(prior_list)[i],
        emitted_allocations = ordered_allocation_keys
      )
      ordered_allocation_keys <- unique(c(
        ordered_allocation_keys,
        evaluator[["allocation_keys"]]
      ))
      evaluators[[i]] <- evaluator[["evaluator"]]
    }else{
      evaluators[[i]] <- .bt_JAGS_bridge_compile_prior_evaluator(
        prior_list[[i]],
        names(prior_list)[i]
      )
    }
  }

  list(
    log_prior = function(samples){
      marglik <- 0
      for(evaluator in evaluators){
        marglik <- marglik + evaluator$log_prior(samples)
      }
      marglik
    },
    parameters = function(samples){
      parameters <- list()
      for(evaluator in evaluators){
        parameters <- c(parameters, evaluator$parameters(samples))
      }
      parameters
    },
    allocation_keys = ordered_allocation_keys
  )
}

.bt_JAGS_bridge_compile_prior_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  if(is.prior.weightfunction(prior_object)){
    return(.bt_JAGS_bridge_compile_weightfunction_evaluator(prior_object))
  }else if(is.prior.none(prior_object)){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list()
    ))
  }else if(is_prior_phacking(prior_object)){
    return(.bt_JAGS_bridge_compile_phacking_evaluator(prior_object))
  }else if(is_prior_bias(prior_object)){
    return(.bt_JAGS_bridge_compile_bias_evaluator(prior_object))
  }else if(is.prior.mixture(prior_object)){
    .JAGS_marglik_stop_unsupported_mixture(prior_object)
  }else if(is.prior.PET(prior_object) | is.prior.PEESE(prior_object)){
    return(.bt_JAGS_bridge_compile_PP_evaluator(prior_object))
  }else if(is.prior.ordered(prior_object)){
    return(.bt_JAGS_bridge_compile_ordered_evaluator(prior_object, parameter_name)[["evaluator"]])
  }else if(is.prior.factor(prior_object)){
    return(.bt_JAGS_bridge_compile_factor_evaluator(prior_object, parameter_name))
  }else if(is.prior.vector(prior_object)){
    return(.bt_JAGS_bridge_compile_vector_evaluator(prior_object, parameter_name))
  }else if(is.prior.simple(prior_object)){
    return(.bt_JAGS_bridge_compile_simple_evaluator(prior_object, parameter_name))
  }

  stop("Unsupported prior object.", call. = FALSE)
}

.bt_JAGS_bridge_compile_ordered_evaluator <- function(prior_object, parameter_name, emitted_allocations = character()){

  force(prior_object)
  force(parameter_name)
  force(emitted_allocations)

  .prior_ordered_bridge_check(prior_object)
  total_names <- .prior_ordered_total_monitor_names(prior_object, parameter_name)
  total_prior <- prior_object$total

  emitted_now <- character()
  dirichlet_records <- .prior_ordered_dirichlet_records(prior_object)
  dirichlet_records <- dirichlet_records[
    !vapply(dirichlet_records, function(record){
      record$key %in% emitted_allocations
    }, logical(1))
  ]

  for(record in dirichlet_records){
    emitted_now <- c(emitted_now, record$key)
  }

  evaluator <- list(
    log_prior = function(samples){
      marglik <- 0
      if(!is.prior.point(total_prior)){
        total_values <- unname(unlist(samples[total_names], use.names = FALSE))
        marglik <- marglik + sum(lpdf(total_prior, total_values))
      }
      for(record in dirichlet_records){
        eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(record$node), "[", seq_len(record$dim), "]")
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored ordered Dirichlet allocation parameters."
        )
        if(is.null(eta)){
          return(-Inf)
        }
        marglik <- marglik + sum(stats::dgamma(eta, shape = record$spec$alpha, rate = 1, log = TRUE))
      }
      marglik
    },
    parameters = function(samples){
      parameter_names <- .JAGS_prior_factor_names(parameter_name, prior_object)
      parameter <- list()
      parameter[[parameter_name]] <- .bt_JAGS_bridge_compile_parameter_values(
        prior_object,
        parameter_names
      )(samples)
      parameter
    }
  )

  list(evaluator = evaluator, allocation_keys = emitted_now)
}

.bt_JAGS_bridge_compile_simple_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  distribution <- prior_object[["distribution"]]
  log_density  <- .prior_simple_lpdf_evaluator(prior_object)

  if(identical(distribution, "invgamma")){
    return(list(
      log_prior = function(samples){
        value <- .bt_JAGS_marglik_invgamma_values(
          samples = samples,
          parameter_names = parameter_name,
          missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters."
        )
        if(is.null(value)){
          return(-Inf)
        }
        log_density(value)
      },
      parameters = function(samples){
        value <- .bt_JAGS_marglik_invgamma_values(
          samples = samples,
          parameter_names = parameter_name,
          missing_message = "'samples' does not contain all monitored inverse-gamma prior parameters.",
          signal = TRUE
        )
        parameter <- list()
        parameter[[parameter_name]] <- value
        parameter
      }
    ))
  }

  if(identical(distribution, "point")){
    location <- prior_object$parameters[["location"]]
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples){
        parameter <- list()
        parameter[[parameter_name]] <- location
        parameter
      }
    ))
  }

  list(
    log_prior = function(samples){
      log_density(samples[[parameter_name]])
    },
    parameters = function(samples){
      parameter <- list()
      parameter[[parameter_name]] <- samples[[parameter_name]]
      parameter
    }
  )
}

.bt_JAGS_bridge_compile_parameter_values <- function(prior_object, parameter_names){

  force(prior_object)
  force(parameter_names)

  if(is.prior.ordered(prior_object)){
    return(.bt_JAGS_marglik_compile_ordered_parameter_values(
      prior = prior_object,
      parameter_names = parameter_names
    ))
  }

  if(is.prior.point(prior_object)){
    location <- prior_object$parameters[["location"]]
    return(function(samples) rep(location, length(parameter_names)))
  }

  if(identical(prior_object[["distribution"]], "invgamma")){
    return(function(samples){
      .bt_JAGS_marglik_invgamma_values(
        samples = samples,
        parameter_names = parameter_names,
        missing_message = "'samples' does not contain all monitored formula prior parameters.",
        signal = TRUE
      )
    })
  }

  sample_names <- parameter_names
  function(samples){
    if(!all(sample_names %in% names(samples))){
      stop("'samples' does not contain all monitored formula prior parameters.", call. = FALSE)
    }
    unname(unlist(samples[sample_names], use.names = FALSE))
  }
}

.bt_JAGS_bridge_compile_vector_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  if(identical(prior_object[["distribution"]], "dirichlet")){
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(prior_object$parameters[["K"]]), "]")
    alpha <- prior_object$parameters[["alpha"]]
    return(list(
      log_prior = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored Dirichlet prior parameters."
        )
        if(is.null(eta)){
          return(-Inf)
        }
        sum(stats::dgamma(eta, shape = alpha, rate = 1, log = TRUE))
      },
      parameters = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored Dirichlet prior parameters.",
          signal = TRUE
        )
        parameter <- list()
        parameter[[parameter_name]] <- eta / sum(eta)
        parameter
      }
    ))
  }

  if(prior_object$parameters[["K"]] == 1){
    parameter_monitor_names <- parameter_name
  }else{
    parameter_monitor_names <- paste0(parameter_name, "[", seq_len(prior_object$parameters[["K"]]), "]")
  }

  if(identical(prior_object[["distribution"]], "mpoint")){
    location <- prior_object$parameters[["location"]]
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples){
        parameter <- list()
        parameter[[parameter_name]] <- rep(location, length(parameter_monitor_names))
        parameter
      }
    ))
  }

  list(
    log_prior = function(samples){
      if(length(parameter_monitor_names) == 1L){
        lpdf(prior_object, samples[[parameter_monitor_names]])
      }else{
        lpdf(prior_object, samples[parameter_monitor_names])
      }
    },
    parameters = function(samples){
      parameter <- list()
      parameter[[parameter_name]] <- samples[parameter_monitor_names]
      parameter
    }
  )
}

.bt_JAGS_bridge_compile_factor_evaluator <- function(prior_object, parameter_name){

  force(prior_object)
  force(parameter_name)

  if(is.prior.treatment(prior_object) | is.prior.independent(prior_object)){
    levels <- .get_prior_factor_levels(prior_object)
    if(levels == 1L){
      return(.bt_JAGS_bridge_compile_simple_evaluator(prior_object, parameter_name))
    }

    parameter_names <- paste0(parameter_name, "[", seq_len(levels), "]")
    density_evaluators <- lapply(parameter_names, function(name){
      .bt_JAGS_bridge_compile_simple_evaluator(prior_object, name)$log_prior
    })
    parameter_values <- .bt_JAGS_bridge_compile_parameter_values(prior_object, parameter_names)
    return(list(
      log_prior = function(samples){
        marglik <- 0
        for(evaluator in density_evaluators){
          marglik <- marglik + evaluator(samples)
        }
        marglik
      },
      parameters = function(samples){
        parameter <- list()
        parameter[[parameter_name]] <- parameter_values(samples)
        parameter
      }
    ))
  }

  factor_prior <- prior_object
  factor_prior$parameters[["K"]] <- .get_prior_factor_levels(prior_object)
  .bt_JAGS_bridge_compile_vector_evaluator(factor_prior, parameter_name)
}

.bt_JAGS_bridge_compile_PP_evaluator <- function(prior_object){

  force(prior_object)

  if(is.prior.PET(prior_object)){
    return(.bt_JAGS_bridge_compile_simple_evaluator(prior_object, "PET"))
  }
  .bt_JAGS_bridge_compile_simple_evaluator(prior_object, "PEESE")
}

.bt_JAGS_bridge_compile_weightfunction_evaluator <- function(prior_object){

  force(prior_object)

  J <- .weightfunction_n_bins(prior_object)
  expansion <- .weightfunction_mapping_expansion(prior_object, force_one_sided = TRUE)

  if(prior_object$weights$type == "fixed"){
    omega_fixed <- unname(prior_object$weights$omega[expansion$index])
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list(omega = omega_fixed)
    ))
  }

  if(prior_object$weights$type == "cumulative"){
    if(J == 2L){
      beta_parameters <- .weightfunction_alpha_marginal(
        prior_object$weights$alpha,
        2L
      )
      return(list(
        log_prior = function(samples){
          omega <- .bt_JAGS_marglik_binary_cumulative_weight(samples)
          if(is.null(omega)){
            return(-Inf)
          }
          stats::dbeta(
            omega,
            shape1 = beta_parameters$alpha,
            shape2 = beta_parameters$beta,
            log = TRUE
          )
        },
        parameters = function(samples){
          omega <- c(
            1,
            .bt_JAGS_marglik_binary_cumulative_weight(samples, signal = TRUE)
          )
          list(omega = unname(omega[expansion$index]))
        }
      ))
    }

    eta_names <- paste0("eta[", seq_len(J), "]")
    alpha <- prior_object$weights$alpha
    return(list(
      log_prior = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored cumulative weightfunction parameters."
        )
        if(is.null(eta)){
          return(-Inf)
        }
        sum(stats::dgamma(eta, shape = alpha, rate = 1, log = TRUE))
      },
      parameters = function(samples){
        eta <- .bt_JAGS_marglik_positive_auxiliary_values(
          samples = samples,
          parameter_names = eta_names,
          missing_message = "'samples' does not contain all monitored cumulative weightfunction parameters.",
          signal = TRUE
        )
        std_eta <- eta / sum(eta)
        omega <- unname(rev(cumsum(rev(std_eta))))
        list(omega = unname(omega[expansion$index]))
      }
    ))
  }

  omega <- rep(1, J)
  if(J == 1L){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list(omega = unname(omega[expansion$index]))
    ))
  }

  if(prior_object$weights$scale == "omega"){
    omega_names <- paste0("omega[", 2:J, "]")
    weight_prior <- prior_object$weights$prior
    return(list(
      log_prior = function(samples){
        sum(mlpdf(weight_prior, samples[omega_names]))
      },
      parameters = function(samples){
        omega[2:J] <- samples[omega_names]
        list(omega = unname(omega[expansion$index]))
      }
    ))
  }

  log_omega_names <- paste0("log_omega[", 2:J, "]")
  weight_prior <- prior_object$weights$prior
  list(
    log_prior = function(samples){
      sum(mlpdf(weight_prior, samples[log_omega_names]))
    },
    parameters = function(samples){
      omega[2:J] <- exp(samples[log_omega_names])
      list(omega = unname(omega[expansion$index]))
    }
  )
}

.bt_JAGS_bridge_compile_phacking_evaluator <- function(prior_object){

  force(prior_object)

  alpha_evaluator <- .bt_JAGS_bridge_compile_simple_evaluator(prior_object$alpha, "alpha")
  constants <- phack_backend_constants(
    prior_object$form,
    prior_object$source,
    prior_object$destination,
    target = prior_object$target
  )

  list(
    log_prior = alpha_evaluator$log_prior,
    parameters = function(samples){
      alpha <- alpha_evaluator$parameters(samples)[["alpha"]]
      list(
        alpha     = alpha,
        pi_null   = alpha * constants$pi_null_per_alpha,
        beta_null = alpha * constants$beta_null_per_alpha
      )
    }
  )
}

.bt_JAGS_bridge_compile_bias_evaluator <- function(prior_object){

  force(prior_object)

  selection_backend_spec(prior_object)

  selection_evaluator <- if(!is.null(prior_object$selection)){
    .bt_JAGS_bridge_compile_weightfunction_evaluator(prior_object$selection)
  }else{
    NULL
  }
  phacking_evaluator <- if(!is.null(prior_object$phacking)){
    .bt_JAGS_bridge_compile_phacking_evaluator(prior_object$phacking)
  }else{
    NULL
  }

  list(
    log_prior = function(samples){
      marglik <- 0
      if(!is.null(selection_evaluator)){
        marglik <- marglik + selection_evaluator$log_prior(samples)
      }
      if(!is.null(phacking_evaluator)){
        marglik <- marglik + phacking_evaluator$log_prior(samples)
      }
      marglik
    },
    parameters = function(samples){
      parameter <- list()
      if(!is.null(selection_evaluator)){
        parameter <- c(parameter, selection_evaluator$parameters(samples))
      }
      if(!is.null(phacking_evaluator)){
        parameter <- c(parameter, phacking_evaluator$parameters(samples))
      }
      parameter
    }
  )
}

.bt_JAGS_bridge_compile_formula_prior_evaluator <- function(formula_prior_list,
                                                            emitted_allocations = character()){

  if(length(formula_prior_list) == 0L){
    return(list(
      log_prior = function(samples) 0,
      parameters = function(samples) list(),
      allocation_keys = unique(emitted_allocations)
    ))
  }

  prior_list <- do.call(c, unname(formula_prior_list))
  .bt_JAGS_bridge_compile_prior_list_evaluator(
    prior_list = prior_list,
    emitted_allocations = emitted_allocations
  )
}

.bt_JAGS_bridge_compile_model_prior_evaluators <- function(prior_list,
                                                           formula_prior_list){

  prior_evaluator <- .bt_JAGS_bridge_compile_prior_list_evaluator(prior_list)
  formula_prior_evaluator <- .bt_JAGS_bridge_compile_formula_prior_evaluator(
    formula_prior_list = formula_prior_list,
    emitted_allocations = prior_evaluator$allocation_keys
  )

  list(
    prior = prior_evaluator,
    formula = formula_prior_evaluator
  )
}

.bt_JAGS_bridge_compile_formula_random_prior_evaluator <- function(
    formula_design_list,
    omitted_latent = character()){

  .bt_JAGS_bridge_check_no_allocation_inclusion(formula_design_list)

  evaluators <- list()
  if(length(formula_design_list) > 0L){
    design_names <- names(formula_design_list)
    if(is.null(design_names)){
      design_names <- rep("", length(formula_design_list))
    }
    for(design_i in seq_along(formula_design_list)){
      design <- formula_design_list[[design_i]]
      if(.bt_formula_design_has_any_random_effects(design)){
        .bt_JAGS_bridge_validate_formula_random_compile(
          parameter = .bt_JAGS_bridge_design_parameter_name(
            design,
            fallback = design_names[[design_i]]
          ),
          design = design,
          label = "stored"
        )
        for(random_term in .bt_formula_design_random_effects(design)){
          evaluators[[length(evaluators) + 1L]] <-
            .bt_JAGS_bridge_compile_random_effect_prior_evaluator(
              random_term = random_term,
              omitted_latent = omitted_latent
            )
        }
      }
    }
  }

  list(
    uses_posterior_row = length(evaluators) > 0L,
    log_prior = function(samples){
      marglik <- 0
      for(evaluator in evaluators){
        contribution <- evaluator$log_prior(samples)
        if(is.null(contribution)){
          return(-Inf)
        }
        marglik <- marglik + contribution
      }
      marglik
    }
  )
}

.bt_JAGS_bridge_compile_random_effect_prior_evaluator <- function(
    random_term,
    omitted_latent = character()){

  force(random_term)
  force(omitted_latent)

  n_columns <- random_term$n_columns
  sampled_random_effect <- identical(
    .bt_random_effect_term_compile_mode(random_term),
    "sampled"
  )
  latent_evaluator <- .bt_JAGS_bridge_compile_random_effect_latent_log_density(
    random_term = random_term,
    omitted_latent = omitted_latent
  )
  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  scalar_rho_support_evaluator <- .bt_JAGS_bridge_compile_random_effect_scalar_rho_support(
    random_term = random_term,
    structure = structure
  )
  lkj_evaluator <- .bt_JAGS_bridge_compile_random_effect_lkj_prior(
    random_term = random_term,
    structure = structure,
    n_columns = n_columns
  )

  list(
    log_prior = function(samples){
      marglik <- 0
      if(isTRUE(sampled_random_effect)){
        marglik <- latent_evaluator(samples)
      }

      scalar_rho_support <- scalar_rho_support_evaluator(samples)
      if(!is.finite(scalar_rho_support)){
        return(scalar_rho_support)
      }
      marglik <- marglik + scalar_rho_support

      if(!is.null(lkj_evaluator)){
        marglik <- marglik + lkj_evaluator(samples)
        if(is.na(marglik)){
          return(-Inf)
        }
      }

      marglik
    }
  )
}

.bt_JAGS_bridge_compile_random_effect_latent_log_density <- function(
    random_term,
    omitted_latent = character()){

  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  active_names <- setdiff(z_names, omitted_latent)
  if(length(active_names) == 0L){
    return(function(samples) 0)
  }

  complete <- length(active_names) == length(z_names)
  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  known_group <- .bt_random_effect_has_known_group_covariance(random_term)
  if(!complete && (!structure %in% c("diag", "id") || known_group)){
    stop(
      "Bridge sampling can omit only a complete correlated random-effect ",
      "latent block or independent fixed-zero latent components.",
      call. = FALSE
    )
  }
  group_covariance <- if(known_group){
    .bt_random_effect_known_group_covariance(
      random_term,
      context = "Bridge sampling"
    )
  }else{
    NULL
  }
  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns

  function(samples){
    if(!all(active_names %in% names(samples))){
      stop(
        "Bridge samples are missing standardized latent random effects for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    z_values <- samples[active_names]
    if(any(is.na(z_values))){
      return(-Inf)
    }
    if(!is.null(group_covariance)){
      z_values <- matrix(
        as.numeric(z_values),
        nrow = n_groups,
        ncol = n_columns
      )
      out <- 0
      for(column in seq_len(n_columns)){
        out <- out + .bt_mvn_zero_log_density(
          z = z_values[, column],
          precision = group_covariance$precision,
          log_det = group_covariance$log_det
        )
      }
      return(out)
    }

    out <- sum(stats::dnorm(z_values, mean = 0, sd = 1, log = TRUE))
    if(is.na(out)) -Inf else out
  }
}

.bt_JAGS_bridge_compile_random_effect_scalar_rho_support <- function(random_term,
                                                                    structure){

  force(random_term)
  force(structure)

  if(!structure %in% c("cs", "hcs", "ar1", "car", "har") ||
     random_term$n_columns <= 1L){
    return(function(samples) 0)
  }

  support_spec <- .bt_JAGS_random_effect_scalar_rho_support_spec(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  rho_evaluator <- .bt_random_effect_compile_rho_draw_evaluator(
    random_term = random_term,
    missing = "error",
    out_of_support = "null",
    context = "Bridge sampling random-effect metadata"
  )
  force(rho_evaluator)

  function(samples){
    rho <- rho_evaluator(
      .bt_JAGS_marglik_random_effect_posterior_row(samples)
    )
    if(is.null(rho) || any(is.na(rho) | !is.finite(rho))){
      return(-Inf)
    }
    if(any(.bt_random_effect_rho_outside_support(
      rho,
      bounds = support_spec$bounds,
      structure = structure
    ))){
      return(-Inf)
    }

    0
  }
}

.bt_JAGS_bridge_compile_random_effect_lkj_prior <- function(random_term,
                                                            structure,
                                                            n_columns){

  force(random_term)
  force(structure)
  force(n_columns)

  if(!identical(structure, "us") || n_columns <= 1L){
    return(NULL)
  }

  block_name <- random_term$block_name
  u_names <- .bt_random_effect_lkj_primitive_names(
    random_term,
    n_columns,
    context = "Bridge sampling random-effect metadata"
  )
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = "us",
    context = "Bridge sampling random-effect metadata"
  )
  eta <- correlation$eta
  if(!is.numeric(eta) || length(eta) != 1L || is.na(eta)){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical 'random_term$correlation$eta'.",
      call. = FALSE
    )
  }

  function(samples){
    if(!all(u_names %in% names(samples))){
      stop(
        "Bridge samples are missing LKJ primitive coordinates for block '",
        block_name,
        "'.",
        call. = FALSE
      )
    }
    u_values <- unname(samples[u_names])
    if(any(is.na(u_values) | u_values <= 0 | u_values >= 1)){
      return(-Inf)
    }

    .bt_lkj_cholesky_cpc_u_log_prior(
      u_values,
      K = n_columns,
      eta = eta
    )
  }
}
