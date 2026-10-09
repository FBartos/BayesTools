bayestools_reference_formula_random_log_prior <- function(samples,
                                                          formula_design_list){

  if(length(formula_design_list) == 0L){
    return(0)
  }
  BayesTools:::.bt_JAGS_bridge_check_no_allocation_inclusion(
    formula_design_list
  )

  marglik <- 0
  design_names <- names(formula_design_list)
  if(is.null(design_names)){
    design_names <- rep("", length(formula_design_list))
  }
  for(design_i in seq_along(formula_design_list)){
    design <- formula_design_list[[design_i]]
    if(!BayesTools:::.bt_formula_design_has_any_random_effects(design)){
      next
    }
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_compile(
      parameter = BayesTools:::.bt_JAGS_bridge_design_parameter_name(
        design,
        fallback = design_names[[design_i]]
      ),
      design = design,
      label = "stored"
    )
    random_effects <- BayesTools:::.bt_formula_design_random_effects(design)
    for(random_term in random_effects){
      contribution <- bayestools_reference_random_effect_log_prior(
        samples,
        random_term
      )
      if(is.na(contribution)){
        return(-Inf)
      }
      marglik <- marglik + contribution
      if(is.na(marglik)){
        return(-Inf)
      }
    }
  }

  marglik
}

bayestools_reference_random_effect_log_prior <- function(samples,
                                                         random_term){

  n_columns <- random_term$n_columns
  sampled_random_effect <- identical(
    BayesTools:::.bt_random_effect_term_compile_mode(random_term),
    "sampled"
  )
  marglik <- 0
  if(isTRUE(sampled_random_effect)){
    marglik <- bayestools_reference_random_effect_latent_log_density(
      random_term,
      samples
    )
  }

  scalar_rho_support <- bayestools_reference_random_effect_scalar_rho_support(
    samples = samples,
    random_term = random_term
  )
  if(!is.finite(scalar_rho_support)){
    return(scalar_rho_support)
  }
  marglik <- marglik + scalar_rho_support

  structure <- BayesTools:::.bt_JAGS_bridge_random_term_structure(random_term)
  if(identical(structure, "us") && n_columns > 1L){
    u_names <- BayesTools:::.bt_random_effect_lkj_primitive_names(
      random_term,
      n_columns,
      context = "Bridge sampling random-effect metadata"
    )
    if(!all(u_names %in% names(samples))){
      stop(
        "Bridge samples are missing LKJ primitive coordinates for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    u_values <- unname(samples[u_names])
    if(any(is.na(u_values) | u_values <= 0 | u_values >= 1)){
      return(-Inf)
    }
    correlation <- BayesTools:::.bt_random_effect_correlation_metadata(
      random_term,
      structure = "us",
      context = "Bridge sampling random-effect metadata"
    )
    eta <- correlation$eta
    if(!is.numeric(eta) || length(eta) != 1L || is.na(eta)){
      stop(
        "Bridge sampling random-effect metadata",
        BayesTools:::.bt_random_effect_metadata_block_detail(random_term),
        " is missing canonical 'random_term$correlation$eta'.",
        call. = FALSE
      )
    }
    marglik <- marglik + BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(
      u_values,
      K = n_columns,
      eta = eta
    )
    if(is.na(marglik)){
      return(-Inf)
    }
  }

  marglik
}

bayestools_reference_random_effect_scalar_rho_support <- function(samples,
                                                                  random_term){

  structure <- BayesTools:::.bt_JAGS_bridge_random_term_structure(random_term)
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har") ||
     random_term$n_columns <= 1L){
    return(0)
  }

  support_spec <- BayesTools:::.bt_JAGS_random_effect_scalar_rho_support_spec(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )

  posterior <- BayesTools:::.bt_JAGS_marglik_random_effect_posterior_row(samples)
  correlation <- support_spec$correlation
  rho_scale <- BayesTools:::.bt_random_effect_rho_scale_metadata(
    correlation,
    random_term = random_term,
    context = "Bridge sampling random-effect metadata"
  )
  source_name <- if(!identical(rho_scale, "rho") &&
                    correlation$sample_name %in% colnames(posterior)){
    correlation$sample_name
  }else if(correlation$rho_name %in% colnames(posterior)){
    correlation$rho_name
  }else{
    correlation$sample_name
  }
  if(source_name %in% colnames(posterior) &&
     any(!is.finite(posterior[, source_name]))){
    return(-Inf)
  }

  rho <- BayesTools:::.bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "error",
    out_of_support = "null",
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(rho) || any(is.na(rho) | !is.finite(rho))){
    return(-Inf)
  }
  if(any(BayesTools:::.bt_random_effect_rho_outside_support(
    rho,
    bounds = support_spec$bounds,
    structure = structure
  ))){
    return(-Inf)
  }

  0
}

bayestools_reference_random_effect_latent_log_density <- function(random_term,
                                                                  samples){

  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns
  z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  ))
  if(!all(z_names %in% names(samples))){
    stop(
      "Bridge samples are missing standardized latent random effects for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  z_values <- samples[z_names]
  if(any(is.na(z_values))){
    return(-Inf)
  }

  if(BayesTools:::.bt_random_effect_has_known_group_covariance(random_term)){
    group_covariance <- BayesTools:::.bt_random_effect_known_group_covariance(
      random_term,
      context = "Bridge sampling"
    )
    z_values <- matrix(
      as.numeric(z_values),
      nrow = n_groups,
      ncol = n_columns
    )
    out <- 0
    for(column in seq_len(n_columns)){
      out <- out + BayesTools:::.bt_mvn_zero_log_density(
        z = z_values[, column],
        precision = group_covariance$precision,
        log_det = group_covariance$log_det
      )
    }
    return(out)
  }

  marglik <- sum(stats::dnorm(z_values, mean = 0, sd = 1, log = TRUE))
  if(is.na(marglik)){
    return(-Inf)
  }
  marglik
}
