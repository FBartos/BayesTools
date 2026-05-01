#' @title Helper functions for extracting and formatting posterior distributions
#'
#' @description Internal helper functions to extract posterior samples from JAGS
#' fits and reformat them for further processing (summary tables, diagnostics, plots).
#' These functions consolidate common logic that was duplicated across
#' \code{runjags_estimates_table}, \code{.diagnostics_plot_data}, and plotting functions.
#'
#' @name posterior_extraction_helpers
#' @keywords internal
NULL


#' @rdname posterior_extraction_helpers
#' @param fit a JAGS model fit object
#' @param as_list whether to return samples as mcmc.list (TRUE) or merged matrix (FALSE)
#' @return matrix or mcmc.list of posterior samples
.extract_posterior_samples <- function(fit, as_list = FALSE) {
  
  if (as_list) {
    # Use generic function to allow S3 method dispatch (runjags has its own as.mcmc.list method)
    model_samples <- coda::as.mcmc.list(fit)
  } else {
    # Use generic function to allow S3 method dispatch (runjags has its own as.mcmc method)
    model_samples <- suppressWarnings(coda::as.mcmc(fit))
  }
  
  return(model_samples)
}


#' @rdname posterior_extraction_helpers
#' @param model_samples matrix of posterior samples
#' @param prior_list list of prior objects
#' @param remove_parameters character vector of parameter names to remove
#' @return list with cleaned model_samples and updated prior_list
.remove_auxiliary_parameters <- function(model_samples, prior_list, remove_parameters = NULL) {
  
  for (i in rev(seq_along(prior_list))) {
    
    par_name <- names(prior_list)[i]
    
    # invgamma support parameter
    if (is.prior.simple(prior_list[[i]]) && prior_list[[i]][["distribution"]] == "invgamma") {
      model_samples <- model_samples[, colnames(model_samples) != paste0("inv_", par_name), drop = FALSE]
    }
    
    # weightfunction parameters
    if (is.prior.weightfunction(prior_list[[i]])) {
      private_parameters <- .JAGS_monitor_private.weightfunction(prior_list[[i]])
      if(length(private_parameters) > 0){
        private_pattern <- paste0("^(", paste(private_parameters, collapse = "|"), ")(\\[|$)")
        model_samples <- model_samples[, !grepl(private_pattern, colnames(model_samples)), drop = FALSE]
      }
      
      # rename the omegas
      omega_cuts      <- weightfunctions_mapping(prior_list[i], cuts_only = TRUE)
      omega_names_old <- paste0("omega[", 1:(length(omega_cuts) - 1), "]")
      omega_names     <- sapply(1:(length(omega_cuts) - 1), function(j) paste0("omega[", omega_cuts[j], ",", omega_cuts[j + 1], "]"))
      omega_raw_names <- grep("^omega\\[[0-9]+\\]$", colnames(model_samples), value = TRUE)

      colnames(model_samples)[which(colnames(model_samples) %in% omega_names_old)] <- omega_names

      # Two-sided weightfunctions are expanded internally to mirrored one-sided
      # coefficients for JAGS. Only the public local bins belong in summaries.
      omega_extra_names <- setdiff(omega_raw_names, omega_names_old)
      if(length(omega_extra_names) > 0L){
        model_samples <- model_samples[, !colnames(model_samples) %in% omega_extra_names, drop = FALSE]
      }

      # remove omegas if requested
      if ("omega" %in% remove_parameters) {
        model_samples <- model_samples[, !colnames(model_samples) %in% omega_names, drop = FALSE]
        prior_list[i] <- NULL
      }

    } else if (is_prior_phacking(prior_list[[i]])) {
      if ("omega" %in% remove_parameters) {
        omega_names <- colnames(model_samples)[grepl("^omega\\[", colnames(model_samples))]
        model_samples <- model_samples[, !colnames(model_samples) %in% omega_names, drop = FALSE]
      }

      if (par_name %in% remove_parameters) {
        model_samples <- .remove_parameter_columns(model_samples, prior_list[[i]], par_name)
        prior_list[i] <- NULL
      }

    } else if (is_prior_bias(prior_list[[i]])) {
      if (.selection_prior_has_selection(prior_list[[i]])) {
        omega_cuts      <- weightfunctions_mapping(.selection_prior_selection_priors(prior_list[[i]]), cuts_only = TRUE, one_sided = TRUE)
        omega_names_old <- paste0("omega[", 1:(length(omega_cuts) - 1), "]")
        omega_names     <- sapply(1:(length(omega_cuts) - 1), function(j) paste0("omega[", omega_cuts[j], ",", omega_cuts[j + 1], "]"))
        colnames(model_samples)[which(colnames(model_samples) %in% omega_names_old)] <- omega_names
      } else {
        omega_names <- colnames(model_samples)[grepl("^omega\\[", colnames(model_samples))]
      }

      if ("omega" %in% remove_parameters) {
        model_samples <- model_samples[, !colnames(model_samples) %in% omega_names, drop = FALSE]
      }

      if (par_name %in% remove_parameters) {
        model_samples <- .remove_parameter_columns(model_samples, prior_list[[i]], par_name)
        prior_list[i] <- NULL
      }
      
    } else if (par_name %in% remove_parameters) {
      # remove parameters to be excluded (note: spike_0 removal is handled by caller)
      model_samples <- .remove_parameter_columns(model_samples, prior_list[[i]], par_name)
      prior_list[i] <- NULL
    }
  }
  
  return(list(model_samples = model_samples, prior_list = prior_list))
}


#' @rdname posterior_extraction_helpers
#' @description Helper to remove all columns associated with a parameter
#' @param model_samples matrix of posterior samples
#' @param prior prior object for the parameter
#' @param par_name name of the parameter
#' @return updated model_samples matrix
.remove_parameter_columns <- function(model_samples, prior, par_name) {
  
  # collect all column patterns to remove
  cols_to_remove <- character(0)
  
  if (is.prior.spike_and_slab(prior)) {
    # spike and slab: remove main parameter, indicator, inclusion, variable
    cols_to_remove <- c(
      par_name,
      paste0(par_name, "_indicator"),
      paste0(par_name, "_inclusion"),
      paste0(par_name, "_variable")
    )
    # also handle factor spike and slab with indexed columns
    cols_to_remove <- c(cols_to_remove, 
      colnames(model_samples)[grepl(paste0("^", par_name, "\\["), colnames(model_samples))],
      colnames(model_samples)[grepl(paste0("^", par_name, "_variable\\["), colnames(model_samples))]
    )
    
  } else if (is.prior.mixture(prior)) {
    # mixture: remove main parameter, indicator, and component-specific columns
    cols_to_remove <- c(
      par_name,
      paste0(par_name, "_indicator")
    )
    # handle factor mixture with indexed columns
    cols_to_remove <- c(cols_to_remove,
      colnames(model_samples)[grepl(paste0("^", par_name, "\\["), colnames(model_samples))]
    )
    
    # check for bias mixture (PET, PEESE, omega, p-hacking)
    if (inherits(prior, "prior.bias_mixture")) {
      cols_to_remove <- c(cols_to_remove, .selection_bias_parameter_names(prior))
      cols_to_remove <- c(cols_to_remove,
        colnames(model_samples)[grepl("^omega\\[", colnames(model_samples))]
      )
    }
    
  } else if (is_prior_phacking(prior)) {
    cols_to_remove <- c(par_name, "omega", "alpha", "pi_null", "phack_kind")
    cols_to_remove <- c(cols_to_remove,
      colnames(model_samples)[grepl("^omega\\[", colnames(model_samples))]
    )

  } else if (is_prior_bias(prior)) {
    cols_to_remove <- c(par_name, .selection_bias_parameter_names(prior))
    cols_to_remove <- c(cols_to_remove,
      colnames(model_samples)[grepl("^omega\\[", colnames(model_samples))]
    )

  } else if (is.prior.factor(prior)) {
    # factor prior: remove all indexed columns
    cols_to_remove <- .JAGS_prior_factor_names(par_name, prior)
    
  } else if (is.prior.PET(prior)) {
    # PET prior: remove the PET column (samples are stored as "PET", not par_name)
    cols_to_remove <- c(par_name, "PET")
    
  } else if (is.prior.PEESE(prior)) {
    # PEESE prior: remove the PEESE column (samples are stored as "PEESE", not par_name)
    cols_to_remove <- c(par_name, "PEESE")
    
  } else {
    # simple prior: just remove the main column
    cols_to_remove <- par_name
  }
  
  # remove duplicates and filter to existing columns
  cols_to_remove <- unique(cols_to_remove)
  cols_to_remove <- cols_to_remove[cols_to_remove %in% colnames(model_samples)]
  
  model_samples <- model_samples[, !colnames(model_samples) %in% cols_to_remove, drop = FALSE]
  
  return(model_samples)
}


#' @rdname posterior_extraction_helpers
#' @param remove_parameters character vector of parameter names to remove, or TRUE to remove all non-formula parameters.
#' If "bias" is specified and the bias prior contains PET, PEESE, or weightfunction priors,
#' the corresponding parameters (PET, PEESE, omega) are also added to the removal list.
#' @param remove_formulas character vector of formula names whose parameters should be removed
#' @param keep_parameters character vector of parameter names to keep (all others removed unless in keep_formulas).
#' If "bias" is specified and the bias prior contains PET, PEESE, or weightfunction priors,
#' the corresponding parameters (PET, PEESE, omega) are also added to the keep list.
#' @param keep_formulas character vector of formula names whose parameters should be kept (all others removed unless in keep_parameters)
#' @param remove_spike_0 whether to remove spike at 0 priors
#' @return list with filtered model_samples and prior_list
.filter_parameters <- function(prior_list, remove_parameters = NULL, remove_formulas = NULL,
                               keep_parameters = NULL, keep_formulas = NULL, remove_spike_0 = TRUE) {
  
  # get formula parameter for each prior
  prior_formulas <- sapply(prior_list, function(p) {
    form <- attr(p, "parameter")
    if (is.null(form)) "__none" else form

  })
  
  # helper function to get bias-related parameters from a bias prior
  .get_bias_params <- function(prior_list, bias_name = "bias") {
    bias_params <- character(0)
    if (bias_name %in% names(prior_list)) {
      bias_params <- .selection_bias_parameter_names(prior_list[[bias_name]])
    }
    return(bias_params)
  }
  
  # initialize parameters to remove
  params_to_remove <- character(0)
  
  # handle remove_spike_0
  if (remove_spike_0) {
    spike_0_params <- names(prior_list)[sapply(seq_along(prior_list), function(i) {
      is.prior.point(prior_list[[i]]) && prior_list[[i]][["parameters"]][["location"]] == 0
    })]
    params_to_remove <- c(params_to_remove, spike_0_params)
  }
  
  # handle remove_parameters
  if (is.logical(remove_parameters) && isTRUE(remove_parameters)) {
    # remove all non-formula parameters
    non_formula_params <- names(prior_list)[prior_formulas == "__none"]
    params_to_remove <- c(params_to_remove, non_formula_params)
  } else if (is.character(remove_parameters)) {
    params_to_remove <- c(params_to_remove, remove_parameters)
    # if "bias" is in remove_parameters, also add corresponding bias-related parameters
    if ("bias" %in% remove_parameters) {
      params_to_remove <- c(params_to_remove, .get_bias_params(prior_list, "bias"))
    }
  }
  
  # handle remove_formulas
  if (!is.null(remove_formulas)) {
    formula_params <- names(prior_list)[prior_formulas %in% remove_formulas]
    params_to_remove <- c(params_to_remove, formula_params)
  }
  
  # handle keep_parameters and keep_formulas (these define what to keep, everything else is removed)
  if (!is.null(keep_parameters) || !is.null(keep_formulas)) {
    # start with all parameters as candidates for removal
    all_params <- names(prior_list)
    
    # determine which parameters to keep
    params_to_keep <- character(0)
    
    if (!is.null(keep_parameters)) {
      params_to_keep <- c(params_to_keep, keep_parameters)
      # if "bias" is in keep_parameters, also add corresponding bias-related parameters
      if ("bias" %in% keep_parameters) {
        params_to_keep <- c(params_to_keep, .get_bias_params(prior_list, "bias"))
      }
    }
    
    if (!is.null(keep_formulas)) {
      formula_params_to_keep <- names(prior_list)[prior_formulas %in% keep_formulas]
      params_to_keep <- c(params_to_keep, formula_params_to_keep)
    }
    
    # add parameters not in keep list to removal list
    params_not_kept <- all_params[!all_params %in% params_to_keep]
    params_to_remove <- c(params_to_remove, params_not_kept)
    
    # if "bias" is in params_not_kept, also add corresponding bias-related parameters
    if ("bias" %in% params_not_kept) {
      params_to_remove <- c(params_to_remove, .get_bias_params(prior_list, "bias"))
    }
  }
  
  # remove duplicates

  params_to_remove <- unique(params_to_remove)
  
  return(params_to_remove)
}


#' @rdname posterior_extraction_helpers
#' @param par parameter name
#' @param conditional whether to compute conditional summary
#' @param remove_inclusion whether to remove inclusion indicators
#' @param warnings character vector for collecting warnings
#' @return list with updated model_samples, prior_list, and warnings
.process_spike_and_slab <- function(model_samples, prior_list, par, conditional = FALSE, remove_inclusion = FALSE, warnings = NULL) {
  
  # prepare parameter names
  if (is.prior.factor(.get_spike_and_slab_variable(prior_list[[par]]))) {
    if (.get_prior_factor_levels(.get_spike_and_slab_variable(prior_list[[par]])) == 1) {
      par_names <- par
    } else {
      par_names <- paste0(par, "[", 1:.get_prior_factor_levels(.get_spike_and_slab_variable(prior_list[[par]])), "]")
    }
  } else {
    par_names <- par
  }
  
  # change the samples between conditional/averaged based on the preferences
  if (conditional) {
    # compute the number of conditional samples
    n_conditional_samples <- sum(model_samples[, colnames(model_samples) == paste0(par, "_indicator")] == 1)
    
    # replace null samples with NAs (important for later transformations)
    model_samples[model_samples[, colnames(model_samples) == paste0(par, "_indicator")] != 1, par_names] <- NA
    
    # add warnings about conditional summary
    warnings <- c(warnings, .runjags_conditional_warning(par_names, n_conditional_samples))
  }
  
  # remove the inclusion
  model_samples <- model_samples[, colnames(model_samples) != paste0(par, "_inclusion"), drop = FALSE]
  
  # remove the latent variable
  model_samples <- model_samples[, !colnames(model_samples) %in% gsub(par, paste0(par, "_variable"), par_names), drop = FALSE]
  
  # remove/rename the inclusions probabilities
  if (remove_inclusion) {
    model_samples <- model_samples[, colnames(model_samples) != paste0(par, "_indicator"), drop = FALSE]
  } else {
    colnames(model_samples)[colnames(model_samples) == paste0(par, "_indicator")] <- paste0(par, " (inclusion)")
  }
  
  # modify the parameter list (forward the parameter attribute)
  variable_component <- .get_spike_and_slab_variable(prior_list[[par]])
  attr(variable_component, "parameter") <- attr(prior_list[[par]], "parameter")
  prior_list[[par]] <- variable_component
  
  return(list(model_samples = model_samples, prior_list = prior_list, warnings = warnings))
}


#' @rdname posterior_extraction_helpers
#' @param transformations list of transformations to apply
#' @param transform_factors whether orthonormal/meandif will be transformed later
#' @return updated model_samples matrix
.apply_parameter_transformations <- function(model_samples, transformations, prior_list, transform_factors = FALSE) {
  
  if (is.null(transformations)) {
    return(model_samples)
  }
  
  for (par in names(transformations)) {
    if (!is.prior.factor(prior_list[[par]])) {
      # non-factor priors
      model_samples[, par] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[, par]), transformations[[par]][["arg"]]))
    } else if ((!transform_factors && (is.prior.orthonormal(prior_list[[par]]) || is.prior.meandif(prior_list[[par]]))) || is.prior.treatment(prior_list[[par]])) {
      # treatment priors, or orthonormal/meandif that won't be transformed to differences
      par_names <- .JAGS_prior_factor_names(par, prior_list[[par]])
      
      for (i in seq_along(par_names)) {
        model_samples[, par_names[i]] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[, par_names[i]]), transformations[[par]][["arg"]]))
      }
    }
  }
  
  return(model_samples)
}


#' @rdname posterior_extraction_helpers
#' @param transform_factors whether to transform orthonormal/meandif to differences
#' @return updated model_samples matrix
.transform_factor_contrasts <- function(model_samples, prior_list, transform_factors = FALSE, transformations = NULL) {
  
  factor_parameters <- names(prior_list)[sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x))]

  if (!transform_factors || length(factor_parameters) == 0) {
    return(model_samples)
  }

  if (any(factor_parameters %in% names(transformations))) {
    message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the orthonormal/meandif contrasts to the differences from the mean.")
  }

  for (par in factor_parameters) {
    
    prior_list[[par]] <- .add_factor_metadata_from_named_objects(prior_list[[par]], par, prior_list)
    par_names <- .JAGS_prior_factor_names(par, prior_list[[par]])
    
    temp_position <- min(which(colnames(model_samples) %in% par_names))
    temp_samples  <- model_samples[, colnames(model_samples) %in% par_names, drop = FALSE]
    model_samples <- model_samples[, !colnames(model_samples) %in% par_names, drop = FALSE]
    
    transformed_samples <- .transform_factor_contrast_samples(
      coefficient_samples = temp_samples,
      metadata            = prior_list[[par]],
      parameter           = par,
      transformed_class   = if(is.prior.orthonormal(prior_list[[par]])) {
        "mixed_posteriors.orthonormal_transformed"
      }else{
        "mixed_posteriors.meandif_transformed"
      }
    )
    
    # apply transformation if specified
    if (!is.null(transformations[[par]])) {
      for (i in seq_len(ncol(transformed_samples))) {
        transformed_samples[, i] <- do.call(transformations[[par]][["fun"]], c(list(transformed_samples[, i]), transformations[[par]][["arg"]]))
      }
    }
    
    # place the transformed samples back
    model_samples <- cbind(
      if (temp_position > 1) model_samples[, 1:(temp_position - 1), drop = FALSE],
      transformed_samples,
      if (temp_position <= ncol(model_samples)) model_samples[, temp_position:ncol(model_samples), drop = FALSE]
    )
  }
  
  return(model_samples)
}


# Helper: Attach factor levels to the factor components inside an interaction
# term instead of appending the level to the full interaction string.
.format_factor_level_parameter_names <- function(parameter, level_names, n_parameters = NULL) {

  if (!is.list(level_names)) {
    return(paste0(parameter, "[", level_names, "]"))
  }

  parameter_terms <- strsplit(parameter, "__xXx__", fixed = TRUE)[[1]]
  factor_terms <- names(level_names)
  factor_positions <- vapply(factor_terms, function(factor_term) {
    factor_position <- which(
      parameter_terms == factor_term |
        endsWith(parameter_terms, paste0("_", factor_term))
    )

    if (length(factor_position) == 1) {
      return(factor_position)
    }

    return(NA_integer_)
  }, integer(1))

  if (all(!is.na(factor_positions))) {
    level_grid <- expand.grid(level_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    formatted_names <- vapply(seq_len(nrow(level_grid)), function(i) {
      formatted_terms <- parameter_terms

      for (factor_term in factor_terms) {
        formatted_terms[factor_positions[[factor_term]]] <- sub(
          paste0(factor_term, "$"),
          paste0(factor_term, "[", level_grid[[factor_term]][i], "]"),
          formatted_terms[factor_positions[[factor_term]]]
        )
      }

      paste0(formatted_terms, collapse = "__xXx__")
    }, character(1))

    if (is.null(n_parameters) || length(formatted_names) == n_parameters) {
      return(formatted_names)
    }
  }

  fallback_names <- paste0(parameter, "[", unlist(level_names, use.names = FALSE), "]")
  if (!is.null(n_parameters) && length(fallback_names) != n_parameters) {
    return(paste0(parameter, "[", seq_len(n_parameters), "]"))
  }

  return(fallback_names)
}


#' @rdname posterior_extraction_helpers
#' @return updated model_samples matrix with renamed columns
.rename_factor_levels <- function(model_samples, prior_list) {
  
  # rename treatment factor levels
  if (any(sapply(prior_list, is.prior.treatment))) {
    for (par in names(prior_list)[sapply(prior_list, is.prior.treatment)]) {
      if (!.is_prior_interaction(prior_list[[par]])) {
        renamed_levels <- .format_factor_level_parameter_names(
          par,
          .get_prior_factor_level_names(prior_list[[par]])[-1],
          .get_prior_factor_levels(prior_list[[par]])
        )
        if (.get_prior_factor_levels(prior_list[[par]]) == 1) {
          colnames(model_samples)[colnames(model_samples) == par] <-
            renamed_levels[1]
        } else {
          colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
            renamed_levels
        }
      } else if (length(attr(prior_list[[par]], "levels")) == 1) {
        interaction_level_names <- .get_prior_factor_level_names(prior_list[[par]])
        interaction_level_names <- lapply(interaction_level_names, function(level_name) level_name[-1])
        colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
          .format_factor_level_parameter_names(par, interaction_level_names, .get_prior_factor_levels(prior_list[[par]]))
      }
    }
  }
  
  # rename independent factor levels
  if (any(sapply(prior_list, is.prior.independent))) {
    for (par in names(prior_list)[sapply(prior_list, is.prior.independent)]) {
      if (!.is_prior_interaction(prior_list[[par]])) {
        renamed_levels <- .format_factor_level_parameter_names(
          par,
          .get_prior_factor_level_names(prior_list[[par]]),
          .get_prior_factor_levels(prior_list[[par]])
        )
        if (.get_prior_factor_levels(prior_list[[par]]) == 1) {
          colnames(model_samples)[colnames(model_samples) == par] <-
            renamed_levels[1]
        } else {
          colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
            renamed_levels
        }
      } else if (length(attr(prior_list[[par]], "levels")) == 1) {
        colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
          .format_factor_level_parameter_names(par, .get_prior_factor_level_names(prior_list[[par]]), .get_prior_factor_levels(prior_list[[par]]))
      }
    }
  }
  
  return(model_samples)
}
